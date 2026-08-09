// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.archaeopteryx.tools;

import java.util.List;

import org.forester.archaeopteryx.tools.NodeDataImporter.ImportResult;
import org.forester.archaeopteryx.tools.NodeDataImporter.MatchBy;
import org.forester.archaeopteryx.tools.NodeDataImporter.MatchReport;
import org.forester.archaeopteryx.tools.NodeDataImporter.Table;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headless tests for {@link NodeDataImporter}: the export -> blank -> import round trip, custom columns
 * becoming node properties, the overwrite / never-clobber-with-blank merge policy, the node_id key path,
 * unmatched-row / tip-not-in-table accounting, the deliberately-ignored branch_length column, and
 * per-cell validation warnings.
 */
public final class NodeDataImporterTest {

    public static void main( final String[] args ) {
        System.out.println( "NodeDataImporter: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return roundTrip() && customColumns() && mergePolicy() && nodeIdKey() && accounting()
                    && branchLengthIgnored() && validationWarnings() && structuralErrors() && duplicateNames()
                    && normalizationAndTrim() && blankKeys() && resultImmutability()
                    && csvParsing() && delimiterDetect() && matchByAccession() && matchByTaxonomy() && dryRunCounts()
                    && embeddedNewlineCsv() && explicitDelimiter() && columnPlan()
                    && nodeIdColumnSkipped() && leadingBlankLine() && emptyRenameSkipped() && dryRunValidatesKeyCol()
                    && userMatchOptionsAndSummary() && defaultKeyPrefersName() && importProfile();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected exception: " + e );
        }
    }

    /** Export a fully-populated tree, then import that TSV onto a bare copy and assert the data is restored. */
    private static boolean roundTrip() throws Exception {
        final Phylogeny source = richTree();
        final String tsv = NodeDataExporter.toNodeDataTsv( source );
        final Phylogeny target = bareTree(); // tips A, B, C with names only
        final ImportResult res = NodeDataImporter.apply( target, tsv );
        if ( res.getTipsAnnotated() != 3 ) {
            return fail( "round trip should annotate all 3 tips, got " + res.getTipsAnnotated() );
        }
        if ( ( res.getRowsMatched() != 3 ) || ( res.getTipsNotInTable() != 0 ) || !res.getUnmatchedRowKeys().isEmpty() ) {
            return fail( "round trip accounting wrong: " + res.summary() );
        }
        final PhylogenyNode a = tip( target, "A" );
        if ( !"Orthoflavivirus zikae".equals( a.getNodeData().getTaxonomy().getScientificName() )
                || !"Zika virus".equals( a.getNodeData().getTaxonomy().getCommonName() )
                || !"9FLAV".equals( a.getNodeData().getTaxonomy().getTaxonomyCode() )
                || !"64320".equals( a.getNodeData().getTaxonomy().getIdentifier().getValue() )
                || !"species".equals( a.getNodeData().getTaxonomy().getRank() ) ) {
            return fail( "taxonomy not round-tripped onto A: " + a.getNodeData().getTaxonomy() );
        }
        if ( !"polyprotein".equals( a.getNodeData().getSequence().getName() )
                || !"POLY".equals( a.getNodeData().getSequence().getGeneName() )
                || !"PP".equals( a.getNodeData().getSequence().getSymbol() )
                || !"PQ48392".equals( a.getNodeData().getSequence().getAccession().getValue() )
                || !"protein".equals( a.getNodeData().getSequence().getType() ) ) {
            return fail( "sequence not round-tripped onto A: " + a.getNodeData().getSequence() );
        }
        if ( !"Aedes".equals( propertyValue( a, "data:host" ) ) ) {
            return fail( "property not round-tripped onto A: " + propertyValue( a, "data:host" ) );
        }
        // branch_length is exported but, by design, never applied -> target keeps its default branch length
        if ( a.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
            return fail( "branch_length must NOT be applied by an annotation import: " + a.getDistanceToParent() );
        }
        return true;
    }

    /** A non-reserved column becomes a node property; a colon-less header is namespaced with data:. */
    private static boolean customColumns() throws Exception {
        final Phylogeny phy = bareTree();
        final String tsv = "name\tcountry\tdata:host\n" + "A\tBrazil\tAedes\n" + "B\tUSA\t\n";
        final ImportResult res = NodeDataImporter.apply( phy, tsv );
        if ( !"Brazil".equals( propertyValue( tip( phy, "A" ), "data:country" ) ) ) {
            return fail( "colon-less custom column should become data:country" );
        }
        if ( !"Aedes".equals( propertyValue( tip( phy, "A" ), "data:host" ) ) ) {
            return fail( "namespaced custom column should be used verbatim" );
        }
        // B's blank data:host cell must not create an (empty) property
        if ( propertyValue( tip( phy, "B" ), "data:host" ) != null ) {
            return fail( "a blank custom cell must not create a property" );
        }
        final List<String> cols = res.getPropertyColumns();
        if ( !cols.contains( "data:country" ) || !cols.contains( "data:host" ) ) {
            return fail( "property columns should be reported for coloring: " + cols );
        }
        if ( !res.summary().contains( "coloring" ) ) {
            return fail( "summary should advertise the color-able columns: " + res.summary() );
        }
        return true;
    }

    /** A non-empty cell overwrites; a blank cell leaves the existing value alone. */
    private static boolean mergePolicy() throws Exception {
        final Phylogeny phy = bareTree();
        final PhylogenyNode a = tip( phy, "A" );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( "Old name" );
        a.getNodeData().addTaxonomy( t );
        final Sequence s = new Sequence();
        s.setGeneName( "KEEPME" );
        a.getNodeData().addSequence( s );
        // taxonomy_scientific_name carries a value (overwrite), gene_name is blank (must not clobber)
        final String tsv = "name\ttaxonomy_scientific_name\tgene_name\n" + "A\tNew name\t\n";
        NodeDataImporter.apply( phy, tsv );
        if ( !"New name".equals( a.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "non-empty cell should overwrite: " + a.getNodeData().getTaxonomy().getScientificName() );
        }
        if ( !"KEEPME".equals( a.getNodeData().getSequence().getGeneName() ) ) {
            return fail( "blank cell must not clobber an existing value: " + a.getNodeData().getSequence().getGeneName() );
        }
        // overwriting a property replaces (not appends) the value for that ref
        final String tsv1 = "name\tdata:clade\n" + "A\tI\n";
        final String tsv2 = "name\tdata:clade\n" + "A\tII\n";
        NodeDataImporter.apply( phy, tsv1 );
        NodeDataImporter.apply( phy, tsv2 );
        if ( a.getNodeData().getProperties().getPropertiesWithGivenRef( "data:clade" ).size() != 1
                || !"II".equals( propertyValue( a, "data:clade" ) ) ) {
            return fail( "re-importing a property should overwrite, not duplicate: "
                    + a.getNodeData().getProperties().getPropertiesWithGivenRef( "data:clade" ) );
        }
        return true;
    }

    /** When a node_id column is present it is the key, and the name column is then an applied field. */
    private static boolean nodeIdKey() throws Exception {
        final Phylogeny phy = bareTree();
        final PhylogenyNode b = tip( phy, "B" );
        final String tsv = "node_id\tname\ttaxonomy_code\n" + b.getId() + "\tRENAMED\tVZ9SC\n";
        final ImportResult res = NodeDataImporter.apply( phy, tsv );
        if ( !"RENAMED".equals( b.getName() ) ) {
            return fail( "name column should be applied when node_id is the key: " + b.getName() );
        }
        if ( !"VZ9SC".equals( b.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
            return fail( "taxonomy_code not applied via node_id key" );
        }
        if ( res.getRowsMatched() != 1 ) {
            return fail( "node_id row should match exactly one tip" );
        }
        return true;
    }

    /** Rows with no matching tip, and tips absent from the table, are counted. */
    private static boolean accounting() throws Exception {
        final Phylogeny phy = bareTree(); // A, B, C
        final String tsv = "name\tdata:x\n" + "A\t1\n" + "ZZZ\t9\n"; // B, C absent; ZZZ matches nothing
        final ImportResult res = NodeDataImporter.apply( phy, tsv );
        if ( res.getTipsAnnotated() != 1 ) {
            return fail( "only A should be annotated: " + res.getTipsAnnotated() );
        }
        if ( ( res.getUnmatchedRowKeys().size() != 1 ) || !res.getUnmatchedRowKeys().contains( "ZZZ" ) ) {
            return fail( "ZZZ row should be unmatched: " + res.getUnmatchedRowKeys() );
        }
        if ( res.getTipsNotInTable() != 2 ) {
            return fail( "B and C should be reported as not-in-table: " + res.getTipsNotInTable() );
        }
        return true;
    }

    private static boolean branchLengthIgnored() throws Exception {
        final Phylogeny phy = bareTree();
        final PhylogenyNode a = tip( phy, "A" );
        NodeDataImporter.apply( phy, "name\tbranch_length\n" + "A\t0.75\n" );
        if ( a.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
            return fail( "branch_length column must be ignored: " + a.getDistanceToParent() );
        }
        return true;
    }

    /** A value the model rejects (illegal rank) is recorded as a warning, not thrown, and the row still applies. */
    private static boolean validationWarnings() throws Exception {
        final Phylogeny phy = bareTree();
        final PhylogenyNode a = tip( phy, "A" );
        final String tsv = "name\ttaxonomy_rank\ttaxonomy_scientific_name\n"
                + "A\tnotarank\tFoo bar\n"  // bad rank, good sci name: A keeps a taxonomy (with the sci name)
                + "B\tbadrank\t\n";          // bad rank is B's ONLY taxonomy cell: B must get no taxonomy at all
        final ImportResult res = NodeDataImporter.apply( phy, tsv );
        if ( res.getWarnings().size() != 2 ) {
            return fail( "each illegal rank should produce a warning: " + res.getWarnings() );
        }
        if ( !"Foo bar".equals( a.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "a bad cell must not prevent the other cells in the row from applying" );
        }
        // the rejected-only cell must NOT leave a phantom empty taxonomy (which would flip isHasTaxonomy on)
        if ( tip( phy, "B" ).getNodeData().isHasTaxonomy() ) {
            return fail( "a rejected-only taxonomy cell must not attach an empty taxonomy" );
        }
        return true;
    }

    private static boolean structuralErrors() throws Exception {
        if ( !throwsIAE( "foo\tbar\nA\t1\n" ) ) {
            return fail( "a table with no name/node_id column should throw" );
        }
        if ( !throwsIAE( "" ) || !throwsIAE( "   \t  \n" ) ) {
            return fail( "an empty or all-whitespace header should throw" );
        }
        // a UTF-8 BOM on the header must not hide the key column -- and the data must actually apply through it
        final Phylogeny phy = bareTree();
        NodeDataImporter.apply( phy, "﻿name\tdata:x\nA\t1\n" );
        if ( !"1".equals( propertyValue( tip( phy, "A" ), "data:x" ) ) ) {
            return fail( "a BOM-prefixed header should be parsed and its data applied" );
        }
        return true;
    }

    /** Headers are trimmed and matched case-insensitively; a tip name is trimmed before matching. */
    private static boolean normalizationAndTrim() throws Exception {
        // surrounding whitespace + mixed case on the headers; a custom column keeps its original ref
        final Phylogeny phy = bareTree();
        final String tsv = " NAME \t Taxonomy_Code \t Data:Region \n" + "A\t9FLAV\tWest\n";
        NodeDataImporter.apply( phy, tsv );
        final PhylogenyNode a = tip( phy, "A" );
        if ( !"9FLAV".equals( a.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
            return fail( "a whitespaced/mixed-case reserved header should still map to its field" );
        }
        if ( !"West".equals( propertyValue( a, "Data:Region" ) ) ) {
            return fail( "a custom column ref should keep its original (un-lowercased) header" );
        }
        // a tip whose name carries surrounding whitespace still matches a trimmed row key
        final Phylogeny spaced = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode sp = new PhylogenyNode();
        sp.setName( " sp " );
        root.addAsChild( sp );
        spaced.setRoot( root );
        spaced.externalNodesHaveChanged();
        final ImportResult res = NodeDataImporter.apply( spaced, "name\tdata:x\nsp\t9\n" );
        if ( !"9".equals( propertyValue( sp, "data:x" ) ) || ( res.getTipsNotInTable() != 0 ) ) {
            return fail( "a whitespace-padded tip name should match a trimmed row key: " + res.summary() );
        }
        // a short row (fewer cells than the header) applies what it has and does not crash
        final Phylogeny shortrow = bareTree();
        NodeDataImporter.apply( shortrow, "name\ttaxonomy_code\tdata:x\nA\n" ); // only the key cell present
        if ( tip( shortrow, "A" ).getNodeData().isHasTaxonomy() ) {
            return fail( "a short row must not fabricate data for its missing cells" );
        }
        return true;
    }

    /** Blank-keyed rows and blank-named tips are accounted, not matched. */
    private static boolean blankKeys() throws Exception {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode named = new PhylogenyNode();
        named.setName( "A" );
        final PhylogenyNode blank = new PhylogenyNode(); // no name
        root.addAsChild( named );
        root.addAsChild( blank );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        final ImportResult res = NodeDataImporter.apply( phy, "name\tdata:x\nA\t1\n\t9\n" ); // 2nd row has a blank key
        if ( !res.getUnmatchedRowKeys().contains( "(blank)" ) ) {
            return fail( "a blank-key row should be reported as unmatched (blank): " + res.getUnmatchedRowKeys() );
        }
        if ( res.getTipsNotInTable() != 1 ) { // the blank-named tip is never keyed
            return fail( "the blank-named tip should count as not-in-table: " + res.getTipsNotInTable() );
        }
        return true;
    }

    /** The result's lists are immutable, defensive copies. */
    private static boolean resultImmutability() throws Exception {
        final ImportResult res = NodeDataImporter.apply( bareTree(), "name\tdata:x\nA\t1\nZZZ\t9\n" );
        for ( final List<String> list : List.of( res.getPropertyColumns(), res.getUnmatchedRowKeys(),
                res.getWarnings() ) ) {
            try {
                list.add( "x" );
                return fail( "result lists should be immutable" );
            }
            catch ( final UnsupportedOperationException expected ) {
                // good
            }
        }
        return true;
    }

    /** One row keyed on a name shared by two tips annotates both; the row is matched once. */
    private static boolean duplicateNames() throws Exception {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a1 = new PhylogenyNode();
        a1.setName( "dup" );
        final PhylogenyNode a2 = new PhylogenyNode();
        a2.setName( "dup" );
        root.addAsChild( a1 );
        root.addAsChild( a2 );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        final ImportResult res = NodeDataImporter.apply( phy, "name\tdata:x\ndup\thit\n" );
        if ( !"hit".equals( propertyValue( a1, "data:x" ) ) || !"hit".equals( propertyValue( a2, "data:x" ) ) ) {
            return fail( "a duplicate-name row should annotate every matching tip" );
        }
        if ( ( res.getRowsMatched() != 1 ) || ( res.getTipsAnnotated() != 2 ) ) {
            return fail( "rows matched=" + res.getRowsMatched() + " tips annotated=" + res.getTipsAnnotated() );
        }
        return true;
    }

    /** A comma-delimited table is detected as CSV and its quoted fields (embedded comma, doubled-quote escape) unquote. */
    private static boolean csvParsing() throws Exception {
        final Phylogeny phy = bareTree();
        final String csv = "name,data:city,data:note\n"
                + "A,\"Rio, RJ\",\"she said \"\"hi\"\"\"\n"
                + "B,Lima,plain\n";
        final Table t = NodeDataImporter.parseTable( csv );
        if ( t.getDelimiter() != ',' ) {
            return fail( "a comma header should be detected as CSV: " + t.getDelimiter() );
        }
        final ImportResult res = NodeDataImporter.apply( phy, t, t.defaultKeyColumn(), MatchBy.TIP_NAME );
        if ( !"Rio, RJ".equals( propertyValue( tip( phy, "A" ), "data:city" ) ) ) {
            return fail( "a quoted CSV field with an embedded comma should be preserved: "
                    + propertyValue( tip( phy, "A" ), "data:city" ) );
        }
        if ( !"she said \"hi\"".equals( propertyValue( tip( phy, "A" ), "data:note" ) ) ) {
            return fail( "a doubled-quote CSV escape should unquote to one quote: "
                    + propertyValue( tip( phy, "A" ), "data:note" ) );
        }
        if ( !"Lima".equals( propertyValue( tip( phy, "B" ), "data:city" ) ) || ( res.getTipsAnnotated() != 2 ) ) {
            return fail( "a plain CSV row should import: " + res.summary() );
        }
        return true;
    }

    /** Delimiter auto-detect: tab header -> TSV, comma header -> CSV, and tab wins when the header has both. */
    private static boolean delimiterDetect() throws Exception {
        if ( NodeDataImporter.parseTable( "name\tdata:x\nA\t1\n" ).getDelimiter() != '\t' ) {
            return fail( "a tab in the header should select TAB" );
        }
        if ( NodeDataImporter.parseTable( "name,data:x\nA,1\n" ).getDelimiter() != ',' ) {
            return fail( "a comma header (no tab) should select comma" );
        }
        // a header carrying both a tab AND a comma prefers tab (so a comma inside a header value is not a delimiter)
        final Table both = NodeDataImporter.parseTable( "name\tdata:x,y\nA\t1\n" );
        if ( ( both.getDelimiter() != '\t' ) || ( both.getColumnCount() != 2 )
                || !"data:x,y".equals( both.getHeaders()[ 1 ] ) ) {
            return fail( "tab should win over comma: cols=" + both.getColumnCount() );
        }
        return true;
    }

    /** Keyed join on the sequence accession attribute (tip names differ from the table key). */
    private static boolean matchByAccession() throws Exception {
        final Phylogeny phy = accessionTree(); // tips n1/ACC1, n2/ACC2
        final Table t = NodeDataImporter.parseTable( "acc,data:host\nACC1,cat\nACC2,dog\n" );
        final ImportResult res = NodeDataImporter.apply( phy, t, 0, MatchBy.SEQUENCE_ACCESSION );
        if ( res.getTipsAnnotated() != 2 ) {
            return fail( "an accession join should annotate both tips: " + res.summary() );
        }
        if ( !"cat".equals( propertyValue( tipByAccession( phy, "ACC1" ), "data:host" ) ) ) {
            return fail( "accession join wrote the wrong value" );
        }
        return true;
    }

    /** Keyed join on taxonomy id and on taxonomy scientific name. */
    private static boolean matchByTaxonomy() throws Exception {
        final Phylogeny phy = taxonomyTree(); // 9606/Homo sapiens, 10090/Mus musculus
        final ImportResult by_id = NodeDataImporter.apply( phy,
                NodeDataImporter.parseTable( "taxid,data:group\n9606,primate\n10090,rodent\n" ), 0, MatchBy.TAXONOMY_ID );
        if ( ( by_id.getTipsAnnotated() != 2 )
                || !"primate".equals( propertyValue( tipByTaxId( phy, "9606" ), "data:group" ) ) ) {
            return fail( "taxonomy-id join wrong: " + by_id.summary() );
        }
        final ImportResult by_sci = NodeDataImporter.apply( phy,
                NodeDataImporter.parseTable( "sci,data:g2\nHomo sapiens,H\nMus musculus,M\n" ), 0,
                MatchBy.TAXONOMY_SCIENTIFIC_NAME );
        if ( ( by_sci.getTipsAnnotated() != 2 )
                || !"H".equals( propertyValue( tipByTaxId( phy, "9606" ), "data:g2" ) ) ) {
            return fail( "taxonomy-scientific-name join wrong: " + by_sci.summary() );
        }
        return true;
    }

    /** dryRun reports matched / ambiguous / unmatched rows, tips-without-row, and the import columns -- no mutation. */
    private static boolean dryRunCounts() throws Exception {
        // tips: A, A (duplicate), B, C
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final String name : new String[] { "A", "A", "B", "C" } ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( name );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        final Table t = NodeDataImporter.parseTable( "name,data:x\nA,1\nB,2\nZZZ,9\n" );
        final MatchReport rep = NodeDataImporter.dryRun( phy, t, 0, MatchBy.TIP_NAME );
        if ( ( rep.getTotalRows() != 3 ) || ( rep.getRowsMatched() != 2 ) || ( rep.getRowsAmbiguous() != 1 )
                || ( rep.getRowsUnmatched() != 1 ) ) {
            return fail( "dry-run row counts wrong: " + rep.summaryLine() );
        }
        if ( ( rep.getTotalTips() != 4 ) || ( rep.getTipsWithoutRow() != 1 ) ) { // C has no row
            return fail( "dry-run tip counts wrong: " + rep.summaryLine() );
        }
        if ( ( rep.getPropertyColumns().size() != 1 ) || !rep.getPropertyColumns().contains( "data:x" ) ) {
            return fail( "dry-run should list the import column: " + rep.getPropertyColumns() );
        }
        // a dry run must NOT mutate the tree
        if ( tip( phy, "B" ).getNodeData().getProperties() != null ) {
            return fail( "dryRun must not write anything onto the tree" );
        }
        return true;
    }

    /** RFC-4180: a double-quoted CSV field may contain an embedded newline (the value spans two physical lines). */
    private static boolean embeddedNewlineCsv() throws Exception {
        final Phylogeny phy = bareTree();
        final String csv = "name,data:note\n" + "A,\"line one\nline two\"\n" + "B,plain\n";
        final Table t = NodeDataImporter.parseTable( csv );
        if ( t.getRowCount() != 2 ) {
            return fail( "an embedded-newline quoted field must NOT split the record: rows=" + t.getRowCount() );
        }
        NodeDataImporter.apply( phy, t, 0, MatchBy.TIP_NAME );
        if ( !"line one\nline two".equals( propertyValue( tip( phy, "A" ), "data:note" ) ) ) {
            return fail( "a quoted field should preserve its embedded newline: "
                    + propertyValue( tip( phy, "A" ), "data:note" ) );
        }
        if ( !"plain".equals( propertyValue( tip( phy, "B" ), "data:note" ) ) ) {
            return fail( "the row after an embedded-newline field should still parse" );
        }
        return true;
    }

    /** An explicit delimiter overrides auto-detect (forcing the "wrong" one collapses the table to a single column). */
    private static boolean explicitDelimiter() throws Exception {
        final Table forced_tab = NodeDataImporter.parseTable( "name,data:x\nA,1\n", Character.valueOf( '\t' ) );
        if ( ( forced_tab.getColumnCount() != 1 ) || !"name,data:x".equals( forced_tab.getHeaders()[ 0 ] ) ) {
            return fail( "forcing TAB on comma data should treat commas as literal (one column): "
                    + forced_tab.getColumnCount() );
        }
        final Table forced_comma = NodeDataImporter.parseTable( "name\tdata:x\nA\t1\n", Character.valueOf( ',' ) );
        if ( forced_comma.getColumnCount() != 1 ) {
            return fail( "forcing comma on a tab file should yield one column: " + forced_comma.getColumnCount() );
        }
        return true;
    }

    /** A ColumnPlan de-selects columns and renames others (a rename to a reserved name routes to that model field). */
    private static boolean columnPlan() throws Exception {
        final Phylogeny phy = bareTree();
        final Table t = NodeDataImporter.parseTable( "name,keepme,dropme,raw\nA,yes,no,7\n" );
        final NodeDataImporter.ColumnPlan plan = NodeDataImporter.ColumnPlan.importAll( t );
        plan.setIncluded( 2, false );   // drop "dropme"
        plan.setHeader( 3, "score" );   // rename "raw" -> a data:score property
        final ImportResult res = NodeDataImporter.apply( phy, t, 0, MatchBy.TIP_NAME, plan );
        final PhylogenyNode a = tip( phy, "A" );
        if ( !"yes".equals( propertyValue( a, "data:keepme" ) ) ) {
            return fail( "an included column should import" );
        }
        if ( propertyValue( a, "data:dropme" ) != null ) {
            return fail( "a de-selected column must NOT import" );
        }
        if ( !"7".equals( propertyValue( a, "data:score" ) ) || ( propertyValue( a, "data:raw" ) != null ) ) {
            return fail( "a renamed column should import under the new ref only" );
        }
        if ( ( res.getPropertyColumns().size() != 2 ) || !res.getPropertyColumns().contains( "data:score" ) ) {
            return fail( "the result should report the included/renamed property columns: " + res.getPropertyColumns() );
        }
        final MatchReport rep = NodeDataImporter.dryRun( phy, t, 0, MatchBy.TIP_NAME, plan );
        if ( rep.getPropertyColumns().contains( "dropme" ) || !rep.getPropertyColumns().contains( "score" ) ) {
            return fail( "dry-run should reflect the plan's include/rename: " + rep.getPropertyColumns() );
        }
        // renaming a column to a RESERVED name routes it to that model field instead of a property
        final Phylogeny phy2 = bareTree();
        final Table t2 = NodeDataImporter.parseTable( "name,code\nA,9FLAV\n" );
        final NodeDataImporter.ColumnPlan plan2 = NodeDataImporter.ColumnPlan.importAll( t2 );
        plan2.setHeader( 1, "taxonomy_code" );
        NodeDataImporter.apply( phy2, t2, 0, MatchBy.TIP_NAME, plan2 );
        if ( !"9FLAV".equals( tip( phy2, "A" ).getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
            return fail( "renaming a column to a reserved name should fill that model field" );
        }
        return true;
    }

    /** A non-key node_id column is recognized but never written, and must NOT count the tip as annotated. */
    private static boolean nodeIdColumnSkipped() throws Exception {
        final Phylogeny phy = bareTree(); // A, B, C
        final Table t = NodeDataImporter.parseTable( "node_id,name\n5,A\n" ); // node_id is a non-key, included column
        final ImportResult res = NodeDataImporter.apply( phy, t, 1, MatchBy.TIP_NAME ); // key = name (col 1)
        if ( res.getTipsAnnotated() != 0 ) {
            return fail( "a node_id column writes nothing and must not count as annotating a tip: "
                    + res.getTipsAnnotated() );
        }
        if ( tip( phy, "A" ).getNodeData().getProperties() != null ) {
            return fail( "a node_id column must not become a property" );
        }
        return true;
    }

    /** A CSV/TSV whose first line is blank still parses: the header is the first NON-blank record. */
    private static boolean leadingBlankLine() throws Exception {
        final Phylogeny phy = bareTree();
        final Table t = NodeDataImporter.parseTable( "\nname,data:x\nA,7\n" ); // leading blank line
        if ( ( t.getColumnCount() != 2 ) || !"name".equals( t.getHeaders()[ 0 ] ) ) {
            return fail( "a leading blank line should not become the header: cols=" + t.getColumnCount() );
        }
        NodeDataImporter.apply( phy, t, 0, MatchBy.TIP_NAME );
        if ( !"7".equals( propertyValue( tip( phy, "A" ), "data:x" ) ) ) {
            return fail( "a CSV with a leading blank line should still import its data" );
        }
        return true;
    }

    /** A blanked rename (empty effective header) is skipped, not written as a malformed "data:" property. */
    private static boolean emptyRenameSkipped() throws Exception {
        final Phylogeny phy = bareTree();
        final Table t = NodeDataImporter.parseTable( "name,keepme,blankme\nA,yes,no\n" );
        final NodeDataImporter.ColumnPlan plan = NodeDataImporter.ColumnPlan.importAll( t );
        plan.setHeader( 2, "" ); // blank the effective header of "blankme"
        NodeDataImporter.apply( phy, t, 0, MatchBy.TIP_NAME, plan );
        final PhylogenyNode a = tip( phy, "A" );
        if ( !"yes".equals( propertyValue( a, "data:keepme" ) ) ) {
            return fail( "an included column should import" );
        }
        if ( ( propertyValue( a, "data:" ) != null ) || ( propertyValue( a, "data:blankme" ) != null ) ) {
            return fail( "a blank effective header must NOT create a property (no malformed \"data:\" ref)" );
        }
        return true;
    }

    /** dryRun validates the key column like apply does (a negative index throws, not AIOOBE). */
    private static boolean dryRunValidatesKeyCol() throws Exception {
        final Phylogeny phy = bareTree();
        final Table t = NodeDataImporter.parseTable( "name,data:x\nA,1\n" );
        try {
            NodeDataImporter.dryRun( phy, t, -1, MatchBy.TIP_NAME );
            return fail( "dryRun should reject a negative key column with IllegalArgumentException" );
        }
        catch ( final IllegalArgumentException expected ) {
            return true;
        }
    }

    /** userMatchOptions offers the four tip attributes (never NODE_ID), and summaryLine formats the join. */
    private static boolean userMatchOptionsAndSummary() throws Exception {
        final MatchBy[] opts = NodeDataImporter.userMatchOptions();
        if ( ( opts.length != 4 ) || ( opts[ 0 ] != MatchBy.TIP_NAME ) ) {
            return fail( "userMatchOptions should list the 4 tip attributes, Tip name first" );
        }
        for( final MatchBy m : opts ) {
            if ( m == MatchBy.NODE_ID ) {
                return fail( "NODE_ID (internal, session-only) must not be user-selectable" );
            }
        }
        final Phylogeny phy = bareTree(); // A, B, C
        final MatchReport rep = NodeDataImporter.dryRun( phy,
                NodeDataImporter.parseTable( "name,data:x\nA,1\nZZZ,9\n" ), 0, MatchBy.TIP_NAME );
        final String line = rep.summaryLine();
        if ( !line.contains( "1/2 rows match" ) || !line.contains( "1 unmatched" ) || !line.contains( "tips have no row" )
                || !line.contains( "data:x" ) ) {
            return fail( "summaryLine format wrong: " + line );
        }
        return true;
    }

    /** The dialog's default key column prefers a "name" column over "node_id" (which it cannot match against). */
    private static boolean defaultKeyPrefersName() throws Exception {
        if ( NodeDataImporter.parseTable( "node_id,name\n5,A\n" ).defaultKeyColumn() != 1 ) {
            return fail( "defaultKeyColumn should prefer the name column (1) over node_id (0)" );
        }
        return true;
    }

    /** An ImportProfile (the re-import "annotation profile") remembers the mapping by HEADER NAME, so re-resolving it
     *  against a source with reordered / added columns keeps the right key, excludes, and renames. */
    private static boolean importProfile() throws Exception {
        final Table t = NodeDataImporter.parseTable( "name,host,country,reads\nA,cat,US,5\n" );
        final NodeDataImporter.ColumnPlan plan = NodeDataImporter.ColumnPlan.importAll( t );
        plan.setIncluded( 2, false ); // exclude "country"
        plan.setHeader( 3, "depth" );  // rename "reads" -> data:depth
        final NodeDataImporter.ImportProfile prof = NodeDataImporter.ImportProfile.from( t, 0, MatchBy.TIP_NAME, plan,
                "/x/data.csv", false );
        if ( prof.isUrl() || !"/x/data.csv".equals( prof.getSource() ) || ( prof.getMatchBy() != MatchBy.TIP_NAME )
                || ( prof.getDelimiter() == null ) || ( prof.getDelimiter().charValue() != ',' ) ) {
            return fail( "profile fields not captured: url=" + prof.isUrl() + " delim=" + prof.getDelimiter() );
        }
        // re-resolve against a source whose columns are REORDERED and gained a NEW column ("extra")
        final Table t2 = NodeDataImporter.parseTable( "reads,name,host,country,extra\n9,B,dog,UK,new\n" );
        if ( prof.keyColumn( t2 ) != 1 ) {
            return fail( "keyColumn should resolve 'name' by header name to index 1, got " + prof.keyColumn( t2 ) );
        }
        final NodeDataImporter.ColumnPlan p2 = prof.columnPlan( t2 );
        if ( p2.isIncluded( 1 ) ) {
            return fail( "the key column must stay excluded after re-resolving" );
        }
        if ( !p2.isIncluded( 2 ) || !"host".equals( p2.header( 2 ) ) ) {
            return fail( "host should stay included under its own header" );
        }
        if ( p2.isIncluded( 3 ) ) {
            return fail( "country should stay excluded across a re-parse (remembered by name)" );
        }
        if ( !p2.isIncluded( 0 ) || !"depth".equals( p2.header( 0 ) ) ) {
            return fail( "reads should stay renamed to depth across a re-parse" );
        }
        if ( !p2.isIncluded( 4 ) || !"extra".equals( p2.header( 4 ) ) ) {
            return fail( "a NEW column not in the profile should default-include" );
        }
        // apply end-to-end (tip B) and confirm the resolved mapping wrote the expected properties
        final Phylogeny phy = bareTree(); // A, B, C
        NodeDataImporter.apply( phy, t2, prof.keyColumn( t2 ), prof.getMatchBy(), p2 );
        final PhylogenyNode b = tip( phy, "B" );
        if ( !"dog".equals( propertyValue( b, "data:host" ) ) || !"9".equals( propertyValue( b, "data:depth" ) )
                || !"new".equals( propertyValue( b, "data:extra" ) ) || ( propertyValue( b, "data:country" ) != null ) ) {
            return fail( "re-resolved profile wrote the wrong properties onto B" );
        }
        // a source that lost the key header falls back to the table's default key column
        final Table t3 = NodeDataImporter.parseTable( "id,host\nZ,x\n" );
        if ( prof.keyColumn( t3 ) != t3.defaultKeyColumn() ) {
            return fail( "a missing key header should fall back to the default key column" );
        }
        return true;
    }

    // ---- helpers ----

    private static boolean throwsIAE( final String tsv ) {
        try {
            NodeDataImporter.apply( bareTree(), tsv );
            return false;
        }
        catch ( final IllegalArgumentException e ) {
            return true;
        }
    }

    private static Phylogeny richTree() throws Exception {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "A" );
        a.setDistanceToParent( 0.25 );
        final Taxonomy tax = new Taxonomy();
        tax.setScientificName( "Orthoflavivirus zikae" );
        tax.setCommonName( "Zika virus" );
        tax.setTaxonomyCode( "9FLAV" );
        tax.setIdentifier( new Identifier( "64320", "ncbi" ) );
        tax.setRank( "species" );
        a.getNodeData().addTaxonomy( tax );
        final Sequence seq = new Sequence();
        seq.setName( "polyprotein" );
        seq.setGeneName( "POLY" );
        seq.setSymbol( "PP" );
        seq.setAccession( new Accession( "PQ48392", "ncbi" ) );
        seq.setType( "protein" );
        a.getNodeData().addSequence( seq );
        final org.forester.phylogeny.data.PropertiesList props = new org.forester.phylogeny.data.PropertiesList();
        props.addProperty( new Property( "data:host", "Aedes", "", "xsd:string", Property.AppliesTo.NODE ) );
        a.getNodeData().setProperties( props );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "B" );
        b.getNodeData().addTaxonomy( taxon( "Dengue virus" ) );
        final PhylogenyNode c = new PhylogenyNode();
        c.setName( "C" );
        c.getNodeData().addTaxonomy( taxon( "West Nile virus" ) );
        root.addAsChild( a );
        root.addAsChild( b );
        root.addAsChild( c );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static Taxonomy taxon( final String sci_name ) {
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci_name );
        return t;
    }

    private static Phylogeny bareTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final String name : new String[] { "A", "B", "C" } ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( name );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode tip( final Phylogeny phy, final String name ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        throw new IllegalStateException( "no tip named " + name );
    }

    /** Two tips whose NAMES (n1/n2) differ from their sequence accessions (ACC1/ACC2) -- for accession-join tests. */
    private static Phylogeny accessionTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final String[] pair : new String[][] { { "n1", "ACC1" }, { "n2", "ACC2" } } ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( pair[ 0 ] );
            final Sequence s = new Sequence();
            s.setAccession( new Accession( pair[ 1 ] ) );
            n.getNodeData().addSequence( s );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode tipByAccession( final Phylogeny phy, final String acc ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( n.getNodeData().isHasSequence() && ( n.getNodeData().getSequence().getAccession() != null )
                    && acc.equals( n.getNodeData().getSequence().getAccession().getValue() ) ) {
                return n;
            }
        }
        throw new IllegalStateException( "no tip with accession " + acc );
    }

    /** Two tips carrying a taxonomy id + scientific name -- for taxonomy-join tests. */
    private static Phylogeny taxonomyTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final String[] pair : new String[][] { { "9606", "Homo sapiens" }, { "10090", "Mus musculus" } } ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "tip_" + pair[ 0 ] );
            final Taxonomy t = new Taxonomy();
            t.setIdentifier( new Identifier( pair[ 0 ], "ncbi" ) );
            t.setScientificName( pair[ 1 ] );
            n.getNodeData().addTaxonomy( t );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode tipByTaxId( final Phylogeny phy, final String id ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( n.getNodeData().isHasTaxonomy() && ( n.getNodeData().getTaxonomy().getIdentifier() != null )
                    && id.equals( n.getNodeData().getTaxonomy().getIdentifier().getValue() ) ) {
                return n;
            }
        }
        throw new IllegalStateException( "no tip with tax id " + id );
    }

    private static String propertyValue( final PhylogenyNode n, final String ref ) {
        if ( n.getNodeData().getProperties() == null ) {
            return null;
        }
        final List<Property> ps = n.getNodeData().getProperties().getPropertiesWithGivenRef( ref );
        return ( ( ps != null ) && !ps.isEmpty() ) ? ps.get( 0 ).getValue() : null;
    }

    private static boolean fail( final String message ) {
        System.out.println( "\nNodeDataImporterTest failed: " + message );
        return false;
    }

    private NodeDataImporterTest() {
    }
}
