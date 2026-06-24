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
                    && normalizationAndTrim() && blankKeys() && resultImmutability();
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
