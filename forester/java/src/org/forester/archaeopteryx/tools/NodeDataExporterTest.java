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

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headless tests for {@link NodeDataExporter}: FASTA of tip molecular sequences, and the tab-separated
 * tip-data table with its dynamic (only-populated) columns.
 */
public final class NodeDataExporterTest {

    public static void main( final String[] args ) {
        System.out.println( "NodeDataExporter: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        try {
            final Phylogeny phy = buildTree();

            // ---- FASTA: a record per tip with a molecular sequence, with a self-describing header ----
            final String fasta = NodeDataExporter.toFasta( phy );
            if ( !fasta.contains( "MKNPK" ) || !fasta.contains( "MGGKL" ) ) {
                return fail( "fasta missing a tip sequence:\n" + fasta );
            }
            // ZIKV_1's rich header: id, then the (differing) accession, organism, and " | " + description
            if ( !fasta.contains( ">ZIKV_1 PQ48392 Orthoflavivirus zikae | polyprotein" ) ) {
                return fail( "ZIKV_1 FASTA header not enriched as expected:\n" + fasta );
            }
            // ZIKV_2 has only a name + sequence -> a bare header (nothing accession/organism/desc appended)
            if ( !fasta.contains( ">ZIKV_2" ) || fasta.contains( ">ZIKV_2 " ) || fasta.contains( ">ZIKV_2 |" ) ) {
                return fail( "ZIKV_2 FASTA header should be bare:\n" + fasta );
            }
            // tip with no molecular sequence is not in the FASTA
            if ( fasta.contains( "ZIKV_3" ) ) {
                return fail( "tip without a molecular sequence must not appear in the FASTA" );
            }
            // record count = one per emitted sequence (ZIKV_1 + ZIKV_2; ZIKV_3 has none)
            if ( ( NodeDataExporter.fastaRecordCount( fasta ) != 2 ) || ( NodeDataExporter.fastaRecordCount( "" ) != 0 ) ) {
                return fail( "fastaRecordCount wrong: " + NodeDataExporter.fastaRecordCount( fasta ) );
            }

            // ---- TSV: header has only the populated columns, rows carry the values ----
            final String[] lines = NodeDataExporter.toNodeDataTsv( phy ).split( "\\R" );
            final List<String> header = Arrays.asList( lines[ 0 ].split( "\t", -1 ) );
            // present (some tip has a value)
            for ( final String col : new String[] { "name", "taxonomy_scientific_name", "taxonomy_code",
                    "taxonomy_id", "sequence_name", "gene_name", "sequence_accession", "branch_length",
                    "data:host" } ) {
                if ( !header.contains( col ) ) {
                    return fail( "expected column '" + col + "' in header: " + header );
                }
            }
            // absent (no tip has a value) -> column omitted
            for ( final String col : new String[] { "taxonomy_common_name", "taxonomy_rank", "sequence_symbol",
                    "sequence_type" } ) {
                if ( header.contains( col ) ) {
                    return fail( "empty column '" + col + "' should have been omitted: " + header );
                }
            }
            // header + 3 tip rows
            if ( lines.length != 4 ) {
                return fail( "expected a header + 3 tip rows; got " + lines.length );
            }
            // ZIKV_1's scientific-name and property cells
            final List<String> row1 = Arrays.asList( lines[ 1 ].split( "\t", -1 ) );
            if ( !"Orthoflavivirus zikae".equals( row1.get( header.indexOf( "taxonomy_scientific_name" ) ) ) ) {
                return fail( "ZIKV_1 scientific name wrong: " + row1 );
            }
            if ( !"Aedes".equals( row1.get( header.indexOf( "data:host" ) ) ) ) {
                return fail( "ZIKV_1 host property wrong: " + row1 );
            }
            // ZIKV_2 has no taxonomy -> empty scientific-name cell (column still present because ZIKV_1 has one)
            final List<String> row2 = Arrays.asList( lines[ 2 ].split( "\t", -1 ) );
            if ( !"".equals( row2.get( header.indexOf( "taxonomy_scientific_name" ) ) ) ) {
                return fail( "ZIKV_2 should have an empty scientific-name cell: " + row2 );
            }

            // ---- TSV key robustness ----
            // clean, distinct, non-empty names form a key -> no node_id column
            if ( header.contains( "node_id" ) || !NodeDataExporter.tipNamesFormUniqueKey( phy ) ) {
                return fail( "clean tip names should be the key, with no node_id column: " + header );
            }
            // a blank or duplicated name breaks the key -> a unique node_id key column is prepended
            final Phylogeny ambig = buildAmbiguousTree();
            if ( NodeDataExporter.tipNamesFormUniqueKey( ambig ) ) {
                return fail( "tree with a duplicate and a blank tip name must not form a unique key" );
            }
            final String[] alines = NodeDataExporter.toNodeDataTsv( ambig ).split( "\\R" );
            final List<String> aheader = Arrays.asList( alines[ 0 ].split( "\t", -1 ) );
            if ( !"node_id".equals( aheader.get( 0 ) ) || !"name".equals( aheader.get( 1 ) ) ) {
                return fail( "expected node_id then name as the first columns: " + aheader );
            }
            final Set<String> ids = new HashSet<>();
            for ( int i = 1; i < alines.length; i++ ) {
                final String id = alines[ i ].split( "\t", -1 )[ 0 ];
                if ( id.isEmpty() || !ids.add( id ) ) {
                    return fail( "node_id values must be present and unique: row " + i + " = " + alines[ i ] );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected exception: " + e );
        }
    }

    private static Phylogeny buildTree() throws Exception {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        // ZIKV_1: full data + a property
        final PhylogenyNode t1 = new PhylogenyNode();
        t1.setName( "ZIKV_1" );
        t1.setDistanceToParent( 0.25 );
        final Taxonomy tax = new Taxonomy();
        tax.setScientificName( "Orthoflavivirus zikae" );
        tax.setTaxonomyCode( "9FLAV" );
        tax.setIdentifier( new Identifier( "64320", "ncbi" ) );
        t1.getNodeData().addTaxonomy( tax );
        final Sequence s1 = new Sequence();
        s1.setName( "polyprotein" );
        s1.setGeneName( "POLY" );
        s1.setAccession( new Accession( "PQ48392", "ncbi" ) );
        s1.setMolecularSequence( "MKNPK" );
        t1.getNodeData().addSequence( s1 );
        final PropertiesList props = new PropertiesList();
        props.addProperty( new Property( "data:host", "Aedes", "", "xsd:string", AppliesTo.NODE ) );
        t1.getNodeData().setProperties( props );
        // ZIKV_2: just a name + molecular sequence (no taxonomy, no property)
        final PhylogenyNode t2 = new PhylogenyNode();
        t2.setName( "ZIKV_2" );
        final Sequence s2 = new Sequence();
        s2.setMolecularSequence( "MGGKL" );
        t2.getNodeData().addSequence( s2 );
        // ZIKV_3: name only, no sequence -> not in FASTA
        final PhylogenyNode t3 = new PhylogenyNode();
        t3.setName( "ZIKV_3" );
        root.addAsChild( t1 );
        root.addAsChild( t2 );
        root.addAsChild( t3 );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /** A tree whose tip names cannot key the table: two share a name and one is blank. */
    private static Phylogeny buildAmbiguousTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "dup" );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "dup" ); // duplicate name
        final PhylogenyNode c = new PhylogenyNode(); // blank name
        root.addAsChild( a );
        root.addAsChild( b );
        root.addAsChild( c );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String message ) {
        System.out.println( "\nNodeDataExporterTest failed: " + message );
        return false;
    }

    private NodeDataExporterTest() {
    }
}
