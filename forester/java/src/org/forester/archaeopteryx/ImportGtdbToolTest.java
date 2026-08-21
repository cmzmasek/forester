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

package org.forester.archaeopteryx;

import java.awt.GraphicsEnvironment;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.tools.NodeDataImporter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;

/**
 * Integration test for the "Import GTDB Taxonomy" tool: drives {@link MainFrame#importGtdbAndRefit} on a real
 * {@link MainFrame}/{@link TreePanel} -- an accession-named tree + a GTDB-Tk-style TSV -- and asserts the classification
 * column + key column are auto-detected, the tips gain gtdb:&lt;rank&gt; properties + a species taxonomy, the import
 * appends a provenance sentence, sets the edited flag, surfaces gtdb:phylum in "Color by", is UNDOABLE, and returns -1
 * for a non-GTDB table. Headful; a green no-op when headless.
 */
public final class ImportGtdbToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ImportGtdbTool: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( !columnDetectionOk() ) {
            return false; // pure, headless-safe
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // the import+refit+undo path needs a realized MainFrame/TreePanel
        }
        try {
            // two Bacteria in the SAME phylum (Pseudomonadota) + one Archaeon, so gtdb:phylum and gtdb:domain each have
            // a repeated value (distinct < leaves) and are OFFERED in Color-by (an all-unique column is filtered out)
            final String tsv = "user_genome\tclassification\n"
                    + "GCF_A\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n"
                    + "GCF_B\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas aeruginosa\n"
                    + "GCF_C\td__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrosopumilales;f__Nitrosopumilaceae;g__Nitrosopumilus;s__Nitrosopumilus maritimus\n";
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { genomeTree() }, conf, "gtdb" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Phylogeny live = tp.getPhylogeny();
                    // precondition: the accession-named tips carry no taxonomy yet
                    if ( tip( live, "GCF_A" ).getNodeData().isHasTaxonomy() ) {
                        fail( ok, "the genome tips should start without a taxonomy" );
                    }
                    final int annotated = frame.importGtdbAndRefit( live, tsv, "gtdb-classifications.tsv" );
                    if ( annotated != 3 ) {
                        fail( ok, "the GTDB import should annotate all 3 accession tips, got " + annotated );
                    }
                    final PhylogenyNode a = tip( live, "GCF_A" );
                    if ( !"Pseudomonadota".equals( propertyValue( a, "gtdb:phylum" ) )
                            || !"Escherichia".equals( propertyValue( a, "gtdb:genus" ) )
                            || !"Escherichia coli".equals( propertyValue( a, "gtdb:species" ) ) ) {
                        fail( ok, "GCF_A should carry the gtdb:* rank properties from its classification" );
                    }
                    if ( !a.getNodeData().isHasTaxonomy()
                            || !"Escherichia coli".equals( a.getNodeData().getTaxonomy().getScientificName() )
                            || !"species".equals( a.getNodeData().getTaxonomy().getRank() ) ) {
                        fail( ok, "GCF_A should gain a species-rank taxonomy from the most specific GTDB rank" );
                    }
                    // the archaeal tip's domain resolves distinctly (so Color-by gtdb:domain = Bacteria vs Archaea)
                    if ( !"Archaea".equals( propertyValue( tip( live, "GCF_C" ), "gtdb:domain" ) ) ) {
                        fail( ok, "GCF_C should be gtdb:domain=Archaea" );
                    }
                    final String desc = live.getDescription();
                    if ( ( desc == null ) || !desc.contains( "Imported GTDB taxonomy" ) ) {
                        fail( ok, "a GTDB import should append a provenance sentence: " + desc );
                    }
                    if ( !tp.isEdited() ) {
                        fail( ok, "a GTDB import that changed the tree should set the edited flag" );
                    }
                    final List<String> refs = tp.getControlPanel().colorByPropertyRefs();
                    if ( !refs.contains( "gtdb:phylum" ) || !refs.contains( "gtdb:domain" ) ) {
                        fail( ok, "the Color-by dropdown should now offer the gtdb:* ranks: " + refs );
                    }
                    // UNDO restores the pre-import tree (accession tips, no taxonomy, no gtdb:* property)
                    frame.undo();
                    final Phylogeny after_undo = frame.getMainPanel().getCurrentTreePanel().getPhylogeny();
                    final PhylogenyNode a2 = tip( after_undo, "GCF_A" );
                    if ( a2.getNodeData().isHasTaxonomy() || ( propertyValue( a2, "gtdb:phylum" ) != null ) ) {
                        fail( ok, "undo should restore the accession tree (no taxonomy / no gtdb:* property)" );
                    }
                    // NEGATIVE: a table with no GTDB classification column returns -1 and changes nothing
                    final Phylogeny cur = frame.getMainPanel().getCurrentTreePanel().getPhylogeny();
                    final int none = frame.importGtdbAndRefit( cur, "id\tvalue\nGCF_A\thello\nGCF_B\tworld\n", "plain.tsv" );
                    if ( none != -1 ) {
                        fail( ok, "a non-GTDB table should return -1, got " + none );
                    }
                    if ( tip( cur, "GCF_A" ).getNodeData().isHasTaxonomy() ) {
                        fail( ok, "a non-GTDB table must not mutate the tree" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Pure (headless-safe): the classification + key column auto-detection. The classification column is the one
     *  whose values look like GTDB, and a GTDB-Tk genome-id header wins the key over an unrelated "name" column. */
    private static boolean columnDetectionOk() {
        // standard GTDB-Tk: [user_genome, classification] -> class col 1, key col 0
        NodeDataImporter.Table t = NodeDataImporter
                .parseTable( "user_genome\tclassification\nGCF_A\td__Bacteria;p__Pseudomonadota;s__Escherichia coli\n" );
        if ( ( MainFrame.gtdbClassificationColumn( t ) != 1 ) || ( MainFrame.gtdbKeyColumn( t, 1 ) != 0 ) ) {
            return note( "GTDB-Tk table should detect class col 1 + key col 0 (user_genome), got class="
                    + MainFrame.gtdbClassificationColumn( t ) + " key=" + MainFrame.gtdbKeyColumn( t, 1 ) );
        }
        // a DECOY unrelated "name" column must NOT hijack the key -> user_genome (col 0) still wins
        t = NodeDataImporter.parseTable(
                "user_genome\tclassification\tname\nGCF_A\td__Bacteria;p__Pseudomonadota;s__Escherichia coli\tfoo\n" );
        if ( MainFrame.gtdbKeyColumn( t, MainFrame.gtdbClassificationColumn( t ) ) != 0 ) {
            return note( "a decoy 'name' column must not hijack the key from user_genome" );
        }
        // no genome-id header: a 'name' column IS the key (via defaultKeyColumn)
        t = NodeDataImporter.parseTable( "name\tclassification\nGCF_A\td__Bacteria;p__Pseudomonadota\n" );
        if ( MainFrame.gtdbKeyColumn( t, MainFrame.gtdbClassificationColumn( t ) ) != 0 ) {
            return note( "a plain name-keyed table should key on the name column" );
        }
        return true;
    }

    /** A tiny tree of three accession-named tips (GCF_A/B/C), no taxonomy -- the input to the GTDB import. */
    private static Phylogeny genomeTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "GCF_A" );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "GCF_B" );
        final PhylogenyNode c = new PhylogenyNode();
        c.setName( "GCF_C" );
        root.addAsChild( a );
        root.addAsChild( b );
        root.addAsChild( c );
        phy.setRoot( root );
        phy.setRooted( true );
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
        if ( ( n.getNodeData() == null ) || ( n.getNodeData().getProperties() == null ) ) {
            return null;
        }
        final List<Property> ps = n.getNodeData().getProperties().getPropertiesWithGivenRef( ref );
        return ( ( ps != null ) && !ps.isEmpty() ) ? ps.get( 0 ).getValue() : null;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ImportGtdbToolTest] " + msg );
        ok[ 0 ] = false;
    }

    private static boolean note( final String msg ) {
        System.out.println( "  [ImportGtdbToolTest] " + msg );
        return false;
    }

    private ImportGtdbToolTest() {
    }
}
