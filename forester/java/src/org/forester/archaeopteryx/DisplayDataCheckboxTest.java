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

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Integration test for the data-driven "Display Data" checkboxes: only checkboxes whose data is
 * actually present in the loaded tree are shown (Node Name is always shown). Loads a tree with
 * taxonomy scientific names but no sequences, asserts the right checkboxes are visible/hidden, then
 * adds a sequence and confirms the Seq Name checkbox appears after a refresh. Needs FlatLaf and a
 * display, so it runs standalone, not in the headless suite.
 */
public final class DisplayDataCheckboxTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DisplayDataCheckbox: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            final Phylogeny phy = taxonomyOnlyTree();
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "dd" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
                // present data -> visible
                check( ok, cp, Configuration.show_node_names, true, "Node Name (always)" );
                check( ok, cp, Configuration.show_taxonomy_scientific_names, true, "Taxonomy Scientific" );
                // absent data -> hidden
                check( ok, cp, Configuration.show_seq_names, false, "Seq Name (no sequences)" );
                check( ok, cp, Configuration.show_gene_names, false, "Gene Name (no sequences)" );
                check( ok, cp, Configuration.show_domain_architectures, false, "Domain Architectures" );
                check( ok, cp, Configuration.show_vector_data, false, "Vector Data" );
                check( ok, cp, Configuration.show_taxonomy_common_names, false, "Taxonomy Common" );
                // now add a sequence name to a leaf; after a refresh the Seq Name checkbox must appear
                final PhylogenyNode leaf = phy.getFirstExternalNode();
                final Sequence seq = new Sequence();
                seq.setName( "myseq" );
                leaf.getNodeData().setSequence( seq );
                cp.displayedPhylogenyMightHaveChanged( true );
                check( ok, cp, Configuration.show_seq_names, true, "Seq Name (after adding a sequence)" );
                ( (JFrame) mf[ 0 ] ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void check( final boolean[] ok,
                               final ControlPanel cp,
                               final int which,
                               final boolean expected_visible,
                               final String label ) {
        if ( cp.isDisplayDataCheckboxVisible( which ) != expected_visible ) {
            ok[ 0 ] = false;
            System.out.println( "  expected " + label + " checkbox visible=" + expected_visible
                    + " but was " + cp.isDisplayDataCheckboxVisible( which ) );
        }
    }

    /** Two leaves, each with a taxonomy scientific name and a node name, but no sequence/other data. */
    private static Phylogeny taxonomyOnlyTree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( final String sn : new String[] { "Homo sapiens", "Pan troglodytes" } ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( sn.replace( ' ', '_' ) );
            final Taxonomy tax = new Taxonomy();
            tax.setScientificName( sn );
            leaf.getNodeData().setTaxonomy( tax );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private DisplayDataCheckboxTest() {
    }
}
