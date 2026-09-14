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
import java.io.File;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * GUI integration test of the re-rooting rules ({@link Rerooting}) in a real window: the controls that re-root are
 * greyed out -- with the reason as their tooltip -- for a tree marked rerootable="false" and for a time tree, and follow
 * the tab; every re-root path refuses such a tree without touching it; the node-data warning is asked only when a
 * re-root changes an annotated node's clade, and Cancel leaves the tree as it was; and a tree declared unrooted, shown
 * unrooted, drops the root-dependent search fields and the internal branch-length field. Guarded to a no-op on a
 * headless box.
 */
public final class RootingGuardsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RootingGuards: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            return lockedTree() && timeTreeFollowsTheTab() && nodeDataWarning() && rootFreeView();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean lockedTree() throws Exception {
        final Phylogeny phy = nhx( "((A:1,B:2)x:1,(C:1,D:3):1)" );
        phy.setRerootable( false );
        final MainFrame mf = open( phy );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
            // greyed out as soon as the tree is shown, with the reason as tooltip
            if ( !greyedWith( mf, Rerooting.NOT_REROOTABLE ) ) {
                ok[ 0 ] = TestFail.here( "the re-rooting menu items must be greyed out with the reason" );
            }
            final ControlPanel cp = tp.getControlPanel();
            final int reroot_index = cp.rerootClickToIndexForTest();
            if ( cp.clickToEntryForTest( reroot_index ).isEnabled() ) {
                ok[ 0 ] = TestFail.here( "the Root/Reroot click option must be drawn greyed out" );
            }
            final JComboBox<String> combo = cp.clickToComboForTest();
            final int chosen = combo.getSelectedIndex();
            combo.setSelectedIndex( reroot_index );
            if ( combo.getSelectedIndex() != chosen ) {
                ok[ 0 ] = TestFail.here( "choosing the greyed-out Root/Reroot option must fall back to the previous one" );
            }
            // every re-root path refuses, leaves the tree as it was, and records no undo step
            final List<String> said = new ArrayList<>();
            tp.setRerootDialogsForTest( said::add, w -> {
                said.add( "WARNING " + w );
                return true;
            } );
            final String before = tp.getPhylogeny().toNewHampshire();
            tp.madRoot();
            tp.midpointRoot();
            tp.reRoot( tp.getPhylogeny().getNode( "C" ) );
            if ( !said.equals( List.of( Rerooting.NOT_REROOTABLE, Rerooting.NOT_REROOTABLE, Rerooting.NOT_REROOTABLE ) )
                    || !before.equals( tp.getPhylogeny().toNewHampshire() ) || tp.canUndo() ) {
                ok[ 0 ] = TestFail.here( "refused: " + said );
            }
            // a subtree view (which shares the full tree's nodes) keeps the flag
            tp.subTree( tp.getPhylogeny().getNode( "x" ) );
            if ( tp.rerootRefusal() != Rerooting.NOT_REROOTABLE ) {
                ok[ 0 ] = TestFail.here( "a subtree view of a not-re-rootable tree must not be re-rootable" );
            }
            ( (JFrame) mf ).dispose();
        } );
        return ok[ 0 ];
    }

    private static boolean timeTreeFollowsTheTab() throws Exception {
        final String nh = "((A:1,B:1)x:1,(C:1,D:1)y:1)";
        final Phylogeny tip_dated = nhx( nh ); // a divergence tree with sampling dates: may be re-rooted
        for( final PhylogenyNode t : tip_dated.getExternalNodes() ) {
            date( t );
        }
        final Phylogeny time_tree = nhx( nh );
        date( time_tree.getRoot() );
        date( time_tree.getNode( "x" ) );
        date( time_tree.getNode( "y" ) );
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { tip_dated, time_tree }, new Configuration(), "rooting" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            // the LAST tab (the time tree) is selected on open
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            if ( ( tp.getPhylogeny() != time_tree ) || !greyedWith( mf[ 0 ], Rerooting.TIME_TREE ) ) {
                ok[ 0 ] = TestFail.here( "a time tree greys the re-rooting controls out" );
            }
            final List<String> said = new ArrayList<>();
            tp.setRerootDialogsForTest( said::add, w -> true );
            tp.madRoot();
            if ( !said.equals( List.of( Rerooting.TIME_TREE ) ) ) {
                ok[ 0 ] = TestFail.here( "a time tree refuses MAD-Root: " + said );
            }
            // switching to the tip-dated tree's tab re-enables them, with their own tooltips back
            mf[ 0 ].getMainPanel().getTabbedPane().setSelectedIndex( 0 );
            final TreePanel tip_tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            if ( ( tip_tp.getPhylogeny() != tip_dated ) || !mf[ 0 ]._mad_root_item.isEnabled()
                    || !mf[ 0 ]._midpoint_root_item.isEnabled()
                    || Rerooting.TIME_TREE.equals( mf[ 0 ]._mad_root_item.getToolTipText() )
                    || !tip_tp.getControlPanel().clickToEntryForTest( tip_tp.getControlPanel().rerootClickToIndexForTest() )
                            .isEnabled() ) {
                ok[ 0 ] = TestFail.here( "tip dates alone keep the tree re-rootable, and the controls follow the tab" );
            }
            ( (JFrame) mf[ 0 ] ).dispose();
        } );
        return ok[ 0 ];
    }

    private static boolean nodeDataWarning() throws Exception {
        // ab, abc and de carry data (names); midpoint-rooting moves the root into C's long branch, which changes
        // abc's clade (1 of the 3)
        final String nh = "(((A:1,B:1)ab:1,C:6)abc:1,(D:1,E:1)de:1)";
        final MainFrame mf = open( nhx( nh ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
            final List<String> asked = new ArrayList<>();
            final boolean[] answer = { false };
            tp.setRerootDialogsForTest( m -> asked.add( "REFUSED " + m ), w -> {
                asked.add( w );
                return answer[ 0 ];
            } );
            final String before = tp.getPhylogeny().toNewHampshire();
            tp.midpointRoot(); // Cancel
            if ( !asked.equals( List.of( Rerooting.dataWarning( 3, 1 ) ) ) || !before.equals( tp.getPhylogeny().toNewHampshire() )
                    || tp.canUndo() ) {
                ok[ 0 ] = TestFail.here( "Cancel must leave the tree and the undo history untouched: " + asked );
            }
            asked.clear();
            answer[ 0 ] = true;
            tp.midpointRoot(); // Re-root
            if ( ( asked.size() != 1 ) || before.equals( tp.getPhylogeny().toNewHampshire() ) || !tp.canUndo() ) {
                ok[ 0 ] = TestFail.here( "Re-root must re-root, with an undo step" );
            }
            ( (JFrame) mf ).dispose();
        } );
        // a manual re-root on the other side of a two-child root changes no clade: no question asked
        final MainFrame mf2 = open( nhx( nh ) );
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf2.getMainPanel().getCurrentTreePanel();
            final List<String> asked = new ArrayList<>();
            tp.setRerootDialogsForTest( asked::add, w -> {
                asked.add( w );
                return false;
            } );
            tp.reRoot( tp.getPhylogeny().getNode( "de" ) );
            if ( !asked.isEmpty() || !tp.canUndo() ) {
                ok[ 0 ] = TestFail.here( "a re-root that changes no annotated clade must not ask: " + asked );
            }
            tp.reRoot( tp.getPhylogeny().getNode( "C" ) );
            if ( asked.size() != 1 ) {
                ok[ 0 ] = TestFail.here( "a manual re-root that changes an annotated clade must ask: " + asked );
            }
            ( (JFrame) mf2 ).dispose();
        } );
        return ok[ 0 ];
    }

    private static boolean rootFreeView() throws Exception {
        final File xml = File.createTempFile( "unrooted", ".xml" );
        xml.deleteOnExit();
        java.nio.file.Files.writeString( xml.toPath(),
                "<?xml version=\"1.0\" encoding=\"UTF-8\"?><phyloxml xmlns=\"http://www.phyloxml.org\">"
                        + "<phylogeny rooted=\"false\"><clade>"
                        + "<clade><name>x</name><branch_length>0.3</branch_length>"
                        + "<clade><name>A</name><branch_length>1</branch_length></clade>"
                        + "<clade><name>B</name><branch_length>1</branch_length></clade></clade>"
                        + "<clade><name>C</name><branch_length>2</branch_length></clade>"
                        + "<clade><name>D</name><branch_length>2</branch_length></clade>"
                        + "</clade></phylogeny></phyloxml>" );
        final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( xml, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
        final MainFrame mf = open( phy );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
            final ControlPanel cp = tp.getControlPanel();
            final PhylogenyNode x = tp.getPhylogeny().getNode( "x" );
            final String depth = "Structure: Depth from Root (edges)";
            cp.rebuildSearchFields( true );
            if ( tp.hidesRootDependentValues() || !cp.searchFieldLabelsForTest( true ).contains( depth )
                    || ( new NodeDataForm( x, tp, NodeDataForm.Mode.EDIT ).fieldForTest( NodeDataDraft.BRANCH_LENGTH ) == null ) ) {
                ok[ 0 ] = TestFail.here( "the rectangular layout keeps everything" );
            }
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
            cp.rebuildSearchFields( false ); // the same tree: rebuilt because the root-free mode changed
            if ( !tp.hidesRootDependentValues() || cp.searchFieldLabelsForTest( true ).contains( depth )
                    || cp.searchFieldLabelsForTest( true ).contains( "Structure: Clade Size (tips)" ) ) {
                ok[ 0 ] = TestFail.here( "unrooted + declared unrooted drops the root-dependent search fields: "
                        + cp.searchFieldLabelsForTest( true ) );
            }
            final NodeDataForm internal_form = new NodeDataForm( x, tp, NodeDataForm.Mode.EDIT );
            if ( ( internal_form.fieldForTest( NodeDataDraft.BRANCH_LENGTH ) != null )
                    || !NodeDataDraft.from( x ).branchLength.equals( internal_form.collect().branchLength )
                    || internal_form.collect().branchLength.isEmpty() ) {
                ok[ 0 ] = TestFail.here( "an internal node's branch length is not offered, and a write keeps it" );
            }
            if ( new NodeDataForm( tp.getPhylogeny().getNode( "A" ), tp, NodeDataForm.Mode.EDIT )
                    .fieldForTest( NodeDataDraft.BRANCH_LENGTH ) == null ) {
                ok[ 0 ] = TestFail.here( "a tip keeps its (single) branch length" );
            }
            ( (JFrame) mf ).dispose();
        } );
        return ok[ 0 ];
    }

    private static boolean greyedWith( final MainFrame mf, final String why ) {
        for( final JMenuItem item : new JMenuItem[] { mf._mad_root_item, mf._midpoint_root_item, mf._gsdir_item,
                mf._gsdir_taxonomy_item } ) {
            if ( ( item == null ) || item.isEnabled() || !why.equals( item.getToolTipText() ) ) {
                return false;
            }
        }
        return true;
    }

    private static MainFrame open( final Phylogeny phy ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy },
                                                                                          new Configuration(),
                                                                                          "rooting" ) );
        return mf[ 0 ];
    }

    private static Phylogeny nhx( final String nh ) throws Exception {
        return ParserBasedPhylogenyFactory.getInstance().create( nh, new NHXParser() )[ 0 ];
    }

    private static void date( final PhylogenyNode n ) {
        n.getNodeData().setDate( new Date( "", BigDecimal.ONE, null, null, "mya" ) );
    }

    private RootingGuardsTest() {
    }
}
