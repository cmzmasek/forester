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

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.GraphicsEnvironment;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.util.ForesterUtil;

/**
 * Headful integration test for {@link NodeStyleDialog} + {@link TreePanel#applyNodeStyleEdit}: a single-node edit
 * builds a spec from the ticked controls, writes exactly those attributes, turns on "Use Visual Styles", and is
 * undoable (without writing provenance); a bulk edit styles every target node; and opening the dialog on an
 * already-styled node PRE-fills (pre-ticks) that attribute. A no-op on a headless box (needs FlatLaf via
 * {@code createInstance}).
 */
public final class NodeStyleDialogTest {

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree( 4 ) }, conf,
                                                                         "node style test" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    final PhylogenyNode[] leaves = tp.getPhylogeny().getExternalNodes()
                            .toArray( new PhylogenyNode[ 0 ] );

                    // ---- single node: build spec from ticked controls, apply, undo ----
                    final NodeStyleDialog d = new NodeStyleDialog( tp, Arrays.asList( leaves[ 0 ] ), true );
                    d.setApplyFontColorForTest( Color.RED );
                    d.setApplyShapeForTest( NodeShape.CIRCLE );
                    final NodeStyleEditor.Spec spec = d.buildSpec();
                    ck( ok, ( spec.fontColor() != null ) && spec.fontColor().equals( Color.RED ),
                        "spec carries the ticked font color" );
                    ck( ok, spec.shape() == NodeShape.CIRCLE, "spec carries the ticked shape" );
                    ck( ok, spec.fill() == null, "spec leaves the un-ticked fill null" );
                    d.applyForTest();
                    final NodeVisualData vis0 = leaves[ 0 ].getNodeData().getNodeVisualData();
                    ck( ok, ( vis0 != null ) && Color.RED.equals( vis0.getFontColor() ), "node got the font color" );
                    ck( ok, ( vis0 != null ) && ( vis0.getShape() == NodeShape.CIRCLE ), "node got the shape" );
                    ck( ok, tp.getControlPanel().isUseVisualStyles(), "\"Use Visual Styles\" was turned on" );
                    // a visual-style edit is undoable but does NOT write provenance to the description
                    ck( ok, ForesterUtil.isEmpty( tp.getPhylogeny().getDescription() ),
                        "no provenance sentence written for a visual-style edit" );
                    // undo restores the pre-edit (un-styled) tree
                    mf[ 0 ].undo();
                    final NodeVisualData after_undo = tp.getPhylogeny().getExternalNodes().get( 0 ).getNodeData()
                            .getNodeVisualData();
                    ck( ok, ( after_undo == null ) || ( after_undo.getFontColor() == null ),
                        "undo removed the style" );

                    // ---- bulk: the spec applies to every target node ----
                    final NodeStyleDialog b = new NodeStyleDialog( tp,
                                                                   Arrays.asList( leaves[ 1 ], leaves[ 2 ] ), false );
                    b.setApplyFontColorForTest( Color.BLUE );
                    b.applyForTest();
                    ck( ok, Color.BLUE.equals( leaves[ 1 ].getNodeData().getNodeVisualData().getFontColor() ),
                        "bulk styled node 1" );
                    ck( ok, Color.BLUE.equals( leaves[ 2 ].getNodeData().getNodeVisualData().getFontColor() ),
                        "bulk styled node 2" );

                    // ---- pre-fill: opening the dialog on an already-styled node pre-ticks that attribute ----
                    final NodeStyleDialog p = new NodeStyleDialog( tp, Arrays.asList( leaves[ 1 ] ), true );
                    ck( ok, ( p.buildSpec().fontColor() != null ),
                        "single-node dialog pre-fills the node's existing font color" );

                    // ---- pre-fill NEVER shrinks an out-of-spinner-range stored value (review #1) ----
                    final NodeVisualData big = new NodeVisualData();
                    big.setFontName( "Source Sans 3" );
                    big.setFontSize( 72 ); // valid in the model (<=127) but beyond the spinner's default max of 64
                    leaves[ 3 ].getNodeData().setNodeVisualData( big );
                    final NodeStyleDialog q = new NodeStyleDialog( tp, Arrays.asList( leaves[ 3 ] ), true );
                    ck( ok, ( q.buildSpec().fontSize() != null ) && ( q.buildSpec().fontSize().intValue() == 72 ),
                        "pre-fill keeps an out-of-range font size (72), does not clamp to 64" );

                    // ---- entry-point wiring (the methods open modal dialogs, untestable headlessly) ----
                    ck( ok, toolsItem( mf[ 0 ].getJMenuBar(), "Node Style for Selected Nodes" ) != null,
                        "Tools menu has the 'Node Style for Selected Nodes' item" );
                    ck( ok, comboHasItem( tp.getControlPanel(), "Node Style" ),
                        "the 'Click on Node to:' dropdown offers 'Node Style'" );

                    // ---- the label offset uses the ACTUAL per-node mark size so a big mark doesn't overlap ----
                    // (a prior edit already turned "Use Visual Styles" on, so the per-node size applies)
                    final PhylogenyNode big_node = new PhylogenyNode();
                    final NodeVisualData big_mark = new NodeVisualData();
                    big_mark.setSize( 30f );
                    big_node.getNodeData().setNodeVisualData( big_mark );
                    ck( ok, tp.effectiveNodeHalfBoxSize( big_node ) == 15,
                        "label offset uses the per-node mark size (30 -> 15)" );
                    ck( ok, tp.effectiveNodeHalfBoxSize( new PhylogenyNode() )
                            == ( tp.getOptions().getDefaultNodeShapeSize() / 2 ),
                        "a plain node uses the default half size" );

                    // ---- diamond shape: the dialog maps it both ways ----
                    final NodeStyleDialog dia = new NodeStyleDialog( tp, Arrays.asList( leaves[ 0 ] ), false );
                    dia.setApplyShapeForTest( NodeShape.DIAMOND );
                    ck( ok, dia.buildSpec().shape() == NodeShape.DIAMOND, "dialog builds a DIAMOND shape spec" );
                    final PhylogenyNode dia_node = new PhylogenyNode();
                    final NodeVisualData dia_vis = new NodeVisualData();
                    dia_vis.setShape( NodeShape.DIAMOND );
                    dia_node.getNodeData().setNodeVisualData( dia_vis );
                    ck( ok, new NodeStyleDialog( tp, Arrays.asList( dia_node ), true ).buildSpec().shape()
                            == NodeShape.DIAMOND, "single-node dialog pre-fills a diamond shape" );

                    mf[ 0 ].dispose();
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** A Tools-menu item whose text contains {@code contains}, or null. */
    private static JMenuItem toolsItem( final JMenuBar bar, final String contains ) {
        for ( int i = 0; i < bar.getMenuCount(); ++i ) {
            final JMenu m = bar.getMenu( i );
            if ( ( m == null ) || !"Tools".equals( m.getText() ) ) {
                continue;
            }
            for ( int j = 0; j < m.getItemCount(); ++j ) {
                final JMenuItem it = m.getItem( j );
                if ( ( it != null ) && ( it.getText() != null ) && it.getText().contains( contains ) ) {
                    return it;
                }
            }
        }
        return null;
    }

    /** Whether some JComboBox under {@code root} carries an item equal to {@code label} (the click-to dropdown). */
    private static boolean comboHasItem( final Container root, final String label ) {
        final List<JComboBox> combos = new ArrayList<>();
        collect( root, JComboBox.class, combos );
        for ( final JComboBox<?> cb : combos ) {
            for ( int i = 0; i < cb.getItemCount(); ++i ) {
                if ( label.equals( String.valueOf( cb.getItemAt( i ) ) ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    @SuppressWarnings( "unchecked" )
    private static <T> void collect( final Container c, final Class<T> type, final List<T> out ) {
        for ( final Component comp : c.getComponents() ) {
            if ( type.isInstance( comp ) ) {
                out.add( (T) comp );
            }
            if ( comp instanceof Container ) {
                collect( (Container) comp, type, out );
            }
        }
    }

    private static void ck( final boolean[] ok, final boolean cond, final String msg ) {
        if ( !cond ) {
            System.out.println( "  [NodeStyleDialogTest] " + msg );
            ok[ 0 ] = false;
        }
    }

    private static Phylogeny tree( final int n ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for ( int i = 0; i < n; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            root.addAsChild( leaf );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private NodeStyleDialogTest() {
    }
}
