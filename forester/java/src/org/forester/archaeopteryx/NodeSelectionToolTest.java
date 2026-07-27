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
import java.awt.Cursor;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Point;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.List;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.ControlPanel.NodeClickAction;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Integration test for branch-click clade selection on a real {@link MainFrame}/{@link TreePanel}: a click on
 * a branch (not a node) selects/deselects all external tips of that subtree, reusing the "Select Node(s)"
 * selection (found set 0). Covers findBranch resolution (horizontal leg + vertical connector), the
 * all-or-nothing subtree toggle, the deselect-to-empty count label, the mouseClicked/shift-click wiring, the
 * hover-cursor affordance (hand over a branch, arrow off it), and the collapsed-clade skip. Guarded to a no-op
 * on a headless box (needs FlatLaf via {@code createInstance}).
 */
public final class NodeSelectionToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeSelectionTool: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { balancedTree() }, conf,
                                                                        "nodesel" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    exercise( frame, ok );
                }
                finally {
                    ( (JFrame) frame ).dispose(); // never leak the window, even on an assertion throw
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void exercise( final MainFrame frame, final boolean[] ok ) {
        final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
        final ControlPanel cp = tp.getControlPanel();
        tp.getOptions().setShowOverview( false ); // keep the overview rectangle out of the click path
        render( tp );

        // structure: root -> {A -> (t0,t1), B -> (t2,t3)}
        final PhylogenyNode root = tp.getPhylogeny().getRoot();
        final PhylogenyNode a = root.getChildNode( 0 );
        final PhylogenyNode b = root.getChildNode( 1 );

        // findBranch resolves a horizontal-branch midpoint to that child; empty space -> null
        final Point a_mid = branchMid( a );
        if ( tp.findNode( a_mid.x, a_mid.y ) != null ) {
            fail( ok, "test setup: A's branch midpoint must be clear of node boxes" );
        }
        if ( tp.findBranch( a_mid.x, a_mid.y ) != a ) {
            fail( ok, "findBranch at A's branch midpoint should return A" );
        }
        if ( tp.findBranch( 2, 2 ) != null ) {
            fail( ok, "findBranch in empty space should be null" );
        }

        // the vertical fork connector (a point between A's children, off both legs) also resolves to A
        final Point a_conn = new Point( Math.round( a.getXcoord() ),
                                        Math.round( ( a.getChildNode( 0 ).getYcoord() + a.getYcoord() ) / 2f ) );
        if ( tp.findBranch( a_conn.x, a_conn.y ) != a ) {
            fail( ok, "clicking A's vertical connector should resolve to A" );
        }

        // selectSubtreeTips: all-or-nothing toggle
        tp.selectSubtreeTips( a );
        if ( ( count( tp ) != 2 ) || !allSelected( tp, a.getAllExternalDescendants() ) ) {
            fail( ok, "selecting A's subtree should select exactly its 2 tips, got " + count( tp ) );
        }
        tp.selectSubtreeTips( a );
        if ( count( tp ) != 0 ) {
            fail( ok, "a second selectSubtreeTips(A) should deselect all of A's tips, got " + count( tp ) );
        }

        // deselecting the last single node must reset the "Found:" count label to 0
        tp.selectNode( a.getChildNode( 0 ) );
        tp.selectNode( a.getChildNode( 0 ) ); // toggle off -> empty
        if ( !"Found: 0".equals( cp.getSearchFoundCountsLabel0().getText() ) ) {
            fail( ok, "deselecting the last node should reset the label to 'Found: 0', got '"
                    + cp.getSearchFoundCountsLabel0().getText() + "'" );
        }

        // a partially-selected subtree fills to all, then a further toggle clears
        tp.selectNode( a.getChildNode( 0 ) );
        tp.selectSubtreeTips( a );
        if ( count( tp ) != 2 ) {
            fail( ok, "toggling a partially-selected subtree should select all its tips, got " + count( tp ) );
        }
        tp.selectSubtreeTips( a );

        // mouseClicked on B's branch in Select-Node(s) mode selects B's tips
        cp.setActionWhenNodeClicked( NodeClickAction.SELECT_NODES );
        final Point b_mid = branchMid( b );
        if ( tp.findNode( b_mid.x, b_mid.y ) != null ) {
            fail( ok, "test setup: B's branch midpoint must be clear of node boxes" );
        }
        tp.mouseClicked( click( tp, b_mid, 0 ) );
        if ( ( count( tp ) != 2 ) || !allSelected( tp, b.getAllExternalDescendants() ) ) {
            fail( ok, "a branch click should select B's 2 tips, got " + count( tp ) );
        }

        // shift-click on a branch selects the subtree in ANY click-to mode (not just Select-Node(s))
        tp.setFoundNodes0( null );
        cp.setActionWhenNodeClicked( NodeClickAction.SHOW_DATA );
        tp.mouseClicked( click( tp, a_mid, InputEvent.SHIFT_DOWN_MASK ) );
        if ( ( count( tp ) != 2 ) || !allSelected( tp, a.getAllExternalDescendants() ) ) {
            fail( ok, "shift-click on a branch should select its subtree regardless of the click-to mode" );
        }

        // hover: over a branch (Select-Node(s) mode) shows the hand cursor; off any branch reverts to the arrow
        cp.setActionWhenNodeClicked( NodeClickAction.SELECT_NODES );
        tp.mouseMoved( move( tp, b_mid ) );
        if ( tp.getCursor().getType() != Cursor.HAND_CURSOR ) {
            fail( ok, "hovering a branch in Select-Node(s) mode should show the hand cursor" );
        }
        tp.mouseMoved( move( tp, new Point( 3, 3 ) ) );
        if ( tp.getCursor().getType() != Cursor.DEFAULT_CURSOR ) {
            fail( ok, "hovering empty space should revert to the arrow cursor" );
        }

        // a collapsed clade's own branch, and any branch hidden under it, are not selectable
        a.setCollapse( true );
        if ( tp.findBranch( a_mid.x, a_mid.y ) != null ) {
            fail( ok, "a collapsed clade's branch should not be selectable" );
        }
        final Point t0_mid = branchMid( a.getChildNode( 0 ) );
        if ( tp.findBranch( t0_mid.x, t0_mid.y ) != null ) {
            fail( ok, "a branch hidden under a collapsed ancestor should not be selectable" );
        }
        a.setCollapse( false );
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [NodeSelectionToolTest] " + msg );
        ok[ 0 ] = false;
    }

    private static MouseEvent click( final TreePanel tp, final Point p, final int modifiers ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(), modifiers, p.x, p.y, 1,
                               false, MouseEvent.BUTTON1 );
    }

    private static MouseEvent move( final TreePanel tp, final Point p ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_MOVED, System.currentTimeMillis(), 0, p.x, p.y, 0, false );
    }

    /** The midpoint of {@code n}'s incoming horizontal branch (from its parent's x to its x, at its y). */
    private static Point branchMid( final PhylogenyNode n ) {
        final PhylogenyNode p = n.getParent();
        return new Point( Math.round( ( p.getXcoord() + n.getXcoord() ) / 2f ), Math.round( n.getYcoord() ) );
    }

    private static int count( final TreePanel tp ) {
        return ( tp.getFoundNodes0() == null ) ? 0 : tp.getFoundNodes0().size();
    }

    private static boolean allSelected( final TreePanel tp, final List<PhylogenyNode> nodes ) {
        final Set<Long> f = tp.getFoundNodes0();
        if ( f == null ) {
            return nodes.isEmpty();
        }
        for( final PhylogenyNode n : nodes ) {
            if ( !f.contains( n.getId() ) ) {
                return false;
            }
        }
        return true;
    }

    private static void render( final TreePanel tp ) {
        final int w = 700, h = 500;
        tp.setSize( w, h );
        tp.validate();
        tp.doLayout();
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.printAll( g );
        g.dispose();
    }

    /** root -> {A -> (t0,t1), B -> (t2,t3)} */
    private static Phylogeny balancedTree() {
        final PhylogenyNode root = new PhylogenyNode();
        int t = 0;
        for( int k = 0; k < 2; ++k ) {
            final PhylogenyNode internal = new PhylogenyNode();
            root.addAsChild( internal );
            for( int j = 0; j < 2; ++j ) {
                final PhylogenyNode leaf = new PhylogenyNode();
                leaf.setName( "t" + ( t++ ) );
                internal.addAsChild( leaf );
            }
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }
}
