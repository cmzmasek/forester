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
 * Integration test for node/branch selection on a real {@link MainFrame}/{@link TreePanel}, all in
 * "Select Node(s)" mode over found set 0. Covers: findBranch resolution (horizontal leg, vertical fork
 * connector, and the ROOT via its connector), the all-or-nothing subtree toggle, the deselect-to-empty count
 * label, the mouseClicked/shift-click wiring, the collapsed-clade skip, and the hand/arrow hover cursor. Also
 * the select/deselect HOVER PREVIEW: the direction (add vs remove) for a hovered branch (subtree) and a hovered
 * single node (leaf or internal -- Option A), that a collapsed triangle is never circled by the preview (only
 * the hand cursor + its own fill/outline) and a select-click on it toggles the whole clade's tips, the post-click
 * suppression (armed by a click, kept across clearHoverPreview/mouse-exit, released on
 * move-off or a tree change), clearing on mouse-exit and tree swap, that a node hidden under a collapsed
 * ancestor is neither found nor selectable, and that a collapsed triangle reflects its hidden tips' selection
 * state (none/partial/all). Guarded to a no-op on a headless box (needs FlatLaf via {@code createInstance}).
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

        // the ROOT is selectable via its own fork connector -> selects the whole tree's tips
        final Point root_conn = new Point( Math.round( root.getXcoord() ),
                                           Math.round( ( a.getYcoord() + root.getYcoord() ) / 2f ) );
        if ( tp.findBranch( root_conn.x, root_conn.y ) != root ) {
            fail( ok, "clicking the root's connector should resolve to the root" );
        }
        tp.setFoundNodes0( null );
        tp.selectSubtreeTips( root );
        if ( count( tp ) != 4 ) {
            fail( ok, "selecting the root should select all 4 tips, got " + count( tp ) );
        }
        tp.setFoundNodes0( null ); // reset for the assertions that follow

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
        // after the click the preview is suppressed on that branch (so it doesn't instantly flip to the "remove"
        // grey over the just-selected tips) until the pointer moves off it
        if ( tp.hoverNodeForTest() != null ) {
            fail( ok, "a branch click should suppress the hover preview on the just-clicked branch" );
        }
        tp.mouseMoved( move( tp, b_mid ) ); // still on B -> stays suppressed
        if ( tp.hoverNodeForTest() != null ) {
            fail( ok, "the preview should stay suppressed while the pointer remains on the just-clicked branch" );
        }
        tp.mouseMoved( move( tp, new Point( 3, 3 ) ) ); // off any branch -> clears the suppression
        tp.mouseMoved( move( tp, b_mid ) ); // back onto B -> preview resumes (B is selected -> a remove)
        if ( ( tp.hoverNodeForTest() != b ) || !tp.hoverWouldDeselectForTest() ) {
            fail( ok, "moving off and back onto the branch should resume the (remove) preview" );
        }

        // shift-click on a branch selects the subtree in ANY click-to mode (not just Select-Node(s))
        tp.setFoundNodes0( null );
        cp.setActionWhenNodeClicked( NodeClickAction.SHOW_DATA );
        tp.mouseClicked( click( tp, a_mid, InputEvent.SHIFT_DOWN_MASK ) );
        if ( ( count( tp ) != 2 ) || !allSelected( tp, a.getAllExternalDescendants() ) ) {
            fail( ok, "shift-click on a branch should select its subtree regardless of the click-to mode" );
        }

        // hover preview: over B's branch (Select-Node(s) mode) -> hand cursor, B is the preview node, and since
        // B is not yet selected the preview is an "add" (not a deselect)
        tp.setFoundNodes0( null );
        cp.setActionWhenNodeClicked( NodeClickAction.SELECT_NODES );
        tp.mouseMoved( move( tp, b_mid ) );
        if ( tp.getCursor().getType() != Cursor.HAND_CURSOR ) {
            fail( ok, "hovering a branch in Select-Node(s) mode should show the hand cursor" );
        }
        if ( tp.hoverNodeForTest() != b ) {
            fail( ok, "hovering B's branch should set B as the hover-preview node" );
        }
        if ( tp.hoverWouldDeselectForTest() ) {
            fail( ok, "hovering an unselected clade should preview an add, not a remove" );
        }
        render( tp ); // exercise paintHoverPreview (add path) -- must not throw

        // once B's subtree is fully selected, hovering it previews a REMOVE (direction-aware)
        tp.selectSubtreeTips( b );
        tp.mouseMoved( move( tp, b_mid ) );
        if ( !tp.hoverWouldDeselectForTest() ) {
            fail( ok, "hovering a fully-selected clade should preview a remove" );
        }
        render( tp ); // exercise paintHoverPreview (remove path)

        // moving off any branch reverts the cursor AND clears the preview; a mouse-exit clears it too
        tp.mouseMoved( move( tp, new Point( 3, 3 ) ) );
        if ( tp.getCursor().getType() != Cursor.DEFAULT_CURSOR ) {
            fail( ok, "hovering empty space should revert to the arrow cursor" );
        }
        if ( tp.hoverNodeForTest() != null ) {
            fail( ok, "moving off a branch should clear the hover preview" );
        }
        tp.mouseMoved( move( tp, b_mid ) );
        tp.dispatchEvent( new MouseEvent( tp, MouseEvent.MOUSE_EXITED, System.currentTimeMillis(), 0, -1, -1, 0,
                                          false ) ); // routes through the registered MouseListener.mouseExited
        if ( tp.hoverNodeForTest() != null ) {
            fail( ok, "a mouse-exit should clear the hover preview" );
        }

        // replacing the tree clears any lingering hover preview (navigation swaps the tree)
        tp.mouseMoved( move( tp, b_mid ) );
        if ( tp.hoverNodeForTest() == null ) {
            fail( ok, "test setup: hover should be set before the tree swap" );
        }
        tp.setTree( tp.getPhylogeny() );
        if ( tp.hoverNodeForTest() != null ) {
            fail( ok, "setTree should clear the hover preview" );
        }

        // Option A: hovering a single node previews just that node (one marker), not a subtree
        tp.setFoundNodes0( null );
        final PhylogenyNode t0 = a.getChildNode( 0 );
        final Point t0_pt = new Point( Math.round( t0.getXcoord() ), Math.round( t0.getYcoord() ) );
        tp.mouseMoved( move( tp, t0_pt ) );
        if ( ( tp.hoverNodeForTest() != t0 ) || tp.hoverIsSubtreeForTest() ) {
            fail( ok, "hovering a leaf should preview that single node (not a subtree)" );
        }
        if ( tp.hoverWouldDeselectForTest() ) {
            fail( ok, "hovering an unselected leaf should preview an add" );
        }
        tp.selectNode( t0 );
        tp.mouseMoved( move( tp, t0_pt ) );
        if ( !tp.hoverWouldDeselectForTest() ) {
            fail( ok, "hovering a selected leaf should preview a remove" );
        }
        // an internal node is also previewed as a single node (its own marker, not its tips)
        final Point a_pt = new Point( Math.round( a.getXcoord() ), Math.round( a.getYcoord() ) );
        tp.mouseMoved( move( tp, a_pt ) );
        if ( ( tp.hoverNodeForTest() != a ) || tp.hoverIsSubtreeForTest() ) {
            fail( ok, "hovering an internal node should preview that single node" );
        }
        render( tp ); // exercise the single-node preview draw path

        // clicking a node suppresses its preview until the pointer moves off it (like a branch click)
        tp.setFoundNodes0( null );
        tp.mouseClicked( click( tp, t0_pt, 0 ) ); // node path -> selects t0
        if ( count( tp ) != 1 ) {
            fail( ok, "clicking a leaf should select it, got " + count( tp ) );
        }
        if ( tp.hoverNodeForTest() != null ) {
            fail( ok, "a node click should suppress its hover preview" );
        }
        tp.mouseMoved( move( tp, t0_pt ) ); // still on t0 -> suppressed
        if ( tp.hoverNodeForTest() != null ) {
            fail( ok, "the preview should stay suppressed while the pointer remains on the just-clicked node" );
        }
        tp.mouseMoved( move( tp, new Point( 3, 3 ) ) ); // off -> clears the suppression
        tp.mouseMoved( move( tp, t0_pt ) ); // back onto t0 -> preview resumes
        if ( tp.hoverNodeForTest() != t0 ) {
            fail( ok, "moving off and back onto the node should resume the preview" );
        }
        tp.setFoundNodes0( null ); // reset for the collapse assertions

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
        render( tp );

        // suppression lifecycle: a click arms hover-suppression on its target; clearHoverPreview (what a
        // mouse-exit and the legend branch call) must NOT release it -- otherwise a spurious exit would let the
        // just-selected node instantly flip to the grey "remove" preview -- but a tree change MUST clear it.
        tp.setFoundNodes0( null );
        cp.setActionWhenNodeClicked( NodeClickAction.SELECT_NODES );
        final PhylogenyNode leaf = a.getChildNode( 0 );
        final Point leaf_pt = new Point( Math.round( leaf.getXcoord() ), Math.round( leaf.getYcoord() ) );
        tp.mouseClicked( click( tp, leaf_pt, 0 ) );
        if ( tp.clickSuppressedForTest() != leaf ) {
            fail( ok, "a node click should arm hover-suppression on the clicked node" );
        }
        tp.clearHoverPreview(); // MouseListener.mouseExited and the legend branch both call this
        if ( tp.clickSuppressedForTest() != leaf ) {
            fail( ok, "clearHoverPreview must NOT release the click suppression (spurious-exit guard)" );
        }
        tp.setTree( tp.getPhylogeny() ); // a structural change routes through setNodeInPreorderToNull
        if ( tp.clickSuppressedForTest() != null ) {
            fail( ok, "a tree change should clear the click suppression" );
        }

        // the hover preview never circles a collapsed triangle -- it is ugly and redundant (the triangle's own
        // fill/outline shows selection). So a collapsed sub-clade inside a hovered subtree contributes NO mark...
        render( tp );
        tp.setFoundNodes0( null );
        cp.setActionWhenNodeClicked( NodeClickAction.SELECT_NODES );
        a.setCollapse( true ); // collapse only A; hovering the root's branch previews the whole tree
        render( tp );
        final Point root_conn2 = new Point( Math.round( root.getXcoord() ),
                                            Math.round( ( a.getYcoord() + root.getYcoord() ) / 2f ) );
        if ( tp.findNode( root_conn2.x, root_conn2.y ) != null ) {
            fail( ok, "test setup: the root connector must be clear of node boxes" );
        }
        tp.mouseMoved( move( tp, root_conn2 ) );
        if ( ( tp.hoverNodeForTest() != root ) || !tp.hoverIsSubtreeForTest() ) {
            fail( ok, "hovering the root's branch should preview the whole subtree" );
        }
        final List<PhylogenyNode> marks = tp.hoverPreviewMarksForTest();
        if ( marks.contains( a ) ) {
            fail( ok, "a collapsed sub-clade must NOT be circled in the hover preview" );
        }
        if ( marks.size() != 2 ) { // only B's two visible tips (A's collapsed triangle is not circled)
            fail( ok, "the preview should mark only the visible tips, got " + marks.size() );
        }
        render( tp ); // exercise the draw path -- must not throw
        // ...and hovering the collapsed triangle ITSELF sets no preview (just the hand cursor)
        final Point a_tri = new Point( Math.round( a.getXcoord() ), Math.round( a.getYcoord() ) );
        if ( tp.findNode( a_tri.x, a_tri.y ) != a ) {
            fail( ok, "test setup: the collapsed clade-root should be findable at its own coords" );
        }
        tp.mouseMoved( move( tp, a_tri ) );
        if ( tp.hoverNodeForTest() != null ) {
            fail( ok, "hovering a collapsed triangle must not set a hover preview (no circle over it)" );
        }
        if ( tp.getCursor().getType() != Cursor.HAND_CURSOR ) {
            fail( ok, "hovering a collapsed triangle should still show the hand cursor (it is clickable)" );
        }
        render( tp );
        // clicking the collapsed triangle toggles its WHOLE clade's tips (like a branch click), so the fill responds
        tp.setFoundNodes0( null );
        tp.mouseClicked( click( tp, a_tri, 0 ) );
        if ( ( count( tp ) != 2 ) || !allSelected( tp, a.getAllExternalDescendants() ) ) {
            fail( ok, "clicking a collapsed triangle should select its clade's 2 tips, got " + count( tp ) );
        }
        final int[] clicked = tp.collapsedCladeFoundCounts( a );
        if ( ( clicked[ 2 ] != 2 ) || ( clicked[ 3 ] != 2 ) ) {
            fail( ok, "clicking a collapsed triangle should fill it (all tips found), got " + clicked[ 2 ] + "/"
                    + clicked[ 3 ] );
        }
        tp.mouseClicked( click( tp, a_tri, 0 ) ); // click again -> deselect the whole clade
        if ( count( tp ) != 0 ) {
            fail( ok, "clicking a fully-selected collapsed triangle should deselect its tips, got " + count( tp ) );
        }
        render( tp );
        a.setCollapse( false );

        // collapsed-clade selection FEEDBACK: the triangle reflects its hidden tips' selection -- all tips
        // selected -> found fill; some -> found outline; none -> neither. collapsedCladeFoundCounts drives it.
        tp.setFoundNodes0( null );
        a.setCollapse( true );
        render( tp );
        final int[] none = tp.collapsedCladeFoundCounts( a );
        if ( ( none[ 2 ] != 0 ) || ( none[ 3 ] != 2 ) ) {
            fail( ok, "an unselected collapsed clade should report 0 of 2 tips, got " + none[ 2 ] + "/" + none[ 3 ] );
        }
        tp.selectNode( a.getChildNode( 0 ) ); // select one of A's two (hidden) tips -> partial
        final int[] partial = tp.collapsedCladeFoundCounts( a );
        if ( ( partial[ 0 ] != 1 ) || ( partial[ 2 ] != 1 ) || ( partial[ 3 ] != 2 ) ) {
            fail( ok, "a partially-selected collapsed clade should report 1 of 2, got " + partial[ 2 ] + "/"
                    + partial[ 3 ] );
        }
        render( tp ); // partial -> found-outline draw path, must not throw
        tp.selectSubtreeTips( a ); // fill A up to all tips selected
        final int[] all = tp.collapsedCladeFoundCounts( a );
        if ( ( all[ 0 ] != 2 ) || ( all[ 2 ] != 2 ) || ( all[ 3 ] != 2 ) ) {
            fail( ok, "a fully-selected collapsed clade should report 2 of 2, got " + all[ 2 ] + "/" + all[ 3 ] );
        }
        render( tp ); // all -> found-fill draw path, must not throw
        a.setCollapse( false );
        tp.setFoundNodes0( null );

        // a node hidden under a collapsed ANCESTOR (its own parent NOT the collapse-root -- the case the direct
        // getParent().isCollapse() guard missed) must not be found at its stale coords, so it gets no phantom
        // hover marker and a click there cannot select it. The collapsed clade-root itself stays selectable.
        final Phylogeny deep = deepTree();
        tp.setTree( deep );
        render( tp );
        final PhylogenyNode x = deep.getRoot().getChildNode( 0 );      // the clade to collapse
        final PhylogenyNode gchild = x.getChildNode( 0 ).getChildNode( 0 ); // a grandchild of X (parent is Y, not X)
        x.setCollapse( true );
        render( tp );
        final Point gchild_pt = new Point( Math.round( gchild.getXcoord() ), Math.round( gchild.getYcoord() ) );
        if ( tp.findNode( gchild_pt.x, gchild_pt.y ) != null ) {
            fail( ok, "a node hidden under a collapsed ancestor must not be found at its stale coords" );
        }
        tp.setFoundNodes0( null );
        cp.setActionWhenNodeClicked( NodeClickAction.SELECT_NODES );
        tp.mouseClicked( click( tp, gchild_pt, 0 ) );
        if ( count( tp ) != 0 ) {
            fail( ok, "clicking at a hidden node's stale coords must not select it, got " + count( tp ) );
        }
        final Point x_pt = new Point( Math.round( x.getXcoord() ), Math.round( x.getYcoord() ) );
        if ( tp.findNode( x_pt.x, x_pt.y ) != x ) {
            fail( ok, "the collapsed clade-root itself must still be findable (it is drawn as the triangle)" );
        }
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

    /** deepRoot -> { X -> (Y -> (g0,g1), w), z } -- a depth-3 clade so collapsing X hides grandchildren whose
     *  own parent (Y) is NOT the collapse-root. */
    private static Phylogeny deepTree() {
        final PhylogenyNode r = new PhylogenyNode();
        final PhylogenyNode x = new PhylogenyNode();
        final PhylogenyNode y = new PhylogenyNode();
        y.addAsChild( leaf( "g0" ) );
        y.addAsChild( leaf( "g1" ) );
        x.addAsChild( y );
        x.addAsChild( leaf( "w" ) );
        r.addAsChild( x );
        r.addAsChild( leaf( "z" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( r );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }
}
