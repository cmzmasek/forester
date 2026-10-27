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
import java.awt.Cursor;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Point;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.SequenceAccessionTools;

/**
 * The node menu (0.11.191, Christian 2026-10-26: "No more old-skool 'Click on Node to:' menu, no more 'status' of
 * click behavior"): a click on a node -- left, right or ctrl -- opens ONE menu for that node, titled with its
 * name; each entry acts once and nothing is remembered. A shift-click selects instead. Pinned here by their
 * words: the entries and their groups in order; which are greyed on a tip, on the root and in a subtree view
 * (Root/Reroot with the reason as its tooltip); the editing entries only on an editable tree; the sequence link
 * named by its accession; the delete entry's wording and warning colour; the title's fallbacks and its cut. And
 * the gestures: Shift held turns the focus glow into the selection preview (branches become targets), with or
 * without a pointer move; the clipboard highlight shows as soon as something is held and fades at the next
 * edit of its own tab (there being no mode to leave). The unrooted and circular layouts are asked too, and the
 * subtree entry at subTree's nesting cap.
 * Guarded to a no-op on a headless box.
 */
public final class NodeMenuTest {

    private static final int W = 640;
    private static final int H = 420;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeMenu: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI test; needs a display toolkit
        }
        try {
            return pureOk() && menuOk() && subtreeCapOk();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** The pure parts: the title (name, else taxonomy, else "Internal node" / "Node"; cut to 34 code points with
     *  an ellipsis), the cut itself, and the delete dialog's choices by node. */
    private static boolean pureOk() {
        boolean ok = true;
        final PhylogenyNode named = new PhylogenyNode();
        named.setName( " BRCA1_HUMAN " );
        ok &= ck( "BRCA1_HUMAN".equals( TreePanelUtil.nodeMenuTitle( named ) ), "the title is the node's name, trimmed" );
        final PhylogenyNode tax = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( "Homo sapiens" );
        tax.getNodeData().setTaxonomy( t );
        ok &= ck( "Homo sapiens".equals( TreePanelUtil.nodeMenuTitle( tax ) ), "a nameless node is titled by its taxonomy" );
        final PhylogenyNode internal = new PhylogenyNode();
        internal.addAsChild( new PhylogenyNode() );
        ok &= ck( "Internal node".equals( TreePanelUtil.nodeMenuTitle( internal ) ), "a nameless internal node: 'Internal node'" );
        ok &= ck( "Node".equals( TreePanelUtil.nodeMenuTitle( new PhylogenyNode() ) ), "a nameless tip: 'Node'" );
        final PhylogenyNode lng = new PhylogenyNode();
        lng.setName( "0123456789012345678901234567890123456789" ); // 40 characters
        final String cut = TreePanelUtil.nodeMenuTitle( lng );
        ok &= ck( ( cut.length() == 34 ) && cut.startsWith( "012345678901234567890123456789012" ) && cut.endsWith( "…" ),
                  "a long name is cut to 33 characters plus an ellipsis, got '" + cut + "'" );
        final PhylogenyNode exact = new PhylogenyNode();
        exact.setName( "0123456789012345678901234567890123" ); // 34: not cut
        ok &= ck( exact.getName().equals( TreePanelUtil.nodeMenuTitle( exact ) ), "a 34-character name is not cut" );
        // the cut counts code points: an emoji at the end is never halved
        final String ascii33 = "012345678901234567890123456789012";
        final PhylogenyNode emoji34 = new PhylogenyNode();
        emoji34.setName( ascii33 + "\uD83D\uDE00" ); // 34 code points in 35 chars: not cut
        ok &= ck( emoji34.getName().equals( TreePanelUtil.nodeMenuTitle( emoji34 ) ),
                  "34 code points are not cut, whatever their char count" );
        final PhylogenyNode emoji35 = new PhylogenyNode();
        emoji35.setName( ascii33 + "\uD83D\uDE00x" ); // 35 code points: cut after 33, before the emoji
        ok &= ck( ( ascii33 + "\u2026" ).equals( TreePanelUtil.nodeMenuTitle( emoji35 ) ),
                  "the cut never leaves half a surrogate pair, got '" + TreePanelUtil.nodeMenuTitle( emoji35 ) + "'" );
        ok &= ck( "abcdefghijklmnopqrs\u2026".equals( TreePanelUtil.ellipsize( "abcdefghijklmnopqrstu", 20 ) )
                  && "abcdefghijklmnopqrst".equals( TreePanelUtil.ellipsize( "abcdefghijklmnopqrst", 20 ) )
                  && ( TreePanelUtil.ellipsize( null, 20 ) == null ),
                  "ellipsize: 20 stays, 21 becomes 19 + an ellipsis, null stays null" );
        // the delete dialog's choices: a tip is deleted or not, a root only itself, an internal node asks which
        ok &= ck( Arrays.equals( new String[] { "Delete", "Cancel" }, TreePanel.deleteOptionsFor( new PhylogenyNode() ) ),
                  "a tip: delete or not" );
        ok &= ck( Arrays.equals( new String[] { "Node only", "Cancel" }, TreePanel.deleteOptionsFor( internal ) ),
                  "a root: itself only, never the whole tree" );
        final PhylogenyNode mid = internal.getChildNode( 0 );
        mid.addAsChild( new PhylogenyNode() );
        ok &= ck( Arrays.equals( new String[] { "Node only", "Entire subtree", "Cancel" }, TreePanel.deleteOptionsFor( mid ) ),
                  "an internal node: node only or its entire subtree" );
        return ok;
    }

    private static boolean menuOk() throws Exception {
        // ((tipA, B[P38398])ab, (Homo sapiens, D)cd)
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode ab = child( root, "ab" );
        final PhylogenyNode cd = child( root, "cd" );
        final PhylogenyNode a = child( ab, "tipA" );
        final PhylogenyNode b = child( ab, "B" );
        final PhylogenyNode c = child( cd, null );
        final PhylogenyNode d = child( cd, "D" );
        final Sequence seq = new Sequence();
        seq.setName( "P38398" ); // named by its accession: the link says it once
        seq.setAccession( new Accession( "P38398", "UniProtKB" ) );
        b.getNodeData().addSequence( seq );
        final Taxonomy tax = new Taxonomy();
        tax.setScientificName( "Homo sapiens" );
        c.getNodeData().setTaxonomy( tax );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        if ( SequenceAccessionTools.obtainAccessorFromSequence( seq ) == null ) {
            return fail( "precondition: the sequence's accession must be recognised" );
        }

        // a second tab, for the clipboard highlight's ownership
        final PhylogenyNode root2 = new PhylogenyNode();
        child( root2, "X" );
        child( root2, "Y" );
        final Phylogeny phy2 = new Phylogeny();
        phy2.setRoot( root2 );
        phy2.setRooted( true );
        phy2.externalNodesHaveChanged();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy, phy2 }, new Configuration(), "node menu" ) );
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> {
                try {
                    ok[ 0 ] = exercise( mf[ 0 ], ab, cd, a, b, c, d );
                }
                catch ( final Throwable e ) {
                    e.printStackTrace();
                    ok[ 0 ] = false;
                }
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ].dispose() );
        }
        return ok[ 0 ];
    }

    private static boolean exercise( final MainFrame mf, final PhylogenyNode ab, final PhylogenyNode cd,
                                     final PhylogenyNode a, final PhylogenyNode b, final PhylogenyNode c,
                                     final PhylogenyNode d ) {
        boolean ok = true;
        final javax.swing.JTabbedPane tabs = mf.getMainPanel().getTabbedPane();
        tabs.setSelectedIndex( 0 ); // the LAST tab is selected on open; the menu tree is the first
        final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
        final PhylogenyNode root = tp.getPhylogeny().getRoot();
        ok &= ck( root == ab.getParent(), "precondition: the first tab shows the menu tree" );
        tp.getOptions().setShowOverview( false ); // keep the overview rectangle out of the click path
        tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        render( tp );

        // ---- (1) the entries, their groups and their order, pinned by their words --------------------------
        final List<String> expected_editable = Arrays.asList( "|", "Display Node Data", "Edit Node Data",
                "Edit Node Style", "|", "Switch to Subtree", "Collapse/Uncollapse", "Uncollapse Subtree",
                "Swap Descendants", "Ladderize Subtree", "Colorize Subtree", "Root/Reroot", "|",
                "Select/Deselect Node", "Select/Deselect Tips", "|", "Open Sequence DB", "Open Taxonomy DB", "|",
                "Cut Subtree", "Copy Subtree", "Paste Subtree", "Add New Node", "|", "Delete Node/Subtree…" );
        ok &= ck( tp.getOptions().isEditable(), "precondition: the tree is editable by default" );
        ok &= ck( expected_editable.equals( words( tp.buildNodeMenu( ab ) ) ),
                  "the menu for an internal node, editable tree: " + words( tp.buildNodeMenu( ab ) ) );
        tp.getOptions().setEditableForTest( false );
        final List<String> expected_read_only = Arrays.asList( "|", "Display Node Data", "Edit Node Style", "|",
                "Switch to Subtree", "Collapse/Uncollapse", "Uncollapse Subtree", "Swap Descendants",
                "Ladderize Subtree", "Colorize Subtree", "Root/Reroot", "|", "Select/Deselect Node",
                "Select/Deselect Tips", "|", "Open Sequence DB", "Open Taxonomy DB" );
        ok &= ck( expected_read_only.equals( words( tp.buildNodeMenu( ab ) ) ),
                  "the editing entries exist only on an editable tree: " + words( tp.buildNodeMenu( ab ) ) );
        tp.getOptions().setEditableForTest( true );

        // ---- (2) the title names the node ------------------------------------------------------------------
        ok &= ck( "tipA".equals( title( tp.buildNodeMenu( a ) ) ), "the menu is titled with the node's name" );
        ok &= ck( "Homo sapiens".equals( title( tp.buildNodeMenu( c ) ) ), "a nameless tip: its taxonomy" );
        ok &= ck( "Internal node".equals( title( tp.buildNodeMenu( root ) ) ), "the nameless root: 'Internal node'" );

        // ---- (3) greyed, never hidden: a tip, the root, a sequence link, a taxonomy link ------------------
        final JPopupMenu on_tip = tp.buildNodeMenu( a );
        for ( final String on : new String[] { "Display Node Data", "Edit Node Data", "Edit Node Style",
                "Switch to Subtree", "Colorize Subtree", "Root/Reroot", "Select/Deselect Node", "Cut Subtree",
                "Copy Subtree", "Add New Node" } ) {
            ok &= ck( enabled( on_tip, on ), "on a tip, '" + on + "' applies" );
        }
        for ( final String off : new String[] { "Collapse/Uncollapse", "Uncollapse Subtree", "Swap Descendants",
                "Ladderize Subtree", "Select/Deselect Tips", "Open Sequence DB", "Open Taxonomy DB" } ) {
            ok &= ck( item( on_tip, off ) != null, "on a tip, '" + off + "' is offered" );
            ok &= ck( !enabled( on_tip, off ), "on a tip, '" + off + "' is greyed" );
        }
        // Paste is greyed only for want of a clipboard: with a clade held it applies on a tip as on any node
        mf.getMainPanel().setCutOrCopiedTree( null );
        ok &= ck( !enabled( tp.buildNodeMenu( a ), "Paste Subtree" ), "with nothing on the clipboard Paste is greyed" );
        tp.runNodeMenuEntry( NodeMenuEntry.COPY_SUBTREE, cd );
        ok &= ck( enabled( tp.buildNodeMenu( a ), "Paste Subtree" ), "with a clade on the clipboard Paste applies, on a tip too" );
        mf.getMainPanel().setCutOrCopiedTree( null );
        mf.getMainPanel().setCopiedAndPastedNodes( null );
        final JMenuItem del_tip = item( on_tip, "Delete Node…" );
        Color danger = UIManager.getColor( "Actions.Red" );
        if ( danger == null ) {
            danger = new Color( 0xC0, 0x39, 0x2B );
        }
        ok &= ck( ( del_tip != null ) && del_tip.isEnabled() && danger.equals( del_tip.getForeground() ),
                  "on a tip the delete entry reads 'Delete Node…', applies, and wears the warning colour" );
        ok &= ck( item( on_tip, "Delete Node/Subtree…" ) == null, "a tip has no subtree to delete" );
        final JPopupMenu on_root = tp.buildNodeMenu( root );
        ok &= ck( !enabled( on_root, "Root/Reroot" )
                  && "this node is the root already".equals( item( on_root, "Root/Reroot" ).getToolTipText() ),
                  "on the root, Root/Reroot is greyed and says why" );
        ok &= ck( !enabled( on_root, "Switch to Subtree" ) && !enabled( on_root, "Cut Subtree" )
                  && !enabled( on_root, "Delete Node/Subtree…" ) && enabled( on_root, "Select/Deselect Tips" )
                  && enabled( on_root, "Swap Descendants" ) && enabled( on_root, "Copy Subtree" ),
                  "on the root: no subtree view, cut or delete; tips, swap and copy apply" );
        ok &= ck( item( on_root, "Root/Reroot" ).getToolTipText() != null && item( on_root, "Copy Subtree" ).getToolTipText() == null,
                  "only a greyed entry carries a reason" );
        ok &= ck( !enabled( on_root, "Collapse/Uncollapse" )
                  && "the root cannot be collapsed".equals( item( on_root, "Collapse/Uncollapse" ).getToolTipText() ),
                  "on the root, Collapse/Uncollapse is greyed and says why (collapse() skips the root)" );
        final JPopupMenu on_b = tp.buildNodeMenu( b );
        ok &= ck( ( item( on_b, "Open Sequence DB [P38398]" ) != null ) && enabled( on_b, "Open Sequence DB [P38398]" )
                  && ( item( on_b, "Open Sequence DB" ) == null ),
                  "a sequence with an accession: one link named by it (once, the name being the accession), in place of the bare entry: " + words( on_b ) );
        ok &= ck( enabled( tp.buildNodeMenu( c ), "Open Taxonomy DB" ), "a taxonomy: the taxonomy link applies" );

        // ---- (4) a subtree view: its root climbs back out, nothing re-roots ------------------------------
        tp.subTree( ab );
        final JPopupMenu on_sub_root = tp.buildNodeMenu( tp.getPhylogeny().getRoot() );
        ok &= ck( ( item( on_sub_root, "Return to Supertree" ) != null ) && enabled( on_sub_root, "Return to Supertree" )
                  && ( item( on_sub_root, "Switch to Subtree" ) == null ),
                  "on a subtree view's root the subtree entry reads 'Return to Supertree': " + words( on_sub_root ) );
        ok &= ck( !enabled( on_sub_root, "Root/Reroot" )
                  && "not in a subtree view: return to the whole tree first".equals( item( on_sub_root, "Root/Reroot" ).getToolTipText() ),
                  "in a subtree view Root/Reroot is greyed and says why" );
        final JPopupMenu on_sub_tip = tp.buildNodeMenu( a );
        ok &= ck( !enabled( on_sub_tip, "Switch to Subtree" ), "a tip under the displayed root has no smaller clade to show" );
        tp.superTree();
        tp.getControlPanel().showWhole();
        render( tp );
        ok &= ck( tp.getPhylogeny().getRoot() == root, "back on the whole tree" );

        // ---- (4b) the layouts: unrooted greys what it cannot offer and says so; circular can collapse ------
        tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
        final JPopupMenu unrooted = tp.buildNodeMenu( cd );
        for ( final String off : new String[] { "Switch to Subtree", "Collapse/Uncollapse", "Uncollapse Subtree",
                "Colorize Subtree", "Root/Reroot", "Cut Subtree", "Copy Subtree", "Paste Subtree", "Add New Node",
                "Delete Node/Subtree\u2026" } ) {
            ok &= ck( !enabled( unrooted, off )
                      && "not in the unrooted layout".equals( item( unrooted, off ).getToolTipText() ),
                      "in the unrooted layout '" + off + "' is greyed and says so" );
        }
        for ( final String on : new String[] { "Display Node Data", "Edit Node Data", "Edit Node Style",
                "Swap Descendants", "Ladderize Subtree", "Select/Deselect Node", "Select/Deselect Tips" } ) {
            ok &= ck( enabled( unrooted, on ), "in the unrooted layout '" + on + "' still applies" );
        }
        tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        final JPopupMenu circular = tp.buildNodeMenu( cd );
        ok &= ck( enabled( circular, "Collapse/Uncollapse" ) && enabled( circular, "Root/Reroot" )
                  && enabled( circular, "Cut Subtree" ) && enabled( circular, "Switch to Subtree" ),
                  "the circular layout offers collapse, re-root, cut and the subtree view" );
        tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        render( tp );

        // ---- (5) every click opens the menu, a shift-click selects, and nothing is sticky ----------------
        final Point at_d = at( d );
        ok &= ck( tp.findNode( at_d.x, at_d.y ) == d, "precondition: D is under its own coordinates" );
        tp.mouseClicked( click( tp, at_d, 0, MouseEvent.BUTTON1 ) );
        final JPopupMenu left = tp.nodeMenuForTest();
        ok &= ck( ( left != null ) && left.isVisible() && "D".equals( title( left ) ), "a plain left click opens D's menu" );
        hide( left );
        tp.mouseClicked( click( tp, at_d, 0, MouseEvent.BUTTON3 ) );
        final JPopupMenu right = tp.nodeMenuForTest();
        ok &= ck( ( right != null ) && ( right != left ) && right.isVisible(), "a right click opens the menu too" );
        hide( right );
        tp.mouseClicked( click( tp, at_d, InputEvent.CTRL_DOWN_MASK, MouseEvent.BUTTON1 ) );
        final JPopupMenu ctrl = tp.nodeMenuForTest();
        ok &= ck( ( ctrl != null ) && ( ctrl != right ) && ctrl.isVisible(), "a ctrl-click opens the menu too" );
        hide( ctrl );
        tp.setFoundNodes0( null );
        tp.mouseClicked( click( tp, at_d, InputEvent.SHIFT_DOWN_MASK, MouseEvent.BUTTON1 ) );
        ok &= ck( ( tp.nodeMenuForTest() == ctrl ) && !ctrl.isVisible() && selected( tp, d ),
                  "a shift-click selects the node and opens no menu" );
        tp.mouseClicked( click( tp, at_d, InputEvent.SHIFT_DOWN_MASK, MouseEvent.BUTTON1 ) );
        ok &= ck( !selected( tp, d ), "a second shift-click deselects it" );
        // nothing sticky: an entry acts once, and the next plain click opens the menu again
        tp.runNodeMenuEntry( NodeMenuEntry.COLLAPSE_UNCOLLAPSE, ab );
        ok &= ck( ab.isCollapse(), "the collapse entry collapses the clade it was chosen on" );
        render( tp );
        final Point at_cd = at( cd );
        ok &= ck( tp.findNode( at_cd.x, at_cd.y ) == cd, "precondition: cd is under its own coordinates" );
        tp.mouseClicked( click( tp, at_cd, 0, MouseEvent.BUTTON1 ) );
        final JPopupMenu again = tp.nodeMenuForTest();
        ok &= ck( ( again != null ) && again.isVisible() && "cd".equals( title( again ) ) && !cd.isCollapse(),
                  "the next plain click opens the menu again -- the previous choice was not kept as a mode" );
        hide( again );
        tp.runNodeMenuEntry( NodeMenuEntry.COLLAPSE_UNCOLLAPSE, ab );
        ok &= ck( !ab.isCollapse(), "chosen again, the entry uncollapses" );
        render( tp );

        // ---- (6) Shift turns the glow into the selection preview; a branch is a target only then ---------
        tp.setFoundNodes0( null );
        final Point d_pt = at( d );
        tp.mouseMoved( move( tp, d_pt, 0 ) );
        ok &= ck( !tp.hoverIsSelectPreviewForTest() && ( tp.hoverNodeForTest() == d )
                  && ( tp.getCursor().getType() == Cursor.HAND_CURSOR ),
                  "without Shift: the neutral focus glow on D and the hand cursor" );
        tp.mouseMoved( move( tp, d_pt, InputEvent.SHIFT_DOWN_MASK ) );
        ok &= ck( tp.hoverIsSelectPreviewForTest() && ( tp.hoverNodeForTest() == d ) && !tp.hoverWouldDeselectForTest(),
                  "with Shift: the selection preview, D unselected so a click adds" );
        final Point cd_mid = branchMid( cd );
        ok &= ck( tp.findNode( cd_mid.x, cd_mid.y ) == null, "precondition: cd's branch midpoint is clear of node boxes" );
        tp.mouseMoved( move( tp, cd_mid, InputEvent.SHIFT_DOWN_MASK ) );
        ok &= ck( ( tp.hoverNodeForTest() == cd ) && tp.hoverIsSubtreeForTest()
                  && ( tp.getCursor().getType() == Cursor.HAND_CURSOR ),
                  "with Shift a branch is a target: its clade is previewed" );
        tp.mouseMoved( move( tp, cd_mid, 0 ) );
        ok &= ck( ( tp.hoverNodeForTest() == null ) && !tp.hoverIsSelectPreviewForTest()
                  && ( tp.getCursor().getType() != Cursor.HAND_CURSOR ),
                  "without Shift a branch is not a target" );
        // the key alone, with the pointer still: pressing Shift starts the preview, releasing it ends it
        for ( final KeyListener kl : tp.getKeyListeners() ) {
            kl.keyPressed( new KeyEvent( tp, KeyEvent.KEY_PRESSED, System.currentTimeMillis(), InputEvent.SHIFT_DOWN_MASK,
                                         KeyEvent.VK_SHIFT, KeyEvent.CHAR_UNDEFINED ) );
        }
        ok &= ck( tp.hoverIsSelectPreviewForTest(), "pressing Shift starts the selection preview without a pointer move" );
        for ( final KeyListener kl : tp.getKeyListeners() ) {
            kl.keyReleased( new KeyEvent( tp, KeyEvent.KEY_RELEASED, System.currentTimeMillis(), 0, KeyEvent.VK_SHIFT,
                                          KeyEvent.CHAR_UNDEFINED ) );
        }
        ok &= ck( !tp.hoverIsSelectPreviewForTest(), "releasing Shift ends it" );
        tp.clearHoverPreview();

        // ---- (7) the clipboard highlight needs no mode: it shows whenever something is held, and fades at the
        //          next edit OF ITS OWN TAB (the clipboard is shared by the tabs; the highlight is the tab's) ---
        mf.getMainPanel().setCutOrCopiedTree( null );
        mf.getMainPanel().setCopiedAndPastedNodes( null );
        final int before = clipboardInk( shot( tp ) );
        tp.runNodeMenuEntry( NodeMenuEntry.COPY_SUBTREE, cd );
        final int after = clipboardInk( shot( tp ) );
        ok &= ck( ( before == 0 ) && ( after > 20 ), "a copied clade is highlighted at once (" + before + " -> " + after + " px)" );
        // an edit in ANOTHER tab is not an edit here
        tabs.setSelectedIndex( 1 );
        final TreePanel other = mf.getMainPanel().getCurrentTreePanel();
        ok &= ck( ( other != tp ) && ( other.getPhylogeny().getNumberOfExternalNodes() == 2 ), "precondition: the second tab" );
        other.runNodeMenuEntry( NodeMenuEntry.SWAP, other.getPhylogeny().getRoot() );
        tabs.setSelectedIndex( 0 );
        ok &= ck( clipboardInk( shot( tp ) ) > 20, "an edit in another tab leaves this tab's highlight" );
        // ... the next edit HERE fades it (every tool checkpoints before it mutates), the clipboard itself kept
        tp.runNodeMenuEntry( NodeMenuEntry.SWAP, root );
        ok &= ck( clipboardInk( shot( tp ) ) == 0, "the next edit fades the highlight" );
        ok &= ck( tp.isCanPaste(), "...but the clipboard still holds the clade: Paste stays possible" );
        // ... and at an undo or redo (the restored tree is another copy)
        tp.runNodeMenuEntry( NodeMenuEntry.COPY_SUBTREE, cd );
        ok &= ck( clipboardInk( shot( tp ) ) > 20, "precondition: copied again, highlighted again" );
        ok &= ck( tp.undo(), "precondition: the swap can be undone" );
        ok &= ck( clipboardInk( shot( tp ) ) == 0, "an undo fades the highlight too" );
        tp.redo();
        // ... and it needs the clipboard: a fresh copy, then the clipboard emptied under it
        tp.runNodeMenuEntry( NodeMenuEntry.COPY_SUBTREE, cd );
        ok &= ck( clipboardInk( shot( tp ) ) > 20, "precondition: copied once more" );
        mf.getMainPanel().setCutOrCopiedTree( null );
        ok &= ck( clipboardInk( shot( tp ) ) == 0, "...and only while the clipboard holds it" );
        mf.getMainPanel().setCopiedAndPastedNodes( null );
        return ok;
    }

    /** subTree() nests at most 100 views deep: at the cap the subtree entry is greyed and says so; one level up
     *  it applies again. */
    private static boolean subtreeCapOk() throws Exception {
        final PhylogenyNode root = new PhylogenyNode();
        PhylogenyNode cur = root;
        for ( int i = 1; i <= 103; ++i ) {
            child( cur, "t" + i );
            cur = child( cur, "n" + i );
        }
        child( cur, "tA" );
        child( cur, "tB" );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "subtree cap" ) );
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    tp.getOptions().setShowOverview( false );
                    tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    render( tp );
                    for ( int i = 1; i <= 100; ++i ) {
                        tp.subTree( tp.getPhylogeny().getNode( "n" + i ) );
                    }
                    ok[ 0 ] &= ck( tp.isCurrentTreeIsSubtree() && "n100".equals( tp.getPhylogeny().getRoot().getName() ),
                                   "precondition: 100 subtree views deep, on " + tp.getPhylogeny().getRoot().getName() );
                    final JPopupMenu at_cap = tp.buildNodeMenu( tp.getPhylogeny().getNode( "n101" ) );
                    ok[ 0 ] &= ck( !enabled( at_cap, "Switch to Subtree" )
                                   && "no deeper subtree view: return to a supertree first".equals( item( at_cap, "Switch to Subtree" ).getToolTipText() ),
                                   "at the cap the subtree entry is greyed and says so" );
                    tp.superTree();
                    ok[ 0 ] &= ck( enabled( tp.buildNodeMenu( tp.getPhylogeny().getNode( "n100" ) ), "Switch to Subtree" ),
                                   "one level up it applies again" );
                }
                catch ( final Throwable e ) {
                    e.printStackTrace();
                    ok[ 0 ] = false;
                }
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ].dispose() );
        }
        return ok[ 0 ];
    }

    // ---- helpers --------------------------------------------------------------------------------------------

    private static PhylogenyNode child( final PhylogenyNode parent, final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        if ( name != null ) {
            n.setName( name );
        }
        n.setDistanceToParent( 1 );
        parent.addAsChild( n );
        return n;
    }

    /** The menu's words in order: "|" for a separator, else the item's text; the title is not among them. */
    private static List<String> words( final JPopupMenu menu ) {
        final List<String> out = new ArrayList<>();
        for ( final Component c : menu.getComponents() ) {
            if ( c instanceof JPopupMenu.Separator ) {
                out.add( "|" );
            }
            else if ( c instanceof JMenuItem ) {
                out.add( ( (JMenuItem) c ).getText() );
            }
        }
        return out;
    }

    private static String title( final JPopupMenu menu ) {
        final Component first = ( menu.getComponentCount() > 0 ) ? menu.getComponent( 0 ) : null;
        return ( first instanceof JLabel ) ? ( (JLabel) first ).getText() : null;
    }

    /** The menu's item with exactly this text, or null (shared with the other tests that read the node menu). */
    static JMenuItem item( final JPopupMenu menu, final String text ) {
        for ( final Component c : menu.getComponents() ) {
            if ( ( c instanceof JMenuItem ) && text.equals( ( (JMenuItem) c ).getText() ) ) {
                return (JMenuItem) c;
            }
        }
        return null;
    }

    private static boolean enabled( final JPopupMenu menu, final String text ) {
        final JMenuItem i = item( menu, text );
        return ( i != null ) && i.isEnabled();
    }

    private static void hide( final JPopupMenu menu ) {
        if ( menu != null ) {
            menu.setVisible( false );
        }
    }

    private static boolean selected( final TreePanel tp, final PhylogenyNode n ) {
        return ( tp.getFoundNodes0() != null ) && tp.getFoundNodes0().contains( n.getId() );
    }

    private static Point at( final PhylogenyNode n ) {
        return new Point( Math.round( n.getXcoord() ), Math.round( n.getYcoord() ) );
    }

    /** The midpoint of {@code n}'s incoming horizontal branch (rectangular, root at the left). */
    private static Point branchMid( final PhylogenyNode n ) {
        final PhylogenyNode p = n.getParent();
        return new Point( Math.round( ( p.getXcoord() + n.getXcoord() ) / 2f ), Math.round( n.getYcoord() ) );
    }

    private static MouseEvent click( final TreePanel tp, final Point p, final int modifiers, final int button ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(), modifiers, p.x, p.y, 1, false,
                               button );
    }

    private static MouseEvent move( final TreePanel tp, final Point p, final int modifiers ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_MOVED, System.currentTimeMillis(), modifiers, p.x, p.y, 0, false );
    }

    private static void render( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        shot( tp ); // the paint assigns the node coordinates
    }

    private static BufferedImage shot( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.printAll( g );
        g.dispose();
        return img;
    }

    /** Pixels within a small distance of the clipboard highlight colour. */
    private static int clipboardInk( final BufferedImage img ) {
        final Color h = TreePanel.CLIPBOARD_HIGHLIGHT;
        int n = 0;
        for ( int y = 0; y < img.getHeight(); ++y ) {
            for ( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int dr = Math.abs( ( ( rgb >> 16 ) & 0xff ) - h.getRed() );
                final int dg = Math.abs( ( ( rgb >> 8 ) & 0xff ) - h.getGreen() );
                final int db = Math.abs( ( rgb & 0xff ) - h.getBlue() );
                if ( ( dr + dg + db ) < 30 ) {
                    n++;
                }
            }
        }
        return n;
    }

    private static boolean ck( final boolean condition, final String what ) {
        if ( !condition ) {
            System.out.println( "  [NodeMenuTest] " + what );
        }
        return condition;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [NodeMenuTest] " + message );
        return false;
    }

    private NodeMenuTest() {
        // not instantiable
    }
}
