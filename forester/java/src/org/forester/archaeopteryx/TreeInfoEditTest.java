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
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/**
 * Integration test with teeth for {@link TreeInfoDialog#apply}: renaming a tree and setting its description
 * must (1) trim and store both fields, (2) rename the current tab to match the name so the two cannot drift,
 * (3) mark the panel edited, (4) be undoable -- undo restoring BOTH the old name and old description (they
 * ride along on {@link Phylogeny#copy()}), and (5) be a strict no-op when the (trimmed) values are unchanged,
 * pushing no undo checkpoint. Headless-guarded (needs FlatLaf via {@code createInstance}); run standalone or
 * in the non-headless suite.
 */
public final class TreeInfoEditTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreeInfoEdit: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        return mainScenario() && emptyNameUndoScenario();
    }

    /** A named tree: rename/undo/redo, whitespace collapse, description preservation, blank-name and sub-tree guards. */
    private static boolean mainScenario() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { threeLevel( "orig" ) }, conf,
                                                                        "orig" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    final int sel = mp.getTabbedPane().getSelectedIndex();

                    // (5) unchanged values (description just whitespace -> trims to the existing "") is a no-op:
                    // no change reported, and crucially no undo checkpoint pushed.
                    final boolean could_undo_before = tp.canUndo();
                    if ( TreeInfoDialog.apply( mf[ 0 ], tp, "orig", "   " ) ) {
                        fail( ok, "unchanged apply should report no change" );
                    }
                    if ( tp.canUndo() != could_undo_before ) {
                        fail( ok, "a no-op apply must not push an undo checkpoint" );
                    }

                    // (1)+(2)+(3) real edit: values are trimmed and stored, tab renamed, panel edited.
                    if ( !TreeInfoDialog.apply( mf[ 0 ], tp, "  renamed  ", "  a helpful description  " ) ) {
                        fail( ok, "a real change should report changed" );
                    }
                    final Phylogeny phy = tp.getPhylogeny();
                    if ( !"renamed".equals( phy.getName() ) ) {
                        fail( ok, "name should be trimmed+stored, got \"" + phy.getName() + "\"" );
                    }
                    if ( !"a helpful description".equals( phy.getDescription() ) ) {
                        fail( ok, "description should be trimmed+stored, got \"" + phy.getDescription() + "\"" );
                    }
                    if ( !"renamed".equals( mp.getTabbedPane().getTitleAt( sel ) ) ) {
                        fail( ok, "tab title should track the name, got \"" + mp.getTabbedPane().getTitleAt( sel )
                                + "\"" );
                    }
                    if ( !tp.isEdited() ) {
                        fail( ok, "panel should be marked edited after a change" );
                    }
                    if ( !tp.canUndo() ) {
                        fail( ok, "the change should be undoable" );
                    }

                    // (4) undo (via MainFrame, as the ⌘Z menu does) restores BOTH old name and old (empty)
                    // description AND re-syncs the tab title -- the invariant that name and tab cannot drift.
                    mf[ 0 ].undo();
                    final Phylogeny restored = tp.getPhylogeny();
                    if ( !"orig".equals( restored.getName() ) ) {
                        fail( ok, "undo should restore the old name, got \"" + restored.getName() + "\"" );
                    }
                    if ( !ForesterUtil.isEmpty( restored.getDescription() ) ) {
                        fail( ok, "undo should restore the old (empty) description, got \""
                                + restored.getDescription() + "\"" );
                    }
                    if ( !"orig".equals( mp.getTabbedPane().getTitleAt( sel ) ) ) {
                        fail( ok, "undo should re-sync the tab title to the restored name, got \""
                                + mp.getTabbedPane().getTitleAt( sel ) + "\"" );
                    }

                    // (4b) redo re-applies the rename AND re-syncs the tab title (redo has the same sync as undo).
                    mf[ 0 ].redo();
                    if ( !"renamed".equals( tp.getPhylogeny().getName() ) ) {
                        fail( ok, "redo should re-apply the name, got \"" + tp.getPhylogeny().getName() + "\"" );
                    }
                    if ( !"renamed".equals( mp.getTabbedPane().getTitleAt( sel ) ) ) {
                        fail( ok, "redo should re-sync the tab title, got \"" + mp.getTabbedPane().getTitleAt( sel )
                                + "\"" );
                    }
                    mf[ 0 ].undo(); // back to "orig" for the following steps

                    // (5) the name is whitespace-cleaned on apply: internal runs (incl. a tab) collapse to single
                    // spaces, since it is a one-line label -- both on the tree and on the tab title.
                    if ( !TreeInfoDialog.apply( mf[ 0 ], tp, "  multi   word\tname  ", "" ) ) {
                        fail( ok, "a whitespace-collapsing rename should report changed" );
                    }
                    if ( !"multi word name".equals( tp.getPhylogeny().getName() ) ) {
                        fail( ok, "internal name whitespace should collapse, got \"" + tp.getPhylogeny().getName()
                                + "\"" );
                    }
                    if ( !"multi word name".equals( mp.getTabbedPane().getTitleAt( sel ) ) ) {
                        fail( ok, "tab title should show the cleaned name, got \""
                                + mp.getTabbedPane().getTitleAt( sel ) + "\"" );
                    }
                    mf[ 0 ].undo(); // back to "orig" for the following steps

                    // (5b) UNLIKE the name, the description is only end-trimmed -- internal runs and newlines are
                    // preserved (it is multi-line free text), so collapsing would be a regression.
                    final String multiline = "para one\n\n  indented   two";
                    if ( !TreeInfoDialog.apply( mf[ 0 ], tp, "orig", "  " + multiline + "  " ) ) {
                        fail( ok, "setting a multi-line description should report changed" );
                    }
                    if ( !multiline.equals( tp.getPhylogeny().getDescription() ) ) {
                        fail( ok, "description whitespace/newlines must be preserved (only ends trimmed), got \""
                                + tp.getPhylogeny().getDescription() + "\"" );
                    }

                    // (6) a blank/whitespace name is not a rename: the existing name and tab title are kept (not
                    // blanked); only the description changes.
                    if ( !TreeInfoDialog.apply( mf[ 0 ], tp, "   ", "desc only" ) ) {
                        fail( ok, "a description-only change (blank name) should report changed" );
                    }
                    if ( !"orig".equals( tp.getPhylogeny().getName() ) ) {
                        fail( ok, "a blank name must keep the existing name, got \"" + tp.getPhylogeny().getName()
                                + "\"" );
                    }
                    if ( !"orig".equals( mp.getTabbedPane().getTitleAt( sel ) ) ) {
                        fail( ok, "a blank name must not blank the tab title, got \""
                                + mp.getTabbedPane().getTitleAt( sel ) + "\"" );
                    }
                    if ( !"desc only".equals( tp.getPhylogeny().getDescription() ) ) {
                        fail( ok, "description should still update on a blank-name edit, got \""
                                + tp.getPhylogeny().getDescription() + "\"" );
                    }

                    // (6b) a blank name with an UNCHANGED description is a pure no-op: no change, no checkpoint.
                    final boolean could_undo_noop = tp.canUndo();
                    if ( TreeInfoDialog.apply( mf[ 0 ], tp, "  ", "desc only" ) ) {
                        fail( ok, "blank name + unchanged description should report no change" );
                    }
                    if ( tp.canUndo() != could_undo_noop ) {
                        fail( ok, "a blank-name no-op must not push an undo checkpoint" );
                    }

                    // (7) editing is refused while a SUB-TREE is displayed (the edit would be discarded on
                    // returning to the whole tree), so apply() is a no-op and the sub-tree name is untouched.
                    final PhylogenyNode inner = innerNode( tp.getPhylogeny() );
                    if ( inner == null ) {
                        fail( ok, "test tree should have an internal non-root node to form a sub-tree" );
                    }
                    else {
                        tp.subTree( inner );
                        if ( !tp.isCurrentTreeIsSubtree() ) {
                            fail( ok, "sub-tree navigation should have entered sub-tree mode" );
                        }
                        if ( TreeInfoDialog.apply( mf[ 0 ], tp, "shouldNotStick", "nope" ) ) {
                            fail( ok, "apply must be a no-op while a sub-tree is displayed" );
                        }
                        if ( "shouldNotStick".equals( tp.getPhylogeny().getName() ) ) {
                            fail( ok, "a refused sub-tree edit must not change the displayed tree's name" );
                        }
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose(); // never leak the frame into the shared suite JVM
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * A previously-UNNAMED tree (tab titled from its file): naming it and then undoing must NOT leave the undone
     * name on the tab, because the save-time backfill would otherwise persist that name -- undo would be defeated
     * on disk. The tab must revert to the original file-derived placeholder.
     */
    private static boolean emptyNameUndoScenario() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { twoTipUnnamed() }, conf, "myfile.xml" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    final int sel = mp.getTabbedPane().getSelectedIndex();
                    if ( !ForesterUtil.isEmpty( tp.getPhylogeny().getName() ) ) {
                        fail( ok, "fixture tree should be unnamed, got \"" + tp.getPhylogeny().getName() + "\"" );
                    }
                    final String placeholder = mp.getTabbedPane().getTitleAt( sel ); // file-derived, e.g. "myfile.xml"
                    if ( ForesterUtil.isEmpty( placeholder ) ) {
                        fail( ok, "an unnamed tree's tab should still have a non-empty placeholder title" );
                    }
                    if ( !TreeInfoDialog.apply( mf[ 0 ], tp, "Primates", "" ) ) {
                        fail( ok, "naming a previously-unnamed tree should report changed" );
                    }
                    if ( !"Primates".equals( mp.getTabbedPane().getTitleAt( sel ) ) ) {
                        fail( ok, "tab should show the new name after naming" );
                    }
                    mf[ 0 ].undo();
                    if ( !ForesterUtil.isEmpty( tp.getPhylogeny().getName() ) ) {
                        fail( ok, "undo should restore the empty name, got \"" + tp.getPhylogeny().getName() + "\"" );
                    }
                    if ( "Primates".equals( mp.getTabbedPane().getTitleAt( sel ) ) ) {
                        fail( ok, "undo must not leave the undone name on the tab (a save would resurrect it)" );
                    }
                    if ( !placeholder.equals( mp.getTabbedPane().getTitleAt( sel ) ) ) {
                        fail( ok, "undo should restore the original placeholder tab title \"" + placeholder
                                + "\", got \"" + mp.getTabbedPane().getTitleAt( sel ) + "\"" );
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [TreeInfoEditTest] " + msg );
        ok[ 0 ] = false;
    }

    // root -> ( X, inner( A, B ) ); the tree carries a name so the rename/undo assertions have a known start,
    // and the internal non-root node `inner` lets the test navigate into a sub-tree
    private static Phylogeny threeLevel( final String name ) {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode x = new PhylogenyNode();
        x.setName( "X" );
        root.addAsChild( x );
        final PhylogenyNode inner = new PhylogenyNode();
        root.addAsChild( inner );
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "A" );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "B" );
        inner.addAsChild( a );
        inner.addAsChild( b );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setName( name );
        phy.externalNodesHaveChanged();
        return phy;
    }

    // root -> ( A, B ) with NO tree name, so the tab title comes from the loaded file placeholder instead
    private static Phylogeny twoTipUnnamed() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "A" );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "B" );
        root.addAsChild( a );
        root.addAsChild( b );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy; // Phylogeny.init() leaves the name as "" -> unnamed
    }

    /** The first internal, non-root node in preorder (the sub-tree seed), or null if none. */
    private static PhylogenyNode innerNode( final Phylogeny phy ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.isRoot() ) {
                return n;
            }
        }
        return null;
    }

    private TreeInfoEditTest() {
    }
}
