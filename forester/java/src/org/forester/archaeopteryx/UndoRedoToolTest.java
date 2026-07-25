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
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;

/**
 * Integration test for Undo/Redo on a real {@link MainFrame}/{@link TreePanel}: checks the Edit menu exists
 * and is enabled/labeled correctly, and that a checkpointed mutation round-trips through undo (tree + edited
 * flag restored) and redo (re-applied), with a new checkpoint clearing the redo history. Guarded to a no-op
 * on a headless box (needs FlatLaf via {@code createInstance}).
 */
public final class UndoRedoToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "UndoRedoTool: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return basicUndoRedo() && undoWithDomainArchitectures();
    }

    private static boolean basicUndoRedo() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree( 6 ) }, conf,
                                                                        "undo" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                final Phylogeny phy = tp.getPhylogeny();

                // Edit menu present; Undo/Redo start disabled
                if ( menuByTitle( frame.getJMenuBar(), "Edit" ) == null ) {
                    fail( ok, "Edit menu not found" );
                }
                frame.updateEditMenu();
                if ( frame._undo_item.isEnabled() || frame._redo_item.isEnabled() ) {
                    fail( ok, "Undo/Redo should start disabled" );
                }
                if ( !"Undo".equals( frame._undo_item.getText() ) ) {
                    fail( ok, "disabled Undo label should be plain 'Undo'" );
                }

                // checkpoint + mutate (delete a tip), like a real Delete Nodes op
                final int before = phy.getNumberOfExternalNodes(); // 6
                tp.pushUndoCheckpoint( "Delete Nodes" );
                phy.deleteSubtree( phy.getExternalNodes().get( 0 ), true );
                phy.externalNodesHaveChanged();
                phy.clearHashIdToNodeMap();
                phy.recalculateNumberOfExternalDescendants( true );
                tp.setEdited( true );
                if ( !tp.canUndo() || tp.canRedo() ) {
                    fail( ok, "after a checkpoint: canUndo true, canRedo false" );
                }
                frame.updateEditMenu();
                if ( !frame._undo_item.isEnabled() || !"Undo Delete Nodes".equals( frame._undo_item.getText() ) ) {
                    fail( ok, "Undo item should be enabled and labeled 'Undo Delete Nodes'" );
                }

                // undo restores the tree AND the pre-edit (clean) flag
                if ( !tp.undo() ) {
                    fail( ok, "undo should report success" );
                }
                if ( tp.getPhylogeny().getNumberOfExternalNodes() != before ) {
                    fail( ok, "undo should restore the original tip count, got "
                            + tp.getPhylogeny().getNumberOfExternalNodes() );
                }
                if ( tp.isEdited() ) {
                    fail( ok, "undo back to the loaded state should restore a clean (not-edited) marker" );
                }
                // tp.undo() notifies the frame, so the Edit menu is already refreshed
                if ( !tp.canRedo() || !frame._redo_item.isEnabled()
                        || !"Redo Delete Nodes".equals( frame._redo_item.getText() ) ) {
                    fail( ok, "after undo: Redo should be enabled and labeled 'Redo Delete Nodes'" );
                }

                // redo re-applies the deletion and the edited flag
                if ( !tp.redo() ) {
                    fail( ok, "redo should report success" );
                }
                if ( tp.getPhylogeny().getNumberOfExternalNodes() != ( before - 1 ) ) {
                    fail( ok, "redo should re-delete the tip, got "
                            + tp.getPhylogeny().getNumberOfExternalNodes() );
                }
                if ( !tp.isEdited() ) {
                    fail( ok, "redo should restore the edited flag" );
                }

                // safety net: an edit made WITHOUT a checkpoint clears any pending redo, so a stale Redo can
                // never install an unrelated tree
                if ( tp.undo() && tp.canRedo() ) {
                    tp.setEdited( true ); // simulate an uncheckpointed mutation
                    if ( tp.canRedo() ) {
                        fail( ok, "an uncheckpointed edit must clear the redo history (safety net)" );
                    }
                }

                // the real Edit-menu handlers (frame.undo()/redo()), not just the panel methods, round-trip
                final Phylogeny cur = tp.getPhylogeny();
                final int n0 = cur.getNumberOfExternalNodes();
                tp.pushUndoCheckpoint( "Delete Nodes" );
                cur.deleteSubtree( cur.getExternalNodes().get( 0 ), true );
                cur.externalNodesHaveChanged();
                cur.clearHashIdToNodeMap();
                cur.recalculateNumberOfExternalDescendants( true );
                tp.setEdited( true );
                frame.undo();
                if ( tp.getPhylogeny().getNumberOfExternalNodes() != n0 ) {
                    fail( ok, "frame.undo() should restore the pre-edit tip count" );
                }
                frame.redo();
                if ( tp.getPhylogeny().getNumberOfExternalNodes() != ( n0 - 1 ) ) {
                    fail( ok, "frame.redo() should re-apply the edit" );
                }

                ( (JFrame) frame ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [UndoRedoToolTest] " + msg );
        ok[ 0 ] = false;
    }

    private static JMenu menuByTitle( final JMenuBar bar, final String title ) {
        for( int i = 0; i < bar.getMenuCount(); ++i ) {
            final JMenu m = bar.getMenu( i );
            if ( ( m != null ) && title.equals( m.getText() ) ) {
                return m;
            }
        }
        return null;
    }

    /**
     * Undoing an edit on a tree that has domain architectures must not crash: Phylogeny.copy() yields plain
     * DomainArchitectures, and the layout's calculateLongestExtNodeInfo casts them to
     * RenderableDomainArchitecture -- so restoreSnapshot must re-wrap them first.
     */
    private static boolean undoWithDomainArchitectures() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { treeWithDomains() }, conf,
                                                                        "domains" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                frame.getMainPanel().getControlPanel().setCheckbox( Configuration.show_domain_architectures, true );
                if ( !frame.getMainPanel().getControlPanel().isShowDomainArchitectures() ) {
                    System.out.println( "  [UndoRedoToolTest] note: domain-architecture display not available; "
                            + "skipping the domain-undo assertion" );
                    ( (JFrame) frame ).dispose();
                    return;
                }
                final Phylogeny phy = tp.getPhylogeny();
                final int before = phy.getNumberOfExternalNodes();
                tp.pushUndoCheckpoint( "Delete Node" );
                phy.deleteSubtree( phy.getExternalNodes().get( 0 ), true );
                phy.externalNodesHaveChanged();
                phy.clearHashIdToNodeMap();
                phy.recalculateNumberOfExternalDescendants( true );
                tp.setEdited( true );
                try {
                    tp.undo(); // restoreSnapshot -> calculateLongestExtNodeInfo casts the domain architectures
                }
                catch ( final Throwable t ) {
                    System.out.println( "  [UndoRedoToolTest] undo with domain architectures threw: " + t );
                    ok[ 0 ] = false;
                }
                if ( tp.getPhylogeny().getNumberOfExternalNodes() != before ) {
                    System.out.println( "  [UndoRedoToolTest] undo did not restore the domain-bearing tree" );
                    ok[ 0 ] = false;
                }
                ( (JFrame) frame ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static Phylogeny treeWithDomains() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 3; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            final Sequence seq = new Sequence();
            final List<PhylogenyData> domains = new ArrayList<PhylogenyData>();
            domains.add( new ProteinDomain( "PF000" + i, 1, 40 ) );
            seq.setDomainArchitecture( new DomainArchitecture( domains, 100 ) );
            leaf.getNodeData().addSequence( seq );
            root.addAsChild( leaf );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static Phylogeny tree( final int n ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < n; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            root.addAsChild( leaf );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }
}
