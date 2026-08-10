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
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.tools.AncestralTaxonomyInferrer;
import org.forester.archaeopteryx.tools.SequenceAndTaxonomyDataObtainer;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;
import org.forester.ws.seqdb.SequenceTaxonomyResolver;

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
        return basicUndoRedo() && undoWithDomainArchitectures() && searchBoxTextUndo() && asyncToolCheckpoints();
    }

    /**
     * The async data tools checkpoint at their EDT commit so their (potentially large) mutations are undoable:
     * Infer always checkpoints "Infer Ancestor Taxonomies"; Fetch checkpoints "Fetch Sequence & Taxonomy Data"
     * only when the run actually wrote something. Drives the extracted commit() paths directly (no network).
     */
    private static boolean asyncToolCheckpoints() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree( 5 ) }, conf,
                                                                        "async" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();

                // Infer: commit(assigned) checkpoints "Infer Ancestor Taxonomies", stamps provenance, sets edited
                new AncestralTaxonomyInferrer( frame, tp, tp.getPhylogeny().copy(), false ).commit( 3 );
                if ( !tp.canUndo() || !"Infer Ancestor Taxonomies".equals( tp.undoLabel() ) ) {
                    fail( ok, "infer commit should checkpoint 'Infer Ancestor Taxonomies', got '" + tp.undoLabel() + "'" );
                }
                if ( !tp.isEdited() ) {
                    fail( ok, "infer commit should mark the tree edited" );
                }
                final String infer_desc = tp.getPhylogeny().getDescription();
                if ( ( infer_desc == null ) || !infer_desc.contains( "ancestral-taxonomy inference" )
                        || !infer_desc.contains( "3 internal nodes" ) ) {
                    fail( ok, "infer commit should append a provenance sentence, got: " + infer_desc );
                }

                final SequenceAndTaxonomyDataObtainer obtainer = new SequenceAndTaxonomyDataObtainer(
                        (MainFrameApplication) frame, tp, tp.getPhylogeny().copy() );
                // a no-write fetch result must NOT checkpoint (undo label stays the infer one)
                obtainer.commit( result( 0, 0 ) );
                if ( !"Infer Ancestor Taxonomies".equals( tp.undoLabel() ) ) {
                    fail( ok, "a no-write fetch must not checkpoint, undo label is now '" + tp.undoLabel() + "'" );
                }
                // a result with writes DOES checkpoint "Fetch Sequence & Taxonomy Data"
                obtainer.commit( result( 2, 1 ) );
                if ( !tp.canUndo() || !"Fetch Sequence & Taxonomy Data".equals( tp.undoLabel() ) ) {
                    fail( ok, "fetch commit with writes should checkpoint 'Fetch Sequence & Taxonomy Data', got '"
                            + tp.undoLabel() + "'" );
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

    private static SequenceTaxonomyResolver.Result result( final int sequences_written, final int taxonomies_written ) {
        return new SequenceTaxonomyResolver.Result( new TreeSet<String>(), sequences_written, taxonomies_written,
                                                    false, null );
    }

    /**
     * The Undo shortcut (Cmd/Ctrl-Z) must edit the search text field when it has focus, not revert the tree:
     * the field carries its own WHEN_FOCUSED text-undo binding that takes precedence over the Edit menu's
     * window-scoped accelerator. Verifies the binding is installed and that firing it reverts the typed text
     * while leaving the tree (tip count + edited flag) untouched.
     */
    private static boolean searchBoxTextUndo() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree( 6 ) }, conf,
                                                                        "search-undo" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                final JTextField sf = frame.getMainPanel().getControlPanel().getSearchTextField0();
                if ( sf == null ) {
                    fail( ok, "search field 0 not found" );
                    ( (JFrame) frame ).dispose();
                    return;
                }
                // the text-undo/redo actions and the Cmd/Ctrl-Z WHEN_FOCUSED binding are installed
                if ( ( sf.getActionMap().get( "text-undo" ) == null )
                        || ( sf.getActionMap().get( "text-redo" ) == null ) ) {
                    fail( ok, "search field is missing the text-undo/redo actions" );
                }
                final int sc = Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx();
                final Object binding = sf.getInputMap( JComponent.WHEN_FOCUSED )
                        .get( KeyStroke.getKeyStroke( KeyEvent.VK_Z, sc ) );
                if ( !"text-undo".equals( binding ) ) {
                    fail( ok, "Cmd/Ctrl-Z is not bound to text-undo on the focused search field" );
                }
                // typing then firing text-undo reverts the text, and does NOT touch the tree
                final int tips = tp.getPhylogeny().getNumberOfExternalNodes();
                final boolean edited_before = tp.isEdited();
                sf.setText( "hello" );
                sf.getActionMap().get( "text-undo" )
                        .actionPerformed( new ActionEvent( sf, ActionEvent.ACTION_PERFORMED, "text-undo" ) );
                if ( !sf.getText().isEmpty() ) {
                    fail( ok, "text-undo should revert the typed search text, got '" + sf.getText() + "'" );
                }
                if ( ( tp.getPhylogeny().getNumberOfExternalNodes() != tips ) || ( tp.isEdited() != edited_before ) ) {
                    fail( ok, "text-undo in the search box must not change the tree" );
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
