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
import javax.swing.JTree;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import org.forester.archaeopteryx.tools.AncestralTaxonomyInferrer;
import org.forester.archaeopteryx.tools.SequenceAndTaxonomyDataObtainer;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
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
        return basicUndoRedo() && undoWithDomainArchitectures() && searchBoxTextUndo() && asyncToolCheckpoints()
                && collapseUndo() && nodeEditUndo();
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
                frame.getMainPanel().getControlPanel().setCheckbox( DisplayOption.SHOW_DOMAIN_ARCHITECTURES, true );
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

    /**
     * Collapsing a clade is undoable. Collapse is DISPLAY state, but it lives on the tree
     * ({@code PhylogenyNode._collapse}, which {@code copyNodeData()} carries), so the snapshot history restores
     * it like any other mutation -- and {@code restoreSnapshot} already refreshes the collapse-derived caches.
     * An uncollapse that would change nothing must NOT push an entry, or merely pressing the button on an
     * expanded tree would bury the user's real edits in a 25-deep history and clear the redo stack.
     */
    private static boolean collapseUndo() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { nestedTree() }, conf, "collapse" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                // Pin the layout: collapse()/uncollapseAll() pop a MODAL dialog in UNROOTED, and a modal opened
                // from inside invokeAndWait never returns -- it would hang the whole suite. A standalone run
                // inherits the developer's persisted display type, so this cannot be left to chance.
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                // pure guard first: it is what keeps a no-op uncollapse out of the history
                if ( TreePanel.hasCollapsedNodeIn( tp.getPhylogeny().getRoot() )
                        || TreePanel.hasCollapsedNodeIn( null ) ) {
                    fail( ok, "a freshly loaded tree has nothing collapsed" );
                }
                // an uncollapse with nothing collapsed must leave the history alone
                mf[ 0 ].getMainPanel().getControlPanel().uncollapseAll( tp );
                if ( tp.canUndo() ) {
                    fail( ok, "uncollapse-all on an expanded tree must not push an undo entry, got '"
                            + tp.undoLabel() + "'" );
                }
                // collapse -> undoable, and the label says what happened
                tp.collapse( named( tp, "cladeA" ) );
                if ( !named( tp, "cladeA" ).isCollapse() ) {
                    fail( ok, "precondition: the clade should be collapsed" );
                }
                if ( !tp.canUndo() || !"Collapse Clade".equals( tp.undoLabel() ) ) {
                    fail( ok, "collapsing should checkpoint 'Collapse Clade', got '" + tp.undoLabel() + "'" );
                }
                if ( !TreePanel.hasCollapsedNodeIn( tp.getPhylogeny().getRoot() ) ) {
                    fail( ok, "hasCollapsedNodeIn should see the collapsed clade" );
                }
                // undo restores the EXPANDED state (the snapshot must carry _collapse through Phylogeny.copy())
                tp.undo();
                if ( named( tp, "cladeA" ).isCollapse() ) {
                    fail( ok, "undo must un-collapse the clade" );
                }
                // ...and redo re-collapses it
                tp.redo();
                if ( !named( tp, "cladeA" ).isCollapse() ) {
                    fail( ok, "redo must re-collapse the clade" );
                }
                // uncollapse-all now HAS something to do -> its own entry, and undo brings the collapse back
                mf[ 0 ].getMainPanel().getControlPanel().uncollapseAll( tp );
                if ( !"Uncollapse All".equals( tp.undoLabel() ) ) {
                    fail( ok, "uncollapse-all should checkpoint 'Uncollapse All', got '" + tp.undoLabel() + "'" );
                }
                if ( named( tp, "cladeA" ).isCollapse() ) {
                    fail( ok, "precondition: uncollapse-all should have expanded the clade" );
                }
                tp.undo();
                if ( !named( tp, "cladeA" ).isCollapse() ) {
                    fail( ok, "undoing an uncollapse-all must bring the collapsed clade back" );
                }
            } );
            SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * A node-data edit is undoable as ONE step, and merely OPENING the editor is not an edit.
     * <p>
     * That second half is why this was deferred for so long: {@code NodeEditPanel.writeBack} commits fields on
     * every selection change and again on close, with no change detection, so a checkpoint taken when the editor
     * opens would push a no-op undo -- and wipe the redo stack -- for someone who only opened a node to read it.
     * The checkpoint is instead taken on the first write that FOLLOWS a committed cell edit.
     */
    private static boolean nodeEditUndo() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { nestedTree() }, conf, "nodeedit" ) );
            final boolean[] ok = { true };
            final NodeFrame[] nf = new NodeFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                // Pin the layout: collapse()/uncollapseAll() pop a MODAL dialog in UNROOTED, and a modal opened
                // from inside invokeAndWait never returns -- it would hang the whole suite. A standalone run
                // inherits the developer's persisted display type, so this cannot be left to chance.
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                // Two collapses then one undo leaves BOTH stacks non-empty, so opening the editor can be shown
                // to disturb neither.
                tp.collapse( named( tp, "cladeA" ) );
                tp.collapse( named( tp, "cladeB" ) );
                tp.undo();
                if ( !tp.canUndo() || !tp.canRedo() ) {
                    fail( ok, "precondition: both stacks should be non-empty before opening the editor" );
                }
                // AFTER the undo: it installs the snapshot as the new live tree, so a node taken before it would
                // be detached from the tree the editor is supposed to edit
                final PhylogenyNode n = tp.getPhylogeny().getFirstExternalNode();
                tp.showNodeEditFrame( n ); // the production path: the panel TRACKS the frame it opens
                nf[ 0 ] = openNodeFrame();
                if ( ( nf[ 0 ] == null ) || ( tp.openNodeFrameCountForTest() != 1 ) ) {
                    fail( ok, "the editor should be open and tracked by the panel" );
                    return;
                }
                // opening alone must leave the history completely untouched
                if ( !"Collapse Clade".equals( tp.undoLabel() ) ) {
                    fail( ok, "merely opening the node editor must not push an undo entry, got '"
                            + tp.undoLabel() + "'" );
                }
                if ( !tp.canRedo() ) {
                    fail( ok, "merely opening the node editor must not touch the history at all" );
                }
            } );
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final JTree jt = findJTree( nf[ 0 ] );
                if ( jt == null ) {
                    fail( ok, "the node editor should contain a JTree" );
                    return;
                }
                for( int i = 0; i < jt.getRowCount(); ++i ) {
                    jt.expandRow( i ); // the editable value rows sit under collapsed category rows
                }
                // INSPECTING: clicking through rows fires writeBack on every selection change. Those write-backs
                // must not be mistaken for edits and must not cost the user an undo step. (They do still call
                // setEdited(true), which clears REDO -- pre-existing, and deliberately left alone: gating the
                // write on the edit flag would risk dropping a real edit if the flag were ever wrong.)
                jt.setSelectionRow( 2 );
                jt.setSelectionRow( 3 );
                if ( "Edit Node Data".equals( tp.undoLabel() ) ) {
                    fail( ok, "clicking through the node editor without editing must not push an undo entry" );
                }
                final String before = tp.getPhylogeny().getFirstExternalNode().getName();
                TreePath name_row = null;
                for( int i = 0; i < jt.getRowCount(); ++i ) {
                    if ( before.equals( jt.getPathForRow( i ).getLastPathComponent().toString() ) ) {
                        name_row = jt.getPathForRow( i );
                        break;
                    }
                }
                if ( name_row == null ) {
                    fail( ok, "could not find the editor row holding the node name" );
                    return;
                }
                // A BARE CLICK is not an edit. editingStopped fires even when nothing was typed -- and this panel
                // opens the inline editor by itself for any empty field -- so opening and closing the editor
                // without changing the text must not count, or reading a node would cost an undo step.
                jt.setSelectionPath( name_row );
                jt.startEditingAtPath( name_row );
                jt.stopEditing();
                jt.setSelectionRow( 0 );
                if ( "Edit Node Data".equals( tp.undoLabel() ) ) {
                    fail( ok, "opening and closing a cell editor WITHOUT changing the text is not an edit" );
                }
                jt.setSelectionPath( name_row );
                jt.startEditingAtPath( name_row );
                // Type into the editor the way a user does: the inline editor is a JTextField inside the tree
                // while editing. Setting the text and stopping the edit is what makes the COMMITTED value differ
                // from the one the edit started at -- which is exactly what the panel watches for.
                final javax.swing.JTextField field = findTextField( jt );
                if ( field == null ) {
                    fail( ok, "the node editor should use an inline text field" );
                    return;
                }
                field.setText( "RENAMED" );
                jt.stopEditing();
                if ( "Edit Node Data".equals( tp.undoLabel() ) ) {
                    fail( ok, "the checkpoint belongs on the WRITE, not on the cell edit itself" );
                }
                jt.setSelectionRow( 0 ); // moving off the row is what commits the field to the phylogeny
                if ( !"RENAMED".equals( tp.getPhylogeny().getFirstExternalNode().getName() ) ) {
                    fail( ok, "precondition: the edit should have reached the tree, got "
                            + tp.getPhylogeny().getFirstExternalNode().getName() );
                }
                if ( !tp.canUndo() || !"Edit Node Data".equals( tp.undoLabel() ) ) {
                    fail( ok, "a node-data edit should checkpoint 'Edit Node Data', got '" + tp.undoLabel() + "'" );
                }
                tp.undo();
                if ( !before.equals( tp.getPhylogeny().getFirstExternalNode().getName() ) ) {
                    fail( ok, "undo should restore the old node name, got "
                            + tp.getPhylogeny().getFirstExternalNode().getName() );
                }
                // An undo installs a DIFFERENT tree, so an editor left open would be editing a node detached from
                // the one on screen -- the edits would vanish silently while still marking the file dirty.
                if ( ( tp.openNodeFrameCountForTest() != 0 ) || nf[ 0 ].isDisplayable() ) {
                    fail( ok, "undo must close the open node editor -- its node belongs to the replaced tree" );
                }
            } );
            // A THEME SWITCH must not turn the editor read-only. setDarkMode runs updateComponentTreeUI over every
            // open window, which reinstalls the tree UI and SWAPS IN A NEW CELL EDITOR; a change-watcher bound only
            // at construction would stop firing from then on.
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final PhylogenyNode n2 = tp.getPhylogeny().getFirstExternalNode();
                tp.showNodeEditFrame( n2 );
                nf[ 0 ] = openNodeFrame();
                if ( nf[ 0 ] == null ) {
                    fail( ok, "the editor should have reopened" );
                    return;
                }
                SwingUtilities.updateComponentTreeUI( nf[ 0 ] ); // what a light/dark switch does to this window
                final JTree jt = findJTree( nf[ 0 ] );
                if ( jt == null ) {
                    fail( ok, "the node editor should still contain a JTree after a UI refresh" );
                    return;
                }
                for( int i = 0; i < jt.getRowCount(); ++i ) {
                    jt.expandRow( i );
                }
                final String before2 = n2.getName();
                TreePath row = null;
                for( int i = 0; i < jt.getRowCount(); ++i ) {
                    if ( before2.equals( jt.getPathForRow( i ).getLastPathComponent().toString() ) ) {
                        row = jt.getPathForRow( i );
                        break;
                    }
                }
                if ( row == null ) {
                    fail( ok, "could not find the name row after the UI refresh" );
                    return;
                }
                jt.setSelectionPath( row );
                jt.startEditingAtPath( row );
                final javax.swing.JTextField field2 = findTextField( jt );
                if ( field2 == null ) {
                    fail( ok, "no inline text field after the UI refresh" );
                    return;
                }
                field2.setText( "AFTER_THEME_SWITCH" );
                jt.stopEditing();
                jt.setSelectionRow( 0 );
                if ( !"AFTER_THEME_SWITCH".equals( tp.getPhylogeny().getFirstExternalNode().getName() ) ) {
                    fail( ok, "a theme switch must not stop the node editor writing edits through, name is "
                            + tp.getPhylogeny().getFirstExternalNode().getName() );
                }
                if ( !"Edit Node Data".equals( tp.undoLabel() ) ) {
                    fail( ok, "an edit after a theme switch should still be undoable, got '" + tp.undoLabel() + "'" );
                }
            } );
            SwingUtilities.invokeAndWait( () -> {
                nf[ 0 ].dispose();
                ( (JFrame) mf[ 0 ] ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** The NodeFrame currently on screen (the panel tracks it internally but does not hand it out). */
    private static NodeFrame openNodeFrame() {
        for( final java.awt.Window w : java.awt.Window.getWindows() ) {
            if ( ( w instanceof NodeFrame ) && w.isDisplayable() ) {
                return (NodeFrame) w;
            }
        }
        return null;
    }

    private static javax.swing.JTextField findTextField( final java.awt.Container c ) {
        for( final java.awt.Component k : c.getComponents() ) {
            if ( k instanceof javax.swing.JTextField ) {
                return (javax.swing.JTextField) k;
            }
            if ( k instanceof java.awt.Container ) {
                final javax.swing.JTextField f = findTextField( (java.awt.Container) k );
                if ( f != null ) {
                    return f;
                }
            }
        }
        return null;
    }

    private static JTree findJTree( final java.awt.Container c ) {
        for( final java.awt.Component k : c.getComponents() ) {
            if ( k instanceof JTree ) {
                return (JTree) k;
            }
            if ( k instanceof java.awt.Container ) {
                final JTree t = findJTree( (java.awt.Container) k );
                if ( t != null ) {
                    return t;
                }
            }
        }
        return null;
    }

    /** The node named {@code name} in the panel's CURRENT tree -- undo/redo replace the tree object, so a node
     *  reference taken before a restore is stale and must be looked up again. */
    private static PhylogenyNode named( final TreePanel tp, final String name ) {
        for( final PhylogenyNodeIterator it = tp.getPhylogeny().iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        throw new IllegalStateException( "no node named " + name );
    }

    /** A tree with a non-root INTERNAL node, which is what collapse() requires (it refuses tips and the root). */
    private static Phylogeny nestedTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final String clade_name : new String[] { "cladeA", "cladeB" } ) {
            final PhylogenyNode clade = new PhylogenyNode();
            clade.setName( clade_name );
            for( int i = 0; i < 2; ++i ) {
                final PhylogenyNode leaf = new PhylogenyNode();
                leaf.setName( clade_name + "_t" + i );
                clade.addAsChild( leaf );
            }
            root.addAsChild( clade );
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
