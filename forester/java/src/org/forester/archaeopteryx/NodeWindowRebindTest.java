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
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Open node windows SURVIVE an undo or redo: each re-attaches to its node in the restored tree (by node id), a
 * clean one reloads from it, a dirty one keeps its unwritten edits (now measured against the restored node) and
 * writes them into the tree on screen, and a window whose node is gone from the tree stays open with a notice
 * instead of quietly discarding what was typed. Before this, every undo closed every node window.
 */
public final class NodeWindowRebindTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeWindowRebind: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // needs a display
        }
        try {
            return cleanAndDirtySurviveUndo() && sequenceCardsRebound() && goneNodeDetaches() && viewAndKindChange();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** A clean editor reloads from the restored node; a dirty one keeps its edit and writes it into the live tree. */
    private static boolean cleanAndDirtySurviveUndo() throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { tree() }, new Configuration(), "rebind1" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            final PhylogenyNode a = named( tp, "A" );
            final long id_a = a.getId();
            tp.showNodeEditFrame( a );
            final NodeFrame fr = openFrame( 0 );
            check( ok, "editor open", ( fr != null ) && ( tp.openNodeFrameCountForTest() == 1 ) );
            // a Write is one undo step; the undo installs a COPY of the tree, so "A" is a different object now
            fr.getForm().setTextForTest( NodeDataDraft.NAME, "A2" );
            check( ok, "write A2", fr.writeNow() && "A2".equals( named( tp, "A2" ).getName() ) );
            tp.undo();
            check( ok, "undo restored the name", "A".equals( tp.getPhylogeny().getNode( id_a ).getName() ) );
            check( ok, "window still open and tracked",
                   fr.isDisplayable() && ( tp.openNodeFrameCountForTest() == 1 ) && tp.isNodeFrameTrackedForTest( fr ) );
            check( ok, "window now holds the RESTORED tree's node",
                   fr.getForm().node() == tp.getPhylogeny().getNode( id_a ) );
            check( ok, "clean window reloaded the field: " + valueOf( fr, NodeDataDraft.NAME ),
                   "A".equals( valueOf( fr, NodeDataDraft.NAME ) ) );
            check( ok, "not dirty after reload", !fr.isDirty() );
            check( ok, "title follows the restored node: " + fr.getTitle(), "Edit Node: A".equals( fr.getTitle() ) );
            check( ok, "no notice", fr.getForm().notice() == null );
            tp.redo();
            check( ok, "redo reloads too: " + valueOf( fr, NodeDataDraft.NAME ),
                   "A2".equals( valueOf( fr, NodeDataDraft.NAME ) ) && !fr.isDirty() );
            // -- dirty: type, do NOT write, then undo --
            fr.getForm().setTextForTest( NodeDataDraft.NAME, "A3" );
            check( ok, "dirty before undo", fr.isDirty() );
            tp.undo();
            check( ok, "tree back to A", "A".equals( tp.getPhylogeny().getNode( id_a ).getName() ) );
            check( ok, "unwritten edit KEPT across the undo: " + valueOf( fr, NodeDataDraft.NAME ),
                   "A3".equals( valueOf( fr, NodeDataDraft.NAME ) ) );
            check( ok, "still dirty (vs the restored node)", fr.isDirty() && fr.getTitle().startsWith( "• " ) );
            check( ok, "status says unsaved", "Unsaved changes".equals( fr.statusTextForTest() ) );
            check( ok, "write enabled", fr.writeButtonForTest().isEnabled() );
            check( ok, "header shows the restored label", "A".equals( headerTitle( fr ) ) );
            check( ok, "write goes to the tree ON SCREEN", fr.writeNow()
                    && "A3".equals( tp.getPhylogeny().getNode( id_a ).getName() ) );
            check( ok, "that write is an undo step", "Edit Node Data".equals( tp.undoLabel() ) );
            check( ok, "clean after write", !fr.isDirty() );
            fr.close();
            ( (JFrame) mf[ 0 ] ).dispose();
        } );
        return ok[ 0 ];
    }

    /** A dirty sequence card is re-bound to the restored node's sequence OBJECT, so the write mutates it in
     *  place (its annotations survive) instead of dragging the replaced tree's sequence into the live one. */
    private static boolean sequenceCardsRebound() throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { tree() }, new Configuration(), "rebind2" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            final PhylogenyNode a = named( tp, "A" );
            final long id_a = a.getId();
            tp.showNodeEditFrame( a );
            final NodeFrame fr = openFrame( 0 );
            check( ok, "one sequence card", fr.getForm().sequenceCardCountForTest() == 1 );
            // an unrelated undoable step, so an undo installs a copy WITHOUT touching A's sequence
            tp.collapse( named( tp, "cladeC" ) );
            fr.getForm().setTextForTest( "sequence.0." + NodeDataDraft.SEQ_NAME, "seqA-renamed" );
            check( ok, "dirty", fr.isDirty() );
            tp.undo();
            final Sequence restored = tp.getPhylogeny().getNode( id_a ).getNodeData().getSequence();
            check( ok, "restored node has its own sequence object", ( restored != null ) && ( restored != a.getNodeData()
                    .getSequence() ) );
            check( ok, "edit kept", "seqA-renamed".equals( valueOf( fr, "sequence.0." + NodeDataDraft.SEQ_NAME ) ) );
            check( ok, "write ok", fr.writeNow() );
            final Sequence after = tp.getPhylogeny().getNode( id_a ).getNodeData().getSequence();
            check( ok, "the RESTORED sequence object was mutated in place", after == restored );
            check( ok, "renamed", "seqA-renamed".equals( after.getName() ) );
            check( ok, "its annotation survived", ( after.getAnnotations() != null )
                    && ( after.getAnnotations().size() == 1 ) && "GO:0005634".equals( after.getAnnotation( 0 ).getRef() ) );
            fr.close();
            ( (JFrame) mf[ 0 ] ).dispose();
        } );
        return ok[ 0 ];
    }

    /** A window whose node is deleted (in place, or by a redo) stays open, says so, cannot write, and re-attaches
     *  -- edits intact -- when an undo brings the node back. */
    private static boolean goneNodeDetaches() throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { tree() }, new Configuration(), "rebind3" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            final PhylogenyNode d = named( tp, "D" );
            final long id_d = d.getId();
            tp.showNodeEditFrame( d );
            final NodeFrame fr = openFrame( 0 );
            fr.getForm().setTextForTest( NodeDataDraft.NAME, "D9" );
            check( ok, "dirty", fr.isDirty() );
            tp.deleteNodeOrSubtreeConfirmed( d, true ); // in place: the live tree loses D
            check( ok, "D gone", tp.getPhylogeny().getNode( id_d ) == null );
            check( ok, "window stays open", fr.isDisplayable() && ( tp.openNodeFrameCountForTest() == 1 ) );
            final String notice = fr.getForm().notice();
            check( ok, "notice set", ( notice != null ) && notice.contains( "no longer in the tree" ) );
            check( ok, "status shows the notice: " + fr.statusTextForTest(), notice.equals( fr.statusTextForTest() ) );
            check( ok, "cannot write", !fr.writeButtonForTest().isEnabled() && !fr.writeNow() );
            check( ok, "not dirty while detached (closing must not offer a Write)", !fr.isDirty() );
            check( ok, "title says so: " + fr.getTitle(), fr.getTitle().endsWith( "(no longer in the tree)" ) );
            check( ok, "typed value still visible", "D9".equals( valueOf( fr, NodeDataDraft.NAME ) ) );
            tp.undo(); // D is back (a copy)
            check( ok, "re-attached", fr.getForm().notice() == null
                    && ( fr.getForm().node() == tp.getPhylogeny().getNode( id_d ) ) );
            check( ok, "edit survived the round trip", "D9".equals( valueOf( fr, NodeDataDraft.NAME ) ) && fr.isDirty() );
            check( ok, "write enabled again", fr.writeButtonForTest().isEnabled() );
            check( ok, "plain title again: " + fr.getTitle(), "• Edit Node: D".equals( fr.getTitle() ) );
            tp.redo(); // D deleted again
            check( ok, "detached again", ( fr.getForm().notice() != null ) && !fr.isDirty() );
            tp.undo();
            check( ok, "and back, edit still there", fr.isDirty() && "D9".equals( valueOf( fr, NodeDataDraft.NAME ) ) );
            check( ok, "writes into the live tree", fr.writeNow()
                    && "D9".equals( tp.getPhylogeny().getNode( id_d ).getName() ) );
            fr.close();
            ( (JFrame) mf[ 0 ] ).dispose();
        } );
        return ok[ 0 ];
    }

    /** A read-only window rebuilds from the restored node; an editor whose node changed kind (a tip that became
     *  internal when a child was added, then undone) keeps its edit through the page rebuild. */
    private static boolean viewAndKindChange() throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { tree() }, new Configuration(), "rebind4" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            final long id_b = named( tp, "B" ).getId();
            // -- VIEW window on B, then rename B through an editor and undo: the view follows --
            tp.showNodeFrameForTest( named( tp, "B" ) );
            final NodeFrame view = openFrame( 0 );
            check( ok, "view open", ( view != null ) && !view.isEditable() );
            tp.showNodeEditFrame( named( tp, "B" ) );
            final NodeFrame edit = openFrame( 1 );
            edit.getForm().setTextForTest( NodeDataDraft.NAME, "B2" );
            check( ok, "write B2", edit.writeNow() );
            tp.undo();
            check( ok, "view still open", view.isDisplayable() && ( tp.openNodeFrameCountForTest() == 2 ) );
            check( ok, "view rebuilt from the restored node: " + view.getTitle(), "Node: B".equals( view.getTitle() )
                    && "B".equals( headerTitle( view ) ) && ( view.getForm().node() == tp.getPhylogeny().getNode( id_b ) ) );
            check( ok, "view never dirty", !view.isDirty() && ( view.getForm().notice() == null ) );
            view.close();
            // -- kind change: B (a tip) gets a child, so it is internal; the editor is dirty; undo makes it a tip --
            edit.getForm().setTextForTest( NodeDataDraft.NAME, "B3" );
            tp.pushUndoCheckpoint( "Add Node" ); // what the click-to Add Node does, minus its dialog
            final PhylogenyNode child = new PhylogenyNode();
            child.setName( "B_child" );
            tp.getPhylogeny().getNode( id_b ).addAsChild( child );
            tp.afterTreeStructureChanged();
            check( ok, "B is internal now", !tp.getPhylogeny().getNode( id_b ).isExternal() );
            check( ok, "edit kept while the tree changed around the node",
                   "B3".equals( valueOf( edit, NodeDataDraft.NAME ) ) && edit.isDirty() );
            tp.undo();
            check( ok, "B is a tip again", tp.getPhylogeny().getNode( id_b ).isExternal() );
            check( ok, "edit kept through the kind change: " + valueOf( edit, NodeDataDraft.NAME ),
                   "B3".equals( valueOf( edit, NodeDataDraft.NAME ) ) && edit.isDirty() );
            check( ok, "no events section on a tip", !edit.getForm().hasSectionForTest( NodeDataDraft.SEC_EVENTS ) );
            check( ok, "writes", edit.writeNow() && "B3".equals( tp.getPhylogeny().getNode( id_b ).getName() ) );
            edit.close();
            ( (JFrame) mf[ 0 ] ).dispose();
        } );
        return ok[ 0 ];
    }

    // ---- fixtures / helpers ----

    /** (A:0.1, B:0.2, (C:0.1, D:0.1)cladeC:0.3); A carries a sequence with one annotation. */
    private static Phylogeny tree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = tip( "A", 0.1 );
        final Sequence seq = new Sequence();
        seq.setName( "seqA" );
        seq.addAnnotation( new Annotation( "GO:0005634" ) );
        a.getNodeData().setSequence( seq );
        root.addAsChild( a );
        root.addAsChild( tip( "B", 0.2 ) );
        final PhylogenyNode clade = new PhylogenyNode();
        clade.setName( "cladeC" );
        clade.setDistanceToParent( 0.3 );
        clade.addAsChild( tip( "C", 0.1 ) );
        clade.addAsChild( tip( "D", 0.1 ) );
        root.addAsChild( clade );
        phy.setRoot( root );
        phy.setRooted( true );
        return phy;
    }

    private static PhylogenyNode tip( final String name, final double bl ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( bl );
        return n;
    }

    private static PhylogenyNode named( final TreePanel tp, final String name ) {
        for( final PhylogenyNodeIterator it = tp.getPhylogeny().iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        throw new IllegalStateException( "no node named " + name );
    }

    /** The i-th node window currently on screen (in opening order). */
    private static NodeFrame openFrame( final int i ) {
        final List<NodeFrame> frames = new ArrayList<>();
        for( final java.awt.Window w : java.awt.Window.getWindows() ) {
            if ( ( w instanceof NodeFrame ) && w.isDisplayable() ) {
                frames.add( (NodeFrame) w );
            }
        }
        return ( i < frames.size() ) ? frames.get( i ) : null;
    }

    private static String valueOf( final NodeFrame fr, final String key ) {
        final javax.swing.JComponent c = fr.getForm().fieldForTest( key );
        return ( c instanceof javax.swing.text.JTextComponent ) ? ( (javax.swing.text.JTextComponent) c ).getText()
                : null;
    }

    private static String headerTitle( final NodeFrame fr ) {
        return fr.getForm().headerTitleForTest();
    }

    private static void check( final boolean[] ok, final String what, final boolean condition ) {
        if ( !condition ) {
            System.out.println( "  [NodeWindowRebindTest] " + what );
            ok[ 0 ] = false;
        }
    }

    private NodeWindowRebindTest() {
        // not instantiable
    }
}
