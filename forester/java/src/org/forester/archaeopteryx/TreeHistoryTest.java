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

import org.forester.archaeopteryx.TreeHistory.Snapshot;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Unit tests for {@link TreeHistory} (the pure snapshot-based undo/redo stack). Same package as the class
 * under test (package-private). Headless; run via the suite or {@link #main(String[])}.
 */
public final class TreeHistoryTest {

    public static void main( final String[] args ) {
        System.out.println( "TreeHistory: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        return roundTrip() && redoClearedOnCheckpoint() && capEviction() && snapshotIsIndependentCopy()
                && clearRedoOnly() && emptyAndLabels();
    }

    // ---- clearRedo() drops the redo history but leaves the undo stack intact (the safety-net primitive) ----
    private static boolean clearRedoOnly() {
        final TreeHistory h = new TreeHistory();
        h.checkpoint( tree( 2 ), "op1", false );
        h.checkpoint( tree( 2 ), "op2", false ); // undo = [op2, op1]
        h.undo( tree( 3 ), true );               // undo = [op1], redo = [op2-post]
        if ( !h.canRedo() || !h.canUndo() ) {
            return fail( "clearRedo setup: both undo and redo should be available" );
        }
        h.clearRedo();
        if ( h.canRedo() ) {
            return fail( "clearRedo() should empty the redo stack" );
        }
        if ( !h.canUndo() || !"op1".equals( h.undoLabel() ) ) {
            return fail( "clearRedo() must leave the undo stack intact" );
        }
        return true;
    }

    // ---- undo restores the previous state; redo re-applies it; labels + edited flag travel ----
    private static boolean roundTrip() {
        final TreeHistory h = new TreeHistory();
        final Phylogeny before = tree( 2 );
        h.checkpoint( before, "Delete Nodes", false ); // capture the pre-mutation state
        final Phylogeny after = tree( 3 ); // simulate: the mutation replaced the live tree
        if ( !h.canUndo() || h.canRedo() ) {
            return fail( "after one checkpoint: canUndo true, canRedo false" );
        }
        if ( !"Delete Nodes".equals( h.undoLabel() ) ) {
            return fail( "undo label" );
        }
        final Snapshot u = h.undo( after, true );
        if ( ( u == null ) || ( u.getPhylogeny().getNumberOfExternalNodes() != 2 ) || !"Delete Nodes".equals( u.getLabel() )
                || u.isEdited() ) {
            return fail( "undo should restore the 2-tip pre-state with its label and not-edited flag" );
        }
        if ( !h.canRedo() || !"Delete Nodes".equals( h.redoLabel() ) ) {
            return fail( "after undo: canRedo true with the same label" );
        }
        // caller installed u.getPhylogeny() as the live tree; now redo
        final Snapshot r = h.redo( u.getPhylogeny(), false );
        if ( ( r == null ) || ( r.getPhylogeny().getNumberOfExternalNodes() != 3 ) || r.isEdited() != true ) {
            return fail( "redo should restore the 3-tip post-state with edited=true" );
        }
        if ( !h.canUndo() || h.canRedo() ) {
            return fail( "after redo: canUndo true, canRedo false" );
        }
        return true;
    }

    // ---- a new checkpoint discards the redo history ----
    private static boolean redoClearedOnCheckpoint() {
        final TreeHistory h = new TreeHistory();
        h.checkpoint( tree( 2 ), "op1", false );
        h.undo( tree( 3 ), true );
        if ( !h.canRedo() ) {
            return fail( "redo should be available right after an undo" );
        }
        h.checkpoint( tree( 4 ), "op2", true );
        if ( h.canRedo() ) {
            return fail( "a new checkpoint must clear the redo history" );
        }
        if ( !"op2".equals( h.undoLabel() ) ) {
            return fail( "undo label after the new checkpoint" );
        }
        return true;
    }

    // ---- the depth cap drops the oldest states ----
    private static boolean capEviction() {
        final TreeHistory h = new TreeHistory( 2 );
        h.checkpoint( tree( 2 ), "op1", false );
        h.checkpoint( tree( 2 ), "op2", false );
        h.checkpoint( tree( 2 ), "op3", false );
        if ( h.undoDepth() != 2 ) {
            return fail( "depth cap 2 should retain only 2 undo entries, got " + h.undoDepth() );
        }
        if ( !"op3".equals( h.undoLabel() ) ) {
            return fail( "newest entry should be on top" );
        }
        // only two undos are possible (op3, op2); op1 was evicted
        if ( ( h.undo( tree( 2 ), false ) == null ) || ( h.undo( tree( 2 ), false ) == null ) ) {
            return fail( "two undos should succeed" );
        }
        if ( h.undo( tree( 2 ), false ) != null ) {
            return fail( "the evicted third undo must return null" );
        }
        return true;
    }

    // ---- checkpoint stores an independent COPY, so later mutating the live tree can't corrupt it ----
    private static boolean snapshotIsIndependentCopy() {
        final TreeHistory h = new TreeHistory();
        final Phylogeny live = tree( 2 );
        h.checkpoint( live, "op1", false );
        // mutate the live tree AFTER the checkpoint
        final PhylogenyNode extra = new PhylogenyNode();
        extra.setName( "extra" );
        live.getRoot().addAsChild( extra );
        live.externalNodesHaveChanged();
        if ( live.getNumberOfExternalNodes() != 3 ) {
            return fail( "test setup: mutated live tree should have 3 tips" );
        }
        final Snapshot u = h.undo( live, true );
        if ( u.getPhylogeny().getNumberOfExternalNodes() != 2 ) {
            return fail( "the snapshot must preserve the pre-mutation 2-tip state (independent copy)" );
        }
        return true;
    }

    // ---- empty history + clear ----
    private static boolean emptyAndLabels() {
        final TreeHistory h = new TreeHistory();
        if ( h.canUndo() || h.canRedo() || ( h.undoLabel() != null ) || ( h.redoLabel() != null ) ) {
            return fail( "a fresh history is empty" );
        }
        if ( ( h.undo( tree( 2 ), false ) != null ) || ( h.redo( tree( 2 ), false ) != null ) ) {
            return fail( "undo/redo on an empty history return null" );
        }
        h.checkpoint( tree( 2 ), "op1", false );
        h.clear();
        if ( h.canUndo() || h.canRedo() ) {
            return fail( "clear() empties both stacks" );
        }
        // a null / empty tree checkpoint is a no-op (must not throw)
        h.checkpoint( null, "op", false );
        h.checkpoint( new Phylogeny(), "op", false );
        if ( h.canUndo() ) {
            return fail( "checkpoint of a null/empty tree should be a no-op" );
        }
        return true;
    }

    private static Phylogeny tree( final int n ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < n; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "L" + i );
            root.addAsChild( leaf );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [TreeHistoryTest] " + message );
        return false;
    }

    private TreeHistoryTest() {
        // not instantiable
    }
}
