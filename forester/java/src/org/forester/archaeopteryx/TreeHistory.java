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

import java.util.ArrayDeque;
import java.util.Deque;

import org.forester.phylogeny.Phylogeny;

/**
 * A bounded per-tree undo/redo history built by <em>snapshotting</em> the whole phylogeny (a deep
 * {@link Phylogeny#copy()}) before each mutating operation — the tree-native, command-free approach: any
 * operation is undoable just by taking a snapshot first, no per-operation inverse needed. This deliberately
 * does NOT use {@code javax.swing.undo.UndoManager} (which is command-based).
 * <p>
 * Ownership rule that keeps the stacks independent of the live tree: {@link #checkpoint} stores a <em>copy</em>
 * of the current tree (the live tree keeps being mutated), and {@link #undo}/{@link #redo} copy the current
 * live tree onto the opposite stack before returning the (de-stacked) snapshot for the caller to install as
 * the new live tree. So a stacked snapshot is never the live tree and is never mutated in place.
 * <p>
 * Each snapshot also records the "edited" (unsaved-changes) flag as it was at capture time, so undoing all the
 * way back to the freshly-loaded state correctly restores a clean (not-edited) marker.
 * <p>
 * Pure and headless-testable: it only reads/copies {@link Phylogeny} objects; the GUI refresh lives in
 * {@code TreePanel}.
 */
final class TreeHistory {

    static final int DEFAULT_MAX_DEPTH = 25;

    /** One recorded tree state: the (copied) tree, the operation label, and the edited flag at capture time. */
    static final class Snapshot {

        private final Phylogeny _phylogeny;
        private final String    _label;
        private final boolean   _edited;

        Snapshot( final Phylogeny phylogeny, final String label, final boolean edited ) {
            _phylogeny = phylogeny;
            _label = label;
            _edited = edited;
        }

        Phylogeny getPhylogeny() {
            return _phylogeny;
        }

        String getLabel() {
            return _label;
        }

        boolean isEdited() {
            return _edited;
        }
    }

    private final int             _max_depth;
    private final Deque<Snapshot> _undo = new ArrayDeque<Snapshot>();
    private final Deque<Snapshot> _redo = new ArrayDeque<Snapshot>();

    TreeHistory() {
        this( DEFAULT_MAX_DEPTH );
    }

    TreeHistory( final int max_depth ) {
        _max_depth = Math.max( 1, max_depth );
    }

    /**
     * Records the current tree state (as it is BEFORE the mutation about to happen) under {@code label}, and
     * clears the redo history (a new action invalidates any redoable future). The oldest entries beyond the
     * depth cap are dropped.
     */
    void checkpoint( final Phylogeny current, final String label, final boolean edited ) {
        if ( ( current == null ) || current.isEmpty() ) {
            return;
        }
        _undo.addFirst( new Snapshot( current.copy(), label, edited ) );
        while ( _undo.size() > _max_depth ) {
            _undo.pollLast();
        }
        _redo.clear();
    }

    boolean canUndo() {
        return !_undo.isEmpty();
    }

    boolean canRedo() {
        return !_redo.isEmpty();
    }

    /** The label of the next undoable operation, or {@code null}. */
    String undoLabel() {
        return _undo.isEmpty() ? null : _undo.peekFirst().getLabel();
    }

    /** The label of the next redoable operation, or {@code null}. */
    String redoLabel() {
        return _redo.isEmpty() ? null : _redo.peekFirst().getLabel();
    }

    /**
     * Undo: pushes a copy of the current live state onto the redo stack, then returns the previous state to
     * install as the new live tree (or {@code null} when there is nothing to undo).
     */
    Snapshot undo( final Phylogeny current, final boolean current_edited ) {
        if ( _undo.isEmpty() ) {
            return null;
        }
        final Snapshot previous = _undo.pollFirst();
        _redo.addFirst( new Snapshot( current.copy(), previous.getLabel(), current_edited ) );
        while ( _redo.size() > _max_depth ) {
            _redo.pollLast();
        }
        return previous;
    }

    /**
     * Redo: pushes a copy of the current live state onto the undo stack, then returns the next state to
     * install as the new live tree (or {@code null} when there is nothing to redo).
     */
    Snapshot redo( final Phylogeny current, final boolean current_edited ) {
        if ( _redo.isEmpty() ) {
            return null;
        }
        final Snapshot next = _redo.pollFirst();
        _undo.addFirst( new Snapshot( current.copy(), next.getLabel(), current_edited ) );
        while ( _undo.size() > _max_depth ) {
            _undo.pollLast();
        }
        return next;
    }

    void clear() {
        _undo.clear();
        _redo.clear();
    }

    /**
     * Clears only the redo history. Called when the tree is edited outside of undo/redo (a mutation that was
     * not checkpointed, or one that came after an undo): any redoable "future" is now stale, and dropping it
     * makes a later Redo impossible rather than letting it install an unrelated tree.
     */
    void clearRedo() {
        _redo.clear();
    }

    int undoDepth() {
        return _undo.size();
    }

    int redoDepth() {
        return _redo.size();
    }
}
