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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The per-node window: a {@link NodeDataForm} in the shared {@link EditorFrame} chrome. {@link NodeDataForm.Mode#VIEW}
 * ("Show node data") is read-only with a single Close button; {@link NodeDataForm.Mode#EDIT} ("Edit node data") adds
 * the status line and the <b>Write to Tree</b> button. The window is non-modal so several can stay open; the tree
 * panel tracks them and, when an undo or redo swaps the tree underneath, {@link #rebind rebinds} each to its node in
 * the restored tree (unwritten edits are kept).
 */
final class NodeFrame extends EditorFrame {

    private static final long  serialVersionUID = -6943510233968557246L;
    private final TreePanel    _tree_panel;
    private final NodeDataForm _form;

    NodeFrame( final PhylogenyNode n, final TreePanel tp, final int index, final NodeDataForm.Mode mode ) {
        this( new NodeDataForm( n, tp, mode ), n, tp, index );
    }

    private NodeFrame( final NodeDataForm form, final PhylogenyNode n, final TreePanel tp, final int index ) {
        super( form, ( form.isEditable() ? "Edit Node: " : "Node: " ) + NodeDataForm.nodeLabel( n ),
               "This node has changes that have not been written to the tree." );
        _tree_panel = tp;
        _form = form;
        sizeAndPlace( tp, index, 44, 24 ); // the slot index only cascades the window position
        setVisible( true );
    }

    NodeDataForm getForm() {
        return _form;
    }

    /**
     * After an undo or redo installed {@code tree}: re-attaches this window to the node with the same id in it
     * (node ids survive the snapshot copy), keeping any unwritten edits. When no such node exists any more -- the
     * undone step was what added it, or a redo deletes it -- the window stays open, says so in its status line and
     * cannot write, so nothing typed disappears unseen; a later redo/undo that brings the node back re-attaches it.
     */
    void rebind( final Phylogeny tree ) {
        final long id = _form.node().getId();
        final PhylogenyNode n = ( ( tree == null ) || tree.isEmpty() ) ? null : tree.getNode( id );
        final String prefix = _form.isEditable() ? "Edit Node: " : "Node: ";
        if ( n == _form.node() ) { // the live tree changed around the node (a delete elsewhere): nothing to re-read
            _form.refreshHeader();
            setBaseTitle( prefix + NodeDataForm.nodeLabel( n ) );
        }
        else if ( n != null ) {
            _form.rebind( n );
            setBaseTitle( prefix + NodeDataForm.nodeLabel( n ) );
        }
        else {
            _form.markDetached();
            setBaseTitle( prefix + NodeDataForm.nodeLabel( _form.node() ) + " (no longer in the tree)" );
        }
    }

    /** Releases the tree panel's slot. */
    @Override
    protected void onClosed() {
        if ( _tree_panel != null ) {
            _tree_panel.removeEditNodeFrame( this );
        }
    }
}
