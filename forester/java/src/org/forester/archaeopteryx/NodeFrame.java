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

import org.forester.phylogeny.PhylogenyNode;

/**
 * The per-node window: a {@link NodeDataForm} in the shared {@link EditorFrame} chrome. {@link NodeDataForm.Mode#VIEW}
 * ("Show node data") is read-only with a single Close button; {@link NodeDataForm.Mode#EDIT} ("Edit node data") adds
 * the status line and the <b>Write to Tree</b> button. The window is non-modal so several can stay open; the tree
 * panel tracks them and closes them all when an undo or redo swaps the tree underneath (their node would belong to
 * the replaced tree).
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

    /** Releases the tree panel's slot. */
    @Override
    protected void onClosed() {
        if ( _tree_panel != null ) {
            _tree_panel.removeEditNodeFrame( this );
        }
    }
}
