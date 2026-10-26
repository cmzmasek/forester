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

/**
 * The entries of the node menu -- the menu a click on a node opens (left, right or ctrl-click alike, see
 * {@code TreePanel.mouseClicked}), in the order they are listed. Each entry acts ONCE on the clicked node and
 * nothing is remembered: there is no "click-to" mode any more (Christian, 2026-10-26: "No more old-skool 'Click
 * on Node to:' menu, no more 'status' of click behavior"), as in Archaeopteryx.js, whose node menu this mirrors.
 * <p>
 * Entries are grouped ({@link Group}; a separator between groups), and an entry that does not apply to the
 * clicked node is GREYED, never hidden -- with the reason as its tooltip where there is one (Root/Reroot). The
 * exception is {@link #editableOnly()}: the editing entries exist only when the tree is editable
 * ({@code Options.isEditable()}), as they always did. The {@link #title()} is the entry's wording; a few entries
 * read differently by node ("Return to Supertree" on a subtree view's root, "Delete Node…" on a tip), see
 * {@code TreePanel.nodeMenuLabel}.
 */
public enum NodeMenuEntry {

    DISPLAY_NODE_DATA( "Display Node Data", Group.DATA, false ),
    EDIT_NODE_DATA( "Edit Node Data", Group.DATA, true ),
    NODE_STYLE( "Edit Node Style", Group.DATA, false ),
    SUBTREE( "Switch to Subtree", Group.SHAPE, false ),
    COLLAPSE_UNCOLLAPSE( "Collapse/Uncollapse", Group.SHAPE, false ),
    UNCOLLAPSE_ALL( "Uncollapse Subtree", Group.SHAPE, false ),
    SWAP( "Swap Descendants", Group.SHAPE, false ),
    ORDER_SUBTREE( "Ladderize Subtree", Group.SHAPE, false ),
    COLOR_SUBTREE( "Colorize Subtree", Group.SHAPE, false ),
    REROOT( "Root/Reroot", Group.SHAPE, false ),
    SELECT_NODE( "Select/Deselect Node", Group.SELECTION, false ),
    SELECT_TIPS( "Select/Deselect Tips", Group.SELECTION, false ),
    OPEN_SEQ_WEB( "Open Sequence DB", Group.LINKS, false ),
    OPEN_TAX_WEB( "Open Taxonomy DB", Group.LINKS, false ),
    CUT_SUBTREE( "Cut Subtree", Group.EDITING, true ),
    COPY_SUBTREE( "Copy Subtree", Group.EDITING, true ),
    PASTE_SUBTREE( "Paste Subtree", Group.EDITING, true ),
    ADD_NEW_NODE( "Add New Node", Group.EDITING, true ),
    DELETE( "Delete Node/Subtree", Group.DELETE, true );

    /** The menu's groups, in menu order; a separator is drawn between two groups that both have an entry. */
    public enum Group {
        /** the node's data: show, edit, style */
        DATA,
        /** the tree's shape around the node: subtree view, collapse, swap, ladderize, colorize, re-root */
        SHAPE,
        /** the manual selection (found set 0) */
        SELECTION,
        /** web links for the node's accession / taxonomy */
        LINKS,
        /** tree editing: cut, copy, paste, add */
        EDITING,
        /** deletion -- last, and drawn in the warning colour */
        DELETE
    }

    private final String  _title;
    private final Group   _group;
    private final boolean _editable_only;

    NodeMenuEntry( final String title, final Group group, final boolean editable_only ) {
        _title = title;
        _group = group;
        _editable_only = editable_only;
    }

    /** The entry's wording (the label when no node-specific label applies). */
    public String title() {
        return _title;
    }

    /** The group the entry is listed in. */
    public Group group() {
        return _group;
    }

    /** Whether the entry exists only on an editable tree ({@code Options.isEditable()}). */
    public boolean editableOnly() {
        return _editable_only;
    }
}
