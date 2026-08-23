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
 * The "click on a node" actions offered in the control panel's click-to dropdown. Replaces the old
 * positional {@code String[][] Configuration.clickto_options} array and its parallel
 * {@code final static int} index constants with a type-safe enum carrying each action's menu title.
 * All click-to actions are always offered (the old "display"/"nodisplay" second column was uniformly
 * "display"); the editable-only subset (cut/copy/paste/delete/add/edit) is still gated separately by
 * {@code Options.isEditable()} at the call site.
 */
public enum ClickToOption {

    DISPLAY_NODE_DATA( "Display Node Data" ),
    SUBTREE( "Switch to Sub/Supertree" ),
    SELECT_NODES( "Select Node/Select Tips" ),
    COLLAPSE_UNCOLLAPSE( "Collapse/Uncollapse" ),
    UNCOLLAPSE_ALL( "Uncollapse Subtree" ),
    SWAP( "Swap Descendants" ),
    ORDER_SUBTREE( "Ladderize Subtree" ),
    OPEN_TAX_WEB( "Open Taxonomy DB" ),
    OPEN_SEQ_WEB( "Open Sequence DB" ),
    NODE_STYLE( "Edit Node Style" ),
    COLOR_SUBTREE( "Colorize Subtree" ),
    EDIT_NODE_DATA( "Edit Node Data" ),
    REROOT( "Root/Reroot" ),
    DELETE_SUBTREE_OR_NODE( "Delete Node/Subtree" ),
    ADD_NEW_NODE( "Add New Node" ),
    CUT_SUBTREE( "Cut Subtree" ),
    COPY_SUBTREE( "Copy Subtree" ),
    PASTE_SUBTREE( "Paste Subtree" );

    private final String _title;

    ClickToOption( final String title ) {
        _title = title;
    }

    /** The dropdown / popup-menu label. */
    public String title() {
        return _title;
    }
}
