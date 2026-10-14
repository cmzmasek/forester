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
 * Where "Annotate Clades by Rank" and "Colorize Subtrees via Taxonomic Rank" take each tip's taxon at a rank from.
 * The user picks ONE source up front, so the figure never depends on an invisible state (what an earlier session
 * happened to look up) and the rank chooser's coverage is the coverage of THIS source.
 * <ul>
 * <li>{@link #TREE}: the tree's own annotations -- a tip's {@code <taxonomy>} at the rank, a per-rank metadata
 *     field, else the nearest ancestor's annotation. Offline; nothing is ever looked up.</li>
 * <li>{@link #DATABASE}: NCBI and UniProt for every tip; a tip's OWN identity (its annotation at the rank, else its
 *     database lineage) overrides an ancestor's annotation, which only fills the gaps.</li>
 * <li>{@link #TREE_THEN_DATABASE}: the tree first; the databases only for the tips the tree cannot place.</li>
 * </ul>
 */
enum TaxonSource {

    TREE( "The tree's own annotations (offline)", false ),
    DATABASE( "NCBI and UniProt for every tip; a tip's own identity overrides the tree (online)", true ),
    TREE_THEN_DATABASE( "The tree first; NCBI and UniProt only for the tips it cannot place (online)", true );

    private final String  _label;
    private final boolean _online;

    TaxonSource( final String label, final boolean online ) {
        _label = label;
        _online = online;
    }

    /** The radio-button text. */
    String getLabel() {
        return _label;
    }

    /** True when this source may look names up (and so may need a network connection). */
    boolean isOnline() {
        return _online;
    }

    /** The source a dialog opens on: the tree's own annotations when the tree carries any rank, else the
     *  databases -- never a network source for a tree that can be drawn offline. */
    static TaxonSource defaultFor( final boolean tree_has_ranks ) {
        return tree_has_ranks ? TREE : DATABASE;
    }
}
