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

import java.awt.Color;

import org.forester.phylogeny.PhylogenyNode;

/**
 * One maximal-monophyletic clade to annotate with a shaded box or a right-edge bar: its taxon (at the
 * chosen rank), its distinct color, and the clade-root node. The clade's tips (the root's external
 * descendants) give the band's vertical extent at paint time, so a {@link CladeBand} survives zoom and
 * scroll. Produced by {@link TreePanelUtil#cladeBands} from the SAME tip&rarr;taxon assignment the rank
 * colorizer uses, so paraphyletic groups become several same-colored bands.
 */
final class CladeBand {

    private final String        _taxon;
    private final Color         _color;
    private final PhylogenyNode _root;

    CladeBand( final String taxon, final Color color, final PhylogenyNode root ) {
        _taxon = taxon;
        _color = color;
        _root = root;
    }

    String getTaxon() {
        return _taxon;
    }

    Color getColor() {
        return _color;
    }

    /** The clade's root node; its external descendants define the band's vertical span. */
    PhylogenyNode getRoot() {
        return _root;
    }
}
