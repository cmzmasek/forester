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

package org.forester.phylogeny.iterators;

import org.forester.phylogeny.PhylogenyNode;

/*
 * @author Christian M. Zmasek
 *
 * @version 1.00 -- last modified: 06/15/00
 */
public final class PostOrderStackObject {

    final private PhylogenyNode _node;
    final private int           _phase;

    public PostOrderStackObject( final PhylogenyNode n, final int i ) {
        _node = n;
        _phase = i;
    }

    final public PhylogenyNode getNode() {
        return _node;
    }

    final public int getPhase() {
        return _phase;
    }
}
