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

package org.forester.surfacing;

import java.util.SortedSet;

final class UniprotNameResult {

    final String            _name;
    final Long              _count;
    final int               _out_of;
    final SortedSet<String> _all_names;
    public UniprotNameResult( final String name,
                              final Long count,
                              final int out_of,
                              final SortedSet<String> all_names ) {
        _name = name;
        _count = count;
        _out_of = out_of;
        _all_names = all_names;
    }

    public String getName() {
        return _name;
    }

    public Long getCount() {
        return _count;
    }

    public int getOutOf() {
        return _out_of;
    }

    public SortedSet<String> getAllNames() {
        return _all_names;
    }
}
