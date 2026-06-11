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

package org.forester.util;

import java.util.Comparator;

public class StringInt {

    private final String _s;
    private final int    _i;

    public StringInt( final String s, final int i ) {
        _s = s;
        _i = i;
    }

    public String getString() {
        return _s;
    }

    public int getInt() {
        return _i;
    }

    public static final class DescendingIntComparator implements Comparator<StringInt> {

        @Override
        public int compare( final StringInt o1, final StringInt o2 ) {
            if ( o1.getInt() > o2.getInt() ) {
                return -1;
            }
            if ( o1.getInt() < o2.getInt() ) {
                return 1;
            }
            return 0;
        }
    }
}
