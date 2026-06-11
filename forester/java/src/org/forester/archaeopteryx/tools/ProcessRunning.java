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

package org.forester.archaeopteryx.tools;

import java.text.SimpleDateFormat;
import java.util.Calendar;

final public class ProcessRunning {

    private static long  count = 0;
    final private long   _id;
    final private String _name;
    final private String _start;

    public long getId() {
        return _id;
    }

    public String getName() {
        return _name;
    }

    public String getStart() {
        return _start;
    }

    @Override
    public String toString() {
        return getName() + " [id=" + getId() + "] [start=" + getStart() + "]";
    }

    synchronized static ProcessRunning createInstance( final String name ) {
        final Calendar cal = Calendar.getInstance();
        final SimpleDateFormat sdf = new SimpleDateFormat( "HH:mm:ss" );
        return new ProcessRunning( count++, name, sdf.format( cal.getTime() ) );
    }

    private ProcessRunning( final long id, final String name, final String start ) {
        if ( id < 0 ) {
            throw new IllegalArgumentException( "process id cannot be negative" );
        }
        _id = id;
        _name = name;
        _start = start;
    }
}
