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

package org.forester.go;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GoId implements Comparable<GoId> {

    private static final int     SIZE       = 10;
    private static final String  GO_PREFIX  = "GO:";
    private static final String  GO_FORMAT  = GO_PREFIX + "\\d{7}";
    private static final Pattern GO_PATTERN = Pattern.compile( GO_FORMAT );
    private final String         _id;

    public GoId( final String id ) {
        if ( id.length() != SIZE ) {
            throw new IllegalArgumentException( "unexpected format for GO id: " + id );
        }
        final Matcher m = GO_PATTERN.matcher( id );
        if ( !m.matches() ) {
            throw new IllegalArgumentException( "unexpected format for GO id: " + id );
        }
        _id = id.substring( 3 );
    }

    @Override
    public int compareTo( final GoId id ) {
        return getId().compareTo( id.getId() );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check go id equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check go id equality to " + o + " [" + o.getClass() + "]" );
        }
        else {
            return getId().equals( ( ( GoId ) o ).getId() );
        }
    }

    public String getId() {
        return GO_PREFIX + _id;
    }

    @Override
    public int hashCode() {
        return getId().hashCode();
    }

    @Override
    public String toString() {
        return getId();
    }
}
