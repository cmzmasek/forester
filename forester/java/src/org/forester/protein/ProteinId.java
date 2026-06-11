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

package org.forester.protein;

import org.forester.util.ForesterUtil;

public class ProteinId implements Comparable<ProteinId> {

    final private String _id;

    public ProteinId( final String id ) {
        if ( ForesterUtil.isEmpty( id ) ) {
            throw new IllegalArgumentException( "attempt to create new protein id from empty or null string" );
        }
        _id = id.trim();
    }

    @Override
    public int compareTo( final ProteinId protein_id ) {
        if ( this == protein_id ) {
            return 0;
        }
        return getId().toLowerCase().compareTo( protein_id.getId().toLowerCase() );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check protein id equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check protein id equality to " + o + " [" + o.getClass()
                                                + "]" );
        }
        else {
            return getId().equals( ( ( ProteinId ) o ).getId() );
        }
    }

    public String getId() {
        return _id;
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
