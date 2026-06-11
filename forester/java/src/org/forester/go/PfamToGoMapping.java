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

public class PfamToGoMapping implements Mapping {

    private final String _pfam_domain_id;
    private final GoId   _go_id;

    public PfamToGoMapping( final String pfam_domain_id, final GoId go_id ) {
        _pfam_domain_id = pfam_domain_id;
        _go_id = go_id;
    }

    @Override
    public int compareTo( final Mapping m ) {
        if ( this == m ) {
            return 0;
        }
        return getKey().compareTo( ( String ) m.getKey() );
    }

    /**
     * Based on key and value.
     *
     *
     */
    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check pfam to go mapping equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check pfam to go mapping equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return getKey().equals( ( ( PfamToGoMapping ) o ).getKey() )
                    && getValue().equals( ( ( PfamToGoMapping ) o ).getValue() );
        }
    }

    @Override
    public String getKey() {
        return _pfam_domain_id;
    }

    @Override
    public GoId getValue() {
        return _go_id;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( getKey().toString() );
        sb.append( " > " );
        sb.append( getValue().toString() );
        return sb.toString();
    }
}
