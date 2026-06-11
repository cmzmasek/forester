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

import java.util.HashMap;
import java.util.Map;

import org.forester.protein.BasicDomain;
import org.forester.protein.BinaryDomainCombination;

final class DirectedBinaryDomainCombination extends BasicBinaryDomainCombination {

    final private static Map<Integer, DirectedBinaryDomainCombination> DDC_POOL = new HashMap<Integer, DirectedBinaryDomainCombination>();

    private DirectedBinaryDomainCombination( final String n_terminal, final String c_terminal ) {
        super();
        if ( ( n_terminal == null ) || ( c_terminal == null ) ) {
            throw new IllegalArgumentException( "attempt to create binary domain combination using null" );
        }
        _id0 = BasicDomain.obtainIdAsShort( n_terminal );
        _id1 = BasicDomain.obtainIdAsShort( c_terminal );
    }

    public final static BinaryDomainCombination obtainInstance( final String ids ) {
        if ( ids.indexOf( BinaryDomainCombination.SEPARATOR ) < 1 ) {
            throw new IllegalArgumentException( "Unexpected format for binary domain combination [" + ids + "]" );
        }
        final String[] ids_ary = ids.split( BinaryDomainCombination.SEPARATOR );
        if ( ids_ary.length != 2 ) {
            throw new IllegalArgumentException( "Unexpected format for binary domain combination [" + ids + "]" );
        }
        return DirectedBinaryDomainCombination.obtainInstance( ids_ary[ 0 ], ids_ary[ 1 ] );
    }

    public final static DirectedBinaryDomainCombination obtainInstance( final String n_terminal, final String c_terminal ) {
        final int code = calcCode( BasicDomain.obtainIdAsShort( n_terminal ), BasicDomain.obtainIdAsShort( c_terminal ) );
        if ( DDC_POOL.containsKey( code ) ) {
            return DDC_POOL.get( code );
        }
        else {
            final DirectedBinaryDomainCombination dc = new DirectedBinaryDomainCombination( n_terminal, c_terminal );
            DDC_POOL.put( code, dc );
            if ( VERBOSE && ( ( DDC_POOL.size() % 100 ) == 0 ) ) {
                System.out.println( " ddc pool size: " + DDC_POOL.size() );
            }
            return dc;
        }
    }
}
