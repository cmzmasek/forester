// $Id:
// Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.surfacing;

import java.util.HashMap;
import java.util.Map;

import org.forester.protein.BasicDomain;
import org.forester.protein.BinaryDomainCombination;
import org.forester.util.ForesterUtil;

public class BasicBinaryDomainCombination implements BinaryDomainCombination {

    final static boolean                                            VERBOSE = false;
    final private static Map<Integer, BasicBinaryDomainCombination> DC_POOL = new HashMap<Integer, BasicBinaryDomainCombination>();
    final private static Map<Integer, String>                       S_POOL  = new HashMap<Integer, String>();
    short                                                           _id0;
    short                                                           _id1;

    BasicBinaryDomainCombination() {
        _id0 = -1;
        _id1 = -1;
    }

    private BasicBinaryDomainCombination( final String id0, final String id1 ) {
        if ( ( id0 == null ) || ( id1 == null ) ) {
            throw new IllegalArgumentException( "attempt to create binary domain combination using null" );
        }
        if ( ( id0.indexOf( SEPARATOR ) != -1 ) || ( id1.indexOf( SEPARATOR ) != -1 ) ) {
            throw new IllegalArgumentException( "ill formatted domain id: " + id0 + ", " + id1 );
        }
        if ( id0.toLowerCase().compareTo( id1.toLowerCase() ) < 0 ) {
            _id0 = BasicDomain.obtainIdAsShort( id0 );
            _id1 = BasicDomain.obtainIdAsShort( id1 );
        }
        else {
            _id0 = BasicDomain.obtainIdAsShort( id1 );
            _id1 = BasicDomain.obtainIdAsShort( id0 );
        }
    }

    @Override
    final public int compareTo( final BinaryDomainCombination binary_domain_combination ) {
        if ( binary_domain_combination.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to compare [" + binary_domain_combination.getClass() + "] to "
                    + "[" + this.getClass() + "]" );
        }
        if ( equals( binary_domain_combination ) ) {
            return 0;
        }
        final int x = getId0().compareTo( binary_domain_combination.getId0() );
        if ( x != 0 ) {
            return x;
        }
        else {
            return getId1().compareTo( binary_domain_combination.getId1() );
        }
    }

    @Override
    final public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to ["
                    + o.getClass() + "]" );
        }
        else {
            return ( getId0Code() == ( ( BinaryDomainCombination ) o ).getId0Code() )
                    && ( getId1Code() == ( ( BinaryDomainCombination ) o ).getId1Code() );
        }
    }

    @Override
    final public String getId0() {
        return BasicDomain.obtainIdFromShort( _id0 );
    }

    @Override
    final public short getId0Code() {
        return _id0;
    }

    @Override
    final public String getId1() {
        return BasicDomain.obtainIdFromShort( _id1 );
    }

    @Override
    final public short getId1Code() {
        return _id1;
    }

    @Override
    final public int hashCode() {
        return calcCode( _id0, _id1 );
    }

    @Override
    final public StringBuffer toGraphDescribingLanguage( final OutputFormat format,
                                                         final String node_attribute,
                                                         final String edge_attribute ) {
        final StringBuffer sb = new StringBuffer();
        switch ( format ) {
            case DOT:
                if ( ForesterUtil.isEmpty( node_attribute ) ) {
                    sb.append( getId0() );
                    sb.append( " -- " );
                    sb.append( getId1() );
                    if ( !ForesterUtil.isEmpty( edge_attribute ) ) {
                        sb.append( " " );
                        sb.append( edge_attribute );
                    }
                    sb.append( ";" );
                }
                else {
                    sb.append( getId0() );
                    sb.append( " " );
                    sb.append( node_attribute );
                    sb.append( ";" );
                    sb.append( ForesterUtil.LINE_SEPARATOR );
                    sb.append( getId1() );
                    sb.append( " " );
                    sb.append( node_attribute );
                    sb.append( ";" );
                    sb.append( ForesterUtil.LINE_SEPARATOR );
                    sb.append( getId0() );
                    sb.append( " -- " );
                    sb.append( getId1() );
                    if ( !ForesterUtil.isEmpty( edge_attribute ) ) {
                        sb.append( " " );
                        sb.append( edge_attribute );
                    }
                    sb.append( ";" );
                }
                break;
            default:
                throw new AssertionError( "unknown format:" + format );
        }
        return sb;
    }

    @Override
    final public String toString() {
        final int code = calcCode( _id0, _id1 );
        if ( S_POOL.containsKey( code ) ) {
            return S_POOL.get( code );
        }
        else {
            final String s = getId0() + SEPARATOR + getId1();
            S_POOL.put( code, s );
            return s;
        }
    }

    public static BinaryDomainCombination obtainInstance( final String ids ) {
        if ( ids.indexOf( BinaryDomainCombination.SEPARATOR ) < 1 ) {
            throw new IllegalArgumentException( "Unexpected format for binary domain combination [" + ids + "]" );
        }
        final String[] ids_ary = ids.split( BinaryDomainCombination.SEPARATOR );
        if ( ids_ary.length != 2 ) {
            throw new IllegalArgumentException( "Unexpected format for binary domain combination [" + ids + "]" );
        }
        return BasicBinaryDomainCombination.obtainInstance( ids_ary[ 0 ], ids_ary[ 1 ] );
    }

    public static BasicBinaryDomainCombination obtainInstance( final String id0, final String id1 ) {
        int code;
        if ( id0.toLowerCase().compareTo( id1.toLowerCase() ) < 0 ) {
            code = calcCode( BasicDomain.obtainIdAsShort( id0 ), BasicDomain.obtainIdAsShort( id1 ) );
        }
        else {
            code = calcCode( BasicDomain.obtainIdAsShort( id1 ), BasicDomain.obtainIdAsShort( id0 ) );
        }
        if ( DC_POOL.containsKey( code ) ) {
            return DC_POOL.get( code );
        }
        else {
            final BasicBinaryDomainCombination dc = new BasicBinaryDomainCombination( id0, id1 );
            DC_POOL.put( code, dc );
            if ( VERBOSE && ( DC_POOL.size() % 100 == 0 ) ) {
                System.out.println( " dc pool size: " + DC_POOL.size() );
            }
            return dc;
        }
    }

    final static int calcCode( final int id0, final int id1 ) {
        return ( id0 * ( Short.MAX_VALUE + 1 ) ) + id1;
    }
}
