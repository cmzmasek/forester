// $Id:
//
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

package org.forester.protein;

import java.util.HashMap;
import java.util.Map;

import org.forester.util.ForesterUtil;

public class BasicDomain implements Domain {

    private static short                    COUNT        = 0;
    private final static Map<Short, String> ID_TO_STRING = new HashMap<Short, String>();
    private final static Map<String, Short> STRING_TO_ID = new HashMap<String, Short>();
    final private int                       _from;
    final private short                     _id;
    final private short                     _number;
    final private double                    _per_domain_evalue;
    final private double                    _per_domain_score;
    final private int                       _to;
    final private short                     _total_count;

    public BasicDomain( final String id ) {
        if ( ForesterUtil.isEmpty( id ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain with null or empty id" );
        }
        _id = obtainIdAsShort( id );
        _from = -1;
        _to = -1;
        _number = -1;
        _total_count = -1;
        _per_domain_evalue = -1;
        _per_domain_score = -1;
    }

    public BasicDomain( final String id,
                        final int from,
                        final int to,
                        final short number,
                        final short total_count,
                        final double per_domain_evalue,
                        final double per_domain_score ) {
        if ( ( from >= to ) || ( from < 0 ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain from " + from + " to " + to );
        }
        if ( ForesterUtil.isEmpty( id ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain with null or empty id" );
        }
        if ( ( number > total_count ) || ( number < 0 ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain number " + number + " out of "
                    + total_count );
        }
        if ( per_domain_evalue < 0.0 ) {
            throw new IllegalArgumentException( "attempt to create protein domain with negative E-value" );
        }
        _id = obtainIdAsShort( id );
        _from = from;
        _to = to;
        _number = number;
        _total_count = total_count;
        _per_domain_evalue = per_domain_evalue;
        _per_domain_score = per_domain_score;
    }

    /**
     * Basic domains are compared/sorted based upon their identifiers (case
     * insensitive) and their numbers.
     *
     */
    @Override
    public int compareTo( final Domain domain ) {
        if ( domain.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to compare [" + domain.getClass() + "] to " + "["
                    + this.getClass() + "]" );
        }
        if ( this == domain ) {
            return 0;
        }
        return getDomainId().compareTo( domain.getDomainId() );
    }

    /**
     * Basic domains are considered equal if they have the same identifier (case
     * sensitive).
     *
     */
    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return getDomainId().equals( ( ( Domain ) o ).getDomainId() );
        }
    }

    @Override
    public String getDomainId() {
        return obtainIdFromShort( _id );
    }

    @Override
    public int getFrom() {
        return _from;
    }

    @Override
    public int getLength() {
        return ( 1 + getTo() ) - getFrom();
    }

    @Override
    public short getNumber() {
        return _number;
    }

    @Override
    public double getPerDomainEvalue() {
        return _per_domain_evalue;
    }

    @Override
    public double getPerDomainScore() {
        return _per_domain_score;
    }

    @Override
    public int getTo() {
        return _to;
    }

    @Override
    public short getTotalCount() {
        return _total_count;
    }

    @Override
    public int hashCode() {
        return getDomainId().hashCode();
    }

    @Override
    public String toString() {
        return getDomainId();
    }

    public StringBuffer toStringBuffer() {
        return new StringBuffer( getDomainId() );
    }

    public final static short obtainIdAsShort( final String id ) {
        if ( !STRING_TO_ID.containsKey( id ) ) {
            if ( COUNT >= ( Short.MAX_VALUE - 2 ) ) {
                throw new RuntimeException( "too many domain ids!" );
            }
            ID_TO_STRING.put( COUNT, id );
            STRING_TO_ID.put( id, COUNT );
            ++COUNT;
        }
        return STRING_TO_ID.get( id );
    }

    public final static String obtainIdFromShort( final short id ) {
        return ID_TO_STRING.get( id );
    }
}
