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

package org.forester.surfacing;

import org.forester.protein.BasicDomain;
import org.forester.protein.Domain;
import org.forester.util.ForesterUtil;

/*
 * A limited implementation of Domain. Its intended use is for when only a
 * domain identifier is needed. Note intended for general use.
 */
public class SimpleDomain implements Domain {

    final private short _id;

    public SimpleDomain( final String id ) {
        if ( ForesterUtil.isEmpty( id ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain with null or empty id" );
        }
        _id = BasicDomain.obtainIdAsShort( id );
    }

    @Override
    public int compareTo( final Domain domain ) {
        if ( this == domain ) {
            return 0;
        }
        return getDomainId().compareTo( domain.getDomainId() );
    }

    @Override
    public String getDomainId() {
        return BasicDomain.obtainIdFromShort( _id );
    }

    @Override
    public int getFrom() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public int getLength() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public short getNumber() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public double getPerDomainEvalue() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public double getPerDomainScore() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public int getTo() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public short getTotalCount() {
        throw new RuntimeException( "method not implemented" );
    }
}
