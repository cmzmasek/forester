// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2014 Christian M. Zmasek
// Copyright (C) 2014 Sanford-Burnham Medical Research Institute
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.msa_compactor;

import org.forester.util.ForesterUtil;

public final class GapContribution implements Comparable<GapContribution> {

    private final String _id;
    private double       _value;

    GapContribution( final String id ) {
        if ( ForesterUtil.isEmpty( id ) ) {
            throw new IllegalArgumentException( "id is empty or null" );
        }
        _id = id;
        _value = 0;
    }

    final String getId() {
        return _id;
    }

    final double getValue() {
        return _value;
    }

    final void addToValue( final double v ) {
        if ( v < 0 ) {
            throw new IllegalArgumentException( "cannot add negative value" );
        }
        _value += v;
    }

    final void divideValue( final double d ) {
        if ( d <= 0 ) {
            throw new IllegalArgumentException( "attempt to divide by non-positive value" );
        }
        _value /= d;
    }

    @Override
    public int compareTo( final GapContribution o ) {
        if ( getValue() < o.getValue() ) {
            return 1;
        }
        else if ( getValue() > o.getValue() ) {
            return -1;
        }
        return 0;
    }
}
