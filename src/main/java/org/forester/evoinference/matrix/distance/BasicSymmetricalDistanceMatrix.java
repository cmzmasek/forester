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

package org.forester.evoinference.matrix.distance;

import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.StringTokenizer;

import org.forester.util.ForesterUtil;
import org.forester.util.IllegalFormatUseException;

public final class BasicSymmetricalDistanceMatrix implements DistanceMatrix {

    // NumberFormat                      nf1              = NumberFormat.getInstance();
    private final static NumberFormat PHYLIP_FORMATTER = new DecimalFormat( "0.000000" );
    final String[]                    _identifiers;
    final double[][]                  _values;

    public BasicSymmetricalDistanceMatrix( final int size ) {
        _values = new double[ size ][ size ];
        _identifiers = new String[ size ];
    }

    @Override
    public final String getIdentifier( final int i ) {
        return _identifiers[ i ];
    }

    @Override
    public final int getIndex( final String identifier ) {
        for( int i = 0; i < _identifiers.length; i++ ) {
            if ( getIdentifier( i ).equals( identifier ) ) {
                return i;
            }
        }
        throw new IllegalArgumentException( "identifier [" + identifier + "] not found in distance matrix" );
    }

    @Override
    public final int getSize() {
        return _values.length;
    }

    @Override
    public final double getValue( final int col, final int row ) {
        if ( col == row ) {
            if ( col >= _values.length ) {
                throw new IndexOutOfBoundsException( "" );
            }
            return 0.0;
        }
        else if ( col > row ) {
            return _values[ row ][ col ];
        }
        return _values[ col ][ row ];
    }

    public final double[][] getValues() {
        return _values;
    }

    public final void randomize( final long seed ) {
        final java.util.Random r = new java.util.Random( seed );
        for( int j = 0; j < getSize(); ++j ) {
            for( int i = 0; i < j; ++i ) {
                setValue( i, j, r.nextDouble() );
            }
        }
    }

    @Override
    public final void setIdentifier( final int i, final String identifier ) {
        _identifiers[ i ] = identifier;
    }

    public final void setRow( final String s, final int row ) {
        final StringTokenizer tk = new StringTokenizer( s );
        int i = 0;
        while ( tk.hasMoreElements() ) {
            setValue( i, row, new Double( tk.nextToken() ).doubleValue() );
            i++;
        }
    }

    @Override
    public final void setValue( final int col, final int row, final double d ) {
        if ( d < 0 ) {
            throw new IllegalArgumentException( "negative distance value" );
        }
        if ( ( col == row ) && ( d != 0.0 ) ) {
            throw new IllegalArgumentException( "attempt to set a non-zero value on the diagonal of a symmetrical distance matrix" );
        }
        else if ( col > row ) {
            _values[ row ][ col ] = d;
        }
        _values[ col ][ row ] = d;
    }

    @Override
    public final String toString() {
        return toPhylip().toString();
    }

    @Override
    public final StringBuffer toStringBuffer( final Format format ) {
        switch ( format ) {
            case PHYLIP:
                return toPhylip();
            default:
                throw new IllegalArgumentException( "Unknown format:" + format );
        }
    }

    public final void write( final Writer w ) throws IOException {
        w.write( "    " );
        w.write( getSize() + "" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        for( int row = 0; row < getSize(); ++row ) {
            if ( !ForesterUtil.isEmpty( getIdentifier( row ) ) ) {
                w.write( ForesterUtil.pad( getIdentifier( row ), 10, ' ', false ).toString() );
                w.write( ' ' );
                w.write( ' ' );
            }
            else {
                throw new IllegalFormatUseException( "Phylip format does not allow empty identifiers" );
            }
            for( int col = 0; col < getSize(); ++col ) {
                w.write( PHYLIP_FORMATTER.format( getValue( col, row ) ) );
                if ( col < ( getSize() - 1 ) ) {
                    w.write( ' ' );
                    w.write( ' ' );
                }
            }
            if ( row < ( getSize() - 1 ) ) {
                w.write( ForesterUtil.LINE_SEPARATOR );
            }
        }
    }

    private final StringBuffer toPhylip() {
        final StringBuffer sb = new StringBuffer();
        sb.append( ' ' );
        sb.append( ' ' );
        sb.append( ' ' );
        sb.append( ' ' );
        sb.append( getSize() );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        for( int row = 0; row < getSize(); ++row ) {
            if ( !ForesterUtil.isEmpty( getIdentifier( row ) ) ) {
                sb.append( ForesterUtil.pad( getIdentifier( row ), 10, ' ', false ) );
                sb.append( ' ' );
                sb.append( ' ' );
            }
            else {
                throw new IllegalFormatUseException( "Phylip format does not allow empty identifiers" );
            }
            //sb.append( "" );
            for( int col = 0; col < getSize(); ++col ) {
                sb.append( PHYLIP_FORMATTER.format( getValue( col, row ) ) );
                if ( col < ( getSize() - 1 ) ) {
                    sb.append( ' ' );
                    sb.append( ' ' );
                }
            }
            if ( row < ( getSize() - 1 ) ) {
                sb.append( ForesterUtil.LINE_SEPARATOR );
            }
        }
        return sb;
    }
}
