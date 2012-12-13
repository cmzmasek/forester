
package org.forester.datastructures;

import java.util.List;

import org.forester.util.ForesterUtil;

public final class IntMatrix {

    private final int    _data[][];
    private final String _labels[];

    public IntMatrix( final int size ) {
        _data = new int[ size ][ size ];
        _labels = new String[ size ];
    }

    public IntMatrix( final List<String> labels ) {
        final int size = labels.size();
        _data = new int[ size ][ size ];
        _labels = new String[ size ];
        for( int i = 0; i < size; ++i ) {
            setLabel( i, labels.get( i ) );
        }
    }

    final public int get( final int x, final int y ) {
        return _data[ x ][ y ];
    }

    final public String getLabel( final int x ) {
        return _labels[ x ];
    }

    final public String getRowAsString( final int x, final char separator ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( getLabel( x ) );
        for( int y = 0; y < size(); ++y ) {
            sb.append( separator );
            sb.append( get( x, y ) );
        }
        return sb.toString();
    }

    final public void inreaseByOne( final int x, final int y ) {
        _data[ x ][ y ]++;
    }

    final public void set( final int x, final int y, final int value ) {
        _data[ x ][ y ] = value;
    }

    final public void setLabel( final int x, final String label ) {
        if ( label == null ) {
            throw new IllegalArgumentException( "matrix label must not be null" );
        }
        _labels[ x ] = label;
    }

    final public int size() {
        return _labels.length;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        for( int x = 0; x < size(); ++x ) {
            if ( getLabel( x ) != null ) {
                sb.append( getLabel( x ) );
                sb.append( "\t" );
            }
            for( int y = 0; y < size(); ++y ) {
                sb.append( get( x, y ) );
                sb.append( "\t" );
            }
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        return sb.toString();
    }
}
