
package org.forester.evoinference.distance;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;

public final class Sarray {

    public final static int                       FACTOR = 1000000;
    private final static boolean                  DEBUG  = true;
    private final List<SortedMap<Integer, int[]>> _data;

    public Sarray() {
        _data = new ArrayList<SortedMap<Integer, int[]>>();
    }

    final public void addPairing( final double key, final int value, final int j ) {
        addPairing( ( int ) ( FACTOR * key ), value, getS( j ) );
    }

    final public void addPairing( final int key, final int value, final int j ) {
        addPairing( key, value, getS( j ) );
    }

    final public SortedMap<Integer, int[]> getS( final int j ) {
        return _data.get( j );
    }

    final public int[] getValues( final int key, final int j ) {
        return getS( j ).get( key );
    }

    final public void initialize( final BasicSymmetricalDistanceMatrix d ) {
        for( int j = 0; j < d.getSize(); ++j ) {
            final TreeMap<Integer, int[]> map = new TreeMap<Integer, int[]>();
            _data.add( map );
            for( int i = 0; i < j; ++i ) {
                addPairing( ( int ) ( FACTOR * d.getValues()[ i ][ j ] ), i, map );
            }
        }
        //System.out.println( toString() );
    }

    final public void initialize( final int size ) {
        for( int j = 0; j < size; ++j ) {
            final TreeMap<Integer, int[]> map = new TreeMap<Integer, int[]>();
            _data.add( map );
        }
    }

    final public void removePairing( final double key, final int value, final int j ) {
        removePairing( ( int ) ( key * FACTOR ), value, j );
    }

    final public void removePairing( final int key, final int value, final int j ) {
        final SortedMap<Integer, int[]> m = _data.get( j );
        final int[] x = m.get( key );
        if ( x == null ) {
            System.out.println();
            System.out
                    .println( "________________________________________________________________________________________" );
            System.out.println( toString() );
            throw new IllegalArgumentException( "key " + key + " (->" + value + ") does not exist for row " + j );
        }
        if ( x.length == 1 ) {
            m.remove( key );
        }
        else {
            int[] xnew = new int[ x.length - 1 ];
            int xc = 0;
            for( int i = 0; ++i < x.length; ++i ) {
                int xv = x[ i ];
                if ( xv != value ) {
                    xnew[ xc++ ] = xv;
                }
            }
            m.put( key, xnew );
        }
    }

    final public int size() {
        return _data.size();
    }

    // Slow, only for testing
    @SuppressWarnings( "unchecked")
    final public Set<Integer>[] toArray( final int j ) {
        return _data.get( j ).values().toArray( new Set[ _data.get( j ).size() ] );
    }

    @Override
    final public String toString() {
        final DecimalFormat df = new DecimalFormat( "0.000000" );
        final StringBuilder sb = new StringBuilder();
        for( int j = 0; j < size(); ++j ) {
            sb.append( j );
            sb.append( ": " );
            for( final Entry<Integer, int[]> entry : getSentrySet( j ) ) {
                final double key = entry.getKey();
                final int[] values = entry.getValue();
                sb.append( df.format( key / FACTOR ) + "->" );
                boolean first = true;
                for( final int v : values ) {
                    if ( !first ) {
                        sb.append( "," );
                    }
                    first = false;
                    sb.append( v );
                }
                sb.append( "  " );
            }
            sb.append( "\n" );
        }
        return sb.toString();
    }

    final Set<Entry<Integer, int[]>> getSentrySet( final int j ) {
        return getS( j ).entrySet();
    }

    final private static void addPairing( final int key, final int value, final SortedMap<Integer, int[]> m ) {
        if ( !m.containsKey( key ) ) {
            final int[] x = new int[ 1 ];
            x[ 0 ] = value;
            m.put( key, x );
        }
        else {
            final int[] x = new int[ m.get( key ).length + 1 ];
            for( int i = 0; i < x.length - 1; i++ ) {
                x[ i ] = m.get( key )[ i ];
            }
            x[ x.length - 1 ] = value;
            m.put( key, x );
        }
    }
}
