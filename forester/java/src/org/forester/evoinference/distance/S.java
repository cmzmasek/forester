
package org.forester.evoinference.distance;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;

public final class S {

    private final static boolean                              DEBUG = true;
    private final List<SortedMap<Double, SortedSet<Integer>>> _data;

    public S() {
        _data = new ArrayList<SortedMap<Double, SortedSet<Integer>>>();
    }

    final public void addPairing( final double key, final int value, final int j ) {
        addPairing( key, value, getS( j ) );
    }

    final public SortedMap<Double, SortedSet<Integer>> getS( final int j ) {
        return _data.get( j );
    }

    final public SortedSet<Integer> getValues( final double key, final int j ) {
        return getS( j ).get( key );
    }

    final public void initialize( final BasicSymmetricalDistanceMatrix d ) {
        for( int j = 0; j < d.getSize(); ++j ) {
            final TreeMap<Double, SortedSet<Integer>> map = new TreeMap<Double, SortedSet<Integer>>();
            _data.add( map );
            for( int i = 0; i < j; ++i ) {
                addPairing( d.getValues()[ i ][ j ], i, map );
            }
        }
        System.out.println( toString() );
    }

    final public void initialize( final int size ) {
        for( int j = 0; j < size; ++j ) {
            final TreeMap<Double, SortedSet<Integer>> map = new TreeMap<Double, SortedSet<Integer>>();
            _data.add( map );
        }
    }

    final public void removePairing( final double key, final int value, final int j ) {
        final SortedMap<Double, SortedSet<Integer>> m = _data.get( j );
        final SortedSet<Integer> x = m.get( key );
        if ( x.size() == 1 ) {
            if ( DEBUG ) {
                if ( !x.contains( value ) ) {
                    throw new IllegalArgumentException( "pairing " + key + " -> " + value + " does not exist for row "
                            + j );
                }
            }
            m.remove( key );
        }
        else if ( x.size() > 1 ) {
            if ( DEBUG ) {
                if ( !x.remove( value ) ) {
                    throw new IllegalArgumentException( "pairing " + key + " -> " + value
                            + " does not exist (could not be removed) for row " + j );
                }
            }
            else {
                x.remove( value );
            }
        }
        else if ( DEBUG ) {
            throw new IllegalStateException();
        }
    }

    final public int size() {
        return _data.size();
    }

    // Slow, only for testing
    @SuppressWarnings( "unchecked")
    final public SortedSet<Integer>[] toArray( final int j ) {
        return _data.get( j ).values().toArray( new SortedSet[ _data.get( j ).size() ] );
    }

    @Override
    final public String toString() {
        final DecimalFormat df = new DecimalFormat( "0.00" );
        final StringBuilder sb = new StringBuilder();
        for( int j = 0; j < size(); ++j ) {
            for( final Entry<Double, SortedSet<Integer>> entry : getSentrySet( j ) ) {
                final double key = entry.getKey();
                final SortedSet<Integer> values = entry.getValue();
                sb.append( df.format( key ) + "->" );
                boolean first = true;
                for( final Integer v : values ) {
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

    final Set<Entry<Double, SortedSet<Integer>>> getSentrySet( final int j ) {
        return getS( j ).entrySet();
    }

    final private static void addPairing( final double key,
                                          final int value,
                                          final SortedMap<Double, SortedSet<Integer>> m ) {
        if ( !m.containsKey( key ) ) {
            final TreeSet<Integer> x = new TreeSet<Integer>();
            x.add( value );
            m.put( key, x );
        }
        else {
            if ( DEBUG ) {
                if ( m.get( key ).contains( value ) ) {
                    throw new IllegalArgumentException( "pairing " + key + " -> " + value + " already exists" );
                }
            }
            m.get( key ).add( value );
        }
    }
}
