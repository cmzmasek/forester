
package org.forester.evoinference.distance;

import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;

public class S {

    private final List<SortedMap<Double, SortedSet<Integer>>> _data;

    S() {
        _data = new ArrayList<SortedMap<Double, SortedSet<Integer>>>();
    }

    void addPairing( final double key, final int value, final int j ) {
        final SortedMap<Double, SortedSet<Integer>> m = _data.get( j );
        addPairing( key, value, m );
    }

    SortedMap<Double, SortedSet<Integer>> getS( final int j ) {
        return _data.get( j );
    }

    Set<Entry<Double, SortedSet<Integer>>> getSentrySet( final int j ) {
        return _data.get( j ).entrySet();
    }

    void initialize( final BasicSymmetricalDistanceMatrix d ) {
        for( int j = 0; j < d.getSize(); ++j ) {
            final TreeMap<Double, SortedSet<Integer>> map = new TreeMap<Double, SortedSet<Integer>>();
            _data.add( map );
            for( int i = 0; i < j; ++i ) {
                addPairing( d.getValues()[ i ][ j ], i, map );
            }
        }
    }

    void removePairing( final double key, final int value, final int j ) {
        final SortedMap<Double, SortedSet<Integer>> m = _data.get( j );
        final SortedSet<Integer> x = m.get( key );
        if ( x.size() == 1 ) {
            if ( !x.contains( value ) ) {
                //TODO remove me later
                throw new IllegalStateException( "!x.contains( value )" );
            }
            m.remove( key );
        }
        else if ( x.size() > 1 ) {
            final boolean removed = x.remove( value );
            if ( !removed ) {
                //TODO remove me later
                throw new IllegalStateException( value + " not found" );
            }
        }
        else {
            //TODO remove me later
            throw new IllegalStateException( "empty" );
        }
    }

    private static void addPairing( final double key, final int value, final SortedMap<Double, SortedSet<Integer>> m ) {
        if ( !m.containsKey( key ) ) {
            final TreeSet<Integer> x = new TreeSet<Integer>();
            x.add( value );
            m.put( key, x );
        }
        else {
            if ( m.get( key ).contains( value ) ) {
                //TODO remove me later
                throw new IllegalStateException( "pairing " + key + " -> " + value + " already exists" );
            }
            m.get( key ).add( value );
        }
    }
}
