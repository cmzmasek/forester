
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

    void initialize( final BasicSymmetricalDistanceMatrix d ) {
        for( int j = 0; j < d.getSize(); ++j ) {
            final TreeMap<Double, SortedSet<Integer>> map = new TreeMap<Double, SortedSet<Integer>>();
            _data.add( map );
            for( int i = 0; i < j; ++i ) {
                addValue( d.getValues()[ i ][ j ], i, map );
            }
        }
    }

    SortedMap<Double, SortedSet<Integer>> getS( final int j ) {
        return _data.get( j );
    }

    Set<Entry<Double, SortedSet<Integer>>> getSentrySet( final int j ) {
        return _data.get( j ).entrySet();
    }

    void addValue( double key, int value, int j ) {
        SortedMap<Double, SortedSet<Integer>> m = _data.get( j );
        addValue( key, value, m );
    }

    static void addValue( double key, int value, SortedMap<Double, SortedSet<Integer>> m ) {
        if ( !m.containsKey( key ) ) {
            TreeSet<Integer> x = new TreeSet<Integer>();
            x.add( value );
            m.put( key, x );
        }
        else {
            m.get( key ).add( value );
        }
    }
}
