
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
