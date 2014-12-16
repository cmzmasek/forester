
package org.forester.util;

import java.util.Comparator;

public class StringInt {

    private final String _s;
    private final int    _i;

    public StringInt( final String s, final int i ) {
        _s = s;
        _i = i;
    }

    public String getString() {
        return _s;
    }

    public int getInt() {
        return _i;
    }

    public static final class DescendingIntComparator implements Comparator<StringInt> {

        @Override
        public int compare( final StringInt o1, final StringInt o2 ) {
            if ( o1.getInt() > o2.getInt() ) {
                return -1;
            }
            if ( o1.getInt() < o2.getInt() ) {
                return 1;
            }
            return 0;
        }
    }
}
