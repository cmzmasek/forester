
package org.forester.surfacing;

import java.util.SortedSet;

final class UniprotNameResult {

    final String            _name;
    final Long              _count;
    final int               _out_of;
    final SortedSet<String> _all_names;
    public UniprotNameResult( final String name,
                              final Long count,
                              final int out_of,
                              final SortedSet<String> all_names ) {
        _name = name;
        _count = count;
        _out_of = out_of;
        _all_names = all_names;
    }

    public String getName() {
        return _name;
    }

    public Long getCount() {
        return _count;
    }

    public int getOutOf() {
        return _out_of;
    }

    public SortedSet<String> getAllNames() {
        return _all_names;
    }
}
