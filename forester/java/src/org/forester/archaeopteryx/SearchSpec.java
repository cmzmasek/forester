// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.archaeopteryx;

import java.util.regex.Pattern;

/**
 * One immutable search predicate for the redesigned search tool: a {@link SearchField} (what to look at), a
 * {@link SearchMode} (how to compare), a query {@code value} (and a second {@code value2} for the numeric
 * {@link SearchMode#RANGE} upper bound), plus the two orthogonal modifiers that stay checkboxes -- case
 * sensitivity and a global {@code inverse} that flips the whole predicate.
 *
 * <p>Replaces the fistful of boolean flags threaded through the legacy {@code searchData(...)} calls with a
 * single value object. Interpreted by {@link SearchMatcher}.
 */
final class SearchSpec {

    private final SearchField _field;
    private final SearchMode  _mode;
    private final String      _value;
    private final String      _value2;
    private final boolean     _case_sensitive;
    private final boolean     _inverse;
    // lazily compiled once for a REGEX/WHOLE_WORD spec, then reused across every node of one search
    private Pattern           _compiled_pattern;
    private boolean           _pattern_compiled;

    SearchSpec( final SearchField field, final SearchMode mode, final String value, final String value2,
                final boolean case_sensitive, final boolean inverse ) {
        if ( field == null ) {
            throw new IllegalArgumentException( "search field is null" );
        }
        if ( mode == null ) {
            throw new IllegalArgumentException( "search mode is null" );
        }
        if ( field.isNumeric() != mode.isNumeric() ) {
            throw new IllegalArgumentException( "mode [" + mode + "] does not apply to field [" + field.label()
                    + "] (numeric mismatch)" );
        }
        _field = field;
        _mode = mode;
        _value = value;
        _value2 = value2;
        _case_sensitive = case_sensitive;
        _inverse = inverse;
    }

    /** Convenience for a string predicate with no range and default modifiers off. */
    SearchSpec( final SearchField field, final SearchMode mode, final String value ) {
        this( field, mode, value, null, false, false );
    }

    SearchField field() {
        return _field;
    }

    SearchMode mode() {
        return _mode;
    }

    /** The query text, or -- for a numeric mode -- the operand (the low bound for {@link SearchMode#RANGE}). */
    String value() {
        return _value;
    }

    /** The upper bound for {@link SearchMode#RANGE}; unused (may be {@code null}) for every other mode. */
    String value2() {
        return _value2;
    }

    boolean isCaseSensitive() {
        return _case_sensitive;
    }

    /** Whether the whole predicate is negated (the global "inverse" modifier). */
    boolean isInverse() {
        return _inverse;
    }

    /** The compiled {@link Pattern} for a REGEX/WHOLE_WORD spec (computed once and cached, so a whole-tree search
     *  reuses it across every node); {@code null} for other modes, an empty query, or an invalid regex. */
    Pattern compiledPattern() {
        if ( !_pattern_compiled ) {
            _compiled_pattern = SearchMatcher.compilePattern( this );
            _pattern_compiled = true;
        }
        return _compiled_pattern;
    }
}
