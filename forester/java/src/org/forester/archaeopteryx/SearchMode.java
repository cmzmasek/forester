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

/**
 * A match mode for the redesigned search tool: how a {@link SearchField}'s value(s) are compared against
 * the query. String fields use the {@link #stringModes()}; numeric fields (branch length, support, numeric
 * custom properties) use the {@link #numericModes()}. The two sets are disjoint and a mode is only valid for
 * a field of its own {@link #isNumeric() kind}.
 *
 * <p>Negation is intentionally NOT expressed as a mode (there is no "does not contain"); the numeric set does
 * carry a first-class {@link #NE} (&ne;) because it reads naturally as an operator, and a single global
 * "inverse" flag on {@link SearchSpec} flips the whole predicate. See the search-redesign notes.
 */
enum SearchMode {

    // string modes
    CONTAINS( "contains", false ),
    STARTS_WITH( "starts with", false ),
    ENDS_WITH( "ends with", false ),
    WHOLE_WORD( "whole word", false ),
    REGEX( "regex", false ),
    // numeric modes -- spelled out (with the symbol as a hint) so they read clearly at small sizes / low
    // resolution and match the word-based string modes above
    EQ( "equals (=)", true ),
    NE( "not equal (!=)", true ),
    LT( "less than (<)", true ),
    LE( "at most (<=)", true ),
    GT( "greater than (>)", true ),
    GE( "at least (>=)", true ),
    RANGE( "range", true );

    private final String  _label;
    private final boolean _numeric;

    SearchMode( final String label, final boolean numeric ) {
        _label = label;
        _numeric = numeric;
    }

    /** A short, user-facing label for the mode selector (e.g. "contains", "&ge;"). */
    String label() {
        return _label;
    }

    /** Whether this mode compares numbers (true) or strings (false). */
    boolean isNumeric() {
        return _numeric;
    }

    /** The modes applicable to a string field, in menu order. */
    static SearchMode[] stringModes() {
        return new SearchMode[] { CONTAINS, STARTS_WITH, ENDS_WITH, WHOLE_WORD, REGEX };
    }

    /** The modes applicable to a numeric field, in menu order. */
    static SearchMode[] numericModes() {
        return new SearchMode[] { EQ, NE, LT, LE, GT, GE, RANGE };
    }

    /** The default mode for a field of the given kind: {@link #CONTAINS} for string, {@link #EQ} for numeric. */
    static SearchMode defaultFor( final boolean numeric ) {
        return numeric ? EQ : CONTAINS;
    }
}
