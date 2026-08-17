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

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/**
 * The matching engine for the redesigned search tool: applies a {@link SearchSpec} (field + mode + value(s) +
 * modifiers) to a node or a whole phylogeny. Pure and Swing-free, so it is fully headless-testable -- the point
 * of the redesign is that the search semantics live here, not scattered as boolean flags through the GUI.
 *
 * <p>The per-node predicate is: extract the field's value(s), test each against the query with the mode (a
 * match against ANY value is a match), then flip the whole result when {@link SearchSpec#isInverse()} is set.
 * With inverse on, a node that has no value for the field matches (there is nothing to match), which mirrors
 * the intuitive "select the nodes that do NOT match"; a driver may further restrict inverse to data-bearing
 * nodes if it wants the legacy behaviour exactly.
 */
final class SearchMatcher {

    private SearchMatcher() {
    }

    /** Whether {@code node} satisfies {@code spec}, including the {@link SearchSpec#isInverse() inverse} flag. */
    static boolean matches( final SearchSpec spec, final PhylogenyNode node ) {
        final boolean base = matchesPositive( spec, node );
        return spec.isInverse() ? !base : base;
    }

    /** The positive predicate (field + mode + value + case sensitivity), ignoring the inverse flag. */
    static boolean matchesPositive( final SearchSpec spec, final PhylogenyNode node ) {
        if ( ( spec == null ) || ( node == null ) ) {
            return false;
        }
        return spec.field().isNumeric() ? matchesNumeric( spec, node ) : matchesString( spec, node );
    }

    /** The ids of the nodes of {@code phy} matching {@code spec} (inverse applied), in preorder. */
    static List<Long> search( final SearchSpec spec, final Phylogeny phy ) {
        final List<Long> ids = new ArrayList<>();
        if ( ( spec == null ) || ( phy == null ) || phy.isEmpty() ) {
            return ids;
        }
        for ( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( matches( spec, node ) ) {
                ids.add( node.getId() );
            }
        }
        return ids;
    }

    // ---- string matching ------------------------------------------------------------------------------------

    private static boolean matchesString( final SearchSpec spec, final PhylogenyNode node ) {
        if ( ForesterUtil.isEmpty( spec.value() ) ) {
            return false;
        }
        for ( final String value : spec.field().stringValues( node ) ) {
            if ( matchOneString( value, spec ) ) {
                return true;
            }
        }
        return false;
    }

    private static boolean matchOneString( final String s, final SearchSpec spec ) {
        if ( ForesterUtil.isEmpty( s ) ) {
            return false;
        }
        final String raw_s = s.trim();
        final SearchMode mode = spec.mode();
        // REGEX / WHOLE_WORD go through a Pattern compiled ONCE per search (cached on the spec), not per value.
        if ( ( mode == SearchMode.REGEX ) || ( mode == SearchMode.WHOLE_WORD ) ) {
            final Pattern p = spec.compiledPattern();
            return ( p != null ) && p.matcher( raw_s ).find();
        }
        final boolean case_sensitive = spec.isCaseSensitive();
        final String v = case_sensitive ? raw_s : raw_s.toLowerCase( Locale.ROOT );
        final String q = case_sensitive ? spec.value().trim() : spec.value().trim().toLowerCase( Locale.ROOT );
        switch ( mode ) {
            case CONTAINS:
                return v.contains( q );
            case STARTS_WITH:
                return v.startsWith( q );
            case ENDS_WITH:
                return v.endsWith( q );
            default:
                return false;
        }
    }

    /** Compiles the {@link Pattern} for a REGEX or WHOLE_WORD spec (case-insensitive unless the spec is
     *  case-sensitive), or {@code null} for any other mode, an empty query, or an invalid regex. Called once per
     *  {@link SearchSpec} (see {@link SearchSpec#compiledPattern()}), so a search over a whole tree compiles once. */
    static Pattern compilePattern( final SearchSpec spec ) {
        final String q = ( spec.value() == null ) ? null : spec.value().trim();
        if ( ForesterUtil.isEmpty( q ) ) {
            return null;
        }
        final int flags = spec.isCaseSensitive() ? 0 : Pattern.CASE_INSENSITIVE;
        try {
            if ( spec.mode() == SearchMode.REGEX ) {
                return Pattern.compile( q, flags );
            }
            if ( spec.mode() == SearchMode.WHOLE_WORD ) {
                return Pattern.compile( "(^|\\s)" + Pattern.quote( q ) + "($|\\s)", flags );
            }
        }
        catch ( final PatternSyntaxException e ) {
            return null;
        }
        return null;
    }

    // ---- numeric matching -----------------------------------------------------------------------------------

    private static boolean matchesNumeric( final SearchSpec spec, final PhylogenyNode node ) {
        final Double a = parseFiniteDouble( spec.value() );
        if ( a == null ) {
            return false;
        }
        double lo = a;
        double hi = a;
        if ( spec.mode() == SearchMode.RANGE ) {
            final Double b = parseFiniteDouble( spec.value2() );
            if ( b == null ) {
                return false;
            }
            lo = Math.min( a, b );
            hi = Math.max( a, b );
        }
        for ( final double v : spec.field().numericValues( node ) ) {
            if ( numericCompare( spec.mode(), v, a, lo, hi ) ) {
                return true;
            }
        }
        return false;
    }

    private static boolean numericCompare( final SearchMode mode, final double v, final double a, final double lo,
                                           final double hi ) {
        switch ( mode ) {
            case EQ:
                return approxEqual( v, a );
            case NE:
                return !approxEqual( v, a );
            case LT:
                return v < a;
            case LE:
                return v <= a;
            case GT:
                return v > a;
            case GE:
                return v >= a;
            case RANGE:
                return ( v >= lo ) && ( v <= hi );
            default:
                return false;
        }
    }

    /** Equality with a small relative tolerance, so a typed "0.5" matches a stored 0.5 despite float noise. */
    private static boolean approxEqual( final double a, final double b ) {
        return Math.abs( a - b ) <= ( 1e-9 * Math.max( 1.0, Math.abs( b ) ) );
    }

    /** Parses a finite {@code double} (trimmed); returns {@code null} for empty, non-numeric, NaN or infinite. */
    static Double parseFiniteDouble( final String s ) {
        if ( ForesterUtil.isEmpty( s ) ) {
            return null;
        }
        try {
            final double d = Double.parseDouble( s.trim() );
            return Double.isFinite( d ) ? Double.valueOf( d ) : null;
        }
        catch ( final NumberFormatException e ) {
            return null;
        }
    }
}
