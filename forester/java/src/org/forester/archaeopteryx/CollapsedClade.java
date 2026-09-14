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

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * How a collapsed clade is drawn and named -- the pure rules, a PORT of Archaeopteryx.js ({@code collapsedRows},
 * {@code collapsedName}, {@code collapsedLabel}, {@code collapsedColor}, {@code collapsedFullyMarked} and the wedge
 * style in {@code drawCollapsedClades}, archaeopteryx.js 407f204). Christian, 2026-09-13: the collapsed-subtree
 * display must be the same in both viewers.
 * <ul>
 * <li>A collapsed clade is a leaf that takes {@link #rows} rows: 1 + log2(tips) / 4, at least 1 and at most 2.5. In the
 * breadth layout it weighs {@link #rowWeight} = 2 rows - 1 against a tip's 1, so the gap to a neighbour is the mean of
 * the two weights.</li>
 * <li>It is drawn as a wedge from its node, one edge to the clade's nearest tip and the other to its farthest, {@link
 * #height} tall, filled in the colour most of its tips wear under Color-by.</li>
 * <li>It is named by the node's own name, else the one Color-by value at least 95% of its tips share, else the tips'
 * common name prefix; the label adds " · N tips" and, while a search hits inside, " [found/total]".</li>
 * </ul>
 */
final class CollapsedClade {

    /** A collapsed clade never draws taller than this many rows, however many tips it hides. */
    static final double MAX_ROWS = 2.5;
    /** The wedge's height is this share of the rows it spans, so neighbouring wedges do not touch. */
    static final double HEIGHT_OF_ROWS = 0.82;
    /** The wedge is never thinner than this (px). */
    static final double MIN_HEIGHT = 6;
    /** The share of a clade's tips that must carry one Color-by value for the value to name the clade. */
    static final double NAME_VALUE_SHARE = 0.95;
    /** The wedge fill's opacity, and its opacity when every hidden tip is a hit. */
    static final double FILL_ALPHA = 0.22;
    static final double FILL_ALPHA_MARKED = 0.45;
    /** The wedge outline's opacity, and its width without and with a hit inside. */
    static final double STROKE_ALPHA = 0.9;
    static final float STROKE_WIDTH = 1f;
    static final float STROKE_WIDTH_MARKED = 1.5f;
    /** The separator between a clade's name and its tip count (a middle dot). */
    static final String NAME_SEPARATOR = " · ";
    /** A common name prefix is used only when at least this long once its trailing separators are dropped. */
    static final int MIN_PREFIX_NAME_LENGTH = 2;

    private CollapsedClade() {
    }

    /** The rows a collapsed clade of {@code tips} tips takes: 1 for a pair, growing with the logarithm of the tip
     *  count, capped at {@link #MAX_ROWS}. A single tip counts as a pair. */
    static double rows( final int tips ) {
        final double log2 = Math.log( Math.max( 2, tips ) ) / Math.log( 2 );
        return Math.max( 1, Math.min( MAX_ROWS, 1 + ( log2 / 4 ) ) );
    }

    /** The clade's weight in the breadth layout, against a tip's 1: {@code 2 rows - 1}. */
    static double rowWeight( final int tips ) {
        return ( 2 * rows( tips ) ) - 1;
    }

    /** The wedge's height (px) for a clade of {@code tips} tips when one row is {@code row_unit} px. */
    static double height( final int tips, final double row_unit ) {
        return Math.max( MIN_HEIGHT, rows( tips ) * row_unit * HEIGHT_OF_ROWS );
    }

    /**
     * The Color-by value shared by at least {@link #NAME_VALUE_SHARE} of a clade's {@code tips} tips, or null.
     * {@code values} holds each tip's value, null for a tip without one -- a tip without a value still counts in the
     * denominator.
     */
    static String dominantValue( final List<String> values, final int tips ) {
        if ( ( values == null ) || ( tips < 1 ) ) {
            return null;
        }
        final Map<String, Integer> counts = new HashMap<>();
        String best = null;
        for( final String v : values ) {
            if ( v == null ) {
                continue;
            }
            final int c = counts.merge( v, 1, Integer::sum );
            if ( ( best == null ) || ( c > counts.get( best ) ) ) {
                best = v;
            }
        }
        return ( ( best != null ) && ( counts.get( best ) >= ( NAME_VALUE_SHARE * tips ) ) ) ? best : null;
    }

    /** A common name prefix as a clade name: its trailing separators (whitespace and {@code _ - . : | /}) dropped,
     *  and "" unless at least {@link #MIN_PREFIX_NAME_LENGTH} characters remain. */
    static String prefixName( final String prefix ) {
        if ( prefix == null ) {
            return "";
        }
        int end = prefix.length();
        while ( ( end > 0 ) && isTrailingSeparator( prefix.charAt( end - 1 ) ) ) {
            --end;
        }
        final String p = prefix.substring( 0, end );
        return ( p.length() >= MIN_PREFIX_NAME_LENGTH ) ? p : "";
    }

    private static boolean isTrailingSeparator( final char c ) {
        return Character.isWhitespace( c ) || Character.isSpaceChar( c ) || ( "_-.:|/".indexOf( c ) >= 0 );
    }

    /** The clade's name: its node's own name (trimmed), else the dominant Color-by value, else the prefix name. */
    static String name( final String node_name, final String dominant_value, final String prefix ) {
        final String own = ( node_name == null ) ? "" : PropertyColorScheme.jsTrim( node_name );
        if ( !own.isEmpty() ) {
            return own;
        }
        if ( dominant_value != null ) {
            return dominant_value;
        }
        return prefixName( prefix );
    }

    /** The clade's label: {@code "name · N tips"} ({@code "1 tip"}), without the name when there is none, followed by
     *  {@code " [found/total]"} while any of its tips is a search hit. */
    static String label( final String name, final int total, final int found ) {
        final StringBuilder sb = new StringBuilder();
        if ( ( name != null ) && !name.isEmpty() ) {
            sb.append( name ).append( NAME_SEPARATOR );
        }
        sb.append( total ).append( ( total == 1 ) ? " tip" : " tips" );
        if ( found > 0 ) {
            sb.append( " [" ).append( found ).append( '/' ).append( total ).append( ']' );
        }
        return sb.toString();
    }

    /** The most frequent non-null vote, or null when there is none. On a tie the vote that reached the top count first
     *  wins, as in the JS running count. */
    static <T> T mostFrequent( final List<T> votes ) {
        final Map<T, Integer> counts = new HashMap<>();
        T best = null;
        for( final T v : votes ) {
            if ( v == null ) {
                continue;
            }
            final int c = counts.merge( v, 1, Integer::sum );
            if ( ( best == null ) || ( c > counts.get( best ) ) ) {
                best = v;
            }
        }
        return best;
    }

    /** Whether every tip the clade hides is a hit -- the wedge is then filled, and its label set, in the hit colour. */
    static boolean fullyMarked( final int found, final int total ) {
        return ( found > 0 ) && ( found == total );
    }
}
