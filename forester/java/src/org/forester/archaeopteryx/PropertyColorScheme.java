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

import java.awt.Color;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.Property;
import org.forester.util.ForesterUtil;

/**
 * Assigns a color to each distinct value of a chosen phyloXML property (e.g.
 * {@code repseq:host}), so leaves can be colored on the fly by that property: the same
 * value always maps to the same color. The categorical palette is cycled when a property
 * has more distinct values than there are palette entries.
 * <p>
 * Trivially different spellings are grouped before coloring -- values are compared after
 * trimming, collapsing whitespace, treating {@code _} as a space, and case-folding. The one
 * deliberate synonym fold is {@code human}/{@code humans} into {@code Homo sapiens}, so e.g.
 * {@code human}/{@code Human}/{@code homo_sapiens}/{@code Homo sapiens} get one color; other
 * semantically-equal but lexically-different values (such as {@code man} vs {@code H. sapiens})
 * are not merged. A predominantly-numeric column (e.g. {@code year}, {@code age}) is colored by
 * a continuous gradient over its numeric range, while {@code country} and {@code host} first drop a
 * trailing qualifier -- everything from the first {@code :} (country, so {@code USA:CA} and
 * {@code USA:IL} share a color) or {@code ;} (host, so {@code Homo sapiens; male 35} and
 * {@code Homo sapiens; female old} share a color) onward -- before that grouping.
 * <p>
 * This is independent of the explicit {@code style:*} visualization properties and
 * the "Visual Styles" feature.
 */
final class PropertyColorScheme {

    // The default categorical palette IS the shared Tableau 10 (AptxUtil.TABLEAU_10), so color strips, symbol
    // columns, and color-by-property match the rank/clade/series coloring. Extended past ten (rather than plainly
    // repeated) by AptxUtil.qualitativeColor when a property has more distinct values than palette entries.
    private static final Color[] DEFAULT_PALETTE = AptxUtil.TABLEAU_10;
    // The Okabe-Ito colorblind-safe set (grey instead of black so it reads on a dark canvas too).
    private static final Color[] COLORBLIND_PALETTE = { new Color( 0xE69F00 ), new Color( 0x56B4E9 ),
            new Color( 0x009E73 ), new Color( 0xF0E442 ), new Color( 0x0072B2 ), new Color( 0xD55E00 ),
            new Color( 0xCC79A7 ), new Color( 0x999999 ) };

    /** The default palette name; selectable palettes are listed by {@link #paletteNames()}. */
    static final String DEFAULT_PALETTE_NAME = "Default";
    private static final java.util.LinkedHashMap<String, Color[]> PALETTES = new java.util.LinkedHashMap<>();
    static {
        PALETTES.put( DEFAULT_PALETTE_NAME, DEFAULT_PALETTE );
        PALETTES.put( "Colorblind-friendly", COLORBLIND_PALETTE );
    }

    /** The names of the selectable categorical palettes, in display order. */
    static List<String> paletteNames() {
        return new ArrayList<String>( PALETTES.keySet() );
    }

    private static Color[] paletteByName( final String name ) {
        final Color[] p = PALETTES.get( name );
        return ( p != null ) ? p : DEFAULT_PALETTE;
    }

    private final String             _ref;
    // Categorical mode: one palette color per distinct value. _value_to_color maps the
    // (display) representative of each group to its color, for the legend; _key_to_color
    // maps the normalized grouping key to the same color, for looking up a node's color.
    private final Map<String, Color> _value_to_color;
    private final Map<String, Color> _key_to_color;
    // representative label -> number of (visible) leaves in that group, for the legend
    private final Map<String, Integer> _value_to_count;
    // representative label -> its (stable, normalized) group key, for keying per-value overrides
    private final Map<String, String>  _value_to_key;
    private final Color[]              _palette; // the categorical palette in use
    // Continuous mode (numeric properties such as "year"): a blue->red gradient spanning
    // [_min, _max] instead of distinct colors. _gradient is false for categorical refs.
    // NOT final: a heat-map MATRIX column forces continuous mode on the range SHARED across the matrix
    // (see setSharedRange), even when this column's own values would auto-detect as categorical.
    private boolean                  _gradient;
    // NOT final: a heat-map MATRIX column overrides its per-ref range with the range SHARED across the matrix
    // (see setSharedRange), so every column of the matrix reads on the same color scale.
    private double                   _min;
    private double                   _max;
    // For refs whose value carries a trailing qualifier, the delimiter at which the value is
    // truncated before grouping: ':' for "country" (USA:CA == USA:IL), ';' for "host"
    // (Homo sapiens; male 35 == Homo sapiens; female old). 0 means keep the whole value.
    private final char               _truncate_at;
    // Coverage, for the legend's "no value" row (JS-parity: total - coverage, pinned last, dashed swatch):
    // how many visible tips the scheme was built over, and how many of those actually carry a value.
    private final int                _visible_tip_count;
    private final int                _tips_with_nonempty_value;
    private final int                _tips_with_numeric_value;

    PropertyColorScheme( final Phylogeny phylogeny, final String ref ) {
        this( phylogeny, ref, null, DEFAULT_PALETTE_NAME );
    }

    PropertyColorScheme( final Phylogeny phylogeny, final String ref, final Map<String, Color> overrides ) {
        this( phylogeny, ref, overrides, DEFAULT_PALETTE_NAME );
    }

    PropertyColorScheme( final Phylogeny phylogeny, final String ref, final Map<String, Color> overrides,
                         final String palette_name ) {
        this( phylogeny, ref, overrides, palette_name, null, null );
    }

    /**
     * @param overrides    optional user-assigned colors, keyed by group key (see {@link #getValueKeys()}),
     *                     that replace the automatic palette color for those values; may be null/empty.
     * @param palette_name the categorical palette to assign colors from (see {@link #paletteNames()}).
     * @param memory       optional value-color IDENTITY memory (group key -> color), MUTATED: a value keeps its
     *                     remembered color across view changes (subtree navigation, collapse, deletion), and a
     *                     value met for the first time takes the next free palette slot and is remembered.
     *                     JS-parity: in Archaeopteryx.js a value's color is an identity for the whole session, so
     *                     diving into a subtree never recolors what stays visible -- without this, the
     *                     frequency-sorted palette assignment RE-SPREADS on every rebuild and a subtree figure's
     *                     colors stop matching the whole-tree figure of the same data. On the FIRST build the
     *                     memory is empty, so assignment equals the plain frequency-indexed palette (unchanged
     *                     first-view behavior). Gradient (numeric) mode deliberately ignores it: a ramp's color
     *                     is position in the VIEW's range. Overrides win over the memory but are never stored in
     *                     it, so clearing an override returns the remembered automatic color. Null = no memory
     *                     (legacy per-view re-spread; annotation-column schemes still use this).
     * @param memory_next  the memory's next-free-palette-slot counter ({@code int[1]}, MUTATED); required when
     *                     {@code memory} is non-null.
     */
    PropertyColorScheme( final Phylogeny phylogeny, final String ref, final Map<String, Color> overrides,
                         final String palette_name, final Map<String, Color> memory, final int[] memory_next ) {
        _ref = ref;
        _palette = paletteByName( palette_name );
        _truncate_at = truncationDelimiter( ref );
        _value_to_color = new LinkedHashMap<String, Color>();
        _key_to_color = new LinkedHashMap<String, Color>();
        _value_to_count = new LinkedHashMap<String, Integer>();
        _value_to_key = new LinkedHashMap<String, String>();
        // Color from the leaves actually on screen (those hidden under a collapsed node are
        // excluded), so the colors and legend track the displayed (sub)tree as the user
        // navigates into subtrees, collapses clades, or deletes nodes.
        final List<PhylogenyNode> leaves = visibleExternalNodes( phylogeny );
        // Coverage counts for the "no value" legend row. BOTH the non-empty and the numeric count are kept,
        // because what "has a value" means depends on the mode -- and the mode can still change after
        // construction (setSharedRange forces a matrix column continuous): categorically, any non-empty value
        // draws a mark; under a gradient, only a PARSEABLE number does (an "n/a" tip draws nothing, so it
        // counts as missing -- the row must agree with what the tree actually shows).
        int nonempty = 0;
        int numeric = 0;
        for( final PhylogenyNode node : leaves ) {
            final String v = valueFor( node, ref );
            if ( !ForesterUtil.isEmpty( v ) && !displayLabel( v ).isEmpty() ) { // fold-to-empty draws no mark
                nonempty++;
                if ( parseNumber( v ) != null ) {
                    numeric++;
                }
            }
        }
        _visible_tip_count = leaves.size();
        _tips_with_nonempty_value = nonempty;
        _tips_with_numeric_value = numeric;
        // Use a continuous gradient when the column is predominantly numeric (year / age / percent-identity /
        // ...); otherwise color by distinct categories. Decided from the visible values, not the ref name.
        _gradient = shouldUseGradient( leaves, ref );
        if ( _gradient ) {
            double min = Double.POSITIVE_INFINITY;
            double max = Double.NEGATIVE_INFINITY;
            for( final PhylogenyNode node : leaves ) {
                final Double d = parseNumber( valueFor( node, ref ) );
                if ( d != null ) {
                    min = Math.min( min, d );
                    max = Math.max( max, d );
                }
            }
            _min = min;
            _max = max;
        }
        else {
            _min = 0;
            _max = 0;
            // Group trivial variants together (case, whitespace, underscores; and, for
            // "country"/"host", the qualifier after ':'/';') so e.g. "Human"/"human"/
            // "homo_sapiens " share one color. Each group's legend label is its most frequent spelling.
            final Map<String, Map<String, Integer>> key_to_label_counts = new HashMap<String, Map<String, Integer>>();
            for( final PhylogenyNode node : leaves ) {
                final String v = valueFor( node, ref );
                if ( !ForesterUtil.isEmpty( v ) ) {
                    final String label = displayLabel( v );
                    if ( label.isEmpty() ) {
                        continue; // a value that folds to nothing (e.g. "_") forms no group (JS rule)
                    }
                    final String key = label.toLowerCase( Locale.ROOT );
                    Map<String, Integer> counts = key_to_label_counts.get( key );
                    if ( counts == null ) {
                        counts = new HashMap<String, Integer>();
                        key_to_label_counts.put( key, counts );
                    }
                    final Integer c = counts.get( label );
                    counts.put( label, ( c == null ) ? 1 : ( c + 1 ) );
                }
            }
            // [ representative label, key ] per group, plus each group's total leaf count
            final List<String[]> groups = new ArrayList<String[]>();
            final Map<String, Integer> key_to_total = new HashMap<String, Integer>();
            for( final Map.Entry<String, Map<String, Integer>> e : key_to_label_counts.entrySet() ) {
                int total = 0;
                for( final int c : e.getValue().values() ) {
                    total += c;
                }
                key_to_total.put( e.getKey(), total );
                groups.add( new String[] { representative( e.getValue() ), e.getKey() } );
            }
            // Most frequent first (ties broken alphabetically), so the most common values get the
            // most distinct palette colors and head the legend; palette cycling, when there are more
            // distinct values than colors, then only affects the rarest values.
            Collections.sort( groups, new Comparator<String[]>() {

                @Override
                public int compare( final String[] a, final String[] b ) {
                    final int by_count = Integer.compare( key_to_total.get( b[ 1 ] ), key_to_total.get( a[ 1 ] ) );
                    return ( by_count != 0 ) ? by_count : String.CASE_INSENSITIVE_ORDER.compare( a[ 0 ], b[ 0 ] );
                }
            } );
            int i = 0;
            for( final String[] g : groups ) {
                Color color;
                if ( memory != null ) { // identity memory: keep a known value's color, remember a new one
                    color = memory.get( g[ 1 ] );
                    if ( color == null ) {
                        color = AptxUtil.qualitativeColor( _palette, memory_next[ 0 ]++ );
                        memory.put( g[ 1 ], color );
                    }
                }
                else {
                    color = AptxUtil.qualitativeColor( _palette, i++ );
                }
                if ( ( overrides != null ) && overrides.containsKey( g[ 1 ] ) ) {
                    color = overrides.get( g[ 1 ] ); // user-assigned color for this value (never stored in memory)
                }
                _value_to_color.put( g[ 0 ], color ); // _value_to_color is now ordered most-frequent first
                _key_to_color.put( g[ 1 ], color );
                _value_to_count.put( g[ 0 ], key_to_total.get( g[ 1 ] ) );
                _value_to_key.put( g[ 0 ], g[ 1 ] );
            }
        }
    }

    String getRef() {
        return _ref;
    }

    boolean isEmpty() {
        return _gradient ? ( _min > _max ) : _value_to_color.isEmpty();
    }

    /** Whether this scheme colors by a continuous range (a gradient) rather than distinct values. */
    boolean isGradient() {
        return _gradient;
    }

    int numberOfValues() {
        return _value_to_color.size();
    }

    /** Number of visible external nodes the scheme was built over (tips hidden under a collapse excluded). */
    int visibleTipCount() {
        return _visible_tip_count;
    }

    /**
     * How many of the visible tips draw NO mark under this scheme -- the count behind the legend's
     * "no value" row. Defined as "would {@link #colorFor} return null": categorically, tips whose value is
     * empty/absent; under a gradient, also tips whose value does not parse as a number (they draw nothing,
     * so a row claiming they were covered would lie).
     */
    int missingCount() {
        return _visible_tip_count - ( _gradient ? _tips_with_numeric_value : _tips_with_nonempty_value );
    }

    /** The color for this node's value of the property, or {@code null} if it has none. */
    Color colorFor( final PhylogenyNode node ) {
        if ( _gradient ) {
            final Double d = parseNumber( valueFor( node, _ref ) );
            if ( d == null ) {
                return null;
            }
            final double t = ( _max > _min ) ? ( ( d - _min ) / ( _max - _min ) ) : 0.0;
            return gradientColorAt( t );
        }
        final String v = valueFor( node, _ref );
        return ForesterUtil.isEmpty( v ) ? null : _key_to_color.get( groupKey( v ) );
    }

    /**
     * For a gradient (numeric) scheme, the node's value as a fraction in {@code [0, 1]} from the minimum to
     * the maximum visible value (used to scale a bar), or {@code null} when this is not a gradient scheme or
     * the node has no numeric value.
     */
    Double gradientFraction( final PhylogenyNode node ) {
        if ( !_gradient ) {
            return null;
        }
        final Double d = parseNumber( valueFor( node, _ref ) );
        if ( d == null ) {
            return null;
        }
        final double t = ( _max > _min ) ? ( ( d - _min ) / ( _max - _min ) ) : 0.0;
        return ( t < 0.0 ) ? 0.0 : ( ( t > 1.0 ) ? 1.0 : t ); // clamp so a bar never over/underflows its column
    }

    /** Representative-label to color map of all distinct values, ordered most-frequent first. */
    Map<String, Color> getValueColors() {
        return _value_to_color;
    }

    /** Representative-label to (visible) leaf count for each value; empty in gradient mode. */
    Map<String, Integer> getValueCounts() {
        return _value_to_count;
    }

    /** Representative-label to its stable group key (the key to use for a per-value color override). */
    Map<String, String> getValueKeys() {
        return _value_to_key;
    }

    /**
     * The display label for a raw property value: trimmed, underscores as spaces, internal
     * whitespace collapsed; for refs that carry a trailing qualifier ("country", "host") only
     * the part before the first ':'/';' (so "USA:CA" reads as "USA" and "Homo sapiens; male 35"
     * reads as "Homo sapiens"). The synonym "human" is folded to "Homo sapiens". Case is otherwise
     * preserved -- this is what the legend shows.
     */
    private String displayLabel( final String v ) {
        String s = v;
        if ( _truncate_at != 0 ) {
            final int idx = s.indexOf( _truncate_at );
            if ( idx >= 0 ) {
                s = s.substring( 0, idx );
                // the cut may have landed INSIDE a parenthetical ("Saimiri boliviensis (squirrel monkey;
                // voucher: X)" cut at ';' leaves an unclosed '(' ): trim back to before the first unmatched '('
                int first_unmatched = -1;
                int depth = 0;
                for( int i = 0; i < s.length(); ++i ) {
                    final char c = s.charAt( i );
                    if ( c == '(' ) {
                        if ( depth == 0 ) {
                            first_unmatched = i;
                        }
                        depth++;
                    }
                    else if ( ( c == ')' ) && ( depth > 0 ) ) {
                        depth--;
                        if ( depth == 0 ) {
                            first_unmatched = -1;
                        }
                    }
                }
                if ( ( depth > 0 ) && ( first_unmatched >= 0 ) ) {
                    s = s.substring( 0, first_unmatched );
                }
            }
        }
        s = s.trim().replace( '_', ' ' );
        // fold runs of whitespace, then trim AGAIN: replacing underscores can create leading/trailing
        // space ("_cat_" -> " cat "), and "_" alone must fold to EMPTY so it is dropped. (The JS spec says
        // exactly this -- "a value that becomes empty is dropped" -- though its code misses the final trim;
        // flagged to the JS side rather than copying the wart.)
        s = s.replaceAll( "\\s+", " " ).trim();
        return canonicalSynonym( s );
    }

    /**
     * Folds a short list of unambiguous common-animal synonyms to a canonical common name (so e.g.
     * {@code swine}/{@code porcine}/{@code Sus scrofa} share one color and one "Pig" legend row). Matching is
     * WHOLE-VALUE only, never substring -- "ferret badger" and "42-day-old pig" keep their own groups. A miss
     * whose value ends in a parenthetical is retried once with that one trailing {@code (...)} removed
     * ("Bos taurus (cattle)" -> "Cow"); the display form of a miss keeps the parenthetical.
     * <p>
     * CROSS-IMPLEMENTATION CONTRACT: this is Archaeopteryx.js's {@code VIS_SYNONYMS} (forester.js), verbatim,
     * so the two viewers group and label a shared tree identically -- extend BOTH or NEITHER. This is
     * deliberately display grouping, not data cleaning: spelling plus a short unambiguous dictionary, nothing
     * semantic beyond it; raw values are untouched everywhere else (search, exports, the node dialog).
     */
    private static String canonicalSynonym( final String label ) {
        final String lower = label.toLowerCase( Locale.ROOT );
        String hit = SYNONYM_LOOKUP.get( lower );
        if ( hit == null ) {
            final String stripped = lower.replaceAll( "\\s*\\([^()]*\\)\\s*$", "" );
            if ( !stripped.equals( lower ) && !stripped.isEmpty() ) {
                hit = SYNONYM_LOOKUP.get( stripped );
            }
        }
        return ( hit != null ) ? hit : label;
    }

    // lowercase synonym (and lowercase canonical) -> canonical display name; see canonicalSynonym
    private static final Map<String, String> SYNONYM_LOOKUP = buildSynonymLookup();

    private static Map<String, String> buildSynonymLookup() {
        final String[][] table = {
            { "Human", "humans", "homo sapiens", "h. sapiens" },
            { "Cow", "bovine", "calf", "cattle", "bull", "heifer", "bos taurus", "b. taurus" },
            { "Chicken", "broiler chicken", "broiler", "hen", "rooster", "gallus gallus", "g. gallus",
              "gallus gallus domesticus" },
            { "Mouse", "house mouse", "murine", "mus musculus", "m. musculus" },
            { "Rat", "brown rat", "norway rat", "black rat", "rattus norvegicus", "r. norvegicus",
              "rattus rattus" },
            { "Ferret", "domestic ferret", "mustela putorius furo", "mustela furo", "m. putorius furo" },
            { "Guinea pig", "cavy", "domestic guinea pig", "cavia porcellus", "c. porcellus" },
            { "Rhesus monkey", "rhesus macaque", "macaca mulatta", "m. mulatta" },
            { "Rabbit", "european rabbit", "oryctolagus cuniculus", "o. cuniculus" },
            { "Dog", "canine", "canis familiaris", "canis lupus familiaris", "c. familiaris" },
            { "Cat", "feline", "domestic cat", "felis catus", "f. catus", "felis silvestris catus" },
            { "Duck", "mallard", "mallard duck", "domestic duck", "anas platyrhynchos", "a. platyrhynchos" },
            { "Pig", "swine", "porcine", "hog", "piglet", "sus scrofa", "s. scrofa", "sus scrofa domesticus" },
            { "Horse", "equine", "mare", "stallion", "equus caballus", "e. caballus" },
            { "Sheep", "ovine", "lamb", "ewe", "ovis aries", "o. aries" },
            { "Goat", "caprine", "capra hircus", "c. hircus" },
            { "Camel", "dromedary", "bactrian camel", "camelus dromedarius", "camelus bactrianus",
              "c. dromedarius" } };
        final Map<String, String> lookup = new LinkedHashMap<String, String>();
        for( final String[] row : table ) {
            lookup.put( row[ 0 ].toLowerCase( Locale.ROOT ), row[ 0 ] ); // a canonical matches itself
            for( int i = 1; i < row.length; ++i ) {
                lookup.put( row[ i ], row[ 0 ] );
            }
        }
        return lookup;
    }

    /** The normalized key a value is grouped/colored by: its display label, case-folded. */
    private String groupKey( final String v ) {
        return displayLabel( v ).toLowerCase( Locale.ROOT );
    }

    /** The most frequent spelling in a group (ties broken by code-point order, ascending), with its first
     *  character uppercased -- so a legend of raw lowercase values still reads cleanly ("cat" -> "Cat").
     *  A dictionary hit's group holds only the canonical name, which already starts uppercase (JS parity). */
    private static String representative( final Map<String, Integer> label_counts ) {
        String best = null;
        int best_count = -1;
        for( final Map.Entry<String, Integer> e : label_counts.entrySet() ) {
            final int n = e.getValue();
            if ( ( n > best_count ) || ( ( n == best_count ) && ( e.getKey().compareTo( best ) < 0 ) ) ) {
                best = e.getKey();
                best_count = n;
            }
        }
        if ( ( best != null ) && !best.isEmpty() && Character.isLowerCase( best.charAt( 0 ) ) ) {
            best = Character.toUpperCase( best.charAt( 0 ) ) + best.substring( 1 );
        }
        return best;
    }

    /** Color at fraction {@code t} (0..1, low value to high value) of the gradient. */
    // Viridis: a perceptually-uniform sequential ramp (dark blue-purple -> teal -> green -> yellow) -- calmer
    // and colorblind-safe versus a blue->red hue sweep, and it reads on both the white and the dark canvas.
    // Ten anchor stops, linearly interpolated in RGB.
    private static final int[] VIRIDIS = { 0x440154, 0x482878, 0x3E4A89, 0x31688E, 0x26828E, 0x1F9E89, 0x35B779,
            0x6DCD59, 0xB4DE2C, 0xFDE725 };

    Color gradientColorAt( final double t ) {
        final double tt = ( t < 0.0 ) ? 0.0 : ( ( t > 1.0 ) ? 1.0 : t );
        final double x = tt * ( VIRIDIS.length - 1 );
        final int i = (int) Math.floor( x );
        if ( i >= ( VIRIDIS.length - 1 ) ) {
            return new Color( VIRIDIS[ VIRIDIS.length - 1 ] );
        }
        final double f = x - i;
        final Color a = new Color( VIRIDIS[ i ] );
        final Color b = new Color( VIRIDIS[ i + 1 ] );
        return new Color( (int) Math.round( a.getRed() + ( f * ( b.getRed() - a.getRed() ) ) ),
                          (int) Math.round( a.getGreen() + ( f * ( b.getGreen() - a.getGreen() ) ) ),
                          (int) Math.round( a.getBlue() + ( f * ( b.getBlue() - a.getBlue() ) ) ) );
    }

    String getGradientMinLabel() {
        return formatNumber( _min );
    }

    String getGradientMaxLabel() {
        return formatNumber( _max );
    }

    /** Make this scheme a continuous gradient over the range SHARED across a heat-map matrix, so all of the
     *  matrix's columns read on one color scale (and one gradient legend). A MATRIX column is by definition a
     *  shared-range numeric gradient, so this FORCES continuous mode even when this column's own values happened
     *  to auto-detect as categorical (e.g. after a collapse left it a single distinct value) -- otherwise that
     *  one column would render on a categorical palette and break the shared scale. */
    void setSharedRange( final double min, final double max ) {
        _gradient = true;
        _min = min;
        _max = max;
    }

    /** Parses a value as a finite number, or {@code null} if empty, non-numeric, or non-finite (NaN/Infinity).
     * Rejecting non-finite values keeps the gradient min/max, the color, and the bar fraction well-defined
     * even when a column carries a "NaN"/"Infinity" sentinel or an overflowing number. */
    static Double parseNumber( final String s ) {
        if ( ForesterUtil.isEmpty( s ) ) {
            return null;
        }
        try {
            final double d = Double.parseDouble( s.trim() );
            return ( Double.isNaN( d ) || Double.isInfinite( d ) ) ? null : Double.valueOf( d );
        }
        catch ( final NumberFormatException e ) {
            return null;
        }
    }

    private static String formatNumber( final double d ) {
        return ( d == Math.rint( d ) ) ? Long.toString( (long) d ) : Double.toString( d );
    }

    /**
     * Whether the values of this ref should be colored as a continuous gradient rather than as distinct
     * categories: the column is <em>predominantly numeric</em> -- a strict majority of the non-empty visible
     * values parse as finite numbers, and there are at least two distinct numbers (so there is a real range).
     * A few non-numeric sentinels ("n/a", "unknown") are tolerated -- they simply get no color -- while a
     * mostly-textual column with a stray number stays categorical.
     */
    static boolean shouldUseGradient( final List<PhylogenyNode> leaves, final String ref ) {
        int total = 0;
        int numeric = 0;
        final Set<Double> distinct = new HashSet<Double>();
        for( final PhylogenyNode node : leaves ) {
            final String v = valueFor( node, ref );
            if ( !ForesterUtil.isEmpty( v ) ) {
                ++total;
                final Double d = parseNumber( v );
                if ( d != null ) {
                    ++numeric;
                    distinct.add( d );
                }
            }
        }
        return ( ( numeric * 2 ) > total ) && ( distinct.size() >= 2 );
    }

    /**
     * The delimiter at which a value of this ref is truncated before grouping (dropping a
     * trailing qualifier), or {@code 0} for refs whose whole value is used. A {@code country}
     * value keeps only the part before the first {@code :} (the subdivision); a {@code host}
     * value keeps only the part before the first {@code ;} (sex/age qualifiers). Matched on the
     * ref name in any namespace.
     */
    private static char truncationDelimiter( final String ref ) {
        if ( refNameEquals( ref, "country" ) ) {
            return ':';
        }
        if ( refNameEquals( ref, "host" ) ) {
            return ';';
        }
        return 0;
    }

    private static boolean refNameEquals( final String ref, final String name ) {
        if ( ForesterUtil.isEmpty( ref ) ) {
            return false;
        }
        final int colon = ref.lastIndexOf( ':' );
        return ( ( colon >= 0 ) ? ref.substring( colon + 1 ) : ref ).equalsIgnoreCase( name );
    }

    /**
     * A user-friendly display label for a property ref: the namespace prefix is dropped,
     * underscores become spaces, and the first letter of each word is capitalized
     * (e.g. {@code repseq:protein_names} becomes {@code Protein Names}). The rest of each
     * word is left untouched so acronyms such as {@code RNA} are preserved. For display
     * only -- the underlying ref and the stored property values are never modified.
     */
    static String displayName( final String ref ) {
        if ( ForesterUtil.isEmpty( ref ) ) {
            return ref;
        }
        final int colon = ref.lastIndexOf( ':' );
        final String name = ( colon >= 0 ) ? ref.substring( colon + 1 ) : ref;
        final StringBuilder sb = new StringBuilder( name.length() );
        boolean start_of_word = true;
        for( int i = 0; i < name.length(); ++i ) {
            char c = name.charAt( i );
            if ( c == '_' ) {
                c = ' ';
            }
            if ( Character.isWhitespace( c ) ) {
                start_of_word = true;
                sb.append( c );
            }
            else if ( start_of_word ) {
                sb.append( Character.toUpperCase( c ) );
                start_of_word = false;
            }
            else {
                sb.append( c );
            }
        }
        return sb.toString();
    }

    /**
     * The external nodes actually visible in the given (sub)tree: leaves hidden underneath a
     * collapsed node are excluded, so collapsing a clade removes its values from the coloring
     * and legend. A collapsed node is itself internal and carries no leaf value, so the
     * collapsed clade contributes nothing. A {@code null}/empty phylogeny yields no leaves.
     */
    static List<PhylogenyNode> visibleExternalNodes( final Phylogeny phylogeny ) {
        final List<PhylogenyNode> leaves = new ArrayList<PhylogenyNode>();
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return leaves;
        }
        final Deque<PhylogenyNode> stack = new ArrayDeque<PhylogenyNode>();
        stack.push( phylogeny.getRoot() );
        while ( !stack.isEmpty() ) {
            final PhylogenyNode n = stack.pop();
            if ( n.isExternal() ) {
                leaves.add( n );
            }
            else if ( !n.isCollapse() ) {
                for( final PhylogenyNode child : n.getDescendants() ) {
                    stack.push( child );
                }
            }
        }
        return leaves;
    }

    static String valueFor( final PhylogenyNode node, final String ref ) {
        if ( ( node.getNodeData() == null ) || ( node.getNodeData().getProperties() == null ) ) {
            return null;
        }
        final List<Property> props = node.getNodeData().getProperties().getProperties( ref );
        return props.isEmpty() ? null : props.get( 0 ).getValue();
    }

    /**
     * The property references worth coloring by: those present on the tree's external
     * nodes with between two distinct values and one less than the number of leaves
     * (constant and per-leaf-unique properties are useless to color by). The explicit
     * {@code style:*} visualization properties are excluded. Sorted alphabetically.
     */
    static List<String> colorableRefs( final Phylogeny phylogeny ) {
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return new ArrayList<String>();
        }
        final Map<String, Set<String>> ref_to_values = new TreeMap<String, Set<String>>();
        final Map<String, int[]> ref_to_counts = new HashMap<String, int[]>();       // [total non-empty, numeric]
        final Map<String, Set<Double>> ref_to_numbers = new HashMap<String, Set<Double>>();
        for( final PhylogenyNode node : phylogeny.getExternalNodes() ) {
            if ( ( node.getNodeData() != null ) && ( node.getNodeData().getProperties() != null ) ) {
                for( final Property p : node.getNodeData().getProperties().getProperties() ) {
                    final String ref = p.getRef();
                    if ( ref.startsWith( NodeVisualData.APTX_VISUALIZATION_REF ) ) {
                        continue;
                    }
                    Set<String> vs = ref_to_values.get( ref );
                    if ( vs == null ) {
                        vs = new HashSet<String>();
                        ref_to_values.put( ref, vs );
                        ref_to_counts.put( ref, new int[ 2 ] );
                        ref_to_numbers.put( ref, new HashSet<Double>() );
                    }
                    vs.add( p.getValue() );
                    if ( !ForesterUtil.isEmpty( p.getValue() ) ) {
                        ref_to_counts.get( ref )[ 0 ]++;
                        final Double d = parseNumber( p.getValue() );
                        if ( d != null ) {
                            ref_to_counts.get( ref )[ 1 ]++;
                            ref_to_numbers.get( ref ).add( d );
                        }
                    }
                }
            }
        }
        final int leaves = phylogeny.getExternalNodes().size();
        final List<String> refs = new ArrayList<String>();
        for( final Map.Entry<String, Set<String>> e : ref_to_values.entrySet() ) {
            final String ref = e.getKey();
            final int distinct = e.getValue().size();
            if ( distinct < 2 ) {
                continue; // constant column: nothing to distinguish
            }
            final int[] c = ref_to_counts.get( ref );
            final boolean gradient = ( ( c[ 1 ] * 2 ) > c[ 0 ] ) && ( ref_to_numbers.get( ref ).size() >= 2 );
            // Drop per-leaf-unique CATEGORICAL columns (one color per tip is useless to color by); a numeric
            // column stays colorable even when every tip has a distinct value, because it renders as a gradient.
            if ( gradient || ( distinct < leaves ) ) {
                refs.add( ref );
            }
        }
        return refs;
    }

    /**
     * The colorable refs whose (visible) values are predominantly numeric -- i.e. the ones that can be SIZED by
     * ("Size by" scales a node symbol by the value). A subset of {@link #colorableRefs(Phylogeny)}; the same
     * "predominantly numeric with a real range" test the gradient coloring uses.
     */
    static List<String> numericRefs( final Phylogeny phylogeny ) {
        return numericRefs( phylogeny, colorableRefs( phylogeny ) );
    }

    /** As {@link #numericRefs(Phylogeny)} but reusing an already-computed {@code colorable} list, so a caller that
     *  also needs {@link #colorableRefs(Phylogeny)} does not scan the whole tree's properties twice. */
    static List<String> numericRefs( final Phylogeny phylogeny, final List<String> colorable ) {
        final List<PhylogenyNode> leaves = visibleExternalNodes( phylogeny );
        final List<String> refs = new ArrayList<String>();
        for( final String ref : colorable ) {
            if ( shouldUseGradient( leaves, ref ) ) {
                refs.add( ref );
            }
        }
        return refs;
    }
}
