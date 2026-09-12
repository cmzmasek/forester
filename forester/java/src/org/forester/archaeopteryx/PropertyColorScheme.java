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
import java.util.regex.Pattern;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeData;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
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
    // whether _ref is a reserved element slot (taxonomy/sequence field): values verbatim, no grouping
    private final boolean            _element_slot;
    // the ref's visualization candidate when the scheme was built for one: values are read as it grouped them
    private final VisCandidate       _candidate;
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
        this( phylogeny, ref, overrides, palette_name, memory, memory_next, null );
    }

    /**
     * @param forced_gradient non-null forces the mode -- TRUE gradient, FALSE categorical -- overriding the
     *                        automatic numeric detection. The Color-by path always passes it (the three-band
     *                        default plus the user's per-field {@code [colors]}/{@code [gradient]} choice);
     *                        null keeps the plain detection (annotation columns, legacy callers).
     */
    PropertyColorScheme( final Phylogeny phylogeny, final String ref, final Map<String, Color> overrides,
                         final String palette_name, final Map<String, Color> memory, final int[] memory_next,
                         final Boolean forced_gradient ) {
        this( phylogeny, ref, overrides, palette_name, memory, memory_next, forced_gradient, null );
    }

    /**
     * @param candidate the ref's visualization CANDIDATE (see {@link #visualizationCandidates(Phylogeny)}), or null.
     *                  When given, a node's value is read exactly as the candidate grouped it
     *                  ({@link #visualizationNodeValue}): folded and mapped to its group's representative, so a
     *                  numeric field's "1" and "1.0" are one legend row and one colour. Null keeps the plain per-ref
     *                  reading (annotation columns, and refs that are not candidates).
     */
    PropertyColorScheme( final Phylogeny phylogeny, final String ref, final Map<String, Color> overrides,
                         final String palette_name, final Map<String, Color> memory, final int[] memory_next,
                         final Boolean forced_gradient, final VisCandidate candidate ) {
        _ref = ref;
        _candidate = candidate;
        _palette = paletteByName( palette_name );
        _truncate_at = truncationDelimiter( ref );
        _element_slot = isElementSlot( ref );
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
            final String v = nodeValue( node );
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
        _gradient = ( forced_gradient != null ) ? forced_gradient.booleanValue()
                : shouldUseGradient( leaves, ref );
        if ( _gradient ) {
            double min = Double.POSITIVE_INFINITY;
            double max = Double.NEGATIVE_INFINITY;
            for( final PhylogenyNode node : leaves ) {
                final Double d = parseNumber( nodeValue( node ) );
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
                final String v = nodeValue( node );
                if ( !ForesterUtil.isEmpty( v ) ) {
                    final String label = displayLabel( v );
                    if ( label.isEmpty() ) {
                        continue; // a value that folds to nothing (e.g. "_") forms no group (JS rule)
                    }
                    final String key = _element_slot ? label : label.toLowerCase( Locale.ROOT );
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
                // element-slot labels stay verbatim (one spelling per group by construction)
                groups.add( new String[] { _element_slot ? e.getKey() : representative( e.getValue() ),
                                           e.getKey() } );
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
            final Double d = parseNumber( nodeValue( node ) );
            if ( d == null ) {
                return null;
            }
            final double t = ( _max > _min ) ? ( ( d - _min ) / ( _max - _min ) ) : 0.0;
            return gradientColorAt( t );
        }
        final String v = nodeValue( node );
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
        final Double d = parseNumber( nodeValue( node ) );
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
        return ( _candidate != null ) ? v : foldLabel( v, _truncate_at, _element_slot ); // a candidate value is folded already
    }

    /** A node's value of this scheme's field: as its candidate grouped it, else the plain first value. */
    private String nodeValue( final PhylogenyNode node ) {
        return ( _candidate != null ) ? visualizationNodeValue( node, _candidate ) : valueFor( node, _ref );
    }

    /** The display form a raw value is grouped under, for a scheme built without a candidate. */
    private static String foldLabel( final String v, final char truncate_at, final boolean element_slot ) {
        if ( element_slot ) {
            return v; // element-slot values are used VERBATIM (no cut/fold/dictionary) -- the JS rule
        }
        String s = v;
        if ( truncate_at != 0 ) {
            final int idx = s.indexOf( truncate_at );
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

    /** The normalized key a value is grouped/colored by: its display label, case-folded -- except for an
     *  element slot, whose values (hence keys) are verbatim (JS parity: "verbatim for element slots"). */
    private String groupKey( final String v ) {
        return _element_slot ? v : displayLabel( v ).toLowerCase( Locale.ROOT );
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
    /** The numeric three-band system (JS parity), deciding the COLOR-BY default mode and switchability:
     *  NOT_NUMERIC; SMALL (<= 10 distinct numbers -- numbers that few are usually codes, like HA/NA
     *  subtypes: individual colors by default, switchable); MEDIUM (11..20: gradient by default,
     *  switchable); LARGE (> 20: gradient, fixed). Only the Color-by path uses it -- the annotation-column
     *  schemes keep the plain numeric default (a BAR/HEATMAP column must stay continuous regardless). */
    enum ModeBand {
        NOT_NUMERIC, SMALL, MEDIUM, LARGE;

        boolean isSwitchable() {
            return ( this == SMALL ) || ( this == MEDIUM );
        }

        boolean defaultsToGradient() {
            return ( this == MEDIUM ) || ( this == LARGE );
        }
    }

    static ModeBand colorModeBand( final List<PhylogenyNode> leaves, final String ref ) {
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
        if ( ( ( numeric * 2 ) <= total ) || ( distinct.size() < 2 ) ) {
            return ModeBand.NOT_NUMERIC;
        }
        if ( distinct.size() <= 10 ) {
            return ModeBand.SMALL;
        }
        return ( distinct.size() <= 20 ) ? ModeBand.MEDIUM : ModeBand.LARGE;
    }

    /** Whether the field would DEFAULT to a gradient in the Color-by band system (candidate tiering). */
    private static boolean defaultsToGradient( final List<PhylogenyNode> leaves, final String ref ) {
        return colorModeBand( leaves, ref ).defaultsToGradient();
    }

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
     * The name the menus and legends show for a ref -- a PORT of Archaeopteryx.js {@code propertyDisplayName} /
     * {@code prettifyVisLabel} (forester.js 79dd9de), character for character, because the exclusion rules match
     * THIS string and the two programs are diffed against a shared fixture. The namespace is everything up to the
     * FIRST colon; underscores become spaces; camelCase splits on lowercase->uppercase and on a digit before an
     * uppercase+lowercase pair (so "H5N1" stays whole and "GlobalH1Clade" reads "Global H1 Clade"); a word whose
     * first character is a lowercase ASCII letter has that letter capitalised, and every other word is left exactly
     * as written. A hyphen is not a separator. The six element slots have fixed names. Display only -- the ref and
     * the stored values are never modified.
     */
    static String displayName( final String ref ) {
        if ( ForesterUtil.isEmpty( ref ) ) {
            return ref;
        }
        // the element slots carry the same display labels as the Archaeopteryx.js Color menu
        if ( TAX_CODE_REF.equals( ref ) ) {
            return "Taxonomy Code";
        }
        if ( TAX_SCI_NAME_REF.equals( ref ) ) {
            return "Scientific Name";
        }
        if ( TAX_COMMON_NAME_REF.equals( ref ) ) {
            return "Common Name";
        }
        if ( SEQ_NAME_REF.equals( ref ) ) {
            return "Sequence Name";
        }
        if ( SEQ_SYMBOL_REF.equals( ref ) ) {
            return "Sequence Symbol";
        }
        if ( SEQ_GENE_NAME_REF.equals( ref ) ) {
            return "Gene Name";
        }
        return prettifyVisLabel( localName( ref ) );
    }

    /** The ref minus its namespace: everything after the FIRST colon, as forester.js reads it
     *  ({@code ref.indexOf(':')}); the whole ref when there is none. */
    static String localName( final String ref ) {
        final int colon = ref.indexOf( ':' );
        return ( colon >= 0 ) ? ref.substring( colon + 1 ) : ref;
    }

    /** forester.js {@code prettifyVisLabel}, verbatim -- see {@link #displayName(String)}. Splits on single spaces
     *  and keeps empty words, exactly like the JS {@code split(' ')}, so runs of spaces survive as written. */
    static String prettifyVisLabel( final String name ) {
        final String s = name.replace( '_', ' ' ).replaceAll( "([a-z])([A-Z])", "$1 $2" )
                .replaceAll( "([0-9])([A-Z][a-z])", "$1 $2" );
        final String[] words = s.split( " ", -1 );
        final StringBuilder sb = new StringBuilder( s.length() );
        for( int i = 0; i < words.length; ++i ) {
            if ( i > 0 ) {
                sb.append( ' ' );
            }
            final String w = words[ i ];
            if ( !w.isEmpty() && ( w.charAt( 0 ) >= 'a' ) && ( w.charAt( 0 ) <= 'z' ) ) {
                sb.append( Character.toUpperCase( w.charAt( 0 ) ) ).append( w, 1, w.length() );
            }
            else {
                sb.append( w );
            }
        }
        return sb.toString();
    }

    /** The displayed name split into lower-case words: every run of non-alphanumerics is one break (forester.js
     *  {@code visNameWords}). The word rules below read this, so a rule and the label it judges cannot drift. */
    private static String visNameWords( final String ref ) {
        return prettifyVisLabel( localName( ref ) ).toLowerCase( Locale.ROOT ).replaceAll( "[^a-z0-9]+", " " ).trim();
    }

    private static final Pattern   VIS_EXCLUDED_LOCAL_NAME_RE = Pattern.compile( "(taxonomy|taxon|tax)id$" );
    private static final Pattern[] VIS_EXCLUDED_WORD_RES      = { Pattern.compile( "(^| )authors?( |$)" ),
            Pattern.compile( "(^| )set( |$)" ), Pattern.compile( "(^| )data use( |$)" ),
            Pattern.compile( "(^| )ids?( |$)" ), Pattern.compile( "accessions?$" ),
            Pattern.compile( "identifiers?$" ), Pattern.compile( "(^| )restricted ?until( |$)" ) };
    private static final String[]  VIS_EXCLUDED_WORD_REASONS  = { "author", "set", "data use", "id", "accession",
            "identifier", "restricted until" };
    private static final Pattern[] VIS_DEPRIORITIZED_WORD_RES = { Pattern.compile( "(^| )(in|out) groups?( |$)" ),
            Pattern.compile( "(^| )(in|out)groups?( |$)" ) };

    /**
     * Why a property ref is never offered for visualization, or {@code null} when it may be -- forester.js
     * {@code visExcludedRef}, verbatim. These describe the RECORD rather than the organism (who deposited it, which
     * collection it belongs to, what may be done with it, what a database calls it): they repeat like categories and
     * so pass every statistical test, but a colour spent on one says nothing about the tree.
     * <ul>
     * <li>the {@code style:} namespace (per-node rendering instructions, not data);</li>
     * <li>a taxon id: the LOCAL name, lower-cased with every non-alphanumeric removed, ends in taxonomyid / taxonid /
     * taxid (so vipr:NCBI_Taxon_Id, ncbi_taxid, taxon_id all count);</li>
     * <li>on the DISPLAYED name split into words: the WORDS author(s), set, data use and id(s) -- word-anchored, so
     * "Dataset", "Authority", "Plasmid", "Hybrid" and "Lipid" survive -- and the SUFFIXES accession(s) and
     * identifier(s), which are long enough that a name ending in those letters is one ("GBAccession").</li>
     * </ul>
     * Matching the displayed name rather than the raw ref is load-bearing: {@code dataUseTerms} reaches the user as
     * "Data Use Terms". CROSS-IMPLEMENTATION CONTRACT, JS-authoritative (Christian, 2026-09-12): pinned by the JS
     * repo's test/fixtures/vis-contract.tsv; never retune it here alone.
     */
    static String excludedRefReason( final String ref ) {
        if ( ForesterUtil.isEmpty( ref ) ) {
            return null;
        }
        if ( ref.startsWith( "style:" ) ) {
            return "style";
        }
        if ( VIS_EXCLUDED_LOCAL_NAME_RE.matcher( localName( ref ).toLowerCase( Locale.ROOT ).replaceAll( "[^a-z0-9]",
                                                                                                          "" ) )
                .find() ) {
            return "taxon id";
        }
        final String words = visNameWords( ref );
        for( int i = 0; i < VIS_EXCLUDED_WORD_RES.length; ++i ) {
            if ( VIS_EXCLUDED_WORD_RES[ i ].matcher( words ).find() ) {
                return VIS_EXCLUDED_WORD_REASONS[ i ];
            }
        }
        return null;
    }

    /** Whether a ref must never be offered for visualization. See {@link #excludedRefReason(String)}. */
    static boolean isExcludedRef( final String ref ) {
        return excludedRefReason( ref ) != null;
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
        if ( isElementSlot( ref ) ) { // reserved refs take precedence over a like-named property
            return elementValue( node, ref );
        }
        if ( ( node.getNodeData() == null ) || ( node.getNodeData().getProperties() == null ) ) {
            return null;
        }
        final List<Property> props = node.getNodeData().getProperties().getProperties( ref );
        return props.isEmpty() ? null : props.get( 0 ).getValue();
    }

    // ---- ELEMENT SLOTS (JS parity): the phyloXML taxonomy/sequence fields offered in "Color by" alongside
    //      the properties, under the SAME reserved ids Archaeopteryx.js uses (VIS_ELEMENT_SLOTS in
    //      forester.js), so a shared tree colors by the same fields in both viewers. Their values are used
    //      VERBATIM -- no qualifier cut, no spelling fold, no synonym dictionary (a common name "swine" must
    //      NOT become a "Pig" legend row when coloring by the taxonomy field itself; the JS rule). ----
    static final String   TAX_CODE_REF        = "tax:code";
    static final String   TAX_SCI_NAME_REF    = "tax:scientific_name";
    static final String   TAX_COMMON_NAME_REF = "tax:common_name";
    static final String   SEQ_NAME_REF        = "seq:name";
    static final String   SEQ_SYMBOL_REF      = "seq:symbol";
    static final String   SEQ_GENE_NAME_REF   = "seq:gene_name";
    static final String[] ELEMENT_SLOT_REFS   = { TAX_CODE_REF, TAX_SCI_NAME_REF, TAX_COMMON_NAME_REF,
                                                  SEQ_NAME_REF, SEQ_SYMBOL_REF, SEQ_GENE_NAME_REF };

    /** Whether {@code ref} is one of the reserved element-slot ids above. */
    static boolean isElementSlot( final String ref ) {
        if ( ref == null ) {
            return false;
        }
        for( final String r : ELEMENT_SLOT_REFS ) {
            if ( r.equals( ref ) ) {
                return true;
            }
        }
        return false;
    }

    /** The node's value for an element slot (null when the element or the field is absent/empty). */
    private static String elementValue( final PhylogenyNode node, final String ref ) {
        if ( node.getNodeData() == null ) {
            return null;
        }
        String v = null;
        if ( ref.startsWith( "tax:" ) ) {
            if ( node.getNodeData().isHasTaxonomy() ) {
                final org.forester.phylogeny.data.Taxonomy t = node.getNodeData().getTaxonomy();
                v = TAX_CODE_REF.equals( ref ) ? t.getTaxonomyCode()
                        : ( TAX_SCI_NAME_REF.equals( ref ) ? t.getScientificName() : t.getCommonName() );
            }
        }
        else if ( node.getNodeData().isHasSequence() ) {
            final org.forester.phylogeny.data.Sequence q = node.getNodeData().getSequence();
            v = SEQ_NAME_REF.equals( ref ) ? q.getName()
                    : ( SEQ_SYMBOL_REF.equals( ref ) ? q.getSymbol() : q.getGeneName() );
        }
        return ForesterUtil.isEmpty( v ) ? null : v;
    }

    // =====================================================================================================
    // AUTOMATIC VISUALIZATION CANDIDATES -- a PORT of Archaeopteryx.js, forester.js "Automatic visualization
    // candidates" (commit 79dd9de). Christian, 2026-09-12: "The JS rules are the authoritative ones. Desktop will
    // have to follow them 100%." The acceptance test is the JS repo's pair of fixtures -- test/fixtures/
    // vis-contract.tsv (names) and vis-trees.tsv + vis-trees/*.xml (data) -- copied into
    // forester/test_data/vis_contract/ and diffed by PropertyColorSchemeTest. Never retune a rule here alone.
    //
    // Candidacy is decided on the TREE: when it is loaded, and again only after it is EDITED. A VIEW -- a subtree,
    // a collapse -- never re-decides it; it only re-summarizes a candidate over the tips on screen
    // (visualizationSummary). A clade is by nature a set of tips sharing a value, and one value is refused, so
    // re-classifying per view would drop the chosen colouring in most clades.
    //
    // THE REFUSAL RULES DECIDE WHAT IS OFFERED, NEVER WHAT IS ALREADY CHOSEN: after an edit, a field the user had
    // chosen that still carries a value on any tip is KEPT, even where the rules would now refuse it
    // (visualizationCandidatesKeeping).
    //
    // The rules, per candidate, in order (external nodes only):
    //   multi-valued  a ref carried more than once by any one tip -> refused: a node cannot be two colours.
    //   numeric       every grouped value matches VIS_NUMERIC_RE -> numeric; spellings of one number fold.
    //   repetition    fewer than 2 distinct values -> refused.
    //   coverage      on fewer than 2/3 of the tips -> SPARSE: offered, ranked low, never refused.
    //   categorical   as many distinct values as the tree has TIPS (all tips, not the covered ones) -> refused as
    //                 an identifier; more than 20 -> WIDE (offered, never opens a tree), and if distinct/covered
    //                 > 3/5 also NEAR-UNIQUE (the very bottom).
    //   numeric       no uniqueness test of any kind: a measurement is naturally one value per sample.
    //   shape         <= 7 distinct values, numeric or not.
    // Tiers, best first: 0 clean categorical, 1 every numeric, 2 wide, 3 in/out-group, 4 sparse, 5 near-unique;
    // within a tier by score = coverage x normalised entropy, then label, then id. The tree OPENS with the first
    // candidate that is not wide.
    // =====================================================================================================

    static final int             VIS_MIN_COVERAGE_NUM     = 2;
    static final int             VIS_MIN_COVERAGE_DEN     = 3;
    static final int             VIS_MAX_COLOR_CATEGORIES = 20;
    static final int             VIS_MAX_SHAPE_CATEGORIES = 7;
    static final int             VIS_NUMERIC_CATEGORY_MAX = 10;
    static final int             VIS_WIDE_REPEAT_NUM      = 3;
    static final int             VIS_WIDE_REPEAT_DEN      = 5;
    /** What "numeric" means, spelled out: an optional sign, decimal digits with an optional fraction, an optional
     *  exponent. PINNED as a grammar because a host language's own parser is not portable -- Java's parseDouble
     *  accepts "Infinity" and "NaN", JavaScript's Number() accepts "0x1A" and "0b101" -- and a field of either would
     *  be a gradient in one program and a category in the other. Matched against the whole (trimmed) value. */
    private static final Pattern VIS_NUMERIC_RE           = Pattern
            .compile( "[+-]?(\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+)?" );

    /** Whether a (trimmed) value is a number under the pinned grammar -- see {@link #VIS_NUMERIC_RE}. */
    static boolean isVisNumber( final String v ) {
        return ( v != null ) && VIS_NUMERIC_RE.matcher( v ).matches();
    }

    /** The number a grammar-valid spelling denotes, as a fold key: "1", "1.0" and "+1" meet here. Minus zero is
     *  zero, as {@code String(Number("-0"))} is "0" in forester.js. */
    private static Double visNumberKey( final String v ) {
        final double d = Double.parseDouble( v );
        return ( d == 0.0 ) ? 0.0 : d;
    }

    /** JavaScript's whitespace set ({@code String.prototype.trim} and {@code \s}), which is wider than Java's: it
     *  includes the no-break space and the other Unicode spaces a spreadsheet export carries. */
    private static final String  JS_WS                     = "\\t\\n\\u000B\\f\\r \\u00A0\\u1680\\u2000-\\u200A"
            + "\\u2028\\u2029\\u202F\\u205F\\u3000\\uFEFF";
    private static final Pattern JS_WS_RUN                 = Pattern.compile( "[" + JS_WS + "]+" );
    private static final Pattern JS_TRAILING_PARENTHETICAL = Pattern
            .compile( "[" + JS_WS + "]*\\([^()]*\\)[" + JS_WS + "]*\\z" );

    private static boolean isJsWhitespace( final char c ) {
        switch ( c ) {
            case '\t':
            case '\n':
            case '':
            case '\f':
            case '\r':
            case ' ':
            case ' ':
            case ' ':
            case ' ':
            case ' ':
            case ' ':
            case ' ':
            case '　':
            case '﻿':
                return true;
            default:
                return ( c >= ' ' ) && ( c <= ' ' );
        }
    }

    /** {@code String.prototype.trim}. */
    static String jsTrim( final String s ) {
        int b = 0;
        int e = s.length();
        while ( ( b < e ) && isJsWhitespace( s.charAt( b ) ) ) {
            ++b;
        }
        while ( ( e > b ) && isJsWhitespace( s.charAt( e - 1 ) ) ) {
            --e;
        }
        return s.substring( b, e );
    }

    /** The trimmed value, or null when there is none -- how forester.js reads a value before anything else. */
    private static String jsClean( final String v ) {
        if ( v == null ) {
            return null;
        }
        final String s = jsTrim( v );
        return s.isEmpty() ? null : s;
    }

    /**
     * The display form a raw PROPERTY value is grouped under -- forester.js {@code visDisplayLabel}, verbatim. For
     * refs named exactly "country" / "host" a trailing qualifier is cut at the first ':' / ';', and a cut that lands
     * inside a parenthetical is trimmed back to before the first unclosed '('. Then trimmed, underscores read as
     * spaces, whitespace runs collapsed, trimmed again (so "_" folds to empty and is dropped), and the synonym
     * dictionary applied as a WHOLE-VALUE match, retried once without a trailing parenthetical. Case is preserved
     * here; grouping lower-cases it. Element-slot values never come through here -- they are verbatim.
     */
    static String visDisplayLabel( final String value, final char cut ) {
        String s = value;
        if ( cut != 0 ) {
            final int at = s.indexOf( cut );
            if ( at >= 0 ) {
                s = s.substring( 0, at );
                final List<Integer> open = new ArrayList<Integer>();
                for( int i = 0; i < s.length(); ++i ) {
                    final char c = s.charAt( i );
                    if ( c == '(' ) {
                        open.add( i );
                    }
                    else if ( ( c == ')' ) && !open.isEmpty() ) {
                        open.remove( open.size() - 1 );
                    }
                }
                if ( !open.isEmpty() ) {
                    s = s.substring( 0, open.get( 0 ) );
                }
            }
        }
        s = jsTrim( JS_WS_RUN.matcher( jsTrim( s ).replace( '_', ' ' ) ).replaceAll( " " ) );
        String hit = SYNONYM_LOOKUP.get( s.toLowerCase( Locale.ROOT ) );
        if ( hit == null ) {
            final String stripped = JS_TRAILING_PARENTHETICAL.matcher( s ).replaceFirst( "" );
            if ( !stripped.equals( s ) && !stripped.isEmpty() ) {
                hit = SYNONYM_LOOKUP.get( stripped.toLowerCase( Locale.ROOT ) );
            }
        }
        return ( hit != null ) ? hit : s;
    }

    /**
     * Whether a property describes the node itself: phyloXML applies_to "node", or "clade" -- the clade rooted at an
     * external node IS that node (the repseq pipeline writes clade for every field). "parent_branch" is out: that is
     * the branch above the node. forester.js {@code isNodeScopedProperty}.
     */
    static boolean isNodeScopedProperty( final Property p ) {
        return ( p.getAppliesTo() == Property.AppliesTo.NODE ) || ( p.getAppliesTo() == Property.AppliesTo.CLADE );
    }

    /**
     * Offered, but ranked after the numerics: "In-Group" and "Out-Group" say which tips were the study set and which
     * were there to root the tree -- a fact about the ANALYSIS, which the person who rooted it already knows -- and
     * they are typically an even two-value split at full coverage, the shape that would otherwise win. Matched as
     * WORDS of the displayed name ("Within Group" contains "in group" and must keep leading), in both hyphenations,
     * both one-word spellings and the plurals. forester.js {@code visDeprioritizedRef}.
     */
    static boolean isDeprioritizedRef( final String ref ) {
        if ( ForesterUtil.isEmpty( ref ) ) {
            return false;
        }
        final String words = visNameWords( ref );
        for( final Pattern p : VIS_DEPRIORITIZED_WORD_RES ) {
            if ( p.matcher( words ).find() ) {
                return true;
            }
        }
        return false;
    }

    /**
     * One visualization candidate: the desktop form of the object forester.js {@code visualizationCandidates}
     * returns. The summary fields (values, counts, coverage, total and a numeric field's colour mode) describe the
     * whole tree when built, and are re-set by {@link #visualizationCandidatesKeeping} for a KEPT field.
     */
    static final class VisCandidate {

        final String              _id;            // "prop:" + ref, or the element-slot id -- as in forester.js
        final String              _ref;           // the property ref or slot id: what menus and schemes key on
        final boolean             _element_slot;
        String                    _label;
        final boolean             _numeric;
        int                       _coverage;
        int                       _total;
        List<String>              _values;
        Map<String, Integer>      _counts;
        final Map<String, String> _canon;         // group key -> representative; null for a non-numeric slot
        final char                _cut;
        final double              _score;
        String                    _color_mode;    // "category" or "range"
        boolean                   _switchable;
        final boolean             _wide;
        final boolean             _near_unique;
        final boolean             _sparse;
        final boolean             _deprioritized;
        final boolean             _shape;
        boolean                   _kept;

        VisCandidate( final String id,
                      final String ref,
                      final boolean element_slot,
                      final String label,
                      final boolean numeric,
                      final int coverage,
                      final int total,
                      final List<String> values,
                      final Map<String, Integer> counts,
                      final Map<String, String> canon,
                      final char cut,
                      final double score,
                      final String color_mode,
                      final boolean switchable,
                      final boolean wide,
                      final boolean near_unique,
                      final boolean sparse,
                      final boolean deprioritized,
                      final boolean shape ) {
            _id = id;
            _ref = ref;
            _element_slot = element_slot;
            _label = label;
            _numeric = numeric;
            _coverage = coverage;
            _total = total;
            _values = values;
            _counts = counts;
            _canon = canon;
            _cut = cut;
            _score = score;
            _color_mode = color_mode;
            _switchable = switchable;
            _wide = wide;
            _near_unique = near_unique;
            _sparse = sparse;
            _deprioritized = deprioritized;
            _shape = shape;
        }

        /** forester.js {@code tierOf}: the FIRST flag that applies, in this precedence. */
        int tier() {
            if ( _near_unique ) {
                return 5;
            }
            if ( _sparse ) {
                return 4;
            }
            if ( _deprioritized ) {
                return 3;
            }
            if ( _numeric ) {
                return 1;
            }
            return _wide ? 2 : 0;
        }
    }

    /** What a view shows of a candidate: forester.js {@code visualizationSummary}. */
    static final class VisSummary {

        final List<String>         _values;
        final Map<String, Integer> _counts;
        final int                  _coverage;
        final int                  _total;
        final int                  _distinct;
        final String               _color_mode;  // a numeric candidate's band in this view; null for a category
        final boolean              _switchable;

        VisSummary( final List<String> values,
                    final Map<String, Integer> counts,
                    final int coverage,
                    final int total,
                    final String color_mode,
                    final boolean switchable ) {
            _values = values;
            _counts = counts;
            _coverage = coverage;
            _total = total;
            _distinct = values.size();
            _color_mode = color_mode;
            _switchable = switchable;
        }
    }

    private static final class VisGroup {

        int                        _count;
        final Map<String, Integer> _spellings = new HashMap<String, Integer>();
    }

    private static final class VisStats {

        final String                _ref;
        final boolean               _element_slot;
        final char                  _cut;
        int                         _nodes;
        boolean                     _multi;
        final Map<String, VisGroup> _keys = new LinkedHashMap<String, VisGroup>();

        VisStats( final String ref, final boolean element_slot, final char cut ) {
            _ref = ref;
            _element_slot = element_slot;
            _cut = cut;
        }
    }

    /** Best first: tier, then score (higher first), then label ignoring case, then id. forester.js's sort. */
    static final Comparator<VisCandidate> VIS_ORDER = new Comparator<VisCandidate>() {

        @Override
        public int compare( final VisCandidate a, final VisCandidate b ) {
            final int ta = a.tier();
            final int tb = b.tier();
            if ( ta != tb ) {
                return ta - tb;
            }
            if ( a._score != b._score ) {
                return Double.compare( b._score, a._score );
            }
            final int by_label = a._label.toLowerCase( Locale.ROOT ).compareTo( b._label.toLowerCase( Locale.ROOT ) );
            return ( by_label != 0 ) ? by_label : a._id.compareTo( b._id );
        }
    };

    /**
     * The visualization candidates of a tree, best first -- forester.js {@code visualizationCandidates}, verbatim.
     * See the rules at the top of this section. Candidates come from the six element slots (verbatim values, never
     * folded) and from node-scoped properties that no name rule excludes. Every external node counts, collapsed or
     * not: this is the TREE, not a view.
     */
    static List<VisCandidate> visualizationCandidates( final Phylogeny tree ) {
        final List<VisCandidate> candidates = new ArrayList<VisCandidate>();
        if ( ( tree == null ) || tree.isEmpty() ) {
            return candidates;
        }
        int total = 0;
        final Map<String, VisStats> stats = new LinkedHashMap<String, VisStats>();
        for( final PhylogenyNode n : allExternalNodes( tree ) ) {
            total++;
            // gather this node's values per candidate first, so carrying the same ref twice is visible as such
            final Map<String, List<String>> per_node = new LinkedHashMap<String, List<String>>();
            final NodeData nd = n.getNodeData();
            if ( nd != null ) {
                for( final String slot : ELEMENT_SLOT_REFS ) {
                    if ( slot.startsWith( "tax:" ) ) {
                        if ( nd.getTaxonomies() != null ) {
                            for( final Taxonomy t : nd.getTaxonomies() ) {
                                addVisValue( per_node, slot, taxonomySlotValue( t, slot ) );
                            }
                        }
                    }
                    else if ( nd.getSequences() != null ) {
                        for( final Sequence q : nd.getSequences() ) {
                            addVisValue( per_node, slot, sequenceSlotValue( q, slot ) );
                        }
                    }
                }
                if ( nd.getProperties() != null ) {
                    for( final Property p : nd.getProperties().getProperties() ) {
                        if ( !ForesterUtil.isEmpty( p.getRef() ) && isNodeScopedProperty( p )
                                && !isExcludedRef( p.getRef() ) ) {
                            addVisValue( per_node, "prop:" + p.getRef(), p.getValue() );
                        }
                    }
                }
            }
            for( final Map.Entry<String, List<String>> e : per_node.entrySet() ) {
                final String id = e.getKey();
                final boolean slot = !id.startsWith( "prop:" );
                final String ref = slot ? id : id.substring( 5 );
                final char cut = slot ? 0 : truncationDelimiter( ref );
                // properties group under their normalized display form; element slots are curated, verbatim. A value
                // that folds to NOTHING ("_", a host that is only its ";" qualifier) is NO VALUE, exactly like an empty
                // string: it does not cover the tip, makes no group, and does not count as "carried twice"
                // (Christian, 2026-09-12; forester.js 635741b)
                final List<String> displays = new ArrayList<String>( e.getValue().size() );
                for( final String v : e.getValue() ) {
                    final String display = slot ? v : visDisplayLabel( v, cut );
                    if ( !display.isEmpty() ) {
                        displays.add( display );
                    }
                }
                if ( displays.isEmpty() ) {
                    continue;
                }
                VisStats s = stats.get( id );
                if ( s == null ) {
                    s = new VisStats( ref, slot, cut );
                    stats.put( id, s );
                }
                s._nodes++;
                if ( displays.size() > 1 ) {
                    s._multi = true;
                }
                for( final String display : displays ) {
                    final String key = s._element_slot ? display : display.toLowerCase( Locale.ROOT );
                    VisGroup g = s._keys.get( key );
                    if ( g == null ) {
                        g = new VisGroup();
                        s._keys.put( key, g );
                    }
                    g._count++;
                    final Integer c = g._spellings.get( display );
                    g._spellings.put( display, ( c == null ) ? 1 : ( c + 1 ) );
                }
            }
        }
        for( final Map.Entry<String, VisStats> e : stats.entrySet() ) {
            final String id = e.getKey();
            final VisStats s = e.getValue();
            if ( s._multi ) {
                continue;
            }
            final int covered = s._nodes;
            // one legend row per group: the most frequent spelling (ties by code point), first letter capitalised
            // for a property
            final Map<String, String> canon = new LinkedHashMap<String, String>();
            Map<String, Integer> counts = new LinkedHashMap<String, Integer>();
            for( final Map.Entry<String, VisGroup> ge : s._keys.entrySet() ) {
                String rep = null;
                int best = -1;
                for( final Map.Entry<String, Integer> sp : ge.getValue()._spellings.entrySet() ) {
                    final int n = sp.getValue();
                    if ( ( n > best ) || ( ( n == best ) && ( sp.getKey().compareTo( rep ) < 0 ) ) ) {
                        rep = sp.getKey();
                        best = n;
                    }
                }
                if ( !s._element_slot ) {
                    rep = rep.substring( 0, 1 ).toUpperCase( Locale.ROOT ) + rep.substring( 1 );
                }
                canon.put( ge.getKey(), rep );
                counts.put( rep, ge.getValue()._count );
            }
            List<String> values = new ArrayList<String>( counts.keySet() );
            boolean numeric = true;
            for( final String v : values ) {
                if ( !isVisNumber( v ) ) {
                    numeric = false;
                    break;
                }
            }
            if ( numeric ) {
                // "1", "1.0" and "+1" are one value: group the spellings by the number they denote; the
                // representative is the SHORTEST spelling, ties by code point (a number has no preferred spelling,
                // so frequency would only pick the exporter's habit). canon then maps every group to it, so
                // visualizationNodeValue folds a node's "1.0" the same way.
                final Map<Double, String> rep_by_number = new LinkedHashMap<Double, String>();
                final Map<Double, Integer> count_by_number = new HashMap<Double, Integer>();
                for( final String v : values ) {
                    final Double k = visNumberKey( v );
                    final String r = rep_by_number.get( k );
                    if ( r == null ) {
                        rep_by_number.put( k, v );
                        count_by_number.put( k, counts.get( v ) );
                    }
                    else {
                        count_by_number.put( k, count_by_number.get( k ) + counts.get( v ) );
                        if ( ( v.length() < r.length() ) || ( ( v.length() == r.length() ) && ( v.compareTo( r ) < 0 ) ) ) {
                            rep_by_number.put( k, v );
                        }
                    }
                }
                final Map<String, Integer> folded = new LinkedHashMap<String, Integer>();
                for( final Map.Entry<Double, String> re : rep_by_number.entrySet() ) {
                    folded.put( re.getValue(), count_by_number.get( re.getKey() ) );
                }
                for( final Map.Entry<String, String> ce : canon.entrySet() ) {
                    ce.setValue( rep_by_number.get( visNumberKey( ce.getValue() ) ) );
                }
                counts = folded;
                values = new ArrayList<String>( counts.keySet() );
            }
            final int distinct = values.size();
            // a sparse field is RANKED LOW, never refused: a half-annotated field is often the interesting one
            final boolean sparse = ( covered * VIS_MIN_COVERAGE_DEN ) < ( total * VIS_MIN_COVERAGE_NUM );
            if ( distinct < 2 ) {
                continue;
            }
            final String color_mode;
            boolean switchable = false;
            boolean wide = false;
            boolean near_unique = false;
            if ( numeric ) {
                // no uniqueness refusal for numbers: a read count, a viral load or a year is naturally one value per
                // sample, which is what a measurement is -- identifiers are caught by the NAME rules
                color_mode = ( distinct <= VIS_NUMERIC_CATEGORY_MAX ) ? "category" : "range";
                switchable = distinct <= VIS_MAX_COLOR_CATEGORIES;
                Collections.sort( values, new Comparator<String>() {

                    @Override
                    public int compare( final String a, final String b ) {
                        return Double.compare( Double.parseDouble( a ), Double.parseDouble( b ) );
                    }
                } );
            }
            else {
                // Counted against TOTAL tips, not the covered ones -- Christian's decision 2026-09-12, left alone
                // until a user complains: a field unique across its annotated subset (1,168 strains over 1,168 of
                // 1,170 tips) is offered at the bottom rather than refused, because counting against covered would
                // also refuse a field carried by two tips with two values, where all-unique is a sample of two.
                if ( distinct >= total ) {
                    continue;
                }
                if ( distinct > VIS_MAX_COLOR_CATEGORIES ) {
                    wide = true;
                    if ( ( distinct * VIS_WIDE_REPEAT_DEN ) > ( covered * VIS_WIDE_REPEAT_NUM ) ) {
                        near_unique = true;
                    }
                }
                color_mode = "category";
                Collections.sort( values );
            }
            // coverage times balance, where balance is the normalised entropy of the value distribution; StrictMath
            // (fdlibm) is what V8 uses too, so equal scores stay equal across the two programs
            double entropy = 0;
            for( final String v : values ) {
                final double p = (double) counts.get( v ) / covered;
                entropy -= p * StrictMath.log( p );
            }
            final double balance = entropy / StrictMath.log( distinct );
            candidates.add( new VisCandidate( id,
                                              s._ref,
                                              s._element_slot,
                                              s._element_slot ? displayName( s._ref )
                                                      : prettifyVisLabel( localName( s._ref ) ),
                                              numeric,
                                              covered,
                                              total,
                                              values,
                                              counts,
                                              ( !s._element_slot || numeric ) ? canon : null,
                                              s._cut,
                                              ( (double) covered / total ) * balance,
                                              color_mode,
                                              switchable,
                                              wide,
                                              near_unique,
                                              sparse,
                                              isDeprioritizedRef( s._ref ),
                                              distinct <= VIS_MAX_SHAPE_CATEGORIES ) );
        }
        // two refs whose prettified labels coincide (different namespaces) are BOTH shown as their full ref
        final Map<String, Integer> label_count = new HashMap<String, Integer>();
        for( final VisCandidate c : candidates ) {
            final Integer n = label_count.get( c._label );
            label_count.put( c._label, ( n == null ) ? 1 : ( n + 1 ) );
        }
        for( final VisCandidate c : candidates ) {
            if ( ( label_count.get( c._label ) > 1 ) && !c._element_slot ) {
                c._label = c._ref;
            }
        }
        Collections.sort( candidates, VIS_ORDER );
        return candidates;
    }

    private static void addVisValue( final Map<String, List<String>> per_node, final String id, final String value ) {
        final String v = jsClean( value );
        if ( v == null ) {
            return;
        }
        List<String> l = per_node.get( id );
        if ( l == null ) {
            l = new ArrayList<String>( 1 );
            per_node.put( id, l );
        }
        l.add( v );
    }

    private static String taxonomySlotValue( final Taxonomy t, final String slot ) {
        return TAX_CODE_REF.equals( slot ) ? t.getTaxonomyCode()
                : ( TAX_SCI_NAME_REF.equals( slot ) ? t.getScientificName() : t.getCommonName() );
    }

    private static String sequenceSlotValue( final Sequence q, final String slot ) {
        return SEQ_NAME_REF.equals( slot ) ? q.getName()
                : ( SEQ_SYMBOL_REF.equals( slot ) ? q.getSymbol() : q.getGeneName() );
    }

    /** Every external node of the tree, collapsed or not (a collapse is a view; candidacy reads the tree). Walks
     *  the topology rather than the cached external-node list, which an in-place edit can leave stale. */
    static List<PhylogenyNode> allExternalNodes( final Phylogeny tree ) {
        final List<PhylogenyNode> tips = new ArrayList<PhylogenyNode>();
        if ( ( tree == null ) || tree.isEmpty() ) {
            return tips;
        }
        final Deque<PhylogenyNode> stack = new ArrayDeque<PhylogenyNode>();
        stack.push( tree.getRoot() );
        while ( !stack.isEmpty() ) {
            final PhylogenyNode n = stack.pop();
            if ( n.isExternal() ) {
                tips.add( n );
            }
            else {
                for( final PhylogenyNode child : n.getDescendants() ) {
                    stack.push( child );
                }
            }
        }
        return tips;
    }

    /**
     * Reads a node's value for a candidate exactly as the classifier grouped it -- forester.js
     * {@code visualizationNodeValue} -- so a node can never carry a value that maps to no colour: a property is
     * folded through {@link #visDisplayLabel} and mapped to its group's representative (a numeric field's "1.0"
     * reads "1"); an element slot is verbatim, except that a numeric slot folds its spellings too. Null when the node
     * has none.
     * <p>
     * A value that folds to NOTHING ("_", a host that is only its ";" qualifier) is no value, and the scan moves on to
     * the next property with the ref, so "_" beside "Human" on one tip reads "Human" (forester.js 635741b).
     */
    static String visualizationNodeValue( final PhylogenyNode node, final VisCandidate c ) {
        final NodeData nd = node.getNodeData();
        if ( nd == null ) {
            return null;
        }
        if ( !c._element_slot ) {
            if ( nd.getProperties() != null ) {
                for( final Property p : nd.getProperties().getProperties() ) {
                    if ( c._ref.equals( p.getRef() ) && isNodeScopedProperty( p ) ) {
                        final String v = jsClean( p.getValue() );
                        if ( v != null ) {
                            if ( c._canon == null ) {
                                return v;
                            }
                            final String display = visDisplayLabel( v, c._cut );
                            if ( display.isEmpty() ) {
                                continue; // folds to nothing: no value, as the classifier counted it -- keep looking
                            }
                            final String r = c._canon.get( display.toLowerCase( Locale.ROOT ) );
                            return ( r != null ) ? r : display;
                        }
                    }
                }
            }
            return null;
        }
        if ( c._ref.startsWith( "tax:" ) ) {
            if ( nd.getTaxonomies() != null ) {
                for( final Taxonomy t : nd.getTaxonomies() ) {
                    final String v = jsClean( taxonomySlotValue( t, c._ref ) );
                    if ( v != null ) {
                        return foldSlotValue( c, v );
                    }
                }
            }
            return null;
        }
        if ( nd.getSequences() != null ) {
            for( final Sequence q : nd.getSequences() ) {
                final String v = jsClean( sequenceSlotValue( q, c._ref ) );
                if ( v != null ) {
                    return foldSlotValue( c, v );
                }
            }
        }
        return null;
    }

    private static String foldSlotValue( final VisCandidate c, final String v ) {
        if ( c._canon == null ) {
            return v;
        }
        final String r = c._canon.get( v );
        return ( r != null ) ? r : v;
    }

    /** The visualization a tree OPENS with: the first candidate that is not wide (21+ values are offered, never
     *  imposed -- and they do not block what sits below them), or null. forester.js {@code openingVisualization}. */
    static VisCandidate openingVisualization( final List<VisCandidate> candidates ) {
        if ( candidates != null ) {
            for( final VisCandidate c : candidates ) {
                if ( !c._wide ) {
                    return c;
                }
            }
        }
        return null;
    }

    /**
     * What a VIEW shows of a candidate: its values, counts and coverage over {@code tips} -- forester.js
     * {@code visualizationSummary}. Never re-decides candidacy. A numeric candidate also gets the view's colour-mode
     * band (up to 10 distinct: colours, switchable; 11-20: range, switchable; above: range, fixed); a category keeps
     * its mode.
     */
    static VisSummary visualizationSummary( final VisCandidate c, final List<PhylogenyNode> tips ) {
        final Map<String, Integer> counts = new LinkedHashMap<String, Integer>();
        int total = 0;
        int coverage = 0;
        for( final PhylogenyNode n : tips ) {
            total++;
            final String v = visualizationNodeValue( n, c );
            if ( v != null ) {
                coverage++;
                final Integer k = counts.get( v );
                counts.put( v, ( k == null ) ? 1 : ( k + 1 ) );
            }
        }
        final List<String> values = new ArrayList<String>( counts.keySet() );
        if ( c._numeric ) {
            Collections.sort( values, new Comparator<String>() {

                @Override
                public int compare( final String a, final String b ) {
                    return Double.compare( Double.parseDouble( a ), Double.parseDouble( b ) );
                }
            } );
            final int d = values.size();
            return new VisSummary( values, counts, coverage, total,
                                   ( d <= VIS_NUMERIC_CATEGORY_MAX ) ? "category" : "range",
                                   d <= VIS_MAX_COLOR_CATEGORIES );
        }
        Collections.sort( values );
        return new VisSummary( values, counts, coverage, total, null, false );
    }

    /**
     * The candidates of a tree the user has EDITED, keeping the fields they had chosen -- forester.js
     * {@code visualizationCandidatesKeeping}. The refusal rules decide what is offered, never what is already
     * chosen: colouring by Host and deleting every clade but one must not silently uncolour the tree because one
     * host is "not a category". A chosen field that is no longer offered but still carries a value on some tip is
     * appended after the offered ones and flagged {@code _kept}, so the menu holds it exactly as long as the user
     * does; the next edit drops it unless it is still chosen. {@code chosen} are the previous candidate objects --
     * their grouping travels with them. A chosen field with no value left anywhere is dropped.
     */
    static List<VisCandidate> visualizationCandidatesKeeping( final Phylogeny tree,
                                                              final java.util.Collection<VisCandidate> chosen ) {
        final List<VisCandidate> candidates = visualizationCandidates( tree );
        final Set<String> ids = new HashSet<String>();
        for( final VisCandidate c : candidates ) {
            ids.add( c._id );
        }
        if ( chosen != null ) {
            for( final VisCandidate c : chosen ) {
                if ( ( c == null ) || ids.contains( c._id ) ) {
                    continue;
                }
                final VisSummary s = visualizationSummary( c, allExternalNodes( tree ) );
                if ( s._coverage == 0 ) {
                    continue;
                }
                c._values = s._values;
                c._counts = s._counts;
                c._coverage = s._coverage;
                c._total = s._total;
                if ( c._numeric ) {
                    c._color_mode = s._color_mode;
                    c._switchable = s._switchable;
                }
                c._kept = true;
                ids.add( c._id );
                candidates.add( c );
            }
        }
        return candidates;
    }

    /** The candidate with this ref in {@code candidates}, or null. */
    static VisCandidate findCandidate( final List<VisCandidate> candidates, final String ref ) {
        if ( ( candidates != null ) && ( ref != null ) ) {
            for( final VisCandidate c : candidates ) {
                if ( ref.equals( c._ref ) ) {
                    return c;
                }
            }
        }
        return null;
    }

    /** The Color-by band of a CANDIDATE in a view: a numeric candidate's band follows its distinct count among the
     *  tips on screen; a categorical candidate is a category in every view, whatever its values happen to look like
     *  there. See {@link #colorModeBand(List, String)} for the legacy detection used by non-candidate refs. */
    static ModeBand colorModeBand( final VisCandidate c, final List<PhylogenyNode> tips ) {
        if ( !c._numeric ) {
            return ModeBand.NOT_NUMERIC;
        }
        final int d = visualizationSummary( c, tips )._distinct;
        if ( d <= VIS_NUMERIC_CATEGORY_MAX ) {
            return ModeBand.SMALL;
        }
        return ( d <= VIS_MAX_COLOR_CATEGORIES ) ? ModeBand.MEDIUM : ModeBand.LARGE;
    }

    /** The refs of {@link #visualizationCandidates(Phylogeny)}, best first -- what "Color by" offers for a tree. */
    static List<String> colorableRefs( final Phylogeny phylogeny ) {
        final List<String> refs = new ArrayList<String>();
        for( final VisCandidate c : visualizationCandidates( phylogeny ) ) {
            refs.add( c._ref );
        }
        return refs;
    }

    /** The field a newly opened tree is coloured by ({@link #openingVisualization}), or null for none. */
    static String autoColorCandidate( final Phylogeny phylogeny ) {
        final VisCandidate c = openingVisualization( visualizationCandidates( phylogeny ) );
        return ( c == null ) ? null : c._ref;
    }

    /** The NUMERIC candidates of a tree (every value a number under the pinned grammar), best first -- the fields
     *  "Size by" can scale a symbol by. */
    static List<String> numericRefs( final Phylogeny phylogeny ) {
        final List<String> refs = new ArrayList<String>();
        for( final VisCandidate c : visualizationCandidates( phylogeny ) ) {
            if ( c._numeric ) {
                refs.add( c._ref );
            }
        }
        return refs;
    }

    /** The numeric refs among {@code colorable}, in its order. */
    static List<String> numericRefs( final Phylogeny phylogeny, final List<String> colorable ) {
        final Set<String> numeric = new HashSet<String>( numericRefs( phylogeny ) );
        final List<String> refs = new ArrayList<String>();
        for( final String ref : colorable ) {
            if ( numeric.contains( ref ) ) {
                refs.add( ref );
            }
        }
        return refs;
    }
}
