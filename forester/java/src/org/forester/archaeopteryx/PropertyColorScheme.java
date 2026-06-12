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
 * are not merged. A few refs are special-cased: {@code year} is colored by a continuous
 * gradient over its numeric range, while {@code country} and {@code host} first drop a
 * trailing qualifier -- everything from the first {@code :} (country, so {@code USA:CA} and
 * {@code USA:IL} share a color) or {@code ;} (host, so {@code Homo sapiens; male 35} and
 * {@code Homo sapiens; female old} share a color) onward -- before that grouping.
 * <p>
 * This is independent of the explicit {@code style:*} visualization properties and
 * the "Visual Styles" feature.
 */
final class PropertyColorScheme {

    // A categorical palette of reasonably distinct colors that read on a white or a
    // dark canvas; cycled when a property has more distinct values than entries.
    private static final Color[] PALETTE = { new Color( 0xE6194B ), new Color( 0x3CB44B ), new Color( 0x4363D8 ),
            new Color( 0xF58231 ), new Color( 0x911EB4 ), new Color( 0x469990 ), new Color( 0xF032E6 ),
            new Color( 0x9A6324 ), new Color( 0x800000 ), new Color( 0x000075 ), new Color( 0x808000 ),
            new Color( 0xE67E22 ), new Color( 0x2980B9 ), new Color( 0x16A085 ), new Color( 0xC0392B ),
            new Color( 0x8E44AD ), new Color( 0xD35400 ), new Color( 0x27AE60 ), new Color( 0x7F8C8D ),
            new Color( 0xB7950B ), new Color( 0x1F618D ), new Color( 0x6C3483 ), new Color( 0xBD1E51 ),
            new Color( 0x117864 ) };

    private final String             _ref;
    // Categorical mode: one palette color per distinct value. _value_to_color maps the
    // (display) representative of each group to its color, for the legend; _key_to_color
    // maps the normalized grouping key to the same color, for looking up a node's color.
    private final Map<String, Color> _value_to_color;
    private final Map<String, Color> _key_to_color;
    // representative label -> number of (visible) leaves in that group, for the legend
    private final Map<String, Integer> _value_to_count;
    // Continuous mode (numeric properties such as "year"): a blue->red gradient spanning
    // [_min, _max] instead of distinct colors. _gradient is false for categorical refs.
    private final boolean            _gradient;
    private final double             _min;
    private final double             _max;
    // For refs whose value carries a trailing qualifier, the delimiter at which the value is
    // truncated before grouping: ':' for "country" (USA:CA == USA:IL), ';' for "host"
    // (Homo sapiens; male 35 == Homo sapiens; female old). 0 means keep the whole value.
    private final char               _truncate_at;

    PropertyColorScheme( final Phylogeny phylogeny, final String ref ) {
        _ref = ref;
        _gradient = isYearRef( ref );
        _truncate_at = truncationDelimiter( ref );
        _value_to_color = new LinkedHashMap<String, Color>();
        _key_to_color = new LinkedHashMap<String, Color>();
        _value_to_count = new LinkedHashMap<String, Integer>();
        // Color from the leaves actually on screen (those hidden under a collapsed node are
        // excluded), so the colors and legend track the displayed (sub)tree as the user
        // navigates into subtrees, collapses clades, or deletes nodes.
        final List<PhylogenyNode> leaves = visibleExternalNodes( phylogeny );
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
                final Color color = PALETTE[ i++ % PALETTE.length ];
                _value_to_color.put( g[ 0 ], color ); // _value_to_color is now ordered most-frequent first
                _key_to_color.put( g[ 1 ], color );
                _value_to_count.put( g[ 0 ], key_to_total.get( g[ 1 ] ) );
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

    /** Representative-label to color map of all distinct values, ordered most-frequent first. */
    Map<String, Color> getValueColors() {
        return _value_to_color;
    }

    /** Representative-label to (visible) leaf count for each value; empty in gradient mode. */
    Map<String, Integer> getValueCounts() {
        return _value_to_count;
    }

    /**
     * The legend entries: the {@code max} most frequent values (so the legend describes what is most
     * visible on the tree), re-sorted alphabetically for readability, mapped to their colors. Returns
     * fewer than {@code max} entries when there are fewer distinct values; empty in gradient mode.
     */
    Map<String, Color> legendValues( final int max ) {
        final List<Map.Entry<String, Color>> top = new ArrayList<Map.Entry<String, Color>>();
        int i = 0;
        for( final Map.Entry<String, Color> e : _value_to_color.entrySet() ) { // already most-frequent first
            if ( i++ >= max ) {
                break;
            }
            top.add( e );
        }
        Collections.sort( top, new Comparator<Map.Entry<String, Color>>() {

            @Override
            public int compare( final Map.Entry<String, Color> a, final Map.Entry<String, Color> b ) {
                return String.CASE_INSENSITIVE_ORDER.compare( a.getKey(), b.getKey() );
            }
        } );
        final Map<String, Color> result = new LinkedHashMap<String, Color>();
        for( final Map.Entry<String, Color> e : top ) {
            result.put( e.getKey(), e.getValue() );
        }
        return result;
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
            }
        }
        s = s.trim().replace( '_', ' ' );
        s = s.replaceAll( "\\s+", " " );
        return canonicalSynonym( s );
    }

    /**
     * Folds a tiny set of unambiguous synonyms to a canonical name so they share a color/legend
     * entry. Deliberately minimal -- this is not general metadata normalization; currently only the
     * common-name "human"/"humans" is folded into the scientific name "Homo sapiens".
     */
    private static String canonicalSynonym( final String label ) {
        final String lower = label.toLowerCase( Locale.ROOT );
        if ( lower.equals( "human" ) || lower.equals( "humans" ) ) {
            return "Homo sapiens";
        }
        return label;
    }

    /** The normalized key a value is grouped/colored by: its display label, case-folded. */
    private String groupKey( final String v ) {
        return displayLabel( v ).toLowerCase( Locale.ROOT );
    }

    /** The most frequent spelling in a group (ties broken alphabetically). */
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
        return best;
    }

    /** Color at fraction {@code t} (0..1, low value to high value) of the gradient. */
    Color gradientColorAt( final double t ) {
        final double tt = ( t < 0.0 ) ? 0.0 : ( ( t > 1.0 ) ? 1.0 : t );
        // hue sweep blue (low) -> red (high), at a fixed saturation/brightness so every
        // step stays legible on both the white and the dark canvas.
        final float hue = (float) ( 0.66 * ( 1.0 - tt ) );
        return Color.getHSBColor( hue, 0.78f, 0.92f );
    }

    String getGradientMinLabel() {
        return formatNumber( _min );
    }

    String getGradientMaxLabel() {
        return formatNumber( _max );
    }

    private static Double parseNumber( final String s ) {
        if ( ForesterUtil.isEmpty( s ) ) {
            return null;
        }
        try {
            return Double.valueOf( s.trim() );
        }
        catch ( final NumberFormatException e ) {
            return null;
        }
    }

    private static String formatNumber( final double d ) {
        return ( d == Math.rint( d ) ) ? Long.toString( (long) d ) : Double.toString( d );
    }

    /** True for the {@code year} property (in any namespace), which is colored as a gradient. */
    private static boolean isYearRef( final String ref ) {
        return refNameEquals( ref, "year" );
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
    private static List<PhylogenyNode> visibleExternalNodes( final Phylogeny phylogeny ) {
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

    private static String valueFor( final PhylogenyNode node, final String ref ) {
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
        for( final PhylogenyNode node : phylogeny.getExternalNodes() ) {
            if ( ( node.getNodeData() != null ) && ( node.getNodeData().getProperties() != null ) ) {
                for( final Property p : node.getNodeData().getProperties().getProperties() ) {
                    if ( p.getRef().startsWith( NodeVisualData.APTX_VISUALIZATION_REF ) ) {
                        continue;
                    }
                    Set<String> vs = ref_to_values.get( p.getRef() );
                    if ( vs == null ) {
                        vs = new HashSet<String>();
                        ref_to_values.put( p.getRef(), vs );
                    }
                    vs.add( p.getValue() );
                }
            }
        }
        final int leaves = phylogeny.getExternalNodes().size();
        final List<String> refs = new ArrayList<String>();
        for( final Map.Entry<String, Set<String>> e : ref_to_values.entrySet() ) {
            final int distinct = e.getValue().size();
            if ( ( distinct >= 2 ) && ( distinct < leaves ) ) {
                refs.add( e.getKey() );
            }
        }
        return refs;
    }
}
