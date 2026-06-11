// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// Contact: phylosoft @ gmail . com

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
 * trimming, collapsing whitespace, treating {@code _} as a space, and case-folding -- so
 * e.g. {@code Human}/{@code human}/{@code homo_sapiens} get one color (it cannot, however,
 * merge semantically-equal but lexically-different values such as {@code man} vs
 * {@code H. sapiens}). Two refs are special-cased: {@code year} is colored by a continuous
 * gradient over its numeric range, and {@code country} is grouped by the part before the
 * first {@code :} (so {@code USA:CA} and {@code USA:IL} share a color).
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
    // Continuous mode (numeric properties such as "year"): a blue->red gradient spanning
    // [_min, _max] instead of distinct colors. _gradient is false for categorical refs.
    private final boolean            _gradient;
    private final double             _min;
    private final double             _max;
    // "country": values are grouped by the part before the first ':' (USA:CA == USA:IL).
    private final boolean            _country;

    PropertyColorScheme( final Phylogeny phylogeny, final String ref ) {
        _ref = ref;
        _gradient = isYearRef( ref );
        _country = isCountryRef( ref );
        _value_to_color = new LinkedHashMap<String, Color>();
        _key_to_color = new LinkedHashMap<String, Color>();
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
            // "country", the subdivision after ':') so e.g. "Human"/"human"/"homo_sapiens "
            // share one color. Each group's legend label is its most frequent spelling.
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
            // [ representative label, key ] per group, ordered by label (case-insensitive)
            final List<String[]> groups = new ArrayList<String[]>();
            for( final Map.Entry<String, Map<String, Integer>> e : key_to_label_counts.entrySet() ) {
                groups.add( new String[] { representative( e.getValue() ), e.getKey() } );
            }
            Collections.sort( groups, new Comparator<String[]>() {

                @Override
                public int compare( final String[] a, final String[] b ) {
                    return String.CASE_INSENSITIVE_ORDER.compare( a[ 0 ], b[ 0 ] );
                }
            } );
            int i = 0;
            for( final String[] g : groups ) {
                final Color color = PALETTE[ i++ % PALETTE.length ];
                _value_to_color.put( g[ 0 ], color );
                _key_to_color.put( g[ 1 ], color );
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

    /** Ordered (alphabetical) representative-label to color map, for building a (categorical) legend. */
    Map<String, Color> getValueColors() {
        return _value_to_color;
    }

    /**
     * The display label for a raw property value: trimmed, underscores as spaces, internal
     * whitespace collapsed; for "country" only the part before the first ':' (so USA:CA and
     * USA:IL both read as "USA"). Case is preserved -- this is what the legend shows.
     */
    private String displayLabel( final String v ) {
        String s = v;
        if ( _country ) {
            final int colon = s.indexOf( ':' );
            if ( colon >= 0 ) {
                s = s.substring( 0, colon );
            }
        }
        s = s.trim().replace( '_', ' ' );
        return s.replaceAll( "\\s+", " " );
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

    /** True for the {@code country} property (in any namespace), which drops the ':' subdivision. */
    private static boolean isCountryRef( final String ref ) {
        return refNameEquals( ref, "country" );
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
