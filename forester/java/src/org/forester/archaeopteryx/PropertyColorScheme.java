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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.Property;
import org.forester.util.ForesterUtil;

/**
 * Assigns a distinct color to each distinct value of a chosen phyloXML property
 * (e.g. {@code repseq:host}), so leaves can be colored on the fly by that property:
 * the same value always maps to the same color. The categorical palette is cycled
 * when a property has more distinct values than there are palette entries.
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
    // Categorical mode: one palette color per distinct value.
    private final Map<String, Color> _value_to_color;
    // Continuous mode (numeric properties such as "year"): a blue->red gradient spanning
    // [_min, _max] instead of distinct colors. _gradient is false for categorical refs.
    private final boolean            _gradient;
    private final double             _min;
    private final double             _max;

    PropertyColorScheme( final Phylogeny phylogeny, final String ref ) {
        _ref = ref;
        _gradient = isYearRef( ref );
        _value_to_color = new LinkedHashMap<String, Color>();
        if ( _gradient ) {
            double min = Double.POSITIVE_INFINITY;
            double max = Double.NEGATIVE_INFINITY;
            if ( ( phylogeny != null ) && !phylogeny.isEmpty() ) {
                for( final PhylogenyNode node : phylogeny.getExternalNodes() ) {
                    final Double d = parseNumber( valueFor( node, ref ) );
                    if ( d != null ) {
                        min = Math.min( min, d );
                        max = Math.max( max, d );
                    }
                }
            }
            _min = min;
            _max = max;
        }
        else {
            _min = 0;
            _max = 0;
            final TreeSet<String> values = new TreeSet<String>();
            if ( ( phylogeny != null ) && !phylogeny.isEmpty() ) {
                for( final PhylogenyNode node : phylogeny.getExternalNodes() ) {
                    final String v = valueFor( node, ref );
                    if ( !ForesterUtil.isEmpty( v ) ) {
                        values.add( v );
                    }
                }
            }
            int i = 0;
            for( final String v : values ) {
                _value_to_color.put( v, PALETTE[ i++ % PALETTE.length ] );
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
        return ( v == null ) ? null : _value_to_color.get( v );
    }

    /** Ordered (alphabetical) value to color map, for building a (categorical) legend. */
    Map<String, Color> getValueColors() {
        return _value_to_color;
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
        if ( ForesterUtil.isEmpty( ref ) ) {
            return false;
        }
        final int colon = ref.lastIndexOf( ':' );
        final String name = ( colon >= 0 ) ? ref.substring( colon + 1 ) : ref;
        return name.equalsIgnoreCase( "year" );
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
