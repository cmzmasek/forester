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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * "Size by property": maps a numeric phyloXML property value at a (tip) node to a symbol diameter, so a node's dot
 * grows with the value -- the size counterpart of {@link PropertyColorScheme}. Combine the two (color = property A,
 * size = property B) for a two-dimensional figure. The value maps to the disc AREA, not the diameter, so the
 * perceived quantity is honest (a linear diameter map over-weights large values).
 *
 * <p>Built over the VISIBLE external nodes (those hidden under a collapsed clade are excluded), so the range tracks
 * the displayed (sub)tree as the user navigates -- exactly like {@link PropertyColorScheme}'s gradient. The
 * numeric-value extraction, the visible-leaf walk and the numeric-ref detection are shared with that class.
 */
final class PropertySizeScale {

    // The symbol diameter spans [MIN_FACTOR, MAX_FACTOR] x the caller's base size; MIN_FACTOR = 1 so the smallest
    // (and any value-less) node keeps the base dot, and the largest grows to MAX_FACTOR x.
    private static final double MIN_FACTOR = 1.0;
    private static final double MAX_FACTOR = 3.0;

    private final String _ref;
    private final double _min;
    private final double _max;

    PropertySizeScale( final Phylogeny phylogeny, final String ref ) {
        _ref = ref;
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for( final PhylogenyNode leaf : PropertyColorScheme.visibleExternalNodes( phylogeny ) ) {
            final Double d = PropertyColorScheme.parseNumber( PropertyColorScheme.valueFor( leaf, ref ) );
            if ( d != null ) {
                min = Math.min( min, d );
                max = Math.max( max, d );
            }
        }
        _min = min;
        _max = max;
    }

    String getRef() {
        return _ref;
    }

    /** True when no visible leaf carries a numeric value for this ref -- nothing to size by. */
    boolean isEmpty() {
        return _min > _max;
    }

    /** True iff {@code node} carries a numeric value for this ref, so it gets a real, value-driven diameter.
     *  A value-less tip has no size to show -- the caller draws no size dot for it (see {@link #diameterFor}),
     *  keeping "no data" visually distinct from the minimum value and matching Color-by's draw-nothing behavior. */
    boolean hasValueFor( final PhylogenyNode node ) {
        return PropertyColorScheme.parseNumber( PropertyColorScheme.valueFor( node, _ref ) ) != null;
    }

    /** The symbol diameter (px) for {@code node}: area-proportional within [MIN_FACTOR,MAX_FACTOR]x{@code base_size}.
     *  A node with no numeric value keeps {@code base_size} (the smallest / unscaled dot). */
    float diameterFor( final PhylogenyNode node, final float base_size ) {
        final Double d = PropertyColorScheme.parseNumber( PropertyColorScheme.valueFor( node, _ref ) );
        if ( d == null ) {
            return base_size;
        }
        return diameter( normalized( d ), base_size );
    }

    /** The node's value mapped to [0,1]; an all-equal range (no spread) maps to 0 (all at the base size). */
    double normalized( final double value ) {
        return ( _max > _min ) ? ( ( value - _min ) / ( _max - _min ) ) : 0.0;
    }

    /**
     * Area-proportional diameter: {@code t} in [0,1] maps LINEARLY to the disc area between (MIN_FACTOR*base)^2 and
     * (MAX_FACTOR*base)^2, so the returned diameter is the square root of the interpolated area. {@code t} is
     * clamped to [0,1]. Pure -- the single source of the size mapping.
     */
    static float diameter( final double t, final float base_size ) {
        final double min_d = MIN_FACTOR * base_size;
        final double max_d = MAX_FACTOR * base_size;
        final double clamped = ( t < 0 ) ? 0 : ( ( t > 1 ) ? 1 : t );
        final double area = ( min_d * min_d ) + ( clamped * ( ( max_d * max_d ) - ( min_d * min_d ) ) );
        return (float) Math.sqrt( area );
    }
}
