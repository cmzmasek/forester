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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.util.ForesterUtil;

/**
 * The tip-aligned annotation columns drawn to the right of the tree: an ordered set of columns, each backed
 * by one node property (ref) and rendered as a color strip, a heat-map cell, a bar, or the value's text.
 * <p>
 * The render type follows the data: a <em>categorical</em> field is a color strip (or text), a
 * <em>numeric</em> field a heat map / bar (or text). Colors and the categorical-vs-numeric decision reuse
 * {@link PropertyColorScheme}, so a column's colors match what "Color by: {field}" would produce and share
 * the same palettes, value grouping, and gradient.
 * <p>
 * Pure data model: it reads the tree and answers per-cell questions (color / bar fraction / text); the
 * rendering and the width/legend layout live in {@code TreePanel}.
 */
final class AnnotationColumns {

    enum Type {
        /** A filled cell colored by the value's categorical color. */
        COLOR_STRIP,
        /** A centered shape glyph colored by the value's categorical color (the same palette a COLOR_STRIP uses):
         *  FILLED when the value is present, a hollow OUTLINE when the value reads as explicitly false/absent
         *  (see {@link #isFalsy}), and nothing at all when the value is missing. One clean tip-aligned mark column
         *  covering the presence/absence (binary), symbol, and external-shape annotation cases. */
        SYMBOL,
        /** A filled cell colored by the value's position in the numeric gradient. */
        HEATMAP,
        /** Like HEATMAP, but all MATRIX columns share ONE color scale (one min/max, one legend) and draw as a
         *  contiguous grid -- the clustergram heat-map matrix (e.g. abundance across samples). */
        MATRIX,
        /** A horizontal bar whose length is the value's fraction of the numeric range. */
        BAR,
        /** The raw value drawn as text. */
        TEXT
    }

    /** How a {@link Type#SYMBOL} cell renders for a given tip. */
    enum Fill {
        /** A solid glyph -- the value is present and not explicitly false. */
        FILLED,
        /** A hollow (outlined) glyph -- the value reads as explicitly false/absent. */
        OUTLINE,
        /** Nothing is drawn -- the value is missing (or the column is not a SYMBOL column). */
        NONE
    }

    /** Values (case-insensitive, trimmed) that a SYMBOL column renders as a hollow OUTLINE rather than a
     *  filled glyph -- the "explicitly false/absent" tokens (a present-but-negative binary state). Anything
     *  else present is FILLED; a missing value draws NONE. */
    private static final Set<String> FALSY = new HashSet<String>( Arrays.asList( "0", "false", "no", "n",
            "absent", "-" ) );

    /** True iff {@code value} reads as explicitly false/absent (see {@link #FALSY}); a hollow SYMBOL glyph. */
    static boolean isFalsy( final String value ) {
        return ( value != null ) && FALSY.contains( value.trim().toLowerCase( Locale.ROOT ) );
    }

    /** A short human-readable name for a render type (for the picker). */
    static String label( final Type type ) {
        switch ( type ) {
            case COLOR_STRIP:
                return "Color strip";
            case SYMBOL:
                return "Symbol";
            case HEATMAP:
                return "Heat map";
            case MATRIX:
                return "Heat map (matrix)";
            case BAR:
                return "Bar";
            case TEXT:
                return "Text";
            default:
                return type.name();
        }
    }

    /** A requested column: which field, rendered how. */
    static final class ColumnSpec {

        final String _ref;
        final Type   _type;

        ColumnSpec( final String ref, final Type type ) {
            _ref = ref;
            _type = type;
        }
    }

    /** A resolved column: its field, type, header, and (for non-text columns) its color scheme. */
    static final class Column {

        private final String              _ref;
        private final Type                _type;
        private final PropertyColorScheme _scheme; // null for TEXT

        Column( final String ref, final Type type, final PropertyColorScheme scheme ) {
            _ref = ref;
            _type = type;
            _scheme = scheme;
        }

        Type getType() {
            return _type;
        }

        /** The human-readable header (the field name, namespace stripped and prettified). */
        String getHeader() {
            return PropertyColorScheme.displayName( _ref );
        }

        PropertyColorScheme getScheme() {
            return _scheme;
        }
    }

    private final List<Column> _columns;

    AnnotationColumns( final Phylogeny phylogeny, final List<ColumnSpec> specs ) {
        _columns = new ArrayList<Column>();
        for( final ColumnSpec spec : specs ) {
            final PropertyColorScheme scheme = ( spec._type == Type.TEXT ) ? null
                    : new PropertyColorScheme( phylogeny, spec._ref );
            _columns.add( new Column( spec._ref, spec._type, scheme ) );
        }
        // heat-map MATRIX: one shared color scale across ALL matrix columns, so a cell's color is comparable across
        // samples. Compute the range over every matrix ref's values and stamp it onto each matrix column's scheme.
        // Deliberately over getExternalNodes() (ALL tips, incl. collapsed) rather than the visible tips a per-column
        // HEATMAP uses: a shared matrix scale should stay STABLE as the user collapses clades (its whole point is
        // cross-sample comparability), not rescale and change what a color means on every collapse.
        double mmin = Double.POSITIVE_INFINITY, mmax = Double.NEGATIVE_INFINITY;
        for( final Column c : _columns ) {
            if ( c._type == Type.MATRIX ) {
                for( final PhylogenyNode leaf : phylogeny.getExternalNodes() ) {
                    final Double d = PropertyColorScheme.parseNumber( PropertyColorScheme.valueFor( leaf, c._ref ) );
                    if ( d != null ) {
                        mmin = Math.min( mmin, d );
                        mmax = Math.max( mmax, d );
                    }
                }
            }
        }
        if ( mmin <= mmax ) {
            for( final Column c : _columns ) {
                if ( ( c._type == Type.MATRIX ) && ( c._scheme != null ) ) {
                    c._scheme.setSharedRange( mmin, mmax );
                }
            }
        }
    }

    int size() {
        return _columns.size();
    }

    List<Column> getColumns() {
        return _columns;
    }

    Column getColumn( final int i ) {
        return _columns.get( i );
    }

    /**
     * The fill color for {@code node}'s cell in column {@code i} for a color-strip or heat-map column, or
     * {@code null} when the node has no value there, or the column is a bar / text column.
     */
    Color cellColor( final PhylogenyNode node, final int i ) {
        final Column c = _columns.get( i );
        if ( ( c._scheme == null ) || ( c._type == Type.BAR ) || ( c._type == Type.TEXT ) ) {
            return null;
        }
        return c._scheme.colorFor( node );
    }

    /**
     * For a SYMBOL column, whether {@code node}'s cell draws a FILLED glyph, a hollow OUTLINE, or NONE: NONE when
     * the node has no value there, OUTLINE when the value reads as explicitly false/absent (see {@link #isFalsy}),
     * else FILLED. Always NONE for a non-SYMBOL column. The glyph color comes from {@link #cellColor}.
     */
    Fill symbolFill( final PhylogenyNode node, final int i ) {
        final Column c = _columns.get( i );
        if ( c._type != Type.SYMBOL ) {
            return Fill.NONE;
        }
        final String v = valueOrEmpty( node, c._ref );
        if ( ForesterUtil.isEmpty( v ) ) {
            return Fill.NONE;
        }
        return isFalsy( v ) ? Fill.OUTLINE : Fill.FILLED;
    }

    /**
     * For a bar column, {@code node}'s value as a fraction of the numeric range in {@code [0, 1]}; {@code NaN}
     * when the column is not a bar or the node has no numeric value (draw nothing).
     */
    double barFraction( final PhylogenyNode node, final int i ) {
        final Column c = _columns.get( i );
        if ( ( c._type != Type.BAR ) || ( c._scheme == null ) ) {
            return Double.NaN;
        }
        final Double f = c._scheme.gradientFraction( node );
        return ( f == null ) ? Double.NaN : f.doubleValue();
    }

    /** The raw value drawn in {@code node}'s cell for a text column, or {@code ""} when it has none. */
    String cellText( final PhylogenyNode node, final int i ) {
        return valueOrEmpty( node, _columns.get( i )._ref );
    }

    private static String valueOrEmpty( final PhylogenyNode node, final String ref ) {
        if ( ( node.getNodeData() == null ) || ( node.getNodeData().getProperties() == null ) ) {
            return "";
        }
        final List<Property> ps = node.getNodeData().getProperties().getProperties( ref );
        if ( ps.isEmpty() || ForesterUtil.isEmpty( ps.get( 0 ).getValue() ) ) {
            return "";
        }
        return ps.get( 0 ).getValue();
    }

    /** The default render type for a field: heat map for a numeric field, color strip for a categorical one. */
    static Type defaultType( final Phylogeny phylogeny, final String ref ) {
        return new PropertyColorScheme( phylogeny, ref ).isGradient() ? Type.HEATMAP : Type.COLOR_STRIP;
    }

    /**
     * The render types offered for a field, given its data: a numeric field can be a heat map, a bar, or text;
     * a categorical field a color strip or text.
     */
    static List<Type> allowedTypes( final Phylogeny phylogeny, final String ref ) {
        final List<Type> types = new ArrayList<Type>();
        if ( new PropertyColorScheme( phylogeny, ref ).isGradient() ) {
            types.add( Type.HEATMAP );
            types.add( Type.MATRIX );
            types.add( Type.BAR );
            types.add( Type.TEXT );
        }
        else {
            types.add( Type.COLOR_STRIP );
            types.add( Type.SYMBOL );
            types.add( Type.TEXT );
        }
        return types;
    }
}
