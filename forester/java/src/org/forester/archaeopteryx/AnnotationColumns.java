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
        /** A composite SEGMENTED bar: several numeric fields (each set to STACKED_BAR) MERGE into one bar per tip,
         *  each field a distinctly-coloured series stacked end-to-end. <em>Absolute</em> by default -- a segment's
         *  length is its value as a fraction of the largest per-tip series total, so the bar's overall length shows
         *  magnitude AND its segments show composition; a <em>normalized</em> column instead fills the whole width and
         *  each segment shows that series' proportion of the tip's own total (composition only). For compositional /
         *  proportional data (iTOL's DATASET_MULTIBAR). */
        STACKED_BAR,
        /** A PIE glyph: several numeric fields (each set to PIE) MERGE into one pie per tip, each field a
         *  distinctly-coloured wedge whose angle is its share of the tip's total. The pie-chart alternative to a
         *  <em>normalized</em> stacked bar -- the same proportional composition, mapped to wedge angles instead of
         *  segment lengths (iTOL's DATASET_PIECHART). A pie is inherently proportional, so it has no absolute mode. */
        PIE,
        /** The raw value drawn as text. */
        TEXT
    }

    /** Whether {@code type} is a MERGED multi-series type -- several numeric fields drawn as one glyph per tip (a
     *  segmented {@link Type#STACKED_BAR} or a {@link Type#PIE}). */
    static boolean isMergedType( final Type type ) {
        return ( type == Type.STACKED_BAR ) || ( type == Type.PIE );
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

    /** The glyph a SYMBOL column draws -- chosen per column by the user (default {@link #CIRCLE}). Color still
     *  encodes the value and fill still encodes presence; the shape just distinguishes one symbol column from
     *  another (e.g. gene A = circle, gene B = square). */
    enum SymbolShape {
        CIRCLE, SQUARE, DIAMOND, TRIANGLE
    }

    /** A short human-readable name for a symbol shape (for the picker). */
    static String shapeLabel( final SymbolShape shape ) {
        switch ( shape ) {
            case CIRCLE:
                return "Circle";
            case SQUARE:
                return "Square";
            case DIAMOND:
                return "Diamond";
            case TRIANGLE:
                return "Triangle";
            default:
                return shape.name();
        }
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
            case STACKED_BAR:
                return "Stacked bar";
            case PIE:
                return "Pie";
            case TEXT:
                return "Text";
            default:
                return type.name();
        }
    }

    /** A requested column: which field, rendered how (for a SYMBOL column which glyph shape; for a STACKED_BAR
     *  field whether the merged bar is normalized). */
    static final class ColumnSpec {

        final String      _ref;
        final Type        _type;
        final SymbolShape _shape;      // the glyph for a SYMBOL column; CIRCLE (ignored) for every other type
        final boolean     _normalized; // STACKED_BAR only: normalize the merged bar to 100% (else absolute); false otherwise

        ColumnSpec( final String ref, final Type type ) {
            this( ref, type, SymbolShape.CIRCLE, false );
        }

        ColumnSpec( final String ref, final Type type, final SymbolShape shape ) {
            this( ref, type, shape, false );
        }

        /** For a STACKED_BAR field: the field plus the merged bar's normalize flag (glyph shape is irrelevant). */
        ColumnSpec( final String ref, final Type type, final boolean normalized ) {
            this( ref, type, SymbolShape.CIRCLE, normalized );
        }

        ColumnSpec( final String ref, final Type type, final SymbolShape shape, final boolean normalized ) {
            _ref = ref;
            _type = type;
            _shape = ( shape == null ) ? SymbolShape.CIRCLE : shape;
            _normalized = normalized;
        }
    }

    /** A resolved column: its field, type, glyph shape, header, and (for non-text columns) its color scheme. A
     *  STACKED_BAR column instead carries an ordered list of series fields with distinct colours + a shared scale. */
    static final class Column {

        private final String              _ref;
        private final Type                _type;
        private final SymbolShape         _shape;
        private final PropertyColorScheme _scheme;           // null for TEXT and a merged (STACKED_BAR / PIE) column
        private final List<String>        _stack_refs;        // merged only (else null): the series fields, in order
        private final List<Color>         _stack_colors;      // merged only: one distinct colour per series (parallel)
        private final double              _stack_max;         // merged only: the largest per-tip series total (STACKED_BAR absolute scale)
        private final boolean             _stack_normalized;  // merged only: divide each tip's series by its OWN total (always true for a PIE)

        /** A normal (single-field) column. */
        Column( final String ref, final Type type, final SymbolShape shape, final PropertyColorScheme scheme ) {
            this( ref, type, shape, scheme, null, null, 0, false );
        }

        /** A merged multi-series column ({@link Type#STACKED_BAR} or {@link Type#PIE}): several numeric fields drawn
         *  as one glyph per tip. */
        Column( final Type type, final List<String> stack_refs, final List<Color> stack_colors, final double stack_max,
                final boolean normalized ) {
            this( stack_refs.get( 0 ), type, SymbolShape.CIRCLE, null, stack_refs, stack_colors, stack_max, normalized );
        }

        private Column( final String ref, final Type type, final SymbolShape shape, final PropertyColorScheme scheme,
                        final List<String> stack_refs, final List<Color> stack_colors, final double stack_max,
                        final boolean normalized ) {
            _ref = ref;
            _type = type;
            _shape = shape;
            _scheme = scheme;
            _stack_refs = stack_refs;
            _stack_colors = stack_colors;
            _stack_max = stack_max;
            _stack_normalized = normalized;
        }

        Type getType() {
            return _type;
        }

        /** The glyph shape for a SYMBOL column (CIRCLE for a non-SYMBOL column). */
        SymbolShape getSymbolShape() {
            return _shape;
        }

        /** The human-readable header: the field name (namespace stripped, prettified), or a generic label for a
         *  merged column (which has several fields -- their names are in the legend, not one header). */
        String getHeader() {
            return isMergedType( _type ) ? label( _type ) : PropertyColorScheme.displayName( _ref );
        }

        PropertyColorScheme getScheme() {
            return _scheme;
        }
    }

    private final List<Column> _columns;

    AnnotationColumns( final Phylogeny phylogeny, final List<ColumnSpec> specs ) {
        _columns = new ArrayList<Column>();
        // every STACKED_BAR field MERGES into one segmented-bar column, and every PIE field into one pie column (each
        // shares its series set; a stacked bar also shares one normalize flag). Collect each group's refs first so the
        // merged column can be built with the whole series set, then emit it at the FIRST field of its type.
        final List<String> stack_refs = new ArrayList<String>();
        boolean stack_normalized = false;
        final List<String> pie_refs = new ArrayList<String>();
        for( final ColumnSpec spec : specs ) {
            if ( spec._type == Type.STACKED_BAR ) {
                if ( stack_refs.isEmpty() ) {
                    stack_normalized = spec._normalized; // the dialog sets every stacked field's flag alike; take the first
                }
                stack_refs.add( spec._ref );
            }
            else if ( spec._type == Type.PIE ) {
                pie_refs.add( spec._ref );
            }
        }
        boolean stack_added = false, pie_added = false;
        for( final ColumnSpec spec : specs ) {
            if ( spec._type == Type.STACKED_BAR ) {
                if ( !stack_added && !stack_refs.isEmpty() ) {
                    // a stacked bar is absolute unless normalized (its own toggle)
                    _columns.add( buildMergedColumn( phylogeny, Type.STACKED_BAR, stack_refs, stack_normalized ) );
                    stack_added = true; // the merged column sits at the FIRST stacked field's position; skip the rest
                }
                continue;
            }
            if ( spec._type == Type.PIE ) {
                if ( !pie_added && !pie_refs.isEmpty() ) {
                    // a pie is inherently proportional -> always "normalized" (each wedge = value / the tip's total)
                    _columns.add( buildMergedColumn( phylogeny, Type.PIE, pie_refs, true ) );
                    pie_added = true;
                }
                continue;
            }
            final PropertyColorScheme scheme = ( spec._type == Type.TEXT ) ? null
                    : new PropertyColorScheme( phylogeny, spec._ref );
            _columns.add( new Column( spec._ref, spec._type, spec._shape, scheme ) );
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

    /** Builds a merged multi-series column ({@code type} = STACKED_BAR or PIE) from its series {@code refs}: assigns
     *  each series a distinct colour and computes the absolute scale = the largest per-tip sum of (non-negative) series
     *  values (used by an absolute stacked bar; a normalized bar and every pie divide by the tip's own total instead).
     *  Deliberately over getExternalNodes() (like the MATRIX scale) so the absolute scale is stable across a collapse;
     *  it still rescales when navigating into a subtree, because the model is rebuilt with that subtree's phylogeny. */
    private static Column buildMergedColumn( final Phylogeny phylogeny, final Type type, final List<String> refs,
                                             final boolean normalized ) {
        double max_total = 0;
        for( final PhylogenyNode leaf : phylogeny.getExternalNodes() ) {
            double sum = 0;
            for( final String ref : refs ) {
                final Double d = PropertyColorScheme.parseNumber( PropertyColorScheme.valueFor( leaf, ref ) );
                if ( ( d != null ) && ( d.doubleValue() > 0 ) ) {
                    sum += d.doubleValue();
                }
            }
            max_total = Math.max( max_total, sum );
        }
        return new Column( type, refs, AptxUtil.distinctColors( refs.size() ), max_total, normalized );
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

    /** The glyph shape of SYMBOL column {@code i} (CIRCLE for a non-SYMBOL column). */
    SymbolShape symbolShape( final int i ) {
        return _columns.get( i ).getSymbolShape();
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

    /**
     * For a MERGED column (STACKED_BAR or PIE), {@code node}'s ordered per-series fractions, parallel to
     * {@link #stackColors(int)}: each series value divided by the scale -- the absolute max per-tip total (an absolute
     * stacked bar), or the tip's OWN total (a normalized stacked bar, and every pie -- so a pie's fractions ARE its
     * wedge angles). Missing / non-positive series values count as 0. The sum is {@code <= 1} (absolute: {@code = 1} at
     * the largest-total tip; normalized/pie: {@code = 1} for any tip with data, {@code 0} for a tip with none). Empty
     * when the column is not a merged type. A stacked bar draws each segment {@code k} from cumulative offset to
     * offset+fraction; a pie sweeps each wedge by {@code fraction * 360}; both skip a 0.
     */
    double[] stackFractions( final PhylogenyNode node, final int i ) {
        final Column c = _columns.get( i );
        if ( !isMergedType( c._type ) || ( c._stack_refs == null ) ) {
            return new double[ 0 ];
        }
        final double[] vals = new double[ c._stack_refs.size() ];
        double tip_total = 0;
        for( int k = 0; k < vals.length; ++k ) {
            final Double d = PropertyColorScheme.parseNumber( PropertyColorScheme.valueFor( node, c._stack_refs.get( k ) ) );
            final double v = ( ( d != null ) && ( d.doubleValue() > 0 ) ) ? d.doubleValue() : 0;
            vals[ k ] = v;
            tip_total += v;
        }
        final double denom = c._stack_normalized ? tip_total : c._stack_max;
        if ( denom <= 0 ) {
            return new double[ vals.length ]; // all zero -> nothing to draw (missing/absent tip, or an empty scale)
        }
        for( int k = 0; k < vals.length; ++k ) {
            vals[ k ] = vals[ k ] / denom;
        }
        return vals;
    }

    /** The series colours of STACKED_BAR column {@code i}, in series order (parallel to {@link #stackFractions} and
     *  {@link #stackHeaders}); empty for a non-stacked column. */
    List<Color> stackColors( final int i ) {
        final List<Color> colors = _columns.get( i )._stack_colors;
        return ( colors == null ) ? java.util.Collections.<Color> emptyList()
                : java.util.Collections.unmodifiableList( colors );
    }

    /** The series field headers (display names) of STACKED_BAR column {@code i}, in series order -- for the legend;
     *  empty for a non-stacked column. */
    List<String> stackHeaders( final int i ) {
        final List<String> refs = _columns.get( i )._stack_refs;
        if ( refs == null ) {
            return java.util.Collections.<String> emptyList();
        }
        final List<String> out = new ArrayList<String>();
        for( final String ref : refs ) {
            out.add( PropertyColorScheme.displayName( ref ) );
        }
        return out;
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
            types.add( Type.STACKED_BAR );
            types.add( Type.PIE );
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
