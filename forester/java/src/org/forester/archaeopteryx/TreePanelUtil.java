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
import java.awt.Component;
import java.awt.geom.AffineTransform;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Function;

import javax.swing.JOptionPane;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;
import org.forester.util.TaxonomyUtil;
import org.forester.ws.seqdb.AccessionAwareLineageService;
import org.forester.ws.seqdb.NcbiTaxonomyLineageService;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonomicLineageService;
import org.forester.ws.seqdb.WebOrganismSource;

public class TreePanelUtil {

    public final static String createUriForSeqWeb( final PhylogenyNode node,
                                                   final Configuration conf,
                                                   final TreePanel tp ) {
        String uri_str = null;
        final String upkb = SequenceAccessionTools.obtainUniProtAccessorFromDataFields( node );
        if ( !ForesterUtil.isEmpty( upkb ) ) {
            try {
                uri_str = ForesterUtil.UNIPROT_KB + URLEncoder.encode( upkb, ForesterConstants.UTF_8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                AptxUtil.showErrorMessage( tp, e.toString() );
                e.printStackTrace();
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = SequenceAccessionTools.obtainGenbankAccessorFromDataFields( node );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    if ( SequenceAccessionTools.isProteinDbQuery( v ) ) {
                        uri_str = ForesterUtil.NCBI_PROTEIN + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                    }
                    else {
                        uri_str = ForesterUtil.NCBI_NUCCORE + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                    }
                }
                catch ( final UnsupportedEncodingException e ) {
                    AptxUtil.showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = SequenceAccessionTools.obtainRefSeqAccessorFromDataFields( node );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    if ( SequenceAccessionTools.isProteinDbQuery( v ) ) {
                        uri_str = ForesterUtil.NCBI_PROTEIN + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                    }
                    else {
                        uri_str = ForesterUtil.NCBI_NUCCORE + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                    }
                }
                catch ( final UnsupportedEncodingException e ) {
                    AptxUtil.showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = SequenceAccessionTools.obtainGiNumberFromDataFields( node );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    uri_str = ForesterUtil.NCBI_GI + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                }
                catch ( final UnsupportedEncodingException e ) {
                    AptxUtil.showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        return uri_str;
    }

    public final static String createUriForSeqWeb( final Sequence seq,
                                                   final Configuration conf,
                                                   final TreePanel tp ) {
        String uri_str = null;
        final String upkb = SequenceAccessionTools.obtainUniProtAccessorFromSequence( seq );
        if ( !ForesterUtil.isEmpty( upkb ) ) {
            try {
                uri_str = ForesterUtil.UNIPROT_KB + URLEncoder.encode( upkb, ForesterConstants.UTF_8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                AptxUtil.showErrorMessage( tp, e.toString() );
                e.printStackTrace();
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = SequenceAccessionTools.obtainGenbankAccessorFromSequence( seq );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    if ( SequenceAccessionTools.isProteinDbQuery( v ) ) {
                        uri_str = ForesterUtil.NCBI_PROTEIN + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                    }
                    else {
                        uri_str = ForesterUtil.NCBI_NUCCORE + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                    }
                }
                catch ( final UnsupportedEncodingException e ) {
                    AptxUtil.showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = SequenceAccessionTools.obtainRefSeqAccessorFromSequence( seq );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    if ( SequenceAccessionTools.isProteinDbQuery( v ) ) {
                        uri_str = ForesterUtil.NCBI_PROTEIN + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                    }
                    else {
                        uri_str = ForesterUtil.NCBI_NUCCORE + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                    }
                }
                catch ( final UnsupportedEncodingException e ) {
                    AptxUtil.showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = SequenceAccessionTools.obtainGiNumberFromSequence( seq );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    uri_str = ForesterUtil.NCBI_GI + URLEncoder.encode( v, ForesterConstants.UTF_8 );
                }
                catch ( final UnsupportedEncodingException e ) {
                    AptxUtil.showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        return uri_str;
    }

    final public static void showInformationMessage( final Component parent, final String title, final String msg ) {
        JOptionPane.showMessageDialog( parent, msg, title, JOptionPane.INFORMATION_MESSAGE );
    }

    final static void collapseSubtree( final PhylogenyNode node, final boolean collapse ) {
        node.setCollapse( collapse );
        if ( node.isExternal() ) {
            return;
        }
        final PhylogenyNodeIterator it = new PreorderTreeIterator( node );
        while ( it.hasNext() ) {
            it.next().setCollapse( collapse );
        }
    }

    final static void uncollapseSubtree( final PhylogenyNode node ) {
        node.setCollapse( false );
        if ( node.isExternal() ) {
            return;
        }
        final PhylogenyNodeIterator it = new PreorderTreeIterator( node );
        while ( it.hasNext() ) {
            it.next().setCollapse( false );
        }
    }

    static void colorizeSubtree( final PhylogenyNode node, final BranchColor c ) {
        node.getBranchData().setBranchColor( c );
        final List<PhylogenyNode> descs = PhylogenyMethods.getAllDescendants( node );
        for( final PhylogenyNode desc : descs ) {
            desc.getBranchData().setBranchColor( c );
        }
    }

    // --- Node-symbol support visualization (see TreePanel.paintNodeSupportSymbol) -----------------
    // Support values come on different absolute scales -- posterior probabilities and aLRT in 0..1,
    // bootstrap and SH-aLRT in 0..100. We pick the scale ceiling from the data (anything above 1
    // implies the 0..100 family) rather than normalizing to the max observed value, so a given
    // symbol size/threshold means the same thing across trees.

    /** The support-scale ceiling implied by the largest value present: 1 (0-1 posteriors/aLRT), 100 (bootstrap
     *  / SH-aLRT), or 1000 -- the smallest of those that is >= the observed maximum. */
    final static double confidenceScaleMaxFor( final double observed_max ) {
        if ( observed_max > 100.0 ) {
            return 1000.0;
        }
        return ( observed_max > 1.0 ) ? 100.0 : 1.0;
    }

    /** Scans a tree's internal-node confidences and returns the implied scale ceiling (1, 100, or 1000). */
    final static double detectConfidenceScaleMax( final Phylogeny tree ) {
        double max = 0.0;
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isInternal() && n.getBranchData().isHasConfidences() ) {
                final double c = PhylogenyMethods.getConfidenceValue( n );
                if ( c > max ) {
                    max = c;
                }
            }
        }
        return confidenceScaleMaxFor( max );
    }

    /** Support as a fraction of the scale ceiling, clamped to 0..1. */
    final static double supportFraction( final double confidence, final double scale_max ) {
        if ( scale_max <= 0.0 ) {
            return 0.0;
        }
        final double f = confidence / scale_max;
        if ( f < 0.0 ) {
            return 0.0;
        }
        if ( f > 1.0 ) {
            return 1.0;
        }
        return f;
    }

    /** SIZE_SCALED diameter: linearly interpolates min..max by the support fraction. */
    final static float supportSymbolSize( final double confidence,
                                          final double scale_max,
                                          final float min_size,
                                          final float max_size ) {
        return (float) ( min_size + ( supportFraction( confidence, scale_max ) * ( max_size - min_size ) ) );
    }

    // The most a weakly-supported branch fades toward the background in COLOR_BRANCHES mode (a fraction-0
    // branch keeps 1 - this of its color, so it stays faint-but-visible rather than vanishing).
    private static final double SUPPORT_COLOR_MAX_FADE = 0.8;

    /**
     * COLOR_BRANCHES branch color: the full {@code strong} (branch) color at support {@code fraction}=1,
     * fading toward the {@code background} as support drops (theme-aware -- "weak support fades into the
     * background"). Pure; clamps the fraction to 0..1.
     */
    final static Color supportColor( final double fraction, final Color strong, final Color background ) {
        final double f = ( fraction < 0.0 ) ? 0.0 : ( ( fraction > 1.0 ) ? 1.0 : fraction );
        return blend( strong, background, SUPPORT_COLOR_MAX_FADE * ( 1.0 - f ) );
    }

    /** {@code a} blended {@code t} (clamped to 0..1) of the way toward {@code b}, per channel. */
    final static Color blend( final Color a, final Color b, final double t ) {
        final double tt = ( t < 0.0 ) ? 0.0 : ( ( t > 1.0 ) ? 1.0 : t );
        return new Color( (int) Math.round( a.getRed() + ( tt * ( b.getRed() - a.getRed() ) ) ),
                          (int) Math.round( a.getGreen() + ( tt * ( b.getGreen() - a.getGreen() ) ) ),
                          (int) Math.round( a.getBlue() + ( tt * ( b.getBlue() - a.getBlue() ) ) ) );
    }

    /**
     * X positions of the vertical distance grid lines: starting one {@code spacing} to the right of
     * {@code origin_x} (the distance-0 root) and stepping by {@code spacing} up to and including
     * {@code max_x} (the deepest tip). Empty when {@code spacing} is non-positive or the tree has no depth.
     */
    final static float[] scaleGridLineXs( final float origin_x, final float spacing, final float max_x ) {
        if ( spacing <= 0.0f ) {
            return new float[ 0 ];
        }
        final int n = (int) Math.floor( ( max_x - origin_x ) / spacing );
        if ( n <= 0 ) {
            return new float[ 0 ];
        }
        final float[] xs = new float[ n ];
        for ( int i = 0; i < n; i++ ) {
            xs[ i ] = origin_x + ( ( i + 1 ) * spacing );
        }
        return xs;
    }

    // A ceiling on the number of scale-axis ticks: on a pathological depth/spacing ratio (corrupt/raw-count branch
    // lengths, e.g. a max distance of billions) draw NO axis rather than allocate a giant array / overflow the int
    // count / hang the paint. Any real tree keeps depth/getScaleDistance() in the low tens.
    private static final int MAX_AXIS_TICKS = 1000;

    /**
     * The distance VALUES at which the labeled scale axis places ticks: 0 (the root), then stepping by
     * {@code spacing} up to and including {@code max_distance} (the deepest tip). Always includes 0 (the origin
     * tick, which {@link #scaleGridLineXs} skips). Empty when {@code spacing} is non-positive, the tree has no
     * depth, or the tick count would be absurd. Pure -- the caller maps each value to an x and a label.
     */
    final static double[] scaleAxisTickValues( final double max_distance, final double spacing ) {
        if ( ( spacing <= 0.0 ) || ( max_distance <= 0.0 ) ) {
            return new double[ 0 ];
        }
        // steps in units of spacing; +1e-9 tolerance (in steps) so a tick that lands on max_distance survives float
        // error (e.g. 3*0.1 = 0.30000000000000004). Checked BEFORE the int cast so a huge ratio can't overflow.
        final double steps = ( max_distance / spacing ) + 1.0e-9;
        if ( steps > MAX_AXIS_TICKS ) {
            return new double[ 0 ];
        }
        final int n = (int) Math.floor( steps );
        final double[] values = new double[ n + 1 ]; // + the 0 tick
        for ( int i = 0; i <= n; i++ ) {
            values[ i ] = i * spacing;
        }
        return values;
    }

    /** "Nice" calendar tick years covering {@code [from_year, to_year]} (each a multiple of a nice step, aiming for
     *  ~6-8 ticks). Steps are whole-year multiples (1, 2, 5, 10, 20, 25, 50, 100, ...) so the labels are calendar
     *  years, not decimals; a sub-year span still gets its enclosing whole years. Empty for a non-positive span. */
    final static double[] calendarTickYears( final double from_year, final double to_year ) {
        final double span = to_year - from_year;
        if ( span <= 0.0 ) {
            return new double[ 0 ];
        }
        final double step = niceYearStep( span / 7.0 );
        if ( step <= 0.0 ) {
            return new double[ 0 ];
        }
        final double first = Math.ceil( ( from_year / step ) - 1.0e-9 ) * step;
        final double count_d = Math.floor( ( ( to_year - first ) / step ) + 1.0e-9 ) + 1;
        if ( ( count_d <= 0 ) || ( count_d > MAX_AXIS_TICKS ) ) {
            return new double[ 0 ];
        }
        final int count = (int) count_d;
        final double[] out = new double[ count ];
        for ( int i = 0; i < count; ++i ) {
            out[ i ] = first + ( i * step );
        }
        return out;
    }

    /** The smallest WHOLE-YEAR "nice" step (1/2/5 x 10^k) that is &gt;= {@code raw}; at least 1 year. Only whole-year
     *  multipliers are used (no 2.5, which would give a fractional 2.5-year step at pow=1 and mislabel half-year ticks). */
    static double niceYearStep( final double raw ) {
        if ( raw <= 1.0 ) {
            return 1.0;
        }
        double pow = 1.0;
        while ( ( pow * 10.0 ) <= raw ) {
            pow *= 10.0;
        }
        for ( final double m : new double[] { 1, 2, 5, 10 } ) {
            if ( ( pow * m ) >= raw ) {
                return pow * m;
            }
        }
        return pow * 10.0;
    }

    /** The smallest "nice" step (1/2/5 x 10^k, k ANY integer incl. NEGATIVE for sub-1 steps) that is &gt;= {@code raw}.
     *  Unlike {@link #niceYearStep} this allows fractional steps (0.5, 0.2, 0.1 ...), so a shallow geologic tree (root a
     *  few Ma old) gets fine ticks and a deep one gets 50/100/500 Ma ticks. 0 for a non-positive raw. */
    static double niceAxisStep( final double raw ) {
        if ( !( raw > 0.0 ) ) {
            return 0.0;
        }
        final double pow = Math.pow( 10.0, Math.floor( Math.log10( raw ) ) );
        for ( final double m : new double[] { 1, 2, 5, 10 } ) {
            if ( ( pow * m ) >= ( raw - ( raw * 1.0e-9 ) ) ) {
                return pow * m;
            }
        }
        return pow * 10.0;
    }

    /** "Nice" age ticks (Ma before present) for a numeric geologic axis over {@code [0, root_age]}: 0 (the present) and
     *  multiples of a nice step (aiming for ~7 ticks) up to {@code root_age}. Reuses {@link #scaleAxisTickValues} with a
     *  {@link #niceAxisStep} spacing (so sub-1 Ma steps are allowed). Empty for a non-positive root age. */
    final static double[] maAxisTickValues( final double root_age ) {
        return scaleAxisTickValues( root_age, niceAxisStep( root_age / 8.0 ) );
    }

    /**
     * The device-y the horizontal scale axis line is drawn at (its top). On SCREEN the axis FLOATS at the viewport
     * bottom so it never scrolls out of view when zoomed in (PearTree-style), exactly like the viewport-fixed scale
     * bar. A FILE export stays anchored to the tree/export extent bottom so figures remain WYSIWYG; the direct
     * File&gt;Print path (an export flag set but {@code graphics_file_height == 0}) anchors to the whole canvas.
     */
    final static int scaleAxisFloatingBottom( final boolean to_pdf,
                                              final boolean to_graphics_file,
                                              final int graphics_file_y,
                                              final int graphics_file_height,
                                              final int canvas_height,
                                              final int viewport_bottom ) {
        if ( to_pdf || to_graphics_file ) {
            return ( graphics_file_height > 0 ) ? ( graphics_file_y + graphics_file_height ) : canvas_height;
        }
        return viewport_bottom;
    }

    /**
     * The device-x the VERTICAL-orientation scale ruler is pinned to. On SCREEN the ruler FLOATS to the viewport
     * breadth EDGE on its own side -- the side AWAY from the tree, given by {@code in} (&gt;0 = tree to the right, so
     * the ruler sits at the left edge; &lt;=0 = tree to the left, ruler at the right edge) -- so it stays visible when
     * the breadth is zoomed/scrolled. A FILE export or File&gt;Print keeps the tree-anchored breadth position so
     * figures remain WYSIWYG.
     */
    final static int scaleAxisRulerX( final boolean to_pdf,
                                      final boolean to_graphics_file,
                                      final int tree_anchored_x,
                                      final int in,
                                      final int viewport_x,
                                      final int viewport_width ) {
        if ( to_pdf || to_graphics_file ) {
            return tree_anchored_x;
        }
        return ( in > 0 ) ? ( viewport_x + 1 ) : ( ( viewport_x + viewport_width ) - 1 );
    }

    /**
     * The horizontal x-range {@code [left, right]} of a node-age (HPD) bar, anchored to the node's OWN drawn x
     * ({@code node_x}) and offset by the signed age deltas: an older bound ({@code max}) sits to the LEFT of the node,
     * a younger bound ({@code min}) to the RIGHT, each scaled by {@code corr} (px per branch-length/time unit). Using
     * the node's real position (not a tree-height-derived age->x map) keeps the bar centred on the node it annotates
     * even when the tree is NOT strictly ultrametric or the root carries a branch length. {@code value} is the node's
     * point age. Pure/testable.
     */
    final static float[] hpdBarXRange( final float node_x, final double value, final double min, final double max,
                                       final double corr ) {
        final float left = (float) ( node_x - ( ( max - value ) * corr ) );
        final float right = (float) ( node_x + ( ( value - min ) * corr ) );
        return new float[] { left, right };
    }

    /**
     * The half-thickness of a node-age SPINDLE at position {@code p} along its span {@code [lo, hi]} (either device x
     * for the rectangular shape or radius for the circular one): 0 at both ends {@code lo}/{@code hi}, rising smoothly
     * (a quarter-sine) to {@code h_max} at the point-estimate position {@code peak}. Asymmetric when the estimate is
     * off-centre in the HPD interval -- so the spindle shows WHERE the point age sits within its 95% HPD, unlike the
     * flat bar. This is a schematic of the SUMMARISED uncertainty (point + 95% HPD), not the raw posterior density.
     * Pure/testable.
     */
    /**
     * The pixel height of a protein-domain-architecture box for a leaf row of spacing {@code row_spacing}
     * (= {@code TreePanel.getYdistance()}): track the vertical zoom, but clamp into a legible band
     * [{@code min}, {@code max}] so the boxes stay readable when zoomed out (floor) and don't become giant blocks
     * when zoomed in (ceiling). Extracted + shared by BOTH domain-height sites in {@code TreePanel} so the clamp is
     * testable and cannot drift between them. Pure.
     */
    final static int domainBoxHeight( final float row_spacing, final int min, final int max ) {
        return Math.max( min, Math.min( ForesterUtil.roundToInt( row_spacing ), max ) );
    }

    /**
     * Trim {@code s} (binary search) so it renders within {@code max_px} pixels in {@code fm}, appending an ellipsis;
     * returns the whole string if it already fits, "" if not even the ellipsis fits, and s for null/empty. Used to
     * truncate radial (circular/unrooted) tip labels so tree + labels + domains stay on-canvas. Pure.
     */
    final static String truncateToPixelWidth( final java.awt.FontMetrics fm, final String s, final int max_px ) {
        if ( ( s == null ) || ( s.length() == 0 ) || ( fm.stringWidth( s ) <= max_px ) ) {
            return s;
        }
        final String ell = "…";
        final int ell_w = fm.stringWidth( ell );
        if ( max_px <= ell_w ) {
            return "";
        }
        int lo = 0, hi = s.length();
        while ( lo < hi ) {
            final int mid = ( lo + hi + 1 ) / 2;
            if ( ( fm.stringWidth( s.substring( 0, mid ) ) + ell_w ) <= max_px ) {
                lo = mid;
            }
            else {
                hi = mid - 1;
            }
        }
        return s.substring( 0, lo ) + ell;
    }

    final static double spindleHalfHeightAt( final double p, final double lo, final double hi, final double peak,
                                             final double h_max ) {
        if ( ( hi <= lo ) || ( h_max <= 0 ) ) {
            return 0;
        }
        final double pk = Math.min( Math.max( peak, lo ), hi ); // clamp the estimate into the interval
        double t;
        if ( p <= pk ) {
            final double span = pk - lo;
            t = ( span > 0 ) ? ( ( p - lo ) / span ) : 1.0; // estimate at the low end -> that lobe is a step to full height
        }
        else {
            final double span = hi - pk;
            t = ( span > 0 ) ? ( ( hi - p ) / span ) : 1.0;
        }
        t = Math.min( Math.max( t, 0 ), 1 );
        return h_max * Math.sin( ( t * Math.PI ) / 2.0 ); // 0 at the tips, h_max at the peak
    }

    /**
     * Formats a number for a compact figure label (scale-axis ticks, size-legend samples): a whole number as an
     * integer, otherwise with enough decimals to stay legible across magnitudes -- 2 decimals for values &gt;= 1
     * (years, distances, counts) but MORE for small magnitudes, so a 0..1 property/distance does not collapse to
     * "0"/duplicate labels. Trailing zeros dropped. US-locale (reproducible across locales; a fresh formatter per
     * call -- DecimalFormat is not thread-safe, so a shared static would be an off-EDT hazard).
     */
    final static String formatCompactNumber( final double v ) {
        if ( ( v == Math.rint( v ) ) && !Double.isInfinite( v ) && ( Math.abs( v ) < 1e15 ) ) {
            return Long.toString( (long) v );
        }
        final double abs = Math.abs( v );
        int decimals = ( ( abs >= 1 ) || ( abs == 0 ) ) ? 2 : ( 2 - (int) Math.floor( Math.log10( abs ) ) );
        decimals = Math.max( 0, Math.min( decimals, 10 ) );
        final StringBuilder pattern = new StringBuilder( "0." );
        for ( int i = 0; i < decimals; ++i ) {
            pattern.append( '#' ); // '#' drops trailing zeros
        }
        return new java.text.DecimalFormat( pattern.toString(),
                                            java.text.DecimalFormatSymbols.getInstance( java.util.Locale.US ) )
                .format( v );
    }

    /** THRESHOLD_MARKS test: is the support at or above the cutoff (a fraction 0..1 of the scale)? */
    final static boolean isSupportAtOrAboveThreshold( final double confidence,
                                                      final double scale_max,
                                                      final double threshold_fraction ) {
        return supportFraction( confidence, scale_max ) >= threshold_fraction;
    }

    /**
     * The {@code {x, y}} center at which a branch-support symbol is drawn: the middle of the branch
     * (parent&rarr;node), since support is a branch property. The horizontal x is always the branch
     * midpoint. For {@code radial} (unrooted/circular) layouts the branch is a slanted segment, so the y
     * is the segment midpoint too; for the rectangular layouts the branch is a horizontal segment at the
     * node's y, so the y is simply {@code node_y}.
     */
    final static float[] supportSymbolCenter( final float parent_x,
                                              final float node_x,
                                              final float parent_y,
                                              final float node_y,
                                              final boolean radial ) {
        final float cx = ( parent_x + node_x ) / 2.0f;
        final float cy = radial ? ( ( parent_y + node_y ) / 2.0f ) : node_y;
        return new float[] { cx, cy };
    }

    /**
     * The support-symbol centre for a CIRCULAR layout: the branch to a node is drawn as a RADIAL leg along the node's
     * own spoke (from the node's radius inward to the parent's radius) plus an arc, so the straight-line
     * (Cartesian) midpoint of node and parent lies OFF the branch. This puts the symbol on the node's spoke at the
     * mid-radius between node and parent -- i.e. the midpoint of the radial leg -- where the branch actually is.
     * {@code (root_x, root_y)} is the ring centre (the circular layout places the root there). Pure.
     */
    final static float[] circularSupportSymbolCenter( final float root_x,
                                                      final float root_y,
                                                      final float node_x,
                                                      final float node_y,
                                                      final float parent_x,
                                                      final float parent_y ) {
        final double dnx = node_x - root_x;
        final double dny = node_y - root_y;
        final double r_node = Math.sqrt( ( dnx * dnx ) + ( dny * dny ) );
        if ( r_node <= 0.0 ) {
            return new float[] { node_x, node_y }; // node sits at the centre (no spoke to place it on)
        }
        final double r_parent = Math.sqrt( ( ( parent_x - root_x ) * ( parent_x - root_x ) )
                + ( ( parent_y - root_y ) * ( parent_y - root_y ) ) );
        final double mid_r = ( r_node + r_parent ) / 2.0;
        return new float[] { (float) ( root_x + ( ( dnx / r_node ) * mid_r ) ),
                             (float) ( root_y + ( ( dny / r_node ) * mid_r ) ) };
    }

    /**
     * Draw positions for an internal node's label placed to the LEFT of the node, right-aligned so it
     * ends just left of the node and sits on top of the incoming branch (the publication-style
     * placement). The label is two adjacent segments read left-to-right: an optional taxonomy segment
     * then an optional node-data segment, with the node-data segment's right edge at the node. Returns
     * {@code {taxo_x, data_x, baseline_y}}: the left x at which to draw each segment and the shared text
     * baseline. The inter-segment {@code gap} is only applied when both segments are present.
     *
     * <p>If right-alignment would push the label's leftmost glyph left of {@code min_x} (a long label on
     * an internal node near the root), the whole label is shifted right to start at {@code min_x} so it
     * stays on-canvas rather than being clipped -- it then extends rightward from {@code min_x} instead
     * of ending exactly at the node.
     */
    final static float[] internalLabelAboveBranchLayout( final float node_x,
                                                         final float node_y,
                                                         final int half_box_size,
                                                         final int taxo_width,
                                                         final int data_width,
                                                         final int gap,
                                                         final int font_descent,
                                                         final float min_x ) {
        // "- 2" is the small gap between the node and the label's right edge (mirrors the classic
        // right-of-node path's "+ 2 + half_box_size"); "- 1" on the baseline lifts the glyph bottoms
        // just clear of the horizontal branch line at node_y (screen y grows downward).
        final float right = node_x - half_box_size - 2;
        float data_x = right - data_width;
        final int effective_gap = ( ( taxo_width > 0 ) && ( data_width > 0 ) ) ? gap : 0;
        float taxo_x = data_x - effective_gap - taxo_width;
        final float leftmost = ( taxo_width > 0 ) ? taxo_x : data_x;
        if ( leftmost < min_x ) {
            final float shift = min_x - leftmost;
            data_x += shift;
            taxo_x += shift;
        }
        final float baseline_y = node_y - font_descent - 1;
        return new float[] { taxo_x, data_x, baseline_y };
    }

    /**
     * Abbreviates a binomial scientific name to the genus initial + ". " + the full species epithet, per
     * the standard convention (e.g. {@code "Homo sapiens"} &rarr; {@code "H. sapiens"}); any further
     * epithets are kept verbatim ({@code "Homo sapiens neanderthalensis"} &rarr;
     * {@code "H. sapiens neanderthalensis"}). Display-only: the caller never writes this back to the
     * taxonomy. A name that is not an abbreviatable binomial -- fewer than two whitespace-separated tokens
     * or an empty first token (leading whitespace) -- is returned unchanged rather than throwing.
     */
    final static String abbreviateScientificName( final String scientific_name ) {
        final String[] a = scientific_name.split( "\\s+" );
        if ( ( a.length < 2 ) || a[ 0 ].isEmpty() ) {
            return scientific_name;
        }
        final StringBuilder sb = new StringBuilder();
        sb.append( a[ 0 ].charAt( 0 ) );
        sb.append( ". " );
        sb.append( a[ 1 ] );
        for( int i = 2; i < a.length; i++ ) {
            sb.append( " " );
            sb.append( a[ i ] );
        }
        return sb.toString();
    }

    // Ranks whose scientific name is italicized by convention (genus and below). Family/order/etc. are roman.
    private static final java.util.Set<String> ITALIC_RANKS = new java.util.HashSet<String>( java.util.Arrays.asList(
            "genus", "subgenus", "species", "subspecies", "species group", "species subgroup", "varietas", "variety",
            "subvariety", "forma", "form", "form species", "strain", "isolate", "biotype", "serotype", "pathovar" ) );
    // Suffixes of higher-rank names (family/order/...), to keep them roman when a tree carries no rank field.
    private static final String[] HIGHER_RANK_SUFFIXES = { "idae", "aceae", "inae", "oideae", "acea", "ineae", "ales",
            "oidea", "ini", "formes", "phyta", "mycota", "phyceae", "opsida", "viridae", "virales", "viricetes",
            "viricota", "virinae" };

    /**
     * Whether a taxonomy scientific name is italicized by taxonomic convention -- true only for genus/species-level
     * names, so family/order/higher names stay upright. Uses {@code rank} when present; otherwise a shape/suffix
     * heuristic (a binomial is a species; a single word ending in a known higher-rank suffix is not italic, else it
     * is assumed to be a genus). Does NOT consider the display option (the caller gates on that).
     */
    final static boolean scientificNameIsItalic( final String rank, final String name ) {
        if ( !ForesterUtil.isEmpty( rank ) ) {
            return ITALIC_RANKS.contains( rank.trim().toLowerCase( java.util.Locale.ROOT ) );
        }
        if ( ForesterUtil.isEmpty( name ) ) {
            return false;
        }
        final String n = name.trim();
        final int sp = n.indexOf( ' ' );
        if ( sp > 0 ) {
            return Character.isLowerCase( n.charAt( sp + 1 ) ); // "Genus epithet" -> a species name
        }
        final String lc = n.toLowerCase( java.util.Locale.ROOT );
        for( final String suffix : HIGHER_RANK_SUFFIXES ) {
            if ( lc.endsWith( suffix ) ) {
                return false;
            }
        }
        return true; // a single word that is not a known higher-rank name -> assume a genus
    }

    /**
     * Whether a node's own name merely repeats a taxonomy label already being shown (so the same word is not drawn
     * twice). Compares case-insensitively against the scientific / common name, each only when that field is shown.
     */
    final static boolean nodeNameDuplicatesTaxonomy( final String node_name, final String scientific_name,
                                                     final String common_name, final boolean show_scientific,
                                                     final boolean show_common ) {
        if ( ForesterUtil.isEmpty( node_name ) ) {
            return false;
        }
        final String n = node_name.trim();
        if ( show_scientific && !ForesterUtil.isEmpty( scientific_name )
                && n.equalsIgnoreCase( scientific_name.trim() ) ) {
            return true;
        }
        if ( show_common && !ForesterUtil.isEmpty( common_name ) && n.equalsIgnoreCase( common_name.trim() ) ) {
            return true;
        }
        return false;
    }

    /** The best display label for a taxonomy: scientific name, else common name, else taxonomy code, else "". */
    final static String taxonomyLabel( final Taxonomy tax ) {
        if ( tax != null ) {
            if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                return tax.getScientificName();
            }
            if ( !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
                return tax.getCommonName();
            }
            if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                return tax.getTaxonomyCode();
            }
        }
        return "";
    }

    /**
     * True when {@code node} is a non-root INTERNAL node whose <i>visible</i> taxonomy label equals that of the
     * nearest ancestor with a visible label -- i.e. its label would just repeat an enclosing clade's. Used to
     * suppress the redundant label at draw time (a clade is marked once, at its topmost node), which declutters
     * nested same-taxon clades (common after ancestral-taxonomy inference, e.g. a Boreoeutheria node inside a
     * Boreoeutheria node) WITHOUT mutating the tree. Tips, collapsed nodes and the root never qualify.
     * <p>
     * {@code labeler} yields the exact string the paint path would draw for a node (which depends on the active
     * taxonomy checkboxes -- rank/code/scientific/common -- so two nodes sharing a scientific name but rendering a
     * different code/rank are NOT judged equal). Ancestors that render an empty label are skipped so the walk
     * reaches the nearest VISIBLE ancestor label.
     */
    final static boolean isDuplicateOfAncestorTaxon( final PhylogenyNode node,
                                                     final Function<PhylogenyNode, String> labeler ) {
        if ( ( node == null ) || node.isExternal() || node.isCollapse() || node.isRoot()
                || !node.getNodeData().isHasTaxonomy() ) {
            return false;
        }
        final String own = labeler.apply( node );
        if ( ForesterUtil.isEmpty( own ) ) {
            return false;
        }
        for( PhylogenyNode a = node.getParent(); a != null; a = a.getParent() ) {
            if ( !a.getNodeData().isHasTaxonomy() ) {
                continue;
            }
            final String anc = labeler.apply( a );
            if ( !ForesterUtil.isEmpty( anc ) ) {
                return own.equalsIgnoreCase( anc ); // nearest ancestor with a visible label
            }
            // an ancestor with a Taxonomy but no visible label (e.g. only a tax-id) -- keep walking
        }
        return false;
    }

    // ===== Ancestral-state pie charts (BEAST discrete/geographic traits) =====

    private final static String BEAST_PREFIX    = "beast:";
    private final static String NEXTSTRAIN_PREFIX = "nextstrain:";
    // discrete ancestral-state distributions come from BEAST/TreeAnnotator (beast:) or Auspice/Nextstrain (nextstrain:),
    // both stored as a <prefix><trait>_set + _set_prob brace-list pair; the pies read either.
    private final static String[] STATE_PREFIXES = { BEAST_PREFIX, NEXTSTRAIN_PREFIX };
    private final static String SET_SUFFIX      = "_set";
    private final static String SET_PROB_SUFFIX = "_set_prob";

    /** A discrete state and its posterior probability at a node (a wedge of the ancestral-state pie). Immutable. */
    public final static class StateProbability {

        private final String _state;
        private final double _probability;

        StateProbability( final String state, final double probability ) {
            _state = state;
            _probability = probability;
        }

        public String getState() {
            return _state;
        }

        public double getProbability() {
            return _probability;
        }
    }

    /** The discrete-trait names for which the tree carries an ancestral-state DISTRIBUTION -- a node with both a
     *  {@code beast:<trait>_set} and a {@code beast:<trait>_set_prob} property (BEAST/TreeAnnotator). These are the
     *  traits offerable as ancestral-state pies. */
    final static SortedSet<String> ancestralStateTraits( final Phylogeny tree ) {
        final SortedSet<String> traits = new TreeSet<String>();
        if ( ( tree == null ) || tree.isEmpty() ) {
            return traits;
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().getProperties() == null ) {
                continue;
            }
            for( final Property p : n.getNodeData().getProperties().getProperties() ) {
                final String ref = p.getRef();
                for( final String prefix : STATE_PREFIXES ) {
                    if ( ref.startsWith( prefix ) && ref.endsWith( SET_PROB_SUFFIX ) ) {
                        final String trait = ref.substring( prefix.length(),
                                                            ref.length() - SET_PROB_SUFFIX.length() );
                        if ( !ForesterUtil.isEmpty( trait )
                                && !firstProperty( n, prefix + trait + SET_SUFFIX ).isEmpty() ) {
                            traits.add( trait ); // require the matching states list too (a distribution is set + set_prob)
                        }
                    }
                }
            }
        }
        return traits;
    }

    /** The ordered {@code {state, probability}} distribution of {@code trait} at {@code node}: from
     *  {@code beast:<trait>_set} + {@code beast:<trait>_set_prob} (zipped, probabilities normalized to sum 1), or a
     *  single {@code {state, 1.0}} from {@code beast:<trait>} (a tip's observed state), or empty. Pure. */
    final static List<StateProbability> stateDistribution( final PhylogenyNode node, final String trait ) {
        final List<StateProbability> out = new ArrayList<StateProbability>();
        if ( ( node == null ) || ForesterUtil.isEmpty( trait ) || ( node.getNodeData().getProperties() == null ) ) {
            return out;
        }
        String set_raw = "";
        String prob_raw = "";
        for( final String prefix : STATE_PREFIXES ) { // beast: or nextstrain:, whichever carries this trait's distribution
            final String s = firstProperty( node, prefix + trait + SET_SUFFIX );
            final String p = firstProperty( node, prefix + trait + SET_PROB_SUFFIX );
            if ( !s.isEmpty() || !p.isEmpty() ) {
                set_raw = s;
                prob_raw = p;
                break;
            }
        }
        // If the node carries EITHER set property it is MEANT to have a distribution: require a well-formed
        // (equal-length, all-finite, all-non-negative) pair, else return empty (no pie) -- never fall through to
        // the single-state disc below, which would misrepresent an uncertain distribution as a confident 100%.
        if ( !set_raw.isEmpty() || !prob_raw.isEmpty() ) {
            final List<String> states = parseBraceList( set_raw );
            final List<String> probs = parseBraceList( prob_raw );
            if ( states.isEmpty() || ( states.size() != probs.size() ) ) {
                return out; // malformed / length-mismatched -> no pie (never half-drawn, never overstated)
            }
            final double[] p = new double[ probs.size() ];
            double sum = 0.0;
            for( int i = 0; i < probs.size(); i++ ) {
                final Double d = parseFinite( probs.get( i ) );
                if ( ( d == null ) || ( d.doubleValue() < 0.0 ) ) {
                    return out; // non-numeric / NaN / Infinity / negative probability -> no pie
                }
                p[ i ] = d.doubleValue();
                sum += p[ i ];
            }
            if ( sum <= 0.0 ) {
                return out;
            }
            for( int i = 0; i < states.size(); i++ ) {
                if ( p[ i ] > 0.0 ) {
                    out.add( new StateProbability( states.get( i ), p[ i ] / sum ) );
                }
            }
            return out;
        }
        for( final String prefix : STATE_PREFIXES ) {
            final String single = firstProperty( node, prefix + trait ); // a tip's single observed state
            if ( !ForesterUtil.isEmpty( single ) ) {
                out.add( new StateProbability( single, 1.0 ) );
                return out;
            }
        }
        return out;
    }

    /** Every distinct state of {@code trait} across the tree (from the per-node distributions + single states), for a
     *  stable state&rarr;color assignment shared by all pies and the legend. */
    final static SortedSet<String> collectAncestralStates( final Phylogeny tree, final String trait ) {
        final SortedSet<String> states = new TreeSet<String>();
        if ( ( tree == null ) || tree.isEmpty() || ForesterUtil.isEmpty( trait ) ) {
            return states;
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            for( final StateProbability sp : stateDistribution( it.next(), trait ) ) {
                states.add( sp.getState() );
            }
        }
        return states;
    }

    /** Split a BEAST {@code {a,b,c}} (or {@code [a,b,c]}) list into trimmed, unquoted elements -- top-level commas
     *  only (a nested brace/bracket or a DOUBLE-quoted value is not split). A single-quote/apostrophe is treated as a
     *  LITERAL character (not a delimiter), so a state name like {@code Xi'an} or {@code Côte d'Ivoire} -- quoted by
     *  BEAST with double quotes, or unquoted -- is not swallowed. Empty for a non-list / empty input. */
    final static List<String> parseBraceList( final String raw ) {
        final List<String> out = new ArrayList<String>();
        if ( ForesterUtil.isEmpty( raw ) ) {
            return out;
        }
        String s = raw.trim();
        if ( ( s.length() >= 2 ) && ( ( s.charAt( 0 ) == '{' ) || ( s.charAt( 0 ) == '[' ) )
                && ( ( s.charAt( s.length() - 1 ) == '}' ) || ( s.charAt( s.length() - 1 ) == ']' ) ) ) {
            s = s.substring( 1, s.length() - 1 );
        }
        int depth = 0;
        boolean in_dq = false;
        final StringBuilder cur = new StringBuilder();
        for( int i = 0; i < s.length(); i++ ) {
            final char c = s.charAt( i );
            if ( in_dq ) {
                if ( c == '"' ) {
                    in_dq = false;
                }
                else {
                    cur.append( c );
                }
            }
            else if ( c == '"' ) {
                in_dq = true;
            }
            else if ( ( c == '{' ) || ( c == '[' ) ) {
                depth++;
                cur.append( c );
            }
            else if ( ( c == '}' ) || ( c == ']' ) ) {
                if ( depth > 0 ) {
                    depth--;
                }
                cur.append( c );
            }
            else if ( ( c == ',' ) && ( depth == 0 ) ) {
                addTrimmed( out, cur );
                cur.setLength( 0 );
            }
            else {
                cur.append( c );
            }
        }
        addTrimmed( out, cur );
        return out;
    }

    private final static void addTrimmed( final List<String> out, final StringBuilder sb ) {
        final String v = sb.toString().trim();
        if ( !v.isEmpty() ) {
            out.add( v );
        }
    }

    /** The value of the first node property with ref {@code ref}, or "". */
    private final static String firstProperty( final PhylogenyNode node, final String ref ) {
        if ( node.getNodeData().getProperties() == null ) {
            return "";
        }
        final List<Property> ps = node.getNodeData().getProperties().getProperties( ref );
        return ps.isEmpty() ? "" : ps.get( 0 ).getValue();
    }

    private final static Double parseFinite( final String s ) {
        try {
            final double d = Double.parseDouble( s );
            return ( Double.isNaN( d ) || Double.isInfinite( d ) ) ? null : Double.valueOf( d );
        }
        catch ( final NumberFormatException e ) {
            return null;
        }
    }

    /** An evenly-spaced gray for element {@code i} of {@code n} (0-based), for a monochrome (black-and-white export)
     *  rendering of a categorical color set -- e.g. ancestral-pie wedges + their legend swatches -- so the elements
     *  stay distinguishable and the two match. A single element (or n&le;1) is a neutral mid-gray. */
    final static Color grayShade( final int i, final int n ) {
        final int shade = ( n <= 1 ) ? 128 : ( 40 + Math.round( ( i / (float) ( n - 1 ) ) * 180f ) );
        return new Color( shade, shade, shade );
    }

    /**
     * A taxon at a rank, identified by its NCBI tax-id when available (else by name) -- the grouping key (Spine B)
     * shared by the rank branch-colorizer and the clade bands. Keying on the tax-id (Spine A's
     * {@link TaxonLineage#taxIdAt}) rather than the scientific name distinguishes HOMONYMS (two distinct taxa that
     * share a name at a rank), while an id-less source (a name-only tree, or a resolver that supplies no id) falls
     * back to name equality -- the pre-Spine-B behavior. The {@link #getName() name} is always kept for display, so
     * the legend / color / override chain stays name-keyed.
     */
    final static class RankTaxon {

        private final String _id;   // NCBI tax-id, or null
        private final String _name; // display name (non-empty for a real taxon)

        RankTaxon( final String id, final String name ) {
            _id = ForesterUtil.isEmpty( id ) ? null : id;
            _name = ( name == null ) ? "" : name;
        }

        String getId() {
            return _id;
        }

        String getName() {
            return _name;
        }

        /** The equality / grouping key: the tax-id when present, else the lower-cased name. */
        private String key() {
            return ( _id != null ) ? ( "id:" + _id ) : ( "nm:" + _name.toLowerCase( Locale.ROOT ) );
        }

        @Override
        public boolean equals( final Object o ) {
            return ( o instanceof RankTaxon ) && key().equals( ( (RankTaxon) o ).key() );
        }

        @Override
        public int hashCode() {
            return key().hashCode();
        }

        @Override
        public String toString() {
            return _name + ( ( _id != null ) ? ( " [" + _id + "]" ) : "" );
        }
    }

    /** Sentinel for {@link #maximalMonochromaticRoots}: a subtree whose descendants are not all one rank taxon
     *  (compared by reference identity, never equal to a real {@link RankTaxon}). */
    private final static RankTaxon MIXED = new RankTaxon( null, "<MIXED-sentinel>" );

    private static TaxonomicLineageService _default_lineage_service;

    /**
     * The process-wide {@link TaxonomicLineageService} used by the rank colorizer and "Annotate Clades by
     * Rank". It wraps the shared NCBI taxonomy singleton (whose in-memory + persistent caches it shares
     * with the Fetch tool and the Settings cache panel) in an {@link AccessionAwareLineageService}, so
     * tips identified by a UniProt/SwissProt/RefSeq/GenBank/GI <i>sequence</i> accession -- which a bare
     * taxonomy database cannot place -- are resolved to their organism (taxonomy-only; the full protein
     * record is never cached) first. Trees with UniProt and/or mixed NCBI/UniProt identifiers are very
     * common.
     */
    final static synchronized TaxonomicLineageService getDefaultLineageService() {
        if ( _default_lineage_service == null ) {
            _default_lineage_service = new AccessionAwareLineageService( NcbiTaxonomyLineageService.getShared(),
                                                                         new WebOrganismSource() );
        }
        return _default_lineage_service;
    }

    /**
     * Colorizes the tree by taxonomic {@code rank}: every external node is assigned to the taxon it
     * belongs to at {@code rank} (from in-tree rank annotations first, then the {@code service}'s
     * cached lineages), then each maximal clade whose tips all share one such taxon is colored with a
     * distinct color. Unlike the old "color the subtree of any node literally annotated at the rank"
     * approach this places a genus-only tip (e.g. <i>Felis</i>) under its order (Carnivora) and
     * colors paraphyletic groups as several same-colored runs. When {@code legend_out} is non-null it
     * is filled with the taxon&rarr;color pairs used. Returns the number of colored clades.
     *
     * <p>Network-pure: it only reads {@code service}'s cache ({@link TaxonomicLineageService#lineageOf})
     * and never fetches, so it is safe on the EDT and unit-testable with an in-memory service. Callers
     * fetch unresolved taxa (see {@link #unresolvedTipTaxa}) off the EDT first, then call this again.
     */
    final static int colorPhylogenyAccordingToRanks( final Phylogeny tree,
                                                     final String rank,
                                                     final TaxonomicLineageService service,
                                                     final Map<String, Color> legend_out ) {
        return colorPhylogenyAccordingToRanks( tree, rank, service, legend_out, null );
    }

    /** {@code overrides} (taxon -&gt; user-chosen color) replaces the auto-assigned color for those taxa. */
    final static int colorPhylogenyAccordingToRanks( final Phylogeny tree,
                                                     final String rank,
                                                     final TaxonomicLineageService service,
                                                     final Map<String, Color> legend_out,
                                                     final Map<String, Color> overrides ) {
        return colorPhylogenyAccordingToRanks( tree, rank, service, legend_out, overrides, null );
    }

    /**
     * {@code counts_out}, when non-null, receives each legend taxon's tip count (how many tips were assigned
     * to it at {@code rank}), so the legend can show "(N)" and sort by count.
     */
    final static int colorPhylogenyAccordingToRanks( final Phylogeny tree,
                                                     final String rank,
                                                     final TaxonomicLineageService service,
                                                     final Map<String, Color> legend_out,
                                                     final Map<String, Color> overrides,
                                                     final Map<String, Integer> counts_out ) {
        final Map<PhylogenyNode, RankTaxon> assignment = assignNodesToRankTaxon( tree, rank, service );
        final Map<PhylogenyNode, RankTaxon> roots = maximalMonochromaticRoots( tree, assignment );
        final Map<String, Color> colors = AptxUtil.assignDistinctColors( legendTaxa( assignment, roots ) );
        applyColorOverrides( colors, overrides );
        int colorizations = 0;
        for( final Entry<PhylogenyNode, RankTaxon> e : roots.entrySet() ) {
            final Color c = colors.get( e.getValue().getName() );
            if ( c != null ) {
                TreePanelUtil.colorizeSubtree( e.getKey(), new BranchColor( c ) );
                ++colorizations;
            }
        }
        if ( legend_out != null ) {
            legend_out.putAll( colors );
        }
        countTipsPerTaxon( assignment, counts_out );
        return colorizations;
    }

    /** The distinct taxon NAMES for the legend / color domain: every TIP's taxon (so the key shows what is present,
     *  as before Spine B) PLUS every colored ROOT's taxon (so a gap-filled internal-only clade still gets a color). */
    private static SortedSet<String> legendTaxa( final Map<PhylogenyNode, RankTaxon> assignment,
                                                 final Map<PhylogenyNode, RankTaxon> roots ) {
        final SortedSet<String> taxa = new TreeSet<String>();
        for( final Entry<PhylogenyNode, RankTaxon> e : assignment.entrySet() ) {
            if ( e.getKey().isExternal() ) {
                taxa.add( e.getValue().getName() );
            }
        }
        for( final RankTaxon rt : roots.values() ) {
            taxa.add( rt.getName() );
        }
        return taxa;
    }

    /** Fills {@code counts_out} (when non-null) with the number of TIPS assigned to each taxon name (internal-node
     *  assignments feed the grouping, not the tip count, so the legend "(N)" is unchanged). */
    private static void countTipsPerTaxon( final Map<PhylogenyNode, RankTaxon> assignment,
                                           final Map<String, Integer> counts_out ) {
        if ( counts_out == null ) {
            return;
        }
        for( final Entry<PhylogenyNode, RankTaxon> e : assignment.entrySet() ) {
            if ( e.getKey().isExternal() ) {
                counts_out.merge( e.getValue().getName(), 1, Integer::sum );
            }
        }
    }

    /**
     * The clade bands for annotating {@code tree} at {@code rank} with shaded boxes or right-edge bars:
     * one {@link CladeBand} (taxon + distinct color + clade-root) per maximal-monophyletic clade, from
     * the SAME assignment the rank colorizer uses (so paraphyletic groups yield several same-colored
     * bands). Network-pure (cache-only via {@code service}); the band geometry is computed later, at
     * paint time, from each clade's tip coordinates. Unit-testable with an in-memory service.
     */
    final static List<CladeBand> cladeBands( final Phylogeny tree,
                                             final String rank,
                                             final TaxonomicLineageService service ) {
        return cladeBands( tree, rank, service, null );
    }

    /** {@code overrides} (taxon -&gt; user-chosen color) replaces the auto-assigned color for those taxa. */
    final static List<CladeBand> cladeBands( final Phylogeny tree,
                                             final String rank,
                                             final TaxonomicLineageService service,
                                             final Map<String, Color> overrides ) {
        return cladeBands( tree, rank, service, overrides, null );
    }

    /**
     * {@code counts_out}, when non-null, receives each taxon's tip count (tips assigned to it at {@code rank}),
     * so the clade-band legend can show "(N)" and sort by count.
     */
    final static List<CladeBand> cladeBands( final Phylogeny tree,
                                             final String rank,
                                             final TaxonomicLineageService service,
                                             final Map<String, Color> overrides,
                                             final Map<String, Integer> counts_out ) {
        final List<CladeBand> bands = new ArrayList<CladeBand>();
        if ( ( tree == null ) || tree.isEmpty() || ForesterUtil.isEmpty( rank ) ) {
            return bands;
        }
        final Map<PhylogenyNode, RankTaxon> assignment = assignNodesToRankTaxon( tree, rank, service );
        final Map<PhylogenyNode, RankTaxon> roots = maximalMonochromaticRoots( tree, assignment );
        final Map<String, Color> colors = AptxUtil.assignDistinctColors( legendTaxa( assignment, roots ) );
        applyColorOverrides( colors, overrides );
        for( final Entry<PhylogenyNode, RankTaxon> e : roots.entrySet() ) {
            final String taxon = e.getValue().getName();
            final Color c = colors.get( taxon );
            if ( c != null ) {
                bands.add( new CladeBand( taxon, c, e.getKey() ) );
            }
        }
        countTipsPerTaxon( assignment, counts_out );
        return bands;
    }

    /**
     * WRITES each maximal-monochromatic clade's taxon at {@code rank} onto that clade root's internal {@code
     * <taxonomy>} (scientific name + rank + NCBI tax-id) -- a real tree enrichment that round-trips to phyloXML,
     * the persistent counterpart of the display-only {@link #cladeBands}. Uses the SAME assignment + grouping (Spine
     * B), so it annotates exactly the clades the bands would draw. Only INTERNAL clade roots are written (a tip that
     * is its own clade keeps its more specific own taxonomy); a root that already carries a taxonomy is left in place
     * unless {@code overwrite}. Network-pure (cache-only via {@code service}). Returns the number of nodes written.
     */
    final static int writeCladeTaxonomies( final Phylogeny tree, final String rank,
                                           final TaxonomicLineageService service, final boolean overwrite ) {
        if ( ( tree == null ) || tree.isEmpty() || ForesterUtil.isEmpty( rank ) ) {
            return 0;
        }
        final Map<PhylogenyNode, RankTaxon> assignment = assignNodesToRankTaxon( tree, rank, service );
        int written = 0;
        for( final Entry<PhylogenyNode, RankTaxon> e : maximalMonochromaticRoots( tree, assignment ).entrySet() ) {
            final PhylogenyNode root = e.getKey();
            if ( root.isExternal() ) {
                continue; // a single-tip clade keeps its own (more specific) taxonomy -- never downgrade a tip
            }
            if ( root.getNodeData().isHasTaxonomy() && !overwrite ) {
                continue; // preserve a hand-curated / previously-inferred internal taxon unless overwriting
            }
            final RankTaxon rt = e.getValue();
            root.getNodeData().setTaxonomy( TaxonomyUtil.buildNcbiTaxonomy( rt.getName(), rank, rt.getId() ) );
            written++;
        }
        return written;
    }

    /** The provenance sentence for "Annotate Clades by Rank -&gt; write into the tree" (pure/testable; the caller
     *  APPENDS it to the description, never overwrites). */
    final static String cladeTaxonomyProvenance( final Phylogeny phy, final String rank, final int count,
                                                 final boolean overwrite ) {
        final String name = ForesterUtil.isEmpty( phy.getName() ) ? "" : phy.getName();
        return "Used annotate-clades-by-rank (rank " + rank + ( overwrite ? ", overwriting existing internal taxa" : "" )
                + ") to write a taxonomy onto " + count + " internal clade node" + ( count == 1 ? "" : "s" )
                + " in tree named \"" + name + "\" with " + phy.getNumberOfExternalNodes() + " tips.";
    }

    /** Canonical Linnaean rank order for the internal-taxonomy key (ranks not listed are appended alphabetically). */
    private static final String[] CANONICAL_RANKS = { "domain", "superkingdom", "kingdom", "subkingdom",
            "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "superorder",
            "order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe",
            "genus", "subgenus", "species group", "species", "subspecies", "varietas", "forma" };

    /**
     * The DISTINCT internal-node taxa {@code phy} carries, grouped by rank &rarr; (taxon name &rarr; count of internal
     * nodes carrying it), with the ranks in canonical Linnaean order (unknown ranks appended) and the taxa within a
     * rank by count-desc then name. Only INTERNAL nodes with a rank + name are counted (the tips are the source data,
     * excluded). Feeds the draggable "Internal Taxonomy Key" -- the taxonomic backbone an inference / curation /
     * clade-annotation left on the tree.
     */
    final static LinkedHashMap<String, LinkedHashMap<String, Integer>> internalTaxaByRank( final Phylogeny phy ) {
        final LinkedHashMap<String, LinkedHashMap<String, Integer>> out = new LinkedHashMap<String, LinkedHashMap<String, Integer>>();
        if ( ( phy == null ) || phy.isEmpty() ) {
            return out;
        }
        final Map<String, Map<String, Integer>> by_rank = new HashMap<String, Map<String, Integer>>();
        for( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isExternal() || !n.getNodeData().isHasTaxonomy() ) {
                continue;
            }
            final Taxonomy tax = n.getNodeData().getTaxonomy();
            final String rank = tax.getRank();
            final String name = taxonomyLabel( tax );
            if ( ForesterUtil.isEmpty( rank ) || ForesterUtil.isEmpty( name ) ) {
                continue;
            }
            by_rank.computeIfAbsent( rank.toLowerCase( Locale.ROOT ), k -> new HashMap<String, Integer>() )
                    .merge( name, 1, Integer::sum );
        }
        final List<String> ranks = new ArrayList<String>();
        for( final String r : CANONICAL_RANKS ) {
            if ( by_rank.containsKey( r ) ) {
                ranks.add( r );
            }
        }
        final SortedSet<String> extra = new TreeSet<String>( by_rank.keySet() );
        extra.removeAll( ranks );
        ranks.addAll( extra ); // any non-canonical ranks, alphabetically
        for( final String r : ranks ) {
            out.put( r, sortTaxaByCountThenName( by_rank.get( r ) ) );
        }
        return out;
    }

    /** {@code taxon -> count} ordered by count DESC then name ASC (case-insensitive), into a {@link LinkedHashMap}. */
    private static LinkedHashMap<String, Integer> sortTaxaByCountThenName( final Map<String, Integer> m ) {
        final List<Entry<String, Integer>> es = new ArrayList<Entry<String, Integer>>( m.entrySet() );
        es.sort( ( a, b ) -> {
            final int c = b.getValue().compareTo( a.getValue() );
            return ( c != 0 ) ? c : a.getKey().compareToIgnoreCase( b.getKey() );
        } );
        final LinkedHashMap<String, Integer> out = new LinkedHashMap<String, Integer>();
        for( final Entry<String, Integer> e : es ) {
            out.put( e.getKey(), e.getValue() );
        }
        return out;
    }

    /** Replaces the auto-assigned color with the user's override for each taxon that has one. */
    private static void applyColorOverrides( final Map<String, Color> colors, final Map<String, Color> overrides ) {
        if ( ( colors == null ) || ( overrides == null ) || overrides.isEmpty() ) {
            return;
        }
        for( final String taxon : colors.keySet() ) {
            final Color o = overrides.get( taxon );
            if ( o != null ) {
                colors.put( taxon, o ); // value-only update of an existing key is safe during keySet iteration
            }
        }
    }

    /**
     * Every node -- tip AND internal -- mapped to its OWN taxon at {@code rank} (absent when it has none): the one
     * node&rarr;taxon assignment (Spine B) that both the branch colorizer and the clade bands read. A node's own
     * taxon is resolved, in order: (1) its own {@code <taxonomy>} annotated at exactly {@code rank}; else (2) its
     * own taxon resolved to a lineage via the cache-only {@code service} (looked up by its NCBI tax-id and/or name);
     * else -- for a TIP only -- (3) the nearest PROPER-ancestor annotation at
     * {@code rank}. An INTERNAL node stops at (2): it contributes its OWN identity (so a curated / inferred internal
     * taxon becomes visible to the visualization) but never an ancestor's.
     * <p>
     * Values are keyed by tax-id where available ({@link RankTaxon}), so HOMONYMS do not merge; a closing pass
     * ({@link #canonicalizeTaxonIds}) rewrites an id-less entry to a same-named id-bearing one, so the SAME taxon
     * arriving from mixed sources (some with an id, some without) stays ONE group. Network-pure (cache-only).
     * <p>
     * "Tip identity wins" still holds: this per-node own taxon feeds {@link #maximalMonochromaticRoots}, which lets a
     * clade's resolvable tips win over a conflicting ancestor annotation and uses an internal annotation only to fill
     * gaps. The escape hatch to override a bad DB hit is to annotate the node itself at {@code rank}.
     */
    final static Map<PhylogenyNode, RankTaxon> assignNodesToRankTaxon( final Phylogeny tree,
                                                                       final String rank,
                                                                       final TaxonomicLineageService service ) {
        final Map<PhylogenyNode, RankTaxon> assignment = new HashMap<PhylogenyNode, RankTaxon>();
        if ( ( tree == null ) || tree.isEmpty() || ForesterUtil.isEmpty( rank ) ) {
            return assignment;
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            final RankTaxon rt = node.isExternal() ? tipRankTaxon( node, rank, service )
                                                   : ownRankTaxon( node, rank, service );
            if ( rt != null ) {
                assignment.put( node, rt );
            }
        }
        canonicalizeTaxonIds( assignment );
        return assignment;
    }

    /** Backward-compatible tip-only, NAME-valued view of {@link #assignNodesToRankTaxon} (each {@link RankTaxon}
     *  projected to its display name) -- the historical signature the tests read. */
    final static Map<PhylogenyNode, String> assignTipsToRankTaxon( final Phylogeny tree,
                                                                   final String rank,
                                                                   final TaxonomicLineageService service ) {
        final Map<PhylogenyNode, String> out = new HashMap<PhylogenyNode, String>();
        for( final Entry<PhylogenyNode, RankTaxon> e : assignNodesToRankTaxon( tree, rank, service ).entrySet() ) {
            if ( e.getKey().isExternal() ) {
                out.put( e.getKey(), e.getValue().getName() );
            }
        }
        return out;
    }

    /** A node's OWN taxon at {@code rank}: (1) its own {@code <taxonomy>} annotated at exactly {@code rank}; else
     *  (2) its own taxon resolved to a cached lineage, tax-id first then name. No ancestor / node-name fallback. */
    final static RankTaxon ownRankTaxon( final PhylogenyNode node, final String rank,
                                         final TaxonomicLineageService service ) {
        // (1) own taxonomy annotated at exactly this rank -- the most specific identity, no network
        final RankTaxon self = selfRankTaxonRT( node, rank );
        if ( self != null ) {
            return self;
        }
        // (2) resolve the own taxon's lineage (cache-only): try the NCBI tax-id AND the own name -- either may be the
        // key the cache was primed under (the colorize fetch flow primes by NAME via unresolvedTipTaxa/tipQueryName,
        // while a Fetch-tool run primes by tax-id), so an id-only lookup would miss a name-primed cache
        if ( ( service != null ) && node.getNodeData().isHasTaxonomy() ) {
            final Taxonomy tax = node.getNodeData().getTaxonomy();
            return resolveRankViaService( service, rank, ncbiId( tax ), ownTaxonName( tax ) );
        }
        return null;
    }

    /** Reads {@code rank} off the lineage the cache-only {@code service} holds for a taxon, trying the NCBI tax-id
     *  first and then the name (either may be the key the cache was primed under), or null. The {@link RankTaxon}
     *  carries the per-rank tax-id ({@link TaxonLineage#taxIdAt}) whichever key resolved it. */
    private static RankTaxon resolveRankViaService( final TaxonomicLineageService service, final String rank,
                                                    final String id, final String name ) {
        if ( service == null ) {
            return null;
        }
        TaxonLineage rl = null;
        if ( !ForesterUtil.isEmpty( id ) ) {
            rl = service.lineageOf( id );
        }
        if ( ( rl == null ) && !ForesterUtil.isEmpty( name ) ) {
            rl = service.lineageOf( name );
        }
        if ( rl != null ) {
            final String at = rl.at( rank );
            if ( !ForesterUtil.isEmpty( at ) ) {
                return new RankTaxon( rl.taxIdAt( rank ), at );
            }
        }
        return null;
    }

    /** A TIP's taxon at {@code rank}: {@link #ownRankTaxon}, else (2b) resolve by the tip's node NAME (a bare-named
     *  tip with no {@code <taxonomy>}; {@link #tipQueryName} also walks to an ancestor scientific name -- the
     *  long-standing tip behavior), else (3) the nearest PROPER-ancestor annotation at {@code rank}. */
    private final static RankTaxon tipRankTaxon( final PhylogenyNode tip, final String rank,
                                                 final TaxonomicLineageService service ) {
        final RankTaxon own = ownRankTaxon( tip, rank, service );
        if ( own != null ) {
            return own;
        }
        // (2b) resolve by the tip's query name -- its own name, else the NODE name, else an ancestor scientific name
        // (tipQueryName): the long-standing tip fallback. NOT gated on isHasTaxonomy, so a tip carrying an effectively
        // nameless taxonomy (rank-only / a non-NCBI id) is still resolved by its node name, as before Spine B.
        final RankTaxon by_name = resolveRankViaService( service, rank, null, tipQueryName( tip ) );
        if ( by_name != null ) {
            return by_name;
        }
        // (3) nearest ancestor annotation at rank
        return ancestorRankTaxonRT( tip, rank );
    }

    /** The queryable NAME of a taxon's OWN identity (scientific &rarr; code &rarr; common), matching
     *  {@link #tipQueryName}'s own-identity order; null when it carries none. */
    private final static String ownTaxonName( final Taxonomy tax ) {
        if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            return tax.getScientificName();
        }
        if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
            return tax.getTaxonomyCode();
        }
        if ( !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
            return tax.getCommonName();
        }
        return null;
    }

    /** Canonicalizes an assignment so the SAME taxon from mixed sources is one group: when a taxon NAME appears
     *  WITH an NCBI id in some entry and id-less in another, the id-less entries adopt that id. A clean homonym
     *  (one name, two DIFFERENT ids) is left split -- the whole point of tax-id keying. */
    private static void canonicalizeTaxonIds( final Map<PhylogenyNode, RankTaxon> assignment ) {
        final Map<String, String> name_to_id = new HashMap<String, String>();
        for( final RankTaxon rt : assignment.values() ) {
            if ( rt.getId() != null ) {
                name_to_id.putIfAbsent( rt.getName().toLowerCase( Locale.ROOT ), rt.getId() );
            }
        }
        if ( name_to_id.isEmpty() ) {
            return;
        }
        for( final Entry<PhylogenyNode, RankTaxon> e : assignment.entrySet() ) {
            final RankTaxon rt = e.getValue();
            if ( rt.getId() == null ) {
                final String id = name_to_id.get( rt.getName().toLowerCase( Locale.ROOT ) );
                if ( id != null ) {
                    e.setValue( new RankTaxon( id, rt.getName() ) );
                }
            }
        }
    }

    /** The node's OWN {@code <taxonomy>} as a {@link RankTaxon} (name + NCBI id) iff its rank equals {@code rank}
     *  (case-insensitive) and the label is non-empty, else null. */
    final static RankTaxon selfRankTaxonRT( final PhylogenyNode node, final String rank ) {
        if ( node.getNodeData().isHasTaxonomy() ) {
            final Taxonomy tax = node.getNodeData().getTaxonomy();
            if ( !ForesterUtil.isEmpty( tax.getRank() ) && tax.getRank().equalsIgnoreCase( rank ) ) {
                final String label = taxonomyLabel( tax );
                if ( !ForesterUtil.isEmpty( label ) ) {
                    return new RankTaxon( ncbiId( tax ), label );
                }
            }
        }
        return null;
    }

    /** The taxon label of {@code node}'s OWN taxonomy iff its rank equals {@code rank}, else null -- used by
     *  {@link #unresolvedTipTaxa} as the "the tip resolves its own identity in-tree" test. */
    final static String selfRankTaxon( final PhylogenyNode node, final String rank ) {
        final RankTaxon rt = selfRankTaxonRT( node, rank );
        return ( rt == null ) ? null : rt.getName();
    }

    /** The nearest PROPER ancestor's own {@code <taxonomy>} at exactly {@code rank} as a {@link RankTaxon}, or null
     *  -- the tip fallback in {@link #tipRankTaxon} for a tip that cannot resolve its own identity. */
    private final static RankTaxon ancestorRankTaxonRT( final PhylogenyNode tip, final String rank ) {
        for( PhylogenyNode n = tip.getParent(); n != null; n = n.getParent() ) {
            final RankTaxon rt = selfRankTaxonRT( n, rank );
            if ( rt != null ) {
                return rt;
            }
        }
        return null;
    }

    /**
     * The most specific name to query a taxonomy DB with for {@code tip}. The tip's OWN identity is
     * the most specific, so it is preferred (scientific name, else code, else common name, else node
     * name); only when the tip carries no identity at all do we fall back to the nearest ancestor's
     * scientific name (which can still place the tip at a rank at/above that ancestor). Querying an
     * ancestor's name before the tip's own code/common name would lose specificity and mis-resolve.
     */
    final static String tipQueryName( final PhylogenyNode tip ) {
        if ( tip.getNodeData().isHasTaxonomy() ) {
            final Taxonomy tax = tip.getNodeData().getTaxonomy();
            if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                return tax.getScientificName();
            }
            if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                return tax.getTaxonomyCode();
            }
            if ( !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
                return tax.getCommonName();
            }
        }
        if ( !ForesterUtil.isEmpty( tip.getName() ) ) {
            return tip.getName();
        }
        for( PhylogenyNode n = tip.getParent(); n != null; n = n.getParent() ) {
            if ( n.getNodeData().isHasTaxonomy()
                    && !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                return n.getNodeData().getTaxonomy().getScientificName();
            }
        }
        return null;
    }

    /**
     * The distinct taxon query-names of tips that cannot resolve their OWN identity at {@code rank} yet -- i.e.
     * exactly the names a caller must {@link TaxonomicLineageService#fetch} (off the EDT) to place more tips.
     * A tip is excluded only when it self-resolves (its own taxonomy is annotated at {@code rank}) or its
     * query-name is already in the cache; an ANCESTOR annotation does NOT suppress the fetch, because under
     * "tip identity wins" the tip's own DB resolution must be available to override a wrong/partial ancestor.
     * Cache hits are excluded even when the cache lacks {@code rank} (refetching would not help), so a second
     * call after a fetch pass returns an empty set (no repeated prompts).
     */
    final static SortedSet<String> unresolvedTipTaxa( final Phylogeny tree,
                                                      final String rank,
                                                      final TaxonomicLineageService service ) {
        final SortedSet<String> names = new TreeSet<String>();
        if ( ( tree == null ) || tree.isEmpty() || ForesterUtil.isEmpty( rank ) ) {
            return names;
        }
        for( final PhylogenyNodeIterator it = tree.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode tip = it.next();
            if ( selfRankTaxon( tip, rank ) != null ) {
                continue; // the tip resolves its OWN identity in-tree; an ancestor annotation does NOT count here
            }
            final String q = tipQueryName( tip );
            if ( ForesterUtil.isEmpty( q ) ) {
                continue;
            }
            if ( ( service != null ) && ( service.lineageOf( q ) != null ) ) {
                continue; // already attempted/cached -- refetching would not help
            }
            names.add( q );
        }
        return names;
    }

    /**
     * Each external node mapped to its resolved {@link TaxonLineage} for ancestral-taxonomy inference: the tip's
     * cached lineage from {@code service} when present, else a minimal lineage reconstructed from the tip's own
     * stored {@link Taxonomy} (its scientific name / rank / NCBI id as the deepest level, plus any in-memory
     * {@code getLineage()} names as ancestors). Tips with nothing usable map to {@link TaxonLineage#EMPTY}.
     * Pure -- reads the tree + the lineage cache, no network.
     */
    final static Map<PhylogenyNode, TaxonLineage> tipLineages( final Phylogeny tree,
                                                               final TaxonomicLineageService service ) {
        final Map<PhylogenyNode, TaxonLineage> out = new HashMap<PhylogenyNode, TaxonLineage>();
        if ( ( tree == null ) || tree.isEmpty() ) {
            return out;
        }
        for( final PhylogenyNodeIterator it = tree.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode tip = it.next();
            TaxonLineage tl = null;
            if ( service != null ) {
                final String q = tipQueryName( tip );
                if ( !ForesterUtil.isEmpty( q ) ) {
                    tl = service.lineageOf( q );
                }
            }
            if ( ( tl == null ) || tl.isEmpty() ) {
                tl = lineageFromStoredTaxonomy( tip );
            }
            out.put( tip, ( tl == null ) ? TaxonLineage.EMPTY : tl );
        }
        return out;
    }

    /** Reconstruct a minimal {@link TaxonLineage} from a tip's own stored {@link Taxonomy}. Per-ancestor ranks/ids
     *  do not survive a phyloXML round-trip, so the ancestors (from the in-memory {@code getLineage()} names) carry
     *  names only; the own taxon keeps its rank + NCBI id. {@link TaxonLineage#EMPTY} when nothing is usable. */
    private final static TaxonLineage lineageFromStoredTaxonomy( final PhylogenyNode tip ) {
        if ( !tip.getNodeData().isHasTaxonomy() ) {
            return TaxonLineage.EMPTY;
        }
        final Taxonomy tax = tip.getNodeData().getTaxonomy();
        final String own = taxonomyLabel( tax );
        if ( ForesterUtil.isEmpty( own ) ) {
            return TaxonLineage.EMPTY;
        }
        final List<TaxonLineage.Ancestor> anc = new ArrayList<TaxonLineage.Ancestor>();
        if ( tax.getLineage() != null ) {
            for( final String name : tax.getLineage() ) {
                if ( !ForesterUtil.isEmpty( name ) && !name.equalsIgnoreCase( own ) ) {
                    anc.add( new TaxonLineage.Ancestor( name, null, null ) );
                }
            }
        }
        return new TaxonLineage( ncbiId( tax ), tax.getRank(), own, null, anc );
    }

    /** The NCBI tax-id from a taxonomy's identifier (only the "ncbi" provider), or null. */
    private final static String ncbiId( final Taxonomy tax ) {
        if ( ( tax.getIdentifier() != null ) && !ForesterUtil.isEmpty( tax.getIdentifier().getValue() )
                && "ncbi".equalsIgnoreCase( tax.getIdentifier().getProvider() ) ) {
            return tax.getIdentifier().getValue();
        }
        return null;
    }

    /**
     * The distinct query-names of tips that have NO usable lineage yet for ancestral-taxonomy inference -- i.e.
     * exactly the names a caller must {@link TaxonomicLineageService#fetch} to place more internal nodes. A tip is
     * excluded when it already carries a multi-level lineage in-tree (in-memory {@code getLineage()} names) or when
     * its query-name is cached (refetching would not help), so a second call after a fetch pass returns an empty
     * set (no repeated prompts). Rank-agnostic sibling of {@link #unresolvedTipTaxa}.
     */
    final static SortedSet<String> tipsWithoutLineage( final Phylogeny tree,
                                                       final TaxonomicLineageService service ) {
        final SortedSet<String> names = new TreeSet<String>();
        if ( ( tree == null ) || tree.isEmpty() ) {
            return names;
        }
        for( final PhylogenyNodeIterator it = tree.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode tip = it.next();
            if ( hasStoredLineage( tip ) ) {
                continue; // the tip already carries an ancestor lineage in-tree -- no fetch needed
            }
            final String q = tipQueryName( tip );
            if ( ForesterUtil.isEmpty( q ) ) {
                continue;
            }
            if ( ( service != null ) && ( service.lineageOf( q ) != null ) ) {
                continue; // already attempted/cached
            }
            names.add( q );
        }
        return names;
    }

    /** True when a tip carries a usable ANCESTOR lineage in the tree itself: an in-memory {@code getLineage()} name
     *  other than the tip's own taxon. An own taxon alone is NOT enough -- inference needs ancestors to intersect,
     *  so it must fetch. Mirrors {@link #lineageFromStoredTaxonomy}, which strips own-name entries, so the two agree
     *  on what counts as an ancestor (else an own-only lineage would be neither fetched nor usable). */
    private final static boolean hasStoredLineage( final PhylogenyNode tip ) {
        if ( !tip.getNodeData().isHasTaxonomy() ) {
            return false;
        }
        final Taxonomy tax = tip.getNodeData().getTaxonomy();
        final List<String> lineage = tax.getLineage();
        if ( ( lineage == null ) || lineage.isEmpty() ) {
            return false;
        }
        final String own = taxonomyLabel( tax );
        for( final String name : lineage ) {
            if ( !ForesterUtil.isEmpty( name ) && !name.equalsIgnoreCase( own ) ) {
                return true; // a real ancestor (not just the tip's own taxon)
            }
        }
        return false;
    }

    /**
     * Each node that roots a <i>maximal</i> clade resolving to one rank taxon, mapped to that taxon (Spine B:
     * tax-id-aware, respecting internal annotations as gap-fillers). Per node the subtree taxon is reconciled from
     * the children and the node's OWN {@code assignment} entry:
     * <ul>
     * <li>children genuinely CONFLICT (&ge;2 distinct placed child taxa, or a mixed child) &rarr; mixed (never
     *     filled; the node's own annotation is ignored here);</li>
     * <li>the placed children AGREE on a taxon T &rarr; T (so a clade's resolvable tips WIN over a differing
     *     ancestor annotation -- "tip identity wins"); an UNPLACED sibling is swept into T only when the node's own
     *     annotation authorizes that same T (else the unplaced tip forces mixed, as before);</li>
     * <li>ALL children unplaced &rarr; the node's OWN annotation fills the gap (defines the clade), else unplaced.</li>
     * </ul>
     * Then the maximal (topmost) node of each taxon is returned. Paraphyly yields several same-taxon roots.
     * Equality is tax-id-aware ({@link RankTaxon}), so homonyms are not merged into one clade.
     */
    final static Map<PhylogenyNode, RankTaxon> maximalMonochromaticRoots( final Phylogeny tree,
                                                                          final Map<PhylogenyNode, RankTaxon> assignment ) {
        final Map<PhylogenyNode, RankTaxon> subtree = new HashMap<PhylogenyNode, RankTaxon>();
        final Map<PhylogenyNode, RankTaxon> roots = new LinkedHashMap<PhylogenyNode, RankTaxon>();
        if ( ( tree == null ) || tree.isEmpty() ) {
            return roots;
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final RankTaxon own = assignment.get( n ); // the node's OWN taxon (null when it has none)
            if ( n.isExternal() ) {
                subtree.put( n, own ); // a placed tip, or null (unplaced); the leaf never "conflicts"
                continue;
            }
            RankTaxon consensus = null;   // the single taxon the PLACED children agree on
            boolean conflict = false;     // >= 2 distinct placed child taxa, or a mixed child
            boolean any_unplaced = false; // a child whose subtree pins no taxon
            for( final PhylogenyNode c : n.getDescendants() ) {
                final RankTaxon cs = subtree.get( c );
                if ( cs == MIXED ) {
                    conflict = true;
                    break;
                }
                if ( cs == null ) {
                    any_unplaced = true;
                }
                else if ( consensus == null ) {
                    consensus = cs;
                }
                else if ( !consensus.equals( cs ) ) {
                    conflict = true;
                    break;
                }
            }
            final RankTaxon result;
            if ( conflict ) {
                result = MIXED;
            }
            else if ( consensus != null ) {
                // placed children agree on `consensus`; an unplaced sibling is swept in only if the node's own
                // annotation names that same taxon (authorizing the fill), else the unplaced tip keeps it mixed
                result = ( !any_unplaced || consensus.equals( own ) ) ? consensus : MIXED;
            }
            else {
                result = own; // all children unplaced -> the node's own annotation fills the gap (may be null)
            }
            subtree.put( n, result );
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final RankTaxon t = subtree.get( n );
            if ( ( t != null ) && ( t != MIXED ) ) {
                final PhylogenyNode p = n.getParent();
                final RankTaxon pt = ( p == null ) ? null : subtree.get( p );
                if ( ( pt == null ) || ( pt == MIXED ) || !t.equals( pt ) ) {
                    roots.put( n, t );
                }
            }
        }
        return roots;
    }

    final static boolean isHasAssignedEvent( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasEvent() ) {
            return false;
        }
        if ( ( node.getNodeData().getEvent() ).isUnassigned() ) {
            return false;
        }
        return true;
    }

    final static boolean isSequenceEmpty( final Sequence seq ) {
        return ( seq.getAccession() == null ) && ForesterUtil.isEmpty( seq.getName() )
                && ForesterUtil.isEmpty( seq.getGeneName() ) && ForesterUtil.isEmpty( seq.getSymbol() );
    }

    final static boolean isTaxonomyEmpty( final Taxonomy tax ) {
        return ( ( tax.getIdentifier() == null ) && ForesterUtil.isEmpty( tax.getTaxonomyCode() )
                && ForesterUtil.isEmpty( tax.getCommonName() ) && ForesterUtil.isEmpty( tax.getScientificName() )
                && tax.getSynonyms().isEmpty() );
    }

    final static Phylogeny subTree( final PhylogenyNode new_root, final Phylogeny source_phy ) {
        final Phylogeny new_phy = new Phylogeny();
        new_phy.setRooted( true );
        new_phy.setName( source_phy.getName() );
        new_phy.setDescription( source_phy.getDescription() );
        new_phy.setType( source_phy.getType() );
        new_phy.setDistanceUnit( source_phy.getDistanceUnit() );
        new_phy.setConfidence( source_phy.getConfidence() );
        new_phy.setIdentifier( source_phy.getIdentifier() );
        new_phy.setRoot( new_root.copyNodeDataShallow() );
        int i = 0;
        for( final PhylogenyNode n : new_root.getDescendants() ) {
            new_phy.getRoot().setChildNode( i++, n );
        }
        return new_phy;
    }

    /**
     * The minimum vertical leaf-to-leaf spacing -- expressed as a y-distance -- at which leaf
     * labels of the given pixel height stop overlapping. Adjacent leaf rows are spaced
     * {@code 2 * y-distance} apart (see {@code TreePanel.resetPreferredSize} /
     * {@code calcParametersForPainting}), so labels of height {@code h} no longer overlap once
     * {@code 2 * y-distance >= h}, i.e. {@code y-distance >= h / 2}. A small margin is added for
     * breathing room; it also keeps the dynamic-hiding factor
     * ({@code round( h / (1.5 * y-distance) )}, see {@code TreePanel.calcDynamicHidingFactor}) at
     * {@code <= 1}, so the "Dyna Hide" indicator clears.
     */
    final static float yDistanceToAvoidLabelOverlap( final int label_height_px ) {
        return ( label_height_px / 2.0f ) * 1.1f;
    }

    /**
     * The logical-&gt;screen rotation for a vertical (root-top / root-bottom) tree orientation, given the tree's
     * LOGICAL width {@code w} (depth extent, root x=0 .. tip x=w) and height {@code h} (breadth/tip-spread extent,
     * y=0 .. y=h). Pure rotations (determinant +1, so nothing mirrors); the translate keeps the rotated tree in the
     * positive quadrant:
     * <ul>
     *   <li>ROOT_TOP turns the page 90&deg; clockwise: {@code (x,y) -> (h - y, x)} (root to the top, tips to bottom)</li>
     *   <li>ROOT_BOTTOM turns it 90&deg; counter-clockwise: {@code (x,y) -> (y, w - x)} (root to the bottom)</li>
     * </ul>
     * Any other value (ROOT_LEFT) returns the identity. Pure math, no toolkit -&gt; headless-testable.
     */
    static AffineTransform orientationTransformFor( final Options.TREE_ORIENTATION orientation, final double w,
                                                    final double h ) {
        if ( orientation == Options.TREE_ORIENTATION.ROOT_TOP ) {
            final AffineTransform r = AffineTransform.getTranslateInstance( h, 0 );
            r.rotate( Math.PI / 2.0 );
            return r;
        }
        if ( orientation == Options.TREE_ORIENTATION.ROOT_BOTTOM ) {
            final AffineTransform r = AffineTransform.getTranslateInstance( 0, w );
            r.rotate( -Math.PI / 2.0 );
            return r;
        }
        return new AffineTransform(); // identity for ROOT_LEFT
    }

    /**
     * The width to right-align a vertical-orientation internal-node label on. The label is drawn as an (optional)
     * taxonomy segment followed by the node-data segment; the taxonomy label always ends with a trailing part-separator
     * space, so when it is the RIGHTMOST drawn element (no node data follows) that trailing space is excluded from the
     * alignment width -- otherwise the visible label would right-align one space-width left of the branch. Clamped to
     * {@code >= 0}.
     *
     * @param tax_w        measured width of the taxonomy segment (includes its trailing space), 0 if none
     * @param data_w       measured width of the node-data segment, 0 if none
     * @param space_width  width of a space in the label font
     * @param show_tax     whether a taxonomy segment is drawn
     * @param data_empty   whether the node-data segment is empty
     */
    static int internalLabelAlignWidth( final int tax_w, final int data_w, final int space_width,
                                        final boolean show_tax, final boolean data_empty ) {
        final int total = tax_w + data_w;
        if ( show_tax && data_empty ) {
            return Math.max( 0, total - space_width );
        }
        return total;
    }

    /**
     * Auto-pick the tip-label angle for a vertical dendrogram from how much room each tip has: upright (0°) when the
     * longest label fits between adjacent tips, else diagonal (45°) while its horizontal projection fits, else
     * vertical (90°). {@code tip_spacing} is the distance between two adjacent tips along the breadth axis (i.e.
     * {@code 2 * getYdistance()}); {@code longest_label_width} the widest tip label. A degenerate/absent layout
     * (non-positive inputs) falls back to vertical (always fits). A centred horizontal label reaches ±width/2 per
     * side, so it needs width ≤ spacing; a 45° label's horizontal footprint is ~width·cos45 = width/√2.
     */
    static Options.TIP_LABEL_DIRECTION autoTipLabelDirection( final double tip_spacing,
                                                              final double longest_label_width ) {
        if ( ( tip_spacing <= 0 ) || ( longest_label_width <= 0 ) ) {
            return Options.TIP_LABEL_DIRECTION.VERTICAL;
        }
        if ( longest_label_width <= tip_spacing ) {
            return Options.TIP_LABEL_DIRECTION.HORIZONTAL;
        }
        if ( ( longest_label_width * 0.70710678 ) <= tip_spacing ) {
            return Options.TIP_LABEL_DIRECTION.ANGLED;
        }
        return Options.TIP_LABEL_DIRECTION.VERTICAL;
    }

    /** ref-namespace prefix for INTERNAL Aptx metadata properties (e.g. the persisted Re-import annotation profile on
     *  the root). These are machinery for the save/reload round-trip, not user content, so they are hidden from the
     *  user-facing node-data displays (rollover popup, Display Node Data) by {@link #userVisiblePropertiesText}. */
    final static String INTERNAL_PROPERTY_REF_PREFIX = "aptx:";

    static boolean isInternalPropertyRef( final String ref ) {
        return ( ref != null ) && ref.startsWith( INTERNAL_PROPERTY_REF_PREFIX );
    }

    /** The property list as newline-joined display text, EXCLUDING internal {@code aptx:*} metadata -- mirrors
     *  {@link PropertiesList#asSimpleText()} but drops the properties an end user should never see. */
    static StringBuffer userVisiblePropertiesText( final PropertiesList props ) {
        final StringBuffer sb = new StringBuffer();
        if ( props != null ) {
            for( final Property p : props.getProperties() ) {
                if ( isInternalPropertyRef( p.getRef() ) || isVisualStylePropertyRef( p.getRef() ) ) {
                    continue;
                }
                if ( sb.length() > 0 ) {
                    sb.append( "\n" );
                }
                sb.append( p.asText() );
            }
        }
        return sb;
    }

    /** A {@code style:*} property ref -- Archaeopteryx's own node-visual-style metadata. A NODE-level one is
     *  consumed into {@link NodeVisualData} at parse time and never reaches the property list, but a stray one
     *  with a different {@code applies_to} can, and it is machinery rather than user content either way. */
    static boolean isVisualStylePropertyRef( final String ref ) {
        return ( ref != null ) && ref.startsWith( NodeVisualData.APTX_VISUALIZATION_REF );
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Node properties in the tip label
    // ---------------------------------------------------------------------------------------------------------------
    // A node easily carries ten properties, and the "Properties" display option used to append EVERY one of them to
    // the label as "<full ref>: <value>" -- NEWLINE-joined, into a string that is drawn as a single line. That gave an
    // unreadably long label ("seq1 nextstrain:country: Brazil" and worse) with literal newlines embedded in it.
    // The label now shows VALUES ONLY, comma-joined, on one line; WHICH fields appear and in WHICH order is chosen in
    // Tools > Annotation Fields (see TreePanel#setLabelPropertyRefs). The full "ref: value" list is still one hover
    // away in the rollover popup, and in Display Node Data -- so narrowing the label loses nothing.

    /** Separator between property values in a tip label. A comma reads as a list, keeping the attributes visibly
     *  distinct from the identity part of the label (name / accession), which is space-joined. */
    private static final String LABEL_PROPERTY_SEPARATOR = ", ";

    /**
     * The node's properties as ONE-LINE tip-label text: values only (no {@code ref:} prefix), comma-joined.
     * <p>
     * {@code refs_in_order} selects and ORDERS the fields; a ref the node does not carry is simply skipped. Pass
     * {@code null} for the default -- every user-visible property, in the property list's own (ref-sorted) order,
     * which is what the rollover popup and the node panel show too. Internal {@code aptx:*} / {@code style:*}
     * metadata and empty values never appear. A property's unit, if it has one, follows its value.
     */
    static String labelPropertiesText( final PropertiesList props, final List<String> refs_in_order ) {
        return labelPropertiesText( props, refs_in_order, null );
    }

    /**
     * As {@link #labelPropertiesText(PropertiesList, List)}, but with the refs that are already drawn as annotation
     * COLUMNS excluded from the default selection -- so "a field has one display role" holds from the moment a tree
     * is opened, not only once the user has visited the Annotation Fields chooser. {@code excluded} is ignored when
     * {@code refs_in_order} is given: an explicit choice is the user's, and the chooser already keeps it exclusive.
     */
    static String labelPropertiesText( final PropertiesList props, final List<String> refs_in_order,
                                       final Set<String> excluded ) {
        final StringBuilder sb = new StringBuilder();
        if ( props == null ) {
            return "";
        }
        if ( refs_in_order == null ) {
            for( final Property p : props.getProperties() ) {
                if ( ( excluded != null ) && excluded.contains( p.getRef() ) ) {
                    continue;
                }
                appendLabelPropertyValue( p, sb );
            }
        }
        else {
            // the chosen order wins over the list's own order, so the caller controls how the label reads
            for( final String ref : refs_in_order ) {
                for( final Property p : props.getProperties() ) {
                    if ( ( ref != null ) && ref.equals( p.getRef() ) ) {
                        appendLabelPropertyValue( p, sb );
                    }
                }
            }
        }
        return sb.toString();
    }

    private static void appendLabelPropertyValue( final Property p, final StringBuilder sb ) {
        // isEmptyTrimmed, not isEmpty: a whitespace-only value contributes no text, so emitting the ", " separator
        // for it would leave a doubled separator mid-label ("Brazil, , cat") or a trailing one
        if ( isInternalPropertyRef( p.getRef() ) || isVisualStylePropertyRef( p.getRef() )
                || ForesterUtil.isEmptyTrimmed( p.getValue() ) ) {
            return;
        }
        if ( sb.length() > 0 ) {
            sb.append( LABEL_PROPERTY_SEPARATOR );
        }
        sb.append( p.getValue().trim() );
        if ( !ForesterUtil.isEmpty( p.getUnit() ) ) {
            // a phyloXML unit is namespaced ("METRIC:m"); the label wants the unit, not the namespace
            sb.append( " " );
            sb.append( p.getUnit().substring( p.getUnit().lastIndexOf( ':' ) + 1 ) );
        }
    }

    /**
     * Every user-visible property ref carried by ANY node of the tree -- internal nodes included -- sorted by the
     * name the user actually sees ({@link PropertyColorScheme#displayName(String)}), ties broken by ref.
     * <p>
     * This is the field inventory for the Annotation Fields chooser, and it is deliberately BROADER than
     * {@link PropertyColorScheme#colorableRefs(Phylogeny)}: that one drops constant fields, per-tip-unique
     * categorical fields and internal-node-only fields, because none of the three can COLOR anything -- but all
     * three are perfectly good label text (a study tag, an accession, an ancestral-state annotation). Reusing the
     * colorable list here would silently hide the most label-worthy fields.
     */
    static List<String> userVisiblePropertyRefs( final Phylogeny phy ) {
        final List<String> refs = new ArrayList<String>();
        if ( ( phy == null ) || phy.isEmpty() ) {
            return refs;
        }
        final SortedSet<String> seen = new TreeSet<String>();
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( ( n.getNodeData() == null ) || ( n.getNodeData().getProperties() == null ) ) {
                continue;
            }
            for( final Property p : n.getNodeData().getProperties().getProperties() ) {
                if ( !isInternalPropertyRef( p.getRef() ) && !isVisualStylePropertyRef( p.getRef() ) ) {
                    seen.add( p.getRef() );
                }
            }
        }
        refs.addAll( seen );
        refs.sort( ( a, b ) -> {
            final int c = PropertyColorScheme.displayName( a ).compareToIgnoreCase( PropertyColorScheme.displayName( b ) );
            return ( c != 0 ) ? c : a.compareTo( b );
        } );
        return refs;
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Broken / truncated long branches (display only)
    // ---------------------------------------------------------------------------------------------------------------
    // A branch far longer than the rest (a distant outgroup, a fast-evolving lineage) squashes the informative part of
    // a phylogram to an unreadable sliver. The "Break Long Branches" display option caps how long such a branch is
    // DRAWN (marking it with an axis-break glyph) and re-derives the depth scale from the CAPPED height, so the rest of
    // the tree reclaims the freed width. The tree data is never touched -- the true branch length is still shown as the
    // branch-length label. A branch is "long" when its length exceeds LONG_BRANCH_BREAK_MULTIPLIER times the median of
    // the tree's strictly-positive branch lengths -- a robust reference: the median is unaffected by the one huge
    // outlier we are trying to detect and by the many zero-length branches of a polytomy-heavy tree.

    /** Median of the strictly-positive branch lengths (distance-to-parent) over every node, or 0 if the tree has none
     *  (a cladogram, or an all-zero-length tree). Pure. */
    static double medianPositiveBranchLength( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return 0;
        }
        final List<Double> pos = new ArrayList<>();
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final double dtp = it.next().getDistanceToParent();
            if ( dtp > 0 ) {
                pos.add( dtp );
            }
        }
        if ( pos.isEmpty() ) {
            return 0;
        }
        final double[] a = new double[ pos.size() ];
        for( int i = 0; i < a.length; i++ ) {
            a[ i ] = pos.get( i );
        }
        Arrays.sort( a );
        final int n = a.length;
        return ( ( n % 2 ) == 1 ) ? a[ n / 2 ] : ( ( a[ ( n / 2 ) - 1 ] + a[ n / 2 ] ) / 2.0 );
    }

    /** The model-length cap above which a branch is drawn broken: {@code multiplier * median(positive branch length)},
     *  or 0 when there is no positive branch length to reference (so capping is inactive). Pure. */
    static double longBranchBreakCap( final Phylogeny phy, final double multiplier ) {
        final double median = medianPositiveBranchLength( phy );
        return ( median > 0 ) ? ( multiplier * median ) : 0;
    }

    /** Deepest root-to-tip path length after each branch is capped at {@code cap}, taking display-collapse into account.
     *  Mirrors {@code Phylogeny.calculateSubtreeHeight} (from which the non-break depth scale is derived) EXACTLY -- the
     *  root's own incoming branch is included, and a display-collapsed clade is measured only to its collapsed root (the
     *  height accumulation resets when the walk passes a collapsed node) -- but with {@code min(dtp, cap)} for the length
     *  of each branch. So with a cap larger than every branch this equals the ordinary tree height. Pure; returns 0 for
     *  an empty / branch-length-less tree or a non-positive cap. */
    static double cappedTreeHeight( final Phylogeny phy, final double cap, final boolean take_collapse_into_account ) {
        if ( ( phy == null ) || phy.isEmpty() || ( cap <= 0 ) ) {
            return 0;
        }
        final PhylogenyNode root = phy.getRoot();
        if ( root.isExternal() || ( take_collapse_into_account && root.isCollapse() ) ) {
            return cappedDistance( root, cap );
        }
        double max = 0;
        for( final PhylogenyNode ext : root.getAllExternalDescendants() ) {
            double h = 0;
            PhylogenyNode d = ext;
            while ( d != root ) {
                if ( take_collapse_into_account && d.isCollapse() ) {
                    h = 0; // a collapsed clade is drawn only to its root -> the depth below it does not count
                }
                h += cappedDistance( d, cap );
                d = d.getParent();
            }
            h += cappedDistance( root, cap ); // include the (capped) root branch, like calculateHeight
            if ( h > max ) {
                max = h;
            }
        }
        return max;
    }

    /** {@code cappedTreeHeight} ignoring collapse (take_collapse_into_account = true is the usual non-average case). */
    static double cappedTreeHeight( final Phylogeny phy, final double cap ) {
        return cappedTreeHeight( phy, cap, true );
    }

    /** Sum of the CAPPED branch lengths ({@code min(dtp, cap)}) from {@code node} up to (not including) the root -- the
     *  capped analogue of {@code PhylogenyNode.calculateDistanceToRoot}, used to place a node's radius/spoke in the
     *  radial (circular / unrooted) layouts while Break Long Branches is active. Pure. */
    static double cappedDistanceToRoot( final PhylogenyNode node, final double cap ) {
        double d = 0;
        PhylogenyNode n = node;
        while ( ( n != null ) && !n.isRoot() ) {
            final double dtp = n.getDistanceToParent();
            if ( dtp > 0 ) {
                d += ( dtp > cap ) ? cap : dtp;
            }
            n = n.getParent();
        }
        return d;
    }

    /** The largest {@link #cappedDistanceToRoot} over all external tips (root branch EXCLUDED) -- the capped analogue of
     *  {@code PhylogenyMethods.calculateMaxDistanceToRoot}, used as the RADIAL (circular / unrooted) normalizer so the
     *  deepest capped tip lands on the outer ring/diameter exactly as the uncapped deepest tip does. Differs from
     *  {@link #cappedTreeHeight}, which INCLUDES the capped root branch (parity with calculateHeight -- the RECTANGULAR
     *  depth) and is collapse-aware; this mirrors getMaxDistanceToRoot (root-excluded, collapse-unaware). Pure. */
    static double cappedMaxDistanceToRoot( final Phylogeny phy, final double cap ) {
        if ( ( phy == null ) || phy.isEmpty() || ( cap <= 0 ) ) {
            return 0;
        }
        double max = 0;
        for( final PhylogenyNode ext : phy.getExternalNodes() ) {
            final double d = cappedDistanceToRoot( ext, cap );
            if ( d > max ) {
                max = d;
            }
        }
        return max;
    }

    private static double cappedDistance( final PhylogenyNode n, final double cap ) {
        final double dtp = n.getDistanceToParent();
        if ( dtp <= 0 ) {
            return 0;
        }
        return ( dtp > cap ) ? cap : dtp;
    }

    /** The "nice" scale-bar distance (0.01 / 0.1 / 1 / 10 / 100) for a tree whose deepest root-to-tip distance is
     *  {@code height}, so the bar is a readable fraction of the tree width; 0 for a non-positive height. Pure. Used
     *  both for the normal scale bar (from the true max distance) and, while capping, from the DRAWN (capped) extent
     *  so the bar reflects the un-broken (ingroup) scale rather than the outlier-inflated one. */
    static double niceScaleBarDistance( final double height ) {
        if ( height <= 0 ) {
            return 0.0;
        }
        if ( height <= 0.5 ) {
            return 0.01;
        }
        if ( height <= 5.0 ) {
            return 0.1;
        }
        if ( height <= 50.0 ) {
            return 1;
        }
        if ( height <= 500.0 ) {
            return 10;
        }
        return 100;
    }

    /**
     * Provenance sentence for a ladderize (order-appearance) operation -- pure/testable; the caller appends it
     * to the phylogeny description. {@code larger_first} is {@code null} when the direction is not tracked (the
     * per-node "Ladderize Subtree" action auto-picks a direction); {@code TRUE} means the larger clade is drawn
     * first at each node (the whole-tree "order all" toggle in its {@code order=true} state), {@code FALSE} the
     * smaller clade first. See {@link #ladderizeStateOf}.
     */
    final static String ladderizeProvenanceSentence( final boolean whole_tree, final Boolean larger_first,
                                                      final String tree_name, final int num_ext_nodes ) {
        final String scope = whole_tree ? "the whole tree" : "a subtree";
        final String dir = ( larger_first == null ) ? ""
                : ( " so that " + ( larger_first ? "larger" : "smaller" ) + " clades are drawn first" );
        final String name = ForesterUtil.isEmpty( tree_name ) ? "" : ( " named \"" + tree_name + "\"" );
        return "Ladderized " + scope + dir + " (tree" + name + " with " + num_ext_nodes + " tips).";
    }

    /**
     * The direction the subtree rooted at {@code root} is currently ladderized in -- used to sync the "order all"
     * toggle icon to the tree after an Undo/Redo or navigation. {@code TRUE} = larger clade drawn first at every
     * 2-child node ({@code orderAppearance}'s {@code order=true} result); {@code FALSE} = smaller first; {@code
     * null} = mixed, or not ladderized (no 2-child node with unequal-sized children). Requires the external-node
     * counts to be current. Pure/testable.
     */
    final static Boolean ladderizeStateOf( final PhylogenyNode root ) {
        final boolean[] flags = { false, true, true }; // { any-unequal, all-larger-first, all-smaller-first }
        collectLadderizeState( root, flags );
        if ( !flags[ 0 ] ) {
            return null;
        }
        if ( flags[ 1 ] ) {
            return Boolean.TRUE;
        }
        if ( flags[ 2 ] ) {
            return Boolean.FALSE;
        }
        return null;
    }

    private static void collectLadderizeState( final PhylogenyNode n, final boolean[] flags ) {
        // Every ADJACENT child pair with unequal clade sizes is evidence -- over however many children the
        // node has (a sorted sequence is exactly one whose adjacent pairs are all ordered), so polytomies
        // participate in the direction detection just as the generalized ladderize now reorders them.
        // Identical to the old logic for a 2-child node.
        for ( int i = 0; i + 1 < n.getNumberOfDescendants(); i++ ) {
            final int a = n.getChildNode( i ).getNumberOfExternalNodes();
            final int b = n.getChildNode( i + 1 ).getNumberOfExternalNodes();
            if ( a != b ) {
                flags[ 0 ] = true;
                if ( a > b ) {
                    flags[ 2 ] = false; // an earlier child larger -> not smaller-first
                }
                else {
                    flags[ 1 ] = false; // an earlier child smaller -> not larger-first
                }
            }
        }
        for ( int i = 0; i < n.getNumberOfDescendants(); i++ ) {
            collectLadderizeState( n.getChildNode( i ), flags );
        }
    }
}
