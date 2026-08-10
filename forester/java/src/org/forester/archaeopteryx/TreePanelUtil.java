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
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Function;

import javax.swing.JOptionPane;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;
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

    public static List<String> createUrisForPdbWeb( final PhylogenyNode node,
                                                    final List<Accession> pdb_accs,
                                                    final Configuration configuration,
                                                    final TreePanel treePanel ) {
        final List<String> uris = new ArrayList<String>();
        if ( !ForesterUtil.isEmpty( pdb_accs ) ) {
            for( final Accession pdb_acc : pdb_accs ) {
                if ( !ForesterUtil.isEmpty( pdb_acc.getValue() ) ) {
                    uris.add( ForesterUtil.PDB + pdb_acc.getValue() );
                }
            }
        }
        return uris;
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

    /** Sentinel for {@link #maximalMonochromaticRoots}: a subtree whose tips are not all one rank taxon. */
    private final static String MIXED_TAXON = "<<MIXED>>";

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
        final Map<PhylogenyNode, String> assignment = assignTipsToRankTaxon( tree, rank, service );
        final SortedSet<String> taxa = new TreeSet<String>( assignment.values() );
        final Map<String, Color> colors = AptxUtil.assignDistinctColors( taxa );
        applyColorOverrides( colors, overrides );
        final Map<PhylogenyNode, String> roots = maximalMonochromaticRoots( tree, assignment );
        int colorizations = 0;
        for( final Entry<PhylogenyNode, String> e : roots.entrySet() ) {
            final Color c = colors.get( e.getValue() );
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

    /** Fills {@code counts_out} (when non-null) with the number of tips assigned to each taxon. */
    private static void countTipsPerTaxon( final Map<PhylogenyNode, String> assignment,
                                           final Map<String, Integer> counts_out ) {
        if ( counts_out == null ) {
            return;
        }
        for( final String taxon : assignment.values() ) {
            counts_out.merge( taxon, 1, Integer::sum );
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
        final Map<PhylogenyNode, String> assignment = assignTipsToRankTaxon( tree, rank, service );
        final Map<String, Color> colors = AptxUtil.assignDistinctColors( new TreeSet<String>( assignment.values() ) );
        applyColorOverrides( colors, overrides );
        for( final Entry<PhylogenyNode, String> e : maximalMonochromaticRoots( tree, assignment ).entrySet() ) {
            final Color c = colors.get( e.getValue() );
            if ( c != null ) {
                bands.add( new CladeBand( e.getValue(), c, e.getKey() ) );
            }
        }
        countTipsPerTaxon( assignment, counts_out );
        return bands;
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
     * Maps each external node to its taxon at {@code rank}, omitting tips that cannot be placed.
     * Resolution order per tip -- <i>tip identity wins, internal annotations only fill gaps</i>: (1) the tip's
     * OWN taxonomy annotated at exactly that rank; else (2) the tip's cached {@link TaxonLineage} from
     * {@code service} (no network here -- a cache miss just skips this step); else (3) the nearest
     * PROPER-ancestor annotation at that rank -- a fallback only for tips that cannot resolve their own
     * identity, so a wrong/partial internal-node annotation can no longer override a tip that resolves correctly.
     * <p>
     * Trade-off of this precedence: a wrong DB resolution of an ambiguous tip NAME (step 2) now beats a correct
     * in-tree ANCESTOR annotation (step 3). The escape hatch is step 1 -- annotate the TIP itself at {@code rank}
     * (its own taxonomy always wins), which is exactly how a curator overrides a bad DB hit.
     */
    final static Map<PhylogenyNode, String> assignTipsToRankTaxon( final Phylogeny tree,
                                                                   final String rank,
                                                                   final TaxonomicLineageService service ) {
        final Map<PhylogenyNode, String> assignment = new HashMap<PhylogenyNode, String>();
        if ( ( tree == null ) || tree.isEmpty() || ForesterUtil.isEmpty( rank ) ) {
            return assignment;
        }
        for( final PhylogenyNodeIterator it = tree.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode tip = it.next();
            // (1) the tip's OWN taxonomy at this rank -- the most specific identity, always trusted
            String taxon = selfRankTaxon( tip, rank );
            // (2) else resolve the tip's OWN name against the (cached) lineage DB
            if ( ForesterUtil.isEmpty( taxon ) && ( service != null ) ) {
                final String q = tipQueryName( tip );
                if ( !ForesterUtil.isEmpty( q ) ) {
                    final TaxonLineage rl = service.lineageOf( q );
                    if ( rl != null ) {
                        taxon = rl.at( rank );
                    }
                }
            }
            // (3) else fall back to the nearest ANCESTOR's manual annotation at this rank
            if ( ForesterUtil.isEmpty( taxon ) ) {
                taxon = ancestorRankTaxon( tip, rank );
            }
            if ( !ForesterUtil.isEmpty( taxon ) ) {
                assignment.put( tip, taxon );
            }
        }
        return assignment;
    }

    /** The taxon label of {@code node}'s OWN taxonomy iff its rank equals {@code rank} (case-insensitive) and the
     *  label is non-empty, else null. The tip's own identity -- highest priority in {@link #assignTipsToRankTaxon}. */
    final static String selfRankTaxon( final PhylogenyNode node, final String rank ) {
        if ( node.getNodeData().isHasTaxonomy() ) {
            final Taxonomy tax = node.getNodeData().getTaxonomy();
            if ( !ForesterUtil.isEmpty( tax.getRank() ) && tax.getRank().equalsIgnoreCase( rank ) ) {
                final String label = taxonomyLabel( tax );
                if ( !ForesterUtil.isEmpty( label ) ) {
                    return label;
                }
            }
        }
        return null;
    }

    /** The taxon label on the nearest PROPER ancestor carrying exactly {@code rank}, or null -- the fallback used
     *  in {@link #assignTipsToRankTaxon} only for a tip that cannot resolve its own identity. */
    final static String ancestorRankTaxon( final PhylogenyNode tip, final String rank ) {
        for( PhylogenyNode n = tip.getParent(); n != null; n = n.getParent() ) {
            final String label = selfRankTaxon( n, rank );
            if ( !ForesterUtil.isEmpty( label ) ) {
                return label;
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
     * Each node that roots a <i>maximal</i> clade whose external descendants all share one rank
     * taxon, mapped to that taxon. A node qualifies iff its whole subtree is uniform in
     * {@code assignment} and its parent's subtree is not the same taxon (so only the topmost such
     * node is returned). Handles paraphyly: a taxon split across the tree yields several roots, all
     * mapping to the same taxon (hence the same color). A tip with no assignment breaks uniformity,
     * so an unplaced tip is never swept into a neighboring clade's color.
     */
    final static Map<PhylogenyNode, String> maximalMonochromaticRoots( final Phylogeny tree,
                                                                       final Map<PhylogenyNode, String> assignment ) {
        final Map<PhylogenyNode, String> subtree = new HashMap<PhylogenyNode, String>();
        final Map<PhylogenyNode, String> roots = new LinkedHashMap<PhylogenyNode, String>();
        if ( ( tree == null ) || tree.isEmpty() ) {
            return roots;
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isExternal() ) {
                final String t = assignment.get( n );
                subtree.put( n, ( t != null ) ? t : MIXED_TAXON );
            }
            else {
                String uniform = null;
                boolean mixed = false;
                for( final PhylogenyNode c : n.getDescendants() ) {
                    final String cs = subtree.get( c );
                    if ( ( cs == null ) || cs.equals( MIXED_TAXON ) ) {
                        mixed = true;
                        break;
                    }
                    if ( uniform == null ) {
                        uniform = cs;
                    }
                    else if ( !uniform.equals( cs ) ) {
                        mixed = true;
                        break;
                    }
                }
                subtree.put( n, ( !mixed && ( uniform != null ) ) ? uniform : MIXED_TAXON );
            }
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final String t = subtree.get( n );
            if ( ( t != null ) && !t.equals( MIXED_TAXON ) ) {
                final PhylogenyNode p = n.getParent();
                if ( ( p == null ) || !t.equals( subtree.get( p ) ) ) {
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

    final static String pdbAccToString( final List<Accession> accs, final int i ) {
        if ( ForesterUtil.isEmpty( accs.get( i ).getComment() ) ) {
            return accs.get( i ).getValue();
        }
        return accs.get( i ).getValue() + " (" + accs.get( i ).getComment().toLowerCase() + ")";
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
                if ( isInternalPropertyRef( p.getRef() ) ) {
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
}
