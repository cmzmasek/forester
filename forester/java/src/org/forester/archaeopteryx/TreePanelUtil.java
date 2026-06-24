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
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.JOptionPane;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.NodeDataField;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;
import org.forester.util.StringInt;
import org.forester.ws.seqdb.AccessionAwareLineageService;
import org.forester.ws.seqdb.NcbiTaxonomyLineageService;
import org.forester.ws.seqdb.RankedLineage;
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

   

    public final static void showExtDescNodeDataUserSelectedHelper( final ControlPanel cp,
                                                                    final PhylogenyNode node,
                                                                    final List<String> data ) {
        final StringBuilder sb = new StringBuilder();
        if ( cp.isShowNodeNames() && !ForesterUtil.isEmpty( node.getName() ) ) {
            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getName(), sb );
        }
        if ( cp.isShowSeqNames() && node.getNodeData().isHasSequence()
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getName() ) ) {
            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence().getName(), sb );
        }
        if ( cp.isShowSeqSymbols() && node.getNodeData().isHasSequence()
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getSymbol() ) ) {
            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence().getSymbol(),
                                                                       sb );
        }
        if ( cp.isShowGeneNames() && node.getNodeData().isHasSequence()
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getGeneName() ) ) {
            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence().getGeneName(),
                                                                       sb );
        }
        if ( cp.isShowSequenceAcc() && node.getNodeData().isHasSequence()
                && ( node.getNodeData().getSequence().getAccession() != null )
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().toString() ) ) {
            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence().getAccession()
                    .toString(), sb );
        }
        if ( cp.isShowTaxonomyCode() && node.getNodeData().isHasTaxonomy()
                && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
            TreePanelUtil
                    .showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getTaxonomy().getTaxonomyCode(),
                                                                  sb );
        }
        if ( cp.isShowTaxonomyScientificNames() && node.getNodeData().isHasTaxonomy()
                && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) ) {
            TreePanelUtil
                    .showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getTaxonomy().getScientificName(),
                                                                  sb );
        }
        if ( cp.isShowTaxonomyCommonNames() && node.getNodeData().isHasTaxonomy()
                && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getCommonName() ) ) {
            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getTaxonomy().getCommonName(),
                                                                       sb );
        }
        //        if ( ( cp.isShowSeqNames() || cp.isShowSeqSymbols() || cp.isShowSequenceAcc() )
        //                && node.getNodeData().isHasSequence()
        //                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() ) ) {
        //            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence()
        //                    .getMolecularSequence(), sb );
        //        }
        final String s = sb.toString().trim();
        if ( !ForesterUtil.isEmpty( s ) ) {
            data.add( s );
        }
    }

    public final static void showExtDescNodeDataUserSelectedHelperHelper( final String s, final StringBuilder sb ) {
        if ( sb.length() > 0 ) {
            sb.append( "\t" );
        }
        sb.append( s );
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

    /** The support-scale ceiling implied by the largest value present: 100 if any value exceeds 1, else 1. */
    final static double confidenceScaleMaxFor( final double observed_max ) {
        return ( observed_max > 1.0 ) ? 100.0 : 1.0;
    }

    /** Scans a tree's internal-node confidences and returns the implied scale ceiling (1 or 100). */
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
        final double fade = SUPPORT_COLOR_MAX_FADE * ( 1.0 - f );
        return new Color( blend( strong.getRed(), background.getRed(), fade ),
                          blend( strong.getGreen(), background.getGreen(), fade ),
                          blend( strong.getBlue(), background.getBlue(), fade ) );
    }

    /** One 8-bit channel blended {@code t} (0..1) of the way from {@code a} toward {@code b}. */
    private static int blend( final int a, final int b, final double t ) {
        return (int) Math.round( a + ( t * ( b - a ) ) );
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
        return colorizations;
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
     * Resolution order per tip: (a) the nearest self-or-ancestor node annotated with exactly that
     * rank (free, in-tree); (b) the tip's cached {@link RankedLineage} from {@code service} (no
     * network here -- a cache miss simply leaves the tip unplaced).
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
            String taxon = inTreeRankTaxon( tip, rank );
            if ( ( taxon == null ) && ( service != null ) ) {
                final String q = tipQueryName( tip );
                if ( !ForesterUtil.isEmpty( q ) ) {
                    final RankedLineage rl = service.lineageOf( q );
                    if ( rl != null ) {
                        taxon = rl.at( rank );
                    }
                }
            }
            if ( !ForesterUtil.isEmpty( taxon ) ) {
                assignment.put( tip, taxon );
            }
        }
        return assignment;
    }

    /** The taxon label on the nearest self-or-ancestor node carrying exactly {@code rank}, or null. */
    final static String inTreeRankTaxon( final PhylogenyNode tip, final String rank ) {
        for( PhylogenyNode n = tip; n != null; n = n.getParent() ) {
            if ( n.getNodeData().isHasTaxonomy() ) {
                final Taxonomy tax = n.getNodeData().getTaxonomy();
                if ( !ForesterUtil.isEmpty( tax.getRank() ) && tax.getRank().equalsIgnoreCase( rank ) ) {
                    final String label = taxonomyLabel( tax );
                    if ( !ForesterUtil.isEmpty( label ) ) {
                        return label;
                    }
                }
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
     * The distinct taxon query-names of tips that are neither placeable from in-tree rank annotations
     * nor already in {@code service}'s cache -- i.e. exactly the names a caller must
     * {@link TaxonomicLineageService#fetch} (off the EDT) to place more tips at {@code rank}. Taxa
     * already in the cache are excluded even if the cache lacks {@code rank} (refetching would not
     * help), so a second call after a fetch pass returns an empty set (no repeated prompts).
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
            if ( inTreeRankTaxon( tip, rank ) != null ) {
                continue;
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

    /**
     * The first {@code max} entries of {@code in} in iteration order, for capping a legend; an
     * unmodifiable view is not needed -- callers treat it read-only. Used to bound the rank legend.
     */
    final static Map<String, Color> capEntries( final Map<String, Color> in, final int max ) {
        final Map<String, Color> out = new LinkedHashMap<String, Color>();
        if ( in != null ) {
            int i = 0;
            for( final Entry<String, Color> e : in.entrySet() ) {
                if ( i++ >= max ) {
                    break;
                }
                out.put( e.getKey(), e.getValue() );
            }
        }
        return out;
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

    static final int nodeDataIntoStringBuffer( final List<String> data, final Options optz, final StringBuilder sb ) {
        final SortedMap<String, Integer> map = new TreeMap<String, Integer>();
        int size = 0;
        if ( ( optz.getExtDescNodeDataToReturn() != NodeDataField.SEQUENCE_MOL_SEQ_FASTA )
                && ( optz.getExtDescNodeDataToReturn() != NodeDataField.GO_TERM_IDS ) ) {
            for( final String d : data ) {
                if ( !ForesterUtil.isEmpty( d ) ) {
                    if ( map.containsKey( d ) ) {
                        map.put( d, map.get( d ) + 1 );
                    }
                    else {
                        map.put( d, 1 );
                    }
                }
            }
            if ( ( optz.getExtDescNodeDataToReturn() == NodeDataField.DOMAINS_ALL )
                    || ( optz.getExtDescNodeDataToReturn() == NodeDataField.DOMAINS_COLLAPSED_PER_PROTEIN )
                    || ( optz.getExtDescNodeDataToReturn() == NodeDataField.SEQ_ANNOTATIONS ) ) {
                final ArrayList<StringInt> sis = new ArrayList<StringInt>();
                for( final Entry<String, Integer> e : map.entrySet() ) {
                    sis.add( new StringInt( e.getKey(), e.getValue() ) );
                }
                Collections.sort( sis, new StringInt.DescendingIntComparator() );
                for( final StringInt si : sis ) {
                    sb.append( si.getString() );
                    sb.append( "\t" );
                    sb.append( si.getInt() );
                    sb.append( ForesterUtil.LINE_SEPARATOR );
                }
            }
            else {
                for( final Entry<String, Integer> e : map.entrySet() ) {
                    final String v = e.getKey();
                    final Object c = e.getValue();
                    sb.append( v );
                    sb.append( "\t" );
                    sb.append( c );
                    sb.append( ForesterUtil.LINE_SEPARATOR );
                }
            }
            size = map.size();
        }
        else {
            for( final String d : data ) {
                if ( !ForesterUtil.isEmpty( d ) ) {
                    sb.append( d );
                    sb.append( ForesterUtil.LINE_SEPARATOR );
                }
            }
            size = data.size();
        }
        return size;
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
}
