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
import java.util.TreeMap;

import javax.swing.JOptionPane;

import org.forester.analysis.TaxonomyDataManager;
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
import org.forester.ws.seqdb.UniProtTaxonomy;

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

    /** THRESHOLD_MARKS test: is the support at or above the cutoff (a fraction 0..1 of the scale)? */
    final static boolean isSupportAtOrAboveThreshold( final double confidence,
                                                      final double scale_max,
                                                      final double threshold_fraction ) {
        return supportFraction( confidence, scale_max ) >= threshold_fraction;
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
     * Colorizes subtrees by the taxonomy at the given {@code rank}. When {@code legend_out} is
     * non-null it is filled with the taxon-name &rarr; color pairs used, so the caller can show a
     * legend mapping colors to taxa.
     */
    final static int colorPhylogenyAccordingToRanks( final Phylogeny tree,
                                                     final String rank,
                                                     final TreePanel tree_panel,
                                                     final Map<String, Color> legend_out ) {
        final Map<String, Color> true_lineage_to_color_map = new HashMap<String, Color>();
        int colorizations = 0;
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy()
                    && ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() )
                            || !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getCommonName() )
                            || !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getTaxonomyCode() ) ) ) {
                if ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getRank() )
                        && n.getNodeData().getTaxonomy().getRank().equalsIgnoreCase( rank ) ) {
                    final BranchColor c = new BranchColor( tree_panel
                            .calculateTaxonomyBasedColor( n.getNodeData().getTaxonomy() ) );
                    TreePanelUtil.colorizeSubtree( n, c );
                    ++colorizations;
                    // legend label: scientific name, else common name, else taxonomy code -- so taxa
                    // identified only by a code/common name still get a legend row (for scientific
                    // names this key also doubles as the lineage-match key used in the next pass)
                    final String label = taxonomyLabel( n.getNodeData().getTaxonomy() );
                    if ( !ForesterUtil.isEmpty( label ) ) {
                        true_lineage_to_color_map.put( label, c.getValue() );
                    }
                }
            }
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( ( node.getBranchData().getBranchColor() == null ) && node.getNodeData().isHasTaxonomy()
                    && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getLineage() ) ) {
                boolean success = false;
                if ( !true_lineage_to_color_map.isEmpty() ) {
                    for( final String lin : node.getNodeData().getTaxonomy().getLineage() ) {
                        if ( true_lineage_to_color_map.containsKey( lin ) ) {
                            TreePanelUtil.colorizeSubtree( node,
                                                           new BranchColor( true_lineage_to_color_map.get( lin ) ) );
                            ++colorizations;
                            success = true;
                            break;
                        }
                    }
                }
                if ( !success ) {
                    final Map<String, String> lineage_to_rank_map = MainPanel.getLineageToRankMap();
                    for( final String lin : node.getNodeData().getTaxonomy().getLineage() ) {
                        final Taxonomy temp_tax = new Taxonomy();
                        temp_tax.setScientificName( lin );
                        if ( lineage_to_rank_map.containsKey( lin )
                                && !ForesterUtil.isEmpty( lineage_to_rank_map.get( lin ) )
                                && lineage_to_rank_map.get( lin ).equalsIgnoreCase( rank ) ) {
                            final BranchColor c = new BranchColor( tree_panel.calculateTaxonomyBasedColor( temp_tax ) );
                            TreePanelUtil.colorizeSubtree( node, c );
                            ++colorizations;
                            true_lineage_to_color_map.put( lin, c.getValue() );
                            break;
                        }
                        else {
                            UniProtTaxonomy up = null;
                            try {
                                up = TaxonomyDataManager.obtainUniProtTaxonomy( temp_tax, null, null );
                            }
                            catch ( final Exception e ) {
                                e.printStackTrace();
                            }
                            if ( ( up != null ) && !ForesterUtil.isEmpty( up.getRank() ) ) {
                                lineage_to_rank_map.put( lin, up.getRank() );
                                System.out.println( lin + "->" + up.getRank() );
                                if ( up.getRank().equalsIgnoreCase( rank ) ) {
                                    final BranchColor c = new BranchColor( tree_panel
                                            .calculateTaxonomyBasedColor( temp_tax ) );
                                    TreePanelUtil.colorizeSubtree( node, c );
                                    ++colorizations;
                                    true_lineage_to_color_map.put( lin, c.getValue() );
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        if ( legend_out != null ) {
            legend_out.putAll( true_lineage_to_color_map );
        }
        return colorizations;
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
