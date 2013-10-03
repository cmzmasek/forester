
package org.forester.archaeopteryx;

import java.awt.Color;
import java.awt.Component;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;

import javax.swing.JOptionPane;

import org.forester.analysis.TaxonomyDataManager;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.NodeData.NODE_DATA;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;
import org.forester.ws.seqdb.UniProtTaxonomy;

public class TreePanelUtil {

    static int makeSB( final List<String> data, final Options optz, final StringBuilder sb ) {
        final SortedMap<String, Integer> map = new TreeMap<String, Integer>();
        if ( ( optz.getExtDescNodeDataToReturn() != NODE_DATA.SEQUENCE_MOL_SEQ )
                && ( optz.getExtDescNodeDataToReturn() != NODE_DATA.SEQUENCE_MOL_SEQ_FASTA ) ) {
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
        }
        int size = 0;
        if ( ( optz.getExtDescNodeDataToReturn() != NODE_DATA.SEQUENCE_MOL_SEQ )
                && ( optz.getExtDescNodeDataToReturn() != NODE_DATA.SEQUENCE_MOL_SEQ_FASTA ) ) {
            for( final Entry<String, Integer> e : map.entrySet() ) {
                final String v = e.getKey();
                final Object c = e.getValue();
                sb.append( v );
                sb.append( "\t" );
                sb.append( c );
                sb.append( ForesterUtil.LINE_SEPARATOR );
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

    public final static String createUriForSeqWeb( final PhylogenyNode node,
                                                   final Configuration conf,
                                                   final TreePanel tp ) {
        String uri_str = null;
        final String upkb = SequenceAccessionTools.extractUniProtKbProteinSeqIdentifier( node );
        if ( !ForesterUtil.isEmpty( upkb ) ) {
            try {
                uri_str = ForesterUtil.UNIPROT_KB + URLEncoder.encode( upkb, ForesterConstants.UTF8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                AptxUtil.showErrorMessage( tp, e.toString() );
                e.printStackTrace();
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = SequenceAccessionTools.extractGenbankAccessor( node );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    if ( SequenceAccessionTools.isProtein( v ) ) {
                        uri_str = ForesterUtil.NCBI_PROTEIN + URLEncoder.encode( v, ForesterConstants.UTF8 );
                    }
                    else {
                        uri_str = ForesterUtil.NCBI_NUCCORE + URLEncoder.encode( v, ForesterConstants.UTF8 );
                    }
                }
                catch ( final UnsupportedEncodingException e ) {
                    AptxUtil.showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = SequenceAccessionTools.extractRefSeqAccessor( node );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    if ( SequenceAccessionTools.isProtein( v ) ) {
                        uri_str = ForesterUtil.NCBI_PROTEIN + URLEncoder.encode( v, ForesterConstants.UTF8 );
                    }
                    else {
                        uri_str = ForesterUtil.NCBI_NUCCORE + URLEncoder.encode( v, ForesterConstants.UTF8 );
                    }
                }
                catch ( final UnsupportedEncodingException e ) {
                    AptxUtil.showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = SequenceAccessionTools.extractGInumber( node );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    uri_str = ForesterUtil.NCBI_GI + URLEncoder.encode( v, ForesterConstants.UTF8 );
                }
                catch ( final UnsupportedEncodingException e ) {
                    AptxUtil.showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        return uri_str;
    }

    /**
     * Returns the set of distinct taxonomies of
     * all external nodes of node.
     * If at least one the external nodes has no taxonomy,
     * null is returned.
     * 
     */
    public static Set<Taxonomy> obtainDistinctTaxonomies( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        final Set<Taxonomy> tax_set = new HashSet<Taxonomy>();
        for( final PhylogenyNode n : descs ) {
            if ( !n.getNodeData().isHasTaxonomy() || n.getNodeData().getTaxonomy().isEmpty() ) {
                return null;
            }
            tax_set.add( n.getNodeData().getTaxonomy() );
        }
        return tax_set;
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
            TreePanelUtil
                    .showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence().getSymbol(), sb );
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
            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getTaxonomy()
                    .getTaxonomyCode(), sb );
        }
        if ( cp.isShowTaxonomyScientificNames() && node.getNodeData().isHasTaxonomy()
                && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) ) {
            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getTaxonomy()
                    .getScientificName(), sb );
        }
        if ( cp.isShowTaxonomyCommonNames() && node.getNodeData().isHasTaxonomy()
                && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getCommonName() ) ) {
            TreePanelUtil
                    .showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getTaxonomy().getCommonName(), sb );
        }
        if ( ( cp.isShowSeqNames() || cp.isShowSeqSymbols() || cp.isShowSequenceAcc() )
                && node.getNodeData().isHasSequence()
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() ) ) {
            TreePanelUtil.showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence()
                    .getMolecularSequence(), sb );
        }
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

    final static Color calculateColorFromString( final String str, final boolean is_taxonomy ) {
        final String my_str = str.toUpperCase();
        char first = my_str.charAt( 0 );
        char second = ' ';
        char third = ' ';
        if ( my_str.length() > 1 ) {
            if ( is_taxonomy ) {
                second = my_str.charAt( 1 );
            }
            else {
                second = my_str.charAt( my_str.length() - 1 );
            }
            if ( is_taxonomy ) {
                if ( my_str.length() > 2 ) {
                    if ( my_str.indexOf( " " ) > 0 ) {
                        third = my_str.charAt( my_str.indexOf( " " ) + 1 );
                    }
                    else {
                        third = my_str.charAt( 2 );
                    }
                }
            }
            else if ( my_str.length() > 2 ) {
                third = my_str.charAt( ( my_str.length() - 1 ) / 2 );
            }
        }
        first = TreePanelUtil.normalizeCharForRGB( first );
        second = TreePanelUtil.normalizeCharForRGB( second );
        third = TreePanelUtil.normalizeCharForRGB( third );
        if ( ( first > 235 ) && ( second > 235 ) && ( third > 235 ) ) {
            first = 0;
        }
        else if ( ( first < 60 ) && ( second < 60 ) && ( third < 60 ) ) {
            second = 255;
        }
        return new Color( first, second, third );
    }

    final static void collapseSpeciesSpecificSubtrees( final Phylogeny phy ) {
        boolean inferred = false;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.isCollapse() && ( n.getNumberOfDescendants() > 1 ) ) {
                final Set<Taxonomy> taxs = TreePanelUtil.obtainDistinctTaxonomies( n );
                if ( ( taxs != null ) && ( taxs.size() == 1 ) ) {
                    TreePanelUtil.collapseSubtree( n, true );
                    if ( !n.getNodeData().isHasTaxonomy() ) {
                        n.getNodeData().setTaxonomy( ( Taxonomy ) n.getAllExternalDescendants().get( 0 ).getNodeData()
                                .getTaxonomy().copy() );
                    }
                    inferred = true;
                }
                else {
                    n.setCollapse( false );
                }
            }
        }
        if ( inferred ) {
            phy.setRerootable( false );
        }
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

    static void colorizeSubtree( final PhylogenyNode node, final BranchColor c ) {
        node.getBranchData().setBranchColor( c );
        final List<PhylogenyNode> descs = PhylogenyMethods.getAllDescendants( node );
        for( final PhylogenyNode desc : descs ) {
            desc.getBranchData().setBranchColor( c );
        }
    }

    final static void colorPhylogenyAccordingToConfidenceValues( final Phylogeny tree, final TreePanel tree_panel ) {
        double max_conf = 0.0;
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.getBranchData().setBranchColor( null );
            if ( n.getBranchData().isHasConfidences() ) {
                final double conf = PhylogenyMethods.getConfidenceValue( n );
                if ( conf > max_conf ) {
                    max_conf = conf;
                }
            }
        }
        if ( max_conf > 0.0 ) {
            final Color bg = tree_panel.getTreeColorSet().getBackgroundColor();
            final Color br = tree_panel.getTreeColorSet().getBranchColor();
            for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( n.getBranchData().isHasConfidences() ) {
                    final double conf = PhylogenyMethods.getConfidenceValue( n );
                    final BranchColor c = new BranchColor( ForesterUtil.calcColor( conf, 0.0, max_conf, bg, br ) );
                    TreePanelUtil.colorizeSubtree( n, c );
                }
            }
        }
    }

    final static void colorPhylogenyAccordingToExternalTaxonomy( final Phylogeny tree, final TreePanel tree_panel ) {
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            it.next().getBranchData().setBranchColor( null );
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.getBranchData().isHasBranchColor() ) {
                final Taxonomy tax = PhylogenyMethods.getExternalDescendantsTaxonomy( n );
                if ( tax != null ) {
                    n.getBranchData().setBranchColor( new BranchColor( tree_panel.calculateTaxonomyBasedColor( tax ) ) );
                    final List<PhylogenyNode> descs = PhylogenyMethods.getAllDescendants( n );
                    for( final PhylogenyNode desc : descs ) {
                        desc.getBranchData()
                                .setBranchColor( new BranchColor( tree_panel.calculateTaxonomyBasedColor( tax ) ) );
                    }
                }
            }
        }
    }

    final static int colorPhylogenyAccordingToRanks( final Phylogeny tree, final String rank, final TreePanel tree_panel ) {
        final Map<String, Color> true_lineage_to_color_map = new HashMap<String, Color>();
        int colorizations = 0;
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy()
                    && ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() )
                            || !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getCommonName() ) || !ForesterUtil
                            .isEmpty( n.getNodeData().getTaxonomy().getTaxonomyCode() ) ) ) {
                if ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getRank() )
                        && n.getNodeData().getTaxonomy().getRank().equalsIgnoreCase( rank ) ) {
                    final BranchColor c = new BranchColor( tree_panel.calculateTaxonomyBasedColor( n.getNodeData()
                            .getTaxonomy() ) );
                    TreePanelUtil.colorizeSubtree( n, c );
                    ++colorizations;
                    if ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                        true_lineage_to_color_map.put( n.getNodeData().getTaxonomy().getScientificName(), c.getValue() );
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
                            TreePanelUtil
                                    .colorizeSubtree( node, new BranchColor( true_lineage_to_color_map.get( lin ) ) );
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
                                if ( up.getRank().equalsIgnoreCase( rank ) ) {
                                    final BranchColor c = new BranchColor( tree_panel.calculateTaxonomyBasedColor( temp_tax ) );
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
        return colorizations;
    }

    final static String createAnnotationString( final SortedSet<Annotation> annotations, final boolean show_ref_sources ) {
        final SortedMap<String, List<Annotation>> m = new TreeMap<String, List<Annotation>>();
        for( final Annotation an : annotations ) {
            final String ref_source = ForesterUtil.isEmpty( an.getRefSource() ) ? "?" : an.getRefSource();
            if ( !m.containsKey( ref_source ) ) {
                m.put( ref_source, new ArrayList<Annotation>() );
            }
            m.get( ref_source ).add( an );
        }
        final StringBuilder sb = new StringBuilder();
        for( final Entry<String, List<Annotation>> e : m.entrySet() ) {
            final String ref_source = e.getKey();
            final List<Annotation> ans = e.getValue();
            if ( m.size() > 1 ) {
                sb.append( "[" );
            }
            if ( show_ref_sources && !ref_source.equals( "?" ) ) {
                sb.append( ref_source );
                sb.append( ": " );
            }
            for( int i = 0; i < ans.size(); ++i ) {
                final Annotation an = ans.get( i );
                if ( !ForesterUtil.isEmpty( an.getRefValue() ) ) {
                    sb.append( an.getRefValue() );
                    sb.append( " " );
                }
                if ( !ForesterUtil.isEmpty( an.getDesc() ) ) {
                    sb.append( an.getDesc() );
                }
                if ( sb.charAt( sb.length() - 1 ) == ' ' ) {
                    sb.deleteCharAt( sb.length() - 1 );
                }
                if ( i < ans.size() - 1 ) {
                    sb.append( ", " );
                }
            }
            if ( m.size() > 1 ) {
                sb.append( "] " );
            }
        }
        return sb.toString();
    }

    final static String getPartAfterColon( final String s ) {
        final int i = s.indexOf( ':' );
        if ( ( i < 1 ) || ( i == ( s.length() - 1 ) ) ) {
            return s;
        }
        return s.substring( i + 1, s.length() );
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
                && ForesterUtil.isEmpty( tax.getCommonName() ) && ForesterUtil.isEmpty( tax.getScientificName() ) && tax
                .getSynonyms().isEmpty() );
    }

    final static char normalizeCharForRGB( char c ) {
        c -= 65;
        c *= 10.2;
        c = c > 255 ? 255 : c;
        c = c < 0 ? 0 : c;
        return c;
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
}
