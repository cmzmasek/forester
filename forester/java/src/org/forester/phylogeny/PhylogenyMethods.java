// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.phylogeny;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.io.parsers.util.PhylogenyParserException;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class PhylogenyMethods {

    private PhylogenyMethods() {
        // Hidden constructor.
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();
    }

    public static DescriptiveStatistics calculatBranchLengthStatistics( final Phylogeny phy ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( !n.isRoot() && ( n.getDistanceToParent() >= 0.0 ) ) {
                stats.addValue( n.getDistanceToParent() );
            }
        }
        return stats;
    }

    public static List<DescriptiveStatistics> calculatConfidenceStatistics( final Phylogeny phy ) {
        final List<DescriptiveStatistics> stats = new ArrayList<DescriptiveStatistics>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( !n.isExternal() && !n.isRoot() ) {
                if ( n.getBranchData().isHasConfidences() ) {
                    for( int i = 0; i < n.getBranchData().getConfidences().size(); ++i ) {
                        final Confidence c = n.getBranchData().getConfidences().get( i );
                        if ( ( i > ( stats.size() - 1 ) ) || ( stats.get( i ) == null ) ) {
                            stats.add( i, new BasicDescriptiveStatistics() );
                        }
                        if ( !ForesterUtil.isEmpty( c.getType() ) ) {
                            if ( !ForesterUtil.isEmpty( stats.get( i ).getDescription() ) ) {
                                if ( !stats.get( i ).getDescription().equalsIgnoreCase( c.getType() ) ) {
                                    throw new IllegalArgumentException( "support values in node [" + n.toString()
                                            + "] appear inconsistently ordered" );
                                }
                            }
                            stats.get( i ).setDescription( c.getType() );
                        }
                        stats.get( i ).addValue( ( ( c != null ) && ( c.getValue() >= 0 ) ) ? c.getValue() : 0 );
                    }
                }
            }
        }
        return stats;
    }

    /**
     * Calculates the distance between PhylogenyNodes node1 and node2.
     * 
     * 
     * @param node1
     * @param node2
     * @return distance between node1 and node2
     */
    public static double calculateDistance( final PhylogenyNode node1, final PhylogenyNode node2 ) {
        final PhylogenyNode lca = calculateLCA( node1, node2 );
        final PhylogenyNode n1 = node1;
        final PhylogenyNode n2 = node2;
        return ( PhylogenyMethods.getDistance( n1, lca ) + PhylogenyMethods.getDistance( n2, lca ) );
    }

    /**
     * Returns the LCA of PhylogenyNodes node1 and node2.
     * 
     * 
     * @param node1
     * @param node2
     * @return LCA of node1 and node2
     */
    public final static PhylogenyNode calculateLCA( PhylogenyNode node1, PhylogenyNode node2 ) {
        if ( node1 == null ) {
            throw new IllegalArgumentException( "first argument (node) is null" );
        }
        if ( node2 == null ) {
            throw new IllegalArgumentException( "second argument (node) is null" );
        }
        if ( node1 == node2 ) {
            return node1;
        }
        if ( ( node1.getParent() == node2.getParent() ) ) {
            return node1.getParent();
        }
        int depth1 = node1.calculateDepth();
        int depth2 = node2.calculateDepth();
        while ( ( depth1 > -1 ) && ( depth2 > -1 ) ) {
            if ( depth1 > depth2 ) {
                node1 = node1.getParent();
                depth1--;
            }
            else if ( depth2 > depth1 ) {
                node2 = node2.getParent();
                depth2--;
            }
            else {
                if ( node1 == node2 ) {
                    return node1;
                }
                node1 = node1.getParent();
                node2 = node2.getParent();
                depth1--;
                depth2--;
            }
        }
        throw new IllegalArgumentException( "illegal attempt to calculate LCA of two nodes which do not share a common root" );
    }

    /**
     * Returns the LCA of PhylogenyNodes node1 and node2.
     * Precondition: ids are in pre-order (or level-order).
     * 
     * 
     * @param node1
     * @param node2
     * @return LCA of node1 and node2
     */
    public final static PhylogenyNode calculateLCAonTreeWithIdsInPreOrder( PhylogenyNode node1, PhylogenyNode node2 ) {
        if ( node1 == null ) {
            throw new IllegalArgumentException( "first argument (node) is null" );
        }
        if ( node2 == null ) {
            throw new IllegalArgumentException( "second argument (node) is null" );
        }
        while ( node1 != node2 ) {
            if ( node1.getId() > node2.getId() ) {
                node1 = node1.getParent();
            }
            else {
                node2 = node2.getParent();
            }
        }
        return node1;
    }

    public static short calculateMaxBranchesToLeaf( final PhylogenyNode node ) {
        if ( node.isExternal() ) {
            return 0;
        }
        short max = 0;
        for( PhylogenyNode d : node.getAllExternalDescendants() ) {
            short steps = 0;
            while ( d != node ) {
                if ( d.isCollapse() ) {
                    steps = 0;
                }
                else {
                    steps++;
                }
                d = d.getParent();
            }
            if ( max < steps ) {
                max = steps;
            }
        }
        return max;
    }

    public static int calculateMaxDepth( final Phylogeny phy ) {
        int max = 0;
        for( final PhylogenyNodeIterator iter = phy.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final int steps = node.calculateDepth();
            if ( steps > max ) {
                max = steps;
            }
        }
        return max;
    }

    public static double calculateMaxDistanceToRoot( final Phylogeny phy ) {
        double max = 0.0;
        for( final PhylogenyNodeIterator iter = phy.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final double d = node.calculateDistanceToRoot();
            if ( d > max ) {
                max = d;
            }
        }
        return max;
    }

    public static int calculateNumberOfExternalNodesWithoutTaxonomy( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        int x = 0;
        for( final PhylogenyNode n : descs ) {
            if ( !n.getNodeData().isHasTaxonomy() || n.getNodeData().getTaxonomy().isEmpty() ) {
                x++;
            }
        }
        return x;
    }

    public static DescriptiveStatistics calculatNumberOfDescendantsPerNodeStatistics( final Phylogeny phy ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( !n.isExternal() ) {
                stats.addValue( n.getNumberOfDescendants() );
            }
        }
        return stats;
    }

    public final static void collapseSubtreeStructure( final PhylogenyNode n ) {
        final List<PhylogenyNode> eds = n.getAllExternalDescendants();
        final List<Double> d = new ArrayList<Double>();
        for( final PhylogenyNode ed : eds ) {
            d.add( calculateDistanceToAncestor( n, ed ) );
        }
        for( int i = 0; i < eds.size(); ++i ) {
            n.setChildNode( i, eds.get( i ) );
            eds.get( i ).setDistanceToParent( d.get( i ) );
        }
    }

    public static int countNumberOfOneDescendantNodes( final Phylogeny phy ) {
        int count = 0;
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( !n.isExternal() && ( n.getNumberOfDescendants() == 1 ) ) {
                count++;
            }
        }
        return count;
    }

    public static int countNumberOfPolytomies( final Phylogeny phy ) {
        int count = 0;
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( !n.isExternal() && ( n.getNumberOfDescendants() > 2 ) ) {
                count++;
            }
        }
        return count;
    }

    public static final HashMap<String, PhylogenyNode> createNameToExtNodeMap( final Phylogeny phy ) {
        final HashMap<String, PhylogenyNode> nodes = new HashMap<String, PhylogenyNode>();
        final List<PhylogenyNode> ext = phy.getExternalNodes();
        for( final PhylogenyNode n : ext ) {
            nodes.put( n.getName(), n );
        }
        return nodes;
    }

    public static void deleteExternalNodesNegativeSelection( final Set<Long> to_delete, final Phylogeny phy ) {
        for( final Long id : to_delete ) {
            phy.deleteSubtree( phy.getNode( id ), true );
        }
        phy.clearHashIdToNodeMap();
        phy.externalNodesHaveChanged();
    }

    public static void deleteExternalNodesNegativeSelection( final String[] node_names_to_delete, final Phylogeny p )
            throws IllegalArgumentException {
        for( final String element : node_names_to_delete ) {
            if ( ForesterUtil.isEmpty( element ) ) {
                continue;
            }
            List<PhylogenyNode> nodes = null;
            nodes = p.getNodes( element );
            final Iterator<PhylogenyNode> it = nodes.iterator();
            while ( it.hasNext() ) {
                final PhylogenyNode n = it.next();
                if ( !n.isExternal() ) {
                    throw new IllegalArgumentException( "attempt to delete non-external node \"" + element + "\"" );
                }
                p.deleteSubtree( n, true );
            }
        }
        p.clearHashIdToNodeMap();
        p.externalNodesHaveChanged();
    }

    public static List<String> deleteExternalNodesPositiveSelection( final String[] node_names_to_keep,
                                                                     final Phylogeny p ) {
        final PhylogenyNodeIterator it = p.iteratorExternalForward();
        final String[] to_delete = new String[ p.getNumberOfExternalNodes() ];
        int i = 0;
        Arrays.sort( node_names_to_keep );
        while ( it.hasNext() ) {
            final String curent_name = it.next().getName();
            if ( Arrays.binarySearch( node_names_to_keep, curent_name ) < 0 ) {
                to_delete[ i++ ] = curent_name;
            }
        }
        PhylogenyMethods.deleteExternalNodesNegativeSelection( to_delete, p );
        final List<String> deleted = new ArrayList<String>();
        for( final String n : to_delete ) {
            if ( !ForesterUtil.isEmpty( n ) ) {
                deleted.add( n );
            }
        }
        return deleted;
    }

    public static void deleteExternalNodesPositiveSelectionT( final List<Taxonomy> species_to_keep, final Phylogeny phy ) {
        final Set<Long> to_delete = new HashSet<Long>();
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                if ( !species_to_keep.contains( n.getNodeData().getTaxonomy() ) ) {
                    to_delete.add( n.getId() );
                }
            }
            else {
                throw new IllegalArgumentException( "node " + n.getId() + " has no taxonomic data" );
            }
        }
        deleteExternalNodesNegativeSelection( to_delete, phy );
    }

    final public static void deleteInternalNodesWithOnlyOneDescendent( final Phylogeny phy ) {
        final ArrayList<PhylogenyNode> to_delete = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( ( !n.isExternal() ) && ( n.getNumberOfDescendants() == 1 ) ) {
                to_delete.add( n );
            }
        }
        for( final PhylogenyNode d : to_delete ) {
            PhylogenyMethods.removeNode( d, phy );
        }
        phy.clearHashIdToNodeMap();
        phy.externalNodesHaveChanged();
    }

    final public static void deleteNonOrthologousExternalNodes( final Phylogeny phy, final PhylogenyNode n ) {
        if ( n.isInternal() ) {
            throw new IllegalArgumentException( "node is not external" );
        }
        final ArrayList<PhylogenyNode> to_delete = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode i = it.next();
            if ( !PhylogenyMethods.getEventAtLCA( n, i ).isSpeciation() ) {
                to_delete.add( i );
            }
        }
        for( final PhylogenyNode d : to_delete ) {
            phy.deleteSubtree( d, true );
        }
        phy.clearHashIdToNodeMap();
        phy.externalNodesHaveChanged();
    }

    public static List<PhylogenyNode> getAllDescendants( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = new ArrayList<PhylogenyNode>();
        final Set<Long> encountered = new HashSet<Long>();
        if ( !node.isExternal() ) {
            final List<PhylogenyNode> exts = node.getAllExternalDescendants();
            for( PhylogenyNode current : exts ) {
                descs.add( current );
                while ( current != node ) {
                    current = current.getParent();
                    if ( encountered.contains( current.getId() ) ) {
                        continue;
                    }
                    descs.add( current );
                    encountered.add( current.getId() );
                }
            }
        }
        return descs;
    }

    /**
     * 
     * Convenience method
     * 
     * @param node
     * @return
     */
    public static Color getBranchColorValue( final PhylogenyNode node ) {
        if ( node.getBranchData().getBranchColor() == null ) {
            return null;
        }
        return node.getBranchData().getBranchColor().getValue();
    }

    /**
     * Convenience method
     */
    public static double getBranchWidthValue( final PhylogenyNode node ) {
        if ( !node.getBranchData().isHasBranchWidth() ) {
            return BranchWidth.BRANCH_WIDTH_DEFAULT_VALUE;
        }
        return node.getBranchData().getBranchWidth().getValue();
    }

    /**
     * Convenience method
     */
    public static double getConfidenceValue( final PhylogenyNode node ) {
        if ( !node.getBranchData().isHasConfidences() ) {
            return Confidence.CONFIDENCE_DEFAULT_VALUE;
        }
        return node.getBranchData().getConfidence( 0 ).getValue();
    }

    /**
     * Convenience method
     */
    public static double[] getConfidenceValuesAsArray( final PhylogenyNode node ) {
        if ( !node.getBranchData().isHasConfidences() ) {
            return new double[ 0 ];
        }
        final double[] values = new double[ node.getBranchData().getConfidences().size() ];
        int i = 0;
        for( final Confidence c : node.getBranchData().getConfidences() ) {
            values[ i++ ] = c.getValue();
        }
        return values;
    }

    final public static Event getEventAtLCA( final PhylogenyNode n1, final PhylogenyNode n2 ) {
        return calculateLCA( n1, n2 ).getNodeData().getEvent();
    }

    /**
     * Returns taxonomy t if all external descendants have 
     * the same taxonomy t, null otherwise.
     * 
     */
    public static Taxonomy getExternalDescendantsTaxonomy( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        Taxonomy tax = null;
        for( final PhylogenyNode n : descs ) {
            if ( !n.getNodeData().isHasTaxonomy() || n.getNodeData().getTaxonomy().isEmpty() ) {
                return null;
            }
            else if ( tax == null ) {
                tax = n.getNodeData().getTaxonomy();
            }
            else if ( n.getNodeData().getTaxonomy().isEmpty() || !tax.isEqual( n.getNodeData().getTaxonomy() ) ) {
                return null;
            }
        }
        return tax;
    }

    public static PhylogenyNode getFurthestDescendant( final PhylogenyNode node ) {
        final List<PhylogenyNode> children = node.getAllExternalDescendants();
        PhylogenyNode farthest = null;
        double longest = -Double.MAX_VALUE;
        for( final PhylogenyNode child : children ) {
            if ( PhylogenyMethods.getDistance( child, node ) > longest ) {
                farthest = child;
                longest = PhylogenyMethods.getDistance( child, node );
            }
        }
        return farthest;
    }

    // public static PhylogenyMethods getInstance() {
    //     if ( PhylogenyMethods._instance == null ) {
    //         PhylogenyMethods._instance = new PhylogenyMethods();
    //    }
    //    return PhylogenyMethods._instance;
    //  }
    /**
     * Returns the largest confidence value found on phy.
     */
    static public double getMaximumConfidenceValue( final Phylogeny phy ) {
        double max = -Double.MAX_VALUE;
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final double s = PhylogenyMethods.getConfidenceValue( iter.next() );
            if ( ( s != Confidence.CONFIDENCE_DEFAULT_VALUE ) && ( s > max ) ) {
                max = s;
            }
        }
        return max;
    }

    static public int getMinimumDescendentsPerInternalNodes( final Phylogeny phy ) {
        int min = Integer.MAX_VALUE;
        int d = 0;
        PhylogenyNode n;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            n = it.next();
            if ( n.isInternal() ) {
                d = n.getNumberOfDescendants();
                if ( d < min ) {
                    min = d;
                }
            }
        }
        return min;
    }

    /**
     * Convenience method for display purposes.
     * Not intended for algorithms.
     */
    public static String getSpecies( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasTaxonomy() ) {
            return "";
        }
        else if ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) ) {
            return node.getNodeData().getTaxonomy().getScientificName();
        }
        if ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
            return node.getNodeData().getTaxonomy().getTaxonomyCode();
        }
        else {
            return node.getNodeData().getTaxonomy().getCommonName();
        }
    }

    /**
     * Convenience method for display purposes.
     * Not intended for algorithms.
     */
    public static String getTaxonomyIdentifier( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasTaxonomy() || ( node.getNodeData().getTaxonomy().getIdentifier() == null ) ) {
            return "";
        }
        return node.getNodeData().getTaxonomy().getIdentifier().getValue();
    }

    public final static boolean isAllDecendentsAreDuplications( final PhylogenyNode n ) {
        if ( n.isExternal() ) {
            return true;
        }
        else {
            if ( n.isDuplication() ) {
                for( final PhylogenyNode desc : n.getDescendants() ) {
                    if ( !isAllDecendentsAreDuplications( desc ) ) {
                        return false;
                    }
                }
                return true;
            }
            else {
                return false;
            }
        }
    }

    public static boolean isHasExternalDescendant( final PhylogenyNode node ) {
        for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            if ( node.getChildNode( i ).isExternal() ) {
                return true;
            }
        }
        return false;
    }

    /*
     * This is case insensitive.
     * 
     */
    public synchronized static boolean isTaxonomyHasIdentifierOfGivenProvider( final Taxonomy tax,
                                                                               final String[] providers ) {
        if ( ( tax.getIdentifier() != null ) && !ForesterUtil.isEmpty( tax.getIdentifier().getProvider() ) ) {
            final String my_tax_prov = tax.getIdentifier().getProvider();
            for( final String provider : providers ) {
                if ( provider.equalsIgnoreCase( my_tax_prov ) ) {
                    return true;
                }
            }
            return false;
        }
        else {
            return false;
        }
    }

    public static void midpointRoot( final Phylogeny phylogeny ) {
        if ( ( phylogeny.getNumberOfExternalNodes() < 2 ) || ( calculateMaxDistanceToRoot( phylogeny ) <= 0 ) ) {
            return;
        }
        int counter = 0;
        final int total_nodes = phylogeny.getNodeCount();
        while ( true ) {
            if ( ++counter > total_nodes ) {
                throw new RuntimeException( "this should not have happened: midpoint rooting does not converge" );
            }
            PhylogenyNode a = null;
            double da = 0;
            double db = 0;
            for( int i = 0; i < phylogeny.getRoot().getNumberOfDescendants(); ++i ) {
                final PhylogenyNode f = getFurthestDescendant( phylogeny.getRoot().getChildNode( i ) );
                final double df = getDistance( f, phylogeny.getRoot() );
                if ( df > 0 ) {
                    if ( df > da ) {
                        db = da;
                        da = df;
                        a = f;
                    }
                    else if ( df > db ) {
                        db = df;
                    }
                }
            }
            final double diff = da - db;
            if ( diff < 0.000001 ) {
                break;
            }
            double x = da - ( diff / 2.0 );
            while ( ( x > a.getDistanceToParent() ) && !a.isRoot() ) {
                x -= ( a.getDistanceToParent() > 0 ? a.getDistanceToParent() : 0 );
                a = a.getParent();
            }
            phylogeny.reRoot( a, x );
        }
        phylogeny.recalculateNumberOfExternalDescendants( true );
    }

    public static void normalizeBootstrapValues( final Phylogeny phylogeny,
                                                 final double max_bootstrap_value,
                                                 final double max_normalized_value ) {
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isInternal() ) {
                final double confidence = getConfidenceValue( node );
                if ( confidence != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                    if ( confidence >= max_bootstrap_value ) {
                        setBootstrapConfidence( node, max_normalized_value );
                    }
                    else {
                        setBootstrapConfidence( node, ( confidence * max_normalized_value ) / max_bootstrap_value );
                    }
                }
            }
        }
    }

    public static List<PhylogenyNode> obtainAllNodesAsList( final Phylogeny phy ) {
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        if ( phy.isEmpty() ) {
            return nodes;
        }
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            nodes.add( iter.next() );
        }
        return nodes;
    }

    /**
     * Returns a map of distinct taxonomies of
     * all external nodes of node.
     * If at least one of the external nodes has no taxonomy,
     * null is returned.
     * 
     */
    public static Map<Taxonomy, Integer> obtainDistinctTaxonomyCounts( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        final Map<Taxonomy, Integer> tax_map = new HashMap<Taxonomy, Integer>();
        for( final PhylogenyNode n : descs ) {
            if ( !n.getNodeData().isHasTaxonomy() || n.getNodeData().getTaxonomy().isEmpty() ) {
                return null;
            }
            final Taxonomy t = n.getNodeData().getTaxonomy();
            if ( tax_map.containsKey( t ) ) {
                tax_map.put( t, tax_map.get( t ) + 1 );
            }
            else {
                tax_map.put( t, 1 );
            }
        }
        return tax_map;
    }

    /**
     * Arranges the order of childern for each node of this Phylogeny in such a
     * way that either the branch with more children is on top (right) or on
     * bottom (left), dependent on the value of boolean order.
     * 
     * @param order
     *            decides in which direction to order
     * @param pri 
     */
    public static void orderAppearance( final PhylogenyNode n,
                                        final boolean order,
                                        final boolean order_ext_alphabetically,
                                        final DESCENDANT_SORT_PRIORITY pri ) {
        if ( n.isExternal() ) {
            return;
        }
        else {
            PhylogenyNode temp = null;
            if ( ( n.getNumberOfDescendants() == 2 )
                    && ( n.getChildNode1().getNumberOfExternalNodes() != n.getChildNode2().getNumberOfExternalNodes() )
                    && ( ( n.getChildNode1().getNumberOfExternalNodes() < n.getChildNode2().getNumberOfExternalNodes() ) == order ) ) {
                temp = n.getChildNode1();
                n.setChild1( n.getChildNode2() );
                n.setChild2( temp );
            }
            else if ( order_ext_alphabetically ) {
                boolean all_ext = true;
                for( final PhylogenyNode i : n.getDescendants() ) {
                    if ( !i.isExternal() ) {
                        all_ext = false;
                        break;
                    }
                }
                if ( all_ext ) {
                    PhylogenyMethods.sortNodeDescendents( n, pri );
                }
            }
            for( int i = 0; i < n.getNumberOfDescendants(); ++i ) {
                orderAppearance( n.getChildNode( i ), order, order_ext_alphabetically, pri );
            }
        }
    }

    public static void postorderBranchColorAveragingExternalNodeBased( final Phylogeny p ) {
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            double red = 0.0;
            double green = 0.0;
            double blue = 0.0;
            int n = 0;
            if ( node.isInternal() ) {
                //for( final PhylogenyNodeIterator iterator = node.iterateChildNodesForward(); iterator.hasNext(); ) {
                for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
                    final PhylogenyNode child_node = node.getChildNode( i );
                    final Color child_color = getBranchColorValue( child_node );
                    if ( child_color != null ) {
                        ++n;
                        red += child_color.getRed();
                        green += child_color.getGreen();
                        blue += child_color.getBlue();
                    }
                }
                setBranchColorValue( node,
                                     new Color( ForesterUtil.roundToInt( red / n ),
                                                ForesterUtil.roundToInt( green / n ),
                                                ForesterUtil.roundToInt( blue / n ) ) );
            }
        }
    }

    public static final void preOrderReId( final Phylogeny phy ) {
        if ( phy.isEmpty() ) {
            return;
        }
        phy.setIdToNodeMap( null );
        long i = PhylogenyNode.getNodeCount();
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            it.next().setId( i++ );
        }
        PhylogenyNode.setNodeCount( i );
    }

    public final static Phylogeny[] readPhylogenies( final PhylogenyParser parser, final File file ) throws IOException {
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final Phylogeny[] trees = factory.create( file, parser );
        if ( ( trees == null ) || ( trees.length == 0 ) ) {
            throw new PhylogenyParserException( "Unable to parse phylogeny from file: " + file );
        }
        return trees;
    }

    public final static Phylogeny[] readPhylogenies( final PhylogenyParser parser, final List<File> files )
            throws IOException {
        final List<Phylogeny> tree_list = new ArrayList<Phylogeny>();
        for( final File file : files ) {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny[] trees = factory.create( file, parser );
            if ( ( trees == null ) || ( trees.length == 0 ) ) {
                throw new PhylogenyParserException( "Unable to parse phylogeny from file: " + file );
            }
            tree_list.addAll( Arrays.asList( trees ) );
        }
        return tree_list.toArray( new Phylogeny[ tree_list.size() ] );
    }

    public static void removeNode( final PhylogenyNode remove_me, final Phylogeny phylogeny ) {
        if ( remove_me.isRoot() ) {
            if ( remove_me.getNumberOfDescendants() == 1 ) {
                final PhylogenyNode desc = remove_me.getDescendants().get( 0 );
                desc.setDistanceToParent( addPhylogenyDistances( remove_me.getDistanceToParent(),
                                                                 desc.getDistanceToParent() ) );
                desc.setParent( null );
                phylogeny.setRoot( desc );
                phylogeny.clearHashIdToNodeMap();
            }
            else {
                throw new IllegalArgumentException( "attempt to remove a root node with more than one descendants" );
            }
        }
        else if ( remove_me.isExternal() ) {
            phylogeny.deleteSubtree( remove_me, false );
            phylogeny.clearHashIdToNodeMap();
            phylogeny.externalNodesHaveChanged();
        }
        else {
            final PhylogenyNode parent = remove_me.getParent();
            final List<PhylogenyNode> descs = remove_me.getDescendants();
            parent.removeChildNode( remove_me );
            for( final PhylogenyNode desc : descs ) {
                parent.addAsChild( desc );
                desc.setDistanceToParent( addPhylogenyDistances( remove_me.getDistanceToParent(),
                                                                 desc.getDistanceToParent() ) );
            }
            remove_me.setParent( null );
            phylogeny.clearHashIdToNodeMap();
            phylogeny.externalNodesHaveChanged();
        }
    }

    public static List<PhylogenyNode> searchData( final String query,
                                                  final Phylogeny phy,
                                                  final boolean case_sensitive,
                                                  final boolean partial,
                                                  final boolean search_domains ) {
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        if ( phy.isEmpty() || ( query == null ) ) {
            return nodes;
        }
        if ( ForesterUtil.isEmpty( query ) ) {
            return nodes;
        }
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            boolean match = false;
            if ( match( node.getName(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasTaxonomy()
                    && match( node.getNodeData().getTaxonomy().getTaxonomyCode(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasTaxonomy()
                    && match( node.getNodeData().getTaxonomy().getCommonName(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasTaxonomy()
                    && match( node.getNodeData().getTaxonomy().getScientificName(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasTaxonomy()
                    && ( node.getNodeData().getTaxonomy().getIdentifier() != null )
                    && match( node.getNodeData().getTaxonomy().getIdentifier().getValue(),
                              query,
                              case_sensitive,
                              partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasTaxonomy() && !node.getNodeData().getTaxonomy().getSynonyms().isEmpty() ) {
                final List<String> syns = node.getNodeData().getTaxonomy().getSynonyms();
                I: for( final String syn : syns ) {
                    if ( match( syn, query, case_sensitive, partial ) ) {
                        match = true;
                        break I;
                    }
                }
            }
            if ( !match && node.getNodeData().isHasSequence()
                    && match( node.getNodeData().getSequence().getName(), query, case_sensitive, partial ) ) {
                match = true;
            }
            if ( !match && node.getNodeData().isHasSequence()
                    && match( node.getNodeData().getSequence().getGeneName(), query, case_sensitive, partial ) ) {
                match = true;
            }
            if ( !match && node.getNodeData().isHasSequence()
                    && match( node.getNodeData().getSequence().getSymbol(), query, case_sensitive, partial ) ) {
                match = true;
            }
            if ( !match
                    && node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getAccession() != null )
                    && match( node.getNodeData().getSequence().getAccession().getValue(),
                              query,
                              case_sensitive,
                              partial ) ) {
                match = true;
            }
            if ( search_domains && !match && node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                final DomainArchitecture da = node.getNodeData().getSequence().getDomainArchitecture();
                I: for( int i = 0; i < da.getNumberOfDomains(); ++i ) {
                    if ( match( da.getDomain( i ).getName(), query, case_sensitive, partial ) ) {
                        match = true;
                        break I;
                    }
                }
            }
            //
            if ( !match && node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getAnnotations() != null ) ) {
                for( final Annotation ann : node.getNodeData().getSequence().getAnnotations() ) {
                    if ( match( ann.getDesc(), query, case_sensitive, partial ) ) {
                        match = true;
                        break;
                    }
                    if ( match( ann.getRef(), query, case_sensitive, partial ) ) {
                        match = true;
                        break;
                    }
                }
            }
            if ( !match && node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getCrossReferences() != null ) ) {
                for( final Accession x : node.getNodeData().getSequence().getCrossReferences() ) {
                    if ( match( x.getComment(), query, case_sensitive, partial ) ) {
                        match = true;
                        break;
                    }
                    if ( match( x.getSource(), query, case_sensitive, partial ) ) {
                        match = true;
                        break;
                    }
                    if ( match( x.getValue(), query, case_sensitive, partial ) ) {
                        match = true;
                        break;
                    }
                }
            }
            //
            if ( !match && ( node.getNodeData().getBinaryCharacters() != null ) ) {
                Iterator<String> it = node.getNodeData().getBinaryCharacters().getPresentCharacters().iterator();
                I: while ( it.hasNext() ) {
                    if ( match( it.next(), query, case_sensitive, partial ) ) {
                        match = true;
                        break I;
                    }
                }
                it = node.getNodeData().getBinaryCharacters().getGainedCharacters().iterator();
                I: while ( it.hasNext() ) {
                    if ( match( it.next(), query, case_sensitive, partial ) ) {
                        match = true;
                        break I;
                    }
                }
            }
            if ( match ) {
                nodes.add( node );
            }
        }
        return nodes;
    }

    public static List<PhylogenyNode> searchDataLogicalAnd( final String[] queries,
                                                            final Phylogeny phy,
                                                            final boolean case_sensitive,
                                                            final boolean partial,
                                                            final boolean search_domains ) {
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        if ( phy.isEmpty() || ( queries == null ) || ( queries.length < 1 ) ) {
            return nodes;
        }
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            boolean all_matched = true;
            for( final String query : queries ) {
                boolean match = false;
                if ( ForesterUtil.isEmpty( query ) ) {
                    continue;
                }
                if ( match( node.getName(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasTaxonomy()
                        && match( node.getNodeData().getTaxonomy().getTaxonomyCode(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasTaxonomy()
                        && match( node.getNodeData().getTaxonomy().getCommonName(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasTaxonomy()
                        && match( node.getNodeData().getTaxonomy().getScientificName(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasTaxonomy()
                        && ( node.getNodeData().getTaxonomy().getIdentifier() != null )
                        && match( node.getNodeData().getTaxonomy().getIdentifier().getValue(),
                                  query,
                                  case_sensitive,
                                  partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasTaxonomy()
                        && !node.getNodeData().getTaxonomy().getSynonyms().isEmpty() ) {
                    final List<String> syns = node.getNodeData().getTaxonomy().getSynonyms();
                    I: for( final String syn : syns ) {
                        if ( match( syn, query, case_sensitive, partial ) ) {
                            match = true;
                            break I;
                        }
                    }
                }
                if ( !match && node.getNodeData().isHasSequence()
                        && match( node.getNodeData().getSequence().getName(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                if ( !match && node.getNodeData().isHasSequence()
                        && match( node.getNodeData().getSequence().getGeneName(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                if ( !match && node.getNodeData().isHasSequence()
                        && match( node.getNodeData().getSequence().getSymbol(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                if ( !match
                        && node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getAccession() != null )
                        && match( node.getNodeData().getSequence().getAccession().getValue(),
                                  query,
                                  case_sensitive,
                                  partial ) ) {
                    match = true;
                }
                if ( search_domains && !match && node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                    final DomainArchitecture da = node.getNodeData().getSequence().getDomainArchitecture();
                    I: for( int i = 0; i < da.getNumberOfDomains(); ++i ) {
                        if ( match( da.getDomain( i ).getName(), query, case_sensitive, partial ) ) {
                            match = true;
                            break I;
                        }
                    }
                }
                //
                if ( !match && node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getAnnotations() != null ) ) {
                    for( final Annotation ann : node.getNodeData().getSequence().getAnnotations() ) {
                        if ( match( ann.getDesc(), query, case_sensitive, partial ) ) {
                            match = true;
                            break;
                        }
                        if ( match( ann.getRef(), query, case_sensitive, partial ) ) {
                            match = true;
                            break;
                        }
                    }
                }
                if ( !match && node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getCrossReferences() != null ) ) {
                    for( final Accession x : node.getNodeData().getSequence().getCrossReferences() ) {
                        if ( match( x.getComment(), query, case_sensitive, partial ) ) {
                            match = true;
                            break;
                        }
                        if ( match( x.getSource(), query, case_sensitive, partial ) ) {
                            match = true;
                            break;
                        }
                        if ( match( x.getValue(), query, case_sensitive, partial ) ) {
                            match = true;
                            break;
                        }
                    }
                }
                //
                if ( !match && ( node.getNodeData().getBinaryCharacters() != null ) ) {
                    Iterator<String> it = node.getNodeData().getBinaryCharacters().getPresentCharacters().iterator();
                    I: while ( it.hasNext() ) {
                        if ( match( it.next(), query, case_sensitive, partial ) ) {
                            match = true;
                            break I;
                        }
                    }
                    it = node.getNodeData().getBinaryCharacters().getGainedCharacters().iterator();
                    I: while ( it.hasNext() ) {
                        if ( match( it.next(), query, case_sensitive, partial ) ) {
                            match = true;
                            break I;
                        }
                    }
                }
                if ( !match ) {
                    all_matched = false;
                    break;
                }
            }
            if ( all_matched ) {
                nodes.add( node );
            }
        }
        return nodes;
    }

    /**
     * Convenience method.
     * Sets value for the first confidence value (created if not present, values overwritten otherwise). 
     */
    public static void setBootstrapConfidence( final PhylogenyNode node, final double bootstrap_confidence_value ) {
        setConfidence( node, bootstrap_confidence_value, "bootstrap" );
    }

    public static void setBranchColorValue( final PhylogenyNode node, final Color color ) {
        if ( node.getBranchData().getBranchColor() == null ) {
            node.getBranchData().setBranchColor( new BranchColor() );
        }
        node.getBranchData().getBranchColor().setValue( color );
    }

    /**
     * Convenience method
     */
    public static void setBranchWidthValue( final PhylogenyNode node, final double branch_width_value ) {
        node.getBranchData().setBranchWidth( new BranchWidth( branch_width_value ) );
    }

    /**
     * Convenience method.
     * Sets value for the first confidence value (created if not present, values overwritten otherwise). 
     */
    public static void setConfidence( final PhylogenyNode node, final double confidence_value ) {
        setConfidence( node, confidence_value, "" );
    }

    /**
     * Convenience method.
     * Sets value for the first confidence value (created if not present, values overwritten otherwise). 
     */
    public static void setConfidence( final PhylogenyNode node, final double confidence_value, final String type ) {
        Confidence c = null;
        if ( node.getBranchData().getNumberOfConfidences() > 0 ) {
            c = node.getBranchData().getConfidence( 0 );
        }
        else {
            c = new Confidence();
            node.getBranchData().addConfidence( c );
        }
        c.setType( type );
        c.setValue( confidence_value );
    }

    public static void setScientificName( final PhylogenyNode node, final String scientific_name ) {
        if ( !node.getNodeData().isHasTaxonomy() ) {
            node.getNodeData().setTaxonomy( new Taxonomy() );
        }
        node.getNodeData().getTaxonomy().setScientificName( scientific_name );
    }

    /**
     * Convenience method to set the taxonomy code of a phylogeny node.
     * 
     * 
     * @param node
     * @param taxonomy_code
     * @throws PhyloXmlDataFormatException 
     */
    public static void setTaxonomyCode( final PhylogenyNode node, final String taxonomy_code )
            throws PhyloXmlDataFormatException {
        if ( !node.getNodeData().isHasTaxonomy() ) {
            node.getNodeData().setTaxonomy( new Taxonomy() );
        }
        node.getNodeData().getTaxonomy().setTaxonomyCode( taxonomy_code );
    }

    final static public void sortNodeDescendents( final PhylogenyNode node, final DESCENDANT_SORT_PRIORITY pri ) {
        class PhylogenyNodeSortTaxonomyPriority implements Comparator<PhylogenyNode> {

            @Override
            public int compare( final PhylogenyNode n1, final PhylogenyNode n2 ) {
                if ( n1.getNodeData().isHasTaxonomy() && n2.getNodeData().isHasTaxonomy() ) {
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getTaxonomy().getScientificName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getTaxonomy().getScientificName() ) ) ) {
                        return n1.getNodeData().getTaxonomy().getScientificName().toLowerCase()
                                .compareTo( n2.getNodeData().getTaxonomy().getScientificName().toLowerCase() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getTaxonomy().getTaxonomyCode() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getTaxonomy().getTaxonomyCode() ) ) ) {
                        return n1.getNodeData().getTaxonomy().getTaxonomyCode()
                                .compareTo( n2.getNodeData().getTaxonomy().getTaxonomyCode() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getTaxonomy().getCommonName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getTaxonomy().getCommonName() ) ) ) {
                        return n1.getNodeData().getTaxonomy().getCommonName().toLowerCase()
                                .compareTo( n2.getNodeData().getTaxonomy().getCommonName().toLowerCase() );
                    }
                }
                if ( n1.getNodeData().isHasSequence() && n2.getNodeData().isHasSequence() ) {
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getName() ) ) ) {
                        return n1.getNodeData().getSequence().getName().toLowerCase()
                                .compareTo( n2.getNodeData().getSequence().getName().toLowerCase() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getSymbol() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getSymbol() ) ) ) {
                        return n1.getNodeData().getSequence().getSymbol()
                                .compareTo( n2.getNodeData().getSequence().getSymbol() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getGeneName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getGeneName() ) ) ) {
                        return n1.getNodeData().getSequence().getGeneName()
                                .compareTo( n2.getNodeData().getSequence().getGeneName() );
                    }
                    if ( ( n1.getNodeData().getSequence().getAccession() != null )
                            && ( n2.getNodeData().getSequence().getAccession() != null )
                            && !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getAccession().getValue() )
                            && !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getAccession().getValue() ) ) {
                        return n1.getNodeData().getSequence().getAccession().getValue()
                                .compareTo( n2.getNodeData().getSequence().getAccession().getValue() );
                    }
                }
                if ( ( !ForesterUtil.isEmpty( n1.getName() ) ) && ( !ForesterUtil.isEmpty( n2.getName() ) ) ) {
                    return n1.getName().toLowerCase().compareTo( n2.getName().toLowerCase() );
                }
                return 0;
            }
        }
        class PhylogenyNodeSortSequencePriority implements Comparator<PhylogenyNode> {

            @Override
            public int compare( final PhylogenyNode n1, final PhylogenyNode n2 ) {
                if ( n1.getNodeData().isHasSequence() && n2.getNodeData().isHasSequence() ) {
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getName() ) ) ) {
                        return n1.getNodeData().getSequence().getName().toLowerCase()
                                .compareTo( n2.getNodeData().getSequence().getName().toLowerCase() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getSymbol() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getSymbol() ) ) ) {
                        return n1.getNodeData().getSequence().getSymbol()
                                .compareTo( n2.getNodeData().getSequence().getSymbol() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getGeneName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getGeneName() ) ) ) {
                        return n1.getNodeData().getSequence().getGeneName()
                                .compareTo( n2.getNodeData().getSequence().getGeneName() );
                    }
                    if ( ( n1.getNodeData().getSequence().getAccession() != null )
                            && ( n2.getNodeData().getSequence().getAccession() != null )
                            && !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getAccession().getValue() )
                            && !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getAccession().getValue() ) ) {
                        return n1.getNodeData().getSequence().getAccession().getValue()
                                .compareTo( n2.getNodeData().getSequence().getAccession().getValue() );
                    }
                }
                if ( n1.getNodeData().isHasTaxonomy() && n2.getNodeData().isHasTaxonomy() ) {
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getTaxonomy().getScientificName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getTaxonomy().getScientificName() ) ) ) {
                        return n1.getNodeData().getTaxonomy().getScientificName().toLowerCase()
                                .compareTo( n2.getNodeData().getTaxonomy().getScientificName().toLowerCase() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getTaxonomy().getTaxonomyCode() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getTaxonomy().getTaxonomyCode() ) ) ) {
                        return n1.getNodeData().getTaxonomy().getTaxonomyCode()
                                .compareTo( n2.getNodeData().getTaxonomy().getTaxonomyCode() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getTaxonomy().getCommonName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getTaxonomy().getCommonName() ) ) ) {
                        return n1.getNodeData().getTaxonomy().getCommonName().toLowerCase()
                                .compareTo( n2.getNodeData().getTaxonomy().getCommonName().toLowerCase() );
                    }
                }
                if ( ( !ForesterUtil.isEmpty( n1.getName() ) ) && ( !ForesterUtil.isEmpty( n2.getName() ) ) ) {
                    return n1.getName().toLowerCase().compareTo( n2.getName().toLowerCase() );
                }
                return 0;
            }
        }
        class PhylogenyNodeSortNodeNamePriority implements Comparator<PhylogenyNode> {

            @Override
            public int compare( final PhylogenyNode n1, final PhylogenyNode n2 ) {
                if ( ( !ForesterUtil.isEmpty( n1.getName() ) ) && ( !ForesterUtil.isEmpty( n2.getName() ) ) ) {
                    return n1.getName().toLowerCase().compareTo( n2.getName().toLowerCase() );
                }
                if ( n1.getNodeData().isHasTaxonomy() && n2.getNodeData().isHasTaxonomy() ) {
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getTaxonomy().getScientificName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getTaxonomy().getScientificName() ) ) ) {
                        return n1.getNodeData().getTaxonomy().getScientificName().toLowerCase()
                                .compareTo( n2.getNodeData().getTaxonomy().getScientificName().toLowerCase() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getTaxonomy().getTaxonomyCode() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getTaxonomy().getTaxonomyCode() ) ) ) {
                        return n1.getNodeData().getTaxonomy().getTaxonomyCode()
                                .compareTo( n2.getNodeData().getTaxonomy().getTaxonomyCode() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getTaxonomy().getCommonName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getTaxonomy().getCommonName() ) ) ) {
                        return n1.getNodeData().getTaxonomy().getCommonName().toLowerCase()
                                .compareTo( n2.getNodeData().getTaxonomy().getCommonName().toLowerCase() );
                    }
                }
                if ( n1.getNodeData().isHasSequence() && n2.getNodeData().isHasSequence() ) {
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getName() ) ) ) {
                        return n1.getNodeData().getSequence().getName().toLowerCase()
                                .compareTo( n2.getNodeData().getSequence().getName().toLowerCase() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getSymbol() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getSymbol() ) ) ) {
                        return n1.getNodeData().getSequence().getSymbol()
                                .compareTo( n2.getNodeData().getSequence().getSymbol() );
                    }
                    if ( ( !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getGeneName() ) )
                            && ( !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getGeneName() ) ) ) {
                        return n1.getNodeData().getSequence().getGeneName()
                                .compareTo( n2.getNodeData().getSequence().getGeneName() );
                    }
                    if ( ( n1.getNodeData().getSequence().getAccession() != null )
                            && ( n2.getNodeData().getSequence().getAccession() != null )
                            && !ForesterUtil.isEmpty( n1.getNodeData().getSequence().getAccession().getValue() )
                            && !ForesterUtil.isEmpty( n2.getNodeData().getSequence().getAccession().getValue() ) ) {
                        return n1.getNodeData().getSequence().getAccession().getValue()
                                .compareTo( n2.getNodeData().getSequence().getAccession().getValue() );
                    }
                }
                return 0;
            }
        }
        Comparator<PhylogenyNode> c;
        switch ( pri ) {
            case SEQUENCE:
                c = new PhylogenyNodeSortSequencePriority();
                break;
            case NODE_NAME:
                c = new PhylogenyNodeSortNodeNamePriority();
                break;
            default:
                c = new PhylogenyNodeSortTaxonomyPriority();
        }
        final List<PhylogenyNode> descs = node.getDescendants();
        Collections.sort( descs, c );
        int i = 0;
        for( final PhylogenyNode desc : descs ) {
            node.setChildNode( i++, desc );
        }
    }

    /**
     * Removes from Phylogeny to_be_stripped all external Nodes which are
     * associated with a species NOT found in Phylogeny reference.
     * 
     * @param reference
     *            a reference Phylogeny
     * @param to_be_stripped
     *            Phylogeny to be stripped
     * @return nodes removed from to_be_stripped
     */
    public static List<PhylogenyNode> taxonomyBasedDeletionOfExternalNodes( final Phylogeny reference,
                                                                            final Phylogeny to_be_stripped ) {
        final Set<String> ref_ext_taxo = new HashSet<String>();
        for( final PhylogenyNodeIterator it = reference.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.getNodeData().isHasTaxonomy() ) {
                throw new IllegalArgumentException( "no taxonomic data in node: " + n );
            }
            if ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                ref_ext_taxo.add( n.getNodeData().getTaxonomy().getScientificName() );
            }
            if ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
                ref_ext_taxo.add( n.getNodeData().getTaxonomy().getTaxonomyCode() );
            }
            if ( ( n.getNodeData().getTaxonomy().getIdentifier() != null )
                    && !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getIdentifier().getValue() ) ) {
                ref_ext_taxo.add( n.getNodeData().getTaxonomy().getIdentifier().getValuePlusProvider() );
            }
        }
        final ArrayList<PhylogenyNode> nodes_to_delete = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator it = to_be_stripped.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.getNodeData().isHasTaxonomy() ) {
                nodes_to_delete.add( n );
            }
            else if ( !( ref_ext_taxo.contains( n.getNodeData().getTaxonomy().getScientificName() ) )
                    && !( ref_ext_taxo.contains( n.getNodeData().getTaxonomy().getTaxonomyCode() ) )
                    && !( ( n.getNodeData().getTaxonomy().getIdentifier() != null ) && ref_ext_taxo.contains( n
                            .getNodeData().getTaxonomy().getIdentifier().getValuePlusProvider() ) ) ) {
                nodes_to_delete.add( n );
            }
        }
        for( final PhylogenyNode n : nodes_to_delete ) {
            to_be_stripped.deleteSubtree( n, true );
        }
        to_be_stripped.clearHashIdToNodeMap();
        to_be_stripped.externalNodesHaveChanged();
        return nodes_to_delete;
    }

    final static public void transferInternalNamesToBootstrapSupport( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !ForesterUtil.isEmpty( n.getName() ) ) {
                double value = -1;
                try {
                    value = Double.parseDouble( n.getName() );
                }
                catch ( final NumberFormatException e ) {
                    throw new IllegalArgumentException( "failed to parse number from [" + n.getName() + "]: "
                            + e.getLocalizedMessage() );
                }
                if ( value >= 0.0 ) {
                    n.getBranchData().addConfidence( new Confidence( value, "bootstrap" ) );
                    n.setName( "" );
                }
            }
        }
    }

    final static public void transferInternalNodeNamesToConfidence( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.getBranchData().isHasConfidences() ) {
                if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                    double d = -1.0;
                    try {
                        d = Double.parseDouble( n.getName() );
                    }
                    catch ( final Exception e ) {
                        d = -1.0;
                    }
                    if ( d >= 0.0 ) {
                        n.getBranchData().addConfidence( new Confidence( d, "" ) );
                        n.setName( "" );
                    }
                }
            }
        }
    }

    final static public void transferNodeNameToField( final Phylogeny phy,
                                                      final PhylogenyNodeField field,
                                                      final boolean external_only ) throws PhyloXmlDataFormatException {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( external_only && n.isInternal() ) {
                continue;
            }
            final String name = n.getName().trim();
            if ( !ForesterUtil.isEmpty( name ) ) {
                switch ( field ) {
                    case TAXONOMY_CODE:
                        n.setName( "" );
                        setTaxonomyCode( n, name );
                        break;
                    case TAXONOMY_SCIENTIFIC_NAME:
                        n.setName( "" );
                        if ( !n.getNodeData().isHasTaxonomy() ) {
                            n.getNodeData().setTaxonomy( new Taxonomy() );
                        }
                        n.getNodeData().getTaxonomy().setScientificName( name );
                        break;
                    case TAXONOMY_COMMON_NAME:
                        n.setName( "" );
                        if ( !n.getNodeData().isHasTaxonomy() ) {
                            n.getNodeData().setTaxonomy( new Taxonomy() );
                        }
                        n.getNodeData().getTaxonomy().setCommonName( name );
                        break;
                    case SEQUENCE_SYMBOL:
                        n.setName( "" );
                        if ( !n.getNodeData().isHasSequence() ) {
                            n.getNodeData().setSequence( new Sequence() );
                        }
                        n.getNodeData().getSequence().setSymbol( name );
                        break;
                    case SEQUENCE_NAME:
                        n.setName( "" );
                        if ( !n.getNodeData().isHasSequence() ) {
                            n.getNodeData().setSequence( new Sequence() );
                        }
                        n.getNodeData().getSequence().setName( name );
                        break;
                    case TAXONOMY_ID_UNIPROT_1: {
                        if ( !n.getNodeData().isHasTaxonomy() ) {
                            n.getNodeData().setTaxonomy( new Taxonomy() );
                        }
                        String id = name;
                        final int i = name.indexOf( '_' );
                        if ( i > 0 ) {
                            id = name.substring( 0, i );
                        }
                        else {
                            n.setName( "" );
                        }
                        n.getNodeData().getTaxonomy()
                                .setIdentifier( new Identifier( id, PhyloXmlUtil.UNIPROT_TAX_PROVIDER ) );
                        break;
                    }
                    case TAXONOMY_ID_UNIPROT_2: {
                        if ( !n.getNodeData().isHasTaxonomy() ) {
                            n.getNodeData().setTaxonomy( new Taxonomy() );
                        }
                        String id = name;
                        final int i = name.indexOf( '_' );
                        if ( i > 0 ) {
                            id = name.substring( i + 1, name.length() );
                        }
                        else {
                            n.setName( "" );
                        }
                        n.getNodeData().getTaxonomy()
                                .setIdentifier( new Identifier( id, PhyloXmlUtil.UNIPROT_TAX_PROVIDER ) );
                        break;
                    }
                    case TAXONOMY_ID: {
                        if ( !n.getNodeData().isHasTaxonomy() ) {
                            n.getNodeData().setTaxonomy( new Taxonomy() );
                        }
                        n.getNodeData().getTaxonomy().setIdentifier( new Identifier( name ) );
                        break;
                    }
                }
            }
        }
    }

    static double addPhylogenyDistances( final double a, final double b ) {
        if ( ( a >= 0.0 ) && ( b >= 0.0 ) ) {
            return a + b;
        }
        else if ( a >= 0.0 ) {
            return a;
        }
        else if ( b >= 0.0 ) {
            return b;
        }
        return PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT;
    }

    static double calculateDistanceToAncestor( final PhylogenyNode anc, PhylogenyNode desc ) {
        double d = 0;
        boolean all_default = true;
        while ( anc != desc ) {
            if ( desc.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                d += desc.getDistanceToParent();
                if ( all_default ) {
                    all_default = false;
                }
            }
            desc = desc.getParent();
        }
        if ( all_default ) {
            return PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT;
        }
        return d;
    }

    /**
     * Deep copies the phylogeny originating from this node.
     */
    static PhylogenyNode copySubTree( final PhylogenyNode source ) {
        if ( source == null ) {
            return null;
        }
        else {
            final PhylogenyNode newnode = source.copyNodeData();
            if ( !source.isExternal() ) {
                for( int i = 0; i < source.getNumberOfDescendants(); ++i ) {
                    newnode.setChildNode( i, PhylogenyMethods.copySubTree( source.getChildNode( i ) ) );
                }
            }
            return newnode;
        }
    }

    /**
     * Shallow copies the phylogeny originating from this node.
     */
    static PhylogenyNode copySubTreeShallow( final PhylogenyNode source ) {
        if ( source == null ) {
            return null;
        }
        else {
            final PhylogenyNode newnode = source.copyNodeDataShallow();
            if ( !source.isExternal() ) {
                for( int i = 0; i < source.getNumberOfDescendants(); ++i ) {
                    newnode.setChildNode( i, PhylogenyMethods.copySubTreeShallow( source.getChildNode( i ) ) );
                }
            }
            return newnode;
        }
    }

    /**
     * Calculates the distance between PhylogenyNodes n1 and n2.
     * PRECONDITION: n1 is a descendant of n2.
     * 
     * @param n1
     *            a descendant of n2
     * @param n2
     * @return distance between n1 and n2
     */
    private static double getDistance( PhylogenyNode n1, final PhylogenyNode n2 ) {
        double d = 0.0;
        while ( n1 != n2 ) {
            if ( n1.getDistanceToParent() > 0.0 ) {
                d += n1.getDistanceToParent();
            }
            n1 = n1.getParent();
        }
        return d;
    }

    private static boolean match( final String s,
                                  final String query,
                                  final boolean case_sensitive,
                                  final boolean partial ) {
        if ( ForesterUtil.isEmpty( s ) || ForesterUtil.isEmpty( query ) ) {
            return false;
        }
        String my_s = s.trim();
        String my_query = query.trim();
        if ( !case_sensitive ) {
            my_s = my_s.toLowerCase();
            my_query = my_query.toLowerCase();
        }
        if ( partial ) {
            return my_s.indexOf( my_query ) >= 0;
        }
        else {
            return my_s.equals( my_query );
        }
    }

    public static enum DESCENDANT_SORT_PRIORITY {
        TAXONOMY, SEQUENCE, NODE_NAME;
    }

    public static enum PhylogenyNodeField {
        CLADE_NAME,
        TAXONOMY_CODE,
        TAXONOMY_SCIENTIFIC_NAME,
        TAXONOMY_COMMON_NAME,
        SEQUENCE_SYMBOL,
        SEQUENCE_NAME,
        TAXONOMY_ID_UNIPROT_1,
        TAXONOMY_ID_UNIPROT_2,
        TAXONOMY_ID;
    }
}
