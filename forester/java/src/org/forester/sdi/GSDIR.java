// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2013 Christian M. Zmasek
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
// WWW: www.phylosoft.org

package org.forester.sdi;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyBranch;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sdi.SDIutil.TaxonomyComparisonBase;
import org.forester.util.BasicDescriptiveStatistics;

public class GSDIR implements GSDII {

    private final int                        _min_duplications_sum;
    private final int                        _speciations_sum;
    private final BasicDescriptiveStatistics _duplications_sum_stats;
    private Phylogeny                        _min_duplications_sum_gene_tree;
    private final List<PhylogenyNode>        _stripped_gene_tree_nodes;
    private final List<PhylogenyNode>        _stripped_species_tree_nodes;
    private final Set<PhylogenyNode>         _mapped_species_tree_nodes;
    private final TaxonomyComparisonBase     _tax_comp_base;
    private final SortedSet<String>          _scientific_names_mapped_to_reduced_specificity;

    public GSDIR( final Phylogeny gene_tree,
                  final Phylogeny species_tree,
                  final boolean strip_gene_tree,
                  final boolean strip_species_tree ) throws SDIException {
        final NodesLinkingResult nodes_linking_result = GSDI.linkNodesOfG( gene_tree,
                                                                           species_tree,
                                                                           strip_gene_tree,
                                                                           strip_species_tree );
        _stripped_gene_tree_nodes = nodes_linking_result.getStrippedGeneTreeNodes();
        _stripped_species_tree_nodes = nodes_linking_result.getStrippedSpeciesTreeNodes();
        _mapped_species_tree_nodes = nodes_linking_result.getMappedSpeciesTreeNodes();
        _scientific_names_mapped_to_reduced_specificity = nodes_linking_result
                .getScientificNamesMappedToReducedSpecificity();
        _tax_comp_base = nodes_linking_result.getTaxCompBase();
        final List<PhylogenyBranch> gene_tree_branches_post_order = new ArrayList<PhylogenyBranch>();
        for( final PhylogenyNodeIterator it = gene_tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isRoot() /*&& !( n.getParent().isRoot() && n.isFirstChildNode() )*/) {
                gene_tree_branches_post_order.add( new PhylogenyBranch( n, n.getParent() ) );
            }
        }
        int min_duplications_sum = Integer.MAX_VALUE;
        int speciations_sum = 0;
        _duplications_sum_stats = new BasicDescriptiveStatistics();
        for( final PhylogenyBranch branch : gene_tree_branches_post_order ) {
            gene_tree.reRoot( branch );
            PhylogenyMethods.preOrderReId( species_tree );
            //TEST, remove later
            //            for( final PhylogenyNodeIterator it = _gene_tree.iteratorPostorder(); it.hasNext(); ) {
            //                final PhylogenyNode g = it.next();
            //                if ( g.isInternal() ) {
            //                    g.setLink( null );
            //                }
            //            }
            final GSDIsummaryResult gsdi_result = GSDI.geneTreePostOrderTraversal( gene_tree,
                                                                                   true,
                                                                                   min_duplications_sum );
            if ( gsdi_result == null ) {
                continue;
            }
            if ( gsdi_result.getDuplicationsSum() < min_duplications_sum ) {
                min_duplications_sum = gsdi_result.getDuplicationsSum();
                speciations_sum = gsdi_result.getSpeciationsSum();
                _min_duplications_sum_gene_tree = gene_tree.copy();
            }
            else if ( gsdi_result.getDuplicationsSum() == min_duplications_sum ) {
                final List<Phylogeny> l = new ArrayList<Phylogeny>();
                l.add( _min_duplications_sum_gene_tree );
                l.add( gene_tree );
                final int index = getIndexesOfShortestTree( l ).get( 0 );
                if ( index == 1 ) {
                    _min_duplications_sum_gene_tree = gene_tree.copy();
                }
            }
            _duplications_sum_stats.addValue( gsdi_result.getDuplicationsSum() );
        }
        _min_duplications_sum = min_duplications_sum;
        _speciations_sum = speciations_sum;
    }

    public BasicDescriptiveStatistics getDuplicationsSumStats() {
        return _duplications_sum_stats;
    }

    @Override
    public Set<PhylogenyNode> getMappedExternalSpeciesTreeNodes() {
        return _mapped_species_tree_nodes;
    }

    public int getMinDuplicationsSum() {
        return _min_duplications_sum;
    }

    public Phylogeny getMinDuplicationsSumGeneTree() {
        return _min_duplications_sum_gene_tree;
    }

    @Override
    public final SortedSet<String> getReMappedScientificNamesFromGeneTree() {
        return _scientific_names_mapped_to_reduced_specificity;
    }

    @Override
    public int getSpeciationsSum() {
        return _speciations_sum;
    }

    @Override
    public List<PhylogenyNode> getStrippedExternalGeneTreeNodes() {
        return _stripped_gene_tree_nodes;
    }

    @Override
    public List<PhylogenyNode> getStrippedSpeciesTreeNodes() {
        return _stripped_species_tree_nodes;
    }

    @Override
    public TaxonomyComparisonBase getTaxCompBase() {
        return _tax_comp_base;
    }

    public final static List<Integer> getIndexesOfShortestTree( final List<Phylogeny> assigned_trees ) {
        final List<Integer> shortests = new ArrayList<Integer>();
        boolean depth = true;
        double x = Double.MAX_VALUE;
        for( int i = 0; i < assigned_trees.size(); ++i ) {
            final Phylogeny phy = assigned_trees.get( i );
            if ( i == 0 ) {
                if ( PhylogenyMethods.calculateMaxDistanceToRoot( phy ) > 0 ) {
                    depth = false;
                }
            }
            final double d;
            if ( depth ) {
                d = PhylogenyMethods.calculateMaxDepth( phy );
            }
            else {
                d = PhylogenyMethods.calculateMaxDistanceToRoot( phy );
            }
            if ( d < x ) {
                x = d;
                shortests.clear();
                shortests.add( i );
            }
            else if ( d == x ) {
                shortests.add( i );
            }
        }
        return shortests;
    }
}
