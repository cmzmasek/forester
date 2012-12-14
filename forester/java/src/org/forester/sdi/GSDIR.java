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

public class GSDIR {

    private int                              _min_duplications_sum;
    private final BasicDescriptiveStatistics _duplications_sum_stats;
    private final List<Phylogeny>            _min_duplications_sum_gene_trees;
    protected int                            _speciations_sum;
    protected int                            _duplications_sum;
    private final List<PhylogenyNode>        _stripped_gene_tree_nodes;
    private final List<PhylogenyNode>        _stripped_species_tree_nodes;
    private final Set<PhylogenyNode>         _mapped_species_tree_nodes;
    private final TaxonomyComparisonBase     _tax_comp_base;
    private final SortedSet<String>          _scientific_names_mapped_to_reduced_specificity;

    public GSDIR( final Phylogeny gene_tree,
                  final Phylogeny species_tree,
                  final boolean strip_gene_tree,
                  final boolean strip_species_tree ) throws SDIException {
        _speciations_sum = 0;
        _duplications_sum = 0;
        final NodesLinkingResult nodes_linking_result = GSDI.linkNodesOfG( gene_tree,
                                                                           species_tree,
                                                                           null,
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
            if ( !n.isRoot() && !( n.getParent().isRoot() && n.isFirstChildNode() ) ) {
                gene_tree_branches_post_order.add( new PhylogenyBranch( n, n.getParent() ) );
            }
        }
        _min_duplications_sum = Integer.MAX_VALUE;
        _min_duplications_sum_gene_trees = new ArrayList<Phylogeny>();
        _duplications_sum_stats = new BasicDescriptiveStatistics();
        for( final PhylogenyBranch branch : gene_tree_branches_post_order ) {
            _duplications_sum = 0;
            _speciations_sum = 0;
            gene_tree.reRoot( branch );
            PhylogenyMethods.preOrderReId( species_tree );
            //TEST, remove later
            //            for( final PhylogenyNodeIterator it = _gene_tree.iteratorPostorder(); it.hasNext(); ) {
            //                final PhylogenyNode g = it.next();
            //                if ( g.isInternal() ) {
            //                    g.setLink( null );
            //                }
            //            }
            final GSDIsummaryResult gsdi_summary_result = new GSDIsummaryResult();
            GSDI.geneTreePostOrderTraversal( gene_tree, true, gsdi_summary_result );
            if ( _duplications_sum < _min_duplications_sum ) {
                _min_duplications_sum = _duplications_sum;
                _min_duplications_sum_gene_trees.clear();
                _min_duplications_sum_gene_trees.add( gene_tree.copy() );
                //_speciations_sum
            }
            else if ( _duplications_sum == _min_duplications_sum ) {
                _min_duplications_sum_gene_trees.add( gene_tree.copy() );
            }
            _duplications_sum_stats.addValue( _duplications_sum );
        }
    }

    public int getMinDuplicationsSum() {
        return _min_duplications_sum;
    }

    public List<Phylogeny> getMinDuplicationsSumGeneTrees() {
        return _min_duplications_sum_gene_trees;
    }

    public BasicDescriptiveStatistics getDuplicationsSumStats() {
        return _duplications_sum_stats;
    }

    public Set<PhylogenyNode> getMappedExternalSpeciesTreeNodes() {
        return _mapped_species_tree_nodes;
    }

    public final SortedSet<String> getReMappedScientificNamesFromGeneTree() {
        return _scientific_names_mapped_to_reduced_specificity;
    }

    public List<PhylogenyNode> getStrippedExternalGeneTreeNodes() {
        return _stripped_gene_tree_nodes;
    }

    public List<PhylogenyNode> getStrippedSpeciesTreeNodes() {
        return _stripped_species_tree_nodes;
    }

    public TaxonomyComparisonBase getTaxCompBase() {
        return _tax_comp_base;
    }
}
