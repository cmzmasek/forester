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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicDescriptiveStatistics;

public class GSDIR extends GSDI {

    private int                   _min_duplications_sum;
    private final BasicDescriptiveStatistics                  _duplications_sum_stats;
    private final List<Phylogeny> _min_duplications_sum_gene_trees;

    public GSDIR( final Phylogeny gene_tree, final Phylogeny species_tree, final boolean strip_gene_tree, final int x )
            throws SDIException {
        super( gene_tree, species_tree, true, strip_gene_tree, true, 1 );
        _min_duplications_sum = Integer.MAX_VALUE;
        _min_duplications_sum_gene_trees = new ArrayList<Phylogeny>();
        _duplications_sum_stats = new  BasicDescriptiveStatistics();
        linkNodesOfG();
        final List<PhylogenyNode> gene_tree_nodes_post_order = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator it = gene_tree.iteratorPostorder(); it.hasNext(); ) {
            gene_tree_nodes_post_order.add( it.next() );
        }
        for( final PhylogenyNode root : gene_tree_nodes_post_order ) {
            _gene_tree.reRoot( root );
            PhylogenyMethods.preOrderReId( getSpeciesTree() );
            //TEST, remove later
            for( final PhylogenyNodeIterator it = getGeneTree().iteratorPostorder(); it.hasNext(); ) {
                final PhylogenyNode g = it.next();
                if ( g.isInternal() ) {
                    g.setLink( null );
                }
            }
            geneTreePostOrderTraversal();
            _duplications_sum_stats.addValue( getMinDuplicationsSum() );
            if ( getDuplicationsSum() < getMinDuplicationsSum() ) {
                _min_duplications_sum = getDuplicationsSum();
                _min_duplications_sum_gene_trees.clear();
                _min_duplications_sum_gene_trees.add( getGeneTree().copy() );
            }
            else if ( getDuplicationsSum() == getMinDuplicationsSum() ) {
                _min_duplications_sum_gene_trees.add( getGeneTree().copy() );
            }
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
}
