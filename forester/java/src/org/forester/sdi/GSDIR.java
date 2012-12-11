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
import org.forester.phylogeny.PhylogenyBranch;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicDescriptiveStatistics;

public class GSDIR extends GSDI {

    private int                              _min_duplications_sum;
    private final BasicDescriptiveStatistics _duplications_sum_stats;
    private final List<Phylogeny>            _min_duplications_sum_gene_trees;

    public GSDIR( final Phylogeny gene_tree, final Phylogeny species_tree, final boolean strip_gene_tree, final int x )
            throws SDIException {
        super( gene_tree.copy(), species_tree, true, strip_gene_tree, true, 1 );
        _min_duplications_sum = Integer.MAX_VALUE;
        _min_duplications_sum_gene_trees = new ArrayList<Phylogeny>();
        _duplications_sum_stats = new BasicDescriptiveStatistics();
        linkNodesOfG();
        final List<PhylogenyBranch> gene_tree_branches_post_order = new ArrayList<PhylogenyBranch>();
        for( final PhylogenyNodeIterator it = _gene_tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isRoot() ) {
                gene_tree_branches_post_order.add( new PhylogenyBranch( n, n.getParent() ) );
            }
        }
        for( final PhylogenyBranch branch : gene_tree_branches_post_order ) {
            _duplications_sum = 0;
            _speciation_or_duplication_events_sum = 0;
            _speciations_sum = 0;
            _gene_tree.reRoot( branch );
            PhylogenyMethods.preOrderReId( getSpeciesTree() );
            //TEST, remove later
            for( final PhylogenyNodeIterator it = _gene_tree.iteratorPostorder(); it.hasNext(); ) {
                final PhylogenyNode g = it.next();
                if ( g.isInternal() ) {
                    g.setLink( null );
                }
            }
            geneTreePostOrderTraversal();
            if ( _duplications_sum < _min_duplications_sum ) {
                _min_duplications_sum = _duplications_sum;
                _min_duplications_sum_gene_trees.clear();
                _min_duplications_sum_gene_trees.add( getGeneTree().copy() );
            }
            else if ( _duplications_sum == _min_duplications_sum ) {
                _min_duplications_sum_gene_trees.add( getGeneTree().copy() );
            }
            _duplications_sum_stats.addValue( _duplications_sum );
        }
        System.out.println( _duplications_sum_stats.getSummaryAsString() );
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
