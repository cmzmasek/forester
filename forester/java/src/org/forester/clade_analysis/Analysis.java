// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2017 Christian M. Zmasek
// Copyright (C) 2017 J. Craig Venter Institute
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
// Contact: phyloxml @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester
// --------------------
// TODO
// * Multiple "hits" with different "M" values
// * More tests (including multiple children per node), especially on edge cases
// * Utilize relevant support values for warnings
// * Better system for "clade label creation" (e.g. 1.3.4 + 1.3.6 -> 1.3), use
// specific separator (eg . | _ )

package org.forester.clade_analysis;

import java.util.ArrayList;
import java.util.List;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

public final class Analysis {

    public static Result execute( final Phylogeny p, final String query ) {
        final PhylogenyNode qnode = p.getNode( query );
        if ( qnode.isRoot() ) {
            throw new IllegalStateException( "Unexpected error: Query " + query
                    + " is root. This should have never happened" );
        }
        if ( qnode.getParent().isRoot() ) {
            throw new IllegalStateException( "Unexpected error: Parent of query " + query
                    + " is root. This should have never happened" );
        }
        final PhylogenyNode qnode_p = qnode.getParent();
        final PhylogenyNode qnode_pp = qnode.getParent().getParent();
        final List<PhylogenyNode> qnode_ext_nodes = qnode_pp.getAllExternalDescendants();
        final int lec_ext_nodes = qnode_ext_nodes.size() - 1;
        final int p_ext_nodes = p.getNumberOfExternalNodes() - 1;
        final List<String> qnode_ext_nodes_names = new ArrayList<>();
        for( final PhylogenyNode qnode_ext_node : qnode_ext_nodes ) {
            String name = qnode_ext_node.getName();
            if ( ForesterUtil.isEmptyTrimmed( name ) ) {
                throw new IllegalArgumentException( "external node(s) with empty names found" );
            }
            name = name.trim();
            if ( !name.equals( query ) ) {
                qnode_ext_nodes_names.add( name );
            }
        }
        final String greatest_common_prefix = ForesterUtil.greatestCommonPrefix( qnode_ext_nodes_names );
        final Result res = new Result();
        if ( greatest_common_prefix.length() < 1 ) {
            res.addWarning( "No greatest common prefix" );
            res.setGreatestCommonPrefix( "" );
        }
        else {
            res.setGreatestCommonPrefix( greatest_common_prefix );
        }
        if ( qnode_pp.isRoot() ) {
            res.addWarning( "Least Encompassing Clade is entire tree" );
        }
        res.setLeastEncompassingCladeSize( lec_ext_nodes );
        res.setTreeSize( p_ext_nodes );
        final String greatest_common_prefix_a = analyzeSiblings( qnode_p, qnode_pp );
        res.setGreatestCommonPrefixUp( greatest_common_prefix_a );
        final String greatest_common_prefix_b = analyzeSiblings( qnode, qnode_p );
        res.setGreatestCommonPrefixDown( greatest_common_prefix_b );
        return res;
    }

    private final static String analyzeSiblings( final PhylogenyNode child, final PhylogenyNode parent ) {
        final int child_index = child.getChildNodeIndex();
        final List<String> ext_nodes_names = new ArrayList<>();
        final List<PhylogenyNode> descs = parent.getDescendants();
        for( int i = 0; i < descs.size(); ++i ) {
            if ( i != child_index ) {
                final PhylogenyNode d = descs.get( i );
                for( final PhylogenyNode n : d.getAllExternalDescendants() ) {
                    final String name = n.getName();
                    if ( ForesterUtil.isEmptyTrimmed( name ) ) {
                        throw new IllegalArgumentException( "external node(s) with empty names found" );
                    }
                    ext_nodes_names.add( name.trim() );
                }
            }
        }
        final String greatest_common_prefix = ForesterUtil.greatestCommonPrefix( ext_nodes_names );
        return greatest_common_prefix;
    }
}
