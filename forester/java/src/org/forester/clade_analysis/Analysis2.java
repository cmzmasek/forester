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

package org.forester.clade_analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.util.ForesterUtil;

public final class Analysis2 {

    public static Result2 execute( final Phylogeny p, final Pattern query, final String separator ) {
        final List<PhylogenyNode> qnodes = p.getNodes( query );
        final Result2 res = new Result2();
        for( int i = 0; i < qnodes.size(); ++i ) {
            final PhylogenyNode qnode = qnodes.get( i );
            System.out.println( ">>" + qnode.getName() );
            if ( qnode.isRoot() ) {
                throw new IllegalArgumentException( "Query " + query + " is root." );
            }
            if ( qnode.getParent().isRoot() ) {
                throw new IllegalArgumentException( "Parent of query " + query + " is root." );
            }
            PhylogenyNode qnode_p = qnode.getParent();
            PhylogenyNode qnode_pp = qnode.getParent().getParent();
            //This is to deal with internal nodes with 1 descendant.
            while ( qnode_p.getNumberOfDescendants() == 1 ) {
                qnode_p = qnode_p.getParent();
            }
            while ( qnode_pp.getNumberOfDescendants() == 1 ) {
                qnode_pp = qnode_pp.getParent();
            }
            // final List<PhylogenyNode> qnode_ext_nodes = new ArrayList<PhylogenyNode>();
            final List<String> qnode_ext_nodes_names = new ArrayList<>();
            for( final PhylogenyNode qnode_ext_node : qnode_pp.getAllExternalDescendants() ) {
                final String name = qnode_ext_node.getName();
                if ( ForesterUtil.isEmptyTrimmed( name ) ) {
                    throw new IllegalArgumentException( "external node(s) with empty names found" );
                }
                final Matcher m = query.matcher( name );
                if ( !m.find() ) {
                    qnode_ext_nodes_names.add( name );
                }
            }
            final int lec_ext_nodes = qnode_ext_nodes_names.size();
            final int p_ext_nodes = p.getNumberOfExternalNodes() - 1;
            final String greatest_common_prefix = ForesterUtil.greatestCommonPrefix( qnode_ext_nodes_names, separator );
            System.out.println( greatest_common_prefix );
            Matcher matcher = query.matcher( qnode.getName() );
            String conf_str = null;
            if ( matcher.find() ) {
                conf_str = matcher.group( 1 );
            }
            else {
                throw new IllegalStateException( "pattern did not match -- this should have never happened!" );
            }
            res.setLeastEncompassingCladeSize( lec_ext_nodes );
            res.setTreeSize( p_ext_nodes );
            final double conf = Double.parseDouble( conf_str );
            if ( !ForesterUtil.isEmpty( greatest_common_prefix ) ) {
                res.addGreatestCommonPrefix( greatest_common_prefix, conf );
            }
            else {
                res.addGreatestCommonPrefix( "?", conf );
            }
        }
        /* for( final PhylogenyNode qnode_ext_node : qnode_ext_nodes ) {
            String name = qnode_ext_node.getName();
            if ( ForesterUtil.isEmptyTrimmed( name ) ) {
                throw new IllegalArgumentException( "external node(s) with empty names found" );
            }
            name = name.trim();
            if ( !name.equals( query ) ) {
                qnode_ext_nodes_names.add( name );
            }
        }*/
        //   if ( greatest_common_prefix.length() < 1 ) {
        //       res.addWarning( "No greatest common prefix" );
        //res.setGreatestCommonPrefix( "" );
        //  }
        // else {
        //    //  res.setGreatestCommonPrefix( greatest_common_prefix );
        // res.addGreatestCommonPrefix( prefix, confidence, separator ); //TODO
        //   }
        // if ( qnode_pp.isRoot() ) {
        //     res.addWarning( "Least Encompassing Clade is entire tree" );
        // }
        /*    final String conf = obtainConfidence( qnode_pp );
        if ( conf != null ) {
            res.setGreatestCommonCladeSubtreeConfidence(conf);
        }*/
        /*  final String greatest_common_prefix_up[] = analyzeSiblings( qnode_p, qnode_pp, separator );
        res.setGreatestCommonPrefixUp( greatest_common_prefix_up[ 0 ] );
        if ( greatest_common_prefix_up[ 1 ] != null ) {
            res.setGreatestCommonCladeUpSubtreeConfidence( greatest_common_prefix_up[ 1 ] );
        }
        final String greatest_common_prefix_down[] = analyzeSiblings( qnode, qnode_p, separator );
        res.setGreatestCommonPrefixDown( greatest_common_prefix_down[ 0 ] );
        if ( greatest_common_prefix_down[ 1 ] != null ) {
            res.setGreatestCommonCladeDownSubtreeConfidence( greatest_common_prefix_down[ 1 ] );
        }*/
        return res;
    }

    private final static String[] analyzeSiblings( final PhylogenyNode child,
                                                   final PhylogenyNode parent,
                                                   final String separator ) {
        final int child_index = child.getChildNodeIndex();
        final List<String> ext_nodes_names = new ArrayList<>();
        final List<PhylogenyNode> descs = parent.getDescendants();
        String conf = null;
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
                if ( descs.size() == 2 ) {
                    conf = obtainConfidence( d );
                }
            }
        }
        final String greatest_common_prefix = ForesterUtil.greatestCommonPrefix( ext_nodes_names, separator );
        return new String[] { greatest_common_prefix, conf };
    }

    private final static String obtainConfidence( final PhylogenyNode n ) {
        if ( n.getBranchData().getConfidences() != null && n.getBranchData().getConfidences().size() > 0 ) {
            final List<Confidence> confidences = n.getBranchData().getConfidences();
            boolean not_first = false;
            Collections.sort( confidences );
            final StringBuilder sb = new StringBuilder();
            for( final Confidence confidence : confidences ) {
                final double value = confidence.getValue();
                if ( value != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                    if ( not_first ) {
                        sb.append( " / " );
                    }
                    else {
                        not_first = true;
                    }
                    sb.append( ( ForesterUtil.isEmpty( confidence.getType() ) ? "confidence: "
                            : confidence.getType() + ": " ) + value );
                }
            }
            return sb.toString();
        }
        return null;
    }
}
