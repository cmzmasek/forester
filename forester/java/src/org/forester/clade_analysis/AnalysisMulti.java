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

public final class AnalysisMulti {

    private final static String UNKNOWN = "?";
    public final static double DEFAULT_CUTOFF_FOR_SPECIFICS = 0.5;
    public final static String DEFAULT_SEPARATOR = ".";
    public final static Pattern DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE = Pattern.compile( ".+#\\d+_M=(.+)" );
    

    public static ResultMulti execute( final Phylogeny p ) {
        return execute( p, DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE, DEFAULT_SEPARATOR, DEFAULT_CUTOFF_FOR_SPECIFICS );
    }
    
    public static ResultMulti execute( final Phylogeny p, final String separator ) {
        return execute( p, DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE, separator , DEFAULT_CUTOFF_FOR_SPECIFICS );
    }
    
    public static ResultMulti execute( final Phylogeny p, final String separator,  final double cutoff_for_specifics ) {
        return execute( p, DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE, separator , cutoff_for_specifics );
    }
    
    public static ResultMulti execute( final Phylogeny p, final double cutoff_for_specifics ) {
        return execute( p, DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE, DEFAULT_SEPARATOR , cutoff_for_specifics );
    }

    public static ResultMulti execute( final Phylogeny p,
                                   final Pattern query,
                                   final String separator,
                                   final double cutoff_for_specifics ) {
        final List<PhylogenyNode> qnodes = p.getNodes( query );
        final ResultMulti res = new ResultMulti();
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
            final String greatest_common_prefix = ForesterUtil.greatestCommonPrefix( qnode_ext_nodes_names, separator );
            //  System.out.println( greatest_common_prefix );
            Matcher matcher = query.matcher( qnode.getName() );
            String conf_str = null;
            if ( matcher.find() ) {
                conf_str = matcher.group( 1 );
            }
            else {
                throw new IllegalStateException( "pattern did not match -- this should have never happened!" );
            }
            final double conf = Double.parseDouble( conf_str );
            if ( !ForesterUtil.isEmpty( greatest_common_prefix ) ) {
                res.addGreatestCommonPrefix( greatest_common_prefix, conf );
            }
            else {
                res.addGreatestCommonPrefix( UNKNOWN, conf );
            }
            //final String greatest_common_prefix_up[] = analyzeSiblings( qnode_p, qnode_pp, separator, query, res );
            final String greatest_common_prefix_up = analyzeSiblings( qnode_p, qnode_pp, separator, query );
            System.out.println( "greatest_common_prefix_up=" + greatest_common_prefix_up + " " + conf );
            if ( !ForesterUtil.isEmpty( greatest_common_prefix_up ) ) {
                res.addGreatestCommonPrefixUp( greatest_common_prefix_up, conf );
            }
            else {
                res.addGreatestCommonPrefixUp( UNKNOWN, conf );
            }
            final String greatest_common_prefix_down = analyzeSiblings( qnode, qnode_p, separator, query );
            System.out.println( "greatest_common_prefix_down=" + greatest_common_prefix_down + " " + conf );
            if ( !ForesterUtil.isEmpty( greatest_common_prefix_down ) ) {
                res.addGreatestCommonPrefixDown( greatest_common_prefix_down, conf );
            }
            else {
                res.addGreatestCommonPrefixDown( UNKNOWN, conf );
            }
        }
        res.analyze( cutoff_for_specifics );
        return res;
    }

    private final static String analyzeSiblings( final PhylogenyNode child,
                                                 final PhylogenyNode parent,
                                                 final String separator,
                                                 final Pattern query ) {
        final int child_index = child.getChildNodeIndex();
        final List<String> ext_nodes_names = new ArrayList<>();
        final List<PhylogenyNode> descs = parent.getDescendants();
        // String conf = null;
        for( int i = 0; i < descs.size(); ++i ) {
            if ( i != child_index ) {
                final PhylogenyNode d = descs.get( i );
                for( final PhylogenyNode n : d.getAllExternalDescendants() ) {
                    final String name = n.getName();
                    if ( ForesterUtil.isEmptyTrimmed( name ) ) {
                        throw new IllegalArgumentException( "external node(s) with empty names found" );
                    }
                    final Matcher m = query.matcher( name );
                    if ( !m.find() ) {
                        ext_nodes_names.add( name );
                    }
                }
                // if ( descs.size() == 2 ) {
                //     conf = obtainConfidence( d );
                // }
            }
        }
        final String greatest_common_prefix = ForesterUtil.greatestCommonPrefix( ext_nodes_names, separator );
        return greatest_common_prefix;
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
