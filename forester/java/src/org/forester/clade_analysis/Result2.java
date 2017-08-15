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

package org.forester.clade_analysis;

import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.util.ForesterUtil;

public final class Result2 {

    private List<Prefix>        _greatest_common_prefix                = new ArrayList<Prefix>();
    private String             _greatest_common_prefix_up             = "";
    private String             _greatest_common_prefix_down           = "";
    private final List<String> _warnings                              = new ArrayList<>();
    private int                _lec_ext_nodes                         = 0;
    private int                _p_ext_nodes                           = 0;
    private String             _greatest_common_clade_subtree_confidence      = "";
    private String             _greatest_common_clade_subtree_confidence_up   = "";
    private String             _greatest_common_clade_subtree_confidence_down = "";

    void addWarning( final String warning ) {
        _warnings.add( warning );
    }

    void addGreatestCommonPrefix( final String prefix, final double confidence ) {
        _greatest_common_prefix.add( new Prefix(prefix, confidence) );
    }

    void setGreatestCommonPrefixUp( final String greatest_common_prefix_up ) {
        _greatest_common_prefix_up = greatest_common_prefix_up;
    }

    void setGreatestCommonPrefixDown( final String greatest_common_prefix_down ) {
        _greatest_common_prefix_down = greatest_common_prefix_down;
    }

    void setGreatestCommonCladeSubtreeConfidence( final String greatest_common_clade_confidence ) {
        _greatest_common_clade_subtree_confidence = greatest_common_clade_confidence;
    }

    void setGreatestCommonCladeUpSubtreeConfidence( final String greatest_common_clade_confidence_up ) {
        _greatest_common_clade_subtree_confidence_up = greatest_common_clade_confidence_up;
    }

    void setGreatestCommonCladeDownSubtreeConfidence( final String greatest_common_clade_confidence_down ) {
        _greatest_common_clade_subtree_confidence_down = greatest_common_clade_confidence_down;
    }

  //  public String getGreatestCommonPrefix() {
  //      return _greatest_common_prefix;
  //  }

    public String getGreatestCommonPrefixUp() {
        return _greatest_common_prefix_up;
    }

    public String getGreatestCommonPrefixDown() {
        return _greatest_common_prefix_down;
    }

    public String getGreatestCommonCladeSubtreeConfidence() {
        return _greatest_common_clade_subtree_confidence;
    }

    public String getGreatestCommonCladeUpSubtreeConfidence() {
        return _greatest_common_clade_subtree_confidence_up;
    }

    public String getGreatestCommonCladeDownSubtreeConfidence() {
        return _greatest_common_clade_subtree_confidence_down;
    }

    public List<String> getWarnings() {
        return _warnings;
    }

    void setLeastEncompassingCladeSize( final int lec_ext_nodes ) {
        _lec_ext_nodes = lec_ext_nodes;
    }

    void setTreeSize( final int p_ext_nodes ) {
        _p_ext_nodes = p_ext_nodes;
    }

    public int getLeastEncompassingCladeSize() {
        return _lec_ext_nodes;
    }

    public int getTreeSize() {
        return _p_ext_nodes;
    }
    
    public void analyzeGreatestCommonPrefixes(final String separator ) {
        final SortedMap<String,Double> map = new TreeMap<String,Double>();
        for( final Prefix prefix : _greatest_common_prefix ) {
            List<String> prefixes = ForesterUtil.spliIntoPrefixes( prefix.getPrefix(), separator );
            for( final String p : prefixes ) {
                map.put( p, 0.0 );
            }
        }
       // System.out.println( map );
        for (final String key : map.keySet()) {
            //System.out.println(key);
            for( final Prefix prefix : _greatest_common_prefix ) {
                if ( prefix.getPrefix().startsWith( key ) ) {
                    map.put( key, map.get( key ) + prefix.getConfidence()  );
                }
            }
        }
        System.out.println( map );
    }
    
    
}
