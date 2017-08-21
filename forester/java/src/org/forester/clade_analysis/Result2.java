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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.util.ForesterUtil;

public final class Result2 {

    private final String       _separator;
    private final List<Prefix> _greatest_common_prefixes                      = new ArrayList<>();
    private String             _greatest_common_prefix_up                     = "";
    private String             _greatest_common_prefix_down                   = "";
    private final List<String> _warnings                                      = new ArrayList<>();
    private int                _lec_ext_nodes                                 = 0;
    private int                _p_ext_nodes                                   = 0;
    private String             _greatest_common_clade_subtree_confidence      = "";
    private String             _greatest_common_clade_subtree_confidence_up   = "";
    private String             _greatest_common_clade_subtree_confidence_down = "";
    private List<Prefix>       _all                                           = null;
    private List<Prefix>       _collapsed                                     = null;
    private List<Prefix>       _cleaned_spec                                  = null;
    private boolean            _has_specifics;

    public Result2( final String separator ) {
        _separator = separator;
    }

    public Result2() {
        _separator = ".";//TODO make const somewhere
    }

    public List<Prefix> getAllMultiHitPrefixes() {
        return _all;
    }
    
    public List<Prefix> getCollapsedMultiHitPrefixes() {
        return _collapsed;
    }
    
    public List<Prefix> getSpecificMultiHitPrefixes() {
        return _cleaned_spec;
    }
    
    public boolean isHasSpecificMultiHitsPrefixes() {
        return _has_specifics;
    }
    
    
    void addWarning( final String warning ) {
        _warnings.add( warning );
    }

    void addGreatestCommonPrefix( final String prefix, final double confidence ) {
        _greatest_common_prefixes.add( new Prefix( prefix, confidence, _separator ) );
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

    public void analyzeGreatestCommonPrefixes( final double cutoff ) {
        analyzeGreatestCommonPrefixes( _greatest_common_prefixes, _separator, cutoff );
    }

    public void analyzeGreatestCommonPrefixes() {
        analyzeGreatestCommonPrefixes( _greatest_common_prefixes, _separator, -1 );
    }

    private final void analyzeGreatestCommonPrefixes( final List<Prefix> greatest_common_prefixes,
                                                      final String separator,
                                                      final double cutoff ) {
        final List<Prefix> l = obtainAllPrefixes( greatest_common_prefixes, separator );
        sortPrefixesAccordingToConfidence( l );
        _all = removeLessSpecificPrefixes( l );
        _collapsed = collapse( _all );
        _has_specifics = false;
        if ( cutoff >= 0 ) {
            _cleaned_spec = obtainSpecifics( cutoff, _all, _collapsed );
            if ( _cleaned_spec != null && _cleaned_spec.size() > 0 ) {
                _has_specifics = true;
            }
        }
        else {
            _cleaned_spec = null;
        }
    }

    private final static List<Prefix> obtainSpecifics( final double cutoff,
                                                       final List<Prefix> cleaned,
                                                       final List<Prefix> collapsed ) {
        final List<Prefix> cleaned_spec = new ArrayList<>();
        final Set<String> collapsed_set = new HashSet<>();
        for( final Prefix prefix : collapsed ) {
            collapsed_set.add( prefix.getPrefix() );
        }
        final List<Prefix> spec = new ArrayList<>();
        for( final Prefix prefix : cleaned ) {
            if ( ( prefix.getConfidence() >= cutoff ) && !collapsed_set.contains( prefix.getPrefix() ) ) {
                spec.add( prefix );
            }
        }
        for( final Prefix o : spec ) {
            boolean ok = true;
            for( final Prefix i : spec ) {
                if ( ( !o.getPrefix().equals( i.getPrefix() ) ) && ( i.getPrefix().startsWith( o.getPrefix() ) ) ) {
                    ok = false;
                    break;
                }
            }
            if ( ok ) {
                cleaned_spec.add( o );
            }
        }
        return cleaned_spec;
    }

    private final static List<Prefix> collapse( final List<Prefix> cleaned ) {
        final List<Prefix> collapsed = new ArrayList<>();
        final Set<String> firsts = new HashSet<>();
        double confidence_sum = 0;
        for( final Prefix prefix : cleaned ) {
            final String f = prefix.getPrefixFirstElement();
            if ( !firsts.contains( f ) ) {
                firsts.add( f );
                collapsed.add( prefix );
                confidence_sum += prefix.getConfidence();
            }
        }
        if ( !ForesterUtil.isEqual( confidence_sum, 1.0, 1E-5 ) ) {
            throw new IllegalArgumentException( "Confidences add up to " + confidence_sum + " instead of 1.0" );
        }
        return collapsed;
    }

    /*
     * This replaces (by way of example)
     * A.1.1 0.9
     * A.1   0.9
     * with
     * A.1.1 0.9
     *
     * I.e. it removes less specific prefixes.
     *
     */
    private final static List<Prefix> removeLessSpecificPrefixes( final List<Prefix> l ) {
        final List<Prefix> cleaned = new ArrayList<>();
        for( final Prefix o : l ) {
            boolean ok = true;
            for( final Prefix i : l ) {
                if ( ( !o.getPrefix().equals( i.getPrefix() ) ) && ( i.getPrefix().startsWith( o.getPrefix() ) )
                        && ForesterUtil.isEqual( i.getConfidence(),
                                                 o.getConfidence() ) ) {
                    ok = false;
                    break;
                }
            }
            if ( ok ) {
                cleaned.add( o );
            }
        }
        return cleaned;
    }

    private static void sortPrefixesAccordingToConfidence( final List<Prefix> l ) {
        Collections.sort( l, new Comparator<Prefix>() {

            @Override
            public int compare( final Prefix x, final Prefix y ) {
                return compare( x.getConfidence(), y.getConfidence() );
            }

            private int compare( final double a, final double b ) {
                return a > b ? -1 : a > b ? 1 : 0;
            }
        } );
    }

    private final static List<Prefix> obtainAllPrefixes( final List<Prefix> greatest_common_prefixes,
                                                         final String separator ) {
        final SortedMap<String, Double> map = new TreeMap<>();
        for( final Prefix prefix : greatest_common_prefixes ) {
            final List<String> prefixes = ForesterUtil.spliIntoPrefixes( prefix.getPrefix(), separator );
            for( final String p : prefixes ) {
                map.put( p, 0.0 );
            }
        }
        for( final String key : map.keySet() ) {
            for( final Prefix prefix : greatest_common_prefixes ) {
                if ( prefix.getPrefix().startsWith( key ) ) {
                    map.put( key, map.get( key ) + prefix.getConfidence() );
                }
            }
        }
        final List<Prefix> l = new ArrayList<>();
        for( final Entry<String, Double> entry : map.entrySet() ) {
            l.add( new Prefix( entry.getKey(), entry.getValue(), separator ) );
        }
        return l;
    }

    public final String toString() {
        final StringBuilder sb = new StringBuilder();
        //TODO add all other stuff
        sb.append( "Cleaned:" );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        for( final Prefix prefix : _all ) {
            sb.append( prefix );
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        sb.append( ForesterUtil.LINE_SEPARATOR );
        sb.append( "Collapsed:" );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        for( final Prefix prefix : _collapsed ) {
            sb.append( prefix );
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        if ( _has_specifics ) {
            sb.append( ForesterUtil.LINE_SEPARATOR );
            sb.append( "Specifics:" );
            sb.append( ForesterUtil.LINE_SEPARATOR );
            for( final Prefix prefix : _cleaned_spec ) {
                sb.append( prefix );
                sb.append( ForesterUtil.LINE_SEPARATOR );
            }
            sb.append( ForesterUtil.LINE_SEPARATOR );
            sb.append( "Collapsed with specifics:" );
            sb.append( ForesterUtil.LINE_SEPARATOR );
            for( final Prefix prefix : _collapsed ) {
                sb.append( prefix );
                sb.append( ForesterUtil.LINE_SEPARATOR );
                for( final Prefix spec : _cleaned_spec ) {
                    if ( spec.getPrefix().startsWith( prefix.getPrefix() ) ) {
                        sb.append( "    " + spec );
                        sb.append( ForesterUtil.LINE_SEPARATOR );
                    }
                }
            }
        }
        return sb.toString();
    }
}
