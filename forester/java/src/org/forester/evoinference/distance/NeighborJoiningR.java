// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2014 Christian M. Zmasek
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

package org.forester.evoinference.distance;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.SortedSet;

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

public final class NeighborJoiningR {

    private final static DecimalFormat     DF = new DecimalFormat( "0.00" );
    private BasicSymmetricalDistanceMatrix _d;
    private double[][]                     _d_values;
    private final DecimalFormat            _df;
    private PhylogenyNode[]                _external_nodes;
    private int[]                          _mappings;
    private int                            _n;
    private double[]                       _r;
    private final boolean                  _verbose;
    private int                            _min_i;
    private int                            _min_j;
    private S                              _s;
    private double                         _d_min;                          //TODO remove me

    private NeighborJoiningR() {
        _verbose = false;
        _df = null;
    }

    private NeighborJoiningR( final boolean verbose, final int maximum_fraction_digits_for_distances ) {
        if ( ( maximum_fraction_digits_for_distances < 1 ) || ( maximum_fraction_digits_for_distances > 9 ) ) {
            throw new IllegalArgumentException( "maximum fraction digits for distances is out of range: "
                    + maximum_fraction_digits_for_distances );
        }
        _verbose = verbose;
        _df = new DecimalFormat();
        _df.setMaximumFractionDigits( maximum_fraction_digits_for_distances );
        _df.setRoundingMode( RoundingMode.HALF_UP );
    }

    public final Phylogeny execute( final BasicSymmetricalDistanceMatrix distance ) {
        reset( distance );
        final Phylogeny phylogeny = new Phylogeny();
        while ( _n > 2 ) {
            System.out.println( "N=" + _n );
            System.out.println();
            // Calculates the minimal distance.
            // If more than one minimal distances, always the first found is used
            final double m = updateM();
            final int otu1 = _min_i;
            final int otu2 = _min_j;
            System.out.println( _min_i + " " + _min_j + " => " + DF.format( m ) + " (" + DF.format( _d_min ) + ")" );
            // It is a condition that otu1 < otu2.
            final PhylogenyNode node = new PhylogenyNode();
            final double d = getDvalue( otu1, otu2 );
            final double d1 = ( d / 2 ) + ( ( _r[ otu1 ] - _r[ otu2 ] ) / ( 2 * ( _n - 2 ) ) );
            final double d2 = d - d1;
            if ( _df == null ) {
                getExternalPhylogenyNode( otu1 ).setDistanceToParent( d1 );
                getExternalPhylogenyNode( otu2 ).setDistanceToParent( d2 );
            }
            else {
                // yes, yes, slow but only grows with n (and not n^2 or worse)...
                getExternalPhylogenyNode( otu1 ).setDistanceToParent( Double.parseDouble( _df.format( d1 ) ) );
                getExternalPhylogenyNode( otu2 ).setDistanceToParent( Double.parseDouble( _df.format( d2 ) ) );
            }
            node.addAsChild( getExternalPhylogenyNode( otu1 ) );
            node.addAsChild( getExternalPhylogenyNode( otu2 ) );
            if ( _verbose ) {
                printProgress( otu1, otu2 );
            }
            calculateDistancesFromNewNode( otu1, otu2, d );
            _external_nodes[ _mappings[ otu1 ] ] = node;
            updateMappings( otu2 );
            --_n;
            System.out.println( "-------------------------------------------------------------" );
            System.out.println( "" );
        }
        final double d = getDvalue( 0, 1 ) / 2;
        if ( _df == null ) {
            getExternalPhylogenyNode( 0 ).setDistanceToParent( d );
            getExternalPhylogenyNode( 1 ).setDistanceToParent( d );
        }
        else {
            final double dd = Double.parseDouble( _df.format( d ) );
            getExternalPhylogenyNode( 0 ).setDistanceToParent( dd );
            getExternalPhylogenyNode( 1 ).setDistanceToParent( dd );
        }
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( getExternalPhylogenyNode( 0 ) );
        root.addAsChild( getExternalPhylogenyNode( 1 ) );
        if ( _verbose ) {
            printProgress( 0, 1 );
        }
        phylogeny.setRoot( root );
        phylogeny.setRooted( false );
        return phylogeny;
    }

    public final List<Phylogeny> execute( final List<BasicSymmetricalDistanceMatrix> distances_list ) {
        final List<Phylogeny> pl = new ArrayList<Phylogeny>();
        for( final BasicSymmetricalDistanceMatrix distances : distances_list ) {
            pl.add( execute( distances ) );
        }
        return pl;
    }

    private final void calculateDistancesFromNewNode( final int otu1, final int otu2, final double d ) {
        for( int i = 0; i < _n; ++i ) {
            if ( ( i == otu1 ) || ( i == otu2 ) ) {
                continue;
            }
            updateDvalue( otu1, otu2, i, d );
        }
    }

    private final void updateDvalue( final int otu1, final int otu2, final int i, final double d ) {
        setDvalue( otu1, i, ( getDvalue( otu1, i ) + getDvalue( i, otu2 ) - d ) / 2 );
    }

    private void setDvalue( final int i, final int j, final double d ) {
        if ( i < j ) {
            _d_values[ _mappings[ i ] ][ _mappings[ j ] ] = d;
        }
        _d_values[ _mappings[ j ] ][ _mappings[ i ] ] = d;
    }

    private double getDvalue( final int i, final int j ) {
        if ( i < j ) {
            return _d_values[ _mappings[ i ] ][ _mappings[ j ] ];
        }
        return _d_values[ _mappings[ j ] ][ _mappings[ i ] ];
    }

    private double getDvalueUnmapped( final int i, final int j ) {
        if ( i < j ) {
            return _d_values[ i ][ j ];
        }
        return _d_values[ j ][ i ];
    }

    private final void calculateNetDivergences() {
        for( int i = 0; i < _n; ++i ) {
            _r[ i ] = calculateNetDivergence( i );
        }
    }

    private double calculateNetDivergence( final int i ) {
        double d = 0;
        for( int n = 0; n < _n; ++n ) {
            if ( i != n ) {
                d += getDvalue( n, i );
            }
        }
        return d;
    }

    private final PhylogenyNode getExternalPhylogenyNode( final int i ) {
        return _external_nodes[ _mappings[ i ] ];
    }

    private final void initExternalNodes() {
        _external_nodes = new PhylogenyNode[ _n ];
        String id;
        for( int i = 0; i < _n; ++i ) {
            _external_nodes[ i ] = new PhylogenyNode();
            id = _d.getIdentifier( i );
            if ( id != null ) {
                _external_nodes[ i ].setName( id );
            }
            else {
                _external_nodes[ i ].setName( Integer.toString( i ) );
            }
            _mappings[ i ] = i;
        }
    }

    private final void printProgress( final int otu1, final int otu2 ) {
        System.out.println( "Node " + printProgressNodeToString( getExternalPhylogenyNode( otu1 ) ) + " joins "
                + ( printProgressNodeToString( getExternalPhylogenyNode( otu2 ) ) ) );
    }

    private final String printProgressNodeToString( final PhylogenyNode n ) {
        if ( n.isExternal() ) {
            if ( ForesterUtil.isEmpty( n.getName() ) ) {
                return Long.toString( n.getId() );
            }
            return n.getName();
        }
        return n.getId()
                + " ("
                + ( ForesterUtil.isEmpty( n.getChildNode1().getName() ) ? n.getChildNode1().getId() : n.getChildNode1()
                        .getName() )
                + "+"
                + ( ForesterUtil.isEmpty( n.getChildNode2().getName() ) ? n.getChildNode2().getId() : n.getChildNode2()
                        .getName() ) + ")";
    }

    // only the values in the lower triangle are used.
    // !matrix values will be changed!
    private final void reset( final BasicSymmetricalDistanceMatrix distances ) {
        _n = distances.getSize();
        _d = distances;
        _r = new double[ _n ];
        _mappings = new int[ _n ];
        _d_values = _d.getValues();
        _s = new S();
        _s.initialize( distances );
        initExternalNodes();
        printM();
    }

    final private void printM() {
        for( int j = 1; j < _n; ++j ) {
            for( int i = 0; i < _n; ++i ) {
                System.out.print( DF.format( _d_values[ _mappings[ i ] ][ _mappings[ j ] ] ) );
                System.out.print( " " );
            }
            System.out.print( "    " );
            for( final Entry<Double, SortedSet<Integer>> entry : _s.getSentrySet( _mappings[ j ] ) ) {
                final double key = entry.getKey();
                final SortedSet<Integer> value = entry.getValue();
                System.out.print( DF.format( key ) + "=" );
                boolean first = true;
                for( final Integer v : value ) {
                    if ( !first ) {
                        System.out.print( "," );
                    }
                    first = false;
                    System.out.print( v );
                }
                System.out.print( "  " );
            }
            System.out.println();
        }
    }

    private final double updateM() {
        printM();
        calculateNetDivergences();
        Double min = Double.MAX_VALUE;
        _min_i = -1;
        _min_j = -1;
        final int n_minus_2 = _n - 2;
        for( int j = 1; j < _n; ++j ) {
            final double r_j = _r[ j ];
            final int m_j = _mappings[ j ];
            int counter = 0;
            int counter_all = 0;
            X: for( final Entry<Double, SortedSet<Integer>> entry : _s.getSentrySet( m_j ) ) {
                for( final int sorted_i : entry.getValue() ) {
                    //if ( counter_all >= j ) {
                    //    break X;
                    //}
                    if ( _mappings[ counter ] == counter_all ) {
                        System.out.print( sorted_i + " " );
                        System.out.print( "(" + DF.format( getDvalueUnmapped( sorted_i, m_j ) ) + ") " );
                        final double m = getDvalueUnmapped( sorted_i, m_j ) - ( ( _r[ sorted_i ] + r_j ) / n_minus_2 );
                        if ( m < min ) {
                            _d_min = getDvalueUnmapped( sorted_i, m_j );
                            min = m;
                            _min_i = sorted_i;
                            _min_j = j;
                        }
                        ++counter;
                    }
                    ++counter_all;
                }
            }
            System.out.println();
            /*
            for( int i = 0; i < j; ++i ) {
                final double m = getDvalue( i, j ) - ( ( _r[ i ] + r_j ) / n_minus_2 );
                if ( m < min ) {
                    min = m;
                    _d_min = getDvalue( i, j );
                    _min_i = i;
                    _min_j = j;
                }
            }*/
        }
        System.out.println();
        return min;
    }

    // otu2 will, in effect, be "deleted" from the matrix.
    private final void updateMappings( final int otu2 ) {
        for( int i = otu2; i < ( _mappings.length - 1 ); ++i ) {
            _mappings[ i ] = _mappings[ i + 1 ];
        }
    }

    public final static NeighborJoiningR createInstance() {
        return new NeighborJoiningR();
    }

    public final static NeighborJoiningR createInstance( final boolean verbose,
                                                         final int maximum_fraction_digits_for_distances ) {
        return new NeighborJoiningR( verbose, maximum_fraction_digits_for_distances );
    }
}
