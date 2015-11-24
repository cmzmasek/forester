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

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

public final class NeighborJoining {

    private final static DecimalFormat     DF = new DecimalFormat( "0.00000" );
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

    private NeighborJoining() {
        _verbose = false;
        _df = null;
    }

    private NeighborJoining( final boolean verbose, final int maximum_fraction_digits_for_distances ) {
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
            // Calculates the minimal distance.
            // If more than one minimal distances, always the first found is used
            updateM();
            final int otu1 = _min_i;
            final int otu2 = _min_j;
            //System.out.println( _min_i + " " + _min_j );
            // It is a condition that otu1 < otu2.
            final PhylogenyNode node = new PhylogenyNode();
            final double d = _d_values[ _mappings[ otu1 ] ][ _mappings[ otu2 ] ];
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
        }
        final double d = _d_values[ _mappings[ 0 ] ][ _mappings[ 1 ] ] / 2;
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
        final int m_otu1 = _mappings[ otu1 ];
        final int m_otu2 = _mappings[ otu2 ];
        for( int i = 0; i < _n; ++i ) {
            if ( ( i == otu1 ) || ( i == otu2 ) ) {
                continue;
            }
            final int m_i = _mappings[ i ];
            if ( otu1 < i ) {
                if ( otu2 > i ) {
                    _d_values[ m_otu1 ][ m_i ] = ( ( _d_values[ m_otu1 ][ m_i ] + _d_values[ m_i ][ m_otu2 ] ) - d ) / 2;
                    //System.out.print( DF.format( _d_values[ m_otu1 ][ m_i ] ) );
                }
                else {
                    _d_values[ m_otu1 ][ m_i ] = ( ( _d_values[ m_otu1 ][ m_i ] + _d_values[ m_otu2 ][ m_i ] ) - d ) / 2;
                    //System.out.print( DF.format( _d_values[ m_otu1 ][ m_i ] ) );
                }
            }
            else {
                if ( otu2 > i ) {
                    _d_values[ m_i ][ m_otu1 ] = ( ( _d_values[ m_i ][ m_otu1 ] + _d_values[ m_i ][ m_otu2 ] ) - d ) / 2;
                    //System.out.print( DF.format( _d_values[ m_i ][ m_otu1 ] ) );
                }
                else {
                    _d_values[ m_i ][ m_otu1 ] = ( ( _d_values[ m_i ][ m_otu1 ] + _d_values[ m_otu2 ][ m_i ] ) - d ) / 2;
                    // System.out.print( DF.format( _d_values[ m_otu1 ][ m_i ] ) );
                }
            }
            //System.out.print( " " );
        }
    }

    private final void calculateNetDivergences() {
        double d;
        for( int i = 0; i < _n; ++i ) {
            d = 0;
            final int m_i = _mappings[ i ];
            for( int n = 0; n < _n; ++n ) {
                if ( i != n ) {
                    if ( i > n ) {
                        d += _d_values[ _mappings[ n ] ][ m_i ];
                    }
                    else {
                        d += _d_values[ m_i ][ _mappings[ n ] ];
                    }
                }
            }
            _r[ i ] = d;
        }
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
        initExternalNodes();
    }

    private final void updateM() {
        calculateNetDivergences();
        Double min = Double.MAX_VALUE;
        _min_i = -1;
        _min_j = -1;
        final int n_minus_2 = _n - 2;
        for( int j = 1; j < _n; ++j ) {
            final double r_j = _r[ j ];
            final int m_j = _mappings[ j ];
            for( int i = 0; i < j; ++i ) {
                final double m = _d_values[ _mappings[ i ] ][ m_j ] - ( ( _r[ i ] + r_j ) / n_minus_2 );
                if ( m < min ) {
                    min = m;
                    _min_i = i;
                    _min_j = j;
                }
            }
        }
        //        for( int j = 1; j < _n; ++j ) {
        //            final double r_j = _r[ j ];
        //            final int m_j = _mappings[ j ];
        //            for( int i = 0; i < j; ++i ) {
        //                System.out.print( i );
        //                System.out.print( "->" );
        //                System.out.print( DF.format( _r[ i ] ) );
        //                System.out.print( "  " );
        //            }
        //            System.out.println();
        //        }
    }

    // otu2 will, in effect, be "deleted" from the matrix.
    private final void updateMappings( final int otu2 ) {
        for( int i = otu2; i < ( _mappings.length - 1 ); ++i ) {
            _mappings[ i ] = _mappings[ i + 1 ];
        }
    }

    public final static NeighborJoining createInstance() {
        return new NeighborJoining();
    }

    public final static NeighborJoining createInstance( final boolean verbose,
                                                        final int maximum_fraction_digits_for_distances ) {
        return new NeighborJoining( verbose, maximum_fraction_digits_for_distances );
    }
}
