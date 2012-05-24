// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
// WWW: www.phylosoft.org/forester

package org.forester.evoinference.distance;

import java.util.ArrayList;
import java.util.List;

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

public final class NeighborJoining {

    private BasicSymmetricalDistanceMatrix _d;
    private BasicSymmetricalDistanceMatrix _m;
    private double[]                       _r;
    private int                            _n;
    private PhylogenyNode[]                _external_nodes;
    private int[]                          _mappings;
    private final boolean                  _verbose;
    private final static boolean           DEBUG = false;

    private NeighborJoining( final boolean verbose ) {
        _verbose = verbose;
    }

    private final void calculateDistancesFromNewNode( final int otu1, final int otu2, final double d ) {
        for( int i = 0; i < _n; ++i ) {
            if ( ( i == otu1 ) || ( i == otu2 ) ) {
                continue;
            }
            //final double nd = ;
            //  setValueInD( nd, otu1, i );
            //   _d.setValue( _mappings[ otu1 ],
            //                _mappings[ i ],
            //                ( getValueFromD( otu1, i ) + getValueFromD( i, otu2 ) - d ) / 2 );
            _d._values[ _mappings[ otu1 ] ][ _mappings[ i ] ] = ( getValueFromD( otu1, i ) + getValueFromD( i, otu2 ) - d ) / 2;
        }
    }

    private final void calculateNetDivergences() {
        double d;
        for( int i = 0; i < _n; ++i ) {
            d = 0;
            for( int n = 0; n < _n; ++n ) {
                d += getValueFromD( i, n );
            }
            _r[ i ] = d;
        }
    }

    public final Phylogeny execute( final BasicSymmetricalDistanceMatrix distance ) {
        reset( distance );
        final Phylogeny phylogeny = new Phylogeny();
        while ( _n > 2 ) {
            updateM();
            //final int[] s = findMinimalDistance();
            // Calculates the minimal distance.
            // If more than one minimal distances, always the first found is used
            // could randomize this, so that any would be returned in a randomized fashion...
            double minimum = Double.MAX_VALUE;
            int otu1 = -1;
            int otu2 = -1;
            for( int j = 1; j < _n; ++j ) {
                for( int i = 0; i < j; ++i ) {
                    //   if ( _m.getValue( i, j ) < minimum ) {
                    if ( _m._values[ i ][ j ] < minimum ) {
                        //minimum = _m.getValue( i, j );
                        minimum = _m._values[ i ][ j ];
                        otu1 = i;
                        otu2 = j;
                    }
                }
            }
            //
            // final int otu1 = s[ 0 ];
            // final int otu2 = s[ 1 ];
            // It is a condition that otu1 < otu2.
            if ( DEBUG ) {
                if ( otu1 > otu2 ) {
                    throw new RuntimeException( "NJ code is faulty: otu1 > otu2" );
                }
            }
            final PhylogenyNode node = new PhylogenyNode();
            final double d = getValueFromD( otu1, otu2 );
            final double d1 = ( d / 2 ) + ( ( _r[ otu1 ] - _r[ otu2 ] ) / ( 2 * ( _n - 2 ) ) );
            final double d2 = d - d1;
            getExternalPhylogenyNode( otu1 ).setDistanceToParent( d1 );
            getExternalPhylogenyNode( otu2 ).setDistanceToParent( d2 );
            node.addAsChild( getExternalPhylogenyNode( otu1 ) );
            node.addAsChild( getExternalPhylogenyNode( otu2 ) );
            if ( _verbose ) {
                printProgress( otu1, otu2 );
            }
            calculateDistancesFromNewNode( otu1, otu2, d );
            //setExternalPhylogenyNode( node, otu1 );
            _external_nodes[ _mappings[ otu1 ] ] = node;
            updateMappings( otu2 );
            --_n;
        }
        final double d = getValueFromD( 0, 1 ) / 2;
        getExternalPhylogenyNode( 0 ).setDistanceToParent( d );
        getExternalPhylogenyNode( 1 ).setDistanceToParent( d );
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

    //    private int[] findMinimalDistance() {
    //        // if more than one minimal distances, always the first found is
    //        // returned
    //        // i could randomize this, so that any would be returned in a randomized
    //        // fashion...
    //        double minimum = Double.MAX_VALUE;
    //        int otu_1 = -1;
    //        int otu_2 = -1;
    //        for( int j = 1; j < _n; ++j ) {
    //            for( int i = 0; i < j; ++i ) {
    //                if ( _m.getValue( i, j ) < minimum ) {
    //                    minimum = _m.getValue( i, j );
    //                    otu_1 = i;
    //                    otu_2 = j;
    //                }
    //            }
    //        }
    //        return new int[] { otu_1, otu_2 };
    //    }
    private final PhylogenyNode getExternalPhylogenyNode( final int i ) {
        return _external_nodes[ _mappings[ i ] ];
    }

    private final double getValueFromD( final int otu1, final int otu2 ) {
        //return _d.getValue( _mappings[ otu1 ], _mappings[ otu2 ] );
        return _d._values[ _mappings[ otu1 ] ][ _mappings[ otu2 ] ];
    }

    private final void initExternalNodes() {
        _external_nodes = new PhylogenyNode[ _n ];
        for( int i = 0; i < _n; ++i ) {
            _external_nodes[ i ] = new PhylogenyNode();
            final String id = _d.getIdentifier( i );
            if ( id != null ) {
                _external_nodes[ i ].setName( id );
            }
            else {
                _external_nodes[ i ].setName( "" + i );
            }
            _mappings[ i ] = i;
        }
    }

    private final void printProgress( final int otu1, final int otu2 ) {
        final PhylogenyNode n1 = getExternalPhylogenyNode( otu1 );
        final PhylogenyNode n2 = getExternalPhylogenyNode( otu2 );
        System.out.println( "Node " + ( ForesterUtil.isEmpty( n1.getName() ) ? n1.getId() : n1.getName() ) + " joins "
                + ( ForesterUtil.isEmpty( n2.getName() ) ? n2.getId() : n2.getName() ) );
    }

    // only the values in the lower triangle are used.
    // !matrix values will be changed!
    private final void reset( final BasicSymmetricalDistanceMatrix distances ) {
        _n = distances.getSize();
        _d = distances;
        _m = new BasicSymmetricalDistanceMatrix( _n );
        _r = new double[ _n ];
        _mappings = new int[ _n ];
        initExternalNodes();
    }

    //  private final void setExternalPhylogenyNode( final PhylogenyNode node, final int i ) {
    //      _external_nodes[ _mappings[ i ] ] = node;
    //  }
    //  private final void setValueInD( final double d, final int otu1, final int otu2 ) {
    //      _d.setValue( _mappings[ otu1 ], _mappings[ otu2 ], d );
    //  }
    private final void updateM() {
        calculateNetDivergences();
        for( int j = 1; j < _n; ++j ) {
            for( int i = 0; i < j; ++i ) {
                //_m.setValue( i, j, calculateM( i, j ) );
                //_m.setValue( i, j, getValueFromD( i, j ) - ( _r[ i ] + _r[ j ] ) / ( _n - 2 ) );
                _m._values[ i ][ j ] = getValueFromD( i, j ) - ( _r[ i ] + _r[ j ] ) / ( _n - 2 );
            }
        }
    }

    //private double calculateM( final int i, final int j ) {
    //    return getValueFromD( i, j ) - ( _r[ i ] + _r[ j ] ) / ( _n - 2 );
    //}
    // otu2 will, in effect, be "deleted" from the matrix.
    private final void updateMappings( final int otu2 ) {
        for( int i = otu2; i < _mappings.length - 1; ++i ) {
            _mappings[ i ] = _mappings[ i + 1 ];
        }
    }

    public final static NeighborJoining createInstance() {
        return new NeighborJoining( false );
    }

    public final static NeighborJoining createInstance( final boolean verbose ) {
        return new NeighborJoining( verbose );
    }
}
