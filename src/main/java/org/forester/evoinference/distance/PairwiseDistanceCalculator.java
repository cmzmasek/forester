// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.msa.Msa;

public final class PairwiseDistanceCalculator {

    public static final double DEFAULT_VALUE_FOR_TOO_LARGE_DISTANCE_FOR_KIMURA_FORMULA = 10; // Felsenstein uses -1
    private final Msa          _msa;
    private final double       _value_for_too_large_distance_for_kimura_formula;

    private PairwiseDistanceCalculator( final Msa msa, final double value_for_too_large_distance_for_kimura_formula ) {
        _msa = msa;
        _value_for_too_large_distance_for_kimura_formula = value_for_too_large_distance_for_kimura_formula;
    }

    private double calcFractionalDissimilarity( final int row_1, final int row_2 ) {
        final int length = _msa.getLength();
        int nd = 0;
        for( int col = 0; col < length; ++col ) {
            if ( _msa.getResidueAt( row_1, col ) != _msa.getResidueAt( row_2, col ) ) {
                ++nd;
            }
        }
        return ( double ) nd / length;
    }

    /**
     * "Kimura Distance"
     * Kimura, 1983
     *
     * @param row_1
     * @param row_2
     * @return
     */
    private double calcKimuraDistance( final int row_1, final int row_2 ) {
        final double p = calcFractionalDissimilarity( row_1, row_2 );
        final double dp = 1 - p - ( 0.2 * p * p );
        if ( dp <= 0.0 ) {
            return _value_for_too_large_distance_for_kimura_formula;
        }
        if ( dp == 1 ) {
            return 0; // Too avoid -0.
        }
        return -Math.log( dp );
    }

    private double calcPoissonDistance( final int row_1, final int row_2 ) {
        final double p = calcFractionalDissimilarity( row_1, row_2 );
        final double dp = 1 - p;
        if ( dp <= 0.0 ) {
            return _value_for_too_large_distance_for_kimura_formula;
        }
        if ( dp == 1 ) {
            return 0; // Too avoid -0.
        }
        return -Math.log( dp );
    }

    private BasicSymmetricalDistanceMatrix calcKimuraDistances() {
        final int s = _msa.getNumberOfSequences();
        final BasicSymmetricalDistanceMatrix d = new BasicSymmetricalDistanceMatrix( s );
        copyIdentifiers( s, d );
        calcKimuraDistances( s, d );
        return d;
    }

    private BasicSymmetricalDistanceMatrix calcPoissonDistances() {
        final int s = _msa.getNumberOfSequences();
        final BasicSymmetricalDistanceMatrix d = new BasicSymmetricalDistanceMatrix( s );
        copyIdentifiers( s, d );
        calcPoissonDistances( s, d );
        return d;
    }

    private BasicSymmetricalDistanceMatrix calcFractionalDissimilarities() {
        final int s = _msa.getNumberOfSequences();
        final BasicSymmetricalDistanceMatrix d = new BasicSymmetricalDistanceMatrix( s );
        copyIdentifiers( s, d );
        calcFractionalDissimilarities( s, d );
        return d;
    }

    private void calcKimuraDistances( final int s, final BasicSymmetricalDistanceMatrix d ) {
        for( int i = 1; i < s; i++ ) {
            for( int j = 0; j < i; j++ ) {
                d.setValue( i, j, calcKimuraDistance( i, j ) );
            }
        }
    }

    private void calcPoissonDistances( final int s, final BasicSymmetricalDistanceMatrix d ) {
        for( int i = 1; i < s; i++ ) {
            for( int j = 0; j < i; j++ ) {
                d.setValue( i, j, calcPoissonDistance( i, j ) );
            }
        }
    }

    private void calcFractionalDissimilarities( final int s, final BasicSymmetricalDistanceMatrix d ) {
        for( int i = 1; i < s; i++ ) {
            for( int j = 0; j < i; j++ ) {
                d.setValue( i, j, calcFractionalDissimilarity( i, j ) );
            }
        }
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();
    }

    private void copyIdentifiers( final int s, final BasicSymmetricalDistanceMatrix d ) {
        for( int i = 0; i < s; i++ ) {
            d.setIdentifier( i, _msa.getIdentifier( i ) );
        }
    }

    public static BasicSymmetricalDistanceMatrix calcFractionalDissimilarities( final Msa msa ) {
        return new PairwiseDistanceCalculator( msa, DEFAULT_VALUE_FOR_TOO_LARGE_DISTANCE_FOR_KIMURA_FORMULA )
        .calcFractionalDissimilarities();
    }

    public static BasicSymmetricalDistanceMatrix calcPoissonDistances( final Msa msa ) {
        return new PairwiseDistanceCalculator( msa, DEFAULT_VALUE_FOR_TOO_LARGE_DISTANCE_FOR_KIMURA_FORMULA )
        .calcPoissonDistances();
    }

    public static BasicSymmetricalDistanceMatrix calcKimuraDistances( final Msa msa ) {
        return new PairwiseDistanceCalculator( msa, DEFAULT_VALUE_FOR_TOO_LARGE_DISTANCE_FOR_KIMURA_FORMULA )
        .calcKimuraDistances();
    }

    public static BasicSymmetricalDistanceMatrix calcKimuraDistances( final Msa msa,
                                                                      final double value_for_too_large_distance_for_kimura_formula ) {
        return new PairwiseDistanceCalculator( msa, value_for_too_large_distance_for_kimura_formula )
        .calcKimuraDistances();
    }

    public enum PWD_DISTANCE_METHOD {
        KIMURA_DISTANCE, POISSON_DISTANCE, FRACTIONAL_DISSIMILARITY;
    }
}
