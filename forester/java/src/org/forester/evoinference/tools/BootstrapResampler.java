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

package org.forester.evoinference.tools;

import java.util.Random;

import org.forester.msa.BasicMsa;
import org.forester.msa.Msa;

public class BootstrapResampler {

    private static void copyIdentifiers( final Msa msa, final Msa new_msa ) {
        for( int i = 0; i < msa.getNumberOfSequences(); ++i ) {
            new_msa.setIdentifier( i, msa.getIdentifier( i ) );
        }
    }

    private static void preconditionCheck( final Msa msa, final int n ) {
        if ( msa.getLength() < 2 ) {
            throw new IllegalArgumentException( "Msa length cannot be smaller than two for bootstrap resampling" );
        }
        if ( msa.getNumberOfSequences() < 1 ) {
            throw new IllegalArgumentException( "Attempt to bootstrap resample empty multiple sequence alignment" );
        }
        if ( n < 1 ) {
            throw new IllegalArgumentException( "Number of bootstrap resamples cannot be zero or negative" );
        }
    }

    private static void preconditionCheck( final int length, final int n ) {
        if ( length < 2 ) {
            throw new IllegalArgumentException( "Msa length cannot be smaller than two for bootstrap resampling" );
        }
        if ( n < 1 ) {
            throw new IllegalArgumentException( "Number of bootstrap resamples cannot be zero or negative" );
        }
    }

    public static Msa[] resample( final Msa msa, final int n, final long seed ) {
        preconditionCheck( msa, n );
        final Random random = new Random( seed );
        final Msa[] msas = new Msa[ n ];
        for( int i = 0; i < n; ++i ) {
            final Msa new_msa = new BasicMsa( msa.getNumberOfSequences(), msa.getLength(), msa.getType() );
            msas[ i ] = new_msa;
            copyIdentifiers( msa, new_msa );
            for( int col = 0; col < msa.getLength(); ++col ) {
                final int random_col = random.nextInt( msa.getLength() );
                for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                    new_msa.setResidueAt( row, col, msa.getResidueAt( row, random_col ) );
                }
            }
        }
        return msas;
    }

    public static int[][] createResampledColumnPositions( final int length, final int n, final long seed ) {
        preconditionCheck( length, n );
        final Random random = new Random( seed );
        final int[][] columns = new int[ n ][ length ];
        for( int i = 0; i < n; ++i ) {
            for( int col = 0; col < length; ++col ) {
                columns[ i ][ col ] = random.nextInt( length );
            }
        }
        return columns;
    }
}
