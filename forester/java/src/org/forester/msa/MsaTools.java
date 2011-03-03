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
// WWW: www.phylosoft.org/forester

package org.forester.msa;

import java.util.ArrayList;
import java.util.List;

import org.forester.sequence.BasicSequence;
import org.forester.sequence.Sequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;

public final class MsaTools {

    private ArrayList<String> _ignored_seqs_ids;

    synchronized public ArrayList<String> getIgnoredSequenceIds() {
        return _ignored_seqs_ids;
    }

    synchronized public static MsaTools createInstance() {
        return new MsaTools();
    }

    private MsaTools() {
        init();
    }

    synchronized private void init() {
        _ignored_seqs_ids = new ArrayList<String>();
    }

    @Override
    public Object clone() {
        throw new NoSuchMethodError();
    }

    public static int calcGapSumPerColumn( final Msa msa, final int col ) {
        int gap_rows = 0;
        for( int j = 0; j < msa.getNumberOfSequences(); ++j ) {
            if ( msa.getResidueAt( j, col ) == Sequence.GAP ) {
                gap_rows++;
            }
        }
        return gap_rows;
    }

    synchronized public Msa removeGapColumns( final double max_allowed_gap_ratio,
                                              final int min_allowed_length,
                                              final Msa msa ) {
        init();
        if ( ( max_allowed_gap_ratio < 0 ) || ( max_allowed_gap_ratio > 1 ) ) {
            throw new IllegalArgumentException( "max allowed gap ration is out of range: " + max_allowed_gap_ratio );
        }
        final boolean ignore_too_short_seqs = min_allowed_length > 0;
        final boolean[] delete_cols = new boolean[ msa.getLength() ];
        int new_length = 0;
        for( int col = 0; col < msa.getLength(); ++col ) {
            delete_cols[ col ] = ( ( double ) calcGapSumPerColumn( msa, col ) / msa.getNumberOfSequences() ) > max_allowed_gap_ratio;
            if ( !delete_cols[ col ] ) {
                ++new_length;
            }
        }
        final List<Sequence> seqs = new ArrayList<Sequence>( msa.getNumberOfSequences() );
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            final char[] mol_seq = new char[ new_length ];
            int new_col = 0;
            int non_gap_cols_sum = 0;
            for( int col = 0; col < msa.getLength(); ++col ) {
                if ( !delete_cols[ col ] ) {
                    final char residue = msa.getResidueAt( row, col );
                    mol_seq[ new_col++ ] = ( residue );
                    if ( residue != Sequence.GAP ) {
                        ++non_gap_cols_sum;
                    }
                }
            }
            if ( ignore_too_short_seqs ) {
                if ( non_gap_cols_sum >= min_allowed_length ) {
                    seqs.add( new BasicSequence( msa.getIdentifier( row ), mol_seq, msa.getType() ) );
                }
                else {
                    _ignored_seqs_ids.add( msa.getIdentifier( row ).toString() );
                }
            }
            else {
                seqs.add( new BasicSequence( msa.getIdentifier( row ), mol_seq, msa.getType() ) );
            }
        }
        if ( seqs.size() < 1 ) {
            return null;
        }
        return BasicMsa.createInstance( seqs );
    }

    public static DescriptiveStatistics calcBasicGapinessStatistics( final Msa msa ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int i = 0; i < msa.getLength(); ++i ) {
            stats.addValue( ( double ) calcGapSumPerColumn( msa, i ) / msa.getNumberOfSequences() );
        }
        return stats;
    }
}
