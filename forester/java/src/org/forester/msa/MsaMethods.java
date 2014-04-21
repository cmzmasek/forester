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

package org.forester.msa;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.sequence.BasicSequence;
import org.forester.sequence.Sequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;

public final class MsaMethods {

    private ArrayList<String> _ignored_seqs_ids;

    synchronized public ArrayList<String> getIgnoredSequenceIds() {
        return _ignored_seqs_ids;
    }

    synchronized public static MsaMethods createInstance() {
        return new MsaMethods();
    }

    private MsaMethods() {
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
            if ( msa.isGapAt( j, col ) ) {
                gap_rows++;
            }
        }
        return gap_rows;
    }

    final public static Msa removeSequence( final Msa msa, final String to_remove_id ) {
        final List<Sequence> seqs = new ArrayList<Sequence>();
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            if ( !to_remove_id.equals( msa.getIdentifier( row ) ) ) {
                seqs.add( msa.getSequence( row ) );
            }
        }
        if ( seqs.size() < 1 ) {
            return null;
        }
        return BasicMsa.createInstance( seqs );
    }

    final public static Msa removeSequences( final Msa msa, final List<String> to_remove_ids ) {
        final List<Sequence> seqs = new ArrayList<Sequence>();
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            if ( !to_remove_ids.contains( msa.getIdentifier( row ) ) ) {
                seqs.add( msa.getSequence( row ) );
            }
        }
        if ( seqs.size() < 1 ) {
            return null;
        }
        return BasicMsa.createInstance( seqs );
    }

    final public static Msa removeSequencesByRow( final Msa msa, final List<Integer> to_remove_rows ) {
        final List<Sequence> seqs = new ArrayList<Sequence>();
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            if ( !to_remove_rows.contains( row ) ) {
                seqs.add( msa.getSequence( row ) );
            }
        }
        if ( seqs.size() < 1 ) {
            return null;
        }
        return BasicMsa.createInstance( seqs );
    }

    synchronized final public Msa removeGapColumns( final double max_allowed_gap_ratio,
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
            delete_cols[ col ] = ( ( double ) calcGapSumPerColumn( msa, col ) / msa.getNumberOfSequences() ) >= max_allowed_gap_ratio;
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

    synchronized final public static void removeGapColumns( final double max_allowed_gap_ratio, final DeleteableMsa msa ) {
        if ( ( max_allowed_gap_ratio < 0 ) || ( max_allowed_gap_ratio > 1 ) ) {
            throw new IllegalArgumentException( "max allowed gap ration is out of range: " + max_allowed_gap_ratio );
        }
        //   final boolean ignore_too_short_seqs = min_allowed_length > 0;
        for( int col = 0; col < msa.getLength(); ++col ) {
            final boolean delete = ( ( double ) calcGapSumPerColumn( msa, col ) / msa.getNumberOfSequences() ) >= max_allowed_gap_ratio;
            if ( delete ) {
                msa.deleteColumn( col );
            }
        }
    }

    public static DescriptiveStatistics calculateIdentityRatio( final int from, final int to, final Msa msa ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int c = from; c <= to; ++c ) {
            stats.addValue( calculateIdentityRatio( msa, c ) );
        }
        return stats;
    }

    public static double calculateIdentityRatio( final Msa msa, final int column ) {
        final SortedMap<Character, Integer> dist = calculateResidueDestributionPerColumn( msa, column );
        int majority_count = 0;
        final Iterator<Map.Entry<Character, Integer>> it = dist.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<Character, Integer> pair = it.next();
            if ( pair.getValue() > majority_count ) {
                majority_count = pair.getValue();
            }
        }
        return ( double ) majority_count / msa.getNumberOfSequences();
    }

    public static SortedMap<Character, Integer> calculateResidueDestributionPerColumn( final Msa msa, final int column ) {
        final SortedMap<Character, Integer> map = new TreeMap<Character, Integer>();
        for( final Character r : msa.getColumnAt( column ) ) {
            if ( r != Sequence.GAP ) {
                if ( !map.containsKey( r ) ) {
                    map.put( r, 1 );
                }
                else {
                    map.put( r, map.get( r ) + 1 );
                }
            }
        }
        return map;
    }

    public static DescriptiveStatistics calcBasicGapinessStatistics( final Msa msa ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int i = 0; i < msa.getLength(); ++i ) {
            stats.addValue( ( double ) calcGapSumPerColumn( msa, i ) / msa.getNumberOfSequences() );
        }
        return stats;
    }

    public static Msa removeSequencesByMinimalLength( final Msa msa, final int min_effective_length ) {
        final List<Integer> to_remove_rows = new ArrayList<Integer>();
        for( int seq = 0; seq < msa.getNumberOfSequences(); ++seq ) {
            int eff_length = 0;
            for( int i = 0; i < msa.getLength(); ++i ) {
                if ( msa.getResidueAt( seq, i ) != Sequence.GAP ) {
                    eff_length++;
                }
            }
            if ( eff_length < min_effective_length ) {
                to_remove_rows.add( seq );
            }
        }
        return removeSequencesByRow( msa, to_remove_rows );
    }

    public static double calcGapRatio( final Msa msa ) {
        int gaps = 0;
        for( int seq = 0; seq < msa.getNumberOfSequences(); ++seq ) {
            for( int i = 0; i < msa.getLength(); ++i ) {
                if ( msa.getResidueAt( seq, i ) == Sequence.GAP ) {
                    gaps++;
                }
            }
        }
        return ( double ) gaps / ( msa.getLength() * msa.getNumberOfSequences() );
    }
}
