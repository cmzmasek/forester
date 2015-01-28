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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;

public final class MsaMethods {

    private ArrayList<String> _ignored_seqs_ids;

    private MsaMethods() {
        init();
    }

    @Override
    public Object clone() {
        throw new NoSuchMethodError();
    }

    synchronized final public Msa deleteGapColumns( final double max_allowed_gap_ratio,
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
        final List<MolecularSequence> seqs = new ArrayList<MolecularSequence>( msa.getNumberOfSequences() );
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            final char[] mol_seq = new char[ new_length ];
            int new_col = 0;
            int non_gap_cols_sum = 0;
            for( int col = 0; col < msa.getLength(); ++col ) {
                if ( !delete_cols[ col ] ) {
                    final char residue = msa.getResidueAt( row, col );
                    mol_seq[ new_col++ ] = ( residue );
                    if ( residue != MolecularSequence.GAP ) {
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

    synchronized public ArrayList<String> getIgnoredSequenceIds() {
        return _ignored_seqs_ids;
    }

    synchronized private void init() {
        _ignored_seqs_ids = new ArrayList<String>();
    }

    public static final DescriptiveStatistics calcNumberOfGapsPer100Stats( final Msa msa ) {
        final int[] gaps = calcNumberOfGapsInMsa( msa );
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        final double n = 100.0 / msa.getLength();
        for( final int gap : gaps ) {
            stats.addValue( n * gap );
        }
        return stats;
    }

    public static final int[] calcNumberOfGapsInMsa( final Msa msa ) {
        final int seqs = msa.getNumberOfSequences();
        final int[]  gaps= new int[ seqs ];
        for( int i = 0; i < seqs; ++i ) {
            gaps[ i ] =  calcNumberOfGaps( msa.getSequence( i ) );
        }
        return gaps;
    }
    
    

    public final static int calcNumberOfGaps( final MolecularSequence seq  ) {
        int gaps = 0;
        boolean was_gap = false;
        for( int i = 0; i < seq.getLength(); ++i ) {
            if ( seq.isGapAt( i ) ) {
               if ( !was_gap ) {
                   ++gaps;
                   was_gap = true;
               }
            }
            else {
                was_gap = false;
            }
        }
        return gaps;
    }

    public static DescriptiveStatistics calcBasicGapinessStatistics( final Msa msa ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int i = 0; i < msa.getLength(); ++i ) {
            stats.addValue( ( double ) calcGapSumPerColumn( msa, i ) / msa.getNumberOfSequences() );
        }
        return stats;
    }

    public static double calcGapRatio( final Msa msa ) {
        int gaps = 0;
        for( int seq = 0; seq < msa.getNumberOfSequences(); ++seq ) {
            for( int i = 0; i < msa.getLength(); ++i ) {
                if ( msa.getResidueAt( seq, i ) == MolecularSequence.GAP ) {
                    gaps++;
                }
            }
        }
        return ( double ) gaps / ( msa.getLength() * msa.getNumberOfSequences() );
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

    final public static double calcNormalizedShannonsEntropy( final int k, final Msa msa ) {
        double s = 0;
        for( int col = 0; col < msa.getLength(); ++col ) {
            s += calcNormalizedShannonsEntropy( k, msa, col );
        }
        return s / msa.getLength();
    }

    final public static double calcNormalizedShannonsEntropy( final int k, final Msa msa, final int col ) {
        // http://www.ebi.ac.uk/thornton-srv/databases/valdarprograms/scorecons_server_help.html
        // n: number of residues in column
        // k: number of residue types
        // na: number of residues of type a
        // pa = na/n
        // S=-sum pa log2 pa
        double s = 0;
        final double n = msa.getNumberOfSequences();
        HashMap<Character, Integer> dist = null;
        if ( k == 6 ) {
            dist = calcResidueDistribution6( msa, col );
        }
        else if ( k == 7 ) {
            dist = calcResidueDistribution7( msa, col );
        }
        else if ( k == 20 ) {
            dist = calcResidueDistribution20( msa, col );
        }
        else if ( k == 21 ) {
            dist = calcResidueDistribution21( msa, col );
        }
        else {
            throw new IllegalArgumentException( "illegal value for k: " + k );
        }
        if ( dist.size() == 1 ) {
            return 0;
        }
        //        if ( dist.size() == n ) {
        //            return 0;
        //        }
        for( final int na : dist.values() ) {
            final double pa = na / n;
            s += pa * Math.log( pa );
        }
        if ( n < k ) {
            return -( s / ( Math.log( n ) ) );
        }
        else {
            return -( s / ( Math.log( k ) ) );
        }
    }

    final public static DescriptiveStatistics calculateEffectiveLengthStatistics( final Msa msa ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            final MolecularSequence s = msa.getSequence( row );
            stats.addValue( s.getLength() - s.getNumberOfGapResidues() );
        }
        return stats;
    }

    final public static DescriptiveStatistics calculateIdentityRatio( final int from, final int to, final Msa msa ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int c = from; c <= to; ++c ) {
            stats.addValue( calculateIdentityRatio( msa, c ) );
        }
        return stats;
    }

    final public static double calculateIdentityRatio( final Msa msa, final int column ) {
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
            if ( r != MolecularSequence.GAP ) {
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

    synchronized public static MsaMethods createInstance() {
        return new MsaMethods();
    }

    final public static Msa removeSequence( final Msa msa, final String to_remove_id ) {
        final List<MolecularSequence> seqs = new ArrayList<MolecularSequence>();
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
        final List<MolecularSequence> seqs = new ArrayList<MolecularSequence>();
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

    public static Msa removeSequencesByMinimalLength( final Msa msa, final int min_effective_length ) {
        final List<Integer> to_remove_rows = new ArrayList<Integer>();
        for( int seq = 0; seq < msa.getNumberOfSequences(); ++seq ) {
            int eff_length = 0;
            for( int i = 0; i < msa.getLength(); ++i ) {
                if ( msa.getResidueAt( seq, i ) != MolecularSequence.GAP ) {
                    eff_length++;
                }
            }
            if ( eff_length < min_effective_length ) {
                to_remove_rows.add( seq );
            }
        }
        return removeSequencesByRow( msa, to_remove_rows );
    }

    final public static Msa removeSequencesByRow( final Msa msa, final List<Integer> to_remove_rows ) {
        final List<MolecularSequence> seqs = new ArrayList<MolecularSequence>();
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

    final private static HashMap<Character, Integer> calcResidueDistribution20( final Msa msa, final int col ) {
        final HashMap<Character, Integer> counts = new HashMap<Character, Integer>();
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            final char c = msa.getResidueAt( row, col );
            if ( c != MolecularSequence.GAP ) {
                if ( !counts.containsKey( c ) ) {
                    counts.put( c, 1 );
                }
                else {
                    counts.put( c, 1 + counts.get( c ) );
                }
            }
        }
        return counts;
    }

    final private static HashMap<Character, Integer> calcResidueDistribution21( final Msa msa, final int col ) {
        final HashMap<Character, Integer> counts = new HashMap<Character, Integer>();
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            final char c = msa.getResidueAt( row, col );
            if ( !counts.containsKey( c ) ) {
                counts.put( c, 1 );
            }
            else {
                counts.put( c, 1 + counts.get( c ) );
            }
        }
        return counts;
    }

    final private static HashMap<Character, Integer> calcResidueDistribution6( final Msa msa, final int col ) {
        // Residues are classified into one of tex2html_wrap199 types:
        // aliphatic [AVLIMC], aromatic [FWYH], polar [STNQ], positive [KR], negative [DE],
        // special conformations [GP] and gaps. This convention follows that
        // of Mirny & Shakhnovich (1999, J Mol Biol 291:177-196).
        final HashMap<Character, Integer> counts = new HashMap<Character, Integer>();
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            final char c = msa.getResidueAt( row, col );
            char x;
            if ( ( c == 'A' ) || ( c == 'V' ) || ( c == 'L' ) || ( c == 'I' ) || ( c == 'M' ) || ( c == 'C' ) ) {
                // aliphatic
                x = 'a';
            }
            else if ( ( c == 'F' ) || ( c == 'W' ) || ( c == 'Y' ) || ( c == 'H' ) ) {
                // aromatic
                x = 'r';
            }
            else if ( ( c == 'S' ) || ( c == 'T' ) || ( c == 'N' ) || ( c == 'Q' ) ) {
                // polar
                x = 'p';
            }
            else if ( ( c == 'K' ) || ( c == 'R' ) ) {
                // positive
                x = 'o';
            }
            else if ( ( c == 'D' ) || ( c == 'E' ) ) {
                // negative
                x = 'e';
            }
            else if ( ( c == 'G' ) || ( c == 'P' ) ) {
                // aliphatic - special conformation
                x = 's';
            }
            else {
                continue;
            }
            if ( !counts.containsKey( x ) ) {
                counts.put( x, 1 );
            }
            else {
                counts.put( x, 1 + counts.get( x ) );
            }
        }
        return counts;
    }

    final private static HashMap<Character, Integer> calcResidueDistribution7( final Msa msa, final int col ) {
        // Residues are classified into one of tex2html_wrap199 types:
        // aliphatic [AVLIMC], aromatic [FWYH], polar [STNQ], positive [KR], negative [DE],
        // special conformations [GP] and gaps. This convention follows that
        // of Mirny & Shakhnovich (1999, J Mol Biol 291:177-196).
        final HashMap<Character, Integer> counts = new HashMap<Character, Integer>();
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            final char c = msa.getResidueAt( row, col );
            char x = '-';
            if ( ( c == 'A' ) || ( c == 'V' ) || ( c == 'L' ) || ( c == 'I' ) || ( c == 'M' ) || ( c == 'C' ) ) {
                // aliphatic
                x = 'a';
            }
            else if ( ( c == 'F' ) || ( c == 'W' ) || ( c == 'Y' ) || ( c == 'H' ) ) {
                // aromatic
                x = 'r';
            }
            else if ( ( c == 'S' ) || ( c == 'T' ) || ( c == 'N' ) || ( c == 'Q' ) ) {
                // polar
                x = 'p';
            }
            else if ( ( c == 'K' ) || ( c == 'R' ) ) {
                // positive
                x = 'o';
            }
            else if ( ( c == 'D' ) || ( c == 'E' ) ) {
                // negative
                x = 'e';
            }
            else if ( ( c == 'G' ) || ( c == 'P' ) ) {
                // aliphatic - special conformation
                x = 's';
            }
            if ( !counts.containsKey( x ) ) {
                counts.put( x, 1 );
            }
            else {
                counts.put( x, 1 + counts.get( x ) );
            }
        }
        return counts;
    }
}
