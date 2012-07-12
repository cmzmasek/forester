
package org.forester.msa;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.sequence.Sequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class MsaCompactor {

    private static final boolean    VERBOSE = true;
    private Msa                     _msa;
    private final SortedSet<String> _removed_seq_ids;

    private MsaCompactor( final Msa msa ) {
        _msa = msa;
        _removed_seq_ids = new TreeSet<String>();
    }

    final public SortedSet<String> getRemovedSeqIds() {
        return _removed_seq_ids;
    }

    final public Msa getMsa() {
        return _msa;
    }

    public final static MsaCompactor removeWorstOffenders( final Msa msa,
                                                           final int worst_offenders_to_remove,
                                                           final boolean realign ) throws IOException,
            InterruptedException {
        final MsaCompactor mc = new MsaCompactor( msa );
        mc.removeWorstOffenders( worst_offenders_to_remove, 1, realign );
        return mc;
    }

    public final static MsaCompactor reduceGapAverage( final Msa msa,
                                                       final double max_gap_average,
                                                       final int step,
                                                       final boolean realign ) throws IOException, InterruptedException {
        final MsaCompactor mc = new MsaCompactor( msa );
        mc.removeViaGapAverage( max_gap_average, step, realign );
        return mc;
    }

    public final static MsaCompactor reduceLength( final Msa msa,
                                                   final int length,
                                                   final int step,
                                                   final boolean realign ) throws IOException, InterruptedException {
        final MsaCompactor mc = new MsaCompactor( msa );
        mc.removeViaLength( length, step, realign );
        return mc;
    }

    final private void removeGapColumns() {
        _msa = MsaMethods.createInstance().removeGapColumns( 1, 0, _msa );
    }

    final private void removeWorstOffenders( final int to_remove, final int step, final boolean realign )
            throws IOException, InterruptedException {
        final DescriptiveStatistics stats[] = calcStats();
        final List<String> to_remove_ids = new ArrayList<String>();
        for( int j = 0; j < to_remove; ++j ) {
            to_remove_ids.add( stats[ j ].getDescription() );
            _removed_seq_ids.add( stats[ j ].getDescription() );
        }
        _msa = MsaMethods.removeSequences( _msa, to_remove_ids );
        removeGapColumns();
        if ( realign ) {
            mafft();
        }
    }

    final private void mafft() throws IOException, InterruptedException {
        final MsaInferrer mafft = Mafft.createInstance( "/home/czmasek/bin/mafft" );
        final List<String> opts = new ArrayList<String>();
        opts.add( "--maxiterate" );
        opts.add( "1000" );
        opts.add( "--localpair" );
        opts.add( "--quiet" );
        _msa = mafft.infer( _msa.asSequenceList(), opts );
    }

    final private void removeViaGapAverage( final double mean_gapiness, final int step, final boolean realign )
            throws IOException, InterruptedException {
        if ( step < 1 ) {
            throw new IllegalArgumentException( "step cannot be less than 1" );
        }
        if ( mean_gapiness < 0 ) {
            throw new IllegalArgumentException( "target average gap ratio cannot be less than 0" );
        }
        if ( VERBOSE ) {
            System.out.println( "start: " + _msa.getLength() + " "
                    + ForesterUtil.round( MsaMethods.calcBasicGapinessStatistics( _msa ).arithmeticMean(), 3 ) );
        }
        int counter = step;
        while ( MsaMethods.calcBasicGapinessStatistics( _msa ).arithmeticMean() > mean_gapiness ) {
            removeWorstOffenders( step, 1, false );
            if ( realign ) {
                mafft();
            }
            if ( VERBOSE ) {
                System.out.println( counter + ": " + _msa.getLength() + " "
                        + ForesterUtil.round( MsaMethods.calcBasicGapinessStatistics( _msa ).arithmeticMean(), 3 ) );
            }
            counter += step;
        }
    }

    final private void removeViaLength( final int length, final int step, final boolean realign ) throws IOException,
            InterruptedException {
        if ( step < 1 ) {
            throw new IllegalArgumentException( "step cannot be less than 1" );
        }
        if ( length < 11 ) {
            throw new IllegalArgumentException( "target length cannot be less than 1" );
        }
        if ( VERBOSE ) {
            System.out.println( "start: " + _msa.getLength() + " "
                    + ForesterUtil.round( MsaMethods.calcBasicGapinessStatistics( _msa ).arithmeticMean(), 3 ) );
        }
        int counter = step;
        while ( _msa.getLength() > length ) {
            removeWorstOffenders( step, 1, false );
            if ( realign ) {
                mafft();
            }
            if ( VERBOSE ) {
                System.out.println( counter + ": " + _msa.getLength() + " "
                        + ForesterUtil.round( MsaMethods.calcBasicGapinessStatistics( _msa ).arithmeticMean(), 3 ) );
            }
            counter += step;
        }
    }

    final private DescriptiveStatistics[] calcStats() {
        final DescriptiveStatistics stats[] = calc();
        sort( stats );
        for( int i = 0; i < stats.length; ++i ) {
            final DescriptiveStatistics s = stats[ i ];
            //            System.out.print( s.getDescription() );
            //            System.out.print( "\t" );
            //            System.out.print( s.arithmeticMean() );
            //            System.out.print( "\t(" );
            //            System.out.print( s.arithmeticMean() );
            //            System.out.print( ")" );
            //            System.out.print( "\t" );
            //            System.out.print( s.getMin() );
            //            System.out.print( "\t" );
            //            System.out.print( s.getMax() );
            //            System.out.println();
        }
        return stats;
    }

    private final static void sort( final DescriptiveStatistics stats[] ) {
        Arrays.sort( stats, new DescriptiveStatisticsComparator( false ) );
    }

    private final DescriptiveStatistics[] calc() {
        final double gappiness[] = calcGappiness();
        final DescriptiveStatistics stats[] = new DescriptiveStatistics[ _msa.getNumberOfSequences() ];
        for( int row = 0; row < _msa.getNumberOfSequences(); ++row ) {
            stats[ row ] = new BasicDescriptiveStatistics( _msa.getIdentifier( row ) );
            for( int col = 0; col < _msa.getLength(); ++col ) {
                if ( _msa.getResidueAt( row, col ) != Sequence.GAP ) {
                    stats[ row ].addValue( gappiness[ col ] );
                }
            }
        }
        return stats;
    }

    private final double[] calcGappiness() {
        final double gappiness[] = new double[ _msa.getLength() ];
        final int seqs = _msa.getNumberOfSequences();
        for( int i = 0; i < gappiness.length; ++i ) {
            gappiness[ i ] = ( double ) MsaMethods.calcGapSumPerColumn( _msa, i ) / seqs;
        }
        return gappiness;
    }

    final static class DescriptiveStatisticsComparator implements Comparator<DescriptiveStatistics> {

        final boolean _ascending;

        public DescriptiveStatisticsComparator( final boolean ascending ) {
            _ascending = ascending;
        }

        @Override
        public final int compare( final DescriptiveStatistics s0, final DescriptiveStatistics s1 ) {
            if ( s0.arithmeticMean() < s1.arithmeticMean() ) {
                return _ascending ? -1 : 1;
            }
            else if ( s0.arithmeticMean() > s1.arithmeticMean() ) {
                return _ascending ? 1 : -1;
            }
            return 0;
        }
    }
}
