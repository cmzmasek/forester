
package org.forester.msa_compactor;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.msa.Mafft;
import org.forester.msa.Msa;
import org.forester.msa.Msa.MSA_FORMAT;
import org.forester.msa.MsaInferrer;
import org.forester.msa.MsaMethods;
import org.forester.sequence.Sequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class MsaCompactor {

    final private static NumberFormat NF_3    = new DecimalFormat( "#.###" );
    final private static NumberFormat NF_4    = new DecimalFormat( "#.####" );
    private static final boolean      VERBOSE = true;
    private Msa                       _msa;
    private final SortedSet<String>   _removed_seq_ids;
    static {
        NF_4.setRoundingMode( RoundingMode.HALF_UP );
        NF_3.setRoundingMode( RoundingMode.HALF_UP );
    }

    private MsaCompactor( final Msa msa ) {
        _msa = msa;
        _removed_seq_ids = new TreeSet<String>();
    }

    final public Msa getMsa() {
        return _msa;
    }

    final public SortedSet<String> getRemovedSeqIds() {
        return _removed_seq_ids;
    }

    final public void writeMsa( final File outfile, final MSA_FORMAT format, final String suffix ) throws IOException {
        final Double gr = MsaMethods.calcGapRatio( _msa );
        writeMsa( outfile + "_" + _msa.getNumberOfSequences() + "_" + _msa.getLength() + "_"
                          + ForesterUtil.roundToInt( gr * 100 ) + suffix,
                  format );
    }

    final int calcNonGapResidues( final Sequence seq ) {
        int ng = 0;
        for( int i = 0; i < seq.getLength(); ++i ) {
            if ( !seq.isGapAt( i ) ) {
                ++ng;
            }
        }
        return ng;
    }

    private final DescriptiveStatistics[] calcGapContribtionsX( final boolean normalize_for_effective_seq_length ) {
        final double gappiness[] = calcGappiness();
        final DescriptiveStatistics stats[] = new DescriptiveStatistics[ _msa.getNumberOfSequences() ];
        for( int row = 0; row < _msa.getNumberOfSequences(); ++row ) {
            stats[ row ] = new BasicDescriptiveStatistics( _msa.getIdentifier( row ) );
            final double l = calculateEffectiveLengthRatio( row );
            for( int col = 0; col < _msa.getLength(); ++col ) {
                if ( !_msa.isGapAt( row, col ) ) {
                    if ( normalize_for_effective_seq_length ) {
                        stats[ row ].addValue( gappiness[ col ] / l );
                    }
                    else {
                        stats[ row ].addValue( gappiness[ col ] );
                    }
                }
            }
        }
        return stats;
    }

    private final GapContribution[] calcGapContribtions( final boolean normalize_for_effective_seq_length ) {
        final double gappiness[] = calcGappiness();
        final GapContribution stats[] = new GapContribution[ _msa.getNumberOfSequences() ];
        for( int row = 0; row < _msa.getNumberOfSequences(); ++row ) {
            stats[ row ] = new GapContribution( _msa.getIdentifier( row ) );
            for( int col = 0; col < _msa.getLength(); ++col ) {
                if ( !_msa.isGapAt( row, col ) ) {
                    stats[ row ].addToValue( gappiness[ col ] );
                }
            }
            if ( normalize_for_effective_seq_length ) {
                stats[ row ].divideValue( calculateEffectiveLengthRatio( row ) );
            }
            else {
                // 
            }
        }
        return stats;
    }

    final private GapContribution[] calcGapContribtionsStats( final boolean norm ) {
        final GapContribution stats[] = calcGapContribtions( norm );
        Arrays.sort( stats );
        for( final GapContribution stat : stats ) {
            final StringBuilder sb = new StringBuilder();
            sb.append( stat.getId() );
            sb.append( "\t" );
            sb.append( NF_4.format( stat.getValue() ) );
            sb.append( "\t" );
            //            sb.append( NF_4.format( stat.median() ) );
            //            sb.append( "\t" );
            //            sb.append( NF_4.format( stat.getMin() ) );
            //            sb.append( "\t" );
            //            sb.append( NF_4.format( stat.getMax() ) );
            //sb.append( "\t" );
            System.out.println( sb );
        }
        return stats;
    }

    private final double[] calcGappiness() {
        final int l = _msa.getLength();
        final double gappiness[] = new double[ l ];
        final int seqs = _msa.getNumberOfSequences();
        for( int i = 0; i < l; ++i ) {
            gappiness[ i ] = ( double ) MsaMethods.calcGapSumPerColumn( _msa, i ) / seqs;
        }
        return gappiness;
    }

    private double calculateEffectiveLengthRatio( final int row ) {
        return ( double ) calcNonGapResidues( _msa.getSequence( row ) ) / _msa.getLength();
    }

    final private void mafft() throws IOException, InterruptedException {
        final MsaInferrer mafft = Mafft
                .createInstance( "/home/czmasek/SOFTWARE/MSA/MAFFT/mafft-7.130-without-extensions/scripts/mafft" );
        final List<String> opts = new ArrayList<String>();
        // opts.add( "--maxiterate" );
        // opts.add( "1000" );
        // opts.add( "--localpair" );
        opts.add( "--quiet" );
        _msa = mafft.infer( _msa.asSequenceList(), opts );
    }

    private StringBuilder msaStatsAsSB() {
        final StringBuilder sb = new StringBuilder();
        sb.append( _msa.getNumberOfSequences() );
        sb.append( "\t" );
        sb.append( _msa.getLength() );
        sb.append( "\t" );
        sb.append( NF_3.format( MsaMethods.calcGapRatio( _msa ) ) );
        return sb;
    }

    final private void removeGapColumns() {
        _msa = MsaMethods.createInstance().removeGapColumns( 1, 0, _msa );
    }

    final private void removeViaGapAverage( final double mean_gapiness,
                                            final int step,
                                            final boolean realign,
                                            final File outfile,
                                            final int minimal_effective_length ) throws IOException,
            InterruptedException {
        if ( step < 1 ) {
            throw new IllegalArgumentException( "step cannot be less than 1" );
        }
        if ( mean_gapiness < 0 ) {
            throw new IllegalArgumentException( "target average gap ratio cannot be less than 0" );
        }
        if ( VERBOSE ) {
            System.out.println( "orig: " + msaStatsAsSB() );
        }
        if ( minimal_effective_length > 1 ) {
            _msa = MsaMethods.removeSequencesByMinimalLength( _msa, minimal_effective_length );
            if ( VERBOSE ) {
                System.out.println( "short seq removal: " + msaStatsAsSB() );
            }
        }
        int counter = step;
        double gr;
        do {
            removeWorstOffenders( step, 1, false, false );
            if ( realign ) {
                mafft();
            }
            gr = MsaMethods.calcGapRatio( _msa );
            if ( VERBOSE ) {
                System.out.println( counter + ": " + msaStatsAsSB() );
            }
            //   write( outfile, gr );
            counter += step;
        } while ( gr > mean_gapiness );
        if ( VERBOSE ) {
            System.out.println( "final: " + msaStatsAsSB() );
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
            System.out.println( "orig: " + msaStatsAsSB() );
        }
        int counter = step;
        while ( _msa.getLength() > length ) {
            removeWorstOffenders( step, 1, false, false );
            if ( realign ) {
                mafft();
            }
            if ( VERBOSE ) {
                System.out.println( counter + ": " + msaStatsAsSB() );
            }
            counter += step;
        }
    }

    final private void removeWorstOffenders( final int to_remove,
                                             final int step,
                                             final boolean realign,
                                             final boolean norm ) throws IOException, InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        for( int j = 0; j < to_remove; ++j ) {
            to_remove_ids.add( stats[ j ].getId() );
            _removed_seq_ids.add( stats[ j ].getId() );
        }
        //TODO if verbose/interactve
        for( final String id : to_remove_ids ) {
            _msa = MsaMethods.removeSequence( _msa, id );
            removeGapColumns();
            System.out.print( id );
            System.out.print( "\t" );
            final StringBuilder sb = msaStatsAsSB();
            System.out.println( sb );
        }
        //TODO else:
        //_msa = MsaMethods.removeSequences( _msa, to_remove_ids );
        //removeGapColumns();
        if ( realign ) {
            mafft();
        }
    }

    final private void writeMsa( final String outfile, final MSA_FORMAT format ) throws IOException {
        final Writer w = ForesterUtil.createBufferedWriter( outfile );
        _msa.write( w, format );
        w.close();
    }

    public final static MsaCompactor reduceGapAverage( final Msa msa,
                                                       final double max_gap_average,
                                                       final int step,
                                                       final boolean realign,
                                                       final File out,
                                                       final int minimal_effective_length ) throws IOException,
            InterruptedException {
        final MsaCompactor mc = new MsaCompactor( msa );
        mc.removeViaGapAverage( max_gap_average, step, realign, out, minimal_effective_length );
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

    public final static MsaCompactor removeWorstOffenders( final Msa msa,
                                                           final int worst_offenders_to_remove,
                                                           final boolean realign,
                                                           final boolean norm ) throws IOException,
            InterruptedException {
        final MsaCompactor mc = new MsaCompactor( msa );
        mc.removeWorstOffenders( worst_offenders_to_remove, 1, realign, norm );
        return mc;
    }

    public static enum SORT_BY {
        MAX, MEAN, MEDIAN;
    }

    final static class DescriptiveStatisticsComparator implements Comparator<DescriptiveStatistics> {

        final private boolean _ascending;
        final private SORT_BY _sort_by;

        public DescriptiveStatisticsComparator( final boolean ascending, final SORT_BY sort_by ) {
            _ascending = ascending;
            _sort_by = sort_by;
        }

        @Override
        public final int compare( final DescriptiveStatistics s0, final DescriptiveStatistics s1 ) {
            switch ( _sort_by ) {
                case MAX:
                    if ( s0.getMax() < s1.getMax() ) {
                        return _ascending ? -1 : 1;
                    }
                    else if ( s0.getMax() > s1.getMax() ) {
                        return _ascending ? 1 : -1;
                    }
                    return 0;
                case MEAN:
                    if ( s0.arithmeticMean() < s1.arithmeticMean() ) {
                        return _ascending ? -1 : 1;
                    }
                    else if ( s0.arithmeticMean() > s1.arithmeticMean() ) {
                        return _ascending ? 1 : -1;
                    }
                    return 0;
                case MEDIAN:
                    if ( s0.median() < s1.median() ) {
                        return _ascending ? -1 : 1;
                    }
                    else if ( s0.median() > s1.median() ) {
                        return _ascending ? 1 : -1;
                    }
                    return 0;
                default:
                    return 0;
            }
        }
    }
}
