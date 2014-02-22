
package org.forester.msa;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.msa.Msa.MSA_FORMAT;
import org.forester.sequence.Sequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class MsaCompactor {

    private static final boolean VERBOSE = true;

    public static enum SORT_BY {
        MAX, MEAN, MEDIAN;
    }
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
        final MsaInferrer mafft = Mafft.createInstance( "mafft" );
        final List<String> opts = new ArrayList<String>();
        // opts.add( "--maxiterate" );
        // opts.add( "1000" );
        // opts.add( "--localpair" );
        opts.add( "--quiet" );
        _msa = mafft.infer( _msa.asSequenceList(), opts );
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
            removeWorstOffenders( step, 1, false );
            if ( realign ) {
                mafft();
            }
            gr = MsaMethods.calcGapRatio( _msa );
            if ( VERBOSE ) {
                System.out.println( counter + ": " + msaStatsAsSB() );
            }
            write( outfile, gr );
            counter += step;
        } while ( gr > mean_gapiness );
        if ( VERBOSE ) {
            System.out.println( "final: " + msaStatsAsSB() );
        }
    }

    final private void write( final File outfile, final double gr ) throws IOException {
        writeMsa( outfile + "_" + _msa.getNumberOfSequences() + "_" + _msa.getLength() + "_"
                + ForesterUtil.roundToInt( gr * 100 ) + ".fasta" );
    }

    final private void writeMsa( final String outfile ) throws IOException {
        final Writer w = ForesterUtil.createBufferedWriter( outfile );
        _msa.write( w, MSA_FORMAT.FASTA );
        w.close();
    }

    final private StringBuilder msaStatsAsSB() {
        final StringBuilder sb = new StringBuilder();
        sb.append( _msa.getLength() );
        sb.append( "\t" );
        sb.append( _msa.getNumberOfSequences() );
        sb.append( "\t" );
        sb.append( ForesterUtil.round( MsaMethods.calcGapRatio( _msa ), 4 ) );
        sb.append( "\t" );
        return sb;
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
            removeWorstOffenders( step, 1, false );
            if ( realign ) {
                mafft();
            }
            if ( VERBOSE ) {
                System.out.println( counter + ": " + msaStatsAsSB() );
            }
            counter += step;
        }
    }

    final private DescriptiveStatistics[] calcStats() {
        final DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator( '.' );
        final NumberFormat f = new DecimalFormat( "#.####", dfs );
        f.setRoundingMode( RoundingMode.HALF_UP );
        final DescriptiveStatistics stats[] = calcGapContribtions();
        Arrays.sort( stats, new DescriptiveStatisticsComparator( false, SORT_BY.MEAN ) );
        for( final DescriptiveStatistics stat : stats ) {
            final StringBuilder sb = new StringBuilder();
            sb.append( stat.getDescription() );
            sb.append( "\t" );
            sb.append( f.format( stat.arithmeticMean() ) );
            sb.append( "\t" );
            sb.append( f.format( stat.median() ) );
            sb.append( "\t" );
            sb.append( f.format( stat.getMin() ) );
            sb.append( "\t" );
            sb.append( f.format( stat.getMax() ) );
            sb.append( "\t" );
            System.out.println( sb );
        }
        return stats;
    }

    private final DescriptiveStatistics[] calcGapContribtions() {
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
        final int l = _msa.getLength();
        final double gappiness[] = new double[ l ];
        final int seqs = _msa.getNumberOfSequences();
        for( int i = 0; i < l; ++i ) {
            gappiness[ i ] = ( double ) MsaMethods.calcGapSumPerColumn( _msa, i ) / seqs;
        }
        return gappiness;
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
