
package org.forester.msa_compactor;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.evoinference.distance.NeighborJoiningF;
import org.forester.evoinference.distance.PairwiseDistanceCalculator;
import org.forester.evoinference.distance.PairwiseDistanceCalculator.PWD_DISTANCE_METHOD;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.evoinference.tools.BootstrapResampler;
import org.forester.msa.BasicMsa;
import org.forester.msa.Mafft;
import org.forester.msa.Msa;
import org.forester.msa.Msa.MSA_FORMAT;
import org.forester.msa.MsaInferrer;
import org.forester.msa.MsaMethods;
import org.forester.msa.ResampleableMsa;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.sequence.Sequence;
import org.forester.tools.ConfidenceAssessor;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class MsaCompactor {

    final private static NumberFormat NF_3    = new DecimalFormat( "#.###" );
    final private static NumberFormat NF_4    = new DecimalFormat( "#.####" );
    private static final boolean      VERBOSE = false;
    private Msa                       _msa;
    private File                      _out_file_base;
    private String                    _path_to_mafft;
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

    final public void setOutFileBase( final File out_file_base ) {
        _out_file_base = out_file_base;
    }

    final public String writeMsa( final File outfile, final MSA_FORMAT format, final String suffix ) throws IOException {
        final Double gr = MsaMethods.calcGapRatio( _msa );
        final String s = outfile + "_" + _msa.getNumberOfSequences() + "_" + _msa.getLength() + "_"
                + ForesterUtil.roundToInt( gr * 100 );
        writeMsa( s + suffix, format );
        return s;
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

    Phylogeny pi( final String matrix ) {
        final Phylogeny master_phy = inferNJphylogeny( PWD_DISTANCE_METHOD.KIMURA_DISTANCE, _msa, true, matrix );
        final int seed = 15;
        final int n = 100;
        final ResampleableMsa resampleable_msa = new ResampleableMsa( ( BasicMsa ) _msa );
        final int[][] resampled_column_positions = BootstrapResampler.createResampledColumnPositions( _msa.getLength(),
                                                                                                      n,
                                                                                                      seed );
        final Phylogeny[] eval_phys = new Phylogeny[ n ];
        for( int i = 0; i < n; ++i ) {
            resampleable_msa.resample( resampled_column_positions[ i ] );
            eval_phys[ i ] = inferNJphylogeny( PWD_DISTANCE_METHOD.KIMURA_DISTANCE, resampleable_msa, false, null );
        }
        ConfidenceAssessor.evaluate( "bootstrap", eval_phys, master_phy, true, 1 );
        PhylogenyMethods.extractFastaInformation( master_phy );
        return master_phy;
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
                stats[ row ].divideValue( calcNonGapResidues( _msa.getSequence( row ) ) );
            }
            else {
                stats[ row ].divideValue( _msa.getLength() );
            }
        }
        return stats;
    }

    final private GapContribution[] calcGapContribtionsStats( final boolean norm ) {
        final GapContribution stats[] = calcGapContribtions( norm );
        Arrays.sort( stats );
        // for( final GapContribution stat : stats ) {
        //  final StringBuilder sb = new StringBuilder();
        //  sb.append( stat.getId() );
        //  sb.append( "\t" );
        //  sb.append( NF_4.format( stat.getValue() ) );
        //  sb.append( "\t" );
        //            sb.append( NF_4.format( stat.median() ) );
        //            sb.append( "\t" );
        //            sb.append( NF_4.format( stat.getMin() ) );
        //            sb.append( "\t" );
        //            sb.append( NF_4.format( stat.getMax() ) );
        //sb.append( "\t" );
        //System.out.println( sb );
        // }
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

    private Phylogeny inferNJphylogeny( final PWD_DISTANCE_METHOD pwd_distance_method,
                                        final Msa msa,
                                        final boolean write_matrix,
                                        final String matrix_name ) {
        BasicSymmetricalDistanceMatrix m = null;
        switch ( pwd_distance_method ) {
            case KIMURA_DISTANCE:
                m = PairwiseDistanceCalculator.calcKimuraDistances( msa );
                break;
            case POISSON_DISTANCE:
                m = PairwiseDistanceCalculator.calcPoissonDistances( msa );
                break;
            case FRACTIONAL_DISSIMILARITY:
                m = PairwiseDistanceCalculator.calcFractionalDissimilarities( msa );
                break;
            default:
                throw new IllegalArgumentException( "invalid pwd method" );
        }
        if ( write_matrix ) {
            try {
                m.write( ForesterUtil.createBufferedWriter( matrix_name ) );
            }
            catch ( final IOException e ) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        final NeighborJoiningF nj = NeighborJoiningF.createInstance( false, 5 );
        final Phylogeny phy = nj.execute( m );
        return phy;
    }

    private StringBuilder msaStatsAsSB() {
        final StringBuilder sb = new StringBuilder();
        sb.append( _msa.getNumberOfSequences() );
        sb.append( "\t" );
        sb.append( _msa.getLength() );
        sb.append( "\t" );
        sb.append( NF_3.format( MsaMethods.calcGapRatio( _msa ) ) );
        sb.append( "\t" );
        sb.append( NF_3.format( calculateIdentityRatio( 0, _msa.getLength() - 1, _msa ).arithmeticMean() ) );
        return sb;
    }

    final private void realignWithMafft() throws IOException, InterruptedException {
        //  final MsaInferrer mafft = Mafft
        //       .createInstance( "/home/czmasek/SOFTWARE/MSA/MAFFT/mafft-7.130-without-extensions/scripts/mafft" );
        final MsaInferrer mafft = Mafft.createInstance( _path_to_mafft );
        final List<String> opts = new ArrayList<String>();
        opts.add( "--maxiterate" );
        opts.add( "1000" );
        opts.add( "--localpair" );
        opts.add( "--quiet" );
        _msa = mafft.infer( _msa.asSequenceList(), opts );
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
            removeWorstOffenders( step, 1, false, false, false );
            if ( realign ) {
                realignWithMafft();
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

    final private void removeWorstOffenders( final int to_remove,
                                             final int step,
                                             final boolean realign,
                                             final boolean norm,
                                             final boolean verbose ) throws IOException, InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        for( int j = 0; j < to_remove; ++j ) {
            to_remove_ids.add( stats[ j ].getId() );
            _removed_seq_ids.add( stats[ j ].getId() );
        }
        for( int i = 0; i < to_remove_ids.size(); ++i ) {
            final String id = to_remove_ids.get( i );
            _msa = MsaMethods.removeSequence( _msa, id );
            removeGapColumns();
            if ( verbose ) {
                System.out.print( ForesterUtil.pad( id, 20, ' ', false ) );
                System.out.print( "\t" );
                final StringBuilder sb = msaStatsAsSB();
                System.out.print( sb );
                System.out.print( "\t" );
            }
            if ( ( ( ( i + 1 ) % step ) == 0 ) || ( i == ( to_remove_ids.size() - 1 ) ) ) {
                if ( realign ) {
                    realignWithMafft();
                }
                final String s = writeOutfile();
                if ( verbose ) {
                    System.out.print( "-> " + s );
                }
            }
            if ( verbose ) {
                System.out.println();
            }
        }
    }

    final private void removeViaLength( final int length,
                                        final int step,
                                        final boolean realign,
                                        final boolean norm,
                                        final boolean verbose ) throws IOException, InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        for( final GapContribution gap_gontribution : stats ) {
            to_remove_ids.add( gap_gontribution.getId() );
        }
        int i = 0;
        while ( _msa.getLength() > length ) {
            final String id = to_remove_ids.get( i );
            _msa = MsaMethods.removeSequence( _msa, id );
            removeGapColumns();
            if ( verbose ) {
                System.out.print( ForesterUtil.pad( id, 20, ' ', false ) );
                System.out.print( "\t" );
                final StringBuilder sb = msaStatsAsSB();
                System.out.print( sb );
                System.out.print( "\t" );
            }
            if ( ( ( ( i + 1 ) % step ) == 0 ) || ( _msa.getLength() <= length ) ) {
                if ( realign ) {
                    realignWithMafft();
                }
                final String s = writeOutfile();
                if ( verbose ) {
                    System.out.print( "-> " + s );
                }
            }
            if ( verbose ) {
                System.out.println();
            }
            ++i;
        }
    }

    private void setPathToMafft( final String path_to_mafft ) {
        _path_to_mafft = path_to_mafft;
    }

    final private void writeMsa( final String outfile, final MSA_FORMAT format ) throws IOException {
        final Writer w = ForesterUtil.createBufferedWriter( outfile );
        _msa.write( w, format );
        w.close();
    }

    private String writeOutfile() throws IOException {
        final String s = writeMsa( _out_file_base, MSA_FORMAT.PHYLIP, ".aln" );
        //writeMsa( _out_file_base, MSA_FORMAT.FASTA, ".fasta" );
        return s;
    }

    // Returns null if not path found.
    final public static String guessPathToMafft() {
        String path;
        if ( ForesterUtil.OS_NAME.toLowerCase().indexOf( "win" ) >= 0 ) {
            path = "C:\\Program Files\\mafft-win\\mafft.bat";
            if ( MsaInferrer.isInstalled( path ) ) {
                return path;
            }
        }
        path = "/usr/local/bin/mafft";
        if ( MsaInferrer.isInstalled( path ) ) {
            return path;
        }
        path = "/usr/bin/mafft";
        if ( MsaInferrer.isInstalled( path ) ) {
            return path;
        }
        path = "/bin/mafft";
        if ( MsaInferrer.isInstalled( path ) ) {
            return path;
        }
        path = "mafft";
        if ( MsaInferrer.isInstalled( path ) ) {
            return path;
        }
        return null;
    }

    public final static MsaCompactor reduceGapAverage( final Msa msa,
                                                       final double max_gap_average,
                                                       final int step,
                                                       final boolean realign,
                                                       final int minimal_effective_length,
                                                       final String path_to_mafft,
                                                       final File out ) throws IOException, InterruptedException {
        final MsaCompactor mc = new MsaCompactor( msa );
        if ( realign ) {
            mc.setPathToMafft( path_to_mafft );
        }
        mc.setOutFileBase( out );
        mc.removeViaGapAverage( max_gap_average, step, realign, out, minimal_effective_length );
        return mc;
    }

    public final static MsaCompactor reduceLength( final Msa msa,
                                                   final int length,
                                                   final int step,
                                                   final boolean realign,
                                                   final boolean norm,
                                                   final String path_to_mafft,
                                                   final File out ) throws IOException, InterruptedException {
        final MsaCompactor mc = new MsaCompactor( msa );
        if ( realign ) {
            mc.setPathToMafft( path_to_mafft );
        }
        mc.setOutFileBase( out );
        mc.removeViaLength( length, step, realign, norm, true );
        return mc;
    }

    public final static MsaCompactor removeWorstOffenders( final Msa msa,
                                                           final int worst_offenders_to_remove,
                                                           final int step,
                                                           final boolean realign,
                                                           final boolean norm,
                                                           final String path_to_mafft,
                                                           final File out ) throws IOException, InterruptedException {
        final MsaCompactor mc = new MsaCompactor( msa );
        if ( realign ) {
            mc.setPathToMafft( path_to_mafft );
        }
        mc.setOutFileBase( out );
        mc.removeWorstOffenders( worst_offenders_to_remove, step, realign, norm, true );
        return mc;
    }

    private static DescriptiveStatistics calculateIdentityRatio( final int from, final int to, final Msa msa ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int c = from; c <= to; ++c ) {
            stats.addValue( MsaMethods.calculateIdentityRatio( msa, c ) );
        }
        return stats;
    }
}
