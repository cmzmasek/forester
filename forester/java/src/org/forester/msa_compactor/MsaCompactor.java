// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2014 Christian M. Zmasek
// Copyright (C) 2014 Sanford-Burnham Medical Research Institute
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

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
import org.forester.msa.DeleteableMsa;
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
import org.forester.util.ForesterUtil;

public class MsaCompactor {

    final private static NumberFormat NF_3                      = new DecimalFormat( "#.###" );
    final private static NumberFormat NF_4                      = new DecimalFormat( "#.####" );
    //
    private final String              _maffts_opts              = "--auto";
    private int                       _step                     = 1;
    //
    private boolean                   _realign                  = false;
    private boolean                   _norm                     = true;
    private int                       _step_for_diagnostics     = 1;
    private int                       _min_length               = -1;
    private double                    _gap_ratio                = -1;
    private final boolean             _report_aln_mean_identity = false;
    private MSA_FORMAT                _output_format            = MSA_FORMAT.FASTA;
    private final File                _removed_seqs_out_base    = null;
    //
    private DeleteableMsa             _msa;
    private File                      _out_file_base;
    private String                    _path_to_mafft;
    private final SortedSet<String>   _removed_seq_ids;
    static {
        NF_4.setRoundingMode( RoundingMode.HALF_UP );
        NF_3.setRoundingMode( RoundingMode.HALF_UP );
    }

    public MsaCompactor( final DeleteableMsa msa ) {
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
        final ResampleableMsa resampleable_msa = new ResampleableMsa( _msa );
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

    public final List<MsaProperties> chart( final int step,
                                            final boolean realign,
                                            final boolean norm,
                                            final boolean verbose ) throws IOException, InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        final List<MsaProperties> msa_props = new ArrayList<MsaProperties>();
        for( final GapContribution gap_gontribution : stats ) {
            to_remove_ids.add( gap_gontribution.getId() );
        }
        if ( verbose ) {
            printTableHeader();
        }
        int i = 0;
        final int s = _msa.getNumberOfSequences();
        final int x = ForesterUtil.roundToInt( s / 20.0 );
        while ( _msa.getNumberOfSequences() > x ) {
            final String id = to_remove_ids.get( i );
            //~_msa = MsaMethods.removeSequence( _msa, id );
            _msa.deleteRow( id );
            if ( ( s < 500 ) || ( ( step > 0 ) && ( ( ( i + 1 ) % step ) == 0 ) ) ) {
                removeGapColumns();
                if ( realign && ( ( step > 0 ) && ( ( ( i + 1 ) % step ) == 0 ) ) ) {
                    realignWithMafft();
                    msa_props.add( new MsaProperties( _msa ) );
                    if ( verbose ) {
                        printMsaStats( id );
                    }
                    if ( verbose ) {
                        System.out.print( "(realigned)" );
                    }
                }
                else {
                    msa_props.add( new MsaProperties( _msa ) );
                    if ( verbose ) {
                        printMsaStats( id );
                    }
                }
                if ( verbose ) {
                    System.out.println();
                }
            }
            ++i;
        }
        return msa_props;
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
        sb.append( NF_4.format( MsaMethods.calcGapRatio( _msa ) ) );
        sb.append( "\t" );
        sb.append( NF_4.format( MsaMethods.calculateIdentityRatio( 0, _msa.getLength() - 1, _msa ).arithmeticMean() ) );
        return sb;
    }

    private final void printMsaStats( final String id ) {
        System.out.print( ForesterUtil.pad( id, 20, ' ', false ) );
        System.out.print( "\t" );
        final StringBuilder sb = msaStatsAsSB();
        System.out.print( sb );
        System.out.print( "\t" );
    }

    final private void printMsaStatsWriteOutfileAndRealign( final boolean realign,
                                                            final boolean verbose,
                                                            final String id ) throws IOException, InterruptedException {
        if ( realign ) {
            realignWithMafft();
        }
        if ( verbose ) {
            printMsaStats( id );
        }
        final String s = writeOutfile();
        if ( verbose ) {
            System.out.print( "-> " + s + ( realign ? "\t(realigned)" : "" ) );
        }
    }

    final private void realignWithMafft() throws IOException, InterruptedException {
        //  final MsaInferrer mafft = Mafft
        //       .createInstance( "/home/czmasek/SOFTWARE/MSA/MAFFT/mafft-7.130-without-extensions/scripts/mafft" );
        final MsaInferrer mafft = Mafft.createInstance( _path_to_mafft );
        final List<String> opts = new ArrayList<String>();
        for( final String o : _maffts_opts.split( "\\s" ) ) {
            opts.add( o );
        }
        _msa = DeleteableMsa.createInstance( mafft.infer( _msa.asSequenceList(), opts ) );
    }

    final private void removeGapColumns() {
        _msa.deleteGapOnlyColumns();
    }

    public final void removeViaGapAverage( final double mean_gapiness, final boolean verbose ) throws IOException,
            InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( _norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        for( final GapContribution gap_gontribution : stats ) {
            to_remove_ids.add( gap_gontribution.getId() );
        }
        if ( verbose ) {
            printTableHeader();
        }
        int i = 0;
        while ( MsaMethods.calcGapRatio( _msa ) > mean_gapiness ) {
            final String id = to_remove_ids.get( i );
            _msa.deleteRow( id );
            removeGapColumns();
            if ( ( ( _step > 0 ) && ( ( ( i + 1 ) % _step ) == 0 ) )
                    || ( MsaMethods.calcGapRatio( _msa ) <= mean_gapiness ) ) {
                printMsaStatsWriteOutfileAndRealign( _realign, verbose, id );
            }
            else if ( verbose ) {
                printMsaStats( id );
            }
            if ( verbose ) {
                System.out.println();
            }
            ++i;
        }
    }

    public void removeViaLength( final int length, final boolean verbose ) throws IOException, InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( _norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        for( final GapContribution gap_gontribution : stats ) {
            to_remove_ids.add( gap_gontribution.getId() );
        }
        if ( verbose ) {
            printTableHeader();
        }
        int i = 0;
        while ( _msa.getLength() > length ) {
            final String id = to_remove_ids.get( i );
            _msa.deleteRow( id );
            removeGapColumns();
            if ( ( ( _step > 0 ) && ( ( ( i + 1 ) % _step ) == 0 ) ) || ( _msa.getLength() <= length ) ) {
                printMsaStatsWriteOutfileAndRealign( _realign, verbose, id );
            }
            else if ( verbose ) {
                printMsaStats( id );
            }
            if ( verbose ) {
                System.out.println();
            }
            ++i;
        }
    }

    public final void removeWorstOffenders( final int to_remove, final boolean verbose ) throws IOException,
            InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( _norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        for( int j = 0; j < to_remove; ++j ) {
            to_remove_ids.add( stats[ j ].getId() );
            _removed_seq_ids.add( stats[ j ].getId() );
        }
        if ( verbose ) {
            printTableHeader();
        }
        for( int i = 0; i < to_remove_ids.size(); ++i ) {
            final String id = to_remove_ids.get( i );
            _msa.deleteRow( id );
            removeGapColumns();
            if ( ( ( _step > 0 ) && ( ( ( i + 1 ) % _step ) == 0 ) ) || ( i == ( to_remove_ids.size() - 1 ) ) ) {
                printMsaStatsWriteOutfileAndRealign( _realign, verbose, id );
            }
            else if ( verbose ) {
                printMsaStats( id );
            }
            if ( verbose ) {
                System.out.println();
            }
        }
    }

    public void setPathToMafft( final String path_to_mafft ) {
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
        path = "/home/czmasek/SOFTWARE/MSA/MAFFT/mafft-7.130-without-extensions/scripts/mafft";
        if ( MsaInferrer.isInstalled( path ) ) {
            return path;
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

    public final void setStep( int step ) {
        _step = step;
    }

    public final void setNorm( boolean norm ) {
        _norm = norm;
    }

    public final void setStepForDiagnostics( int step_for_diagnostics ) {
        _step_for_diagnostics = step_for_diagnostics;
    }

    public final void setMinLength( int min_length ) {
        _min_length = min_length;
    }

    public final void setGapRatio( double gap_ratio ) {
        _gap_ratio = gap_ratio;
    }

    public final void setOutputFormat( MSA_FORMAT output_format ) {
        _output_format = output_format;
    }

    public final void setRealign( boolean realign ) {
        _realign = realign;
    }

    private final static void printTableHeader() {
        System.out.print( ForesterUtil.pad( "Id", 20, ' ', false ) );
        System.out.print( "\t" );
        System.out.print( "Seqs" );
        System.out.print( "\t" );
        System.out.print( "Length" );
        System.out.print( "\t" );
        System.out.print( "Gaps" );
        System.out.print( "\t" );
        System.out.print( "MSA qual" );
        System.out.print( "\t" );
        System.out.println();
    }
}
