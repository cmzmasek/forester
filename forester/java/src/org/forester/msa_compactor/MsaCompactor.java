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
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
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

    final private static NumberFormat NF_3                   = new DecimalFormat( "#.###" );
    final private static NumberFormat NF_4                   = new DecimalFormat( "#.####" );
    private double                    _gap_ratio             = -1;
    //
    private String                    _maffts_opts           = "--auto";
    private int                       _min_length            = -1;
    //
    private DeleteableMsa             _msa                   = null;
    private boolean                   _norm                  = true;
    private File                      _out_file_base         = null;
    private MSA_FORMAT                _output_format         = MSA_FORMAT.FASTA;
    private String                    _path_to_mafft         = null;
    //
    private boolean                   _realign               = false;
    private final SortedSet<String>   _removed_seq_ids;
    private File                      _removed_seqs_out_base = null;

    public final void setRemovedSeqsOutBase( final File removed_seqs_out_base ) {
        _removed_seqs_out_base = removed_seqs_out_base;
    }
    private boolean                   _report_aln_mean_identity = false;
    private int                       _step                     = -1;
    private int                       _step_for_diagnostics     = -1;
    private final short               _longest_id_length;
    private final ArrayList<Sequence> _removed_seqs;
    static {
        NF_4.setRoundingMode( RoundingMode.HALF_UP );
        NF_3.setRoundingMode( RoundingMode.HALF_UP );
    }

    public MsaCompactor( final DeleteableMsa msa ) {
        _msa = msa;
        _removed_seq_ids = new TreeSet<String>();
        _longest_id_length = _msa.determineMaxIdLength();
        _removed_seqs = new ArrayList<Sequence>();
    }

    final public Msa getMsa() {
        return _msa;
    }

    final public SortedSet<String> getRemovedSeqIds() {
        return _removed_seq_ids;
    }

    public final List<MsaProperties> removeViaGapAverage( final double mean_gapiness ) throws IOException,
            InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( _norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        final List<MsaProperties> msa_props = new ArrayList<MsaProperties>();
        for( final GapContribution gap_gontribution : stats ) {
            to_remove_ids.add( gap_gontribution.getId() );
        }
        printTableHeader();
        MsaProperties msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
        msa_props.add( msa_prop );
        printMsaProperties( "", msa_prop );
        System.out.println();
        int i = 0;
        while ( MsaMethods.calcGapRatio( _msa ) > mean_gapiness ) {
            final String id = to_remove_ids.get( i );
            _removed_seq_ids.add( id );
            final Sequence deleted = _msa.deleteRow( id, true );
            _removed_seqs.add( deleted );
            removeGapColumns();
            if ( isPrintMsaStatsWriteOutfileAndRealign( i ) || ( MsaMethods.calcGapRatio( _msa ) <= mean_gapiness ) ) {
                msa_prop = printMsaStatsWriteOutfileAndRealign( _realign, id );
                msa_props.add( msa_prop );
                System.out.println();
            }
            else if ( isPrintMsaStats( i ) ) {
                msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
                msa_props.add( msa_prop );
                printMsaProperties( id, msa_prop );
                System.out.println();
            }
            ++i;
        }
        if ( _removed_seqs_out_base != null ) {
            final String msg = writeAndAlignRemovedSeqs();
            System.out.println();
            System.out.println( msg );
        }
        return msa_props;
    }

    public List<MsaProperties> removeViaLength( final int length ) throws IOException, InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( _norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        final List<MsaProperties> msa_props = new ArrayList<MsaProperties>();
        for( final GapContribution gap_gontribution : stats ) {
            to_remove_ids.add( gap_gontribution.getId() );
        }
        printTableHeader();
        MsaProperties msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
        msa_props.add( msa_prop );
        printMsaProperties( "", msa_prop );
        System.out.println();
        int i = 0;
        while ( _msa.getLength() > length ) {
            final String id = to_remove_ids.get( i );
            _removed_seq_ids.add( id );
            final Sequence deleted = _msa.deleteRow( id, true );
            _removed_seqs.add( deleted );
            removeGapColumns();
            if ( isPrintMsaStatsWriteOutfileAndRealign( i ) || ( _msa.getLength() <= length ) ) {
                msa_prop = printMsaStatsWriteOutfileAndRealign( _realign, id );
                msa_props.add( msa_prop );
                System.out.println();
            }
            else if ( isPrintMsaStats( i ) ) {
                msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
                printMsaProperties( id, msa_prop );
                msa_props.add( msa_prop );
                System.out.println();
            }
            ++i;
        }
        if ( _removed_seqs_out_base != null ) {
            final String msg = writeAndAlignRemovedSeqs();
            System.out.println();
            System.out.println( msg );
        }
        return msa_props;
    }

    public final List<MsaProperties> removeWorstOffenders( final int to_remove ) throws IOException,
            InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( _norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        final List<MsaProperties> msa_props = new ArrayList<MsaProperties>();
        for( int j = 0; j < to_remove; ++j ) {
            to_remove_ids.add( stats[ j ].getId() );
            _removed_seq_ids.add( stats[ j ].getId() );
        }
        printTableHeader();
        MsaProperties msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
        msa_props.add( msa_prop );
        printMsaProperties( "", msa_prop );
        System.out.println();
        for( int i = 0; i < to_remove_ids.size(); ++i ) {
            final String id = to_remove_ids.get( i );
            _removed_seq_ids.add( id );
            final Sequence deleted = _msa.deleteRow( id, true );
            _removed_seqs.add( deleted );
            removeGapColumns();
            if ( isPrintMsaStatsWriteOutfileAndRealign( i ) || ( i == ( to_remove_ids.size() - 1 ) ) ) {
                msa_prop = printMsaStatsWriteOutfileAndRealign( _realign, id );
                msa_props.add( msa_prop );
                System.out.println();
            }
            else if ( isPrintMsaStats( i ) ) {
                msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
                msa_props.add( msa_prop );
                printMsaProperties( id, msa_prop );
                System.out.println();
            }
        }
        if ( _removed_seqs_out_base != null ) {
            final String msg = writeAndAlignRemovedSeqs();
            System.out.println();
            System.out.println( msg );
        }
        return msa_props;
    }

    public final List<MsaProperties> chart( final int step, final boolean realign, final boolean norm )
            throws IOException, InterruptedException {
        final GapContribution stats[] = calcGapContribtionsStats( norm );
        final List<String> to_remove_ids = new ArrayList<String>();
        final List<MsaProperties> msa_props = new ArrayList<MsaProperties>();
        for( final GapContribution gap_gontribution : stats ) {
            to_remove_ids.add( gap_gontribution.getId() );
        }
        printTableHeader();
        final int x = ForesterUtil.roundToInt( _msa.getNumberOfSequences() / 20.0 );
        MsaProperties msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
        msa_props.add( msa_prop );
        printMsaProperties( "", msa_prop );
        System.out.println();
        int i = 0;
        while ( _msa.getNumberOfSequences() > x ) {
            final String id = to_remove_ids.get( i );
            _msa.deleteRow( id, false );
            if ( realign && isPrintMsaStatsWriteOutfileAndRealign( i ) ) {
                removeGapColumns();
                realignWithMafft();
                msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
                msa_props.add( msa_prop );
                printMsaProperties( id, msa_prop );
                System.out.print( "(realigned)" );
                System.out.println();
            }
            else if ( isPrintMsaStats( i ) ) {
                removeGapColumns();
                msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
                msa_props.add( msa_prop );
                printMsaProperties( id, msa_prop );
                System.out.println();
            }
            ++i;
        }
        return msa_props;
    }

    private final boolean isPrintMsaStats( final int i ) {
        return ( ( ( _step < 2 ) && ( _step_for_diagnostics < 2 ) ) || ( ( _step_for_diagnostics > 0 ) && ( ( ( i + 1 ) % _step_for_diagnostics ) == 0 ) ) );
    }

    private final boolean isPrintMsaStatsWriteOutfileAndRealign( final int i ) {
        return ( ( ( _step < 2 ) && ( _step_for_diagnostics < 2 ) ) || ( ( _step > 0 ) && ( ( ( i + 1 ) % _step ) == 0 ) ) );
    }

    public final void setGapRatio( final double gap_ratio ) {
        _gap_ratio = gap_ratio;
    }

    public final void setMinLength( final int min_length ) {
        _min_length = min_length;
    }

    public final void setNorm( final boolean norm ) {
        _norm = norm;
    }

    final public void setOutFileBase( final File out_file_base ) {
        _out_file_base = out_file_base;
    }

    public final void setOutputFormat( final MSA_FORMAT output_format ) {
        _output_format = output_format;
    }

    public void setPathToMafft( final String path_to_mafft ) {
        _path_to_mafft = path_to_mafft;
    }

    public final void setRealign( final boolean realign ) {
        _realign = realign;
    }

    public final void setStep( final int step ) {
        _step = step;
    }

    public final void setStepForDiagnostics( final int step_for_diagnostics ) {
        _step_for_diagnostics = step_for_diagnostics;
    }

    public final void setReportAlnMeanIdentity( final boolean report_aln_mean_identity ) {
        _report_aln_mean_identity = report_aln_mean_identity;
    }

    final public String writeMsa( final File outfile ) throws IOException {
        final Double gr = MsaMethods.calcGapRatio( _msa );
        final String s = outfile + "_" + _msa.getNumberOfSequences() + "_" + _msa.getLength() + "_"
                + ForesterUtil.roundToInt( gr * 100 );
        writeMsa( _msa, s + obtainSuffix(), _output_format );
        return s;
    }

    final public String writeAndAlignRemovedSeqs() throws IOException, InterruptedException {
        final StringBuilder msg = new StringBuilder();
        final String n = _removed_seqs_out_base + "_" + _removed_seqs.size() + ".fasta";
        SequenceWriter.writeSeqs( _removed_seqs, new File( n ), SEQ_FORMAT.FASTA, 100 );
        msg.append( "wrote " + _removed_seqs.size() + " removed sequences to " + "\"" + n + "\"" );
        if ( _realign ) {
            final MsaInferrer mafft = Mafft.createInstance( _path_to_mafft );
            final List<String> opts = new ArrayList<String>();
            for( final String o : _maffts_opts.split( "\\s" ) ) {
                opts.add( o );
            }
            final Msa removed_msa = mafft.infer( _removed_seqs, opts );
            final Double gr = MsaMethods.calcGapRatio( removed_msa );
            String s = _removed_seqs_out_base + "_" + removed_msa.getNumberOfSequences() + "_"
                    + removed_msa.getLength() + "_" + ForesterUtil.roundToInt( gr * 100 );
            final String suffix = obtainSuffix();
            s += suffix;
            writeMsa( removed_msa, s, _output_format );
            msg.append( ", and as MSA of length " + removed_msa.getLength() + " to \"" + s + "\"" );
        }
        return msg.toString();
    }

    private String obtainSuffix() {
        if ( _output_format == MSA_FORMAT.FASTA ) {
            return ".fasta";
        }
        else if ( _output_format == MSA_FORMAT.PHYLIP ) {
            return ".aln";
        }
        return "";
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

    private final Phylogeny inferNJphylogeny( final PWD_DISTANCE_METHOD pwd_distance_method,
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

    private final Phylogeny pi( final String matrix ) {
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

    private final void printMsaProperties( final String id, final MsaProperties msa_properties ) {
        if ( ( _step < 2 ) || ( _step_for_diagnostics < 2 ) ) {
            System.out.print( ForesterUtil.pad( id, _longest_id_length, ' ', false ) );
            System.out.print( "\t" );
        }
        System.out.print( msaPropertiesAsSB( msa_properties ) );
        System.out.print( "\t" );
    }

    private final StringBuilder msaPropertiesAsSB( final MsaProperties msa_properties ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( msa_properties.getNumberOfSequences() );
        sb.append( "\t" );
        sb.append( msa_properties.getLength() );
        sb.append( "\t" );
        sb.append( NF_4.format( msa_properties.getGapRatio() ) );
        if ( _report_aln_mean_identity /*msa_properties.getAverageIdentityRatio() >= 0*/) {
            sb.append( "\t" );
            sb.append( NF_4.format( msa_properties.getAverageIdentityRatio() ) );
        }
        return sb;
    }

    final private MsaProperties printMsaStatsWriteOutfileAndRealign( final boolean realign, final String id )
            throws IOException, InterruptedException {
        if ( realign ) {
            realignWithMafft();
        }
        final MsaProperties msa_prop = new MsaProperties( _msa, _report_aln_mean_identity );
        printMsaProperties( id, msa_prop );
        final String s = writeOutfile();
        System.out.print( "-> " + s + ( realign ? "\t(realigned)" : "" ) );
        return msa_prop;
    }

    final private void realignWithMafft() throws IOException, InterruptedException {
        final MsaInferrer mafft = Mafft.createInstance( _path_to_mafft );
        final List<String> opts = new ArrayList<String>();
        for( final String o : _maffts_opts.split( "\\s" ) ) {
            opts.add( o );
        }
        _msa = DeleteableMsa.createInstance( mafft.infer( _msa.asSequenceList(), opts ) );
    }

    public final void setMafftOptions( final String maffts_opts ) {
        _maffts_opts = maffts_opts;
    }

    final private void removeGapColumns() {
        _msa.deleteGapOnlyColumns();
    }

    final private static void writeMsa( final Msa msa, final String outfile, final MSA_FORMAT format )
            throws IOException {
        final Writer w = ForesterUtil.createBufferedWriter( outfile );
        msa.write( w, format );
        w.close();
    }

    private final String writeOutfile() throws IOException {
        final String s = writeMsa( _out_file_base );
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

    private final void printTableHeader() {
        if ( ( _step < 2 ) || ( _step_for_diagnostics < 2 ) ) {
            System.out.print( ForesterUtil.pad( "Id", _longest_id_length, ' ', false ) );
            System.out.print( "\t" );
        }
        System.out.print( "Seqs" );
        System.out.print( "\t" );
        System.out.print( "Length" );
        System.out.print( "\t" );
        System.out.print( "Gaps" );
        System.out.print( "\t" );
        if ( _report_aln_mean_identity ) {
            System.out.print( "MSA qual" );
            System.out.print( "\t" );
        }
        System.out.println();
    }
}
