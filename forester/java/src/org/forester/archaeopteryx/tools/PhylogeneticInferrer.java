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

package org.forester.archaeopteryx.tools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JOptionPane;

import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.evoinference.distance.NeighborJoiningF;
import org.forester.evoinference.distance.PairwiseDistanceCalculator;
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
import org.forester.util.ForesterUtil;

public class PhylogeneticInferrer extends RunnableProcess {

    private Msa                                _msa;
    private final MainFrameApplication         _mf;
    private final PhylogeneticInferenceOptions _options;
    private final List<Sequence>               _seqs;
    private final boolean                      DEBUG           = true;
    public final static String                 MSA_FILE_SUFFIX = ".aln";
    public final static String                 PWD_FILE_SUFFIX = ".pwd";

    public PhylogeneticInferrer( final List<Sequence> seqs,
                                 final PhylogeneticInferenceOptions options,
                                 final MainFrameApplication mf ) {
        _msa = null;
        _seqs = seqs;
        _mf = mf;
        _options = options;
    }

    public PhylogeneticInferrer( final Msa msa,
                                 final PhylogeneticInferenceOptions options,
                                 final MainFrameApplication mf ) {
        _msa = msa;
        _seqs = null;
        _mf = mf;
        _options = options;
    }

    private Msa inferMsa( final MSA_PRG msa_prg ) throws IOException, InterruptedException {
        //        final File temp_seqs_file = File.createTempFile( "__msa__temp__", ".fasta" );
        //        if ( DEBUG ) {
        //            System.out.println();
        //            System.out.println( "temp file: " + temp_seqs_file );
        //            System.out.println();
        //        }
        //        //final File temp_seqs_file = new File( _options.getTempDir() + ForesterUtil.FILE_SEPARATOR + "s.fasta" );
        //        final BufferedWriter writer = new BufferedWriter( new FileWriter( temp_seqs_file ) );
        //        SequenceWriter.writeSeqs( _seqs, writer, SEQ_FORMAT.FASTA, 100 );
        //        writer.close();
        switch ( msa_prg ) {
            case MAFFT:
                return runMAFFT( _seqs, processMafftOptions() );
            default:
                return null;
        }
    }

    private List<String> processMafftOptions() {
        final String opts_str = _options.getMsaPrgParameters().trim().toLowerCase();
        final String[] opts_ary = opts_str.split( " " );
        final List<String> opts = new ArrayList<String>();
        boolean saw_quiet = false;
        for( final String opt : opts_ary ) {
            opts.add( opt );
            if ( opt.equals( "--quiet" ) ) {
                saw_quiet = true;
            }
        }
        if ( !saw_quiet ) {
            opts.add( "--quiet" );
        }
        return opts;
    }

    private Phylogeny inferPhylogeny( final Msa msa ) {
        BasicSymmetricalDistanceMatrix m = null;
        switch ( _options.getPwdDistanceMethod() ) {
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
                throw new RuntimeException( "invalid pwd method" );
        }
        if ( !ForesterUtil.isEmpty( _options.getIntermediateFilesBase() ) ) {
            BufferedWriter pwd_writer;
            try {
                pwd_writer = new BufferedWriter( new FileWriter( _options.getIntermediateFilesBase() + PWD_FILE_SUFFIX ) );
                m.write( pwd_writer );
                pwd_writer.close();
            }
            catch ( final IOException e ) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        final NeighborJoiningF nj = NeighborJoiningF.createInstance( false, 5 );
        final Phylogeny phy = nj.execute( m );
        PhylogenyMethods.extractFastaInformation( phy );
        return phy;
    }

    private void infer() throws InterruptedException {
        //_mf.getMainPanel().getCurrentTreePanel().setWaitCursor();
        if ( ( _msa == null ) && ( _seqs == null ) ) {
            throw new IllegalArgumentException( "cannot run phylogenetic analysis with null msa and seq array" );
        }
        start( _mf, "phylogenetic inference" );
        if ( _msa == null ) {
            Msa msa = null;
            try {
                msa = inferMsa( MSA_PRG.MAFFT );
            }
            catch ( final IOException e ) {
                end( _mf );
                JOptionPane.showMessageDialog( _mf,
                                               "Could not create multiple sequence alignment with \""
                                                       + _options.getMsaPrg() + "\" and the following parameters:\n\""
                                                       + _options.getMsaPrgParameters() + "\"\nError: "
                                                       + e.getLocalizedMessage(),
                                               "Failed to Calculate MSA",
                                               JOptionPane.ERROR_MESSAGE );
                if ( DEBUG ) {
                    e.printStackTrace();
                }
                return;
            }
            catch ( final Exception e ) {
                end( _mf );
                JOptionPane.showMessageDialog( _mf,
                                               "Could not create multiple sequence alignment with \""
                                                       + _options.getMsaPrg() + "\" and the following parameters:\n\""
                                                       + _options.getMsaPrgParameters() + "\"\nError: "
                                                       + e.getLocalizedMessage(),
                                               "Unexpected Exception During MSA Calculation",
                                               JOptionPane.ERROR_MESSAGE );
                if ( DEBUG ) {
                    e.printStackTrace();
                }
                return;
            }
            if ( msa == null ) {
                end( _mf );
                JOptionPane.showMessageDialog( _mf,
                                               "Could not create multiple sequence alignment with "
                                                       + _options.getMsaPrg() + "\nand the following parameters:\n\""
                                                       + _options.getMsaPrgParameters() + "\"",
                                               "Failed to Calculate MSA",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            if ( DEBUG ) {
                System.out.println( msa.toString() );
                System.out.println( MsaMethods.calcGapRatio( msa ) );
            }
            final MsaMethods msa_tools = MsaMethods.createInstance();
            if ( _options.isExecuteMsaProcessing() ) {
                msa = msa_tools.removeGapColumns( _options.getMsaProcessingMaxAllowedGapRatio(),
                                                  _options.getMsaProcessingMinAllowedLength(),
                                                  msa );
                if ( msa == null ) {
                    end( _mf );
                    JOptionPane.showMessageDialog( _mf,
                                                   "Less than two sequences longer than "
                                                           + _options.getMsaProcessingMinAllowedLength()
                                                           + " residues left after MSA processing",
                                                   "MSA Processing Settings Too Stringent",
                                                   JOptionPane.ERROR_MESSAGE );
                    return;
                }
            }
            if ( DEBUG ) {
                System.out.println( msa_tools.getIgnoredSequenceIds() );
                System.out.println( msa.toString() );
                System.out.println( MsaMethods.calcGapRatio( msa ) );
            }
            _msa = msa;
        }
        final int n = _options.getBootstrapSamples();
        final long seed = _options.getRandomNumberGeneratorSeed();
        final Phylogeny master_phy = inferPhylogeny( _msa );
        if ( _options.isPerformBootstrapResampling() && ( n > 0 ) ) {
            final ResampleableMsa resampleable_msa = new ResampleableMsa( ( BasicMsa ) _msa );
            final int[][] resampled_column_positions = BootstrapResampler.createResampledColumnPositions( _msa
                    .getLength(), n, seed );
            final Phylogeny[] eval_phys = new Phylogeny[ n ];
            for( int i = 0; i < n; ++i ) {
                resampleable_msa.resample( resampled_column_positions[ i ] );
                eval_phys[ i ] = inferPhylogeny( resampleable_msa );
            }
            ConfidenceAssessor.evaluate( "bootstrap", eval_phys, master_phy, true, 1 );
        }
        _mf.getMainPanel().addPhylogenyInNewTab( master_phy, _mf.getConfiguration(), "nj", "njpath" );
        //  _mf.getMainPanel().getCurrentTreePanel().setArrowCursor();
        end( _mf );
        JOptionPane.showMessageDialog( _mf,
                                       "Inference successfully completed",
                                       "Inference Completed",
                                       JOptionPane.INFORMATION_MESSAGE );
    }

    @Override
    public void run() {
        try {
            infer();
        }
        catch ( final InterruptedException e ) {
            // TODO need to handle this exception SOMEHOW!
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private Msa runMAFFT( final List<Sequence> seqs, final List<String> opts ) throws IOException, InterruptedException {
        Msa msa = null;
        final MsaInferrer mafft = Mafft.createInstance( _mf.getInferenceManager().getPathToLocalMafft()
                .getCanonicalPath() );
        try {
            msa = mafft.infer( seqs, opts );
        }
        catch ( final IOException e ) {
            System.out.println( mafft.getErrorDescription() );
        }
        return msa;
    }

    private void writeToFiles( final BasicSymmetricalDistanceMatrix m ) {
        if ( !ForesterUtil.isEmpty( _options.getIntermediateFilesBase() ) ) {
            try {
                final BufferedWriter msa_writer = new BufferedWriter( new FileWriter( _options.getIntermediateFilesBase()
                        + MSA_FILE_SUFFIX ) );
                _msa.write( msa_writer, MSA_FORMAT.PHYLIP );
                msa_writer.close();
                final BufferedWriter pwd_writer = new BufferedWriter( new FileWriter( _options.getIntermediateFilesBase()
                        + PWD_FILE_SUFFIX ) );
                m.write( pwd_writer );
                pwd_writer.close();
            }
            catch ( final Exception e ) {
                System.out.println( "Error: " + e.getMessage() );
            }
        }
    }

    public enum MSA_PRG {
        MAFFT;
    }
}
