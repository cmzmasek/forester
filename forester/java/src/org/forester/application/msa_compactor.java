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

package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.DeleteableMsa;
import org.forester.msa.Msa.MSA_FORMAT;
import org.forester.msa.MsaInferrer;
import org.forester.msa.MsaMethods;
import org.forester.msa_compactor.Chart;
import org.forester.msa_compactor.MsaCompactor;
import org.forester.msa_compactor.MsaProperties;
import org.forester.util.CommandLineArguments;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class msa_compactor {

    final private static NumberFormat NF_1                                   = new DecimalFormat( "#.0" );
    static {
        NF_1.setRoundingMode( RoundingMode.HALF_UP );
    }
    final static private String       HELP_OPTION_1                          = "help";
    final static private String       HELP_OPTION_2                          = "h";
    final static private String       REMOVE_WORST_OFFENDERS_OPTION          = "r";
    final static private String       AV_GAPINESS_OPTION                     = "g";
    final static private String       STEP_OPTION                            = "s";
    final static private String       LENGTH_OPTION                          = "l";
    final static private String       REALIGN_OPTION                         = "a";
    //
    final static private String       STEP_FOR_DIAGNOSTICS_OPTION            = "sd";
    final static private String       MIN_LENGTH_OPTION                      = "ml";
    final static private String       GAP_RATIO_LENGTH_OPTION                = "gr";
    final static private String       REPORT_ALN_MEAN_IDENTITY               = "q";
    final static private String       OUTPUT_FORMAT_PHYLIP_OPTION            = "p";
    final static private String       OUTPUT_REMOVED_SEQS_OPTION             = "ro";
    final static private String       MAFFT_OPTIONS                          = "mo";
    final static private String       PERFORM_PHYLOGENETIC_INFERENCE         = "t";
    //        
    final static private String       PATH_TO_MAFFT_OPTION                   = "mafft";
    final static private String       DO_NOT_NORMALIZE_FOR_EFF_LENGTH_OPTION = "nn";
    final static private String       PRG_NAME                               = "msa_compactor";
    final static private String       PRG_DESC                               = "multiple sequence aligment compactor";
    final static private String       PRG_VERSION                            = "0.2";
    final static private String       PRG_DATE                               = "140428";
    final static private String       E_MAIL                                 = "czmasek@sanfordburham.org";
    final static private String       WWW                                    = "https://sites.google.com/site/cmzmasek/home/software/forester";

    public static void main( final String args[] ) {
        try {
            final CommandLineArguments cla = new CommandLineArguments( args );
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 )
                    || ( ( cla.getNumberOfNames() < 1 ) || ( cla.getNumberOfNames() > 2 ) ) ) {
                printHelp();
                System.exit( 0 );
            }
            final File in = cla.getFile( 0 );
            File out = null;
            if ( cla.getNumberOfNames() > 1 ) {
                out = cla.getFile( 1 );
            }
            int worst_remove = -1;
            double av_gap = -1;
            int length = -1;
            int step = 1;
            boolean realign = false;
            boolean norm = true;
            String path_to_mafft = null;
            int step_for_diagnostics = 1;
            int min_length = -1;
            double gap_ratio = -1;
            boolean report_aln_mean_identity = false;
            MSA_FORMAT output_format = MSA_FORMAT.FASTA;
            File removed_seqs_out_base = null;
            String mafft_options = "--auto";
            boolean perform_phylogenetic_inference = false;
            final List<String> allowed_options = new ArrayList<String>();
            allowed_options.add( REMOVE_WORST_OFFENDERS_OPTION );
            allowed_options.add( AV_GAPINESS_OPTION );
            allowed_options.add( LENGTH_OPTION );
            allowed_options.add( REALIGN_OPTION );
            allowed_options.add( DO_NOT_NORMALIZE_FOR_EFF_LENGTH_OPTION );
            allowed_options.add( STEP_OPTION );
            allowed_options.add( PATH_TO_MAFFT_OPTION );
            allowed_options.add( STEP_FOR_DIAGNOSTICS_OPTION );
            allowed_options.add( MIN_LENGTH_OPTION );
            allowed_options.add( GAP_RATIO_LENGTH_OPTION );
            allowed_options.add( REPORT_ALN_MEAN_IDENTITY );
            allowed_options.add( OUTPUT_FORMAT_PHYLIP_OPTION );
            allowed_options.add( OUTPUT_REMOVED_SEQS_OPTION );
            allowed_options.add( MAFFT_OPTIONS );
            allowed_options.add( PERFORM_PHYLOGENETIC_INFERENCE );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            DeleteableMsa msa = null;
            final FileInputStream is = new FileInputStream( in );
            if ( FastaParser.isLikelyFasta( in ) ) {
                msa = DeleteableMsa.createInstance( FastaParser.parseMsa( is ) );
            }
            else {
                msa = DeleteableMsa.createInstance( GeneralMsaParser.parse( is ) );
            }
            final DescriptiveStatistics initial_msa_stats = MsaMethods.calculateEffectiveLengthStatistics( msa );
            final boolean chart_only = ( !cla.isOptionSet( LENGTH_OPTION ) )
                    && ( !cla.isOptionSet( REMOVE_WORST_OFFENDERS_OPTION ) )
                    && ( !cla.isOptionSet( AV_GAPINESS_OPTION ) );
            if ( !chart_only && ( out == null ) ) {
                ForesterUtil.fatalError( PRG_NAME, "outfile file missing" );
            }
            if ( cla.isOptionSet( REMOVE_WORST_OFFENDERS_OPTION ) ) {
                worst_remove = cla.getOptionValueAsInt( REMOVE_WORST_OFFENDERS_OPTION );
                if ( ( worst_remove < 1 ) || ( worst_remove >= msa.getNumberOfSequences() - 1 ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "number of worst offender sequences to remove is out of range: "
                            + worst_remove );
                }
            }
            if ( cla.isOptionSet( AV_GAPINESS_OPTION ) ) {
                if ( cla.isOptionSet( REMOVE_WORST_OFFENDERS_OPTION ) ) {
                    printHelp();
                    System.exit( 0 );
                }
                av_gap = cla.getOptionValueAsDouble( AV_GAPINESS_OPTION );
                if ( ( av_gap < 0 ) || ( av_gap >= 1 ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "target gap-ratio is out of range: " + av_gap );
                }
            }
            if ( cla.isOptionSet( LENGTH_OPTION ) ) {
                if ( cla.isOptionSet( REMOVE_WORST_OFFENDERS_OPTION ) || cla.isOptionSet( AV_GAPINESS_OPTION ) ) {
                    printHelp();
                    System.exit( 0 );
                }
                length = cla.getOptionValueAsInt( LENGTH_OPTION );
                if ( length >= msa.getLength() ) {
                    ForesterUtil.fatalError( PRG_NAME,
                                             "target length is out of range [longer than MSA (" + msa.getLength()
                                                     + ")]: " + length );
                }
                else if ( length < initial_msa_stats.getMin() ) {
                    ForesterUtil.fatalError( PRG_NAME,
                                             "target length is out of range [shorter than the shortest sequence ("
                                                     + initial_msa_stats.getMin() + ") ]: " + length );
                }
            }
            if ( cla.isOptionSet( STEP_OPTION ) ) {
                step = cla.getOptionValueAsInt( STEP_OPTION );
                if ( ( step < 1 )
                        || ( ( step > msa.getNumberOfSequences() ) || ( ( worst_remove > 0 ) && ( step > worst_remove ) ) ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "value for step is out of range: " + step );
                }
            }
            if ( cla.isOptionSet( REALIGN_OPTION ) ) {
                realign = true;
            }
            if ( cla.isOptionSet( PATH_TO_MAFFT_OPTION ) ) {
                if ( !realign ) {
                    ForesterUtil.fatalError( PRG_NAME, "no need to indicate path to MAFFT without realigning" );
                }
                path_to_mafft = cla.getOptionValueAsCleanString( PATH_TO_MAFFT_OPTION );
            }
            if ( cla.isOptionSet( DO_NOT_NORMALIZE_FOR_EFF_LENGTH_OPTION ) ) {
                norm = false;
            }
            if ( cla.isOptionSet( STEP_FOR_DIAGNOSTICS_OPTION ) ) {
                step_for_diagnostics = cla.getOptionValueAsInt( STEP_FOR_DIAGNOSTICS_OPTION );
                if ( ( step_for_diagnostics < 1 )
                        || ( ( step_for_diagnostics > msa.getNumberOfSequences() ) || ( ( worst_remove > 0 ) && ( step_for_diagnostics > worst_remove ) ) ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "value for diagnostic step is out of range: "
                            + step_for_diagnostics );
                }
            }
            if ( cla.isOptionSet( MIN_LENGTH_OPTION ) ) {
                min_length = cla.getOptionValueAsInt( MIN_LENGTH_OPTION );
                if ( ( min_length < 2 ) || ( min_length > initial_msa_stats.getMax() ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "value for minimal sequence length is out of range: "
                            + min_length );
                }
            }
            if ( cla.isOptionSet( GAP_RATIO_LENGTH_OPTION ) ) {
                gap_ratio = cla.getOptionValueAsDouble( GAP_RATIO_LENGTH_OPTION );
                if ( ( gap_ratio < 0 ) || ( gap_ratio > 1 ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "gap ratio is out of range: " + gap_ratio );
                }
            }
            if ( cla.isOptionSet( REPORT_ALN_MEAN_IDENTITY ) ) {
                report_aln_mean_identity = true;
            }
            if ( cla.isOptionSet( OUTPUT_FORMAT_PHYLIP_OPTION ) ) {
                output_format = MSA_FORMAT.PHYLIP;
            }
            if ( cla.isOptionSet( OUTPUT_REMOVED_SEQS_OPTION ) ) {
                final String s = cla.getOptionValueAsCleanString( OUTPUT_REMOVED_SEQS_OPTION );
                removed_seqs_out_base = new File( s );
            }
            if ( realign ) {
                if ( ForesterUtil.isEmpty( path_to_mafft ) ) {
                    path_to_mafft = MsaCompactor.guessPathToMafft();
                }
                checkPathToMafft( path_to_mafft );
                if ( cla.isOptionSet( MAFFT_OPTIONS ) ) {
                    mafft_options = cla.getOptionValueAsCleanString( MAFFT_OPTIONS );
                    if ( ForesterUtil.isEmpty( mafft_options ) || ( mafft_options.length() < 3 ) ) {
                        ForesterUtil.fatalError( PRG_NAME, "illegal or empty MAFFT options: " + mafft_options );
                    }
                }
            }
            else if ( cla.isOptionSet( MAFFT_OPTIONS ) ) {
                ForesterUtil.fatalError( PRG_NAME, "no need to indicate MAFFT options without realigning" );
            }
            if ( cla.isOptionSet( PERFORM_PHYLOGENETIC_INFERENCE ) ) {
                perform_phylogenetic_inference = true;
            }
            if ( chart_only ) {
                if ( ( out != null ) || ( removed_seqs_out_base != null ) ) {
                    ForesterUtil
                            .fatalError( PRG_NAME,
                                         "chart only, no outfile(s) produced, thus no need to indicate output file(s)" );
                }
                if ( !realign && cla.isOptionSet( STEP_OPTION ) ) {
                    ForesterUtil.fatalError( PRG_NAME,
                                             "chart only, no re-aligning, thus no need to use step for output and re-aligning; use -"
                                                     + STEP_FOR_DIAGNOSTICS_OPTION + " instead" );
                }
            }
            ForesterUtil.printProgramInformation( PRG_NAME,
                                                  PRG_DESC,
                                                  PRG_VERSION,
                                                  PRG_DATE,
                                                  E_MAIL,
                                                  WWW,
                                                  ForesterUtil.getForesterLibraryInformation() );
            System.out.println( "Input MSA                            : " + in );
            System.out.println( "  MSA length                         : " + msa.getLength() );
            System.out.println( "  Number of sequences                : " + msa.getNumberOfSequences() );
            System.out.println( "  Median sequence length             : " + NF_1.format( initial_msa_stats.median() ) );
            System.out.println( "  Mean sequence length               : "
                    + NF_1.format( initial_msa_stats.arithmeticMean() ) );
            System.out.println( "  Max sequence length                : " + ( ( int ) initial_msa_stats.getMax() ) );
            System.out.println( "  Min sequence length                : " + ( ( int ) initial_msa_stats.getMin() ) );
            if ( !chart_only ) {
                System.out.println( "Output                               : " + out );
            }
            else {
                System.out.println( "Output                               : n/a" );
            }
            if ( removed_seqs_out_base != null ) {
                System.out.println( "Write removed sequences to           : " + removed_seqs_out_base );
            }
            if ( worst_remove > 0 ) {
                System.out.println( "Number of worst offenders to remove  : " + worst_remove );
            }
            if ( av_gap > 0 ) {
                System.out.println( "Target gap-ratio                     : " + av_gap );
            }
            if ( length > 0 ) {
                System.out.println( "Target MSA length                    : " + length );
            }
            if ( min_length > 1 ) {
                System.out.println( "Minimal effective sequence length    : " + min_length );
            }
            if ( gap_ratio > -1 ) {
                System.out.println( "Maximum allowed gap ratio per column : " + gap_ratio );
            }
            if ( ( out != null ) || ( removed_seqs_out_base != null ) ) {
                System.out.println( "Output format                        : "
                        + ( output_format == MSA_FORMAT.FASTA ? "fasta" : "phylip" ) );
            }
            if ( chart_only && !realign ) {
                System.out.println( "Step for output and re-aligning      : n/a" );
            }
            else {
                if ( chart_only ) {
                    System.out.println( "Step for re-aligning                 : " + step );
                }
                else {
                    System.out.println( "Step for output and re-aligning      : " + step );
                }
            }
            System.out.println( "Step for diagnostics reports         : " + step_for_diagnostics );
            System.out.println( "Calculate mean identity              : " + report_aln_mean_identity );
            if ( !norm ) {
                System.out.println( "Normalize                            : " + norm );
            }
            System.out.println( "Realign with MAFFT                   : " + realign );
            if ( realign ) {
                System.out.println( "MAFFT options                        : " + mafft_options );
            }
            System.out.println( "Simple tree (Kimura distances, NJ)   : " + perform_phylogenetic_inference );
            System.out.println();
            final int initial_number_of_seqs = msa.getNumberOfSequences();
            List<MsaProperties> msa_props = null;
            final MsaCompactor mc = new MsaCompactor( msa );
            mc.setInfileName( in.getName() );
            mc.setNorm( norm );
            mc.setRealign( realign );
            if ( realign ) {
                mc.setPathToMafft( path_to_mafft );
                mc.setMafftOptions( mafft_options );
            }
            mc.setStep( step );
            mc.setStepForDiagnostics( step_for_diagnostics );
            mc.setReportAlnMeanIdentity( report_aln_mean_identity );
            mc.setPeformPhylogenticInference( perform_phylogenetic_inference );
            if ( ( worst_remove > 0 ) || ( av_gap > 0 ) || ( length > 0 ) ) {
                mc.setOutputFormat( output_format );
                mc.setOutFileBase( out );
                if ( removed_seqs_out_base != null ) {
                    mc.setRemovedSeqsOutBase( removed_seqs_out_base );
                }
            }
            if ( min_length > 1 ) {
                mc.removeSequencesByMinimalLength( min_length );
                mc.writeMsa( new File( "removed" ) );
            }
            if ( worst_remove > 0 ) {
                msa_props = mc.removeWorstOffenders( worst_remove );
            }
            else if ( av_gap > 0 ) {
                msa_props = mc.removeViaGapAverage( av_gap );
            }
            else if ( length > 0 ) {
                msa_props = mc.removeViaLength( length );
            }
            else {
                msa_props = mc.chart( step, realign, norm );
            }
            Chart.display( msa_props, initial_number_of_seqs, report_aln_mean_identity, in.getName() );
        }
        catch ( final IllegalArgumentException iae ) {
            //  iae.printStackTrace(); //TODO remove me
            ForesterUtil.fatalError( PRG_NAME, iae.getMessage() );
        }
        catch ( final IOException ioe ) {
            // ioe.printStackTrace(); //TODO remove me
            ForesterUtil.fatalError( PRG_NAME, ioe.getMessage() );
        }
        catch ( final Exception e ) {
            ForesterUtil.unexpectedFatalError( PRG_NAME, e );
        }
    }

    private static void checkPathToMafft( final String path_to_mafft ) {
        if ( !ForesterUtil.isEmpty( path_to_mafft ) && MsaInferrer.isInstalled( path_to_mafft ) ) {
        }
        else {
            if ( ForesterUtil.isEmpty( path_to_mafft ) ) {
                ForesterUtil.fatalError( PRG_NAME, "no MAFFT executable found, use -\"" + PATH_TO_MAFFT_OPTION
                        + "=<path to MAFFT>\" option" );
            }
            else {
                ForesterUtil.fatalError( PRG_NAME, "no MAFFT executable at \"" + path_to_mafft + "\"" );
            }
        }
    }

    private static void printHelp() {
        ForesterUtil.printProgramInformation( PRG_NAME,
                                              PRG_DESC,
                                              PRG_VERSION,
                                              PRG_DATE,
                                              E_MAIL,
                                              WWW,
                                              ForesterUtil.getForesterLibraryInformation() );
        final String path_to_mafft = MsaCompactor.guessPathToMafft();
        String mafft_comment;
        if ( !ForesterUtil.isEmpty( path_to_mafft ) ) {
            mafft_comment = " (using " + path_to_mafft + ")";
        }
        else {
            mafft_comment = " (no path to MAFFT found, use -\"" + PATH_TO_MAFFT_OPTION + "=<path to MAFFT>\" option";
        }
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " <options> <msa input file> <output file base>" );
        System.out.println();
        System.out.println( " options: " );
        System.out.println();
        System.out.println( "   -" + REMOVE_WORST_OFFENDERS_OPTION
                + "=<integer>   number of worst offender sequences to remove" );
        System.out.println( "   -" + LENGTH_OPTION + "=<integer>   target MSA length" );
        System.out.println( "   -" + AV_GAPINESS_OPTION + "=<decimal>   target gap-ratio (0.0-1.0)" );
        System.out.println( "   -" + REALIGN_OPTION + "             to realign using MAFFT" + mafft_comment );
        System.out.println( "   -" + MAFFT_OPTIONS + "=<string>   options for MAFFT (default: --auto)" );
        System.out.println( "   -" + STEP_OPTION + "=<integer>   step for output and re-aligning (default: 1)" );
        System.out.println( "   -" + STEP_FOR_DIAGNOSTICS_OPTION
                + "=<integer>  step for diagnostics reports (default: 1)" );
        System.out
                .println( "   -"
                        + REPORT_ALN_MEAN_IDENTITY
                        + "             to calculate mean MSA column identity (\"MSA quality\") (not recommended for very large alignments)" );
        System.out.println( "   -" + OUTPUT_FORMAT_PHYLIP_OPTION
                + "             to write output alignments in phylip format instead of fasta" );
        System.out.println( "   -" + OUTPUT_REMOVED_SEQS_OPTION + "=<file>     to output the removed sequences" );
        System.out.println( "   -" + MIN_LENGTH_OPTION
                + "=<integer>  minimal effecive sequence length (for deleting of shorter sequences)" );
        System.out.println( "   -" + GAP_RATIO_LENGTH_OPTION
                + "=<decimal>  maximal allowed gap ratio per column (for deleting of columms) (0.0-1.0)" );
        System.out.println( "   -" + PERFORM_PHYLOGENETIC_INFERENCE
                + "             to calculate a simple phylogenetic tree (Kimura distances, NJ)" );
        System.out.println();
        System.out.println();
        System.out.println();
    }
}
