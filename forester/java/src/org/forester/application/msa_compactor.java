
package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.Msa;
import org.forester.msa.MsaInferrer;
import org.forester.msa_compactor.MsaCompactor;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class msa_compactor {

    final static private String HELP_OPTION_1                          = "help";
    final static private String HELP_OPTION_2                          = "h";
    final static private String REMOVE_WORST_OFFENDERS_OPTION          = "r";
    final static private String AV_GAPINESS_OPTION                     = "g";
    final static private String STEP_OPTION                            = "s";
    final static private String LENGTH_OPTION                          = "l";
    final static private String REALIGN_OPTION                         = "a";
    final static private String PATH_TO_MAFFT_OPTION                   = "mafft";
    final static private String DO_NOT_NORMALIZE_FOR_EFF_LENGTH_OPTION = "nn";
    final static private String PRG_NAME                               = "msa_compactor";
    final static private String PRG_DESC                               = "multiple sequence aligment compactor";
    final static private String PRG_VERSION                            = "0.01";
    final static private String PRG_DATE                               = "140316";
    final static private String E_MAIL                                 = "phylosoft@gmail.com";
    final static private String WWW                                    = "https://sites.google.com/site/cmzmasek/home/software/forester";

    public static void main( final String args[] ) {
        try {
            final CommandLineArguments cla = new CommandLineArguments( args );
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( cla.getNumberOfNames() != 2 ) ) {
                printHelp();
                System.exit( 0 );
            }
            final File in = cla.getFile( 0 );
            final File out = cla.getFile( 1 );
            int worst_remove = -1;
            double av_gap = -1;
            int length = -1;
            int step = -1;
            boolean realign = false;
            boolean norm = true;
            String path_to_mafft = null;
            final List<String> allowed_options = new ArrayList<String>();
            allowed_options.add( REMOVE_WORST_OFFENDERS_OPTION );
            allowed_options.add( AV_GAPINESS_OPTION );
            allowed_options.add( LENGTH_OPTION );
            allowed_options.add( REALIGN_OPTION );
            allowed_options.add( DO_NOT_NORMALIZE_FOR_EFF_LENGTH_OPTION );
            allowed_options.add( STEP_OPTION );
            allowed_options.add( PATH_TO_MAFFT_OPTION );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            Msa msa = null;
            final FileInputStream is = new FileInputStream( in );
            if ( FastaParser.isLikelyFasta( in ) ) {
                msa = FastaParser.parseMsa( is );
            }
            else {
                msa = GeneralMsaParser.parse( is );
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
                if ( ( length < 2 ) || ( length >= msa.getLength() ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "target length is out of range: " + length );
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
            if ( realign ) {
                if ( ForesterUtil.isEmpty( path_to_mafft ) ) {
                    path_to_mafft = MsaCompactor.guessPathToMafft();
                }
                checkPathToMafft( path_to_mafft );
            }
            if ( worst_remove > 0 ) {
                MsaCompactor.removeWorstOffenders( msa, worst_remove, step, realign, norm, path_to_mafft, out );
            }
            else if ( av_gap > 0 ) {
                MsaCompactor.reduceGapAverage( msa, av_gap, step, realign, norm, path_to_mafft, out );
            }
            else if ( length > 0 ) {
                // TODO if < shortest seq -> error
                MsaCompactor.reduceLength( msa, length, step, realign, norm, path_to_mafft, out );
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    private static void checkPathToMafft( final String path_to_mafft ) {
        if ( !ForesterUtil.isEmpty( path_to_mafft ) && MsaInferrer.isInstalled( path_to_mafft ) ) {
            ForesterUtil.programMessage( PRG_NAME, "using MAFFT at \"" + path_to_mafft + "\"" );
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
        System.out.println( PRG_NAME + " <options> <msa input file> <output file>" );
        System.out.println();
        System.out.println( " options: " );
        System.out.println();
        System.out.println( "   -" + REMOVE_WORST_OFFENDERS_OPTION
                + "=<integer>  number of worst offender sequences to remove" );
        System.out.println( "   -" + LENGTH_OPTION + "=<integer>  target MSA length" );
        System.out.println( "   -" + AV_GAPINESS_OPTION + "=<decimal>  target gap-ratio (0.0-1.0)" );
        System.out.println( "   -" + STEP_OPTION + "=<integer>  step (for output and re-aligning)" );
        System.out.println( "   -" + REALIGN_OPTION + "            to realign using MAFFT" + mafft_comment );
        System.out.println();
        System.out.println();
        System.out.println();
    }
}
