
package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.Msa;
import org.forester.msa.Msa.MSA_FORMAT;
import org.forester.msa_compactor.MsaCompactor;
import org.forester.msa.MsaMethods;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class msa_compactor {

    final static private String HELP_OPTION_1                          = "help";
    final static private String HELP_OPTION_2                          = "h";
    final static private String REMOVE_WORST_OFFENDERS_OPTION          = "w";
    final static private String AV_GAPINESS_OPTION                     = "a";
    final static private String STEP_OPTION                            = "s";
    final static private String LENGTH_OPTION                          = "l";
    final static private String REALIGN_OPTION                         = "r";
    final static private String DO_NOT_NORMALIZE_FOR_EFF_LENGTH_OPTION = "nn";
    final static private String PRG_NAME                               = "msa_compactor";
    final static private String PRG_DESC                               = "multiple sequnce aligment compactor";
    final static private String PRG_VERSION                            = "0.01";
    final static private String PRG_DATE                               = "140221";
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
            double av = -1;
            int length = -1;
            int step = 1;
            boolean realign = false;
            boolean norm = true;
            final List<String> allowed_options = new ArrayList<String>();
            allowed_options.add( REMOVE_WORST_OFFENDERS_OPTION );
            allowed_options.add( AV_GAPINESS_OPTION );
            allowed_options.add( LENGTH_OPTION );
            allowed_options.add( REALIGN_OPTION );
            allowed_options.add( DO_NOT_NORMALIZE_FOR_EFF_LENGTH_OPTION );
            allowed_options.add( STEP_OPTION );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            if ( cla.isOptionSet( REMOVE_WORST_OFFENDERS_OPTION ) ) {
                worst_remove = cla.getOptionValueAsInt( REMOVE_WORST_OFFENDERS_OPTION );
            }
            if ( cla.isOptionSet( AV_GAPINESS_OPTION ) ) {
                av = cla.getOptionValueAsDouble( AV_GAPINESS_OPTION );
            }
            if ( cla.isOptionSet( LENGTH_OPTION ) ) {
                length = cla.getOptionValueAsInt( LENGTH_OPTION );
            }
            if ( cla.isOptionSet( STEP_OPTION ) ) {
                step = cla.getOptionValueAsInt( STEP_OPTION );
            }
            if ( cla.isOptionSet( REALIGN_OPTION ) ) {
                realign = true;
            }
            if ( cla.isOptionSet( DO_NOT_NORMALIZE_FOR_EFF_LENGTH_OPTION ) ) {
                norm = false;
            }
            //            else if ( cla.isOptionSet( STEP_OPTION ) && cla.isOptionSet( WINDOW_OPTION ) ) {
            //                step = cla.getOptionValueAsInt( STEP_OPTION );
            //                window = cla.getOptionValueAsInt( WINDOW_OPTION );
            //            }
            //            else {
            //                printHelp();
            //                System.exit( 0 );
            //            }
            Msa msa = null;
            final FileInputStream is = new FileInputStream( in );
            if ( FastaParser.isLikelyFasta( in ) ) {
                msa = FastaParser.parseMsa( is );
            }
            else {
                msa = GeneralMsaParser.parse( is );
            }
            MsaCompactor mc = null;
            if ( worst_remove > 0 ) {
                mc = MsaCompactor.removeWorstOffenders( msa, worst_remove, realign, norm );
            }
            else if ( av > 0 ) {
                mc = MsaCompactor.reduceGapAverage( msa, av, step, realign, out, 50 );
            }
            else if ( length > 0 ) {
                mc = MsaCompactor.reduceLength( msa, length, step, realign );
            }
            System.out.println( MsaMethods.calcGapRatio( mc.getMsa() ) );
            for( final String id : mc.getRemovedSeqIds() ) {
                System.out.println( id );
            }
            mc.writeMsa( out, MSA_FORMAT.PHYLIP, ".aln" );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
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
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " <options> <msa input file>" );
        System.out.println();
        System.out.println( " options: " );
        System.out.println();
        System.out.println( "   -" + REMOVE_WORST_OFFENDERS_OPTION + "=<integer>: step size (msa columns)" );
        System.out.println();
        System.out.println();
        System.out.println();
    }
}
