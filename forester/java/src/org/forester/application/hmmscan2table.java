
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.HmmscanPerDomainTableParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser.INDIVIDUAL_SCORE_CUTOFF;
import org.forester.protein.Protein;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class hmmscan2table {

    final static private String PRG_NAME                        = "hmmscan2table";
    final static private String PRG_VERSION                     = "1.0.0";
    final static private String PRG_DATE                        = "20220519";
    final static private String PRG_DESC                        = "hmmscan to table";
    final static private String E_MAIL                          = "phyloxml@gmail.com";
    final static private String WWW                             = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String HELP_OPTION_1                   = "help";
    final static private String HELP_OPTION_2                   = "h";
    private static final String MAX_I_E_VALUE_OPTION            = "ie";
    private static final String MIN_REL_ENV_LENGTH_RATIO_OPTION = "mrel";
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME,
                                              PRG_DESC,
                                              PRG_VERSION,
                                              PRG_DATE,
                                              E_MAIL,
                                              WWW,
                                              ForesterUtil.getForesterLibraryInformation() );
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) ) {
            System.out.println();
            print_help();
            System.exit( 0 );
        }
        if ( ( cla.getNumberOfNames() != 1 ) && ( cla.getNumberOfNames() != 2 ) ) {
            print_help();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<>();
        allowed_options.add( MAX_I_E_VALUE_OPTION );
        allowed_options.add( MIN_REL_ENV_LENGTH_RATIO_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        double rel_env_length_ratio_cutoff = -1;
        double ie_value_max = -1;
        if ( cla.isOptionSet( MIN_REL_ENV_LENGTH_RATIO_OPTION ) ) {
            try {
                rel_env_length_ratio_cutoff = cla.getOptionValueAsDouble( MIN_REL_ENV_LENGTH_RATIO_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, "no acceptable value for min rel env length ratio" );
            }
        }
        if ( cla.isOptionSet( MAX_I_E_VALUE_OPTION ) ) {
            try {
                ie_value_max = cla.getOptionValueAsDouble( MAX_I_E_VALUE_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, "no acceptable value for E-value maximum" );
            }
        }
        final File infile = cla.getFile( 0 );
        final HmmscanPerDomainTableParser parser = new HmmscanPerDomainTableParser( infile,
                                                                                    "unknown",
                                                                                    INDIVIDUAL_SCORE_CUTOFF.NONE );
        if ( ie_value_max > 0.0 ) {
            parser.setIEValueMaximum( ie_value_max );
        }
        if ( rel_env_length_ratio_cutoff > 0.0 ) {
            parser.setRelEnvLengthRatioCutoff( rel_env_length_ratio_cutoff );
        }
        List<Protein> protein_list = null;
        try {
            protein_list = parser.parse();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        catch ( final Exception e ) {
            ForesterUtil.unexpectedFatalError( PRG_NAME, e.getMessage(), e );
        }
        for( final Protein protein : protein_list ) {
            System.out.print( protein.getProteinId() );
            System.out.print( '\t' );
            System.out.print( protein.toDomainArchitectureString( "--", 1 ) );
            System.out.println();
        }
    }

    private static void print_help() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " [options] <hmmscan output> " );
        System.out.println();
        System.out.println( " options:" );
        System.out.println( MAX_I_E_VALUE_OPTION + ": max (inclusive) iE-value" );
        System.out.println( MIN_REL_ENV_LENGTH_RATIO_OPTION + ": min (inclusive) relative envelope length ratio" );
    }
}
