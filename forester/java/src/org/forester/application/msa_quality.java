
package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.Msa;
import org.forester.msa.MsaMethods;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.CommandLineArguments;
import org.forester.util.DescriptiveStatistics;

public class msa_quality {

    public static void main( final String args[] ) {
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            // ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        // if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length == 0 ) ) {
        //     printHelp();
        //     System.exit( 0 );
        // }
        final File in = cla.getFile( 0 );
        Msa msa = null;
        try {
            msa = GeneralMsaParser.parse( new FileInputStream( in ) );
        }
        catch ( final FileNotFoundException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        catch ( final IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        final int end = 2;
        final int start = 6;
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int c = start; c <= end; ++c ) {
            stats.addValue( MsaMethods.calculateIdentityRatio( msa, c ) );
        }
    }
}
