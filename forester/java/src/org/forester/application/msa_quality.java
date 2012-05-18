
package org.forester.application;

import java.io.File;
import java.io.FileInputStream;

import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.Msa;
import org.forester.msa.MsaMethods;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.CommandLineArguments;
import org.forester.util.DescriptiveStatistics;

public class msa_quality {

    public static void main( final String args[] ) {
        try {
            final CommandLineArguments cla = new CommandLineArguments( args );
            // if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length == 0 ) ) {
            //     printHelp();
            //     System.exit( 0 );
            // }
            final File in = cla.getFile( 0 );
            int from = 0;
            int to = 0;
            int window = 0;
            int step = 0;
            if ( cla.isOptionSet( "f" ) && cla.isOptionSet( "t" ) ) {
                from = cla.getOptionValueAsInt( "f" );
                to = cla.getOptionValueAsInt( "t" );
            }
            else if ( cla.isOptionSet( "s" ) && cla.isOptionSet( "w" ) ) {
                step = cla.getOptionValueAsInt( "s" );
                window = cla.getOptionValueAsInt( "w" );
            }
            else {
            }
            Msa msa = null;
            msa = GeneralMsaParser.parse( new FileInputStream( in ) );
            if ( cla.isOptionSet( "f" ) && cla.isOptionSet( "t" ) ) {
                singleCalc( in, from, to, msa );
            }
            else {
                windowedCalcs( window, step, msa );
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            // ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    private static void windowedCalcs( int window, int step, final Msa msa ) {
        if ( window < 1 ) {
            window = 1;
        }
        if ( step < 1 ) {
            step = 1;
        }
        String min_pos = "";
        String max_pos = "";
        double min = 1;
        double max = 0;
        for( int i = 0; i <= msa.getLength() - 1; i += step ) {
            int to = i + window - 1;
            if ( to > ( msa.getLength() - 1 ) ) {
                to = msa.getLength() - 1;
            }
            final DescriptiveStatistics stats = calc( i, to, msa );
            final double mean = stats.arithmeticMean();
            final String pos = i + "-" + to;
            System.out.print( pos );
            System.out.print( ":\t" );
            System.out.print( mean );
            if ( window > 2 ) {
                System.out.print( "\t" );
                System.out.print( stats.median() );
                System.out.print( "\t" );
                System.out.print( stats.sampleStandardDeviation() );
            }
            System.out.println();
            if ( mean > max ) {
                max = mean;
                max_pos = pos;
            }
            if ( mean < min ) {
                min = mean;
                min_pos = pos;
            }
        }
        System.out.println( "Min: " + min_pos + ": " + min );
        System.out.println( "Max: " + max_pos + ": " + max );
    }

    private static void singleCalc( final File in, int from, int to, final Msa msa ) {
        if ( from < 0 ) {
            from = 0;
        }
        if ( to > ( msa.getLength() - 1 ) ) {
            to = msa.getLength() - 1;
        }
        final DescriptiveStatistics stats = calc( from, to, msa );
        System.out.println( in.toString() + ": " + from + "-" + to + ":" );
        System.out.println();
        System.out.println( stats.toString() );
    }

    private static DescriptiveStatistics calc( final int from, final int to, final Msa msa ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int c = from; c <= to; ++c ) {
            stats.addValue( MsaMethods.calculateIdentityRatio( msa, c ) );
        }
        return stats;
    }
}
