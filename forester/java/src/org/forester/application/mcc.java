// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2012 Christian M. Zmasek
// Copyright (C) 2012 Sanford Burnham Medical Research Institute
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

package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.Msa;
import org.forester.msa.MsaMethods;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.CommandLineArguments;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class mcc {

    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";
    final static private String FROM_OPTION   = "f";
    final static private String TO_OPTION     = "t";
    final static private String STEP_OPTION   = "s";
    final static private String WINDOW_OPTION = "w";
    final static private String PRG_NAME      = "mcc";
    final static private String PRG_DESC      = "msa consensus conservation";
    final static private String PRG_VERSION   = "1.00";
    final static private String PRG_DATE      = "2012.05.18";
    final static private String E_MAIL        = "phylosoft@gmail.com";
    final static private String WWW           = "www.phylosoft.org/forester/";

    public static void main( final String args[] ) {
        try {
            final CommandLineArguments cla = new CommandLineArguments( args );
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length != 3 ) ) {
                printHelp();
                System.exit( 0 );
            }
            final File in = cla.getFile( 0 );
            int from = 0;
            int to = 0;
            int window = 0;
            int step = 0;
            if ( cla.isOptionSet( FROM_OPTION ) && cla.isOptionSet( TO_OPTION ) ) {
                from = cla.getOptionValueAsInt( FROM_OPTION );
                to = cla.getOptionValueAsInt( TO_OPTION );
            }
            else if ( cla.isOptionSet( STEP_OPTION ) && cla.isOptionSet( WINDOW_OPTION ) ) {
                step = cla.getOptionValueAsInt( STEP_OPTION );
                window = cla.getOptionValueAsInt( WINDOW_OPTION );
            }
            else {
                printHelp();
                System.exit( 0 );
            }
            Msa msa = null;
            final InputStream is = new FileInputStream( in );
            if ( FastaParser.isLikelyFasta( in ) ) {
                msa = FastaParser.parseMsa( is );
            }
            else {
                msa = GeneralMsaParser.parse( is );
            }
            if ( cla.isOptionSet( FROM_OPTION ) ) {
                singleCalc( in, from, to, msa );
            }
            else {
                windowedCalcs( window, step, msa );
            }
        }
        catch ( final Exception e ) {
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
        System.out.println( "   -" + FROM_OPTION + "=<integer>: from (msa column)" );
        System.out.println( "   -" + TO_OPTION + "=<integer>: to (msa column)" );
        System.out.println( "    or" );
        System.out.println( "   -" + WINDOW_OPTION + "=<integer>: window size (msa columns)" );
        System.out.println( "   -" + STEP_OPTION + "=<integer>: step size (msa columns)" );
        System.out.println();
        System.out.println();
        System.out.println();
    }

    private static void windowedCalcs( int window, int step, final Msa msa ) {
        if ( window < 1 ) {
            window = 1;
        }
        if ( step < 1 ) {
            step = 1;
        }
        final double id_ratios[] = new double[ msa.getLength() ];
        for( int i = 0; i <= ( msa.getLength() - 1 ); ++i ) {
            id_ratios[ i ] = MsaMethods.calculateIdentityRatio( msa, i );
        }
        String min_pos = "";
        String max_pos = "";
        double min = 1;
        double max = 0;
        for( int i = 0; i <= ( msa.getLength() - 1 ); i += step ) {
            int to = ( i + window ) - 1;
            if ( to > ( msa.getLength() - 1 ) ) {
                to = msa.getLength() - 1;
            }
            final DescriptiveStatistics stats = calc( i, to, id_ratios );
            final double mean = stats.arithmeticMean();
            final String pos = i + "-" + to;
            System.out.print( pos );
            System.out.print( ":\t" );
            System.out.print( mean );
            if ( stats.getN() > 2 ) {
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

    private static DescriptiveStatistics calc( final int from, final int to, final double id_ratios[] ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( int c = from; c <= to; ++c ) {
            stats.addValue( id_ratios[ c ] );
        }
        return stats;
    }
}
