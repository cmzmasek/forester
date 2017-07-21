// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2017 Christian M. Zmasek
// Copyright (C) 2017 J. Craig Venter Institute
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
// Contact: phyloxml @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.application;

import java.io.File;
import java.text.DecimalFormat;

import org.forester.clade_analysis.Analysis;
import org.forester.clade_analysis.Result;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class cladinator {

    final static private String        PRG_NAME      = "cladinator";
    final static private String        PRG_VERSION   = "0.100";
    final static private String        PRG_DATE      = "170721";
    final static private String        PRG_DESC      = "clades within clades -- analysis of pplacer type outputs";
    final static private String        E_MAIL        = "phyloxml@gmail.com";
    final static private String        WWW           = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String        HELP_OPTION_1 = "help";
    final static private String        HELP_OPTION_2 = "h";
    private final static DecimalFormat df2           = new DecimalFormat( ".##" );

    public static void main( final String args[] ) {
        try {
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
            else if ( ( args.length != 2 ) ) {
                System.out.println();
                System.out.println( "Wrong number of arguments." );
                System.out.println();
                print_help();
                System.exit( -1 );
            }
            //final List<String> allowed_options = new ArrayList<>();
            final File intreefile = cla.getFile( 0 );
            final String query = cla.getName( 1 );
            System.out.println( "Input tree: " + intreefile );
            System.out.println( "Query:      " + query );
            Phylogeny p = null;
            try {
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intreefile, true );
                p = factory.create( intreefile, pp )[ 0 ];
            }
            catch ( final Exception e ) {
                System.out.println( "\nCould not read \"" + intreefile + "\" [" + e.getMessage() + "]\n" );
                System.exit( -1 );
            }
            final Result res = Analysis.execute( p, query );
            System.out.println();
            System.out.println( "Result:" );
            System.out.println( "Greatest common prefix     : " + res.getGreatestCommonPrefix() );
            System.out.println( "Greatest common prefix up  : " + res.getGreatestCommonPrefixUp() );
            System.out.println( "Greatest common prefix down: " + res.getGreatestCommonPrefixDown() );
            final double lec_ratio = ( 100.0 * res.getLeastEncompassingCladeSize() ) / res.getTreeSize();
            System.out.println( "Least Encompassing Clade has " + res.getLeastEncompassingCladeSize()
                    + " external nodes (" + df2.format( lec_ratio ) + "% of a total of " + res.getTreeSize() + ")" );
            if ( res.getWarnings().size() > 0 ) {
                System.out.println( "Warnings:" );
                for( final String s : res.getWarnings() ) {
                    System.out.println( s );
                }
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    private final static void print_help() {
        System.out.println( "Usage: " + PRG_NAME + " <gene tree file> <query>" );
        System.out.println();
    }
}
