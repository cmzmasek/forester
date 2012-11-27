// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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
// WWW: www.phylosoft.org/forester

package org.forester.application;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sdi.RIO;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class rio {

    final static private String PRG_NAME              = "rio";
    final static private String PRG_VERSION           = "3.00 beta 1";
    final static private String PRG_DATE              = "2010.01.15";
    final static private String E_MAIL                = "czmasek@burnham.org";
    final static private String WWW                   = "www.phylosoft.org/forester/";
    final static private String HELP_OPTION_1         = "help";
    final static private String HELP_OPTION_2         = "h";
    final static private String QUERY_OPTION          = "q";
    final static private String SORT_OPTION           = "s";
    final static private String OUTPUT_ULTRA_P_OPTION = "u";
    final static private String CUTOFF_ULTRA_P_OPTION = "cu";
    final static private String CUTOFF_ORTHO_OPTION   = "co";
    final static private String TABLE_OUTPUT_OPTION   = "t";

    public static void main( final String[] args ) {
        ForesterUtil.printProgramInformation( PRG_NAME,
                                              "resampled inference of orthologs",
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
        if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length == 0 ) ) {
            printHelp();
            System.exit( 0 );
        }
        if ( ( args.length < 3 ) || ( args.length > 10 ) ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( QUERY_OPTION );
        allowed_options.add( SORT_OPTION );
        allowed_options.add( CUTOFF_ULTRA_P_OPTION );
        allowed_options.add( CUTOFF_ORTHO_OPTION );
        allowed_options.add( TABLE_OUTPUT_OPTION );
        allowed_options.add( OUTPUT_ULTRA_P_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final File multiple_trees_file = cla.getFile( 0 );
        final File species_tree_file = cla.getFile( 1 );
        final File outfile = cla.getFile( 2 );
        ForesterUtil.fatalErrorIfFileNotReadable( PRG_NAME, multiple_trees_file );
        ForesterUtil.fatalErrorIfFileNotReadable( PRG_NAME, species_tree_file );
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        String seq_name = null;
        if ( cla.isOptionSet( QUERY_OPTION ) ) {
            seq_name = cla.getOptionValue( QUERY_OPTION );
        }
        File table_outfile = null;
        if ( cla.isOptionSet( TABLE_OUTPUT_OPTION ) ) {
            table_outfile = new File( cla.getOptionValue( TABLE_OUTPUT_OPTION ) );
            if ( table_outfile.exists() ) {
                ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
            }
        }
        boolean output_ultraparalogs = false;
        if ( cla.isOptionSet( OUTPUT_ULTRA_P_OPTION ) ) {
            output_ultraparalogs = true;
        }
        double t_orthologs = 0.0;
        double threshold_ultra_paralogs = 0.0;
        int sort = 2;
        try {
            if ( cla.isOptionSet( CUTOFF_ORTHO_OPTION ) ) {
                t_orthologs = cla.getOptionValueAsDouble( CUTOFF_ORTHO_OPTION );
            }
            if ( cla.isOptionSet( CUTOFF_ULTRA_P_OPTION ) ) {
                threshold_ultra_paralogs = cla.getOptionValueAsDouble( CUTOFF_ULTRA_P_OPTION );
            }
            if ( cla.isOptionSet( SORT_OPTION ) ) {
                sort = cla.getOptionValueAsInt( SORT_OPTION );
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "error in command line: " + e.getLocalizedMessage() );
        }
        if ( sort < 0 ) {
            sort = 0;
        }
        else if ( sort > 2 ) {
            sort = 2;
        }
        long time = 0;
        System.out.println( "\n" );
        System.out.println( "Gene trees:                              " + multiple_trees_file );
        System.out.println( "Species tree:                            " + species_tree_file );
        System.out.println( "Query:                                     " + seq_name );
        System.out.println( "Outfile:                                      " + outfile );
        System.out.println( "Outfile:                                      " + table_outfile );
        System.out.println( "Sort:                                         " + sort );
        System.out.println( "Threshold orthologs:                          " + t_orthologs );
        if ( output_ultraparalogs ) {
            System.out.println( "Threshold ultra paralogs:                     " + threshold_ultra_paralogs );
        }
        time = System.currentTimeMillis();
        Phylogeny species_tree = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            species_tree = factory.create( species_tree_file, new PhyloXmlParser() )[ 0 ];
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
        if ( !species_tree.isRooted() ) {
            ForesterUtil.printErrorMessage( PRG_NAME, "Species tree is not rooted" );
            System.exit( -1 );
        }
        if ( !species_tree.isCompletelyBinary() ) {
            ForesterUtil.printErrorMessage( PRG_NAME, "Species tree is not completely binary" );
            System.exit( -1 );
        }
        final RIO rio_instance = new RIO();
        final StringBuffer output = new StringBuffer();
        PrintWriter out = null;
        try {
            rio_instance.inferOrthologs( multiple_trees_file, species_tree.copy(), seq_name );
            output.append( rio_instance.inferredOrthologsToString( seq_name, sort, t_orthologs ) );
            if ( output_ultraparalogs ) {
                output.append( "\n\nUltra paralogs:\n" );
                output.append( rio_instance.inferredUltraParalogsToString( seq_name, threshold_ultra_paralogs ) );
            }
            output.append( "\n\nSort priority: " + RIO.getOrder( sort ) );
            output.append( "\nExt nodes    : " + rio_instance.getExtNodesOfAnalyzedGeneTrees() );
            output.append( "\nSamples      : " + rio_instance.getNumberOfSamples() + "\n" );
            out = new PrintWriter( new FileWriter( outfile ), true );
        }
        catch ( final Exception e ) {
            ForesterUtil.printErrorMessage( PRG_NAME, e.getLocalizedMessage() );
            e.printStackTrace();
            System.exit( -1 );
        }
        out.println( output );
        out.close();
        ForesterUtil.programMessage( PRG_NAME, "wrote results to \"" + outfile + "\"" );
        time = System.currentTimeMillis() - time;
        ForesterUtil.programMessage( PRG_NAME, "time: " + time + "ms" );
        ForesterUtil.programMessage( PRG_NAME, "OK." );
        System.exit( 0 );
    }

    private final static void printHelp() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " [options] <gene trees file> <species tree file> [outfile]" );
        System.out.println();
        System.out.println( "options:" );
        System.out.println();
        //        System.out.println( " -" + STRICT_OPTION
        //                + "    : strict [default: non-strict]: all nodes between 'target' and 'evaluators' must match" );
        //        System.out.println( " -" + NORMALIZE_OPTION
        //                + "=<d>: normalize to this value (e.g. 100 for most bootstrap analyses) [default: no normalization]" );
        //        System.out.println( " -" + FIRST_OPTION + "=<i>: first evaluator topology to use (0-based) [default: 0]" );
        //        System.out.println( " -" + LAST_OPTION
        //                + "=<i>: last evaluator topology to use (0-based) [default: use all until final topology]" );
        //        System.out.println();
        //        System.out.println( "M= (String) Multiple gene tree file (mandatory)" );
        //        System.out.println( "N= (String) Query sequence name (mandatory)" );
        //        System.out.println( "S= (String) Species tree file (mandatory)" );
        //        System.out.println( "O= (String) Output file name -- overwritten without warning! (mandatory)" );
        //        System.out.println( "P= (int)    Sort priority" );
        //        System.out.println( "L= (double) Threshold orthologs for output" );
        //        System.out.println( " Sort priority (\"P=\"):" );
        System.out.println( RIO.getOrderHelp().toString() );
        System.out.println();
        System.out
                .println( " Example: \"rio -q=D_NEMVE -s=1 -t=out -u Bcl-2_e1_20_mafft_05_40_fme.mlt species.xml out\"" );
        System.out.println();
    }
}
