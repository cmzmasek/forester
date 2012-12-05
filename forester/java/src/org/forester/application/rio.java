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
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.forester.datastructures.IntMatrix;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sdi.RIO;
import org.forester.sdi.RioException;
import org.forester.sdi.SDIException;
import org.forester.util.CommandLineArguments;
import org.forester.util.EasyWriter;
import org.forester.util.ForesterUtil;

public class rio {

    final static private String PRG_NAME              = "rio";
    final static private String PRG_VERSION           = "3.00 beta 3";
    final static private String PRG_DATE              = "2012.12.05";
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
        }
        if ( ( args.length < 2 ) || ( args.length > 10 ) ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
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
        final File gene_trees_file = cla.getFile( 0 );
        final File species_tree_file = cla.getFile( 1 );
        File outfile = null;
        if ( cla.getNumberOfNames() > 2 ) {
            outfile = cla.getFile( 2 );
        }
        ForesterUtil.fatalErrorIfFileNotReadable( PRG_NAME, gene_trees_file );
        ForesterUtil.fatalErrorIfFileNotReadable( PRG_NAME, species_tree_file );
        if ( ( outfile != null ) && outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        String query = null;
        if ( cla.isOptionSet( QUERY_OPTION ) ) {
            query = cla.getOptionValue( QUERY_OPTION );
        }
        File table_outfile = null;
        if ( cla.isOptionSet( TABLE_OUTPUT_OPTION ) ) {
            table_outfile = new File( cla.getOptionValue( TABLE_OUTPUT_OPTION ) );
            if ( table_outfile.exists() ) {
                ForesterUtil.fatalError( PRG_NAME, "[" + table_outfile + "] already exists" );
            }
        }
        boolean output_ultraparalogs = false;
        if ( cla.isOptionSet( OUTPUT_ULTRA_P_OPTION ) ) {
            output_ultraparalogs = true;
        }
        double cutoff_for_orthologs = 50;
        double cutoff_for_ultra_paralogs = 50;
        int sort = 1;
        try {
            if ( cla.isOptionSet( CUTOFF_ORTHO_OPTION ) ) {
                cutoff_for_orthologs = cla.getOptionValueAsDouble( CUTOFF_ORTHO_OPTION );
                if ( query == null ) {
                    ForesterUtil.fatalError( PRG_NAME, "missing query name, type \"rio -h\" for help" );
                }
                if ( outfile == null ) {
                    ForesterUtil.fatalError( PRG_NAME, "missing outfile, type \"rio -h\" for help" );
                }
            }
            if ( cla.isOptionSet( CUTOFF_ULTRA_P_OPTION ) ) {
                cutoff_for_ultra_paralogs = cla.getOptionValueAsDouble( CUTOFF_ULTRA_P_OPTION );
                output_ultraparalogs = true;
            }
            if ( cla.isOptionSet( SORT_OPTION ) ) {
                sort = cla.getOptionValueAsInt( SORT_OPTION );
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "error in command line: " + e.getLocalizedMessage() );
        }
        if ( ( cutoff_for_orthologs < 0 ) || ( cutoff_for_ultra_paralogs < 0 ) || ( sort < 0 ) || ( sort > 2 ) ) {
            ForesterUtil.fatalError( PRG_NAME, "numberical option out of range, type \"rio -h\" for help" );
        }
        if ( ( ( query == null ) && ( ( outfile != null ) || output_ultraparalogs ) ) ) {
            ForesterUtil.fatalError( PRG_NAME, "missing query name, type \"rio -h\" for help" );
        }
        if ( ( output_ultraparalogs && ( outfile == null ) ) || ( ( query != null ) && ( outfile == null ) ) ) {
            ForesterUtil.fatalError( PRG_NAME, "missing outfile, type \"rio -h\" for help" );
        }
        long time = 0;
        System.out.println( "Gene trees                : " + gene_trees_file );
        System.out.println( "Species tree              : " + species_tree_file );
        if ( query != null ) {
            System.out.println( "Query                     : " + query );
            System.out.println( "Outfile                   : " + outfile );
            System.out.println( "Sort                      : " + sort );
            System.out.println( "Cutoff for  orthologs     : " + cutoff_for_orthologs );
            if ( output_ultraparalogs ) {
                System.out.println( "Cutoff for ultra paralogs : " + cutoff_for_ultra_paralogs );
            }
        }
        if ( table_outfile != null ) {
            System.out.println( "Table output              : " + table_outfile );
        }
        System.out.println();
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
            ForesterUtil.fatalError( PRG_NAME, "species tree is not rooted" );
        }
        if ( !species_tree.isCompletelyBinary() ) {
            ForesterUtil.fatalError( PRG_NAME, "species tree is not completely binary" );
        }
        try {
            final RIO rio;
            if ( ForesterUtil.isEmpty( query ) ) {
                rio = new RIO( gene_trees_file, species_tree );
            }
            else {
                rio = new RIO( gene_trees_file, species_tree, query );
            }
            if ( outfile != null ) {
                final StringBuilder output = new StringBuilder();
                output.append( rio.inferredOrthologsToString( query, sort, cutoff_for_orthologs ) );
                if ( output_ultraparalogs ) {
                    output.append( "\n\nUltra paralogs:\n" );
                    output.append( rio.inferredUltraParalogsToString( query, cutoff_for_ultra_paralogs ) );
                }
                output.append( "\n\nSort priority: " + RIO.getOrder( sort ) );
                output.append( "\nExt nodes    : " + rio.getExtNodesOfAnalyzedGeneTrees() );
                output.append( "\nSamples      : " + rio.getNumberOfSamples() + "\n" );
                final PrintWriter out = new PrintWriter( new FileWriter( outfile ), true );
                out.println( output );
                out.close();
            }
            if ( table_outfile != null ) {
                tableOutput( table_outfile, rio );
            }
        }
        catch ( final RioException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
        }
        catch ( final SDIException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
        }
        catch ( final Exception e ) {
            ForesterUtil.unexpectedFatalError( PRG_NAME, e );
        }
        if ( outfile != null ) {
            ForesterUtil.programMessage( PRG_NAME, "wrote results to \"" + outfile + "\"" );
        }
        time = System.currentTimeMillis() - time;
        ForesterUtil.programMessage( PRG_NAME, "time: " + time + "ms" );
        ForesterUtil.programMessage( PRG_NAME, "OK" );
        System.exit( 0 );
    }

    private static void tableOutput( final File table_outfile, final RIO rio ) throws IOException, RioException {
        final IntMatrix m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees() );
        writeTable( table_outfile, rio, m );
    }

    private static void writeTable( final File table_outfile, final RIO rio, final IntMatrix m ) throws IOException {
        final EasyWriter w = ForesterUtil.createEasyWriter( table_outfile );
        final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.###" );
        df.setDecimalSeparatorAlwaysShown( false );
        for( int i = 0; i < m.size(); ++i ) {
            w.print( "\t" );
            w.print( m.getLabel( i ) );
        }
        w.println();
        for( int x = 0; x < m.size(); ++x ) {
            w.print( m.getLabel( x ) );
            for( int y = 0; y < m.size(); ++y ) {
                w.print( "\t" );
                if ( x == y ) {
                    if ( m.get( x, y ) != rio.getNumberOfSamples() ) {
                        ForesterUtil.unexpectedFatalError( PRG_NAME, "diagonal value is off" );
                    }
                    w.print( "-" );
                }
                else {
                    w.print( df.format( ( ( double ) m.get( x, y ) ) / rio.getNumberOfSamples() ) );
                }
            }
            w.println();
        }
        w.close();
        ForesterUtil.programMessage( PRG_NAME, "wrote table to \"" + table_outfile + "\"" );
    }

    private final static void printHelp() {
        System.out.println( "Usage" );
        System.out.println();
        System.out.println( PRG_NAME + " [options] <gene trees file> <species tree file> [outfile]" );
        System.out.println();
        System.out.println( " Options" );
        System.out.println( "  -" + CUTOFF_ORTHO_OPTION + " : cutoff for ortholog output (default: 50)" );
        System.out.println( "  -" + TABLE_OUTPUT_OPTION
                + "  : file-name for output table of all vs. all ortholgy support" );
        System.out.println( "  -" + QUERY_OPTION
                + "  : name for query (sequence/node), if this is used, [outfile] is required as well" );
        System.out.println( "  -" + SORT_OPTION + "  : sort (default: 1)" );
        System.out.println( "  -" + OUTPUT_ULTRA_P_OPTION
                + "  : to output ultra-paralogs (species specific expansions/paralogs)" );
        System.out.println( "  -" + CUTOFF_ULTRA_P_OPTION + " : cutoff for ultra-paralog output (default: 50)" );
        System.out.println();
        System.out.println( " Note" );
        System.out.println( "  Either output of all vs. all ortholgy support with -t=<output table> and/or output for" );
        System.out.println( "  one query sequence with -q=<query name> and a [outfile] are required." );
        System.out.println();
        System.out.println( " Sort" );
        System.out.println( RIO.getOrderHelp().toString() );
        System.out.println( " Formats" );
        System.out.println( "  The species tree is expected to be in phyloXML format." );
        System.out
                .println( "  The gene trees ideally are in phyloXML as well, but can also be in New Hamphshire (Newick)" );
        System.out.println( "  or Nexus format as long as species information can be extracted from the gene names" );
        System.out.println( "  (e.g. \"HUMAN\" from \"BCL2_HUMAN\")." );
        System.out.println();
        System.out.println( " Examples" );
        System.out.println( "  \"rio gene_trees.nh species.xml outfile -q=BCL2_HUMAN -t=outtable -u -cu=60 -co=60\"" );
        System.out.println( "  \"rio gene_trees.nh species.xml -t=outtable\"" );
        System.out.println();
        System.out.println( " More information: http://code.google.com/p/forester/wiki/RIO" );
        System.out.println();
        System.exit( -1 );
    }
}
