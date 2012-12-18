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
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.datastructures.IntMatrix;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.rio.RIO;
import org.forester.rio.RIO.REROOTING;
import org.forester.rio.RIOException;
import org.forester.sdi.SDIException;
import org.forester.sdi.SDIutil.ALGORITHM;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.CommandLineArguments;
import org.forester.util.EasyWriter;
import org.forester.util.ForesterUtil;

public class rio {

    final static private String PRG_NAME      = "rio";
    final static private String PRG_VERSION   = "4.000 beta 3";
    final static private String PRG_DATE      = "2012.12.17";
    final static private String E_MAIL        = "czmasek@burnham.org";
    final static private String WWW           = "www.phylosoft.org/forester/";
    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";
    final static private String USE_SDIR      = "b";

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
        if ( ( args.length < 3 ) || ( args.length > 8 ) ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( USE_SDIR );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final File gene_trees_file = cla.getFile( 0 );
        final File species_tree_file = cla.getFile( 1 );
        final File orthology_outtable = cla.getFile( 2 );
        final File logfile;
        if ( cla.getNumberOfNames() > 3 ) {
            logfile = cla.getFile( 3 );
            if ( logfile.exists() ) {
                ForesterUtil.fatalError( PRG_NAME, "\"" + logfile + "\" already exists" );
            }
        }
        else {
            logfile = null;
        }
        final String outgroup = "";
        ForesterUtil.fatalErrorIfFileNotReadable( PRG_NAME, gene_trees_file );
        ForesterUtil.fatalErrorIfFileNotReadable( PRG_NAME, species_tree_file );
        if ( orthology_outtable.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "\"" + orthology_outtable + "\" already exists" );
        }
        boolean sdir = false;
        if ( cla.isOptionSet( USE_SDIR ) ) {
            sdir = true;
            if ( logfile != null ) {
                ForesterUtil.fatalError( PRG_NAME, "logfile output only for GSDIR algorithm" );
            }
        }
        long time = 0;
        System.out.println( "Gene trees                : " + gene_trees_file );
        System.out.println( "Species tree              : " + species_tree_file );
        System.out.println( "All vs all orthology table: " + orthology_outtable );
        if ( !sdir ) {
            if ( logfile != null ) {
                System.out.println( "Logfile                   : " + logfile );
            }
            System.out.println( "Non binary species tree   : allowed (GSDIR algorithm)" );
        }
        else {
            System.out.println( "Non binary species tree   : disallowed (SDIR algorithm)" );
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
            ForesterUtil.fatalError( PRG_NAME, "species tree is not rooted" );
        }
        final int o = PhylogenyMethods.countNumberOfOneDescendantNodes( species_tree );
        if ( o > 0 ) {
            ForesterUtil.printWarningMessage( PRG_NAME, "species tree has " + o
                    + " internal nodes with only one descendent! Going to strip them." );
            PhylogenyMethods.deleteInternalNodesWithOnlyOneDescendent( species_tree );
            if ( PhylogenyMethods.countNumberOfOneDescendantNodes( species_tree ) > 0 ) {
                ForesterUtil.unexpectedFatalError( PRG_NAME, "stripping of one-desc nodes failed" );
            }
        }
        final ALGORITHM algorithm;
        if ( sdir ) {
            algorithm = ALGORITHM.SDIR;
        }
        else {
            algorithm = ALGORITHM.GSDIR;
        }
        try {
            final RIO rio = RIO.executeAnalysis( gene_trees_file,
                                                 species_tree,
                                                 algorithm,
                                                 REROOTING.BY_ALGORITHM,
                                                 outgroup,
                                                 logfile != null,
                                                 true );
            if ( algorithm == ALGORITHM.GSDIR ) {
                ForesterUtil.programMessage( PRG_NAME, "taxonomy linking based on: " + rio.getGSDIRtaxCompBase() );
            }
            tableOutput( orthology_outtable, rio );
            if ( ( algorithm != ALGORITHM.SDIR ) && ( logfile != null ) ) {
                writeLogFile( logfile,
                              rio,
                              species_tree_file,
                              gene_trees_file,
                              orthology_outtable,
                              PRG_NAME,
                              PRG_VERSION,
                              PRG_DATE,
                              ForesterUtil.getForesterLibraryInformation() );
            }
            final BasicDescriptiveStatistics stats = rio.getDuplicationsStatistics();
            final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.#" );
            ForesterUtil.programMessage( PRG_NAME,
                                         "Mean number of duplications  : " + df.format( stats.arithmeticMean() )
                                                 + " (sd: " + df.format( stats.sampleStandardDeviation() ) + ")" );
            if ( stats.getN() > 3 ) {
                ForesterUtil.programMessage( PRG_NAME, "Median number of duplications: " + df.format( stats.median() ) );
            }
            ForesterUtil.programMessage( PRG_NAME, "Minimum duplications         : " + ( int ) stats.getMin() );
            ForesterUtil.programMessage( PRG_NAME, "Maximum duplications         : " + ( int ) stats.getMax() );
        }
        catch ( final RIOException e ) {
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
        time = System.currentTimeMillis() - time;
        ForesterUtil.programMessage( PRG_NAME, "time: " + time + "ms" );
        ForesterUtil.programMessage( PRG_NAME, "OK" );
        System.exit( 0 );
    }

    private final static void printHelp() {
        System.out.println( "Usage" );
        System.out.println();
        System.out
                .println( PRG_NAME
                        + " [options] <gene trees infile> <species tree infile> <all vs all orthology table outfile> [logfile]" );
        System.out.println();
        System.out.println( " Options" );
        System.out.println( "  -" + USE_SDIR
                + " : to use SDIR instead of GSDIR (faster, but non-binary species trees are disallowed)" );
        System.out.println();
        System.out.println( " Formats" );
        System.out.println( "  The species tree is expected to be in phyloXML format." );
        System.out
                .println( "  The gene trees ideally are in phyloXML as well, but can also be in New Hamphshire (Newick)" );
        System.out.println( "  or Nexus format as long as species information can be extracted from the gene names" );
        System.out.println( "  (e.g. \"HUMAN\" from \"BCL2_HUMAN\")." );
        System.out.println();
        System.out.println( " Examples" );
        System.out.println( "  \"rio gene_trees.nh species.xml outtable.tsv log.txt\"" );
        System.out.println();
        System.out.println( " More information: http://code.google.com/p/forester/wiki/RIO" );
        System.out.println();
        System.exit( -1 );
    }

    private static void tableOutput( final File table_outfile, final RIO rio ) throws IOException, RIOException {
        final IntMatrix m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
        writeTable( table_outfile, rio, m );
    }

    private static void writeLogFile( final File logfile,
                                      final RIO rio,
                                      final File species_tree_file,
                                      final File gene_trees_file,
                                      final File outtable,
                                      final String prg_name,
                                      final String prg_v,
                                      final String prg_date,
                                      final String f ) throws IOException {
        final EasyWriter out = ForesterUtil.createEasyWriter( logfile );
        out.println( prg_name );
        out.println( "version : " + prg_v );
        out.println( "date    : " + prg_date );
        out.println( "based on: " + f );
        out.println( "----------------------------------" );
        out.println( "Gene trees                                      : " + gene_trees_file );
        out.println( "Species tree                                    : " + species_tree_file );
        out.println( "All vs all orthology table                      : " + outtable );
        out.flush();
        out.println( rio.getLog().toString() );
        out.close();
        ForesterUtil.programMessage( PRG_NAME, "wrote log to \"" + logfile + "\"" );
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
                    if ( m.get( x, y ) != rio.getAnalyzedGeneTrees().length ) {
                        ForesterUtil.unexpectedFatalError( PRG_NAME, "diagonal value is off" );
                    }
                    w.print( "-" );
                }
                else {
                    w.print( df.format( ( ( double ) m.get( x, y ) ) / rio.getAnalyzedGeneTrees().length ) );
                }
            }
            w.println();
        }
        w.close();
        ForesterUtil.programMessage( PRG_NAME, "wrote table to \"" + table_outfile + "\"" );
    }
}
