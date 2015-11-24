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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;

import org.forester.datastructures.IntMatrix;
import org.forester.io.parsers.IteratingPhylogenyParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
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

    final static private String PRG_NAME                 = "rio";
    final static private String PRG_VERSION              = "4.000 beta 10";
    final static private String PRG_DATE                 = "140211";
    final static private String E_MAIL                   = "phyloxml@gmail.com";
    final static private String WWW                      = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String HELP_OPTION_1            = "help";
    final static private String HELP_OPTION_2            = "h";
    final static private String GT_FIRST                 = "f";
    final static private String GT_LAST                  = "l";
    final static private String REROOTING_OPT            = "r";
    final static private String OUTGROUP                 = "o";
    final static private String RETURN_SPECIES_TREE      = "s";
    final static private String RETURN_BEST_GENE_TREE    = "g";
    final static private String USE_SDIR                 = "b";
    final static private String TRANSFER_TAXONOMY_OPTION = "t";

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
            ForesterUtil.fatalError( e.getMessage() );
        }
        if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length == 0 ) ) {
            printHelp();
        }
        if ( ( args.length < 3 ) || ( args.length > 11 ) || ( cla.getNumberOfNames() < 3 ) ) {
            System.out.println();
            System.out.println( "error: incorrect number of arguments" );
            System.out.println();
            printHelp();
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( GT_FIRST );
        allowed_options.add( GT_LAST );
        allowed_options.add( REROOTING_OPT );
        allowed_options.add( OUTGROUP );
        allowed_options.add( USE_SDIR );
        allowed_options.add( RETURN_SPECIES_TREE );
        allowed_options.add( RETURN_BEST_GENE_TREE );
        allowed_options.add( TRANSFER_TAXONOMY_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( "unknown option(s): " + dissallowed_options );
        }
        final File gene_trees_file = cla.getFile( 0 );
        final File species_tree_file = cla.getFile( 1 );
        final File orthology_outtable = cla.getFile( 2 );
        final File logfile;
        if ( cla.getNumberOfNames() > 3 ) {
            logfile = cla.getFile( 3 );
            if ( logfile.exists() ) {
                ForesterUtil.fatalError( "\"" + logfile + "\" already exists" );
            }
        }
        else {
            logfile = null;
        }
        boolean sdir = false;
        if ( cla.isOptionSet( USE_SDIR ) ) {
            if ( cla.isOptionHasAValue( USE_SDIR ) ) {
                ForesterUtil.fatalError( "no value allowed for -" + USE_SDIR );
            }
            sdir = true;
            if ( logfile != null ) {
                ForesterUtil.fatalError( "no logfile output for SDIR algorithm" );
            }
        }
        String outgroup = null;
        if ( cla.isOptionSet( OUTGROUP ) ) {
            if ( !cla.isOptionHasAValue( OUTGROUP ) ) {
                ForesterUtil.fatalError( "no value for -" + OUTGROUP );
            }
            if ( sdir ) {
                ForesterUtil.fatalError( "no outgroup option for SDIR algorithm" );
            }
            outgroup = cla.getOptionValueAsCleanString( OUTGROUP );
        }
        REROOTING rerooting = REROOTING.BY_ALGORITHM;
        if ( cla.isOptionSet( REROOTING_OPT ) ) {
            if ( !cla.isOptionHasAValue( REROOTING_OPT ) ) {
                ForesterUtil.fatalError( "no value for -" + REROOTING_OPT );
            }
            if ( sdir ) {
                ForesterUtil.fatalError( "no re-rooting option for SDIR algorithm" );
            }
            final String rerooting_str = cla.getOptionValueAsCleanString( REROOTING_OPT ).toLowerCase();
            if ( rerooting_str.equals( "none" ) ) {
                rerooting = REROOTING.NONE;
            }
            else if ( rerooting_str.equals( "midpoint" ) ) {
                rerooting = REROOTING.MIDPOINT;
            }
            else if ( rerooting_str.equals( "outgroup" ) ) {
                rerooting = REROOTING.OUTGROUP;
            }
            else {
                ForesterUtil
                .fatalError( "values for re-rooting are: 'none', 'midpoint', or 'outgroup' (minizming duplications is default)" );
            }
        }
        if ( ForesterUtil.isEmpty( outgroup ) && ( rerooting == REROOTING.OUTGROUP ) ) {
            ForesterUtil.fatalError( "selected re-rooting by outgroup, but outgroup not set" );
        }
        if ( !ForesterUtil.isEmpty( outgroup ) && ( rerooting != REROOTING.OUTGROUP ) ) {
            ForesterUtil.fatalError( "outgroup set, but selected re-rooting by other approach" );
        }
        int gt_first = RIO.DEFAULT_RANGE;
        int gt_last = RIO.DEFAULT_RANGE;
        if ( cla.isOptionSet( GT_FIRST ) ) {
            if ( !cla.isOptionHasAValue( GT_FIRST ) ) {
                ForesterUtil.fatalError( "no value for -" + GT_FIRST );
            }
            if ( sdir ) {
                ForesterUtil.fatalError( "no gene tree range option for SDIR algorithm" );
            }
            try {
                gt_first = cla.getOptionValueAsInt( GT_FIRST );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( "could not parse integer for -" + GT_FIRST + " option" );
            }
            if ( gt_first < 0 ) {
                ForesterUtil.fatalError( "attempt to set index of first tree to analyze to: " + gt_first );
            }
        }
        if ( cla.isOptionSet( GT_LAST ) ) {
            if ( !cla.isOptionHasAValue( GT_LAST ) ) {
                ForesterUtil.fatalError( "no value for -" + GT_LAST );
            }
            if ( sdir ) {
                ForesterUtil.fatalError( "no gene tree range option for SDIR algorithm" );
            }
            try {
                gt_last = cla.getOptionValueAsInt( GT_LAST );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( "could not parse integer for -" + GT_LAST + " option" );
            }
            if ( gt_last < 0 ) {
                ForesterUtil.fatalError( "attempt to set index of last tree to analyze to: " + gt_last );
            }
        }
        if ( ( ( gt_last != RIO.DEFAULT_RANGE ) && ( gt_first != RIO.DEFAULT_RANGE ) ) && ( ( gt_last < gt_first ) ) ) {
            ForesterUtil.fatalError( "attempt to set range (0-based) of gene to analyze to: from " + gt_first + " to "
                    + gt_last );
        }
        File return_species_tree = null;
        if ( !sdir && cla.isOptionSet( RETURN_SPECIES_TREE ) ) {
            if ( !cla.isOptionHasAValue( RETURN_SPECIES_TREE ) ) {
                ForesterUtil.fatalError( "no value for -" + RETURN_SPECIES_TREE );
            }
            final String s = cla.getOptionValueAsCleanString( RETURN_SPECIES_TREE );
            return_species_tree = new File( s );
            if ( return_species_tree.exists() ) {
                ForesterUtil.fatalError( "\"" + return_species_tree + "\" already exists" );
            }
        }
        File return_gene_tree = null;
        if ( !sdir && cla.isOptionSet( RETURN_BEST_GENE_TREE ) ) {
            if ( !cla.isOptionHasAValue( RETURN_BEST_GENE_TREE ) ) {
                ForesterUtil.fatalError( "no value for -" + RETURN_BEST_GENE_TREE );
            }
            final String s = cla.getOptionValueAsCleanString( RETURN_BEST_GENE_TREE );
            return_gene_tree = new File( s );
            if ( return_gene_tree.exists() ) {
                ForesterUtil.fatalError( "\"" + return_gene_tree + "\" already exists" );
            }
        }
        boolean transfer_taxonomy = false;
        if ( !sdir && cla.isOptionSet( TRANSFER_TAXONOMY_OPTION ) ) {
            if ( return_gene_tree == null ) {
                ForesterUtil.fatalError( "no point in transferring taxonomy data without returning best gene tree" );
            }
            transfer_taxonomy = true;
        }
        ForesterUtil.fatalErrorIfFileNotReadable( gene_trees_file );
        ForesterUtil.fatalErrorIfFileNotReadable( species_tree_file );
        if ( orthology_outtable.exists() ) {
            ForesterUtil.fatalError( "\"" + orthology_outtable + "\" already exists" );
        }
        long time = 0;
        System.out.println( "Gene trees                : " + gene_trees_file );
        System.out.println( "Species tree              : " + species_tree_file );
        System.out.println( "All vs all orthology table: " + orthology_outtable );
        if ( logfile != null ) {
            System.out.println( "Logfile                   : " + logfile );
        }
        if ( gt_first != RIO.DEFAULT_RANGE ) {
            System.out.println( "First gene tree to analyze: " + gt_first );
        }
        if ( gt_last != RIO.DEFAULT_RANGE ) {
            System.out.println( "Last gene tree to analyze : " + gt_last );
        }
        String rerooting_str = "";
        switch ( rerooting ) {
            case BY_ALGORITHM: {
                rerooting_str = "by minimizing duplications";
                break;
            }
            case MIDPOINT: {
                rerooting_str = "by midpoint method";
                break;
            }
            case OUTGROUP: {
                rerooting_str = "by outgroup: " + outgroup;
                break;
            }
            case NONE: {
                rerooting_str = "none";
                break;
            }
        }
        System.out.println( "Re-rooting                : " + rerooting_str );
        if ( !sdir ) {
            System.out.println( "Non binary species tree   : allowed" );
        }
        else {
            System.out.println( "Non binary species tree   : disallowed" );
        }
        if ( return_species_tree != null ) {
            System.out.println( "Write used species tree to: " + return_species_tree );
        }
        if ( return_gene_tree != null ) {
            System.out.println( "Write best gene tree to   : " + return_gene_tree );
            System.out.println( "Transfer taxonomic data   : " + transfer_taxonomy );
        }
        time = System.currentTimeMillis();
        final ALGORITHM algorithm;
        if ( sdir ) {
            algorithm = ALGORITHM.SDIR;
        }
        else {
            algorithm = ALGORITHM.GSDIR;
        }
        try {
            final RIO rio;
            boolean iterating = false;
            final PhylogenyParser p = ParserUtils.createParserDependingOnFileType( gene_trees_file, true );
            if ( p instanceof PhyloXmlParser ) {
                rio = RIO.executeAnalysis( gene_trees_file,
                                           species_tree_file,
                                           algorithm,
                                           rerooting,
                                           outgroup,
                                           gt_first,
                                           gt_last,
                                           logfile != null,
                                           true,
                                           transfer_taxonomy );
            }
            else {
                iterating = true;
                if ( p instanceof NHXParser ) {
                    final NHXParser nhx = ( NHXParser ) p;
                    nhx.setReplaceUnderscores( false );
                    nhx.setIgnoreQuotes( true );
                    nhx.setTaxonomyExtraction( TAXONOMY_EXTRACTION.AGGRESSIVE );
                }
                else if ( p instanceof NexusPhylogeniesParser ) {
                    final NexusPhylogeniesParser nex = ( NexusPhylogeniesParser ) p;
                    nex.setReplaceUnderscores( false );
                    nex.setIgnoreQuotes( true );
                    nex.setTaxonomyExtraction( TAXONOMY_EXTRACTION.AGGRESSIVE );
                }
                else {
                    throw new RuntimeException( "unknown parser type: " + p );
                }
                final IteratingPhylogenyParser ip = ( IteratingPhylogenyParser ) p;
                ip.setSource( gene_trees_file );
                rio = RIO.executeAnalysis( ip,
                                           species_tree_file,
                                           algorithm,
                                           rerooting,
                                           outgroup,
                                           gt_first,
                                           gt_last,
                                           logfile != null,
                                           true,
                                           transfer_taxonomy );
            }
            if ( algorithm == ALGORITHM.GSDIR ) {
                System.out.println( "Taxonomy linking based on : " + rio.getGSDIRtaxCompBase() );
            }
            final IntMatrix m;
            if ( iterating ) {
                m = rio.getOrthologTable();
            }
            else {
                m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            }
            final BasicDescriptiveStatistics stats = rio.getDuplicationsStatistics();
            writeTable( orthology_outtable, stats.getN(), m );
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
            if ( return_species_tree != null ) {
                writeTree( rio.getSpeciesTree(), return_species_tree, "Wrote (stripped) species tree to" );
            }
            if ( return_gene_tree != null ) {
                String tt = "";
                if ( transfer_taxonomy ) {
                    tt = "(with transferred taxonomic data) ";
                }
                writeTree( rio.getMinDuplicationsGeneTree(),
                           return_gene_tree,
                           "Wrote (one) minimal duplication gene tree " + tt + "to" );
            }
            final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.#" );
            System.out.println( "Mean number of duplications  : " + df.format( stats.arithmeticMean() ) + " (sd: "
                    + df.format( stats.sampleStandardDeviation() ) + ") ("
                    + df.format( ( 100.0 * stats.arithmeticMean() ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%)" );
            if ( stats.getN() > 3 ) {
                System.out.println( "Median number of duplications: " + df.format( stats.median() ) + " ("
                        + df.format( ( 100.0 * stats.median() ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%)" );
            }
            System.out.println( "Minimum duplications         : " + ( int ) stats.getMin() + " ("
                    + df.format( ( 100.0 * stats.getMin() ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%)" );
            System.out.println( "Maximum duplications         : " + ( int ) stats.getMax() + " ("
                    + df.format( ( 100.0 * stats.getMax() ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%)" );
            System.out.println( "Gene tree internal nodes     : " + rio.getIntNodesOfAnalyzedGeneTrees() );
            System.out.println( "Gene tree external nodes     : " + rio.getExtNodesOfAnalyzedGeneTrees() );
        }
        catch ( final RIOException e ) {
            ForesterUtil.fatalError( e.getLocalizedMessage() );
        }
        catch ( final SDIException e ) {
            ForesterUtil.fatalError( e.getLocalizedMessage() );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( e.getLocalizedMessage() );
        }
        catch ( final OutOfMemoryError e ) {
            ForesterUtil.outOfMemoryError( e );
        }
        catch ( final Exception e ) {
            ForesterUtil.unexpectedFatalError( e );
        }
        catch ( final Error e ) {
            ForesterUtil.unexpectedFatalError( e );
        }
        time = System.currentTimeMillis() - time;
        System.out.println( "Time: " + time + "ms" );
        System.out.println( "OK" );
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
        System.out.println( "  -" + GT_FIRST + "=<first>     : first gene tree to analyze (0-based index)" );
        System.out.println( "  -" + GT_LAST + "=<last>      : last gene tree to analyze (0-based index)" );
        System.out.println( "  -" + REROOTING_OPT
                            + "=<re-rooting>: re-rooting method for gene trees, possible values or 'none', 'midpoint'," );
        System.out.println( "                   or 'outgroup' (default: by minizming duplications)" );
        System.out.println( "  -" + OUTGROUP
                            + "=<outgroup>  : for rooting by outgroup, name of outgroup (external gene tree node)" );
        System.out
        .println( "  -" + RETURN_SPECIES_TREE + "=<outfile>   : to write the (stripped) species tree to file" );
        System.out.println( "  -" + RETURN_BEST_GENE_TREE
                            + "=<outfile>   : to write (one) minimal duplication gene tree to file" );
        System.out
        .println( "  -"
                + TRANSFER_TAXONOMY_OPTION
                + "             : to transfer taxonomic data from species tree to returned minimal duplication gene tree\n"
                + "                   (if -" + RETURN_BEST_GENE_TREE + " option is used)" );
        System.out.println( "  -" + USE_SDIR
                            + "             : to use SDIR instead of GSDIR (faster, but non-binary species trees are" );
        System.out.println( "                   disallowed, as are most options)" );
        System.out.println();
        System.out.println( " Formats" );
        System.out
        .println( "  The gene trees, as well as the species tree, ideally are in phyloXML (www.phyloxml.org) format," );
        System.out
        .println( "  but can also be in New Hamphshire (Newick) or Nexus format as long as species information can be" );
        System.out
        .println( "  extracted from the gene names (e.g. \"HUMAN\" from \"BCL2_HUMAN\") and matched to a single species" );
        System.out.println( "  in the species tree." );
        System.out.println();
        System.out.println( " Examples" );
        System.out.println( "  \"rio gene_trees.nh species.xml outtable.tsv log.txt\"" );
        System.out.println();
        System.out.println( " More information: http://code.google.com/p/forester/wiki/RIO" );
        System.out.println();
        System.exit( -1 );
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
        System.out.println( "Wrote log to \"" + logfile + "\"" );
    }

    private static void writeTable( final File table_outfile, final int gene_trees_analyzed, final IntMatrix m )
            throws IOException {
        final EasyWriter w = ForesterUtil.createEasyWriter( table_outfile );
        final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.####" );
        df.setDecimalSeparatorAlwaysShown( false );
        df.setRoundingMode( RoundingMode.HALF_UP );
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
                    if ( m.get( x, y ) != gene_trees_analyzed ) {
                        ForesterUtil.unexpectedFatalError( "diagonal value is off" );
                    }
                    w.print( "-" );
                }
                else {
                    w.print( df.format( ( ( double ) m.get( x, y ) ) / gene_trees_analyzed ) );
                }
            }
            w.println();
        }
        w.close();
        System.out.println( "Wrote table to \"" + table_outfile + "\"" );
    }

    private static void writeTree( final Phylogeny p, final File f, final String comment ) throws IOException {
        final PhylogenyWriter writer = new PhylogenyWriter();
        writer.toPhyloXML( f, p, 0 );
        System.out.println( comment + " \"" + f + "\"" );
    }
}
