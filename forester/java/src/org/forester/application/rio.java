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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.application;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
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

    final static private String PRG_NAME                     = "rio";
    final static private String PRG_VERSION                  = "4.000 beta 11";
    final static private String PRG_DATE                     = "170410";
    final static private String E_MAIL                       = "phyloxml@gmail.com";
    final static private String WWW                          = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String HELP_OPTION_1                = "help";
    final static private String LOGFILE_SUFFIX               = "_RIO_log.tsv";
    final static private String STRIPPED_SPECIES_TREE_SUFFIX = "_RIO_sst.xml";
    final static private String ORTHO_OUTTABLE_SUFFIX        = "_RIO_o_table.tsv";
    final static private String OUT_GENE_TREE_SUFFIX         = "_RIO_gene_tree.xml";
    final static private String HELP_OPTION_2                = "h";
    final static private String GT_FIRST                     = "f";
    final static private String GT_LAST                      = "l";
    final static private String REROOTING_OPT                = "r";
    final static private String OUTGROUP                     = "o";
    final static private String RETURN_SPECIES_TREE          = "s";
    final static private String RETURN_BEST_GENE_TREE        = "g";
    final static private String USE_SDIR                     = "b";
    final static private String TRANSFER_TAXONOMY_OPTION     = "t";
    final static private String GENE_TREES_SUFFIX_OPTION     = "u";

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
        allowed_options.add( GENE_TREES_SUFFIX_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( "unknown option(s): " + dissallowed_options );
        }
        final File gene_trees_file = cla.getFile( 0 );
        final boolean use_dir;
        File indir = null;
        File outdir = null;
        if ( gene_trees_file.isDirectory() ) {
            if ( !gene_trees_file.exists() ) {
                ForesterUtil.fatalError( "gene trees directory \"" + gene_trees_file + "\" does not exist" );
            }
            use_dir = true;
            indir = gene_trees_file;
        }
        else {
            use_dir = false;
        }
        final File species_tree_file = cla.getFile( 1 );
        File orthology_outtable = null;
        if ( use_dir ) {
            outdir = cla.getFile( 2 );
        }
        else {
            orthology_outtable = cla.getFile( 2 );
        }
        File logfile;
        if ( use_dir ) {
            if ( ( cla.getNumberOfNames() < 4 ) ) {
                System.out.println();
                System.out.println( "error: incorrect number of arguments" );
                System.out.println();
                printHelp();
            }
            logfile = cla.getFile( 3 );
            if ( logfile.exists() ) {
                ForesterUtil.fatalError( "\"" + logfile + "\" already exists" );
            }
        }
        else {
            if ( cla.getNumberOfNames() > 3 ) {
                logfile = cla.getFile( 3 );
                if ( logfile.exists() ) {
                    ForesterUtil.fatalError( "\"" + logfile + "\" already exists" );
                }
            }
            else {
                logfile = null;
            }
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
            if ( sdir ) {
                ForesterUtil.fatalError( "no outgroup option for SDIR algorithm" );
            }
            if ( use_dir ) {
                ForesterUtil.fatalError( "no outgroup option for operating on gene trees directory" );
            }
            if ( !cla.isOptionHasAValue( OUTGROUP ) ) {
                ForesterUtil.fatalError( "no value for -" + OUTGROUP );
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
                if ( use_dir ) {
                    ForesterUtil.fatalError( "no outgroup option for operating on gene trees directory" );
                }
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
            if ( use_dir ) {
                ForesterUtil.fatalError( "no return species tree option when operating on gene trees directory" );
            }
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
            if ( use_dir ) {
                ForesterUtil.fatalError( "no best gene tree return option when operating on gene trees directory" );
            }
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
            if ( use_dir ) {
                ForesterUtil.fatalError( "no transferring taxonomy option when operating on gene trees directory" );
            }
            if ( return_gene_tree == null ) {
                ForesterUtil.fatalError( "no point in transferring taxonomy data without returning best gene tree" );
            }
            transfer_taxonomy = true;
        }
        if ( !use_dir ) {
            ForesterUtil.fatalErrorIfFileNotReadable( gene_trees_file );
        }
        else {
            transfer_taxonomy = true;
        }
        final String gene_trees_suffix;
        if ( cla.isOptionSet( GENE_TREES_SUFFIX_OPTION ) ) {
            if ( !use_dir ) {
                ForesterUtil.fatalError( "no gene tree suffix option when operating on indivual gene trees" );
            }
            if ( !cla.isOptionHasAValue( GENE_TREES_SUFFIX_OPTION ) ) {
                ForesterUtil.fatalError( "no value for -" + GENE_TREES_SUFFIX_OPTION );
            }
            gene_trees_suffix = cla.getOptionValueAsCleanString( GENE_TREES_SUFFIX_OPTION );
        }
        else {
            gene_trees_suffix = ".mlt";
        }
        ForesterUtil.fatalErrorIfFileNotReadable( species_tree_file );
        if ( !use_dir && orthology_outtable.exists() ) {
            ForesterUtil.fatalError( "\"" + orthology_outtable + "\" already exists" );
        }
        long time = 0;
        try {
            if ( use_dir ) {
                System.out.println( "Gene trees in-dir                   :\t" + indir.getCanonicalPath() );
                System.out.println( "Gene trees suffix                   :\t" + gene_trees_suffix );
            }
            else {
                System.out.println( "Gene trees                          :\t" + gene_trees_file.getCanonicalPath() );
            }
            System.out.println( "Species tree                        :\t" + species_tree_file.getCanonicalPath() );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( e.getLocalizedMessage() );
        }
        if ( use_dir ) {
            System.out.println( "Out-dir                             :\t" + outdir );
        }
        else {
            System.out.println( "All vs all orthology results table  :\t" + orthology_outtable );
        }
        if ( logfile != null ) {
            System.out.println( "Logfile                             :\t" + logfile );
        }
        if ( gt_first != RIO.DEFAULT_RANGE ) {
            System.out.println( "First gene tree to analyze          :\t" + gt_first );
        }
        if ( gt_last != RIO.DEFAULT_RANGE ) {
            System.out.println( "Last gene tree to analyze           :\t" + gt_last );
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
        System.out.println( "Re-rooting                          : \t" + rerooting_str );
        if ( !sdir ) {
            System.out.println( "Non binary species tree             :\tallowed" );
        }
        else {
            System.out.println( "Non binary species tree             :\tdisallowed" );
        }
        if ( return_species_tree != null ) {
            System.out.println( "Write used species tree to          :\t" + return_species_tree );
        }
        if ( return_gene_tree != null ) {
            System.out.println( "Write best gene tree to             :\t" + return_gene_tree );
            System.out.println( "Transfer taxonomic data             :\t" + transfer_taxonomy );
        }
        time = System.currentTimeMillis();
        final ALGORITHM algorithm;
        if ( sdir ) {
            algorithm = ALGORITHM.SDIR;
        }
        else {
            algorithm = ALGORITHM.GSDIR;
        }
        EasyWriter log = null;
        if ( use_dir ) {
            if ( outdir.exists() ) {
                if ( !outdir.isDirectory() ) {
                    ForesterUtil.fatalError( PRG_NAME,
                                             "out-directory [" + outdir + "] already exists but is not a directory" );
                }
            }
            else {
                final boolean success = outdir.mkdirs();
                if ( !success ) {
                    ForesterUtil.fatalError( PRG_NAME, "could not create out-directory [" + outdir + "]" );
                }
            }
            final String species_tree_file_name = species_tree_file.getName();
            final File gene_trees_files[] = indir.listFiles( new FilenameFilter() {

                @Override
                public boolean accept( final File dir, final String name ) {
                    return ( ( name.endsWith( gene_trees_suffix ) ) && !( name.equals( species_tree_file_name ) ) );
                }
            } );
            if ( gene_trees_files.length < 1 ) {
                ForesterUtil.fatalError( PRG_NAME,
                                         "in-directory [" + indir
                                                 + "] does not contain any gene tree files with suffix "
                                                 + gene_trees_suffix );
            }
            try {
                log = ForesterUtil.createEasyWriter( logfile );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, "could not create [" + logfile + "]" );
            }
            Arrays.sort( gene_trees_files );
            try {
                log.print( "# program" );
                log.print( "\t" );
                log.print( PRG_NAME );
                log.println();
                log.print( "# version" );
                log.print( "\t" );
                log.print( PRG_VERSION );
                log.println();
                log.print( "# date" );
                log.print( "\t" );
                log.print( PRG_DATE );
                log.println();
                log.print( "# Algorithm " );
                log.print( "\t" );
                log.print( algorithm.toString() );
                log.println();
                log.print( "# Gene trees in-dir" );
                log.print( "\t" );
                log.print( indir.getCanonicalPath() );
                log.println();
                log.print( "# Gene trees suffix" );
                log.print( "\t" );
                log.print( gene_trees_suffix );
                log.println();
                log.print( "# Species tree" );
                log.print( "\t" );
                log.print( species_tree_file.getCanonicalPath() );
                log.println();
                log.print( "# Out-dir" );
                log.print( "\t" );
                log.print( outdir.getCanonicalPath() );
                log.println();
                log.print( "# Logfile" );
                log.print( "\t" );
                log.print( logfile.getCanonicalPath() );
                log.println();
                if ( gt_first != RIO.DEFAULT_RANGE ) {
                    log.print( "# First gene tree to analyze" );
                    log.print( "\t" );
                    log.print( Integer.toString( gt_first ) );
                    log.println();
                }
                if ( gt_last != RIO.DEFAULT_RANGE ) {
                    log.print( "# Last gene tree to analyze" );
                    log.print( "\t" );
                    log.print( Integer.toString( gt_last ) );
                    log.println();
                }
                log.print( "# Re-rooting" );
                log.print( "\t" );
                log.print( rerooting_str );
                log.println();
                log.print( "# Non binary species tree" );
                log.print( "\t" );
                if ( !sdir ) {
                    log.print( "allowed" );
                }
                else {
                    log.print( "disallowed" );
                }
                log.println();
                log.println();
                log.print( "NAME" );
                log.print( "\t" );
                log.print( "EXT NODES" );
                log.print( "\t" );
                log.print( "MEAN DUP" );
                log.print( "\t" );
                log.print( "MEAN DUP SD" );
                log.print( "\t" );
                log.print( "MEDIAN DUP" );
                log.print( "\t" );
                log.print( "MIN DUP" );
                log.print( "\t" );
                log.print( "MAX DUP" );
                log.print( "\t" );
                log.print( "REMOVED EXT NODES" );
                log.print( "\t" );
                log.print( "N" );
                log.println();
            }
            catch ( IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
            }
            int counter = 1;
            for( final File gf : gene_trees_files ) {
                String outname = gf.getName();
                System.out
                        .print( "\r                                                                                            " );
                System.out.print( "\r" + counter + "/" + gene_trees_files.length + ": " + outname );
                counter++;
                if ( outname.indexOf( "." ) > 0 ) {
                    outname = outname.substring( 0, outname.lastIndexOf( "." ) );
                }
                try {
                    executeAnalysis( gf,
                                     species_tree_file,
                                     new File( outdir.getCanonicalFile() + "/" + outname + ORTHO_OUTTABLE_SUFFIX ),
                                     new File( outdir.getCanonicalFile() + "/" + outname + LOGFILE_SUFFIX ),
                                     outgroup,
                                     rerooting,
                                     gt_first,
                                     gt_last,
                                     new File( outdir.getCanonicalFile() + "/" + outname
                                             + STRIPPED_SPECIES_TREE_SUFFIX ),
                                     new File( outdir.getCanonicalFile() + "/" + outname + OUT_GENE_TREE_SUFFIX ),
                                     transfer_taxonomy,
                                     algorithm,
                                     true,
                                     log );
                }
                catch ( IOException e ) {
                    ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
                }
            }
            System.out
                    .print( "\r                                                                                        " );
            System.out.println();
        }
        else {
            executeAnalysis( gene_trees_file,
                             species_tree_file,
                             orthology_outtable,
                             logfile,
                             outgroup,
                             rerooting,
                             gt_first,
                             gt_last,
                             return_species_tree,
                             return_gene_tree,
                             transfer_taxonomy,
                             algorithm,
                             false,
                             null );
        }
        if ( !use_dir ) {
            time = System.currentTimeMillis() - time;
            System.out.println( "Time                                :\t" + time + "ms" );
        }
        else {
            try {
                log.close();
            }
            catch ( IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
            }
            time = System.currentTimeMillis() - time;
            System.out.println( "Time                                :\t" + time + "ms" );
        }
        System.exit( 0 );
    }

    private static final void executeAnalysis( final File gene_trees_file,
                                               final File species_tree_file,
                                               final File orthology_outtable,
                                               final File logfile,
                                               final String outgroup,
                                               final REROOTING rerooting,
                                               final int gt_first,
                                               final int gt_last,
                                               final File return_species_tree,
                                               final File return_gene_tree,
                                               final boolean transfer_taxonomy,
                                               final ALGORITHM algorithm,
                                               final boolean use_gene_trees_dir,
                                               final EasyWriter log ) {
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
                                           !use_gene_trees_dir,
                                           transfer_taxonomy );
            }
            if ( !use_gene_trees_dir ) {
                if ( algorithm == ALGORITHM.GSDIR ) {
                    System.out.println( "Taxonomy linking based on           :\t" + rio.getGSDIRtaxCompBase() );
                }
            }
            final IntMatrix m;
            if ( iterating ) {
                m = rio.getOrthologTable();
            }
            else {
                m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            }
            final BasicDescriptiveStatistics stats = rio.getDuplicationsStatistics();
            writeTable( orthology_outtable, stats.getN(), m, !use_gene_trees_dir );
            if ( ( algorithm != ALGORITHM.SDIR ) && ( logfile != null ) ) {
                writeLogFile( logfile,
                              rio,
                              species_tree_file,
                              gene_trees_file,
                              orthology_outtable,
                              PRG_NAME,
                              PRG_VERSION,
                              PRG_DATE,
                              ForesterUtil.getForesterLibraryInformation(),
                              !use_gene_trees_dir );
            }
            if ( return_species_tree != null ) {
                writeTree( rio.getSpeciesTree(),
                           return_species_tree,
                           use_gene_trees_dir ? null : "Wrote (stripped) species tree to    :\t" );
            }
            if ( return_gene_tree != null ) {
                writeTree( rio.getMinDuplicationsGeneTree(),
                           return_gene_tree,
                           use_gene_trees_dir ? null : "Wrote one min duplication gene tree :\t" );
            }
            final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.##" );
            final int min = ( int ) stats.getMin();
            final int max = ( int ) stats.getMax();
            final int median = ( int ) stats.median();
            int min_count = 0;
            int max_count = 0;
            int median_count = 0;
            for( double d : stats.getData() ) {
                if ( ( ( int ) d ) == min ) {
                    ++min_count;
                }
                if ( ( ( int ) d ) == max ) {
                    ++max_count;
                }
                if ( ( ( int ) d ) == median ) {
                    ++median_count;
                }
            }
            final double min_count_percentage = ( 100.0 * min_count ) / stats.getN();
            final double max_count_percentage = ( 100.0 * max_count ) / stats.getN();
            final double median_count_percentage = ( 100.0 * median_count ) / stats.getN();
            if ( use_gene_trees_dir ) {
                String name = gene_trees_file.getName();
                if ( name.indexOf( "." ) > 0 ) {
                    name = name.substring( 0, name.lastIndexOf( "." ) );
                }
                log.print( name );
                log.print( "\t" );
                log.print( Integer.toString( rio.getExtNodesOfAnalyzedGeneTrees() ) );
                log.print( "\t" );
                log.print( df.format( stats.arithmeticMean() ) );
                log.print( "\t" );
                log.print( df.format( stats.sampleStandardDeviation() ) );
                log.print( "\t" );
                if ( stats.getN() > 3 ) {
                    log.print( df.format( median ) );
                }
                else {
                    log.print( "" );
                }
                log.print( "\t" );
                log.print( Integer.toString( min ) );
                log.print( "\t" );
                log.print( Integer.toString( max ) );
                log.print( "\t" );
                log.print( Integer.toString( rio.getRemovedGeneTreeNodes().size() ) );
                log.print( "\t" );
                log.print( Integer.toString( stats.getN() ) );
                log.println();
            }
            else {
                System.out.println( "Gene tree internal nodes            :\t" + rio.getIntNodesOfAnalyzedGeneTrees() );
                System.out.println( "Gene tree external nodes            :\t" + rio.getExtNodesOfAnalyzedGeneTrees() );
                System.out.println( "Mean number of duplications         :\t" + df.format( stats.arithmeticMean() )
                        + "\t" + df.format( ( 100.0 * stats.arithmeticMean() ) / rio.getIntNodesOfAnalyzedGeneTrees() )
                        + "%\t(sd: " + df.format( stats.sampleStandardDeviation() ) + ")" );
                if ( stats.getN() > 3 ) {
                    System.out.println( "Median number of duplications       :\t" + df.format( median ) + "\t"
                            + df.format( ( 100.0 * median ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%" );
                }
                System.out.println( "Minimum duplications                :\t" + min + "\t"
                        + df.format( ( 100.0 * min ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%" );
                System.out.println( "Maximum duplications                :\t" + ( int ) max + "\t"
                        + df.format( ( 100.0 * max ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%" );
                System.out.println( "Gene trees with median duplications :\t" + median_count + "\t"
                        + df.format( median_count_percentage ) + "%" );
                System.out.println( "Gene trees with minimum duplications:\t" + min_count + "\t"
                        + df.format( min_count_percentage ) + "%" );
                System.out.println( "Gene trees with maximum duplications:\t" + max_count + "\t"
                        + df.format( max_count_percentage ) + "%" );
                System.out.println( "Removed ext gene tree nodes:\t" + rio.getRemovedGeneTreeNodes().size() );
            }
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
    }

    private final static void printHelp() {
        System.out.println( "Usage" );
        System.out.println();
        System.out.println( PRG_NAME
                + " [options] <gene trees infile> <species tree infile> <all vs all orthology table outfile> [logfile]" );
        System.out.println();
        System.out.println( PRG_NAME + " [options] <gene trees indir> <species tree infile> <outdir> <logfile>" );
        System.out.println();
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
        System.out.println( "  -" + TRANSFER_TAXONOMY_OPTION
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
        System.out.println( "  rio gene_trees.nh species.xml outtable.tsv log.txt" );
        System.out
                .println( "  rio -t -f=10 -l=100 -r=none -g=out_gene_tree.xml -s=stripped_species.xml gene_trees.xml species.xml outtable.tsv log.txt" );
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
                                      final String f,
                                      final boolean verbose )
            throws IOException {
        final EasyWriter out = ForesterUtil.createEasyWriter( logfile );
        out.println( "# " + prg_name );
        out.println( "# version : " + prg_v );
        out.println( "# date    : " + prg_date );
        out.println( "# based on: " + f );
        out.println( "# ----------------------------------" );
        out.println( "Gene trees                          :\t" + gene_trees_file.getCanonicalPath() );
        out.println( "Species tree                        :\t" + species_tree_file.getCanonicalPath() );
        out.println( "All vs all orthology table          :\t" + outtable.getCanonicalPath() );
        out.flush();
        out.println( rio.getLog().toString() );
        out.close();
        if ( verbose ) {
            System.out.println( "Wrote log to                        :\t" + logfile.getCanonicalPath() );
        }
    }

    private static void writeTable( final File table_outfile,
                                    final int gene_trees_analyzed,
                                    final IntMatrix m,
                                    final boolean verbose )
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
        if ( verbose ) {
            System.out.println( "Wrote table to                      :\t" + table_outfile.getCanonicalPath() );
        }
    }

    private static void writeTree( final Phylogeny p, final File f, final String comment ) throws IOException {
        final PhylogenyWriter writer = new PhylogenyWriter();
        writer.toPhyloXML( f, p, 0 );
        if ( comment != null ) {
            System.out.println( comment + f.getCanonicalPath() );
        }
    }
}
