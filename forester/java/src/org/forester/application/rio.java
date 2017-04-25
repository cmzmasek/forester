// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2017 Christian M. Zmasek
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.forester.rio.RIO;
import org.forester.rio.RIO.REROOTING;
import org.forester.rio.RIOUtil;
import org.forester.sdi.SDIutil.ALGORITHM;
import org.forester.util.CommandLineArguments;
import org.forester.util.EasyWriter;
import org.forester.util.ForesterUtil;

public class rio {
    //

    public final static String  PRG_NAME                       = "rio";
    public final static String  PRG_VERSION                    = "5.900";
    public final static String  PRG_DATE                       = "170420";
    final static private String E_MAIL                         = "phyloxml@gmail.com";
    final static private String WWW                            = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String HELP_OPTION_1                  = "help";
    final static private String HELP_OPTION_2                  = "h";
    final static private String GT_FIRST                       = "f";
    final static private String GT_LAST                        = "l";
    final static private String REROOTING_OPT                  = "r";
    final static private String OUTGROUP                       = "o";
    final static private String USE_SDIR                       = "s";
    final static private String GENE_TREES_SUFFIX_OPTION       = "g";
    final static private String MAPPINGS_DIR_OPTION            = "m";
    final static private String MAPPINGS_SUFFIX_OPTION         = "ms";
    final static private String CONSENSUS_TREES_DIR_OPTION     = "co";
    final static private String CONSENSUS_TREES_SUFFIX_OPTION  = "cos";
    final static private String MAPPINGS_SUFFIX_DEFAULT        = ".nim";
    final static private String CONSENSUS_TREE_SUFFIX_DEFAULT  = ".xml";
    final static private String ORTHOLOG_GROUPS_CUTOFF_OPTION  = "c";
    final static private String GENE_TREES_SUFFIX_DEFAULT      = ".mlt";
    final static private double ORTHOLOG_GROUPS_CUTOFF_DEFAULT = 0.5;

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
        allowed_options.add( GENE_TREES_SUFFIX_OPTION );
        allowed_options.add( ORTHOLOG_GROUPS_CUTOFF_OPTION );
        allowed_options.add( MAPPINGS_DIR_OPTION );
        allowed_options.add( MAPPINGS_SUFFIX_OPTION );
        allowed_options.add( CONSENSUS_TREES_DIR_OPTION );
        allowed_options.add( CONSENSUS_TREES_SUFFIX_OPTION );
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
            if ( gene_trees_file.listFiles().length < 1 ) {
                ForesterUtil.fatalError( "gene trees directory \"" + gene_trees_file + "\" is empty" );
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
            if ( !use_dir && logfile != null ) {
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
            if ( sdir ) {
                ForesterUtil.fatalError( "no gene tree range option for SDIR algorithm" );
            }
            if ( !cla.isOptionHasAValue( GT_FIRST ) ) {
                ForesterUtil.fatalError( "no value for -" + GT_FIRST );
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
            if ( sdir ) {
                ForesterUtil.fatalError( "no gene tree range option for SDIR algorithm" );
            }
            if ( !cla.isOptionHasAValue( GT_LAST ) ) {
                ForesterUtil.fatalError( "no value for -" + GT_LAST );
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
        double ortholog_group_cutoff = ORTHOLOG_GROUPS_CUTOFF_DEFAULT;
        if ( cla.isOptionSet( ORTHOLOG_GROUPS_CUTOFF_OPTION ) ) {
            if ( sdir ) {
                ForesterUtil.fatalError( "ortholog groups cutoff for SDIR algorithm" );
            }
            if ( !cla.isOptionHasAValue( ORTHOLOG_GROUPS_CUTOFF_OPTION ) ) {
                ForesterUtil.fatalError( "no value for -" + ORTHOLOG_GROUPS_CUTOFF_OPTION );
            }
            try {
                ortholog_group_cutoff = cla.getOptionValueAsDouble( ORTHOLOG_GROUPS_CUTOFF_OPTION );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( "could not parse double for -" + ORTHOLOG_GROUPS_CUTOFF_OPTION + " option" );
            }
            if ( ortholog_group_cutoff < 0 ) {
                ForesterUtil.fatalError( "attempt to set ortholog groups cutoff to: " + ortholog_group_cutoff );
            }
            if ( ortholog_group_cutoff > 1 ) {
                ForesterUtil.fatalError( "attempt to set ortholog groups cutoff to: " + ortholog_group_cutoff );
            }
        }
        if ( !use_dir ) {
            ForesterUtil.fatalErrorIfFileNotReadable( gene_trees_file );
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
            gene_trees_suffix = GENE_TREES_SUFFIX_DEFAULT;
        }
        final boolean perform_id_mapping;
        final File id_mapping_dir;
        if ( cla.isOptionSet( MAPPINGS_DIR_OPTION ) ) {
            id_mapping_dir = new File( cla.getOptionValue( MAPPINGS_DIR_OPTION ) );
            perform_id_mapping = true;
            if ( !use_dir ) {
                ForesterUtil.fatalError( "no id mapping when operating on indivual gene trees" );
            }
            if ( !id_mapping_dir.exists() ) {
                ForesterUtil.fatalError( "id mappings directory \"" + id_mapping_dir + "\" does not exist" );
            }
            if ( !id_mapping_dir.isDirectory() ) {
                ForesterUtil.fatalError( "id mappings directory \"" + id_mapping_dir + "\" is not a directory" );
            }
            if ( id_mapping_dir.listFiles().length < 1 ) {
                ForesterUtil.fatalError( "id mappings directory \"" + id_mapping_dir + "\" is empty" );
            }
        }
        else {
            id_mapping_dir = null;
            perform_id_mapping = false;
        }
        final String id_mapping_suffix;
        if ( cla.isOptionSet( MAPPINGS_SUFFIX_OPTION ) ) {
            if ( !use_dir ) {
                ForesterUtil.fatalError( "no id mapping file suffix option when operating on indivual gene trees" );
            }
            if ( !perform_id_mapping ) {
                ForesterUtil.fatalError( "no id mapping directory given" );
            }
            if ( !cla.isOptionHasAValue( MAPPINGS_SUFFIX_OPTION ) ) {
                ForesterUtil.fatalError( "no value for -" + MAPPINGS_SUFFIX_OPTION );
            }
            id_mapping_suffix = cla.getOptionValueAsCleanString( MAPPINGS_SUFFIX_OPTION );
        }
        else {
            id_mapping_suffix = MAPPINGS_SUFFIX_DEFAULT;
        }
        boolean perform_gsdir_on_best_tree;
        final File best_trees_indir;
        if ( cla.isOptionSet( CONSENSUS_TREES_DIR_OPTION ) ) {
            best_trees_indir = new File( cla.getOptionValue( CONSENSUS_TREES_DIR_OPTION ) );
            perform_gsdir_on_best_tree = true;
            if ( !use_dir ) {
                ForesterUtil
                        .fatalError( "no consensus (\"best\") gene tree GSDIR analysis when operating on individual gene trees" );
            }
            if ( !best_trees_indir.exists() ) {
                ForesterUtil.fatalError( "consensus (\"best\") gene tree directory \"" + best_trees_indir
                        + "\" does not exist" );
            }
            if ( !best_trees_indir.isDirectory() ) {
                ForesterUtil.fatalError( "consensus (\"best\") gene tree directory \"" + best_trees_indir
                        + "\" is not a directory" );
            }
            if ( best_trees_indir.listFiles().length < 1 ) {
                ForesterUtil
                        .fatalError( "consensus (\"best\") gene tree directory \"" + best_trees_indir + "\" is empty" );
            }
        }
        else {
            best_trees_indir = null;
            perform_gsdir_on_best_tree = false;
        }
        final String best_trees_suffix;
        if ( cla.isOptionSet( CONSENSUS_TREES_SUFFIX_OPTION ) ) {
            if ( !use_dir ) {
                ForesterUtil
                        .fatalError( "no consensus (\"best\") gene tree suffix option when operating on individual gene trees" );
            }
            if ( !perform_gsdir_on_best_tree ) {
                ForesterUtil.fatalError( "no consensus (\"best\") gene tree directory given" );
            }
            if ( !cla.isOptionHasAValue( CONSENSUS_TREES_SUFFIX_OPTION ) ) {
                ForesterUtil.fatalError( "no value for -" + CONSENSUS_TREES_SUFFIX_OPTION );
            }
            best_trees_suffix = cla.getOptionValueAsCleanString( CONSENSUS_TREES_SUFFIX_OPTION );
        }
        else {
            best_trees_suffix = CONSENSUS_TREE_SUFFIX_DEFAULT;
        }
        ////////////////////////////////
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
        if ( perform_id_mapping ) {
            try {
                System.out.println( "Id mappings in-dir                  :\t" + id_mapping_dir.getCanonicalPath() );
            }
            catch ( IOException e ) {
                ForesterUtil.fatalError( e.getLocalizedMessage() );
            }
            System.out.println( "Id mappings suffix                  :\t" + id_mapping_suffix );
        }
        if ( perform_gsdir_on_best_tree ) {
            try {
                System.out.println( "Consensus (\"best\") gene tree dir    :\t" + best_trees_indir.getCanonicalPath() );
            }
            catch ( IOException e ) {
                ForesterUtil.fatalError( e.getLocalizedMessage() );
            }
            System.out.println( "Consensus (\"best\") gene tree suffix :\t" + best_trees_suffix );
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
        System.out.println( "Ortholog groups cutoff              :\t" + ortholog_group_cutoff );
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
                log.print( "# Ortholog groups cutoff" );
                log.print( "\t" );
                log.print( Double.toString( ortholog_group_cutoff ) );
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
                log.print( ortholog_group_cutoff + " O GROUPS" );
                log.print( "\t" );
                log.print( "0.05 O GROUPS" );
                log.print( "\t" );
                log.print( "0.25 O GROUPS" );
                log.print( "\t" );
                log.print( "0.5 O GROUPS" );
                log.print( "\t" );
                log.print( "0.75 O GROUPS" );
                log.print( "\t" );
                log.print( "0.95 O GROUPS" );
                log.print( "\t" );
                if ( perform_gsdir_on_best_tree ) {
                    log.print( "BEST TREE DUP" );
                    log.print( "\t" );
                    log.print( "MEDIAN DUP - BEST TREE DUP" );
                    log.print( "\t" );
                }
                log.print( "MEDIAN DUP" );
                log.print( "\t" );
                log.print( "MEAN DUP" );
                log.print( "\t" );
                log.print( "MEAN DUP SD" );
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
                    RIOUtil.executeAnalysis( gf,
                                             species_tree_file,
                                             new File( outdir.getCanonicalFile() + "/" + outname
                                                     + RIOUtil.ORTHO_OUTTABLE_SUFFIX ),
                                             new File( outdir.getCanonicalFile() + "/" + outname
                                                     + RIOUtil.ORTHO_OUTTABLE_WITH_MAP_SUFFIX ),
                                             new File( outdir.getCanonicalFile() + "/" + outname
                                                     + RIOUtil.ORTHOLOG_GROUPS_SUFFIX ),
                                             new File( outdir.getCanonicalFile() + "/" + outname
                                                     + RIOUtil.LOGFILE_SUFFIX ),
                                             outgroup,
                                             rerooting,
                                             gt_first,
                                             gt_last,
                                             new File( outdir.getCanonicalFile() + "/" + outname
                                                     + RIOUtil.STRIPPED_SPECIES_TREE_SUFFIX ),
                                             new File( outdir.getCanonicalFile() + "/" + outname
                                                     + RIOUtil.OUT_MIN_DUP_GENE_TREE_SUFFIX ),
                                             new File( outdir.getCanonicalFile() + "/" + outname
                                                     + RIOUtil.OUT_MED_DUP_GENE_TREE_SUFFIX ),
                                             true,
                                             algorithm,
                                             true,
                                             log,
                                             ortholog_group_cutoff,
                                             perform_id_mapping,
                                             id_mapping_dir,
                                             id_mapping_suffix,
                                             perform_gsdir_on_best_tree,
                                             outdir,
                                             best_trees_indir,
                                             best_trees_suffix );
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
            String outname = ForesterUtil.removeFileExtension( orthology_outtable.toString() );
            RIOUtil.executeAnalysis( gene_trees_file,
                                     species_tree_file,
                                     orthology_outtable,
                                     null,
                                     new File( outname + RIOUtil.ORTHOLOG_GROUPS_SUFFIX ),
                                     logfile,
                                     outgroup,
                                     rerooting,
                                     gt_first,
                                     gt_last,
                                     new File( outname + RIOUtil.STRIPPED_SPECIES_TREE_SUFFIX ),
                                     new File( outname + RIOUtil.OUT_MIN_DUP_GENE_TREE_SUFFIX ),
                                     new File( outname + RIOUtil.OUT_MED_DUP_GENE_TREE_SUFFIX ),
                                     algorithm == ALGORITHM.GSDIR,
                                     algorithm,
                                     false,
                                     null,
                                     ortholog_group_cutoff,
                                     false,
                                     null,
                                     null,
                                     false,
                                     null,
                                     null,
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
        System.out.println( "  -" + ORTHOLOG_GROUPS_CUTOFF_OPTION
                + "=<cutoff>    : cutoff value for ortholog groups (default: " + ORTHOLOG_GROUPS_CUTOFF_DEFAULT + ")" );
        System.out.println( "  -" + REROOTING_OPT
                + "=<re-rooting>: re-rooting method for gene trees, possible values or 'none', 'midpoint'," );
        System.out.println( "                   or 'outgroup' (default: by minizming duplications)" );
        System.out.println( "  -" + OUTGROUP
                + "=<outgroup>  : for rooting by outgroup, name of outgroup (external gene tree node)" );
        System.out.println( "  -" + USE_SDIR
                + "             : to use SDIR instead of GSDIR (faster, but non-binary species trees are" );
        System.out.println( "                   disallowed, as are most options)" );
        System.out.println( "  -" + GENE_TREES_SUFFIX_OPTION
                + "=<suffix>    : suffix for gene trees when operating on gene tree directories (default: "
                + GENE_TREES_SUFFIX_DEFAULT + ")" );
        System.out.println( "  -" + MAPPINGS_DIR_OPTION + "=<dir>       : directory for id mapping files" );
        System.out.println( "  -" + MAPPINGS_SUFFIX_OPTION + "=<suffix>   : suffix for id mapping files (default: "
                + MAPPINGS_SUFFIX_DEFAULT + ")" );
        System.out.println( "  -" + CONSENSUS_TREES_DIR_OPTION
                + "=<dir>      : directory with consenus (\"best\") gene trees to be analyzed with GSDIR" );
        System.out.println( "  -" + CONSENSUS_TREES_SUFFIX_OPTION
                + "=<suffix>  : suffix for consenus (\"best\") gene trees (default: " + CONSENSUS_TREE_SUFFIX_DEFAULT
                + ")" );
        ///
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
        System.out.println( "  rio -s gene_trees.nh species.xml outtable.tsv" );
        System.out.println( "  rio gene_trees.nh species.xml outtable.tsv log.txt" );
        System.out.println( "  rio -c=0.9 -f=10 -l=100 -r=none gene_trees.xml species.xml outtable.tsv log.txt" );
        System.out.println( "  rio -g=.xml gene_trees_dir species.xml out_dir log.tsv" );
        System.out.println( "  rio -g=.mlt -m=id_maps_dir -ms=.nim -c=0.8 gene_trees_dir species.xml out_dir log.tsv" );
        System.out.println( "  rio -m=id_maps_dir -c=0.8 gene_trees_dir species.xml out_dir log.tsv" );
        System.out
                .println( "  rio -m=id_maps_dir -co=consensus_dir -cos=.xml -c=0.8 gene_trees_dir species.xml out_dir log.tsv" );
        System.out.println();
        System.exit( -1 );
    }
}
