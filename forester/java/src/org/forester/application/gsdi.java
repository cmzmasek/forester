// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
import java.io.FilenameFilter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sdi.GSDI;
import org.forester.sdi.GSDII;
import org.forester.sdi.GSDIR;
import org.forester.sdi.SDIException;
import org.forester.sdi.SDIutil;
import org.forester.sdi.SDIutil.ALGORITHM;
import org.forester.util.CommandLineArguments;
import org.forester.util.EasyWriter;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public final class gsdi {

    final static public boolean REPLACE_UNDERSCORES_IN_NH_SPECIES_TREE = true;
    final static private String ALLOW_STRIPPING_OF_GENE_TREE_OPTION    = "g";
    final static private String GSDIR_OPTION                           = "r";
    final static private String MOST_PARSIMONIOUS_OPTION               = "m";
    final static private String SUFFIX_FOR_DIR_OPTION                  = "s";
    final static private String GUESS_FORMAT_OF_SPECIES_TREE           = "q";
    final static private String TRANSFER_TAXONOMY_OPTION               = "t";
    final static private String HELP_OPTION_1                          = "help";
    final static private String HELP_OPTION_2                          = "h";
    final static private String SUFFIX_FOR_SPECIES_TREE_USED           = "_species_tree_used.xml";
    final static private String OUTTREE_SUFFIX                         = "_gsdir.xml";
    final static private String LOGFILE_NAME                           = "00_gsdi_log.tsv";
    final static private String LOGFILE_SUFFIX                         = "_gsdi_log.txt";
    final static private String REMAPPED_SUFFIX                        = "_gsdi_remapped.txt";
    final static private String PRG_NAME                               = "gsdi";
    final static private String PRG_VERSION                            = "1.100";
    final static private String PRG_DATE                               = "170403";
    final static private String PRG_DESC                               = "general speciation duplication inference";
    final static private String E_MAIL                                 = "phyloxml@gmail.com";
    final static private String WWW                                    = "https://sites.google.com/site/cmzmasek/home/software/forester";

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
                gsdi.print_help();
                System.exit( 0 );
            }
            else if ( ( args.length < 2 ) || ( cla.getNumberOfNames() != 3 ) ) {
                System.out.println();
                System.out.println( "Wrong number of arguments." );
                System.out.println();
                gsdi.print_help();
                System.exit( -1 );
            }
            final List<String> allowed_options = new ArrayList<String>();
            allowed_options.add( GSDIR_OPTION );
            allowed_options.add( GUESS_FORMAT_OF_SPECIES_TREE );
            allowed_options.add( MOST_PARSIMONIOUS_OPTION );
            allowed_options.add( ALLOW_STRIPPING_OF_GENE_TREE_OPTION );
            allowed_options.add( TRANSFER_TAXONOMY_OPTION );
            allowed_options.add( SUFFIX_FOR_DIR_OPTION );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            execute( cla );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, e.getMessage() );
        }
    }

    private final static void execute( final CommandLineArguments cla ) throws IOException {
        ALGORITHM base_algorithm = ALGORITHM.GSDI;
        boolean most_parsimonous_duplication_model = false;
        boolean allow_stripping_of_gene_tree = false;
        if ( cla.isOptionSet( GSDIR_OPTION ) ) {
            base_algorithm = ALGORITHM.GSDIR;
        }
        if ( cla.isOptionSet( MOST_PARSIMONIOUS_OPTION ) ) {
            if ( base_algorithm == ALGORITHM.SDI ) {
                ForesterUtil.fatalError( PRG_NAME, "Cannot use most parsimonious duplication mode with SDI" );
            }
            most_parsimonous_duplication_model = true;
        }
        if ( cla.isOptionSet( ALLOW_STRIPPING_OF_GENE_TREE_OPTION ) ) {
            if ( base_algorithm == ALGORITHM.SDI ) {
                ForesterUtil.fatalError( PRG_NAME, "Cannot allow stripping of gene tree with SDI" );
            }
            allow_stripping_of_gene_tree = true;
        }
        boolean transfer_taxonomy = false;
        if ( cla.isOptionSet( TRANSFER_TAXONOMY_OPTION ) ) {
            transfer_taxonomy = true;
        }
        boolean use_gene_tree_dir = false;
        final String gene_tree_suffix;
        if ( cla.isOptionSet( SUFFIX_FOR_DIR_OPTION ) ) {
            gene_tree_suffix = cla.getOptionValue( SUFFIX_FOR_DIR_OPTION );
            use_gene_tree_dir = true;
        }
        else {
            gene_tree_suffix = null;
        }
        File gene_tree_file = null;
        File species_tree_file = null;
        File out_file = null;
        File log_file = null;
        File out_dir = null;
        try {
            gene_tree_file = cla.getFile( 0 );
            species_tree_file = cla.getFile( 1 );
            if ( use_gene_tree_dir ) {
                out_dir = cla.getFile( 2 );
                if ( out_dir.exists() ) {
                    if ( !out_dir.isDirectory() ) {
                        ForesterUtil
                                .fatalError( gsdi.PRG_NAME,
                                             "out-directory [" + out_dir + "] already exists but is not a directory" );
                    }
                }
                else {
                    final boolean success = out_dir.mkdirs();
                    if ( !success ) {
                        ForesterUtil.fatalError( gsdi.PRG_NAME, "could not create out-directory [" + out_dir + "]" );
                    }
                }
            }
            else {
                out_file = cla.getFile( 2 );
                log_file = new File( ForesterUtil.removeSuffix( out_file.toString() ) + LOGFILE_SUFFIX );
            }
        }
        catch ( final IllegalArgumentException e ) {
            ForesterUtil.fatalError( PRG_NAME, "error in command line: " + e.getMessage() );
        }
        if ( use_gene_tree_dir ) {
            final File indir = new File( gene_tree_file.toString() );
            if ( !indir.exists() ) {
                ForesterUtil.fatalError( gsdi.PRG_NAME, "in-directory [" + indir + "] does not exist" );
            }
            if ( !indir.isDirectory() ) {
                ForesterUtil.fatalError( gsdi.PRG_NAME, "in-directory [" + indir + "] is not a directory" );
            }
            final String species_tree_file_name = species_tree_file.getName();
            final File gene_tree_files[] = indir.listFiles( new FilenameFilter() {

                @Override
                public boolean accept( final File dir, final String name ) {
                    return ( ( name.endsWith( gene_tree_suffix ) ) && !( name.equals( species_tree_file_name ) ) );
                }
            } );
            if ( gene_tree_files.length < 1 ) {
                ForesterUtil.fatalError( gsdi.PRG_NAME,
                                         "in-directory [" + indir
                                                 + "] does not contain any gene tree files with suffix "
                                                 + gene_tree_suffix );
            }
            executeDir( base_algorithm,
                        most_parsimonous_duplication_model,
                        allow_stripping_of_gene_tree,
                        transfer_taxonomy,
                        gene_tree_files,
                        species_tree_file,
                        out_dir );
        }
        else {
            execute( base_algorithm,
                     most_parsimonous_duplication_model,
                     allow_stripping_of_gene_tree,
                     transfer_taxonomy,
                     gene_tree_file,
                     species_tree_file,
                     out_file,
                     log_file );
        }
    }

    private final static void executeDir( final ALGORITHM base_algorithm,
                                          final boolean most_parsimonous_duplication_model,
                                          final boolean allow_stripping_of_gene_tree,
                                          final boolean transfer_taxonomy,
                                          final File gene_tree_files[],
                                          final File species_tree_file,
                                          final File outdir )
            throws IOException {
        final File log_file = new File( outdir, LOGFILE_NAME );
        if ( ForesterUtil.isWritableFile( log_file ) != null ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, ForesterUtil.isWritableFile( log_file ) );
        }
        EasyWriter log_writer = null;
        try {
            log_writer = ForesterUtil.createEasyWriter( log_file );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, "Failed to create [" + log_file + "]: " + e.getMessage() );
        }
        log_writer.println( "# " + PRG_NAME );
        log_writer.println( "# Version\t" + PRG_VERSION );
        log_writer.println( "# Date\t" + PRG_DATE );
        log_writer.println( "# Forester version\t" + ForesterConstants.FORESTER_VERSION );
        log_writer.println( "# Species tree\t" + species_tree_file.getCanonicalPath() );
        if ( base_algorithm == ALGORITHM.GSDI ) {
            log_writer.println( "# Algorithm\tGSDI" );
        }
        else if ( base_algorithm == ALGORITHM.GSDIR ) {
            log_writer.println( "# Algorithm\tGSDIR" );
        }
        log_writer.println( "# Use most parsimonous duplication model\t" + most_parsimonous_duplication_model );
        log_writer.println( "# Allow stripping of gene tree nodes\t" + allow_stripping_of_gene_tree );
        log_writer.println( "# Start time\t" + new SimpleDateFormat( "yyyyMMdd HH:mm:ss" ).format( new Date() ) );
        log_writer.println();
        log_writer.print( "Gene-tree file\t" );
        log_writer.print( "Gene-tree name/#\t" );
        log_writer.print( "Ext. nodes\t" );
        log_writer.print( "Speciations\t" );
        log_writer.print( "Duplications\t" );
        if ( !most_parsimonous_duplication_model ) {
            log_writer.print( "Spec. or Dup.\t" );
        }
        if ( allow_stripping_of_gene_tree ) {
            log_writer.print( "Stripped gene-tree ext. nodes\t" );
        }
        log_writer.print( "Taxonomy mapping" );
        log_writer.println();
        int counter = 0;
        Arrays.sort( gene_tree_files );
        for( final File gene_tree_file : gene_tree_files ) {
            String outname = gene_tree_file.getName();
            if ( outname.indexOf( "." ) > 0 ) {
                outname = outname.substring( 0, outname.lastIndexOf( "." ) );
            }
            outname = outname + OUTTREE_SUFFIX;
            counter += executeOneTreeInDir( base_algorithm,
                                            most_parsimonous_duplication_model,
                                            allow_stripping_of_gene_tree,
                                            transfer_taxonomy,
                                            gene_tree_file,
                                            species_tree_file,
                                            new File( outdir, outname ),
                                            log_writer );
            log_writer.flush();
            System.out.print( "\r" + counter );
        }
        System.out.print( "\r" );
        log_writer.close();
        System.out.println( "Analyzed " + counter + " gene trees" );
        System.out.println();
        System.out.println( "Wrote log to: " + log_file.getCanonicalPath() );
        System.out.println();
    }

    private final static int executeOneTreeInDir( final ALGORITHM base_algorithm,
                                                  final boolean most_parsimonous_duplication_model,
                                                  final boolean allow_stripping_of_gene_tree,
                                                  final boolean transfer_taxonomy,
                                                  final File gene_tree_file,
                                                  final File species_tree_file,
                                                  final File out_file,
                                                  final EasyWriter log_writer )
            throws IOException {
        if ( ForesterUtil.isReadableFile( gene_tree_file ) != null ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, ForesterUtil.isReadableFile( gene_tree_file ) );
        }
        if ( ForesterUtil.isReadableFile( species_tree_file ) != null ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, ForesterUtil.isReadableFile( species_tree_file ) );
        }
        if ( ForesterUtil.isWritableFile( out_file ) != null ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, ForesterUtil.isWritableFile( out_file ) );
        }
        Phylogeny gene_trees[] = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            gene_trees = factory.create( gene_tree_file, PhyloXmlParser.createPhyloXmlParserXsdValidating() );
        }
        catch ( final IOException e ) {
            fatalError( "error",
                        "failed to read gene tree from [" + gene_tree_file + "]: " + e.getMessage(),
                        log_writer );
        }
        int counter = 0;
        final List<Phylogeny> out_trees = new ArrayList<Phylogeny>();
        for( final Phylogeny gene_tree : gene_trees ) {
            if ( !gene_tree.isEmpty() && gene_tree.getNumberOfExternalNodes() > 1 ) {
                Phylogeny species_tree = null;
                try {
                    species_tree = SDIutil.parseSpeciesTree( gene_tree,
                                                             species_tree_file,
                                                             REPLACE_UNDERSCORES_IN_NH_SPECIES_TREE,
                                                             true,
                                                             TAXONOMY_EXTRACTION.NO );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    fatalError( "user error",
                                "failed to transfer general node name, in [" + species_tree_file + "]: "
                                        + e.getMessage(),
                                log_writer );
                }
                catch ( final SDIException e ) {
                    fatalError( "user error", e.getMessage(), log_writer );
                }
                catch ( final IOException e ) {
                    fatalError( "error",
                                "Failed to read species tree from [" + species_tree_file + "]: " + e.getMessage(),
                                log_writer );
                }
                gene_tree.setRooted( true );
                species_tree.setRooted( true );
                if ( !gene_tree.isCompletelyBinary() ) {
                    fatalError( "user error",
                                "gene tree [" + gene_tree_file + "] is not completely binary",
                                log_writer );
                }
                if ( base_algorithm == ALGORITHM.SDI ) {
                    if ( !species_tree.isCompletelyBinary() ) {
                        fatalError( "user error",
                                    "species tree is not completely binary, use GSDI or GSDIR instead",
                                    log_writer );
                    }
                }
                log_writer.print( gene_tree_file.getName() );
                log_writer.print( "\t" );
                log_writer.print( ( ForesterUtil.isEmpty( gene_tree.getName() ) ? "" : gene_tree.getName() ) );
                if ( gene_trees.length > 1 ) {
                    log_writer.print( ( ForesterUtil.isEmpty( gene_tree.getName() ) ? Integer.toString( counter )
                            : ( ":" + Integer.toString( counter ) ) ) );
                }
                log_writer.print( "\t" );
                GSDII gsdii = null;
                try {
                    if ( base_algorithm == ALGORITHM.GSDI ) {
                        gsdii = new GSDI( gene_tree,
                                          species_tree,
                                          most_parsimonous_duplication_model,
                                          allow_stripping_of_gene_tree,
                                          true,
                                          transfer_taxonomy );
                    }
                    else if ( base_algorithm == ALGORITHM.GSDIR ) {
                        gsdii = new GSDIR( gene_tree,
                                           species_tree,
                                           allow_stripping_of_gene_tree,
                                           true,
                                           transfer_taxonomy );
                    }
                }
                catch ( final SDIException e ) {
                    fatalError( "user error", e.getLocalizedMessage(), log_writer );
                }
                catch ( final OutOfMemoryError e ) {
                    ForesterUtil.outOfMemoryError( e );
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                    fatalError( "unexpected error", e.toString(), log_writer );
                }
                if ( base_algorithm == ALGORITHM.GSDIR ) {
                    final Phylogeny gt = ( ( GSDIR ) gsdii ).getMinDuplicationsSumGeneTree();
                    gt.setRerootable( false );
                    out_trees.add( gt );
                }
                else {
                    gene_tree.setRerootable( false );
                    out_trees.add( gene_tree );
                }
                log_writer.print( gene_tree.getNumberOfExternalNodes() + "\t" );
                log_writer.print( gsdii.getSpeciationsSum() + "\t" );
                if ( ( base_algorithm == ALGORITHM.GSDIR ) ) {
                    final GSDIR gsdir = ( GSDIR ) gsdii;
                    log_writer.print( gsdir.getMinDuplicationsSum() + "\t" );
                }
                else if ( ( base_algorithm == ALGORITHM.GSDI ) ) {
                    final GSDI gsdi = ( GSDI ) gsdii;
                    log_writer.print( gsdi.getDuplicationsSum() + "\t" );
                    if ( !most_parsimonous_duplication_model ) {
                        log_writer.print( gsdi.getSpeciationOrDuplicationEventsSum() + "\t" );
                    }
                }
                if ( allow_stripping_of_gene_tree ) {
                    log_writer.print( gsdii.getStrippedExternalGeneTreeNodes().size() + "\t" );
                }
                log_writer.print( gsdii.getTaxCompBase().toString() );
                log_writer.println();
                ++counter;
            }
        }
        if ( counter > 0 ) {
            try {
                final PhylogenyWriter writer = new PhylogenyWriter();
                writer.toPhyloXML( out_file, out_trees, 0, ForesterUtil.LINE_SEPARATOR );
            }
            catch ( final IOException e ) {
                ForesterUtil
                        .fatalError( PRG_NAME,
                                     "Failed to write to [" + out_file.getCanonicalPath() + "]: " + e.getMessage() );
            }
        }
        return counter;
    }

    private final static void execute( final ALGORITHM base_algorithm,
                                       final boolean most_parsimonous_duplication_model,
                                       final boolean allow_stripping_of_gene_tree,
                                       final boolean transfer_taxonomy,
                                       final File gene_tree_file,
                                       final File species_tree_file,
                                       final File out_file,
                                       final File log_file )
            throws IOException {
        if ( ForesterUtil.isReadableFile( gene_tree_file ) != null ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, ForesterUtil.isReadableFile( gene_tree_file ) );
        }
        if ( ForesterUtil.isReadableFile( species_tree_file ) != null ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, ForesterUtil.isReadableFile( species_tree_file ) );
        }
        if ( ForesterUtil.isWritableFile( out_file ) != null ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, ForesterUtil.isWritableFile( out_file ) );
        }
        if ( ForesterUtil.isWritableFile( log_file ) != null ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, ForesterUtil.isWritableFile( log_file ) );
        }
        EasyWriter log_writer = null;
        try {
            log_writer = ForesterUtil.createEasyWriter( log_file );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( gsdi.PRG_NAME, "Failed to create [" + log_file + "]: " + e.getMessage() );
        }
        Phylogeny species_tree = null;
        Phylogeny gene_tree = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            gene_tree = factory.create( gene_tree_file, PhyloXmlParser.createPhyloXmlParserXsdValidating() )[ 0 ];
        }
        catch ( final IOException e ) {
            fatalError( "error",
                        "failed to read gene tree from [" + gene_tree_file + "]: " + e.getMessage(),
                        log_writer );
        }
        try {
            species_tree = SDIutil.parseSpeciesTree( gene_tree,
                                                     species_tree_file,
                                                     REPLACE_UNDERSCORES_IN_NH_SPECIES_TREE,
                                                     true,
                                                     TAXONOMY_EXTRACTION.NO );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            fatalError( "user error",
                        "failed to transfer general node name, in [" + species_tree_file + "]: " + e.getMessage(),
                        log_writer );
        }
        catch ( final SDIException e ) {
            fatalError( "user error", e.getMessage(), log_writer );
        }
        catch ( final IOException e ) {
            fatalError( "error",
                        "Failed to read species tree from [" + species_tree_file + "]: " + e.getMessage(),
                        log_writer );
        }
        gene_tree.setRooted( true );
        species_tree.setRooted( true );
        if ( !gene_tree.isCompletelyBinary() ) {
            fatalError( "user error", "gene tree [" + gene_tree_file + "] is not completely binary", log_writer );
        }
        if ( base_algorithm == ALGORITHM.SDI ) {
            if ( !species_tree.isCompletelyBinary() ) {
                fatalError( "user error",
                            "species tree is not completely binary, use GSDI or GSDIR instead",
                            log_writer );
            }
        }
        log_writer.println( PRG_NAME + " - " + PRG_DESC );
        log_writer.println( "  version         : " + PRG_VERSION );
        log_writer.println( "  date            : " + PRG_DATE );
        log_writer.println( "  forester version: " + ForesterConstants.FORESTER_VERSION );
        log_writer.println();
        log_writer.println( "Start time                               : "
                + new SimpleDateFormat( "yyyyMMdd HH:mm:ss" ).format( new Date() ) );
        System.out.println( "Start time                               : "
                + new SimpleDateFormat( "yyyyMMdd HH:mm:ss" ).format( new Date() ) );
        log_writer.println( "Gene tree file                           : " + gene_tree_file.getCanonicalPath() );
        System.out.println( "Gene tree file                           : " + gene_tree_file.getCanonicalPath() );
        log_writer.println( "Gene tree name                           : "
                + ( ForesterUtil.isEmpty( gene_tree.getName() ) ? "" : gene_tree.getName() ) );
        System.out.println( "Gene tree name                           : "
                + ( ForesterUtil.isEmpty( gene_tree.getName() ) ? "" : gene_tree.getName() ) );
        log_writer.println( "Species tree file                        : " + species_tree_file.getCanonicalPath() );
        System.out.println( "Species tree file                        : " + species_tree_file.getCanonicalPath() );
        log_writer.println( "Species tree name                        : "
                + ( ForesterUtil.isEmpty( species_tree.getName() ) ? "" : gene_tree.getName() ) );
        System.out.println( "Species tree name                        : "
                + ( ForesterUtil.isEmpty( species_tree.getName() ) ? "" : gene_tree.getName() ) );
        System.out.println( "Transfer taxonomy                        : " + transfer_taxonomy );
        GSDII gsdii = null;
        final long start_time = new Date().getTime();
        try {
            if ( base_algorithm == ALGORITHM.GSDI ) {
                System.out.println( "Algorithm                                : GSDI" );
                log_writer.println( "Algorithm                                : GSDI" );
            }
            else if ( base_algorithm == ALGORITHM.GSDIR ) {
                System.out.println( "Algorithm                                : GSDIR" );
                log_writer.println( "Algorithm                                : GSDIR" );
            }
            System.out.println( "Use most parsimonous duplication model   : " + most_parsimonous_duplication_model );
            System.out.println( "Allow stripping of gene tree nodes       : " + allow_stripping_of_gene_tree );
            log_writer.println( "Use most parsimonous duplication model   : " + most_parsimonous_duplication_model );
            log_writer.println( "Allow stripping of gene tree nodes       : " + allow_stripping_of_gene_tree );
            log_writer.flush();
            if ( base_algorithm == ALGORITHM.GSDI ) {
                gsdii = new GSDI( gene_tree,
                                  species_tree,
                                  most_parsimonous_duplication_model,
                                  allow_stripping_of_gene_tree,
                                  true,
                                  transfer_taxonomy );
            }
            else if ( base_algorithm == ALGORITHM.GSDIR ) {
                gsdii = new GSDIR( gene_tree, species_tree, allow_stripping_of_gene_tree, true, transfer_taxonomy );
            }
        }
        catch ( final SDIException e ) {
            fatalError( "user error", e.getLocalizedMessage(), log_writer );
        }
        catch ( final IOException e ) {
            fatalError( "error", e.toString(), log_writer );
        }
        catch ( final OutOfMemoryError e ) {
            ForesterUtil.outOfMemoryError( e );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            fatalError( "unexpected error", e.toString(), log_writer );
        }
        System.out.println( "Running time (excluding I/O)             : " + ( new Date().getTime() - start_time )
                + "ms" );
        log_writer.println( "Running time (excluding I/O)             : " + ( new Date().getTime() - start_time )
                + "ms" );
        System.out.println( "Mapping based on                         : " + gsdii.getTaxCompBase() );
        log_writer.println( "Mapping based on                         : " + gsdii.getTaxCompBase() );
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            if ( base_algorithm == ALGORITHM.GSDIR ) {
                final Phylogeny gt = ( ( GSDIR ) gsdii ).getMinDuplicationsSumGeneTree();
                gt.setRerootable( false );
                writer.toPhyloXML( out_file, gt, 0 );
            }
            else {
                gene_tree.setRerootable( false );
                writer.toPhyloXML( out_file, gene_tree, 0 );
            }
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME,
                                     "Failed to write to [" + out_file.getCanonicalPath() + "]: " + e.getMessage() );
        }
        System.out.println( "Wrote resulting gene tree to             : " + out_file.getCanonicalPath() );
        log_writer.println( "Wrote resulting gene tree to             : " + out_file.getCanonicalPath() );
        final File species_tree_used_file = new File( ForesterUtil.removeSuffix( out_file.toString() )
                + SUFFIX_FOR_SPECIES_TREE_USED );
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( species_tree_used_file, species_tree, 0 );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME,
                                     "Failed to write to [" + species_tree_used_file.getCanonicalPath() + "]: "
                                             + e.getMessage() );
        }
        System.out.println( "Wrote (stripped) species tree to         : " + species_tree_used_file.getCanonicalPath() );
        log_writer.println( "Wrote (stripped) species tree to         : " + species_tree_used_file.getCanonicalPath() );
        if ( ( gsdii.getReMappedScientificNamesFromGeneTree() != null )
                && !gsdii.getReMappedScientificNamesFromGeneTree().isEmpty() ) {
            System.out.println( "Number of gene tree species remapped     : "
                    + gsdii.getReMappedScientificNamesFromGeneTree().size() );
            log_writer.println( "Number of gene tree species remapped     : "
                    + gsdii.getReMappedScientificNamesFromGeneTree().size() );
            writeToRemappedFile( out_file, gsdii.getReMappedScientificNamesFromGeneTree(), log_writer );
        }
        System.out.println( "Number of external nodes in gene tree    : " + gene_tree.getNumberOfExternalNodes() );
        log_writer.println( "Number of external nodes in gene tree    : " + gene_tree.getNumberOfExternalNodes() );
        System.out.println( "Number of external nodes in species tree : " + species_tree.getNumberOfExternalNodes() );
        log_writer.println( "Number of external nodes in species tree : " + species_tree.getNumberOfExternalNodes() );
        final int poly = PhylogenyMethods.countNumberOfPolytomies( species_tree );
        System.out.println( "Number of polytomies in species tree     : " + poly );
        log_writer.println( "Number of polytomies in species tree     : " + poly );
        System.out.println( "External nodes stripped from gene tree   : "
                + gsdii.getStrippedExternalGeneTreeNodes().size() );
        log_writer.println( "External nodes stripped from gene tree   : "
                + gsdii.getStrippedExternalGeneTreeNodes().size() );
        System.out
                .println( "External nodes stripped from species tree: " + gsdii.getStrippedSpeciesTreeNodes().size() );
        log_writer
                .println( "External nodes stripped from species tree: " + gsdii.getStrippedSpeciesTreeNodes().size() );
        System.out.println();
        System.out.println( "Number of speciations                    : " + gsdii.getSpeciationsSum() );
        log_writer.println( "Number of speciations                    : " + gsdii.getSpeciationsSum() );
        if ( ( base_algorithm == ALGORITHM.GSDIR ) ) {
            final GSDIR gsdir = ( GSDIR ) gsdii;
            System.out.println( "Minimal number of duplications           : " + gsdir.getMinDuplicationsSum() );
            log_writer.println( "Minimal number of duplications           : " + gsdir.getMinDuplicationsSum() );
        }
        else if ( ( base_algorithm == ALGORITHM.GSDI ) ) {
            final GSDI gsdi = ( GSDI ) gsdii;
            System.out.println( "Number of duplications                   : " + gsdi.getDuplicationsSum() );
            log_writer.println( "Number of duplications                   : " + gsdi.getDuplicationsSum() );
            if ( !most_parsimonous_duplication_model ) {
                final int u = gsdi.getSpeciationOrDuplicationEventsSum();
                System.out.println( "Number of potential duplications         : " + u );
                log_writer.println( "Number of potential duplications         : " + u );
            }
        }
        log_writer.println();
        printMappedNodesToLog( log_writer, gsdii );
        log_writer.println();
        printStrippedGeneTreeNodesToLog( log_writer, gsdii );
        System.out.println();
        System.out.println( "Wrote log to                             : " + log_file.getCanonicalPath() );
        System.out.println();
        log_writer.close();
    }

    private final static void fatalError( final String type, final String msg, final EasyWriter log_writer ) {
        try {
            log_writer.flush();
            log_writer.println();
            log_writer.print( type.toUpperCase() + ": " );
            log_writer.println( msg );
            log_writer.close();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        ForesterUtil.fatalError( gsdi.PRG_NAME, msg );
    }

    private final static void print_help() {
        System.out.println( "Usage: " + PRG_NAME
                + " [-options] <gene tree file, or gene trees in-directory> <species tree> <outfile, or out-directory>" );
        System.out.println();
        System.out.println( "Options:" );
        System.out.println( " -" + ALLOW_STRIPPING_OF_GENE_TREE_OPTION
                + "         : to allow stripping of gene tree nodes without a matching species" );
        System.out.println( " -" + MOST_PARSIMONIOUS_OPTION
                + "         : use most parimonious duplication model for GSDI: " );
        System.out.println( "              assign nodes as speciations which would otherwise be assiged" );
        System.out.println( "              as potential duplications due to polytomies in the species tree" );
        System.out.println( " -" + GUESS_FORMAT_OF_SPECIES_TREE
                + "         : to allow species tree in other formats than phyloXML (i.e. Newick, NHX, Nexus)" );
        System.out.println( " -" + GSDIR_OPTION
                + "         : to use GSDIR algorithm instead of GSDI algorithm (re-rooting)" );
        System.out.println( " -" + TRANSFER_TAXONOMY_OPTION
                + "         : to transfer taxonomic data from species tree to gene tree" );
        System.out.println( " -" + SUFFIX_FOR_DIR_OPTION
                + "=<suffix>: suffix for gene trees for analyzing entire directory of trees" );
        System.out.println();
        System.out.println();
        System.out.println( "Gene tree(s):" );
        System.out.println( " in phyloXM format, with taxonomy and sequence data in appropriate fields" );
        System.out.println();
        System.out.println( "Species tree:" );
        System.out.println( " in phyloXML format (unless option -" + GUESS_FORMAT_OF_SPECIES_TREE + " is used)" );
        System.out.println();
        System.out.println( "Examples: gsdi -" + ALLOW_STRIPPING_OF_GENE_TREE_OPTION
                + " gene_tree.xml tree_of_life.xml out.xml" );
        System.out.println( "          gsdi -" + ALLOW_STRIPPING_OF_GENE_TREE_OPTION + " -" + SUFFIX_FOR_DIR_OPTION
                + "=.xml" + " gene_tree_dir tree_of_life.xml out_dir" );
        System.out.println( "          gsdi -" + ALLOW_STRIPPING_OF_GENE_TREE_OPTION + " -" + MOST_PARSIMONIOUS_OPTION
                + " -" + GSDIR_OPTION + " -" + TRANSFER_TAXONOMY_OPTION + " -" + SUFFIX_FOR_DIR_OPTION + "=.xml"
                + " gene_tree_dir tree_of_life.xml out_dir" );
        System.out.println();
    }

    private final static void printMappedNodesToLog( final EasyWriter log_writer, final GSDII gsdi )
            throws IOException {
        final SortedSet<String> ss = new TreeSet<String>();
        for( final PhylogenyNode n : gsdi.getMappedExternalSpeciesTreeNodes() ) {
            ss.add( n.toString() );
        }
        log_writer.println( "The following " + ss.size() + " species were used: " );
        for( final String s : ss ) {
            log_writer.println( "  " + s );
        }
    }

    private final static void printStrippedGeneTreeNodesToLog( final EasyWriter log_writer, final GSDII gsdi )
            throws IOException {
        final SortedMap<String, Integer> sm = new TreeMap<String, Integer>();
        for( final PhylogenyNode n : gsdi.getStrippedExternalGeneTreeNodes() ) {
            final String s = n.toString();
            if ( sm.containsKey( s ) ) {
                sm.put( s, sm.get( s ) + 1 );
            }
            else {
                sm.put( s, 1 );
            }
        }
        log_writer.println( "The following " + sm.size() + " nodes were stripped from the gene tree: " );
        for( final String s : sm.keySet() ) {
            final int count = sm.get( s );
            if ( count == 1 ) {
                log_writer.println( "  " + s );
            }
            else {
                log_writer.println( "  " + s + " [" + count + "]" );
            }
        }
    }

    private final static void writeToRemappedFile( final File out_file,
                                                   final SortedSet<String> remapped,
                                                   final EasyWriter log_writer )
            throws IOException {
        final File file = new File( ForesterUtil.removeSuffix( out_file.toString() ) + REMAPPED_SUFFIX );
        final EasyWriter remapped_writer = ForesterUtil.createEasyWriter( file );
        for( final String s : remapped ) {
            remapped_writer.println( s );
        }
        remapped_writer.close();
        System.out.println( "Wrote remapped gene tree species to      : " + file.getCanonicalPath() );
        log_writer.println( "Wrote remapped gene tree species to      : " + file.getCanonicalPath() );
    }
}
