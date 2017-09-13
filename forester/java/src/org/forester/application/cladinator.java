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
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.forester.clade_analysis.AnalysisMulti;
import org.forester.clade_analysis.Prefix;
import org.forester.clade_analysis.ResultMulti;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.EasyWriter;
import org.forester.util.ForesterUtil;
import org.forester.util.UserException;

public final class cladinator {

    final static private String        PRG_NAME                             = "cladinator";
    final static private String        PRG_VERSION                          = "1.03";
    final static private String        PRG_DATE                             = "170913";
    final static private String        PRG_DESC                             = "clades within clades of annotated labels -- analysis of pplacer-type outputs";
    final static private String        E_MAIL                               = "phyloxml@gmail.com";
    final static private String        WWW                                  = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String        HELP_OPTION_1                        = "help";
    final static private String        HELP_OPTION_2                        = "h";
    final static private String        SEP_OPTION                           = "s";
    final static private String        QUERY_PATTERN_OPTION                 = "q";
    final static private String        SPECIFICS_CUTOFF_OPTION              = "c";
    final static private String        MAPPING_FILE_OPTION                  = "m";
    final static private String        EXTRA_PROCESSING_OPTION1             = "x";
    final static private String        EXTRA_PROCESSING1_SEP_OPTION         = "xs";
    final static private String        EXTRA_PROCESSING1_KEEP_EXTRA_OPTION  = "xk";
    final static private String        VERBOSE_OPTION                       = "v";
    final static private double        SPECIFICS_CUTOFF_DEFAULT             = 0.8;
    final static private String        SEP_DEFAULT                          = ".";
    final static private Pattern       QUERY_PATTERN_DEFAULT                = AnalysisMulti.DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE;
    final static private String        EXTRA_PROCESSING1_SEP_DEFAULT        = "|";
    final static private boolean       EXTRA_PROCESSING1_KEEP_EXTRA_DEFAULT = false;
    private final static DecimalFormat df                                   = new DecimalFormat( "0.0###" );

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
            if ( ( cla.getNumberOfNames() != 1 ) && ( cla.getNumberOfNames() != 2 ) ) {
                print_help();
                System.exit( -1 );
            }
            final List<String> allowed_options = new ArrayList<>();
            allowed_options.add( SEP_OPTION );
            allowed_options.add( QUERY_PATTERN_OPTION );
            allowed_options.add( SPECIFICS_CUTOFF_OPTION );
            allowed_options.add( MAPPING_FILE_OPTION );
            allowed_options.add( EXTRA_PROCESSING_OPTION1 );
            allowed_options.add( EXTRA_PROCESSING1_SEP_OPTION );
            allowed_options.add( EXTRA_PROCESSING1_KEEP_EXTRA_OPTION );
            allowed_options.add( VERBOSE_OPTION );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            double cutoff_specifics = SPECIFICS_CUTOFF_DEFAULT;
            if ( cla.isOptionSet( SPECIFICS_CUTOFF_OPTION ) ) {
                if ( cla.isOptionValueSet( SPECIFICS_CUTOFF_OPTION ) ) {
                    cutoff_specifics = cla.getOptionValueAsDouble( SPECIFICS_CUTOFF_OPTION );
                    if ( cutoff_specifics < 0 ) {
                        ForesterUtil.fatalError( PRG_NAME, "cutoff cannot be negative" );
                    }
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "no value for cutoff for specifics" );
                }
            }
            String separator = SEP_DEFAULT;
            if ( cla.isOptionSet( SEP_OPTION ) ) {
                if ( cla.isOptionValueSet( SEP_OPTION ) ) {
                    separator = cla.getOptionValue( SEP_OPTION );
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "no value for separator option" );
                }
            }
            Pattern compiled_query_str = null;
            if ( cla.isOptionSet( QUERY_PATTERN_OPTION ) ) {
                if ( cla.isOptionValueSet( QUERY_PATTERN_OPTION ) ) {
                    final String query_str = cla.getOptionValue( QUERY_PATTERN_OPTION );
                    try {
                        compiled_query_str = Pattern.compile( query_str );
                    }
                    catch ( final PatternSyntaxException e ) {
                        ForesterUtil.fatalError( PRG_NAME, "error in regular expression: " + e.getMessage() );
                    }
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "no value for query pattern option" );
                }
            }
            File mapping_file = null;
            if ( cla.isOptionSet( MAPPING_FILE_OPTION ) ) {
                if ( cla.isOptionValueSet( MAPPING_FILE_OPTION ) ) {
                    final String mapping_file_str = cla.getOptionValue( MAPPING_FILE_OPTION );
                    final String error = ForesterUtil.isReadableFile( mapping_file_str );
                    if ( !ForesterUtil.isEmpty( error ) ) {
                        ForesterUtil.fatalError( PRG_NAME, error );
                    }
                    mapping_file = new File( mapping_file_str );
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "no value for mapping file" );
                }
            }
            final Pattern pattern = ( compiled_query_str != null ) ? compiled_query_str : QUERY_PATTERN_DEFAULT;
            final File intreefile = cla.getFile( 0 );
            final String error_intreefile = ForesterUtil.isReadableFile( intreefile );
            if ( !ForesterUtil.isEmpty( error_intreefile ) ) {
                ForesterUtil.fatalError( PRG_NAME, error_intreefile );
            }
            final File outtablefile;
            if ( cla.getNumberOfNames() > 1 ) {
                outtablefile = cla.getFile( 1 );
                final String error_outtablefile = ForesterUtil.isWritableFile( outtablefile );
                if ( !ForesterUtil.isEmpty( error_outtablefile ) ) {
                    ForesterUtil.fatalError( PRG_NAME, error_outtablefile );
                }
            }
            else {
                outtablefile = null;
            }
            final BasicTable<String> t;
            final SortedMap<String, String> map;
            if ( mapping_file != null ) {
                t = BasicTableParser.parse( mapping_file, '\t' );
                if ( t.getNumberOfColumns() != 2 ) {
                    ForesterUtil.fatalError( PRG_NAME,
                                             "mapping file needs to have 2 tab-separated columns, not "
                                                     + t.getNumberOfColumns() );
                }
                map = t.getColumnsAsMap( 0, 1 );
            }
            else {
                t = null;
                map = null;
            }
            final boolean extra_processing1;
            if ( cla.isOptionSet( EXTRA_PROCESSING_OPTION1 ) ) {
                extra_processing1 = true;
            }
            else {
                extra_processing1 = false;
            }
            String extra_processing1_sep = EXTRA_PROCESSING1_SEP_DEFAULT;
            if ( cla.isOptionSet( EXTRA_PROCESSING1_SEP_OPTION ) ) {
                if ( !extra_processing1 ) {
                    ForesterUtil.fatalError( PRG_NAME,
                                             "extra processing is not enabled, cannot set -"
                                                     + EXTRA_PROCESSING1_SEP_OPTION + " option" );
                }
                if ( cla.isOptionValueSet( EXTRA_PROCESSING1_SEP_OPTION ) ) {
                    extra_processing1_sep = cla.getOptionValue( EXTRA_PROCESSING1_SEP_OPTION );
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "no value for extra processing separator" );
                }
            }
            if ( ( extra_processing1_sep != null ) && extra_processing1_sep.equals( separator ) ) {
                ForesterUtil.fatalError( PRG_NAME,
                                         "extra processing separator must not be the same the annotation-separator" );
            }
            boolean extra_processing1_keep = EXTRA_PROCESSING1_KEEP_EXTRA_DEFAULT;
            if ( cla.isOptionSet( EXTRA_PROCESSING1_KEEP_EXTRA_OPTION ) ) {
                if ( !extra_processing1 ) {
                    ForesterUtil.fatalError( PRG_NAME,
                                             "extra processing is not enabled, cannot set -"
                                                     + EXTRA_PROCESSING1_KEEP_EXTRA_OPTION + " option" );
                }
                extra_processing1_keep = true;
            }
            final boolean verbose;
            if ( cla.isOptionSet( VERBOSE_OPTION ) ) {
                verbose = true;
            }
            else {
                verbose = false;
            }
            System.out.println( "Input tree                 : " + intreefile );
            System.out.println( "Specific-hit support cutoff: " + cutoff_specifics );
            if ( mapping_file != null ) {
                System.out.println( "Mapping file               : " + mapping_file + " (" + t.getNumberOfRows()
                        + " rows)" );
            }
            System.out.println( "Annotation-separator       : " + separator );
            System.out.println( "Query pattern              : " + pattern );
            System.out.println( "Extra processing           : " + extra_processing1 );
            if ( extra_processing1 ) {
                System.out.println( "Extra processing separator : " + extra_processing1_sep );
                System.out.println( "Keep extra annotations     : " + extra_processing1_keep );
            }
            if ( outtablefile != null ) {
                System.out.println( "Output table               : " + outtablefile );
            }
            Phylogeny phys[] = null;
            try {
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intreefile, true );
                phys = factory.create( intreefile, pp );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, "Could not read \"" + intreefile + "\" [" + e.getMessage() + "]" );
            }
            if ( phys.length == 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "\"" + intreefile + "\" does not contain any trees" );
            }
            System.out.println( "Number of input trees      : " + phys.length );
            if ( phys.length == 1 ) {
                System.out.println( "Ext. nodes in input tree 1 : " + phys[ 0 ].getNumberOfExternalNodes() );
            }
            else {
                System.out.println( "Ext. nodes in input tree   : " + phys[ 0 ].getNumberOfExternalNodes() );
            }
            final EasyWriter outtable_writer;
            if ( outtablefile != null ) {
                outtable_writer = ForesterUtil.createEasyWriter( outtablefile );
            }
            else {
                outtable_writer = null;
            }
            for( final Phylogeny phy : phys ) {
                if ( map != null ) {
                    AnalysisMulti.performMapping( pattern, map, phy, verbose );
                }
                if ( extra_processing1 ) {
                    AnalysisMulti.performExtraProcessing1( pattern,
                                                           phy,
                                                           extra_processing1_sep,
                                                           extra_processing1_keep,
                                                           separator,
                                                           verbose );
                }
                final ResultMulti res = AnalysisMulti.execute( phy, pattern, separator, cutoff_specifics );
                printResult( res );
                if ( outtable_writer != null ) {
                    writeResultToTable( res, outtable_writer );
                    outtable_writer.flush();
                }
            }
            if ( outtable_writer != null ) {
                outtable_writer.flush();
                outtable_writer.close();
            }
        }
        catch ( final UserException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ForesterUtil.fatalError( PRG_NAME, "Unexpected errror!" );
        }
    }

    private final static void printResult( final ResultMulti res ) {
        System.out.println();
        System.out.println( "Result for " + res.getQueryNamePrefix() );
        System.out.println();
        if ( ( res.getAllMultiHitPrefixes() == null ) | ( res.getAllMultiHitPrefixes().size() < 1 ) ) {
            System.out.println( " No match to query pattern!" );
        }
        else {
            System.out.println( " Matching Clade(s):" );
            for( final Prefix prefix : res.getCollapsedMultiHitPrefixes() ) {
                System.out.println( " " + prefix );
            }
            if ( res.isHasSpecificMultiHitsPrefixes() ) {
                System.out.println();
                System.out.println( " Specific-hit(s):" );
                for( final Prefix prefix : res.getSpecificMultiHitPrefixes() ) {
                    System.out.println( " " + prefix );
                }
                System.out.println();
                System.out.println( " Matching Clade(s) with Specific-hit(s):" );
                for( final Prefix prefix : res.getCollapsedMultiHitPrefixes() ) {
                    System.out.println( " " + prefix );
                    for( final Prefix spec : res.getSpecificMultiHitPrefixes() ) {
                        if ( spec.getPrefix().startsWith( prefix.getPrefix() ) ) {
                            System.out.println( "     " + spec );
                        }
                    }
                }
            }
            if ( !ForesterUtil.isEmpty( res.getAllMultiHitPrefixesDown() ) ) {
                System.out.println();
                System.out.println( " Matching Down-tree Bracketing Clade(s):" );
                for( final Prefix prefix : res.getCollapsedMultiHitPrefixesDown() ) {
                    System.out.println( " " + prefix );
                }
            }
            if ( !ForesterUtil.isEmpty( res.getAllMultiHitPrefixesUp() ) ) {
                System.out.println();
                System.out.println( " Matching Up-tree Bracketing Clade(s):" );
                for( final Prefix prefix : res.getCollapsedMultiHitPrefixesUp() ) {
                    System.out.println( " " + prefix );
                }
            }
        }
        System.out.println();
    }

    private final static void writeResultToTable( final ResultMulti res, final EasyWriter w ) throws IOException {
        if ( ( res.getAllMultiHitPrefixes() == null ) | ( res.getAllMultiHitPrefixes().size() < 1 ) ) {
            w.print( res.getQueryNamePrefix() );
            w.print( "\t" );
            w.println( "No match to query pattern!" );
        }
        else {
            for( final Prefix prefix : res.getCollapsedMultiHitPrefixes() ) {
                w.print( res.getQueryNamePrefix() );
                w.print( "\t" );
                w.print( "Matching Clades" );
                w.print( "\t" );
                w.print( prefix.getPrefix() );
                w.print( "\t" );
                w.print( df.format( prefix.getConfidence() ) );
                w.println();
            }
            if ( res.isHasSpecificMultiHitsPrefixes() ) {
                for( final Prefix prefix : res.getSpecificMultiHitPrefixes() ) {
                    w.print( res.getQueryNamePrefix() );
                    w.print( "\t" );
                    w.print( "Specific-hits" );
                    w.print( "\t" );
                    w.print( prefix.getPrefix() );
                    w.print( "\t" );
                    w.print( df.format( prefix.getConfidence() ) );
                    w.println();
                }
            }
            if ( !ForesterUtil.isEmpty( res.getAllMultiHitPrefixesDown() ) ) {
                for( final Prefix prefix : res.getCollapsedMultiHitPrefixesDown() ) {
                    w.print( res.getQueryNamePrefix() );
                    w.print( "\t" );
                    w.print( "Matching Down-tree Bracketing Clades" );
                    w.print( "\t" );
                    w.print( prefix.getPrefix() );
                    w.print( "\t" );
                    w.print( df.format( prefix.getConfidence() ) );
                    w.println();
                }
            }
            if ( !ForesterUtil.isEmpty( res.getAllMultiHitPrefixesUp() ) ) {
                for( final Prefix prefix : res.getCollapsedMultiHitPrefixesUp() ) {
                    w.print( res.getQueryNamePrefix() );
                    w.print( "\t" );
                    w.print( "Matching Up-tree Bracketing Clades" );
                    w.print( "\t" );
                    w.print( prefix.getPrefix() );
                    w.print( "\t" );
                    w.print( df.format( prefix.getConfidence() ) );
                    w.println();
                }
            }
        }
    }

    private final static void print_help() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " [options] <input tree(s) file> [output table file]" );
        System.out.println();
        System.out.println( " options:" );
        System.out.println( "  -" + SPECIFICS_CUTOFF_OPTION
                + "=<double>        : the cutoff for \"specific-hit\" support values (default: "
                + SPECIFICS_CUTOFF_DEFAULT + ")" );
        System.out.println( "  -" + SEP_OPTION + "=<separator>     : the annotation-separator to be used (default: "
                + SEP_DEFAULT + ")" );
        System.out.println( "  -" + MAPPING_FILE_OPTION
                + "=<mapping table> : to map node names to appropriate annotations (tab-separated, two columns) (default: no mapping)" );
        System.out.println( "  -" + EXTRA_PROCESSING_OPTION1
                + "                 : to enable extra processing of annotations (e.g. \"Q16611|A.1.1\" becomes \"A.1.1\")" );
        System.out.println( "  -" + EXTRA_PROCESSING1_SEP_OPTION
                + "=<separator>    : the separator for extra annotations (default: \"" + EXTRA_PROCESSING1_SEP_DEFAULT
                + "\")" );
        System.out.println( "  -" + EXTRA_PROCESSING1_KEEP_EXTRA_OPTION
                + "                : to keep extra annotations (e.g. \"Q16611|A.1.1\" becomes \"A.1.1.Q16611\")" );
        System.out.println( "  -" + VERBOSE_OPTION + "                 : verbose" );
        System.out.println( "  --" + QUERY_PATTERN_OPTION
                + "=<query pattern>: the regular expression for the query (default: \"" + QUERY_PATTERN_DEFAULT
                + "\" for pplacer output)" );
        System.out.println();
        System.out.println( "Examples:" );
        System.out.println();
        System.out.println( " " + PRG_NAME + " my_tree.nh result.tsv" );
        System.out.println( " " + PRG_NAME + " -c=0.5 -s=. my_tree.nh result.tsv" );
        System.out.println( " " + PRG_NAME + " -c=0.9 -s=_ -m=map.tsv my_tree.nh result.tsv" );
        System.out.println( " " + PRG_NAME + " -x -xs=& -xk my_tree.nh result.tsv" );
        System.out.println( " " + PRG_NAME + " -x -xs=\"|\" my_tree.nh result.tsv" );
        System.out.println( " " + PRG_NAME + " -x -xk -m=map.tsv pplacer_out_trees.sing.tre result.tsv" );
        System.out.println();
    }
}
