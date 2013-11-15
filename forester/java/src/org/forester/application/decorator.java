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
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sequence.Sequence;
import org.forester.tools.PhylogenyDecorator;
import org.forester.tools.PhylogenyDecorator.FIELD;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class decorator {

    private static final String SEQUENCE_NAME_FIELD                     = "s";
    private static final String MOL_SEQ                                 = "m";
    private static final String TAXONOMY_CODE_FIELD                     = "c";
    private static final String TAXONOMY_SCIENTIFIC_NAME_FIELD          = "sn";
    private static final String DS_FILED                                = "d";
    private static final String SEQUENCE_ANNOTATION_DESC                = "a";
    private static final String NODE_NAME_FIELD                         = "n";
    final static private String PICKY_OPTION                            = "p";
    final static private String FIELD_OPTION                            = "f";
    final static private String TRIM_AFTER_TILDE_OPTION                 = "t";
    final static private String VERBOSE_OPTION                          = "ve";
    final static private String TREE_NAME_OPTION                        = "pn";
    final static private String TREE_ID_OPTION                          = "pi";
    final static private String TREE_DESC_OPTION                        = "pd";
    final static private String MIDPOINT_ROOT_OPTION                    = "mp";
    final static private String ORDER_TREE_OPTION                       = "or";
    final static private String EXTRACT_BRACKETED_SCIENTIC_NAME_OPTION  = "sn";
    final static private String EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION = "tc";
    final static private String CUT_NAME_AFTER_FIRST_SPACE_OPTION       = "c";
    final static private String ADVANCED_TABLE_OPTION                   = "table";
    final static private String KEY_COLUMN                              = "k";
    final static private String VALUE_COLUMN                            = "v";
    final static private String MAPPING_FILE_SEPARATOR_OPTION           = "s";
    final static private char   MAPPING_FILE_SEPARATOR_DEFAULT          = '\t';
    final static private String PRG_NAME                                = "decorator";
    final static private String PRG_VERSION                             = "1.16";
    final static private String PRG_DATE                                = "131113";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( decorator.PRG_NAME, decorator.PRG_VERSION, decorator.PRG_DATE );
        System.out.println();
        if ( ( args.length < 4 ) || ( args.length > 13 ) ) {
            decorator.argumentsError();
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( ( cla.getNumberOfNames() < 3 ) || ( cla.getNumberOfNames() > 4 ) ) {
            decorator.argumentsError();
        }
        final File phylogenies_infile = cla.getFile( 0 );
        final File mapping_infile = cla.getFile( 1 );
        final File phylogenies_outfile = cla.getFile( 2 );
        if ( phylogenies_outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + phylogenies_outfile + "] already exists" );
        }
        String err = ForesterUtil.isReadableFile( phylogenies_infile );
        if ( !ForesterUtil.isEmpty( err ) ) {
            ForesterUtil.fatalError( PRG_NAME, err );
        }
        err = ForesterUtil.isReadableFile( mapping_infile );
        if ( !ForesterUtil.isEmpty( err ) ) {
            ForesterUtil.fatalError( PRG_NAME, err );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( decorator.ADVANCED_TABLE_OPTION );
        allowed_options.add( decorator.PICKY_OPTION );
        allowed_options.add( decorator.FIELD_OPTION );
        allowed_options.add( decorator.CUT_NAME_AFTER_FIRST_SPACE_OPTION );
        allowed_options.add( decorator.KEY_COLUMN );
        allowed_options.add( decorator.VALUE_COLUMN );
        allowed_options.add( decorator.MAPPING_FILE_SEPARATOR_OPTION );
        allowed_options.add( decorator.EXTRACT_BRACKETED_SCIENTIC_NAME_OPTION );
        allowed_options.add( decorator.EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION );
        allowed_options.add( decorator.TREE_NAME_OPTION );
        allowed_options.add( decorator.TREE_ID_OPTION );
        allowed_options.add( decorator.TREE_DESC_OPTION );
        allowed_options.add( decorator.TRIM_AFTER_TILDE_OPTION );
        allowed_options.add( decorator.ORDER_TREE_OPTION );
        allowed_options.add( decorator.MIDPOINT_ROOT_OPTION );
        allowed_options.add( decorator.VERBOSE_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final boolean advanced_table = cla.isOptionSet( decorator.ADVANCED_TABLE_OPTION );
        if ( !advanced_table ) {
            final List<String> mandatory_options = new ArrayList<String>();
            mandatory_options.add( decorator.FIELD_OPTION );
            final String missing_options = cla.validateMandatoryOptionsAsString( mandatory_options );
            if ( missing_options.length() > 0 ) {
                ForesterUtil.fatalError( decorator.PRG_NAME, "missing option(s): " + missing_options );
            }
        }
        final boolean picky = cla.isOptionSet( decorator.PICKY_OPTION );
        char separator = decorator.MAPPING_FILE_SEPARATOR_DEFAULT;
        if ( cla.isOptionSet( decorator.MAPPING_FILE_SEPARATOR_OPTION ) ) {
            if ( advanced_table ) {
                argumentsError();
            }
            separator = cla.getOptionValueAsChar( decorator.MAPPING_FILE_SEPARATOR_OPTION );
        }
        int key_column = 0;
        int value_column = 1;
        String field_str = "";
        FIELD field = FIELD.NODE_NAME;
        boolean cut_name_after_space = false;
        boolean extract_bracketed_scientific_name = false;
        boolean extract_bracketed_tax_code = false;
        boolean trim_after_tilde = false;
        boolean order_tree = false;
        boolean midpoint_root = false;
        boolean verbose = false;
        String tree_name = "";
        String tree_id = "";
        String tree_desc = "";
        try {
            if ( cla.isOptionSet( decorator.TREE_NAME_OPTION ) ) {
                tree_name = cla.getOptionValueAsCleanString( decorator.TREE_NAME_OPTION );
            }
            if ( cla.isOptionSet( decorator.TREE_ID_OPTION ) ) {
                tree_id = cla.getOptionValueAsCleanString( decorator.TREE_ID_OPTION );
            }
            if ( cla.isOptionSet( decorator.TREE_DESC_OPTION ) ) {
                tree_desc = cla.getOptionValueAsCleanString( decorator.TREE_DESC_OPTION );
            }
            if ( cla.isOptionSet( decorator.EXTRACT_BRACKETED_SCIENTIC_NAME_OPTION ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                extract_bracketed_scientific_name = true;
            }
            if ( cla.isOptionSet( decorator.EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                extract_bracketed_tax_code = true;
            }
            if ( cla.isOptionSet( decorator.KEY_COLUMN ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                key_column = cla.getOptionValueAsInt( decorator.KEY_COLUMN );
            }
            if ( cla.isOptionSet( decorator.VALUE_COLUMN ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                value_column = cla.getOptionValueAsInt( decorator.VALUE_COLUMN );
            }
            if ( cla.isOptionSet( decorator.CUT_NAME_AFTER_FIRST_SPACE_OPTION ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                cut_name_after_space = true;
            }
            if ( cla.isOptionSet( decorator.TRIM_AFTER_TILDE_OPTION ) ) {
                if ( advanced_table ) {
                    argumentsError();
                }
                trim_after_tilde = true;
            }
            if ( cla.isOptionSet( decorator.MIDPOINT_ROOT_OPTION ) ) {
                midpoint_root = true;
            }
            if ( cla.isOptionSet( decorator.ORDER_TREE_OPTION ) ) {
                order_tree = true;
            }
            if ( cla.isOptionSet( decorator.VERBOSE_OPTION ) ) {
                verbose = true;
            }
            if ( cla.isOptionSet( decorator.FIELD_OPTION ) ) {
                field_str = cla.getOptionValue( decorator.FIELD_OPTION );
                if ( field_str.equals( NODE_NAME_FIELD ) ) {
                    field = FIELD.NODE_NAME;
                }
                else if ( field_str.equals( SEQUENCE_ANNOTATION_DESC ) ) {
                    field = FIELD.SEQUENCE_ANNOTATION_DESC;
                }
                else if ( field_str.equals( DS_FILED ) ) {
                    field = FIELD.DOMAIN_STRUCTURE;
                    extract_bracketed_scientific_name = false;
                    extract_bracketed_tax_code = false;
                }
                else if ( field_str.equals( TAXONOMY_CODE_FIELD ) ) {
                    field = FIELD.TAXONOMY_CODE;
                }
                else if ( field_str.equals( SEQUENCE_NAME_FIELD ) ) {
                    field = FIELD.SEQUENCE_NAME;
                }
                else if ( field_str.equals( MOL_SEQ ) ) {
                    field = FIELD.MOL_SEQ;
                }
                else if ( field_str.equals( TAXONOMY_SCIENTIFIC_NAME_FIELD ) ) {
                    field = FIELD.TAXONOMY_SCIENTIFIC_NAME;
                    extract_bracketed_scientific_name = false;
                    extract_bracketed_tax_code = false;
                }
                else {
                    ForesterUtil.fatalError( decorator.PRG_NAME, "unknown value for \"" + decorator.FIELD_OPTION
                            + "\" option: \"" + field_str + "\"" );
                }
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "error in command line: " + e.getMessage() );
        }
        if ( extract_bracketed_scientific_name && extract_bracketed_tax_code ) {
            argumentsError();
        }
        ForesterUtil.programMessage( PRG_NAME, "input tree(s) : " + phylogenies_infile );
        ForesterUtil.programMessage( PRG_NAME, "map           : " + mapping_infile );
        ForesterUtil.programMessage( PRG_NAME, "output tree(s): " + phylogenies_outfile );
        System.out.println();
        Phylogeny[] phylogenies = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( phylogenies_infile, true );
            phylogenies = factory.create( phylogenies_infile, pp );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "failed to read phylgenies from [" + phylogenies_infile
                    + "] [" + e.getMessage() + "]" );
        }
        Map<String, String> map = null;
        if ( !advanced_table ) {
            if ( field != FIELD.MOL_SEQ ) {
                BasicTable<String> mapping_table = null;
                try {
                    mapping_table = BasicTableParser.parse( mapping_infile, separator, true, false );
                }
                catch ( final Exception e ) {
                    ForesterUtil.fatalError( decorator.PRG_NAME,
                                             "failed to read [" + mapping_infile + "] [" + e.getMessage() + "]" );
                }
                if ( ( key_column < 0 ) || ( key_column >= mapping_table.getNumberOfColumns() ) ) {
                    ForesterUtil.fatalError( decorator.PRG_NAME, "illegal value for key column" );
                }
                if ( ( value_column < 0 ) || ( value_column >= mapping_table.getNumberOfColumns() ) ) {
                    ForesterUtil.fatalError( decorator.PRG_NAME, "illegal value for value column" );
                }
                if ( mapping_table.isEmpty() || ( mapping_table.getNumberOfColumns() < 1 ) ) {
                    ForesterUtil.fatalError( decorator.PRG_NAME, "mapping table is empty" );
                }
                if ( mapping_table.getNumberOfColumns() == 1 ) {
                    ForesterUtil.fatalError( decorator.PRG_NAME, "mapping table has only one column" );
                }
                map = mapping_table.getColumnsAsMap( key_column, value_column );
                final Iterator<Entry<String, String>> iter = map.entrySet().iterator();
                if ( verbose ) {
                    System.out.println();
                }
                while ( iter.hasNext() ) {
                    final Entry<String, String> e = iter.next();
                    if ( ForesterUtil.isEmpty( e.getKey() ) ) {
                        ForesterUtil.fatalError( decorator.PRG_NAME, "mapping table contains empty key" );
                    }
                    if ( ForesterUtil.isEmpty( e.getValue() ) ) {
                        ForesterUtil.fatalError( decorator.PRG_NAME, "mapping table contains empty value" );
                    }
                    if ( verbose ) {
                        System.out.println( e.getKey() + " => " + e.getValue() );
                    }
                }
                if ( verbose ) {
                    System.out.println();
                }
            }
            else {
                map = readFastaFileIntoMap( mapping_infile, verbose );
            }
        }
        if ( !ForesterUtil.isEmpty( tree_name ) || !ForesterUtil.isEmpty( tree_id )
                || !ForesterUtil.isEmpty( tree_desc ) ) {
            if ( ( phylogenies.length > 1 )
                    && ( !ForesterUtil.isEmpty( tree_name ) || !ForesterUtil.isEmpty( tree_id ) ) ) {
                ForesterUtil.fatalError( decorator.PRG_NAME,
                                         "attempt to set same name or id on more than one phylogeny" );
            }
            if ( !ForesterUtil.isEmpty( tree_name ) ) {
                phylogenies[ 0 ].setName( tree_name );
            }
            if ( !ForesterUtil.isEmpty( tree_id ) ) {
                final String[] s_ary = tree_id.split( ":" );
                phylogenies[ 0 ].setIdentifier( new Identifier( s_ary[ 1 ], s_ary[ 0 ] ) );
            }
            if ( !ForesterUtil.isEmpty( tree_desc ) ) {
                for( final Phylogeny phylogenie : phylogenies ) {
                    phylogenie.setDescription( tree_desc );
                }
            }
        }
        try {
            if ( advanced_table ) {
                Map<String, Map<String, String>> table = null;
                try {
                    table = PhylogenyDecorator.parseMappingTable( mapping_infile );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( decorator.PRG_NAME,
                                             "failed to read \"" + mapping_infile + "\" [" + e.getMessage() + "]" );
                }
                for( final Phylogeny phylogenie : phylogenies ) {
                    PhylogenyDecorator.decorate( phylogenie, table, picky );
                }
            }
            else {
                for( final Phylogeny phylogenie : phylogenies ) {
                    final String msg = PhylogenyDecorator.decorate( phylogenie,
                                                                    map,
                                                                    field,
                                                                    extract_bracketed_scientific_name,
                                                                    extract_bracketed_tax_code,
                                                                    picky,
                                                                    cut_name_after_space,
                                                                    trim_after_tilde,
                                                                    verbose );
                    ForesterUtil.programMessage( PRG_NAME, msg );
                }
            }
        }
        catch ( final NullPointerException e ) {
            ForesterUtil.unexpectedFatalError( decorator.PRG_NAME, e );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, e.getLocalizedMessage() );
        }
        if ( midpoint_root || order_tree ) {
            for( final Phylogeny phy : phylogenies ) {
                if ( midpoint_root ) {
                    PhylogenyMethods.midpointRoot( phy );
                }
                if ( order_tree ) {
                    PhylogenyMethods.orderAppearance( phy.getRoot(), true, true, DESCENDANT_SORT_PRIORITY.TAXONOMY );
                }
            }
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( phylogenies, 0, phylogenies_outfile, ForesterUtil.getLineSeparator() );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "failed to write output [" + e.getMessage() + "]" );
        }
        System.out.println();
        ForesterUtil.programMessage( PRG_NAME, "wrote: " + phylogenies_outfile );
        ForesterUtil.programMessage( PRG_NAME, "OK." );
    }

    private static Map<String, String> readFastaFileIntoMap( final File mapping_infile, final boolean verbose ) {
        List<Sequence> seqs = null;
        try {
            seqs = FastaParser.parse( new FileInputStream( mapping_infile ) );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "failed to read fasta-file from [" + mapping_infile + "] ["
                    + e.getMessage() + "]" );
        }
        if ( ForesterUtil.isEmpty( seqs ) ) {
            ForesterUtil.fatalError( decorator.PRG_NAME, "fasta-file [" + mapping_infile
                    + "] is devoid of fasta-formatted sequences" );
        }
        final Map<String, String> map = new HashMap<String, String>();
        for( final Sequence seq : seqs ) {
            if ( ForesterUtil.isEmpty( seq.getIdentifier() ) ) {
                ForesterUtil.fatalError( decorator.PRG_NAME, "fasta-file [" + mapping_infile
                        + "] contains sequence with empty identifier" );
            }
            if ( map.containsKey( seq.getIdentifier() ) ) {
                ForesterUtil.fatalError( decorator.PRG_NAME, "sequence identifier [" + seq.getIdentifier()
                        + "] is not unique" );
            }
            if ( seq.getLength() < 1 ) {
                ForesterUtil.fatalError( decorator.PRG_NAME, "sequence [" + seq.getIdentifier() + "] is empty" );
            }
            map.put( seq.getIdentifier(), seq.getMolecularSequenceAsString() );
            if ( verbose ) {
                System.out.println( seq.getIdentifier() + " => " + seq.getMolecularSequenceAsString() );
            }
        }
        return map;
    }

    private static void argumentsError() {
        System.out.println();
        System.out.println( decorator.PRG_NAME + " -" + ADVANCED_TABLE_OPTION + " | -f=<c> <phylogenies infile> "
                + "<mapping table file|fasta-file> <phylogenies outfile>" );
        System.out.println();
        System.out.println( "options:" );
        System.out.println();
        System.out.println( " -" + ADVANCED_TABLE_OPTION + " : table instead of one to one map (-f=<c>)" );
        System.out.println( " -p     : picky, fails if node name not found in mapping table" );
        System.out.println( " -" + TREE_NAME_OPTION + "=<s>: name for the phylogeny" );
        System.out.println( " -" + TREE_ID_OPTION + "=<s>: identifier for the phylogeny (in the form provider:value)" );
        System.out.println( " -" + TREE_DESC_OPTION + "=<s>: description for phylogenies" );
        System.out.println();
        System.out.println();
        System.out.println( "advanced options, only available if -" + ADVANCED_TABLE_OPTION + " is not used:" );
        System.out.println();
        System.out.println( " -f=<c> : field to be replaced: " + NODE_NAME_FIELD + " : node name" );
        System.out.println( "                                " + SEQUENCE_ANNOTATION_DESC
                + " : sequence annotation description" );
        System.out.println( "                                " + DS_FILED + " : domain structure" );
        System.out.println( "                                " + TAXONOMY_CODE_FIELD + " : taxonomy code" );
        System.out.println( "                                " + TAXONOMY_SCIENTIFIC_NAME_FIELD
                + ": taxonomy scientific name" );
        System.out.println( "                                " + SEQUENCE_NAME_FIELD + " : sequence name" );
        System.out.println( "                                " + MOL_SEQ + " : molecular sequence" );
        System.out.println( " -k=<n> : key column in mapping table (0 based)," );
        System.out.println( "          names of the node to be decorated - default is 0" );
        System.out.println( " -v=<n> : value column in mapping table (0 based)," );
        System.out.println( "          data which with to decorate - default is 1" );
        System.out.println( " -" + EXTRACT_BRACKETED_SCIENTIC_NAME_OPTION
                + "    : to extract bracketed scientific names, e.g. [Nematostella vectensis]" );
        System.out.println( " -" + EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION
                + "    : to extract bracketed taxonomic codes, e.g. [NEMVE]" );
        System.out.println( " -s=<c> : column separator in mapping file, default is tab" );
        System.out.println( " -c     : cut name after first space (only for -f=n)" );
        System.out.println( " -" + decorator.TRIM_AFTER_TILDE_OPTION
                + "     : trim node name to be replaced after tilde" );
        System.out.println( " -" + decorator.MIDPOINT_ROOT_OPTION + "    : to midpoint-root the tree" );
        System.out.println( " -" + decorator.ORDER_TREE_OPTION + "    : to order tree branches" );
        System.out.println( " -" + decorator.VERBOSE_OPTION + "    : verbose" );
        System.out.println();
        System.exit( -1 );
    }
}
