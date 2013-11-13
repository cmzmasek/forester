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

package org.forester.tools;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;

import org.forester.io.parsers.nhx.NHXFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

public final class PhylogenyDecorator {

    public final static boolean SANITIZE                = false;
    final private static String TP_NODE_NAME            = "NODE_NAME";
    final private static String TP_SEQ_ACCESSION        = "SEQ_ACCESSION";
    final private static String TP_SEQ_ACCESSION_SOURCE = "SEQ_ACCESSION_SOURCE";
    final private static String TP_SEQ_ANNOTATION_DESC  = "SEQ_ANNOTATION_DESC";
    final private static String TP_SEQ_ANNOTATION_REF   = "SEQ_ANNOTATION_REF";
    final private static String TP_SEQ_MOL_SEQ          = "SEQ_MOL_SEQ";
    final private static String TP_SEQ_NAME             = "SEQ_NAME";
    final private static String TP_SEQ_SYMBOL           = "SEQ_SYMBOL";
    final private static String TP_TAXONOMY_CN          = "TAXONOMY_CN";
    // From evoruby/lib/evo/apps/tseq_taxonomy_processor.rb:
    final private static String TP_TAXONOMY_CODE        = "TAXONOMY_CODE";
    final private static String TP_TAXONOMY_ID          = "TAXONOMY_ID";
    final private static String TP_TAXONOMY_ID_PROVIDER = "TAXONOMY_ID_PROVIDER";
    final private static String TP_TAXONOMY_SN          = "TAXONOMY_SN";
    final private static String TP_TAXONOMY_SYN         = "TAXONOMY_SYN";

    private PhylogenyDecorator() {
        // Not needed.
    }

    public static void decorate( final Phylogeny phylogeny,
                                 final Map<String, Map<String, String>> map,
                                 final boolean picky,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map )
            throws IllegalArgumentException, PhyloXmlDataFormatException {
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final String name = node.getName();
            if ( !ForesterUtil.isEmpty( name ) ) {
                if ( map.containsKey( name ) || ( numbers_of_chars_allowed_to_remove_if_not_found_in_map > 0 ) ) {
                    Map<String, String> new_values = map.get( name );
                    int x = 0;
                    while ( ( new_values == null ) && ( numbers_of_chars_allowed_to_remove_if_not_found_in_map > 0 )
                            && ( x <= numbers_of_chars_allowed_to_remove_if_not_found_in_map ) ) {
                        new_values = map.get( name.substring( 0, name.length() - x ) );
                        ++x;
                    }
                    if ( new_values != null ) {
                        if ( new_values.containsKey( TP_TAXONOMY_CODE ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().setTaxonomyCode( new_values.get( TP_TAXONOMY_CODE ) );
                        }
                        if ( new_values.containsKey( TP_TAXONOMY_ID )
                                && new_values.containsKey( TP_TAXONOMY_ID_PROVIDER ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData()
                                    .getTaxonomy()
                                    .setIdentifier( new Identifier( new_values.get( TP_TAXONOMY_ID ),
                                                                    new_values.get( TP_TAXONOMY_ID_PROVIDER ) ) );
                        }
                        else if ( new_values.containsKey( TP_TAXONOMY_ID ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy()
                                    .setIdentifier( new Identifier( new_values.get( TP_TAXONOMY_ID ) ) );
                        }
                        if ( new_values.containsKey( TP_TAXONOMY_SN ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().setScientificName( new_values.get( TP_TAXONOMY_SN ) );
                        }
                        if ( new_values.containsKey( TP_TAXONOMY_CN ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().setCommonName( new_values.get( TP_TAXONOMY_CN ) );
                        }
                        if ( new_values.containsKey( TP_TAXONOMY_SYN ) ) {
                            ForesterUtil.ensurePresenceOfTaxonomy( node );
                            node.getNodeData().getTaxonomy().getSynonyms().add( new_values.get( TP_TAXONOMY_SYN ) );
                        }
                        if ( new_values.containsKey( TP_SEQ_ACCESSION )
                                && new_values.containsKey( TP_SEQ_ACCESSION_SOURCE ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            node.getNodeData()
                                    .getSequence()
                                    .setAccession( new Accession( new_values.get( TP_SEQ_ACCESSION ),
                                                                  new_values.get( TP_SEQ_ACCESSION_SOURCE ) ) );
                        }
                        if ( new_values.containsKey( TP_SEQ_ANNOTATION_DESC ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            final Annotation ann = new Annotation();
                            ann.setDesc( new_values.get( TP_SEQ_ANNOTATION_DESC ) );
                            node.getNodeData().getSequence().addAnnotation( ann );
                        }
                        if ( new_values.containsKey( TP_SEQ_ANNOTATION_REF ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            final Annotation ann = new Annotation( new_values.get( TP_SEQ_ANNOTATION_REF ) );
                            node.getNodeData().getSequence().addAnnotation( ann );
                        }
                        if ( new_values.containsKey( TP_SEQ_SYMBOL ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            node.getNodeData().getSequence().setSymbol( new_values.get( TP_SEQ_SYMBOL ) );
                        }
                        if ( new_values.containsKey( TP_SEQ_NAME ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            node.getNodeData().getSequence().setName( new_values.get( TP_SEQ_NAME ) );
                        }
                        if ( new_values.containsKey( TP_SEQ_MOL_SEQ ) ) {
                            ForesterUtil.ensurePresenceOfSequence( node );
                            node.getNodeData().getSequence().setMolecularSequence( new_values.get( TP_SEQ_MOL_SEQ ) );
                        }
                        if ( new_values.containsKey( TP_NODE_NAME ) ) {
                            node.setName( new_values.get( TP_NODE_NAME ) );
                        }
                    } // if ( new_values != null ) 
                } // if ( map.containsKey( name ) || ( numbers_of_chars_allowed_to_remove_if_not_found_in_map > 0 ) )
                else if ( picky ) {
                    throw new IllegalArgumentException( "\"" + name + "\" not found in name map" );
                }
            }
        }
    }

    public static void decorate( final Phylogeny phylogeny,
                                 final Map<String, String> map,
                                 final FIELD field,
                                 final boolean extract_bracketed_scientific_name,
                                 final boolean extract_bracketed_tax_code,
                                 final boolean picky,
                                 final boolean cut_name_after_space,
                                 final boolean process_name_intelligently,
                                 final boolean process_similar_to,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                 final boolean trim_after_tilde,
                                 final boolean verbose ) throws IllegalArgumentException, NHXFormatException,
            PhyloXmlDataFormatException {
        PhylogenyDecorator.decorate( phylogeny,
                                     map,
                                     field,
                                     extract_bracketed_scientific_name,
                                     extract_bracketed_tax_code,
                                     picky,
                                     null,
                                     cut_name_after_space,
                                     process_name_intelligently,
                                     process_similar_to,
                                     numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                     trim_after_tilde,
                                     verbose );
    }

    /**
     * 
     * 
     * 
     * @param phylogeny
     * @param map
     *            maps names (in phylogeny) to new values if intermediate_map is
     *            null otherwise maps intermediate value to new value
     * @param field
     * @param picky
     * @param intermediate_map
     *            maps name (in phylogeny) to a intermediate value
     * @throws IllegalArgumentException
     * @throws PhyloXmlDataFormatException 
     */
    public static void decorate( final Phylogeny phylogeny,
                                 final Map<String, String> map,
                                 final FIELD field,
                                 final boolean extract_bracketed_scientific_name,
                                 final boolean extract_bracketed_tax_code,
                                 final boolean picky,
                                 final Map<String, String> intermediate_map,
                                 final boolean cut_name_after_space,
                                 final boolean process_name_intelligently,
                                 final boolean process_similar_to,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                 final boolean trim_after_tilde,
                                 final boolean verbose ) throws IllegalArgumentException, PhyloXmlDataFormatException {
        if ( extract_bracketed_scientific_name && ( field == FIELD.TAXONOMY_SCIENTIFIC_NAME ) ) {
            throw new IllegalArgumentException( "attempt to extract bracketed scientific name together with data field pointing to scientific name" );
        }
        if ( map.isEmpty() ) {
            throw new IllegalArgumentException( "map is empty" );
        }
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            String name = node.getName();
            String tilde_annotation = null;
            if ( trim_after_tilde && ( name.indexOf( '~' ) > 0 ) ) {
                final int ti = name.indexOf( '~' );
                tilde_annotation = name.substring( ti );
                name = name.substring( 0, ti );
            }
            if ( !ForesterUtil.isEmpty( name ) ) {
                if ( intermediate_map != null ) {
                    name = PhylogenyDecorator.extractIntermediate( intermediate_map, name, verbose );
                }
                if ( map.containsKey( name ) || ( numbers_of_chars_allowed_to_remove_if_not_found_in_map > 0 ) ) {
                    String new_value = map.get( name );
                    int x = 0;
                    while ( ( new_value == null ) && ( numbers_of_chars_allowed_to_remove_if_not_found_in_map > 0 )
                            && ( x <= numbers_of_chars_allowed_to_remove_if_not_found_in_map ) ) {
                        new_value = map.get( name.substring( 0, name.length() - x ) );
                        ++x;
                    }
                    if ( new_value != null ) {
                        new_value = new_value.trim();
                        new_value.replaceAll( "/\\s+/", " " );
                        if ( extract_bracketed_scientific_name && new_value.endsWith( "]" ) ) {
                            new_value = extractBracketedScientificNames( node, new_value );
                        }
                        else if ( extract_bracketed_tax_code ) {
                            if ( ParserUtils.TAXOMONY_CODE_PATTERN_BRACKETED.matcher( new_value ).find() ) {
                                new_value = extractBracketedTaxCodes( node, new_value );
                            }
                            else if ( picky ) {
                                throw new IllegalArgumentException( " could not get taxonomy from \"" + new_value
                                        + "\"" );
                            }
                        }
                        switch ( field ) {
                            case MOL_SEQ:
                                if ( verbose ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                if ( !node.getNodeData().isHasSequence() ) {
                                    node.getNodeData().setSequence( new Sequence() );
                                }
                                node.getNodeData().getSequence().setMolecularSequence( new_value );
                                break;
                            case SEQUENCE_ANNOTATION_DESC:
                                if ( verbose ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                if ( !node.getNodeData().isHasSequence() ) {
                                    node.getNodeData().setSequence( new Sequence() );
                                }
                                final Annotation annotation = new Annotation();
                                annotation.setDesc( new_value );
                                node.getNodeData().getSequence().addAnnotation( annotation );
                                break;
                            case DOMAIN_STRUCTURE:
                                if ( verbose ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                if ( !node.getNodeData().isHasSequence() ) {
                                    node.getNodeData().setSequence( new Sequence() );
                                }
                                node.getNodeData().getSequence()
                                        .setDomainArchitecture( new DomainArchitecture( new_value ) );
                                break;
                            case TAXONOMY_CODE:
                                if ( verbose ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                ForesterUtil.ensurePresenceOfTaxonomy( node );
                                node.getNodeData().getTaxonomy().setTaxonomyCode( new_value );
                                break;
                            case TAXONOMY_SCIENTIFIC_NAME:
                                if ( verbose ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                ForesterUtil.ensurePresenceOfTaxonomy( node );
                                node.getNodeData().getTaxonomy().setScientificName( new_value );
                                break;
                            case SEQUENCE_NAME:
                                if ( trim_after_tilde ) {
                                    new_value = addTildeAnnotation( tilde_annotation, new_value );
                                }
                                if ( verbose ) {
                                    System.out.println( name + ": " + new_value );
                                }
                                if ( !node.getNodeData().isHasSequence() ) {
                                    node.getNodeData().setSequence( new Sequence() );
                                }
                                node.getNodeData().getSequence().setName( new_value );
                                break;
                            case NODE_NAME:
                                if ( verbose ) {
                                    System.out.print( name + " -> " );
                                }
                                if ( cut_name_after_space ) {
                                    if ( verbose ) {
                                        System.out.print( new_value + " -> " );
                                    }
                                    new_value = PhylogenyDecorator.deleteAtFirstSpace( new_value );
                                }
                                else if ( process_name_intelligently ) {
                                    if ( verbose ) {
                                        System.out.print( new_value + " -> " );
                                    }
                                    new_value = PhylogenyDecorator.processNameIntelligently( new_value );
                                }
                                else if ( process_similar_to ) {
                                    if ( verbose ) {
                                        System.out.print( new_value + " -> " );
                                    }
                                    new_value = PhylogenyDecorator.processSimilarTo( new_value );
                                }
                                if ( PhylogenyDecorator.SANITIZE ) {
                                    new_value = PhylogenyDecorator.sanitize( new_value );
                                }
                                if ( trim_after_tilde ) {
                                    new_value = addTildeAnnotation( tilde_annotation, new_value );
                                }
                                if ( verbose ) {
                                    System.out.println( new_value );
                                }
                                node.setName( new_value );
                                break;
                            default:
                                throw new RuntimeException( "unknown field \"" + field + "\"" );
                        }
                    }
                }
                else if ( picky ) {
                    throw new IllegalArgumentException( "\"" + name + "\" not found in name map" );
                }
            }
        }
    }

    public static void decorate( final Phylogeny[] phylogenies,
                                 final Map<String, Map<String, String>> map,
                                 final boolean picky,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map )
            throws IllegalArgumentException, NHXFormatException, PhyloXmlDataFormatException {
        for( final Phylogeny phylogenie : phylogenies ) {
            PhylogenyDecorator
                    .decorate( phylogenie, map, picky, numbers_of_chars_allowed_to_remove_if_not_found_in_map );
        }
    }

    public static void decorate( final Phylogeny[] phylogenies,
                                 final Map<String, String> map,
                                 final FIELD field,
                                 final boolean extract_bracketed_scientific_name,
                                 final boolean extract_bracketed_tax_code,
                                 final boolean picky,
                                 final boolean cut_name_after_space,
                                 final boolean process_name_intelligently,
                                 final boolean process_similar_to,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                 final boolean trim_after_tilde,
                                 final boolean verbose ) throws IllegalArgumentException, NHXFormatException,
            PhyloXmlDataFormatException {
        for( final Phylogeny phylogenie : phylogenies ) {
            PhylogenyDecorator.decorate( phylogenie,
                                         map,
                                         field,
                                         extract_bracketed_scientific_name,
                                         extract_bracketed_tax_code,
                                         picky,
                                         cut_name_after_space,
                                         process_name_intelligently,
                                         process_similar_to,
                                         numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                         trim_after_tilde,
                                         verbose );
        }
    }

    public static void decorate( final Phylogeny[] phylogenies,
                                 final Map<String, String> map,
                                 final FIELD field,
                                 final boolean extract_bracketed_scientific_name,
                                 final boolean extract_bracketed_tax_code,
                                 final boolean picky,
                                 final Map<String, String> intermediate_map,
                                 final boolean cut_name_after_space,
                                 final boolean process_name_intelligently,
                                 final boolean process_similar_to,
                                 final int numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                 final boolean trim_after_tilde,
                                 final boolean verbose ) throws IllegalArgumentException, NHXFormatException,
            PhyloXmlDataFormatException {
        for( final Phylogeny phylogenie : phylogenies ) {
            PhylogenyDecorator.decorate( phylogenie,
                                         map,
                                         field,
                                         extract_bracketed_scientific_name,
                                         extract_bracketed_tax_code,
                                         picky,
                                         intermediate_map,
                                         cut_name_after_space,
                                         process_name_intelligently,
                                         process_similar_to,
                                         numbers_of_chars_allowed_to_remove_if_not_found_in_map,
                                         trim_after_tilde,
                                         verbose );
        }
    }

    public static Map<String, Map<String, String>> parseMappingTable( final File mapping_table_file )
            throws IOException {
        final Map<String, Map<String, String>> map = new HashMap<String, Map<String, String>>();
        BasicTable<String> mapping_table = null;
        mapping_table = BasicTableParser.parse( mapping_table_file, '\t', false, false );
        for( int row = 0; row < mapping_table.getNumberOfRows(); ++row ) {
            final Map<String, String> row_map = new HashMap<String, String>();
            String name = null;
            for( int col = 0; col < mapping_table.getNumberOfColumns(); ++col ) {
                final String table_cell = mapping_table.getValue( col, row );
                if ( col == 0 ) {
                    name = table_cell;
                }
                else if ( table_cell != null ) {
                    final String key = table_cell.substring( 0, table_cell.indexOf( ':' ) );
                    final String val = table_cell.substring( table_cell.indexOf( ':' ) + 1, table_cell.length() );
                    row_map.put( key, val );
                }
            }
            map.put( name, row_map );
        }
        return map;
    }

    private final static String addTildeAnnotation( final String tilde_annotation, final String new_value ) {
        if ( ForesterUtil.isEmpty( tilde_annotation ) ) {
            return new_value;
        }
        return new_value + tilde_annotation;
    }

    private static String deleteAtFirstSpace( final String name ) {
        final int first_space = name.indexOf( " " );
        if ( first_space > 1 ) {
            return name.substring( 0, first_space ).trim();
        }
        return name;
    }

    private static String extractBracketedScientificNames( final PhylogenyNode node, final String new_value ) {
        final int i = new_value.lastIndexOf( "[" );
        final String scientific_name = new_value.substring( i + 1, new_value.length() - 1 );
        ForesterUtil.ensurePresenceOfTaxonomy( node );
        node.getNodeData().getTaxonomy().setScientificName( scientific_name );
        return new_value.substring( 0, i - 1 ).trim();
    }

    private static String extractBracketedTaxCodes( final PhylogenyNode node, final String new_value ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( new_value );
        final String tc = extractBracketedTaxCodes( sb );
        if ( !ForesterUtil.isEmpty( tc ) ) {
            ForesterUtil.ensurePresenceOfTaxonomy( node );
            try {
                node.getNodeData().getTaxonomy().setTaxonomyCode( tc );
            }
            catch ( final PhyloXmlDataFormatException e ) {
                throw new IllegalArgumentException( "illegal format for taxonomy code: " + tc );
            }
            return sb.toString().trim();
        }
        return new_value;
    }

    private static String extractBracketedTaxCodes( final StringBuilder sb ) {
        final Matcher m = ParserUtils.TAXOMONY_CODE_PATTERN_BRACKETED.matcher( sb );
        if ( m.find() ) {
            final String tc = m.group( 1 );
            sb.delete( m.start( 1 ) - 1, m.end( 1 ) + 1 );
            return tc;
        }
        return null;
    }

    private static String extractIntermediate( final Map<String, String> intermediate_map,
                                               final String name,
                                               final boolean verbose ) {
        String new_name = null;
        if ( verbose ) {
            System.out.print( name + " => " );
        }
        if ( intermediate_map.containsKey( name ) ) {
            new_name = intermediate_map.get( name );
            if ( ForesterUtil.isEmpty( new_name ) ) {
                throw new IllegalArgumentException( "\"" + name + "\" maps to null or empty string in secondary map" );
            }
        }
        else {
            throw new IllegalArgumentException( "\"" + name + "\" not found in name secondary map" );
        }
        if ( verbose ) {
            System.out.println( new_name + "  " );
        }
        return new_name;
    }

    private static String processNameIntelligently( final String name ) {
        final String[] s = name.split( " " );
        if ( s.length < 2 ) {
            return name;
        }
        else if ( ( s[ 0 ].indexOf( "_" ) > 0 ) && ( s[ 0 ].indexOf( "|" ) > 0 ) ) {
            return s[ 0 ];
        }
        else if ( ( s[ 1 ].indexOf( "_" ) > 0 ) && ( s[ 1 ].indexOf( "|" ) > 0 ) ) {
            return s[ 1 ];
        }
        else if ( ( s[ 0 ].indexOf( "_" ) > 0 ) && ( s[ 0 ].indexOf( "." ) > 0 ) ) {
            return s[ 0 ];
        }
        else if ( ( s[ 1 ].indexOf( "_" ) > 0 ) && ( s[ 1 ].indexOf( "." ) > 0 ) ) {
            return s[ 1 ];
        }
        else if ( s[ 0 ].indexOf( "_" ) > 0 ) {
            return s[ 0 ];
        }
        else if ( s[ 1 ].indexOf( "_" ) > 0 ) {
            return s[ 1 ];
        }
        else {
            return s[ 0 ];
        }
    }

    private static String processSimilarTo( final String name ) {
        final int i = name.toLowerCase().indexOf( "similar to" );
        String similar_to = "";
        if ( i >= 0 ) {
            similar_to = " similarity=" + name.substring( i + 10 ).trim();
        }
        final String pi = processNameIntelligently( name );
        return pi + similar_to;
    }

    private static String sanitize( String s ) {
        s = s.replace( ' ', '_' );
        s = s.replace( '(', '{' );
        s = s.replace( ')', '}' );
        s = s.replace( '[', '{' );
        s = s.replace( ']', '}' );
        s = s.replace( ',', '_' );
        return s;
    }

    public static enum FIELD {
        DOMAIN_STRUCTURE,
        MOL_SEQ,
        NODE_NAME,
        SEQUENCE_ANNOTATION_DESC,
        SEQUENCE_NAME,
        TAXONOMY_CODE,
        TAXONOMY_SCIENTIFIC_NAME;
    }
}
