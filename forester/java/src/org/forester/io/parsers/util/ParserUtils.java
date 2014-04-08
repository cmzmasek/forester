// $Id:
//
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
// WWW: www.phylosoft.org/

package org.forester.io.parsers.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.net.URL;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.tol.TolParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public final class ParserUtils {

    final public static String   TAX_CODE                        = "(?:[A-Z9][A-Z]{2}[A-Z0-9]{2})|RAT|PIG|PEA";
    final public static Pattern  TAXOMONY_CODE_PATTERN_A         = Pattern.compile( "(?:\\b|_)(" + TAX_CODE + ")\\b" );
    final public static Pattern  TAXOMONY_CODE_PATTERN_BRACKETED = Pattern.compile( "\\[(" + TAX_CODE + ")\\]" );
    final public static Pattern  TAXOMONY_CODE_PATTERN_PFR       = Pattern.compile( "(?:\\b|_)[a-zA-Z0-9]{3,}_("
                                                                         + TAX_CODE + ")\\b" );
    final public static Pattern  TAXOMONY_SN_PATTERN             = Pattern
                                                                         .compile( "(?:\\b|_)[a-zA-Z0-9]{3,}_([A-Z][a-z]+_[a-z]{2,}(?:_[a-z][a-z0-9_]+)?)\\b" );
    final public static Pattern  TAXOMONY_SN_PATTERN_SN          = Pattern
                                                                         .compile( "\\b([A-Z][a-z]+[_ ][a-z]{2,}(?:[_ ][a-z]+)?)(?:\\b|_)" );
    final private static Pattern TAXOMONY_CODE_PATTERN_PFS       = Pattern.compile( "(?:\\b|_)[A-Z0-9]{4,}_("
                                                                         + TAX_CODE + ")/\\d+-\\d+\\b" );
    final private static Pattern TAXOMONY_UNIPROT_ID_PATTERN_PFR = Pattern
                                                                         .compile( "(?:\\b|_)[A-Z0-9]{1,}_(\\d{1,7})\\b" );
    final private static Pattern TAXOMONY_UNIPROT_ID_PATTERN_PFS = Pattern
                                                                         .compile( "(?:\\b|_)[A-Z0-9]{4,}_(\\d{1,7})/\\d+-\\d+\\b" );

    final public static PhylogenyParser createParserDependingFileContents( final File file,
                                                                           final boolean phyloxml_validate_against_xsd )
            throws FileNotFoundException, IOException {
        PhylogenyParser parser = null;
        final String first_line = ForesterUtil.getFirstLine( file ).trim().toLowerCase();
        if ( first_line.startsWith( "<" ) ) {
            parser = PhyloXmlParser.createPhyloXmlParser();
            if ( phyloxml_validate_against_xsd ) {
                final ClassLoader cl = PhyloXmlParser.class.getClassLoader();
                final URL xsd_url = cl.getResource( ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE );
                if ( xsd_url != null ) {
                    ( ( PhyloXmlParser ) parser ).setValidateAgainstSchema( xsd_url.toString() );
                }
                else {
                    if ( ForesterConstants.RELEASE ) {
                        throw new RuntimeException( "failed to get URL for phyloXML XSD from jar file from ["
                                + ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE + "]" );
                    }
                }
            }
        }
        else if ( ( first_line.startsWith( "nexus" ) ) || ( first_line.startsWith( "#nexus" ) )
                || ( first_line.startsWith( "# nexus" ) ) || ( first_line.startsWith( "begin" ) ) ) {
            parser = new NexusPhylogeniesParser();
        }
        else {
            parser = new NHXParser();
        }
        return parser;
    }

    final public static PhylogenyParser createParserDependingOnFileType( final File file,
                                                                         final boolean phyloxml_validate_against_xsd )
            throws FileNotFoundException, IOException {
        PhylogenyParser parser = null;
        parser = ParserUtils.createParserDependingOnSuffix( file.getName(), phyloxml_validate_against_xsd );
        if ( parser == null ) {
            parser = createParserDependingFileContents( file, phyloxml_validate_against_xsd );
        }
        if ( ( parser != null ) && file.toString().toLowerCase().endsWith( ".zip" ) ) {
            if ( parser instanceof PhyloXmlParser ) {
                ( ( PhyloXmlParser ) parser ).setZippedInputstream( true );
            }
            else if ( parser instanceof TolParser ) {
                ( ( TolParser ) parser ).setZippedInputstream( true );
            }
        }
        return parser;
    }

    final public static PhylogenyParser createParserDependingOnUrlContents( final URL url,
                                                                            final boolean phyloxml_validate_against_xsd )
            throws FileNotFoundException, IOException {
        final String lc_filename = url.getFile().toString().toLowerCase();
        PhylogenyParser parser = createParserDependingOnSuffix( lc_filename, phyloxml_validate_against_xsd );
        if ( parser == null ) {
            final String first_line = ForesterUtil.getFirstLine( url ).trim().toLowerCase();
            if ( first_line.startsWith( "<" ) ) {
                parser = PhyloXmlParser.createPhyloXmlParser();
                if ( phyloxml_validate_against_xsd ) {
                    final ClassLoader cl = PhyloXmlParser.class.getClassLoader();
                    final URL xsd_url = cl.getResource( ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE );
                    if ( xsd_url != null ) {
                        ( ( PhyloXmlParser ) parser ).setValidateAgainstSchema( xsd_url.toString() );
                    }
                    else {
                        throw new RuntimeException( "failed to get URL for phyloXML XSD from jar file from ["
                                + ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE + "]" );
                    }
                }
            }
            else if ( ( first_line.startsWith( "nexus" ) ) || ( first_line.startsWith( "#nexus" ) )
                    || ( first_line.startsWith( "# nexus" ) ) || ( first_line.startsWith( "begin" ) ) ) {
                parser = new NexusPhylogeniesParser();
            }
            else {
                parser = new NHXParser();
            }
        }
        if ( ( parser != null ) && lc_filename.endsWith( ".zip" ) ) {
            if ( parser instanceof PhyloXmlParser ) {
                ( ( PhyloXmlParser ) parser ).setZippedInputstream( true );
            }
            else if ( parser instanceof TolParser ) {
                ( ( TolParser ) parser ).setZippedInputstream( true );
            }
        }
        return parser;
    }

    public static BufferedReader createReader( final Object source ) throws IOException, FileNotFoundException {
        BufferedReader reader = null;
        if ( ( source instanceof File ) || ( source instanceof String ) ) {
            File f = null;
            if ( source instanceof File ) {
                f = ( File ) source;
            }
            else {
                f = new File( ( String ) source );
            }
            if ( !f.exists() ) {
                throw new IOException( "[" + f.getAbsolutePath() + "] does not exist" );
            }
            else if ( !f.isFile() ) {
                throw new IOException( "[" + f.getAbsolutePath() + "] is not a file" );
            }
            else if ( !f.canRead() ) {
                throw new IOException( "[" + f.getAbsolutePath() + "] is not a readable" );
            }
            reader = new BufferedReader( new FileReader( f ) );
        }
        else if ( source instanceof InputStream ) {
            reader = new BufferedReader( new InputStreamReader( ( InputStream ) source ) );
        }
        else if ( ( source instanceof StringBuffer ) || ( source instanceof StringBuilder ) ) {
            reader = new BufferedReader( new StringReader( source.toString() ) );
        }
        else {
            throw new IllegalArgumentException( "attempt to parse object of type [" + source.getClass()
                    + "] (can only parse objects of type File/String, InputStream, StringBuffer, or StringBuilder)" );
        }
        return reader;
    }

    public final static String extractScientificNameFromNodeName( final String name ) {
        final Matcher m = TAXOMONY_SN_PATTERN.matcher( name );
        if ( m.find() ) {
            return m.group( 1 ).replace( '_', ' ' );
        }
        final Matcher m_sn = TAXOMONY_SN_PATTERN_SN.matcher( name );
        if ( m_sn.find() ) {
            return m_sn.group( 1 ).replace( '_', ' ' );
        }
        return null;
    }

    public final static String extractTaxonomyCodeFromNodeName( final String name,
                                                                final TAXONOMY_EXTRACTION taxonomy_extraction ) {
        Matcher m = TAXOMONY_CODE_PATTERN_PFS.matcher( name );
        if ( m.find() ) {
            return m.group( 1 );
        }
        else if ( ( taxonomy_extraction == TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED )
                || ( taxonomy_extraction == TAXONOMY_EXTRACTION.AGGRESSIVE ) ) {
            m = TAXOMONY_CODE_PATTERN_PFR.matcher( name );
            if ( m.find() ) {
                return m.group( 1 );
            }
            else if ( taxonomy_extraction == TAXONOMY_EXTRACTION.AGGRESSIVE ) {
                m = TAXOMONY_CODE_PATTERN_A.matcher( name );
                if ( m.find() ) {
                    return m.group( 1 );
                }
            }
        }
        return null;
    }

    public final static String extractTaxonomyDataFromNodeName( final PhylogenyNode node,
                                                                final NHXParser.TAXONOMY_EXTRACTION taxonomy_extraction )
            throws PhyloXmlDataFormatException {
        if ( taxonomy_extraction == TAXONOMY_EXTRACTION.NO ) {
            throw new IllegalArgumentException();
        }
        final String id = extractUniprotTaxonomyIdFromNodeName( node.getName(), taxonomy_extraction );
        if ( !ForesterUtil.isEmpty( id ) ) {
            if ( !node.getNodeData().isHasTaxonomy() ) {
                node.getNodeData().setTaxonomy( new Taxonomy() );
            }
            node.getNodeData().getTaxonomy().setIdentifier( new Identifier( id, "uniprot" ) );
            return id;
        }
        else {
            final String code = extractTaxonomyCodeFromNodeName( node.getName(), taxonomy_extraction );
            if ( !ForesterUtil.isEmpty( code ) ) {
                if ( !node.getNodeData().isHasTaxonomy() ) {
                    node.getNodeData().setTaxonomy( new Taxonomy() );
                }
                node.getNodeData().getTaxonomy().setTaxonomyCode( code );
                return code;
            }
            else if ( taxonomy_extraction == TAXONOMY_EXTRACTION.AGGRESSIVE ) {
                final String sn = extractScientificNameFromNodeName( node.getName() );
                if ( !ForesterUtil.isEmpty( sn ) ) {
                    if ( !node.getNodeData().isHasTaxonomy() ) {
                        node.getNodeData().setTaxonomy( new Taxonomy() );
                    }
                    node.getNodeData().getTaxonomy().setScientificName( sn );
                    return sn;
                }
            }
        }
        return null;
    }

    public final static String extractUniprotTaxonomyIdFromNodeName( final String name,
                                                                     final TAXONOMY_EXTRACTION taxonomy_extraction ) {
        Matcher m = TAXOMONY_UNIPROT_ID_PATTERN_PFS.matcher( name );
        if ( m.find() ) {
            return m.group( 1 );
        }
        else if ( ( taxonomy_extraction == TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED )
                || ( taxonomy_extraction == TAXONOMY_EXTRACTION.AGGRESSIVE ) ) {
            m = TAXOMONY_UNIPROT_ID_PATTERN_PFR.matcher( name );
            if ( m.find() ) {
                return m.group( 1 );
            }
            //else if ( taxonomy_extraction == TAXONOMY_EXTRACTION.AGGRESSIVE ) {
            //    m = TAXOMONY_UNIPROT_ID_PATTERN_A.matcher( name );
            //    if ( m.find() ) {
            //        return m.group( 1 );
            //    }
            //}
        }
        return null;
    }

    public final static Phylogeny[] readPhylogenies( final File file ) throws FileNotFoundException, IOException {
        return PhylogenyMethods.readPhylogenies( ParserUtils.createParserDependingOnFileType( file, true ), file );
    }

    public final static Phylogeny[] readPhylogenies( final String file_name ) throws FileNotFoundException, IOException {
        return readPhylogenies( new File( file_name ) );
    }

    /**
     * Return null if it can not guess the parser to use based on name suffix.
     * 
     * @param filename
     * @return
     */
    final private static PhylogenyParser createParserDependingOnSuffix( final String filename,
                                                                        final boolean phyloxml_validate_against_xsd ) {
        PhylogenyParser parser = null;
        final String filename_lc = filename.toLowerCase();
        if ( filename_lc.endsWith( ".tol" ) || filename_lc.endsWith( ".tolxml" ) || filename_lc.endsWith( ".tol.zip" ) ) {
            parser = new TolParser();
        }
        else if ( filename_lc.endsWith( ".xml" ) || filename_lc.endsWith( "phyloxml" ) || filename_lc.endsWith( ".zip" ) ) {
            parser = PhyloXmlParser.createPhyloXmlParser();
            if ( phyloxml_validate_against_xsd ) {
                final ClassLoader cl = PhyloXmlParser.class.getClassLoader();
                final URL xsd_url = cl.getResource( ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE );
                if ( xsd_url != null ) {
                    ( ( PhyloXmlParser ) parser ).setValidateAgainstSchema( xsd_url.toString() );
                }
                else {
                    if ( ForesterConstants.RELEASE ) {
                        throw new RuntimeException( "failed to get URL for phyloXML XSD from jar file from ["
                                + ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE + "]" );
                    }
                }
            }
        }
        else if ( filename_lc.endsWith( ".nexus" ) || filename_lc.endsWith( ".nex" ) || filename_lc.endsWith( ".nx" ) ) {
            parser = new NexusPhylogeniesParser();
        }
        else if ( filename_lc.endsWith( ".nhx" ) || filename_lc.endsWith( ".nh" ) || filename_lc.endsWith( ".newick" )
                || filename_lc.endsWith( ".nwk" ) ) {
            parser = new NHXParser();
        }
        return parser;
    }
}
