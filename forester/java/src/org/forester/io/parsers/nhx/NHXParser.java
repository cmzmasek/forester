// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2013 Christian M. Zmasek
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

package org.forester.io.parsers.nhx;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.IteratingPhylogenyParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.parsers.util.PhylogenyParserException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public final class NHXParser implements PhylogenyParser, IteratingPhylogenyParser {

    public final static Pattern             MB_BL_PATTERN                              = Pattern
                                                                                               .compile( "length_median=([^,]+)" );
    public final static Pattern             MB_PROB_PATTERN                            = Pattern
                                                                                               .compile( "prob=([^,]+)" );
    public final static Pattern             MB_PROB_SD_PATTERN                         = Pattern
                                                                                               .compile( "prob_stddev=([^,]+)" );
    public final static Pattern             NUMBERS_ONLY_PATTERN                       = Pattern
                                                                                               .compile( "^[0-9\\.]+$" );
    final static public boolean             REPLACE_UNDERSCORES_DEFAULT                = false;
    public static final TAXONOMY_EXTRACTION TAXONOMY_EXTRACTION_DEFAULT                = TAXONOMY_EXTRACTION.NO;
    private static final boolean            ALLOW_ERRORS_IN_DISTANCE_TO_PARENT_DEFAULT = false;
    final static private byte               BUFFERED_READER                            = 3;
    final static private byte               CHAR_ARRAY                                 = 2;
    final static private boolean            GUESS_IF_SUPPORT_VALUES                    = true;
    final static private boolean            GUESS_ROOTEDNESS_DEFAULT                   = true;
    final static private boolean            IGNORE_QUOTES_DEFAULT                      = false;
    final static private byte               STRING                                     = 0;
    final static private byte               STRING_BUFFER                              = 1;
    final static private byte               STRING_BUILDER                             = 4;
    private boolean                         _allow_errors_in_distance_to_parent;
    private int                             _clade_level;
    private StringBuilder                   _current_anotation;
    private PhylogenyNode                   _current_node;
    private Phylogeny                       _current_phylogeny;
    private boolean                         _guess_rootedness;
    private int                             _i;
    private boolean                         _ignore_quotes;
    private boolean                         _in_comment                                = false;
    private boolean                         _in_double_quote                           = false;
    private boolean                         _in_open_bracket                           = false;
    private boolean                         _in_single_quote                           = false;
    private byte                            _input_type;
    private BufferedReader                  _my_source_br                              = null;
    private char[]                          _my_source_charary                         = null;
    private StringBuffer                    _my_source_sbuff                           = null;
    private StringBuilder                   _my_source_sbuil                           = null;
    private String                          _my_source_str                             = null;
    private Phylogeny                       _next;
    private Object                          _nhx_source;
    private boolean                         _replace_underscores;
    private boolean                         _saw_closing_paren;
    private boolean                         _saw_colon                                 = false;
    private boolean                         _saw_open_bracket                          = false;
    private Object                          _source;
    private int                             _source_length;
    private TAXONOMY_EXTRACTION             _taxonomy_extraction;

    public NHXParser() {
        init();
    }

    @Override
    public String getName() {
        return "NN/NHX Parser";
    }

    public final TAXONOMY_EXTRACTION getTaxonomyExtraction() {
        return _taxonomy_extraction;
    }

    @Override
    public final boolean hasNext() {
        return _next != null;
    }

    @Override
    public final Phylogeny next() throws NHXFormatException, IOException {
        final Phylogeny phy = _next;
        parseNext();
        return phy;
    }

    @Override
    public final Phylogeny[] parse() throws IOException {
        final List<Phylogeny> l = new ArrayList<Phylogeny>();
        int c = 0;
        while ( hasNext() ) {
            l.add( next() );
            c++;
        }
        final Phylogeny[] p = new Phylogeny[ l.size() ];
        for( int i = 0; i < l.size(); ++i ) {
            p[ i ] = l.get( i );
        }
        reset();
        return p;
    }

    @Override
    public final void reset() throws NHXFormatException, IOException {
        _i = 0;
        _next = null;
        _in_comment = false;
        _saw_colon = false;
        _saw_open_bracket = false;
        _in_open_bracket = false;
        _in_double_quote = false;
        _in_single_quote = false;
        _clade_level = 0;
        _current_anotation = new StringBuilder();
        _current_phylogeny = null;
        _current_node = null;
        _my_source_str = null;
        _my_source_sbuff = null;
        _my_source_sbuil = null;
        _my_source_charary = null;
        determineSourceType( _source );
        switch ( _input_type ) {
            case STRING:
                _my_source_br = null;
                _my_source_str = ( String ) _nhx_source;
                break;
            case STRING_BUFFER:
                _my_source_br = null;
                _my_source_sbuff = ( StringBuffer ) _nhx_source;
                break;
            case STRING_BUILDER:
                _my_source_br = null;
                _my_source_sbuil = ( StringBuilder ) _nhx_source;
                break;
            case CHAR_ARRAY:
                _my_source_br = null;
                _my_source_charary = ( char[] ) _nhx_source;
                break;
            case BUFFERED_READER:
                _my_source_br = ( BufferedReader ) _nhx_source;
                break;
            default:
                throw new RuntimeException( "unknown input type" );
        }
        parseNext();
    }

    public final void setGuessRootedness( final boolean guess_rootedness ) {
        _guess_rootedness = guess_rootedness;
    }

    public final void setIgnoreQuotes( final boolean ignore_quotes ) {
        _ignore_quotes = ignore_quotes;
    }

    public final void setReplaceUnderscores( final boolean replace_underscores ) {
        _replace_underscores = replace_underscores;
    }

    @Override
    public final void setSource( final Object nhx_source ) throws NHXFormatException, IOException {
        _source = nhx_source;
        reset();
    }

    public final void setTaxonomyExtraction( final TAXONOMY_EXTRACTION taxonomy_extraction ) {
        _taxonomy_extraction = taxonomy_extraction;
    }

    public final void setAllowErrorsInDistanceToParent( final boolean allow_errors_in_distance_to_parent ) {
        _allow_errors_in_distance_to_parent = allow_errors_in_distance_to_parent;
    }

    private final void determineSourceType( final Object nhx_source ) throws IOException {
        if ( nhx_source == null ) {
            throw new PhylogenyParserException( getClass() + ": attempt to parse null object." );
        }
        else if ( nhx_source instanceof String ) {
            _input_type = NHXParser.STRING;
            _source_length = ( ( String ) nhx_source ).length();
            _nhx_source = nhx_source;
        }
        else if ( nhx_source instanceof StringBuilder ) {
            _input_type = NHXParser.STRING_BUILDER;
            _source_length = ( ( StringBuilder ) nhx_source ).length();
            _nhx_source = nhx_source;
        }
        else if ( nhx_source instanceof StringBuffer ) {
            _input_type = NHXParser.STRING_BUFFER;
            _source_length = ( ( StringBuffer ) nhx_source ).length();
            _nhx_source = nhx_source;
        }
        else if ( nhx_source instanceof StringBuilder ) {
            _input_type = NHXParser.STRING_BUILDER;
            _source_length = ( ( StringBuilder ) nhx_source ).length();
            _nhx_source = nhx_source;
        }
        else if ( nhx_source instanceof char[] ) {
            _input_type = NHXParser.CHAR_ARRAY;
            _source_length = ( ( char[] ) nhx_source ).length;
            _nhx_source = nhx_source;
        }
        else if ( nhx_source instanceof File ) {
            _input_type = NHXParser.BUFFERED_READER;
            _source_length = 0;
            if ( _my_source_br != null ) {
                try {
                    _my_source_br.close();
                }
                catch ( final IOException e ) {
                }
            }
            final File f = ( File ) nhx_source;
            final String error = ForesterUtil.isReadableFile( f );
            if ( !ForesterUtil.isEmpty( error ) ) {
                throw new PhylogenyParserException( error );
            }
            _nhx_source = new BufferedReader( new FileReader( f ) );
        }
        else if ( nhx_source instanceof URL ) {
            _input_type = NHXParser.BUFFERED_READER;
            _source_length = 0;
            if ( _my_source_br != null ) {
                try {
                    _my_source_br.close();
                }
                catch ( final IOException e ) {
                }
            }
            final InputStreamReader isr = new InputStreamReader( ( ( URL ) nhx_source ).openStream() );
            _nhx_source = new BufferedReader( isr );
        }
        else if ( nhx_source instanceof InputStream ) {
            _input_type = NHXParser.BUFFERED_READER;
            _source_length = 0;
            if ( _my_source_br != null ) {
                try {
                    _my_source_br.close();
                }
                catch ( final IOException e ) {
                }
            }
            final InputStreamReader isr = new InputStreamReader( ( InputStream ) nhx_source );
            _nhx_source = new BufferedReader( isr );
        }
        else {
            throw new IllegalArgumentException( getClass() + " can only parse objects of type String,"
                    + " StringBuffer, StringBuilder, char[], File, InputStream, or URL "
                    + " [attempt to parse object of " + nhx_source.getClass() + "]." );
        }
    }

    private final Phylogeny finishPhylogeny() throws PhylogenyParserException, NHXFormatException,
            PhyloXmlDataFormatException {
        if ( _current_phylogeny != null ) {
            parseNHX( _current_anotation != null ? _current_anotation.toString() : "",
                      _current_phylogeny.getRoot(),
                      getTaxonomyExtraction(),
                      isReplaceUnderscores(),
                      isAllowErrorsInDistanceToParent() );
            if ( GUESS_IF_SUPPORT_VALUES ) {
                if ( isBranchLengthsLikeBootstrapValues( _current_phylogeny ) ) {
                    moveBranchLengthsToConfidenceValues( _current_phylogeny );
                }
            }
            if ( isGuessRootedness() ) {
                final PhylogenyNode root = _current_phylogeny.getRoot();
                if ( ( root.getDistanceToParent() >= 0.0 ) || !ForesterUtil.isEmpty( root.getName() )
                        || !ForesterUtil.isEmpty( PhylogenyMethods.getSpecies( root ) ) || root.isHasAssignedEvent() ) {
                    _current_phylogeny.setRooted( true );
                }
            }
            return _current_phylogeny;
        }
        return null;
    }

    private final Phylogeny finishSingleNodePhylogeny() throws PhylogenyParserException, NHXFormatException,
            PhyloXmlDataFormatException {
        final PhylogenyNode new_node = new PhylogenyNode();
        parseNHX( _current_anotation.toString(),
                  new_node,
                  getTaxonomyExtraction(),
                  isReplaceUnderscores(),
                  isAllowErrorsInDistanceToParent() );
        _current_phylogeny = new Phylogeny();
        _current_phylogeny.setRoot( new_node );
        return _current_phylogeny;
    }

    private final void init() {
        setTaxonomyExtraction( TAXONOMY_EXTRACTION_DEFAULT );
        setReplaceUnderscores( REPLACE_UNDERSCORES_DEFAULT );
        setGuessRootedness( GUESS_ROOTEDNESS_DEFAULT );
        setIgnoreQuotes( IGNORE_QUOTES_DEFAULT );
        setAllowErrorsInDistanceToParent( ALLOW_ERRORS_IN_DISTANCE_TO_PARENT_DEFAULT );
    }

    private final boolean isAllowErrorsInDistanceToParent() {
        return _allow_errors_in_distance_to_parent;
    }

    private final boolean isGuessRootedness() {
        return _guess_rootedness;
    }

    private final boolean isIgnoreQuotes() {
        return _ignore_quotes;
    }

    private final boolean isReplaceUnderscores() {
        return _replace_underscores;
    }

    private final void parseNext() throws IOException, NHXFormatException {
        if ( _source == null ) {
            throw new IOException( "source is not set" );
        }
        while ( true ) {
            char c = '\b';
            if ( _input_type == BUFFERED_READER ) {
                final int ci = _my_source_br.read();
                if ( ci >= 0 ) {
                    c = ( char ) ci;
                }
                else {
                    break;
                }
            }
            else {
                if ( _i >= _source_length ) {
                    break;
                }
                else {
                    switch ( _input_type ) {
                        case STRING:
                            c = _my_source_str.charAt( _i );
                            break;
                        case STRING_BUFFER:
                            c = _my_source_sbuff.charAt( _i );
                            break;
                        case STRING_BUILDER:
                            c = _my_source_sbuil.charAt( _i );
                            break;
                        case CHAR_ARRAY:
                            c = _my_source_charary[ _i ];
                            break;
                    }
                }
            }
            if ( !_in_single_quote && !_in_double_quote ) {
                if ( c == ':' ) {
                    _saw_colon = true;
                }
                else if ( !( ( c < 33 ) || ( c > 126 ) ) && _saw_colon
                        && ( ( c != '[' ) && ( c != '.' ) && ( ( c < 48 ) || ( c > 57 ) ) ) ) {
                    _saw_colon = false;
                }
                if ( _in_open_bracket && ( c == ']' ) ) {
                    _in_open_bracket = false;
                }
            }
            // \n\t is always ignored,
            // as is " (34) and ' (39) (space is 32):
            if ( ( isIgnoreQuotes() && ( ( c < 33 ) || ( c > 126 ) || ( c == 34 ) || ( c == 39 ) || ( ( _clade_level == 0 ) && ( c == ';' ) ) ) )
                    || ( !isIgnoreQuotes() && ( ( c < 32 ) || ( c > 126 ) || ( ( _clade_level == 0 ) && ( c == ';' ) ) ) ) ) {
                //do nothing
            }
            else if ( ( c == 32 ) && ( !_in_single_quote && !_in_double_quote ) ) {
                //do nothing
            }
            else if ( _in_comment ) {
                if ( c == ']' ) {
                    _in_comment = false;
                }
            }
            else if ( _in_double_quote ) {
                if ( c == '"' ) {
                    _in_double_quote = false;
                }
                else {
                    _current_anotation.append( c );
                }
            }
            else if ( c == '"' ) {
                _in_double_quote = true;
            }
            else if ( _in_single_quote ) {
                if ( c == 39 ) {
                    _in_single_quote = false;
                }
                else {
                    _current_anotation.append( c );
                }
            }
            else if ( c == 39 ) {
                _in_single_quote = true;
            }
            else if ( c == '[' ) {
                _saw_open_bracket = true;
                _in_open_bracket = true;
            }
            else if ( _saw_open_bracket ) {
                if ( c != ']' ) {
                    // everything not starting with "[&" is considered a comment
                    // unless ":digits and/or . [bootstrap]":
                    if ( c == '&' ) {
                        _current_anotation.append( "[&" );
                    }
                    else if ( _saw_colon ) {
                        _current_anotation.append( "[" + c );
                    }
                    else {
                        _in_comment = true;
                    }
                }
                // comment consisting just of "[]":
                _saw_open_bracket = false;
            }
            else if ( ( c == '(' ) && !_in_open_bracket ) {
                final Phylogeny phy = processOpenParen();
                if ( phy != null ) {
                    ++_i;
                    //  return phy;
                    _next = phy;
                    return;
                }
            }
            else if ( ( c == ')' ) && !_in_open_bracket ) {
                processCloseParen();
            }
            else if ( ( c == ',' ) && !_in_open_bracket ) {
                processComma();
            }
            else {
                _current_anotation.append( c );
            }
            ++_i;
        } //  while ( true ) 
        if ( _clade_level != 0 ) {
            throw new PhylogenyParserException( "error in NH (Newick) formatted data: most likely cause: number of open parens does not equal number of close parens" );
        }
        if ( _current_phylogeny != null ) {
            _next = finishPhylogeny();
            _current_phylogeny = null;
            _current_anotation = null;
        }
        else if ( ( _current_anotation != null ) && ( _current_anotation.length() > 0 ) ) {
            _next = finishSingleNodePhylogeny();
            _current_anotation = null;
        }
        else {
            _next = null;
        }
    }

    private final void processCloseParen() throws PhylogenyParserException, NHXFormatException,
            PhyloXmlDataFormatException {
        if ( _clade_level < 0 ) {
            throw new PhylogenyParserException( "error in NH (Newick)/NHX formatted data: most likely cause: number of close parens is larger than number of open parens" );
        }
        --_clade_level;
        if ( !_saw_closing_paren ) {
            final PhylogenyNode new_node = new PhylogenyNode();
            parseNHX( _current_anotation.toString(),
                      new_node,
                      getTaxonomyExtraction(),
                      isReplaceUnderscores(),
                      isAllowErrorsInDistanceToParent() );
            _current_anotation = new StringBuilder();
            _current_node.addAsChild( new_node );
        }
        else {
            parseNHX( _current_anotation.toString(),
                      _current_node.getLastChildNode(),
                      getTaxonomyExtraction(),
                      isReplaceUnderscores(),
                      isAllowErrorsInDistanceToParent() );
            _current_anotation = new StringBuilder();
        }
        if ( !_current_node.isRoot() ) {
            _current_node = _current_node.getParent();
        }
        _saw_closing_paren = true;
    }

    private final void processComma() throws PhylogenyParserException, NHXFormatException, PhyloXmlDataFormatException {
        if ( !_saw_closing_paren ) {
            final PhylogenyNode new_node = new PhylogenyNode();
            parseNHX( _current_anotation.toString(),
                      new_node,
                      getTaxonomyExtraction(),
                      isReplaceUnderscores(),
                      isAllowErrorsInDistanceToParent() );
            if ( _current_node == null ) {
                throw new NHXFormatException( "format might not be NH or NHX" );
            }
            _current_node.addAsChild( new_node );
        }
        else {
            parseNHX( _current_anotation.toString(),
                      _current_node.getLastChildNode(),
                      getTaxonomyExtraction(),
                      isReplaceUnderscores(),
                      isAllowErrorsInDistanceToParent() );
        }
        _current_anotation = new StringBuilder();
        _saw_closing_paren = false;
    }

    private final Phylogeny processOpenParen() throws PhylogenyParserException, NHXFormatException,
            PhyloXmlDataFormatException {
        Phylogeny phy = null;
        final PhylogenyNode new_node = new PhylogenyNode();
        if ( _clade_level == 0 ) {
            if ( _current_phylogeny != null ) {
                phy = finishPhylogeny();
            }
            _clade_level = 1;
            _current_anotation = new StringBuilder();
            _current_phylogeny = new Phylogeny();
            _current_phylogeny.setRoot( new_node );
        }
        else {
            ++_clade_level;
            _current_node.addAsChild( new_node );
        }
        _current_node = new_node;
        _saw_closing_paren = false;
        return phy;
    }

    public final static NHXParser createInstance( final Object nhx_source ) throws NHXFormatException, IOException {
        final NHXParser parser = new NHXParser();
        parser.setSource( nhx_source );
        return parser;
    }

    public final static Phylogeny[] parse( final Object nhx_source ) throws NHXFormatException, IOException {
        return NHXParser.createInstance( nhx_source ).parse();
    }

    public final static void parseNHX( String s,
                                       final PhylogenyNode node_to_annotate,
                                       final TAXONOMY_EXTRACTION taxonomy_extraction,
                                       final boolean replace_underscores,
                                       final boolean allow_errors_in_distance_to_parent ) throws NHXFormatException,
            PhyloXmlDataFormatException {
        if ( ( taxonomy_extraction != TAXONOMY_EXTRACTION.NO ) && replace_underscores ) {
            throw new IllegalArgumentException( "cannot extract taxonomies and replace under scores at the same time" );
        }
        if ( ( s != null ) && ( s.length() > 0 ) ) {
            if ( replace_underscores ) {
                s = s.replaceAll( "_+", " " );
            }
            boolean is_nhx = false;
            final int ob = s.indexOf( "[" );
            if ( ob > -1 ) {
                String b = "";
                is_nhx = true;
                final int cb = s.indexOf( "]" );
                if ( cb < 0 ) {
                    throw new NHXFormatException( "error in NHX formatted data: no closing \"]\" in \"" + s + "\"" );
                }
                if ( s.indexOf( "&&NHX" ) == ( ob + 1 ) ) {
                    b = s.substring( ob + 6, cb );
                }
                else {
                    // No &&NHX and digits only: is likely to be a support value.
                    final String bracketed = s.substring( ob + 1, cb );
                    final Matcher numbers_only = NUMBERS_ONLY_PATTERN.matcher( bracketed );
                    if ( numbers_only.matches() ) {
                        b = ":" + NHXtags.SUPPORT + bracketed;
                    }
                    else if ( s.indexOf( "prob=" ) > -1 ) {
                        processMrBayes3Data( s, node_to_annotate );
                    }
                }
                s = s.substring( 0, ob ) + b;
                if ( ( s.indexOf( "[" ) > -1 ) || ( s.indexOf( "]" ) > -1 ) ) {
                    throw new NHXFormatException( "error in NHX formatted data: more than one \"]\" or \"[\"" );
                }
            }
            final StringTokenizer t = new StringTokenizer( s, ":" );
            if ( t.countTokens() > 0 ) {
                if ( !s.startsWith( ":" ) ) {
                    node_to_annotate.setName( t.nextToken() );
                    if ( !replace_underscores && ( !is_nhx && ( taxonomy_extraction != TAXONOMY_EXTRACTION.NO ) ) ) {
                        ParserUtils.extractTaxonomyDataFromNodeName( node_to_annotate, taxonomy_extraction );
                    }
                }
                while ( t.hasMoreTokens() ) {
                    s = t.nextToken();
                    if ( s.startsWith( NHXtags.SPECIES_NAME ) ) {
                        if ( !node_to_annotate.getNodeData().isHasTaxonomy() ) {
                            node_to_annotate.getNodeData().setTaxonomy( new Taxonomy() );
                        }
                        node_to_annotate.getNodeData().getTaxonomy().setScientificName( s.substring( 2 ) );
                    }
                    else if ( s.startsWith( NHXtags.IS_DUPLICATION ) ) {
                        if ( ( s.charAt( 2 ) == 'Y' ) || ( s.charAt( 2 ) == 'T' ) ) {
                            node_to_annotate.getNodeData().setEvent( Event.createSingleDuplicationEvent() );
                        }
                        else if ( ( s.charAt( 2 ) == 'N' ) || ( s.charAt( 2 ) == 'F' ) ) {
                            node_to_annotate.getNodeData().setEvent( Event.createSingleSpeciationEvent() );
                        }
                        else if ( s.charAt( 2 ) == '?' ) {
                            node_to_annotate.getNodeData().setEvent( Event.createSingleSpeciationOrDuplicationEvent() );
                        }
                        else {
                            throw new NHXFormatException( "error in NHX formatted data: :D=Y or :D=N or :D=?" );
                        }
                    }
                    else if ( s.startsWith( NHXtags.SUPPORT ) ) {
                        PhylogenyMethods.setConfidence( node_to_annotate, doubleValue( s.substring( 2 ), false ) );
                    }
                    else if ( s.startsWith( NHXtags.TAXONOMY_ID ) ) {
                        if ( !node_to_annotate.getNodeData().isHasTaxonomy() ) {
                            node_to_annotate.getNodeData().setTaxonomy( new Taxonomy() );
                        }
                        node_to_annotate.getNodeData().getTaxonomy().setIdentifier( new Identifier( s.substring( 2 ) ) );
                    }
                    else if ( s.startsWith( NHXtags.SEQUENCE_ACCESSION ) ) {
                        if ( !node_to_annotate.getNodeData().isHasSequence() ) {
                            node_to_annotate.getNodeData().setSequence( new Sequence() );
                        }
                        node_to_annotate.getNodeData().getSequence()
                                .setAccession( new Accession( s.substring( 3 ), "?" ) );
                    }
                    else if ( s.startsWith( NHXtags.GENE_NAME ) ) {
                        if ( !node_to_annotate.getNodeData().isHasSequence() ) {
                            node_to_annotate.getNodeData().setSequence( new Sequence() );
                        }
                        node_to_annotate.getNodeData().getSequence().setName( s.substring( 3 ) );
                    }
                    else if ( s.indexOf( '=' ) < 0 ) {
                        if ( node_to_annotate.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT
                                && !allow_errors_in_distance_to_parent ) {
                            throw new NHXFormatException( "error in NHX formatted data: more than one distance to parent:"
                                    + "\"" + s + "\"" );
                        }
                        node_to_annotate.setDistanceToParent( doubleValue( s, allow_errors_in_distance_to_parent ) );
                    }
                } // while ( t.hasMoreTokens() ) 
            }
        }
    }

    private final static double doubleValue( final String str, final boolean allow_errors ) throws NHXFormatException {
        try {
            return Double.valueOf( str ).doubleValue();
        }
        catch ( final NumberFormatException ex ) {
            if ( !allow_errors ) {
                throw new NHXFormatException( "error in NH/NHX formatted data: failed to parse number from " + "\""
                        + str + "\"" );
            }
        }
        return 0.0;
    }

    private final static boolean isBranchLengthsLikeBootstrapValues( final Phylogeny p ) {
        final PhylogenyNodeIterator it = p.iteratorExternalForward();
        final double d0 = it.next().getDistanceToParent();
        if ( ( d0 < 10 ) || !it.hasNext() ) {
            return false;
        }
        while ( it.hasNext() ) {
            final double d = it.next().getDistanceToParent();
            if ( ( d != d0 ) || ( d < 10 ) ) {
                return false;
            }
        }
        return true;
    }

    private final static void moveBranchLengthsToConfidenceValues( final Phylogeny p ) {
        final PhylogenyNodeIterator it = p.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            PhylogenyMethods.setBootstrapConfidence( n, n.getDistanceToParent() );
            n.setDistanceToParent( PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
        }
    }

    private final static void processMrBayes3Data( final String s, final PhylogenyNode node_to_annotate )
            throws NHXFormatException {
        double sd = -1;
        final Matcher mb_prob_sd_matcher = MB_PROB_SD_PATTERN.matcher( s );
        if ( mb_prob_sd_matcher.find() ) {
            try {
                sd = Double.parseDouble( mb_prob_sd_matcher.group( 1 ) );
            }
            catch ( final NumberFormatException e ) {
                throw new NHXFormatException( "failed to parse probability standard deviation (Mr Bayes output) from \""
                        + s + "\"" );
            }
        }
        final Matcher mb_prob_matcher = MB_PROB_PATTERN.matcher( s );
        if ( mb_prob_matcher.find() ) {
            double prob = -1;
            try {
                prob = Double.parseDouble( mb_prob_matcher.group( 1 ) );
            }
            catch ( final NumberFormatException e ) {
                throw new NHXFormatException( "failed to parse probability (Mr Bayes output) from \"" + s + "\"" );
            }
            if ( prob >= 0.0 ) {
                if ( sd >= 0.0 ) {
                    node_to_annotate.getBranchData()
                            .addConfidence( new Confidence( prob, "posterior probability", sd ) );
                }
                else {
                    node_to_annotate.getBranchData().addConfidence( new Confidence( prob, "posterior probability" ) );
                }
            }
        }
        final Matcher mb_bl_matcher = MB_BL_PATTERN.matcher( s );
        if ( mb_bl_matcher.find() ) {
            double bl = -1;
            try {
                bl = Double.parseDouble( mb_bl_matcher.group( 1 ) );
            }
            catch ( final NumberFormatException e ) {
                throw new NHXFormatException( "failed to parse median branch length (Mr Bayes output) from \"" + s
                        + "\"" );
            }
            if ( bl >= 0.0 ) {
                node_to_annotate.setDistanceToParent( bl );
            }
        }
    }

    public static enum TAXONOMY_EXTRACTION {
        AGGRESSIVE, NO, PFAM_STYLE_RELAXED, PFAM_STYLE_STRICT;
    }
}
