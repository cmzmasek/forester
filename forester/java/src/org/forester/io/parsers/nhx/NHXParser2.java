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

package org.forester.io.parsers.nhx;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.parsers.util.PhylogenyParserException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public final class NHXParser2 implements PhylogenyParser {

    public static final TAXONOMY_EXTRACTION TAXONOMY_EXTRACTION_DEFAULT = TAXONOMY_EXTRACTION.NO;
    final static private boolean            GUESS_ROOTEDNESS_DEFAULT    = true;
    final static private boolean            GUESS_IF_SUPPORT_VALUES     = true;
    final static private boolean            IGNORE_QUOTES_DEFAULT       = false;
    final static public boolean             REPLACE_UNDERSCORES_DEFAULT = false;
    private boolean                         _saw_closing_paren;
    final static private byte               STRING                      = 0;
    final static private byte               STRING_BUFFER               = 1;
    final static private byte               CHAR_ARRAY                  = 2;
    final static private byte               BUFFERED_READER             = 3;
    final static private byte               STRING_BUILDER              = 4;
    private boolean                         _guess_rootedness;
    private boolean                         _ignore_quotes;
    private byte                            _input_type;
    private int                             _source_length;
    private PhylogenyNode                   _current_node;
    private StringBuilder                   _current_anotation;
    private Object                          _nhx_source;
    private int                             _clade_level;
    private Phylogeny                       _current_phylogeny;
    private TAXONOMY_EXTRACTION             _taxonomy_extraction;
    private boolean                         _replace_underscores;
    public final static Pattern             UC_LETTERS_NUMBERS_PATTERN  = Pattern.compile( "^[A-Z0-9]+$" );
    public final static Pattern             NUMBERS_ONLY_PATTERN        = Pattern.compile( "^[0-9\\.]+$" );
    public final static Pattern             MB_PROB_PATTERN             = Pattern.compile( "prob=([^,]+)" );
    public final static Pattern             MB_PROB_SD_PATTERN          = Pattern.compile( "prob_stddev=([^,]+)" );
    public final static Pattern             MB_BL_PATTERN               = Pattern.compile( "length_median=([^,]+)" );
    boolean                                 _in_comment                 = false;
    boolean                                 _saw_colon                  = false;
    boolean                                 _saw_open_bracket           = false;
    boolean                                 _in_open_bracket            = false;
    boolean                                 _in_double_quote            = false;
    boolean                                 _in_single_quote            = false;
    String                                  _my_source_str              = null;
    StringBuffer                            _my_source_sbuff            = null;
    StringBuilder                           _my_source_sbuil            = null;
    char[]                                  _my_source_charary          = null;
    BufferedReader                          _my_source_br               = null;
    int                                     _i;
    private Phylogeny                       _next;
    private Object                          _source;

    public NHXParser2() {
        init();
    }

    public TAXONOMY_EXTRACTION getTaxonomyExtraction() {
        return _taxonomy_extraction;
    }

    public boolean hasNext() {
        return _next != null;
    }

    public Phylogeny next() throws NHXFormatException, IOException {
        final Phylogeny phy = _next;
        getNext();
        return phy;
    }

    @Override
    public Phylogeny[] parse() throws IOException {
        reset();
        List<Phylogeny> l = new ArrayList<Phylogeny>();
        System.out.println( ">> _next=" + _next );
        while ( hasNext() ) {
            Phylogeny n = next();
            System.out.println( ">> going to add " + n );
            l.add( n );
        }
        final Phylogeny[] p = new Phylogeny[ l.size() ];
        for( int i = 0; i < l.size(); ++i ) {
            p[ i ] = l.get( i );
        }
        return p;
    }

    public void reset() throws NHXFormatException, IOException {
        _i = 0;
        _next = null;
        _in_comment = false;
        _saw_colon = false;
        _saw_open_bracket = false;
        _in_open_bracket = false;
        _in_double_quote = false;
        _in_single_quote = false;
        setCladeLevel( 0 );
        newCurrentAnotation();
        setCurrentPhylogeny( null );
        setCurrentNode( null );
        _my_source_str = null;
        _my_source_sbuff = null;
        _my_source_sbuil = null;
        _my_source_charary = null;
        _my_source_br = null;
        determineSourceType( _source );
        switch ( getInputType() ) {
            case STRING:
                _my_source_str = ( String ) getNhxSource();
                break;
            case STRING_BUFFER:
                _my_source_sbuff = ( StringBuffer ) getNhxSource();
                break;
            case STRING_BUILDER:
                _my_source_sbuil = ( StringBuilder ) getNhxSource();
                break;
            case CHAR_ARRAY:
                _my_source_charary = ( char[] ) getNhxSource();
                break;
            case BUFFERED_READER:
                if ( _my_source_br != null ) {
                    try {
                        _my_source_br.close();
                    }
                    catch ( IOException e ) {
                        //do nothing
                    }
                }
                _my_source_br = ( BufferedReader ) getNhxSource();
                break;
            default:
                throw new RuntimeException( "unknown input type" );
        }
        getNext();
    }

    public void setGuessRootedness( final boolean guess_rootedness ) {
        _guess_rootedness = guess_rootedness;
    }

    public void setIgnoreQuotes( final boolean ignore_quotes ) {
        _ignore_quotes = ignore_quotes;
    }

    public void setReplaceUnderscores( final boolean replace_underscores ) {
        _replace_underscores = replace_underscores;
    }

    /**
     * This sets the source to be parsed. The source can be: String,
     * StringBuffer, char[], File, or InputStream. The source can contain more
     * than one phylogenies in either New Hamphshire (NH) or New Hamphshire
     * Extended (NHX) format. There is no need to separate phylogenies with any
     * special character. White space is always ignored, as are semicolons
     * inbetween phylogenies. Example of a source describing two phylogenies
     * (source is a String, in this example): "(A,(B,(C,(D,E)de)cde)bcde)abcde
     * ((((A,B)ab,C)abc,D)abcd,E)abcde". Everything between a '[' followed by any
     * character other than '&' and ']' is considered a comment and ignored
     * (example: "[this is a comment]"). NHX tags are surrounded by '[&&NHX' and
     * ']' (example: "[&&NHX:S=Varanus_storri]"). A sequence like "[& some
     * info]" is ignored, too (at the PhylogenyNode level, though).
     * Exception: numbers only between [ and ] (e.g. [90]) are interpreted as support values.
     * 
     * @see #parse()
     * @see org.forester.io.parsers.PhylogenyParser#setSource(java.lang.Object)
     * @param nhx_source
     *            the source to be parsed (String, StringBuffer, char[], File,
     *            or InputStream)
     * @throws NHXFormatException 
     * @throws IOException
     * @throws PhylogenyParserException
     */
    @Override
    public void setSource( final Object nhx_source ) throws NHXFormatException, IOException {
        _source = nhx_source;
        reset();
    }

    private void determineSourceType( final Object nhx_source ) throws PhylogenyParserException, FileNotFoundException {
        if ( nhx_source == null ) {
            throw new PhylogenyParserException( getClass() + ": attempt to parse null object." );
        }
        else if ( nhx_source instanceof String ) {
            setInputType( NHXParser2.STRING );
            setSourceLength( ( ( String ) nhx_source ).length() );
            setNhxSource( nhx_source );
        }
        else if ( nhx_source instanceof StringBuilder ) {
            setInputType( NHXParser2.STRING_BUILDER );
            setSourceLength( ( ( StringBuilder ) nhx_source ).length() );
            setNhxSource( nhx_source );
        }
        else if ( nhx_source instanceof StringBuffer ) {
            setInputType( NHXParser2.STRING_BUFFER );
            setSourceLength( ( ( StringBuffer ) nhx_source ).length() );
            setNhxSource( nhx_source );
        }
        else if ( nhx_source instanceof StringBuilder ) {
            setInputType( NHXParser2.STRING_BUILDER );
            setSourceLength( ( ( StringBuilder ) nhx_source ).length() );
            setNhxSource( nhx_source );
        }
        else if ( nhx_source instanceof char[] ) {
            setInputType( NHXParser2.CHAR_ARRAY );
            setSourceLength( ( ( char[] ) nhx_source ).length );
            setNhxSource( nhx_source );
        }
        else if ( nhx_source instanceof File ) {
            setInputType( NHXParser2.BUFFERED_READER );
            setSourceLength( 0 );
            final File f = ( File ) nhx_source;
            final String error = ForesterUtil.isReadableFile( f );
            if ( !ForesterUtil.isEmpty( error ) ) {
                throw new PhylogenyParserException( error );
            }
            setNhxSource( new BufferedReader( new FileReader( f ) ) );
        }
        else if ( nhx_source instanceof InputStream ) {
            setInputType( NHXParser2.BUFFERED_READER );
            setSourceLength( 0 );
            final InputStreamReader isr = new InputStreamReader( ( InputStream ) nhx_source );
            setNhxSource( new BufferedReader( isr ) );
        }
        else {
            throw new IllegalArgumentException( getClass() + " can only parse objects of type String,"
                    + " StringBuffer, StringBuilder, char[], File," + " or InputStream "
                    + " [attempt to parse object of " + nhx_source.getClass() + "]." );
        }
    }

    public void setTaxonomyExtraction( final TAXONOMY_EXTRACTION taxonomy_extraction ) {
        _taxonomy_extraction = taxonomy_extraction;
    }

    /**
     * Decreases the clade level by one.
     * 
     * @throws PhylogenyParserException
     *             if level goes below zero.
     */
    private void decreaseCladeLevel() throws PhylogenyParserException {
        if ( getCladeLevel() < 0 ) {
            throw new PhylogenyParserException( "error in NH (Newick)/NHX formatted data: most likely cause: number of close parens is larger than number of open parens" );
        }
        --_clade_level;
    }

    private Phylogeny finishPhylogeny2() throws PhylogenyParserException, NHXFormatException,
            PhyloXmlDataFormatException {
        //setCladeLevel( 0 );
        if ( getCurrentPhylogeny() != null ) {
            System.out.println( "fp: cp=" + getCurrentPhylogeny() );
            if ( getCurrentAnotation() != null ) {
                System.out.println( "fp: ca=" + getCurrentAnotation().toString() );
            }
            else {
                System.out.println( "fp: ca=null" );
            }
            parseNHX( getCurrentAnotation() != null ? getCurrentAnotation().toString() : "", getCurrentPhylogeny()
                    .getRoot(), getTaxonomyExtraction(), isReplaceUnderscores() );
            if ( GUESS_IF_SUPPORT_VALUES ) {
                if ( isBranchLengthsLikeBootstrapValues( getCurrentPhylogeny() ) ) {
                    moveBranchLengthsToConfidenceValues( getCurrentPhylogeny() );
                }
            }
            if ( isGuessRootedness() ) {
                final PhylogenyNode root = getCurrentPhylogeny().getRoot();
                if ( ( root.getDistanceToParent() >= 0.0 ) || !ForesterUtil.isEmpty( root.getName() )
                        || !ForesterUtil.isEmpty( PhylogenyMethods.getSpecies( root ) ) || root.isHasAssignedEvent() ) {
                    getCurrentPhylogeny().setRooted( true );
                }
            }
            return getCurrentPhylogeny();
        }
        return null;
    }

    private Phylogeny finishSingleNodePhylogeny2() throws PhylogenyParserException, NHXFormatException,
            PhyloXmlDataFormatException {
        // setCladeLevel( 0 );
        final PhylogenyNode new_node = new PhylogenyNode();
        parseNHX( getCurrentAnotation().toString(), new_node, getTaxonomyExtraction(), isReplaceUnderscores() );
        setCurrentPhylogeny( new Phylogeny() );
        getCurrentPhylogeny().setRoot( new_node );
        return getCurrentPhylogeny();
    }

    private int getCladeLevel() {
        return _clade_level;
    }

    private StringBuilder getCurrentAnotation() {
        return _current_anotation;
    }

    private PhylogenyNode getCurrentNode() {
        return _current_node;
    }

    private Phylogeny getCurrentPhylogeny() {
        return _current_phylogeny;
    }

    private byte getInputType() {
        return _input_type;
    }

    private void getNext() throws IOException, NHXFormatException {
        while ( true ) {
            char c = '\b';
            if ( getInputType() == BUFFERED_READER ) {
                final int ci = _my_source_br.read();
                if ( ci >= 0 ) {
                    c = ( char ) ci;
                }
                else {
                    break;
                }
            }
            else {
                if ( _i >= getSourceLength() ) {
                    break;
                }
                else {
                    switch ( getInputType() ) {
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
            if ( ( isIgnoreQuotes() && ( ( c < 33 ) || ( c > 126 ) || ( c == 34 ) || ( c == 39 ) || ( ( getCladeLevel() == 0 ) && ( c == ';' ) ) ) )
                    || ( !isIgnoreQuotes() && ( ( c < 32 ) || ( c > 126 ) || ( ( getCladeLevel() == 0 ) && ( c == ';' ) ) ) ) ) {
                // Do nothing.
            }
            else if ( ( c == 32 ) && ( !_in_single_quote && !_in_double_quote ) ) {
                // Do nothing.
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
                    getCurrentAnotation().append( c );
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
                    getCurrentAnotation().append( c );
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
                        getCurrentAnotation().append( "[&" );
                    }
                    else if ( _saw_colon ) {
                        getCurrentAnotation().append( "[" + c );
                    }
                    else {
                        _in_comment = true;
                    }
                }
                // comment consisting just of "[]":
                _saw_open_bracket = false;
            }
            else if ( ( c == '(' ) && !_in_open_bracket ) {
                final Phylogeny phy = processOpenParen2();
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
                getCurrentAnotation().append( c );
            }
            ++_i;
        } //  while ( true ) 
        System.out.println( "done with loop" );
        if ( getCurrentPhylogeny() == null ) {
            System.out.println( "... but is null" );
        }
        if ( getCladeLevel() != 0 ) {
            throw new PhylogenyParserException( "error in NH (Newick) formatted data: most likely cause: number of open parens does not equal number of close parens" );
        }
        if ( getCurrentPhylogeny() != null ) {
            System.out.println( "... and current=" + getCurrentPhylogeny() );
            _next = finishPhylogeny2();
            System.out.println( "... _next=" + _next );
            setCurrentPhylogeny( null );
            setCurrentAnotation( null );
            //return finishPhylogeny2();
        }
        else if ( ( getCurrentAnotation() != null ) && ( getCurrentAnotation().length() > 0 ) ) {
            System.out.println( "1node=" + getCurrentAnotation() );
            _next = finishSingleNodePhylogeny2();
            setCurrentAnotation( null );
            //return finishSingleNodePhylogeny2();
        }
        else {
            _next = null;
            //return null;
        }
    }

    private Object getNhxSource() {
        return _nhx_source;
    }

    private int getSourceLength() {
        return _source_length;
    }

    private void increaseCladeLevel() {
        ++_clade_level;
    }

    private void init() {
        setTaxonomyExtraction( TAXONOMY_EXTRACTION_DEFAULT );
        setReplaceUnderscores( REPLACE_UNDERSCORES_DEFAULT );
        setGuessRootedness( GUESS_ROOTEDNESS_DEFAULT );
        setIgnoreQuotes( IGNORE_QUOTES_DEFAULT );
    }

    private boolean isGuessRootedness() {
        return _guess_rootedness;
    }

    private boolean isIgnoreQuotes() {
        return _ignore_quotes;
    }

    private boolean isReplaceUnderscores() {
        return _replace_underscores;
    }

    private boolean isSawClosingParen() {
        return _saw_closing_paren;
    }

    /**
     * Replaces the current annotation with a new StringBuffer.
     */
    private void newCurrentAnotation() {
        setCurrentAnotation( new StringBuilder() );
    }

    /**
     * Called if a closing paren is encountered.
     * 
     * @throws PhylogenyParserException
     * @throws NHXFormatException
     * @throws PhyloXmlDataFormatException 
     */
    private void processCloseParen() throws PhylogenyParserException, NHXFormatException, PhyloXmlDataFormatException {
        decreaseCladeLevel();
        if ( !isSawClosingParen() ) {
            final PhylogenyNode new_node = new PhylogenyNode();
            parseNHX( getCurrentAnotation().toString(), new_node, getTaxonomyExtraction(), isReplaceUnderscores() );
            newCurrentAnotation();
            getCurrentNode().addAsChild( new_node );
        }
        else {
            parseNHX( getCurrentAnotation().toString(),
                      getCurrentNode().getLastChildNode(),
                      getTaxonomyExtraction(),
                      isReplaceUnderscores() );
            newCurrentAnotation();
        }
        if ( !getCurrentNode().isRoot() ) {
            setCurrentNode( getCurrentNode().getParent() );
        }
        setSawClosingParen( true );
    }

    /**
     * Called if a comma is encountered.
     * 
     * @throws PhylogenyParserException
     * @throws NHXFormatException
     * @throws PhyloXmlDataFormatException 
     */
    private void processComma() throws PhylogenyParserException, NHXFormatException, PhyloXmlDataFormatException {
        if ( !isSawClosingParen() ) {
            final PhylogenyNode new_node = new PhylogenyNode();
            parseNHX( getCurrentAnotation().toString(), new_node, getTaxonomyExtraction(), isReplaceUnderscores() );
            if ( getCurrentNode() == null ) {
                throw new NHXFormatException( "format might not be NH or NHX" );
            }
            getCurrentNode().addAsChild( new_node );
        }
        else {
            parseNHX( getCurrentAnotation().toString(),
                      getCurrentNode().getLastChildNode(),
                      getTaxonomyExtraction(),
                      isReplaceUnderscores() );
        }
        newCurrentAnotation();
        setSawClosingParen( false );
    }

    private Phylogeny processOpenParen2() throws PhylogenyParserException, NHXFormatException,
            PhyloXmlDataFormatException {
        Phylogeny phy = null;
        final PhylogenyNode new_node = new PhylogenyNode();
        System.out.println( "level=" + getCladeLevel() );
        if ( getCladeLevel() == 0 ) {
            if ( getCurrentPhylogeny() != null ) {
                phy = finishPhylogeny2();
            }
            setCladeLevel( 1 );
            newCurrentAnotation();
            setCurrentPhylogeny( new Phylogeny() );
            getCurrentPhylogeny().setRoot( new_node );
        }
        else {
            increaseCladeLevel();
            getCurrentNode().addAsChild( new_node );
        }
        setCurrentNode( new_node );
        setSawClosingParen( false );
        if ( phy != null ) {
            System.out.println( "processOpenParen2 returns " + phy.toString() );
        }
        else {
            System.out.println( "processOpenParen2 returns null" );
        }
        return phy;
    }

    private void setCladeLevel( final int clade_level ) {
        if ( clade_level < 0 ) {
            throw new IllegalArgumentException( "attempt to set clade level to a number smaller than zero" );
        }
        _clade_level = clade_level;
    }

    private void setCurrentAnotation( final StringBuilder current_anotation ) {
        _current_anotation = current_anotation;
    }

    private void setCurrentNode( final PhylogenyNode current_node ) {
        _current_node = current_node;
    }

    private void setCurrentPhylogeny( final Phylogeny current_phylogeny ) {
        _current_phylogeny = current_phylogeny;
    }

    private void setInputType( final byte input_type ) {
        _input_type = input_type;
    }

    private void setNhxSource( final Object nhx_source ) {
        _nhx_source = nhx_source;
    }

    private void setSawClosingParen( final boolean saw_closing_paren ) {
        _saw_closing_paren = saw_closing_paren;
    }

    private void setSourceLength( final int source_length ) {
        _source_length = source_length;
    }

    public static void parseNHX( String s,
                                 final PhylogenyNode node_to_annotate,
                                 final TAXONOMY_EXTRACTION taxonomy_extraction,
                                 final boolean replace_underscores ) throws NHXFormatException,
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
                    if ( s.startsWith( org.forester.io.parsers.nhx.NHXtags.SPECIES_NAME ) ) {
                        if ( !node_to_annotate.getNodeData().isHasTaxonomy() ) {
                            node_to_annotate.getNodeData().setTaxonomy( new Taxonomy() );
                        }
                        node_to_annotate.getNodeData().getTaxonomy().setScientificName( s.substring( 2 ) );
                    }
                    else if ( s.startsWith( org.forester.io.parsers.nhx.NHXtags.IS_DUPLICATION ) ) {
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
                        PhylogenyMethods.setConfidence( node_to_annotate, doubleValue( s.substring( 2 ) ) );
                    }
                    else if ( s.startsWith( NHXtags.TAXONOMY_ID ) ) {
                        if ( !node_to_annotate.getNodeData().isHasTaxonomy() ) {
                            node_to_annotate.getNodeData().setTaxonomy( new Taxonomy() );
                        }
                        node_to_annotate.getNodeData().getTaxonomy().setIdentifier( new Identifier( s.substring( 2 ) ) );
                    }
                    else if ( s.startsWith( NHXtags.PARENT_BRANCH_WIDTH ) ) {
                        PhylogenyMethods.setBranchWidthValue( node_to_annotate, Integer.parseInt( s.substring( 2 ) ) );
                    }
                    else if ( s.startsWith( NHXtags.COLOR ) ) {
                        final Color c = NHXParser2.stringToColor( s.substring( 2 ) );
                        if ( c != null ) {
                            PhylogenyMethods.setBranchColorValue( node_to_annotate, c );
                        }
                    }
                    else if ( s.startsWith( NHXtags.DOMAIN_STRUCTURE ) ) {
                        if ( !node_to_annotate.getNodeData().isHasSequence() ) {
                            node_to_annotate.getNodeData().setSequence( new Sequence() );
                        }
                        node_to_annotate.getNodeData().getSequence()
                                .setDomainArchitecture( new DomainArchitecture( s.substring( 3 ) ) );
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
                        if ( node_to_annotate.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                            throw new NHXFormatException( "error in NHX formatted data: more than one distance to parent:"
                                    + "\"" + s + "\"" );
                        }
                        node_to_annotate.setDistanceToParent( doubleValue( s ) );
                    }
                } // while ( t.hasMoreTokens() ) 
            }
        }
    }

    private static double doubleValue( final String str ) throws NHXFormatException {
        try {
            return Double.valueOf( str ).doubleValue();
        }
        catch ( final NumberFormatException ex ) {
            throw new NHXFormatException( "error in NH/NHX formatted data: failed to parse number from " + "\"" + str
                    + "\"" );
        }
    }

    private static boolean isBranchLengthsLikeBootstrapValues( final Phylogeny p ) {
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

    private static void moveBranchLengthsToConfidenceValues( final Phylogeny p ) {
        final PhylogenyNodeIterator it = p.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            PhylogenyMethods.setBootstrapConfidence( n, n.getDistanceToParent() );
            n.setDistanceToParent( PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
        }
    }

    private static void processMrBayes3Data( final String s, final PhylogenyNode node_to_annotate )
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

    /**
     * Parses String s in the format r.g.b (e.g. "12.34.234" ) into red, green,
     * and blue and returns the corresponding Color.
     */
    private static Color stringToColor( final String s ) {
        final StringTokenizer st = new StringTokenizer( s, "." );
        if ( st.countTokens() != 3 ) {
            throw new IllegalArgumentException( "illegal format for color: " + s );
        }
        final int red = ForesterUtil.limitRangeForColor( Integer.parseInt( st.nextToken() ) );
        final int green = ForesterUtil.limitRangeForColor( Integer.parseInt( st.nextToken() ) );
        final int blu = ForesterUtil.limitRangeForColor( Integer.parseInt( st.nextToken() ) );
        return new Color( red, green, blu );
    }
}
