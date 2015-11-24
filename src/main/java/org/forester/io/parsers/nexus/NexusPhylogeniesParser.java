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

package org.forester.io.parsers.nexus;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.archaeopteryx.Constants;
import org.forester.io.parsers.IteratingPhylogenyParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nhx.NHXFormatException;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.parsers.util.PhylogenyParserException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterUtil;

public final class NexusPhylogeniesParser implements IteratingPhylogenyParser, PhylogenyParser {

    final private static String            begin_trees               = NexusConstants.BEGIN_TREES.toLowerCase();
    final private static String            end                       = NexusConstants.END.toLowerCase();
    final private static String            endblock                  = "endblock";
    final private static Pattern           ROOTEDNESS_PATTERN        = Pattern.compile( ".+=\\s*\\[&([R|U])\\].*" );
    final private static String            taxlabels                 = NexusConstants.TAXLABELS.toLowerCase();
    final private static Pattern           TITLE_PATTERN             = Pattern.compile( "TITLE.?\\s+([^;]+)",
                                                                                        Pattern.CASE_INSENSITIVE );
    final private static String            translate                 = NexusConstants.TRANSLATE.toLowerCase();
    final private static String            data                      = NexusConstants.BEGIN_CHARACTERS.toLowerCase();
    final private static String            characters                = NexusConstants.BEGIN_DATA.toLowerCase();
    final private static String            tree                      = NexusConstants.TREE.toLowerCase();
    final private static Pattern           TREE_NAME_PATTERN         = Pattern.compile( "\\s*.?Tree\\s+(.+?)\\s*=.+",
                                                                                        Pattern.CASE_INSENSITIVE );
    final private static Pattern           TRANSLATE_PATTERN         = Pattern.compile( "([0-9A-Za-z]+)\\s+(.+)" );
    final private static Pattern           ALN_PATTERN               = Pattern.compile( "(.+)\\s+([A-Za-z-_\\*\\?]+)" );
    final private static Pattern           DATATYPE_PATTERN          = Pattern.compile( "datatype\\s?.\\s?([a-z]+)" );
    final private static Pattern           LINK_TAXA_PATTERN         = Pattern.compile( "link\\s+taxa\\s?.\\s?([^;]+)",
                                                                                        Pattern.CASE_INSENSITIVE );
    final private static String            utree                     = NexusConstants.UTREE.toLowerCase();
    private BufferedReader                 _br;
    private boolean                        _ignore_quotes_in_nh_data = Constants.NH_PARSING_IGNORE_QUOTES_DEFAULT;
    private boolean                        _in_taxalabels;
    private boolean                        _in_translate;
    private boolean                        _in_tree;
    private boolean                        _in_trees_block;
    private boolean                        _in_data_block;
    private boolean                        _is_rooted;
    private String                         _datatype;
    private String                         _name;
    private Phylogeny                      _next;
    private Object                         _nexus_source;
    private StringBuilder                  _nh;
    private boolean                        _replace_underscores      = NHXParser.REPLACE_UNDERSCORES_DEFAULT;
    private boolean                        _rooted_info_present;
    private List<String>                   _taxlabels;
    private TAXONOMY_EXTRACTION            _taxonomy_extraction      = TAXONOMY_EXTRACTION.NO;
    private String                         _title;
    private Map<String, String>            _translate_map;
    private StringBuilder                  _translate_sb;
    private Map<String, MolecularSequence> _seqs;
    private final boolean                  _add_sequences            = true;

    @Override
    public String getName() {
        return "Nexus Phylogenies Parser";
    }

    @Override
    public final boolean hasNext() {
        return _next != null;
    }

    @Override
    public final Phylogeny next() throws NHXFormatException, IOException {
        final Phylogeny phy = _next;
        getNext();
        return phy;
    }

    @Override
    public final Phylogeny[] parse() throws IOException {
        final List<Phylogeny> l = new ArrayList<Phylogeny>();
        while ( hasNext() ) {
            l.add( next() );
        }
        final Phylogeny[] p = new Phylogeny[ l.size() ];
        for( int i = 0; i < l.size(); ++i ) {
            p[ i ] = l.get( i );
        }
        reset();
        return p;
    }

    @Override
    public final void reset() throws FileNotFoundException, IOException {
        _taxlabels = new ArrayList<String>();
        _translate_map = new HashMap<String, String>();
        _nh = new StringBuilder();
        _name = "";
        _title = "";
        _translate_sb = null;
        _next = null;
        _in_trees_block = false;
        _in_taxalabels = false;
        _in_translate = false;
        _in_tree = false;
        _rooted_info_present = false;
        _is_rooted = false;
        _seqs = new HashMap<String, MolecularSequence>();
        _br = ParserUtils.createReader( _nexus_source );
        getNext();
    }

    public final void setIgnoreQuotes( final boolean ignore_quotes_in_nh_data ) {
        _ignore_quotes_in_nh_data = ignore_quotes_in_nh_data;
    }

    public final void setReplaceUnderscores( final boolean replace_underscores ) {
        _replace_underscores = replace_underscores;
    }

    @Override
    public final void setSource( final Object nexus_source ) throws PhylogenyParserException, IOException {
        if ( nexus_source == null ) {
            throw new PhylogenyParserException( "attempt to parse null object" );
        }
        _nexus_source = nexus_source;
        reset();
    }

    public final void setTaxonomyExtraction( final TAXONOMY_EXTRACTION taxonomy_extraction ) {
        _taxonomy_extraction = taxonomy_extraction;
    }

    private final void createPhylogeny( final String title,
                                        final String name,
                                        final StringBuilder nhx,
                                        final boolean rooted_info_present,
                                        final boolean is_rooted ) throws IOException {
        _next = null;
        final NHXParser pars = new NHXParser();
        pars.setTaxonomyExtraction( _taxonomy_extraction );
        pars.setReplaceUnderscores( _replace_underscores );
        pars.setIgnoreQuotes( _ignore_quotes_in_nh_data );
        if ( rooted_info_present ) {
            pars.setGuessRootedness( false );
        }
        pars.setSource( nhx );
        final Phylogeny p = pars.next();
        if ( p == null ) {
            throw new PhylogenyParserException( "failed to create phylogeny" );
        }
        String myname = null;
        if ( !ForesterUtil.isEmpty( title ) && !ForesterUtil.isEmpty( name ) ) {
            myname = title.replace( '_', ' ' ).trim() + " (" + name.trim() + ")";
        }
        else if ( !ForesterUtil.isEmpty( title ) ) {
            myname = title.replace( '_', ' ' ).trim();
        }
        else if ( !ForesterUtil.isEmpty( name ) ) {
            myname = name.trim();
        }
        if ( !ForesterUtil.isEmpty( myname ) ) {
            p.setName( myname );
        }
        if ( rooted_info_present ) {
            p.setRooted( is_rooted );
        }
        if ( ( _taxlabels.size() > 0 ) || ( _translate_map.size() > 0 ) ) {
            final PhylogenyNodeIterator it = p.iteratorExternalForward();
            while ( it.hasNext() ) {
                final PhylogenyNode node = it.next();
                if ( ( _translate_map.size() > 0 ) && _translate_map.containsKey( node.getName() ) ) {
                    node.setName( _translate_map.get( node.getName() ).replaceAll( "['\"]+", "" ) );
                }
                else if ( _taxlabels.size() > 0 ) {
                    int i = -1;
                    try {
                        i = Integer.parseInt( node.getName() );
                    }
                    catch ( final NumberFormatException e ) {
                        // Ignore.
                    }
                    if ( i > 0 ) {
                        node.setName( _taxlabels.get( i - 1 ).replaceAll( "['\"]+", "" ) );
                    }
                }
                if ( !_replace_underscores && ( ( _taxonomy_extraction != TAXONOMY_EXTRACTION.NO ) ) ) {
                    ParserUtils.extractTaxonomyDataFromNodeName( node, _taxonomy_extraction );
                }
                else if ( _replace_underscores ) {
                    if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                        node.setName( node.getName().replace( '_', ' ' ).trim() );
                    }
                }
                if ( _add_sequences ) {
                    if ( _seqs.containsKey( node.getName() ) ) {
                        final MolecularSequence s = _seqs.get( node.getName() );
                        //TODO need to check for uniqueness when adding seqs....
                        final Sequence ns = new Sequence( s );
                        ns.setMolecularSequenceAligned( true ); //TODO need to check if all same length
                        node.getNodeData().addSequence( ns );
                    }
                }
            }
        }
        _next = p;
    }

    private final void getNext() throws IOException, NHXFormatException {
        _next = null;
        String line;
        while ( ( line = _br.readLine() ) != null ) {
            line = line.trim();
            if ( ( line.length() > 0 ) && !line.startsWith( "#" ) && !line.startsWith( ">" ) ) {
                line = ForesterUtil.collapseWhiteSpace( line );
                line = removeWhiteSpaceBeforeSemicolon( line );
                final String line_lc = line.toLowerCase();
                if ( line_lc.startsWith( begin_trees ) ) {
                    _in_trees_block = true;
                    _in_taxalabels = false;
                    _in_translate = false;
                    _in_data_block = false;
                    _datatype = null;
                    _title = "";
                }
                else if ( line_lc.startsWith( taxlabels ) ) {
                    //TODO need to be taxa block instead
                    _in_trees_block = false;
                    _in_taxalabels = true;
                    _in_translate = false;
                    _in_data_block = false;
                    _datatype = null;
                }
                else if ( line_lc.startsWith( translate ) ) {
                    _translate_sb = new StringBuilder();
                    _in_taxalabels = false;
                    _in_translate = true;
                    _in_data_block = false;
                    _datatype = null;
                }
                else if ( line_lc.startsWith( characters ) || line_lc.startsWith( data ) ) {
                    _in_taxalabels = false;
                    _in_trees_block = false;
                    _in_translate = false;
                    _in_data_block = true;
                    _datatype = null;
                }
                else if ( _in_trees_block ) {
                    if ( line_lc.startsWith( "title" ) ) {
                        final Matcher title_m = TITLE_PATTERN.matcher( line );
                        if ( title_m.lookingAt() ) {
                            _title = title_m.group( 1 );
                        }
                    }
                    else if ( line_lc.startsWith( "link" ) ) {
                        final Matcher link_m = LINK_TAXA_PATTERN.matcher( line );
                        if ( link_m.lookingAt() ) {
                            final String link = link_m.group( 1 );
                            //System.out.println( "link taxa:" + link );
                        }
                    }
                    else if ( line_lc.startsWith( end ) || line_lc.startsWith( endblock ) ) {
                        _in_trees_block = false;
                        _in_tree = false;
                        _in_translate = false;
                        if ( _nh.length() > 0 ) {
                            createPhylogeny( _title, _name, _nh, _rooted_info_present, _is_rooted );
                            _nh = new StringBuilder();
                            _name = "";
                            _rooted_info_present = false;
                            _is_rooted = false;
                            if ( _next != null ) {
                                return;
                            }
                        }
                    }
                    else if ( line_lc.startsWith( tree ) || ( line_lc.startsWith( utree ) ) ) {
                        boolean might = false;
                        if ( _nh.length() > 0 ) {
                            might = true;
                            createPhylogeny( _title, _name, _nh, _rooted_info_present, _is_rooted );
                            _nh = new StringBuilder();
                            _name = "";
                            _rooted_info_present = false;
                            _is_rooted = false;
                        }
                        _in_tree = true;
                        _nh.append( line.substring( line.indexOf( '=' ) ) );
                        final Matcher name_matcher = TREE_NAME_PATTERN.matcher( line );
                        if ( name_matcher.matches() ) {
                            _name = name_matcher.group( 1 );
                            _name = _name.replaceAll( "['\"]+", "" );
                        }
                        final Matcher rootedness_matcher = ROOTEDNESS_PATTERN.matcher( line );
                        if ( rootedness_matcher.matches() ) {
                            final String s = rootedness_matcher.group( 1 );
                            line = line.replaceAll( "\\[\\&.\\]", "" );
                            _rooted_info_present = true;
                            if ( s.toUpperCase().equals( "R" ) ) {
                                _is_rooted = true;
                            }
                        }
                        if ( might && ( _next != null ) ) {
                            return;
                        }
                    }
                    else if ( _in_tree && !_in_translate ) {
                        _nh.append( line );
                    }
                    if ( !line_lc.startsWith( "title" ) && !line_lc.startsWith( "link" ) && !_in_translate
                            && !line_lc.startsWith( end ) && !line_lc.startsWith( endblock ) && line_lc.endsWith( ";" ) ) {
                        _in_tree = false;
                        _in_translate = false;
                        createPhylogeny( _title, _name, _nh, _rooted_info_present, _is_rooted );
                        _nh = new StringBuilder();
                        _name = "";
                        _rooted_info_present = false;
                        _is_rooted = false;
                        if ( _next != null ) {
                            return;
                        }
                    }
                }
                if ( _in_taxalabels ) {
                    if ( line_lc.startsWith( end ) || line_lc.startsWith( endblock ) ) {
                        _in_taxalabels = false;
                    }
                    else {
                        final String[] labels = line.split( "\\s+" );
                        for( String label : labels ) {
                            if ( !label.toLowerCase().equals( taxlabels ) ) {
                                if ( label.endsWith( ";" ) ) {
                                    _in_taxalabels = false;
                                    label = label.substring( 0, label.length() - 1 );
                                }
                                if ( label.length() > 0 ) {
                                    _taxlabels.add( label );
                                }
                            }
                        }
                    }
                }
                if ( _in_translate ) {
                    if ( line_lc.startsWith( end ) || line_lc.startsWith( endblock ) ) {
                        _in_translate = false;
                    }
                    else {
                        _translate_sb.append( " " );
                        _translate_sb.append( line.trim() );
                        if ( line.endsWith( ";" ) ) {
                            _in_translate = false;
                            setTranslateKeyValuePairs( _translate_sb );
                        }
                    }
                }
                if ( _in_data_block ) {
                    if ( line_lc.startsWith( end ) || line_lc.startsWith( endblock ) ) {
                        _in_data_block = false;
                        _datatype = null;
                    }
                    else if ( line_lc.startsWith( "link" ) ) {
                        final Matcher link_m = LINK_TAXA_PATTERN.matcher( line );
                        if ( link_m.lookingAt() ) {
                            final String link = link_m.group( 1 );
                            //System.out.println( "link taxa:" + link );
                        }
                    }
                    else {
                        final Matcher datatype_matcher = DATATYPE_PATTERN.matcher( line_lc );
                        if ( datatype_matcher.find() ) {
                            _datatype = datatype_matcher.group( 1 );
                            //System.out.println( _datatype );
                        }
                        else {
                            if ( ( _datatype != null )
                                    && ( _datatype.equals( "protein" ) || _datatype.equals( "dna" ) || _datatype
                                            .equals( "rna" ) ) ) {
                                if ( line.endsWith( ";" ) ) {
                                    _in_data_block = false;
                                    line = line.substring( 0, line.length() - 1 );
                                }
                                final Matcher aln_matcher = ALN_PATTERN.matcher( line );
                                if ( aln_matcher.matches() ) {
                                    final String id = aln_matcher.group( 1 );
                                    final String seq = aln_matcher.group( 2 );
                                    MolecularSequence s = null;
                                    if ( _datatype.equals( "protein" ) ) {
                                        s = BasicSequence.createAaSequence( id, seq );
                                    }
                                    else if ( _datatype.equals( "dna" ) ) {
                                        s = BasicSequence.createDnaSequence( id, seq );
                                    }
                                    else {
                                        s = BasicSequence.createRnaSequence( id, seq );
                                    }
                                    _seqs.put( id, s );
                                    //System.out.println( s );
                                }
                            }
                        }
                    }
                }
            }
        }
        if ( _nh.length() > 0 ) {
            createPhylogeny( _title, _name, _nh, _rooted_info_present, _is_rooted );
            if ( _next != null ) {
                return;
            }
        }
    }

    private final void setTranslateKeyValuePairs( final StringBuilder translate_sb ) throws IOException {
        String s = translate_sb.toString().trim();
        if ( s.endsWith( ";" ) ) {
            s = s.substring( 0, s.length() - 1 ).trim();
        }
        for( String pair : s.split( "," ) ) {
            String key = "";
            String value = "";
            final int ti = pair.toLowerCase().indexOf( "translate" );
            if ( ti > -1 ) {
                pair = pair.substring( ti + 9 );
            }
            final Matcher m = TRANSLATE_PATTERN.matcher( pair );
            if ( m.find() ) {
                key = m.group( 1 );
                value = m.group( 2 ).replaceAll( "\'", "" ).replaceAll( "\"", "" ).trim();
            }
            else {
                throw new IOException( "ill-formatted translate values: " + pair );
            }
            if ( value.endsWith( ";" ) ) {
                value = value.substring( 0, value.length() - 1 );
            }
            _translate_map.put( key, value );
        }
    }

    private final static String removeWhiteSpaceBeforeSemicolon( final String s ) {
        return s.replaceAll( "\\s+;", ";" );
    }
}
