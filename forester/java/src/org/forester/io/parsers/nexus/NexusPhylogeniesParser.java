// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

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

import org.forester.archaeopteryx.AptxConstants;
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
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public final class NexusPhylogeniesParser implements IteratingPhylogenyParser, PhylogenyParser {

   
    final private static boolean DEBUG                               = false;
    
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
    final private static Pattern           RESIDUES_PATTERN          = Pattern.compile( "[A-Za-z\\-_\\*\\?\\.]+" );
    final private static Pattern           MATCHCHAR_PATTERN         = Pattern.compile( "matchchar\\s?.\\s?(\\S)" );
    final private static String            matrix                    = NexusConstants.MATRIX.toLowerCase();
    final private static Pattern           DATATYPE_PATTERN          = Pattern.compile( "datatype\\s?.\\s?([a-z]+)" );
    //final private static Pattern           LINK_TAXA_PATTERN         = Pattern.compile( "link\\s+taxa\\s?.\\s?([^;]+)",
    //                                                                                    Pattern.CASE_INSENSITIVE );
    final private static String            utree                     = NexusConstants.UTREE.toLowerCase();
    private BufferedReader                 _br;
    private boolean                        _ignore_quotes_in_nh_data = AptxConstants.NH_PARSING_IGNORE_QUOTES_DEFAULT;
    private boolean                        _in_taxalabels;
    private boolean                        _in_translate;
    private boolean                        _in_tree;
    private boolean                        _in_trees_block;
    private boolean                        _in_data_block;
    private boolean                        _in_matrix;
    private boolean                        _in_data_comment;
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
    private char                           _matchchar;
    private String                         _matrix_reference_id;
    private final boolean                  _add_sequences            = true;
    private boolean                       _parse_beast_style_extended_tags           = false;
           

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
        _in_matrix = false;
        _in_data_block = false;
        _in_data_comment = false;
        _datatype = null;
        _rooted_info_present = false;
        _is_rooted = false;
        _matchchar = 0;
        _matrix_reference_id = null;
        _seqs = new HashMap<String, MolecularSequence>();
        _br = ParserUtils.createReader( _nexus_source, ForesterConstants.UTF_8 );
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
        pars.setParseBeastStyleExtendedTags( _parse_beast_style_extended_tags );
        if ( rooted_info_present ) {
            pars.setGuessRootedness( false );
        }
        pars.setSource( nhx.toString() );
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
        if ( ( _taxlabels.size() > 0 ) || ( _translate_map.size() > 0 )
                || ( _add_sequences && !_seqs.isEmpty() ) ) {
            // Match sequence rows to tips by a CANONICAL key (see joinKey): a matrix often capitalizes
            // taxon names differently from the tree, Nexus treats '_' and ' ' as equivalent (and the
            // tree tip may already be underscore-replaced while the matrix id keeps its '_'), and labels
            // may be quoted. Last-wins on a key collision (taxa should not collide under this key).
            final Map<String, MolecularSequence> seqs_by_key = new HashMap<String, MolecularSequence>();
            for( final String seq_id : _seqs.keySet() ) {
                seqs_by_key.put( joinKey( seq_id ), _seqs.get( seq_id ) );
            }
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
                if ( _add_sequences && ( node.getName() != null ) ) {
                    final MolecularSequence s = seqs_by_key.get( joinKey( node.getName() ) );
                    if ( s != null ) {
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
            if ( DEBUG ) {
                System.out.println( line );
            }
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
                    _in_matrix = false;
                    _in_data_comment = false;
                    _datatype = null;
                    _matchchar = 0;
                    _matrix_reference_id = null;
                    // Scope the sequence rows to THIS matrix block, so an interleaved matrix
                    // concatenates within its own block and a later DATA/CHARACTERS block does
                    // not cross-contaminate an earlier one (multi-block Nexus).
                    _seqs = new HashMap<String, MolecularSequence>();
                }
                else if ( _in_trees_block ) {
                    if ( line_lc.startsWith( "title" ) ) {
                        final Matcher title_m = TITLE_PATTERN.matcher( line );
                        if ( title_m.lookingAt() ) {
                            _title = title_m.group( 1 );
                        }
                    }
                    else if ( line_lc.startsWith( "link" ) ) {
                        //final Matcher link_m = LINK_TAXA_PATTERN.matcher( line );
                        //if ( link_m.lookingAt() ) {
                            //final String link = link_m.group( 1 );  //TODO why?
                       // }
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
                    // Strip Nexus [ ... ] comments here (inside the matrix they are pure comments;
                    // unlike the trees block, where [&R]/[&...] are semantic and must be kept).
                    line = stripNexusComments( line );
                    final String dlc = line.toLowerCase();
                    if ( line.length() == 0 ) {
                        // comment-only (or now-empty) line: ignore
                    }
                    else if ( dlc.startsWith( end ) || dlc.startsWith( endblock ) ) {
                        _in_data_block = false;
                        _in_matrix = false;
                        _datatype = null;
                    }
                    else if ( dlc.startsWith( "link " ) ) {
                        // a LINK sub-command (e.g. "LINK TAXA=...") -- ignore. Note the trailing space:
                        // a bare "link" prefix would swallow a taxon row whose name starts with "link".
                    }
                    else if ( !_in_matrix ) {
                        // Block header: DIMENSIONS / FORMAT / CHARLABELS / CHARSTATELABELS / OPTIONS / ...
                        // Read the datatype (and MATCHCHAR) off FORMAT; enter the matrix on the MATRIX
                        // keyword; ignore the rest (a sub-command ending in ';' must NOT be mistaken for
                        // the end of the block -- that dropped the whole matrix when a CHARLABELS line
                        // preceded it).
                        final Matcher datatype_matcher = DATATYPE_PATTERN.matcher( dlc );
                        if ( datatype_matcher.find() ) {
                            _datatype = datatype_matcher.group( 1 );
                        }
                        final Matcher matchchar_matcher = MATCHCHAR_PATTERN.matcher( dlc );
                        if ( matchchar_matcher.find() ) {
                            _matchchar = matchchar_matcher.group( 1 ).charAt( 0 );
                        }
                        if ( dlc.equals( matrix ) || dlc.startsWith( matrix + " " ) ) {
                            _in_matrix = true;
                            String after = line.substring( matrix.length() ).trim();
                            boolean matrix_ends = false;
                            if ( after.endsWith( ";" ) ) {
                                matrix_ends = true;
                                after = after.substring( 0, after.length() - 1 ).trim();
                            }
                            if ( after.length() > 0 ) {
                                addMatrixRow( after );
                            }
                            if ( matrix_ends ) {
                                _in_matrix = false;
                                _in_data_block = false;
                                _datatype = null;
                            }
                        }
                    }
                    else {
                        // Inside the MATRIX: one taxon row per line until the terminating ';' (or END).
                        boolean matrix_ends = false;
                        if ( line.endsWith( ";" ) ) {
                            matrix_ends = true;
                            line = line.substring( 0, line.length() - 1 ).trim();
                        }
                        if ( line.length() > 0 ) {
                            addMatrixRow( line );
                        }
                        if ( matrix_ends ) {
                            _in_matrix = false;
                            _in_data_block = false;
                            _datatype = null;
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

    // Parse one MATRIX row ("taxon residues...") into an aligned molecular sequence. The id is the
    // first whitespace-delimited token; the residues are everything after it with all internal
    // whitespace removed (Nexus/PHYLIP commonly break the residues into space-separated blocks).
    // Only protein/dna/rna matrices become sequences. In an interleaved matrix each id reappears in
    // a later block, so a repeated id is CONCATENATED onto its existing row (never overwritten).
    private final void addMatrixRow( final String row ) {
        if ( ( _datatype == null )
                || !( _datatype.equals( "protein" ) || _datatype.equals( "dna" ) || _datatype.equals( "rna" ) ) ) {
            return;
        }
        if ( row.length() == 0 ) {
            return;
        }
        final String id;
        final String rest;
        final char c0 = row.charAt( 0 );
        if ( ( c0 == '\'' ) || ( c0 == '"' ) ) {
            // a quoted taxon label may contain spaces: the id runs to the matching closing quote
            final int close = row.indexOf( c0, 1 );
            if ( close < 1 ) {
                return;
            }
            id = row.substring( 0, close + 1 );
            rest = row.substring( close + 1 );
        }
        else {
            final int sp = row.indexOf( ' ' );
            if ( sp < 1 ) {
                return;
            }
            id = row.substring( 0, sp );
            rest = row.substring( sp + 1 );
        }
        String block = rest.replaceAll( "\\s+", "" );
        if ( ( block.length() == 0 ) || !RESIDUES_PATTERN.matcher( block ).matches() ) {
            return;
        }
        // MATCHCHAR (e.g. '.') means "identical to the first taxon (the reference) at this position".
        // Resolve against the reference at the same ABSOLUTE positions -- the reference row is read
        // first (sequential) or ahead of the others in each block (interleaved), so it always covers
        // the positions being resolved. Without this, createDnaSequence would turn '.' into a gap.
        if ( ( _matchchar != 0 ) && ( _matrix_reference_id != null ) && !id.equals( _matrix_reference_id )
                && _seqs.containsKey( _matrix_reference_id ) ) {
            final String ref = _seqs.get( _matrix_reference_id ).getMolecularSequenceAsString();
            final int offset = _seqs.containsKey( id ) ? _seqs.get( id ).getLength() : 0;
            final StringBuilder sb = new StringBuilder( block.length() );
            for( int j = 0; j < block.length(); ++j ) {
                final char c = block.charAt( j );
                if ( ( c == _matchchar ) && ( ( offset + j ) < ref.length() ) ) {
                    sb.append( ref.charAt( offset + j ) );
                }
                else {
                    sb.append( c );
                }
            }
            block = sb.toString();
        }
        String seq = block;
        if ( _seqs.containsKey( id ) ) {
            seq = _seqs.get( id ).getMolecularSequenceAsString() + block;
        }
        final MolecularSequence s;
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
        if ( _matrix_reference_id == null ) {
            // The first taxon of the matrix is the MATCHCHAR reference.
            _matrix_reference_id = id;
        }
    }

    // Canonical key for joining a matrix row id to a tree tip name. Nexus treats '_' and ' ' as
    // equivalent and labels may be quoted, a matrix often capitalizes names differently from the tree,
    // and the tree tip may already have been underscore-replaced while the matrix id keeps its '_' --
    // so normalize underscores to spaces, strip quotes, trim, and lower-case both sides.
    private final static String joinKey( final String name ) {
        return name.replace( '_', ' ' ).replaceAll( "['\"]+", "" ).trim().toLowerCase();
    }

    // Strip Nexus [ ... ] comments, tracking an OPEN comment across lines (via _in_data_comment) so a
    // MULTI-LINE comment inside the matrix cannot leak prose as a spurious taxon row (which could even
    // become the MATCHCHAR reference and corrupt other rows). Single-level: a nested "[[" closes on the
    // first ']' (nested Nexus comments are rare). Called only inside the data block -- the trees block
    // keeps [&R]/[&...], which are semantic. A dangling '[' with no ']' safely swallows to end of line.
    private final String stripNexusComments( final String s ) {
        if ( !_in_data_comment && ( s.indexOf( '[' ) < 0 ) ) {
            return s;
        }
        final StringBuilder out = new StringBuilder( s.length() );
        for( int i = 0; i < s.length(); ++i ) {
            final char c = s.charAt( i );
            if ( _in_data_comment ) {
                if ( c == ']' ) {
                    _in_data_comment = false;
                }
            }
            else if ( c == '[' ) {
                _in_data_comment = true;
            }
            else {
                out.append( c );
            }
        }
        return out.toString().trim();
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
    
    public final void setParseBeastStyleExtendedTags( final boolean parse_beast_style_extended_tags ) {
        _parse_beast_style_extended_tags = parse_beast_style_extended_tags;
    }
    
    private final static String removeWhiteSpaceBeforeSemicolon( final String s ) {
        return s.replaceAll( "\\s+;", ";" );
    }
}
