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

package org.forester.io.parsers.nhx;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
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
import org.forester.phylogeny.data.*;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public final class NHXParser implements PhylogenyParser, IteratingPhylogenyParser {

    private final static Pattern MB_BL_PATTERN = Pattern
            .compile("length.median=([-+eE0-9\\.]+)");
    private final static Pattern MB_PROB_PATTERN = Pattern.compile("prob=([-+eE0-9\\.]+)");
    private final static Pattern MB_PROB_SD_PATTERN = Pattern
            .compile("prob.stddev=([-+eE0-9\\.]+)");
    // A branch length may be signed and may start with its decimal point: ":-0.1" (a negative length is what some
    // distance methods write, and a parser reports what the file says) and ":.5" (legal Newick) used to vanish.
    private final static Pattern BRANCH_LENGTH_START_PATTERN = Pattern.compile("^[-+]?(\\d|\\.\\d)");
    private final static Pattern NUMBERS_ONLY_PATTERN = Pattern.compile("^[-+]?[0-9\\.]+$");
    private final static Pattern ENDS_WITH_NUMBER_PATTERN = Pattern.compile("(:[-+eE0-9\\.]+$)");
    public final static boolean REPLACE_UNDERSCORES_DEFAULT = false;
    private final static boolean ALLOW_ERRORS_IN_DISTANCE_TO_PARENT_DEFAULT = false;
    private final static byte BUFFERED_READER = 3;
    private final static byte CHAR_ARRAY = 2;
    private final static boolean GUESS_IF_SUPPORT_VALUES = true;
    private final static boolean GUESS_ROOTEDNESS_DEFAULT = true;
    private final static boolean IGNORE_QUOTES_DEFAULT = false;
    private final static char BELL = 7;
    final static char BLOB_OPEN_BRACKET = 1;
    final static char BLOB_CLOSE_BRACKET = 2;
    private final static int BLOB_QUOTE_LOOKAHEAD = 1 << 16;
    private final static String ENCODING_DEFAULT = ForesterConstants.UTF_8;
    private boolean _allow_errors_in_distance_to_parent;
    private int _clade_level;
    private StringBuilder _current_anotation;
    private PhylogenyNode _current_node;
    private Phylogeny _current_phylogeny;
    private boolean _guess_rootedness;
    private int _i;
    private boolean _ignore_quotes;
    private boolean _in_comment = false;
    private boolean _in_double_quote = false;
    // The quote character that closed a quoted label on the PREVIOUS character, or 0.
    // An immediately following identical quote is the Nexus/Newick escape for a literal one.
    private char _just_closed_quote = 0;
    private boolean _in_open_bracket = false;
    // True only once a bracket has been committed to as a real "[&...]" extended annotation (the character right
    // after '[' was '&'). Narrower than _in_open_bracket on purpose: _in_open_bracket also covers a plain "[...]"
    // comment and the leading "[91]"-style numeric confidence, neither of which is written to _current_anotation
    // character-by-character the way an annotation blob is, and a comment's content is discarded regardless of
    // whitespace -- so exempting spaces/quotes for THOSE from the usual drop rule would only risk misreading the
    // single character right after '[' that decides comment-vs-annotation-vs-numeric in the first place (a literal
    // space there must still be swallowed so e.g. "[ 91 ]" is still recognised as the numeric confidence "91").
    private boolean _in_kept_annotation = false;
    // True once the second character of a kept annotation is also '&', i.e. the legacy "[&&NHX:key=value:...]"
    // tag scheme (as opposed to a single-'&' "[&key=value,...]" BEAST/FigTree/TreeTime/Auspice-style blob). NHX
    // tag values are programmatic identifiers (species/gene names, normally underscore-separated) and whitespace
    // in and around them has always been pure formatting noise -- testNHXParsingQuotes pins "mo\tnkey !" reducing
    // to "monkey!". A single-'&' blob's values are free text where whitespace (and quoting) is real data (Auspice
    // writes "country=Democratic Republic of the Congo" with no quotes at all). Set once, right after the '&&' is
    // seen, so it must be checked before deciding whether to keep or drop the very next character.
    private boolean _in_legacy_nhx_tag = false;
    // True for exactly one character: the one right after a kept annotation's opening '&', which decides
    // _in_legacy_nhx_tag. Consumption skips characters that are unconditionally dropped everywhere (control
    // characters), so e.g. the tab in "[&\t&NHX:...]" is not mistaken for the deciding character.
    private boolean _awaiting_second_annotation_char = false;
    // > 0 while inside a quoted VALUE of a kept blob: the characters left up to and including its closing quote
    private int _blob_quote_remaining = 0;
    private boolean _in_single_quote = false;
    private byte _input_type;
    private BufferedReader _my_source_br = null;
    private char[] _my_source_charary = null;
    private Phylogeny _next;
    private Object _nhx_source;
    private boolean _replace_underscores;
    private boolean _saw_closing_paren;
    private boolean _saw_colon = false;
    private boolean _saw_open_bracket = false;
    private boolean _after_close_paren = false;
    private Object _source;
    private int _source_length;
    private TAXONOMY_EXTRACTION _taxonomy_extraction;
    // ON by default: "[&...]" bracket annotations are what BEAST, MrBayes, FigTree, TreeTime and Auspice write, and
    // a library that silently drops posteriors, node ages and traits unless told otherwise is a trap (aptx_render
    // fell into it). Off, a blob is kept whole as an nh:comment, exactly as before.
    public final static boolean PARSE_BRACKET_ANNOTATIONS_DEFAULT = true;
    private boolean _parse_beast_style_extended_tags = PARSE_BRACKET_ANNOTATIONS_DEFAULT;
    private boolean _normalize_bracket_annotations = true;
    private final String _encoding;

    public NHXParser() {
        _encoding = ENCODING_DEFAULT;
        init();
    }

    public NHXParser(final String encoding) {
        _encoding = encoding;
        init();
    }

    @Override
    public String getName() {
        return "NH/NHX Parser";
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
        while (hasNext()) {
            l.add(next());
        }
        final Phylogeny[] p = new Phylogeny[l.size()];
        for (int i = 0; i < l.size(); ++i) {
            p[i] = l.get(i);
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
        _in_kept_annotation = false;
        _in_legacy_nhx_tag = false;
        _awaiting_second_annotation_char = false;
        _blob_quote_remaining = 0;
        _in_double_quote = false;
        _in_single_quote = false;
        _just_closed_quote = 0;
        _after_close_paren = false;
        _clade_level = 0;
        _current_anotation = new StringBuilder();
        _current_phylogeny = null;
        _current_node = null;
        _my_source_charary = null;
        determineAndProcessSourceType(_source);
        switch (_input_type) {
            case CHAR_ARRAY:
                _my_source_br = null;
                _my_source_charary = (char[]) _nhx_source;
                break;
            case BUFFERED_READER:
                _my_source_br = (BufferedReader) _nhx_source;
                break;
            default:
                throw new RuntimeException("unknown input type");
        }
        parseNext();
    }

    public final void setGuessRootedness(final boolean guess_rootedness) {
        _guess_rootedness = guess_rootedness;
    }

    public final void setIgnoreQuotes(final boolean ignore_quotes) {
        _ignore_quotes = ignore_quotes;
    }

    public final void setReplaceUnderscores(final boolean replace_underscores) {
        _replace_underscores = replace_underscores;
    }

    @Override
    public final void setSource(final Object nhx_source) throws NHXFormatException, IOException {
        _source = nhx_source;
        reset();
    }

    public final void setTaxonomyExtraction(final TAXONOMY_EXTRACTION taxonomy_extraction) {
        _taxonomy_extraction = taxonomy_extraction;
    }

    public final void setAllowErrorsInDistanceToParent(final boolean allow_errors_in_distance_to_parent) {
        _allow_errors_in_distance_to_parent = allow_errors_in_distance_to_parent;
    }

    private final void determineAndProcessSourceType(final Object nhx_source) throws IOException {
        if (nhx_source == null) {
            throw new PhylogenyParserException(getClass() + ": attempt to parse null object.");
        } else if (nhx_source instanceof String) {
            _nhx_source = nhx_source;
            _input_type = NHXParser.BUFFERED_READER;
            _source_length = 0;
            final InputStream is = new ByteArrayInputStream(((String) nhx_source).getBytes(getEncoding()));
            final InputStreamReader isr = new InputStreamReader(is, getEncoding());
            _nhx_source = new BufferedReader(isr);
        } else if (nhx_source instanceof char[]) {
            _input_type = NHXParser.CHAR_ARRAY;
            _source_length = ((char[]) nhx_source).length;
            _nhx_source = nhx_source;
        } else if (nhx_source instanceof File) {
            _input_type = NHXParser.BUFFERED_READER;
            _source_length = 0;
            if (_my_source_br != null) {
                //I am REALLY not sure if it is a "good" idea NOT to close the stream...
                //                try {
                //                    _my_source_br.close();
                //                }
                //                catch ( final IOException e ) {
                //                }
            }
            final File f = (File) nhx_source;
            final String error = ForesterUtil.isReadableFile(f);
            if (!ForesterUtil.isEmpty(error)) {
                throw new PhylogenyParserException(error);
            }
            final InputStream is = new FileInputStream(f);
            final InputStreamReader isr = new InputStreamReader(is, getEncoding());
            _nhx_source = new BufferedReader(isr);
        } else if (nhx_source instanceof URL) {
            _input_type = NHXParser.BUFFERED_READER;
            _source_length = 0;
            if (_my_source_br != null) {
                //I am REALLY not sure if it is a "good" idea NOT to close the stream...
                //                try {
                //                    _my_source_br.close();
                //                }
                //                catch ( final IOException e ) {
                //                }
            }
            final InputStream is = ((URL) nhx_source).openStream();
            final InputStreamReader isr = new InputStreamReader(is, getEncoding());
            _nhx_source = new BufferedReader(isr);
        } else if (nhx_source instanceof InputStream) {
            _input_type = NHXParser.BUFFERED_READER;
            _source_length = 0;
            if (_my_source_br != null) {
                //I am REALLY not sure if it is a "good" idea NOT to close the stream...
                //                try {
                //                    _my_source_br.close();
                //                }
                //                catch ( final IOException e ) {
                //                }
            }
            final InputStream is = (InputStream) nhx_source;
            final InputStreamReader isr = new InputStreamReader(is, getEncoding());
            _nhx_source = new BufferedReader(isr);
        } else {
            throw new IllegalArgumentException(getClass() + " can only parse objects of type String,"
                    + " char[], File, InputStream, or URL " + " [attempt to parse object of " + nhx_source.getClass()
                    + "].");
        }
    }

    private final Phylogeny finishPhylogeny()
            throws PhylogenyParserException, NHXFormatException, PhyloXmlDataFormatException {
        if (_current_phylogeny != null) {
            parseNHX(_current_anotation != null ? _current_anotation.toString() : "",
                    _current_phylogeny.getRoot(),
                    getTaxonomyExtraction(),
                    isReplaceUnderscores(),
                    isAllowErrorsInDistanceToParent(),
                    true,
                    isParseBeastStyleExtendedTags());
            if (GUESS_IF_SUPPORT_VALUES) {
                if (isBranchLengthsLikeBootstrapValues(_current_phylogeny)) {
                    moveBranchLengthsToConfidenceValues(_current_phylogeny);
                }
            }
            if (isGuessRootedness()) {
                final PhylogenyNode root = _current_phylogeny.getRoot();
                if ((root.getDistanceToParent() >= 0.0) || !ForesterUtil.isEmpty(root.getName())
                        || !ForesterUtil.isEmpty(PhylogenyMethods.getSpecies(root)) || root.isHasAssignedEvent()) {
                    _current_phylogeny.setRooted(true);
                }
            }
            if (isParseBeastStyleExtendedTags() && _normalize_bracket_annotations) {
                // what no single annotation can say: whose they are, and whether the tree confirms its dates
                BracketAnnotationNormalizer.normalize(_current_phylogeny);
            }
            return _current_phylogeny;
        }
        return null;
    }

    private final Phylogeny finishSingleNodePhylogeny()
            throws PhylogenyParserException, NHXFormatException, PhyloXmlDataFormatException {
        final PhylogenyNode new_node = new PhylogenyNode();
        parseNHX(_current_anotation.toString(),
                new_node,
                getTaxonomyExtraction(),
                isReplaceUnderscores(),
                isAllowErrorsInDistanceToParent(),
                true,
                isParseBeastStyleExtendedTags());
        _current_phylogeny = new Phylogeny();
        _current_phylogeny.setRoot(new_node);
        return _current_phylogeny;
    }

    private final void init() {
        setTaxonomyExtraction(TAXONOMY_EXTRACTION.NO);
        setReplaceUnderscores(REPLACE_UNDERSCORES_DEFAULT);
        setGuessRootedness(GUESS_ROOTEDNESS_DEFAULT);
        setIgnoreQuotes(IGNORE_QUOTES_DEFAULT);
        setAllowErrorsInDistanceToParent(ALLOW_ERRORS_IN_DISTANCE_TO_PARENT_DEFAULT);
        setParseBeastStyleExtendedTags(PARSE_BRACKET_ANNOTATIONS_DEFAULT);
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
        if (_source == null) {
            throw new IOException("source is not set");
        }
        while (true) {
            char c = '\b';
            if (_input_type == BUFFERED_READER) {
                final int ci = _my_source_br.read();
                if (ci >= 0) {
                    c = (char) ci;
                } else {
                    break;
                }
            } else {
                if (_i >= _source_length) {
                    break;
                }
                c = _my_source_charary[_i];
            }
            // '' (or "") is the escape for a literal quote INSIDE a quoted label -- 'Seba''s bat' is
            // one label reading "Seba's bat", not a label ending and another starting. Only a pair on
            // ADJACENT characters counts, so this is consumed and cleared for every character.
            final char closed_quote = _just_closed_quote;
            _just_closed_quote = 0;
            if (_blob_quote_remaining > 0) {
                // inside a quoted VALUE of a kept blob, up to its closing quote (found when it opened): all data
                --_blob_quote_remaining;
                if (!((c < 32) || (c == 127))) {
                    _current_anotation.append((_blob_quote_remaining == 0) ? c : blobQuotedChar(c));
                }
                ++_i;
                continue;
            }
            // The character right after a kept annotation's opening '&' decides whether this is the legacy
            // "&&NHX" tag scheme (see _in_legacy_nhx_tag) -- but white space there decides nothing: a control
            // character is always dropped (the "\n\t is always ignored" rule below), and "[ & & NHX : S = x ]" has
            // always been read as the NHX tag it is.
            if (_awaiting_second_annotation_char && !((c < 33) || (c == 127))) {
                _awaiting_second_annotation_char = false;
                if (c == '&') {
                    _in_legacy_nhx_tag = true;
                }
            }
            // Data, not the "&&NHX" formatting-noise scheme: whitespace/quotes inside it are real and must survive.
            final boolean raw_annotation_data = _in_kept_annotation && !_in_legacy_nhx_tag;
            if (!_in_single_quote && !_in_double_quote) {
                if (c == ':') {
                    _saw_colon = true;
                } else if (!((c < 33) || (c == 127)) && _saw_colon
                        && ((c != '[') && (c != '.') && ((c < 48) || (c > 57)))) {
                    _saw_colon = false;
                }
                if (_in_open_bracket && (c == ']')) {
                    _in_open_bracket = false;
                    _in_kept_annotation = false;
                    _in_legacy_nhx_tag = false;
                }
            }
            // \n\t is always ignored,
            // "=34  '=39 space=32
            // A plain space is dropped everywhere EXCEPT inside a kept, non-legacy "[&...]" blob, where it is
            // DATA: real Auspice/Nextstrain Nexus writes unquoted spaces in values ("country=Democratic Republic
            // of the Congo", "outbreak_geo=Kikwit 1995"), and squashing them silently produced a different value
            // here than the JSON path gives for the same dataset. Gated on raw_annotation_data, NOT
            // _in_open_bracket: a plain "[comment]", the leading character of a "[91]" numeric confidence, and a
            // legacy "[&&NHX:...]" tag must all still have their spaces swallowed exactly as before (a comment's
            // content is discarded either way, but the single character right after '[' decides
            // comment-vs-annotation-vs-numeric and must not be a stray space; NHX tag values are underscore-based
            // identifiers where whitespace has always been pure formatting noise -- testNHXParsingQuotes pins
            // "mo\tnkey !" reducing to "monkey!").
            if ((c < 32) || (c == 127) || (isIgnoreQuotes() && ((c == 32) || (c == 34) || (c == 39)))
                    || ((c == 32) && (!_in_single_quote && !_in_double_quote) && !raw_annotation_data)
                    || ((_clade_level == 0) && (c == ';') && (!_in_single_quote && !_in_double_quote))) {
                //do nothing
            } else if (_in_comment) {
                if (c == ']') {
                    _in_comment = false;
                }
            } else if (_in_double_quote) {
                if (c == '"') {
                    _in_double_quote = false;
                    _just_closed_quote = '"';
                } else {
                    _current_anotation.append(changeCharInParens(c));
                }
            } else if ((c == '"') && !_in_single_quote) {
                if (raw_annotation_data) {
                    openBlobQuote(c);
                } else {
                    if (closed_quote == '"') {
                        _current_anotation.append('"');
                    }
                    _in_double_quote = true;
                }
            } else if (_in_single_quote) {
                if (c == 39) {
                    _in_single_quote = false;
                    _just_closed_quote = 39;
                } else {
                    _current_anotation.append(changeCharInParens(c));
                }
            } else if (c == 39) {
                if (raw_annotation_data) {
                    openBlobQuote(c);
                } else {
                    if (closed_quote == 39) {
                        _current_anotation.append('\'');
                    }
                    _in_single_quote = true;
                }
            } else if (c == '[') {
                _saw_open_bracket = true;
                _in_open_bracket = true;
            } else if (_saw_open_bracket) {
                if (c != ']') {
                    // everything not starting with "[&" is considered a comment
                    // unless ":digits and/or . [bootstrap]":
                    if (c == '&') {
                        _current_anotation.append("[&");
                        _in_kept_annotation = true;
                        _awaiting_second_annotation_char = true;
                    } else if ((_saw_colon || _after_close_paren)
                            && (((c > 47) && (c < 58)) || (c == 46) || (c == 45) || (c == 43))) {
                        _current_anotation.append("[" + c);
                    } else {
                        _in_comment = true;
                    }
                }
                // comment consisting just of "[]":
                _saw_open_bracket = false;
            } else if ((c == '(') && !_in_open_bracket) {
                _after_close_paren = false;
                final Phylogeny phy = processOpenParen();
                if (phy != null) {
                    ++_i;
                    _next = phy;
                    return;
                }
            } else if ((c == ')') && !_in_open_bracket) {
                _after_close_paren = true;
                processCloseParen();
            } else if ((c == ',') && !_in_open_bracket) {
                _after_close_paren = false;
                processComma();
            } else {
                _current_anotation.append(c);
            }
            ++_i;
        } //  while ( true )
        if (_clade_level != 0) {
            throw new PhylogenyParserException("error in NH (Newick) formatted data: most likely cause: number of open parens does not equal number of close parens");
        }
        if (_current_phylogeny != null) {
            _next = finishPhylogeny();
            _current_phylogeny = null;
            _current_anotation = null;
        } else if ((_current_anotation != null) && (_current_anotation.length() > 0)) {
            _next = finishSingleNodePhylogeny();
            _current_anotation = null;
        } else {
            _next = null;
        }
    }

    /**
     * A quote character inside a kept "[&...]" blob. It is DATA -- "country=Côte d'Ivoire", a real Nextstrain export
     * -- unless the same character closes it standing where a value can END (before ',', '}' or ']'); then the run
     * between them is one quoted stretch whose ']' does not end the blob. The search for that partner gives up at a
     * ']' followed by Newick structure, so a bare value that merely BEGINS with an apostrophe
     * ("division='s-Hertogenbosch") cannot reach into the next node's blob for one. Without this an odd count of
     * apostrophes fails the parse, and an even count is worse: two of them pair up across tips and silently swallow
     * every bracket, comma and paren -- whole tips -- in between. The accepted cost of the bound: a QUOTED value
     * containing "]," ends its blob early ("a]b" is fine).
     * <p>
     * Whether such a run is a quoted VALUE (it must also START where a value can: after '=', ',', '{') is decided
     * where it matters, in {@link BeastAnnotationParser#splitTopLevel}. Here it would change nothing: inside a blob
     * the only structural character is ']', a well-formed blob's own ']' is always followed by Newick structure, so
     * the bound already keeps every run inside its blob -- and what a run does to the text is undone downstream.
     * JOINT with Archaeopteryx.js (opensBlobQuote / blobQuoteClose in forester.js).
     */
    private final void openBlobQuote(final char q) throws IOException {
        final int close = blobQuoteCloseOffset(q);
        _current_anotation.append(q);
        _blob_quote_remaining = Math.max(close, 0);
    }

    /** How many characters ahead the quote closing a run opened HERE stands, or -1 when nothing closes it. Looks
     *  ahead without consuming (the reader is marked and reset). */
    private final int blobQuoteCloseOffset(final char q) throws IOException {
        final boolean from_reader = _input_type == BUFFERED_READER;
        if (from_reader) {
            _my_source_br.mark(BLOB_QUOTE_LOOKAHEAD);
        }
        try {
            int pending = 0; // what waits for the next non-space character: 1 a candidate closing quote, 2 a ']'
            int pending_offset = -1;
            for (int offset = 1; offset < BLOB_QUOTE_LOOKAHEAD; ++offset) {
                final int ci = from_reader ? _my_source_br.read()
                        : (((_i + offset) < _source_length) ? _my_source_charary[_i + offset] : -1);
                if (ci < 0) {
                    return -1; // the input ends inside the blob: nothing closes here
                }
                final char c = (char) ci;
                if (pending != 0) {
                    if (Character.isWhitespace(c)) {
                        continue;
                    }
                    if (pending == 1) {
                        if ((c == ',') || (c == '}') || (c == ']')) {
                            return pending_offset;
                        }
                    } else if ((c == ',') || (c == ')') || (c == ':') || (c == ';') || (c == '(') || (c == '[')) {
                        return -1; // the blob has ended: this quote never was an opening one
                    }
                    pending = 0;
                }
                if (c == q) {
                    pending = 1;
                    pending_offset = offset;
                } else if (c == ']') {
                    pending = 2;
                }
            }
            return -1;
        } finally {
            if (from_reader) {
                _my_source_br.reset();
            }
        }
    }

    /** Inside a blob's quoted value ':' '[' ']' are data, but the node-level parse still splits on them: they
     *  travel as placeholders, which {@link BeastAnnotationParser} turns back. */
    private final static char blobQuotedChar(final char c) {
        if (c == ':') {
            return BELL;
        } else if (c == '[') {
            return BLOB_OPEN_BRACKET;
        } else if (c == ']') {
            return BLOB_CLOSE_BRACKET;
        }
        return c;
    }

    private final static char changeCharInParens(char c) {
        if (c == ':') {
            c = BELL;
        } else if (c == '[') {
            c = '{';
        } else if (c == ']') {
            c = '}';
        }
        return c;
    }

    private final void processCloseParen()
            throws PhylogenyParserException, NHXFormatException, PhyloXmlDataFormatException {
        if (_clade_level < 0) {
            throw new PhylogenyParserException("error in NH (Newick)/NHX formatted data: most likely cause: number of close parens is larger than number of open parens");
        }
        --_clade_level;
        if (!_saw_closing_paren) {
            final PhylogenyNode new_node = new PhylogenyNode();
            parseNHX(_current_anotation.toString(),
                    new_node,
                    getTaxonomyExtraction(),
                    isReplaceUnderscores(),
                    isAllowErrorsInDistanceToParent(),
                    true,
                    isParseBeastStyleExtendedTags());
            _current_anotation = new StringBuilder();
            _current_node.addAsChild(new_node);
        } else {
            parseNHX(_current_anotation.toString(),
                    _current_node.getLastChildNode(),
                    getTaxonomyExtraction(),
                    isReplaceUnderscores(),
                    isAllowErrorsInDistanceToParent(),
                    true,
                    isParseBeastStyleExtendedTags());
            _current_anotation = new StringBuilder();
        }
        if (!_current_node.isRoot()) {
            _current_node = _current_node.getParent();
        }
        _saw_closing_paren = true;
    }

    private final void processComma() throws PhylogenyParserException, NHXFormatException, PhyloXmlDataFormatException {
        if (!_saw_closing_paren) {
            final PhylogenyNode new_node = new PhylogenyNode();
            parseNHX(_current_anotation.toString(),
                    new_node,
                    getTaxonomyExtraction(),
                    isReplaceUnderscores(),
                    isAllowErrorsInDistanceToParent(),
                    true,
                    isParseBeastStyleExtendedTags());
            if (_current_node == null) {
                throw new NHXFormatException("format might not be NH or NHX");
            }
            _current_node.addAsChild(new_node);
        } else {
            parseNHX(_current_anotation.toString(),
                    _current_node.getLastChildNode(),
                    getTaxonomyExtraction(),
                    isReplaceUnderscores(),
                    isAllowErrorsInDistanceToParent(),
                    true,
                    isParseBeastStyleExtendedTags());
        }
        _current_anotation = new StringBuilder();
        _saw_closing_paren = false;
    }

    private final Phylogeny processOpenParen()
            throws PhylogenyParserException, NHXFormatException, PhyloXmlDataFormatException {
        Phylogeny phy = null;
        final PhylogenyNode new_node = new PhylogenyNode();
        if (_clade_level == 0) {
            if (_current_phylogeny != null) {
                phy = finishPhylogeny();
            }
            _clade_level = 1;
            _current_anotation = new StringBuilder();
            _current_phylogeny = new Phylogeny();
            _current_phylogeny.setRoot(new_node);
        } else {
            ++_clade_level;
            _current_node.addAsChild(new_node);
        }
        _current_node = new_node;
        _saw_closing_paren = false;
        return phy;
    }

    private final static NHXParser createInstance(final Object nhx_source) throws NHXFormatException, IOException {
        final NHXParser parser = new NHXParser();
        parser.setSource(nhx_source);
        return parser;
    }

    public final static Phylogeny[] parse(final Object nhx_source) throws NHXFormatException, IOException {
        return NHXParser.createInstance(nhx_source).parse();
    }

    public final static void parseNHX(String s,
                                      final PhylogenyNode node_to_annotate,
                                      final TAXONOMY_EXTRACTION taxonomy_extraction,
                                      final boolean replace_underscores,
                                      final boolean allow_errors_in_distance_to_parent,
                                      final boolean replace_bell,
                                      final boolean parse_beast_style_extended_tags)
            throws NHXFormatException, PhyloXmlDataFormatException {
        if ((taxonomy_extraction != TAXONOMY_EXTRACTION.NO) && replace_underscores) {
            throw new IllegalArgumentException("cannot extract taxonomies and replace under scores at the same time");
        }
        if ((s != null) && (s.length() > 0)) {
            boolean is_nhx = false;
            if (parse_beast_style_extended_tags && (s.indexOf('[') > -1)) {
                // EVERY bracket group of the node, not only the first, and the text around them: MrBayes writes
                // "[&prob=...]:0.04129[&length_mean=...]" -- a length BETWEEN two groups -- and FigTree leads a
                // BEAST blob with "!color". The keys belong to the structured parser whatever leads the blob, and
                // underscores are replaced in the label and in legacy NHX tag values (S=Homo_sapiens), as they always
                // were -- never inside a key=value blob ("height_median" is a key, not a name).
                final StringBuilder outside = new StringBuilder();
                final List<String> groups = bracketGroups(s, outside);
                final StringBuilder blob = new StringBuilder();
                final StringBuilder b = new StringBuilder();
                for (final String group : groups) {
                    if (isNhxTagSyntax(group)) {
                        is_nhx = true;
                        b.append(group.substring(group.indexOf(':')));
                    } else if (group.startsWith("&")) {
                        if (blob.length() > 0) {
                            blob.append(',');
                        }
                        blob.append(group.substring(1));
                    } else if (NUMBERS_ONLY_PATTERN.matcher(group.trim()).matches()) {
                        // No &&NHX and digits only: is likely to be a support value.
                        b.append(":" + NHXtags.SUPPORT + group.trim());
                    }
                }
                if (blob.length() > 0) {
                    BeastAnnotationParser.apply(blob.toString(), node_to_annotate);
                }
                s = outside.toString() + b; // the label and the legacy NHX tags -- never the key=value blobs
                if (replace_underscores) {
                    s = s.replaceAll("_+", " ");
                }
                s = s.replaceAll("\\s+", " ").trim();
            } else {
                if (replace_underscores) {
                    s = s.replaceAll("_+", " ");
                }
                s = s.replaceAll("\\s+", " ").trim();
            }
            final int ob = s.indexOf("[");
            if (ob > -1) {
                String b = "";
                // Only a genuine "[&&NHX" tag states the node's taxonomy itself; a BEAST blob, a "[91]" support or a
                // kept comment says nothing about it, and must not switch the name-based extraction off.
                is_nhx = s.indexOf("&&NHX") == (ob + 1);
                final int cb = s.indexOf("]");
                if (cb < 0) {
                    throw new NHXFormatException("error in NHX formatted data: no closing \"]\" in \"" + s + "\"");
                }
                // Only the option-OFF path gets here with a bracket (the option-ON path above has consumed them
                // all): ONE group per node, kept bit for bit as it always was.
                if (s.indexOf("&&NHX") == (ob + 1)) {
                    b = s.substring(ob + 6, cb);
                } else if (s.indexOf("&") == (ob + 1) && s.indexOf("[&prob") == -1 &&
                        s.indexOf("[&boot") == -1 &&
                        s.indexOf("[&!colo") == -1) {
                    final String bracketed = s.substring(ob + 1, cb);
                    b = ":" + NHXtags.COMMENT + bracketed; // keep the raw comment when the option is off
                    final Matcher ewn_matcher = ENDS_WITH_NUMBER_PATTERN.matcher(s);
                    if (ewn_matcher.find()) {
                        b = b + ewn_matcher.group(1);
                    }
                } else {
                    // No &&NHX and digits only: is likely to be a support value.
                    final String bracketed = s.substring(ob + 1, cb);
                    final Matcher numbers_only = NUMBERS_ONLY_PATTERN.matcher(bracketed);
                    if (numbers_only.matches()) {
                        b = ":" + NHXtags.SUPPORT + bracketed;
                    } else if (s.indexOf("prob=") > -1) {
                        processMrBayes3Data(s, node_to_annotate);
                    }
                    final Matcher ewn_matcher = ENDS_WITH_NUMBER_PATTERN.matcher(s);
                    if (ewn_matcher.find()) {
                        b = b + ewn_matcher.group(1);
                    }
                }
                s = s.substring(0, ob) + b;
                if ((s.indexOf("[") > -1) || (s.indexOf("]") > -1)) {
                    throw new NHXFormatException("error in NHX formatted data: more than one \"]\" or \"[\"");
                }
            }

            final StringTokenizer t = new StringTokenizer(s, ":");
            if (t.countTokens() > 0) {
                if (!s.startsWith(":")) {
                    if ((s.indexOf(BELL) <= -1) || !replace_bell) {
                        node_to_annotate.setName(t.nextToken());
                    } else {
                        node_to_annotate.setName(t.nextToken().replace(BELL, ':'));
                    }
                    if (!replace_underscores && (!is_nhx && (taxonomy_extraction != TAXONOMY_EXTRACTION.NO))) {
                        ParserUtils.extractTaxonomyDataFromNodeName(node_to_annotate, taxonomy_extraction);
                    }
                }
                while (t.hasMoreTokens()) {
                    s = trimNhxField(t.nextToken());
                    if ((s.indexOf(BELL) > -1) && replace_bell) {
                        s = s.replace(BELL, ':');
                    }
                    if (BRANCH_LENGTH_START_PATTERN.matcher(s).find()) {
                        if ((node_to_annotate.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT)
                                && !allow_errors_in_distance_to_parent) {
                            throw new NHXFormatException("error in NHX formatted data: more than one distance to parent:"
                                    + "\"" + s + "\"");
                        }
                        node_to_annotate.setDistanceToParent(doubleValue(s, allow_errors_in_distance_to_parent));
                    } else if (s.startsWith(NHXtags.SPECIES_NAME)) {
                        if (!node_to_annotate.getNodeData().isHasTaxonomy()) {
                            node_to_annotate.getNodeData().setTaxonomy(new Taxonomy());
                        }
                        node_to_annotate.getNodeData().getTaxonomy().setScientificName(s.substring(2));
                    } else if (s.startsWith(NHXtags.IS_DUPLICATION)) {
                        if ((s.charAt(2) == 'Y') || (s.charAt(2) == 'T')) {
                            node_to_annotate.getNodeData().setEvent(Event.createSingleDuplicationEvent());
                        } else if ((s.charAt(2) == 'N') || (s.charAt(2) == 'F')) {
                            node_to_annotate.getNodeData().setEvent(Event.createSingleSpeciationEvent());
                        } else if (s.charAt(2) == '?') {
                            node_to_annotate.getNodeData().setEvent(Event.createSingleSpeciationOrDuplicationEvent());
                        } else {
                            throw new NHXFormatException("error in NHX formatted data: :D=Y or :D=N or :D=?");
                        }
                    } else if (s.startsWith(NHXtags.SUPPORT)) {
                        PhylogenyMethods.setConfidence(node_to_annotate, doubleValue(s.substring(2), false));
                    } else if (s.startsWith(NHXtags.TAXONOMY_ID)) {
                        if (!node_to_annotate.getNodeData().isHasTaxonomy()) {
                            node_to_annotate.getNodeData().setTaxonomy(new Taxonomy());
                        }
                        node_to_annotate.getNodeData().getTaxonomy()
                                .setIdentifier(new Identifier(s.substring(2)));
                    } else if (s.startsWith(NHXtags.SEQUENCE_ACCESSION)) {
                        if (!node_to_annotate.getNodeData().isHasSequence()) {
                            node_to_annotate.getNodeData().setSequence(new Sequence());
                        }
                        node_to_annotate.getNodeData().getSequence()
                                .setAccession(new Accession(s.substring(3), "?"));
                    } else if (s.startsWith(NHXtags.GENE_NAME)) {
                        if (!node_to_annotate.getNodeData().isHasSequence()) {
                            node_to_annotate.getNodeData().setSequence(new Sequence());
                        }
                        node_to_annotate.getNodeData().getSequence().setName(s.substring(3));
                    } else if (s.startsWith(NHXtags.COMMENT)) {
                        // "C=" + the comment; a bracket blob kept as a comment (option off) arrives as "C=&..."
                        String comment = s.substring(NHXtags.COMMENT.length());
                        if (comment.startsWith("&")) {
                            comment = comment.substring(1);
                        }
                        comment = comment.trim().replace(BLOB_OPEN_BRACKET, '[').replace(BLOB_CLOSE_BRACKET, ']');
                        if (!ForesterUtil.isEmpty(comment)) {
                            // MERGED into what the node already carries: a bracket group parsed before this tag
                            // (the option-ON path reads them all) has put its properties there
                            PropertiesList custom_data = node_to_annotate.getNodeData().getProperties();
                            if (custom_data == null) {
                                custom_data = new PropertiesList();
                                node_to_annotate.getNodeData().setProperties(custom_data);
                            }
                            custom_data.addProperty(new Property(ForesterConstants.NH_COMMENT, comment, "", "xsd:string", Property.AppliesTo.NODE));
                        }
                    }
                } // while ( t.hasMoreTokens() )
            }
            if (parse_beast_style_extended_tags) {
                useLengthMedianWhenNoLength(node_to_annotate);
            }
        }
    }

    /** The file's literal ":length" is the branch length (Christian, 2026-09-16). Only when a node states NONE does
     *  MrBayes' / TreeAnnotator's length_median stand in -- what the option-OFF path has always done for a MrBayes
     *  node without a literal length. */
    private final static void useLengthMedianWhenNoLength(final PhylogenyNode node) {
        if ((node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT)
                || (node.getNodeData().getProperties() == null)) {
            return;
        }
        final List<Property> median = node.getNodeData().getProperties().getProperties("beast:length_median");
        if (median.isEmpty()) {
            return;
        }
        final Double d = BeastAnnotationParser.parseNumber(median.get(0).getValue());
        if ((d != null) && (d.doubleValue() >= 0.0)) {
            node.setDistanceToParent(d.doubleValue());
        }
    }

    private final static double doubleValue(final String str, final boolean allow_errors) throws NHXFormatException {
        try {
            return Double.parseDouble(str);
        } catch (final NumberFormatException ex) {
            if (!allow_errors) {
                throw new NHXFormatException("error in NH/NHX formatted data: failed to parse number from " + "\""
                        + str + "\"");
            }
        }
        return 0.0;
    }

    private final static boolean isBranchLengthsLikeBootstrapValues(final Phylogeny p) {
        final PhylogenyNodeIterator it = p.iteratorExternalForward();
        final double d0 = it.next().getDistanceToParent();
        if ((d0 < 10) || !it.hasNext()) {
            return false;
        }
        while (it.hasNext()) {
            final double d = it.next().getDistanceToParent();
            if ((d != d0) || (d < 10)) {
                return false;
            }
        }
        return true;
    }

    private final static void moveBranchLengthsToConfidenceValues(final Phylogeny p) {
        final PhylogenyNodeIterator it = p.iteratorPostorder();
        while (it.hasNext()) {
            final PhylogenyNode n = it.next();
            PhylogenyMethods.setBootstrapConfidence(n, n.getDistanceToParent());
            n.setDistanceToParent(PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT);
        }
    }

    /** One "tag=value" field of an NHX tag with the padding around its tag and around its value removed: a quoted
     *  run keeps the white space INSIDE it ("homo sapiens"), but padding at its ends (S=" homo ") is no more part
     *  of the value than the quotes are. JOINT with Archaeopteryx.js, which trims each field. */
    final static String trimNhxField(final String field) {
        final int eq = field.indexOf('=');
        if (eq < 0) {
            return field.trim();
        }
        return field.substring(0, eq).trim() + "=" + field.substring(eq + 1).trim();
    }

    /** NHX tags -- ":S=species:B=91" -- rather than a "key=value,key=value" blob: an '&'-led group with a ':' before
     *  its first '=' (a BEAST / FigTree / MrBayes key never contains a colon). That is the genuine "&&NHX:..." and
     *  also the sloppy spellings this parser has always forgiven -- "&NHX:S=x", "&:S=x", "&&NH:S=x", "&&:S=x". */
    final static boolean isNhxTagSyntax(final String group) {
        if (!group.startsWith("&")) {
            return false;
        }
        final int colon = group.indexOf(':');
        final int eq = group.indexOf('=');
        return (colon > -1) && ((eq < 0) || (colon < eq));
    }

    /** The contents of every top-level {@code [...]} group of one node's annotation text, in order (without their
     *  brackets), with everything OUTSIDE the groups appended to {@code outside} -- so a branch length sitting
     *  between two groups ({@code "A[&prob=1]:0.04[&length_mean=0.05]"}) survives. The streaming scanner has already
     *  turned a bracket inside a quoted value into a brace, so every bracket here is a real one. */
    final static List<String> bracketGroups(final String s, final StringBuilder outside) throws NHXFormatException {
        final List<String> groups = new ArrayList<String>();
        int pos = 0;
        while (pos < s.length()) {
            final int ob = s.indexOf('[', pos);
            if (ob < 0) {
                break;
            }
            final int cb = s.indexOf(']', ob);
            if (cb < 0) {
                throw new NHXFormatException("error in NHX formatted data: no closing \"]\" in \"" + s + "\"");
            }
            outside.append(s, pos, ob);
            groups.add(s.substring(ob + 1, cb));
            pos = cb + 1;
        }
        outside.append(s, pos, s.length());
        if (outside.indexOf("]") > -1) {
            throw new NHXFormatException("error in NHX formatted data: a \"]\" without its \"[\" in \"" + s + "\"");
        }
        return groups;
    }

    private final static void processMrBayes3Data(final String s, final PhylogenyNode node_to_annotate)
            throws NHXFormatException {
        double sd = -1;
        final Matcher mb_prob_sd_matcher = MB_PROB_SD_PATTERN.matcher(s);
        if (mb_prob_sd_matcher.find()) {
            try {
                sd = Double.parseDouble(mb_prob_sd_matcher.group(1));
            } catch (final NumberFormatException e) {
                throw new NHXFormatException("failed to parse probability standard deviation (Mr Bayes output) from \""
                        + s + "\"");
            }
        }
        final Matcher mb_prob_matcher = MB_PROB_PATTERN.matcher(s);
        if (mb_prob_matcher.find()) {
            double prob = -1;
            try {
                prob = Double.parseDouble(mb_prob_matcher.group(1));
            } catch (final NumberFormatException e) {
                throw new NHXFormatException("failed to parse probability (Mr Bayes output) from \"" + s + "\"");
            }
            if (prob >= 0.0) {
                if (sd >= 0.0) {
                    node_to_annotate.getBranchData()
                            .addConfidence(new Confidence(prob, "posterior probability", sd));
                } else {
                    node_to_annotate.getBranchData().addConfidence(new Confidence(prob, "posterior probability"));
                }
            }
        }
        final Matcher mb_bl_matcher = MB_BL_PATTERN.matcher(s);
        if (mb_bl_matcher.find()) {
            double bl = -1;
            try {
                bl = Double.parseDouble(mb_bl_matcher.group(1));
            } catch (final NumberFormatException e) {
                throw new NHXFormatException("failed to parse median branch length (Mr Bayes output) from \"" + s
                        + "\"");
            }
            if (bl >= 0.0) {
                node_to_annotate.setDistanceToParent(bl);
            }
        }
    }

    public String getEncoding() {
        return _encoding;
    }

    private final boolean isParseBeastStyleExtendedTags() {
        return _parse_beast_style_extended_tags;
    }

    /** Whether the per-tree pass ({@link BracketAnnotationNormalizer}) runs when a tree is finished (default true).
     *  Off ONLY for a caller that adds annotations of its own afterwards and runs the pass itself -- the Nexus
     *  reader, whose TAXLABELS annotations reach the tips after the tree string: normalizing before them would leave
     *  one tree with both beast: and treetime: refs. */
    public final void setNormalizeBracketAnnotations(final boolean normalize_bracket_annotations) {
        _normalize_bracket_annotations = normalize_bracket_annotations;
    }

    public final void setParseBeastStyleExtendedTags(final boolean parse_beast_style_extended_tags) {
        _parse_beast_style_extended_tags = parse_beast_style_extended_tags;
    }

    public static enum TAXONOMY_EXTRACTION {
        AGGRESSIVE,
        NO,
        PFAM_STYLE_RELAXED,
        PFAM_STYLE_STRICT;
    }
}
