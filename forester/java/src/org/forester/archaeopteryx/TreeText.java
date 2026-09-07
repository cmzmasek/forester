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

package org.forester.archaeopteryx;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.util.ForesterUtil;

/**
 * The text side of the "Tree as Text" window, free of Swing: which {@link Format}s exist, how a tree is rendered
 * to each, and a light tokenizer ({@link #spans}) that says which stretches of the text are markup rather than
 * data -- so the window can mute tags, punctuation and numbers and let the names stand out.
 */
final class TreeText {

    /** The three text formats a tree can be shown in. */
    enum Format {
        PHYLOXML( "phyloXML", ".xml" ), NEWICK( "Newick", ".nwk" ), NEXUS( "Nexus", ".nex" );

        final String label;
        final String suffix;

        Format( final String label, final String suffix ) {
            this.label = label;
            this.suffix = suffix;
        }

        /** Whether this format is one long line by nature (so wrapping is on by default). */
        boolean wrapsByDefault() {
            return this != PHYLOXML;
        }
    }

    /** What a span of the text is, for tinting. Everything not covered by a span is data (plain). */
    enum Kind {
        /** An XML tag (its brackets and name), a Newick/Nexus bracket, comma, colon or semicolon. */
        MARKUP,
        /** An XML attribute name (with its "="). */
        ATTRIBUTE,
        /** A number that is markup-ish data: a branch length or a support value. */
        NUMBER,
        /** A Nexus keyword (#NEXUS, BEGIN, END, TREE, ...). */
        KEYWORD,
        /** An XML comment / declaration, or a Newick "[...]" comment (incl. NHX blocks). */
        COMMENT
    }

    /** A half-open [start, end) stretch of the text of one {@link Kind}. */
    static final class Span {

        final int  start;
        final int  end;
        final Kind kind;

        Span( final int start, final int end, final Kind kind ) {
            this.start = start;
            this.end = end;
            this.kind = kind;
        }

        @Override
        public String toString() {
            return kind + "[" + start + "," + end + ")";
        }
    }

    /** Texts longer than this are shown plain (tinting a multi-megabyte document is not worth the wait). */
    static final int TINT_LIMIT = 512 * 1024;

    /** The words that can open a Nexus COMMAND (a statement terminated by ";"). A word is only matched against
     *  this in command position -- at the start of a statement -- so a taxon that happens to be called "Matrix"
     *  stays plain data. */
    private static final Set<String> NEXUS_COMMANDS = new HashSet<>( Arrays.asList( "#nexus", "begin", "end",
            "endblock", "tree", "utree", "translate", "taxlabels", "dimensions", "format", "matrix", "options",
            "charlabels", "charstatelabels", "statelabels", "link", "title" ) );
    /** Commands whose body is a list of {@code Name=Value} settings, so the names in it are tinted like XML
     *  attributes ({@code Dimensions NTax=9;}, {@code Format DataType=DNA;}). */
    private static final Set<String> NEXUS_SETTING_COMMANDS = new HashSet<>( Arrays.asList( "dimensions", "format",
            "options" ) );

    private TreeText() {
    }

    /** The tree in {@code format} (what File &gt; Save As writes, modulo options); "" for no tree. */
    static String render( final Phylogeny phy, final Format format, final NH_CONVERSION_SUPPORT_VALUE_STYLE style ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return "";
        }
        switch ( format ) {
            case PHYLOXML:
                return phy.toPhyloXML( 0 );
            case NEWICK:
                return phy.toNewHampshire( style );
            case NEXUS:
                return phy.toNexus( style );
            default:
                return "";
        }
    }

    /** "12,345 characters · 210 lines" (a one-line Newick reports "1 line"). */
    static String sizeText( final String text ) {
        final int chars = text.length();
        int lines = text.isEmpty() ? 0 : 1;
        for( int i = 0; i < chars; ++i ) {
            if ( text.charAt( i ) == '\n' ) {
                ++lines;
            }
        }
        if ( text.endsWith( "\n" ) ) {
            --lines;
        }
        return TreeFacts.count( chars ) + ( chars == 1 ? " character" : " characters" ) + " · "
                + TreeFacts.count( lines ) + ( lines == 1 ? " line" : " lines" );
    }

    /** The tint spans of {@code text} in {@code format}, in order, non-overlapping. Empty above {@link #TINT_LIMIT}. */
    static List<Span> spans( final String text, final Format format ) {
        if ( ( text == null ) || ( text.length() > TINT_LIMIT ) ) {
            return new ArrayList<>();
        }
        return ( format == Format.PHYLOXML ) ? xmlSpans( text ) : newickSpans( text, format == Format.NEXUS );
    }

    // ------------------------------------------------------------------ XML
    private static List<Span> xmlSpans( final String t ) {
        final List<Span> out = new ArrayList<>();
        final int n = t.length();
        int i = 0;
        while( i < n ) {
            final int lt = t.indexOf( '<', i );
            if ( lt < 0 ) {
                break;
            }
            if ( t.startsWith( "<!--", lt ) ) {
                final int e = t.indexOf( "-->", lt + 4 );
                final int end = ( e < 0 ) ? n : e + 3;
                out.add( new Span( lt, end, Kind.COMMENT ) );
                i = end;
                continue;
            }
            if ( t.startsWith( "<?", lt ) || t.startsWith( "<!", lt ) ) {
                final int e = t.indexOf( '>', lt + 2 );
                final int end = ( e < 0 ) ? n : e + 1;
                out.add( new Span( lt, end, Kind.COMMENT ) );
                i = end;
                continue;
            }
            // a tag: "<" or "</" + name, then attributes, then ">" or "/>"
            int j = lt + 1;
            if ( ( j < n ) && ( t.charAt( j ) == '/' ) ) {
                ++j;
            }
            while( ( j < n ) && isNameChar( t.charAt( j ) ) ) {
                ++j;
            }
            out.add( new Span( lt, j, Kind.MARKUP ) );
            // attributes up to the closing bracket (quotes may contain '>')
            while( j < n ) {
                final char c = t.charAt( j );
                if ( c == '>' ) {
                    out.add( new Span( j, j + 1, Kind.MARKUP ) );
                    ++j;
                    break;
                }
                if ( ( c == '/' ) && ( j + 1 < n ) && ( t.charAt( j + 1 ) == '>' ) ) {
                    out.add( new Span( j, j + 2, Kind.MARKUP ) );
                    j += 2;
                    break;
                }
                if ( Character.isWhitespace( c ) ) {
                    ++j;
                    continue;
                }
                if ( isNameChar( c ) ) {
                    final int s = j;
                    while( ( j < n ) && isNameChar( t.charAt( j ) ) ) {
                        ++j;
                    }
                    while( ( j < n ) && Character.isWhitespace( t.charAt( j ) ) ) {
                        ++j;
                    }
                    if ( ( j < n ) && ( t.charAt( j ) == '=' ) ) {
                        ++j;
                        out.add( new Span( s, j, Kind.ATTRIBUTE ) );
                        while( ( j < n ) && Character.isWhitespace( t.charAt( j ) ) ) {
                            ++j;
                        }
                        if ( ( j < n ) && ( ( t.charAt( j ) == '"' ) || ( t.charAt( j ) == '\'' ) ) ) {
                            final char q = t.charAt( j );
                            final int e = t.indexOf( q, j + 1 );
                            j = ( e < 0 ) ? n : e + 1; // the value stays plain
                        }
                    }
                    else {
                        out.add( new Span( s, j, Kind.ATTRIBUTE ) );
                    }
                    continue;
                }
                ++j; // anything else inside a tag: skip
            }
            i = j;
        }
        return out;
    }

    private static boolean isNameChar( final char c ) {
        return Character.isLetterOrDigit( c ) || ( c == '_' ) || ( c == ':' ) || ( c == '-' ) || ( c == '.' );
    }

    // ------------------------------------------------------------------ Newick / Nexus
    /**
     * Newick, and Nexus as Newick plus a small command tokenizer. Nexus is a sequence of {@code ";"}-terminated
     * COMMANDS ({@code Begin Taxa;}, {@code Dimensions NTax=9;}, {@code Tree x=(...);}), so a word is only read as
     * a keyword when it stands where a command can stand -- which is what keeps a taxon called "End" plain and
     * still tints the block name in {@code Begin Trees;} and the setting name in {@code Dimensions NTax=9;}.
     */
    private static List<Span> newickSpans( final String t, final boolean nexus ) {
        final List<Span> out = new ArrayList<>();
        final int n = t.length();
        int i = 0;
        boolean at_command_start = true;  // a command word may begin here (start of text, or after a ";")
        boolean expect_block_name = false; // the word after "Begin" is the block's name, whatever it is
        boolean in_tree_body = false;      // inside the Newick of a "Tree name=..." command
        String cmd = null;                 // the command being read, lower-cased ("dimensions", "tree", ...)
        while( i < n ) {
            final char c = t.charAt( i );
            if ( c == '\n' ) {
                // Only a finished statement lets the next line start a command: this way a wrapped TaxLabels list
                // or a multi-line tree cannot have one of its labels read as a keyword.
                at_command_start = ( cmd == null );
                ++i;
                continue;
            }
            if ( ( c == '\'' ) || ( c == '"' ) ) {
                // A quoted label: plain. Single quotes are the usual form ('' inside is an escaped quote); the
                // writer switches to DOUBLE quotes for a name that itself contains an apostrophe ("Seba's
                // short-tailed bat"), and a stray apostrophe read outside a quote would open one that swallows
                // everything up to the next apostrophe -- possibly the rest of the file.
                int j = i + 1;
                while( j < n ) {
                    if ( t.charAt( j ) == c ) {
                        if ( ( c == '\'' ) && ( j + 1 < n ) && ( t.charAt( j + 1 ) == '\'' ) ) {
                            j += 2;
                            continue;
                        }
                        break;
                    }
                    ++j;
                }
                i = Math.min( n, j + 1 );
                at_command_start = false;
                continue;
            }
            if ( c == '[' ) { // a comment / NHX block
                final int e = t.indexOf( ']', i + 1 );
                final int end = ( e < 0 ) ? n : e + 1;
                out.add( new Span( i, end, Kind.COMMENT ) );
                i = end;
                at_command_start = false;
                continue;
            }
            if ( ( c == '(' ) || ( c == ')' ) || ( c == ',' ) || ( c == ';' ) || ( c == ':' ) || ( c == '=' ) ) {
                out.add( new Span( i, i + 1, Kind.MARKUP ) );
                if ( c == ';' ) { // end of a Nexus command: the next word may open a new one
                    at_command_start = true;
                    expect_block_name = false;
                    in_tree_body = false;
                    cmd = null;
                    ++i;
                    continue;
                }
                if ( nexus && ( c == '=' ) && ( ( "tree".equals( cmd ) ) || ( "utree".equals( cmd ) ) ) ) {
                    in_tree_body = true; // everything up to the ";" is Newick, not Nexus
                }
                // a number right after ':' (a branch length), ')' (a support value) or, in Nexus, '=' (NTax=9)
                if ( ( c == ':' ) || ( c == ')' ) || ( nexus && ( c == '=' ) ) ) {
                    final int e = numberEnd( t, i + 1 );
                    if ( e > i + 1 ) {
                        out.add( new Span( i + 1, e, Kind.NUMBER ) );
                        i = e;
                        at_command_start = false;
                        continue;
                    }
                }
                ++i;
                at_command_start = false;
                continue;
            }
            if ( Character.isWhitespace( c ) ) {
                ++i;
                continue;
            }
            // a word (a label, or in Nexus possibly a keyword). A quote INSIDE a word (Foo"bar, O'Neil written
            // unquoted by some tool) is just a character of it: only a quote at a token boundary -- the branch
            // above -- opens a quoted label. Otherwise one stray quote would un-tint the rest of the file.
            int j = i;
            while( ( j < n ) && !Character.isWhitespace( t.charAt( j ) ) && "(),;:=[]".indexOf( t.charAt( j ) ) < 0 ) {
                ++j;
            }
            if ( nexus && !in_tree_body && ( j > i ) ) {
                final String w = t.substring( i, j ).toLowerCase( Locale.US );
                if ( expect_block_name ) { // "Begin Taxa;" -- the name IS the block, tint it whatever it is
                    out.add( new Span( i, j, Kind.KEYWORD ) );
                    expect_block_name = false;
                }
                else if ( at_command_start ) {
                    if ( NEXUS_COMMANDS.contains( w ) ) {
                        out.add( new Span( i, j, Kind.KEYWORD ) );
                    }
                    cmd = w;
                    at_command_start = false;
                    if ( "begin".equals( w ) ) {
                        expect_block_name = true;
                    }
                    else if ( "#nexus".equals( w ) ) {
                        cmd = null; // the file header stands alone: it carries no ";"
                        at_command_start = true;
                    }
                }
                else if ( NEXUS_SETTING_COMMANDS.contains( cmd ) && isFollowedByEquals( t, j ) ) {
                    out.add( new Span( i, j, Kind.ATTRIBUTE ) ); // "NTax" in "Dimensions NTax=9;"
                }
            }
            i = Math.max( j, i + 1 );
            at_command_start = false;
        }
        return out;
    }

    /** Whether the next non-space character at or after {@code from} is an "=" (a Nexus setting name). */
    private static boolean isFollowedByEquals( final String t, final int from ) {
        int j = from;
        while( ( j < t.length() ) && ( ( t.charAt( j ) == ' ' ) || ( t.charAt( j ) == '\t' ) ) ) {
            ++j;
        }
        return ( j < t.length() ) && ( t.charAt( j ) == '=' );
    }

    /** The end of a number (digits, sign, dot, exponent) starting at {@code from}, or {@code from} if none. */
    static int numberEnd( final String t, final int from ) {
        final int n = t.length();
        int j = from;
        if ( ( j < n ) && ( ( t.charAt( j ) == '-' ) || ( t.charAt( j ) == '+' ) ) ) {
            ++j;
        }
        int digits = 0;
        while( ( j < n ) && ( Character.isDigit( t.charAt( j ) ) || ( t.charAt( j ) == '.' ) ) ) {
            if ( Character.isDigit( t.charAt( j ) ) ) {
                ++digits;
            }
            ++j;
        }
        if ( digits == 0 ) {
            return from;
        }
        if ( ( j < n ) && ( ( t.charAt( j ) == 'e' ) || ( t.charAt( j ) == 'E' ) ) ) {
            int k = j + 1;
            if ( ( k < n ) && ( ( t.charAt( k ) == '-' ) || ( t.charAt( k ) == '+' ) ) ) {
                ++k;
            }
            int ed = 0;
            while( ( k < n ) && Character.isDigit( t.charAt( k ) ) ) {
                ++k;
                ++ed;
            }
            if ( ed > 0 ) {
                j = k;
            }
        }
        return j;
    }

    /** A file name for "Save As": the tree file's base name (else the tree name, else "tree") + the suffix. */
    static String suggestedFileName( final Phylogeny phy, final java.io.File tree_file, final Format format ) {
        String base = null;
        if ( tree_file != null ) {
            base = tree_file.getName();
            final int dot = base.lastIndexOf( '.' );
            if ( dot > 0 ) {
                base = base.substring( 0, dot );
            }
        }
        if ( ForesterUtil.isEmpty( base ) && ( phy != null ) && !ForesterUtil.isEmpty( phy.getName() ) ) {
            base = phy.getName().replaceAll( "[\\\\/:*?\"<>|\\s]+", "_" );
        }
        if ( ForesterUtil.isEmpty( base ) ) {
            base = "tree";
        }
        return base + format.suffix;
    }
}
