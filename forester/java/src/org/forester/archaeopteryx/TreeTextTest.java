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

import java.io.File;
import java.util.List;

import org.forester.archaeopteryx.TreeText.Format;
import org.forester.archaeopteryx.TreeText.Kind;
import org.forester.archaeopteryx.TreeText.Span;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;

/**
 * Headless tests for {@link TreeText}: rendering a tree to each format, the size line, the tint tokenizer for
 * phyloXML (tags, attributes, comments, declarations; attribute values and element text stay plain), Newick
 * (brackets and separators, branch lengths and support values, comments / NHX blocks, quoted labels untouched)
 * and Nexus (keywords), the number scanner, the tint size limit, and the suggested file name.
 */
public final class TreeTextTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreeText: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return render() && sizeText() && xmlSpans() && newickSpans() && nexusSpans() && numbers() && limit()
                    && fileName();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean render() {
        final Phylogeny phy = TreeFactsTest.fixture();
        phy.setName( "fixture" );
        final String xml = TreeText.render( phy, Format.PHYLOXML, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
        final String nwk = TreeText.render( phy, Format.NEWICK, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
        final String nex = TreeText.render( phy, Format.NEXUS, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
        if ( !xml.contains( "<phyloxml" ) || !xml.contains( "<name>fixture</name>" ) ) {
            return TestFail.here( xml.substring( 0, Math.min( 200, xml.length() ) ) );
        }
        if ( !nwk.startsWith( "(" ) || !nwk.contains( "B:0.3" ) || !nwk.endsWith( ";" ) ) {
            return TestFail.here( nwk );
        }
        if ( !nex.startsWith( "#NEXUS" ) || !nex.contains( "Begin Trees;" ) ) {
            return TestFail.here( nex );
        }
        // the support-value style is honoured (the unnamed node Y's bootstrap 70 becomes its label)
        final String named = TreeText.render( phy, Format.NEWICK,
                                              NH_CONVERSION_SUPPORT_VALUE_STYLE.AS_INTERNAL_NODE_NAMES );
        if ( !named.contains( ")70:" ) && !named.contains( ")70.0:" ) ) {
            return TestFail.here( named );
        }
        // no tree: ""
        if ( !"".equals( TreeText.render( null, Format.NEWICK, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ) )
                || !"".equals( TreeText.render( new Phylogeny(), Format.NEXUS, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ) ) ) {
            return TestFail.here();
        }
        if ( Format.PHYLOXML.wrapsByDefault() || !Format.NEWICK.wrapsByDefault() || !Format.NEXUS.wrapsByDefault() ) {
            return TestFail.here();
        }
        return true;
    }

    private static boolean sizeText() {
        if ( !"0 characters · 0 lines".equals( TreeText.sizeText( "" ) )
                || !"1 character · 1 line".equals( TreeText.sizeText( ";" ) )
                || !"5 characters · 1 line".equals( TreeText.sizeText( "abcd\n" ) )
                || !"7 characters · 3 lines".equals( TreeText.sizeText( "ab\ncd\n\n" ) )
                || !"1,200 characters · 3 lines".equals( TreeText.sizeText( "a\nb\n" + "x".repeat( 1196 ) ) ) ) {
            return TestFail.here( TreeText.sizeText( "ab\ncd\n\n" ) );
        }
        return true;
    }

    private static Span at( final List<Span> spans, final int offset ) {
        for( final Span s : spans ) {
            if ( ( offset >= s.start ) && ( offset < s.end ) ) {
                return s;
            }
        }
        return null;
    }

    private static boolean kindAt( final List<Span> spans, final String text, final String needle, final Kind expected ) {
        final int i = text.indexOf( needle );
        if ( i < 0 ) {
            return TestFail.here( "needle not in text: " + needle );
        }
        final Span s = at( spans, i );
        final Kind got = ( s == null ) ? null : s.kind;
        if ( got != expected ) {
            return TestFail.here( "at \"" + needle + "\": expected " + expected + ", got " + got );
        }
        return true;
    }

    private static boolean xmlSpans() {
        final String t = "<?xml version=\"1.0\"?>\n<!-- a comment -->\n<phyloxml xmlns=\"http://www.phyloxml.org\">\n"
                + "  <phylogeny rooted=\"true\" type='gene>tree'>\n    <name>Tree &amp; more</name>\n"
                + "    <clade branch_length=\"0.5\"/>\n  </phylogeny>\n</phyloxml>";
        final List<Span> spans = TreeText.spans( t, Format.PHYLOXML );
        boolean ok = true;
        ok &= kindAt( spans, t, "<?xml", Kind.COMMENT );
        ok &= kindAt( spans, t, "<!-- a", Kind.COMMENT );
        ok &= kindAt( spans, t, "<phyloxml", Kind.MARKUP );
        ok &= kindAt( spans, t, "xmlns=", Kind.ATTRIBUTE );
        ok &= kindAt( spans, t, "rooted=", Kind.ATTRIBUTE );
        ok &= kindAt( spans, t, "type=", Kind.ATTRIBUTE );
        ok &= kindAt( spans, t, "<name>", Kind.MARKUP );
        ok &= kindAt( spans, t, "</name>", Kind.MARKUP );
        ok &= kindAt( spans, t, "/>", Kind.MARKUP );
        ok &= kindAt( spans, t, "</phyloxml>", Kind.MARKUP );
        // plain: attribute values (even with a '>' inside a quote), element text
        ok &= kindAt( spans, t, "http://www.phyloxml.org", null );
        ok &= kindAt( spans, t, "gene>tree", null );
        ok &= kindAt( spans, t, "Tree &amp; more", null );
        ok &= kindAt( spans, t, "0.5", null );
        // the closing '>' of a tag with attributes is markup, and the quoted '>' did not end the tag early
        final int close = t.indexOf( "'>" ) + 1;
        if ( ( at( spans, close ) == null ) || ( at( spans, close ).kind != Kind.MARKUP ) ) {
            return TestFail.here( "closing bracket after a quoted '>'" );
        }
        // spans are ordered and non-overlapping
        for( int i = 1; i < spans.size(); ++i ) {
            if ( spans.get( i ).start < spans.get( i - 1 ).end ) {
                return TestFail.here( spans.get( i - 1 ) + " / " + spans.get( i ) );
            }
        }
        // an unterminated comment / tag does not throw
        if ( TreeText.spans( "<!-- open", Format.PHYLOXML ).isEmpty() || TreeText.spans( "<a b=\"x", Format.PHYLOXML ).isEmpty() ) {
            return TestFail.here();
        }
        return ok;
    }

    private static boolean newickSpans() {
        final String t = "((Homo_sapiens:0.12,'quoted (label):1'[&&NHX:S=x]:1e-3)90:0.5,Mus:-0.1)0.99:3;";
        final List<Span> spans = TreeText.spans( t, Format.NEWICK );
        boolean ok = true;
        ok &= kindAt( spans, t, "((", Kind.MARKUP );
        ok &= kindAt( spans, t, ":0.12", Kind.MARKUP );
        ok &= kindAt( spans, t, "0.12", Kind.NUMBER );
        ok &= kindAt( spans, t, "Homo_sapiens", null );
        ok &= kindAt( spans, t, "quoted (label):1", null ); // inside quotes: plain, brackets and colon ignored
        ok &= kindAt( spans, t, "[&&NHX", Kind.COMMENT );
        ok &= kindAt( spans, t, "1e-3", Kind.NUMBER );
        ok &= kindAt( spans, t, "90:", Kind.NUMBER ); // a support value right after ')'
        ok &= kindAt( spans, t, "Mus", null );
        ok &= kindAt( spans, t, "-0.1", Kind.NUMBER );
        ok &= kindAt( spans, t, "0.99", Kind.NUMBER );
        ok &= kindAt( spans, t, ";", Kind.MARKUP );
        // a label right after ')' (an internal node NAME) is not a number
        final String named = "(A,B)inner:1;";
        final List<Span> ns = TreeText.spans( named, Format.NEWICK );
        ok &= kindAt( ns, named, "inner", null );
        ok &= kindAt( ns, named, "1;", Kind.NUMBER );
        // Newick has no keywords, even for a word like 'tree'
        final String w = "(tree,begin);";
        ok &= kindAt( TreeText.spans( w, Format.NEWICK ), w, "tree", null );
        return ok;
    }

    private static boolean nexusSpans() {
        final String t = "#NEXUS\nBegin Taxa;\n Dimensions NTax=2;\n TaxLabels A B;\nEnd;\nBegin Trees;\n Tree "
                + "'my tree'=[&R](A:1,B:2);\nEnd;\n";
        final List<Span> spans = TreeText.spans( t, Format.NEXUS );
        boolean ok = true;
        ok &= kindAt( spans, t, "#NEXUS", Kind.KEYWORD );
        ok &= kindAt( spans, t, "Begin Taxa", Kind.KEYWORD );
        // the block name is part of the statement that opens the block, so it is tinted with it (it used to be
        // left plain, which made "Begin Taxa;" read as half markup, half data)
        ok &= kindAt( spans, t, "Taxa;", Kind.KEYWORD );
        ok &= kindAt( spans, t, "Trees;", Kind.KEYWORD );
        ok &= kindAt( spans, t, "Dimensions", Kind.KEYWORD );
        ok &= kindAt( spans, t, "NTax=2", Kind.ATTRIBUTE ); // a setting name, like an XML attribute
        ok &= kindAt( spans, t, "=2", Kind.MARKUP );
        ok &= kindAt( spans, t, "2;", Kind.NUMBER );        // its value
        ok &= kindAt( spans, t, "TaxLabels", Kind.KEYWORD );
        ok &= kindAt( spans, t, "A B", null );
        ok &= kindAt( spans, t, "End;", Kind.KEYWORD );
        ok &= kindAt( spans, t, "Tree '", Kind.KEYWORD );
        ok &= kindAt( spans, t, "my tree", null );
        ok &= kindAt( spans, t, "[&R]", Kind.COMMENT );
        ok &= kindAt( spans, t, "A:1", null );
        ok &= kindAt( spans, t, "1,B", Kind.NUMBER );
        // ---- a keyword is only a keyword where a COMMAND can stand ----
        // a taxon called "End" or "Matrix" in a label list is data, not markup
        final String labels = "Begin Taxa;\n TaxLabels End Matrix Format;\nEnd;\n";
        final List<Span> ls = TreeText.spans( labels, Format.NEXUS );
        ok &= kindAt( ls, labels, "End Matrix", null );
        ok &= kindAt( ls, labels, "Matrix Format", null );
        ok &= kindAt( ls, labels, "Format;", null );
        ok &= kindAt( ls, labels, "End;", Kind.KEYWORD ); // ... but the real one still is
        // a wrapped statement: the continuation line does not start a command, so a taxon there stays data even
        // when it is spelled like one (only a finished ";" statement lets the next line open a command)
        final String wrapped = "Begin Taxa;\n TaxLabels A\n Tree Matrix;\nEnd;\n";
        final List<Span> ws = TreeText.spans( wrapped, Format.NEXUS );
        ok &= kindAt( ws, wrapped, "Tree Matrix", null );
        ok &= kindAt( ws, wrapped, "Matrix;", null );
        ok &= kindAt( ws, wrapped, "End;", Kind.KEYWORD );
        // a label inside a tree is never a keyword either
        final String intree = "Begin Trees;\n Tree t1=(Matrix:1,End:2)Format;\nEnd;\n";
        final List<Span> is = TreeText.spans( intree, Format.NEXUS );
        ok &= kindAt( is, intree, "Matrix:1", null );
        ok &= kindAt( is, intree, "End:2", null );
        ok &= kindAt( is, intree, "Format;", null );
        // several commands on ONE line: the ";" is what opens the next one
        final String one = "Begin Taxa; Dimensions NTax=3; End;";
        final List<Span> os = TreeText.spans( one, Format.NEXUS );
        ok &= kindAt( os, one, "Taxa;", Kind.KEYWORD );
        ok &= kindAt( os, one, "Dimensions", Kind.KEYWORD );
        ok &= kindAt( os, one, "NTax=3", Kind.ATTRIBUTE );
        ok &= kindAt( os, one, "3;", Kind.NUMBER );
        ok &= kindAt( os, one, "End;", Kind.KEYWORD );
        // A name with an apostrophe is written in DOUBLE quotes ("Seba's short-tailed bat"). Read as a plain word
        // plus a stray apostrophe, it opened a single-quoted "label" that ran to the next apostrophe -- in the bat
        // demo that swallowed "End;", "Begin Trees;" and "Tree", which is exactly what the user saw.
        final String dq = "Begin Taxa;\n TaxLabels \"Seba's bat\" 'Large fox' \"Pallas's bat\";\nEnd;\nBegin Trees;\n"
                + " Tree t=[&R](\"Seba's bat\":1,('Large fox':2,\"Pallas's bat\":3)80:1);\nEnd;\n";
        final List<Span> ds = TreeText.spans( dq, Format.NEXUS );
        ok &= kindAt( ds, dq, "Seba's bat\" 'Large", null );
        ok &= kindAt( ds, dq, "Pallas's bat\";", null );
        ok &= kindAt( ds, dq, "End;\nBegin Trees", Kind.KEYWORD );
        ok &= kindAt( ds, dq, "Trees;", Kind.KEYWORD );
        ok &= kindAt( ds, dq, "Tree t=", Kind.KEYWORD );
        ok &= kindAt( ds, dq, "[&R]", Kind.COMMENT );
        ok &= kindAt( ds, dq, "Seba's bat\":1", null );
        ok &= kindAt( ds, dq, "1,(", Kind.NUMBER );
        ok &= kindAt( ds, dq, "80:1", Kind.NUMBER );
        ok &= kindAt( ds, dq, "End;\n", Kind.KEYWORD );
        if ( ds.isEmpty() || ( ds.get( ds.size() - 1 ).end < dq.lastIndexOf( "End" ) ) ) {
            return TestFail.here( "the last End; must still be reached: " + ds );
        }
        // a quote INSIDE a word is not a quote: a stray one must not un-tint everything after it
        final String stray = "Begin Taxa;\n TaxLabels Foo\"bar O'Neil baz;\nEnd;\nBegin Trees;\n Tree t=(a:1,b:2);\nEnd;\n";
        final List<Span> ss = TreeText.spans( stray, Format.NEXUS );
        ok &= kindAt( ss, stray, "Foo\"bar", null );
        ok &= kindAt( ss, stray, "O'Neil", null );
        ok &= kindAt( ss, stray, "End;\nBegin Trees", Kind.KEYWORD );
        ok &= kindAt( ss, stray, "Trees;", Kind.KEYWORD );
        ok &= kindAt( ss, stray, "Tree t=", Kind.KEYWORD );
        ok &= kindAt( ss, stray, "2)", Kind.NUMBER );
        // ... and the same in plain Newick, where the tree's own quote style is all there is
        final String dn = "(\"Seba's bat\":0.1,'it''s':0.2)90:1;";
        final List<Span> dns = TreeText.spans( dn, Format.NEWICK );
        ok &= kindAt( dns, dn, "Seba's", null );
        ok &= kindAt( dns, dn, "0.1", Kind.NUMBER );
        ok &= kindAt( dns, dn, "it''s", null );
        ok &= kindAt( dns, dn, "0.2", Kind.NUMBER );
        ok &= kindAt( dns, dn, "90:1", Kind.NUMBER );
        // a Format command's settings, and a non-numeric value left plain
        final String fmt = "Begin Data;\n Format DataType=DNA Gap=-;\nEnd;\n";
        final List<Span> fs = TreeText.spans( fmt, Format.NEXUS );
        ok &= kindAt( fs, fmt, "Data;", Kind.KEYWORD );
        ok &= kindAt( fs, fmt, "Format Data", Kind.KEYWORD );
        ok &= kindAt( fs, fmt, "DataType=", Kind.ATTRIBUTE );
        ok &= kindAt( fs, fmt, "DNA", null );
        ok &= kindAt( fs, fmt, "Gap=", Kind.ATTRIBUTE );
        return ok;
    }

    private static boolean numbers() {
        if ( ( TreeText.numberEnd( "0.12,", 0 ) != 4 ) || ( TreeText.numberEnd( "-1e-3)", 0 ) != 5 )
                || ( TreeText.numberEnd( "1E+10;", 0 ) != 5 ) || ( TreeText.numberEnd( "abc", 0 ) != 0 )
                || ( TreeText.numberEnd( "-", 0 ) != 0 ) || ( TreeText.numberEnd( ".", 0 ) != 0 )
                || ( TreeText.numberEnd( "5e", 0 ) != 1 ) || ( TreeText.numberEnd( "x:7", 2 ) != 3 )
                || ( TreeText.numberEnd( "", 0 ) != 0 ) ) {
            return TestFail.here( String.valueOf( TreeText.numberEnd( "-1e-3)", 0 ) ) );
        }
        return true;
    }

    private static boolean limit() {
        final String big = "<a>" + "x".repeat( TreeText.TINT_LIMIT ) + "</a>";
        if ( !TreeText.spans( big, Format.PHYLOXML ).isEmpty() || !TreeText.spans( null, Format.NEWICK ).isEmpty() ) {
            return TestFail.here( "over the limit / null: no spans" );
        }
        if ( TreeText.spans( "<a/>", Format.PHYLOXML ).isEmpty() || !TreeText.spans( "", Format.NEXUS ).isEmpty() ) {
            return TestFail.here();
        }
        return true;
    }

    private static boolean fileName() {
        final Phylogeny phy = new Phylogeny();
        phy.setName( "My tree / v2" );
        if ( !"bats.nwk".equals( TreeText.suggestedFileName( phy, new File( "/x/bats.xml" ), Format.NEWICK ) )
                || !"bats.tar.nex".equals( TreeText.suggestedFileName( phy, new File( "bats.tar.gz" ), Format.NEXUS ) )
                || !"My_tree_v2.xml".equals( TreeText.suggestedFileName( phy, null, Format.PHYLOXML ) )
                || !"tree.nex".equals( TreeText.suggestedFileName( new Phylogeny(), null, Format.NEXUS ) )
                || !"tree.nwk".equals( TreeText.suggestedFileName( null, null, Format.NEWICK ) )
                || !".xml".equals( Format.PHYLOXML.suffix ) ) {
            return TestFail.here( TreeText.suggestedFileName( phy, null, Format.PHYLOXML ) );
        }
        return true;
    }

    private TreeTextTest() {
    }
}
