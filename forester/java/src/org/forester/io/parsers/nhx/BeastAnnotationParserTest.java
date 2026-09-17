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

import java.io.File;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Headless unit tests for {@link BeastAnnotationParser}: the field mapping (posterior&rarr;confidence,
 * height/HPD&rarr;date, rate/traits&rarr;properties), the top-level tokenizer + interval/number helpers, and an
 * end-to-end parse through {@link NHXParser} (the wiring seam).
 */
public final class BeastAnnotationParserTest {

    public static boolean test() {
        return testFieldMapping() && testHelpers() && testMalformedTolerance() && testEndToEnd() && testMoreFields()
                && testLengthLedBlob() && testOptionOffKeepsComment() && testRoundTrip()
                && testQuotedCommaValueSurvivesScanner() && testUnquotedSpaceSurvivesInBracket()
                && testUnquotedSpaceStillDroppedInLabel() && testLegacyNhxTagDropsBothQuoteStyles()
                && testAuspiceVocabulary() && testBeastDatesUntouchedByAuspiceVocabulary()
                && testNumDateOutranksHeight() && testAuspiceVocabularyEdgeCases() && testBracketGroups()
                && testMrBayesLiteralLengthWins() && testFigTreeColouredBeastBlob() && testParseColor()
                && testHpdSpellingVariants() && testBootstrapAndProbKeys() && testUnderscoreReplacementSparesKeys()
                && testCommentTagMergesWithBlobProperties() && testSupportGroupBesideBlob()
                && testApostropheInsideBlobValueIsData() && testTaxonomyExtractionSurvivesABracket()
                && testSignedAndLeadingDotLengths() && testNumberGrammar() && testBlobQuoteRulesJoint()
                && testNhxTagSyntaxEitherPath() && testLibraryDefaultReadsAnnotations();
    }

    /** NHX tags are NHX tags whichever way the bracket-annotation option stands -- the genuine "[&&NHX:..]" and the
     *  sloppy spellings this parser has always forgiven (Test.testNHXParsing pins them; with the option on they used
     *  to reach the key=value parser and lose their species). A colon before the first '=' is what tells them from a
     *  blob: a BEAST key never contains one -- while a colon INSIDE a value does not make a blob NHX. */
    private static boolean testNhxTagSyntaxEitherPath() {
        final String sloppy = "(((((((A:0.2[&NHX:S=qw,erty]):0.2[&:S=u(io)p]):0.3[&NHX:S=asdf]):0.4[S=zxc]):0.5[]):0.6[&&NH:S=asd]):0.7[&&HX:S=za]):0.8[&&:S=zaq]";
        final String clean = "(((((((A:0.2[&&NHX:S=qw,erty]):0.2[&&NHX:S=u(io)p]):0.3[&&NHX:S=asdf]):0.4):0.5):0.6[&&NHX:S=asd]):0.7[&&NHX:S=za]):0.8[&&NHX:S=zaq]";
        try {
            for( final boolean on : new boolean[] { true, false } ) {
                final NHXParser p = new NHXParser();
                p.setParseBeastStyleExtendedTags( on );
                p.setSource( sloppy );
                final String got = p.parse()[ 0 ].toNewHampshireX();
                if ( !clean.equals( got ) ) {
                    return fail( "option " + ( on ? "ON" : "OFF" ) + ": sloppy NHX tags must still read as NHX, got " + got );
                }
            }
            // inside an NHX tag: unquoted white space is noise, a QUOTED run keeps its own (quotes protect, as in a
            // label) and may carry a colon -- and spaces BETWEEN the two ampersands change nothing
            final String[][] spaced = { { "(A[&&NHX:S=Homo sapiens]:1,B:1);", "Homosapiens" },
                    { "(A[&&NHX:S=\"homo sapiens\"]:1,B:1);", "homo sapiens" }, { "(A[&&NHX:S='homo sapiens']:1,B:1);", "homo sapiens" },
                    { "(A[&&NHX:S=\"a:b c\":D=Y]:1,B:1);", "a:b c" }, { "(A[ & & NHX : S = \"homo sapiens\" ]:1,B:1);", "homo sapiens" },
                    { "(A[& &NHX:S=x y]:1,B:1);", "xy" },
                    // padding at the ENDS of a quoted value is no more part of it than the quotes are (as Archaeopteryx.js)
                    { "(A[&&NHX:S=\" homo \"]:1,B:1);", "homo" }, { "(A[&&NHX:S=\" homo sapiens \":D=Y]:1,B:1);", "homo sapiens" },
                    { "(A[&&NHX:D=Y:S=' homo sapiens ']:1,B:1);", "homo sapiens" } };
            for( final String[] c : spaced ) {
                for( final boolean on : new boolean[] { true, false } ) {
                    final NHXParser sp = new NHXParser();
                    sp.setParseBeastStyleExtendedTags( on );
                    sp.setSource( c[ 0 ] );
                    final PhylogenyNode n = sp.parse()[ 0 ].getNode( "A" );
                    final String sn = n.getNodeData().isHasTaxonomy() ? n.getNodeData().getTaxonomy().getScientificName() : null;
                    if ( !c[ 1 ].equals( sn ) ) {
                        return fail( "option " + ( on ? "ON" : "OFF" ) + ": " + c[ 0 ] + " must give the species '" + c[ 1 ]
                                + "', got '" + sn + "'" );
                    }
                }
            }
            final String[][] fields = { { "S= homo ", "S=homo" }, { " S = a b ", "S=a b" }, { "D=Y", "D=Y" }, { " 0.5 ", "0.5" },
                    { "C= a=b ", "C=a=b" }, { "", "" } };
            for( final String[] f : fields ) {
                if ( !f[ 1 ].equals( NHXParser.trimNhxField( f[ 0 ] ) ) ) {
                    return fail( "trimNhxField('" + f[ 0 ] + "') must be '" + f[ 1 ] + "', got '" + NHXParser.trimNhxField( f[ 0 ] ) + "'" );
                }
            }
            final String[] nhx = { "&&NHX:S=x", "&NHX:S=x", "&:S=x", "&&NH:S=x", "&&:S=x", "&&NHX:D=Y:S=x" };
            for( final String g : nhx ) {
                if ( !NHXParser.isNhxTagSyntax( g ) ) {
                    return fail( "'" + g + "' is NHX tag syntax" );
                }
            }
            final String[] blobs = { "&rate=0.5", "&note=\"a:b\",x=1", "&date=2014-03-17,t=12:30", "&&NHX", "&", "91", "S=x:B=1",
                    "" };
            for( final String g : blobs ) {
                if ( NHXParser.isNhxTagSyntax( g ) ) {
                    return fail( "'" + g + "' is NOT NHX tag syntax" );
                }
            }
            // a colon inside a VALUE: still a blob, and the value keeps its colon
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A[&sampled=12:30,rate=0.5]:1.0,B:1.0);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            if ( ( prop( a, "beast:sampled" ) == null ) || !"12:30".equals( prop( a, "beast:sampled" ).getValue() )
                    || ( prop( a, "beast:rate" ) == null ) || ( a.getDistanceToParent() != 1.0 ) ) {
                return fail( "a colon inside a blob VALUE does not make the blob NHX" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "NHX-tag-syntax parse threw: " + e );
        }
    }

    /** The LIBRARY default: a parser nobody configured reads bracket annotations. It used to drop them -- which is
     *  how aptx_render came to draw BEAST trees without their posteriors -- and the Nexus parser follows the same
     *  constant. Off is still one call away, and still keeps the blob whole as a comment. */
    private static boolean testLibraryDefaultReadsAnnotations() {
        try {
            if ( !NHXParser.PARSE_BRACKET_ANNOTATIONS_DEFAULT ) {
                return fail( "the library default for bracket annotations is ON" );
            }
            final NHXParser p = new NHXParser();
            p.setSource( "(A[&posterior=0.9,rate=0.5]:1.0,B:1.0);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            if ( ( prop( a, "beast:rate" ) == null ) || !a.getBranchData().isHasConfidences() ) {
                return fail( "an unconfigured NHXParser must read bracket annotations" );
            }
            final org.forester.io.parsers.nexus.NexusPhylogeniesParser nex = new org.forester.io.parsers.nexus.NexusPhylogeniesParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                    .create( new java.io.ByteArrayInputStream( "#NEXUS\nbegin trees;\n\ttree t = [&R] (A[&rate=0.5]:1.0,B:1.0);\nend;\n"
                            .getBytes( java.nio.charset.StandardCharsets.UTF_8 ) ), nex )[ 0 ];
            if ( prop( phy.getNode( "A" ), "beast:rate" ) == null ) {
                return fail( "an unconfigured NexusPhylogeniesParser must read bracket annotations" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "library-default parse threw: " + e );
        }
    }

    /** The JOINT blob-quote cases, as Archaeopteryx.js runs them (its testBlobQuotes), through BOTH input paths of
     *  the streaming scanner (a char[] is scanned in place; a String, like a file, through a marked reader). The three halves
     *  of the rule -- where a quote may OPEN, where it may CLOSE, and the BOUND on the search for its partner --
     *  mask one another on the everyday cases, so each also gets a case where it is the only thing standing. */
    private static boolean testBlobQuoteRulesJoint() {
        final String ci = "C\u00f4te d'Ivoire";
        final String[][] cases = {
                // the odd count: used to throw
                { "(A:1[&country=" + ci + ",region=Africa],B:1);", "A=beast:country=" + ci + " | beast:region=Africa :1.0", "B= :1.0" },
                // the even count: used to read WRONG, silently, and lose tip B
                { "(A:1[&country=" + ci + "],B:2[&country=" + ci + "]);", "A=beast:country=" + ci + " :1.0", "B=beast:country=" + ci + " :2.0" },
                // a bare value that BEGINS with an apostrophe must not reach into the next tip's blob for a partner
                { "(A:1[&division='s-Hertogenbosch,region=Europe],B:2[&division='s-Gravenhage]);",
                        "A=beast:division='s-Hertogenbosch | beast:region=Europe :1.0", "B=beast:division='s-Gravenhage :2.0" },
                // OPEN alone: the second apostrophe stands where a value can end, so only "not after '='" keeps
                // these two fields apart
                { "(A:1[&a=rock'n,b=roll'],B:1);", "A=beast:a=rock'n | beast:b=roll' :1.0", "B= :1.0" },
                // CLOSE alone: the inner apostrophe is followed by a letter, so it is data and the comma after it
                // still belongs to the value
                { "(A:1[&note='a'b,c',x=1],B:1);", "A=beast:note=a'b,c | beast:x=1 :1.0", "B= :1.0" },
                // the BOUND alone: this value opens legitimately and never closes, and the next tip's blob ends in
                // an apostrophe that would make a fine partner
                { "(A:1[&division='s-Hertogenbosch],B:2[&k=v'],C:3);", "A=beast:division='s-Hertogenbosch :1.0", "B=beast:k=v' :2.0",
                        "C= :3.0" },
                // quoted values still work, either quote, the other one inside
                { "(A:1[&country=\"" + ci + "\"],B:1);", "A=beast:country=" + ci + " :1.0", "B= :1.0" },
                { "(A:1[&country='" + ci + "',region=Africa],B:1);", "A=beast:country=" + ci + " | beast:region=Africa :1.0", "B= :1.0" },
                // ... and what a quoted value protects is still protected
                { "(A:1[&k=\"a,b\",x=1],B:1);", "A=beast:k=a,b | beast:x=1 :1.0", "B= :1.0" },
                { "(A:1[&note=\"a]b\",x=1],B:1);", "A=beast:note=a]b | beast:x=1 :1.0", "B= :1.0" },
                // the scanner's CLOSE rule alone: closing at the inner apostrophe would let the ']' end the blob
                { "(A:1[&note='a'b]c',x=1],B:1);", "A=beast:note=a'b]c | beast:x=1 :1.0", "B= :1.0" },
                { "(A:1[&note=\"a[b:c\",x=1],B:1);", "A=beast:note=a[b:c | beast:x=1 :1.0", "B= :1.0" },
                { "(A:1[&loc.set={\"Hong Kong\",\"Korea, Republic of\"},x=1],B:1);",
                        "A=beast:loc_set={\"Hong Kong\",\"Korea, Republic of\"} | beast:x=1 :1.0", "B= :1.0" } };
        try {
            for( final String[] c : cases ) {
                for( final boolean from_file : new boolean[] { false, true } ) {
                    final NHXParser p = new NHXParser();
                    p.setParseBeastStyleExtendedTags( true );
                    p.setSource( from_file ? ( Object ) c[ 0 ] : ( Object ) c[ 0 ].toCharArray() );
                    final Phylogeny phy = p.parse()[ 0 ];
                    if ( phy.getNumberOfExternalNodes() != ( c.length - 1 ) ) {
                        return fail( ( from_file ? "[reader] " : "[char array] " ) + c[ 0 ] + ": expected " + ( c.length - 1 )
                                + " tips, got " + phy.getNumberOfExternalNodes() );
                    }
                    for( int i = 1; i < c.length; ++i ) {
                        final String name = c[ i ].substring( 0, c[ i ].indexOf( '=' ) );
                        final String got = name + "=" + describe( phy.getNode( name ) );
                        if ( !c[ i ].equals( got ) ) {
                            return fail( ( from_file ? "[reader] " : "[char array] " ) + c[ 0 ] + "\n      expected " + c[ i ]
                                    + "\n      got      " + got );
                        }
                    }
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "joint blob-quote parse threw: " + e );
        }
    }

    /** "ref=value | ref=value :length", in the order the properties were read. */
    private static String describe( final PhylogenyNode n ) {
        final StringBuilder sb = new StringBuilder();
        if ( n.getNodeData().getProperties() != null ) {
            for( final Property pr : n.getNodeData().getProperties().getProperties() ) {
                sb.append( sb.length() > 0 ? " | " : "" ).append( pr.getRef() ).append( '=' ).append( pr.getValue() );
            }
        }
        return sb + " :" + n.getDistanceToParent();
    }

    /** A number is a plain decimal with an optional exponent -- what Archaeopteryx.js reads as one. Java alone also
     *  takes "3f", "1d" and hex floats, which typed a clade called "3f" as xsd:decimal: invalid phyloXML. */
    private static boolean testNumberGrammar() {
        final String[] numbers = { "5", "5.", ".5", "-1e-3", "+2", "2.5E10", "0.0431319", "2006.349" };
        for( final String v : numbers ) {
            if ( BeastAnnotationParser.parseNumber( v ) == null ) {
                return fail( "'" + v + "' is a number" );
            }
        }
        final String[] not_numbers = { "3f", "1d", "2F", "0x1p3", "0x1F", " 5", "5 ", "NaN", "Infinity", "-Infinity", ".",
                "1e", "e5", "1,5", "", null };
        for( final String v : not_numbers ) {
            if ( BeastAnnotationParser.parseNumber( v ) != null ) {
                return fail( "'" + v + "' is not a number" );
            }
        }
        final PhylogenyNode n = new PhylogenyNode();
        BeastAnnotationParser.apply( "&clade=3f,rate=1e-3", n );
        if ( !"xsd:string".equals( prop( n, "beast:clade" ).getDataType() )
                || !"xsd:decimal".equals( prop( n, "beast:rate" ).getDataType() ) ) {
            return fail( "a clade called 3f is text; 1e-3 is a decimal" );
        }
        return true;
    }

    /** Increment D -- name-based taxonomy extraction used to be switched off by ANY bracket on the node. Only a
     *  genuine "[&&NHX" tag states the taxonomy itself; a BEAST blob, a "[91]" support or a kept comment does not. */
    private static boolean testTaxonomyExtractionSurvivesABracket() {
        try {
            for( final boolean option_on : new boolean[] { true, false } ) {
                final NHXParser p = new NHXParser();
                p.setParseBeastStyleExtendedTags( option_on );
                p.setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
                p.setSource( "((BCL2_HUMAN[&rate=0.5]:0.1,BCL2_MOUSE:0.2)BCL2_RAT:0.3[91],BCL2_CHICK:0.4[&&NHX:D=N],"
                        + "BCL2_XENLA:0.5);" );
                final Phylogeny phy = p.parse()[ 0 ];
                final String[][] expected = { { "BCL2_HUMAN", "HUMAN" }, { "BCL2_MOUSE", "MOUSE" },
                        { "BCL2_RAT", "RAT" }, { "BCL2_XENLA", "XENLA" } };
                for( final String[] e : expected ) {
                    final PhylogenyNode n = phy.getNode( e[ 0 ] );
                    if ( !n.getNodeData().isHasTaxonomy()
                            || !e[ 1 ].equals( n.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
                        return fail( "option " + ( option_on ? "ON" : "OFF" ) + ": " + e[ 0 ]
                                + " must get its taxonomy code from its name, bracket or not" );
                    }
                }
                // non-behaviour: a genuine &&NHX tag still switches the name-based extraction off
                if ( phy.getNode( "BCL2_CHICK" ).getNodeData().isHasTaxonomy() ) {
                    return fail( "option " + ( option_on ? "ON" : "OFF" )
                            + ": an &&NHX node states its own taxonomy; its name must not be mined" );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "taxonomy-extraction parse threw: " + e );
        }
    }

    /** Increment D -- a branch length may be signed and may start with its decimal point. ":-0.1" and the legal
     *  Newick ":.5" both used to vanish into "no length". Not clamped: a parser reports what the file says. */
    private static boolean testSignedAndLeadingDotLengths() {
        try {
            final NHXParser p = new NHXParser();
            p.setSource( "((A:-0.1,B:.5)AB:-1,(C:+0.25,D:-.5)CD:1e-3,E:-2.5E-2,F:7);" );
            final Phylogeny phy = p.parse()[ 0 ];
            final Object[][] expected = { { "A", -0.1 }, { "B", 0.5 }, { "AB", -1.0 }, { "C", 0.25 }, { "D", -0.5 },
                    { "CD", 0.001 }, { "E", -0.025 }, { "F", 7.0 } };
            for( final Object[] e : expected ) {
                final double d = phy.getNode( ( String ) e[ 0 ] ).getDistanceToParent();
                if ( Math.abs( d - ( ( Double ) e[ 1 ] ).doubleValue() ) > 1e-12 ) {
                    return fail( "length of " + e[ 0 ] + " must be " + e[ 1 ] + ", got " + d );
                }
            }
            // a negative length survives a Newick round trip
            final Phylogeny again = NHXParser.parse( phy.toNewHampshire() )[ 0 ];
            if ( again.getNode( "A" ).getDistanceToParent() != -0.1 ) {
                return fail( "a negative length must survive write + re-read, got "
                        + again.getNode( "A" ).getDistanceToParent() );
            }
            // non-behaviour: a sign with no number behind it is still not a length
            final NHXParser q = new NHXParser();
            q.setSource( "(A:-,B:+x,C:.,D:1);" );
            final Phylogeny junk = q.parse()[ 0 ];
            for( final String name : new String[] { "A", "B", "C" } ) {
                if ( junk.getNode( name ).getDistanceToParent() != org.forester.phylogeny.data.PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                    return fail( "'" + name + "' carries no number: no length, got "
                            + junk.getNode( name ).getDistanceToParent() );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "signed-length parse threw: " + e );
        }
    }

    /** An apostrophe INSIDE a blob value is data: "country=Côte d'Ivoire", 16 times in a real Nextstrain export
     *  (nextstrain_chikv_global_timetree.nexus), which could not be opened. Were it to open a quoted string, an ODD
     *  count fails the parse loudly and an EVEN count is worse -- two apostrophes pair up ACROSS tips and silently
     *  swallow every bracket, comma and paren between them, so a tip vanishes. A quote opens a string only where a
     *  value starts, which is where BEAST / FigTree / TreeTime put theirs. Found by the Archaeopteryx.js session. */
    private static boolean testApostropheInsideBlobValueIsData() {
        try {
            final String ci = "C\u00f4te d'Ivoire";
            // odd count
            NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A:1[&country=" + ci + ",region=Africa],B:1);" );
            Phylogeny phy = p.parse()[ 0 ];
            if ( phy.getNumberOfExternalNodes() != 2 ) {
                return fail( "odd apostrophe count: both tips must survive" );
            }
            Property country = prop( phy.getNode( "A" ), "beast:country" );
            Property region = prop( phy.getNode( "A" ), "beast:region" );
            if ( ( country == null ) || !ci.equals( country.getValue() ) || ( region == null )
                    || !"Africa".equals( region.getValue() ) ) {
                return fail( "an apostrophe inside a value must not swallow the next field, got country="
                        + ( country == null ? "null" : country.getValue() ) + " region="
                        + ( region == null ? "null" : region.getValue() ) );
            }
            // even count: the silent one
            p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A:1[&country=" + ci + "],B:2[&country=" + ci + "]);" );
            phy = p.parse()[ 0 ];
            if ( phy.getNumberOfExternalNodes() != 2 ) {
                return fail( "even apostrophe count: two apostrophes must not pair up across tips and swallow tip B" );
            }
            for( final String tip : new String[] { "A", "B" } ) {
                country = prop( phy.getNode( tip ), "beast:country" );
                if ( ( country == null ) || !ci.equals( country.getValue() ) ) {
                    return fail( "even apostrophe count: tip " + tip + " must carry country=" + ci + ", got "
                            + ( country == null ? "null" : country.getValue() ) );
                }
            }
            if ( phy.getNode( "B" ).getDistanceToParent() != 2 ) {
                return fail( "even apostrophe count: tip B keeps its length" );
            }
            // a stray double quote inside a value is data too; a set element may carry an apostrophe
            p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A:1[&size=5\" pipe,loc.set={" + ci + ",Ghana},n=1],B:1);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            if ( ( prop( a, "beast:size" ) == null ) || !"5\" pipe".equals( prop( a, "beast:size" ).getValue() )
                    || ( prop( a, "beast:loc_set" ) == null )
                    || !( "{" + ci + ",Ghana}" ).equals( prop( a, "beast:loc_set" ).getValue() )
                    || ( prop( a, "beast:n" ) == null ) ) {
                return fail( "a stray double quote / an apostrophe in a set element must stay data" );
            }
            // non-behaviour: a quote WHERE A VALUE STARTS still opens a string that protects , ] ( and an apostrophe
            p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A:1[&country=\"" + ci + "\",note='x,y(z)',m={'a,b',c}],B:1);" );
            final PhylogenyNode q = p.parse()[ 0 ].getNode( "A" );
            if ( ( prop( q, "beast:country" ) == null ) || !ci.equals( prop( q, "beast:country" ).getValue() ) ) {
                return fail( "a double-quoted value may carry an apostrophe" );
            }
            if ( ( prop( q, "beast:note" ) == null ) || !"x,y(z)".equals( prop( q, "beast:note" ).getValue() ) ) {
                return fail( "a single-quoted value still protects its comma and parens, got "
                        + ( prop( q, "beast:note" ) == null ? "null" : prop( q, "beast:note" ).getValue() ) );
            }
            if ( ( prop( q, "beast:m" ) == null ) || !"{'a,b',c}".equals( prop( q, "beast:m" ).getValue() ) ) {
                return fail( "a quoted set element still protects its comma" );
            }
            // non-behaviour: in a LABEL a quote still opens a string, as ever
            p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "('A (x,y)':1,B:1);" );
            if ( p.parse()[ 0 ].getNode( "A (x,y)" ) == null ) {
                return fail( "a quoted label is untouched by the blob rule" );
            }
            // the pure tokenizer obeys the same rule
            final List<String> t = BeastAnnotationParser.splitTopLevel( "country=" + ci + ",region=Africa" );
            if ( ( t.size() != 2 ) || !"region=Africa".equals( t.get( 1 ) ) ) {
                return fail( "splitTopLevel: an apostrophe inside a value must not hide the next comma, got " + t );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "apostrophe-in-blob parse threw: " + e );
        }
    }

    private static final String MRBAYES_NODE = "[&prob=0.95,prob_stddev=0.11,prob_range={0.9,1.0},"
            + "prob(percent)=\"95\",prob+-sd=\"95+-11\"]:0.05[&length_mean=0.0415,length_median=0.04129,"
            + "length_95%HPD={0.032,0.0503}]";

    /** Increment C -- the pure splitter: every top-level group in order, and the text around them kept together (a
     *  length BETWEEN two groups is the MrBayes case). An unbalanced bracket is an error, never a silent guess. */
    private static boolean testBracketGroups() {
        try {
            final StringBuilder outside = new StringBuilder();
            final List<String> g = NHXParser.bracketGroups( "A[&prob=1]:0.04[&length_mean=0.05]", outside );
            if ( ( g.size() != 2 ) || !"&prob=1".equals( g.get( 0 ) ) || !"&length_mean=0.05".equals( g.get( 1 ) ) ) {
                return fail( "bracketGroups must return both groups in order, got " + g );
            }
            if ( !"A:0.04".equals( outside.toString() ) ) {
                return fail( "the text outside the groups must survive together, got '" + outside + "'" );
            }
            final StringBuilder plain = new StringBuilder();
            if ( !NHXParser.bracketGroups( "A:0.1", plain ).isEmpty() || !"A:0.1".equals( plain.toString() ) ) {
                return fail( "no bracket: no groups, the text untouched" );
            }
            final StringBuilder trailing = new StringBuilder();
            if ( ( NHXParser.bracketGroups( "[&a=1]A:0.1[]", trailing ).size() != 2 )
                    || !"A:0.1".equals( trailing.toString() ) ) {
                return fail( "a leading group and an empty trailing group are both groups, got '" + trailing + "'" );
            }
        }
        catch ( final Exception e ) {
            return fail( "bracketGroups threw on well-formed input: " + e );
        }
        final String[] unbalanced = { "A[&prob=1:0.04", "A&prob=1]:0.04", "A[&a=1]:0.1]" };
        for( final String u : unbalanced ) {
            try {
                NHXParser.bracketGroups( u, new StringBuilder() );
                return fail( "an unbalanced bracket must be an error: " + u );
            }
            catch ( final NHXFormatException expected ) {
                // expected
            }
        }
        return true;
    }

    /** Decision 2 -- the file's literal ":length" is the branch length; MrBayes' length_mean / length_median /
     *  length_95%HPD are DATA about it. The fixture's literal (0.05) deliberately differs from its length_median
     *  (0.04129): only then can the two be told apart. With the option OFF the old single-group path is untouched
     *  (it never sees the literal, and still falls back on length_median). */
    private static boolean testMrBayesLiteralLengthWins() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(1" + MRBAYES_NODE + ",2" + MRBAYES_NODE.replace( ":0.05", ":0.07" ) + ");" );
            final Phylogeny phy = p.parse()[ 0 ];
            final PhylogenyNode n1 = phy.getNode( "1" );
            if ( n1.getDistanceToParent() != 0.05 ) {
                return fail( "the literal length between the two groups must win, got " + n1.getDistanceToParent() );
            }
            if ( phy.getNode( "2" ).getDistanceToParent() != 0.07 ) {
                return fail( "the second node's literal length must win too, got "
                        + phy.getNode( "2" ).getDistanceToParent() );
            }
            if ( n1.getBranchData().getNumberOfConfidences() != 1 ) {
                return fail( "prob must become exactly one confidence, got "
                        + n1.getBranchData().getNumberOfConfidences() );
            }
            final org.forester.phylogeny.data.Confidence c = n1.getBranchData().getConfidence( 0 );
            if ( ( c.getValue() != 0.95 ) || ( c.getStandardDeviation() != 0.11 )
                    || !"posterior probability".equals( c.getType() ) ) {
                return fail( "prob + prob_stddev -> posterior probability 0.95 sd 0.11, got " + c.getValue() + " sd "
                        + c.getStandardDeviation() + " type " + c.getType() );
            }
            final String[][] recovered = { { "beast:prob_range", "{0.9,1.0}" }, { "beast:prob_percent", "95" },
                    { "beast:prob_sd", "95+-11" }, { "beast:length_mean", "0.0415" },
                    { "beast:length_median", "0.04129" }, { "beast:length_95_HPD", "{0.032,0.0503}" } };
            for( final String[] r : recovered ) {
                final Property pr = prop( n1, r[ 0 ] );
                if ( ( pr == null ) || !r[ 1 ].equals( pr.getValue() ) ) {
                    return fail( "MrBayes field " + r[ 0 ] + " must be recovered as " + r[ 1 ] + ", got "
                            + ( pr == null ? "null" : pr.getValue() ) );
                }
            }
            if ( ( prop( n1, "beast:prob" ) != null ) || ( prop( n1, "beast:prob_stddev" ) != null ) ) {
                return fail( "prob / prob_stddev are the confidence, not ALSO properties" );
            }
            if ( n1.getNodeData().isHasDate() ) {
                return fail( "length_95%HPD is a BRANCH-LENGTH interval: it must never become a date" );
            }
            // a node that states NO literal length: length_median stands in, as the option-OFF path always did
            final NHXParser nolen = new NHXParser();
            nolen.setParseBeastStyleExtendedTags( true );
            nolen.setSource( "(1[&prob=1.0,prob_stddev=0,length_mean=0.05,length_median=0.04],2[&length_median=x]:0.1,3);" );
            final Phylogeny nl = nolen.parse()[ 0 ];
            if ( nl.getNode( "1" ).getDistanceToParent() != 0.04 ) {
                return fail( "no literal length: length_median must stand in, got " + nl.getNode( "1" ).getDistanceToParent() );
            }
            if ( ( nl.getNode( "2" ).getDistanceToParent() != 0.1 )
                    || ( nl.getNode( "3" ).getDistanceToParent() != org.forester.phylogeny.data.PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) ) {
                return fail( "a literal length is never replaced; a node with neither gets none" );
            }
            // option OFF: the old path, bit for bit
            final NHXParser off = new NHXParser();
            off.setParseBeastStyleExtendedTags( false );
            off.setSource( "(1" + MRBAYES_NODE + ",2" + MRBAYES_NODE + ");" );
            final PhylogenyNode o1 = off.parse()[ 0 ].getNode( "1" );
            if ( o1.getDistanceToParent() != 0.04129 ) {
                return fail( "option OFF keeps the old length_median fallback, got " + o1.getDistanceToParent() );
            }
            if ( ( o1.getBranchData().getNumberOfConfidences() != 1 )
                    || ( o1.getBranchData().getConfidence( 0 ).getValue() != 0.95 ) ) {
                return fail( "option OFF still reads MrBayes' prob" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "MrBayes parse threw: " + e );
        }
    }

    /** A FigTree-coloured BEAST tree LEADS its blob with "!color" -- which used to send the whole blob to a regex
     *  that recovered the colour alone (and not even that: FigTree writes a signed int, which the regex could not
     *  match). Everything in it must be read, wherever it stands. */
    private static boolean testFigTreeColouredBeastBlob() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "((A:1.0,B:1.0)[&!color=#-8381639,posterior=0.98,height=3.1,height_95%_HPD={2.9,3.4},"
                    + "rate=0.002]:0.5,C[&!color=#ff0000]:1.5);" );
            final Phylogeny phy = p.parse()[ 0 ];
            final PhylogenyNode ab = phy.getNode( "A" ).getParent();
            if ( ( ab.getBranchData().getBranchColor() == null )
                    || ( ab.getBranchData().getBranchColor().getValue().getRGB() != 0xFF801B39 ) ) {
                return fail( "!color=#-8381639 must decode to opaque 0x801B39, got "
                        + ab.getBranchData().getBranchColor() );
            }
            if ( !ab.getBranchData().isHasConfidences() || ( ab.getBranchData().getConfidence( 0 ).getValue() != 0.98 ) ) {
                return fail( "a !color-led blob must keep its posterior" );
            }
            final Date d = ab.getNodeData().getDate();
            if ( ( d == null ) || ( d.getValue() == null ) || ( d.getValue().doubleValue() != 3.1 ) || ( d.getMin() == null )
                    || ( d.getMin().doubleValue() != 2.9 ) || ( d.getMax().doubleValue() != 3.4 ) ) {
                return fail( "a !color-led blob must keep its height and HPD, got " + d );
            }
            if ( prop( ab, "beast:rate" ) == null ) {
                return fail( "a !color-led blob must keep its rate" );
            }
            if ( ab.getDistanceToParent() != 0.5 ) {
                return fail( "the length after a !color-led blob must survive, got " + ab.getDistanceToParent() );
            }
            final PhylogenyNode c = phy.getNode( "C" );
            if ( ( c.getBranchData().getBranchColor() == null )
                    || ( c.getBranchData().getBranchColor().getValue().getRGB() != 0xFFFF0000 ) ) {
                return fail( "!color=#ff0000 must still decode as hex red" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "FigTree-coloured parse threw: " + e );
        }
    }

    private static boolean testParseColor() {
        final java.awt.Color signed = BeastAnnotationParser.parseColor( "#-8381639" );
        if ( ( signed == null ) || ( signed.getRed() != 0x80 ) || ( signed.getGreen() != 0x1B )
                || ( signed.getBlue() != 0x39 ) || ( signed.getAlpha() != 255 ) ) {
            return fail( "#-8381639 is 0xFF801B39 -> opaque (128,27,57), got " + signed );
        }
        final java.awt.Color white = BeastAnnotationParser.parseColor( "#-1" );
        if ( ( white == null ) || ( white.getRGB() != 0xFFFFFFFF ) ) {
            return fail( "#-1 is opaque white, got " + white );
        }
        final java.awt.Color translucent = BeastAnnotationParser.parseColor( "#99999999999" ); // beyond an int
        if ( translucent != null ) {
            return fail( "a number beyond an int is not a colour" );
        }
        final java.awt.Color hex = BeastAnnotationParser.parseColor( "#00Ff7f" );
        if ( ( hex == null ) || ( hex.getRGB() != 0xFF00FF7F ) ) {
            return fail( "#00Ff7f is hex, either case, got " + hex );
        }
        final String[] not_colours = { "-8381639", "ff0000", "#", "#red", "#12345g", "#ff00", "#-", "#1.5", "", null };
        for( final String nc : not_colours ) {
            if ( BeastAnnotationParser.parseColor( nc ) != null ) {
                return fail( "'" + nc + "' is not a colour" );
            }
        }
        // what is not a colour is kept as data, never thrown away -- and never aborts the parse
        final PhylogenyNode n = new PhylogenyNode();
        BeastAnnotationParser.apply( "&!color=crimson,rate=0.5", n );
        if ( ( n.getBranchData().getBranchColor() != null ) || ( prop( n, "beast:_color" ) == null )
                || ( prop( n, "beast:rate" ) == null ) ) {
            return fail( "an unreadable colour stays a property and the rest of the blob is still read" );
        }
        // ...as TEXT, even when it looks like a number (a colour with no '#' is refused -- Christian, 2026-09-17): a
        // "!" key is a FigTree display directive, never a measurement for Color-by to draw a gradient of. A user's
        // own trait called "color" is a trait like any other.
        final PhylogenyNode bare = new PhylogenyNode();
        BeastAnnotationParser.apply( "&!color=-8381639,!rotate=1,color=3", bare );
        if ( ( bare.getBranchData().getBranchColor() != null ) || ( prop( bare, "beast:_color" ) == null )
                || !"-8381639".equals( prop( bare, "beast:_color" ).getValue() )
                || !"xsd:string".equals( prop( bare, "beast:_color" ).getDataType() )
                || !"xsd:string".equals( prop( bare, "beast:_rotate" ).getDataType() ) ) {
            return fail( "a refused !color (and any ! directive) is kept as text, never typed numeric" );
        }
        if ( !"xsd:decimal".equals( prop( bare, "beast:color" ).getDataType() ) ) {
            return fail( "a plain trait called color=3 is an ordinary numeric trait" );
        }
        return true;
    }

    /** The HPD key is spelled three ways across BEAST 1 / BEAST 2 / MrBayes, in any case; all are the node-age
     *  interval. A length_... interval is NOT: it describes the branch, and must never bracket a date. */
    private static boolean testHpdSpellingVariants() {
        final String[] spellings = { "height_95%_HPD", "height_95%HPD", "height95%HPD", "HEIGHT_95%_hpd" };
        for( final String sp : spellings ) {
            final PhylogenyNode n = new PhylogenyNode();
            BeastAnnotationParser.apply( "&height=3.1," + sp + "={2.9,3.4}", n );
            final Date d = n.getNodeData().getDate();
            if ( ( d == null ) || ( d.getMin() == null ) || ( d.getMin().doubleValue() != 2.9 ) || ( d.getMax() == null )
                    || ( d.getMax().doubleValue() != 3.4 ) ) {
                return fail( sp + " must fill the date interval, got " + d );
            }
        }
        final String[] length_spellings = { "length_95%HPD", "length_95%_HPD", "length_range" };
        for( final String sp : length_spellings ) {
            final PhylogenyNode n = new PhylogenyNode();
            BeastAnnotationParser.apply( "&height=3.1," + sp + "={0.1,0.2}", n );
            final Date d = n.getNodeData().getDate();
            if ( ( d == null ) || ( d.getMin() != null ) || ( d.getMax() != null ) ) {
                return fail( sp + " is a branch-length interval and must NOT fill the date interval, got " + d );
            }
            if ( n.getNodeData().getProperties() == null ) {
                return fail( sp + " must be kept as a property" );
            }
        }
        if ( !"height95hpd".equals( BeastAnnotationParser.normalizedKey( "Height_95%_HPD" ) )
                || !"length95hpd".equals( BeastAnnotationParser.normalizedKey( "length_95%HPD" ) )
                || !"location195hpd1".equals( BeastAnnotationParser.normalizedKey( "location1_95%HPD_1" ) ) ) {
            return fail( "normalizedKey: lower-case, no % and no _" );
        }
        return true;
    }

    private static boolean testBootstrapAndProbKeys() {
        // the last one leads EVERY field with '&' (Test.testNHXNodeParsing2 has carried that shape for years)
        final String[] boots = { "&bootstrap=87", "&boot=87", "&rate=1,Bootstrap=87", "&rate=1,&bootstrap=87" };
        for( final String b : boots ) {
            final PhylogenyNode n = new PhylogenyNode();
            BeastAnnotationParser.apply( b, n );
            if ( ( n.getBranchData().getNumberOfConfidences() != 1 )
                    || ( n.getBranchData().getConfidence( 0 ).getValue() != 87 )
                    || !"bootstrap".equals( n.getBranchData().getConfidence( 0 ).getType() ) ) {
                return fail( b + " must become one confidence of type bootstrap" );
            }
        }
        final PhylogenyNode amp = new PhylogenyNode();
        BeastAnnotationParser.apply( "&bootstrap=69,&!color=#FFFFFF,&rate=0.5", amp );
        if ( ( amp.getBranchData().getBranchColor() == null ) || ( prop( amp, "beast:rate" ) == null )
                || ( amp.getNodeData().getProperties().size() != 1 ) ) {
            return fail( "a field's own leading '&' is punctuation, not part of its key" );
        }
        final PhylogenyNode neg = new PhylogenyNode();
        BeastAnnotationParser.apply( "&bootstrap=-1,prob=-0.5", neg );
        if ( neg.getBranchData().isHasConfidences() ) {
            return fail( "a negative bootstrap / prob is not a confidence" );
        }
        if ( prop( neg, "beast:bootstrap" ) == null ) {
            return fail( "a negative bootstrap is kept as data" );
        }
        // prob alone: no standard deviation invented; the type is the one MrBayes' prob has always had here
        final PhylogenyNode alone = new PhylogenyNode();
        BeastAnnotationParser.apply( "&prob=0.81", alone );
        if ( ( alone.getBranchData().getNumberOfConfidences() != 1 )
                || !BeastAnnotationParser.MRBAYES_CONFIDENCE_TYPE.equals( alone.getBranchData().getConfidence( 0 ).getType() )
                || ( alone.getBranchData().getConfidence( 0 ).getStandardDeviation() == 0.11 ) ) {
            return fail( "prob alone -> one posterior-probability confidence" );
        }
        // the deviation may come BEFORE the probability
        final PhylogenyNode reversed = new PhylogenyNode();
        BeastAnnotationParser.apply( "&prob_stddev=0.11,prob=0.95", reversed );
        if ( ( reversed.getBranchData().getNumberOfConfidences() != 1 )
                || ( reversed.getBranchData().getConfidence( 0 ).getStandardDeviation() != 0.11 ) ) {
            return fail( "prob_stddev before prob must still qualify it" );
        }
        final PhylogenyNode sd_only = new PhylogenyNode();
        BeastAnnotationParser.apply( "&prob_stddev=0.11", sd_only );
        if ( sd_only.getBranchData().isHasConfidences() || ( prop( sd_only, "beast:prob_stddev" ) == null ) ) {
            return fail( "a deviation with no probability is kept as data, not made a confidence" );
        }
        return true;
    }

    /** "Replace underscores" is about LABELS. With the option on it used to run over the whole annotation text, so
     *  "height_median" arrived as "height median" and the age was lost. */
    private static boolean testUnderscoreReplacementSparesKeys() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setReplaceUnderscores( true );
            p.setSource( "(Homo_sapiens[&height_median=1.2,my_trait=a_b]:1.0,B:1.0);" );
            final Phylogeny phy = p.parse()[ 0 ];
            final PhylogenyNode h = phy.getNode( "Homo sapiens" );
            final Date d = h.getNodeData().getDate();
            if ( ( d == null ) || ( d.getValue() == null ) || ( d.getValue().doubleValue() != 1.2 ) ) {
                return fail( "replacing underscores in labels must not break the key height_median, got " + d );
            }
            final Property t = prop( h, "beast:my_trait" );
            if ( ( t == null ) || !"a_b".equals( t.getValue() ) ) {
                return fail( "a value inside the blob keeps its underscores, got " + ( t == null ? "null" : t.getValue() ) );
            }
            // ...while a legacy NHX tag value is replaced like the label, on BOTH option paths, as it always was
            for( final boolean on : new boolean[] { true, false } ) {
                final NHXParser q = new NHXParser();
                q.setParseBeastStyleExtendedTags( on );
                q.setReplaceUnderscores( true );
                q.setSource( "(A_x[&&NHX:S=Homo_sapiens]:1,B:1);" );
                final PhylogenyNode ax = q.parse()[ 0 ].getNode( "A x" );
                if ( !ax.getNodeData().isHasTaxonomy() || !"Homo sapiens".equals( ax.getNodeData().getTaxonomy().getScientificName() ) ) {
                    return fail( "option " + ( on ? "ON" : "OFF" ) + ": Replace Underscores must reach an NHX tag value (S=Homo_sapiens)" );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "replace-underscores parse threw: " + e );
        }
    }

    /** A legacy &&NHX group beside a single-& blob on one node: both are read, and the C= comment is ADDED to the
     *  properties the blob produced (it used to replace the whole list) -- with its first letter. */
    private static boolean testCommentTagMergesWithBlobProperties() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A[&rate=0.5]:0.1[&&NHX:C=hello:S=homo],B:0.2);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            if ( prop( a, "beast:rate" ) == null ) {
                return fail( "the C= comment must not wipe the blob's properties" );
            }
            final Property c = prop( a, "nh:comment" );
            if ( ( c == null ) || !"hello".equals( c.getValue() ) ) {
                return fail( "C=hello is the comment 'hello', got " + ( c == null ? "null" : c.getValue() ) );
            }
            if ( !a.getNodeData().isHasTaxonomy() || !"homo".equals( a.getNodeData().getTaxonomy().getScientificName() ) ) {
                return fail( "the &&NHX group after the length must be read too" );
            }
            if ( a.getDistanceToParent() != 0.1 ) {
                return fail( "the length between the blob and the &&NHX group must survive" );
            }
            // option OFF: a bracket inside a quoted value comes back as itself in the kept comment
            final NHXParser off_b = new NHXParser();
            off_b.setParseBeastStyleExtendedTags( false );
            off_b.setSource( "(A[&note=\"a]b[c:d\"]:0.1,B:0.2);" );
            final Property ob = prop( off_b.parse()[ 0 ].getNode( "A" ), "nh:comment" );
            if ( ( ob == null ) || !"note=\"a]b[c:d\"".equals( ob.getValue() ) ) {
                return fail( "option OFF: brackets and a colon inside a quoted value survive in the comment, got "
                        + ( ob == null ? "null" : ob.getValue() ) );
            }
            // option OFF: the kept blob comment is what it always was
            final NHXParser off = new NHXParser();
            off.setParseBeastStyleExtendedTags( false );
            off.setSource( "(A[&rate=0.5]:0.1,B:0.2);" );
            final Property oc = prop( off.parse()[ 0 ].getNode( "A" ), "nh:comment" );
            if ( ( oc == null ) || !"rate=0.5".equals( oc.getValue() ) ) {
                return fail( "option OFF: the blob is kept as the comment 'rate=0.5', got "
                        + ( oc == null ? "null" : oc.getValue() ) );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "comment-merge parse threw: " + e );
        }
    }

    /** A bare numeric "[91]" support group beside a blob, either side of the length: both are read. */
    private static boolean testSupportGroupBesideBlob() {
        try {
            final String[] shapes = { "((A:1,B:1):0.5[91][&rate=2],C:1);", "((A:1,B:1)[&rate=2]:0.5[91],C:1);" };
            for( final String shape : shapes ) {
                final NHXParser p = new NHXParser();
                p.setParseBeastStyleExtendedTags( true );
                p.setSource( shape );
                final PhylogenyNode ab = p.parse()[ 0 ].getNode( "A" ).getParent();
                if ( !ab.getBranchData().isHasConfidences() || ( ab.getBranchData().getConfidence( 0 ).getValue() != 91 ) ) {
                    return fail( shape + ": the [91] support must be read" );
                }
                if ( prop( ab, "beast:rate" ) == null ) {
                    return fail( shape + ": the blob beside the support must be read" );
                }
                if ( ab.getDistanceToParent() != 0.5 ) {
                    return fail( shape + ": the length must survive, got " + ab.getDistanceToParent() );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "support-beside-blob parse threw: " + e );
        }
    }

    /** Increment B -- the Auspice "download Nexus" vocabulary, in the real file's shape (the blob BEFORE the length,
     *  NODE_ internals), must land where AuspiceJsonParser puts the same dataset's JSON: a date with unit "year",
     *  nextstrain:num_date, nextstrain:div. The calendar-axis assertion is the one that proves the actual bug fixed:
     *  a date with no unit derives NO time axis, however right its value is. */
    private static boolean testAuspiceVocabulary() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "((A[&num_date=2006.349,num_date_CI={2006.1,2006.6},clade=H1,country=China,div=0.0431319]:0.5,"
                    + "B[&num_date=2006.849,num_date_CI={2006.7,2006.9},clade=H1,div=0.05]:1.0)"
                    + "NODE_0000001[&num_date=2005.849,num_date_CI={2005.2,2006.0},div=0.04]:1.0,"
                    + "C[&num_date=2005.5,div=0.01]:0.651)NODE_0000000[&num_date=2004.849,div=0];" );
            final Phylogeny phy = p.parse()[ 0 ];
            final PhylogenyNode a = phy.getNode( "A" );
            final Date d = a.getNodeData().getDate();
            if ( ( d == null ) || ( d.getValue() == null ) || ( d.getValue().doubleValue() != 2006.349 ) ) {
                return fail( "num_date must become the date value, got " + d );
            }
            if ( !BeastAnnotationParser.YEAR_UNIT.equals( d.getUnit() ) || !"year".equals( d.getUnit() ) ) {
                return fail( "a num_date is in calendar years: unit must be \"year\", got '" + d.getUnit() + "'" );
            }
            // the interval: on the tip and on the internal node alike
            if ( ( d.getMin() == null ) || ( d.getMin().doubleValue() != 2006.1 ) || ( d.getMax() == null )
                    || ( d.getMax().doubleValue() != 2006.6 ) ) {
                return fail( "num_date_CI must become the date's min/max, got " + d.getMin() + " .. " + d.getMax() );
            }
            final Date di = a.getParent().getNodeData().getDate();
            if ( ( di == null ) || ( di.getMin() == null ) || ( di.getMin().doubleValue() != 2005.2 ) || ( di.getMax() == null )
                    || ( di.getMax().doubleValue() != 2006.0 ) ) {
                return fail( "an internal node's num_date_CI must become its date's min/max, got " + di );
            }
            if ( !isEmpty( d.getDesc() ) ) {
                return fail( "no date= field, so no desc, got '" + d.getDesc() + "'" );
            }
            final Property div = prop( a, "nextstrain:div" );
            if ( ( div == null ) || !"0.0431319".equals( div.getValue() ) || !"xsd:decimal".equals( div.getDataType() ) ) {
                return fail( "div must become the numeric nextstrain:div (the Time|Div ref), got " + div );
            }
            final Property nd = prop( a, "nextstrain:num_date" );
            if ( ( nd == null ) || !"2006.349".equals( nd.getValue() ) || !"xsd:decimal".equals( nd.getDataType() ) ) {
                return fail( "num_date must ALSO be the numeric nextstrain:num_date property, got " + nd );
            }
            if ( ( prop( a, "beast:num_date" ) != null ) || ( prop( a, "beast:div" ) != null )
                    || ( prop( a, "beast:num_date_CI" ) != null ) ) {
                return fail( "the Auspice vocabulary must not ALSO be left behind as opaque beast: properties" );
            }
            final Property clade = prop( a, "beast:clade" );
            if ( ( clade == null ) || !"H1".equals( clade.getValue() ) ) {
                return fail( "an ordinary trait beside the Auspice vocabulary still becomes a beast: property" );
            }
            if ( a.getDistanceToParent() != 0.5 ) {
                return fail( "the length after the blob must survive, got " + a.getDistanceToParent() );
            }
            // a node with a num_date but no CI: a value and a unit, no interval
            final Date dc = phy.getNode( "C" ).getNodeData().getDate();
            if ( ( dc == null ) || ( dc.getMin() != null ) || ( dc.getMax() != null ) || !"year".equals( dc.getUnit() ) ) {
                return fail( "a num_date without a CI is a bare value with unit year, got " + dc );
            }
            if ( org.forester.archaeopteryx.AptxUtil
                    .deriveTimeAxisType( phy ) != org.forester.archaeopteryx.Options.TIME_AXIS_TYPE.CALENDAR ) {
                return fail( "an Auspice Nexus tree must derive the CALENDAR axis, got "
                        + org.forester.archaeopteryx.AptxUtil.deriveTimeAxisType( phy ) );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "Auspice vocabulary parse threw: " + e );
        }
    }

    /** Non-behaviour: what BEAST trees did before Increment B they still do. A date= beside a height* stays the
     *  desc, the height stays the value, and an age gets NO unit -- whether the date is a calendar string
     *  (test_trees/peartree_example.nexus) or a bare decimal year (a tip-dated BEAST tree: promoting it would put the
     *  tips on calendar years and the internals on ages, one axis over two scales). */
    private static boolean testBeastDatesUntouchedByAuspiceVocabulary() {
        final String[][] cases = { { "&height_median=1.2,date=\"2014-03-17\"", "1.2", "2014-03-17" },
                { "&height=3.1,date=2014.21", "3.1", "2014.21" } };
        for( final String[] c : cases ) {
            final PhylogenyNode n = new PhylogenyNode();
            BeastAnnotationParser.apply( c[ 0 ], n );
            final Date d = n.getNodeData().getDate();
            if ( ( d == null ) || ( d.getValue() == null ) || !c[ 1 ].equals( d.getValue().toPlainString() ) ) {
                return fail( c[ 0 ] + ": the height must stay the date value, got " + d );
            }
            if ( !c[ 2 ].equals( d.getDesc() ) ) {
                return fail( c[ 0 ] + ": date= must stay the desc, got '" + d.getDesc() + "'" );
            }
            if ( !isEmpty( d.getUnit() ) ) {
                return fail( c[ 0 ] + ": an age has no unit, got '" + d.getUnit() + "'" );
            }
        }
        return true;
    }

    /** Precedence: a num_date outranks every height as the date value AND brings its own interval -- a height HPD is
     *  on the age scale, so it must never bracket a calendar year. The age fields are consumed, as before. */
    private static boolean testNumDateOutranksHeight() {
        final PhylogenyNode n = new PhylogenyNode();
        BeastAnnotationParser.apply( "&height=3.1,height_95%_HPD={2.9,3.4},num_date=2011.5,num_date_CI={2011.2,2011.8}",
                                     n );
        final Date d = n.getNodeData().getDate();
        if ( ( d == null ) || ( d.getValue() == null ) || ( d.getValue().doubleValue() != 2011.5 ) ) {
            return fail( "num_date must outrank height as the date value, got " + d );
        }
        if ( ( d.getMin() == null ) || ( d.getMin().doubleValue() != 2011.2 ) || ( d.getMax() == null )
                || ( d.getMax().doubleValue() != 2011.8 ) ) {
            return fail( "a num_date's interval is its own CI, never the height HPD, got " + d.getMin() + " .. "
                    + d.getMax() );
        }
        if ( !"year".equals( d.getUnit() ) ) {
            return fail( "num_date brings unit year even beside a height, got '" + d.getUnit() + "'" );
        }
        // the heights it outranked are kept as what they were written as (as Archaeopteryx.js keeps them)
        if ( ( prop( n, "beast:height" ) == null ) || !"3.1".equals( prop( n, "beast:height" ).getValue() )
                || ( prop( n, "beast:height_95_HPD" ) == null )
                || !"{2.9,3.4}".equals( prop( n, "beast:height_95_HPD" ).getValue() ) ) {
            return fail( "a height outranked by a num_date is kept as a beast: property, not dropped" );
        }
        final PhylogenyNode plain = new PhylogenyNode();
        BeastAnnotationParser.apply( "&height=3.1,height_95%_HPD={2.9,3.4}", plain );
        if ( plain.getNodeData().getProperties() != null ) {
            return fail( "a height that IS the date is consumed by it, as ever: no property" );
        }
        // no CI of its own: the height HPD must still not be borrowed
        final PhylogenyNode m = new PhylogenyNode();
        BeastAnnotationParser.apply( "&height=3.1,height_95%_HPD={2.9,3.4},num_date=2011.5", m );
        final Date dm = m.getNodeData().getDate();
        if ( ( dm == null ) || ( dm.getMin() != null ) || ( dm.getMax() != null ) ) {
            return fail( "a num_date without its own CI must not borrow the height HPD, got " + dm );
        }
        return true;
    }

    /** Edge cases: an unparseable num_date / div is kept as plain data, never a date or a numeric ref; a CI with no
     *  point date to hang on is kept as a property; TreeTime's mutations / mcc are text even when they look numeric. */
    private static boolean testAuspiceVocabularyEdgeCases() {
        final PhylogenyNode bad = new PhylogenyNode();
        BeastAnnotationParser.apply( "&num_date=unknown,div=NaN", bad );
        if ( bad.getNodeData().isHasDate() ) {
            return fail( "an unparseable num_date must not create a date" );
        }
        final Property bad_nd = prop( bad, "beast:num_date" );
        if ( ( bad_nd == null ) || !"unknown".equals( bad_nd.getValue() ) || !"xsd:string".equals( bad_nd.getDataType() ) ) {
            return fail( "an unparseable num_date is kept as a plain text property, got " + bad_nd );
        }
        if ( ( prop( bad, "nextstrain:div" ) != null ) || ( prop( bad, "nextstrain:num_date" ) != null ) ) {
            return fail( "an unparseable value must never reach the numeric nextstrain: refs" );
        }
        final PhylogenyNode ci_only = new PhylogenyNode();
        BeastAnnotationParser.apply( "&num_date_CI={2011.2,2011.8}", ci_only );
        if ( ci_only.getNodeData().isHasDate() ) {
            return fail( "a CI with no num_date must not create a date" );
        }
        final Property ci = prop( ci_only, "nextstrain:num_date_CI" );
        if ( ( ci == null ) || !"{2011.2,2011.8}".equals( ci.getValue() ) ) {
            return fail( "a CI with no num_date is kept as data, got " + ci );
        }
        final PhylogenyNode tt = new PhylogenyNode();
        BeastAnnotationParser.apply( "&mutations=\"123\",mcc=4,rate=0.5", tt );
        final Property mut = prop( tt, "beast:mutations" );
        final Property mcc = prop( tt, "beast:mcc" );
        final Property rate = prop( tt, "beast:rate" );
        if ( ( mut == null ) || !"xsd:string".equals( mut.getDataType() ) || ( mcc == null )
                || !"xsd:string".equals( mcc.getDataType() ) ) {
            return fail( "mutations / mcc are text even when they look numeric" );
        }
        if ( ( rate == null ) || !"xsd:decimal".equals( rate.getDataType() ) ) {
            return fail( "an ordinary numeric field stays xsd:decimal" );
        }
        return true;
    }

    private static boolean isEmpty( final String s ) {
        return org.forester.util.ForesterUtil.isEmpty( s );
    }

    /** Non-behaviour: quote preservation is for single-'&' data blobs only. Inside a legacy "[&&NHX:...]" tag a quote
     *  character stays formatting noise, for BOTH quote styles and for the opening AND the closing quote alike -- an
     *  opening single quote was once kept while its closing partner was dropped, leaving S='homo as the name. */
    private static boolean testLegacyNhxTagDropsBothQuoteStyles() {
        try {
            final String[] quoted = { "'homo'", "\"homo\"" };
            for( final String q : quoted ) {
                final NHXParser p = new NHXParser();
                p.setParseBeastStyleExtendedTags( true );
                p.setSource( "(A:1.0[&&NHX:S=" + q + "],B:1.0);" );
                final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
                if ( !a.getNodeData().isHasTaxonomy() ) {
                    return fail( "legacy tag S=" + q + " must still set the taxonomy" );
                }
                final String sn = a.getNodeData().getTaxonomy().getScientificName();
                if ( !"homo".equals( sn ) ) {
                    return fail( "quotes inside a legacy &&NHX tag must be dropped (S=" + q + "), got '" + sn + "'" );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "legacy-tag quote parse threw: " + e );
        }
    }

    /** Increment A -- scanner bracket fidelity: a quoted value containing a top-level comma survives whole through
     *  NHXParser's streaming scanner into BeastAnnotationParser. Before this fix the scanner itself consumed the
     *  quote characters, so by the time BeastAnnotationParser.splitTopLevel saw the blob the embedded commas looked
     *  like ordinary field separators -- every comma-separated piece after the first (here T92C and G100A) was
     *  silently dropped because it has no '='. */
    private static boolean testQuotedCommaValueSurvivesScanner() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A[&mutations=\"A54G,T92C,G100A\"]:1.0,B:1.0);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            // treetime: -- a tree with mutations and no node age is TreeTime's (BracketAnnotationNormalizer)
            final Property mutations = prop( a, "treetime:mutations" );
            if ( ( mutations == null ) || !"A54G,T92C,G100A".equals( mutations.getValue() ) ) {
                return fail( "a quoted comma-separated value must survive the scanner whole, got "
                        + ( mutations == null ? "null" : mutations.getValue() ) );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "quoted-comma end-to-end parse threw: " + e );
        }
    }

    /** Increment A -- an UNQUOTED value inside a kept bracket may legitimately contain spaces: real
     *  Auspice/Nextstrain Nexus exports write "country=Democratic Republic of the Congo" and
     *  "outbreak_geo=Kikwit 1995" with no quotes at all. Before this fix the scanner dropped every space
     *  unconditionally, silently producing a different (concatenated) value here than the same dataset's JSON
     *  download gives for the same field. */
    private static boolean testUnquotedSpaceSurvivesInBracket() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A[&country=Democratic Republic of the Congo,outbreak_geo=Kikwit 1995]:1.0,B:1.0);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            final Property country = prop( a, "beast:country" );
            final Property geo = prop( a, "beast:outbreak_geo" );
            if ( ( country == null ) || !"Democratic Republic of the Congo".equals( country.getValue() ) ) {
                return fail( "an unquoted space-containing value inside a bracket must survive, got "
                        + ( country == null ? "null" : country.getValue() ) );
            }
            if ( ( geo == null ) || !"Kikwit 1995".equals( geo.getValue() ) ) {
                return fail( "a second unquoted space-containing value on the same blob must also survive, got "
                        + ( geo == null ? "null" : geo.getValue() ) );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "unquoted-space end-to-end parse threw: " + e );
        }
    }

    /** Non-behaviour: the bracket exemption must NOT leak into ordinary label parsing. An unquoted space in a node
     *  label OUTSIDE any bracket is still silently dropped, exactly as before Increment A -- this is not valid
     *  Newick, but changing it was never part of the request and must not happen as a side effect. */
    private static boolean testUnquotedSpaceStillDroppedInLabel() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(New York:1.0,B:1.0);" );
            final Phylogeny phy = p.parse()[ 0 ];
            boolean found_new_york = false;
            for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
                final String name = it.next().getName();
                if ( "NewYork".equals( name ) ) {
                    found_new_york = true;
                }
                else if ( name.indexOf( ' ' ) >= 0 ) {
                    return fail( "a space in a label outside any bracket must be dropped, not kept, got '" + name
                            + "'" );
                }
            }
            if ( !found_new_york ) {
                return fail( "expected the space-less label 'NewYork' from unquoted \"New York\"" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "label-space end-to-end parse threw: " + e );
        }
    }

    /** A blob whose FIRST field is 'length' (TreeAnnotator emits field order varies) must still route to the
     *  structured parser -- guards the removed '[&length' routing exclusion in NHXParser. */
    private static boolean testLengthLedBlob() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "(A[&length=1.0,posterior=0.9,height=1.2,height_95%_HPD={0.8,1.5},rate=0.01]:1.0,B:1.0);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            if ( !a.getNodeData().isHasDate() || ( a.getNodeData().getDate().getMin() == null )
                    || !a.getBranchData().isHasConfidences() || ( prop( a, "beast:rate" ) == null ) ) {
                return fail( "a length-led BEAST blob must still be parsed (date interval + posterior + rate)" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "length-led parse threw: " + e );
        }
    }

    /** refKey sanitization ('%' -> '_'), the height_range interval fallback (and 95%-HPD precedence), and the
     *  'date' desc + plain 'height' value branches. */
    private static boolean testMoreFields() {
        final PhylogenyNode n = new PhylogenyNode();
        BeastAnnotationParser.apply( "&rate_95%_HPD={0.001,0.004}", n );
        final Property hpd = prop( n, "beast:rate_95_HPD" );
        if ( ( hpd == null ) || !"{0.001,0.004}".equals( hpd.getValue() ) || !"xsd:string".equals( hpd.getDataType() ) ) {
            return fail( "a 'rate_95%_HPD' key must sanitize to a valid beast:rate_95_HPD string property" );
        }
        final PhylogenyNode n2 = new PhylogenyNode();
        BeastAnnotationParser.apply( "&height=1.5,height_range={1.0,2.0}", n2 );
        final Date d2 = n2.getNodeData().getDate();
        if ( ( d2 == null ) || ( d2.getMin() == null ) || ( Math.abs( d2.getMin().doubleValue() - 1.0 ) > 1e-9 )
                || ( d2.getMax() == null ) || ( Math.abs( d2.getMax().doubleValue() - 2.0 ) > 1e-9 ) ) {
            return fail( "height_range must fill the date interval when no 95% HPD is present" );
        }
        final PhylogenyNode n3 = new PhylogenyNode();
        BeastAnnotationParser.apply( "&height_95%_HPD={0.8,1.2},height_range={0.5,1.5}", n3 );
        if ( Math.abs( n3.getNodeData().getDate().getMin().doubleValue() - 0.8 ) > 1e-9 ) {
            return fail( "height_95%_HPD must take precedence over height_range" );
        }
        final PhylogenyNode n4 = new PhylogenyNode();
        BeastAnnotationParser.apply( "&date=2014-06-10,height=3.0", n4 );
        final Date d4 = n4.getNodeData().getDate();
        if ( ( d4 == null ) || !"2014-06-10".equals( d4.getDesc() ) || ( d4.getValue() == null )
                || ( Math.abs( d4.getValue().doubleValue() - 3.0 ) > 1e-9 ) ) {
            return fail( "date -> desc and plain height -> value; got desc='" + ( d4 == null ? "?" : d4.getDesc() )
                    + "'" );
        }
        return true;
    }

    /** With the option OFF the raw [&...] blob is kept as an nh:comment (backward compatible), not structured. */
    private static boolean testOptionOffKeepsComment() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( false );
            p.setSource( "(A[&posterior=0.9,height=1.2]:1.0,B:1.0);" );
            final PhylogenyNode a = p.parse()[ 0 ].getNode( "A" );
            if ( a.getNodeData().isHasDate() || a.getBranchData().isHasConfidences() ) {
                return fail( "with the option OFF, BEAST fields must NOT be parsed into structured data" );
            }
            if ( prop( a, "nh:comment" ) == null ) {
                return fail( "with the option OFF, the raw [&...] blob must be kept as an nh:comment" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "option-off parse threw: " + e );
        }
    }

    /** A parsed BEAST tree, saved as phyloXML and reloaded, reproduces the date intervals, posteriors and rate
     *  properties (the <date>/<confidence>/<property> writers + parsers round-trip). */
    private static boolean testRoundTrip() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "((A[&height=0.0,rate=0.003]:1.0,B[&height=0.0,rate=0.002]:1.0)"
                    + "[&posterior=0.95,height=1.0,height_95%_HPD={0.8,1.3},rate=0.0025]:0.5,"
                    + "C[&height=0.0,rate=0.004]:1.5)[&posterior=1.0,height=1.5,height_95%_HPD={1.2,1.9},rate=0.003];" );
            final Phylogeny phy = p.parse()[ 0 ];
            final File tmp = File.createTempFile( "beast_roundtrip", ".xml" );
            try {
                new PhylogenyWriter().toPhyloXML( phy, 0, tmp );
                final Phylogeny back = ParserBasedPhylogenyFactory.getInstance()
                        .create( tmp, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
                int intervals = 0;
                int posteriors = 0;
                int rates = 0;
                for( final PhylogenyNodeIterator it = back.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if ( n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getMin() != null )
                            && ( n.getNodeData().getDate().getMax() != null ) ) {
                        intervals++;
                    }
                    if ( n.getBranchData().isHasConfidences() ) {
                        posteriors++;
                    }
                    if ( prop( n, "beast:rate" ) != null ) {
                        rates++;
                    }
                }
                if ( ( intervals != 2 ) || ( posteriors != 2 ) || ( rates != 5 ) ) {
                    return fail( "phyloXML round-trip lost data: intervals=" + intervals + " posteriors=" + posteriors
                            + " rates=" + rates );
                }
                return true;
            }
            finally {
                tmp.delete();
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "round-trip threw: " + e );
        }
    }

    private static boolean testFieldMapping() {
        final PhylogenyNode n = new PhylogenyNode();
        BeastAnnotationParser.apply(
                "&posterior=0.99,height_median=1.44,height_mean=1.40,height_95%_HPD={1.435,1.465},"
                        + "rate=0.0031,location=\"Africa\"",
                n );
        if ( !n.getBranchData().isHasConfidences()
                || ( Math.abs( n.getBranchData().getConfidence( 0 ).getValue() - 0.99 ) > 1e-9 )
                || !"posterior".equals( n.getBranchData().getConfidence( 0 ).getType() ) ) {
            return fail( "posterior must become a Confidence of type 'posterior'" );
        }
        final Date d = n.getNodeData().getDate();
        if ( d == null ) {
            return fail( "a height field must produce a <date>" );
        }
        if ( ( d.getValue() == null ) || ( Math.abs( d.getValue().doubleValue() - 1.44 ) > 1e-9 ) ) {
            return fail( "date value must be height_median (1.44), got " + d.getValue() );
        }
        if ( ( d.getMin() == null ) || ( Math.abs( d.getMin().doubleValue() - 1.435 ) > 1e-9 )
                || ( d.getMax() == null ) || ( Math.abs( d.getMax().doubleValue() - 1.465 ) > 1e-9 ) ) {
            return fail( "date min/max must be the 95% HPD bounds {1.435,1.465}, got " + d.getMin() + "/" + d.getMax() );
        }
        final Property rate = prop( n, "beast:rate" );
        if ( ( rate == null ) || !"xsd:decimal".equals( rate.getDataType() ) || !"0.0031".equals( rate.getValue() ) ) {
            return fail( "rate must become a numeric (xsd:decimal) beast:rate property" );
        }
        final Property loc = prop( n, "beast:location" );
        if ( ( loc == null ) || !"xsd:string".equals( loc.getDataType() ) || !"Africa".equals( loc.getValue() ) ) {
            return fail( "location must become a categorical property with quotes stripped, got "
                    + ( loc == null ? "null" : loc.getValue() ) );
        }
        // height_mean is used only when there is no height_median
        final PhylogenyNode n2 = new PhylogenyNode();
        BeastAnnotationParser.apply( "&height_mean=2.5", n2 );
        if ( ( n2.getNodeData().getDate() == null ) || ( n2.getNodeData().getDate().getValue() == null )
                || ( Math.abs( n2.getNodeData().getDate().getValue().doubleValue() - 2.5 ) > 1e-9 ) ) {
            return fail( "height_mean must be the date value when no median is present" );
        }
        return true;
    }

    private static boolean testHelpers() {
        final List<String> toks = BeastAnnotationParser.splitTopLevel( "a=1,b={2,3},c=\"x,y\",d=4" );
        if ( ( toks.size() != 4 ) || !toks.get( 1 ).equals( "b={2,3}" ) || !toks.get( 2 ).equals( "c=\"x,y\"" ) ) {
            return fail( "splitTopLevel must not split commas inside {} or \"\": " + toks );
        }
        final String[] iv = BeastAnnotationParser.parseInterval( "{1.2,3.4}" );
        if ( ( iv == null ) || !"1.2".equals( iv[ 0 ] ) || !"3.4".equals( iv[ 1 ] ) ) {
            return fail( "parseInterval must return the two raw numbers of {1.2,3.4}" );
        }
        if ( ( BeastAnnotationParser.parseInterval( "{1.2}" ) != null )
                || ( BeastAnnotationParser.parseInterval( "5" ) != null )
                || ( BeastAnnotationParser.parseInterval( "{a,b}" ) != null ) ) {
            return fail( "parseInterval must reject non-two-number sets" );
        }
        if ( ( BeastAnnotationParser.parseNumber( "0.003" ) == null )
                || ( BeastAnnotationParser.parseNumber( "abc" ) != null )
                || ( BeastAnnotationParser.parseNumber( "NaN" ) != null )
                || ( BeastAnnotationParser.parseNumber( "{1,2}" ) != null ) ) {
            return fail( "parseNumber must accept finite numbers only" );
        }
        return true;
    }

    private static boolean testMalformedTolerance() {
        final PhylogenyNode n = new PhylogenyNode();
        // a bare token, an empty value, and a good field mixed together -- the good fields must still parse
        BeastAnnotationParser.apply( "&garbage,empty=,posterior=0.9,rate=0.01", n );
        if ( !n.getBranchData().isHasConfidences()
                || ( Math.abs( n.getBranchData().getConfidence( 0 ).getValue() - 0.9 ) > 1e-9 ) ) {
            return fail( "malformed fields must be skipped while good fields (posterior) still parse" );
        }
        if ( prop( n, "beast:rate" ) == null ) {
            return fail( "rate must still parse alongside malformed fields" );
        }
        // null / empty input must be a safe no-op
        BeastAnnotationParser.apply( null, n );
        BeastAnnotationParser.apply( "&", new PhylogenyNode() );
        return true;
    }

    private static boolean testEndToEnd() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( true );
            p.setSource( "((A[&height=0.0,rate=0.003]:1.0,B[&height=0.0,rate=0.002]:1.0)"
                    + "[&posterior=0.95,height=1.0,height_95%_HPD={0.8,1.3},rate=0.0025]:0.5,"
                    + "C[&height=0.0,rate=0.004]:1.5)[&posterior=1.0,height=1.5,height_95%_HPD={1.2,1.9},rate=0.003];" );
            final Phylogeny[] phys = p.parse();
            if ( ( phys.length != 1 ) || ( phys[ 0 ].getNumberOfExternalNodes() != 3 ) ) {
                return fail( "should parse one 3-tip tree, got " + phys.length + " tree(s)" );
            }
            int internal_intervals = 0;
            int internal_posteriors = 0;
            int rate_props = 0;
            for( final PhylogenyNodeIterator it = phys[ 0 ].iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isExternal() && n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getMin() != null )
                        && ( n.getNodeData().getDate().getMax() != null ) ) {
                    internal_intervals++;
                }
                if ( !n.isExternal() && n.getBranchData().isHasConfidences() ) {
                    internal_posteriors++;
                }
                if ( prop( n, "beast:rate" ) != null ) {
                    rate_props++;
                }
            }
            if ( internal_intervals != 2 ) {
                return fail( "both internal nodes must carry an HPD interval, got " + internal_intervals );
            }
            if ( internal_posteriors != 2 ) {
                return fail( "both internal nodes must carry a posterior confidence, got " + internal_posteriors );
            }
            if ( rate_props != 5 ) {
                return fail( "all 5 nodes carry a beast:rate property, got " + rate_props );
            }
            // and NO opaque nh:comment blob is left behind
            for( final PhylogenyNodeIterator it = phys[ 0 ].iteratorPreorder(); it.hasNext(); ) {
                if ( prop( it.next(), "nh:comment" ) != null ) {
                    return fail( "the structured parse must leave no opaque nh:comment blob" );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "end-to-end parse threw: " + e );
        }
    }

    private static Property prop( final PhylogenyNode node, final String ref ) {
        if ( node.getNodeData().getProperties() == null ) {
            return null;
        }
        final List<Property> ps = node.getNodeData().getProperties().getProperties( ref );
        return ps.isEmpty() ? null : ps.get( 0 );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "BeastAnnotationParser test failed: " + msg );
        return false;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }
}
