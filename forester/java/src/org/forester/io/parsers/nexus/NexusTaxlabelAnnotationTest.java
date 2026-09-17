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

import java.awt.Color;
import java.io.ByteArrayInputStream;
import java.nio.charset.StandardCharsets;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Headless tests for the annotation a Nexus TAXLABELS entry may carry. FigTree writes a coloured TAXON there
 * ({@code 'NewYork_454_1999.05'[&!color=#-8381639]}), never in the tree string -- all 17 colours of a real FigTree
 * file (test_trees/influenza.tree) sit in that block. A coloured taxon is a coloured LABEL.
 */
public final class NexusTaxlabelAnnotationTest {

    private static final String HEAD = "#NEXUS\nbegin taxa;\n\tdimensions ntax=4;\n\ttaxlabels\n";
    private static final String TAIL = ";\nend;\n\nbegin trees;\n\ttree TREE1 = [&R] (('New_York':1.0,B:1.0)[&!color=#ff0000]:1.0,"
            + "('Seba''s bat':1.0,D:1.0):1.0);\nend;\n";
    private static final String FIGTREE = HEAD + "\t'New_York'[&!color=#-8381639]\n\tB\n\t'Seba''s bat'[&!color=#0000ff]\n"
            + "\tD[a plain comment]\n" + TAIL;

    public static boolean test() {
        return testSplit() && testLabelColours() && testOptionOffReadsNoColour() && testIndexedTipsGetCleanNames()
                && testOneNamespacePerTree() && testJoinKeyNeverCrossesTaxa() && testAnnotationWithSpaces();
    }

    private static boolean testSplit() {
        final String[][] cases = { { "'New_York'[&!color=#-8381639]", "'New_York'", "&!color=#-8381639" },
                { "B[&!color=#ff0000,!name=\"x\"]", "B", "&!color=#ff0000,!name=\"x\"" }, { "B", "B", null },
                { "D[a plain comment]", "D", null }, { "[&!color=#ff0000]", "", "&!color=#ff0000" },
                { "'a[b]c'", "'a[b]c'", null }, { "'it''s[x]'[&k=v]", "'it''s[x]'", "&k=v" },
                { "B[&unclosed", "B", "&unclosed" } };
        for( final String[] c : cases ) {
            final String[] r = NexusPhylogeniesParser.splitTaxlabelAnnotation( c[ 0 ] );
            if ( !c[ 1 ].equals( r[ 0 ] ) || ( ( c[ 2 ] == null ) ? ( r[ 1 ] != null ) : !c[ 2 ].equals( r[ 1 ] ) ) ) {
                return fail( c[ 0 ] + " must split into label '" + c[ 1 ] + "' + annotation '" + c[ 2 ] + "', got '"
                        + r[ 0 ] + "' + '" + r[ 1 ] + "'" );
            }
        }
        return true;
    }

    /** The taxon's colour is the tip's LABEL colour; the tree string's "!color" stays the BRANCH colour. */
    private static boolean testLabelColours() {
        try {
            final Phylogeny phy = parse( FIGTREE, true );
            final PhylogenyNode ny = phy.getNode( "New_York" );
            if ( ( fontColor( ny ) == null ) || ( fontColor( ny ).getRGB() != 0xFF801B39 ) ) {
                return fail( "'New_York'[&!color=#-8381639] must give the tip the label colour 0x801B39, got "
                        + fontColor( ny ) );
            }
            if ( ny.getBranchData().getBranchColor() != null ) {
                return fail( "a coloured TAXON is a coloured label: it must not colour the tip's branch" );
            }
            final PhylogenyNode seba = phy.getNode( "Seba's bat" );
            if ( ( seba == null ) || ( fontColor( seba ) == null ) || ( fontColor( seba ).getRGB() != 0xFF0000FF ) ) {
                return fail( "a quoted label with an escaped quote must still find its annotation" );
            }
            if ( ( fontColor( phy.getNode( "B" ) ) != null ) || ( fontColor( phy.getNode( "D" ) ) != null ) ) {
                return fail( "a taxon without a colour gets none" );
            }
            final PhylogenyNode clade = ny.getParent();
            if ( ( clade.getBranchData().getBranchColor() == null )
                    || ( clade.getBranchData().getBranchColor().getValue().getRGB() != 0xFFFF0000 )
                    || ( fontColor( clade ) != null ) ) {
                return fail( "a !color in the TREE STRING stays the branch colour, never a label colour" );
            }
            // the tree may spell a taxon with a space where TAXLABELS has an underscore (Nexus treats them alike)
            final Phylogeny spaced = parse( FIGTREE.replace( "('New_York':1.0", "('New York':1.0" ), true );
            if ( fontColor( spaced.getNode( "New York" ) ) == null ) {
                return fail( "'New York' in the tree must find the annotation of the taxlabel 'New_York'" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "taxlabel-colour parse threw: " + e );
        }
    }

    private static boolean testOptionOffReadsNoColour() {
        try {
            final Phylogeny phy = parse( FIGTREE, false );
            if ( fontColor( phy.getNode( "New_York" ) ) != null ) {
                return fail( "with the bracket-annotation option OFF no taxlabel annotation is applied" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "option-off parse threw: " + e );
        }
    }

    /** When the tips are TAXLABELS indices they are NAMED from the block -- and used to be named
     *  "New_York[&!color=#-8381639]" / "D[a plain comment]", the annotation leaking into the name. */
    private static boolean testIndexedTipsGetCleanNames() {
        try {
            final String indexed = FIGTREE.substring( 0, FIGTREE.indexOf( "begin trees;" ) )
                    + "begin trees;\n\ttree TREE1 = [&R] ((1:1.0,2:1.0):1.0,(3:1.0,4:1.0):1.0);\nend;\n";
            for( final boolean option_on : new boolean[] { true, false } ) {
                final Phylogeny phy = parse( indexed, option_on );
                for( final String name : new String[] { "New_York", "B", "Seba's bat", "D" } ) {
                    if ( phy.getNode( name ) == null ) {
                        return fail( "an indexed tip must be named by its CLEAN taxlabel '" + name + "'" );
                    }
                }
                if ( option_on && ( fontColor( phy.getNode( "New_York" ) ) == null ) ) {
                    return fail( "an indexed tip gets its taxlabel's colour too" );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "indexed-tips parse threw: " + e );
        }
    }

    /** The per-tree pass sees the TAXLABELS annotations too: on a TreeTime tree (mutations, no node ages) EVERY
     *  bracket-annotation ref is treetime:, whichever block it came from -- it used to run before the taxlabels were
     *  applied, leaving one tree with both beast: and treetime: refs. And the pass runs exactly once. */
    private static boolean testOneNamespacePerTree() {
        try {
            final String treetime = "#NEXUS\nbegin taxa;\n\tdimensions ntax=2;\n\ttaxlabels\n\tA[&!color=-8381639]\n\tB[&host=duck]\n;\nend;\n"
                    + "begin trees;\n\ttree t = [&R] (A:1[&mutations=\"A1G\"],B:1[&mutations=\"C2T\"]);\nend;\n";
            final Phylogeny phy = parse( treetime, true );
            for( final String name : new String[] { "A", "B" } ) {
                for( final org.forester.phylogeny.data.Property p : phy.getNode( name ).getNodeData().getProperties().getProperties() ) {
                    if ( p.getRef().startsWith( "beast:" ) ) {
                        return fail( "a TreeTime tree must carry treetime: refs only, tip " + name + " has " + p.getRef() );
                    }
                }
            }
            if ( phy.getNode( "A" ).getNodeData().getProperties().getProperties( "treetime:_color" ).isEmpty()
                    || phy.getNode( "B" ).getNodeData().getProperties().getProperties( "treetime:host" ).isEmpty()
                    || phy.getNode( "A" ).getNodeData().getProperties().getProperties( "treetime:mutations" ).isEmpty() ) {
                return fail( "the taxlabel annotations must be renamed with the rest (treetime:_color, treetime:host)" );
            }
            // a TIME-scaled TreeTime tree: the dates are promoted ONCE (a second pass would see date values, call the
            // tree "not TreeTime's own" and skip the rename -- the ordering this test exists for)
            final String timetree = "#NEXUS\nbegin taxa;\n\tdimensions ntax=3;\n\ttaxlabels\n\tA[&!color=-8381639]\n\tB\n\tC\n;\nend;\n"
                    + "begin trees;\n\ttree t = [&R] ((A:1.0[&mutations=\"A1G\",date=2003.00],B:2.5[&date=2004.50]):1.0[&date=2002.00],"
                    + "C:2.25[&date=2003.25]):0.0[&date=2001.00];\nend;\n";
            final Phylogeny tt = parse( timetree, true );
            if ( ( tt.getNode( "A" ).getNodeData().getDate().getValue() == null )
                    || tt.getNode( "A" ).getNodeData().getProperties().getProperties( "treetime:_color" ).isEmpty() ) {
                return fail( "a time-scaled TreeTime tree: dates promoted AND the taxlabel colour renamed" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "one-namespace parse threw: " + e );
        }
    }

    /** "Taxon_A" and "taxon_a" are two taxa: one's annotation must never reach the other through the loose join key.
     *  The key serves only a tip that is not itself a taxlabel, and only when it names exactly one taxon. */
    private static boolean testJoinKeyNeverCrossesTaxa() {
        try {
            final String head = "#NEXUS\nbegin taxa;\n\tdimensions ntax=2;\n\ttaxlabels\n";
            // one annotated, one not
            Phylogeny phy = parse( head + "\tTaxon_A[&!color=#-65536]\n\ttaxon_a\n;\nend;\nbegin trees;\n\ttree t = [&R] (Taxon_A:1,taxon_a:1);\nend;\n",
                                   true );
            if ( ( fontColor( phy.getNode( "Taxon_A" ) ) == null ) || ( fontColor( phy.getNode( "taxon_a" ) ) != null ) ) {
                return fail( "only Taxon_A is coloured; taxon_a, a taxon of its own, must not borrow its colour" );
            }
            // both annotated: each keeps its own, whatever order they came in
            phy = parse( head + "\tTaxon_A[&!color=#ff0000]\n\ttaxon_a[&!color=#0000ff]\n;\nend;\nbegin trees;\n\ttree t = [&R] (Taxon_A:1,taxon_a:1);\nend;\n",
                         true );
            if ( ( fontColor( phy.getNode( "Taxon_A" ) ).getRGB() != 0xFFFF0000 )
                    || ( fontColor( phy.getNode( "taxon_a" ) ).getRGB() != 0xFF0000FF ) ) {
                return fail( "two annotated taxa that share a join key each keep their own colour" );
            }
            // a tip spelled neither way, whose key names TWO taxa: ambiguous, no colour
            phy = parse( "#NEXUS\nbegin taxa;\n\tdimensions ntax=3;\n\ttaxlabels\n\tTaxon_A[&!color=#ff0000]\n\ttaxon_a[&!color=#0000ff]\n\tB\n;\nend;\n"
                    + "begin trees;\n\ttree t = [&R] ('Taxon a':1,B:1);\nend;\n", true );
            if ( fontColor( phy.getNode( "Taxon a" ) ) != null ) {
                return fail( "a tip whose join key names two taxa is ambiguous and gets no colour" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "join-key parse threw: " + e );
        }
    }

    /** An annotation with a space inside stays with its taxon: split apart, "y\"]" would become a taxon of its own and
     *  every later taxlabel index -- how a tree with numbered tips is named -- would be off by one. */
    private static boolean testAnnotationWithSpaces() {
        try {
            final String nexus = "#NEXUS\nbegin taxa;\n\tdimensions ntax=3;\n\ttaxlabels\n\t'A'[&!color=#ff0000,note=\"x y  z\"] 'B[1'\n\tC\n;\nend;\n"
                    + "begin trees;\n\ttree t = [&R] ((1:1,2:1):1,3:1);\nend;\n";
            final Phylogeny phy = parse( nexus, true );
            if ( ( phy.getNodes( "A" ).size() != 1 ) || ( phy.getNodes( "B[1" ).size() != 1 ) || ( phy.getNodes( "C" ).size() != 1 ) ) {
                return fail( "numbered tips must be named A, B[1, C -- the annotation's inner spaces must not make extra taxa" );
            }
            final PhylogenyNode a = phy.getNode( "A" );
            if ( ( fontColor( a ) == null ) || a.getNodeData().getProperties().getProperties( "beast:note" ).isEmpty() ) {
                return fail( "A keeps its whole annotation: colour AND note" );
            }
            final java.util.List<String> joined = NexusPhylogeniesParser.joinBracketedTokens( java.util.Arrays.asList( "'a[b'", "c[&k=\"1", "2\"]", "d" ) );
            if ( !java.util.Arrays.asList( "'a[b'", "c[&k=\"1 2\"]", "d" ).equals( joined ) ) {
                return fail( "joinBracketedTokens: a bracket inside a quoted label does not count; a split annotation is rejoined, got " + joined );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "annotation-with-spaces parse threw: " + e );
        }
    }

    private static Phylogeny parse( final String nexus, final boolean option_on ) throws Exception {
        final NexusPhylogeniesParser p = new NexusPhylogeniesParser();
        p.setParseBeastStyleExtendedTags( option_on );
        // a String source is a file NAME to this parser: hand it the text as a stream
        return ParserBasedPhylogenyFactory.getInstance()
                .create( new ByteArrayInputStream( nexus.getBytes( StandardCharsets.UTF_8 ) ), p )[ 0 ];
    }

    private static Color fontColor( final PhylogenyNode n ) {
        final NodeVisualData v = n.getNodeData().getNodeVisualData();
        return ( v == null ) ? null : v.getFontColor();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "NexusTaxlabelAnnotation test failed: " + msg );
        return false;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }
}
