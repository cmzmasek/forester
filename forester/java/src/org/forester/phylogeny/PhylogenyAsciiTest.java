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

// Phylogeny.toAscii(): the tree drawn as plain ASCII for a terminal, and the size it refuses.
//
// Two things this pins beyond the drawing itself. A node's name is the one toNewHampshire would write
// (newHampshireLabel), NOT node.getName() -- much of real phyloXML carries its tips' identity in <taxonomy>
// or <sequence>, and under getName() such a tree would draw as a column of blank lines, which is exactly the
// tree you reached for toAscii to look at. And the drawing carries structure and names only: a branch length,
// a support value or a molecular sequence appearing here would be a defect, so each is built into a fixture
// and asserted ABSENT rather than assumed absent.

package org.forester.phylogeny;

import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.FailedConditionCheckException;
import org.forester.util.ForesterUtil;

public class PhylogenyAsciiTest {

    private static final String NL = ForesterUtil.LINE_SEPARATOR;

    private static boolean fail( final String msg ) {
        System.out.println( "  [PhylogenyAsciiTest] " + msg );
        return false;
    }

    /** The drawing with every line separator normalised to \n, so an expectation reads as one block here. */
    private static String drawn( final Phylogeny p ) {
        return p.toAscii().replace( NL, "\n" );
    }

    public static boolean test() {
        try {
            // (1) the shape, in full. Covers what is easiest to get subtly wrong: the last child takes a
            // backtick, every earlier one a plus, and an ancestor that still has siblings below it keeps a
            // vertical bar running down through its subtree while the last one does not.
            // Note what is NOT here: a line of its own for the unnamed root. A node earns a line only by
            // being a tip or carrying a label, so the root contributes its elbow and Primates continues on
            // the same line. Drawn the other way this six-line tree was eleven lines, five of them bare
            // connectors -- which is what Christian objected to, and rightly.
            final String expected = "+--+-- Primates\n"
                    + "   |  +-- Human\n"
                    + "   |  `-- Chimp\n"
                    + "   `-- Mouse\n";
            final String got = drawn( Phylogeny
                    .createInstanceFromNhxString( "((Human,Chimp)Primates,Mouse)" ) );
            if ( !expected.equals( got ) ) {
                return fail( "the drawing changed:\n" + got + "--- expected ---\n" + expected );
            }
            // (2) a deeper nesting: the bar must stop at the last child at EVERY level, not just the first
            final String deep = drawn( Phylogeny.createInstanceFromNhxString( "(A,(B,(C,D)))" ) );
            if ( !deep.equals( "+--+-- A\n"
                    + "   `--+-- B\n"
                    + "      `--+-- C\n"
                    + "         `-- D\n" ) ) {
                return fail( "nested last-children draw wrong:\n" + deep );
            }
            // (2b) THE LINE RULE, stated as a count. A tree whose internal nodes are all unnamed -- which is
            // what a gene tree looks like -- must draw exactly one line per tip: nothing is spent saying
            // "there is a branch point here", because the elbows already say it.
            final String genes = drawn( Phylogeny.createInstanceFromNhxString(
                    "(((BCL2_human,BCL2_chimp),BCL2_mouse),((BCLX_human,BCLX_chimp),BCLX_mouse))" ) );
            final int gene_lines = genes.split( "\n", -1 ).length - 1; // the trailing separator is not a line
            if ( gene_lines != 6 ) {
                return fail( "a 6-tip tree with no internal names must draw 6 lines, drew " + gene_lines
                        + ":\n" + genes );
            }
            // the NEIGHBOURING case, differing by exactly the thing under test: name ONE internal node and
            // the drawing grows by exactly that one line. Without this, "never draw internal nodes at all"
            // would pass the count above just as well -- and would throw away the clade names that make the
            // bat and dinosaur trees readable.
            final Phylogeny one_named = Phylogeny.createInstanceFromNhxString(
                    "(((BCL2_human,BCL2_chimp),BCL2_mouse),((BCLX_human,BCLX_chimp),BCLX_mouse))" );
            one_named.getNode( "BCL2_human" ).getParent().getParent().setName( "BCL2_clade" );
            final String named_drawn = drawn( one_named );
            final int named_lines = named_drawn.split( "\n", -1 ).length - 1;
            if ( named_lines != 7 ) {
                return fail( "naming one internal node must add exactly one line (6 -> 7), got " + named_lines
                        + ":\n" + named_drawn );
            }
            if ( !named_drawn.contains( "BCL2_clade" ) ) {
                return fail( "a NAMED internal node must still draw its own line:\n" + named_drawn );
            }
            // (3) a multifurcation: every child but the last takes a plus
            final String multi = drawn( Phylogeny.createInstanceFromNhxString( "(A,B,C,D)" ) );
            if ( !multi.equals( "+--+-- A\n   +-- B\n   +-- C\n   `-- D\n" ) ) {
                return fail( "a multifurcation draws wrong:\n" + multi );
            }
            // (4) THE LABEL RULE: the name toNewHampshire would write, not node.getName().
            final Phylogeny named = Phylogeny.createInstanceFromNhxString( "(a,b,c,d,e,f)" );
            // By NAME, never by index: getExternalNodes()' javadoc says the order "is random -- and hence
            // cannot be relied on". Indexing worked only because today's implementation happens to walk in
            // tree order, and the bare-connector assertion below depended on the nameless tip landing LAST.
            final PhylogenyNode by_tax = named.getNode( "a" );
            by_tax.setName( "" );
            final Taxonomy tax = new Taxonomy();
            tax.setScientificName( "Nematostella vectensis" );
            by_tax.getNodeData().setTaxonomy( tax );
            final PhylogenyNode by_seq = named.getNode( "b" );
            by_seq.setName( "" );
            final Sequence seq = new Sequence();
            seq.setName( "Apoptosis regulator Bcl-2" );
            by_seq.getNodeData().addSequence( seq );
            final PhylogenyNode by_acc = named.getNode( "c" );
            by_acc.setName( "" );
            final Sequence acc_seq = new Sequence();
            acc_seq.setAccession( new Accession( "P10415", "uniprot" ) );
            by_acc.getNodeData().addSequence( acc_seq );
            // the NEIGHBOURING case, differing by exactly the thing under test: this one HAS a name, and a
            // taxonomy too, so it must print the NAME -- otherwise the test would pass for a chain that
            // simply always preferred the taxonomy.
            final PhylogenyNode both = named.getNode( "d" );
            both.setName( "MyOwnName" );
            final Taxonomy other = new Taxonomy();
            other.setScientificName( "Homo sapiens" );
            both.getNodeData().setTaxonomy( other );
            // a tip with BOTH a name and an accession: the name wins. Without this the fixture cannot tell
            // newHampshireLabel(false,..) from newHampshireLabel(true,..) -- every other node here has only
            // one of the two, so forcing sequence ids would change nothing and the mutation survived.
            final PhylogenyNode name_and_acc = named.getNode( "e" );
            name_and_acc.setName( "HasBoth" );
            final Sequence both_seq = new Sequence();
            both_seq.setAccession( new Accession( "Q99999", "uniprot" ) );
            name_and_acc.getNodeData().addSequence( both_seq );
            // a tip that NOTHING names: it must draw as a bare connector, never as an invented placeholder.
            // Every other nameless node in this fixture is internal, and the placeholder only ever applies to
            // an external one -- so without this tip, inventing placeholders changed nothing either.
            named.getNode( "f" ).setName( "" );
            final String labels = drawn( named );
            for( final String want : new String[] { "Nematostella vectensis", "Apoptosis regulator Bcl-2",
                    "P10415", "MyOwnName" } ) {
                if ( !labels.contains( want ) ) {
                    return fail( "the label chain must reach [" + want + "]:\n" + labels );
                }
            }
            if ( labels.contains( "Homo sapiens" ) ) {
                return fail( "a node with a name must draw the NAME, not its taxonomy:\n" + labels );
            }
            if ( labels.contains( "Q99999" ) ) {
                return fail( "a node with a name must draw the NAME, not its accession:\n" + labels );
            }
            boolean bare = false;
            for( final String line : labels.split( "\n", -1 ) ) {
                // an elbow and nothing after it, wherever in the drawing it fell
                if ( line.endsWith( "+--" ) || line.endsWith( "`--" ) ) {
                    bare = true;
                }
            }
            if ( !bare ) {
                return fail( "a tip nothing names must draw as a bare connector:\n" + labels );
            }
            // the PLACEHOLDER's own shape, not the word "node": a label like "nodulin-26" or a taxonomy
            // "Nodosaurus" would trip a bare contains("node") against entirely correct code
            if ( java.util.regex.Pattern.compile( "--\\s+node\\d+" ).matcher( labels ).find() ) {
                return fail( "no placeholder may be invented for a nameless tip:\n" + labels );
            }
            // (5) DELIBERATE NON-BEHAVIOUR: structure and names only. Each of these is put ON the tree and
            // then required to be absent -- assuming absence would pass on a tree that never had them.
            // The support value goes on as a CONFIDENCE, not as a Newick internal label: an internal label
            // is a NAME until AptxUtil.applyInternalLabelPolicy turns it into one, so writing "0.97" in the
            // string would have put the node's name at 0.97 and the assertion below would have failed
            // against correct code. It did, on this test's first run.
            final Phylogeny rich = Phylogeny
                    .createInstanceFromNhxString( "((Human:0.12345,Chimp:0.6789):0.5,Mouse:0.4)" );
            rich.getNode( "Human" ).getParent().getBranchData()
                    .addConfidence( new Confidence( 97.0, "bootstrap" ) );
            final PhylogenyNode h = rich.getNode( "Human" );
            final Sequence mol = new Sequence();
            mol.setName( "Human" );
            mol.setMolecularSequence( "MKALIVWQNP" );
            h.getNodeData().addSequence( mol );
            final String plain = drawn( rich );
            for( final String unwanted : new String[] { "0.12345", "0.6789", "97", "0.4", "MKALIVWQNP" } ) {
                if ( plain.contains( unwanted ) ) {
                    return fail( "[" + unwanted + "] must not appear in a structure-and-names drawing:\n"
                            + plain );
                }
            }
            // (6) no line may carry trailing whitespace, and the drawing ends in a newline
            if ( !plain.endsWith( "\n" ) ) {
                return fail( "the drawing must end in a line separator" );
            }
            for( final String line : plain.split( "\n", -1 ) ) {
                if ( ( line.length() > 0 ) && Character.isWhitespace( line.charAt( line.length() - 1 ) ) ) {
                    return fail( "a line carries trailing whitespace: [" + line + "]" );
                }
            }
            // (7) drawing a tree must not change it. On a FRESHLY parsed tree, because `rich` has already
            // been drawn above -- a first-call side effect (writing back a fallback name, say) would already
            // have been applied, and the baseline would record the mutated tree and agree with itself.
            final Phylogeny pure = Phylogeny.createInstanceFromNhxString( "((Human,Chimp)Primates,Mouse)" );
            final String before = pure.toNewHampshire();
            pure.toAscii();
            if ( !before.equals( pure.toNewHampshire() ) ) {
                return fail( "toAscii changed the tree it drew" );
            }
            // (8) edges: an empty tree, and a single node
            if ( !new Phylogeny().toAscii().isEmpty() ) {
                return fail( "an empty tree draws nothing" );
            }
            if ( !drawn( Phylogeny.createInstanceFromNhxString( "A" ) ).equals( "+-- A\n" ) ) {
                return fail( "a one-node tree draws wrong: "
                        + drawn( Phylogeny.createInstanceFromNhxString( "A" ) ) );
            }
            // (9) the size limit, on both sides of it
            if ( Phylogeny.MAX_ASCII_EXTERNAL_NODES != 1000 ) {
                return fail( "the documented limit is 1000, not " + Phylogeny.MAX_ASCII_EXTERNAL_NODES );
            }
            final Phylogeny at_limit = star( Phylogeny.MAX_ASCII_EXTERNAL_NODES );
            if ( at_limit.getNumberOfExternalNodes() != Phylogeny.MAX_ASCII_EXTERNAL_NODES ) {
                return fail( "the fixture is not the size it claims ("
                        + at_limit.getNumberOfExternalNodes() + ")" );
            }
            final String big = at_limit.toAscii(); // exactly at the limit: drawn, not refused
            if ( !big.contains( "tip999" ) ) {
                return fail( "a tree AT the limit must be drawn in full" );
            }
            final Phylogeny over = star( Phylogeny.MAX_ASCII_EXTERNAL_NODES + 1 );
            try {
                over.toAscii();
                return fail( "a tree over the limit must be refused, not drawn" );
            }
            catch ( final FailedConditionCheckException e ) {
                final String m = String.valueOf( e.getMessage() );
                // the message has to say what was too big and what the limit is, or the caller cannot act
                if ( !m.contains( "1001" ) || !m.contains( "1000" ) ) {
                    return fail( "the refusal must name the count and the limit: " + m );
                }
            }
            // (10) the guard must not be fooled by the EXTERNAL-NODE CACHE. getNumberOfExternalNodes() is
            // served from a cache that setRoot() does not invalidate, so a small tree re-rooted onto a large
            // one reports the old count; measured, a two-tip tree given a 5000-tip root drew all 5000 lines.
            final Phylogeny stale = Phylogeny.createInstanceFromNhxString( "(A,B)" );
            stale.getNumberOfExternalNodes(); // warm the cache at 2
            stale.setRoot( star( Phylogeny.MAX_ASCII_EXTERNAL_NODES + 1 ).getRoot() );
            try {
                stale.toAscii();
                return fail( "the size guard was fooled by a stale external-node count ("
                        + stale.getNumberOfExternalNodes() + " reported)" );
            }
            catch ( final FailedConditionCheckException refused_stale ) {
                // refused, as it must be
            }
            // (11) DEPTH, which the tip count cannot cover: a chain of unnamed nodes ending in ONE tip passes
            // any tip limit and used to die with a StackOverflowError rather than the documented refusal.
            final PhylogenyNode chain_root = new PhylogenyNode();
            PhylogenyNode cur = chain_root;
            for( int i = 0; i < ( Phylogeny.MAX_ASCII_EXTERNAL_NODES + 50 ); ++i ) {
                final PhylogenyNode next = new PhylogenyNode();
                cur.addAsChild( next );
                cur = next;
            }
            cur.setName( "the_one_tip" );
            final Phylogeny deep_chain = new Phylogeny();
            deep_chain.setRoot( chain_root );
            deep_chain.setRooted( true );
            if ( deep_chain.getNumberOfExternalNodes() != 1 ) {
                return fail( "the deep fixture must have exactly one tip, so only DEPTH can refuse it" );
            }
            try {
                deep_chain.toAscii();
                return fail( "a tree deeper than the limit must be refused" );
            }
            catch ( final FailedConditionCheckException refused_deep ) {
                if ( !String.valueOf( refused_deep.getMessage() ).contains( "deep" ) ) {
                    return fail( "the refusal must say it was the DEPTH: " + refused_deep.getMessage() );
                }
            }
            catch ( final StackOverflowError e ) {
                return fail( "a deep tree must be refused, not overflow the stack" );
            }
            // a tree just inside the depth limit still draws
            final PhylogenyNode ok_root = new PhylogenyNode();
            PhylogenyNode c2 = ok_root;
            for( int i = 0; i < ( Phylogeny.MAX_ASCII_EXTERNAL_NODES - 10 ); ++i ) {
                final PhylogenyNode next = new PhylogenyNode();
                c2.addAsChild( next );
                c2 = next;
            }
            c2.setName( "deep_but_legal" );
            final Phylogeny ok_deep = new Phylogeny();
            ok_deep.setRoot( ok_root );
            ok_deep.setRooted( true );
            if ( !ok_deep.toAscii().contains( "deep_but_legal" ) ) {
                return fail( "a tree inside the depth limit must still be drawn" );
            }
            // (12) a name is free text: a line break in one must not split a node across two lines, and a
            // trailing space must not survive. Newick answers this by quoting; a drawing cannot.
            final Phylogeny dirty = Phylogeny.createInstanceFromNhxString( "(x,y)" );
            dirty.getNode( "x" ).setName( "has\nnewline" );
            dirty.getNode( "y" ).setName( "trailing " );
            final String cleaned = drawn( dirty );
            final String[] lines = cleaned.split( "\n", -1 );
            if ( ( lines.length - 1 ) != 2 ) {
                return fail( "a newline in a name split the drawing into " + ( lines.length - 1 )
                        + " lines:\n" + cleaned );
            }
            for( final String line : lines ) {
                if ( ( line.length() > 0 ) && Character.isWhitespace( line.charAt( line.length() - 1 ) ) ) {
                    return fail( "a name's trailing space survived into the drawing: [" + line + "]" );
                }
                if ( ( line.length() > 0 ) && !line.startsWith( "+" ) && !line.startsWith( " " ) ) {
                    return fail( "a line of the drawing does not start with the tree: [" + line + "]" );
                }
            }
            if ( !cleaned.contains( "has newline" ) ) {
                return fail( "a control character in a name should become a space:\n" + cleaned );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
    }

    /** A root with {@code tips} children -- {@code tips} external nodes, cheaply. */
    private static Phylogeny star( final int tips ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < tips; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "tip" + i );
            root.addAsChild( tip );
        }
        final Phylogeny p = new Phylogeny();
        p.setRoot( root );
        p.setRooted( true );
        return p;
    }

    public static void main( final String[] args ) {
        if ( test() ) {
            System.out.println( "PhylogenyAsciiTest: OK." );
        }
        else {
            System.out.println( "PhylogenyAsciiTest: FAILED." );
        }
    }
}
