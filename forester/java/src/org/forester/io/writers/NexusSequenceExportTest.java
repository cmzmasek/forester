// Tests for writing molecular sequences as a Nexus Characters block.
//
// See PhylogenyWriter.writeNexusCharactersBlock.

package org.forester.io.writers;

import java.io.File;

import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

public class NexusSequenceExportTest {

    private static Phylogeny tree( final String nh ) throws Exception {
        return Phylogeny.createInstanceFromNhxString( nh );
    }

    private static void setSeq( final Phylogeny phy, final int i, final String mol ) {
        final PhylogenyNode n = phy.getExternalNodes().get( i );
        final Sequence s = new Sequence();
        s.setMolecularSequence( mol );
        n.getNodeData().setSequence( s );
    }

    private static String nexus( final Phylogeny phy ) throws Exception {
        return new PhylogenyWriter().toNexus( phy, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ).toString();
    }

    private static String charactersBlock( final String nexus ) {
        final int start = nexus.indexOf( "Begin Characters;" );
        if ( start < 0 ) {
            return null;
        }
        return nexus.substring( start, nexus.indexOf( "Begin Trees;" ) );
    }

    private static String taxaBlock( final String nexus ) {
        return nexus.substring( nexus.indexOf( "Begin Taxa;" ), nexus.indexOf( "Begin Characters;" ) < 0
                ? nexus.indexOf( "Begin Trees;" ) : nexus.indexOf( "Begin Characters;" ) );
    }

    public static boolean test() {
        try {
            // (1) a tree whose tips carry an alignment gets a Characters block
            final Phylogeny p1 = tree( "((A,B),C)" );
            setSeq( p1, 0, "MKAL-IV" );
            setSeq( p1, 1, "MKAL-IW" );
            setSeq( p1, 2, "MKAL-IY" );
            final String n1 = nexus( p1 );
            final String b1 = charactersBlock( n1 );
            if ( b1 == null ) {
                System.out.println( "no Characters block for a tree with an alignment" );
                return false;
            }
            if ( !b1.contains( "NChar=7" ) ) {
                System.out.println( "wrong NChar, expected 7: " + b1 );
                return false;
            }
            if ( !b1.contains( "DataType=Protein" ) ) {
                System.out.println( "protein sequences not typed as Protein: " + b1 );
                return false;
            }
            if ( !b1.contains( "MKAL-IV" ) || !b1.contains( "MKAL-IW" ) || !b1.contains( "MKAL-IY" ) ) {
                System.out.println( "a sequence is missing from the matrix: " + b1 );
                return false;
            }
            // NTax in a Characters block's Dimensions is illegal without NEWTAXA. The taxon count comes from the
            // Taxa block, which is why only NChar may appear here.
            if ( b1.contains( "NTax" ) ) {
                System.out.println( "NTax must not appear in a Characters block: " + b1 );
                return false;
            }
            // block order must be Taxa, Characters, Trees
            if ( !( n1.indexOf( "Begin Taxa;" ) < n1.indexOf( "Begin Characters;" )
                    && n1.indexOf( "Begin Characters;" ) < n1.indexOf( "Begin Trees;" ) ) ) {
                System.out.println( "Nexus blocks out of order" );
                return false;
            }
            // (2) deliberate non-behaviour: no sequences, no block at all
            final String n2 = nexus( tree( "((A,B),C)" ) );
            if ( charactersBlock( n2 ) != null ) {
                System.out.println( "a tree with no sequences must not get a Characters block" );
                return false;
            }
            // (3) deliberate non-behaviour: unaligned sequences are not a character matrix. Padding them would
            // state an alignment that does not exist, so nothing is written -- but a comment says why.
            final Phylogeny p3 = tree( "((A,B),C)" );
            setSeq( p3, 0, "MKAL" );
            setSeq( p3, 1, "MKALIVGD" );
            final String n3 = nexus( p3 );
            if ( charactersBlock( n3 ) != null ) {
                System.out.println( "unequal-length sequences must not be written as a matrix" );
                return false;
            }
            if ( !n3.contains( "unequal length" ) ) {
                System.out.println( "no comment explaining why the sequences were dropped" );
                return false;
            }
            // the comment must be a legal Nexus comment, and the file must still parse
            if ( !n3.contains( "[ Molecular sequences were not written" ) || !n3.contains( "]" ) ) {
                System.out.println( "the explanation is not a bracketed Nexus comment" );
                return false;
            }
            // (4) partial coverage: tips without a sequence become rows of the missing symbol, so the matrix
            // still covers every taxon in the Taxa block
            final Phylogeny p4 = tree( "((A,B),C)" );
            setSeq( p4, 0, "MKAL" );
            setSeq( p4, 2, "MKIV" );
            final String b4 = charactersBlock( nexus( p4 ) );
            if ( b4 == null ) {
                System.out.println( "partial sequence coverage should still produce a matrix" );
                return false;
            }
            if ( !b4.contains( "????" ) ) {
                System.out.println( "a tip with no sequence must get a missing-data row: " + b4 );
                return false;
            }
            if ( !b4.contains( "Missing=?" ) ) {
                System.out.println( "the missing symbol must be declared in Format" );
                return false;
            }
            int rows = 0;
            for( final String line : b4.split( "\\R" ) ) {
                final String t = line.trim();
                if ( t.startsWith( "A " ) || t.startsWith( "B " ) || t.startsWith( "C " ) ) {
                    ++rows;
                }
            }
            if ( rows != 3 ) {
                System.out.println( "expected one matrix row per taxon, got " + rows );
                return false;
            }
            // (5) nucleotide sequences are typed as DNA, not silently called Protein
            final Phylogeny p5 = tree( "(A,B)" );
            setSeq( p5, 0, "ACGTACGTACGTACGTACGT" );
            setSeq( p5, 1, "ACGTACGTACGTACGTACGA" );
            if ( !charactersBlock( nexus( p5 ) ).contains( "DataType=DNA" ) ) {
                System.out.println( "DNA not detected" );
                return false;
            }
            // (6) a label needing quoting must come out IDENTICAL in both blocks, or the file is unreadable
            final Phylogeny p6 = tree( "(X,Y)" );
            p6.getExternalNodes().get( 0 ).setName( "Seba's bat" );
            p6.getExternalNodes().get( 1 ).setName( "a b" );
            setSeq( p6, 0, "MKAL" );
            setSeq( p6, 1, "MKIV" );
            final String n6 = nexus( p6 );
            final String b6 = charactersBlock( n6 );
            final String t6 = taxaBlock( n6 );
            // The quoting style is ForesterUtil.santitizeStringForNH's business, and it is shared with the Newick
            // writer -- this test must not pin it. What MUST hold is that the matrix uses the very same token as
            // TaxLabels: a row whose label is not a declared taxon makes the file unreadable.
            if ( t6.contains( "Seba" ) == false ) {
                System.out.println( "TaxLabels lost the name entirely: " + t6 );
                return false;
            }
            final String labels6 = t6.substring( t6.indexOf( "TaxLabels" ) + "TaxLabels".length(),
                                                 t6.indexOf( ";", t6.indexOf( "TaxLabels" ) ) ).trim();
            // both names need quoting, so each declared label is a quoted token
            int quoted = 0;
            for( final String tok : labels6.split( "(?<=[\"'])\\s+(?=[\"'])" ) ) {
                final String label = tok.trim();
                if ( label.length() == 0 ) {
                    continue;
                }
                ++quoted;
                if ( !b6.contains( label ) ) {
                    System.out.println( "the matrix label does not match TaxLabels for [" + label + "]: " + b6 );
                    return false;
                }
            }
            if ( quoted != 2 ) {
                System.out.println( "expected two quoted labels, parsed " + quoted + " from: " + labels6 );
                return false;
            }
            // and a name needing quotes must still survive a round trip through our own parser
            final File tmp6 = File.createTempFile( "aptx_nexus_q_", ".nex" );
            tmp6.deleteOnExit();
            new PhylogenyWriter().toNexus( tmp6, p6, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
            final NexusPhylogeniesParser pq = new NexusPhylogeniesParser();
            pq.setSource( tmp6 );
            final Phylogeny back6 = pq.parse()[ 0 ];
            boolean found_seba = false;
            for( final PhylogenyNode n : back6.getExternalNodes() ) {
                if ( "Seba's bat".equals( n.getName() ) ) {
                    found_seba = true;
                }
            }
            if ( !found_seba ) {
                System.out.println( "a quoted name did not survive the round trip with a Characters block" );
                return false;
            }
            // (7) the label falls back the same way in both blocks: no name, but a taxonomy code
            final Phylogeny p7 = tree( "(,)" );
            final Taxonomy tax = new Taxonomy();
            tax.setTaxonomyCode( "HUMAN" );
            p7.getExternalNodes().get( 0 ).getNodeData().setTaxonomy( tax );
            setSeq( p7, 0, "MKAL" );
            setSeq( p7, 1, "MKIV" );
            final String n7 = nexus( p7 );
            if ( !taxaBlock( n7 ).contains( "HUMAN" ) || !charactersBlock( n7 ).contains( "HUMAN" ) ) {
                System.out.println( "the taxonomy-code fallback did not reach both blocks: " + n7 );
                return false;
            }
            // (8) deliberate non-behaviour: an internal node's sequence is not a taxon and is not written
            final Phylogeny p8 = tree( "((A,B)INNER,C)" );
            setSeq( p8, 0, "MKAL" );
            setSeq( p8, 1, "MKIV" );
            setSeq( p8, 2, "MKLL" );
            for( final PhylogenyNode n : p8.getNodes( "INNER" ) ) {
                final Sequence s = new Sequence();
                s.setMolecularSequence( "WWWW" );
                n.getNodeData().setSequence( s );
            }
            final String b8 = charactersBlock( nexus( p8 ) );
            if ( b8.contains( "WWWW" ) ) {
                System.out.println( "an internal node's sequence must not appear in the matrix: " + b8 );
                return false;
            }
            if ( !b8.contains( "NChar=4" ) ) {
                System.out.println( "internal sequence changed NChar: " + b8 );
                return false;
            }
            // (9) the whole file must still be readable as a tree: adding data must not break the trees block
            final File tmp = File.createTempFile( "aptx_nexus_seq_", ".nex" );
            tmp.deleteOnExit();
            new PhylogenyWriter().toNexus( tmp, p1, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
            final NexusPhylogeniesParser parser = new NexusPhylogeniesParser();
            parser.setSource( tmp );
            final Phylogeny[] back = parser.parse();
            if ( ( back.length != 1 ) || ( back[ 0 ].getNumberOfExternalNodes() != 3 ) ) {
                System.out.println( "the tree no longer round trips with a Characters block present" );
                return false;
            }
            // (10) the sequences themselves must survive the round trip, and come back marked as aligned
            int back_with_seq = 0;
            for( final PhylogenyNode n : back[ 0 ].getExternalNodes() ) {
                if ( n.getNodeData().isHasSequence()
                        && !org.forester.util.ForesterUtil
                                .isEmpty( n.getNodeData().getSequence().getMolecularSequence() ) ) {
                    ++back_with_seq;
                    if ( !n.getNodeData().getSequence().isMolecularSequenceAligned() ) {
                        System.out.println( "a sequence from a Nexus matrix must be marked aligned" );
                        return false;
                    }
                    if ( n.getNodeData().getSequence().getMolecularSequence().length() != 7 ) {
                        System.out.println( "wrong sequence length after the round trip: "
                                + n.getNodeData().getSequence().getMolecularSequence() );
                        return false;
                    }
                }
            }
            if ( back_with_seq != 3 ) {
                System.out.println( "expected 3 sequences back from the matrix, got " + back_with_seq );
                return false;
            }
            // (11) Christian, 2026-09-23: the NAME wins over a sequence accession, and the label must be
            // the same in all three blocks (a matrix row naming a taxon the tree does not contain is
            // unreadable). Archaeopteryx.js writes the name too.
            final Phylogeny p11 = tree( "(A,B)" );
            final Sequence acc_seq = new Sequence();
            acc_seq.setAccession( new org.forester.phylogeny.data.Accession( "P12345", "UniProt" ) );
            acc_seq.setMolecularSequence( "MKAL" );
            p11.getExternalNodes().get( 0 ).getNodeData().setSequence( acc_seq );
            setSeq( p11, 1, "MKIV" );
            final String n11 = nexus( p11 );
            if ( n11.contains( "P12345" ) ) {
                System.out.println( "the accession beat the node name: " + n11 );
                return false;
            }
            for( final String block : new String[] { taxaBlock( n11 ), charactersBlock( n11 ),
                                                     n11.substring( n11.indexOf( "Begin Trees;" ) ) } ) {
                if ( !block.contains( "A" ) ) {
                    System.out.println( "a block does not carry the node name: " + block );
                    return false;
                }
            }
            // and that tree must still round trip its sequences
            final File t11 = File.createTempFile( "aptx_nexus_acc_", ".nex" );
            t11.deleteOnExit();
            new PhylogenyWriter().toNexus( t11, p11, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
            final NexusPhylogeniesParser pa = new NexusPhylogeniesParser();
            pa.setSource( t11 );
            final Phylogeny b11 = pa.parse()[ 0 ];
            int acc_seqs = 0;
            for( final PhylogenyNode n : b11.getExternalNodes() ) {
                if ( n.getNodeData().isHasSequence() && !org.forester.util.ForesterUtil
                        .isEmpty( n.getNodeData().getSequence().getMolecularSequence() ) ) {
                    ++acc_seqs;
                }
            }
            if ( acc_seqs != 2 ) {
                System.out.println( "a name-labelled tree lost its sequences on re-read: " + acc_seqs );
                return false;
            }
            // (11b) but an accession is still a real identifier: a tip that NOTHING else names must be
            // labelled with it rather than replaced by an opaque placeholder
            final Phylogeny p11b = tree( "(,B)" );
            final Sequence only_acc = new Sequence();
            only_acc.setAccession( new org.forester.phylogeny.data.Accession( "Q99999", "UniProt" ) );
            only_acc.setMolecularSequence( "MKAL" );
            p11b.getExternalNodes().get( 0 ).getNodeData().setSequence( only_acc );
            setSeq( p11b, 1, "MKIV" );
            final String n11b = nexus( p11b );
            if ( !taxaBlock( n11b ).contains( "Q99999" ) || !charactersBlock( n11b ).contains( "Q99999" ) ) {
                System.out.println( "an accession-only tip was not labelled with its accession: " + n11b );
                return false;
            }
            if ( n11b.contains( "node1" ) ) {
                System.out.println( "a placeholder displaced a real accession: " + n11b );
                return false;
            }
            // (11c) a tip that nothing names at all gets "nodeN" by tip index -- in every block and in the
            // trees string, or the file cannot be read back. Before this, such a tip produced an empty
            // label, TaxLabels shorter than NTax, and a matrix row that was silently dropped on re-read.
            final Phylogeny p11c = tree( "(,)" );
            final Taxonomy tx = new Taxonomy();
            tx.setTaxonomyCode( "HUMAN" );
            p11c.getExternalNodes().get( 0 ).getNodeData().setTaxonomy( tx );
            setSeq( p11c, 0, "MKAL" );
            setSeq( p11c, 1, "MKIV" );
            final String n11c = nexus( p11c );
            for( final String block : new String[] { taxaBlock( n11c ), charactersBlock( n11c ),
                                                     n11c.substring( n11c.indexOf( "Begin Trees;" ) ) } ) {
                if ( !block.contains( "node2" ) ) {
                    System.out.println( "the placeholder is missing from a block: " + block );
                    return false;
                }
            }
            if ( n11c.contains( "node1" ) ) {
                System.out.println( "a placeholder displaced the HUMAN taxonomy code: " + n11c );
                return false;
            }
            final File t11c = File.createTempFile( "aptx_nexus_ph_", ".nex" );
            t11c.deleteOnExit();
            new PhylogenyWriter().toNexus( t11c, p11c, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
            final NexusPhylogeniesParser pc = new NexusPhylogeniesParser();
            pc.setSource( t11c );
            final Phylogeny b11c = pc.parse()[ 0 ];
            int kept = 0;
            for( final PhylogenyNode n : b11c.getExternalNodes() ) {
                if ( n.getNodeData().isHasSequence() && !org.forester.util.ForesterUtil
                        .isEmpty( n.getNodeData().getSequence().getMolecularSequence() ) ) {
                    ++kept;
                }
            }
            if ( kept != 2 ) {
                System.out.println( "a nameless tip lost its sequence on re-read: " + kept );
                return false;
            }
            // (11d) the same rule reaches plain New Hampshire, which is where the trees block comes from
            if ( !p11c.toNewHampshire( NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ).contains( "node2" ) ) {
                System.out.println( "plain Newick still writes a nameless tip as an empty label: "
                        + p11c.toNewHampshire( NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ) );
                return false;
            }
            // (11f) the placeholder numbers TIP POSITIONS, not placeholders. The two rules agree on every
            // tree whose nameless tips come first, so only a NAMED tip before a nameless one separates
            // them: "(a,)" is node2 under the real rule and node1 under the other. Found by the
            // Archaeopteryx.js session, whose own sabotage survived for exactly this reason.
            final Phylogeny p11f = tree( "(a,)" );
            final String nh11f = p11f.toNewHampshire( NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
            if ( !nh11f.contains( "node2" ) || nh11f.contains( "node1" ) ) {
                System.out.println( "the placeholder is numbered by placeholder, not by tip position: "
                        + nh11f );
                return false;
            }
            final String n11f = nexus( p11f );
            if ( !taxaBlock( n11f ).contains( "node2" ) || taxaBlock( n11f ).contains( "node1" ) ) {
                System.out.println( "TaxLabels numbers the placeholder differently from the tree: " + n11f );
                return false;
            }
            // the same shape the other way round, where both rules DO agree, as the neighbouring case
            if ( !tree( "(,a)" ).toNewHampshire( NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ).contains( "node1" ) ) {
                System.out.println( "a leading nameless tip should be node1" );
                return false;
            }
            // (11j) NHX gives a nameless TIP the same placeholder as New Hampshire, so the two writers no
            // longer disagree about the same tree's tip names (Christian, 2026-09-23).
            for( final String nh : new String[] { "(,)", "(a,)", "(((,),),)" } ) {
                final Phylogeny t = tree( nh );
                final String plain = t.toNewHampshire( NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
                final String nhx = t.toNewHampshireX();
                if ( !plain.equals( nhx + ";" ) ) {
                    System.out.println( "NH and NHX disagree about " + nh + ": [" + plain + "] vs [" + nhx
                            + "]" );
                    return false;
                }
            }
            // (11k) but NHX takes ONLY the placeholder, not the rest of the chain: its tags already carry a
            // node's taxonomy and sequence, so labelling from them would duplicate the tags AND give a name
            // to a node that never had one. A tagged tree must come back byte for byte.
            final String tagged = "((((a,b),c),d)[&&NHX:S=lizards],e[&&NHX:S=reptiles])r[&&NHX:S=animals]";
            final String rewritten = Phylogeny.createInstanceFromNhxString( tagged ).toNewHampshireX();
            if ( !tagged.equals( rewritten ) ) {
                System.out.println( "NHX is no longer byte-faithful for a tagged tree:\n  in  " + tagged
                        + "\n  out " + rewritten );
                return false;
            }
            // and that same tree in NH DOES take the chain, which is the difference being kept on purpose
            if ( !Phylogeny.createInstanceFromNhxString( tagged )
                    .toNewHampshire( NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ).contains( "lizards" ) ) {
                System.out.println( "New Hampshire stopped labelling a node from its taxonomy" );
                return false;
            }
            // (11e) deliberate non-behaviour: an unlabeled INTERNAL node is normal and keeps no
            // placeholder. Our own writer only ever offers one for external nodes, so this is checked
            // through the public four-argument call, which is the only way to reach the guard.
            final Phylogeny p11e = tree( "((A,B),C)" );
            if ( p11e.toNewHampshire( NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ).contains( "node" ) ) {
                System.out.println( "an internal node was given a placeholder: "
                        + p11e.toNewHampshire( NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ) );
                return false;
            }
            final PhylogenyNode inner = p11e.getRoot().getChildNode( 0 );
            if ( inner.isExternal() ) {
                System.out.println( "the fixture's first child is not an internal node" );
                return false;
            }
            if ( inner.toNewHampshire( false, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE, false, "PLACEHOLDER" )
                    .contains( "PLACEHOLDER" ) ) {
                System.out.println( "an internal node took a placeholder it was offered" );
                return false;
            }
            // and an external one does take it
            if ( !tree( "(,B)" ).getExternalNodes().get( 0 )
                    .toNewHampshire( false, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE, false, "PLACEHOLDER" )
                    .contains( "PLACEHOLDER" ) ) {
                System.out.println( "a nameless external node refused the placeholder" );
                return false;
            }
            // (11g) the datatype is decided by the WHOLE matrix. guessMolecularSequenceType looks for
            // residues only protein has, so a short protein built from nucleotide letters guesses DNA --
            // and a matrix wrongly declared DNA is read back with every other residue replaced by N.
            final Phylogeny p11g = tree( "(A,B)" );
            setSeq( p11g, 0, "MKATSWNP" );
            setSeq( p11g, 1, "MKLTSWNP" );
            if ( !charactersBlock( nexus( p11g ) ).contains( "DataType=Protein" ) ) {
                System.out.println( "one sequence's guess decided the datatype for the matrix: "
                        + charactersBlock( nexus( p11g ) ) );
                return false;
            }
            final File t11g = File.createTempFile( "aptx_nexus_dt_", ".nex" );
            t11g.deleteOnExit();
            new PhylogenyWriter().toNexus( t11g, p11g, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
            final NexusPhylogeniesParser pg = new NexusPhylogeniesParser();
            pg.setSource( t11g );
            // compared against the ORIGINAL text, not against a residue letter: MKATSWNP contains an N of
            // its own (asparagine), so "did an N appear" cannot tell a destroyed residue from a real one
            final String[] dt_expected = { "MKATSWNP", "MKLTSWNP" };
            int dt_i = 0;
            for( final PhylogenyNode n : pg.parse()[ 0 ].getExternalNodes() ) {
                final String dt_back = n.getNodeData().getSequence().getMolecularSequence();
                if ( !dt_expected[ dt_i++ ].equals( dt_back ) ) {
                    System.out.println( "a residue was destroyed by a wrong datatype: wrote "
                            + dt_expected[ dt_i - 1 ] + ", read back " + dt_back );
                    return false;
                }
            }
            // the rule this fix actually implements is "protein wins a DISAGREEMENT", which needs a matrix
            // where one row guesses DNA and another guesses protein -- two protein rows cannot test it.
            // (Before F/P/V were added, MKATSWNP guessed DNA and this case was exercised by accident; it
            // stopped being exercised the moment P began marking protein.)
            final Phylogeny p11g_mixed = tree( "(A,B)" );
            setSeq( p11g_mixed, 0, "ACGTACGTACGT" );
            setSeq( p11g_mixed, 1, "MKLTSWNPMKLT" );
            if ( !charactersBlock( nexus( p11g_mixed ) ).contains( "DataType=Protein" ) ) {
                System.out.println( "a matrix holding both a nucleotide-looking and a protein row was not "
                        + "typed Protein: " + charactersBlock( nexus( p11g_mixed ) ) );
                return false;
            }
            // the neighbouring case: a genuine nucleotide matrix must still be typed DNA
            final Phylogeny p11g2 = tree( "(A,B)" );
            setSeq( p11g2, 0, "ACGTACGTACGTACGT" );
            setSeq( p11g2, 1, "ACGTACGTACGTACGA" );
            if ( !charactersBlock( nexus( p11g2 ) ).contains( "DataType=DNA" ) ) {
                System.out.println( "a nucleotide matrix is no longer typed DNA" );
                return false;
            }

            // (11l) the datatype alphabet, pinned letter by letter and jointly with Archaeopteryx.js: a
            // letter added on one side only types the same file two ways, which is worse than a blind spot
            // both share. Each letter is tested as "<letter>T": a bare letter that gets no verdict falls
            // through to Protein and would look identical to one tested AS protein, so the T gives the
            // undecided letters a DNA answer and makes the two distinguishable.
            final String protein_letters = "LIEHDQFPV";
            for( final char c : protein_letters.toCharArray() ) {
                final Phylogeny t = tree( "(A,B)" );
                setSeq( t, 0, c + "TTT" );
                setSeq( t, 1, c + "TTT" );
                if ( !charactersBlock( nexus( t ) ).contains( "DataType=Protein" ) ) {
                    System.out.println( "'" + c + "' no longer marks a sequence as protein" );
                    return false;
                }
            }
            // and the letters that must NOT decide protein, or a nucleotide alignment gets called one
            for( final char c : "ACGRYMKWSN".toCharArray() ) {
                final Phylogeny t = tree( "(A,B)" );
                setSeq( t, 0, c + "TTT" );
                setSeq( t, 1, c + "TTT" );
                if ( !charactersBlock( nexus( t ) ).contains( "DataType=DNA" ) ) {
                    System.out.println( "'" + c + "' now marks a sequence as protein, but it is a "
                            + "nucleotide code" );
                    return false;
                }
            }
            // the case that motivated adding P: its only protein-exclusive residue
            final Phylogeny p11l = tree( "(A,B)" );
            setSeq( p11l, 0, "MKATSWNP" );
            setSeq( p11l, 1, "MKATSWNP" );
            if ( !charactersBlock( nexus( p11l ) ).contains( "DataType=Protein" ) ) {
                System.out.println( "a protein whose only exclusive residue is P is still typed DNA" );
                return false;
            }

            // (11h) a tip literally named "node2" must not be given away to a different tip: two taxa under
            // one label is illegal Nexus, and the parser reads the repeat as an interleaved continuation,
            // so BOTH tips come back with the two sequences concatenated.
            final Phylogeny p11h = tree( "(node2,)" );
            setSeq( p11h, 0, "MKAL" );
            setSeq( p11h, 1, "MKIV" );
            final File t11h = File.createTempFile( "aptx_nexus_coll_", ".nex" );
            t11h.deleteOnExit();
            new PhylogenyWriter().toNexus( t11h, p11h, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
            final NexusPhylogeniesParser ph = new NexusPhylogeniesParser();
            ph.setSource( t11h );
            for( final PhylogenyNode n : ph.parse()[ 0 ].getExternalNodes() ) {
                final String coll_back = n.getNodeData().getSequence().getMolecularSequence();
                if ( coll_back.length() != 4 ) {
                    System.out.println( "a placeholder collided with a real tip name, concatenating "
                            + "sequences: " + coll_back );
                    return false;
                }
            }

            // (11i) when two tips genuinely share a label the matrix cannot be keyed on it, so none is
            // written -- the Taxa and Trees blocks have always written such a tree, but a matrix would turn
            // a cosmetic problem into a corrupting one.
            final Phylogeny p11i = tree( "((A,A),B)" );
            setSeq( p11i, 0, "MKAL" );
            setSeq( p11i, 1, "MKIV" );
            setSeq( p11i, 2, "MKLL" );
            final String n11i = nexus( p11i );
            if ( charactersBlock( n11i ) != null ) {
                System.out.println( "a matrix was written for a tree with duplicate taxon labels" );
                return false;
            }
            if ( !n11i.contains( "share the taxon label" ) ) {
                System.out.println( "no comment explaining the missing matrix" );
                return false;
            }

            // (12) a row of nothing but the missing symbol is absence of data, not a sequence: it must not
            // come back as one, or a round trip would invent sequences for tips that never had any
            final Phylogeny p12 = tree( "(A,B)" );
            setSeq( p12, 0, "MKAL" );
            final File t12 = File.createTempFile( "aptx_nexus_miss_", ".nex" );
            t12.deleteOnExit();
            new PhylogenyWriter().toNexus( t12, p12, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
            final NexusPhylogeniesParser pm = new NexusPhylogeniesParser();
            pm.setSource( t12 );
            final Phylogeny b12 = pm.parse()[ 0 ];
            for( final PhylogenyNode n : b12.getExternalNodes() ) {
                final boolean has = n.getNodeData().isHasSequence() && !org.forester.util.ForesterUtil
                        .isEmpty( n.getNodeData().getSequence().getMolecularSequence() );
                if ( "B".equals( n.getName() ) && has ) {
                    System.out.println( "an all-missing row came back as a sequence: ["
                            + n.getNodeData().getSequence().getMolecularSequence() + "]" );
                    return false;
                }
                if ( "A".equals( n.getName() ) && !has ) {
                    System.out.println( "the real sequence was lost" );
                    return false;
                }
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    public static void main( final String[] args ) {
        if ( test() ) {
            System.out.println( "NexusSequenceExportTest: OK." );
        }
        else {
            System.out.println( "NexusSequenceExportTest: FAILED." );
        }
    }
}
