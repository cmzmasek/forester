// Round trips between phyloXML and Nexus, for trees carrying node names, branch lengths, support values
// and molecular sequences.
//
//   (1) phyloXML -> Nexus -> phyloXML
//   (2) Nexus    -> phyloXML -> Nexus
//
// Both must come back with the same tree and the same sequences.

package org.forester.io.writers;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class NexusPhyloXmlRoundTripTest {

    // Support only survives a Newick-based format when it is actually written. IN_SQUARE_BRACKETS is the one
    // style that both writes it and reads it back as a confidence: NONE never writes it, and
    // AS_INTERNAL_NODE_NAMES writes it where an internal node's NAME goes, so it returns as a name.
    private static final NH_CONVERSION_SUPPORT_VALUE_STYLE SVS = NH_CONVERSION_SUPPORT_VALUE_STYLE.IN_SQUARE_BRACKETS;

    private static Phylogeny fixture() throws Exception {
        final Phylogeny p = Phylogeny
                .createInstanceFromNhxString( "((A:0.1,B:0.2):0.3,(C:0.15,D:0.25):0.35)" );
        // A NAMED tree: an unnamed one is written as the generated "tree1", which comes back as a real name
        // and is then quoted on the next write -- a difference in the file that is not a difference in data.
        p.setName( "round trip fixture" );
        final String[] seqs = { "MKAL-IVG", "MKAL-IVA", "MKAL-IVC", "MKAL-IVD" };
        int i = 0;
        for( final PhylogenyNode n : p.getExternalNodes() ) {
            final Sequence s = new Sequence();
            s.setMolecularSequence( seqs[ i++ ] );
            s.setMolecularSequenceAligned( true );
            n.getNodeData().addSequence( s );
        }
        final double[] support = { 90, 75 };
        int j = 0;
        for( final PhylogenyNodeIterator it = p.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.isRoot() ) {
                n.getBranchData().addConfidence( new Confidence( support[ j++ % 2 ], "bootstrap" ) );
            }
        }
        return p;
    }

    /** Names, branch lengths, support values and molecular sequences -- the four things under test. */
    private static String canonical( final Phylogeny p ) {
        final List<String> tips = new ArrayList<String>();
        for( final PhylogenyNode n : p.getExternalNodes() ) {
            final String mol = n.getNodeData().isHasSequence()
                    ? n.getNodeData().getSequence().getMolecularSequence() : null;
            tips.add( n.getName() + "|" + n.getDistanceToParent() + "|" + ( mol == null ? "-" : mol ) );
        }
        Collections.sort( tips );
        final List<String> internals = new ArrayList<String>();
        for( final PhylogenyNodeIterator it = p.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isExternal() ) {
                continue;
            }
            String c = "-";
            if ( n.getBranchData().isHasConfidences() && ( n.getBranchData().getConfidences().size() > 0 ) ) {
                c = String.valueOf( n.getBranchData().getConfidence( 0 ).getValue() );
            }
            internals.add( "[" + n.getName() + "|" + n.getDistanceToParent() + "|" + c + "]" );
        }
        return "tips=" + tips + " internal=" + internals;
    }

    private static File writePhyloXml( final Phylogeny p ) throws Exception {
        final File f = File.createTempFile( "aptx_rt_px_", ".xml" );
        f.deleteOnExit();
        new PhylogenyWriter().toPhyloXML( f, p, 0 );
        return f;
    }

    private static Phylogeny readPhyloXml( final File f ) throws Exception {
        return PhylogenyMethods.readPhylogenies( PhyloXmlParser.createPhyloXmlParserXsdValidating(), f )[ 0 ];
    }

    private static File writeNexus( final Phylogeny p, final NH_CONVERSION_SUPPORT_VALUE_STYLE svs )
            throws Exception {
        final File f = File.createTempFile( "aptx_rt_nex_", ".nex" );
        f.deleteOnExit();
        new PhylogenyWriter().toNexus( f, p, svs );
        return f;
    }

    private static Phylogeny readNexus( final File f ) throws Exception {
        final NexusPhylogeniesParser parser = new NexusPhylogeniesParser();
        parser.setSource( f );
        return parser.parse()[ 0 ];
    }

    private static String text( final File f ) throws Exception {
        return new String( Files.readAllBytes( f.toPath() ), StandardCharsets.UTF_8 );
    }

    public static boolean test() {
        try {
            // (0) the fixture must really carry all four kinds of data, or the round trips below would be
            // comparing nothing to nothing and would pass however broken the writers were.
            final Phylogeny fix = fixture();
            int named = 0, with_len = 0, with_seq = 0, with_support = 0;
            for( final PhylogenyNode n : fix.getExternalNodes() ) {
                if ( ( n.getName() != null ) && ( n.getName().length() > 0 ) ) {
                    ++named;
                }
                if ( n.getDistanceToParent() > 0 ) {
                    ++with_len;
                }
                if ( n.getNodeData().isHasSequence()
                        && ( n.getNodeData().getSequence().getMolecularSequence() != null ) ) {
                    ++with_seq;
                }
            }
            for( final PhylogenyNodeIterator it = fix.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isExternal() && n.getBranchData().isHasConfidences() ) {
                    ++with_support;
                }
            }
            if ( ( named != 4 ) || ( with_len != 4 ) || ( with_seq != 4 ) || ( with_support != 2 ) ) {
                System.out.println( "the fixture does not carry what these tests claim to check: names="
                        + named + " lengths=" + with_len + " seqs=" + with_seq + " support=" + with_support );
                return false;
            }

            // (1) phyloXML -> Nexus -> phyloXML
            final Phylogeny start1 = readPhyloXml( writePhyloXml( fixture() ) );
            final Phylogeny mid1 = readNexus( writeNexus( start1, SVS ) );
            final Phylogeny end1 = readPhyloXml( writePhyloXml( mid1 ) );
            if ( !canonical( start1 ).equals( canonical( end1 ) ) ) {
                System.out.println( "phyloXML -> Nexus -> phyloXML changed the tree" );
                System.out.println( "  before: " + canonical( start1 ) );
                System.out.println( "  after : " + canonical( end1 ) );
                return false;
            }
            // and the sequences must still be marked aligned at the end
            for( final PhylogenyNode n : end1.getExternalNodes() ) {
                if ( !n.getNodeData().getSequence().isMolecularSequenceAligned() ) {
                    System.out.println( "the aligned flag was lost on the way round" );
                    return false;
                }
            }

            // (2) Nexus -> phyloXML -> Nexus, compared as BYTES: the strongest form this direction allows
            final File nex1 = writeNexus( fixture(), SVS );
            final Phylogeny mid2 = readPhyloXml( writePhyloXml( readNexus( nex1 ) ) );
            final File nex2 = writeNexus( mid2, SVS );
            if ( !text( nex1 ).equals( text( nex2 ) ) ) {
                System.out.println( "Nexus -> phyloXML -> Nexus is not byte identical" );
                final String[] a = text( nex1 ).split( "\\R" );
                final String[] b = text( nex2 ).split( "\\R" );
                for( int i = 0; i < Math.max( a.length, b.length ); ++i ) {
                    final String x = i < a.length ? a[ i ] : "<none>";
                    final String y = i < b.length ? b[ i ] : "<none>";
                    if ( !x.equals( y ) ) {
                        System.out.println( "  line " + ( i + 1 ) + "\n    A: " + x + "\n    B: " + y );
                    }
                }
                return false;
            }

            // (3) what Nexus cannot carry, pinned so that a change is noticed rather than discovered later:
            // a square-bracket support value has nowhere to record its TYPE, so "bootstrap" comes back as
            // "unknown". The VALUE survives, which is what (1) and (2) check.
            String type_after = null;
            for( final PhylogenyNodeIterator it = end1.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isExternal() && n.getBranchData().isHasConfidences() ) {
                    type_after = n.getBranchData().getConfidence( 0 ).getType();
                    break;
                }
            }
            if ( "bootstrap".equals( type_after ) ) {
                System.out.println( "the confidence type now survives Nexus -- good, but this test and the "
                        + "comment above it need updating" );
                return false;
            }

            // (4) deliberate non-behaviour: with NONE the support is not written at all, so it cannot come
            // back. This was the Options.init default until 2026-09-23, which made a plain Save As Nexus
            // drop support values silently; the default is IN_SQUARE_BRACKETS now, so the round trip above
            // is what a user actually gets. NONE remains selectable, and still loses support.
            final Phylogeny none_rt = readNexus( writeNexus( fixture(), NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ) );
            for( final PhylogenyNodeIterator it = none_rt.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isExternal() && n.getBranchData().isHasConfidences()
                        && ( n.getBranchData().getConfidences().size() > 0 ) ) {
                    System.out.println( "support survived style NONE -- the default now keeps it, so this "
                            + "test needs updating" );
                    return false;
                }
            }

            // (5) a tip with NO sequence must not acquire one: the writer gives it a row of the missing
            // symbol to keep the matrix rectangular, and that row is not data.
            final Phylogeny partial = fixture();
            partial.getExternalNodes().get( 1 ).getNodeData().setSequence( null );
            final Phylogeny partial_rt = readPhyloXml( writePhyloXml( readNexus( writeNexus( partial, SVS ) ) ) );
            int seqs_back = 0;
            for( final PhylogenyNode n : partial_rt.getExternalNodes() ) {
                if ( n.getNodeData().isHasSequence()
                        && ( n.getNodeData().getSequence().getMolecularSequence() != null )
                        && ( n.getNodeData().getSequence().getMolecularSequence().length() > 0 ) ) {
                    ++seqs_back;
                }
            }
            if ( seqs_back != 3 ) {
                System.out.println( "expected 3 sequences after a round trip with one tip lacking one, got "
                        + seqs_back );
                return false;
            }

            // (6) gaps are part of an alignment and must survive verbatim
            for( final PhylogenyNode n : end1.getExternalNodes() ) {
                if ( n.getNodeData().getSequence().getMolecularSequence().indexOf( '-' ) < 0 ) {
                    System.out.println( "a gap character was lost from the alignment" );
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
            System.out.println( "NexusPhyloXmlRoundTripTest: OK." );
        }
        else {
            System.out.println( "NexusPhyloXmlRoundTripTest: FAILED." );
        }
    }
}
