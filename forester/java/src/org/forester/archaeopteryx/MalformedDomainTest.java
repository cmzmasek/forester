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

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * A malformed protein domain must never cost the user the tree.
 * <p>
 * {@code ProteinDomain}'s constructor rejects a zero-or-negative length ({@code from >= to}) with an UNCHECKED
 * IllegalArgumentException, which used to sail past the phyloXML parser's checked-exception handling and out of the
 * factory -- so one degenerate domain aborted the ENTIRE file load: the tree, every other node, and every other
 * domain beside it. Worse, the phyloXML schema constrains {@code from}/{@code to} only to nonNegativeInteger, so such
 * a file is SCHEMA-VALID and the XSD-validating parser met exactly the same fate: Archaeopteryx refused a document
 * that conforms to phyloXML. A domain is a display annotation and the tree is the data, so now the annotation is what
 * gets dropped -- skipped, counted, and reported once at load.
 */
public final class MalformedDomainTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MalformedDomain: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    /** Two tips; node A carries {@code domains}, node B one always-valid domain. */
    private static String xml( final String domains ) {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                + "<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
                + "xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" "
                + "xmlns=\"http://www.phyloxml.org\"><phylogeny rooted=\"true\"><clade>"
                + "<clade><name>A</name><branch_length>0.1</branch_length><sequence>"
                + "<domain_architecture length=\"100\">" + domains + "</domain_architecture>"
                + "</sequence></clade>"
                + "<clade><name>B</name><branch_length>0.2</branch_length><sequence>"
                + "<domain_architecture length=\"100\">"
                + "<domain from=\"10\" to=\"70\" confidence=\"1e-8\">Fine</domain>"
                + "</domain_architecture></sequence></clade>"
                + "</clade></phylogeny></phyloxml>";
    }

    private static final String GOOD_A = "<domain from=\"5\" to=\"40\" confidence=\"1e-8\">Good</domain>";
    private static final String GOOD_B = "<domain from=\"60\" to=\"90\" confidence=\"1e-7\">AlsoGood</domain>";

    private static Phylogeny[] load( final String xml, final boolean validating ) throws Exception {
        final PhyloXmlParser p = validating ? PhyloXmlParser.createPhyloXmlParserXsdValidating()
                : PhyloXmlParser.createPhyloXmlParser();
        return ParserBasedPhylogenyFactory.getInstance().create( new StringBuffer( xml ), p );
    }

    private static DomainArchitecture architectureOf( final Phylogeny phy, final String node_name ) {
        final PhylogenyNode n = phy.getNode( node_name );
        return n.getNodeData().getSequence().getDomainArchitecture();
    }

    public static boolean test() {
        try {
            // ---- the reported case: a zero-length domain (from == to) --------------------------------------
            // Both parsers, because the file is schema-valid: the XSD-validating one must not be the strict one.
            for ( final boolean validating : new boolean[] { false, true } ) {
                final Phylogeny[] phys = load( xml( GOOD_A
                        + "<domain from=\"50\" to=\"50\" confidence=\"1e-9\">ZeroLen</domain>" + GOOD_B ),
                                               validating );
                final String who = validating ? "XSD-validating parser" : "parser";
                if ( ( phys == null ) || ( phys.length != 1 ) ) {
                    return fail( who + ": a zero-length domain must not abort the load" );
                }
                if ( phys[ 0 ].getNumberOfExternalNodes() != 2 ) {
                    return fail( who + ": the whole tree must survive one bad domain, got "
                            + phys[ 0 ].getNumberOfExternalNodes() + " tips" );
                }
                final DomainArchitecture a = architectureOf( phys[ 0 ], "A" );
                // the GOOD domains beside it must be kept -- dropping the architecture wholesale would be its own bug
                if ( a.getDomains().size() != 2 ) {
                    return fail( who + ": the two valid domains beside the bad one must be kept, got "
                            + a.getDomains().size() );
                }
                if ( a.getIgnoredDomainCount() != 1 ) {
                    return fail( who + ": the skipped domain must be counted, got " + a.getIgnoredDomainCount() );
                }
                // an untouched node must not inherit the count
                if ( architectureOf( phys[ 0 ], "B" ).getIgnoredDomainCount() != 0 ) {
                    return fail( who + ": a node with only valid domains must report none ignored" );
                }
                if ( MainFrame.countIgnoredDomains( phys ) != 1 ) {
                    return fail( who + ": countIgnoredDomains must sum across the tree" );
                }
            }
            // ---- the other shapes of "cannot be built" -----------------------------------------------------
            // inverted (from > to): no sensible reading, must be dropped just the same
            final Phylogeny[] inverted = load( xml( GOOD_A
                    + "<domain from=\"80\" to=\"60\" confidence=\"1e-9\">Inverted</domain>" ), false );
            if ( MainFrame.countIgnoredDomains( inverted ) != 1 ) {
                return fail( "an inverted domain (from > to) must be skipped and counted" );
            }
            if ( architectureOf( inverted[ 0 ], "A" ).getDomains().size() != 1 ) {
                return fail( "the valid domain beside an inverted one must be kept" );
            }
            // a missing coordinate takes the parser's OWN checked-exception path -- also non-fatal now
            final Phylogeny[] no_to = load( xml( GOOD_A
                    + "<domain from=\"50\" confidence=\"1e-9\">NoTo</domain>" ), false );
            if ( MainFrame.countIgnoredDomains( no_to ) != 1 ) {
                return fail( "a domain with a missing \"to\" must be skipped and counted" );
            }
            // ---- and NO false positives on a clean file ----------------------------------------------------
            final Phylogeny[] clean = load( xml( GOOD_A + GOOD_B ), false );
            if ( MainFrame.countIgnoredDomains( clean ) != 0 ) {
                return fail( "a clean file must report nothing ignored, got "
                        + MainFrame.countIgnoredDomains( clean ) );
            }
            if ( architectureOf( clean[ 0 ], "A" ).getDomains().size() != 2 ) {
                return fail( "a clean file must keep all its domains" );
            }
            // ---- the accounting helper's own edges ---------------------------------------------------------
            if ( MainFrame.countIgnoredDomains( null ) != 0 ) {
                return fail( "countIgnoredDomains(null) must be 0, not a crash" );
            }
            if ( MainFrame.countIgnoredDomains( new Phylogeny[] { null } ) != 0 ) {
                return fail( "countIgnoredDomains must tolerate a null tree in the array" );
            }
            // An EMPTY tree is a reachable state -- a phyloXML file may hold one alongside real trees, and
            // addPhylogeniesToTabs explicitly skips them -- and iteratorPreorder() THROWS on one. Since the
            // report is the FIRST thing called after the tabs are built, letting that escape would take out the
            // whole post-load chain (label extraction, tip dates, the time-tree offer): the very class of
            // failure this change exists to remove.
            final Phylogeny[] with_empty = load( xml( GOOD_A
                    + "<domain from=\"50\" to=\"50\" confidence=\"1e-9\">ZeroLen</domain>" ), false );
            if ( MainFrame.countIgnoredDomains( new Phylogeny[] { new Phylogeny(), with_empty[ 0 ] } ) != 1 ) {
                return fail( "countIgnoredDomains must skip an EMPTY tree and still count the rest" );
            }
            // a domain-less tree has no architecture to ask
            if ( MainFrame.countIgnoredDomains( load(
                    "<?xml version=\"1.0\" encoding=\"UTF-8\"?><phyloxml xmlns=\"http://www.phyloxml.org\">"
                            + "<phylogeny rooted=\"true\"><clade><clade><name>X</name></clade>"
                            + "<clade><name>Y</name></clade></clade></phylogeny></phyloxml>",
                    false ) ) != 0 ) {
                return fail( "a tree without domains must report nothing ignored" );
            }
            // ---- the counter is a PARSE-time artifact, not data --------------------------------------------
            // copy() builds a new architecture from the domain list, so a copied tree carries no phantom count
            final DomainArchitecture copied = ( DomainArchitecture ) architectureOf( inverted[ 0 ], "A" ).copy();
            if ( copied.getIgnoredDomainCount() != 0 ) {
                return fail( "a copied architecture must not carry a parse-time ignored count" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "unexpected exception: " + e );
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [MalformedDomainTest] " + msg );
        return false;
    }
}
