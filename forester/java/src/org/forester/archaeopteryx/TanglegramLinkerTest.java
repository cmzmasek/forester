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
import java.util.List;
import java.util.Random;

import org.forester.archaeopteryx.TanglegramLinker.LinkField;
import org.forester.archaeopteryx.TanglegramLinker.Result;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headless tests for {@link TanglegramLinker}: node-name / scientific-name keying, per-tree link fields (linking two
 * trees on different fields holding the same value), all-matches (many:many) linking, unmatched-tip accounting, and
 * the crossing count.
 */
public final class TanglegramLinkerTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TanglegramLinker: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return nodeNameMatchingOk() && manyToManyOk() && unmatchedOk() && scientificNameOk() && crossFieldLinkOk()
                && crossingsOk() && crossingsMatchBruteForceOk() && crossingConnectorsOk();
    }

    /** Two trees that store the same value in DIFFERENT fields link when a field is chosen per tree (e.g. a gene
     *  tree's taxonomy scientific name joined to a species tree's node name), but not with one shared field. */
    private static boolean crossFieldLinkOk() {
        // gene tree: node name = accession, the species lives in the taxonomy scientific name
        final Phylogeny gene = tree( clade( taxonLeaf( "acc1", "Felis catus" ), taxonLeaf( "acc2", "Canis lupus" ) ) );
        // species tree: the species IS the node name
        final Phylogeny species = tree( clade( leaf( "Canis lupus" ), leaf( "Felis catus" ) ) );
        final Result r = TanglegramLinker.link( gene, species, LinkField.SCIENTIFIC_NAME, LinkField.NODE_NAME );
        if ( ( r.getLinks().size() != 2 ) || !r.getUnmatchedA().isEmpty() || !r.getUnmatchedB().isEmpty() ) {
            return fail( "cross-field link (sci-name <-> node-name) should match both species with none unmatched; got "
                    + r.getLinks().size() + " links, " + r.getUnmatchedA().size() + "/" + r.getUnmatchedB().size()
                    + " unmatched" );
        }
        if ( ( r.getLeftField() != LinkField.SCIENTIFIC_NAME ) || ( r.getRightField() != LinkField.NODE_NAME ) ) {
            return fail( "the result should carry both link fields" );
        }
        // the same field on both sides must NOT link these trees -- which is exactly why per-tree fields are needed:
        // node name (accession vs species) and scientific name (the species tree has no taxonomy) each match nothing
        if ( !TanglegramLinker.link( gene, species, LinkField.NODE_NAME ).getLinks().isEmpty()
                || !TanglegramLinker.link( gene, species, LinkField.SCIENTIFIC_NAME ).getLinks().isEmpty() ) {
            return fail( "a single shared field must not link the accession/species trees" );
        }
        // the single-field overload reports the same field on both sides
        final Result same = TanglegramLinker.link( gene, gene, LinkField.NODE_NAME );
        if ( same.getLeftField() != same.getRightField() ) {
            return fail( "single-field link should report equal left/right fields" );
        }
        return true;
    }

    private static boolean crossingConnectorsOk() {
        // {0,0} and {3,3} cross nobody; {1,2} and {2,1} cross each other
        final List<int[]> pairs = new ArrayList<>();
        pairs.add( new int[] { 0, 0 } );
        pairs.add( new int[] { 1, 2 } );
        pairs.add( new int[] { 2, 1 } );
        pairs.add( new int[] { 3, 3 } );
        final boolean[] flags = TanglegramLinker.crossingConnectors( pairs );
        final boolean[] expected = { false, true, true, false };
        for( int i = 0; i < expected.length; i++ ) {
            if ( flags[ i ] != expected[ i ] ) {
                return fail( "crossingConnectors[" + i + "] = " + flags[ i ] + ", expected " + expected[ i ] );
            }
        }
        return true;
    }

    /** The O(n log n) inversion-count crossings must equal the O(n^2) definition on random pair-lists, including
     *  duplicate indices (a tip linked to several partners -- shared indices must NOT count as a crossing). */
    private static boolean crossingsMatchBruteForceOk() {
        final Random rng = new Random( 3 );
        for( int trial = 0; trial < 60; trial++ ) {
            final int n = rng.nextInt( 40 );
            final List<int[]> pairs = new ArrayList<>();
            for( int i = 0; i < n; i++ ) {
                pairs.add( new int[] { rng.nextInt( 12 ), rng.nextInt( 12 ) } ); // small range -> forces duplicates
            }
            int brute = 0;
            for( int i = 0; i < pairs.size(); i++ ) {
                for( int j = i + 1; j < pairs.size(); j++ ) {
                    final int[] p = pairs.get( i );
                    final int[] q = pairs.get( j );
                    if ( ( ( p[ 0 ] - q[ 0 ] ) * ( p[ 1 ] - q[ 1 ] ) ) < 0 ) {
                        brute++;
                    }
                }
            }
            final int fast = TanglegramLinker.countCrossings( pairs );
            if ( fast != brute ) {
                return fail( "trial " + trial + " (n=" + n + "): inversion count " + fast + " != brute-force " + brute );
            }
        }
        return true;
    }

    private static boolean nodeNameMatchingOk() {
        // A over (A,B,C); B over (C,B,A) -- same names, reversed order -> 3 links, none unmatched, fully crossing
        final Phylogeny a = tree( clade( leaf( "A" ), leaf( "B" ), leaf( "C" ) ) );
        final Phylogeny b = tree( clade( leaf( "C" ), leaf( "B" ), leaf( "A" ) ) );
        final Result r = TanglegramLinker.link( a, b, LinkField.NODE_NAME );
        if ( r.getLinks().size() != 3 ) {
            return fail( "expected 3 node-name links, got " + r.getLinks().size() );
        }
        if ( !r.getUnmatchedA().isEmpty() || !r.getUnmatchedB().isEmpty() ) {
            return fail( "expected no unmatched tips for the reversed same-name pair" );
        }
        return true;
    }

    private static boolean manyToManyOk() {
        // two "X" tips on each side -> ALL four cross pairs (draw-all-matches), plus a clean 1:1 for "Y"
        final Phylogeny a = tree( clade( leaf( "X" ), leaf( "X" ), leaf( "Y" ) ) );
        final Phylogeny b = tree( clade( leaf( "X" ), leaf( "X" ), leaf( "Y" ) ) );
        final Result r = TanglegramLinker.link( a, b, LinkField.NODE_NAME );
        if ( r.getLinks().size() != 5 ) {
            return fail( "expected 5 links (2x2 for X + 1 for Y), got " + r.getLinks().size() );
        }
        int x_links = 0;
        for( final TanglegramLinker.Link link : r.getLinks() ) {
            if ( "X".equals( link.getKey() ) ) {
                ++x_links;
            }
        }
        if ( x_links != 4 ) {
            return fail( "expected 4 X links (all matches), got " + x_links );
        }
        return true;
    }

    private static boolean unmatchedOk() {
        // A has Q + an unnamed tip; B has Q + Z -- Q links, the unnamed A tip and Z (B) are unmatched
        final Phylogeny a = tree( clade( leaf( "Q" ), leaf( "" ) ) );
        final Phylogeny b = tree( clade( leaf( "Q" ), leaf( "Z" ) ) );
        final Result r = TanglegramLinker.link( a, b, LinkField.NODE_NAME );
        if ( r.getLinks().size() != 1 ) {
            return fail( "expected 1 link (Q), got " + r.getLinks().size() );
        }
        if ( ( r.getUnmatchedA().size() != 1 ) || ( r.getUnmatchedB().size() != 1 ) ) {
            return fail( "expected 1 unmatched on each side, got a=" + r.getUnmatchedA().size() + " b="
                    + r.getUnmatchedB().size() );
        }
        return true;
    }

    private static boolean scientificNameOk() {
        // node names differ, scientific names match -> links on SCIENTIFIC_NAME but NOT on NODE_NAME
        final Phylogeny a = tree( clade( taxonLeaf( "a1", "Felis catus" ), taxonLeaf( "a2", "Canis lupus" ) ) );
        final Phylogeny b = tree( clade( taxonLeaf( "b1", "Canis lupus" ), taxonLeaf( "b2", "Felis catus" ) ) );
        final Result by_sci = TanglegramLinker.link( a, b, LinkField.SCIENTIFIC_NAME );
        if ( ( by_sci.getLinks().size() != 2 ) || !by_sci.getUnmatchedA().isEmpty() ) {
            return fail( "expected 2 scientific-name links, none unmatched; got " + by_sci.getLinks().size() );
        }
        final Result by_name = TanglegramLinker.link( a, b, LinkField.NODE_NAME );
        if ( !by_name.getLinks().isEmpty() ) {
            return fail( "the distinct node names must NOT link, got " + by_name.getLinks().size() );
        }
        return true;
    }

    private static boolean crossingsOk() {
        final List<int[]> parallel = new ArrayList<>();
        parallel.add( new int[] { 0, 0 } );
        parallel.add( new int[] { 1, 1 } );
        parallel.add( new int[] { 2, 2 } );
        if ( TanglegramLinker.countCrossings( parallel ) != 0 ) {
            return fail( "parallel connectors should have 0 crossings" );
        }
        final List<int[]> crossed = new ArrayList<>();
        crossed.add( new int[] { 0, 1 } );
        crossed.add( new int[] { 1, 0 } );
        if ( TanglegramLinker.countCrossings( crossed ) != 1 ) {
            return fail( "one swapped pair should have 1 crossing" );
        }
        final List<int[]> reversed = new ArrayList<>();
        for( int i = 0; i < 4; i++ ) {
            reversed.add( new int[] { i, 3 - i } );
        }
        if ( TanglegramLinker.countCrossings( reversed ) != 6 ) {
            return fail( "a fully reversed 4-tip mapping should have 6 crossings, got "
                    + TanglegramLinker.countCrossings( reversed ) );
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "TanglegramLinker test failed: " + message );
        return false;
    }

    // ---- fixtures ----------------------------------------------------------------------------------------------

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode taxonLeaf( final String name, final String scientific_name ) {
        final PhylogenyNode n = leaf( name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode clade( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode child : children ) {
            n.addAsChild( child );
        }
        return n;
    }

    private static Phylogeny tree( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }
}
