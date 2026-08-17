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
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.forester.archaeopteryx.TanglegramLinker.Link;
import org.forester.archaeopteryx.TanglegramLinker.LinkField;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Headless tests for {@link TanglegramUntangler}: it lowers the crossing count (to 0 when the two trees have the same
 * unordered topology, i.e. are reachable from one another by flips), leaves an already-untangled pair alone, and the
 * net flip set it returns is a true inverse (reversing it restores the original ordering).
 */
public final class TanglegramUntanglerTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TanglegramUntangler: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return tangledToZeroOk() && alreadyCleanOk() && eightTipUntangleOk() && undoRestoresOk() && reachesOptimumOk()
                && largeTreeUntangleOk();
    }

    /** Performance / scale: a large (300-tip) heavily-scrambled SAME-topology pair (so a flip-reachable 0 exists) must
     *  untangle to 0 crossings, quickly. Guards the O(n log n) crossing count + barycentre passes against a
     *  re-quadratic regression -- the untangler evaluates the crossing count hundreds of times, so this is where a
     *  regression would bite. Deterministic (fixed seeds + the untangler's own fixed seed). */
    private static boolean largeTreeUntangleOk() {
        final Random rng = new Random( 99 );
        final int n = 300;
        final Phylogeny left = tree( randomBinary( n, rng ) );
        final Phylogeny right = left.copy(); // same unordered topology -> untangle can reach exactly 0 by flips
        for( final PhylogenyNode node : internalNodes( right ) ) {
            if ( rng.nextBoolean() ) {
                TanglegramUntangler.reverse( node );
            }
        }
        final List<Link> links = links( left, right );
        if ( links.size() != n ) {
            return fail( "expected " + n + " 1:1 links on the large pair, got " + links.size() );
        }
        final int before = crossings( left, right, links );
        if ( before < 100 ) {
            return fail( "the scrambled " + n + "-tip pair should start heavily tangled, got " + before );
        }
        final long t0 = System.nanoTime();
        TanglegramUntangler.untangle( left, right, links );
        final long ms = ( System.nanoTime() - t0 ) / 1_000_000L;
        final int after = crossings( left, right, links );
        if ( after != 0 ) {
            return fail( "untangle should reach 0 on the flip-reachable " + n + "-tip pair, got " + after + " (was "
                    + before + ")" );
        }
        // very generous bound: this runs in a few ms; it only trips on a catastrophic (e.g. O(n^2)) regression
        if ( ms > 10_000 ) {
            return fail( "untangle took " + ms + " ms on " + n + " tips -- likely a performance regression" );
        }
        // the HARDER, realistic case: two INDEPENDENT large topologies (the random restarts + net-parity flips
        // genuinely engage here, so this is the true performance path). Untangle must never worsen, must reduce, and
        // stay fast even with all the restarts running.
        final Random rng2 = new Random( 7 );
        final Phylogeny l2 = tree( randomBinary( n, rng2 ) );
        final Phylogeny r2 = tree( randomBinary( n, rng2 ) );
        final List<Link> links2 = links( l2, r2 );
        final int before2 = crossings( l2, r2, links2 );
        final long t2 = System.nanoTime();
        TanglegramUntangler.untangle( l2, r2, links2 );
        final long ms2 = ( System.nanoTime() - t2 ) / 1_000_000L;
        final int after2 = crossings( l2, r2, links2 );
        if ( after2 >= before2 ) { // never-worse guarantees after2 <= before2, so this means "did not reduce"
            return fail( "untangle should reduce (and never worsen) the " + n + "-tip different-topology pair, before="
                    + before2 + " after=" + after2 );
        }
        if ( ms2 > 10_000 ) {
            return fail( "untangle (with restarts) took " + ms2 + " ms on " + n + " different-topology tips -- "
                    + "likely a performance regression" );
        }
        return true;
    }

    private static boolean reachesOptimumOk() {
        // scrambled random pairs where barycentre-from-the-original is often SUBOPTIMAL, so the random restarts (and
        // the net-parity flip toggle they exercise) must engage; untangle must still hit the brute-force minimum
        final Random rng = new Random( 7 );
        for( int trial = 0; trial < 20; trial++ ) {
            final Phylogeny left = tree( randomBinary( 7, rng ) );
            final Phylogeny right = tree( randomBinary( 7, rng ) );
            final List<Link> links = links( left, right );
            final int optimum = bruteForceMin( left, right, links );
            TanglegramUntangler.untangle( left, right, links );
            final int result = crossings( left, right, links );
            if ( result != optimum ) {
                return fail( "trial " + trial + ": untangle reached " + result + " crossings, brute-force optimum is "
                        + optimum );
            }
        }
        return true;
    }

    private static PhylogenyNode randomBinary( final int n, final Random rng ) {
        final List<PhylogenyNode> nodes = new ArrayList<>();
        for( int i = 0; i < n; i++ ) {
            nodes.add( leaf( "t" + i ) );
        }
        while ( nodes.size() > 1 ) {
            final PhylogenyNode a = nodes.remove( rng.nextInt( nodes.size() ) );
            final PhylogenyNode b = nodes.remove( rng.nextInt( nodes.size() ) );
            nodes.add( clade( a, b ) );
        }
        return nodes.get( 0 );
    }

    /** Exhaustively tries every flip combination of both trees' internal nodes and returns the fewest crossings. */
    private static int bruteForceMin( final Phylogeny left, final Phylogeny right, final List<Link> links ) {
        final List<PhylogenyNode> li = internalNodes( left );
        final List<PhylogenyNode> ri = internalNodes( right );
        int min = Integer.MAX_VALUE;
        for( int lm = 0; lm < ( 1 << li.size() ); lm++ ) {
            applyMask( li, lm );
            for( int rm = 0; rm < ( 1 << ri.size() ); rm++ ) {
                applyMask( ri, rm );
                min = Math.min( min, crossings( left, right, links ) );
                applyMask( ri, rm ); // revert (reverse is involutive)
            }
            applyMask( li, lm ); // revert
        }
        return min;
    }

    private static void applyMask( final List<PhylogenyNode> internals, final int mask ) {
        for( int i = 0; i < internals.size(); i++ ) {
            if ( ( ( mask >> i ) & 1 ) == 1 ) {
                TanglegramUntangler.reverse( internals.get( i ) );
            }
        }
    }

    private static List<PhylogenyNode> internalNodes( final Phylogeny phy ) {
        final List<PhylogenyNode> out = new ArrayList<>();
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() ) {
                out.add( n );
            }
        }
        return out;
    }

    private static boolean tangledToZeroOk() {
        final Phylogeny left = tree( clade( clade( leaf( "A" ), leaf( "B" ) ), clade( leaf( "C" ), leaf( "D" ) ) ) );
        final Phylogeny right = tree( clade( clade( leaf( "C" ), leaf( "D" ) ), clade( leaf( "A" ), leaf( "B" ) ) ) );
        final List<Link> links = links( left, right );
        if ( crossings( left, right, links ) != 4 ) {
            return fail( "the tangled 4-tip pair should start at 4 crossings" );
        }
        final List<PhylogenyNode> flips = TanglegramUntangler.untangle( left, right, links );
        if ( crossings( left, right, links ) != 0 ) {
            return fail( "untangle should reach 0 crossings on a flip-reachable pair, got "
                    + crossings( left, right, links ) );
        }
        if ( flips.isEmpty() ) {
            return fail( "untangle changed the layout but reported no flips" );
        }
        return true;
    }

    private static boolean alreadyCleanOk() {
        final Phylogeny left = tree( clade( clade( leaf( "A" ), leaf( "B" ) ), clade( leaf( "C" ), leaf( "D" ) ) ) );
        final Phylogeny right = tree( clade( clade( leaf( "A" ), leaf( "B" ) ), clade( leaf( "C" ), leaf( "D" ) ) ) );
        final List<Link> links = links( left, right );
        final List<PhylogenyNode> flips = TanglegramUntangler.untangle( left, right, links );
        if ( !flips.isEmpty() ) {
            return fail( "an already-untangled pair should need no flips, got " + flips.size() );
        }
        if ( crossings( left, right, links ) != 0 ) {
            return fail( "an already-untangled pair should stay at 0 crossings" );
        }
        return true;
    }

    private static boolean eightTipUntangleOk() {
        // right has the SAME unordered topology as left (every cherry + some pairs flipped), so untangle can reach 0
        final Phylogeny left = tree( clade( clade( clade( leaf( "A" ), leaf( "B" ) ), clade( leaf( "C" ), leaf( "D" ) ) ),
                                            clade( clade( leaf( "E" ), leaf( "F" ) ),
                                                   clade( leaf( "G" ), leaf( "H" ) ) ) ) );
        final Phylogeny right = tree( clade( clade( clade( leaf( "D" ), leaf( "C" ) ), clade( leaf( "B" ), leaf( "A" ) ) ),
                                             clade( clade( leaf( "H" ), leaf( "G" ) ),
                                                    clade( leaf( "F" ), leaf( "E" ) ) ) ) );
        final List<Link> links = links( left, right );
        final int before = crossings( left, right, links );
        if ( before == 0 ) {
            return fail( "the scrambled 8-tip pair should start with crossings" );
        }
        TanglegramUntangler.untangle( left, right, links );
        if ( crossings( left, right, links ) != 0 ) {
            return fail( "untangle should reach 0 on the flip-reachable 8-tip pair, got "
                    + crossings( left, right, links ) + " (was " + before + ")" );
        }
        return true;
    }

    private static boolean undoRestoresOk() {
        final Phylogeny left = tree( clade( clade( leaf( "A" ), leaf( "B" ) ), clade( leaf( "C" ), leaf( "D" ) ) ) );
        final Phylogeny right = tree( clade( clade( leaf( "C" ), leaf( "D" ) ), clade( leaf( "A" ), leaf( "B" ) ) ) );
        final List<Link> links = links( left, right );
        final int before = crossings( left, right, links );
        final List<PhylogenyNode> flips = TanglegramUntangler.untangle( left, right, links );
        // reversing every node in the returned action must restore the original crossing count (undo correctness)
        for( final PhylogenyNode node : flips ) {
            TanglegramUntangler.reverse( node );
        }
        if ( crossings( left, right, links ) != before ) {
            return fail( "reversing the flip set should restore the original " + before + " crossings, got "
                    + crossings( left, right, links ) );
        }
        return true;
    }

    // ---- helpers -----------------------------------------------------------------------------------------------

    private static List<Link> links( final Phylogeny left, final Phylogeny right ) {
        return TanglegramLinker.link( left, right, LinkField.NODE_NAME ).getLinks();
    }

    private static int crossings( final Phylogeny left, final Phylogeny right, final List<Link> links ) {
        final Map<PhylogenyNode, Integer> li = index( left );
        final Map<PhylogenyNode, Integer> ri = index( right );
        final List<int[]> pairs = new ArrayList<>();
        for( final Link link : links ) {
            final Integer a = li.get( link.getA() );
            final Integer b = ri.get( link.getB() );
            if ( ( a != null ) && ( b != null ) ) {
                pairs.add( new int[] { a, b } );
            }
        }
        return TanglegramLinker.countCrossings( pairs );
    }

    private static Map<PhylogenyNode, Integer> index( final Phylogeny phy ) {
        final List<PhylogenyNode> tips = TanglegramLinker.externalTipsInDisplayOrder( phy );
        final Map<PhylogenyNode, Integer> m = new IdentityHashMap<>();
        for( int i = 0; i < tips.size(); i++ ) {
            m.put( tips.get( i ), i );
        }
        return m;
    }

    private static boolean fail( final String message ) {
        System.out.println( "TanglegramUntangler test failed: " + message );
        return false;
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
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
