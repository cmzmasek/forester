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
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Headless crossing-minimisation for a tanglegram. Reorders BOTH trees using flip-only moves (reverse a clade's child
 * order -- the same topology-preserving operation as a manual rotation), so the result is reachable by, and undoable
 * as, ordinary flips.
 *
 * <p>Heuristic: alternating <b>barycentre</b> passes (order each clade's children so the mean position of the tips
 * they link to in the OTHER tree ascends -- top = smallest -- iterated until the crossing count stops dropping), run
 * from the current orientation AND from several random-restart orientations, keeping whichever layout has the fewest
 * crossings. The random restarts escape the local minima a single barycentre pass gets stuck in (empirically it then
 * reaches the brute-force optimum on all small pairs). Optimal untangling is NP-hard, so this is a strong heuristic,
 * not a guaranteed minimum; it never leaves the tanglegram worse than it found it.
 *
 * <p>Returns the NET set of flipped nodes (already applied to the trees), so the caller records one undoable action.
 */
final class TanglegramUntangler {

    private static final int  MAX_ITERATIONS = 8;
    // barycentre alone sticks in local minima on tangled pairs (escaping needs coordinated multi-node flips), so
    // restart from random flip states and keep the best layout. Fixed seed => deterministic / reproducible.
    private static final int  MAX_RESTARTS   = 40;
    private static final long RANDOM_SEED    = 42L;

    /** Flip a node's child order in place -- the tanglegram's rotate primitive. Its own inverse (reversing twice
     *  restores the original order), which is what makes undo/redo of both manual rotation and auto-untangle trivial.
     *  Works for any node arity (unlike {@code PhylogenyNode.swapChildren}, which requires exactly two children). */
    static void reverse( final PhylogenyNode node ) {
        final List<PhylogenyNode> children = node.getDescendants();
        final int n = children.size();
        final PhylogenyNode[] snapshot = children.toArray( new PhylogenyNode[ 0 ] );
        for( int i = 0; i < n; i++ ) {
            node.setChildNode( i, snapshot[ n - 1 - i ] );
        }
    }

    /** Reorders both trees (flips only) to reduce connector crossings; returns the net set of flipped nodes. Runs the
     *  barycentre local optimiser from the current orientation AND from several random-restart orientations, keeping
     *  whichever gives the fewest crossings. Never leaves the tanglegram worse than it found it (the best-kept layout
     *  starts as the original). */
    static List<PhylogenyNode> untangle( final Phylogeny left, final Phylogeny right,
                                         final List<TanglegramLinker.Link> links ) {
        final int initial = crossings( left, right, links );
        if ( initial == 0 ) {
            return new ArrayList<>();
        }
        final List<PhylogenyNode> internals = new ArrayList<>();
        collectInternals( left, internals );
        collectInternals( right, internals );
        final Set<PhylogenyNode> flipped = newFlipSet();
        Set<PhylogenyNode> best = newFlipSet(); // empty == the original layout, so the result can never be worse
        int best_crossings = initial;
        final Random rng = new Random( RANDOM_SEED );
        for( int restart = 0; ( restart <= MAX_RESTARTS ) && ( best_crossings > 0 ); restart++ ) {
            resetToOriginal( flipped );
            if ( restart > 0 ) {
                for( final PhylogenyNode node : internals ) {
                    if ( rng.nextBoolean() ) {
                        flip( node, flipped );
                    }
                }
            }
            localOptimize( left, right, links, flipped );
            final int c = crossings( left, right, links );
            if ( c < best_crossings ) {
                best_crossings = c;
                best = copyOf( flipped );
            }
        }
        resetToOriginal( flipped );
        for( final PhylogenyNode node : best ) {
            flip( node, flipped );
        }
        return new ArrayList<>( flipped );
    }

    /** Barycentre iteration (both sides) from the current orientation, returning the BEST layout it reaches -- a
     *  double-pass is not monotonic in crossings, so we keep the best seen rather than the (possibly overshot) last. */
    private static void localOptimize( final Phylogeny left, final Phylogeny right,
                                       final List<TanglegramLinker.Link> links, final Set<PhylogenyNode> flipped ) {
        int prev = crossings( left, right, links );
        int best = prev;
        Set<PhylogenyNode> best_state = copyOf( flipped );
        for( int iter = 0; ( iter < MAX_ITERATIONS ) && ( prev > 0 ); iter++ ) {
            barycenterPass( left, right, links, true, flipped );
            barycenterPass( right, left, links, false, flipped );
            final int now = crossings( left, right, links );
            if ( now < best ) {
                best = now;
                best_state = copyOf( flipped );
            }
            if ( now >= prev ) {
                break;
            }
            prev = now;
        }
        restoreTo( flipped, best_state );
    }

    /** Moves the trees from the current flip state to {@code target} by reversing the symmetric difference of the two
     *  flip sets (each such node's parity must change), leaving {@code flipped} equal to {@code target}. */
    private static void restoreTo( final Set<PhylogenyNode> flipped, final Set<PhylogenyNode> target ) {
        for( final PhylogenyNode node : flipped ) {
            if ( !target.contains( node ) ) {
                reverse( node );
            }
        }
        for( final PhylogenyNode node : target ) {
            if ( !flipped.contains( node ) ) {
                reverse( node );
            }
        }
        flipped.clear();
        flipped.addAll( target );
    }

    private static void resetToOriginal( final Set<PhylogenyNode> flipped ) {
        for( final PhylogenyNode node : flipped ) {
            reverse( node );
        }
        flipped.clear();
    }

    private static Set<PhylogenyNode> newFlipSet() {
        return Collections.newSetFromMap( new IdentityHashMap<PhylogenyNode, Boolean>() );
    }

    private static Set<PhylogenyNode> copyOf( final Set<PhylogenyNode> set ) {
        final Set<PhylogenyNode> copy = newFlipSet();
        copy.addAll( set );
        return copy;
    }

    private static void collectInternals( final Phylogeny phy, final List<PhylogenyNode> out ) {
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( !node.isExternal() ) {
                out.add( node );
            }
        }
    }

    /** Reverses a node and toggles its membership in {@code flipped}, so the set always holds the NET parity of flips
     *  (a node flipped an even number of times is not in the set, and its subtree is back to the original order). */
    private static void flip( final PhylogenyNode node, final Set<PhylogenyNode> flipped ) {
        reverse( node );
        if ( !flipped.add( node ) ) {
            flipped.remove( node );
        }
    }

    /** One barycentre pass over {@code tree} (holding {@code other} fixed): postorder, ordering each clade's children
     *  so their subtree barycentres ascend, achieved by flipping a node when its first child's mean target exceeds
     *  its last child's (for a binary node this is exactly the crossing-reducing choice). */
    private static void barycenterPass( final Phylogeny tree, final Phylogeny other,
                                        final List<TanglegramLinker.Link> links, final boolean tree_is_left,
                                        final Set<PhylogenyNode> flipped ) {
        final Map<PhylogenyNode, Integer> other_index = indexOf( other );
        final Map<PhylogenyNode, double[]> tip_target = new IdentityHashMap<>(); // tip -> {sum, count} of partner indices
        for( final TanglegramLinker.Link link : links ) {
            final PhylogenyNode tip = tree_is_left ? link.getA() : link.getB();
            final PhylogenyNode partner = tree_is_left ? link.getB() : link.getA();
            final Integer idx = other_index.get( partner );
            if ( idx != null ) {
                final double[] acc = tip_target.computeIfAbsent( tip, k -> new double[ 2 ] );
                acc[ 0 ] += idx;
                acc[ 1 ] += 1;
            }
        }
        subtreeMeanAndFlip( tree.getRoot(), tip_target, flipped );
    }

    /** Postorder: returns {sum,count} of tip targets in the subtree, and flips this node if its first child's mean
     *  target is greater than its last child's (a flip-only barycentre ordering). Each child is visited exactly once.
     *  For a polytomy (&gt;2 children) the whole child list is only reversed-or-not -- a flip can't reorder a node's
     *  interior -- so a multifurcating clade's middle children keep their relative order, exactly as a manual rotation
     *  would leave them; binary nodes (the common case) are ordered exactly. */
    private static double[] subtreeMeanAndFlip( final PhylogenyNode node, final Map<PhylogenyNode, double[]> tip_target,
                                                final Set<PhylogenyNode> flipped ) {
        if ( node.isExternal() ) {
            final double[] t = tip_target.get( node );
            return ( t == null ) ? new double[] { 0, 0 } : new double[] { t[ 0 ], t[ 1 ] };
        }
        final List<PhylogenyNode> children = node.getDescendants();
        final int nc = children.size();
        final double[][] child_sc = new double[ nc ][];
        double total_sum = 0;
        double total_count = 0;
        for( int i = 0; i < nc; i++ ) {
            child_sc[ i ] = subtreeMeanAndFlip( children.get( i ), tip_target, flipped );
            total_sum += child_sc[ i ][ 0 ];
            total_count += child_sc[ i ][ 1 ];
        }
        final Double first_mean = mean( child_sc[ 0 ] );
        final Double last_mean = mean( child_sc[ nc - 1 ] );
        if ( ( first_mean != null ) && ( last_mean != null ) && ( first_mean > last_mean ) ) {
            flip( node, flipped );
        }
        return new double[] { total_sum, total_count };
    }

    private static Double mean( final double[] sum_count ) {
        return ( sum_count[ 1 ] > 0 ) ? ( sum_count[ 0 ] / sum_count[ 1 ] ) : null;
    }

    private static int crossings( final Phylogeny left, final Phylogeny right,
                                  final List<TanglegramLinker.Link> links ) {
        final Map<PhylogenyNode, Integer> left_index = indexOf( left );
        final Map<PhylogenyNode, Integer> right_index = indexOf( right );
        final List<int[]> pairs = new ArrayList<>();
        for( final TanglegramLinker.Link link : links ) {
            final Integer a = left_index.get( link.getA() );
            final Integer b = right_index.get( link.getB() );
            if ( ( a != null ) && ( b != null ) ) {
                pairs.add( new int[] { a, b } );
            }
        }
        return TanglegramLinker.countCrossings( pairs );
    }

    private static Map<PhylogenyNode, Integer> indexOf( final Phylogeny phy ) {
        final List<PhylogenyNode> tips = TanglegramLinker.externalTipsInDisplayOrder( phy );
        final Map<PhylogenyNode, Integer> index = new IdentityHashMap<>();
        for( int i = 0; i < tips.size(); i++ ) {
            index.put( tips.get( i ), i );
        }
        return index;
    }

    private TanglegramUntangler() {
    }
}
