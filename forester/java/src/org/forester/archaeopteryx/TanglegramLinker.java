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
import java.util.Arrays;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/**
 * Headless, Swing-free core of the tanglegram feature: given two phylogenies and a link field, it pairs up their
 * external tips by a chosen identity key. Matching is ALL-matches -- a key shared by several tips on either side
 * yields every cross pair (so gene-vs-species and host-vs-parasite many:many links are drawn, not collapsed). It
 * also reports the tips left unmatched on each side, and counts connector crossings for a given tip ordering (a
 * cheap congruence read-out). Kept independent of {@link TanglegramPanel} so it can be unit-tested headless.
 */
final class TanglegramLinker {

    /** The tip attribute two trees are linked on. */
    enum LinkField {
        NODE_NAME( "Node Name" ),
        SCIENTIFIC_NAME( "Taxonomy: Scientific Name" ),
        TAXONOMY_CODE( "Taxonomy: Code" ),
        TAXONOMY_ID( "Taxonomy: NCBI ID" ),
        SEQUENCE_NAME( "Sequence: Name" );

        private final String _label;

        LinkField( final String label ) {
            _label = label;
        }

        String label() {
            return _label;
        }

        /** The link key for a tip, or "" when the tip carries no value for this field (an empty key never matches). */
        String keyFor( final PhylogenyNode node ) {
            switch ( this ) {
                case NODE_NAME:
                    return trimOrEmpty( node.getName() );
                case SCIENTIFIC_NAME:
                    return node.getNodeData().isHasTaxonomy()
                            ? trimOrEmpty( node.getNodeData().getTaxonomy().getScientificName() ) : "";
                case TAXONOMY_CODE:
                    return node.getNodeData().isHasTaxonomy()
                            ? trimOrEmpty( node.getNodeData().getTaxonomy().getTaxonomyCode() ) : "";
                case TAXONOMY_ID:
                    return trimOrEmpty( PhylogenyMethods.getTaxonomyIdentifier( node ) );
                case SEQUENCE_NAME:
                    return node.getNodeData().isHasSequence()
                            ? trimOrEmpty( node.getNodeData().getSequence().getName() ) : "";
                default:
                    return "";
            }
        }

        private static String trimOrEmpty( final String s ) {
            return ForesterUtil.isEmpty( s ) ? "" : s.trim();
        }
    }

    /** One connector: tip {@code a} in the left tree linked to tip {@code b} in the right tree by the shared key. */
    static final class Link {

        private final PhylogenyNode _a;
        private final PhylogenyNode _b;
        private final String        _key;

        Link( final PhylogenyNode a, final PhylogenyNode b, final String key ) {
            _a = a;
            _b = b;
            _key = key;
        }

        PhylogenyNode getA() {
            return _a;
        }

        PhylogenyNode getB() {
            return _b;
        }

        String getKey() {
            return _key;
        }
    }

    /**
     * The links + the tips each side could not match. The two trees may be linked on DIFFERENT fields that hold the
     * same information (e.g. a gene tree's {@code Taxonomy: Scientific Name} joined to a species tree's
     * {@code Node Name}); {@link #getLeftField()}/{@link #getRightField()} are equal in the common single-field case.
     * All lists are unmodifiable.
     */
    static final class Result {

        private final List<Link>          _links;
        private final List<PhylogenyNode> _unmatched_a;
        private final List<PhylogenyNode> _unmatched_b;
        private final LinkField           _left_field;
        private final LinkField           _right_field;

        Result( final List<Link> links, final List<PhylogenyNode> unmatched_a, final List<PhylogenyNode> unmatched_b,
                final LinkField left_field, final LinkField right_field ) {
            _links = Collections.unmodifiableList( links );
            _unmatched_a = Collections.unmodifiableList( unmatched_a );
            _unmatched_b = Collections.unmodifiableList( unmatched_b );
            _left_field = left_field;
            _right_field = right_field;
        }

        List<Link> getLinks() {
            return _links;
        }

        List<PhylogenyNode> getUnmatchedA() {
            return _unmatched_a;
        }

        List<PhylogenyNode> getUnmatchedB() {
            return _unmatched_b;
        }

        /** The field the LEFT tree is linked on. */
        LinkField getLeftField() {
            return _left_field;
        }

        /** The field the RIGHT tree is linked on (equal to {@link #getLeftField()} in the single-field case). */
        LinkField getRightField() {
            return _right_field;
        }
    }

    /** Pair every left tip with every right tip that shares its non-empty {@code field} key (all matches). */
    static Result link( final Phylogeny left, final Phylogeny right, final LinkField field ) {
        return link( left, right, field, field );
    }

    /**
     * Pair every left tip with every right tip that shares its non-empty key, where the left tip's key comes from
     * {@code left_field} and the right tip's from {@code right_field}. This links two trees that store the same
     * information in DIFFERENT fields -- e.g. a gene tree carrying the species in {@code Taxonomy: Scientific Name}
     * joined to a species tree carrying it as {@code Node Name}. When {@code left_field == right_field} this is the
     * ordinary single-field link. The join is still on equal string VALUE, so both fields must hold the same value
     * (this is not an external parasite-&gt;host mapping table).
     */
    static Result link( final Phylogeny left, final Phylogeny right, final LinkField left_field,
                        final LinkField right_field ) {
        final List<PhylogenyNode> left_tips = externalTipsInDisplayOrder( left );
        final List<PhylogenyNode> right_tips = externalTipsInDisplayOrder( right );
        final Map<String, List<PhylogenyNode>> right_by_key = keyIndex( right_tips, right_field );
        final Map<String, List<PhylogenyNode>> left_by_key = keyIndex( left_tips, left_field );
        final List<Link> links = new ArrayList<>();
        for( final PhylogenyNode a : left_tips ) {
            final String key = left_field.keyFor( a );
            if ( key.isEmpty() ) {
                continue;
            }
            final List<PhylogenyNode> matches = right_by_key.get( key );
            if ( matches != null ) {
                for( final PhylogenyNode b : matches ) {
                    links.add( new Link( a, b, key ) );
                }
            }
        }
        final List<PhylogenyNode> unmatched_a = new ArrayList<>();
        for( final PhylogenyNode a : left_tips ) {
            final String key = left_field.keyFor( a );
            if ( key.isEmpty() || !right_by_key.containsKey( key ) ) {
                unmatched_a.add( a );
            }
        }
        final List<PhylogenyNode> unmatched_b = new ArrayList<>();
        for( final PhylogenyNode b : right_tips ) {
            final String key = right_field.keyFor( b );
            if ( key.isEmpty() || !left_by_key.containsKey( key ) ) {
                unmatched_b.add( b );
            }
        }
        return new Result( links, unmatched_a, unmatched_b, left_field, right_field );
    }

    /**
     * Pair tips through an external ASSOCIATION table, for two trees whose linking tips carry DIFFERENT names/values
     * (the classic parasite-vs-host or gene-vs-species-with-different-ids case that a value join cannot handle). A
     * left tip {@code a} is linked to a right tip {@code b} when {@code associations} maps {@code left_field.keyFor(a)}
     * to a key equal to {@code right_field.keyFor(b)}. The mapping is many:many (a left key may list several right
     * keys, and several left keys may share a right key), so every implied cross pair is drawn. A tip is reported
     * unmatched when it takes part in no link. This is the general form of {@link #link}: an identity value join is
     * just the association where every value maps to itself.
     */
    static Result linkByAssociation( final Phylogeny left, final Phylogeny right, final LinkField left_field,
                                     final LinkField right_field, final Map<String, List<String>> associations ) {
        final List<PhylogenyNode> left_tips = externalTipsInDisplayOrder( left );
        final List<PhylogenyNode> right_tips = externalTipsInDisplayOrder( right );
        final Map<String, List<PhylogenyNode>> right_by_key = keyIndex( right_tips, right_field );
        final List<Link> links = new ArrayList<>();
        final Set<PhylogenyNode> matched_a = newIdentitySet();
        final Set<PhylogenyNode> matched_b = newIdentitySet();
        for( final PhylogenyNode a : left_tips ) {
            final String key = left_field.keyFor( a );
            if ( key.isEmpty() ) {
                continue;
            }
            final List<String> right_keys = associations.get( key );
            if ( right_keys == null ) {
                continue;
            }
            for( final String right_key : right_keys ) {
                final List<PhylogenyNode> matches = right_by_key.get( right_key );
                if ( matches != null ) {
                    for( final PhylogenyNode b : matches ) {
                        links.add( new Link( a, b, key ) );
                        matched_a.add( a );
                        matched_b.add( b );
                    }
                }
            }
        }
        final List<PhylogenyNode> unmatched_a = new ArrayList<>();
        for( final PhylogenyNode a : left_tips ) {
            if ( !matched_a.contains( a ) ) {
                unmatched_a.add( a );
            }
        }
        final List<PhylogenyNode> unmatched_b = new ArrayList<>();
        for( final PhylogenyNode b : right_tips ) {
            if ( !matched_b.contains( b ) ) {
                unmatched_b.add( b );
            }
        }
        return new Result( links, unmatched_a, unmatched_b, left_field, right_field );
    }

    private static Set<PhylogenyNode> newIdentitySet() {
        return Collections.newSetFromMap( new IdentityHashMap<PhylogenyNode, Boolean>() );
    }

    private static Map<String, List<PhylogenyNode>> keyIndex( final List<PhylogenyNode> tips, final LinkField field ) {
        final Map<String, List<PhylogenyNode>> index = new LinkedHashMap<>();
        for( final PhylogenyNode tip : tips ) {
            final String key = field.keyFor( tip );
            if ( !key.isEmpty() ) {
                index.computeIfAbsent( key, k -> new ArrayList<>() ).add( tip );
            }
        }
        return index;
    }

    /** External tips top-to-bottom in the tree's current display order (reflects child order / any rotation). */
    static List<PhylogenyNode> externalTipsInDisplayOrder( final Phylogeny phy ) {
        final List<PhylogenyNode> tips = new ArrayList<>();
        if ( ( phy != null ) && !phy.isEmpty() ) {
            for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
                tips.add( it.next() );
            }
        }
        return tips;
    }

    /**
     * Number of connector crossings for the given {left-tip index, right-tip index} pairs: two connectors cross
     * exactly when their left and right endpoints are in the opposite vertical order (a rough congruence measure --
     * 0 == perfectly concordant). Computed as an inversion count in O(n log n): sort by (left asc, right asc) -- so
     * two connectors sharing a left index never count -- then count inversions of the right-index sequence. The
     * auto-untangler evaluates this hundreds of times, so the sub-quadratic cost matters on large trees.
     */
    static int countCrossings( final List<int[]> pairs ) {
        final int n = pairs.size();
        if ( n < 2 ) {
            return 0;
        }
        final int[][] sorted = pairs.toArray( new int[ 0 ][] );
        Arrays.sort( sorted, ( a, b ) -> ( a[ 0 ] != b[ 0 ] ) ? Integer.compare( a[ 0 ], b[ 0 ] )
                : Integer.compare( a[ 1 ], b[ 1 ] ) );
        final int[] right = new int[ n ];
        for( int i = 0; i < n; i++ ) {
            right[ i ] = sorted[ i ][ 1 ];
        }
        final long inversions = countInversions( right, new int[ n ], 0, n - 1 );
        return (int) Math.min( inversions, Integer.MAX_VALUE );
    }

    /**
     * A size-normalised "entanglement" score in [0,1] for the current tip ordering: the connector-crossing (inversion)
     * count divided by the maximum possible for that many connectors ({@code n*(n-1)/2}). 0 == perfectly concordant
     * (fully untangled), 1 == maximally discordant. Unlike the raw crossing count this is comparable across
     * tanglegrams of different sizes -- the standard way to report topological (dis)agreement. Returns 0 for fewer
     * than two connectors (nothing can cross). The crossing count is an inversion count, which never exceeds
     * {@code n*(n-1)/2}, so the result is always in range (the clamp is defensive).
     */
    static double entanglement( final int crossings, final int connector_count ) {
        if ( connector_count < 2 ) {
            return 0.0;
        }
        final long max = ( (long) connector_count * ( connector_count - 1 ) ) / 2L;
        return ( max <= 0 ) ? 0.0 : Math.min( 1.0, crossings / (double) max );
    }

    /** Inversions in a[lo..hi] (pairs i&lt;j with a[i] &gt; a[j]) via a counting merge sort; sorts a[lo..hi] in place. */
    private static long countInversions( final int[] a, final int[] tmp, final int lo, final int hi ) {
        if ( lo >= hi ) {
            return 0;
        }
        final int mid = ( lo + hi ) >>> 1;
        long inv = countInversions( a, tmp, lo, mid ) + countInversions( a, tmp, mid + 1, hi );
        int i = lo;
        int j = mid + 1;
        int k = lo;
        while ( ( i <= mid ) && ( j <= hi ) ) {
            if ( a[ i ] <= a[ j ] ) {
                tmp[ k++ ] = a[ i++ ];
            }
            else {
                tmp[ k++ ] = a[ j++ ];
                inv += ( mid - i + 1 ); // a[i..mid] are all > a[j]
            }
        }
        while ( i <= mid ) {
            tmp[ k++ ] = a[ i++ ];
        }
        while ( j <= hi ) {
            tmp[ k++ ] = a[ j++ ];
        }
        System.arraycopy( tmp, lo, a, lo, ( hi - lo ) + 1 );
        return inv;
    }

    /**
     * For each connector (in the input order), whether it crosses at least one other connector. Two connectors cross
     * when their left and right endpoints are in the opposite vertical order (the same test as {@link #countCrossings}).
     * O(n^2), computed on demand (not per frame) -- used only to colour the crossing connectors.
     */
    static boolean[] crossingConnectors( final List<int[]> pairs ) {
        final int n = pairs.size();
        final boolean[] crossing = new boolean[ n ];
        for( int i = 0; i < n; i++ ) {
            final int[] p = pairs.get( i );
            for( int j = i + 1; j < n; j++ ) {
                final int[] q = pairs.get( j );
                // long product: index differences can each exceed 46k on a huge tanglegram, overflowing an int
                if ( ( (long) ( p[ 0 ] - q[ 0 ] ) * ( p[ 1 ] - q[ 1 ] ) ) < 0 ) {
                    crossing[ i ] = true;
                    crossing[ j ] = true;
                }
            }
        }
        return crossing;
    }

    private TanglegramLinker() {
    }
}
