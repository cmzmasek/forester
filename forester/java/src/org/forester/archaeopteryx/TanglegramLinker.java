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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

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

    /** The links + the tips each side could not match, for one field. All lists are unmodifiable. */
    static final class Result {

        private final List<Link>          _links;
        private final List<PhylogenyNode> _unmatched_a;
        private final List<PhylogenyNode> _unmatched_b;
        private final LinkField           _field;

        Result( final List<Link> links, final List<PhylogenyNode> unmatched_a, final List<PhylogenyNode> unmatched_b,
                final LinkField field ) {
            _links = Collections.unmodifiableList( links );
            _unmatched_a = Collections.unmodifiableList( unmatched_a );
            _unmatched_b = Collections.unmodifiableList( unmatched_b );
            _field = field;
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

        LinkField getField() {
            return _field;
        }
    }

    /** Pair every left tip with every right tip that shares its non-empty {@code field} key (all matches). */
    static Result link( final Phylogeny left, final Phylogeny right, final LinkField field ) {
        final List<PhylogenyNode> left_tips = externalTipsInDisplayOrder( left );
        final List<PhylogenyNode> right_tips = externalTipsInDisplayOrder( right );
        final Map<String, List<PhylogenyNode>> right_by_key = keyIndex( right_tips, field );
        final Map<String, List<PhylogenyNode>> left_by_key = keyIndex( left_tips, field );
        final List<Link> links = new ArrayList<>();
        for( final PhylogenyNode a : left_tips ) {
            final String key = field.keyFor( a );
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
            final String key = field.keyFor( a );
            if ( key.isEmpty() || !right_by_key.containsKey( key ) ) {
                unmatched_a.add( a );
            }
        }
        final List<PhylogenyNode> unmatched_b = new ArrayList<>();
        for( final PhylogenyNode b : right_tips ) {
            final String key = field.keyFor( b );
            if ( key.isEmpty() || !left_by_key.containsKey( key ) ) {
                unmatched_b.add( b );
            }
        }
        return new Result( links, unmatched_a, unmatched_b, field );
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

    private TanglegramLinker() {
    }
}
