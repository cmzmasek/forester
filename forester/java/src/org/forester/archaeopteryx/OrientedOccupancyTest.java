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

/**
 * {@link OrientedOccupancy}: the separating-axis test itself (the case that motivates the class -- two rotated boxes
 * whose axis-aligned bounds overlap while the boxes do not), first-come placing (ask, then place -- the two calls
 * production makes), placing without asking, the cell grid, and the non-behaviours (a reset forgets; a read records
 * nothing; a zero-size box reserves nothing). Since 0.11.181 also the COST -- a box is threaded through the cells it
 * covers, linear in its length, not through its bounds -- that the index then still finds every box a query
 * overlaps, checked against testing every box, and that an abnormal pass's memory is given back.
 * No GUI: this is geometry.
 */
public final class OrientedOccupancyTest {

    private static final double Q = Math.PI / 4;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "OrientedOccupancy: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        separatingAxis( ok );
        claiming( ok );
        grid( ok );
        nonBehaviour( ok );
        thinBoxCost( ok );
        indexEqualsEveryBox( ok );
        givesBackMemory( ok );
        return ok[ 0 ];
    }

    /** A long label along a 45-degree spoke covers about two cells per column; its bounds hold a square of them.
     *  A 10,000 x 14 px box at 16 px cells: under 3,000 entries against the bounds' 196,000. */
    private static void thinBoxCost( final boolean[] ok ) {
        final OrientedOccupancy o = new OrientedOccupancy();
        o.reset( 16 );
        o.place( 5000, 5000, 5000, 7, Q );
        final int entries = o.cellEntryCountForTest();
        check( ok, entries <= 3000, "a 10,000 x 14 px box at 45 degrees threads a few cells per column, got " + entries );
        check( ok, entries >= 440, "precondition: at least one cell per column, got " + entries );
        final double c = Math.cos( Q ), s = Math.sin( Q );
        check( ok, o.overlaps( 5000 + ( 4900 * c ), 5000 + ( 4900 * s ), 10, 7, Q ), "found at its far end" );
        check( ok, !o.overlaps( 5000 + ( 4900 * c ) - ( 40 * s ), 5000 + ( 4900 * s ) + ( 40 * c ), 10, 7, Q ),
                "a query 40 px beside the far end, inside the box's bounds, is clear" );
        // the END edge of a box reaches into a column on its own: a 200 x 40 box at 30 degrees has its far end
        // edge crossing the column edges at x = 80 and 96, and in the column [80, 96) its highest point (y = 61.4)
        // lies on that end edge alone -- the other three edges reach only y = 32 there. The point (81, 58) is inside
        // the box (99.1 along it, 9.7 across) and in the cell the end edge alone reaches; a query that fits in that
        // one cell must find the box.
        o.reset( 16 );
        o.place( 0, 0, 100, 20, Math.PI / 6 );
        check( ok, o.overlaps( 81, 58, 1, 1, 0 ), "a box is threaded through the cell its end edge alone reaches" );
        check( ok, !o.overlaps( 81, 66, 1, 1, 0 ), "...and a query just past that end is clear" );
    }

    /** The grid is an index, never a filter: for random rotated boxes placed, and random rotated queries, its answer
     *  equals the answer of testing EVERY placed box with the separating-axis test, both ways. Boxes centred on
     *  cell edges and boxes thinner than a cell are included. */
    private static void indexEqualsEveryBox( final boolean[] ok ) {
        final java.util.Random rnd = new java.util.Random( 20261016 );
        final OrientedOccupancy o = new OrientedOccupancy();
        o.reset( 16 );
        final int n = 300;
        final double[][] boxes = new double[ n ][ 5 ];
        for( int i = 0; i < n; ++i ) {
            final double[] b = boxes[ i ];
            b[ 0 ] = ( rnd.nextInt( 4 ) == 0 ) ? ( ( rnd.nextInt( 60 ) - 10 ) * 16.0 ) : ( ( rnd.nextDouble() * 1000 ) - 150 );
            b[ 1 ] = ( rnd.nextInt( 4 ) == 0 ) ? ( ( rnd.nextInt( 60 ) - 10 ) * 16.0 ) : ( ( rnd.nextDouble() * 1000 ) - 150 );
            b[ 2 ] = 1 + ( rnd.nextDouble() * ( ( rnd.nextInt( 3 ) == 0 ) ? 300 : 40 ) ); // half width, some long
            b[ 3 ] = 0.5 + ( rnd.nextDouble() * 8 );                                       // half height
            b[ 4 ] = rnd.nextDouble() * 2 * Math.PI;
            o.place( b[ 0 ], b[ 1 ], b[ 2 ], b[ 3 ], b[ 4 ] ); // placed regardless of overlap, as a found label is
        }
        check( ok, o.claimedCountForTest() == n, "precondition: every box was placed" );
        int hits = 0, misses = 0;
        for( int q = 0; q < 3000; ++q ) {
            final double cx = ( rnd.nextDouble() * 1000 ) - 150, cy = ( rnd.nextDouble() * 1000 ) - 150;
            final double hw = 1 + ( rnd.nextDouble() * 60 ), hh = 0.5 + ( rnd.nextDouble() * 8 );
            final double theta = rnd.nextDouble() * 2 * Math.PI;
            final double qc = Math.cos( theta ), qs = Math.sin( theta );
            boolean expected = false;
            for( int i = 0; ( i < n ) && !expected; ++i ) {
                final double[] b = boxes[ i ];
                // the stored floats, as the grid holds them
                expected = OrientedOccupancy.overlap( cx, cy, hw, hh, qc, qs, (float) b[ 0 ], (float) b[ 1 ],
                                                      (float) b[ 2 ], (float) b[ 3 ], (float) Math.cos( b[ 4 ] ),
                                                      (float) Math.sin( b[ 4 ] ) );
            }
            final boolean got = o.overlaps( cx, cy, hw, hh, theta );
            if ( got != expected ) {
                check( ok, false, "query " + q + " at (" + cx + ", " + cy + "): the index says " + got
                        + ", testing every box says " + expected );
                return;
            }
            if ( expected ) {
                ++hits;
            }
            else {
                ++misses;
            }
        }
        check( ok, ( hits > 200 ) && ( misses > 200 ), "precondition: both answers occur (hits " + hits
                + ", misses " + misses + ")" );
    }

    /** A normal pass keeps its tables across resets; an abnormal one, past SHRINK_ABOVE slots, has them given
     *  back at the next reset -- and the grid works on from there. */
    private static void givesBackMemory( final boolean[] ok ) {
        final OrientedOccupancy o = new OrientedOccupancy();
        o.reset( 8 );
        for( int i = 0; i < 5000; ++i ) {
            o.place( ( i * 8 ) + 4, 4, 2, 2, 0 ); // 5000 cells
        }
        final int normal = o.slotTableSizeForTest();
        check( ok, ( normal >= 8192 ) && ( normal <= OrientedOccupancy.SHRINK_ABOVE ), "precondition: 5000 cells grew the table past its initial size, got " + normal );
        o.reset( 8 );
        check( ok, o.slotTableSizeForTest() == normal, "a normal pass keeps its table across a reset, got " + o.slotTableSizeForTest() );
        for( int i = 0; i < 300000; ++i ) {
            o.place( ( i * 8 ) + 4, 4, 2, 2, 0 ); // 300,000 cells: an abnormal pass
        }
        check( ok, o.slotTableSizeForTest() > OrientedOccupancy.SHRINK_ABOVE, "precondition: 300,000 cells grew the table past the threshold, got " + o.slotTableSizeForTest() );
        o.reset( 8 );
        check( ok, o.slotTableSizeForTest() == 1024, "an abnormal pass's table is given back at the next reset, got " + o.slotTableSizeForTest() );
        check( ok, o.cellEntryCountForTest() == 0, "...and the pass starts empty" );
        check( ok, claim( o, 50, 50, 20, 7, 0 ), "...and the grid works on from the small table" );
        check( ok, !claim( o, 50, 50, 20, 7, 0 ), "...and remembers what it placed" );
    }

    /** The static test, on pairs where the axis-aligned answer and the oriented answer differ. */
    private static void separatingAxis( final boolean[] ok ) {
        // two 100 x 14 labels along the same 45-degree spoke, 20 px apart across it: their bounds (80 x 80 each)
        // overlap almost entirely, the labels themselves do not
        final double c = Math.cos( Q ), s = Math.sin( Q );
        check( ok, !OrientedOccupancy.overlap( 0, 0, 50, 7, c, s, -20 * s, 20 * c, 50, 7, c, s ),
                "two parallel rotated labels 20 px apart across the spoke do not overlap (their bounds do)" );
        check( ok, OrientedOccupancy.overlap( 0, 0, 50, 7, c, s, -10 * s, 10 * c, 50, 7, c, s ),
                "...10 px apart, they do (14 px tall)" );
        check( ok, OrientedOccupancy.overlap( 0, 0, 50, 7, c, s, 60 * c, 60 * s, 50, 7, c, s ),
                "end to end along the spoke with 40 px of overlap, they do" );
        check( ok, !OrientedOccupancy.overlap( 0, 0, 50, 7, c, s, 110 * c, 110 * s, 50, 7, c, s ),
                "end to end with a 10 px gap, they do not" );
        // crossing at right angles: overlap iff the crossing point lies inside both
        check( ok, OrientedOccupancy.overlap( 0, 0, 50, 7, 1, 0, 0, 0, 50, 7, 0, 1 ), "a cross of two labels overlaps" );
        // (a 100 px vertical label centred 90 px up has its near end at 40, 33 px above the horizontal one's top edge
        // at 7; centred at 40 its near end would be at -10, INSIDE the other -- the first draft of this line)
        check( ok, !OrientedOccupancy.overlap( 0, 0, 50, 7, 1, 0, 0, 90, 50, 7, 0, 1 ),
                "a vertical label whose near end stops 33 px short of a horizontal one does not" );
        // touching along an edge is NOT overlap (strict, like LabelOccupancy)
        check( ok, !OrientedOccupancy.overlap( 0, 0, 50, 7, 1, 0, 0, 14, 50, 7, 1, 0 ), "edge to edge, touching, is not overlap" );
        check( ok, OrientedOccupancy.overlap( 0, 0, 50, 7, 1, 0, 0, 13.9, 50, 7, 1, 0 ), "...a tenth of a pixel closer is" );
        // an axis of the SECOND box must be tried too. A short 45-degree box off the flat box's corner: both of the
        // flat box's axes see overlapping projections (x 48.8..67.2 against ..50, y 5.8..24.2 against ..7), and only
        // the tilted box's own direction axis separates them (41.6..61.6 against ..40.3). Two pixels closer and the
        // tilted box's near end pokes into the flat one.
        final double h = Math.sqrt( 0.5 );
        check( ok, !OrientedOccupancy.overlap( 0, 0, 50, 7, 1, 0, 58, 15, 10, 3, h, h ),
                "a tilted box off the corner is separated by ITS OWN axis alone" );
        check( ok, OrientedOccupancy.overlap( 0, 0, 50, 7, 1, 0, 56, 13, 10, 3, h, h ),
                "...two pixels closer its near end is inside the flat box" );
    }

    private static void claiming( final boolean[] ok ) {
        final OrientedOccupancy o = new OrientedOccupancy();
        o.reset( 16 );
        check( ok, claim( o, 100, 100, 50, 7, Q ), "the first label is free" );
        check( ok, o.claimedCountForTest() == 1, "...and recorded" );
        final double c = Math.cos( Q ), s = Math.sin( Q );
        check( ok, !claim( o, 100 - ( 10 * s ), 100 + ( 10 * c ), 50, 7, Q ), "a label 10 px beside it is refused" );
        check( ok, o.claimedCountForTest() == 1, "a refused label records nothing, so it never blocks a later one" );
        check( ok, claim( o, 100 - ( 20 * s ), 100 + ( 20 * c ), 50, 7, Q ), "a label 20 px beside it is free although the bounds overlap" );
        check( ok, o.overlaps( 100, 100, 50, 7, Q ), "a read sees the claimed label" );
        check( ok, o.claimedCountForTest() == 2, "...and records nothing" );
        // placing without asking: a found node's label is drawn regardless, and still blocks what follows
        o.place( 300, 300, 50, 7, 0 );
        check( ok, o.claimedCountForTest() == 3, "a placed label is recorded" );
        check( ok, !claim( o, 320, 305, 50, 7, 0 ), "...and what comes after keeps clear of it" );
    }

    private static void grid( final boolean[] ok ) {
        final OrientedOccupancy o = new OrientedOccupancy();
        o.reset( 16 );
        // a long rotated label spans many cells: a query at its far end must find it
        check( ok, claim( o, 0, 0, 400, 7, Q ), "a long label claims" );
        final double c = Math.cos( Q ), s = Math.sin( Q );
        check( ok, o.overlaps( 380 * c, 380 * s, 10, 7, Q ), "found at its far end, cells away from its centre" );
        check( ok, !o.overlaps( ( 380 * c ) - ( 30 * s ), ( 380 * s ) + ( 30 * c ), 10, 7, Q ),
                "a query 30 px beside the far end, inside the label's bounds, is clear" );
        // many labels: the slot table grows and nothing is lost
        o.reset( 16 );
        for( int i = 0; i < 5000; ++i ) {
            check( ok, claim( o, i * 30, 0, 10, 5, 0 ), "label " + i + " of 5000 in a row claims" );
        }
        check( ok, o.claimedCountForTest() == 5000, "5000 recorded, got " + o.claimedCountForTest() );
        check( ok, o.overlaps( 4999 * 30, 0, 1, 1, 0 ), "the last of 5000 is still found after the table grew" );
        check( ok, o.overlaps( 0, 0, 1, 1, 0 ), "...and so is the first" );
    }

    private static void nonBehaviour( final boolean[] ok ) {
        final OrientedOccupancy o = new OrientedOccupancy();
        o.reset( 16 );
        check( ok, claim( o, 10, 10, 0, 7, 0 ), "a zero-width box is granted" );
        check( ok, o.claimedCountForTest() == 0, "...and reserves nothing" );
        o.place( 10, 10, 5, 0, 0 );
        check( ok, o.claimedCountForTest() == 0, "a zero-height placement reserves nothing" );
        check( ok, claim( o, 50, 50, 20, 7, 0 ), "precondition: a label is placed" );
        // a degenerate QUERY (an empty label: no width, or no height standing across) still looks in its own
        // column -- its edges are all vertical, and the index must answer as testing every box does
        check( ok, o.overlaps( 50, 50, 0, 7, 0 ), "a zero-width query through a placed label still finds it" );
        check( ok, o.overlaps( 50, 50, 20, 0, Math.PI / 2 ), "...and so does a zero-height one standing across it" );
        check( ok, !o.overlaps( 50, 80, 0, 7, 0 ), "...and a zero-width query beside it is clear" );
        o.reset( 16 );
        check( ok, claim( o, 50, 50, 20, 7, 0 ), "a reset forgets every label" );
    }

    /** What production does with a mark: ask, and place only if the way is clear (TreePanel asks for a name and
     *  its image together, then places both -- the reason the class has no single ask-and-record call). */
    private static boolean claim( final OrientedOccupancy o, final double cx, final double cy, final double hw,
                                  final double hh, final double theta ) {
        if ( o.overlaps( cx, cy, hw, hh, theta ) ) {
            return false;
        }
        o.place( cx, cy, hw, hh, theta );
        return true;
    }

    private static void check( final boolean[] ok, final boolean condition, final String what ) {
        if ( !condition ) {
            System.out.println( "  [OrientedOccupancyTest] " + what );
            ok[ 0 ] = false;
        }
    }

    private OrientedOccupancyTest() {
    }
}
