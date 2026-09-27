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
 * nothing; a zero-size box reserves nothing).
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
        return ok[ 0 ];
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
