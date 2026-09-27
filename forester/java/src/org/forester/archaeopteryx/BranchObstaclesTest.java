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
 * {@link BranchObstacles}: the clip itself, the rotated frame a number is tested in, the owner exclusion, the cell
 * grid (a line far from where it was registered is still found), arcs as chords that hug the circle, and the two
 * deliberate non-behaviours -- a zero-length line registers nothing, and a reset forgets everything. No GUI: this
 * is geometry, and a wrong answer here would be invisible in a rendered figure until the one tree that hits it.
 */
public final class BranchObstaclesTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "BranchObstacles: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        clip( ok );
        frameAndOwner( ok );
        farCells( ok );
        arcs( ok );
        nonBehaviour( ok );
        return ok[ 0 ];
    }

    /** Liang-Barsky against [0,10]x[0,10]: every way a segment can relate to a box. */
    private static void clip( final boolean[] ok ) {
        check( ok, BranchObstacles.segmentMeetsBox( -5, -5, 15, 15, 0, 0, 10, 10 ), "a diagonal through the box meets it" );
        check( ok, BranchObstacles.segmentMeetsBox( 2, 2, 4, 4, 0, 0, 10, 10 ), "a segment wholly inside meets it" );
        check( ok, !BranchObstacles.segmentMeetsBox( -5, -5, -1, -1, 0, 0, 10, 10 ), "a segment short of the corner misses" );
        check( ok, !BranchObstacles.segmentMeetsBox( 20, -5, 20, 15, 0, 0, 10, 10 ), "a vertical line to the right misses" );
        check( ok, !BranchObstacles.segmentMeetsBox( -5, 12, 15, 12, 0, 0, 10, 10 ), "a horizontal line above misses" );
        check( ok, BranchObstacles.segmentMeetsBox( -5, 10, 15, 10, 0, 0, 10, 10 ), "a line along the top edge TOUCHES, and touching counts" );
        check( ok, BranchObstacles.segmentMeetsBox( 10, -5, 10, 15, 0, 0, 10, 10 ), "a line along the right edge touches" );
        // a segment that meets the box at ONE point -- its corner -- is the case where 'touching counts' is decided:
        // the clip's entry and exit parameters coincide there, and a strict comparison would call it a miss
        check( ok, BranchObstacles.segmentMeetsBox( 0, 20, 20, 0, 0, 0, 10, 10 ), "a segment through the corner alone touches, and touching counts" );
        check( ok, BranchObstacles.segmentMeetsBox( 5, 5, 5, 5, 0, 0, 10, 10 ), "a point inside meets" );
        check( ok, !BranchObstacles.segmentMeetsBox( 15, 5, 15, 5, 0, 0, 10, 10 ), "a point outside misses" );
        check( ok, !BranchObstacles.segmentMeetsBox( -20, 0, -1, 30, 0, 0, 10, 10 ), "a diagonal that passes the box's corner region but never enters misses" );
    }

    /**
     * The number's box is ROTATED: a segment that meets it in one orientation and not another pins that the test
     * runs in the box's own frame, and its owner's own line is never an obstacle.
     */
    private static void frameAndOwner( final boolean[] ok ) {
        final BranchObstacles o = new BranchObstacles();
        o.reset( 16 );
        o.add( 1, 50, -100, 50, 100 ); // a vertical line at x = 50
        // a 20 x 6 box centred 6 px right of the line: lying along x it reaches back over the line...
        check( ok, o.crosses( 2, 56, 20, 0, 10, -3, 3, 0 ), "a box across the line, unrotated, crosses it" );
        // ...standing along y (rotated a quarter turn) it is 3 px wide and stops 3 px short of it
        check( ok, !o.crosses( 2, 56, 20, Math.PI / 2, 10, -3, 3, 0 ), "the same box rotated a quarter turn clears the line" );
        check( ok, o.crosses( 2, 56, 20, Math.PI / 2, 10, -3, 3, 3.5 ), "...until the margin reaches it" );
        check( ok, !o.crosses( 1, 56, 20, 0, 10, -3, 3, 0 ), "its own line is never an obstacle" );
        check( ok, !o.crosses( 2, 56, 300, 0, 10, -3, 3, 0 ), "a box far along y, past the line's end, is clear" );
        // the box's OFF-CENTRE band, rotated: local +y is the canvas's +x here (m = a quarter turn maps local
        // (u, v) to device (-v, u)), so a band from 2 to 8 on local +y lies at device x 48..54 and meets the line
        // at x = 50 -- while its mirror image, from -8 to -2, lies at x 58..64 and misses. A sign slip in the frame
        // transform swaps the two answers; a band centred on the line cannot tell (it meets both ways).
        check( ok, o.crosses( 2, 56, 20, Math.PI / 2, 10, 2, 8, 0 ), "a band on the local +y side, rotated, meets the line that lies on that side" );
        check( ok, !o.crosses( 2, 56, 20, Math.PI / 2, 10, -8, -2, 0 ), "...and the band on the other side misses it" );
        // the local box is off-centre (top..bottom need not straddle 0), and the frame's +y is the canvas's
        // when unrotated: a box from 4 to 10 below a horizontal line at y = 0 clears it, one from -2 to 4 does not
        o.reset( 16 );
        o.add( 7, -100, 0, 100, 0 );
        check( ok, !o.crosses( 8, 0, 0, 0, 10, 4, 10, 0 ), "a box wholly below the line clears it" );
        check( ok, o.crosses( 8, 0, 0, 0, 10, -2, 4, 0 ), "a box straddling the line crosses it" );
        check( ok, !o.crosses( 8, 0, 0, 0, 10, -10, -4, 0 ), "a box wholly above the line clears it" );
        // rotated half a turn the same local box lands on the other side of the line
        check( ok, !o.crosses( 8, 0, 0, Math.PI, 10, 4, 10, 0 ), "half a turn round, 'below' is above, and still clear" );
        check( ok, o.crosses( 8, 0, 0, Math.PI, 10, -2, 4, 0 ), "...and a straddling box still crosses" );
    }

    /** A long line lives in many cells; a query at its far end must find it, and one beside it must not. */
    private static void farCells( final boolean[] ok ) {
        final BranchObstacles o = new BranchObstacles();
        o.reset( 16 );
        o.add( 1, 0, 0, 1000, 0 );
        check( ok, o.crosses( 2, 950, 1, 0, 10, -3, 3, 0 ), "a line is found in a cell far from where it starts" );
        check( ok, !o.crosses( 2, 950, 40, 0, 10, -3, 3, 0 ), "a box a few cells off the line is clear" );
        // a diagonal line: its cells are those of its bounds, and a box beside it inside those bounds is clear
        o.add( 3, 0, 100, 400, 500 );
        check( ok, o.crosses( 2, 200, 300, 0, 10, -3, 3, 0 ), "a diagonal is met on the diagonal" );
        check( ok, !o.crosses( 2, 350, 150, 0, 10, -3, 3, 0 ), "a box inside a diagonal's bounds but off it is clear" );
        // many lines: the slot table grows and nothing is lost
        o.reset( 16 );
        for( int i = 0; i < 5000; ++i ) {
            o.add( i, i * 3, 0, ( i * 3 ) + 2, 0 );
        }
        check( ok, o.segmentCountForTest() == 5000, "5000 lines registered, got " + o.segmentCountForTest() );
        check( ok, o.crosses( -1, 14999, 0, 0, 1, -1, 1, 0 ), "the last of 5000 lines is still found after the table grew" );
        check( ok, o.crosses( -1, 3, 0, 0, 1, -1, 1, 0 ), "...and so is the first" );
    }

    /** An arc is registered as chords that never sag more than half a pixel from the circle. */
    private static void arcs( final boolean[] ok ) {
        final BranchObstacles o = new BranchObstacles();
        o.reset( 16 );
        final double r = 100;
        o.addArc( 0, 0, r, 0, Math.PI / 2 ); // a quarter circle in the +x/+y quadrant
        check( ok, o.segmentCountForTest() > 1, "an arc of radius 100 needs several chords, got " + o.segmentCountForTest() );
        // a thin box ON the circle, at an angle no chord ends at, must be met: that is the sag bound at work. The box
        // reaches 0.6 px either side of the circle -- just past the half-pixel a chord may sag inward (measured: at
        // r = 100 the 8 chords of a quarter circle sag 0.48 px, so a 0.45 px box missed by three hundredths).
        final double a = 0.1;
        check( ok, o.crosses( 2, r * Math.cos( a ), r * Math.sin( a ), a + ( Math.PI / 2 ), 2, -0.6, 0.6, 0 ),
                "a 1.2 px thin box lying on the arc between chord ends is met (chords sag under half a pixel)" );
        check( ok, !o.crosses( 2, 80 * Math.cos( a ), 80 * Math.sin( a ), a + ( Math.PI / 2 ), 2, -0.6, 0.6, 0 ),
                "the same box 20 px inside the arc is clear" );
        check( ok, o.crosses( 1, r * Math.cos( a ), r * Math.sin( a ), a + ( Math.PI / 2 ), 2, -0.6, 0.6, 0 ),
                "an arc is nobody's own: the branch that drew it is not excused from it either" );
        check( ok, !o.crosses( 2, r * Math.cos( -0.3 ), r * Math.sin( -0.3 ), 0, 2, -2, 2, 0 ),
                "a box on the circle but outside the arc's sweep is clear" );
        // sweeping the other way registers the same circle
        o.reset( 16 );
        o.addArc( 0, 0, r, Math.PI / 2, 0 );
        check( ok, o.crosses( 2, r * Math.cos( a ), r * Math.sin( a ), a + ( Math.PI / 2 ), 2, -0.6, 0.6, 0 ),
                "an arc swept backwards is the same obstacle" );
    }

    /** What is deliberately NOT an obstacle. */
    private static void nonBehaviour( final boolean[] ok ) {
        final BranchObstacles o = new BranchObstacles();
        o.reset( 16 );
        o.add( 1, 40, 40, 40, 40 ); // a zero-length branch: the painter draws nothing for it
        check( ok, o.segmentCountForTest() == 0, "a zero-length line registers nothing" );
        check( ok, !o.crosses( 2, 40, 40, 0, 10, -3, 3, 1 ), "...so a number on top of it is clear" );
        o.addArc( 0, 0, 0, 0, 1 );
        check( ok, o.segmentCountForTest() == 0, "an arc of zero radius registers nothing" );
        o.add( 1, 0, 0, 100, 0 );
        check( ok, o.crosses( 2, 50, 0, 0, 10, -3, 3, 0 ), "precondition: a line is in the way" );
        o.reset( 16 );
        check( ok, !o.crosses( 2, 50, 0, 0, 10, -3, 3, 0 ), "a reset forgets every line" );
        check( ok, o.segmentCountForTest() == 0, "...and the count says so" );
    }

    private static void check( final boolean[] ok, final boolean condition, final String what ) {
        if ( !condition ) {
            System.out.println( "  [BranchObstaclesTest] " + what );
            ok[ 0 ] = false;
        }
    }

    private BranchObstaclesTest() {
    }
}
