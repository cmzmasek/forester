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
 * deliberate non-behaviours -- a zero-length line registers nothing, and a reset forgets everything. Since
 * 0.11.181 also the COST: a line is threaded through the cells it passes, linear in its length (the bounding box
 * was quadratic and a circular paint ran out of heap in this class's table); that the index then still finds
 * every line a box meets, checked against testing every line; and that an abnormal pass's memory is given back.
 * No GUI: this is geometry, and a wrong answer here would be invisible in a rendered figure until the one tree
 * that hits it.
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
        linearCost( ok );
        indexEqualsEveryLine( ok );
        givesBackMemory( ok );
        return ok[ 0 ];
    }

    /** A line costs cells in proportion to its LENGTH, not to the area of its bounds: the 10,000 px diagonal leg a
     *  radial layout draws at its maximum zoom threads about two cells per column at 16 px cells (~1,250 entries),
     *  where its bounding box holds 390,625 -- and it is still met along its whole length, and not beside it. */
    private static void linearCost( final boolean[] ok ) {
        final BranchObstacles o = new BranchObstacles();
        o.reset( 16 );
        o.add( 1, 0, 0, 10000, 10000 );
        final int entries = o.cellEntryCountForTest();
        check( ok, entries <= 2600, "a 10,000 px diagonal threads at most a few cells per column, got " + entries );
        check( ok, entries >= 625, "precondition: it was threaded through at least one cell per column, got " + entries );
        check( ok, o.crosses( 2, 7000, 7000, 0, 10, -3, 3, 0 ), "the diagonal is met on itself, far from its start" );
        check( ok, !o.crosses( 2, 7000, 3000, 0, 10, -3, 3, 0 ), "...and not inside its bounds, off it" );
        // a horizontal and a vertical line: one cell per column / per row
        o.reset( 16 );
        o.add( 1, 0, 100, 10000, 100 );
        check( ok, o.cellEntryCountForTest() <= 1300, "a horizontal 10,000 px line threads about one cell per column, got " + o.cellEntryCountForTest() );
        o.reset( 16 );
        o.add( 1, 100, 0, 100, 10000 );
        check( ok, o.cellEntryCountForTest() <= 1300, "a vertical 10,000 px line threads about one cell per row, got " + o.cellEntryCountForTest() );
        check( ok, o.crosses( 2, 100, 9990, 0, 10, -3, 3, 0 ), "...and is met at its far end" );
        // a NEAR-vertical line stays in one column: the piece inside the column is the line's own extent, not the
        // line extrapolated to the column's edges (which for a slope of 20,000 would reach 240,000 px past its end)
        o.reset( 16 );
        o.add( 1, 100, 0, 100.5f, 10000 );
        check( ok, o.cellEntryCountForTest() <= 1300, "a near-vertical 10,000 px line threads its own rows only, got " + o.cellEntryCountForTest() );
        check( ok, o.crosses( 2, 100.25, 5000, 0, 10, -3, 3, 0 ), "...and is met in the middle" );
    }

    /** The grid is an index, never a filter: for random lines and random rotated boxes its answer equals the answer
     *  of testing EVERY line, both ways -- lines ending exactly on cell edges and crossing cell corners included,
     *  since those are where a threading rule can leave out the cell a line touches. Same arithmetic as crosses,
     *  so the only thing under test is which lines the index offers. */
    private static void indexEqualsEveryLine( final boolean[] ok ) {
        final java.util.Random rnd = new java.util.Random( 20261016 );
        final BranchObstacles o = new BranchObstacles();
        o.reset( 16 );
        final int n = 300;
        final float[][] lines = new float[ n ][ 4 ];
        for( int i = 0; i < n; ++i ) {
            for( int k = 0; k < 4; ++k ) {
                // one coordinate in four lands exactly on a cell edge
                lines[ i ][ k ] = ( rnd.nextInt( 4 ) == 0 ) ? ( ( rnd.nextInt( 60 ) - 10 ) * 16f )
                        : (float) ( ( rnd.nextDouble() * 1000 ) - 150 );
            }
            o.add( i, lines[ i ][ 0 ], lines[ i ][ 1 ], lines[ i ][ 2 ], lines[ i ][ 3 ] );
        }
        int hits = 0, misses = 0;
        for( int q = 0; q < 3000; ++q ) {
            final double mx = ( rnd.nextDouble() * 1000 ) - 150, my = ( rnd.nextDouble() * 1000 ) - 150;
            final double m = rnd.nextDouble() * 2 * Math.PI;
            final double hw = 1 + ( rnd.nextDouble() * 60 ), top = -( 1 + ( rnd.nextDouble() * 10 ) );
            final double bottom = 1 + ( rnd.nextDouble() * 10 ), margin = rnd.nextDouble() * 2;
            final long owner = rnd.nextInt( n + 5 ) - 2; // sometimes a real line's owner, excused
            boolean expected = false;
            final double cos = Math.cos( m ), sin = Math.sin( m );
            for( int i = 0; ( i < n ) && !expected; ++i ) {
                final float[] l = lines[ i ];
                if ( ( i == owner ) || ( ( l[ 0 ] == l[ 2 ] ) && ( l[ 1 ] == l[ 3 ] ) ) ) {
                    continue; // its own line, or one that registers nothing
                }
                final double dx0 = l[ 0 ] - mx, dy0 = l[ 1 ] - my, dx1 = l[ 2 ] - mx, dy1 = l[ 3 ] - my;
                final double u0 = ( dx0 * cos ) + ( dy0 * sin ), v0 = ( dy0 * cos ) - ( dx0 * sin );
                final double u1 = ( dx1 * cos ) + ( dy1 * sin ), v1 = ( dy1 * cos ) - ( dx1 * sin );
                expected = BranchObstacles.segmentMeetsBox( u0, v0, u1, v1, -( hw + margin ), top - margin,
                                                            hw + margin, bottom + margin );
            }
            final boolean got = o.crosses( owner, mx, my, m, hw, top, bottom, margin );
            if ( got != expected ) {
                check( ok, false, "query " + q + " at (" + mx + ", " + my + "): the index says " + got
                        + ", testing every line says " + expected );
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

    /** A normal pass keeps its tables across resets (no allocation in the steady state); an abnormal one, past
     *  SHRINK_ABOVE slots, has them given back at the next reset -- and the grid works on from there. */
    private static void givesBackMemory( final boolean[] ok ) {
        final BranchObstacles o = new BranchObstacles();
        o.reset( 8 );
        for( int i = 0; i < 5000; ++i ) {
            o.add( i, i * 8, 0, ( i * 8 ) + 4, 0 ); // 5000 cells
        }
        final int normal = o.slotTableSizeForTest();
        check( ok, ( normal >= 8192 ) && ( normal <= BranchObstacles.SHRINK_ABOVE ), "precondition: 5000 cells grew the table past its initial size, got " + normal );
        o.reset( 8 );
        check( ok, o.slotTableSizeForTest() == normal, "a normal pass keeps its table across a reset, got " + o.slotTableSizeForTest() );
        for( int i = 0; i < 300000; ++i ) {
            o.add( i, i * 8, 0, ( i * 8 ) + 4, 0 ); // 300,000 cells: an abnormal pass
        }
        check( ok, o.slotTableSizeForTest() > BranchObstacles.SHRINK_ABOVE, "precondition: 300,000 cells grew the table past the threshold, got " + o.slotTableSizeForTest() );
        o.reset( 8 );
        check( ok, o.slotTableSizeForTest() == 1024, "an abnormal pass's table is given back at the next reset, got " + o.slotTableSizeForTest() );
        check( ok, o.cellEntryCountForTest() == 0, "...and the pass starts empty" );
        o.add( 1, 0, 0, 100, 0 );
        check( ok, o.crosses( 2, 50, 0, 0, 10, -3, 3, 0 ), "...and the grid works on from the small table" );
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
        o.add( 1, Float.NaN, 0, 10, 10 );
        o.add( 1, 0, 0, Float.POSITIVE_INFINITY, 10 );
        check( ok, o.segmentCountForTest() == 0, "a line with a NaN or infinite end draws as nothing and registers nothing" );
        offCanvas( ok );
        o.add( 1, 0, 0, 100, 0 );
        check( ok, o.crosses( 2, 50, 0, 0, 10, -3, 3, 0 ), "precondition: a line is in the way" );
        o.reset( 16 );
        check( ok, !o.crosses( 2, 50, 0, 0, 10, -3, 3, 0 ), "a reset forgets every line" );
        check( ok, o.segmentCountForTest() == 0, "...and the count says so" );
    }

    /** What is off every canvas registers nothing, and costs nothing: an internal node laid out at infinity
     *  (2026-10-16, a tree opened straight into the circular layout) handed addArc an infinite radius, whose
     *  chord count -- the sweep over an arc piece of 0 -- came out as two billion; each chord at infinity, so the
     *  old bounding-box threading filed each in one far cell until the heap was gone (0.11.180), and the linear
     *  threading refused each, two billion times per arc, for a hang (0.11.181). The bound is MAX_EXTENT: past it
     *  a float places no pixel and no canvas reaches. A sweep past a full turn is the same circle. */
    private static void offCanvas( final boolean[] ok ) {
        final BranchObstacles o = new BranchObstacles();
        o.reset( 16 );
        final long t0 = System.nanoTime();
        o.addArc( 370, 406, Double.POSITIVE_INFINITY, 4.71, 4.72 ); // the arc of 2026-10-16, as it was called
        o.addArc( 0, 0, Double.NaN, 0, 1 );
        o.addArc( 0, 0, 100, Double.NaN, 1 );
        o.addArc( 0, 0, 100, 0, Double.POSITIVE_INFINITY );
        o.addArc( Double.POSITIVE_INFINITY, 0, 100, 0, 1 ); // every chord off the canvas: add refuses each
        o.addArc( 0, 0, 2.0 * BranchObstacles.MAX_EXTENT, 0, 1 );
        o.addArc( 0, 0, Double.MAX_VALUE, 0, 1 ); // finite, and an arc piece of 0: two billion chords without the bound
        o.add( 1, 0, 0, 2.0f * BranchObstacles.MAX_EXTENT, 10 );
        o.add( 1, -2.0f * BranchObstacles.MAX_EXTENT, 0, 10, 10 );
        final long ms = ( System.nanoTime() - t0 ) / 1_000_000;
        check( ok, o.segmentCountForTest() == 0, "an arc or line with an infinite, NaN or beyond-MAX_EXTENT part registers nothing" );
        check( ok, o.cellEntryCountForTest() == 0, "...and threads no cell" );
        check( ok, ms < 2000, "...and returns at once, not after two billion chords (took " + ms + " ms)" );
        // the bound is what the canvas can hold, not a round number below it: an arc AT the bound is registered
        o.addArc( 0, 0, BranchObstacles.MAX_EXTENT, 0, 1e-6 );
        check( ok, o.segmentCountForTest() > 0, "an arc at MAX_EXTENT is still a drawable arc and registers" );
        // a sweep past a full turn is the same circle: as many chords as one turn, not as many as the angle says
        o.reset( 16 );
        o.addArc( 0, 0, 100, 0, 2 * Math.PI );
        final int one_turn = o.segmentCountForTest();
        o.reset( 16 );
        o.addArc( 0, 0, 100, 0, 200 * Math.PI );
        check( ok, one_turn > 0, "precondition: a full turn registers chords" );
        check( ok, o.segmentCountForTest() == one_turn, "a sweep of a hundred turns registers one turn of chords, " + o.segmentCountForTest() + " vs " + one_turn );
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
