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
 * Tests the occupancy grid behind "auto-hide crowded data". Pure arithmetic, so it runs everywhere.
 */
public final class LabelOccupancyTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LabelOccupancy: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        final LabelOccupancy o = new LabelOccupancy();
        o.reset( 16 );

        // (1) free space is granted; the same space a second time is not
        if ( !o.claim( 100, 100, 30, 12 ) ) {
            fail( ok, "the first claim on empty space must be granted" );
        }
        if ( o.claim( 100, 100, 30, 12 ) ) {
            fail( ok, "the identical box must be refused the second time" );
        }
        if ( o.claim( 120, 105, 30, 12 ) ) {
            fail( ok, "a partly overlapping box must be refused" );
        }

        // (2) ...but a box that merely ABUTS is not an overlap -- otherwise everything tightly packed would vanish
        if ( !o.claim( 130, 100, 30, 12 ) ) {
            fail( ok, "a box starting exactly where the last one ends must be granted" );
        }
        if ( !o.claim( 100, 112, 30, 12 ) ) {
            fail( ok, "a box directly below, sharing only an edge, must be granted" );
        }

        // (3) a REFUSED box must not be recorded -- else a rejected mark would go on blocking later ones
        final LabelOccupancy o3 = new LabelOccupancy();
        o3.reset( 16 );
        o3.claim( 0, 0, 10, 10 );
        if ( o3.claim( 5, 0, 10, 10 ) ) {
            fail( ok, "precondition: that box overlaps and must be refused" );
        }
        // the refused box covered x 5..15; a box at 12..22 overlaps only IT, not the accepted one
        if ( !o3.claim( 12, 0, 10, 10 ) ) {
            fail( ok, "a refused box must not be recorded, so it cannot block a later claim" );
        }

        // (4) boxes wider than a cell are still found (the grid must test every cell a box spans)
        final LabelOccupancy o4 = new LabelOccupancy();
        o4.reset( 4 ); // deliberately much smaller than the boxes
        if ( !o4.claim( 0, 0, 100, 10 ) ) {
            fail( ok, "the first wide box must be granted" );
        }
        if ( o4.claim( 90, 0, 100, 10 ) ) {
            fail( ok, "an overlap at the FAR end of a multi-cell box must still be detected" );
        }
        if ( !o4.claim( 100, 0, 100, 10 ) ) {
            fail( ok, "...and a box just past its end must be granted" );
        }

        // (4b) a cell holding SEVERAL boxes must keep them ALL reachable. Threading the list through a
        // per-BOX link cannot express a box that spans several cells -- each extra cell overwrites the previous
        // cell's link and truncates its chain -- and the truncation only shows when a probe lands in a shared
        // cell and overlaps the EARLIER occupant, which is what this builds. (An A/B against the map-based grid
        // it replaced is what found the bug: 12000 placed against 12070 over the same 20k claims.)
        // cell = 4: A covers rows 0-1, B covers rows 1-3, so row 1 holds BOTH, with B linked second.
        final LabelOccupancy multi = new LabelOccupancy();
        multi.reset( 4 );
        if ( !multi.claim( 0, 0, 20, 6 ) ) {       // A: y 0..6
            fail( ok, "A must be granted" );
        }
        if ( !multi.claim( 0, 7, 20, 6 ) ) {       // B: y 7..13, clear of A, but shares row 1 with it
            fail( ok, "B is clear of A in y and must be granted" );
        }
        // The probe lands in row 1 ONLY, and overlaps A ONLY. Reachable solely by walking past B to A.
        if ( multi.claim( 0, 4, 20, 2 ) ) {        // y 4..6: overlaps A (0..6), not B (7..13)
            fail( ok, "a probe overlapping the EARLIER occupant of a shared cell must be refused -- the cell's "
                    + "chain was truncated at its newest box" );
        }
        // ...and the mirror, so this is not simply "everything is refused"
        if ( !multi.claim( 200, 4, 20, 2 ) ) {
            fail( ok, "the same probe in empty space must be granted" );
        }

        // (5) a degenerate box reserves nothing and is always granted (nothing is drawn, nothing to protect)
        final LabelOccupancy o5 = new LabelOccupancy();
        o5.reset( 16 );
        if ( !o5.claim( 10, 10, 0, 12 ) || !o5.claim( 10, 10, 30, 0 ) || !o5.claim( 10, 10, -5, 12 ) ) {
            fail( ok, "a zero/negative-size box must be granted and reserve nothing" );
        }
        if ( o5.claimedCountForTest() != 0 ) {
            fail( ok, "a degenerate box must reserve NOTHING, but the map holds "
                    + o5.claimedCountForTest() + " mark(s)" );
        }
        if ( !o5.claim( 10, 10, 30, 12 ) ) {
            fail( ok, "...so a real box in that same space is still free to take it" );
        }

        // (6) reset really clears
        final LabelOccupancy o6 = new LabelOccupancy();
        o6.reset( 16 );
        o6.claim( 50, 50, 20, 20 );
        o6.reset( 16 );
        if ( !o6.claim( 50, 50, 20, 20 ) ) {
            fail( ok, "reset must clear the map, or the next paint inherits the last one's marks" );
        }

        // (7) negative coordinates work (a scrolled/rotated layout hands them out)
        final LabelOccupancy o7 = new LabelOccupancy();
        o7.reset( 16 );
        if ( !o7.claim( -100, -100, 20, 10 ) ) {
            fail( ok, "a box at negative coordinates must be granted" );
        }
        if ( o7.claim( -95, -95, 20, 10 ) ) {
            fail( ok, "...and overlaps there must still be detected" );
        }

        // (8) DETERMINISM: the same sequence must accept exactly the same boxes every time -- this is what makes
        // a figure reproducible, and what lets screen and export agree
        final int[] first = acceptedPattern();
        for( int rep = 0; rep < 3; ++rep ) {
            final int[] again = acceptedPattern();
            for( int i = 0; i < first.length; ++i ) {
                if ( first[ i ] != again[ i ] ) {
                    fail( ok, "the same sequence of claims must give the same answers every time (box " + i + ")" );
                    break;
                }
            }
        }

        // (9) counting: a box spanning many cells is still ONE placed mark
        final LabelOccupancy o9 = new LabelOccupancy();
        o9.reset( 4 );
        o9.claim( 0, 0, 100, 10 );
        o9.claim( 200, 0, 100, 10 );
        if ( o9.claimedCountForTest() != 2 ) {
            fail( ok, "two multi-cell boxes must count as two marks, got " + o9.claimedCountForTest() );
        }

        // (10) it must stay LINEAR: 20k claims is the scale of a big tree, and a quadratic map would crawl here
        final LabelOccupancy big = new LabelOccupancy();
        big.reset( 16 );
        final long started = System.nanoTime();
        int placed = 0;
        for( int i = 0; i < 20000; ++i ) {
            if ( big.claim( ( i * 37 ) % 4000, ( i * 53 ) % 3000, 28, 12 ) ) {
                ++placed;
            }
        }
        final long ms = ( System.nanoTime() - started ) / 1_000_000;
        if ( ms > 500 ) {
            fail( ok, "20000 claims took " + ms + " ms -- the grid is not doing its job" );
        }
        if ( placed < 1 ) {
            fail( ok, "the sweep placed nothing, so it measured nothing" );
        }
        return ok[ 0 ];
    }

    /** A fixed sequence of overlapping and non-overlapping claims; 1 = granted, 0 = refused. */
    private static int[] acceptedPattern() {
        final LabelOccupancy o = new LabelOccupancy();
        o.reset( 16 );
        final int[] out = new int[ 12 ];
        for( int i = 0; i < out.length; ++i ) {
            out[ i ] = o.claim( ( i % 4 ) * 20, ( i / 4 ) * 6, 30, 12 ) ? 1 : 0;
        }
        return out;
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [LabelOccupancyTest] " + message );
        ok[ 0 ] = false;
    }

    private LabelOccupancyTest() {
        // not instantiable
    }
}
