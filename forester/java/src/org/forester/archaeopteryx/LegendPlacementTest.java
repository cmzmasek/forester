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

import java.awt.Point;
import java.awt.Rectangle;

/**
 * Unit tests for {@link TreePanel#legendTopLeftFor}: the dragged property-color legend position must
 * be honored not only on screen but also in PDF/PNG/etc. exports. The bug it guards against: exports
 * snapped the legend back to the default top-right corner, ignoring the user's move. The dragged
 * position (a pixel offset within the on-screen viewport) is mapped onto the export target by its
 * fractional position, so e.g. "bottom-right of my view" exports as "bottom-right of the image". In
 * {@code org.forester.archaeopteryx} so it can reach the package-private method; run via the suite
 * ({@link #test()}).
 */
public final class LegendPlacementTest {

    public static void main( final String[] args ) {
        System.out.println( "LegendPlacement: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        final int bw = 100;
        final int bh = 60;
        final Rectangle vp = new Rectangle( 0, 0, 800, 600 );

        // 1. no drag -> default top-right corner of the target, inset 10px
        Point p = TreePanel.legendTopLeftFor( vp, vp, null, bw, bh );
        if ( ( p.x != ( vp.width - bw - 10 ) ) || ( p.y != 10 ) ) {
            return fail( "default corner wrong: " + p );
        }

        // 2. on screen (target == viewport): exact dragged pixel position, no drift
        final Point off = new Point( 250, 400 );
        p = TreePanel.legendTopLeftFor( vp, vp, off, bw, bh );
        if ( ( p.x != 250 ) || ( p.y != 400 ) ) {
            return fail( "on-screen exact placement wrong: " + p );
        }

        // 2b. "visible region only" export: target is the (scrolled) viewport; placement matches the
        //     same absolute coordinates the user sees
        final Rectangle vis = new Rectangle( 120, 90, 800, 600 );
        p = TreePanel.legendTopLeftFor( vis, vis, off, bw, bh );
        if ( ( p.x != ( 120 + 250 ) ) || ( p.y != ( 90 + 400 ) ) ) {
            return fail( "visible-only export placement wrong: " + p );
        }

        // 3. full-tree export onto a much larger image: proportional mapping.
        //    This is the case the old code got wrong (it used the default corner).
        final Rectangle full = new Rectangle( 0, 0, 4000, 3000 );
        // 3a. far corner of the viewport -> far corner of the image
        p = TreePanel.legendTopLeftFor( full, vp, new Point( vp.width - bw, vp.height - bh ), bw, bh );
        if ( ( p.x != ( full.width - bw ) ) || ( p.y != ( full.height - bh ) ) ) {
            return fail( "full-export bottom-right wrong: " + p );
        }
        // 3b. origin stays at the origin
        p = TreePanel.legendTopLeftFor( full, vp, new Point( 0, 0 ), bw, bh );
        if ( ( p.x != 0 ) || ( p.y != 0 ) ) {
            return fail( "full-export top-left wrong: " + p );
        }
        // 3c. centered in the viewport -> centered in the image
        p = TreePanel.legendTopLeftFor( full, vp, new Point( ( vp.width - bw ) / 2, ( vp.height - bh ) / 2 ), bw, bh );
        if ( ( Math.abs( p.x - ( ( full.width - bw ) / 2 ) ) > 2 )
                || ( Math.abs( p.y - ( ( full.height - bh ) / 2 ) ) > 2 ) ) {
            return fail( "full-export center wrong: " + p );
        }

        // 4. an over-range offset is clamped so the legend never lands off the image
        p = TreePanel.legendTopLeftFor( full, vp, new Point( 999999, 999999 ), bw, bh );
        if ( ( p.x != ( full.width - bw ) ) || ( p.y != ( full.height - bh ) ) ) {
            return fail( "over-range offset not clamped: " + p );
        }

        // 5. a degenerate viewport narrower/shorter than the box must not divide by zero
        p = TreePanel.legendTopLeftFor( full, new Rectangle( 0, 0, 50, 40 ), off, bw, bh );
        if ( ( p.x != 0 ) || ( p.y != 0 ) ) {
            return fail( "tiny viewport should fall back to the origin: " + p );
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [LegendPlacementTest] " + message );
        return false;
    }

    private LegendPlacementTest() {
    }
}
