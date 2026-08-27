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

import org.forester.archaeopteryx.ExportSizeSpec.Unit;

/**
 * Headless unit test for {@link ExportSizeSpec}: the pure mm/inch/pixel -> inches -> layout-points / raster-pixels
 * / raster-scale conversions that drive the "export at a fixed size" feature.
 */
public final class ExportSizeSpecTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ExportSizeSpec: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        boolean ok = true;

        // 170 x 120 mm at 300 DPI: inches = 6.6929 x 4.7244; layout pt = 482 x 340; image px = 2008 x 1417
        final ExportSizeSpec mm = new ExportSizeSpec( Unit.MILLIMETERS, 170, 120, 300 );
        ok &= close( "mm widthInches", mm.widthInches(), 170.0 / 25.4 );
        ok &= eq( "mm layoutWidthPt", mm.layoutWidthPt(), 482 );
        ok &= eq( "mm layoutHeightPt", mm.layoutHeightPt(), 340 );
        ok &= eq( "mm imageWidthPx", mm.imageWidthPx(), 2008 );
        ok &= eq( "mm imageHeightPx", mm.imageHeightPx(), 1417 );
        ok &= close( "mm rasterScale", mm.rasterScale(), 300.0 / 72.0 );

        // 6 x 4 inch at 150 DPI: layout pt = 432 x 288; image px = 900 x 600
        final ExportSizeSpec in = new ExportSizeSpec( Unit.INCHES, 6, 4, 150 );
        ok &= eq( "in layoutWidthPt", in.layoutWidthPt(), 432 );
        ok &= eq( "in imageWidthPx", in.imageWidthPx(), 900 );
        ok &= eq( "in imageHeightPx", in.imageHeightPx(), 600 );

        // PIXELS unit: the pixel count round-trips exactly at any DPI (this is the whole point of px mode)
        final ExportSizeSpec px96 = new ExportSizeSpec( Unit.PIXELS, 800, 600, 96 );
        ok &= eq( "px@96 imageWidthPx", px96.imageWidthPx(), 800 );
        ok &= eq( "px@96 imageHeightPx", px96.imageHeightPx(), 600 );
        final ExportSizeSpec px300 = new ExportSizeSpec( Unit.PIXELS, 2400, 1800, 300 );
        ok &= eq( "px@300 imageWidthPx", px300.imageWidthPx(), 2400 );
        ok &= eq( "px@300 imageHeightPx", px300.imageHeightPx(), 1800 );
        // ... but the physical size (hence layout points + embedded DPI) DOES change with the DPI: 2400 px is
        // 8 in at 300 DPI -> 576 pt, whereas at 96 DPI it would be 25 in -> 1800 pt.
        ok &= eq( "px@300 layoutWidthPt", px300.layoutWidthPt(), 576 );

        // 72 DPI -> raster scale 1.0 (raster == layout point size); 288 DPI -> 4.0
        ok &= close( "scale@72", new ExportSizeSpec( Unit.INCHES, 5, 5, 72 ).rasterScale(), 1.0 );
        ok &= close( "scale@288", new ExportSizeSpec( Unit.INCHES, 5, 5, 288 ).rasterScale(), 4.0 );

        // fromInches inverts toInches (preserve physical size across a unit switch)
        ok &= close( "fromInches mm", ExportSizeSpec.fromInches( 170.0 / 25.4, Unit.MILLIMETERS, 300 ), 170.0 );
        ok &= close( "fromInches px", ExportSizeSpec.fromInches( 8.0, Unit.PIXELS, 300 ), 2400.0 );
        ok &= close( "fromInches in", ExportSizeSpec.fromInches( 6.0, Unit.INCHES, 300 ), 6.0 );

        // defensive clamping: a non-positive dimension / DPI never yields a zero/negative canvas
        final ExportSizeSpec zero = new ExportSizeSpec( Unit.MILLIMETERS, 0, -5, 0 );
        ok &= eq( "zero-dim layoutWidthPt >= 1", Math.max( 1, zero.layoutWidthPt() ), zero.layoutWidthPt() );
        ok &= ( zero.layoutWidthPt() >= 1 ) && ( zero.imageHeightPx() >= 1 ) && ( zero.dpi() >= 1 );
        if ( !( ( zero.layoutWidthPt() >= 1 ) && ( zero.imageHeightPx() >= 1 ) ) ) {
            System.out.println( "  [ExportSizeSpecTest] zero/negative inputs must clamp to a positive canvas" );
        }

        // summary mentions the raster px, DPI and the physical mm size
        final String s = mm.summary();
        if ( !s.contains( "2008 × 1417 px" ) || !s.contains( "300 DPI" ) || !s.contains( "170 × 120 mm" )
                || !s.contains( "482 × 340 pt" ) ) {
            ok = false;
            System.out.println( "  [ExportSizeSpecTest] summary missing expected parts: " + s );
        }

        return ok;
    }

    private static boolean eq( final String name, final int a, final int b ) {
        if ( a != b ) {
            System.out.println( "  [ExportSizeSpecTest] " + name + ": got " + a + ", expected " + b );
            return false;
        }
        return true;
    }

    private static boolean close( final String name, final double a, final double b ) {
        if ( Math.abs( a - b ) > 1e-6 ) {
            System.out.println( "  [ExportSizeSpecTest] " + name + ": got " + a + ", expected " + b );
            return false;
        }
        return true;
    }

    private ExportSizeSpecTest() {
    }
}
