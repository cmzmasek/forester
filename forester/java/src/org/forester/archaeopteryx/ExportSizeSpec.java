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

import java.util.Locale;

/**
 * A resolved target size for a fixed-size ("export at exactly this size") graphics export. Pure and GUI-free
 * (no Swing / {@link TreePanel} dependency), so its arithmetic is unit-testable headlessly.
 *
 * <p><b>Everything reduces to a physical size in inches plus a DPI.</b> A width/height is given in millimetres,
 * inches, or pixels; pixels are interpreted as a pixel count at the chosen DPI (so a pixel value round-trips to
 * exactly that many raster pixels). From the physical size two independent quantities are derived:
 * <ul>
 *   <li>the <b>layout size in points</b> ({@code inches * 72}) -- the tree is laid out into this point-space
 *       canvas, so a font at its point size is physically correct on the page. This drives the vector page size
 *       directly (SVG / EPS / PDF are resolution-independent, so the DPI does not apply to them); and</li>
 *   <li>the <b>raster scale</b> ({@code dpi / 72}) -- the point-space layout is rasterised at this scale, so the
 *       output image is {@code round(layoutPt * dpi/72)} px and carries the requested DPI. This keeps the raster
 *       and vector renders consistent (both lay the tree out in the same point-space).</li>
 * </ul>
 *
 * <p>Instances are immutable and defensively clamp their inputs (a non-positive dimension or DPI can never
 * produce a zero/negative canvas). The very large end is bounded downstream: the raster path caps the image at
 * ~100 MP (see {@link AptxUtil#renderPhylogenyToImage}), and a vector page is resolution-independent.
 */
final class ExportSizeSpec {

    /** The unit a fixed export size is expressed in. {@code toString()} is the user-facing combo label. */
    enum Unit {
        MILLIMETERS( "Millimeters (mm)", "mm" ),
        INCHES( "Inches (in)", "in" ),
        PIXELS( "Pixels (px)", "px" );

        private final String _label;
        private final String _abbrev;

        Unit( final String label, final String abbrev ) {
            _label = label;
            _abbrev = abbrev;
        }

        String abbrev() {
            return _abbrev;
        }

        @Override
        public String toString() {
            return _label;
        }
    }

    static final double MM_PER_INCH  = 25.4;
    static final double POINTS_PER_INCH = 72.0;
    /** Defensive upper bound on any derived dimension, so a pathological input can't overflow an int canvas. */
    private static final int MAX_DIM = 200_000;

    private final Unit   _unit;
    private final double _width;
    private final double _height;
    private final int    _dpi;

    ExportSizeSpec( final Unit unit, final double width, final double height, final int dpi ) {
        _unit = ( unit == null ) ? Unit.MILLIMETERS : unit;
        _width = ( width > 0 ) ? width : 0;
        _height = ( height > 0 ) ? height : 0;
        _dpi = ( dpi > 0 ) ? dpi : 1;
    }

    Unit unit() {
        return _unit;
    }

    int dpi() {
        return _dpi;
    }

    /** The width converted to inches (mm/25.4, inches as-is, pixels/DPI). */
    double widthInches() {
        return toInches( _width );
    }

    /** The height converted to inches. */
    double heightInches() {
        return toInches( _height );
    }

    private double toInches( final double v ) {
        switch ( _unit ) {
            case MILLIMETERS:
                return v / MM_PER_INCH;
            case PIXELS:
                return v / _dpi;
            case INCHES:
            default:
                return v;
        }
    }

    /** Converts an inches value BACK into {@code unit} at {@code dpi} -- used when the user switches units and the
     *  numeric value should follow to preserve the physical size. */
    static double fromInches( final double inches, final Unit unit, final int dpi ) {
        final int d = ( dpi > 0 ) ? dpi : 1;
        switch ( unit ) {
            case MILLIMETERS:
                return inches * MM_PER_INCH;
            case PIXELS:
                return inches * d;
            case INCHES:
            default:
                return inches;
        }
    }

    /** The layout width in points (1 pt = 1/72"): the tree is laid out into this so fonts are physically correct.
     *  This IS the vector (SVG/EPS/PDF) page width. */
    int layoutWidthPt() {
        return clamp( (int) Math.round( widthInches() * POINTS_PER_INCH ) );
    }

    /** The layout / vector-page height in points. */
    int layoutHeightPt() {
        return clamp( (int) Math.round( heightInches() * POINTS_PER_INCH ) );
    }

    /** The factor the point-space layout is rasterised by ({@code dpi/72}); the raster image is
     *  {@code round(layoutPt * rasterScale())} px. 1.0 at 72 DPI. */
    double rasterScale() {
        return _dpi / POINTS_PER_INCH;
    }

    /** The raster image width in pixels actually produced ({@code round(layoutWidthPt * dpi/72)}), so a summary
     *  and a completion message match the file exactly. For a {@code PIXELS}-unit spec this equals the value the
     *  user typed (a pixel count round-trips). */
    int imageWidthPx() {
        return clamp( (int) Math.round( layoutWidthPt() * rasterScale() ) );
    }

    /** The raster image height in pixels actually produced. */
    int imageHeightPx() {
        return clamp( (int) Math.round( layoutHeightPt() * rasterScale() ) );
    }

    private static int clamp( final int v ) {
        if ( v < 1 ) {
            return 1;
        }
        return ( v > MAX_DIM ) ? MAX_DIM : v;
    }

    /** A one-line, unit-independent description of exactly what will be produced: the raster pixel size + DPI, the
     *  physical size (mm and inches), and the vector page size in points. */
    String summary() {
        final double mm_w = widthInches() * MM_PER_INCH;
        final double mm_h = heightInches() * MM_PER_INCH;
        return String.format( Locale.ROOT,
                              "%d × %d px @ %d DPI  ·  %s × %s mm (%s × %s in)  ·  "
                                      + "vector page %d × %d pt",
                              imageWidthPx(), imageHeightPx(), _dpi,
                              trim( mm_w ), trim( mm_h ), trim( widthInches() ), trim( heightInches() ),
                              layoutWidthPt(), layoutHeightPt() );
    }

    /** Compact number: an integer prints without a decimal, else one decimal place. */
    private static String trim( final double v ) {
        if ( Math.abs( v - Math.rint( v ) ) < 0.05 ) {
            return Long.toString( Math.round( v ) );
        }
        return String.format( Locale.ROOT, "%.1f", v );
    }
}
