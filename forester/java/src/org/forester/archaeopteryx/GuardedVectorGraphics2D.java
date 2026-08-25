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

import java.awt.Color;
import java.awt.GradientPaint;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.geom.Rectangle2D;

import de.erichseifert.vectorgraphics2d.VectorGraphics2D;

/**
 * A {@link VectorGraphics2D} for SVG / EPS export that FLATTENS any non-{@link Color} paint to a solid colour.
 *
 * <p>VectorGraphics2D cannot express a {@link GradientPaint} / {@link java.awt.TexturePaint}: setting one arms an
 * internal "paint the shape as an image" filter that then rasterizes EVERY subsequent fill into a
 * {@code BufferedImage} the size of the shape's rounded bounds -- and that filter is cleared only by a
 * {@code dispose()}, which Archaeopteryx's single-context paint pass never issues. Two failures follow, both
 * reproduced on real domain trees (the domain box is the only {@code GradientPaint} the tree paint draws):
 * <ul>
 *   <li><b>Crash:</b> once a domain's gradient body is drawn, a later sub-pixel fill (a second domain's shadow, an
 *       empty glyph outline) rounds to a zero dimension, so {@code new BufferedImage(0, ..)} throws
 *       {@code IllegalArgumentException: Width (0) and height (0) cannot be <= 0}, aborting the whole export at
 *       document assembly (after the paint pass, so the exporter's paint try/catch never sees it).</li>
 *   <li><b>Corruption:</b> normal solid content drawn after a gradient (tip labels, boxes) is rasterized AND tinted
 *       with the stale gradient colour, embedding as raster images instead of clean vector.</li>
 * </ul>
 *
 * <p>The fix records a solid {@link Color} for any gradient up front (via {@code setColor}, never a gradient
 * {@code setPaint}), so the filter is never armed: no rasterization, no crash, no tint, and the domain box exports
 * as a true flat vector rectangle. The mild on-screen domain gradient is subtle, and PDF (OpenPDF, a separate
 * exporter with native shading) and the on-screen / raster paths are unaffected -- only SVG/EPS render the box flat.
 * {@link #fill(Shape)} additionally drops a would-be zero-size fill as a belt-and-suspenders guard.
 */
class GuardedVectorGraphics2D extends VectorGraphics2D {

    @Override
    public void setPaint( final Paint paint ) {
        if ( ( paint != null ) && !( paint instanceof Color ) ) {
            // route through setColor (a SetColorCommand) so the image-rasterization filter is never armed
            setColor( flatten( paint ) );
            return;
        }
        super.setPaint( paint );
    }

    /** A representative solid colour for a paint VectorGraphics2D can't express: the two-stop average of a
     *  {@link GradientPaint} (~ the domain box's base tone), else neutral grey (Archaeopteryx draws no other
     *  non-{@code Color} paint). */
    private static Color flatten( final Paint paint ) {
        if ( paint instanceof GradientPaint ) {
            final Color a = ( (GradientPaint) paint ).getColor1();
            final Color b = ( (GradientPaint) paint ).getColor2();
            return new Color( ( a.getRed() + b.getRed() ) / 2, ( a.getGreen() + b.getGreen() ) / 2,
                              ( a.getBlue() + b.getBlue() ) / 2, ( a.getAlpha() + b.getAlpha() ) / 2 );
        }
        return Color.GRAY;
    }

    @Override
    public void fill( final Shape s ) {
        if ( s == null ) {
            return;
        }
        // Secondary net: a fill whose bounds round to a zero dimension would crash IF it were ever rasterized (it
        // isn't, now that gradients are flattened before they can arm the filter) and is invisible regardless, so
        // drop it rather than risk a BufferedImage(0, ..).
        final Rectangle2D b = s.getBounds2D();
        if ( ( Math.round( b.getWidth() ) <= 0 ) || ( Math.round( b.getHeight() ) <= 0 ) ) {
            return;
        }
        super.fill( s );
    }
}
