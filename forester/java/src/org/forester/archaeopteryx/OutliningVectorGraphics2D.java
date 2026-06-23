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

import java.awt.Font;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.font.TextLayout;
import java.awt.geom.AffineTransform;
import java.text.AttributedCharacterIterator;

import de.erichseifert.vectorgraphics2d.VectorGraphics2D;

/**
 * A {@link VectorGraphics2D} that renders every text string as filled glyph OUTLINES instead of
 * font-referenced text.
 *
 * <p>VectorGraphics2D embeds no fonts: it emits font-family names and lets the viewer supply the glyphs.
 * For a non-standard family (Archaeopteryx's bundled Source Sans 3) the viewer substitutes -- EPS /
 * PostScript falls back to a Times-family serif, SVG to a generic sans -- and italic faces are lost or
 * mis-sized, so the exported figure no longer matches the screen. Outlining bakes the exact glyph
 * geometry (family, weight, slant, advance) into the vector file, so SVG/EPS reproduce the on-screen
 * figure regardless of the viewer's installed fonts -- the WYSIWYG guarantee {@link VectorGraphicsExporter}
 * promises.
 *
 * <p>Glyphs are laid out with a fractional {@link FontRenderContext} so an outline's advance matches the
 * width Archaeopteryx measured to place the next label (via {@code fractionalAdvanceWidth}); without this
 * the per-glyph integer rounding would accumulate into label overlap on the vector path.
 *
 * <p>Cost: larger files and text that is paths, not selectable characters -- the right trade for a
 * non-embedding backend. The PDF path (OpenPDF, which embeds fonts) is unaffected and keeps selectable
 * text. The on-screen paint path never calls {@code create()}, so a single instance covers a whole export.
 */
final class OutliningVectorGraphics2D extends VectorGraphics2D {

    private static final FontRenderContext FRC = new FontRenderContext( null, true, true );

    @Override
    public void drawString( final String str, final float x, final float y ) {
        if ( ( str == null ) || str.isEmpty() ) {
            return;
        }
        final Font font = getFont();
        final GlyphVector gv = font.createGlyphVector( FRC, str );
        fill( gv.getOutline( x, y ) );
    }

    @Override
    public void drawString( final String str, final int x, final int y ) {
        drawString( str, (float) x, (float) y );
    }

    @Override
    public void drawString( final AttributedCharacterIterator iterator, final float x, final float y ) {
        if ( ( iterator == null ) || ( iterator.getBeginIndex() == iterator.getEndIndex() ) ) {
            return;
        }
        // Honors per-run attributes (font, posture) carried by the iterator.
        final TextLayout tl = new TextLayout( iterator, FRC );
        fill( tl.getOutline( AffineTransform.getTranslateInstance( x, y ) ) );
    }

    @Override
    public void drawString( final AttributedCharacterIterator iterator, final int x, final int y ) {
        drawString( iterator, (float) x, (float) y );
    }

    @Override
    public void drawChars( final char[] data, final int offset, final int length, final int x, final int y ) {
        if ( ( data == null ) || ( length <= 0 ) ) {
            return;
        }
        drawString( new String( data, offset, length ), (float) x, (float) y );
    }

    @Override
    public void drawGlyphVector( final GlyphVector gv, final float x, final float y ) {
        if ( gv == null ) {
            return;
        }
        fill( gv.getOutline( x, y ) );
    }
}
