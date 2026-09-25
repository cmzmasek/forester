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
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.Map;

/**
 * The INK box of a single character -- how far its drawn shape actually extends above and below the baseline, and
 * how wide it is -- measured once per (font, character) and cached.
 * <p>
 * <b>Why not FontMetrics.</b> {@code FontMetrics.getAscent()/getDescent()/getHeight()} describe the font's LINE box:
 * the room the face reserves for any character, the same numbers for {@code A} as for {@code g}. A sequence logo
 * scales every letter to a height of its own, so the divisor has to be the letter's own ink. Use the line box and
 * every letter comes out slightly short, a letter with a descender badly so -- and the baseline lands in the wrong
 * place, leaving {@code G}, {@code Q} and {@code J} hanging through whatever is drawn beneath the stack. (The
 * archaeopteryx.js session lost an afternoon to the same distinction in its own toolkit, where
 * {@code getBoundingClientRect} on an SVG text node likewise returns the layout box.)
 * <p>
 * A glyph with no ink at all -- a space, or a character the face cannot render -- reports zero size; callers skip it
 * rather than dividing by zero.
 */
final class MsaLogoGlyphs {

    /** Antialiased, no fractional metrics: a measuring context, never used to draw. */
    private final static FontRenderContext FRC = new FontRenderContext( null, true, false );

    private final static Map<Font, Map<Character, float[]>> CACHE = new HashMap<Font, Map<Character, float[]>>();

    /**
     * {@code { advance, inkAscent, inkDescent }} of {@code ch} in {@code font}, in pixels -- deliberately MIXED, and
     * the mix is the contract.
     * <ul>
     * <li>{@code advance} is the layout width, what the next character would be offset by. The horizontal scale uses
     * it so that scaling a monospace glyph by {@code columnWidth / advance} lands its advance box exactly on the
     * column, leaving the face's own side bearings as the gap between neighbouring letters. Scaling by the INK width
     * instead would butt every letter against its neighbours and stretch a narrow character like {@code I} to the
     * full column.</li>
     * <li>{@code inkAscent} / {@code inkDescent} come from the glyph's VISUAL bounds -- the ink above and below the
     * baseline -- because the vertical scale has to make this letter's own shape a given number of pixels tall.</li>
     * </ul>
     * Both {@code inkAscent} and {@code inkDescent} are 0 for a glyph that paints nothing (a space, or a character
     * the face cannot render); {@code advance} may still be positive.
     */
    static float[] metrics( final Font font, final char ch ) {
        Map<Character, float[]> per_font;
        synchronized ( CACHE ) {
            per_font = CACHE.get( font );
            if ( per_font == null ) {
                per_font = new HashMap<Character, float[]>();
                CACHE.put( font, per_font );
            }
            final float[] hit = per_font.get( Character.valueOf( ch ) );
            if ( hit != null ) {
                return hit;
            }
        }
        final GlyphVector gv = font.createGlyphVector( FRC, new char[] { ch } );
        final float advance = gv.getGlyphMetrics( 0 ).getAdvance();
        final Rectangle2D ink = gv.getVisualBounds();
        final float[] box;
        if ( ( ink == null ) || ( ink.getWidth() <= 0 ) || ( ink.getHeight() <= 0 ) ) {
            box = new float[] { Math.max( 0, advance ), 0, 0 };
        }
        else {
            // getVisualBounds is relative to the BASELINE, with y increasing downward: a box from y=-7 to y=2 is
            // 7 px of ink above the baseline and 2 below.
            box = new float[] { Math.max( 0, advance ), (float) Math.max( 0, -ink.getY() ),
                                (float) Math.max( 0, ink.getMaxY() ) };
        }
        synchronized ( CACHE ) {
            final Map<Character, float[]> m = CACHE.get( font );
            if ( m != null ) {
                m.put( Character.valueOf( ch ), box );
            }
        }
        return box;
    }

    /** The full ink height of {@code ch}: what a letter scaled to a given pixel height must be divided by. Never
     *  the font's line height, which is the same for every character and so scales them all slightly wrong. */
    static float inkHeight( final Font font, final char ch ) {
        final float[] b = metrics( font, ch );
        return b[ 1 ] + b[ 2 ];
    }

    /** Test hook: forget everything measured so far, so a test can measure a cold cache. */
    static void clearCacheForTest() {
        synchronized ( CACHE ) {
            CACHE.clear();
        }
    }

    private MsaLogoGlyphs() {
        // not instantiable
    }
}
