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
import java.awt.Graphics2D;
import java.awt.font.FontRenderContext;
import java.awt.image.BufferedImage;

/**
 * Tests for {@link TreePanel#fractionalAdvanceWidth(Graphics2D, String)}, the helper that positions
 * the label following a taxonomy/scientific name on the vector-export path.
 *
 * <p>Background (the bug this guards): on screen and in raster exports AWT grid-fits each glyph
 * advance to an integer, so {@code FontMetrics.stringWidth} matches what is drawn. The vector
 * backends (PDF via OpenPDF glyph outlines, SVG/EPS via VectorGraphics2D {@code <text>}) advance
 * glyphs by the font's true *fractional* widths instead; at the small leaf-label font the per-glyph
 * rounding deficit accumulates over a long label until the next label (placed at the integer width)
 * overlaps it by a character or two. The fix measures the advance with a fractional
 * {@link FontRenderContext}; this test pins that contract so a revert to integer
 * {@code stringWidth} is caught, and confirms the fractional advance is not narrower than the
 * integer one at a small font (the direction that prevents the overlap).
 *
 * <p>Pure font math on an offscreen {@link BufferedImage} graphics -- no display -- so it runs in
 * the headless suite.
 */
public final class LabelAdvanceWidthTest {

    // matches the constant inside TreePanel (private there); fractional metrics on.
    private static final FontRenderContext FRC = new FontRenderContext( null, true, true );

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LabelAdvanceWidth: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        Graphics2D g = null;
        try {
            g = new BufferedImage( 8, 8, BufferedImage.TYPE_INT_RGB ).createGraphics();

            // empty string advances nothing
            g.setFont( new Font( Font.SANS_SERIF, Font.PLAIN, 6 ) );
            if ( TreePanel.fractionalAdvanceWidth( g, "" ) != 0 ) {
                return fail( "empty string must advance 0" );
            }

            // contract: the method returns ceil() of the font's *fractional* advance for g's current
            // font -- NOT FontMetrics.stringWidth. Checked across sizes/styles and a realistic label.
            final String[] strings = { "A", "flavivirus ", "Ochlerotatus scapularis flavivirus ",
                                       "Orthoflavivirus_encephalitidis|KT224355.1" };
            final Font[] fonts = { new Font( Font.SANS_SERIF, Font.PLAIN, 6 ),
                                   new Font( Font.SANS_SERIF, Font.PLAIN, 12 ),
                                   new Font( Font.SERIF, Font.BOLD, 6 ),
                                   new Font( Font.MONOSPACED, Font.ITALIC, 9 ) };
            for( final Font f : fonts ) {
                g.setFont( f );
                for( final String s : strings ) {
                    final int got = TreePanel.fractionalAdvanceWidth( g, s );
                    final int expected = (int) Math.ceil( f.getStringBounds( s, FRC ).getWidth() );
                    if ( got != expected ) {
                        return fail( "advance for [" + s + "] @ " + f.getSize() + "px: got " + got + ", expected "
                                + expected + " (ceil of fractional bounds) -- did it revert to stringWidth()?" );
                    }
                    if ( ( s.length() > 0 ) && ( got <= 0 ) ) {
                        return fail( "non-empty string must advance > 0: [" + s + "]" );
                    }
                }
            }

            // the overlap-preventing direction: at the small label font, the fractional advance for a
            // long label is at least the integer stringWidth, so the next label is never placed left of
            // where the (wider) vector glyphs actually end.
            final Font small = new Font( Font.SANS_SERIF, Font.PLAIN, 6 );
            g.setFont( small );
            final String label = "Ochlerotatus scapularis flavivirus Ochlerotatus_scapularis_flavivirus";
            final int frac = TreePanel.fractionalAdvanceWidth( g, label );
            final int integer = g.getFontMetrics().stringWidth( label );
            if ( frac < integer ) {
                return fail( "fractional advance (" + frac + ") < integer stringWidth (" + integer
                        + ") for a long small-font label -- would still let the next label overlap" );
            }

            // monotonicity sanity: a longer string never advances less than its prefix
            if ( TreePanel.fractionalAdvanceWidth( g, "flavivirus" ) > TreePanel
                    .fractionalAdvanceWidth( g, "flavivirus virus" ) ) {
                return fail( "advance should be monotonic in string length" );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( e.toString() );
        }
        finally {
            if ( g != null ) {
                g.dispose();
            }
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "\nLabelAdvanceWidthTest failed: " + msg );
        return false;
    }

    private LabelAdvanceWidthTest() {
    }
}
