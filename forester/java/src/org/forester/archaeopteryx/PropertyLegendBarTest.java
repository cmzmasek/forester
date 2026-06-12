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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;

/**
 * Unit tests for {@link TreePanel#paintGradientBar}: the "Color by:" gradient legend must render at
 * full saturation regardless of the {@link Graphics2D}'s incoming stroke. The bug it guards against:
 * the legend used to inherit the zoom-dependent branch stroke (set by {@code setupStroke}); when the
 * tree was zoomed out that stroke shrank to sub-pixel widths, so the bar -- then drawn as 1px-wide
 * {@code drawLine}s -- was antialiased into a washed-out version of itself. Painting each column with
 * {@code fillRect} makes it stroke- and antialias-independent. In {@code org.forester.archaeopteryx}
 * so it can reach the package-private method; run via the suite ({@link #test()}).
 */
public final class PropertyLegendBarTest {

    public static void main( final String[] args ) {
        System.out.println( "PropertyLegendBar: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        final Color bar_color = new Color( 0, 128, 255 ); // a saturated, clearly-non-background color
        final int w = 40;
        final int h = 12;
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        // start from a black background, antialiasing ON, and -- crucially -- a deliberately
        // sub-pixel stroke, exactly the state a zoomed-out tree leaves on the Graphics2D
        g.setColor( Color.BLACK );
        g.fillRect( 0, 0, w, h );
        g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        g.setStroke( new BasicStroke( 0.05f ) );
        TreePanel.paintGradientBar( g, 0, 0, w, h, t -> bar_color );
        g.dispose();
        // every column must be the full color, not blended toward the black background
        for ( int x = 0; x < w; ++x ) {
            final Color got = new Color( img.getRGB( x, h / 2 ) );
            if ( !closeTo( got, bar_color, 6 ) ) {
                return fail( "column " + x + " washed out: expected " + bar_color + " but got " + got
                        + " (a thin inherited stroke is blending the bar into the background)" );
            }
        }
        // a left-to-right gradient must actually advance: left end != right end, and t spans [0,1]
        final BufferedImage grad = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g2 = grad.createGraphics();
        TreePanel.paintGradientBar( g2, 0, 0, w, h, t -> new Color( (int) Math.round( 255 * t ), 0, 0 ) );
        g2.dispose();
        final Color left = new Color( grad.getRGB( 0, h / 2 ) );
        final Color right = new Color( grad.getRGB( w - 1, h / 2 ) );
        if ( left.getRed() != 0 ) {
            return fail( "gradient left end should be t=0 (red 0) but was " + left.getRed() );
        }
        if ( right.getRed() != 255 ) {
            return fail( "gradient right end should be t=1 (red 255) but was " + right.getRed() );
        }
        // a single-column bar must not divide by zero and uses t=0
        final BufferedImage one = new BufferedImage( 1, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g3 = one.createGraphics();
        TreePanel.paintGradientBar( g3, 0, 0, 1, h, t -> ( t == 0.0 ) ? bar_color : Color.BLACK );
        g3.dispose();
        if ( !closeTo( new Color( one.getRGB( 0, h / 2 ) ), bar_color, 6 ) ) {
            return fail( "single-column bar should use t=0" );
        }
        return true;
    }

    private static boolean closeTo( final Color a, final Color b, final int tol ) {
        return ( Math.abs( a.getRed() - b.getRed() ) <= tol ) && ( Math.abs( a.getGreen() - b.getGreen() ) <= tol )
                && ( Math.abs( a.getBlue() - b.getBlue() ) <= tol );
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [PropertyLegendBarTest] " + message );
        return false;
    }

    private PropertyLegendBarTest() {
    }
}
