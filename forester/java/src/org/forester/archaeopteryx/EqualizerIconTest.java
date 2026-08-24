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

import java.awt.image.BufferedImage;

import javax.swing.JLabel;

/**
 * The menu-bar "processes running" activity indicator paints a small set of equalizer bars that must (a) render
 * visible ink and (b) actually MOVE as the animation advances. This is a pure paint test -- no GUI/look-and-feel
 * needed (the icon falls back to a default accent color) -- so it runs in the headless suite.
 */
public final class EqualizerIconTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "EqualizerIcon: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final EqualizerIcon icon = new EqualizerIcon( 16, 12, 4 );
        final BufferedImage f0 = render( icon );
        if ( inkPixels( f0 ) < 4 ) {
            return fail( "the equalizer icon must paint visible bars" );
        }
        // a half-cycle phase shift must change the bar heights (the animation)
        icon.setPhaseForTest( Math.PI );
        if ( diffPixels( f0, render( icon ) ) < 4 ) {
            return fail( "the bars must differ between phases" );
        }
        // advance() (what the timer calls each tick) must also change the frame
        final EqualizerIcon ticking = new EqualizerIcon( 16, 12, 4 );
        final BufferedImage a0 = render( ticking );
        for ( int i = 0; i < 4; i++ ) {
            ticking.advance();
        }
        if ( diffPixels( a0, render( ticking ) ) < 4 ) {
            return fail( "advance() must move the animation forward" );
        }
        return true;
    }

    private static BufferedImage render( final EqualizerIcon icon ) {
        final BufferedImage img = new BufferedImage( icon.getIconWidth() + 2, icon.getIconHeight() + 2,
                                                     BufferedImage.TYPE_INT_ARGB );
        icon.paintIcon( new JLabel(), img.getGraphics(), 1, 1 );
        return img;
    }

    private static int inkPixels( final BufferedImage im ) {
        int n = 0;
        for ( int y = 0; y < im.getHeight(); y++ ) {
            for ( int x = 0; x < im.getWidth(); x++ ) {
                if ( ( ( im.getRGB( x, y ) >> 24 ) & 0xff ) > 20 ) {
                    n++;
                }
            }
        }
        return n;
    }

    private static int diffPixels( final BufferedImage a, final BufferedImage b ) {
        int n = 0;
        for ( int y = 0; y < a.getHeight(); y++ ) {
            for ( int x = 0; x < a.getWidth(); x++ ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    n++;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [EqualizerIconTest] " + msg );
        return false;
    }

    private EqualizerIconTest() {
    }
}
