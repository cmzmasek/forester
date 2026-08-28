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
import java.awt.image.BufferedImage;

import javax.swing.JLabel;

/**
 * The control panel's single light/dark switch is icon-only, so the glyph IS the affordance: the sun and the
 * moon must both paint visible ink, must be clearly DIFFERENT from each other, and each must be recognizable --
 * the sun symmetric about its center with ink out at the rim (the rays), the crescent moon strongly asymmetric
 * (all its ink piled on one side) and hollow at the center of its "bite". Pure paint tests -- no GUI or
 * look-and-feel needed -- so this runs in the headless suite.
 */
public final class ThemeToggleIconTest {

    private static final int SIZE = 16;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ThemeToggleIcon: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final ThemeToggleIcon sun = new ThemeToggleIcon( true, SIZE );
        final ThemeToggleIcon moon = new ThemeToggleIcon( false, SIZE );
        // the two flavors must identify themselves consistently (ControlPanel reads isMoon() to know what the
        // button is currently offering)
        if ( !sun.isSun() || sun.isMoon() || moon.isMoon() == false || moon.isSun() ) {
            return fail( "isSun()/isMoon() must be exact opposites and match the constructor flag" );
        }
        if ( ( sun.getIconWidth() != SIZE ) || ( sun.getIconHeight() != SIZE ) ) {
            return fail( "the icon must report the requested size, got " + sun.getIconWidth() + "x"
                    + sun.getIconHeight() );
        }
        final BufferedImage sun_img = render( sun );
        final BufferedImage moon_img = render( moon );
        if ( ink( sun_img ) < 10 ) {
            return fail( "the sun icon must paint visible ink, got " + ink( sun_img ) + " px" );
        }
        if ( ink( moon_img ) < 10 ) {
            return fail( "the moon icon must paint visible ink, got " + ink( moon_img ) + " px" );
        }
        // the whole point of one button instead of two: the two states must be unmistakably different
        if ( diff( sun_img, moon_img ) < ( ( SIZE * SIZE ) / 8 ) ) {
            return fail( "the sun and moon glyphs must look clearly different, differing px = "
                    + diff( sun_img, moon_img ) );
        }
        // the sun is a disc with rays: ink reaches the rim on BOTH sides, and is left/right balanced
        if ( ( leftInk( sun_img ) == 0 ) || ( rightInk( sun_img ) == 0 ) ) {
            return fail( "the sun must paint on both sides of its center" );
        }
        final int sun_bias = Math.abs( leftInk( sun_img ) - rightInk( sun_img ) );
        if ( sun_bias > ( ink( sun_img ) / 4 ) ) {
            return fail( "the sun must be roughly symmetric about its center, left/right imbalance = " + sun_bias );
        }
        // the crescent is strongly one-sided -- that asymmetry is what makes it read as a moon and not a blob
        final int moon_bias = Math.abs( leftInk( moon_img ) - rightInk( moon_img ) );
        if ( moon_bias < ( ink( moon_img ) / 3 ) ) {
            return fail( "the crescent must be clearly lopsided, left/right imbalance = " + moon_bias );
        }
        // ... and the subtracted disc must really be carved out: the icon's own center pixel stays empty
        if ( opaque( moon_img, moon_img.getWidth() / 2, moon_img.getHeight() / 2 ) ) {
            return fail( "the moon's center must be hollow (the crescent's bite), not filled" );
        }
        // the glyph follows the host component's foreground, so it tracks the FlatLaf light/dark theme
        final JLabel red = new JLabel();
        red.setForeground( Color.RED );
        if ( !paintsIn( moon, red, Color.RED ) || !paintsIn( sun, red, Color.RED ) ) {
            return fail( "the icons must paint in the host component's foreground color" );
        }
        // a tiny requested size must still produce a drawable icon (the size is clamped, never zero/negative)
        final ThemeToggleIcon tiny = new ThemeToggleIcon( false, 0 );
        if ( ( tiny.getIconWidth() < 1 ) || ( ink( render( tiny ) ) < 1 ) ) {
            return fail( "a degenerate size must still yield a paintable icon" );
        }
        return true;
    }

    private static BufferedImage render( final ThemeToggleIcon icon ) {
        return render( icon, new JLabel() );
    }

    private static BufferedImage render( final ThemeToggleIcon icon, final JLabel host ) {
        final BufferedImage img = new BufferedImage( icon.getIconWidth(), icon.getIconHeight(),
                                                     BufferedImage.TYPE_INT_ARGB );
        icon.paintIcon( host, img.getGraphics(), 0, 0 );
        return img;
    }

    private static boolean paintsIn( final ThemeToggleIcon icon, final JLabel host, final Color want ) {
        final BufferedImage img = render( icon, host );
        for ( int y = 0; y < img.getHeight(); y++ ) {
            for ( int x = 0; x < img.getWidth(); x++ ) {
                final int rgb = img.getRGB( x, y );
                if ( ( ( ( rgb >> 24 ) & 0xff ) > 200 ) && ( ( rgb & 0xffffff ) == ( want.getRGB() & 0xffffff ) ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    private static boolean opaque( final BufferedImage img, final int x, final int y ) {
        return ( ( img.getRGB( x, y ) >> 24 ) & 0xff ) > 20;
    }

    private static int ink( final BufferedImage img ) {
        return countInk( img, 0, img.getWidth() );
    }

    private static int leftInk( final BufferedImage img ) {
        return countInk( img, 0, img.getWidth() / 2 );
    }

    private static int rightInk( final BufferedImage img ) {
        return countInk( img, ( img.getWidth() + 1 ) / 2, img.getWidth() );
    }

    private static int countInk( final BufferedImage img, final int from_x, final int to_x ) {
        int n = 0;
        for ( int y = 0; y < img.getHeight(); y++ ) {
            for ( int x = from_x; x < to_x; x++ ) {
                if ( opaque( img, x, y ) ) {
                    n++;
                }
            }
        }
        return n;
    }

    private static int diff( final BufferedImage a, final BufferedImage b ) {
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
        System.out.println( "  [ThemeToggleIconTest] " + msg );
        return false;
    }

    private ThemeToggleIconTest() {
    }
}
