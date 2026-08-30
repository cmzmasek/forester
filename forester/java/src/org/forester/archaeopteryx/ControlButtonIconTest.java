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

import javax.swing.Icon;
import javax.swing.JLabel;

/**
 * The control panel's action glyphs replaced bare letters, so the picture now IS the label. Pure paint tests --
 * no GUI needed -- covering the things that would quietly make one of them wrong:
 * <ul>
 * <li>every glyph paints, and no two of the five look alike;</li>
 * <li>the two rotate arcs are true MIRROR images of each other -- clockwise and counter-clockwise are the pair a
 * reader is most easily shown the wrong way round, and nothing else in the figure would betray it;</li>
 * <li>neither rotate arc is confusable with the circular LAYOUT glyph a row above it (the reason they carry no
 * centre hub) nor with the theme toggle's sun;</li>
 * <li>"return to the whole tree" is "one level back" plus the bar it stops against -- strictly more ink, and ink at
 * the left edge -- so the pair reads as all-the-way versus one-step rather than as two unrelated arrows.</li>
 * </ul>
 */
public final class ControlButtonIconTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ControlButtonIcon: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    /**
     * Runs the whole battery at every size the panel can actually ship a glyph at, NOT just at a comfortable
     * one. {@link ControlPanel#glyphIconSize()} scales with the GUI font and bottoms out at 13 px, so a design
     * that only reads at 24 px would sail through a single-size test and be mush in the shipped product -- the
     * exact failure mode this increment kept running into. Thresholds are therefore expressed as a fraction of
     * the glyph's AREA rather than as pixel counts, so they mean the same thing at 13 px as at 24.
     */
    public static boolean test() {
        final int shipped = ControlPanel.glyphIconSize();
        final int[] sizes = ( shipped == GLYPH_SIZE_FLOOR ) ? new int[] { GLYPH_SIZE_FLOOR, 24 }
                : new int[] { GLYPH_SIZE_FLOOR, shipped, 24 };
        for ( final int size : sizes ) {
            if ( !testAtSize( size ) ) {
                System.out.println( "   (at glyph size " + size + " px)" );
                return false;
            }
        }
        return true;
    }

    /** The smallest glyph the panel can ship (ControlPanel.glyphIconSize() clamps here). */
    private static final int GLYPH_SIZE_FLOOR = 13;

    private static boolean testAtSize( final int SIZE ) {
        final ControlButtonIcon.Kind[] kinds = ControlButtonIcon.Kind.values();
        final BufferedImage[] imgs = new BufferedImage[ kinds.length ];
        for ( int i = 0; i < kinds.length; i++ ) {
            final ControlButtonIcon icon = new ControlButtonIcon( kinds[ i ], SIZE );
            if ( icon.getKind() != kinds[ i ] ) {
                return fail( "the icon must report the kind it was built with" );
            }
            if ( ( icon.getIconWidth() != SIZE ) || ( icon.getIconHeight() != SIZE ) ) {
                return fail( kinds[ i ] + " must report the requested size" );
            }
            imgs[ i ] = render( icon );
            if ( ink( imgs[ i ] ) < minimumInk( SIZE ) ) {
                return fail( kinds[ i ] + " paints almost nothing (" + ink( imgs[ i ] ) + " px)" );
            }
        }
        for ( int i = 0; i < kinds.length; i++ ) {
            for ( int j = i + 1; j < kinds.length; j++ ) {
                // The rotate pair is EXEMPT here on purpose: the two are the same ring mirrored, so of course
                // they share most of their ink -- that is the design, not a defect. A generic "no two alike"
                // rule cannot express it, so the pair gets its own, stronger assertions below (they must be
                // true mirrors, mirroring must not be a no-op, and the difference must land where the eye
                // looks). Loosening this threshold for everyone would have hidden a real collision elsewhere.
                if ( isRotatePair( kinds[ i ], kinds[ j ] ) ) {
                    continue;
                }
                if ( diff( imgs[ i ], imgs[ j ] ) < ( ( SIZE * SIZE ) / 12 ) ) {
                    return fail( kinds[ i ] + " and " + kinds[ j ] + " look too much alike ("
                            + diff( imgs[ i ], imgs[ j ] ) + " differing px)" );
                }
            }
        }
        // ---- the rotate pair are mirror images -------------------------------------------------------
        final BufferedImage cw = imgs[ ControlButtonIcon.Kind.ROTATE_CW.ordinal() ];
        final BufferedImage ccw = imgs[ ControlButtonIcon.Kind.ROTATE_CCW.ordinal() ];
        final int mirror_diff = diff( mirrored( cw ), ccw );
        if ( mirror_diff > ( ( SIZE * SIZE ) / 12 ) ) {
            return fail( "clockwise and counter-clockwise must be mirror images; flipping one leaves "
                    + mirror_diff + " differing px" );
        }
        // ... and a mirror is a real difference, not a no-op (a symmetric glyph would pass the check above).
        // The ring itself IS symmetric about the vertical, so all of the asymmetry is the head changing sides --
        // roughly twice the head's own area (~38 px at 24). A 24th of the glyph's area sits well above the
        // antialiasing noise and well below that at EVERY shipped size, so the check bites on a genuinely
        // symmetric glyph without being tuned to the current head dimensions.
        if ( diff( mirrored( cw ), cw ) < ( ( SIZE * SIZE ) / 24 ) ) {
            return fail( "a rotate glyph must not be left-right symmetric, or the two directions are the same" );
        }
        // The two must actually be TELLABLE APART, and where it counts: the pivot and most of the orbit are
        // identical, so every bit of the difference has to be down at the heads, either side of the gap at the
        // BOTTOM of the glyph. If it drifted elsewhere on the ring the pair would be technically distinct and
        // practically indistinguishable -- which is exactly what a wrong-way arc sweep produces.
        final int total_diff = diff( cw, ccw );
        if ( total_diff < ( ( SIZE * SIZE ) / 24 ) ) {
            return fail( "the two rotate directions must be visibly different, only " + total_diff + " px apart" );
        }
        if ( bottomHalfDiff( cw, ccw ) < ( ( total_diff * 3 ) / 4 ) ) {
            return fail( "the difference between the rotate directions must sit at the heads by the gap, not "
                    + "elsewhere on the ring (" + bottomHalfDiff( cw, ccw ) + " of " + total_diff
                    + " px in the bottom half)" );
        }
        // ---- and not confusable with the round glyphs elsewhere on the panel --------------------------
        final BufferedImage circular_layout = render( new LayoutIcon( LayoutIcon.Kind.CIRCULAR, SIZE ) );
        final BufferedImage sun = render( new ThemeToggleIcon( true, SIZE ) );
        for ( final BufferedImage other : new BufferedImage[] { circular_layout, sun } ) {
            if ( diff( cw, other ) < ( ( SIZE * SIZE ) / 12 ) ) {
                return fail( "a rotate glyph must not look like the circular-layout button or the theme sun" );
            }
        }
        // ---- the fit family: the FRAME's shape carries the axis --------------------------------------
        final BufferedImage fit_w = imgs[ ControlButtonIcon.Kind.FIT_WIDTH.ordinal() ];
        final BufferedImage fit_h = imgs[ ControlButtonIcon.Kind.FIT_HEIGHT.ordinal() ];
        final BufferedImage fit_all = imgs[ ControlButtonIcon.Kind.FIT_ALL.ordinal() ];
        if ( aspect( fit_w ) <= 1.15 ) {
            return fail( "fit-width must be a LANDSCAPE frame (clearly wider than tall), aspect " + aspect( fit_w ) );
        }
        if ( aspect( fit_h ) >= 0.87 ) {
            return fail( "fit-height must be a PORTRAIT frame (clearly taller than wide), aspect " + aspect( fit_h ) );
        }
        if ( ( aspect( fit_all ) < 0.9 ) || ( aspect( fit_all ) > 1.1 ) ) {
            return fail( "fit-everything must be a roughly SQUARE frame, aspect " + aspect( fit_all ) );
        }
        // ... and width/height are the same design rotated, so their ink amounts must be close
        if ( Math.abs( ink( fit_w ) - ink( fit_h ) ) > ( ink( fit_w ) / 4 ) ) {
            return fail( "fit-width and fit-height should be the same design turned 90 degrees, got " + ink( fit_w )
                    + " vs " + ink( fit_h ) + " px" );
        }

        // ---- a DISABLED button must fade its glyph -----------------------------------------------------
        // FlatLaf synthesizes disabled icons only for ImageIcons; for a custom Icon it paints the normal one
        // full-strength, so without inkFor() a disabled R/R1 at startup looks exactly like an enabled button.
        // Black-on-white host: the enabled glyph must contain near-black ink, the disabled one must contain NONE
        // (a solid mid grey), yet must still be visible.
        final JLabel disabled_host = new JLabel();
        disabled_host.setForeground( Color.BLACK );
        disabled_host.setBackground( Color.WHITE );
        disabled_host.setEnabled( false );
        final BufferedImage faded = new BufferedImage( SIZE, SIZE, BufferedImage.TYPE_INT_ARGB );
        new ControlButtonIcon( ControlButtonIcon.Kind.WHOLE_TREE, SIZE ).paintIcon( disabled_host,
                                                                                   faded.getGraphics(), 0, 0 );
        if ( darkestLuma( imgs[ ControlButtonIcon.Kind.WHOLE_TREE.ordinal() ] ) > 60 ) {
            return fail( "the enabled glyph should paint near-black on this host" );
        }
        if ( darkestLuma( faded ) < 110 ) {
            return fail( "a disabled button's glyph must fade -- its darkest ink is still luma "
                    + darkestLuma( faded ) );
        }
        if ( ink( faded ) < minimumInk( SIZE ) ) {
            return fail( "a disabled glyph should still be visible, got " + ink( faded ) + " px of ink" );
        }

        // ---- whole-tree is up-one-level plus the bar it stops against --------------------------------
        final BufferedImage whole = imgs[ ControlButtonIcon.Kind.WHOLE_TREE.ordinal() ];
        final BufferedImage one = imgs[ ControlButtonIcon.Kind.UP_ONE_LEVEL.ordinal() ];
        if ( ink( whole ) <= ink( one ) ) {
            return fail( "'whole tree' should add the stop bar to the plain arrow, got " + ink( whole ) + " vs "
                    + ink( one ) + " px" );
        }
        if ( leftEdgeInk( whole ) <= leftEdgeInk( one ) ) {
            return fail( "'whole tree' must draw the bar it stops against, at the left edge of the glyph" );
        }
        return true;
    }

    /** A glyph has to paint SOMETHING; scales with the area so it means the same at 13 px as at 24. */
    private static int minimumInk( final int size ) {
        return Math.max( 4, ( size * size ) / 48 );
    }

    private static boolean isRotatePair( final ControlButtonIcon.Kind a, final ControlButtonIcon.Kind b ) {
        return ( ( a == ControlButtonIcon.Kind.ROTATE_CW ) && ( b == ControlButtonIcon.Kind.ROTATE_CCW ) )
                || ( ( a == ControlButtonIcon.Kind.ROTATE_CCW ) && ( b == ControlButtonIcon.Kind.ROTATE_CW ) );
    }

    /** Differing pixels in the bottom half only (where the rotate heads and their gap live). */
    private static int bottomHalfDiff( final BufferedImage a, final BufferedImage b ) {
        int n = 0;
        for ( int y = a.getHeight() / 2; y < a.getHeight(); y++ ) {
            for ( int x = 0; x < a.getWidth(); x++ ) {
                if ( opaque( a, x, y ) != opaque( b, x, y ) ) {
                    n++;
                }
            }
        }
        return n;
    }

    /** Width over height of the glyph's ink bounding box. */
    private static double aspect( final BufferedImage img ) {
        int minx = Integer.MAX_VALUE, miny = Integer.MAX_VALUE, maxx = -1, maxy = -1;
        for ( int y = 0; y < img.getHeight(); y++ ) {
            for ( int x = 0; x < img.getWidth(); x++ ) {
                if ( opaque( img, x, y ) ) {
                    minx = Math.min( minx, x );
                    miny = Math.min( miny, y );
                    maxx = Math.max( maxx, x );
                    maxy = Math.max( maxy, y );
                }
            }
        }
        return ( maxy < miny ) ? 1.0 : ( ( maxx - minx ) + 1 ) / ( double ) ( ( maxy - miny ) + 1 );
    }

    /** The darkest opaque pixel's luminance -- near 0 for full-strength black ink, well above for faded grey. */
    private static int darkestLuma( final BufferedImage img ) {
        int darkest = 255;
        for ( int y = 0; y < img.getHeight(); y++ ) {
            for ( int x = 0; x < img.getWidth(); x++ ) {
                if ( opaque( img, x, y ) ) {
                    final int rgb = img.getRGB( x, y );
                    final int luma = ( ( ( rgb >> 16 ) & 0xff ) + ( ( rgb >> 8 ) & 0xff ) + ( rgb & 0xff ) ) / 3;
                    darkest = Math.min( darkest, luma );
                }
            }
        }
        return darkest;
    }

    /** The image flipped left-to-right. */
    private static BufferedImage mirrored( final BufferedImage src ) {
        final BufferedImage out = new BufferedImage( src.getWidth(), src.getHeight(), BufferedImage.TYPE_INT_ARGB );
        for ( int y = 0; y < src.getHeight(); y++ ) {
            for ( int x = 0; x < src.getWidth(); x++ ) {
                out.setRGB( ( src.getWidth() - 1 ) - x, y, src.getRGB( x, y ) );
            }
        }
        return out;
    }

    /** Ink in the left fifth of the glyph, where the "all the way back" stop bar sits. */
    private static int leftEdgeInk( final BufferedImage img ) {
        int n = 0;
        for ( int y = 0; y < img.getHeight(); y++ ) {
            for ( int x = 0; x < ( img.getWidth() / 5 ); x++ ) {
                if ( opaque( img, x, y ) ) {
                    n++;
                }
            }
        }
        return n;
    }

    private static BufferedImage render( final Icon icon ) {
        final BufferedImage img = new BufferedImage( icon.getIconWidth(), icon.getIconHeight(),
                                                     BufferedImage.TYPE_INT_ARGB );
        final JLabel host = new JLabel();
        host.setForeground( Color.BLACK );
        icon.paintIcon( host, img.getGraphics(), 0, 0 );
        return img;
    }

    private static boolean opaque( final BufferedImage img, final int x, final int y ) {
        return ( ( img.getRGB( x, y ) >> 24 ) & 0xff ) > 40;
    }

    private static int ink( final BufferedImage img ) {
        int n = 0;
        for ( int y = 0; y < img.getHeight(); y++ ) {
            for ( int x = 0; x < img.getWidth(); x++ ) {
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
                if ( opaque( a, x, y ) != opaque( b, x, y ) ) {
                    n++;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ControlButtonIconTest] " + msg );
        return false;
    }

    private ControlButtonIconTest() {
    }
}
