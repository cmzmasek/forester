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
 * The control panel's layout row ({@link LayoutIcon}) and phylogram/cladogram row ({@link DisplayTypeIcon}) are
 * icon-only, so the glyph IS the label: if two of them look alike, the control is unusable. These are pure paint
 * tests -- no GUI or look-and-feel needed -- asserting that each glyph paints ink, that every glyph in a row is
 * distinguishable from every other, and that each one actually depicts what it claims:
 * <ul>
 * <li>the three rectangular layout glyphs put their three TIPS on the edge the layout's tips point at (right for
 * root-left, bottom for root-top, top for root-bottom) -- the orientation is the entire meaning of those three;</li>
 * <li>the circular glyph paints a rim all the way round, and is not the theme toggle's sun;</li>
 * <li>the phylogram glyph stops SHORT of the tip column while the cladogram is flush with it, and the aligned
 * phylogram has both -- solid branches short of the column plus a solid rule ON it.</li>
 * </ul>
 */
public final class TreeGlyphIconTest {

    // large enough that a 3-tip glyph has separable tips at 1 px of antialiasing slack; the control panel draws
    // these at ~17 px, and the structural assertions below hold across that range
    private static final int SIZE = 24;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreeGlyphIcon: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return layoutIcons() && displayTypeIcons() && disabledFades();
    }

    /** A DISABLED layout / P-A-C button must fade its glyph (Swing/FlatLaf will not do it for a custom Icon);
     *  these buttons really are disabled in practice -- all of P/A/C for a tree without branch lengths, "A" in
     *  unrooted -- and a full-strength glyph on a dead button looks clickable. */
    private static boolean disabledFades() {
        final JLabel disabled = new JLabel();
        disabled.setForeground( java.awt.Color.BLACK );
        disabled.setBackground( java.awt.Color.WHITE );
        disabled.setEnabled( false );
        final Icon[] icons = { new DisplayTypeIcon( DisplayTypeIcon.Kind.ALIGNED_PHYLOGRAM, SIZE ),
                new LayoutIcon( LayoutIcon.Kind.CIRCULAR, SIZE ) };
        for ( final Icon icon : icons ) {
            final BufferedImage down = new BufferedImage( icon.getIconWidth(), icon.getIconHeight(),
                                                          BufferedImage.TYPE_INT_ARGB );
            icon.paintIcon( disabled, down.getGraphics(), 0, 0 );
            if ( darkestLuma( down ) < 110 ) {
                return fail( icon.getClass().getSimpleName()
                        + " must fade to grey on a disabled black-on-white button, darkest luma "
                        + darkestLuma( down ) );
            }
            if ( ink( down ) < 12 ) {
                return fail( icon.getClass().getSimpleName() + " should still be visible when disabled" );
            }
        }
        return true;
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

    // ---- LayoutIcon ------------------------------------------------------------------------------------

    private static boolean layoutIcons() {
        final LayoutIcon.Kind[] kinds = LayoutIcon.Kind.values();
        if ( kinds.length != 5 ) {
            return fail( "there must be exactly five primary display types, got " + kinds.length );
        }
        final BufferedImage[] imgs = new BufferedImage[ kinds.length ];
        for ( int i = 0; i < kinds.length; i++ ) {
            final LayoutIcon icon = new LayoutIcon( kinds[ i ], SIZE );
            if ( icon.getKind() != kinds[ i ] ) {
                return fail( "LayoutIcon must report the kind it was built with" );
            }
            if ( ( icon.getIconWidth() != SIZE ) || ( icon.getIconHeight() != SIZE ) ) {
                return fail( "LayoutIcon must report the requested size" );
            }
            imgs[ i ] = render( icon );
            if ( ink( imgs[ i ] ) < 12 ) {
                return fail( "the " + kinds[ i ] + " glyph paints almost nothing (" + ink( imgs[ i ] ) + " px)" );
            }
        }
        for ( int i = 0; i < kinds.length; i++ ) {
            for ( int j = i + 1; j < kinds.length; j++ ) {
                if ( diff( imgs[ i ], imgs[ j ] ) < ( ( SIZE * SIZE ) / 12 ) ) {
                    return fail( "the " + kinds[ i ] + " and " + kinds[ j ] + " glyphs look too much alike ("
                            + diff( imgs[ i ], imgs[ j ] ) + " differing px)" );
                }
            }
        }
        // The three rectangular glyphs are the same tree rotated, so each must land its three tips on the edge
        // its layout's tips point at -- that is what tells root-left from root-top from root-bottom at a glance.
        final BufferedImage left = imgs[ LayoutIcon.Kind.ROOT_LEFT.ordinal() ];
        final BufferedImage top = imgs[ LayoutIcon.Kind.ROOT_TOP.ordinal() ];
        final BufferedImage bottom = imgs[ LayoutIcon.Kind.ROOT_BOTTOM.ordinal() ];
        if ( runsOnEdge( left, Edge.RIGHT ) != 3 ) {
            return fail( "root-left must show 3 tips on the RIGHT edge, got " + runsOnEdge( left, Edge.RIGHT ) );
        }
        if ( runsOnEdge( top, Edge.BOTTOM ) != 3 ) {
            return fail( "root-top must show 3 tips on the BOTTOM edge, got " + runsOnEdge( top, Edge.BOTTOM ) );
        }
        if ( runsOnEdge( bottom, Edge.TOP ) != 3 ) {
            return fail( "root-bottom must show 3 tips on the TOP edge, got " + runsOnEdge( bottom, Edge.TOP ) );
        }
        // ... and the root end must be a single stub on the opposite edge, not a spray of tips
        if ( runsOnEdge( left, Edge.LEFT ) != 1 ) {
            return fail( "root-left must show a single root stub on the LEFT edge, got "
                    + runsOnEdge( left, Edge.LEFT ) );
        }
        if ( runsOnEdge( top, Edge.TOP ) != 1 ) {
            return fail( "root-top must show a single root stub on the TOP edge, got " + runsOnEdge( top, Edge.TOP ) );
        }
        if ( runsOnEdge( bottom, Edge.BOTTOM ) != 1 ) {
            return fail( "root-bottom must show a single root stub on the BOTTOM edge, got "
                    + runsOnEdge( bottom, Edge.BOTTOM ) );
        }
        // The circular glyph's rim has to actually go round -- ink out near the edge on THREE of the four sides.
        // Not four: the rim is deliberately left open on the right, where a real circular tree's root wedge is,
        // and that gap is also part of what keeps it from reading as the sun. (A glyph reaching only two sides
        // would read as a bracket rather than a circle.)
        final BufferedImage circular = imgs[ LayoutIcon.Kind.CIRCULAR.ordinal() ];
        int sides_reached = 0;
        for ( final Edge e : Edge.values() ) {
            if ( inkNearEdge( circular, e ) > 0 ) {
                sides_reached++;
            }
        }
        if ( sides_reached < 3 ) {
            return fail( "the circular glyph's rim must reach at least 3 sides, got " + sides_reached );
        }
        if ( inkNearEdge( circular, Edge.RIGHT ) > 0 ) {
            return fail( "the circular glyph must leave its root wedge open on the right" );
        }
        // ... and it must not be confusable with the sun on the theme toggle, which sits two rows above it
        if ( diff( circular, render( new ThemeToggleIcon( true, SIZE ) ) ) < ( ( SIZE * SIZE ) / 12 ) ) {
            return fail( "the circular glyph must not look like the theme toggle's sun" );
        }
        return true;
    }

    // ---- DisplayTypeIcon -------------------------------------------------------------------------------

    /**
     * The three phylogram/cladogram glyphs encode two independent cues -- are the BRANCHES ragged, and are the
     * LABELS lined up -- so that adjacent buttons differ in exactly one of them:
     * <pre>
     *                        branches   labels
     *   PHYLOGRAM            ragged     ragged
     *   ALIGNED_PHYLOGRAM    ragged     aligned
     *   CLADOGRAM            flush      aligned
     * </pre>
     * Crucially, PHYLOGRAM and ALIGNED_PHYLOGRAM must draw the very SAME tree -- if their branches drifted apart
     * the pair would be telling the reader two things at once and the "only the labels moved" story would be
     * lost. This reads the rendered pixels row by row: in each tip row the last ink run is the label tick and
     * the one before it is the branch.
     */
    private static boolean displayTypeIcons() {
        final DisplayTypeIcon.Kind[] kinds = DisplayTypeIcon.Kind.values();
        final BufferedImage[] imgs = new BufferedImage[ kinds.length ];
        for ( int i = 0; i < kinds.length; i++ ) {
            final DisplayTypeIcon icon = new DisplayTypeIcon( kinds[ i ], SIZE );
            if ( icon.getKind() != kinds[ i ] ) {
                return fail( "DisplayTypeIcon must report the kind it was built with" );
            }
            if ( icon.getIconHeight() != SIZE ) {
                return fail( "DisplayTypeIcon must report the requested height" );
            }
            // a tree PLUS a label column is a wide thing, and these are the widest buttons on the panel
            if ( icon.getIconWidth() <= icon.getIconHeight() ) {
                return fail( "the display-type glyph must be wider than it is tall, got " + icon.getIconWidth()
                        + "x" + icon.getIconHeight() );
            }
            imgs[ i ] = render( icon );
            if ( ink( imgs[ i ] ) < 12 ) {
                return fail( "the " + kinds[ i ] + " glyph paints almost nothing (" + ink( imgs[ i ] ) + " px)" );
            }
        }
        for ( int i = 0; i < kinds.length; i++ ) {
            for ( int j = i + 1; j < kinds.length; j++ ) {
                if ( diff( imgs[ i ], imgs[ j ] ) < ( SIZE / 2 ) ) {
                    return fail( "the " + kinds[ i ] + " and " + kinds[ j ] + " glyphs look too much alike ("
                            + diff( imgs[ i ], imgs[ j ] ) + " differing px)" );
                }
            }
        }
        final int[] rows = DisplayTypeIcon.tipRowsForTest( SIZE );
        final int[][] branch = new int[ kinds.length ][];
        final int[][] label = new int[ kinds.length ][];
        for ( int i = 0; i < kinds.length; i++ ) {
            branch[ i ] = new int[ rows.length ];
            label[ i ] = new int[ rows.length ];
            for ( int r = 0; r < rows.length; r++ ) {
                final int[][] runs = runsNearRow( imgs[ i ], rows[ r ] );
                if ( runs.length < 2 ) {
                    return fail( "the " + kinds[ i ] + " glyph's tip row " + r
                            + " must show a branch AND its label tick, found " + runs.length + " ink run(s)" );
                }
                branch[ i ][ r ] = runs[ runs.length - 2 ][ 1 ]; // where the branch ends
                label[ i ][ r ] = runs[ runs.length - 1 ][ 0 ];  // where the label starts
            }
        }
        final int p = DisplayTypeIcon.Kind.PHYLOGRAM.ordinal();
        final int a = DisplayTypeIcon.Kind.ALIGNED_PHYLOGRAM.ordinal();
        final int c = DisplayTypeIcon.Kind.CLADOGRAM.ordinal();
        // branch lengths: ragged for both phylogram flavors, flush for the cladogram
        if ( allEqual( branch[ p ] ) ) {
            return fail( "the phylogram glyph's branches must be ragged, got " + str( branch[ p ] ) );
        }
        if ( allEqual( branch[ a ] ) ) {
            return fail( "the aligned-phylogram glyph's branches must be ragged, got " + str( branch[ a ] ) );
        }
        if ( !allEqual( branch[ c ] ) ) {
            return fail( "the cladogram glyph's branches must be flush, got " + str( branch[ c ] ) );
        }
        // ... and the two phylogram flavors must draw the SAME tree -- only the labels may move between them
        for ( int r = 0; r < rows.length; r++ ) {
            if ( branch[ p ][ r ] != branch[ a ][ r ] ) {
                return fail( "phylogram and aligned phylogram must draw the identical tree; branch " + r
                        + " ends at " + branch[ p ][ r ] + " vs " + branch[ a ][ r ] );
            }
        }
        // labels: at each own tip for the plain phylogram, in one column for the other two
        if ( allEqual( label[ p ] ) ) {
            return fail( "the phylogram glyph's labels must sit at their own tips, got " + str( label[ p ] ) );
        }
        if ( !allEqual( label[ a ] ) ) {
            return fail( "the aligned-phylogram glyph's labels must line up in one column, got "
                    + str( label[ a ] ) );
        }
        if ( !allEqual( label[ c ] ) ) {
            return fail( "the cladogram glyph's labels must line up in one column, got " + str( label[ c ] ) );
        }
        // every label must clear its own branch, or the tick reads as more branch rather than as a label
        for ( int i = 0; i < kinds.length; i++ ) {
            for ( int r = 0; r < rows.length; r++ ) {
                if ( label[ i ][ r ] <= branch[ i ][ r ] ) {
                    return fail( "the " + kinds[ i ] + " glyph's label " + r + " must start past its branch ("
                            + label[ i ][ r ] + " vs " + branch[ i ][ r ] + ")" );
                }
            }
        }
        return true;
    }

    private static boolean allEqual( final int[] v ) {
        for ( int i = 1; i < v.length; i++ ) {
            if ( v[ i ] != v[ 0 ] ) {
                return false;
            }
        }
        return true;
    }

    private static String str( final int[] v ) {
        return java.util.Arrays.toString( v );
    }

    /** The runs at {@code y}, tolerating the one-pixel snap STROKE_NORMALIZE applies to thin axis-aligned
     *  lines: of the rows y-1..y+1, the one carrying the most runs answers (the exact row wins a tie). */
    private static int[][] runsNearRow( final BufferedImage img, final int y ) {
        int[][] best = runsInRow( img, y );
        for ( final int dy : new int[] { -1, 1 } ) {
            if ( ( ( y + dy ) >= 0 ) && ( ( y + dy ) < img.getHeight() ) ) {
                final int[][] cand = runsInRow( img, y + dy );
                if ( cand.length > best.length ) {
                    best = cand;
                }
            }
        }
        return best;
    }

    /** The contiguous ink runs in one pixel row, each as {startX, endX}. */
    private static int[][] runsInRow( final BufferedImage img, final int y ) {
        final java.util.List<int[]> runs = new java.util.ArrayList<>();
        int start = -1;
        for ( int x = 0; x < img.getWidth(); x++ ) {
            if ( opaque( img, x, y ) ) {
                if ( start < 0 ) {
                    start = x;
                }
            }
            else if ( start >= 0 ) {
                runs.add( new int[] { start, x - 1 } );
                start = -1;
            }
        }
        if ( start >= 0 ) {
            runs.add( new int[] { start, img.getWidth() - 1 } );
        }
        return runs.toArray( new int[ 0 ][] );
    }

    // ---- pixel helpers ---------------------------------------------------------------------------------

    private enum Edge {
        LEFT, RIGHT, TOP, BOTTOM
    }

    /** How many separate ink RUNS touch the given edge -- i.e. how many branch ends arrive there. */
    private static int runsOnEdge( final BufferedImage img, final Edge edge ) {
        final int n = ( ( edge == Edge.LEFT ) || ( edge == Edge.RIGHT ) ) ? img.getHeight() : img.getWidth();
        int runs = 0;
        boolean in_run = false;
        for ( int i = 0; i < n; i++ ) {
            if ( inkAtEdge( img, edge, i ) ) {
                if ( !in_run ) {
                    runs++;
                    in_run = true;
                }
            }
            else {
                in_run = false;
            }
        }
        return runs;
    }

    /** Total ink in the two-pixel band along the given edge. */
    private static int inkNearEdge( final BufferedImage img, final Edge edge ) {
        final int n = ( ( edge == Edge.LEFT ) || ( edge == Edge.RIGHT ) ) ? img.getHeight() : img.getWidth();
        int total = 0;
        for ( int i = 0; i < n; i++ ) {
            if ( inkAtEdge( img, edge, i ) ) {
                total++;
            }
        }
        return total;
    }

    /** Ink anywhere in the outermost two pixels of {@code edge}, at position {@code i} along that edge. */
    private static boolean inkAtEdge( final BufferedImage img, final Edge edge, final int i ) {
        for ( int d = 0; d < 2; d++ ) {
            switch ( edge ) {
                case LEFT:
                    if ( opaque( img, d, i ) ) {
                        return true;
                    }
                    break;
                case RIGHT:
                    if ( opaque( img, img.getWidth() - 1 - d, i ) ) {
                        return true;
                    }
                    break;
                case TOP:
                    if ( opaque( img, i, d ) ) {
                        return true;
                    }
                    break;
                default:
                    if ( opaque( img, i, img.getHeight() - 1 - d ) ) {
                        return true;
                    }
                    break;
            }
        }
        return false;
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
        System.out.println( "  [TreeGlyphIconTest] " + msg );
        return false;
    }

    private TreeGlyphIconTest() {
    }
}
