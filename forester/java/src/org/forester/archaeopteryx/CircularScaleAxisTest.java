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
import java.awt.Point;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The distance scale in the CIRCULAR layout (Christian, 2026-09-30, for both programs): the Scale Axis is a ruler from
 * the root out to the deepest tip along the SEAM (the gap between the last tip and the first), numbered on one side and
 * upright, the unit past its end; the Scale Grid is the concentric distance rings; the Scale bar is the ordinary bar,
 * one unit = radius / tree height. None of the three on a cladogram, with capped long branches, or when time rings are
 * the tree's scale. Measured on the paint (export and screen), against geometry derived independently: the seam from
 * the tips' own angles, the ruler's end and the ring radii from the tree's numbers.
 */
public final class CircularScaleAxisTest {

    private static int       W   = 820;
    private static int       H   = 820;
    private static final int GAP = 5;  // TreePanel.SCALE_AXIS_UNIT_GAP

    public static void main( final String[] args ) {
        System.out.println( test() ? "CircularScaleAxisTest: OK." : "CircularScaleAxisTest: FAILED." );
        System.exit( 0 );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        withDemo( "scale-axis.xml", ( frame, tp, o ) -> {
            prepare( frame, tp, o );
            seamOk( tp, o, ok );
            axisInSeamOk( tp, o, ok );
            uprightOk( tp, o, ok );
            gridOk( tp, o, ok );
            barOk( tp, o, ok );
            backdropOk( tp, ok );
            // a cladogram states no distance: no bar, axis or grid
            tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM );
            if ( ( diff( render( tp, o, false, false, false ), render( tp, o, true, true, true ) ) ).size() != 0 ) {
                fail( ok, "a circular CLADOGRAM must draw no distance scale" );
            }
        }, ok );
        // a root with a branch of its own: the tree height (what a radius is normalised by) then exceeds the deepest
        // tip's distance, so the ruler and the rings must stop short of the ring, at the deepest tip
        withDemo( "scale-axis.xml", phy -> phy.getRoot().setDistanceToParent( 0.25 ), ( frame, tp, o ) -> {
            prepare( frame, tp, o );
            final double[] n = tp.circularScaleNumbersForTest();
            if ( !( n[ 1 ] < ( n[ 0 ] - ( 2 * n[ 2 ] ) ) ) ) {
                fail( ok, "precondition: the root branch must put the deepest tip two scale distances inside the "
                        + "tree height, got " + java.util.Arrays.toString( n ) );
            }
            axisInSeamOk( tp, o, ok );
            gridOk( tp, o, ok );
        }, ok );
        // a small canvas: the ticks crowd, and the numbers must be thinned so none overlap
        withDemo( "scale-axis.xml", ( frame, tp, o ) -> {
            W = 300;
            H = 300;
            try {
                prepare( frame, tp, o );
                thinningOk( tp, o, ok );
            }
            finally {
                W = 820;
                H = 820;
            }
        }, ok );
        withDemo( "long-branch-break.xml", ( frame, tp, o ) -> {
            prepare( frame, tp, o );
            // precondition: uncapped, the axis draws
            if ( diff( render( tp, o, false, false, false ), render( tp, o, true, false, false ) ).isEmpty() ) {
                fail( ok, "precondition: the long-branch tree must draw its axis while uncapped" );
            }
            o.setBreakLongBranches( true );
            if ( !tp.breakLongBranchesActiveCircular() ) {
                fail( ok, "precondition: Break Long Branches must cap this circular tree" );
            }
            else {
                // no linear ruler spans a capped branch: no axis, no grid
                if ( !diff( render( tp, o, false, false, false ), render( tp, o, true, true, false ), 4 ).isEmpty() ) {
                    fail( ok, "a CAPPED circular tree draws no Scale Axis or Grid" );
                }
                // the bar stays, sized to the unbroken scale -- radius / the CAPPED height per unit (as the
                // rectangular bar stays when capped)
                final double[] n = tp.circularScaleNumbersForTest();
                if ( !( n[ 3 ] < ( n[ 0 ] * 0.8 ) ) ) {
                    fail( ok, "precondition: capping must shrink the radial normaliser, got " + n[ 3 ] + " vs " + n[ 0 ] );
                }
                checkBar( "capped", diff( render( tp, o, false, false, false ), render( tp, o, false, false, true ) ),
                          ( n[ 2 ] / n[ 3 ] ) * tp.circularRadiusForTest(), W, H, ok );
            }
            o.setBreakLongBranches( false );
        }, ok );
        withDemo( "dinosaur-time-tree.xml", ( frame, tp, o ) -> {
            prepare( frame, tp, o );
            if ( !tp.geologicRingsApplyCircular() && !tp.calendarRingsApplyCircular() ) {
                fail( ok, "precondition: the dinosaur tree must show its time rings in circular" );
            }
            else if ( !diff( render( tp, o, false, false, false ), render( tp, o, true, true, true ) ).isEmpty() ) {
                fail( ok, "time rings ARE a dated circular tree's scale: no distance bar, axis or grid beside them" );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** The seam lies half a row past the last tip and half a row before the first; no tip is nearer. */
    private static void seamOk( final TreePanel tp, final Options o, final boolean[] ok ) {
        for ( final double start : new double[] { tp.getStartingAngle(), 1.0, 4.0 } ) {
            tp.setStartingAngle( start );
            render( tp, o, true, false, false );
            final double seam = tp.circularSeamAngleForTest();
            final List<PhylogenyNode> tips = tp.getPhylogeny().getExternalNodes();
            final double row = ( 2 * Math.PI ) / tips.size();
            double nearest = Double.MAX_VALUE;
            int at_half_row = 0;
            for ( final PhylogenyNode t : tips ) {
                final double d = Math.abs( wrap( tp.circularAngleForTest( t ) - seam ) );
                nearest = Math.min( nearest, d );
                if ( Math.abs( d - ( row / 2 ) ) < 1e-9 ) {
                    ++at_half_row;
                }
            }
            if ( ( Math.abs( nearest - ( row / 2 ) ) > 1e-9 ) || ( at_half_row != 2 ) ) {
                fail( ok, "start " + start + ": the seam must sit half a row from exactly the first and last tips, "
                        + "nearest " + nearest + " vs " + ( row / 2 ) + ", at half a row: " + at_half_row );
            }
        }
    }

    /** All axis ink lies along the seam ray, the ruler reaches the deepest tip, and the unit comes past its end only
     *  when the tree has one. */
    private static void axisInSeamOk( final TreePanel tp, final Options o, final boolean[] ok ) {
        final List<int[]> ink = diff( render( tp, o, false, false, false ), render( tp, o, true, false, false ) );
        final double[] n = tp.circularScaleNumbersForTest();
        final double r_end = ( n[ 1 ] / n[ 0 ] ) * tp.circularRadiusForTest();
        if ( ink.isEmpty() ) {
            fail( ok, "Scale Axis on must draw a ruler in circular" );
            return;
        }
        final double[][] ap = alongPerp( tp, ink );
        double line_end = 0;
        double unit_end = 0;
        int numbers = 0;
        int strays = 0;
        for ( final double[] p : ap ) {
            if ( ( Math.abs( p[ 1 ] ) > 40 ) || ( p[ 0 ] < -12 ) ) {
                ++strays;
            }
            if ( Math.abs( p[ 1 ] ) <= 1.5 ) {
                if ( p[ 0 ] <= ( r_end + 2 ) ) {
                    line_end = Math.max( line_end, p[ 0 ] );
                }
                else {
                    unit_end = Math.max( unit_end, p[ 0 ] );
                }
            }
            else if ( Math.abs( p[ 1 ] ) > 6 ) {
                ++numbers;
            }
        }
        if ( strays > 0 ) {
            fail( ok, strays + " axis pixels lie off the seam ray" );
        }
        if ( Math.abs( line_end - r_end ) > 2 ) {
            fail( ok, "the ruler must end at the deepest tip (r " + r_end + "), ends at " + line_end );
        }
        if ( numbers < 60 ) {
            fail( ok, "the ruler must carry its numbers beside it, got " + numbers + " px" );
        }
        if ( unit_end < ( r_end + GAP + 40 ) ) {
            fail( ok, "the unit [substitutions/site] must follow past the ruler's end, reaches " + unit_end + " vs r "
                    + r_end );
        }
        // no unit: nothing past the end on the ruler's line
        final String unit = tp.getPhylogeny().getDistanceUnit();
        tp.getPhylogeny().setDistanceUnit( null );
        for ( final double[] p : alongPerp( tp, diff( render( tp, o, false, false, false ),
                                                      render( tp, o, true, false, false ) ) ) ) {
            if ( ( Math.abs( p[ 1 ] ) <= 1.5 ) && ( p[ 0 ] > ( r_end + 4 ) ) ) {
                fail( ok, "a tree without a unit must draw nothing past the ruler's end, found at " + p[ 0 ] );
                break;
            }
        }
        tp.getPhylogeny().setDistanceUnit( unit );
    }

    /** On a ruler whose ticks stand closer than a number is wide, every other number gives way: along the ruler the
     *  numbers' ink falls into separate runs, none longer than the widest number. */
    private static void thinningOk( final TreePanel tp, final Options o, final boolean[] ok ) {
        final BufferedImage off = render( tp, o, false, false, false );
        final BufferedImage on = render( tp, o, true, false, false );
        final double[] n = tp.circularScaleNumbersForTest();
        final double step = ( n[ 2 ] / n[ 0 ] ) * tp.circularRadiusForTest();
        final java.awt.FontMetrics fm = tp.getMainPanel().getTreeFontSet().getFontMetricsSmall();
        int widest = 0;
        for ( final double v : TreePanelUtil.scaleAxisTickValues( n[ 1 ], n[ 2 ] ) ) {
            widest = Math.max( widest, fm.stringWidth( TreePanelUtil.formatCompactNumber( v ) ) );
        }
        if ( step >= ( widest + 4 ) ) {
            fail( ok, "precondition: the ticks must stand closer (" + step + " px) than a number is wide (" + widest
                    + ")" );
            return;
        }
        final java.util.TreeSet<Integer> along = new java.util.TreeSet<>();
        // the numbers' own INK: dark in the axis paint (a backdrop also changes pixels, where it covers a branch, but
        // lightens them)
        final List<int[]> ink = new ArrayList<>();
        for ( final int[] p : diff( off, on ) ) {
            if ( dark( on, p[ 0 ], p[ 1 ] ) ) {
                ink.add( p );
            }
        }
        for ( final double[] p : alongPerp( tp, ink ) ) {
            if ( Math.abs( p[ 1 ] ) > 6 ) {
                along.add( (int) Math.round( p[ 0 ] ) );
            }
        }
        int runs = 0;
        int run = 0;
        int longest = 0;
        Integer prev = null;
        for ( final int a : along ) {
            // one number's characters stand up to ~2 px apart; thinned numbers at least SCALE_AXIS_LABEL_GAP (4)
            run = ( ( prev != null ) && ( a <= ( prev + 3 ) ) ) ? ( ( run + a ) - prev ) : 1;
            if ( run == 1 ) {
                ++runs;
            }
            longest = Math.max( longest, run );
            prev = a;
        }
        if ( ( runs < 2 ) || ( longest > ( widest + 2 ) ) ) {
            fail( ok, "crowded numbers must be thinned apart: " + runs + " runs, the longest " + longest
                    + " px against a widest number of " + widest );
        }
    }

    /** The numbers read upright on a ruler pointing right AND on one pointing left: below a horizontal ruler either
     *  way, and the same glyphs, not the same glyphs turned over. */
    private static void uprightOk( final TreePanel tp, final Options o, final boolean[] ok ) {
        final double row = ( 2 * Math.PI ) / tp.getPhylogeny().getNumberOfExternalNodes();
        final BufferedImage[] img = new BufferedImage[ 2 ];
        final Point[] c = new Point[ 2 ];
        for ( int i = 0; i < 2; ++i ) {
            tp.setStartingAngle( ( i * Math.PI ) + ( row / 2 ) ); // seam at 0 (right), then at pi (left)
            final BufferedImage off = render( tp, o, false, false, false );
            img[ i ] = render( tp, o, true, false, false );
            c[ i ] = tp.circularCenterForTest();
            if ( Math.abs( wrap( tp.circularSeamAngleForTest() - ( i * Math.PI ) ) ) > 1e-9 ) {
                fail( ok, "precondition: the seam must point " + ( i == 0 ? "right" : "left" ) );
                return;
            }
            int above = 0;
            int below = 0;
            for ( final int[] p : diff( off, img[ i ] ) ) {
                if ( p[ 1 ] < ( c[ i ].y - 6 ) ) {
                    ++above;
                }
                else if ( p[ 1 ] > ( c[ i ].y + 6 ) ) {
                    ++below;
                }
            }
            if ( ( below < 60 ) || ( above > 0 ) ) {
                fail( ok, "the numbers of a " + ( i == 0 ? "right" : "left" ) + "-pointing ruler belong below it, got "
                        + below + " below, " + above + " above" );
            }
        }
        // each number's patch on the right equals its twin's on the left (upright both), not its twin turned over
        final double[] n = tp.circularScaleNumbersForTest();
        final double radius = tp.circularRadiusForTest();
        long same = 0;
        long turned = 0;
        for ( double v = n[ 2 ]; v <= ( n[ 1 ] + 1e-9 ); v += n[ 2 ] ) {
            final int r = (int) Math.round( ( v / n[ 0 ] ) * radius );
            for ( int dy = 7; dy < 20; ++dy ) {
                for ( int dx = -14; dx <= 14; ++dx ) {
                    final boolean right = dark( img[ 0 ], c[ 0 ].x + r + dx, c[ 0 ].y + dy );
                    same += ( right == dark( img[ 1 ], ( c[ 1 ].x - r ) + dx, c[ 1 ].y + dy ) ) ? 0 : 1;
                    turned += ( right == dark( img[ 1 ], ( c[ 1 ].x - r ) - dx, ( c[ 1 ].y + 26 ) - dy ) ) ? 0 : 1;
                }
            }
        }
        if ( !( ( same * 3 ) < turned ) ) {
            fail( ok, "the left-pointing ruler's numbers must read upright (mismatch " + same + " vs turned over "
                    + turned + ")" );
        }
        tp.setStartingAngle( 0 );
    }

    /** Scale Grid draws a ring at each scale distance out to the deepest tip, and nothing else. */
    private static void gridOk( final TreePanel tp, final Options o, final boolean[] ok ) {
        // the grid is FAINT (the background nudged toward the branch colour), so a finer difference sees it
        final List<int[]> ink = diff( render( tp, o, false, false, false ), render( tp, o, false, true, false ), 4 );
        final double[] n = tp.circularScaleNumbersForTest();
        final double radius = tp.circularRadiusForTest();
        final Point c = tp.circularCenterForTest();
        if ( ink.isEmpty() ) {
            fail( ok, "Scale Grid on must draw rings in circular" );
            return;
        }
        final int rings = (int) Math.floor( ( n[ 1 ] / n[ 2 ] ) + 1e-9 );
        if ( rings < 2 ) {
            fail( ok, "precondition: the demo must span at least two scale distances, got " + rings );
            return;
        }
        final int[] hits = new int[ rings + 1 ];
        int off_ring = 0;
        for ( final int[] p : ink ) {
            final double rr = Math.hypot( p[ 0 ] - c.x, p[ 1 ] - c.y );
            final int k = (int) Math.round( rr / ( ( n[ 2 ] / n[ 0 ] ) * radius ) );
            if ( ( k >= 1 ) && ( k <= rings ) && ( Math.abs( rr - ( ( ( k * n[ 2 ] ) / n[ 0 ] ) * radius ) ) <= 2 ) ) {
                ++hits[ k ];
            }
            else {
                ++off_ring;
            }
        }
        if ( off_ring > ( ink.size() / 50 ) ) {
            fail( ok, off_ring + " of " + ink.size() + " grid pixels lie off the distance rings" );
        }
        for ( int k = 1; k <= rings; ++k ) {
            final double r = ( ( k * n[ 2 ] ) / n[ 0 ] ) * radius;
            if ( hits[ k ] < ( Math.PI * r ) ) { // at least half a circumference of ink
                fail( ok, "ring " + k + " (r " + Math.round( r ) + ") must be drawn, got " + hits[ k ] + " px" );
            }
        }
    }

    /** Scale draws the bar, bottom left, one scale distance long in the circle's own units -- and no rings (it drew
     *  them until 2026-09-30). On screen and in an export. */
    private static void barOk( final TreePanel tp, final Options o, final boolean[] ok ) {
        final double[] n = tp.circularScaleNumbersForTest();
        final List<int[]> ink = diff( render( tp, o, false, false, false ), render( tp, o, false, false, true ) );
        final double expected = ( n[ 2 ] / n[ 3 ] ) * tp.circularRadiusForTest();
        checkBar( "export", ink, expected, W, H, ok );
        final List<int[]> screen = diff( screen( tp, o, false ), screen( tp, o, true ) );
        checkBar( "screen", screen, ( n[ 2 ] / n[ 3 ] ) * tp.circularRadiusForTest(), tp.getVisibleRect().width,
                  tp.getVisibleRect().height, ok );
    }

    private static void checkBar( final String where, final List<int[]> ink, final double expected, final int w,
                                  final int h, final boolean[] ok ) {
        if ( ink.isEmpty() ) {
            fail( ok, where + ": Scale must draw a bar in circular" );
            return;
        }
        // the bar is measured between its two END TICKS: they reach 4 px below the bar line, below the label text
        // (baseline 2 px above the line), so the lowest two rows of the ink hold nothing else. (The longest
        // horizontal run is no measure: on a short bar it is a run inside the label.)
        int bottom = -1;
        for ( final int[] p : ink ) {
            if ( ( p[ 0 ] > ( w / 3 ) ) || ( p[ 1 ] < ( h - 60 ) ) ) {
                fail( ok, where + ": Scale alone draws only the bar at the bottom left, found ink at " + p[ 0 ] + ","
                        + p[ 1 ] );
                return;
            }
            bottom = Math.max( bottom, p[ 1 ] );
        }
        int left = Integer.MAX_VALUE;
        int right = Integer.MIN_VALUE;
        for ( final int[] p : ink ) {
            if ( p[ 1 ] >= ( bottom - 1 ) ) {
                left = Math.min( left, p[ 0 ] );
                right = Math.max( right, p[ 0 ] );
            }
        }
        final int longest = right - left;
        if ( Math.abs( longest - expected ) > 2 ) {
            fail( ok, where + ": the bar must be one scale distance long (" + expected + " px), got " + longest );
        }
    }

    private static void backdropOk( final TreePanel tp, final boolean[] ok ) {
        final Color b = tp.scaleAxisBackdropForTest();
        final Color bg = tp.getTreeColorSet().getBackgroundColor();
        if ( ( b == null ) || ( b.getAlpha() != 217 ) || ( b.getRGB() & 0xFFFFFF ) != ( bg.getRGB() & 0xFFFFFF ) ) {
            fail( ok, "a scale-axis number sits on the canvas colour at 0.85, got " + b );
        }
        tp.setExportTransparentBackground( true );
        if ( tp.scaleAxisBackdropForTest() != null ) {
            fail( ok, "a transparent export draws no backdrop box" );
        }
        tp.setExportTransparentBackground( false );
    }

    // ---- helpers ----

    private static void prepare( final MainFrame frame, final TreePanel tp, final Options o ) {
        o.setGraphicsExportWhiteBackground( true );
        o.setShowOverview( false );
        tp.setOvOn( false );
        tp.setColorByPropertyRef( null );
        tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
        frame.showWhole();
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        tp.setPreferredSize( new java.awt.Dimension( W, H ) );
        tp.setSize( W, H );
    }

    private static BufferedImage render( final TreePanel tp, final Options o, final boolean axis, final boolean grid,
                                         final boolean scale ) {
        o.setShowScaleAxis( axis );
        o.setShowScaleGrid( grid );
        o.setShowScale( scale );
        tp.calcParametersForPainting( W, H );
        return AptxUtil.renderPhylogenyToImage( W, H, tp, o, false, 1, false );
    }

    private static BufferedImage screen( final TreePanel tp, final Options o, final boolean scale ) {
        o.setShowScaleAxis( false );
        o.setShowScaleGrid( false );
        o.setShowScale( scale );
        final BufferedImage img = new BufferedImage( tp.getWidth(), tp.getHeight(), BufferedImage.TYPE_INT_RGB );
        tp.printAll( img.getGraphics() );
        return img;
    }

    /** Pixels that differ clearly between two paints of the same layout. */
    private static List<int[]> diff( final BufferedImage a, final BufferedImage b ) {
        return diff( a, b, 40 );
    }

    private static List<int[]> diff( final BufferedImage a, final BufferedImage b, final int by ) {
        final List<int[]> out = new ArrayList<>();
        for ( int y = 0; y < Math.min( a.getHeight(), b.getHeight() ); ++y ) {
            for ( int x = 0; x < Math.min( a.getWidth(), b.getWidth() ); ++x ) {
                final int p = a.getRGB( x, y );
                final int q = b.getRGB( x, y );
                if ( ( Math.abs( ( ( p >> 16 ) & 255 ) - ( ( q >> 16 ) & 255 ) ) > by )
                        || ( Math.abs( ( ( p >> 8 ) & 255 ) - ( ( q >> 8 ) & 255 ) ) > by )
                        || ( Math.abs( ( p & 255 ) - ( q & 255 ) ) > by ) ) {
                    out.add( new int[] { x, y } );
                }
            }
        }
        return out;
    }

    /** Each pixel as {along the seam ray from the centre, across it}. */
    private static double[][] alongPerp( final TreePanel tp, final List<int[]> ink ) {
        final Point c = tp.circularCenterForTest();
        final double ux = Math.cos( tp.circularSeamAngleForTest() );
        final double uy = Math.sin( tp.circularSeamAngleForTest() );
        final double[][] out = new double[ ink.size() ][];
        for ( int i = 0; i < ink.size(); ++i ) {
            final double dx = ink.get( i )[ 0 ] - c.x;
            final double dy = ink.get( i )[ 1 ] - c.y;
            out[ i ] = new double[] { ( dx * ux ) + ( dy * uy ), ( dx * -uy ) + ( dy * ux ) };
        }
        return out;
    }

    private static boolean dark( final BufferedImage img, final int x, final int y ) {
        if ( ( x < 0 ) || ( y < 0 ) || ( x >= img.getWidth() ) || ( y >= img.getHeight() ) ) {
            return false;
        }
        final int p = img.getRGB( x, y );
        return ( ( p >> 16 ) & 255 ) < 0xA0;
    }

    private static double wrap( final double a ) {
        double d = a % ( 2 * Math.PI );
        if ( d > Math.PI ) {
            d -= 2 * Math.PI;
        }
        else if ( d <= -Math.PI ) {
            d += 2 * Math.PI;
        }
        return d;
    }

    private interface Body {

        void run( MainFrame frame, TreePanel tp, Options o );
    }

    private static void withDemo( final String demo, final Body body, final boolean[] ok ) {
        withDemo( demo, phy -> {
        }, body, ok );
    }

    private static void withDemo( final String demo, final java.util.function.Consumer<Phylogeny> tweak,
                                  final Body body, final boolean[] ok ) {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
            if ( !file.exists() ) {
                fail( ok, "demo tree missing: " + file.getAbsolutePath() );
                return;
            }
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                    .create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
            tweak.accept( phy );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "circular-scale" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    body.run( mf[ 0 ], mf[ 0 ].getMainPanel().getCurrentTreePanel(), mf[ 0 ].getOptions() );
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    fail( ok, demo + ": unexpected " + t );
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            fail( ok, demo + ": unexpected " + e );
        }
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [CircularScaleAxisTest] " + message );
        ok[ 0 ] = false;
    }

    private CircularScaleAxisTest() {
    }
}
