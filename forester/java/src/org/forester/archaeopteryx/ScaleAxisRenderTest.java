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
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the scale-axis demo phylogram (forester/demo/scale-axis.xml) with "Scale Axis" ON and asserts a wide
 * horizontal axis line appears along the bottom -- and NONE when it is off. Headful; a green no-op when headless.
 * Dogfoods the demo tree (the "feature test loads its demo tree" convention).
 */
public final class ScaleAxisRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ScaleAxisRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return axisRendersOk();
    }

    private static boolean axisRendersOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/scale-axis.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "axis" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( false ); // keep our forced colors (no light-theme switch)
                    final int magenta = 0xFF00FF; // the scale axis is drawn in the BRANCH_LENGTH color -- force it unique
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, new Color( 255, 255, 255 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BRANCH_LENGTH, new Color( magenta ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BRANCH, new Color( 0, 0, 0 ) ); // branches black
                    tp.getTreeColorSet().setColorSchema( 0 );
                    final int w = 800, h = 480;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    o.setShowScaleAxis( false );
                    final int off = widestBottomRun(
                            AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ), magenta );
                    o.setShowScaleAxis( true );
                    final BufferedImage on_img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int on = widestBottomRun( on_img, magenta );
                    if ( on < ( w / 3 ) ) {
                        fail( ok, "the scale axis line should span much of the width, got " + on + "px" );
                    }
                    if ( off >= ( w / 4 ) ) {
                        fail( ok, "no wide axis line should appear when Scale Axis is off, got " + off + "px" );
                    }
                    // the tree name (default ON, ~120px "Scale axis (demo)") must be RAISED above the axis line: below
                    // the axis line, the left region holds only narrow tick labels (~<20px), never the wide name.
                    final int axis_row = widestBottomRunRow( on_img, magenta );
                    final int below = widestRunInRegion( on_img, magenta, axis_row + 2, h, 0, ( w * 2 ) / 5 );
                    if ( below >= 40 ) {
                        fail( ok, "the tree name must be raised above the axis, not left below it (run " + below + "px)" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Widest horizontal run (px) of {@code color} in the bottom ~28% of the image, where the axis line sits. */
    private static int widestBottomRun( final BufferedImage img, final int color ) {
        final int y0 = (int) ( img.getHeight() * 0.72 );
        int best = 0;
        for( int y = y0; y < img.getHeight(); ++y ) {
            int run = 0;
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) == color ) {
                    if ( ++run > best ) {
                        best = run;
                    }
                }
                else {
                    run = 0;
                }
            }
        }
        return best;
    }

    /** The y of the row (bottom ~28%) with the widest {@code color} run -- the axis line. */
    private static int widestBottomRunRow( final BufferedImage img, final int color ) {
        final int y0 = (int) ( img.getHeight() * 0.72 );
        int best = -1, best_y = img.getHeight() - 1;
        for( int y = y0; y < img.getHeight(); ++y ) {
            int run = 0, row_best = 0;
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) == color ) {
                    if ( ++run > row_best ) {
                        row_best = run;
                    }
                }
                else {
                    run = 0;
                }
            }
            if ( row_best > best ) {
                best = row_best;
                best_y = y;
            }
        }
        return best_y;
    }

    /** Widest horizontal run (px) of {@code color} within the rectangle [x0,x1) x [y0,y1). */
    private static int widestRunInRegion( final BufferedImage img, final int color, final int y0, final int y1,
                                          final int x0, final int x1 ) {
        int best = 0;
        for( int y = Math.max( 0, y0 ); y < Math.min( img.getHeight(), y1 ); ++y ) {
            int run = 0;
            for( int x = Math.max( 0, x0 ); x < Math.min( img.getWidth(), x1 ); ++x ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) == color ) {
                    if ( ++run > best ) {
                        best = run;
                    }
                }
                else {
                    run = 0;
                }
            }
        }
        return best;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ScaleAxisRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ScaleAxisRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private ScaleAxisRenderTest() {
    }
}
