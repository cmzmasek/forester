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

import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.phylogeny.Phylogeny;

/**
 * End-to-end render test for high-resolution raster export ({@link Options#getRasterExportScale()}) and
 * transparent-background PNG export ({@link Options#isTransparentExportBackground()}). Drives a real
 * {@link TreePanel} through the PNG/JPG export path and checks that the scale multiplies the image
 * dimensions (a true re-render), and that a transparent PNG has an unfilled (alpha 0) corner while the
 * same setting on JPG -- which cannot carry alpha -- stays opaque.
 *
 * <p>Needs FlatLaf + a display, so {@link #test()} is a no-op (returns true) when headless -- the same
 * pattern the other GUI tests in this package use; it runs for real from {@link #main} or any
 * non-headless invocation.
 */
public final class ExportResolutionTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ExportResolution: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = Phylogeny.createInstanceFromNhxString( "(((A,B),(C,D)),(E,F));" );
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "export test" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final JFrame f = (JFrame) mf[ 0 ];
                try {
                    f.setSize( 600, 360 );
                    f.validate();
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    mp.getControlPanel().showWhole();
                    final int w = tp.getWidth();
                    final int h = tp.getHeight();
                    if ( ( w < 100 ) || ( h < 100 ) ) {
                        return; // no usable viewport in this environment; nothing to assert
                    }
                    // 1x: image == panel size, opaque corner
                    mp.getOptions().setRasterExportScale( 1 );
                    mp.getOptions().setTransparentExportBackground( false );
                    final BufferedImage s1 = export( tp, mp, GraphicsExportType.PNG );
                    if ( ( s1.getWidth() != w ) || ( s1.getHeight() != h ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  1x export should equal panel size; got " + s1.getWidth() + "x" + s1.getHeight() );
                    }
                    if ( alpha( s1, 2, 2 ) != 255 ) {
                        ok[ 0 ] = false;
                        System.out.println( "  opaque export should have an opaque corner" );
                    }
                    // 3x: image dimensions tripled (true re-render onto a larger canvas)
                    mp.getOptions().setRasterExportScale( 3 );
                    final BufferedImage s3 = export( tp, mp, GraphicsExportType.PNG );
                    if ( ( s3.getWidth() != ( 3 * w ) ) || ( s3.getHeight() != ( 3 * h ) ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  3x export should triple the dimensions; got " + s3.getWidth() + "x"
                                + s3.getHeight() + " vs " + ( 3 * w ) + "x" + ( 3 * h ) );
                    }
                    // transparent PNG: unfilled (alpha 0) corner
                    mp.getOptions().setRasterExportScale( 1 );
                    mp.getOptions().setTransparentExportBackground( true );
                    final BufferedImage tr = export( tp, mp, GraphicsExportType.PNG );
                    if ( alpha( tr, 2, 2 ) != 0 ) {
                        ok[ 0 ] = false;
                        System.out.println( "  transparent PNG should have an alpha-0 corner; got alpha " + alpha( tr, 2, 2 ) );
                    }
                    // same setting on JPG (no alpha channel) must stay opaque
                    final BufferedImage trj = export( tp, mp, GraphicsExportType.JPG );
                    if ( alpha( trj, 2, 2 ) != 255 ) {
                        ok[ 0 ] = false;
                        System.out.println( "  JPG cannot be transparent; corner should be opaque" );
                    }
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                    ok[ 0 ] = false;
                }
                finally {
                    f.dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    private static BufferedImage export( final TreePanel tp, final MainPanel mp, final GraphicsExportType type )
            throws Exception {
        final File out = File.createTempFile( "aptx_export_" + type, "." + type.toString() );
        try {
            AptxUtil.writePhylogenyToGraphicsFile( out.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                   mp.getControlPanel(), type, mp.getOptions() );
            return ImageIO.read( out );
        }
        finally {
            out.delete();
        }
    }

    private static int alpha( final BufferedImage img, final int x, final int y ) {
        return ( img.getRGB( x, y ) >>> 24 ) & 0xFF;
    }

    private ExportResolutionTest() {
    }
}
