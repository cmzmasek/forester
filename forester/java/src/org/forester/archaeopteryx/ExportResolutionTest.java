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
import javax.imageio.ImageReader;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.metadata.IIOMetadataNode;
import javax.imageio.stream.ImageInputStream;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.phylogeny.Phylogeny;
import org.w3c.dom.NodeList;

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
        if ( !dpiRoundTripOk() || !truncateOnOverwriteOk() ) {
            return false;
        }
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
                    // the 3x PNG must carry a physical resolution of 3*72 = 216 DPI (so it drops into a document
                    // at its intended print size instead of a nominal 72); the 1x PNG must carry 72
                    final double dpi3 = exportedDpi( tp, mp, 3 );
                    if ( Math.abs( dpi3 - 216.0 ) > 1.0 ) {
                        ok[ 0 ] = false;
                        System.out.println( "  3x PNG should embed ~216 DPI; got " + dpi3 );
                    }
                    final double dpi1 = exportedDpi( tp, mp, 1 );
                    if ( Math.abs( dpi1 - 72.0 ) > 1.0 ) {
                        ok[ 0 ] = false;
                        System.out.println( "  1x PNG should embed ~72 DPI; got " + dpi1 );
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

    /** Exports a PNG at {@code scale} through the real file path and returns the horizontal DPI it embedded. */
    private static double exportedDpi( final TreePanel tp, final MainPanel mp, final int scale ) throws Exception {
        mp.getOptions().setRasterExportScale( scale );
        mp.getOptions().setTransparentExportBackground( false );
        final File out = File.createTempFile( "aptx_dpi_" + scale, ".png" );
        try {
            AptxUtil.writePhylogenyToGraphicsFile( out.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                   mp.getControlPanel(), GraphicsExportType.PNG, mp.getOptions() );
            return readHorizontalDpi( out );
        }
        finally {
            out.delete();
        }
    }

    /** Pure, headless-safe check of {@link AptxUtil#writePngWithDpi}: a written pHYs resolution reads back. */
    private static boolean dpiRoundTripOk() {
        try {
            final BufferedImage img = new BufferedImage( 12, 8, BufferedImage.TYPE_INT_RGB );
            final File out = File.createTempFile( "aptx_dpi_pure", ".png" );
            try {
                AptxUtil.writePngWithDpi( img, out, 288.0 );
                final BufferedImage back = ImageIO.read( out );
                if ( ( back == null ) || ( back.getWidth() != 12 ) || ( back.getHeight() != 8 ) ) {
                    System.out.println( "  [ExportResolutionTest] writePngWithDpi produced an unreadable/wrong-size PNG" );
                    return false;
                }
                final double dpi = readHorizontalDpi( out );
                if ( Math.abs( dpi - 288.0 ) > 1.0 ) {
                    System.out.println( "  [ExportResolutionTest] writePngWithDpi should embed ~288 DPI; got " + dpi );
                    return false;
                }
                return true;
            }
            finally {
                out.delete();
            }
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    /** Overwriting a larger PNG with a smaller one must TRUNCATE the file: writePngWithDpi writes through a
     *  truncating stream, so no stale trailing bytes from the previous (larger) export survive. Headless-safe. */
    private static boolean truncateOnOverwriteOk() {
        try {
            final File out = File.createTempFile( "aptx_trunc", ".png" );
            try {
                AptxUtil.writePngWithDpi( new BufferedImage( 400, 400, BufferedImage.TYPE_INT_RGB ), out, 72.0 );
                final long large_len = out.length();
                AptxUtil.writePngWithDpi( new BufferedImage( 6, 6, BufferedImage.TYPE_INT_RGB ), out, 72.0 );
                final long small_len = out.length();
                if ( small_len >= large_len ) {
                    System.out.println( "  [ExportResolutionTest] overwrite did not truncate: large=" + large_len
                            + " small=" + small_len + " (stale trailing bytes -> corrupt PNG)" );
                    return false;
                }
                final BufferedImage back = ImageIO.read( out );
                if ( ( back == null ) || ( back.getWidth() != 6 ) || ( back.getHeight() != 6 ) ) {
                    System.out.println( "  [ExportResolutionTest] overwritten PNG is not a clean 6x6 image" );
                    return false;
                }
                return true;
            }
            finally {
                out.delete();
            }
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    /** The horizontal DPI recorded in a PNG's standard metadata (pHYs), or -1 if absent. */
    private static double readHorizontalDpi( final File png ) throws Exception {
        final ImageReader reader = ImageIO.getImageReadersByFormatName( "png" ).next();
        try ( final ImageInputStream iis = ImageIO.createImageInputStream( png ) ) {
            reader.setInput( iis );
            final IIOMetadata meta = reader.getImageMetadata( 0 );
            final IIOMetadataNode root = (IIOMetadataNode) meta.getAsTree( "javax_imageio_1.0" );
            final NodeList sizes = root.getElementsByTagName( "HorizontalPixelSize" );
            if ( sizes.getLength() == 0 ) {
                return -1;
            }
            final double mm_per_pixel = Double.parseDouble( ( (IIOMetadataNode) sizes.item( 0 ) ).getAttribute( "value" ) );
            return 25.4 / mm_per_pixel;
        }
        finally {
            reader.dispose();
        }
    }

    private ExportResolutionTest() {
    }
}
