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
import java.nio.file.Files;

import javax.imageio.ImageIO;
import javax.imageio.ImageReader;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.metadata.IIOMetadataNode;
import javax.imageio.stream.ImageInputStream;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.w3c.dom.NodeList;

/**
 * End-to-end render test for "export at a fixed figure size" ({@link AptxUtil#writePhylogenyToGraphicsFileAtSize}
 * driven by an {@link ExportSizeSpec}). Drives a real {@link TreePanel} and checks that:
 * <ul>
 *   <li>a raster (PNG) export produces an image of EXACTLY the spec's pixel size and embeds the requested DPI;</li>
 *   <li>a vector (SVG) export produces a page of EXACTLY the spec's point size;</li>
 *   <li>the on-screen layout is RESTORED afterwards (a normal export before and after a fixed-size export renders
 *       identically -- the layout is not left at the export size); and</li>
 *   <li>a CIRCULAR (radial) tree exports at the target size, renders real ink, and restores its radial zoom.</li>
 * </ul>
 *
 * <p>Needs FlatLaf + a display, so {@link #test()} is a no-op (returns true) when headless -- the same pattern the
 * other GUI tests in this package use; it runs for real from {@link #main} or any non-headless invocation.
 */
public final class FixedExportSizeTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "FixedExportSize: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = Phylogeny.createInstanceFromNhxString( "(((A:0.2,B:0.2):0.2,(C:0.2,D:0.2):0.2):0.2,(E:0.2,F:0.2):0.2);" );
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "fixed size" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final JFrame f = (JFrame) mf[ 0 ];
                try {
                    f.setSize( 640, 420 );
                    f.validate();
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    final Options o = mp.getOptions();
                    mp.getControlPanel().showWhole();
                    if ( ( tp.getWidth() < 100 ) || ( tp.getHeight() < 100 ) ) {
                        return; // no usable viewport in this environment; nothing to assert
                    }
                    o.setTransparentExportBackground( false );

                    // 1. RASTER: 1200 x 900 px at 300 DPI. In PIXELS mode the pixel count round-trips exactly.
                    final ExportSizeSpec px = new ExportSizeSpec( ExportSizeSpec.Unit.PIXELS, 1200, 900, 300 );
                    final File png = File.createTempFile( "aptx_fixed_", ".png" );
                    try {
                        final String msg = AptxUtil.writePhylogenyToGraphicsFileAtSize( png.getAbsolutePath(), tp,
                                GraphicsExportType.PNG, o, px );
                        final BufferedImage img = ImageIO.read( png );
                        if ( ( img.getWidth() != px.imageWidthPx() ) || ( img.getHeight() != px.imageHeightPx() ) ) {
                            fail( ok, "fixed PNG must be exactly " + px.imageWidthPx() + "x" + px.imageHeightPx()
                                    + " px; got " + img.getWidth() + "x" + img.getHeight() );
                        }
                        if ( ( img.getWidth() != 1200 ) || ( img.getHeight() != 900 ) ) {
                            fail( ok, "a PIXELS-unit spec must round-trip to exactly the pixel count; got "
                                    + img.getWidth() + "x" + img.getHeight() );
                        }
                        final double dpi = readHorizontalDpi( png );
                        if ( Math.abs( dpi - 300.0 ) > 1.0 ) {
                            fail( ok, "fixed PNG must embed the requested ~300 DPI; got " + dpi );
                        }
                        if ( ( msg == null ) || !msg.contains( "1200 × 900 px" ) || !msg.contains( "300 DPI" )
                                || !msg.contains( "Fixed size: yes" ) || !msg.contains( "Font size:" )
                                || !msg.contains( "Physical:" ) || !msg.contains( " mm" ) ) {
                            fail( ok, "fixed PNG report should state px size + DPI + physical mm + fixed + font; got "
                                    + msg );
                        }
                    }
                    finally {
                        png.delete();
                    }

                    // 2. VECTOR: 4 x 3 inch at 72 DPI -> page EXACTLY 288 x 216 pt (the returned message states it).
                    final ExportSizeSpec in = new ExportSizeSpec( ExportSizeSpec.Unit.INCHES, 4, 3, 72 );
                    if ( ( in.layoutWidthPt() != 288 ) || ( in.layoutHeightPt() != 216 ) ) {
                        fail( ok, "precondition: 4x3in should be 288x216 pt; got " + in.layoutWidthPt() + "x"
                                + in.layoutHeightPt() );
                    }
                    final File svg = File.createTempFile( "aptx_fixed_", ".svg" );
                    try {
                        final String msg = AptxUtil.writePhylogenyToGraphicsFileAtSize( svg.getAbsolutePath(), tp,
                                GraphicsExportType.SVG, o, in );
                        if ( ( msg == null ) || !msg.contains( "288 × 216 pt" ) || !msg.contains( "Fixed size: yes" ) ) {
                            fail( ok, "fixed SVG page must be exactly 288x216 pt (fixed); report was " + msg );
                        }
                        final String body = new String( Files.readAllBytes( svg.toPath() ), "UTF-8" );
                        if ( !body.contains( "288" ) || !body.contains( "216" ) || !body.contains( "svg" ) ) {
                            fail( ok, "fixed SVG document should carry the 288x216 page size" );
                        }
                    }
                    finally {
                        svg.delete();
                    }

                    // 3. RESTORE: a normal (on-screen size) export before and after a fixed-size export must render
                    //    identically -- the fixed export must not leave the tree laid out at the export size.
                    final BufferedImage before = normalExport( tp, mp );
                    AptxUtil.writePhylogenyToGraphicsFileAtSize( File.createTempFile( "aptx_x_", ".png" ).getAbsolutePath(),
                            tp, GraphicsExportType.PNG,
                            o, new ExportSizeSpec( ExportSizeSpec.Unit.MILLIMETERS, 60, 45, 300 ) );
                    final BufferedImage after = normalExport( tp, mp );
                    if ( !nearlyIdentical( before, after ) ) {
                        fail( ok, "the on-screen layout was not restored after a fixed-size export (a normal export "
                                + "before/after differs)" );
                    }

                    // 3b. DYNA-HIDE report warning: with "Dyna Hide" on, a cramped export layout hides labels
                    //     (labelsDynamicallyHidden -> the report warns); a roomy one does not; off -> never.
                    final javax.swing.JCheckBox dh = mp.getControlPanel().getDynamicallyHideData();
                    if ( dh != null ) {
                        dh.setSelected( false );
                        int[] t = tp.layoutForExportSize( 300, 24 );
                        final boolean off_cramped = tp.labelsDynamicallyHidden();
                        tp.restoreLayoutAfterExport( t );
                        dh.setSelected( true );
                        t = tp.layoutForExportSize( 300, 24 );
                        final boolean on_cramped = tp.labelsDynamicallyHidden();
                        tp.restoreLayoutAfterExport( t );
                        t = tp.layoutForExportSize( 400, 1200 );
                        final boolean on_roomy = tp.labelsDynamicallyHidden();
                        tp.restoreLayoutAfterExport( t );
                        // a REAL cramped fixed export (Dyna Hide still on) must carry the warning line in the report
                        final File whpng = File.createTempFile( "aptx_dh_", ".png" );
                        try {
                            final String whmsg = AptxUtil.writePhylogenyToGraphicsFileAtSize( whpng.getAbsolutePath(),
                                    tp, GraphicsExportType.PNG, o,
                                    new ExportSizeSpec( ExportSizeSpec.Unit.PIXELS, 300, 24, 72 ) );
                            if ( ( whmsg == null ) || !whmsg.contains( "Warning:" ) || !whmsg.contains( "Auto-hide Labels" ) ) {
                                fail( ok, "a cramped fixed export with Dyna Hide on must warn in the report; got "
                                        + whmsg );
                            }
                        }
                        finally {
                            whpng.delete();
                        }
                        dh.setSelected( false );
                        if ( off_cramped || !on_cramped || on_roomy ) {
                            fail( ok, "Dyna-Hide warning wrong: off_cramped=" + off_cramped + " on_cramped="
                                    + on_cramped + " on_roomy=" + on_roomy );
                        }
                    }

                    // 4. RADIAL (circular): the fixed export produces the target size, real ink, and restores the zoom.
                    tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    mp.getControlPanel().showWhole();
                    final int radial_before = tp.radialDiameter();
                    final ExportSizeSpec sq = new ExportSizeSpec( ExportSizeSpec.Unit.PIXELS, 700, 700, 150 );
                    final File cpng = File.createTempFile( "aptx_circ_", ".png" );
                    try {
                        AptxUtil.writePhylogenyToGraphicsFileAtSize( cpng.getAbsolutePath(), tp, GraphicsExportType.PNG,
                                o, sq );
                        final BufferedImage cimg = ImageIO.read( cpng );
                        if ( ( cimg.getWidth() != sq.imageWidthPx() ) || ( cimg.getHeight() != sq.imageHeightPx() ) ) {
                            fail( ok, "fixed circular PNG must be exactly " + sq.imageWidthPx() + "x"
                                    + sq.imageHeightPx() + "; got " + cimg.getWidth() + "x" + cimg.getHeight() );
                        }
                        if ( nonBackgroundPixels( cimg ) < 200 ) {
                            fail( ok, "fixed circular export drew (almost) no ink -- the radial tree did not render" );
                        }
                        // CENTERING: the ring must be centred in the EXPORT canvas (not the on-screen panel). The bug
                        // centred at getWidth()/2 (screen size != export size for a fixed export) -> the ring is
                        // pushed off-canvas / clipped. Assert the ink's bounding-box centre is near the image centre.
                        final int[] b = inkBounds( cimg );
                        if ( b != null ) {
                            final double cx = ( b[ 0 ] + b[ 2 ] ) / 2.0;
                            final double cy = ( b[ 1 ] + b[ 3 ] ) / 2.0;
                            if ( ( Math.abs( cx - ( cimg.getWidth() / 2.0 ) ) > ( cimg.getWidth() * 0.15 ) )
                                    || ( Math.abs( cy - ( cimg.getHeight() / 2.0 ) ) > ( cimg.getHeight() * 0.15 ) ) ) {
                                fail( ok, "fixed circular export is not centred in the canvas: ink centre (" + cx + ","
                                        + cy + ") vs image centre (" + ( cimg.getWidth() / 2 ) + ","
                                        + ( cimg.getHeight() / 2 ) + ")" );
                            }
                        }
                    }
                    finally {
                        cpng.delete();
                    }
                    if ( tp.radialDiameter() != radial_before ) {
                        fail( ok, "the radial zoom must be restored after a fixed-size export; was " + radial_before
                                + ", now " + tp.radialDiameter() );
                    }

                    // 5. PDF fixed-size (exact page): reach the PdfExporter.writePhylogenyToPdfExactSize path the way
                    //    MainFrame does (layout -> write -> restore), and check the report states the exact page size.
                    tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    mp.getControlPanel().showWhole();
                    final ExportSizeSpec pspec = new ExportSizeSpec( ExportSizeSpec.Unit.INCHES, 4, 3, 72 );
                    final File pdf = File.createTempFile( "aptx_fixed_", ".pdf" );
                    final int[] ptok = tp.layoutForExportSize( pspec.layoutWidthPt(), pspec.layoutHeightPt() );
                    try {
                        final String pmsg = PdfExporter.writePhylogenyToPdfExactSize( pdf.getAbsolutePath(), tp,
                                pspec.layoutWidthPt(), pspec.layoutHeightPt(), o.isGraphicsExportWhiteBackground() );
                        if ( ( pmsg == null ) || !pmsg.contains( "288 × 216 pt" ) || !pmsg.contains( "Fixed size: yes" )
                                || ( pdf.length() < 400 ) ) {
                            fail( ok, "fixed PDF page must be exactly 288x216 pt (fixed) + a real file; report was "
                                    + pmsg + ", len " + pdf.length() );
                        }
                    }
                    finally {
                        tp.restoreLayoutAfterExport( ptok );
                        pdf.delete();
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
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

    /** A plain, on-screen-size raster export at scale 1 using the CURRENT panel layout (no re-layout), so a broken
     *  restore in the fixed-size path shows up as a different render here. */
    private static BufferedImage normalExport( final TreePanel tp, final MainPanel mp ) throws Exception {
        mp.getOptions().setRasterExportScale( 1 );
        final File out = File.createTempFile( "aptx_normal_", ".png" );
        try {
            AptxUtil.writePhylogenyToGraphicsFile( out.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                   mp.getControlPanel(), GraphicsExportType.PNG, mp.getOptions() );
            return ImageIO.read( out );
        }
        finally {
            out.delete();
        }
    }

    /** Two renders of the same layout are byte-identical in principle; allow a tiny number of differing pixels for
     *  robustness. A broken restore mislays the whole tree, so it differs in far more than the threshold. */
    private static boolean nearlyIdentical( final BufferedImage a, final BufferedImage b ) {
        if ( ( a.getWidth() != b.getWidth() ) || ( a.getHeight() != b.getHeight() ) ) {
            return false;
        }
        int diff = 0;
        for ( int y = 0; y < a.getHeight(); y += 4 ) {
            for ( int x = 0; x < a.getWidth(); x += 4 ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    ++diff;
                }
            }
        }
        final int sampled = ( ( a.getWidth() + 3 ) / 4 ) * ( ( a.getHeight() + 3 ) / 4 );
        return diff <= ( sampled / 200 ); // <= 0.5% of sampled pixels may differ
    }

    /** Bounding box [minX, minY, maxX, maxY] of the non-(white)-background pixels, or null if the image is blank. */
    private static int[] inkBounds( final BufferedImage img ) {
        int minx = Integer.MAX_VALUE, miny = Integer.MAX_VALUE, maxx = -1, maxy = -1;
        for ( int y = 0; y < img.getHeight(); y += 2 ) {
            for ( int x = 0; x < img.getWidth(); x += 2 ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) != 0xFFFFFF ) {
                    if ( x < minx ) { minx = x; }
                    if ( x > maxx ) { maxx = x; }
                    if ( y < miny ) { miny = y; }
                    if ( y > maxy ) { maxy = y; }
                }
            }
        }
        return ( maxx < 0 ) ? null : new int[] { minx, miny, maxx, maxy };
    }

    /** Count of pixels that are not the (white) export background, on an opaque light-theme export. */
    private static int nonBackgroundPixels( final BufferedImage img ) {
        int n = 0;
        for ( int y = 0; y < img.getHeight(); y += 3 ) {
            for ( int x = 0; x < img.getWidth(); x += 3 ) {
                final int rgb = img.getRGB( x, y ) & 0xFFFFFF;
                if ( rgb != 0xFFFFFF ) {
                    ++n;
                }
            }
        }
        return n;
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
            final double mm_per_pixel = Double
                    .parseDouble( ( (IIOMetadataNode) sizes.item( 0 ) ).getAttribute( "value" ) );
            return 25.4 / mm_per_pixel;
        }
        finally {
            reader.dispose();
        }
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [FixedExportSizeTest] " + msg );
        ok[ 0 ] = false;
    }

    private FixedExportSizeTest() {
    }
}
