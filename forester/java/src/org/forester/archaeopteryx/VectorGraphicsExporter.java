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

import java.awt.Graphics2D;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.function.Consumer;

import org.forester.phylogeny.Phylogeny;

import de.erichseifert.vectorgraphics2d.Processor;
import de.erichseifert.vectorgraphics2d.Processors;
import de.erichseifert.vectorgraphics2d.VectorGraphics2D;
import de.erichseifert.vectorgraphics2d.util.PageSize;

/**
 * Vector-graphics (SVG / EPS) export for Archaeopteryx, built on the VectorGraphics2D library.
 *
 * <p>Like {@link PdfExporter}, this reuses the <em>exact</em> on-screen paint routine
 * ({@link TreePanel#paintPhylogeny}): a {@link VectorGraphics2D} (which is itself a
 * {@link Graphics2D} that records every Java2D call) is handed to {@code paintPhylogeny}, so the
 * vector file is produced from the same call stream that paints the screen. The resulting figure is
 * true vector output that scales losslessly and can be finished in Illustrator / Inkscape -- which
 * raster PNG/TIFF cannot. WYSIWYG by default, except when {@code white_background} requests a
 * document-ready (light theme, white page) render -- see {@link ExportTheme}.
 */
final class VectorGraphicsExporter {

    /** A vector output format understood by VectorGraphics2D's processors. */
    enum Format {
        SVG( "svg" ),
        EPS( "eps" );

        private final String _id;

        Format( final String id ) {
            _id = id;
        }

        String id() {
            return _id;
        }
    }

    private static final int    HEIGHT_LIMIT = 100;
    private static final int    WIDTH_LIMIT  = 60;
    private static final int    MARGIN_X     = 20;
    private static final int    MARGIN_Y     = 10;
    // VectorGraphics2D sizes a page in millimetres and treats the DRAWING's coordinates as millimetres too,
    // scaling them to PostScript points (1 mm = 72/25.4 pt) for EPS output. SVG uses the page numbers directly as
    // the pixel viewBox (drawing maps 1:1), so it needs no conversion. For EPS we therefore both size the page in
    // px->mm AND pre-scale the drawing px->mm (see render): that keeps the point bounding box numerically equal to
    // the pixel size while the content still fills it. (Sizing the page alone -- the old approach -- left the
    // pixel-valued drawing ~2.83x too big for the page, so only a corner of the figure appeared.)
    private static final double PX_TO_MM     = 25.4 / 72.0;

    private VectorGraphicsExporter() {
        // not instantiable
    }

    /**
     * Renders a paint callback into vector-graphics bytes of the requested format. Pure and
     * GUI-free (needs no {@link TreePanel}), hence unit-testable headlessly. The painter receives a
     * {@link VectorGraphics2D}; whatever it draws is what the file contains.
     *
     * @param width  intended width in pixels (becomes the SVG viewBox width / EPS point bbox width)
     * @param height intended height in pixels
     * @return the encoded SVG or EPS document
     */
    static byte[] render( final int width, final int height, final Format fmt, final boolean outline_text,
                          final Consumer<Graphics2D> painter )
            throws IOException {
        final ByteArrayOutputStream out = new ByteArrayOutputStream( 1 << 16 );
        render( width, height, fmt, outline_text, painter, out );
        return out.toByteArray();
    }

    /** Core renderer: streams the encoded document straight to {@code out}, no intermediate copy. */
    private static void render( final int width,
                                final int height,
                                final Format fmt,
                                final boolean outline_text,
                                final Consumer<Graphics2D> painter,
                                final OutputStream out )
            throws IOException {
        // When outline_text, every text string is rendered as glyph outlines: VectorGraphics2D embeds no
        // fonts, so font-referenced text would be substituted by the viewer (EPS -> Times serif, SVG ->
        // generic sans, italics lost). See OutliningVectorGraphics2D. Off keeps selectable text.
        final VectorGraphics2D g = outline_text ? new OutliningVectorGraphics2D() : new VectorGraphics2D();
        final PageSize page;
        if ( fmt == Format.SVG ) {
            page = new PageSize( width, height );
        }
        else {
            // EPS: page in px->mm, AND draw in mm (pre-scale px->mm) so the mm-interpreted commands fill the mm page
            page = new PageSize( width * PX_TO_MM, height * PX_TO_MM );
            g.scale( PX_TO_MM, PX_TO_MM );
        }
        painter.accept( g );
        final Processor processor = Processors.get( fmt.id() );
        processor.getDocument( g.getCommands(), page ).writeTo( out );
    }

    /**
     * Exports the tree currently shown in {@code tree_panel} to a vector-graphics file. Mirrors
     * {@link PdfExporter#writePhylogenyToPdf}: same size clamping/margins and the same
     * {@code paintPhylogeny(g, true, false, ...)} call.
     *
     * @return a human-readable description of what was written, or "" if the tree was empty
     */
    static String writePhylogenyToVectorGraphicsFile( final String file_name,
                                                      final TreePanel tree_panel,
                                                      final int width,
                                                      final int height,
                                                      final Format fmt,
                                                      final boolean outline_text,
                                                      final boolean white_background )
            throws IOException {
        final Phylogeny phylogeny = tree_panel.getPhylogeny();
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return "";
        }
        if ( tree_panel.getMainPanel().getTreeFontSet().getSmallFont().getSize() < 1 ) {
            throw new IOException( "fonts are too small for vector graphics export" );
        }
        if ( tree_panel.getMainPanel().getTreeFontSet().getLargeFont().getSize() < 1 ) {
            throw new IOException( "fonts are too small for vector graphics export" );
        }
        final File file = new File( file_name );
        if ( file.isDirectory() ) {
            throw new IOException( "[" + file_name + "] is a directory" );
        }
        final int my_width = ( width < WIDTH_LIMIT ? WIDTH_LIMIT : width ) + ( 2 * MARGIN_X );
        final int my_height = ( height < HEIGHT_LIMIT ? HEIGHT_LIMIT : height ) + ( 2 * MARGIN_Y );
        // Document-ready (light theme) export when requested -- same behavior as raster/clipboard, so a
        // dark-theme vector figure isn't light-on-(page-)white; restored afterwards. The paint exception (if any)
        // is captured and reported AFTER restore(), so the modal error dialog never pumps the EDT while the
        // transient light theme is applied (which would flash the on-screen tree white).
        final ExportTheme export_theme = ExportTheme.applyIf( tree_panel, white_background );
        final Exception[] paint_error = { null };
        try ( final OutputStream out = new BufferedOutputStream( new FileOutputStream( file ) ) ) {
            render( my_width, my_height, fmt, outline_text, g -> {
                if ( white_background ) {
                    // an explicit white page (SVG/EPS are otherwise transparent), so the figure reads on a
                    // white document -- matching the raster/PDF background
                    g.setColor( java.awt.Color.WHITE );
                    g.fillRect( 0, 0, my_width, my_height );
                }
                try {
                    tree_panel.paintPhylogeny( g, true, false, my_width, my_height, 0, 0 );
                }
                catch ( final Exception e ) {
                    paint_error[ 0 ] = e;
                }
            }, out );
        }
        finally {
            export_theme.restore();
        }
        if ( paint_error[ 0 ] != null ) {
            AptxUtil.unexpectedException( paint_error[ 0 ] );
        }
        return file.toString() + " [size: " + my_width + ", " + my_height + "]";
    }
}
