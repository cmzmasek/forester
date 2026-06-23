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
 * vector file is produced from the same call stream that paints the screen == WYSIWYG. The
 * resulting figure is true vector output that scales losslessly and can be finished in
 * Illustrator / Inkscape -- which raster PNG/TIFF cannot.
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
    // VectorGraphics2D sizes a page in millimeters. SVG output uses the page numbers directly as the
    // pixel viewBox, but the EPS/PDF processors convert mm -> PostScript points (1 mm = 72/25.4 pt).
    // Converting pixels -> mm for those formats makes the point bounding box numerically equal to the
    // pixel size, so 1 drawing-pixel maps to 1 point and the page tightly bounds the content.
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
        painter.accept( g );
        final Processor processor = Processors.get( fmt.id() );
        final PageSize page = ( fmt == Format.SVG ) ? new PageSize( width, height )
                : new PageSize( width * PX_TO_MM, height * PX_TO_MM );
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
                                                      final boolean outline_text )
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
        try ( final OutputStream out = new BufferedOutputStream( new FileOutputStream( file ) ) ) {
            render( my_width, my_height, fmt, outline_text, g -> {
                try {
                    tree_panel.paintPhylogeny( g, true, false, my_width, my_height, 0, 0 );
                }
                catch ( final Exception e ) {
                    AptxUtil.unexpectedException( e );
                }
            }, out );
        }
        return file.toString() + " [size: " + my_width + ", " + my_height + "]";
    }
}
