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
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

/**
 * Writes a {@link TanglegramPanel} figure to a file. Vector formats (PDF / SVG / EPS) are the publication-quality
 * choices; PNG is a raster fallback. Every format renders through {@link TanglegramPanel#paintForExport} (white
 * background, black ink, WYSIWYG for the connector colouring / labels / legend), reusing the app's existing
 * exporters -- {@link VectorGraphicsExporter} for SVG/EPS, {@link PdfExporter} for PDF -- so the figures match the
 * tree exports.
 */
final class TanglegramExporter {

    /** The supported output formats (a label for the file-chooser + the file extension). */
    enum Format {
        PDF( "pdf", "PDF (vector)" ),
        SVG( "svg", "SVG (vector)" ),
        EPS( "eps", "EPS (vector)" ),
        PNG( "png", "PNG (image)" );

        private final String _extension;
        private final String _label;

        Format( final String extension, final String label ) {
            _extension = extension;
            _label = label;
        }

        String extension() {
            return _extension;
        }

        String label() {
            return _label;
        }
    }

    private static final int  PNG_SCALE         = 2;             // render the raster at 2x for a crisper image
    private static final long MAX_RASTER_PIXELS = 100_000_000L;  // guard against an OutOfMemoryError on a huge tree

    /** Writes the tanglegram in {@code panel} to {@code file} in {@code format} at (w x h). Returns a short message. */
    static String write( final File file, final Format format, final TanglegramPanel panel, final int w, final int h )
            throws IOException {
        switch ( format ) {
            case SVG:
            case EPS: {
                final VectorGraphicsExporter.Format vector = ( format == Format.SVG ) ? VectorGraphicsExporter.Format.SVG
                        : VectorGraphicsExporter.Format.EPS;
                Files.write( file.toPath(),
                             VectorGraphicsExporter.render( w, h, vector, true, g -> panel.paintForExport( g, w, h ) ) );
                return file.toString() + " [" + w + " x " + h + "]";
            }
            case PDF:
                return PdfExporter.writeToPdf( file.getAbsolutePath(), w, h, g -> panel.paintForExport( g, w, h ) );
            case PNG:
                return writePng( file, panel, w, h );
            default:
                throw new IOException( "unsupported tanglegram export format: " + format );
        }
    }

    private static String writePng( final File file, final TanglegramPanel panel, final int w, final int h )
            throws IOException {
        // drop to 1x if 2x would exceed the raster cap; refuse (cleanly) if even 1x is too big for a raster
        int scale = PNG_SCALE;
        while ( ( scale > 1 ) && ( ( (long) w * scale ) * ( (long) h * scale ) > MAX_RASTER_PIXELS ) ) {
            --scale;
        }
        if ( ( ( (long) w * scale ) * ( (long) h * scale ) ) > MAX_RASTER_PIXELS ) {
            throw new IOException( "the tanglegram is too large to export as PNG (" + w + " x " + h
                    + " px) -- use PDF, SVG, or EPS instead" );
        }
        final BufferedImage img = new BufferedImage( w * scale, h * scale, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            g.scale( scale, scale );
            panel.paintForExport( g, w, h );
        }
        finally {
            g.dispose();
        }
        AptxUtil.writePngWithDpi( img, file, 72.0 * scale );
        return file.toString() + " [" + ( w * scale ) + " x " + ( h * scale ) + "]";
    }

    private TanglegramExporter() {
    }
}
