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
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.function.Consumer;

import org.forester.phylogeny.Phylogeny;

import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.Rectangle;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfWriter;

/*
 * PDF export via OpenPDF (https://github.com/LibrePDF/OpenPDF), the LGPL/MPL community fork of
 * iText 4 -- API-compatible (com.lowagie.text.*) but without iText 5's AGPL license.
 */
final class PdfExporter {

    private static final int HEIGHT_LIMIT = 100;
    private static final int WIDTH_LIMIT  = 60;
    private static final int MARGIN_X = 20;
    private static final int MARGIN_Y = 10;
    
    private PdfExporter() {
        // Empty constructor.
    }

    /**
     * Generic vector PDF export for any figure (not just a TreePanel): hands the {@code painter} a shape-based
     * (outline-text, so no font embedding) {@link Graphics2D} sized (w x h) to draw into. The painter is responsible
     * for its own background fill. Text is drawn as vector outlines via {@code createGraphicsShapes}, exactly as the
     * TreePanel PDF path, so the figure is fully portable.
     */
    static String writeToPdf( final String file_name, final int width, final int height,
                              final Consumer<Graphics2D> painter )
            throws IOException {
        final int my_width = Math.max( width, WIDTH_LIMIT );
        final int my_height = Math.max( height, HEIGHT_LIMIT );
        final File file = new File( file_name );
        if ( file.isDirectory() ) {
            throw new IOException( "[" + file_name + "] is a directory" );
        }
        final Document document = new Document();
        document.setPageSize( new Rectangle( my_width, my_height ) );
        document.setMargins( 0, 0, 0, 0 ); // the figure supplies its own margins
        PdfWriter writer;
        try {
            writer = PdfWriter.getInstance( document, new FileOutputStream( file_name ) );
        }
        catch ( final DocumentException e ) {
            throw new IOException( e );
        }
        document.open();
        final PdfContentByte cb = writer.getDirectContent();
        final Graphics2D g2 = cb.createGraphicsShapes( my_width, my_height );
        Exception paint_error = null;
        try {
            painter.accept( g2 );
        }
        catch ( final Exception e ) {
            paint_error = e;
        }
        finally {
            try {
                g2.dispose();
                document.close();
            }
            catch ( final Exception e ) {
                // Do nothing.
            }
        }
        // Unlike the TreePanel path (which reports via a dialog), propagate to the caller so it reports ONE clean
        // failure -- and don't leave a truncated PDF behind.
        if ( paint_error != null ) {
            file.delete();
            throw new IOException( "PDF rendering failed", paint_error );
        }
        return file.toString() + " [size: " + my_width + ", " + my_height + "]";
    }

    /**
     * WYSIWYG tree PDF export: the tree is laid out at {@code width} x {@code height} (its on-screen size), and the
     * page is that plus a small margin. This is the default File -> Export to PDF path.
     */
    static String writePhylogenyToPdf( final String file_name, final TreePanel tree_panel, final int width,
                                       final int height, final boolean white_background )
            throws IOException {
        final int page_h = ( ( height < HEIGHT_LIMIT ) ? HEIGHT_LIMIT : height ) + ( 2 * MARGIN_Y );
        final int page_w = ( ( width < WIDTH_LIMIT ) ? WIDTH_LIMIT : width ) + ( 2 * MARGIN_X );
        return renderPhylogenyToPdf( file_name, tree_panel, page_w, page_h, white_background, false );
    }

    /**
     * Fixed-size ("export at exactly this size") PDF export: the PDF page is EXACTLY {@code page_w} x {@code page_h}
     * points (no added margin), and the tree was already laid out to fill that frame by the caller (see
     * {@link TreePanel#layoutForExportSize}). A small floor keeps a degenerate size renderable.
     */
    static String writePhylogenyToPdfExactSize( final String file_name, final TreePanel tree_panel,
                                                final int page_w, final int page_h, final boolean white_background )
            throws IOException {
        return renderPhylogenyToPdf( file_name, tree_panel, Math.max( page_w, WIDTH_LIMIT ),
                                     Math.max( page_h, HEIGHT_LIMIT ), white_background, true );
    }

    /** Shared core: renders the tree onto a PDF page of EXACTLY {@code page_w} x {@code page_h} points. Text is drawn
     *  as vector outlines (createGraphicsShapes), so the figure needs no font embedding/mapping and is fully
     *  portable; it also sidesteps the bold-glyph stroke-color bleed that the glyph-font path has in iText/OpenPDF. */
    private static String renderPhylogenyToPdf( final String file_name, final TreePanel tree_panel,
                                                final int page_w, final int page_h, final boolean white_background,
                                                final boolean fixed_size )
            throws IOException {
        final Phylogeny phylogeny = tree_panel.getPhylogeny();
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return "";
        }
        if ( tree_panel.getMainPanel().getTreeFontSet().getSmallFont().getSize() < 1 ) {
            throw new IOException( "fonts are too small for PDF export" );
        }
        if ( tree_panel.getMainPanel().getTreeFontSet().getLargeFont().getSize() < 1 ) {
            throw new IOException( "fonts are too small for PDF export" );
        }
        final File file = new File( file_name );
        if ( file.isDirectory() ) {
            throw new IOException( "[" + file_name + "] is a directory" );
        }
        final Document document = new Document();
        document.setPageSize( new Rectangle( page_w, page_h ) );
        // the figure lays out into the full page (with its own internal margin); createGraphicsShapes draws to the
        // content byte directly, so the document margins are cosmetic -- 0, like the generic writeToPdf path
        document.setMargins( 0, 0, 0, 0 );
        PdfWriter writer;
        try {
            writer = PdfWriter.getInstance( document, new FileOutputStream( file_name ) );
        }
        catch ( final DocumentException e ) {
            throw new IOException( e );
        }
        document.open();
        final PdfContentByte cb = writer.getDirectContent();
        final Graphics2D g2 = cb.createGraphicsShapes( page_w, page_h );
        // Document-ready (light theme) export when requested, so a dark-theme figure isn't light-on-white on
        // the white PDF page; restored afterwards. Same behavior as raster/clipboard/SVG export. The paint
        // exception (if any) is reported AFTER restore(), so its modal dialog never pumps the EDT while the
        // transient light theme is applied.
        final ExportTheme export_theme = ExportTheme.applyIf( tree_panel, white_background );
        Exception paint_error = null;
        try {
            if ( white_background ) {
                g2.setColor( java.awt.Color.WHITE );
                g2.fillRect( 0, 0, page_w, page_h );
            }
            tree_panel.paintPhylogeny( g2, true, false, page_w, page_h, 0, 0 );
        }
        catch ( final Exception e ) {
            paint_error = e;
        }
        finally {
            // restore() is guarded on its own so a (near-impossible) failure there cannot skip document.close(),
            // which would leave a truncated/unreadable PDF and a leaked stream
            try {
                export_theme.restore();
            }
            catch ( final Exception e ) {
                //Do nothing.
            }
            try {
                g2.dispose();
                document.close();
            }
            catch ( final Exception e ) {
                //Do nothing.
            }
        }
        if ( paint_error != null ) {
            AptxUtil.unexpectedException( paint_error );
        }
        // A fixed-size caller restores the layout AFTER this returns, so the font size + dyna-hide state read for
        // the report still reflect the export layout.
        return AptxUtil.formatExportReport( file.toString(), tree_panel, page_w, page_h, true, 72, fixed_size );
    }
}
