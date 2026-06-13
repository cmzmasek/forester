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

    static String writePhylogenyToPdf( final String file_name, final TreePanel tree_panel, final int width, final int height )
            throws IOException {
        final int my_height;
        final int my_width;
        if ( height < HEIGHT_LIMIT ) {
            my_height = HEIGHT_LIMIT + 2 * MARGIN_Y;
        }
        else {
            my_height = height + 2 * MARGIN_Y;
        }
        if ( width < WIDTH_LIMIT ) {
            my_width = WIDTH_LIMIT +  2 * MARGIN_X;
        }
        else {
            my_width = width +  2 * MARGIN_X;
        }
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
        document.setPageSize( new Rectangle( my_width, my_height ) );
        document.setMargins( MARGIN_X, MARGIN_X, MARGIN_Y, MARGIN_Y );
        PdfWriter writer = null;
        try {
            writer = PdfWriter.getInstance( document, new FileOutputStream( file_name ) );
           
        }
        catch ( final DocumentException e ) {
            throw new IOException( e );
        }
        document.open();
        // Text is rendered as vector outlines (createGraphicsShapes), so the figure needs no font
        // embedding/mapping and is fully portable; it also sidesteps the bold-glyph stroke-color
        // bleed that the glyph-font path has in iText/OpenPDF.
        final PdfContentByte cb = writer.getDirectContent();
        final Graphics2D g2 = cb.createGraphicsShapes(my_width, my_height);
    
        try {
            tree_panel.paintPhylogeny( g2, true, false, my_width, my_height, 0, 0 );
        }
        catch ( final Exception e ) {
            AptxUtil.unexpectedException( e );
        }
        finally {
            try {
                g2.dispose();
                document.close();
            }
            catch ( final Exception e ) {
                //Do nothing.
            }
        }
        final String msg = file.toString() +  " [size: " + my_width + ", " + my_height + "]";
        return msg;
    }
}
