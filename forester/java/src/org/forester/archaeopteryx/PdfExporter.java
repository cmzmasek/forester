// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// Copyright (C) 2003-2007 Ethalinda K.S. Cannon
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx;

import java.awt.Graphics2D;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;

import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.awt.DefaultFontMapper.BaseFontParameters;
import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.FontFactory;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.BaseFont;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;

/*
 * 
 * This uses iText.
 * 
 * See: http://www.lowagie.com/iText/
 * 
 * Current version: iText-5.5.9
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
        final DefaultFontMapper mapper = new DefaultFontMapper();
        FontFactory.registerDirectories();
        if ( ForesterUtil.isWindows() ) {
            mapper.insertDirectory( "c:/windows/fonts" );
        }
        else if ( ForesterUtil.isMac() ) {
            mapper.insertDirectory( "/Library/Fonts/" );
            mapper.insertDirectory( "/System/Library/Fonts/" );
        }
        else {
            mapper.insertDirectory( "/usr/X/lib/X11/fonts/TrueType/" );
            mapper.insertDirectory( "/usr/X/lib/X11/fonts/Type1/" );
            mapper.insertDirectory( "/usr/share/fonts/default/TrueType/" );
            mapper.insertDirectory( "/usr/share/fonts/default/Type1/" );
        }
        enableUnicode( mapper );
        final PdfContentByte cb = writer.getDirectContent();
        
        final Graphics2D g2 = new PdfGraphics2D(cb, my_width, my_height, mapper); 
    
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

    private final static void enableUnicode( final DefaultFontMapper mapper ) {
        final Map<String, DefaultFontMapper.BaseFontParameters> map = mapper.getMapper();
        for (final Iterator<String> i = map.keySet().iterator(); i.hasNext();) {
            final String name = i.next();
            final String name_lc = name.toLowerCase();
            if ( name_lc.contains( "unicode" ) || name_lc.equals( "dialog" ) ) {
                final BaseFontParameters pfps = map.get(name);
                try {
                    pfps.encoding = BaseFont.IDENTITY_H;
                    pfps.embedded = true;
                }
                catch ( Exception e )  {
                    //Ignore.
                }
            }
        }
    }
    
    /* not used currently 
    static FontMapper arial_uni = new FontMapper() {
        public BaseFont awtToPdf(Font font) {
            System.out.println( font.toString() );
            try {
                return BaseFont.createFont(
                        "c:/windows/fonts/arialuni.ttf",
                        BaseFont.IDENTITY_H, BaseFont.EMBEDDED);
            }
            catch (DocumentException e) {
                e.printStackTrace();
            }
            catch (IOException e) {
                e.printStackTrace();
            }
            return null;
        }
       
        @Override
        public Font pdfToAwt( BaseFont arg0, int arg1 ) {
            return null;
        }
    };
    */
}
