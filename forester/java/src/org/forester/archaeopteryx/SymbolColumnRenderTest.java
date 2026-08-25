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
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the symbol-columns demo (forester/demo/symbol-columns.xml) with SYMBOL annotation columns and asserts:
 * (a) the shape glyphs draw (an all-present categorical column adds colored ink over the plain tree); (b) a
 * mixed binary column (some FILLED, some hollow OUTLINE, one missing/NONE) draws clearly LESS ink than the
 * all-filled column -- so a hollow or missing value is NOT painted as a filled disc; and (c) focusing a SYMBOL
 * column exposes its categorical color-key legend. Headful; a green no-op when headless. Dogfoods the demo.
 */
public final class SymbolColumnRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SymbolColumnRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return symbolRendersOk();
    }

    private static boolean symbolRendersOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/symbol-columns.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "sym" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true ); // gray tree on white -> only the marks are colored
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM );
                    final int w = 700, h = 500;

                    // baseline: no annotation columns
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final int base = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );

                    // 'host' as a SYMBOL column: all 8 tips carry a value -> all FILLED discs
                    tp.setAnnotationColumns(
                            Arrays.asList( new AnnotationColumns.ColumnSpec( "data:host", AnnotationColumns.Type.SYMBOL ) ) );
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final int host_ink = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( host_ink <= ( base + 50 ) ) {
                        fail( ok, "the symbol column drew no shape ink (host=" + host_ink + ", base=" + base + ")" );
                    }
                    if ( tp.annotationColumnWidthForTest( 0 ) <= 0 ) {
                        fail( ok, "the symbol column reserved no width" );
                    }
                    // focusing a SYMBOL column exposes its categorical color key
                    tp.setFocusedAnnotationColumn( 0 );
                    if ( !tp.hasFocusedAnnotationColumn() ) {
                        fail( ok, "focusing a SYMBOL column must expose its legend (columnHasLegend missing SYMBOL?)" );
                    }
                    final BufferedImage leg = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
                    final Graphics2D lg = leg.createGraphics();
                    lg.setColor( Color.WHITE );
                    lg.fillRect( 0, 0, w, h );
                    tp.drawLegendForTest( lg, new Rectangle( 0, 0, w, h ), true );
                    lg.dispose();
                    if ( tp.getPropertyLegendBounds() == null ) {
                        fail( ok, "the focused SYMBOL legend did not draw into the legend slot" );
                    }
                    tp.setFocusedAnnotationColumn( -1 );

                    // 'resistant' as a SYMBOL column: FILLED (yes) + hollow OUTLINE (no) + one missing (NONE) tip.
                    // It must draw SOME ink (the yes/no marks) but clearly LESS than the all-filled host column --
                    // if a hollow or missing value were painted as a filled disc, its ink would climb toward host's.
                    tp.setAnnotationColumns(
                            Arrays.asList( new AnnotationColumns.ColumnSpec( "data:resistant", AnnotationColumns.Type.SYMBOL ) ) );
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final int res_ink = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( res_ink <= ( base + 20 ) ) {
                        fail( ok, "the binary symbol column drew no marks (res=" + res_ink + ", base=" + base + ")" );
                    }
                    if ( host_ink <= ( res_ink * 1.3 ) ) {
                        fail( ok, "a hollow/missing symbol should draw less ink than an all-filled column: host="
                                + host_ink + " res=" + res_ink );
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

    /** Saturated (non-gray, colored) pixels -- the annotation marks; the tree + labels render gray on white. */
    private static int countSaturated( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( ( Math.max( r, Math.max( g, b ) ) - Math.min( r, Math.min( g, b ) ) ) > 40 )
                        && ( Math.max( r, Math.max( g, b ) ) > 80 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [SymbolColumnRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [SymbolColumnRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private SymbolColumnRenderTest() {
    }
}
