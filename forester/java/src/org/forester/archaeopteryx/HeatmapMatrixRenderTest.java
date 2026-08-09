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
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the heat-map-matrix demo (forester/demo/heatmap-matrix.xml) with s1..s6 added as MATRIX annotation columns
 * and asserts the heat-map grid draws (many colored cells) and spans cool (low) to warm (high), and that a matrix
 * column exposes its color-scale legend. This is a render smoke: the SHARED-scale invariant itself (equal values in
 * different columns get equal colors) is asserted headlessly with real teeth in AnnotationColumnsTest.testMatrix --
 * warm+cool cells alone would also appear under per-column scaling, so they do not prove sharing. Headful; a green
 * no-op when headless. Dogfoods the demo.
 */
public final class HeatmapMatrixRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "HeatmapMatrixRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return matrixRendersOk();
    }

    private static boolean matrixRendersOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/heatmap-matrix.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "mtx" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM );
                    final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<>();
                    for( int i = 1; i <= 6; ++i ) {
                        specs.add( new AnnotationColumns.ColumnSpec( "data:s" + i, AnnotationColumns.Type.MATRIX ) );
                    }
                    tp.setAnnotationColumns( specs );
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    final int w = 640, h = 680;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int colored = countColored( img );
                    if ( colored < 2000 ) {
                        fail( ok, "the heat-map matrix should paint many colored cells, got " + colored );
                    }
                    // the grid spans the gradient: the global MAX value (9) reads warm (r>b) and the global MIN (0)
                    // cool (b>r). Count clearly-warm vs clearly-cool cells; a matrix on a real gradient must have both.
                    final int warm = countWarm( img ), cool = countCool( img );
                    if ( ( warm < 200 ) || ( cool < 200 ) ) {
                        fail( ok, "the matrix should span cool (low) to warm (high) cells, warm=" + warm
                                + " cool=" + cool );
                    }
                    // ALWAYS-ON matrix legend: a heat-map matrix shows its shared-scale legend with NO header click.
                    if ( tp.hasFocusedAnnotationColumn() ) {
                        fail( ok, "no annotation column was focused yet" );
                    }
                    if ( !tp.hasAnnotationColumnLegend() ) {
                        fail( ok, "a heat-map matrix must show its shared-scale legend WITHOUT a click (always-on)" );
                    }
                    // ...and it actually draws into the legend slot (records its draggable bounds + paints ink)
                    final BufferedImage leg = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
                    final java.awt.Graphics2D lg = leg.createGraphics();
                    lg.setColor( java.awt.Color.WHITE );
                    lg.fillRect( 0, 0, w, h );
                    tp.drawLegendForTest( lg, new java.awt.Rectangle( 0, 0, w, h ), true );
                    lg.dispose();
                    if ( tp.getPropertyLegendBounds() == null ) {
                        fail( ok, "the always-on matrix legend did not draw into the legend slot" );
                    }
                    // the always-on matrix legend DEFERS to an active "Color by" legend (they share the one slot, and
                    // the user explicitly turned Color-by on) -- it must not silently usurp it
                    tp.setColorByPropertyRef( "data:s1" );
                    if ( tp.hasAnnotationColumnLegend() ) {
                        fail( ok, "the matrix legend must defer to an active Color-by legend, not usurp the slot" );
                    }
                    tp.setColorByPropertyRef( null );
                    if ( !tp.hasAnnotationColumnLegend() ) {
                        fail( ok, "with Color-by off, the always-on matrix legend should return" );
                    }
                    // explicit focus still overrides to a specific column's legend (guards columnHasLegend(MATRIX))
                    tp.setFocusedAnnotationColumn( 0 );
                    if ( !tp.hasFocusedAnnotationColumn() ) {
                        fail( ok, "focusing a MATRIX column header must expose its legend (columnHasLegend missing MATRIX)" );
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

    private static int countColored( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int max = Math.max( r, Math.max( g, b ) ), min = Math.min( r, Math.min( g, b ) );
                if ( ( ( max - min ) > 60 ) && ( max > 100 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Clearly warm (red-dominant) cells -- the high end of the shared gradient. */
    private static int countWarm( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( r > 150 ) && ( r > ( b + 60 ) ) && ( r > ( g + 40 ) ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Clearly cool (blue-dominant) cells -- the low end of the shared gradient. */
    private static int countCool( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( b > 150 ) && ( b > ( r + 60 ) ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [HeatmapMatrixRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [HeatmapMatrixRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private HeatmapMatrixRenderTest() {
    }
}
