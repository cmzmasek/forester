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
 * Verifies the clustergram "Tip Labels Below Columns" layout on the annotation-columns demo: in a root-top vertical
 * orientation with tip-aligned columns, turning the option ON moves the tip labels from between the tree and the
 * columns down to BELOW the columns (so the dendrogram sits directly on the grid). Asserts label ink appears past the
 * colored column band only when the option is on. Headful; a green no-op when headless. Dogfoods the demo.
 */
public final class ClustergramRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ClustergramRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return clustergramOk();
    }

    private static boolean clustergramOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/annotation-columns.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "clust" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    tp.setColorByPropertyRef( null ); // establish this test's premise: no Color-by active (undo the load-time auto-color)
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true );
                    o.setShowTreeName( false ); // else the lower-left tree name is dark text below the columns too
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<>();
                    specs.add( new AnnotationColumns.ColumnSpec( "data:host", AnnotationColumns.Type.COLOR_STRIP ) );
                    specs.add( new AnnotationColumns.ColumnSpec( "data:viral_load", AnnotationColumns.Type.HEATMAP ) );
                    tp.setAnnotationColumns( specs );
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    final int w = 720, h = 720;

                    o.setTipLabelsBelowColumns( false );
                    final BufferedImage off = render( frame, tp, o, w, h );
                    o.setTipLabelsBelowColumns( true );
                    final BufferedImage on = render( frame, tp, o, w, h );

                    // the sample labels must appear BELOW the colored column band only when the option is on. Find the
                    // lowest colored (column) row in each image, then count dark (text) ink in the strip just below it.
                    final int off_bottom = lowestColoredRow( off );
                    final int on_bottom = lowestColoredRow( on );
                    if ( ( off_bottom < 0 ) || ( on_bottom < 0 ) ) {
                        fail( ok, "the annotation columns did not render (no colored band found)" );
                    }
                    else {
                        final int below_on = darkText( on, on_bottom + 3, on_bottom + 60 );
                        final int below_off = darkText( off, off_bottom + 3, off_bottom + 60 );
                        if ( below_on < 300 ) {
                            fail( ok, "clustergram: tip labels should draw BELOW the columns when the option is on, "
                                    + "got " + below_on + " dark px" );
                        }
                        if ( below_off > ( below_on / 3 ) ) {
                            fail( ok, "without the clustergram layout, the tip labels should NOT be below the columns "
                                    + "(off=" + below_off + " on=" + below_on + ")" );
                        }
                    }
                    // the option must be inert in the HORIZONTAL orientation (no columns axis to sit below)
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    if ( tp.tipLabelsBelowColumns() ) {
                        fail( ok, "Tip Labels Below Columns must be a no-op in the horizontal orientation" );
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

    private static BufferedImage render( final MainFrame frame, final TreePanel tp, final Options o, final int w,
                                         final int h ) {
        frame.showWhole();
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        return AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
    }

    /** The lowest image row that still carries a run of clearly-colored (non-gray, saturated) pixels -- the bottom of
     *  the annotation-column band. -1 if none. */
    private static int lowestColoredRow( final BufferedImage img ) {
        for( int y = img.getHeight() - 1; y >= 0; --y ) {
            int colored = 0;
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int max = Math.max( r, Math.max( g, b ) ), min = Math.min( r, Math.min( g, b ) );
                if ( ( ( max - min ) > 60 ) && ( max > 100 ) ) {
                    ++colored;
                }
            }
            if ( colored > 40 ) {
                return y;
            }
        }
        return -1;
    }

    /** Count of clearly-dark (text ink) pixels in the row band [y0, y1). */
    private static int darkText( final BufferedImage img, final int y0, final int y1 ) {
        int n = 0;
        for( int y = Math.max( 0, y0 ); y < Math.min( img.getHeight(), y1 ); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                if ( ( ( ( rgb >> 16 ) & 0xFF ) < 120 ) && ( ( ( rgb >> 8 ) & 0xFF ) < 120 )
                        && ( ( rgb & 0xFF ) < 120 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ClustergramRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ClustergramRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private ClustergramRenderTest() {
    }
}
