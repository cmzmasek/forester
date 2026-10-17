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

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.archaeopteryx.Options.SUPPORT_VISUALIZATION;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;

/**
 * End-to-end render test for the node-symbol support visualization
 * ({@link TreePanel#paintNodeSupportSymbol}). The per-value math is unit-tested headlessly in
 * {@link TreePanelUtilTest}; this drives a real {@link TreePanel} through the PNG export path and
 * checks that the symbols actually appear on internal nodes and that the modes differ as designed.
 *
 * <p>Needs FlatLaf + a display, so {@link #test()} is a no-op (returns true) when headless -- the
 * same pattern the other GUI tests in this package use; it runs for real from {@link #main} or any
 * non-headless invocation.
 *
 * <p>Tree: two internal nodes under the root carry bootstrap support 99 and 90; a third (60) sits
 * deeper. With the default 0.95 cutoff, THRESHOLD_MARKS marks only the 99 node, SIZE_SCALED draws a
 * dot at all three (sized by support), and NONE draws none -- so dark-pixel counts must satisfy
 * none &lt; threshold &lt; size-scaled.
 */
public final class NodeSupportSymbolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeSupportSymbol: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = supportTree();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "support test" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final JFrame f = (JFrame) mf[ 0 ];
                try {
                    f.setSize( 900, 600 );
                    f.validate();
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    mp.getOptions().setRasterExportScale( 1 ); // pixel checks below are in logical coordinates
                    mp.getControlPanel().showWhole();
                    if ( ( tp.getWidth() < 200 ) || ( tp.getHeight() < 200 ) ) {
                        return; // no usable viewport in this environment; nothing to assert
                    }
                    mp.getOptions().setSupportThreshold( 0.95 );
                    final long none = renderDarkPixels( tp, mp, SUPPORT_VISUALIZATION.NONE );
                    final long threshold = renderDarkPixels( tp, mp, SUPPORT_VISUALIZATION.THRESHOLD_MARKS );
                    final long size = renderDarkPixels( tp, mp, SUPPORT_VISUALIZATION.SIZE_SCALED );
                    if ( !( none < threshold ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  threshold marks added no symbol (none=" + none + ", threshold="
                                + threshold + ")" );
                    }
                    if ( !( threshold < size ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  size-scaled should mark more nodes than the 0.95 threshold (threshold="
                                + threshold + ", size=" + size + ")" );
                    }
                    // the support dot belongs on the middle of the branch (a branch property), not the node
                    PhylogenyNode i90 = null;
                    for( final PhylogenyNode ext : phy.getExternalNodes() ) {
                        if ( "D".equals( ext.getName() ) ) {
                            i90 = ext.getParent(); // ((A,(B,C)60)99,(D,E)90)root -> parent of D is the support-90 node
                            break;
                        }
                    }
                    if ( ( i90 != null )
                            && !dotIsMidBranch( renderImage( tp, mp, SUPPORT_VISUALIZATION.SIZE_SCALED ), i90 ) ) {
                        ok[ 0 ] = TestFail.here();
                    }
                    // the threshold mark is a yes/no mark at the default node size, smaller than the size-scaled
                    // dot of the same (near-full) support, whose diameter IS the value (8 px until 0.11.184 for both;
                    // Christian, 2026-10-16: "too large")
                    PhylogenyNode i99 = null;
                    for( final PhylogenyNode ext : phy.getExternalNodes() ) {
                        if ( "A".equals( ext.getName() ) ) {
                            i99 = ext.getParent();
                            break;
                        }
                    }
                    if ( i99 != null ) {
                        final int mark = verticalExtent( renderImage( tp, mp, SUPPORT_VISUALIZATION.THRESHOLD_MARKS ), i99 );
                        final int scaled = verticalExtent( renderImage( tp, mp, SUPPORT_VISUALIZATION.SIZE_SCALED ), i99 );
                        // what a disc of the nominal diameter measures when drawn the painter's way at the same
                        // fractional position (an antialiased disc measures a pixel under its diameter: a +-1 check
                        // let a 7 px mutant pass, 2026-10-16), so the mark must measure exactly that
                        final int want = verticalExtent( referenceDisc( AptxConstants.SUPPORT_THRESHOLD_MARK_DIAMETER, i99 ), i99 );
                        if ( mark != want ) {
                            ok[ 0 ] = false;
                            System.out.println( "  the threshold mark must measure as a " + AptxConstants.SUPPORT_THRESHOLD_MARK_DIAMETER
                                    + " px disc does (" + want + " px), measured " + mark );
                        }
                        if ( !( mark < scaled ) ) {
                            ok[ 0 ] = false;
                            System.out.println( "  the threshold mark (" + mark + " px) must be smaller than the size-scaled dot of the same support ("
                                    + scaled + " px)" );
                        }
                    }
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                    ok[ 0 ] = TestFail.here();
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

    private static long renderDarkPixels( final TreePanel tp, final MainPanel mp, final SUPPORT_VISUALIZATION mode )
            throws Exception {
        final BufferedImage img = renderImage( tp, mp, mode );
        long dark = 0;
        for( int y = 0; y < img.getHeight(); y++ ) {
            for( int x = 0; x < img.getWidth(); x++ ) {
                if ( isDark( img, x, y ) ) {
                    dark++;
                }
            }
        }
        return dark;
    }

    private static BufferedImage renderImage( final TreePanel tp, final MainPanel mp,
                                              final SUPPORT_VISUALIZATION mode ) throws Exception {
        mp.getOptions().setSupportVisualization( mode );
        final File out = File.createTempFile( "aptx_support_" + mode, ".png" );
        try {
            AptxUtil.writePhylogenyToGraphicsFile( out.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                   mp.getControlPanel(), GraphicsExportType.PNG, mp.getOptions() );
            return ImageIO.read( out );
        }
        finally {
            out.delete();
        }
    }

    private static boolean isDark( final BufferedImage img, final int x, final int y ) {
        if ( ( x < 0 ) || ( y < 0 ) || ( x >= img.getWidth() ) || ( y >= img.getHeight() ) ) {
            return false;
        }
        final int rgb = img.getRGB( x, y );
        return ( ( ( ( rgb >> 16 ) & 0xFF ) + ( ( rgb >> 8 ) & 0xFF ) + ( rgb & 0xFF ) ) < 200 );
    }

    /** Dark pixels in the (2r+1)x(2r+1) window centered at (cx, cy). */
    private static long darkInWindow( final BufferedImage img, final int cx, final int cy, final int r ) {
        long n = 0;
        for( int y = cy - r; y <= ( cy + r ); ++y ) {
            for( int x = cx - r; x <= ( cx + r ); ++x ) {
                if ( isDark( img, x, y ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /**
     * Verifies the support dot is drawn on the MIDDLE of an internal branch, not at the node: a window at
     * the branch midpoint must hold clearly more dark pixels (the filled dot, on top of the branch line)
     * than a same-size window on the branch just left of the node (the bare line). Returns true if OK.
     */
    private static boolean dotIsMidBranch( final BufferedImage img, final PhylogenyNode internal ) {
        final int by = Math.round( internal.getYcoord() );
        final int parent_x = Math.round( internal.getParent().getXcoord() );
        final int node_x = Math.round( internal.getXcoord() );
        final int mid_x = ( parent_x + node_x ) / 2;
        final long mid_dark = darkInWindow( img, mid_x, by, 6 );
        final long node_dark = darkInWindow( img, node_x - 9, by, 6 ); // on the branch, clear of the node connector
        final boolean ok = mid_dark > ( node_dark + 15 );
        if ( !ok ) {
            System.out.println( "  support dot not at branch middle (mid_dark=" + mid_dark + ", node_dark=" + node_dark
                    + ", parent_x=" + parent_x + ", node_x=" + node_x + ")" );
        }
        return ok;
    }

    /** A disc of {@code diameter} drawn as the painter draws a support symbol -- an antialiased filled ellipse framed
     *  at (centre - d/2) -- on the branch centre of {@code internal}, alone on a white image of the panel's size. */
    private static BufferedImage referenceDisc( final float diameter, final PhylogenyNode internal ) {
        final PhylogenyNode parent = internal.getParent();
        final float[] c = TreePanelUtil.supportSymbolCenter( parent.getXcoord(), internal.getXcoord(), parent.getYcoord(),
                                                            internal.getYcoord(), false );
        final BufferedImage img = new BufferedImage( Math.round( c[ 0 ] ) + 40, Math.round( c[ 1 ] ) + 40, BufferedImage.TYPE_INT_RGB );
        final java.awt.Graphics2D g = img.createGraphics();
        g.setColor( java.awt.Color.WHITE );
        g.fillRect( 0, 0, img.getWidth(), img.getHeight() );
        g.setRenderingHint( java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON );
        g.setColor( java.awt.Color.BLACK );
        final float half = diameter / 2.0f;
        g.fill( new java.awt.geom.Ellipse2D.Double( c[ 0 ] - half, c[ 1 ] - half, diameter, diameter ) );
        g.dispose();
        return img;
    }

    /** The height of the dark run through the support symbol's centre on {@code internal}'s branch: the symbol's
     *  diameter (the branch line under it is horizontal, so it adds nothing to a VERTICAL run). */
    private static int verticalExtent( final BufferedImage img, final PhylogenyNode internal ) {
        final PhylogenyNode parent = internal.getParent();
        final float[] c = TreePanelUtil.supportSymbolCenter( parent.getXcoord(), internal.getXcoord(), parent.getYcoord(),
                                                            internal.getYcoord(), false );
        final int cx = Math.round( c[ 0 ] ), cy = Math.round( c[ 1 ] );
        if ( !isDark( img, cx, cy ) ) {
            System.out.println( "  no symbol at the branch centre " + cx + "," + cy );
            return 0;
        }
        int top = cy, bottom = cy;
        while ( isDark( img, cx, top - 1 ) ) {
            --top;
        }
        while ( isDark( img, cx, bottom + 1 ) ) {
            ++bottom;
        }
        return ( bottom - top ) + 1;
    }

    /** ((A,(B,C)60)99,(D,E)90)root -- bootstrap support on the internal branches. */
    private static Phylogeny supportTree() {
        final Phylogeny p = new Phylogeny();
        final PhylogenyNode i60 = internal( 60.0, leaf( "B" ), leaf( "C" ) );
        final PhylogenyNode i99 = internal( 99.0, leaf( "A" ), i60 );
        final PhylogenyNode i90 = internal( 90.0, leaf( "D" ), leaf( "E" ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( i99 );
        root.addAsChild( i90 );
        p.setRoot( root );
        p.externalNodesHaveChanged();
        return p;
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode internal( final double conf, final PhylogenyNode a, final PhylogenyNode b ) {
        final PhylogenyNode i = new PhylogenyNode();
        i.getBranchData().addConfidence( new Confidence( conf, "bootstrap" ) );
        i.addAsChild( a );
        i.addAsChild( b );
        return i;
    }

    private NodeSupportSymbolTest() {
    }
}
