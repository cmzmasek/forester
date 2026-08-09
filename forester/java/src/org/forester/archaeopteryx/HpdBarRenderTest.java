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

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the dated demo phylogram (forester/demo/node-hpd-bars.xml) as a phylogram with "Node Age Bars (HPD)" ON and
 * asserts translucent-blue bars appear (and none when off). Headful; a green no-op when headless. Dogfoods the demo.
 */
public final class HpdBarRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "HpdBarRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return barsRenderOk();
    }

    private static boolean barsRenderOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/node-hpd-bars.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "hpd" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true ); // predictable white background for the blue composite
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final int w = 900, h = 460;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    o.setShowHpdBars( false );
                    final int off = countBluish( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    o.setShowHpdBars( true );
                    final java.awt.image.BufferedImage on_img =
                            AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int on = countBluish( on_img );
                    if ( on <= ( off + 300 ) ) {
                        fail( ok, "Node Age Bars should add many blue pixels (on=" + on + " off=" + off + ")" );
                    }
                    if ( off >= 100 ) {
                        fail( ok, "no blue bars should appear when Node Age Bars is off, got " + off );
                    }
                    // the bar must STRADDLE its node (anchored to the node's own x): pick an internal, non-root dated
                    // node and assert blue pixels appear BOTH left and right of its x at its y
                    org.forester.phylogeny.PhylogenyNode inode = null;
                    for ( final java.util.Iterator<org.forester.phylogeny.PhylogenyNode> it = phy.iteratorPreorder();
                            it.hasNext(); ) {
                        final org.forester.phylogeny.PhylogenyNode n = it.next();
                        if ( !n.isExternal() && !n.isRoot() && n.getNodeData().isHasDate() ) {
                            inode = n;
                            break;
                        }
                    }
                    if ( inode == null ) {
                        fail( ok, "expected an internal dated node in the demo" );
                    }
                    else {
                        final int nx = Math.round( inode.getXcoord() ), ny = Math.round( inode.getYcoord() );
                        if ( !bluishNear( on_img, nx - 12, nx - 2, ny ) || !bluishNear( on_img, nx + 2, nx + 12, ny ) ) {
                            fail( ok, "the HPD bar must straddle its node's x (" + nx + "," + ny + ")" );
                        }
                    }
                    // VERTICAL PARITY: HPD bars are plain rects, so they ride the rotation R into vertical bars at the
                    // internal nodes in a root-top/bottom orientation. Render in ROOT_TOP and confirm the blue draws.
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    final int vertical_on = countBluish( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( vertical_on <= ( off + 300 ) ) {
                        fail( ok, "Node Age Bars should draw in a vertical orientation (on=" + vertical_on + " off="
                                + off + ")" );
                    }
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );

                    // collapsing a clade must REMOVE its hidden internal descendants' bars (not draw them at stale
                    // coords): find a dated clade with an internal child, collapse it, and assert the blue drops
                    org.forester.phylogeny.PhylogenyNode clade = null;
                    for ( final java.util.Iterator<org.forester.phylogeny.PhylogenyNode> it = phy.iteratorPreorder();
                            it.hasNext(); ) {
                        final org.forester.phylogeny.PhylogenyNode n = it.next();
                        if ( n.isExternal() || n.isRoot() || !n.getNodeData().isHasDate() ) {
                            continue;
                        }
                        for( int i = 0; i < n.getNumberOfDescendants(); ++i ) {
                            if ( !n.getChildNode( i ).isExternal() ) {
                                clade = n;
                                break;
                            }
                        }
                        if ( clade != null ) {
                            break;
                        }
                    }
                    if ( clade != null ) {
                        tp.collapse( clade );
                        tp.calcParametersForPainting( w, h );
                        final int collapsed = countBluish(
                                AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                        if ( collapsed >= ( ( on * 3 ) / 4 ) ) {
                            fail( ok, "collapsing a clade must remove its hidden descendants' HPD bars (collapsed="
                                    + collapsed + " expanded=" + on + ")" );
                        }
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

    /** Pixels where blue clearly dominates red and green -- the translucent-blue HPD bars (over white OR over the
     *  black branches), but not white background, black branches, or gray antialiased text. */
    private static int countBluish( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( b >= ( r + 20 ) ) && ( b >= ( g + 15 ) ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Is there any bluish pixel in the band [x0,x1] x [y-4,y+4]? (the bar is HPD_BAR_HEIGHT=7 tall, centred on y) */
    private static boolean bluishNear( final BufferedImage img, final int x0, final int x1, final int y ) {
        for( int yy = Math.max( 0, y - 4 ); yy < Math.min( img.getHeight(), y + 5 ); ++yy ) {
            for( int xx = Math.max( 0, x0 ); xx < Math.min( img.getWidth(), x1 + 1 ); ++xx ) {
                final int rgb = img.getRGB( xx, yy );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( b >= ( r + 20 ) ) && ( b >= ( g + 15 ) ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [HpdBarRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [HpdBarRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private HpdBarRenderTest() {
    }
}
