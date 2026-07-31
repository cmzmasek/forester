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
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Integration test with teeth for the "Tree Name" canvas display toggle ({@code Options.isShowTreeName},
 * default on). Renders a named phylogram offscreen and pixel-scans the bottom band for the tree name, drawn
 * in the branch-length color (set to a unique magenta): (1) by default the name is drawn in the lower-LEFT;
 * (2) with the option off nothing is drawn there; (3) when the scale is shown (which occupies the lower-left)
 * the name slides to the lower-RIGHT so the two never overlap. Fails if the toggle stops gating the draw or if
 * the name no longer relocates when the scale is on. Headless-guarded (needs FlatLaf via createInstance).
 */
public final class ShowTreeNameRenderTest {

    private static final int    W    = 900;
    private static final int    H    = 600;
    private static final String NAME = "MYTREENAME";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ShowTreeNameRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phylogram( NAME ) }, conf,
                                                                        "n" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, Color.WHITE );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BRANCH, Color.BLACK );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BRANCH_LENGTH, new Color( 255, 0, 255 ) );
                    tp.getTreeColorSet().setColorSchema( 0 ); // propagate the scheme into the live color fields

                    if ( !mf[ 0 ].getOptions().isShowTreeName() ) {
                        fail( ok, "isShowTreeName should default to ON" );
                    }
                    // precondition: branch-length VALUE labels use the same (branch-length) color as the name,
                    // so keep them off, otherwise the magenta scan would not isolate the name
                    if ( tp.getControlPanel().isWriteBranchLengthValues() ) {
                        fail( ok, "test precondition: branch-length values should be off" );
                    }

                    // (1) default: name shown, scale off -> lower-LEFT (and NOT on the right)
                    mf[ 0 ].getOptions().setShowScale( false );
                    mf[ 0 ].getOptions().setShowTreeName( true );
                    final BufferedImage on = render( tp );
                    if ( magenta( on, 0, ( W / 2 ) - 1, H - 45, H - 6 ) < 4 ) {
                        fail( ok, "tree name should be drawn in the lower-left by default" );
                    }
                    if ( magenta( on, W / 2, W - 1, H - 45, H - 6 ) != 0 ) {
                        fail( ok, "with the scale off the name should not appear on the right" );
                    }

                    // (2) toggled off -> nothing in the name band
                    mf[ 0 ].getOptions().setShowTreeName( false );
                    final BufferedImage off = render( tp );
                    if ( magenta( off, 0, W - 1, H - 45, H - 6 ) != 0 ) {
                        fail( ok, "with the tree-name option off, nothing should be drawn in the name band" );
                    }

                    // (3) scale shown -> the name relocates to the lower-RIGHT while the scale keeps the
                    // lower-left. Rendered via the EXPORT path with an explicit width, so the right-alignment
                    // reference is a known W (deterministic, and it exercises the graphics-file branch a figure
                    // export actually uses) rather than the viewport-dependent screen width.
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    mf[ 0 ].getOptions().setPrintBlackAndWhite( false ); // keep the export name in the branch-length color
                    mf[ 0 ].getOptions().setShowScale( true );
                    mf[ 0 ].getOptions().setShowTreeName( true );
                    final BufferedImage exported = renderExport( tp );
                    if ( magenta( exported, W / 2, W - 1, H - 45, H - 6 ) < 4 ) {
                        fail( ok, "with the scale shown, the tree name should move to the lower-right (export path)" );
                    }
                    // positively confirm the scale actually rendered lower-left (so step 3's premise holds and the
                    // right-hand signal is the relocated name, not a missing scale)
                    if ( magenta( exported, 0, ( W / 2 ) - 1, H - 45, H - 6 ) < 4 ) {
                        fail( ok, "the scale should occupy the lower-left when shown" );
                    }

                    // (4) the menu/Settings checkbox actually drives the option (actionPerformed -> updateOptions)
                    final boolean before = mf[ 0 ].getOptions().isShowTreeName();
                    mf[ 0 ]._show_tree_name_cbmi.doClick();
                    if ( mf[ 0 ].getOptions().isShowTreeName() == before ) {
                        fail( ok, "clicking the Tree Name checkbox should flip the option via updateOptions" );
                    }
                    mf[ 0 ]._show_tree_name_cbmi.doClick();
                    if ( mf[ 0 ].getOptions().isShowTreeName() != before ) {
                        fail( ok, "clicking the Tree Name checkbox again should restore the option" );
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ShowTreeNameRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private static BufferedImage render( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, W, H );
        tp.paintPhylogeny( g, false, false, 0, 0, 0, 0 ); // screen path (getVisibleRect-based region)
        g.dispose();
        return img;
    }

    // the graphics-file export path: region is the explicit W x H rectangle at the origin, so the tree name
    // right-aligns to a known W (independent of the on-screen viewport width)
    private static BufferedImage renderExport( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, W, H );
        tp.paintPhylogeny( g, false, true, W, H, 0, 0 ); // to_graphics_file=true, region = full W x H
        g.dispose();
        return img;
    }

    private static int magenta( final BufferedImage img, final int x0, final int x1, final int y0, final int y1 ) {
        int count = 0;
        for( int y = Math.max( 0, y0 ); y <= Math.min( img.getHeight() - 1, y1 ); ++y ) {
            for( int x = Math.max( 0, x0 ); x <= Math.min( img.getWidth() - 1, x1 ); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, gg = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( r > 180 ) && ( b > 180 ) && ( gg < 110 ) ) {
                    ++count;
                }
            }
        }
        return count;
    }

    // root -> ( (A,B), C ) with branch lengths (so it is a phylogram with a positive scale distance)
    private static Phylogeny phylogram( final String name ) {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode ab = child( root, null, 0.3 );
        child( ab, "A", 0.2 );
        child( ab, "B", 0.2 );
        child( root, "C", 0.5 );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setName( name );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode child( final PhylogenyNode parent, final String n, final double dist ) {
        final PhylogenyNode node = new PhylogenyNode();
        if ( n != null ) {
            node.setName( n );
        }
        parent.addAsChild( node );
        node.setDistanceToParent( dist );
        return node;
    }

    private ShowTreeNameRenderTest() {
    }
}
