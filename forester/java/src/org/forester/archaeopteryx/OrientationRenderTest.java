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

import org.forester.archaeopteryx.Options.TREE_ORIENTATION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Verifies the root-on-top / root-on-bottom orientation on the demo tree (forester/demo/root-on-top.xml): the tree
 * is genuinely rotated (in ROOT_TOP the root is ABOVE every tip and the tips fan out HORIZONTALLY; ROOT_BOTTOM is
 * the mirror), hit-testing follows the rotation (findNode at a tip's on-screen position returns that tip -- proving
 * the R-inverse path), and a scaled raster export still renders non-empty at the right size (proving the orientation
 * transform COMPOSES with the export's own scale rather than replacing it). Headful; a green no-op when headless.
 * Dogfoods the demo.
 */
public final class OrientationRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "OrientationRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return orientationRendersOk();
    }

    private static boolean orientationRendersOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/root-on-top.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "orient" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final PhylogenyNode root = phy.getRoot();
                    final PhylogenyNode[] tips = phy.getExternalNodes().toArray( new PhylogenyNode[ 0 ] );
                    if ( tips.length < 4 ) {
                        fail( ok, "demo must have several tips" );
                        return;
                    }
                    final int w = 620, h = 520;

                    // the "Tip label angle" setting drives the tip-label tilt (independent of the display type -- it was
                    // auto before): VERTICAL -> +/-90deg, ANGLED -> +/-45deg, with the sign following the orientation so
                    // the labels always extend AWAY from the tree.
                    o.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    o.setTipLabelDirection( Options.TIP_LABEL_DIRECTION.VERTICAL );
                    if ( Math.abs( tp.tipLabelAngleForTest() - ( Math.PI / 2.0 ) ) > 1e-6 ) {
                        fail( ok, "VERTICAL tip labels must be 90deg in ROOT_TOP, got " + tp.tipLabelAngleForTest() );
                    }
                    o.setTipLabelDirection( Options.TIP_LABEL_DIRECTION.ANGLED );
                    if ( Math.abs( tp.tipLabelAngleForTest() - ( Math.PI / 4.0 ) ) > 1e-6 ) {
                        fail( ok, "ANGLED tip labels must be 45deg in ROOT_TOP, got " + tp.tipLabelAngleForTest() );
                    }
                    o.setTreeOrientation( TREE_ORIENTATION.ROOT_BOTTOM );
                    if ( Math.abs( tp.tipLabelAngleForTest() - ( -Math.PI / 4.0 ) ) > 1e-6 ) {
                        fail( ok, "ANGLED tip labels must be -45deg in ROOT_BOTTOM, got " + tp.tipLabelAngleForTest() );
                    }
                    o.setTipLabelDirection( Options.TIP_LABEL_DIRECTION.VERTICAL );
                    if ( Math.abs( tp.tipLabelAngleForTest() - ( -Math.PI / 2.0 ) ) > 1e-6 ) {
                        fail( ok, "VERTICAL tip labels must be -90deg in ROOT_BOTTOM, got " + tp.tipLabelAngleForTest() );
                    }
                    o.setTipLabelDirection( Options.TIP_LABEL_DIRECTION.VERTICAL ); // restore default for the rest

                    // ROOT_TOP: the root sits ABOVE every tip (smaller screen y), and the tips fan out HORIZONTALLY
                    layout( frame, tp, o, TREE_ORIENTATION.ROOT_TOP, w, h );
                    final double root_y_top = tp.screenPointFor( root ).y;
                    double min_tip_y = Double.MAX_VALUE, min_tip_x = Double.MAX_VALUE, max_tip_x = -Double.MAX_VALUE;
                    for ( final PhylogenyNode t : tips ) {
                        final java.awt.geom.Point2D.Double p = tp.screenPointFor( t );
                        min_tip_y = Math.min( min_tip_y, p.y );
                        min_tip_x = Math.min( min_tip_x, p.x );
                        max_tip_x = Math.max( max_tip_x, p.x );
                    }
                    if ( !( root_y_top < min_tip_y ) ) {
                        fail( ok, "ROOT_TOP: root (y=" + root_y_top + ") must be above every tip (min tip y="
                                + min_tip_y + ")" );
                    }
                    if ( !( ( max_tip_x - min_tip_x ) > ( w / 4.0 ) ) ) {
                        fail( ok, "ROOT_TOP: tips must fan out horizontally (x span=" + ( max_tip_x - min_tip_x )
                                + " over width " + w + ")" );
                    }

                    // hit-testing must follow the rotation: findNode at a tip's on-screen position returns that tip
                    final PhylogenyNode tip0 = tips[ 0 ];
                    final java.awt.geom.Point2D.Double sp = tp.screenPointFor( tip0 );
                    final PhylogenyNode hit = tp.findNode( (int) Math.round( sp.x ), (int) Math.round( sp.y ) );
                    if ( hit != tip0 ) {
                        fail( ok, "ROOT_TOP: hit-testing must follow the rotation -- findNode returned "
                                + ( hit == null ? "null" : hit.getName() ) + " (expected " + tip0.getName() + ")" );
                    }

                    // internal-node data placement: the internal label sits LEFT of the (vertical) branch and the
                    // support + length sit RIGHT of it. Toggling this text doesn't move the tip layout, so a data-on
                    // vs data-off render diff at an internal branch's midpoint isolates the label (left) and the
                    // support/length (right) from the (unchanged) tree structure.
                    final ControlPanel cp0 = tp.getControlPanel();
                    setInternalNames( phy, true );
                    cp0.setCheckbox( Configuration.display_internal_data, true ); // the internal label needs this on
                    if ( cp0.getWriteConfidenceCb() != null ) {
                        cp0.getWriteConfidenceCb().setSelected( true );
                    }
                    cp0.setCheckbox( Configuration.write_branch_length_values, true );
                    final BufferedImage data_on = render( frame, tp, o, TREE_ORIENTATION.ROOT_TOP, w, h );
                    final PhylogenyNode inode = centralInternalNode( tp, phy, w, h );
                    if ( inode == null ) {
                        fail( ok, "test setup: no central internal node for the data-placement check" );
                    }
                    else {
                        final double imx = ( inode.getParent().getXcoord() + inode.getXcoord() ) / 2.0;
                        final java.awt.geom.Point2D.Double bm = tp.screenPoint( imx, inode.getYcoord() );
                        final int bx = (int) Math.round( bm.x );
                        final int by = (int) Math.round( bm.y );
                        final int left_on = darkInRegion( data_on, bx - 60, by - 8, bx - 4, by + 8 );
                        final int right_on = darkInRegion( data_on, bx + 4, by - 8, bx + 60, by + 8 );
                        // data OFF (same layout): clear internal names + support + length
                        setInternalNames( phy, false );
                        if ( cp0.getWriteConfidenceCb() != null ) {
                            cp0.getWriteConfidenceCb().setSelected( false );
                        }
                        cp0.setCheckbox( Configuration.write_branch_length_values, false );
                        final BufferedImage data_off = render( frame, tp, o, TREE_ORIENTATION.ROOT_TOP, w, h );
                        final int left_off = darkInRegion( data_off, bx - 60, by - 8, bx - 4, by + 8 );
                        final int right_off = darkInRegion( data_off, bx + 4, by - 8, bx + 60, by + 8 );
                        if ( ( left_on - left_off ) < 15 ) {
                            fail( ok, "internal label not drawn to the LEFT of the vertical branch (dark px "
                                    + left_off + " -> " + left_on + ")" );
                        }
                        if ( ( right_on - right_off ) < 15 ) {
                            fail( ok, "support/length not drawn to the RIGHT of the vertical branch (dark px "
                                    + right_off + " -> " + right_on + ")" );
                        }
                    }

                    // pulse-halo dirty-rects must be recorded in DEVICE space (not logical) in a vertical orientation,
                    // else the pulse timer invalidates the wrong region and the halo never breathes. Set a found node,
                    // paint the SCREEN path (which records _found_halo_bounds), and assert a region covers the node's
                    // ON-SCREEN position.
                    o.setPulseFoundNodes( true );
                    final PhylogenyNode found = tips[ tips.length / 2 ];
                    final java.util.HashSet<Long> fset = new java.util.HashSet<>();
                    fset.add( found.getId() );
                    tp.setFoundNodes0( fset );
                    o.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage pulse_img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
                    final java.awt.Graphics2D pg = pulse_img.createGraphics();
                    tp.paintComponent( pg ); // screen path -> records _found_halo_bounds
                    pg.dispose();
                    final java.awt.geom.Point2D.Double fsp = tp.screenPointFor( found );
                    if ( tp.getFoundHaloBoundsCountForTest() == 0 ) {
                        fail( ok, "test setup: no pulse-halo regions recorded" );
                    }
                    else if ( !tp.foundHaloBoundContainsForTest( (int) Math.round( fsp.x ),
                            (int) Math.round( fsp.y ) ) ) {
                        fail( ok, "pulse-halo dirty-rect is not at the found node's on-screen position in vertical "
                                + "mode (recorded in logical space?)" );
                    }
                    tp.setFoundNodes0( null );

                    // ROOT_BOTTOM: the mirror -- the root sits BELOW every tip (larger screen y)
                    layout( frame, tp, o, TREE_ORIENTATION.ROOT_BOTTOM, w, h );
                    final double root_y_bot = tp.screenPointFor( root ).y;
                    double max_tip_y = -Double.MAX_VALUE;
                    for ( final PhylogenyNode t : tips ) {
                        max_tip_y = Math.max( max_tip_y, tp.screenPointFor( t ).y );
                    }
                    if ( !( root_y_bot > max_tip_y ) ) {
                        fail( ok, "ROOT_BOTTOM: root (y=" + root_y_bot + ") must be below every tip (max tip y="
                                + max_tip_y + ")" );
                    }

                    // the canvas transform R must ACTUALLY be applied to the BRANCH geometry (the checks above use
                    // screenPointFor / findNode, and labels ride withNodeTextFrame -- all of which apply R
                    // independently of the geometry paint). Compare ROOT_TOP vs ROOT_BOTTOM: both are vertical, so
                    // they share the identical layout (calcParametersForPainting swaps the budget for both), and they
                    // differ ONLY in R. Turn off overlays + all node text so the ONLY thing that can differ is the
                    // rotated branches: if R were built but never applied to g, both would render the same horizontal
                    // tree -> byte-identical (diff 0). Correct rotation -> the tree is flipped top<->bottom.
                    o.setShowZebraStripes( false );
                    o.setShowScaleAxis( false );
                    o.setShowScaleGrid( false );
                    o.setShowScale( false );
                    final ControlPanel cp = tp.getControlPanel();
                    cp.setCheckbox( Configuration.show_node_names, false );
                    cp.setCheckbox( Configuration.display_external_data, false );
                    cp.setCheckbox( Configuration.display_internal_data, false );
                    final BufferedImage top = render( frame, tp, o, TREE_ORIENTATION.ROOT_TOP, w, h );
                    final BufferedImage bottom = render( frame, tp, o, TREE_ORIENTATION.ROOT_BOTTOM, w, h );
                    if ( ( top == null ) || ( bottom == null ) || !hasInk( top ) || !hasInk( bottom ) ) {
                        fail( ok, "a comparison render produced no ink" );
                    }
                    else {
                        final double diff = pixelDiffFraction( top, bottom );
                        if ( diff < 0.01 ) {
                            fail( ok, "ROOT_TOP and ROOT_BOTTOM branch renders are nearly identical (diff=" + diff
                                    + ") -- the rotation transform was not applied to the canvas" );
                        }
                    }
                    // export composition: a 2x raster export in ROOT_TOP must render non-empty at the doubled size
                    // (proving R composes with the export scale via g.transform, not setTransform)
                    final BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 2, false );
                    if ( ( img == null ) || ( img.getWidth() != ( 2 * w ) ) || ( img.getHeight() != ( 2 * h ) ) ) {
                        fail( ok, "2x export should be " + ( 2 * w ) + "x" + ( 2 * h ) + ", got "
                                + ( img == null ? "null" : img.getWidth() + "x" + img.getHeight() ) );
                    }
                    else if ( !hasInk( img ) ) {
                        fail( ok, "2x ROOT_TOP export rendered no ink (blank image)" );
                    }

                    // the on-screen paint must NOT wrong-cull vertical branches when the tip-spread is zoomed past the
                    // viewport: the screen path applies a y-based cull that is device-space, but in a vertical
                    // orientation logical-y is the HORIZONTAL axis, so that cull drops on-screen branches (they vanish
                    // on X+ zoom-in). Compare the SCREEN paint to the never-culling EXPORT at the same size after
                    // zooming the tip-spread well past the viewport.
                    o.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    for ( int i = 0; i < 8; ++i ) {
                        cp.zoomInScreenX( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                    }
                    tp.setSize( w, h );
                    final BufferedImage screen = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
                    final java.awt.Graphics2D sg = screen.createGraphics();
                    sg.setColor( java.awt.Color.WHITE );
                    sg.fillRect( 0, 0, w, h );
                    tp.paintComponent( sg ); // SCREEN path (to_graphics_file=false) -> the cull is active
                    sg.dispose();
                    final BufferedImage exp = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int screen_dark = darkPixels( screen );
                    final int export_dark = darkPixels( exp );
                    if ( screen_dark < ( export_dark * 0.8 ) ) {
                        fail( ok, "screen paint wrong-culled vertical branches on zoom-in (screen dark=" + screen_dark
                                + " vs export dark=" + export_dark + ")" );
                    }

                    // aligned phylogram in a vertical orientation: the tip->label leaders must draw as vertical
                    // geometry (not slanted inside the 45deg label frame). Smoke: it renders non-empty and a tip is
                    // still hit-testable (a regression here -- e.g. drawing the leader in the tilted frame -- broke
                    // the layout badly enough that this exercises the aligned path end to end).
                    cp.setCheckbox( Configuration.show_node_names, true );
                    cp.setCheckbox( Configuration.display_external_data, true );
                    cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM );
                    final BufferedImage aligned = render( frame, tp, o, TREE_ORIENTATION.ROOT_TOP, w, h );
                    if ( ( aligned == null ) || !hasInk( aligned ) ) {
                        fail( ok, "aligned-phylogram ROOT_TOP rendered no ink" );
                    }
                    else {
                        final java.awt.geom.Point2D.Double asp = tp.screenPointFor( tips[ 0 ] );
                        if ( tp.findNode( (int) Math.round( asp.x ), (int) Math.round( asp.y ) ) != tips[ 0 ] ) {
                            fail( ok, "aligned-phylogram ROOT_TOP: tip hit-testing broke" );
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

    /** Applies the orientation, fits + sizes the panel, and renders once (discarding the image) so the node coords
     *  are assigned and the orientation transform R is (re)built during the paint. */
    private static void layout( final MainFrame frame, final TreePanel tp, final Options o,
                                final TREE_ORIENTATION orient, final int w, final int h ) {
        o.setTreeOrientation( orient );
        frame.showWhole();
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
    }

    /** Like {@link #layout} but returns the rendered image (scale 1). */
    private static BufferedImage render( final MainFrame frame, final TreePanel tp, final Options o,
                                         final TREE_ORIENTATION orient, final int w, final int h ) {
        o.setTreeOrientation( orient );
        frame.showWhole();
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        return AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
    }

    /** Fraction of sampled pixels that differ between two equally-sized images (1.0 if the sizes differ). */
    private static double pixelDiffFraction( final BufferedImage a, final BufferedImage b ) {
        if ( ( a.getWidth() != b.getWidth() ) || ( a.getHeight() != b.getHeight() ) ) {
            return 1.0;
        }
        int diff = 0, total = 0;
        for( int y = 0; y < a.getHeight(); ++y ) {
            for( int x = 0; x < a.getWidth(); ++x ) {
                ++total;
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    ++diff;
                }
            }
        }
        return ( total == 0 ) ? 0.0 : ( (double) diff / total );
    }

    /** Names (or clears) every non-root internal node, so the internal-label placement can be exercised. */
    private static void setInternalNames( final Phylogeny phy, final boolean set ) {
        int c = 0;
        for ( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = phy.iteratorPreorder(); it
                .hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.isRoot() ) {
                n.setName( set ? ( "clade" + ( ++c ) ) : "" );
            }
        }
    }

    /** The first non-root internal node whose incoming-branch midpoint (screen) is comfortably inside the image. */
    private static PhylogenyNode centralInternalNode( final TreePanel tp, final Phylogeny phy, final int w,
                                                      final int h ) {
        for ( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = phy.iteratorPreorder(); it
                .hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isExternal() || n.isRoot() ) {
                continue;
            }
            final double mx = ( n.getParent().getXcoord() + n.getXcoord() ) / 2.0;
            final java.awt.geom.Point2D.Double bm = tp.screenPoint( mx, n.getYcoord() ); // the branch's leg midpoint
            if ( ( bm.x > 90 ) && ( bm.x < ( w - 90 ) ) && ( bm.y > 40 ) && ( bm.y < ( h - 40 ) ) ) {
                return n;
            }
        }
        return null;
    }

    /** Count of clearly-dark pixels in the (clamped) rectangle. */
    private static int darkInRegion( final BufferedImage img, int x0, int y0, int x1, int y1 ) {
        x0 = Math.max( 0, x0 );
        y0 = Math.max( 0, y0 );
        x1 = Math.min( img.getWidth() - 1, x1 );
        y1 = Math.min( img.getHeight() - 1, y1 );
        int n = 0;
        for( int y = y0; y <= y1; ++y ) {
            for( int x = x0; x <= x1; ++x ) {
                final int rgb = img.getRGB( x, y );
                if ( ( ( ( rgb >> 16 ) & 0xFF ) < 150 ) && ( ( ( rgb >> 8 ) & 0xFF ) < 150 )
                        && ( ( rgb & 0xFF ) < 150 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Count of clearly-dark (branch/text ink) pixels -- a proxy for how much of the tree actually drew. */
    private static int darkPixels( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
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

    /** True if the image has any non-near-white pixel (i.e. the tree actually drew ink on the white export bg). */
    private static boolean hasInk( final BufferedImage img ) {
        for( int y = 0; y < img.getHeight(); y += 3 ) {
            for( int x = 0; x < img.getWidth(); x += 3 ) {
                final int rgb = img.getRGB( x, y );
                if ( ( ( ( rgb >> 16 ) & 0xFF ) < 230 ) || ( ( ( rgb >> 8 ) & 0xFF ) < 230 )
                        || ( ( rgb & 0xFF ) < 230 ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [OrientationRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [OrientationRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private OrientationRenderTest() {
    }
}
