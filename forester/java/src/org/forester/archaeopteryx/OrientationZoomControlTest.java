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
import java.io.File;

import javax.swing.JFrame;
import javax.swing.JScrollBar;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.TREE_ORIENTATION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Verifies the orientation-sensitive zoom controls: the fit button relabels W&lt;-&gt;H with the orientation, and the
 * screen-oriented zoom flips depth&lt;-&gt;breadth in a vertical orientation (so the "X"/horizontal action drives the
 * tip-spread and "Y"/vertical drives the depth axis). Headful; a green no-op when headless.
 */
public final class OrientationZoomControlTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "OrientationZoomControl: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return zoomControlsOk();
    }

    private static boolean zoomControlsOk() {
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
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "zoom" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final ControlPanel cp = frame.getMainPanel().getControlPanel();
                    final Options o = frame.getOptions();
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final int w = 620, h = 520;

                    // the fit button's GLYPH tracks the orientation: a landscape frame with horizontal arrows
                    // (fit-width) in root-left, a portrait frame with vertical arrows (fit-height) in root-top
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_LEFT );
                    cp.updateZoomButtonsForLayout();
                    if ( cp.getFitButtonIconKind() != ControlButtonIcon.Kind.FIT_WIDTH ) {
                        fail( ok, "ROOT_LEFT fit button should carry FIT_WIDTH, got " + cp.getFitButtonIconKind() );
                    }
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    cp.updateZoomButtonsForLayout();
                    if ( cp.getFitButtonIconKind() != ControlButtonIcon.Kind.FIT_HEIGHT ) {
                        fail( ok, "ROOT_TOP fit button should carry FIT_HEIGHT, got " + cp.getFitButtonIconKind() );
                    }
                    // ... and the EXPAND glyph turns with it: the label rows it spreads run horizontally in
                    // root-left (spread vertically) and vertically in root-top (spread horizontally)
                    if ( !( cp.getExpandButtonForTest().getIcon() instanceof ControlButtonIcon )
                            || ( ( ( ControlButtonIcon ) cp.getExpandButtonForTest().getIcon() )
                                    .getKind() != ControlButtonIcon.Kind.EXPAND_HORIZONTAL ) ) {
                        fail( ok, "ROOT_TOP expand button should carry EXPAND_HORIZONTAL" );
                    }
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_LEFT );
                    cp.updateZoomButtonsForLayout();
                    if ( !( cp.getExpandButtonForTest().getIcon() instanceof ControlButtonIcon )
                            || ( ( ( ControlButtonIcon ) cp.getExpandButtonForTest().getIcon() )
                                    .getKind() != ControlButtonIcon.Kind.EXPAND_VERTICAL ) ) {
                        fail( ok, "ROOT_LEFT expand button should carry EXPAND_VERTICAL" );
                    }
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    cp.updateZoomButtonsForLayout();

                    // screen-oriented zoom flips depth<->breadth in a vertical orientation: the "X" (horizontal-screen)
                    // action must grow the tip-spread (y-distance), and "Y" (vertical-screen) must grow the depth
                    // (x-correction factor) -- the OPPOSITE of the horizontal layout
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    layout( frame, tp, w, h );
                    final float yd0 = tp.getYdistance();
                    cp.zoomInScreenX( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                            AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                    if ( !( tp.getYdistance() > yd0 ) ) {
                        fail( ok, "vertical: zoomInScreenX must grow the tip-spread (y-distance " + yd0 + " -> "
                                + tp.getYdistance() + ")" );
                    }
                    final float xc0 = tp.getXcorrectionFactor();
                    cp.zoomInScreenY( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                            AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                    if ( !( tp.getXcorrectionFactor() > xc0 ) ) {
                        fail( ok, "vertical: zoomInScreenY must grow the depth (x-correction " + xc0 + " -> "
                                + tp.getXcorrectionFactor() + ")" );
                    }

                    // in the horizontal layout it is the other way round: zoomInScreenX grows the depth
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_LEFT );
                    layout( frame, tp, w, h );
                    final float xcL = tp.getXcorrectionFactor();
                    cp.zoomInScreenX( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                            AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                    if ( !( tp.getXcorrectionFactor() > xcL ) ) {
                        fail( ok, "horizontal: zoomInScreenX must grow the depth (x-correction " + xcL + " -> "
                                + tp.getXcorrectionFactor() + ")" );
                    }

                    // realize the frame so the scroll bars / viewport are live (needed for re-center + fitHeight)
                    ( (JFrame) frame ).setSize( 500, 400 );
                    ( (JFrame) frame ).setVisible( true );

                    // fitHeight ("H") must be IDEMPOTENT: repeated presses fit the depth to the window height and keep
                    // the breadth (tip-spread) zoom -- they must not drift the y-distance (they did when the padded
                    // preferred width was fed back as the breadth budget).
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    frame.showWhole();
                    cp.fitHeight();
                    frame.getMainPanel().getCurrentScrollPane().validate();
                    final float yd_fit = tp.getYdistance();
                    cp.fitHeight();
                    cp.fitHeight();
                    frame.getMainPanel().getCurrentScrollPane().validate();
                    if ( Math.abs( tp.getYdistance() - yd_fit ) > ( ( yd_fit * 0.05f ) + 0.5f ) ) {
                        fail( ok, "fitHeight must not drift the breadth zoom (y-distance " + yd_fit + " -> "
                                + tp.getYdistance() + ")" );
                    }

                    // re-centering: X+ (breadth) zoom grows the HORIZONTAL extent, so it must re-center the HORIZONTAL
                    // scroll bar (keep the tree centered) -- in BOTH vertical orientations, not just ROOT_TOP.
                    for ( final TREE_ORIENTATION orient : new TREE_ORIENTATION[] { TREE_ORIENTATION.ROOT_TOP,
                            TREE_ORIENTATION.ROOT_BOTTOM } ) {
                        tp.setTreeOrientation( orient );
                        frame.showWhole();
                        for ( int i = 0; i < 7; ++i ) {
                            cp.zoomInScreenX( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                    AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                        }
                        frame.getMainPanel().getCurrentScrollPane().validate();
                        final JScrollBar hb = frame.getMainPanel().getCurrentScrollPane().getHorizontalScrollBar();
                        if ( ( hb.getMaximum() - hb.getMinimum() ) <= hb.getVisibleAmount() ) {
                            fail( ok, orient + ": horizontal scrollbar not active after zooming the tip-spread" );
                        }
                        else {
                            final double cf = centerFraction( hb );
                            if ( Math.abs( cf - 0.5 ) > 0.2 ) {
                                fail( ok, orient + ": X+ zoom did not re-center on the tree (horizontal center "
                                        + "fraction " + cf + ")" );
                            }
                        }
                    }

                    // centerOnNode must use the node's ON-SCREEN (rotated) position in a vertical orientation, so
                    // centering actually brings the node INTO the viewport near its middle -- not to a wrong spot
                    // derived from the raw logical coords.
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    frame.showWhole();
                    for ( int i = 0; i < 8; ++i ) {
                        cp.zoomInScreenX( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                        cp.zoomInScreenY( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                    }
                    frame.getMainPanel().getCurrentScrollPane().validate();
                    final java.util.List<org.forester.phylogeny.PhylogenyNode> exts = phy.getExternalNodes();
                    final org.forester.phylogeny.PhylogenyNode mid = exts.get( exts.size() / 2 );
                    tp.centerOnNodeForTest( mid );
                    frame.getMainPanel().getCurrentScrollPane().validate();
                    final java.awt.Rectangle vr = frame.getMainPanel().getCurrentScrollPane().getViewport()
                            .getViewRect();
                    final java.awt.geom.Point2D.Double sp = tp.screenPointFor( mid );
                    final int nx = (int) Math.round( sp.x ), ny = (int) Math.round( sp.y );
                    if ( !vr.contains( nx, ny ) ) {
                        fail( ok, "vertical centerOnNode did not bring the node into the viewport (node at " + nx + ","
                                + ny + " outside " + vr + ")" );
                    }
                    else {
                        // on the breadth (horizontal) axis a mid-tree tip is interior, so unless clamped at an edge it
                        // should land near the viewport's horizontal center
                        final boolean clamped_x = ( nx < ( vr.width / 2 ) )
                                || ( nx > ( tp.getWidth() - ( vr.width / 2 ) ) );
                        final int off_center = Math.abs( ( nx - vr.x ) - ( vr.width / 2 ) );
                        if ( !clamped_x && ( off_center > ( vr.width / 4 ) ) ) {
                            fail( ok, "vertical centerOnNode did not horizontally center the node (off-center "
                                    + off_center + "px of " + vr.width + ")" );
                        }
                    }

                    // the orientation transform R is cached across plain repaints and rebuilt only when the layout or
                    // orientation changes -- so a hover/pulse/scroll repaint does not re-walk the tree just to build R.
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    tp.calcParametersForPainting( w, h );
                    tp.rebuildOrientationTransform();
                    final java.awt.geom.AffineTransform r1 = tp.getOrientationRForTest();
                    tp.rebuildOrientationTransform(); // no layout change -> the cached instance must be reused
                    if ( tp.getOrientationRForTest() != r1 ) {
                        fail( ok, "orientation transform R was rebuilt without a layout change (cache ineffective)" );
                    }
                    tp.calcParametersForPainting( w, h ); // a layout pass must invalidate the cache
                    tp.rebuildOrientationTransform();
                    if ( tp.getOrientationRForTest() == r1 ) {
                        fail( ok, "orientation transform R was not rebuilt after a layout change (stale R)" );
                    }
                    final java.awt.geom.AffineTransform r2 = tp.getOrientationRForTest();
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_BOTTOM ); // an orientation change must also rebuild R
                    tp.rebuildOrientationTransform();
                    if ( tp.getOrientationRForTest() == r2 ) {
                        fail( ok, "orientation transform R was not rebuilt after an orientation change" );
                    }

                    // screen-culling in a vertical orientation (via the viewport mapped back into logical space): with
                    // the tree zoomed well past the REALIZED viewport, a node centered in view must NOT be culled and a
                    // node far outside MUST be -- i.e. the optimization is both correct and actually engaging.
                    tp.setTreeOrientation( TREE_ORIENTATION.ROOT_TOP );
                    frame.showWhole();
                    for ( int i = 0; i < 12; ++i ) {
                        cp.zoomInScreenX( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                        cp.zoomInScreenY( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                    }
                    frame.getMainPanel().getCurrentScrollPane().validate();
                    // a real paint assigns every node's coords at the zoomed scale (calcParametersForPainting sets only
                    // the scale factors, not per-node coords) and rebuilds R, so screenPointFor + the cull read true
                    // positions
                    final java.awt.image.BufferedImage scratch = new java.awt.image.BufferedImage( 60, 60,
                            java.awt.image.BufferedImage.TYPE_INT_RGB );
                    final java.awt.Graphics2D scg = scratch.createGraphics();
                    tp.paintComponent( scg );
                    scg.dispose();
                    final java.util.List<org.forester.phylogeny.PhylogenyNode> tips = phy.getExternalNodes();
                    final org.forester.phylogeny.PhylogenyNode centered = tips.get( tips.size() / 2 );
                    tp.centerOnNodeForTest( centered );
                    frame.getMainPanel().getCurrentScrollPane().validate();
                    if ( tp.isNodeDataInvisibleForTest( centered ) ) {
                        fail( ok, "vertical culling wrongly culled a node centered in the viewport" );
                    }
                    // the tip whose on-screen position is farthest outside the (scrolled) viewport must be culled
                    final java.awt.Rectangle vw = frame.getMainPanel().getCurrentScrollPane().getViewport()
                            .getViewRect();
                    org.forester.phylogeny.PhylogenyNode far = null;
                    double far_dist = 0;
                    for ( final org.forester.phylogeny.PhylogenyNode t : tips ) {
                        final java.awt.geom.Point2D.Double p = tp.screenPointFor( t );
                        final double dx = Math.max( 0, Math.max( vw.x - p.x, p.x - ( vw.x + vw.width ) ) );
                        final double dy = Math.max( 0, Math.max( vw.y - p.y, p.y - ( vw.y + vw.height ) ) );
                        final double d = Math.max( dx, dy );
                        if ( d > far_dist ) {
                            far_dist = d;
                            far = t;
                        }
                    }
                    if ( ( far == null ) || ( far_dist < 600 ) ) {
                        fail( ok, "test setup: no tip far enough outside the viewport after zoom (max dist " + far_dist
                                + ")" );
                    }
                    else if ( !tp.isNodeDataInvisibleForTest( far ) ) {
                        fail( ok, "vertical culling did not cull a node " + (int) far_dist
                                + "px outside the viewport (optimization not engaging)" );
                    }

                    // "F"/showWhole must FIT the tree to the window in a vertical orientation: the tip-label footprint
                    // is reserved on the axis it actually extends along (depth for the vertical/45deg labels), so the
                    // breadth does NOT overflow the viewport width. Long labels make the (old) breadth over-reservation
                    // overflow clearly. Checked for both vertical orientations.
                    for ( final org.forester.phylogeny.PhylogenyNode t : tips ) {
                        t.setName( t.getName() + "_sapiens_longlabel" );
                    }
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    for ( final TREE_ORIENTATION orient : new TREE_ORIENTATION[] { TREE_ORIENTATION.ROOT_TOP,
                            TREE_ORIENTATION.ROOT_BOTTOM } ) {
                        tp.setTreeOrientation( orient );
                        frame.showWhole();
                        frame.getMainPanel().getCurrentScrollPane().validate();
                        final java.awt.Dimension vpz = frame.getMainPanel().getCurrentScrollPane().getViewport()
                                .getExtentSize();
                        final java.awt.Dimension pref = tp.getPreferredSize();
                        if ( pref.width > ( vpz.width + 8 ) ) {
                            fail( ok, orient + ": showWhole overflows the viewport width (pref " + pref.width + " > vp "
                                    + vpz.width + ") -- breadth over-reserved" );
                        }
                        if ( pref.height < ( vpz.height * 0.7 ) ) {
                            fail( ok, orient + ": showWhole under-fills the viewport height (pref " + pref.height
                                    + " vs vp " + vpz.height + ") -- depth over-reserved" );
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

    /** The fraction of the scroll range at the viewport's center (0 = start, 1 = end); ~0.5 when centered. */
    private static double centerFraction( final JScrollBar sb ) {
        final int range = sb.getMaximum() - sb.getMinimum();
        return ( range <= 0 ) ? 0.5 : ( ( sb.getValue() + ( sb.getVisibleAmount() / 2.0 ) ) / range );
    }

    private static void layout( final MainFrame frame, final TreePanel tp, final int w, final int h ) {
        frame.showWhole();
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        AptxUtil.renderPhylogenyToImage( w, h, tp, frame.getOptions(), false, 1, false );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [OrientationZoomControlTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [OrientationZoomControlTest] " + msg );
        ok[ 0 ] = false;
    }

    private OrientationZoomControlTest() {
    }
}
