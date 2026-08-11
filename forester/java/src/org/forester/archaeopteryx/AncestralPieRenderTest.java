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
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Iterator;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the ancestral-pie demo (forester/demo/ancestral-pie-charts.xml) with "Ancestral pie: location" selected and
 * asserts: an internal (multi-state) node shows &ge;2 distinct wedge colors (a real pie); a tip shows one solid state
 * color; the pies appear in a vertical (root-top) orientation too; the state legend renders; and selecting "None"
 * removes them. Headful; a green no-op when headless. Dogfoods the demo.
 */
public final class AncestralPieRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "AncestralPieRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return piesRenderOk();
    }

    private static boolean piesRenderOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/ancestral-pie-charts.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "pie" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final ControlPanel cp = tp.getControlPanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true ); // predictable white background for vivid-color detection
                    final int w = 900, h = 520;
                    frame.showWhole();
                    tp.setSize( w, h );

                    // the dropdown appears only for a tree with discrete traits, and offers "location"
                    if ( !cp.isAncestralPieControlVisible() ) {
                        fail( ok, "the 'Ancestral pie' control must be visible for a tree with discrete traits" );
                    }
                    if ( !cp.ancestralPieTraitRefs().contains( "location" ) ) {
                        fail( ok, "the 'Ancestral pie' dropdown must offer 'location', got "
                                + cp.ancestralPieTraitRefs() );
                    }

                    // OFF baseline: no trait selected -> no vivid pie ink
                    tp.setAncestralPieTrait( null );
                    tp.calcParametersForPainting( w, h );
                    final int off = countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );

                    // ON: select the trait
                    tp.setAncestralPieTrait( "location" );
                    cp.populateAncestralPieBox();
                    if ( !tp.isShowAncestralPies() ) {
                        fail( ok, "ancestral pies should be active after selecting the trait" );
                    }
                    if ( !"location".equals( cp.getAncestralPieSelection() ) ) {
                        fail( ok, "the dropdown should show 'location', got " + cp.getAncestralPieSelection() );
                    }
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage on_img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int on = countVivid( on_img );
                    if ( on <= ( off + 500 ) ) {
                        fail( ok, "selecting a trait must add many colored pie pixels (on=" + on + " off=" + off + ")" );
                    }
                    if ( off >= 100 ) {
                        fail( ok, "no colored pies should appear when no trait is selected, got " + off );
                    }

                    // an internal (multi-state) node must show >= 2 DISTINCT wedge colors -> a real pie
                    final int r = (int) Math.ceil( tp.ancestralPieDiameterForTest() / 2.0 ) + 1;
                    final PhylogenyNode inode = firstMultiStateInternalNode( phy );
                    if ( inode == null ) {
                        fail( ok, "expected a multi-state internal node in the demo" );
                    }
                    else {
                        final int colors = distinctVividColors( on_img, Math.round( inode.getXcoord() ),
                                Math.round( inode.getYcoord() ), r );
                        if ( colors < 2 ) {
                            fail( ok, "an internal node's pie must show >= 2 distinct wedge colors, got " + colors );
                        }
                    }

                    // a tip must show exactly ONE solid state color (a single-state disc)
                    final PhylogenyNode tip = phy.getFirstExternalNode();
                    final int tip_colors = distinctVividColors( on_img, Math.round( tip.getXcoord() ),
                            Math.round( tip.getYcoord() ), r );
                    if ( tip_colors != 1 ) {
                        fail( ok, "a tip's pie must be a single solid color, got " + tip_colors );
                    }

                    // STATE->COLOR STABILITY: the shared map assigns each of the 4 states a distinct, stable color
                    // (a regression to per-node coloring would break the legend's meaning)
                    final String[] all_states = { "Africa", "Americas", "Asia", "Europe" };
                    final java.util.Set<java.awt.Color> state_colors = new java.util.HashSet<>();
                    for ( final String st : all_states ) {
                        final java.awt.Color c = tp.ancestralPieColor( st );
                        if ( c == null ) {
                            fail( ok, "every state must have a stable pie color; '" + st + "' had none" );
                        }
                        else {
                            state_colors.add( c );
                        }
                    }
                    if ( state_colors.size() != all_states.length ) {
                        fail( ok, "the " + all_states.length + " states must get DISTINCT colors, got "
                                + state_colors.size() );
                    }

                    // VERTICAL PARITY: pies ride the rotation in a root-top orientation (the disc stays a disc). Sample
                    // the internal node's ROTATED device position (screenPointFor) and require a real multi-wedge pie
                    // there -- so this asserts parity, not merely "some vivid ink somewhere".
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    final BufferedImage vimg = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    if ( countVivid( vimg ) < ( on / 3 ) ) {
                        fail( ok, "pies must draw in a vertical (root-top) orientation" );
                    }
                    if ( inode != null ) {
                        final java.awt.geom.Point2D.Double vp = tp.screenPointFor( inode );
                        final int vcolors = distinctVividColors( vimg, (int) Math.round( vp.x ),
                                (int) Math.round( vp.y ), r );
                        if ( vcolors < 2 ) {
                            fail( ok, "a vertical-orientation internal pie must show >= 2 wedge colors at its rotated "
                                    + "position, got " + vcolors );
                        }
                    }
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );

                    // the state legend renders (>= 2 distinct swatch colors for a multi-state trait); and in B&W it
                    // grays the swatches (0 vivid colors) so the key matches the gray wedges
                    final BufferedImage legend = new BufferedImage( 320, 260, BufferedImage.TYPE_INT_ARGB );
                    Graphics2D lg = legend.createGraphics();
                    tp.drawAncestralPieLegendForTest( lg, new Rectangle( 0, 0, 320, 260 ), true, false );
                    lg.dispose();
                    if ( distinctVividColors( legend, 160, 130, 200 ) < 2 ) {
                        fail( ok, "the ancestral-pie legend must render >= 2 distinct state swatch colors" );
                    }
                    final BufferedImage bwlegend = new BufferedImage( 320, 260, BufferedImage.TYPE_INT_ARGB );
                    lg = bwlegend.createGraphics();
                    tp.drawAncestralPieLegendForTest( lg, new Rectangle( 0, 0, 320, 260 ), false, true );
                    lg.dispose();
                    if ( countVivid( bwlegend ) != 0 ) {
                        fail( ok, "a black-and-white pie legend must gray its swatches (0 vivid colors), got "
                                + countVivid( bwlegend ) );
                    }

                    // LEGEND INTERACTION: hit-test, drag, and double-click-to-reset the pie legend
                    final BufferedImage limg = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
                    drawPieLegend( tp, limg, w, h );
                    final Rectangle lb = tp.getAncestralPieLegendBounds();
                    if ( lb == null ) {
                        fail( ok, "pie legend bounds must be recorded when drawn draggable" );
                    }
                    else {
                        final int cx = lb.x + ( lb.width / 2 ), cy = lb.y + ( lb.height / 2 );
                        if ( !tp.isOnAncestralPieLegend( mouseAt( tp, cx, cy ) ) ) {
                            fail( ok, "isOnAncestralPieLegend must be true inside the box" );
                        }
                        if ( tp.isOnAncestralPieLegend( mouseAt( tp, lb.x - 40, cy ) ) ) {
                            fail( ok, "isOnAncestralPieLegend must be false outside the box" );
                        }
                        if ( !tp.isOnAnyLegend( mouseAt( tp, cx, cy ) ) ) {
                            fail( ok, "isOnAnyLegend must be true over the pie legend" );
                        }
                        // drag it to the center, then double-click to reset it back to its default corner
                        tp.startLegendDrag( mouseAt( tp, cx, cy ) );
                        tp.dragLegend( mouseAt( tp, w / 2, h / 2 ) );
                        drawPieLegend( tp, limg, w, h );
                        final Rectangle moved = tp.getAncestralPieLegendBounds();
                        if ( moved.equals( lb ) ) {
                            fail( ok, "dragging the pie legend must move it" );
                        }
                        new MouseListener( tp )
                                .mouseClicked( mouseAtClicks( tp, moved.x + ( moved.width / 2 ),
                                        moved.y + ( moved.height / 2 ), 2 ) );
                        drawPieLegend( tp, limg, w, h );
                        if ( !tp.getAncestralPieLegendBounds().equals( lb ) ) {
                            fail( ok, "double-clicking the pie legend must RESET it to its default corner" );
                        }
                    }

                    // ORPHAN-LEGEND GUARD: in a circular/unrooted layout the pies are not drawn, so neither the legend
                    // nor its hit region may appear (the layout gate)
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    final int circular = countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( circular >= 100 ) {
                        fail( ok, "no pies OR pie legend may draw in a circular layout, got " + circular
                                + " colored pixels" );
                    }
                    if ( tp.getAncestralPieLegendBounds() != null
                            && tp.isOnAncestralPieLegend( mouseAt( tp, tp.getAncestralPieLegendBounds().x + 2,
                                    tp.getAncestralPieLegendBounds().y + 2 ) ) ) {
                        fail( ok, "the pie legend must not be interactive in a circular layout (layout gate)" );
                    }
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );

                    // selecting "None" removes the pies (mutation guard for the whole feature)
                    tp.setAncestralPieTrait( null );
                    tp.calcParametersForPainting( w, h );
                    final int none = countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( none >= 100 ) {
                        fail( ok, "selecting 'None' must remove the pies, got " + none + " colored pixels" );
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

    /** Draws the pie legend (draggable, so its bounds are recorded) into {@code img} at (0,0,w,h). */
    private static void drawPieLegend( final TreePanel tp, final BufferedImage img, final int w, final int h ) {
        final Graphics2D g = img.createGraphics();
        tp.drawAncestralPieLegendForTest( g, new Rectangle( 0, 0, w, h ), true, false );
        g.dispose();
    }

    private static MouseEvent mouseAt( final TreePanel tp, final int x, final int y ) {
        return mouseAtClicks( tp, x, y, 1 );
    }

    private static MouseEvent mouseAtClicks( final TreePanel tp, final int x, final int y, final int clicks ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(), 0, x, y, clicks, false );
    }

    private static PhylogenyNode firstMultiStateInternalNode( final Phylogeny phy ) {
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && ( TreePanelUtil.stateDistribution( n, "location" ).size() >= 2 ) ) {
                return n;
            }
        }
        return null;
    }

    private static final int HUE_BINS = 12;   // 30-degree hue bins -- enough to separate the 4 evenly-spread states
    private static final int MIN_HUE_PIXELS = 8; // a "dominant" wedge/disc, so AA edge-blend stragglers don't count

    /** A "vivid" pixel: clearly chromatic (not white/black/gray, near-white AA, or antialiased text) -- pie ink.
     *  Gated on SATURATION (not raw RGB chroma) so the lightening ramp at a disc's antialiased edge -- same hue,
     *  falling saturation -- is excluded rather than split into many buckets. */
    private static boolean isVivid( final int rgb ) {
        if ( ( ( rgb >>> 24 ) & 0xFF ) < 128 ) {
            return false; // transparent (the legend image is ARGB)
        }
        final float[] hsb = java.awt.Color.RGBtoHSB( ( rgb >> 16 ) & 0xFF, ( rgb >> 8 ) & 0xFF, rgb & 0xFF, null );
        return ( hsb[ 1 ] >= 0.35f ) && ( hsb[ 2 ] >= 0.30f );
    }

    private static int countVivid( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( isVivid( img.getRGB( x, y ) ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Number of DOMINANT hues (hue bins with &ge; {@link #MIN_HUE_PIXELS} vivid pixels) in a box of half-size
     *  {@code rad} centred on {@code (cx,cy)} -- so a solid 1-state disc yields 1 (its AA lightening ramp shares
     *  the one hue) and a multi-wedge pie yields &ge;2 (one bin per substantial wedge). */
    private static int distinctVividColors( final BufferedImage img, final int cx, final int cy, final int rad ) {
        final int[] counts = new int[ HUE_BINS ];
        for( int y = Math.max( 0, cy - rad ); y < Math.min( img.getHeight(), cy + rad + 1 ); ++y ) {
            for( int x = Math.max( 0, cx - rad ); x < Math.min( img.getWidth(), cx + rad + 1 ); ++x ) {
                final int rgb = img.getRGB( x, y );
                if ( isVivid( rgb ) ) {
                    final float[] hsb = java.awt.Color.RGBtoHSB( ( rgb >> 16 ) & 0xFF, ( rgb >> 8 ) & 0xFF,
                            rgb & 0xFF, null );
                    counts[ Math.min( HUE_BINS - 1, (int) ( hsb[ 0 ] * HUE_BINS ) ) ]++;
                }
            }
        }
        int dominant = 0;
        for( final int c : counts ) {
            if ( c >= MIN_HUE_PIXELS ) {
                ++dominant;
            }
        }
        return dominant;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [AncestralPieRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [AncestralPieRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private AncestralPieRenderTest() {
    }
}
