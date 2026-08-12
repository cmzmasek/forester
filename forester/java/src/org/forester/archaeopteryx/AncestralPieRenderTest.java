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

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
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
            // add a numeric tip property so Size-by can be activated for the Size-by radial orphan-legend check
            int szi = 1;
            for ( final org.forester.phylogeny.PhylogenyNode leaf : phy.getExternalNodes() ) {
                org.forester.phylogeny.data.PropertiesList pl = leaf.getNodeData().getProperties();
                if ( pl == null ) {
                    pl = new org.forester.phylogeny.data.PropertiesList();
                    leaf.getNodeData().setProperties( pl );
                }
                pl.addProperty( new org.forester.phylogeny.data.Property( "data:sz", Integer.toString( szi++ ), "",
                        "xsd:decimal", org.forester.phylogeny.data.Property.AppliesTo.NODE ) );
            }
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
                    // a REAL multi-wedge pie must render: some small (8px) region holds >= 2 distinct state hues
                    // (adjacent wedges of one pie). The window is smaller than the tip/pie spacing and the legend's
                    // row pitch, so it is NOT triggered by two separate discs or two legend swatches -- only by wedges
                    // meeting inside one pie. Transform-independent (scans the image, not node coords).
                    if ( !hasMultiWedgeRegion( on_img ) ) {
                        fail( ok, "a real multi-wedge ancestral pie must render (>= 2 state hues within one pie)" );
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

                    // VERTICAL PARITY: pies ride the rotation in a root-top orientation (the disc stays a disc); assert
                    // the colored pie ink is still there (comparable to horizontal, above the legend-only floor).
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    if ( countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) ) < ( on / 3 ) ) {
                        fail( ok, "pies must draw in a vertical (root-top) orientation" );
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

                    // RADIAL PARITY: ancestral pies also render in circular AND unrooted layouts (node-centered discs),
                    // with their legend interactive there.
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    tp.calcParametersForPainting( w, h );
                    if ( countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) ) < ( on / 3 ) ) {
                        fail( ok, "pies must draw in a CIRCULAR layout" );
                    }
                    final Rectangle clb = tp.getAncestralPieLegendBounds();
                    if ( ( clb == null ) || !tp.isOnAncestralPieLegend(
                            mouseAt( tp, clb.x + ( clb.width / 2 ), clb.y + ( clb.height / 2 ) ) ) ) {
                        fail( ok, "the pie legend must be interactive in a circular layout" );
                    }
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                    tp.calcParametersForPainting( w, h );
                    if ( countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) ) < ( on / 3 ) ) {
                        fail( ok, "pies must draw in an UNROOTED layout" );
                    }

                    // ORPHAN-LEGEND FIX: the Color-by tip-dot legend draws in rectangular but is SUPPRESSED (its bounds
                    // nulled) in a radial layout, where the tip dots it keys are not drawn. Uses the SCREEN paint path
                    // (printAll), which is where legend bounds are recorded.
                    tp.setColorByPropertyRef( "beast:location" );
                    if ( tp.isColorByProperty() ) {
                        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                        tp.calcParametersForPainting( w, h );
                        paintScreen( tp, w, h );
                        if ( tp.getPropertyLegendBounds() == null ) {
                            fail( ok, "the Color-by legend must draw (record bounds) in a rectangular layout" );
                        }
                        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                        tp.calcParametersForPainting( w, h );
                        paintScreen( tp, w, h );
                        if ( tp.getPropertyLegendBounds() != null ) {
                            fail( ok, "the Color-by legend must be suppressed (bounds nulled) in a radial layout" );
                        }
                    }
                    tp.setColorByPropertyRef( null );
                    // ORPHAN-LEGEND FIX (Size-by): symmetric to Color-by -- the size-dot legend is suppressed (bounds
                    // nulled) in a radial layout too, where the size dots it keys are not drawn.
                    tp.setSizeByPropertyRef( "data:sz" );
                    if ( tp.isSizeByProperty() ) {
                        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                        tp.calcParametersForPainting( w, h );
                        paintScreen( tp, w, h );
                        if ( tp.getSizeLegendBounds() == null ) {
                            fail( ok, "the Size-by legend must draw (record bounds) in a rectangular layout" );
                        }
                        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                        tp.calcParametersForPainting( w, h );
                        paintScreen( tp, w, h );
                        if ( tp.getSizeLegendBounds() != null ) {
                            fail( ok, "the Size-by legend must be suppressed (bounds nulled) in a radial layout" );
                        }
                    }
                    tp.setSizeByPropertyRef( null );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.calcParametersForPainting( w, h );

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

    /** Drives the SCREEN paint path (printAll -> paintPhylogeny with to_screen=true), which is where the legends
     *  record their draggable bounds -- unlike renderPhylogenyToImage, which is an export (bounds not recorded). */
    private static void paintScreen( final TreePanel tp, final int w, final int h ) {
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.printAll( g );
        g.dispose();
    }

    private static MouseEvent mouseAt( final TreePanel tp, final int x, final int y ) {
        return mouseAtClicks( tp, x, y, 1 );
    }

    private static MouseEvent mouseAtClicks( final TreePanel tp, final int x, final int y, final int clicks ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(), 0, x, y, clicks, false );
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

    /** Whether some small (WIN x WIN) region of the image holds &ge;2 distinct dominant state hues -- a real
     *  multi-wedge pie (adjacent wedges of ONE pie). WIN is smaller than the tip/pie spacing and the legend's row
     *  pitch, so two separate discs or two stacked legend swatches never both land in one window; only wedges meeting
     *  inside a pie do. Precomputes a per-pixel hue bin once, then window-scans -- transform/scale-independent. */
    private static boolean hasMultiWedgeRegion( final BufferedImage img ) {
        final int w = img.getWidth(), h = img.getHeight();
        final int[] hue = new int[ w * h ];
        for( int y = 0; y < h; ++y ) {
            for( int x = 0; x < w; ++x ) {
                final int rgb = img.getRGB( x, y );
                if ( !isVivid( rgb ) ) {
                    hue[ ( y * w ) + x ] = -1;
                }
                else {
                    final float[] hsb = java.awt.Color.RGBtoHSB( ( rgb >> 16 ) & 0xFF, ( rgb >> 8 ) & 0xFF,
                            rgb & 0xFF, null );
                    hue[ ( y * w ) + x ] = Math.min( HUE_BINS - 1, (int) ( hsb[ 0 ] * HUE_BINS ) );
                }
            }
        }
        final int win = 8, step = 3, min_px = 4; // small window + low per-hue floor to catch a thin second wedge
        final int[] counts = new int[ HUE_BINS ];
        final int[] dom = new int[ HUE_BINS ];
        for( int y0 = 0; ( y0 + win ) <= h; y0 += step ) {
            for( int x0 = 0; ( x0 + win ) <= w; x0 += step ) {
                java.util.Arrays.fill( counts, 0 );
                for( int y = y0; y < ( y0 + win ); ++y ) {
                    for( int x = x0; x < ( x0 + win ); ++x ) {
                        final int hb = hue[ ( y * w ) + x ];
                        if ( hb >= 0 ) {
                            counts[ hb ]++;
                        }
                    }
                }
                int nd = 0;
                for( int k = 0; k < HUE_BINS; ++k ) {
                    if ( counts[ k ] >= min_px ) {
                        dom[ nd++ ] = k;
                    }
                }
                // a real second WEDGE is a DIFFERENT state color, which sits >= 2 hue bins away (the state palette is
                // evenly spread around the wheel). Two ADJACENT dominant bins are just ONE color's antialiasing
                // straddling a bin boundary -- not a second wedge -- so require a circular bin distance >= 2.
                for( int p = 0; p < nd; ++p ) {
                    for( int q = p + 1; q < nd; ++q ) {
                        final int diff = Math.abs( dom[ p ] - dom[ q ] );
                        if ( Math.min( diff, HUE_BINS - diff ) >= 2 ) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
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
