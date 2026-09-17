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
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the fossil-range demo phylogram (forester/demo/fossil-range-bars.xml) with "Fossil Range Bars (FAD/LAD)" ON
 * and asserts sepia range bars appear on the tips (and none when off), that loading the demo auto-enables the toggle,
 * and that the bars draw in the vertical (root-top, riding R) AND circular layouts -- display-type parity. Headful; a
 * green no-op when headless. Dogfoods the demo.
 */
public final class FossilRangeBarRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "FossilRangeBarRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return barsRenderOk() && noiseIntervalsDrawNothingOk();
    }

    /**
     * TreeAnnotator writes an exactly dated tip's height as a subtraction, so its "range" comes out one float ULP
     * wide ({@code height_95%_HPD={9.0,9.000000000000004}}); influenza.tree carries 686 of them. Before 2026-09-17
     * that switched the fossil-range overlay on and drew a bracketed 1 px range under every tip of a virus tree.
     * Neither the auto-enable nor the painter may treat it as a range -- and the same fixture with a REAL range must
     * still draw, or this would pass on a tree that simply cannot reach the painter.
     */
    private static boolean noiseIntervalsDrawNothingOk() {
        try {
            final Phylogeny phy = noiseIntervalTree();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "noise" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    if ( o.isShowFossilRangeBars() ) {
                        fail( ok, "float-noise tip intervals must NOT auto-enable the Fossil Range Bars" );
                    }
                    if ( o.isShowHpdBars() ) {
                        fail( ok, "float-noise intervals must not auto-enable the Node Age Bars either" );
                    }
                    o.setGraphicsExportWhiteBackground( true );
                    tp.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final int w = 900, h = 460;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    // forced ON, BOTH overlays: even then, noise must draw nothing -- in either layout
                    o.setShowFossilRangeBars( true );
                    o.setShowHpdBars( true );
                    final BufferedImage noise_img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    if ( ( countSepia( noise_img ) > 0 ) || ( countHpdBlue( noise_img ) > 0 ) ) {
                        fail( ok, "rectangular: float noise must draw no bars, got " + countSepia( noise_img )
                                + " sepia + " + countHpdBlue( noise_img ) + " blue px" );
                    }
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    o.setShowOverview( false );
                    tp.setOvOn( false );
                    frame.showWhole();
                    tp.setPreferredSize( new java.awt.Dimension( w, h ) );
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage circ_noise = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    if ( ( countSepia( circ_noise ) > 0 ) || ( countHpdBlue( circ_noise ) > 0 ) ) {
                        fail( ok, "circular: float noise must draw no bars, got " + countSepia( circ_noise )
                                + " sepia + " + countHpdBlue( circ_noise ) + " blue px" );
                    }
                    // REACHABILITY: give one tip a real range and one internal node a real HPD -- the same fixture,
                    // in both layouts, must now draw both. Without this the checks above would also pass on a tree
                    // that simply cannot reach the painters.
                    final PhylogenyNode tip = phy.getFirstExternalNode();
                    tip.getNodeData().getDate().setMin( new java.math.BigDecimal( "4.0" ) );
                    tip.getNodeData().getDate().setMax( new java.math.BigDecimal( "7.0" ) );
                    final PhylogenyNode inner = phy.getRoot().getChildNode( 0 );
                    inner.getNodeData().getDate().setMin( new java.math.BigDecimal( "5.0" ) );
                    inner.getNodeData().getDate().setMax( new java.math.BigDecimal( "8.0" ) );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage circ_real = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    if ( ( countSepia( circ_real ) < 20 ) || ( countHpdBlue( circ_real ) < 20 ) ) {
                        fail( ok, "circular: a REAL range and HPD must draw, got " + countSepia( circ_real )
                                + " sepia + " + countHpdBlue( circ_real ) + " blue px" );
                    }
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage real_img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    if ( ( countSepia( real_img ) < 20 ) || ( countHpdBlue( real_img ) < 20 ) ) {
                        fail( ok, "rectangular: a REAL range and HPD must draw, got " + countSepia( real_img )
                                + " sepia + " + countHpdBlue( real_img ) + " blue px" );
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

    /** Four dated tips with branch lengths, unit-less ages (so the tree is not on calendar time), every tip carrying
     *  a height and its own float noise as the "interval" -- a BEAST MCC tree that did not convert. */
    private static Phylogeny noiseIntervalTree() {
        final PhylogenyNode root = new PhylogenyNode();
        noiseDate( root, "10" );
        final PhylogenyNode left = new PhylogenyNode();
        noiseDate( left, "6" );
        final PhylogenyNode right = new PhylogenyNode();
        noiseDate( right, "4" );
        left.setDistanceToParent( 4 );
        right.setDistanceToParent( 6 );
        root.addAsChild( left );
        root.addAsChild( right );
        final String[] names = { "tip_a", "tip_b", "tip_c", "tip_d" };
        final String[] heights = { "5", "0", "3", "0" };
        for( int i = 0; i < names.length; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( names[ i ] );
            noiseDate( tip, heights[ i ] );
            final PhylogenyNode parent = ( i < 2 ) ? left : right;
            tip.setDistanceToParent( Double.parseDouble( parent.getNodeData().getDate().getValue().toPlainString() )
                    - Double.parseDouble( heights[ i ] ) );
            parent.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /** A height with TreeAnnotator's own float noise as its 95% HPD: {h, h + 1e-14}. */
    private static void noiseDate( final PhylogenyNode n, final String height ) {
        final java.math.BigDecimal h = new java.math.BigDecimal( height );
        n.getNodeData().setDate( new org.forester.phylogeny.data.Date( "", h, h,
                                                                       h.add( new java.math.BigDecimal( "0.00000000000001" ) ),
                                                                       "" ) );
    }

    private static boolean barsRenderOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/fossil-range-bars.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "fossil" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    // fossil-range-bars.xml has EXTERNAL-tip date intervals (FAD/LAD), so loading it must AUTO-ENABLE
                    // "Fossil Range Bars" -- both the Option and its kept-in-sync menu item -- so the ranges show
                    // without the user hunting for the toggle. Assert this before the test overrides the option below.
                    if ( !o.isShowFossilRangeBars() ) {
                        fail( ok, "Fossil Range Bars must auto-enable when the loaded tree has tip date intervals" );
                    }
                    if ( ( frame._show_fossil_range_bars_cbmi != null )
                            && !frame._show_fossil_range_bars_cbmi.isSelected() ) {
                        fail( ok, "the Fossil-Range-Bars menu item must be checked after the auto-enable" );
                    }
                    o.setGraphicsExportWhiteBackground( true ); // predictable white background for the sepia composite
                    // the demo auto-derives a GEOLOGIC axis (unit "mya"); its coloured ICS bands include sepia/brown
                    // tones, so turn the time axis OFF here to isolate the fossil bars as the only sepia content
                    tp.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final int w = 900, h = 460;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    o.setShowFossilRangeBars( false );
                    final int off = countSepia( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    o.setShowFossilRangeBars( true );
                    final BufferedImage on_img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int on = countSepia( on_img );
                    if ( on <= ( off + 300 ) ) {
                        fail( ok, "Fossil Range Bars should add many sepia pixels (on=" + on + " off=" + off + ")" );
                    }
                    if ( off >= 100 ) {
                        fail( ok, "no sepia bars should appear when Fossil Range Bars is off, got " + off );
                    }
                    // the bar must sit at a fossil TIP's row: pick an external tip with a FAD/LAD interval and assert
                    // sepia pixels appear across its y (the bar runs from FAD to LAD along that row)
                    PhylogenyNode tip = null;
                    for ( final java.util.Iterator<PhylogenyNode> it = phy.iteratorExternalForward(); it.hasNext(); ) {
                        final PhylogenyNode n = it.next();
                        if ( n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getMin() != null )
                                && ( n.getNodeData().getDate().getMax() != null ) ) {
                            tip = n;
                            break;
                        }
                    }
                    if ( tip == null ) {
                        fail( ok, "expected a fossil tip with a date interval in the demo" );
                    }
                    else if ( !sepiaInRow( on_img, Math.round( tip.getYcoord() ) ) ) {
                        fail( ok, "a fossil range bar must appear at the tip's row y=" + Math.round( tip.getYcoord() ) );
                    }
                    // VERTICAL PARITY: fossil bars are plain rects, so they ride the rotation R into vertical bars in a
                    // root-top orientation. Render in ROOT_TOP and confirm the sepia draws.
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    final int vertical_on = countSepia( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( vertical_on <= ( off + 300 ) ) {
                        fail( ok, "Fossil Range Bars should draw in a vertical orientation (on=" + vertical_on + " off="
                                + off + ")" );
                    }
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );

                    // CIRCULAR PARITY: on a circular phylogram the fossil ranges draw as radial sepia segments on the
                    // tips (paintFossilRangeBarsCircular). Confirm on >> off in the circular layout.
                    final int cw = 820, ch = 820;
                    o.setShowOverview( false );
                    tp.setOvOn( false );
                    frame.showWhole();
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    tp.setPreferredSize( new java.awt.Dimension( cw, ch ) );
                    tp.setSize( cw, ch );
                    o.setShowFossilRangeBars( false );
                    tp.calcParametersForPainting( cw, ch );
                    final int circ_off = countSepia( AptxUtil.renderPhylogenyToImage( cw, ch, tp, o, false, 1, false ) );
                    o.setShowFossilRangeBars( true );
                    tp.calcParametersForPainting( cw, ch );
                    final int circ_on = countSepia( AptxUtil.renderPhylogenyToImage( cw, ch, tp, o, false, 1, false ) );
                    if ( circ_on <= ( circ_off + 150 ) ) {
                        fail( ok, "circular Fossil Range Bars must add sepia radial segments (on=" + circ_on + " off="
                                + circ_off + ")" );
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

    /** Pixels of the sepia fossil-range bar (FOSSIL_BAR_COLOR ~= 150,100,55 composited over white or the black
     *  branches): red clearly dominates green dominates blue -- excludes white/black/gray (r~=g~=b) and any bluish
     *  or pure-red content. */
    private static int countSepia( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( isSepia( img.getRGB( x, y ) ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Any sepia pixel anywhere in the row band [y-3, y+3]? (the bar is FOSSIL_BAR_HEIGHT=5 tall, centred on y). */
    private static boolean sepiaInRow( final BufferedImage img, final int y ) {
        for( int yy = Math.max( 0, y - 3 ); yy < Math.min( img.getHeight(), y + 4 ); ++yy ) {
            for( int xx = 0; xx < img.getWidth(); ++xx ) {
                if ( isSepia( img.getRGB( xx, yy ) ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    /** Pixels of the translucent blue node-age (HPD) bar (HPD_BAR_COLOR = 70,130,220 at alpha 90, composited over
     *  white): blue clearly dominates, which no branch, label or sepia pixel does. */
    private static int countHpdBlue( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( b >= ( r + 25 ) ) && ( b >= ( g + 10 ) ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean isSepia( final int rgb ) {
        final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
        return ( r >= ( g + 20 ) ) && ( g >= ( b + 20 ) ) && ( r >= ( b + 55 ) );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [FossilRangeBarRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [FossilRangeBarRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private FossilRangeBarRenderTest() {
    }
}
