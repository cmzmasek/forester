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
        return barsRenderOk();
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
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    final int vertical_on = countSepia( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( vertical_on <= ( off + 300 ) ) {
                        fail( ok, "Fossil Range Bars should draw in a vertical orientation (on=" + vertical_on + " off="
                                + off + ")" );
                    }
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );

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
