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
        return barsRenderOk() && calendarDirectionOk() && sampledTipBarsOk();
    }

    /** On CALENDAR time a tip is a dated SAMPLE. A tip whose sampling date is known only to the year carries a genuine
     *  date interval, and the Node Age Bars draw it (in the calendar direction); a tip dated exactly ({d,d}) draws
     *  nothing, a tip with no interval draws nothing -- and the geologic Fossil Range Bars draw NOTHING on such a tree
     *  even when their global toggle is on (another open tab may hold fossils). Both Nextstrain readers used to throw
     *  tip intervals away just to keep those bars off. */
    private static boolean sampledTipBarsOk() {
        final boolean[] ok = { true };
        try {
            // YEAR: dated only to 2020 -> value 2020.9 inside [2020.0, 2020.999]; EXACT: {2021.0, 2021.0}; PLAIN: no interval
            final String json = "{\"version\":\"v2\",\"tree\":{\"name\":\"R\",\"node_attrs\":{\"num_date\":{\"value\":2019.0}},"
                    + "\"children\":[{\"name\":\"PLAIN\",\"node_attrs\":{\"num_date\":{\"value\":2019.5}}},"
                    + "{\"name\":\"MID\",\"node_attrs\":{\"num_date\":{\"value\":2019.6}},"
                    + "\"children\":[{\"name\":\"YEAR\",\"node_attrs\":{\"num_date\":{\"value\":2020.9,\"confidence\":[2020.0,2020.999]}}},"
                    + "{\"name\":\"EXACT\",\"node_attrs\":{\"num_date\":{\"value\":2021.0,\"confidence\":[2021.0,2021.0]}}}]}]}}";
            final org.forester.io.parsers.json.AuspiceJsonParser parser = new org.forester.io.parsers.json.AuspiceJsonParser();
            parser.setSource( new StringBuffer( json ) );
            final Phylogeny phy = parser.parse()[ 0 ];
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(), "tips" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    // the load-time auto-enable: node-age bars ON for the uncertain tip, fossil bars OFF
                    if ( !o.isShowHpdBars() ) {
                        fail( ok, "a calendar tree with an uncertain sampling date must auto-enable the Node Age Bars" );
                    }
                    if ( o.isShowFossilRangeBars() ) {
                        fail( ok, "a calendar tree's tip intervals must NOT auto-enable the Fossil Range Bars" );
                    }
                    o.setGraphicsExportWhiteBackground( true );
                    o.setShowHpdBars( true );
                    o.setShowFossilRangeBars( true ); // as if another tab held fossils: still nothing to draw here
                    o.setNodeAgeShape( Options.NODE_AGE_SHAPE.BAR );
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final int w = 900, h = 460;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final org.forester.phylogeny.PhylogenyNode year = phy.getNode( "YEAR" );
                    final org.forester.phylogeny.PhylogenyNode exact = phy.getNode( "EXACT" );
                    final org.forester.phylogeny.PhylogenyNode plain = phy.getNode( "PLAIN" );
                    // YEAR: value 2020.9 in [2020.0, 2020.999] -> the bar reaches far to the EARLIER side (left)
                    final int yx = Math.round( year.getXcoord() ), yy = Math.round( year.getYcoord() );
                    int left_ext = 0;
                    for ( int dx = 8; dx <= 400; ++dx ) { // from 8 px out: clear of the node's own mark
                        if ( bluishNear( img, yx - dx, yx - dx, yy ) ) {
                            left_ext = dx;
                        }
                    }
                    if ( left_ext < 40 ) {
                        fail( ok, "a tip dated only to the year must draw its sampling-date uncertainty (bar reach " + left_ext + " px)" );
                    }
                    for ( final org.forester.phylogeny.PhylogenyNode n : new org.forester.phylogeny.PhylogenyNode[] { exact, plain } ) {
                        final int nx = Math.round( n.getXcoord() ), ny = Math.round( n.getYcoord() );
                        if ( bluishNear( img, nx - 60, nx - 8, ny ) ) {
                            fail( ok, "tip " + n.getName() + " has no date uncertainty and must draw no age bar" );
                        }
                    }
                    // colour-independent (this tree auto-colours by date, and that palette has sepia-like oranges): with
                    // the fossil toggle OFF the figure must be the very same, pixel for pixel
                    o.setShowFossilRangeBars( false );
                    final BufferedImage without_fossil_toggle = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    int differing = 0;
                    for ( int y = 0; y < img.getHeight(); ++y ) {
                        for ( int x = 0; x < img.getWidth(); ++x ) {
                            if ( img.getRGB( x, y ) != without_fossil_toggle.getRGB( x, y ) ) {
                                ++differing;
                            }
                        }
                    }
                    if ( differing > 0 ) {
                        fail( ok, "the Fossil Range Bars must draw nothing on a calendar tree, even with their toggle on ("
                                + differing + " px differ)" );
                    }
                    // CIRCULAR parity. The only interval in this tree is the uncertain tip's, so any pixel the Node Age Bars
                    // toggle changes IS that tip's bar; and the fossil toggle must change nothing here either.
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    frame.showWhole();
                    tp.setSize( 700, 700 );
                    tp.calcParametersForPainting( 700, 700 );
                    o.setShowFossilRangeBars( false );
                    o.setShowHpdBars( true );
                    final BufferedImage circ_on = AptxUtil.renderPhylogenyToImage( 700, 700, tp, o, false, 1, false );
                    o.setShowHpdBars( false );
                    final BufferedImage circ_off = AptxUtil.renderPhylogenyToImage( 700, 700, tp, o, false, 1, false );
                    o.setShowFossilRangeBars( true );
                    final BufferedImage circ_fossil = AptxUtil.renderPhylogenyToImage( 700, 700, tp, o, false, 1, false );
                    int bar_px = 0, fossil_px = 0;
                    for ( int y = 0; y < circ_on.getHeight(); ++y ) {
                        for ( int x = 0; x < circ_on.getWidth(); ++x ) {
                            bar_px += ( circ_on.getRGB( x, y ) != circ_off.getRGB( x, y ) ) ? 1 : 0;
                            fossil_px += ( circ_fossil.getRGB( x, y ) != circ_off.getRGB( x, y ) ) ? 1 : 0;
                        }
                    }
                    if ( bar_px < 20 ) {
                        fail( ok, "circular: the uncertain tip must draw its age bar too (" + bar_px + " px)" );
                    }
                    if ( fossil_px > 0 ) {
                        fail( ok, "circular: the Fossil Range Bars must draw nothing on a calendar tree (" + fossil_px + " px)" );
                    }
                    o.setShowHpdBars( true );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    // the Time Axis the user SHOWS is a display choice: switching it Off or to Geologic must not turn the
                    // sampling-date uncertainty into a fossil range (the fossil bars would then draw it)
                    for ( final Options.TIME_AXIS_TYPE shown : new Options.TIME_AXIS_TYPE[] { Options.TIME_AXIS_TYPE.NONE,
                            Options.TIME_AXIS_TYPE.GEOLOGIC } ) {
                        tp.setTimeAxisType( shown );
                        if ( !tp.isSampledTipWithDateUncertainty( year ) || !tp.isOnCalendarTime() ) {
                            fail( ok, "showing the " + shown + " axis must not change what the tip interval means" );
                        }
                        final BufferedImage fossil_on = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                        o.setShowFossilRangeBars( false );
                        final BufferedImage fossil_off = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                        o.setShowFossilRangeBars( true );
                        int diff = 0;
                        for ( int y = 0; y < fossil_on.getHeight(); ++y ) {
                            for ( int x = 0; x < fossil_on.getWidth(); ++x ) {
                                diff += ( fossil_on.getRGB( x, y ) != fossil_off.getRGB( x, y ) ) ? 1 : 0;
                            }
                        }
                        if ( diff > 0 ) {
                            fail( ok, "with the " + shown + " axis shown, the Fossil Range Bars must still draw nothing (" + diff + " px)" );
                        }
                    }
                    tp.setTimeAxisType( null );
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

    /** On a CALENDAR tree (dates increase toward the tips, opposite of geologic age), the node-age bar for an internal
     *  node with a right-skewed date confidence must extend FARTHER toward the later date (right) than the earlier date
     *  (left) -- i.e. the signed-corr flip is applied. A geologic (age) mapping would mirror it. */
    private static boolean calendarDirectionOk() {
        final boolean[] ok = { true };
        try {
            // MID (internal, mid-tree) has value 2020.0 with confidence [2019.9, 2020.5]: strongly skewed to LATER
            final String json = "{\"version\":\"v2\",\"tree\":{\"name\":\"R\",\"node_attrs\":{\"num_date\":{\"value\":2019.0,\"confidence\":[2018.9,2019.1]}},"
                    + "\"children\":[{\"name\":\"x\",\"node_attrs\":{\"num_date\":{\"value\":2019.5}}},"
                    + "{\"name\":\"MID\",\"node_attrs\":{\"num_date\":{\"value\":2020.0,\"confidence\":[2019.9,2020.5]}},"
                    + "\"children\":[{\"name\":\"a\",\"node_attrs\":{\"num_date\":{\"value\":2020.9}}},"
                    + "{\"name\":\"b\",\"node_attrs\":{\"num_date\":{\"value\":2021.0}}}]}]}}";
            final org.forester.io.parsers.json.AuspiceJsonParser parser = new org.forester.io.parsers.json.AuspiceJsonParser();
            parser.setSource( new StringBuffer( json ) );
            final Phylogeny phy = parser.parse()[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "cal" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true );
                    o.setShowHpdBars( true );
                    o.setNodeAgeShape( Options.NODE_AGE_SHAPE.BAR ); // measure the flat bar's extent, not a spindle
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    tp.setTimeAxisType( Options.TIME_AXIS_TYPE.CALENDAR );
                    final int w = 900, h = 460;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    org.forester.phylogeny.PhylogenyNode mid = null;
                    for ( final java.util.Iterator<org.forester.phylogeny.PhylogenyNode> it = phy.iteratorPreorder();
                            it.hasNext(); ) {
                        final org.forester.phylogeny.PhylogenyNode n = it.next();
                        if ( "MID".equals( n.getName() ) ) {
                            mid = n;
                            break;
                        }
                    }
                    if ( mid == null ) {
                        fail( ok, "expected the MID node" );
                        return;
                    }
                    final int nx = Math.round( mid.getXcoord() ), ny = Math.round( mid.getYcoord() );
                    int left_ext = 0, right_ext = 0;
                    for ( int dx = 1; dx <= 200; ++dx ) {
                        if ( bluishNear( img, nx - dx, nx - dx, ny ) ) {
                            left_ext = dx;
                        }
                        if ( bluishNear( img, nx + dx, nx + dx, ny ) ) {
                            right_ext = dx;
                        }
                    }
                    // later (right) reach must clearly exceed earlier (left) reach -- the calendar direction
                    if ( right_ext <= ( left_ext + 5 ) ) {
                        fail( ok, "a calendar tree's node-age bar must reach farther toward the LATER date (right "
                                + right_ext + " vs left " + left_ext + ")" );
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
                    // node-hpd-bars.xml has internal-node date intervals, so loading it must AUTO-ENABLE "Node Age
                    // Bars (HPD)" -- both the Option and its (kept-in-sync) menu item -- so the bars show without the
                    // user hunting for the toggle. Assert this before the test overrides the option below.
                    if ( !o.isShowHpdBars() ) {
                        ok[ 0 ] = false;
                        System.out.println( "  HPD bars must auto-enable when the loaded tree has date intervals" );
                    }
                    if ( ( frame._show_hpd_bars_cbmi != null ) && !frame._show_hpd_bars_cbmi.isSelected() ) {
                        ok[ 0 ] = false;
                        System.out.println( "  the HPD-bars menu item must be checked after the auto-enable" );
                    }
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
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    final int vertical_on = countBluish( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( vertical_on <= ( off + 300 ) ) {
                        fail( ok, "Node Age Bars should draw in a vertical orientation (on=" + vertical_on + " off="
                                + off + ")" );
                    }
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );

                    // NODE AGE SHAPE = SPINDLE: the tapered lens is TALLER at the node's point-estimate x (its peak)
                    // than a flat HPD_BAR_HEIGHT (7px) bar, so its vertical blue extent there exceeds 7 -- proving the
                    // spindle shape rendered (not a fallback to the flat bar)
                    if ( inode != null ) {
                        o.setNodeAgeShape( Options.NODE_AGE_SHAPE.SPINDLE );
                        tp.calcParametersForPainting( w, h );
                        final java.awt.image.BufferedImage spindle_img =
                                AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                        final int nx = Math.round( inode.getXcoord() ), ny = Math.round( inode.getYcoord() );
                        final int span = bluishVerticalSpan( spindle_img, nx, ny );
                        if ( span <= 7 ) {
                            fail( ok, "the SPINDLE must be taller than a flat 7px bar at the peak (blue span " + span
                                    + ")" );
                        }
                        o.setNodeAgeShape( Options.NODE_AGE_SHAPE.BAR );
                        tp.calcParametersForPainting( w, h );
                    }

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

    /** The vertical extent (px) of bluish overlay pixels in a thin x-window around {@code x}, near row {@code y} -- how
     *  tall the node-age overlay is at the point estimate (a flat bar is 7px; the spindle peaks taller). */
    private static int bluishVerticalSpan( final BufferedImage img, final int x, final int y ) {
        int min_y = Integer.MAX_VALUE, max_y = Integer.MIN_VALUE;
        for ( int yy = Math.max( 0, y - 12 ); yy < Math.min( img.getHeight(), y + 13 ); ++yy ) {
            for ( int xx = Math.max( 0, x - 1 ); xx <= Math.min( img.getWidth() - 1, x + 1 ); ++xx ) {
                final int rgb = img.getRGB( xx, yy );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( b >= ( r + 20 ) ) && ( b >= ( g + 15 ) ) ) {
                    min_y = Math.min( min_y, yy );
                    max_y = Math.max( max_y, yy );
                }
            }
        }
        return ( max_y >= min_y ) ? ( ( max_y - min_y ) + 1 ) : 0;
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
