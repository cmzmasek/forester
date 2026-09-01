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
import java.io.File;

import javax.swing.JFrame;
import javax.swing.JViewport;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Verifies the overview thumbnail is drawn in a VERTICAL orientation. Two things are checked: (1) the navigator
 * rectangle draws (the UI accent -- the navigator is chrome, so it no longer borrows the search-hit colour, and with
 * no search active the accent is the only thing painting those pixels); and (2) the ROTATED MINI-TREE itself
 * draws branches -- the navigator alone draws even if R is degenerate, so the mini-tree is detected by rendering the
 * viewport with the overview ON and OFF and counting the pixels that DIFFER but are not the navigator colour (the main
 * tree is identical between the two, so those differing pixels are the mini-tree's branches). Headful; a green no-op
 * when headless.
 */
public final class OrientationOverviewTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "OrientationOverview: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return overviewShownInVerticalView();
    }

    private static boolean overviewShownInVerticalView() {
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
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "ov" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    final ControlPanel cp = frame.getMainPanel().getControlPanel();
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    o.setShowOverview( true );
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    ( (JFrame) frame ).setSize( 520, 420 );
                    ( (JFrame) frame ).setVisible( true );
                    frame.showWhole();
                    // zoom both axes so the tree exceeds the viewport -> the overview turns on
                    for ( int i = 0; i < 5; ++i ) {
                        cp.zoomInScreenX( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                        cp.zoomInScreenY( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                    }
                    tp.updateOvSizes();
                    if ( !tp.isOvOn() ) {
                        fail( ok, "test setup: overview did not turn on after zooming past the viewport" );
                        return;
                    }
                    final JViewport vp = frame.getMainPanel().getCurrentScrollPane().getViewport();
                    vp.validate();
                    final BufferedImage img = renderViewport( vp ); // the visible tree + the overview at the corner
                    final Color nav = TreePanel.uiAccentColor(); // the navigator rectangle colour (UI chrome)
                    final int n = countColor( img, nav );
                    if ( n < 8 ) {
                        fail( ok, "the overview navigator rectangle was not drawn in the vertical view (accent px="
                                + n + ")" );
                    }
                    // ...and it must NOT be the search-hit colour: nothing is found here, so none may be painted
                    if ( countColor( img, tp.getTreeColorSet().getFoundColor0() ) > 0 ) {
                        fail( ok, "the navigator must not be drawn in the search-hit colour" );
                    }
                    // render again with the overview OFF: the main tree is identical, so every pixel that DIFFERS is
                    // overview ink; excluding the navigator colour leaves the rotated mini-tree's own branches, proving
                    // it actually drew (the navigator alone would draw even with a broken R -- finding #10).
                    o.setShowOverview( false );
                    vp.validate();
                    final BufferedImage img_off = renderViewport( vp );
                    final int mini = countMiniTreeInk( img, img_off, nav );
                    if ( mini < 20 ) {
                        fail( ok, "the rotated overview mini-tree drew (almost) no branches (mini-tree ink px=" + mini
                                + ")" );
                    }
                    o.setShowOverview( true );

                    // finding #9: a COLLAPSED clade's OWN incoming branch must still draw in the mini-tree. Collapse a
                    // clade with an internal child (its hidden descendants' branches vanish, but the apex's incoming
                    // limb must remain) and check for mini-tree ink at the collapsed apex's overview position.
                    org.forester.phylogeny.PhylogenyNode clade = null;
                    for ( final java.util.Iterator<org.forester.phylogeny.PhylogenyNode> it = phy.iteratorPreorder();
                            it.hasNext(); ) {
                        final org.forester.phylogeny.PhylogenyNode nd = it.next();
                        if ( nd.isExternal() || nd.isRoot() ) {
                            continue;
                        }
                        for ( int i = 0; i < nd.getNumberOfDescendants(); ++i ) {
                            if ( !nd.getChildNode( i ).isExternal() ) {
                                clade = nd;
                                break;
                            }
                        }
                        if ( clade != null ) {
                            break;
                        }
                    }
                    if ( clade == null ) {
                        fail( ok, "test setup: no internal clade with an internal child to collapse" );
                    }
                    else {
                        tp.collapse( clade );
                        for ( int i = 0; i < 3; ++i ) { // keep the (now smaller) tree past the viewport so the ov stays on
                            cp.zoomInScreenX( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                    AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                            cp.zoomInScreenY( AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                    AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                        }
                        tp.updateOvSizes();
                        vp.validate();
                        if ( !tp.isOvOn() ) {
                            fail( ok, "test setup: overview turned off after collapsing" );
                        }
                        else {
                            final BufferedImage con = renderViewport( vp ); // overview on, clade collapsed
                            final java.awt.geom.Point2D.Double apex = tp.overviewPointForTest( clade );
                            o.setShowOverview( false );
                            vp.validate();
                            final BufferedImage coff = renderViewport( vp ); // overview off
                            o.setShowOverview( true );
                            final int apex_ink = ( apex == null ) ? 0
                                    : countInkInBox( con, coff, nav, (int) Math.round( apex.x ),
                                            (int) Math.round( apex.y ), 5 );
                            if ( apex_ink < 1 ) {
                                fail( ok, "the collapsed clade's incoming branch is missing from the mini-tree "
                                        + "(finding #9): no ink at the apex overview point (" + apex + ")" );
                            }
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

    private static BufferedImage renderViewport( final JViewport vp ) {
        final BufferedImage img = new BufferedImage( Math.max( 1, vp.getWidth() ), Math.max( 1, vp.getHeight() ),
                BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, img.getWidth(), img.getHeight() );
        vp.printAll( g );
        g.dispose();
        return img;
    }

    /** Non-navigator differing pixels within a small box centred at (cx,cy) -- i.e. mini-tree ink at that spot. */
    private static int countInkInBox( final BufferedImage on, final BufferedImage off, final Color navColor,
                                      final int cx, final int cy, final int r ) {
        final int nr = navColor.getRed(), ng = navColor.getGreen(), nb = navColor.getBlue();
        int n = 0;
        for( int y = Math.max( 0, cy - r ); y <= Math.min( on.getHeight() - 1, cy + r ); ++y ) {
            for( int x = Math.max( 0, cx - r ); x <= Math.min( on.getWidth() - 1, cx + r ); ++x ) {
                final int a = on.getRGB( x, y ), b = off.getRGB( x, y );
                final int ar = ( a >> 16 ) & 0xFF, ag = ( a >> 8 ) & 0xFF, ab = a & 0xFF;
                final int br = ( b >> 16 ) & 0xFF, bg = ( b >> 8 ) & 0xFF, bb = b & 0xFF;
                final boolean differs = ( Math.abs( ar - br ) > 10 ) || ( Math.abs( ag - bg ) > 10 )
                        || ( Math.abs( ab - bb ) > 10 );
                final boolean is_nav = ( Math.abs( ar - nr ) <= 40 ) && ( Math.abs( ag - ng ) <= 40 )
                        && ( Math.abs( ab - nb ) <= 40 );
                if ( differs && !is_nav ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static int countColor( final BufferedImage img, final Color c ) {
        final int target = c.getRGB() & 0xFFFFFF;
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) == target ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Pixels that DIFFER between the overview-on and overview-off renders but are not the navigator colour -- i.e. the
     *  mini-tree's own (possibly faint, anti-aliased) branch ink. */
    private static int countMiniTreeInk( final BufferedImage on, final BufferedImage off, final Color navColor ) {
        final int nr = navColor.getRed(), ng = navColor.getGreen(), nb = navColor.getBlue();
        int n = 0;
        for( int y = 0; y < on.getHeight(); ++y ) {
            for( int x = 0; x < on.getWidth(); ++x ) {
                final int a = on.getRGB( x, y ), b = off.getRGB( x, y );
                final int ar = ( a >> 16 ) & 0xFF, ag = ( a >> 8 ) & 0xFF, ab = a & 0xFF;
                final int br = ( b >> 16 ) & 0xFF, bg = ( b >> 8 ) & 0xFF, bb = b & 0xFF;
                final boolean differs = ( Math.abs( ar - br ) > 10 ) || ( Math.abs( ag - bg ) > 10 )
                        || ( Math.abs( ab - bb ) > 10 );
                if ( !differs ) {
                    continue;
                }
                final boolean is_nav = ( Math.abs( ar - nr ) <= 40 ) && ( Math.abs( ag - ng ) <= 40 )
                        && ( Math.abs( ab - nb ) <= 40 );
                if ( !is_nav ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [OrientationOverviewTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [OrientationOverviewTest] " + msg );
        ok[ 0 ] = false;
    }

    private OrientationOverviewTest() {
    }
}
