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
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Verifies "Flip Vertically" on the ladder demo (forester/demo/flip-vertically.xml): with the toggle off tip_01 is
 * above tip_08; with it on their vertical order reverses (a true reflection about the tree centre), and -- crucially --
 * hit-testing follows, i.e. findNode at each tip's on-screen position still returns that tip. Headful; a green no-op
 * when headless. Dogfoods the demo.
 */
public final class FlipVerticallyRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "FlipVerticallyRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return flipRendersOk();
    }

    private static boolean flipRendersOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/flip-vertically.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "flip" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final PhylogenyNode t1 = tipNamed( phy, "tip_01" );
                    final PhylogenyNode t8 = tipNamed( phy, "tip_08" );
                    if ( ( t1 == null ) || ( t8 == null ) ) {
                        fail( ok, "demo must contain tip_01 and tip_08" );
                        return;
                    }
                    final int w = 560, h = 460;
                    // OFF: normal order -- tip_01 above tip_08
                    o.setFlipVertically( false );
                    layout( frame, tp, o, w, h );
                    final float y1_off = t1.getYcoord(), y8_off = t8.getYcoord();
                    if ( !( y1_off < y8_off ) ) {
                        fail( ok, "without flip, tip_01 (" + y1_off + ") must be above tip_08 (" + y8_off + ")" );
                    }
                    // ON: reversed -- tip_01 below tip_08, and a true reflection (each swaps to the other's row)
                    o.setFlipVertically( true );
                    layout( frame, tp, o, w, h );
                    final float y1_on = t1.getYcoord(), y8_on = t8.getYcoord();
                    if ( !( y1_on > y8_on ) ) {
                        fail( ok, "with flip, tip_01 (" + y1_on + ") must be below tip_08 (" + y8_on + ")" );
                    }
                    if ( ( Math.abs( y1_on - y8_off ) > 3f ) || ( Math.abs( y8_on - y1_off ) > 3f ) ) {
                        fail( ok, "flip must be a reflection: tip_01/tip_08 should swap rows (y1_on=" + y1_on
                                + " y8_off=" + y8_off + ", y8_on=" + y8_on + " y1_off=" + y1_off + ")" );
                    }
                    // hit-testing must follow the flip: a click at a tip's drawn position still returns that tip
                    final PhylogenyNode h1 = tp.findNode( Math.round( t1.getXcoord() ), Math.round( y1_on ) );
                    final PhylogenyNode h8 = tp.findNode( Math.round( t8.getXcoord() ), Math.round( y8_on ) );
                    if ( ( h1 != t1 ) || ( h8 != t8 ) ) {
                        fail( ok, "hit-testing must follow the flip: findNode returned "
                                + ( h1 == null ? "null" : h1.getName() ) + " / "
                                + ( h8 == null ? "null" : h8.getName() ) + " (expected tip_01 / tip_08)" );
                    }
                    // the OVERVIEW thumbnail lays out its OWN (YSecondary) coords, so it must be flipped too -- else
                    // the mini-map mirrors the canvas and its navigator maps to the wrong region. The overview only
                    // activates when the tree is larger than the viewport, so realize the frame + zoom in first.
                    ( (JFrame) frame ).setSize( 760, 620 );
                    ( (JFrame) frame ).validate();
                    o.setShowOverview( true );
                    final int bw = 1600, bh = 1600;
                    final float ov1_off = overviewTipY( tp, o, t1, false, bw, bh );
                    final float ov8_off = t8.getYSecondary();
                    final float ov1_on = overviewTipY( tp, o, t1, true, bw, bh );
                    final float ov8_on = t8.getYSecondary();
                    if ( tp.isOvOn() && ( ( ov1_off >= ov8_off ) || ( ov1_on <= ov8_on ) ) ) {
                        fail( ok, "the overview thumbnail must follow the flip (off: t1=" + ov1_off + " t8=" + ov8_off
                                + ", on: t1=" + ov1_on + " t8=" + ov8_on + ")" );
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

    /** Sets flip, zooms the tree larger than the viewport, and paints via the SCREEN path so the overview
     *  ({@code paintPhylogenyLite}) runs and assigns each tip's YSecondary; returns the given tip's YSecondary. */
    private static float overviewTipY( final TreePanel tp, final Options o, final PhylogenyNode tip,
                                       final boolean flip, final int w, final int h ) {
        o.setFlipVertically( flip );
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        tp.updateOvSizes(); // turns the overview on when the tree exceeds the viewport + computes its layout
        final java.awt.image.BufferedImage img = new java.awt.image.BufferedImage( w, h,
                java.awt.image.BufferedImage.TYPE_INT_RGB );
        final java.awt.Graphics2D g = img.createGraphics();
        try {
            tp.paintComponent( g ); // screen path (not to_graphics_file) -> paintPhylogenyLite assigns YSecondary
        }
        finally {
            g.dispose();
        }
        return tip.getYSecondary();
    }

    private static void layout( final MainFrame frame, final TreePanel tp, final Options o, final int w, final int h ) {
        frame.showWhole();
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        // render once (discarding the image) so the node coords are assigned -- they are set during the paint loop
        AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
    }

    private static PhylogenyNode tipNamed( final Phylogeny phy, final String name ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        return null;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [FlipVerticallyRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [FlipVerticallyRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private FlipVerticallyRenderTest() {
    }
}
