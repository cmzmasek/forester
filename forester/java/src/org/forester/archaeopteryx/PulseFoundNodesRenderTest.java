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
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Verifies "Pulse Found Nodes" on the search-emphasis demo: an EXPORT renders a static found-color glow around a hit
 * (absent when off); an on-screen paint draws the same halo and it BREATHES across phases; and the animation timer
 * has a sane lifecycle (stopped after the panel leaves the hierarchy). Headful; a green no-op when headless.
 */
public final class PulseFoundNodesRenderTest {

    private static final int W = 640, H = 560;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "PulseFoundNodesRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return pulseOk();
    }

    private static boolean pulseOk() {
        final File file = new File( System.getProperty( "user.dir" ), "forester/demo/search-emphasis.xml" );
        if ( !file.exists() ) {
            return fail( "demo tree missing: " + file.getAbsolutePath() );
        }
        final MainFrame[] mf = new MainFrame[ 1 ];
        final boolean[] ok = { true };
        try {
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "pulse" ) );
            final PhylogenyNode[] akt = new PhylogenyNode[ 1 ];
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                final Options o = frame.getOptions();
                o.setGraphicsExportWhiteBackground( true ); // LIGHT theme -> found color is red -> pink halo on white
                tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                final Set<Long> found = new HashSet<>();
                for( final PhylogenyNode n : phy.getExternalNodes() ) {
                    if ( n.getName().equals( "AKT1_kinase" ) ) {
                        akt[ 0 ] = n;
                    }
                    if ( n.getName().contains( "kinase" ) ) {
                        found.add( n.getId() );
                    }
                }
                tp.setFoundNodes0( found );
                o.setPulseFoundNodes( true );
                frame.showWhole();
                tp.setSize( W, H );
                tp.calcParametersForPainting( W, H );
                // (1) EXPORT static glow: a found-color halo appears just left of the hit's node (absent when off)
                final int glow_on = haloInk( AptxUtil.renderPhylogenyToImage( W, H, tp, o, false, 1, false ), akt[ 0 ] );
                o.setPulseFoundNodes( false );
                final int glow_off = haloInk( AptxUtil.renderPhylogenyToImage( W, H, tp, o, false, 1, false ),
                        akt[ 0 ] );
                o.setPulseFoundNodes( true );
                // off still has a few px from the base found-highlight dot's antialiased edge -- the halo must ADD
                // clearly more on top; the meaningful assertion is on >> off.
                if ( glow_off > 40 ) {
                    fail( ok, "no halo should render when Pulse Found Nodes is off, got " + glow_off );
                }
                if ( glow_on <= ( glow_off + 25 ) ) {
                    fail( ok, "an export should render a static found-node glow (halo px on=" + glow_on + " off="
                            + glow_off + ")" );
                }
            } );
            if ( akt[ 0 ] == null ) {
                fail( ok, "demo must contain AKT1_kinase" );
                return ok[ 0 ];
            }
            // (2) ON-SCREEN animation: paint at 3 phases across the pulse period; the halo must BREATHE (vary)
            final int[] counts = new int[ 3 ];
            for( int i = 0; i < 3; i++ ) {
                final int fi = i;
                SwingUtilities.invokeAndWait(
                        () -> counts[ fi ] = haloInk( paint( mf[ 0 ].getMainPanel().getCurrentTreePanel() ), akt[ 0 ] ) );
                Thread.sleep( 430 ); // ~1/3 of PULSE_PERIOD, off the EDT, so the next paint lands at a new phase
            }
            int lo = counts[ 0 ], hi = counts[ 0 ];
            for( final int c : counts ) {
                lo = Math.min( lo, c );
                hi = Math.max( hi, c );
            }
            if ( ( hi - lo ) < 10 ) {
                fail( ok, "the on-screen halo must breathe across phases (counts " + counts[ 0 ] + "," + counts[ 1 ]
                        + "," + counts[ 2 ] + ")" );
            }
            // (3) LEAK-FIX: a rectangular screen paint records halo repaint regions (the timer would run); switching
            // to a CIRCULAR layout clears the halo state on the next paint, so the timer reconciles to stopped --
            // the fix for the "timer keeps ticking on stale regions after a layout switch" bug.
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.calcParametersForPainting( W, H );
                paint( tp );
                final int rect_bounds = tp.getFoundHaloBoundsCountForTest();
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                tp.calcParametersForPainting( W, H );
                paint( tp );
                final int circ_bounds = tp.getFoundHaloBoundsCountForTest();
                if ( rect_bounds <= 0 ) {
                    fail( ok, "a rectangular paint should record found-node halo regions, got " + rect_bounds );
                }
                if ( circ_bounds != 0 ) {
                    fail( ok, "switching to a circular layout must clear the halo state so the pulse timer stops "
                            + "(leftover regions=" + circ_bounds + ")" );
                }
            } );
            // (4) LIFECYCLE smoke: after the panel is disposed the animation timer must not be running (no leak)
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                ( (JFrame) mf[ 0 ] ).dispose(); // -> removeNotify() stops the timer
                if ( tp.isPulseTimerRunning() ) {
                    fail( ok, "the pulse timer must be stopped once the panel/window is disposed" );
                }
            } );
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            fail( ok, "unexpected: " + t );
        }
        finally {
            if ( mf[ 0 ] != null ) {
                try {
                    SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() ); // best-effort teardown
                }
                catch ( final Exception ignore ) {
                    // ignore teardown failure
                }
            }
        }
        return ok[ 0 ];
    }

    private static BufferedImage paint( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            tp.paintComponent( g ); // screen path -> time-based halo phase + halo animation state recorded
        }
        finally {
            g.dispose();
        }
        return img;
    }

    /** Count pink-red halo pixels just LEFT of a hit's node (over the white/branch background, away from the red
     *  label to the right AND clear of the node dot's own antialiased edge, which would otherwise count when off):
     *  translucent red over white -> (255, 255-a, 255-a), i.e. r high and r-g positive. */
    private static int haloInk( final BufferedImage img, final PhylogenyNode node ) {
        final int x = Math.round( node.getXcoord() ), y = Math.round( node.getYcoord() );
        int n = 0;
        for( int yy = Math.max( 0, y - 11 ); yy < Math.min( img.getHeight(), y + 12 ); ++yy ) {
            for( int xx = Math.max( 0, x - 13 ); xx < Math.min( img.getWidth(), x - 3 ); ++xx ) {
                final int rgb = img.getRGB( xx, yy );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( r >= 230 ) && ( ( r - g ) >= 15 ) && ( Math.abs( g - b ) <= 12 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [PulseFoundNodesRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [PulseFoundNodesRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private PulseFoundNodesRenderTest() {
    }
}
