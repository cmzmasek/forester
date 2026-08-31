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

import org.forester.archaeopteryx.ControlPanel.NodeClickAction;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The FOCUS GLOW: a halo marking the node under the pointer.
 * <p>
 * The point of the feature is that it appears in EVERY click-to mode. Before it, on-canvas hover feedback existed
 * only in Select-Node(s) mode, so in the DEFAULT mode -- and in the destructive ones (Reroot, Delete, Cut) -- the
 * hand cursor said "clickable" without ever saying WHICH node. So the load-bearing assertion here is the one for a
 * NON-select mode; the rest guard the colour semantics that let one indicator serve every mode.
 * <p>
 * Headful; a green no-op when headless.
 */
public final class FocusGlowRenderTest {

    private static final int W = 700, H = 420;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "FocusGlowRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        // Run the whole body ON THE EDT. Off it, a real repaint can interleave between two exports and re-run
        // calcParametersForPainting, which drifts the layout -- so the export-invariant check below compared two
        // differently-laid-out images and failed inside the suite while passing standalone.
        final boolean[] ok = { false };
        try {
            javax.swing.SwingUtilities.invokeAndWait( () -> ok[ 0 ] = glowOk() );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    private static BufferedImage shot( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.printAll( g );
        g.dispose();
        return img;
    }

    /** Pixels near the node that CHANGED versus the un-hovered render, split by dominant channel. */
    private static int[] glowInk( final BufferedImage off, final BufferedImage on, final int cx, final int cy ) {
        int bluish = 0, reddish = 0, neutral = 0;
        final int r = 40;
        for ( int y = Math.max( 0, cy - r ); y < Math.min( on.getHeight(), cy + r ); y++ ) {
            for ( int x = Math.max( 0, cx - r ); x < Math.min( on.getWidth(), cx + r ); x++ ) {
                if ( off.getRGB( x, y ) == on.getRGB( x, y ) ) {
                    continue;
                }
                final int rgb = on.getRGB( x, y );
                final int rr = ( rgb >> 16 ) & 0xff, gg = ( rgb >> 8 ) & 0xff, bb = rgb & 0xff;
                if ( ( bb > ( rr + 12 ) ) && ( bb > ( gg + 6 ) ) ) {
                    bluish++;
                }
                else if ( ( rr > ( bb + 12 ) ) && ( rr > ( gg + 6 ) ) ) {
                    reddish++;
                }
                else {
                    neutral++;
                }
            }
        }
        return new int[] { bluish, reddish, neutral };
    }

    /**
     * Changed pixels in an ANNULUS around the node, classified by dominant channel. The glow has to be judged out
     * here rather than over the node: an already-selected node carries its own solid red found marker, and counting
     * that in swamps the muted grey the glow is supposed to be showing.
     */
    private static int[] ringInk( final BufferedImage off, final BufferedImage on, final int cx, final int cy ) {
        int bluish = 0, reddish = 0, neutral = 0;
        final int r_min = 12, r_max = 24;
        for ( int y = Math.max( 0, cy - r_max ); y < Math.min( on.getHeight(), cy + r_max ); y++ ) {
            for ( int x = Math.max( 0, cx - r_max ); x < Math.min( on.getWidth(), cx + r_max ); x++ ) {
                final int dx = x - cx, dy = y - cy;
                final int d2 = ( dx * dx ) + ( dy * dy );
                if ( ( d2 < ( r_min * r_min ) ) || ( d2 > ( r_max * r_max ) ) ) {
                    continue;
                }
                if ( off.getRGB( x, y ) == on.getRGB( x, y ) ) {
                    continue;
                }
                final int rgb = on.getRGB( x, y );
                final int rr = ( rgb >> 16 ) & 0xff, gg = ( rgb >> 8 ) & 0xff, bb = rgb & 0xff;
                if ( ( bb > ( rr + 8 ) ) && ( bb > ( gg + 4 ) ) ) {
                    bluish++;
                }
                else if ( ( rr > ( bb + 8 ) ) && ( rr > ( gg + 4 ) ) ) {
                    reddish++;
                }
                else {
                    neutral++;
                }
            }
        }
        return new int[] { bluish, reddish, neutral };
    }

    private static int changed( final BufferedImage off, final BufferedImage on, final int cx, final int cy ) {
        final int[] c = glowInk( off, on, cx, cy );
        return c[ 0 ] + c[ 1 ] + c[ 2 ];
    }

    private static boolean glowOk() {
        MainFrameApplication mf = null;
        File nh = null;
        try {
            nh = File.createTempFile( "focusglow", ".nh" );
            java.nio.file.Files.writeString( nh.toPath(),
                                             "((AAAA:0.1,BBBB:0.1)N1:0.1,(CCCC:0.1,DDDD:0.1)N2:0.1)root;\n" );
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( nh, new NHXParser() )[ 0 ];
            mf = ( MainFrameApplication ) MainFrameApplication.createInstance( new Phylogeny[] { phy },
                                                                              new Configuration(), "focusglow" );
            mf.setSize( W, H );
            // NOTE: deliberately NOT pinning the theme with setDarkMode() -- reinstalling the look-and-feel
            // mid-suite leaves FlatLaf's tabbed pane in a broken state (a null tabInsets NPE, then a
            // "missing initial moveto" from the next paint). It is not needed: the found colour is red and the
            // accent blue in BOTH themes (the hues were unified in 0.11.42), so the channel comparisons below
            // hold either way, and the grey de-select mark is theme-independent.
            final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
            final ControlPanel cp = mf.getMainPanel().getControlPanel();
            tp.setSize( W, H );
            // The rollover node-description popup is HEAVYWEIGHT near a screen edge, and PopupFactory CACHES
            // its native window instead of disposing it -- a displayable window keeps AWT's non-daemon threads
            // alive, so the suite JVM prints all its results and then never exits. This test drives mouseMoved
            // below, which would open one. (Same trap as the standalone-quit hang fixed in 0.9.68.)
            if ( cp.getNodeDescPopupCb() != null ) {
                cp.getNodeDescPopupCb().setSelected( false );
            }
            cp.showWhole();
            final PhylogenyNode target = phy.getNode( "BBBB" );
            shot( tp ); // one render so the node has coordinates at all (they are assigned during PAINT)

            // ---- (1) hover is transient screen state: it must NEVER reach an export ----------------------
            tp.setHoverForTest( null, false );
            cp.showWhole();
            final PhylogenyNode t2 = target;
            final File png_off = File.createTempFile( "glowexpoff", ".png" );
            final File png_on = File.createTempFile( "glowexpon", ".png" );
            AptxUtil.writePhylogenyToGraphicsFile( png_off.getAbsolutePath(), W, H, tp, cp,
                                                   AptxUtil.GraphicsExportType.PNG, tp.getOptions() );
            tp.setHoverForTest( t2, false );
            AptxUtil.writePhylogenyToGraphicsFile( png_on.getAbsolutePath(), W, H, tp, cp,
                                                   AptxUtil.GraphicsExportType.PNG, tp.getOptions() );
            final BufferedImage e_off = javax.imageio.ImageIO.read( png_off );
            final BufferedImage e_on = javax.imageio.ImageIO.read( png_on );
            int diff = 0;
            for ( int y = 0; y < Math.min( e_off.getHeight(), e_on.getHeight() ); y++ ) {
                for ( int x = 0; x < Math.min( e_off.getWidth(), e_on.getWidth() ); x++ ) {
                    if ( e_off.getRGB( x, y ) != e_on.getRGB( x, y ) ) {
                        diff++;
                    }
                }
            }
            png_off.delete();
            png_on.delete();
            if ( diff != 0 ) {
                return fail( "the focus glow must not appear in an export, " + diff + " px differ" );
            }
            // ---- (2) THE POINT OF THE FEATURE: a non-select mode glows at all ----------------------------
            cp.setActionWhenNodeClicked( NodeClickAction.SHOW_DATA );
            tp.setHoverForTest( null, false );
            final BufferedImage off = shot( tp );
            // read the coordinates HERE, not earlier: an export re-lays the tree out at the export size, so a
            // position captured before the export block above would be stale by now
            final int cx = Math.round( target.getXcoord() ), cy = Math.round( target.getYcoord() );
            tp.setHoverForTest( target, false );
            final BufferedImage plain = shot( tp );
            final int[] accent = glowInk( off, plain, cx, cy );
            if ( changed( off, plain, cx, cy ) < 60 ) {
                return fail( "hovering a node in a NON-select mode must draw a focus glow, changed only "
                        + changed( off, plain, cx, cy ) + " px" );
            }
            // ... in the neutral accent, NOT the found colour -- a focus ring must not read as a selection state
            if ( accent[ 0 ] <= accent[ 1 ] ) {
                return fail( "the glow outside Select mode must be the accent colour, got bluish " + accent[ 0 ]
                        + " vs reddish " + accent[ 1 ] );
            }
            // ... and it must reach BEYOND the old flat preview disc (max(6, node size)) -- it is a halo
            if ( !inkAt( off, plain, cx + 13, cy ) && !inkAt( off, plain, cx - 13, cy ) ) {
                return fail( "the glow must extend well past the node, like a halo" );
            }

            // ---- (3) in Select mode the SAME glow carries the selection meaning --------------------------
            cp.setActionWhenNodeClicked( NodeClickAction.SELECT_NODES );
            tp.setHoverForTest( null, false );
            final BufferedImage sel_off = shot( tp );
            tp.setHoverForTest( target, false );
            final int[] add = ringInk( sel_off, shot( tp ), cx, cy );
            if ( add[ 1 ] <= ( add[ 0 ] + add[ 2 ] ) ) {
                return fail( "in Select mode a click that ADDS must glow in the found colour, got reddish "
                        + add[ 1 ] + " vs bluish " + add[ 0 ] + " / neutral " + add[ 2 ] );
            }

            // ---- (3b) a shift-click outside Select mode must NOT blank the focus ring ---------------------
            // shift-click arms the just-clicked suppression in EVERY mode, but that exists only to stop a
            // selection preview flipping to grey under a stationary pointer; a neutral focus mark has nothing to
            // flip, so blanking it there would just make the ring vanish under the pointer.
            cp.setActionWhenNodeClicked( NodeClickAction.SHOW_DATA );
            tp.setHoverForTest( null, false );
            tp.mouseClicked( new java.awt.event.MouseEvent( tp, java.awt.event.MouseEvent.MOUSE_CLICKED, 0,
                                                            java.awt.event.InputEvent.SHIFT_DOWN_MASK, cx, cy, 1,
                                                            false ) );
            tp.mouseMoved( new java.awt.event.MouseEvent( tp, java.awt.event.MouseEvent.MOUSE_MOVED, 0, 0, cx, cy, 0,
                                                          false ) );
            if ( tp.hoverNodeForTest() != target ) {
                return fail( "a shift-click outside Select mode must leave the focus glow on the node" );
            }
            tp.setFoundNodes0( null );

            // ---- (4) hovering an ALREADY-SELECTED node glows grey ("a click removes it") -----------------
            cp.setActionWhenNodeClicked( NodeClickAction.SELECT_NODES ); // (3b) left the panel in SHOW_DATA
            tp.setHoverForTest( null, false );
            final java.util.Set<Long> picked_set = new java.util.HashSet<Long>();
            picked_set.add( target.getId() );
            tp.setFoundNodes0( picked_set );
            final BufferedImage picked = shot( tp );
            tp.setHoverForTest( target, false );
            final int[] remove = ringInk( picked, shot( tp ), cx, cy );
            if ( remove[ 2 ] <= ( remove[ 1 ] * 2 ) ) {
                return fail( "hovering a SELECTED node must glow in the muted grey, got neutral " + remove[ 2 ]
                        + " vs reddish " + remove[ 1 ] + " (ring)" );
            }
            tp.setFoundNodes0( null );

            // ---- (4) a COLLAPSED clade glows too (it is clickable in every mode) -------------------------
            cp.setActionWhenNodeClicked( NodeClickAction.SHOW_DATA );
            final PhylogenyNode n2 = phy.getNode( "N2" );
            n2.setCollapse( true );
            cp.displayedPhylogenyMightHaveChanged( true );
            tp.calcParametersForPainting( W, H ); // explicit: showWhole() needs a realised viewport
            tp.setHoverForTest( null, false );
            final BufferedImage coll_off = shot( tp );
            tp.setHoverForTest( n2, false );
            final BufferedImage coll_on = shot( tp ); // paint first: it is what assigns the collapsed node's coords
            final int ccx = Math.round( n2.getXcoord() ), ccy = Math.round( n2.getYcoord() );
            if ( changed( coll_off, coll_on, ccx, ccy ) < 60 ) {
                return fail( "a collapsed clade must glow on hover as well" );
            }
            // ... and its glow must read the CLADE's state, not the clade root's own id: a selection click on a
            // triangle routes to selectSubtreeTips, so a fully-selected collapsed clade is about to be CLEARED.
            // Reading the root (whose id is never in the found set) would show "a click adds" and invert the cue.
            cp.setActionWhenNodeClicked( NodeClickAction.SELECT_NODES );
            tp.setFoundNodes0( null );
            tp.selectSubtreeTips( n2 );
            tp.setHoverForTest( null, false );
            tp.setHoverForTest( n2, false );
            if ( !tp.hoverWouldDeselectForTest() ) {
                return fail( "a fully-selected COLLAPSED clade must glow as 'a click removes', not 'adds'" );
            }
            tp.setFoundNodes0( null );
            tp.setHoverForTest( null, false );

            n2.setCollapse( false );
            cp.displayedPhylogenyMightHaveChanged( true );
            tp.calcParametersForPainting( W, H );

            return true;
        }
        catch ( final Exception e ) {
            return fail( "unexpected exception: " + e );
        }
        finally {
            if ( mf != null ) {
                mf.dispose();
            }
            if ( nh != null ) {
                nh.delete();
            }
        }
    }

    private static boolean inkAt( final BufferedImage off, final BufferedImage on, final int x, final int y ) {
        return ( x >= 0 ) && ( y >= 0 ) && ( x < on.getWidth() ) && ( y < on.getHeight() )
                && ( off.getRGB( x, y ) != on.getRGB( x, y ) );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [FocusGlowRenderTest] " + msg );
        return false;
    }
}
