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
import java.awt.Window;
import java.awt.event.ComponentEvent;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The node-description rollover popup must NEVER be stranded.
 * <p>
 * The popup can be a HEAVYWEIGHT window -- its own native window -- so any path that leaves the canvas without
 * hiding it leaves a phantom popup floating on the desktop, even outside the application window (the
 * long-standing user-reported annoyance). This drives every exit path through the REAL registered listeners
 * (dispatchEvent, not internal calls): mouse exit, mouse press, mouse wheel, key press, window deactivation,
 * window move, the display option being switched off while a popup is showing, and the panel being removed.
 * <p>
 * SUITE-HYGIENE NOTE: this test deliberately opens real popups -- the exact machinery documented to hang the
 * suite JVM (PopupFactory CACHES a heavyweight popup's native window instead of disposing it, keeping AWT's
 * non-daemon threads alive after all tests pass). The finally block therefore disposes every leftover popup
 * window by hand. Do not remove that sweep.
 */
public final class NodePopupLifecycleTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "Node popup lifecycle: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        try {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { tree() }, new Configuration(), "popup-lifecycle" ) );
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                // deterministic layout regardless of the developer's persisted display settings, and no
                // overview thumbnail sitting over the node we aim mouse events at
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.getOptions().setShowOverview( false );
                tp.setOvOn( false );
                mf[ 0 ].getMainPanel().getControlPanel().showWhole();
                if ( !tp.shows( DisplayOption.NODE_DATA_POPUP ) ) {
                    ok[ 0 ] = fail( "precondition: the Rollover display option must default ON" );
                    return;
                }
                // node coordinates are assigned during PAINT, so render once before aiming events at a node
                final BufferedImage img = new BufferedImage( 600, 460, BufferedImage.TYPE_INT_RGB );
                final Graphics2D g = img.createGraphics();
                tp.printAll( g );
                g.dispose();
                final PhylogenyNode tip = tp.getPhylogeny().getFirstExternalNode();
                final int nx = ( int ) tip.getXcoord();
                final int ny = ( int ) tip.getYcoord();

                // (1) leaving the canvas hides it -- THE everyday leak
                if ( !show( tp, nx, ny ) ) {
                    ok[ 0 ] = fail( "precondition: a mouse move over a named tip must show the popup" );
                    return;
                }
                tp.dispatchEvent( new MouseEvent( tp, MouseEvent.MOUSE_EXITED, now(), 0, -5, -5, 0, false ) );
                if ( tp.isNodeDescPopupShowingForTest() ) {
                    ok[ 0 ] = fail( "leaving the canvas must hide the rollover popup" );
                }

                // (2) a mouse press hides it (a dialog/menu is about to come up)
                show( tp, nx, ny );
                tp.dispatchEvent( new MouseEvent( tp, MouseEvent.MOUSE_PRESSED, now(),
                                                  InputEvent.BUTTON1_DOWN_MASK, nx, ny, 1, false,
                                                  MouseEvent.BUTTON1 ) );
                if ( tp.isNodeDescPopupShowingForTest() ) {
                    ok[ 0 ] = fail( "a mouse press must hide the rollover popup" );
                }

                // (3) the wheel scrolls the tree under a stationary pointer -- popup must not go stale
                show( tp, nx, ny );
                tp.dispatchEvent( new MouseWheelEvent( tp, MouseEvent.MOUSE_WHEEL, now(), 0, nx, ny, 0, false,
                                                       MouseWheelEvent.WHEEL_UNIT_SCROLL, 1, 1 ) );
                if ( tp.isNodeDescPopupShowingForTest() ) {
                    ok[ 0 ] = fail( "a mouse-wheel scroll must hide the rollover popup" );
                }

                // (4) so does keyboard scrolling (and Esc, and any other key)
                show( tp, nx, ny );
                tp.dispatchEvent( new KeyEvent( tp, KeyEvent.KEY_PRESSED, now(), 0, KeyEvent.VK_DOWN,
                                                KeyEvent.CHAR_UNDEFINED ) );
                if ( tp.isNodeDescPopupShowingForTest() ) {
                    ok[ 0 ] = fail( "a key press must hide the rollover popup" );
                }

                // (5) Cmd-Tab / another window taking focus: no mouse event fires, the WINDOW must hide it.
                // A synthetic WINDOW_DEACTIVATED through dispatchEvent is CONSUMED by the KeyboardFocusManager
                // whenever the frame is not the OS's actual active window (always, in a test run) -- so invoke
                // the registered window listeners directly, exactly as processWindowEvent would on the real event.
                show( tp, nx, ny );
                final WindowEvent deact = new WindowEvent( mf[ 0 ], WindowEvent.WINDOW_DEACTIVATED );
                for( final java.awt.event.WindowListener wl : mf[ 0 ].getWindowListeners() ) {
                    wl.windowDeactivated( deact );
                }
                if ( tp.isNodeDescPopupShowingForTest() ) {
                    ok[ 0 ] = fail( "deactivating the window must hide the rollover popup" );
                }

                // (5b) minimizing the window (WINDOW_ICONIFIED passes dispatchEvent uninspected by the KFM)
                show( tp, nx, ny );
                mf[ 0 ].dispatchEvent( new WindowEvent( mf[ 0 ], WindowEvent.WINDOW_ICONIFIED ) );
                if ( tp.isNodeDescPopupShowingForTest() ) {
                    ok[ 0 ] = fail( "iconifying the window must hide the rollover popup" );
                }

                // (6) dragging the window: the popup must not stay behind at the old spot
                show( tp, nx, ny );
                mf[ 0 ].dispatchEvent( new ComponentEvent( mf[ 0 ], ComponentEvent.COMPONENT_MOVED ) );
                if ( tp.isNodeDescPopupShowingForTest() ) {
                    ok[ 0 ] = fail( "moving the window must hide the rollover popup" );
                }

                // (7) THE GATE BUG: switching the option off while a popup is showing used to strand it
                // permanently (the only hide was itself gated on the option being ON)
                show( tp, nx, ny );
                tp.setShows( DisplayOption.NODE_DATA_POPUP, false );
                tp.dispatchEvent( new MouseEvent( tp, MouseEvent.MOUSE_MOVED, now(), 0, nx, ny, 0, false ) );
                if ( tp.isNodeDescPopupShowingForTest() ) {
                    ok[ 0 ] = fail( "a popup showing when the Rollover option is switched OFF must still hide" );
                }
                tp.setShows( DisplayOption.NODE_DATA_POPUP, true );

                // (8) closing the tab/window (removeNotify) hides it
                show( tp, nx, ny );
                mf[ 0 ].dispose();
                if ( tp.isNodeDescPopupShowingForTest() ) {
                    ok[ 0 ] = fail( "disposing the window must hide the rollover popup" );
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        finally {
            try {
                SwingUtilities.invokeAndWait( () -> {
                    if ( mf[ 0 ] != null ) {
                        mf[ 0 ].dispose();
                    }
                    // PopupFactory CACHES heavyweight popup windows instead of disposing them; a cached native
                    // window keeps AWT's non-daemon threads alive and hangs the suite JVM after all tests pass
                    // (the documented gotcha). Dispose them by hand.
                    for( final Window w : Window.getWindows() ) {
                        if ( w.getClass().getName().contains( "HeavyWeightWindow" ) ) {
                            w.dispose();
                        }
                    }
                } );
            }
            catch ( final Exception e ) {
                // best effort
            }
        }
    }

    /** Shows the popup by dispatching a REAL mouse move over the tip; true iff it is then showing. */
    private static boolean show( final TreePanel tp, final int x, final int y ) {
        tp.dispatchEvent( new MouseEvent( tp, MouseEvent.MOUSE_MOVED, now(), 0, x, y, 0, false ) );
        return tp.isNodeDescPopupShowingForTest();
    }

    private static long now() {
        return System.currentTimeMillis();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [NodePopupLifecycleTest] " + msg );
        return false;
    }

    private static Phylogeny tree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 4; ++i ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "tip_" + i ); // a non-empty name, so the popup has content to show
            n.setDistanceToParent( 0.2 );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private NodePopupLifecycleTest() {
    }
}
