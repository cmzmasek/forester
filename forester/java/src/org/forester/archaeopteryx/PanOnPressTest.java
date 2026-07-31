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
import java.awt.Point;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Regression guard for the removal of the press-time pan block in {@link MouseListener}: a plain left-button
 * press (with no movement) must NOT start a pan/overview drag and must NOT scroll the viewport -- panning begins
 * only on the first {@code mouseDragged} (which initializes the drag anchor). Guarded to a no-op when headless
 * (needs FlatLaf via {@code createInstance}).
 */
public final class PanOnPressTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "PanOnPress: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, conf, "pan" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final int w = 600, h = 400;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );

                    final JScrollPane sp = frame.getMainPanel().getCurrentScrollPane();
                    final Point before = ( sp != null ) ? sp.getViewport().getViewPosition() : null;
                    final MouseListener ml = new MouseListener( tp );
                    final int px = w / 2, py = h / 2; // middle of the tree, away from the legend/overview corners

                    // A plain left-button press must NOT initiate a pan/drag. The isDraggingForTest() check is the
                    // discriminating one (the removed block set _being_dragged on press); the viewport check is a
                    // belt-and-suspenders guard against a future press-time scroll (the removed block scrolled by 0).
                    ml.mousePressed( ev( tp, MouseEvent.MOUSE_PRESSED, InputEvent.BUTTON1_DOWN_MASK, px, py ) );
                    if ( ml.isDraggingForTest() ) {
                        fail( ok, "a press without movement must not start a pan/drag" );
                    }
                    if ( ( before != null ) && !sp.getViewport().getViewPosition().equals( before ) ) {
                        fail( ok, "a press without movement must not scroll the viewport" );
                    }

                    // a real drag DOES engage the pan path (so panning still works after the removal)
                    ml.mouseDragged( ev( tp, MouseEvent.MOUSE_DRAGGED, InputEvent.BUTTON1_DOWN_MASK, px + 40, py + 25 ) );
                    if ( !ml.isDraggingForTest() ) {
                        fail( ok, "a drag must engage the pan/drag path" );
                    }

                    // release ends the drag
                    ml.mouseReleased( ev( tp, MouseEvent.MOUSE_RELEASED, 0, px + 40, py + 25 ) );
                    if ( ml.isDraggingForTest() ) {
                        fail( ok, "release must end the drag" );
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

    private static MouseEvent ev( final TreePanel tp, final int id, final int mods, final int x, final int y ) {
        return new MouseEvent( tp, id, System.currentTimeMillis(), mods, x, y, 1, false );
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [PanOnPressTest] " + msg );
        ok[ 0 ] = false;
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 6; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private PanOnPressTest() {
    }
}
