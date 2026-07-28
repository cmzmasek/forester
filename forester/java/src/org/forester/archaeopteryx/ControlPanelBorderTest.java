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

import java.awt.Component;
import java.awt.GraphicsEnvironment;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;

import org.forester.phylogeny.Phylogeny;

/**
 * Regression guard for the "control panel has no border on a fresh, tree-less startup" bug. The control panel is
 * wrapped in a {@link JScrollPane} placed in the WEST region; the active look-and-feel installs a thin scroll-pane
 * border (a separator between the control panel and the tree canvas) at construction. {@code MainPanel} used to
 * null that border out, so it reappeared only after the first runtime theme switch (which reinstalls it via
 * {@code updateComponentTreeUI}) or a tree load -- leaving the initial window looking "incomplete". The fix keeps
 * the L&F border. This test asserts the scroll pane already has a (non-null) border on a freshly-created,
 * tree-less {@link MainFrame}, and that simulating a theme reinstall does not change it (startup already matches
 * the themed look). Guarded to a no-op on a headless box (needs FlatLaf via {@code createInstance}).
 */
public final class ControlPanelBorderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ControlPanelBorder: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] {}, new Configuration(),
                                                                        "cp-border" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    exercise( frame, ok );
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

    private static void exercise( final MainFrame frame, final boolean[] ok ) {
        frame.setSize( 1000, 600 );
        ( (JFrame) frame ).validate();
        final JScrollPane scroller = controlPanelScroller( frame );
        if ( scroller == null ) {
            fail( ok, "could not locate the control-panel scroll pane" );
            return;
        }
        // THE REGRESSION: on a fresh, tree-less startup the scroll pane must already carry the L&F's border.
        // A null border here is exactly the "incomplete"-looking startup the fix addresses.
        final Border fresh = scroller.getBorder();
        if ( fresh == null ) {
            fail( ok, "control-panel scroll pane has NO border at startup (window looks incomplete)" );
        }
        // A runtime theme switch reinstalls component UIs via updateComponentTreeUI; that is what used to make the
        // border appear. It must NOT change the border now -- startup already matches the themed look.
        SwingUtilities.updateComponentTreeUI( (JFrame) frame );
        final Border after = scroller.getBorder();
        if ( after == null ) {
            fail( ok, "control-panel scroll pane border went missing after a theme reinstall" );
        }
        if ( ( fresh != null ) && ( after != null ) && !fresh.getClass().equals( after.getClass() ) ) {
            fail( ok, "startup border (" + fresh.getClass().getSimpleName() + ") differs from the themed border ("
                    + after.getClass().getSimpleName() + ")" );
        }
    }

    // the WEST scroll pane whose viewport view is the ControlPanel (the CENTER child is the JTabbedPane)
    private static JScrollPane controlPanelScroller( final MainFrame frame ) {
        final MainPanel mp = frame.getMainPanel();
        final ControlPanel cp = mp.getControlPanel();
        for( final Component c : mp.getComponents() ) {
            if ( c instanceof JScrollPane ) {
                final JScrollPane sp = (JScrollPane) c;
                if ( ( sp.getViewport() != null ) && ( sp.getViewport().getView() == cp ) ) {
                    return sp;
                }
            }
        }
        return null;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ControlPanelBorderTest] " + msg );
        ok[ 0 ] = false;
    }

    private ControlPanelBorderTest() {
    }
}
