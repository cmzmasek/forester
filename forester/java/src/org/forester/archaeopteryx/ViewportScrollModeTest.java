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

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JViewport;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Locks in the fix for the "fixed items jitter when scrolling with the scrollbars" bug: the tree's scroll pane
 * must use {@link JViewport#SIMPLE_SCROLL_MODE}. The overview, scale and tree name are painted at
 * getVisibleRect-relative (viewport-fixed) positions, and the default BLIT_SCROLL_MODE blits them to the wrong
 * place on a scrollbar drag, flickering for a frame before a full repaint corrects them. The visual jitter
 * itself is a paint-timing artifact that cannot be asserted headlessly, so this asserts the scroll mode that
 * eliminates it. Headless-guarded (needs FlatLaf via createInstance).
 */
public final class ViewportScrollModeTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ViewportScrollMode: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { twoTip() }, conf, "n" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final JScrollPane sp = mf[ 0 ].getMainPanel().getCurrentScrollPane();
                    if ( sp == null ) {
                        System.out.println( "  [ViewportScrollModeTest] no current scroll pane" );
                        ok[ 0 ] = false;
                    }
                    else if ( sp.getViewport().getScrollMode() != JViewport.SIMPLE_SCROLL_MODE ) {
                        System.out.println( "  [ViewportScrollModeTest] tree scroll pane should use SIMPLE_SCROLL_MODE "
                                + "(to stop viewport-fixed overlays jittering on scrollbar drags), got scrollMode="
                                + sp.getViewport().getScrollMode() );
                        ok[ 0 ] = false;
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = TestFail.here();
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static Phylogeny twoTip() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "A" );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "B" );
        root.addAsChild( a );
        root.addAsChild( b );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setName( "scroll" );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private ViewportScrollModeTest() {
    }
}
