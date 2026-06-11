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
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Integration test for the tree-tab right-click "Close Tab" support: {@code closeTabAt(index)} must close the
 * tab at the given index (the one that was right-clicked), not merely the previously selected tab. Guarded to a
 * no-op on a headless box.
 */
public final class TabContextMenuTest {

    private static Phylogeny make( final String name ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( new PhylogenyNode() );
        root.addAsChild( new PhylogenyNode() );
        phy.setRoot( root );
        phy.setName( name );
        phy.externalNodesHaveChanged();
        return phy;
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration( null, false, false, true );
            final Phylogeny a = make( "A" );
            final Phylogeny b = make( "B" );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { a, b }, conf, "tabs" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final JTabbedPane tabs = mf[ 0 ].getMainPanel().getTabbedPane();
                if ( tabs.getTabCount() != 2 ) {
                    ok[ 0 ] = false;
                }
                // build the right-click popup for tab 0 and invoke its "Close Tab" item
                final JPopupMenu popup = ( (MainFrameApplication) mf[ 0 ] ).createTabPopupMenu( 0 );
                if ( ( popup.getComponentCount() != 1 ) || !( popup.getComponent( 0 ) instanceof JMenuItem )
                        || !"Close Tab".equals( ( (JMenuItem) popup.getComponent( 0 ) ).getText() ) ) {
                    ok[ 0 ] = false;
                }
                ( (JMenuItem) popup.getComponent( 0 ) ).doClick();
                if ( tabs.getTabCount() != 1 ) {
                    ok[ 0 ] = false;
                }
                // the clicked tab ("A") must be gone, leaving "B"
                final Phylogeny remaining = mf[ 0 ].getMainPanel().getCurrentTreePanel().getPhylogeny();
                if ( ( remaining == null ) || !"B".equals( remaining.getName() ) ) {
                    ok[ 0 ] = false;
                }
                ( (JFrame) mf[ 0 ] ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private TabContextMenuTest() {
    }
}
