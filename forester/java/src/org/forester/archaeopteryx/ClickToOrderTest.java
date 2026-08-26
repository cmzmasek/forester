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
import java.util.List;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Verifies the order of the "Click on Node to:" actions -- which feed BOTH the control-panel dropdown AND the
 * node right-click popup from one list ({@link ControlPanel#getSingleClickToNames()}). In particular "Root/Reroot"
 * must come before "Edit Node Data". Guarded to a no-op on a headless box (needs FlatLaf via {@code createInstance}).
 */
public final class ClickToOrderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ClickToOrder: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI test; needs a display toolkit
        }
        try {
            final PhylogenyNode root = new PhylogenyNode();
            root.addAsChild( new PhylogenyNode() );
            root.addAsChild( new PhylogenyNode() );
            final Phylogeny phy = new Phylogeny();
            phy.setRoot( root );
            phy.setRooted( true );
            phy.externalNodesHaveChanged();

            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "clickto" ) );
            final List<String> names = mf[ 0 ].getMainPanel().getControlPanel().getSingleClickToNames();
            final int reroot = names.indexOf( ClickToOption.REROOT.title() );
            final int edit = names.indexOf( ClickToOption.EDIT_NODE_DATA.title() );

            boolean ok = true;
            if ( ( reroot < 0 ) || ( edit < 0 ) ) {
                ok = fail( "both 'Root/Reroot' and 'Edit Node Data' must be offered (reroot=" + reroot + ", edit="
                        + edit + ")" );
            }
            else if ( reroot > edit ) {
                ok = fail( "'Root/Reroot' must come before 'Edit Node Data' (reroot=" + reroot + ", edit=" + edit
                        + ")" );
            }
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            return ok;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [ClickToOrderTest] " + message );
        return false;
    }

    private ClickToOrderTest() {
        // not instantiable
    }
}
