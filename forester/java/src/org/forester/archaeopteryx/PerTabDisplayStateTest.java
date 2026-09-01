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

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Every Display option now has a PER-TAB value, not just the two that were migrated earlier.
 * <p>
 * The checkboxes are shared by the window and show whichever tab is in front; the tab's own copy is what a figure
 * is made of, and what a saved figure has to restore. Two things must hold for that to work: reading an option and
 * writing it must agree for ALL of them (they are two directions over one widget mapping, and two switches would
 * eventually disagree), and a tab must be able to carry a value the widgets do not currently show.
 */
public final class PerTabDisplayStateTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "PerTabDisplayState: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { tree() }, new Configuration(), "pertabdisp" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();

                // (1) SYMMETRY over every option -- the guard against the read and the write drifting apart
                for( final DisplayOption opt : DisplayOption.values() ) {
                    cp.setCheckbox( opt, true );
                    if ( !cp.isCheckboxSelected( opt ) ) {
                        ok[ 0 ] = fail( opt + ": set to true but reads back false" );
                        return;
                    }
                    cp.setCheckbox( opt, false );
                    if ( cp.isCheckboxSelected( opt ) ) {
                        ok[ 0 ] = fail( opt + ": set to false but reads back true" );
                        return;
                    }
                }

                // (2) the whole set is copied onto the tab, not just the two Display-Data toggles
                for( final DisplayOption opt : DisplayOption.values() ) {
                    cp.setCheckbox( opt, true );
                }
                cp.pushDisplayDataToCurrentTab();
                for( final DisplayOption opt : DisplayOption.values() ) {
                    if ( opt == DisplayOption.DISPLAY_AS_PHYLOGRAM ) {
                        // a THREE-way radio choice, not a toggle: forcing it through a boolean turned "aligned
                        // phylogram" into "cladogram" on the way back (LayoutButtonsTest caught exactly that),
                        // so it is excluded from the boolean per-tab state and keeps its own representation
                        if ( tp.hasOwnDisplayState( opt ) ) {
                            ok[ 0 ] = fail( "DISPLAY_AS_PHYLOGRAM must NOT be stored as a per-tab boolean" );
                            return;
                        }
                        continue;
                    }
                    if ( !tp.hasOwnDisplayState( opt ) ) {
                        ok[ 0 ] = fail( opt + " was not copied onto the tab -- the tab must own the whole figure" );
                        return;
                    }
                    if ( !tp.shows( opt ) ) {
                        ok[ 0 ] = fail( opt + " should be on for this tab" );
                        return;
                    }
                }

                // (3) a tab can hold a value the shared widgets are NOT currently showing -- which is the entire
                // point: a second tab, or a restored figure, must not be overwritten by the front tab's widgets
                // (setCheckbox deliberately pushes to the current tab, so set the WIDGET first and the tab after)
                cp.setCheckbox( DisplayOption.SHOW_TAX_RANK, true );
                tp.setShows( DisplayOption.SHOW_TAX_RANK, false );
                if ( tp.shows( DisplayOption.SHOW_TAX_RANK ) ) {
                    ok[ 0 ] = fail( "the tab must be able to hold a value the shared widget is not showing" );
                    return;
                }

                // (4) ...and re-seeding puts the tab's value back into the widgets
                cp.reseedDisplayDataFromCurrentTab();
                if ( cp.isCheckboxSelected( DisplayOption.SHOW_TAX_RANK ) ) {
                    ok[ 0 ] = fail( "re-seeding must show the TAB's value in the shared checkbox" );
                }
            } );
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [PerTabDisplayStateTest] " + msg );
        return false;
    }

    private static Phylogeny tree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 4; ++i ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "t" + i );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private PerTabDisplayStateTest() {
    }
}
