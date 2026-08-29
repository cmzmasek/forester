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
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Three pieces of view state are PER TAB, not per window: the tree's orientation and the two Display-Data
 * toggles ("Show Internal Data" / "Show External Data"). They used to live on the shared {@code Options} /
 * {@code ControlPanel}, so choosing root-top for one tree silently re-oriented every other open tree, and
 * opening a demo that wants internal labels off turned them off for everything.
 * <p>
 * What is asserted here is the isolation itself -- change a value in one tab, and the OTHER tab must be
 * untouched -- plus the two halves that make the shared widgets honest: a tab switch re-seeds the controls from
 * the incoming tab, and a deliberate choice still updates the DEFAULT that a newly opened tab inherits (which is
 * what persists across restarts). Guarded to a no-op on a headless box.
 */
public final class PerTabViewStateTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "PerTabViewState: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { tree(), tree() }, new Configuration(), "per-tab" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    exercise( mf[ 0 ], ok );
                }
                finally {
                    ( ( JFrame ) mf[ 0 ] ).dispose();
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
        frame.setSize( 1000, 700 );
        ( ( JFrame ) frame ).validate();
        final ControlPanel cp = frame.getMainPanel().getControlPanel();
        final JTabbedPane tabs = frame.getMainPanel().getTabbedPane();
        if ( tabs.getTabCount() != 2 ) {
            fail( ok, "precondition: two tabs, got " + tabs.getTabCount() );
            return;
        }
        final TreePanel tp0 = frame.getMainPanel().getTreePanels().get( 0 );
        final TreePanel tp1 = frame.getMainPanel().getTreePanels().get( 1 );

        // ---- ORIENTATION is per tab -------------------------------------------------------------------
        tabs.setSelectedIndex( 0 );
        cp.getLayoutButton( LayoutIcon.Kind.ROOT_TOP ).doClick();
        if ( tp0.getTreeOrientation() != Options.TREE_ORIENTATION.ROOT_TOP ) {
            fail( ok, "tab 0 should be root-top, got " + tp0.getTreeOrientation() );
        }
        if ( tp1.getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_TOP ) {
            fail( ok, "choosing root-top in tab 0 must NOT re-orient tab 1" );
        }
        // ... and the choice still sets the default a NEW tab would inherit (what persists across restarts)
        if ( frame.getOptions().getTreeOrientation() != Options.TREE_ORIENTATION.ROOT_TOP ) {
            fail( ok, "a deliberate orientation choice should update the default for new tabs, got "
                    + frame.getOptions().getTreeOrientation() );
        }
        // switching to tab 1 must show TAB 1's orientation on the layout row, not tab 0's
        tabs.setSelectedIndex( 1 );
        if ( cp.selectedLayoutKind() == LayoutIcon.Kind.ROOT_TOP ) {
            fail( ok, "the layout row must show tab 1's own orientation after the switch, got "
                    + cp.selectedLayoutKind() );
        }
        tabs.setSelectedIndex( 0 );
        if ( cp.selectedLayoutKind() != LayoutIcon.Kind.ROOT_TOP ) {
            fail( ok, "returning to tab 0 must show its root-top again, got " + cp.selectedLayoutKind() );
        }

        // ---- the two Display-Data toggles are per tab -------------------------------------------------
        tabs.setSelectedIndex( 0 );
        cp.setShowInternalData( false );
        cp.setShowExternalData( false );
        if ( tp0.isShowInternalDataForThisTab() || tp0.isShowExternalDataForThisTab() ) {
            fail( ok, "tab 0's Display-Data toggles should both be off" );
        }
        if ( !tp1.isShowInternalDataForThisTab() || !tp1.isShowExternalDataForThisTab() ) {
            fail( ok, "turning Display Data off in tab 0 must NOT reach tab 1 -- this is what opening a demo "
                    + "used to do to every tree already open" );
        }
        // the shared checkboxes follow whichever tab is current
        tabs.setSelectedIndex( 1 );
        if ( !cp.isShowInternalData() || !cp.isShowExternalData() ) {
            fail( ok, "the Display-Data checkboxes must re-seed from tab 1 (both on)" );
        }
        tabs.setSelectedIndex( 0 );
        if ( cp.isShowInternalData() || cp.isShowExternalData() ) {
            fail( ok, "returning to tab 0 must show its own (off) Display-Data state" );
        }

        // ---- a REAL user click on the checkbox must reach the tab -------------------------------------
        // The tests above drive setShowInternalData/setCheckbox, which push to the tab by construction. A user
        // clicks the widget itself, which goes through ControlPanel.actionPerformed -- and with no branch there
        // the checkbox toggled, the tree repainted, and NOTHING changed. Drive the actual widget.
        tabs.setSelectedIndex( 0 );
        final javax.swing.AbstractButton internal_cb = findButton( cp, "Show Internal Data" );
        final javax.swing.AbstractButton external_cb = findButton( cp, "Show External Data" );
        if ( ( internal_cb == null ) || ( external_cb == null ) ) {
            fail( ok, "could not find the Display-Data checkboxes on the control panel" );
        }
        else {
            internal_cb.doClick(); // was off (set above) -> on
            if ( !tp0.isShowInternalDataForThisTab() ) {
                fail( ok, "clicking 'Show Internal Data' must reach the tab -- the checkbox must not be inert" );
            }
            external_cb.doClick();
            if ( !tp0.isShowExternalDataForThisTab() ) {
                fail( ok, "clicking 'Show External Data' must reach the tab" );
            }
            // ... and it must SURVIVE a tab round trip, rather than being undone by the re-seed
            tabs.setSelectedIndex( 1 );
            tabs.setSelectedIndex( 0 );
            if ( !tp0.isShowInternalDataForThisTab() || !cp.isShowInternalData() ) {
                fail( ok, "a clicked Display-Data change must survive switching away and back" );
            }
            internal_cb.doClick(); // back off, so the reset check below has something to restore
        }

        // ---- Reset to Defaults puts EVERY tab back ----------------------------------------------------
        frame.resetToDefaults();
        for ( final TreePanel tp : new TreePanel[] { tp0, tp1 } ) {
            if ( tp.getTreeOrientation() != Options.TREE_ORIENTATION.ROOT_LEFT ) {
                fail( ok, "reset must return every tab to root-left, got " + tp.getTreeOrientation() );
            }
            if ( !tp.isShowInternalDataForThisTab() || !tp.isShowExternalDataForThisTab() ) {
                fail( ok, "reset must turn both Display-Data toggles back on in every tab" );
            }
        }
        // ... and the SHARED checkboxes must say so too. A stale "off" here is not just cosmetic: the next push
        // writes both values from the widget, which would put the stale off straight back onto the tab.
        if ( !cp.isShowInternalData() || !cp.isShowExternalData() ) {
            fail( ok, "reset must re-seed the Display-Data checkboxes, not just the per-tab flags" );
        }

        // ---- a NEW tab inherits the current Display-Data state, as orientation does --------------------
        cp.setShowExternalData( false );
        if ( frame._new_item == null ) {
            fail( ok, "could not reach File -> New" );
            return;
        }
        frame._new_item.doClick();
        final TreePanel fresh = frame.getCurrentTreePanel();
        if ( ( fresh == null ) || fresh.isShowExternalDataForThisTab() ) {
            fail( ok, "a new tab should open with the Display-Data state you are currently looking at" );
        }
    }

    /** Depth-first search for the first AbstractButton with the given text, or null. */
    private static javax.swing.AbstractButton findButton( final java.awt.Container c, final String text ) {
        for ( final java.awt.Component comp : c.getComponents() ) {
            if ( ( comp instanceof javax.swing.AbstractButton )
                    && text.equals( ( ( javax.swing.AbstractButton ) comp ).getText() ) ) {
                return ( javax.swing.AbstractButton ) comp;
            }
            if ( comp instanceof java.awt.Container ) {
                final javax.swing.AbstractButton found = findButton( ( java.awt.Container ) comp, text );
                if ( found != null ) {
                    return found;
                }
            }
        }
        return null;
    }

    /** A small tree with branch lengths, so every layout and both Display-Data toggles are meaningful. */
    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode inner = new PhylogenyNode();
        inner.setDistanceToParent( 0.2 );
        for ( final String name : new String[] { "a", "b" } ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( name );
            tip.setDistanceToParent( 0.3 );
            inner.addAsChild( tip );
        }
        final PhylogenyNode c = new PhylogenyNode();
        c.setName( "c" );
        c.setDistanceToParent( 0.7 );
        root.addAsChild( inner );
        root.addAsChild( c );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        return phy;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [PerTabViewStateTest] " + msg );
        ok[ 0 ] = false;
    }

    private PerTabViewStateTest() {
    }
}
