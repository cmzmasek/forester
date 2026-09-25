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
import java.util.function.BiConsumer;

import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;

import com.formdev.flatlaf.FlatClientProperties;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The "x" on a tree tab. FlatLaf draws and hit-tests the button (a Swing {@code JTabbedPane} has none of its own)
 * and hands the callback the tab's index; the callback must route to {@code closeTabAt}, which is the same path as
 * File &gt; Close Tab and the tab right-click menu -- so the unsaved-changes confirmation and the per-tab teardown
 * stay in one place rather than being reimplemented behind the button.
 * <p>
 * The check that matters is that it closes the tab whose "x" was pressed, NOT the selected one: closeTabAt selects
 * first and then closes, and a button that skipped that step would quietly close the wrong tree.
 */
public final class TabCloseButtonTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TabCloseButton: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication.createInstance(
                    new Phylogeny[] { tree( "alpha" ), tree( "beta" ), tree( "gamma" ) }, new Configuration(),
                    "tabclose" ) );
            final JTabbedPane tabs = mf[ 0 ].getMainPanel().getTabbedPane();
            if ( tabs.getTabCount() != 3 ) {
                fail( ok, "precondition: three tabs, got " + tabs.getTabCount() );
                dispose( mf );
                return ok[ 0 ];
            }
            // (1) the tabs really are closable, and carry a callback
            if ( !Boolean.TRUE.equals( tabs.getClientProperty( FlatClientProperties.TABBED_PANE_TAB_CLOSABLE ) ) ) {
                fail( ok, "the tabbed pane must declare its tabs closable, or no 'x' is drawn at all" );
            }
            final Object cb = tabs.getClientProperty( FlatClientProperties.TABBED_PANE_TAB_CLOSE_CALLBACK );
            if ( !( cb instanceof BiConsumer ) ) {
                fail( ok, "there must be a close callback; got " + cb );
                dispose( mf );
                return ok[ 0 ];
            }
            // (2) pressing the "x" of a tab that is NOT selected closes THAT tab. Select the last one first, so a
            // callback that simply closed the current pane would remove "gamma" and pass an index-blind check.
            SwingUtilities.invokeAndWait( () -> tabs.setSelectedIndex( 2 ) );
            @SuppressWarnings( "unchecked" )
            final BiConsumer<JTabbedPane, Integer> close = (BiConsumer<JTabbedPane, Integer>) cb;
            SwingUtilities.invokeAndWait( () -> close.accept( tabs, Integer.valueOf( 1 ) ) ); // the MIDDLE tab
            if ( tabs.getTabCount() != 2 ) {
                fail( ok, "closing one tab must leave two, got " + tabs.getTabCount() );
            }
            if ( titles( tabs ).contains( "beta" ) ) {
                fail( ok, "the tab whose 'x' was pressed must be the one closed; still have: " + titles( tabs ) );
            }
            if ( !titles( tabs ).contains( "gamma" ) ) {
                fail( ok, "the SELECTED tab must not be closed instead of the pressed one; have: "
                        + titles( tabs ) );
            }
            if ( !titles( tabs ).contains( "alpha" ) ) {
                fail( ok, "the untouched tab must survive; have: " + titles( tabs ) );
            }
            // (3) an out-of-range index must be ignored rather than close something arbitrary -- the button is
            // drawn per tab and an index can go stale between the draw and the click
            SwingUtilities.invokeAndWait( () -> close.accept( tabs, Integer.valueOf( 99 ) ) );
            if ( tabs.getTabCount() != 2 ) {
                fail( ok, "a stale/out-of-range index must close nothing, got " + tabs.getTabCount() + " tabs" );
            }
            dispose( mf );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
        return ok[ 0 ];
    }

    private static String titles( final JTabbedPane t ) {
        final StringBuilder sb = new StringBuilder();
        for( int i = 0; i < t.getTabCount(); ++i ) {
            sb.append( i > 0 ? ", " : "" ).append( t.getTitleAt( i ) );
        }
        return sb.toString();
    }

    private static Phylogeny tree( final String name ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( final String s : new String[] { name + "_a", name + "_b" } ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( s );
            tip.setDistanceToParent( 0.2 );
            root.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        phy.setName( name );
        return phy;
    }

    private static void dispose( final MainFrame[] mf ) throws Exception {
        SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [TabCloseButtonTest] " + message );
        ok[ 0 ] = false;
    }

    private TabCloseButtonTest() {
        // not instantiable
    }
}
