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

import java.awt.Color;
import java.awt.GraphicsEnvironment;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * When a background task is registered, the menu bar's "processes running" menu gains the animated
 * {@link EqualizerIcon} activity indicator (in the theme accent, not red) and its animation timer runs; when the
 * last task finishes the indicator is removed and the timer stops. Needs a display (via
 * {@code MainFrameApplication.createInstance}), so it is a no-op (returns true) when headless.
 */
public final class ProcessMenuAnimationTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ProcessMenuAnimation: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, new Configuration(),
                                                                         "proc" ) );
            final long[] id = { -1 };
            // register a task and refresh the process menu (updateProcessMenu posts doUpdateProcessMenu to the EDT)
            SwingUtilities.invokeAndWait( () -> {
                id[ 0 ] = mf[ 0 ].getProcessPool().addProcess( "probe-task" );
                mf[ 0 ].updateProcessMenu();
            } );
            SwingUtilities.invokeAndWait( () -> { /* flush the queued doUpdateProcessMenu */ } );
            SwingUtilities.invokeAndWait( () -> {
                if ( ( mf[ 0 ]._process_menu == null ) || !( mf[ 0 ]._process_menu.getIcon() instanceof EqualizerIcon ) ) {
                    ok[ 0 ] = fail( "a running task must show the equalizer activity icon on the process menu" );
                }
                else if ( !mf[ 0 ].isProcessAnimationRunningForTest() ) {
                    ok[ 0 ] = fail( "the activity animation timer must run while a task is running" );
                }
                else {
                    final Color fg = mf[ 0 ]._process_menu.getForeground();
                    if ( ( fg == null ) || Color.RED.equals( fg ) ) {
                        ok[ 0 ] = fail( "the process menu must use the theme accent color, not red" );
                    }
                }
            } );
            // finish the task -> the indicator + its timer go away
            SwingUtilities.invokeAndWait( () -> {
                mf[ 0 ].getProcessPool().removeProcess( id[ 0 ] );
                mf[ 0 ].updateProcessMenu();
            } );
            SwingUtilities.invokeAndWait( () -> { /* flush */ } );
            SwingUtilities.invokeAndWait( () -> {
                if ( mf[ 0 ].isProcessAnimationRunningForTest() ) {
                    ok[ 0 ] = fail( "the animation timer must stop when no task is running" );
                }
                ( (JFrame) mf[ 0 ] ).dispose();
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ProcessMenuAnimationTest] " + msg );
        return false;
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        for ( int i = 0; i < 3; i++ ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "t" + i );
            root.addAsChild( n );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private ProcessMenuAnimationTest() {
    }
}
