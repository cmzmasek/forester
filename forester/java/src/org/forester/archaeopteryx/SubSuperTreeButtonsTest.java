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
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.util.ArrayDeque;
import java.util.Deque;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Integration test for the two sub-tree navigation buttons in the control panel: <b>R</b> must
 * jump all the way back to the complete tree, while <b>R1</b> must move up by exactly one level
 * towards it. Both are exercised through {@link ControlPanel#returnedToWholeTreePressed()} and
 * {@link ControlPanel#returnedToSuperTreePressed()}; depth is checked via the number of external
 * nodes of the displayed (sub)tree. Guarded to a no-op on a headless box. Needs FlatLaf on the
 * classpath (via {@code MainFrameApplication.createInstance}), so it is run standalone, not as
 * part of the headless suite.
 */
public final class SubSuperTreeButtonsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SubSuperTreeButtons: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration( null, false, false, true );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { nestedTree() }, conf,
                                                                         "subsuper" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainPanel mp = mf[ 0 ].getMainPanel();
                final TreePanel tp = mp.getCurrentTreePanel();
                final ControlPanel cp = mp.getControlPanel();
                // the complete tree has 4 leaves and is not a sub-tree
                if ( ( leaves( tp ) != 4 ) || tp.isCurrentTreeIsSubtree() ) {
                    ok[ 0 ] = false;
                }
                // descend into A (3 leaves) and then the nested B (2 leaves): now two levels deep
                tp.subTree( internalNamed( tp.getPhylogeny(), "A" ) );
                if ( leaves( tp ) != 3 ) {
                    ok[ 0 ] = false;
                }
                tp.subTree( internalNamed( tp.getPhylogeny(), "B" ) );
                if ( ( leaves( tp ) != 2 ) || !tp.isCurrentTreeIsSubtree() ) {
                    ok[ 0 ] = false;
                }
                // R1 moves up exactly one level -> back at A (3 leaves), still a sub-tree
                cp.returnedToSuperTreePressed();
                if ( ( leaves( tp ) != 3 ) || !tp.isCurrentTreeIsSubtree() ) {
                    ok[ 0 ] = false;
                }
                // R1 again -> back at the complete tree (4 leaves), no longer a sub-tree
                cp.returnedToSuperTreePressed();
                if ( ( leaves( tp ) != 4 ) || tp.isCurrentTreeIsSubtree() ) {
                    ok[ 0 ] = false;
                }
                // descend two levels again, then R jumps straight to the root in a single press
                tp.subTree( internalNamed( tp.getPhylogeny(), "A" ) );
                tp.subTree( internalNamed( tp.getPhylogeny(), "B" ) );
                if ( leaves( tp ) != 2 ) {
                    ok[ 0 ] = false;
                }
                cp.returnedToWholeTreePressed();
                if ( ( leaves( tp ) != 4 ) || tp.isCurrentTreeIsSubtree() ) {
                    ok[ 0 ] = false;
                }
                // R on the complete tree is a harmless no-op
                cp.returnedToWholeTreePressed();
                if ( ( leaves( tp ) != 4 ) || tp.isCurrentTreeIsSubtree() ) {
                    ok[ 0 ] = false;
                }
                // keyboard shortcuts (dispatched through the panel's own key listener):
                // Alt+R moves up one level, Alt+Shift+R returns all the way to the root.
                tp.subTree( internalNamed( tp.getPhylogeny(), "A" ) );
                tp.subTree( internalNamed( tp.getPhylogeny(), "B" ) );
                pressR( tp, false ); // Alt+R: up one -> A (3 leaves)
                if ( ( leaves( tp ) != 3 ) || !tp.isCurrentTreeIsSubtree() ) {
                    ok[ 0 ] = false;
                }
                pressR( tp, false ); // Alt+R again: up one -> whole tree (4 leaves)
                if ( ( leaves( tp ) != 4 ) || tp.isCurrentTreeIsSubtree() ) {
                    ok[ 0 ] = false;
                }
                tp.subTree( internalNamed( tp.getPhylogeny(), "A" ) );
                tp.subTree( internalNamed( tp.getPhylogeny(), "B" ) );
                pressR( tp, true ); // Alt+Shift+R: straight to the root (4 leaves)
                if ( ( leaves( tp ) != 4 ) || tp.isCurrentTreeIsSubtree() ) {
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

    private static int leaves( final TreePanel tp ) {
        return tp.getPhylogeny().getNumberOfExternalNodes();
    }

    /** Press Alt+R (or Alt+Shift+R when {@code shift}) on the panel via its registered key listener. */
    private static void pressR( final TreePanel tp, final boolean shift ) {
        int mods = InputEvent.ALT_DOWN_MASK;
        if ( shift ) {
            mods |= InputEvent.SHIFT_DOWN_MASK;
        }
        tp.dispatchEvent( new KeyEvent( tp, KeyEvent.KEY_PRESSED, System.currentTimeMillis(), mods, KeyEvent.VK_R,
                                        'r' ) );
    }

    /**
     *      root
     *     /    \
     *    A      r2          leaves: b1, b2, a3, r2  (whole = 4)
     *   / \                 sub-tree A: b1, b2, a3  (= 3)
     *  B   a3               sub-tree B: b1, b2      (= 2)
     * / \
     * b1 b2
     */
    private static Phylogeny nestedTree() {
        final PhylogenyNode root = named( "root" );
        final PhylogenyNode a = named( "A" );
        final PhylogenyNode b = named( "B" );
        b.addAsChild( named( "b1" ) );
        b.addAsChild( named( "b2" ) );
        a.addAsChild( b );
        a.addAsChild( named( "a3" ) );
        root.addAsChild( a );
        root.addAsChild( named( "r2" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode named( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    /** The first internal (non-external) node with the given name in the (sub)tree, or null. */
    private static PhylogenyNode internalNamed( final Phylogeny phy, final String name ) {
        final Deque<PhylogenyNode> stack = new ArrayDeque<PhylogenyNode>();
        stack.push( phy.getRoot() );
        while ( !stack.isEmpty() ) {
            final PhylogenyNode n = stack.pop();
            if ( !n.isExternal() && name.equals( n.getName() ) ) {
                return n;
            }
            for( final PhylogenyNode c : n.getDescendants() ) {
                stack.push( c );
            }
        }
        return null;
    }

    private SubSuperTreeButtonsTest() {
    }
}
