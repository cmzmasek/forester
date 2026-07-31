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
import java.awt.Toolkit;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import javax.swing.JFrame;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Verifies the File-menu keyboard accelerators: Open, Save, Close Tab, Copy Image, and (when editable) New carry
 * the expected platform-shortcut keystrokes. Copy Image must be Cmd/Ctrl+SHIFT+C -- deliberately NOT plain
 * Cmd/Ctrl+C, so it cannot hijack text copy. Also verifies the focused-canvas key listener does not CONSUME a
 * menu-shortcut key (which would shadow the accelerator) while still handling plain keys. Guarded to a no-op when
 * headless (needs FlatLaf via createInstance).
 */
public final class MenuAcceleratorTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MenuAccelerator: " + ( ok ? "OK." : "FAILED." ) );
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
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, conf, "acc" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final int sc = Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx();
                    expect( ok, "Read Tree from File", frame._open_item, KeyStroke.getKeyStroke( KeyEvent.VK_O, sc ) );
                    expect( ok, "Save Tree As", frame._save_item, KeyStroke.getKeyStroke( KeyEvent.VK_S, sc ) );
                    expect( ok, "Close Tab", frame._close_item, KeyStroke.getKeyStroke( KeyEvent.VK_W, sc ) );
                    expect( ok, "Copy Image", frame._copy_image_to_clipboard_item,
                            KeyStroke.getKeyStroke( KeyEvent.VK_C, sc | InputEvent.SHIFT_DOWN_MASK ) );
                    // Copy Image must NOT be plain shortcut+C (that would shadow text copy)
                    if ( KeyStroke.getKeyStroke( KeyEvent.VK_C, sc )
                            .equals( frame._copy_image_to_clipboard_item.getAccelerator() ) ) {
                        fail( ok, "Copy Image must not use plain Cmd/Ctrl+C" );
                    }
                    // New exists only in an editable configuration
                    if ( frame._new_item != null ) {
                        expect( ok, "New", frame._new_item, KeyStroke.getKeyStroke( KeyEvent.VK_N, sc ) );
                    }

                    // The focused-canvas key listener must NOT consume a menu-shortcut (Cmd/Ctrl) key -- otherwise
                    // it would shadow the accelerator entirely (a KeyListener runs before WHEN_IN_FOCUSED_WINDOW
                    // menu accelerators). A plain letter must still be handled/consumed as before.
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final long when = System.currentTimeMillis();
                    final KeyEvent with_shortcut = new KeyEvent( tp, KeyEvent.KEY_PRESSED, when, sc, KeyEvent.VK_W,
                                                                 'W' );
                    for( final KeyListener kl : tp.getKeyListeners() ) {
                        kl.keyPressed( with_shortcut );
                    }
                    if ( with_shortcut.isConsumed() ) {
                        fail( ok, "the canvas key listener must not consume a menu-shortcut key (it would shadow the "
                                + "menu accelerator)" );
                    }
                    final KeyEvent plain = new KeyEvent( tp, KeyEvent.KEY_PRESSED, when, 0, KeyEvent.VK_B, 'B' );
                    for( final KeyListener kl : tp.getKeyListeners() ) {
                        kl.keyPressed( plain );
                    }
                    if ( !plain.isConsumed() ) {
                        fail( ok, "the canvas key listener must still consume a plain (non-shortcut) key" );
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

    private static void expect( final boolean[] ok, final String name, final javax.swing.JMenuItem item,
                                final KeyStroke expected ) {
        if ( item == null ) {
            fail( ok, name + " menu item is missing" );
        }
        else if ( !expected.equals( item.getAccelerator() ) ) {
            fail( ok, name + " accelerator should be " + expected + ", got " + item.getAccelerator() );
        }
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [MenuAcceleratorTest] " + msg );
        ok[ 0 ] = false;
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 3; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private MenuAcceleratorTest() {
    }
}
