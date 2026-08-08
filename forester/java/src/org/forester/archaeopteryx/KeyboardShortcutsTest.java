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
import java.awt.Window;
import java.util.List;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.KeyboardShortcuts.Shortcut;
import org.forester.archaeopteryx.KeyboardShortcuts.ShortcutGroup;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Tests the "Help &gt; Keyboard Shortcuts" reference. The content ({@link KeyboardShortcuts#groups()} /
 * {@link KeyboardShortcuts#toHtml}) is checked purely (runs headless); headful, the dialog builds and a real
 * {@link MainFrame}'s Help menu carries the launcher item.
 */
public final class KeyboardShortcutsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "KeyboardShortcuts: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( !contentOk() ) {
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return guiOk();
    }

    private static boolean contentOk() {
        final List<ShortcutGroup> groups = KeyboardShortcuts.groups();
        if ( ( groups == null ) || groups.isEmpty() ) {
            return fail( "groups() must not be empty" );
        }
        for ( final ShortcutGroup g : groups ) {
            if ( ( g.title() == null ) || g.title().isEmpty() || ( g.shortcuts() == null ) || g.shortcuts().isEmpty() ) {
                return fail( "each group needs a title and at least one shortcut: " + g.title() );
            }
            for ( final Shortcut s : g.shortcuts() ) {
                if ( ( s.keys() == null ) || s.keys().isEmpty() || ( s.action() == null ) || s.action().isEmpty() ) {
                    return fail( "each shortcut needs non-empty keys and action (group " + g.title() + ")" );
                }
            }
        }
        final String html = KeyboardShortcuts.toHtml( groups );
        if ( !html.startsWith( "<html>" ) || !html.endsWith( "</html>" ) ) {
            return fail( "toHtml must produce a full <html>...</html> fragment" );
        }
        // a couple of representative shortcuts must be present, so the sheet is not silently empty of content
        if ( !html.contains( "Fit the whole tree to the window" )
                || !html.contains( "Copy the tree image to the clipboard" ) ) {
            return fail( "toHtml is missing expected shortcut text" );
        }
        // the clumsy keystroke-cycles for tree style (X) and label direction (D) were removed -- both settings live
        // in the Settings dialog / menu now; guard that they do not creep back into the shortcut sheet
        if ( html.contains( "Switch the tree style" ) || html.contains( "Toggle label direction" ) ) {
            return fail( "the removed D/X display keystrokes must no longer be listed" );
        }
        return true;
    }

    private static boolean guiOk() {
        try {
            // the dialog builds (with a null owner) and holds content
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final JDialog dialog = KeyboardShortcuts.buildDialog( null );
                if ( ( dialog == null ) || ( dialog.getContentPane().getComponentCount() == 0 ) ) {
                    fail( ok, "buildDialog must return a dialog with content" );
                }
                if ( dialog != null ) {
                    dialog.dispose();
                }
            } );
            // a real MainFrame's Help menu must carry the launcher, WIRED to an action, and show() must open a dialog
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, conf, "ks" ) );
            SwingUtilities.invokeAndWait( () -> {
                final JFrame frame = (JFrame) mf[ 0 ];
                try {
                    final JMenuItem item = findMenuItem( frame.getJMenuBar(), "Help", "Keyboard Shortcuts" );
                    if ( item == null ) {
                        fail( ok, "the Help menu must contain a \"Keyboard Shortcuts\" item" );
                    }
                    else {
                        // ACTIVATING the item must open the dialog -- this exercises the whole chain (menu item ->
                        // MainFrame.actionPerformed dispatch -> KeyboardShortcuts.show), so a dropped dispatch branch
                        // is caught (merely checking the item has a listener would not: the frame is always attached).
                        item.doClick();
                        boolean shown = false;
                        for ( final Window w : frame.getOwnedWindows() ) {
                            if ( ( w instanceof JDialog ) && "Keyboard Shortcuts".equals( ( (JDialog) w ).getTitle() ) ) {
                                shown = true;
                                w.dispose(); // don't leak the non-modal dialog
                            }
                        }
                        if ( !shown ) {
                            fail( ok, "activating the \"Keyboard Shortcuts\" item must open the dialog" );
                        }
                    }
                }
                finally {
                    frame.dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static JMenuItem findMenuItem( final JMenuBar bar, final String menu, final String item ) {
        if ( bar == null ) {
            return null;
        }
        for ( int i = 0; i < bar.getMenuCount(); ++i ) {
            final JMenu m = bar.getMenu( i );
            if ( ( m != null ) && menu.equals( m.getText() ) ) {
                for ( int j = 0; j < m.getItemCount(); ++j ) {
                    final JMenuItem mi = m.getItem( j );
                    if ( ( mi != null ) && item.equals( mi.getText() ) ) {
                        return mi;
                    }
                }
            }
        }
        return null;
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( new PhylogenyNode() );
        root.addAsChild( new PhylogenyNode() );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [KeyboardShortcutsTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [KeyboardShortcutsTest] " + msg );
        ok[ 0 ] = false;
    }

    private KeyboardShortcutsTest() {
    }
}
