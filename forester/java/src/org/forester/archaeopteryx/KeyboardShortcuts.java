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

import java.awt.Component;
import java.awt.Window;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import org.forester.util.ForesterUtil;

/**
 * The "Help &gt; Keyboard Shortcuts" reference. The canvas shortcuts are dispatched by hand in
 * {@link TreePanel#keyPressedCalls} (not by Swing mnemonics), so they are otherwise invisible to the user; this
 * gathers them -- plus the menu accelerators -- into one panel. Modifier symbols are platform-aware (⌘/⌥ on macOS,
 * Ctrl/Alt elsewhere). {@link #groups()} and {@link #toHtml(List)} are pure (testable headless); {@link #show} and
 * {@link #buildDialog} construct the (non-modal) dialog.
 */
final class KeyboardShortcuts {

    private static final boolean MAC = ForesterUtil.isMac();
    private static final String  CMD = MAC ? "⌘" : "Ctrl+";       // ⌘ vs Ctrl+
    private static final String  ALT = MAC ? "⌥" : "Alt+";        // ⌥ (Option) vs Alt+
    private static final String  SHIFT = MAC ? "⇧" : "Shift+";    // ⇧ vs Shift+

    record Shortcut(String keys, String action) {
    }

    record ShortcutGroup(String title, List<Shortcut> shortcuts) {
    }

    /** The grouped shortcut reference, in display order. Pure -- no toolkit needed. */
    static List<ShortcutGroup> groups() {
        final List<ShortcutGroup> groups = new ArrayList<>();
        groups.add( new ShortcutGroup( "View & navigation", List.of(
                new Shortcut( "Arrow keys", "Scroll the view" ),
                new Shortcut( "Home  /  Esc", "Fit the whole tree to the window" ),
                new Shortcut( "Page Up  /  Page Down", "Increase / decrease the font size" ) ) ) );
        groups.add( new ShortcutGroup( "Zoom", List.of(
                new Shortcut( ALT + "↑  /  " + ALT + "↓", "Zoom in / out vertically" ),
                new Shortcut( ALT + "→  /  " + ALT + "←", "Zoom in / out horizontally" ),
                new Shortcut( ALT + "+  /  " + ALT + "−", "Zoom in / out (both axes)" ),
                new Shortcut( ALT + SHIFT + "+  /  " + ALT + SHIFT + "−", "Increase / decrease the font size" ) ) ) );
        groups.add( new ShortcutGroup( "Layout", List.of(
                new Shortcut( ALT + "W", "Fit the tree width to the window" ),
                new Shortcut( ALT + "E", "Expand vertically so labels stop overlapping" ),
                new Shortcut( ALT + "O", "Order (ladderize) the subtrees" ),
                new Shortcut( ALT + "U", "Uncollapse all collapsed clades" ),
                new Shortcut( ALT + "R  /  " + ALT + SHIFT + "R", "Return to the parent tree / to the whole tree" ) ) ) );
        groups.add( new ShortcutGroup( "Circular / unrooted style only", List.of(
                new Shortcut( "S  /  A", "Rotate clockwise / counter-clockwise" ) ) ) );
        groups.add( new ShortcutGroup( "Overview (when shown)", List.of(
                new Shortcut( "O", "Move the overview to the next corner" ),
                new Shortcut( "I  /  U", "Enlarge / shrink the overview" ) ) ) );
        groups.add( new ShortcutGroup( "Menu commands", List.of(
                new Shortcut( CMD + "O", "Open a tree from a file" ),
                new Shortcut( CMD + "S", "Save the tree" ),
                new Shortcut( CMD + "W", "Close the current tab" ),
                new Shortcut( CMD + SHIFT + "C", "Copy the tree image to the clipboard" ),
                new Shortcut( CMD + "0", "Fit the tree to the window" ),
                new Shortcut( CMD + "Z  /  " + CMD + SHIFT + "Z", "Undo / redo" ),
                new Shortcut( CMD + "N", "New tree (when editing)" ) ) ) );
        return groups;
    }

    /** Renders {@link #groups()} as a self-contained HTML fragment for a Swing {@link JLabel}. Pure. */
    static String toHtml( final List<ShortcutGroup> groups ) {
        final StringBuilder sb = new StringBuilder( "<html><body style='width:360px'>" );
        for( final ShortcutGroup group : groups ) {
            sb.append( "<h3>" ).append( escape( group.title() ) ).append( "</h3>" );
            sb.append( "<table cellspacing='0' cellpadding='2'>" );
            for( final Shortcut s : group.shortcuts() ) {
                sb.append( "<tr><td valign='top'><b>" ).append( escape( s.keys() ) )
                        .append( "</b></td><td valign='top'>&nbsp;&nbsp;&nbsp;" ).append( escape( s.action() ) )
                        .append( "</td></tr>" );
            }
            sb.append( "</table>" );
        }
        return sb.append( "</body></html>" ).toString();
    }

    private static String escape( final String s ) {
        return s.replace( "&", "&amp;" ).replace( "<", "&lt;" ).replace( ">", "&gt;" );
    }

    /** Builds the (non-modal) shortcuts dialog owned by {@code owner}, sized and centered, but not yet shown. */
    static JDialog buildDialog( final Window owner ) {
        final JDialog dialog = new JDialog( owner, "Keyboard Shortcuts" );
        dialog.setDefaultCloseOperation( WindowConstants.DISPOSE_ON_CLOSE );
        final JLabel content = new JLabel( toHtml( groups() ) );
        content.setVerticalAlignment( JLabel.TOP );
        content.setBorder( BorderFactory.createEmptyBorder( 8, 14, 12, 14 ) );
        final JScrollPane scroll = new JScrollPane( content, ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
                ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER );
        scroll.setBorder( null );
        scroll.getVerticalScrollBar().setUnitIncrement( 16 );
        dialog.getContentPane().add( scroll );
        // Size to the content's natural width (so nothing is clipped horizontally -- there is no h-scrollbar) but
        // cap the height so a long list scrolls vertically instead of running off a short screen. Add the vertical
        // scrollbar's width when that scrollbar will be shown, so it does not eat into the content width.
        dialog.pack();
        final int max_height = 640;
        if ( dialog.getHeight() > max_height ) {
            dialog.setSize( dialog.getWidth() + scroll.getVerticalScrollBar().getPreferredSize().width, max_height );
        }
        dialog.setLocationRelativeTo( owner );
        return dialog;
    }

    /** Shows the shortcuts reference (non-modal) over {@code parent}. */
    static void show( final Component parent ) {
        final Window owner = ( parent instanceof Window ) ? (Window) parent
                : SwingUtilities.getWindowAncestor( parent );
        buildDialog( owner ).setVisible( true );
    }

    private KeyboardShortcuts() {
    }
}
