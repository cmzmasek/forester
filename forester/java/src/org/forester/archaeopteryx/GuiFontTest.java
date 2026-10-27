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
import java.awt.Container;
import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.util.ArrayList;
import java.util.List;

import javax.swing.AbstractButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import com.formdev.flatlaf.util.UIScale;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * One GUI font (Christian, 2026-10-27: "the same font for all menus and control panels", at "12"): the
 * look-and-feel's own face at 12 points is FlatLaf's default font ({@link MainFrame#applyGuiFont}), and every
 * menu, menu item, control-panel widget and the node menu take theirs from it -- the panel's headers and the
 * node menu's title in bold, nothing else. Pins the LITERAL 12 (not {@code Configuration.getGuiFontSize()}, which
 * the font is set from -- the two must not drift together; at FlatLaf's user scale factor of 1), that a theme
 * switch keeps the font (it re-installs the look-and-feel, which replaces the look-and-feel's defaults; a
 * developer default outlives that AND the font is re-applied after every setup -- each safeguard is killed on
 * its own), and the glyph size that follows it. Needs a display:
 * a no-op when headless.
 */
public final class GuiFontTest {

    private static final int SIZE = 12;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "GuiFont: " + ( ok ? "OK." : "FAILED." ) );
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
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, conf, "gui font" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    ok[ 0 ] = liveOk( mf[ 0 ] );
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
            SwingUtilities.invokeAndWait( () -> ok[ 0 ] &= themeSwitchKeepsItOk( conf ) );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    /** The window as built: the look-and-feel's font everywhere, bold only where a header or a title is. */
    private static boolean liveOk( final MainFrame mf ) {
        boolean ok = true;
        final Font base = UIManager.getFont( "defaultFont" );
        ok &= ck( base != null, "the look-and-feel has a default font" );
        if ( base == null ) {
            return false;
        }
        // the literal, at FlatLaf's user scale factor of 1 (the suite's; a flatlaf.uiScale property or a Windows
        // text-size setting scales it, and the font must follow -- proven standalone with -Dflatlaf.uiScale=1.5)
        final int expected = Math.round( UIScale.scale( (float) SIZE ) );
        ok &= ck( ( UIScale.getUserScaleFactor() != 1f ) || ( expected == SIZE ), "precondition: scale 1 means " + SIZE );
        ok &= ck( base.getSize() == expected, "the GUI font is " + expected + " points (" + SIZE + " at scale "
                + UIScale.getUserScaleFactor() + "), got " + base.getSize() );
        ok &= ck( !base.isBold() && !base.isItalic(), "the GUI font is plain: " + base );
        for( final String key : new String[] { "Label.font", "Menu.font", "MenuItem.font", "CheckBoxMenuItem.font",
                "RadioButtonMenuItem.font", "PopupMenu.font", "CheckBox.font", "RadioButton.font", "Button.font",
                "ToggleButton.font", "ComboBox.font", "TextField.font", "ToolTip.font", "TabbedPane.font",
                "TitledBorder.font" } ) {
            ok &= ck( sameFace( base, UIManager.getFont( key ) ), key + " is the GUI font: " + UIManager.getFont( key ) );
        }
        // ---- the control panel: every widget with text, headers bold ------------------------------------
        final ControlPanel cp = mf.getMainPanel().getControlPanel();
        final List<Component> all = new ArrayList<>();
        collect( cp, all );
        int checked = 0, bold = 0;
        for( final Component c : all ) {
            if ( ( c instanceof JLabel ) || ( c instanceof AbstractButton ) || ( c instanceof JComboBox )
                    || ( c instanceof JTextField ) ) {
                checked++;
                ok &= ck( sameFace( base, c.getFont() ),
                          "control '" + textOf( c ) + "' (" + c.getClass().getSimpleName() + ") has the GUI font, got "
                                  + c.getFont() );
                if ( ( c.getFont() != null ) && c.getFont().isBold() ) {
                    bold++;
                }
            }
        }
        ok &= ck( checked >= 60, "precondition: the control panel was walked (" + checked + " widgets)" );
        for( final String header : new String[] { "Display Data:", "Search Options:", "Search (A):", "Search (B):" } ) {
            final JLabel label = label( all, header );
            ok &= ck( label != null, "precondition: the header '" + header + "' is there" );
            ok &= ck( ( label != null ) && label.getFont().isBold(), "the header '" + header + "' is bold" );
        }
        final AbstractButton checkbox = firstCheckbox( all );
        ok &= ck( checkbox != null, "precondition: the panel has a checkbox" );
        ok &= ck( ( checkbox != null ) && !checkbox.getFont().isBold(),
                  "a checkbox ('" + textOf( checkbox ) + "') is not bold" );
        ok &= ck( bold < checked / 2, "bold is the exception on the panel: " + bold + " of " + checked );
        // ---- the menu bar: every menu and every item, submenus included --------------------------------
        final JMenuBar bar = ( (JFrame) mf ).getJMenuBar();
        final int[] counts = { 0, 0 }; // items, submenus
        for( int i = 0; i < bar.getMenuCount(); i++ ) {
            final JMenu menu = bar.getMenu( i );
            if ( menu == null ) {
                continue;
            }
            ok &= ck( sameFace( base, menu.getFont() ), "the menu '" + menu.getText() + "' has the GUI font: "
                    + menu.getFont() );
            ok &= menuOk( base, menu, counts );
        }
        ok &= ck( counts[ 0 ] >= 80, "precondition: the menu bar was walked (" + counts[ 0 ] + " items)" );
        ok &= ck( counts[ 1 ] >= 3, "precondition: submenus were walked (" + counts[ 1 ] + ")" );
        // ---- the node menu: a bold title, plain entries ----------------------------------------------
        final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
        final JPopupMenu node_menu = tp.buildNodeMenu( tp.getPhylogeny().getRoot() );
        final Component title = node_menu.getComponent( 0 );
        ok &= ck( ( title instanceof JLabel ) && sameFace( base, title.getFont() ) && title.getFont().isBold(),
                  "the node menu's title is the GUI font in bold: " + title.getFont() );
        int entries = 0;
        for( final Component c : node_menu.getComponents() ) {
            if ( c instanceof JMenuItem ) {
                entries++;
                ok &= ck( sameFace( base, c.getFont() ) && !c.getFont().isBold(),
                          "the node menu entry '" + ( (JMenuItem) c ).getText() + "' is the GUI font, plain: " + c.getFont() );
            }
        }
        ok &= ck( entries >= 15, "precondition: the node menu was walked (" + entries + " entries)" );
        // ---- what follows the size ------------------------------------------------------------------
        ok &= ck( ControlPanel.glyphIconSize() == 18, "the layout glyphs are 18 px at this font, got "
                + ControlPanel.glyphIconSize() );
        return ok;
    }

    /**
     * A theme switch re-installs the look-and-feel, which replaces the look-and-feel's defaults. The GUI font
     * survives that two ways, and each is killed on its own here: (1) a BARE {@code setup()} with nothing
     * re-applied -- only a developer default ({@code UIManager.put}) carries the font across it, one put into the
     * look-and-feel defaults is gone; (2) the developer default CLEARED, then {@code installLookAndFeel} -- only
     * the re-apply after the setup restores it. Then the switch as the application does it, both ways round.
     */
    private static boolean themeSwitchKeepsItOk( final Configuration conf ) {
        boolean ok = true;
        final Configuration.UI original = conf.getUi();
        final int expected = Math.round( UIScale.scale( (float) SIZE ) );
        // (1) the developer default alone
        com.formdev.flatlaf.FlatDarkLaf.setup();
        ok &= ck( sizeIs( expected ), "after a bare dark setup() the developer default carries the GUI font: "
                + UIManager.getFont( "defaultFont" ) );
        // (2) the re-apply alone
        UIManager.put( "defaultFont", null );
        ok &= ck( !sizeIs( expected ), "precondition: the developer default is cleared: " + UIManager.getFont( "defaultFont" ) );
        MainFrame.installLookAndFeel( Configuration.UI.FLAT_LIGHT );
        ok &= ck( sizeIs( expected ), "after installLookAndFeel with no developer default the re-apply sets the GUI font: "
                + UIManager.getFont( "defaultFont" ) );
        // the switch as the application makes it
        for( final Configuration.UI ui : new Configuration.UI[] { Configuration.UI.FLAT_DARK, Configuration.UI.FLAT_LIGHT } ) {
            MainFrame.installLookAndFeel( ui );
            ok &= ck( sizeIs( expected ), "after installing " + ui + " the GUI font is " + expected + ": "
                    + UIManager.getFont( "defaultFont" ) );
            ok &= ck( sameFace( UIManager.getFont( "defaultFont" ), UIManager.getFont( "MenuItem.font" ) ),
                      "after installing " + ui + " the menu items follow: " + UIManager.getFont( "MenuItem.font" ) );
        }
        MainFrame.installLookAndFeel( original );
        return ok;
    }

    private static boolean sizeIs( final int size ) {
        final Font f = UIManager.getFont( "defaultFont" );
        return ( f != null ) && ( f.getSize() == size );
    }

    private static boolean menuOk( final Font base, final JMenu menu, final int[] counts ) {
        boolean ok = true;
        for( final Component c : menu.getMenuComponents() ) {
            if ( c instanceof JMenuItem ) {
                counts[ 0 ]++;
                ok &= ck( sameFace( base, c.getFont() ), "the item '" + ( (JMenuItem) c ).getText() + "' under '"
                        + menu.getText() + "' has the GUI font: " + c.getFont() );
                if ( c instanceof JMenu ) {
                    counts[ 1 ]++;
                    ok &= menuOk( base, (JMenu) c, counts );
                }
            }
        }
        return ok;
    }

    /** Same family and size; style is the caller's business. */
    static boolean sameFace( final Font base, final Font f ) {
        return ( f != null ) && base.getFamily().equals( f.getFamily() ) && ( base.getSize() == f.getSize() );
    }

    private static void collect( final Container c, final List<Component> out ) {
        for( final Component k : c.getComponents() ) {
            out.add( k );
            if ( k instanceof Container ) {
                collect( (Container) k, out );
            }
        }
    }

    private static JLabel label( final List<Component> all, final String text ) {
        for( final Component c : all ) {
            if ( ( c instanceof JLabel ) && text.equals( ( (JLabel) c ).getText() ) ) {
                return (JLabel) c;
            }
        }
        return null;
    }

    private static AbstractButton firstCheckbox( final List<Component> all ) {
        for( final Component c : all ) {
            if ( c instanceof javax.swing.JCheckBox ) {
                return (AbstractButton) c;
            }
        }
        return null;
    }

    private static String textOf( final Component c ) {
        if ( c instanceof JLabel ) {
            return ( (JLabel) c ).getText();
        }
        if ( c instanceof AbstractButton ) {
            return ( (AbstractButton) c ).getText();
        }
        return c.getName() == null ? "" : c.getName();
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode ab = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "a" );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "b" );
        final PhylogenyNode c = new PhylogenyNode();
        c.setName( "c" );
        ab.addAsChild( a );
        ab.addAsChild( b );
        root.addAsChild( ab );
        root.addAsChild( c );
        final Phylogeny p = new Phylogeny();
        p.setRoot( root );
        p.setRooted( true );
        return p;
    }

    private static boolean ck( final boolean cond, final String what ) {
        if ( !cond ) {
            System.out.println( "  [GuiFontTest] " + what );
        }
        return cond;
    }
}
