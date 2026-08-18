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
import java.awt.GraphicsEnvironment;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JTabbedPane;
import javax.swing.ListCellRenderer;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;

/**
 * Integration test for {@link SettingsDialog}: builds the dialog over a real {@link MainFrame} and checks
 * that it has the expected tabs and that a bound checkbox drives its backing menu item (the {@code doClick}
 * binding that all the apply logic depends on). Also checks {@link SettingsDialog#prettyEnumName} (the
 * Node-shape/Node-fill display fix: "RECTANGLE" -> "Rectangle") purely, plus -- headful -- that those combos
 * actually render through it. The pure part runs everywhere; the GUI part is a no-op on a headless box.
 */
public final class SettingsDialogTest {

    public static boolean test() {
        if ( !prettyEnumNameOk() ) {
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; nothing meaningful to do without a display toolkit
        }
        try {
            final Phylogeny phy = new Phylogeny();
            final PhylogenyNode root = new PhylogenyNode();
            root.addAsChild( new PhylogenyNode() );
            root.addAsChild( new PhylogenyNode() );
            phy.setRoot( root );
            phy.externalNodesHaveChanged();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, conf, "settings test" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                // the menu bar must carry a "Settings" launcher and no longer the old "Options"/"Type" menus
                final JMenuBar bar = ( (JFrame) mf[ 0 ] ).getJMenuBar();
                boolean has_settings = false;
                boolean has_old = false;
                for ( int i = 0; i < bar.getMenuCount(); ++i ) {
                    final JMenu m = bar.getMenu( i );
                    if ( m == null ) {
                        continue;
                    }
                    if ( "Settings".equals( m.getText() ) ) {
                        has_settings = true;
                    }
                    if ( "Options".equals( m.getText() ) || "Type".equals( m.getText() ) ) {
                        has_old = true;
                    }
                }
                if ( !has_settings || has_old ) {
                    ok[ 0 ] = false;
                }
                final SettingsDialog dlg = new SettingsDialog( mf[ 0 ] );
                dlg.pack();
                final List<JTabbedPane> tabs = new ArrayList<>();
                collect( dlg.getContentPane(), JTabbedPane.class, tabs );
                // the former single "Display" tab was split into "Layout" / "Labels & Colors" / "Overlays", and the
                // former "Search" tab (its one clumsy colorize-all-found setting) was removed entirely -> 8 tabs
                if ( tabs.isEmpty() || ( tabs.get( 0 ).getTabCount() != 8 ) ) {
                    ok[ 0 ] = false;
                }
                else {
                    final List<String> titles = new ArrayList<>();
                    for ( int i = 0; i < tabs.get( 0 ).getTabCount(); ++i ) {
                        titles.add( tabs.get( 0 ).getTitleAt( i ) );
                    }
                    // the three split tabs are present, the retired combined "Display" and standalone "Search" tabs
                    // are gone, and the persistent-taxonomy-cache tab (with its on/off checkbox) is still there
                    if ( !titles.contains( "Layout" ) || !titles.contains( "Labels & Colors" )
                            || !titles.contains( "Overlays" ) || titles.contains( "Display" )
                            || titles.contains( "Search" ) || !titles.contains( "Taxonomy Cache" )
                            || ( findCheckBox( dlg.getContentPane(), "Use persistent cache" ) == null ) ) {
                        ok[ 0 ] = false;
                    }
                }
                // the default width is widened so the whole tab-header row fits on ONE row (pack() alone sizes to
                // the narrow tab content and would wrap the header onto two cramped rows)
                if ( dlg.getWidth() < 900 ) {
                    ok[ 0 ] = false;
                }
                // the retired "Behavior" section and its "Data returned on copy:" combo (the old "List Node
                // Data" config) must no longer appear anywhere in the dialog
                if ( ( findLabel( dlg.getContentPane(), "Behavior" ) != null )
                        || ( findLabel( dlg.getContentPane(), "Data returned on copy:" ) != null ) ) {
                    ok[ 0 ] = false;
                }
                // the "Reset to Defaults" button must be present in the footer (its action is covered by
                // ResetToDefaultsTest; the modal confirm can't be clicked through here)
                boolean has_reset = false;
                final List<javax.swing.JButton> buttons = new ArrayList<>();
                collect( dlg.getContentPane(), javax.swing.JButton.class, buttons );
                for ( final javax.swing.JButton b : buttons ) {
                    if ( ( b.getText() != null ) && b.getText().startsWith( "Reset to Defaults" ) ) {
                        has_reset = true;
                    }
                }
                if ( !has_reset ) {
                    ok[ 0 ] = false;
                }
                // a dialog checkbox must drive its backing menu item (the doClick binding)
                final boolean before = mf[ 0 ]._show_scale_cbmi.isSelected();
                final JCheckBox cb = findCheckBox( dlg.getContentPane(), mf[ 0 ]._show_scale_cbmi.getText() );
                if ( cb == null ) {
                    ok[ 0 ] = false;
                }
                else {
                    cb.doClick();
                    if ( mf[ 0 ]._show_scale_cbmi.isSelected() == before ) {
                        ok[ 0 ] = false; // toggling the dialog control did not flip the menu item
                    }
                }
                // the Node-shape / Node-fill combos must render enum constants proper-cased, not ALL-CAPS: find
                // the actual combos in the dialog and confirm their installed renderer turns RECTANGLE/SOLID into
                // "Rectangle"/"Solid" (guards that the pretty renderer is really wired to these two combos)
                final List<JComboBox> combos = new ArrayList<>();
                collect( dlg.getContentPane(), JComboBox.class, combos );
                boolean saw_shape = false;
                boolean saw_fill = false;
                for ( final JComboBox<?> combo : combos ) {
                    if ( ( combo.getItemCount() > 0 ) && ( combo.getItemAt( 0 ) instanceof NodeShape ) ) {
                        saw_shape = true;
                        if ( !"Rectangle".equals( renderedText( combo, NodeShape.RECTANGLE ) ) ) {
                            ok[ 0 ] = false;
                        }
                    }
                    else if ( ( combo.getItemCount() > 0 ) && ( combo.getItemAt( 0 ) instanceof NodeFill ) ) {
                        saw_fill = true;
                        if ( !"Solid".equals( renderedText( combo, NodeFill.SOLID ) ) ) {
                            ok[ 0 ] = false;
                        }
                    }
                }
                if ( !saw_shape || !saw_fill ) {
                    ok[ 0 ] = false; // the combos we mean to check must actually be present
                }
                dlg.dispose();
                ( (JFrame) mf[ 0 ] ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Pure check of {@link SettingsDialog#prettyEnumName}: single words are title-cased and underscores
     *  become spaces. Runs even when headless -- the display-string logic needs no toolkit. */
    private static boolean prettyEnumNameOk() {
        boolean ok = true;
        ok &= expectPretty( "Rectangle", SettingsDialog.prettyEnumName( NodeShape.RECTANGLE ) );
        ok &= expectPretty( "Circle", SettingsDialog.prettyEnumName( NodeShape.CIRCLE ) );
        ok &= expectPretty( "Default", SettingsDialog.prettyEnumName( NodeShape.DEFAULT ) );
        ok &= expectPretty( "None", SettingsDialog.prettyEnumName( NodeFill.NONE ) );
        ok &= expectPretty( "Solid", SettingsDialog.prettyEnumName( NodeFill.SOLID ) );
        // multi-word (underscore) constant: each word title-cased, underscore -> space
        ok &= expectPretty( "Lower Left",
                            SettingsDialog.prettyEnumName( Options.OVERVIEW_PLACEMENT_TYPE.LOWER_LEFT ) );
        return ok;
    }

    private static boolean expectPretty( final String expected, final String actual ) {
        if ( !expected.equals( actual ) ) {
            System.out.println( "  [SettingsDialogTest] prettyEnumName expected \"" + expected + "\", got \""
                    + actual + "\"" );
            return false;
        }
        return true;
    }

    /** The text the combo's installed renderer produces for {@code value} (how the item actually appears). */
    @SuppressWarnings( { "unchecked", "rawtypes" } )
    private static String renderedText( final JComboBox combo, final Object value ) {
        final ListCellRenderer r = combo.getRenderer();
        final Component c = r.getListCellRendererComponent( new JList(), value, -1, false, false );
        return ( c instanceof JLabel ) ? ( (JLabel) c ).getText() : null;
    }

    @SuppressWarnings( "unchecked" )
    private static <T> void collect( final Container c, final Class<T> type, final List<T> out ) {
        for ( final Component comp : c.getComponents() ) {
            if ( type.isInstance( comp ) ) {
                out.add( (T) comp );
            }
            if ( comp instanceof Container ) {
                collect( (Container) comp, type, out );
            }
        }
    }

    private static JCheckBox findCheckBox( final Container c, final String text ) {
        final List<JCheckBox> all = new ArrayList<>();
        collect( c, JCheckBox.class, all );
        for ( final JCheckBox cb : all ) {
            if ( text.equals( cb.getText() ) ) {
                return cb;
            }
        }
        return null;
    }

    private static JLabel findLabel( final Container c, final String text ) {
        final List<JLabel> all = new ArrayList<>();
        collect( c, JLabel.class, all );
        for ( final JLabel l : all ) {
            if ( text.equals( l.getText() ) ) {
                return l;
            }
        }
        return null;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private SettingsDialogTest() {
    }
}
