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
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Integration test asserting that every top-level menu carries a mouse-over tooltip (as the
 * "Settings" menu does). Checks for a non-empty {@code getToolTipText()} on each menu rather than
 * exact wording, so it survives copy edits. Guarded to a no-op on a headless box. Needs FlatLaf on
 * the classpath (via {@code MainFrameApplication.createInstance}), so it is run standalone, not as
 * part of the headless suite.
 */
public final class MenuTooltipsTest {

    // "Font Size" was retired as a top-level menu (it is now a slider in the control panel).
    private static final String[] EXPECTED_MENUS = { "File", "Analysis", "Tools", "View", "Settings", "Help" };

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MenuTooltips: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { twoLeafTree() }, conf,
                                                                         "menus" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final JMenuBar bar = mf[ 0 ].getJMenuBar();
                for( final String title : EXPECTED_MENUS ) {
                    final JMenu menu = menuByTitle( bar, title );
                    if ( menu == null ) {
                        System.out.println( "  [MenuTooltipsTest] menu not found: " + title );
                        ok[ 0 ] = false;
                    }
                    else {
                        final String tip = menu.getToolTipText();
                        if ( ( tip == null ) || tip.trim().isEmpty() ) {
                            System.out.println( "  [MenuTooltipsTest] missing tooltip on menu: " + title );
                            ok[ 0 ] = false;
                        }
                    }
                }
                // the two data-export items must advertise how to scope the export (discoverability hint)
                final JMenu file = menuByTitle( bar, "File" );
                for( final String label : new String[] { "Export Sequences (FASTA)...", "Export Node Data (TSV)..." } ) {
                    final JMenuItem item = itemByText( file, label );
                    if ( item == null ) {
                        System.out.println( "  [MenuTooltipsTest] export item not found: " + label );
                        ok[ 0 ] = false;
                    }
                    else {
                        final String tip = item.getToolTipText();
                        if ( ( tip == null ) || !tip.contains( "restrict the export" ) ) {
                            System.out.println( "  [MenuTooltipsTest] export item missing scope hint: " + label );
                            ok[ 0 ] = false;
                        }
                    }
                }
                // the import item must advertise that custom columns become color-able properties
                final JMenuItem import_item = itemByText( file, "Import Node Data (TSV)..." );
                if ( import_item == null ) {
                    System.out.println( "  [MenuTooltipsTest] import item not found" );
                    ok[ 0 ] = false;
                }
                else {
                    final String tip = import_item.getToolTipText();
                    if ( ( tip == null ) || !tip.contains( "color by" ) ) {
                        System.out.println( "  [MenuTooltipsTest] import item missing color-by hint" );
                        ok[ 0 ] = false;
                    }
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

    private static JMenu menuByTitle( final JMenuBar bar, final String title ) {
        for( int i = 0; i < bar.getMenuCount(); ++i ) {
            final JMenu m = bar.getMenu( i );
            if ( ( m != null ) && title.equals( m.getText() ) ) {
                return m;
            }
        }
        return null;
    }

    private static JMenuItem itemByText( final JMenu menu, final String text ) {
        if ( menu == null ) {
            return null;
        }
        for( int i = 0; i < menu.getItemCount(); ++i ) {
            final JMenuItem it = menu.getItem( i );
            if ( ( it != null ) && text.equals( it.getText() ) ) {
                return it;
            }
        }
        return null;
    }

    private static Phylogeny twoLeafTree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( new PhylogenyNode() );
        root.addAsChild( new PhylogenyNode() );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private MenuTooltipsTest() {
    }
}
