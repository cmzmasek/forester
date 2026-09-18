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
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.tools.NodeDataImporter;
import org.forester.phylogeny.Phylogeny;

/**
 * The wiring of <b>View &rarr; Order Matrix Columns</b> (the orders themselves are {@link MatrixColumnOrderTest}'s):
 * the submenu and its five radio items; the Clustergram preset ordering by the tab's mode, Clustered by default; each
 * radio item re-ordering the CURRENT tab's matrix and Manual never re-sorting; the mode being PER TAB with the radios
 * following a tab switch; a restored figure setting its tab to Manual; and Reset to Defaults putting every tab back
 * to Clustered. Headful; a green no-op when headless.
 * <p>
 * The fixture is the pangenome demo (40 genes), and the test first asserts that its four data-driven orders are
 * pairwise DIFFERENT -- otherwise a radio that did nothing could pass for one that worked.
 */
public final class MatrixColumnOrderMenuTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MatrixColumnOrderMenu: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return menuOk();
    }

    private static boolean fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [MatrixColumnOrderMenuTest] " + msg );
        ok[ 0 ] = false;
        return false;
    }

    /** The demo tree with its table joined on, as the import would leave it. */
    private static Phylogeny demo() throws Exception {
        final File dir = new File( System.getProperty( "user.dir" ), "forester/demo" );
        final Phylogeny phy = FigureRenderer.readTrees( new File( dir, "pangenome-presence-absence.xml" ) )[ 0 ];
        final NodeDataImporter.Table table = NodeDataImporter
                .parseTable( Files.readString( new File( dir, "pangenome-presence-absence.tsv" ).toPath() ) );
        NodeDataImporter.apply( phy, table, table.defaultKeyColumn(), NodeDataImporter.MatchBy.TIP_NAME );
        return phy;
    }

    private static boolean menuOk() {
        final boolean[] ok = { true };
        try {
            final Phylogeny first = demo();
            final Phylogeny second = demo();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { first }, new Configuration(), "order" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    ( ( JFrame ) frame ).setSize( 1000, 700 );
                    checkMenu( ok, frame );
                    final TreePanel tp0 = frame.getMainPanel().getCurrentTreePanel();
                    final Phylogeny phy = tp0.getPhylogeny();
                    final List<String> table = new ArrayList<String>(
                            MatrixColumnOrder.matrixRefs( MainFrame.clustergramColumnSpecs( phy ) ) );
                    final List<String> alpha = MatrixColumnOrder.order( table, MatrixColumnOrder.Mode.ALPHABETICAL, phy );
                    final List<String> freq = MatrixColumnOrder.order( table, MatrixColumnOrder.Mode.FREQUENCY, phy );
                    final List<String> clus = MatrixColumnOrder.order( table, MatrixColumnOrder.Mode.CLUSTERED, phy );
                    // GUARD: four orders that differ pairwise, or a radio that did nothing could pass
                    final List<List<String>> all = java.util.Arrays.asList( table, alpha, freq, clus );
                    for( int i = 0; i < all.size(); ++i ) {
                        for( int j = i + 1; j < all.size(); ++j ) {
                            if ( all.get( i ).equals( all.get( j ) ) ) {
                                fail( ok, "the fixture cannot tell two modes apart (orders " + i + " and " + j + ")" );
                            }
                        }
                    }
                    if ( table.size() != 40 ) {
                        fail( ok, "the demo should give 40 matrix columns, got " + table.size() );
                    }
                    // the preset orders by the tab's mode: Clustered by default
                    frame.applyClustergramPreset();
                    expect( ok, frame, tp0, clus, MatrixColumnOrder.Mode.CLUSTERED, "View > Clustergram (default)" );
                    // each radio re-orders the current tab
                    click( frame, MatrixColumnOrder.Mode.TABLE );
                    expect( ok, frame, tp0, table, MatrixColumnOrder.Mode.TABLE, "Same as Table" );
                    click( frame, MatrixColumnOrder.Mode.ALPHABETICAL );
                    expect( ok, frame, tp0, alpha, MatrixColumnOrder.Mode.ALPHABETICAL, "Alphabetical" );
                    click( frame, MatrixColumnOrder.Mode.FREQUENCY );
                    expect( ok, frame, tp0, freq, MatrixColumnOrder.Mode.FREQUENCY, "Frequency" );
                    // Manual keeps what is on screen -- here the frequency order -- and re-sorts nothing
                    click( frame, MatrixColumnOrder.Mode.MANUAL );
                    expect( ok, frame, tp0, freq, MatrixColumnOrder.Mode.MANUAL, "Manual (keeps the current order)" );
                    // PER TAB: a new tab starts Clustered, and the radios follow a switch both ways
                    AptxUtil.addPhylogeniesToTabs( new Phylogeny[] { second }, "", "", frame.getConfiguration(),
                                                   frame.getMainPanel() );
                    final TreePanel tp1 = frame.getMainPanel().getCurrentTreePanel();
                    if ( ( tp1 == tp0 ) || ( tp1.getMatrixColumnOrder() != MatrixColumnOrder.DEFAULT ) ) {
                        fail( ok, "a new tab must start in the default mode, got " + tp1.getMatrixColumnOrder() );
                    }
                    radioShows( ok, frame, MatrixColumnOrder.Mode.CLUSTERED, "on the new tab" );
                    frame.getMainPanel().getTabbedPane().setSelectedIndex( 0 );
                    radioShows( ok, frame, MatrixColumnOrder.Mode.MANUAL, "after switching back to the first tab" );
                    if ( tp1.getMatrixColumnOrder() != MatrixColumnOrder.Mode.CLUSTERED ) {
                        fail( ok, "the first tab's Manual must not leak into the second" );
                    }
                    frame.getMainPanel().getTabbedPane().setSelectedIndex( 1 );
                    radioShows( ok, frame, MatrixColumnOrder.Mode.CLUSTERED, "after switching to the second tab" );
                    // a restored figure restores an EXPLICIT order: its tab becomes Manual and nothing is re-sorted
                    click( frame, MatrixColumnOrder.Mode.TABLE ); // tab 1 has no columns yet -- only its mode changes
                    FigureSpec.capture( tp0 ).applyTo( tp1 );
                    if ( tp1.getMatrixColumnOrder() != MatrixColumnOrder.Mode.MANUAL ) {
                        fail( ok, "restoring a figure with a matrix must set its tab to Manual, got "
                                + tp1.getMatrixColumnOrder() );
                    }
                    if ( !MatrixColumnOrder.matrixRefs( tp1.getAnnotationColumnSpecs() ).equals( freq ) ) {
                        fail( ok, "a restored figure must keep the saved column order exactly" );
                    }
                    // Reset to Defaults: every tab back to Clustered, and the radios with it
                    frame.resetToDefaults();
                    for( final TreePanel tp : frame.getMainPanel().getTreePanels() ) {
                        if ( tp.getMatrixColumnOrder() != MatrixColumnOrder.DEFAULT ) {
                            fail( ok, "Reset to Defaults must put every tab back to Clustered, got "
                                    + tp.getMatrixColumnOrder() );
                        }
                    }
                    radioShows( ok, frame, MatrixColumnOrder.Mode.CLUSTERED, "after Reset to Defaults" );
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                    t.printStackTrace();
                }
                finally {
                    ( ( JFrame ) frame ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    /** The submenu exists in View, is named, uses the menu font, and holds one radio per mode, default selected. */
    private static void checkMenu( final boolean[] ok, final MainFrame frame ) {
        if ( frame._matrix_order_items.size() != MatrixColumnOrder.Mode.values().length ) {
            fail( ok, "one radio item per mode expected, got " + frame._matrix_order_items.size() );
            return;
        }
        JMenu submenu = null;
        for( final MatrixColumnOrder.Mode m : MatrixColumnOrder.Mode.values() ) {
            final JRadioButtonMenuItem item = frame._matrix_order_items.get( m );
            if ( !m.label().equals( item.getText() ) ) {
                fail( ok, "radio for " + m + " must read \"" + m.label() + "\", reads \"" + item.getText() + "\"" );
            }
            if ( !MainFrame.menu_font.equals( item.getFont() ) ) {
                fail( ok, "radio \"" + item.getText() + "\" must use the menu font" );
            }
            if ( item.isSelected() != ( m == MatrixColumnOrder.DEFAULT ) ) {
                fail( ok, "only the default mode's radio may be selected at start; " + m + " is " + item.isSelected() );
            }
            final Component parent = item.getParent();
            if ( parent instanceof JPopupMenu ) {
                submenu = ( JMenu ) ( ( JPopupMenu ) parent ).getInvoker();
            }
        }
        if ( ( submenu == null ) || !"Order Matrix Columns".equals( submenu.getText() ) ) {
            fail( ok, "the radios must sit in an \"Order Matrix Columns\" submenu" );
            return;
        }
        if ( !MainFrame.menu_font.equals( submenu.getFont() ) ) {
            fail( ok, "the submenu must use the menu font (createMenu only sets it in custom-colors mode)" );
        }
        boolean in_view = false;
        for( final Component c : frame._view_jmenu.getMenuComponents() ) {
            in_view |= ( c == submenu );
        }
        if ( !in_view ) {
            fail( ok, "the submenu must be in the View menu" );
        }
    }

    private static void click( final MainFrame frame, final MatrixColumnOrder.Mode mode ) {
        frame._matrix_order_items.get( mode ).doClick(); // a real click, through actionPerformed's dispatch
    }

    private static void expect( final boolean[] ok, final MainFrame frame, final TreePanel tp, final List<String> order,
                                final MatrixColumnOrder.Mode mode, final String what ) {
        final List<String> got = MatrixColumnOrder.matrixRefs( tp.getAnnotationColumnSpecs() );
        if ( !got.equals( order ) ) {
            fail( ok, what + ": the matrix is not in that order" );
        }
        if ( tp.getMatrixColumnOrder() != mode ) {
            fail( ok, what + ": the tab's mode should be " + mode + ", is " + tp.getMatrixColumnOrder() );
        }
        radioShows( ok, frame, mode, what );
    }

    private static void radioShows( final boolean[] ok, final MainFrame frame, final MatrixColumnOrder.Mode mode,
                                    final String when ) {
        if ( !frame._matrix_order_items.get( mode ).isSelected() ) {
            fail( ok, when + ": the " + mode + " radio should be selected" );
        }
    }
}
