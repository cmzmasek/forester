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
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterConstants;

/**
 * Tests the Help &rarr; References content ({@link AlgorithmReferences}) and its wiring: the reference table is
 * complete (every algorithm has a non-blank label and citation), the rendered text and the dialog view carry
 * every algorithm, and the Help menu actually exposes a "References" item. The content checks are headless-safe;
 * the menu-presence check needs a display (FlatLaf via {@code createInstance}) and is skipped when headless.
 */
public final class AlgorithmReferencesTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "AlgorithmReferences: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        boolean ok = checkContent();
        if ( !GraphicsEnvironment.isHeadless() ) {
            ok = checkMenu() && ok;
        }
        return ok;
    }

    private static boolean checkContent() {
        boolean ok = true;
        final java.util.List<AlgorithmReferences.Reference> all = AlgorithmReferences.all();
        if ( all.size() < 5 ) {
            System.out.println( "  [AlgorithmReferencesTest] reference list is suspiciously short: " + all.size() );
            ok = false;
        }
        // every algorithm must have a non-blank label AND a non-blank citation (no silently unlabeled/empty rows)
        for ( final AlgorithmReferences.Reference r : all ) {
            if ( ( r.algorithm() == null ) || r.algorithm().isBlank() || ( r.citation() == null )
                    || r.citation().isBlank() ) {
                System.out.println( "  [AlgorithmReferencesTest] a reference has a blank algorithm or citation: " + r );
                ok = false;
            }
        }
        final String text = AlgorithmReferences.asText();
        // the rendered text must mention every algorithm label
        for ( final AlgorithmReferences.Reference r : all ) {
            if ( !text.contains( r.algorithm() ) ) {
                System.out.println( "  [AlgorithmReferencesTest] rendered text is missing an algorithm: " + r.algorithm() );
                ok = false;
            }
        }
        // spot-check known citations
        if ( !text.contains( "minimal ancestor deviation" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] MAD-rooting citation is missing" );
            ok = false;
        }
        boolean phyloxml_present = false;
        for ( final AlgorithmReferences.Reference r : all ) {
            if ( ForesterConstants.PHYLO_XML_REFERENCE.equals( r.citation() ) ) {
                phyloxml_present = true;
            }
        }
        if ( !phyloxml_present ) {
            System.out.println( "  [AlgorithmReferencesTest] phyloXML reference (ForesterConstants) not reused verbatim" );
            ok = false;
        }
        // the dialog view (headless-safe to build) must show exactly the reference text
        final JScrollPane view = MainFrame.buildReferencesView();
        final JTextArea ta = (JTextArea) view.getViewport().getView();
        if ( !text.equals( ta.getText() ) ) {
            System.out.println( "  [AlgorithmReferencesTest] references view text does not match AlgorithmReferences.asText()" );
            ok = false;
        }
        if ( ta.isEditable() ) {
            System.out.println( "  [AlgorithmReferencesTest] references view should be read-only (but still selectable)" );
            ok = false;
        }
        return ok;
    }

    private static boolean checkMenu() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] {}, conf, "refs" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                if ( helpItem( mf[ 0 ].getJMenuBar(), "References" ) == null ) {
                    System.out.println( "  [AlgorithmReferencesTest] Help menu has no \"References\" item" );
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

    /** Finds a menu item with the given label inside the "Help" menu of the menu bar. */
    private static JMenuItem helpItem( final JMenuBar bar, final String label ) {
        if ( bar == null ) {
            return null;
        }
        for ( int i = 0; i < bar.getMenuCount(); i++ ) {
            final JMenu menu = bar.getMenu( i );
            if ( ( menu != null ) && "Help".equals( menu.getText() ) ) {
                for ( int j = 0; j < menu.getItemCount(); j++ ) {
                    final JMenuItem item = menu.getItem( j );
                    if ( ( item != null ) && label.equals( item.getText() ) ) {
                        return item;
                    }
                }
            }
        }
        return null;
    }

    private AlgorithmReferencesTest() {
    }
}
