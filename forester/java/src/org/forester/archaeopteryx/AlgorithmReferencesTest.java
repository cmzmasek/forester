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
        // the tanglegram auto-untangle heuristic must be described, with its barycentre + tanglegram citations
        if ( !text.contains( "Auto-untangle" ) || !text.contains( "barycentre" ) || !text.contains( "Sugiyama" )
                || !text.contains( "Tanglegrams for rooted phylogenetic trees" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] tanglegram auto-untangle reference/citation is missing" );
            ok = false;
        }
        // the tanglegram entanglement / congruence score must be described, with its rank-correlation basis
        if ( !text.contains( "entanglement" ) || !text.contains( "discordant" )
                || !text.contains( "A New Measure of Rank Correlation" ) || !text.contains( "Kendall" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] tanglegram entanglement reference/citation is missing" );
            ok = false;
        }
        // the geologic time axis must be described, with its ICS / Cohen et al. chart citation
        if ( !text.contains( "Geologic time axis" ) || !text.contains( "International Chronostratigraphic Chart" )
                || !text.contains( "International Commission on Stratigraphy" )
                || !text.contains( "The ICS International Chronostratigraphic Chart this decade" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] geologic time axis reference/citation is missing" );
            ok = false;
        }
        // the fossil range bars (FAD/LAD) must be described, with the strap (Bell & Lloyd 2015) citation
        if ( !text.contains( "Fossil range bars" ) || !text.contains( "First Appearance Datum" )
                || !text.contains( "Last Appearance Datum" ) || !text.contains( "Bell MA, Lloyd GT (2015)" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] fossil range bars (FAD/LAD) reference/citation is missing" );
            ok = false;
        }
        // the sequence-alignment display must be described, with the Jalview (Waterhouse et al. 2009) colouring citation
        if ( !text.contains( "Sequence-alignment display" ) || !text.contains( "Zappo" )
                || !text.contains( "Jalview Version 2" ) || !text.contains( "Bioinformatics 25(9):1189" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] sequence-alignment display reference/citation is missing" );
            ok = false;
        }
        // the node-age bars/spindles must be described (incl. the honest "schematic, not raw posterior" caveat) + BEAST
        if ( !text.contains( "Node age bars / spindles" ) || !text.contains( "Highest Posterior Density" )
                || !text.contains( "SCHEMATIC" ) || !text.contains( "BEAST 1.10" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] node age bars/spindles reference/citation is missing" );
            ok = false;
        }
        // the broken-long-branch display convention must be described, with its median-multiple threshold
        if ( !text.contains( "Break Long Branches" ) || !text.contains( "axis-break glyph" )
                || !text.contains( "MEDIAN of the tree's strictly-positive branch lengths" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] broken long-branch reference is missing" );
            ok = false;
        }
        // the Auspice/Nextstrain JSON import must be described, with the Nextstrain (Hadfield et al. 2018) citation
        if ( !text.contains( "Auspice / Nextstrain JSON import" ) || !text.contains( "dataset.json" )
                || !text.contains( "Hadfield J" ) || !text.contains( "real-time tracking of pathogen evolution" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] Auspice/Nextstrain import reference/citation is missing" );
            ok = false;
        }
        // the NCBI-taxonomy-as-proxy reconciliation must be described (honest "approximate/classification not a phylogeny"
        // caveat) with the NCBI Taxonomy (Schoch et al. 2020) citation
        if ( !text.contains( "Reconciliation using the NCBI taxonomy" ) || !text.contains( "APPROXIMATE" )
                || !text.contains( "not a phylogeny" ) || !text.contains( "Schoch CL, Ciufo S" )
                || !text.contains( "NCBI Taxonomy: a comprehensive update" ) || !text.contains( "baaa062" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] NCBI-taxonomy reconciliation reference/citation is missing" );
            ok = false;
        }
        // the GTDB taxonomy import must be described, with the GTDB (Parks et al. 2022) + GTDB-Tk (Chaumeil et al. 2020) citations
        if ( !text.contains( "GTDB taxonomy import" ) || !text.contains( "genome-based" )
                || !text.contains( "Parks DH, Chuvochina M" ) || !text.contains( "GTDB-Tk" )
                || !text.contains( "Chaumeil P-A, Mussig AJ" ) ) {
            System.out.println( "  [AlgorithmReferencesTest] GTDB taxonomy import reference/citation is missing" );
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
