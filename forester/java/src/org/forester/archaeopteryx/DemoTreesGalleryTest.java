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
import java.io.InputStream;

import java.awt.Component;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;

/**
 * Guards the File -> Demo Trees gallery ({@link DemoTrees}): (1) headless -- every bundled demo resource loads + parses
 * from the CLASSPATH (i.e. it ships in the jar, not just the demo folder), and the catalog is the expected size; (2)
 * headful -- every catalog entry OPENS without throwing and actually turns on the overlay it advertises (color-by set,
 * annotation columns present, ancestral-pie trait set, the time axis applied), and the tanglegram window is created.
 */
public final class DemoTreesGalleryTest {

    // the resources the catalog opens (a fossil-only ammonite tree, a tanglegram pair + its association file, etc.)
    private static final String[] TREE_RESOURCES = { "color-by-property.xml", "annotation-columns.xml",
            "symbol-columns.xml", "stacked-bar-columns.xml", "domain-architectures.xml", "alignment.xml",
            "bat-phylogeny.xml",
            "ancestral-pie-charts.xml",
            "node-hpd-bars.xml",
            "long-branch-break.xml", "sars-cov-2-time-tree.xml", "nextstrain-ncov.json", "filoviridae-tree.xml",
            "dinosaur-time-tree.xml", "lagomorph-time-tree.xml", "ammonite-time-tree.xml", "tree-of-life-deep-time.xml",
            "tanglegram-host-tree.xml", "tanglegram-parasite-tree.xml", "gtdb-genomes.xml" };

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DemoTreesGallery: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( !resourcesLoadOk() ) {
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // the open-and-configure checks need a realized GUI
        }
        return galleryOpensOk();
    }

    /** Headless: every bundled demo resource loads + parses from the classpath, and the catalog is the expected size. */
    private static boolean resourcesLoadOk() {
        try {
            if ( DemoTrees.catalog().size() != 20 ) {
                return fail( "the demo catalog should have 20 curated entries, has " + DemoTrees.catalog().size() );
            }
            // if the demo resources were not staged onto the classpath (a raw IDE compile that skipped the Ant
            // copy_resources step), skip the load checks rather than fail the whole suite -- the authoritative
            // 'ant compile'/'ant all' path always stages them, so this still guards the jar bundling there
            try ( final InputStream probe = DemoTrees.class.getResourceAsStream( "/resources/demo/color-by-property.xml" ) ) {
                if ( probe == null ) {
                    System.out.println( "  [DemoTreesGalleryTest] demo resources not on the classpath (run 'ant "
                            + "copy_resources'); skipping the bundle checks" );
                    return true;
                }
            }
            for ( final String r : TREE_RESOURCES ) {
                final Phylogeny phy = DemoTrees.loadTree( r );
                if ( ( phy == null ) || phy.isEmpty() || ( phy.getNumberOfExternalNodes() < 2 ) ) {
                    return fail( "bundled demo resource did not load a non-trivial tree: " + r );
                }
            }
            final String assoc = DemoTrees.loadText( "tanglegram-association.tsv" );
            if ( ForesterUtilShim.isBlank( assoc ) ) {
                return fail( "the bundled tanglegram association TSV is empty" );
            }
            final String gtdb = DemoTrees.loadText( "gtdb-classifications.tsv" );
            if ( ForesterUtilShim.isBlank( gtdb ) || !gtdb.contains( "d__Bacteria" ) ) {
                return fail( "the bundled GTDB classifications TSV is empty or missing a d__ classification" );
            }
            return true;
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    /** Headful: each catalog entry opens without throwing AND turns on the overlay it advertises. */
    private static boolean galleryOpensOk() {
        final boolean[] ok = { true };
        try {
            final Phylogeny seed = DemoTrees.loadTree( "color-by-property.xml" );
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { seed }, conf, "demo" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrameApplication frame = (MainFrameApplication) mf[ 0 ];
                try {
                    ( (JFrame) frame ).setSize( 1000, 640 ); // give showWhole a viewport to fit into
                    assertSubmenu( ok, frame ); // the File -> Demo Trees submenu is built with one item per catalog entry
                    for ( final DemoTrees.Demo demo : DemoTrees.catalog() ) {
                        try {
                            demo.open( frame );
                        }
                        catch ( final Throwable t ) {
                            fail( ok, "demo \"" + demo.label() + "\" threw on open: " + t );
                            continue;
                        }
                        assertConfigured( ok, frame, demo.label() );
                    }
                    // the tanglegram opens its own window (not a tab) -- confirm a TanglegramFrame was created
                    boolean tanglegram_shown = false;
                    for ( final Window w : Window.getWindows() ) {
                        if ( w instanceof TanglegramFrame ) {
                            tanglegram_shown = true;
                        }
                    }
                    if ( !tanglegram_shown ) {
                        fail( ok, "the Tanglegram demo must open a TanglegramFrame window" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    for ( final Window w : Window.getWindows() ) {
                        w.dispose(); // dispose the main frame AND the tanglegram window
                    }
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** The File -> Demo Trees submenu must be built with exactly one item per catalog entry, in order, labelled to
     *  match (covers MainFrameApplication.buildDemoTreesSubmenu). */
    private static void assertSubmenu( final boolean[] ok, final MainFrameApplication frame ) {
        JMenu demo_menu = null;
        if ( frame._file_jmenu != null ) {
            for ( final Component c : frame._file_jmenu.getMenuComponents() ) {
                if ( ( c instanceof JMenu ) && "Demo Trees".equals( ( (JMenu) c ).getText() ) ) {
                    demo_menu = (JMenu) c;
                    break;
                }
            }
        }
        if ( demo_menu == null ) {
            fail( ok, "the File menu must contain a \"Demo Trees\" submenu" );
            return;
        }
        final java.util.List<DemoTrees.Demo> catalog = DemoTrees.catalog();
        if ( demo_menu.getItemCount() != catalog.size() ) {
            fail( ok, "the Demo Trees submenu must have " + catalog.size() + " items, has " + demo_menu.getItemCount() );
            return;
        }
        for ( int i = 0; i < catalog.size(); ++i ) {
            final JMenuItem item = demo_menu.getItem( i );
            if ( ( item == null ) || !catalog.get( i ).label().equals( item.getText() ) ) {
                fail( ok, "the Demo Trees submenu item " + i + " must be labelled \"" + catalog.get( i ).label()
                        + "\", was \"" + ( ( item == null ) ? "null" : item.getText() ) + "\"" );
            }
        }
    }

    /** After a demo opens, assert the overlay it advertises is actually on (for the tab-based, configurable demos). */
    private static void assertConfigured( final boolean[] ok, final MainFrameApplication frame, final String label ) {
        final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
        if ( tp == null ) {
            return; // the tanglegram demo doesn't add a tab; nothing to check here
        }
        if ( label.startsWith( "Color Tips" ) ) {
            if ( ( tp.getPropertyColorScheme() == null )
                    || !"data:host".equals( tp.getPropertyColorScheme().getRef() ) ) {
                fail( ok, "the color-by demo must colour by data:host" );
            }
        }
        else if ( label.startsWith( "Filoviridae" ) ) {
            if ( ( tp.getPropertyColorScheme() == null )
                    || !"repseq:species".equals( tp.getPropertyColorScheme().getRef() ) ) {
                fail( ok, "the Filoviridae demo must colour by repseq:species" );
            }
        }
        else if ( label.startsWith( "Annotation Columns" ) ) {
            if ( !tp.hasAnnotationColumns() ) {
                fail( ok, "the annotation-columns demo must have annotation columns" );
            }
        }
        else if ( label.startsWith( "Alignment" ) ) {
            if ( !tp.getOptions().isShowMsa() || !AptxUtil.hasAlignedSequences( tp.getPhylogeny() ) ) {
                fail( ok, "the alignment demo must auto-enable the alignment display over aligned tip sequences" );
            }
        }
        else if ( label.startsWith( "Symbol Columns" ) ) {
            if ( !tp.hasAnnotationColumns() ) {
                fail( ok, "the symbol-columns demo must have annotation columns" );
            }
            else {
                boolean has_symbol = false;
                for ( final AnnotationColumns.ColumnSpec s : tp.getAnnotationColumnSpecs() ) {
                    if ( s._type == AnnotationColumns.Type.SYMBOL ) {
                        has_symbol = true;
                    }
                }
                if ( !has_symbol ) {
                    fail( ok, "the symbol-columns demo must configure a SYMBOL column" );
                }
            }
        }
        else if ( label.startsWith( "Stacked Bar Columns" ) ) {
            if ( !tp.hasAnnotationColumns() ) {
                fail( ok, "the stacked-bar demo must have annotation columns" );
            }
            else {
                boolean has_stacked = false;
                for ( final AnnotationColumns.ColumnSpec s : tp.getAnnotationColumnSpecs() ) {
                    if ( s._type == AnnotationColumns.Type.STACKED_BAR ) {
                        has_stacked = true;
                    }
                }
                if ( !has_stacked ) {
                    fail( ok, "the stacked-bar demo must configure a STACKED_BAR column" );
                }
            }
        }
        else if ( label.startsWith( "Pie Columns" ) ) {
            if ( !tp.hasAnnotationColumns() ) {
                fail( ok, "the pie demo must have annotation columns" );
            }
            else {
                boolean has_pie = false;
                for ( final AnnotationColumns.ColumnSpec s : tp.getAnnotationColumnSpecs() ) {
                    if ( s._type == AnnotationColumns.Type.PIE ) {
                        has_pie = true;
                    }
                }
                if ( !has_pie ) {
                    fail( ok, "the pie demo must configure a PIE column" );
                }
            }
        }
        else if ( label.startsWith( "Ancestral State Pies" ) ) {
            if ( tp.getAncestralPieTrait() == null ) {
                fail( ok, "the ancestral-pie demo must select a discrete trait" );
            }
        }
        else if ( label.startsWith( "Break Long Branches" ) ) {
            if ( !frame.getOptions().isBreakLongBranches() ) {
                fail( ok, "the break-long-branches demo must turn on the Break Long Branches option" );
            }
        }
        else if ( label.startsWith( "Node Age Spindles" ) ) {
            if ( !frame.getOptions().isShowHpdBars()
                    || ( frame.getOptions().getNodeAgeShape() != Options.NODE_AGE_SHAPE.SPINDLE ) ) {
                fail( ok, "the node-age spindles demo must turn on the node-age bars in SPINDLE shape (bars="
                        + frame.getOptions().isShowHpdBars() + " shape=" + frame.getOptions().getNodeAgeShape() + ")" );
            }
        }
        else if ( label.startsWith( "SARS-CoV-2" ) ) {
            if ( tp.effectiveTimeAxisType() != Options.TIME_AXIS_TYPE.CALENDAR ) {
                fail( ok, "the SARS-CoV-2 demo must show the CALENDAR axis, got " + tp.effectiveTimeAxisType() );
            }
        }
        else if ( label.startsWith( "Phylodynamics" ) ) {
            if ( tp.effectiveTimeAxisType() != Options.TIME_AXIS_TYPE.CALENDAR ) {
                fail( ok, "the Nextstrain JSON demo must show the CALENDAR axis, got " + tp.effectiveTimeAxisType() );
            }
            if ( !"region".equals( tp.getAncestralPieTrait() ) ) {
                fail( ok, "the Nextstrain JSON demo must show region ancestral-state pies, got "
                        + tp.getAncestralPieTrait() );
            }
        }
        else if ( label.startsWith( "Dinosaur" ) || label.startsWith( "Ammonite" ) || label.startsWith( "Tree of Life" )
                || label.startsWith( "Lagomorph" ) ) {
            if ( tp.effectiveTimeAxisType() != Options.TIME_AXIS_TYPE.GEOLOGIC ) {
                fail( ok, "the " + label + " demo must show the GEOLOGIC axis, got " + tp.effectiveTimeAxisType() );
            }
        }
        else if ( label.startsWith( "Bat Phylogeny" ) ) {
            if ( !tp.hasRankLegend() ) {
                fail( ok, "the bat-phylogeny demo must colorize by taxonomic rank (family) from the in-tree clade "
                        + "annotations" );
            }
        }
        else if ( label.startsWith( "GTDB Taxonomy" ) ) {
            if ( ( tp.getPropertyColorScheme() == null )
                    || !"gtdb:phylum".equals( tp.getPropertyColorScheme().getRef() ) ) {
                fail( ok, "the GTDB demo must import the classifications and colour by gtdb:phylum, got "
                        + ( tp.getPropertyColorScheme() == null ? "no color scheme"
                                                                : tp.getPropertyColorScheme().getRef() ) );
            }
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [DemoTreesGalleryTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [DemoTreesGalleryTest] " + msg );
        ok[ 0 ] = false;
    }

    /** Tiny local blank check (avoids taking a dependency on a specific ForesterUtil signature). */
    private static final class ForesterUtilShim {
        static boolean isBlank( final String s ) {
            return ( s == null ) || ( s.trim().length() == 0 );
        }
    }

    private DemoTreesGalleryTest() {
    }
}
