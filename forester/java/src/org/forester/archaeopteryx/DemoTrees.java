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

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The in-app "demo gallery": a curated, ordered set of the {@code forester/demo} trees, bundled into the jar as
 * resources and opened PRE-CONFIGURED so each one immediately shows off a distinct capability -- like Cytoscape's
 * pre-loaded demo networks, but hardcoded into the (tiny) jar rather than unpacked onto the user's disk. Wired into
 * <b>File &rarr; Demo Trees</b>. The catalog is deliberately a SUBSET of the published demo folder (the showcases, not
 * the regression fixtures); add entries here as new headline features land.
 */
final class DemoTrees {

    /** Where the demo resources live inside the jar (see build.xml copy_resources). */
    private static final String RESOURCE_DIR = "/resources/demo/";

    /** The action that opens (and pre-configures) one demo; may fail on a bad resource -> surfaced as a dialog. */
    @FunctionalInterface
    interface DemoAction {
        void run( MainFrameApplication mf ) throws Exception;
    }

    /** One demo-gallery entry: a menu label, a tooltip, and the open-and-configure action. */
    static final class Demo {

        private final String     _label;
        private final String     _tooltip;
        private final DemoAction  _action;

        Demo( final String label, final String tooltip, final DemoAction action ) {
            _label = label;
            _tooltip = tooltip;
            _action = action;
        }

        String label() {
            return _label;
        }

        String tooltip() {
            return _tooltip;
        }

        void open( final MainFrameApplication mf ) throws Exception {
            _action.run( mf );
        }
    }

    /**
     * The curated, ordered catalog -- approachable first, flagship last. Each entry opens its bundled tree and turns on
     * the overlay it showcases (so the demo is impressive on open, not a plain tree the user must then configure).
     */
    static List<Demo> catalog() {
        final List<Demo> demos = new ArrayList<>();
        demos.add( new Demo( "Color Tips by Metadata",
                             "Colour a tree's tips by a categorical property (host species) -- your data, on the tree.",
                             mf -> {
                                 openTree( mf, "color-by-property.xml" );
                                 colorBy( mf, "data:host" );
                             } ) );
        demos.add( new Demo( "Annotation Columns",
                             "Tip-aligned data columns beside the tree: categorical colour strips + a numeric heat-map.",
                             mf -> {
                                 final TreePanel tp = openTree( mf, "annotation-columns.xml" );
                                 tp.setAnnotationColumns( Arrays.asList(
                                         new AnnotationColumns.ColumnSpec( "data:host", AnnotationColumns.Type.COLOR_STRIP ),
                                         new AnnotationColumns.ColumnSpec( "data:segment", AnnotationColumns.Type.COLOR_STRIP ),
                                         new AnnotationColumns.ColumnSpec( "data:viral_load", AnnotationColumns.Type.HEATMAP ) ) );
                                 refit( mf );
                             } ) );
        demos.add( new Demo( "Protein Domain Architectures",
                             "Multi-domain protein sequences drawn to scale at each tip (auto-enabled on load).",
                             mf -> openTree( mf, "domain-architectures.xml" ) ) ); // domains auto-enable + fit on load
        demos.add( new Demo( "Ancestral State Pies (Phylogeography)",
                             "A discrete geographic trait as posterior-probability pie charts at the internal nodes.",
                             mf -> {
                                 final TreePanel tp = openTree( mf, "ancestral-pie-charts.xml" );
                                 firstAncestralPie( mf, tp );
                             } ) );
        demos.add( new Demo( "SARS-CoV-2 Time Tree (Calendar Axis)",
                             "A tip-dated viral tree read against a labelled calendar-year axis -- molecular epidemiology.",
                             mf -> {
                                 final TreePanel tp = openTree( mf, "sars-cov-2-time-tree.xml" );
                                 timeAxis( mf, tp, Options.TIME_AXIS_TYPE.CALENDAR );
                             } ) );
        demos.add( new Demo( "Dinosaur Time Tree (Geologic Axis)",
                             "A dated archosaur tree -- with Archaeopteryx! -- against the ICS geologic time scale.",
                             mf -> {
                                 final TreePanel tp = openTree( mf, "dinosaur-time-tree.xml" );
                                 timeAxis( mf, tp, Options.TIME_AXIS_TYPE.GEOLOGIC );
                             } ) );
        demos.add( new Demo( "Ammonite Time Tree (Fossil Ranges)",
                             "An all-extinct fossil clade with FAD/LAD stratigraphic-range bars against the geologic axis.",
                             mf -> {
                                 final TreePanel tp = openTree( mf, "ammonite-time-tree.xml" );
                                 timeAxis( mf, tp, Options.TIME_AXIS_TYPE.GEOLOGIC ); // fossil range bars auto-enable
                             } ) );
        demos.add( new Demo( "Tree of Life (Deep Time)",
                             "A time-calibrated tree of life back to LUCA (~3.8 Ga), with Eon/Era geologic bands.",
                             mf -> {
                                 final TreePanel tp = openTree( mf, "tree-of-life-deep-time.xml" );
                                 timeAxis( mf, tp, Options.TIME_AXIS_TYPE.GEOLOGIC );
                             } ) );
        demos.add( new Demo( "Tanglegram (Cophylogeny)",
                             "Two trees -- pocket gophers vs their chewing lice -- linked by an association file and "
                                     + "compared side by side.",
                             DemoTrees::openTanglegram ) );
        return demos;
    }

    /** Parse a bundled phyloXML demo tree from the jar. */
    static Phylogeny loadTree( final String resource ) throws IOException {
        try ( final InputStream in = DemoTrees.class.getResourceAsStream( RESOURCE_DIR + resource ) ) {
            if ( in == null ) {
                throw new IOException( "bundled demo tree not found on the classpath: " + RESOURCE_DIR + resource );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny[] phys = ParserBasedPhylogenyFactory.getInstance().create( in, parser );
            if ( ( phys == null ) || ( phys.length == 0 ) || ( phys[ 0 ] == null ) || phys[ 0 ].isEmpty() ) {
                throw new IOException( "bundled demo tree is empty: " + resource );
            }
            return phys[ 0 ];
        }
    }

    /** Read a bundled text resource (e.g. the tanglegram association TSV) from the jar. */
    static String loadText( final String resource ) throws IOException {
        try ( final InputStream in = DemoTrees.class.getResourceAsStream( RESOURCE_DIR + resource ) ) {
            if ( in == null ) {
                throw new IOException( "bundled demo resource not found on the classpath: " + RESOURCE_DIR + resource );
            }
            return new String( in.readAllBytes(), StandardCharsets.UTF_8 );
        }
    }

    /** Open a bundled demo tree in a new tab (the same load path as File -> Open, so the per-load auto-detect fires),
     *  and return the newly-current TreePanel for post-load configuration. */
    private static TreePanel openTree( final MainFrameApplication mf, final String resource ) throws IOException {
        final Phylogeny phy = loadTree( resource );
        AptxUtil.addPhylogeniesToTabs( new Phylogeny[] { phy }, "", "", mf.getConfiguration(), mf.getMainPanel() );
        return mf.getMainPanel().getCurrentTreePanel();
    }

    /** Colour the current tree's tips by a property ref, reflected in the "Color by" dropdown. */
    private static void colorBy( final MainFrameApplication mf, final String ref ) {
        final ControlPanel cp = mf.getMainPanel().getControlPanel();
        if ( cp != null ) {
            cp.demoSelectColorByProperty( ref );
        }
    }

    /** Turn on ancestral-state pies for the tree's first discrete trait (robust to the exact trait name). */
    private static void firstAncestralPie( final MainFrameApplication mf, final TreePanel tp ) {
        final ControlPanel cp = mf.getMainPanel().getControlPanel();
        if ( ( cp == null ) || ( tp == null ) || ( tp.getPhylogeny() == null ) ) {
            return;
        }
        final java.util.SortedSet<String> traits = TreePanelUtil.ancestralStateTraits( tp.getPhylogeny() );
        if ( !traits.isEmpty() ) {
            cp.demoSelectAncestralPie( traits.first() );
        }
    }

    /** Show the given time axis on the current tree: force a phylogram (the axis needs branch lengths = time) and set
     *  the per-tree axis type, then re-fit so its reserved band is applied. */
    private static void timeAxis( final MainFrameApplication mf, final TreePanel tp,
                                  final Options.TIME_AXIS_TYPE type ) {
        final ControlPanel cp = mf.getMainPanel().getControlPanel();
        if ( ( cp == null ) || ( tp == null ) ) {
            return;
        }
        cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM ); // "Draw Phylogram"
        tp.setTimeAxisType( type );
        refit( mf );
    }

    /** Re-fit the current tree to the window so a newly-enabled overlay's reserved space is applied. */
    private static void refit( final MainFrameApplication mf ) {
        final ControlPanel cp = mf.getMainPanel().getControlPanel();
        if ( cp != null ) {
            cp.showWhole();
        }
    }

    /** Open the cophylogeny tanglegram demo in its own window: two trees linked by a bundled association file. Built
     *  directly from the bundled resources (no tabs needed), mirroring MainFrameApplication.createTanglegram. */
    private static void openTanglegram( final MainFrameApplication mf ) throws IOException {
        final Phylogeny host = loadTree( "tanglegram-host-tree.xml" );
        final Phylogeny parasite = loadTree( "tanglegram-parasite-tree.xml" );
        final TanglegramAssociation assoc = TanglegramAssociation.parse( loadText( "tanglegram-association.tsv" ) );
        final TanglegramFrame frame = new TanglegramFrame( host, parasite, TanglegramLinker.LinkField.NODE_NAME,
                TanglegramLinker.LinkField.NODE_NAME, assoc.leftToRight(), "Pocket gophers (hosts)",
                "Chewing lice (parasites)" );
        final TreeColorSet colorset = mf.getMainPanel().getTreeColorSet();
        if ( colorset != null ) {
            frame.getTanglegramPanel().applyThemeColors( colorset.getBackgroundColor(), colorset.getBranchColor() );
        }
        frame.setVisible( true );
    }

    private DemoTrees() {
    }
}
