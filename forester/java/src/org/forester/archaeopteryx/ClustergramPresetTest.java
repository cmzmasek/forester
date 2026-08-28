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
import java.io.File;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The one-click "Clustergram" preset (View menu). Pure part: {@code MainFrame.clustergramColumnSpecs} maps a tree's
 * numeric fields to a shared-scale heat-map MATRIX (first) and its categorical fields to color strips. Headful part:
 * {@code applyClustergramPreset} on the heat-map-matrix demo sets a rectangular, root-at-top, tip-aligned,
 * labels-below layout with the numeric columns added as a matrix. A green no-op when headless.
 */
public final class ClustergramPresetTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ClustergramPreset: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( !columnSpecsOk() ) {
            return false; // pure, headless-safe
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return presetOk();
    }

    /** Pure: numeric fields become a MATRIX (before) categorical color strips; a data-less tree yields no columns. */
    private static boolean columnSpecsOk() {
        final Phylogeny phy = mixedTree();
        final List<AnnotationColumns.ColumnSpec> specs = MainFrame.clustergramColumnSpecs( phy );
        if ( specs.size() != 2 ) {
            return fail( "a tree with one numeric + one categorical field should yield 2 columns, got " + specs.size() );
        }
        if ( ( specs.get( 0 )._type != AnnotationColumns.Type.MATRIX ) || !"data:score".equals( specs.get( 0 )._ref ) ) {
            return fail( "the numeric field should be a MATRIX column, first: " + specs.get( 0 )._ref );
        }
        if ( ( specs.get( 1 )._type != AnnotationColumns.Type.COLOR_STRIP )
                || !"data:host".equals( specs.get( 1 )._ref ) ) {
            return fail( "the categorical field should be a COLOR_STRIP column: " + specs.get( 1 )._ref );
        }
        if ( !MainFrame.clustergramColumnSpecs( barePhy() ).isEmpty() ) {
            return fail( "a tree with no per-tip data should yield no clustergram columns" );
        }
        return true;
    }

    /** Headful: the one click flips the whole layout and adds the matrix columns, and it renders as a colored figure. */
    private static boolean presetOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/heatmap-matrix.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                    .create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(),
                            "cg" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    // precondition: a plain horizontal tree with no columns
                    if ( tp.hasAnnotationColumns() || tp.isVerticalOrientation() ) {
                        fail( ok, "precondition: the tree should start as a plain horizontal tree" );
                    }
                    frame.applyClustergramPreset();
                    final Options o = frame.getOptions();
                    if ( o.getTreeOrientation() != Options.TREE_ORIENTATION.ROOT_TOP ) {
                        fail( ok, "Clustergram should set root-at-top orientation" );
                    }
                    if ( !o.isTipLabelsBelowColumns() ) {
                        fail( ok, "Clustergram should enable Tip Labels Below Columns" );
                    }
                    if ( tp.getPhylogenyGraphicsType() != Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) {
                        fail( ok, "Clustergram should force a rectangular graphics type" );
                    }
                    if ( tp.getControlPanel().getTreeDisplayType() != Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM ) {
                        fail( ok, "Clustergram should align the tips" );
                    }
                    if ( !tp.isVerticalOrientation() || !tp.hasAnnotationColumns() ) {
                        fail( ok, "Clustergram should be a vertical tree WITH annotation columns" );
                    }
                    if ( !tp.tipLabelsBelowColumns() ) {
                        fail( ok, "the labels-below clustergram layout should be active after the preset" );
                    }
                    // The preset sets the style + orientation directly rather than through MainFrame.typeChanged,
                    // so the control panel's five-way layout row has to be re-seeded by hand -- otherwise it stays
                    // lit on root-left while the clustergram is drawn root-at-top.
                    if ( tp.getControlPanel().selectedLayoutKind() != LayoutIcon.Kind.ROOT_TOP ) {
                        fail( ok, "Clustergram must leave the control-panel layout row on root-top, got "
                                + tp.getControlPanel().selectedLayoutKind() );
                    }
                    if ( tp.isEdited() ) {
                        fail( ok, "Clustergram is display-only and must NOT mark the tree edited (would pop a save "
                                + "prompt + clear the redo stack)" );
                    }
                    // it renders as a colored clustergram (the matrix cells paint)
                    final int w = 640, h = 680;
                    o.setGraphicsExportWhiteBackground( true );
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final java.awt.image.BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1,
                            false );
                    if ( countColored( img ) < 1500 ) {
                        fail( ok, "the preset should render a colored heat-map grid, got " + countColored( img ) );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static int countColored( final java.awt.image.BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int max = Math.max( r, Math.max( g, b ) ), min = Math.min( r, Math.min( g, b ) );
                if ( ( ( max - min ) > 60 ) && ( max > 100 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    // ---- fixtures ----

    /** Four tips with a numeric field (data:score) and a categorical one (data:host); host repeats (2 of 4 distinct)
     *  so it is color-able (a per-tip-unique categorical is intentionally excluded by colorableRefs). */
    private static Phylogeny mixedTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( mixedLeaf( "a", "1", "cat" ) );
        root.addAsChild( mixedLeaf( "b", "2", "dog" ) );
        root.addAsChild( mixedLeaf( "c", "3", "cat" ) );
        root.addAsChild( mixedLeaf( "d", "4", "dog" ) );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode mixedLeaf( final String name, final String score, final String host ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( "data:score", score, "", "xsd:string", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:host", host, "", "xsd:string", AppliesTo.NODE ) );
        n.getNodeData().setProperties( pl );
        return n;
    }

    private static Phylogeny barePhy() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final String name : new String[] { "a", "b", "c" } ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( name );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ClustergramPresetTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ClustergramPresetTest] " + msg );
        ok[ 0 ] = false;
    }

    private ClustergramPresetTest() {
    }
}
