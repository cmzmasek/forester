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
import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Radial (circular/unrooted) node-content parity, increment 1a: asserts that in a CIRCULAR layout the enriched tip
 * labels and the INTERNAL-node labels now render (they went through an impoverished, external-only path before), and
 * that Color-by adds colored ink radially (the tip dots + property-colored labels). Uses ON-vs-OFF ink deltas so it is
 * independent of the display-density-dependent node->device transform. Headful; a green no-op when headless.
 */
public final class RadialNodeContentRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RadialNodeContentRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return labelsOk() && dotsOk() && numbersOk() && collapseOk();
    }

    /** A collapsed clade renders in circular AND unrooted WITHOUT crashing (collapsed clade-roots are given a ring
     *  angle + coords -- reading their absent angle used to NPE in circular), AND unrooted now HONORS collapse (its
     *  hidden subtree is not drawn -- paintUnrooted used to recurse through collapsed clades, drawing them expanded). */
    private static boolean collapseOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 760, h = 760;
            o.setGraphicsExportWhiteBackground( true );
            frame.showWhole();
            tp.setSize( w, h );
            // pick the LARGEST non-root internal clade so hiding it produces an unambiguous ink drop
            org.forester.phylogeny.PhylogenyNode target = null;
            for ( final java.util.Iterator<org.forester.phylogeny.PhylogenyNode> it =
                    tp.getPhylogeny().iteratorPreorder(); it.hasNext(); ) {
                final org.forester.phylogeny.PhylogenyNode n = it.next();
                if ( !n.isExternal() && !n.isRoot() && ( n.getNumberOfExternalNodes() >= 2 )
                        && ( ( target == null ) || ( n.getNumberOfExternalNodes() > target.getNumberOfExternalNodes() ) ) ) {
                    target = n;
                }
            }
            if ( ( target == null ) || target.getAllExternalDescendants().isEmpty() ) {
                fail( ok, "precondition: expected an internal clade to collapse" );
                return;
            }
            // lay the tree out UNROOTED with the clade EXPANDED, and capture a hidden-to-be descendant's coords
            final org.forester.phylogeny.PhylogenyNode hidden = target.getAllExternalDescendants().get( 0 );
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
            tp.calcParametersForPainting( w, h );
            AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
            final float hx = hidden.getXcoord(), hy = hidden.getYcoord();
            // collapse must be done in a NON-unrooted layout (the app refuses "Cannot collapse in unrooted display
            // type"); the user collapses in rectangular/circular, then VIEWS unrooted
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.calcParametersForPainting( w, h );
            tp.collapse( target );
            // renders in both radial layouts (the withFrame wrapper turns any thrown exception -- e.g. the old circular
            // collapse NPE -- into a failure) and still draws content
            for ( final Options.PHYLOGENY_GRAPHICS_TYPE gt : new Options.PHYLOGENY_GRAPHICS_TYPE[] {
                    Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                tp.setPhylogenyGraphicsType( gt );
                tp.calcParametersForPainting( w, h );
                if ( countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) ) < 200 ) {
                    fail( ok, "a collapsed-clade tree must still render its branches/labels in " + gt );
                }
            }
            // UNROOTED now HONORS collapse: paintUnrooted stops at the collapsed clade, so its hidden descendants are
            // NOT re-laid-out -- the hidden leaf keeps its pre-collapse coords (without the fix it recurses and
            // re-positions them for the reflowed layout).
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
            tp.calcParametersForPainting( w, h );
            AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
            if ( ( hidden.getXcoord() != hx ) || ( hidden.getYcoord() != hy ) ) {
                fail( ok, "unrooted must not lay out a collapsed clade's hidden descendants (leaf moved from (" + hx
                        + "," + hy + ") to (" + hidden.getXcoord() + "," + hidden.getYcoord() + "))" );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** Support (confidence) + branch-length NUMBERS render on the branches in circular AND unrooted layouts
     *  (root-on-top.xml carries bootstrap support on internals + branch lengths on every branch). */
    private static boolean numbersOk() {
        final boolean[] ok = { true };
        withFrame( "root-on-top.xml", ( frame, tp, o ) -> {
            final ControlPanel cp = tp.getControlPanel();
            final int w = 840, h = 840;
            o.setGraphicsExportWhiteBackground( true );
            frame.showWhole();
            tp.setSize( w, h );
            for ( final Options.PHYLOGENY_GRAPHICS_TYPE gt : new Options.PHYLOGENY_GRAPHICS_TYPE[] {
                    Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                tp.setPhylogenyGraphicsType( gt );
                cp.setCheckbox( Configuration.write_confidence_values, false );
                cp.setCheckbox( Configuration.write_branch_length_values, false );
                tp.calcParametersForPainting( w, h );
                final int no_num = countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                cp.setCheckbox( Configuration.write_confidence_values, true );
                cp.setCheckbox( Configuration.write_branch_length_values, true );
                tp.calcParametersForPainting( w, h );
                final int with_num = countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                if ( with_num <= ( no_num + 150 ) ) {
                    fail( ok, "support + branch-length numbers must render on the branches in " + gt + " (dark ink "
                            + with_num + " vs " + no_num + ")" );
                }
            }
        }, ok );
        return ok[ 0 ];
    }

    /** Internal-node labels + enriched tip labels render in a circular layout (colorize-by-rank.xml: tips carry an
     *  'order' taxonomy + a species node name; internal clade roots carry the order too). */
    private static boolean labelsOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final ControlPanel cp = tp.getControlPanel();
            final int w = 820, h = 820;
            o.setGraphicsExportWhiteBackground( true );
            frame.showWhole();
            tp.setSize( w, h );
            cp.setCheckbox( Configuration.show_taxonomy_scientific_names, true );
            cp.setCheckbox( Configuration.show_tax_rank, true );
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );

            // TIP labels render radially: showing taxonomy adds dark text ink over a no-labels baseline
            cp.setCheckbox( Configuration.show_taxonomy_scientific_names, false );
            cp.setCheckbox( Configuration.display_node_data, false );
            tp.calcParametersForPainting( w, h );
            final int no_labels = countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            cp.setCheckbox( Configuration.show_taxonomy_scientific_names, true );
            cp.setCheckbox( Configuration.display_node_data, true );
            cp.setCheckbox( Configuration.display_internal_data, false );
            tp.calcParametersForPainting( w, h );
            final int tips_only = countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            if ( tips_only <= ( no_labels + 200 ) ) {
                fail( ok, "tip labels must render in a circular layout (dark ink " + tips_only + " vs " + no_labels
                        + ")" );
            }

            // INTERNAL-node labels render radially: turning "Show Internal Data" on adds MORE dark ink (the clade
            // roots' "[order] <taxon>" labels) over the tips-only render
            cp.setCheckbox( Configuration.display_internal_data, true );
            tp.calcParametersForPainting( w, h );
            final int with_internal = countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            if ( with_internal <= ( tips_only + 100 ) ) {
                fail( ok, "internal-node labels must render radially with Show Internal Data on (dark ink "
                        + with_internal + " vs tips-only " + tips_only + ")" );
            }

            // UNROOTED: internal-node labels ride the branch there too (added dark ink over internal-data-off)
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
            cp.setCheckbox( Configuration.display_internal_data, false );
            tp.calcParametersForPainting( w, h );
            final int u_ext = countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            cp.setCheckbox( Configuration.display_internal_data, true );
            tp.calcParametersForPainting( w, h );
            final int u_all = countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            if ( u_all <= ( u_ext + 100 ) ) {
                fail( ok, "internal-node labels must render in an UNROOTED layout too (dark ink " + u_all + " vs "
                        + u_ext + ")" );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** Color-by adds colored ink in a circular layout (the tip dots + property-colored labels), where it was a no-op
     *  before. color-by-property.xml carries a categorical 'data:host'. */
    private static boolean dotsOk() {
        final boolean[] ok = { true };
        withFrame( "color-by-property.xml", ( frame, tp, o ) -> {
            final int w = 780, h = 780;
            o.setGraphicsExportWhiteBackground( true );
            frame.showWhole();
            tp.setSize( w, h );
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.calcParametersForPainting( w, h );
            final int off = countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            tp.setColorByPropertyRef( "data:host" );
            tp.getControlPanel().populateColorByPropertyBox();
            if ( !tp.isColorByProperty() ) {
                fail( ok, "precondition: data:host should be colorable" );
                return;
            }
            tp.calcParametersForPainting( w, h );
            final int on = countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            if ( on <= ( off + 300 ) ) {
                fail( ok, "Color-by must add colored ink (dots + labels) in a circular layout (on=" + on + " off="
                        + off + ")" );
            }
        }, ok );
        return ok[ 0 ];
    }

    private interface FrameBody {
        void run( MainFrame frame, TreePanel tp, Options o ) throws Exception;
    }

    private static void withFrame( final String demo, final FrameBody body, final boolean[] ok ) {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
            if ( !file.exists() ) {
                fail( ok, "demo tree missing: " + file.getAbsolutePath() );
                return;
            }
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                    .create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "radial" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    body.run( frame, frame.getMainPanel().getCurrentTreePanel(), frame.getOptions() );
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
    }

    /** Count of near-black pixels (text + branches) -- an ON-vs-OFF delta isolates newly-drawn label text. */
    private static int countDark( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                if ( ( ( ( rgb >> 16 ) & 0xFF ) < 90 ) && ( ( ( rgb >> 8 ) & 0xFF ) < 90 )
                        && ( ( rgb & 0xFF ) < 90 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Count of vividly-chromatic pixels (property colors), excluding white/black/gray. */
    private static int countVivid( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final float[] hsb = java.awt.Color.RGBtoHSB( ( rgb >> 16 ) & 0xFF, ( rgb >> 8 ) & 0xFF, rgb & 0xFF,
                        null );
                if ( ( hsb[ 1 ] >= 0.35f ) && ( hsb[ 2 ] >= 0.30f ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [RadialNodeContentRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private RadialNodeContentRenderTest() {
    }
}
