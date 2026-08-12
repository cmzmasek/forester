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
        return labelsOk() && dotsOk() && numbersOk() && collapseOk() && collapseMarkerOk() && collapseOverviewOk()
                && circularCenteredOk() && circularPhylogramOk() && overviewCollapseLayoutOk()
                && circularCladeBandsOk() && circularAnnotationRingsOk() && circularZebraWedgesOk();
    }

    /** Zebra stripes render as alternating angular WEDGES in the circular layout (the radial analogue of the
     *  rectangular alternating row bands): turning zebra on shades every other tip's slice a faint grey, so a lot of
     *  faint-grey pixels appear over the (otherwise white) figure. */
    private static boolean circularZebraWedgesOk() {
        final boolean[] ok = { true };
        withFrame( "zebra-stripes.xml", ( frame, tp, o ) -> {
            final int w = 820, h = 820;
            o.setGraphicsExportWhiteBackground( true );
            o.setShowOverview( false );
            tp.setOvOn( false );
            frame.showWhole();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.setPreferredSize( new java.awt.Dimension( w, h ) );
            tp.setSize( w, h );
            o.setShowZebraStripes( false );
            tp.calcParametersForPainting( w, h );
            final int off = countFaintGray( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            o.setShowZebraStripes( true );
            tp.calcParametersForPainting( w, h );
            final int on = countFaintGray( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            if ( on <= ( off + 5000 ) ) {
                fail( ok, "circular zebra wedges must add faint shaded wedges (faint-grey on " + on + " vs off " + off
                        + ")" );
            }
            o.setShowZebraStripes( false );
        }, ok );
        return ok[ 0 ];
    }

    /** Count of very-light near-neutral GREY pixels (all channels in [230,248] and near-equal) -- the faint
     *  zebra-over-white shade (~239); pure white (255) and darker branch/label ink fall outside the band. */
    private static int countFaintGray( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( r >= 230 ) && ( r <= 248 ) && ( Math.abs( r - g ) <= 4 ) && ( Math.abs( g - b ) <= 4 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Annotation columns render as concentric RINGS in the circular layout (the radial analogue of the rectangular
     *  tip-aligned columns): adding colour-strip + heat-map columns adds a lot of coloured ink over the columns-off
     *  baseline, AND a click on a ring focuses that column's colour key (the polar analogue of a header click). */
    private static boolean circularAnnotationRingsOk() {
        final boolean[] ok = { true };
        withFrame( "annotation-columns.xml", ( frame, tp, o ) -> {
            final int w = 900, h = 900;
            o.setGraphicsExportWhiteBackground( true );
            o.setShowOverview( false );
            tp.setOvOn( false );
            frame.showWhole();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.setPreferredSize( new java.awt.Dimension( w, h ) );
            tp.setSize( w, h );
            tp.clearAnnotationColumns();
            tp.calcParametersForPainting( w, h );
            final int v_off = countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            final java.util.List<AnnotationColumns.ColumnSpec> specs = new java.util.ArrayList<>();
            specs.add( new AnnotationColumns.ColumnSpec( "data:host", AnnotationColumns.Type.COLOR_STRIP ) );
            specs.add( new AnnotationColumns.ColumnSpec( "data:segment", AnnotationColumns.Type.COLOR_STRIP ) );
            specs.add( new AnnotationColumns.ColumnSpec( "data:viral_load", AnnotationColumns.Type.HEATMAP ) );
            tp.setAnnotationColumns( specs );
            tp.calcParametersForPainting( w, h );
            final int v_on = countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            if ( v_on <= ( v_off + 2000 ) ) {
                fail( ok, "circular annotation rings must add coloured ink (vivid on " + v_on + " vs off " + v_off
                        + ")" );
            }
            // a click ON a ring focuses that column's legend -- sweep the 3-o'clock spoke outward until one lands in a
            // legend-bearing ring band (radius-only hit test, so any tip angle works); exercises circularAnnotationRingAt
            boolean focused = false;
            for ( int r = 1; r < ( Math.min( w, h ) / 2 ); r += 2 ) {
                final java.awt.event.MouseEvent click = new java.awt.event.MouseEvent( tp,
                        java.awt.event.MouseEvent.MOUSE_CLICKED, 0L, 0, ( w / 2 ) + r, h / 2, 1, false );
                if ( tp.handleAnnotationHeaderClick( click ) && tp.hasFocusedAnnotationColumn() ) {
                    focused = true;
                    break;
                }
            }
            if ( !focused ) {
                fail( ok, "a click on an annotation ring must focus that column's legend (circular)" );
            }
            tp.clearAnnotationColumns();
        }, ok );
        return ok[ 0 ];
    }

    /** Clade bands render as POLAR arcs/sectors in the circular layout (the radial analogue of the rectangular clade
     *  wash). Uses the BARS mode (solid clade colours = a robust vivid-ink signal): turning the bands on adds a lot of
     *  coloured ink over the bands-off baseline. */
    private static boolean circularCladeBandsOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 820, h = 820;
            o.setGraphicsExportWhiteBackground( true );
            o.setShowOverview( false );
            tp.setOvOn( false );
            frame.showWhole();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.setPreferredSize( new java.awt.Dimension( w, h ) );
            tp.setSize( w, h );
            tp.clearCladeBands();
            tp.calcParametersForPainting( w, h );
            final int v_off = countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            final int n = tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS );
            if ( n < 2 ) {
                fail( ok, "precondition: expected >= 2 clade bands for rank 'order' (got " + n + ")" );
                return;
            }
            tp.calcParametersForPainting( w, h );
            final int v_on = countVivid( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
            if ( v_on <= ( v_off + 800 ) ) {
                fail( ok, "circular clade bands (bars) must add coloured arcs (vivid ink on " + v_on + " vs off "
                        + v_off + ")" );
            }
            tp.clearCladeBands();
        }, ok );
        return ok[ 0 ];
    }

    /** The circular OVERVIEW thumbnail HONORS collapse: it reuses the main paint's displayed-tip angles and SKIPS tips
     *  hidden under a collapse, instead of spreading every structural tip by the full external count (which left a
     *  collapsed thumbnail out of scale). Checked via coord-stability: a hidden tip's thumbnail coords must not be
     *  re-laid-out after its clade is collapsed. */
    private static boolean overviewCollapseLayoutOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 1200, h = 1200;
            o.setShowOverview( true );
            frame.showWhole();
            tp.setSize( w, h ); // larger than the frame viewport -> the overview turns on
            org.forester.phylogeny.PhylogenyNode target = null;
            for ( final java.util.Iterator<org.forester.phylogeny.PhylogenyNode> it =
                    tp.getPhylogeny().iteratorPreorder(); it.hasNext(); ) {
                final org.forester.phylogeny.PhylogenyNode n = it.next();
                if ( !n.isExternal() && !n.isRoot() && ( n.getNumberOfExternalNodes() >= 2 ) && ( ( target == null )
                        || ( n.getNumberOfExternalNodes() > target.getNumberOfExternalNodes() ) ) ) {
                    target = n;
                }
            }
            if ( ( target == null ) || target.getAllExternalDescendants().isEmpty() ) {
                fail( ok, "precondition: expected an internal clade with tips to collapse" );
                return;
            }
            // a VISIBLE tip -- one NOT under the collapsed clade
            final java.util.List<org.forester.phylogeny.PhylogenyNode> hidden_tips = target.getAllExternalDescendants();
            org.forester.phylogeny.PhylogenyNode vis = null;
            for ( final org.forester.phylogeny.PhylogenyNode n : tp.getPhylogeny().getExternalNodes() ) {
                if ( !hidden_tips.contains( n ) ) {
                    vis = n;
                    break;
                }
            }
            if ( vis == null ) {
                fail( ok, "precondition: need a tip outside the collapsed clade" );
                return;
            }
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.collapse( target );
            tp.calcParametersForPainting( w, h );
            paintScreen( tp, w, h ); // main circular paint + the overview thumbnail
            if ( !tp.isOvOn() ) {
                fail( ok, "precondition: the overview must be ON so the thumbnail layout is exercised" );
                return;
            }
            // the overview reuses the main paint's DISPLAYED-tip angles, so a visible tip sits at the SAME angle in the
            // thumbnail as in the main tree. The old full-external-count spacing would put it at a different angle once
            // a clade is collapsed (displayed count < full count), leaving the thumbnail out of scale.
            final org.forester.phylogeny.PhylogenyNode root = tp.getPhylogeny().getRoot();
            final double main_angle = Math.atan2( vis.getYcoord() - root.getYcoord(), vis.getXcoord() - root.getXcoord() );
            final double ov_angle = Math.atan2( vis.getYSecondary() - root.getYSecondary(),
                    vis.getXSecondary() - root.getXSecondary() );
            double d = Math.abs( main_angle - ov_angle ) % ( 2 * Math.PI );
            if ( d > Math.PI ) {
                d = ( 2 * Math.PI ) - d;
            }
            if ( d > 0.05 ) {
                fail( ok, "the circular overview must place a visible tip at the SAME angle as the main tree once a "
                        + "clade is collapsed (main " + main_angle + " vs overview " + ov_angle + " rad)" );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** Paint the panel through its SCREEN path (paintComponent), which -- unlike the export path -- draws the overview
     *  thumbnail. Discards the image; used to exercise the overview layout. */
    private static void paintScreen( final TreePanel tp, final int w, final int h ) {
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final java.awt.Graphics2D g = img.createGraphics();
        try {
            tp.paintComponent( g );
        }
        finally {
            g.dispose();
        }
    }

    /** A circular PHYLOGRAM (Draw Phylogram + branch lengths) positions each tip by its DISTANCE-to-root (not all on
     *  one outer ring like a cladogram), and draws concentric distance RINGS. Checked transform-independently: a tip's
     *  distance from the ring centre (root) must scale with its distance-to-root; and the faint light-grey ring circles
     *  must be present on the white theme. */
    private static boolean circularPhylogramOk() {
        final boolean[] ok = { true };
        withFrame( "scale-axis.xml", ( frame, tp, o ) -> {
            final int w = 820, h = 820;
            o.setGraphicsExportWhiteBackground( true ); // light theme: rings = light grey on white
            o.setShowOverview( false );
            tp.setOvOn( false );
            tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM ); // Draw Phylogram
            frame.showWhole();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.setPreferredSize( new java.awt.Dimension( w, h ) );
            tp.setSize( w, h );
            // render the phylogram with the SCALE OFF (no rings) -- positions the tips + the rings-off reference
            o.setShowScale( false );
            tp.calcParametersForPainting( w, h );
            final BufferedImage no_rings = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
            final org.forester.phylogeny.PhylogenyNode root = tp.getPhylogeny().getRoot();
            // (1) tips are placed by distance-to-root: find the farthest tip and a clearly-closer one, and assert
            // their radii-from-centre are in the same ratio as their distances-to-root (a cladogram would put both
            // on the outer ring -> ratio 1).
            org.forester.phylogeny.PhylogenyNode far = null, near = null;
            for ( final org.forester.phylogeny.PhylogenyNode n : tp.getPhylogeny().getExternalNodes() ) {
                if ( ( far == null ) || ( n.calculateDistanceToRoot() > far.calculateDistanceToRoot() ) ) {
                    far = n;
                }
            }
            for ( final org.forester.phylogeny.PhylogenyNode n : tp.getPhylogeny().getExternalNodes() ) {
                if ( ( n.calculateDistanceToRoot() < ( 0.7 * far.calculateDistanceToRoot() ) )
                        && ( ( near == null ) || ( n.calculateDistanceToRoot() < near.calculateDistanceToRoot() ) ) ) {
                    near = n;
                }
            }
            if ( ( near == null ) || ( near.calculateDistanceToRoot() <= 0 ) ) {
                fail( ok, "precondition: need tips at clearly different distances-to-root" );
                return;
            }
            final double r_far = Math.hypot( far.getXcoord() - root.getXcoord(), far.getYcoord() - root.getYcoord() );
            final double r_near = Math.hypot( near.getXcoord() - root.getXcoord(), near.getYcoord() - root.getYcoord() );
            final double expected = far.calculateDistanceToRoot() / near.calculateDistanceToRoot();
            final double actual = r_near > 0 ? ( r_far / r_near ) : 0;
            if ( Math.abs( actual - expected ) > ( 0.2 * expected ) ) {
                fail( ok, "circular phylogram tip radius must scale with distance-to-root (radius ratio " + actual
                        + " vs distance ratio " + expected + ")" );
            }
            // (2) turning the SCALE on adds the concentric distance rings on the SAME phylogram layout, so the grey
            // anti-aliasing ink from branches/labels cancels and the delta is PURELY the rings -- platform-robust
            // (no absolute threshold). The rings are gated on the "Scale" option, like the rectangular scale bar.
            o.setShowScale( true );
            tp.calcParametersForPainting( w, h );
            final BufferedImage with_rings = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
            final int g_off = countGrayRings( no_rings );
            final int g_on = countGrayRings( with_rings );
            if ( g_on <= ( g_off + 1500 ) ) {
                fail( ok, "a circular phylogram must draw concentric distance rings when the scale is shown (grey px "
                        + g_on + " with vs " + g_off + " without)" );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** The CIRCULAR ring is centred in the drawing area (so an export/view is centred + fills the canvas, not pushed
     *  into a corner). paintCircular puts the ROOT at the ring centre; on a deliberately NON-square canvas the centre
     *  must be (width/2, height/2), not the old min(w,h)/2 (which biased it toward the top-left). */
    private static boolean circularCenteredOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 1000, h = 640; // non-square on purpose: min/2 = 320 vs width/2 = 500
            o.setGraphicsExportWhiteBackground( true );
            frame.showWhole();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.setPreferredSize( new java.awt.Dimension( w, h ) ); // the circle sizes/centres from the preferred size
            tp.setSize( w, h );
            tp.calcParametersForPainting( w, h );
            AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
            final org.forester.phylogeny.PhylogenyNode root = tp.getPhylogeny().getRoot();
            if ( ( Math.abs( root.getXcoord() - ( w / 2.0 ) ) > 5 ) || ( Math.abs( root.getYcoord() - ( h / 2.0 ) ) > 5 ) ) {
                fail( ok, "the circular ring must be centred in the canvas: root (centre) at (" + root.getXcoord() + ","
                        + root.getYcoord() + "), expected ~(" + ( w / 2 ) + "," + ( h / 2 ) + ")" );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** Collapsing a clade in CIRCULAR with the OVERVIEW shown must not crash the on-screen paint: the overview's
     *  paintCircularsLite recursed into a collapsed clade's hidden INTERNAL nodes, which have no angle -> NPE (the
     *  overview counterpart of the 0.11.3 main-path fix). The export render path skips the overview, so this drives
     *  the SCREEN paint (paintComponent) with a tp larger than its viewport (so updateOvSizes turns the overview on).
     *  CIRCULAR-only: the fix is entirely in paintCircularsLite; the unrooted overview (paintUnrootedLite) never had
     *  this NPE. The collapsed clade is required to hold a hidden INTERNAL node, which is what actually triggers the
     *  recursion into the null-angle lookup -- a cherry of all-leaves would not exercise the fix. */
    private static boolean collapseOverviewOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 1200, h = 1200;
            o.setShowOverview( true );
            frame.showWhole();
            tp.setSize( w, h ); // larger than the frame viewport -> updateOvSizes() sets the overview on
            // largest non-root internal clade that CONTAINS a hidden internal node (so the NPE path is exercised)
            org.forester.phylogeny.PhylogenyNode target = null;
            for ( final java.util.Iterator<org.forester.phylogeny.PhylogenyNode> it =
                    tp.getPhylogeny().iteratorPreorder(); it.hasNext(); ) {
                final org.forester.phylogeny.PhylogenyNode n = it.next();
                if ( n.isExternal() || n.isRoot() || ( n.getNumberOfExternalNodes() < 2 ) || !hasInternalDescendant( n )
                        || ( ( target != null )
                                && ( n.getNumberOfExternalNodes() <= target.getNumberOfExternalNodes() ) ) ) {
                    continue;
                }
                target = n;
            }
            if ( target == null ) {
                fail( ok, "precondition: expected an internal clade with a hidden internal node to collapse" );
                return;
            }
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.calcParametersForPainting( w, h );
            tp.collapse( target ); // calls updateOvSizes() -> overview on (tree bigger than viewport)
            if ( !tp.isOvOn() ) {
                fail( ok, "precondition: the overview must be ON so the collapsed-circular overview paint is exercised" );
                return;
            }
            // a SCREEN paint draws the overview; the withFrame wrapper turns the old NPE into a failure
            tp.calcParametersForPainting( w, h );
            final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
            final java.awt.Graphics2D g = img.createGraphics();
            try {
                tp.paintComponent( g );
            }
            finally {
                g.dispose();
            }
        }, ok );
        return ok[ 0 ];
    }

    /** True if {@code node}'s subtree contains an internal (non-external) node other than {@code node} itself. */
    private static boolean hasInternalDescendant( final org.forester.phylogeny.PhylogenyNode node ) {
        for ( final org.forester.phylogeny.PhylogenyNode d : org.forester.phylogeny.PhylogenyMethods
                .getAllDescendants( node ) ) {
            if ( !d.isExternal() ) {
                return true;
            }
        }
        return false;
    }

    /** A collapsed clade draws a COLLAPSE MARKER (a filled triangle + "(N)" count) in circular AND unrooted -- the
     *  radial analogue of the rectangular collapsed-clade triangle, replacing the bare branch stub. The collapse-fill
     *  colour equals the branch/ink colour, so the triangle can't be isolated by colour directly. Instead the clade is
     *  rendered twice in the LIGHT theme: UNSELECTED (triangle in the black collapse-fill = dark ink) vs ALL TIPS
     *  SELECTED (triangle in the red found colour = NOT dark). The unselected render must carry MORE dark ink -- the
     *  triangle's fill -- which exercises the everyday collapse-fill path AND the found-fill response, and fails if the
     *  marker is gone (both renders then match). */
    private static boolean collapseMarkerOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 760, h = 760;
            o.setGraphicsExportWhiteBackground( true ); // LIGHT theme: collapse-fill = black (dark), found = red (not dark)
            o.setPulseFoundNodes( false );
            frame.showWhole();
            tp.setSize( w, h );
            org.forester.phylogeny.PhylogenyNode target = null;
            for ( final java.util.Iterator<org.forester.phylogeny.PhylogenyNode> it =
                    tp.getPhylogeny().iteratorPreorder(); it.hasNext(); ) {
                final org.forester.phylogeny.PhylogenyNode n = it.next();
                if ( !n.isExternal() && !n.isRoot() && ( n.getNumberOfExternalNodes() >= 2 ) && ( ( target == null )
                        || ( n.getNumberOfExternalNodes() > target.getNumberOfExternalNodes() ) ) ) {
                    target = n;
                }
            }
            if ( target == null ) {
                fail( ok, "precondition: expected an internal clade to collapse" );
                return;
            }
            final java.util.Set<Long> found = new java.util.HashSet<Long>();
            for ( final org.forester.phylogeny.PhylogenyNode t : target.getAllExternalDescendants() ) {
                found.add( t.getId() );
            }
            // collapse must be done in a non-unrooted layout (the app refuses it in unrooted)
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.calcParametersForPainting( w, h );
            tp.collapse( target );
            for ( final Options.PHYLOGENY_GRAPHICS_TYPE gt : new Options.PHYLOGENY_GRAPHICS_TYPE[] {
                    Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                tp.setPhylogenyGraphicsType( gt );
                tp.calcParametersForPainting( w, h );
                tp.setFoundNodes0( null );
                final int d_unsel = countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                tp.setFoundNodes0( found );
                final int d_found = countDark( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                tp.setFoundNodes0( null );
                if ( d_unsel <= ( d_found + 20 ) ) {
                    fail( ok, "the collapse-marker triangle must render (dark collapse-fill when unselected, found "
                            + "colour when selected) in " + gt + " (dark unsel=" + d_unsel + " found=" + d_found + ")" );
                }
            }
        }, ok );
        return ok[ 0 ];
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

    /** Count of light-medium GREY pixels (all channels in [175,235] and near-equal) -- the faint distance-ring
     *  circles on a white theme; branches/labels (near-black) and the white background fall outside the band. */
    private static int countGrayRings( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( r >= 175 ) && ( r <= 235 ) && ( g >= 175 ) && ( g <= 235 ) && ( b >= 175 ) && ( b <= 235 )
                        && ( Math.abs( r - g ) <= 8 ) && ( Math.abs( g - b ) <= 8 ) ) {
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
