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

import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Integration test for the "Colorize Subtrees via Taxonomic Rank" legend: colorizing by rank builds
 * a taxon-&gt;color legend, the legend is drawn (its on-screen bounds recorded), it shares the
 * property-legend drag machinery (move + reset), and clearing branch colors removes it. Guarded to a
 * no-op when headless or when the panel has no usable viewport. Needs FlatLaf on the classpath (via
 * {@code MainFrameApplication.createInstance}), so it runs standalone, not in the headless suite.
 */
public final class RankLegendTest {

    private static final String[] ORDERS = { "Coleoptera", "Diptera", "Hymenoptera" };

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RankLegend: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return testRankColorizeLegend() && testCladeBandLegend() && testLegendHandBack() && testSingleMemberSkip()
                && testLabelAngle();
    }

    /** Clade-bar/bracket label ANGLE (root-left): HORIZONTAL / DIAGONAL labels extend rightward and reserve more
     *  right-margin than the compact VERTICAL labels, so the labels of small clades close together don't overlap. */
    private static boolean testLabelAngle() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { orderTree() }, conf, "ang" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT ); // the angle applies here
                tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS, false, TreePanel.CLADE_LABEL_ANGLE.VERTICAL );
                final int r_vert = tp.cladeBandRightReserveForTest();
                tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS, false, TreePanel.CLADE_LABEL_ANGLE.HORIZONTAL );
                final int r_horiz = tp.cladeBandRightReserveForTest();
                tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS, false, TreePanel.CLADE_LABEL_ANGLE.DIAGONAL );
                final int r_diag = tp.cladeBandRightReserveForTest();
                // "Coleoptera"/"Diptera"/"Hymenoptera" are much wider than one line, so HORIZONTAL reserves far more
                if ( r_horiz <= ( r_vert + 20 ) ) {
                    ok[ 0 ] = false;
                    System.out.println( "  HORIZONTAL labels must reserve more right-margin than VERTICAL: vert="
                            + r_vert + " horiz=" + r_horiz );
                }
                if ( ( r_diag <= r_vert ) || ( r_diag > r_horiz ) ) {
                    ok[ 0 ] = false;
                    System.out.println( "  DIAGONAL reserve should sit between VERTICAL and HORIZONTAL: vert=" + r_vert
                            + " diag=" + r_diag + " horiz=" + r_horiz );
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

    /** "Skip single-member clades" (Bars/Brackets, on by default): a bar/bracket over a single-tip clade is a
     *  degenerate stub, so it is not drawn (the drawn count excludes it); OFF draws it; BOXES ignore the option. */
    private static boolean testSingleMemberSkip() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { mixedSizeOrderTree() }, conf, "sm" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                // two clades at rank 'order': Grupo (3 tips, multi-member) + Solitaria (1 tip, single-member)
                final int all = tp.setCladeBands( "order", TreePanel.CLADE_VIS.BRACKETS, false );
                if ( all != 2 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  expected 2 clade bands (Grupo + Solitaria), got " + all );
                }
                final int brackets = tp.setCladeBands( "order", TreePanel.CLADE_VIS.BRACKETS, true );
                if ( brackets != 1 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  skip-singletons: only the 3-tip clade should get a bracket (1), got " + brackets );
                }
                final int bars = tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS, true );
                if ( bars != 1 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  skip-singletons applies to BARS too (expected 1), got " + bars );
                }
                final int boxes = tp.setCladeBands( "order", TreePanel.CLADE_VIS.BOXES, true );
                if ( boxes != 2 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  BOXES must ignore skip-singletons (draw all 2), got " + boxes );
                }
                // cladeBandCount() reports ALL placed clades (2) even when a bar/bracket run skips one -- the total the
                // "all placed but all single-member" message keys on (distinct from the drawn count).
                tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS, true ); // drawn 1, but 2 clades were placed
                if ( tp.cladeBandCount() != 2 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  cladeBandCount() must report all 2 placed clades, got " + tp.cladeBandCount() );
                }
                // COLLAPSE-INDEPENDENCE: collapse the 3-tip Grupo clade so its DISPLAY count becomes 1; a real multi-tip
                // clade must STILL draw its bar (single-member is the STRUCTURAL tip count, not the collapsed count), so
                // skip does NOT drop it. (Mutation: keying isSingleMemberClade on getNumberOfExternalNodes() -> 0 drawn.)
                final Phylogeny phy = tp.getPhylogeny();
                final PhylogenyNode grupo = phy.getRoot().getChildNode( 0 );
                grupo.setCollapse( true );
                phy.recalculateNumberOfExternalDescendants( true );
                if ( grupo.getNumberOfExternalNodes() != 1 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  collapse should drop the display count to 1, got " + grupo.getNumberOfExternalNodes() );
                }
                final int after_collapse = tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS, true );
                if ( after_collapse != 1 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  a COLLAPSED multi-tip clade must still draw its bar (expected 1), got " + after_collapse );
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

    /** Two clades at rank 'order': "Grupo" over 3 leaves (multi-member) and a lone "Solitaria" leaf (single-member).
     *  Each is resolvable OFFLINE (in-tree order taxonomy), so no network / no "resolve online" prompt. */
    private static Phylogeny mixedSizeOrderTree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode grupo = new PhylogenyNode();
        grupo.getNodeData().setTaxonomy( orderTaxon( "Grupo" ) );
        for( int c = 0; c < 3; ++c ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "g" + c );
            grupo.addAsChild( leaf ); // inherits the ancestor 'order' annotation
        }
        root.addAsChild( grupo );
        final PhylogenyNode solo = new PhylogenyNode(); // a single leaf carrying its OWN order (a 1-member clade)
        solo.setName( "s0" );
        solo.getNodeData().setTaxonomy( orderTaxon( "Solitaria" ) );
        root.addAsChild( solo );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        phy.recalculateNumberOfExternalDescendants( false );
        return phy;
    }

    private static Taxonomy orderTaxon( final String name ) {
        final Taxonomy tax = new Taxonomy();
        tax.setScientificName( name );
        try {
            tax.setRank( "order" );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        return tax;
    }

    /** Clade bands ("Annotate Clades by Rank") must reuse the same draggable taxon-&gt;color legend, so the
     *  boxes/bars are labeled; clearing the bands must remove it. */
    private static boolean testCladeBandLegend() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { orderTree() }, conf, "clb" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final int bands = tp.setCladeBands( "order", TreePanel.CLADE_VIS.BOXES );
                if ( bands < 1 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  no clade bands at rank 'order'" );
                }
                if ( !tp.hasRankLegend() ) {
                    ok[ 0 ] = false;
                    System.out.println( "  clade bands did not build the legend" );
                }
                final Rectangle vp = tp.getVisibleRect();
                if ( ( vp.width >= 300 ) && ( vp.height >= 300 ) ) {
                    paint( tp, vp ); // on-screen path records the legend bounds
                    if ( tp.getPropertyLegendBounds() == null ) {
                        ok[ 0 ] = false;
                        System.out.println( "  clade-band legend was not drawn (no bounds)" );
                    }
                }
                // unified click-to-recolor: a legend recolor updates the band/legend color and persists,
                // and "Use Automatic Color" reverts it
                tp.setRankColorOverride( "order", ORDERS[ 0 ], java.awt.Color.MAGENTA );
                if ( !java.awt.Color.MAGENTA.equals( tp.rankLegendColor( ORDERS[ 0 ] ) ) ) {
                    ok[ 0 ] = false;
                    System.out.println( "  clade-band legend recolor did not take effect" );
                }
                tp.clearRankColorOverride( "order", ORDERS[ 0 ] );
                if ( java.awt.Color.MAGENTA.equals( tp.rankLegendColor( ORDERS[ 0 ] ) ) ) {
                    ok[ 0 ] = false;
                    System.out.println( "  'automatic' did not revert the clade-band recolor" );
                }
                // VERTICAL PARITY: the clade boxes/bars/brackets must draw in a vertical orientation (geometry rides R,
                // labels re-anchored upright). Render BOXES in ROOT_TOP and count TINTED pixels: the translucent box
                // wash (~46 saturation) is distinct from the OPAQUE legend chips (~255) and the grayscale tree (~0), so
                // this isolates the boxes from the legend -- a disabled vertical clade paint drops it to ~0.
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                tp.setCladeBands( "order", TreePanel.CLADE_VIS.BOXES );
                if ( cladeBoxTintedPixels( tp ) < 200 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  vertical clade boxes drew (almost) no tinted pixels" );
                }
                // BARS mode exercises drawCladeBandLabel's vertical branch (upright taxon labels re-anchored to the base
                // frame). The colored bars sit below the tips and the upright labels below the bars, so grayscale ink
                // BELOW the deepest bar row is the labels; a broken label branch (absent / sideways / unrestored
                // transform) leaves that region empty.
                tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS );
                if ( cladeBarLabelInk( tp ) < 30 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  vertical clade BARS drew (almost) no upright taxon-label ink below the bars" );
                }
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.setCladeBands( "order", TreePanel.CLADE_VIS.BOXES ); // restore for the checks below
                // BRACKETS mode is monochrome -> it must NOT create a color legend
                tp.setCladeBands( "order", TreePanel.CLADE_VIS.BRACKETS );
                if ( tp.hasRankLegend() ) {
                    ok[ 0 ] = false;
                    System.out.println( "  bracket (B&W) mode must not create a color legend" );
                }
                // clearing the bands must remove the legend
                tp.clearCladeBands();
                if ( tp.hasRankLegend() ) {
                    ok[ 0 ] = false;
                    System.out.println( "  clade-band legend was not cleared with the bands" );
                }
                // ...but when a BRANCH rank-colorization still owns the legend, clearing the bands must hand it
                // back rather than leave the clade rows on screen. The bands overwrite the legend contents, so
                // keeping them would show a key at a rank the branches are NOT colored by -- and a colour picked
                // on one of those rows would be stored against the branch rank instead.
                ( (JFrame) mf[ 0 ] ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean testRankColorizeLegend() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { orderTree() }, conf, "rl" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                // #4: activate color-by-property first; a rank colorization must turn it back off so its
                // (property) legend is not drawn over the rank-colored branches
                tp.setColorByPropertyRef( "test:grp" );
                if ( !tp.isColorByProperty() ) {
                    ok[ 0 ] = false;
                    System.out.println( "  color-by-property did not activate (test setup problem)" );
                }
                final int colorized = tp.colorByRank( "order" ); // UI-free variant (no result dialog)
                if ( tp.isColorByProperty() ) {
                    ok[ 0 ] = false;
                    System.out.println( "  colorByRank did not turn off color-by-property" );
                }
                mf[ 0 ].getMainPanel().getControlPanel().showWhole();
                if ( colorized < 1 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  no subtrees colorized by rank 'order'" );
                }
                if ( !tp.hasRankLegend() ) {
                    ok[ 0 ] = false;
                    System.out.println( "  rank legend was not built" );
                }
                // each of the 3 orders covers 2 tips, so the rank legend must record a count of 2 per taxon
                for( final String order : ORDERS ) {
                    final Integer c = tp.rankLegendCount( order );
                    if ( ( c == null ) || ( c.intValue() != 2 ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  rank legend count for " + order + " expected 2, got " + c );
                    }
                }
                final Rectangle vp = tp.getVisibleRect();
                if ( ( vp.width < 300 ) || ( vp.height < 300 ) ) {
                    ( (JFrame) mf[ 0 ] ).dispose();
                    return; // no usable viewport; nothing more to assert
                }
                paint( tp, vp );
                final Rectangle home = tp.getPropertyLegendBounds();
                if ( home == null ) {
                    ok[ 0 ] = false;
                    System.out.println( "  rank legend was not drawn (no bounds)" );
                }
                else {
                    // a legend row maps to one of the colorized taxa
                    final Set<String> expected = new HashSet<>( Arrays.asList( ORDERS ) );
                    String found = null;
                    for( int yy = home.y + 1; yy < ( home.y + home.height ); ++yy ) {
                        final String v = tp.legendValueAt( at( tp, home.x + 12, yy ) );
                        if ( v != null ) {
                            found = v;
                            break;
                        }
                    }
                    if ( ( found == null ) || !expected.contains( found ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  legend row did not map to a colorized taxon: " + found );
                    }
                    // a click in the title row must NOT map to a value (guards the negative-index-truncates-to-0 bug)
                    if ( tp.legendValueAt( at( tp, home.x + 12, home.y + 8 ) ) != null ) {
                        ok[ 0 ] = false;
                        System.out.println( "  title-row click wrongly mapped to a value row" );
                    }
                    // hit-test inside the legend vs. the opposite corner
                    if ( !tp.isOnPropertyLegend( at( tp, home.x + ( home.width / 2 ), home.y + 3 ) )
                            || tp.isOnPropertyLegend( at( tp, vp.x + 5, vp.y + 5 ) ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  legend hit-test failed" );
                    }
                    // drag it left 80 / down 120; the box must follow
                    tp.startLegendDrag( at( tp, home.x + 10, home.y + 5 ) );
                    tp.dragLegend( at( tp, home.x + 10 - 80, home.y + 5 + 120 ) );
                    tp.endLegendDrag();
                    paint( tp, vp );
                    final Rectangle moved = tp.getPropertyLegendBounds();
                    if ( ( Math.abs( moved.x - ( home.x - 80 ) ) > 1 ) || ( Math.abs( moved.y - ( home.y + 120 ) ) > 1 ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  legend did not follow the drag: " + moved + " vs home " + home );
                    }
                    // reset returns it to the default corner
                    tp.resetLegendPosition();
                    paint( tp, vp );
                    final Rectangle back = tp.getPropertyLegendBounds();
                    if ( ( Math.abs( back.x - home.x ) > 1 ) || ( Math.abs( back.y - home.y ) > 1 ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  legend did not reset to its corner" );
                    }
                    // unified click-to-recolor also works on the branch rank-colorize legend (re-applies
                    // to the branches), and "Use Automatic Color" reverts it
                    if ( found != null ) {
                        tp.setRankColorOverride( "order", found, java.awt.Color.PINK );
                        if ( !java.awt.Color.PINK.equals( tp.rankLegendColor( found ) ) ) {
                            ok[ 0 ] = false;
                            System.out.println( "  rank-colorize legend recolor did not take effect" );
                        }
                        tp.clearRankColorOverride( "order", found );
                        if ( java.awt.Color.PINK.equals( tp.rankLegendColor( found ) ) ) {
                            ok[ 0 ] = false;
                            System.out.println( "  'automatic' did not revert the rank-colorize recolor" );
                        }
                    }
                    // clearing branch colors removes the legend
                    tp.clearRankLegend();
                    if ( tp.hasRankLegend() ) {
                        ok[ 0 ] = false;
                        System.out.println( "  rank legend was not cleared" );
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

    private static void paint( final TreePanel tp, final Rectangle vp ) {
        final BufferedImage img = new BufferedImage( Math.max( 1, vp.x + vp.width ), Math.max( 1, vp.y + vp.height ),
                                                     BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        tp.paintPhylogeny( g, false, false, 0, 0, 0, 0 ); // on-screen path: records the legend bounds
        g.dispose();
    }

    /** Fits the tree to a fixed size, paints it, and counts TINTED pixels (a translucent clade-box wash, distinct
     *  from the opaque legend chips and the grayscale tree). */
    private static int cladeBoxTintedPixels( final TreePanel tp ) {
        final int w = 640, h = 520;
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        tp.resetPreferredSize();
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( java.awt.Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.paintPhylogeny( g, false, false, 0, 0, 0, 0 );
        g.dispose();
        int n = 0;
        for( int x = 0; x < w; ++x ) {
            for( int y = 0; y < h; ++y ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, gc = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int d = Math.max( r, Math.max( gc, b ) ) - Math.min( r, Math.min( gc, b ) );
                if ( ( d > 25 ) && ( d < 120 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Renders CLADE_VIS.BARS in a vertical orientation and counts GRAYSCALE non-white ink BELOW the deepest colored
     *  bar row -- i.e. the upright taxon labels drawn by drawCladeBandLabel's vertical branch (the bars are saturated
     *  and above; the tree/tip-labels are above the bars). Isolates the labels from the bars AND the legend. */
    private static int cladeBarLabelInk( final TreePanel tp ) {
        final int w = 640, h = 520;
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        tp.resetPreferredSize();
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( java.awt.Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.paintPhylogeny( g, false, false, 0, 0, 0, 0 );
        g.dispose();
        int max_bar_y = -1; // deepest row that is mostly a colored BAR (many saturated px, unlike the small legend chips)
        for( int y = 0; y < h; ++y ) {
            int sat = 0;
            for( int x = 0; x < w; ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, gc = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( Math.max( r, Math.max( gc, b ) ) - Math.min( r, Math.min( gc, b ) ) ) > 60 ) {
                    ++sat;
                }
            }
            if ( sat > 40 ) {
                max_bar_y = y;
            }
        }
        if ( max_bar_y < 0 ) {
            return 0; // no bars found
        }
        int n = 0; // grayscale non-white ink below the bars = the upright taxon labels
        for( int y = max_bar_y + 2; y < h; ++y ) {
            for( int x = 0; x < w; ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, gc = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final boolean grayscale = ( Math.max( r, Math.max( gc, b ) ) - Math.min( r, Math.min( gc, b ) ) ) < 40;
                if ( ( ( rgb & 0xFFFFFF ) != 0xFFFFFF ) && grayscale ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static MouseEvent at( final TreePanel tp, final int x, final int y ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_PRESSED, System.currentTimeMillis(), 0, x, y, 1, false );
    }

    /** Three internal "order" nodes, each over two leaves, so colorizing by rank "order" yields 3 colors. */
    /**
     * Removing the clade marks while a BRANCH rank-colorization is active must hand the colour key back to that
     * colorization, not leave the clade rows on screen.
     * <p>
     * The bands do not merely add a legend, they OVERWRITE it ({@code updateCladeBandLegend} replaces the rows and
     * the title). So a naive "keep the legend if something else owns it" leaves a key describing marks that are
     * gone, at a rank the branches are not coloured by -- and because {@code currentRankLegendRank()} falls back to
     * the branch rank once the bands are gone, a colour picked on one of those stale rows is stored against the
     * WRONG rank. Unreachable until the "stop drawing the clade marks" entry gave {@code clearCladeBands()} its
     * first GUI caller.
     */
    private static boolean testLegendHandBack() {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { twoRankTree() }, conf, "handback" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                if ( tp.colorByRank( "class" ) < 1 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  precondition: the fixture must colorize branches at rank 'class'" );
                    return;
                }
                if ( !"Taxonomy: class".equals( tp.rankLegendTitleForTest() ) ) {
                    ok[ 0 ] = false;
                    System.out.println( "  precondition: branch colorization should own a 'class' legend, got "
                            + tp.rankLegendTitleForTest() );
                    return;
                }
                if ( tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS ) < 1 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  precondition: the fixture must place clade bands at rank 'order'" );
                    return;
                }
                if ( !"Taxonomy: order".equals( tp.rankLegendTitleForTest() ) ) {
                    ok[ 0 ] = false;
                    System.out.println( "  precondition: the clade bands should have taken the legend over, got "
                            + tp.rankLegendTitleForTest() );
                    return;
                }
                tp.clearCladeBands();
                if ( !tp.hasRankLegend() ) {
                    ok[ 0 ] = false;
                    System.out.println( "  the branch colorization still needs its legend after the marks go" );
                }
                else if ( !"Taxonomy: class".equals( tp.rankLegendTitleForTest() ) ) {
                    ok[ 0 ] = false;
                    System.out.println( "  removing the clade marks must hand the legend back to the BRANCH rank, "
                            + "got " + tp.rankLegendTitleForTest() );
                }
            } );
            SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Two resolvable ranks: orders nested inside classes, so one rank can colour the branches while the other
     *  draws the clade marks -- which is what makes the legend hand-back observable at all. */
    private static Phylogeny twoRankTree() {
        final PhylogenyNode root = new PhylogenyNode();
        int id = 0;
        final String[][] classes = { { "Mammalia", "Chiroptera", "Rodentia" },
                                     { "Aves", "Passeriformes", "Falconiformes" } };
        for( final String[] grp : classes ) {
            final PhylogenyNode class_node = new PhylogenyNode();
            class_node.getNodeData().setTaxonomy( rankedTaxon( grp[ 0 ], "class" ) );
            for( int o = 1; o < grp.length; ++o ) {
                final PhylogenyNode order_node = new PhylogenyNode();
                order_node.getNodeData().setTaxonomy( rankedTaxon( grp[ o ], "order" ) );
                for( int c = 0; c < 2; ++c ) {
                    final PhylogenyNode leaf = new PhylogenyNode();
                    leaf.setName( "n" + ( id++ ) );
                    order_node.addAsChild( leaf );
                }
                class_node.addAsChild( order_node );
            }
            root.addAsChild( class_node );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static Taxonomy rankedTaxon( final String name, final String rank ) {
        final Taxonomy tax = new Taxonomy();
        tax.setScientificName( name );
        try {
            tax.setRank( rank ); // validated against the controlled rank vocabulary
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        return tax;
    }

    private static Phylogeny orderTree() {
        final PhylogenyNode root = new PhylogenyNode();
        int id = 0;
        for( final String order : ORDERS ) {
            final PhylogenyNode internal = new PhylogenyNode();
            final Taxonomy tax = new Taxonomy();
            tax.setScientificName( order );
            try {
                tax.setRank( "order" ); // validated against the controlled rank vocabulary
            }
            catch ( final Exception e ) {
                throw new RuntimeException( e );
            }
            internal.getNodeData().setTaxonomy( tax );
            for( int c = 0; c < 2; ++c ) {
                final PhylogenyNode leaf = new PhylogenyNode();
                leaf.setName( "n" + ( id++ ) );
                // a colorable property so the test can activate color-by-property (see the #4 assertion)
                final PropertiesList pl = new PropertiesList();
                pl.addProperty( new Property( "test:grp", order, "", "xsd:string", Property.AppliesTo.NODE ) );
                leaf.getNodeData().setProperties( pl );
                internal.addAsChild( leaf );
            }
            root.addAsChild( internal );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private RankLegendTest() {
    }
}
