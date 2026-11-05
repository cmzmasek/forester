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

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Window;
import java.awt.event.MouseEvent;
import java.awt.geom.Area;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JRadioButton;
import javax.swing.SwingUtilities;
import javax.swing.Timer;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.Options.TREE_ORIENTATION;
import org.forester.archaeopteryx.TreePanel.CLADE_VIS;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Shaded clade boxes ("Annotate Clades by Rank", Christian 2026-10-31), on the two demos made for them:
 * <ul>
 * <li><b>Two legends.</b> The boxes' legend stays on screen beside the Color-by legend -- through Color by on, off and
 * on again, and with Size by -- in all five display types; each legend is clicked, dragged and reset on its own.
 * (Until 0.11.198 Color by took the one legend slot and wiped the boxes' legend for good.)</li>
 * <li><b>Names beside the boxes.</b> "Named beside each box" draws each clade's name past its box instead of a legend:
 * in the rectangular layout only where the box is a line tall (wide enough, root at top/bottom), so no two names
 * meet -- the crowded twin's six one-tip genera carry none, the coronavirus demo's one-tip genus does. Circular and
 * unrooted place the names largest clade first: unrooted moves a name outward where it would touch another, circular
 * leaves it out -- and only then (each name left out overlaps a drawn one). Names never sit on a tip label, never
 * leave the canvas, and travel with a saved figure.</li>
 * <li><b>Unrooted sectors.</b> Each node's equal-angle wedge tiles its parent's; distinct clades' wedges never
 * overlap; every drawn tip label of a clade lies inside its sector (the far corners within the radius, the text's
 * axis within the angle) -- also where the tips sit on zero-length branches; and the sector is painted where it is
 * claimed to be.</li>
 * <li><b>Unrooted bars and brackets</b> (Christian, 2026-11-01): no sector; at each wedge's mouth a square on the
 * tree's side and the name past it, placed by the same rules -- a bar's square filled in its clade's colour, a
 * bracket's an outline in the text colour; the one-tip genus skipped when single-member clades are.</li>
 * <li><b>The dialog</b> sets the choice, greyed for bars and brackets.</li>
 * </ul>
 * Guarded to a no-op on a headless box.
 */
public final class CladeBoxesTest {

    // the layout divides arcs in float (a tip count over a float), so its wedges tile to ~1e-7 rad: 0.005 px at 500 px
    private static final double      EPS           = 1e-5;
    private static final String      DEMO          = "forester/demo/clade-boxes.xml";
    private static final String      CROWDED       = "forester/demo/clade-boxes-crowded.xml";
    private static final Set<String> GENERA        = set( "Alphacoronavirus", "Betacoronavirus", "Gammacoronavirus",
                                                          "Deltacoronavirus", "Alphaletovirus" );
    private static final Set<String> HOSTS         = set( "Human", "Bird", "Pig", "Bat", "Camel", "Cat", "Chicken",
                                                          "Frog", "Mouse", "Turkey", "Whale" );
    private static final Set<String> CROWDED_TALL  = set( "Enterovirus", "Kobuvirus", "Aphthovirus", "Sapelovirus",
                                                          "Teschovirus", "Avihepatovirus", "Salivirus" );
    private static final Set<String> CROWDED_ALL   = set( "Enterovirus", "Hepatovirus", "Kobuvirus", "Parechovirus",
                                                          "Aphthovirus", "Cardiovirus", "Sapelovirus", "Senecavirus",
                                                          "Teschovirus", "Tremovirus", "Avihepatovirus", "Cosavirus",
                                                          "Salivirus" );
    private static final Object[][]  LAYOUTS       = {
            { "rectangular", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, TREE_ORIENTATION.ROOT_LEFT },
            { "root at top", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, TREE_ORIENTATION.ROOT_TOP },
            { "root at bottom", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, TREE_ORIENTATION.ROOT_BOTTOM },
            { "circular", PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, TREE_ORIENTATION.ROOT_LEFT },
            { "unrooted", PHYLOGENY_GRAPHICS_TYPE.UNROOTED, TREE_ORIENTATION.ROOT_LEFT } };

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "Shaded clade boxes (two legends, names, unrooted): " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            boolean ok = sectorShapes();
            ok &= sideBySide();
            ok &= withDemo( DEMO, CladeBoxesTest::onDemo );
            ok &= withDemo( CROWDED, CladeBoxesTest::onCrowded );
            ok &= withDemo( "forester/demo/domain-architectures.xml", CladeBoxesTest::groupedDomains );
            ok &= crowdedDomains();
            ok &= bigCollapsed();
            ok &= dialog();
            return ok;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    // ---- the sector shape, alone: the angle convention (y down) is where a sign error would hide -----------------

    private static boolean sectorShapes() {
        boolean ok = true;
        final double a0 = 0.3, a1 = 1.1; // below-right of the centre on screen (y grows downward)
        final Shape s = TreePanel.unrootedSectorForTest( 100, 100, 50, a0, a1 );
        ok &= ck( s.contains( 100 + ( 40 * Math.cos( 0.7 ) ), 100 + ( 40 * Math.sin( 0.7 ) ) ),
                  "a sector holds a point on its middle spoke" );
        ok &= ck( !s.contains( 100 + ( 40 * Math.cos( -0.7 ) ), 100 + ( 40 * Math.sin( -0.7 ) ) ),
                  "a sector is not its mirror image across the x axis" );
        ok &= ck( !s.contains( 100 + ( 40 * Math.cos( a1 + 0.15 ) ), 100 + ( 40 * Math.sin( a1 + 0.15 ) ) ),
                  "a sector ends at its upper angle" );
        ok &= ck( !s.contains( 100 + ( 40 * Math.cos( a0 - 0.15 ) ), 100 + ( 40 * Math.sin( a0 - 0.15 ) ) ),
                  "a sector starts at its lower angle" );
        ok &= ck( !s.contains( 100 + ( 55 * Math.cos( 0.7 ) ), 100 + ( 55 * Math.sin( 0.7 ) ) ),
                  "a sector ends at its radius" );
        // a wedge that runs past 2 pi (the layout's angles do: they start at an offset and add)
        final Shape w = TreePanel.unrootedSectorForTest( 0, 0, 50, 6.0, 6.8 );
        ok &= ck( w.contains( 40 * Math.cos( 6.4 ), 40 * Math.sin( 6.4 ) ), "a sector across 2 pi holds its middle" );
        ok &= ck( TreePanel.unrootedSectorForTest( 0, 0, 50, 1.0, ( 1.0 + ( 2 * Math.PI ) ) - 2.9e-7 )
                instanceof java.awt.geom.Ellipse2D, "the layout's float full turn (2.9e-7 rad short) is a disc, no seam" );
        final Shape disc = TreePanel.unrootedSectorForTest( 0, 0, 50, 1.0, 1.0 + ( 2 * Math.PI ) );
        for ( double a = 0; a < 6.28; a += 0.5 ) {
            ok &= ck( disc.contains( 40 * Math.cos( a ), 40 * Math.sin( a ) ), "a full turn is the whole disc (angle "
                    + a + ")" );
        }
        return ok;
    }

    /** The legend column's width for the two right-edge legends: one above the other when they fit the frame's
     *  height (each 10 px from the edge and apart), else side by side, 10 px apart. */
    private static boolean sideBySide() {
        boolean ok = true;
        final java.awt.Dimension shared = new java.awt.Dimension( 200, 300 ), rank = new java.awt.Dimension( 150, 100 );
        ok &= ck( TreePanel.sideBySideLegendsWidth( shared, rank, 420 ) == 0, "stacked: 10 + 300 + 10 + 100 fits 420" );
        ok &= ck( TreePanel.sideBySideLegendsWidth( shared, rank, 419 ) == 360, "side by side below that: 200 + 10 + "
                + "150; got " + TreePanel.sideBySideLegendsWidth( shared, rank, 419 ) );
        ok &= ck( ( TreePanel.sideBySideLegendsWidth( null, rank, 10 ) == 0 )
                && ( TreePanel.sideBySideLegendsWidth( shared, null, 10 ) == 0 ), "one legend: nothing beside it" );
        return ok;
    }

    // ---- the coronavirus demo --------------------------------------------------------------------------------

    private static boolean onDemo( final MainFrame mf ) {
        final TreePanel tp = mf.getCurrentTreePanel();
        boolean ok = true;
        // as it opens: the saved figure
        paint( tp );
        ok &= ck( ( tp.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && tp.isCladeBoxLabels()
                && "vipr:Host".equals( tp.getColorByPropertyRef() ) && ( tp.getCladeBandsMode() == CLADE_VIS.BOXES ),
                  "the demo opens unrooted, its boxes named, coloured by host; got " + tp.getPhylogenyGraphicsType()
                          + " names " + tp.isCladeBoxLabels() + " colour " + tp.getColorByPropertyRef() );
        ok &= ck( new HashSet<String>( tp.cladeBoxNamesDrawnForTest() ).equals( GENERA ),
                  "as it opens, the five genera are named; got " + tp.cladeBoxNamesDrawnForTest() );
        ok &= twoLegends( mf, tp );
        ok &= namesBesideBoxes( mf, tp );
        ok &= unrootedSectors( mf, tp );
        ok &= unrootedBarsAndBrackets( mf, tp );
        ok &= unrootedLevels( mf, tp );
        ok &= unrootedViews( mf, tp );
        ok &= figureRoundTrip( tp );
        setLayout( mf, tp, LAYOUTS[ 0 ] );
        return ok;
    }

    private static boolean twoLegends( final MainFrame mf, final TreePanel tp ) {
        boolean ok = true;
        tp.setCladeBoxLabels( false );
        ok &= ck( tp.setCladeBands( "genus", CLADE_VIS.BOXES ) == 5, "the demo places 5 genus clades" );
        Rectangle alone_home = null;
        for ( final Object[] layout : LAYOUTS ) {
            final String name = ( String ) layout[ 0 ];
            setLayout( mf, tp, layout );
            tp.setSizeByPropertyRef( null );
            tp.setColorByPropertyRef( null );
            paint( tp );
            final Rectangle alone = tp.getRankLegendBoundsForTest();
            ok &= ck( ( tp.getColorByPropertyRef() == null ) && ( alone != null ), name
                    + ": the boxes' legend is drawn with no Color by" );
            if ( layout == LAYOUTS[ 0 ] ) {
                alone_home = alone;
            }
            for ( int round = 1; round <= 2; ++round ) { // on, off, ON AGAIN: the legend must come back
                tp.setColorByPropertyRef( "vipr:Host" );
                paint( tp );
                final Rectangle shared = tp.getPropertyLegendBounds();
                final Rectangle rank = tp.getRankLegendBoundsForTest();
                ok &= ck( "vipr:Host".equals( tp.getColorByPropertyRef() ) && ( shared != null ) && ( rank != null )
                        && !shared.intersects( rank ), name + ", Color by on (" + round
                                + "): both legends drawn, apart; got " + shared + " / " + rank );
                tp.setColorByPropertyRef( null );
                paint( tp );
                ok &= ck( ( tp.getColorByPropertyRef() == null ) && ( tp.getRankLegendBoundsForTest() != null ), name
                        + ", Color by off (" + round + "): the boxes' legend stays" );
            }
            tp.setColorByPropertyRef( "vipr:Host" );
            tp.setSizeByPropertyRef( "vipr:Genome_kb" );
            paint( tp );
            final Rectangle shared = tp.getPropertyLegendBounds(), rank = tp.getRankLegendBoundsForTest(),
                    size = tp.getSizeLegendBounds();
            ok &= ck( tp.isSizeByProperty() && ( shared != null ) && ( rank != null ) && ( size != null )
                    && !shared.intersects( rank ) && !shared.intersects( size ) && !rank.intersects( size ), name
                            + ", Color by + Size by: three legends drawn, apart; got " + shared + " / " + rank + " / "
                            + size );
            tp.setSizeByPropertyRef( null );
        }
        // the legend column is reserved for the boxes' legend alone too: it sits past the boxes, not over them
        setLayout( mf, tp, LAYOUTS[ 0 ] );
        tp.setColorByPropertyRef( null );
        if ( !mf.getOptions().isReserveLegendColumn() ) {
            mf._legend_column_cbmi.doClick(); // the shipped default, through its menu item
        }
        mf.showWhole();
        paint( tp );
        ok &= ck( mf.getOptions().isReserveLegendColumn() && ( tp.getRankLegendBoundsForTest() != null )
                && ( tp.getRankLegendBoundsForTest().x >= tp.cladeBandRightEdgeForTest() ),
                  "the boxes' legend alone has its own column, past the boxes; legend " + tp.getRankLegendBoundsForTest()
                          + ", boxes end at " + tp.cladeBandRightEdgeForTest() );
        // a window too short to stack the two: the boxes' legend stands beside the Color-by legend, and the column
        // holds both -- not the wider one alone, which left the other over the tree
        final java.awt.Dimension frame_size = mf.getSize();
        tp.setColorByPropertyRef( "vipr:Host" );
        mf.setSize( 1300, 330 ); // wide enough that both stay under the column's 40% cap
        mf.validate();
        mf.showWhole();
        paint( tp );
        final Rectangle beside_shared = tp.getPropertyLegendBounds(), beside_rank = tp.getRankLegendBoundsForTest();
        final Rectangle view_now = tp.getVisibleRect();
        ok &= ck( ( beside_shared != null ) && ( beside_rank != null )
                && ( ( beside_rank.x + beside_rank.width ) <= beside_shared.x ), "precondition: in a short window the "
                        + "boxes' legend stands beside the Color-by legend; got " + beside_shared + " / " + beside_rank
                        + " in " + view_now );
        ok &= ck( ( beside_rank != null ) && ( tp.legendColumnReserve() > 0 ) && ( beside_rank.x >= ( ( view_now.x
                + view_now.width ) - tp.legendColumnReserve() ) ), "...and the legend column holds both: column "
                        + tp.legendColumnReserve() + " px at the right of " + view_now + ", the boxes' legend at "
                        + beside_rank );
        // a fixed-size export decides by ITS canvas while laid out for it
        final int[] token = tp.layoutForExportSize( 900, 2000 );
        final int during = tp.legendPlacementHeightForTest();
        tp.restoreLayoutAfterExport( token );
        ok &= ck( ( during == 2000 ) && ( tp.legendPlacementHeightForTest() == tp.getVisibleRect().height ),
                  "the legends are placed in a fixed-size export's height while laid out for it (" + during + "), the "
                          + "view's after (" + tp.legendPlacementHeightForTest() + ")" );
        mf.setSize( frame_size );
        mf.validate();
        tp.setColorByPropertyRef( null );
        mf.showWhole();
        // the displayed tree gone empty: the bands go, and their legend with them
        {
            tp.setCladeBoxLabels( false );
            tp.setCladeBands( "genus", CLADE_VIS.BOXES );
            final Phylogeny shown = tp.getPhylogeny();
            ok &= ck( tp.hasRankLegend(), "precondition: the boxes' legend" );
            tp.setTree( new Phylogeny() );
            tp.rebuildCladeBands();
            ok &= ck( !tp.hasRankLegend(), "an empty tree keeps no bands' legend" );
            tp.setTree( shown );
            tp.rebuildCladeBands();
            ok &= ck( tp.hasRankLegend(), "...and it comes back with the tree" );
        }
        // clicked, dragged and reset each on its own (rectangular)
        setLayout( mf, tp, LAYOUTS[ 0 ] );
        tp.setColorByPropertyRef( "vipr:Host" );
        paint( tp );
        final Rectangle shared = tp.getPropertyLegendBounds(), rank = tp.getRankLegendBoundsForTest();
        if ( ( shared == null ) || ( rank == null ) || ( alone_home == null ) ) {
            return ck( false, "both legends must be drawn before they are clicked" );
        }
        ok &= ck( GENERA.contains( valueIn( tp, rank ) ), "a row of the boxes' legend reads a genus, got "
                + valueIn( tp, rank ) );
        ok &= ck( HOSTS.contains( valueIn( tp, shared ) ), "a row of the Color-by legend reads a host, got "
                + valueIn( tp, shared ) );
        ok &= ck( tp.isOnPropertyLegend( at( tp, rank.x + 5, rank.y + 5 ) )
                && tp.isOnPropertyLegend( at( tp, shared.x + 5, shared.y + 5 ) ), "a click on either legend is a "
                        + "legend click" );
        // measuring a legend (the legend column asks after the draw) leaves both where the draw put them -- their
        // clicks, drags and rows read these bounds
        tp.rankLegendSizeForTest();
        tp.sharedLegendSizeForTest();
        ok &= ck( rank.equals( tp.getRankLegendBoundsForTest() ) && shared.equals( tp.getPropertyLegendBounds() )
                && GENERA.contains( valueIn( tp, rank ) ), "measuring the legends leaves them where they were drawn; "
                        + "got " + tp.getRankLegendBoundsForTest() + " / " + tp.getPropertyLegendBounds() );
        // a click on a row of the boxes' legend, with Color by on, opens its colour and recolours that genus there --
        // through the click's own route (the row's slot decides: the rank legend's rows are not Color by's values)
        ok &= ck( pickOnRow( tp, rank, new Color( 0x13579B ) ), "a click on a row of the boxes' legend with Color by "
                + "on recolours that genus in the boxes' legend" );
        drag( tp, rank, -60, 40 );
        paint( tp );
        final Rectangle rank_moved = tp.getRankLegendBoundsForTest();
        ok &= ck( ( Math.abs( rank_moved.x - ( rank.x - 60 ) ) <= 1 ) && ( Math.abs( rank_moved.y - ( rank.y + 40 ) ) <= 1 )
                && tp.getPropertyLegendBounds().getLocation().equals( shared.getLocation() ),
                  "dragging the boxes' legend moves it alone; got " + rank_moved + " (from " + rank + "), Color by "
                          + tp.getPropertyLegendBounds() );
        drag( tp, shared, -40, 30 );
        paint( tp );
        ok &= ck( ( Math.abs( tp.getPropertyLegendBounds().x - ( shared.x - 40 ) ) <= 1 )
                && tp.getRankLegendBoundsForTest().getLocation().equals( rank_moved.getLocation() ),
                  "dragging the Color-by legend moves it alone" );
        // a double-click on the boxes' legend's title, through the panel's own mouse listener, returns it to its
        // corner (the top one: the Color-by legend was dragged away)
        final Rectangle title = tp.getRankLegendBoundsForTest();
        new MouseListener( tp ).mouseClicked( new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(),
                                                              0, title.x + 12, title.y + 6, 2, false ) );
        paint( tp );
        ok &= ck( tp.getRankLegendBoundsForTest().getLocation().equals( alone_home.getLocation() ),
                  "a double-click on the boxes' legend returns it to its corner; got " + tp.getRankLegendBoundsForTest()
                          + ", home " + alone_home );
        // a paint while a click on the boxes' legend is handled (a row's colour chooser is modal, and paints run in its
        // event loop) records each legend in its own slot: the Color-by legend, dragged first, is found where it went
        final Rectangle shared_before = tp.getPropertyLegendBounds();
        final Rectangle rank_before = tp.getRankLegendBoundsForTest();
        drag( tp, shared_before, 0, 50 );
        tp.inRankSlotForTest( () -> paint( tp ) );
        ok &= ck( ( Math.abs( tp.getPropertyLegendBounds().y - ( shared_before.y + 50 ) ) <= 1 )
                && tp.getRankLegendBoundsForTest().equals( rank_before )
                && HOSTS.contains( valueIn( tp, tp.getPropertyLegendBounds() ) )
                && GENERA.contains( valueIn( tp, tp.getRankLegendBoundsForTest() ) ),
                  "a paint during a click on the boxes' legend keeps each legend in its own place; Color by "
                          + tp.getPropertyLegendBounds() + " (dragged from " + shared_before + "), boxes "
                          + tp.getRankLegendBoundsForTest() + " (was " + rank_before + ")" );
        drag( tp, tp.getRankLegendBoundsForTest(), -60, 40 ); // away from home, so the reset has something to undo
        paint( tp );
        ok &= ck( !tp.getRankLegendBoundsForTest().getLocation().equals( alone_home.getLocation() ),
                  "precondition: the boxes' legend is away from its corner before the reset" );
        tp.resetColorStateToDefaults(); // Color by off, both legends back to their corners
        paint( tp );
        ok &= ck( tp.getRankLegendBoundsForTest() != null
                && tp.getRankLegendBoundsForTest().getLocation().equals( alone_home.getLocation() ),
                  "Reset returns the boxes' legend to its corner; got " + tp.getRankLegendBoundsForTest() + ", home "
                          + alone_home );
        return ok;
    }

    private static boolean namesBesideBoxes( final MainFrame mf, final TreePanel tp ) {
        boolean ok = true;
        // circular names in the large font, as the other layouts name their boxes (and as the ring reserved them)
        {
            final java.awt.FontMetrics large = mf.getMainPanel().getTreeFontSet().getFontMetricsLarge();
            final java.awt.FontMetrics small = mf.getMainPanel().getTreeFontSet().getFontMetricsSmall();
            tp.setCladeBoxLabels( true );
            tp.setCladeBands( "genus", CLADE_VIS.BOXES );
            setLayout( mf, tp, LAYOUTS[ 3 ] );
            paint( tp );
            final int want = large.getAscent() + large.getDescent();
            ok &= ck( want != ( small.getAscent() + small.getDescent() ), "precondition: the two fonts differ" );
            for ( final Shape name : tp.cladeBoxNameShapesForTest() ) {
                ok &= ck( Math.abs( thickness( name ) - want ) < 0.5, "circular: a box's name is in the large font ("
                        + thickness( name ) + " px thick, the large font " + want + ")" );
            }
            ok &= ck( !tp.cladeBoxNameShapesForTest().isEmpty(), "circular: precondition: names drawn" );
        }
        tp.setCladeBoxLabels( true );
        tp.setCladeBands( "genus", CLADE_VIS.BOXES );
        tp.setColorByPropertyRef( "vipr:Host" );
        for ( final Object[] layout : LAYOUTS ) {
            final String name = ( String ) layout[ 0 ];
            setLayout( mf, tp, layout );
            tp.setRecordLabelBoxesForTest( true );
            paint( tp );
            ok &= ck( tp.getRankLegendBoundsForTest() == null, name + ": named boxes draw no legend of their own" );
            final List<String> names = tp.cladeBoxNamesDrawnForTest();
            final List<Shape> shapes = tp.cladeBoxNameShapesForTest();
            ok &= ck( new HashSet<String>( names ).size() == names.size(), name + ": no name twice " + names );
            ok &= ck( disjoint( shapes ), name + ": no two names meet" );
            final PHYLOGENY_GRAPHICS_TYPE type = ( PHYLOGENY_GRAPHICS_TYPE ) layout[ 1 ];
            if ( type == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) {
                final float right = tp.cladeBandRightEdgeForTest();
                for ( final Shape s : shapes ) {
                    ok &= ck( s.getBounds2D().getMinX() > right, name + ": a name starts past the boxes' edge ("
                            + s.getBounds2D().getMinX() + " vs " + right + ")" );
                }
                if ( layout[ 2 ] == TREE_ORIENTATION.ROOT_LEFT ) {
                    // the precondition that makes "all five" mean something: the one-tip genus's box is a line tall
                    ok &= ck( ( 2 * tp.getYdistance() ) >= tp.getLargeFontHeight(), name + ": precondition: a one-tip "
                            + "box is a line tall (" + ( 2 * tp.getYdistance() ) + " vs " + tp.getLargeFontHeight() + ")" );
                    ok &= ck( new HashSet<String>( names ).equals( GENERA ), name + ": every box named, the one-tip "
                            + "Alphaletovirus too; got " + names );
                    final Rectangle legend = tp.getPropertyLegendBounds();
                    final double limit = ( legend != null ) ? legend.x : tp.getWidth();
                    for ( final Shape s : shapes ) {
                        ok &= ck( s.getBounds2D().getMaxX() <= limit, name + ": a name ends before the legend's "
                                + "column / the canvas edge (" + s.getBounds2D().getMaxX() + " vs " + limit + ")" );
                    }
                }
                else {
                    // upright names are as wide as their text: each lies within its own box, or is not drawn
                    final double pitch = 2 * tp.getYdistance();
                    int refused = 0;
                    for ( final CladeBand band : tp.cladeBandsForTest() ) {
                        final int i = names.indexOf( band.getTaxon() );
                        double lo = Double.MAX_VALUE, hi = -Double.MAX_VALUE;
                        for ( final PhylogenyNode tip : band.getRoot().getAllExternalDescendants() ) {
                            lo = Math.min( lo, tip.getYcoord() );
                            hi = Math.max( hi, tip.getYcoord() );
                        }
                        if ( i < 0 ) {
                            ++refused;
                            continue;
                        }
                        final Rectangle2D b = shapes.get( i ).getBounds2D();
                        ok &= ck( ( b.getMinY() >= ( lo - ( pitch / 2 ) - 0.5 ) ) && ( b.getMaxY() <= ( hi + ( pitch / 2 )
                                + 0.5 ) ), name + ": " + band.getTaxon() + "'s name lies within its box's width" );
                    }
                    ok &= ck( ( refused > 0 ) && ( refused < tp.cladeBandsForTest().size() ), name
                            + ": precondition: some boxes are named and some too narrow; " + refused + " refused" );
                }
            }
            else {
                ok &= ck( new HashSet<String>( names ).equals( GENERA ) && tp.cladeBoxNamesLeftOutForTest().isEmpty(),
                          name + ": every box named, none left out; got " + names + " / "
                                  + tp.cladeBoxNamesLeftOutForTest() );
                final Rectangle canvas = new Rectangle( 0, 0, tp.getWidth(), tp.getHeight() );
                for ( int i = 0; i < shapes.size(); ++i ) {
                    ok &= ck( canvas.contains( shapes.get( i ).getBounds2D() ), name + ": " + names.get( i )
                            + " stays on the canvas (" + shapes.get( i ).getBounds2D() + " in " + canvas + ")" );
                }
                int labels = 0;
                for ( final double[] b : tp.labelBoxesForTest() ) {
                    ++labels;
                    final Area label = new Area( orientedBox( b ) );
                    for ( int i = 0; i < shapes.size(); ++i ) {
                        final Area a = new Area( shapes.get( i ) );
                        a.intersect( label );
                        ok &= ck( a.isEmpty(), name + ": " + names.get( i ) + " sits on no tip label" );
                    }
                }
                ok &= ck( labels >= 8, name + ": precondition: tip labels drawn to compare with (" + labels + ")" );
                if ( type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
                    ok &= onTheirSpokes( tp, name );
                }
            }
            tp.setRecordLabelBoxesForTest( false );
        }
        return ok;
    }

    /** Each unrooted name starts on its clade's spoke (the wedge's middle), past the sector -- or, where the fan's
     *  reserve falls short (a name sits at the clade root's offset plus the sector's radius, which can pass the
     *  longest tip label the fan keeps room for), it was pulled back onto the canvas and touches its edge. The
     *  reserve keeps that rare: measured 1 name of 5 on this demo (Gammacoronavirus, 11 px), 2026-10-31. */
    private static boolean onTheirSpokes( final TreePanel tp, final String name ) {
        boolean ok = true;
        final List<String> names = tp.cladeBoxNamesDrawnForTest();
        int pulled = 0;
        for ( int i = 0; i < names.size(); ++i ) {
            final CladeBand b = band( tp.cladeBandsForTest(), names.get( i ) );
            final double[] w = tp.unrootedWedgeForTest( b.getRoot() );
            final double mid = ( w[ 0 ] + w[ 1 ] ) / 2.0;
            final Rectangle2D box = tp.cladeBoxNameShapesForTest().get( i ).getBounds2D();
            final double ax = ( ( Math.cos( mid ) >= 0 ) ? box.getMinX() : box.getMaxX() ) - b.getRoot().getXcoord();
            final double ay = box.getCenterY() - b.getRoot().getYcoord();
            final double along = ( ax * Math.cos( mid ) ) + ( ay * Math.sin( mid ) );
            final double across = Math.abs( ( ax * Math.sin( mid ) ) - ( ay * Math.cos( mid ) ) );
            final boolean on_spoke = ( across <= 1 ) && ( along >= tp.unrootedCladeReachForTest( b.getRoot() ) );
            final boolean at_edge = ( box.getMinX() <= 0.5 ) || ( box.getMinY() <= 0.5 )
                    || ( box.getMaxX() >= ( tp.getWidth() - 0.5 ) ) || ( box.getMaxY() >= ( tp.getHeight() - 0.5 ) );
            pulled += on_spoke ? 0 : 1;
            ok &= ck( on_spoke || at_edge, name + ": " + names.get( i ) + " starts on its spoke past its sector, or "
                    + "was pulled back to the canvas edge (" + across + " px off the spoke, " + along + " px out, the "
                    + "sector " + tp.unrootedCladeReachForTest( b.getRoot() ) + ", box " + box + ")" );
        }
        ok &= ck( pulled <= 1, name + ": the fan keeps room for the names, so at most one is pulled back (measured 1); "
                + pulled + " were" );
        return ok;
    }

    private static boolean unrootedSectors( final MainFrame mf, final TreePanel tp ) {
        boolean ok = true;
        tp.setCladeBoxLabels( false );
        tp.setCladeBands( "genus", CLADE_VIS.BOXES );
        tp.setColorByPropertyRef( null );
        setLayout( mf, tp, LAYOUTS[ 4 ] );
        tp.setRecordLabelBoxesForTest( true );
        paint( tp );
        final Phylogeny phy = tp.getPhylogeny();
        final double[] whole = tp.unrootedWedgeForTest( phy.getRoot() );
        ok &= ck( ( whole != null ) && ( Math.abs( ( whole[ 1 ] - whole[ 0 ] ) - ( 2 * Math.PI ) ) < EPS ),
                  "the root's wedge is the full turn (" + ( ( whole[ 1 ] - whole[ 0 ] ) - ( 2 * Math.PI ) ) + ")" );
        for ( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isExternal() ) {
                continue;
            }
            final double[] w = tp.unrootedWedgeForTest( n );
            double at = w[ 0 ];
            for ( int i = 0; i < n.getNumberOfDescendants(); ++i ) {
                final double[] c = tp.unrootedWedgeForTest( n.getChildNode( i ) );
                ok &= ck( Math.abs( c[ 0 ] - at ) < EPS, "a child's wedge starts where its sibling's ends ("
                        + ( c[ 0 ] - at ) + ")" );
                at = c[ 1 ];
            }
            ok &= ck( Math.abs( at - w[ 1 ] ) < EPS, "the children's wedges fill their parent's (" + ( at - w[ 1 ] )
                    + ")" );
        }
        final List<CladeBand> bands = tp.cladeBandsForTest();
        for ( int i = 0; i < bands.size(); ++i ) {
            for ( int j = i + 1; j < bands.size(); ++j ) {
                final double[] a = tp.unrootedWedgeForTest( bands.get( i ).getRoot() );
                final double[] b = tp.unrootedWedgeForTest( bands.get( j ).getRoot() );
                ok &= ck( ( Math.min( a[ 1 ], b[ 1 ] ) - Math.max( a[ 0 ], b[ 0 ] ) ) <= EPS, bands.get( i ).getTaxon()
                        + " and " + bands.get( j ).getTaxon() + ": distinct clades' wedges do not overlap ("
                        + ( Math.min( a[ 1 ], b[ 1 ] ) - Math.max( a[ 0 ], b[ 0 ] ) ) + ")" );
            }
        }
        // every drawn tip label of a clade inside its sector: all four corners within the radius, the axis in the angle
        int checked = 0;
        for ( final double[] box : tp.labelBoxesForTest() ) {
            final CladeBand band = bandOf( bands, phy.getNode( ( long ) box[ 5 ] ) );
            if ( band == null ) {
                continue;
            }
            ++checked;
            final PhylogenyNode r = band.getRoot();
            final double[] w = tp.unrootedWedgeForTest( r );
            final double reach = tp.unrootedCladeReachForTest( r );
            final double c = Math.cos( box[ 4 ] ), s = Math.sin( box[ 4 ] );
            for ( int i = -1; i <= 1; i += 2 ) {
                for ( int j = -1; j <= 1; j += 2 ) {
                    final double d = Math.hypot( ( box[ 0 ] + ( i * box[ 2 ] * c ) ) - ( j * box[ 3 ] * s ) - r.getXcoord(),
                                                 ( box[ 1 ] + ( i * box[ 2 ] * s ) + ( j * box[ 3 ] * c ) ) - r.getYcoord() );
                    ok &= ck( d <= ( reach + 0.01 ), band.getTaxon() + ": a label corner lies within the sector's "
                            + "radius (" + d + " vs " + reach + ")" );
                }
                final double x = box[ 0 ] + ( i * box[ 2 ] * c ) - r.getXcoord();
                final double y = box[ 1 ] + ( i * box[ 2 ] * s ) - r.getYcoord();
                if ( Math.hypot( x, y ) > 1 ) {
                    double a = Math.atan2( y, x );
                    while ( a < w[ 0 ] ) {
                        a += 2 * Math.PI;
                    }
                    ok &= ck( a <= ( w[ 1 ] + EPS ), band.getTaxon() + ": a label's axis lies within the sector's "
                            + "angle" );
                }
            }
        }
        ok &= ck( checked >= 8, "precondition: drawn tip labels inside clades to check (" + checked + ")" );
        // the deltacoronaviruses sit ON their clade's root, and the sector still reaches past their labels
        final CladeBand delta = band( bands, "Deltacoronavirus" );
        double on_root = 0;
        for ( final PhylogenyNode tip : delta.getRoot().getAllExternalDescendants() ) {
            on_root = Math.max( on_root, Math.hypot( tip.getXcoord() - delta.getRoot().getXcoord(), tip.getYcoord()
                    - delta.getRoot().getYcoord() ) );
        }
        ok &= ck( on_root < 0.5, "precondition: the deltacoronaviruses sit on their clade's root (" + on_root + ")" );
        ok &= ck( tp.unrootedCladeReachForTest( delta.getRoot() ) > 40, "the zero-length clade's sector reaches past "
                + "its labels (" + tp.unrootedCladeReachForTest( delta.getRoot() ) + " px)" );
        tp.setRecordLabelBoxesForTest( false );
        // painted where it is claimed: just inside the radius near the wedge's edge, the box's tint; just past, none
        final BufferedImage img = render( tp );
        for ( final String taxon : new String[] { "Betacoronavirus", "Deltacoronavirus" } ) {
            final CladeBand b = band( bands, taxon );
            final double[] w = tp.unrootedWedgeForTest( b.getRoot() );
            final double a = w[ 0 ] + ( 0.1 * ( w[ 1 ] - w[ 0 ] ) );
            final double reach = tp.unrootedCladeReachForTest( b.getRoot() );
            final Color in = pixel( img, tp, b.getRoot(), a, reach - 2 );
            final Color out = pixel( img, tp, b.getRoot(), a, reach + 12 );
            final Color k = b.getColor();
            final Color expected = new Color( blend( k.getRed() ), blend( k.getGreen() ), blend( k.getBlue() ) );
            ok &= ck( near( in, expected, 8 ), taxon + ": the sector is painted inside its radius (" + in + " vs "
                    + expected + ")" );
            ok &= ck( near( out, Color.WHITE, 3 ) || near( out, tp.getBackground(), 3 ), taxon
                    + ": nothing painted just past the radius (" + out + ")" );
        }
        return ok;
    }

    private static boolean unrootedBarsAndBrackets( final MainFrame mf, final TreePanel tp ) {
        boolean ok = true;
        tp.setColorByPropertyRef( "vipr:Host" );
        setLayout( mf, tp, LAYOUTS[ 4 ] );
        final Set<String> multi = new HashSet<String>( GENERA );
        multi.remove( "Alphaletovirus" );
        final Color text = tp.getTreeColorSet().getSequenceColor();
        // the fan keeps room past the tip labels for what sits at the mouths: named boxes their names, bars and
        // brackets the 9 px square and its 4 px gap as well; boxes named in the legend nothing
        tp.setCladeBoxLabels( false );
        tp.setCladeBands( "genus", CLADE_VIS.BOXES, false );
        final int legend_reserve = tp.unrootedCladeNameReserveForTest();
        tp.setCladeBoxLabels( true );
        tp.setCladeBands( "genus", CLADE_VIS.BOXES, false );
        final int named_reserve = tp.unrootedCladeNameReserveForTest();
        tp.setCladeBands( "genus", CLADE_VIS.BARS, false );
        final int bar_reserve = tp.unrootedCladeNameReserveForTest();
        tp.setCladeBands( "genus", CLADE_VIS.BRACKETS, false );
        final int bracket_reserve = tp.unrootedCladeNameReserveForTest();
        tp.setCladeBoxLabels( false );
        ok &= ck( ( legend_reserve == 0 ) && ( named_reserve > 0 ) && ( bar_reserve == ( named_reserve + 13 ) )
                && ( bracket_reserve == bar_reserve ), "unrooted reserve: boxes in the legend 0, named boxes > 0, bars "
                        + "and brackets 13 px more; got " + legend_reserve + " / " + named_reserve + " / " + bar_reserve
                        + " / " + bracket_reserve );
        // a bigger tree font, the same marks: the room follows the names' new width (the widths are kept per level,
        // and a font change replaces no level)
        final TreeFontSet tfs = mf.getMainPanel().getTreeFontSet();
        final java.awt.Font base = tfs.getBaseFont();
        tfs.setBaseFont( base.deriveFont( base.getSize2D() + 8f ) );
        final int bigger = tp.unrootedCladeNameReserveForTest();
        tfs.setBaseFont( base );
        ok &= ck( ( bigger > bracket_reserve ) && ( tp.unrootedCladeNameReserveForTest() == bracket_reserve ),
                  "unrooted reserve follows the tree font: " + bracket_reserve + " -> " + bigger + " -> "
                          + tp.unrootedCladeNameReserveForTest() );
        for ( final CLADE_VIS mode : new CLADE_VIS[] { CLADE_VIS.BARS, CLADE_VIS.BRACKETS } ) {
            for ( final boolean skip : new boolean[] { true, false } ) {
                final String name = "unrooted " + mode + ( skip ? ", single-member clades skipped" : "" );
                tp.setCladeBands( "genus", mode, skip );
                mf.showWhole();
                tp.setRecordLabelBoxesForTest( true );
                paint( tp );
                final List<String> names = tp.cladeBoxNamesDrawnForTest();
                final List<Shape> shapes = tp.cladeBoxNameShapesForTest();
                ok &= ck( new HashSet<String>( names ).equals( skip ? multi : GENERA )
                        && ( names.size() == ( skip ? 4 : 5 ) ) && tp.cladeBoxNamesLeftOutForTest().isEmpty(), name
                                + ": every " + ( skip ? "multi-tip " : "" ) + "genus named once, none left out; got "
                                + names + " / " + tp.cladeBoxNamesLeftOutForTest() );
                ok &= ck( disjoint( shapes ), name + ": no two names meet" );
                final Rectangle canvas = new Rectangle( 0, 0, tp.getWidth(), tp.getHeight() );
                for ( int i = 0; i < shapes.size(); ++i ) {
                    ok &= ck( canvas.contains( shapes.get( i ).getBounds2D() ), name + ": " + names.get( i )
                            + " stays on the canvas" );
                    for ( final double[] b : tp.labelBoxesForTest() ) {
                        final Area a = new Area( shapes.get( i ) );
                        a.intersect( new Area( orientedBox( b ) ) );
                        ok &= ck( a.isEmpty(), name + ": " + names.get( i ) + " sits on no tip label" );
                    }
                }
                ok &= ck( tp.labelBoxesForTest().size() >= 8, name + ": precondition: tip labels drawn to compare with" );
                ok &= onTheirSpokes( tp, name );
                ok &= ck( ( mode == CLADE_VIS.BARS ) == ( tp.getRankLegendBoundsForTest() != null ), name
                        + ": bars keep their legend, brackets have none" );
                tp.setRecordLabelBoxesForTest( false );
                final BufferedImage img = render( tp );
                // no sector: where the box would be tinted (unrootedSectors), the canvas
                final CladeBand beta = band( tp.cladeBandsForTest(), "Betacoronavirus" );
                final double[] w = tp.unrootedWedgeForTest( beta.getRoot() );
                final Color wash = pixel( img, tp, beta.getRoot(), w[ 0 ] + ( 0.1 * ( w[ 1 ] - w[ 0 ] ) ),
                                          tp.unrootedCladeReachForTest( beta.getRoot() ) - 2 );
                ok &= ck( near( wash, Color.WHITE, 3 ) || near( wash, tp.getBackground(), 3 ), name
                        + ": no sector is painted (" + wash + ")" );
                // the square: at the name's end nearer the tree -- where no legend is drawn over it (the radial
                // layouts' legends float over the canvas: here the boxes' legend covers Alphacoronavirus's)
                int squares = 0;
                for ( int i = 0; i < names.size(); ++i ) {
                    if ( underALegend( tp, shapes.get( i ).getBounds2D() ) ) {
                        continue;
                    }
                    ++squares;
                    final CladeBand b = band( tp.cladeBandsForTest(), names.get( i ) );
                    final double[] wb = tp.unrootedWedgeForTest( b.getRoot() );
                    final boolean right = Math.cos( ( wb[ 0 ] + wb[ 1 ] ) / 2 ) >= 0;
                    final Rectangle2D box = shapes.get( i ).getBounds2D();
                    final double sx = right ? box.getMinX() : ( box.getMaxX() - 9 ); // the square is 9 px
                    final double sy = box.getCenterY() - 4.5;
                    final Color centre = new Color( img.getRGB( ( int ) Math.floor( sx + 4.5 ),
                                                                ( int ) Math.floor( sy + 4.5 ) ) );
                    // the name starts past the square: the gap between them (4 px) is clear, top to bottom
                    final int gx = ( int ) Math.floor( right ? ( sx + 9 + 2 ) : ( sx - 2 ) );
                    boolean clear = true;
                    for ( int y = ( int ) Math.ceil( box.getMinY() ); y < Math.floor( box.getMaxY() ); ++y ) {
                        final Color c = new Color( img.getRGB( gx, y ) );
                        clear &= near( c, Color.WHITE, 40 ) || near( c, tp.getBackground(), 40 );
                    }
                    ok &= ck( clear, name + ": " + names.get( i ) + " starts past its square (the gap is clear)" );
                    if ( mode == CLADE_VIS.BARS ) {
                        ok &= ck( near( centre, b.getColor(), 2 ), name + ": " + names.get( i ) + "'s square, on the "
                                + "tree's side, is its clade's colour (" + centre + " vs " + b.getColor() + ")" );
                    }
                    else {
                        ok &= ck( near( centre, Color.WHITE, 3 ) || near( centre, tp.getBackground(), 3 ), name + ": "
                                + names.get( i ) + "'s square, on the tree's side, is hollow (" + centre + ")" );
                        int darkest = 255 * 3;
                        Color edge = null;
                        for ( int dx = 0; dx <= 1; ++dx ) {
                            final Color e = new Color( img.getRGB( ( int ) Math.floor( sx ) + dx,
                                                                   ( int ) Math.floor( sy + 4.5 ) ) );
                            if ( ( e.getRed() + e.getGreen() + e.getBlue() ) < darkest ) {
                                darkest = e.getRed() + e.getGreen() + e.getBlue();
                                edge = e;
                            }
                        }
                        ok &= ck( near( edge, text, 110 ), name + ": " + names.get( i ) + "'s square is outlined in the "
                                + "text colour (" + edge + " vs " + text + ")" );
                    }
                }
                ok &= ck( squares >= 3, name + ": precondition: squares clear of the legends to look at (" + squares
                        + ")" );
            }
        }
        return ok;
    }

    /** Unrooted bars at TWO levels: every level drawn, a subfamily's name past the names of the genera inside it (the
     *  finest rank nearest the tips, as the other layouts nest their columns), the fan's room the sum of the levels'. */
    private static boolean unrootedLevels( final MainFrame mf, final TreePanel tp ) {
        boolean ok = true;
        tp.setColorByPropertyRef( "vipr:Host" );
        setLayout( mf, tp, LAYOUTS[ 4 ] );
        for ( final PhylogenyNode tip : tp.getPhylogeny().getExternalNodes() ) {
            final String g = tip.getNodeData().getProperties().getProperties( "vipr:Genus" ).get( 0 ).getValue();
            if ( tip.getNodeData().getProperties().getProperties( "test:subfamily" ).isEmpty() ) {
                tip.getNodeData().getProperties().addProperty( new org.forester.phylogeny.data.Property(
                        "test:subfamily", g.equals( "Alphaletovirus" ) ? "Letovirinae" : "Orthocoronavirinae", "",
                        "xsd:string", org.forester.phylogeny.data.Property.AppliesTo.NODE ) );
            }
        }
        tp.setCladeLevels( Arrays.asList( new CladeLevel.Spec( "genus", TreePanel.CLADE_LABEL_ANGLE.VERTICAL ) ),
                           CLADE_VIS.BARS, false, TaxonSource.TREE );
        final int one = tp.unrootedCladeNameReserveForTest();
        tp.setCladeLevels( Arrays.asList( new CladeLevel.Spec( "genus", TreePanel.CLADE_LABEL_ANGLE.VERTICAL ),
                                          new CladeLevel.Spec( "subfamily", TreePanel.CLADE_LABEL_ANGLE.VERTICAL ) ),
                           CLADE_VIS.BARS, false, TaxonSource.TREE );
        final int two = tp.unrootedCladeNameReserveForTest();
        final int sub_w = mf.getMainPanel().getTreeFontSet().getFontMetricsLarge().stringWidth( "Orthocoronavirinae" );
        ok &= ck( ( tp.cladeBandLevelsForTest().size() == 2 ) && ( two >= ( one + sub_w + 13 ) ), "unrooted, two "
                + "levels: the fan keeps room for both (" + one + " for one, " + two + " for two, the wider subfamily "
                + "name " + sub_w + " px)" );
        final java.awt.Dimension frame_size = mf.getSize();
        mf.setSize( 1600, 1250 ); // room for two levels: a smaller window holds the outer names at the canvas edge
        mf.validate();
        mf.showWhole();
        tp.setRecordLabelBoxesForTest( true );
        paint( tp );
        final List<String> names = tp.cladeBoxNamesDrawnForTest();
        final List<Shape> shapes = tp.cladeBoxNameShapesForTest();
        final Set<String> all = new HashSet<String>( GENERA );
        all.add( "Letovirinae" );
        all.add( "Orthocoronavirinae" );
        ok &= ck( new HashSet<String>( names ).equals( all ) && ( names.size() == 7 )
                && tp.cladeBoxNamesLeftOutForTest().isEmpty(), "unrooted, two levels: every genus and subfamily named "
                        + "once; got " + names + " / " + tp.cladeBoxNamesLeftOutForTest() );
        ok &= ck( disjoint( shapes ), "unrooted, two levels: no two names meet" );
        final Rectangle canvas = new Rectangle( 0, 0, tp.getWidth(), tp.getHeight() );
        for ( int i = 0; i < shapes.size(); ++i ) {
            ok &= ck( canvas.contains( shapes.get( i ).getBounds2D() ), "unrooted, two levels: " + names.get( i )
                    + " stays on the canvas" );
            for ( final double[] b : tp.labelBoxesForTest() ) {
                final Area a = new Area( shapes.get( i ) );
                a.intersect( new Area( orientedBox( b ) ) );
                ok &= ck( a.isEmpty(), "unrooted, two levels: " + names.get( i ) + " sits on no tip label" );
            }
        }
        // a subfamily's mark starts a whole genus level further out than its own labels reach -- where the genus marks
        // of its clades end, as the next ring does in circular -- its square's edge nearer the tree on its spoke
        for ( final CladeBand sub : tp.cladeBandLevelsForTest().get( 1 ) ) {
            final PhylogenyNode s_root = sub.getRoot();
            final double[] w = tp.unrootedWedgeForTest( s_root );
            final boolean right = Math.cos( ( w[ 0 ] + w[ 1 ] ) / 2 ) >= 0;
            final Rectangle2D s_box = shapes.get( names.indexOf( sub.getTaxon() ) ).getBounds2D();
            final double at = Math.hypot( ( right ? s_box.getMinX() : s_box.getMaxX() ) - s_root.getXcoord(),
                                          s_box.getCenterY() - s_root.getYcoord() );
            final double own = tp.unrootedCladeReachForTest( s_root );
            ok &= ck( at >= ( ( own + one ) - 1 ), "unrooted, two levels: " + sub.getTaxon() + "'s mark starts a genus "
                    + "level past its labels (" + at + " vs its labels' reach " + own + " + the genus level " + one
                    + "), in " + tp.getWidth() + "x" + tp.getHeight() );
        }
        mf.setSize( frame_size );
        mf.validate();
        mf.showWhole();
        tp.setRecordLabelBoxesForTest( false );
        tp.setCladeBands( "genus", CLADE_VIS.BOXES );
        return ok;
    }

    /** Unrooted views of the same marks: a visible-region export keeps the names where the screen has them (the
     *  canvas they keep to is the panel's, not the crop's); scrolled, every sector reaches as far as unscrolled (a
     *  label scrolled out of view still counts); a collapsed clade's wedge and label are inside its own sector. */
    private static boolean unrootedViews( final MainFrame mf, final TreePanel tp ) {
        boolean ok = true;
        tp.setColorByPropertyRef( "vipr:Host" );
        tp.setCladeBoxLabels( true );
        tp.setCladeBands( "genus", CLADE_VIS.BOXES );
        setLayout( mf, tp, LAYOUTS[ 4 ] );
        paint( tp );
        final java.util.Map<String, Double> reach = new java.util.HashMap<String, Double>();
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            reach.put( b.getTaxon(), tp.unrootedCladeReachForTest( b.getRoot() ) );
        }
        // zoomed in, so the panel is larger than the view, and scrolled
        for ( int i = 0; i < 3; ++i ) {
            tp.getControlPanel().zoomInX( 1.25f, 1.2f );
        }
        mf.getMainPanel().validate();
        final javax.swing.JViewport port = mf.getMainPanel().getCurrentScrollPane().getViewport();
        final java.awt.Dimension view = port.getExtentSize();
        ok &= ck( ( tp.getWidth() > ( view.width + 200 ) ) && ( tp.getHeight() > ( view.height + 200 ) ),
                  "precondition: zoomed in, the panel is larger than the view (" + tp.getWidth() + "x" + tp.getHeight()
                          + " vs " + view + ")" );
        port.setViewPosition( new java.awt.Point( ( tp.getWidth() - view.width ) / 2 + 150, 60 ) );
        paint( tp );
        boolean some_off = false;
        final Rectangle vis = tp.getVisibleRect();
        for ( final PhylogenyNode tip : tp.getPhylogeny().getExternalNodes() ) {
            some_off |= !vis.contains( tip.getXcoord(), tip.getYcoord() );
        }
        ok &= ck( some_off, "precondition: scrolled, some tips (and their labels) are out of view" );
        // zoom scales every reach alike: compare each clade's reach with the others', not with the unzoomed numbers
        final java.util.Map<String, Double> scrolled = new java.util.HashMap<String, Double>();
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            scrolled.put( b.getTaxon(), tp.unrootedCladeReachForTest( b.getRoot() ) );
        }
        port.setViewPosition( new java.awt.Point( 0, 0 ) );
        paint( tp );
        port.setViewPosition( new java.awt.Point( tp.getWidth() - view.width, tp.getHeight() - view.height ) );
        paint( tp );
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            ok &= ck( Math.abs( tp.unrootedCladeReachForTest( b.getRoot() ) - scrolled.get( b.getTaxon() ) ) < 0.5,
                      b.getTaxon() + "'s sector reaches as far whichever way the view is scrolled ("
                              + scrolled.get( b.getTaxon() ) + " vs " + tp.unrootedCladeReachForTest( b.getRoot() ) + ")" );
        }
        // what is out of view is laid out for the marks but not drawn (nor claimed, nor counted)
        tp.setRecordLabelBoxesForTest( true );
        paint( tp );
        final Rectangle seen = tp.getVisibleRect();
        int drawn_out_of_view = 0;
        for ( final double[] b : tp.labelBoxesForTest() ) {
            final Rectangle2D r = orientedBox( b ).getBounds2D();
            if ( !r.intersects( seen.x - 20, seen.y - 20, seen.width + 40, seen.height + 40 ) ) {
                ++drawn_out_of_view;
            }
        }
        tp.setRecordLabelBoxesForTest( false );
        ok &= ck( drawn_out_of_view == 0, "scrolled, no label out of view is drawn (" + drawn_out_of_view + ")" );
        // a visible-region export: the names where the screen has them
        port.setViewPosition( new java.awt.Point( ( tp.getWidth() - view.width ) / 2 + 150, 60 ) );
        paint( tp );
        final List<Rectangle2D> on_screen = new java.util.ArrayList<Rectangle2D>();
        for ( final Shape sh : tp.cladeBoxNameShapesForTest() ) {
            on_screen.add( sh.getBounds2D() );
        }
        boolean outside_crop_frame = false;
        for ( final Rectangle2D r : on_screen ) {
            outside_crop_frame |= !new Rectangle( 0, 0, view.width, view.height ).contains( r );
        }
        ok &= ck( !on_screen.isEmpty() && outside_crop_frame, "precondition: names lie outside the crop's own frame "
                + "(0, 0, " + view.width + ", " + view.height + ")" );
        final boolean was = mf.getOptions().isGraphicsExportVisibleOnly();
        mf.getOptions().setGraphicsExportVisibleOnly( true );
        try {
            AptxUtil.renderPhylogenyToImage( tp.getVisibleRect().width, tp.getVisibleRect().height, tp, mf.getOptions(),
                                             false, 1, true );
        }
        finally {
            mf.getOptions().setGraphicsExportVisibleOnly( was );
        }
        final List<Shape> exported = tp.cladeBoxNameShapesForTest();
        boolean same = exported.size() == on_screen.size();
        for ( int i = 0; same && ( i < exported.size() ); ++i ) {
            final Rectangle2D a = exported.get( i ).getBounds2D(), b = on_screen.get( i );
            same = ( Math.abs( a.getX() - b.getX() ) < 0.5 ) && ( Math.abs( a.getY() - b.getY() ) < 0.5 );
        }
        ok &= ck( same, "a visible-region export keeps the names where the screen has them; screen " + on_screen
                + ", export " + exported );
        mf.showWhole();
        // a collapsed clade (collapsed in the rectangular layout: the unrooted one refuses to collapse) is inside the
        // sector of the band rooted at it -- its wedge and its label
        setLayout( mf, tp, LAYOUTS[ 0 ] );
        final PhylogenyNode beta = band( tp.cladeBandsForTest(), "Betacoronavirus" ).getRoot();
        tp.collapse( beta );
        setLayout( mf, tp, LAYOUTS[ 4 ] );
        tp.setRecordLabelBoxesForTest( true );
        paint( tp );
        final double[] wedge = tp.collapsedWedgeRadialForTest( beta );
        double label_far = 0;
        for ( final double[] b : tp.labelBoxesForTest() ) {
            if ( ( long ) b[ 5 ] == beta.getId() ) {
                label_far = Math.max( label_far, farCorner( b, beta.getXcoord(), beta.getYcoord() ) );
            }
        }
        tp.setRecordLabelBoxesForTest( false );
        final double r = tp.unrootedCladeReachForTest( beta );
        ok &= ck( beta.isCollapse() && ( label_far > wedge[ 1 ] ), "precondition: Betacoronavirus collapsed, its label "
                + "drawn past its wedge (" + label_far + ", wedge " + wedge[ 1 ] + ")" );
        ok &= ck( r >= Math.max( wedge[ 1 ], label_far ), "a collapsed clade's sector reaches past its wedge and its "
                + "label (" + r + " vs wedge " + wedge[ 1 ] + ", label " + label_far + ")" );
        setLayout( mf, tp, LAYOUTS[ 0 ] );
        tp.collapse( beta ); // open again
        setLayout( mf, tp, LAYOUTS[ 4 ] );
        return ok;
    }

    /** How far from ({@code x},{@code y}) the farthest corner of the oriented box {@code b} lies (its own corners,
     *  not its bounds': a slanted box's bounds reach farther than it does). */
    private static double farCorner( final double[] b, final double x, final double y ) {
        double far = 0;
        final java.awt.geom.PathIterator it = orientedBox( b ).getPathIterator( null );
        final double[] c = new double[ 6 ];
        while ( !it.isDone() ) {
            if ( it.currentSegment( c ) != java.awt.geom.PathIterator.SEG_CLOSE ) {
                far = Math.max( far, Math.hypot( c[ 0 ] - x, c[ 1 ] - y ) );
            }
            it.next();
        }
        return far;
    }

    /** Unrooted bars over tips with domain architectures (drawn past the labels, along the spokes), no name on a
     *  NEIGHBOUR's domain track: a crowded generated tree, one-tip genera between large ones, every tip with a track,
     *  the fan turned in steps (the demo below covers a clade's own tracks). */
    private static boolean crowdedDomains() throws Exception {
        final int[] sizes = { 9, 1, 7, 1, 1, 8, 1, 6, 1, 1 };
        final StringBuilder x = new StringBuilder( "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<phyloxml "
                + "xmlns=\"http://www.phyloxml.org\"><phylogeny rooted=\"true\"><clade>" );
        int n = 0;
        for ( int g = 0; g < sizes.length; ++g ) {
            x.append( "<clade><branch_length>" ).append( 0.05 + ( 0.07 * ( g % 3 ) ) ).append( "</branch_length>" );
            for ( int i = 0; i < sizes[ g ]; ++i ) {
                ++n;
                x.append( "<clade><name>Isolate_" ).append( g ).append( '_' ).append( i ).append( "</name>" )
                        .append( "<branch_length>" ).append( 0.1 + ( 0.03 * ( i % 4 ) ) ).append( "</branch_length>" )
                        .append( "<sequence><domain_architecture length=\"300\"><domain from=\"20\" to=\"140\" "
                                + "confidence=\"0.001\">Kinase</domain><domain from=\"180\" to=\"280\" "
                                + "confidence=\"0.01\">SH2</domain></domain_architecture></sequence>" )
                        .append( "<property ref=\"test:genus\" datatype=\"xsd:string\" applies_to=\"node\">Genus"
                                + ( char ) ( 'A' + g ) + "virus</property></clade>" );
            }
            x.append( "</clade>" );
        }
        x.append( "</clade></phylogeny></phyloxml>" );
        final File f = File.createTempFile( "crowded-domains", ".xml" );
        f.deleteOnExit();
        java.nio.file.Files.write( f.toPath(), x.toString().getBytes( java.nio.charset.StandardCharsets.UTF_8 ) );
        final Phylogeny phy = FigureRenderer.readTrees( f )[ 0 ];
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "crowded domains" ) );
        final boolean[] ok = { true };
        final int tips = n;
        try {
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getCurrentTreePanel();
                setLayout( mf[ 0 ], tp, LAYOUTS[ 4 ] );
                tp.setCladeLevels( Arrays.asList( new CladeLevel.Spec( "genus", TreePanel.CLADE_LABEL_ANGLE.VERTICAL ) ),
                                   CLADE_VIS.BARS, false, TaxonSource.TREE );
                int moved = 0, checked = 0;
                final double start = tp.getStartingAngle();
                for ( double a = 0; a < ( 2 * Math.PI ); a += Math.PI / 6 ) {
                    tp.setStartingAngle( a );
                    mf[ 0 ].showWhole();
                    paint( tp );
                    // (a tip Auto-hide drops draws no track either)
                    ok[ 0 ] &= ck( tp.radialDomainBoxesForTest().size() >= ( tips / 2 ), "precondition: domain tracks "
                            + "drawn (" + tp.radialDomainBoxesForTest().size() + " of " + tips + ")" );
                    moved += tp.cladeBoxNamesMovedForTest() + tp.cladeBoxNamesLeftOutForTest().size();
                    final List<Shape> shapes = tp.cladeBoxNameShapesForTest();
                    for ( int i = 0; i < shapes.size(); ++i ) {
                        ++checked;
                        for ( final double[] d : tp.radialDomainBoxesForTest().values() ) {
                            final Area o = new Area( shapes.get( i ) );
                            o.intersect( new Area( orientedBox( d ) ) );
                            ok[ 0 ] &= ck( o.isEmpty(), "crowded, turned " + Math.round( Math.toDegrees( a ) )
                                    + " degrees: " + tp.cladeBoxNamesDrawnForTest().get( i ) + " sits on no domain track" );
                        }
                    }
                }
                tp.setStartingAngle( start );
                ok[ 0 ] &= ck( ( checked >= 60 ) && ( moved > 0 ), "precondition: names drawn (" + checked + ") and some "
                        + "had to move or give way (" + moved + ")" );
                // circular: the names ride past the domain ring (it paints its domains after the clade marks, which
                // ring them by the measure that counts the track)
                setLayout( mf[ 0 ], tp, LAYOUTS[ 3 ] );
                tp.setCladeBoxLabels( true );
                tp.setCladeLevels( Arrays.asList( new CladeLevel.Spec( "genus", TreePanel.CLADE_LABEL_ANGLE.VERTICAL ) ),
                                   CLADE_VIS.BOXES, false, TaxonSource.TREE );
                mf[ 0 ].showWhole();
                paint( tp );
                final List<Shape> circ = tp.cladeBoxNameShapesForTest();
                ok[ 0 ] &= ck( ( circ.size() >= 5 ) && ( tp.radialDomainBoxesForTest().size() >= ( tips / 2 ) ),
                               "circular: precondition: names (" + circ.size() + ") and domain tracks ("
                                       + tp.radialDomainBoxesForTest().size() + ") drawn" );
                for ( int i = 0; i < circ.size(); ++i ) {
                    for ( final double[] d : tp.radialDomainBoxesForTest().values() ) {
                        final Area o = new Area( circ.get( i ) );
                        o.intersect( new Area( orientedBox( d ) ) );
                        ok[ 0 ] &= ck( o.isEmpty(), "circular: " + tp.cladeBoxNamesDrawnForTest().get( i )
                                + " sits on no domain track" );
                    }
                }
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( ( JFrame ) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    /** A LARGE collapsed clade (60 tips, short names) in unrooted: its wedge is tall enough that its corners reach
     *  past its own label, and its band's mark reaches past them. */
    private static boolean bigCollapsed() throws Exception {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode big = new PhylogenyNode();
        final PhylogenyNode small = new PhylogenyNode();
        root.addAsChild( big );
        root.addAsChild( small );
        big.setDistanceToParent( 1 );
        small.setDistanceToParent( 1 );
        for ( int i = 0; i < 63; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( ( i < 60 ? "b" : "s" ) + i );
            tip.setDistanceToParent( 0.1 );
            tip.getNodeData().setProperties( new org.forester.phylogeny.data.PropertiesList() );
            tip.getNodeData().getProperties().addProperty( new org.forester.phylogeny.data.Property( "test:genus",
                    ( i < 60 ) ? "Big" : "Small", "", "xsd:string", org.forester.phylogeny.data.Property.AppliesTo.NODE ) );
            ( ( i < 60 ) ? big : small ).addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "big collapsed" ) );
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getCurrentTreePanel();
                tp.setColorByPropertyRef( null );
                setLayout( mf[ 0 ], tp, LAYOUTS[ 0 ] );
                tp.setCladeLevels( Arrays.asList( new CladeLevel.Spec( "genus", TreePanel.CLADE_LABEL_ANGLE.VERTICAL ) ),
                                   CLADE_VIS.BARS, false, TaxonSource.TREE );
                tp.collapse( big );
                setLayout( mf[ 0 ], tp, LAYOUTS[ 4 ] );
                tp.setRecordLabelBoxesForTest( true );
                paint( tp );
                final double[] wedge = tp.collapsedWedgeRadialForTest( big );
                final double corner = Math.hypot( wedge[ 1 ], wedge[ 2 ] );
                double label_far = 0;
                for ( final double[] b : tp.labelBoxesForTest() ) {
                    if ( ( long ) b[ 5 ] == big.getId() ) {
                        label_far = Math.max( label_far, farCorner( b, big.getXcoord(), big.getYcoord() ) );
                    }
                }
                tp.setRecordLabelBoxesForTest( false );
                ok[ 0 ] &= ck( big.isCollapse() && ( corner > ( label_far + 2 ) ), "precondition: the collapsed clade's "
                        + "wedge corners (" + corner + ") reach past its label (" + label_far + ")" );
                ok[ 0 ] &= ck( tp.unrootedCladeReachForTest( big ) >= ( corner - 0.5 ), "a large collapsed clade's mark "
                        + "reaches past its wedge's corners (" + tp.unrootedCladeReachForTest( big ) + " vs " + corner
                        + ")" );
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( ( JFrame ) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    /** Unrooted bars on the domain demo: each clade's mark reaches past its own tips' domain tracks, and no name sits
     *  on one. */
    private static boolean groupedDomains( final MainFrame mf ) {
        final TreePanel tp = mf.getCurrentTreePanel();
        boolean ok = true;
        final PhylogenyNode root = tp.getPhylogeny().getRoot();
        for ( int c = 0; c < root.getNumberOfDescendants(); ++c ) {
            for ( final PhylogenyNode tip : root.getChildNode( c ).getAllExternalDescendants() ) {
                if ( tip.getNodeData().getProperties() == null ) {
                    tip.getNodeData().setProperties( new org.forester.phylogeny.data.PropertiesList() );
                }
                tip.getNodeData().getProperties().addProperty( new org.forester.phylogeny.data.Property( "test:genus",
                        "Group" + ( char ) ( 'A' + c ), "", "xsd:string",
                        org.forester.phylogeny.data.Property.AppliesTo.NODE ) );
            }
        }
        setLayout( mf, tp, LAYOUTS[ 4 ] );
        tp.setCladeLevels( Arrays.asList( new CladeLevel.Spec( "genus", TreePanel.CLADE_LABEL_ANGLE.VERTICAL ) ),
                           CLADE_VIS.BARS, false, TaxonSource.TREE );
        mf.showWhole();
        paint( tp );
        final java.util.Map<Long, double[]> domains = tp.radialDomainBoxesForTest();
        ok &= ck( ( domains.size() >= 4 ) && ( tp.cladeBandsForTest().size() >= 2 ), "precondition: domains drawn ("
                + domains.size() + ") under at least two clades (" + tp.cladeBandsForTest().size() + ")" );
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            double far = 0;
            for ( final PhylogenyNode tip : b.getRoot().getAllExternalDescendants() ) {
                final double[] d = domains.get( tip.getId() );
                if ( d != null ) {
                    far = Math.max( far, farCorner( d, b.getRoot().getXcoord(), b.getRoot().getYcoord() ) );
                }
            }
            ok &= ck( tp.unrootedCladeReachForTest( b.getRoot() ) >= ( far - 1 ), b.getTaxon() + "'s mark reaches past "
                    + "its tips' domains (" + tp.unrootedCladeReachForTest( b.getRoot() ) + " vs " + far + ")" );
        }
        final List<String> names = tp.cladeBoxNamesDrawnForTest();
        final List<Shape> shapes = tp.cladeBoxNameShapesForTest();
        ok &= ck( !names.isEmpty(), "precondition: names drawn" );
        for ( int i = 0; i < shapes.size(); ++i ) {
            for ( final double[] d : domains.values() ) {
                final Area a = new Area( shapes.get( i ) );
                a.intersect( new Area( orientedBox( d ) ) );
                ok &= ck( a.isEmpty(), names.get( i ) + " sits on no domain track" );
            }
        }
        return ok;
    }

    private static boolean underALegend( final TreePanel tp, final Rectangle2D r ) {
        for ( final Rectangle legend : new Rectangle[] { tp.getPropertyLegendBounds(), tp.getRankLegendBoundsForTest(),
                tp.getSizeLegendBounds() } ) {
            if ( ( legend != null ) && legend.intersects( r ) ) {
                return true;
            }
        }
        return false;
    }

    private static boolean figureRoundTrip( final TreePanel tp ) {
        boolean ok = true;
        tp.setCladeBoxLabels( true );
        tp.setCladeBands( "genus", CLADE_VIS.BOXES );
        final String saved = FigureSpec.capture( tp ).toPropertyValue();
        ok &= ck( saved.contains( ";cladenames=labels" ), "a figure saves the names choice: " + saved );
        tp.setCladeBoxLabels( false );
        tp.setCladeBands( "genus", CLADE_VIS.BOXES );
        FigureSpec.parse( saved ).applyTo( tp );
        ok &= ck( tp.isCladeBoxLabels() && !tp.hasRankLegend(), "a saved figure brings the names back (and no legend)" );
        // a figure from before the choice existed: the boxes' legend, as it always was
        FigureSpec.parse( saved.replace( ";cladenames=labels", "" ) ).applyTo( tp );
        ok &= ck( !tp.isCladeBoxLabels() && tp.hasRankLegend(), "an older figure opens with the boxes' legend" );
        return ok;
    }

    // ---- the crowded twin: the name the rectangular layout refuses ----------------------------------------------

    private static boolean onCrowded( final MainFrame mf ) {
        final TreePanel tp = mf.getCurrentTreePanel();
        boolean ok = true;
        paint( tp );
        ok &= ck( ( tp.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) && tp.isCladeBoxLabels()
                && ( tp.cladeBandsForTest().size() == 13 ), "the twin opens rectangular, 13 genus boxes named" );
        // the precondition its README row quotes: in the window as it opens, a one-tip box is shorter than a line
        ok &= ck( ( 2 * tp.getYdistance() ) < tp.getLargeFontHeight(), "precondition: a one-tip box ("
                + ( 2 * tp.getYdistance() ) + " px) is shorter than a line (" + tp.getLargeFontHeight() + " px) in a "
                + ( ( JFrame ) mf ).getHeight() + " px window" );
        ok &= ck( new HashSet<String>( tp.cladeBoxNamesDrawnForTest() ).equals( CROWDED_TALL ),
                  "rectangular: the seven large genera are named, the six one-tip genera are not; got "
                          + tp.cladeBoxNamesDrawnForTest() );
        ok &= ck( disjoint( tp.cladeBoxNameShapesForTest() ), "rectangular: no two names meet" );
        for ( final Object[] layout : new Object[][] { LAYOUTS[ 3 ], LAYOUTS[ 4 ] } ) {
            setLayout( mf, tp, layout );
            tp.setRecordLabelBoxesForTest( true );
            paint( tp );
            final List<double[]> labels = new java.util.ArrayList<double[]>( tp.labelBoxesForTest() );
            tp.setRecordLabelBoxesForTest( false );
            ok &= ck( labels.size() >= 15, layout[ 0 ] + ": precondition: tip labels drawn to compare with ("
                    + labels.size() + ")" );
            final List<Shape> label_shapes = new java.util.ArrayList<Shape>();
            for ( final double[] b : labels ) {
                label_shapes.add( orientedBox( b ) );
            }
            // no name on a tip label: the circular names ride past every label as drawn (the label measure is the
            // SHORTENED name, the radial label the full one), the unrooted ones keep clear of them
            for ( final Shape name : drawn( tp ) ) {
                ok &= ck( !touchesAny( name, label_shapes ), layout[ 0 ] + ": a name sits on no tip label ("
                        + name.getBounds2D() + ")" );
            }
            final List<String> drawn = tp.cladeBoxNamesDrawnForTest();
            final List<String> left_out = tp.cladeBoxNamesLeftOutForTest();
            final Set<String> all = new HashSet<String>( drawn );
            all.addAll( left_out );
            ok &= ck( new HashSet<String>( drawn ).containsAll( CROWDED_TALL ), layout[ 0 ]
                    + ": the large genera are always named (placed first); got " + drawn );
            ok &= ck( all.equals( CROWDED_ALL ) && ( ( drawn.size() + left_out.size() ) == CROWDED_ALL.size() ),
                      layout[ 0 ] + ": every box named or left out, once; got " + drawn + " / " + left_out );
            ok &= ck( disjoint( drawn( tp ) ), layout[ 0 ] + ": no two names meet" );
            int previous = Integer.MAX_VALUE;
            for ( final String taxon : drawn ) { // placed largest clade first, so a small clade never takes a large one's place
                final int size = band( tp.cladeBandsForTest(), taxon ).getRoot().getNumberOfExternalNodes();
                ok &= ck( size <= previous, layout[ 0 ] + ": names placed largest clade first; " + taxon + " (" + size
                        + " tips) after a clade of " + previous );
                previous = size;
            }
            for ( int i = 0; i < left_out.size(); ++i ) {
                final Shape out = tp.cladeBoxLeftOutShapesForTest().get( i );
                ok &= ck( !disjoint( with( drawn( tp ), out ) ) || ( ( layout == LAYOUTS[ 4 ] )
                        && touchesAny( out, label_shapes ) ), layout[ 0 ] + ": " + left_out.get( i )
                                + " was left out only because it would touch a drawn name"
                                + ( ( layout == LAYOUTS[ 4 ] ) ? " or a label" : "" ) );
            }
            if ( layout == LAYOUTS[ 3 ] ) {
                ok &= ck( !left_out.isEmpty(), "circular: precondition: neighbouring one-tip genera leave a name out" );
            }
            else {
                ok &= ck( tp.cladeBoxNamesMovedForTest() > 0, "unrooted: precondition: some names had to move "
                        + "outward (" + tp.cladeBoxNamesMovedForTest() + ")" );
            }
        }
        // circular with LARGE node marks: each label starts past its mark, so the labels reach past the label measure
        // the ring is laid out by -- the names still ride past every label as drawn
        {
            final short size = mf.getOptions().getDefaultNodeShapeSize();
            mf.getOptions().setDefaultNodeShapeSize( ( short ) 40 );
            try {
                setLayout( mf, tp, LAYOUTS[ 3 ] );
                tp.setRecordLabelBoxesForTest( true );
                paint( tp );
                final List<Shape> labels = new java.util.ArrayList<Shape>();
                for ( final double[] b : tp.labelBoxesForTest() ) {
                    labels.add( orientedBox( b ) );
                }
                tp.setRecordLabelBoxesForTest( false );
                ok &= ck( !labels.isEmpty() && !drawn( tp ).isEmpty(), "circular, large node marks: precondition: "
                        + "labels and names drawn" );
                for ( final Shape name : drawn( tp ) ) {
                    ok &= ck( !touchesAny( name, labels ), "circular, large node marks: a name sits on no tip label ("
                            + name.getBounds2D() + ")" );
                }
            }
            finally {
                mf.getOptions().setDefaultNodeShapeSize( size );
            }
        }
        setLayout( mf, tp, LAYOUTS[ 0 ] );
        return ok;
    }

    // ---- the dialog --------------------------------------------------------------------------------------------

    private static boolean dialog() throws Exception {
        final Phylogeny phy = FigureRenderer.readTrees( demoFile( DEMO ) )[ 0 ];
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "clade boxes dialog" ) );
        final boolean[] ok = { true };
        try {
            // the demo opens with its boxes named beside each box: named in the legend, beside again, then removed --
            // each a change of the room the marks take past the labels, which the dialog must give back or take
            // ...and "boxbracket": the boxes plus a bracket column at their OWN rank (the only other rank this demo
            // places is none), which names each box -- the names choice greyed, the boxes' colours the key
            for ( final String choice : new String[] { "legend", "beside", "boxbracket", "bars", "remove" } ) {
                final boolean box_bracket = "boxbracket".equals( choice );
                // "boxbracket" picks "named beside each box" too: the choice is kept for next time, and the names are
                // still not drawn, because the brackets at the boxes' rank name them
                final boolean beside = "beside".equals( choice ) || box_bracket;
                final boolean bars_run = "bars".equals( choice );
                final boolean remove = "remove".equals( choice );
                final boolean[] seen = { false };
                SwingUtilities.invokeAndWait( () -> {
                    final Timer t = new Timer( 150, null );
                    final int[] ticks = { 0 };
                    t.addActionListener( e -> {
                        for ( final Window w : Window.getWindows() ) {
                            if ( !( w instanceof JDialog ) || !w.isShowing() ) {
                                continue;
                            }
                            final String title = ( ( JDialog ) w ).getTitle();
                            final JOptionPane op = find( ( ( JDialog ) w ).getContentPane(), JOptionPane.class );
                            if ( "Annotate Clades by Rank".equals( title ) && ( op != null ) ) {
                                seen[ 0 ] = true;
                                final Component form = ( Component ) op.getMessage();
                                final JRadioButton boxes = radio( form, "Shaded boxes" );
                                final JRadioButton bars = radio( form, "Bars" );
                                final JRadioButton legend = radio( form, "named in the legend" );
                                final JRadioButton named = radio( form, "named beside each box" );
                                bars.doClick();
                                ok[ 0 ] &= ck( !legend.isEnabled() && !named.isEnabled(), "the names choice is greyed "
                                        + "for bars" );
                                boxes.doClick();
                                ok[ 0 ] &= ck( legend.isEnabled() && named.isEnabled(), "the names choice is offered "
                                        + "for shaded boxes" );
                                ( beside ? named : legend ).doClick();
                                if ( bars_run ) { // bars, with "named in the legend" left selected (and greyed)
                                    legend.doClick();
                                    bars.doClick();
                                }
                                if ( box_bracket ) {
                                    final List<javax.swing.JComboBox<?>> combos = combos( form );
                                    @SuppressWarnings( "unchecked" )
                                    final javax.swing.JComboBox<String> ranks2 = ( javax.swing.JComboBox<String> ) combos.get( 2 );
                                    final JRadioButton as_brackets = radio( form, "...drawn as brackets" );
                                    ok[ 0 ] &= ck( !as_brackets.isEnabled(), "the extra levels' style is greyed until "
                                            + "an extra level is set" );
                                    for ( int i = 0; i < ranks2.getItemCount(); ++i ) {
                                        if ( "genus".equals( RankSourcePanel.bareRank( ranks2.getItemAt( i ) ) ) ) {
                                            ranks2.setSelectedIndex( i );
                                        }
                                    }
                                    ok[ 0 ] &= ck( "genus".equals( RankSourcePanel.bareRank( ( String ) ranks2.getSelectedItem() ) ),
                                                   "precondition: the second level set to genus, the boxes' own rank" );
                                    ok[ 0 ] &= ck( as_brackets.isEnabled(), "with an extra level past shaded boxes its "
                                            + "style is offered" );
                                    as_brackets.doClick();
                                    ok[ 0 ] &= ck( !legend.isEnabled() && !named.isEnabled(), "a column at the boxes' "
                                            + "own rank names each box: the names choice is greyed" );
                                }
                                if ( remove ) { // the rank chooser's first entry: stop drawing the clade marks
                                    @SuppressWarnings( "unchecked" )
                                    final javax.swing.JComboBox<String> ranks = find( form, javax.swing.JComboBox.class );
                                    ranks.setSelectedIndex( 0 );
                                    ok[ 0 ] &= ck( MainFrame.CLADE_LEVEL_CLEAR.equals( ranks.getSelectedItem() ),
                                                   "precondition: the first entry removes the marks" );
                                    t.stop(); // no result dialog follows
                                }
                                op.setValue( JOptionPane.OK_OPTION );
                                w.dispose();
                                return;
                            }
                            if ( ( title != null ) && title.startsWith( "Annotate Clades by Rank (" ) ) {
                                t.stop(); // the result: close it, and the flow returns
                                if ( op != null ) {
                                    op.setValue( JOptionPane.CLOSED_OPTION );
                                }
                                w.dispose();
                                return;
                            }
                        }
                        if ( ++ticks[ 0 ] > 60 ) { // ~9 s: close whatever modal is up, fail, never hang
                            t.stop();
                            for ( final Window w : Window.getWindows() ) {
                                if ( ( w instanceof JDialog ) && w.isShowing() ) {
                                    w.dispose();
                                }
                            }
                        }
                    } );
                    t.start();
                    mf[ 0 ].setLastTaxonSourceForTest( TaxonSource.TREE );
                    mf[ 0 ].labelCladesByRank();
                } );
                final TreePanel tp = mf[ 0 ].getCurrentTreePanel();
                final boolean[] state = new boolean[ 5 ];
                final double[] spread = new double[ 2 ];
                SwingUtilities.invokeAndWait( () -> {
                    paint( tp );
                    state[ 0 ] = tp.isCladeBoxLabels();
                    state[ 1 ] = tp.getRankLegendBoundsForTest() != null;
                    state[ 2 ] = tp.hasCladeBands();
                    state[ 3 ] = ( tp.getCladeBandsMode() == CLADE_VIS.BRACKETS ) && "genus".equals( tp.cladeBoxRank() )
                            && tp.cladeColumnRanks().equals( java.util.Arrays.asList( "genus" ) )
                            && ( tp.drawnCladeColumnCount() > 0 );
                    state[ 4 ] = !tp.cladeBoxNamesShownForTest(); // (the names at the mouths are the brackets')
                    spread[ 0 ] = tipSpread( tp );
                    tp.getControlPanel().fitWidth(); // what the dialog should have done already
                    paint( tp );
                    spread[ 1 ] = tipSpread( tp );
                } );
                ok[ 0 ] &= ck( seen[ 0 ], "the dialog opened" );
                if ( remove ) {
                    ok[ 0 ] &= ck( !state[ 2 ], "the dialog's remove entry removes the marks" );
                }
                else if ( bars_run ) {
                    ok[ 0 ] &= ck( state[ 2 ] && state[ 0 ], "a Bars run keeps the boxes' \"named beside each box\" "
                            + "for next time (names " + state[ 0 ] + ")" );
                }
                else if ( box_bracket ) {
                    ok[ 0 ] &= ck( state[ 2 ] && state[ 3 ], "genus boxes + a genus bracket column reach the tab" );
                    ok[ 0 ] &= ck( state[ 4 ] && state[ 1 ] && state[ 0 ], "the brackets name the boxes (no names beside "
                            + "them, though \"beside\" stays chosen) and the boxes' colours are the key (a legend); names "
                            + state[ 4 ] + ", legend " + state[ 1 ] + ", chosen " + state[ 0 ] );
                }
                else {
                    ok[ 0 ] &= ck( state[ 2 ] && ( state[ 0 ] == beside ) && ( state[ 1 ] == !beside ), "the dialog's \""
                            + ( beside ? "named beside each box" : "named in the legend" ) + "\" reaches the tab; names "
                            + state[ 0 ] + ", legend " + state[ 1 ] );
                }
                ok[ 0 ] &= ck( Math.abs( spread[ 0 ] - spread[ 1 ] ) < 1, "unrooted, " + choice + ": the dialog re-fits the "
                        + "fan to the room the marks now take (spread " + spread[ 0 ] + ", re-fitted " + spread[ 1 ] + ")" );
            }
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( ( JFrame ) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    // ---- helpers -------------------------------------------------------------------------------------------------

    private interface Exercise {

        boolean run( MainFrame mf );
    }

    private static boolean withDemo( final String path, final Exercise ex ) throws Exception {
        final Phylogeny phy = FigureRenderer.readTrees( demoFile( path ) )[ 0 ];
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "clade boxes" ) );
        final boolean[] ok = { false };
        try {
            SwingUtilities.invokeAndWait( () -> ok[ 0 ] = ex.run( mf[ 0 ] ) );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( ( JFrame ) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    private static File demoFile( final String path ) {
        return new File( System.getProperty( "user.dir" ), path );
    }

    private static void setLayout( final MainFrame mf, final TreePanel tp, final Object[] layout ) {
        tp.setPhylogenyGraphicsType( ( PHYLOGENY_GRAPHICS_TYPE ) layout[ 1 ] );
        tp.setTreeOrientation( ( TREE_ORIENTATION ) layout[ 2 ] );
        mf.showWhole();
    }

    private static void paint( final TreePanel tp ) {
        tp.paintImmediately( 0, 0, tp.getWidth(), tp.getHeight() );
        tp.paintImmediately( 0, 0, tp.getWidth(), tp.getHeight() );
    }

    private static BufferedImage render( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( tp.getWidth(), tp.getHeight(), BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            tp.printAll( g );
        }
        finally {
            g.dispose();
        }
        return img;
    }

    private static Color pixel( final BufferedImage img, final TreePanel tp, final PhylogenyNode root, final double a,
                                final double r ) {
        final int x = ( int ) Math.round( root.getXcoord() + ( r * Math.cos( a ) ) );
        final int y = ( int ) Math.round( root.getYcoord() + ( r * Math.sin( a ) ) );
        return new Color( img.getRGB( x, y ) );
    }

    private static int blend( final int c ) { // the box's colour at its alpha (46 of 255) over white
        return ( int ) Math.round( 255 - ( ( 255 - c ) * ( 46 / 255.0 ) ) );
    }

    private static boolean near( final Color a, final Color b, final int tol ) {
        return ( Math.abs( a.getRed() - b.getRed() ) <= tol ) && ( Math.abs( a.getGreen() - b.getGreen() ) <= tol )
                && ( Math.abs( a.getBlue() - b.getBlue() ) <= tol );
    }

    private static Shape orientedBox( final double[] b ) {
        final java.awt.geom.Path2D.Double p = new java.awt.geom.Path2D.Double();
        final double c = Math.cos( b[ 4 ] ), s = Math.sin( b[ 4 ] );
        final double[][] k = { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } };
        for ( int i = 0; i < 4; ++i ) {
            final double x = ( b[ 0 ] + ( k[ i ][ 0 ] * b[ 2 ] * c ) ) - ( k[ i ][ 1 ] * b[ 3 ] * s );
            final double y = b[ 1 ] + ( k[ i ][ 0 ] * b[ 2 ] * s ) + ( k[ i ][ 1 ] * b[ 3 ] * c );
            if ( i == 0 ) {
                p.moveTo( x, y );
            }
            else {
                p.lineTo( x, y );
            }
        }
        p.closePath();
        return p;
    }

    private static List<Shape> drawn( final TreePanel tp ) {
        return tp.cladeBoxNameShapesForTest();
    }

    /** {@code shapes} and {@code extra}: !disjoint(with(drawn, s)) asks whether s touches a drawn name. */
    private static List<Shape> with( final List<Shape> shapes, final Shape extra ) {
        final List<Shape> out = new java.util.ArrayList<Shape>();
        out.add( extra ); // first: disjoint() then compares it with every drawn name (drawn names are apart)
        out.addAll( shapes );
        return out;
    }

    private static boolean disjoint( final List<Shape> shapes ) {
        for ( int i = 0; i < shapes.size(); ++i ) {
            for ( int j = i + 1; j < shapes.size(); ++j ) {
                final Area a = new Area( shapes.get( i ) );
                a.intersect( new Area( shapes.get( j ) ) );
                if ( !a.isEmpty() && ( a.getBounds2D().getWidth() > 0.5 ) && ( a.getBounds2D().getHeight() > 0.5 ) ) {
                    return false;
                }
            }
        }
        return true;
    }

    private static CladeBand bandOf( final List<CladeBand> bands, final PhylogenyNode n ) {
        for ( final CladeBand b : bands ) {
            for ( PhylogenyNode x = n; x != null; x = x.getParent() ) {
                if ( x == b.getRoot() ) {
                    return b;
                }
            }
        }
        return null;
    }

    private static CladeBand band( final List<CladeBand> bands, final String taxon ) {
        for ( final CladeBand b : bands ) {
            if ( taxon.equals( b.getTaxon() ) ) {
                return b;
            }
        }
        throw new IllegalStateException( "no " + taxon + " clade" );
    }

    /** The value under the first value row of a legend, scanning down its left column. */
    private static String valueIn( final TreePanel tp, final Rectangle legend ) {
        for ( int y = legend.y + 1; y < ( legend.y + legend.height ); ++y ) {
            final String v = tp.legendValueAt( at( tp, legend.x + 12, y ) );
            if ( v != null ) {
                return v;
            }
        }
        return null;
    }

    /** Clicks the first value row of the rank legend at {@code legend} through the panel's legend-click route, picks
     *  {@code chosen} in the colour dialog that opens, and says whether that row's taxon took it on the rank legend
     *  (and nothing else did on the Color-by legend). The override is cleared again. */
    private static boolean pickOnRow( final TreePanel tp, final Rectangle legend, final Color chosen ) {
        int row_y = -1;
        String taxon = null;
        for ( int y = legend.y + 1; ( y < ( legend.y + legend.height ) ) && ( taxon == null ); ++y ) {
            taxon = tp.legendValueAt( at( tp, legend.x + 12, y ) );
            row_y = y;
        }
        if ( taxon == null ) {
            return false;
        }
        final java.util.Map<String, Color> color_by_before = new java.util.HashMap<String, Color>(
                tp.getPropertyColorScheme().getValueColors() );
        final String[] title = { null };
        final long give_up = System.currentTimeMillis() + 9000;
        final Timer t = new Timer( 100, null );
        t.addActionListener( ev -> {
            for ( final Window w : Window.getWindows() ) {
                if ( ( w instanceof JDialog ) && w.isShowing() ) {
                    final JDialog d = ( JDialog ) w;
                    final javax.swing.JColorChooser cc = find( d, javax.swing.JColorChooser.class, null );
                    final javax.swing.JButton okb = find( d, javax.swing.JButton.class, "OK" );
                    if ( ( cc != null ) && ( okb != null ) ) {
                        title[ 0 ] = d.getTitle();
                        cc.setColor( chosen );
                        okb.doClick();
                        t.stop();
                        return;
                    }
                }
            }
            if ( System.currentTimeMillis() > give_up ) { // a missing OK fails the check rather than hang the suite
                for ( final Window w : Window.getWindows() ) {
                    if ( ( w instanceof JDialog ) && w.isShowing() ) {
                        w.dispose();
                    }
                }
                t.stop();
            }
        } );
        t.start();
        tp.handleLegendClick( at( tp, legend.x + 12, row_y ) );
        t.stop();
        final boolean ok = ( title[ 0 ] != null ) && title[ 0 ].contains( taxon ) && chosen.equals( tp.rankLegendColor(
                taxon ) ) && tp.getPropertyColorScheme().getValueColors().equals( color_by_before );
        tp.clearRankColorOverride( "genus", taxon );
        return ok;
    }

    private static <T> T find( final Component c, final Class<T> type, final String text ) {
        if ( type.isInstance( c ) && ( ( text == null ) || ( ( c instanceof javax.swing.AbstractButton )
                && text.equals( ( ( javax.swing.AbstractButton ) c ).getText() ) ) ) ) {
            return type.cast( c );
        }
        if ( c instanceof java.awt.Container ) {
            for ( final Component x : ( ( java.awt.Container ) c ).getComponents() ) {
                final T found = find( x, type, text );
                if ( found != null ) {
                    return found;
                }
            }
        }
        return null;
    }

    private static void drag( final TreePanel tp, final Rectangle legend, final int dx, final int dy ) {
        tp.startLegendDrag( at( tp, legend.x + 10, legend.y + 5 ) );
        tp.dragLegend( at( tp, legend.x + 10 + dx, legend.y + 5 + dy ) );
        tp.endLegendDrag();
    }

    private static MouseEvent at( final TreePanel tp, final int x, final int y ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_PRESSED, 0, 0, x, y, 1, false );
    }

    /** Every combo box under {@code c}, in component order. */
    private static List<javax.swing.JComboBox<?>> combos( final Component c ) {
        final List<javax.swing.JComboBox<?>> out = new java.util.ArrayList<>();
        if ( c instanceof javax.swing.JComboBox ) {
            out.add( ( javax.swing.JComboBox<?> ) c );
        }
        if ( c instanceof java.awt.Container ) {
            for ( final Component x : ( ( java.awt.Container ) c ).getComponents() ) {
                out.addAll( combos( x ) );
            }
        }
        return out;
    }

    private static JRadioButton radio( final Component c, final String prefix ) {
        if ( ( c instanceof JRadioButton ) && ( ( JRadioButton ) c ).getText().trim().startsWith( prefix ) ) {
            return ( JRadioButton ) c;
        }
        if ( c instanceof java.awt.Container ) {
            for ( final Component x : ( ( java.awt.Container ) c ).getComponents() ) {
                final JRadioButton r = radio( x, prefix );
                if ( r != null ) {
                    return r;
                }
            }
        }
        return null;
    }

    /** Whether {@code s} meets any of {@code others} (by more than half a pixel each way, as {@link #disjoint}). */
    private static boolean touchesAny( final Shape s, final List<Shape> others ) {
        for ( final Shape o : others ) {
            final Area a = new Area( s );
            a.intersect( new Area( o ) );
            if ( !a.isEmpty() && ( a.getBounds2D().getWidth() > 0.5 ) && ( a.getBounds2D().getHeight() > 0.5 ) ) {
                return true;
            }
        }
        return false;
    }

    /** The shorter side of a turned rectangle (a radial name's box). */
    private static double thickness( final Shape s ) {
        final java.awt.geom.PathIterator it = s.getPathIterator( null );
        final double[] c = new double[ 6 ];
        final List<double[]> pts = new java.util.ArrayList<double[]>();
        while ( !it.isDone() && ( pts.size() < 3 ) ) {
            if ( it.currentSegment( c ) != java.awt.geom.PathIterator.SEG_CLOSE ) {
                pts.add( new double[] { c[ 0 ], c[ 1 ] } );
            }
            it.next();
        }
        final double a = Math.hypot( pts.get( 1 )[ 0 ] - pts.get( 0 )[ 0 ], pts.get( 1 )[ 1 ] - pts.get( 0 )[ 1 ] );
        final double b = Math.hypot( pts.get( 2 )[ 0 ] - pts.get( 1 )[ 0 ], pts.get( 2 )[ 1 ] - pts.get( 1 )[ 1 ] );
        return Math.min( a, b );
    }

    /** How far the farthest tip of the last paint lies from the root (the unrooted fan's spread). */
    private static double tipSpread( final TreePanel tp ) {
        final PhylogenyNode root = tp.getPhylogeny().getRoot();
        double far = 0;
        for ( final PhylogenyNode tip : tp.getPhylogeny().getExternalNodes() ) {
            far = Math.max( far, Math.hypot( tip.getXcoord() - root.getXcoord(), tip.getYcoord() - root.getYcoord() ) );
        }
        return far;
    }

    private static <T> T find( final Component c, final Class<T> k ) {
        if ( k.isInstance( c ) ) {
            return k.cast( c );
        }
        if ( c instanceof java.awt.Container ) {
            for ( final Component x : ( ( java.awt.Container ) c ).getComponents() ) {
                final T r = find( x, k );
                if ( r != null ) {
                    return r;
                }
            }
        }
        return null;
    }

    private static Set<String> set( final String... s ) {
        return new HashSet<String>( Arrays.asList( s ) );
    }

    private static boolean ck( final boolean b, final String msg ) {
        if ( !b ) {
            System.out.println( "  [CladeBoxesTest] " + msg );
        }
        return b;
    }
}
