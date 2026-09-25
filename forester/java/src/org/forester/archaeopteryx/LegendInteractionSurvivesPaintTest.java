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

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;

/**
 * A legend is only a legend if you can still USE it after it has been drawn: move it, re-order it ([by count] /
 * [A-Z]), expand it ([show all] / [show fewer]) and click a colour square to recolour a value.
 * <p>
 * All of that is hit-tested against state the paint RECORDS while drawing -- the box, the chips, and the value-row
 * layout. The legend is also MEASURED, by running that same drawing code into a scratch image inside a sentinel
 * {@code (0, 0, 100000, 100000)} rectangle, and the measurement used to leave its off-screen rectangles behind. Any
 * caller that measured AFTER the real draw -- the paint-time FPS readout did, because it asks for the
 * legend-column reserve to place itself -- therefore replaced every hit region with coordinates ~99843 px to the
 * right, and every legend went dead: no drag, no chips, no recolour. This pins that measuring leaves no trace.
 */
public final class LegendInteractionSurvivesPaintTest {

    private static final int W = 1200;
    private static final int H = 800;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LegendInteractionSurvivesPaint: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { coloredTree() }, new Configuration(), "legendlive" ) );
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
            SwingUtilities.invokeAndWait( () -> {
                tp.getOptions().setShowOverview( false );
                tp.setOvOn( false );
                cp.demoSelectColorByProperty( "test:grp" ); // a categorical field -> a legend with chips
                // The trigger: the readout asks for the legend-column reserve to place itself, and it is painted
                // LAST. With it off the measurement still runs (the time-tree badge asks too) but only BEFORE the
                // legend is drawn, so the damage was invisible.
                tp.getOptions().setShowFps( true );
                paint( tp );
            } );
            final Rectangle box = tp.getPropertyLegendBounds();
            if ( box == null ) {
                fail( ok, "precondition: the tree must draw a property-color legend" );
                dispose( mf );
                return ok[ 0 ];
            }
            // (1) the legend box is ON SCREEN -- the bug left it at x ~99843 on a 1200 px panel
            if ( ( box.x < 0 ) || ( box.y < 0 ) || ( box.x > W ) || ( box.y > H ) ) {
                fail( ok, "the legend's hit region must be on screen after a paint, got " + box );
            }
            // (2) ...and a click on it is actually recognised, which is what a drag and a recolour both start from
            if ( !tp.isOnPropertyLegend( at( tp, box.x + ( box.width / 2 ), box.y + 3 ) ) ) {
                fail( ok, "a click on the drawn legend must hit it" );
            }
            if ( tp.isOnPropertyLegend( at( tp, box.x - 60, box.y + 3 ) ) ) {
                fail( ok, "a click well outside the legend must NOT hit it" );
            }
            // (3) the in-legend CONTROLS are on screen too: [by count]/[A-Z] and, when the legend is capped,
            // "+N more". These are separate fields from the box, and the bug left them behind at the sentinel.
            final Rectangle sort = tp.legendSortToggleBoundsForTest();
            if ( sort == null ) {
                fail( ok, "a categorical legend must offer the [by count] / [A-Z] control" );
            }
            else if ( !box.contains( sort.x + 1, sort.y + 1 ) ) {
                fail( ok, "the sort control must sit INSIDE the legend box; box=" + box + " chip=" + sort );
            }
            final Rectangle more = tp.legendMoreBoundsForTest();
            if ( more == null ) {
                fail( ok, "precondition: this fixture has more values than the cap, so '+N more' must be offered" );
            }
            else if ( !box.contains( more.x + 1, more.y + 1 ) ) {
                fail( ok, "the '+N more' control must sit INSIDE the legend box; box=" + box + " chip=" + more );
            }
            // (4) the controls WORK: clicking them changes what the next paint draws
            if ( more != null ) {
                final int[] rows_before = { tp.legendRowCountForTest() };
                SwingUtilities.invokeAndWait( () -> {
                    tp.handleLegendClick( at( tp, more.x + ( more.width / 2 ), more.y + ( more.height / 2 ) ) );
                    paint( tp );
                } );
                if ( tp.legendRowCountForTest() <= rows_before[ 0 ] ) {
                    fail( ok, "clicking '+N more' must show more rows, got " + rows_before[ 0 ] + " -> "
                            + tp.legendRowCountForTest() );
                }
            }
            // (5) a value ROW resolves to its value -- the click that opens the colour chooser for a swatch.
            // The point is derived from the BOX, deliberately, not from the recorded row layout: computing it
            // from _legend_rows_top/_legend_row_height would move the probe along with the very corruption this
            // is meant to detect, and the check would pass however broken the layout was.
            final Rectangle box_now = tp.getPropertyLegendBounds();
            if ( box_now.height < 60 ) {
                fail( ok, "precondition: the legend must be tall enough to have value rows, got " + box_now );
            }
            else {
                boolean any_row_resolves = false;
                for( int y = box_now.y + 20; y < ( box_now.y + box_now.height - 4 ); y += 4 ) {
                    if ( tp.legendValueAtForTest( at( tp, box_now.x + 12, y ) ) != null ) {
                        any_row_resolves = true;
                        break;
                    }
                }
                if ( !any_row_resolves ) {
                    fail( ok, "no point inside the drawn legend resolves to a value -- the recolour gesture is "
                            + "dead (row layout left at the measuring pass's coordinates?)" );
                }
            }
            // (6) and the legend still MOVES
            final Rectangle before_drag = tp.getPropertyLegendBounds();
            SwingUtilities.invokeAndWait( () -> {
                tp.startLegendDrag( at( tp, before_drag.x + ( before_drag.width / 2 ), before_drag.y + 3 ) );
                tp.dragLegend( at( tp, before_drag.x + ( before_drag.width / 2 ) + 90, before_drag.y + 70 ) );
                paint( tp );
            } );
            final Rectangle after = tp.getPropertyLegendBounds();
            if ( ( after == null ) || after.getLocation().equals( before_drag.getLocation() ) ) {
                fail( ok, "the legend must move when dragged: " + before_drag + " -> " + after );
            }
            dispose( mf );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
        return ok[ 0 ];
    }

    private static void paint( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.paintPhylogeny( g, false, false, W, H, 0, 0 );
        g.dispose();
    }

    private static MouseEvent at( final TreePanel tp, final int x, final int y ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, 0, 0, x, y, 1, false );
    }

    /** 32 distinct values, comfortably over DEFAULT_LEGEND_MAX_ENTRIES (20), so the legend is capped and
     *  really does offer "+N more" -- otherwise the control check below would pass vacuously. */
    private static Phylogeny coloredTree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 64; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "tip_" + i );
            tip.setDistanceToParent( 0.1 + ( i * 0.01 ) );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "test:grp", "group_" + ( i % 32 ), "", "xsd:string", AppliesTo.NODE ) );
            tip.getNodeData().setProperties( pl );
            root.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static void dispose( final MainFrame[] mf ) throws Exception {
        SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [LegendInteractionSurvivesPaintTest] " + message );
        ok[ 0 ] = false;
    }

    private LegendInteractionSurvivesPaintTest() {
        // not instantiable
    }
}
