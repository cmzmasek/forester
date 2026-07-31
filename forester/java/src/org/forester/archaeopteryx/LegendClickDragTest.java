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
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;

/**
 * Integration test for the legend drag/click discrimination in {@link MouseListener}: a press on a
 * legend chip followed by only sub-threshold jitter must stay a pure click (the chip acts, the legend
 * does NOT move); a press followed by real movement is a drag (the legend follows and the trailing
 * click is swallowed, so it does not ALSO fire the chip). Regression guard for the misfire where a
 * tiny press+jitter both nudged the legend and toggled/recolored a chip. Guarded to a no-op when
 * headless or when the environment gives no usable viewport (needs FlatLaf via {@code createInstance}).
 */
public final class LegendClickDragTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LegendClickDrag: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree( 35 ) }, conf,
                                                                         "lgclick" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                tp.setColorByPropertyRef( "data:region" ); // 35 distinct values -> a categorical legend with a sort chip
                frame.showWhole();
                final Rectangle vp = tp.getVisibleRect();
                if ( ( vp.width < 300 ) || ( vp.height < 300 ) ) {
                    ( (JFrame) frame ).dispose(); // no usable viewport in this environment; nothing to assert
                    return;
                }
                paint( tp, vp );
                final Rectangle toggle = tp.legendSortToggleBoundsForTest();
                final Rectangle home = tp.getPropertyLegendBounds();
                if ( ( toggle == null ) || ( home == null ) ) {
                    fail( ok, "setup: the sort chip / legend bounds were not drawn (35 values over the cap)" );
                    ( (JFrame) frame ).dispose();
                    return;
                }
                final MouseListener ml = new MouseListener( tp );

                // --- Scenario A: press on the sort chip + only sub-threshold jitter = a pure click.
                // The chip must toggle (click honored) and the legend must NOT move (no nudge).
                final boolean sort_before_a = tp.isLegendSortByCount();
                final int ax = toggle.x + ( toggle.width / 2 );
                final int ay = toggle.y + ( toggle.height / 2 );
                ml.mousePressed( ev( tp, MouseEvent.MOUSE_PRESSED, InputEvent.BUTTON1_DOWN_MASK, ax, ay, 1 ) );
                ml.mouseDragged( ev( tp, MouseEvent.MOUSE_DRAGGED, InputEvent.BUTTON1_DOWN_MASK, ax + 2, ay + 2, 1 ) );
                ml.mouseReleased( ev( tp, MouseEvent.MOUSE_RELEASED, 0, ax + 2, ay + 2, 1 ) );
                ml.mouseClicked( ev( tp, MouseEvent.MOUSE_CLICKED, 0, ax, ay, 1 ) );
                if ( tp.isLegendSortByCount() == sort_before_a ) {
                    fail( ok, "a sub-threshold press+jitter on the sort chip should have toggled it (a click)" );
                }
                paint( tp, vp );
                final Rectangle after_a = tp.getPropertyLegendBounds();
                // The toggle relabels the legend, so its WIDTH changes and (being anchored to the top-right
                // corner) its x with it -- that is not a nudge. The invariants of a still-un-nudged legend are
                // its top edge (y) and its right edge (x + width); a real jitter-nudge would move both.
                if ( ( Math.abs( after_a.y - home.y ) > 1 )
                        || ( Math.abs( ( after_a.x + after_a.width ) - ( home.x + home.width ) ) > 1 ) ) {
                    fail( ok, "a click on a legend chip must not nudge the legend: home=" + home + " after=" + after_a );
                }

                // --- Scenario B: press on the sort chip + real movement = a drag.
                // The legend must follow the drag and the trailing click must be swallowed (chip NOT toggled).
                paint( tp, vp ); // re-record the (label-changed) chip bounds and the current home
                final Rectangle toggle_b = tp.legendSortToggleBoundsForTest();
                final Rectangle home_b = tp.getPropertyLegendBounds();
                final boolean sort_before_b = tp.isLegendSortByCount();
                if ( ( toggle_b == null ) || ( home_b == null ) ) {
                    fail( ok, "setup: the sort chip / legend bounds vanished before scenario B" );
                }
                else {
                    final int bx = toggle_b.x + ( toggle_b.width / 2 );
                    final int by = toggle_b.y + ( toggle_b.height / 2 );
                    ml.mousePressed( ev( tp, MouseEvent.MOUSE_PRESSED, InputEvent.BUTTON1_DOWN_MASK, bx, by, 1 ) );
                    ml.mouseDragged( ev( tp, MouseEvent.MOUSE_DRAGGED, InputEvent.BUTTON1_DOWN_MASK, bx - 40, by + 40, 1 ) );
                    ml.mouseReleased( ev( tp, MouseEvent.MOUSE_RELEASED, 0, bx - 40, by + 40, 1 ) );
                    // the legend must have followed the drag
                    paint( tp, vp ); // a live drag repaints continuously; refresh so the chip bounds track the move
                    final Rectangle moved = tp.getPropertyLegendBounds();
                    if ( ( Math.abs( moved.x - ( home_b.x - 40 ) ) > 2 )
                            || ( Math.abs( moved.y - ( home_b.y + 40 ) ) > 2 ) ) {
                        fail( ok, "the legend did not follow the drag: home=" + home_b + " moved=" + moved );
                    }
                    // the trailing click -- delivered onto the chip at its NEW (moved) location -- must be
                    // swallowed as the tail of the drag, NOT fire the toggle
                    final Rectangle moved_toggle = tp.legendSortToggleBoundsForTest();
                    if ( moved_toggle == null ) {
                        fail( ok, "the sort chip vanished after the drag" );
                    }
                    else {
                        ml.mouseClicked( ev( tp, MouseEvent.MOUSE_CLICKED, 0, moved_toggle.x + ( moved_toggle.width / 2 ),
                                             moved_toggle.y + ( moved_toggle.height / 2 ), 1 ) );
                        if ( tp.isLegendSortByCount() != sort_before_b ) {
                            fail( ok, "a drag's trailing click on the (moved) chip must NOT also toggle the sort" );
                        }
                    }
                }
                ( (JFrame) frame ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [LegendClickDragTest] " + msg );
        ok[ 0 ] = false;
    }

    private static MouseEvent ev( final TreePanel tp, final int id, final int mods, final int x, final int y,
                                  final int clicks ) {
        return new MouseEvent( tp, id, System.currentTimeMillis(), mods, x, y, clicks, false );
    }

    /** On-screen paint into a vp-sized buffer so the legend box and its chip bounds are recorded for hit-testing. */
    private static void paint( final TreePanel tp, final Rectangle vp ) {
        final BufferedImage img = new BufferedImage( Math.max( 1, vp.x + vp.width ), Math.max( 1, vp.y + vp.height ),
                                                     BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        tp.paintPhylogeny( g, false, false, 0, 0, 0, 0 );
        g.dispose();
    }

    private static Phylogeny tree( final int n ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < n; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:region", String.format( "v%02d", i ), "", "xsd:string",
                                          AppliesTo.NODE ) );
            leaf.getNodeData().setProperties( pl );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private LegendClickDragTest() {
    }
}
