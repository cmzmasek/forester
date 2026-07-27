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
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
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
 * Integration test for the in-legend controls on a real {@link MainFrame}/{@link TreePanel}: a categorical
 * color-by-property legend with more distinct values than the default cap must draw a clickable sort toggle
 * and a clickable "+N more" footer; clicking them flips the sort order and expands/collapses the legend.
 * Guarded to a no-op on a headless box (needs FlatLaf via {@code createInstance}).
 */
public final class LegendControlsToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LegendControlsTool: " + ( ok ? "OK." : "FAILED." ) );
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
                                                                        "legend" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                tp.setColorByPropertyRef( "data:region" ); // 35 distinct categorical values
                frame.showWhole();
                render( tp );

                // the two controls are drawn (35 values > the default 20 cap -> "+15 more")
                if ( tp.legendSortToggleBoundsForTest() == null ) {
                    fail( ok, "sort toggle was not drawn" );
                }
                if ( tp.legendMoreBoundsForTest() == null ) {
                    fail( ok, "'+N more' control was not drawn for 35 values over the cap" );
                }
                if ( !tp.isLegendSortByCount() ) {
                    fail( ok, "default sort should be by count" );
                }

                // a static export omits the interactive chips ([by count]/[show all]) but keeps the "+N more"
                // line, so it draws fewer pixels than the on-screen legend (yet still draws the legend)
                final int screen_px = legendPixels( tp, true );
                final int export_px = legendPixels( tp, false );
                if ( export_px <= 0 ) {
                    fail( ok, "export legend drew nothing" );
                }
                if ( export_px >= screen_px ) {
                    fail( ok, "export legend should drop the interactive chips (fewer pixels than on-screen): screen="
                            + screen_px + " export=" + export_px );
                }

                // clicking the sort toggle flips count <-> A-Z (re-render between clicks: the chip label, hence
                // its bounds, changes)
                clickCenter( tp, tp.legendSortToggleBoundsForTest() );
                if ( tp.isLegendSortByCount() ) {
                    fail( ok, "clicking the sort toggle should switch to A-Z" );
                }
                render( tp );
                clickCenter( tp, tp.legendSortToggleBoundsForTest() );
                if ( !tp.isLegendSortByCount() ) {
                    fail( ok, "clicking the sort toggle again should switch back to by-count" );
                }

                // clicking "+N more" expands to show all; then a "show fewer" collapses back to the default cap
                render( tp );
                clickCenter( tp, tp.legendMoreBoundsForTest() );
                if ( tp.legendMaxEntries() != Integer.MAX_VALUE ) {
                    fail( ok, "clicking '+N more' should show all entries" );
                }
                render( tp );
                if ( tp.legendMoreBoundsForTest() == null ) {
                    fail( ok, "an expanded legend must offer a 'show fewer' control" );
                }
                clickCenter( tp, tp.legendMoreBoundsForTest() );
                if ( tp.legendMaxEntries() != 20 ) {
                    fail( ok, "clicking 'show fewer' should collapse to the default cap (20), got "
                            + tp.legendMaxEntries() );
                }

                // switching to a DIFFERENT legend subject must reset the expand state (not leak "show all")
                clickCenter( tp, tp.legendMoreBoundsForTest() ); // expand the region legend again
                render( tp );
                if ( tp.legendMaxEntries() != Integer.MAX_VALUE ) {
                    fail( ok, "test setup: legend should be expanded before the subject switch" );
                }
                tp.setColorByPropertyRef( "data:kind" ); // a different property = a different legend subject
                render( tp );
                if ( tp.legendMaxEntries() != 20 ) {
                    fail( ok, "switching legend subject must reset the expand cap to 20, got " + tp.legendMaxEntries() );
                }

                // a non-null-but-EMPTY color scheme draws no legend, so its stale bounds must not stay clickable
                final Rectangle stale = tp.getPropertyLegendBounds();
                tp.setColorByPropertyRef( "no:such:field" ); // no tip has this -> empty scheme, no legend drawn
                render( tp );
                if ( stale != null ) {
                    final MouseEvent in_old = new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(),
                                                              0, stale.x + ( stale.width / 2 ),
                                                              stale.y + ( stale.height / 2 ), 1, false );
                    if ( tp.isOnPropertyLegend( in_old ) ) {
                        fail( ok, "an empty color scheme must not leave a clickable phantom legend" );
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
        System.out.println( "  [LegendControlsToolTest] " + msg );
        ok[ 0 ] = false;
    }

    private static void clickCenter( final TreePanel tp, final Rectangle r ) {
        if ( r == null ) {
            return;
        }
        final int cx = r.x + ( r.width / 2 );
        final int cy = r.y + ( r.height / 2 );
        tp.handleLegendClick( new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(), 0, cx, cy, 1,
                                              false ) );
    }

    /** Non-white pixel count of the legend alone, drawn on-screen (draggable) or in static-export mode. */
    private static int legendPixels( final TreePanel tp, final boolean draggable ) {
        final int w = 700, h = 500;
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.drawLegendForTest( g, new Rectangle( 0, 0, w, h ), draggable );
        g.dispose();
        int n = 0;
        for( int x = 0; x < w; ++x ) {
            for( int y = 0; y < h; ++y ) {
                if ( ( img.getRGB( x, y ) & 0x00FFFFFF ) != 0x00FFFFFF ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Renders offscreen so the legend (and its control bounds) are laid out for the next hit-test. */
    private static void render( final TreePanel tp ) {
        final int w = 700, h = 500;
        tp.setSize( w, h );
        tp.validate();
        tp.doLayout();
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.printAll( g );
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
            pl.addProperty( new Property( "data:kind", "k" + ( i % 4 ), "", "xsd:string", AppliesTo.NODE ) );
            leaf.getNodeData().setProperties( pl );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }
}
