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

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;

/**
 * Integration test for dragging the "Color by:" legend: hit-testing the legend box, moving it by a
 * mouse drag (the box follows, clamped into the viewport), and double-click / reset returning it to
 * the default corner. Guarded to a no-op on a headless box (or when the panel has no usable
 * viewport). Needs FlatLaf on the classpath (via {@code MainFrameApplication.createInstance}), so it
 * is run standalone, not as part of the headless suite.
 */
public final class LegendDragTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LegendDrag: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Configuration conf = new Configuration( null, false, false, true );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { hostTree() }, conf, "lg" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                tp.setColorByPropertyRef( "repseq:host" );
                mf[ 0 ].getMainPanel().getControlPanel().showWhole();
                final Rectangle vp = tp.getVisibleRect();
                if ( ( vp.width < 300 ) || ( vp.height < 300 ) ) {
                    ok[ 0 ] = true; // no usable viewport in this environment; nothing to assert
                    ( (JFrame) mf[ 0 ] ).dispose();
                    return;
                }
                paint( tp, vp );
                final Rectangle home = tp.getPropertyLegendBounds();
                if ( home == null ) {
                    ok[ 0 ] = false;
                }
                else {
                    // hit-test: a point inside the legend vs. one in the opposite corner
                    if ( !tp.isOnPropertyLegend( at( tp, home.x + ( home.width / 2 ), home.y + 3 ) )
                            || tp.isOnPropertyLegend( at( tp, vp.x + 5, vp.y + 5 ) ) ) {
                        ok[ 0 ] = false;
                    }
                    // drag it left 80 / down 120; the box must follow
                    tp.startLegendDrag( at( tp, home.x + 10, home.y + 5 ) );
                    tp.dragLegend( at( tp, home.x + 10 - 80, home.y + 5 + 120 ) );
                    tp.endLegendDrag();
                    paint( tp, vp );
                    final Rectangle moved = tp.getPropertyLegendBounds();
                    if ( ( Math.abs( moved.x - ( home.x - 80 ) ) > 1 ) || ( Math.abs( moved.y - ( home.y + 120 ) ) > 1 ) ) {
                        ok[ 0 ] = false;
                    }
                    // reset returns it to the default corner
                    tp.resetLegendPosition();
                    paint( tp, vp );
                    final Rectangle back = tp.getPropertyLegendBounds();
                    if ( ( Math.abs( back.x - home.x ) > 1 ) || ( Math.abs( back.y - home.y ) > 1 ) ) {
                        ok[ 0 ] = false;
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

    private static MouseEvent at( final TreePanel tp, final int x, final int y ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_PRESSED, System.currentTimeMillis(), 0, x, y, 1, false );
    }

    private static Phylogeny hostTree() {
        final PhylogenyNode root = new PhylogenyNode();
        final String[] hosts = { "Homo sapiens", "Mus musculus", "Gallus gallus" };
        int id = 0;
        for( final String host : hosts ) {
            for( int c = 0; c < 4; ++c ) {
                final PhylogenyNode leaf = new PhylogenyNode();
                leaf.setName( "n" + ( id++ ) );
                final PropertiesList pl = new PropertiesList();
                pl.addProperty( new Property( "repseq:host", host, "", "xsd:string", Property.AppliesTo.NODE ) );
                leaf.getNodeData().setProperties( pl );
                root.addAsChild( leaf );
            }
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private LegendDragTest() {
    }
}
