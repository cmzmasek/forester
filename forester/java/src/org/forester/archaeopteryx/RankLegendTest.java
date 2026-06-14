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

    private static MouseEvent at( final TreePanel tp, final int x, final int y ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_PRESSED, System.currentTimeMillis(), 0, x, y, 1, false );
    }

    /** Three internal "order" nodes, each over two leaves, so colorizing by rank "order" yields 3 colors. */
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
