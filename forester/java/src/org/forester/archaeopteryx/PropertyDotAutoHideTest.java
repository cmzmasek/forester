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
import java.awt.event.MouseEvent;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.Options.TREE_ORIENTATION;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;

/**
 * The Color-by tip dot is drawn for EVERY tip of a rectangular tree, whether or not "Auto-hide Labels" draws the
 * tip's label -- as the radial layouts always did. Until 0.11.184 the rectangular dot rode the label routine, so a
 * crowded tree, its rows closer than the font: the labels are thinned to every k-th) showed a dot only beside a
 * drawn label (Christian, 2026-10-16: "always show the colored circle even if the label is hidden"). Held in all
 * three rectangular orientations, by sampling the pixel at every tip -- on 5 px rows at the smallest node size, so
 * the labels are thinned (that starts under 7 px) while the 4 px dots stay apart (on tighter rows a later tip's dot
 * covers its neighbour's centre, and the sample would read the overlap order, not the rule); what is deliberately
 * NOT drawn: no dot with "Show External Data" off, none with Color-by off.
 */
public final class PropertyDotAutoHideTest {

    private static final int    TIPS = 90;
    private static final float  ROW  = 5f; // px between tips: the labels are thinned (under 7 px), the 4 px dots stay apart
    private static final String REF  = "demo:group";

    public static void main( final String[] args ) {
        System.out.println( test() ? "PropertyDotAutoHideTest: OK." : "PropertyDotAutoHideTest: FAILED." );
    }

    public static boolean test() {
        if ( java.awt.GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final Phylogeny phy = star();
        final MainFrame[] mf = new MainFrame[ 1 ];
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "dots under auto-hide" ) );
            SwingUtilities.invokeAndWait( () -> {
                final JFrame f = mf[ 0 ]._window;
                f.setSize( 1300, 900 ); // room for the whole panel in every orientation: a tip off the viewport is culled
                f.validate();
                final MainPanel mp = mf[ 0 ].getMainPanel();
                final TreePanel tp = mp.getCurrentTreePanel();
                mp.getOptions().setDefaultNodeShapeSize( (short) 1 ); // 4 px dots: apart on 5 px rows
                tp.setColorByPropertyRef( REF );
                layRows( tp );
                if ( tp.baseDotSize() != 4 ) {
                    ok[ 0 ] = fail( "fixture: the dot must be 4 px at node size 1, got " + tp.baseDotSize() );
                }
                // the fixture: a rectangular tree whose labels ARE thinned (the factor the paint uses is above 1
                // and the checkbox is on), Color-by active, the panel sized to paint
                if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) {
                    ok[ 0 ] = fail( "fixture: the tree must open rectangular" );
                }
                if ( !mp.getControlPanel().checkboxForTest( DisplayOption.DYNAMICALLY_HIDE_DATA ).isSelected() ) {
                    ok[ 0 ] = fail( "fixture: Auto-hide Labels must be on (the default)" );
                }
                if ( tp.dynamicHidingFactorForTest() <= 1 ) {
                    ok[ 0 ] = fail( "fixture: " + TIPS + " tips in a 600 px window must thin the labels, factor "
                            + tp.dynamicHidingFactorForTest() );
                }
                if ( ( tp.getWidth() < 200 ) || ( tp.getHeight() < 200 ) ) {
                    ok[ 0 ] = fail( "fixture: the panel must have a size to paint, got " + tp.getSize() );
                }
                if ( !ok[ 0 ] ) {
                    return;
                }
                final Color[] colors = new Color[ TIPS ];
                int i = 0;
                for( final PhylogenyNode tip : phy.getExternalNodes() ) {
                    colors[ i++ ] = tp.getPropertyBasedColor( tip );
                }
                int distinct = 0;
                for( int a = 0; a < colors.length; ++a ) {
                    if ( colors[ a ].equals( Color.BLACK ) || colors[ a ].equals( Color.WHITE ) ) {
                        ok[ 0 ] = fail( "fixture: a property colour must be neither black nor white (the ink and the paper)" );
                    }
                    if ( ( a == 0 ) || !colors[ a ].equals( colors[ a - 1 ] ) ) {
                        ++distinct;
                    }
                }
                if ( distinct < 4 ) {
                    ok[ 0 ] = fail( "fixture: the tips must cycle through several colours" );
                }
                for( final TREE_ORIENTATION o : new TREE_ORIENTATION[] { TREE_ORIENTATION.ROOT_LEFT,
                        TREE_ORIENTATION.ROOT_TOP, TREE_ORIENTATION.ROOT_BOTTOM } ) {
                    tp.setTreeOrientation( o );
                    layRows( tp );
                    final int[] c = countDots( tp, phy, colors );
                    final int missing = c[ 0 ];
                    if ( c[ 1 ] < ( ( TIPS * 9 ) / 10 ) ) {
                        ok[ 0 ] = fail( "fixture: only " + c[ 1 ] + " of " + TIPS + " tips are clear of the legend in " + o );
                    }
                    if ( tp.dynamicHidingFactorForTest() <= 1 ) {
                        ok[ 0 ] = fail( "fixture: the labels must be thinned in " + o + " too, factor "
                                + tp.dynamicHidingFactorForTest() );
                    }
                    if ( tp.getYdistance() != ROW ) {
                        ok[ 0 ] = fail( "fixture: the rows must stay " + ROW + " px apart in " + o + ", got " + tp.getYdistance() );
                    }
                    if ( missing > 0 ) {
                        ok[ 0 ] = fail( missing + " of " + TIPS + " tips have no dot in " + o + " (labels thinned to every "
                                + tp.dynamicHidingFactorForTest() + "th; rows " + tp.getYdistance() + " px, panel " + tp.getWidth()
                                + "x" + tp.getHeight() + " in a viewport " + tp.getVisibleRect().width + "x"
                                + tp.getVisibleRect().height + ", dot " + tp.baseDotSize() + " px)" );
                    }
                }
                tp.setTreeOrientation( TREE_ORIENTATION.ROOT_LEFT );
                layRows( tp );
                // what is deliberately not drawn
                mp.getControlPanel().checkboxForTest( DisplayOption.DISPLAY_EXTERNAL_DATA ).doClick();
                if ( tp.isShowExternalDataForThisTab() ) {
                    ok[ 0 ] = fail( "fixture: the click must turn Show External Data off" );
                }
                else if ( !noneDrawn( countDots( tp, phy, colors ) ) ) {
                    ok[ 0 ] = fail( "with Show External Data off no tip draws a dot" );
                }
                mp.getControlPanel().checkboxForTest( DisplayOption.DISPLAY_EXTERNAL_DATA ).doClick();
                if ( countDots( tp, phy, colors )[ 0 ] != 0 ) {
                    ok[ 0 ] = fail( "precondition: with Show External Data back on every tip has its dot again" );
                }
                tp.setColorByPropertyRef( null );
                if ( !noneDrawn( countDots( tp, phy, colors ) ) ) {
                    ok[ 0 ] = fail( "with Color-by off no tip draws a dot" );
                }
            } );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ok[ 0 ] = fail( "unexpected: " + e );
        }
        finally {
            try {
                if ( mf[ 0 ] != null ) {
                    SwingUtilities.invokeAndWait( () -> mf[ 0 ].dispose() );
                }
            }
            catch ( final Exception e ) {
                e.printStackTrace();
            }
        }
        return ok[ 0 ];
    }

    /** Lays the tips ROW px apart (a fit to the window would pack them far tighter) and sizes the panel to the tree,
     *  within the viewport, so no tip is culled as off-screen. */
    private static void layRows( final TreePanel tp ) {
        tp.setYdistance( ROW );
        tp.resetPreferredSize();
        tp.setSize( tp.getPreferredSize() );
        tp.validate();
    }

    /** Paints the panel as the screen does and counts the tips with no pixel of their property colour within one
     *  pixel of the tip (the dot is placed on truncated logical coordinates and, in a vertical orientation, rotated,
     *  so the exact centre can land a pixel off the mapped point; the rows are 5 px, so the window is inside the
     *  4 px dot and clear of both neighbours). A tip under the legend overlay is not sampled (the legend paints
     *  over it, in the top right -- where a vertical orientation's first tips run); at most a tenth may be. */
    private static boolean noneDrawn( final int[] c ) {
        return ( c[ 0 ] == c[ 1 ] ) && ( c[ 1 ] >= ( ( TIPS * 9 ) / 10 ) );
    }

    /** {tips sampled and found without a dot, tips sampled}. */
    private static int[] countDots( final TreePanel tp, final Phylogeny phy, final Color[] colors ) {
        final BufferedImage img = new BufferedImage( tp.getWidth(), tp.getHeight(), BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, img.getWidth(), img.getHeight() );
        tp.printAll( g );
        g.dispose();
        int missing = 0;
        int under_legend = 0;
        int i = 0;
        for( final PhylogenyNode tip : phy.getExternalNodes() ) {
            final Point2D.Double p = tp.screenPoint( tip.getXcoord(), tip.getYcoord() );
            final int cx = (int) Math.round( p.x ), cy = (int) Math.round( p.y );
            if ( tp.isOnAnyLegend( new MouseEvent( tp, MouseEvent.MOUSE_MOVED, 0, 0, cx, cy, 0, false ) ) ) {
                ++under_legend;
                ++i;
                continue;
            }
            boolean found = false;
            for( int y = cy - 1; ( y <= ( cy + 1 ) ) && !found; ++y ) {
                for( int x = cx - 1; ( x <= ( cx + 1 ) ) && !found; ++x ) {
                    found = ( x >= 0 ) && ( y >= 0 ) && ( x < img.getWidth() ) && ( y < img.getHeight() )
                            && close( new Color( img.getRGB( x, y ) ), colors[ i ] );
                }
            }
            if ( !found ) {
                ++missing;
            }
            ++i;
        }
        return new int[] { missing, TIPS - under_legend };
    }

    private static boolean close( final Color a, final Color b ) {
        return ( Math.abs( a.getRed() - b.getRed() ) <= 48 ) && ( Math.abs( a.getGreen() - b.getGreen() ) <= 48 )
                && ( Math.abs( a.getBlue() - b.getBlue() ) <= 48 );
    }

    /** 240 tips under one root, each carrying one of four group values. */
    private static Phylogeny star() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < TIPS; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "tip" + i );
            tip.setDistanceToParent( 1.0 + ( i % 7 ) );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( REF, "group " + ( i % 4 ), "", "xsd:string", AppliesTo.NODE ) );
            tip.getNodeData().setProperties( pl );
            root.addAsChild( tip );
        }
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String what ) {
        System.out.println( "  [PropertyDotAutoHideTest] " + what );
        return false;
    }

    private PropertyDotAutoHideTest() {
    }
}
