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
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Taxonomy;

/**
 * The rollover that reads an annotation-column CELL: {@link AnnotationColumns#cellRows} (what it says) and
 * {@link TreePanel#annotationCellAt} (which cell the pointer is on), plus the gesture that ties them together.
 * <p>
 * The content half runs anywhere. The geometry half is headful and covers every layout that draws cells -- the
 * rectangular root-left, the vertical clustergram (where the point has to go back through the rotation) and the
 * circular one (where it is a radius and an angle, probed on each tip's OWN spoke rather than a convenient one) --
 * and pins the NON-hits too: the gaps between columns, the margins either side, above the first row and below the
 * last, and the unrooted layout, which draws no cells at all.
 * <p>
 * Every probe point is built from production geometry that already had to be right for something else -- the drawn
 * header's x, the tip row bands, the ring mid-radius, the tips' own drawn coordinates -- never from this test's own
 * copy of the layout arithmetic.
 */
public final class AnnotationCellHoverTest {

    private static final int W = 1100;
    private static final int H = 800;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "AnnotationCellHover: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        cellRowsOk( ok );
        if ( GraphicsEnvironment.isHeadless() ) {
            return ok[ 0 ];
        }
        return geometryAndGestureOk( ok ) && ok[ 0 ];
    }

    private static boolean fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [AnnotationCellHoverTest] " + msg );
        ok[ 0 ] = false;
        return false;
    }

    // ---- fixture -------------------------------------------------------------------------------------------------

    /**
     * 8 tips: a categorical host, four numeric matrix fields, two stacked-bar fields. tip_3 is missing {@code m2}
     * (the blank cell), and tip_7 has NO NAME but a taxonomy (the heading's fallback).
     */
    private static Phylogeny tree() {
        final String[] host = { "cat", "dog", "cat", "dog", "cow", "cat", "dog", "cow" };
        final int[][] m = { { 0, 1, 2, 3, 4, 4, 3, 2 }, { 4, 4, 3, 1, 0, 0, 1, 2 }, { 1, 3, 1, 3, 1, 3, 1, 3 },
                { 2, 2, 4, 4, 0, 0, 2, 2 } };
        final int[][] st = { { 5, 9, 2, 7, 4, 6, 8, 3 }, { 3, 1, 6, 2, 5, 4, 2, 7 } };
        final PhylogenyNode root = new PhylogenyNode();
        final List<PhylogenyNode> tips = new ArrayList<PhylogenyNode>();
        for( int i = 0; i < 8; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            if ( i == 7 ) {
                final Taxonomy tax = new Taxonomy(); // no name: the heading must fall back to something
                tax.setScientificName( "Bos taurus" );
                tip.getNodeData().setTaxonomy( tax );
            }
            else {
                tip.setName( "tip_" + i );
            }
            tip.setDistanceToParent( 0.1 + ( 0.02 * i ) );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:host", host[ i ], "", "xsd:string", AppliesTo.NODE ) );
            for( int k = 0; k < 4; ++k ) {
                if ( ( i == 3 ) && ( k == 1 ) ) {
                    continue; // tip_3 was never assessed for m2
                }
                pl.addProperty( new Property( "data:m" + ( k + 1 ), String.valueOf( m[ k ][ i ] ), "", "xsd:decimal",
                                              AppliesTo.NODE ) );
            }
            pl.addProperty( new Property( "data:x", String.valueOf( st[ 0 ][ i ] ), "", "xsd:decimal", AppliesTo.NODE ) );
            pl.addProperty( new Property( "data:y", String.valueOf( st[ 1 ][ i ] ), "", "xsd:decimal", AppliesTo.NODE ) );
            tip.getNodeData().setProperties( pl );
            tips.add( tip );
        }
        for( int g = 0; g < 4; ++g ) {
            final PhylogenyNode pair = new PhylogenyNode();
            pair.setDistanceToParent( 0.2 );
            pair.addAsChild( tips.get( 2 * g ) );
            pair.addAsChild( tips.get( ( 2 * g ) + 1 ) );
            root.addAsChild( pair );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static AnnotationColumns.ColumnSpec spec( final String field, final AnnotationColumns.Type t ) {
        return new AnnotationColumns.ColumnSpec( "data:" + field, t );
    }

    /** host | m1 m2 m3 | stacked(x, y) | m4 -- six DRAWN columns from seven specs. */
    private static List<AnnotationColumns.ColumnSpec> specs() {
        return new ArrayList<AnnotationColumns.ColumnSpec>( Arrays.asList(
                spec( "host", AnnotationColumns.Type.COLOR_STRIP ), spec( "m1", AnnotationColumns.Type.MATRIX ),
                spec( "m2", AnnotationColumns.Type.MATRIX ), spec( "m3", AnnotationColumns.Type.MATRIX ),
                spec( "x", AnnotationColumns.Type.STACKED_BAR ), spec( "y", AnnotationColumns.Type.STACKED_BAR ),
                spec( "m4", AnnotationColumns.Type.MATRIX ) ) );
    }

    // ---- what the card says --------------------------------------------------------------------------------------

    private static String render( final List<NodeHoverText.Row> rows ) {
        final StringBuilder sb = new StringBuilder();
        for( final NodeHoverText.Row r : rows ) {
            if ( sb.length() > 0 ) {
                sb.append( " | " );
            }
            sb.append( r.toString() ); // the Row's OWN rendering, so this test cannot drift from it
        }
        return sb.toString();
    }

    private static void eq( final boolean[] ok, final String what, final String got, final String want ) {
        if ( !want.equals( got ) ) {
            fail( ok, what + ": expected \"" + want + "\", got \"" + got + "\"" );
        }
    }

    private static void cellRowsOk( final boolean[] ok ) {
        final Phylogeny phy = tree();
        final AnnotationColumns cols = new AnnotationColumns( phy, specs() );
        final List<PhylogenyNode> tips = phy.getExternalNodes();
        if ( cols.size() != 6 ) {
            fail( ok, "the fixture should draw 6 columns, draws " + cols.size() );
            return;
        }
        // a MATRIX cell: which tip, the field and its raw value, and the SHARED scale the colour came from
        eq( ok, "matrix cell", render( cols.cellRows( tips.get( 0 ), 1 ) ), "Tip: tip_0 | M1: 0 | Scale: 0 – 4" );
        // the neighbouring case, differing by exactly one thing: the same column on a tip that HAS no value there.
        // "0" would be a lie -- a blank cell is absence of evidence, and every other order/colour rule says so too.
        eq( ok, "blank matrix cell", render( cols.cellRows( tips.get( 3 ), 2 ) ),
            "Tip: tip_3 | M2: not assessed | Scale: 0 – 4" );
        eq( ok, "assessed neighbour", render( cols.cellRows( tips.get( 2 ), 2 ) ),
            "Tip: tip_2 | M2: 3 | Scale: 0 – 4" );
        // a CATEGORICAL strip has no numeric scale to read a colour back from, so no Scale row
        eq( ok, "colour-strip cell", render( cols.cellRows( tips.get( 1 ), 0 ) ), "Tip: tip_1 | Host: dog" );
        // a MERGED column is several fields in one cell: every series, each with its own value
        eq( ok, "stacked-bar cell", render( cols.cellRows( tips.get( 1 ), 4 ) ), "Tip: tip_1 | X: 9 | Y: 1" );
        // an unnamed tip still has to say WHICH row this is
        eq( ok, "unnamed tip", render( cols.cellRows( tips.get( 7 ), 1 ) ),
            "Tip: Bos taurus | M1: 2 | Scale: 0 – 4" );
        final PhylogenyNode bare = new PhylogenyNode();
        eq( ok, "nothing to name it by", render( cols.cellRows( bare, 1 ) ),
            "Tip: (unnamed tip) | M1: not assessed | Scale: 0 – 4" );
    }

    // ---- which cell the pointer is on ----------------------------------------------------------------------------

    private static boolean geometryAndGestureOk( final boolean[] ok ) {
        final Phylogeny phy = tree();
        final MainFrame[] mf = new MainFrame[ 1 ];
        try {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "cell hover" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    ( ( JFrame ) frame ).setSize( W, H );
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    for( final String l : new String[] { "root-left", "clustergram", "circular" } ) {
                        layout( tp, l );
                        geometry( ok, tp, l );
                    }
                    layout( tp, "root-left" );
                    gesture( ok, tp );
                    unrooted( ok, tp );
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                    t.printStackTrace();
                }
                finally {
                    ( ( JFrame ) frame ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    private static void layout( final TreePanel tp, final String layout ) {
        tp.setAnnotationColumns( specs() );
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.TABLE );
        tp.setTreeOrientation( "clustergram".equals( layout ) ? Options.TREE_ORIENTATION.ROOT_TOP
                : Options.TREE_ORIENTATION.ROOT_LEFT );
        tp.setPhylogenyGraphicsType( "circular".equals( layout ) ? PHYLOGENY_GRAPHICS_TYPE.CIRCULAR
                : "unrooted".equals( layout ) ? PHYLOGENY_GRAPHICS_TYPE.UNROOTED
                        : PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        screenPaint( tp ); // lays out the rings / the rotation the hit-tests read
    }

    /** A paint through the SCREEN path, which is what sets the geometry the hit-tests read. */
    private static void screenPaint( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.paint( g );
        g.dispose();
    }

    /** The device point at the centre of the cell of drawn column {@code col} on tip {@code row}, or null. */
    private static Point cellPoint( final TreePanel tp, final int col, final int row ) {
        final Point on_column = tp.annotationColumnGrabPointForTest( col );
        if ( on_column == null ) {
            return null;
        }
        if ( tp.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
            // annotationColumnGrabPointForTest lands on the ring's mid-radius at 3 o'clock; swing it round to THIS
            // tip's own spoke, whose angle comes from where the tip itself was drawn
            final Point c = tp.circularCenterForTest();
            final PhylogenyNode t = tp.getPhylogeny().getExternalNodes().get( row );
            if ( c == null ) {
                return null;
            }
            final double r = Math.hypot( on_column.x - c.x, on_column.y - c.y );
            final double a = Math.atan2( t.getYcoord() - c.y, t.getXcoord() - c.x );
            return new Point( ( int ) Math.round( c.x + ( r * Math.cos( a ) ) ),
                              ( int ) Math.round( c.y + ( r * Math.sin( a ) ) ) );
        }
        if ( !tp.isVerticalOrientation() ) {
            final int[] band = tp.tipRowBandsForTest()[ row ];
            return new Point( on_column.x, band[ 0 ] + ( band[ 1 ] / 2 ) );
        }
        // A vertical orientation turns the picture a quarter turn, so the two axes swap in DEVICE space: the
        // column's device y is the one its upright header anchor already sits at, and the tip's device x is where
        // the tip itself was drawn. Both come from production geometry that is not the hit-test's.
        final java.awt.geom.Point2D.Double t = tp.screenPointFor( tp.getPhylogeny().getExternalNodes().get( row ) );
        return new Point( ( int ) Math.round( t.x ), on_column.y );
    }

    /** Every cell reads back as itself; the gaps, the margins and the rows past the ends read back as nothing. */
    private static void geometry( final boolean[] ok, final TreePanel tp, final String layout ) {
        final List<PhylogenyNode> tips = tp.getPhylogeny().getExternalNodes();
        int probed = 0;
        for( int c = 0; c < 6; ++c ) {
            for( int r = 0; r < tips.size(); ++r ) {
                final Point p = cellPoint( tp, c, r );
                if ( p == null ) {
                    fail( ok, layout + ": no probe point for cell (" + c + ", " + r + ")" );
                    return;
                }
                final TreePanel.AnnotationCell hit = tp.annotationCellAt( p.x, p.y );
                if ( hit == null ) {
                    fail( ok, layout + ": cell (" + c + ", " + r + ") at " + p + " reads back as nothing" );
                    return;
                }
                if ( ( hit.column() != c ) || ( hit.tip() != tips.get( r ) ) ) {
                    fail( ok, layout + ": the cell at " + p + " should be (" + c + ", " + r + "), reads back as ("
                            + hit.column() + ", " + hit.tip().getName() + ")" );
                    return;
                }
                ++probed;
            }
        }
        if ( probed != ( 6 * tips.size() ) ) {
            fail( ok, layout + ": " + probed + " cells probed, expected " + ( 6 * tips.size() ) );
        }
        if ( !"circular".equals( layout ) ) {
            // The grid's two axes are DEVICE axes, and a vertical orientation swaps them, so both "off the end"
            // probes step along a direction measured from the cells themselves rather than assuming x or y.
            final Point r0 = cellPoint( tp, 3, 0 );
            final Point rn = cellPoint( tp, 3, tips.size() - 1 );
            if ( tp.annotationCellAt( rn.x + ( rn.x - r0.x ), rn.y + ( rn.y - r0.y ) ) != null ) {
                fail( ok, layout + ": there is no cell a whole grid past the last tip's row" );
            }
            final Point c0 = cellPoint( tp, 0, 0 );
            final Point c1 = cellPoint( tp, 1, 0 );
            final Point c5 = cellPoint( tp, 5, 0 );
            final int ux = c1.x - c0.x;
            final int uy = c1.y - c0.y;
            if ( tp.annotationCellAt( c0.x - ( 6 * ux ), c0.y - ( 6 * uy ) ) != null ) {
                fail( ok, layout + ": there is no cell before the first column" );
            }
            if ( tp.annotationCellAt( c5.x + ( 6 * ux ), c5.y + ( 6 * uy ) ) != null ) {
                fail( ok, layout + ": there is no cell past the last column" );
            }
            // the GAP between two columns belongs to neither: walking from one column's centre to the next, some
            // point in between must read as nothing (the colour strip and the matrix beside it are drawn apart)
            boolean gap = false;
            for( int k = 1; k < 100; ++k ) {
                final int gx = c0.x + ( ( ( c1.x - c0.x ) * k ) / 100 );
                final int gy = c0.y + ( ( ( c1.y - c0.y ) * k ) / 100 );
                gap |= tp.annotationCellAt( gx, gy ) == null;
            }
            if ( !gap ) {
                fail( ok, layout + ": the gap between two columns must belong to neither of them" );
            }
        }
    }

    /** The whole gesture, through a real {@link MouseListener}: moving onto a cell shows that cell's card. */
    private static void gesture( final boolean[] ok, final TreePanel tp ) {
        final MouseListener ml = new MouseListener( tp );
        final List<PhylogenyNode> tips = tp.getPhylogeny().getExternalNodes();
        final Point cell = cellPoint( tp, 2, 4 );
        ml.mouseMoved( move( tp, cell ) );
        final TreePanel.AnnotationCell showing = tp.hoverCardCellForTest();
        if ( ( showing == null ) || ( showing.column() != 2 ) || ( showing.tip() != tips.get( 4 ) ) ) {
            fail( ok, "moving onto a cell must show THAT cell's rollover, got " + showing );
        }
        // a HEADER is for dragging, not for reading a value: the header hit wins and no cell card appears
        final Point header = tp.annotationColumnGrabPointForTest( 2 );
        ml.mouseMoved( move( tp, header ) );
        if ( tp.hoverCardCellForTest() != null ) {
            fail( ok, "a column header must not show a cell rollover, got " + tp.hoverCardCellForTest() );
        }
        // off the columns entirely: the card goes. A card has to be SHOWING first -- the header move above already
        // took one down, so without this the check would pass on a card that was never up (measured: it did).
        ml.mouseMoved( move( tp, cell ) );
        if ( !tp.isNodeDescPopupShowingForTest() ) {
            fail( ok, "fixture: a cell rollover must be up before the check that moving off takes it down" );
        }
        ml.mouseMoved( move( tp, new Point( 2, H - 2 ) ) );
        if ( tp.isNodeDescPopupShowingForTest() ) {
            fail( ok, "moving off the columns must take the rollover down" );
        }
        // ...and with Rollover switched off, a cell shows nothing at all
        tp.setShows( DisplayOption.NODE_DATA_POPUP, false );
        ml.mouseMoved( move( tp, cell ) );
        if ( tp.isNodeDescPopupShowingForTest() ) {
            fail( ok, "with Rollover off, a cell must show no card" );
        }
        tp.setShows( DisplayOption.NODE_DATA_POPUP, true );
        ml.mouseMoved( move( tp, cell ) );
        if ( tp.hoverCardCellForTest() == null ) {
            fail( ok, "...and switching Rollover back on must bring it back" );
        }
    }

    private static MouseEvent move( final TreePanel tp, final Point p ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_MOVED, System.currentTimeMillis(), 0, p.x, p.y, 0, false );
    }

    /** The unrooted layout draws no cells, so nothing in it can be read out -- probed at the points the rectangular
     *  hit-test WOULD use, and then over the whole canvas. */
    private static void unrooted( final boolean[] ok, final TreePanel tp ) {
        layout( tp, "unrooted" );
        // The notional points are computed IN the unrooted layout, from the same column and row hooks the
        // rectangular hit-test reads -- points taken from another layout mean nothing here, and a scan alone left
        // the guard removable without anything noticing (measured).
        final List<Point> notional = new ArrayList<Point>();
        for( int c = 0; c < 6; ++c ) {
            for( int r = 0; r < tp.getPhylogeny().getExternalNodes().size(); ++r ) {
                final Point p = cellPoint( tp, c, r );
                if ( p != null ) {
                    notional.add( p );
                }
            }
        }
        if ( notional.isEmpty() ) {
            fail( ok, "the unrooted probe found no points to test the guard with" );
            return;
        }
        for( final Point p : notional ) {
            if ( tp.annotationCellAt( p.x, p.y ) != null ) {
                fail( ok, "unrooted draws no cells, yet one reads out at " + p );
                return;
            }
        }
        for( int x = 0; x < W; x += 7 ) {
            for( int y = 0; y < H; y += 7 ) {
                if ( tp.annotationCellAt( x, y ) != null ) {
                    fail( ok, "unrooted draws no cells, yet one reads out at " + x + "," + y );
                    return;
                }
            }
        }
    }
}
