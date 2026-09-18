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
import java.awt.Point;
import java.awt.event.InputEvent;
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

/**
 * Dragging an annotation column to a new place: press its header (rectangular, vertical) or its ring (circular) and
 * drag. Pinned here:
 * <ul>
 * <li>{@link AnnotationColumns#specIndicesByColumn} against the REAL constructor -- a merged stacked-bar or pie
 * column is several specs, so a drawn column must move as its whole group;</li>
 * <li>the move itself: a group moves whole, a column dropped where it was does not move, a matrix column moved by
 * hand makes the tab Manual while a colour strip moved around does not, and a shown legend follows its column;</li>
 * <li>the geometry in all three layouts that draw columns (root-left, the clustergram, circular) and none in
 * unrooted, which draws no columns;</li>
 * <li>the whole gesture through a real {@link MouseListener}: a sub-threshold wiggle is still a click (it toggles the
 * legend, moves nothing, pans nothing), a real drag moves the column and swallows the trailing click, marks the tree
 * edited and syncs the Order Matrix Columns radios;</li>
 * <li>the drop marker on screen and never in an export.</li>
 * </ul>
 * Headful; a green no-op when headless.
 */
public final class AnnotationColumnDragTest {

    private static final int W = 1100;
    private static final int H = 800;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "AnnotationColumnDrag: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return dragOk();
    }

    private static boolean fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [AnnotationColumnDragTest] " + msg );
        ok[ 0 ] = false;
        return false;
    }

    // ---- fixture -------------------------------------------------------------------------------------------------

    /** 8 tips: a categorical host, four numeric matrix fields with different patterns, two stacked-bar fields. */
    private static Phylogeny tree() {
        final String[] host = { "cat", "dog", "cat", "dog", "cow", "cat", "dog", "cow" };
        final int[][] m = { { 0, 1, 2, 3, 4, 4, 3, 2 }, { 4, 4, 3, 1, 0, 0, 1, 2 }, { 1, 3, 1, 3, 1, 3, 1, 3 },
                { 2, 2, 4, 4, 0, 0, 2, 2 } };
        final int[][] st = { { 5, 9, 2, 7, 4, 6, 8, 3 }, { 3, 1, 6, 2, 5, 4, 2, 7 } };
        final PhylogenyNode root = new PhylogenyNode();
        final List<PhylogenyNode> tips = new ArrayList<PhylogenyNode>();
        for( int i = 0; i < 8; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "tip_" + i );
            tip.setDistanceToParent( 0.1 + ( 0.02 * i ) );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:host", host[ i ], "", "xsd:string", AppliesTo.NODE ) );
            for( int k = 0; k < 4; ++k ) {
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

    /** What the user sees: each drawn column's header, a merged column as the list of its series. */
    private static List<String> drawn( final TreePanel tp ) {
        final AnnotationColumns cols = new AnnotationColumns( tp.getPhylogeny(), tp.getAnnotationColumnSpecs() );
        final List<String> out = new ArrayList<String>();
        for( int i = 0; i < cols.size(); ++i ) {
            final List<String> stack = cols.stackHeaders( i );
            out.add( stack.isEmpty() ? cols.getColumn( i ).getHeader() : String.join( "+", stack ) );
        }
        return out;
    }

    private static List<String> moved( final List<String> before, final int from, final int slot ) {
        final List<String> out = new ArrayList<String>( before );
        final String c = out.remove( from );
        out.add( ( slot > from ) ? ( slot - 1 ) : slot, c );
        return out;
    }

    // ---- the test ------------------------------------------------------------------------------------------------

    private static boolean dragOk() {
        final boolean[] ok = { true };
        try {
            groupingMatchesConstructor( ok, tree() );
            final Phylogeny phy = tree();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "drag" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    ( ( JFrame ) frame ).setSize( W, H );
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    moveModel( ok, tp );
                    for( final String layout : new String[] { "root-left", "clustergram", "circular" } ) {
                        layout( tp, layout );
                        geometry( ok, tp, layout );
                        gesture( ok, frame, tp, layout );
                    }
                    layout( tp, "unrooted" );
                    // where the rectangular hit-test WOULD put each header -- the scan below never reaches those
                    // points (unrooted lays out elsewhere), so without this the unrooted guard could be deleted unseen
                    for( int c = 0; c < 6; ++c ) {
                        final Point p = tp.annotationColumnGrabPointForTest( c );
                        if ( ( p != null ) && ( tp.annotationColumnGrabbedAt( p.x, p.y ) >= 0 ) ) {
                            fail( ok, "unrooted draws no columns, yet column " + c + " can be grabbed at " + p );
                        }
                    }
                    for( int x = 0; x < W; x += 7 ) {
                        for( int y = 0; y < H; y += 7 ) {
                            if ( tp.annotationColumnGrabbedAt( x, y ) >= 0 ) {
                                fail( ok, "unrooted draws no columns, yet a column can be grabbed at " + x + "," + y );
                                x = W;
                                break;
                            }
                        }
                    }
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

    /**
     * specIndicesByColumn must agree with what the constructor actually draws -- checked on several spec lists: the
     * merged columns split around others, a pie AND a stacked bar together, and a LABEL spec, which draws nothing.
     */
    private static void groupingMatchesConstructor( final boolean[] ok, final Phylogeny phy ) {
        final List<List<AnnotationColumns.ColumnSpec>> cases = new ArrayList<List<AnnotationColumns.ColumnSpec>>();
        cases.add( specs() );
        cases.add( Arrays.asList( spec( "x", AnnotationColumns.Type.STACKED_BAR ), spec( "host", AnnotationColumns.Type.COLOR_STRIP ),
                                  spec( "y", AnnotationColumns.Type.STACKED_BAR ), spec( "m1", AnnotationColumns.Type.MATRIX ) ) );
        cases.add( Arrays.asList( spec( "m1", AnnotationColumns.Type.PIE ), spec( "x", AnnotationColumns.Type.STACKED_BAR ),
                                  spec( "m2", AnnotationColumns.Type.PIE ), spec( "y", AnnotationColumns.Type.STACKED_BAR ),
                                  spec( "host", AnnotationColumns.Type.LABEL ), spec( "m3", AnnotationColumns.Type.MATRIX ) ) );
        for( final List<AnnotationColumns.ColumnSpec> specs : cases ) {
            final AnnotationColumns cols = new AnnotationColumns( phy, specs );
            final List<List<Integer>> groups = AnnotationColumns.specIndicesByColumn( specs );
            if ( groups.size() != cols.size() ) {
                fail( ok, "specIndicesByColumn gives " + groups.size() + " columns, the constructor draws " + cols.size() );
                continue;
            }
            for( int i = 0; i < groups.size(); ++i ) {
                final List<String> want = new ArrayList<String>();
                for( final int k : groups.get( i ) ) {
                    want.add( PropertyColorScheme.displayName( specs.get( k )._ref ) );
                }
                final List<String> got = cols.stackHeaders( i ).isEmpty()
                        ? Arrays.asList( cols.getColumn( i ).getHeader() ) : cols.stackHeaders( i );
                if ( !want.equals( got ) ) {
                    fail( ok, "drawn column " + i + " is built from " + got + ", specIndicesByColumn says " + want );
                }
            }
        }
    }

    /** The move, on the model: groups move whole, in-place drops do nothing, Manual only when the matrix order changes,
     *  the legend follows its column. */
    private static void moveModel( final boolean[] ok, final TreePanel tp ) {
        tp.setAnnotationColumns( specs() );
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.TABLE );
        final List<String> start = drawn( tp );
        if ( start.size() != 6 ) {
            fail( ok, "the fixture should draw 6 columns (7 specs, one merged pair), draws " + start );
            return;
        }
        // the merged stacked bar (drawn column 4) moves to the front AS A WHOLE; the matrix order is untouched
        if ( !tp.moveAnnotationColumn( 4, 0 ) || !drawn( tp ).equals( moved( start, 4, 0 ) ) ) {
            fail( ok, "a merged column must move as a whole: " + drawn( tp ) );
        }
        if ( tp.getMatrixColumnOrder() != MatrixColumnOrder.Mode.TABLE ) {
            fail( ok, "moving a non-matrix column must not make the tab Manual, is " + tp.getMatrixColumnOrder() );
        }
        // dropped where it was: slot == from, and slot == from + 1, move nothing
        final List<String> before = drawn( tp );
        if ( tp.moveAnnotationColumn( 2, 2 ) || tp.moveAnnotationColumn( 2, 3 ) || !drawn( tp ).equals( before ) ) {
            fail( ok, "a column dropped where it was must not move" );
        }
        // a MATRIX column moved past its neighbours: the matrix order changes -> Manual
        final int m1 = before.indexOf( "M1" );
        final int after_m3 = before.indexOf( "M3" ) + 1;
        if ( ( m1 < 0 ) || ( after_m3 < 1 ) ) {
            fail( ok, "fixture headers not as expected: " + before );
            return;
        }
        if ( !tp.moveAnnotationColumn( m1, after_m3 ) || !drawn( tp ).equals( moved( before, m1, after_m3 ) ) ) {
            fail( ok, "the matrix column did not land where it was dropped: " + drawn( tp ) );
        }
        if ( tp.getMatrixColumnOrder() != MatrixColumnOrder.Mode.MANUAL ) {
            fail( ok, "a matrix column moved by hand must make the tab Manual, is " + tp.getMatrixColumnOrder() );
        }
        // the legend shown for a column follows it, when it moves and when a neighbour moves past it
        final int host = drawn( tp ).indexOf( "Host" );
        tp.setFocusedAnnotationColumn( host );
        tp.moveAnnotationColumn( drawn( tp ).size() - 1, 0 ); // the last column jumps in front of Host
        if ( ( tp.focusedAnnotationColumnForTest() < 0 )
                || !"Host".equals( drawn( tp ).get( tp.focusedAnnotationColumnForTest() ) ) ) {
            fail( ok, "the shown legend must stay on Host when another column moves past it" );
        }
        tp.moveAnnotationColumn( drawn( tp ).indexOf( "Host" ), drawn( tp ).size() ); // Host itself to the end
        if ( ( tp.focusedAnnotationColumnForTest() < 0 )
                || !"Host".equals( drawn( tp ).get( tp.focusedAnnotationColumnForTest() ) ) ) {
            fail( ok, "the shown legend must follow Host when Host moves" );
        }
        if ( tp.focusedAnnotationColumnForTest() >= 0 ) {
            tp.setFocusedAnnotationColumn( tp.focusedAnnotationColumnForTest() ); // toggle it off again
        }
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

    /** A paint through the SCREEN path (paintComponent). */
    private static void screenPaint( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.paint( g );
        g.dispose();
    }

    /** Every drawn column can be grabbed where it is drawn, and every slot position maps back to its own slot. */
    private static void geometry( final boolean[] ok, final TreePanel tp, final String layout ) {
        final int n = drawn( tp ).size();
        for( int c = 0; c < n; ++c ) {
            final Point p = tp.annotationColumnGrabPointForTest( c );
            if ( ( p == null ) || ( tp.annotationColumnGrabbedAt( p.x, p.y ) != c ) ) {
                fail( ok, layout + ": column " + c + " cannot be grabbed where it is drawn (" + p + ")" );
            }
        }
        if ( "circular".equals( layout ) ) {
            // off the 3-o'clock spoke, where x-offset and radius are NOT the same number
            for( final double angle : new double[] { 0.7, 1.9, 3.0, 4.4, 5.6 } ) {
                for( int s = 0; s <= n; ++s ) {
                    final Point p = tp.annotationColumnSlotPointForTest( s, angle );
                    if ( ( p == null ) || ( tp.annotationColumnInsertionSlotAt( p.x, p.y ) != s ) ) {
                        fail( ok, "circular: at angle " + angle + " the point at slot " + s + " maps elsewhere" );
                    }
                }
            }
        }
        Point prev = null;
        for( int s = 0; s <= n; ++s ) {
            final Point p = tp.annotationColumnSlotPointForTest( s );
            if ( ( p == null ) || ( tp.annotationColumnInsertionSlotAt( p.x, p.y ) != s ) ) {
                fail( ok, layout + ": the point at slot " + s + " maps to slot "
                        + ( ( p == null ) ? "none" : tp.annotationColumnInsertionSlotAt( p.x, p.y ) ) );
            }
            if ( ( prev != null ) && ( p != null ) && prev.equals( p ) ) {
                fail( ok, layout + ": slots " + ( s - 1 ) + " and " + s + " sit at the same point" );
            }
            prev = p;
        }
    }

    private static MouseEvent ev( final TreePanel tp, final int id, final Point p, final boolean button_down ) {
        return new MouseEvent( tp, id, System.currentTimeMillis(), button_down ? InputEvent.BUTTON1_DOWN_MASK : 0, p.x,
                               p.y, 1, false, MouseEvent.BUTTON1 );
    }

    /** The whole gesture through a real MouseListener, in one layout. */
    private static void gesture( final boolean[] ok, final MainFrame frame, final TreePanel tp, final String layout ) {
        final MouseListener ml = new MouseListener( tp );
        final List<String> start = drawn( tp );
        // (1) a press + sub-threshold wiggle + release is a CLICK: the legend toggles, nothing moves, nothing pans
        final int host = start.indexOf( "Host" );
        final Point hp = tp.annotationColumnGrabPointForTest( host );
        ml.mousePressed( ev( tp, MouseEvent.MOUSE_PRESSED, hp, true ) );
        ml.mouseDragged( ev( tp, MouseEvent.MOUSE_DRAGGED, new Point( hp.x + 2, hp.y + 1 ), true ) );
        if ( ml.isDraggingForTest() || ( tp.columnDragSlotForTest() >= 0 ) ) {
            fail( ok, layout + ": a sub-threshold wiggle on a header must neither pan nor start a column drag" );
        }
        ml.mouseReleased( ev( tp, MouseEvent.MOUSE_RELEASED, hp, false ) );
        ml.mouseClicked( ev( tp, MouseEvent.MOUSE_CLICKED, hp, false ) );
        if ( !tp.hasFocusedAnnotationColumn() || !drawn( tp ).equals( start ) ) {
            fail( ok, layout + ": a click on a header must still show its legend and move nothing" );
        }
        tp.setFocusedAnnotationColumn( tp.focusedAnnotationColumnForTest() ); // toggle the legend off again
        // (2) a real drag: the matrix column M2 goes to the very end
        final int from = start.indexOf( "M2" );
        final int slot = start.size();
        final Point grab = tp.annotationColumnGrabPointForTest( from );
        final Point drop = tp.annotationColumnSlotPointForTest( slot );
        tp.setEdited( false );
        ml.mousePressed( ev( tp, MouseEvent.MOUSE_PRESSED, grab, true ) );
        ml.mouseDragged( ev( tp, MouseEvent.MOUSE_DRAGGED, drop, true ) );
        if ( tp.columnDragSlotForTest() != slot ) {
            fail( ok, layout + ": dragging to the end should mark slot " + slot + ", marks " + tp.columnDragSlotForTest() );
        }
        if ( ml.isDraggingForTest() ) {
            fail( ok, layout + ": a column drag must not also pan the view" );
        }
        // the marker: drawn on screen, never into an export
        final int before_export = tp._column_drag_marker_paints;
        AptxUtil.renderPhylogenyToImage( W, H, tp, frame.getOptions(), false, 1, false );
        if ( tp._column_drag_marker_paints != before_export ) {
            fail( ok, layout + ": the drop marker must never be drawn into an exported image" );
        }
        screenPaint( tp );
        if ( tp._column_drag_marker_paints <= before_export ) {
            fail( ok, layout + ": the drop marker must be drawn on screen during a drag" );
        }
        ml.mouseReleased( ev( tp, MouseEvent.MOUSE_RELEASED, drop, false ) );
        // the click that trails the drag, landing ON a header (a drag along the header row ends on one): it must be
        // swallowed, not read as a click on that header -- at the drop point, inside the cells, it would prove nothing
        final Point on_header = tp.annotationColumnGrabPointForTest( 0 );
        ml.mouseClicked( ev( tp, MouseEvent.MOUSE_CLICKED, on_header, false ) );
        if ( !drawn( tp ).equals( moved( start, from, slot ) ) ) {
            fail( ok, layout + ": after the drop the columns should read " + moved( start, from, slot ) + ", read "
                    + drawn( tp ) );
        }
        if ( tp.hasFocusedAnnotationColumn() ) {
            fail( ok, layout + ": the click that trails a real drag must be swallowed, not toggle a legend" );
        }
        if ( tp.isDraggingAnnotationColumn() || ( tp.columnDragSlotForTest() >= 0 ) ) {
            fail( ok, layout + ": the drag must be over after the release" );
        }
        if ( !tp.isEdited() ) {
            fail( ok, layout + ": a moved column must mark the tree edited (the arrangement is saved with it)" );
        }
        if ( ( tp.getMatrixColumnOrder() != MatrixColumnOrder.Mode.MANUAL )
                || !frame._matrix_order_items.get( MatrixColumnOrder.Mode.MANUAL ).isSelected() ) {
            fail( ok, layout + ": a matrix column dragged by hand must make the tab Manual, radios included" );
        }
    }
}
