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
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * The node click / hover target ({@link TreePanel#findNode}, {@link TreePanel#nodeHitHalfSize}): a square around
 * the node whose half side is the drawn mark's half plus WIGGLE, never under {@link TreePanel#MIN_NODE_HIT_HALF_SIZE}
 * (4: 9 x 9 px), with no ceiling; where squares overlap the NEAREST node wins, a tie going to the first in postorder.
 * {@link TreePanel#findBranch} shares the floor ({@code hitReach}). Three parts: hand-placed coordinates (every
 * boundary to the pixel); the demo node-click-targets.xml in a real window at the size it opens (the layout a user
 * gets); and the same demo in all five layouts (rectangular root-left / top / bottom, circular, unrooted), clicked at
 * DEVICE points (a vertical orientation rotates the canvas). Each part pins its layout and orientation: a saved
 * Root at Top would otherwise rotate every click.
 */
public final class NodeHitTargetTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeHitTarget: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // needs a display
        }
        try {
            return placed() & demo() & parity();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    // ---- hand-placed: a (100,100) b (100,112) c (100,200) d (300,300) g (312,300); e collapsed at (300,100) ----
    private static boolean placed() throws Exception {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = tip( root, "a" ), b = tip( root, "b" ), c = tip( root, "c" ), d = tip( root, "d" );
        final PhylogenyNode g = tip( root, "g" ); // beside d on the same row: the horizontal distance decides
        final PhylogenyNode e = new PhylogenyNode();
        e.setName( "e" );
        root.addAsChild( e );
        final PhylogenyNode e1 = tip( e, "e1" ), e2 = tip( e, "e2" );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "hit targets" ) );
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final Options o = tp.getOptions();
                final Runnable place = () -> {
                    at( root, 20, 400 );
                    at( a, 100, 100 );
                    at( b, 100, 112 );
                    at( c, 100, 200 );
                    at( d, 300, 300 );
                    at( g, 312, 300 );
                    at( e, 300, 100 );
                    at( e1, 100, 118 ); // stale coords of a tip hidden under the collapsed e
                    at( e2, 500, 500 );
                };
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.getControlPanel().setCheckbox( DisplayOption.USE_STYLE, false );
                place.run();
                // ---- the floor: the default size's 9 x 9 px target, also at node size 0 and 1 ----
                o.setDefaultNodeShapeSize( (short) 5 );
                ok[ 0 ] &= ck( tp.nodeHitHalfSize( c ) == TreePanel.MIN_NODE_HIT_HALF_SIZE,
                               "control: the default node size (5) gives the 9 x 9 px target" );
                ok[ 0 ] &= ck( ( tp.findNode( 104, 204 ) == c ) && ( tp.findNode( 105, 200 ) == null )
                               && ( tp.findNode( 100, 205 ) == null ), "control: at size 5 the corner hits, a pixel beyond misses" );
                for( final short size : new short[] { 0, 1 } ) {
                    o.setDefaultNodeShapeSize( size );
                    ok[ 0 ] &= ck( tp.effectiveNodeHalfBoxSize( c ) == 0, "precondition: at size " + size + " the mark's half is 0" );
                    ok[ 0 ] &= ck( ( tp.nodeHitHalfSize( c ) == 4 ) && ( tp.findNode( 104, 204 ) == c )
                                   && ( tp.findNode( 105, 200 ) == null ),
                                   "node size " + size + ": still the 9 x 9 px target (the corner hits, a pixel beyond misses)" );
                }
                // ---- a branch shares the floor (its reach is strict: closer than 4 px), off every node ----
                o.setDefaultNodeShapeSize( (short) 0 );
                ok[ 0 ] &= ck( ( tp.findNode( 60, 203 ) == null ) && ( tp.findBranch( 60, 203 ) == c )
                               && ( tp.findBranch( 60, 204 ) == null ),
                               "node size 0: c's branch still answers 3 px off its leg, not 4" );
                o.setDefaultNodeShapeSize( (short) 5 );
                ok[ 0 ] &= ck( ( tp.findBranch( 60, 204 ) == c ) && ( tp.findBranch( 60, 205 ) == null ),
                               "control: at size 5 a branch answers within 4.5 px" );
                // ---- no ceiling: the mark plus the wiggle, to the pixel ----
                o.setDefaultNodeShapeSize( (short) 40 );
                ok[ 0 ] &= ck( ( tp.nodeHitHalfSize( c ) == 22 ) && ( tp.findNode( 122, 200 ) == c )
                               && ( tp.findNode( 123, 200 ) == null ), "node size 40: a 45 x 45 px target" );
                o.setDefaultNodeShapeSize( (short) 100 );
                ok[ 0 ] &= ck( ( tp.nodeHitHalfSize( c ) == 52 ) && ( tp.findNode( 152, 200 ) == c )
                               && ( tp.findNode( 153, 200 ) == null ), "node size 100: a 105 x 105 px target, no ceiling" );
                // a node's own style size, with Use Visual Styles on (off: the global size)
                o.setDefaultNodeShapeSize( (short) 5 );
                final NodeVisualData big = new NodeVisualData();
                big.setSize( 60f );
                d.getNodeData().setNodeVisualData( big );
                ok[ 0 ] &= ck( tp.findNode( 332, 300 ) == null, "control: Use Visual Styles off, d's style size is not used" );
                tp.getControlPanel().setCheckbox( DisplayOption.USE_STYLE, true );
                place.run();
                ok[ 0 ] &= ck( ( tp.findNode( 332, 300 ) == d ) && ( tp.findNode( 333, 300 ) == null ),
                               "a style size of 60: a 65 x 65 px target" );
                tp.getControlPanel().setCheckbox( DisplayOption.USE_STYLE, false );
                place.run();
                // ---- overlapping squares: the nearest node wins ----
                o.setDefaultNodeShapeSize( (short) 40 ); // half 22: a's square covers b's node (12 px away)
                ok[ 0 ] &= ck( postorderIndex( phy, a ) < postorderIndex( phy, b ),
                               "precondition: a comes before b in postorder (the old rule's winner)" );
                ok[ 0 ] &= ck( Math.abs( b.getYcoord() - a.getYcoord() ) <= tp.nodeHitHalfSize( a ),
                               "precondition: a's square covers b's node" );
                ok[ 0 ] &= ck( tp.findNode( 100, 112 ) == b, "a click right on b opens b, not a (whose square covers it)" );
                ok[ 0 ] &= ck( tp.findNode( 100, 104 ) == a, "4 px from a, 8 from b: a" );
                ok[ 0 ] &= ck( tp.findNode( 100, 107 ) == b, "7 px from a, 5 from b: b" );
                ok[ 0 ] &= ck( tp.findNode( 100, 106 ) == a, "6 px from each: the first in postorder, a" );
                ok[ 0 ] &= ck( ( postorderIndex( phy, d ) < postorderIndex( phy, g ) ) && ( tp.findNode( 310, 300 ) == g )
                               && ( tp.findNode( 304, 300 ) == d ), "one row, two nodes: 2 px from g opens g, 4 px from d opens d" );
                // ---- a node hidden under a collapsed clade is never a target, even the nearest ----
                e.setCollapse( true );
                ok[ 0 ] &= ck( tp.findNode( 100, 118 ) == b, "e1 (hidden, 0 px away) is skipped: b (6 px)" );
                ok[ 0 ] &= ck( tp.findNode( 300, 100 ) == e, "the collapsed clade's own node is a target" );
                e.setCollapse( false );
                ok[ 0 ] &= ck( tp.findNode( 100, 118 ) == e1, "control: not collapsed, e1 is the nearest" );
                o.setDefaultNodeShapeSize( (short) 5 );
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    // ---- the demo, laid out by the app at the size its window opens ----
    private static boolean demo() throws Exception {
        final File f = new File( System.getProperty( "user.dir" ), "forester/demo/node-click-targets.xml" );
        final Phylogeny phy = FigureRenderer.readTrees( f )[ 0 ];
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "click targets" ) );
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> {
                final JFrame frame = (JFrame) mf[ 0 ];
                ok[ 0 ] &= ck( ( frame.getWidth() == 900 ) && ( frame.getHeight() <= 900 ),
                               "the window opens 900 px wide and at most 900 px tall, as the demo says ("
                                       + frame.getWidth() + " x " + frame.getHeight() + ")" );
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                ok[ 0 ] &= ck( tp.shows( DisplayOption.USE_STYLE ), "Use Visual Styles is on as the file opens" );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                mf[ 0 ].showWhole();
                tp.paintImmediately( 0, 0, tp.getWidth(), tp.getHeight() );
                final Phylogeny p = tp.getPhylogeny();
                int covered = 0, right = 0, ring = 0, pairs = 0, old_wrong = 0;
                for( int k = 1; k <= 20; ++k ) {
                    final String n = String.format( "%02d", k );
                    final PhylogenyNode large = p.getNodes( "large_" + n ).get( 0 );
                    final PhylogenyNode plain = p.getNodes( "plain_" + n ).get( 0 );
                    ++pairs;
                    final int half = tp.nodeHitHalfSize( large );
                    if ( ( Math.abs( plain.getYcoord() - large.getYcoord() ) <= half )
                            && ( Math.abs( plain.getXcoord() - large.getXcoord() ) <= half ) ) {
                        ++covered;
                    }
                    if ( postorderIndex( p, large ) < postorderIndex( p, plain ) ) {
                        ++old_wrong; // the old rule (first in postorder) would have answered the large one
                    }
                    if ( ( tp.findNode( Math.round( plain.getXcoord() ), Math.round( plain.getYcoord() ) ) == plain )
                            && ( tp.findNode( Math.round( large.getXcoord() ), Math.round( large.getYcoord() ) ) == large ) ) {
                        ++right;
                    }
                    // inside the ring, sideways from its node: nearer the large one than any other
                    if ( tp.findNode( Math.round( large.getXcoord() ) + ( half - 4 ), Math.round( large.getYcoord() ) ) == large ) {
                        ++ring;
                    }
                }
                ok[ 0 ] &= ck( tp.nodeHitHalfSize( p.getNodes( "large_01" ).get( 0 ) ) == 26,
                               "a 48 px ring: a 53 x 53 px target (half 26)" );
                ok[ 0 ] &= ck( covered == pairs, "precondition: every ring's square covers its plain neighbour's node ("
                        + covered + " of " + pairs + ", rows "
                        + Math.abs( p.getNodes( "plain_01" ).get( 0 ).getYcoord() - p.getNodes( "large_01" ).get( 0 ).getYcoord() )
                        + " px apart)" );
                ok[ 0 ] &= ck( old_wrong == 10, "precondition: in 10 of the 20 pairs the large tip comes first in postorder ("
                        + old_wrong + ")" );
                ok[ 0 ] &= ck( right == pairs, "a click right on each tip opens that tip (" + right + " of " + pairs + ")" );
                ok[ 0 ] &= ck( ring == pairs, "a click inside a ring, away from other nodes, opens the large one (" + ring
                        + " of " + pairs + ")" );
                for( final String z : new String[] { "no_mark_1", "no_mark_2" } ) {
                    final PhylogenyNode t = p.getNodes( z ).get( 0 );
                    // node coordinates are fractional: the last pixel inside the square, and the first beyond it
                    final int x4 = (int) Math.floor( t.getXcoord() + 4 ), y4 = (int) Math.floor( t.getYcoord() + 4 );
                    final int y = Math.round( t.getYcoord() );
                    ok[ 0 ] &= ck( ( tp.effectiveNodeHalfBoxSize( t ) == 0 ) && ( tp.nodeHitHalfSize( t ) == 4 )
                                   && ( tp.findNode( x4, y4 ) == t ) && ( tp.findNode( x4 + 1, y ) != t ),
                                   z + ": node size 0, still a 9 x 9 px target" );
                }
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    // ---- all five layouts, at device points: the nearest wins, and a mark of size 0 keeps the 9 x 9 px floor ----
    private static boolean parity() throws Exception {
        final File f = new File( System.getProperty( "user.dir" ), "forester/demo/node-click-targets.xml" );
        final Phylogeny phy = FigureRenderer.readTrees( f )[ 0 ];
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "click targets, every layout" ) );
        final boolean[] ok = { true };
        final Object[][] layouts = {
                { "rectangular root-left", Options.TREE_ORIENTATION.ROOT_LEFT, Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR },
                { "rectangular root-top", Options.TREE_ORIENTATION.ROOT_TOP, Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR },
                { "rectangular root-bottom", Options.TREE_ORIENTATION.ROOT_BOTTOM, Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR },
                { "circular", Options.TREE_ORIENTATION.ROOT_LEFT, Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR },
                { "unrooted", Options.TREE_ORIENTATION.ROOT_LEFT, Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED } };
        try {
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final Phylogeny p = tp.getPhylogeny();
                for( final Object[] l : layouts ) {
                    final String name = (String) l[ 0 ];
                    tp.setTreeOrientation( (Options.TREE_ORIENTATION) l[ 1 ] );
                    tp.setPhylogenyGraphicsType( (Options.PHYLOGENY_GRAPHICS_TYPE) l[ 2 ] );
                    final int w = 450, h = 450; // small enough that the rings overlap their neighbours in every layout
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final java.awt.image.BufferedImage img = new java.awt.image.BufferedImage( w, h,
                            java.awt.image.BufferedImage.TYPE_INT_RGB );
                    final java.awt.Graphics2D g = img.createGraphics();
                    tp.paint( g ); // a screen paint: lays the tree out for the panel, and the rotation the hit test reads
                    g.dispose();
                    int covered = 0, right = 0;
                    for( int k = 1; k <= 20; ++k ) {
                        final String n = String.format( "%02d", k );
                        final PhylogenyNode large = p.getNodes( "large_" + n ).get( 0 );
                        final PhylogenyNode plain = p.getNodes( "plain_" + n ).get( 0 );
                        final int half = tp.nodeHitHalfSize( large );
                        if ( ( Math.abs( plain.getYcoord() - large.getYcoord() ) <= half )
                                && ( Math.abs( plain.getXcoord() - large.getXcoord() ) <= half ) ) {
                            ++covered;
                        }
                        if ( ( hit( tp, plain, 0, 0 ) == plain ) && ( hit( tp, large, 0, 0 ) == large ) ) {
                            ++right;
                        }
                    }
                    ok[ 0 ] &= ck( covered >= 10, name + ": precondition: rings cover their neighbours' nodes (" + covered
                            + " of 20 pairs)" );
                    ok[ 0 ] &= ck( right == 20, name + ": a click right on each tip opens that tip (" + right + " of 20)" );
                    for( final String z : new String[] { "no_mark_1", "no_mark_2" } ) {
                        final PhylogenyNode t = p.getNodes( z ).get( 0 );
                        // outward, from its parent through the tip: no other node lies that way
                        final double ux = t.getXcoord() - t.getParent().getXcoord();
                        final double uy = t.getYcoord() - t.getParent().getYcoord();
                        final double len = Math.sqrt( ( ux * ux ) + ( uy * uy ) );
                        ok[ 0 ] &= ck( ( len > 0 ) && ( hit( tp, t, 3 * ux / len, 3 * uy / len ) == t )
                                       && ( hit( tp, t, 6 * ux / len, 6 * uy / len ) != t ),
                                       name + ", " + z + ": node size 0 answers 3 px out, not 6" );
                    }
                }
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    /** findNode at the DEVICE point of {@code n}'s logical position moved by ({@code dx},{@code dy}). */
    private static PhylogenyNode hit( final TreePanel tp, final PhylogenyNode n, final double dx, final double dy ) {
        final java.awt.geom.Point2D.Double d = tp.toDevicePointForTest( n.getXcoord() + dx, n.getYcoord() + dy );
        return tp.findNode( (int) Math.round( d.x ), (int) Math.round( d.y ) );
    }

    private static PhylogenyNode tip( final PhylogenyNode parent, final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        parent.addAsChild( n );
        return n;
    }

    private static void at( final PhylogenyNode n, final float x, final float y ) {
        n.setXcoord( x );
        n.setYcoord( y );
    }

    private static int postorderIndex( final Phylogeny phy, final PhylogenyNode target ) {
        final List<PhylogenyNode> order = new ArrayList<>();
        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
            order.add( it.next() );
        }
        return order.indexOf( target );
    }

    private static boolean ck( final boolean condition, final String what ) {
        if ( !condition ) {
            System.out.println( "  [NodeHitTargetTest] " + what );
        }
        return condition;
    }

    private NodeHitTargetTest() {
        // not instantiable
    }
}
