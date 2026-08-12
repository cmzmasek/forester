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
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Radial (circular/unrooted) interaction parity, increment 2: asserts that branch-click clade selection works radially
 * ({@link TreePanel#findBranch} now hit-tests the radial branch geometry), that selecting off a found radial branch
 * populates the selection, and that the "Pulse Found Nodes" halo now renders in circular AND unrooted. The findBranch
 * checks compute a point ON a node's drawn branch from that node's coords and hit-test against the same coords, so they
 * are independent of the display-density-dependent node->device transform. Headful; a green no-op when headless.
 */
public final class RadialInteractionTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RadialInteraction: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return findBranchOk() && haloOk();
    }

    /** findBranch returns the node whose incoming branch is under a point placed ON that branch -- for a LEAF and an
     *  INTERNAL clade root, in CIRCULAR and UNROOTED; and selecting off a found internal branch selects its tips. */
    private static boolean findBranchOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 820, h = 820;
            o.setGraphicsExportWhiteBackground( true );
            frame.showWhole();
            tp.setSize( w, h );
            tp.getControlPanel().setActionWhenNodeClicked( ControlPanel.NodeClickAction.SELECT_NODES );
            final PhylogenyNode leaf = firstMatch( tp.getPhylogeny(), true );
            final PhylogenyNode clade = firstMatch( tp.getPhylogeny(), false );
            if ( ( leaf == null ) || ( clade == null ) ) {
                fail( ok, "precondition: need a leaf and an internal clade with a non-root parent" );
                return;
            }
            for ( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                    PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                tp.setPhylogenyGraphicsType( gt );
                tp.calcParametersForPainting( w, h );
                AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ); // assigns radial node coords
                final PhylogenyNode root = tp.getPhylogeny().getRoot();
                if ( ( Math.abs( leaf.getXcoord() - root.getXcoord() )
                        + Math.abs( leaf.getYcoord() - root.getYcoord() ) ) < 5 ) {
                    fail( ok, "precondition: radial layout did not assign meaningful node coords in " + gt );
                    return;
                }
                // a click ON the leaf's own incoming branch resolves to the leaf
                final double[] lp = branchPoint( leaf, root, gt );
                final PhylogenyNode leaf_hit = tp.findBranch( (int) Math.round( lp[ 0 ] ), (int) Math.round( lp[ 1 ] ) );
                if ( leaf_hit != leaf ) {
                    fail( ok, "findBranch on a leaf's branch must return that leaf in " + gt + " (got " + leaf_hit
                            + ")" );
                }
                // a click ON the clade root's branch resolves to the clade root, and selecting it selects its tips
                final double[] cp = branchPoint( clade, root, gt );
                final PhylogenyNode clade_hit = tp.findBranch( (int) Math.round( cp[ 0 ] ), (int) Math.round( cp[ 1 ] ) );
                if ( clade_hit != clade ) {
                    fail( ok, "findBranch on an internal clade's branch must return that clade in " + gt + " (got "
                            + clade_hit + ")" );
                    continue;
                }
                // hovering that branch (Select-Node(s) mode) sets the hover-preview target -- the source of the
                // translucent on-screen preview, which is now dispatched in the radial branches (was dead before).
                // Assert it resolves to THAT clade as a subtree (not just any node) -- proving the radial findBranch
                // branch-hover path, not the node-hover fallthrough (mouseMoved tries findNode first).
                tp.clearHoverPreview();
                tp.mouseMoved( new java.awt.event.MouseEvent( tp, java.awt.event.MouseEvent.MOUSE_MOVED,
                        System.currentTimeMillis(), 0, (int) Math.round( cp[ 0 ] ), (int) Math.round( cp[ 1 ] ), 0,
                        false ) );
                if ( ( tp.hoverNodeForTest() != clade ) || !tp.hoverIsSubtreeForTest() ) {
                    fail( ok, "hovering a radial branch must set that clade as the subtree hover-preview target in "
                            + gt + " (got " + tp.hoverNodeForTest() + ", subtree=" + tp.hoverIsSubtreeForTest() + ")" );
                }
                tp.clearHoverPreview();
                tp.setFoundNodes0( null );
                tp.selectSubtreeTips( clade_hit );
                final Set<Long> found = tp.getFoundNodes0();
                if ( ( found == null ) || ( found.size() != clade.getNumberOfExternalNodes() ) ) {
                    fail( ok, "branch-click selecting a clade must select all its tips in " + gt + " (selected "
                            + ( found == null ? "null" : found.size() ) + " of " + clade.getNumberOfExternalNodes()
                            + ")" );
                }
                tp.setFoundNodes0( null );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** The "Pulse Found Nodes" halo now renders in circular AND unrooted: with a found node set, turning the pulse ON
     *  adds a translucent found-colour disc (a static glow in an export) -> more tinted ink than with it OFF. */
    private static boolean haloOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 780, h = 780;
            o.setGraphicsExportWhiteBackground( true );
            frame.showWhole();
            tp.setSize( w, h );
            final PhylogenyNode leaf = firstMatch( tp.getPhylogeny(), true );
            if ( leaf == null ) {
                fail( ok, "precondition: need a leaf to mark as found" );
                return;
            }
            final Set<Long> found = new HashSet<Long>();
            found.add( leaf.getId() );
            tp.setFoundNodes0( found );
            for ( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                    PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                tp.setPhylogenyGraphicsType( gt );
                tp.calcParametersForPainting( w, h );
                o.setPulseFoundNodes( false );
                final int no_halo = countTinted( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                o.setPulseFoundNodes( true );
                final int with_halo = countTinted( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                if ( with_halo <= ( no_halo + 60 ) ) {
                    fail( ok, "the found-node pulse halo must render in " + gt + " (tinted ink " + with_halo + " vs "
                            + no_halo + ")" );
                }
            }
            o.setPulseFoundNodes( false );
            tp.setFoundNodes0( null );
        }, ok );
        return ok[ 0 ];
    }

    /** A point on {@code node}'s drawn incoming branch. UNROOTED = the midpoint of the straight parent->node line,
     *  which is EXACTLY the segment paintUnrooted draws (drawLine(parent, node)) -- an independent paint-anchored
     *  oracle. CIRCULAR = the midpoint of the radial leg (node -> the point on the parent's radius at the node's
     *  angle); this mirrors both findBranch's reconstruction AND paintBranchCircular's drawn leg (verified equal:
     *  the node sits at root + r*(cos,sin), so atan2(node-root) recovers its stored angle and the inward point matches
     *  drawLine's), so it exercises dispatch + coordinate-space + tolerance + collapse-skip rather than being a fully
     *  independent geometry oracle. */
    private static double[] branchPoint( final PhylogenyNode node, final PhylogenyNode root,
            final PHYLOGENY_GRAPHICS_TYPE gt ) {
        final PhylogenyNode p = node.getParent();
        if ( gt == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
            final double angle = Math.atan2( node.getYcoord() - root.getYcoord(), node.getXcoord() - root.getXcoord() );
            final double pdx = p.getXcoord() - root.getXcoord(), pdy = p.getYcoord() - root.getYcoord();
            final double parent_radius = Math.sqrt( ( pdx * pdx ) + ( pdy * pdy ) );
            final double inward_x = root.getXcoord() + ( Math.cos( angle ) * parent_radius );
            final double inward_y = root.getYcoord() + ( Math.sin( angle ) * parent_radius );
            return new double[] { ( node.getXcoord() + inward_x ) / 2.0, ( node.getYcoord() + inward_y ) / 2.0 };
        }
        return new double[] { ( node.getXcoord() + p.getXcoord() ) / 2.0, ( node.getYcoord() + p.getYcoord() ) / 2.0 };
    }

    /** First external node ({@code leaf}=true) or first internal node ({@code leaf}=false) with a non-root, non-collapsed
     *  parent, itself not collapsed; internal must have >= 2 tips (so its branch is a real clade selection). */
    private static PhylogenyNode firstMatch( final Phylogeny phy, final boolean leaf ) {
        for ( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || n.isCollapse() || n.getParent().isRoot() || n.getParent().isCollapse() ) {
                continue;
            }
            if ( leaf && n.isExternal() ) {
                return n;
            }
            if ( !leaf && !n.isExternal() && ( n.getNumberOfExternalNodes() >= 2 ) ) {
                return n;
            }
        }
        return null;
    }

    private interface FrameBody {
        void run( MainFrame frame, TreePanel tp, Options o ) throws Exception;
    }

    private static void withFrame( final String demo, final FrameBody body, final boolean[] ok ) {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
            if ( !file.exists() ) {
                fail( ok, "demo tree missing: " + file.getAbsolutePath() );
                return;
            }
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                    .create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "radial" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    body.run( frame, frame.getMainPanel().getCurrentTreePanel(), frame.getOptions() );
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
    }

    /** Count of tinted (chromatic) pixels -- a low saturation threshold so the translucent halo disc registers. */
    private static int countTinted( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final float[] hsb = java.awt.Color.RGBtoHSB( ( rgb >> 16 ) & 0xFF, ( rgb >> 8 ) & 0xFF, rgb & 0xFF,
                        null );
                if ( ( hsb[ 1 ] >= 0.12f ) && ( hsb[ 2 ] >= 0.30f ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [RadialInteractionTest] " + msg );
        ok[ 0 ] = false;
    }

    private RadialInteractionTest() {
    }
}
