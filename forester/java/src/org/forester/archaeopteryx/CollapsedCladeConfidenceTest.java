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
import java.awt.image.BufferedImage;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Integration test that the support/confidence value on the branch leading to a <em>collapsed</em> clade is
 * drawn -- a collapsed clade is still an internal node, but the paint path used to {@code return} before the
 * confidence label was drawn, so collapsing a subtree silently dropped its branch support.
 * <p>
 * The tree is rendered offscreen with a deliberately unique confidence color (magenta) on a white background
 * with black branches, so any magenta pixel found on a node's incoming-branch region can only be a confidence
 * label. The test asserts magenta appears on the collapsed clade's incoming branch (the fix) and on a
 * non-collapsed sibling's branch (a baseline that the render path works at all), and that neither appears when
 * "Confidence Values" is turned off (so the magenta really is the gated confidence label).
 * <p>
 * Guarded to a no-op on a headless box (needs FlatLaf via {@code createInstance}); run standalone or as part of
 * the non-headless suite.
 */
public final class CollapsedCladeConfidenceTest {

    private static final int W = 800;
    private static final int H = 600;
    // bootstrap-style support (scale 0-100 -> ceiling 100); both clear the default cutoff (fraction 0.5 -> 50)
    private static final double SUPPORT_AB = 95.0;
    private static final double SUPPORT_CD = 80.0;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CollapsedCladeConfidence: " + ( ok ? "OK." : "FAILED." ) );
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
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { twoClades() }, conf,
                                                                        "conf" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    if ( mp.getControlPanel().getWriteConfidenceCb() == null ) {
                        System.out.println( "  [CollapsedCladeConfidenceTest] no confidence control -- cannot test" );
                        ok[ 0 ] = false;
                        ( (JFrame) mf[ 0 ] ).dispose();
                        return;
                    }
                    // unique colors so a confidence label is the only possible magenta on the canvas
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, Color.WHITE );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BRANCH, Color.BLACK );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.CONFIDENCE, new Color( 255, 0, 255 ) );
                    tp.getTreeColorSet().setColorSchema( 0 ); // propagate scheme[0] overrides into the live fields

                    final PhylogenyNode ab = internalWithFirstChild( tp.getPhylogeny(), "A" );
                    final PhylogenyNode cd = internalWithFirstChild( tp.getPhylogeny(), "C" );
                    if ( ( ab == null ) || ( cd == null ) ) {
                        System.out.println( "  [CollapsedCladeConfidenceTest] test tree shape unexpected" );
                        ok[ 0 ] = false;
                        ( (JFrame) mf[ 0 ] ).dispose();
                        return;
                    }
                    tp.collapse( ab );
                    if ( !ab.isCollapse() ) {
                        System.out.println( "  [CollapsedCladeConfidenceTest] failed to collapse the AB clade" );
                        ok[ 0 ] = false;
                    }

                    // confidence ON: paint and read node coordinates (set during the paint pass)
                    mp.getControlPanel().setCheckbox( DisplayOption.WRITE_CONFIDENCE_VALUES, true );
                    final BufferedImage on = paint( tp );
                    final int ab_on = branchMagenta( on, ab );
                    final int cd_on = branchMagenta( on, cd );
                    if ( ab_on < 5 ) {
                        System.out.println( "  [CollapsedCladeConfidenceTest] no support drawn on the collapsed clade's "
                                + "incoming branch (magenta px = " + ab_on + ")" );
                        ok[ 0 ] = false;
                    }
                    if ( cd_on < 5 ) {
                        System.out.println( "  [CollapsedCladeConfidenceTest] baseline failed: no support drawn on the "
                                + "non-collapsed sibling's branch (magenta px = " + cd_on + ")" );
                        ok[ 0 ] = false;
                    }

                    // confidence OFF: the magenta label must disappear from the collapsed clade's branch, proving
                    // the magenta really is the (toggle-gated) confidence label and not some other artifact
                    mp.getControlPanel().setCheckbox( DisplayOption.WRITE_CONFIDENCE_VALUES, false );
                    final BufferedImage off = paint( tp );
                    final int ab_off = branchMagenta( off, ab );
                    if ( ab_off != 0 ) {
                        System.out.println( "  [CollapsedCladeConfidenceTest] confidence label still drawn on the "
                                + "collapsed branch with the toggle OFF (magenta px = " + ab_off + ")" );
                        ok[ 0 ] = false;
                    }
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = TestFail.here();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Renders the panel offscreen onto a white canvas and returns the image (node coords are set during paint). */
    private static BufferedImage paint( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, W, H );
        tp.paintPhylogeny( g, false, false, 0, 0, 0, 0 );
        g.dispose();
        return img;
    }

    /** Counts magenta-ish (confidence-colored) pixels in the box over {@code node}'s incoming branch. */
    private static int branchMagenta( final BufferedImage img, final PhylogenyNode node ) {
        final int x0 = Math.max( 0, (int) node.getParent().getXcoord() );
        final int x1 = Math.min( img.getWidth() - 1, (int) node.getXcoord() );
        final int y0 = Math.max( 0, (int) ( node.getYcoord() - 22 ) );
        final int y1 = Math.min( img.getHeight() - 1, (int) ( node.getYcoord() + 22 ) );
        int count = 0;
        for( int y = y0; y <= y1; ++y ) {
            for( int x = x0; x <= x1; ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, gg = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( r > 180 ) && ( b > 180 ) && ( gg < 110 ) ) {
                    ++count;
                }
            }
        }
        return count;
    }

    /** First non-root internal node with exactly two descendants whose first child has the given name. */
    private static PhylogenyNode internalWithFirstChild( final Phylogeny phy, final String first_child_name ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.isRoot() && ( n.getNumberOfDescendants() == 2 )
                    && first_child_name.equals( n.getChildNode( 0 ).getName() ) ) {
                return n;
            }
        }
        return null;
    }

    // root -> [ (A,B) with SUPPORT_AB on its branch, (C,D) with SUPPORT_CD on its branch ]
    private static Phylogeny twoClades() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode ab = child( root, null, 0.3 );
        ab.getBranchData().addConfidence( new Confidence( SUPPORT_AB, "bootstrap" ) );
        child( ab, "A", 0.1 );
        child( ab, "B", 0.1 );
        final PhylogenyNode cd = child( root, null, 0.3 );
        cd.getBranchData().addConfidence( new Confidence( SUPPORT_CD, "bootstrap" ) );
        child( cd, "C", 0.1 );
        child( cd, "D", 0.1 );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode child( final PhylogenyNode parent, final String name, final double dist ) {
        final PhylogenyNode n = new PhylogenyNode();
        if ( name != null ) {
            n.setName( name );
        }
        parent.addAsChild( n );
        n.setDistanceToParent( dist );
        return n;
    }

    private CollapsedCladeConfidenceTest() {
    }
}
