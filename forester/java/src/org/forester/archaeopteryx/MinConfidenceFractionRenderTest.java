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
 * Integration test with teeth for the scale-aware "Min. confidence shown" default: the setting is a FRACTION of
 * the tree's detected support scale (default 0.5), so on a 0-100 bootstrap tree the effective cutoff is 50. This
 * renders a tree with one clade above the cutoff (95) and one below it (40) at the default setting and asserts
 * the first is drawn and the second is NOT -- which fails if {@code paintConfidenceValues} ever regresses to
 * using the fraction as an absolute threshold (0.5, under which both 95 and 40 would show). Uses the same
 * unique-confidence-color pixel probe as {@link CollapsedCladeConfidenceTest}. Headless-guarded (needs FlatLaf
 * via {@code createInstance}); run standalone or in the non-headless suite.
 */
public final class MinConfidenceFractionRenderTest {

    private static final int W = 800;
    private static final int H = 600;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MinConfidenceFractionRender: " + ( ok ? "OK." : "FAILED." ) );
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
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { bootstrapTree() }, conf,
                                                                        "frac" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    if ( mp.getControlPanel().getWriteConfidenceCb() == null ) {
                        System.out.println( "  [MinConfidenceFractionRenderTest] no confidence control -- cannot test" );
                        ok[ 0 ] = false;
                        ( (JFrame) mf[ 0 ] ).dispose();
                        return;
                    }
                    // default "min. confidence shown" fraction (0.5), NOT overridden -> cutoff = 0.5 * 100 = 50
                    if ( mf[ 0 ].getOptions().getMinConfidenceFraction() != 0.5 ) {
                        System.out.println( "  [MinConfidenceFractionRenderTest] expected default fraction 0.5" );
                        ok[ 0 ] = false;
                    }
                    mp.getControlPanel().getWriteConfidenceCb().setSelected( true );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, Color.WHITE );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BRANCH, Color.BLACK );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.CONFIDENCE, new Color( 255, 0, 255 ) );
                    tp.getTreeColorSet().setColorSchema( 0 );

                    final PhylogenyNode high = internalWithFirstChild( tp.getPhylogeny(), "A" ); // support 95 -> shown
                    final PhylogenyNode low = internalWithFirstChild( tp.getPhylogeny(), "C" );  // support 40 -> hidden
                    tp.setSize( W, H );
                    tp.calcParametersForPainting( W, H );
                    final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
                    final Graphics2D g = img.createGraphics();
                    g.setColor( Color.WHITE );
                    g.fillRect( 0, 0, W, H );
                    tp.paintPhylogeny( g, false, false, 0, 0, 0, 0 );
                    g.dispose();

                    if ( branchMagenta( img, high ) < 5 ) {
                        System.out.println( "  [MinConfidenceFractionRenderTest] support 95 (>= cutoff 50) should be shown" );
                        ok[ 0 ] = false;
                    }
                    // The teeth: 40 is above the fraction (0.5) but below the scale-aware cutoff (50), so a correct
                    // fraction*scale filter hides it; a fraction-as-absolute regression would wrongly show it.
                    if ( branchMagenta( img, low ) != 0 ) {
                        System.out.println( "  [MinConfidenceFractionRenderTest] support 40 (< cutoff 50) should be "
                                + "HIDDEN by the fraction-of-scale default (magenta px = " + branchMagenta( img, low ) + ")" );
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

    // root -> [ (A,B) bootstrap 95, (C,D) bootstrap 40 ] -- scale 0-100 (ceiling 100)
    private static Phylogeny bootstrapTree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode ab = child( root, null, 0.3 );
        ab.getBranchData().addConfidence( new Confidence( 95.0, "bootstrap" ) );
        child( ab, "A", 0.1 );
        child( ab, "B", 0.1 );
        final PhylogenyNode cd = child( root, null, 0.3 );
        cd.getBranchData().addConfidence( new Confidence( 40.0, "bootstrap" ) );
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

    private MinConfidenceFractionRenderTest() {
    }
}
