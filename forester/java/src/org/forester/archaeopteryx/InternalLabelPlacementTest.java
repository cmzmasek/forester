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

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * End-to-end render test for the publication-style internal-label placement (the option
 * {@link Options#isInternalLabelsAboveBranch()}, drawn by {@link TreePanel#paintInternalLabelAboveBranch}).
 * The placement arithmetic is unit-tested headlessly in {@link TreePanelUtilTest}; this drives a real
 * {@link TreePanel} through the PNG export path and checks that toggling the option actually moves an
 * internal node's label from the RIGHT of the node to ABOVE-and-LEFT of it.
 *
 * <p>Needs FlatLaf + a display, so {@link #test()} is a no-op (returns true) when headless -- the same
 * pattern the other GUI tests in this package use; it runs for real from {@link #main} or any
 * non-headless invocation.
 *
 * <p>Tree: {@code ((A,B)CladeAlpha,(C,D)CladeBeta)root} with internal-node names shown. The named
 * internal node CladeAlpha carries a label; with the option ON it must sit above-and-left of the node,
 * with it OFF it must sit to the right -- so the dark-pixel counts in the two regions must swap.
 */
public final class InternalLabelPlacementTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "InternalLabelPlacement: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = namedInternalTree();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, conf, "internal label test" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final JFrame f = (JFrame) mf[ 0 ];
                try {
                    f.setSize( 900, 600 );
                    f.validate();
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    final ControlPanel cp = mp.getControlPanel();
                    cp.setCheckbox( Configuration.display_internal_data, true );
                    cp.setCheckbox( Configuration.show_node_names, true );
                    cp.showWhole();
                    if ( ( tp.getWidth() < 200 ) || ( tp.getHeight() < 200 ) || !cp.isShowNodeNames() ) {
                        return; // no usable viewport or node-name display in this environment; nothing to assert
                    }
                    final PhylogenyNode clade = phy.getNode( "A" ).getParent(); // the named internal node "CladeAlpha"
                    // ON: label above-and-left of the node. Read coords from this render (layout is identical
                    // either way, since the label grows leftward over the branch and changes no geometry).
                    final BufferedImage on = render( tp, mp, true );
                    final int nx = Math.round( clade.getXcoord() );
                    final int ny = Math.round( clade.getYcoord() );
                    final BufferedImage off = render( tp, mp, false );
                    // windows clear of the branch line (which is at y = ny in both images): the label only
                    final long aboveleft_on = darkInWindow( on, nx - 28, ny - 9, 9 );
                    final long aboveleft_off = darkInWindow( off, nx - 28, ny - 9, 9 );
                    final long right_on = darkInWindow( on, nx + 32, ny - 3, 9 );
                    final long right_off = darkInWindow( off, nx + 32, ny - 3, 9 );
                    if ( !( aboveleft_on > ( aboveleft_off + 15 ) ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  label did not move above-and-left when option ON (aboveleft_on="
                                + aboveleft_on + ", aboveleft_off=" + aboveleft_off + ")" );
                    }
                    if ( !( right_off > ( right_on + 15 ) ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  label did not sit to the right when option OFF (right_off=" + right_off
                                + ", right_on=" + right_on + ")" );
                    }
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                    ok[ 0 ] = false;
                }
                finally {
                    f.dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    private static BufferedImage render( final TreePanel tp, final MainPanel mp, final boolean above_branch )
            throws Exception {
        mp.getOptions().setInternalLabelsAboveBranch( above_branch );
        final File out = File.createTempFile( "aptx_intlabel_" + above_branch, ".png" );
        try {
            AptxUtil.writePhylogenyToGraphicsFile( out.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                   mp.getControlPanel(), GraphicsExportType.PNG, mp.getOptions() );
            return ImageIO.read( out );
        }
        finally {
            out.delete();
        }
    }

    private static boolean isDark( final BufferedImage img, final int x, final int y ) {
        if ( ( x < 0 ) || ( y < 0 ) || ( x >= img.getWidth() ) || ( y >= img.getHeight() ) ) {
            return false;
        }
        final int rgb = img.getRGB( x, y );
        return ( ( ( ( rgb >> 16 ) & 0xFF ) + ( ( rgb >> 8 ) & 0xFF ) + ( rgb & 0xFF ) ) < 200 );
    }

    /** Dark pixels in the (2r+1)x(2r+1) window centered at (cx, cy). */
    private static long darkInWindow( final BufferedImage img, final int cx, final int cy, final int r ) {
        long n = 0;
        for( int y = cy - r; y <= ( cy + r ); ++y ) {
            for( int x = cx - r; x <= ( cx + r ); ++x ) {
                if ( isDark( img, x, y ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** {@code ((A,B)CladeAlpha,(C,D)CladeBeta)root} -- names on the internal nodes. */
    private static Phylogeny namedInternalTree() {
        final Phylogeny p = new Phylogeny();
        final PhylogenyNode alpha = internal( "CladeAlpha", leaf( "A" ), leaf( "B" ) );
        final PhylogenyNode beta = internal( "CladeBeta", leaf( "C" ), leaf( "D" ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( alpha );
        root.addAsChild( beta );
        p.setRoot( root );
        p.externalNodesHaveChanged();
        return p;
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode internal( final String name, final PhylogenyNode a, final PhylogenyNode b ) {
        final PhylogenyNode i = new PhylogenyNode();
        i.setName( name );
        i.addAsChild( a );
        i.addAsChild( b );
        return i;
    }

    private InternalLabelPlacementTest() {
    }
}
