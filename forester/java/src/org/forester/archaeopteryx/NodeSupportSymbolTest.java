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
import org.forester.archaeopteryx.Options.SUPPORT_VISUALIZATION;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;

/**
 * End-to-end render test for the node-symbol support visualization
 * ({@link TreePanel#paintNodeSupportSymbol}). The per-value math is unit-tested headlessly in
 * {@link TreePanelUtilTest}; this drives a real {@link TreePanel} through the PNG export path and
 * checks that the symbols actually appear on internal nodes and that the modes differ as designed.
 *
 * <p>Needs FlatLaf + a display, so {@link #test()} is a no-op (returns true) when headless -- the
 * same pattern the other GUI tests in this package use; it runs for real from {@link #main} or any
 * non-headless invocation.
 *
 * <p>Tree: two internal nodes under the root carry bootstrap support 99 and 90; a third (60) sits
 * deeper. With the default 0.95 cutoff, THRESHOLD_MARKS marks only the 99 node, SIZE_SCALED draws a
 * dot at all three (sized by support), and NONE draws none -- so dark-pixel counts must satisfy
 * none &lt; threshold &lt; size-scaled.
 */
public final class NodeSupportSymbolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeSupportSymbol: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = supportTree();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "support test" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final JFrame f = (JFrame) mf[ 0 ];
                try {
                    f.setSize( 900, 600 );
                    f.validate();
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    mp.getControlPanel().showWhole();
                    if ( ( tp.getWidth() < 200 ) || ( tp.getHeight() < 200 ) ) {
                        return; // no usable viewport in this environment; nothing to assert
                    }
                    mp.getOptions().setSupportThreshold( 0.95 );
                    final long none = renderDarkPixels( tp, mp, SUPPORT_VISUALIZATION.NONE );
                    final long threshold = renderDarkPixels( tp, mp, SUPPORT_VISUALIZATION.THRESHOLD_MARKS );
                    final long size = renderDarkPixels( tp, mp, SUPPORT_VISUALIZATION.SIZE_SCALED );
                    if ( !( none < threshold ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  threshold marks added no symbol (none=" + none + ", threshold="
                                + threshold + ")" );
                    }
                    if ( !( threshold < size ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  size-scaled should mark more nodes than the 0.95 threshold (threshold="
                                + threshold + ", size=" + size + ")" );
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

    private static long renderDarkPixels( final TreePanel tp, final MainPanel mp, final SUPPORT_VISUALIZATION mode )
            throws Exception {
        mp.getOptions().setSupportVisualization( mode );
        final File out = File.createTempFile( "aptx_support_" + mode, ".png" );
        try {
            AptxUtil.writePhylogenyToGraphicsFile( out.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                   mp.getControlPanel(), GraphicsExportType.PNG, mp.getOptions() );
            final BufferedImage img = ImageIO.read( out );
            long dark = 0;
            for( int y = 0; y < img.getHeight(); y++ ) {
                for( int x = 0; x < img.getWidth(); x++ ) {
                    final int rgb = img.getRGB( x, y );
                    if ( ( ( ( rgb >> 16 ) & 0xFF ) + ( ( rgb >> 8 ) & 0xFF ) + ( rgb & 0xFF ) ) < 200 ) {
                        dark++;
                    }
                }
            }
            return dark;
        }
        finally {
            out.delete();
        }
    }

    /** ((A,(B,C)60)99,(D,E)90)root -- bootstrap support on the internal branches. */
    private static Phylogeny supportTree() {
        final Phylogeny p = new Phylogeny();
        final PhylogenyNode i60 = internal( 60.0, leaf( "B" ), leaf( "C" ) );
        final PhylogenyNode i99 = internal( 99.0, leaf( "A" ), i60 );
        final PhylogenyNode i90 = internal( 90.0, leaf( "D" ), leaf( "E" ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( i99 );
        root.addAsChild( i90 );
        p.setRoot( root );
        p.externalNodesHaveChanged();
        return p;
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode internal( final double conf, final PhylogenyNode a, final PhylogenyNode b ) {
        final PhylogenyNode i = new PhylogenyNode();
        i.getBranchData().addConfidence( new Confidence( conf, "bootstrap" ) );
        i.addAsChild( a );
        i.addAsChild( b );
        return i;
    }

    private NodeSupportSymbolTest() {
    }
}
