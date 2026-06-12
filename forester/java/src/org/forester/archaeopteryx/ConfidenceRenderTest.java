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
import java.awt.image.BufferedImage;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;

/**
 * Integration test that a non-finite (NaN/Infinity) branch confidence value does not crash the tree
 * rendering: with "Confidence Values" turned on, painting a tree that carries a NaN confidence used
 * to throw {@code NumberFormatException} from {@code ForesterUtil.round}. Guarded to a no-op on a
 * headless box. Needs FlatLaf on the classpath (via {@code MainFrameApplication.createInstance}), so
 * it is run standalone, not as part of the headless suite.
 */
public final class ConfidenceRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ConfidenceRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration( null, false, false, true );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { naNConfidenceTree() }, conf,
                                                                        "conf" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    if ( mp.getControlPanel().getWriteConfidenceCb() == null ) {
                        ok[ 0 ] = false; // the test can only exercise the render path with this control present
                    }
                    else {
                        mp.getControlPanel().getWriteConfidenceCb().setSelected( true );
                        tp.setSize( 800, 600 );
                        final BufferedImage img = new BufferedImage( 800, 600, BufferedImage.TYPE_INT_ARGB );
                        final Graphics2D g = img.createGraphics();
                        tp.printAll( g ); // threw NumberFormatException("Infinite or NaN") before the fix
                        g.dispose();
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
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

    // root -> [ internal(A,B) with a NaN confidence on its branch, leaf C ]
    private static Phylogeny naNConfidenceTree() {
        final PhylogenyNode internal = new PhylogenyNode();
        internal.addAsChild( leaf( "A" ) );
        internal.addAsChild( leaf( "B" ) );
        internal.setDistanceToParent( 1.0 );
        internal.getBranchData().addConfidence( new Confidence( Double.NaN, "MAD" ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( internal );
        root.addAsChild( leaf( "C" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( 1.0 );
        return n;
    }

    private ConfidenceRenderTest() {
    }
}
