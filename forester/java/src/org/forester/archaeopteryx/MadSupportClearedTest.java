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

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;

/**
 * Integration test that the per-branch MAD support added by {@link TreePanel#madRoot()} is removed
 * when the tree is rerooted by another means: {@link TreePanel#midpointRoot()} (Tools menu) and
 * {@link TreePanel#reRoot(PhylogenyNode)} (clicking a non-root node). Guarded to a no-op on a
 * headless box. Needs FlatLaf on the classpath (via {@code MainFrameApplication.createInstance}), so
 * it is run standalone, not as part of the headless suite.
 */
public final class MadSupportClearedTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MadSupportCleared: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny phy = factory.create( "(A:1,B:1,(C:1,D:3):1)", new NHXParser() )[ 0 ];
            final Configuration conf = new Configuration( null, false, false, true );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "mad" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                tp.madRoot();
                if ( !hasMad( tp.getPhylogeny() ) ) {
                    ok[ 0 ] = false; // MAD rooting adds the support
                }
                tp.midpointRoot();
                if ( hasMad( tp.getPhylogeny() ) ) {
                    ok[ 0 ] = false; // midpoint rooting clears it
                }
                tp.madRoot();
                if ( !hasMad( tp.getPhylogeny() ) ) {
                    ok[ 0 ] = false;
                }
                final PhylogenyNode internal = firstInternalNonRoot( tp.getPhylogeny() );
                if ( internal == null ) {
                    ok[ 0 ] = false;
                }
                else {
                    tp.reRoot( internal ); // a manual reroot on a non-root node clears it
                    if ( hasMad( tp.getPhylogeny() ) ) {
                        ok[ 0 ] = false;
                    }
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

    private static boolean hasMad( final Phylogeny phy ) {
        for ( final PhylogenyNode nd : PhylogenyMethods.obtainAllNodesAsList( phy ) ) {
            for ( final Confidence c : nd.getBranchData().getConfidences() ) {
                if ( PhylogenyMethods.MAD_CONFIDENCE_TYPE.equals( c.getType() ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    private static PhylogenyNode firstInternalNonRoot( final Phylogeny phy ) {
        for ( final PhylogenyNode nd : PhylogenyMethods.obtainAllNodesAsList( phy ) ) {
            if ( nd.isInternal() && !nd.isRoot() ) {
                return nd;
            }
        }
        return null;
    }

    private MadSupportClearedTest() {
    }
}
