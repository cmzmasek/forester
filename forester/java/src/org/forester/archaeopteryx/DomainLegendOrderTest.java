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
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;

/**
 * The protein-domain legend lists names in first-appearance order with the tips taken in DISPLAY order, top to bottom
 * -- so with "Reverse Tip Order" on it follows what is on screen, not the stored child order (Christian, 2026-09-12;
 * the rule Archaeopteryx.js uses). Before this, the reversed display kept the forward legend and read bottom to top.
 */
public final class DomainLegendOrderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DomainLegendOrder: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // needs a display
        }
        try {
            return legendFollowsReverseTipOrder();
        }
        catch ( final Throwable e ) {
            e.printStackTrace( System.out );
            return false;
        }
    }

    private static boolean legendFollowsReverseTipOrder() throws Exception {
        // three tips, top to bottom: t0 [Alpha, Shared], t1 [Beta], t2 [Gamma, Shared]. The forward and the reversed
        // first-appearance orders differ, so the test can see which one the legend used.
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( tip( "t0", "Alpha", "Shared" ) );
        root.addAsChild( tip( "t1", "Beta", null ) );
        root.addAsChild( tip( "t2", "Gamma", "Shared" ) );
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        final List<String> forward = Arrays.asList( "Alpha (1)", "Shared (2)", "Beta (1)", "Gamma (1)" );
        final List<String> reversed = Arrays.asList( "Gamma (1)", "Shared (2)", "Beta (1)", "Alpha (1)" );
        final MainFrame[] mf = new MainFrame[ 1 ];
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "domainlegend" ) );
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                mf[ 0 ].getOptions().setReverseTipOrder( false );
                final List<String> off = tp.domainLegendRowsForTest();
                if ( !forward.equals( off ) ) {
                    ok[ 0 ] = fail( "tips drawn in stored order: the legend should read " + forward + ", got " + off );
                }
                mf[ 0 ].getOptions().setReverseTipOrder( true );
                final List<String> on = tp.domainLegendRowsForTest();
                if ( !reversed.equals( on ) ) {
                    ok[ 0 ] = fail( "with Reverse Tip Order on, the legend must follow the display " + reversed + ", got "
                            + on );
                }
                mf[ 0 ].getOptions().setReverseTipOrder( false );
                if ( !forward.equals( tp.domainLegendRowsForTest() ) ) {
                    ok[ 0 ] = fail( "turning Reverse Tip Order off again must restore the forward legend" );
                }
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> {
                if ( mf[ 0 ] != null ) {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
        }
        return ok[ 0 ];
    }

    /** A tip carrying one sequence whose architecture holds {@code first} (and {@code second}, when not null). */
    private static PhylogenyNode tip( final String name, final String first, final String second ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final DomainArchitecture da = new DomainArchitecture();
        da.setTotalLength( 100 );
        da.addDomain( new ProteinDomain( first, 1, 40, 1e-10 ) );
        if ( second != null ) {
            da.addDomain( new ProteinDomain( second, 50, 90, 1e-10 ) );
        }
        final Sequence s = new Sequence();
        s.setDomainArchitecture( da );
        n.getNodeData().setSequence( s );
        return n;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [DomainLegendOrderTest] " + message );
        return false;
    }

    private DomainLegendOrderTest() {
        // not instantiable
    }
}
