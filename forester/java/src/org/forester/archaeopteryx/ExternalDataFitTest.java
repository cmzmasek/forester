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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Guards the "Show External Data" fit bug: when that toggle is OFF, no external-node labels are drawn, so the tip-label
 * reach must NOT be reserved by the fit ("F"/"W") -- otherwise there is a large empty right margin. Asserts the
 * single-source reservation ({@link TreePanel#calculateLongestExtNodeInfo()} -> {@code getLongestExtNodeInfo()})
 * collapses to 0 with the toggle off, and that the depth scale then expands into the freed width. A green no-op when
 * headless.
 */
public final class ExternalDataFitTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ExternalDataFit: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        try {
            final Phylogeny phy = buildPhylogram();
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(),
                                                                        "extfit" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    final ControlPanel cp = tp.getControlPanel();
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );

                    // Show External Data ON: the long tip labels reserve a real right margin
                    cp.setCheckbox( DisplayOption.DISPLAY_EXTERNAL_DATA, true );
                    tp.calcParametersForPainting( 800, 600 );
                    final int reserve_on = tp.getLongestExtNodeInfo();
                    final double corr_on = tp.getXcorrectionFactor();
                    ck( ok, reserve_on > 0, "with external data ON the tip-label reach must be reserved (>0), got "
                            + reserve_on );

                    // Show External Data OFF: labels are not drawn -> reserve NOTHING, branches expand into the space
                    cp.setCheckbox( DisplayOption.DISPLAY_EXTERNAL_DATA, false );
                    tp.calculateLongestExtNodeInfo();
                    ck( ok, tp.getLongestExtNodeInfo() == 0,
                        "with external data OFF the label reservation must collapse to 0, got "
                                + tp.getLongestExtNodeInfo() );
                    tp.calcParametersForPainting( 800, 600 );
                    final double corr_off = tp.getXcorrectionFactor();
                    ck( ok, corr_off > corr_on,
                        "with external data OFF the depth scale must grow (branches use the freed width): off "
                                + corr_off + " vs on " + corr_on );

                    ( (JFrame) mf[ 0 ] ).dispose();
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = TestFail.here();
                }
            } );
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    private static Phylogeny buildPhylogram() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode mid = new PhylogenyNode();
        mid.setDistanceToParent( 0.3 );
        mid.addAsChild( leaf( "a_reasonably_long_tip_label_XXXXXXXXXXXXXXXXXXXX", 0.5 ) );
        mid.addAsChild( leaf( "b_reasonably_long_tip_label_XXXXXXXXXXXXXXXXXXXX", 0.5 ) );
        root.addAsChild( mid );
        root.addAsChild( leaf( "c_reasonably_long_tip_label_XXXXXXXXXXXXXXXXXXXX", 0.9 ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode leaf( final String name, final double bl ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( bl );
        return n;
    }

    private static void ck( final boolean[] ok, final boolean cond, final String msg ) {
        if ( !cond ) {
            System.out.println( "  [ExternalDataFitTest] " + msg );
            ok[ 0 ] = false;
        }
    }

    private ExternalDataFitTest() {
    }
}
