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
import java.math.BigDecimal;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.data.Date;

/**
 * Headful test for the "Time tree" badge on a real {@link TreePanel}: a DATED tree is auto-labeled (with its unit),
 * a plain phylogram is not, and an ULTRAMETRIC tree is labeled only after {@code setConfirmedTimeTree}. A no-op on a
 * headless box (needs FlatLaf via {@code createInstance}).
 */
public final class TimeTreeLabelTest {

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final boolean[] ok = { true };
            // DATED -> auto-labelled with the unit
            withPanel( datedTree(), tp -> {
                final String label = tp.timeTreeBadgeLabel();
                ck( ok, ( label != null ) && label.contains( "Time tree" ) && label.contains( "mya" ),
                    "a dated tree is auto-labelled a time tree (with unit): " + label );
            } );
            // plain ragged phylogram -> no label
            withPanel( raggedTree(), tp -> ck( ok, tp.timeTreeBadgeLabel() == null,
                                               "a ragged phylogram is NOT labelled a time tree" ) );
            // ULTRAMETRIC -> no label until the user confirms
            withPanel( ultrametricTree(), tp -> {
                ck( ok, tp.timeTreeBadgeLabel() == null, "an ultrametric tree is NOT auto-labelled (it is offered)" );
                tp.setConfirmedTimeTree( true );
                ck( ok, "Time tree".equals( tp.timeTreeBadgeLabel() ) && tp.isConfirmedTimeTreeForTest(),
                    "a confirmed ultrametric tree IS labelled a time tree" );
                // the confirmation is bound to the tree object: replacing the tree drops the stale badge
                tp.setTree( raggedTree() );
                ck( ok, tp.timeTreeBadgeLabel() == null,
                    "the time-tree confirmation lapses when the panel's tree is replaced" );
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private interface PanelCheck {

        void run( TreePanel tp );
    }

    private static void withPanel( final Phylogeny phy, final PanelCheck check ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(),
                                                                     "time tree test" ) );
        SwingUtilities.invokeAndWait( () -> {
            try {
                check.run( mf[ 0 ].getMainPanel().getCurrentTreePanel() );
            }
            finally {
                ( (JFrame) mf[ 0 ] ).dispose();
            }
        } );
    }

    private static Phylogeny datedTree() {
        final Phylogeny phy = AptxUtilTest.balancedFourTip( 0, 0, 0 );
        phy.getRoot().getNodeData().setDate( new Date( "", new BigDecimal( "90" ), null, null, "mya" ) );
        phy.getRoot().getChildNode( 0 ).getNodeData()
                .setDate( new Date( "", new BigDecimal( "30" ), null, null, "mya" ) );
        return phy;
    }

    private static Phylogeny ultrametricTree() {
        return AptxUtilTest.balancedFourTip( 1, 1, 1 );
    }

    private static Phylogeny raggedTree() {
        final Phylogeny phy = AptxUtilTest.balancedFourTip( 1, 1, 1 );
        phy.getFirstExternalNode().setDistanceToParent( 5 );
        return phy;
    }

    private static void ck( final boolean[] ok, final boolean cond, final String msg ) {
        if ( !cond ) {
            System.out.println( "  [TimeTreeLabelTest] " + msg );
            ok[ 0 ] = false;
        }
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private TimeTreeLabelTest() {
    }
}
