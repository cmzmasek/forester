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
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;

import javax.swing.Icon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The "order all" (ladderize) control-panel button is a two-state toggle: each press flips the ladderize
 * direction and shows the just-applied direction on a custom {@code LadderizeIcon}. This drives that toggle
 * on an asymmetric tree and asserts (a) the tip order alternates between the two states, (b) the button icon
 * flips with it, and (c) the two icon states actually paint differently. Needs a display (via
 * {@code MainFrameApplication.createInstance}), so it is a no-op (returns true) when headless; run standalone
 * or with the full suite non-headless.
 */
public final class OrderButtonToggleTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "OrderButtonToggle: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { asymmetricTree() }, new Configuration(), "order" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    ok[ 0 ] = run( mf[ 0 ] );
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
                }
                ( (JFrame) mf[ 0 ] ).dispose();
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    private static boolean run( final MainFrame mf ) {
        final ControlPanel cp = mf.getMainPanel().getControlPanel();
        final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();

        if ( cp.getOrderButtonIconForTest() == null ) {
            return fail( "the order button must carry a (ladderize) icon, not text" );
        }

        cp.orderPressed( tp );
        final String tip1 = firstTip( tp );
        final Icon icon1 = cp.getOrderButtonIconForTest();
        // the mutation appends a ladderize provenance sentence to the tree description
        final String desc1 = tp.getPhylogeny().getDescription();
        if ( ( desc1 == null ) || !desc1.contains( "Ladderized" ) ) {
            return fail( "orderPressed must append a ladderize provenance sentence" );
        }

        cp.orderPressed( tp );
        final String tip2 = firstTip( tp );
        final Icon icon2 = cp.getOrderButtonIconForTest();

        cp.orderPressed( tp );
        final String tip3 = firstTip( tp );
        final Icon icon3 = cp.getOrderButtonIconForTest();

        // the two ladderize states alternate: press 1 and 3 agree, press 2 is the opposite
        if ( tip1.equals( tip2 ) || tip2.equals( tip3 ) || !tip1.equals( tip3 ) ) {
            return fail( "ladderize toggle must alternate the tip order: got " + tip1 + " / " + tip2 + " / " + tip3 );
        }
        // the button icon flips with the state (two distinct icon instances) and returns on the third press
        if ( ( icon1 == icon2 ) || ( icon2 == icon3 ) || ( icon1 != icon3 ) ) {
            return fail( "the order button icon must flip with the ladderize direction" );
        }
        // and the two icon states must actually render differently (mirrored staircase), each with real ink
        final int ink1 = darkPixels( render( icon1 ) );
        final int ink2 = darkPixels( render( icon2 ) );
        if ( ( ink1 < 4 ) || ( ink2 < 4 ) ) {
            return fail( "each ladderize icon must paint a visible staircase (ink1=" + ink1 + ", ink2=" + ink2 + ")" );
        }
        if ( diffPixels( render( icon1 ), render( icon2 ) ) < 4 ) {
            return fail( "the ascending and descending ladderize icons must paint differently" );
        }
        // undo restores the order from before the last press (proves the ladderize is undoable)
        mf.undo();
        final String after_undo = firstTip( tp );
        if ( !after_undo.equals( tip2 ) ) {
            return fail( "undo must restore the pre-ladderize order: expected " + tip2 + " got " + after_undo );
        }
        // and the toggle icon re-syncs to the restored tree's direction (post-press2 = smaller-first here),
        // instead of staying stale on the pre-undo (press3) direction -- matches icon2, and not ascending
        if ( ( cp.getOrderButtonIconForTest() != icon2 ) || cp.isOrderAscendingForTest() ) {
            return fail( "undo must re-sync the ladderize icon to the restored tree's direction" );
        }
        return true;
    }

    private static String firstTip( final TreePanel tp ) {
        return tp.getPhylogeny().getFirstExternalNode().getName();
    }

    private static BufferedImage render( final Icon icon ) {
        final BufferedImage img = new BufferedImage( icon.getIconWidth() + 2, icon.getIconHeight() + 2,
                                                     BufferedImage.TYPE_INT_ARGB );
        final JLabel dummy = new JLabel();
        dummy.setForeground( Color.BLACK );
        icon.paintIcon( dummy, img.getGraphics(), 1, 1 );
        return img;
    }

    private static int darkPixels( final BufferedImage im ) {
        int n = 0;
        for ( int y = 0; y < im.getHeight(); y++ ) {
            for ( int x = 0; x < im.getWidth(); x++ ) {
                final int p = im.getRGB( x, y );
                if ( ( ( ( p >> 24 ) & 0xff ) > 20 ) && ( ( ( p >> 16 ) & 0xff ) < 128 ) ) {
                    n++;
                }
            }
        }
        return n;
    }

    private static int diffPixels( final BufferedImage a, final BufferedImage b ) {
        int n = 0;
        for ( int y = 0; y < a.getHeight(); y++ ) {
            for ( int x = 0; x < a.getWidth(); x++ ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    n++;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [OrderButtonToggleTest] " + msg );
        return false;
    }

    // root -> [ leaf A , clade(B,C,D) ]; ladderizing flips the small leaf vs the larger clade to the top
    private static Phylogeny asymmetricTree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "A" );
        final PhylogenyNode clade = new PhylogenyNode();
        for ( final String s : new String[] { "B", "C", "D" } ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( s );
            clade.addAsChild( n );
        }
        root.addAsChild( a );
        root.addAsChild( clade );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        phy.recalculateNumberOfExternalDescendants( true );
        return phy;
    }

    private OrderButtonToggleTest() {
    }
}
