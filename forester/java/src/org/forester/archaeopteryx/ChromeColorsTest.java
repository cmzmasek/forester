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

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Interactive CHROME must not wear a DATA colour.
 * <p>
 * Red means "this matched your search" and nothing else. The overview's viewport box used to be drawn in that very
 * colour -- inside the overview, where the found-node marks legitimately use it -- so a navigation affordance and a
 * search hit were indistinguishable. A subtree held by Cut/Copy borrowed it too. Both now have colours of their own,
 * and this pins that, so the shortcut of "just reuse the found colour" cannot come back.
 */
public final class ChromeColorsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ChromeColors: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return distinctFromData() && overviewBoxNotFoundColour();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ChromeColorsTest] " + msg );
        return false;
    }

    /** The chrome colours differ from every data colour they could be confused with. */
    private static boolean distinctFromData() {
        final TreeColorSet tcs = TreeColorSet.createInstance();
        for( final int scheme : new int[] { TreeColorSet.DARK_COLOR_SCHEME, TreeColorSet.LIGHT_COLOR_SCHEME } ) {
            tcs.setColorSchema( scheme );
            final Color accent = TreePanel.uiAccentColor();
            final Color clipboard = TreePanel.clipboardHighlightColorForTest();
            for( final Color found : new Color[] { tcs.getFoundColor0(), tcs.getFoundColor1(),
                                                   tcs.getFoundColor0and1() } ) {
                if ( accent.equals( found ) ) {
                    return fail( "the navigation accent must not be a search-hit colour: " + accent );
                }
                if ( clipboard.equals( found ) ) {
                    return fail( "the Cut/Copy highlight must not be a search-hit colour: " + clipboard );
                }
            }
            // ...nor may the held-subtree amber be mistaken for the duplication-or-speciation event box
            if ( clipboard.equals( tcs.getDuplicationOrSpeciationColor() ) ) {
                return fail( "the Cut/Copy highlight must differ from the event amber" );
            }
        }
        // and the error marker is a warning tone, not raw red
        if ( Color.RED.equals( MainFrame.errorMarkerColorForTest() ) ) {
            return fail( "the error marker must not be raw red" );
        }
        return true;
    }

    /** The drawn proof: with the overview showing a viewport box, no search-hit red is painted in it. */
    private static boolean overviewBoxNotFoundColour() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final boolean[] ok = { true };
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { tree() }, new Configuration(), "chrome" ) );
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final int w = 600;
                final int h = 460;
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.getOptions().setShowOverview( true );
                tp.getOptions().setGraphicsExportWhiteBackground( true );
                mf[ 0 ].getMainPanel().getControlPanel().showWhole();
                // the tree must be LARGER than the viewport, or there is no viewport box to draw
                tp.setPreferredSize( new java.awt.Dimension( w * 3, h * 3 ) );
                tp.setSize( w, h );
                tp.validate();
                tp.doLayout();
                tp.calcParametersForPainting( w, h );
                tp.setOvOn( true );
                final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
                final Graphics2D g = img.createGraphics();
                g.setColor( Color.WHITE );
                g.fillRect( 0, 0, w, h );
                tp.printAll( g );
                g.dispose();
                if ( !tp.isOvOn() ) {
                    ok[ 0 ] = fail( "precondition: the overview must be on, or this proves nothing" );
                    return;
                }
                // the viewport box is the only thing drawn in the accent here, so it must be present...
                final Color accent = TreePanel.uiAccentColor();
                int accent_px = 0;
                int found_px = 0;
                final Color found = tp.getTreeColorSet().getFoundColor0();
                for( int x = 0; x < w; ++x ) {
                    for( int y = 0; y < h; ++y ) {
                        final int rgb = img.getRGB( x, y );
                        if ( near( rgb, accent ) ) {
                            accent_px++;
                        }
                        if ( near( rgb, found ) ) {
                            found_px++;
                        }
                    }
                }
                if ( accent_px < 10 ) {
                    ok[ 0 ] = fail( "the overview's viewport box must be drawn in the accent, found " + accent_px
                            + " such pixels" );
                }
                // ...and NOTHING may be drawn in the search-hit colour: nothing is selected in this tree
                if ( found_px > 0 ) {
                    ok[ 0 ] = fail( "no search-hit colour may appear when nothing is found: " + found_px + " pixels" );
                }
            } );
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean near( final int rgb, final Color c ) {
        return ( Math.abs( ( ( rgb >> 16 ) & 0xff ) - c.getRed() ) < 12 )
                && ( Math.abs( ( ( rgb >> 8 ) & 0xff ) - c.getGreen() ) < 12 )
                && ( Math.abs( ( rgb & 0xff ) - c.getBlue() ) < 12 );
    }

    private static Phylogeny tree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 8; ++i ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "tip_" + i );
            n.setDistanceToParent( 0.1 );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private ChromeColorsTest() {
    }
}
