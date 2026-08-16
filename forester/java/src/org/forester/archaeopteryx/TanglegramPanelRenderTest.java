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
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;

import org.forester.archaeopteryx.TanglegramLinker.LinkField;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headful render test for {@link TanglegramPanel}: renders a reversed 4-tip pair offscreen and asserts that both
 * trees + the tip-to-tip connectors are actually drawn (ink in the central connector band, which is clear of the
 * two trees and their labels), and that the link/crossing counts are correct. A green no-op when headless.
 */
public final class TanglegramPanelRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TanglegramPanelRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        // constructing a panel does not paint, so the count / label checks are headless-safe; only renderOk paints
        if ( !countsOk() || !labelFallbackOk() ) {
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return renderOk();
    }

    private static boolean labelFallbackOk() {
        // tips carry ONLY a taxonomy code (no node name, no scientific name) and are linked by that code -- the
        // drawn label must still resolve to the code, not a blank row (regression for the label-fallback fix)
        try {
            final TanglegramPanel panel = new TanglegramPanel( tree( clade( codeLeaf( "FELCA" ), codeLeaf( "CANLU" ) ) ),
                                                               tree( clade( codeLeaf( "CANLU" ), codeLeaf( "FELCA" ) ) ),
                                                               LinkField.TAXONOMY_CODE );
            if ( panel.getResult().getLinks().size() != 2 ) {
                return fail( "code-linked tips should form 2 links, got " + panel.getResult().getLinks().size() );
            }
            if ( !"FELCA".equals( panel.leftTipLabelForTest( 0 ) )
                    || !"CANLU".equals( panel.leftTipLabelForTest( 1 ) ) ) {
                return fail( "code-only tips must label with their taxonomy code, got '"
                        + panel.leftTipLabelForTest( 0 ) + "'/'" + panel.leftTipLabelForTest( 1 ) + "'" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "code-leaf setup failed: " + e.getMessage() );
        }
    }

    private static boolean countsOk() {
        final TanglegramPanel panel = new TanglegramPanel( treeABCD(), treeDCBA(), LinkField.NODE_NAME );
        if ( panel.getResult().getLinks().size() != 4 ) {
            return fail( "expected 4 connectors, got " + panel.getResult().getLinks().size() );
        }
        if ( panel.getUnmatchedCount() != 0 ) {
            return fail( "expected 0 unmatched, got " + panel.getUnmatchedCount() );
        }
        if ( panel.getCrossingCount() != 6 ) {
            return fail( "a fully reversed 4-tip tanglegram should have 6 crossings, got " + panel.getCrossingCount() );
        }
        return true;
    }

    private static boolean renderOk() {
        final TanglegramPanel panel = new TanglegramPanel( treeABCD(), treeDCBA(), LinkField.NODE_NAME );
        final int w = 800;
        final int h = 400;
        panel.setBackground( Color.WHITE );
        panel.setForeground( Color.BLACK );
        panel.setFont( new Font( "SansSerif", Font.PLAIN, 12 ) );
        panel.setSize( w, h );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            panel.printAll( g );
        }
        finally {
            g.dispose();
        }
        final int total_ink = countNonWhite( img, 0, w );
        if ( total_ink < 200 ) {
            return fail( "the tanglegram rendered almost no ink (" + total_ink + " px) -- trees not drawn?" );
        }
        // a narrow central band lies BETWEEN the two trees + their labels, so only connectors cross it
        final int band_ink = countNonWhite( img, ( w / 2 ) - 5, ( w / 2 ) + 5 );
        if ( band_ink <= 0 ) {
            return fail( "no ink in the central connector band -- connectors not drawn" );
        }
        return true;
    }

    private static int countNonWhite( final BufferedImage img, final int x_from, final int x_to ) {
        int count = 0;
        for( int x = Math.max( 0, x_from ); x < Math.min( img.getWidth(), x_to ); x++ ) {
            for( int y = 0; y < img.getHeight(); y++ ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) != 0xFFFFFF ) {
                    ++count;
                }
            }
        }
        return count;
    }

    private static boolean fail( final String message ) {
        System.out.println( "TanglegramPanelRender test failed: " + message );
        return false;
    }

    private static Phylogeny treeABCD() {
        return tree( clade( leaf( "A" ), leaf( "B" ), leaf( "C" ), leaf( "D" ) ) );
    }

    private static Phylogeny treeDCBA() {
        return tree( clade( leaf( "D" ), leaf( "C" ), leaf( "B" ), leaf( "A" ) ) );
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    /** A tip with ONLY a taxonomy code (no node name, no scientific name). */
    private static PhylogenyNode codeLeaf( final String code ) throws Exception {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setTaxonomyCode( code );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode clade( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode child : children ) {
            n.addAsChild( child );
        }
        return n;
    }

    private static Phylogeny tree( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }
}
