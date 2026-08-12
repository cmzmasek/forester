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
import java.awt.Rectangle;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * The draggable "Internal Taxonomy Key" draws (records its bounds) when the option is on AND the tree carries internal
 * taxa, and draws NOTHING (bounds nulled) for a tree with no internal taxa. Headful; a green no-op when headless.
 */
public final class InternalTaxaKeyRenderTest {

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            // internal taxa at order + family; the tips carry none
            final Phylogeny withTaxa = phylogeny( node(
                    tax( "order", "Carnivora", tax( "family", "Felidae", tip( "a" ), tip( "b" ) ), tip( "c" ) ),
                    tax( "order", "Rodentia", tip( "d" ), tip( "e" ) ) ) );
            final Phylogeny plain = phylogeny( node( node( tip( "p" ), tip( "q" ) ), tip( "r" ) ) );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { withTaxa, plain }, new Configuration(), "internal taxa key" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    frame.getOptions().setShowInternalTaxonomyKey( true );
                    // tab 0 carries internal taxa -> the key draws (records bounds)
                    frame.getMainPanel().getTabbedPane().setSelectedIndex( 0 );
                    final TreePanel tp0 = frame.getMainPanel().getCurrentTreePanel();
                    paint( tp0 );
                    final Rectangle b = tp0.getInternalTaxaKeyBoundsForTest();
                    if ( ( b == null ) || (b.width <= 0) || (b.height <= 0) ) {
                        fail( ok, "the internal-taxonomy key must draw (record bounds) for a tree with internal taxa" );
                    }
                    // tab 1 has NO internal taxa -> the key draws nothing (bounds nulled) even with the option on
                    frame.getMainPanel().getTabbedPane().setSelectedIndex( 1 );
                    final TreePanel tp1 = frame.getMainPanel().getCurrentTreePanel();
                    paint( tp1 );
                    if ( tp1.getInternalTaxaKeyBoundsForTest() != null ) {
                        fail( ok, "the internal-taxonomy key must be empty (no bounds) for a tree with no internal taxa" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
        return ok[ 0 ];
    }

    private static void paint( final TreePanel tp ) {
        tp.setSize( 760, 460 );
        tp.calcParametersForPainting( 760, 460 );
        final BufferedImage img = new BufferedImage( 760, 460, BufferedImage.TYPE_INT_ARGB );
        tp.paintPhylogeny( img.createGraphics(), false, false, 0, 0, 0, 0 );
    }

    private static PhylogenyNode tip( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode tax( final String rank, final String sci, final PhylogenyNode... kids ) {
        final PhylogenyNode n = node( kids );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        try {
            t.setRank( rank );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode node( final PhylogenyNode... kids ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode k : kids ) {
            n.addAsChild( k );
        }
        return n;
    }

    private static Phylogeny phylogeny( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [InternalTaxaKeyRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    public static void main( final String[] a ) {
        final boolean ok = test();
        System.out.println( "InternalTaxaKeyRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    private InternalTaxaKeyRenderTest() {
    }
}
