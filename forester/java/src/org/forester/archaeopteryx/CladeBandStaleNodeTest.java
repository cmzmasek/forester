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
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Deleting nodes must not leave the clade bands pointing at nodes that are no longer in the tree.
 * <p>
 * A {@link CladeBand} holds the node its mark spans. Deleting detaches that node, and the next paint walked it with
 * {@code getAllExternalDescendants()}, which climbs parent links -- off a detached subtree those run out and it
 * threw. In the EDT paint loop that repeats forever, so the window fills with NPEs and stops drawing.
 * <p>
 * Both delete paths are covered: the click-to "Delete Node / Delete Subtree" on the panel, and Tools' prune of the
 * selected external nodes. The subtree / return-to-whole-tree paths already rebuilt the bands; the deletes did not.
 */
public final class CladeBandStaleNodeTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CladeBandStaleNode: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return deleteSubtree() && deleteNodeOnly() && prune();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [CladeBandStaleNodeTest] " + msg );
        return false;
    }

    /** Click-to "Delete Subtree" on a clade that a band spans. */
    private static boolean deleteSubtree() {
        return withBandedTree( ( tp, ok ) -> {
            final int before = tp.cladeBandCount();
            tp.deleteNodeOrSubtreeConfirmed( named( tp, "Chiroptera" ), false );
            if ( !paintsWithoutThrowing( tp ) ) {
                ok[ 0 ] = fail( "painting after deleting a banded subtree threw" );
            }
            // NOT just "it did not crash": the band for the clade that is gone must be gone too, or the figure
            // shows a mark for a taxon that is no longer in the tree
            else if ( tp.cladeBandCount() != ( before - 1 ) ) {
                ok[ 0 ] = fail( "deleting a banded clade should drop its band: " + before + " -> "
                        + tp.cladeBandCount() + " (expected " + ( before - 1 ) + ")" );
            }
        } );
    }

    /** Click-to "Node only" on a banded clade's root -- removeNode splices it out, which also invalidates a band. */
    private static boolean deleteNodeOnly() {
        return withBandedTree( ( tp, ok ) -> {
            tp.deleteNodeOrSubtreeConfirmed( named( tp, "Rodentia" ), true );
            if ( !paintsWithoutThrowing( tp ) ) {
                ok[ 0 ] = fail( "painting after deleting a banded clade's root node threw" );
            }
        } );
    }

    /** The Tools prune path, which deletes several external nodes at once. */
    private static boolean prune() {
        return withBandedTree( ( tp, ok ) -> {
            final int before = tp.cladeBandCount();
            for( final PhylogenyNode tip : tp.getPhylogeny().getExternalNodes()
                    .toArray( new PhylogenyNode[ 0 ] ) ) {
                if ( tip.getName().startsWith( "Chiroptera" ) ) {
                    tp.getPhylogeny().deleteSubtree( tip, true );
                }
            }
            tp.afterTreeStructureChanged(); // what the Tools prune handler runs
            if ( !paintsWithoutThrowing( tp ) ) {
                ok[ 0 ] = fail( "painting after pruning the tips of a banded clade threw" );
            }
            else if ( tp.cladeBandCount() >= before ) {
                ok[ 0 ] = fail( "pruning a banded clade's tips should drop its band: " + before + " -> "
                        + tp.cladeBandCount() );
            }
        } );
    }

    private interface Check {

        void run( TreePanel tp, boolean[] ok );
    }

    private static boolean withBandedTree( final Check check ) {
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { rankedTree() }, new Configuration(), "cladestale" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                if ( tp.setCladeBands( "order", TreePanel.CLADE_VIS.BARS ) < 1 ) {
                    ok[ 0 ] = fail( "precondition: the fixture must place clade bands at rank 'order'" );
                    return;
                }
                if ( !paintsWithoutThrowing( tp ) ) {
                    ok[ 0 ] = fail( "precondition: the banded tree must paint before anything is deleted" );
                    return;
                }
                check.run( tp, ok );
            } );
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Paints offscreen, reporting whether the paint threw (the real symptom -- an NPE in the EDT paint loop). */
    private static boolean paintsWithoutThrowing( final TreePanel tp ) {
        try {
            final int w = 700, h = 500;
            tp.setSize( w, h );
            tp.validate();
            tp.doLayout();
            final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
            final Graphics2D g = img.createGraphics();
            g.setColor( Color.WHITE );
            g.fillRect( 0, 0, w, h );
            tp.paintPhylogeny( g, false, false, w, h, 0, 0 );
            g.dispose();
            return true;
        }
        catch ( final Throwable t ) {
            System.out.println( "  [CladeBandStaleNodeTest] paint threw: " + t );
            return false;
        }
    }

    private static PhylogenyNode named( final TreePanel tp, final String name ) {
        for( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = tp.getPhylogeny()
                .iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( ( n.getNodeData().isHasTaxonomy() )
                    && name.equals( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                return n;
            }
        }
        throw new IllegalStateException( "no clade " + name );
    }

    /** Orders inside a class, so 'order' places real clade bands with several tips each. */
    private static Phylogeny rankedTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode mammalia = new PhylogenyNode();
        mammalia.getNodeData().setTaxonomy( ranked( "Mammalia", "class" ) );
        for( final String order : new String[] { "Chiroptera", "Rodentia", "Primates" } ) {
            final PhylogenyNode on = new PhylogenyNode();
            on.getNodeData().setTaxonomy( ranked( order, "order" ) );
            for( int i = 0; i < 3; ++i ) {
                final PhylogenyNode leaf = new PhylogenyNode();
                leaf.setName( order + "_" + i );
                on.addAsChild( leaf );
            }
            mammalia.addAsChild( on );
        }
        root.addAsChild( mammalia );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static Taxonomy ranked( final String name, final String rank ) {
        final Taxonomy t = new Taxonomy();
        t.setScientificName( name );
        try {
            t.setRank( rank );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        return t;
    }

    private CladeBandStaleNodeTest() {
    }
}
