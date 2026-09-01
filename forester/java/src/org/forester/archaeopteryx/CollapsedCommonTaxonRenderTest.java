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
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * A collapsed clade shows its COMPUTED COMMON TAXON (the deepest taxon shared by all its tips) as its label, instead
 * of the boundary tip names. Headful (builds a live {@code MainFrame}); a green no-op when headless. Uses invented
 * taxa carrying an in-tree {@code <lineage>} so the computation is offline + deterministic (never hits a cache).
 */
public final class CollapsedCommonTaxonRenderTest {

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            // ((Aaa, Bbb), Ccc) all share order "Feketiformia" (divergent families) under an UN-annotated clade root
            final PhylogenyNode clade = internal(
                    internal( leaf( "Aaa_one", "Aaa one", "Feketia", "Feketiformia", "Aaidae", "Aaa one" ),
                              leaf( "Bbb_two", "Bbb two", "Feketia", "Feketiformia", "Bbidae", "Bbb two" ) ),
                    leaf( "Ccc_three", "Ccc three", "Feketia", "Feketiformia", "Ccidae", "Ccc three" ) );
            final PhylogenyNode root = internal( clade,
                    leaf( "Zzz_out", "Zzz out", "Feketia", "Zzzformia", "Zzidae", "Zzz out" ) );
            final Phylogeny phy = new Phylogeny();
            phy.setRoot( root );
            phy.externalNodesHaveChanged();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "collapsed common taxon" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    tp.collapse( clade );
                    final String common = tp.collapsedCommonTaxonForTest( clade );
                    if ( !"Feketiformia".equals( common ) ) {
                        fail( ok, "a collapsed clade must show its common taxon (Feketiformia), got \"" + common + "\"" );
                    }
                    // an un-collapsed single node has no "clade" of its own, but the helper still computes over its
                    // tips -- the outgroup tip alone resolves to its own deepest taxon, not shared with the clade
                    final String out = tp.collapsedCommonTaxonForTest( root.getChildNode( 1 ) );
                    if ( "Feketiformia".equals( out ) ) {
                        fail( ok, "the outgroup tip must not share the clade's taxon" );
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
            ok[ 0 ] = TestFail.here();
        }
        return ok[ 0 ];
    }

    private static PhylogenyNode leaf( final String node_name, final String sci, final String... lineage ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( node_name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        t.setLineage( new ArrayList<String>( Arrays.asList( lineage ) ) );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode internal( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode c : children ) {
            n.addAsChild( c );
        }
        return n;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [CollapsedCommonTaxonRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    public static void main( final String[] a ) {
        final boolean ok = test();
        System.out.println( "CollapsedCommonTaxonRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    private CollapsedCommonTaxonRenderTest() {
    }
}
