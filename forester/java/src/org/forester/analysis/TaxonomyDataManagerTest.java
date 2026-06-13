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

package org.forester.analysis;

import java.util.SortedSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Test for the offline behavior of {@link TaxonomyDataManager#obtainDetailedTaxonomicInformation}.
 *
 * <p>With {@code allow_to_use_basic_node_names == false} and nodes that carry no {@link
 * org.forester.phylogeny.data.Taxonomy} object, the method never contacts a web service: every
 * external node is simply reported as not-found and internal nodes are ignored. That lets us
 * deterministically test the iteration / collection logic (and the absence of the old delete
 * path) without any network. The actual UniProt resolution path needs live services and is not
 * unit-tested here.
 */
public final class TaxonomyDataManagerTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TaxonomyDataManager: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            // Only external nodes are reported; named internal nodes are NOT.
            final SortedSet<String> nf = TaxonomyDataManager.obtainDetailedTaxonomicInformation( nestedTree(), false );
            if ( nf.size() != 3 ) {
                return fail( "expected 3 not-found external nodes, got " + nf.size() + ": " + nf );
            }
            if ( !nf.contains( "A" ) || !nf.contains( "B" ) || !nf.contains( "C" ) ) {
                return fail( "not-found set is missing an external leaf name: " + nf );
            }
            if ( nf.contains( "ROOT" ) || nf.contains( "INNER" ) ) {
                return fail( "internal nodes must not be reported as not-found: " + nf );
            }

            // A tree whose leaves all resolve trivially would contact the network; instead use an
            // unnamed leaf, which is reported (by node.toString()) without any lookup.
            final SortedSet<String> nf2 = TaxonomyDataManager.obtainDetailedTaxonomicInformation( oneUnnamedLeaf(),
                                                                                                  false );
            if ( nf2.size() != 1 ) {
                return fail( "an unnamed external leaf must yield exactly one not-found entry, got " + nf2.size() );
            }

            // The tree is not mutated/pruned (the old "delete unresolved nodes" behavior is gone).
            final Phylogeny phy = nestedTree();
            TaxonomyDataManager.obtainDetailedTaxonomicInformation( phy, false );
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return fail( "no node may be deleted; expected 3 external nodes, got "
                        + phy.getNumberOfExternalNodes() );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( e.toString() );
        }
    }

    // root[ROOT] -> { inner[INNER] -> { A, B }, C }
    private static Phylogeny nestedTree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.setName( "ROOT" );
        final PhylogenyNode inner = new PhylogenyNode();
        inner.setName( "INNER" );
        inner.addAsChild( leaf( "A" ) );
        inner.addAsChild( leaf( "B" ) );
        root.addAsChild( inner );
        root.addAsChild( leaf( "C" ) );
        return phylogeny( root );
    }

    private static Phylogeny oneUnnamedLeaf() {
        final PhylogenyNode root = new PhylogenyNode();
        root.setName( "ROOT" );
        root.addAsChild( new PhylogenyNode() ); // single unnamed external leaf; root stays internal
        return phylogeny( root );
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static Phylogeny phylogeny( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "\nTaxonomyDataManagerTest failed: " + msg );
        return false;
    }

    private TaxonomyDataManagerTest() {
    }
}
