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

package org.forester.phylogeny;

import java.util.LinkedHashSet;
import java.util.Set;

/**
 * Headless test that pruning external nodes gives the same tree whatever order the tips are deleted in:
 * {@link PhylogenyMethods#deleteExternalNodesNegativeSelection(Set, Phylogeny)}, its by-name overload, and
 * {@link PhylogenyMethods#deleteExternalNodesPositiveSelection(String[], Phylogeny)} (which delegates to it). The pruned
 * root keeps the ORIGINAL root's own branch length; every other length is the sum along the collapsed path. Before the
 * fix, keeping only C of {@link #TREE} left C as the root with length 0.05, 0.25 or 0.45 depending on the order (the
 * id overload iterates a HashSet, i.e. the node ids assigned at load). Joint rule with Archaeopteryx.js (Select
 * Representative Tips), Christian 2026-09-15.
 */
public final class ExternalNodeDeletionTest {

    private static final String     TREE             = "(A:0.05,(B:0.1,(C:0.05,D:0.1):0.2):0.2);";
    private static final String     TREE_ROOT_LENGTH = "(A:0.05,(B:0.1,(C:0.05,D:0.1):0.2):0.2):0.7;";
    // every order of deleting A, B and D (keeps only C)
    private static final String[][] ORDERS_ABD       = { { "A", "B", "D" }, { "A", "D", "B" }, { "B", "A", "D" },
            { "B", "D", "A" }, { "D", "A", "B" }, { "D", "B", "A" } };
    private static final double     EPS              = 1e-12;

    public static void main( final String[] args ) {
        System.out.println( "ExternalNodeDeletion: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        try {
            boolean ok = true;
            ok &= testRootKeepsOriginalRootDistanceById( TREE );
            ok &= testRootKeepsOriginalRootDistanceById( TREE_ROOT_LENGTH );
            ok &= testOriginalRootLengthIsCarried();
            ok &= testInteriorLengthsAreSummed();
            ok &= testByNameOverload();
            ok &= testPositiveSelection();
            ok &= testNothingDeleted();
            ok &= testEverythingDeleted();
            return ok;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected exception: " + e.getMessage() );
        }
    }

    private static boolean testRootKeepsOriginalRootDistanceById( final String nhx ) throws Exception {
        boolean ok = true;
        for( final String[] order : ORDERS_ABD ) {
            final Phylogeny phy = Phylogeny.createInstanceFromNhxString( nhx );
            final double original = phy.getRoot().getDistanceToParent();
            final Set<Long> to_delete = new LinkedHashSet<>();
            for( final String name : order ) {
                to_delete.add( phy.getNode( name ).getId() );
            }
            PhylogenyMethods.deleteExternalNodesNegativeSelection( to_delete, phy );
            ok &= checkOnlyC( phy, original, "ids " + String.join( ",", order ) + " from " + nhx );
        }
        return ok;
    }

    // the fixture must actually carry a root length, or TREE_ROOT_LENGTH would test the same thing as TREE
    private static boolean testOriginalRootLengthIsCarried() throws Exception {
        final Phylogeny phy = Phylogeny.createInstanceFromNhxString( TREE_ROOT_LENGTH );
        final Set<Long> to_delete = new LinkedHashSet<>();
        for( final String name : new String[] { "D", "B", "A" } ) {
            to_delete.add( phy.getNode( name ).getId() );
        }
        PhylogenyMethods.deleteExternalNodesNegativeSelection( to_delete, phy );
        if ( Math.abs( phy.getRoot().getDistanceToParent() - 0.7 ) > EPS ) {
            return fail( "root length should stay 0.7, is " + phy.getRoot().getDistanceToParent() );
        }
        return true;
    }

    // keeping B and C: the root has B (0.1) and C (0.05 + 0.2) in either deletion order
    private static boolean testInteriorLengthsAreSummed() throws Exception {
        boolean ok = true;
        for( final String[] order : new String[][] { { "A", "D" }, { "D", "A" } } ) {
            final Phylogeny phy = Phylogeny.createInstanceFromNhxString( TREE );
            final double original = phy.getRoot().getDistanceToParent();
            final Set<Long> to_delete = new LinkedHashSet<>();
            for( final String name : order ) {
                to_delete.add( phy.getNode( name ).getId() );
            }
            PhylogenyMethods.deleteExternalNodesNegativeSelection( to_delete, phy );
            final String what = "ids " + String.join( ",", order ) + ": ";
            if ( ( phy.getNumberOfExternalNodes() != 2 ) || ( phy.getRoot().getNumberOfDescendants() != 2 ) ) {
                ok = fail( what + "expected a root with the two tips B and C" );
                continue;
            }
            if ( Math.abs( phy.getNode( "B" ).getDistanceToParent() - 0.1 ) > EPS ) {
                ok = fail( what + "B should keep 0.1, has " + phy.getNode( "B" ).getDistanceToParent() );
            }
            if ( Math.abs( phy.getNode( "C" ).getDistanceToParent() - 0.25 ) > EPS ) {
                ok = fail( what + "C should have 0.05 + 0.2, has " + phy.getNode( "C" ).getDistanceToParent() );
            }
            if ( phy.getRoot().getDistanceToParent() != original ) {
                ok = fail( what + "root length " + phy.getRoot().getDistanceToParent() + ", expected " + original );
            }
        }
        return ok;
    }

    private static boolean testByNameOverload() throws Exception {
        boolean ok = true;
        for( final String nhx : new String[] { TREE, TREE_ROOT_LENGTH } ) {
            for( final String[] order : ORDERS_ABD ) {
                final Phylogeny phy = Phylogeny.createInstanceFromNhxString( nhx );
                final double original = phy.getRoot().getDistanceToParent();
                PhylogenyMethods.deleteExternalNodesNegativeSelection( order.clone(), phy );
                ok &= checkOnlyC( phy, original, "names " + String.join( ",", order ) + " from " + nhx );
            }
        }
        return ok;
    }

    private static boolean testPositiveSelection() throws Exception {
        final Phylogeny phy = Phylogeny.createInstanceFromNhxString( TREE_ROOT_LENGTH );
        final double original = phy.getRoot().getDistanceToParent();
        PhylogenyMethods.deleteExternalNodesPositiveSelection( new String[] { "C" }, phy );
        return checkOnlyC( phy, original, "positive selection of C from " + TREE_ROOT_LENGTH );
    }

    // deliberate non-behaviour: deleting nothing changes nothing
    private static boolean testNothingDeleted() throws Exception {
        final Phylogeny phy = Phylogeny.createInstanceFromNhxString( TREE_ROOT_LENGTH );
        PhylogenyMethods.deleteExternalNodesNegativeSelection( new LinkedHashSet<Long>(), phy );
        if ( ( phy.getNumberOfExternalNodes() != 4 ) || ( Math.abs( phy.getRoot().getDistanceToParent() - 0.7 ) > EPS )
                || ( Math.abs( phy.getNode( "C" ).getDistanceToParent() - 0.05 ) > EPS ) ) {
            return fail( "deleting nothing changed the tree" );
        }
        return true;
    }

    // edge case: deleting every tip empties the tree without an exception (no root to restore a length on)
    private static boolean testEverythingDeleted() throws Exception {
        final Phylogeny phy = Phylogeny.createInstanceFromNhxString( "(A:0.1,B:0.2):0.3;" );
        final Set<Long> to_delete = new LinkedHashSet<>();
        to_delete.add( phy.getNode( "A" ).getId() );
        to_delete.add( phy.getNode( "B" ).getId() );
        PhylogenyMethods.deleteExternalNodesNegativeSelection( to_delete, phy );
        if ( !phy.isEmpty() ) {
            return fail( "deleting every tip should leave an empty tree" );
        }
        return true;
    }

    private static boolean checkOnlyC( final Phylogeny phy, final double original, final String what ) {
        if ( ( phy.getNumberOfExternalNodes() != 1 ) || !"C".equals( phy.getRoot().getName() ) ) {
            return fail( what + ": expected C alone as the root" );
        }
        if ( phy.getRoot().getDistanceToParent() != original ) {
            return fail( what + ": root length " + phy.getRoot().getDistanceToParent() + ", expected the original root's "
                    + original );
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ExternalNodeDeletionTest] " + msg );
        return false;
    }

    private ExternalNodeDeletionTest() {
    }
}
