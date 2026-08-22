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

package org.forester.archaeopteryx.tools;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Headless unit test for {@link TreeAppProperty} (the phylogeny-level app-state read/write helper): a new write lands as
 * a single {@code applies_to="phylogeny"} property and on NO node; a re-write replaces (no accumulation); a
 * {@code null} value clears both locations; {@link TreeAppProperty#read} prefers the phylogeny level but falls back to
 * a pre-0.11.x root-NODE copy (backward-compat); a re-save MIGRATES a node-level copy onto the phylogeny; and
 * null/empty trees are safe no-ops.
 */
public final class TreeAppPropertyTest {

    private static final String REF = "aptx:time_axis";

    public static void main( final String[] args ) {
        System.out.println( "TreeAppProperty: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        try {
            // --- write-new: a single phylogeny-level property, on no node ---
            final Phylogeny phy = tree();
            TreeAppProperty.write( phy, REF, "hello" );
            if ( !"hello".equals( TreeAppProperty.read( phy, REF ) ) ) {
                return fail( "read should return the written value" );
            }
            if ( !phy.isHasProperties() || ( phyCount( phy ) != 1 )
                    || ( AppliesTo.PHYLOGENY != phy.getProperties().getProperties( REF ).get( 0 ).getAppliesTo() ) ) {
                return fail( "a write should create exactly one applies_to=phylogeny property" );
            }
            if ( nodeCount( phy ) != 0 ) {
                return fail( "a write must not put the property on any node" );
            }

            // --- re-write REPLACES (no accumulation) ---
            TreeAppProperty.write( phy, REF, "world" );
            if ( ( phyCount( phy ) != 1 ) || !"world".equals( TreeAppProperty.read( phy, REF ) ) ) {
                return fail( "a re-write must replace, not accumulate: count=" + phyCount( phy ) );
            }

            // --- write(null) clears it ---
            TreeAppProperty.write( phy, REF, null );
            if ( ( TreeAppProperty.read( phy, REF ) != null ) || ( phyCount( phy ) != 0 ) ) {
                return fail( "write(null) should clear the property" );
            }

            // --- backward-compat: read falls back to a root-NODE copy (pre-0.11.x location) ---
            final Phylogeny old = tree();
            addNodeProperty( old, "legacy" );
            if ( old.isHasProperties() ) {
                return fail( "precondition: no phylogeny-level property yet" );
            }
            if ( !"legacy".equals( TreeAppProperty.read( old, REF ) ) ) {
                return fail( "read must fall back to a node-level (pre-0.11.x) copy" );
            }

            // --- a re-save MIGRATES it: phylogeny-level set, node-level gone ---
            TreeAppProperty.write( old, REF, TreeAppProperty.read( old, REF ) );
            if ( ( phyCount( old ) != 1 ) || ( nodeCount( old ) != 0 ) ) {
                return fail( "a re-save should migrate the node-level copy onto the phylogeny and drop the node copy" );
            }

            // --- read precedence: phylogeny-level wins over a lingering node-level copy ---
            final Phylogeny both = tree();
            addNodeProperty( both, "NODE" );
            TreeAppProperty.write( both, REF, "PHY" ); // strips the node copy AND adds the phylogeny copy...
            addNodeProperty( both, "NODE2" );          // ...then a stray node copy is re-introduced
            if ( !"PHY".equals( TreeAppProperty.read( both, REF ) ) ) {
                return fail( "read must prefer the phylogeny-level value over a node-level one" );
            }

            // --- null / empty trees are safe no-ops ---
            if ( TreeAppProperty.read( null, REF ) != null ) {
                return fail( "read(null) must be null" );
            }
            final Phylogeny empty = new Phylogeny();
            TreeAppProperty.write( empty, REF, "x" ); // no-op, no throw
            if ( TreeAppProperty.read( empty, REF ) != null ) {
                return fail( "read on an empty tree must be null" );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected exception: " + e.getMessage() );
        }
    }

    private static Phylogeny tree() throws Exception {
        return Phylogeny.createInstanceFromNhxString( "((A,B),C);" );
    }

    /** Attach a NODE-level property with {@link #REF} onto the root (the pre-0.11.x location). */
    private static void addNodeProperty( final Phylogeny phy, final String value ) {
        PropertiesList pl = phy.getRoot().getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            phy.getRoot().getNodeData().setProperties( pl );
        }
        pl.addProperty( new Property( REF, value, "", "xsd:string", AppliesTo.NODE ) );
    }

    private static int phyCount( final Phylogeny phy ) {
        return ( phy.getProperties() == null ) ? 0 : phy.getProperties().getProperties( REF ).size();
    }

    private static int nodeCount( final Phylogeny phy ) {
        int n = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PropertiesList pl = it.next().getNodeData().getProperties();
            if ( pl != null ) {
                n += pl.getProperties( REF ).size();
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [TreeAppPropertyTest] " + msg );
        return false;
    }

    private TreeAppPropertyTest() {
    }
}
