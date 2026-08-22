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

import java.util.Iterator;

import org.forester.archaeopteryx.Options.TIME_AXIS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/** Headless unit test for {@link TimeAxisConfig}: the serialize/deserialize codec (round-trip + malformed -> null) and
 *  the write/read of the {@code aptx:time_axis} property -- stored at the PHYLOGENY level (via {@link org.forester.archaeopteryx.tools.TreeAppProperty}),
 *  read backward-compatibly from the pre-0.11.76 root-node location, and migrated onto {@code <phylogeny>} on re-save. */
public final class TimeAxisConfigTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TimeAxisConfig: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        // round-trip every field (including the two calibration overrides and the flags)
        final TimeAxisConfig a = new TimeAxisConfig( TIME_AXIS_TYPE.CALENDAR, 250.5, 2021.25, true, false );
        final TimeAxisConfig b = TimeAxisConfig.deserialize( a.serialize() );
        if ( b == null ) {
            return fail( "round-trip returned null for: " + a.serialize() );
        }
        if ( ( b.getType() != TIME_AXIS_TYPE.CALENDAR ) || ( Math.abs( b.getRootAgeOverride() - 250.5 ) > 1e-9 )
                || ( Math.abs( b.getPresentDateOverride() - 2021.25 ) > 1e-9 ) || !b.isGrid() || b.isAges() ) {
            return fail( "round-trip lost a field: " + b.serialize() );
        }
        // a null (follow-auto-derive) type round-trips as null -- the config can represent "no type override"
        final TimeAxisConfig auto_grid = TimeAxisConfig.deserialize(
                new TimeAxisConfig( null, 0, 0, true, false ).serialize() );
        if ( ( auto_grid == null ) || ( auto_grid.getType() != null ) || !auto_grid.isGrid() ) {
            return fail( "a null (auto-derive) type with a grid deviation must round-trip as {type=null, grid=true}" );
        }
        // isDefault(): pure auto-derive == null type + no flags + no overrides (this is what does NOT get written)
        if ( !new TimeAxisConfig( null, 0, 0, false, false ).isDefault() ) {
            return fail( "a null type with no flags/overrides is the pure auto-derive default" );
        }
        if ( new TimeAxisConfig( TIME_AXIS_TYPE.GEOLOGIC, 0, 0, false, false ).isDefault() ) {
            return fail( "an EXPLICIT type is a deviation, not the auto-derive default" );
        }
        if ( new TimeAxisConfig( null, 0, 0, true, false ).isDefault() ) {
            return fail( "a grid-only deviation (with auto type) is NOT the default" );
        }
        final TimeAxisConfig def = new TimeAxisConfig( TIME_AXIS_TYPE.GEOLOGIC, 0, 0, false, false ); // explicit type
        // malformed inputs all yield null (never throw)
        final String[] bad = { null, "", "junk", "v2;GEOLOGIC;0.0;0.0;0;0", // wrong version
                "v1;GEOLOGIC;0.0;0.0;0", // 5 fields
                "v1;GEOLOGIC;0.0;0.0;0;0;0", // 7 fields
                "v1;BOGUS;0.0;0.0;0;0", // bad enum
                "v1;GEOLOGIC;abc;0.0;0;0" }; // non-numeric
        for ( final String s : bad ) {
            if ( TimeAxisConfig.deserialize( s ) != null ) {
                return fail( "malformed input must deserialize to null: " + s );
            }
        }
        // write/read on a real tree: strip-then-add at the PHYLOGENY level, round-trips; null clears
        final Phylogeny phy = twoTipTree();
        TimeAxisConfig.writeToTree( phy, a );
        final TimeAxisConfig read = TimeAxisConfig.readFromTree( phy );
        if ( ( read == null ) || ( read.getType() != TIME_AXIS_TYPE.CALENDAR ) || !read.isGrid() ) {
            return fail( "writeToTree/readFromTree did not round-trip the config" );
        }
        // the property lives at the PHYLOGENY level (<property applies_to="phylogeny">), NOT on any node
        if ( !hasTimeAxisPropAtPhylogenyLevel( phy ) ) {
            return fail( "the time-axis property must be written at the phylogeny level" );
        }
        if ( anyNodeHasTimeAxisProp( phy ) ) {
            return fail( "the time-axis property must NOT ride on any node" );
        }
        // a second write REPLACES (no duplicate at the phylogeny level); null clears entirely
        TimeAxisConfig.writeToTree( phy, def );
        final TimeAxisConfig read2 = TimeAxisConfig.readFromTree( phy );
        if ( ( read2 == null ) || ( read2.getType() != TIME_AXIS_TYPE.GEOLOGIC ) ) {
            return fail( "a second writeToTree must replace the config, got " + ( read2 == null ? "null" : read2.serialize() ) );
        }
        if ( countTimeAxisPropsAtPhylogenyLevel( phy ) != 1 ) {
            return fail( "a second writeToTree must not leave a duplicate phylogeny-level property, found "
                    + countTimeAxisPropsAtPhylogenyLevel( phy ) );
        }
        TimeAxisConfig.writeToTree( phy, null );
        if ( TimeAxisConfig.readFromTree( phy ) != null ) {
            return fail( "writeToTree(null) must strip the property" );
        }
        if ( hasTimeAxisPropAtPhylogenyLevel( phy ) ) {
            return fail( "writeToTree(null) must leave no phylogeny-level property" );
        }
        // backward-compat + migration: a pre-0.11.76 tree stored the config on the ROOT node. readFromTree must still
        // read it, and re-saving must MIGRATE it to the phylogeny level and drop the stale node copy.
        final Phylogeny legacy = twoTipTree();
        writeLegacyRootNodeProperty( legacy, a.serialize() );
        final TimeAxisConfig read_legacy = TimeAxisConfig.readFromTree( legacy );
        if ( ( read_legacy == null ) || ( read_legacy.getType() != TIME_AXIS_TYPE.CALENDAR ) ) {
            return fail( "readFromTree must read a pre-0.11.76 root-node property (backward-compat)" );
        }
        TimeAxisConfig.writeToTree( legacy, def ); // re-save
        if ( !hasTimeAxisPropAtPhylogenyLevel( legacy ) || anyNodeHasTimeAxisProp( legacy ) ) {
            return fail( "re-saving a legacy file must migrate the property to the phylogeny level and clear the node copy" );
        }
        return true;
    }

    /** True when a {@code TIME_AXIS_REF} property with {@code applies_to="phylogeny"} sits at the phylogeny level. */
    private static boolean hasTimeAxisPropAtPhylogenyLevel( final Phylogeny phy ) {
        return countTimeAxisPropsAtPhylogenyLevel( phy ) > 0;
    }

    private static int countTimeAxisPropsAtPhylogenyLevel( final Phylogeny phy ) {
        if ( phy.getProperties() == null ) {
            return 0;
        }
        int count = 0;
        for ( final Property p : phy.getProperties().getProperties() ) {
            if ( TimeAxisConfig.TIME_AXIS_REF.equals( p.getRef() ) && ( p.getAppliesTo() == AppliesTo.PHYLOGENY ) ) {
                count++;
            }
        }
        return count;
    }

    private static boolean anyNodeHasTimeAxisProp( final Phylogeny phy ) {
        for ( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            if ( hasTimeAxisProp( it.next() ) ) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasTimeAxisProp( final PhylogenyNode n ) {
        if ( n.getNodeData().getProperties() == null ) {
            return false;
        }
        for ( final Property p : n.getNodeData().getProperties().getProperties() ) {
            if ( TimeAxisConfig.TIME_AXIS_REF.equals( p.getRef() ) ) {
                return true;
            }
        }
        return false;
    }

    /** Write a config the pre-0.11.76 way: an {@code applies_to="node"} property on the ROOT node (the old location
     *  the backward-compat read + re-save migration must handle). */
    private static void writeLegacyRootNodeProperty( final Phylogeny phy, final String value ) {
        final PhylogenyNode root = phy.getRoot();
        PropertiesList pl = root.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            root.getNodeData().setProperties( pl );
        }
        pl.addProperty( new Property( TimeAxisConfig.TIME_AXIS_REF, value, "", "xsd:string", AppliesTo.NODE ) );
    }

    private static Phylogeny twoTipTree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( new PhylogenyNode() );
        root.addAsChild( new PhylogenyNode() );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String m ) {
        System.out.println( "  TimeAxisConfigTest: " + m );
        return false;
    }
}
