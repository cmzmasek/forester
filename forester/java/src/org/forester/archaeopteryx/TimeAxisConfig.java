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

/**
 * The per-tree Time-Axis configuration (type + the two refinement toggles + the calibration overrides), an immutable
 * value object that can be embedded in a saved phyloXML tree as a single {@code aptx:time_axis} property on the ROOT
 * node -- so a deliberately-chosen axis rides along with the tree and is restored on reload. The property is written
 * only when the config DEVIATES from what {@link AptxUtil#deriveTimeAxisType(Phylogeny)} would auto-derive, so files
 * stay clean; auto-derive is the always-on baseline for unsaved / never-configured trees.
 *
 * <p>Serialization is a versioned, {@code ;}-separated string (mirroring {@code NodeDataImporter.ImportProfile}) -- but
 * with NO escaping, deliberately: the value carries only an enum name, two numbers and two digit flags, none of which
 * can contain whitespace / {@code ;} / {@code |} / {@code ~} / a backslash, so {@code XmlElement}'s
 * whitespace-collapse-on-reload cannot corrupt it (unlike the free-text import profile, which must escape).
 */
final class TimeAxisConfig {

    /** phyloXML node-property ref under which the config is persisted, on the ROOT node (internal, so it stays out of
     *  the tip-facing features and is auto-hidden from displays by {@code TreePanelUtil.isInternalPropertyRef}). */
    static final String            TIME_AXIS_REF = "aptx:time_axis";
    private static final String    VERSION       = "v1";
    private static final String    AUTO_TOKEN    = "AUTO"; // serialized sentinel for a null (follow-auto-derive) type
    private final TIME_AXIS_TYPE   _type;                  // null == follow auto-derive (no explicit type override)
    private final double           _root_age_override;      // Ma; 0 == none (derive from the tree's dates)
    private final double           _present_date_override;  // calendar year; 0 == none
    private final boolean          _grid;
    private final boolean          _ages;

    TimeAxisConfig( final TIME_AXIS_TYPE type, final double root_age_override, final double present_date_override,
                    final boolean grid, final boolean ages ) {
        _type = type; // null is meaningful: the type is not overridden, it follows auto-derive
        _root_age_override = Math.max( 0, root_age_override );
        _present_date_override = Math.max( 0, present_date_override );
        _grid = grid;
        _ages = ages;
    }

    /** The explicit per-tab type override, or {@code null} when the type follows auto-derive (no override). */
    TIME_AXIS_TYPE getType() {
        return _type;
    }

    double getRootAgeOverride() {
        return _root_age_override;
    }

    double getPresentDateOverride() {
        return _present_date_override;
    }

    boolean isGrid() {
        return _grid;
    }

    boolean isAges() {
        return _ages;
    }

    /** True when this config represents the pure auto-derive state: no explicit type override, no refinements, no
     *  calibration overrides. Such a config carries NO deviation from what the tree's own dates would produce, so it is
     *  NOT written to the tree (files stay clean). Deviation is thus decided WITHOUT re-deriving -- so an in-place date
     *  edit (which the identity-keyed derive cache doesn't see) can never make a follow-derive tree persist a type. */
    boolean isDefault() {
        return ( _type == null ) && !_grid && !_ages && ( _root_age_override == 0 ) && ( _present_date_override == 0 );
    }

    String serialize() {
        return VERSION + ";" + ( ( _type == null ) ? AUTO_TOKEN : _type.name() ) + ";"
                + Double.toString( _root_age_override ) + ";" + Double.toString( _present_date_override ) + ";"
                + ( _grid ? "1" : "0" ) + ";" + ( _ages ? "1" : "0" );
    }

    /** Parse a serialized config, or {@code null} on ANY malformation (wrong version / field count, bad enum,
     *  non-numeric) -- mirrors {@code NodeDataImporter.ImportProfile.deserialize}; never throws. */
    static TimeAxisConfig deserialize( final String s ) {
        if ( s == null ) {
            return null;
        }
        final String[] f = s.split( ";", -1 );
        if ( ( f.length != 6 ) || !VERSION.equals( f[ 0 ] ) ) {
            return null;
        }
        try {
            final TIME_AXIS_TYPE type = AUTO_TOKEN.equals( f[ 1 ] ) ? null : TIME_AXIS_TYPE.valueOf( f[ 1 ] );
            final double root_age = Double.parseDouble( f[ 2 ] );
            final double present = Double.parseDouble( f[ 3 ] );
            final boolean grid = "1".equals( f[ 4 ] );
            final boolean ages = "1".equals( f[ 5 ] );
            return new TimeAxisConfig( type, root_age, present, grid, ages );
        }
        catch ( final RuntimeException e ) { // IllegalArgumentException (bad enum) / NumberFormatException
            return null;
        }
    }

    /** Persist {@code cfg} onto {@code phy} as a property on the root node (so it saves with the tree); any prior
     *  time-axis property anywhere in the tree is removed first. A {@code null} cfg just clears it. Structural copy of
     *  {@code NodeDataImporter.writeProfileToTree}. */
    static void writeToTree( final Phylogeny phy, final TimeAxisConfig cfg ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return;
        }
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PropertiesList pl = it.next().getNodeData().getProperties();
            if ( pl != null ) {
                final Iterator<Property> pi = pl.getProperties().iterator();
                while ( pi.hasNext() ) {
                    if ( TIME_AXIS_REF.equals( pi.next().getRef() ) ) {
                        pi.remove();
                    }
                }
            }
        }
        if ( cfg == null ) {
            return;
        }
        final PhylogenyNode root = phy.getRoot();
        PropertiesList pl = root.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            root.getNodeData().setProperties( pl );
        }
        pl.addProperty( new Property( TIME_AXIS_REF, cfg.serialize(), "", "xsd:string", AppliesTo.NODE ) );
    }

    /** Read a persisted config from {@code phy} (scans preorder, robust to a reroot), or {@code null} if none /
     *  unparsable. */
    static TimeAxisConfig readFromTree( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return null;
        }
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PropertiesList pl = it.next().getNodeData().getProperties();
            if ( pl != null ) {
                for( final Property p : pl.getProperties() ) {
                    if ( TIME_AXIS_REF.equals( p.getRef() ) ) {
                        return deserialize( p.getValue() );
                    }
                }
            }
        }
        return null;
    }
}
