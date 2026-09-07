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

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import org.forester.archaeopteryx.NodeDataDraft.Problem;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.data.Identifier;
import org.forester.util.ForesterUtil;

/**
 * The editable, tree-level metadata of a {@link Phylogeny} as plain text -- the six phyloXML fields a user may set
 * by hand: name, description, identifier (value + provider), type, and the branch-length unit. Like
 * {@link NodeDataDraft} it is a side-effect-free view-model: {@link #from} reads, {@link #validate} checks,
 * {@link #writeTo} applies, and equality of the {@link #normalized} drafts decides whether anything changed.
 * <p>
 * The name is a one-line label (tab title, figure title), so it is normalized by collapsing every run of
 * whitespace and trimming; the description keeps its line breaks and is only trimmed. A tree that HAS a name
 * cannot be renamed to nothing (it would blank the tab and the on-disk name), which {@link #validate} reports.
 */
final class TreePropertiesDraft {

    static final String NAME          = "name";
    static final String DESCRIPTION   = "description";
    static final String ID_VALUE      = "id.value";
    static final String ID_PROVIDER   = "id.provider";
    static final String TYPE          = "type";
    static final String DISTANCE_UNIT = "distance_unit";

    String name          = "";
    String description   = "";
    String idValue       = "";
    String idProvider    = "";
    String type          = "";
    String distanceUnit  = "";

    static TreePropertiesDraft from( final Phylogeny phy ) {
        final TreePropertiesDraft d = new TreePropertiesDraft();
        if ( phy == null ) {
            return d;
        }
        d.name = nn( phy.getName() );
        d.description = nn( phy.getDescription() );
        if ( phy.getIdentifier() != null ) {
            d.idValue = nn( phy.getIdentifier().getValue() );
            d.idProvider = nn( phy.getIdentifier().getProvider() );
        }
        d.type = nn( phy.getType() );
        d.distanceUnit = nn( phy.getDistanceUnit() );
        return d;
    }

    TreePropertiesDraft copy() {
        final TreePropertiesDraft d = new TreePropertiesDraft();
        d.name = name;
        d.description = description;
        d.idValue = idValue;
        d.idProvider = idProvider;
        d.type = type;
        d.distanceUnit = distanceUnit;
        return d;
    }

    /** The draft as it would be written: name whitespace collapsed, everything trimmed. */
    TreePropertiesDraft normalized() {
        final TreePropertiesDraft d = new TreePropertiesDraft();
        d.name = ForesterUtil.collapseWhiteSpace( name ).trim();
        d.description = description.trim();
        d.idValue = idValue.trim();
        d.idProvider = idProvider.trim();
        d.type = ForesterUtil.collapseWhiteSpace( type ).trim();
        d.distanceUnit = ForesterUtil.collapseWhiteSpace( distanceUnit ).trim();
        return d;
    }

    /** Validation problems against {@code baseline} (what the tree currently has); empty means writable. */
    List<Problem> validate( final TreePropertiesDraft baseline ) {
        final List<Problem> out = new ArrayList<>();
        final TreePropertiesDraft n = normalized();
        if ( n.name.isEmpty() && ( baseline != null ) && !baseline.normalized().name.isEmpty() ) {
            out.add( new Problem( NAME, "The name cannot be blank (the tab and the saved tree are named by it)." ) );
        }
        if ( n.idValue.isEmpty() && !n.idProvider.isEmpty() ) {
            out.add( new Problem( ID_PROVIDER, "An identifier provider needs an identifier value." ) );
        }
        return out;
    }

    /** Which fields differ from {@code baseline} after normalization, as user-facing labels (in field order). */
    Set<String> changedFields( final TreePropertiesDraft baseline ) {
        final Set<String> out = new LinkedHashSet<>();
        final TreePropertiesDraft a = normalized();
        final TreePropertiesDraft b = baseline.normalized();
        if ( !a.name.equals( b.name ) ) {
            out.add( "name" );
        }
        if ( !a.description.equals( b.description ) ) {
            out.add( "description" );
        }
        if ( !a.idValue.equals( b.idValue ) || !a.idProvider.equals( b.idProvider ) ) {
            out.add( "identifier" );
        }
        if ( !a.type.equals( b.type ) ) {
            out.add( "type" );
        }
        if ( !a.distanceUnit.equals( b.distanceUnit ) ) {
            out.add( "branch-length unit" );
        }
        return out;
    }

    /** Applies the normalized values to {@code phy} (every field, so a blanked one is cleared). */
    void writeTo( final Phylogeny phy ) {
        final TreePropertiesDraft n = normalized();
        phy.setName( n.name );
        phy.setDescription( n.description );
        phy.setIdentifier( n.idValue.isEmpty() ? null
                : new Identifier( n.idValue, n.idProvider.isEmpty() ? null : n.idProvider ) );
        phy.setType( n.type );
        phy.setDistanceUnit( n.distanceUnit );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( !( o instanceof TreePropertiesDraft ) ) {
            return false;
        }
        final TreePropertiesDraft d = (TreePropertiesDraft) o;
        return name.equals( d.name ) && description.equals( d.description ) && idValue.equals( d.idValue )
                && idProvider.equals( d.idProvider ) && type.equals( d.type ) && distanceUnit.equals( d.distanceUnit );
    }

    @Override
    public int hashCode() {
        return Objects.hash( name, description, idValue, idProvider, type, distanceUnit );
    }

    @Override
    public String toString() {
        return "TreePropertiesDraft[name=" + name + ", id=" + idValue + "/" + idProvider + ", type=" + type
                + ", unit=" + distanceUnit + "]";
    }

    private static String nn( final String s ) {
        return ( s == null ) ? "" : s;
    }
}
