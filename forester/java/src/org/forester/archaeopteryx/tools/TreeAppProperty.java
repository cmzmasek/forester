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

import java.util.Iterator;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * Reads and writes a single Archaeopteryx app-state string property keyed by an {@code aptx:<name>} ref at the
 * PHYLOGENY level -- a phyloXML {@code <property applies_to="phylogeny">} that is a direct child of {@code <phylogeny>}
 * (tree-level viewer/app metadata, NOT node data). Used to persist app state (the annotation-import profile, the
 * time-axis config, ...) inside a saved tree file, cleanly separated from the biological tree: the {@code <clade>} data
 * is left untouched, and an {@code aptx:} ref is auto-hidden from node displays by
 * {@code TreePanelUtil.isInternalPropertyRef}.
 *
 * <p>{@link #read} is backward-compatible: it first looks for the property at the phylogeny level (the current
 * location), then falls back to scanning the node properties (the pre-0.11.x location, where such state rode on the
 * ROOT node). {@link #write} clears BOTH locations before writing, so re-saving an old file migrates the property from
 * the node onto the {@code <phylogeny>} element and leaves no stale node-level copy. Pure and Swing-free.</p>
 *
 * <p><b>Undo invariant:</b> {@link #write} always REPLACES (strip + add a fresh {@link Property}); it never mutates a
 * Property in place. {@link org.forester.phylogeny.Phylogeny#copy() phy.copy()} (the undo-snapshot path) shares the
 * Property objects, so mutating one in place would change the value inside every snapshot -- writers of tree-level app
 * state MUST go through {@code write} (strip-and-add), never {@code getProperties().get(i).setValue(...)}.</p>
 */
public final class TreeAppProperty {

    private TreeAppProperty() {
    }

    /**
     * Write {@code value} under {@code ref} as a PHYLOGENY-level property, first removing any existing property with
     * that {@code ref} at BOTH the phylogeny level AND every node (so a re-save migrates the pre-0.11.x node-level
     * copy and never duplicates). A {@code null} {@code value} just strips it.
     */
    public static void write( final Phylogeny phy, final String ref, final String value ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return;
        }
        stripFromTree( phy, ref );
        if ( value == null ) {
            return;
        }
        if ( phy.getProperties() == null ) {
            phy.setProperties( new PropertiesList() );
        }
        phy.getProperties().addProperty( new Property( ref, value, "", "xsd:string", AppliesTo.PHYLOGENY ) );
    }

    /**
     * Read the value stored under {@code ref}: the phylogeny-level property if present, else (backward-compat) the
     * first node-level property with that {@code ref} anywhere in the tree. {@code null} if neither is found.
     */
    public static String read( final Phylogeny phy, final String ref ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return null;
        }
        if ( phy.getProperties() != null ) {
            for( final Property p : phy.getProperties().getProperties() ) {
                if ( ref.equals( p.getRef() ) ) {
                    return p.getValue();
                }
            }
        }
        // backward-compat: pre-0.11.x, tree-level app state rode on the root NODE (robust to a reroot -> full scan)
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PropertiesList pl = it.next().getNodeData().getProperties();
            if ( pl != null ) {
                for( final Property p : pl.getProperties() ) {
                    if ( ref.equals( p.getRef() ) ) {
                        return p.getValue();
                    }
                }
            }
        }
        return null;
    }

    /** Remove every property with {@code ref} from both the phylogeny level and every node. */
    private static void stripFromTree( final Phylogeny phy, final String ref ) {
        if ( phy.getProperties() != null ) {
            phy.getProperties().getProperties().removeIf( p -> ref.equals( p.getRef() ) );
        }
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PropertiesList pl = it.next().getNodeData().getProperties();
            if ( pl != null ) {
                pl.getProperties().removeIf( p -> ref.equals( p.getRef() ) );
            }
        }
    }
}
