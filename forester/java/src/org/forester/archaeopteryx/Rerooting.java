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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeData;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/**
 * The rules every way of re-rooting a tree follows -- manual Re-Root, Midpoint-Root, MAD-Root, GSDI and GSDIR --
 * agreed with Archaeopteryx.js (Christian, 2026-09-14): whether the tree may be re-rooted at all, and which internal
 * nodes carry data whose meaning a re-root would change. Also the root-free view of a node, used when a tree declared
 * unrooted is shown in the unrooted layout. Pure/headless.
 */
final class Rerooting {

    static final String NOT_REROOTABLE = "This tree is marked as not re-rootable (rerootable=\"false\").";
    static final String TIME_TREE      = "Time trees can't be re-rooted: their branch lengths are times measured from"
            + " this root.";

    private Rerooting() {
    }

    /**
     * Why {@code phy} must not be re-rooted ({@link #NOT_REROOTABLE}, {@link #TIME_TREE}), or null when it may be.
     */
    static String refusal( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return null;
        }
        if ( !phy.isRerootable() ) {
            return NOT_REROOTABLE;
        }
        if ( AptxUtil.isTimeTree( phy ) ) {
            return TIME_TREE;
        }
        return null;
    }

    /**
     * Whether a node carries data whose meaning depends on its clade: a name, a non-empty taxonomy, a sequence, an
     * event, a distribution, a date, a reference, binary characters, or a node property (not an internal aptx: one, a
     * visual-style one, or one that applies to the parent branch). Branch lengths and support values belong to the
     * branch, and visual styling (colours, node styles, collapsed state) is not data.
     */
    static boolean hasNodeData( final PhylogenyNode n ) {
        if ( !ForesterUtil.isEmpty( n.getName() ) ) {
            return true;
        }
        final NodeData nd = n.getNodeData();
        if ( nd.isHasTaxonomy() ) {
            for( final Taxonomy t : nd.getTaxonomies() ) {
                if ( ( t != null ) && !t.isEmpty() ) {
                    return true;
                }
            }
        }
        if ( nd.isHasSequence() || nd.isHasEvent() || nd.isHasDistribution() || nd.isHasDate() || nd.isHasReference()
                || nd.isHasBinaryCharacters() ) {
            return true;
        }
        if ( nd.isHasProperties() ) {
            for( final Property p : nd.getProperties().getProperties() ) {
                if ( ( p != null ) && !TreePanelUtil.isInternalPropertyRef( p.getRef() )
                        && !TreePanelUtil.isVisualStylePropertyRef( p.getRef() )
                        && ( p.getAppliesTo() != Property.AppliesTo.PARENT_BRANCH ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    /** The number of internal nodes (the root included) that carry data ({@link #hasNodeData}). */
    static int internalNodesWithData( final Phylogeny phy ) {
        int count = 0;
        if ( ( phy != null ) && !phy.isEmpty() ) {
            for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( n.isInternal() && hasNodeData( n ) ) {
                    ++count;
                }
            }
        }
        return count;
    }

    /**
     * The number of data-carrying internal nodes of {@code before} whose clade differs in {@code after}, a re-rooted
     * copy of it (node ids survive {@link Phylogeny#copy()}). A re-root flips parent and child only along the path
     * between the old and the new root, so a node's clade changes exactly when its set of children does -- or when
     * the node is gone (a two-child old root is dissolved).
     */
    static int dataNodesWhoseCladeChanges( final Phylogeny before, final Phylogeny after ) {
        if ( ( before == null ) || before.isEmpty() || ( after == null ) || after.isEmpty() ) {
            return 0;
        }
        final Map<Long, PhylogenyNode> after_by_id = new HashMap<>();
        for( final PhylogenyNodeIterator it = after.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            after_by_id.put( n.getId(), n );
        }
        int count = 0;
        for( final PhylogenyNodeIterator it = before.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isInternal() && hasNodeData( n ) ) {
                final PhylogenyNode m = after_by_id.get( n.getId() );
                if ( ( m == null ) || !childIds( n ).equals( childIds( m ) ) ) {
                    ++count;
                }
            }
        }
        return count;
    }

    private static Set<Long> childIds( final PhylogenyNode n ) {
        final Set<Long> ids = new HashSet<>();
        for( final PhylogenyNode c : n.getDescendants() ) {
            ids.add( c.getId() );
        }
        return ids;
    }

    /** The warning shown before a re-root that changes the clade of {@code affected} of the {@code with_data} nodes. */
    static String dataWarning( final int with_data, final int affected ) {
        final String nodes = with_data + ( ( with_data == 1 ) ? " internal node" : " internal nodes" );
        final String change;
        if ( with_data == 1 ) {
            change = "Re-rooting changes its clade, so its data may no longer describe it.";
        }
        else if ( affected == 1 ) {
            change = "Re-rooting changes the clade of 1 of them, so its data may no longer describe it.";
        }
        else {
            change = "Re-rooting changes the clade of " + affected + " of them, so their data may no longer describe"
                    + " them.";
        }
        return "This tree has data on " + nodes + ". " + change;
    }

    /**
     * Whether values that only mean something relative to a root (distance to parent, depth, tips below, distance
     * from the root, height) are hidden: the tree is shown in the unrooted layout AND its file declared it unrooted.
     */
    static boolean hidesRootDependentValues( final Options.PHYLOGENY_GRAPHICS_TYPE type, final Phylogeny phy ) {
        return ( type == Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && ( phy != null ) && phy.isDeclaredUnrooted();
    }

    /**
     * The number of tips on each side of an internal node -- one count per neighbour (each child's clade, plus the
     * rest of the tree beyond its parent), ascending. Unlike "tips below", this does not depend on where the tree is
     * stored as rooted. Empty for a tip.
     */
    static List<Integer> tipsAround( final PhylogenyNode n ) {
        final List<Integer> sides = new ArrayList<>();
        if ( n.isExternal() ) {
            return sides;
        }
        int below = 0;
        for( final PhylogenyNode c : n.getDescendants() ) {
            final int tips = c.getAllExternalDescendants().size();
            sides.add( tips );
            below += tips;
        }
        if ( !n.isRoot() ) {
            PhylogenyNode root = n;
            while ( !root.isRoot() ) {
                root = root.getParent();
            }
            sides.add( root.getAllExternalDescendants().size() - below );
        }
        Collections.sort( sides );
        return sides;
    }

    /** {@link #tipsAround} as text, e.g. "2 · 3 · 5". */
    static String tipsAroundText( final PhylogenyNode n ) {
        final StringBuilder sb = new StringBuilder();
        for( final int tips : tipsAround( n ) ) {
            if ( sb.length() > 0 ) {
                sb.append( " · " );
            }
            sb.append( tips );
        }
        return sb.toString();
    }
}
