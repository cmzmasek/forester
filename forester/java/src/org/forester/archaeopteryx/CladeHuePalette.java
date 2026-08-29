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
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.phylogeny.PhylogenyNode;

/**
 * The colour scheme for NESTED clade levels: each clade at the broadest annotated rank owns a slice of the hue
 * wheel, and the finer clades inside it are coloured from within that slice. A reader then sees at a glance that
 * these genera belong to that family -- the point of drawing several ranks at once. Give each level its own
 * independent palette instead and the figure is two or three unrelated rainbows stacked side by side, implying
 * relationships that are not there.
 * <p>
 * Containment is read from the TREE, not from the names: a fine clade belongs to the coarse clade whose root is
 * its ancestor. That keeps the colouring honest when the taxonomy is patchy -- a clade whose broader rank could
 * not be resolved has no containing clade, and is deliberately drawn DESATURATED ({@link #ORPHAN_SATURATION})
 * rather than being given a confident hue that would claim a parent it does not have.
 * <p>
 * Only applies from two levels up; a single level keeps the plain distinct-colour palette, since with no
 * hierarchy to express there is nothing to band.
 * <p>
 * Pure and deterministic: taxa are processed in name order, so the same tree always yields the same figure.
 */
final class CladeHuePalette {

    /** How much of a parent's hue slice its children may spread across. Below 1 so neighbouring families' children
     *  cannot run into each other; well above 0 so sibling genera are actually told apart by hue. */
    private static final float CHILD_SPREAD     = 0.86f;
    private static final float BASE_SATURATION  = 0.62f;
    private static final float BASE_BRIGHTNESS  = 0.78f;
    /** Each level inward loses a little saturation, so an outer mark reads as the grouping and an inner one as
     *  its subdivision rather than the two competing. */
    private static final float SATURATION_STEP  = 0.05f;
    /** How far a parent's children ramp in lightness across the group. Hue alone is not enough when a clade has
     *  many children -- a big family's genera would otherwise all land within a few degrees of each other. */
    private static final float SIBLING_LIGHTNESS_RANGE = 0.26f;
    /** A clade with no containing clade at the next-broader rank: honest about the gap, quiet on the page. */
    private static final float ORPHAN_SATURATION = 0.22f;

    /** One taxon's claim on the hue wheel: the centre of its slice and how wide that slice is. Carrying the WIDTH
     *  (not just the colour) is what lets the next level subdivide inside its parent instead of guessing. */
    private static final class Arc {

        private final float _hue;
        private final float _width;

        Arc( final float hue, final float width ) {
            _hue = hue;
            _width = width;
        }
    }

    /**
     * Re-colours nested clade levels so each level's colours sit inside their containing clade's hue slice.
     *
     * @param levels the bands per level, FINEST FIRST (the order they are drawn in, innermost nearest the tips)
     * @return the same levels, same order, with hue-banded colours; the input is returned unchanged when there
     *         is only one level (or none), which is the case that has no hierarchy to show
     */
    static List<List<CladeBand>> hueBand( final List<List<CladeBand>> levels ) {
        if ( ( levels == null ) || ( levels.size() < 2 ) ) {
            return levels;
        }
        // work broadest-first: the broadest rank allocates the wheel, each finer rank subdivides what it inherits
        final List<List<CladeBand>> broad_first = new ArrayList<>( levels );
        Collections.reverse( broad_first );
        final int n = broad_first.size();
        // which coarse clade each fine clade belongs to, read from the tree, one map per subdivision
        final List<Map<String, String>> parent_of = new ArrayList<>();
        parent_of.add( Collections.<String, String> emptyMap() ); // the broadest rank has no parent
        for ( int d = 1; d < n; d++ ) {
            parent_of.add( parentTaxa( broad_first.get( d ), broad_first.get( d - 1 ) ) );
        }
        final List<Map<String, Arc>> arcs = new ArrayList<>();
        final List<Map<String, Color>> colors = new ArrayList<>();
        for ( int d = 0; d < n; d++ ) {
            // a clade's slice of the wheel is proportional to how many clades it will have to hold at the next
            // level down -- otherwise a family with eight genera and a family with one get the same arc, and the
            // eight genera end up within a few degrees of each other
            final Map<String, Integer> weights = ( d + 1 < n ) ? childCounts( broad_first.get( d ),
                                                                             parent_of.get( d + 1 ) )
                                                              : Collections.<String, Integer> emptyMap();
            final Map<String, Arc> level_arcs = ( d == 0 ) ? topLevelArcs( broad_first.get( 0 ), weights )
                    : nestedArcs( broad_first.get( d ), parent_of.get( d ), arcs.get( d - 1 ), weights );
            arcs.add( level_arcs );
            colors.add( toColors( broad_first.get( d ), level_arcs, parent_of.get( d ), d ) );
        }
        final List<List<CladeBand>> out = new ArrayList<>();
        for ( int d = 0; d < n; d++ ) {
            out.add( recolor( broad_first.get( d ), colors.get( d ) ) );
        }
        Collections.reverse( out ); // back to finest-first, the order the caller draws in
        return out;
    }

    /** How many distinct finer clades each coarse clade contains (at least 1, so a childless clade still gets a slice). */
    private static Map<String, Integer> childCounts( final List<CladeBand> bands,
                                                     final Map<String, String> child_parent ) {
        final Map<String, Integer> out = new TreeMap<>();
        for ( final String taxon : taxaOf( bands ) ) {
            out.put( taxon, 0 );
        }
        for ( final String parent : child_parent.values() ) {
            if ( out.containsKey( parent ) ) {
                out.put( parent, out.get( parent ) + 1 );
            }
        }
        for ( final Map.Entry<String, Integer> e : out.entrySet() ) {
            e.setValue( Math.max( 1, e.getValue() ) );
        }
        return out;
    }

    /** The broadest rank: its clades divide the whole wheel, each taking a share of it proportional to its weight.
     *  Name order throughout, so the same tree always produces the same figure. */
    private static Map<String, Arc> topLevelArcs( final List<CladeBand> bands, final Map<String, Integer> weights ) {
        final TreeSet<String> taxa = taxaOf( bands );
        final Map<String, Arc> out = new LinkedHashMap<>();
        int total = 0;
        for ( final String taxon : taxa ) {
            total += weight( weights, taxon );
        }
        final float scale = ( total <= 0 ) ? 0f : ( 1.0f / total );
        float cursor = 0f;
        for ( final String taxon : taxa ) {
            final float width = weight( weights, taxon ) * scale;
            out.put( taxon, new Arc( wrap( cursor + ( width / 2f ) ), width ) );
            cursor += width;
        }
        return out;
    }

    /** A finer rank: every clade is placed INSIDE the arc of the clade that contains it, sharing that arc with its
     *  siblings in proportion to what each of them will itself have to hold further in. */
    private static Map<String, Arc> nestedArcs( final List<CladeBand> bands, final Map<String, String> parent_of,
                                                final Map<String, Arc> parent_arcs,
                                                final Map<String, Integer> weights ) {
        final Map<String, TreeSet<String>> by_parent = new TreeMap<>();
        for ( final String taxon : taxaOf( bands ) ) {
            final String parent = parent_of.get( taxon );
            if ( ( parent != null ) && parent_arcs.containsKey( parent ) ) {
                by_parent.computeIfAbsent( parent, k -> new TreeSet<>() ).add( taxon );
            }
        }
        final Map<String, Arc> out = new LinkedHashMap<>();
        for ( final Map.Entry<String, TreeSet<String>> e : by_parent.entrySet() ) {
            final Arc parent_arc = parent_arcs.get( e.getKey() );
            final float inner = parent_arc._width * CHILD_SPREAD;
            int total = 0;
            for ( final String child : e.getValue() ) {
                total += weight( weights, child );
            }
            final float scale = ( total <= 0 ) ? 0f : ( inner / total );
            float cursor = parent_arc._hue - ( inner / 2f );
            for ( final String child : e.getValue() ) {
                final float width = weight( weights, child ) * scale;
                out.put( child, new Arc( wrap( cursor + ( width / 2f ) ), width ) );
                cursor += width;
            }
        }
        return out;
    }

    /**
     * Turns each clade's arc into its colour. Siblings additionally ramp in lightness across their group, because
     * hue alone runs out when a clade has many children. A clade with no containing clade at the broader rank has
     * no arc, and is drawn desaturated rather than being given a hue that would claim a parent it does not have.
     *
     * @param depth how many levels in from the broadest (0 = the broadest itself)
     */
    private static Map<String, Color> toColors( final List<CladeBand> bands, final Map<String, Arc> arcs,
                                                final Map<String, String> parent_of, final int depth ) {
        final float saturation = clamp( BASE_SATURATION - ( SATURATION_STEP * depth ) );
        final Map<String, TreeSet<String>> siblings = new TreeMap<>();
        for ( final String taxon : taxaOf( bands ) ) {
            final String parent = parent_of.get( taxon );
            siblings.computeIfAbsent( ( parent == null ) ? "" : parent, k -> new TreeSet<>() ).add( taxon );
        }
        final Map<String, Color> out = new LinkedHashMap<>();
        final TreeSet<String> orphans = new TreeSet<>();
        for ( final TreeSet<String> group : siblings.values() ) {
            final List<String> members = new ArrayList<>( group );
            for ( int j = 0; j < members.size(); j++ ) {
                final String taxon = members.get( j );
                final Arc arc = arcs.get( taxon );
                if ( arc == null ) {
                    orphans.add( taxon );
                    continue;
                }
                final float t = ( members.size() == 1 ) ? 0.5f : ( j / ( float ) ( members.size() - 1 ) );
                final float bright = clamp( ( BASE_BRIGHTNESS - ( SIBLING_LIGHTNESS_RANGE / 2f ) )
                        + ( t * SIBLING_LIGHTNESS_RANGE ) );
                out.put( taxon, Color.getHSBColor( arc._hue, saturation, bright ) );
            }
        }
        int i = 0;
        for ( final String orphan : orphans ) { // no containing clade at the broader rank -- say so, quietly
            out.put( orphan, Color.getHSBColor( ( float ) i / Math.max( 1, orphans.size() ), ORPHAN_SATURATION,
                                                BASE_BRIGHTNESS ) );
            i++;
        }
        return out;
    }

    private static int weight( final Map<String, Integer> weights, final String taxon ) {
        final Integer w = weights.get( taxon );
        return ( w == null ) ? 1 : Math.max( 1, w.intValue() );
    }

    /**
     * Which coarse clade contains each fine one, read from the tree: a fine band belongs to the coarse band whose
     * root is an ancestor of (or is) its own root. Where a taxon appears in several places (a paraphyletic group
     * yields several bands) the first containing parent in name order wins, so one taxon still gets one colour.
     */
    private static Map<String, String> parentTaxa( final List<CladeBand> bands, final List<CladeBand> parents ) {
        final Map<String, String> out = new TreeMap<>();
        for ( final CladeBand band : bands ) {
            if ( out.containsKey( band.getTaxon() ) ) {
                continue;
            }
            for ( final CladeBand parent : parents ) {
                if ( isDescendantOrSelf( band.getRoot(), parent.getRoot() ) ) {
                    out.put( band.getTaxon(), parent.getTaxon() );
                    break;
                }
            }
        }
        return out;
    }

    /** Whether {@code node} is {@code ancestor} or sits somewhere below it. */
    static boolean isDescendantOrSelf( final PhylogenyNode node, final PhylogenyNode ancestor ) {
        if ( ( node == null ) || ( ancestor == null ) ) {
            return false;
        }
        PhylogenyNode n = node;
        while ( n != null ) {
            if ( n == ancestor ) {
                return true;
            }
            n = n.isRoot() ? null : n.getParent();
        }
        return false;
    }

    private static List<CladeBand> recolor( final List<CladeBand> bands, final Map<String, Color> colors ) {
        final List<CladeBand> out = new ArrayList<>();
        for ( final CladeBand band : bands ) {
            final Color c = colors.get( band.getTaxon() );
            out.add( ( c == null ) ? band : new CladeBand( band.getTaxon(), c, band.getRoot() ) );
        }
        return out;
    }

    private static TreeSet<String> taxaOf( final List<CladeBand> bands ) {
        final TreeSet<String> taxa = new TreeSet<>();
        for ( final CladeBand band : bands ) {
            taxa.add( band.getTaxon() );
        }
        return taxa;
    }

    static float hueOf( final Color c ) {
        return Color.RGBtoHSB( c.getRed(), c.getGreen(), c.getBlue(), null )[ 0 ];
    }

    private static float clamp( final float v ) {
        return Math.max( 0f, Math.min( 1f, v ) );
    }

    private static float wrap( final float hue ) {
        final float h = hue % 1.0f;
        return ( h < 0 ) ? ( h + 1.0f ) : h;
    }

    private CladeHuePalette() {
    }
}
