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
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.forester.phylogeny.PhylogenyNode;

/**
 * The whole reason for drawing several ranks at once is that a reader should see WHICH family a genus belongs to
 * without consulting a legend. That works only if a genus's colour is visibly a variation of its own family's --
 * so the property pinned down here is comparative: every fine clade's hue must be closer to its own containing
 * clade's hue than to any other containing clade's. Independent per-level palettes fail that instantly.
 * <p>
 * Also guards the honesty case: a clade whose broader rank could not be resolved has no parent to take a hue
 * from, and must be drawn desaturated rather than handed a confident colour that implies a parent it lacks.
 * Pure + headless.
 */
public final class CladeHuePaletteTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CladeHuePalette: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return singleLevelUntouched() && childrenSitInParentHue() && orphansAreDesaturated() && deterministic();
    }

    /** One level has no hierarchy to express, so it must keep the plain distinct-colour palette untouched. */
    private static boolean singleLevelUntouched() {
        final Tree t = new Tree();
        final List<List<CladeBand>> one = Arrays.asList( t.genusBands() );
        if ( CladeHuePalette.hueBand( one ) != one ) {
            return fail( "a single level must be returned unchanged" );
        }
        if ( CladeHuePalette.hueBand( null ) != null ) {
            return fail( "a null level list must be returned unchanged" );
        }
        return true;
    }

    /** The headline property: each genus is nearer its OWN family's hue than any other family's. */
    private static boolean childrenSitInParentHue() {
        final Tree t = new Tree();
        // finest first, which is the order the panel stores and draws them in
        final List<List<CladeBand>> out = CladeHuePalette.hueBand( Arrays.asList( t.genusBands(), t.familyBands() ) );
        final Map<String, Color> genus = colors( out.get( 0 ) );
        final Map<String, Color> family = colors( out.get( 1 ) );
        final Map<String, String> belongs_to = new HashMap<>();
        belongs_to.put( "GenusA1", "FamilyA" );
        belongs_to.put( "GenusA2", "FamilyA" );
        belongs_to.put( "GenusB1", "FamilyB" );
        for ( final Map.Entry<String, String> e : belongs_to.entrySet() ) {
            final float own = hueDistance( genus.get( e.getKey() ), family.get( e.getValue() ) );
            for ( final String other : family.keySet() ) {
                if ( other.equals( e.getValue() ) ) {
                    continue;
                }
                final float alien = hueDistance( genus.get( e.getKey() ), family.get( other ) );
                if ( own >= alien ) {
                    return fail( e.getKey() + " must sit nearer its own " + e.getValue() + " (hue distance " + own
                            + ") than " + other + " (" + alien + ")" );
                }
            }
        }
        // siblings must still be told apart -- a family's genera may share its hue band but not its exact colour
        if ( genus.get( "GenusA1" ).equals( genus.get( "GenusA2" ) ) ) {
            return fail( "two genera of the same family must not get the identical colour" );
        }
        // ... and the finer level is drawn slightly less saturated, so an outer mark reads as the grouping
        if ( saturation( genus.get( "GenusA1" ) ) >= saturation( family.get( "FamilyA" ) ) ) {
            return fail( "the finer level should be no more saturated than the level containing it" );
        }
        return true;
    }

    /** A genus hanging outside every annotated family: no parent hue exists, so it must not pretend to have one. */
    private static boolean orphansAreDesaturated() {
        final Tree t = new Tree();
        final List<CladeBand> genera = t.genusBands();
        genera.add( new CladeBand( "GenusOrphan", Color.BLACK, t.orphan ) ); // sits under no annotated family
        final List<List<CladeBand>> out = CladeHuePalette.hueBand( Arrays.asList( genera, t.familyBands() ) );
        final Map<String, Color> genus = colors( out.get( 0 ) );
        final float orphan_sat = saturation( genus.get( "GenusOrphan" ) );
        if ( orphan_sat >= saturation( genus.get( "GenusA1" ) ) ) {
            return fail( "a clade with no containing clade must be drawn DESATURATED relative to a placed one, got "
                    + orphan_sat + " vs " + saturation( genus.get( "GenusA1" ) ) );
        }
        if ( orphan_sat > 0.35f ) {
            return fail( "an orphan clade should read as unaffiliated, got saturation " + orphan_sat );
        }
        // the placed genera must be unaffected by an orphan joining them
        final Map<String, Color> without = colors(
                CladeHuePalette.hueBand( Arrays.asList( t.genusBands(), t.familyBands() ) ).get( 0 ) );
        if ( !genus.get( "GenusA1" ).equals( without.get( "GenusA1" ) ) ) {
            return fail( "an unaffiliated clade must not shift the colours of the placed ones" );
        }
        return true;
    }

    /** Same tree in, same figure out -- colours are keyed on taxon NAME order, never on map iteration order. */
    private static boolean deterministic() {
        final Tree t = new Tree();
        final Map<String, Color> a = colors(
                CladeHuePalette.hueBand( Arrays.asList( t.genusBands(), t.familyBands() ) ).get( 0 ) );
        // the same bands, offered in a different order (as a HashMap-backed producer well might)
        final List<CladeBand> shuffled = t.genusBands();
        java.util.Collections.reverse( shuffled );
        final Map<String, Color> b = colors(
                CladeHuePalette.hueBand( Arrays.asList( shuffled, t.familyBands() ) ).get( 0 ) );
        if ( !a.equals( b ) ) {
            return fail( "the colours must not depend on the order the bands arrive in: " + a + " vs " + b );
        }
        return true;
    }

    // ---- fixture ---------------------------------------------------------------------------------------

    /** root -&gt; (FamilyA -&gt; GenusA1, GenusA2), (FamilyB -&gt; GenusB1), plus a clade under no family at all. */
    private static final class Tree {

        private final PhylogenyNode famA = new PhylogenyNode();
        private final PhylogenyNode famB = new PhylogenyNode();
        private final PhylogenyNode genA1 = new PhylogenyNode();
        private final PhylogenyNode genA2 = new PhylogenyNode();
        private final PhylogenyNode genB1 = new PhylogenyNode();
        private final PhylogenyNode orphan = new PhylogenyNode();

        Tree() {
            final PhylogenyNode root = new PhylogenyNode();
            famA.addAsChild( genA1 );
            famA.addAsChild( genA2 );
            famB.addAsChild( genB1 );
            root.addAsChild( famA );
            root.addAsChild( famB );
            root.addAsChild( orphan );
        }

        List<CladeBand> familyBands() {
            final List<CladeBand> out = new ArrayList<>();
            out.add( new CladeBand( "FamilyA", Color.BLACK, famA ) );
            out.add( new CladeBand( "FamilyB", Color.BLACK, famB ) );
            return out;
        }

        List<CladeBand> genusBands() {
            final List<CladeBand> out = new ArrayList<>();
            out.add( new CladeBand( "GenusA1", Color.BLACK, genA1 ) );
            out.add( new CladeBand( "GenusA2", Color.BLACK, genA2 ) );
            out.add( new CladeBand( "GenusB1", Color.BLACK, genB1 ) );
            return out;
        }
    }

    private static Map<String, Color> colors( final List<CladeBand> bands ) {
        final Map<String, Color> out = new java.util.TreeMap<>();
        for ( final CladeBand b : bands ) {
            out.put( b.getTaxon(), b.getColor() );
        }
        return out;
    }

    /** Distance around the hue wheel (0..0.5), so red at 0.99 and red at 0.01 count as neighbours. */
    private static float hueDistance( final Color a, final Color b ) {
        final float d = Math.abs( CladeHuePalette.hueOf( a ) - CladeHuePalette.hueOf( b ) );
        return Math.min( d, 1.0f - d );
    }

    private static float saturation( final Color c ) {
        return Color.RGBtoHSB( c.getRed(), c.getGreen(), c.getBlue(), null )[ 1 ];
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [CladeHuePaletteTest] " + msg );
        return false;
    }

    private CladeHuePaletteTest() {
    }
}
