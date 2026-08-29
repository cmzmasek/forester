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
import java.util.Arrays;
import java.util.List;

import org.forester.archaeopteryx.TreePanel.CLADE_LABEL_ANGLE;

/**
 * "Annotate Clades by Rank" can draw up to three ranks at once as nested bar/bracket columns. Which rank ends up
 * where is DERIVED, never chosen: {@link CladeLevel#order} sorts finest-rank-first so the finest is drawn nearest
 * the tips and the broadest outermost, which is the only order in which the nesting reads correctly -- a broad
 * clade's mark has to span the several finer marks inside it. If this ordering breaks, the figure silently lies
 * about the hierarchy, so it is pinned down here. Pure + headless.
 */
public final class CladeLevelTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CladeLevel: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return rankDepth() && ordering();
    }

    private static boolean rankDepth() {
        // finer rank => strictly greater depth, right down the chain a user actually picks
        final String[] broad_to_fine = { "phylum", "class", "order", "family", "subfamily", "genus", "species" };
        for ( int i = 1; i < broad_to_fine.length; i++ ) {
            final int prev = CladeLevel.rankDepth( broad_to_fine[ i - 1 ] );
            final int cur = CladeLevel.rankDepth( broad_to_fine[ i ] );
            if ( cur <= prev ) {
                return fail( broad_to_fine[ i ] + " must be FINER than " + broad_to_fine[ i - 1 ] + ", got depths "
                        + cur + " vs " + prev );
            }
        }
        if ( CladeLevel.rankDepth( "GENUS" ) != CladeLevel.rankDepth( "genus" ) ) {
            return fail( "rank depth must be case-insensitive" );
        }
        // "unknown"/"other" sit at the END of the underlying rank list but carry no real depth; treating them as
        // the finest ranks would wedge a rank of unknown breadth up against the tips
        for ( final String no_depth : new String[] { "unknown", "other", "not-a-rank", "", null } ) {
            if ( CladeLevel.rankDepth( no_depth ) != -1 ) {
                return fail( "\"" + no_depth + "\" has no meaningful depth and must report -1, got "
                        + CladeLevel.rankDepth( no_depth ) );
            }
        }
        return true;
    }

    private static boolean ordering() {
        // the headline behavior: whatever order the user picks, the FINEST rank is drawn first (nearest the tips)
        for ( final List<String> picked : Arrays.asList( Arrays.asList( "family", "genus", "order" ),
                                                         Arrays.asList( "order", "family", "genus" ),
                                                         Arrays.asList( "genus", "order", "family" ) ) ) {
            final List<String> got = ranks( CladeLevel.order( specs( picked ) ) );
            if ( !got.equals( Arrays.asList( "genus", "family", "order" ) ) ) {
                return fail( "picked " + picked + " must draw as [genus, family, order], got " + got );
            }
        }
        // blanks are dropped, and a repeated rank cannot claim two columns (it would draw the same marks twice)
        final List<CladeLevel.Spec> messy = new ArrayList<>();
        messy.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL ) );
        messy.add( null );
        messy.add( new CladeLevel.Spec( "", CLADE_LABEL_ANGLE.VERTICAL ) );
        messy.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL ) );
        messy.add( new CladeLevel.Spec( "FAMILY", CLADE_LABEL_ANGLE.VERTICAL ) ); // same rank, other case
        final List<String> cleaned = ranks( CladeLevel.order( messy ) );
        if ( !cleaned.equals( Arrays.asList( "genus", "family" ) ) ) {
            return fail( "blanks/nulls/duplicates must be dropped, got " + cleaned );
        }
        // never more columns than the layout can carry
        final List<String> many = Arrays.asList( "species", "genus", "subfamily", "family", "order", "class" );
        final List<CladeLevel.Spec> capped = CladeLevel.order( specs( many ) );
        if ( capped.size() != CladeLevel.MAX_LEVELS ) {
            return fail( "at most " + CladeLevel.MAX_LEVELS + " levels may survive, got " + capped.size() );
        }
        if ( !ranks( capped ).equals( Arrays.asList( "species", "genus", "subfamily" ) ) ) {
            return fail( "the cap must keep the FINEST ranks (they are the ones nearest the tips), got "
                    + ranks( capped ) );
        }
        // a rank of unknown depth goes OUTERMOST (the safe end) rather than being taken for the finest
        final List<String> with_unknown = ranks( CladeLevel.order( specs( Arrays.asList( "mystery", "genus" ) ) ) );
        if ( !with_unknown.equals( Arrays.asList( "genus", "mystery" ) ) ) {
            return fail( "an unrecognized rank must sort outermost, got " + with_unknown );
        }
        // equal-depth (here: both unknown) specs keep the order the user gave them -- a stable sort
        final List<String> stable = ranks( CladeLevel.order( specs( Arrays.asList( "zeta", "alpha" ) ) ) );
        if ( !stable.equals( Arrays.asList( "zeta", "alpha" ) ) ) {
            return fail( "equal-depth specs must keep the user's order, got " + stable );
        }
        // the per-level label angle survives ordering (it is the whole point of the Spec)
        final List<CladeLevel.Spec> angled = new ArrayList<>();
        angled.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.HORIZONTAL ) );
        angled.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL ) );
        final List<CladeLevel.Spec> ordered = CladeLevel.order( angled );
        if ( ( ordered.get( 0 ).getLabelAngle() != CLADE_LABEL_ANGLE.VERTICAL )
                || ( ordered.get( 1 ).getLabelAngle() != CLADE_LABEL_ANGLE.HORIZONTAL ) ) {
            return fail( "each level's label angle must travel with its rank through the re-ordering" );
        }
        if ( !CladeLevel.order( null ).isEmpty() ) {
            return fail( "a null spec list must yield no levels" );
        }
        return true;
    }

    private static List<CladeLevel.Spec> specs( final List<String> ranks ) {
        final List<CladeLevel.Spec> out = new ArrayList<>();
        for ( final String r : ranks ) {
            out.add( new CladeLevel.Spec( r, CLADE_LABEL_ANGLE.VERTICAL ) );
        }
        return out;
    }

    private static List<String> ranks( final List<CladeLevel.Spec> specs ) {
        final List<String> out = new ArrayList<>();
        for ( final CladeLevel.Spec s : specs ) {
            out.add( s.getRank() );
        }
        return out;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [CladeLevelTest] " + msg );
        return false;
    }

    private CladeLevelTest() {
    }
}
