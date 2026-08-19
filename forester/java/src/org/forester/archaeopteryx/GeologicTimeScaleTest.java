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

import java.util.List;

import org.forester.archaeopteryx.GeologicTimeScale.Interval;
import org.forester.archaeopteryx.GeologicTimeScale.Rank;

/**
 * Sanity checks on the embedded ICS geologic time scale: the four populated ranks (Eon / Era / Period / Epoch) have
 * the expected counts, their intervals TILE their range with no gaps or overlaps, boundaries are monotonic, the age
 * lookups / overlap queries return the right units (incl. the thematic ~150 Ma Late Jurassic where <em>Archaeopteryx</em>
 * lived), and {@link GeologicTimeScale#bandRanks} adapts the two-band pair to the tree's depth (Period/Epoch ->
 * Era/Period -> Eon/Era) so deep-time (Precambrian) trees are fully banded.
 */
public final class GeologicTimeScaleTest {

    public static boolean test() {
        try {
            return countsOk() && contiguousOk( Rank.EON ) && contiguousOk( Rank.ERA ) && contiguousOk( Rank.PERIOD )
                    && contiguousOk( Rank.EPOCH ) && lookupOk() && overlapOk() && bandRanksOk();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean countsOk() {
        ck( GeologicTimeScale.intervals( Rank.EON ).size() == 3, "3 eons" );
        ck( GeologicTimeScale.intervals( Rank.ERA ).size() == 10, "10 eras" );
        ck( GeologicTimeScale.intervals( Rank.PERIOD ).size() == 22, "22 periods" );
        ck( GeologicTimeScale.intervals( Rank.EPOCH ).size() == 34, "34 epochs" );
        return true;
    }

    /** The two-band axis picks a coarse+fine rank pair that adapts to the tree's depth, so BOTH bands always fully
     *  cover [0, oldest] -- Period/Epoch for a Phanerozoic tree, Era/Period into the Proterozoic (epochs run out at
     *  ~539 Ma), Eon/Era for a deep Archean tree (periods run out at ~2500 Ma). */
    private static boolean bandRanksOk() {
        ck( java.util.Arrays.equals( GeologicTimeScale.bandRanks( 150 ), new Rank[] { Rank.PERIOD, Rank.EPOCH } ),
            "a 150 Ma (Phanerozoic) tree bands Period over Epoch" );
        ck( java.util.Arrays.equals( GeologicTimeScale.bandRanks( 538.8 ), new Rank[] { Rank.PERIOD, Rank.EPOCH } ),
            "at exactly 538.8 Ma (base of the Phanerozoic) epochs still cover it -> Period over Epoch" );
        ck( java.util.Arrays.equals( GeologicTimeScale.bandRanks( 1000 ), new Rank[] { Rank.ERA, Rank.PERIOD } ),
            "a 1000 Ma (Proterozoic) tree bands Era over Period (epochs run out at ~539 Ma)" );
        ck( java.util.Arrays.equals( GeologicTimeScale.bandRanks( 3500 ), new Rank[] { Rank.EON, Rank.ERA } ),
            "a 3500 Ma (Archean) tree bands Eon over Era (periods run out at ~2500 Ma)" );
        // the finer band of each pair fully covers the depth (no blank deep segment)
        ck( GeologicTimeScale.coverageMa( Rank.EPOCH ) >= 538.0, "epochs cover to the base of the Phanerozoic" );
        ck( GeologicTimeScale.coverageMa( Rank.PERIOD ) >= 2500.0, "periods cover into the Paleoproterozoic" );
        ck( GeologicTimeScale.coverageMa( Rank.ERA ) >= 4031.0, "eras cover the whole Archean" );
        // deep-time lookups
        ck( "Archean".equals( named( GeologicTimeScale.at( Rank.EON, 3000 ) ) ), "3000 Ma is in the Archean eon" );
        ck( "Mesoarchean".equals( named( GeologicTimeScale.at( Rank.ERA, 3000 ) ) ), "3000 Ma is the Mesoarchean era" );
        ck( "Proterozoic".equals( named( GeologicTimeScale.at( Rank.EON, 1000 ) ) ), "1000 Ma is the Proterozoic eon" );
        return true;
    }

    /** Within a rank, intervals are ordered youngest-first, each has young < old, and each interval's OLD boundary
     *  equals the next interval's YOUNG boundary (contiguous, no gaps/overlaps). */
    private static boolean contiguousOk( final Rank rank ) {
        final List<Interval> ivs = GeologicTimeScale.intervals( rank );
        for ( int i = 0; i < ivs.size(); ++i ) {
            final Interval iv = ivs.get( i );
            ck( iv.youngMa() < iv.oldMa(), rank + " " + iv.name() + " must have young < old" );
            if ( i > 0 ) {
                ck( Math.abs( ivs.get( i - 1 ).oldMa() - iv.youngMa() ) < 1e-6,
                    rank + ": " + ivs.get( i - 1 ).name() + " (old " + ivs.get( i - 1 ).oldMa() + ") must abut "
                            + iv.name() + " (young " + iv.youngMa() + ")" );
            }
        }
        ck( ivs.get( 0 ).youngMa() == 0.0, rank + " must start at the present (0 Ma)" );
        return true;
    }

    private static boolean lookupOk() {
        // ~150 Ma: Jurassic period, Late Jurassic epoch -- when Archaeopteryx lived
        ck( "Jurassic".equals( named( GeologicTimeScale.at( Rank.PERIOD, 150 ) ) ), "150 Ma is in the Jurassic" );
        ck( "Late Jurassic".equals( named( GeologicTimeScale.at( Rank.EPOCH, 150 ) ) ),
            "150 Ma is in the Late Jurassic (Archaeopteryx)" );
        ck( "Cretaceous".equals( named( GeologicTimeScale.at( Rank.PERIOD, 66 ) ) ),
            "66 Ma (the K-Pg boundary) is the young edge of the Cretaceous" );
        ck( "Quaternary".equals( named( GeologicTimeScale.at( Rank.PERIOD, 0 ) ) ), "0 Ma is the Quaternary" );
        ck( GeologicTimeScale.at( Rank.PERIOD, 99999 ) == null, "an out-of-range age has no interval" );
        return true;
    }

    private static boolean overlapOk() {
        // a "dinosaur" window 66-201 Ma spans Cretaceous + Jurassic + (just) Triassic at period rank
        final List<Interval> periods = GeologicTimeScale.overlapping( Rank.PERIOD, 66, 200 );
        ck( containsName( periods, "Cretaceous" ) && containsName( periods, "Jurassic" ),
            "a 66-200 Ma window overlaps the Cretaceous and Jurassic periods" );
        ck( !containsName( periods, "Paleogene" ), "a 66-200 Ma window does NOT overlap the (younger) Paleogene" );
        // at epoch rank, a 66-100 Ma window is just the Late Cretaceous
        final List<Interval> epochs = GeologicTimeScale.overlapping( Rank.EPOCH, 66, 100 );
        ck( containsName( epochs, "Late Cretaceous" ), "a 66-100 Ma window overlaps the Late Cretaceous epoch" );
        return true;
    }

    private static boolean containsName( final List<Interval> ivs, final String name ) {
        for ( final Interval iv : ivs ) {
            if ( name.equals( iv.name() ) ) {
                return true;
            }
        }
        return false;
    }

    private static String named( final Interval iv ) {
        return ( iv == null ) ? null : iv.name();
    }

    private static void ck( final boolean cond, final String msg ) {
        if ( !cond ) {
            throw new RuntimeException( "GeologicTimeScaleTest: " + msg );
        }
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private GeologicTimeScaleTest() {
    }
}
