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
import java.util.List;

/**
 * The International Chronostratigraphic Chart (the official geologic time scale) as embedded data, for drawing a
 * colored geologic time axis under a dated ("time tree") phylogeny. Each interval carries its rank, its boundary
 * ages in millions of years before present (Ma), and its official ICS colour.
 * <p>
 * Data source / reference (put in Help &gt; References and the README):
 * Cohen, K.M., Harper, D.A.T., Gibbard, P.L. &amp; Car, N. (2025, updated) "The ICS International Chronostratigraphic
 * Chart this decade." Episodes 48: 105-115. Boundary ages and unit colours from the International Commission on
 * Stratigraphy (ICS / IUGS), www.stratigraphy.org.
 * <p>
 * Populated ranks: EON, ERA, PERIOD, EPOCH and AGE (Stage) -- so the two-band axis can adapt its rank pair to the
 * window the tree actually spans: Epoch/Stage for a tree inside one or two Series, Period/Epoch for a Phanerozoic
 * tree, Era/Period into the Proterozoic, Eon/Era for a deep Archean tree (see {@link #bandRanks(double, double)}).
 * Ages/colours follow the ICS scheme (via Macrostrat); every value here was checked against that source, and the
 * tests pin the coverage and contiguity of each rank.
 */
final class GeologicTimeScale {

    /** ICS reference for the citation dialogs. */
    static final String REFERENCE =
            "Cohen, K.M., Harper, D.A.T., Gibbard, P.L. & Car, N. (2025, updated) "
                    + "The ICS International Chronostratigraphic Chart this decade. Episodes 48: 105-115.";

    /** The chronostratigraphic ranks, coarse -> fine. EON, ERA, PERIOD, EPOCH are populated (AGE is a follow-on). */
    enum Rank {
        EON, ERA, PERIOD, EPOCH, AGE
    }

    /** One named time interval: [young, old] in Ma (young &lt; old), plus its rank and official colour. */
    static final class Interval {

        private final String _name;
        private final Rank   _rank;
        private final double _young_ma; // younger (smaller-Ma) boundary
        private final double _old_ma;   // older (larger-Ma) boundary
        private final Color  _color;

        Interval( final String name, final Rank rank, final double young_ma, final double old_ma, final Color color ) {
            _name = name;
            _rank = rank;
            _young_ma = young_ma;
            _old_ma = old_ma;
            _color = color;
        }

        String name() {
            return _name;
        }

        Rank rank() {
            return _rank;
        }

        double youngMa() {
            return _young_ma;
        }

        double oldMa() {
            return _old_ma;
        }

        Color color() {
            return _color;
        }
    }

    private static final List<Interval> EONS    = new ArrayList<>();
    private static final List<Interval> ERAS    = new ArrayList<>();
    private static final List<Interval> PERIODS = new ArrayList<>();
    private static final List<Interval> EPOCHS  = new ArrayList<>();
    private static final List<Interval> AGES    = new ArrayList<>();

    static {
        // ---- Eons (Eonothem/Eon), youngest -> oldest ----
        eon( "Phanerozoic", 0, 538.8, 0x9AD9DD );
        eon( "Proterozoic", 538.8, 2500, 0xFF70B8 );
        eon( "Archean", 2500, 4031, 0xFF3399 );
        // ---- Eras (Erathem/Era), youngest -> oldest ----
        era( "Cenozoic", 0, 66, 0xF2F91D );
        era( "Mesozoic", 66, 251.902, 0x67C5CA );
        era( "Paleozoic", 251.902, 538.8, 0x99C08D );
        era( "Neoproterozoic", 538.8, 1000, 0xFF9BCD );
        era( "Mesoproterozoic", 1000, 1600, 0xFF7EBF );
        era( "Paleoproterozoic", 1600, 2500, 0xE665A6 );
        era( "Neoarchean", 2500, 2800, 0xFF5CAD );
        era( "Mesoarchean", 2800, 3200, 0xE62E8A );
        era( "Paleoarchean", 3200, 3600, 0xCC297A );
        era( "Eoarchean", 3600, 4031, 0xB2246B );
        // ---- Periods (System/Period), youngest -> oldest ----
        p( "Quaternary", 0, 2.58, 0xF9F97F );
        p( "Neogene", 2.58, 23.04, 0xFFE619 );
        p( "Paleogene", 23.04, 66, 0xFD9A52 );
        p( "Cretaceous", 66, 143.1, 0x7FC64E );
        p( "Jurassic", 143.1, 201.4, 0x34B2C9 );
        p( "Triassic", 201.4, 251.902, 0x812B92 );
        p( "Permian", 251.902, 298.9, 0xF04028 );
        p( "Carboniferous", 298.9, 358.86, 0x67A599 );
        p( "Devonian", 358.86, 419.62, 0xCB8C37 );
        p( "Silurian", 419.62, 443.1, 0xB3E1B6 );
        p( "Ordovician", 443.1, 486.85, 0x009270 );
        p( "Cambrian", 486.85, 538.8, 0x7FA056 );
        p( "Ediacaran", 538.8, 635, 0xFFC3E1 );
        p( "Cryogenian", 635, 720, 0xFFAFD7 );
        p( "Tonian", 720, 1000, 0xFFA5D2 );
        p( "Stenian", 1000, 1200, 0xFFA5D2 );
        p( "Ectasian", 1200, 1400, 0xFF98CC );
        p( "Calymmian", 1400, 1600, 0xFF8BC5 );
        p( "Statherian", 1600, 1800, 0xEE93C1 );
        p( "Orosirian", 1800, 2050, 0xE874AF );
        p( "Rhyacian", 2050, 2300, 0xEB84B8 );
        p( "Siderian", 2300, 2500, 0xE874AF );
        // ---- Epochs (Series/Epoch), youngest -> oldest ----
        e( "Holocene", 0, 0.0117, 0xFEF2E0 );
        e( "Pleistocene", 0.0117, 2.58, 0xFFF2AE );
        e( "Pliocene", 2.58, 5.333, 0xFFFF99 );
        e( "Miocene", 5.333, 23.04, 0xFFFF00 );
        e( "Oligocene", 23.04, 33.9, 0xFDC07A );
        e( "Eocene", 33.9, 56, 0xFDB46C );
        e( "Paleocene", 56, 66, 0xFDA75F );
        e( "Late Cretaceous", 66, 100.5, 0xA6D84A );
        e( "Early Cretaceous", 100.5, 143.1, 0x8CCD57 );
        e( "Late Jurassic", 143.1, 161.5, 0xB3E3EE );
        e( "Middle Jurassic", 161.5, 174.7, 0x80CFD8 );
        e( "Early Jurassic", 174.7, 201.4, 0x42AED0 );
        e( "Late Triassic", 201.4, 237, 0xBD8CC3 );
        e( "Middle Triassic", 237, 246.7, 0xB168B1 );
        e( "Early Triassic", 246.7, 251.902, 0x983999 );
        e( "Lopingian", 251.902, 259.51, 0xFBA794 );
        e( "Guadalupian", 259.51, 274.4, 0xFB745C );
        e( "Cisuralian", 274.4, 298.9, 0xEF5845 );
        e( "Pennsylvanian", 298.9, 323.4, 0x99C2B5 );
        e( "Mississippian", 323.4, 358.86, 0x678F66 );
        e( "Late Devonian", 358.86, 382.31, 0xF1E19D );
        e( "Middle Devonian", 382.31, 393.47, 0xF1C868 );
        e( "Early Devonian", 393.47, 419.62, 0xE5AC4D );
        e( "Pridoli", 419.62, 422.7, 0xE6F5E1 );
        e( "Ludlow", 422.7, 426.7, 0xBFE6CF );
        e( "Wenlock", 426.7, 432.9, 0xB3E1C2 );
        e( "Llandovery", 432.9, 443.1, 0x99D7B3 );
        e( "Late Ordovician", 443.1, 458.2, 0x7FCA93 );
        e( "Middle Ordovician", 458.2, 471.3, 0x4DB47E );
        e( "Early Ordovician", 471.3, 486.85, 0x1A9D6F );
        e( "Furongian", 486.85, 497, 0xB3E095 );
        e( "Miaolingian", 497, 506.5, 0xA6CF86 );
        e( "Series 2", 506.5, 521, 0x99C078 );
        e( "Terreneuvian", 521, 538.8, 0x8CB06C );
        // ---- Ages (Stage/Age), youngest -> oldest: the finest ICS rank, and the finer band for a tree that
        // spans only one or two Series. They exist for the Phanerozoic only (0 - 538.8 Ma) -- the Precambrian
        // has no ratified stages -- which is exactly the depth at which bandRanks can reach for them. The
        // PRIDOLI is a Series with no stages of its own, so (as on the ICS chart itself) it stands in the stage
        // row for its own span; without it the band would have a hole at 419.62-422.7 Ma.
        a( "Meghalayan", 0, 0.0042, 0xFEF2E0 );
        a( "Northgrippian", 0.0042, 0.0082, 0xFEF2E0 );
        a( "Greenlandian", 0.0082, 0.0117, 0xFEF2E0 );
        a( "Late Pleistocene", 0.0117, 0.129, 0xFFF2C7 );
        a( "Chibanian", 0.129, 0.774, 0xFFF2C7 );
        a( "Calabrian", 0.774, 1.8, 0xFFF2C7 );
        a( "Gelasian", 1.8, 2.58, 0xFFEDB3 );
        a( "Piacenzian", 2.58, 3.6, 0xFFFFBF );
        a( "Zanclean", 3.6, 5.333, 0xFFFFB3 );
        a( "Messinian", 5.333, 7.246, 0xFFFF73 );
        a( "Tortonian", 7.246, 11.63, 0xFFFF66 );
        a( "Serravallian", 11.63, 13.82, 0xFFFF59 );
        a( "Langhian", 13.82, 15.98, 0xFFFF4D );
        a( "Burdigalian", 15.98, 20.45, 0xFFFF41 );
        a( "Aquitanian", 20.45, 23.04, 0xFFFF33 );
        a( "Chattian", 23.04, 27.3, 0xFEE6AA );
        a( "Rupelian", 27.3, 33.9, 0xFED99A );
        a( "Priabonian", 33.9, 37.71, 0xFDCDA1 );
        a( "Bartonian", 37.71, 41.03, 0xFDC091 );
        a( "Lutetian", 41.03, 48.07, 0xFCB482 );
        a( "Ypresian", 48.07, 56, 0xFCA773 );
        a( "Thanetian", 56, 59.24, 0xFDBF6F );
        a( "Selandian", 59.24, 61.66, 0xFEBF65 );
        a( "Danian", 61.66, 66, 0xFDB462 );
        a( "Maastrichtian", 66, 72.2, 0xF2FA8C );
        a( "Campanian", 72.2, 83.6, 0xE6F47F );
        a( "Santonian", 83.6, 85.7, 0xD9EF74 );
        a( "Coniacian", 85.7, 89.8, 0xCCE968 );
        a( "Turonian", 89.8, 93.9, 0xBFE35D );
        a( "Cenomanian", 93.9, 100.5, 0xB3DE53 );
        a( "Albian", 100.5, 113.2, 0xCCEA97 );
        a( "Aptian", 113.2, 121.4, 0xBFE48A );
        a( "Barremian", 121.4, 125.77, 0xB3DF7F );
        a( "Hauterivian", 125.77, 132.6, 0xA6D975 );
        a( "Valanginian", 132.6, 137.05, 0x99D36A );
        a( "Berriasian", 137.05, 143.1, 0x8CCD60 );
        a( "Tithonian", 143.1, 149.2, 0xD9F1F7 );
        a( "Kimmeridgian", 149.2, 154.8, 0xCCECF4 );
        a( "Oxfordian", 154.8, 161.5, 0xBFE7F1 );
        a( "Callovian", 161.5, 165.3, 0xBFE7E5 );
        a( "Bathonian", 165.3, 168.2, 0xB3E2E3 );
        a( "Bajocian", 168.2, 170.9, 0xA6DDE0 );
        a( "Aalenian", 170.9, 174.7, 0x9AD9DD );
        a( "Toarcian", 174.7, 184.2, 0x99CEE3 );
        a( "Pliensbachian", 184.2, 192.9, 0x80C5DD );
        a( "Sinemurian", 192.9, 199.5, 0x67BCD8 );
        a( "Hettangian", 199.5, 201.4, 0x4EB3D3 );
        a( "Rhaetian", 201.4, 205.7, 0xE3B9DB );
        a( "Norian", 205.7, 227.3, 0xD6AAD3 );
        a( "Carnian", 227.3, 237, 0xC99BCB );
        a( "Ladinian", 237, 241.464, 0xC983BF );
        a( "Anisian", 241.464, 246.7, 0xBC75B7 );
        a( "Olenekian", 246.7, 249.9, 0xB051A5 );
        a( "Induan", 249.9, 251.902, 0xA4469F );
        a( "Changhsingian", 251.902, 254.14, 0xFCC0B2 );
        a( "Wuchiapingian", 254.14, 259.51, 0xFCB4A2 );
        a( "Capitanian", 259.51, 264.28, 0xFB9A85 );
        a( "Wordian", 264.28, 266.9, 0xFB8D76 );
        a( "Roadian", 266.9, 274.4, 0xFB8069 );
        a( "Kungurian", 274.4, 283.3, 0xE38776 );
        a( "Artinskian", 283.3, 290.1, 0xE37B68 );
        a( "Sakmarian", 290.1, 293.52, 0xE36F5C );
        a( "Asselian", 293.52, 298.9, 0xE36350 );
        a( "Gzhelian", 298.9, 303.7, 0xCCD4C7 );
        a( "Kasimovian", 303.7, 307, 0xBFD0C5 );
        a( "Moscovian", 307, 315.2, 0xC7CBB9 );
        a( "Bashkirian", 315.2, 323.4, 0x99C2B5 );
        a( "Serpukhovian", 323.4, 330.3, 0xBFC26B );
        a( "Visean", 330.3, 346.7, 0xA6B96C );
        a( "Tournaisian", 346.7, 358.86, 0x8CB06C );
        a( "Famennian", 358.86, 372.15, 0xF2EDC5 );
        a( "Frasnian", 372.15, 382.31, 0xF2EDAD );
        a( "Givetian", 382.31, 387.95, 0xF1E185 );
        a( "Eifelian", 387.95, 393.47, 0xF1D576 );
        a( "Emsian", 393.47, 410.62, 0xE5D075 );
        a( "Pragian", 410.62, 413.02, 0xE5C468 );
        a( "Lochkovian", 413.02, 419.62, 0xE5B75A );
        a( "Pridoli", 419.62, 422.7, 0xE6F5E1 );
        a( "Ludfordian", 422.7, 425, 0xD9F0DF );
        a( "Gorstian", 425, 426.7, 0xCCECDD );
        a( "Homerian", 426.7, 430.6, 0xCCEBD1 );
        a( "Sheinwoodian", 430.6, 432.9, 0xBFE6C3 );
        a( "Telychian", 432.9, 438.6, 0xBFE6CF );
        a( "Aeronian", 438.6, 440.5, 0xB3E1C2 );
        a( "Rhuddanian", 440.5, 443.1, 0xA6DCB5 );
        a( "Hirnantian", 443.1, 445.2, 0xA6DBAB );
        a( "Katian", 445.2, 452.8, 0x99D69F );
        a( "Sandbian", 452.8, 458.2, 0x8CD094 );
        a( "Darriwilian", 458.2, 469.4, 0x74C69C );
        a( "Dapingian", 469.4, 471.3, 0x66C092 );
        a( "Floian", 471.3, 477.1, 0x41B087 );
        a( "Tremadocian", 477.1, 486.85, 0x33A97E );
        a( "Stage 10", 486.85, 491, 0xE6F5C9 );
        a( "Jiangshanian", 491, 494.2, 0xD9F0BB );
        a( "Paibian", 494.2, 497, 0xCCEBAE );
        a( "Guzhangian", 497, 500.5, 0xCCDFAA );
        a( "Drumian", 500.5, 504.5, 0xBFD99D );
        a( "Wuliuan", 504.5, 506.5, 0xB3D492 );
        a( "Stage 4", 506.5, 514.5, 0xB3CA8E );
        a( "Stage 3", 514.5, 521, 0xA6C583 );
        a( "Stage 2", 521, 529, 0xA6BA80 );
        a( "Fortunian", 529, 538.8, 0x99B575 );
    }

    private static void eon( final String name, final double young, final double old, final int rgb ) {
        EONS.add( new Interval( name, Rank.EON, young, old, new Color( rgb ) ) );
    }

    private static void era( final String name, final double young, final double old, final int rgb ) {
        ERAS.add( new Interval( name, Rank.ERA, young, old, new Color( rgb ) ) );
    }

    private static void p( final String name, final double young, final double old, final int rgb ) {
        PERIODS.add( new Interval( name, Rank.PERIOD, young, old, new Color( rgb ) ) );
    }

    private static void e( final String name, final double young, final double old, final int rgb ) {
        EPOCHS.add( new Interval( name, Rank.EPOCH, young, old, new Color( rgb ) ) );
    }

    private static void a( final String name, final double young, final double old, final int rgb ) {
        AGES.add( new Interval( name, Rank.AGE, young, old, new Color( rgb ) ) );
    }

    /** All intervals of a rank (unmodifiable-ish; callers must not mutate), youngest first. */
    static List<Interval> intervals( final Rank rank ) {
        switch ( rank ) {
            case EON:
                return EONS;
            case ERA:
                return ERAS;
            case PERIOD:
                return PERIODS;
            case EPOCH:
                return EPOCHS;
            case AGE:
                return AGES;
            default:
                return new ArrayList<>();
        }
    }

    /** The oldest boundary (Ma) covered by a rank's data -- i.e. how deep that rank can band. */
    static double coverageMa( final Rank rank ) {
        double m = 0;
        for ( final Interval iv : intervals( rank ) ) {
            m = Math.max( m, iv.oldMa() );
        }
        return m;
    }

    /** A window this narrow -- one or two Series -- is banded Epoch over Stage rather than Period over Epoch. */
    static final int MAX_SERIES_FOR_STAGE_BANDS = 2;

    /** The pair for a tree measured from the present; see {@link #bandRanks(double, double)}. */
    static Rank[] bandRanks( final double old_ma ) {
        return bandRanks( 0, old_ma );
    }

    /** The coarse+fine rank pair for the two-band axis, adapted to a tree spanning {@code [young_ma, old_ma]}:
     *  Epoch over Stage for a window narrow enough to sit in one or two Series (else the axis would be two huge
     *  blocks and no scale at all), Period over Epoch for a Phanerozoic tree, Era over Period once it reaches into
     *  the Proterozoic (epochs run out at ~539 Ma), and Eon over Era for a deep Archean tree (periods run out at
     *  ~2500 Ma) -- so BOTH bands always fully cover the range (no blank segment at either end). Returns
     *  {@code [upper (coarser), lower (finer)]}. */
    static Rank[] bandRanks( final double young_ma, final double old_ma ) {
        if ( old_ma <= coverageMa( Rank.AGE ) ) { // stages exist for the Phanerozoic only
            final int series = overlapping( Rank.EPOCH, young_ma, old_ma ).size();
            if ( ( series > 0 ) && ( series <= MAX_SERIES_FOR_STAGE_BANDS ) ) {
                return new Rank[] { Rank.EPOCH, Rank.AGE };
            }
        }
        if ( old_ma <= coverageMa( Rank.EPOCH ) ) {
            return new Rank[] { Rank.PERIOD, Rank.EPOCH };
        }
        if ( old_ma <= coverageMa( Rank.PERIOD ) ) {
            return new Rank[] { Rank.ERA, Rank.PERIOD };
        }
        return new Rank[] { Rank.EON, Rank.ERA };
    }

    /** The intervals of {@code rank} that overlap the age window [young_ma, old_ma] (Ma), youngest first -- i.e. the
     *  bands to draw for a tree spanning that window. A zero/negative window yields the single interval at young_ma. */
    static List<Interval> overlapping( final Rank rank, final double young_ma, final double old_ma ) {
        final double lo = Math.min( young_ma, old_ma );
        final double hi = Math.max( young_ma, old_ma );
        final List<Interval> out = new ArrayList<>();
        for ( final Interval iv : intervals( rank ) ) {
            // overlap of [iv.young, iv.old] with [lo, hi], treating a zero-width window as a point query
            if ( ( iv.oldMa() > lo ) && ( iv.youngMa() < hi ) ) {
                out.add( iv );
            }
            else if ( ( lo == hi ) && ( iv.youngMa() <= lo ) && ( lo < iv.oldMa() ) ) {
                out.add( iv );
            }
        }
        return out;
    }

    /** The interval of {@code rank} containing {@code age_ma} (young &lt;= age &lt; old), or null if out of range. */
    static Interval at( final Rank rank, final double age_ma ) {
        for ( final Interval iv : intervals( rank ) ) {
            if ( ( iv.youngMa() <= age_ma ) && ( age_ma < iv.oldMa() ) ) {
                return iv;
            }
        }
        return null;
    }

    private GeologicTimeScale() {
    }
}
