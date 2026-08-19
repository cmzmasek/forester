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
 * Populated ranks: EON, ERA, PERIOD, EPOCH -- so the two-band axis can adapt its rank pair to the tree's depth
 * (Period/Epoch for a Phanerozoic tree, Era/Period into the Proterozoic, Eon/Era for a deep Archean tree; see
 * {@link #bandRanks(double)}). STAGE/AGE (the finer band for very recent trees) is a follow-on. Ages/colours follow
 * the ICS scheme (via Macrostrat) and should be spot-checked against the official chart PDF before shipping.
 */
final class GeologicTimeScale {

    /** ICS reference for the citation dialogs. */
    static final String REFERENCE =
            "Cohen, K.M., Harper, D.A.T., Gibbard, P.L. & Car, N. (2025, updated) "
                    + "The ICS International Chronostratigraphic Chart this decade. Episodes 48: 105-115.";
    static final String SOURCE = "International Commission on Stratigraphy (ICS / IUGS), www.stratigraphy.org";

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

        double durationMa() {
            return _old_ma - _young_ma;
        }

        Color color() {
            return _color;
        }
    }

    private static final List<Interval> EONS    = new ArrayList<>();
    private static final List<Interval> ERAS    = new ArrayList<>();
    private static final List<Interval> PERIODS = new ArrayList<>();
    private static final List<Interval> EPOCHS  = new ArrayList<>();

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

    /** All intervals of a rank (unmodifiable-ish; callers must not mutate), youngest first. Empty for AGE (a follow-on). */
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

    /** The coarse+fine rank pair for the two-band axis, adapted to a tree spanning {@code [0, old_ma]}: Period over
     *  Epoch for a Phanerozoic tree (epochs cover it), Era over Period once it reaches into the Proterozoic (epochs
     *  run out at ~539 Ma), and Eon over Era for a deep Archean tree (periods run out at ~2500 Ma) -- so BOTH bands
     *  always fully cover the range (no blank deep segment). Returns {@code [upper (coarser), lower (finer)]}. */
    static Rank[] bandRanks( final double old_ma ) {
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

    /** The oldest boundary we have data for (Ma) -- the deep end of the embedded scale (the Eon/Era ranks reach the
     *  base of the Archean, ~4031 Ma). */
    static double oldestMa() {
        return Math.max( coverageMa( Rank.EON ), Math.max( coverageMa( Rank.ERA ),
                Math.max( coverageMa( Rank.PERIOD ), coverageMa( Rank.EPOCH ) ) ) );
    }

    private GeologicTimeScale() {
    }
}
