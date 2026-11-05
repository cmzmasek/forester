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
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.forester.util.ForesterUtil;
import org.forester.util.TaxonomyUtil;

/**
 * One RANK's worth of clade annotation -- the rank, how its marks are drawn and how its taxon labels should read,
 * and the {@link CladeBand}s computed for it. "Annotate Clades by Rank" draws up to {@link #MAX_LEVELS} of these
 * at once: at most ONE level of shaded boxes (a translucent wash behind each clade; of several box ranks the
 * finest that places anything is the one drawn), and the others as nested columns of bars or brackets past the
 * tip labels (e.g. genus inside family inside order), which is how the nested rank brackets of a published
 * phylogeny figure are built -- or, with a box level, the shaded clade with the bracket that names it (Christian,
 * 2026-11-04, "go, box + two").
 * <p>
 * Split, like {@code AnnotationColumns}, into the user's {@link Spec} (which rank, which style, which label angle
 * -- kept across navigation) and the bands built from it for the current view (recomputed whenever the displayed
 * tree changes).
 * <p>
 * <b>Placement order is derived, never chosen.</b> {@link #order} sorts the levels FINEST-RANK-FIRST, so the
 * finest rank is drawn nearest the tips and the broadest outermost -- the convention of the published figures
 * this feature exists to reproduce, and the only order in which the nesting reads correctly (a broad clade's
 * mark must span the several finer marks inside it). The user therefore cannot pick an order that draws wrong.
 * The box level is not a column: it sits behind the tree whatever its rank, and the paint draws it first.
 */
final class CladeLevel {

    /** The most levels that can be annotated at once, the box level included. Three fits with compact labels;
     *  more does not read. */
    static final int MAX_LEVELS = 3;

    /**
     * The user's choice for one level: the taxonomic rank, the style of its marks (shaded boxes, bars or
     * brackets; null = whatever the tab's mode says, resolved by {@link #order}), and the on-screen angle of that
     * level's taxon labels. The angle is per level on purpose -- an outermost rank has few, long names worth
     * reading straight, while inner levels must stay narrow or three columns will not fit beside the tree.
     */
    static final class Spec {

        private final String                    _rank;
        private final TreePanel.CLADE_LABEL_ANGLE _label_angle;
        private final TreePanel.CLADE_VIS         _style;

        /** A spec whose style the tab's mode decides (every caller before styles were per level). */
        Spec( final String rank, final TreePanel.CLADE_LABEL_ANGLE label_angle ) {
            this( rank, label_angle, null );
        }

        Spec( final String rank, final TreePanel.CLADE_LABEL_ANGLE label_angle, final TreePanel.CLADE_VIS style ) {
            _rank = rank;
            _label_angle = ( label_angle == null ) ? TreePanel.CLADE_LABEL_ANGLE.VERTICAL : label_angle;
            _style = style;
        }

        String getRank() {
            return _rank;
        }

        TreePanel.CLADE_LABEL_ANGLE getLabelAngle() {
            return _label_angle;
        }

        /** How this level's marks are drawn; null until {@link #order} resolved it against the tab's mode. */
        TreePanel.CLADE_VIS getStyle() {
            return _style;
        }

        /** Whether this level is the shaded boxes (else a bar or bracket column). */
        boolean isBox() {
            return _style == TreePanel.CLADE_VIS.BOXES;
        }

        /** This spec with {@code style} where it has none. */
        Spec resolved( final TreePanel.CLADE_VIS style ) {
            return ( _style != null ) ? this : new Spec( _rank, _label_angle, style );
        }

        @Override
        public String toString() {
            return _rank + "/" + _label_angle + "/" + _style;
        }
    }

    private final Spec            _spec;
    private final List<CladeBand> _bands;

    /** @throws IllegalArgumentException for a spec without a style: every level comes through {@link #order}, which
     *          resolves it, and the paints switch on it */
    CladeLevel( final Spec spec, final List<CladeBand> bands ) {
        if ( ( spec == null ) || ( spec.getStyle() == null ) ) {
            throw new IllegalArgumentException( "a clade level's spec must carry a style (resolve it through CladeLevel.order)" );
        }
        _spec = spec;
        _bands = ( bands == null ) ? Collections.<CladeBand> emptyList() : bands;
    }

    String getRank() {
        return _spec.getRank();
    }

    TreePanel.CLADE_LABEL_ANGLE getLabelAngle() {
        return _spec.getLabelAngle();
    }

    /** How this level's marks are drawn (never null: the specs a level is built from are resolved). */
    TreePanel.CLADE_VIS getStyle() {
        return _spec.getStyle();
    }

    /** Whether this level is the shaded boxes (else a bar or bracket column). */
    boolean isBox() {
        return _spec.isBox();
    }

    List<CladeBand> getBands() {
        return _bands;
    }

    /**
     * How deep a rank sits in the taxonomic hierarchy: bigger = finer (species &gt; genus &gt; family). Reads
     * {@link TaxonomyUtil#RANK_TO_INT}, whose list runs broad to fine.
     * <p>
     * An unrecognized rank -- and the placeholder "unknown"/"other", which carry no real depth despite sitting
     * at the end of that list -- gets {@code -1}, which {@link #order} places OUTERMOST. That is the safe end: a
     * rank of unknown breadth wedged up against the tips would claim to be the finest thing on the figure.
     */
    static int rankDepth( final String rank ) {
        if ( ForesterUtil.isEmpty( rank ) || TaxonomyUtil.UNKNOWN.equalsIgnoreCase( rank )
                || TaxonomyUtil.OTHER.equalsIgnoreCase( rank ) ) {
            return -1;
        }
        final Integer d = TaxonomyUtil.RANK_TO_INT.get( rank.toLowerCase() );
        return ( d == null ) ? -1 : d.intValue();
    }

    /**
     * The specs in DRAWING order -- finest rank first, i.e. nearest the tips -- every style resolved
     * ({@code default_style} where a spec names none), with blanks dropped and a rank kept once as a box and once
     * as a column (a column at the box level's own rank is the bracket that names the shaded clade). Every box
     * rank is kept: only the finest that places anything is DRAWN (two translucent washes multiply, so a clade
     * inside another would paint darker than either and no longer the colour its legend row claims), decided
     * where the bands are known ({@code TreePanel.drawableCladeLevels}), so a figure naming a box rank the tree
     * cannot place still draws its next one -- as several ranks in the old BOXES mode always did. The COLUMNS are
     * capped so the drawn levels never exceed {@link #MAX_LEVELS}, the broadest column going first and the box
     * never. A stable sort, so two ranks of equal (or unknown) depth keep the order the user gave them; a box
     * comes before a column at its own rank. Pure; the single source of the placement order for the rectangular
     * columns, the circular rings and every space reservation derived from them.
     */
    static List<Spec> order( final List<Spec> specs, final TreePanel.CLADE_VIS default_style ) {
        final List<Spec> out = new ArrayList<>();
        if ( specs == null ) {
            return out;
        }
        final TreePanel.CLADE_VIS fallback = ( default_style == null ) ? TreePanel.CLADE_VIS.BOXES : default_style;
        final Set<String> seen = new LinkedHashSet<>();
        for ( final Spec raw : specs ) {
            if ( ( raw == null ) || ForesterUtil.isEmpty( raw.getRank() ) ) {
                continue;
            }
            final Spec s = raw.resolved( fallback );
            // a repeated rank cannot claim two columns (it would draw the same marks twice); a box and a column at
            // one rank are two different marks of the same clades
            if ( seen.add( s.getRank().toLowerCase() + ( s.isBox() ? "|box" : "|column" ) ) ) {
                out.add( s );
            }
        }
        // stable: List.sort is guaranteed stable, so equal-depth specs keep the user's relative order -- except that
        // the box goes before a column at its own rank (one order for the figure text, whichever was picked first)
        out.sort( Comparator.comparingInt( ( final Spec s ) -> rankDepth( s.getRank() ) ).reversed()
                .thenComparingInt( ( final Spec s ) -> s.isBox() ? 0 : 1 ) );
        // the cap: the box level (one is drawn) plus at most MAX_LEVELS - 1 columns, the broadest columns dropped
        boolean has_box = false;
        int columns = 0;
        for ( final Spec s : out ) {
            if ( s.isBox() ) {
                has_box = true;
            }
            else {
                ++columns;
            }
        }
        final int max_columns = MAX_LEVELS - ( has_box ? 1 : 0 );
        for ( int i = out.size() - 1; ( i >= 0 ) && ( columns > max_columns ); --i ) {
            if ( !out.get( i ).isBox() ) {
                out.remove( i );
                --columns;
            }
        }
        return out;
    }
}
