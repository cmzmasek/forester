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
 * One RANK's worth of clade annotation -- the rank, how its taxon labels should read, and the
 * {@link CladeBand}s computed for it. "Annotate Clades by Rank" draws up to {@link #MAX_LEVELS} of these at
 * once as nested columns of bars or brackets (e.g. genus inside family inside order), which is how the nested
 * rank brackets of a published phylogeny figure are built.
 * <p>
 * Split, like {@code AnnotationColumns}, into the user's {@link Spec} (which rank, which label angle -- kept
 * across navigation) and the bands built from it for the current view (recomputed whenever the displayed tree
 * changes).
 * <p>
 * <b>Placement order is derived, never chosen.</b> {@link #order} sorts the levels FINEST-RANK-FIRST, so the
 * finest rank is drawn nearest the tips and the broadest outermost -- the convention of the published figures
 * this feature exists to reproduce, and the only order in which the nesting reads correctly (a broad clade's
 * mark must span the several finer marks inside it). The user therefore cannot pick an order that draws wrong.
 */
final class CladeLevel {

    /** The most levels that can be annotated at once. Three fits with compact labels; more does not read. */
    static final int MAX_LEVELS = 3;

    /**
     * The user's choice for one level: the taxonomic rank, and the on-screen angle of that level's taxon
     * labels. The angle is per level on purpose -- an outermost rank has few, long names worth reading
     * straight, while inner levels must stay narrow or three columns will not fit beside the tree.
     */
    static final class Spec {

        private final String                    _rank;
        private final TreePanel.CLADE_LABEL_ANGLE _label_angle;

        Spec( final String rank, final TreePanel.CLADE_LABEL_ANGLE label_angle ) {
            _rank = rank;
            _label_angle = ( label_angle == null ) ? TreePanel.CLADE_LABEL_ANGLE.VERTICAL : label_angle;
        }

        String getRank() {
            return _rank;
        }

        TreePanel.CLADE_LABEL_ANGLE getLabelAngle() {
            return _label_angle;
        }

        @Override
        public String toString() {
            return _rank + "/" + _label_angle;
        }
    }

    private final Spec            _spec;
    private final List<CladeBand> _bands;

    CladeLevel( final Spec spec, final List<CladeBand> bands ) {
        _spec = spec;
        _bands = ( bands == null ) ? Collections.<CladeBand> emptyList() : bands;
    }

    Spec getSpec() {
        return _spec;
    }

    String getRank() {
        return _spec.getRank();
    }

    TreePanel.CLADE_LABEL_ANGLE getLabelAngle() {
        return _spec.getLabelAngle();
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
     * The specs in DRAWING order -- finest rank first, i.e. nearest the tips -- with blanks and duplicate ranks
     * dropped and the list capped at {@link #MAX_LEVELS}. A stable sort, so two ranks of equal (or unknown)
     * depth keep the order the user gave them. Pure; the single source of the placement order for the
     * rectangular columns, the circular rings and every space reservation derived from them.
     */
    static List<Spec> order( final List<Spec> specs ) {
        final List<Spec> out = new ArrayList<>();
        if ( specs == null ) {
            return out;
        }
        final Set<String> seen = new LinkedHashSet<>();
        for ( final Spec s : specs ) {
            if ( ( s == null ) || ForesterUtil.isEmpty( s.getRank() ) ) {
                continue;
            }
            if ( seen.add( s.getRank().toLowerCase() ) ) {
                out.add( s );
            }
        }
        // stable: List.sort is guaranteed stable, so equal-depth specs keep the user's relative order
        out.sort( Comparator.comparingInt( ( final Spec s ) -> rankDepth( s.getRank() ) ).reversed() );
        return ( out.size() > MAX_LEVELS ) ? new ArrayList<>( out.subList( 0, MAX_LEVELS ) ) : out;
    }
}
