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

import org.forester.archaeopteryx.TreePanelUtil.RankTaxon;

/**
 * One node's taxon at a rank TOGETHER WITH WHERE IT CAME FROM -- the record "Annotate Clades by Rank" keeps for
 * every tip so its report can say, per tip, which source placed it, why it stayed unplaced, and when the database
 * contradicted the tree. Immutable. The taxon is {@code null} exactly when the {@link Source} is an UNPLACED one.
 */
final class TaxonPlacement {

    /** Where a placement came from, or why there is none. */
    enum Source {
        /** The node's own {@code <taxonomy>} is annotated at exactly this rank. */
        OWN_ANNOTATION( "own annotation", true ),
        /** The node's own per-rank metadata field, in any namespace ({@code vipr:Genus}, {@code ncbi:family},
         *  {@code gtdb:family}). */
        OWN_METADATA( "own metadata field", true ),
        /** The nearest ancestor annotated at this rank (the clade the tree says it sits in). */
        ANCESTOR_ANNOTATION( "ancestor's annotation", true ),
        /** The database lineage of the node's own NCBI tax-id. */
        DATABASE_BY_ID( "database, by its NCBI id", true ),
        /** The database lineage of the node's own taxonomy name (scientific, code or common). */
        DATABASE_BY_NAME( "database, by its taxonomy name", true ),
        /** The database lineage found from the node's label (its name, an accession, an ancestor's name). */
        DATABASE_BY_LABEL( "database, by its label", true ),
        /** Nothing in the tree places it, and the source does not look names up. */
        UNPLACED_NO_ANNOTATION( "unplaced: nothing in the tree at this rank", false ),
        /** It would be looked up, but has not been yet. */
        UNPLACED_NOT_LOOKED_UP( "unplaced: not looked up yet", false ),
        /** The database does not know the name. */
        UNPLACED_UNKNOWN( "unplaced: not in the database", false ),
        /** Several taxa carry the name; none was picked. */
        UNPLACED_AMBIGUOUS( "unplaced: ambiguous name", false ),
        /** The database knows it, but its lineage has no taxon at this rank (NCBI: "clade"/"no rank"). */
        UNPLACED_NO_SUCH_RANK( "unplaced: its lineage has no taxon at this rank", false );

        private final String  _label;
        private final boolean _placed;

        Source( final String label, final boolean placed ) {
            _label = label;
            _placed = placed;
        }

        String getLabel() {
            return _label;
        }

        boolean isPlaced() {
            return _placed;
        }

        /** True for the three database sources. */
        boolean isDatabase() {
            return ( this == DATABASE_BY_ID ) || ( this == DATABASE_BY_NAME ) || ( this == DATABASE_BY_LABEL );
        }
    }

    private final RankTaxon _taxon;
    private final Source    _source;
    private final String    _detail;
    private final RankTaxon _tree_says;

    /**
     * @param taxon     the taxon, or null when unplaced
     * @param source    where it came from / why there is none
     * @param detail    what was consulted: the ancestor's label, the query the database was asked, the
     *                  candidates of an ambiguous name; may be null
     * @param tree_says the ancestor annotation the DATABASE source overrode (a conflict), else null
     */
    TaxonPlacement( final RankTaxon taxon, final Source source, final String detail, final RankTaxon tree_says ) {
        if ( ( taxon == null ) == source.isPlaced() ) {
            throw new IllegalArgumentException( "a placed source needs a taxon and an unplaced one must not have one: "
                    + source + " / " + taxon );
        }
        _taxon = taxon;
        _source = source;
        _detail = detail;
        _tree_says = tree_says;
    }

    static TaxonPlacement unplaced( final Source why, final String detail ) {
        return new TaxonPlacement( null, why, detail, null );
    }

    RankTaxon getTaxon() {
        return _taxon;
    }

    Source getSource() {
        return _source;
    }

    String getDetail() {
        return _detail;
    }

    /** The tree's ancestor annotation at this rank that the database overrode, or null when they agree or there
     *  was none. */
    RankTaxon getTreeSays() {
        return _tree_says;
    }

    boolean isPlaced() {
        return _source.isPlaced();
    }

    /** True when the database placed this node somewhere other than where the tree's ancestor annotation says. */
    boolean isConflict() {
        return _tree_says != null;
    }

    @Override
    public String toString() {
        return ( ( _taxon == null ) ? "-" : _taxon.toString() ) + " <" + _source + ( ( _detail == null ) ? "" : ": " + _detail )
                + ( ( _tree_says == null ) ? "" : "; tree says " + _tree_says ) + ">";
    }
}
