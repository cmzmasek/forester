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

package org.forester.ws.seqdb;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The full, raw taxonomy detail for one taxon as parsed from an NCBI efetch response, kept in a single
 * shape so that <b>both</b> taxonomy consumers can be served from one cache entry: the rank colorizer
 * (which needs a rank&rarr;name map &mdash; {@link RankedLineage}) and the "Fetch Sequence &amp;
 * Taxonomic Data" tool (which needs the flat name lineage &mdash; {@link ResolvedTaxonomy}). It is the
 * unit persisted by {@link TaxonomyDiskCache}; the two views are derived in
 * {@link NcbiTaxonomyLineageService}. Immutable.
 *
 * <p>Ancestors are stored root&rarr;parent (NOT including the taxon itself), each as a
 * {@code {scientific-name, rank}} pair where the rank may be empty/"no rank"/"clade" (kept verbatim so
 * the derivations can filter as each consumer requires). Any field may be {@code null}/empty.
 */
final class CachedTaxon {

    private final String         _tax_id;
    private final String         _rank;
    private final String         _scientific_name;
    private final String         _common_name;
    private final List<String[]> _ancestors;

    CachedTaxon( final String tax_id,
                 final String rank,
                 final String scientific_name,
                 final String common_name,
                 final List<String[]> ancestors ) {
        _tax_id = tax_id;
        _rank = rank;
        _scientific_name = scientific_name;
        _common_name = common_name;
        _ancestors = ( ancestors == null ) ? Collections.<String[]> emptyList()
                : Collections.unmodifiableList( new ArrayList<String[]>( ancestors ) );
    }

    String getTaxId() {
        return _tax_id;
    }

    /** The taxon's own raw NCBI rank (may be null / "no rank" / "clade"). */
    String getRank() {
        return _rank;
    }

    String getScientificName() {
        return _scientific_name;
    }

    String getCommonName() {
        return _common_name;
    }

    /** Ancestors root&rarr;parent (excluding the taxon itself); each entry is {@code {name, rank}}. */
    List<String[]> getAncestors() {
        return _ancestors;
    }

    /** True when nothing usable was parsed (no scientific name and no ancestors). */
    boolean isEmpty() {
        return ( ( _scientific_name == null ) || _scientific_name.isEmpty() ) && _ancestors.isEmpty();
    }
}
