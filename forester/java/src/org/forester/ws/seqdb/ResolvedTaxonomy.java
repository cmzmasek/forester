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
 * The full taxonomy detail for one taxon as returned by a taxonomy database: scientific name, rank,
 * tax-id, common name, and the lineage of scientific names (root&rarr;taxon). Immutable. Produced by
 * {@code NcbiTaxonomyLineageService.resolveTaxonomy} and consumed by the "Fetch Sequence &amp;
 * Taxonomic Data" tool to enrich a node's {@code Taxonomy}.
 *
 * <p>Any field may be {@code null}/empty (the database may not supply it); {@link #isEmpty()} is true
 * when nothing usable was resolved.
 */
public final class ResolvedTaxonomy {

    /** An empty result (nothing resolved). */
    public static final ResolvedTaxonomy EMPTY = new ResolvedTaxonomy( null, null, null, null, null );
    private final String                 _scientific_name;
    private final String                 _rank;
    private final String                 _tax_id;
    private final String                 _common_name;
    private final List<String>           _lineage;

    /** @param lineage root&rarr;taxon scientific names; copied defensively, never null after construction. */
    public ResolvedTaxonomy( final String scientific_name,
                             final String rank,
                             final String tax_id,
                             final String common_name,
                             final List<String> lineage ) {
        _scientific_name = scientific_name;
        _rank = rank;
        _tax_id = tax_id;
        _common_name = common_name;
        _lineage = ( lineage == null ) ? Collections.<String> emptyList()
                : Collections.unmodifiableList( new ArrayList<String>( lineage ) );
    }

    public String getScientificName() {
        return _scientific_name;
    }

    /** The taxon's own rank (lower-cased Linnaean rank), or {@code null} if none / "no rank"/"clade". */
    public String getRank() {
        return _rank;
    }

    /** The NCBI tax-id (a numeric string), or {@code null}. */
    public String getTaxId() {
        return _tax_id;
    }

    public String getCommonName() {
        return _common_name;
    }

    /** Unmodifiable lineage of scientific names, root&rarr;taxon (never null; may be empty). */
    public List<String> getLineage() {
        return _lineage;
    }

    /** True when nothing usable was resolved (no scientific name and no lineage). */
    public boolean isEmpty() {
        return ( ( _scientific_name == null ) || _scientific_name.isEmpty() ) && _lineage.isEmpty();
    }

    @Override
    public String toString() {
        return "ResolvedTaxonomy[" + _scientific_name + ", rank=" + _rank + ", taxId=" + _tax_id + ", common="
                + _common_name + ", lineage=" + _lineage + "]";
    }
}
