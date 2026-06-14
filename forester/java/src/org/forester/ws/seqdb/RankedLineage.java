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

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

/**
 * An immutable taxonomic lineage with a scientific name for each (standard) rank -- e.g. for
 * <i>Felis catus</i>: {@code superkingdom=Eukaryota, ..., class=Mammalia, order=Carnivora,
 * family=Felidae, genus=Felis, species=Felis catus}. The point is constant-time
 * {@link #at(String)} so a caller can ask "what is this taxon's <i>order</i>?" without walking a
 * tree or making a second query.
 *
 * <p>Ranks are compared case-insensitively. Built from a single taxonomy-database response (see
 * {@code NcbiTaxonomyLineageService}); a clean replacement for resolving ranks one ancestor at a
 * time against the old UniProt code.
 */
public final class RankedLineage {

    /** A lineage with no ranked entries (e.g. an unresolved or rank-less taxon). */
    public static final RankedLineage EMPTY = new RankedLineage( new LinkedHashMap<String, String>() );
    private final Map<String, String> _by_rank;

    /**
     * @param rank_to_name ordered rank&rarr;scientific-name pairs (root&rarr;tip); keys are
     *            lower-cased on the way in, so lookups are case-insensitive. A defensive copy is made.
     */
    public RankedLineage( final Map<String, String> rank_to_name ) {
        final Map<String, String> m = new LinkedHashMap<String, String>();
        if ( rank_to_name != null ) {
            for( final Map.Entry<String, String> e : rank_to_name.entrySet() ) {
                if ( ( e.getKey() != null ) && ( e.getValue() != null ) ) {
                    m.put( e.getKey().toLowerCase( Locale.ROOT ), e.getValue() );
                }
            }
        }
        _by_rank = Collections.unmodifiableMap( m );
    }

    /**
     * The scientific name at {@code rank} (case-insensitive), or {@code null} if this lineage has
     * none. "domain" and "superkingdom" are treated as aliases: NCBI labels the domain level
     * "superkingdom", so a tree annotated with either spelling resolves against this lineage.
     */
    public String at( final String rank ) {
        if ( rank == null ) {
            return null;
        }
        final String r = rank.toLowerCase( Locale.ROOT );
        final String v = _by_rank.get( r );
        if ( v != null ) {
            return v;
        }
        if ( r.equals( "domain" ) ) {
            return _by_rank.get( "superkingdom" );
        }
        if ( r.equals( "superkingdom" ) ) {
            return _by_rank.get( "domain" );
        }
        return null;
    }

    /** True if no ranked entries are present. */
    public boolean isEmpty() {
        return _by_rank.isEmpty();
    }

    /** Unmodifiable rank&rarr;name view (ranks lower-cased), in insertion (root&rarr;tip) order. */
    public Map<String, String> asMap() {
        return _by_rank;
    }

    @Override
    public String toString() {
        return _by_rank.toString();
    }
}
