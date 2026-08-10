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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * The single, immutable representation of a resolved taxon and its lineage -- the one shape returned by
 * {@link TaxonomicLineageService} and persisted by {@link TaxonomyDiskCache}, serving every consumer:
 * <ul>
 * <li>the rank colorizer / clade-bands, via {@link #at(String)} (constant-time "what is this taxon's
 *     <i>order</i>?", case-insensitive, aliasing {@code domain}&harr;{@code superkingdom});</li>
 * <li>Spine B (assignment keyed by tax-id), via {@link #taxIdAt(String)};</li>
 * <li>the "Fetch Sequence &amp; Taxonomic Data" tool, via the queried taxon's own
 *     {@link #getScientificName()}/{@link #ownRankOrNull()}/{@link #getTaxId()}/{@link #getCommonName()}
 *     plus the flat {@link #lineageNames()}.</li>
 * </ul>
 *
 * <p>Ancestors are stored root&rarr;parent (NOT including the taxon itself), each an {@link Ancestor}
 * {@code {name, rank, tax-id}} where the rank may be empty/"no rank"/"clade" (kept verbatim) and the
 * per-ancestor tax-id may be {@code null} (absent on pre-tax-id cached entries). A single rank&rarr;
 * {@code {name, tax-id}} lookup map and the flat name lineage are built ONCE in the constructor
 * (non-Linnaean "no rank"/"clade"/empty ranks skipped, keys lower-cased, the taxon's own rank folded in
 * last so it wins over a same-rank ancestor), so an instance is deeply immutable and safe to publish
 * across the EDT/worker boundary. Any field may be {@code null}/empty.
 */
public final class TaxonLineage {

    static final String                       NO_RANK = "no rank";
    static final String                       CLADE   = "clade";

    /** An unresolved / empty taxon (no scientific name, no ancestors). */
    public static final TaxonLineage          EMPTY   = new TaxonLineage( null, null, null, null, null );

    private final String                      _tax_id;
    private final String                      _rank;
    private final String                      _scientific_name;
    private final String                      _common_name;
    private final List<Ancestor>              _ancestors;
    // rank (lower-cased Linnaean) -> {name, tax-id}. ONE map (not two parallel ones) so at(rank) and
    // taxIdAt(rank) can never disagree: folding the own taxon in last replaces a same-rank entry WHOLE.
    private final Map<String, String[]>       _by_rank;
    private final List<String>                _lineage_names;

    public TaxonLineage( final String tax_id,
                         final String rank,
                         final String scientific_name,
                         final String common_name,
                         final List<Ancestor> ancestors ) {
        _tax_id = tax_id;
        _rank = rank;
        _scientific_name = scientific_name;
        _common_name = common_name;
        _ancestors = ( ancestors == null ) ? Collections.<Ancestor> emptyList()
                : Collections.unmodifiableList( new ArrayList<Ancestor>( ancestors ) );
        final Map<String, String[]> by_rank = new LinkedHashMap<String, String[]>();
        final List<String> lineage = new ArrayList<String>();
        for( final Ancestor a : _ancestors ) {
            addRankedEntry( by_rank, a.getRank(), a.getName(), a.getTaxId() );
            lineage.add( a.getName() );
        }
        // the queried taxon itself (its own rank/name/tax-id -- e.g. genus "Felis" for a Felis query); folded in
        // LAST so its {name, tax-id} replaces an ancestor entry of the same rank as ONE unit (no name/id mismatch)
        addRankedEntry( by_rank, _rank, _scientific_name, _tax_id );
        if ( _scientific_name != null ) {
            lineage.add( _scientific_name );
        }
        _by_rank = Collections.unmodifiableMap( by_rank );
        _lineage_names = Collections.unmodifiableList( lineage );
    }

    private static void addRankedEntry( final Map<String, String[]> by_rank, final String rank, final String name,
                                        final String tax_id ) {
        if ( isLinnaeanRank( rank ) && !isEmpty( name ) ) {
            by_rank.put( rank.toLowerCase( Locale.ROOT ), new String[] { name, isEmpty( tax_id ) ? null : tax_id } );
        }
    }

    /** A rank that a colorization can key on: present, and not NCBI's non-Linnaean "no rank"/"clade". */
    private static boolean isLinnaeanRank( final String rank ) {
        return !isEmpty( rank ) && !NO_RANK.equalsIgnoreCase( rank ) && !CLADE.equalsIgnoreCase( rank );
    }

    /** The queried taxon's own NCBI tax-id (may be {@code null}). */
    public String getTaxId() {
        return _tax_id;
    }

    /** The taxon's own <b>raw</b> NCBI rank (may be null / "no rank" / "clade"). Use for serialization; consumers
     *  that want a Linnaean-or-nothing rank should use {@link #ownRankOrNull()}. */
    public String getRank() {
        return _rank;
    }

    /** The taxon's own rank, or {@code null} for a non-Linnaean rank ("no rank"/"clade") -- for the Fetch writer. */
    public String ownRankOrNull() {
        return isLinnaeanRank( _rank ) ? _rank : null;
    }

    public String getScientificName() {
        return _scientific_name;
    }

    public String getCommonName() {
        return _common_name;
    }

    /** Ancestors root&rarr;parent (excluding the taxon itself). */
    public List<Ancestor> getAncestors() {
        return _ancestors;
    }

    /**
     * The scientific name at {@code rank} (case-insensitive), or {@code null} if this lineage has none.
     * "domain" and "superkingdom" are aliases: NCBI labels the domain level "superkingdom", so a tree
     * annotated with either spelling resolves.
     */
    public String at( final String rank ) {
        final String[] e = lookup( rank );
        return ( e == null ) ? null : e[ 0 ];
    }

    /** The NCBI tax-id at {@code rank} (case-insensitive, same aliasing as {@link #at}), or {@code null} when the
     *  rank is absent OR present without a tax-id. Always consistent with {@link #at} -- they read one entry. */
    public String taxIdAt( final String rank ) {
        final String[] e = lookup( rank );
        return ( e == null ) ? null : e[ 1 ];
    }

    private String[] lookup( final String rank ) {
        if ( rank == null ) {
            return null;
        }
        final String r = rank.toLowerCase( Locale.ROOT );
        final String[] e = _by_rank.get( r );
        if ( e != null ) {
            return e;
        }
        if ( r.equals( "domain" ) ) {
            return _by_rank.get( "superkingdom" );
        }
        if ( r.equals( "superkingdom" ) ) {
            return _by_rank.get( "domain" );
        }
        return null;
    }

    /** Unmodifiable rank&rarr;name view (non-Linnaean ranks skipped, keys lower-cased), root&rarr;tip order. */
    public Map<String, String> namesByRank() {
        final Map<String, String> m = new LinkedHashMap<String, String>();
        for( final Map.Entry<String, String[]> e : _by_rank.entrySet() ) {
            m.put( e.getKey(), e.getValue()[ 0 ] );
        }
        return Collections.unmodifiableMap( m );
    }

    /** The flat lineage of scientific names root&rarr;taxon: all ancestor names (incl. no-rank ones), then the
     *  taxon itself (when it has a scientific name). Precomputed. Used by the Fetch tool. */
    public List<String> lineageNames() {
        return _lineage_names;
    }

    /** True when nothing usable was resolved (no scientific name and no ancestors). */
    public boolean isEmpty() {
        return isEmpty( _scientific_name ) && _ancestors.isEmpty();
    }

    private static boolean isEmpty( final String s ) {
        return ( s == null ) || s.isEmpty();
    }

    @Override
    public String toString() {
        return "TaxonLineage[" + _scientific_name + " (" + _rank + ", id=" + _tax_id + "), " + namesByRank() + "]";
    }

    /**
     * One ancestor in the lineage. Immutable. Any field may be {@code null}/empty; the rank is kept verbatim
     * (may be "no rank"/"clade") and the tax-id may be {@code null} (absent on pre-tax-id cached entries).
     */
    public static final class Ancestor {

        private final String _name;
        private final String _rank;
        private final String _tax_id;

        public Ancestor( final String name, final String rank, final String tax_id ) {
            _name = name;
            _rank = rank;
            _tax_id = tax_id;
        }

        public String getName() {
            return _name;
        }

        public String getRank() {
            return _rank;
        }

        public String getTaxId() {
            return _tax_id;
        }
    }
}
