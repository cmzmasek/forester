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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.util.ForesterUtil;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A clean {@link TaxonomicLineageService} backed by NCBI Taxonomy via the E-utilities REST API:
 * <i>esearch</i> maps a name to a tax-id, one <i>efetch</i> returns the full
 * lineage in a single XML response (so we never resolve ranks one ancestor at a time).
 *
 * <p><b>A name is never guessed.</b> The search runs in three narrowing stages -- the quoted name as a
 * {@code [Scientific Name]}, then as any name ({@code [All Names]}: synonyms, common names, "human"), then
 * unrestricted, as NCBI's own search box does -- and the FIRST stage with a hit decides. One hit resolves the
 * name; several hits make it AMBIGUOUS ({@link #ambiguousCandidates}): it is left unplaced, because the first hit
 * of a homonym is a coin toss ("Bacteria" used to resolve to the stick-insect genus <i>Bacteria</i> Latreille,
 * "Bacillus" still would). Quoting keeps a name's own brackets and parentheses ("[Haemophilus] ducreyi",
 * "Buchnera aphidicola (Baizongia pistaciae)") from being read as search syntax. A tax-id carried IN the label
 * ({@code OX=47803}, the UniProt FASTA header) is used directly ({@link #embeddedTaxId}). It serves every
 * consumer -- the rank colorizer / clade-bands ({@link #lineageOf}/{@link #fetch} &rarr;
 * {@link TaxonLineage#at}) and the "Fetch Sequence &amp; Taxonomic Data" tool (the taxon's own fields +
 * {@link TaxonLineage#lineageNames}) -- from one {@link TaxonLineage} per taxon. Results -- including
 * "not found" -- are cached for the life of the process.
 *
 * <p>Intentionally new, self-contained code: HTTP/throttle/timeout/XML plumbing comes from
 * {@link WsHttp}, not the legacy taxonomy classes. The XML parser is split from the network call so it
 * can be unit-tested from captured fixtures.
 */
public final class NcbiTaxonomyLineageService implements TaxonomicLineageService {

    private static final String                  ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&tool=Archaeopteryx&term=";
    private static final String                  EFETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&tool=Archaeopteryx&retmode=xml&id=";
    /** The esearch field restrictions, tried in this order; the first stage with any hit decides. */
    static final String[]                        SEARCH_STAGES = { "[Scientific Name]", "[All Names]", "" };
    // "OX=47803" as a UniProt FASTA header writes it: at the start or after whitespace / a '|' separator
    private static final Pattern                 EMBEDDED_TAX_ID = Pattern.compile( "(?:^|[\\s|])OX=(\\d+)(?![\\w])" );
    // EMPTY entries double as a negative cache (taxon was queried, nothing found -- or found ambiguous).
    private final Map<String, TaxonLineage>      _cache        = Collections
            .synchronizedMap( new HashMap<String, TaxonLineage>() );
    // name -> the several taxa that carry it (an EMPTY in _cache with an entry here = ambiguous, not unknown)
    private final Map<String, List<TaxonLineage>> _ambiguous   = Collections
            .synchronizedMap( new HashMap<String, List<TaxonLineage>>() );
    // persistent (cross-session) cache of positive resolutions; best-effort, never a dependency.
    private final TaxonomyDiskCache              _disk;
    private final Http                           _http;
    private volatile boolean                     _loaded;

    private static NcbiTaxonomyLineageService    _shared;

    /** The transport seam: a URL to its response body. Production is {@link WsHttp#httpGet}; tests answer from
     *  canned fixtures, so the staged search and the ambiguity rule run without a network. */
    interface Http {
        String get( String url ) throws IOException;
    }

    public NcbiTaxonomyLineageService() {
        this( WsHttp::httpGet, new TaxonomyDiskCache() );
    }

    /** For tests: an in-memory transport and an explicit (temporary) disk cache. */
    NcbiTaxonomyLineageService( final Http http, final TaxonomyDiskCache disk ) {
        _http = http;
        _disk = disk;
    }

    /**
     * The process-wide shared instance, so the rank colorizer, the "Fetch Sequence &amp; Taxonomic
     * Data" tool and the Settings dialog all share one in-memory cache <i>and</i> one disk cache.
     */
    public static synchronized NcbiTaxonomyLineageService getShared() {
        if ( _shared == null ) {
            _shared = new NcbiTaxonomyLineageService();
        }
        return _shared;
    }

    @Override
    public TaxonLineage lineageOf( final String taxon ) {
        if ( ForesterUtil.isEmpty( taxon ) ) {
            return null;
        }
        ensureLoaded();
        return _cache.get( key( taxon ) );
    }

    /**
     * Warms the in-memory caches from disk on a background daemon thread, so the first (EDT) rank
     * colorize doesn't pay the disk read + compaction. Safe to call repeatedly; a no-op once loaded.
     * Intended for GUI startup -- tests use the constructor directly and load lazily.
     */
    public void primeAsync() {
        if ( _loaded ) {
            return;
        }
        final Thread t = new Thread( () -> {
            try {
                ensureLoaded();
            }
            catch ( final Throwable th ) {
                // best-effort warm-up; the lazy load on the first real lookup still covers correctness
            }
        }, "aptx-taxonomy-cache-prime" );
        t.setDaemon( true );
        t.start();
    }

    /** Seeds the in-memory cache from disk exactly once (lazily, on first use). */
    private void ensureLoaded() {
        if ( _loaded ) {
            return;
        }
        synchronized ( this ) {
            if ( _loaded ) {
                return;
            }
            for( final Map.Entry<String, TaxonLineage> e : _disk.load().entrySet() ) {
                if ( !e.getValue().isEmpty() ) {
                    _cache.put( e.getKey(), e.getValue() );
                }
            }
            _loaded = true;
        }
    }

    // ---- persistent-cache control (used by the Settings dialog) --------------------------------

    /** A snapshot of the on-disk cache (location, size, age, availability) for display. */
    public TaxonomyCacheStatus getCacheStatus() {
        return _disk.status();
    }

    /** Deletes the on-disk cache and clears the in-memory cache. */
    public synchronized void clearPersistentCache() {
        _disk.clear();
        _cache.clear();
    }

    public boolean isPersistentCacheEnabled() {
        return _disk.isEnabled();
    }

    /** Turns the persistent cache on/off; enabling re-seeds the in-memory caches from disk. */
    public synchronized void setPersistentCacheEnabled( final boolean enabled ) {
        final boolean was = _disk.isEnabled();
        _disk.setEnabled( enabled );
        if ( enabled && !was ) {
            _loaded = false; // re-seed on the next lookup now that the disk cache is in play again
        }
    }

    @Override
    public TaxonLineage fetch( final String taxon ) throws IOException {
        if ( ForesterUtil.isEmpty( taxon ) ) {
            return TaxonLineage.EMPTY;
        }
        ensureLoaded();
        final String k = key( taxon );
        synchronized ( _cache ) {
            if ( _cache.containsKey( k ) ) {
                return _cache.get( k );
            }
        }
        // A bare NCBI tax-id (e.g. the organism id read from a sequence entry), or one carried in the label
        // ("OX=47803"), is authoritative: efetch it directly, skipping esearch and its name-ambiguity. Otherwise map
        // the name to a tax-id via the staged esearch.
        final String id;
        if ( isTaxId( taxon ) ) {
            id = taxon.trim();
        }
        else {
            final String embedded = embeddedTaxId( taxon );
            if ( embedded != null ) {
                id = embedded;
            }
            else {
                final List<String> ids = searchName( taxon );
                if ( ids.size() > 1 ) {
                    // several taxa carry this name: never pick one. Remember the candidates (for the report) and
                    // cache the name as unplaced, in memory only, so it is neither re-queried nor re-prompted.
                    _ambiguous.put( k, parseEfetchAll( _http.get( EFETCH + WsHttp.encode( String.join( ",", ids ) ) ) ) );
                    _cache.put( k, TaxonLineage.EMPTY );
                    return TaxonLineage.EMPTY;
                }
                id = ids.isEmpty() ? null : ids.get( 0 );
            }
        }
        if ( ForesterUtil.isEmpty( id ) ) {
            // esearch found nothing -- a definitive negative; cache it (in memory only) so we never re-query.
            _cache.put( k, TaxonLineage.EMPTY );
            return TaxonLineage.EMPTY;
        }
        final TaxonLineage parsed = parseEfetchFull( _http.get( EFETCH + WsHttp.encode( id ) ) );
        // "usable" = TaxonLineage.isEmpty() is false: a scientific name OR any ancestors. This is deliberately more
        // lenient than the pre-Spine-A colorizer gate ("has a Linnaean rank"): a genuinely rank-less taxon (a
        // top-level taxon, or one NCBI has no Linnaean ranks for) is now cached positive instead of re-prompted /
        // re-fetched forever -- a net win. Trade-off: a well-formed-but-partial efetch (name, no ranks) is also
        // cached (till the 30-day TTL); truncation usually yields malformed XML -> empty -> uncached, so this is rare.
        if ( ( parsed != null ) && !parsed.isEmpty() ) {
            _cache.put( k, parsed );
            _disk.put( k, parsed );
            return parsed;
        }
        // esearch found the taxon but efetch returned nothing usable (empty) -- possibly a truncated/transient
        // response. Do NOT cache, so a later attempt can retry.
        return TaxonLineage.EMPTY;
    }

    @Override
    public List<TaxonLineage> ambiguousCandidates( final String taxon ) {
        if ( ForesterUtil.isEmpty( taxon ) ) {
            return Collections.emptyList();
        }
        final List<TaxonLineage> c = _ambiguous.get( key( taxon ) );
        return ( c == null ) ? Collections.<TaxonLineage> emptyList() : Collections.unmodifiableList( c );
    }

    /**
     * The tax-ids NCBI returns for {@code name}, from the first {@link #SEARCH_STAGES stage} that returns any:
     * empty = unknown, one = resolved, several = ambiguous. The name is sent in double quotes so its own brackets,
     * parentheses and spaces are data, not Entrez syntax; the unrestricted last stage is NCBI's own search-box
     * behaviour and catches forms NCBI lists differently ("Bacillus subtilis (strain 168)", "Oryza sativa subsp.
     * japonica"). Network.
     */
    private List<String> searchName( final String name ) throws IOException {
        final String quoted = '"' + name.trim().replace( '"', ' ' ) + '"';
        for( final String stage : SEARCH_STAGES ) {
            final String term = stage.isEmpty() ? name.trim() : ( quoted + stage );
            final List<String> ids = parseEsearchIds( _http.get( ESEARCH + WsHttp.encode( term ) ) );
            if ( !ids.isEmpty() ) {
                return ids;
            }
        }
        return Collections.emptyList();
    }

    /**
     * The NCBI tax-id a label carries as a UniProt-style {@code OX=} qualifier ("... OS=Fusarium phyllophilum
     * OX=47803 GN=..."), or {@code null}. Pure.
     */
    public static String embeddedTaxId( final String label ) {
        if ( ForesterUtil.isEmpty( label ) ) {
            return null;
        }
        final Matcher m = EMBEDDED_TAX_ID.matcher( label );
        return m.find() ? m.group( 1 ) : null;
    }

    /** True if {@code s} is a non-empty run of digits (an NCBI tax-id). */
    static boolean isTaxId( final String s ) {
        final String t = s.trim();
        if ( t.isEmpty() ) {
            return false;
        }
        for( int i = 0; i < t.length(); ++i ) {
            if ( !Character.isDigit( t.charAt( i ) ) ) {
                return false;
            }
        }
        return true;
    }

    private static String key( final String taxon ) {
        return taxon.trim().toLowerCase( Locale.ROOT );
    }

    /** The first {@code <Id>} under {@code <IdList>} in an esearch response, or {@code null}. Pure. */
    static String parseEsearchFirstId( final String xml ) {
        final List<String> ids = parseEsearchIds( xml );
        return ids.isEmpty() ? null : ids.get( 0 );
    }

    /** Every {@code <Id>} under {@code <IdList>} in an esearch response, in NCBI's order; empty for none or for a
     *  malformed response. Pure. */
    static List<String> parseEsearchIds( final String xml ) {
        final List<String> out = new ArrayList<String>();
        final Document doc = WsHttp.parseXml( xml );
        if ( doc == null ) {
            return out;
        }
        final Element list = WsHttp.firstChildElement( doc.getDocumentElement(), "IdList" );
        if ( list == null ) {
            return out;
        }
        for( final Element e : WsHttp.childElements( list, "Id" ) ) {
            final String id = WsHttp.text( e );
            if ( !ForesterUtil.isEmpty( id ) ) {
                out.add( id.trim() );
            }
        }
        return out;
    }

    /**
     * Parses an NCBI efetch taxonomy XML response into a {@link TaxonLineage}: the queried taxon's
     * scientific name, rank, tax-id and common name, plus its ancestors ({@code <LineageEx>}
     * root&rarr;parent) each as {@code {name, rank, tax-id}} (ranks kept verbatim, including "no
     * rank"/"clade"; the ranked lookup maps are derived inside {@link TaxonLineage}). Pure -- no I/O;
     * never throws.
     */
    static TaxonLineage parseEfetchFull( final String xml ) {
        final List<TaxonLineage> all = parseEfetchAll( xml );
        return all.isEmpty() ? TaxonLineage.EMPTY : all.get( 0 );
    }

    /** Every {@code <Taxon>} of an efetch response (one per id of a multi-id request), each parsed as
     *  {@link #parseEfetchFull} parses a single one; empty for none or a malformed response. Pure. */
    static List<TaxonLineage> parseEfetchAll( final String xml ) {
        final List<TaxonLineage> out = new ArrayList<TaxonLineage>();
        final Document doc = WsHttp.parseXml( xml );
        if ( doc == null ) {
            return out;
        }
        for( final Element taxon : WsHttp.childElements( doc.getDocumentElement(), "Taxon" ) ) {
            final TaxonLineage tl = parseTaxon( taxon );
            if ( !tl.isEmpty() ) {
                out.add( tl );
            }
        }
        return out;
    }

    private static TaxonLineage parseTaxon( final Element taxon ) {
        final String sci = trimOrNull( WsHttp.text( WsHttp.firstChildElement( taxon, "ScientificName" ) ) );
        final String tax_id = trimOrNull( WsHttp.text( WsHttp.firstChildElement( taxon, "TaxId" ) ) );
        final String rank = trimOrNull( WsHttp.text( WsHttp.firstChildElement( taxon, "Rank" ) ) );
        // common name: prefer the GenBank common name, else any common name
        String common = null;
        final Element other = WsHttp.firstChildElement( taxon, "OtherNames" );
        if ( other != null ) {
            common = trimOrNull( WsHttp.text( WsHttp.firstChildElement( other, "GenbankCommonName" ) ) );
            if ( common == null ) {
                common = trimOrNull( WsHttp.text( WsHttp.firstChildElement( other, "CommonName" ) ) );
            }
        }
        final List<TaxonLineage.Ancestor> ancestors = new ArrayList<TaxonLineage.Ancestor>();
        final Element lineage_ex = WsHttp.firstChildElement( taxon, "LineageEx" );
        if ( lineage_ex != null ) {
            for( final Element anc : WsHttp.childElements( lineage_ex, "Taxon" ) ) {
                final String n = trimOrNull( WsHttp.text( WsHttp.firstChildElement( anc, "ScientificName" ) ) );
                if ( n != null ) {
                    final String r = trimOrNull( WsHttp.text( WsHttp.firstChildElement( anc, "Rank" ) ) );
                    final String id = trimOrNull( WsHttp.text( WsHttp.firstChildElement( anc, "TaxId" ) ) );
                    ancestors.add( new TaxonLineage.Ancestor( n, ( r == null ) ? "" : r, id ) );
                }
            }
        }
        return new TaxonLineage( tax_id, rank, sci, common, ancestors );
    }

    private static String trimOrNull( final String s ) {
        if ( s == null ) {
            return null;
        }
        final String t = s.trim();
        return t.isEmpty() ? null : t;
    }
}
