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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.forester.util.ForesterUtil;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * A clean {@link TaxonomicLineageService} backed by NCBI Taxonomy via the E-utilities REST API:
 * one <i>esearch</i> maps a scientific name to a tax-id, one <i>efetch</i> returns the full
 * lineage in a single XML response (so we never resolve ranks one ancestor at a time). It serves both
 * the rank colorizer ({@link #fetch}/{@link #lineageOf} &rarr; {@link RankedLineage}) and the
 * "Fetch Sequence &amp; Taxonomic Data" tool ({@link #resolveTaxonomy} &rarr; full
 * {@link ResolvedTaxonomy}). Results -- including "not found" -- are cached for the life of the process.
 *
 * <p>Intentionally new, self-contained code: HTTP/throttle/timeout/XML plumbing comes from
 * {@link WsHttp}, not the legacy taxonomy classes. The XML parsers are split from the network call so
 * they can be unit-tested from captured fixtures.
 */
public final class NcbiTaxonomyLineageService implements TaxonomicLineageService, TaxonomyResolver {

    private static final String                  ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&tool=Archaeopteryx&term=";
    private static final String                  EFETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&tool=Archaeopteryx&retmode=xml&id=";
    // NCBI ranks that are not Linnaean levels a colorization can key on; skipped when building a ranked lineage.
    private static final String                  NO_RANK = "no rank";
    private static final String                  CLADE   = "clade";
    private static final CachedTaxon            EMPTY_TAXON   = new CachedTaxon( null, null, null, null, null );
    // null-valued entries double as a negative cache (taxon was queried, nothing found).
    private final Map<String, RankedLineage>     _cache        = Collections
            .synchronizedMap( new HashMap<String, RankedLineage>() );
    private final Map<String, ResolvedTaxonomy>  _detail_cache = Collections
            .synchronizedMap( new HashMap<String, ResolvedTaxonomy>() );
    // persistent (cross-session) cache of positive resolutions; best-effort, never a dependency.
    private final TaxonomyDiskCache              _disk         = new TaxonomyDiskCache();
    private volatile boolean                     _loaded;

    private static NcbiTaxonomyLineageService    _shared;

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
    public RankedLineage lineageOf( final String taxon ) {
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

    /** Seeds the in-memory caches from disk exactly once (lazily, on first use). */
    private void ensureLoaded() {
        if ( _loaded ) {
            return;
        }
        synchronized ( this ) {
            if ( _loaded ) {
                return;
            }
            for( final Map.Entry<String, CachedTaxon> e : _disk.load().entrySet() ) {
                final RankedLineage rl = toRankedLineage( e.getValue() );
                if ( !rl.isEmpty() ) {
                    _cache.put( e.getKey(), rl );
                }
                final ResolvedTaxonomy rt = toResolvedTaxonomy( e.getValue() );
                if ( !rt.isEmpty() ) {
                    _detail_cache.put( e.getKey(), rt );
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

    /** Deletes the on-disk cache and clears the in-memory caches. */
    public synchronized void clearPersistentCache() {
        _disk.clear();
        _cache.clear();
        _detail_cache.clear();
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
    public RankedLineage fetch( final String taxon ) throws IOException {
        if ( ForesterUtil.isEmpty( taxon ) ) {
            return RankedLineage.EMPTY;
        }
        ensureLoaded();
        final String k = key( taxon );
        synchronized ( _cache ) {
            if ( _cache.containsKey( k ) ) {
                return _cache.get( k );
            }
        }
        // A bare NCBI tax-id (e.g. the organism id read from a sequence entry) is authoritative: efetch it
        // directly, skipping esearch and its name-ambiguity. Otherwise map the name to a tax-id via esearch.
        final String id;
        if ( isTaxId( taxon ) ) {
            id = taxon.trim();
        }
        else {
            id = parseEsearchFirstId( WsHttp.httpGet( ESEARCH + WsHttp.encode( taxon ) ) );
        }
        if ( ForesterUtil.isEmpty( id ) ) {
            // esearch found nothing -- a definitive negative; cache it (in memory only) so we never re-query.
            _cache.put( k, RankedLineage.EMPTY );
            return RankedLineage.EMPTY;
        }
        final CachedTaxon ct = parseEfetchFull( WsHttp.httpGet( EFETCH + WsHttp.encode( id ) ) );
        final RankedLineage parsed = toRankedLineage( ct );
        if ( ( parsed != null ) && !parsed.isEmpty() ) {
            _cache.put( k, parsed );
            _disk.put( k, ct );
            return parsed;
        }
        // esearch found the taxon but efetch returned no usable ranked lineage -- possibly a
        // truncated/transient response. Do NOT cache this as a negative, so a later attempt can retry.
        return RankedLineage.EMPTY;
    }

    /**
     * Resolves the full taxonomy detail (scientific name, rank, NCBI tax-id, common name, lineage) for
     * {@code query} (a scientific name, code, or tax-id), caching the result. Returns
     * {@link ResolvedTaxonomy#EMPTY} (cached) when esearch finds nothing; an esearch-hit whose efetch
     * yields nothing usable is NOT cached (likely transient) so a later attempt can retry. Does network
     * I/O -- call off the EDT.
     *
     * @throws IOException on a connection/transport failure.
     */
    public ResolvedTaxonomy resolveTaxonomy( final String query ) throws IOException {
        if ( ForesterUtil.isEmpty( query ) ) {
            return ResolvedTaxonomy.EMPTY;
        }
        ensureLoaded();
        final String k = key( query );
        synchronized ( _detail_cache ) {
            if ( _detail_cache.containsKey( k ) ) {
                return _detail_cache.get( k );
            }
        }
        // A bare NCBI tax-id is the authoritative key: efetch it directly (skip esearch), which also avoids
        // the name-ambiguity clobber when the caller already has the organism's tax-id from a sequence entry.
        final String id;
        if ( isTaxId( query ) ) {
            id = query.trim();
        }
        else {
            id = parseEsearchFirstId( WsHttp.httpGet( ESEARCH + WsHttp.encode( query ) ) );
        }
        if ( ForesterUtil.isEmpty( id ) ) {
            _detail_cache.put( k, ResolvedTaxonomy.EMPTY );
            return ResolvedTaxonomy.EMPTY;
        }
        final CachedTaxon ct = parseEfetchFull( WsHttp.httpGet( EFETCH + WsHttp.encode( id ) ) );
        final ResolvedTaxonomy detail = toResolvedTaxonomy( ct );
        if ( ( detail != null ) && !detail.isEmpty() ) {
            _detail_cache.put( k, detail );
            _disk.put( k, ct );
            return detail;
        }
        return ResolvedTaxonomy.EMPTY;
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
        final Document doc = WsHttp.parseXml( xml );
        if ( doc == null ) {
            return null;
        }
        final NodeList ids = doc.getElementsByTagName( "Id" );
        if ( ids.getLength() < 1 ) {
            return null;
        }
        final String id = ids.item( 0 ).getTextContent();
        return ForesterUtil.isEmpty( id ) ? null : id.trim();
    }

    /**
     * Builds a {@link RankedLineage} (rank&rarr;name) for the rank colorizer. Thin derivation of
     * {@link #parseEfetchFull}: entries with a non-Linnaean rank ("no rank"/"clade"/empty) are
     * skipped. Pure -- no I/O.
     */
    static RankedLineage parseEfetchTaxonomyXml( final String xml ) {
        return toRankedLineage( parseEfetchFull( xml ) );
    }

    /**
     * Builds a full {@link ResolvedTaxonomy} (scientific name, rank, tax-id, common name, full name
     * lineage) for the "Fetch Sequence &amp; Taxonomic Data" tool. Thin derivation of
     * {@link #parseEfetchFull}. Pure -- no I/O.
     */
    static ResolvedTaxonomy parseEfetchTaxonomyDetail( final String xml ) {
        return toResolvedTaxonomy( parseEfetchFull( xml ) );
    }

    /**
     * Parses an NCBI efetch taxonomy XML response into the single, raw {@link CachedTaxon} shape that
     * feeds both consumers (and the disk cache): the queried taxon's scientific name, rank, tax-id and
     * common name, plus its ancestors ({@code <LineageEx>} root&rarr;parent) as {@code {name, rank}}
     * pairs (ranks kept verbatim, including "no rank"/"clade"). Pure -- no I/O; never throws.
     */
    static CachedTaxon parseEfetchFull( final String xml ) {
        final Document doc = WsHttp.parseXml( xml );
        if ( doc == null ) {
            return EMPTY_TAXON;
        }
        final Element taxon = WsHttp.firstChildElement( doc.getDocumentElement(), "Taxon" );
        if ( taxon == null ) {
            return EMPTY_TAXON;
        }
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
        final List<String[]> ancestors = new ArrayList<String[]>();
        final Element lineage_ex = WsHttp.firstChildElement( taxon, "LineageEx" );
        if ( lineage_ex != null ) {
            for( final Element anc : WsHttp.childElements( lineage_ex, "Taxon" ) ) {
                final String n = trimOrNull( WsHttp.text( WsHttp.firstChildElement( anc, "ScientificName" ) ) );
                if ( n != null ) {
                    final String r = trimOrNull( WsHttp.text( WsHttp.firstChildElement( anc, "Rank" ) ) );
                    ancestors.add( new String[] { n, ( r == null ) ? "" : r } );
                }
            }
        }
        return new CachedTaxon( tax_id, rank, sci, common, ancestors );
    }

    /** Derives the rank colorizer's {@link RankedLineage} from a {@link CachedTaxon}. Pure. */
    static RankedLineage toRankedLineage( final CachedTaxon ct ) {
        if ( ( ct == null ) || ct.isEmpty() ) {
            return RankedLineage.EMPTY;
        }
        final Map<String, String> rank_to_name = new LinkedHashMap<String, String>();
        for( final String[] a : ct.getAncestors() ) {
            addRankedEntry( rank_to_name, a[ 1 ], a[ 0 ] );
        }
        // the queried taxon itself (its own rank/name -- e.g. genus "Felis" for a Felis query)
        addRankedEntry( rank_to_name, ct.getRank(), ct.getScientificName() );
        return rank_to_name.isEmpty() ? RankedLineage.EMPTY : new RankedLineage( rank_to_name );
    }

    private static void addRankedEntry( final Map<String, String> rank_to_name, final String rank, final String name ) {
        if ( !ForesterUtil.isEmpty( rank ) && !ForesterUtil.isEmpty( name ) && !NO_RANK.equalsIgnoreCase( rank )
                && !CLADE.equalsIgnoreCase( rank ) ) {
            rank_to_name.put( rank, name );
        }
    }

    /** Derives the Fetch tool's {@link ResolvedTaxonomy} from a {@link CachedTaxon}. Pure. */
    static ResolvedTaxonomy toResolvedTaxonomy( final CachedTaxon ct ) {
        if ( ( ct == null ) || ct.isEmpty() ) {
            return ResolvedTaxonomy.EMPTY;
        }
        String rank = ct.getRank();
        if ( ( rank != null ) && ( NO_RANK.equalsIgnoreCase( rank ) || CLADE.equalsIgnoreCase( rank ) ) ) {
            rank = null;
        }
        // full lineage of scientific names (all ancestors, including no-rank ones), then the taxon itself
        final List<String> lineage = new ArrayList<String>();
        for( final String[] a : ct.getAncestors() ) {
            lineage.add( a[ 0 ] );
        }
        if ( ct.getScientificName() != null ) {
            lineage.add( ct.getScientificName() );
        }
        final ResolvedTaxonomy result = new ResolvedTaxonomy( ct.getScientificName(), rank, ct.getTaxId(),
                                                              ct.getCommonName(), lineage );
        return result.isEmpty() ? ResolvedTaxonomy.EMPTY : result;
    }

    private static String trimOrNull( final String s ) {
        if ( s == null ) {
            return null;
        }
        final String t = s.trim();
        return t.isEmpty() ? null : t;
    }
}
