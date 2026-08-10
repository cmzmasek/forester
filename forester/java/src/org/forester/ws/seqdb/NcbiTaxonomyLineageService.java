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

import org.forester.util.ForesterUtil;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * A clean {@link TaxonomicLineageService} backed by NCBI Taxonomy via the E-utilities REST API:
 * one <i>esearch</i> maps a scientific name to a tax-id, one <i>efetch</i> returns the full
 * lineage in a single XML response (so we never resolve ranks one ancestor at a time). It serves every
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
    // EMPTY entries double as a negative cache (taxon was queried, nothing found).
    private final Map<String, TaxonLineage>      _cache        = Collections
            .synchronizedMap( new HashMap<String, TaxonLineage>() );
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
            _cache.put( k, TaxonLineage.EMPTY );
            return TaxonLineage.EMPTY;
        }
        final TaxonLineage parsed = parseEfetchFull( WsHttp.httpGet( EFETCH + WsHttp.encode( id ) ) );
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
     * Parses an NCBI efetch taxonomy XML response into a {@link TaxonLineage}: the queried taxon's
     * scientific name, rank, tax-id and common name, plus its ancestors ({@code <LineageEx>}
     * root&rarr;parent) each as {@code {name, rank, tax-id}} (ranks kept verbatim, including "no
     * rank"/"clade"; the ranked lookup maps are derived inside {@link TaxonLineage}). Pure -- no I/O;
     * never throws.
     */
    static TaxonLineage parseEfetchFull( final String xml ) {
        final Document doc = WsHttp.parseXml( xml );
        if ( doc == null ) {
            return TaxonLineage.EMPTY;
        }
        final Element taxon = WsHttp.firstChildElement( doc.getDocumentElement(), "Taxon" );
        if ( taxon == null ) {
            return TaxonLineage.EMPTY;
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
