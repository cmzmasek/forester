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
    // null-valued entries double as a negative cache (taxon was queried, nothing found).
    private final Map<String, RankedLineage>     _cache        = Collections
            .synchronizedMap( new HashMap<String, RankedLineage>() );
    private final Map<String, ResolvedTaxonomy>  _detail_cache = Collections
            .synchronizedMap( new HashMap<String, ResolvedTaxonomy>() );

    @Override
    public RankedLineage lineageOf( final String taxon ) {
        if ( ForesterUtil.isEmpty( taxon ) ) {
            return null;
        }
        return _cache.get( key( taxon ) );
    }

    @Override
    public RankedLineage fetch( final String taxon ) throws IOException {
        if ( ForesterUtil.isEmpty( taxon ) ) {
            return RankedLineage.EMPTY;
        }
        final String k = key( taxon );
        synchronized ( _cache ) {
            if ( _cache.containsKey( k ) ) {
                return _cache.get( k );
            }
        }
        final String id = parseEsearchFirstId( WsHttp.httpGet( ESEARCH + WsHttp.encode( taxon ) ) );
        if ( ForesterUtil.isEmpty( id ) ) {
            // esearch found nothing -- a definitive negative; cache it so we never re-query.
            _cache.put( k, RankedLineage.EMPTY );
            return RankedLineage.EMPTY;
        }
        final RankedLineage parsed = parseEfetchTaxonomyXml( WsHttp.httpGet( EFETCH + WsHttp.encode( id ) ) );
        if ( ( parsed != null ) && !parsed.isEmpty() ) {
            _cache.put( k, parsed );
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
        final ResolvedTaxonomy detail = parseEfetchTaxonomyDetail( WsHttp.httpGet( EFETCH + WsHttp.encode( id ) ) );
        if ( ( detail != null ) && !detail.isEmpty() ) {
            _detail_cache.put( k, detail );
            return detail;
        }
        return ResolvedTaxonomy.EMPTY;
    }

    /** True if {@code s} is a non-empty run of digits (an NCBI tax-id). */
    private static boolean isTaxId( final String s ) {
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
     * Builds a {@link RankedLineage} from an NCBI efetch taxonomy XML response: the ranked ancestors
     * in {@code <LineageEx>} (root&rarr;parent) plus the queried taxon's own rank and scientific
     * name. Entries with a non-Linnaean rank ("no rank"/"clade"/empty) are skipped. Pure -- no I/O.
     */
    static RankedLineage parseEfetchTaxonomyXml( final String xml ) {
        final Document doc = WsHttp.parseXml( xml );
        if ( doc == null ) {
            return RankedLineage.EMPTY;
        }
        final Element taxon = WsHttp.firstChildElement( doc.getDocumentElement(), "Taxon" );
        if ( taxon == null ) {
            return RankedLineage.EMPTY;
        }
        final Map<String, String> rank_to_name = new LinkedHashMap<String, String>();
        final Element lineage_ex = WsHttp.firstChildElement( taxon, "LineageEx" );
        if ( lineage_ex != null ) {
            for( final Element anc : WsHttp.childElements( lineage_ex, "Taxon" ) ) {
                addRankedEntry( rank_to_name, anc );
            }
        }
        // the queried taxon itself (its own rank/name -- e.g. genus "Felis" for a Felis query)
        addRankedEntry( rank_to_name, taxon );
        return rank_to_name.isEmpty() ? RankedLineage.EMPTY : new RankedLineage( rank_to_name );
    }

    private static void addRankedEntry( final Map<String, String> rank_to_name, final Element taxon_el ) {
        final String rank = WsHttp.text( WsHttp.firstChildElement( taxon_el, "Rank" ) );
        final String name = WsHttp.text( WsHttp.firstChildElement( taxon_el, "ScientificName" ) );
        if ( !ForesterUtil.isEmpty( rank ) && !ForesterUtil.isEmpty( name ) && !NO_RANK.equalsIgnoreCase( rank )
                && !CLADE.equalsIgnoreCase( rank ) ) {
            rank_to_name.put( rank, name );
        }
    }

    /**
     * Builds a full {@link ResolvedTaxonomy} from an NCBI efetch taxonomy XML response: the queried
     * taxon's scientific name, rank, tax-id and common name, plus the complete lineage of scientific
     * names ({@code <LineageEx>} root&rarr;parent, then the taxon itself). Pure -- no I/O.
     */
    static ResolvedTaxonomy parseEfetchTaxonomyDetail( final String xml ) {
        final Document doc = WsHttp.parseXml( xml );
        if ( doc == null ) {
            return ResolvedTaxonomy.EMPTY;
        }
        final Element taxon = WsHttp.firstChildElement( doc.getDocumentElement(), "Taxon" );
        if ( taxon == null ) {
            return ResolvedTaxonomy.EMPTY;
        }
        final String sci = trimOrNull( WsHttp.text( WsHttp.firstChildElement( taxon, "ScientificName" ) ) );
        final String tax_id = trimOrNull( WsHttp.text( WsHttp.firstChildElement( taxon, "TaxId" ) ) );
        String rank = trimOrNull( WsHttp.text( WsHttp.firstChildElement( taxon, "Rank" ) ) );
        if ( ( rank != null ) && ( NO_RANK.equalsIgnoreCase( rank ) || CLADE.equalsIgnoreCase( rank ) ) ) {
            rank = null;
        }
        // common name: prefer the GenBank common name, else any common name
        String common = null;
        final Element other = WsHttp.firstChildElement( taxon, "OtherNames" );
        if ( other != null ) {
            common = trimOrNull( WsHttp.text( WsHttp.firstChildElement( other, "GenbankCommonName" ) ) );
            if ( common == null ) {
                common = trimOrNull( WsHttp.text( WsHttp.firstChildElement( other, "CommonName" ) ) );
            }
        }
        // full lineage of scientific names (all ancestors, including no-rank ones), then the taxon itself
        final List<String> lineage = new ArrayList<String>();
        final Element lineage_ex = WsHttp.firstChildElement( taxon, "LineageEx" );
        if ( lineage_ex != null ) {
            for( final Element anc : WsHttp.childElements( lineage_ex, "Taxon" ) ) {
                final String n = trimOrNull( WsHttp.text( WsHttp.firstChildElement( anc, "ScientificName" ) ) );
                if ( n != null ) {
                    lineage.add( n );
                }
            }
        }
        if ( sci != null ) {
            lineage.add( sci );
        }
        final ResolvedTaxonomy result = new ResolvedTaxonomy( sci, rank, tax_id, common, lineage );
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
