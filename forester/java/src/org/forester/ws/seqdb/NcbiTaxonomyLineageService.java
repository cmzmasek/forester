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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.forester.util.ForesterUtil;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

/**
 * A clean {@link TaxonomicLineageService} backed by NCBI Taxonomy via the E-utilities REST API:
 * one <i>esearch</i> maps a scientific name to a tax-id, one <i>efetch</i> returns the full
 * <b>ranked</b> lineage in a single XML response (so we never resolve ranks one ancestor at a time).
 * Results -- including "not found" -- are cached for the life of the process.
 *
 * <p>This is intentionally new, self-contained code: it reuses only generic plumbing
 * ({@link ForesterUtil#readUrl(String)} and the JDK XML parser), not the legacy taxonomy classes.
 * The XML parsers ({@link #parseEfetchTaxonomyXml(String)} / {@link #parseEsearchFirstId(String)})
 * are split out from the network call so they can be unit-tested from captured fixtures.
 */
public final class NcbiTaxonomyLineageService implements TaxonomicLineageService {

    private static final String                ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&tool=Archaeopteryx&term=";
    private static final String                EFETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&tool=Archaeopteryx&retmode=xml&id=";
    // NCBI ranks that are not Linnaean levels a colorization can key on; skipped when building a lineage.
    private static final String                NO_RANK             = "no rank";
    private static final String                CLADE               = "clade";
    // NCBI E-utilities allow ~3 requests/sec without an API key; space requests to stay under that
    // (each fetch makes two: esearch + efetch), or the server answers HTTP 429.
    private static final long                  MIN_REQUEST_INTERVAL_MS = 350L;
    private static final int                   CONNECT_TIMEOUT_MS  = 10000;
    private static final int                   READ_TIMEOUT_MS     = 20000;
    private static final Object                REQUEST_LOCK        = new Object();
    private static long                        _last_request_ms    = 0L;
    // null-valued entries double as a negative cache (taxon was queried, nothing found).
    private final Map<String, RankedLineage>   _cache  = Collections.synchronizedMap( new HashMap<String, RankedLineage>() );

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
        final String id = parseEsearchFirstId( httpGet( ESEARCH + encode( taxon ) ) );
        if ( ForesterUtil.isEmpty( id ) ) {
            // esearch found nothing -- a definitive negative; cache it so we never re-query.
            _cache.put( k, RankedLineage.EMPTY );
            return RankedLineage.EMPTY;
        }
        final RankedLineage parsed = parseEfetchTaxonomyXml( httpGet( EFETCH + encode( id ) ) );
        if ( ( parsed != null ) && !parsed.isEmpty() ) {
            _cache.put( k, parsed );
            return parsed;
        }
        // esearch found the taxon but efetch returned no usable ranked lineage -- possibly a
        // truncated/transient response. Do NOT cache this as a negative, so a later attempt can
        // retry (caching EMPTY here would strand a resolvable taxon for the life of the process).
        return RankedLineage.EMPTY;
    }

    private static String key( final String taxon ) {
        return taxon.trim().toLowerCase( Locale.ROOT );
    }

    private static String encode( final String s ) {
        try {
            return URLEncoder.encode( s.trim(), "UTF-8" );
        }
        catch ( final Exception e ) {
            return s.trim();
        }
    }

    /** Milliseconds to wait before a request to keep at least {@code interval_ms} between requests. Pure. */
    static long throttleDelayMs( final long last_request_ms, final long now_ms, final long interval_ms ) {
        final long wait = ( last_request_ms + interval_ms ) - now_ms;
        return ( wait > 0 ) ? wait : 0;
    }

    /** Blocks just long enough to honor {@link #MIN_REQUEST_INTERVAL_MS} between successive requests. */
    private static void throttle() {
        synchronized ( REQUEST_LOCK ) {
            final long wait = throttleDelayMs( _last_request_ms, System.currentTimeMillis(), MIN_REQUEST_INTERVAL_MS );
            if ( wait > 0 ) {
                try {
                    Thread.sleep( wait );
                }
                catch ( final InterruptedException e ) {
                    Thread.currentThread().interrupt();
                }
            }
            _last_request_ms = System.currentTimeMillis();
        }
    }

    private static String httpGet( final String url ) throws IOException {
        throttle();
        try {
            return read( url );
        }
        catch ( final IOException e ) {
            // a rate-limit (429) or transient server error (503): back off once and retry
            final String m = String.valueOf( e.getMessage() );
            if ( m.contains( "429" ) || m.contains( "503" ) ) {
                try {
                    Thread.sleep( 1500L );
                }
                catch ( final InterruptedException ie ) {
                    Thread.currentThread().interrupt();
                }
                throttle(); // keep the retry within the request-spacing budget too
                return read( url );
            }
            throw e;
        }
    }

    /**
     * A timeout-bounded GET (unlike {@link ForesterUtil#readUrl}, which sets no timeouts) so a stalled
     * NCBI socket cannot hang the background resolver thread forever.
     */
    private static String read( final String url_str ) throws IOException {
        final URLConnection c;
        try {
            c = new URI( url_str ).toURL().openConnection();
        }
        catch ( final URISyntaxException e ) {
            throw new IOException( "bad URL: " + url_str, e );
        }
        c.setConnectTimeout( CONNECT_TIMEOUT_MS );
        c.setReadTimeout( READ_TIMEOUT_MS );
        final StringBuilder sb = new StringBuilder();
        try ( BufferedReader in = new BufferedReader( new InputStreamReader( c.getInputStream() ) ) ) {
            String line;
            while ( ( line = in.readLine() ) != null ) {
                sb.append( line ).append( '\n' );
            }
        }
        return sb.toString();
    }

    /** The first {@code <Id>} under {@code <IdList>} in an esearch response, or {@code null}. Pure. */
    static String parseEsearchFirstId( final String xml ) {
        final Document doc = parse( xml );
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
     * Returns {@link RankedLineage#EMPTY} when nothing usable is present.
     */
    static RankedLineage parseEfetchTaxonomyXml( final String xml ) {
        final Document doc = parse( xml );
        if ( doc == null ) {
            return RankedLineage.EMPTY;
        }
        final Element taxon = firstChildElement( doc.getDocumentElement(), "Taxon" );
        if ( taxon == null ) {
            return RankedLineage.EMPTY;
        }
        final Map<String, String> rank_to_name = new LinkedHashMap<String, String>();
        final Element lineage_ex = firstChildElement( taxon, "LineageEx" );
        if ( lineage_ex != null ) {
            for( final Element anc : childElements( lineage_ex, "Taxon" ) ) {
                addRankedEntry( rank_to_name, anc );
            }
        }
        // the queried taxon itself (its own rank/name -- e.g. genus "Felis" for a Felis query)
        addRankedEntry( rank_to_name, taxon );
        return rank_to_name.isEmpty() ? RankedLineage.EMPTY : new RankedLineage( rank_to_name );
    }

    private static void addRankedEntry( final Map<String, String> rank_to_name, final Element taxon_el ) {
        final String rank = text( firstChildElement( taxon_el, "Rank" ) );
        final String name = text( firstChildElement( taxon_el, "ScientificName" ) );
        if ( !ForesterUtil.isEmpty( rank ) && !ForesterUtil.isEmpty( name ) && !NO_RANK.equalsIgnoreCase( rank )
                && !CLADE.equalsIgnoreCase( rank ) ) {
            rank_to_name.put( rank, name );
        }
    }

    private static Document parse( final String xml ) {
        if ( ForesterUtil.isEmpty( xml ) ) {
            return null;
        }
        try {
            final DocumentBuilderFactory f = DocumentBuilderFactory.newInstance();
            f.setValidating( false );
            f.setNamespaceAware( false );
            // never reach out for the NCBI DTD referenced in the response's DOCTYPE
            setFeatureQuietly( f, "http://apache.org/xml/features/nonvalidating/load-external-dtd", false );
            setFeatureQuietly( f, "http://xml.org/sax/features/external-general-entities", false );
            setFeatureQuietly( f, "http://xml.org/sax/features/external-parameter-entities", false );
            final DocumentBuilder b = f.newDocumentBuilder();
            // swallow parser warnings/errors -- a malformed response just yields a null document here
            b.setErrorHandler( new org.xml.sax.helpers.DefaultHandler() );
            return b.parse( new InputSource( new StringReader( xml ) ) );
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    private static void setFeatureQuietly( final DocumentBuilderFactory f, final String feature, final boolean value ) {
        try {
            f.setFeature( feature, value );
        }
        catch ( final Exception e ) {
            // feature unsupported by this parser -- best effort
        }
    }

    private static Element firstChildElement( final Element parent, final String tag ) {
        if ( parent == null ) {
            return null;
        }
        for( Node n = parent.getFirstChild(); n != null; n = n.getNextSibling() ) {
            if ( ( n.getNodeType() == Node.ELEMENT_NODE ) && tag.equals( n.getNodeName() ) ) {
                return (Element) n;
            }
        }
        return null;
    }

    private static java.util.List<Element> childElements( final Element parent, final String tag ) {
        final java.util.List<Element> result = new java.util.ArrayList<Element>();
        if ( parent != null ) {
            for( Node n = parent.getFirstChild(); n != null; n = n.getNextSibling() ) {
                if ( ( n.getNodeType() == Node.ELEMENT_NODE ) && tag.equals( n.getNodeName() ) ) {
                    result.add( (Element) n );
                }
            }
        }
        return result;
    }

    private static String text( final Element e ) {
        return ( e == null ) ? null : e.getTextContent();
    }
}
