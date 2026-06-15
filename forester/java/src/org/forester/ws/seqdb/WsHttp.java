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
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.forester.util.ForesterUtil;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.xml.sax.InputSource;

/**
 * Shared, self-contained HTTP + XML plumbing for the web-service clients
 * ({@link NcbiTaxonomyLineageService}, {@code UniProtKbClient}, {@code NcbiSequenceClient}). Every
 * request is <b>timeout bounded</b> (so a stalled socket cannot hang a worker thread forever) and
 * <b>throttled</b> to a polite rate (NCBI E-utilities allow ~3 requests/sec without an API key, or
 * answer HTTP 429). The throttle is process-global on purpose: it serializes all of these clients'
 * requests, which keeps us under NCBI's limit even when several run.
 */
final class WsHttp {

    private static final long   MIN_REQUEST_INTERVAL_MS = 350L;
    private static final int    CONNECT_TIMEOUT_MS      = 10000;
    private static final int    READ_TIMEOUT_MS         = 20000;
    private static final long   RETRY_BACKOFF_MS        = 1500L;
    private static final Object REQUEST_LOCK            = new Object();
    private static long         _last_request_ms        = 0L;

    /** URL-encodes {@code s} (UTF-8); on the (impossible-for-UTF-8) failure returns the trimmed input. */
    static String encode( final String s ) {
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

    /** A throttled, timeout-bounded GET returning the body as a String, with one 429/503 back-off retry. */
    static String httpGet( final String url ) throws IOException {
        throttle();
        try {
            return read( url );
        }
        catch ( final IOException e ) {
            // a rate-limit (429) or transient server error (503): back off once and retry
            final String m = String.valueOf( e.getMessage() );
            if ( m.contains( "429" ) || m.contains( "503" ) ) {
                try {
                    Thread.sleep( RETRY_BACKOFF_MS );
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

    /** A timeout-bounded GET (unlike {@link ForesterUtil#readUrl}, which sets no timeouts). */
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
        // decode as UTF-8 (NCBI/UniProt responses are UTF-8) so accented organism/common names survive
        try ( BufferedReader in = new BufferedReader( new InputStreamReader( c.getInputStream(), StandardCharsets.UTF_8 ) ) ) {
            String line;
            while ( ( line = in.readLine() ) != null ) {
                sb.append( line ).append( '\n' );
            }
        }
        return sb.toString();
    }

    // --- XML helpers (namespace-unaware, external-DTD/entity loading disabled, errors swallowed) ------

    /** Parses {@code xml} into a DOM document, or {@code null} if empty/malformed (never throws). */
    static Document parseXml( final String xml ) {
        if ( ForesterUtil.isEmpty( xml ) ) {
            return null;
        }
        try {
            final DocumentBuilderFactory f = DocumentBuilderFactory.newInstance();
            f.setValidating( false );
            f.setNamespaceAware( false );
            // never reach out for a DTD referenced in the response's DOCTYPE
            setFeatureQuietly( f, "http://apache.org/xml/features/nonvalidating/load-external-dtd", false );
            setFeatureQuietly( f, "http://xml.org/sax/features/external-general-entities", false );
            setFeatureQuietly( f, "http://xml.org/sax/features/external-parameter-entities", false );
            final DocumentBuilder b = f.newDocumentBuilder();
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

    /** The first direct child element named {@code tag} of {@code parent}, or {@code null}. */
    static Element firstChildElement( final Element parent, final String tag ) {
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

    /** All direct child elements named {@code tag} of {@code parent} (in document order). */
    static List<Element> childElements( final Element parent, final String tag ) {
        final List<Element> result = new ArrayList<Element>();
        if ( parent != null ) {
            for( Node n = parent.getFirstChild(); n != null; n = n.getNextSibling() ) {
                if ( ( n.getNodeType() == Node.ELEMENT_NODE ) && tag.equals( n.getNodeName() ) ) {
                    result.add( (Element) n );
                }
            }
        }
        return result;
    }

    /** The text content of {@code e}, or {@code null} if {@code e} is null. */
    static String text( final Element e ) {
        return ( e == null ) ? null : e.getTextContent();
    }

    private WsHttp() {
        // not instantiable
    }
}
