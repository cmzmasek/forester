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

package org.forester.archaeopteryx;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URI;
import java.nio.file.Files;
import java.security.MessageDigest;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.imageio.ImageIO;
import javax.swing.SwingUtilities;

import org.forester.util.ForesterUtil;

/**
 * Loads and caches the tip images (see {@link TipImages}). Two sources, auto-detected from the reference: a LOCAL
 * file (resolved against the loaded tree file's directory, or an absolute path) or an http(s) URL (fetched once and
 * kept in a small on-disk cache, so re-opening a tree of web images is instant and offline). A reference that fails
 * to load is remembered so it is not retried every repaint.
 * <p>
 * {@link #get(String, File)} is the paint-time entry point: it returns the decoded image if already in memory, or
 * {@code null} while it schedules a background load (off the EDT). When a load finishes it runs the repaint callback
 * on the EDT, so the freshly-loaded image is drawn on the next frame. The in-memory cache is bounded (LRU).
 */
final class TipImageCache {

    /** A negative-result marker kept in {@link #_cache} so a broken reference is not reloaded on every repaint. */
    private static final BufferedImage FAILED           = new BufferedImage( 1, 1, BufferedImage.TYPE_INT_ARGB );

    private static final int           DEFAULT_MAX      = 300;    // in-memory images (LRU)
    private static final int           CONNECT_TIMEOUT  = 10_000; // ms
    private static final int           READ_TIMEOUT     = 20_000; // ms
    private static final int           MAX_IMAGE_BYTES  = 25 * 1024 * 1024; // reject a runaway/huge download
    private static final int           MAX_REDIRECTS    = 5;      // Location followed manually (incl. http<->https)
    private static final String        USER_AGENT       = "Archaeopteryx/" + AptxConstants.VERSION
            + " (https://github.com/cmzmasek/archaeopteryx)";

    /** ONE shared daemon pool for every cache instance (a per-tab cache would otherwise leak its threads on tab churn). */
    private static final ExecutorService LOADER = Executors.newFixedThreadPool( 4, r -> {
                                                     final Thread t = new Thread( r, "tip-image-loader" );
                                                     t.setDaemon( true ); // never block JVM exit
                                                     return t;
                                                 } );

    /** access-ordered LRU: the eldest (least-recently drawn) image is evicted past the cap. Guarded by {@code this}. */
    private final Map<String, BufferedImage> _cache;
    private final Set<String>                _inflight = ConcurrentHashMap.newKeySet();
    private volatile Runnable                _on_loaded;

    TipImageCache() {
        this( DEFAULT_MAX );
    }

    TipImageCache( final int max_entries ) {
        final int cap = Math.max( 1, max_entries );
        _cache = new LinkedHashMap<String, BufferedImage>( 16, 0.75f, true ) {

            @Override
            protected boolean removeEldestEntry( final Map.Entry<String, BufferedImage> eldest ) {
                return size() > cap;
            }
        };
    }

    /** Sets the callback run on the EDT after a background load completes (so the tree repaints with the new image). */
    void setRepaintCallback( final Runnable on_loaded ) {
        _on_loaded = on_loaded;
    }

    /**
     * The decoded image for {@code ref} (a local path resolved against {@code base_dir}, or an http(s) URL), or
     * {@code null} if it is not yet loaded (a background load is scheduled) or the reference is broken. Call from the
     * EDT paint path.
     */
    BufferedImage get( final String ref, final File base_dir ) {
        if ( ( ref == null ) || ref.isEmpty() ) {
            return null;
        }
        final String key = resolveKey( ref, base_dir );
        synchronized ( this ) {
            final BufferedImage cached = _cache.get( key );
            if ( cached != null ) {
                return ( cached == FAILED ) ? null : cached;
            }
        }
        if ( _inflight.add( key ) ) { // schedule exactly one load per key
            LOADER.submit( () -> {
                BufferedImage img = null;
                try {
                    img = loadImage( ref, base_dir );
                }
                catch ( final Throwable ignored ) {
                    img = null; // any failure -> negative-cache it below
                }
                try {
                    final BufferedImage result = ( img == null ) ? FAILED : img;
                    synchronized ( TipImageCache.this ) {
                        _cache.put( key, result );
                    }
                }
                finally {
                    _inflight.remove( key ); // even if the put throws (e.g. OOM), never strand the key in-flight forever
                }
                final Runnable cb = _on_loaded;
                if ( cb != null ) {
                    SwingUtilities.invokeLater( cb );
                }
            } );
        }
        return null;
    }

    /** Whether {@code ref} has been tried and failed to load (so the caller can draw a small placeholder). */
    boolean isFailed( final String ref, final File base_dir ) {
        if ( ( ref == null ) || ref.isEmpty() ) {
            return false;
        }
        synchronized ( this ) {
            return _cache.get( resolveKey( ref, base_dir ) ) == FAILED;
        }
    }

    /** Drops all cached images (on tree change / for tests). Does not clear the on-disk URL cache. */
    synchronized void clear() {
        _cache.clear();
    }

    /** The in-memory cache key: the URL as-is, or a local path made absolute against {@code base_dir}. */
    private static String resolveKey( final String ref, final File base_dir ) {
        if ( TipImages.isUrl( ref ) ) {
            return ref.trim();
        }
        return localFile( ref, base_dir ).getAbsolutePath();
    }

    private static File localFile( final String ref, final File base_dir ) {
        final File f = new File( ref.trim() );
        if ( f.isAbsolute() || ( base_dir == null ) ) {
            return f;
        }
        return new File( base_dir, ref.trim() );
    }

    /**
     * Synchronously decodes {@code ref}: a local file (resolved against {@code base_dir}) via {@code ImageIO}, or an
     * http(s) URL fetched through the on-disk URL cache. Returns null on any failure (missing file, network error,
     * unsupported/undecodable content). Package-visible so tests can exercise the load logic directly.
     */
    static BufferedImage loadImage( final String ref, final File base_dir ) throws IOException {
        if ( ( ref == null ) || ref.trim().isEmpty() ) {
            return null;
        }
        if ( TipImages.isUrl( ref ) ) {
            return loadUrl( ref.trim() );
        }
        final File f = localFile( ref, base_dir );
        if ( !f.isFile() ) {
            return null;
        }
        return ImageIO.read( f ); // null if the bytes are not a decodable image
    }

    /**
     * Loads a URL image via the on-disk cache: a cache hit reads the file; a miss fetches, DECODES, and only then (for a
     * genuine image) stores it. Decoding before writing keeps HTML error pages (served 200) and truncated bodies out of
     * the disk cache, so a bad response can never permanently poison it.
     */
    private static BufferedImage loadUrl( final String url ) throws IOException {
        final File cached = urlCacheFile( url );
        if ( ( cached != null ) && cached.isFile() ) {
            final BufferedImage img = readImageQuietly( cached ); // null if it returns null OR throws (corrupt/partial)
            if ( img != null ) {
                return img;
            }
            cached.delete(); // a poisoned cache file -> drop it and refetch
        }
        final byte[] bytes = fetch( url );
        if ( ( bytes == null ) || ( bytes.length == 0 ) ) {
            return null;
        }
        final BufferedImage img = ImageIO.read( new java.io.ByteArrayInputStream( bytes ) );
        if ( img == null ) {
            return null; // a non-image body (e.g. an HTML page served 200) -> do not cache the junk
        }
        if ( cached != null ) {
            writeQuietly( cached, bytes ); // only a decodable image reaches the disk cache
        }
        return img;
    }

    /** {@link ImageIO#read(File)} that maps a THROW (a truncated/corrupt file makes it throw, not return null) to null,
     *  so the caller treats it the same as an undecodable file: delete + refetch, never permanently poisoned. */
    private static BufferedImage readImageQuietly( final File f ) {
        try {
            return ImageIO.read( f );
        }
        catch ( final Throwable t ) {
            return null;
        }
    }

    /** Fetches the bytes at {@code url}, following Location redirects MANUALLY (so an http-&gt;https cross-protocol
     *  redirect works -- {@link HttpURLConnection}'s auto-follow refuses those), bounded by {@link #MAX_IMAGE_BYTES}. */
    private static byte[] fetch( final String url ) throws IOException {
        String current = url;
        for( int hop = 0; hop <= MAX_REDIRECTS; ++hop ) {
            final HttpURLConnection c = openConnection( current );
            try {
                final int code = c.getResponseCode();
                if ( ( code >= 300 ) && ( code < 400 ) ) { // a redirect we follow ourselves
                    final String location = c.getHeaderField( "Location" );
                    if ( ForesterUtil.isEmpty( location ) ) {
                        return null;
                    }
                    current = resolveRedirect( current, location );
                    continue;
                }
                if ( code != HttpURLConnection.HTTP_OK ) {
                    return null;
                }
                try ( final InputStream in = c.getInputStream() ) {
                    final java.io.ByteArrayOutputStream out = new java.io.ByteArrayOutputStream();
                    final byte[] buf = new byte[ 8192 ];
                    long total = 0;
                    int n;
                    while ( ( n = in.read( buf ) ) > 0 ) {
                        total += n;
                        if ( total > MAX_IMAGE_BYTES ) {
                            return null; // runaway/huge download -> give up rather than buffer it all in memory
                        }
                        out.write( buf, 0, n );
                    }
                    return out.toByteArray();
                }
            }
            finally {
                c.disconnect();
            }
        }
        return null; // too many redirects
    }

    private static HttpURLConnection openConnection( final String url ) throws IOException {
        final HttpURLConnection c = (HttpURLConnection) toUrl( url ).openConnection();
        c.setConnectTimeout( CONNECT_TIMEOUT );
        c.setReadTimeout( READ_TIMEOUT );
        // Wikimedia (the flagship Special:FilePath source) rejects a bare/generic User-Agent, so send a descriptive one.
        c.setRequestProperty( "User-Agent", USER_AGENT );
        c.setInstanceFollowRedirects( false ); // we follow Location ourselves so http<->https redirects work
        return c;
    }

    /** Lenient URL parsing: {@link URI#create} is RFC-strict and throws on a raw space (common in a Special:FilePath
     *  file name), so fall back to percent-encoding spaces before giving up. */
    private static java.net.URL toUrl( final String url ) throws IOException {
        try {
            return URI.create( url ).toURL();
        }
        catch ( final Exception e ) {
            try {
                return URI.create( url.replace( " ", "%20" ) ).toURL();
            }
            catch ( final Exception e2 ) {
                throw new IOException( "malformed image URL: " + url );
            }
        }
    }

    private static String resolveRedirect( final String base, final String location ) {
        if ( location.startsWith( "http://" ) || location.startsWith( "https://" ) ) {
            return location; // absolute (what real servers send)
        }
        try {
            return URI.create( base.replace( " ", "%20" ) ).resolve( location.replace( " ", "%20" ) ).toString();
        }
        catch ( final Exception e ) {
            return location; // best effort for a relative Location
        }
    }

    /** The on-disk cache file for a URL: {@code <cache-dir>/image-cache/<md5(url)>.img}, or null if unavailable. */
    private static File urlCacheFile( final String url ) {
        try {
            final File dir = new File( cacheDir(), "image-cache" );
            if ( !dir.isDirectory() && !dir.mkdirs() ) {
                return null;
            }
            return new File( dir, md5Hex( url ) + ".img" );
        }
        catch ( final Throwable t ) {
            return null; // best-effort: no cache dir -> just refetch each session
        }
    }

    private static File cacheDir() {
        final String override = System.getProperty( "archaeopteryx.cache.dir" );
        if ( ( override != null ) && !override.isEmpty() ) {
            return new File( override );
        }
        return new File( System.getProperty( "user.home", "." ), ".archaeopteryx" );
    }

    private static String md5Hex( final String s ) throws Exception {
        final byte[] d = MessageDigest.getInstance( "MD5" ).digest( s.getBytes( "UTF-8" ) );
        final StringBuilder sb = new StringBuilder( d.length * 2 );
        for( final byte b : d ) {
            sb.append( Character.forDigit( ( b >> 4 ) & 0xF, 16 ) ).append( Character.forDigit( b & 0xF, 16 ) );
        }
        return sb.toString().toLowerCase( Locale.ROOT );
    }

    private static void writeQuietly( final File f, final byte[] bytes ) {
        File tmp = null;
        try {
            tmp = File.createTempFile( "img", ".tmp", f.getParentFile() );
            try ( final OutputStream os = Files.newOutputStream( tmp.toPath() ) ) {
                os.write( bytes );
            }
            if ( tmp.renameTo( f ) ) { // atomic-ish; if the target already exists on some FS, just drop the temp
                tmp = null; // renamed into place -> nothing to clean up
            }
        }
        catch ( final Throwable ignored ) {
            // best-effort caching only
        }
        finally {
            if ( tmp != null ) {
                tmp.delete(); // a failed write/rename must not orphan the temp file in the cache dir
            }
        }
    }
}
