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

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

/**
 * A small, best-effort persistent cache for resolved NCBI taxonomy ({@link CachedTaxon}). It lets the
 * rank colorizer and the "Fetch Sequence &amp; Taxonomic Data" tool remember taxa <i>between</i>
 * sessions, so re-opening trees of an organism group you have already seen does not re-query NCBI.
 *
 * <p>Design priorities, in order:
 * <ol>
 *   <li><b>Never a dependency.</b> Every disk operation is wrapped; on any failure the cache flips to
 *       {@link #status() unavailable} and all methods become no-ops &mdash; Archaeopteryx then behaves
 *       exactly as before (a miss just means "go fetch"), only slower on the affected lookups.</li>
 *   <li><b>Dumb format.</b> One line per entry, tab-separated, positionally parsed (same discipline as
 *       the UniProt TSV client) &mdash; no JSON library, no database. The file is append-only, which
 *       makes concurrent appends from two Archaeopteryx instances safe; {@link #load()} compacts it
 *       (last-wins per key, expired entries dropped), bounding growth without an eviction policy.</li>
 *   <li><b>Positive only.</b> "Not found" is never persisted; a transient miss must not become sticky
 *       for the whole TTL.</li>
 * </ol>
 *
 * <p>Location: {@code ${user.home}/.archaeopteryx/taxonomy-cache.tsv}, overridable for tests via the
 * {@value #DIR_PROPERTY} system property. The user toggle is persisted as the presence of a sibling
 * {@code .disabled} marker file. Not instantiated until first use; the constructor does no I/O.
 */
final class TaxonomyDiskCache {

    static final String         DIR_PROPERTY  = "archaeopteryx.cache.dir";
    static final String         DEFAULT_DIR   = ".archaeopteryx";
    static final String         CACHE_FILE    = "taxonomy-cache.tsv";
    static final String         DISABLED_FILE = "taxonomy-cache.disabled";
    static final long           TTL_MS        = TimeUnit.DAYS.toMillis( 30 );
    private static final int    FIELD_COUNT   = 7;

    private final Path          _dir;
    private final Path          _file;
    private final Path          _disabled_marker;
    private boolean             _enabled;
    private boolean             _enabled_known;
    private boolean             _available    = true;
    private String              _reason;

    TaxonomyDiskCache() {
        this( resolveDir() );
    }

    /** For tests: an explicit cache directory. The constructor performs no I/O. */
    TaxonomyDiskCache( final Path dir ) {
        _dir = dir;
        _file = dir.resolve( CACHE_FILE );
        _disabled_marker = dir.resolve( DISABLED_FILE );
    }

    private static Path resolveDir() {
        final String override = System.getProperty( DIR_PROPERTY );
        if ( ( override != null ) && !override.trim().isEmpty() ) {
            return Paths.get( override.trim() );
        }
        return Paths.get( System.getProperty( "user.home", "." ), DEFAULT_DIR );
    }

    // ---- enable toggle (persisted as a marker file) --------------------------------------------

    synchronized boolean isEnabled() {
        if ( !_enabled_known ) {
            boolean disabled = false;
            try {
                disabled = Files.exists( _disabled_marker );
            }
            catch ( final Exception e ) {
                // can't tell -> default to enabled
            }
            _enabled = !disabled;
            _enabled_known = true;
        }
        return _enabled;
    }

    synchronized void setEnabled( final boolean enabled ) {
        _enabled = enabled;
        _enabled_known = true;
        try {
            if ( enabled ) {
                Files.deleteIfExists( _disabled_marker );
            }
            else {
                ensureDir();
                Files.write( _disabled_marker, new byte[ 0 ] );
            }
        }
        catch ( final Exception e ) {
            // best-effort: the toggle still holds for this session even if the marker can't be written
        }
    }

    // ---- core operations (all best-effort) -----------------------------------------------------

    /**
     * Reads the cache, dropping expired entries and keeping the last line per key, then compacts the
     * file to that surviving set. Returns the surviving {@code key -> CachedTaxon} map (empty if the
     * cache is disabled or unavailable). Does disk I/O -- call off the EDT where possible.
     */
    synchronized Map<String, CachedTaxon> load() {
        final Map<String, CachedTaxon> result = new LinkedHashMap<String, CachedTaxon>();
        if ( !isEnabled() ) {
            return result;
        }
        try {
            if ( !Files.exists( _file ) ) { // a pure read never creates the cache directory
                return result;
            }
            final List<String> lines = Files.readAllLines( _file, StandardCharsets.UTF_8 );
            final long now = System.currentTimeMillis();
            final Map<String, String> surviving = new LinkedHashMap<String, String>();
            for( final String line : lines ) {
                final Parsed p = parseLine( line );
                if ( ( p == null ) || ( ( now - p._epoch_ms ) > TTL_MS ) ) {
                    continue;
                }
                surviving.put( p._key, line ); // last line for a key wins
                result.put( p._key, p._taxon );
            }
            if ( surviving.size() != lines.size() ) {
                compact( surviving.values() ); // only rewrite when it would actually shrink the file
            }
            _available = true;
            _reason = null;
        }
        catch ( final Exception e ) {
            markUnavailable( e );
        }
        return result;
    }

    /** Appends one positive entry. No-op when disabled, unavailable, or the record is empty. */
    synchronized void put( final String key, final CachedTaxon taxon ) {
        if ( !isEnabled() || ( key == null ) || ( taxon == null ) || taxon.isEmpty() || !probe() ) {
            return;
        }
        try {
            Files.write( _file, Collections.singletonList( toLine( key, taxon, System.currentTimeMillis() ) ),
                         StandardCharsets.UTF_8, StandardOpenOption.CREATE, StandardOpenOption.APPEND );
        }
        catch ( final Exception e ) {
            markUnavailable( e );
        }
    }

    /** Deletes the cache file (the in-memory maps are cleared by the caller). */
    synchronized void clear() {
        try {
            Files.deleteIfExists( _file );
        }
        catch ( final Exception e ) {
            markUnavailable( e );
        }
    }

    synchronized TaxonomyCacheStatus status() {
        final boolean enabled = isEnabled();
        boolean available = checkAvailableWithoutCreating();
        long bytes = 0L;
        int entries = 0;
        long oldest = 0L;
        if ( available && Files.exists( _file ) ) {
            try {
                bytes = Files.size( _file );
                final long now = System.currentTimeMillis();
                final Map<String, Long> by_key = new LinkedHashMap<String, Long>();
                for( final String line : Files.readAllLines( _file, StandardCharsets.UTF_8 ) ) {
                    final Parsed p = parseLine( line );
                    if ( ( p == null ) || ( ( now - p._epoch_ms ) > TTL_MS ) ) {
                        continue;
                    }
                    by_key.put( p._key, Long.valueOf( p._epoch_ms ) ); // distinct keys; last timestamp wins
                }
                entries = by_key.size();
                for( final Long ts : by_key.values() ) {
                    if ( ( oldest == 0L ) || ( ts.longValue() < oldest ) ) {
                        oldest = ts.longValue();
                    }
                }
            }
            catch ( final Exception e ) {
                markUnavailable( e );
                available = false;
            }
        }
        return new TaxonomyCacheStatus( available, enabled, _file.toString(), bytes, entries, oldest,
                                        available ? null : _reason );
    }

    // ---- filesystem helpers --------------------------------------------------------------------

    /** Can we write here (creating the directory if needed) -- the check used before an actual write. */
    private boolean probe() {
        try {
            ensureDir();
            if ( !Files.isWritable( _dir ) ) {
                _available = false;
                _reason = "cache directory is not writable: " + _dir;
                return false;
            }
            _available = true;
            _reason = null;
            return true;
        }
        catch ( final Exception e ) {
            markUnavailable( e );
            return false;
        }
    }

    /** Availability check for {@link #status()} that never creates anything (so reads have no side effects). */
    private boolean checkAvailableWithoutCreating() {
        try {
            Path d = _dir;
            while ( ( d != null ) && !Files.exists( d ) ) {
                d = d.getParent(); // the cache dir may not exist yet; a writable ancestor means we could create it
            }
            if ( ( d == null ) || !Files.isDirectory( d ) || !Files.isWritable( d ) ) {
                _reason = "cache directory is not writable: " + _dir;
                return false;
            }
            _reason = null;
            return true;
        }
        catch ( final Exception e ) {
            markUnavailable( e );
            return false;
        }
    }

    private void ensureDir() throws java.io.IOException {
        if ( !Files.isDirectory( _dir ) ) {
            Files.createDirectories( _dir );
        }
    }

    private void markUnavailable( final Exception e ) {
        _available = false;
        _reason = ( e.getMessage() == null ) ? e.getClass().getSimpleName() : e.getMessage();
    }

    private void compact( final Collection<String> lines ) {
        Path tmp = null;
        try {
            tmp = Files.createTempFile( _dir, "taxcache", ".tmp" );
            Files.write( tmp, lines, StandardCharsets.UTF_8 );
            try {
                Files.move( tmp, _file, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE );
            }
            catch ( final Exception atomic_unsupported ) {
                Files.move( tmp, _file, StandardCopyOption.REPLACE_EXISTING );
            }
            tmp = null; // moved into place -- nothing left to clean up
        }
        catch ( final Exception e ) {
            // best-effort: leaving the un-compacted file in place is fine -- it still loads next time
        }
        finally {
            if ( tmp != null ) {
                try {
                    Files.deleteIfExists( tmp ); // don't orphan the temp file if write/move failed
                }
                catch ( final Exception e ) {
                    // best-effort cleanup
                }
            }
        }
    }

    // ---- serialization (pure; package-static for unit testing) ---------------------------------

    static String toLine( final String key, final CachedTaxon t, final long epoch_ms ) {
        final StringBuilder anc = new StringBuilder();
        final List<String[]> ancestors = t.getAncestors();
        for( int i = 0; i < ancestors.size(); ++i ) {
            if ( i > 0 ) {
                anc.append( '|' );
            }
            anc.append( esc( ancestors.get( i )[ 0 ] ) ).append( '~' ).append( esc( ancestors.get( i )[ 1 ] ) );
        }
        return esc( key ) + "\t" + epoch_ms + "\t" + esc( t.getTaxId() ) + "\t" + esc( t.getRank() ) + "\t"
                + esc( t.getScientificName() ) + "\t" + esc( t.getCommonName() ) + "\t" + anc;
    }

    static Parsed parseLine( final String line ) {
        if ( line == null ) {
            return null;
        }
        final String[] f = line.split( "\t", -1 );
        if ( f.length < FIELD_COUNT ) {
            return null;
        }
        final long epoch;
        try {
            epoch = Long.parseLong( f[ 1 ].trim() );
        }
        catch ( final NumberFormatException e ) {
            return null;
        }
        final String key = unesc( f[ 0 ] );
        if ( key.isEmpty() ) {
            return null;
        }
        final List<String[]> ancestors = new ArrayList<String[]>();
        if ( !f[ 6 ].isEmpty() ) {
            for( final String pair : f[ 6 ].split( "\\|", -1 ) ) {
                final int sep = pair.indexOf( '~' ); // the only literal '~' is the separator (real ones are escaped)
                if ( sep < 0 ) {
                    ancestors.add( new String[] { unesc( pair ), "" } );
                }
                else {
                    ancestors.add( new String[] { unesc( pair.substring( 0, sep ) ), unesc( pair.substring( sep + 1 ) ) } );
                }
            }
        }
        final CachedTaxon taxon = new CachedTaxon( emptyToNull( unesc( f[ 2 ] ) ), emptyToNull( unesc( f[ 3 ] ) ),
                                                   emptyToNull( unesc( f[ 4 ] ) ), emptyToNull( unesc( f[ 5 ] ) ),
                                                   ancestors );
        return new Parsed( key, epoch, taxon );
    }

    /** Escapes the field/sub-delimiter characters (backslash first) so a value is safe on one line. */
    static String esc( final String s ) {
        if ( s == null ) {
            return "";
        }
        final StringBuilder b = new StringBuilder( s.length() );
        for( int i = 0; i < s.length(); ++i ) {
            final char c = s.charAt( i );
            switch ( c ) {
                case '\\':
                    b.append( "\\\\" );
                    break;
                case '\t':
                    b.append( "\\t" );
                    break;
                case '\n':
                    b.append( "\\n" );
                    break;
                case '\r':
                    b.append( "\\r" );
                    break;
                case '|':
                    b.append( "\\p" );
                    break;
                case '~':
                    b.append( "\\s" );
                    break;
                default:
                    b.append( c );
            }
        }
        return b.toString();
    }

    static String unesc( final String s ) {
        if ( s == null ) {
            return "";
        }
        final StringBuilder b = new StringBuilder( s.length() );
        for( int i = 0; i < s.length(); ++i ) {
            final char c = s.charAt( i );
            if ( ( c == '\\' ) && ( ( i + 1 ) < s.length() ) ) {
                final char n = s.charAt( ++i );
                switch ( n ) {
                    case '\\':
                        b.append( '\\' );
                        break;
                    case 't':
                        b.append( '\t' );
                        break;
                    case 'n':
                        b.append( '\n' );
                        break;
                    case 'r':
                        b.append( '\r' );
                        break;
                    case 'p':
                        b.append( '|' );
                        break;
                    case 's':
                        b.append( '~' );
                        break;
                    default:
                        b.append( n );
                }
            }
            else {
                b.append( c );
            }
        }
        return b.toString();
    }

    private static String emptyToNull( final String s ) {
        return ( ( s == null ) || s.isEmpty() ) ? null : s;
    }

    /** A parsed cache line: its lookup key, write timestamp, and the taxon record. */
    static final class Parsed {

        final String      _key;
        final long        _epoch_ms;
        final CachedTaxon _taxon;

        Parsed( final String key, final long epoch_ms, final CachedTaxon taxon ) {
            _key = key;
            _epoch_ms = epoch_ms;
            _taxon = taxon;
        }
    }
}
