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
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Unit tests for {@link TaxonomyDiskCache} and {@link TaxonomyCacheStatus} -- entirely local (no
 * network): the line serialization (incl. escaping nasty characters), TTL expiry + last-wins
 * compaction, positive-only persistence, the enable toggle, graceful degradation when the cache
 * directory cannot be created, the display-formatting helpers, and that
 * {@link NcbiTaxonomyLineageService} seeds its in-memory caches from a pre-populated disk cache.
 *
 * <p>Each filesystem case runs against a temporary directory selected via the
 * {@value TaxonomyDiskCache#DIR_PROPERTY} system property (restored at the end).
 */
public final class TaxonomyDiskCacheTest {

    private static final long DAY_MS = 24L * 60L * 60L * 1000L;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TaxonomyDiskCache: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return testSerialization() && testFormatting() && testFilesystem();
    }

    private static boolean testSerialization() {
        // a record whose ancestor carries every structural delimiter -- must survive the round-trip on one line
        final List<String[]> anc = Arrays.asList( new String[] { "cellular organisms", "no rank" },
                                                  new String[] { "Eukaryota", "superkingdom" },
                                                  new String[] { "weird|name~with\ttabs", "" } );
        final CachedTaxon ct = new CachedTaxon( "9682", "genus", "Felis", "cats", anc );
        final String line = TaxonomyDiskCache.toLine( "felis", ct, 1234567890L );
        if ( line.indexOf( '\n' ) >= 0 ) {
            return fail( "a serialized entry must be a single line" );
        }
        final TaxonomyDiskCache.Parsed p = TaxonomyDiskCache.parseLine( line );
        if ( ( p == null ) || !"felis".equals( p._key ) || ( p._epoch_ms != 1234567890L ) ) {
            return fail( "key/epoch did not round-trip: " + p );
        }
        final CachedTaxon back = p._taxon;
        if ( !"9682".equals( back.getTaxId() ) || !"genus".equals( back.getRank() )
                || !"Felis".equals( back.getScientificName() ) || !"cats".equals( back.getCommonName() ) ) {
            return fail( "scalar fields did not round-trip" );
        }
        if ( ( back.getAncestors().size() != 3 )
                || !"cellular organisms".equals( back.getAncestors().get( 0 )[ 0 ] )
                || !"no rank".equals( back.getAncestors().get( 0 )[ 1 ] )
                || !"weird|name~with\ttabs".equals( back.getAncestors().get( 2 )[ 0 ] )
                || !"".equals( back.getAncestors().get( 2 )[ 1 ] ) ) {
            return fail( "ancestors (incl. special characters) did not round-trip: " + back.getAncestors().get( 2 )[ 0 ] );
        }
        // esc must remove every structural character; unesc must invert it exactly
        final String nasty = "a\\b\tc\nd|e~f";
        if ( !nasty.equals( TaxonomyDiskCache.unesc( TaxonomyDiskCache.esc( nasty ) ) ) ) {
            return fail( "esc/unesc is not invertible" );
        }
        final String escaped = TaxonomyDiskCache.esc( nasty );
        if ( ( escaped.indexOf( '\t' ) >= 0 ) || ( escaped.indexOf( '\n' ) >= 0 ) || ( escaped.indexOf( '|' ) >= 0 )
                || ( escaped.indexOf( '~' ) >= 0 ) ) {
            return fail( "esc must not leave any literal structural character" );
        }
        // malformed lines must parse to null, never throw
        if ( ( TaxonomyDiskCache.parseLine( "only\ttwo" ) != null )
                || ( TaxonomyDiskCache.parseLine( "k\tNOTANUMBER\t\t\t\t\t" ) != null )
                || ( TaxonomyDiskCache.parseLine( null ) != null ) ) {
            return fail( "malformed/short/null lines must parse to null" );
        }
        return true;
    }

    private static boolean testFormatting() {
        if ( !"0 B".equals( TaxonomyCacheStatus.formatBytes( 0 ) ) || !"512 B".equals( TaxonomyCacheStatus.formatBytes( 512 ) )
                || !"2 KB".equals( TaxonomyCacheStatus.formatBytes( 2048 ) )
                || !TaxonomyCacheStatus.formatBytes( 5L * 1024L * 1024L ).endsWith( "MB" ) ) {
            return fail( "formatBytes wrong" );
        }
        final long now = 1000L * DAY_MS;
        if ( !"".equals( TaxonomyCacheStatus.describeAge( 0L, now ) )
                || !"Oldest entry: today".equals( TaxonomyCacheStatus.describeAge( now, now ) )
                || !"Oldest entry: 1 day ago".equals( TaxonomyCacheStatus.describeAge( now - DAY_MS, now ) )
                || !"Oldest entry: 3 days ago".equals( TaxonomyCacheStatus.describeAge( now - ( 3L * DAY_MS ), now ) ) ) {
            return fail( "describeAge wrong" );
        }
        return true;
    }

    private static boolean testFilesystem() {
        final String saved = System.getProperty( TaxonomyDiskCache.DIR_PROPERTY );
        Path dir = null;
        Path blocker = null;
        try {
            dir = Files.createTempDirectory( "aptx-taxcache-test" );
            System.setProperty( TaxonomyDiskCache.DIR_PROPERTY, dir.toString() );
            final Path file = dir.resolve( TaxonomyDiskCache.CACHE_FILE );
            final CachedTaxon ct = felis( "Felis" );
            final CachedTaxon ctV2 = felis( "FelisV2" );

            // TTL + last-wins compaction: a stale "stale" key is dropped, a duplicated "fresh" key keeps the last write
            final long now = System.currentTimeMillis();
            final List<String> seeded = new ArrayList<String>();
            seeded.add( TaxonomyDiskCache.toLine( "fresh", ct, now ) );
            seeded.add( TaxonomyDiskCache.toLine( "stale", ct, now - TaxonomyDiskCache.TTL_MS - 100000L ) );
            seeded.add( TaxonomyDiskCache.toLine( "fresh", ctV2, now ) );
            Files.write( file, seeded, StandardCharsets.UTF_8 );

            final TaxonomyDiskCache cache = new TaxonomyDiskCache();
            final Map<String, CachedTaxon> loaded = cache.load();
            if ( loaded.containsKey( "stale" ) ) {
                return fail( "an entry older than the TTL must be dropped on load" );
            }
            if ( !loaded.containsKey( "fresh" ) || !"FelisV2".equals( loaded.get( "fresh" ).getScientificName() ) ) {
                return fail( "the last write for a key must win" );
            }
            if ( Files.readAllLines( file, StandardCharsets.UTF_8 ).size() != 1 ) {
                return fail( "load must compact to one line per surviving key" );
            }

            // status reflects the compacted, non-expired set
            final TaxonomyCacheStatus st = cache.status();
            if ( !st.isAvailable() || ( st.getEntries() != 1 ) || ( st.getOldestEpochMs() <= 0L )
                    || !st.getPath().contains( TaxonomyDiskCache.CACHE_FILE ) ) {
                return fail( "status after compaction wrong: entries=" + st.getEntries() );
            }

            // positive-only: an empty record must not be appended
            final long size_before = Files.size( file );
            cache.put( "empty", new CachedTaxon( null, null, null, null, null ) );
            if ( Files.size( file ) != size_before ) {
                return fail( "an empty/not-found record must never be persisted" );
            }
            // a real put appends and is readable back
            cache.put( "mus", felis( "Mus" ) );
            if ( !cache.load().containsKey( "mus" ) ) {
                return fail( "a positive record must be persisted and reloadable" );
            }

            // enable toggle persists via the marker file, across instances
            cache.setEnabled( false );
            if ( cache.isEnabled() || new TaxonomyDiskCache().isEnabled() ) {
                return fail( "disabling must persist (marker file) and be visible to a fresh instance" );
            }
            if ( !cache.load().isEmpty() ) {
                return fail( "a disabled cache must load nothing" );
            }
            cache.setEnabled( true );
            if ( !cache.isEnabled() || cache.load().isEmpty() ) {
                return fail( "re-enabling must restore the cache" );
            }

            // the service seeds BOTH in-memory caches from disk (no network for a seeded taxon)
            Files.write( file, Collections.singletonList( TaxonomyDiskCache.toLine( "felis", ct, now ) ),
                         StandardCharsets.UTF_8 );
            final NcbiTaxonomyLineageService svc = new NcbiTaxonomyLineageService();
            final RankedLineage rl = svc.lineageOf( "Felis" ); // cache-only path -> must be a warm disk hit
            if ( ( rl == null ) || rl.isEmpty() || !"Felis".equals( rl.at( "genus" ) )
                    || !"Carnivora".equals( rl.at( "order" ) ) ) {
                return fail( "service did not seed the rank cache from disk" );
            }
            ResolvedTaxonomy rt;
            try {
                rt = svc.resolveTaxonomy( "felis" ); // seeded -> returns before any network call
            }
            catch ( final Exception e ) {
                return fail( "resolveTaxonomy on a seeded key must not reach the network: " + e );
            }
            if ( rt.isEmpty() || !"Felis".equals( rt.getScientificName() ) || !"9682".equals( rt.getTaxId() ) ) {
                return fail( "service did not seed the detail cache from disk: " + rt );
            }

            // graceful degradation: a cache dir that cannot be created -> unavailable, no throws, no caching
            blocker = Files.createTempFile( "aptx-taxcache-block", ".txt" ); // a regular file...
            System.setProperty( TaxonomyDiskCache.DIR_PROPERTY, blocker.resolve( "sub" ).toString() ); // ...used as a parent
            final TaxonomyDiskCache broken = new TaxonomyDiskCache();
            final TaxonomyCacheStatus bs = broken.status();
            if ( bs.isAvailable() || ( bs.getUnavailableReason() == null ) ) {
                return fail( "an uncreatable cache directory must report unavailable with a reason" );
            }
            broken.put( "x", ct ); // must not throw
            if ( !broken.load().isEmpty() ) {
                return fail( "an unavailable cache must load nothing (and not throw)" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "unexpected exception: " + e );
        }
        finally {
            if ( saved == null ) {
                System.clearProperty( TaxonomyDiskCache.DIR_PROPERTY );
            }
            else {
                System.setProperty( TaxonomyDiskCache.DIR_PROPERTY, saved );
            }
            deleteQuietly( dir );
            deleteQuietly( blocker );
        }
    }

    /** A small Felis-like record with a Linnaean order so seeded lineages can be checked. */
    private static CachedTaxon felis( final String sci_name ) {
        final List<String[]> anc = Arrays.asList( new String[] { "Eukaryota", "superkingdom" },
                                                  new String[] { "Mammalia", "class" },
                                                  new String[] { "Carnivora", "order" } );
        return new CachedTaxon( "9682", "genus", sci_name, "cats", anc );
    }

    private static void deleteQuietly( final Path p ) {
        if ( p == null ) {
            return;
        }
        try {
            if ( Files.isDirectory( p ) ) {
                try ( DirectoryStream<Path> ds = Files.newDirectoryStream( p ) ) {
                    for( final Path c : ds ) {
                        Files.deleteIfExists( c );
                    }
                }
            }
            Files.deleteIfExists( p );
        }
        catch ( final Exception e ) {
            // best-effort cleanup
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  TaxonomyDiskCache test failed: " + msg );
        return false;
    }

    private TaxonomyDiskCacheTest() {
    }
}
