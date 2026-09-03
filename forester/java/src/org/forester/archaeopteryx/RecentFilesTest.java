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

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

/**
 * The File &rarr; Open Recent history: ordering, capping, and that it survives a restart.
 * <p>
 * Headless -- the store is deliberately free of Swing so its ordering rules can be pinned without a frame.
 */
public final class RecentFilesTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RecentFiles: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return merging() && persistence() && degradesQuietly();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [RecentFilesTest] " + msg );
        return false;
    }

    /** The ordering rules, pure: newest first, no duplicates, capped. */
    private static boolean merging() {
        List<String> l = RecentFiles.merge( null, "/a" );
        if ( !l.equals( Arrays.asList( "/a" ) ) ) {
            return fail( "first entry: " + l );
        }
        l = RecentFiles.merge( Arrays.asList( "/a" ), "/b" );
        if ( !l.equals( Arrays.asList( "/b", "/a" ) ) ) {
            return fail( "newest must come first: " + l );
        }
        // re-opening a known file MOVES it up rather than listing it twice
        l = RecentFiles.merge( Arrays.asList( "/b", "/a", "/c" ), "/a" );
        if ( !l.equals( Arrays.asList( "/a", "/b", "/c" ) ) ) {
            return fail( "a re-opened file must move to the top, not duplicate: " + l );
        }
        // ...and the list is capped, dropping the OLDEST
        final String[] many = new String[ RecentFiles.MAX_ENTRIES ];
        for( int i = 0; i < many.length; i++ ) {
            many[ i ] = "/f" + i;
        }
        l = RecentFiles.merge( Arrays.asList( many ), "/new" );
        if ( l.size() != RecentFiles.MAX_ENTRIES ) {
            return fail( "the list must cap at " + RecentFiles.MAX_ENTRIES + ", got " + l.size() );
        }
        if ( !"/new".equals( l.get( 0 ) ) ) {
            return fail( "the newest must survive the cap" );
        }
        if ( l.contains( "/f" + ( many.length - 1 ) ) ) {
            return fail( "the OLDEST must be the one dropped: " + l );
        }
        return true;
    }

    /** It survives a restart -- a new instance over the same file sees the history. */
    private static boolean persistence() {
        Path dir = null;
        try {
            dir = Files.createTempDirectory( "aptx-recent" );
            final Path store = dir.resolve( RecentFiles.RECENT_FILES_FILE );
            final File f1 = Files.createFile( dir.resolve( "one.xml" ) ).toFile();
            final File f2 = Files.createFile( dir.resolve( "two.xml" ) ).toFile();
            new RecentFiles( store ).add( f1 );
            new RecentFiles( store ).add( f2 );
            final List<String> loaded = new RecentFiles( store ).load(); // a "restart"
            if ( loaded.size() != 2 ) {
                return fail( "expected 2 remembered files, got " + loaded );
            }
            if ( !loaded.get( 0 ).endsWith( "two.xml" ) ) {
                return fail( "the most recently opened file must be first, got " + loaded.get( 0 ) );
            }
            // stored paths are absolute, so the entry still resolves from any working directory
            if ( !new File( loaded.get( 0 ) ).isAbsolute() ) {
                return fail( "paths must be stored absolute, got " + loaded.get( 0 ) );
            }
            new RecentFiles( store ).clear();
            if ( !new RecentFiles( store ).load().isEmpty() ) {
                return fail( "clear() must empty the history" );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        finally {
            deleteQuietly( dir );
        }
    }

    /** Every failure degrades to "no history", never to an exception on a load or a launch. */
    private static boolean degradesQuietly() {
        try {
            // no file yet
            final Path missing = Files.createTempDirectory( "aptx-recent-none" ).resolve( "nope.properties" );
            if ( !new RecentFiles( missing ).load().isEmpty() ) {
                return fail( "a missing store must read as empty" );
            }
            deleteQuietly( missing.getParent() );
            // an unwritable location: add() must not throw
            final RecentFiles broken = new RecentFiles( new File( "/" ).toPath().resolve( "aptx-nope" )
                    .resolve( "x.properties" ) );
            broken.add( new File( "/tmp" ) );
            if ( !broken.load().isEmpty() ) {
                return fail( "an unwritable store must stay empty rather than pretend" );
            }
            // a null file is ignored rather than recorded
            final Path dir = Files.createTempDirectory( "aptx-recent-null" );
            final Path store = dir.resolve( RecentFiles.RECENT_FILES_FILE );
            new RecentFiles( store ).add( null );
            if ( !new RecentFiles( store ).load().isEmpty() ) {
                return fail( "a null file must be ignored" );
            }
            deleteQuietly( dir );
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void deleteQuietly( final Path dir ) {
        if ( dir == null ) {
            return;
        }
        try ( final java.util.stream.Stream<Path> s = Files.walk( dir ) ) {
            s.sorted( java.util.Comparator.reverseOrder() ).forEach( p -> {
                try {
                    Files.deleteIfExists( p );
                }
                catch ( final Exception e ) {
                    // best effort
                }
            } );
        }
        catch ( final Exception e ) {
            // best effort
        }
    }

    private RecentFilesTest() {
    }
}
