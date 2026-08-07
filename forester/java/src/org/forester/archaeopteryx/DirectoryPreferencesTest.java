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
import java.util.EnumMap;
import java.util.Map;

import org.forester.archaeopteryx.DirectoryPreferences.Category;

/**
 * Headless unit test for {@link DirectoryPreferences}: the three last-used directories (OPEN / SAVE / EXPORT)
 * round-trip through the properties file; a directory that no longer exists is dropped on load (so the dialog
 * falls back to home rather than a dead path); only the categories that were set are restored; and a corrupt /
 * empty store never throws.
 */
public final class DirectoryPreferencesTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DirectoryPreferences: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return roundTrip() && merges() && skipsEmptyAndMissing();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected: " + e );
        }
    }

    /** All three categories save and reload to the same absolute paths. */
    private static boolean roundTrip() throws Exception {
        final Path store = Files.createTempFile( "aptx-dirs", ".properties" );
        final File open = Files.createTempDirectory( "aptx-open" ).toFile();
        final File save = Files.createTempDirectory( "aptx-save" ).toFile();
        final File export = Files.createTempDirectory( "aptx-export" ).toFile();
        final Map<Category, File> out = new EnumMap<>( Category.class );
        out.put( Category.OPEN, open );
        out.put( Category.SAVE, save );
        out.put( Category.EXPORT, export );
        new DirectoryPreferences( store ).saveFrom( out );

        final Map<Category, File> in = new EnumMap<>( Category.class );
        new DirectoryPreferences( store ).applyTo( in );
        if ( !samePath( in.get( Category.OPEN ), open ) ) {
            return fail( "OPEN did not round-trip: " + in.get( Category.OPEN ) );
        }
        if ( !samePath( in.get( Category.SAVE ), save ) ) {
            return fail( "SAVE did not round-trip: " + in.get( Category.SAVE ) );
        }
        if ( !samePath( in.get( Category.EXPORT ), export ) ) {
            return fail( "EXPORT did not round-trip: " + in.get( Category.EXPORT ) );
        }
        return true;
    }

    /** saveFrom MERGES onto the existing file: a category written in an earlier save (or by another window) is kept
     *  when a later save touches a different category; touching the same category overrides it. */
    private static boolean merges() throws Exception {
        final Path store = Files.createTempFile( "aptx-dirs", ".properties" );
        final File open = Files.createTempDirectory( "aptx-open2" ).toFile();
        final File save = Files.createTempDirectory( "aptx-save2" ).toFile();
        final File save2 = Files.createTempDirectory( "aptx-save3" ).toFile();

        final Map<Category, File> first = new EnumMap<>( Category.class );
        first.put( Category.OPEN, open );
        new DirectoryPreferences( store ).saveFrom( first ); // OPEN=open

        final Map<Category, File> second = new EnumMap<>( Category.class );
        second.put( Category.SAVE, save );
        new DirectoryPreferences( store ).saveFrom( second ); // adds SAVE=save, must KEEP OPEN=open

        final Map<Category, File> in = new EnumMap<>( Category.class );
        new DirectoryPreferences( store ).applyTo( in );
        if ( !samePath( in.get( Category.OPEN ), open ) || !samePath( in.get( Category.SAVE ), save ) ) {
            return fail( "merge should keep OPEN and add SAVE, got " + in );
        }

        final Map<Category, File> third = new EnumMap<>( Category.class );
        third.put( Category.SAVE, save2 );
        new DirectoryPreferences( store ).saveFrom( third ); // overrides SAVE
        final Map<Category, File> in2 = new EnumMap<>( Category.class );
        new DirectoryPreferences( store ).applyTo( in2 );
        if ( !samePath( in2.get( Category.SAVE ), save2 ) || !samePath( in2.get( Category.OPEN ), open ) ) {
            return fail( "re-saving SAVE should override it and keep OPEN, got " + in2 );
        }
        return true;
    }

    /** An empty session writes no file (no stray dotfile on first quit); paths load raw (validation is the caller's
     *  job, deferred off the startup path); a missing file loads to nothing without throwing. */
    private static boolean skipsEmptyAndMissing() throws Exception {
        // (a) saving an empty map creates no file
        final Path store = Files.createTempDirectory( "aptx-empty" ).resolve( "dirs.properties" );
        new DirectoryPreferences( store ).saveFrom( new EnumMap<>( Category.class ) );
        if ( Files.exists( store ) ) {
            return fail( "an empty save must not create a file, but " + store + " exists" );
        }

        // (b) a stored path is loaded RAW -- even one that no longer exists (getCurrentDir validates later, lazily)
        final Path store2 = Files.createTempFile( "aptx-dirs", ".properties" );
        final File gone = Files.createTempDirectory( "aptx-gone" ).toFile();
        final Map<Category, File> out = new EnumMap<>( Category.class );
        out.put( Category.SAVE, gone );
        new DirectoryPreferences( store2 ).saveFrom( out );
        gone.delete(); // remove after saving the path
        final Map<Category, File> in = new EnumMap<>( Category.class );
        new DirectoryPreferences( store2 ).applyTo( in );
        if ( !in.containsKey( Category.SAVE ) || in.containsKey( Category.OPEN ) ) {
            return fail( "applyTo should load the stored path raw (no filesystem stat) and only set categories, got " + in );
        }

        // (c) a missing file loads to nothing without throwing
        final Map<Category, File> in3 = new EnumMap<>( Category.class );
        new DirectoryPreferences( Files.createTempDirectory( "aptx-none" ).resolve( "nope.properties" ) )
                .applyTo( in3 );
        if ( !in3.isEmpty() ) {
            return fail( "a missing store should load to nothing, got " + in3 );
        }
        return true;
    }

    private static boolean samePath( final File a, final File b ) {
        return ( a != null ) && ( b != null ) && a.getAbsolutePath().equals( b.getAbsolutePath() );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [DirectoryPreferencesTest] " + msg );
        return false;
    }

    private DirectoryPreferencesTest() {
    }
}
