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
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

/**
 * The files most recently opened, for File &rarr; Open Recent.
 * <p>
 * Stored beside the other cross-session state in {@code ${user.home}/.archaeopteryx/recent-files.properties} --
 * not in {@link GuiPreferences}' display settings, because a recent-file list is history rather than a preference
 * and must survive "Reset to Defaults" (resetting how the tree LOOKS should not forget what you were working on).
 * Honours the same {@code archaeopteryx.cache.dir} override, so all of Archaeopteryx's state relocates together.
 * <p>
 * Best-effort throughout, like {@link DirectoryPreferences}: every failure degrades to "no recent files" rather
 * than interrupting a load or a launch.
 */
final class RecentFiles {

    static final String RECENT_FILES_FILE = "recent-files.properties";
    /** How many to remember. Ten is the usual length -- long enough to span a working session, short enough that
     *  the submenu stays scannable. */
    static final int    MAX_ENTRIES       = 10;
    private static final String KEY_PREFIX = "recent.";

    private final Path _file;

    RecentFiles() {
        this( resolveDir().resolve( RECENT_FILES_FILE ) );
    }

    /** For tests: an explicit file. The constructor performs no I/O. */
    RecentFiles( final Path file ) {
        _file = file;
    }

    private static Path resolveDir() {
        final String override = System.getProperty( GuiPreferences.DIR_PROPERTY );
        if ( ( override != null ) && !override.trim().isEmpty() ) {
            return Paths.get( override.trim() );
        }
        return Paths.get( System.getProperty( "user.home", "." ), GuiPreferences.DEFAULT_DIR );
    }

    /**
     * The remembered paths, most recent first.
     * <p>
     * Returned WITHOUT touching the filesystem: statting each path here would block the launch if one of them
     * lives on a network mount that has gone away (the same reason {@link DirectoryPreferences} defers
     * validation). Whether a file still exists is the menu builder's question, not this one's.
     */
    List<String> load() {
        final List<String> out = new ArrayList<String>();
        final Properties p = new Properties();
        try {
            if ( !Files.exists( _file ) ) {
                return out;
            }
            try ( final java.io.InputStream in = Files.newInputStream( _file ) ) {
                p.load( in );
            }
        }
        catch ( final Exception e ) {
            return out; // best-effort: no history rather than a failed launch
        }
        for( int i = 0; i < MAX_ENTRIES; i++ ) {
            final String v = p.getProperty( KEY_PREFIX + i );
            if ( ( v != null ) && !v.trim().isEmpty() ) {
                out.add( v.trim() );
            }
        }
        return out;
    }

    /** Records {@code file} as the most recent, moving it up if it was already known and dropping the oldest past
     *  {@link #MAX_ENTRIES}. A null/unusable file is ignored. */
    void add( final File file ) {
        if ( file == null ) {
            return;
        }
        String path;
        try {
            path = file.getCanonicalPath();
        }
        catch ( final Exception e ) {
            path = file.getAbsolutePath(); // a canonical path is nicer, but any stable path will do
        }
        if ( ( path == null ) || path.trim().isEmpty() ) {
            return;
        }
        save( merge( load(), path ) );
    }

    /** Pure: {@code path} first, previous entries after it with any duplicate of it removed, capped. Comparison is
     *  exact -- two spellings of one file would both be listed, which is harmless and better than guessing that
     *  two paths are the same file across symlinks and network mounts. */
    static List<String> merge( final List<String> existing, final String path ) {
        final List<String> out = new ArrayList<String>();
        out.add( path );
        if ( existing != null ) {
            for( final String e : existing ) {
                if ( ( e != null ) && !e.equals( path ) && ( out.size() < MAX_ENTRIES ) ) {
                    out.add( e );
                }
            }
        }
        return out;
    }

    /** Forgets everything (File &rarr; Open Recent &rarr; Clear Menu). */
    void clear() {
        save( new ArrayList<String>() );
    }

    private void save( final List<String> paths ) {
        final Properties p = new Properties();
        for( int i = 0; i < paths.size(); i++ ) {
            p.setProperty( KEY_PREFIX + i, paths.get( i ) );
        }
        try {
            final Path dir = _file.getParent();
            if ( dir != null ) {
                Files.createDirectories( dir );
            }
            try ( final OutputStream out = Files.newOutputStream( _file ) ) {
                p.store( out, "Archaeopteryx recently opened files (most recent first)" );
            }
        }
        catch ( final Exception e ) {
            // best-effort: the list just won't persist this time
        }
    }
}
