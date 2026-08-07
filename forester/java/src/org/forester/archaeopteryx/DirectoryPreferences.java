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
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Properties;

/**
 * Persists the last-used file-dialog directories -- separately for reading ({@link Category#OPEN}), saving
 * ({@link Category#SAVE}), and exporting ({@link Category#EXPORT}) -- so each dialog reopens where the user last
 * left it, across restarts.
 *
 * <p>This is app <em>state</em>, not a display setting, so it lives in its own
 * {@code ${user.home}/.archaeopteryx/directories.properties} rather than in {@link GuiPreferences}' display-settings
 * file. That keeps "Reset to Defaults" (which deletes the display-settings file and resets {@link Options}) from
 * forgetting where the user's files are. It shares the {@code .archaeopteryx} directory and the
 * {@link GuiPreferences#DIR_PROPERTY} relocation override, so the whole of Archaeopteryx's on-disk state still moves
 * (and isolates in tests) together.
 *
 * <p>Best-effort, never a dependency: every disk operation is wrapped and any failure is swallowed, so a
 * missing/corrupt/unwritable file simply means "start in the default directory". A stored directory that no longer
 * exists (deleted, or a path from another machine) is ignored on load, so the dialog falls back to the home directory.
 */
final class DirectoryPreferences {

    /** Which file dialog a directory belongs to. The name is the stored key, so never rename an existing one. */
    enum Category {
        OPEN, SAVE, EXPORT
    }

    static final String DIRECTORIES_FILE = "directories.properties";

    private final Path _file;

    DirectoryPreferences() {
        this( resolveDir().resolve( DIRECTORIES_FILE ) );
    }

    /** For tests: an explicit file. The constructor performs no I/O. */
    DirectoryPreferences( final Path file ) {
        _file = file;
    }

    /** Same resolution as {@link GuiPreferences}, reusing its override property and default directory, so all of
     *  Archaeopteryx's state relocates together. */
    private static Path resolveDir() {
        final String override = System.getProperty( GuiPreferences.DIR_PROPERTY );
        if ( ( override != null ) && !override.trim().isEmpty() ) {
            return Paths.get( override.trim() );
        }
        return Paths.get( System.getProperty( "user.home", "." ), GuiPreferences.DEFAULT_DIR );
    }

    /** Overlays each persisted directory onto {@code dirs} (keyed by category). Paths are loaded WITHOUT touching
     *  the filesystem -- validation (does it still exist / is it readable?) is deferred to the caller's
     *  {@code getCurrentDir}, which is off the launch path. Statting a stored path here would block startup if it
     *  is on a slow/unmounted network share. Never throws. */
    void applyTo( final Map<Category, File> dirs ) {
        if ( dirs == null ) {
            return;
        }
        final Properties p = load();
        for( final Category c : Category.values() ) {
            final String v = p.getProperty( c.name() );
            if ( ( v != null ) && !v.trim().isEmpty() ) {
                dirs.put( c, new File( v.trim() ) );
            }
        }
    }

    /** Persists each set directory in {@code dirs}, MERGED onto whatever is already on disk (so a category this
     *  session did not touch -- or a second window's directory -- is kept, not clobbered). No filesystem stat (which
     *  could block on a stale share at quit); a since-deleted path is written but harmlessly ignored on next load.
     *  An empty result is NOT written, so a session that never used a file dialog creates no file. Never throws. */
    void saveFrom( final Map<Category, File> dirs ) {
        if ( dirs == null ) {
            return;
        }
        final Properties p = load(); // merge, don't overwrite: keep other categories'/windows' entries
        for( final Category c : Category.values() ) {
            final File dir = dirs.get( c );
            if ( dir != null ) {
                p.setProperty( c.name(), dir.getAbsolutePath() );
            }
        }
        if ( !p.isEmpty() ) {
            save( p );
        }
    }

    /** Best-effort read; returns an empty Properties on any problem (missing file, unreadable, corrupt). */
    private Properties load() {
        final Properties p = new Properties();
        try {
            if ( Files.exists( _file ) ) {
                try ( final InputStream in = Files.newInputStream( _file ) ) {
                    p.load( in );
                }
            }
        }
        catch ( final Exception e ) {
            // best-effort: fall back to the default directory
        }
        return p;
    }

    /** Best-effort write; creates the directory if needed and swallows any failure. */
    private void save( final Properties p ) {
        try {
            final Path dir = _file.getParent();
            if ( dir != null ) {
                Files.createDirectories( dir );
            }
            try ( final OutputStream out = Files.newOutputStream( _file ) ) {
                p.store( out, "Archaeopteryx last-used file-dialog directories" );
            }
        }
        catch ( final Exception e ) {
            // best-effort: the directories just won't persist this time
        }
    }
}
