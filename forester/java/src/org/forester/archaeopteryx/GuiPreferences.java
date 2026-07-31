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

import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Properties;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Persists a curated set of boolean <em>display</em> toggles across restarts, so the view a user prefers
 * (tree name on/off, scale, italic species names, ...) is the view they get next time. Stored as a small
 * properties file at {@code ${user.home}/.archaeopteryx/display-settings.properties}.
 *
 * <p>The {@code .archaeopteryx} directory and the {@value #DIR_PROPERTY} override are shared, on purpose, with
 * the taxonomy disk cache: one well-known Archaeopteryx state directory, and one override that relocates ALL of
 * that state together (which is exactly what the test harness wants -- see Test.main). A properties FILE (rather
 * than {@link java.util.prefs.Preferences}, which the theme setting in {@link Configuration} uses) is deliberate:
 * it keeps the growing set of display settings in one human-readable, test-isolatable place and is the backing
 * store the deferred "persist all GUI settings" goal is meant to extend.
 *
 * <p>Best-effort, never a dependency: every disk operation is wrapped and any failure is swallowed, so a
 * missing/corrupt/unwritable file simply means "use the built-in defaults". Only the keys in {@link #PREFS}
 * are read or written; a key absent from the file leaves that option at whatever value it already had (so the
 * config-file / built-in default still applies). Loaded once at startup (before the menus are built from the
 * options) and written once on exit.
 */
final class GuiPreferences {

    static final String DIR_PROPERTY  = "archaeopteryx.cache.dir";
    static final String DEFAULT_DIR   = ".archaeopteryx";
    static final String SETTINGS_FILE = "display-settings.properties";

    /** The persisted toggles: a stable string key plus the Options getter/setter it maps to. Add to this list
     *  (never renumber existing keys) to remember more display options. NOTE: graphics_export_white_background is
     *  the one export-appearance toggle here; its interacting siblings (transparent-PNG, raster export scale) are
     *  not persisted yet -- part of the deferred "persist all settings" goal. */
    private record Pref(String key, Function<Options, Boolean> getter, BiConsumer<Options, Boolean> setter) {
    }

    private static final List<Pref> PREFS = List.of(
            new Pref( "show_tree_name", Options::isShowTreeName, Options::setShowTreeName ),
            new Pref( "show_scale", Options::isShowScale, Options::setShowScale ),
            new Pref( "show_scale_grid", Options::isShowScaleGrid, Options::setShowScaleGrid ),
            new Pref( "show_overview", Options::isShowOverview, Options::setShowOverview ),
            new Pref( "abbreviate_scientific_names", Options::isAbbreviateScientificTaxonNames,
                      Options::setAbbreviateScientificTaxonNames ),
            new Pref( "use_italic_scientific_names", Options::isUseItalicScientificNames,
                      Options::setUseItalicScientificNames ),
            new Pref( "antialias_print", Options::isAntialiasPrint, Options::setAntialiasPrint ),
            new Pref( "graphics_export_white_background", Options::isGraphicsExportWhiteBackground,
                      Options::setGraphicsExportWhiteBackground ) );

    private final Path _file;

    GuiPreferences() {
        this( resolveDir().resolve( SETTINGS_FILE ) );
    }

    /** For tests: an explicit settings file. The constructor performs no I/O. */
    GuiPreferences( final Path file ) {
        _file = file;
    }

    private static Path resolveDir() {
        final String override = System.getProperty( DIR_PROPERTY );
        if ( ( override != null ) && !override.trim().isEmpty() ) {
            return Paths.get( override.trim() );
        }
        return Paths.get( System.getProperty( "user.home", "." ), DEFAULT_DIR );
    }

    /** Overlays the persisted toggles onto {@code options}. Options whose key is not in the file are left as
     *  they are (their config/built-in default stands). Never throws. */
    void applyTo( final Options options ) {
        if ( options == null ) {
            return;
        }
        final Properties p = load();
        for( final Pref pref : PREFS ) {
            final String v = p.getProperty( pref.key() );
            if ( v != null ) {
                pref.setter().accept( options, Boolean.parseBoolean( v ) );
            }
        }
    }

    /** Writes the current value of every persisted toggle from {@code options} to disk. Never throws. */
    void saveFrom( final Options options ) {
        if ( options == null ) {
            return;
        }
        final Properties p = new Properties();
        for( final Pref pref : PREFS ) {
            p.setProperty( pref.key(), Boolean.toString( pref.getter().apply( options ) ) );
        }
        save( p );
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
            // best-effort: fall back to built-in defaults
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
                p.store( out, "Archaeopteryx display settings" );
            }
        }
        catch ( final Exception e ) {
            // best-effort: the settings just won't persist this time
        }
    }
}
