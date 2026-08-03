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

import org.forester.phylogeny.data.NodeVisualData.NodeShape;

/**
 * Persists a curated set of <em>display</em> settings across restarts, so the view a user prefers (tree name
 * on/off, scale, italic species names, default node shape and size, ...) is the view they get next time. Each
 * setting is stored as a string; booleans, enums (e.g. node shape) and numbers (e.g. node size) all round-trip
 * through the same file at {@code ${user.home}/.archaeopteryx/display-settings.properties}.
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

    /** One persisted setting: a stable string key plus a codec that reads it from / writes it to {@link Options}
     *  as a string. {@code reader} turns the option's current value into the stored string; {@code writer} parses
     *  a stored string back onto the option (and silently ignores an unparseable value, so a corrupt file never
     *  breaks startup). Add to {@link #PREFS} via the typed factories below -- never renumber existing keys. */
    private record Pref(String key, Function<Options, String> reader, BiConsumer<Options, String> writer) {
    }

    private static final List<Pref> PREFS = List.of(
            boolPref( "show_tree_name", Options::isShowTreeName, Options::setShowTreeName ),
            boolPref( "show_scale", Options::isShowScale, Options::setShowScale ),
            boolPref( "show_scale_grid", Options::isShowScaleGrid, Options::setShowScaleGrid ),
            boolPref( "show_overview", Options::isShowOverview, Options::setShowOverview ),
            boolPref( "abbreviate_scientific_names", Options::isAbbreviateScientificTaxonNames,
                      Options::setAbbreviateScientificTaxonNames ),
            boolPref( "use_italic_scientific_names", Options::isUseItalicScientificNames,
                      Options::setUseItalicScientificNames ),
            boolPref( "antialias_print", Options::isAntialiasPrint, Options::setAntialiasPrint ),
            boolPref( "graphics_export_white_background", Options::isGraphicsExportWhiteBackground,
                      Options::setGraphicsExportWhiteBackground ),
            enumPref( "default_node_shape", Options::getDefaultNodeShape, Options::setDefaultNodeShape,
                      NodeShape::valueOf ),
            shortPref( "default_node_shape_size", Options::getDefaultNodeShapeSize,
                       Options::setDefaultNodeShapeSize ) );

    /** A boolean setting, stored as {@code "true"}/{@code "false"}. */
    private static Pref boolPref( final String key, final Function<Options, Boolean> getter,
                                  final BiConsumer<Options, Boolean> setter ) {
        return new Pref( key,
                         o -> Boolean.toString( getter.apply( o ) ),
                         ( o, v ) -> setter.accept( o, Boolean.parseBoolean( v ) ) );
    }

    /** An enum setting, stored as its {@link Enum#name()}. A value not in the enum (renamed/removed constant, or a
     *  corrupt file) is ignored, so the option keeps its default. */
    private static <E extends Enum<E>> Pref enumPref( final String key, final Function<Options, E> getter,
                                                      final BiConsumer<Options, E> setter,
                                                      final Function<String, E> parser ) {
        return new Pref( key,
                         o -> getter.apply( o ).name(),
                         ( o, v ) -> {
                             final E parsed = parseOrNull( parser, v );
                             if ( parsed != null ) {
                                 setter.accept( o, parsed );
                             }
                         } );
    }

    /** A short-integer setting (e.g. node size). A non-numeric stored value is ignored, keeping the default. */
    private static Pref shortPref( final String key, final Function<Options, Short> getter,
                                   final BiConsumer<Options, Short> setter ) {
        return new Pref( key,
                         o -> Short.toString( getter.apply( o ) ),
                         ( o, v ) -> {
                             final Short parsed = parseShortOrNull( v );
                             if ( parsed != null ) {
                                 setter.accept( o, parsed );
                             }
                         } );
    }

    private static <E> E parseOrNull( final Function<String, E> parser, final String v ) {
        try {
            return parser.apply( v );
        }
        catch ( final RuntimeException ex ) {
            return null; // e.g. IllegalArgumentException for an unknown enum constant
        }
    }

    private static Short parseShortOrNull( final String v ) {
        try {
            return Short.valueOf( v.trim() );
        }
        catch ( final RuntimeException ex ) {
            return null;
        }
    }

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

    /** Overlays the persisted settings onto {@code options}. Options whose key is not in the file (or whose stored
     *  value does not parse) are left as they are (their config/built-in default stands). Never throws. */
    void applyTo( final Options options ) {
        if ( options == null ) {
            return;
        }
        final Properties p = load();
        for( final Pref pref : PREFS ) {
            final String v = p.getProperty( pref.key() );
            if ( v != null ) {
                pref.writer().accept( options, v );
            }
        }
    }

    /** Writes the current value of every persisted setting from {@code options} to disk. Never throws. */
    void saveFrom( final Options options ) {
        if ( options == null ) {
            return;
        }
        final Properties p = new Properties();
        for( final Pref pref : PREFS ) {
            p.setProperty( pref.key(), pref.reader().apply( options ) );
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
