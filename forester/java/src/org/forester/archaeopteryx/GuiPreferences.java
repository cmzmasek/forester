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
import java.util.function.Predicate;

import org.forester.archaeopteryx.Options.OVERVIEW_PLACEMENT_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.Options.SUPPORT_VISUALIZATION;
import org.forester.archaeopteryx.Options.TIP_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.TREE_ORIENTATION;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
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
            boolPref( "show_scale_axis", Options::isShowScaleAxis, Options::setShowScaleAxis ),
            boolPref( "show_hpd_bars", Options::isShowHpdBars, Options::setShowHpdBars ),
            boolPref( "show_zebra_stripes", Options::isShowZebraStripes, Options::setShowZebraStripes ),
            // key kept as "flip_vertically" (its original name) so existing saved settings still load after the
            // feature was renamed to "Reverse Tip Order"; do NOT rename the key.
            boolPref( "flip_vertically", Options::isReverseTipOrder, Options::setReverseTipOrder ),
            boolPref( "bold_found_labels", Options::isBoldFoundLabels, Options::setBoldFoundLabels ),
            boolPref( "dim_non_matches", Options::isDimNonMatches, Options::setDimNonMatches ),
            boolPref( "pulse_found_nodes", Options::isPulseFoundNodes, Options::setPulseFoundNodes ),
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
            // Numeric prefs pass the SAME [min,max] their Settings spinner enforces: a stored value is clamped into
            // range on load, so a corrupt/older/hand-edited file can never seed a SpinnerNumberModel out of range
            // (which THROWS) or feed a NaN/huge value into paint/export (NaN/±Infinity order outside any finite range,
            // so they clamp to a bound too).
            shortPref( "default_node_shape_size", Options::getDefaultNodeShapeSize,
                       Options::setDefaultNodeShapeSize, (short) 0, (short) 100 ),
            enumPref( "default_node_fill", Options::getDefaultNodeFill, Options::setDefaultNodeFill,
                      NodeFill::valueOf ),
            floatPref( "default_branch_width", Options::getDefaultBranchWidth, Options::setDefaultBranchWidth,
                       0.5f, 20f ),
            enumPref( "support_visualization", Options::getSupportVisualization, Options::setSupportVisualization,
                      SUPPORT_VISUALIZATION::valueOf ),
            doublePref( "support_threshold", Options::getSupportThreshold, Options::setSupportThreshold, 0.0, 1.0 ),
            doublePref( "min_confidence_fraction", Options::getMinConfidenceFraction,
                        Options::setMinConfidenceFraction, 0.0, 1.0 ),
            // Display-tab layout: tree style + overview placement (both clean Options enums). Tree style deliberately
            // does NOT restore the (alpha) CIRCULAR/UNROOTED radial modes: on load, the graphics-type-blind
            // lookAtSomeTreePropertiesForAptxControlSettings would re-enable the phylogram axis for a branch-length
            // tree, drawing a forbidden "circular phylogram". Restoring only the rectangular-family types avoids that;
            // a session left in a radial mode simply reopens rectangular (radial parity is deferred anyway).
            enumPref( "phylogeny_graphics_type", Options::getPhylogenyGraphicsType,
                      Options::setPhylogenyGraphicsType, PHYLOGENY_GRAPHICS_TYPE::valueOf,
                      t -> ( t != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) && ( t != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) ),
            enumPref( "overview_placement", Options::getOvPlacement, Options::setOvPlacement,
                      OVERVIEW_PLACEMENT_TYPE::valueOf ),
            // Root orientation (rectangular family): ROOT_LEFT / ROOT_TOP / ROOT_BOTTOM. All three are safe to
            // restore -- unlike the radial graphics types above, orientation is simply a no-op under CIRCULAR/
            // UNROOTED (the paint path gates it), so no exclusion predicate is needed.
            enumPref( "tree_orientation", Options::getTreeOrientation, Options::setTreeOrientation,
                      TREE_ORIENTATION::valueOf ),
            // the vertical-orientation tip-label angle (45deg vs vertical); a no-op under any non-vertical layout
            enumPref( "tip_label_direction", Options::getTipLabelDirection, Options::setTipLabelDirection,
                      TIP_LABEL_DIRECTION::valueOf ),
            // "Color by" palette: only apply a stored name that is still a known palette (a renamed/removed
            // palette is ignored, keeping the default) so a stale file can't select a non-existent palette
            stringPref( "color_palette", Options::getColorPaletteName, Options::setColorPaletteName,
                        v -> PropertyColorScheme.paletteNames().contains( v ) ),
            // Export appearance: raster scale + the two background toggles (white background already persisted above)
            intPref( "raster_export_scale", Options::getRasterExportScale, Options::setRasterExportScale, 1, 8 ),
            boolPref( "transparent_export_background", Options::isTransparentExportBackground,
                      Options::setTransparentExportBackground ),
            boolPref( "graphics_export_visible_only", Options::isGraphicsExportVisibleOnly,
                      Options::setGraphicsExportVisibleOnly ) );

    /** A boolean setting, stored as {@code "true"}/{@code "false"}. Trimmed so a hand-edited {@code "true "} still
     *  parses (Properties keeps trailing whitespace). */
    private static Pref boolPref( final String key, final Function<Options, Boolean> getter,
                                  final BiConsumer<Options, Boolean> setter ) {
        return new Pref( key,
                         o -> Boolean.toString( getter.apply( o ) ),
                         ( o, v ) -> setter.accept( o, Boolean.parseBoolean( v.trim() ) ) );
    }

    /** An enum setting, stored as its {@link Enum#name()}, applying every valid constant. */
    private static <E extends Enum<E>> Pref enumPref( final String key, final Function<Options, E> getter,
                                                      final BiConsumer<Options, E> setter,
                                                      final Function<String, E> parser ) {
        return enumPref( key, getter, setter, parser, e -> true );
    }

    /** An enum setting, stored as its {@link Enum#name()}. A value not in the enum (renamed/removed constant, or a
     *  corrupt file) OR one {@code accept} rejects is ignored, so the option keeps its default. */
    private static <E extends Enum<E>> Pref enumPref( final String key, final Function<Options, E> getter,
                                                      final BiConsumer<Options, E> setter,
                                                      final Function<String, E> parser, final Predicate<E> accept ) {
        return new Pref( key,
                         o -> getter.apply( o ).name(),
                         ( o, v ) -> {
                             final E parsed = parseOrNull( parser, v );
                             if ( ( parsed != null ) && accept.test( parsed ) ) {
                                 setter.accept( o, parsed );
                             }
                         } );
    }

    /** A numeric setting parsed with {@code parser} then CLAMPED into [{@code min},{@code max}]. A non-parseable
     *  value is ignored (default kept); a finite out-of-range value -- and NaN/±Infinity, which order outside any
     *  finite range -- is clamped to the nearest bound. The single writer shared by all the numeric factories. */
    private static <T extends Comparable<T>> Pref numericPref( final String key,
            final Function<Options, String> reader, final BiConsumer<Options, T> setter,
            final Function<String, T> parser, final T min, final T max ) {
        return new Pref( key, reader, ( o, v ) -> {
            final T parsed = parseOrNull( parser, v );
            if ( parsed != null ) {
                setter.accept( o, ( parsed.compareTo( min ) < 0 ) ? min : ( parsed.compareTo( max ) > 0 ) ? max : parsed );
            }
        } );
    }

    /** A short-integer setting (e.g. node size), clamped into [{@code min},{@code max}]. */
    private static Pref shortPref( final String key, final Function<Options, Short> getter,
                                   final BiConsumer<Options, Short> setter, final short min, final short max ) {
        return numericPref( key, o -> Short.toString( getter.apply( o ) ), setter, Short::valueOf, min, max );
    }

    /** A float setting (e.g. branch width), clamped into [{@code min},{@code max}]. */
    private static Pref floatPref( final String key, final Function<Options, Float> getter,
                                   final BiConsumer<Options, Float> setter, final float min, final float max ) {
        return numericPref( key, o -> Float.toString( getter.apply( o ) ), setter, Float::valueOf, min, max );
    }

    /** An int setting (e.g. raster export scale), clamped into [{@code min},{@code max}]. */
    private static Pref intPref( final String key, final Function<Options, Integer> getter,
                                 final BiConsumer<Options, Integer> setter, final int min, final int max ) {
        return numericPref( key, o -> Integer.toString( getter.apply( o ) ), setter, Integer::valueOf, min, max );
    }

    /** A double setting (e.g. support threshold), clamped into [{@code min},{@code max}]. */
    private static Pref doublePref( final String key, final Function<Options, Double> getter,
                                    final BiConsumer<Options, Double> setter, final double min, final double max ) {
        return numericPref( key, o -> Double.toString( getter.apply( o ) ), setter, Double::valueOf, min, max );
    }

    /** A string setting, applied (trimmed) only when {@code valid} accepts it (so a stale/unknown value -- e.g. a
     *  palette that no longer exists -- is ignored and the default stands). */
    private static Pref stringPref( final String key, final Function<Options, String> getter,
                                    final BiConsumer<Options, String> setter, final Predicate<String> valid ) {
        return new Pref( key,
                         getter,
                         ( o, v ) -> {
                             final String trimmed = v.trim();
                             if ( valid.test( trimmed ) ) {
                                 setter.accept( o, trimmed );
                             }
                         } );
    }

    private static <E> E parseOrNull( final Function<String, E> parser, final String v ) {
        try {
            return parser.apply( v.trim() );
        }
        catch ( final RuntimeException ex ) {
            return null; // e.g. IllegalArgumentException for an unknown enum constant, NumberFormatException for a number
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
            final String value = pref.reader().apply( options );
            if ( value != null ) { // Properties forbids null values; skip rather than abort the whole save (never throw)
                p.setProperty( pref.key(), value );
            }
        }
        save( p );
    }

    /** Best-effort delete of the persisted settings file so a "Reset to Defaults" survives the next launch (the
     *  next startup then finds no file and uses the built-in defaults). Never throws; a missing file is fine. */
    void deleteSettingsFile() {
        try {
            Files.deleteIfExists( _file );
        }
        catch ( final Exception e ) {
            // best-effort: if it can't be removed, the reset still applied to the live session
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
