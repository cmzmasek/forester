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

import java.awt.Font;
import java.awt.FontFormatException;
import java.awt.GraphicsEnvironment;
import java.io.IOException;
import java.io.InputStream;

import org.forester.util.ForesterUtil;

/**
 * Registers the bundled, openly-licensed figure fonts (<b>Source Sans 3</b> -- the default -- plus
 * <b>Liberation Sans</b> and <b>Noto Sans</b>; all SIL OFL, see {@link #PREFERRED}) with the JVM so
 * Archaeopteryx renders identical, reproducible type on every platform regardless of which fonts the
 * user has installed -- instead of falling through proprietary, frequently-absent choices (Arial
 * Unicode MS, Arial, Helvetica) and landing on whatever the platform default happens to be.
 *
 * <p>All four faces (regular / bold / italic / bold-italic) are registered so {@code Font.BOLD} and
 * species-name italics resolve to real faces, not algorithmic slants. The font files ship in the jar
 * under {@value #FONT_DIR} (staged by the build's {@code copy_resources}).
 *
 * <p>Best-effort: a missing or unreadable face is skipped and the existing
 * {@link AptxConstants#DEFAULT_FONT_CHOICES} fallback chain still applies. Idempotent. Registration
 * invalidates {@link AptxUtil}'s cached available-font list, so the new family is visible to the
 * default-font selection in {@link Configuration} no matter the class-load order.
 */
public final class FontResources {

    /** The bundled default figure/UI font family (first of {@link #PREFERRED}). */
    public static final String  DEFAULT_FIGURE_FONT = "Source Sans 3";
    static final String         FONT_DIR            = "/resources/fonts/";

    /** A bundled, always-available figure font: family name, a one-line blurb, and its four face files. */
    public static final class Preferred {

        public final String    family;
        public final String    description;
        final String[]         faces;

        Preferred( final String family, final String description, final String... faces ) {
            this.family = family;
            this.description = description;
            this.faces = faces;
        }
    }

    /** The bundled fonts, in preference order (the first is the default). All are OFL-licensed. */
    public static final Preferred[] PREFERRED = {
            new Preferred( "Source Sans 3", "clean and highly legible at small sizes (recommended default)",
                           "SourceSans3-Regular.ttf", "SourceSans3-Bold.ttf", "SourceSans3-It.ttf",
                           "SourceSans3-BoldIt.ttf" ),
            new Preferred( "Liberation Sans", "metric-compatible with Arial / Helvetica (the journal-familiar look)",
                           "LiberationSans-Regular.ttf", "LiberationSans-Bold.ttf", "LiberationSans-Italic.ttf",
                           "LiberationSans-BoldItalic.ttf" ),
            new Preferred( "Noto Sans", "the widest character coverage (accents, diacritics, symbols)",
                           "NotoSans-Regular.ttf", "NotoSans-Bold.ttf", "NotoSans-Italic.ttf",
                           "NotoSans-BoldItalic.ttf" ), };

    private static boolean      _registered;

    /**
     * Registers every bundled face with the local {@link GraphicsEnvironment}. Idempotent and
     * best-effort; safe in headless mode. Call once at startup, before the default font is chosen.
     */
    public static synchronized void registerBundledFonts() {
        if ( _registered ) {
            return;
        }
        final GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
        boolean any = false;
        for( final Preferred p : PREFERRED ) {
            for( final String face : p.faces ) {
                try ( InputStream in = FontResources.class.getResourceAsStream( FONT_DIR + face ) ) {
                    if ( in == null ) {
                        continue; // not bundled in this build -- fall back to the chain
                    }
                    ge.registerFont( Font.createFont( Font.TRUETYPE_FONT, in ) );
                    any = true;
                }
                catch ( final IOException | FontFormatException e ) {
                    // best-effort: skip a bad face, keep the others and the fallback chain
                }
            }
        }
        if ( any ) {
            _registered = true; // only mark done once at least one face actually registered (so a misbuild can retry)
            // the just-registered families must be visible to the default-font selection
            AptxUtil.refreshAvailableFontFamilies();
        }
        else {
            // a misbuild (resources/fonts/ not staged) -- diagnosable, but Aptx still works via the fallback chain
            ForesterUtil.printWarningMessage( AptxConstants.PRG_NAME,
                                              "no bundled fonts found on the classpath (" + FONT_DIR
                                                      + "); falling back to installed fonts" );
        }
    }

    private FontResources() {
    }
}
