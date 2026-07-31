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

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Unit test for {@link GuiPreferences}: display toggles must round-trip through the settings file, a key
 * absent from the file must leave that option at its default (not clobber it), and a missing file must be a
 * silent no-op. Pure disk I/O against a throwaway temp file -- runs headless, never touches the real
 * {@code ~/.archaeopteryx}.
 */
public final class GuiPreferencesTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "GuiPreferences: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        Path dir = null;
        try {
            dir = Files.createTempDirectory( "aptx-guiprefs" );
            final Path file = dir.resolve( "display-settings.properties" );

            // round-trip: flip a representative subset to the OPPOSITE of its default, save, reload into a
            // fresh default instance, and require the flipped values to come back
            final Options src = Options.createDefaultInstance();
            final boolean tree_name = !src.isShowTreeName();
            final boolean scale = !src.isShowScale();
            final boolean italic = !src.isUseItalicScientificNames();
            final boolean antialias = !src.isAntialiasPrint();
            src.setShowTreeName( tree_name );
            src.setShowScale( scale );
            src.setUseItalicScientificNames( italic );
            src.setAntialiasPrint( antialias );
            new GuiPreferences( file ).saveFrom( src );
            if ( !Files.exists( file ) ) {
                return fail( "saveFrom did not write the settings file" );
            }
            final Options dst = Options.createDefaultInstance();
            new GuiPreferences( file ).applyTo( dst );
            if ( dst.isShowTreeName() != tree_name ) {
                return fail( "show_tree_name did not round-trip" );
            }
            if ( dst.isShowScale() != scale ) {
                return fail( "show_scale did not round-trip" );
            }
            if ( dst.isUseItalicScientificNames() != italic ) {
                return fail( "use_italic_scientific_names did not round-trip" );
            }
            if ( dst.isAntialiasPrint() != antialias ) {
                return fail( "antialias_print did not round-trip" );
            }

            // a key absent from the file must leave that option at its current (default) value
            final Path partial = dir.resolve( "partial.properties" );
            Files.write( partial, "show_scale=true\n".getBytes( StandardCharsets.UTF_8 ) );
            final Options e = Options.createDefaultInstance();
            final boolean name_default = e.isShowTreeName(); // NOT present in the partial file
            new GuiPreferences( partial ).applyTo( e );
            if ( !e.isShowScale() ) {
                return fail( "show_scale=true in the file should have been applied" );
            }
            if ( e.isShowTreeName() != name_default ) {
                return fail( "a key absent from the file must not change that option" );
            }

            // a missing file must be a silent no-op (defaults untouched, no exception)
            final Options m = Options.createDefaultInstance();
            final boolean m_name = m.isShowTreeName();
            new GuiPreferences( dir.resolve( "does-not-exist.properties" ) ).applyTo( m );
            if ( m.isShowTreeName() != m_name ) {
                return fail( "a missing settings file must leave defaults untouched" );
            }
            return true;
        }
        catch ( final Throwable ex ) {
            ex.printStackTrace();
            return false;
        }
        finally {
            deleteQuietly( dir );
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [GuiPreferencesTest] " + msg );
        return false;
    }

    private static void deleteQuietly( final Path dir ) {
        if ( dir == null ) {
            return;
        }
        try ( final java.util.stream.Stream<Path> paths = Files.walk( dir ) ) {
            paths.sorted( java.util.Comparator.reverseOrder() ).forEach( p -> {
                try {
                    Files.deleteIfExists( p );
                }
                catch ( final Exception ignored ) {
                }
            } );
        }
        catch ( final Exception ignored ) {
        }
    }

    private GuiPreferencesTest() {
    }
}
