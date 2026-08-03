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

import org.forester.phylogeny.data.NodeVisualData.NodeShape;

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
            final boolean white_bg = !src.isGraphicsExportWhiteBackground(); // the only default-TRUE key -> flips to false
            // non-boolean settings round-trip too: an enum (node shape) and a number (node size)
            final NodeShape shape = ( src.getDefaultNodeShape() == NodeShape.RECTANGLE ) ? NodeShape.CIRCLE
                    : NodeShape.RECTANGLE;
            final short node_size = (short) ( src.getDefaultNodeShapeSize() + 3 );
            src.setShowTreeName( tree_name );
            src.setShowScale( scale );
            src.setUseItalicScientificNames( italic );
            src.setAntialiasPrint( antialias );
            src.setGraphicsExportWhiteBackground( white_bg );
            src.setDefaultNodeShape( shape );
            src.setDefaultNodeShapeSize( node_size );
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
            if ( dst.isGraphicsExportWhiteBackground() != white_bg ) {
                return fail( "graphics_export_white_background did not round-trip" );
            }
            if ( dst.getDefaultNodeShape() != shape ) {
                return fail( "default_node_shape did not round-trip" );
            }
            if ( dst.getDefaultNodeShapeSize() != node_size ) {
                return fail( "default_node_shape_size did not round-trip" );
            }

            // a key absent from the file must leave that option at its current (default) value
            final Path partial = dir.resolve( "partial.properties" );
            Files.write( partial, "show_scale=true\n".getBytes( StandardCharsets.UTF_8 ) );
            final Options e = Options.createDefaultInstance();
            final boolean name_default = e.isShowTreeName(); // NOT present in the partial file
            // upgrade path: a file written by an older version has no graphics_export_white_background key; that
            // default-TRUE option must stay TRUE (not silently flip to false) when the key is absent
            final boolean white_default = e.isGraphicsExportWhiteBackground();
            // the node shape/size keys are likewise absent -> must keep their built-in defaults
            final NodeShape shape_default = e.getDefaultNodeShape();
            final short size_default = e.getDefaultNodeShapeSize();
            new GuiPreferences( partial ).applyTo( e );
            if ( !e.isShowScale() ) {
                return fail( "show_scale=true in the file should have been applied" );
            }
            if ( e.isShowTreeName() != name_default ) {
                return fail( "a key absent from the file must not change that option" );
            }
            if ( e.isGraphicsExportWhiteBackground() != white_default ) {
                return fail( "an absent default-TRUE key (upgrade path) must keep its default (TRUE)" );
            }
            if ( ( e.getDefaultNodeShape() != shape_default ) || ( e.getDefaultNodeShapeSize() != size_default ) ) {
                return fail( "absent node shape/size keys must keep their defaults" );
            }

            // a corrupt/unknown stored value must be ignored (never throw), leaving the option at its default
            final Path bad = dir.resolve( "bad.properties" );
            Files.write( bad, "default_node_shape=NOT_A_SHAPE\ndefault_node_shape_size=xyz\n"
                    .getBytes( StandardCharsets.UTF_8 ) );
            final Options b = Options.createDefaultInstance();
            final NodeShape b_shape = b.getDefaultNodeShape();
            final short b_size = b.getDefaultNodeShapeSize();
            new GuiPreferences( bad ).applyTo( b );
            if ( ( b.getDefaultNodeShape() != b_shape ) || ( b.getDefaultNodeShapeSize() != b_size ) ) {
                return fail( "a corrupt node shape/size value must be ignored (default kept)" );
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
