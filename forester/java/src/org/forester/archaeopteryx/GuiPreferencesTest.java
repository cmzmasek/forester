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

import org.forester.archaeopteryx.Options.OVERVIEW_PLACEMENT_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_DISPLAY_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.Options.SUPPORT_VISUALIZATION;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
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

            // the palette round-trip below is only meaningful with >1 palette to switch between; guard so this
            // test fails loudly (rather than silently going vacuous) if the palette set is ever reduced to one
            if ( PropertyColorScheme.paletteNames().size() < 2 ) {
                return fail( "expected at least 2 palettes for a meaningful color_palette round-trip" );
            }

            // round-trip: flip a representative subset to the OPPOSITE of its default, save, reload into a
            // fresh default instance, and require the flipped values to come back
            final Options src = Options.createDefaultInstance();
            final boolean tree_name = !src.isShowTreeName();
            final boolean scale = !src.isShowScale();
            final boolean italic = !src.isUseItalicScientificNames();
            final boolean antialias = !src.isAntialiasPrint();
            final boolean white_bg = !src.isGraphicsExportWhiteBackground(); // the only default-TRUE key -> flips to false
            // non-boolean settings round-trip too: enums (node shape/fill, support viz), a short (node size),
            // a float (branch width) and doubles (support threshold, min-confidence fraction)
            final NodeShape shape = ( src.getDefaultNodeShape() == NodeShape.RECTANGLE ) ? NodeShape.CIRCLE
                    : NodeShape.RECTANGLE;
            final NodeFill fill = ( src.getDefaultNodeFill() == NodeFill.NONE ) ? NodeFill.SOLID : NodeFill.NONE;
            final SUPPORT_VISUALIZATION support_viz = ( src.getSupportVisualization() == SUPPORT_VISUALIZATION.NONE )
                    ? SUPPORT_VISUALIZATION.SIZE_SCALED : SUPPORT_VISUALIZATION.NONE;
            final Options.NODE_AGE_SHAPE age_shape = ( src.getNodeAgeShape() == Options.NODE_AGE_SHAPE.BAR )
                    ? Options.NODE_AGE_SHAPE.SPINDLE : Options.NODE_AGE_SHAPE.BAR;
            final short node_size = (short) ( src.getDefaultNodeShapeSize() + 3 );
            // all numeric flips stay INSIDE the persisted [min,max] (branch width [0.5,20], support/min-conf [0,1],
            // raster scale [1,8]) so they survive the clamp on load; out-of-range clamping is checked separately below
            final float branch_width = src.getDefaultBranchWidth() + 2.5f; // 1 -> 3.5
            final double support_threshold = ( src.getSupportThreshold() > 0.5 ) ? 0.25 : 0.75;
            final double min_conf = ( src.getMinConfidenceFraction() > 0.5 ) ? 0.1 : 0.7;
            // Display-tab layout + export settings: enums, an int, booleans and a validated string (palette)
            final PHYLOGENY_GRAPHICS_TYPE gtype = ( src.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE )
                    ? PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR : PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE;
            final PHYLOGENY_DISPLAY_TYPE dtype = ( src.getPhylogenyDisplayType() == PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM )
                    ? PHYLOGENY_DISPLAY_TYPE.CLADOGRAM : PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM;
            final OVERVIEW_PLACEMENT_TYPE ov = ( src.getOvPlacement() == OVERVIEW_PLACEMENT_TYPE.LOWER_RIGHT )
                    ? OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT : OVERVIEW_PLACEMENT_TYPE.LOWER_RIGHT;
            final Options.TREE_ORIENTATION orient = ( src.getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_TOP )
                    ? Options.TREE_ORIENTATION.ROOT_BOTTOM : Options.TREE_ORIENTATION.ROOT_TOP;
            final Options.TIP_LABEL_DIRECTION tip_dir = ( src
                    .getTipLabelDirection() == Options.TIP_LABEL_DIRECTION.VERTICAL )
                            ? Options.TIP_LABEL_DIRECTION.ANGLED : Options.TIP_LABEL_DIRECTION.VERTICAL;
            final Options.NODE_LABEL_DIRECTION node_dir = ( src
                    .getNodeLabelDirection() == Options.NODE_LABEL_DIRECTION.HORIZONTAL )
                            ? Options.NODE_LABEL_DIRECTION.RADIAL : Options.NODE_LABEL_DIRECTION.HORIZONTAL;
            final Options.FOUND_COLOR found_color = ( src.getFoundColor() == Options.FOUND_COLOR.ELECTRIC_VIOLET )
                    ? Options.FOUND_COLOR.NEON_MAGENTA : Options.FOUND_COLOR.ELECTRIC_VIOLET;
            final String palette = altPalette( src.getColorPaletteName() ); // a real, non-current palette name
            final int raster_scale = src.getRasterExportScale() + 2;
            final boolean transparent = !src.isTransparentExportBackground();
            final boolean visible_only = !src.isGraphicsExportVisibleOnly();
            src.setShowTreeName( tree_name );
            src.setShowScale( scale );
            src.setUseItalicScientificNames( italic );
            src.setAntialiasPrint( antialias );
            src.setGraphicsExportWhiteBackground( white_bg );
            src.setDefaultNodeShape( shape );
            src.setDefaultNodeShapeSize( node_size );
            src.setDefaultNodeFill( fill );
            src.setDefaultBranchWidth( branch_width );
            src.setSupportVisualization( support_viz );
            src.setNodeAgeShape( age_shape );
            src.setSupportThreshold( support_threshold );
            src.setMinConfidenceFraction( min_conf );
            src.setPhylogenyGraphicsType( gtype );
            src.setPhylogenyDisplayType( dtype );
            src.setOvPlacement( ov );
            src.setTreeOrientation( orient );
            src.setTipLabelDirection( tip_dir );
            src.setNodeLabelDirection( node_dir );
            src.setFoundColor( found_color );
            src.setColorPaletteName( palette );
            src.setRasterExportScale( raster_scale );
            src.setTransparentExportBackground( transparent );
            src.setGraphicsExportVisibleOnly( visible_only );
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
            if ( dst.getDefaultNodeFill() != fill ) {
                return fail( "default_node_fill did not round-trip" );
            }
            if ( dst.getDefaultBranchWidth() != branch_width ) {
                return fail( "default_branch_width did not round-trip" );
            }
            if ( dst.getSupportVisualization() != support_viz ) {
                return fail( "support_visualization did not round-trip" );
            }
            if ( dst.getNodeAgeShape() != age_shape ) {
                return fail( "node_age_shape did not round-trip" );
            }
            if ( dst.getSupportThreshold() != support_threshold ) {
                return fail( "support_threshold did not round-trip" );
            }
            if ( dst.getMinConfidenceFraction() != min_conf ) {
                return fail( "min_confidence_fraction did not round-trip" );
            }
            if ( dst.getPhylogenyGraphicsType() != gtype ) {
                return fail( "phylogeny_graphics_type did not round-trip" );
            }
            if ( dst.getPhylogenyDisplayType() != dtype ) {
                return fail( "phylogeny_display_type did not round-trip" );
            }
            if ( dst.getOvPlacement() != ov ) {
                return fail( "overview_placement did not round-trip" );
            }
            if ( dst.getTipLabelDirection() != tip_dir ) {
                return fail( "tip_label_direction did not round-trip" );
            }
            if ( dst.getFoundColor() != found_color ) {
                return fail( "found_color did not round-trip" );
            }
            if ( dst.getNodeLabelDirection() != node_dir ) {
                return fail( "node_label_direction did not round-trip" );
            }
            if ( dst.getTreeOrientation() != orient ) {
                return fail( "tree_orientation did not round-trip" );
            }
            if ( !palette.equals( dst.getColorPaletteName() ) ) {
                return fail( "color_palette did not round-trip" );
            }
            if ( dst.getRasterExportScale() != raster_scale ) {
                return fail( "raster_export_scale did not round-trip" );
            }
            if ( dst.isTransparentExportBackground() != transparent ) {
                return fail( "transparent_export_background did not round-trip" );
            }
            if ( dst.isGraphicsExportVisibleOnly() != visible_only ) {
                return fail( "graphics_export_visible_only did not round-trip" );
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

            // a stored palette name that is no longer a known palette must be ignored (validated stringPref),
            // leaving the default palette rather than selecting a non-existent one
            final Path stale = dir.resolve( "stale-palette.properties" );
            Files.write( stale, "color_palette=NoSuchPalette\n".getBytes( StandardCharsets.UTF_8 ) );
            final Options sp = Options.createDefaultInstance();
            final String palette_default = sp.getColorPaletteName();
            new GuiPreferences( stale ).applyTo( sp );
            if ( !palette_default.equals( sp.getColorPaletteName() ) ) {
                return fail( "an unknown palette name must be ignored (default kept)" );
            }

            // an OUT-OF-RANGE or non-finite numeric value must be CLAMPED to the setting's [min,max], never applied
            // verbatim -- otherwise it would later throw when it seeds a range-limited SpinnerNumberModel, or feed
            // NaN/huge values into paint/export. (branch width [0.5,20], support/min-conf [0,1], raster scale [1,8],
            // node size [0,100].)
            final Path oor = dir.resolve( "out-of-range.properties" );
            Files.write( oor, ( "default_branch_width=1000\nsupport_threshold=5.0\nmin_confidence_fraction=-2\n"
                    + "raster_export_scale=999\ndefault_node_shape_size=-7\nsupport_threshold_nan=ignored\n" )
                    .getBytes( StandardCharsets.UTF_8 ) );
            final Options c = Options.createDefaultInstance();
            new GuiPreferences( oor ).applyTo( c );
            if ( c.getDefaultBranchWidth() != 20f ) {
                return fail( "branch width 1000 must clamp to 20, got " + c.getDefaultBranchWidth() );
            }
            if ( c.getSupportThreshold() != 1.0 ) {
                return fail( "support threshold 5.0 must clamp to 1.0, got " + c.getSupportThreshold() );
            }
            if ( c.getMinConfidenceFraction() != 0.0 ) {
                return fail( "min-confidence -2 must clamp to 0.0, got " + c.getMinConfidenceFraction() );
            }
            if ( c.getRasterExportScale() != 8 ) {
                return fail( "raster scale 999 must clamp to 8, got " + c.getRasterExportScale() );
            }
            if ( c.getDefaultNodeShapeSize() != 0 ) {
                return fail( "node size -7 must clamp to 0, got " + c.getDefaultNodeShapeSize() );
            }
            // NaN/Infinity must not reach Options either: they order outside any finite range, so they clamp to a bound
            final Path nan = dir.resolve( "nan.properties" );
            Files.write( nan, "support_threshold=NaN\ndefault_branch_width=Infinity\n"
                    .getBytes( StandardCharsets.UTF_8 ) );
            final Options n = Options.createDefaultInstance();
            new GuiPreferences( nan ).applyTo( n );
            if ( !Double.isFinite( n.getSupportThreshold() ) || !Float.isFinite( n.getDefaultBranchWidth() ) ) {
                return fail( "NaN/Infinity must not reach Options: threshold=" + n.getSupportThreshold()
                        + " branchWidth=" + n.getDefaultBranchWidth() );
            }

            // the radial tree styles ARE restored now: since 0.11.7 a circular branch-length tree is a real
            // circular phylogram, so a persisted CIRCULAR (or UNROOTED) reopens exactly as the user left it
            final Path radial = dir.resolve( "radial.properties" );
            Files.write( radial, "phylogeny_graphics_type=CIRCULAR\n".getBytes( StandardCharsets.UTF_8 ) );
            final Options r = Options.createDefaultInstance();
            new GuiPreferences( radial ).applyTo( r );
            if ( r.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                return fail( "a persisted CIRCULAR graphics type must be restored, got "
                        + r.getPhylogenyGraphicsType() );
            }
            final Path unrooted = dir.resolve( "unrooted.properties" );
            Files.write( unrooted, "phylogeny_graphics_type=UNROOTED\n".getBytes( StandardCharsets.UTF_8 ) );
            final Options ur = Options.createDefaultInstance();
            new GuiPreferences( unrooted ).applyTo( ur );
            if ( ur.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
                return fail( "a persisted UNROOTED graphics type must be restored, got "
                        + ur.getPhylogenyGraphicsType() );
            }
            // ...as is a rectangular-family style
            final Path euro = dir.resolve( "euro.properties" );
            Files.write( euro, "phylogeny_graphics_type=EURO_STYLE\n".getBytes( StandardCharsets.UTF_8 ) );
            final Options eu = Options.createDefaultInstance();
            new GuiPreferences( euro ).applyTo( eu );
            if ( eu.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                return fail( "a persisted EURO_STYLE must be restored, got " + eu.getPhylogenyGraphicsType() );
            }

            // deleteSettingsFile removes the persisted file (so a "Reset to Defaults" survives a restart) and is a
            // silent no-op when the file is already gone
            final Path to_delete = dir.resolve( "to-delete.properties" );
            new GuiPreferences( to_delete ).saveFrom( Options.createDefaultInstance() );
            if ( !Files.exists( to_delete ) ) {
                return fail( "saveFrom should have written the file to delete" );
            }
            new GuiPreferences( to_delete ).deleteSettingsFile();
            if ( Files.exists( to_delete ) ) {
                return fail( "deleteSettingsFile must remove the settings file" );
            }
            new GuiPreferences( to_delete ).deleteSettingsFile(); // second delete: silent no-op, must not throw

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

    /** A real palette name (see {@link PropertyColorScheme#paletteNames()}) different from {@code current}, so a
     *  round-trip actually changes the value. Falls back to {@code current} if only one palette exists. */
    private static String altPalette( final String current ) {
        for ( final String n : PropertyColorScheme.paletteNames() ) {
            if ( !n.equals( current ) ) {
                return n;
            }
        }
        return current;
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
