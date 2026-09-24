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

import java.awt.Component;
import java.awt.Container;
import java.awt.GraphicsEnvironment;
import java.nio.file.Files;
import java.nio.file.Path;

import javax.swing.AbstractButton;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.OVERVIEW_PLACEMENT_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_DISPLAY_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.Options.SUPPORT_VISUALIZATION;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * Tests "Reset to Defaults". Headless: {@link Options#resetToDefaults} returns a mutated Options to exactly the
 * built-in defaults (guards {@code init()} completeness). Headful: {@link MainFrame#resetToDefaults} resets the
 * live Options in place, RE-SEEDS the menu controls (so a later {@code updateOptions} reads defaults, not the
 * flipped state), reverts the theme to light, deletes the persisted settings file, and resets the per-tab palette.
 * The GUI part isolates the settings directory to a temp dir, so it never touches the real {@code ~/.archaeopteryx}.
 */
public final class ResetToDefaultsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ResetToDefaults: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( !optionsResetOk() ) {
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return guiResetOk();
    }

    /** {@link Options#resetToDefaults} must return a mutated Options to the SAME state a fresh
     *  {@link Options#createInstance} produces -- across the persisted + display fields. A field that init() forgets
     *  would come back non-default here. */
    private static boolean optionsResetOk() {
        // Christian, 2026-09-23: support values are DATA. With NONE a Save As Nexus or Newick dropped them
        // without saying so, and the phyloXML <-> Nexus round trips are lossless only with brackets. Pinned
        // here because the reflective comparison below cannot catch it: it only checks that a reset RESTORES
        // the default, so a silent change OF the default would pass unnoticed.
        if ( Options.createInstance()
                .getNhConversionSupportValueStyle() != PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE.IN_SQUARE_BRACKETS ) {
            System.out.println( "  [ResetToDefaultsTest] the default NH/Nexus support style is no longer "
                    + "IN_SQUARE_BRACKETS, so saving a tree drops its support values: "
                    + Options.createInstance().getNhConversionSupportValueStyle() );
            return false;
        }
        final Options o = Options.createInstance();
        // drive a broad set of fields away from their defaults
        o.setShowScale( true );
        o.setShowTreeName( false );
        o.setAntialiasExport( false );
        o.setShowOverview( false );
        o.setUseItalicScientificNames( false );
        o.setGraphicsExportWhiteBackground( false );
        o.setTransparentExportBackground( true );
        o.setGraphicsExportVisibleOnly( true );
        o.setRasterExportScale( 8 );
        o.setSupportThreshold( 0.1 );
        o.setMinConfidenceFraction( 0.9 );
        o.setDefaultBranchWidth( 7f );
        o.setDefaultNodeShape( NodeShape.RECTANGLE );
        o.setDefaultNodeFill( NodeFill.NONE );
        o.setDefaultNodeShapeSize( (short) 42 );
        o.setOvPlacement( OVERVIEW_PLACEMENT_TYPE.LOWER_RIGHT );
        o.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE );
        o.setPhylogenyDisplayType( PHYLOGENY_DISPLAY_TYPE.CLADOGRAM ); // default is UNALIGNED_PHYLOGRAM
        o.setColorPaletteName( "Colorblind-friendly" );
        o.setShowScaleGrid( true );
        o.setShowScaleAxis( true );
        o.setShowHpdBars( true );
        o.setShowFossilRangeBars( true );
        o.setShowZebraStripes( true );
        o.setBreakLongBranches( true );
        o.setShowInternalTaxonomyKey( true );
        o.setTipLabelsBelowColumns( true );
        o.setReserveLegendColumn( false ); // default is ON, so drive to OFF to prove reset restores it
        o.setReverseTipOrder( true );
        o.setBoldFoundLabels( true );
        o.setDimNonMatches( false );   // default is ON, so drive to OFF to prove reset restores it
        o.setPulseFoundNodes( false ); // default is ON, so drive to OFF to prove reset restores it
        o.setCheckForUpdatesAtLaunch( true ); // default OFF (opt-in)
        o.setAutoColorNewTrees( false ); // default is ON (JS-parity auto-color), so drive OFF
        o.setAbbreviateScientificTaxonNames( true );
        o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP ); // default is ROOT_LEFT (the DEFAULT new tabs get)
        o.setTipLabelDirection( Options.TIP_LABEL_DIRECTION.HORIZONTAL ); // default is VERTICAL
        o.setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.RADIAL ); // "Radial Labels"; default is HORIZONTAL
        o.setFoundColor( Options.FOUND_COLOR.NEON_MAGENTA ); // default is ELECTRIC_VIOLET
        o.setSupportVisualization( SUPPORT_VISUALIZATION.SIZE_SCALED );
        o.setNodeAgeShape( Options.NODE_AGE_SHAPE.SPINDLE ); // default is BAR
        // default is AUTO (silently promote numeric internal labels that all look like support)
        o.setConfidenceFromInternalLabels( Options.CONFIDENCE_FROM_INTERNAL_LABELS.NEVER );
        o.setDomainLabelMode( Options.DOMAIN_LABEL_MODE.LEGEND ); // default is ON_DOMAINS
        o.setShowDomainGlow( true ); // default is false
        o.setShowMsa( true ); // default is false
        o.setMsaColumnWidth( AptxConstants.MSA_COLUMN_WIDTH_MAX ); // default is MSA_COLUMN_WIDTH_DEFAULT
        o.setShowMsaConservation( false ); // default is true
        o.setMsaConservationMeasure( MsaConservation.Measure.INFORMATION ); // default is IDENTITY
        // "Export at a fixed size" (default off / MILLIMETERS / 170 x 120 mm / 300 DPI)
        o.setExportUseFixedSize( true );
        o.setExportSizeUnit( ExportSizeSpec.Unit.PIXELS );
        o.setExportWidth( 3000 );
        o.setExportHeight( 2000 );
        o.setExportDpi( 600 );
        // search options (reset by init, resynced onto the control-panel checkboxes by Reset to Defaults)
        o.setSearchCaseSensitive( true );
        o.setMatchWholeTermsOnly( true );
        o.setSearchWithRegex( true );
        o.setInverseSearchResult( true );
        o.setSearchProperties( false );

        // the 18 the fixture used to leave at their defaults, so the comparison after the reset said nothing about
        // them (found by everyFieldMoved below, 2026-09-21)
        o.setOutlineFontsInVectorExport( false );
        o.setInternalLabelsAboveBranch( false );
        o.setColorLabelsSameAsParentBranch( true );
        o.setShowDefaultNodeShapesInternal( true );
        o.setShowDefaultNodeShapesExternal( true );
        o.setShowDefaultNodeShapesForMarkedNodes( true );
        o.setReplaceUnderscoresInNhParsing( true );
        o.setAllowErrorsInDistanceToParent( true );
        o.setParseBeastStyleExtendedNexusTags( false );
        // must differ from the default, which is IN_SQUARE_BRACKETS since 2026-09-23
        o.setNhConversionSupportValueStyle( PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE.AS_INTERNAL_NODE_NAMES );
        o.setCladogramType( Options.CLADOGRAM_TYPE.NON_LINED_UP );
        o.setShowConfidenceStddev( true );
        o.setShowMadConfidence( true );
        o.setExportBlackAndWhite( true );
        o.setPdfLineWidth( 3.5f );
        o.setShowTipImages( true );
        o.setTipImageSize( 99 );
        o.setBaseFont( new java.awt.Font( "Serif", java.awt.Font.BOLD, 21 ) );
        // ...and the three that have NO mutator at all -- init() is their only writer, so a setter cannot move them
        // and only reflection can. They are the fields this test most needs to drive: nothing else in the program
        // would notice if init() stopped assigning them.
        if ( !driveFieldsWithoutSetters( o ) ) {
            return false;
        }

        if ( !everyFieldMoved( o, Options.createInstance() ) ) {
            return false;
        }

        o.resetToDefaults();

        final Options def = Options.createInstance();
        // both halves: every field by reflection (complete, and covers a field added later), and the named list
        // below, which says what each one IS and carries the few value-specific invariants
        return everyFieldReset( o, def ) & sameDefaults( o, def );
    }

    /**
     * The three Options fields with no mutator -- {@code _editable},
     * {@code _number_of_digits_after_comma_for_confidence_values} and {@code _taxonomy_extraction} are assigned in
     * {@code init()} and nowhere else. Driven here by reflection so that {@link #sameDefaults} actually proves
     * something about them; skipping them instead would leave exactly the hole this guard exists to close. If one
     * ever gains a setter, use it -- the field write below would still pass, but the setter would go untested.
     */
    private static boolean driveFieldsWithoutSetters( final Options o ) {
        try {
            set( o, "_editable", Boolean.FALSE );
            set( o, "_number_of_digits_after_comma_for_confidence_values", Short.valueOf( (short) 5 ) );
            set( o, "_taxonomy_extraction", org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION.AGGRESSIVE );
            return true;
        }
        catch ( final Exception e ) {
            System.out.println( "  [ResetToDefaultsTest] could not drive a setter-less field: " + e
                    + " -- if the field was renamed or retyped, fix the name here rather than dropping the check" );
            return false;
        }
    }

    private static void set( final Options o, final String field, final Object value ) throws Exception {
        final java.lang.reflect.Field f = Options.class.getDeclaredField( field );
        f.setAccessible( true );
        f.set( o, value );
    }

    /**
     * Every field of {@link Options} must actually DIFFER from a fresh default at this point -- i.e. the drive-away
     * above really took. Without this the comparison after the reset passes on a value that never moved: a setter
     * that ignores its argument, a getter that reads a different field, or a fixture line that was never typed all
     * look identical to "reset worked". (Measured 2026-09-19: with {@code setReserveLegendColumn} stubbed to ignore
     * its argument, this whole test still passed.)
     * <p>
     * Reflection over Options' own readers, not a hand-kept list, so a field added LATER is covered without anyone
     * remembering to come here: it fails until the fixture drives it or it is named in {@link #NOT_DRIVEN} with a
     * reason. That also makes THIS the check that gives the one after the reset its teeth -- "init() forgot a field"
     * can only be caught for a field that was moved off its default first.
     */
    private static boolean everyFieldMoved( final Options driven, final Options fresh ) {
        boolean ok = true;
        int checked = 0;
        final java.util.List<String> names = new java.util.ArrayList<String>();
        for( final java.lang.reflect.Method m : Options.class.getDeclaredMethods() ) {
            if ( !isReader( m ) ) {
                continue;
            }
            final String n = m.getName();
            names.add( n );
            try {
                m.setAccessible( true );
                final Object a = m.invoke( driven );
                final Object b = m.invoke( fresh );
                final boolean moved = ( a == null ) ? ( b != null ) : !a.equals( b );
                ++checked;
                if ( !moved ) {
                    System.out.println( "  [ResetToDefaultsTest] the fixture never moves " + n + "() off its default ("
                            + a + "), so the reset check below proves nothing about it -- drive it above (via its"
                            + " setter, or driveFieldsWithoutSetters if it has none)" );
                    ok = false;
                }
            }
            catch ( final Exception e ) {
                System.out.println( "  [ResetToDefaultsTest] could not read " + n + "(): " + e );
                ok = false; // a reader that could not run has proved nothing, and must not pass for one that did
            }
        }
        if ( checked < 40 ) {
            System.out.println( "  [ResetToDefaultsTest] only " + checked + " Options readers were found (" + names
                    + ") -- the reflection scan is not seeing the class" );
            ok = false;
        }
        return ok;
    }

    /**
     * The mirror of {@link #everyFieldMoved}: after the reset EVERY field must be back at a fresh default's value.
     * {@link #sameDefaults} below compares a hand-kept list of names, which is good documentation and carries a few
     * value-specific invariants, but a field it forgets is a field {@code init()} may forget too -- so the complete
     * comparison is done here, by reflection, and cannot fall behind the class.
     */
    private static boolean everyFieldReset( final Options reset, final Options fresh ) {
        boolean ok = true;
        for( final java.lang.reflect.Method m : Options.class.getDeclaredMethods() ) {
            if ( !isReader( m ) ) {
                continue;
            }
            try {
                m.setAccessible( true );
                final Object a = m.invoke( reset );
                final Object b = m.invoke( fresh );
                if ( ( a == null ) ? ( b != null ) : !a.equals( b ) ) {
                    System.out.println( "  [ResetToDefaultsTest] " + m.getName() + "() was not reset: got " + a
                            + ", a fresh default has " + b + " -- init() does not restore it" );
                    ok = false;
                }
            }
            catch ( final Exception e ) {
                System.out.println( "  [ResetToDefaultsTest] could not read " + m.getName() + "(): " + e );
                ok = false;
            }
        }
        return ok;
    }

    /** A zero-argument {@code isXxx()} / {@code getXxx()} value reader on Options. */
    private static boolean isReader( final java.lang.reflect.Method m ) {
        return ( m.getParameterCount() == 0 ) && ( m.getReturnType() != void.class )
                && !java.lang.reflect.Modifier.isStatic( m.getModifiers() )
                && ( m.getName().startsWith( "is" ) || m.getName().startsWith( "get" ) );
    }

    /** Compares the persisted + display fields of a reset Options against a fresh default. */
    private static boolean sameDefaults( final Options o, final Options def ) {
        boolean ok = true;
        ok &= eq( "autoColorNewTrees", o.isAutoColorNewTrees(), def.isAutoColorNewTrees() );
        ok &= eq( "showScale", o.isShowScale(), def.isShowScale() );
        ok &= eq( "showTreeName", o.isShowTreeName(), def.isShowTreeName() );
        ok &= eq( "antialiasExport", o.isAntialiasExport(), def.isAntialiasExport() );
        ok &= eq( "showOverview", o.isShowOverview(), def.isShowOverview() );
        ok &= eq( "italicNames", o.isUseItalicScientificNames(), def.isUseItalicScientificNames() );
        ok &= eq( "whiteBg", o.isGraphicsExportWhiteBackground(), def.isGraphicsExportWhiteBackground() );
        ok &= eq( "transparentBg", o.isTransparentExportBackground(), def.isTransparentExportBackground() );
        ok &= eq( "visibleOnly", o.isGraphicsExportVisibleOnly(), def.isGraphicsExportVisibleOnly() );
        ok &= eq( "rasterScale", o.getRasterExportScale(), def.getRasterExportScale() );
        ok &= eq( "supportThreshold", o.getSupportThreshold(), def.getSupportThreshold() );
        ok &= eq( "minConfidence", o.getMinConfidenceFraction(), def.getMinConfidenceFraction() );
        ok &= eq( "branchWidth", o.getDefaultBranchWidth(), def.getDefaultBranchWidth() );
        ok &= eq( "nodeShape", o.getDefaultNodeShape(), def.getDefaultNodeShape() );
        ok &= eq( "nodeFill", o.getDefaultNodeFill(), def.getDefaultNodeFill() );
        ok &= eq( "nodeSize", o.getDefaultNodeShapeSize(), def.getDefaultNodeShapeSize() );
        ok &= eq( "ovPlacement", o.getOvPlacement(), def.getOvPlacement() );
        ok &= eq( "graphicsType", o.getPhylogenyGraphicsType(), def.getPhylogenyGraphicsType() );
        ok &= eq( "displayType", o.getPhylogenyDisplayType(), def.getPhylogenyDisplayType() );
        ok &= eq( "treeOrientation", o.getTreeOrientation(), def.getTreeOrientation() );
        ok &= eq( "tipLabelDirection", o.getTipLabelDirection(), def.getTipLabelDirection() );
        ok &= eq( "nodeLabelDirection", o.getNodeLabelDirection(), def.getNodeLabelDirection() );
        ok &= eq( "foundColor", o.getFoundColor(), def.getFoundColor() );
        ok &= eq( "foundColor default VIOLET", def.getFoundColor(), Options.FOUND_COLOR.ELECTRIC_VIOLET );
        ok &= eq( "palette", o.getColorPaletteName(), def.getColorPaletteName() );
        ok &= eq( "showScaleGrid", o.isShowScaleGrid(), def.isShowScaleGrid() );
        ok &= eq( "showScaleAxis", o.isShowScaleAxis(), def.isShowScaleAxis() );
        ok &= eq( "showHpdBars", o.isShowHpdBars(), def.isShowHpdBars() );
        ok &= eq( "showFossilRangeBars", o.isShowFossilRangeBars(), def.isShowFossilRangeBars() );
        ok &= eq( "showZebraStripes", o.isShowZebraStripes(), def.isShowZebraStripes() );
        ok &= eq( "breakLongBranches", o.isBreakLongBranches(), def.isBreakLongBranches() );
        ok &= eq( "showInternalTaxonomyKey", o.isShowInternalTaxonomyKey(), def.isShowInternalTaxonomyKey() );
        ok &= eq( "tipLabelsBelowColumns", o.isTipLabelsBelowColumns(), def.isTipLabelsBelowColumns() );
        ok &= eq( "reserveLegendColumn", o.isReserveLegendColumn(), def.isReserveLegendColumn() );
        ok &= eq( "reverseTipOrder", o.isReverseTipOrder(), def.isReverseTipOrder() );
        ok &= eq( "boldFoundLabels", o.isBoldFoundLabels(), def.isBoldFoundLabels() );
        ok &= eq( "dimNonMatches", o.isDimNonMatches(), def.isDimNonMatches() );
        ok &= eq( "pulseFoundNodes", o.isPulseFoundNodes(), def.isPulseFoundNodes() );
        ok &= eq( "checkForUpdatesAtLaunch", o.isCheckForUpdatesAtLaunch(), def.isCheckForUpdatesAtLaunch() );
        ok &= eq( "checkForUpdatesAtLaunch default OFF (opt-in)", def.isCheckForUpdatesAtLaunch(), false );
        // pin the shipped defaults so an accidental revert to OFF is caught (the reset comparison above only checks
        // reset == fresh-default, which would still pass if BOTH flipped)
        ok &= eq( "dimNonMatches default ON", def.isDimNonMatches(), true );
        ok &= eq( "pulseFoundNodes default ON", def.isPulseFoundNodes(), true );
        ok &= eq( "treeOrientation default LEFT", def.getTreeOrientation(), Options.TREE_ORIENTATION.ROOT_LEFT );
        ok &= eq( "tipLabelDirection default VERTICAL", def.getTipLabelDirection(),
                  Options.TIP_LABEL_DIRECTION.VERTICAL );
        ok &= eq( "abbreviateNames", o.isAbbreviateScientificTaxonNames(), def.isAbbreviateScientificTaxonNames() );
        ok &= eq( "supportViz", o.getSupportVisualization(), def.getSupportVisualization() );
        ok &= eq( "nodeAgeShape", o.getNodeAgeShape(), def.getNodeAgeShape() );
        ok &= eq( "confidenceFromInternalLabels", o.getConfidenceFromInternalLabels(),
                  def.getConfidenceFromInternalLabels() );
        ok &= eq( "domainLabelMode", o.getDomainLabelMode(), def.getDomainLabelMode() );
        ok &= eq( "showDomainGlow", o.isShowDomainGlow(), def.isShowDomainGlow() );
        ok &= eq( "showMsa", o.isShowMsa(), def.isShowMsa() );
        ok &= eq( "msaColumnWidth", o.getMsaColumnWidth(), def.getMsaColumnWidth() );
        ok &= eq( "showMsaConservation", o.isShowMsaConservation(), def.isShowMsaConservation() );
        ok &= eq( "msaConservationMeasure", o.getMsaConservationMeasure(), def.getMsaConservationMeasure() );
        ok &= eq( "exportUseFixedSize", o.isExportUseFixedSize(), def.isExportUseFixedSize() );
        ok &= eq( "exportSizeUnit", o.getExportSizeUnit(), def.getExportSizeUnit() );
        ok &= eq( "exportWidth", o.getExportWidth(), def.getExportWidth() );
        ok &= eq( "exportHeight", o.getExportHeight(), def.getExportHeight() );
        ok &= eq( "exportDpi", o.getExportDpi(), def.getExportDpi() );
        // pin the shipped default (off) so an accidental default flip to ON is caught
        ok &= eq( "exportUseFixedSize default OFF", def.isExportUseFixedSize(), false );
        ok &= eq( "searchCase", o.isSearchCaseSensitive(), def.isSearchCaseSensitive() );
        ok &= eq( "searchWholeWords", o.isMatchWholeTermsOnly(), def.isMatchWholeTermsOnly() );
        ok &= eq( "searchRegex", o.isSearchWithRegex(), def.isSearchWithRegex() );
        ok &= eq( "searchInverse", o.isInverseSearchResult(), def.isInverseSearchResult() );
        ok &= eq( "searchProperties", o.isSearchProperties(), def.isSearchProperties() );
        return ok;
    }

    private static boolean guiResetOk() {
        final String prev_dir = System.getProperty( GuiPreferences.DIR_PROPERTY );
        // the theme lives in java.util.prefs (NOT isolated by DIR_PROPERTY); setDarkMode below persists it, so
        // snapshot and restore the developer's real theme choice
        final Configuration.UI prev_ui = new Configuration().getUi();
        Path temp = null;
        try {
            temp = Files.createTempDirectory( "aptx-reset" );
            System.setProperty( GuiPreferences.DIR_PROPERTY, temp.toString() ); // isolate from the real ~/.archaeopteryx
            final Path settings = temp.resolve( GuiPreferences.SETTINGS_FILE );
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, conf, "reset" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final Options o = frame.getOptions();
                    final ControlPanel cp = frame.getMainPanel().getControlPanel();
                    // 1. drive the LIVE state away from defaults across every surface Reset must cover:
                    //    (a) EVERY menu control updateOptions reads, flipped to the opposite of default, so the
                    //        Options->control resync is genuinely completeness-checked (not spot-checked): a control
                    //        missed in applyOptionsToMenuStates stays flipped and fails the 3b probe below
                    flip( frame._abbreviate_scientific_names,
                          frame._use_italic_scientific_names_cbmi, frame._outline_fonts_in_vector_export_cbmi,
                          frame._transparent_export_background_cbmi, frame._graphics_export_white_background_cbmi,
                          frame._color_labels_same_as_parent_branch, frame._show_default_node_shapes_internal_cbmi,
                          frame._show_default_node_shapes_external_cbmi, frame._show_default_node_shapes_for_marked_cbmi,
                          frame._show_scale_grid_cbmi,
                          frame._show_scale_axis_cbmi, frame._show_hpd_bars_cbmi,
                          frame._show_fossil_range_bars_cbmi,
                          frame._show_zebra_stripes_cbmi, frame._break_long_branches_cbmi,
                          frame._tip_labels_below_columns_cbmi,
                          frame._reverse_tip_order_cbmi,
                          frame._bold_found_labels_cbmi, frame._dim_non_matches_cbmi, frame._pulse_found_nodes_cbmi,
                          frame._check_for_updates_cbmi,
                          frame._show_scale_cbmi, frame._show_tree_name_cbmi,
                          frame._show_overview_cbmi, frame._show_confidence_stddev_cbmi, frame._show_mad_confidence_cbmi,
                          frame._antialias_export_cbmi, frame._export_black_and_white_cbmi,
                          frame._replace_underscores_cbmi,
                          frame._allow_errors_in_distance_to_parent_cbmi, frame._graphics_export_visible_only_cbmi,
                          frame._parse_beast_style_extended_nexus_tags_cbmi, frame._label_direction_cbmi,
                          frame._internal_labels_above_branch_rbmi, frame._internal_labels_right_of_node_rbmi,
                          frame._non_lined_up_cladograms_rbmi, frame._ext_node_dependent_cladogram_rbmi,
                          frame._use_brackets_for_conf_in_nh_export_cbmi,
                          frame._use_internal_names_for_conf_in_nh_export_cbmi );
                    frame.setSelectedTypeInTypeMenu( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ); // type menu stale (non-default)
                    //    (b) enum/number Options + a non-default tree style ON THE LIVE TREE + palette + theme
                    o.setRasterExportScale( 7 );
                    o.setSupportThreshold( 0.1 );
                    o.setDefaultNodeShape( NodeShape.RECTANGLE );
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ); // live tree drawn circular
                    cp.setTreeDisplayType( PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM ); // drive P/A/C off the default
                    tp.setColorPaletteName( "Colorblind-friendly" );
                    tp.setTimeAxisType( Options.TIME_AXIS_TYPE.GEOLOGIC ); // per-tab Time-Axis override + toggle
                    tp.setTimeAxisGrid( true );
                    tp.setSizeByPropertyRef( "data:sz" ); // turn "Size by" ON (a per-tab display setting)
                    cp.populateSizeByPropertyBox(); // re-seed the dropdown from the tree, so it shows the active ref
                    if ( !tp.isSizeByProperty() ) {
                        fail( ok, "precondition: Size by should be active before reset" );
                    }
                    if ( !"data:sz".equals( cp.getSizeByPropertySelection() ) ) {
                        fail( ok, "precondition: the 'Size by' dropdown should show data:sz before reset, got "
                                + cp.getSizeByPropertySelection() );
                    }
                    tp.setAncestralPieTrait( "location" ); // turn ancestral-state pies ON (a per-tab display setting)
                    cp.populateAncestralPieBox(); // re-seed the dropdown from the tree, so it shows the active trait
                    if ( !tp.isShowAncestralPies() ) {
                        fail( ok, "precondition: ancestral pies should be active before reset" );
                    }
                    if ( !"location".equals( cp.getAncestralPieSelection() ) ) {
                        fail( ok, "precondition: the 'Ancestral pie' dropdown should show location before reset, got "
                                + cp.getAncestralPieSelection() );
                    }
                    tp.setAnnotationColumns( java.util.Arrays.asList(
                            new AnnotationColumns.ColumnSpec( "data:sz", AnnotationColumns.Type.HEATMAP ) ) );
                    if ( !tp.hasAnnotationColumns() ) {
                        fail( ok, "precondition: annotation columns should be active before reset" );
                    }
                    // clade annotations are per-tab display state exactly like the columns, so a reset owes them too
                    tp.setCladeBands( "genus", TreePanel.CLADE_VIS.BARS, false,
                                      TreePanel.CLADE_LABEL_ANGLE.VERTICAL );
                    // assert on the SPECS, not on placed bands: this fixture's taxonomy does not resolve offline,
                    // so hasCladeBands() would be false and the check below would prove nothing
                    if ( tp.cladeLevelSpecCountForTest() != 1 ) {
                        fail( ok, "precondition: one clade level should be configured before reset, got "
                                + tp.cladeLevelSpecCountForTest() );
                    }
                    // the same chooser assigns which properties go in the tip LABEL, so a reset owes them too
                    tp.setLabelPropertyRefs( java.util.Arrays.asList( "data:sz" ) );
                    tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.MANUAL ); // View > Order Matrix Columns, per tab
                    if ( tp.getMatrixColumnOrder() != MatrixColumnOrder.Mode.MANUAL ) {
                        fail( ok, "precondition: the tab's matrix-column order must have been changed before reset" );
                    }
                    if ( tp.getLabelPropertyRefs() == null ) {
                        fail( ok, "precondition: label property fields should be set before reset" );
                    }
                    frame.setDarkMode( true );
                    //    (c) always-visible ControlPanel controls held STALE: the sun/moon theme toggle + the
                    //        Inverse search checkbox (plus its Options field), which would clobber the reset on the
                    //        next click. In the dark theme the toggle must offer the LIGHT theme (the sun).
                    if ( cp.getThemeToggleButton() == null ) {
                        fail( ok, "the control panel should have a theme toggle button" );
                    }
                    if ( cp.themeToggleOffersDark() ) {
                        fail( ok, "in the dark theme the toggle must offer the light theme (a sun), not the moon" );
                    }
                    final javax.swing.AbstractButton inverse_cb = findButton( cp, MainFrame.INVERSE_SEARCH_RESULT_LABEL );
                    if ( inverse_cb == null ) {
                        fail( ok, "could not find the control-panel inverse search checkbox" );
                    }
                    else {
                        inverse_cb.setSelected( true );
                        o.setInverseSearchResult( true );
                    }
                    // a settings file exists (so we can prove reset deletes it)
                    new GuiPreferences().saveFrom( o );
                    if ( !Files.exists( settings ) ) {
                        fail( ok, "precondition: a settings file should exist before reset" );
                    }

                    // 2. RESET
                    frame.resetToDefaults();

                    // 3a. the live Options reset in place
                    if ( o.getRasterExportScale() != 4 ) {
                        fail( ok, "raster scale should reset to 4, got " + o.getRasterExportScale() );
                    }
                    if ( o.getSupportThreshold() != Options.SUPPORT_THRESHOLD_DEFAULT ) {
                        fail( ok, "support threshold should reset to default, got " + o.getSupportThreshold() );
                    }
                    if ( o.getDefaultNodeShape() != NodeShape.CIRCLE ) {
                        fail( ok, "node shape should reset to CIRCLE, got " + o.getDefaultNodeShape() );
                    }
                    // the per-tab Time-Axis override is cleared -> back to auto-derive (this dateless tree derives NONE)
                    if ( ( tp.effectiveTimeAxisType() != Options.TIME_AXIS_TYPE.NONE ) || tp.isTimeAxisGrid() ) {
                        fail( ok, "Reset must clear the per-tab Time-Axis override (back to auto-derive), got type="
                                + tp.effectiveTimeAxisType() + " grid=" + tp.isTimeAxisGrid() );
                    }
                    // 3b. EVERY flipped menu control was re-seeded: reading them all back via updateOptions must equal
                    //     a fresh default across the whole boolean surface (a control the resync missed stays flipped
                    //     and shows up here)
                    final Options probe = Options.createInstance(); // starts == a fresh default
                    frame.updateOptions( probe ); // overwrites the menu-backed fields FROM the (reset) controls
                    if ( !sameDefaults( probe, Options.createInstance() ) ) {
                        fail( ok, "menu controls were not fully re-seeded to defaults (see field above)" );
                    }
                    // 3c. theme reverted to light AND the ControlPanel sun/moon toggle follows (back to the moon,
                    //     i.e. offering the dark theme again)
                    if ( frame.getConfiguration().getUi() != Configuration.UI.FLAT_LIGHT ) {
                        fail( ok, "theme should reset to FLAT_LIGHT, got " + frame.getConfiguration().getUi() );
                    }
                    if ( !cp.themeToggleOffersDark() ) {
                        fail( ok, "after reset the control-panel theme toggle must offer the dark theme (a moon)" );
                    }
                    // 3d. the always-visible Inverse search checkbox was re-seeded (else the next click re-applies it)
                    if ( ( inverse_cb != null ) && inverse_cb.isSelected() ) {
                        fail( ok, "the 'Inverse' search checkbox must be cleared after reset" );
                    }
                    if ( o.isInverseSearchResult() ) {
                        fail( ok, "Options.inverseSearchResult must reset to false" );
                    }
                    // 3e. the LIVE tree's graphics type was reset (a circular tree must not stay circular)
                    if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) {
                        fail( ok, "the live tree style must reset to RECTANGULAR, got " + tp.getPhylogenyGraphicsType() );
                    }
                    // 3e-bis. the LIVE P/A/C display type was reset too: this (mostly-branch-lengthed) tree resets to
                    // the default UNALIGNED_PHYLOGRAM via the same preferredDisplayTypeForBranchLengthTree policy the
                    // load path uses -- the ALIGNED drive above must not survive
                    if ( cp.getTreeDisplayType() != PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM ) {
                        fail( ok, "the live P/A/C display type must reset to UNALIGNED_PHYLOGRAM for a branch-length "
                                + "tree, got " + cp.getTreeDisplayType() );
                    }
                    // 3f. per-tab palette reset AND "Size by" turned off
                    if ( !PropertyColorScheme.DEFAULT_PALETTE_NAME.equals( tp.getColorPaletteName() ) ) {
                        fail( ok, "per-tab palette should reset to default, got " + tp.getColorPaletteName() );
                    }
                    if ( tp.isSizeByProperty() ) {
                        fail( ok, "Size by must be turned off by reset" );
                    }
                    // the always-visible "Size by" dropdown control must be re-seeded to None too (else it stays
                    // showing "data:sz" while the model is cleared -- the stale-control class this DoD guards against)
                    if ( !"None".equals( cp.getSizeByPropertySelection() ) ) {
                        fail( ok, "the 'Size by' dropdown must be re-seeded to None after reset, got "
                                + cp.getSizeByPropertySelection() );
                    }
                    // 3f-bis. per-tab ancestral pies turned off AND the dropdown re-seeded to None
                    if ( tp.isShowAncestralPies() ) {
                        fail( ok, "ancestral pies must be turned off by reset" );
                    }
                    if ( !"None".equals( cp.getAncestralPieSelection() ) ) {
                        fail( ok, "the 'Ancestral pie' dropdown must be re-seeded to None after reset, got "
                                + cp.getAncestralPieSelection() );
                    }
                    // 3f-ter. per-tab annotation columns cleared (a fresh install has none), and with them the
                    //         per-tab label-field selection -- back to the default of showing every property
                    if ( tp.hasAnnotationColumns() ) {
                        fail( ok, "annotation columns must be cleared by reset" );
                    }
                    if ( tp.getLabelPropertyRefs() != null ) {
                        fail( ok, "the label-field selection must be cleared by reset, got "
                                + tp.getLabelPropertyRefs() );
                    }
                    if ( tp.getMatrixColumnOrder() != MatrixColumnOrder.DEFAULT ) {
                        fail( ok, "the tab's matrix-column order must be reset to Clustered, got "
                                + tp.getMatrixColumnOrder() );
                    }
                    // 3f-quater. per-tab clade annotations cleared. Only assertable when the fixture tree could
                    // actually place bands at that rank -- otherwise there was nothing to clear and nothing to prove.
                    if ( ( tp.cladeLevelSpecCountForTest() != 0 ) || tp.hasCladeBands() ) {
                        fail( ok, "clade annotations must be cleared by reset, still "
                                + tp.cladeLevelSpecCountForTest() + " level(s), ranks " + tp.cladeLevelRanks() );
                    }
                    // 3g. persisted settings file deleted (so the reset survives a restart)
                    if ( Files.exists( settings ) ) {
                        fail( ok, "reset must delete the persisted settings file" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        finally {
            if ( prev_dir == null ) {
                System.clearProperty( GuiPreferences.DIR_PROPERTY );
            }
            else {
                System.setProperty( GuiPreferences.DIR_PROPERTY, prev_dir );
            }
            Configuration.saveUiPreference( prev_ui ); // restore the developer's real theme
            deleteQuietly( temp );
        }
    }

    /** Flips each non-null button to the opposite of its current (default) state, for the reset completeness check. */
    private static void flip( final AbstractButton... buttons ) {
        for ( final AbstractButton b : buttons ) {
            if ( b != null ) {
                b.setSelected( !b.isSelected() );
            }
        }
    }

    /** Depth-first search for the first AbstractButton (checkbox/radio) with the given text, or null. */
    private static AbstractButton findButton( final Container c, final String text ) {
        for ( final Component comp : c.getComponents() ) {
            if ( ( comp instanceof AbstractButton ) && text.equals( ( (AbstractButton) comp ).getText() ) ) {
                return (AbstractButton) comp;
            }
            if ( comp instanceof Container ) {
                final AbstractButton found = findButton( (Container) comp, text );
                if ( found != null ) {
                    return found;
                }
            }
        }
        return null;
    }

    private static boolean eq( final String name, final Object a, final Object b ) {
        if ( ( a == null ) ? ( b != null ) : !a.equals( b ) ) {
            System.out.println( "  [ResetToDefaultsTest] " + name + " not reset to default: got " + a + ", expected "
                    + b );
            return false;
        }
        return true;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ResetToDefaultsTest] " + msg );
        ok[ 0 ] = false;
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        for ( int i = 0; i < 3; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            leaf.setDistanceToParent( 0.1 * ( i + 1 ) ); // branch lengths -> a phylogram, so the reset exercises the
                                                         // branch-length (helper) path of the P/A/C display reset
            // a numeric property (distinct values) so "Size by" can be turned on and its reset verified
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:sz", Integer.toString( i + 1 ), "", "xsd:string", AppliesTo.NODE ) );
            leaf.getNodeData().setProperties( pl );
            root.addAsChild( leaf );
        }
        // a BEAST discrete-trait distribution on the internal root, so the "Ancestral pie" control has a trait
        final PropertiesList root_pl = new PropertiesList();
        root_pl.addProperty(
                new Property( "beast:location_set", "{Africa,Europe,Asia}", "", "xsd:string", AppliesTo.NODE ) );
        root_pl.addProperty(
                new Property( "beast:location_set_prob", "{0.5,0.3,0.2}", "", "xsd:string", AppliesTo.NODE ) );
        root.getNodeData().setProperties( root_pl );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
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

    private ResetToDefaultsTest() {
    }
}
