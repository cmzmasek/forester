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

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.util.Locale;
import java.util.function.Function;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTabbedPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.SpinnerNumberModel;

import org.forester.archaeopteryx.Options.OVERVIEW_PLACEMENT_TYPE;
import org.forester.archaeopteryx.Options.SUPPORT_VISUALIZATION;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.ws.seqdb.NcbiTaxonomyLineageService;
import org.forester.ws.seqdb.TaxonomyCacheStatus;

/**
 * A modeless, live-apply "Settings" dialog that replaces the old Options and Type menus. (The light/dark
 * theme is deliberately NOT here -- it lives on the control panel as a one-click sun/moon toggle, which is
 * faster than a dialog trip for a setting people flip often.) Every control is visible at once -- no cycling. Checkboxes,
 * radio buttons and the tree-style dropdown simply drive the existing menu-item fields via
 * {@code doClick()}, so all of the established apply logic runs unchanged; the few former
 * cycle/chooser items become dropdowns or spinners that set the value on {@link Options} and repaint.
 */
final class SettingsDialog extends JDialog {

    private static final long serialVersionUID = 1L;
    private final MainFrame   _mf;
    // "Taxonomy Cache" tab widgets, refreshed whenever that tab is shown
    private JCheckBox         _cache_enabled_cb;
    private JLabel            _cache_location_label;
    private JLabel            _cache_size_label;
    private JLabel            _cache_status_label;
    private int               _cache_tab_index = -1;
    private JLabel            _current_font_label; // "Fonts, Nodes and Branches" tab
    // per-tab controls (tree style, "Color by" palette, Time Axis) reflect the CURRENT tab; on a tab switch a
    // left-open modeless dialog re-seeds them via these reseeders, under a guard that suppresses their own listeners
    private final java.util.List<Runnable> _tab_reseeders = new java.util.ArrayList<>();
    private boolean           _refreshing = false;
    private boolean           _export_size_adjusting = false; // guards the programmatic W/H update on a unit switch
    private JComboBox<String> _axis_combo; // the per-tab Time-Axis combo ("Auto" + Off/Geologic/Calendar; test hook)

    SettingsDialog( final MainFrame mf ) {
        super( mf, "Settings", false );
        _mf = mf;
        final JTabbedPane tabs = new JTabbedPane();
        // the former single "Display" tab had grown too long, so it is split into three focused tabs
        tabs.addTab( "Layout", scroll( layoutTab() ) );
        tabs.addTab( "Labels & Colors", scroll( labelsColorsTab() ) );
        tabs.addTab( "Overlays", scroll( overlaysTab() ) );
        tabs.addTab( "Fonts, Nodes and Branches", scroll( nodesTab() ) );
        tabs.addTab( "Graphics Export", scroll( exportTab() ) );
        tabs.addTab( "File Reading", scroll( readTab() ) );
        tabs.addTab( "File Saving", scroll( saveTab() ) );
        _cache_tab_index = tabs.getTabCount();
        tabs.addTab( "Taxonomy Cache", scroll( cacheTab() ) );
        // the cache stats are read from disk; refresh them each time the tab is brought to the front
        tabs.addChangeListener( e -> {
            if ( tabs.getSelectedIndex() == _cache_tab_index ) {
                refreshCacheTab();
            }
            refreshFontInfo(); // cheap; keep the "Current font" line in sync with any font change
        } );
        final JButton reset = button( "Reset to Defaults…", this::confirmAndResetToDefaults );
        reset.setToolTipText( "Restore all display settings, the theme, and tree colors to the built-in defaults" );
        final JButton close = new JButton( "Close" );
        close.addActionListener( e -> setVisible( false ) );
        final JPanel south = new JPanel( new BorderLayout() );
        final JPanel south_left = new JPanel( new FlowLayout( FlowLayout.LEFT ) );
        south_left.add( reset );
        final JPanel south_right = new JPanel( new FlowLayout( FlowLayout.RIGHT ) );
        south_right.add( close );
        south.add( south_left, BorderLayout.WEST );
        south.add( south_right, BorderLayout.EAST );
        setLayout( new BorderLayout() );
        add( tabs, BorderLayout.CENTER );
        add( south, BorderLayout.SOUTH );
        pack();
        // pack() sizes to the (narrow) tab CONTENT, which wraps the tab-header row onto two rows and feels cramped;
        // widen the default so all tabs sit on a single row with room to spare (never shrink a naturally-wider pack).
        // The minimum WIDTH keeps it from being dragged back into the wrapped, cramped state; height stays free
        // (each tab already scrolls).
        final int min_width = Math.max( 900, getWidth() );
        setSize( min_width, getHeight() );
        setMinimumSize( new Dimension( min_width, 300 ) );
        setLocationRelativeTo( mf );
    }

    /** Re-seed the per-tab controls (tree style, "Color by" palette, Time Axis type + grid/ages) from the now-current
     *  tab -- called on a main-window tab switch so a left-open modeless dialog reflects the current tree. The guard
     *  suppresses the controls' own listeners so re-seeding doesn't write the value back onto the (new) current tab. */
    void refreshCurrentTabControls() {
        if ( _refreshing ) {
            return;
        }
        _refreshing = true;
        try {
            for ( final Runnable r : _tab_reseeders ) {
                r.run();
            }
        }
        finally {
            _refreshing = false;
        }
    }

    /** Tooltip on the rectangular sub-style dropdown -- it points at where the FAMILY is chosen instead. */
    static final String RECTANGULAR_STYLE_TIP =
            "How a branch joint is drawn in the rectangular layouts. Can be set at any time: while a circular or "
                    + "unrooted tree is showing it takes effect the next time you pick a rectangular layout. "
                    + "Rectangular / circular / unrooted and the root's position are picked on the control panel, "
                    + "top left.";

    /** Label for the per-tab Time-Axis combo's "follow auto-derive" entry (index 0). */
    static final String AXIS_AUTO_LABEL = "Auto (from dates)";

    /** Combo index for a raw type override: 0 == "Auto" (null override), else the enum ordinal + 1 (Off/Geo/Cal). */
    private static int axisComboIndexFor( final Options.TIME_AXIS_TYPE override ) {
        return ( override == null ) ? 0 : ( override.ordinal() + 1 );
    }

    /** The type override a combo index maps to: index 0 -> null ("Auto"), else the enum value at ordinal index-1. */
    private static Options.TIME_AXIS_TYPE axisTypeForComboIndex( final int i ) {
        return ( i <= 0 ) ? null : Options.TIME_AXIS_TYPE.values()[ i - 1 ];
    }

    /** Test hook: the type override the per-tab Time-Axis combo currently shows ({@code null} == the "Auto" entry). */
    Options.TIME_AXIS_TYPE axisComboTypeForTest() {
        return ( _axis_combo == null ) ? null : axisTypeForComboIndex( _axis_combo.getSelectedIndex() );
    }

    /** Test hook: simulate a user selecting a Time-Axis entry ({@code null} == "Auto"), firing the real listener. */
    void userSelectAxisForTest( final Options.TIME_AXIS_TYPE type ) {
        if ( _axis_combo != null ) {
            _axis_combo.setSelectedIndex( axisComboIndexFor( type ) );
        }
    }

    // ---- tabs ----------------------------------------------------------------------------------

    /** Tab 1 of the split former "Display" tab: the tree's layout/orientation and collapsed-subtree/domain
     *  handling -- i.e. the tree's shape, not the decorations on it. (The light/dark theme is NOT here: it is the
     *  control panel's one-click sun/moon toggle, so it needs no dialog trip.) */
    private JPanel layoutTab() {
        final JPanel c = column();
        c.add( header( "Layout" ) );
        // The tip-label angle applies only in a vertical (root-top/bottom) orientation, so grey it out otherwise.
        // Orientation itself is no longer set here -- it is one of the control panel's five layout buttons -- so
        // this combo is re-seeded from the live orientation on every tab refresh instead.
        final JComboBox<Options.TIP_LABEL_DIRECTION> tipLabelCombo = enumCombo( Options.TIP_LABEL_DIRECTION.values(),
                _mf.getOptions().getTipLabelDirection(),
                v -> { _mf.getOptions().setTipLabelDirection( v ); reFitCurrentTree(); } );
        tipLabelCombo.setEnabled( tipLabelAngleApplies() );
        _tab_reseeders.add( () -> tipLabelCombo.setEnabled( tipLabelAngleApplies() ) );
        add( c, labeled( "Rectangular style:", treeStyleCombo() ) );
        add( c, labeled( "Tip label angle:", tipLabelCombo ) );
        addRadioGroup( c, _mf._ext_node_dependent_cladogram_rbmi, _mf._non_lined_up_cladograms_rbmi );
        add( c, cb( _mf._label_direction_cbmi ) );
        add( c, cb( _mf._reverse_tip_order_cbmi ) );
        add( c, cb( _mf._break_long_branches_cbmi ) );
        c.add( header( "Collapsed Subtrees & Domains" ) );
        add( c, cb( _mf._collapsed_with_average_height_cbmi ) );
        add( c, cb( _mf._show_abbreviated_labels_for_collapsed_nodes_cbmi ) );
        add( c, labeled( "Domain labels:", enumCombo( Options.DOMAIN_LABEL_MODE.values(),
                                                      _mf.getOptions().getDomainLabelMode(),
                                                      v -> { _mf.getOptions().setDomainLabelMode( v ); repaintTree(); } ) ) );
        final JCheckBox domain_glow = new JCheckBox( "Domain glow", _mf.getOptions().isShowDomainGlow() );
        domain_glow.addActionListener( e -> { _mf.getOptions().setShowDomainGlow( domain_glow.isSelected() ); repaintTree(); } );
        add( c, domain_glow );
        return c;
    }

    /** Tab 2 of the split former "Display" tab: how node/label TEXT reads and its COLOR, including the
     *  found/selected highlight colour + emphasis. */
    private JPanel labelsColorsTab() {
        final JPanel c = column();
        c.add( header( "Labels" ) );
        c.add( new JLabel( "Internal label placement:" ) );
        addRadioGroup( c, _mf._internal_labels_above_branch_rbmi, _mf._internal_labels_right_of_node_rbmi );
        add( c, cb( _mf._use_italic_scientific_names_cbmi ) );
        add( c, cb( _mf._abbreviate_scientific_names ) );
        add( c, cb( _mf._show_confidence_stddev_cbmi ) );
        add( c, cb( _mf._show_mad_confidence_cbmi ) );
        add( c, cb( _mf._tip_labels_below_columns_cbmi ) );
        c.add( header( "Colors" ) );
        add( c, cb( _mf._color_labels_same_as_parent_branch ) );
        add( c, labeled( "\"Color by\" palette:", paletteCombo() ) );
        c.add( header( "Found / Selected" ) );
        add( c, labeled( "Found/Selected colors:", enumCombo( Options.FOUND_COLOR.values(),
                                                              _mf.getOptions().getFoundColor(),
                                                              v -> { _mf.getOptions().setFoundColor( v );
                                                                     applyFoundColor( v ); } ) ) );
        add( c, cb( _mf._bold_found_labels_cbmi ) );
        add( c, cb( _mf._dim_non_matches_cbmi ) );
        add( c, cb( _mf._pulse_found_nodes_cbmi ) );
        return c;
    }

    /** Tab 3 of the split former "Display" tab: extra graphical elements drawn over/around the tree (scale,
     *  data overlays) plus the viewport chrome (tree name, overview). */
    private JPanel overlaysTab() {
        final JPanel c = column();
        c.add( header( "Scale & Grid" ) );
        add( c, cb( _mf._show_scale_cbmi ) );
        add( c, cb( _mf._show_scale_grid_cbmi ) );
        add( c, cb( _mf._show_scale_axis_cbmi ) );
        c.add( header( "Data Overlays" ) );
        add( c, cb( _mf._show_hpd_bars_cbmi ) );
        add( c, labeled( "Node age shape:", enumCombo( Options.NODE_AGE_SHAPE.values(),
                                                       _mf.getOptions().getNodeAgeShape(),
                                                       v -> { _mf.getOptions().setNodeAgeShape( v ); repaintTree(); } ) ) );
        add( c, cb( _mf._show_fossil_range_bars_cbmi ) );
        add( c, cb( _mf._show_zebra_stripes_cbmi ) );
        add( c, cb( _mf._show_internal_taxonomy_key_cbmi ) );
        // Tip images (raster; PNG/JPG): a picture at each tip, from a local path or http(s) URL in a node property
        // (the Import Annotations / CSV path) or a taxonomy <uri>. Options-direct (no menu item); re-layouts so the
        // tip labels shift to make room. Renders in all five display types (rectangular x3 + circular + unrooted).
        final JCheckBox tip_images = new JCheckBox( "Tip Images", _mf.getOptions().isShowTipImages() );
        tip_images.setToolTipText( "Draw an image at each tip (from a local path or URL in a node property / <uri>). "
                + "Add the paths with Tools → Import Annotations." );
        tip_images.addActionListener( e -> {
            _mf.getOptions().setShowTipImages( tip_images.isSelected() );
            _mf.getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
        } );
        add( c, tip_images );
        add( c, labeled( "Tip image size:", intSpinner( _mf.getOptions().getTipImageSize(), 12, 200, 4, v -> {
            _mf.getOptions().setTipImageSize( v );
            _mf.getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
        } ) ) );
        // Sequence alignment: colored residue cells beside the tips (rectangular root-left only). Options-direct;
        // auto-enabled when the tree carries an aligned molecular sequence (loaded FASTA / phyloXML <mol_seq>).
        final JCheckBox show_msa = new JCheckBox( "Sequence Alignment", _mf.getOptions().isShowMsa() );
        show_msa.setToolTipText( "Show a multiple sequence alignment next to the tree (rectangular root-left). "
                + "Load it with File → Load Alignment (FASTA)." );
        show_msa.addActionListener( e -> {
            _mf.getOptions().setShowMsa( show_msa.isSelected() );
            _mf.getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
        } );
        add( c, show_msa );
        add( c, labeled( "Alignment column width:",
                intSpinner( _mf.getOptions().getMsaColumnWidth(), AptxConstants.MSA_COLUMN_WIDTH_MIN,
                        AptxConstants.MSA_COLUMN_WIDTH_MAX, 1, v -> {
                            _mf.getOptions().setMsaColumnWidth( v );
                            _mf.getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
                        } ) ) );
        // The Time Axis is PER-TREE (per-tab): a SARS-CoV-2 (calendar) tab and a Dinosaur (geologic) tab show
        // different axes at once. These controls act on the CURRENT tab and are seeded from it at open (like the
        // Tree style / palette combos above -- so they reflect the tab that was current when the dialog opened).
        c.add( header( "Time Axis (per tree)" ) );
        // "Auto" (index 0) = follow auto-derive (the tab's type override is null); the other entries are explicit
        // overrides (Off / Geologic / Calendar). The combo shows "Auto" until the user picks an explicit type, and
        // selecting "Auto" again returns the tab to auto-derive.
        final TreePanel cur_ta = _mf.getCurrentTreePanel();
        final String[] axis_labels = { AXIS_AUTO_LABEL, Options.TIME_AXIS_TYPE.NONE.toString(),
                Options.TIME_AXIS_TYPE.GEOLOGIC.toString(), Options.TIME_AXIS_TYPE.CALENDAR.toString() };
        final JComboBox<String> axis_combo = new JComboBox<>( axis_labels );
        axis_combo.setSelectedIndex( axisComboIndexFor( ( cur_ta != null ) ? cur_ta.getTimeAxisTypeOverride() : null ) );
        axis_combo.addActionListener( e -> {
            if ( _refreshing ) {
                return; // a programmatic re-seed on a tab switch -- don't write the value back onto the tab
            }
            final TreePanel cur = _mf.getCurrentTreePanel();
            final Options.TIME_AXIS_TYPE sel = axisTypeForComboIndex( axis_combo.getSelectedIndex() );
            if ( cur != null ) {
                cur.setTimeAxisType( sel ); // sel == null -> "Auto" (follow auto-derive)
            }
            if ( sel == null ) {
                reFitCurrentTree(); // "Auto" needs no manual calibration prompt (a dated tree derives its own)
            }
            else {
                maybeCalibrateAndFit( sel );
            }
        } );
        _axis_combo = axis_combo;
        add( c, labeled( "Time axis:", axis_combo ) );
        _tab_reseeders.add( () -> {
            final TreePanel cur = _mf.getCurrentTreePanel();
            axis_combo.setSelectedIndex( axisComboIndexFor( ( cur != null ) ? cur.getTimeAxisTypeOverride() : null ) );
        } );
        add( c, labeled( "", button( "Set root age… (geologic)", this::setTimeAxisRootAgeDialog ) ) );
        add( c, labeled( "", button( "Set most-recent-tip date… (calendar)", this::setTimeAxisPresentDateDialog ) ) );
        add( c, panelCheckbox( MainFrame.GEOLOGIC_GRID_LABEL, MainFrame.GEOLOGIC_GRID_TIP,
                               TreePanel::isTimeAxisGrid, ( p, b ) -> p.setTimeAxisGrid( b ) ) );
        add( c, panelCheckbox( MainFrame.GEOLOGIC_AGES_LABEL, MainFrame.GEOLOGIC_AGES_TIP,
                               TreePanel::isTimeAxisAges, ( p, b ) -> p.setTimeAxisAges( b ) ) );
        c.add( header( "Tree Name & Overview" ) );
        add( c, cb( _mf._show_tree_name_cbmi ) );
        add( c, cb( _mf._show_overview_cbmi ) );
        add( c, labeled( "Overview placement:", enumCombo( OVERVIEW_PLACEMENT_TYPE.values(),
                                                           _mf.getOptions().getOvPlacement(),
                                                           v -> { _mf.getOptions().setOvPlacement( v );
                                                                  updateOverview(); } ) ) );
        return c;
    }

    private JPanel nodesTab() {
        final JPanel c = column();
        c.add( header( "Node Shapes" ) );
        add( c, cb( _mf._show_default_node_shapes_external_cbmi ) );
        add( c, cb( _mf._show_default_node_shapes_internal_cbmi ) );
        add( c, cb( _mf._show_default_node_shapes_for_marked_cbmi ) );
        add( c, labeled( "Node shape:", enumCombo( NodeShape.values(), _mf.getOptions().getDefaultNodeShape(),
                                                   v -> { _mf.getOptions().setDefaultNodeShape( v ); repaintTree(); },
                                                   SettingsDialog::prettyEnumName ) ) );
        // GRADIENT node fill is retired as a user choice (kept in the NodeFill enum only for phyloXML round-tripping).
        add( c, labeled( "Node fill:", enumCombo( new NodeFill[] { NodeFill.DEFAULT, NodeFill.NONE, NodeFill.SOLID },
                                                  _mf.getOptions().getDefaultNodeFill(),
                                                  v -> { _mf.getOptions().setDefaultNodeFill( v ); repaintTree(); },
                                                  SettingsDialog::prettyEnumName ) ) );
        add( c, labeled( "Node size:", intSpinner( _mf.getOptions().getDefaultNodeShapeSize(), 0, 100, 1,
                                                   v -> { _mf.getOptions().setDefaultNodeShapeSize( v.shortValue() ); repaintTree(); } ) ) );
        c.add( header( "Branches & Confidence" ) );
        add( c, labeled( "Branch width:", doubleSpinner( _mf.getOptions().getDefaultBranchWidth(), 0.5, 20, 0.5,
                                                         v -> { _mf.getOptions().setDefaultBranchWidth( v.floatValue() ); repaintTree(); } ) ) );
        add( c, labeled( "Min. support shown (fraction of scale):",
                         // clamp the seed into [0,1]: SpinnerNumberModel THROWS (not clamps) on an out-of-range
                         // value, so a legacy/persisted absolute value (e.g. 50) must not reach the constructor
                         doubleSpinner( Math.max( 0.0, Math.min( 1.0, _mf.getOptions().getMinConfidenceFraction() ) ),
                                        0, 1.0, 0.05,
                                        v -> { _mf.getOptions().setMinConfidenceFraction( v ); repaintTree(); } ) ) );
        add( c, labeled( "Show support as:", enumCombo( SUPPORT_VISUALIZATION.values(),
                                                        _mf.getOptions().getSupportVisualization(),
                                                        v -> { _mf.getOptions().setSupportVisualization( v ); repaintTree(); } ) ) );
        add( c, labeled( "Support threshold (0–1):", doubleSpinner( _mf.getOptions().getSupportThreshold(), 0, 1.0, 0.05,
                                                                       v -> { _mf.getOptions().setSupportThreshold( v ); repaintTree(); } ) ) );
        c.add( header( "Fonts" ) );
        add( c, labeled( "Tree font:", button( "Choose…", () -> { _mf.chooseFont(); refreshFontInfo(); } ) ) );
        _current_font_label = new JLabel();
        add( c, _current_font_label );
        final JLabel font_blurb = new JLabel( buildPreferredFontsHtml() );
        font_blurb.setAlignmentX( Component.LEFT_ALIGNMENT );
        c.add( font_blurb );
        refreshFontInfo();
        return c;
    }

    private JPanel exportTab() {
        final JPanel c = column();
        c.add( header( "Graphics Export" ) );
        final JSpinner raster_scale = intSpinner( _mf.getOptions().getRasterExportScale(), 1, 8, 1,
                                                  v -> _mf.getOptions().setRasterExportScale( v ) );
        raster_scale.setToolTipText( "<html>Resolution multiplier for raster (PNG/JPG/TIFF) export: the figure is "
                + "re-rendered onto an N&times;-larger canvas for crisp, publication-quality output<br>"
                + "(a true re-render, not pixel doubling). 1 = on-screen size. Higher = larger files/memory; "
                + "very large figures are capped automatically.<br>Does not affect vector (SVG/EPS/PDF) export.</html>" );
        add( c, labeled( "Raster export scale (×):", raster_scale ) );
        add( c, cb( _mf._transparent_export_background_cbmi ) );
        add( c, cb( _mf._graphics_export_white_background_cbmi ) );
        add( c, cb( _mf._outline_fonts_in_vector_export_cbmi ) );
        add( c, cb( _mf._antialias_export_cbmi ) );
        add( c, cb( _mf._export_black_and_white_cbmi ) );
        add( c, cb( _mf._graphics_export_visible_only_cbmi ) );
        add( c, labeled( "PDF line width:", doubleSpinner( _mf.getOptions().getPdfLineWidth(), 0.5, 20, 0.5,
                                                           v -> _mf.getOptions().setPdfLineWidth( v.floatValue() ) ) ) );
        addFixedExportSizeSection( c );
        return c;
    }

    /**
     * "Fixed Export Size": a checkbox that makes every graphics export render at EXACTLY the size below (see
     * {@link ExportSizeSpec}), plus a unit combo, width/height/DPI spinners, journal-column presets and a live
     * "Output:" summary. All Options-direct (no menu item). Switching units converts the numeric width/height so
     * the physical size is preserved; a live summary shows exactly what the file will be.
     */
    private void addFixedExportSizeSection( final JPanel c ) {
        c.add( header( "Fixed Export Size" ) );
        final Options o = _mf.getOptions();
        final JLabel summary = new JLabel();
        final Runnable refresh = () -> summary.setText(
                "<html><div style='width:520px'>Output: " + o.exportSizeSpec().summary() + "</div></html>" );
        // Seed clamped into [MIN,MAX]: a unit conversion (below) can push the stored value out of range, and a
        // SpinnerNumberModel THROWS (not clamps) on an out-of-range seed -- which would break reopening the dialog.
        final JSpinner width = doubleSpinner( clampExportDim( o.getExportWidth() ), AptxConstants.EXPORT_SIZE_DIM_MIN,
                                              AptxConstants.EXPORT_SIZE_DIM_MAX, 1.0, v -> {
                                                  if ( _export_size_adjusting ) {
                                                      return; // a programmatic unit conversion -- Options already set
                                                  }
                                                  o.setExportWidth( v );
                                                  refresh.run();
                                              } );
        final JSpinner height = doubleSpinner( clampExportDim( o.getExportHeight() ), AptxConstants.EXPORT_SIZE_DIM_MIN,
                                               AptxConstants.EXPORT_SIZE_DIM_MAX, 1.0, v -> {
                                                   if ( _export_size_adjusting ) {
                                                       return;
                                                   }
                                                   o.setExportHeight( v );
                                                   refresh.run();
                                               } );
        final JSpinner dpi = intSpinner( o.getExportDpi(), AptxConstants.EXPORT_SIZE_DPI_MIN,
                                         AptxConstants.EXPORT_SIZE_DPI_MAX, 10, v -> {
                                             o.setExportDpi( v );
                                             refresh.run();
                                         } );
        final JComboBox<ExportSizeSpec.Unit> unit = enumCombo( ExportSizeSpec.Unit.values(), o.getExportSizeUnit(),
                new_unit -> {
                    final ExportSizeSpec.Unit old_unit = o.getExportSizeUnit();
                    if ( new_unit != old_unit ) {
                        // preserve the physical size across the unit switch: convert the numeric width/height
                        final ExportSizeSpec cur = new ExportSizeSpec( old_unit, o.getExportWidth(),
                                                                       o.getExportHeight(), o.getExportDpi() );
                        // clamp into the spinner range so a conversion to a higher-multiplier unit (e.g. in->px at a
                        // high DPI) can't store an out-of-range value that would later throw when the dialog reopens
                        final double nw = clampExportDim( ExportSizeSpec.fromInches( cur.widthInches(), new_unit,
                                                                                     o.getExportDpi() ) );
                        final double nh = clampExportDim( ExportSizeSpec.fromInches( cur.heightInches(), new_unit,
                                                                                     o.getExportDpi() ) );
                        o.setExportSizeUnit( new_unit );
                        o.setExportWidth( nw );
                        o.setExportHeight( nh );
                        _export_size_adjusting = true;
                        width.setValue( nw );
                        height.setValue( nh );
                        _export_size_adjusting = false;
                    }
                    refresh.run();
                } );
        final JPanel presets = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 1 ) );
        presets.add( new JLabel( "Presets:" ) );
        presets.add( button( "Single column (85 mm)", () -> applySizePreset( unit, width, 85 ) ) );
        presets.add( button( "Double column (170 mm)", () -> applySizePreset( unit, width, 170 ) ) );
        presets.setMaximumSize( new Dimension( Integer.MAX_VALUE, presets.getPreferredSize().height ) );
        final JComponent[] sub = { labeled( "Unit:", unit ), labeled( "Width:", width ), labeled( "Height:", height ),
                labeled( "Resolution (DPI, raster only):", dpi ), presets, summary };
        final JCheckBox fixed = new JCheckBox( "Export at a fixed figure size", o.isExportUseFixedSize() );
        fixed.setToolTipText( "<html>When on, every graphics export (PDF / SVG / EPS / PNG / JPG / TIFF) is rendered at "
                + "EXACTLY the size below,<br>laying the tree out to fill that frame -- so the file matches your target "
                + "figure size (e.g. a journal column width).<br>The whole tree is exported (the raster scale and "
                + "\"visible region only\" options above do not apply to a fixed-size export).</html>" );
        fixed.addActionListener( e -> {
            o.setExportUseFixedSize( fixed.isSelected() );
            setEnabledDeep( sub, fixed.isSelected() );
        } );
        add( c, fixed );
        for ( final JComponent s : sub ) {
            add( c, s );
        }
        setEnabledDeep( sub, o.isExportUseFixedSize() );
        refresh.run();
    }

    /** Clamps a fixed-export width/height into the spinner's [MIN,MAX] range. A {@link SpinnerNumberModel} THROWS
     *  (not clamps) on an out-of-range seed, and a unit conversion can overshoot the max, so both the seed and the
     *  conversion clamp through here -- otherwise reopening the dialog after such a conversion would fail. */
    private static double clampExportDim( final double v ) {
        return Math.max( AptxConstants.EXPORT_SIZE_DIM_MIN, Math.min( AptxConstants.EXPORT_SIZE_DIM_MAX, v ) );
    }

    /** A journal-column preset: sets the unit to millimetres (converting the current height) and the width to
     *  {@code mm}. Both {@code setSelectedItem}/{@code setValue} fire the real listeners, which write Options. */
    private void applySizePreset( final JComboBox<ExportSizeSpec.Unit> unit, final JSpinner width, final double mm ) {
        unit.setSelectedItem( ExportSizeSpec.Unit.MILLIMETERS );
        width.setValue( mm );
    }

    /** Enables/disables each component AND all its descendants (a disabled Swing container does not grey its
     *  children, so the labels/controls inside the {@code labeled(...)} rows must be toggled directly). */
    private static void setEnabledDeep( final JComponent[] comps, final boolean enabled ) {
        for ( final JComponent comp : comps ) {
            setEnabledDeep( comp, enabled );
        }
    }

    private static void setEnabledDeep( final Component comp, final boolean enabled ) {
        comp.setEnabled( enabled );
        if ( comp instanceof java.awt.Container ) {
            for ( final Component child : ( (java.awt.Container) comp ).getComponents() ) {
                setEnabledDeep( child, enabled );
            }
        }
    }

    private JPanel readTab() {
        final JPanel c = column();
        c.add( header( "Newick / NHX / Nexus Reading" ) );
        add( c, cb( _mf._internal_number_are_confidence_for_nh_parsing_cbmi ) );
        add( c, cb( _mf._replace_underscores_cbmi ) );
        add( c, cb( _mf._parse_beast_style_extended_nexus_tags_cbmi ) );
        add( c, cb( _mf._allow_errors_in_distance_to_parent_cbmi ) );
        return c;
    }

    private JPanel saveTab() {
        final JPanel c = column();
        c.add( header( "Newick / Nexus Saving" ) );
        c.add( new JLabel( "Write confidence values as:" ) );
        add( c, cb( _mf._use_brackets_for_conf_in_nh_export_cbmi ) );
        add( c, cb( _mf._use_internal_names_for_conf_in_nh_export_cbmi ) );
        return c;
    }

    private JPanel cacheTab() {
        final JPanel c = column();
        c.add( header( "Persistent Taxonomy Cache" ) );
        final JLabel intro = new JLabel( "<html>NCBI taxonomy lookups are remembered on disk for 30 days, so trees of"
                + " organisms you've already seen<br>load without re-querying. If the cache can't be written,"
                + " Archaeopteryx still works &mdash; lookups<br>just aren't remembered between sessions.</html>" );
        intro.setAlignmentX( Component.LEFT_ALIGNMENT );
        c.add( intro );
        c.add( Box.createVerticalStrut( 8 ) );
        _cache_enabled_cb = new JCheckBox( "Use persistent cache" );
        _cache_enabled_cb.addActionListener( e -> {
            NcbiTaxonomyLineageService.getShared().setPersistentCacheEnabled( _cache_enabled_cb.isSelected() );
            refreshCacheTab();
        } );
        add( c, _cache_enabled_cb );
        c.add( Box.createVerticalStrut( 4 ) );
        _cache_location_label = new JLabel();
        _cache_size_label = new JLabel();
        _cache_status_label = new JLabel();
        add( c, _cache_location_label );
        add( c, _cache_size_label );
        add( c, _cache_status_label );
        c.add( Box.createVerticalStrut( 8 ) );
        add( c, button( "Clear Cache", () -> {
            final int ok = JOptionPane.showConfirmDialog( this, "Delete all cached taxonomy data?",
                                                          "Clear Taxonomy Cache", JOptionPane.OK_CANCEL_OPTION,
                                                          JOptionPane.QUESTION_MESSAGE );
            if ( ok == JOptionPane.OK_OPTION ) {
                NcbiTaxonomyLineageService.getShared().clearPersistentCache();
                refreshCacheTab();
            }
        } ) );
        refreshCacheTab();
        return c;
    }

    /** Reads the current cache status from disk and updates the tab's labels/checkbox. */
    private void refreshCacheTab() {
        final TaxonomyCacheStatus s = NcbiTaxonomyLineageService.getShared().getCacheStatus();
        _cache_enabled_cb.setSelected( s.isEnabled() );
        _cache_location_label.setText( "Location: " + s.getPath() );
        if ( !s.isAvailable() ) {
            _cache_size_label.setText( " " );
            final String why = ( s.getUnavailableReason() == null ) ? "" : ( s.getUnavailableReason() + " " );
            _cache_status_label.setText( "Cache unavailable: " + why + "— lookups still work, just slower." );
        }
        else if ( !s.isEnabled() ) {
            // available, but switched off: make clear nothing is being read/written right now
            _cache_size_label.setText( "Disabled — " + TaxonomyCacheStatus.formatBytes( s.getBytes() ) + " ("
                    + s.getEntries() + " taxa) retained on disk." );
            _cache_status_label.setText( "Re-check \"Use persistent cache\" to use it again." );
        }
        else {
            _cache_size_label.setText( "Size: " + TaxonomyCacheStatus.formatBytes( s.getBytes() ) + " — "
                    + s.getEntries() + " taxa" );
            _cache_status_label.setText( TaxonomyCacheStatus.describeAge( s.getOldestEpochMs(),
                                                                         System.currentTimeMillis() ) );
        }
    }

    /** Updates the "Current font" line from the active tree's base font. */
    private void refreshFontInfo() {
        if ( _current_font_label == null ) {
            return;
        }
        try {
            final Font f = _mf.getMainPanel().getTreeFontSet().getBaseFont();
            _current_font_label.setText( "Current font:   " + f.getFamily() + ",   " + f.getSize() + " pt,   "
                    + styleName( f ) );
        }
        catch ( final Exception e ) {
            _current_font_label.setText( " " );
        }
    }

    private static String styleName( final Font f ) {
        if ( f.isBold() && f.isItalic() ) {
            return "Bold Italic";
        }
        if ( f.isBold() ) {
            return "Bold";
        }
        if ( f.isItalic() ) {
            return "Italic";
        }
        return "Plain";
    }

    /** A short HTML blurb describing the three bundled, always-available figure fonts. */
    private static String buildPreferredFontsHtml() {
        // single-line items: the dialog's width is set by its tab-header row (7 tabs), which comfortably
        // exceeds this blurb, so no width constraint is needed (and constraining only forces taller wrapping)
        final StringBuilder sb = new StringBuilder(
                "<html><i>Three publication-quality fonts are bundled and always available:</i>"
                        + "<ul style='margin-top:3px;margin-bottom:0'>" );
        for ( final FontResources.Preferred p : FontResources.PREFERRED ) {
            sb.append( "<li><b>" ).append( p.family ).append( "</b> &mdash; " ).append( p.description ).append( "</li>" );
        }
        return sb.append( "</ul></html>" ).toString();
    }

    // ---- bindings ------------------------------------------------------------------------------

    /** A checkbox bound to an existing menu item; toggling it does the menu item's full action. */
    private JCheckBox cb( final JCheckBoxMenuItem mi ) {
        if ( mi == null ) {
            return null;
        }
        final JCheckBox c = new JCheckBox( mi.getText(), mi.isSelected() );
        c.addActionListener( e -> {
            if ( mi.isSelected() != c.isSelected() ) {
                mi.doClick();
            }
        } );
        return c;
    }

    /** A checkbox bound to a PER-TREE (per-panel) boolean on the CURRENT TreePanel: seeded from it at construction,
     *  writes to it on toggle (a no-op when no tree is loaded). Used for the per-tab Time-Axis refinement toggles. */
    private JCheckBox panelCheckbox( final String label, final String tip,
                                     final java.util.function.Predicate<TreePanel> getter,
                                     final java.util.function.BiConsumer<TreePanel, Boolean> setter ) {
        final TreePanel cur = _mf.getCurrentTreePanel();
        final JCheckBox c = new JCheckBox( label, ( cur != null ) && getter.test( cur ) );
        c.setToolTipText( tip );
        c.addActionListener( e -> {
            if ( _refreshing ) {
                return; // re-seed on a tab switch -- reflect the tab, don't re-apply the toggle
            }
            final TreePanel p = _mf.getCurrentTreePanel();
            if ( p != null ) {
                setter.accept( p, c.isSelected() );
            }
        } );
        _tab_reseeders.add( () -> {
            final TreePanel p = _mf.getCurrentTreePanel();
            c.setSelected( ( p != null ) && getter.test( p ) );
        } );
        return c;
    }

    /** Lays out the given radio menu items as a mutually-exclusive group of bound radio buttons. */
    private void addRadioGroup( final JPanel col, final JRadioButtonMenuItem... items ) {
        final ButtonGroup group = new ButtonGroup();
        for ( final JRadioButtonMenuItem mi : items ) {
            if ( mi == null ) {
                continue;
            }
            final JRadioButton r = new JRadioButton( mi.getText(), mi.isSelected() );
            group.add( r );
            r.addActionListener( e -> {
                if ( r.isSelected() && !mi.isSelected() ) {
                    mi.doClick();
                }
            } );
            add( col, r );
        }
    }

    /** The four rectangular sub-styles, in dropdown order. They differ only in how a branch JOINT is drawn. */
    private static final Options.PHYLOGENY_GRAPHICS_TYPE[] RECTANGULAR_STYLES = {
            Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, Options.PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE,
            Options.PHYLOGENY_GRAPHICS_TYPE.ROUNDED, Options.PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR };

    /**
     * The RECTANGULAR sub-style dropdown: Square / Euro Type / Rounded / Triangular.
     * <p>
     * Circular and unrooted are deliberately NOT in this list -- which of the five primary display types you are
     * looking at is a control-panel decision now, since it is changed far too often to live behind a dialog.
     * <p>
     * The dropdown stays usable in EVERY layout, including the radial ones. Picking a style while a circular or
     * unrooted tree is on screen does not switch the layout (changing a joint style must not yank you out of the
     * view you are working in); it is remembered and applies the moment a rectangular layout is picked. Both
     * paths go through {@link ControlPanel#setRectangularStyle}, which owns that distinction.
     */
    private JComboBox<String> treeStyleCombo() {
        final String[] labels = new String[ RECTANGULAR_STYLES.length ];
        for ( int i = 0; i < RECTANGULAR_STYLES.length; ++i ) {
            final JCheckBoxMenuItem item = typeMenuItemFor( RECTANGULAR_STYLES[ i ] );
            labels[ i ] = ( item != null ) ? item.getText() : RECTANGULAR_STYLES[ i ].toString();
        }
        final JComboBox<String> combo = new JComboBox<>( labels );
        combo.setSelectedIndex( indexOfCurrentStyle() );
        combo.setToolTipText( RECTANGULAR_STYLE_TIP );
        combo.addActionListener( e -> {
            if ( _refreshing ) {
                return; // re-seed on a tab switch -- reflect the tab, don't re-apply the style
            }
            final int i = combo.getSelectedIndex();
            if ( ( i >= 0 ) && ( i < RECTANGULAR_STYLES.length ) ) {
                controlPanel().setRectangularStyle( RECTANGULAR_STYLES[ i ] );
            }
        } );
        // The style in force is per-tab in a rectangular layout and pending in a radial one; the control panel
        // holds both in one place, so re-seeding is just "show whatever it says". Runs on a tab switch AND after
        // any layout change (MainFrame.typeChanged -> refreshOpenSettingsDialog).
        _tab_reseeders.add( () -> combo.setSelectedIndex( indexOfCurrentStyle() ) );
        return combo;
    }

    /** The dropdown index of the rectangular sub-style currently in force (live or pending). */
    private int indexOfCurrentStyle() {
        final Options.PHYLOGENY_GRAPHICS_TYPE style = controlPanel().getRectangularStyle();
        for ( int i = 0; i < RECTANGULAR_STYLES.length; ++i ) {
            if ( RECTANGULAR_STYLES[ i ] == style ) {
                return i;
            }
        }
        return 0;
    }

    private JCheckBoxMenuItem typeMenuItemFor( final Options.PHYLOGENY_GRAPHICS_TYPE style ) {
        switch ( style ) {
            case EURO_STYLE:
                return _mf._euro_type_cbmi;
            case ROUNDED:
                return _mf._rounded_type_cbmi;
            case TRIANGULAR:
                return _mf._triangular_type_cbmi;
            default:
                return _mf._rectangular_type_cbmi;
        }
    }

    private ControlPanel controlPanel() {
        return _mf.getMainPanel().getControlPanel();
    }

    /** The tip-label angle is a vertical-orientation setting, so it is live only for a rectangular root-top/bottom. */
    private boolean tipLabelAngleApplies() {
        // the CURRENT TAB's orientation: the tip-label angle applies to the tree being looked at, not to the
        // default a future tab would open with
        final TreePanel tp = _mf.getCurrentTreePanel();
        final Options.TREE_ORIENTATION o = ( tp != null ) ? tp.getTreeOrientation()
                : _mf.getOptions().getTreeOrientation();
        return isVerticalOrientation( o ) && !isRadialStyleSelected();
    }

    private boolean isRadialStyleSelected() {
        return ( ( _mf._circular_type_cbmi != null ) && _mf._circular_type_cbmi.isSelected() )
                || ( ( _mf._unrooted_type_cbmi != null ) && _mf._unrooted_type_cbmi.isSelected() );
    }

    private static boolean isVerticalOrientation( final Options.TREE_ORIENTATION o ) {
        return ( o == Options.TREE_ORIENTATION.ROOT_TOP ) || ( o == Options.TREE_ORIENTATION.ROOT_BOTTOM );
    }

    /** Orientation swaps the layout's width and height, so re-fit the whole tree (recomputes the scaling constants
     *  and preferred size) before repainting -- a plain repaint would leave the old scroll extent. */
    private void reFitCurrentTree() {
        if ( _mf.getCurrentTreePanel() != null ) {
            _mf.getMainPanel().getControlPanel().updateZoomButtonsForLayout(); // relabel the zoom cluster for the layout
            _mf.getMainPanel().getControlPanel().showWhole();
            _mf.getCurrentTreePanel().repaint();
        }
    }

    /** When the geologic time axis is switched on for a tree that has NO absolute dates, prompt for the root age so
     *  the axis has an absolute calibration to map; then re-fit (the axis reserves a bottom band). */
    private void maybeCalibrateAndFit( final Options.TIME_AXIS_TYPE type ) {
        final TreePanel tp = _mf.getCurrentTreePanel();
        if ( ( type == Options.TIME_AXIS_TYPE.GEOLOGIC ) && ( tp != null ) && ( tp.timeAxisRootAgeMa() <= 0 ) ) {
            setTimeAxisRootAgeDialog();
        }
        else if ( ( type == Options.TIME_AXIS_TYPE.CALENDAR ) && ( tp != null ) && ( tp.timeAxisPresentDate() <= 0 ) ) {
            setTimeAxisPresentDateDialog();
        }
        reFitCurrentTree();
    }

    /** Asks for the most-recent-tip calendar date (decimal year) for the CALENDAR axis, for a tree WITHOUT calendar
     *  {@code <date>} elements; a tip-dated tree derives it automatically. */
    private void setTimeAxisPresentDateDialog() {
        final TreePanel tp = _mf.getCurrentTreePanel();
        if ( tp == null ) {
            return;
        }
        final double current = tp.timeAxisPresentDate();
        final String in = JOptionPane.showInputDialog( this,
                "Calendar date of the MOST-RECENT tip, as a decimal year (e.g. 2021.5 for mid-2021).\n"
                        + "This calibrates the calendar time axis for a tree WITHOUT calendar <date> elements.\n"
                        + "(A tip-dated tree whose dates are calendar years derives this automatically.)",
                ( current > 0 ) ? Double.toString( current ) : "" );
        if ( in == null ) {
            return; // cancelled
        }
        try {
            final double year = Double.parseDouble( in.trim() );
            if ( year > 0 ) {
                tp.setTimeAxisPresentDate( year );
                reFitCurrentTree();
            }
        }
        catch ( final NumberFormatException e ) {
            JOptionPane.showMessageDialog( this, "Please enter a positive decimal year (e.g. 2021.5).",
                                           "Invalid date", JOptionPane.WARNING_MESSAGE );
        }
    }

    /** Asks for the root age in Ma (for a tree without {@code <date>} elements) and applies it to the current tree. */
    private void setTimeAxisRootAgeDialog() {
        final TreePanel tp = _mf.getCurrentTreePanel();
        if ( tp == null ) {
            return;
        }
        final double current = tp.timeAxisRootAgeMa();
        final String in = JOptionPane.showInputDialog( this,
                "Age of the ROOT, in millions of years before present (Ma).\n"
                        + "This calibrates the geologic time axis for a tree WITHOUT <date> elements.\n"
                        + "(A dated tree derives this automatically.)",
                ( current > 0 ) ? Double.toString( current ) : "" );
        if ( in == null ) {
            return; // cancelled
        }
        try {
            final double ma = Double.parseDouble( in.trim() );
            if ( ma > 0 ) {
                tp.setTimeAxisRootAge( ma );
                reFitCurrentTree();
            }
        }
        catch ( final NumberFormatException e ) {
            JOptionPane.showMessageDialog( this, "Please enter a positive number (Ma).", "Invalid age",
                                           JOptionPane.WARNING_MESSAGE );
        }
    }

    // The categorical palette used by the "Color by:" feature for the current tree.
    private JComboBox<String> paletteCombo() {
        final JComboBox<String> combo = new JComboBox<>( PropertyColorScheme.paletteNames().toArray( new String[ 0 ] ) );
        final TreePanel tp = _mf.getCurrentTreePanel();
        if ( tp != null ) {
            combo.setSelectedItem( tp.getColorPaletteName() );
        }
        combo.addActionListener( e -> {
            if ( _refreshing ) {
                return; // re-seed on a tab switch -- reflect the tab, don't re-apply the palette
            }
            final TreePanel cur = _mf.getCurrentTreePanel();
            if ( ( cur != null ) && ( combo.getSelectedItem() != null ) ) {
                cur.setColorPaletteName( combo.getSelectedItem().toString() );
            }
        } );
        _tab_reseeders.add( () -> { // "Color by" palette is per-tab
            final TreePanel cur = _mf.getCurrentTreePanel();
            if ( cur != null ) {
                combo.setSelectedItem( cur.getColorPaletteName() );
            }
        } );
        return combo;
    }

    private interface Setter<T> {
        void set( T value );
    }

    private <T> JComboBox<T> enumCombo( final T[] values, final T current, final Setter<T> setter ) {
        return enumCombo( values, current, setter, null );
    }

    /**
     * A combo over {@code values}. When {@code labeller} is non-null it renders each item through it -- used for
     * enums whose own name is ALL-CAPS (NodeShape, NodeFill) so the user sees "Rectangle", not "RECTANGLE". The
     * model still holds the real enum constants, so selection and the setter are unaffected. Enums that already
     * define a friendly {@code toString()} (e.g. OVERVIEW_PLACEMENT_TYPE) pass null and render via that.
     */
    private <T> JComboBox<T> enumCombo( final T[] values, final T current, final Setter<T> setter,
                                        final Function<? super T, String> labeller ) {
        final JComboBox<T> combo = new JComboBox<>( values );
        if ( labeller != null ) {
            combo.setRenderer( new DefaultListCellRenderer() {
                @Override
                public Component getListCellRendererComponent( final JList<?> list, final Object value,
                                                               final int index, final boolean isSelected,
                                                               final boolean cellHasFocus ) {
                    super.getListCellRendererComponent( list, value, index, isSelected, cellHasFocus );
                    if ( value != null ) {
                        @SuppressWarnings( "unchecked" )
                        final T v = (T) value;
                        setText( labeller.apply( v ) );
                    }
                    return this;
                }
            } );
        }
        if ( current != null ) {
            combo.setSelectedItem( current );
        }
        combo.addActionListener( e -> {
            if ( _refreshing ) {
                return; // a programmatic re-seed on a tab switch -- do not write the value back onto the tab
            }
            @SuppressWarnings( "unchecked" )
            final T v = (T) combo.getSelectedItem();
            setter.set( v );
        } );
        return combo;
    }

    /** Display form of an ALL-CAPS enum constant: "RECTANGLE" -> "Rectangle", "LOWER_LEFT" -> "Lower Left". */
    static String prettyEnumName( final Enum<?> value ) {
        final StringBuilder sb = new StringBuilder();
        for ( final String word : value.name().toLowerCase( Locale.ROOT ).split( "_" ) ) {
            if ( word.isEmpty() ) {
                continue;
            }
            if ( sb.length() > 0 ) {
                sb.append( ' ' );
            }
            sb.append( Character.toUpperCase( word.charAt( 0 ) ) ).append( word.substring( 1 ) );
        }
        return sb.toString();
    }

    private JSpinner intSpinner( final int value, final int min, final int max, final int step,
                                 final Setter<Integer> setter ) {
        final JSpinner s = new JSpinner( new SpinnerNumberModel( value, min, max, step ) );
        s.addChangeListener( e -> setter.set( ( (Number) s.getValue() ).intValue() ) );
        return s;
    }

    private JSpinner doubleSpinner( final double value, final double min, final double max, final double step,
                                    final Setter<Double> setter ) {
        final JSpinner s = new JSpinner( new SpinnerNumberModel( value, min, max, step ) );
        s.addChangeListener( e -> setter.set( ( (Number) s.getValue() ).doubleValue() ) );
        return s;
    }

    private JButton button( final String label, final Runnable action ) {
        final JButton b = new JButton( label );
        b.addActionListener( e -> action.run() );
        return b;
    }

    /** "Reset to Defaults": confirm, then reset all display settings + theme + tree colors, and close the dialog
     *  (its controls were built from the pre-reset options, so leaving it open would show stale values -- reopen
     *  Settings to see the defaults). */
    private void confirmAndResetToDefaults() {
        final int choice = JOptionPane.showConfirmDialog( this,
                "<html>Reset all display settings to their defaults?<br><br>"
                        + "This also switches the theme to <b>Light</b>, resets the search options, and turns off "
                        + "property-based <b>&quot;Color by&quot;</b> (back to the default palette) on "
                        + "<b>all open trees</b>.<br>"
                        + "Manually applied branch/clade colors and your loaded trees are not changed.</html>",
                "Reset to Defaults", JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE );
        if ( choice == JOptionPane.OK_OPTION ) {
            _mf.resetToDefaults();
            setVisible( false );
            dispose();
        }
    }

    // ---- apply / layout helpers ----------------------------------------------------------------

    private void repaintTree() {
        if ( _mf.getCurrentTreePanel() != null ) {
            _mf.getCurrentTreePanel().repaint();
        }
    }

    /** Pushes the chosen "Found/Selected Colors" hue onto the live tree color set and repaints. */
    private void applyFoundColor( final Options.FOUND_COLOR v ) {
        if ( _mf.getMainPanel() != null ) {
            final TreeColorSet cs = _mf.getMainPanel().getTreeColorSet();
            if ( cs != null ) {
                cs.setFoundColor( v );
            }
        }
        repaintTree();
    }

    private void updateOverview() {
        if ( _mf.getCurrentTreePanel() != null ) {
            _mf.getCurrentTreePanel().updateOvSettings();
            _mf.getCurrentTreePanel().repaint();
        }
    }

    private static JScrollPane scroll( final JPanel content ) {
        final JScrollPane sp = new JScrollPane( content, ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
                                                ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER );
        sp.setBorder( null );
        sp.getVerticalScrollBar().setUnitIncrement( 16 );
        return sp;
    }

    private static JPanel column() {
        final JPanel p = new JPanel();
        p.setLayout( new BoxLayout( p, BoxLayout.Y_AXIS ) );
        p.setBorder( BorderFactory.createEmptyBorder( 10, 14, 14, 14 ) );
        return p;
    }

    private static JLabel header( final String text ) {
        final JLabel l = new JLabel( text );
        l.setFont( l.getFont().deriveFont( Font.BOLD ) );
        l.setAlignmentX( Component.LEFT_ALIGNMENT );
        l.setBorder( BorderFactory.createEmptyBorder( 10, 0, 4, 0 ) );
        return l;
    }

    private static JComponent labeled( final String label, final JComponent control ) {
        final JPanel p = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 1 ) );
        p.add( new JLabel( label ) );
        p.add( control );
        // Without this the row's max height is unbounded, so the Y-axis BoxLayout stretches it and the
        // FlowLayout centers the label/control inside, leaving a large empty gap below the row.
        p.setMaximumSize( new Dimension( Integer.MAX_VALUE, p.getPreferredSize().height ) );
        return p;
    }

    /** Adds a control to the column (skipping nulls), left-aligned, with a little spacing. */
    private static void add( final JPanel col, final JComponent c ) {
        if ( c == null ) {
            return;
        }
        c.setAlignmentX( Component.LEFT_ALIGNMENT );
        col.add( c );
        col.add( Box.createVerticalStrut( 2 ) );
    }
}
