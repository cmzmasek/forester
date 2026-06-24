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

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
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
import org.forester.phylogeny.data.NodeDataField;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.ws.seqdb.NcbiTaxonomyLineageService;
import org.forester.ws.seqdb.TaxonomyCacheStatus;

/**
 * A modeless, live-apply "Settings" dialog that replaces the old Options and Type menus (and the
 * View/control-panel Light-Dark toggle). Every control is visible at once -- no cycling. Checkboxes,
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

    SettingsDialog( final MainFrame mf ) {
        super( mf, "Settings", false );
        _mf = mf;
        final JTabbedPane tabs = new JTabbedPane();
        tabs.addTab( "Display", scroll( displayTab() ) );
        tabs.addTab( "Fonts, Nodes and Branches", scroll( nodesTab() ) );
        tabs.addTab( "Search", scroll( searchTab() ) );
        tabs.addTab( "Export & Print", scroll( exportTab() ) );
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
        final JButton close = new JButton( "Close" );
        close.addActionListener( e -> setVisible( false ) );
        final JPanel south = new JPanel( new FlowLayout( FlowLayout.RIGHT ) );
        south.add( close );
        setLayout( new BorderLayout() );
        add( tabs, BorderLayout.CENTER );
        add( south, BorderLayout.SOUTH );
        pack();
        setLocationRelativeTo( mf );
    }

    // ---- tabs ----------------------------------------------------------------------------------

    private JPanel displayTab() {
        final JPanel c = column();
        c.add( header( "Theme" ) );
        addThemeControls( c );
        c.add( header( "Layout" ) );
        add( c, labeled( "Tree style:", treeStyleCombo() ) );
        addRadioGroup( c, _mf._ext_node_dependent_cladogram_rbmi, _mf._non_lined_up_cladograms_rbmi );
        add( c, cb( _mf._label_direction_cbmi ) );
        add( c, cb( _mf._show_scale_cbmi ) );
        add( c, cb( _mf._show_scale_grid_cbmi ) );
        add( c, cb( _mf._show_overview_cbmi ) );
        add( c, labeled( "Overview placement:", enumCombo( OVERVIEW_PLACEMENT_TYPE.values(),
                                                           _mf.getOptions().getOvPlacement(),
                                                           v -> { _mf.getOptions().setOvPlacement( v );
                                                                  updateOverview(); } ) ) );
        c.add( header( "Labels & Colors" ) );
        c.add( new JLabel( "Internal label placement:" ) );
        addRadioGroup( c, _mf._internal_labels_above_branch_rbmi, _mf._internal_labels_right_of_node_rbmi );
        add( c, cb( _mf._color_labels_same_as_parent_branch ) );
        add( c, cb( _mf._use_italic_scientific_names_cbmi ) );
        add( c, cb( _mf._abbreviate_scientific_names ) );
        add( c, cb( _mf._show_confidence_stddev_cbmi ) );
        add( c, cb( _mf._show_mad_confidence_cbmi ) );
        add( c, labeled( "\"Color by\" palette:", paletteCombo() ) );
        c.add( header( "Collapsed Subtrees & Domains" ) );
        add( c, cb( _mf._collapsed_with_average_height_cbmi ) );
        add( c, cb( _mf._show_abbreviated_labels_for_collapsed_nodes_cbmi ) );
        add( c, cb( _mf._line_up_renderable_data_cbmi ) );
        add( c, cb( _mf._right_line_up_domains_cbmi ) );
        add( c, cb( _mf._show_domain_labels ) );
        return c;
    }

    private JPanel nodesTab() {
        final JPanel c = column();
        c.add( header( "Node Shapes" ) );
        add( c, cb( _mf._show_default_node_shapes_external_cbmi ) );
        add( c, cb( _mf._show_default_node_shapes_internal_cbmi ) );
        add( c, cb( _mf._show_default_node_shapes_for_marked_cbmi ) );
        add( c, labeled( "Node shape:", enumCombo( NodeShape.values(), _mf.getOptions().getDefaultNodeShape(),
                                                   v -> { _mf.getOptions().setDefaultNodeShape( v ); repaintTree(); } ) ) );
        // GRADIENT node fill is retired as a user choice (kept in the NodeFill enum only for phyloXML round-tripping).
        add( c, labeled( "Node fill:", enumCombo( new NodeFill[] { NodeFill.DEFAULT, NodeFill.NONE, NodeFill.SOLID },
                                                  _mf.getOptions().getDefaultNodeFill(),
                                                  v -> { _mf.getOptions().setDefaultNodeFill( v ); repaintTree(); } ) ) );
        add( c, labeled( "Node size:", intSpinner( _mf.getOptions().getDefaultNodeShapeSize(), 0, 100, 1,
                                                   v -> { _mf.getOptions().setDefaultNodeShapeSize( v.shortValue() ); repaintTree(); } ) ) );
        c.add( header( "Branches & Confidence" ) );
        add( c, labeled( "Branch width:", doubleSpinner( _mf.getOptions().getDefaultBranchWidth(), 0.5, 20, 0.5,
                                                         v -> { _mf.getOptions().setDefaultBranchWidth( v.floatValue() ); repaintTree(); } ) ) );
        add( c, labeled( "Min. confidence shown:", doubleSpinner( _mf.getOptions().getMinConfidenceValue(), 0, 1.0E8, 0.1,
                                                                  v -> { _mf.getOptions().setMinConfidenceValue( v ); repaintTree(); } ) ) );
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
        c.add( header( "Behavior" ) );
        add( c, labeled( "Data returned on copy:", enumCombo( NodeDataField.values(),
                                                              _mf.getOptions().getExtDescNodeDataToReturn(),
                                                              v -> _mf.getOptions().setExtDescNodeDataToReturn( v ) ) ) );
        return c;
    }

    private JPanel searchTab() {
        final JPanel c = column();
        c.add( header( "Search & Selection" ) );
        add( c, cb( _mf._color_all_found_nodes_when_coloring_subtree_cbmi ) );
        c.add( Box.createVerticalStrut( 4 ) );
        final JLabel note = new JLabel( "<html><i>The per-search options (Match Case, Words, Regex, Inverse,"
                + " Properties, Visible) are on the left control panel, next to the search boxes.</i></html>" );
        note.setAlignmentX( Component.LEFT_ALIGNMENT );
        c.add( note );
        return c;
    }

    private JPanel exportTab() {
        final JPanel c = column();
        c.add( header( "Graphics Export & Printing" ) );
        final JSpinner raster_scale = intSpinner( _mf.getOptions().getRasterExportScale(), 1, 8, 1,
                                                  v -> _mf.getOptions().setRasterExportScale( v ) );
        raster_scale.setToolTipText( "<html>Resolution multiplier for raster (PNG/JPG/TIFF) export: the figure is "
                + "re-rendered onto an N&times;-larger canvas for crisp, print-quality output<br>"
                + "(a true re-render, not pixel doubling). 1 = on-screen size. Higher = larger files/memory; "
                + "very large figures are capped automatically.<br>Does not affect vector (SVG/EPS/PDF) export.</html>" );
        add( c, labeled( "Raster export scale (×):", raster_scale ) );
        add( c, cb( _mf._transparent_export_background_cbmi ) );
        add( c, cb( _mf._outline_fonts_in_vector_export_cbmi ) );
        add( c, cb( _mf._antialias_print_cbmi ) );
        add( c, cb( _mf._print_black_and_white_cbmi ) );
        add( c, cb( _mf._graphics_export_visible_only_cbmi ) );
        add( c, labeled( "PDF line width:", doubleSpinner( _mf.getOptions().getPrintLineWidth(), 0.5, 20, 0.5,
                                                           v -> _mf.getOptions().setPrintLineWidth( v.floatValue() ) ) ) );
        return c;
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

    private void addThemeControls( final JPanel col ) {
        final boolean dark = _mf.getConfiguration().getUi() == Configuration.UI.FLAT_DARK;
        final ButtonGroup group = new ButtonGroup();
        final JRadioButton light = new JRadioButton( "Light", !dark );
        final JRadioButton dark_rb = new JRadioButton( "Dark", dark );
        group.add( light );
        group.add( dark_rb );
        light.addActionListener( e -> { if ( light.isSelected() ) { _mf.setDarkMode( false ); } } );
        dark_rb.addActionListener( e -> { if ( dark_rb.isSelected() ) { _mf.setDarkMode( true ); } } );
        final JPanel row = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 0 ) );
        row.add( light );
        row.add( dark_rb );
        add( col, row );
    }

    /** The tree-style dropdown (former "Type" menu); selecting an entry clicks its menu item. */
    private JComboBox<String> treeStyleCombo() {
        final JCheckBoxMenuItem[] items = { _mf._rectangular_type_cbmi, _mf._euro_type_cbmi, _mf._rounded_type_cbmi,
                _mf._curved_type_cbmi, _mf._triangular_type_cbmi, _mf._convex_type_cbmi, _mf._circular_type_cbmi,
                _mf._unrooted_type_cbmi };
        final String[] labels = new String[ items.length ];
        int selected = 0;
        for ( int i = 0; i < items.length; ++i ) {
            labels[ i ] = ( items[ i ] != null ) ? items[ i ].getText() : "?";
            if ( ( items[ i ] != null ) && items[ i ].isSelected() ) {
                selected = i;
            }
        }
        final JComboBox<String> combo = new JComboBox<>( labels );
        combo.setSelectedIndex( selected );
        combo.addActionListener( e -> {
            final int i = combo.getSelectedIndex();
            if ( ( i >= 0 ) && ( items[ i ] != null ) && !items[ i ].isSelected() ) {
                items[ i ].doClick();
            }
        } );
        return combo;
    }

    // The categorical palette used by the "Color by:" feature for the current tree.
    private JComboBox<String> paletteCombo() {
        final JComboBox<String> combo = new JComboBox<>( PropertyColorScheme.paletteNames().toArray( new String[ 0 ] ) );
        final TreePanel tp = _mf.getCurrentTreePanel();
        if ( tp != null ) {
            combo.setSelectedItem( tp.getColorPaletteName() );
        }
        combo.addActionListener( e -> {
            final TreePanel cur = _mf.getCurrentTreePanel();
            if ( ( cur != null ) && ( combo.getSelectedItem() != null ) ) {
                cur.setColorPaletteName( combo.getSelectedItem().toString() );
            }
        } );
        return combo;
    }

    private interface Setter<T> {
        void set( T value );
    }

    private <T> JComboBox<T> enumCombo( final T[] values, final T current, final Setter<T> setter ) {
        final JComboBox<T> combo = new JComboBox<>( values );
        if ( current != null ) {
            combo.setSelectedItem( current );
        }
        combo.addActionListener( e -> {
            @SuppressWarnings( "unchecked" )
            final T v = (T) combo.getSelectedItem();
            setter.set( v );
        } );
        return combo;
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

    // ---- apply / layout helpers ----------------------------------------------------------------

    private void repaintTree() {
        if ( _mf.getCurrentTreePanel() != null ) {
            _mf.getCurrentTreePanel().repaint();
        }
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
