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

    SettingsDialog( final MainFrame mf ) {
        super( mf, "Settings", false );
        _mf = mf;
        final JTabbedPane tabs = new JTabbedPane();
        tabs.addTab( "Display", scroll( displayTab() ) );
        tabs.addTab( "Nodes & Branches", scroll( nodesTab() ) );
        tabs.addTab( "Search", scroll( searchTab() ) );
        tabs.addTab( "Export & Print", scroll( exportTab() ) );
        tabs.addTab( "File Reading", scroll( readTab() ) );
        tabs.addTab( "File Saving", scroll( saveTab() ) );
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
        add( c, cb( _mf._show_overview_cbmi ) );
        add( c, labeled( "Overview placement:", enumCombo( OVERVIEW_PLACEMENT_TYPE.values(),
                                                           _mf.getOptions().getOvPlacement(),
                                                           v -> { _mf.getOptions().setOvPlacement( v );
                                                                  updateOverview(); } ) ) );
        c.add( header( "Labels & Colors" ) );
        add( c, cb( _mf._color_labels_same_as_parent_branch ) );
        add( c, cb( _mf._abbreviate_scientific_names ) );
        add( c, cb( _mf._show_confidence_stddev_cbmi ) );
        add( c, cb( _mf._show_mad_confidence_cbmi ) );
        add( c, cb( _mf._screen_antialias_cbmi ) );
        add( c, cb( _mf._background_gradient_cbmi ) );
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
        add( c, labeled( "Node fill:", enumCombo( NodeFill.values(), _mf.getOptions().getDefaultNodeFill(),
                                                  v -> { _mf.getOptions().setDefaultNodeFill( v ); repaintTree(); } ) ) );
        add( c, labeled( "Node size:", intSpinner( _mf.getOptions().getDefaultNodeShapeSize(), 0, 100, 1,
                                                   v -> { _mf.getOptions().setDefaultNodeShapeSize( v.shortValue() ); repaintTree(); } ) ) );
        c.add( header( "Branches & Confidence" ) );
        add( c, labeled( "Branch width:", doubleSpinner( _mf.getOptions().getDefaultBranchWidth(), 0.5, 20, 0.5,
                                                         v -> { _mf.getOptions().setDefaultBranchWidth( v.floatValue() ); repaintTree(); } ) ) );
        add( c, labeled( "Min. confidence shown:", doubleSpinner( _mf.getOptions().getMinConfidenceValue(), 0, 1.0E8, 0.1,
                                                                  v -> { _mf.getOptions().setMinConfidenceValue( v ); repaintTree(); } ) ) );
        add( c, labeled( "Support symbols:", enumCombo( SUPPORT_VISUALIZATION.values(),
                                                        _mf.getOptions().getSupportVisualization(),
                                                        v -> { _mf.getOptions().setSupportVisualization( v ); repaintTree(); } ) ) );
        add( c, labeled( "Support threshold (0–1):", doubleSpinner( _mf.getOptions().getSupportThreshold(), 0, 1.0, 0.05,
                                                                       v -> { _mf.getOptions().setSupportThreshold( v ); repaintTree(); } ) ) );
        c.add( header( "Colors & Font" ) );
        add( c, labeled( "Color scheme:", button( "Choose…", () -> _mf.switchColors() ) ) );
        add( c, labeled( "Tree font:", button( "Choose…", () -> _mf.chooseFont() ) ) );
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
        c.add( header( "Taxonomy Extraction from Node Names" ) );
        addRadioGroup( c, _mf._extract_taxonomy_no_rbmi, _mf._extract_taxonomy_pfam_strict_rbmi,
                       _mf._extract_taxonomy_pfam_relaxed_rbmi, _mf._extract_taxonomy_agressive_rbmi );
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
