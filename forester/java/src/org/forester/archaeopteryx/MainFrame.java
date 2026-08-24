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

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Container;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.SortedSet;
import java.util.NoSuchElementException;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JList;
import javax.swing.BorderFactory;
import javax.swing.JDialog;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import javax.swing.JMenu;
import javax.swing.KeyStroke;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.filechooser.FileFilter;

import com.formdev.flatlaf.FlatDarkLaf;
import com.formdev.flatlaf.FlatLightLaf;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.analysis.AncestralTaxonomyInference;
import org.forester.archaeopteryx.tools.AncestralTaxonomyInferrer;
import org.forester.archaeopteryx.tools.LabelDataExtractor;
import org.forester.archaeopteryx.tools.TipDateExtractor;
import org.forester.archaeopteryx.tools.TipDateExtractor.DayMonthOrder;
import org.forester.archaeopteryx.tools.TipDateExtractor.Summary;
import org.forester.archaeopteryx.tools.TipDateExtractor.TipDate;
import org.forester.archaeopteryx.tools.NodeDataExporter;
import org.forester.archaeopteryx.tools.NodeDataImporter;
import org.forester.archaeopteryx.tools.ProcessPool;
import org.forester.archaeopteryx.tools.ProcessRunning;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonomicLineageService;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.NodeDataField;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sdi.GSDI;
import org.forester.sdi.GSDIR;
import org.forester.sdi.SDIException;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.WindowsUtils;

public abstract class MainFrame extends JFrame implements ActionListener {

    /**
     * Installs the FlatLaf look-and-feel (light or dark) -- the only look-and-feel Archaeopteryx uses. It is
     * pure Java, so it renders identically on macOS, Windows and Linux.
     */
    static void installLookAndFeel(final Configuration.UI ui) {
        try {
            switch (ui) {
                case FLAT_DARK:
                    FlatDarkLaf.setup();
                    break;
                case FLAT_LIGHT:
                default:
                    FlatLightLaf.setup();
                    break;
            }
        }
        catch (final Exception e) {
            ForesterUtil.printWarningMessage(AptxConstants.PRG_NAME,
                    "could not set look and feel [" + ui + "]: " + e.getMessage());
            FlatLightLaf.setup();
        }
    }

    final static NHFilter nhfilter = new NHFilter();
    final static NHXFilter nhxfilter = new NHXFilter();
    final static XMLFilter xmlfilter = new XMLFilter();
    final static JsonFilter jsonfilter = new JsonFilter();
    final static TolFilter tolfilter = new TolFilter();
    final static NexusFilter nexusfilter = new NexusFilter();
    final static PdfFilter pdffilter = new PdfFilter();
    final static GraphicsFileFilter graphicsfilefilter = new GraphicsFileFilter();
    final static DefaultFilter defaultfilter = new DefaultFilter();
    static final String USE_MOUSEWHEEL_SHIFT_TO_ROTATE = "rotate with mousewheel + Shift (or A and S), D toggles between horizontal and radial labels";
    private static final long serialVersionUID = 3655000897845508358L;
    final static Font menu_font = new Font(Configuration.getDefaultFontFamilyName(),
            Font.PLAIN,
            Configuration.getGuiFontSize());
    static final String TYPE_MENU_HEADER = "Type";
    static final String RECTANGULAR_TYPE_CBMI_LABEL = "Rectangular";
    static final String EURO_TYPE_CBMI_LABEL = "Euro Type";
    static final String TRIANGULAR_TYPE_CBMI_LABEL = "Triangular";
    static final String ROUNDED_TYPE_CBMI_LABEL = "Rounded";
    static final String UNROOTED_TYPE_CBMI_LABEL = "Unrooted (alpha)";                                                                                                                                                          //TODO
    static final String CIRCULAR_TYPE_CBMI_LABEL = "Circular (alpha)";                                                                                                                                                          //TODO
    static final String SEARCH_CASE_SENSITIVE_LABEL = "Match Case";
    static final String INVERSE_SEARCH_RESULT_LABEL = "Inverse";
    static final String DISPLAY_SCALE_LABEL = "Scale";
    static final String SHOW_TREE_NAME_LABEL = "Tree Name";
    static final String DISPLAY_SCALE_GRID_LABEL = "Scale Grid Lines";
    static final String DISPLAY_SCALE_GRID_TIP = "Draw faint vertical reference lines at scale intervals across the tree (phylograms only), so branch depths are easy to compare.";
    static final String GEOLOGIC_GRID_LABEL = "Time-Axis Grid Lines";
    static final String GEOLOGIC_GRID_TIP = "When a Time Axis is on, draw faint reference lines across the tree at its tick positions (the geologic band boundaries, or the calendar year ticks).";
    static final String GEOLOGIC_AGES_LABEL = "Geologic Boundary Ages";
    static final String GEOLOGIC_AGES_TIP = "When the Geologic time axis is on, label the coarser band's interval boundaries with their age (e.g. \"201.4 Ma\" at the base of the Jurassic).";
    static final String DISPLAY_SCALE_AXIS_LABEL = "Scale Axis";
    static final String DISPLAY_SCALE_AXIS_TIP = "Draw a labeled distance axis with tick marks along the bottom (phylograms only), so branch lengths can be read off directly.";
    static final String DISPLAY_HPD_BARS_LABEL = "Node Age Bars (HPD)";
    static final String DISPLAY_HPD_BARS_TIP = "On a dated (time-scaled) phylogram, draw a bar at each internal node spanning its age uncertainty (the node's phyloXML date min/max).";
    static final String DISPLAY_FOSSIL_RANGE_BARS_LABEL = "Fossil Range Bars (FAD/LAD)";
    static final String DISPLAY_FOSSIL_RANGE_BARS_TIP = "On a dated (time-scaled) phylogram, draw a stratigraphic-range bar at each fossil tip spanning its first-to-last appearance (the tip's phyloXML date min/max).";
    static final String DISPLAY_ZEBRA_STRIPES_LABEL = "Zebra Stripes";
    static final String DISPLAY_ZEBRA_STRIPES_TIP = "Shade every other tip row with a faint band, so a label is easy to track across a wide tree to its annotation columns.";
    static final String DISPLAY_BREAK_LONG_BRANCHES_LABEL = "Break Long Branches";
    static final String DISPLAY_BREAK_LONG_BRANCHES_TIP = "On a phylogram, draw an outlier-long branch (e.g. a distant outgroup) shortened with a break mark and give the rest of the tree the freed width. The true branch length is unchanged (still shown as its label).";
    static final String DISPLAY_INTERNAL_TAXONOMY_KEY_LABEL = "Internal Taxonomy Key";
    static final String DISPLAY_INTERNAL_TAXONOMY_KEY_TIP = "Show a draggable key of the distinct internal-node taxa (from inference / curation / clade annotation), grouped by rank with counts.";
    static final String DISPLAY_REVERSE_TIP_ORDER_LABEL = "Reverse Tip Order";
    static final String DISPLAY_REVERSE_TIP_ORDER_TIP = "Reverse the order of the tips (mirror the tree across the tip axis). In a root-top/bottom orientation the tips run sideways, so this flips them left-to-right. Rectangular layouts only.";
    static final String DISPLAY_TIP_LABELS_BELOW_COLUMNS_LABEL = "Tip Labels Below Columns";
    static final String DISPLAY_TIP_LABELS_BELOW_COLUMNS_TIP = "Clustergram layout: in a root-top/bottom orientation with annotation columns, draw the tip labels BELOW the columns (so the dendrogram sits directly on the tip-aligned grid) instead of between the tree and the columns.";
    static final String DISPLAY_BOLD_FOUND_LABELS_LABEL = "Bold Found Labels";
    static final String DISPLAY_BOLD_FOUND_LABELS_TIP = "Render the labels of found/selected nodes in bold, so search hits stand out. Works on screen and in exports.";
    static final String DISPLAY_DIM_NON_MATCHES_LABEL = "Dim Non-Matches";
    static final String DISPLAY_DIM_NON_MATCHES_TIP = "While a search or selection is active, fade non-matching labels toward the background so the hits stand out. Works on screen and in exports.";
    static final String DISPLAY_PULSE_FOUND_NODES_LABEL = "Pulse Found Nodes";
    static final String DISPLAY_PULSE_FOUND_NODES_TIP = "Draw a gently pulsing halo around found/selected nodes to draw the eye (a static glow in exports; not in black-and-white; rectangular layouts only).";
    static final String NON_LINED_UP_CLADOGRAMS_LABEL = "Non-Lined Up Cladogram";
    static final String LABEL_DIRECTION_LABEL = "Radial Labels";
    static final String LABEL_DIRECTION_TIP = "To use radial node labels in radial and unrooted display types";
    static final String COLOR_LABELS_LABEL = "Colorize Labels Same as Parent Branch";
    static final String DISPLAY_NODE_BOXES_LABEL_EXT = "Shapes for External Nodes";
    static final String DISPLAY_NODE_BOXES_LABEL_INT = "Shapes for Internal Nodes";
    static final String DISPLAY_NODE_BOXES_LABEL_MARKED = "Shapes for Nodes with Visual Data";
    static final String INTERNAL_LABELS_ABOVE_BRANCH_LABEL = "Above Branch (Left of Node)";
    static final String INTERNAL_LABELS_RIGHT_OF_NODE_LABEL = "Right of Node";
    static final String INTERNAL_LABELS_ABOVE_BRANCH_TIP = "Place internal-node labels to the left of the node, right-aligned, on top of the branch";
    static final String INTERNAL_LABELS_RIGHT_OF_NODE_TIP = "Place internal-node labels to the right of the node (classic placement)";
    static final String SHOW_OVERVIEW_LABEL = "Overview";
    static final String NONUNIFORM_CLADOGRAMS_LABEL = "Lined Up Cladogram";
    static final String COLOR_LABELS_TIP = "To use parent branch colors for node labels as well, need to turn off taxonomy dependent colorization and turn on branch colorization for this to become apparent";
    static final String ABBREV_SN_LABEL = "Abbreviate Scientific Taxonomic Names";
    static final String ITALIC_SN_LABEL = "Italic Scientific Taxonomic Names";
    static final String ITALIC_SN_TIP = "To set scientific (Latin) taxonomic names in italics (e.g. Homo sapiens), per the publication convention; only the scientific-name part is italicized, not the code, common name or rank";
    static final String OUTLINE_FONTS_VECTOR_LABEL = "Outline Fonts (SVG/EPS export)";
    static final String OUTLINE_FONTS_VECTOR_TIP = "SVG/EPS export embeds no fonts; outlining renders text as vector shapes so viewers cannot substitute the figure font (recommended). Turn off to keep selectable/searchable text, at the risk of font substitution.";
    static final String TRANSPARENT_BG_LABEL = "Transparent Background (PNG)";
    static final String TRANSPARENT_BG_TIP = "Export PNG with a transparent (alpha) background instead of the solid figure background. Ignored for formats that cannot carry transparency (JPG, etc.).";
    static final String WHITE_BG_LABEL = "White Image Background";
    static final String WHITE_BG_TIP = "Export with the LIGHT theme -- white background and dark, legible labels -- regardless of the on-screen light/dark theme, so a figure is document-ready (no dark box, and not light-on-white). Applies to raster images (PNG/JPG/TIFF), 'Copy Image to Clipboard', and vector (SVG/EPS/PDF). Turn off to export exactly what is on screen (WYSIWYG). Transparent-PNG is unaffected.";
    static final String SHOW_CONF_STDDEV_LABEL = "Confidence Standard Deviations";
    static final String SHOW_MAD_CONF_LABEL    = "MAD Confidence Values (MAD/regular)";
    static final String USE_BRACKETS_FOR_CONF_IN_NH_LABEL = "Use Brackets for Confidence Values";
    static final String USE_INTERNAL_NAMES_FOR_CONF_IN_NH_LABEL = "Use Internal Node Names for Confidence Values";
    static final String SHOW_BASIC_TREE_INFORMATION_LABEL = "Basic Tree Information";
    static final String EDIT_TREE_INFO_LABEL = "Edit Tree Name and Description...";
    static final String INFER_ANCESTOR_TAXONOMIES = "Infer Ancestor Taxonomies";
    static final String OBTAIN_SEQUENCE_AND_TAXONOMIC_INFORMATION = "Fetch Sequence & Taxonomic Data";
    JMenuBar _jmenubar;
    FoundSelectedCounter _found_selected_counter; // right-aligned menu-bar counter of highlighted (found+selected) nodes
    JMenu _file_jmenu;
    JMenu _tools_menu;
    JMenu _view_jmenu;
    JMenu _settings_jmenu;
    JMenu _help_jmenu;
    // Analysis menu
    JMenu _analysis_menu;
    JMenuItem _load_species_tree_item;
    JMenuItem _gsdi_item;
    JMenuItem _gsdir_item;
    JMenuItem _create_tanglegram_item;
    JMenuItem _lineage_inference;
    // file menu:
    JMenuItem _open_item;
    JMenuItem _save_item;
    JMenuItem _save_all_item;
    JMenuItem _close_item;
    JMenuItem _exit_item;
    JMenuItem _new_item;
    JMenuItem _write_to_pdf_item;
    JMenuItem _write_to_jpg_item;
    JMenuItem _write_to_tif_item;
    JMenuItem _write_to_png_item;
    JMenuItem _copy_image_to_clipboard_item;
    JMenuItem _write_to_svg_item;
    JMenuItem _write_to_eps_item;
    JMenuItem _export_seqs_fasta_item;
    JMenuItem _export_node_data_item;
    JMenuItem _import_annotations_item;
    JMenuItem _import_annotations_url_item;
    JMenuItem _import_gtdb_item;
    JMenuItem _reimport_annotations_item;
    // tools menu:
    JMenuItem _midpoint_root_item;
    JMenuItem _mad_root_item;
    JMenuItem _color_rank_jmi;
    JMenuItem _node_style_selected_jmi;
    JMenuItem _clade_bands_jmi;
    JMenuItem _annotation_columns_jmi;
    JMenu _edit_jmenu;
    JMenuItem _undo_item;
    JMenuItem _redo_item;

    JMenuItem _obtain_seq_and_tax_information_jmi;
    JMenuItem _extract_label_data_jmi;
    JMenuItem _extract_dates_jmi;
    JMenuItem _remove_branch_color_item;
    JMenuItem _remove_visual_styles_item;
    JMenuItem _delete_selected_nodes_item;
    JMenuItem _delete_not_selected_nodes_item;
    // font size menu:
    // options menu:
    // _  screen and print
    JCheckBoxMenuItem _label_direction_cbmi;
    // _  screen display
    JRadioButtonMenuItem _non_lined_up_cladograms_rbmi;
    JRadioButtonMenuItem _ext_node_dependent_cladogram_rbmi;
    JCheckBoxMenuItem _show_scale_cbmi;                                                                                                                                                                                                      //TODO fix me
    JCheckBoxMenuItem _show_tree_name_cbmi;
    JCheckBoxMenuItem _show_scale_grid_cbmi;
    JCheckBoxMenuItem _show_scale_axis_cbmi;
    JCheckBoxMenuItem _show_hpd_bars_cbmi;
    JCheckBoxMenuItem _show_fossil_range_bars_cbmi;
    JCheckBoxMenuItem _show_zebra_stripes_cbmi;
    JCheckBoxMenuItem _break_long_branches_cbmi;
    JCheckBoxMenuItem _show_internal_taxonomy_key_cbmi;
    JCheckBoxMenuItem _tip_labels_below_columns_cbmi;
    JCheckBoxMenuItem _reverse_tip_order_cbmi;
    JCheckBoxMenuItem _bold_found_labels_cbmi;
    JCheckBoxMenuItem _dim_non_matches_cbmi;
    JCheckBoxMenuItem _pulse_found_nodes_cbmi;
    JCheckBoxMenuItem _show_overview_cbmi;
    JCheckBoxMenuItem _abbreviate_scientific_names;
    JCheckBoxMenuItem _use_italic_scientific_names_cbmi;
    JCheckBoxMenuItem _outline_fonts_in_vector_export_cbmi;
    JCheckBoxMenuItem _transparent_export_background_cbmi;
    JCheckBoxMenuItem _graphics_export_white_background_cbmi;
    JCheckBoxMenuItem _color_labels_same_as_parent_branch;
    JCheckBoxMenuItem _show_default_node_shapes_internal_cbmi;
    JRadioButtonMenuItem _internal_labels_above_branch_rbmi;
    JRadioButtonMenuItem _internal_labels_right_of_node_rbmi;
    JCheckBoxMenuItem _show_default_node_shapes_external_cbmi;
    JCheckBoxMenuItem _show_default_node_shapes_for_marked_cbmi;
    JCheckBoxMenuItem _show_confidence_stddev_cbmi;
    JCheckBoxMenuItem _show_mad_confidence_cbmi;
    JCheckBoxMenuItem _collapsed_with_average_height_cbmi;
    JCheckBoxMenuItem _show_abbreviated_labels_for_collapsed_nodes_cbmi;
    // _  print
    JCheckBoxMenuItem _graphics_export_visible_only_cbmi;
    JCheckBoxMenuItem _antialias_print_cbmi;
    JCheckBoxMenuItem _print_black_and_white_cbmi;
    //JMenuItem                        _print_size_mi;
    // _  parsing
    JCheckBoxMenuItem _internal_number_are_confidence_for_nh_parsing_cbmi;
    JCheckBoxMenuItem _replace_underscores_cbmi;
    JCheckBoxMenuItem _allow_errors_in_distance_to_parent_cbmi;
    JCheckBoxMenuItem _use_brackets_for_conf_in_nh_export_cbmi;
    JCheckBoxMenuItem _use_internal_names_for_conf_in_nh_export_cbmi;
    JCheckBoxMenuItem _parse_beast_style_extended_nexus_tags_cbmi;
    // type menu:
    JMenu _type_menu;
    JCheckBoxMenuItem _rectangular_type_cbmi;
    JCheckBoxMenuItem _triangular_type_cbmi;
    JCheckBoxMenuItem _euro_type_cbmi;
    JCheckBoxMenuItem _rounded_type_cbmi;
    JCheckBoxMenuItem _unrooted_type_cbmi;
    JCheckBoxMenuItem _circular_type_cbmi;
    // view as text menu:
    JMenuItem _view_as_NH_item;
    JMenuItem _view_as_XML_item;
    JMenuItem _view_as_nexus_item;
    JMenuItem _display_basic_information_item;
    JMenuItem _edit_tree_info_item;
    JMenuItem _fit_to_window_item;
    JMenuItem _clustergram_item;
    JMenuItem _find_next_hit_item;
    JMenuItem _find_prev_hit_item;
    // help menu:
    JMenuItem _about_item;
    JMenuItem _keyboard_shortcuts_item;
    JMenuItem _help_item;
    JMenuItem _website_item;
    JMenuItem _aptxjs_website_item;
    JMenuItem _references_item;

    //
    // Last-used file-dialog directory, kept separately per purpose (read / save / export) so each dialog reopens
    // where it was last left. Persisted across restarts by DirectoryPreferences; a purpose with no chosen (or a
    // no-longer-readable) directory falls back to the home/Desktop directory via getCurrentDir(...).
    final EnumMap<DirectoryPreferences.Category, File> _current_dirs = new EnumMap<>(DirectoryPreferences.Category.class);
    // A transient launch-time default (e.g. the working directory when opened from the command line with a tree):
    // it seeds where dialogs open before the user has chosen anything, but is NEVER persisted -- only directories
    // the user actually navigated to (stored in _current_dirs) are remembered across restarts.
    File _startup_dir;
    JFileChooser _writetopdf_filechooser;
    JFileChooser _save_filechooser;
    JFileChooser _writetographics_filechooser;
    // process menu:
    JMenu _process_menu;
    MainPanel _mainpanel;
    Container _contentpane;
    final LinkedList<TextFrame> _textframes = new LinkedList<>();
    Configuration _configuration;
    Options _options;
    private Phylogeny _species_tree;
    // the rank last chosen in "Annotate Clades by Rank", pre-selected next time (per session); null = first use
    private String _last_clade_rank;
    final ProcessPool _process_pool;

    MainFrame() {
        _process_pool = ProcessPool.createInstance();
        _writetopdf_filechooser = new JFileChooser();
        _writetopdf_filechooser.setMultiSelectionEnabled(false);
        _writetopdf_filechooser.addChoosableFileFilter(pdffilter);
        _writetographics_filechooser = new JFileChooser();
        _writetographics_filechooser.setMultiSelectionEnabled(false);
        _writetographics_filechooser.addChoosableFileFilter(graphicsfilefilter);
        _save_filechooser = new JFileChooser();
        _save_filechooser.setMultiSelectionEnabled(false);
        _save_filechooser.setFileFilter(xmlfilter);
        _save_filechooser.addChoosableFileFilter(nhfilter);
        _save_filechooser.addChoosableFileFilter(nexusfilter);
        _save_filechooser.addChoosableFileFilter(_save_filechooser.getAcceptAllFileFilter());
        try {
            final String home_dir = System.getProperty("user.home");
            _save_filechooser.setCurrentDirectory(new File(home_dir));
            _writetopdf_filechooser.setCurrentDirectory(new File(home_dir));
            _writetographics_filechooser.setCurrentDirectory(new File(home_dir));
        } catch (final Exception e) {
            e.printStackTrace();
            // Do nothing. Not important.
        }
    }

    /**
     * Action performed.
     */
    @Override
    public void actionPerformed(final ActionEvent e) {
        final Object o = e.getSource();

        if (o == _exit_item) {
            close();
        } else if (o == _gsdi_item) {
            if (isSubtreeDisplayed()) {
                return;
            }
            executeGSDI();
        } else if (o == _gsdir_item) {
            if (isSubtreeDisplayed()) {
                return;
            }
            executeGSDIR();
        } else if (o == _color_rank_jmi) {
            colorRank();
        } else if (o == _node_style_selected_jmi) {
            nodeStyleForSelectedNodes();
        } else if (o == _clade_bands_jmi) {
            labelCladesByRank();
        } else if (o == _annotation_columns_jmi) {
            chooseAnnotationColumns();
        } else if (o == _undo_item) {
            undo();
        } else if (o == _redo_item) {
            redo();
        } else if (o == _remove_branch_color_item) {
            if (isSubtreeDisplayed()) {
                return;
            }
            removeBranchColors();
        } else if (o == _remove_visual_styles_item) {
            if (isSubtreeDisplayed()) {
                return;
            }
            removeVisualStyles();
        } else if (o == _midpoint_root_item) {
            if (isSubtreeDisplayed()) {
                return;
            }
            midpointRoot();
        } else if (o == _mad_root_item) {
            if (isSubtreeDisplayed()) {
                return;
            }
            madRoot();
        } else if (o == _delete_selected_nodes_item) {
            if (isSubtreeDisplayed()) {
                return;
            }
            deleteSelectedNodes(true);
        } else if (o == _delete_not_selected_nodes_item) {
            if (isSubtreeDisplayed()) {
                return;
            }
            deleteSelectedNodes(false);
        } else if (o == _edit_tree_info_item) {
            showTreeInfoDialog();
        } else if (o == _display_basic_information_item) {
            if (getCurrentTreePanel() != null) {
                displayBasicInformation(getCurrentTreePanel().getTreeFile());
            }
        } else if (o == _view_as_NH_item) {
            viewAsNH();
        } else if (o == _view_as_XML_item) {
            viewAsXML();
        } else if (o == _view_as_nexus_item) {
            viewAsNexus();
        } else if (o == _fit_to_window_item) {
            showWhole();
        } else if (o == _clustergram_item) {
            applyClustergramPreset();
        } else if (o == _find_next_hit_item) {
            if (getCurrentTreePanel() != null) {
                getCurrentTreePanel().stepToFoundNode(1);
            }
        } else if (o == _find_prev_hit_item) {
            if (getCurrentTreePanel() != null) {
                getCurrentTreePanel().stepToFoundNode(-1);
            }
        } else if (o == _abbreviate_scientific_names) {
            updateOptions(getOptions());
        } else if (o == _use_italic_scientific_names_cbmi) {
            updateOptions(getOptions());
        } else if (o == _outline_fonts_in_vector_export_cbmi) {
            updateOptions(getOptions());
        } else if (o == _transparent_export_background_cbmi) {
            updateOptions(getOptions());
        } else if (o == _graphics_export_white_background_cbmi) {
            updateOptions(getOptions());
        } else if (o == _color_labels_same_as_parent_branch) {
            updateOptions(getOptions());
        } else if (o == _show_default_node_shapes_internal_cbmi) {
            updateOptions(getOptions());
        } else if (o == _internal_labels_above_branch_rbmi) {
            updateOptions(getOptions());
        } else if (o == _internal_labels_right_of_node_rbmi) {
            updateOptions(getOptions());
        } else if (o == _show_default_node_shapes_external_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_default_node_shapes_for_marked_cbmi) {
            updateOptions(getOptions());
        } else if (o == _non_lined_up_cladograms_rbmi) {
            updateOptions(getOptions());
            showWhole();
        } else if (o == _ext_node_dependent_cladogram_rbmi) {
            updateOptions(getOptions());
            showWhole();
        } else if (o == _parse_beast_style_extended_nexus_tags_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_scale_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_tree_name_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_scale_grid_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_scale_axis_cbmi) {
            updateOptions(getOptions());
            // the scale axis reserves a tip-spread band (a side ruler in a vertical orientation, a bottom band in a
            // horizontal one); toggling it must re-fit so the reserve is applied/reclaimed and the axis clears the
            // tips -- but ONLY where the axis is actually drawn (a rectangular-family phylogram with a scale). For a
            // cladogram or the radial circular/unrooted layouts the toggle is inert, so leave the view (and any radial
            // zoom) untouched. showWhole re-fits to the viewport in BOTH orientations -- no preferred-size feedback,
            // so no depth-zoom drift (a plain fitWidth/fitHeight feeds the extent back and creeps the zoom by MOVE).
            final TreePanel tp = getCurrentTreePanel();
            if ((tp != null) && tp.scaleAxisAppliesToLayout()) {
                showWhole();
            }
        } else if (o == _show_hpd_bars_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_fossil_range_bars_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_zebra_stripes_cbmi) {
            updateOptions(getOptions());
        } else if (o == _break_long_branches_cbmi) {
            updateOptions(getOptions());
            // capping an outlier branch changes the depth scale / radial spread (the informative part reclaims the freed
            // width or disc), so re-fit whenever the toggle FLIPS -- ON to apply the capping, OFF to restore the uncapped
            // layout (the guard is option-INDEPENDENT: it asks whether this layout is one the feature affects, so OFF
            // re-fits too). Fires for a rectangular-family phylogram (unaligned "P" or aligned "A") AND for a radial
            // (circular / unrooted) phylogram -- the latter re-fits the radial zoom via showWhole -- but NOT for a
            // cladogram. showWhole re-fits to the viewport -> no preferred-size feedback / drift.
            final TreePanel tp = getCurrentTreePanel();
            if ((tp != null)
                    && (tp.breakLongBranchesRelevantToLayout() || tp.breakLongBranchesRelevantToRadialLayout())) {
                showWhole();
            }
        } else if (o == _show_internal_taxonomy_key_cbmi) {
            updateOptions(getOptions());
        } else if (o == _tip_labels_below_columns_cbmi) {
            updateOptions(getOptions());
        } else if (o == _reverse_tip_order_cbmi) {
            updateOptions(getOptions());
        } else if (o == _bold_found_labels_cbmi) {
            updateOptions(getOptions());
        } else if (o == _dim_non_matches_cbmi) {
            updateOptions(getOptions());
        } else if (o == _pulse_found_nodes_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_confidence_stddev_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_mad_confidence_cbmi) {
            updateOptions(getOptions());
        } else if (o == _use_brackets_for_conf_in_nh_export_cbmi) {
            if (_use_brackets_for_conf_in_nh_export_cbmi.isSelected()) {
                _use_internal_names_for_conf_in_nh_export_cbmi.setSelected(false);
            }
            updateOptions(getOptions());
        } else if (o == _use_internal_names_for_conf_in_nh_export_cbmi) {
            if (_use_internal_names_for_conf_in_nh_export_cbmi.isSelected()) {
                _use_brackets_for_conf_in_nh_export_cbmi.setSelected(false);
            }
            updateOptions(getOptions());
        } else if (o == _label_direction_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_overview_cbmi) {
            updateOptions(getOptions());
            if (getCurrentTreePanel() != null) {
                getCurrentTreePanel().updateOvSizes();
            }
        } else if (o == _collapsed_with_average_height_cbmi) {
            if (_collapsed_with_average_height_cbmi.isSelected()) {
                _collapsed_with_average_height_cbmi.setSelected(true);
            }
            updateOptions(getOptions());
        } else if (o == _show_abbreviated_labels_for_collapsed_nodes_cbmi) {
            if (_show_abbreviated_labels_for_collapsed_nodes_cbmi.isSelected()) {
                _show_abbreviated_labels_for_collapsed_nodes_cbmi.setSelected(true);
            }
            updateOptions(getOptions());
        } else if ((o == _rectangular_type_cbmi) || (o == _triangular_type_cbmi)
                || (o == _euro_type_cbmi) || (o == _rounded_type_cbmi)
                || (o == _unrooted_type_cbmi) || (o == _circular_type_cbmi)) {
            typeChanged(o);
        } else if (o == _about_item) {
            about();
        } else if (o == _keyboard_shortcuts_item) {
            KeyboardShortcuts.show(this);
        } else if (o == _references_item) {
            showReferences();
        } else if (o == _help_item) {
            try {
                AptxUtil.openWebsite(AptxConstants.APTX_DOC_SITE);
            } catch (final IOException e1) {
                ForesterUtil.printErrorMessage(AptxConstants.PRG_NAME, e1.toString());
            }
        } else if (o == _website_item) {
            try {
                AptxUtil.openWebsite(AptxConstants.APTX_WEB_SITE);
            } catch (final IOException e1) {
                ForesterUtil.printErrorMessage(AptxConstants.PRG_NAME, e1.toString());
            }
        } else if (o == _aptxjs_website_item) {
            try {
                AptxUtil.openWebsite(AptxConstants.APTX_JS_WEB_SITE);
            } catch (final IOException e1) {
                ForesterUtil.printErrorMessage(AptxConstants.PRG_NAME, e1.toString());
            }
        } else if (o == _write_to_pdf_item) {
            final File curr_dir = writeToPdf(_mainpanel.getCurrentPhylogeny(),
                    getMainPanel(),
                    _writetopdf_filechooser,
                    getCurrentDir(DirectoryPreferences.Category.EXPORT),
                    getContentPane(),
                    this);
            if (curr_dir != null) {
                setCurrentDir(DirectoryPreferences.Category.EXPORT, curr_dir);
            }
        } else if (o == _save_all_item) {
            writeAllToFile();
        } else if (o == _write_to_jpg_item) {
            final File new_dir = writeToGraphicsFile(_mainpanel.getCurrentPhylogeny(),
                    GraphicsExportType.JPG,
                    _mainpanel,
                    _writetographics_filechooser,
                    this,
                    getContentPane(),
                    getCurrentDir(DirectoryPreferences.Category.EXPORT));
            if (new_dir != null) {
                setCurrentDir(DirectoryPreferences.Category.EXPORT, new_dir);
            }
        } else if (o == _write_to_tif_item) {
            final File new_dir = writeToGraphicsFile(_mainpanel.getCurrentPhylogeny(),
                    GraphicsExportType.TIFF,
                    _mainpanel,
                    _writetographics_filechooser,
                    this,
                    getContentPane(),
                    getCurrentDir(DirectoryPreferences.Category.EXPORT));
            if (new_dir != null) {
                setCurrentDir(DirectoryPreferences.Category.EXPORT, new_dir);
            }
        } else if (o == _write_to_png_item) {
            final File new_dir = writeToGraphicsFile(_mainpanel.getCurrentPhylogeny(),
                    GraphicsExportType.PNG,
                    _mainpanel,
                    _writetographics_filechooser,
                    this,
                    getContentPane(),
                    getCurrentDir(DirectoryPreferences.Category.EXPORT));
            if (new_dir != null) {
                setCurrentDir(DirectoryPreferences.Category.EXPORT, new_dir);
            }
        } else if (o == _copy_image_to_clipboard_item) {
            copyImageToClipboard();
        } else if (o == _write_to_svg_item) {
            final File new_dir = writeToGraphicsFile(_mainpanel.getCurrentPhylogeny(),
                    GraphicsExportType.SVG,
                    _mainpanel,
                    _writetographics_filechooser,
                    this,
                    getContentPane(),
                    getCurrentDir(DirectoryPreferences.Category.EXPORT));
            if (new_dir != null) {
                setCurrentDir(DirectoryPreferences.Category.EXPORT, new_dir);
            }
        } else if (o == _write_to_eps_item) {
            final File new_dir = writeToGraphicsFile(_mainpanel.getCurrentPhylogeny(),
                    GraphicsExportType.EPS,
                    _mainpanel,
                    _writetographics_filechooser,
                    this,
                    getContentPane(),
                    getCurrentDir(DirectoryPreferences.Category.EXPORT));
            if (new_dir != null) {
                setCurrentDir(DirectoryPreferences.Category.EXPORT, new_dir);
            }
        } else if (o == _export_seqs_fasta_item) {
            exportSequencesAsFasta();
        } else if (o == _export_node_data_item) {
            exportNodeDataAsTsv();
        } else if (o == _import_annotations_item) {
            importAnnotations();
        } else if (o == _import_gtdb_item) {
            importGtdbTaxonomy();
        } else if (o == _import_annotations_url_item) {
            importAnnotationsFromUrl();
        } else if (o == _reimport_annotations_item) {
            reimportAnnotations();
        }  else if (o == _save_item) {
            final File new_dir = writeToFile(_mainpanel.getCurrentPhylogeny(),
                    getMainPanel(),
                    _save_filechooser,
                    getCurrentDir(DirectoryPreferences.Category.SAVE),
                    getContentPane(),
                    this);
            if (new_dir != null) {
                setCurrentDir(DirectoryPreferences.Category.SAVE, new_dir);
            }
        } else if (o == _graphics_export_visible_only_cbmi) {
            updateOptions(getOptions());
        } else if (o == _antialias_print_cbmi) {
            updateOptions(getOptions());
        } else if (o == _print_black_and_white_cbmi) {
            updateOptions(getOptions());
        } else if (o == _lineage_inference) {
            if (isSubtreeDisplayed()) {
                JOptionPane.showMessageDialog(this,
                        "Subtree is shown.",
                        "Cannot infer ancestral taxonomies",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }
            executeLineageInference();
        }
        _contentpane.repaint();
    }

    public Configuration getConfiguration() {
        return _configuration;
    }

    /**
     * This method returns the current external node data which
     * has been selected by the user by clicking the "Return ..."
     * menu item. This method is expected to be called from Javascript or
     * something like it.
     *
     * @return current external node data as String
     */

    public MainPanel getMainPanel() {
        return _mainpanel;
    }

    public Options getOptions() {
        return _options;
    }

    public ProcessPool getProcessPool() {
        return _process_pool;
    }

    public void showTextFrame(final String s, final String title) {
        checkTextFrames();
        _textframes.addLast(TextFrame.instantiate(s, title, _textframes));
    }

    public void showWhole() {
        _mainpanel.getControlPanel().showWhole();
    }

    public void updateProcessMenu() {
        // In general Swing is not thread safe.
        // See "Swing's Threading Policy".
        SwingUtilities.invokeLater(this::doUpdateProcessMenu);
    }

    void chooseFont() {
        final FontChooser fc = new FontChooser();
        fc.setFont(getMainPanel().getTreeFontSet().getBaseFont()); // the user size, not the transient auto-shrunk one
        fc.showDialog(this, "Select the Base Font");
        getMainPanel().getTreeFontSet().setBaseFont(fc.getFont());
        getControlPanel().displayedPhylogenyMightHaveChanged(true);
        if (getMainPanel().getCurrentTreePanel() != null) {
            getMainPanel().getCurrentTreePanel().resetPreferredSize();
            getMainPanel().getCurrentTreePanel().updateOvSizes();
        }

        repaint();
    }

    /**
     * Applies a few figure-display preferences programmatically to the currently displayed tree. The MSA
     * compactor opens a tree window from another package and wants a compact font and a visible scale bar for
     * the compacted-alignment tree; it can no longer push these through {@link Configuration} (which no longer
     * seeds {@link Options}), so it calls this right after {@link Archaeopteryx#createApplication}.
     *
     * @param base_font_family the tree font family, or {@code null} to leave the current font unchanged
     * @param base_font_size   the tree font size in points; ignored when {@code <= 0}
     * @param show_scale       whether to show the scale bar
     */
    public void applyTreeDisplayPreferences(final String base_font_family, final int base_font_size,
                                            final boolean show_scale) {
        getOptions().setShowScale(show_scale);
        if ((base_font_family != null) && (base_font_size > 0)) {
            final Font base_font = new Font(base_font_family, Font.PLAIN, base_font_size);
            getOptions().setBaseFont(base_font);
            getMainPanel().getTreeFontSet().setBaseFont(base_font);
        }
        // re-derive the label-width-driven layout at the new font (mirrors chooseFont)
        getControlPanel().displayedPhylogenyMightHaveChanged(true);
        final TreePanel tp = getMainPanel().getCurrentTreePanel();
        if (tp != null) {
            tp.resetPreferredSize();
            tp.updateOvSizes();
        }
        repaint();
    }

    private void deleteSelectedNodes(final boolean delete) {
        final String function = delete ? "Delete" : "Retain";
        final Phylogeny phy = getMainPanel().getCurrentPhylogeny();
        final int ext = (phy == null) ? 0 : phy.getNumberOfExternalNodes();
        // the external tips the user has selected (search a/b + manual clicks; branch-click a clade to select its
        // tips) -- shared with the export/protect scope so all three treat a selection identically (no internal
        // expansion). selectedExternalTips tolerates a null phy / empty selection, returning no tips.
        final TreePanel tp = getCurrentTreePanel();
        final List<PhylogenyNode> nodes = (tp == null) ? new ArrayList<>()
                : NodeDataExporter.selectedExternalTips(phy, tp.getFoundNodesAsListOfPhylogenyNodes());
        switch (AptxUtil.nodePruningOutcome(ext, nodes.size(), delete)) {
            case NO_TREE:
                JOptionPane.showMessageDialog(this,
                        "Load a tree with at least two external nodes before using \"" + function
                                + " Selected Nodes\".",
                        "Cannot " + function.toLowerCase() + " nodes",
                        JOptionPane.ERROR_MESSAGE);
                return;
            case NO_SELECTION:
                JOptionPane.showMessageDialog(this,
                        "Select one or more external nodes first — click them in the tree, or use the \"Search\" "
                                + "field to find and highlight them — then choose \"" + function
                                + " Selected Nodes\" again.",
                        "No external nodes selected to " + function.toLowerCase(),
                        JOptionPane.ERROR_MESSAGE);
                return;
            case WOULD_REMOVE_ALL:
                JOptionPane.showMessageDialog(this,
                        "That would remove every external node, leaving an empty tree.",
                        "Cannot " + function.toLowerCase() + " all nodes",
                        JOptionPane.ERROR_MESSAGE);
                return;
            default:
                break;
        }
        final int todo = nodes.size();
        final int res = delete ? (ext - todo) : todo;
        final int result = JOptionPane.showConfirmDialog(null, function + " " + todo
                + " external node(s), from a total of " + ext + " external nodes," + "\nresulting in tree with " + res
                + " nodes?", function + " external nodes", JOptionPane.OK_CANCEL_OPTION);
        if (result == JOptionPane.OK_OPTION) {
            getCurrentTreePanel().pushUndoCheckpoint(delete ? "Delete Nodes" : "Retain Nodes");
            if (!delete) {
                final List<PhylogenyNode> to_delete = new ArrayList<>();
                for (final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if (!nodes.contains(n)) {
                        to_delete.add(n);
                    }
                }
                for (final PhylogenyNode n : to_delete) {
                    phy.deleteSubtree(n, true);
                }
            } else {
                for (final PhylogenyNode n : nodes) {
                    phy.deleteSubtree(n, true);
                }
            }
            resetSearch();
            getCurrentTreePanel().setNodeInPreorderToNull();
            phy.externalNodesHaveChanged();
            phy.clearHashIdToNodeMap();
            phy.recalculateNumberOfExternalDescendants(true);
            getCurrentTreePanel().resetNodeIdToDistToLeafMap();
            // pruning changes the visible tips, so recompute the display schemes that summarize them
            // (this path does not go through displayedPhylogenyMightHaveChanged)
            getCurrentTreePanel().rebuildPropertyDisplays();
            getCurrentTreePanel().rebuildAnnotationColumns();
            getCurrentTreePanel().setEdited(true);
            repaint();
        }
    }

    private void doUpdateProcessMenu() {
        if (_process_pool.size() > 0) {
            if (_process_menu == null) {
                _process_menu = createMenu("", getConfiguration());
                _process_menu.setForeground(Color.RED);
            }
            _process_menu.removeAll();
            final String text = "processes running: " + _process_pool.size();
            _process_menu.setText(text);
            _jmenubar.add(_process_menu);
            for (int i = 0; i < _process_pool.size(); ++i) {
                final ProcessRunning p = _process_pool.getProcessByIndex(i);
                final boolean cancellable = (p.getProcess() != null) && !p.getProcess().isCancelled();
                final JMenuItem item = customizeJMenuItem(new JMenuItem(
                        p.getName() + " [" + p.getStart() + "]" + (cancellable ? "  — click to cancel" : "")));
                if (cancellable) {
                    item.addActionListener(e -> p.getProcess().requestCancel());
                }
                _process_menu.add(item);
            }
        } else {
            if (_process_menu != null) {
                _process_menu.removeAll();
                _jmenubar.remove(_process_menu);
            }
        }
        _jmenubar.validate();
        _jmenubar.repaint();
        repaint();
    }

    private void removeBranchColors() {
        if (getMainPanel().getCurrentPhylogeny() != null) {
            getMainPanel().getCurrentTreePanel().pushUndoCheckpoint("Delete Branch Colors");
            AptxUtil.removeBranchColors(getMainPanel().getCurrentPhylogeny());
            if (getMainPanel().getCurrentTreePanel() != null) {
                getMainPanel().getCurrentTreePanel().clearRankLegend(); // the rank legend no longer reflects the tree
            }
        }
    }

    private void removeVisualStyles() {
        if (getMainPanel().getCurrentPhylogeny() != null) {
            getMainPanel().getCurrentTreePanel().pushUndoCheckpoint("Delete Visual Styles");
            AptxUtil.removeVisualStyles(getMainPanel().getCurrentPhylogeny());
        }
    }

    private void writeAllToFile() {
        if ((getMainPanel().getTabbedPane() == null) || (getMainPanel().getTabbedPane().getTabCount() < 1)) {
            return;
        }
        final File my_dir = getCurrentDir(DirectoryPreferences.Category.SAVE);
        if (my_dir != null) {
            _save_filechooser.setCurrentDirectory(my_dir);
        }
        _save_filechooser.setSelectedFile(new File(""));
        final int result = _save_filechooser.showSaveDialog(_contentpane);
        final File file = _save_filechooser.getSelectedFile();
        setCurrentDir(DirectoryPreferences.Category.SAVE, _save_filechooser.getCurrentDirectory());
        if ((file != null) && (result == JFileChooser.APPROVE_OPTION)) {
            if (file.exists()) {
                final int i = JOptionPane.showConfirmDialog(this,
                        file + " already exists. Overwrite?",
                        "Warning",
                        JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.WARNING_MESSAGE);
                if (i != JOptionPane.OK_OPTION) {
                    return;
                } else {
                    try {
                        file.delete();
                    } catch (final Exception e) {
                        JOptionPane.showMessageDialog(this,
                                "Failed to delete: " + file,
                                "Error",
                                JOptionPane.WARNING_MESSAGE);
                    }
                }
            }
            final int count = getMainPanel().getTabbedPane().getTabCount();
            final List<Phylogeny> trees = new ArrayList<>();
            for (int i = 0; i < count; ++i) {
                final Phylogeny phy = getMainPanel().getPhylogeny(i);
                if (ForesterUtil.isEmpty(phy.getName())
                        && !ForesterUtil.isEmpty(getMainPanel().getTabbedPane().getTitleAt(i))) {
                    phy.setName(getMainPanel().getTabbedPane().getTitleAt(i));
                }
                trees.add(phy);
                getMainPanel().getTreePanels().get(i).syncTimeAxisConfigToTree(); // embed each tab's Time-Axis config
                getMainPanel().getTreePanels().get(i).setEdited(false);
            }
            final PhylogenyWriter writer = new PhylogenyWriter();
            try {
                writer.toPhyloXML(file, trees, 0, ForesterUtil.LINE_SEPARATOR);
            } catch (final IOException e) {
                JOptionPane.showMessageDialog(this,
                        "Failed to write to: " + file,
                        "Error",
                        JOptionPane.WARNING_MESSAGE);
            }
        }
    }

    void activateSaveAllIfNeeded() {
        final boolean multiple_trees = (getMainPanel().getTabbedPane() != null)
                && (getMainPanel().getTabbedPane().getTabCount() > 1);
        _save_all_item.setEnabled(multiple_trees);
        // a tanglegram compares two trees, so the item is only enabled once at least two are loaded
        if (_create_tanglegram_item != null) {
            _create_tanglegram_item.setEnabled(multiple_trees);
        }
    }

    void buildHelpMenu() {
        _help_jmenu = createMenu("Help", getConfiguration());
        _help_jmenu.setToolTipText("Documentation, web links, and program information");
        _help_jmenu.add(_help_item = new JMenuItem("Documentation"));
        _help_jmenu.addSeparator();
        _help_jmenu.add(_website_item = new JMenuItem("Archaeopteryx Home"));
        _help_jmenu.add(_aptxjs_website_item = new JMenuItem("Archaeopteryx online version: Archaeopteryx.js"));
        _help_jmenu.add(_references_item = new JMenuItem("References"));
        _help_jmenu.addSeparator();
        _help_jmenu.add(_keyboard_shortcuts_item = new JMenuItem("Keyboard Shortcuts"));
        _help_jmenu.addSeparator();
        _help_jmenu.add(_about_item = new JMenuItem("About"));
        customizeJMenuItem(_help_item);
        customizeJMenuItem(_website_item);
        customizeJMenuItem(_aptxjs_website_item);
        customizeJMenuItem(_references_item);
        customizeJMenuItem(_keyboard_shortcuts_item);
        customizeJMenuItem(_about_item);
        _keyboard_shortcuts_item.setToolTipText("List the keyboard shortcuts for navigating and viewing trees");
        _references_item.setToolTipText("Main literature citations for the algorithms implemented in Archaeopteryx");
        _jmenubar.add(_help_jmenu);
    }

    void buildTypeMenu() {
        _type_menu = createMenu(TYPE_MENU_HEADER, getConfiguration());
        _type_menu.add(_rectangular_type_cbmi = new JCheckBoxMenuItem(MainFrame.RECTANGULAR_TYPE_CBMI_LABEL));
        _type_menu.add(_euro_type_cbmi = new JCheckBoxMenuItem(MainFrame.EURO_TYPE_CBMI_LABEL));
        _type_menu.add(_rounded_type_cbmi = new JCheckBoxMenuItem(MainFrame.ROUNDED_TYPE_CBMI_LABEL));
        _type_menu.add(_triangular_type_cbmi = new JCheckBoxMenuItem(MainFrame.TRIANGULAR_TYPE_CBMI_LABEL));
        _type_menu.add(_unrooted_type_cbmi = new JCheckBoxMenuItem(MainFrame.UNROOTED_TYPE_CBMI_LABEL));
        _type_menu.add(_circular_type_cbmi = new JCheckBoxMenuItem(MainFrame.CIRCULAR_TYPE_CBMI_LABEL));
        customizeCheckBoxMenuItem(_rectangular_type_cbmi, false);
        customizeCheckBoxMenuItem(_triangular_type_cbmi, false);
        customizeCheckBoxMenuItem(_euro_type_cbmi, false);
        customizeCheckBoxMenuItem(_rounded_type_cbmi, false);
        customizeCheckBoxMenuItem(_unrooted_type_cbmi, false);
        customizeCheckBoxMenuItem(_circular_type_cbmi, false);
        _triangular_type_cbmi.setToolTipText(
                "Draws clades as triangles; reads cleanly only with aligned tips, so selecting it switches this tab to "
                        + "Cladogram (C). Switch back to P for a near-clock tree.");
        _unrooted_type_cbmi.setToolTipText(MainFrame.USE_MOUSEWHEEL_SHIFT_TO_ROTATE);
        _circular_type_cbmi.setToolTipText(MainFrame.USE_MOUSEWHEEL_SHIFT_TO_ROTATE);
        initializeTypeMenu(getOptions());
        // _type_menu is not added to the menu bar; its items are folded into the Settings dialog.
    }

    void buildViewMenu() {
        _view_jmenu = createMenu("View", getConfiguration());
        _view_jmenu.setToolTipText("Show tree information, or the tree as phyloXML, Newick, or Nexus");
        _view_jmenu.add(_edit_tree_info_item = new JMenuItem(EDIT_TREE_INFO_LABEL));
        _view_jmenu.add(_display_basic_information_item = new JMenuItem(SHOW_BASIC_TREE_INFORMATION_LABEL));
        _view_jmenu.addSeparator();
        _view_jmenu.add(_view_as_XML_item = new JMenuItem("as phyloXML"));
        _view_jmenu.add(_view_as_NH_item = new JMenuItem("as Newick"));
        _view_jmenu.add(_view_as_nexus_item = new JMenuItem("as Nexus"));
        _view_jmenu.addSeparator();
        _view_jmenu.add(_fit_to_window_item = new JMenuItem("Fit to Window"));
        _fit_to_window_item.setToolTipText("Zoom the tree to fit the window (also HOME / ESC)");
        _fit_to_window_item.setAccelerator(
                KeyStroke.getKeyStroke(KeyEvent.VK_0, Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx()));
        _view_jmenu.add(_clustergram_item = new JMenuItem("Clustergram"));
        _clustergram_item.setToolTipText("<html>One click: turn the tree into a vertical dendrogram over a tip-aligned "
                + "heat map — root at top, tips aligned, sample labels along the bottom, and each numeric per-tip "
                + "field as a shared-scale heat-map column (categorical fields as color strips).<br><i>The figure iTOL "
                + "does clunkily and FigTree/PearTree can't do at all. Import Annotations (CSV/TSV) first if the tree "
                + "has no per-tip data yet.</i></html>");
        _view_jmenu.addSeparator();
        _view_jmenu.add(_find_next_hit_item = new JMenuItem("Find Next"));
        _find_next_hit_item.setToolTipText("Center the next search hit in the view");
        _find_next_hit_item.setAccelerator(
                KeyStroke.getKeyStroke(KeyEvent.VK_G, Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx()));
        _view_jmenu.add(_find_prev_hit_item = new JMenuItem("Find Previous"));
        _find_prev_hit_item.setToolTipText("Center the previous search hit in the view");
        _find_prev_hit_item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_G,
                Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx() | InputEvent.SHIFT_DOWN_MASK));
        customizeJMenuItem(_edit_tree_info_item);
        customizeJMenuItem(_display_basic_information_item);
        customizeJMenuItem(_view_as_NH_item);
        customizeJMenuItem(_view_as_XML_item);
        customizeJMenuItem(_view_as_nexus_item);
        customizeJMenuItem(_fit_to_window_item);
        customizeJMenuItem(_clustergram_item);
        customizeJMenuItem(_find_next_hit_item);
        customizeJMenuItem(_find_prev_hit_item);
        _jmenubar.add(_view_jmenu);
    }

    void checkTextFrames() {
        if (_textframes.size() > 5) {
            try {
                if (_textframes.getFirst() != null) {
                    _textframes.getFirst().removeMe();
                } else {
                    _textframes.removeFirst();
                }
            } catch (final NoSuchElementException e) {
                // Ignore.
            }
        }
    }

    void close() {
        removeAllTextFrames();
        if (_mainpanel != null) {
            _mainpanel.terminate();
        }
        if (_contentpane != null) {
            _contentpane.removeAll();
        }
        setVisible(false);
        dispose();
    }

    /** Tools -> "Node Style for Selected Nodes…": opens the {@link NodeStyleDialog} over the current found/selected
     *  set (search hits A/B + manual selection). Warns when nothing is selected. */
    void nodeStyleForSelectedNodes() {
        final TreePanel tp = _mainpanel.getCurrentTreePanel();
        if (tp == null) {
            return;
        }
        final Phylogeny phy = tp.getPhylogeny();
        if ((phy == null) || phy.isEmpty()) {
            return;
        }
        final List<PhylogenyNode> selected = tp.getFoundNodesAsListOfPhylogenyNodes();
        if ((selected == null) || selected.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "No nodes are selected. Click nodes in the tree, or use the Search field(s), then choose this again.",
                    "No Selection", JOptionPane.WARNING_MESSAGE);
            return;
        }
        new NodeStyleDialog(tp, selected, false).setVisible(true);
    }

    void colorRank() {
        if (_mainpanel.getCurrentTreePanel() == null) {
            return;
        }
        final TreePanel tp = _mainpanel.getCurrentTreePanel();
        final Phylogeny phy = tp.getPhylogeny();
        if ((phy == null) || phy.isEmpty() || (phy.getNumberOfExternalNodes() < 2)) {
            return;
        }
        final Map<String, Integer> present_ranks = AptxUtil.getRankCounts(phy);
        final Map<String, Integer> coverage_counts = AptxUtil.getRankCoverageCounts(phy);
        final String[] ranks = AptxUtil.getRankChoices(present_ranks, coverage_counts, phy.getNumberOfExternalNodes());
        String rank = (String) JOptionPane.showInputDialog(this,
                "What rank should the colorization be based on?",
                "Rank Selection",
                JOptionPane.QUESTION_MESSAGE,
                null,
                ranks,
                ranks[0]);
        if (ForesterUtil.isEmpty(rank)) {
            return;
        }
        if (rank.indexOf('(') > 0) {
            rank = rank.substring(0, rank.indexOf('(')).trim();
        }
        // Decide whether the DB is needed BEFORE colorizing (cache-only, never blocks the EDT), so we
        // colorize the tree exactly once -- either here, or (on the online path) in the resolver after
        // the fetch -- rather than colorizing now and again, which flashed a partial result.
        final String r = rank;
        final SortedSet<String> unresolved = TreePanelUtil.unresolvedTipTaxa(phy, r,
                TreePanelUtil.getDefaultLineageService());
        if (!unresolved.isEmpty()) {
            final int choice = JOptionPane.showConfirmDialog(this,
                    unresolved.size() + " tip " + ((unresolved.size() == 1) ? "taxon lacks" : "taxa lack")
                            + " a \"" + r + "\"-rank identity in the tree itself.\n"
                            + "Resolve online via the NCBI and UniProt databases so each tip's own identity sets its color"
                            + " (overriding any internal-node annotation)?\n"
                            + "Decline to colorize from the tree's existing annotations. (Requires an internet connection.)",
                    "Resolve Taxa Online?",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.QUESTION_MESSAGE);
            if (choice == JOptionPane.YES_OPTION) {
                new Thread(new OnlineTaxonResolver(this, "rank colorization (" + r + ")", unresolved, err -> {
                    tp.pushUndoCheckpoint("Colorize by Rank"); // checkpoint at the actual (deferred) colorization
                    final int colorized = tp.colorByRank(r);
                    if (err != null) {
                        // colorByRank mutated branch colors but does not set the edited flag; on the error
                        // path we skip reportRankColorization, so mark the tree edited here.
                        if (colorized > 0) {
                            tp.setEdited(true);
                        }
                        JOptionPane.showMessageDialog(this,
                                "Colorized " + colorized + " clade(s), but some taxa could not be resolved:\n" + err,
                                "Taxonomy Rank-Colorization (" + r + ")", JOptionPane.WARNING_MESSAGE);
                    } else {
                        tp.reportRankColorization(r, colorized);
                    }
                })).start();
                return; // the background resolver colorizes and reports when done
            }
        }
        // no online resolution -- colorize once from what the tree (and cache) already know
        tp.pushUndoCheckpoint("Colorize by Rank");
        tp.reportRankColorization(r, tp.colorByRank(r));
    }

    /**
     * The Tools "Annotate Clades by Rank…" operation: pick a rank and a mode (shaded boxes or right-edge
     * bars), resolve any unplaced tips online if the user agrees (off the EDT), then draw the bands.
     * Mirrors {@link #colorRank()} but renders {@link CladeBand}s instead of coloring branches.
     */
    void buildEditMenu() {
        _edit_jmenu = createMenu("Edit", getConfiguration()); // font handled like the sibling top-level menus
        _edit_jmenu.setToolTipText("Undo or redo the last tree edit");
        final int shortcut = Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx(); // Cmd on macOS, Ctrl elsewhere
        _edit_jmenu.add(_undo_item = new JMenuItem("Undo"));
        _undo_item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Z, shortcut));
        customizeJMenuItem(_undo_item);
        _edit_jmenu.add(_redo_item = new JMenuItem("Redo"));
        _redo_item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Z, shortcut | InputEvent.SHIFT_DOWN_MASK));
        customizeJMenuItem(_redo_item);
        _jmenubar.add(_edit_jmenu);
        updateEditMenu();
    }

    /** Undo the current tab's last tree edit, then refresh the view + the Edit menu. */
    void undo() {
        if ((getCurrentTreePanel() != null) && getCurrentTreePanel().undo()) {
            syncSelectedTabTitleToTreeName();
            showWhole();
            getCurrentTreePanel().repaint();
            repaint();
        }
        updateEditMenu();
    }

    /** Re-apply the last undone edit, then refresh the view + the Edit menu. */
    void redo() {
        if ((getCurrentTreePanel() != null) && getCurrentTreePanel().redo()) {
            syncSelectedTabTitleToTreeName();
            showWhole();
            getCurrentTreePanel().repaint();
            repaint();
        }
        updateEditMenu();
    }

    /**
     * After an undo/redo restore, keep the selected tab's label in step with the (possibly reverted) tree name,
     * so a "Edit Tree Name and Description" rename and its undo cannot leave the tab and the tree name drifted.
     * When the restored name is empty (undoing the first naming of a previously-unnamed tree) the tab is
     * re-derived to the same never-empty placeholder it had at load -- the tree file's base name, else "[n]" --
     * rather than left showing the just-undone name (which the save-time backfill would otherwise persist).
     */
    private void syncSelectedTabTitleToTreeName() {
        final TreePanel tp = getCurrentTreePanel();
        if ((tp == null) || (tp.getPhylogeny() == null)) {
            return;
        }
        final MainPanel mp = getMainPanel();
        final String fallback = (tp.getTreeFile() != null) ? tp.getTreeFile().getName() : null;
        mp.setTitleOfSelectedTab(MainPanel.tabTitleFor(tp.getPhylogeny(), fallback,
                                                       mp.getTabbedPane().getSelectedIndex() + 1));
    }

    /** Syncs the Edit menu's enabled state + "Undo <op>"/"Redo <op>" labels to the current tab's history. */
    void updateEditMenu() {
        if (_undo_item == null) {
            return;
        }
        final TreePanel tp = (getMainPanel() != null) ? getMainPanel().getCurrentTreePanel() : null;
        final boolean can_undo = (tp != null) && tp.canUndo();
        final boolean can_redo = (tp != null) && tp.canRedo();
        _undo_item.setEnabled(can_undo);
        _undo_item.setText(can_undo ? ("Undo " + tp.undoLabel()) : "Undo");
        _redo_item.setEnabled(can_redo);
        _redo_item.setText(can_redo ? ("Redo " + tp.redoLabel()) : "Redo");
    }

    /** Refreshes the menu-bar "Found / Selected: N" counter from the current tree's highlighted (found + selected)
     *  nodes; hides it when nothing is highlighted. Called from the found-set choke point
     *  ({@link ControlPanel#updateSearchHitNavigation()}), so it stays in sync with search / selection / prune / undo /
     *  tab change. Null-safe during construction/teardown. */
    void updateFoundSelectedCounter() {
        if (_found_selected_counter == null) {
            return;
        }
        final TreePanel tp = (getMainPanel() != null) ? getMainPanel().getCurrentTreePanel() : null;
        if (tp == null) {
            _found_selected_counter.setCounts(0, 0, 0, 0, false, null);
            return;
        }
        // distinct highlighted nodes = the union of the two found sets = a + b - both; derive it from the one
        // breakdown walk instead of a second getSearchHitCount() preorder pass.
        final int[] br = tp.foundSelectedBreakdown();
        final int total = (br[0] + br[1]) - br[2];
        // found-set 0 doubles as "Search A" and the manual selection; it holds Search-A hits only when box A has a
        // query -- otherwise its nodes are a manual selection (a search would have cleared them), so label them so.
        final ControlPanel cp = tp.getControlPanel();
        final boolean a_is_search = (cp != null) && (cp.getSearchTextField0() != null)
                && (cp.getSearchTextField0().getText() != null)
                && !cp.getSearchTextField0().getText().trim().isEmpty();
        // when the two boxes are combined (A AND/OR B), the result is one set in found-set 0 -- label it by the combine
        // mode instead of the per-box A/B breakdown (which would misreport A∩B/A∪B as box-A-only hits)
        final String combine_label = (cp != null) ? cp.searchCombineDescription() : null;
        _found_selected_counter.setCounts(total, br[0], br[1], br[2], a_is_search, combine_label);
    }

    /**
     * Opens the "Annotation Columns" chooser: pick which annotation fields to show as tip-aligned columns and
     * how to render each (color strip / heat map / bar / text), then applies them to the current tree.
     */
    void chooseAnnotationColumns() {
        if (_mainpanel.getCurrentTreePanel() == null) {
            return;
        }
        final TreePanel tp = _mainpanel.getCurrentTreePanel();
        final Phylogeny phy = tp.getPhylogeny();
        if ((phy == null) || phy.isEmpty()) {
            return;
        }
        final List<String> refs = PropertyColorScheme.colorableRefs(phy);
        if (refs.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "This tree has no annotation fields to show as columns.\n"
                            + "Import a table (File → Import Annotations) or load a tree with node properties first.",
                    "No Annotation Fields", JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        // pre-select whatever columns are shown now
        final Map<String, AnnotationColumns.Type> current = new HashMap<String, AnnotationColumns.Type>();
        if (tp.getAnnotationColumnSpecs() != null) {
            for (final AnnotationColumns.ColumnSpec s : tp.getAnnotationColumnSpecs()) {
                current.put(s._ref, s._type);
            }
        }
        final DefaultListCellRenderer type_renderer = new DefaultListCellRenderer() {

            @Override
            public java.awt.Component getListCellRendererComponent(final JList<?> list, final Object value,
                    final int index, final boolean sel, final boolean focus) {
                super.getListCellRendererComponent(list, value, index, sel, focus);
                if (value instanceof AnnotationColumns.Type) {
                    setText(AnnotationColumns.label((AnnotationColumns.Type) value));
                }
                return this;
            }
        };
        final JPanel panel = new JPanel(new GridLayout(0, 2, 10, 3));
        panel.add(new JLabel("Field"));
        panel.add(new JLabel("Show as"));
        final List<JCheckBox> checks = new ArrayList<JCheckBox>();
        final List<JComboBox<AnnotationColumns.Type>> combos = new ArrayList<JComboBox<AnnotationColumns.Type>>();
        for (final String ref : refs) {
            final JCheckBox cb = new JCheckBox(PropertyColorScheme.displayName(ref), current.containsKey(ref));
            final List<AnnotationColumns.Type> types = AnnotationColumns.allowedTypes(phy, ref);
            final JComboBox<AnnotationColumns.Type> combo = new JComboBox<AnnotationColumns.Type>(
                    types.toArray(new AnnotationColumns.Type[0]));
            combo.setRenderer(type_renderer);
            combo.setSelectedItem(current.containsKey(ref) ? current.get(ref)
                    : AnnotationColumns.defaultType(phy, ref));
            panel.add(cb);
            panel.add(combo);
            checks.add(cb);
            combos.add(combo);
        }
        final JScrollPane sp = new JScrollPane(panel);
        sp.setPreferredSize(new java.awt.Dimension(380, Math.min(440, 50 + (refs.size() * 30))));
        if (JOptionPane.showConfirmDialog(this, sp, "Annotation Columns", JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.PLAIN_MESSAGE) != JOptionPane.OK_OPTION) {
            return;
        }
        final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<AnnotationColumns.ColumnSpec>();
        for (int i = 0; i < refs.size(); ++i) {
            if (checks.get(i).isSelected()) {
                specs.add(new AnnotationColumns.ColumnSpec(refs.get(i),
                        (AnnotationColumns.Type) combos.get(i).getSelectedItem()));
            }
        }
        tp.setAnnotationColumns(specs);
        tp.setEdited(true);
        getControlPanel().fitWidth(); // reveal the columns even when they extend past the current width
        tp.repaint();
        repaint();
    }

    void labelCladesByRank() {
        if (_mainpanel.getCurrentTreePanel() == null) {
            return;
        }
        final TreePanel tp = _mainpanel.getCurrentTreePanel();
        final Phylogeny phy = tp.getPhylogeny();
        if ((phy == null) || phy.isEmpty() || (phy.getNumberOfExternalNodes() < 2)) {
            return;
        }
        final String[] ranks = AptxUtil.getRankChoices(AptxUtil.getRankCounts(phy),
                AptxUtil.getRankCoverageCounts(phy), phy.getNumberOfExternalNodes());
        final JComboBox<String> rank_box = new JComboBox<>(ranks);
        preselectLastCladeRank(rank_box, ranks);
        final JRadioButton boxes_rb = new JRadioButton("Shaded boxes (behind the clades)", true);
        final JRadioButton bars_rb = new JRadioButton("Bars + labels (at the right edge)");
        final JRadioButton brackets_rb = new JRadioButton("Brackets ] + labels (black & white, no legend)");
        final ButtonGroup bg = new ButtonGroup();
        bg.add(boxes_rb);
        bg.add(bars_rb);
        bg.add(brackets_rb);
        final JCheckBox write_cb = new JCheckBox("Also write the clade taxa into the tree (rank + NCBI id; undoable)",
                false);
        final JCheckBox overwrite_cb = new JCheckBox("    ...overwriting existing internal-node taxonomies", false);
        final JPanel panel = new JPanel(new GridLayout(0, 1, 0, 2));
        panel.add(new JLabel("Annotate clades by rank:"));
        panel.add(rank_box);
        panel.add(new JLabel(" "));
        panel.add(new JLabel("Show as:"));
        panel.add(boxes_rb);
        panel.add(bars_rb);
        panel.add(brackets_rb);
        panel.add(new JLabel(" "));
        panel.add(write_cb);
        panel.add(overwrite_cb);
        if (JOptionPane.showConfirmDialog(this, panel, "Annotate Clades by Rank", JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.PLAIN_MESSAGE) != JOptionPane.OK_OPTION) {
            return;
        }
        final boolean write = write_cb.isSelected();
        final boolean overwrite = write && overwrite_cb.isSelected();
        String rank = (String) rank_box.getSelectedItem();
        if (ForesterUtil.isEmpty(rank)) {
            return;
        }
        if (rank.indexOf('(') > 0) {
            rank = rank.substring(0, rank.indexOf('(')).trim();
        }
        final String r = rank;
        _last_clade_rank = r; // remember for the next invocation's pre-selection
        final TreePanel.CLADE_VIS mode = bars_rb.isSelected() ? TreePanel.CLADE_VIS.BARS
                : brackets_rb.isSelected() ? TreePanel.CLADE_VIS.BRACKETS : TreePanel.CLADE_VIS.BOXES;
        final SortedSet<String> unresolved = TreePanelUtil.unresolvedTipTaxa(phy, r,
                TreePanelUtil.getDefaultLineageService());
        if (!unresolved.isEmpty()) {
            final int choice = JOptionPane.showConfirmDialog(this,
                    unresolved.size() + " tip " + ((unresolved.size() == 1) ? "taxon lacks" : "taxa lack")
                            + " a \"" + r + "\"-rank identity in the tree itself.\n"
                            + "Resolve online via the NCBI and UniProt databases so each tip's own identity marks its clade"
                            + " (overriding any internal-node annotation)?\n"
                            + "Decline to annotate from the tree's existing annotations. (Requires an internet connection.)",
                    "Resolve Taxa Online?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
            if (choice == JOptionPane.YES_OPTION) {
                new Thread(new OnlineTaxonResolver(this, "clade bands (" + r + ")", unresolved,
                        err -> reportCladeBands(tp, r, mode, write, overwrite, err))).start();
                return;
            }
        }
        reportCladeBands(tp, r, mode, write, overwrite, null);
    }

    /**
     * Pre-selects in {@code rank_box} the rank last used in "Annotate Clades by Rank" (this session), so a
     * repeat invocation defaults to the same rank. A no-op on first use, or when the remembered rank is
     * absent from the current tree's choices.
     */
    private void preselectLastCladeRank(final JComboBox<String> rank_box, final String[] ranks) {
        final int idx = AptxUtil.indexOfRank(ranks, _last_clade_rank);
        if (idx >= 0) {
            rank_box.setSelectedIndex(idx);
        }
    }

    private void reportCladeBands(final TreePanel tp, final String rank, final TreePanel.CLADE_VIS mode,
                                  final boolean write, final boolean overwrite, final String error) {
        int wrote = 0;
        if (write) {
            tp.pushUndoCheckpoint("Annotate Clade Taxa"); // a tree-data mutation -> checkpoint before writing
            wrote = tp.writeCladeTaxonomiesByRank(rank, overwrite);
        }
        final int n = tp.setCladeBands(rank, mode);
        if (n > 0) {
            tp.setEdited(true);
            // bars/brackets extend to the right of the labels; fit the width so they are immediately
            // visible without the user having to press "W" (boxes sit within the existing tree width)
            if ((mode == TreePanel.CLADE_VIS.BARS) || (mode == TreePanel.CLADE_VIS.BRACKETS)) {
                tp.getControlPanel().fitWidth();
            }
        }
        final String kind = (mode == TreePanel.CLADE_VIS.BARS) ? "bar(s)"
                : (mode == TreePanel.CLADE_VIS.BRACKETS) ? "bracket(s)" : "box(es)";
        final StringBuilder msg = new StringBuilder("Drew ").append(n).append(" clade ").append(kind)
                .append(" at rank \"").append(rank).append("\".");
        if (write) {
            msg.append("\nWrote a taxonomy onto ").append(wrote).append(" internal clade node")
                    .append(wrote == 1 ? "" : "s").append(" (rank + NCBI id; saved with the tree).");
        }
        if (error != null) {
            msg.append("\nSome taxa could not be resolved:\n").append(error);
            JOptionPane.showMessageDialog(this, msg.toString(), "Annotate Clades by Rank (" + rank + ")",
                    JOptionPane.WARNING_MESSAGE);
        } else if ((n > 0) || (wrote > 0)) {
            JOptionPane.showMessageDialog(this, msg.toString(), "Annotate Clades by Rank (" + rank + ")",
                    JOptionPane.INFORMATION_MESSAGE);
        } else {
            JOptionPane.showMessageDialog(this, "Could not place any tip at rank \"" + rank + "\".\n"
                    + "Try a different rank, or check that the tips carry resolvable taxonomic names.",
                    "Annotate Clades by Rank (" + rank + ")", JOptionPane.WARNING_MESSAGE);
        }
    }

    /**
     * Just after a tree is loaded, if most of its labels are UniProt FASTA headers, offer once to extract
     * their data (the proactive half of the feature; the Tools menu item is the on-demand half). A no-op
     * for ordinary trees, so it never nags.
     */
    void offerLabelExtraction(final Phylogeny[] phys) {
        if (phys == null) {
            return;
        }
        boolean offer = false;
        for (final Phylogeny p : phys) {
            if (LabelDataExtractor.mostLabelsParsable(p)) {
                offer = true;
                break;
            }
        }
        if (!offer) {
            return;
        }
        final int choice = JOptionPane.showConfirmDialog(this,
                "These node labels look like UniProt or GenBank FASTA headers.\n"
                        + "Extract their accession, description, gene and taxonomy into proper fields?\n"
                        + "(Node names are shortened to the accession; only empty fields are filled.)",
                "Extract Data from Labels?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
        if (choice == JOptionPane.YES_OPTION) {
            extractLabelData();
        }
    }

    /**
     * Load-time OFFER (mirrors {@link #offerLabelExtraction}): when the CURRENT tree's tip labels look like they carry
     * sampling dates and it has no structured {@code <date>} yet, ask whether to extract them. Yes opens the preview
     * dialog ({@link #extractTipDates}). Runs before the ultrametric time-tree offer, so accepting it (which dates the
     * tree) preempts a redundant "treat as time tree" prompt.
     */
    void offerTipDateExtraction() {
        if ((_mainpanel == null) || (_mainpanel.getCurrentTreePanel() == null)) {
            return;
        }
        if (!TipDateExtractor.shouldOffer(_mainpanel.getCurrentTreePanel().getPhylogeny())) {
            return;
        }
        final int choice = JOptionPane.showConfirmDialog(this,
                "The tip labels look like they contain sampling dates (e.g. \"…/2021-03-15\").\n"
                        + "Extract them into node dates? You'll see a preview first, and it's undoable.",
                "Extract Dates from Labels?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
        if (choice == JOptionPane.YES_OPTION) {
            extractTipDates();
        }
    }

    /**
     * Load-time offer (NH/NHX/Nexus only, and only when the manual "Internal Node Names are Confidence
     * Values" option is off): if the internal node labels look like support values
     * ({@link AptxUtil#internalNamesLookLikeConfidenceValues}), ask whether to treat them as confidence
     * values rather than node names.
     */
    /** After a load, OFFER to label any merely-ultrametric (undated) loaded tree a time tree. A DATED tree is
     *  auto-labeled by the badge itself and needs no dialog; this handles only the ambiguous case (ultrametric
     *  could also be a UPGMA distance tree, so we ask rather than assert). One dialog covers all loaded tabs. */
    void offerTreatAsTimeTree() {
        if (_mainpanel == null) {
            return;
        }
        final java.util.List<TreePanel> panels = _mainpanel.getTreePanels();
        if ((panels == null) || panels.isEmpty()) {
            return;
        }
        final java.util.List<TreePanel> ultrametric = new java.util.ArrayList<>();
        for (final TreePanel tp : panels) {
            if ((tp != null)
                    && (AptxUtil.detectTimeTree(tp.getPhylogeny()) == AptxUtil.TIME_TREE_KIND.ULTRAMETRIC)) {
                ultrametric.add(tp);
            }
        }
        if (ultrametric.isEmpty()) {
            return;
        }
        final String what = (ultrametric.size() == 1) ? "This tree is"
                : (ultrametric.size() + " of the loaded trees are");
        final int choice = JOptionPane.showConfirmDialog(this,
                what + " ultrametric (all tips the same distance from the root), which is typical of a\n"
                        + "dated / time tree. Label and treat "
                        + ((ultrametric.size() == 1) ? "it" : "them") + " as " + ((ultrametric.size() == 1) ? "a" : "")
                        + " time " + ((ultrametric.size() == 1) ? "tree" : "trees") + "?\n\n"
                        + "(If the branch lengths are genetic distance rather than time, choose No.)",
                "Tree Looks Like a Time Tree", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
        if (choice == JOptionPane.YES_OPTION) {
            for (final TreePanel tp : ultrametric) {
                tp.setConfirmedTimeTree(true);
            }
        }
    }

    void offerInternalNamesAsConfidence(final Phylogeny[] phys) {
        if ((phys == null) || getOptions().isInternalNumberAreConfidenceForNhParsing()) {
            return;
        }
        boolean offer = false;
        for (final Phylogeny p : phys) {
            if (AptxUtil.internalNamesLookLikeConfidenceValues(p)) {
                offer = true;
                break;
            }
        }
        if (!offer) {
            return;
        }
        final int choice = JOptionPane.showConfirmDialog(this,
                "The internal node labels look like support / confidence values (e.g. bootstrap).\n"
                        + "Treat them as confidence values instead of node names?",
                "Internal Labels Look Like Support Values", JOptionPane.YES_NO_OPTION,
                JOptionPane.QUESTION_MESSAGE);
        if (choice == JOptionPane.YES_OPTION) {
            for (final Phylogeny p : phys) {
                AptxUtil.stripBracketsFromInternalNames(p);
                PhylogenyMethods.transferInternalNodeNamesToConfidence(p, "");
            }
            if (_mainpanel.getCurrentTreePanel() != null) {
                _mainpanel.getCurrentTreePanel().setEdited(true);
                // reveal the freshly-created confidence values so the conversion is immediately visible
                _mainpanel.getControlPanel().setCheckbox(DisplayOption.WRITE_CONFIDENCE_VALUES, true);
                _mainpanel.getControlPanel().displayedPhylogenyMightHaveChanged(true);
                _mainpanel.getCurrentTreePanel().repaint();
            }
        }
    }

    /**
     * The Tools "Extract Data from Labels…" operation: parse UniProt FASTA-header node names into proper
     * sequence + taxonomy fields (offline, only filling empties), shorten the names to their accession,
     * then reveal "Seq Name" + "Taxonomy Scientific" so the cleaned-up labels are visible.
     */
    void extractLabelData() {
        if (_mainpanel.getCurrentTreePanel() == null) {
            return;
        }
        final TreePanel tp = _mainpanel.getCurrentTreePanel();
        final Phylogeny phy = tp.getPhylogeny();
        if ((phy == null) || phy.isEmpty()) {
            return;
        }
        final int n = extractLabelDataAndRefit(tp);
        if (n > 0) {
            JOptionPane.showMessageDialog(this,
                    "Extracted accession, description and taxonomy from " + n + " label" + ((n == 1) ? "" : "s")
                            + ".\nNode names were shortened to their accession; \"Seq Name\" and "
                            + "\"Taxonomy Scientific\" are now shown.",
                    "Extract Data from Labels", JOptionPane.INFORMATION_MESSAGE);
        } else {
            JOptionPane.showMessageDialog(this,
                    "No UniProt or GenBank FASTA-header labels (e.g. \"tr|ACC|ENTRY … OS=… OX=…\" or "
                            + "\"NP_000537.1 … [Homo sapiens]\") were found to extract from.",
                    "Extract Data from Labels", JOptionPane.INFORMATION_MESSAGE);
        }
    }

    /**
     * The dialog-free core of {@link #extractLabelData()} (so it is unit-testable): parses the panel's tree
     * labels, and -- if anything was extracted -- reveals the freshly-populated "Seq Name" + "Taxonomy
     * Scientific" columns and re-fits the tree width. Returns the number of labels changed.
     * <p>
     * The re-fit ({@link ControlPanel#fitWidth()}) is the point: shortening the names and revealing those two
     * columns changes the horizontal extent (and, on domain-bearing trees, shifts the lined-up domain
     * architectures rightward). {@code fitWidth()} recalcs the longest-ext-node info AND resets the preferred
     * size -- which the plain {@link ControlPanel#displayedPhylogenyMightHaveChanged(boolean)} path does not,
     * leaving the panel width stale so the rightmost domain lands at/just past the panel edge and is clipped
     * until the user manually zooms.
     */
    int extractLabelDataAndRefit(final TreePanel tp) {
        final int n = LabelDataExtractor.extract(tp.getPhylogeny());
        if (n > 0) {
            tp.setEdited(true);
            tp.getControlPanel().setCheckbox(DisplayOption.SHOW_SEQ_NAMES, true);
            tp.getControlPanel().setCheckbox(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES, true);
            // the extraction populated new fields (seq name/accession/gene name, taxonomy) on the SAME tree, so force
            // a data-presence refresh: this reveals the newly-populated Display-Data checkboxes AND (via the shared
            // choke point) rebuilds the search field selectors so the new fields are searchable -- on the manual Tools
            // path too, not only the load auto-offer path.
            tp.getControlPanel().updateDataCheckboxVisibility(true);
            tp.getControlPanel().fitWidth();
        }
        return n;
    }

    /**
     * The Tools "Extract Dates from Labels…" operation: parse a sampling date out of each tip label (ISO / numeric /
     * month-name / decimal-year / bare-year) via a PREVIEW dialog, then -- on Apply -- write it into each tip's
     * {@code <date>} + a numeric {@code data:date} property, so a tip-dated Newick/phyloXML/Nexus tree lights up the
     * Calendar axis and Color-by-date. A no-op if no label carries a recognisable date.
     */
    void extractTipDates() {
        if (_mainpanel.getCurrentTreePanel() == null) {
            return;
        }
        final TreePanel tp = _mainpanel.getCurrentTreePanel();
        final Phylogeny phy = tp.getPhylogeny();
        if ((phy == null) || phy.isEmpty()) {
            return;
        }
        boolean any = false;
        for (final TipDate t : TipDateExtractor.preview(phy, DayMonthOrder.DAY_FIRST)) {
            if (t.match() != null) {
                any = true;
                break;
            }
        }
        if (!any) {
            JOptionPane.showMessageDialog(this,
                    "No sampling dates were recognised in the tip labels.\nTried ISO (2021-03-15), numeric "
                            + "(15/03/2021), month-name (01-Dec-2015), decimal-year (2021.37) and bare-year (…/2012) "
                            + "formats.",
                    "Extract Dates from Labels", JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        final TipDateExtractionDialog dlg = new TipDateExtractionDialog(this, phy);
        dlg.setVisible(true); // modal
        if (dlg.isApplied()) {
            final int n = applyTipDatesAndRefit(tp, dlg.getOrder(), dlg.isSkipExisting());
            if (n > 0) {
                JOptionPane.showMessageDialog(this,
                        "Wrote sampling dates to " + n + " tip" + ((n == 1) ? "" : "s")
                                + ".\nThe tree is now on the Calendar axis, and \"Color by\" → data:date shows the "
                                + "date gradient.",
                        "Extract Dates from Labels", JOptionPane.INFORMATION_MESSAGE);
            } else {
                JOptionPane.showMessageDialog(this,
                        "No new dates were written — every tip with a recognised date already has one.\n"
                                + "Uncheck \"Skip tips that already have a date\" to overwrite them.",
                        "Extract Dates from Labels", JOptionPane.INFORMATION_MESSAGE);
            }
        }
    }

    /**
     * The dialog-free core of {@link #extractTipDates()} (so it is unit-testable): writes the parsed date onto each
     * matched tip (skipping already-dated tips when requested), appends a provenance sentence, and re-fits so the
     * auto-derived Calendar axis + the {@code data:date} Color-by option appear. Undoable (a checkpoint is pushed
     * before mutating). Returns the number of tips dated.
     */
    int applyTipDatesAndRefit(final TreePanel tp, final DayMonthOrder order, final boolean skip_existing) {
        final Phylogeny phy = tp.getPhylogeny();
        if ((phy == null) || phy.isEmpty()) {
            return 0;
        }
        final java.util.List<TipDate> rows = TipDateExtractor.preview(phy, order);
        int to_write = 0;
        for (final TipDate t : rows) {
            if ((t.match() != null) && !(skip_existing && t.alreadyDated())) {
                to_write++;
            }
        }
        if (to_write == 0) {
            return 0;
        }
        tp.pushUndoCheckpoint("Extract dates from labels"); // reversible: a tree-data mutation
        int written = 0;
        for (final TipDate t : rows) {
            if ((t.match() != null) && !(skip_existing && t.alreadyDated())
                    && TipDateExtractor.applyToNode(t.node(), t.match())) {
                written++;
            }
        }
        final Summary sum = TipDateExtractor.summarize(rows);
        final String prov = TipDateExtractor.provenanceSentence(written, rows.size(), sum.dominantFormat());
        final String existing = phy.getDescription();
        phy.setDescription(ForesterUtil.isEmpty(existing) ? prov : existing + " " + prov);
        tp.setEdited(true);
        phy.externalNodesHaveChanged();
        tp.invalidateTimeAxisDerivation(); // dates were added in place on the same tree -> re-derive the Calendar axis
        // reveal the data:date Color-by option + let the Calendar axis auto-derive from the new "year" dates; re-fit
        // (null-guarded: a minimal/embedded TreePanel may have no control panel -- the tree is still validly dated)
        if (tp.getControlPanel() != null) {
            tp.getControlPanel().updateDataCheckboxVisibility(true);
            tp.getControlPanel().populateColorByPropertyBox();
            tp.getControlPanel().showWhole();
        }
        return written;
    }


    void customizeCheckBoxMenuItem(final JCheckBoxMenuItem item, final boolean is_selected) {
        if (item != null) {
            item.setFont(MainFrame.menu_font);
            item.setSelected(is_selected);
            item.addActionListener(this);
        }
    }

    JMenuItem customizeJMenuItem(final JMenuItem jmi) {
        if (jmi != null) {
            jmi.setFont(MainFrame.menu_font);
            jmi.addActionListener(this);
        }
        return jmi;
    }

    void customizeRadioButtonMenuItem(final JRadioButtonMenuItem item, final boolean is_selected) {
        if (item != null) {
            item.setFont(MainFrame.menu_font);
            item.setSelected(is_selected);
            item.addActionListener(this);
        }
    }

    /** Opens the small modal editor for the current tree's name and description (View menu / tab menu). */
    void showTreeInfoDialog() {
        TreeInfoDialog.showFor(this);
    }

    void displayBasicInformation(final File treefile) {
        if ((_mainpanel.getCurrentPhylogeny() != null) && !_mainpanel.getCurrentPhylogeny().isEmpty()) {
            String title = "Basic Information";
            if (!ForesterUtil.isEmpty(_mainpanel.getCurrentPhylogeny().getName())) {
                title = title + " for \"" + _mainpanel.getCurrentPhylogeny().getName() + "\"";
            }
            showTextFrame(AptxUtil.createBasicInformation(_mainpanel.getCurrentPhylogeny(), treefile), title);
        }
    }

    void exceptionOccuredDuringOpenFile(final Exception e) {
        try {
            _mainpanel.getCurrentTreePanel().setArrowCursor();
        } catch (final Exception ex) {
            // Do nothing.
        }
        JOptionPane.showMessageDialog(this,
                ForesterUtil.wordWrap(e.getLocalizedMessage(), 80),
                "Error during File|Open",
                JOptionPane.ERROR_MESSAGE);
    }

    void executeGSDI() {
        if (!isOKforSDI(false, true)) {
            return;
        }
        if (!_mainpanel.getCurrentPhylogeny().isRooted()) {
            JOptionPane.showMessageDialog(this,
                    "Gene tree is not rooted.",
                    "Cannot execute GSDI",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        final Phylogeny gene_tree = _mainpanel.getCurrentPhylogeny().copy();
        gene_tree.setAllNodesToNotCollapse();
        gene_tree.recalculateNumberOfExternalDescendants(false);
        GSDI gsdi = null;
        final Phylogeny species_tree = getSpeciesTree().copy();
        try {
            gsdi = new GSDI(gene_tree, species_tree, false, true, true, true);
        } catch (final SDIException e) {
            JOptionPane.showMessageDialog(this,
                    e.getLocalizedMessage(),
                    "Error during GSDI",
                    JOptionPane.ERROR_MESSAGE);
            return;
        } catch (final Exception e) {
            AptxUtil.unexpectedException(e);
            return;
        }
        gene_tree.setRerootable(false);
        gene_tree.clearHashIdToNodeMap();
        gene_tree.recalculateNumberOfExternalDescendants(true);
        PhylogenyMethods.removeMadConfidences(gene_tree); // the result is a reconciliation, not a MAD rooting
        _mainpanel.addPhylogenyInNewTab(gene_tree, getConfiguration(), "gene tree", null);
        getMainPanel().getControlPanel().setShowEvents(true);
        showWhole();
        final int selected = _mainpanel.getTabbedPane().getSelectedIndex();
        _mainpanel.addPhylogenyInNewTab(species_tree, getConfiguration(), "species tree", null);
        showWhole();
        _mainpanel.getTabbedPane().setSelectedIndex(selected);
        showWhole();
        _mainpanel.getCurrentTreePanel().setEdited(true);
        final int poly = PhylogenyMethods.countNumberOfPolytomies(species_tree);
        if (gsdi.getStrippedExternalGeneTreeNodes().size() > 0) {
            JOptionPane.showMessageDialog(this,
                    "Duplications: " + gsdi.getDuplicationsSum() + "\n"
                            + "Potential duplications: "
                            + gsdi.getSpeciationOrDuplicationEventsSum() + "\n"
                            + "Speciations: " + gsdi.getSpeciationsSum() + "\n"
                            + "Stripped gene tree nodes: "
                            + gsdi.getStrippedExternalGeneTreeNodes().size() + "\n"
                            + "Taxonomy linkage based on: " + gsdi.getTaxCompBase() + "\n"
                            + "Number of polytomies in species tree used: " + poly + "\n",
                    "GSDI successfully completed",
                    JOptionPane.WARNING_MESSAGE);
        } else {
            JOptionPane.showMessageDialog(this,
                    "Duplications: " + gsdi.getDuplicationsSum() + "\n"
                            + "Potential duplications: "
                            + gsdi.getSpeciationOrDuplicationEventsSum() + "\n"
                            + "Speciations: " + gsdi.getSpeciationsSum() + "\n"
                            + "Stripped gene tree nodes: "
                            + gsdi.getStrippedExternalGeneTreeNodes().size() + "\n"
                            + "Taxonomy linkage based on: " + gsdi.getTaxCompBase() + "\n"
                            + "Number of polytomies in species tree used: " + poly + "\n",
                    "GSDI successfully completed",
                    JOptionPane.INFORMATION_MESSAGE);
        }
    }

    void executeGSDIR() {
        if (!isOKforSDI(false, false)) {
            return;
        }
        final int p = PhylogenyMethods.countNumberOfPolytomies(_mainpanel.getCurrentPhylogeny());
        if ((p > 0)
                && !((p == 1) && (_mainpanel.getCurrentPhylogeny().getRoot().getNumberOfDescendants() == 3))) {
            JOptionPane.showMessageDialog(this,
                    "Gene tree is not completely binary",
                    "Cannot execute GSDI",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        final Phylogeny gene_tree = _mainpanel.getCurrentPhylogeny().copy();
        gene_tree.setAllNodesToNotCollapse();
        gene_tree.recalculateNumberOfExternalDescendants(false);
        GSDIR gsdir = null;
        final Phylogeny species_tree = getSpeciesTree().copy();
        try {
            gsdir = new GSDIR(gene_tree, species_tree, true, true, true);
        } catch (final SDIException e) {
            JOptionPane.showMessageDialog(this,
                    e.getLocalizedMessage(),
                    "Error during GSDIR",
                    JOptionPane.ERROR_MESSAGE);
            return;
        } catch (final Exception e) {
            AptxUtil.unexpectedException(e);
            return;
        }
        final Phylogeny result_gene_tree = gsdir.getMinDuplicationsSumGeneTree();
        result_gene_tree.setRerootable(false);
        result_gene_tree.clearHashIdToNodeMap();
        result_gene_tree.recalculateNumberOfExternalDescendants(true);
        PhylogenyMethods.removeMadConfidences(result_gene_tree); // GSDIR rerooted the tree; MAD support is stale
        PhylogenyMethods.orderAppearance(result_gene_tree.getRoot(), true, true, DESCENDANT_SORT_PRIORITY.NODE_NAME);
        _mainpanel.addPhylogenyInNewTab(result_gene_tree, getConfiguration(), "gene tree", null);
        getMainPanel().getControlPanel().setShowEvents(true);
        showWhole();
        final int selected = _mainpanel.getTabbedPane().getSelectedIndex();
        _mainpanel.addPhylogenyInNewTab(species_tree, getConfiguration(), "species tree", null);
        showWhole();
        _mainpanel.getTabbedPane().setSelectedIndex(selected);
        showWhole();
        _mainpanel.getCurrentTreePanel().setEdited(true);
        final int poly = PhylogenyMethods.countNumberOfPolytomies(species_tree);
        if (gsdir.getStrippedExternalGeneTreeNodes().size() > 0) {
            JOptionPane.showMessageDialog(this,
                    "Minimal duplications: " + gsdir.getMinDuplicationsSum() + "\n"
                            + "Speciations: " + gsdir.getSpeciationsSum() + "\n"
                            + "Stripped gene tree nodes: "
                            + gsdir.getStrippedExternalGeneTreeNodes().size() + "\n"
                            + "Taxonomy linkage based on: " + gsdir.getTaxCompBase() + "\n"
                            + "Number of polytomies in species tree used: " + poly + "\n",
                    "GSDIR successfully completed",
                    JOptionPane.WARNING_MESSAGE);
        } else {
            JOptionPane.showMessageDialog(this,
                    "Minimal duplications: " + gsdir.getMinDuplicationsSum() + "\n"
                            + "Speciations: " + gsdir.getSpeciationsSum() + "\n"
                            + "Stripped gene tree nodes: "
                            + gsdir.getStrippedExternalGeneTreeNodes().size() + "\n"
                            + "Taxonomy linkage based on: " + gsdir.getTaxCompBase() + "\n"
                            + "Number of polytomies in species tree used: " + poly + "\n",
                    "GSDIR successfully completed",
                    JOptionPane.INFORMATION_MESSAGE);
        }
    }

    void executeLineageInference() {
        if ((_mainpanel.getCurrentPhylogeny() == null) || (_mainpanel.getCurrentPhylogeny().isEmpty())) {
            return;
        }
        if (!_mainpanel.getCurrentPhylogeny().isRooted()) {
            JOptionPane.showMessageDialog(this,
                    "Phylogeny is not rooted.",
                    "Cannot infer ancestral taxonomies",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        final JCheckBox overwrite_cb = new JCheckBox("Overwrite existing internal-node taxonomies", false);
        final Object[] message = {
                "Assign each internal node the deepest taxon shared by all of its descendant tips.",
                "By default, internal nodes that already carry a taxonomy are kept.",
                overwrite_cb };
        final int opt = JOptionPane.showConfirmDialog(this, message, "Infer Ancestor Taxonomies",
                JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
        if (opt != JOptionPane.OK_OPTION) {
            return;
        }
        final boolean overwrite = overwrite_cb.isSelected();
        final TreePanel tp = _mainpanel.getCurrentTreePanel();
        // Work on a copy: the pure engine mutates internal-node taxa, and the undo checkpoint (taken in
        // AncestralTaxonomyInferrer.commit) snapshots the still-untouched LIVE tree before the copy is installed.
        final Phylogeny phy = _mainpanel.getCurrentPhylogeny().copy();
        final TaxonomicLineageService service = TreePanelUtil.getDefaultLineageService();
        final SortedSet<String> unresolved = TreePanelUtil.tipsWithoutLineage(phy, service);
        if (!unresolved.isEmpty()) {
            final int choice = JOptionPane.showConfirmDialog(this,
                    unresolved.size() + " tip " + ((unresolved.size() == 1) ? "taxon has" : "taxa have")
                            + " no lineage in the tree itself.\n"
                            + "Resolve their lineages online via the NCBI and UniProt databases so their ancestors"
                            + " can be inferred?\n"
                            + "Decline to infer only from the lineages the tree already carries."
                            + " (Requires an internet connection.)",
                    "Resolve Taxa Online?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
            if (choice == JOptionPane.YES_OPTION) {
                new Thread(new OnlineTaxonResolver(this, "ancestral taxonomy", unresolved,
                        err -> runAncestralInference(tp, phy, service, overwrite, err))).start();
                return; // the background resolver infers and reports when done
            }
        }
        runAncestralInference(tp, phy, service, overwrite, null);
    }

    /** EDT-only: build the tip lineages (cache + on-tree), run the pure inference engine on {@code phy}, and (when
     *  anything was assigned) commit it undoably; then report. {@code fetch_error} is non-null if an online resolve
     *  partially failed. */
    private void runAncestralInference(final TreePanel tp,
                                       final Phylogeny phy,
                                       final TaxonomicLineageService service,
                                       final boolean overwrite,
                                       final String fetch_error) {
        final Map<PhylogenyNode, TaxonLineage> tip_lineages = TreePanelUtil.tipLineages(phy, service);
        final AncestralTaxonomyInference.InferenceResult result = AncestralTaxonomyInference
                .inferInternalTaxonomies(phy, tip_lineages, overwrite);
        if (result.getAssigned() > 0) {
            new AncestralTaxonomyInferrer(this, tp, phy, overwrite).commit(result.getAssigned());
        }
        reportAncestralInference(result, fetch_error);
    }

    private void reportAncestralInference(final AncestralTaxonomyInference.InferenceResult result,
                                          final String fetch_error) {
        String msg = AncestralTaxonomyInferrer.inferenceSummary(result.getAssigned(), result.getSkippedExisting(),
                result.getUnresolvedTips());
        if (fetch_error != null) {
            msg += "\n\nSome taxa could not be resolved online:\n" + fetch_error;
            JOptionPane.showMessageDialog(this, msg, "Ancestral Taxonomy Inference", JOptionPane.WARNING_MESSAGE);
        }
        else {
            JOptionPane.showMessageDialog(this, msg, "Ancestral Taxonomy Inference Completed",
                    JOptionPane.INFORMATION_MESSAGE);
        }
    }

    boolean GAndSDoHaveMoreThanOneSpeciesInComman(final Phylogeny gene_tree) {
        if ((gene_tree == null) || gene_tree.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "Gene tree and species tree have no species in common.",
                    "Error during SDI",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        } else if (gene_tree.getNumberOfExternalNodes() < 2) {
            JOptionPane.showMessageDialog(this,
                    "Gene tree and species tree have only one species in common.",
                    "Error during SDI",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        } else {
            return true;
        }
    }

    ControlPanel getControlPanel() {
        return getMainPanel().getControlPanel();
    }

    /** The directory the given dialog category should open in: the last one used for that purpose if it is still
     *  readable, otherwise the home/Desktop fallback. The fallback is NOT stored, so only user-chosen directories
     *  are ever remembered/persisted. */
    File getCurrentDir(final DirectoryPreferences.Category cat) {
        final File dir = _current_dirs.get(cat);
        if ((dir != null) && dir.canRead()) { // a persisted / just-used dir that still resolves (validated here,
            return dir;                       // lazily, not at startup -- so a stale network path can't hang launch)
        }
        if ((_startup_dir != null) && _startup_dir.canRead()) {
            return _startup_dir; // transient launch default (never persisted)
        }
        return fallbackDir();
    }

    private static File fallbackDir() {
        File dir = null;
        if (ForesterUtil.isWindows()) {
            try {
                dir = new File(WindowsUtils.getCurrentUserDesktopPath());
            } catch (final Exception e) {
                dir = null;
            }
        }
        if ((dir == null) || !dir.canRead()) {
            if (System.getProperty("user.home") != null) {
                dir = new File(System.getProperty("user.home"));
            } else if (System.getProperty("user.dir") != null) {
                dir = new File(System.getProperty("user.dir"));
            }
        }
        return dir;
    }

    TreePanel getCurrentTreePanel() {
        return getMainPanel().getCurrentTreePanel();
    }

    /** Test hook: the View-menu "Find Next" item (⌘G), so a test can doClick() its accelerator dispatch. */
    JMenuItem getFindNextHitItemForTest() {
        return _find_next_hit_item;
    }

    /** Test hook: the View-menu "Find Previous" item (⌘⇧G). */
    JMenuItem getFindPreviousHitItemForTest() {
        return _find_prev_hit_item;
    }

    JMenu getHelpMenu() {
        return _help_jmenu;
    }


    final Phylogeny getSpeciesTree() {
        return _species_tree;
    }

    void initializeTypeMenu(final Options options) {
        setTypeMenuToAllUnselected();
        switch (options.getPhylogenyGraphicsType()) {
            case EURO_STYLE:
                _euro_type_cbmi.setSelected(true);
                break;
            case ROUNDED:
                _rounded_type_cbmi.setSelected(true);
                break;
            case TRIANGULAR:
                _triangular_type_cbmi.setSelected(true);
                break;
            case UNROOTED:
                _unrooted_type_cbmi.setSelected(true);
                break;
            case CIRCULAR:
                _circular_type_cbmi.setSelected(true);
                break;
            default:
                _rectangular_type_cbmi.setSelected(true);
                break;
        }
    }

    boolean isOKforSDI(final boolean species_tree_has_to_binary, final boolean gene_tree_has_to_binary) {
        if ((_mainpanel.getCurrentPhylogeny() == null) || _mainpanel.getCurrentPhylogeny().isEmpty()) {
            return false;
        } else if ((getSpeciesTree() == null) || getSpeciesTree().isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "No species tree loaded",
                    "Cannot execute GSDI",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        } else if (species_tree_has_to_binary && !getSpeciesTree().isCompletelyBinary()) {
            JOptionPane.showMessageDialog(this,
                    "Species tree is not completely binary",
                    "Cannot execute GSDI",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        } else if (gene_tree_has_to_binary && !_mainpanel.getCurrentPhylogeny().isCompletelyBinary()) {
            JOptionPane.showMessageDialog(this,
                    "Gene tree is not completely binary",
                    "Cannot execute GSDI",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        } else {
            return true;
        }
    }

    boolean isSubtreeDisplayed() {
        if (getCurrentTreePanel() != null) {
            if (getCurrentTreePanel().isCurrentTreeIsSubtree()) {
                JOptionPane
                        .showMessageDialog(this,
                                "This operation can only be performed on a complete tree, not on the currently displayed sub-tree only.",
                                "Operation can not be exectuted on a sub-tree",
                                JOptionPane.WARNING_MESSAGE);
                return true;
            }
        }
        return false;
    }

    void midpointRoot() {
        if (_mainpanel.getCurrentTreePanel() != null) {
            _mainpanel.getCurrentTreePanel().midpointRoot();
        }
    }

    void madRoot() {
        if (_mainpanel.getCurrentTreePanel() != null) {
            _mainpanel.getCurrentTreePanel().madRoot();
        }
    }

    void removeAllTextFrames() {
        for (final TextFrame tf : _textframes) {
            if (tf != null) {
                tf.close();
            }
        }
        _textframes.clear();
    }

    void resetSearch() {
        getMainPanel().getCurrentTreePanel().setFoundNodes0(null);
        getMainPanel().getCurrentTreePanel().setFoundNodes1(null);
        getMainPanel().getControlPanel().setSearchFoundCountsOnLabel0(0);
        getMainPanel().getControlPanel().getSearchFoundCountsLabel0().setVisible(false);
        getMainPanel().getControlPanel().getSearchTextField0().setText("");
        getMainPanel().getControlPanel().getSearchResetButton0().setEnabled(false);
        getMainPanel().getControlPanel().getSearchResetButton0().setVisible(false);
        getMainPanel().getControlPanel().setSearchFoundCountsOnLabel1(0);
        getMainPanel().getControlPanel().getSearchFoundCountsLabel1().setVisible(false);
        getMainPanel().getControlPanel().getSearchTextField1().setText("");
        getMainPanel().getControlPanel().getSearchResetButton1().setEnabled(false);
        getMainPanel().getControlPanel().getSearchResetButton1().setVisible(false);
    }

    void setConfiguration(final Configuration configuration) {
        _configuration = configuration;
    }

    void setCurrentDir(final DirectoryPreferences.Category cat, final File current_dir) {
        if (current_dir != null) {
            _current_dirs.put(cat, current_dir);
        }
    }

    void setOptions(final Options options) {
        _options = options;
    }

    /** The currently-open modeless Settings dialog (set/cleared by the launcher), or null. */
    SettingsDialog _open_settings_dialog;

    /** Track the open Settings dialog so a main-window tab switch can re-seed its per-tab controls; auto-clears when
     *  the dialog is disposed. */
    void setOpenSettingsDialog(final SettingsDialog dialog) {
        _open_settings_dialog = dialog;
    }

    /** Re-seed the open (modeless) Settings dialog's per-tab controls (tree style, palette, Time Axis) after a tab
     *  switch, so it reflects the now-current tree. No-op when the dialog isn't open. */
    void refreshOpenSettingsDialog() {
        if ((_open_settings_dialog != null) && _open_settings_dialog.isShowing()) {
            _open_settings_dialog.refreshCurrentTabControls();
        }
    }

    void setSelectedTypeInTypeMenu(final PHYLOGENY_GRAPHICS_TYPE type) {
        setTypeMenuToAllUnselected();
        switch (type) {
            case CIRCULAR:
                _circular_type_cbmi.setSelected(true);
                break;
            case EURO_STYLE:
                _euro_type_cbmi.setSelected(true);
                break;
            case ROUNDED:
                _rounded_type_cbmi.setSelected(true);
                break;
            case RECTANGULAR:
                _rectangular_type_cbmi.setSelected(true);
                break;
            case TRIANGULAR:
                _triangular_type_cbmi.setSelected(true);
                break;
            case UNROOTED:
                _unrooted_type_cbmi.setSelected(true);
                break;
            default:
                throw new IllegalArgumentException("unknown type: " + type);
        }
    }

    final void setSpeciesTree(final Phylogeny species_tree) {
        _species_tree = species_tree;
    }

    void setTypeMenuToAllUnselected() {
        _euro_type_cbmi.setSelected(false);
        _rounded_type_cbmi.setSelected(false);
        _triangular_type_cbmi.setSelected(false);
        _rectangular_type_cbmi.setSelected(false);
        _unrooted_type_cbmi.setSelected(false);
        _circular_type_cbmi.setSelected(false);
    }

    /**
     * The export/save file choosers are created in the {@link MainFrame} constructor,
     * before the look-and-feel is installed, so on macOS they pick up the native Aqua
     * file dialog and keep it. They are also standalone (never part of a window's
     * component tree), so a runtime theme switch does not reach them. Refresh their UI
     * explicitly so they always match the current FlatLaf theme.
     */
    void refreshFileChoosersLookAndFeel() {
        for (final JFileChooser fc : new JFileChooser[] { _writetopdf_filechooser, _writetographics_filechooser,
                _save_filechooser }) {
            if (fc != null) {
                SwingUtilities.updateComponentTreeUI(fc);
            }
        }
    }

    void setDarkMode(final boolean dark) {
        final Configuration.UI ui = dark ? Configuration.UI.FLAT_DARK : Configuration.UI.FLAT_LIGHT;
        getConfiguration().setUi(ui);
        Configuration.saveUiPreference(ui);
        installLookAndFeel(ui);
        // restyle every open window with the new look-and-feel
        for (final Window window : Window.getWindows()) {
            SwingUtilities.updateComponentTreeUI(window);
        }
        // standalone file choosers are not part of any window, so refresh them too
        refreshFileChoosersLookAndFeel();
        // make the tree canvas follow the light/dark theme
        updateTreeCanvasColors(ui);
    }

    void updateTreeCanvasColors(final Configuration.UI ui) {
        if (getMainPanel() == null) {
            return;
        }
        final TreeColorSet colorset = getMainPanel().getTreeColorSet();
        if (colorset == null) {
            return;
        }
        // keep the found-node palette ("Found/Selected Colors" setting) in sync, then set the light/dark scheme --
        // both re-derive the found colors, so they end up correct for the new scheme + choice
        if (getOptions() != null) {
            colorset.setFoundColor(getOptions().getFoundColor());
        }
        // scheme 0 = Dark, scheme 1 = Light (TreeColorSet has only these two)
        colorset.setColorSchema(ui == Configuration.UI.FLAT_DARK ? 0 : 1);
        for (final TreePanel tree_panel : getMainPanel().getTreePanels()) {
            tree_panel.setBackground(colorset.getBackgroundColor());
        }
        if (getMainPanel().getCurrentTreePanel() != null) {
            getMainPanel().getCurrentTreePanel().repaint();
        }
    }


    void typeChanged(final Object o) {
        updateTypeCheckboxes(getOptions(), o);
        updateOptions(getOptions());
        if (getCurrentTreePanel() != null) {
            final PHYLOGENY_GRAPHICS_TYPE previous_type = getCurrentTreePanel().getPhylogenyGraphicsType();
            final PHYLOGENY_GRAPHICS_TYPE new_type = getOptions().getPhylogenyGraphicsType();
            // Triangular draws straight chords to each clade's extreme tips, so it only reads cleanly when the tips are
            // aligned (a cladogram) or the tree is near-ultrametric -- on a ragged phylogram the chords cross the
            // branches. So on a TRANSITION to Triangular, nudge THIS tab to a CLADOGRAM (the user stays free to switch
            // back to P, e.g. for a near-clock tree). Gated on a real transition (previous_type != TRIANGULAR) so
            // re-selecting the already-current Triangular style does NOT re-clobber a deliberate P choice; only when the
            // tree has branch lengths and is currently a phylogram. The nudge does its own re-fit below, so the
            // radial-exit re-fit is skipped when it fires (that fit would be on the stale pre-cladogram layout).
            final boolean nudge_to_cladogram = (previous_type != PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR)
                    && (new_type == PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR) && getCurrentTreePanel().isPhyHasBranchLengths()
                    && getCurrentTreePanel().getControlPanel().isDrawPhylogram();
            if (!nudge_to_cladogram
                    && (((previous_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && (new_type != PHYLOGENY_GRAPHICS_TYPE.UNROOTED))
                            || ((previous_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) && (new_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR))
                            || ((previous_type != PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && (new_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED))
                            || ((previous_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) && (new_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)))) {
                getCurrentTreePanel().getControlPanel().showWhole();
            }
            // The phylogram/cladogram (P/A/C) radios apply in EVERY layout that can honor branch lengths -- since
            // 0.11.7 the CIRCULAR layout renders a real phylogram (isCircularPhylogram) and UNROOTED always has, so the
            // radios are enabled purely on branch-length presence (the old "&& new_type != CIRCULAR" force-disable was a
            // pre-0.11.7 leftover that greyed the radios out in circular even though the paint now responds to them).
            getCurrentTreePanel().getControlPanel().setDrawPhylogramEnabled(getCurrentTreePanel().isPhyHasBranchLengths());
            getCurrentTreePanel().setPhylogenyGraphicsType(getOptions().getPhylogenyGraphicsType());
            // switching TO a radial layout with domains on: auto-enable "Radial Labels" so the domains show (a
            // spoke-riding domain bar needs radial labels; horizontal labels suppress it)
            enableRadialLabelsIfDomainsInRadialLayout();
            // Apply the Triangular->Cladogram nudge computed above. setTreeDisplayType does NOT fire the P/A/C
            // actionPerformed, so the persisted global display-type default is untouched -- a per-tab nudge only.
            if (nudge_to_cladogram) {
                getCurrentTreePanel().getControlPanel().setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM);
                getCurrentTreePanel().getControlPanel().showWhole();
            }
            if ((new_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) || (new_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
                // the showWhole above ran while the panel was still the OLD (rectangular) type, so it laid out a
                // non-square preferred size; now that the panel IS radial, re-fit its SQUARE canvas to the viewport
                // (otherwise the first radial frame draws in the stale rectangular canvas -- off-centre until a re-fit).
                getCurrentTreePanel().getControlPanel().showWhole();
            }
            // Relabel the zoom cluster for the new layout (rectangular <-> radial): the X-/X+ zoom buttons become
            // rotate controls, Y+/Y- become a plain +/- zoom, E greys out, and W becomes the label-direction flip.
            getCurrentTreePanel().getControlPanel().updateZoomButtonsForLayout();
            updateScreenTextAntialias(getMainPanel().getTreePanels());
            if (getCurrentTreePanel().getControlPanel().getDynamicallyHideData() != null) {
                if (new_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
                    getCurrentTreePanel().getControlPanel().getDynamicallyHideData().setEnabled(false);
                } else {
                    getCurrentTreePanel().getControlPanel().getDynamicallyHideData().setEnabled(true);
                }
            }
        }
    }

    /** Flips the node-label direction (horizontal &lt;-&gt; radial) from the radial "L" zoom-cluster button. Drives
     *  the SAME "Radial Labels" checkbox path a Settings toggle uses -- flip the checkbox, then {@code updateOptions}
     *  writes {@code Options._node_label_direction} -- so the option, its menu/dialog state, GuiPreferences
     *  persistence and Reset-to-Defaults all stay consistent (setSelected does not fire the item's action, hence the
     *  explicit updateOptions). Then repaints the current tree. */
    /** When domain architectures are shown in a radial (circular/unrooted) layout, node labels must ride the spoke
     *  for the domain bars to sit cleanly past them (a horizontal label + a spoke-riding bar clash, and the bars are
     *  suppressed under horizontal labels). So the moment domains + a radial layout coincide -- on load, when the user
     *  turns domains on, or when the user switches to a radial layout -- auto-enable "Radial Labels" so the domains are
     *  immediately visible. Sets the option AND syncs the menu checkbox (kept consistent for updateOptions / Reset /
     *  GuiPreferences). No-op if already radial-labelled, not a radial layout, or domains are off. */
    void enableRadialLabelsIfDomainsInRadialLayout() {
        final TreePanel tp = getCurrentTreePanel();
        if (tp == null) {
            return;
        }
        if (tp.isRadialLayout() && (tp.getControlPanel() != null) && tp.getControlPanel().isShowDomainArchitectures()
                && (getOptions().getNodeLabelDirection() != NODE_LABEL_DIRECTION.RADIAL)) {
            getOptions().setNodeLabelDirection(NODE_LABEL_DIRECTION.RADIAL);
            if (_label_direction_cbmi != null) {
                _label_direction_cbmi.setSelected(true);
            }
            tp.repaint();
        }
    }

    void toggleRadialLabelDirection() {
        if (_label_direction_cbmi != null) {
            _label_direction_cbmi.setSelected(!_label_direction_cbmi.isSelected());
            updateOptions(getOptions());
            if (getCurrentTreePanel() != null) {
                getCurrentTreePanel().repaint();
            }
        }
    }

    void updateOptions(final Options options) {
        options.setAbbreviateScientificTaxonNames((_abbreviate_scientific_names != null)
                && _abbreviate_scientific_names.isSelected());
        options.setUseItalicScientificNames((_use_italic_scientific_names_cbmi != null)
                && _use_italic_scientific_names_cbmi.isSelected());
        options.setOutlineFontsInVectorExport((_outline_fonts_in_vector_export_cbmi != null)
                && _outline_fonts_in_vector_export_cbmi.isSelected());
        options.setTransparentExportBackground((_transparent_export_background_cbmi != null)
                && _transparent_export_background_cbmi.isSelected());
        // guarded (not the `!= null && isSelected` pattern): defaults ON, so a null checkbox must keep the
        // default/persisted value rather than clear it (cf. show_tree_name)
        if (_graphics_export_white_background_cbmi != null) {
            options.setGraphicsExportWhiteBackground(_graphics_export_white_background_cbmi.isSelected());
        }
        options.setColorLabelsSameAsParentBranch((_color_labels_same_as_parent_branch != null)
                && _color_labels_same_as_parent_branch.isSelected());
        options.setShowDefaultNodeShapesInternal((_show_default_node_shapes_internal_cbmi != null)
                && _show_default_node_shapes_internal_cbmi.isSelected());
        if ((_internal_labels_above_branch_rbmi != null) && (_internal_labels_right_of_node_rbmi != null)) {
            options.setInternalLabelsAboveBranch(_internal_labels_above_branch_rbmi.isSelected());
        }
        options.setShowDefaultNodeShapesExternal((_show_default_node_shapes_external_cbmi != null)
                && _show_default_node_shapes_external_cbmi.isSelected());
        options.setShowDefaultNodeShapesForMarkedNodes((_show_default_node_shapes_for_marked_cbmi != null)
                && _show_default_node_shapes_for_marked_cbmi.isSelected());
        if ((_non_lined_up_cladograms_rbmi != null) && (_non_lined_up_cladograms_rbmi.isSelected())) {
            options.setCladogramType(CLADOGRAM_TYPE.NON_LINED_UP);
        } else if ((_ext_node_dependent_cladogram_rbmi != null) && (_ext_node_dependent_cladogram_rbmi.isSelected())) {
            options.setCladogramType(CLADOGRAM_TYPE.LINED_UP);
        }
        // The search options (case/words/regex/inverse/properties) are now set directly on Options
        // by the control-panel checkboxes, so they are intentionally not read from menu items here.
        options.setShowScaleGrid((_show_scale_grid_cbmi != null) && _show_scale_grid_cbmi.isSelected());
        options.setShowScaleAxis((_show_scale_axis_cbmi != null) && _show_scale_axis_cbmi.isSelected());
        options.setShowHpdBars((_show_hpd_bars_cbmi != null) && _show_hpd_bars_cbmi.isSelected());
        options.setShowFossilRangeBars((_show_fossil_range_bars_cbmi != null) && _show_fossil_range_bars_cbmi.isSelected());
        options.setShowZebraStripes((_show_zebra_stripes_cbmi != null) && _show_zebra_stripes_cbmi.isSelected());
        options.setBreakLongBranches((_break_long_branches_cbmi != null) && _break_long_branches_cbmi.isSelected());
        options.setShowInternalTaxonomyKey(
                (_show_internal_taxonomy_key_cbmi != null) && _show_internal_taxonomy_key_cbmi.isSelected());
        options.setTipLabelsBelowColumns((_tip_labels_below_columns_cbmi != null)
                && _tip_labels_below_columns_cbmi.isSelected());
        options.setReverseTipOrder((_reverse_tip_order_cbmi != null) && _reverse_tip_order_cbmi.isSelected());
        options.setBoldFoundLabels((_bold_found_labels_cbmi != null) && _bold_found_labels_cbmi.isSelected());
        options.setDimNonMatches((_dim_non_matches_cbmi != null) && _dim_non_matches_cbmi.isSelected());
        options.setPulseFoundNodes((_pulse_found_nodes_cbmi != null) && _pulse_found_nodes_cbmi.isSelected());
        if ((_show_scale_cbmi != null) && _show_scale_cbmi.isEnabled()) {
            options.setShowScale(_show_scale_cbmi.isSelected());
        }
        if (_show_tree_name_cbmi != null) {
            options.setShowTreeName(_show_tree_name_cbmi.isSelected());
        }
        if (_label_direction_cbmi != null) {
            if (_label_direction_cbmi.isSelected()) {
                options.setNodeLabelDirection(NODE_LABEL_DIRECTION.RADIAL);
            } else {
                options.setNodeLabelDirection(NODE_LABEL_DIRECTION.HORIZONTAL);
            }
        }
        options.setShowOverview((_show_overview_cbmi != null) && _show_overview_cbmi.isSelected());
        options.setShowConfidenceStddev((_show_confidence_stddev_cbmi != null)
                && _show_confidence_stddev_cbmi.isSelected());
        options.setShowMadConfidence((_show_mad_confidence_cbmi != null) && _show_mad_confidence_cbmi.isSelected());
        options.setAntialiasPrint((_antialias_print_cbmi != null) && _antialias_print_cbmi.isSelected());
        if ((_use_brackets_for_conf_in_nh_export_cbmi != null)
                && _use_brackets_for_conf_in_nh_export_cbmi.isSelected()) {
            options.setNhConversionSupportValueStyle(NH_CONVERSION_SUPPORT_VALUE_STYLE.IN_SQUARE_BRACKETS);
        } else if ((_use_internal_names_for_conf_in_nh_export_cbmi != null)
                && _use_internal_names_for_conf_in_nh_export_cbmi.isSelected()) {
            options.setNhConversionSupportValueStyle(NH_CONVERSION_SUPPORT_VALUE_STYLE.AS_INTERNAL_NODE_NAMES);
        } else {
            options.setNhConversionSupportValueStyle(NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE);
        }
        options.setPrintBlackAndWhite((_print_black_and_white_cbmi != null)
                && _print_black_and_white_cbmi.isSelected());
        options.setInternalNumberAreConfidenceForNhParsing((_internal_number_are_confidence_for_nh_parsing_cbmi != null)
                && _internal_number_are_confidence_for_nh_parsing_cbmi.isSelected());
        // Taxonomy extraction from node names has no GUI control any more; Options keeps its default
        // (TAXONOMY_EXTRACTION.NO) so the NHX/Nexus parsers read with no extraction.
        options.setReplaceUnderscoresInNhParsing((_replace_underscores_cbmi != null)
                && _replace_underscores_cbmi.isSelected());
        options.setAllowErrorsInDistanceToParent((_allow_errors_in_distance_to_parent_cbmi != null)
                && _allow_errors_in_distance_to_parent_cbmi.isSelected());
        if (_graphics_export_visible_only_cbmi != null) {
            options.setGraphicsExportVisibleOnly(_graphics_export_visible_only_cbmi.isSelected());
        }
        if ((_rectangular_type_cbmi != null) && _rectangular_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
        } else if ((_triangular_type_cbmi != null) && _triangular_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR);
        } else if ((_euro_type_cbmi != null) && _euro_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE);
        } else if ((_rounded_type_cbmi != null) && _rounded_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.ROUNDED);
        } else if ((_unrooted_type_cbmi != null) && _unrooted_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.UNROOTED);
        } else if ((_circular_type_cbmi != null) && _circular_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.CIRCULAR);
        }
        if ((_parse_beast_style_extended_nexus_tags_cbmi != null) && _parse_beast_style_extended_nexus_tags_cbmi.isEnabled()) {
            options.setParseBeastStyleExtendedNexusTags(_parse_beast_style_extended_nexus_tags_cbmi.isSelected());
        }
        if ((_collapsed_with_average_height_cbmi != null) && _collapsed_with_average_height_cbmi.isEnabled()) {
            options.setCollapsedWithAverageHeigh(_collapsed_with_average_height_cbmi.isSelected());
        }
        if ((_show_abbreviated_labels_for_collapsed_nodes_cbmi != null) && _show_abbreviated_labels_for_collapsed_nodes_cbmi.isEnabled()) {
            options.setShowAbbreviatedLabelsForCollapsedNodes(_show_abbreviated_labels_for_collapsed_nodes_cbmi.isSelected());
        }

    }

    private static void setSelected(final JCheckBoxMenuItem item, final boolean selected) {
        if (item != null) {
            item.setSelected(selected);
        }
    }

    /**
     * The inverse of {@link #updateOptions(Options)}: pushes {@code options}'s values OUT to the menu
     * checkboxes/radios so the controls reflect the options (updateOptions reads controls -&gt; options; this
     * writes options -&gt; controls). Uses setSelected (not doClick), so no actions fire. Used by
     * {@link #resetToDefaults()} so the controls match the reset Options instead of clobbering it on the next
     * menu interaction. Must mirror updateOptions's control set exactly -- the reset test flips every control,
     * resets, then runs updateOptions and asserts a default result, so a control missed here fails that test.
     */
    void applyOptionsToMenuStates(final Options options) {
        setSelected(_abbreviate_scientific_names, options.isAbbreviateScientificTaxonNames());
        setSelected(_use_italic_scientific_names_cbmi, options.isUseItalicScientificNames());
        setSelected(_outline_fonts_in_vector_export_cbmi, options.isOutlineFontsInVectorExport());
        setSelected(_transparent_export_background_cbmi, options.isTransparentExportBackground());
        setSelected(_graphics_export_white_background_cbmi, options.isGraphicsExportWhiteBackground());
        setSelected(_color_labels_same_as_parent_branch, options.isColorLabelsSameAsParentBranch());
        setSelected(_show_default_node_shapes_internal_cbmi, options.isShowDefaultNodeShapesInternal());
        setSelected(_show_default_node_shapes_external_cbmi, options.isShowDefaultNodeShapesExternal());
        setSelected(_show_default_node_shapes_for_marked_cbmi, options.isShowDefaultNodeShapesForMarkedNodes());
        setSelected(_show_scale_grid_cbmi, options.isShowScaleGrid());
        setSelected(_show_scale_axis_cbmi, options.isShowScaleAxis());
        setSelected(_show_hpd_bars_cbmi, options.isShowHpdBars());
        setSelected(_show_fossil_range_bars_cbmi, options.isShowFossilRangeBars());
        setSelected(_show_zebra_stripes_cbmi, options.isShowZebraStripes());
        setSelected(_break_long_branches_cbmi, options.isBreakLongBranches());
        setSelected(_show_internal_taxonomy_key_cbmi, options.isShowInternalTaxonomyKey());
        setSelected(_tip_labels_below_columns_cbmi, options.isTipLabelsBelowColumns());
        setSelected(_reverse_tip_order_cbmi, options.isReverseTipOrder());
        setSelected(_bold_found_labels_cbmi, options.isBoldFoundLabels());
        setSelected(_dim_non_matches_cbmi, options.isDimNonMatches());
        setSelected(_pulse_found_nodes_cbmi, options.isPulseFoundNodes());
        setSelected(_show_scale_cbmi, options.isShowScale());
        setSelected(_show_tree_name_cbmi, options.isShowTreeName());
        setSelected(_show_overview_cbmi, options.isShowOverview());
        setSelected(_show_confidence_stddev_cbmi, options.isShowConfidenceStddev());
        setSelected(_show_mad_confidence_cbmi, options.isShowMadConfidence());
        setSelected(_antialias_print_cbmi, options.isAntialiasPrint());
        setSelected(_print_black_and_white_cbmi, options.isPrintBlackAndWhite());
        setSelected(_internal_number_are_confidence_for_nh_parsing_cbmi,
                options.isInternalNumberAreConfidenceForNhParsing());
        setSelected(_replace_underscores_cbmi, options.isReplaceUnderscoresInNhParsing());
        setSelected(_allow_errors_in_distance_to_parent_cbmi, options.isAllowErrorsInDistanceToParent());
        setSelected(_graphics_export_visible_only_cbmi, options.isGraphicsExportVisibleOnly());
        setSelected(_parse_beast_style_extended_nexus_tags_cbmi, options.isParseBeastStyleExtendedNexusTags());
        setSelected(_collapsed_with_average_height_cbmi, options.isCollapsedWithAverageHeigh());
        setSelected(_show_abbreviated_labels_for_collapsed_nodes_cbmi,
                options.isShowAbbreviatedLabelsForCollapsedNodes());
        setSelected(_label_direction_cbmi, options.getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL);
        // radio groups / tri-states (not simple checkboxes)
        if (_internal_labels_above_branch_rbmi != null) {
            _internal_labels_above_branch_rbmi.setSelected(options.isInternalLabelsAboveBranch());
        }
        if (_internal_labels_right_of_node_rbmi != null) {
            _internal_labels_right_of_node_rbmi.setSelected(!options.isInternalLabelsAboveBranch());
        }
        final boolean lined_up = options.getCladogramType() != CLADOGRAM_TYPE.NON_LINED_UP;
        if (_ext_node_dependent_cladogram_rbmi != null) {
            _ext_node_dependent_cladogram_rbmi.setSelected(lined_up);
        }
        if (_non_lined_up_cladograms_rbmi != null) {
            _non_lined_up_cladograms_rbmi.setSelected(!lined_up);
        }
        final NH_CONVERSION_SUPPORT_VALUE_STYLE nh = options.getNhConversionSupportValueStyle();
        setSelected(_use_brackets_for_conf_in_nh_export_cbmi,
                nh == NH_CONVERSION_SUPPORT_VALUE_STYLE.IN_SQUARE_BRACKETS);
        setSelected(_use_internal_names_for_conf_in_nh_export_cbmi,
                nh == NH_CONVERSION_SUPPORT_VALUE_STYLE.AS_INTERNAL_NODE_NAMES);
        // setSelectedTypeInTypeMenu -> setTypeMenuToAllUnselected dereferences the 6 type items unguarded; match the
        // null-tolerance of the rest of this method (and updateOptions) so a frame without a Type menu can't NPE
        if (_rectangular_type_cbmi != null) {
            setSelectedTypeInTypeMenu(options.getPhylogenyGraphicsType());
        }
    }

    /**
     * Reset to Defaults: returns all display/UI settings to the built-in defaults, the theme to the shipped
     * default (light), the search options and tree style to default, forgets the persisted settings file (so the
     * reset survives a restart), and turns off property-based "Color by" (default palette) on every open tree.
     * Does NOT reload the trees or strip manually applied branch/clade colors (those are tree data / Undo territory).
     */
    void resetToDefaults() {
        // 1. settings: reset the LIVE Options in place against the frame's own configuration (== the baseline
        //    createInstance produced at launch) -- each TreePanel caches this same reference, so an in-place reset
        //    propagates everywhere without a swap
        getOptions().resetToDefaults();
        // 2. push the defaults out to the menu controls, so a later updateOptions() reads defaults, not stale state
        applyOptionsToMenuStates(getOptions());
        // 3. forget the persisted settings so the reset survives the next launch
        new GuiPreferences().deleteSettingsFile();
        // 4. theme -> shipped default (FlatLaf light); re-inits the L&F and restyles every open window live
        setDarkMode(false);
        if (getMainPanel() != null) {
            final Options.PHYLOGENY_GRAPHICS_TYPE type = getOptions().getPhylogenyGraphicsType();
            // 5. per-tab: turn property "Color by" OFF (default palette) AND apply the reset tree style to the live
            //    tree -- the TreePanel caches its own graphics type, so resetting Options alone would leave e.g. a
            //    Circular tree still drawn circular
            for (final TreePanel tp : getMainPanel().getTreePanels()) {
                tp.resetColorStateToDefaults(); // also turns ancestral-state pies OFF (per-tab)
                tp.setPhylogenyGraphicsType(type);
                tp.resetTimeAxisToAutoDerive(); // per-tab: drop any Time-Axis override -> back to auto-derive
                tp.resetNextstrainBranchModeToDefault(); // per-tab: back to the TIME branch-length view (Auspice trees)
            }
            final ControlPanel cp = getMainPanel().getControlPanel();
            if (cp != null) {
                cp.setColorByPropertySelectionToNone();
                cp.setSizeByPropertySelectionToNone();
                cp.setAncestralPieSelectionToNone();
                cp.setBranchLengthsSelectionToTime();
                // re-seed the always-visible control-panel controls (theme radios + search checkboxes) that hold
                // their own state -- else they stay stale and the search checkboxes clobber the reset on next click
                cp.resyncFromOptions();
                cp.updateZoomButtonsForLayout(); // layout reset to RECTANGULAR/ROOT_LEFT -> H reverts to W, radial labels revert
            }
            final TreePanel current = getMainPanel().getCurrentTreePanel();
            if ((current != null) && (cp != null)) {
                // phylogram is available for the default (rectangular) style when the tree has branch lengths; also
                // reset the P/A/C display type to the (now-default) preference for a branch-length tree, else to
                // cladogram -- mirroring the load-time auto-detect so Reset returns the tree shape to a fresh install
                if (current.isPhyHasBranchLengths()) {
                    cp.setDrawPhylogramEnabled(true);
                    // use the SAME preference-aware policy as the load-time auto-detect, so a SPARSE-branch-length
                    // tree resets to a cladogram (not a degenerate phylogram) exactly as a fresh load would
                    cp.setTreeDisplayType(AptxUtil.preferredDisplayTypeForBranchLengthTree(
                            AptxUtil.isHasAtLeast50PercentBranchLengthLargerThanZero(current.getPhylogeny()),
                            getOptions().getPhylogenyDisplayType()));
                }
                else {
                    cp.setDrawPhylogramEnabled(false);
                    cp.setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM);
                }
            }
            // 6. apply the reset options to the live view (base font + relayout + re-fit + repaint), mirroring
            //    chooseFont()/typeChanged()
            getMainPanel().getTreeFontSet().setBaseFont(getOptions().getBaseFont());
            if (cp != null) {
                cp.displayedPhylogenyMightHaveChanged(true);
                cp.showWhole(); // re-fit: the tree style may have changed (e.g. Circular -> Rectangular)
            }
            if (current != null) {
                current.resetPreferredSize();
                current.updateOvSizes();
            }
        }
        repaint();
    }

    void updateTypeCheckboxes(final Options options, final Object o) {
        setTypeMenuToAllUnselected();
        ((JCheckBoxMenuItem) o).setSelected(true);
    }

    void viewAsNexus() {
        if ((_mainpanel.getCurrentPhylogeny() != null) && !_mainpanel.getCurrentPhylogeny().isEmpty()) {
            String title = "Nexus";
            if (!ForesterUtil.isEmpty(_mainpanel.getCurrentPhylogeny().getName())) {
                title = "\"" + getMainPanel().getCurrentPhylogeny().getName() + "\" in " + title;
            }
            showTextFrame(_mainpanel.getCurrentPhylogeny().toNexus(getOptions().getNhConversionSupportValueStyle()),
                    title);
        }
    }

    void viewAsNH() {
        if ((_mainpanel.getCurrentPhylogeny() != null) && !_mainpanel.getCurrentPhylogeny().isEmpty()) {
            String title = "New Hampshire";
            if (!ForesterUtil.isEmpty(_mainpanel.getCurrentPhylogeny().getName())) {
                title = "\"" + getMainPanel().getCurrentPhylogeny().getName() + "\" in " + title;
            }
            showTextFrame(_mainpanel.getCurrentPhylogeny().toNewHampshire(getOptions()
                            .getNhConversionSupportValueStyle()),
                    title);
        }
    }

    void viewAsXML() {
        if ((_mainpanel.getCurrentPhylogeny() != null) && !_mainpanel.getCurrentPhylogeny().isEmpty()) {
            String title = "phyloXML";
            if (!ForesterUtil.isEmpty(_mainpanel.getCurrentPhylogeny().getName())) {
                title = "\"" + getMainPanel().getCurrentPhylogeny().getName() + "\" in " + title;
            }
            showTextFrame(_mainpanel.getCurrentPhylogeny().toPhyloXML(0), title);
        }
    }

    /**
     * Display the about box.
     */
    void about() {
        JOptionPane.showMessageDialog(null, buildAboutText(), AptxConstants.PRG_NAME, JOptionPane.PLAIN_MESSAGE,
                aboutLogo());
    }

    /** The About-box text (program info, environment, and the main reference). Extracted so it is testable, and so
     *  its content stays independent of the dialog chrome. */
    static String buildAboutText() {
        final StringBuilder about = new StringBuilder("Archaeopteryx\nVersion " + AptxConstants.VERSION + "\n");
        about.append("Copyright (C) 2026 Christian M Zmasek\n");
        about.append("All Rights Reserved\n");
        about.append("License: GNU General Public License version 3 (GPL3)\n");
        about.append("Last modified: " + AptxConstants.PRG_DATE + "\n");
        about.append("Based on: ").append(ForesterUtil.getForesterLibraryInformation()).append("\n");
        about.append("phyloXML version : " + ForesterConstants.PHYLO_XML_VERSION + "\n");
        about.append("phyloXML location: " + ForesterConstants.PHYLO_XML_LOCATION + "\n");
        if (!ForesterUtil.isEmpty(ForesterUtil.JAVA_VERSION) && !ForesterUtil.isEmpty(ForesterUtil.JAVA_VENDOR)) {
            about.append("[your Java version: ").append(ForesterUtil.JAVA_VERSION).append(" ").append(ForesterUtil.JAVA_VENDOR).append("]\n");
        }
        if (!ForesterUtil.isEmpty(ForesterUtil.OS_NAME) && !ForesterUtil.isEmpty(ForesterUtil.OS_ARCH)
                && !ForesterUtil.isEmpty(ForesterUtil.OS_VERSION)) {
            about.append("[your OS: ").append(ForesterUtil.OS_NAME).append(" ").append(ForesterUtil.OS_ARCH).append(" ").append(ForesterUtil.OS_VERSION).append("]\n");
        }
        final Runtime rt = java.lang.Runtime.getRuntime();
        final long free_memory = rt.freeMemory() / 1000000;
        final long total_memory = rt.totalMemory() / 1000000;
        about.append("[free memory: ").append(free_memory).append("MB, total memory: ").append(total_memory).append("MB]\n");
        about.append("[locale: ").append(Locale.getDefault()).append("]\n");
        about.append("References:\n");
        about.append(AptxConstants.PHYLOXML_REFERENCE_SHORT + "\n");
        about.append("Comments: " + AptxConstants.AUTHOR_EMAIL);
        return about.toString();
    }

    /** The bundled Archaeopteryx logo for the About box, or {@code null} if the image is missing (best-effort: the
     *  dialog then simply shows without it). Serves a HiDPI-aware image: the icon stays ~128pt but AWT picks the
     *  @2x (256px) variant on a Retina/2x display, so it stays crisp. */
    static javax.swing.Icon aboutLogo() {
        try {
            final java.awt.image.BufferedImage base = readImage("/resources/images/archaeopteryx-logo.png");
            if (base == null) {
                return null;
            }
            final java.awt.image.BufferedImage x2 = readImage("/resources/images/archaeopteryx-logo@2x.png");
            final java.awt.Image image = (x2 != null)
                    ? new java.awt.image.BaseMultiResolutionImage(base, x2) // base = logical (128pt) size; x2 on HiDPI
                    : base;
            return new javax.swing.ImageIcon(image);
        }
        catch (final Exception e) {
            return null; // best-effort: the About box just shows without a logo
        }
    }

    private static java.awt.image.BufferedImage readImage(final String resource) throws java.io.IOException {
        try (final java.io.InputStream in = MainFrame.class.getResourceAsStream(resource)) {
            return (in == null) ? null : javax.imageio.ImageIO.read(in);
        }
    }

    /** Read-only, selectable, scrolling view of {@link AlgorithmReferences} -- extracted so it is testable. */
    static JScrollPane buildReferencesView() {
        final JTextArea ta = new JTextArea(AlgorithmReferences.asText());
        ta.setEditable(false);
        ta.setLineWrap(true);
        ta.setWrapStyleWord(true);
        ta.setCaretPosition(0);
        ta.setBorder(BorderFactory.createEmptyBorder(8, 10, 8, 10));
        final JScrollPane scroller = new JScrollPane(ta, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        // A comfortable, conservative size (well within any usable screen); the dialog is resizable and the
        // list scrolls, so citations are always reachable and copyable.
        scroller.setPreferredSize(new Dimension(620,
                Math.min(520, ta.getPreferredSize().height + 24)));
        return scroller;
    }

    void showReferences() {
        final JOptionPane pane = new JOptionPane(buildReferencesView(), JOptionPane.PLAIN_MESSAGE,
                JOptionPane.DEFAULT_OPTION);
        final JDialog dialog = pane.createDialog(this, "References");
        dialog.setResizable(true);
        dialog.setVisible(true);
        dialog.dispose();
    }

    static JMenu createMenu(final String title, final Configuration conf) {
        return new JMenu(title);
    }

    static void cycleOverview(final Options op, final TreePanel tree_panel) {
        switch (op.getOvPlacement()) {
            case LOWER_LEFT:
                op.setOvPlacement(Options.OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT);
                break;
            case LOWER_RIGHT:
                op.setOvPlacement(Options.OVERVIEW_PLACEMENT_TYPE.LOWER_LEFT);
                break;
            case UPPER_LEFT:
                op.setOvPlacement(Options.OVERVIEW_PLACEMENT_TYPE.UPPER_RIGHT);
                break;
            case UPPER_RIGHT:
                op.setOvPlacement(Options.OVERVIEW_PLACEMENT_TYPE.LOWER_RIGHT);
                break;
            default:
                throw new RuntimeException("unknown placement: " + op.getOvPlacement());
        }
        if (tree_panel != null) {
            tree_panel.updateOvSettings();
        }
    }

    static void exceptionOccuredDuringSaveAs(final Exception e, final TreePanel tp, final Component comp) {
        try {
            tp.setArrowCursor();
        } catch (final Exception ex) {
            // Do nothing.
        }
        JOptionPane.showMessageDialog(comp, "Exception" + e, "Error during File|SaveAs", JOptionPane.ERROR_MESSAGE);
    }



    static void exportPhylogenyToPdf(final String file_name,
                                     final Options opts,
                                     final TreePanel tp,
                                     final Component comp) {

        String pdf_written_to = "";
        boolean error = false;
        try {
            if (opts.isPrintUsingActualSize()) {
                pdf_written_to = PdfExporter.writePhylogenyToPdf(file_name, tp, tp.getWidth(), tp.getHeight(),
                        opts.isGraphicsExportWhiteBackground());
            } else {
                // Never false.
            }
        } catch (final IOException e) {
            error = true;
            JOptionPane.showMessageDialog(comp, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
        }
        if (!error) {
            if (!ForesterUtil.isEmpty(pdf_written_to)) {
                JOptionPane.showMessageDialog(comp,
                        "Wrote PDF to: " + pdf_written_to,
                        "Information",
                        JOptionPane.INFORMATION_MESSAGE);
            } else {
                JOptionPane.showMessageDialog(comp,
                        "There was an unknown problem when attempting to write to PDF file: \""
                                + file_name + "\"",
                        "Error",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
        if (!opts.isPrintUsingActualSize()) {
            tp.getControlPanel().showWhole();
        }
    }

    static void updateScreenTextAntialias(final List<TreePanel> treepanels) {
        for (final TreePanel tree_panel : treepanels) {
            tree_panel.setTextAntialias();
        }
    }

    static boolean writeAsNewHampshire(final TreePanel tp, final Options op, boolean exception, final File file) {
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toNewHampshire(tp.getPhylogeny(), true, op.getNhConversionSupportValueStyle(), file);
        } catch (final Exception e) {
            exception = true;
            exceptionOccuredDuringSaveAs(e, tp, tp);
        }
        return exception;
    }

    static boolean writeAsNexus(final TreePanel tp, final Options op, boolean exception, final File file) {
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toNexus(file, tp.getPhylogeny(), op.getNhConversionSupportValueStyle());
        } catch (final Exception e) {
            exception = true;
            exceptionOccuredDuringSaveAs(e, tp, tp);
        }
        return exception;
    }

    static boolean writeAsPhyloXml(final TreePanel tp, final Options op, boolean exception, final File file) {
        try {
            tp.syncTimeAxisConfigToTree(); // embed the per-tree Time-Axis config (only if it deviates from auto-derive)
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML(file, tp.getPhylogeny(), 0);
        } catch (final Exception e) {
            exception = true;
            exceptionOccuredDuringSaveAs(e, tp, tp);
        }
        return exception;
    }

    static void writePhylogenyToGraphicsFile(final String file_name,
                                             final GraphicsExportType type,
                                             final MainPanel mp,
                                             final Component comp,
                                             final Container contentpane) {
        mp.getCurrentTreePanel().calcParametersForPainting(mp.getCurrentTreePanel().getWidth(),
                mp.getCurrentTreePanel().getHeight());
        String file_written_to = "";
        boolean error = false;
        try {
            file_written_to = AptxUtil.writePhylogenyToGraphicsFile(file_name,
                    mp.getCurrentTreePanel().getWidth(),
                    mp.getCurrentTreePanel().getHeight(),
                    mp.getCurrentTreePanel(),
                    mp.getControlPanel(),
                    type,
                    mp.getOptions());
        } catch (final IOException e) {
            error = true;
            JOptionPane.showMessageDialog(comp, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
        }
        if (!error) {
            if ((file_written_to != null) && (file_written_to.length() > 0)) {
                JOptionPane.showMessageDialog(comp,
                        "Wrote image to: " + file_written_to,
                        "Graphics Export",
                        JOptionPane.INFORMATION_MESSAGE);
            } else {
                JOptionPane.showMessageDialog(comp,
                        "There was an unknown problem when attempting to write to an image file: \""
                                + file_name + "\"",
                        "Error",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
        contentpane.repaint();
    }

    /** Outcome of a Copy-Image-to-Clipboard attempt; kept UI-free so the decision logic is unit-testable. */
    enum ClipboardCopyResult {
        COPIED, NO_TREE, CLIPBOARD_UNAVAILABLE
    }

    /** Testable core of File -> Copy Image to Clipboard: renders the current tree onto {@code clipboard} and
     *  reports what happened. No dialogs, no reference to the system clipboard, so a test can drive it with a
     *  private clipboard. */
    ClipboardCopyResult copyImageToClipboard(final Clipboard clipboard) {
        final TreePanel tp = getCurrentTreePanel();
        if ((tp == null) || (tp.getPhylogeny() == null) || tp.getPhylogeny().isEmpty()) {
            return ClipboardCopyResult.NO_TREE;
        }
        try {
            return AptxUtil.copyPhylogenyImageToClipboard(tp, getOptions(), clipboard, null)
                    ? ClipboardCopyResult.COPIED : ClipboardCopyResult.NO_TREE;
        } catch (final RuntimeException | OutOfMemoryError ex) {
            // best-effort: clipboard busy/owned (IllegalStateException), headless/denied (Headless/Security), an
            // unexpected rendering failure, or an OutOfMemoryError -- rendering now uses the raster-export scale
            // (default 4x), so a large tree can request a big bitmap (bounded by renderPhylogenyToImage's 100 MP
            // cap) and OOM on a small heap. A single failed allocation is recoverable once dropped, so catch it too
            // rather than let it escape onto the EDT; the oversized image is unreferenced on return.
            return ClipboardCopyResult.CLIPBOARD_UNAVAILABLE;
        }
    }

    /** File -> Copy Image to Clipboard: render the current tree as an (opaque) image and put it on the system
     *  clipboard, ready to paste into a document or slide. Silent on success (the paste is the feedback);
     *  warns only when there is no tree or the clipboard is unavailable. */
    void copyImageToClipboard() {
        final Clipboard clipboard;
        try {
            clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
        } catch (final RuntimeException ex) {
            // getSystemClipboard can throw HeadlessException / SecurityException
            JOptionPane.showMessageDialog(this, "The clipboard is currently unavailable; please try again.",
                    "Copy Image", JOptionPane.WARNING_MESSAGE);
            return;
        }
        switch (copyImageToClipboard(clipboard)) {
            case NO_TREE:
                JOptionPane.showMessageDialog(this, "There is no tree to copy.", "Copy Image",
                        JOptionPane.WARNING_MESSAGE);
                break;
            case CLIPBOARD_UNAVAILABLE:
                JOptionPane.showMessageDialog(this,
                        "The image could not be copied; the clipboard may be busy, or the figure may be too large "
                                + "to copy at the current export resolution (lower it in Options). Please try again.",
                        "Copy Image", JOptionPane.WARNING_MESSAGE);
                break;
            case COPIED:
                break; // silent success
        }
    }

    /** File -> Export Sequences (FASTA): write the current tree's tip molecular sequences as FASTA. */
    void exportSequencesAsFasta() {
        final Phylogeny phy = currentPhylogenyForExport();
        if (phy != null) {
            final ExportScope scope = chooseExportTips(phy);
            if (scope == null) {
                return; // cancelled
            }
            final String fasta = NodeDataExporter.toFasta(scope.tips);
            final int records = NodeDataExporter.fastaRecordCount(fasta);
            writeDataExportToFile(fasta, "molecular sequences (FASTA)", records, scope.description,
                    suggestedExportName(phy, ".fasta", scope.restricted ? records : -1), null);
        }
    }

    /** File -> Export Node Data (TSV): write the current tree's tip data as a tab-separated table. */
    void exportNodeDataAsTsv() {
        final Phylogeny phy = currentPhylogenyForExport();
        if (phy != null) {
            final ExportScope scope = chooseExportTips(phy);
            if (scope == null) {
                return; // cancelled
            }
            final String note = NodeDataExporter.tipNamesFormUniqueKey(scope.tips) ? null
                    : "Some tip names are blank or duplicated, so a \"node_id\" column was added as the unique key.";
            // one row per tip, so the row count is the size of the chosen scope
            final int rows = scope.tips.size();
            writeDataExportToFile(NodeDataExporter.toNodeDataTsv(scope.tips), "node-data rows (TSV)", rows,
                    scope.description, suggestedExportName(phy, ".tsv", scope.restricted ? rows : -1), note);
        }
    }

    /**
     * File -> Import Annotations (CSV/TSV): read a tip-keyed CSV or TSV table and write its columns onto the
     * matching external nodes. The user picks the key column and the tip attribute to match it against (tip name /
     * sequence accession / taxonomy id / taxonomy scientific name), with a live dry-run preview of the join before
     * committing. Recognized columns fill the taxonomy/sequence fields; any other column becomes a node property
     * usable for column coloring / annotation columns.
     */
    void importAnnotations() {
        final Phylogeny phy = currentPhylogenyForExport();
        if (phy == null) {
            return;
        }
        final JFileChooser fc = new JFileChooser();
        fc.setMultiSelectionEnabled(false);
        fc.setDialogTitle("Import Annotations (CSV/TSV)");
        fc.setFileFilter(new FileNameExtensionFilter("Annotation tables (*.csv, *.tsv, *.txt)", "csv", "tsv", "txt"));
        if (getCurrentDir(DirectoryPreferences.Category.OPEN) != null) {
            fc.setCurrentDirectory(getCurrentDir(DirectoryPreferences.Category.OPEN));
        }
        if (fc.showOpenDialog(this) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        final File file = fc.getSelectedFile();
        if (file == null) {
            return;
        }
        setCurrentDir(DirectoryPreferences.Category.OPEN, fc.getCurrentDirectory());
        final String text;
        try {
            text = Files.readString(file.toPath());
        }
        catch (final IOException e) {
            JOptionPane.showMessageDialog(this, "Failed to read " + file + ":\n" + e.getMessage(), "Read Failed",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        importAnnotationsFromText(phy, text, file.getName(), file.getAbsolutePath(), false);
    }

    /**
     * File -> Import Annotations from URL: fetch a CSV/TSV from a URL (e.g. a Google Sheet "published to the web as
     * CSV") and run the same import dialog. Reads the whole body (so a quoted field may span lines).
     */
    void importAnnotationsFromUrl() {
        final Phylogeny phy = currentPhylogenyForExport();
        if (phy == null) {
            return;
        }
        final String url = JOptionPane.showInputDialog(this,
                "Enter a CSV/TSV URL (e.g. a Google Sheet published to the web as CSV):",
                "Import Annotations from URL", JOptionPane.PLAIN_MESSAGE);
        if (ForesterUtil.isEmpty(url)) {
            return;
        }
        final String text;
        try {
            text = ForesterUtil.readUrlToString(url.trim()); // full body (preserves newlines inside quoted fields)
        }
        catch (final Exception e) {
            JOptionPane.showMessageDialog(this, "Failed to fetch " + url.trim() + ":\n" + e.getMessage(),
                    "Fetch Failed", JOptionPane.ERROR_MESSAGE);
            return;
        }
        importAnnotationsFromText(phy, text, url.trim(), url.trim(), true);
    }

    /**
     * Shared back-end for the file and URL import sources: validate, show the config dialog, apply, report, and
     * remember the import (source {@code reimport_locator} + column mapping) on the tree for one-click Re-import.
     */
    private void importAnnotationsFromText(final Phylogeny phy, final String text, final String display_name,
            final String reimport_locator, final boolean is_url) {
        try {
            NodeDataImporter.parseTable(text); // validate up front so a bad source shows an error, not an empty dialog
        }
        catch (final IllegalArgumentException e) {
            JOptionPane.showMessageDialog(this, "This does not look like a CSV/TSV table:\n" + e.getMessage(),
                    "Import Failed", JOptionPane.ERROR_MESSAGE);
            return;
        }
        final ImportChoice choice = showImportChoice(phy, text);
        if (choice == null) {
            return;
        }
        final NodeDataImporter.ImportResult result = importAnnotationsAndRefit(phy, choice._table, choice._key_col,
                choice._match_by, choice._plan, display_name);
        final TreePanel tp = getCurrentTreePanel();
        // Only remember/persist the profile when the import ACTUALLY annotated tips -- mirrors the undo + setEdited
        // gating in importAnnotationsAndRefit. A confirmed 0-match import (wrong key/attribute) must not clobber a good
        // prior profile (un-undoably, since no checkpoint was taken) nor mark the tree dirty for a no-op change.
        if ((tp != null) && (result.getTipsAnnotated() > 0)) {
            final NodeDataImporter.ImportProfile profile = NodeDataImporter.ImportProfile.from(choice._table,
                    choice._key_col, choice._match_by, choice._plan, reimport_locator, is_url);
            tp.setLastImportProfile(profile); // enable one-click Re-import this session
            NodeDataImporter.writeProfileToTree(phy, profile); // ...and persist it with the tree (survives save/reload)
        }
        final int type = result.getWarnings().isEmpty() ? JOptionPane.INFORMATION_MESSAGE : JOptionPane.WARNING_MESSAGE;
        JOptionPane.showMessageDialog(this, "Imported from " + display_name + ":\n\n" + result.summary(),
                "Import Complete", type);
    }

    /**
     * File -> Re-import Annotations: re-fetch the tree's last import source (file or URL) and re-apply the same key
     * column, match attribute, and column include/rename mapping with one click -- so you can edit your sheet/file
     * and pull the changes in without walking the dialog again. The mapping is keyed by header name, so it survives
     * the source gaining or reordering rows/columns.
     */
    void reimportAnnotations() {
        final Phylogeny phy = currentPhylogenyForExport();
        if (phy == null) {
            return;
        }
        final TreePanel tp = getCurrentTreePanel();
        final NodeDataImporter.ImportProfile profile = (tp != null) ? tp.getLastImportProfile() : null;
        if (profile == null) {
            JOptionPane.showMessageDialog(this,
                    "No annotation import to repeat for this tree.\nUse \"Import Annotations (CSV/TSV)…\" first.",
                    "Nothing to Re-import", JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        final NodeDataImporter.ImportResult result;
        try {
            result = reimportAnnotationsAndRefit(phy, profile);
        }
        catch (final IOException e) {
            JOptionPane.showMessageDialog(this, "Failed to re-read " + profile.getSource() + ":\n" + e.getMessage(),
                    "Re-import Failed", JOptionPane.ERROR_MESSAGE);
            return;
        }
        catch (final IllegalArgumentException e) {
            JOptionPane.showMessageDialog(this, "The source is no longer a valid CSV/TSV table:\n" + e.getMessage(),
                    "Re-import Failed", JOptionPane.ERROR_MESSAGE);
            return;
        }
        final String name = importSourceDisplayName(profile);
        if (result.getTipsAnnotated() == 0) {
            // a one-click re-import has no dialog/preview, so a 0-match (e.g. the edited source broke the key column)
            // must not look like success -- the tree was left unchanged
            JOptionPane.showMessageDialog(this,
                    "Re-imported from " + name + ", but no rows matched a tip — the tree was not changed.\n\n"
                            + "Check the key column / match attribute against the edited source (use \"Import "
                            + "Annotations (CSV/TSV)…\" to reconfigure).\n\n" + result.summary(),
                    "Re-import: Nothing Matched", JOptionPane.WARNING_MESSAGE);
            return;
        }
        final int type = result.getWarnings().isEmpty() ? JOptionPane.INFORMATION_MESSAGE : JOptionPane.WARNING_MESSAGE;
        JOptionPane.showMessageDialog(this, "Re-imported from " + name + ":\n\n" + result.summary(),
                "Re-import Complete", type);
    }

    /**
     * UI-free re-import (so it is unit-testable without the dialogs): re-fetch the profile's source (file or URL),
     * parse it with the remembered delimiter, resolve the profile's key column + column mapping against the fresh
     * table (by header name), and apply.
     *
     * @throws IOException if the source cannot be re-read; {@code IllegalArgumentException} if it no longer parses
     */
    NodeDataImporter.ImportResult reimportAnnotationsAndRefit(final Phylogeny phy,
            final NodeDataImporter.ImportProfile profile) throws IOException {
        final String text = profile.isUrl() ? ForesterUtil.readUrlToString(profile.getSource())
                : Files.readString(java.nio.file.Path.of(profile.getSource()));
        final NodeDataImporter.Table table = NodeDataImporter.parseTable(text, profile.getDelimiter());
        return importAnnotationsAndRefit(phy, table, profile.keyColumn(table), profile.getMatchBy(),
                profile.columnPlan(table), importSourceDisplayName(profile));
    }

    /** A short display name for a profile's source: the URL as-is, or a file's base name. */
    private static String importSourceDisplayName(final NodeDataImporter.ImportProfile profile) {
        return profile.isUrl() ? profile.getSource() : new File(profile.getSource()).getName();
    }

    /**
     * View -> Clustergram: one click to build the signature vertical-dendrogram-over-a-heat-map figure. Forces a
     * rectangular, root-at-top, tip-aligned layout with the tip labels below the columns, and auto-builds the
     * annotation columns from the tree's own per-tip data (every numeric field -> one shared-scale heat-map MATRIX;
     * every categorical field -> a color strip). Display-only: Undo is N/A, and every setting it flips is already
     * covered by Reset to Defaults (no new Options field).
     */
    void applyClustergramPreset() {
        final TreePanel tp = getCurrentTreePanel();
        if (tp == null) {
            return;
        }
        final Phylogeny phy = tp.getPhylogeny();
        if ((phy == null) || phy.isEmpty()) {
            return;
        }
        final java.util.List<AnnotationColumns.ColumnSpec> specs = clustergramColumnSpecs(phy);
        if (specs.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "This tree has no per-tip data to put in the heat map.\n\nUse File → Import Annotations (CSV/TSV) "
                            + "to add columns, then run Clustergram again.",
                    "Clustergram", JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        tp.setAnnotationColumns(specs);
        // a rectangular, root-at-top, tip-aligned tree with the sample labels below the columns = the clustergram
        getOptions().setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
        tp.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
        getOptions().setTreeOrientation(Options.TREE_ORIENTATION.ROOT_TOP);
        getOptions().setTipLabelsBelowColumns(true);
        getOptions().setTipLabelDirection(Options.TIP_LABEL_DIRECTION.AUTO); // upright short names, tilt as density rises
        applyOptionsToMenuStates(getOptions()); // reflect the labels-below toggle in the menu / Settings dialog
        tp.getControlPanel().setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM); // align tips -> clean grid
        // NB: display-only -- do NOT setEdited(true): nothing here is saved to the tree file, and setEdited would both
        // pop a spurious save prompt and clear the redo stack (the undo safety net). Sibling display toggles don't either.
        final ControlPanel cp = tp.getControlPanel();
        cp.updateZoomButtonsForLayout();
        cp.displayedPhylogenyMightHaveChanged(true); // recompute the label/column extents for the swapped layout
        cp.showWhole(); // re-fit for the vertical extent (a plain repaint would leave the old scroll extent)
        tp.repaint();
    }

    /**
     * The annotation columns a Clustergram builds from a tree's data: every numeric field as a shared-scale heat-map
     * MATRIX (nearest the tips -> the dendrogram sits directly on the grid), then every categorical field as a color
     * strip. Empty when the tree carries no color-able per-tip properties. Pure + testable.
     */
    static java.util.List<AnnotationColumns.ColumnSpec> clustergramColumnSpecs(final Phylogeny phy) {
        final java.util.List<String> colorable = PropertyColorScheme.colorableRefs(phy); // scan the tree ONCE
        final java.util.List<String> numeric = PropertyColorScheme.numericRefs(phy, colorable);
        final java.util.List<AnnotationColumns.ColumnSpec> specs = new java.util.ArrayList<>();
        for (final String ref : numeric) {
            specs.add(new AnnotationColumns.ColumnSpec(ref, AnnotationColumns.Type.MATRIX));
        }
        for (final String ref : colorable) {
            if (!numeric.contains(ref)) {
                specs.add(new AnnotationColumns.ColumnSpec(ref, AnnotationColumns.Type.COLOR_STRIP));
            }
        }
        return specs;
    }

    /** The user's choices from the import config dialog: which table (a given delimiter), key column, match attribute, and column plan. */
    private static final class ImportChoice {

        private final NodeDataImporter.Table      _table;
        private final int                         _key_col;
        private final NodeDataImporter.MatchBy     _match_by;
        private final NodeDataImporter.ColumnPlan _plan;

        private ImportChoice(final NodeDataImporter.Table table, final int key_col,
                final NodeDataImporter.MatchBy match_by, final NodeDataImporter.ColumnPlan plan) {
            _table = table;
            _key_col = key_col;
            _match_by = match_by;
            _plan = plan;
        }
    }

    /**
     * The import configuration dialog: a delimiter override (which re-parses), a preview, the key column + tip
     * attribute to match, a per-column include/rename panel, and a live dry-run summary. Returns the user's choices,
     * or {@code null} on cancel.
     */
    private ImportChoice showImportChoice(final Phylogeny phy, final String text) {
        final javax.swing.JDialog dialog = new javax.swing.JDialog(this, "Import Annotations", true);
        dialog.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE); // 'X' disposes (not just hides)
        final ImportChoice[] result = { null };
        final NodeDataImporter.Table[] cur_table = { null };
        final int[] key_col = { 0 };
        final NodeDataImporter.MatchBy[] match_by = { NodeDataImporter.MatchBy.TIP_NAME };
        final List<JCheckBox> includes = new ArrayList<>();
        final List<JTextField> renames = new ArrayList<>();

        final JComboBox<String> delim_combo = new JComboBox<>(new String[] { "Auto", "Comma", "Tab" });
        final JPanel content = new JPanel(new java.awt.BorderLayout());

        final Runnable rebuild = () -> {
            content.removeAll();
            includes.clear();
            renames.clear();
            final Character forced = (delim_combo.getSelectedIndex() == 1) ? Character.valueOf(',')
                    : ((delim_combo.getSelectedIndex() == 2) ? Character.valueOf('\t') : null);
            final NodeDataImporter.Table table;
            try {
                table = NodeDataImporter.parseTable(text, forced);
            }
            catch (final IllegalArgumentException ex) {
                cur_table[0] = null;
                content.add(new JLabel("Cannot parse with this delimiter: " + ex.getMessage()),
                        java.awt.BorderLayout.NORTH);
                content.revalidate();
                content.repaint();
                return;
            }
            cur_table[0] = table;
            key_col[0] = table.defaultKeyColumn();
            final String[] headers = table.getHeaders();
            final int rows_shown = Math.min(8, table.getRowCount());
            final String[][] preview_data = new String[rows_shown][headers.length];
            for (int r = 0; r < rows_shown; r++) {
                for (int c = 0; c < headers.length; c++) {
                    preview_data[r][c] = table.getCell(r, c);
                }
            }
            final JTable preview = new JTable(preview_data, headers) {
                @Override
                public boolean isCellEditable(final int r, final int c) {
                    return false;
                }
            };
            preview.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
            final JScrollPane preview_scroll = new JScrollPane(preview);
            preview_scroll.setPreferredSize(new Dimension(580, 130));
            final JComboBox<String> key_combo = new JComboBox<>(headers);
            key_combo.setSelectedIndex(key_col[0]);
            final JComboBox<NodeDataImporter.MatchBy> match_combo = new JComboBox<>(NodeDataImporter.userMatchOptions());
            match_combo.setSelectedItem(match_by[0]); // preserve the user's match attribute across a delimiter re-parse
            final JLabel summary = new JLabel();
            // one row per column: [include?] header -> [rename field]; the key column's row is disabled
            final JPanel map = new JPanel(new GridLayout(0, 1, 0, 2));
            for (int c = 0; c < headers.length; c++) {
                final JCheckBox cb = new JCheckBox("", true);
                final JTextField tf = new JTextField(headers[c], 14);
                includes.add(cb);
                renames.add(tf);
                final JPanel row = new JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 4, 0));
                row.add(cb);
                row.add(new JLabel(headers[c] + "  →"));
                row.add(tf);
                map.add(row);
            }
            final JScrollPane map_scroll = new JScrollPane(map);
            map_scroll.setPreferredSize(new Dimension(580, 120));
            final Runnable refresh = () -> summary.setText(NodeDataImporter
                    .dryRun(phy, table, key_col[0], match_by[0], planFromRows(table, key_col[0], includes, renames))
                    .summaryLine());
            final Runnable sync_key_row = () -> {
                for (int c = 0; c < headers.length; c++) {
                    includes.get(c).setEnabled(c != key_col[0]);
                    renames.get(c).setEnabled(c != key_col[0]);
                }
            };
            key_combo.addActionListener(e -> {
                key_col[0] = key_combo.getSelectedIndex();
                sync_key_row.run();
                refresh.run();
            });
            match_combo.addActionListener(e -> {
                match_by[0] = (NodeDataImporter.MatchBy) match_combo.getSelectedItem();
                refresh.run();
            });
            for (final JCheckBox cb : includes) {
                cb.addActionListener(e -> refresh.run());
            }
            sync_key_row.run();
            refresh.run();
            final JPanel controls = new JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 4, 0));
            controls.add(new JLabel("Match table column:"));
            controls.add(key_combo);
            controls.add(new JLabel("against tip:"));
            controls.add(match_combo);
            final JPanel body = new JPanel();
            body.setLayout(new BoxLayout(body, BoxLayout.Y_AXIS));
            addLeft(body, new JLabel("Preview (first " + rows_shown + " of " + table.getRowCount() + " rows):"));
            addLeft(body, preview_scroll);
            addLeft(body, controls);
            addLeft(body, new JLabel("Columns to import (uncheck to skip; edit a name to rename its property):"));
            addLeft(body, map_scroll);
            addLeft(body, summary);
            content.add(body, java.awt.BorderLayout.CENTER);
            content.revalidate();
            content.repaint();
            dialog.pack();
        };
        delim_combo.addActionListener(e -> rebuild.run());
        rebuild.run();

        final JButton ok = new JButton("Import");
        final JButton cancel = new JButton("Cancel");
        ok.addActionListener(e -> {
            final NodeDataImporter.Table table = cur_table[0];
            if (table != null) {
                result[0] = new ImportChoice(table, key_col[0], match_by[0],
                        planFromRows(table, key_col[0], includes, renames));
            }
            dialog.dispose();
        });
        cancel.addActionListener(e -> dialog.dispose());

        final JPanel top = new JPanel(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 4, 0));
        top.add(new JLabel("Delimiter:"));
        top.add(delim_combo);
        final JPanel buttons = new JPanel();
        buttons.add(ok);
        buttons.add(cancel);
        dialog.getContentPane().setLayout(new java.awt.BorderLayout(8, 8));
        dialog.getContentPane().add(top, java.awt.BorderLayout.NORTH);
        dialog.getContentPane().add(content, java.awt.BorderLayout.CENTER);
        dialog.getContentPane().add(buttons, java.awt.BorderLayout.SOUTH);
        dialog.pack();
        dialog.setLocationRelativeTo(this);
        dialog.setVisible(true); // modal: blocks until OK/Cancel
        return result[0];
    }

    /** Build a ColumnPlan from the dialog's per-column checkbox/rename rows (extracts the UI state, then delegates
     *  to the pure {@link #buildColumnPlan}). One source of truth so the live dry-run preview and the committed
     *  import can never disagree. */
    private static NodeDataImporter.ColumnPlan planFromRows(final NodeDataImporter.Table table, final int key_col,
            final List<JCheckBox> includes, final List<JTextField> renames) {
        final boolean[] included = new boolean[table.getColumnCount()];
        final String[] headers = new String[table.getColumnCount()];
        for (int c = 0; c < table.getColumnCount(); c++) {
            included[c] = includes.get(c).isSelected();
            headers[c] = renames.get(c).getText();
        }
        return buildColumnPlan(table, key_col, included, headers);
    }

    /** Pure (testable) ColumnPlan assembly: the key column and any unchecked column are excluded from the import;
     *  each column's effective header is its (possibly renamed) text. */
    static NodeDataImporter.ColumnPlan buildColumnPlan(final NodeDataImporter.Table table, final int key_col,
            final boolean[] included, final String[] headers) {
        final NodeDataImporter.ColumnPlan plan = NodeDataImporter.ColumnPlan.importAll(table);
        for (int c = 0; c < table.getColumnCount(); c++) {
            plan.setIncluded(c, (c != key_col) && included[c]);
            plan.setHeader(c, headers[c]);
        }
        return plan;
    }

    private static void addLeft(final JPanel box, final Component c) {
        ((javax.swing.JComponent) c).setAlignmentX(Component.LEFT_ALIGNMENT);
        box.add(c);
    }

    /**
     * UI-free apply + refresh (so it is unit-testable without the modal chooser): applies the table under the given
     * column plan, and -- only when it actually annotated something -- checkpoints undo, appends a provenance
     * sentence, re-lays-out the tree, and refreshes the Color-by / Size-by dropdowns so the columns are usable.
     */
    NodeDataImporter.ImportResult importAnnotationsAndRefit(final Phylogeny phy, final NodeDataImporter.Table table,
            final int key_col, final NodeDataImporter.MatchBy match_by, final NodeDataImporter.ColumnPlan plan,
            final String source_name) {
        final TreePanel tp = getCurrentTreePanel();
        final Phylogeny before = (tp != null) ? phy.copy() : null;
        final boolean was_edited = (tp != null) && tp.isEdited();
        final int total_tips = phy.getNumberOfExternalNodes();
        final NodeDataImporter.ImportResult result = NodeDataImporter.apply(phy, table, key_col, match_by, plan);
        if ((result.getTipsAnnotated() > 0) && (tp != null)) {
            tp.pushUndoSnapshot(before, was_edited, "Import Annotations"); // now we know it changed the tree
            final String prov = importProvenance(result.getPropertyColumns(), result.getTipsAnnotated(), total_tips,
                    source_name, match_by);
            final String existing = phy.getDescription();
            phy.setDescription(ForesterUtil.isEmpty(existing) ? prov : existing + " " + prov);
            tp.setTree(phy); // recompute the layout so the new labels/data show
            tp.getControlPanel().populateColorByPropertyBox(); // surface any imported columns in "Color by:"
            tp.getControlPanel().populateSizeByPropertyBox(); // and any numeric ones in "Size by:"
            tp.getControlPanel().populateAncestralPieBox(); // and any discrete-trait ones in "Ancestral pie:"
            tp.getControlPanel().rebuildSearchFields(true); // and as searchable fields (forced: same tree, new data)
            showWhole();
            tp.setEdited(true);
        }
        return result;
    }

    /** Pure provenance sentence appended to the tree description on an annotation import (per the repo rule). */
    static String importProvenance(final List<String> property_columns, final int tips_annotated,
            final int total_tips, final String source_name, final NodeDataImporter.MatchBy match_by) {
        final StringBuilder sb = new StringBuilder();
        sb.append("Imported annotations from table \"").append(source_name).append("\" onto ").append(tips_annotated)
          .append(" of ").append(total_tips).append(total_tips == 1 ? " tip" : " tips").append(" (matched by ")
          .append(match_by.toString().toLowerCase()).append(").");
        if (!property_columns.isEmpty()) {
            sb.append(" Columns: ").append(String.join(", ", property_columns)).append(".");
        }
        return sb.toString();
    }

    /**
     * File -> Import GTDB Taxonomy: read a GTDB-Tk-style table (a tip-key column + a GTDB classification column
     * {@code d__…;…;s__…}) and write the genome-based taxonomy onto the matching tips -- each rank as a
     * {@code gtdb:<rank>} property + a taxonomy at the most specific rank -- so GTDB (the bacterial/archaeal standard)
     * drives Color-by / Annotation Columns / search entirely offline. Undoable.
     */
    void importGtdbTaxonomy() {
        final Phylogeny phy = currentPhylogenyForExport();
        if (phy == null) {
            return;
        }
        final JFileChooser fc = new JFileChooser();
        fc.setMultiSelectionEnabled(false);
        fc.setDialogTitle("Import GTDB Taxonomy (GTDB-Tk table)");
        fc.setFileFilter(new FileNameExtensionFilter("GTDB-Tk / GTDB tables (*.tsv, *.csv, *.txt)", "tsv", "csv", "txt"));
        if (getCurrentDir(DirectoryPreferences.Category.OPEN) != null) {
            fc.setCurrentDirectory(getCurrentDir(DirectoryPreferences.Category.OPEN));
        }
        if (fc.showOpenDialog(this) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        final File file = fc.getSelectedFile();
        if (file == null) {
            return;
        }
        setCurrentDir(DirectoryPreferences.Category.OPEN, fc.getCurrentDirectory());
        final String text;
        try {
            text = Files.readString(file.toPath());
        }
        catch (final IOException e) {
            JOptionPane.showMessageDialog(this, "Failed to read " + file + ":\n" + e.getMessage(), "Read Failed",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        final int annotated = importGtdbAndRefit(phy, text, file.getName());
        if (annotated < 0) {
            JOptionPane.showMessageDialog(this,
                    "Could not read a GTDB classification from \"" + file.getName() + "\".\n"
                            + "Expected a table (TSV/CSV) with a tip-name column and a GTDB classification column "
                            + "(values like d__Bacteria;p__…;s__…), e.g. a GTDB-Tk summary.",
                    "Import GTDB Taxonomy", JOptionPane.WARNING_MESSAGE);
        }
        else if (annotated == 0) {
            JOptionPane.showMessageDialog(this,
                    "No tips matched the table's key column by name — nothing was imported.\n"
                            + "The classification column was found, but no tip name equalled a row's key.",
                    "Import GTDB Taxonomy", JOptionPane.WARNING_MESSAGE);
        }
        else {
            JOptionPane.showMessageDialog(this,
                    "Imported GTDB taxonomy onto " + annotated + " of " + phy.getNumberOfExternalNodes()
                            + " tips.\nColor by (or add an Annotation Column for) gtdb:phylum / gtdb:family / gtdb:genus …",
                    "Import GTDB Taxonomy", JOptionPane.INFORMATION_MESSAGE);
        }
    }

    /**
     * Dialog-free testable core of {@link #importGtdbTaxonomy()}: parse the table, find the GTDB classification column
     * (values look like {@code d__…}) and the key column (a genome-id / node_id / name column, else the first
     * non-classification column -- see {@link #gtdbKeyColumn}), apply the classifications to tips matched BY NAME, and
     * -- only when it changed the tree -- checkpoint undo, append a provenance sentence, re-lay-out, and refresh the
     * Color-by / search inputs. Returns the number of tips annotated, or -1 on a parse error / no GTDB classification
     * column.
     */
    int importGtdbAndRefit(final Phylogeny phy, final String text, final String source_name) {
        final NodeDataImporter.Table table;
        try {
            table = NodeDataImporter.parseTable(text);
        }
        catch (final Exception e) {
            return -1;
        }
        final int class_col = gtdbClassificationColumn(table);
        if (class_col < 0) {
            return -1;
        }
        final int key_col = gtdbKeyColumn(table, class_col);
        final java.util.Map<String, String> map = new java.util.HashMap<String, String>();
        for (int r = 0; r < table.getRowCount(); ++r) {
            final String k = table.getCell(r, key_col);
            final String c = table.getCell(r, class_col);
            if ((k != null) && (k.trim().length() > 0) && (c != null)) {
                map.put(k.trim(), c);
            }
        }
        final TreePanel tp = getCurrentTreePanel();
        final Phylogeny before = (tp != null) ? phy.copy() : null;
        final boolean was_edited = (tp != null) && tp.isEdited();
        final int total_tips = phy.getNumberOfExternalNodes();
        final int annotated;
        try {
            annotated = org.forester.archaeopteryx.tools.GtdbTaxonomy.applyByTipName(phy, map);
        }
        catch (final Exception e) {
            return -1;
        }
        if ((annotated > 0) && (tp != null)) {
            tp.pushUndoSnapshot(before, was_edited, "Import GTDB Taxonomy");
            final String prov = "Imported GTDB taxonomy from table \"" + source_name + "\" onto " + annotated + " of "
                    + total_tips + (total_tips == 1 ? " tip" : " tips") + " (gtdb:domain … gtdb:species).";
            final String existing = phy.getDescription();
            phy.setDescription(ForesterUtil.isEmpty(existing) ? prov : existing + " " + prov);
            tp.setTree(phy);
            tp.getControlPanel().populateColorByPropertyBox(); // surface the gtdb:* ranks in "Color by:"
            tp.getControlPanel().populateSizeByPropertyBox(); // (consistency with the annotation-import refresh path)
            tp.getControlPanel().populateAncestralPieBox();
            tp.getControlPanel().rebuildSearchFields(true); // and as searchable fields (forced: same tree, new data)
            showWhole();
            tp.setEdited(true);
        }
        return annotated;
    }

    /** The first table column whose values look like a GTDB classification string ({@code d__…;p__…;s__…}); -1 if
     *  none. Package-visible for {@link ImportGtdbToolTest}. */
    static int gtdbClassificationColumn(final NodeDataImporter.Table table) {
        for (int c = 0; c < table.getColumnCount(); ++c) {
            for (int r = 0; r < table.getRowCount(); ++r) {
                if (org.forester.archaeopteryx.tools.GtdbTaxonomy.looksLikeGtdb(table.getCell(r, c))) {
                    return c;
                }
            }
        }
        return -1;
    }

    /**
     * The column to key the join on (matched against tip names), given the detected {@code class_col}. Prefers a GTDB-Tk
     * genome-id header (user_genome / genome / accession / bin_id / bin), else the table's node_id/name key column, else
     * the first column that is not the classification -- so an unrelated "name" column can't hijack the key on a
     * GTDB-Tk table that also carries one. Package-visible for {@link ImportGtdbToolTest}.
     */
    static int gtdbKeyColumn(final NodeDataImporter.Table table, final int class_col) {
        final String[] headers = table.getHeaders();
        for (int c = 0; c < headers.length; ++c) {
            if (c == class_col) {
                continue;
            }
            final String h = (headers[c] == null) ? "" : headers[c].trim().toLowerCase(java.util.Locale.ROOT);
            if (h.equals("user_genome") || h.equals("genome") || h.equals("accession") || h.equals("bin_id")
                    || h.equals("bin")) {
                return c;
            }
        }
        final int def = table.defaultKeyColumn(); // a node_id/name column if the table has one (else 0)
        if ((def >= 0) && (def != class_col)) {
            return def;
        }
        for (int c = 0; c < table.getColumnCount(); ++c) {
            if (c != class_col) {
                return c;
            }
        }
        return class_col; // single-column table (degenerate; applyByTipName then matches the classification as a name)
    }

    /**
     * A chosen export scope: the tips to export, a human phrase for the completion message, and whether the
     * scope is a restricted subset of the whole loaded tree (sub-tree display or a selection) -- in which case
     * the suggested file name gets a {@code _N} element-count suffix.
     */
    private static final class ExportScope {
        final List<PhylogenyNode> tips;
        final String              description;
        final boolean             restricted;

        ExportScope(final List<PhylogenyNode> tips, final String description, final boolean restricted) {
            this.tips = tips;
            this.description = description;
            this.restricted = restricted;
        }
    }

    /**
     * Decide which tips a data export covers. The scope is the selected EXTERNAL tips (branch-click a clade to
     * select its tips; a selected internal node contributes nothing). With no tips selected the whole displayed
     * (sub)tree is used silently; otherwise the user picks Selected / Not-selected / All displayed / Cancel.
     * Returns null only if the user cancelled.
     */
    private ExportScope chooseExportTips(final Phylogeny phy) {
        final List<PhylogenyNode> all = phy.getExternalNodes();
        final TreePanel tp = _mainpanel.getCurrentTreePanel();
        final boolean from_subtree = (tp != null) && tp.isCurrentTreeIsSubtree();
        final String all_desc = all.size() + " external nodes"
                + (from_subtree ? " in the currently displayed subtree" : "");
        final List<PhylogenyNode> selected = (tp == null) ? new ArrayList<>()
                : NodeDataExporter.selectedExternalTips(phy, tp.getFoundNodesAsListOfPhylogenyNodes());
        if (selected.isEmpty()) {
            // no selection -> whole displayed tree, no prompt; restricted only if that tree is a sub-tree
            return new ExportScope(all, all_desc, from_subtree);
        }
        final List<PhylogenyNode> not_selected = NodeDataExporter.complementExternalTips(all, selected);
        final Object[] options = {"Selected tips (" + selected.size() + ")",
                "Not-selected tips (" + not_selected.size() + ")", "All displayed tips (" + all.size() + ")",
                "Cancel"};
        // FlatLaf/macOS lays option-dialog buttons right-to-left (default on the right). For this multi-choice
        // picker a left-to-right reading order (Selected, Not-selected, All, Cancel) is clearer, so force
        // forward button order for just this dialog and restore the previous setting afterwards.
        final Object prev_yes_last = UIManager.get("OptionPane.isYesLast");
        UIManager.put("OptionPane.isYesLast", Boolean.FALSE);
        final int choice;
        try {
            choice = JOptionPane.showOptionDialog(this,
                    selected.size() + " of " + all.size() + " tips are selected. Which tips do you want to export?",
                    "Export Scope", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, options,
                    options[0]);
        }
        finally {
            UIManager.put("OptionPane.isYesLast", prev_yes_last);
        }
        switch (choice) {
            case 0:
                return new ExportScope(selected, selected.size() + " selected tips", true);
            case 1:
                return new ExportScope(not_selected, not_selected.size() + " not-selected tips", true);
            case 2:
                return new ExportScope(all, all_desc, from_subtree);
            default:
                return null; // Cancel or dialog closed
        }
    }

    private Phylogeny currentPhylogenyForExport() {
        if (_mainpanel.getCurrentTreePanel() == null) {
            return null;
        }
        final Phylogeny phy = _mainpanel.getCurrentPhylogeny();
        return ((phy == null) || phy.isEmpty()) ? null : phy;
    }

    /** A suggested file name for a data export: the loaded tree file's base name (or the tree name) + ext. */
    private String suggestedExportName(final Phylogeny phy, final String ext, final int count) {
        String base = null;
        if ((_mainpanel.getCurrentTreePanel() != null) && (_mainpanel.getCurrentTreePanel().getTreeFile() != null)) {
            base = _mainpanel.getCurrentTreePanel().getTreeFile().getName();
            final int dot = base.lastIndexOf('.');
            if (dot > 0) {
                base = base.substring(0, dot);
            }
        }
        if (ForesterUtil.isEmpty(base)) {
            base = ForesterUtil.isEmpty(phy.getName()) ? "tree_data" : phy.getName();
        }
        // A restricted export (count >= 0) tags the name with how many elements it holds, e.g. flavi_43.tsv.
        final String suffix = (count >= 0) ? ("_" + count) : "";
        return base.replaceAll("[^A-Za-z0-9._-]+", "_") + suffix + ext;
    }

    /** Shared save flow for the read-only data exports: pick a file (overwrite-confirmed) and write the text. */
    private void writeDataExportToFile(final String content, final String what, final int count,
                                       final String scope_desc, final String suggested_name, final String note) {
        if (ForesterUtil.isEmpty(content)) {
            JOptionPane.showMessageDialog(this, "There is no " + what + " in this tree to export.",
                    "Nothing to Export", JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        final JFileChooser fc = new JFileChooser();
        fc.setMultiSelectionEnabled(false);
        if (getCurrentDir(DirectoryPreferences.Category.EXPORT) != null) {
            fc.setCurrentDirectory(getCurrentDir(DirectoryPreferences.Category.EXPORT));
        }
        if (!ForesterUtil.isEmpty(suggested_name)) {
            fc.setSelectedFile(new File(suggested_name));
        }
        if (fc.showSaveDialog(this) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        final File file = fc.getSelectedFile();
        if (file == null) {
            return;
        }
        if (file.exists() && (JOptionPane.showConfirmDialog(this, file + " already exists.\nOverwrite?",
                "Overwrite?", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE) != JOptionPane.OK_OPTION)) {
            return;
        }
        try (final Writer w = new BufferedWriter(new FileWriter(file))) {
            w.write(content);
        }
        catch (final IOException e) {
            JOptionPane.showMessageDialog(this, "Failed to write " + file + ":\n" + e.getMessage(),
                    "Write Failed", JOptionPane.ERROR_MESSAGE);
            return;
        }
        setCurrentDir(DirectoryPreferences.Category.EXPORT, fc.getCurrentDirectory());
        final String msg = "Wrote " + count + " " + what + " from " + scope_desc + " to:\n" + file
                + (ForesterUtil.isEmpty(note) ? "" : "\n\n" + note);
        JOptionPane.showMessageDialog(this, msg, "Export Complete", JOptionPane.INFORMATION_MESSAGE);
    }

    static File writeToFile(final Phylogeny t,
                            final MainPanel mp,
                            final JFileChooser save_filechooser,
                            final File current_dir,
                            final Container contentpane,
                            final Component comp) {
        File new_file = null;
        if (t == null) {
            return null;
        }
        String initial_filename = null;
        if (mp.getCurrentTreePanel().getTreeFile() != null) {
            try {
                initial_filename = mp.getCurrentTreePanel().getTreeFile().getCanonicalPath();
            } catch (final IOException e) {
                initial_filename = null;
            }
        }
        if (!ForesterUtil.isEmpty(initial_filename)) {
            save_filechooser.setSelectedFile(new File(initial_filename));
        } else {
            save_filechooser.setSelectedFile(new File(""));
        }
        final File my_dir = current_dir;
        if (my_dir != null) {
            save_filechooser.setCurrentDirectory(my_dir);
        }
        final int result = save_filechooser.showSaveDialog(contentpane);
        final File file = save_filechooser.getSelectedFile();
        new_file = save_filechooser.getCurrentDirectory();
        boolean exception = false;
        if ((file != null) && (result == JFileChooser.APPROVE_OPTION)) {
            if (file.exists()) {
                final int i = JOptionPane.showConfirmDialog(comp,
                        file + " already exists.\nOverwrite?",
                        "Overwrite?",
                        JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.QUESTION_MESSAGE);
                if (i != JOptionPane.OK_OPTION) {
                    return null;
                } else {
                    final File to = new File(file.getAbsoluteFile().toString() + AptxConstants.BACKUP_FILE_SUFFIX);
                    try {
                        ForesterUtil.copyFile(file, to);
                    } catch (final Exception e) {
                        JOptionPane.showMessageDialog(comp,
                                "Failed to create backup copy " + to,
                                "Failed to Create Backup Copy",
                                JOptionPane.WARNING_MESSAGE);
                    }
                    try {
                        file.delete();
                    } catch (final Exception e) {
                        JOptionPane.showMessageDialog(comp,
                                "Failed to delete: " + file,
                                "Failed to Delete",
                                JOptionPane.WARNING_MESSAGE);
                    }
                }
            }
            if (save_filechooser.getFileFilter() == MainFrame.nhfilter) {
                exception = writeAsNewHampshire(mp.getCurrentTreePanel(), mp.getOptions(), exception, file);
            } else if (save_filechooser.getFileFilter() == MainFrame.xmlfilter) {
                exception = writeAsPhyloXml(mp.getCurrentTreePanel(), mp.getOptions(), exception, file);
            } else if (save_filechooser.getFileFilter() == MainFrame.nexusfilter) {
                exception = writeAsNexus(mp.getCurrentTreePanel(), mp.getOptions(), exception, file);
            }
            // "*.*":
            else {
                final String file_name = file.getName().trim().toLowerCase();
                if (file_name.endsWith(".nh") || file_name.endsWith(".newick") || file_name.endsWith(".phy")
                        || file_name.endsWith(".tree")) {
                    exception = writeAsNewHampshire(mp.getCurrentTreePanel(), mp.getOptions(), exception, file);
                } else if (file_name.endsWith(".nex") || file_name.endsWith(".nexus")) {
                    exception = writeAsNexus(mp.getCurrentTreePanel(), mp.getOptions(), exception, file);
                }
                // XML is default:
                else {
                    exception = writeAsPhyloXml(mp.getCurrentTreePanel(), mp.getOptions(), exception, file);
                }
            }
            if (!exception) {
                mp.setTitleOfSelectedTab(file.getName());
                mp.getCurrentTreePanel().setTreeFile(file);
                mp.getCurrentTreePanel().setEdited(false);
            }
        }
        return new_file;
    }

    static File writeToGraphicsFile(final Phylogeny t,
                                    final GraphicsExportType type,
                                    final MainPanel mp,
                                    final JFileChooser writetographics_filechooser,
                                    final Component component,
                                    final Container contentpane,
                                    final File current_dir) {
        File new_dir = null;
        if ((t == null) || t.isEmpty()) {
            return null;
        }
        String initial_filename = "";
        if (mp.getCurrentTreePanel().getTreeFile() != null) {
            initial_filename = mp.getCurrentTreePanel().getTreeFile().toString();
        }
        if (initial_filename.indexOf('.') > 0) {
            initial_filename = initial_filename.substring(0, initial_filename.lastIndexOf('.'));
        }
        initial_filename = initial_filename + "." + type;
        writetographics_filechooser.setSelectedFile(new File(initial_filename));
        final File my_dir = current_dir;
        if (my_dir != null) {
            writetographics_filechooser.setCurrentDirectory(my_dir);
        }
        final int result = writetographics_filechooser.showSaveDialog(contentpane);
        File file = writetographics_filechooser.getSelectedFile();
        //setCurrentDir( writetographics_filechooser.getCurrentDirectory() );
        new_dir = writetographics_filechooser.getCurrentDirectory();
        if ((file != null) && (result == JFileChooser.APPROVE_OPTION)) {
            if (!file.toString().toLowerCase().endsWith(type.toString())) {
                file = new File(file.toString() + "." + type);
            }
            if (file.exists()) {
                final int i = JOptionPane.showConfirmDialog(component,
                        file + " already exists. Overwrite?",
                        "Warning",
                        JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.WARNING_MESSAGE);
                if (i != JOptionPane.OK_OPTION) {
                    return null;
                } else {
                    try {
                        file.delete();
                    } catch (final Exception e) {
                        JOptionPane.showMessageDialog(component,
                                "Failed to delete: " + file,
                                "Error",
                                JOptionPane.WARNING_MESSAGE);
                    }
                }
            }
            writePhylogenyToGraphicsFile(file.toString(), type, mp, component, contentpane);
        }
        return new_dir;
    }

    static File writeToPdf(final Phylogeny t,
                           final MainPanel mp,
                           final JFileChooser writetopdf_filechooser,
                           final File curr_dir,
                           final Container contentpane,
                           final Component component) {
        if ((t == null) || t.isEmpty()) {
            return null;
        }
        String initial_filename = "";
        if (mp.getCurrentTreePanel().getTreeFile() != null) {
            initial_filename = mp.getCurrentTreePanel().getTreeFile().toString();
        }
        if (initial_filename.indexOf('.') > 0) {
            initial_filename = initial_filename.substring(0, initial_filename.lastIndexOf('.'));
        }
        initial_filename = initial_filename + ".pdf";
        writetopdf_filechooser.setSelectedFile(new File(initial_filename));
        final File my_dir = curr_dir;
        if (my_dir != null) {
            writetopdf_filechooser.setCurrentDirectory(my_dir);
        }
        final int result = writetopdf_filechooser.showSaveDialog(contentpane);
        File file = writetopdf_filechooser.getSelectedFile();
        // setCurrentDir( writetopdf_filechooser.getCurrentDirectory() );
        final File new_current_dir = writetopdf_filechooser.getCurrentDirectory();
        if ((file != null) && (result == JFileChooser.APPROVE_OPTION)) {
            if (!file.toString().toLowerCase().endsWith(".pdf")) {
                file = new File(file.toString() + ".pdf");
            }
            if (file.exists()) {
                final int i = JOptionPane.showConfirmDialog(component,
                        file + " already exists. Overwrite?",
                        "WARNING",
                        JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.WARNING_MESSAGE);
                if (i != JOptionPane.OK_OPTION) {
                    return null;
                }
            }
            exportPhylogenyToPdf(file.toString(), mp.getOptions(), mp.getCurrentTreePanel(), component);
        }
        return new_current_dir;
    }
}

class DefaultFilter extends FileFilter {

    @Override
    public boolean accept(final File f) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith(".nh") || file_name.endsWith(".newick") || file_name.endsWith(".phy")
                || file_name.endsWith(".nwk") || file_name.endsWith(".phb") || file_name.endsWith(".ph")
                || file_name.endsWith(".tr") || file_name.endsWith(".dnd") || file_name.endsWith(".tree")
                || file_name.endsWith(".nhx") || file_name.endsWith(".xml") || file_name.endsWith(".phyloxml")
                || file_name.endsWith("phylo.xml") || file_name.endsWith(".pxml") || file_name.endsWith(".nexus")
                || file_name.endsWith(".nx") || file_name.endsWith(".nex") || file_name.endsWith(".tre")
                || file_name.endsWith(".zip") || file_name.endsWith(".tol") || file_name.endsWith(".tolxml")
                || file_name.endsWith(".con") || file_name.endsWith(".json") || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "All supported files (*.xml, *.phyloxml, *phylo.xml, *.nhx, *.nh, *.newick, *.nex, *.nexus, *.phy, *.tre, *.tree, *.tol, *.json, ...)";
    }
}

class GraphicsFileFilter extends FileFilter {

    @Override
    public boolean accept(final File f) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith(".jpg") || file_name.endsWith(".jpeg") || file_name.endsWith(".png")
                || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "Image files (*.jpg, *.jpeg, *.png)";
    }
}

class NexusFilter extends FileFilter {

    @Override
    public boolean accept(final File f) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith(".nex") || file_name.endsWith(".nexus") || file_name.endsWith(".nx")
                || file_name.endsWith(".tre") || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "Nexus files (*.nex, *.nexus, *.nx, *.tre)";
    }
} // NexusFilter

class NHFilter extends FileFilter {

    @Override
    public boolean accept(final File f) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith(".nh") || file_name.endsWith(".newick") || file_name.endsWith(".phy")
                || file_name.endsWith(".tr") || file_name.endsWith(".tree") || file_name.endsWith(".dnd")
                || file_name.endsWith(".ph") || file_name.endsWith(".phb") || file_name.endsWith(".nwk")
                || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "New Hampshire - Newick files (*.nh, *.newick, *.phy, *.tree, *.dnd, *.tr, *.ph, *.phb, *.nwk)";
    }
} // NHFilter

class NHXFilter extends FileFilter {

    @Override
    public boolean accept(final File f) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith(".nhx") || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "NHX files (*.nhx) [deprecated]";
    }
}

class PdfFilter extends FileFilter {

    @Override
    public boolean accept(final File f) {
        return f.getName().trim().toLowerCase().endsWith(".pdf") || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "PDF files (*.pdf)";
    }
} // PdfFilter

class TolFilter extends FileFilter {

    @Override
    public boolean accept(final File f) {
        final String file_name = f.getName().trim().toLowerCase();
        return (file_name.endsWith(".tol") || file_name.endsWith(".tolxml") || file_name.endsWith(".zip") || f
                .isDirectory()) && (!file_name.endsWith(".xml.zip"));
    }

    @Override
    public String getDescription() {
        return "Tree of Life files (*.tol, *.tolxml)";
    }
} // TolFilter

class XMLFilter extends FileFilter {

    @Override
    public boolean accept(final File f) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith(".xml") || file_name.endsWith(".phyloxml") || file_name.endsWith("phylo.xml")
                || file_name.endsWith(".pxml") || file_name.endsWith(".zip") || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "phyloXML files (*.xml, *.phyloxml, *phylo.xml, *.pxml, *.zip)";
    }
} // XMLFilter

class JsonFilter extends FileFilter {

    @Override
    public boolean accept(final File f) {
        // .auspice.json ends with .json, so the single check covers both
        return f.getName().trim().toLowerCase().endsWith(".json") || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "Auspice / Nextstrain JSON files (*.json, *.auspice.json)";
    }
} // JsonFilter
