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
import java.awt.Container;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.SortedSet;
import java.util.NoSuchElementException;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
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
import org.forester.archaeopteryx.tools.AncestralTaxonomyInferrer;
import org.forester.archaeopteryx.tools.LabelDataExtractor;
import org.forester.archaeopteryx.tools.ProcessPool;
import org.forester.archaeopteryx.tools.ProcessRunning;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.PhylogenyNode;
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
     * Installs the given look-and-feel. FlatLaf (light/dark) is the modern default;
     * the native and cross-platform look-and-feels are kept as alternatives.
     */
    static void installLookAndFeel(final Configuration.UI ui) {
        try {
            switch (ui) {
                case NATIVE:
                    UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
                    break;
                case CROSSPLATFORM:
                    UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
                    break;
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
    final static TolFilter tolfilter = new TolFilter();
    final static NexusFilter nexusfilter = new NexusFilter();
    final static PdfFilter pdffilter = new PdfFilter();
    final static GraphicsFileFilter graphicsfilefilter = new GraphicsFileFilter();
    final static DefaultFilter defaultfilter = new DefaultFilter();
    static final String USE_MOUSEWHEEL_SHIFT_TO_ROTATE = "rotate with mousewheel + Shift (or A and S), D toggles between horizontal and radial labels";
    static final String PHYLOXML_REF_TOOL_TIP = AptxConstants.PHYLOXML_REFERENCE;                                                                                                                                                //TODO //FIXME
    static final String APTX_REF_TOOL_TIP = AptxConstants.APTX_REFERENCE;
    private static final long serialVersionUID = 3655000897845508358L;
    final static Font menu_font = new Font(Configuration.getDefaultFontFamilyName(),
            Font.PLAIN,
            Configuration.getGuiFontSize());
    static final String TYPE_MENU_HEADER = "Type";
    static final String RECTANGULAR_TYPE_CBMI_LABEL = "Rectangular";
    static final String EURO_TYPE_CBMI_LABEL = "Euro Type";
    static final String CURVED_TYPE_CBMI_LABEL = "Curved";
    static final String TRIANGULAR_TYPE_CBMI_LABEL = "Triangular";
    static final String CONVEX_TYPE_CBMI_LABEL = "Convex";
    static final String ROUNDED_TYPE_CBMI_LABEL = "Rounded";
    static final String UNROOTED_TYPE_CBMI_LABEL = "Unrooted (alpha)";                                                                                                                                                          //TODO
    static final String CIRCULAR_TYPE_CBMI_LABEL = "Circular (alpha)";                                                                                                                                                          //TODO
    static final String SEARCH_TERMS_ONLY_LABEL = "Words";
    static final String SEARCH_REGEX_LABEL = "Regex";
    static final String SEARCH_CASE_SENSITIVE_LABEL = "Match Case";
    static final String INVERSE_SEARCH_RESULT_LABEL = "Inverse";
    static final String DISPLAY_SCALE_LABEL = "Scale";
    static final String NON_LINED_UP_CLADOGRAMS_LABEL = "Non-Lined Up Cladogram";
    static final String LABEL_DIRECTION_LABEL = "Radial Labels";
    static final String LABEL_DIRECTION_TIP = "To use radial node labels in radial and unrooted display types";
    static final String SEARCH_WITH_REGEX_TIP = "To search using regular expressions (~Java/Perl syntax). For example, use \"^B.+\\d{2,}$\" to search for everything starting with a B and ending with at least two digits.";
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
    static final String SHOW_DOMAIN_LABELS_LABEL = "Domain Labels";
    static final String COLOR_LABELS_TIP = "To use parent branch colors for node labels as well, need to turn off taxonomy dependent colorization and turn on branch colorization for this to become apparent";
    static final String ABBREV_SN_LABEL = "Abbreviate Scientific Taxonomic Names";
    static final String ITALIC_SN_LABEL = "Italic Scientific Taxonomic Names";
    static final String ITALIC_SN_TIP = "To set scientific (Latin) taxonomic names in italics (e.g. Homo sapiens), per the publication convention; only the scientific-name part is italicized, not the code, common name or rank";
    static final String OUTLINE_FONTS_VECTOR_LABEL = "Outline Fonts (SVG/EPS export)";
    static final String OUTLINE_FONTS_VECTOR_TIP = "SVG/EPS export embeds no fonts; outlining renders text as vector shapes so viewers cannot substitute the figure font (recommended). Turn off to keep selectable/searchable text, at the risk of font substitution.";
    static final String TRANSPARENT_BG_LABEL = "Transparent Background (PNG)";
    static final String TRANSPARENT_BG_TIP = "Export PNG with a transparent (alpha) background instead of the solid figure background. Ignored for formats that cannot carry transparency (JPG, etc.).";
    static final String SHOW_CONF_STDDEV_LABEL = "Confidence Standard Deviations";
    static final String SHOW_MAD_CONF_LABEL    = "MAD Confidence Values (MAD/regular)";
    static final String USE_BRACKETS_FOR_CONF_IN_NH_LABEL = "Use Brackets for Confidence Values";
    static final String USE_INTERNAL_NAMES_FOR_CONF_IN_NH_LABEL = "Use Internal Node Names for Confidence Values";
    static final String SHOW_BASIC_TREE_INFORMATION_LABEL = "Basic Tree Information";
    static final String RIGHT_LINE_UP_DOMAINS = "Right-align Domain Architectures";
    static final String LINE_UP_RENDERABLE_DATA = "Line Up Diagrams (such as Domain Architectures)";
    static final String INFER_ANCESTOR_TAXONOMIES = "Infer Ancestor Taxonomies";
    static final String OBTAIN_SEQUENCE_AND_TAXONOMIC_INFORMATION = "Fetch Sequence & Taxonomic Data";
    JMenuBar _jmenubar;
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
    JMenuItem _write_to_svg_item;
    JMenuItem _write_to_eps_item;
    // tools menu:
    JMenuItem _midpoint_root_item;
    JMenuItem _mad_root_item;
    JMenuItem _color_rank_jmi;
    JMenuItem _clade_bands_jmi;

    JMenuItem _obtain_seq_and_tax_information_jmi;
    JMenuItem _extract_label_data_jmi;
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
    JCheckBoxMenuItem _show_overview_cbmi;
    JCheckBoxMenuItem _show_domain_labels;
    JCheckBoxMenuItem _abbreviate_scientific_names;
    JCheckBoxMenuItem _use_italic_scientific_names_cbmi;
    JCheckBoxMenuItem _outline_fonts_in_vector_export_cbmi;
    JCheckBoxMenuItem _transparent_export_background_cbmi;
    JCheckBoxMenuItem _color_labels_same_as_parent_branch;
    JCheckBoxMenuItem _show_default_node_shapes_internal_cbmi;
    JRadioButtonMenuItem _internal_labels_above_branch_rbmi;
    JRadioButtonMenuItem _internal_labels_right_of_node_rbmi;
    JCheckBoxMenuItem _show_default_node_shapes_external_cbmi;
    JCheckBoxMenuItem _show_default_node_shapes_for_marked_cbmi;
    JCheckBoxMenuItem _show_confidence_stddev_cbmi;
    JCheckBoxMenuItem _show_mad_confidence_cbmi;
    JCheckBoxMenuItem _right_line_up_domains_cbmi;
    JCheckBoxMenuItem _line_up_renderable_data_cbmi;
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
    // search: the case/words/regex/inverse/properties options moved to the left ControlPanel.
    JCheckBoxMenuItem _color_all_found_nodes_when_coloring_subtree_cbmi;
    // type menu:
    JMenu _type_menu;
    JCheckBoxMenuItem _rectangular_type_cbmi;
    JCheckBoxMenuItem _triangular_type_cbmi;
    JCheckBoxMenuItem _curved_type_cbmi;
    JCheckBoxMenuItem _convex_type_cbmi;
    JCheckBoxMenuItem _euro_type_cbmi;
    JCheckBoxMenuItem _rounded_type_cbmi;
    JCheckBoxMenuItem _unrooted_type_cbmi;
    JCheckBoxMenuItem _circular_type_cbmi;
    // view as text menu:
    JMenuItem _view_as_NH_item;
    JMenuItem _view_as_XML_item;
    JMenuItem _view_as_nexus_item;
    JMenuItem _display_basic_information_item;
    // help menu:
    JMenuItem _about_item;
    JMenuItem _help_item;
    JMenuItem _website_item;
    JMenuItem _aptxjs_website_item;
    JMenuItem _phyloxml_ref_item;

    //
    File _current_dir;
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
        } else if (o == _clade_bands_jmi) {
            labelCladesByRank();
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
        } else if (o == _show_domain_labels) {
            updateOptions(getOptions());
        } else if (o == _abbreviate_scientific_names) {
            updateOptions(getOptions());
        } else if (o == _use_italic_scientific_names_cbmi) {
            updateOptions(getOptions());
        } else if (o == _outline_fonts_in_vector_export_cbmi) {
            updateOptions(getOptions());
        } else if (o == _transparent_export_background_cbmi) {
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
        } else if (o == _color_all_found_nodes_when_coloring_subtree_cbmi) {
            updateOptions(getOptions());
        } else if (o == _parse_beast_style_extended_nexus_tags_cbmi) {
            updateOptions(getOptions());
        } else if (o == _show_scale_cbmi) {
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
        } else if (o == _line_up_renderable_data_cbmi) {
            if (!_line_up_renderable_data_cbmi.isSelected()) {
                _right_line_up_domains_cbmi.setSelected(false);
            }
            updateOptions(getOptions());
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
        } else if (o == _right_line_up_domains_cbmi) {
            if (_right_line_up_domains_cbmi.isSelected()) {
                _line_up_renderable_data_cbmi.setSelected(true);
            }
            updateOptions(getOptions());
        } else if ((o == _rectangular_type_cbmi) || (o == _triangular_type_cbmi) || (o == _curved_type_cbmi)
                || (o == _convex_type_cbmi) || (o == _euro_type_cbmi) || (o == _rounded_type_cbmi)
                || (o == _unrooted_type_cbmi) || (o == _circular_type_cbmi)) {
            typeChanged(o);
        } else if (o == _about_item) {
            about();
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
        } else if (o == _phyloxml_ref_item) {
            try {
                AptxUtil.openWebsite(AptxConstants.PHYLOXML_REFERENCE_URL);
            } catch (final IOException e1) {
                ForesterUtil.printErrorMessage(AptxConstants.PRG_NAME, e1.toString());
            }
        } else if (o == _write_to_pdf_item) {
            final File curr_dir = writeToPdf(_mainpanel.getCurrentPhylogeny(),
                    getMainPanel(),
                    _writetopdf_filechooser,
                    _current_dir,
                    getContentPane(),
                    this);
            if (curr_dir != null) {
                setCurrentDir(curr_dir);
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
                    _current_dir);
            if (new_dir != null) {
                setCurrentDir(new_dir);
            }
        } else if (o == _write_to_tif_item) {
            final File new_dir = writeToGraphicsFile(_mainpanel.getCurrentPhylogeny(),
                    GraphicsExportType.TIFF,
                    _mainpanel,
                    _writetographics_filechooser,
                    this,
                    getContentPane(),
                    _current_dir);
            if (new_dir != null) {
                setCurrentDir(new_dir);
            }
        } else if (o == _write_to_png_item) {
            final File new_dir = writeToGraphicsFile(_mainpanel.getCurrentPhylogeny(),
                    GraphicsExportType.PNG,
                    _mainpanel,
                    _writetographics_filechooser,
                    this,
                    getContentPane(),
                    _current_dir);
            if (new_dir != null) {
                setCurrentDir(new_dir);
            }
        } else if (o == _write_to_svg_item) {
            final File new_dir = writeToGraphicsFile(_mainpanel.getCurrentPhylogeny(),
                    GraphicsExportType.SVG,
                    _mainpanel,
                    _writetographics_filechooser,
                    this,
                    getContentPane(),
                    _current_dir);
            if (new_dir != null) {
                setCurrentDir(new_dir);
            }
        } else if (o == _write_to_eps_item) {
            final File new_dir = writeToGraphicsFile(_mainpanel.getCurrentPhylogeny(),
                    GraphicsExportType.EPS,
                    _mainpanel,
                    _writetographics_filechooser,
                    this,
                    getContentPane(),
                    _current_dir);
            if (new_dir != null) {
                setCurrentDir(new_dir);
            }
        }  else if (o == _save_item) {
            final File new_dir = writeToFile(_mainpanel.getCurrentPhylogeny(),
                    getMainPanel(),
                    _save_filechooser,
                    _current_dir,
                    getContentPane(),
                    this);
            if (new_dir != null) {
                setCurrentDir(new_dir);
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

    private void deleteSelectedNodes(final boolean delete) {
        final String function = delete ? "Delete" : "Retain";
        final Phylogeny phy = getMainPanel().getCurrentPhylogeny();
        final int ext = (phy == null) ? 0 : phy.getNumberOfExternalNodes();
        final List<PhylogenyNode> nodes = new ArrayList<>();
        if ((phy != null) && (getCurrentTreePanel() != null)
                && ((getCurrentTreePanel().getFoundNodes0() != null)
                        || (getCurrentTreePanel().getFoundNodes1() != null))) {
            for (final PhylogenyNode n : getCurrentTreePanel().getFoundNodesAsListOfPhylogenyNodes()) {
                if (n.isExternal()) {
                    nodes.add(n);
                }
            }
        }
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
            AptxUtil.removeBranchColors(getMainPanel().getCurrentPhylogeny());
            if (getMainPanel().getCurrentTreePanel() != null) {
                getMainPanel().getCurrentTreePanel().clearRankLegend(); // the rank legend no longer reflects the tree
            }
        }
    }

    private void removeVisualStyles() {
        if (getMainPanel().getCurrentPhylogeny() != null) {
            AptxUtil.removeVisualStyles(getMainPanel().getCurrentPhylogeny());
        }
    }

    private void writeAllToFile() {
        if ((getMainPanel().getTabbedPane() == null) || (getMainPanel().getTabbedPane().getTabCount() < 1)) {
            return;
        }
        final File my_dir = getCurrentDir();
        if (my_dir != null) {
            _save_filechooser.setCurrentDirectory(my_dir);
        }
        _save_filechooser.setSelectedFile(new File(""));
        final int result = _save_filechooser.showSaveDialog(_contentpane);
        final File file = _save_filechooser.getSelectedFile();
        setCurrentDir(_save_filechooser.getCurrentDirectory());
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
        if ((getMainPanel().getTabbedPane() != null) && (getMainPanel().getTabbedPane().getTabCount() > 1)) {
            _save_all_item.setEnabled(true);
        } else {
            _save_all_item.setEnabled(false);
        }
    }

    void buildFileMenu() {
        _file_jmenu = MainFrame.createMenu("File", getConfiguration());
        _file_jmenu.add(_save_item = new JMenuItem("Save Tree As..."));
        _file_jmenu.addSeparator();
        _file_jmenu.add(_write_to_pdf_item = new JMenuItem("Export to PDF file ..."));
        _file_jmenu.add(_write_to_svg_item = new JMenuItem("Export to SVG file..."));
        _write_to_svg_item.setToolTipText("Scalable vector graphics for publication (edit in Illustrator/Inkscape)");
        _file_jmenu.add(_write_to_eps_item = new JMenuItem("Export to EPS file..."));
        _write_to_eps_item.setToolTipText("Encapsulated PostScript vector graphics for publication");
        if (AptxUtil.canWriteFormat("tif") || AptxUtil.canWriteFormat("tiff") || AptxUtil.canWriteFormat("TIF")) {
            _file_jmenu.add(_write_to_tif_item = new JMenuItem("Export to TIFF file..."));
        }
        _file_jmenu.add(_write_to_png_item = new JMenuItem("Export to PNG file..."));
        _file_jmenu.add(_write_to_jpg_item = new JMenuItem("Export to JPG file..."));
        _file_jmenu.addSeparator();

        _file_jmenu.add(_exit_item = new JMenuItem("Exit"));
        customizeJMenuItem(_save_item);
        customizeJMenuItem(_write_to_pdf_item);
        customizeJMenuItem(_write_to_svg_item);
        customizeJMenuItem(_write_to_eps_item);
        customizeJMenuItem(_write_to_png_item);
        customizeJMenuItem(_write_to_jpg_item);
        customizeJMenuItem(_write_to_tif_item);
        customizeJMenuItem(_exit_item);
        _jmenubar.add(_file_jmenu);
    }

    void buildHelpMenu() {
        _help_jmenu = createMenu("Help", getConfiguration());
        _help_jmenu.setToolTipText("Documentation, web links, and program information");
        _help_jmenu.add(_help_item = new JMenuItem("Documentation"));
        _help_jmenu.addSeparator();
        _help_jmenu.add(_website_item = new JMenuItem("Archaeopteryx Home"));
        _help_jmenu.add(_aptxjs_website_item = new JMenuItem("Archaeopteryx online version: Archaeopteryx.js"));
        _help_jmenu.add(_phyloxml_ref_item = new JMenuItem("phyloXML Reference"));
        _help_jmenu.addSeparator();
        _help_jmenu.add(_about_item = new JMenuItem("About"));
        customizeJMenuItem(_help_item);
        customizeJMenuItem(_website_item);
        customizeJMenuItem(_aptxjs_website_item);
        customizeJMenuItem(_phyloxml_ref_item);
        customizeJMenuItem(_about_item);
        _phyloxml_ref_item.setToolTipText(PHYLOXML_REF_TOOL_TIP);
        _jmenubar.add(_help_jmenu);
    }

    void buildTypeMenu() {
        _type_menu = createMenu(TYPE_MENU_HEADER, getConfiguration());
        _type_menu.add(_rectangular_type_cbmi = new JCheckBoxMenuItem(MainFrame.RECTANGULAR_TYPE_CBMI_LABEL));
        _type_menu.add(_euro_type_cbmi = new JCheckBoxMenuItem(MainFrame.EURO_TYPE_CBMI_LABEL));
        _type_menu.add(_rounded_type_cbmi = new JCheckBoxMenuItem(MainFrame.ROUNDED_TYPE_CBMI_LABEL));
        _type_menu.add(_curved_type_cbmi = new JCheckBoxMenuItem(MainFrame.CURVED_TYPE_CBMI_LABEL));
        _type_menu.add(_triangular_type_cbmi = new JCheckBoxMenuItem(MainFrame.TRIANGULAR_TYPE_CBMI_LABEL));
        _type_menu.add(_convex_type_cbmi = new JCheckBoxMenuItem(MainFrame.CONVEX_TYPE_CBMI_LABEL));
        _type_menu.add(_unrooted_type_cbmi = new JCheckBoxMenuItem(MainFrame.UNROOTED_TYPE_CBMI_LABEL));
        _type_menu.add(_circular_type_cbmi = new JCheckBoxMenuItem(MainFrame.CIRCULAR_TYPE_CBMI_LABEL));
        customizeCheckBoxMenuItem(_rectangular_type_cbmi, false);
        customizeCheckBoxMenuItem(_triangular_type_cbmi, false);
        customizeCheckBoxMenuItem(_euro_type_cbmi, false);
        customizeCheckBoxMenuItem(_rounded_type_cbmi, false);
        customizeCheckBoxMenuItem(_curved_type_cbmi, false);
        customizeCheckBoxMenuItem(_convex_type_cbmi, false);
        customizeCheckBoxMenuItem(_unrooted_type_cbmi, false);
        customizeCheckBoxMenuItem(_circular_type_cbmi, false);
        _triangular_type_cbmi.setToolTipText("not suitable for phylograms");
        _curved_type_cbmi.setToolTipText("not suitable for phylograms");
        _unrooted_type_cbmi.setToolTipText(MainFrame.USE_MOUSEWHEEL_SHIFT_TO_ROTATE);
        _circular_type_cbmi.setToolTipText(MainFrame.USE_MOUSEWHEEL_SHIFT_TO_ROTATE);
        initializeTypeMenu(getOptions());
        // _type_menu is not added to the menu bar; its items are folded into the Settings dialog.
    }

    void buildViewMenu() {
        _view_jmenu = createMenu("View", getConfiguration());
        _view_jmenu.setToolTipText("Show tree information, or the tree as phyloXML, Newick, or Nexus");
        _view_jmenu.add(_display_basic_information_item = new JMenuItem(SHOW_BASIC_TREE_INFORMATION_LABEL));
        _view_jmenu.addSeparator();
        _view_jmenu.add(_view_as_XML_item = new JMenuItem("as phyloXML"));
        _view_jmenu.add(_view_as_NH_item = new JMenuItem("as Newick"));
        _view_jmenu.add(_view_as_nexus_item = new JMenuItem("as Nexus"));
        customizeJMenuItem(_display_basic_information_item);
        customizeJMenuItem(_view_as_NH_item);
        customizeJMenuItem(_view_as_XML_item);
        customizeJMenuItem(_view_as_nexus_item);
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
                    unresolved.size() + " tip " + ((unresolved.size() == 1) ? "taxon" : "taxa")
                            + " could not be placed at rank \"" + r + "\" from the tree's own data.\n"
                            + "Resolve online via the NCBI and UniProt databases? (requires an internet connection)",
                    "Resolve Taxa Online?",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.QUESTION_MESSAGE);
            if (choice == JOptionPane.YES_OPTION) {
                new Thread(new OnlineTaxonResolver(this, "rank colorization (" + r + ")", unresolved, err -> {
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
        tp.reportRankColorization(r, tp.colorByRank(r));
    }

    /**
     * The Tools "Annotate Clades by Rank…" operation: pick a rank and a mode (shaded boxes or right-edge
     * bars), resolve any unplaced tips online if the user agrees (off the EDT), then draw the bands.
     * Mirrors {@link #colorRank()} but renders {@link CladeBand}s instead of coloring branches.
     */
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
        final JPanel panel = new JPanel(new GridLayout(0, 1, 0, 2));
        panel.add(new JLabel("Annotate clades by rank:"));
        panel.add(rank_box);
        panel.add(new JLabel(" "));
        panel.add(new JLabel("Show as:"));
        panel.add(boxes_rb);
        panel.add(bars_rb);
        panel.add(brackets_rb);
        if (JOptionPane.showConfirmDialog(this, panel, "Annotate Clades by Rank", JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.PLAIN_MESSAGE) != JOptionPane.OK_OPTION) {
            return;
        }
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
                    unresolved.size() + " tip " + ((unresolved.size() == 1) ? "taxon" : "taxa")
                            + " could not be placed at rank \"" + r + "\" from the tree's own data.\n"
                            + "Resolve online via the NCBI and UniProt databases? (requires an internet connection)",
                    "Resolve Taxa Online?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
            if (choice == JOptionPane.YES_OPTION) {
                new Thread(new OnlineTaxonResolver(this, "clade bands (" + r + ")", unresolved,
                        err -> reportCladeBands(tp, r, mode, err))).start();
                return;
            }
        }
        reportCladeBands(tp, r, mode, null);
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
                                  final String error) {
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
        if (error != null) {
            JOptionPane.showMessageDialog(this, "Drew " + n + " clade " + kind + " at rank \"" + rank
                    + "\", but some taxa could not be resolved:\n" + error, "Annotate Clades by Rank (" + rank + ")",
                    JOptionPane.WARNING_MESSAGE);
        } else if (n > 0) {
            JOptionPane.showMessageDialog(this, "Drew " + n + " clade " + kind + " at rank \"" + rank + "\".",
                    "Annotate Clades by Rank (" + rank + ")", JOptionPane.INFORMATION_MESSAGE);
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
     * Load-time offer (NH/NHX/Nexus only, and only when the manual "Internal Node Names are Confidence
     * Values" option is off): if the internal node labels look like support values
     * ({@link AptxUtil#internalNamesLookLikeConfidenceValues}), ask whether to treat them as confidence
     * values rather than node names.
     */
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
                _mainpanel.getControlPanel().setCheckbox(Configuration.write_confidence_values, true);
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
        final int n = LabelDataExtractor.extract(phy);
        if (n > 0) {
            tp.setEdited(true);
            // surface the freshly-populated fields so the value of the cleanup is immediately visible
            tp.getControlPanel().setCheckbox(Configuration.show_seq_names, true);
            tp.getControlPanel().setCheckbox(Configuration.show_taxonomy_scientific_names, true);
            tp.getControlPanel().displayedPhylogenyMightHaveChanged(true);
            tp.repaint();
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


    void customizeCheckBoxMenuItem(final JCheckBoxMenuItem item, final boolean is_selected) {
        if (item != null) {
            item.setFont(MainFrame.menu_font);
            if (getConfiguration().isApplyCustomGuiColors()) {
                item.setBackground(getConfiguration().getGuiMenuBackgroundColor());
                item.setForeground(getConfiguration().getGuiMenuTextColor());
            }
            item.setSelected(is_selected);
            item.addActionListener(this);
        }
    }

    JMenuItem customizeJMenuItem(final JMenuItem jmi) {
        if (jmi != null) {
            jmi.setFont(MainFrame.menu_font);
            if (getConfiguration().isApplyCustomGuiColors()) {
                jmi.setBackground(getConfiguration().getGuiMenuBackgroundColor());
                jmi.setForeground(getConfiguration().getGuiMenuTextColor());
            }
            jmi.addActionListener(this);
        }
        return jmi;
    }

    void customizeRadioButtonMenuItem(final JRadioButtonMenuItem item, final boolean is_selected) {
        if (item != null) {
            item.setFont(MainFrame.menu_font);
            if (getConfiguration().isApplyCustomGuiColors()) {
                item.setBackground(getConfiguration().getGuiMenuBackgroundColor());
                item.setForeground(getConfiguration().getGuiMenuTextColor());
            }
            item.setSelected(is_selected);
            item.addActionListener(this);
        }
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
        final AncestralTaxonomyInferrer inferrer = new AncestralTaxonomyInferrer(this,
                _mainpanel.getCurrentTreePanel(),
                _mainpanel.getCurrentPhylogeny()
                        .copy());
        new Thread(inferrer).start();
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

    File getCurrentDir() {
        if ((_current_dir == null) || !_current_dir.canRead()) {
            if (ForesterUtil.isWindows()) {
                try {
                    _current_dir = new File(WindowsUtils.getCurrentUserDesktopPath());
                } catch (final Exception e) {
                    _current_dir = null;
                }
            }
        }
        if ((_current_dir == null) || !_current_dir.canRead()) {
            if (System.getProperty("user.home") != null) {
                _current_dir = new File(System.getProperty("user.home"));
            } else if (System.getProperty("user.dir") != null) {
                _current_dir = new File(System.getProperty("user.dir"));
            }
        }
        return _current_dir;
    }

    TreePanel getCurrentTreePanel() {
        return getMainPanel().getCurrentTreePanel();
    }

    JMenu getHelpMenu() {
        return _help_jmenu;
    }

    JCheckBoxMenuItem getlabelDirectionCbmi() {
        return _label_direction_cbmi;
    }

    final Phylogeny getSpeciesTree() {
        return _species_tree;
    }

    void initializeTypeMenu(final Options options) {
        setTypeMenuToAllUnselected();
        switch (options.getPhylogenyGraphicsType()) {
            case CONVEX:
                _convex_type_cbmi.setSelected(true);
                break;
            case CURVED:
                _curved_type_cbmi.setSelected(true);
                break;
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

    void setCurrentDir(final File current_dir) {
        _current_dir = current_dir;
    }

    void setOptions(final Options options) {
        _options = options;
    }

    void setSelectedTypeInTypeMenu(final PHYLOGENY_GRAPHICS_TYPE type) {
        setTypeMenuToAllUnselected();
        switch (type) {
            case CIRCULAR:
                _circular_type_cbmi.setSelected(true);
                break;
            case CONVEX:
                _convex_type_cbmi.setSelected(true);
                break;
            case CURVED:
                _curved_type_cbmi.setSelected(true);
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
        _convex_type_cbmi.setSelected(false);
        _curved_type_cbmi.setSelected(false);
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
            if (((previous_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && (new_type != PHYLOGENY_GRAPHICS_TYPE.UNROOTED))
                    || ((previous_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) && (new_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR))
                    || ((previous_type != PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && (new_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED))
                    || ((previous_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) && (new_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR))) {
                getCurrentTreePanel().getControlPanel().showWhole();
            }
            if (getCurrentTreePanel().isPhyHasBranchLengths() && (new_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)) {
                getCurrentTreePanel().getControlPanel().setDrawPhylogramEnabled(true);
            } else {
                getCurrentTreePanel().getControlPanel().setDrawPhylogramEnabled(false);
            }
            getCurrentTreePanel().setPhylogenyGraphicsType(getOptions().getPhylogenyGraphicsType());
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

    void updateOptions(final Options options) {
        options.setShowDomainLabels((_show_domain_labels != null) && _show_domain_labels.isSelected());
        options.setAbbreviateScientificTaxonNames((_abbreviate_scientific_names != null)
                && _abbreviate_scientific_names.isSelected());
        options.setUseItalicScientificNames((_use_italic_scientific_names_cbmi != null)
                && _use_italic_scientific_names_cbmi.isSelected());
        options.setOutlineFontsInVectorExport((_outline_fonts_in_vector_export_cbmi != null)
                && _outline_fonts_in_vector_export_cbmi.isSelected());
        options.setTransparentExportBackground((_transparent_export_background_cbmi != null)
                && _transparent_export_background_cbmi.isSelected());
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
        if ((_show_scale_cbmi != null) && _show_scale_cbmi.isEnabled()) {
            options.setShowScale(_show_scale_cbmi.isSelected());
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
        } else if ((_curved_type_cbmi != null) && _curved_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.CURVED);
        } else if ((_convex_type_cbmi != null) && _convex_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.CONVEX);
        } else if ((_euro_type_cbmi != null) && _euro_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE);
        } else if ((_rounded_type_cbmi != null) && _rounded_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.ROUNDED);
        } else if ((_unrooted_type_cbmi != null) && _unrooted_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.UNROOTED);
        } else if ((_circular_type_cbmi != null) && _circular_type_cbmi.isSelected()) {
            options.setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.CIRCULAR);
        }
        if ((_right_line_up_domains_cbmi != null) && _right_line_up_domains_cbmi.isEnabled()) {
            options.setRightLineUpDomains(_right_line_up_domains_cbmi.isSelected());
        }
        if ((_line_up_renderable_data_cbmi != null) && _line_up_renderable_data_cbmi.isEnabled()) {
            options.setLineUpRendarableNodeData(_line_up_renderable_data_cbmi.isSelected());
        }
        if ((_color_all_found_nodes_when_coloring_subtree_cbmi != null) && _color_all_found_nodes_when_coloring_subtree_cbmi.isEnabled()) {
            options.setColorAllFoundNodesWhenColoringSubtree(_color_all_found_nodes_when_coloring_subtree_cbmi.isSelected());
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
        final StringBuffer about = new StringBuffer("Archaeopteryx\nVersion " + AptxConstants.VERSION + "\n");
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
        about.append("For more information & download:\n");
        about.append(AptxConstants.APTX_WEB_SITE + "\n");
        about.append("Documentation:\n");
        about.append(AptxConstants.APTX_DOC_SITE + "\n");
        about.append("Comments: " + AptxConstants.AUTHOR_EMAIL);
        JOptionPane.showMessageDialog(null, about, AptxConstants.PRG_NAME, JOptionPane.PLAIN_MESSAGE);
    }

    static JMenu createMenu(final String title, final Configuration conf) {
        final JMenu jmenu = new JMenu(title);
        if (conf.isApplyCustomGuiColors()) {
            jmenu.setFont(MainFrame.menu_font);
            jmenu.setBackground(conf.getGuiMenuBackgroundColor());
            jmenu.setForeground(conf.getGuiMenuTextColor());
        }
        return jmenu;
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



    static void printPhylogenyToPdf(final String file_name,
                                    final Options opts,
                                    final TreePanel tp,
                                    final Component comp) {

        String pdf_written_to = "";
        boolean error = false;
        try {
            if (opts.isPrintUsingActualSize()) {
                pdf_written_to = PdfExporter.writePhylogenyToPdf(file_name, tp, tp.getWidth(), tp.getHeight());
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
            printPhylogenyToPdf(file.toString(), mp.getOptions(), mp.getCurrentTreePanel(), component);
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
                || file_name.endsWith(".con") || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "All supported files (*.xml, *.phyloxml, *phylo.xml, *.nhx, *.nh, *.newick, *.nex, *.nexus, *.phy, *.tre, *.tree, *.tol, ...)";
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
