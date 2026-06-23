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
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.prefs.Preferences;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.OVERVIEW_PLACEMENT_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.phylogeny.data.NodeDataField;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.util.ForesterUtil;

public final class Configuration {

    public enum EXT_NODE_DATA_RETURN_ON {
        BUFFER_ONLY, CONSOLE, WINODW
    }

    public enum UI {
        CROSSPLATFORM, NATIVE, NIMBUS, FLAT_LIGHT, FLAT_DARK, UNKNOWN
    }

    private static final String PREFS_NODE = "org/forester/archaeopteryx";
    private static final String PREFS_KEY_UI = "ui";

    final static String[][] clickto_options = {
            {"Display Node Data", "display"}, {"Collapse/Uncollapse", "display"}, {"Root/Reroot", "display"},
            {"Go to Sub/Supertree", "display"}, {"Swap Descendants", "display"},
            {"Colorize Node(s)", "display"},
            {"Colorize Subtree(s)", "display"}, {"Open Sequence DB", "display"}, {"Open PDB", "display"},
            {"Open Taxonomy DB", "display"}, {"Launch BLAST", "display"}, {"Cut Subtree", "display"},
            {"Copy Subtree", "display"}, {"Paste Subtree", "display"}, {"Delete Subtree/Node", "display"},
            {"Add New Node", "display"}, {"Edit Node Data", "display"}, {"Sort Descendants", "display"},
            {"List Node Data", "display"}, {"Select Node(s)", "display"}, {"Uncollapse All", "display"}, {"Order Subtree", "display"},};

    final static int display_node_data = 0;
    final static int collapse_uncollapse = 1;
    final static int reroot = 2;
    final static int subtree = 3;
    final static int swap = 4;
    final static int color_node_font = 5;
    final static int color_subtree = 6;
    final static int open_seq_web = 7;
    final static int open_pdb_web = 8;
    final static int open_tax_web = 9;
    final static int blast = 10;
    final static int cut_subtree = 11;
    final static int copy_subtree = 12;
    final static int paste_subtree = 13;
    final static int delete_subtree_or_node = 14;
    final static int add_new_node = 15;
    final static int edit_node_data = 16;
    final static int sort_descendents = 17;
    final static int get_ext_desc_data = 18;
    final static int select_nodes = 19;
    final static int uncollapse_all = 20;
    final static int order_subtree = 21;

    // ------------------
    // Click-to options
    // ------------------
    final static String[][] display_options = {
            {"Phylogram", "display", "?"}, {"Node Name", "display", "yes"}, {"Taxonomy Code", "display", "yes"},
            {"Seq Annotations", "display", "no"}, {"Confidence Values", "display", "?"},
            {"Node Events", "display", "?"}, {"Colorize by Taxonomy", "display", "no"},
            {"Colorize by Sequence", "display", "no"}, {"Visual Styles/Branch Colors", "display", "yes"},
            {"Branch Widths", "display", "no"}, {"Domain Architectures", "display", "no"},
            {"Binary Characters", "nodisplay", "no"}, {"Binary Char Counts", "nodisplay", "no"},
            {"Seq Name", "display", "no"}, {"Seq Accession", "display", "no"},
            {"Show Internal Data", "display", "yes"}, {"Dyna Hide", "display", "yes"},
            {"Taxonomy Scientific", "display", "yes"}, {"Taxonomy Common", "display", "no"},
            {"Colorize by Annotation", "nodisplay", "no"}, {"Seq Symbol", "nodisplay", "no"},
            {"Rollover", "display", "yes"}, {"Relation Confidence", "nodisplay", "no"},
            {"Vector Data", "nodisplay", "no"}, {"Taxonomy Images", "nodisplay", "no"},
            {"Properties", "display", "no"}, {"Gene Name", "nodisplay", "no"},
            {"Multiple Seq Alignment", "nodisplay", "no"}, {"Branch Length Values", "display", "no"}
            , {"Taxonomy Rank", "display", "no"}, {"Show External Data", "display", "yes"},
            {"Shorten Labels", "display", "yes"}};
    final static int display_as_phylogram = 0;
    final static int show_node_names = 1;
    final static int show_tax_code = 2;
    final static int write_confidence_values = 4;
    final static int write_events = 5;
    final static int use_style = 8;
    final static int width_branches = 9;
    final static int show_domain_architectures = 10;
    final static int show_binary_characters = 11;
    final static int show_binary_character_counts = 12;
    final static int show_seq_names = 13;
    final static int show_sequence_acc = 14;
    final static int display_internal_data = 15;
    final static int dynamically_hide_data = 16;
    final static int show_taxonomy_scientific_names = 17;
    final static int show_taxonomy_common_names = 18;
    final static int show_seq_symbols = 20;
    final static int node_data_popup = 21;
    final static int show_relation_confidence = 22;
    final static int show_vector_data = 23;
    final static int show_taxonomy_images = 24;
    final static int show_properties = 25;
    final static int show_gene_names = 26;
    final static int show_mol_seqs = 27;
    final static int write_branch_length_values = 28;
    final static int show_tax_rank = 29;
    final static int display_external_data = 30;
    final static int shorten_labels = 31;

    private static Hashtable<String, Color> _domain_colors;
    private static Hashtable<String, Color> _species_colors;
    private static String DEFAULT_FONT_FAMILY = "";

    // ---------------------------
    // Display options for trees
    // ---------------------------
    // ---------------------------------
    // This option is selected in the dropdown
    int default_clickto = Configuration.display_node_data;

    private final boolean _abbreviate_scientific_names = false;
    private String _base_font_family_name = "";
    private int _base_font_size = -1;
    private final CLADOGRAM_TYPE _cladogram_type = AptxConstants.CLADOGRAM_TYPE_DEFAULT;
    private boolean _color_labels_same_as_parent_branch = false;
    private NodeFill _default_node_fill = NodeFill.SOLID;
    private NodeShape _default_node_shape = NodeShape.RECTANGLE;
    private short _default_node_shape_size = AptxConstants.DEFAULT_NODE_SHAPE_SIZE_DEFAULT;
    private SortedMap<String, Color> _display_colors = null;
    private boolean _display_sequence_relations = false;

    private NodeDataField _ext_desc_data_to_return = NodeDataField.UNKNOWN;
    private final EXT_NODE_DATA_RETURN_ON _ext_node_data_return_on = EXT_NODE_DATA_RETURN_ON.WINODW;
    private final Color _gui_background_color = AptxConstants.GUI_BACKGROUND_DEFAULT;
    private final Color _gui_button_background_color = AptxConstants.BUTTON_BACKGROUND_COLOR_DEFAULT;
    private final Color _gui_button_border_color = AptxConstants.BUTTON_BORDER_COLOR_DEFAULT;
    private final Color _gui_button_text_color = AptxConstants.BUTTON_TEXT_COLOR_DEFAULT;
    private final Color _gui_checkbox_and_button_active_color = AptxConstants.CHECKBOX_AND_BUTTON_ACTIVE_COLOR_DEFAULT;
    private final Color _gui_checkbox_text_color = AptxConstants.CHECKBOX_TEXT_COLOR_DEFAULT;
    private final Color _gui_menu_background_color = AptxConstants.MENU_BACKGROUND_COLOR_DEFAULT;
    private final Color _gui_menu_text_color = AptxConstants.MENU_TEXT_COLOR_DEFAULT;
    private final boolean _hide_controls_and_menus = false;
    private final boolean _internal_number_are_confidence_for_nh_parsing = false;
    private String _label_for_get_ext_descendents_data = "";
    private boolean _midpoint_root = false;
    private int _min_base_font_size = 2;
    private double _min_confidence_value = Options.MIN_CONFIDENCE_DEFAULT;
    private boolean _nh_parsing_replace_underscores = false;
    private NODE_LABEL_DIRECTION _node_label_direction = NODE_LABEL_DIRECTION.HORIZONTAL;
    private short _number_of_digits_after_comma_for_branch_length_values = AptxConstants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_BRANCH_LENGTH_VALUES_DEFAULT;
    private short _number_of_digits_after_comma_for_confidence_values = AptxConstants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_CONFIDENCE_VALUES_DEFAULT;
    private short _ov_max_height = 80;
    private short _ov_max_width = 80;
    private OVERVIEW_PLACEMENT_TYPE _ov_placement = OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT;
    private PHYLOGENY_GRAPHICS_TYPE _phylogeny_graphics_type = PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR;
    private float _print_line_width = AptxConstants.PDF_LINE_WIDTH_DEFAULT;
    private boolean _show_default_node_shapes_external = false;
    private boolean _show_default_node_shapes_for_marked_nodes = false;
    private boolean _show_default_node_shapes_internal = false;
    private boolean _show_domain_labels = true;
    private boolean _show_overview = true;
    private boolean _show_scale = false;
    private TAXONOMY_EXTRACTION _taxonomy_extraction = TAXONOMY_EXTRACTION.NO;
    private UI _ui = UI.UNKNOWN;
    private boolean _use_tabbed_display = false;
    private boolean _validate_against_phyloxml_xsd_schema = AptxConstants.VALIDATE_AGAINST_PHYLOXML_XSD_SCJEMA_DEFAULT;
    private Color _vector_data_min_color = Color.BLUE;
    private Color _vector_data_max_color = Color.YELLOW;
    private Color _vector_data_mean_color = Color.WHITE;
    private double _vector_data_height = 12;
    private int _vector_data_width = 120;
    private boolean _line_up_renderable_node_data = true;
    private boolean _right_align_domains = false;

    static {
        for (final String font_name : AptxConstants.DEFAULT_FONT_CHOICES) {
            if (Arrays.binarySearch(AptxUtil.getAvailableFontFamiliesSorted(), font_name) >= 0) {
                DEFAULT_FONT_FAMILY = font_name;
                break;
            }
        }
        if (ForesterUtil.isEmpty(DEFAULT_FONT_FAMILY)) {
            DEFAULT_FONT_FAMILY = AptxConstants.DEFAULT_FONT_CHOICES[AptxConstants.DEFAULT_FONT_CHOICES.length - 1];
        }
    }

    public Configuration() {
        // Archaeopteryx no longer reads configuration files; all settings come from the
        // built-in defaults (see the field initializers and display_options above) and the
        // Settings dialog at runtime.
        setDisplayColors(new TreeMap<String, Color>());
    }

    public String getBaseFontFamilyName() {
        return _base_font_family_name;
    }

    public NodeFill getDefaultNodeFill() {
        return _default_node_fill;
    }

    public NodeShape getDefaultNodeShape() {
        return _default_node_shape;
    }

    public short getDefaultNodeShapeSize() {
        return _default_node_shape_size;
    }

    public NodeDataField getExtDescNodeDataToReturn() {
        return _ext_desc_data_to_return;
    }

    public EXT_NODE_DATA_RETURN_ON getExtNodeDataReturnOn() {
        return _ext_node_data_return_on;
    }


    public String getLabelForGetExtDescendentsData() {
        return _label_for_get_ext_descendents_data;
    }

    public double getVectorDataHeight() {
        return _vector_data_height;
    }

    public Color getVectorDataMaxColor() {
        return _vector_data_max_color;
    }

    public Color getVectorDataMeanColor() {
        return _vector_data_mean_color;
    }

    public Color getVectorDataMinColor() {
        return _vector_data_min_color;
    }

    public int getVectorDataWidth() {
        return _vector_data_width;
    }

    public boolean isAbbreviateScientificTaxonNames() {
        return _abbreviate_scientific_names;
    }

    public boolean isColorByTaxonomicGroup() {
        return false;
    }

    public boolean isColorLabelsSameAsParentBranch() {
        return _color_labels_same_as_parent_branch;
    }

    final public boolean isLineUpRendarableNodeData() {
        return _line_up_renderable_node_data;
    }

    public boolean isMidpointReroot() {
        return _midpoint_root;
    }

    final public boolean isRightLineUpDomains() {
        return _right_align_domains;
    }

    public boolean isShowDefaultNodeShapesExternal() {
        return _show_default_node_shapes_external;
    }

    public boolean isShowDefaultNodeShapesForMarkedNodes() {
        return _show_default_node_shapes_for_marked_nodes;
    }

    public boolean isShowDefaultNodeShapesInternal() {
        return _show_default_node_shapes_internal;
    }

    public boolean isShowDomainLabels() {
        return _show_domain_labels;
    }

    public void putDisplayColors(final String key, final Color color) {
        getDisplayColors().put(key, color);
    }

    public void setAddTaxonomyImagesCB(final boolean b) {
        display_options[show_taxonomy_images][1] = b ? "yes" : "no";
    }

    public void setBaseFontFamilyName(final String base_font_family_name) {
        _base_font_family_name = base_font_family_name;
    }

    public void setBaseFontSize(final int base_font_size) {
        _base_font_size = base_font_size;
    }

    public void setDisplayAsPhylogram(final boolean b) {
        display_options[display_as_phylogram][2] = b ? "yes" : "no";
    }

    public void setDisplayColors(final SortedMap<String, Color> display_colors) {
        _display_colors = display_colors;
    }

    public void setDisplayGeneNames(final boolean b) {
        display_options[show_gene_names][2] = b ? "yes" : "no";
    }

    public void setDisplayMultipleSequenceAlignment(final boolean b) {
        display_options[show_mol_seqs][2] = b ? "yes" : "no";
    }

    public void setDisplaySequenceNames(final boolean b) {
        display_options[show_seq_names][2] = b ? "yes" : "no";
    }


    public void setDisplaySequenceSymbols(final boolean b) {
        display_options[show_seq_symbols][2] = b ? "yes" : "no";
    }

    public void setDisplayTaxonomyCode(final boolean b) {
        display_options[show_tax_code][2] = b ? "yes" : "no";
    }


    public void setDisplayTaxonomyCommonNames(final boolean b) {
        display_options[show_taxonomy_common_names][2] = b ? "yes" : "no";
    }

    public void setDisplayTaxonomyScientificNames(final boolean b) {
        display_options[show_taxonomy_scientific_names][2] = b ? "yes" : "no";
    }


    public void setShowScale(final boolean show_scale) {
        _show_scale = show_scale;
    }

    public void setUseStyle(final boolean b) {
        display_options[use_style][2] = b ? "yes" : "no";
    }

    private final void initSpeciesColors() {
        _species_colors = new Hashtable<String, Color>();
    }


    boolean displaySequenceRelations() {
        return _display_sequence_relations;
    }

    boolean doCheckOption(final int which) {
        return (display_options[which][2].equalsIgnoreCase("yes"))
                || (display_options[which][2].equalsIgnoreCase("true"));
    }

    boolean doDisplayClickToOption(final int which) {
        return clickto_options[which][1].equalsIgnoreCase("display");
    }

    boolean doDisplayOption(final int which) {
        return display_options[which][1].equalsIgnoreCase("display");
    }

    /**
     * Will attempt to use the phylogeny to determine whether to check
     * this or not (e.g. phylogram)
     */
    boolean doGuessCheckOption(final int which) {
        return display_options[which][2].equals("?");
    }

    int getBaseFontSize() {
        return _base_font_size;
    }

    CLADOGRAM_TYPE getCladogramType() {
        return _cladogram_type;
    }

    String getClickToTitle(final int which) {
        return clickto_options[which][0];
    }

    int getDefaultDisplayClicktoOption() {
        return default_clickto;
    }

    SortedMap<String, Color> getDisplayColors() {
        return _display_colors;
    }

    String getDisplayTitle(final int which) {
        return display_options[which][0];
    }

    Map<String, Color> getDomainColors() {
        if (_domain_colors == null) {
            _domain_colors = new Hashtable<String, Color>();
        }
        return _domain_colors;
    }

    Color getGuiBackgroundColor() {
        return _gui_background_color;
    }

    Color getGuiButtonBackgroundColor() {
        return _gui_button_background_color;
    }

    Color getGuiButtonBorderColor() {
        return _gui_button_border_color;
    }

    Color getGuiButtonTextColor() {
        return _gui_button_text_color;
    }

    Color getGuiCheckboxAndButtonActiveColor() {
        return _gui_checkbox_and_button_active_color;
    }

    Color getGuiCheckboxTextColor() {
        return _gui_checkbox_text_color;
    }

    Color getGuiMenuBackgroundColor() {
        return _gui_menu_background_color;
    }

    Color getGuiMenuTextColor() {
        return _gui_menu_text_color;
    }

    static int getGuiFontSize() {
        return 11;
    }

    int getMinBaseFontSize() {
        return _min_base_font_size;
    }

    double getMinConfidenceValue() {
        return _min_confidence_value;
    }

    NODE_LABEL_DIRECTION getNodeLabelDirection() {
        return _node_label_direction;
    }

    short getNumberOfDigitsAfterCommaForBranchLengthValues() {
        return _number_of_digits_after_comma_for_branch_length_values;
    }

    short getNumberOfDigitsAfterCommaForConfidenceValues() {
        return _number_of_digits_after_comma_for_confidence_values;
    }

    short getOvMaxHeight() {
        return _ov_max_height;
    }

    short getOvMaxWidth() {
        return _ov_max_width;
    }

    OVERVIEW_PLACEMENT_TYPE getOvPlacement() {
        return _ov_placement;
    }

    PHYLOGENY_GRAPHICS_TYPE getPhylogenyGraphicsType() {
        return _phylogeny_graphics_type;
    }

    float getPrintLineWidth() {
        return _print_line_width;
    }

    Hashtable<String, Color> getSpeciesColors() {
        if (_species_colors == null) {
            initSpeciesColors();
        }
        return _species_colors;
    }

    final TAXONOMY_EXTRACTION getTaxonomyExtraction() {
        return _taxonomy_extraction;
    }

    /**
     * Convenience method.
     *
     * @return true if the tree should be drawn as a phylogram (vs. cladogram) by default
     */
    boolean isDrawAsPhylogram() {
        return doCheckOption(display_as_phylogram);
    }

    boolean isEditable() {
        return true;
    }


    boolean isInternalNumberAreConfidenceForNhParsing() {
        return _internal_number_are_confidence_for_nh_parsing;
    }

    boolean isReplaceUnderscoresInNhParsing() {
        return _nh_parsing_replace_underscores;
    }

    boolean isShowOverview() {
        return _show_overview;
    }

    boolean isShowScale() {
        return _show_scale;
    }

    /**
     * Returns the resolved look-and-feel selection. When no preference has been set yet
     * ({@code UNKNOWN}), the last theme the user chose at runtime (persisted via
     * {@link java.util.prefs.Preferences}) is used; if none was ever saved, the modern
     * FlatLaf light theme is the default.
     */
    final UI getUi() {
        if (_ui == UI.UNKNOWN) {
            _ui = readUiPreference();
        }
        return _ui;
    }

    final void setUi(final UI ui) {
        _ui = ui;
    }

    private static UI readUiPreference() {
        try {
            final String saved = Preferences.userRoot().node(PREFS_NODE).get(PREFS_KEY_UI, null);
            if (saved != null) {
                final UI ui = UI.valueOf(saved);
                if (ui != UI.UNKNOWN) {
                    return ui;
                }
            }
        } catch (final Exception e) {
            // an invalid or inaccessible preference simply falls through to the default
        }
        return UI.FLAT_LIGHT;
    }

    /**
     * Persists the user's runtime look-and-feel choice so it survives a restart.
     */
    static void saveUiPreference(final UI ui) {
        try {
            Preferences.userRoot().node(PREFS_NODE).put(PREFS_KEY_UI, ui.name());
        } catch (final Exception e) {
            // failing to persist the preference is non-fatal
        }
    }

    final boolean isUseNativeUI() {
        return getUi() == UI.NATIVE;
    }

    /**
     * Whether the legacy, hand-themed cross-platform GUI colors and fonts should be
     * applied to the Swing components. The native and FlatLaf look-and-feels style the
     * components themselves, so custom colors are only applied for {@code CROSSPLATFORM}.
     */
    final boolean isApplyCustomGuiColors() {
        return getUi() == UI.CROSSPLATFORM;
    }


    boolean isValidatePhyloXmlAgainstSchema() {
        return _validate_against_phyloxml_xsd_schema;
    }


    static String getDefaultFontFamilyName() {
        return DEFAULT_FONT_FAMILY;
    }

}
