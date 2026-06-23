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

import java.awt.Font;

import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.NodeDataField;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.util.ForesterUtil;

/*
 * This is to hold changeable options.
 */
final public class Options {


    public static enum CLADOGRAM_TYPE {
        LINED_UP, NON_LINED_UP;
    }

    public static enum NODE_LABEL_DIRECTION {
        HORIZONTAL, RADIAL;
    }

    public static enum PHYLOGENY_GRAPHICS_TYPE {
        CIRCULAR, CONVEX, CURVED, EURO_STYLE, RECTANGULAR, ROUNDED, TRIANGULAR, UNROOTED;
    }

    static enum PHYLOGENY_DISPLAY_TYPE {
        CLADOGRAM,
        ALIGNED_PHYLOGRAM,
        UNALIGNED_PHYLOGRAM
    }

    static enum OVERVIEW_PLACEMENT_TYPE {
        LOWER_LEFT("lower left"),
        LOWER_RIGHT("lower right"),
        UPPER_LEFT("upper left"),
        UPPER_RIGHT("upper right");

        private final String _name;

        private OVERVIEW_PLACEMENT_TYPE(final String name) {
            _name = name;
        }

        @Override
        public String toString() {
            return _name;
        }

        public String toTag() {
            return toString().replaceAll(" ", "_");
        }
    }

    // How branch support (confidence) at internal nodes is drawn as a node symbol. Monochrome by
    // design (uses the branch color) so the color channel stays free for clade/taxonomy identity.
    static enum SUPPORT_VISUALIZATION {
        NONE("none"),
        THRESHOLD_MARKS("threshold marks"),
        SIZE_SCALED("size-scaled");

        private final String _name;

        private SUPPORT_VISUALIZATION(final String name) {
            _name = name;
        }

        @Override
        public String toString() {
            return _name;
        }
    }

    static final double MIN_CONFIDENCE_DEFAULT = 50.0;
    // Default cutoff for THRESHOLD_MARKS, as a fraction (0..1) of the support scale -- 0.95 marks
    // the conventional "well-supported" node (>=95% bootstrap or >=0.95 posterior probability).
    static final double SUPPORT_THRESHOLD_DEFAULT = 0.95;
    private boolean _abbreviate_scientific_names;
    private boolean _allow_errors_in_distance_to_parent;
    private boolean _antialias_print;
    private Font _base_font;
    private CLADOGRAM_TYPE _cladogram_type;
    private boolean _color_by_taxonomic_group;
    private boolean _color_labels_same_as_parent_branch;
    private NodeVisualData.NodeFill _default_node_fill;
    private NodeVisualData.NodeShape _default_node_shape;
    private short _default_node_shape_size;
    private boolean _editable;
    private NodeDataField _ext_desc_data_to_return;
    private final boolean _graphics_export_using_actual_size = true;
    private boolean _graphics_export_visible_only;
    private boolean _internal_number_are_confidence_for_nh_parsing;
    private boolean _inverse_search_result;
    private boolean _match_whole_terms_only;
    private boolean _search_with_regex;
    private double _min_confidence_value;
    private SUPPORT_VISUALIZATION _support_visualization;
    private double _support_threshold;
    private NH_CONVERSION_SUPPORT_VALUE_STYLE _nh_conversion_support_value_style;
    private boolean _nh_parsing_replace_underscores;
    private NODE_LABEL_DIRECTION _node_label_direction;
    private short _number_of_digits_after_comma_for_branch_length_values;
    private short _number_of_digits_after_comma_for_confidence_values;
    private OVERVIEW_PLACEMENT_TYPE _ov_placement;
    private PHYLOGENY_GRAPHICS_TYPE _phylogeny_graphics_type;
    private boolean _print_black_and_white;
    private float _print_line_width;
    private final boolean _print_using_actual_size = true;
    private double _scale_bar_length;
    private boolean _search_case_sensitive;
    private boolean _show_confidence_stddev;
    private boolean _show_mad_confidence;
    private boolean _show_default_node_shapes_for_marked_nodes;
    private boolean _show_default_node_shapes_external;
    private boolean _show_default_node_shapes_internal;
    private boolean _internal_labels_above_branch;
    private boolean _use_italic_scientific_names;
    private boolean _outline_fonts_in_vector_export;
    private boolean _show_domain_labels;
    private boolean _show_overview;
    private boolean _show_scale;
    private TAXONOMY_EXTRACTION _taxonomy_extraction;
    private boolean _line_up_renderable_node_data;
    private boolean _right_align_domains;
    private boolean _color_all_found_nodes_when_coloring_subtree;
    private boolean _parse_beast_style_extended_nexus_tags;
    private boolean _collapsed_with_average_height;
    private boolean _show_abbreviated_labels_for_collapsed_nodes;

    private boolean _search_properties;
    private float _default_branch_width;

    private Options() {
        init();
    }


    public void setDefaultBranchWidth(final float default_branch_width) {
        _default_branch_width = default_branch_width;

    }

    public float getDefaultBranchWidth() {
        return _default_branch_width;
    }

    public NodeDataField getExtDescNodeDataToReturn() {
        return _ext_desc_data_to_return;
    }

    public boolean isAllowErrorsInDistanceToParent() {
        return _allow_errors_in_distance_to_parent;
    }

    final public boolean isLineUpRendarableNodeData() {
        return _line_up_renderable_node_data;
    }

    final public boolean isRightLineUpDomains() {
        return _right_align_domains;
    }

    public final boolean isShowDomainLabels() {
        return _show_domain_labels;
    }

    public final void setAllowErrorsInDistanceToParent(final boolean allow_errors_in_distance_to_parent) {
        _allow_errors_in_distance_to_parent = allow_errors_in_distance_to_parent;
    }

    public void setColorLabelsSameAsParentBranch(final boolean color_labels_same_as_parent_branch) {
        _color_labels_same_as_parent_branch = color_labels_same_as_parent_branch;
    }

    public void setExtDescNodeDataToReturn(final NodeDataField ext_desc_data_to_return) {
        _ext_desc_data_to_return = ext_desc_data_to_return;
    }

    final public void setLineUpRendarableNodeData(final boolean line_up_renderable_node_data) {
        _line_up_renderable_node_data = line_up_renderable_node_data;
    }

    final public void setRightLineUpDomains(final boolean right_align_domains) {
        _right_align_domains = right_align_domains;
    }

    public void setShowDomainLabels(final boolean show_domain_labels) {
        _show_domain_labels = show_domain_labels;
    }

    final private void init() {
        _default_node_shape = NodeShape.CIRCLE;
        // GRADIENT was retired as a user-selectable node fill, so do not default to it.
        _default_node_fill = NodeFill.SOLID;
        _default_node_shape_size = AptxConstants.DEFAULT_NODE_SHAPE_SIZE_DEFAULT;
        _internal_number_are_confidence_for_nh_parsing = false;
        _show_scale = false;
        _antialias_print = true;
        _graphics_export_visible_only = false;
        _editable = true;
        _show_default_node_shapes_internal = false;
        // Publication-style default: internal-node labels (e.g. clade names from "Annotate Clades by
        // Rank") sit to the LEFT of the node, right-aligned, on top of the incoming branch -- where the
        // space is usually empty -- rather than to the right, where the node's own subtree fans out.
        _internal_labels_above_branch = true;
        // Publication convention: scientific (Latin) taxonomic names are set in italics (e.g. Homo
        // sapiens). On by default for figure-ready output; only the scientific-name run is italicized,
        // not the taxonomy code, common name or rank.
        _use_italic_scientific_names = true;
        // SVG/EPS export (VectorGraphics2D) embeds no fonts, so by default all text is rendered as glyph
        // outlines to stop viewers substituting the bundled font (EPS -> Times serif, SVG -> generic
        // sans). Turn off to keep selectable/searchable text at the risk of viewer font substitution.
        _outline_fonts_in_vector_export = true;
        _show_default_node_shapes_external = false;
        _show_default_node_shapes_for_marked_nodes = false;
        _color_all_found_nodes_when_coloring_subtree = false;
        _parse_beast_style_extended_nexus_tags = true;
        _min_confidence_value = MIN_CONFIDENCE_DEFAULT;
        _support_visualization = SUPPORT_VISUALIZATION.NONE;
        _support_threshold = SUPPORT_THRESHOLD_DEFAULT;
        _print_black_and_white = false;
        _phylogeny_graphics_type = PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR;
        _base_font = new Font(Configuration.getDefaultFontFamilyName(), Font.PLAIN, AptxConstants.DEFAULT_TREE_FONT_SIZE);
        _match_whole_terms_only = false;
        _search_with_regex = false;
        _search_case_sensitive = false;
        _print_line_width = AptxConstants.PDF_LINE_WIDTH_DEFAULT;
        _show_overview = true;
        _ov_placement = OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT;
        _node_label_direction = NODE_LABEL_DIRECTION.HORIZONTAL;
        _inverse_search_result = false;
        _scale_bar_length = 0.0;
        _number_of_digits_after_comma_for_branch_length_values = AptxConstants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_BRANCH_LENGTH_VALUES_DEFAULT;
        _number_of_digits_after_comma_for_confidence_values = AptxConstants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_CONFIDENCE_VALUES_DEFAULT;
        _nh_parsing_replace_underscores = false;
        _taxonomy_extraction = TAXONOMY_EXTRACTION.NO;
        _cladogram_type = AptxConstants.CLADOGRAM_TYPE_DEFAULT;
        _show_domain_labels = true;
        setAbbreviateScientificTaxonNames(false);
        _color_labels_same_as_parent_branch = false;
        _show_confidence_stddev = false;
        _show_mad_confidence = false;
        _nh_conversion_support_value_style = NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE;
        _ext_desc_data_to_return = NodeDataField.UNKNOWN;
        _line_up_renderable_node_data = true;
        _right_align_domains = false;
        _collapsed_with_average_height = true;
        _show_abbreviated_labels_for_collapsed_nodes = true;
        _search_properties = true;
        _default_branch_width = 1;
    }

    final private void setNumberOfDigitsAfterCommaForBranchLength(final short number_of_digits_after_comma_for_branch_length_values) {
        _number_of_digits_after_comma_for_branch_length_values = number_of_digits_after_comma_for_branch_length_values;
    }

    final private void setNumberOfDigitsAfterCommaForConfidenceValues(final short number_of_digits_after_comma_for_confidence_values) {
        _number_of_digits_after_comma_for_confidence_values = number_of_digits_after_comma_for_confidence_values;
    }

    final Font getBaseFont() {
        return _base_font;
    }

    final CLADOGRAM_TYPE getCladogramType() {
        return _cladogram_type;
    }

    final NodeFill getDefaultNodeFill() {
        return _default_node_fill;
    }

    final NodeShape getDefaultNodeShape() {
        return _default_node_shape;
    }

    final short getDefaultNodeShapeSize() {
        return _default_node_shape_size;
    }

    final double getMinConfidenceValue() {
        return _min_confidence_value;
    }

    final SUPPORT_VISUALIZATION getSupportVisualization() {
        return _support_visualization;
    }

    final double getSupportThreshold() {
        return _support_threshold;
    }

    NH_CONVERSION_SUPPORT_VALUE_STYLE getNhConversionSupportValueStyle() {
        return _nh_conversion_support_value_style;
    }

    final NODE_LABEL_DIRECTION getNodeLabelDirection() {
        return _node_label_direction;
    }

    final short getNumberOfDigitsAfterCommaForBranchLengthValues() {
        return _number_of_digits_after_comma_for_branch_length_values;
    }

    final short getNumberOfDigitsAfterCommaForConfidenceValues() {
        return _number_of_digits_after_comma_for_confidence_values;
    }

    final OVERVIEW_PLACEMENT_TYPE getOvPlacement() {
        return _ov_placement;
    }

    final PHYLOGENY_GRAPHICS_TYPE getPhylogenyGraphicsType() {
        return _phylogeny_graphics_type;
    }

    final float getPrintLineWidth() {
        return _print_line_width;
    }

    final double getScaleBarLength() {
        return _scale_bar_length;
    }

    final TAXONOMY_EXTRACTION getTaxonomyExtraction() {
        return _taxonomy_extraction;
    }

    final boolean isAbbreviateScientificTaxonNames() {
        return _abbreviate_scientific_names;
    }

    boolean isAllowMagnificationOfTaxonomyImages() {
        return true;
    }

    final boolean isAntialiasPrint() {
        return _antialias_print;
    }

    final boolean isColorByTaxonomicGroup() {
        return _color_by_taxonomic_group;
    }

    final boolean isColorLabelsSameAsParentBranch() {
        return _color_labels_same_as_parent_branch;
    }

    final boolean isEditable() {
        return _editable;
    }

    final boolean isGraphicsExportUsingActualSize() {
        return _graphics_export_using_actual_size;
    }

    final boolean isGraphicsExportVisibleOnly() {
        return _graphics_export_visible_only;
    }

    final boolean isInternalNumberAreConfidenceForNhParsing() {
        return _internal_number_are_confidence_for_nh_parsing;
    }

    final boolean isInverseSearchResult() {
        return _inverse_search_result;
    }

    public boolean isSearchProperties() {
        return _search_properties;
    }


    public void setSearchProperties(final boolean search_properties) {
        _search_properties = search_properties;
    }

    final boolean isMatchWholeTermsOnly() {
        return _match_whole_terms_only;
    }

    final boolean isPrintBlackAndWhite() {
        return _print_black_and_white;
    }

    final boolean isPrintUsingActualSize() {
        return _print_using_actual_size;
    }

    final boolean isReplaceUnderscoresInNhParsing() {
        return _nh_parsing_replace_underscores;
    }

    final boolean isSearchCaseSensitive() {
        return _search_case_sensitive;
    }

    final boolean isSearchWithRegex() {
        return _search_with_regex;
    }

    boolean isShowConfidenceStddev() {
        return _show_confidence_stddev;
    }

    boolean isShowMadConfidence() {
        return _show_mad_confidence;
    }

    void setShowMadConfidence(final boolean show_mad_confidence) {
        _show_mad_confidence = show_mad_confidence;
    }

    boolean isShowDefaultNodeShapesExternal() {
        return _show_default_node_shapes_external;
    }

    boolean isShowDefaultNodeShapesForMarkedNodes() {
        return _show_default_node_shapes_for_marked_nodes;
    }

    boolean isShowDefaultNodeShapesInternal() {
        return _show_default_node_shapes_internal;
    }

    boolean isInternalLabelsAboveBranch() {
        return _internal_labels_above_branch;
    }

    boolean isUseItalicScientificNames() {
        return _use_italic_scientific_names;
    }

    boolean isOutlineFontsInVectorExport() {
        return _outline_fonts_in_vector_export;
    }

    final boolean isShowOverview() {
        return _show_overview;
    }

    final boolean isShowScale() {
        return _show_scale;
    }

    final void setAbbreviateScientificTaxonNames(final boolean abbreviate_scientific_names) {
        _abbreviate_scientific_names = abbreviate_scientific_names;
    }

    final void setAntialiasPrint(final boolean antialias_print) {
        _antialias_print = antialias_print;
    }

    final void setBaseFont(final Font base_font) {
        _base_font = base_font;
    }

    final void setCladogramType(final CLADOGRAM_TYPE cladogram_type) {
        _cladogram_type = cladogram_type;
    }

    final void setColorByTaxonomicGroup(final boolean color_by_taxonomic_group) {
        _color_by_taxonomic_group = color_by_taxonomic_group;
    }

    final void setDefaultNodeFill(final NodeFill default_node_fill) {
        _default_node_fill = default_node_fill;
    }

    final void setDefaultNodeShape(final NodeShape default_node_shape) {
        _default_node_shape = default_node_shape;
    }

    final void setDefaultNodeShapeSize(final short default_node_shape_size) {
        _default_node_shape_size = default_node_shape_size;
    }

    final void setEditable(final boolean editable) {
        _editable = editable;
    }

    final void setGraphicsExportVisibleOnly(final boolean graphics_export_visible_only) {
        _graphics_export_visible_only = graphics_export_visible_only;
    }

    final void setInternalNumberAreConfidenceForNhParsing(final boolean internal_number_are_confidence_for_nh_parsing) {
        _internal_number_are_confidence_for_nh_parsing = internal_number_are_confidence_for_nh_parsing;
    }

    final void setInverseSearchResult(final boolean inverse_search_result) {
        _inverse_search_result = inverse_search_result;
    }

    final void setMatchWholeTermsOnly(final boolean search_whole_words_only) {
        _match_whole_terms_only = search_whole_words_only;
    }

    final void setMinConfidenceValue(final double min_confidence_value) {
        _min_confidence_value = min_confidence_value;
    }

    final void setSupportVisualization(final SUPPORT_VISUALIZATION support_visualization) {
        _support_visualization = support_visualization;
    }

    final void setSupportThreshold(final double support_threshold) {
        _support_threshold = support_threshold;
    }

    void setNhConversionSupportValueStyle(final NH_CONVERSION_SUPPORT_VALUE_STYLE nh_conversion_support_value_style) {
        _nh_conversion_support_value_style = nh_conversion_support_value_style;
    }

    final void setNodeLabelDirection(final NODE_LABEL_DIRECTION node_label_direction) {
        _node_label_direction = node_label_direction;
    }

    final void setOvPlacement(final OVERVIEW_PLACEMENT_TYPE ov_placement) {
        _ov_placement = ov_placement;
    }

    final void setPhylogenyGraphicsType(final PHYLOGENY_GRAPHICS_TYPE phylogeny_graphics_type) {
        _phylogeny_graphics_type = phylogeny_graphics_type;
    }

    final void setPrintBlackAndWhite(final boolean print_black_and_white) {
        _print_black_and_white = print_black_and_white;
    }

    final void setPrintLineWidth(final float print_line_width) {
        _print_line_width = print_line_width;
    }

    final void setReplaceUnderscoresInNhParsing(final boolean nh_parsing_replace_underscores) {
        _nh_parsing_replace_underscores = nh_parsing_replace_underscores;
    }

    final void setScaleBarLength(final double scale_bar_length) {
        _scale_bar_length = scale_bar_length;
    }

    final void setSearchCaseSensitive(final boolean search_case_sensitive) {
        _search_case_sensitive = search_case_sensitive;
    }

    final void setSearchWithRegex(final boolean search_with_regex) {
        _search_with_regex = search_with_regex;
    }

    void setShowConfidenceStddev(final boolean show_confidence_stddev) {
        _show_confidence_stddev = show_confidence_stddev;
    }

    void setShowDefaultNodeShapesExternal(final boolean show_default_node_shapes_external) {
        _show_default_node_shapes_external = show_default_node_shapes_external;
    }

    void setShowDefaultNodeShapesForMarkedNodes(final boolean show_default_node_shapes_for_marked_nodes) {
        _show_default_node_shapes_for_marked_nodes = show_default_node_shapes_for_marked_nodes;
    }

    void setShowDefaultNodeShapesInternal(final boolean show_default_node_shapes_internal) {
        _show_default_node_shapes_internal = show_default_node_shapes_internal;
    }

    void setInternalLabelsAboveBranch(final boolean internal_labels_above_branch) {
        _internal_labels_above_branch = internal_labels_above_branch;
    }

    void setUseItalicScientificNames(final boolean use_italic_scientific_names) {
        _use_italic_scientific_names = use_italic_scientific_names;
    }

    void setOutlineFontsInVectorExport(final boolean outline_fonts_in_vector_export) {
        _outline_fonts_in_vector_export = outline_fonts_in_vector_export;
    }

    final void setShowOverview(final boolean show_overview) {
        _show_overview = show_overview;
    }

    final void setShowScale(final boolean show_scale) {
        _show_scale = show_scale;
    }

    final void setTaxonomyExtraction(final TAXONOMY_EXTRACTION taxonomy_extraction) {
        _taxonomy_extraction = taxonomy_extraction;
    }

    public final static Options createInstance(final Configuration configuration) {
        final Options instance = createDefaultInstance();
        if (configuration != null) {
            instance.setShowScale(configuration.isShowScale());
            instance.setShowOverview(configuration.isShowOverview());
            instance.setColorByTaxonomicGroup(configuration.isColorByTaxonomicGroup());
            instance.setCladogramType(configuration.getCladogramType());
            instance.setOvPlacement(configuration.getOvPlacement());
            instance.setPrintLineWidth(configuration.getPrintLineWidth());
            instance.setNodeLabelDirection(configuration.getNodeLabelDirection());
            if (configuration.getNumberOfDigitsAfterCommaForBranchLengthValues() >= 0) {
                instance.setNumberOfDigitsAfterCommaForBranchLength(configuration
                        .getNumberOfDigitsAfterCommaForBranchLengthValues());
            }
            if (configuration.getNumberOfDigitsAfterCommaForConfidenceValues() >= 0) {
                instance.setNumberOfDigitsAfterCommaForConfidenceValues(configuration
                        .getNumberOfDigitsAfterCommaForConfidenceValues());
            }
            instance.setTaxonomyExtraction(configuration.getTaxonomyExtraction());
            instance.setReplaceUnderscoresInNhParsing(configuration.isReplaceUnderscoresInNhParsing());
            instance.setInternalNumberAreConfidenceForNhParsing(configuration
                    .isInternalNumberAreConfidenceForNhParsing());
            instance.setEditable(configuration.isEditable());
            instance.setColorLabelsSameAsParentBranch(configuration.isColorLabelsSameAsParentBranch());
            instance.setShowDomainLabels(configuration.isShowDomainLabels());
            instance.setAbbreviateScientificTaxonNames(configuration.isAbbreviateScientificTaxonNames());
            if (configuration.getMinConfidenceValue() != MIN_CONFIDENCE_DEFAULT) {
                instance.setMinConfidenceValue(configuration.getMinConfidenceValue());
            }
            if (configuration.getBaseFontSize() > 0) {
                instance.setBaseFont(instance.getBaseFont().deriveFont((float) configuration.getBaseFontSize()));
            }
            if (!ForesterUtil.isEmpty(configuration.getBaseFontFamilyName())) {
                instance.setBaseFont(new Font(configuration.getBaseFontFamilyName(), Font.PLAIN, instance
                        .getBaseFont().getSize()));
            }
            if (configuration.getPhylogenyGraphicsType() != null) {
                instance.setPhylogenyGraphicsType(configuration.getPhylogenyGraphicsType());
            }
            if (configuration.getDefaultNodeFill() != null) {
                instance.setDefaultNodeFill(configuration.getDefaultNodeFill());
            }
            if (configuration.getDefaultNodeShape() != null) {
                instance.setDefaultNodeShape(configuration.getDefaultNodeShape());
            }
            if (configuration.getDefaultNodeShapeSize() >= 0) {
                instance.setDefaultNodeShapeSize(configuration.getDefaultNodeShapeSize());
            }
            instance.setShowDefaultNodeShapesInternal(configuration.isShowDefaultNodeShapesInternal());
            instance.setShowDefaultNodeShapesExternal(configuration.isShowDefaultNodeShapesExternal());
            instance.setShowDefaultNodeShapesForMarkedNodes(configuration.isShowDefaultNodeShapesForMarkedNodes());
            if (configuration.getExtDescNodeDataToReturn() != null) {
                instance.setExtDescNodeDataToReturn(configuration.getExtDescNodeDataToReturn());
            }
            instance.setRightLineUpDomains(configuration.isRightLineUpDomains());
            instance.setLineUpRendarableNodeData(configuration.isLineUpRendarableNodeData());
            instance.setAllowErrorsInDistanceToParent(false);
        }
        return instance;
    }

    final static Options createDefaultInstance() {
        return new Options();
    }

    final boolean isColorAllFoundNodesWhenColoringSubtree() {
        return _color_all_found_nodes_when_coloring_subtree;
    }

    final void setColorAllFoundNodesWhenColoringSubtree(final boolean color_all_found_nodes_when_coloring_subtree) {
        _color_all_found_nodes_when_coloring_subtree = color_all_found_nodes_when_coloring_subtree;
    }

    final boolean isParseBeastStyleExtendedNexusTags() {
        return _parse_beast_style_extended_nexus_tags;
    }

    final void setParseBeastStyleExtendedNexusTags(final boolean parse_beast_style_extended_nexus_tags) {
        _parse_beast_style_extended_nexus_tags = parse_beast_style_extended_nexus_tags;
    }

    final boolean isCollapsedWithAverageHeigh() {
        return _collapsed_with_average_height;
    }

    final void setCollapsedWithAverageHeigh(final boolean collapsed_with_average_height) {
        _collapsed_with_average_height = collapsed_with_average_height;
    }

    final boolean isShowAbbreviatedLabelsForCollapsedNodes() {
        return _show_abbreviated_labels_for_collapsed_nodes;
    }

    final void setShowAbbreviatedLabelsForCollapsedNodes(final boolean show_abbreviated_labels_for_collapsed_nodes) {
        _show_abbreviated_labels_for_collapsed_nodes = show_abbreviated_labels_for_collapsed_nodes;
    }


}
