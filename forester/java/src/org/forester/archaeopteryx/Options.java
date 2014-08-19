// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2009 Christian M. Zmasek
// Copyright (C) 2009 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx;

import java.awt.Font;

import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.NodeData;
import org.forester.phylogeny.data.NodeData.NODE_DATA;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.util.ForesterUtil;
import org.omg.stub.java.rmi._Remote_Stub;

/*
 * This is to hold changeable options.
 */
final public class Options {

    static final double                       MIN_CONFIDENCE_DEFAULT = 0.0;
    private boolean                           _abbreviate_scientific_names;
    private boolean                           _allow_errors_in_distance_to_parent;
    private boolean                           _antialias_print;
    private boolean                           _antialias_screen;
    private boolean                           _background_color_gradient;
    private Font                              _base_font;
    private CLADOGRAM_TYPE                    _cladogram_type;
    private boolean                           _color_by_taxonomic_group;
    private boolean                           _color_labels_same_as_parent_branch;
    private NodeVisualData.NodeFill           _default_node_fill;
    private NodeVisualData.NodeShape          _default_node_shape;
    private short                             _default_node_shape_size;
    private boolean                           _editable;
    private NODE_DATA                         _ext_desc_data_to_return;
    private boolean                           _graphics_export_using_actual_size;
    private boolean                           _graphics_export_visible_only;
    private boolean                           _internal_number_are_confidence_for_nh_parsing;
    private boolean                           _inverse_search_result;
    private boolean                           _match_whole_terms_only;
    private double                            _min_confidence_value;
    private NH_CONVERSION_SUPPORT_VALUE_STYLE _nh_conversion_support_value_style;
    private boolean                           _nh_parsing_replace_underscores;
    private NODE_LABEL_DIRECTION              _node_label_direction;
    private short                             _number_of_digits_after_comma_for_branch_length_values;
    private short                             _number_of_digits_after_comma_for_confidence_values;
    private OVERVIEW_PLACEMENT_TYPE           _ov_placement;
    private PHYLOGENY_GRAPHICS_TYPE           _phylogeny_graphics_type;
    private boolean                           _print_black_and_white;
    private float                             _print_line_width;
    private int                               _print_size_x;
    private int                               _print_size_y;
    private boolean                           _print_using_actual_size;
    private double                            _scale_bar_length;
    private boolean                           _search_case_sensitive;
    private boolean                           _show_annotation_ref_source;
    private boolean                           _show_branch_length_values;
    private boolean                           _show_confidence_stddev;
    private boolean                           _show_default_node_shapes_external;
    private boolean                           _show_default_node_shapes_internal;
    private boolean                           _show_domain_labels;
    private boolean                           _show_overview;
    private boolean                           _show_scale;
    private TAXONOMY_EXTRACTION               _taxonomy_extraction;
    private boolean _line_up_renderable_node_data;
    private boolean _right_align_domains;

    private Options() {
        init();
    }

    public NodeData.NODE_DATA getExtDescNodeDataToReturn() {
        return _ext_desc_data_to_return;
    }

    public boolean isAllowErrorsInDistanceToParent() {
        return _allow_errors_in_distance_to_parent;
    }

    public boolean isAllowFontSizeChange() {
        return true;
    }

    public final boolean isShowAnnotationRefSource() {
        return _show_annotation_ref_source;
    }

    public final boolean isShowDomainLabels() {
        return _show_domain_labels;
    }

    public final void setAllowErrorsInDistanceToParent( final boolean allow_errors_in_distance_to_parent ) {
        _allow_errors_in_distance_to_parent = allow_errors_in_distance_to_parent;
    }

    public void setBackgroundColorGradient( final boolean background_color_gradient ) {
        _background_color_gradient = background_color_gradient;
    }

    public void setColorLabelsSameAsParentBranch( final boolean color_labels_same_as_parent_branch ) {
        _color_labels_same_as_parent_branch = color_labels_same_as_parent_branch;
    }

    public void setExtDescNodeDataToReturn( final NODE_DATA ext_desc_data_to_return ) {
        _ext_desc_data_to_return = ext_desc_data_to_return;
    }

    public final void setShowAnnotationRefSource( final boolean show_annotation_ref_source ) {
        _show_annotation_ref_source = show_annotation_ref_source;
    }

    public void setShowDomainLabels( final boolean show_domain_labels ) {
        _show_domain_labels = show_domain_labels;
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

    final int getPrintSizeX() {
        return _print_size_x;
    }

    final int getPrintSizeY() {
        return _print_size_y;
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

    final boolean isAntialiasScreen() {
        return _antialias_screen;
    }

    final boolean isBackgroundColorGradient() {
        return _background_color_gradient;
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

    final boolean isShowBranchLengthValues() {
        return _show_branch_length_values;
    }

    boolean isShowConfidenceStddev() {
        return _show_confidence_stddev;
    }

    boolean isShowDefaultNodeShapesExternal() {
        return _show_default_node_shapes_external;
    }

    boolean isShowDefaultNodeShapesInternal() {
        return _show_default_node_shapes_internal;
    }

    final boolean isShowOverview() {
        return _show_overview;
    }

    final boolean isShowScale() {
        return _show_scale;
    }

    final void setAbbreviateScientificTaxonNames( final boolean abbreviate_scientific_names ) {
        _abbreviate_scientific_names = abbreviate_scientific_names;
    }

    final void setAntialiasPrint( final boolean antialias_print ) {
        _antialias_print = antialias_print;
    }

    final void setAntialiasScreen( final boolean antialias_screen ) {
        _antialias_screen = antialias_screen;
    }

    final void setBaseFont( final Font base_font ) {
        _base_font = base_font;
    }

    final void setCladogramType( final CLADOGRAM_TYPE cladogram_type ) {
        _cladogram_type = cladogram_type;
    }

    final void setColorByTaxonomicGroup( final boolean color_by_taxonomic_group ) {
        _color_by_taxonomic_group = color_by_taxonomic_group;
    }

    final void setDefaultNodeFill( final NodeFill default_node_fill ) {
        _default_node_fill = default_node_fill;
    }

    final void setDefaultNodeShape( final NodeShape default_node_shape ) {
        _default_node_shape = default_node_shape;
    }

    final void setDefaultNodeShapeSize( final short default_node_shape_size ) {
        _default_node_shape_size = default_node_shape_size;
    }

    final void setEditable( final boolean editable ) {
        _editable = editable;
    }

    final void setGraphicsExportUsingActualSize( final boolean graphics_export_using_actual_size ) {
        _graphics_export_using_actual_size = graphics_export_using_actual_size;
        if ( !graphics_export_using_actual_size ) {
            setGraphicsExportVisibleOnly( false );
        }
    }

    final void setGraphicsExportVisibleOnly( final boolean graphics_export_visible_only ) {
        _graphics_export_visible_only = graphics_export_visible_only;
        if ( graphics_export_visible_only ) {
            setGraphicsExportUsingActualSize( true );
        }
    }

    final void setInternalNumberAreConfidenceForNhParsing( final boolean internal_number_are_confidence_for_nh_parsing ) {
        _internal_number_are_confidence_for_nh_parsing = internal_number_are_confidence_for_nh_parsing;
    }

    final void setInverseSearchResult( final boolean inverse_search_result ) {
        _inverse_search_result = inverse_search_result;
    }

    final void setMatchWholeTermsOnly( final boolean search_whole_words_only ) {
        _match_whole_terms_only = search_whole_words_only;
    }

    final void setMinConfidenceValue( final double min_confidence_value ) {
        _min_confidence_value = min_confidence_value;
    }

    void setNhConversionSupportValueStyle( final NH_CONVERSION_SUPPORT_VALUE_STYLE nh_conversion_support_value_style ) {
        _nh_conversion_support_value_style = nh_conversion_support_value_style;
    }

    final void setNodeLabelDirection( final NODE_LABEL_DIRECTION node_label_direction ) {
        _node_label_direction = node_label_direction;
    }

    final void setOvPlacement( final OVERVIEW_PLACEMENT_TYPE ov_placement ) {
        _ov_placement = ov_placement;
    }

    final void setPhylogenyGraphicsType( final PHYLOGENY_GRAPHICS_TYPE phylogeny_graphics_type ) {
        _phylogeny_graphics_type = phylogeny_graphics_type;
    }

    final void setPrintBlackAndWhite( final boolean print_black_and_white ) {
        _print_black_and_white = print_black_and_white;
    }

    final void setPrintLineWidth( final float print_line_width ) {
        _print_line_width = print_line_width;
    }

    final void setPrintSizeX( final int print_size_x ) {
        _print_size_x = print_size_x;
    }

    final void setPrintSizeY( final int print_size_y ) {
        _print_size_y = print_size_y;
    }

    final void setPrintUsingActualSize( final boolean print_using_actual_size ) {
        _print_using_actual_size = print_using_actual_size;
    }

    final void setReplaceUnderscoresInNhParsing( final boolean nh_parsing_replace_underscores ) {
        _nh_parsing_replace_underscores = nh_parsing_replace_underscores;
    }

    final void setScaleBarLength( final double scale_bar_length ) {
        _scale_bar_length = scale_bar_length;
    }

    final void setSearchCaseSensitive( final boolean search_case_sensitive ) {
        _search_case_sensitive = search_case_sensitive;
    }

    final void setShowBranchLengthValues( final boolean show_branch_length_values ) {
        _show_branch_length_values = show_branch_length_values;
    }

    void setShowConfidenceStddev( final boolean show_confidence_stddev ) {
        _show_confidence_stddev = show_confidence_stddev;
    }

    void setShowDefaultNodeShapesExternal( final boolean show_default_node_shapes_external ) {
        _show_default_node_shapes_external = show_default_node_shapes_external;
    }

    void setShowDefaultNodeShapesInternal( final boolean show_default_node_shapes_internal ) {
        _show_default_node_shapes_internal = show_default_node_shapes_internal;
    }

    final void setShowOverview( final boolean show_overview ) {
        _show_overview = show_overview;
    }

    final void setShowScale( final boolean show_scale ) {
        _show_scale = show_scale;
    }

    final void setTaxonomyExtraction( final TAXONOMY_EXTRACTION taxonomy_extraction ) {
        _taxonomy_extraction = taxonomy_extraction;
    }

    final private void init() {
        _default_node_shape = NodeShape.CIRCLE;
        _default_node_fill = NodeFill.GRADIENT;
        _default_node_shape_size = Constants.DEFAULT_NODE_SHAPE_SIZE_DEFAULT;
        _show_branch_length_values = false;
        _internal_number_are_confidence_for_nh_parsing = false;
        _show_scale = false;
        _antialias_screen = true;
        _antialias_print = true;
        _graphics_export_visible_only = false;
        _editable = true;
        _background_color_gradient = false;
        _show_default_node_shapes_internal = false;
        _show_default_node_shapes_external = false;
        if ( AptxUtil.isUsOrCanada() ) {
            _print_size_x = Constants.US_LETTER_SIZE_X;
            _print_size_y = Constants.US_LETTER_SIZE_Y;
        }
        else {
            _print_size_x = Constants.A4_SIZE_X;
            _print_size_y = Constants.A4_SIZE_Y;
        }
        _min_confidence_value = MIN_CONFIDENCE_DEFAULT;
        _print_black_and_white = false;
        _print_using_actual_size = true;
        _graphics_export_using_actual_size = true;
        _phylogeny_graphics_type = PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR;
        _base_font = new Font( Configuration.getDefaultFontFamilyName(), Font.PLAIN, 10 );
        _match_whole_terms_only = false;
        _search_case_sensitive = false;
        _print_line_width = Constants.PDF_LINE_WIDTH_DEFAULT;
        _show_overview = true;
        _ov_placement = OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT;
        _node_label_direction = NODE_LABEL_DIRECTION.HORIZONTAL;
        _inverse_search_result = false;
        _scale_bar_length = 0.0;
        _number_of_digits_after_comma_for_branch_length_values = Constants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_BRANCH_LENGTH_VALUES_DEFAULT;
        _number_of_digits_after_comma_for_confidence_values = Constants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_CONFIDENCE_VALUES_DEFAULT;
        _nh_parsing_replace_underscores = false;
        _taxonomy_extraction = TAXONOMY_EXTRACTION.NO;
        _cladogram_type = Constants.CLADOGRAM_TYPE_DEFAULT;
        _show_domain_labels = true;
        _show_annotation_ref_source = true;
        setAbbreviateScientificTaxonNames( false );
        _color_labels_same_as_parent_branch = false;
        _show_confidence_stddev = true;
        _nh_conversion_support_value_style = NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE;
        _ext_desc_data_to_return = NODE_DATA.UNKNOWN;
        _line_up_renderable_node_data = false;
        _right_align_domains = false;
    }

    final private void setNumberOfDigitsAfterCommaForBranchLength( final short number_of_digits_after_comma_for_branch_length_values ) {
        _number_of_digits_after_comma_for_branch_length_values = number_of_digits_after_comma_for_branch_length_values;
    }

    final private void setNumberOfDigitsAfterCommaForConfidenceValues( final short number_of_digits_after_comma_for_confidence_values ) {
        _number_of_digits_after_comma_for_confidence_values = number_of_digits_after_comma_for_confidence_values;
    }

    public final static Options createInstance( final Configuration configuration ) {
        final Options instance = createDefaultInstance();
        if ( configuration != null ) {
            instance.setAntialiasScreen( configuration.isAntialiasScreen() );
            instance.setShowScale( configuration.isShowScale() );
            instance.setShowBranchLengthValues( configuration.isShowBranchLengthValues() );
            instance.setShowOverview( configuration.isShowOverview() );
            instance.setColorByTaxonomicGroup( configuration.isColorByTaxonomicGroup() );
            instance.setCladogramType( configuration.getCladogramType() );
            instance.setOvPlacement( configuration.getOvPlacement() );
            instance.setPrintLineWidth( configuration.getPrintLineWidth() );
            instance.setNodeLabelDirection( configuration.getNodeLabelDirection() );
            instance.setBackgroundColorGradient( configuration.isBackgroundColorGradient() );
            if ( configuration.getNumberOfDigitsAfterCommaForBranchLengthValues() >= 0 ) {
                instance.setNumberOfDigitsAfterCommaForBranchLength( configuration
                        .getNumberOfDigitsAfterCommaForBranchLengthValues() );
            }
            if ( configuration.getNumberOfDigitsAfterCommaForConfidenceValues() >= 0 ) {
                instance.setNumberOfDigitsAfterCommaForConfidenceValues( configuration
                        .getNumberOfDigitsAfterCommaForConfidenceValues() );
            }
            instance.setTaxonomyExtraction( configuration.getTaxonomyExtraction() );
            instance.setReplaceUnderscoresInNhParsing( configuration.isReplaceUnderscoresInNhParsing() );
            instance.setInternalNumberAreConfidenceForNhParsing( configuration
                    .isInternalNumberAreConfidenceForNhParsing() );
            instance.setEditable( configuration.isEditable() );
            instance.setColorLabelsSameAsParentBranch( configuration.isColorLabelsSameAsParentBranch() );
            instance.setShowDomainLabels( configuration.isShowDomainLabels() );
            instance.setShowAnnotationRefSource( configuration.isShowAnnotationRefSource() );
            instance.setAbbreviateScientificTaxonNames( configuration.isAbbreviateScientificTaxonNames() );
            if ( configuration.getMinConfidenceValue() != MIN_CONFIDENCE_DEFAULT ) {
                instance.setMinConfidenceValue( configuration.getMinConfidenceValue() );
            }
            if ( configuration.getGraphicsExportX() > 0 ) {
                instance.setPrintSizeX( configuration.getGraphicsExportX() );
            }
            if ( configuration.getGraphicsExportY() > 0 ) {
                instance.setPrintSizeY( configuration.getGraphicsExportY() );
            }
            if ( configuration.getBaseFontSize() > 0 ) {
                instance.setBaseFont( instance.getBaseFont().deriveFont( ( float ) configuration.getBaseFontSize() ) );
            }
            if ( !ForesterUtil.isEmpty( configuration.getBaseFontFamilyName() ) ) {
                instance.setBaseFont( new Font( configuration.getBaseFontFamilyName(), Font.PLAIN, instance
                        .getBaseFont().getSize() ) );
            }
            if ( configuration.getPhylogenyGraphicsType() != null ) {
                instance.setPhylogenyGraphicsType( configuration.getPhylogenyGraphicsType() );
            }
            if ( configuration.getDefaultNodeFill() != null ) {
                instance.setDefaultNodeFill( configuration.getDefaultNodeFill() );
            }
            if ( configuration.getDefaultNodeShape() != null ) {
                instance.setDefaultNodeShape( configuration.getDefaultNodeShape() );
            }
            if ( configuration.getDefaultNodeShapeSize() >= 0 ) {
                instance.setDefaultNodeShapeSize( configuration.getDefaultNodeShapeSize() );
            }
            instance.setShowDefaultNodeShapesInternal( configuration.isShowDefaultNodeShapesInternal() );
            instance.setShowDefaultNodeShapesExternal( configuration.isShowDefaultNodeShapesExternal() );
            if ( configuration.getExtDescNodeDataToReturn() != null ) {
                instance.setExtDescNodeDataToReturn( configuration.getExtDescNodeDataToReturn() );
            }
            instance.setAllowErrorsInDistanceToParent( false );
        }
        return instance;
    }

    final static Options createDefaultInstance() {
        return new Options();
    }

    public static enum CLADOGRAM_TYPE {
        EXT_NODE_SUM_DEP, NON_LINED_UP, TOTAL_NODE_SUM_DEP;
    }

    public static enum NODE_LABEL_DIRECTION {
        HORIZONTAL, RADIAL;
    }

    public static enum PHYLOGENY_GRAPHICS_TYPE {
        CIRCULAR, CONVEX, CURVED, EURO_STYLE, RECTANGULAR, ROUNDED, TRIANGULAR, UNROOTED;
    }

    static enum OVERVIEW_PLACEMENT_TYPE {
        LOWER_LEFT( "lower left" ),
        LOWER_RIGHT( "lower right" ),
        UPPER_LEFT( "upper left" ),
        UPPER_RIGHT( "upper right" );

        private final String _name;

        private OVERVIEW_PLACEMENT_TYPE( final String name ) {
            _name = name;
        }

        @Override
        public String toString() {
            return _name;
        }

        public String toTag() {
            return toString().replaceAll( " ", "_" );
        }
    }

    final public boolean isLineUpRendarableNodeData() {
       
        return _line_up_renderable_node_data;
    }
    
    final public boolean isRightLineUpDomains() {
        
        return _right_align_domains;
    }
    
    final public void setLineUpRendarableNodeData( final boolean line_up_renderable_node_data) {
        
         _line_up_renderable_node_data = line_up_renderable_node_data;
    }
    
    final public void setRightLineUpDomains( final boolean right_align_domains ) {
        
        _right_align_domains = right_align_domains;
    }
    
}
