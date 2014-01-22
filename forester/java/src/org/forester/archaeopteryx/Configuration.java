// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// Copyright (C) 2003-2007 Ethalinda K.S. Cannon
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

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Map;
import java.util.SortedMap;
import java.util.StringTokenizer;
import java.util.TreeMap;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.OVERVIEW_PLACEMENT_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.phylogeny.data.NodeData.NODE_DATA;
import org.forester.phylogeny.data.NodeVisualization;
import org.forester.phylogeny.data.NodeVisualization.NodeFill;
import org.forester.phylogeny.data.NodeVisualization.NodeShape;
import org.forester.util.ForesterUtil;

public final class Configuration {

    final static int                        add_new_node                                           = 14;
    final static int                        blast                                                  = 9;
    final static String                     clickto_options[][]                                    = {
            { "Display Node Data", "display" }, { "Collapse/Uncollapse", "display" }, { "Root/Reroot", "display" },
            { "Sub/Super Tree", "display" }, { "Swap Descendants", "display" },
            { "Colorize Subtree/Node(s)", "display" }, { "Open Sequence DB", "display" }, { "Open PDB", "display" },
            { "Open Taxonomy DB", "display" }, { "Blast", "display" }, { "Cut Subtree", "display" },
            { "Copy Subtree", "display" }, { "Paste Subtree", "display" }, { "Delete Subtree/Node", "display" },
            { "Add New Node", "display" }, { "Edit Node Data", "display" }, { "Sort Descendants", "display" },
            { "Return", "display" }, { "Select Node(s)", "display" }                              };
    final static int                        collapse_uncollapse                                    = 1;
    final static int                        color_according_to_annotation                          = 19;
    final static int                        color_according_to_species                             = 6;
    final static int                        color_branches                                         = 7;
    final static int                        color_subtree                                          = 5;
    final static int                        copy_subtree                                           = 11;
    final static int                        cut_subtree                                            = 10;
    final static int                        delete_subtree_or_node                                 = 13;
    final static int                        display_as_phylogram                                   = 0;
    final static int                        display_internal_data                                  = 15;
    // ------------------
    // Click-to options
    // ------------------
    final static int                        display_node_data                                      = 0;
    final static String                     display_options[][]                                    = {
            { "Phylogram", "display", "?" }, { "Node Name", "display", "yes" }, { "Taxonomy Code", "display", "yes" },
            { "Seq Annotations", "nodisplay", "no" }, { "Confidence Values", "display", "?" },
            { "Node Events", "display", "?" }, { "Colorize by Taxonomy", "display", "no" },
            { "Use Branch Colors", "display", "no" }, { "Use Branch Widths", "display", "no" },
            { "Show Custom Nodes", "display", "yes" }, { "Protein Domains", "nodisplay", "no" },
            { "Binary Characters", "nodisplay", "no" }, { "Binary Char Counts", "nodisplay", "no" },
            { "Seq Name", "display", "yes" }, { "Seq Accession", "display", "no" },
            { "Show Internal Data", "display", "yes" }, { "Dyna Hide", "display", "yes" },
            { "Taxonomy Scientific", "display", "yes" }, { "Taxonomy Common", "display", "no" },
            { "Colorize by Annotation", "nodisplay", "no" }, { "Seq Symbol", "display", "yes" },
            { "Rollover", "display", "yes" }, { "Relation Confidence", "nodisplay", "no" },
            { "Vector Data", "nodisplay", "no" }, { "Taxonomy Images", "display", "no" },
            { "Properties", "nodisplay", "no" }, { "Gene Name", "display", "yes" }                };
    final static int                        dynamically_hide_data                                  = 16;
    final static int                        edit_node_data                                         = 15;
    final static int                        get_ext_desc_data                                      = 17;
    final static int                        node_data_popup                                        = 21;
    final static int                        open_pdb_web                                           = 7;
    final static int                        open_seq_web                                           = 6;
    final static int                        open_tax_web                                           = 8;
    final static int                        paste_subtree                                          = 12;
    final static int                        reroot                                                 = 2;
    final static int                        select_nodes                                           = 18;
    final static int                        show_annotation                                        = 3;
    final static int                        show_binary_character_counts                           = 12;
    final static int                        show_binary_characters                                 = 11;
    final static int                        show_custom_node_shapes                                = 9;
    final static int                        show_domain_architectures                              = 10;
    final static int                        show_gene_names                                        = 26;
    final static int                        show_node_names                                        = 1;
    final static int                        show_properties                                        = 25;
    final static int                        show_relation_confidence                               = 22;
    final static int                        show_seq_names                                         = 13;
    final static int                        show_seq_symbols                                       = 20;
    final static int                        show_sequence_acc                                      = 14;
    final static int                        show_tax_code                                          = 2;
    final static int                        show_taxonomy_common_names                             = 18;
    final static int                        show_taxonomy_images                                   = 24;
    final static int                        show_taxonomy_scientific_names                         = 17;
    final static int                        show_vector_data                                       = 23;
    final static int                        sort_descendents                                       = 16;
    final static int                        subtree                                                = 3;
    final static int                        swap                                                   = 4;
    static final String                     VALIDATE_AGAINST_PHYLOXML_XSD_SCHEMA                   = "validate_against_phyloxml_xsd_schema";
    final static int                        width_branches                                         = 8;
    final static int                        write_confidence_values                                = 4;
    final static int                        write_events                                           = 5;
    // ----------------
    // Function colors
    // ----------------
    private static Hashtable<String, Color> _annotation_colors;
    // ----------------
    // Domain colors
    // ----------------
    private static Hashtable<String, Color> _domain_colors;
    // ----------------
    // Species colors
    // ----------------
    private static Hashtable<String, Color> _species_colors;
    private static String                   DEFAULT_FONT_FAMILY                                    = "";
    private static final int                DEPRECATED                                             = -2;
    private static final String             DISPLAY_COLOR_KEY                                      = "display_color";
    // ---------------------------
    // Display options for trees
    // ---------------------------
    // ---------------------------------
    // Pertaining to the config itself
    // ---------------------------------
    // Full path to config (may be URL)
    String                                  config_filename;
    // This option is selected in the dropdown
    int                                     default_clickto                                        = Configuration.display_node_data;
    String                                  default_config_filename                                = Constants.DEFAULT_CONFIGURATION_FILE_NAME;
    // --------------
    // Color set
    // --------------
    TreeColorSet                            tree_color_set;
    // -------
    // Fonts
    // -------
    TreeFontSet                             tree_font_set;
    boolean                                 verbose                                                = Constants.VERBOSE_DEFAULT;
    private boolean                         _abbreviate_scientific_names                           = false;
    private boolean                         _antialias_screen                                      = true;
    private boolean                         _background_color_gradient                             = false;
    private String                          _base_font_family_name                                 = "";
    private int                             _base_font_size                                        = -1;
    private CLADOGRAM_TYPE                  _cladogram_type                                        = Constants.CLADOGRAM_TYPE_DEFAULT;
    private boolean                         _color_labels_same_as_parent_branch                    = false;
    private int                             _default_bootstrap_samples                             = -1;
    private NodeFill                        _default_node_fill                                     = NodeFill.SOLID;
    private NodeShape                       _default_node_shape                                    = NodeShape.RECTANGLE;
    private short                           _default_node_shape_size                               = Constants.DEFAULT_NODE_SHAPE_SIZE_DEFAULT;
    private SortedMap<String, Color>        _display_colors                                        = null;
    private boolean                         _display_sequence_relations                            = false;
    private Color                           _domain_structure_base_color                           = Constants.DOMAIN_STRUCTURE_BASE_COLOR_DEFAULT;
    private Color                           _domain_structure_font_color                           = Constants.DOMAIN_STRUCTURE_FONT_COLOR_DEFAULT;
    private boolean                         _editable                                              = true;
    private NODE_DATA                       _ext_desc_data_to_return                               = NODE_DATA.UNKNOWN;
    private EXT_NODE_DATA_RETURN_ON         _ext_node_data_return_on                               = EXT_NODE_DATA_RETURN_ON.WINODW;
    private int                             _frame_x_size;
    private int                             _frame_y_size;
    private int                             _graphics_export_x                                     = -1;
    private int                             _graphics_export_y                                     = -1;
    private Color                           _gui_background_color                                  = Constants.GUI_BACKGROUND_DEFAULT;
    private Color                           _gui_button_background_color                           = Constants.BUTTON_BACKGROUND_COLOR_DEFAULT;
    private Color                           _gui_button_border_color                               = Constants.BUTTON_BORDER_COLOR_DEFAULT;
    private Color                           _gui_button_text_color                                 = Constants.BUTTON_TEXT_COLOR_DEFAULT;
    private Color                           _gui_checkbox_and_button_active_color                  = Constants.CHECKBOX_AND_BUTTON_ACTIVE_COLOR_DEFAULT;
    private Color                           _gui_checkbox_text_color                               = Constants.CHECKBOX_TEXT_COLOR_DEFAULT;
    private Color                           _gui_menu_background_color                             = Constants.MENU_BACKGROUND_COLOR_DEFAULT;
    private Color                           _gui_menu_text_color                                   = Constants.MENU_TEXT_COLOR_DEFAULT;
    private boolean                         _hide_controls_and_menus                               = false;
    private boolean                         _internal_number_are_confidence_for_nh_parsing         = false;
    private String                          _label_for_get_ext_descendents_data                    = "";
    private int                             _max_base_font_size                                    = 20;
    private boolean                         _midpoint_root                                         = false;
    private int                             _min_base_font_size                                    = 2;
    private double                          _min_confidence_value                                  = Options.MIN_CONFIDENCE_DEFAULT;
    private boolean                         _nh_parsing_replace_underscores                        = false;
    private NODE_LABEL_DIRECTION            _node_label_direction                                  = NODE_LABEL_DIRECTION.HORIZONTAL;
    private short                           _number_of_digits_after_comma_for_branch_length_values = Constants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_BRANCH_LENGTH_VALUES_DEFAULT;
    private short                           _number_of_digits_after_comma_for_confidence_values    = Constants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_CONFIDENCE_VALUES_DEFAULT;
    private short                           _ov_max_height                                         = 80;
    private short                           _ov_max_width                                          = 80;
    private OVERVIEW_PLACEMENT_TYPE         _ov_placement                                          = OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT;
    private File                            _path_to_local_fastme                                  = null;
    private File                            _path_to_local_mafft                                   = null;
    private File                            _path_to_local_raxml                                   = null;
    private PHYLOGENY_GRAPHICS_TYPE         _phylogeny_graphics_type                               = PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR;
    private float                           _print_line_width                                      = Constants.PDF_LINE_WIDTH_DEFAULT;
    private boolean                         _show_annotation_ref_source                            = true;
    private boolean                         _show_branch_length_values                             = false;
    private boolean                         _show_default_node_shapes_external                     = false;
    private boolean                         _show_default_node_shapes_internal                     = false;
    private boolean                         _show_domain_labels                                    = true;
    private boolean                         _show_overview                                         = true;
    private boolean                         _show_scale                                            = false;
    private boolean                         _taxonomy_colorize_node_shapes                         = false;
    private TAXONOMY_EXTRACTION             _taxonomy_extraction                                   = TAXONOMY_EXTRACTION.NO;
    private UI                              _ui                                                    = UI.UNKNOWN;
    private boolean                         _use_tabbed_display                                    = false;
    private boolean                         _validate_against_phyloxml_xsd_schema                  = Constants.VALIDATE_AGAINST_PHYLOXML_XSD_SCJEMA_DEFAULT;
    static {
        for( final String font_name : Constants.DEFAULT_FONT_CHOICES ) {
            if ( Arrays.binarySearch( AptxUtil.getAvailableFontFamiliesSorted(), font_name ) >= 0 ) {
                DEFAULT_FONT_FAMILY = font_name;
                break;
            }
        }
        if ( ForesterUtil.isEmpty( DEFAULT_FONT_FAMILY ) ) {
            DEFAULT_FONT_FAMILY = Constants.DEFAULT_FONT_CHOICES[ Constants.DEFAULT_FONT_CHOICES.length - 1 ];
        }
    }

    public Configuration() {
        this( null, false, false, false );
    }

    public Configuration( final String cf, final boolean is_url, final boolean is_applet, final boolean verbose ) {
        if ( ForesterUtil.isEmpty( cf ) ) {
            config_filename = default_config_filename;
        }
        else {
            config_filename = cf;
        }
        setDisplayColors( new TreeMap<String, Color>() );
        config_filename = config_filename.trim();
        URL u = null;
        if ( is_url ) {
            // If URL, open accordingly
            try {
                u = new URL( config_filename );
                try {
                    final InputStreamReader isr = new InputStreamReader( u.openStream() );
                    final BufferedReader bf = new BufferedReader( isr );
                    readConfig( bf );
                    bf.close();
                    ForesterUtil.programMessage( Constants.PRG_NAME, "successfully read from configuration url ["
                            + config_filename + "]" );
                }
                catch ( final Exception e ) {
                    ForesterUtil.printWarningMessage( Constants.PRG_NAME, "failed to read configuration from ["
                            + config_filename + "]: " + e.getLocalizedMessage() );
                }
            }
            catch ( final Exception e ) {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "cannot find or open configuration url ["
                        + config_filename + "]" );
            }
        }
        else {
            // Otherwise, open as a file
            File f = new File( config_filename );
            if ( !f.exists() ) {
                f = new File( config_filename + ".txt" );
            }
            if ( f.exists() && f.canRead() ) {
                try {
                    final BufferedReader bf = new BufferedReader( new FileReader( f ) );
                    readConfig( bf );
                    bf.close();
                }
                catch ( final Exception e ) {
                    if ( verbose ) {
                        ForesterUtil.printWarningMessage( Constants.PRG_NAME, "failed to read configuration from ["
                                + config_filename + "]: " + e );
                    }
                }
            }
            else {
                if ( verbose ) {
                    ForesterUtil.printWarningMessage( Constants.PRG_NAME, "cannot find or open configuration file ["
                            + config_filename + "]" );
                }
            }
        }
    }

    public String getBaseFontFamilyName() {
        return _base_font_family_name;
    }

    public int getDefaultBootstrapSamples() {
        return _default_bootstrap_samples;
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

    public Color getDomainStructureBaseColor() {
        return _domain_structure_base_color;
    }

    public Color getDomainStructureFontColor() {
        return _domain_structure_font_color;
    }

    public NODE_DATA getExtDescNodeDataToReturn() {
        return _ext_desc_data_to_return;
    }

    public EXT_NODE_DATA_RETURN_ON getExtNodeDataReturnOn() {
        return _ext_node_data_return_on;
    }

    public int getFrameXSize() {
        return _frame_x_size;
    }

    public int getFrameYSize() {
        return _frame_y_size;
    }

    public String getLabelForGetExtDescendentsData() {
        return _label_for_get_ext_descendents_data;
    }

    public File getPathToLocalFastme() {
        return _path_to_local_fastme;
    }

    public File getpathToLocalMafft() {
        return _path_to_local_mafft;
    }

    public File getPathToLocalRaxml() {
        return _path_to_local_raxml;
    }

    public boolean isAbbreviateScientificTaxonNames() {
        return _abbreviate_scientific_names;
    }

    public boolean isBackgroundColorGradient() {
        return _background_color_gradient;
    }

    public boolean isColorByTaxonomicGroup() {
        return false;
    }

    public boolean isColorLabelsSameAsParentBranch() {
        return _color_labels_same_as_parent_branch;
    }

    public boolean isMidpointReroot() {
        return _midpoint_root;
    }

    public boolean isShowAnnotationRefSource() {
        return _show_annotation_ref_source;
    }

    public boolean isShowDefaultNodeShapesExternal() {
        return _show_default_node_shapes_external;
    }

    public boolean isShowDefaultNodeShapesInternal() {
        return _show_default_node_shapes_internal;
    }

    public boolean isShowDomainLabels() {
        return _show_domain_labels;
    }

    public boolean isTaxonomyColorizeNodeShapes() {
        return _taxonomy_colorize_node_shapes;
    }

    public void putDisplayColors( final String key, final Color color ) {
        getDisplayColors().put( key, color );
    }

    public void setAbbreviateScientificTaxonNames( final boolean abbreviate_scientific_names ) {
        _abbreviate_scientific_names = abbreviate_scientific_names;
    }

    public void setBackgroundColorGradient( final boolean background_color_gradient ) {
        _background_color_gradient = background_color_gradient;
    }

    public void setBaseFontFamilyName( final String base_font_family_name ) {
        _base_font_family_name = base_font_family_name;
    }

    public void setBaseFontSize( final int base_font_size ) {
        _base_font_size = base_font_size;
    }

    public void setColorizeBranches( final boolean b ) {
        display_options[ color_branches ][ 2 ] = b ? "yes" : "no";
    }

    public void setColorLabelsSameAsParentBranch( final boolean color_labels_same_as_parent_branch ) {
        _color_labels_same_as_parent_branch = color_labels_same_as_parent_branch;
    }

    public void setDefaultNodeFill( final NodeFill default_node_fill ) {
        _default_node_fill = default_node_fill;
    }

    public void setDefaultNodeShape( final NodeShape default_node_shape ) {
        _default_node_shape = default_node_shape;
    }

    public void setDefaultNodeShapeSize( final short default_node_shape_size ) {
        _default_node_shape_size = default_node_shape_size;
    }

    public void setDisplayAsPhylogram( final boolean b ) {
        display_options[ display_as_phylogram ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplayColors( final SortedMap<String, Color> display_colors ) {
        _display_colors = display_colors;
    }

    public void setDisplayConfidenceValues( final boolean b ) {
        display_options[ write_confidence_values ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplayInternalData( final boolean b ) {
        display_options[ display_internal_data ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplayNodeNames( final boolean b ) {
        display_options[ show_node_names ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplaySequenceAcc( final boolean b ) {
        display_options[ show_sequence_acc ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplaySequenceNames( final boolean b ) {
        display_options[ show_seq_names ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplaySequenceRelations( final boolean display_sequence_relations ) {
        _display_sequence_relations = display_sequence_relations;
    }

    public void setDisplaySequenceSymbols( final boolean b ) {
        display_options[ show_seq_symbols ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplayTaxonomyCode( final boolean b ) {
        display_options[ show_tax_code ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplayTaxonomyCommonNames( final boolean b ) {
        display_options[ show_taxonomy_common_names ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplayTaxonomyImages( final boolean b ) {
        display_options[ show_taxonomy_images ][ 2 ] = b ? "yes" : "no";
    }

    public void setDisplayTaxonomyScientificNames( final boolean b ) {
        display_options[ show_taxonomy_scientific_names ][ 2 ] = b ? "yes" : "no";
    }

    public void setDynamicallyHideData( final boolean b ) {
        display_options[ dynamically_hide_data ][ 2 ] = b ? "yes" : "no";
    }

    public void setExtDescNodeDataToReturn( final NODE_DATA ext_desc_data_to_return ) {
        _ext_desc_data_to_return = ext_desc_data_to_return;
    }

    public void setFrameXSize( final int frame_x_size ) {
        _frame_x_size = frame_x_size;
    }

    public void setFrameYSize( final int frame_y_size ) {
        _frame_y_size = frame_y_size;
    }

    public void setMidpointReroot( final boolean midpoint_root ) {
        _midpoint_root = midpoint_root;
    }

    public void setMinConfidenceValue( final double min_confidence_value ) {
        _min_confidence_value = min_confidence_value;
    }

    public void setNodeLabelDirection( final NODE_LABEL_DIRECTION node_label_direction ) {
        _node_label_direction = node_label_direction;
    }

    public void setNumberOfDigitsAfterCommaForBranchLengthValue( final short number_of_digits_after_comma_for_branch_length_values ) {
        _number_of_digits_after_comma_for_branch_length_values = number_of_digits_after_comma_for_branch_length_values;
    }

    public void setNumberOfDigitsAfterCommaForConfidenceValues( final short number_of_digits_after_comma_for_confidence_values ) {
        _number_of_digits_after_comma_for_confidence_values = number_of_digits_after_comma_for_confidence_values;
    }

    public void setPhylogenyGraphicsType( final PHYLOGENY_GRAPHICS_TYPE phylogeny_graphics_type ) {
        _phylogeny_graphics_type = phylogeny_graphics_type;
    }

    public void setPrintLineWidth( final float print_line_width ) {
        _print_line_width = print_line_width;
    }

    public void setReplaceUnderscoresInNhParsing( final boolean nh_parsing_replace_underscores ) {
        _nh_parsing_replace_underscores = nh_parsing_replace_underscores;
    }

    public void setShowBranchLengthValues( final boolean show_branch_length_values ) {
        _show_branch_length_values = show_branch_length_values;
    }

    public void setShowDefaultNodeShapesExternal( final boolean show_default_node_shapes_external ) {
        _show_default_node_shapes_external = show_default_node_shapes_external;
    }

    public void setShowDefaultNodeShapesInternal( final boolean show_default_node_shapes_internal ) {
        _show_default_node_shapes_internal = show_default_node_shapes_internal;
    }

    public void setShowDomainLabels( final boolean show_domain_labels ) {
        _show_domain_labels = show_domain_labels;
    }

    public void setShowScale( final boolean show_scale ) {
        _show_scale = show_scale;
    }

    public void setTaxonomyColorize( final boolean b ) {
        display_options[ color_according_to_species ][ 2 ] = b ? "yes" : "no";
    }

    public void setTaxonomyColorizeNodeShapes( final boolean taxonomy_colorize_node_shapes ) {
        _taxonomy_colorize_node_shapes = taxonomy_colorize_node_shapes;
    }

    public void setUseBranchesWidths( final boolean b ) {
        display_options[ width_branches ][ 2 ] = b ? "yes" : "no";
    }

    boolean displaySequenceRelations() {
        return _display_sequence_relations;
    }

    boolean doCheckOption( final int which ) {
        return ( display_options[ which ][ 2 ].equalsIgnoreCase( "yes" ) )
                || ( display_options[ which ][ 2 ].equalsIgnoreCase( "true" ) );
    }

    boolean doDisplayClickToOption( final int which ) {
        return clickto_options[ which ][ 1 ].equalsIgnoreCase( "display" );
    }

    boolean doDisplayOption( final int which ) {
        return display_options[ which ][ 1 ].equalsIgnoreCase( "display" );
    }

    /**
     * Will attempt to use the phylogeny to determine whether to check
     * this or not (e.g. phylogram)
     * 
     */
    boolean doGuessCheckOption( final int which ) {
        return display_options[ which ][ 2 ].equals( "?" );
    }

    Map<String, Color> getAnnotationColors() {
        if ( _annotation_colors == null ) {
            _annotation_colors = new Hashtable<String, Color>();
        }
        return _annotation_colors;
    }

    int getBaseFontSize() {
        return _base_font_size;
    }

    CLADOGRAM_TYPE getCladogramType() {
        return _cladogram_type;
    }

    int getClickToOptionsCount() {
        return clickto_options.length;
    }

    String getClickToTitle( final int which ) {
        return clickto_options[ which ][ 0 ];
    }

    int getDefaultDisplayClicktoOption() {
        return default_clickto;
    }

    SortedMap<String, Color> getDisplayColors() {
        return _display_colors;
    }

    String getDisplayTitle( final int which ) {
        return display_options[ which ][ 0 ];
    }

    Map<String, Color> getDomainColors() {
        if ( _domain_colors == null ) {
            _domain_colors = new Hashtable<String, Color>();
        }
        return _domain_colors;
    }

    int getGraphicsExportX() {
        return _graphics_export_x;
    }

    int getGraphicsExportY() {
        return _graphics_export_y;
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

    int getMaxBaseFontSize() {
        return _max_base_font_size;
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
        if ( _species_colors == null ) {
            _species_colors = new Hashtable<String, Color>();
        }
        return _species_colors;
    }

    final TAXONOMY_EXTRACTION getTaxonomyExtraction() {
        return _taxonomy_extraction;
    }

    boolean isAntialiasScreen() {
        if ( ForesterUtil.isMac() ) {
            //Apple Macintosh graphics are slow, turn off anti-alias.
            return false;
        }
        return _antialias_screen;
    }

    /**
     * Convenience method.
     * 
     * @return true if value in configuration file was 'yes'
     */
    boolean isDrawAsPhylogram() {
        return doCheckOption( display_as_phylogram );
    }

    boolean isEditable() {
        return _editable;
    }

    /**
     * Only used by ArchaeoptryxE.
     *
     */
    boolean isHideControlPanelAndMenubar() {
        return _hide_controls_and_menus;
    }

    boolean isInternalNumberAreConfidenceForNhParsing() {
        return _internal_number_are_confidence_for_nh_parsing;
    }

    boolean isReplaceUnderscoresInNhParsing() {
        return _nh_parsing_replace_underscores;
    }

    boolean isShowBranchLengthValues() {
        return _show_branch_length_values;
    }

    boolean isShowOverview() {
        return _show_overview;
    }

    boolean isShowScale() {
        return _show_scale;
    }

    final boolean isUseNativeUI() {
        if ( ( _ui == UI.UNKNOWN ) && ForesterUtil.isMac() ) {
            _ui = UI.NATIVE;
        }
        return _ui == UI.NATIVE;
    }

    /**
     * Only used by ArchaeoptryxE.
     *
     */
    boolean isUseTabbedDisplay() {
        return _use_tabbed_display;
    }

    boolean isValidatePhyloXmlAgainstSchema() {
        return _validate_against_phyloxml_xsd_schema;
    }

    final void setTaxonomyExtraction( final TAXONOMY_EXTRACTION taxonomy_extraction ) {
        _taxonomy_extraction = taxonomy_extraction;
    }

    private int getClickToIndex( final String name ) {
        int index = -1;
        if ( name.equals( "edit_info" ) ) {
            index = Configuration.display_node_data;
            ForesterUtil
                    .printWarningMessage( Constants.PRG_NAME,
                                          "configuration key [edit_info] is deprecated, use [display node data] instead" );
        }
        else if ( name.equals( "display_node_data" ) ) {
            index = Configuration.display_node_data;
        }
        else if ( name.equals( "collapse_uncollapse" ) ) {
            index = Configuration.collapse_uncollapse;
        }
        else if ( name.equals( "reroot" ) ) {
            index = Configuration.reroot;
        }
        else if ( name.equals( "subtree" ) ) {
            index = Configuration.subtree;
        }
        else if ( name.equals( "swap" ) ) {
            index = Configuration.swap;
        }
        else if ( name.equals( "sort_descendants" ) ) {
            index = Configuration.sort_descendents;
        }
        else if ( name.equals( "get_ext_descendents_data" ) ) {
            index = Configuration.get_ext_desc_data;
        }
        else if ( name.equals( "display_sequences" ) ) {
            ForesterUtil
                    .printWarningMessage( Constants.PRG_NAME, "configuration key [display_sequences] is deprecated" );
            return DEPRECATED;
        }
        else if ( name.equals( "open_seq_web" ) ) {
            index = Configuration.open_seq_web;
        }
        else if ( name.equals( "open_pdb_web" ) ) {
            index = Configuration.open_pdb_web;
        }
        else if ( name.equals( "open_tax_web" ) ) {
            index = Configuration.open_tax_web;
        }
        else if ( name.equals( "blast" ) ) {
            index = Configuration.blast;
        }
        else if ( name.equals( "cut_subtree" ) ) {
            index = Configuration.cut_subtree;
        }
        else if ( name.equals( "copy_subtree" ) ) {
            index = Configuration.copy_subtree;
        }
        else if ( name.equals( "paste_subtree" ) ) {
            index = Configuration.paste_subtree;
        }
        else if ( name.equals( "delete" ) ) {
            index = Configuration.delete_subtree_or_node;
        }
        else if ( name.equals( "add_new_node" ) ) {
            index = Configuration.add_new_node;
        }
        else if ( name.equals( "edit_node_data" ) ) {
            index = Configuration.edit_node_data;
        }
        else if ( name.equals( "select_nodes" ) ) {
            index = Configuration.select_nodes;
        }
        else if ( name.equals( "display_node_popup" ) ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME,
                                              "configuration key [display_node_popup] is deprecated" );
            return DEPRECATED;
        }
        else if ( name.equals( "custom_option" ) ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "configuration key [custom_option] is deprecated" );
            return DEPRECATED;
        }
        else if ( name.equals( "color_subtree" ) ) {
            index = Configuration.color_subtree;
        }
        return index;
    }

    private boolean parseBoolean( final String str ) {
        final String my_str = str.trim().toLowerCase();
        if ( my_str.equals( "yes" ) || my_str.equals( "true" ) ) {
            return true;
        }
        else if ( my_str.equals( "no" ) || my_str.equals( "false" ) ) {
            return false;
        }
        else {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "could not parse boolean value from [" + str + "]" );
            return false;
        }
    }

    private double parseDouble( final String str ) {
        double d = 0.0;
        try {
            d = Double.parseDouble( str.trim() );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "could not parse double from [" + str + "]" );
            d = 0.0;
        }
        return d;
    }

    private float parseFloat( final String str ) {
        float f = 0.0f;
        try {
            f = Float.parseFloat( str.trim() );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "could not parse float from [" + str + "]" );
            f = 0.0f;
        }
        return f;
    }

    private int parseInt( final String str ) {
        int i = -1;
        try {
            i = Integer.parseInt( str.trim() );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "could not parse integer from [" + str + "]" );
            i = -1;
        }
        return i;
    }

    private short parseShort( final String str ) {
        short i = -1;
        try {
            i = Short.parseShort( str.trim() );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "could not parse short from [" + str + "]" );
            i = -1;
        }
        return i;
    }

    private void processFontFamily( final StringTokenizer st ) {
        setBaseFontFamilyName( "" );
        final String font_str = ( ( String ) st.nextElement() ).trim();
        final String[] fonts = font_str.split( ",+" );
        for( String font : fonts ) {
            font = font.replace( '_', ' ' ).trim();
            if ( Arrays.binarySearch( AptxUtil.getAvailableFontFamiliesSorted(), font ) >= 0 ) {
                setBaseFontFamilyName( font );
                break;
            }
        }
    }

    /**
     * read each line of config file, process non-comment lines
     * @throws IOException 
     */
    private void readConfig( final BufferedReader conf_in ) throws IOException {
        String line;
        do {
            line = conf_in.readLine();
            if ( line != null ) {
                line = line.trim();
                // skip comments and blank lines
                if ( !line.startsWith( "#" ) && ( !ForesterUtil.isEmpty( line ) ) ) {
                    // convert runs of spaces to tabs
                    line = line.replaceAll( "\\s+", "\t" );
                    final StringTokenizer st = new StringTokenizer( line, "\t" );
                    setKeyValue( st );
                }
            }
        } while ( line != null );
    }

    private void setAntialiasScreen( final boolean antialias_screen ) {
        _antialias_screen = antialias_screen;
    }

    private void setCladogramType( final CLADOGRAM_TYPE cladogram_type ) {
        _cladogram_type = cladogram_type;
    }

    private void setDefaultBootstrapSamples( final int default_bootstrap_samples ) {
        _default_bootstrap_samples = default_bootstrap_samples;
    }

    private void setEditable( final boolean editable ) {
        _editable = editable;
    }

    private void setExtNodeDataReturnOn( final EXT_NODE_DATA_RETURN_ON ext_node_data_return_on ) {
        _ext_node_data_return_on = ext_node_data_return_on;
    }

    private void setGraphicsExportX( final int graphics_export_x ) {
        _graphics_export_x = graphics_export_x;
    }

    private void setGraphicsExportY( final int graphics_export_y ) {
        _graphics_export_y = graphics_export_y;
    }

    private void setInternalNumberAreConfidenceForNhParsing( final boolean internal_number_are_confidence_for_nh_parsing ) {
        _internal_number_are_confidence_for_nh_parsing = internal_number_are_confidence_for_nh_parsing;
    }

    /**
     * Set a key-value(s) tuple
     */
    private void setKeyValue( final StringTokenizer st ) {
        final String key = ( ( String ) st.nextElement() ).replace( ':', ' ' ).trim().toLowerCase();
        if ( !st.hasMoreElements() ) {
            return;
        }
        // Handle single value settings first:
        if ( key.equals( "default_click_to" ) ) {
            final String clickto_name = ( String ) st.nextElement();
            default_clickto = getClickToIndex( clickto_name );
            if ( default_clickto == -1 ) {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "invalid value [" + clickto_name
                        + "] for [default_click_to]" );
                default_clickto = 0;
            }
            else if ( default_clickto == DEPRECATED ) {
                // Deprecated.
            }
        }
        else if ( key.equals( "native_ui" ) ) {
            final String my_str = ( ( String ) st.nextElement() ).trim().toLowerCase();
            if ( my_str.equals( "yes" ) || my_str.equals( "true" ) ) {
                _ui = UI.NATIVE;
            }
            else if ( my_str.equals( "no" ) || my_str.equals( "false" ) ) {
                _ui = UI.CROSSPLATFORM;
            }
            else if ( my_str.equals( "?" ) ) {
                _ui = UI.UNKNOWN;
            }
            else {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "could not parse yes/no/? value from [" + my_str
                        + "]" );
                _ui = UI.UNKNOWN;
            }
        }
        else if ( key.equals( VALIDATE_AGAINST_PHYLOXML_XSD_SCHEMA ) ) {
            setValidatePhyloXmlAgainstSchema( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "antialias_screen" ) ) {
            setAntialiasScreen( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "phylogeny_graphics_type" ) ) {
            final String type_str = ( ( String ) st.nextElement() ).trim();
            if ( type_str.equalsIgnoreCase( PHYLOGENY_GRAPHICS_TYPE.CONVEX.toString() ) ) {
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CONVEX );
            }
            else if ( type_str.equalsIgnoreCase( PHYLOGENY_GRAPHICS_TYPE.CURVED.toString() ) ) {
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CURVED );
            }
            else if ( type_str.equalsIgnoreCase( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE.toString() ) ) {
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE );
            }
            else if ( type_str.equalsIgnoreCase( PHYLOGENY_GRAPHICS_TYPE.ROUNDED.toString() ) ) {
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.ROUNDED );
            }
            else if ( type_str.equalsIgnoreCase( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR.toString() ) ) {
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            }
            else if ( type_str.equalsIgnoreCase( PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR.toString() ) ) {
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR );
            }
            else if ( type_str.equalsIgnoreCase( PHYLOGENY_GRAPHICS_TYPE.UNROOTED.toString() ) ) {
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
            }
            else if ( type_str.equalsIgnoreCase( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR.toString() ) ) {
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            }
            else {
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown value [" + type_str
                        + "] for [phylogeny_graphics_type]" );
            }
        }
        else if ( key.equals( "min_confidence_value" ) ) {
            final String mcv_str = ( ( String ) st.nextElement() ).trim();
            final double d = parseDouble( mcv_str );
            setMinConfidenceValue( d );
        }
        else if ( key.equals( "font_family" ) ) {
            processFontFamily( st );
        }
        else if ( key.equals( "font_size" ) ) {
            final String size_str = ( ( String ) st.nextElement() ).trim();
            final int i = parseInt( size_str );
            if ( i > 0 ) {
                setBaseFontSize( i );
            }
        }
        else if ( key.equals( "font_size_min" ) ) {
            final String size_str = ( ( String ) st.nextElement() ).trim();
            final int i = parseInt( size_str );
            if ( i > 0 ) {
                setMinBaseFontSize( i );
            }
        }
        else if ( key.equals( "font_size_max" ) ) {
            final String size_str = ( ( String ) st.nextElement() ).trim();
            final int i = parseInt( size_str );
            if ( i > 1 ) {
                setMaxBaseFontSize( i );
            }
        }
        else if ( key.equals( "graphics_export_x" ) ) {
            final String str = ( ( String ) st.nextElement() ).trim();
            final int i = parseInt( str );
            if ( i > 0 ) {
                setGraphicsExportX( i );
            }
        }
        else if ( key.equals( "graphics_export_y" ) ) {
            final String str = ( ( String ) st.nextElement() ).trim();
            final int i = parseInt( str );
            if ( i > 0 ) {
                setGraphicsExportY( i );
            }
        }
        else if ( key.equals( "pdf_export_line_width" ) ) {
            final String str = ( ( String ) st.nextElement() ).trim();
            final float f = parseFloat( str );
            if ( f > 0 ) {
                setPrintLineWidth( f );
            }
            else {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME,
                                                  "value for [pdf_export_line_width] cannot be zero or negative" );
            }
        }
        else if ( key.equals( "window_initial_size_x" ) ) {
            final String str = ( ( String ) st.nextElement() ).trim();
            final int i = parseInt( str );
            if ( i > 0 ) {
                setFrameXSize( i );
            }
        }
        else if ( key.equals( "window_initial_size_y" ) ) {
            final String str = ( ( String ) st.nextElement() ).trim();
            final int i = parseInt( str );
            if ( i > 0 ) {
                setFrameYSize( i );
            }
        }
        else if ( key.equals( "default_number_of_bootstrap_resamples" ) ) {
            final String str = ( ( String ) st.nextElement() ).trim();
            final int i = parseInt( str );
            if ( i >= 0 ) {
                setDefaultBootstrapSamples( i );
            }
            else {
                ForesterUtil
                        .printWarningMessage( Constants.PRG_NAME,
                                              "value for [default_number_of_bootstrap_resamples] cannot be negative" );
            }
        }
        else if ( key.equals( "mafft_local" ) ) {
            final String str = ( ( String ) st.nextElement() ).trim();
            if ( !ForesterUtil.isEmpty( str ) ) {
                setPathToLocalMafft( new File( str ) );
            }
        }
        else if ( key.equals( "fastme_local" ) ) {
            final String str = ( ( String ) st.nextElement() ).trim();
            if ( !ForesterUtil.isEmpty( str ) ) {
                setPathToLocalFastme( new File( str ) );
            }
        }
        else if ( key.equals( "raxml_local" ) ) {
            final String str = ( ( String ) st.nextElement() ).trim();
            if ( !ForesterUtil.isEmpty( str ) ) {
                setPathToLocalRaxml( new File( str ) );
            }
        }
        else if ( key.equals( "show_scale" ) ) {
            setShowScale( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "show_overview" ) ) {
            setShowOverview( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "show_branch_length_values" ) ) {
            setShowBranchLengthValues( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "background_gradient" ) ) {
            setBackgroundColorGradient( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "color_labels_same_as_branch_length_values" ) ) {
            setColorLabelsSameAsParentBranch( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "show_domain_labels" ) ) {
            setShowDomainLabels( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "show_seq_annotation_ref_sources" ) ) {
            setShowAnnotationRefSource( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "abbreviate_scientific_names" ) ) {
            setAbbreviateScientificTaxonNames( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "cladogram_type" ) ) {
            final String type_str = ( ( String ) st.nextElement() ).trim();
            if ( type_str.equalsIgnoreCase( Options.CLADOGRAM_TYPE.NON_LINED_UP.toString() ) ) {
                setCladogramType( Options.CLADOGRAM_TYPE.NON_LINED_UP );
            }
            else if ( type_str.equalsIgnoreCase( Options.CLADOGRAM_TYPE.EXT_NODE_SUM_DEP.toString() ) ) {
                setCladogramType( Options.CLADOGRAM_TYPE.EXT_NODE_SUM_DEP );
            }
            else if ( type_str.equalsIgnoreCase( Options.CLADOGRAM_TYPE.TOTAL_NODE_SUM_DEP.toString() ) ) {
                setCladogramType( Options.CLADOGRAM_TYPE.TOTAL_NODE_SUM_DEP );
            }
            else {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown value [" + type_str
                        + "] for [cladogram_type]" );
            }
        }
        else if ( key.equals( "non_lined_up_cladogram" ) ) {
            ForesterUtil
                    .printWarningMessage( Constants.PRG_NAME,
                                          "configuration key [non_lined_up_cladogram] is deprecated, use [cladogram_type] instead" );
        }
        else if ( key.equals( "hide_controls_and_menus" ) ) {
            _hide_controls_and_menus = parseBoolean( ( String ) st.nextElement() );
        }
        else if ( key.equals( "use_tabbed_display" ) ) {
            _use_tabbed_display = parseBoolean( ( String ) st.nextElement() );
        }
        else if ( key.equals( "overview_width" ) ) {
            final short i = parseShort( ( ( String ) st.nextElement() ) );
            setOvMaxWidth( i );
        }
        else if ( key.equals( "overview_height" ) ) {
            final short i = parseShort( ( ( String ) st.nextElement() ) );
            setOvMaxHeight( i );
        }
        else if ( key.equals( "overview_placement_type" ) ) {
            final String type_str = ( ( String ) st.nextElement() ).trim();
            if ( type_str.equalsIgnoreCase( OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT.toTag() ) ) {
                setOvPlacement( OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT );
            }
            else if ( type_str.equalsIgnoreCase( OVERVIEW_PLACEMENT_TYPE.UPPER_RIGHT.toTag() ) ) {
                setOvPlacement( OVERVIEW_PLACEMENT_TYPE.UPPER_RIGHT );
            }
            else if ( type_str.equalsIgnoreCase( OVERVIEW_PLACEMENT_TYPE.LOWER_LEFT.toTag() ) ) {
                setOvPlacement( OVERVIEW_PLACEMENT_TYPE.LOWER_LEFT );
            }
            else if ( type_str.equalsIgnoreCase( OVERVIEW_PLACEMENT_TYPE.LOWER_RIGHT.toTag() ) ) {
                setOvPlacement( OVERVIEW_PLACEMENT_TYPE.LOWER_RIGHT );
            }
            else {
                setOvPlacement( OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT );
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown value [" + type_str
                        + "] for [overview_placement_type]" );
            }
        }
        else if ( key.equals( "node_label_direction" ) ) {
            final String type_str = ( ( String ) st.nextElement() ).trim();
            if ( type_str.equalsIgnoreCase( NODE_LABEL_DIRECTION.HORIZONTAL.toString() ) ) {
                setNodeLabelDirection( NODE_LABEL_DIRECTION.HORIZONTAL );
            }
            else if ( type_str.equalsIgnoreCase( NODE_LABEL_DIRECTION.RADIAL.toString() ) ) {
                setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
            }
            else {
                setNodeLabelDirection( NODE_LABEL_DIRECTION.HORIZONTAL );
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown value [" + type_str
                        + "] for [node_label_direction]" );
            }
        }
        else if ( key.equals( "branch_length_value_digits" ) ) {
            final short i = parseShort( ( ( String ) st.nextElement() ).trim() );
            if ( i >= 0 ) {
                setNumberOfDigitsAfterCommaForBranchLengthValue( i );
            }
            else {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "illegal value [" + i
                        + "] for [branch_length_value_digits]" );
            }
        }
        else if ( key.equals( "confidence_value_digits" ) ) {
            final short i = parseShort( ( ( String ) st.nextElement() ).trim() );
            if ( i >= 0 ) {
                setNumberOfDigitsAfterCommaForConfidenceValues( i );
            }
            else {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "illegal value [" + i
                        + "] for [confidence_value_digits]" );
            }
        }
        else if ( key.equals( "allow_editing" ) ) {
            setEditable( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "display_sequence_relations" ) ) {
            setDisplaySequenceRelations( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "replace_underscores_in_nh_parsing" ) ) {
            final boolean r = parseBoolean( ( String ) st.nextElement() );
            if ( r && ( getTaxonomyExtraction() != TAXONOMY_EXTRACTION.NO ) ) {
                ForesterUtil
                        .printWarningMessage( Constants.PRG_NAME,
                                              "attempt to extract taxonomies and replace underscores at the same time" );
            }
            else {
                setReplaceUnderscoresInNhParsing( r );
            }
        }
        else if ( key.equals( "taxonomy_extraction_in_nh_parsing" ) ) {
            final String s = ( String ) st.nextElement();
            if ( s.equalsIgnoreCase( "no" ) ) {
                setTaxonomyExtraction( TAXONOMY_EXTRACTION.NO );
            }
            else if ( s.equalsIgnoreCase( "pfam_relaxed" ) ) {
                setTaxonomyExtraction( TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            }
            else if ( s.equalsIgnoreCase( "pfam_strict" ) ) {
                setTaxonomyExtraction( TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            }
            else if ( s.equalsIgnoreCase( "aggressive" ) ) {
                setTaxonomyExtraction( TAXONOMY_EXTRACTION.AGGRESSIVE );
            }
            else {
                ForesterUtil
                        .printWarningMessage( Constants.PRG_NAME,
                                              "unknown value for \"taxonomy_extraction_in_nh_parsing\": "
                                                      + s
                                                      + " (must be either: no, pfam_relaxed, pfam_strict, or aggressive)" );
            }
            if ( ( getTaxonomyExtraction() != TAXONOMY_EXTRACTION.NO ) && isReplaceUnderscoresInNhParsing() ) {
                ForesterUtil
                        .printWarningMessage( Constants.PRG_NAME,
                                              "attempt to extract taxonomies and replace underscores at the same time" );
            }
        }
        else if ( key.equals( "internal_labels_are_confidence_values" ) ) {
            setInternalNumberAreConfidenceForNhParsing( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "gui_background_color" ) ) {
            _gui_background_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "gui_checkbox_text_color" ) ) {
            _gui_checkbox_text_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "gui_checkbox_and_button_active_color" ) ) {
            _gui_checkbox_and_button_active_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "gui_button_text_color" ) ) {
            _gui_button_text_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "gui_button_background_color" ) ) {
            _gui_button_background_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "gui_menu_background_color" ) ) {
            _gui_menu_background_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "gui_menu_text_color" ) ) {
            _gui_menu_text_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "gui_button_border_color" ) ) {
            _gui_button_border_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "domain_structure_font_color" ) ) {
            _domain_structure_font_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "domain_structure_base_color" ) ) {
            _domain_structure_base_color = Color.decode( ( String ) st.nextElement() );
        }
        else if ( key.equals( "show_default_node_shapes" ) ) {
            ForesterUtil
                    .printWarningMessage( Constants.PRG_NAME,
                                          "configuration key [show_default_node_shapes] is deprecated, use [show_default_node_shapes_internal] and [show_default_node_shapes_external] instead" );
            final boolean b = parseBoolean( ( ( String ) st.nextElement() ).trim() );
            setShowDefaultNodeShapesInternal( b );
            setShowDefaultNodeShapesExternal( b );
        }
        else if ( key.equals( "show_default_node_shapes_internal" ) ) {
            setShowDefaultNodeShapesInternal( parseBoolean( ( ( String ) st.nextElement() ).trim() ) );
        }
        else if ( key.equals( "show_default_node_shapes_external" ) ) {
            setShowDefaultNodeShapesExternal( parseBoolean( ( ( String ) st.nextElement() ).trim() ) );
        }
        else if ( key.equals( "default_node_size" ) ) {
            final short i = parseShort( ( ( String ) st.nextElement() ).trim() );
            setDefaultNodeShapeSize( i );
        }
        else if ( key.equals( "default_node_fill" ) ) {
            final String fill_str = ( ( String ) st.nextElement() ).trim();
            if ( fill_str.equalsIgnoreCase( NodeVisualization.NodeFill.NONE.toString() ) ) {
                setDefaultNodeFill( NodeFill.NONE );
            }
            else if ( fill_str.equalsIgnoreCase( NodeVisualization.NodeFill.GRADIENT.toString() ) ) {
                setDefaultNodeFill( NodeFill.GRADIENT );
            }
            else if ( fill_str.equalsIgnoreCase( NodeVisualization.NodeFill.SOLID.toString() ) ) {
                setDefaultNodeFill( NodeFill.SOLID );
            }
            else {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown value [" + fill_str
                        + "] for [default_node_fill]" );
            }
        }
        else if ( key.equals( "default_node_shape" ) ) {
            final String shape_str = ( ( String ) st.nextElement() ).trim();
            if ( shape_str.equalsIgnoreCase( NodeVisualization.NodeShape.CIRCLE.toString() ) ) {
                setDefaultNodeShape( NodeShape.CIRCLE );
            }
            else if ( shape_str.equalsIgnoreCase( NodeVisualization.NodeShape.RECTANGLE.toString() ) ) {
                setDefaultNodeShape( NodeShape.RECTANGLE );
            }
            else {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown value [" + shape_str
                        + "] for [default_node_shape]" );
            }
        }
        else if ( key.equals( "taxonomy_colorize_node_shapes" ) ) {
            setTaxonomyColorizeNodeShapes( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "midpoint_reroot" ) ) {
            setMidpointReroot( parseBoolean( ( String ) st.nextElement() ) );
        }
        else if ( key.equals( "ext_descendents_data_to_return" ) ) {
            final String s = ( ( String ) st.nextElement() ).trim();
            if ( s.equalsIgnoreCase( "node_name" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.NODE_NAME );
            }
            else if ( s.equalsIgnoreCase( "sequence_acc" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.SEQUENCE_ACC );
            }
            else if ( s.equalsIgnoreCase( "sequence_mol_seq_fasta" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.SEQUENCE_MOL_SEQ_FASTA );
            }
            else if ( s.equalsIgnoreCase( "sequence_mol_seq" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.SEQUENCE_MOL_SEQ );
            }
            else if ( s.equalsIgnoreCase( "sequence_name" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.SEQUENCE_NAME );
            }
            else if ( s.equalsIgnoreCase( "gene_name" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.GENE_NAME );
            }
            else if ( s.equalsIgnoreCase( "sequence_symbol" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.SEQUENCE_SYMBOL );
            }
            else if ( s.equalsIgnoreCase( "taxonomy_scientific_name" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.TAXONOMY_SCIENTIFIC_NAME );
            }
            else if ( s.equalsIgnoreCase( "taxonomy_code" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.TAXONOMY_CODE );
            }
            else if ( s.equalsIgnoreCase( "taxonomy_common_name" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.TAXONOMY_COMM0N_NAME );
            }
            else if ( s.equalsIgnoreCase( "user_selected" ) ) {
                setExtDescNodeDataToReturn( NODE_DATA.UNKNOWN );
            }
            else {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown value [" + s
                        + "] for [ext_descendents_data_to_return]" );
            }
        }
        else if ( key.equals( "label_for_get_ext_descendents_data" ) ) {
            final String s = ( ( String ) st.nextElement() ).trim();
            if ( ForesterUtil.isEmpty( s ) || ( s.length() < 2 ) ) {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "illegal value [" + s
                        + "] for [label_for_get_ext_descendents_data]" );
            }
            else {
                setLabelForGetExtDescendentsData( s.replaceAll( "_", " " ) );
            }
        }
        else if ( key.equals( "ext_descendents_data_to_return_on" ) ) {
            final String s = ( ( String ) st.nextElement() ).trim().toLowerCase();
            if ( s.equals( "console" ) ) {
                setExtNodeDataReturnOn( EXT_NODE_DATA_RETURN_ON.CONSOLE );
            }
            else if ( s.equals( "window" ) ) {
                setExtNodeDataReturnOn( EXT_NODE_DATA_RETURN_ON.WINODW );
            }
            else if ( s.equals( "buffer_only" ) ) {
                setExtNodeDataReturnOn( EXT_NODE_DATA_RETURN_ON.BUFFER_ONLY );
            }
            else {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown value [" + s
                        + "] for [ext_descendents_data_to_return_on]" );
            }
        }
        else if ( st.countTokens() >= 2 ) { // counts the tokens that are not
            // yet retrieved!
            int key_index = -1;
            if ( key.equals( "phylogram" ) ) {
                key_index = Configuration.display_as_phylogram;
            }
            else if ( key.equals( "rollover" ) ) {
                key_index = Configuration.node_data_popup;
            }
            else if ( key.equals( "color_according_to_species" ) ) {
                key_index = Configuration.color_according_to_species;
            }
            else if ( key.equals( "show_node_names" ) ) {
                key_index = Configuration.show_node_names;
            }
            else if ( key.equals( "show_taxonomy_code" ) ) {
                key_index = Configuration.show_tax_code;
            }
            else if ( key.equals( "write_confidence_values" ) ) {
                key_index = Configuration.write_confidence_values;
            }
            else if ( key.equals( "write_events" ) ) {
                key_index = Configuration.write_events;
            }
            else if ( key.equals( "color_branches" ) ) {
                key_index = Configuration.color_branches;
            }
            else if ( key.equals( "width_branches" ) ) {
                key_index = Configuration.width_branches;
            }
            else if ( key.equals( "show_domain_architectures" ) ) {
                key_index = Configuration.show_domain_architectures;
            }
            else if ( key.equals( "show_annotations" ) ) {
                key_index = Configuration.show_annotation;
            }
            else if ( key.equals( "show_binary_characters" ) ) {
                key_index = Configuration.show_binary_characters;
            }
            else if ( key.equals( "show_binary_character_counts" ) ) {
                key_index = Configuration.show_binary_character_counts;
            }
            else if ( key.equals( "show_seq_names" ) ) {
                key_index = Configuration.show_seq_names;
            }
            else if ( key.equals( "show_gene_names" ) ) {
                key_index = Configuration.show_gene_names;
            }
            else if ( key.equals( "show_seq_symbols" ) ) {
                key_index = Configuration.show_seq_symbols;
            }
            else if ( key.equals( "show_seq_acc" ) ) {
                key_index = Configuration.show_sequence_acc;
            }
            else if ( key.equals( "display_internal_data" ) ) {
                key_index = Configuration.display_internal_data;
            }
            else if ( key.equals( "dynamically_hide_data" ) ) {
                key_index = Configuration.dynamically_hide_data;
            }
            else if ( key.equals( "show_taxonomy_scientific_names" ) ) {
                key_index = Configuration.show_taxonomy_scientific_names;
            }
            else if ( key.equals( "show_taxonomy_common_names" ) ) {
                key_index = Configuration.show_taxonomy_common_names;
            }
            else if ( key.equals( "show_taxonomy_images" ) ) {
                key_index = Configuration.show_taxonomy_images;
            }
            else if ( key.equals( "color_according_to_annotation" ) ) {
                key_index = Configuration.color_according_to_annotation;
            }
            else if ( key.equals( "show_vector_data" ) ) {
                key_index = Configuration.show_vector_data;
            }
            else if ( key.equals( "show_properties" ) ) {
                key_index = Configuration.show_properties;
            }
            else if ( key.equals( "show_relation_confidence" ) ) {
                key_index = Configuration.show_relation_confidence;
            }
            else if ( key.equals( "show_custom_node_shapes" ) ) {
                key_index = Configuration.show_custom_node_shapes;
            }
            // If we've found the key, set the values
            if ( key_index >= 0 ) {
                display_options[ key_index ][ 1 ] = ( String ) st.nextElement();
                display_options[ key_index ][ 2 ] = ( String ) st.nextElement();
                // otherwise, keep looking
            }
            else {
                if ( key_index == DEPRECATED ) {
                    // Deprecated.
                }
                else if ( key.equals( "click_to" ) ) {
                    final String click_to_name = ( String ) st.nextElement();
                    key_index = getClickToIndex( click_to_name );
                    if ( key_index >= 0 ) {
                        clickto_options[ key_index ][ 1 ] = ( String ) st.nextElement();
                    }
                    else if ( key_index == DEPRECATED ) {
                        // Deprecated.
                    }
                    else {
                        ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown click-to option: "
                                + click_to_name );
                    }
                }
                else if ( key.equals( "species_color" ) ) {
                    getSpeciesColors().put( ( ( String ) st.nextElement() ).replace( '_', ' ' ),
                                            Color.decode( ( String ) st.nextElement() ) );
                }
                else if ( key.equals( "domain_color" ) ) {
                    getDomainColors().put( ( String ) st.nextElement(), Color.decode( ( String ) st.nextElement() ) );
                }
                else if ( key.equals( "annotation_color" ) ) {
                    getAnnotationColors()
                            .put( ( String ) st.nextElement(), Color.decode( ( String ) st.nextElement() ) );
                }
                else if ( key.equals( "function_color" ) ) {
                    ForesterUtil.printWarningMessage( Constants.PRG_NAME,
                                                      "configuration key [function_color] is deprecated" );
                }
                else if ( key.equals( DISPLAY_COLOR_KEY ) ) {
                    putDisplayColors( ( String ) st.nextElement(), Color.decode( ( String ) st.nextElement() ) );
                }
                else {
                    ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown configuration key [" + key
                            + "] in: " + config_filename );
                }
            }
        }
        else {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "unknown configuration key [" + key + "] in: "
                    + config_filename );
        }
    }

    private void setLabelForGetExtDescendentsData( final String label_for_get_ext_descendents_data ) {
        _label_for_get_ext_descendents_data = label_for_get_ext_descendents_data;
    }

    private void setMaxBaseFontSize( final int max_base_font_size ) {
        _max_base_font_size = max_base_font_size;
    }

    private void setMinBaseFontSize( final int min_base_font_size ) {
        _min_base_font_size = min_base_font_size;
    }

    private void setOvMaxHeight( final short ov_max_height ) {
        _ov_max_height = ov_max_height;
    }

    private void setOvMaxWidth( final short ov_max_width ) {
        _ov_max_width = ov_max_width;
    }

    private void setOvPlacement( final OVERVIEW_PLACEMENT_TYPE ov_placement ) {
        _ov_placement = ov_placement;
    }

    private void setPathToLocalFastme( final File path_to_local_fastme ) {
        _path_to_local_fastme = path_to_local_fastme;
    }

    private void setPathToLocalMafft( final File path_to_local_mafft ) {
        _path_to_local_mafft = path_to_local_mafft;
    }

    private void setPathToLocalRaxml( final File path_to_local_raxml ) {
        _path_to_local_raxml = path_to_local_raxml;
    }

    private void setShowAnnotationRefSource( final boolean b ) {
        _show_annotation_ref_source = b;
    }

    private void setShowOverview( final boolean show_overview ) {
        _show_overview = show_overview;
    }

    private void setValidatePhyloXmlAgainstSchema( final boolean validate_against_phyloxml_xsd_schema ) {
        _validate_against_phyloxml_xsd_schema = validate_against_phyloxml_xsd_schema;
    }

    static String getDefaultFontFamilyName() {
        return DEFAULT_FONT_FAMILY;
    }

    public enum EXT_NODE_DATA_RETURN_ON {
        BUFFER_ONLY, CONSOLE, WINODW;
    }

    public enum UI {
        CROSSPLATFORM, NATIVE, NIMBUS, UNKNOWN
    }

    static enum TRIPLET {
        FALSE, TRUE, UNKNOWN
    }
}
