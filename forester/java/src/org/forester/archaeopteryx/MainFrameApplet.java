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

import java.awt.BorderLayout;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

import javax.swing.ButtonGroup;
import javax.swing.JApplet;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;

public final class MainFrameApplet extends MainFrame {

    private static final long    serialVersionUID     = 1941019292746717053L;
    private final static int     DEFAULT_FRAME_X_SIZE = 640;
    private final static int     DEFAULT_FRAME_Y_SIZE = 580;
    private final ArchaeopteryxA _applet;
    private ButtonGroup          _radio_group_1;

    MainFrameApplet( final ArchaeopteryxA parent_applet,
                     final Configuration configuration,
                     final String species_tree_url_str ) {
        setTitle( ArchaeopteryxA.NAME );
        _applet = parent_applet;
        setConfiguration( configuration );
        setOptions( Options.createInstance( configuration ) );
        _mainpanel = new MainPanelApplets( _configuration, this );
        if ( !ForesterUtil.isEmpty( species_tree_url_str ) ) {
            try {
                readSpeciesTree( configuration, species_tree_url_str );
            }
            catch ( final Exception e ) {
                ForesterUtil.printErrorMessage( ArchaeopteryxA.NAME, "failed to read species tree from "
                        + species_tree_url_str );
                ForesterUtil.printErrorMessage( ArchaeopteryxA.NAME, e.toString() );
            }
        }
        // build the menu bar
        _jmenubar = new JMenuBar();
        if ( !_configuration.isUseNativeUI() ) {
            _jmenubar.setBackground( _configuration.getGuiMenuBackgroundColor() );
        }
        if ( getSpeciesTree() != null ) {
            buildAnalysisMenu();
        }
        buildToolsMenu();
        buildViewMenu();
        buildFontSizeMenu();
        buildOptionsMenu();
        buildTypeMenu();
        buildHelpMenu();
        setJMenuBar( _jmenubar );
        _contentpane = getContentPane();
        _contentpane.setLayout( new BorderLayout() );
        _contentpane.add( _mainpanel, BorderLayout.CENTER );
        setSize( getConfiguration().getFrameXSize() > 40 ? getConfiguration().getFrameXSize() : DEFAULT_FRAME_X_SIZE,
                 getConfiguration().getFrameYSize() > 40 ? getConfiguration().getFrameYSize() : DEFAULT_FRAME_Y_SIZE );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                close();
            }
        } );
        addComponentListener( new ComponentAdapter() {

            @Override
            public void componentResized( final ComponentEvent e ) {
                if ( _mainpanel.getCurrentTreePanel() != null ) {
                    _mainpanel.getCurrentTreePanel().calcParametersForPainting( _mainpanel.getCurrentTreePanel()
                                                                                        .getWidth(),
                                                                                _mainpanel.getCurrentTreePanel()
                                                                                        .getHeight() );
                }
            }
        } );
        setFocusable( true );
        requestFocus();
        requestFocusInWindow();
        setVisible( true );
        System.gc();
    }

    private void readSpeciesTree( final Configuration configuration, final String species_tree_url_str )
            throws MalformedURLException, FileNotFoundException, IOException {
        final URL species_tree_url = new URL( species_tree_url_str );
        final Phylogeny[] species_trees = AptxUtil.readPhylogeniesFromUrl( species_tree_url,
                                                                           configuration
                                                                                   .isValidatePhyloXmlAgainstSchema(),
                                                                           configuration
                                                                                   .isReplaceUnderscoresInNhParsing(),
                                                                           false,
                                                                           TAXONOMY_EXTRACTION.NO,
                                                                           false );
        if ( ( species_trees != null ) && ( species_trees.length > 0 ) ) {
            AptxUtil.printAppletMessage( ArchaeopteryxA.NAME, "successfully read species tree" );
            if ( species_trees[ 0 ].isEmpty() ) {
                ForesterUtil.printErrorMessage( ArchaeopteryxA.NAME, "species tree is empty" );
            }
            else if ( !species_trees[ 0 ].isRooted() ) {
                ForesterUtil.printErrorMessage( ArchaeopteryxA.NAME, "species tree is not rooted" );
            }
            else {
                setSpeciesTree( species_trees[ 0 ] );
                AptxUtil.printAppletMessage( ArchaeopteryxA.NAME, "species tree OK" );
            }
        }
        else {
            ForesterUtil.printErrorMessage( ArchaeopteryxA.NAME, "failed to read species tree from "
                    + species_tree_url_str );
        }
    }

    @Override
    public MainPanel getMainPanel() {
        return _mainpanel;
    }

    void buildAnalysisMenu() {
        _analysis_menu = MainFrame.createMenu( "Analysis", getConfiguration() );
        _analysis_menu.add( _gsdi_item = new JMenuItem( "GSDI (Generalized Speciation Duplication Inference)" ) );
        _analysis_menu.add( _gsdir_item = new JMenuItem( "GSDIR (GSDI with re-rooting)" ) );
        customizeJMenuItem( _gsdi_item );
        customizeJMenuItem( _gsdir_item );
        //  _analysis_menu.addSeparator();
        //  _analysis_menu.add( _lineage_inference = new JMenuItem( INFER_ANCESTOR_TAXONOMIES ) );
        //  customizeJMenuItem( _lineage_inference );
        //  _lineage_inference.setToolTipText( "Inference of ancestor taxonomies/lineages" );
        _jmenubar.add( _analysis_menu );
    }

    void buildOptionsMenu() {
        _options_jmenu = MainFrame.createMenu( MainFrame.OPTIONS_HEADER, getConfiguration() );
        _options_jmenu.addChangeListener( new ChangeListener() {

            @Override
            public void stateChanged( final ChangeEvent e ) {
                MainFrame.setOvPlacementColorChooseMenuItem( _overview_placment_mi, getOptions() );
                MainFrame.setTextColorChooseMenuItem( _switch_colors_mi, getCurrentTreePanel() );
                MainFrame
                        .setTextMinSupportMenuItem( _choose_minimal_confidence_mi, getOptions(), getCurrentTreePanel() );
                MainFrame.setTextForFontChooserMenuItem( _choose_font_mi, createCurrentFontDesc( getMainPanel()
                        .getTreeFontSet() ) );
                MainFrame.setCycleNodeFillMenuItem( _cycle_node_fill_mi, getOptions() );
                MainFrame.setCycleNodeShapeMenuItem( _cycle_node_shape_mi, getOptions() );
                MainFrame.setTextNodeSizeMenuItem( _choose_node_size_mi, getOptions() );
                try {
                    getMainPanel().getControlPanel().setVisibilityOfDomainStrucureCB();
                    getMainPanel().getControlPanel().setVisibilityOfX();
                }
                catch ( final Exception ignore ) {
                    // do nothing, not important.
                }
            }
        } );
        _options_jmenu.add( MainFrame.customizeMenuItemAsLabel( new JMenuItem( MainFrame.DISPLAY_SUBHEADER ),
                                                                getConfiguration() ) );
        _options_jmenu
                .add( _ext_node_dependent_cladogram_rbmi = new JRadioButtonMenuItem( MainFrame.NONUNIFORM_CLADOGRAMS_LABEL ) );
        _options_jmenu.add( _uniform_cladograms_rbmi = new JRadioButtonMenuItem( MainFrame.UNIFORM_CLADOGRAMS_LABEL ) );
        _options_jmenu.add( _non_lined_up_cladograms_rbmi = new JRadioButtonMenuItem( NON_LINED_UP_CLADOGRAMS_LABEL ) );
        _radio_group_1 = new ButtonGroup();
        _radio_group_1.add( _ext_node_dependent_cladogram_rbmi );
        _radio_group_1.add( _uniform_cladograms_rbmi );
        _radio_group_1.add( _non_lined_up_cladograms_rbmi );
        _options_jmenu.add( _show_overview_cbmi = new JCheckBoxMenuItem( MainFrame.SHOW_OVERVIEW_LABEL ) );
        _options_jmenu.add( _show_scale_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_SCALE_LABEL ) );
        _options_jmenu
                .add( _show_default_node_shapes_internal_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_NODE_BOXES_LABEL_INT ) );
        _options_jmenu
                .add( _show_default_node_shapes_external_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_NODE_BOXES_LABEL_EXT ) );
        _options_jmenu
                .add( _show_default_node_shapes_for_marked_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_NODE_BOXES_LABEL_MARKED ) );
        _options_jmenu.add( _line_up_renderable_data_cbmi = new JCheckBoxMenuItem( MainFrame.LINE_UP_RENDERABLE_DATA ) );
        if ( getConfiguration().doDisplayOption( Configuration.show_domain_architectures ) ) {
            _options_jmenu.add( _right_line_up_domains_cbmi = new JCheckBoxMenuItem( MainFrame.RIGHT_LINE_UP_DOMAINS ) );
            _options_jmenu.add( _show_domain_labels = new JCheckBoxMenuItem( MainFrame.SHOW_DOMAIN_LABELS_LABEL ) );
        }
        _options_jmenu.add( _show_annotation_ref_source = new JCheckBoxMenuItem( MainFrame.SHOW_ANN_REF_SOURCE_LABEL ) );
        _options_jmenu.add( _show_confidence_stddev_cbmi = new JCheckBoxMenuItem( MainFrame.SHOW_CONF_STDDEV_LABEL ) );
        _options_jmenu
                .add( _color_by_taxonomic_group_cbmi = new JCheckBoxMenuItem( MainFrame.COLOR_BY_TAXONOMIC_GROUP ) );
        _options_jmenu
                .add( _color_labels_same_as_parent_branch = new JCheckBoxMenuItem( MainFrame.COLOR_LABELS_LABEL ) );
        _color_labels_same_as_parent_branch.setToolTipText( MainFrame.COLOR_LABELS_TIP );
        _options_jmenu.add( _abbreviate_scientific_names = new JCheckBoxMenuItem( MainFrame.ABBREV_SN_LABEL ) );
        _options_jmenu.add( _label_direction_cbmi = new JCheckBoxMenuItem( MainFrame.LABEL_DIRECTION_LABEL ) );
        _label_direction_cbmi.setToolTipText( MainFrame.LABEL_DIRECTION_TIP );
        _options_jmenu.add( _screen_antialias_cbmi = new JCheckBoxMenuItem( MainFrame.SCREEN_ANTIALIAS_LABEL ) );
        _options_jmenu.add( _background_gradient_cbmi = new JCheckBoxMenuItem( MainFrame.BG_GRAD_LABEL ) );
        _options_jmenu.add( _cycle_node_shape_mi = new JMenuItem( MainFrame.CYCLE_NODE_SHAPE_LABEL ) );
        _options_jmenu.add( _cycle_node_fill_mi = new JMenuItem( MainFrame.CYCLE_NODE_FILL_LABEL ) );
        _options_jmenu.add( _choose_node_size_mi = new JMenuItem( MainFrame.CHOOSE_NODE_SIZE_LABEL ) );
        _options_jmenu.add( _choose_minimal_confidence_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _overview_placment_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _switch_colors_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _choose_font_mi = new JMenuItem( "" ) );
        _options_jmenu.addSeparator();
        _options_jmenu.add( MainFrame.customizeMenuItemAsLabel( new JMenuItem( MainFrame.SEARCH_SUBHEADER ),
                                                                getConfiguration() ) );
        _options_jmenu
                .add( _search_case_senstive_cbmi = new JCheckBoxMenuItem( MainFrame.SEARCH_CASE_SENSITIVE_LABEL ) );
        _options_jmenu.add( _search_whole_words_only_cbmi = new JCheckBoxMenuItem( MainFrame.SEARCH_TERMS_ONLY_LABEL ) );
        _options_jmenu.add( _search_with_regex_cbmi = new JCheckBoxMenuItem( MainFrame.SEARCH_REGEX_LABEL ) );
        _search_with_regex_cbmi.setToolTipText( MainFrame.SEARCH_WITH_REGEX_TIP );
        _options_jmenu.add( _inverse_search_result_cbmi = new JCheckBoxMenuItem( INVERSE_SEARCH_RESULT_LABEL ) );
        customizeJMenuItem( _choose_font_mi );
        customizeJMenuItem( _switch_colors_mi );
        customizeJMenuItem( _choose_minimal_confidence_mi );
        customizeJMenuItem( _overview_placment_mi );
        customizeCheckBoxMenuItem( _show_default_node_shapes_internal_cbmi, getOptions()
                .isShowDefaultNodeShapesInternal() );
        customizeCheckBoxMenuItem( _show_default_node_shapes_external_cbmi, getOptions()
                .isShowDefaultNodeShapesExternal() );
        customizeCheckBoxMenuItem( _show_default_node_shapes_for_marked_cbmi, getOptions()
                .isShowDefaultNodeShapesForMarkedNodes() );
        customizeJMenuItem( _cycle_node_shape_mi );
        customizeJMenuItem( _cycle_node_fill_mi );
        customizeJMenuItem( _choose_node_size_mi );
        customizeCheckBoxMenuItem( _color_by_taxonomic_group_cbmi, getOptions().isColorByTaxonomicGroup() );
        customizeCheckBoxMenuItem( _color_labels_same_as_parent_branch, getOptions().isColorLabelsSameAsParentBranch() );
        customizeCheckBoxMenuItem( _screen_antialias_cbmi, getOptions().isAntialiasScreen() );
        customizeCheckBoxMenuItem( _background_gradient_cbmi, getOptions().isBackgroundColorGradient() );
        customizeCheckBoxMenuItem( _show_domain_labels, getOptions().isShowDomainLabels() );
        customizeCheckBoxMenuItem( _show_annotation_ref_source, getOptions().isShowAnnotationRefSource() );
        customizeCheckBoxMenuItem( _abbreviate_scientific_names, getOptions().isAbbreviateScientificTaxonNames() );
        customizeCheckBoxMenuItem( _search_case_senstive_cbmi, getOptions().isSearchCaseSensitive() );
        customizeCheckBoxMenuItem( _show_scale_cbmi, getOptions().isShowScale() );
        customizeRadioButtonMenuItem( _non_lined_up_cladograms_rbmi,
                                      getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP );
        customizeRadioButtonMenuItem( _uniform_cladograms_rbmi,
                                      getOptions().getCladogramType() == CLADOGRAM_TYPE.TOTAL_NODE_SUM_DEP );
        customizeRadioButtonMenuItem( _ext_node_dependent_cladogram_rbmi,
                                      getOptions().getCladogramType() == CLADOGRAM_TYPE.EXT_NODE_SUM_DEP );
        customizeCheckBoxMenuItem( _show_overview_cbmi, getOptions().isShowOverview() );
        customizeCheckBoxMenuItem( _label_direction_cbmi,
                                   getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL );
        customizeCheckBoxMenuItem( _search_with_regex_cbmi, getOptions().isSearchWithRegex() );
        customizeCheckBoxMenuItem( _search_whole_words_only_cbmi, getOptions().isMatchWholeTermsOnly() );
        customizeCheckBoxMenuItem( _inverse_search_result_cbmi, getOptions().isInverseSearchResult() );
        customizeCheckBoxMenuItem( _show_confidence_stddev_cbmi, getOptions().isShowConfidenceStddev() );
        customizeCheckBoxMenuItem( _line_up_renderable_data_cbmi, getOptions().isLineUpRendarableNodeData() );
        customizeCheckBoxMenuItem( _right_line_up_domains_cbmi, getOptions().isRightLineUpDomains() );
        _jmenubar.add( _options_jmenu );
    }

    void buildToolsMenu() {
        _tools_menu = MainFrame.createMenu( "Tools", getConfiguration() );
        _tools_menu.add( _confcolor_item = new JMenuItem( "Colorize Branches Depending on Confidence" ) );
        customizeJMenuItem( _confcolor_item );
        _tools_menu.add( _taxcolor_item = new JMenuItem( "Taxonomy Colorize Branches" ) );
        customizeJMenuItem( _taxcolor_item );
        _tools_menu.addSeparator();
        _tools_menu.add( _remove_visual_styles_item = new JMenuItem( "Delete All Visual Styles From Nodes" ) );
        _remove_visual_styles_item
                .setToolTipText( "To remove all node visual styles (fonts, colors) from the current phylogeny." );
        customizeJMenuItem( _remove_visual_styles_item );
        _tools_menu.add( _remove_branch_color_item = new JMenuItem( "Delete All Colors From Branches" ) );
        _remove_branch_color_item.setToolTipText( "To remove all branch color values from the current phylogeny." );
        customizeJMenuItem( _remove_branch_color_item );
        _tools_menu.addSeparator();
        _tools_menu.add( _midpoint_root_item = new JMenuItem( "Midpoint-Root" ) );
        customizeJMenuItem( _midpoint_root_item );
        _tools_menu.addSeparator();
        _tools_menu.add( _collapse_species_specific_subtrees = new JMenuItem( "Collapse Species-Specific Subtrees" ) );
        customizeJMenuItem( _collapse_species_specific_subtrees );
        _jmenubar.add( _tools_menu );
    }

    JApplet getApplet() {
        return _applet;
    }
}
