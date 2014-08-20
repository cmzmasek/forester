
package org.forester.archaeopteryx;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;

import javax.swing.ButtonGroup;
import javax.swing.JApplet;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.apache.commons.codec.binary.Base64;
import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.data.SequenceRelation;
import org.forester.sdi.GSDI;
import org.forester.sdi.GSDIR;
import org.forester.sdi.SDIException;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

// Use like this:
// <applet archive="forester.jar"
// code="org.forester.archaeopteryx.ArchaeopteryxE.class"
// codebase="http://www.myserver.org/path/to/forester"
// width="600"
// height="500"
// alt="ArchaeopteryxE is not working on your system (requires at least Sun Java 1.5)!">
// <param name="url_of_tree_to_load"
// value="http://www.myserver.org/examples/data/apaf.xml">
// <param name="config_file"
// value="http://www.myserver.org/examples/config/config_file.txt">
// </applet>
public class ArchaeopteryxE extends JApplet implements ActionListener {

    private final static String         NAME             = "ArchaeopteryxE";
    private static final long           serialVersionUID = -1220055577935759443L;
    private Configuration               _configuration;
    private MainPanelApplets            _mainpanel;
    private JMenuBar                    _jmenubar;
    private JMenu                       _options_jmenu;
    private JMenu                       _font_size_menu;
    private JMenuItem                   _super_tiny_fonts_mi;
    private JMenuItem                   _tiny_fonts_mi;
    private JMenuItem                   _small_fonts_mi;
    private JMenuItem                   _medium_fonts_mi;
    private JMenuItem                   _large_fonts_mi;
    private JMenu                       _tools_menu;
    private JMenuItem                   _taxcolor_item;
    private JMenuItem                   _confcolor_item;
    private JMenuItem                   _midpoint_root_item;
    private JMenu                       _view_jmenu;
    private JMenuItem                   _view_as_XML_item;
    private JMenuItem                   _view_as_NH_item;
    private JMenuItem                   _view_as_nexus_item;
    private JMenuItem                   _display_basic_information_item;
    private JMenu                       _type_menu;
    private JCheckBoxMenuItem           _rectangular_type_cbmi;
    private JCheckBoxMenuItem           _triangular_type_cbmi;
    private JCheckBoxMenuItem           _curved_type_cbmi;
    private JCheckBoxMenuItem           _convex_type_cbmi;
    private JCheckBoxMenuItem           _euro_type_cbmi;
    private JCheckBoxMenuItem           _rounded_type_cbmi;
    private JCheckBoxMenuItem           _unrooted_type_cbmi;
    private JCheckBoxMenuItem           _circular_type_cbmi;
    private JMenuItem                   _help_item;
    private JMenuItem                   _about_item;
    private JMenu                       _help_jmenu;
    private JMenuItem                   _website_item;
    private JMenuItem                   _phyloxml_website_item;
    private JMenuItem                   _phyloxml_ref_item;
    private JMenuItem                   _aptx_ref_item;
    private JMenuItem                   _remove_branch_color_item;
    private JMenuItem                   _remove_visual_styles_item;
    private JCheckBoxMenuItem           _show_domain_labels;
    private JCheckBoxMenuItem           _show_annotation_ref_source;
    private JCheckBoxMenuItem           _color_labels_same_as_parent_branch;
    private JCheckBoxMenuItem           _abbreviate_scientific_names;
    private JCheckBoxMenuItem           _screen_antialias_cbmi;
    private JCheckBoxMenuItem           _background_gradient_cbmi;
    private JCheckBoxMenuItem           _color_by_taxonomic_group_cbmi;
    private JRadioButtonMenuItem        _non_lined_up_cladograms_rbmi;
    private JRadioButtonMenuItem        _uniform_cladograms_rbmi;
    private JRadioButtonMenuItem        _ext_node_dependent_cladogram_rbmi;
    private Options                     _options;
    private JMenuItem                   _choose_font_mi;
    private JMenuItem                   _switch_colors_mi;
    JCheckBoxMenuItem                   _label_direction_cbmi;
    private JCheckBoxMenuItem           _show_scale_cbmi;
    private JCheckBoxMenuItem           _search_case_senstive_cbmi;
    private JCheckBoxMenuItem           _search_whole_words_only_cbmi;
    private JCheckBoxMenuItem           _inverse_search_result_cbmi;
    private JCheckBoxMenuItem           _show_overview_cbmi;
    private JMenuItem                   _choose_minimal_confidence_mi;
    private JCheckBoxMenuItem           _show_branch_length_values_cbmi;
    private JMenuItem                   _collapse_species_specific_subtrees;
    private JMenuItem                   _overview_placment_mi;
    private ButtonGroup                 _radio_group_1;
    private JCheckBoxMenuItem           _show_default_node_shapes_internal_cbmi;
    private JCheckBoxMenuItem           _show_default_node_shapes_external_cbmi;
    private JMenuItem                   _cycle_node_shape_mi;
    private JMenuItem                   _cycle_node_fill_mi;
    private JMenuItem                   _choose_node_size_mi;
    private JCheckBoxMenuItem           _show_confidence_stddev_cbmi;
    private final LinkedList<TextFrame> _textframes      = new LinkedList<TextFrame>();
    private JMenu                       _analysis_menu;
    private JMenuItem                   _gsdi_item;
    private JMenuItem                   _gsdir_item;
    private Phylogeny                   _species_tree;
    private JCheckBoxMenuItem           _right_line_up_domains_cbmi;
    private JCheckBoxMenuItem           _line_up_renderable_data_cbmi;

    @Override
    public void actionPerformed( final ActionEvent e ) {
        final Object o = e.getSource();
        if ( o == _midpoint_root_item ) {
            getMainPanel().getCurrentTreePanel().midpointRoot();
        }
        else if ( o == _gsdi_item ) {
            if ( isSubtreeDisplayed() ) {
                return;
            }
            executeGSDI();
        }
        else if ( o == _gsdir_item ) {
            if ( isSubtreeDisplayed() ) {
                return;
            }
            executeGSDIR();
        }
        else if ( o == _taxcolor_item ) {
            getMainPanel().getCurrentTreePanel().taxColor();
        }
        else if ( o == _confcolor_item ) {
            getMainPanel().getCurrentTreePanel().confColor();
        }
        else if ( o == _collapse_species_specific_subtrees ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().collapseSpeciesSpecificSubtrees();
            }
        }
        else if ( o == _remove_branch_color_item ) {
            removeBranchColors();
        }
        else if ( o == _remove_visual_styles_item ) {
            removeVisualStyles();
        }
        else if ( o == _switch_colors_mi ) {
            switchColors();
        }
        else if ( o == _display_basic_information_item ) {
            displayBasicInformation();
        }
        else if ( o == _view_as_NH_item ) {
            viewAsNH();
        }
        else if ( o == _view_as_XML_item ) {
            viewAsXML();
        }
        else if ( o == _view_as_nexus_item ) {
            viewAsNexus();
        }
        else if ( o == _super_tiny_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setSuperTinyFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _tiny_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setTinyFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _small_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setSmallFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _medium_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setMediumFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _large_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setLargeFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _choose_font_mi ) {
            chooseFont();
        }
        else if ( o == _choose_minimal_confidence_mi ) {
            chooseMinimalConfidence();
        }
        else if ( o == _choose_node_size_mi ) {
            MainFrame.chooseNodeSize( getOptions(), this );
        }
        else if ( o == _overview_placment_mi ) {
            MainFrame.cycleOverview( getOptions(), getCurrentTreePanel() );
        }
        else if ( o == _cycle_node_fill_mi ) {
            MainFrame.cycleNodeFill( getOptions(), getCurrentTreePanel() );
        }
        else if ( o == _cycle_node_shape_mi ) {
            MainFrame.cycleNodeShape( getOptions(), getCurrentTreePanel() );
        }
        else if ( o == _non_lined_up_cladograms_rbmi ) {
            updateOptions( getOptions() );
            _mainpanel.getControlPanel().showWhole();
        }
        else if ( o == _uniform_cladograms_rbmi ) {
            updateOptions( getOptions() );
            _mainpanel.getControlPanel().showWhole();
        }
        else if ( o == _ext_node_dependent_cladogram_rbmi ) {
            updateOptions( getOptions() );
            _mainpanel.getControlPanel().showWhole();
        }
        else if ( o == _search_case_senstive_cbmi ) {
            updateOptions( getOptions() );
            getMainPanel().getControlPanel().search0();
            getMainPanel().getControlPanel().search1();
        }
        else if ( o == _search_whole_words_only_cbmi ) {
            updateOptions( getOptions() );
            getMainPanel().getControlPanel().search0();
            getMainPanel().getControlPanel().search1();
        }
        else if ( o == _inverse_search_result_cbmi ) {
            updateOptions( getOptions() );
            getMainPanel().getControlPanel().search0();
            getMainPanel().getControlPanel().search1();
        }
        else if ( o == _show_scale_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _show_branch_length_values_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _show_confidence_stddev_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _label_direction_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _abbreviate_scientific_names ) {
            updateOptions( getOptions() );
        }
        else if ( o == _show_overview_cbmi ) {
            updateOptions( getOptions() );
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().updateOvSizes();
            }
        }
        else if ( ( o == _rectangular_type_cbmi ) || ( o == _triangular_type_cbmi ) || ( o == _curved_type_cbmi )
                || ( o == _convex_type_cbmi ) || ( o == _rounded_type_cbmi ) || ( o == _euro_type_cbmi )
                || ( o == _unrooted_type_cbmi ) || ( o == _circular_type_cbmi ) ) {
            typeChanged( o );
        }
        else if ( o == _screen_antialias_cbmi ) {
            updateOptions( getOptions() );
            setupScreenTextAntialias( getMainPanel().getTreePanels(), isScreenAntialias() );
        }
        else if ( o == _background_gradient_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _show_domain_labels ) {
            updateOptions( getOptions() );
        }
        else if ( o == _color_labels_same_as_parent_branch ) {
            updateOptions( getOptions() );
        }
        else if ( o == _show_default_node_shapes_internal_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _show_default_node_shapes_external_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _about_item ) {
            MainFrame.about();
        }
        else if ( o == _help_item ) {
            try {
                AptxUtil.openWebsite( Constants.APTX_DOC_SITE, true, this );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _website_item ) {
            try {
                AptxUtil.openWebsite( Constants.APTX_WEB_SITE, true, this );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _phyloxml_website_item ) {
            try {
                AptxUtil.openWebsite( Constants.PHYLOXML_WEB_SITE, true, this );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _aptx_ref_item ) {
            try {
                AptxUtil.openWebsite( Constants.APTX_REFERENCE_URL, true, this );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _phyloxml_ref_item ) {
            try {
                AptxUtil.openWebsite( Constants.PHYLOXML_REFERENCE_URL, true, this );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _color_by_taxonomic_group_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _line_up_renderable_data_cbmi ) {
            if ( !_line_up_renderable_data_cbmi.isSelected() ) {
                _right_line_up_domains_cbmi.setSelected( false );
            }
            updateOptions( getOptions() );
        }
        else if ( o == _right_line_up_domains_cbmi ) {
            if ( _right_line_up_domains_cbmi.isSelected() ) {
                _line_up_renderable_data_cbmi.setSelected( true );
            }
            updateOptions( getOptions() );
        }
        repaint();
    }

    @Override
    public void destroy() {
        AptxUtil.printAppletMessage( NAME, "going to be destroyed " );
        removeAllTextFrames();
        if ( getMainPanel() != null ) {
            getMainPanel().terminate();
        }
    }

    /**
     * This method returns the current external node data which
     * has been selected by the user by clicking the "Return ..."
     * menu item. This method is expected to be called from Javascript or
     * something like it.
     * 
     * @return current external node data as String
     */
    public String getCurrentExternalNodesDataBuffer() {
        return getCurrentTreePanel().getCurrentExternalNodesDataBufferAsString();
    }

    public int getCurrentExternalNodesDataBufferChangeCounter() {
        return getCurrentTreePanel().getCurrentExternalNodesDataBufferChangeCounter();
    }

    public int getCurrentExternalNodesDataBufferLength() {
        return getCurrentTreePanel().getCurrentExternalNodesDataBufferAsString().length();
    }

    /**
     * This method returns the current phylogeny as a string in the chosen format
     * 
     * @param format must be NH, NHX, NEXUS or PHYLOXML
     * @return the phylogeny string
     * @author Herve Menager
     */
    public String getCurrentPhylogeny( final String format ) {
        removeAllTextFrames();
        if ( ( getMainPanel().getCurrentPhylogeny() == null ) || getMainPanel().getCurrentPhylogeny().isEmpty()
                || ( getMainPanel().getCurrentPhylogeny().getNumberOfExternalNodes() > 10000 ) ) {
            return new String();
        }
        switch ( ForesterConstants.PhylogeneticTreeFormats.valueOf( format ) ) {
            case NH:
                return getMainPanel().getCurrentPhylogeny().toNewHampshire();
            case NHX:
                return getMainPanel().getCurrentPhylogeny().toNewHampshireX();
            case NEXUS:
                return getMainPanel().getCurrentPhylogeny().toNexus();
            case PHYLOXML:
                return getMainPanel().getCurrentPhylogeny().toPhyloXML( -1 );
            default:
                break;
        }
        return new String();
    }

    /**
     * This method returns a view of the current phylogeny in a chosen 
     * graphics format, base64-encoded in a string so that in can be used
     * from javascript.
     * 
     * @param format must be GraphicsExportType (gif, jpg, pdf, png, tif, bmp)
     * @return the phylogeny string
     * @author Herve Menager
     */
    public String getCurrentPhylogenyGraphicsAsBase64EncodedString( final String format ) {
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try {
            AptxUtil.writePhylogenyToGraphicsByteArrayOutputStream( baos,
                                                                    _mainpanel.getWidth(),
                                                                    _mainpanel.getHeight(),
                                                                    getCurrentTreePanel(),
                                                                    getCurrentTreePanel().getControlPanel(),
                                                                    GraphicsExportType.valueOf( format ),
                                                                    getOptions() );
        }
        catch ( final IOException ioe ) {
            ForesterUtil.printErrorMessage( NAME, ioe.toString() );
            ioe.printStackTrace();
            JOptionPane.showMessageDialog( this,
                                           NAME + ": Failed to generate graphics: " + "\nException: " + ioe,
                                           "Failed to generate graphics",
                                           JOptionPane.ERROR_MESSAGE );
            return null;
        }
        final byte[] bytes = baos.toByteArray();
        final String dataImg = Base64.encodeBase64String( bytes );
        return dataImg;
    }

    public Options getOptions() {
        return _options;
    }

    @Override
    public void init() {
        final String config_filename = getParameter( Constants.APPLET_PARAM_NAME_FOR_CONFIG_FILE_URL );
        AptxUtil.printAppletMessage( NAME, "URL for configuration file is: " + config_filename );
        final Configuration configuration = new Configuration( config_filename, true, true, true );
        setConfiguration( configuration );
        setOptions( Options.createInstance( configuration ) );
        setupUI();
        final String tree_url_str = getParameter( Constants.APPLET_PARAM_NAME_FOR_URL_OF_TREE_TO_LOAD );
        if ( ForesterUtil.isEmpty( tree_url_str ) ) {
            ForesterUtil.printErrorMessage( NAME, "could not get tree URL from "
                    + Constants.APPLET_PARAM_NAME_FOR_URL_OF_TREE_TO_LOAD );
            JOptionPane.showMessageDialog( this, NAME + ": could not get tree URL from "
                    + Constants.APPLET_PARAM_NAME_FOR_URL_OF_TREE_TO_LOAD, "Failed get URL", JOptionPane.ERROR_MESSAGE );
            return;
        }
        AptxUtil.printAppletMessage( NAME, "URL for phylogenies is " + tree_url_str );
        // Get URL to tree file
        URL phys_url = null;
        try {
            phys_url = new URL( tree_url_str );
        }
        catch ( final Exception e ) {
            ForesterUtil.printErrorMessage( NAME, "error: " + e );
            e.printStackTrace();
            JOptionPane.showMessageDialog( this, NAME + ": Could not create URL from: \"" + tree_url_str
                    + "\"\nException: " + e, "Failed to create URL", JOptionPane.ERROR_MESSAGE );
        }
        if ( phys_url == null ) {
            ForesterUtil.printErrorMessage( NAME, "failed to get tree URL from "
                    + Constants.APPLET_PARAM_NAME_FOR_URL_OF_TREE_TO_LOAD );
            JOptionPane.showMessageDialog( this,
                                           NAME + ": Could not create URL from: \"" + tree_url_str,
                                           "Failed to create URL",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        // Load the tree from URL
        Phylogeny[] phys = null;
        try {
            phys = AptxUtil.readPhylogeniesFromUrl( phys_url,
                                                    getConfiguration().isValidatePhyloXmlAgainstSchema(),
                                                    getConfiguration().isReplaceUnderscoresInNhParsing(),
                                                    getConfiguration().isInternalNumberAreConfidenceForNhParsing(),
                                                    getConfiguration().getTaxonomyExtraction(),
                                                    getConfiguration().isMidpointReroot() );
        }
        catch ( final Exception e ) {
            ForesterUtil.printErrorMessage( NAME, e.toString() );
            e.printStackTrace();
            JOptionPane.showMessageDialog( this,
                                           NAME + ": Failed to read phylogenies: " + "\nException: " + e,
                                           "Failed to read phylogenies",
                                           JOptionPane.ERROR_MESSAGE );
        }
        if ( phys == null ) {
            ForesterUtil.printErrorMessage( NAME, "phylogenies from [" + phys_url + "] are null" );
            JOptionPane.showMessageDialog( this,
                                           NAME + ": phylogenies from [" + phys_url + "] are null",
                                           "Failed to read phylogenies",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        else if ( phys.length < 1 ) {
            ForesterUtil.printErrorMessage( NAME, "phylogenies from [" + phys_url + "] are empty" );
            JOptionPane.showMessageDialog( this,
                                           NAME + ": phylogenies from [" + phys_url + "] are empty",
                                           "Failed to read phylogenies",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        else {
            AptxUtil.printAppletMessage( NAME, "loaded " + phys.length + " phylogenies from: " + phys_url );
        }
        //
        final String species_tree_url_str = getParameter( Constants.APPLET_PARAM_NAME_FOR_URL_OF_SPECIES_TREE_TO_LOAD );
        if ( !ForesterUtil.isEmpty( species_tree_url_str ) ) {
            AptxUtil.printAppletMessage( NAME, "URL of species tree to load: \"" + species_tree_url_str + "\"" );
            Phylogeny[] species_trees = null;
            try {
                final URL species_tree_url = new URL( species_tree_url_str );
                species_trees = AptxUtil.readPhylogeniesFromUrl( species_tree_url,
                                                                 configuration.isValidatePhyloXmlAgainstSchema(),
                                                                 configuration.isReplaceUnderscoresInNhParsing(),
                                                                 false,
                                                                 TAXONOMY_EXTRACTION.NO,
                                                                 false );
            }
            catch ( final IOException e ) {
                ForesterUtil.printErrorMessage( NAME, "could not read species tree from  [" + species_tree_url_str
                        + "]" );
                JOptionPane.showMessageDialog( this, NAME + ": could not read species tree from  ["
                        + species_tree_url_str + "]", "Failed to read species tree", JOptionPane.ERROR_MESSAGE );
            }
            if ( ( species_trees != null ) && ( species_trees.length > 0 ) ) {
                AptxUtil.printAppletMessage( NAME, "successfully read species tree" );
                if ( species_trees[ 0 ].isEmpty() ) {
                    ForesterUtil.printErrorMessage( NAME, "species tree is empty" );
                }
                else if ( !species_trees[ 0 ].isRooted() ) {
                    ForesterUtil.printErrorMessage( NAME, "species tree is not rooted" );
                }
                else {
                    setSpeciesTree( species_trees[ 0 ] );
                    AptxUtil.printAppletMessage( NAME, "species tree OK" );
                }
            }
        }
        //
        setVisible( false );
        setMainPanel( new MainPanelApplets( getConfiguration(), this ) );
        _jmenubar = new JMenuBar();
        if ( !getConfiguration().isHideControlPanelAndMenubar() ) {
            if ( !getConfiguration().isUseNativeUI() ) {
                _jmenubar.setBackground( getConfiguration().getGuiMenuBackgroundColor() );
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
        }
        final Container contentpane = getContentPane();
        contentpane.setLayout( new BorderLayout() );
        contentpane.add( getMainPanel(), BorderLayout.CENTER );
        addComponentListener( new ComponentAdapter() {

            @Override
            public void componentResized( final ComponentEvent e ) {
                if ( getMainPanel().getCurrentTreePanel() != null ) {
                    getMainPanel().getCurrentTreePanel().calcParametersForPainting( getMainPanel()
                                                                                            .getCurrentTreePanel()
                                                                                            .getWidth(),
                                                                                    getMainPanel()
                                                                                            .getCurrentTreePanel()
                                                                                            .getHeight(),
                                                                                    getOptions()
                                                                                            .isAllowFontSizeChange() );
                }
            }
        } );
        if ( getConfiguration().isUseTabbedDisplay() ) {
            AptxUtil.printAppletMessage( NAME, "using tabbed display" );
            AptxUtil.addPhylogeniesToTabs( phys,
                                           new File( phys_url.getFile() ).getName(),
                                           phys_url.toString(),
                                           getConfiguration(),
                                           getMainPanel() );
        }
        else {
            AptxUtil.printAppletMessage( NAME, "not using tabbed display" );
            if ( getSpeciesTree() != null ) {
                AptxUtil.printAppletMessage( NAME,
                                             "Warning: gsdi (gene duplication inference) only available tabbed display" );
            }
            AptxUtil.addPhylogenyToPanel( phys, getConfiguration(), getMainPanel() );
        }
        validate();
        setName( NAME );
        getMainPanel().getControlPanel().showWholeAll();
        getMainPanel().getControlPanel().showWhole();
        /* GUILHEM_BEG */
        getCurrentTreePanel().getControlPanel().getSequenceRelationTypeBox().removeAllItems();
        for( final SequenceRelation.SEQUENCE_RELATION_TYPE type : getMainPanel().getCurrentPhylogeny()
                .getRelevantSequenceRelationTypes() ) {
            getCurrentTreePanel().getControlPanel().getSequenceRelationTypeBox().addItem( type );
        }
        final String default_relation = getParameter( Constants.APPLET_PARAM_NAME_FOR_DEFAULT_SEQUENCE_RELATION_TYPE );
        if ( default_relation != null ) {
            getCurrentTreePanel().getControlPanel().getSequenceRelationTypeBox().setSelectedItem( default_relation );
        }
        final String default_sequence = getParameter( Constants.APPLET_PARAM_NAME_FOR_DEFAULT_QUERY_SEQUENCE );
        if ( default_sequence != null ) {
            getCurrentTreePanel().getControlPanel().getSequenceRelationBox().setSelectedItem( default_sequence );
        }
        /* GUILHEM_END */
        System.gc();
        AptxUtil.printAppletMessage( NAME, "successfully initialized" );
        setVisible( true );
    }

    public void showTextFrame( final String s, final String title ) {
        checkTextFrames();
        _textframes.addLast( TextFrame.instantiate( s, title, _textframes ) );
    }

    @Override
    public void start() {
        if ( getMainPanel() != null ) {
            getMainPanel().validate();
        }
        requestFocus();
        requestFocusInWindow();
        requestFocus();
        AptxUtil.printAppletMessage( NAME, "started" );
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

    void buildFontSizeMenu() {
        _font_size_menu = MainFrame.createMenu( MainFrame.FONT_SIZE_MENU_LABEL, getConfiguration() );
        _font_size_menu.add( _super_tiny_fonts_mi = new JMenuItem( "Super tiny fonts" ) );
        _font_size_menu.add( _tiny_fonts_mi = new JMenuItem( "Tiny fonts" ) );
        _font_size_menu.add( _small_fonts_mi = new JMenuItem( "Small fonts" ) );
        _font_size_menu.add( _medium_fonts_mi = new JMenuItem( "Medium fonts" ) );
        _font_size_menu.add( _large_fonts_mi = new JMenuItem( "Large fonts" ) );
        customizeJMenuItem( _super_tiny_fonts_mi );
        customizeJMenuItem( _tiny_fonts_mi );
        customizeJMenuItem( _small_fonts_mi );
        customizeJMenuItem( _medium_fonts_mi );
        customizeJMenuItem( _large_fonts_mi );
        _jmenubar.add( _font_size_menu );
    }

    void buildHelpMenu() {
        _help_jmenu = MainFrame.createMenu( "Help", getConfiguration() );
        _help_jmenu.add( _help_item = new JMenuItem( "Documentation" ) );
        _help_jmenu.addSeparator();
        _help_jmenu.add( _website_item = new JMenuItem( "Archaeopteryx Home" ) );
        _aptx_ref_item = new JMenuItem( "Archaeopteryx Reference" );
        _help_jmenu.add( _phyloxml_website_item = new JMenuItem( "phyloXML Home" ) );
        _help_jmenu.add( _phyloxml_ref_item = new JMenuItem( "phyloXML Reference" ) );
        _help_jmenu.addSeparator();
        _help_jmenu.add( _about_item = new JMenuItem( "About" ) );
        customizeJMenuItem( _help_item );
        customizeJMenuItem( _website_item );
        customizeJMenuItem( _phyloxml_website_item );
        customizeJMenuItem( _aptx_ref_item );
        customizeJMenuItem( _phyloxml_ref_item );
        customizeJMenuItem( _about_item );
        _phyloxml_ref_item.setToolTipText( MainFrame.PHYLOXML_REF_TOOL_TIP );
        _aptx_ref_item.setToolTipText( MainFrame.APTX_REF_TOOL_TIP );
        _jmenubar.add( _help_jmenu );
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
                MainFrame.setTextForFontChooserMenuItem( _choose_font_mi, MainFrame
                        .createCurrentFontDesc( getMainPanel().getTreeFontSet() ) );
                MainFrame.updateOptionsMenuDependingOnPhylogenyType( getMainPanel(),
                                                                     _show_scale_cbmi,
                                                                     _show_branch_length_values_cbmi,
                                                                     _non_lined_up_cladograms_rbmi,
                                                                     _uniform_cladograms_rbmi,
                                                                     _ext_node_dependent_cladogram_rbmi,
                                                                     _label_direction_cbmi );
                MainFrame.setCycleNodeFillMenuItem( _cycle_node_fill_mi, getOptions() );
                MainFrame.setCycleNodeShapeMenuItem( _cycle_node_shape_mi, getOptions() );
                MainFrame.setTextNodeSizeMenuItem( _choose_node_size_mi, getOptions() );
            }
        } );
        _options_jmenu.add( MainFrame.customizeMenuItemAsLabel( new JMenuItem( MainFrame.DISPLAY_SUBHEADER ),
                                                                getConfiguration() ) );
        _options_jmenu
                .add( _ext_node_dependent_cladogram_rbmi = new JRadioButtonMenuItem( MainFrame.NONUNIFORM_CLADOGRAMS_LABEL ) );
        _options_jmenu.add( _uniform_cladograms_rbmi = new JRadioButtonMenuItem( MainFrame.UNIFORM_CLADOGRAMS_LABEL ) );
        _options_jmenu
                .add( _non_lined_up_cladograms_rbmi = new JRadioButtonMenuItem( MainFrame.NON_LINED_UP_CLADOGRAMS_LABEL ) );
        _radio_group_1 = new ButtonGroup();
        _radio_group_1.add( _ext_node_dependent_cladogram_rbmi );
        _radio_group_1.add( _uniform_cladograms_rbmi );
        _radio_group_1.add( _non_lined_up_cladograms_rbmi );
        /////
        _options_jmenu.add( _show_overview_cbmi = new JCheckBoxMenuItem( MainFrame.SHOW_OVERVIEW_LABEL ) );
        _options_jmenu.add( _show_scale_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_SCALE_LABEL ) );
        _options_jmenu
                .add( _show_branch_length_values_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_BRANCH_LENGTH_VALUES_LABEL ) );
        _options_jmenu
                .add( _show_default_node_shapes_internal_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_NODE_BOXES_LABEL_INT ) );
        _options_jmenu
                .add( _show_default_node_shapes_external_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_NODE_BOXES_LABEL_EXT ) );
        if ( getConfiguration().doDisplayOption( Configuration.show_domain_architectures ) ) {
            _options_jmenu.add( _show_domain_labels = new JCheckBoxMenuItem( MainFrame.SHOW_DOMAIN_LABELS_LABEL ) );
            _options_jmenu.add( _right_line_up_domains_cbmi = new JCheckBoxMenuItem( MainFrame.RIGHT_LINE_UP_DOMAINS ) );
        }
        _options_jmenu.add( _line_up_renderable_data_cbmi = new JCheckBoxMenuItem( MainFrame.LINE_UP_RENDERABLE_DATA ) );
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
        /////
        _options_jmenu.addSeparator();
        _options_jmenu.add( MainFrame.customizeMenuItemAsLabel( new JMenuItem( MainFrame.SEARCH_SUBHEADER ),
                                                                getConfiguration() ) );
        _options_jmenu
                .add( _search_case_senstive_cbmi = new JCheckBoxMenuItem( MainFrame.SEARCH_CASE_SENSITIVE_LABEL ) );
        _options_jmenu.add( _search_whole_words_only_cbmi = new JCheckBoxMenuItem( MainFrame.SEARCH_TERMS_ONLY_LABEL ) );
        _options_jmenu
                .add( _inverse_search_result_cbmi = new JCheckBoxMenuItem( MainFrame.INVERSE_SEARCH_RESULT_LABEL ) );
        customizeJMenuItem( _choose_font_mi );
        customizeJMenuItem( _choose_minimal_confidence_mi );
        customizeJMenuItem( _switch_colors_mi );
        customizeJMenuItem( _overview_placment_mi );
        customizeCheckBoxMenuItem( _color_by_taxonomic_group_cbmi, getOptions().isColorByTaxonomicGroup() );
        customizeCheckBoxMenuItem( _label_direction_cbmi,
                                   getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL );
        customizeCheckBoxMenuItem( _screen_antialias_cbmi, getOptions().isAntialiasScreen() );
        customizeCheckBoxMenuItem( _background_gradient_cbmi, getOptions().isBackgroundColorGradient() );
        customizeCheckBoxMenuItem( _show_domain_labels, getOptions().isShowDomainLabels() );
        customizeCheckBoxMenuItem( _show_annotation_ref_source, getOptions().isShowAnnotationRefSource() );
        customizeCheckBoxMenuItem( _abbreviate_scientific_names, getOptions().isAbbreviateScientificTaxonNames() );
        customizeCheckBoxMenuItem( _show_default_node_shapes_external_cbmi, getOptions()
                .isShowDefaultNodeShapesExternal() );
        customizeCheckBoxMenuItem( _show_default_node_shapes_internal_cbmi, getOptions()
                .isShowDefaultNodeShapesInternal() );
        customizeJMenuItem( _cycle_node_shape_mi );
        customizeJMenuItem( _cycle_node_fill_mi );
        customizeJMenuItem( _choose_node_size_mi );
        customizeCheckBoxMenuItem( _color_labels_same_as_parent_branch, getOptions().isColorLabelsSameAsParentBranch() );
        customizeCheckBoxMenuItem( _search_case_senstive_cbmi, getOptions().isSearchCaseSensitive() );
        customizeCheckBoxMenuItem( _show_scale_cbmi, getOptions().isShowScale() );
        customizeRadioButtonMenuItem( _non_lined_up_cladograms_rbmi,
                                      getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP );
        customizeRadioButtonMenuItem( _uniform_cladograms_rbmi,
                                      getOptions().getCladogramType() == CLADOGRAM_TYPE.TOTAL_NODE_SUM_DEP );
        customizeRadioButtonMenuItem( _ext_node_dependent_cladogram_rbmi,
                                      getOptions().getCladogramType() == CLADOGRAM_TYPE.EXT_NODE_SUM_DEP );
        customizeCheckBoxMenuItem( _show_branch_length_values_cbmi, getOptions().isShowBranchLengthValues() );
        customizeCheckBoxMenuItem( _show_overview_cbmi, getOptions().isShowOverview() );
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

    void buildTypeMenu() {
        _type_menu = MainFrame.createMenu( MainFrame.TYPE_MENU_HEADER, getConfiguration() );
        _type_menu.add( _rectangular_type_cbmi = new JCheckBoxMenuItem( MainFrame.RECTANGULAR_TYPE_CBMI_LABEL ) );
        _type_menu.add( _euro_type_cbmi = new JCheckBoxMenuItem( MainFrame.EURO_TYPE_CBMI_LABEL ) );
        _type_menu.add( _rounded_type_cbmi = new JCheckBoxMenuItem( MainFrame.ROUNDED_TYPE_CBMI_LABEL ) );
        _type_menu.add( _curved_type_cbmi = new JCheckBoxMenuItem( MainFrame.CURVED_TYPE_CBMI_LABEL ) );
        _type_menu.add( _triangular_type_cbmi = new JCheckBoxMenuItem( MainFrame.TRIANGULAR_TYPE_CBMI_LABEL ) );
        _type_menu.add( _convex_type_cbmi = new JCheckBoxMenuItem( MainFrame.CONVEX_TYPE_CBMI_LABEL ) );
        _type_menu.add( _unrooted_type_cbmi = new JCheckBoxMenuItem( MainFrame.UNROOTED_TYPE_CBMI_LABEL ) );
        _type_menu.add( _circular_type_cbmi = new JCheckBoxMenuItem( MainFrame.CIRCULAR_TYPE_CBMI_LABEL ) );
        customizeCheckBoxMenuItem( _rectangular_type_cbmi, false );
        customizeCheckBoxMenuItem( _triangular_type_cbmi, false );
        customizeCheckBoxMenuItem( _euro_type_cbmi, false );
        customizeCheckBoxMenuItem( _rounded_type_cbmi, false );
        customizeCheckBoxMenuItem( _curved_type_cbmi, false );
        customizeCheckBoxMenuItem( _convex_type_cbmi, false );
        customizeCheckBoxMenuItem( _unrooted_type_cbmi, false );
        customizeCheckBoxMenuItem( _circular_type_cbmi, false );
        _unrooted_type_cbmi.setToolTipText( MainFrame.USE_MOUSEWHEEL_SHIFT_TO_ROTATE );
        _circular_type_cbmi.setToolTipText( MainFrame.USE_MOUSEWHEEL_SHIFT_TO_ROTATE );
        initializeTypeMenu( getOptions() );
        _jmenubar.add( _type_menu );
    }

    void buildViewMenu() {
        _view_jmenu = MainFrame.createMenu( "View", getConfiguration() );
        _view_jmenu
                .add( _display_basic_information_item = new JMenuItem( MainFrame.SHOW_BASIC_TREE_INFORMATION_LABEL ) );
        _view_jmenu.addSeparator();
        _view_jmenu.add( _view_as_XML_item = new JMenuItem( "as phyloXML" ) );
        _view_jmenu.add( _view_as_NH_item = new JMenuItem( "as Newick" ) );
        _view_jmenu.add( _view_as_nexus_item = new JMenuItem( "as Nexus" ) );
        customizeJMenuItem( _display_basic_information_item );
        customizeJMenuItem( _view_as_NH_item );
        customizeJMenuItem( _view_as_XML_item );
        customizeJMenuItem( _view_as_nexus_item );
        _jmenubar.add( _view_jmenu );
    }

    void checkTextFrames() {
        if ( _textframes.size() > 5 ) {
            try {
                if ( _textframes.getFirst() != null ) {
                    _textframes.getFirst().removeMe();
                }
                else {
                    _textframes.removeFirst();
                }
            }
            catch ( final NoSuchElementException e ) {
                // Ignore.
            }
        }
    }

    void clearCurrentExternalNodesDataBuffer() {
        getCurrentTreePanel().clearCurrentExternalNodesDataBuffer();
    }

    void customizeCheckBoxMenuItem( final JCheckBoxMenuItem item, final boolean is_selected ) {
        if ( item != null ) {
            item.setFont( MainFrame.menu_font );
            if ( !getConfiguration().isUseNativeUI() ) {
                item.setBackground( getConfiguration().getGuiMenuBackgroundColor() );
                item.setForeground( getConfiguration().getGuiMenuTextColor() );
            }
            item.setSelected( is_selected );
            item.addActionListener( this );
        }
    }

    void customizeJMenuItem( final JMenuItem jmi ) {
        jmi.setFont( MainFrame.menu_font );
        if ( !getConfiguration().isUseNativeUI() ) {
            jmi.setBackground( getConfiguration().getGuiMenuBackgroundColor() );
            jmi.setForeground( getConfiguration().getGuiMenuTextColor() );
        }
        jmi.addActionListener( this );
    }

    void displayBasicInformation() {
        if ( ( getMainPanel().getCurrentPhylogeny() != null ) && !getMainPanel().getCurrentPhylogeny().isEmpty() ) {
            String title = "Basic Information";
            if ( !ForesterUtil.isEmpty( getMainPanel().getCurrentPhylogeny().getName() ) ) {
                title = title + " for \"" + _mainpanel.getCurrentPhylogeny().getName() + "\"";
            }
            showTextFrame( AptxUtil.createBasicInformation( getMainPanel().getCurrentPhylogeny(), null ), title );
        }
    }

    void executeGSDI() {
        if ( !isOKforSDI( false, true ) ) {
            return;
        }
        if ( !_mainpanel.getCurrentPhylogeny().isRooted() ) {
            JOptionPane.showMessageDialog( this,
                                           "Gene tree is not rooted.",
                                           "Cannot execute GSDI",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        final Phylogeny gene_tree = _mainpanel.getCurrentPhylogeny().copy();
        gene_tree.setAllNodesToNotCollapse();
        gene_tree.recalculateNumberOfExternalDescendants( false );
        GSDI gsdi = null;
        final Phylogeny species_tree = _species_tree.copy();
        try {
            gsdi = new GSDI( gene_tree, species_tree, false, true, true, true );
        }
        catch ( final SDIException e ) {
            JOptionPane.showMessageDialog( this,
                                           e.getLocalizedMessage(),
                                           "Error during GSDI",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final Exception e ) {
            AptxUtil.unexpectedException( e );
            return;
        }
        gene_tree.setRerootable( false );
        gene_tree.clearHashIdToNodeMap();
        gene_tree.recalculateNumberOfExternalDescendants( true );
        _mainpanel.addPhylogenyInNewTab( gene_tree, getConfiguration(), "gene tree", null );
        getMainPanel().getControlPanel().setShowEvents( true );
        showWhole();
        final int selected = _mainpanel.getTabbedPane().getSelectedIndex();
        _mainpanel.addPhylogenyInNewTab( species_tree, getConfiguration(), "species tree", null );
        showWhole();
        _mainpanel.getTabbedPane().setSelectedIndex( selected );
        showWhole();
        _mainpanel.getCurrentTreePanel().setEdited( true );
        final int poly = PhylogenyMethods.countNumberOfPolytomies( species_tree );
        if ( gsdi.getStrippedExternalGeneTreeNodes().size() > 0 ) {
            JOptionPane.showMessageDialog( this,
                                           "Duplications: " + gsdi.getDuplicationsSum() + "\n"
                                                   + "Potential duplications: "
                                                   + gsdi.getSpeciationOrDuplicationEventsSum() + "\n"
                                                   + "Speciations: " + gsdi.getSpeciationsSum() + "\n"
                                                   + "Stripped gene tree nodes: "
                                                   + gsdi.getStrippedExternalGeneTreeNodes().size() + "\n"
                                                   + "Taxonomy linkage based on: " + gsdi.getTaxCompBase() + "\n"
                                                   + "Number of polytomies in species tree used: " + poly + "\n",
                                           "GSDI successfully completed",
                                           JOptionPane.WARNING_MESSAGE );
        }
        else {
            JOptionPane.showMessageDialog( this,
                                           "Duplications: " + gsdi.getDuplicationsSum() + "\n"
                                                   + "Potential duplications: "
                                                   + gsdi.getSpeciationOrDuplicationEventsSum() + "\n"
                                                   + "Speciations: " + gsdi.getSpeciationsSum() + "\n"
                                                   + "Stripped gene tree nodes: "
                                                   + gsdi.getStrippedExternalGeneTreeNodes().size() + "\n"
                                                   + "Taxonomy linkage based on: " + gsdi.getTaxCompBase() + "\n"
                                                   + "Number of polytomies in species tree used: " + poly + "\n",
                                           "GSDI successfully completed",
                                           JOptionPane.INFORMATION_MESSAGE );
        }
    }

    void executeGSDIR() {
        if ( !isOKforSDI( false, false ) ) {
            return;
        }
        final int p = PhylogenyMethods.countNumberOfPolytomies( _mainpanel.getCurrentPhylogeny() );
        if ( ( p > 0 )
                && !( ( p == 1 ) && ( _mainpanel.getCurrentPhylogeny().getRoot().getNumberOfDescendants() == 3 ) ) ) {
            JOptionPane.showMessageDialog( this,
                                           "Gene tree is not completely binary",
                                           "Cannot execute GSDI",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        final Phylogeny gene_tree = _mainpanel.getCurrentPhylogeny().copy();
        gene_tree.setAllNodesToNotCollapse();
        gene_tree.recalculateNumberOfExternalDescendants( false );
        GSDIR gsdir = null;
        final Phylogeny species_tree = _species_tree.copy();
        try {
            gsdir = new GSDIR( gene_tree, species_tree, true, true, true );
        }
        catch ( final SDIException e ) {
            JOptionPane.showMessageDialog( this,
                                           e.getLocalizedMessage(),
                                           "Error during GSDIR",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final Exception e ) {
            AptxUtil.unexpectedException( e );
            return;
        }
        final Phylogeny result_gene_tree = gsdir.getMinDuplicationsSumGeneTree();
        result_gene_tree.setRerootable( false );
        result_gene_tree.clearHashIdToNodeMap();
        result_gene_tree.recalculateNumberOfExternalDescendants( true );
        PhylogenyMethods.orderAppearance( result_gene_tree.getRoot(), true, true, DESCENDANT_SORT_PRIORITY.NODE_NAME );
        _mainpanel.addPhylogenyInNewTab( result_gene_tree, getConfiguration(), "gene tree", null );
        getMainPanel().getControlPanel().setShowEvents( true );
        showWhole();
        final int selected = _mainpanel.getTabbedPane().getSelectedIndex();
        _mainpanel.addPhylogenyInNewTab( species_tree, getConfiguration(), "species tree", null );
        showWhole();
        _mainpanel.getTabbedPane().setSelectedIndex( selected );
        showWhole();
        _mainpanel.getCurrentTreePanel().setEdited( true );
        final int poly = PhylogenyMethods.countNumberOfPolytomies( species_tree );
        if ( gsdir.getStrippedExternalGeneTreeNodes().size() > 0 ) {
            JOptionPane.showMessageDialog( this,
                                           "Minimal duplications: " + gsdir.getMinDuplicationsSum() + "\n"
                                                   + "Speciations: " + gsdir.getSpeciationsSum() + "\n"
                                                   + "Stripped gene tree nodes: "
                                                   + gsdir.getStrippedExternalGeneTreeNodes().size() + "\n"
                                                   + "Taxonomy linkage based on: " + gsdir.getTaxCompBase() + "\n"
                                                   + "Number of polytomies in species tree used: " + poly + "\n",
                                           "GSDIR successfully completed",
                                           JOptionPane.WARNING_MESSAGE );
        }
        else {
            JOptionPane.showMessageDialog( this,
                                           "Minimal duplications: " + gsdir.getMinDuplicationsSum() + "\n"
                                                   + "Speciations: " + gsdir.getSpeciationsSum() + "\n"
                                                   + "Stripped gene tree nodes: "
                                                   + gsdir.getStrippedExternalGeneTreeNodes().size() + "\n"
                                                   + "Taxonomy linkage based on: " + gsdir.getTaxCompBase() + "\n"
                                                   + "Number of polytomies in species tree used: " + poly + "\n",
                                           "GSDIR successfully completed",
                                           JOptionPane.INFORMATION_MESSAGE );
        }
    }

    Configuration getConfiguration() {
        return _configuration;
    }

    TreePanel getCurrentTreePanel() {
        return getMainPanel().getCurrentTreePanel();
    }

    JCheckBoxMenuItem getlabelDirectionCbmi() {
        return _label_direction_cbmi;
    }

    Options getOtions() {
        return _options;
    }

    void initializeTypeMenu( final Options options ) {
        setTypeMenuToAllUnselected();
        try {
            switch ( options.getPhylogenyGraphicsType() ) {
                case CONVEX:
                    _convex_type_cbmi.setSelected( true );
                    break;
                case CURVED:
                    _curved_type_cbmi.setSelected( true );
                    break;
                case EURO_STYLE:
                    _euro_type_cbmi.setSelected( true );
                    break;
                case ROUNDED:
                    _rounded_type_cbmi.setSelected( true );
                    break;
                case TRIANGULAR:
                    _triangular_type_cbmi.setSelected( true );
                    break;
                case UNROOTED:
                    _unrooted_type_cbmi.setSelected( true );
                    break;
                case CIRCULAR:
                    _circular_type_cbmi.setSelected( true );
                    break;
                default:
                    _rectangular_type_cbmi.setSelected( true );
                    break;
            }
        }
        catch ( final NullPointerException np ) {
            // In all likelihood, this is caused by menu-less display.
        }
    }

    boolean isOKforSDI( final boolean species_tree_has_to_binary, final boolean gene_tree_has_to_binary ) {
        if ( ( _mainpanel.getCurrentPhylogeny() == null ) || _mainpanel.getCurrentPhylogeny().isEmpty() ) {
            return false;
        }
        else if ( ( _species_tree == null ) || _species_tree.isEmpty() ) {
            JOptionPane.showMessageDialog( this,
                                           "No species tree loaded",
                                           "Cannot execute GSDI",
                                           JOptionPane.ERROR_MESSAGE );
            return false;
        }
        else if ( species_tree_has_to_binary && !_species_tree.isCompletelyBinary() ) {
            JOptionPane.showMessageDialog( this,
                                           "Species tree is not completely binary",
                                           "Cannot execute GSDI",
                                           JOptionPane.ERROR_MESSAGE );
            return false;
        }
        else if ( gene_tree_has_to_binary && !_mainpanel.getCurrentPhylogeny().isCompletelyBinary() ) {
            JOptionPane.showMessageDialog( this,
                                           "Gene tree is not completely binary",
                                           "Cannot execute GSDI",
                                           JOptionPane.ERROR_MESSAGE );
            return false;
        }
        else {
            return true;
        }
    }

    boolean isSubtreeDisplayed() {
        if ( getCurrentTreePanel() != null ) {
            if ( getCurrentTreePanel().isCurrentTreeIsSubtree() ) {
                JOptionPane
                        .showMessageDialog( this,
                                            "This operation can only be performed on a complete tree, not on the currently displayed sub-tree only.",
                                            "Operation can not be exectuted on a sub-tree",
                                            JOptionPane.WARNING_MESSAGE );
                return true;
            }
        }
        return false;
    }

    void removeAllTextFrames() {
        for( final TextFrame tf : _textframes ) {
            if ( tf != null ) {
                tf.close();
            }
        }
        _textframes.clear();
    }

    void setConfiguration( final Configuration configuration ) {
        _configuration = configuration;
    }

    void setOptions( final Options options ) {
        _options = options;
    }

    void setSelectedTypeInTypeMenu( final PHYLOGENY_GRAPHICS_TYPE type ) {
        setTypeMenuToAllUnselected();
        try {
            switch ( type ) {
                case CIRCULAR:
                    _circular_type_cbmi.setSelected( true );
                    break;
                case CONVEX:
                    _convex_type_cbmi.setSelected( true );
                    break;
                case CURVED:
                    _curved_type_cbmi.setSelected( true );
                    break;
                case EURO_STYLE:
                    _euro_type_cbmi.setSelected( true );
                    break;
                case ROUNDED:
                    _rounded_type_cbmi.setSelected( true );
                    break;
                case RECTANGULAR:
                    _rectangular_type_cbmi.setSelected( true );
                    break;
                case TRIANGULAR:
                    _triangular_type_cbmi.setSelected( true );
                    break;
                case UNROOTED:
                    _unrooted_type_cbmi.setSelected( true );
                    break;
                default:
                    throw new IllegalArgumentException( "unknown type: " + type );
            }
        }
        catch ( final NullPointerException np ) {
            // In all likelihood, this is caused by menu-less display.
        }
    }

    void setTypeMenuToAllUnselected() {
        if ( _convex_type_cbmi != null ) {
            _convex_type_cbmi.setSelected( false );
        }
        if ( _curved_type_cbmi != null ) {
            _curved_type_cbmi.setSelected( false );
        }
        if ( _euro_type_cbmi != null ) {
            _euro_type_cbmi.setSelected( false );
        }
        if ( _rounded_type_cbmi != null ) {
            _rounded_type_cbmi.setSelected( false );
        }
        if ( _triangular_type_cbmi != null ) {
            _triangular_type_cbmi.setSelected( false );
        }
        if ( _rectangular_type_cbmi != null ) {
            _rectangular_type_cbmi.setSelected( false );
        }
        if ( _unrooted_type_cbmi != null ) {
            _unrooted_type_cbmi.setSelected( false );
        }
        if ( _circular_type_cbmi != null ) {
            _circular_type_cbmi.setSelected( false );
        }
    }

    void showWhole() {
        _mainpanel.getControlPanel().showWhole();
    }

    void switchColors() {
        final TreeColorSet colorset = getMainPanel().getCurrentTreePanel().getTreeColorSet();
        final ColorSchemeChooser csc = new ColorSchemeChooser( getMainPanel(), colorset );
        csc.setVisible( true );
        getMainPanel().setTreeColorSet( colorset );
    }

    void typeChanged( final Object o ) {
        updateTypeCheckboxes( getOptions(), o );
        updateOptions( getOptions() );
        if ( getCurrentTreePanel() != null ) {
            final PHYLOGENY_GRAPHICS_TYPE previous_type = getCurrentTreePanel().getPhylogenyGraphicsType();
            final PHYLOGENY_GRAPHICS_TYPE new_type = getOptions().getPhylogenyGraphicsType();
            if ( ( ( previous_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && ( new_type != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) )
                    || ( ( previous_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) && ( new_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) )
                    || ( ( previous_type != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && ( new_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) )
                    || ( ( previous_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) && ( new_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) ) {
                getCurrentTreePanel().getControlPanel().showWhole();
            }
            if ( getCurrentTreePanel().isPhyHasBranchLengths() && ( new_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                getCurrentTreePanel().getControlPanel().setDrawPhylogramEnabled( true );
            }
            else {
                getCurrentTreePanel().getControlPanel().setDrawPhylogramEnabled( false );
            }
            getCurrentTreePanel().setPhylogenyGraphicsType( getOptions().getPhylogenyGraphicsType() );
            MainFrame.updateScreenTextAntialias( getMainPanel().getTreePanels() );
        }
    }

    void updateOptions( final Options options ) {
        options.setAntialiasScreen( ( _screen_antialias_cbmi != null ) && _screen_antialias_cbmi.isSelected() );
        options.setBackgroundColorGradient( ( _background_gradient_cbmi != null )
                && _background_gradient_cbmi.isSelected() );
        options.setShowDomainLabels( ( _show_domain_labels != null ) && _show_domain_labels.isSelected() );
        options.setShowAnnotationRefSource( ( _show_annotation_ref_source != null )
                && _show_annotation_ref_source.isSelected() );
        options.setAbbreviateScientificTaxonNames( ( _abbreviate_scientific_names != null )
                && _abbreviate_scientific_names.isSelected() );
        options.setColorLabelsSameAsParentBranch( ( _color_labels_same_as_parent_branch != null )
                && _color_labels_same_as_parent_branch.isSelected() );
        options.setShowDefaultNodeShapesInternal( ( _show_default_node_shapes_internal_cbmi != null )
                && _show_default_node_shapes_internal_cbmi.isSelected() );
        options.setShowDefaultNodeShapesExternal( ( _show_default_node_shapes_external_cbmi != null )
                && _show_default_node_shapes_external_cbmi.isSelected() );
        if ( ( _non_lined_up_cladograms_rbmi != null ) && ( _non_lined_up_cladograms_rbmi.isSelected() ) ) {
            options.setCladogramType( CLADOGRAM_TYPE.NON_LINED_UP );
        }
        else if ( ( _uniform_cladograms_rbmi != null ) && ( _uniform_cladograms_rbmi.isSelected() ) ) {
            options.setCladogramType( CLADOGRAM_TYPE.TOTAL_NODE_SUM_DEP );
        }
        else if ( ( _ext_node_dependent_cladogram_rbmi != null ) && ( _ext_node_dependent_cladogram_rbmi.isSelected() ) ) {
            options.setCladogramType( CLADOGRAM_TYPE.EXT_NODE_SUM_DEP );
        }
        options.setSearchCaseSensitive( ( _search_case_senstive_cbmi != null )
                && _search_case_senstive_cbmi.isSelected() );
        if ( ( _show_scale_cbmi != null ) && _show_scale_cbmi.isEnabled() ) {
            options.setShowScale( _show_scale_cbmi.isSelected() );
        }
        if ( _label_direction_cbmi != null ) {
            if ( _label_direction_cbmi.isSelected() ) {
                options.setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
            }
            else {
                options.setNodeLabelDirection( NODE_LABEL_DIRECTION.HORIZONTAL );
            }
        }
        options.setShowOverview( ( _show_overview_cbmi != null ) && _show_overview_cbmi.isSelected() );
        options.setShowConfidenceStddev( ( _show_confidence_stddev_cbmi != null )
                && _show_confidence_stddev_cbmi.isSelected() );
        if ( ( _show_branch_length_values_cbmi != null ) && _show_branch_length_values_cbmi.isEnabled() ) {
            options.setShowBranchLengthValues( _show_branch_length_values_cbmi.isSelected() );
        }
        options.setMatchWholeTermsOnly( ( _search_whole_words_only_cbmi != null )
                && _search_whole_words_only_cbmi.isSelected() );
        options.setInverseSearchResult( ( _inverse_search_result_cbmi != null )
                && _inverse_search_result_cbmi.isSelected() );
        if ( ( _rectangular_type_cbmi != null ) && _rectangular_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        }
        else if ( ( _triangular_type_cbmi != null ) && _triangular_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR );
        }
        else if ( ( _curved_type_cbmi != null ) && _curved_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CURVED );
        }
        else if ( ( _convex_type_cbmi != null ) && _convex_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CONVEX );
        }
        else if ( ( _euro_type_cbmi != null ) && _euro_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE );
        }
        else if ( ( _rounded_type_cbmi != null ) && _rounded_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.ROUNDED );
        }
        else if ( ( _unrooted_type_cbmi != null ) && _unrooted_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
        }
        else if ( ( _circular_type_cbmi != null ) && _circular_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        }
        if ( ( _color_by_taxonomic_group_cbmi != null ) && _color_by_taxonomic_group_cbmi.isEnabled() ) {
            options.setColorByTaxonomicGroup( _color_by_taxonomic_group_cbmi.isSelected() );
        }
        if ( ( _right_line_up_domains_cbmi != null ) && _right_line_up_domains_cbmi.isEnabled() ) {
            options.setRightLineUpDomains( _right_line_up_domains_cbmi.isSelected() );
        }
        if ( ( _line_up_renderable_data_cbmi != null ) && _line_up_renderable_data_cbmi.isEnabled() ) {
            options.setLineUpRendarableNodeData( _line_up_renderable_data_cbmi.isSelected() );
        }
    }

    void updateTypeCheckboxes( final Options options, final Object o ) {
        setTypeMenuToAllUnselected();
        ( ( JCheckBoxMenuItem ) o ).setSelected( true );
    }

    void viewAsNexus() {
        if ( ( getMainPanel().getCurrentPhylogeny() != null ) && !getMainPanel().getCurrentPhylogeny().isEmpty() ) {
            String title = "Nexus";
            if ( !ForesterUtil.isEmpty( getMainPanel().getCurrentPhylogeny().getName() ) ) {
                title = "\"" + getMainPanel().getCurrentPhylogeny().getName() + "\" in " + title;
            }
            showTextFrame( getMainPanel().getCurrentPhylogeny().toNexus( getOptions()
                                   .getNhConversionSupportValueStyle() ),
                           title );
        }
    }

    void viewAsNH() {
        if ( ( getMainPanel().getCurrentPhylogeny() != null ) && !getMainPanel().getCurrentPhylogeny().isEmpty() ) {
            String title = "New Hampshire";
            if ( !ForesterUtil.isEmpty( getMainPanel().getCurrentPhylogeny().getName() ) ) {
                title = "\"" + getMainPanel().getCurrentPhylogeny().getName() + "\" in " + title;
            }
            showTextFrame( getMainPanel().getCurrentPhylogeny().toNewHampshire( getOptions()
                                   .getNhConversionSupportValueStyle() ),
                           title );
        }
    }

    void viewAsXML() {
        if ( ( getMainPanel().getCurrentPhylogeny() != null ) && !getMainPanel().getCurrentPhylogeny().isEmpty() ) {
            String title = "phyloXML";
            if ( !ForesterUtil.isEmpty( getMainPanel().getCurrentPhylogeny().getName() ) ) {
                title = "\"" + getMainPanel().getCurrentPhylogeny().getName() + "\" in " + title;
            }
            showTextFrame( getMainPanel().getCurrentPhylogeny().toPhyloXML( 0 ), title );
        }
    }

    private void chooseFont() {
        final FontChooser fc = new FontChooser();
        fc.setFont( getMainPanel().getTreeFontSet().getLargeFont() );
        fc.showDialog( this, "Select the Base Font" );
        getMainPanel().getTreeFontSet().setBaseFont( fc.getFont() );
    }

    private void chooseMinimalConfidence() {
        final String s = ( String ) JOptionPane
                .showInputDialog( this,
                                  "Please the minimum for confidence values to be displayed.\n" + "[current value: "
                                          + getOptions().getMinConfidenceValue() + "]\n",
                                  "Minimal Confidence Value",
                                  JOptionPane.QUESTION_MESSAGE,
                                  null,
                                  null,
                                  getOptions().getMinConfidenceValue() );
        if ( !ForesterUtil.isEmpty( s ) ) {
            boolean success = true;
            double m = 0.0;
            final String m_str = s.trim();
            if ( !ForesterUtil.isEmpty( m_str ) ) {
                try {
                    m = Double.parseDouble( m_str );
                }
                catch ( final Exception ex ) {
                    success = false;
                }
            }
            else {
                success = false;
            }
            if ( success && ( m >= 0.0 ) ) {
                getOptions().setMinConfidenceValue( m );
            }
        }
    }

    private void customizeRadioButtonMenuItem( final JRadioButtonMenuItem item, final boolean is_selected ) {
        if ( item != null ) {
            item.setFont( MainFrame.menu_font );
            if ( !getConfiguration().isUseNativeUI() ) {
                item.setBackground( getConfiguration().getGuiMenuBackgroundColor() );
                item.setForeground( getConfiguration().getGuiMenuTextColor() );
            }
            item.setSelected( is_selected );
            item.addActionListener( this );
        }
    }

    private MainPanel getMainPanel() {
        return _mainpanel;
    }

    private Phylogeny getSpeciesTree() {
        return _species_tree;
    }

    private boolean isScreenAntialias() {
        return true;
    }

    private void removeBranchColors() {
        if ( getMainPanel().getCurrentPhylogeny() != null ) {
            AptxUtil.removeBranchColors( getMainPanel().getCurrentPhylogeny() );
        }
    }

    private void removeVisualStyles() {
        if ( getMainPanel().getCurrentPhylogeny() != null ) {
            AptxUtil.removeVisualStyles( getMainPanel().getCurrentPhylogeny() );
        }
    }

    private void setMainPanel( final MainPanelApplets main_panel ) {
        _mainpanel = main_panel;
    }

    private void setSpeciesTree( final Phylogeny species_tree ) {
        _species_tree = species_tree;
    }

    private void setupUI() {
        try {
            if ( getConfiguration().isUseNativeUI() ) {
                UIManager.setLookAndFeel( UIManager.getSystemLookAndFeelClassName() );
            }
            else {
                UIManager.setLookAndFeel( UIManager.getCrossPlatformLookAndFeelClassName() );
            }
        }
        catch ( final UnsupportedLookAndFeelException e ) {
            AptxUtil.dieWithSystemError( "UnsupportedLookAndFeelException: " + e.toString() );
        }
        catch ( final ClassNotFoundException e ) {
            AptxUtil.dieWithSystemError( "ClassNotFoundException: " + e.toString() );
        }
        catch ( final InstantiationException e ) {
            AptxUtil.dieWithSystemError( "InstantiationException: " + e.toString() );
        }
        catch ( final IllegalAccessException e ) {
            AptxUtil.dieWithSystemError( "IllegalAccessException: " + e.toString() );
        }
        catch ( final Exception e ) {
            AptxUtil.dieWithSystemError( e.toString() );
        }
    }

    static void setupScreenTextAntialias( final List<TreePanel> treepanels, final boolean antialias ) {
        for( final TreePanel tree_panel : treepanels ) {
            tree_panel.setTextAntialias();
        }
    }
}