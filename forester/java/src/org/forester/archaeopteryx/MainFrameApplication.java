// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.WindowConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.plaf.synth.SynthLookAndFeel;

import org.forester.analysis.TaxonomyDataManager;
import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.tools.AncestralTaxonomyInferrer;
import org.forester.archaeopteryx.tools.InferenceManager;
import org.forester.archaeopteryx.tools.PhyloInferenceDialog;
import org.forester.archaeopteryx.tools.PhylogeneticInferenceOptions;
import org.forester.archaeopteryx.tools.PhylogeneticInferrer;
import org.forester.archaeopteryx.tools.SequenceDataRetriver;
import org.forester.archaeopteryx.webservices.PhylogeniesWebserviceClient;
import org.forester.archaeopteryx.webservices.WebservicesManager;
import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.io.parsers.tol.TolParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.io.writers.SequenceWriter;
import org.forester.msa.Msa;
import org.forester.msa.MsaFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sequence.Sequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;
import org.forester.util.WindowsUtils;

public final class MainFrameApplication extends MainFrame {

    static final String                      INFER_ANCESTOR_TAXONOMIES             = "Infer Ancestor Taxonomies";
    static final String                      OBTAIN_DETAILED_TAXONOMIC_INFORMATION = "Obtain Detailed Taxonomic Information";
    private final static int                 FRAME_X_SIZE                          = 800;
    private final static int                 FRAME_Y_SIZE                          = 800;
    // Filters for the file-open dialog (classes defined in this file)
    private final static NHFilter            nhfilter                              = new NHFilter();
    private final static NHXFilter           nhxfilter                             = new NHXFilter();
    private final static XMLFilter           xmlfilter                             = new XMLFilter();
    private final static TolFilter           tolfilter                             = new TolFilter();
    private final static NexusFilter         nexusfilter                           = new NexusFilter();
    private final static PdfFilter           pdffilter                             = new PdfFilter();
    private final static GraphicsFileFilter  graphicsfilefilter                    = new GraphicsFileFilter();
    private final static MsaFileFilter       msafilter                             = new MsaFileFilter();
    private final static SequencesFileFilter seqsfilter                            = new SequencesFileFilter();
    private final static DefaultFilter       defaultfilter                         = new DefaultFilter();
    private static final long                serialVersionUID                      = -799735726778865234L;
    private final JFileChooser               _values_filechooser;
    private final JFileChooser               _sequences_filechooser;
    private final JFileChooser               _open_filechooser;
    private final JFileChooser               _msa_filechooser;
    private final JFileChooser               _seqs_pi_filechooser;
    private final JFileChooser               _open_filechooser_for_species_tree;
    private final JFileChooser               _save_filechooser;
    private final JFileChooser               _writetopdf_filechooser;
    private final JFileChooser               _writetographics_filechooser;
    // Application-only print menu items
    private JMenuItem                        _print_item;
    private JMenuItem                        _write_to_pdf_item;
    private JMenuItem                        _write_to_jpg_item;
    private JMenuItem                        _write_to_gif_item;
    private JMenuItem                        _write_to_tif_item;
    private JMenuItem                        _write_to_png_item;
    private JMenuItem                        _write_to_bmp_item;
    private File                             _current_dir;
    private ButtonGroup                      _radio_group_1;
    private ButtonGroup                      _radio_group_2;
    // Others:
    double                                   _min_not_collapse                     = Constants.MIN_NOT_COLLAPSE_DEFAULT;
    // Phylogeny Inference menu
    private JMenu                            _inference_menu;
    private JMenuItem                        _inference_from_msa_item;
    private JMenuItem                        _inference_from_seqs_item;
    // Phylogeny Inference
    private PhylogeneticInferenceOptions     _phylogenetic_inference_options       = null;
    private Msa                              _msa                                  = null;
    private File                             _msa_file                             = null;
    private List<Sequence>                   _seqs                                 = null;
    private File                             _seqs_file                            = null;
    JMenuItem                                _read_values_jmi;
    JMenuItem                                _read_seqs_jmi;

    private MainFrameApplication( final Phylogeny[] phys, final Configuration config ) {
        _configuration = config;
        if ( _configuration == null ) {
            throw new IllegalArgumentException( "configuration is null" );
        }
        setVisible( false );
        setOptions( Options.createInstance( _configuration ) );
        _mainpanel = new MainPanel( _configuration, this );
        _open_filechooser = null;
        _open_filechooser_for_species_tree = null;
        _save_filechooser = null;
        _writetopdf_filechooser = null;
        _writetographics_filechooser = null;
        _msa_filechooser = null;
        _seqs_pi_filechooser = null;
        _values_filechooser = null;
        _sequences_filechooser = null;
        _jmenubar = new JMenuBar();
        buildFileMenu();
        buildTypeMenu();
        _contentpane = getContentPane();
        _contentpane.setLayout( new BorderLayout() );
        _contentpane.add( _mainpanel, BorderLayout.CENTER );
        // App is this big
        setSize( MainFrameApplication.FRAME_X_SIZE, MainFrameApplication.FRAME_Y_SIZE );
        // The window listener
        setDefaultCloseOperation( WindowConstants.DO_NOTHING_ON_CLOSE );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                exit();
            }
        } );
        //   setVisible( true );
        if ( ( phys != null ) && ( phys.length > 0 ) ) {
            AptxUtil.addPhylogeniesToTabs( phys, "", null, _configuration, _mainpanel );
            validate();
            getMainPanel().getControlPanel().showWholeAll();
            getMainPanel().getControlPanel().showWhole();
        }
        //activateSaveAllIfNeeded();
        // ...and its children
        _contentpane.repaint();
    }

    private MainFrameApplication( final Phylogeny[] phys, final Configuration config, final String title ) {
        this( phys, config, title, null );
    }

    private MainFrameApplication( final Phylogeny[] phys,
                                  final Configuration config,
                                  final String title,
                                  final File current_dir ) {
        super();
        _configuration = config;
        if ( _configuration == null ) {
            throw new IllegalArgumentException( "configuration is null" );
        }
        try {
            boolean synth_exception = false;
            if ( Constants.__SYNTH_LF ) {
                try {
                    final SynthLookAndFeel synth = new SynthLookAndFeel();
                    synth.load( MainFrameApplication.class.getResourceAsStream( "/resources/synth_look_and_feel_1.xml" ),
                                MainFrameApplication.class );
                    UIManager.setLookAndFeel( synth );
                }
                catch ( final Exception ex ) {
                    synth_exception = true;
                    ForesterUtil.printWarningMessage( Constants.PRG_NAME,
                                                      "could not create synth look and feel: "
                                                              + ex.getLocalizedMessage() );
                }
            }
            if ( !Constants.__SYNTH_LF || synth_exception ) {
                if ( _configuration.isUseNativeUI() ) {
                    UIManager.setLookAndFeel( UIManager.getSystemLookAndFeelClassName() );
                }
                else {
                    UIManager.setLookAndFeel( UIManager.getCrossPlatformLookAndFeelClassName() );
                }
            }
            //UIManager.setLookAndFeel( "com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel" );
        }
        catch ( final UnsupportedLookAndFeelException e ) {
            AptxUtil.dieWithSystemError( "unsupported look and feel: " + e.toString() );
        }
        catch ( final ClassNotFoundException e ) {
            AptxUtil.dieWithSystemError( "class not found exception: " + e.toString() );
        }
        catch ( final InstantiationException e ) {
            AptxUtil.dieWithSystemError( "instantiation exception: " + e.toString() );
        }
        catch ( final IllegalAccessException e ) {
            AptxUtil.dieWithSystemError( "illegal access exception: " + e.toString() );
        }
        if ( ( current_dir != null ) && current_dir.canRead() && current_dir.isDirectory() ) {
            setCurrentDir( current_dir );
        }
        // hide until everything is ready
        setVisible( false );
        setOptions( Options.createInstance( _configuration ) );
        setInferenceManager( InferenceManager.createInstance( _configuration ) );
        setPhylogeneticInferenceOptions( PhylogeneticInferenceOptions.createInstance( _configuration ) );
        //     _textframe = null; #~~~~
        // set title
        setTitle( Constants.PRG_NAME + " " + Constants.VERSION + " (" + Constants.PRG_DATE + ")" );
        _mainpanel = new MainPanel( _configuration, this );
        // The file dialogs
        _open_filechooser = new JFileChooser();
        _open_filechooser.setCurrentDirectory( new File( "." ) );
        _open_filechooser.setMultiSelectionEnabled( false );
        _open_filechooser.addChoosableFileFilter( MainFrameApplication.xmlfilter );
        _open_filechooser.addChoosableFileFilter( MainFrameApplication.nhxfilter );
        _open_filechooser.addChoosableFileFilter( MainFrameApplication.nhfilter );
        _open_filechooser.addChoosableFileFilter( MainFrameApplication.nexusfilter );
        _open_filechooser.addChoosableFileFilter( MainFrameApplication.tolfilter );
        _open_filechooser.addChoosableFileFilter( _open_filechooser.getAcceptAllFileFilter() );
        _open_filechooser.setFileFilter( MainFrameApplication.defaultfilter );
        _open_filechooser_for_species_tree = new JFileChooser();
        _open_filechooser_for_species_tree.setCurrentDirectory( new File( "." ) );
        _open_filechooser_for_species_tree.setMultiSelectionEnabled( false );
        _open_filechooser_for_species_tree.addChoosableFileFilter( MainFrameApplication.xmlfilter );
        _open_filechooser_for_species_tree.addChoosableFileFilter( MainFrameApplication.tolfilter );
        _open_filechooser_for_species_tree.setFileFilter( MainFrameApplication.xmlfilter );
        _save_filechooser = new JFileChooser();
        _save_filechooser.setCurrentDirectory( new File( "." ) );
        _save_filechooser.setMultiSelectionEnabled( false );
        _save_filechooser.setFileFilter( MainFrameApplication.xmlfilter );
        _save_filechooser.addChoosableFileFilter( MainFrameApplication.nhfilter );
        _save_filechooser.addChoosableFileFilter( MainFrameApplication.nexusfilter );
        _save_filechooser.addChoosableFileFilter( _save_filechooser.getAcceptAllFileFilter() );
        _writetopdf_filechooser = new JFileChooser();
        _writetopdf_filechooser.addChoosableFileFilter( MainFrameApplication.pdffilter );
        _writetographics_filechooser = new JFileChooser();
        _writetographics_filechooser.addChoosableFileFilter( MainFrameApplication.graphicsfilefilter );
        // Msa:
        _msa_filechooser = new JFileChooser();
        _msa_filechooser.setName( "Read Multiple Sequence Alignment File" );
        _msa_filechooser.setCurrentDirectory( new File( "." ) );
        _msa_filechooser.setMultiSelectionEnabled( false );
        _msa_filechooser.addChoosableFileFilter( _msa_filechooser.getAcceptAllFileFilter() );
        _msa_filechooser.addChoosableFileFilter( MainFrameApplication.msafilter );
        // Seqs:
        _seqs_pi_filechooser = new JFileChooser();
        _seqs_pi_filechooser.setName( "Read Sequences File" );
        _seqs_pi_filechooser.setCurrentDirectory( new File( "." ) );
        _seqs_pi_filechooser.setMultiSelectionEnabled( false );
        _seqs_pi_filechooser.addChoosableFileFilter( _seqs_pi_filechooser.getAcceptAllFileFilter() );
        _seqs_pi_filechooser.addChoosableFileFilter( MainFrameApplication.seqsfilter );
        // Expression
        _values_filechooser = new JFileChooser();
        _values_filechooser.setCurrentDirectory( new File( "." ) );
        _values_filechooser.setMultiSelectionEnabled( false );
        // Sequences
        _sequences_filechooser = new JFileChooser();
        _sequences_filechooser.setCurrentDirectory( new File( "." ) );
        _sequences_filechooser.setMultiSelectionEnabled( false );
        // build the menu bar
        _jmenubar = new JMenuBar();
        if ( !_configuration.isUseNativeUI() ) {
            _jmenubar.setBackground( getConfiguration().getGuiMenuBackgroundColor() );
        }
        buildFileMenu();
        if ( Constants.__ALLOW_PHYLOGENETIC_INFERENCE ) {
            buildPhylogeneticInferenceMenu();
        }
        buildAnalysisMenu();
        buildToolsMenu();
        buildViewMenu();
        buildFontSizeMenu();
        buildOptionsMenu();
        buildTypeMenu();
        buildHelpMenu();
        setJMenuBar( _jmenubar );
        _jmenubar.add( _help_jmenu );
        _contentpane = getContentPane();
        _contentpane.setLayout( new BorderLayout() );
        _contentpane.add( _mainpanel, BorderLayout.CENTER );
        // App is this big
        setSize( MainFrameApplication.FRAME_X_SIZE, MainFrameApplication.FRAME_Y_SIZE );
        //        addWindowFocusListener( new WindowAdapter() {
        //
        //            @Override
        //            public void windowGainedFocus( WindowEvent e ) {
        //                requestFocusInWindow();
        //            }
        //        } );
        // The window listener
        setDefaultCloseOperation( WindowConstants.DO_NOTHING_ON_CLOSE );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                if ( isUnsavedDataPresent() ) {
                    final int r = JOptionPane.showConfirmDialog( null,
                                                                 "Exit despite potentially unsaved changes?",
                                                                 "Exit?",
                                                                 JOptionPane.YES_NO_OPTION );
                    if ( r != JOptionPane.YES_OPTION ) {
                        return;
                    }
                }
                else {
                    final int r = JOptionPane.showConfirmDialog( null,
                                                                 "Exit Archaeopteryx?",
                                                                 "Exit?",
                                                                 JOptionPane.YES_NO_OPTION );
                    if ( r != JOptionPane.YES_OPTION ) {
                        return;
                    }
                }
                exit();
            }
        } );
        // The component listener
        addComponentListener( new ComponentAdapter() {

            @Override
            public void componentResized( final ComponentEvent e ) {
                if ( _mainpanel.getCurrentTreePanel() != null ) {
                    _mainpanel.getCurrentTreePanel().calcParametersForPainting( _mainpanel.getCurrentTreePanel()
                                                                                        .getWidth(),
                                                                                _mainpanel.getCurrentTreePanel()
                                                                                        .getHeight(),
                                                                                getOptions().isAllowFontSizeChange() );
                }
            }
        } );
        requestFocusInWindow();
        // addKeyListener( this );
        setVisible( true );
        if ( ( phys != null ) && ( phys.length > 0 ) ) {
            AptxUtil.addPhylogeniesToTabs( phys, title, null, _configuration, _mainpanel );
            validate();
            getMainPanel().getControlPanel().showWholeAll();
            getMainPanel().getControlPanel().showWhole();
        }
        activateSaveAllIfNeeded();
        // ...and its children
        _contentpane.repaint();
        System.gc();
    }

    private MainFrameApplication( final Phylogeny[] phys, final String config_file, final String title ) {
        // Reads the config file (false, false => not url, not applet):
        this( phys, new Configuration( config_file, false, false, true ), title );
    }

    @Override
    public void actionPerformed( final ActionEvent e ) {
        try {
            super.actionPerformed( e );
            final Object o = e.getSource();
            // Handle app-specific actions here:
            if ( o == _open_item ) {
                readPhylogeniesFromFile();
            }
            else if ( o == _save_item ) {
                writeToFile( _mainpanel.getCurrentPhylogeny() );
                // If subtree currently displayed, save it, instead of complete
                // tree.
            }
            else if ( o == _new_item ) {
                newTree();
            }
            else if ( o == _save_all_item ) {
                writeAllToFile();
            }
            else if ( o == _close_item ) {
                closeCurrentPane();
            }
            else if ( o == _write_to_pdf_item ) {
                writeToPdf( _mainpanel.getCurrentPhylogeny() );
            }
            else if ( o == _write_to_jpg_item ) {
                writeToGraphicsFile( _mainpanel.getCurrentPhylogeny(), GraphicsExportType.JPG );
            }
            else if ( o == _write_to_png_item ) {
                writeToGraphicsFile( _mainpanel.getCurrentPhylogeny(), GraphicsExportType.PNG );
            }
            else if ( o == _write_to_gif_item ) {
                writeToGraphicsFile( _mainpanel.getCurrentPhylogeny(), GraphicsExportType.GIF );
            }
            else if ( o == _write_to_tif_item ) {
                writeToGraphicsFile( _mainpanel.getCurrentPhylogeny(), GraphicsExportType.TIFF );
            }
            else if ( o == _write_to_bmp_item ) {
                writeToGraphicsFile( _mainpanel.getCurrentPhylogeny(), GraphicsExportType.BMP );
            }
            else if ( o == _print_item ) {
                print();
            }
            else if ( o == _load_species_tree_item ) {
                readSpeciesTreeFromFile();
            }
            else if ( o == _lineage_inference ) {
                if ( isSubtreeDisplayed() ) {
                    JOptionPane.showMessageDialog( this,
                                                   "Subtree is shown.",
                                                   "Cannot infer ancestral taxonomies",
                                                   JOptionPane.ERROR_MESSAGE );
                    return;
                }
                executeLineageInference();
            }
            else if ( o == _obtain_detailed_taxonomic_information_jmi ) {
                if ( isSubtreeDisplayed() ) {
                    return;
                }
                obtainDetailedTaxonomicInformation();
            }
            else if ( o == _obtain_detailed_taxonomic_information_deleting_jmi ) {
                if ( isSubtreeDisplayed() ) {
                    return;
                }
                obtainDetailedTaxonomicInformationDelete();
            }
            else if ( o == _obtain_seq_information_jmi ) {
                obtainSequenceInformation();
            }
            else if ( o == _read_values_jmi ) {
                if ( isSubtreeDisplayed() ) {
                    return;
                }
                addExpressionValuesFromFile();
            }
            else if ( o == _read_seqs_jmi ) {
                if ( isSubtreeDisplayed() ) {
                    return;
                }
                addSequencesFromFile();
            }
            else if ( o == _move_node_names_to_tax_sn_jmi ) {
                moveNodeNamesToTaxSn();
            }
            else if ( o == _move_node_names_to_seq_names_jmi ) {
                moveNodeNamesToSeqNames();
            }
            else if ( o == _extract_tax_code_from_node_names_jmi ) {
                extractTaxDataFromNodeNames();
            }
            else if ( o == _graphics_export_visible_only_cbmi ) {
                updateOptions( getOptions() );
            }
            else if ( o == _antialias_print_cbmi ) {
                updateOptions( getOptions() );
            }
            else if ( o == _print_black_and_white_cbmi ) {
                updateOptions( getOptions() );
            }
            else if ( o == _print_using_actual_size_cbmi ) {
                updateOptions( getOptions() );
            }
            else if ( o == _graphics_export_using_actual_size_cbmi ) {
                updateOptions( getOptions() );
            }
            else if ( o == _print_size_mi ) {
                choosePrintSize();
            }
            else if ( o == _choose_pdf_width_mi ) {
                choosePdfWidth();
            }
            else if ( o == _internal_number_are_confidence_for_nh_parsing_cbmi ) {
                updateOptions( getOptions() );
            }
            else if ( o == _replace_underscores_cbmi ) {
                if ( ( _extract_taxonomy_no_rbmi != null ) && !_extract_taxonomy_no_rbmi.isSelected() ) {
                    _extract_taxonomy_no_rbmi.setSelected( true );
                }
                updateOptions( getOptions() );
            }
            else if ( o == _allow_errors_in_distance_to_parent_cbmi ) {
                updateOptions( getOptions() );
            }
            else if ( o == _collapse_below_threshold ) {
                if ( isSubtreeDisplayed() ) {
                    return;
                }
                collapseBelowThreshold();
            }
            else if ( ( o == _extract_taxonomy_pfam_strict_rbmi ) || ( o == _extract_taxonomy_pfam_relaxed_rbmi )
                    || ( o == _extract_taxonomy_agressive_rbmi ) ) {
                if ( _replace_underscores_cbmi != null ) {
                    _replace_underscores_cbmi.setSelected( false );
                }
                updateOptions( getOptions() );
            }
            else if ( o == _extract_taxonomy_no_rbmi ) {
                updateOptions( getOptions() );
            }
            else if ( o == _inference_from_msa_item ) {
                executePhyleneticInference( false );
            }
            else if ( o == _inference_from_seqs_item ) {
                executePhyleneticInference( true );
            }
            _contentpane.repaint();
        }
        catch ( final Exception ex ) {
            AptxUtil.unexpectedException( ex );
        }
        catch ( final Error err ) {
            AptxUtil.unexpectedError( err );
        }
    }

    public void end() {
        _mainpanel.terminate();
        _contentpane.removeAll();
        setVisible( false );
        dispose();
    }

    @Override
    public MainPanel getMainPanel() {
        return _mainpanel;
    }

    public Msa getMsa() {
        return _msa;
    }

    public File getMsaFile() {
        return _msa_file;
    }

    public List<Sequence> getSeqs() {
        return _seqs;
    }

    public File getSeqsFile() {
        return _seqs_file;
    }

    public void readMsaFromFile() {
        // Set an initial directory if none set yet
        final File my_dir = getCurrentDir();
        _msa_filechooser.setMultiSelectionEnabled( false );
        // Open file-open dialog and set current directory
        if ( my_dir != null ) {
            _msa_filechooser.setCurrentDirectory( my_dir );
        }
        final int result = _msa_filechooser.showOpenDialog( _contentpane );
        // All done: get the msa
        final File file = _msa_filechooser.getSelectedFile();
        setCurrentDir( _msa_filechooser.getCurrentDirectory() );
        if ( ( file != null ) && !file.isDirectory() && ( result == JFileChooser.APPROVE_OPTION ) ) {
            setMsaFile( null );
            setMsa( null );
            Msa msa = null;
            try {
                final InputStream is = new FileInputStream( file );
                if ( FastaParser.isLikelyFasta( file ) ) {
                    msa = FastaParser.parseMsa( is );
                }
                else {
                    msa = GeneralMsaParser.parse( is );
                }
            }
            catch ( final MsaFormatException e ) {
                setArrowCursor();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Multiple sequence alignment format error",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            catch ( final IOException e ) {
                setArrowCursor();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Failed to read multiple sequence alignment",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            catch ( final IllegalArgumentException e ) {
                setArrowCursor();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Unexpected error during reading of multiple sequence alignment",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            catch ( final Exception e ) {
                setArrowCursor();
                e.printStackTrace();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Unexpected error during reading of multiple sequence alignment",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            if ( ( msa == null ) || ( msa.getNumberOfSequences() < 1 ) ) {
                JOptionPane.showMessageDialog( this,
                                               "Multiple sequence alignment is empty",
                                               "Illegal Multiple Sequence Alignment",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            if ( msa.getNumberOfSequences() < 4 ) {
                JOptionPane.showMessageDialog( this,
                                               "Multiple sequence alignment needs to contain at least 3 sequences",
                                               "Illegal multiple sequence alignment",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            if ( msa.getLength() < 2 ) {
                JOptionPane.showMessageDialog( this,
                                               "Multiple sequence alignment needs to contain at least 2 residues",
                                               "Illegal multiple sequence alignment",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            System.gc();
            setMsaFile( _msa_filechooser.getSelectedFile() );
            setMsa( msa );
        }
    }

    public void readSeqsFromFileforPI() {
        // Set an initial directory if none set yet
        final File my_dir = getCurrentDir();
        _seqs_pi_filechooser.setMultiSelectionEnabled( false );
        // Open file-open dialog and set current directory
        if ( my_dir != null ) {
            _seqs_pi_filechooser.setCurrentDirectory( my_dir );
        }
        final int result = _seqs_pi_filechooser.showOpenDialog( _contentpane );
        // All done: get the seqs
        final File file = _seqs_pi_filechooser.getSelectedFile();
        setCurrentDir( _seqs_pi_filechooser.getCurrentDirectory() );
        if ( ( file != null ) && !file.isDirectory() && ( result == JFileChooser.APPROVE_OPTION ) ) {
            setSeqsFile( null );
            setSeqs( null );
            List<Sequence> seqs = null;
            try {
                if ( FastaParser.isLikelyFasta( new FileInputStream( file ) ) ) {
                    seqs = FastaParser.parse( new FileInputStream( file ) );
                    for( final Sequence seq : seqs ) {
                        System.out.println( SequenceWriter.toFasta( seq, 60 ) );
                    }
                }
                else {
                    //TODO error
                }
            }
            catch ( final MsaFormatException e ) {
                setArrowCursor();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Multiple sequence file format error",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            catch ( final IOException e ) {
                setArrowCursor();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Failed to read multiple sequence file",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            catch ( final IllegalArgumentException e ) {
                setArrowCursor();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Unexpected error during reading of multiple sequence file",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            catch ( final Exception e ) {
                setArrowCursor();
                e.printStackTrace();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Unexpected error during reading of multiple sequence file",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            if ( ( seqs == null ) || ( seqs.size() < 1 ) ) {
                JOptionPane.showMessageDialog( this,
                                               "Multiple sequence file is empty",
                                               "Illegal multiple sequence file",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            if ( seqs.size() < 4 ) {
                JOptionPane.showMessageDialog( this,
                                               "Multiple sequence file needs to contain at least 3 sequences",
                                               "Illegal multiple sequence file",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            //  if ( msa.getLength() < 2 ) {
            //       JOptionPane.showMessageDialog( this,
            //                                      "Multiple sequence alignment needs to contain at least 2 residues",
            //                                      "Illegal multiple sequence file",
            //                                      JOptionPane.ERROR_MESSAGE );
            //       return;
            //   }
            System.gc();
            setSeqsFile( _seqs_pi_filechooser.getSelectedFile() );
            setSeqs( seqs );
        }
    }

    void buildAnalysisMenu() {
        _analysis_menu = MainFrame.createMenu( "Analysis", getConfiguration() );
        _analysis_menu.add( _gsdi_item = new JMenuItem( "GSDI (Generalized Speciation Duplication Inference)" ) );
        _analysis_menu.add( _gsdir_item = new JMenuItem( "GSDIR (GSDI with re-rooting)" ) );
        _analysis_menu.add( _load_species_tree_item = new JMenuItem( "Load Species Tree..." ) );
        customizeJMenuItem( _gsdi_item );
        customizeJMenuItem( _gsdir_item );
        customizeJMenuItem( _load_species_tree_item );
        _analysis_menu.addSeparator();
        _analysis_menu.add( _lineage_inference = new JMenuItem( INFER_ANCESTOR_TAXONOMIES ) );
        customizeJMenuItem( _lineage_inference );
        _lineage_inference.setToolTipText( "Inference of ancestor taxonomies/lineages" );
        _jmenubar.add( _analysis_menu );
    }

    @Override
    void buildFileMenu() {
        _file_jmenu = MainFrame.createMenu( "File", getConfiguration() );
        _file_jmenu.add( _open_item = new JMenuItem( "Read Tree from File..." ) );
        _file_jmenu.addSeparator();
        final WebservicesManager webservices_manager = WebservicesManager.getInstance();
        _load_phylogeny_from_webservice_menu_items = new JMenuItem[ webservices_manager
                .getAvailablePhylogeniesWebserviceClients().size() ];
        for( int i = 0; i < webservices_manager.getAvailablePhylogeniesWebserviceClients().size(); ++i ) {
            final PhylogeniesWebserviceClient client = webservices_manager.getAvailablePhylogeniesWebserviceClient( i );
            _load_phylogeny_from_webservice_menu_items[ i ] = new JMenuItem( client.getMenuName() );
            _file_jmenu.add( _load_phylogeny_from_webservice_menu_items[ i ] );
        }
        if ( getConfiguration().isEditable() ) {
            _file_jmenu.addSeparator();
            _file_jmenu.add( _new_item = new JMenuItem( "New" ) );
            _new_item.setToolTipText( "to create a new tree with one node, as source for manual tree construction" );
        }
        _file_jmenu.addSeparator();
        _file_jmenu.add( _save_item = new JMenuItem( "Save Tree As..." ) );
        _file_jmenu.add( _save_all_item = new JMenuItem( "Save All Trees As..." ) );
        _save_all_item.setToolTipText( "Write all phylogenies to one file." );
        _save_all_item.setEnabled( false );
        _file_jmenu.addSeparator();
        _file_jmenu.add( _write_to_pdf_item = new JMenuItem( "Export to PDF file ..." ) );
        if ( AptxUtil.canWriteFormat( "tif" ) || AptxUtil.canWriteFormat( "tiff" ) || AptxUtil.canWriteFormat( "TIF" ) ) {
            _file_jmenu.add( _write_to_tif_item = new JMenuItem( "Export to TIFF file..." ) );
        }
        _file_jmenu.add( _write_to_png_item = new JMenuItem( "Export to PNG file..." ) );
        _file_jmenu.add( _write_to_jpg_item = new JMenuItem( "Export to JPG file..." ) );
        if ( AptxUtil.canWriteFormat( "gif" ) ) {
            _file_jmenu.add( _write_to_gif_item = new JMenuItem( "Export to GIF file..." ) );
        }
        if ( AptxUtil.canWriteFormat( "bmp" ) ) {
            _file_jmenu.add( _write_to_bmp_item = new JMenuItem( "Export to BMP file..." ) );
        }
        _file_jmenu.addSeparator();
        _file_jmenu.add( _print_item = new JMenuItem( "Print..." ) );
        _file_jmenu.addSeparator();
        _file_jmenu.add( _close_item = new JMenuItem( "Close Tab" ) );
        _close_item.setToolTipText( "To close the current pane." );
        _close_item.setEnabled( true );
        _file_jmenu.addSeparator();
        _file_jmenu.add( _exit_item = new JMenuItem( "Exit" ) );
        // For print in color option item
        customizeJMenuItem( _open_item );
        _open_item
                .setFont( new Font( _open_item.getFont().getFontName(), Font.BOLD, _open_item.getFont().getSize() + 4 ) );
        for( int i = 0; i < webservices_manager.getAvailablePhylogeniesWebserviceClients().size(); ++i ) {
            customizeJMenuItem( _load_phylogeny_from_webservice_menu_items[ i ] );
        }
        customizeJMenuItem( _save_item );
        if ( getConfiguration().isEditable() ) {
            customizeJMenuItem( _new_item );
        }
        customizeJMenuItem( _close_item );
        customizeJMenuItem( _save_all_item );
        customizeJMenuItem( _write_to_pdf_item );
        customizeJMenuItem( _write_to_png_item );
        customizeJMenuItem( _write_to_jpg_item );
        customizeJMenuItem( _write_to_gif_item );
        customizeJMenuItem( _write_to_tif_item );
        customizeJMenuItem( _write_to_bmp_item );
        customizeJMenuItem( _print_item );
        customizeJMenuItem( _exit_item );
        _jmenubar.add( _file_jmenu );
    }

    void buildOptionsMenu() {
        _options_jmenu = MainFrame.createMenu( OPTIONS_HEADER, getConfiguration() );
        _options_jmenu.addChangeListener( new ChangeListener() {

            @Override
            public void stateChanged( final ChangeEvent e ) {
                MainFrame.setOvPlacementColorChooseMenuItem( _overview_placment_mi, getOptions() );
                MainFrame.setTextColorChooseMenuItem( _switch_colors_mi, getCurrentTreePanel() );
                MainFrame
                        .setTextMinSupportMenuItem( _choose_minimal_confidence_mi, getOptions(), getCurrentTreePanel() );
                MainFrame.setTextForFontChooserMenuItem( _choose_font_mi, MainFrame
                        .createCurrentFontDesc( getMainPanel().getTreeFontSet() ) );
                setTextForGraphicsSizeChooserMenuItem( _print_size_mi, getOptions() );
                setTextForPdfLineWidthChooserMenuItem( _choose_pdf_width_mi, getOptions() );
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
        _options_jmenu.add( customizeMenuItemAsLabel( new JMenuItem( DISPLAY_SUBHEADER ), getConfiguration() ) );
        _options_jmenu
                .add( _ext_node_dependent_cladogram_rbmi = new JRadioButtonMenuItem( MainFrame.NONUNIFORM_CLADOGRAMS_LABEL ) );
        _options_jmenu.add( _uniform_cladograms_rbmi = new JRadioButtonMenuItem( MainFrame.UNIFORM_CLADOGRAMS_LABEL ) );
        _options_jmenu.add( _non_lined_up_cladograms_rbmi = new JRadioButtonMenuItem( NON_LINED_UP_CLADOGRAMS_LABEL ) );
        _radio_group_1 = new ButtonGroup();
        _radio_group_1.add( _ext_node_dependent_cladogram_rbmi );
        _radio_group_1.add( _uniform_cladograms_rbmi );
        _radio_group_1.add( _non_lined_up_cladograms_rbmi );
        ///////
        _options_jmenu.add( _show_overview_cbmi = new JCheckBoxMenuItem( SHOW_OVERVIEW_LABEL ) );
        _options_jmenu.add( _show_scale_cbmi = new JCheckBoxMenuItem( DISPLAY_SCALE_LABEL ) );
        _options_jmenu
                .add( _show_branch_length_values_cbmi = new JCheckBoxMenuItem( DISPLAY_BRANCH_LENGTH_VALUES_LABEL ) );
        _options_jmenu
                .add( _show_default_node_shapes_internal_cbmi = new JCheckBoxMenuItem( DISPLAY_NODE_BOXES_LABEL_INT ) );
        _options_jmenu
                .add( _show_default_node_shapes_external_cbmi = new JCheckBoxMenuItem( DISPLAY_NODE_BOXES_LABEL_EXT ) );
        if ( getConfiguration().doDisplayOption( Configuration.show_domain_architectures ) ) {
            _options_jmenu.add( _show_domain_labels = new JCheckBoxMenuItem( SHOW_DOMAIN_LABELS_LABEL ) );
        }
        _options_jmenu.add( _show_annotation_ref_source = new JCheckBoxMenuItem( SHOW_ANN_REF_SOURCE_LABEL ) );
        _options_jmenu.add( _show_confidence_stddev_cbmi = new JCheckBoxMenuItem( SHOW_CONF_STDDEV_LABEL ) );
        _options_jmenu.add( _color_by_taxonomic_group_cbmi = new JCheckBoxMenuItem( COLOR_BY_TAXONOMIC_GROUP ) );
        _options_jmenu.add( _color_labels_same_as_parent_branch = new JCheckBoxMenuItem( COLOR_LABELS_LABEL ) );
        _color_labels_same_as_parent_branch.setToolTipText( MainFrame.COLOR_LABELS_TIP );
        _options_jmenu.add( _abbreviate_scientific_names = new JCheckBoxMenuItem( ABBREV_SN_LABEL ) );
        _options_jmenu.add( _label_direction_cbmi = new JCheckBoxMenuItem( LABEL_DIRECTION_LABEL ) );
        _label_direction_cbmi.setToolTipText( LABEL_DIRECTION_TIP );
        _options_jmenu.add( _screen_antialias_cbmi = new JCheckBoxMenuItem( SCREEN_ANTIALIAS_LABEL ) );
        _options_jmenu.add( _background_gradient_cbmi = new JCheckBoxMenuItem( BG_GRAD_LABEL ) );
        _options_jmenu.add( _cycle_node_shape_mi = new JMenuItem( MainFrame.CYCLE_NODE_SHAPE_LABEL ) );
        _options_jmenu.add( _cycle_node_fill_mi = new JMenuItem( MainFrame.CYCLE_NODE_FILL_LABEL ) );
        _options_jmenu.add( _choose_node_size_mi = new JMenuItem( MainFrame.CHOOSE_NODE_SIZE_LABEL ) );
        _options_jmenu.add( _choose_minimal_confidence_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _overview_placment_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _switch_colors_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _choose_font_mi = new JMenuItem( "" ) );
        ///////
        _options_jmenu.addSeparator();
        _options_jmenu.add( customizeMenuItemAsLabel( new JMenuItem( SEARCH_SUBHEADER ), getConfiguration() ) );
        _options_jmenu.add( _search_case_senstive_cbmi = new JCheckBoxMenuItem( SEARCH_CASE_SENSITIVE_LABEL ) );
        _options_jmenu.add( _search_whole_words_only_cbmi = new JCheckBoxMenuItem( SEARCH_TERMS_ONLY_LABEL ) );
        _options_jmenu.add( _inverse_search_result_cbmi = new JCheckBoxMenuItem( INVERSE_SEARCH_RESULT_LABEL ) );
        _options_jmenu.addSeparator();
        _options_jmenu.add( customizeMenuItemAsLabel( new JMenuItem( "Graphics Export & Printing:" ),
                                                      getConfiguration() ) );
        _options_jmenu.add( _antialias_print_cbmi = new JCheckBoxMenuItem( "Antialias" ) );
        _options_jmenu.add( _print_black_and_white_cbmi = new JCheckBoxMenuItem( "Export in Black and White" ) );
        _options_jmenu
                .add( _print_using_actual_size_cbmi = new JCheckBoxMenuItem( "Use Current Image Size for PDF export and Printing" ) );
        _options_jmenu
                .add( _graphics_export_using_actual_size_cbmi = new JCheckBoxMenuItem( "Use Current Image Size for PNG, JPG, and GIF export" ) );
        _options_jmenu
                .add( _graphics_export_visible_only_cbmi = new JCheckBoxMenuItem( "Limit to Visible ('Screenshot') for PNG, JPG, and GIF export" ) );
        _options_jmenu.add( _print_size_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _choose_pdf_width_mi = new JMenuItem( "" ) );
        _options_jmenu.addSeparator();
        _options_jmenu.add( customizeMenuItemAsLabel( new JMenuItem( "Newick/NHX/Nexus Input:" ), getConfiguration() ) );
        _options_jmenu
                .add( _internal_number_are_confidence_for_nh_parsing_cbmi = new JCheckBoxMenuItem( "Internal Node Names are Confidence Values" ) );
        _options_jmenu.add( _replace_underscores_cbmi = new JCheckBoxMenuItem( "Replace Underscores with Spaces" ) );
        _options_jmenu
                .add( _allow_errors_in_distance_to_parent_cbmi = new JCheckBoxMenuItem( "Ignore Distance Values Format Errors" ) );
        //
        _options_jmenu.add( _extract_taxonomy_no_rbmi = new JRadioButtonMenuItem( "No Taxonomy Extraction" ) );
        _options_jmenu
                .add( _extract_taxonomy_pfam_strict_rbmi = new JRadioButtonMenuItem( "Extract Taxonomy Codes/Ids from Pfam-style Node Names" ) );
        _options_jmenu
                .add( _extract_taxonomy_pfam_relaxed_rbmi = new JRadioButtonMenuItem( "Extract Taxonomy Codes/Ids from Pfam-style like Node Names" ) );
        _options_jmenu
                .add( _extract_taxonomy_agressive_rbmi = new JRadioButtonMenuItem( "Extract Taxonomy Codes/Ids/Scientific Names from Node Names" ) );
        _extract_taxonomy_pfam_strict_rbmi
                .setToolTipText( "To extract taxonomy codes/ids from node names in the form of e.g. \"BCL2_MOUSE/123-304\" or \"BCL2_10090/123-304\"" );
        _extract_taxonomy_pfam_relaxed_rbmi
                .setToolTipText( "To extract taxonomy codes/ids from node names in the form of e.g. \"bax_MOUSE\" or \"bax_10090\"" );
        _extract_taxonomy_agressive_rbmi
                .setToolTipText( "To extract taxonomy codes/ids or scientific names from node names in the form of e.g. \"MOUSE\" or \"10090\" or \"xyz_Nematostella_vectensis\"" );
        _radio_group_2 = new ButtonGroup();
        _radio_group_2.add( _extract_taxonomy_no_rbmi );
        _radio_group_2.add( _extract_taxonomy_pfam_strict_rbmi );
        _radio_group_2.add( _extract_taxonomy_pfam_relaxed_rbmi );
        _radio_group_2.add( _extract_taxonomy_agressive_rbmi );
        // 
        _options_jmenu.add( customizeMenuItemAsLabel( new JMenuItem( "Newick/Nexus Output:" ), getConfiguration() ) );
        _options_jmenu
                .add( _use_brackets_for_conf_in_nh_export_cbmi = new JCheckBoxMenuItem( USE_BRACKETS_FOR_CONF_IN_NH_LABEL ) );
        _use_brackets_for_conf_in_nh_export_cbmi
                .setToolTipText( "e.g. \"0.1[90]\" for a branch with support 90 and a length of 0.1" );
        _options_jmenu
                .add( _use_internal_names_for_conf_in_nh_export_cbmi = new JCheckBoxMenuItem( USE_INTERNAL_NAMES_FOR_CONF_IN_NH_LABEL ) );
        customizeJMenuItem( _choose_font_mi );
        customizeJMenuItem( _choose_minimal_confidence_mi );
        customizeJMenuItem( _switch_colors_mi );
        customizeJMenuItem( _print_size_mi );
        customizeJMenuItem( _choose_pdf_width_mi );
        customizeJMenuItem( _overview_placment_mi );
        customizeCheckBoxMenuItem( _show_default_node_shapes_external_cbmi, getOptions()
                .isShowDefaultNodeShapesExternal() );
        customizeCheckBoxMenuItem( _show_default_node_shapes_internal_cbmi, getOptions()
                .isShowDefaultNodeShapesInternal() );
        customizeJMenuItem( _cycle_node_shape_mi );
        customizeJMenuItem( _cycle_node_fill_mi );
        customizeJMenuItem( _choose_node_size_mi );
        customizeCheckBoxMenuItem( _color_labels_same_as_parent_branch, getOptions().isColorLabelsSameAsParentBranch() );
        customizeCheckBoxMenuItem( _color_by_taxonomic_group_cbmi, getOptions().isColorByTaxonomicGroup() );
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
        customizeCheckBoxMenuItem( _show_branch_length_values_cbmi, getOptions().isShowBranchLengthValues() );
        customizeCheckBoxMenuItem( _show_overview_cbmi, getOptions().isShowOverview() );
        customizeCheckBoxMenuItem( _label_direction_cbmi,
                                   getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL );
        customizeCheckBoxMenuItem( _antialias_print_cbmi, getOptions().isAntialiasPrint() );
        customizeCheckBoxMenuItem( _print_black_and_white_cbmi, getOptions().isPrintBlackAndWhite() );
        customizeCheckBoxMenuItem( _internal_number_are_confidence_for_nh_parsing_cbmi, getOptions()
                .isInternalNumberAreConfidenceForNhParsing() );
        customizeRadioButtonMenuItem( _extract_taxonomy_no_rbmi,
                                      getOptions().getTaxonomyExtraction() == TAXONOMY_EXTRACTION.NO );
        customizeRadioButtonMenuItem( _extract_taxonomy_pfam_strict_rbmi,
                                      getOptions().getTaxonomyExtraction() == TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
        customizeRadioButtonMenuItem( _extract_taxonomy_pfam_relaxed_rbmi,
                                      getOptions().getTaxonomyExtraction() == TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
        customizeRadioButtonMenuItem( _extract_taxonomy_agressive_rbmi,
                                      getOptions().getTaxonomyExtraction() == TAXONOMY_EXTRACTION.AGGRESSIVE );
        customizeCheckBoxMenuItem( _replace_underscores_cbmi, getOptions().isReplaceUnderscoresInNhParsing() );
        customizeCheckBoxMenuItem( _allow_errors_in_distance_to_parent_cbmi, getOptions()
                .isReplaceUnderscoresInNhParsing() );
        customizeCheckBoxMenuItem( _search_whole_words_only_cbmi, getOptions().isMatchWholeTermsOnly() );
        customizeCheckBoxMenuItem( _inverse_search_result_cbmi, getOptions().isInverseSearchResult() );
        customizeCheckBoxMenuItem( _graphics_export_visible_only_cbmi, getOptions().isGraphicsExportVisibleOnly() );
        customizeCheckBoxMenuItem( _print_using_actual_size_cbmi, getOptions().isPrintUsingActualSize() );
        customizeCheckBoxMenuItem( _graphics_export_using_actual_size_cbmi, getOptions()
                .isGraphicsExportUsingActualSize() );
        customizeCheckBoxMenuItem( _show_confidence_stddev_cbmi, getOptions().isShowConfidenceStddev() );
        customizeCheckBoxMenuItem( _use_brackets_for_conf_in_nh_export_cbmi, getOptions()
                .getNhConversionSupportValueStyle() == NH_CONVERSION_SUPPORT_VALUE_STYLE.IN_SQUARE_BRACKETS );
        customizeCheckBoxMenuItem( _use_internal_names_for_conf_in_nh_export_cbmi, getOptions()
                .getNhConversionSupportValueStyle() == NH_CONVERSION_SUPPORT_VALUE_STYLE.AS_INTERNAL_NODE_NAMES );
        _jmenubar.add( _options_jmenu );
    }

    void buildPhylogeneticInferenceMenu() {
        final InferenceManager im = getInferenceManager();
        _inference_menu = MainFrame.createMenu( "Inference", getConfiguration() );
        _inference_menu.add( _inference_from_msa_item = new JMenuItem( "From Multiple Sequence Alignment..." ) );
        customizeJMenuItem( _inference_from_msa_item );
        _inference_from_msa_item.setToolTipText( "Basic phylogenetic inference from MSA" );
        if ( im.canDoMsa() ) {
            _inference_menu.add( _inference_from_seqs_item = new JMenuItem( "From Unaligned Sequences..." ) );
            customizeJMenuItem( _inference_from_seqs_item );
            _inference_from_seqs_item
                    .setToolTipText( "Basic phylogenetic inference including multiple sequence alignment" );
        }
        else {
            _inference_menu
                    .add( _inference_from_seqs_item = new JMenuItem( "From Unaligned Sequences (no program found)" ) );
            customizeJMenuItem( _inference_from_seqs_item );
            _inference_from_seqs_item.setEnabled( false );
        }
        _jmenubar.add( _inference_menu );
    }

    void buildToolsMenu() {
        _tools_menu = createMenu( "Tools", getConfiguration() );
        _tools_menu.add( _confcolor_item = new JMenuItem( "Colorize Branches Depending on Confidence" ) );
        customizeJMenuItem( _confcolor_item );
        _tools_menu.add( _color_rank_jmi = new JMenuItem( "Colorize Subtrees via Taxonomic Rank" ) );
        customizeJMenuItem( _color_rank_jmi );
        _color_rank_jmi.setToolTipText( "for example, at \"Class\" level, colorize mammal specific subtree red" );
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
        _tools_menu.add( _annotate_item = new JMenuItem( "Annotate Sequences of Selected Nodes" ) );
        customizeJMenuItem( _annotate_item );
        _tools_menu.addSeparator();
        _tools_menu.add( _midpoint_root_item = new JMenuItem( "Midpoint-Root" ) );
        customizeJMenuItem( _midpoint_root_item );
        _tools_menu.addSeparator();
        _tools_menu.add( _collapse_species_specific_subtrees = new JMenuItem( "Collapse Species-Specific Subtrees" ) );
        customizeJMenuItem( _collapse_species_specific_subtrees );
        _tools_menu
                .add( _collapse_below_threshold = new JMenuItem( "Collapse Branches with Confidence Below Threshold into Multifurcations" ) );
        customizeJMenuItem( _collapse_below_threshold );
        _collapse_below_threshold
                .setToolTipText( "To collapse branches with confidence values below a threshold into multifurcations (in the case of multiple confidences per branch: without at least one confidence value above a threshold)" );
        _tools_menu.addSeparator();
        _tools_menu
                .add( _extract_tax_code_from_node_names_jmi = new JMenuItem( "Extract Taxonomic Data from Node Names" ) );
        customizeJMenuItem( _extract_tax_code_from_node_names_jmi );
        _extract_tax_code_from_node_names_jmi
                .setToolTipText( "To extract SwissProt/Uniprot taxonomic codes (mnemonics) from nodes names in the form of 'xyz_CAEEL', Uniprot/NCBI identifiers form of 'xyz_6239', or scientific names form of 'xyz_Caenorhabditis_elegans'" );
        _tools_menu
                .add( _move_node_names_to_tax_sn_jmi = new JMenuItem( "Transfer Node Names to Taxonomic Scientific Names" ) );
        customizeJMenuItem( _move_node_names_to_tax_sn_jmi );
        _move_node_names_to_tax_sn_jmi.setToolTipText( "To interpret node names as taxonomic scientific names" );
        _tools_menu.add( _move_node_names_to_seq_names_jmi = new JMenuItem( "Transfer Node Names to Sequence Names" ) );
        customizeJMenuItem( _move_node_names_to_seq_names_jmi );
        _move_node_names_to_seq_names_jmi.setToolTipText( "To interpret node names as sequence (protein, gene) names" );
        _tools_menu.addSeparator();
        _tools_menu.add( _obtain_seq_information_jmi = new JMenuItem( "Obtain Sequence Information" ) );
        customizeJMenuItem( _obtain_seq_information_jmi );
        _obtain_seq_information_jmi.setToolTipText( "To add additional sequence information" );
        _tools_menu
                .add( _obtain_detailed_taxonomic_information_jmi = new JMenuItem( OBTAIN_DETAILED_TAXONOMIC_INFORMATION ) );
        customizeJMenuItem( _obtain_detailed_taxonomic_information_jmi );
        _obtain_detailed_taxonomic_information_jmi
                .setToolTipText( "To add additional taxonomic information (from UniProt Taxonomy)" );
        _tools_menu
                .add( _obtain_detailed_taxonomic_information_deleting_jmi = new JMenuItem( "Obtain Detailed Taxonomic Information (deletes nodes!)" ) );
        customizeJMenuItem( _obtain_detailed_taxonomic_information_deleting_jmi );
        _obtain_detailed_taxonomic_information_deleting_jmi
                .setToolTipText( "To add additional taxonomic information, deletes nodes for which taxonomy cannot found (from UniProt Taxonomy)" );
        _tools_menu.addSeparator();
        _tools_menu.add( _read_values_jmi = new JMenuItem( "Attach Vector/Expression Values" ) );
        customizeJMenuItem( _read_values_jmi );
        _read_values_jmi.setToolTipText( "To attach vector (e.g. gene expression) values to tree nodes (beta)" );
        _jmenubar.add( _tools_menu );
        _tools_menu.add( _read_seqs_jmi = new JMenuItem( "Attach Molecular Sequences" ) );
        customizeJMenuItem( _read_seqs_jmi );
        _read_seqs_jmi
                .setToolTipText( "To attach molecular sequences to tree nodes (from Fasta-formatted file) (beta)" );
        _jmenubar.add( _tools_menu );
    }

    @Override
    void close() {
        if ( isUnsavedDataPresent() ) {
            final int r = JOptionPane.showConfirmDialog( this,
                                                         "Exit despite potentially unsaved changes?",
                                                         "Exit?",
                                                         JOptionPane.YES_NO_OPTION );
            if ( r != JOptionPane.YES_OPTION ) {
                return;
            }
        }
        exit();
    }

    void executeLineageInference() {
        if ( ( _mainpanel.getCurrentPhylogeny() == null ) || ( _mainpanel.getCurrentPhylogeny().isEmpty() ) ) {
            return;
        }
        if ( !_mainpanel.getCurrentPhylogeny().isRooted() ) {
            JOptionPane.showMessageDialog( this,
                                           "Phylogeny is not rooted.",
                                           "Cannot infer ancestral taxonomies",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        final AncestralTaxonomyInferrer inferrer = new AncestralTaxonomyInferrer( this,
                                                                                  _mainpanel.getCurrentTreePanel(),
                                                                                  _mainpanel.getCurrentPhylogeny()
                                                                                          .copy() );
        new Thread( inferrer ).start();
    }

    void exit() {
        removeAllTextFrames();
        _mainpanel.terminate();
        _contentpane.removeAll();
        setVisible( false );
        dispose();
        System.exit( 0 );
    }

    void setMsa( final Msa msa ) {
        _msa = msa;
    }

    void setMsaFile( final File msa_file ) {
        _msa_file = msa_file;
    }

    void setSeqs( final List<Sequence> seqs ) {
        _seqs = seqs;
    }

    void setSeqsFile( final File seqs_file ) {
        _seqs_file = seqs_file;
    }

    void writePhylogenyToGraphicsFile( final String file_name, final GraphicsExportType type ) {
        _mainpanel.getCurrentTreePanel().calcParametersForPainting( _mainpanel.getCurrentTreePanel().getWidth(),
                                                                    _mainpanel.getCurrentTreePanel().getHeight(),
                                                                    true );
        String file_written_to = "";
        boolean error = false;
        try {
            file_written_to = AptxUtil.writePhylogenyToGraphicsFile( file_name,
                                                                     _mainpanel.getCurrentTreePanel().getWidth(),
                                                                     _mainpanel.getCurrentTreePanel().getHeight(),
                                                                     _mainpanel.getCurrentTreePanel(),
                                                                     _mainpanel.getControlPanel(),
                                                                     type,
                                                                     getOptions() );
        }
        catch ( final IOException e ) {
            error = true;
            JOptionPane.showMessageDialog( this, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE );
        }
        if ( !error ) {
            if ( ( file_written_to != null ) && ( file_written_to.length() > 0 ) ) {
                JOptionPane.showMessageDialog( this,
                                               "Wrote image to: " + file_written_to,
                                               "Graphics Export",
                                               JOptionPane.INFORMATION_MESSAGE );
            }
            else {
                JOptionPane.showMessageDialog( this,
                                               "There was an unknown problem when attempting to write to an image file: \""
                                                       + file_name + "\"",
                                               "Error",
                                               JOptionPane.ERROR_MESSAGE );
            }
        }
        _contentpane.repaint();
    }

    private void addExpressionValuesFromFile() {
        if ( ( getCurrentTreePanel() == null ) || ( getCurrentTreePanel().getPhylogeny() == null ) ) {
            JOptionPane.showMessageDialog( this,
                                           "Need to load evolutionary tree first",
                                           "Can Not Read Expression Values",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        final File my_dir = getCurrentDir();
        if ( my_dir != null ) {
            _values_filechooser.setCurrentDirectory( my_dir );
        }
        final int result = _values_filechooser.showOpenDialog( _contentpane );
        final File file = _values_filechooser.getSelectedFile();
        if ( ( file != null ) && ( file.length() > 0 ) && ( result == JFileChooser.APPROVE_OPTION ) ) {
            BasicTable<String> t = null;
            try {
                t = BasicTableParser.parse( file, '\t' );
                if ( t.getNumberOfColumns() < 2 ) {
                    t = BasicTableParser.parse( file, ',' );
                }
                if ( t.getNumberOfColumns() < 2 ) {
                    t = BasicTableParser.parse( file, ' ' );
                }
            }
            catch ( final IOException e ) {
                JOptionPane.showMessageDialog( this,
                                               e.getMessage(),
                                               "Could Not Read Expression Value Table",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            if ( t.getNumberOfColumns() < 2 ) {
                JOptionPane.showMessageDialog( this,
                                               "Table contains " + t.getNumberOfColumns() + " column(s)",
                                               "Problem with Expression Value Table",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            if ( t.getNumberOfRows() < 1 ) {
                JOptionPane.showMessageDialog( this,
                                               "Table contains zero rows",
                                               "Problem with Expression Value Table",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ( t.getNumberOfRows() != phy.getNumberOfExternalNodes() ) {
                JOptionPane.showMessageDialog( this,
                                               "Table contains " + t.getNumberOfRows() + " rows, but tree contains "
                                                       + phy.getNumberOfExternalNodes() + " external nodes",
                                               "Warning",
                                               JOptionPane.WARNING_MESSAGE );
            }
            final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
            int not_found = 0;
            for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
                final String node_name = node.getName();
                if ( !ForesterUtil.isEmpty( node_name ) ) {
                    int row = -1;
                    try {
                        row = t.findRow( node_name );
                    }
                    catch ( final IllegalArgumentException e ) {
                        JOptionPane
                                .showMessageDialog( this,
                                                    e.getMessage(),
                                                    "Error Mapping Node Identifiers to Expression Value Identifiers",
                                                    JOptionPane.ERROR_MESSAGE );
                        return;
                    }
                    if ( row < 0 ) {
                        if ( node.isExternal() ) {
                            not_found++;
                        }
                        continue;
                    }
                    final List<Double> l = new ArrayList<Double>();
                    for( int col = 1; col < t.getNumberOfColumns(); ++col ) {
                        double d = -100;
                        try {
                            d = Double.parseDouble( t.getValueAsString( col, row ) );
                        }
                        catch ( final NumberFormatException e ) {
                            JOptionPane.showMessageDialog( this,
                                                           "Could not parse \"" + t.getValueAsString( col, row )
                                                                   + "\" into a decimal value",
                                                           "Issue with Expression Value Table",
                                                           JOptionPane.ERROR_MESSAGE );
                            return;
                        }
                        stats.addValue( d );
                        l.add( d );
                    }
                    if ( !l.isEmpty() ) {
                        if ( node.getNodeData().getProperties() != null ) {
                            node.getNodeData().getProperties()
                                    .removePropertiesWithGivenReferencePrefix( PhyloXmlUtil.VECTOR_PROPERTY_REF );
                        }
                        node.getNodeData().setVector( l );
                    }
                }
            }
            if ( not_found > 0 ) {
                JOptionPane.showMessageDialog( this, "Could not fine expression values for " + not_found
                        + " external node(s)", "Warning", JOptionPane.WARNING_MESSAGE );
            }
            getCurrentTreePanel().setStatisticsForExpressionValues( stats );
        }
    }

    private void addSequencesFromFile() {
        if ( ( getCurrentTreePanel() == null ) || ( getCurrentTreePanel().getPhylogeny() == null ) ) {
            JOptionPane.showMessageDialog( this,
                                           "Need to load evolutionary tree first",
                                           "Can Not Read Sequences",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        final File my_dir = getCurrentDir();
        if ( my_dir != null ) {
            _sequences_filechooser.setCurrentDirectory( my_dir );
        }
        final int result = _sequences_filechooser.showOpenDialog( _contentpane );
        final File file = _sequences_filechooser.getSelectedFile();
        List<Sequence> seqs = null;
        if ( ( file != null ) && !file.isDirectory() && ( result == JFileChooser.APPROVE_OPTION ) ) {
            try {
                if ( FastaParser.isLikelyFasta( new FileInputStream( file ) ) ) {
                    seqs = FastaParser.parse( new FileInputStream( file ) );
                }
                else {
                    JOptionPane.showMessageDialog( this,
                                                   "Format does not appear to be Fasta",
                                                   "Multiple sequence file format error",
                                                   JOptionPane.ERROR_MESSAGE );
                    return;
                }
            }
            catch ( final MsaFormatException e ) {
                setArrowCursor();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Multiple sequence file format error",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            catch ( final IOException e ) {
                setArrowCursor();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Failed to read multiple sequence file",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            catch ( final Exception e ) {
                setArrowCursor();
                e.printStackTrace();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Unexpected error during reading of multiple sequence file",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            if ( ( seqs == null ) || ( seqs.size() < 1 ) ) {
                JOptionPane.showMessageDialog( this,
                                               "Multiple sequence file is empty",
                                               "Empty multiple sequence file",
                                               JOptionPane.ERROR_MESSAGE );
                setArrowCursor();
                return;
            }
        }
        if ( seqs != null ) {
            for( final Sequence seq : seqs ) {
                System.out.println( seq.getIdentifier() );
            }
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            int total_counter = 0;
            int attached_counter = 0;
            for( final Sequence seq : seqs ) {
                ++total_counter;
                final String seq_name = seq.getIdentifier();
                if ( !ForesterUtil.isEmpty( seq_name ) ) {
                    List<PhylogenyNode> nodes = phy.getNodesViaSequenceName( seq_name );
                    if ( nodes.isEmpty() ) {
                        nodes = phy.getNodesViaSequenceSymbol( seq_name );
                    }
                    if ( nodes.isEmpty() ) {
                        nodes = phy.getNodesViaGeneName( seq_name );
                    }
                    if ( nodes.isEmpty() ) {
                        nodes = phy.getNodes( seq_name );
                    }
                    if ( nodes.size() > 1 ) {
                        JOptionPane.showMessageDialog( this,
                                                       "Sequence name \"" + seq_name + "\" is not unique",
                                                       "Sequence name not unique",
                                                       JOptionPane.ERROR_MESSAGE );
                        setArrowCursor();
                        return;
                    }
                    final String[] a = seq_name.split( "\\s" );
                    if ( nodes.isEmpty() && ( a.length > 1 ) ) {
                        final String seq_name_split = a[ 0 ];
                        nodes = phy.getNodesViaSequenceName( seq_name_split );
                        if ( nodes.isEmpty() ) {
                            nodes = phy.getNodesViaSequenceSymbol( seq_name_split );
                        }
                        if ( nodes.isEmpty() ) {
                            nodes = phy.getNodes( seq_name_split );
                        }
                        if ( nodes.size() > 1 ) {
                            JOptionPane.showMessageDialog( this, "Split sequence name \"" + seq_name_split
                                    + "\" is not unique", "Sequence name not unique", JOptionPane.ERROR_MESSAGE );
                            setArrowCursor();
                            return;
                        }
                    }
                    if ( nodes.size() == 1 ) {
                        ++attached_counter;
                        final PhylogenyNode n = nodes.get( 0 );
                        if ( !n.getNodeData().isHasSequence() ) {
                            n.getNodeData().addSequence( new org.forester.phylogeny.data.Sequence() );
                        }
                        n.getNodeData().getSequence().setMolecularSequence( seq.getMolecularSequenceAsString() );
                        if ( ForesterUtil.isEmpty( n.getNodeData().getSequence().getName() ) ) {
                            n.getNodeData().getSequence().setName( seq_name );
                        }
                    }
                }
            }
            if ( attached_counter > 0 ) {
                int ext_nodes = 0;
                int ext_nodes_with_seq = 0;
                for( final PhylogenyNodeIterator iter = phy.iteratorExternalForward(); iter.hasNext(); ) {
                    ++ext_nodes;
                    final PhylogenyNode n = iter.next();
                    if ( n.getNodeData().isHasSequence()
                            && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getMolecularSequence() ) ) {
                        ++ext_nodes_with_seq;
                    }
                }
                final String s;
                if ( ext_nodes == ext_nodes_with_seq ) {
                    s = "All " + ext_nodes_with_seq + " external nodes now have a molecular sequence attached to them.";
                }
                else {
                    s = ext_nodes_with_seq + " out of " + ext_nodes
                            + " external nodes now have a molecular sequence attached to them.";
                }
                if ( ( attached_counter == total_counter ) && ( ext_nodes == ext_nodes_with_seq ) ) {
                    JOptionPane.showMessageDialog( this,
                                                   "Attached all " + total_counter + " sequences to tree nodes.\n" + s,
                                                   "All sequences attached",
                                                   JOptionPane.INFORMATION_MESSAGE );
                }
                else {
                    JOptionPane.showMessageDialog( this, "Attached " + attached_counter
                            + " sequences out of a total of " + total_counter + " sequences.\n" + s, attached_counter
                            + " sequences attached", JOptionPane.WARNING_MESSAGE );
                }
            }
            else {
                JOptionPane.showMessageDialog( this, "No maching tree node for any of the " + total_counter
                        + " sequences", "Could not attach any sequences", JOptionPane.ERROR_MESSAGE );
            }
        }
    }

    private void choosePdfWidth() {
        final String s = ( String ) JOptionPane.showInputDialog( this,
                                                                 "Please enter the default line width for PDF export.\n"
                                                                         + "[current value: "
                                                                         + getOptions().getPrintLineWidth() + "]\n",
                                                                 "Line Width for PDF Export",
                                                                 JOptionPane.QUESTION_MESSAGE,
                                                                 null,
                                                                 null,
                                                                 getOptions().getPrintLineWidth() );
        if ( !ForesterUtil.isEmpty( s ) ) {
            boolean success = true;
            float f = 0.0f;
            final String m_str = s.trim();
            if ( !ForesterUtil.isEmpty( m_str ) ) {
                try {
                    f = Float.parseFloat( m_str );
                }
                catch ( final Exception ex ) {
                    success = false;
                }
            }
            else {
                success = false;
            }
            if ( success && ( f > 0.0 ) ) {
                getOptions().setPrintLineWidth( f );
            }
        }
    }

    private void choosePrintSize() {
        final String s = ( String ) JOptionPane.showInputDialog( this,
                                                                 "Please enter values for width and height,\nseparated by a comma.\n"
                                                                         + "[current values: "
                                                                         + getOptions().getPrintSizeX() + ", "
                                                                         + getOptions().getPrintSizeY() + "]\n"
                                                                         + "[A4: " + Constants.A4_SIZE_X + ", "
                                                                         + Constants.A4_SIZE_Y + "]\n" + "[US Letter: "
                                                                         + Constants.US_LETTER_SIZE_X + ", "
                                                                         + Constants.US_LETTER_SIZE_Y + "]",
                                                                 "Default Size for Graphics Export",
                                                                 JOptionPane.QUESTION_MESSAGE,
                                                                 null,
                                                                 null,
                                                                 getOptions().getPrintSizeX() + ", "
                                                                         + getOptions().getPrintSizeY() );
        if ( !ForesterUtil.isEmpty( s ) && ( s.indexOf( ',' ) > 0 ) ) {
            boolean success = true;
            int x = 0;
            int y = 0;
            final String[] str_ary = s.split( "," );
            if ( str_ary.length == 2 ) {
                final String x_str = str_ary[ 0 ].trim();
                final String y_str = str_ary[ 1 ].trim();
                if ( !ForesterUtil.isEmpty( x_str ) && !ForesterUtil.isEmpty( y_str ) ) {
                    try {
                        x = Integer.parseInt( x_str );
                        y = Integer.parseInt( y_str );
                    }
                    catch ( final Exception ex ) {
                        success = false;
                    }
                }
                else {
                    success = false;
                }
            }
            else {
                success = false;
            }
            if ( success && ( x > 1 ) && ( y > 1 ) ) {
                getOptions().setPrintSizeX( x );
                getOptions().setPrintSizeY( y );
            }
        }
    }

    private void closeCurrentPane() {
        if ( getMainPanel().getCurrentTreePanel() != null ) {
            if ( getMainPanel().getCurrentTreePanel().isEdited() ) {
                final int r = JOptionPane.showConfirmDialog( this,
                                                             "Close tab despite potentially unsaved changes?",
                                                             "Close Tab?",
                                                             JOptionPane.YES_NO_OPTION );
                if ( r != JOptionPane.YES_OPTION ) {
                    return;
                }
            }
            getMainPanel().closeCurrentPane();
            activateSaveAllIfNeeded();
        }
    }

    private void collapse( final Phylogeny phy, final double m ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        final List<PhylogenyNode> to_be_removed = new ArrayList<PhylogenyNode>();
        double min_support = Double.MAX_VALUE;
        boolean conf_present = false;
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.isRoot() ) {
                final List<Confidence> c = n.getBranchData().getConfidences();
                if ( ( c != null ) && ( c.size() > 0 ) ) {
                    conf_present = true;
                    double max = 0;
                    for( final Confidence confidence : c ) {
                        if ( confidence.getValue() > max ) {
                            max = confidence.getValue();
                        }
                    }
                    if ( max < getMinNotCollapseConfidenceValue() ) {
                        to_be_removed.add( n );
                    }
                    if ( max < min_support ) {
                        min_support = max;
                    }
                }
            }
        }
        if ( conf_present ) {
            for( final PhylogenyNode node : to_be_removed ) {
                PhylogenyMethods.removeNode( node, phy );
            }
            if ( to_be_removed.size() > 0 ) {
                phy.externalNodesHaveChanged();
                phy.clearHashIdToNodeMap();
                phy.recalculateNumberOfExternalDescendants( true );
                getCurrentTreePanel().resetNodeIdToDistToLeafMap();
                getCurrentTreePanel().updateSetOfCollapsedExternalNodes();
                getCurrentTreePanel().calculateLongestExtNodeInfo();
                getCurrentTreePanel().setNodeInPreorderToNull();
                getCurrentTreePanel().recalculateMaxDistanceToRoot();
                getCurrentTreePanel().resetPreferredSize();
                getCurrentTreePanel().setEdited( true );
                getCurrentTreePanel().repaint();
                repaint();
            }
            if ( to_be_removed.size() > 0 ) {
                JOptionPane.showMessageDialog( this, "Collapsed " + to_be_removed.size()
                        + " branches with\nconfidence values below " + getMinNotCollapseConfidenceValue(), "Collapsed "
                        + to_be_removed.size() + " branches", JOptionPane.INFORMATION_MESSAGE );
            }
            else {
                JOptionPane.showMessageDialog( this, "No branch collapsed,\nminimum confidence value per branch is "
                        + min_support, "No branch collapsed", JOptionPane.INFORMATION_MESSAGE );
            }
        }
        else {
            JOptionPane.showMessageDialog( this,
                                           "No branch collapsed because no confidence values present",
                                           "No confidence values present",
                                           JOptionPane.INFORMATION_MESSAGE );
        }
    }

    private void collapseBelowThreshold() {
        if ( getCurrentTreePanel() != null ) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ( ( phy != null ) && !phy.isEmpty() ) {
                final String s = ( String ) JOptionPane.showInputDialog( this,
                                                                         "Please enter the minimum confidence value\n",
                                                                         "Minimal Confidence Value",
                                                                         JOptionPane.QUESTION_MESSAGE,
                                                                         null,
                                                                         null,
                                                                         getMinNotCollapseConfidenceValue() );
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
                        setMinNotCollapseConfidenceValue( m );
                        collapse( phy, m );
                    }
                }
            }
        }
    }

    private PhyloXmlParser createPhyloXmlParser() {
        PhyloXmlParser xml_parser = null;
        if ( getConfiguration().isValidatePhyloXmlAgainstSchema() ) {
            try {
                xml_parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
            }
            catch ( final Exception e ) {
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "failed to create validating XML parser",
                                               JOptionPane.WARNING_MESSAGE );
            }
        }
        if ( xml_parser == null ) {
            xml_parser = PhyloXmlParser.createPhyloXmlParser();
        }
        return xml_parser;
    }

    private void executePhyleneticInference( final boolean from_unaligned_seqs ) {
        final PhyloInferenceDialog dialog = new PhyloInferenceDialog( this,
                                                                      getPhylogeneticInferenceOptions(),
                                                                      from_unaligned_seqs );
        dialog.activate();
        if ( dialog.getValue() == JOptionPane.OK_OPTION ) {
            if ( !from_unaligned_seqs ) {
                if ( getMsa() != null ) {
                    final PhylogeneticInferrer inferrer = new PhylogeneticInferrer( getMsa(),
                                                                                    getPhylogeneticInferenceOptions()
                                                                                            .copy(), this );
                    new Thread( inferrer ).start();
                }
                else {
                    JOptionPane.showMessageDialog( this,
                                                   "No multiple sequence alignment selected",
                                                   "Phylogenetic Inference Not Launched",
                                                   JOptionPane.WARNING_MESSAGE );
                }
            }
            else {
                if ( getSeqs() != null ) {
                    final PhylogeneticInferrer inferrer = new PhylogeneticInferrer( getSeqs(),
                                                                                    getPhylogeneticInferenceOptions()
                                                                                            .copy(), this );
                    new Thread( inferrer ).start();
                }
                else {
                    JOptionPane.showMessageDialog( this,
                                                   "No input sequences selected",
                                                   "Phylogenetic Inference Not Launched",
                                                   JOptionPane.WARNING_MESSAGE );
                }
            }
        }
    }

    private void extractTaxDataFromNodeNames() throws PhyloXmlDataFormatException {
        final StringBuilder sb = new StringBuilder();
        final StringBuilder sb_failed = new StringBuilder();
        int counter = 0;
        int counter_failed = 0;
        if ( getCurrentTreePanel() != null ) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ( ( phy != null ) && !phy.isEmpty() ) {
                final PhylogenyNodeIterator it = phy.iteratorExternalForward();
                while ( it.hasNext() ) {
                    final PhylogenyNode n = it.next();
                    final String name = n.getName().trim();
                    if ( !ForesterUtil.isEmpty( name ) ) {
                        final String nt = ParserUtils.extractTaxonomyDataFromNodeName( n,
                                                                                       TAXONOMY_EXTRACTION.AGGRESSIVE );
                        if ( !ForesterUtil.isEmpty( nt ) ) {
                            if ( counter < 15 ) {
                                sb.append( name + ": " + nt + "\n" );
                            }
                            else if ( counter == 15 ) {
                                sb.append( "...\n" );
                            }
                            counter++;
                        }
                        else {
                            if ( counter_failed < 15 ) {
                                sb_failed.append( name + "\n" );
                            }
                            else if ( counter_failed == 15 ) {
                                sb_failed.append( "...\n" );
                            }
                            counter_failed++;
                        }
                    }
                }
                if ( counter > 0 ) {
                    String failed = "";
                    String all = "all ";
                    if ( counter_failed > 0 ) {
                        all = "";
                        failed = "\nCould not extract taxonomic data for " + counter_failed
                                + " named external nodes:\n" + sb_failed;
                    }
                    JOptionPane.showMessageDialog( this,
                                                   "Extracted taxonomic data from " + all + counter
                                                           + " named external nodes:\n" + sb.toString() + failed,
                                                   "Taxonomic Data Extraction Completed",
                                                   counter_failed > 0 ? JOptionPane.WARNING_MESSAGE
                                                           : JOptionPane.INFORMATION_MESSAGE );
                }
                else {
                    JOptionPane
                            .showMessageDialog( this,
                                                "Could not extract any taxonomic data.\nMaybe node names are empty\n"
                                                        + "or not in the forms \"XYZ_CAEEL\", \"XYZ_6239\", or \"XYZ_Caenorhabditis_elegans\"\n"
                                                        + "or nodes already have taxonomic data?\n",
                                                "No Taxonomic Data Extracted",
                                                JOptionPane.ERROR_MESSAGE );
                }
            }
        }
    }

    private ControlPanel getControlPanel() {
        return getMainPanel().getControlPanel();
    }

    private File getCurrentDir() {
        if ( ( _current_dir == null ) || !_current_dir.canRead() ) {
            if ( ForesterUtil.isWindows() ) {
                try {
                    _current_dir = new File( WindowsUtils.getCurrentUserDesktopPath() );
                }
                catch ( final Exception e ) {
                    _current_dir = null;
                }
            }
        }
        if ( ( _current_dir == null ) || !_current_dir.canRead() ) {
            if ( System.getProperty( "user.home" ) != null ) {
                _current_dir = new File( System.getProperty( "user.home" ) );
            }
            else if ( System.getProperty( "user.dir" ) != null ) {
                _current_dir = new File( System.getProperty( "user.dir" ) );
            }
        }
        return _current_dir;
    }

    private double getMinNotCollapseConfidenceValue() {
        return _min_not_collapse;
    }

    private PhylogeneticInferenceOptions getPhylogeneticInferenceOptions() {
        if ( _phylogenetic_inference_options == null ) {
            _phylogenetic_inference_options = new PhylogeneticInferenceOptions();
        }
        return _phylogenetic_inference_options;
    }

    private boolean isUnsavedDataPresent() {
        final List<TreePanel> tps = getMainPanel().getTreePanels();
        for( final TreePanel tp : tps ) {
            if ( tp.isEdited() ) {
                return true;
            }
        }
        return false;
    }

    private void moveNodeNamesToSeqNames() throws PhyloXmlDataFormatException {
        if ( getCurrentTreePanel() != null ) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ( ( phy != null ) && !phy.isEmpty() ) {
                PhylogenyMethods
                        .transferNodeNameToField( phy, PhylogenyMethods.PhylogenyNodeField.SEQUENCE_NAME, false );
            }
        }
    }

    private void moveNodeNamesToTaxSn() throws PhyloXmlDataFormatException {
        if ( getCurrentTreePanel() != null ) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ( ( phy != null ) && !phy.isEmpty() ) {
                PhylogenyMethods.transferNodeNameToField( phy,
                                                          PhylogenyMethods.PhylogenyNodeField.TAXONOMY_SCIENTIFIC_NAME,
                                                          false );
            }
        }
    }

    private void newTree() {
        final Phylogeny[] phys = new Phylogeny[ 1 ];
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode node = new PhylogenyNode();
        phy.setRoot( node );
        phy.setRooted( true );
        phys[ 0 ] = phy;
        AptxUtil.addPhylogeniesToTabs( phys, "", "", getConfiguration(), getMainPanel() );
        _mainpanel.getControlPanel().showWhole();
        _mainpanel.getCurrentTreePanel().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        _mainpanel.getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        if ( getMainPanel().getMainFrame() == null ) {
            // Must be "E" applet version.
            ( ( ArchaeopteryxE ) ( ( MainPanelApplets ) getMainPanel() ).getApplet() )
                    .setSelectedTypeInTypeMenu( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        }
        else {
            getMainPanel().getMainFrame().setSelectedTypeInTypeMenu( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        }
        activateSaveAllIfNeeded();
        System.gc();
    }

    private void obtainDetailedTaxonomicInformation() {
        if ( getCurrentTreePanel() != null ) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ( ( phy != null ) && !phy.isEmpty() ) {
                final TaxonomyDataManager t = new TaxonomyDataManager( this,
                                                                       _mainpanel.getCurrentTreePanel(),
                                                                       phy.copy(),
                                                                       false,
                                                                       true );
                new Thread( t ).start();
            }
        }
    }

    private void obtainDetailedTaxonomicInformationDelete() {
        if ( getCurrentTreePanel() != null ) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ( ( phy != null ) && !phy.isEmpty() ) {
                final TaxonomyDataManager t = new TaxonomyDataManager( this,
                                                                       _mainpanel.getCurrentTreePanel(),
                                                                       phy.copy(),
                                                                       true,
                                                                       true );
                new Thread( t ).start();
            }
        }
    }

    private void obtainSequenceInformation() {
        if ( getCurrentTreePanel() != null ) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ( ( phy != null ) && !phy.isEmpty() ) {
                final SequenceDataRetriver u = new SequenceDataRetriver( this,
                                                                         _mainpanel.getCurrentTreePanel(),
                                                                         phy.copy() );
                new Thread( u ).start();
            }
        }
    }

    private void print() {
        if ( ( getCurrentTreePanel() == null ) || ( getCurrentTreePanel().getPhylogeny() == null )
                || getCurrentTreePanel().getPhylogeny().isEmpty() ) {
            return;
        }
        if ( !getOptions().isPrintUsingActualSize() ) {
            getCurrentTreePanel().calcParametersForPainting( getOptions().getPrintSizeX() - 80,
                                                             getOptions().getPrintSizeY() - 140,
                                                             true );
            getCurrentTreePanel().resetPreferredSize();
            getCurrentTreePanel().repaint();
        }
        final String job_name = Constants.PRG_NAME;
        boolean error = false;
        String printer_name = null;
        try {
            printer_name = Printer.print( getCurrentTreePanel(), job_name );
        }
        catch ( final Exception e ) {
            error = true;
            JOptionPane.showMessageDialog( this, e.getMessage(), "Printing Error", JOptionPane.ERROR_MESSAGE );
        }
        if ( !error && ( printer_name != null ) ) {
            String msg = "Printing data sent to printer";
            if ( printer_name.length() > 1 ) {
                msg += " [" + printer_name + "]";
            }
            JOptionPane.showMessageDialog( this, msg, "Printing...", JOptionPane.INFORMATION_MESSAGE );
        }
        if ( !getOptions().isPrintUsingActualSize() ) {
            getControlPanel().showWhole();
        }
    }

    private void printPhylogenyToPdf( final String file_name ) {
        if ( !getOptions().isPrintUsingActualSize() ) {
            getCurrentTreePanel().calcParametersForPainting( getOptions().getPrintSizeX(),
                                                             getOptions().getPrintSizeY(),
                                                             true );
            getCurrentTreePanel().resetPreferredSize();
            getCurrentTreePanel().repaint();
        }
        String pdf_written_to = "";
        boolean error = false;
        try {
            if ( getOptions().isPrintUsingActualSize() ) {
                pdf_written_to = PdfExporter.writePhylogenyToPdf( file_name,
                                                                  getCurrentTreePanel(),
                                                                  getCurrentTreePanel().getWidth(),
                                                                  getCurrentTreePanel().getHeight() );
            }
            else {
                pdf_written_to = PdfExporter.writePhylogenyToPdf( file_name, getCurrentTreePanel(), getOptions()
                        .getPrintSizeX(), getOptions().getPrintSizeY() );
            }
        }
        catch ( final IOException e ) {
            error = true;
            JOptionPane.showMessageDialog( this, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE );
        }
        if ( !error ) {
            if ( !ForesterUtil.isEmpty( pdf_written_to ) ) {
                JOptionPane.showMessageDialog( this,
                                               "Wrote PDF to: " + pdf_written_to,
                                               "Information",
                                               JOptionPane.INFORMATION_MESSAGE );
            }
            else {
                JOptionPane.showMessageDialog( this,
                                               "There was an unknown problem when attempting to write to PDF file: \""
                                                       + file_name + "\"",
                                               "Error",
                                               JOptionPane.ERROR_MESSAGE );
            }
        }
        if ( !getOptions().isPrintUsingActualSize() ) {
            getControlPanel().showWhole();
        }
    }

    private void readPhylogeniesFromFile() {
        boolean exception = false;
        Phylogeny[] phys = null;
        // Set an initial directory if none set yet
        final File my_dir = getCurrentDir();
        _open_filechooser.setMultiSelectionEnabled( true );
        // Open file-open dialog and set current directory
        if ( my_dir != null ) {
            _open_filechooser.setCurrentDirectory( my_dir );
        }
        final int result = _open_filechooser.showOpenDialog( _contentpane );
        // All done: get the file
        final File[] files = _open_filechooser.getSelectedFiles();
        setCurrentDir( _open_filechooser.getCurrentDirectory() );
        boolean nhx_or_nexus = false;
        if ( ( files != null ) && ( files.length > 0 ) && ( result == JFileChooser.APPROVE_OPTION ) ) {
            for( final File file : files ) {
                if ( ( file != null ) && !file.isDirectory() ) {
                    if ( _mainpanel.getCurrentTreePanel() != null ) {
                        _mainpanel.getCurrentTreePanel().setWaitCursor();
                    }
                    else {
                        _mainpanel.setWaitCursor();
                    }
                    if ( ( _open_filechooser.getFileFilter() == MainFrameApplication.nhfilter )
                            || ( _open_filechooser.getFileFilter() == MainFrameApplication.nhxfilter ) ) {
                        try {
                            final NHXParser nhx = new NHXParser();
                            setSpecialOptionsForNhxParser( nhx );
                            phys = PhylogenyMethods.readPhylogenies( nhx, file );
                            nhx_or_nexus = true;
                        }
                        catch ( final Exception e ) {
                            exception = true;
                            exceptionOccuredDuringOpenFile( e );
                        }
                    }
                    else if ( _open_filechooser.getFileFilter() == MainFrameApplication.xmlfilter ) {
                        warnIfNotPhyloXmlValidation( getConfiguration() );
                        try {
                            final PhyloXmlParser xml_parser = createPhyloXmlParser();
                            phys = PhylogenyMethods.readPhylogenies( xml_parser, file );
                        }
                        catch ( final Exception e ) {
                            exception = true;
                            exceptionOccuredDuringOpenFile( e );
                        }
                    }
                    else if ( _open_filechooser.getFileFilter() == MainFrameApplication.tolfilter ) {
                        try {
                            phys = PhylogenyMethods.readPhylogenies( new TolParser(), file );
                        }
                        catch ( final Exception e ) {
                            exception = true;
                            exceptionOccuredDuringOpenFile( e );
                        }
                    }
                    else if ( _open_filechooser.getFileFilter() == MainFrameApplication.nexusfilter ) {
                        try {
                            final NexusPhylogeniesParser nex = new NexusPhylogeniesParser();
                            setSpecialOptionsForNexParser( nex );
                            phys = PhylogenyMethods.readPhylogenies( nex, file );
                            nhx_or_nexus = true;
                        }
                        catch ( final Exception e ) {
                            exception = true;
                            exceptionOccuredDuringOpenFile( e );
                        }
                    }
                    // "*.*":
                    else {
                        try {
                            final PhylogenyParser parser = ParserUtils
                                    .createParserDependingOnFileType( file, getConfiguration()
                                            .isValidatePhyloXmlAgainstSchema() );
                            if ( parser instanceof NexusPhylogeniesParser ) {
                                final NexusPhylogeniesParser nex = ( NexusPhylogeniesParser ) parser;
                                setSpecialOptionsForNexParser( nex );
                                nhx_or_nexus = true;
                            }
                            else if ( parser instanceof NHXParser ) {
                                final NHXParser nhx = ( NHXParser ) parser;
                                setSpecialOptionsForNhxParser( nhx );
                                nhx_or_nexus = true;
                            }
                            else if ( parser instanceof PhyloXmlParser ) {
                                warnIfNotPhyloXmlValidation( getConfiguration() );
                            }
                            phys = PhylogenyMethods.readPhylogenies( parser, file );
                        }
                        catch ( final Exception e ) {
                            exception = true;
                            exceptionOccuredDuringOpenFile( e );
                        }
                    }
                    if ( _mainpanel.getCurrentTreePanel() != null ) {
                        _mainpanel.getCurrentTreePanel().setArrowCursor();
                    }
                    else {
                        _mainpanel.setArrowCursor();
                    }
                    if ( !exception && ( phys != null ) && ( phys.length > 0 ) ) {
                        boolean one_desc = false;
                        if ( nhx_or_nexus ) {
                            for( final Phylogeny phy : phys ) {
                                if ( getOptions().isInternalNumberAreConfidenceForNhParsing() ) {
                                    PhylogenyMethods.transferInternalNodeNamesToConfidence( phy, "" );
                                }
                                if ( PhylogenyMethods.getMinimumDescendentsPerInternalNodes( phy ) == 1 ) {
                                    one_desc = true;
                                    break;
                                }
                            }
                        }
                        AptxUtil.addPhylogeniesToTabs( phys,
                                                       file.getName(),
                                                       file.getAbsolutePath(),
                                                       getConfiguration(),
                                                       getMainPanel() );
                        _mainpanel.getControlPanel().showWhole();
                        if ( nhx_or_nexus && one_desc ) {
                            JOptionPane
                                    .showMessageDialog( this,
                                                        "One or more trees contain (a) node(s) with one descendant, "
                                                                + ForesterUtil.LINE_SEPARATOR
                                                                + "possibly indicating illegal parentheses within node names.",
                                                        "Warning: Possible Error in New Hampshire Formatted Data",
                                                        JOptionPane.WARNING_MESSAGE );
                        }
                    }
                }
            }
        }
        activateSaveAllIfNeeded();
        System.gc();
    }

    private void readSpeciesTreeFromFile() {
        Phylogeny t = null;
        boolean exception = false;
        final File my_dir = getCurrentDir();
        _open_filechooser_for_species_tree.setSelectedFile( new File( "" ) );
        if ( my_dir != null ) {
            _open_filechooser_for_species_tree.setCurrentDirectory( my_dir );
        }
        final int result = _open_filechooser_for_species_tree.showOpenDialog( _contentpane );
        final File file = _open_filechooser_for_species_tree.getSelectedFile();
        if ( ( file != null ) && ( result == JFileChooser.APPROVE_OPTION ) ) {
            if ( _open_filechooser_for_species_tree.getFileFilter() == MainFrameApplication.xmlfilter ) {
                try {
                    final Phylogeny[] trees = PhylogenyMethods.readPhylogenies( PhyloXmlParser
                            .createPhyloXmlParserXsdValidating(), file );
                    t = trees[ 0 ];
                }
                catch ( final Exception e ) {
                    exception = true;
                    exceptionOccuredDuringOpenFile( e );
                }
            }
            else if ( _open_filechooser_for_species_tree.getFileFilter() == MainFrameApplication.tolfilter ) {
                try {
                    final Phylogeny[] trees = PhylogenyMethods.readPhylogenies( new TolParser(), file );
                    t = trees[ 0 ];
                }
                catch ( final Exception e ) {
                    exception = true;
                    exceptionOccuredDuringOpenFile( e );
                }
            }
            // "*.*":
            else {
                try {
                    final Phylogeny[] trees = PhylogenyMethods.readPhylogenies( PhyloXmlParser
                            .createPhyloXmlParserXsdValidating(), file );
                    t = trees[ 0 ];
                }
                catch ( final Exception e ) {
                    exception = true;
                    exceptionOccuredDuringOpenFile( e );
                }
            }
            if ( !exception && ( t != null ) && !t.isRooted() ) {
                exception = true;
                t = null;
                JOptionPane.showMessageDialog( this,
                                               "Species tree is not rooted",
                                               "Species tree not loaded",
                                               JOptionPane.ERROR_MESSAGE );
            }
            if ( !exception && ( t != null ) ) {
                final Set<Taxonomy> tax_set = new HashSet<Taxonomy>();
                for( final PhylogenyNodeIterator it = t.iteratorExternalForward(); it.hasNext(); ) {
                    final PhylogenyNode node = it.next();
                    if ( !node.getNodeData().isHasTaxonomy() ) {
                        exception = true;
                        t = null;
                        JOptionPane
                                .showMessageDialog( this,
                                                    "Species tree contains external node(s) without taxonomy information",
                                                    "Species tree not loaded",
                                                    JOptionPane.ERROR_MESSAGE );
                        break;
                    }
                    else {
                        if ( tax_set.contains( node.getNodeData().getTaxonomy() ) ) {
                            exception = true;
                            t = null;
                            JOptionPane.showMessageDialog( this,
                                                           "Taxonomy ["
                                                                   + node.getNodeData().getTaxonomy().asSimpleText()
                                                                   + "] is not unique in species tree",
                                                           "Species tree not loaded",
                                                           JOptionPane.ERROR_MESSAGE );
                            break;
                        }
                        else {
                            tax_set.add( node.getNodeData().getTaxonomy() );
                        }
                    }
                }
            }
            if ( !exception && ( t != null ) ) {
                setSpeciesTree( t );
                JOptionPane.showMessageDialog( this,
                                               "Species tree successfully loaded",
                                               "Species tree loaded",
                                               JOptionPane.INFORMATION_MESSAGE );
            }
            _contentpane.repaint();
            System.gc();
        }
    }

    private void setArrowCursor() {
        try {
            _mainpanel.getCurrentTreePanel().setArrowCursor();
        }
        catch ( final Exception ex ) {
            // Do nothing.
        }
    }

    private void setCurrentDir( final File current_dir ) {
        _current_dir = current_dir;
    }

    private void setMinNotCollapseConfidenceValue( final double min_not_collapse ) {
        _min_not_collapse = min_not_collapse;
    }

    private void setPhylogeneticInferenceOptions( final PhylogeneticInferenceOptions phylogenetic_inference_options ) {
        _phylogenetic_inference_options = phylogenetic_inference_options;
    }

    private void setSpecialOptionsForNexParser( final NexusPhylogeniesParser nex ) {
        nex.setReplaceUnderscores( getOptions().isReplaceUnderscoresInNhParsing() );
        nex.setTaxonomyExtraction( getOptions().getTaxonomyExtraction() );
    }

    private void setSpecialOptionsForNhxParser( final NHXParser nhx ) {
        nhx.setReplaceUnderscores( getOptions().isReplaceUnderscoresInNhParsing() );
        nhx.setTaxonomyExtraction( getOptions().getTaxonomyExtraction() );
        nhx.setAllowErrorsInDistanceToParent( getOptions().isAllowErrorsInDistanceToParent() );
    }

    private void writeAllToFile() {
        if ( ( getMainPanel().getTabbedPane() == null ) || ( getMainPanel().getTabbedPane().getTabCount() < 1 ) ) {
            return;
        }
        final File my_dir = getCurrentDir();
        if ( my_dir != null ) {
            _save_filechooser.setCurrentDirectory( my_dir );
        }
        _save_filechooser.setSelectedFile( new File( "" ) );
        final int result = _save_filechooser.showSaveDialog( _contentpane );
        final File file = _save_filechooser.getSelectedFile();
        setCurrentDir( _save_filechooser.getCurrentDirectory() );
        if ( ( file != null ) && ( result == JFileChooser.APPROVE_OPTION ) ) {
            if ( file.exists() ) {
                final int i = JOptionPane.showConfirmDialog( this,
                                                             file + " already exists. Overwrite?",
                                                             "Warning",
                                                             JOptionPane.OK_CANCEL_OPTION,
                                                             JOptionPane.WARNING_MESSAGE );
                if ( i != JOptionPane.OK_OPTION ) {
                    return;
                }
                else {
                    try {
                        file.delete();
                    }
                    catch ( final Exception e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "Failed to delete: " + file,
                                                       "Error",
                                                       JOptionPane.WARNING_MESSAGE );
                    }
                }
            }
            final int count = getMainPanel().getTabbedPane().getTabCount();
            final List<Phylogeny> trees = new ArrayList<Phylogeny>();
            for( int i = 0; i < count; ++i ) {
                final Phylogeny phy = getMainPanel().getPhylogeny( i );
                if ( ForesterUtil.isEmpty( phy.getName() )
                        && !ForesterUtil.isEmpty( getMainPanel().getTabbedPane().getTitleAt( i ) ) ) {
                    phy.setName( getMainPanel().getTabbedPane().getTitleAt( i ) );
                }
                trees.add( phy );
                getMainPanel().getTreePanels().get( i ).setEdited( false );
            }
            final PhylogenyWriter writer = new PhylogenyWriter();
            try {
                writer.toPhyloXML( file, trees, 0, ForesterUtil.LINE_SEPARATOR );
            }
            catch ( final IOException e ) {
                JOptionPane.showMessageDialog( this,
                                               "Failed to write to: " + file,
                                               "Error",
                                               JOptionPane.WARNING_MESSAGE );
            }
        }
    }

    private boolean writeAsNewHampshire( final Phylogeny t, boolean exception, final File file ) {
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toNewHampshire( t, false, true, getOptions().getNhConversionSupportValueStyle(), file );
        }
        catch ( final Exception e ) {
            exception = true;
            exceptionOccuredDuringSaveAs( e );
        }
        return exception;
    }

    private boolean writeAsNexus( final Phylogeny t, boolean exception, final File file ) {
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toNexus( file, t, getOptions().getNhConversionSupportValueStyle() );
        }
        catch ( final Exception e ) {
            exception = true;
            exceptionOccuredDuringSaveAs( e );
        }
        return exception;
    }

    private boolean writeAsPhyloXml( final Phylogeny t, boolean exception, final File file ) {
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( file, t, 0 );
        }
        catch ( final Exception e ) {
            exception = true;
            exceptionOccuredDuringSaveAs( e );
        }
        return exception;
    }

    private void writeToFile( final Phylogeny t ) {
        if ( t == null ) {
            return;
        }
        String initial_filename = null;
        if ( getMainPanel().getCurrentTreePanel().getTreeFile() != null ) {
            try {
                initial_filename = getMainPanel().getCurrentTreePanel().getTreeFile().getCanonicalPath();
            }
            catch ( final IOException e ) {
                initial_filename = null;
            }
        }
        if ( !ForesterUtil.isEmpty( initial_filename ) ) {
            _save_filechooser.setSelectedFile( new File( initial_filename ) );
        }
        else {
            _save_filechooser.setSelectedFile( new File( "" ) );
        }
        final File my_dir = getCurrentDir();
        if ( my_dir != null ) {
            _save_filechooser.setCurrentDirectory( my_dir );
        }
        final int result = _save_filechooser.showSaveDialog( _contentpane );
        final File file = _save_filechooser.getSelectedFile();
        setCurrentDir( _save_filechooser.getCurrentDirectory() );
        boolean exception = false;
        if ( ( file != null ) && ( result == JFileChooser.APPROVE_OPTION ) ) {
            if ( file.exists() ) {
                final int i = JOptionPane.showConfirmDialog( this,
                                                             file + " already exists.\nOverwrite?",
                                                             "Overwrite?",
                                                             JOptionPane.OK_CANCEL_OPTION,
                                                             JOptionPane.QUESTION_MESSAGE );
                if ( i != JOptionPane.OK_OPTION ) {
                    return;
                }
                else {
                    final File to = new File( file.getAbsoluteFile().toString() + Constants.BACKUP_FILE_SUFFIX );
                    try {
                        ForesterUtil.copyFile( file, to );
                    }
                    catch ( final Exception e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "Failed to create backup copy " + to,
                                                       "Failed to Create Backup Copy",
                                                       JOptionPane.WARNING_MESSAGE );
                    }
                    try {
                        file.delete();
                    }
                    catch ( final Exception e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "Failed to delete: " + file,
                                                       "Failed to Delete",
                                                       JOptionPane.WARNING_MESSAGE );
                    }
                }
            }
            if ( _save_filechooser.getFileFilter() == MainFrameApplication.nhfilter ) {
                exception = writeAsNewHampshire( t, exception, file );
            }
            else if ( _save_filechooser.getFileFilter() == MainFrameApplication.xmlfilter ) {
                exception = writeAsPhyloXml( t, exception, file );
            }
            else if ( _save_filechooser.getFileFilter() == MainFrameApplication.nexusfilter ) {
                exception = writeAsNexus( t, exception, file );
            }
            // "*.*":
            else {
                final String file_name = file.getName().trim().toLowerCase();
                if ( file_name.endsWith( ".nh" ) || file_name.endsWith( ".newick" ) || file_name.endsWith( ".phy" )
                        || file_name.endsWith( ".tree" ) ) {
                    exception = writeAsNewHampshire( t, exception, file );
                }
                else if ( file_name.endsWith( ".nex" ) || file_name.endsWith( ".nexus" ) ) {
                    exception = writeAsNexus( t, exception, file );
                }
                // XML is default:
                else {
                    exception = writeAsPhyloXml( t, exception, file );
                }
            }
            if ( !exception ) {
                getMainPanel().setTitleOfSelectedTab( file.getName() );
                getMainPanel().getCurrentTreePanel().setTreeFile( file );
                getMainPanel().getCurrentTreePanel().setEdited( false );
            }
        }
    }

    private void writeToGraphicsFile( final Phylogeny t, final GraphicsExportType type ) {
        if ( ( t == null ) || t.isEmpty() ) {
            return;
        }
        String initial_filename = "";
        if ( getMainPanel().getCurrentTreePanel().getTreeFile() != null ) {
            initial_filename = getMainPanel().getCurrentTreePanel().getTreeFile().toString();
        }
        if ( initial_filename.indexOf( '.' ) > 0 ) {
            initial_filename = initial_filename.substring( 0, initial_filename.lastIndexOf( '.' ) );
        }
        initial_filename = initial_filename + "." + type;
        _writetographics_filechooser.setSelectedFile( new File( initial_filename ) );
        final File my_dir = getCurrentDir();
        if ( my_dir != null ) {
            _writetographics_filechooser.setCurrentDirectory( my_dir );
        }
        final int result = _writetographics_filechooser.showSaveDialog( _contentpane );
        File file = _writetographics_filechooser.getSelectedFile();
        setCurrentDir( _writetographics_filechooser.getCurrentDirectory() );
        if ( ( file != null ) && ( result == JFileChooser.APPROVE_OPTION ) ) {
            if ( !file.toString().toLowerCase().endsWith( type.toString() ) ) {
                file = new File( file.toString() + "." + type );
            }
            if ( file.exists() ) {
                final int i = JOptionPane.showConfirmDialog( this,
                                                             file + " already exists. Overwrite?",
                                                             "Warning",
                                                             JOptionPane.OK_CANCEL_OPTION,
                                                             JOptionPane.WARNING_MESSAGE );
                if ( i != JOptionPane.OK_OPTION ) {
                    return;
                }
                else {
                    try {
                        file.delete();
                    }
                    catch ( final Exception e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "Failed to delete: " + file,
                                                       "Error",
                                                       JOptionPane.WARNING_MESSAGE );
                    }
                }
            }
            writePhylogenyToGraphicsFile( file.toString(), type );
        }
    }

    private void writeToPdf( final Phylogeny t ) {
        if ( ( t == null ) || t.isEmpty() ) {
            return;
        }
        String initial_filename = "";
        if ( getMainPanel().getCurrentTreePanel().getTreeFile() != null ) {
            initial_filename = getMainPanel().getCurrentTreePanel().getTreeFile().toString();
        }
        if ( initial_filename.indexOf( '.' ) > 0 ) {
            initial_filename = initial_filename.substring( 0, initial_filename.lastIndexOf( '.' ) );
        }
        initial_filename = initial_filename + ".pdf";
        _writetopdf_filechooser.setSelectedFile( new File( initial_filename ) );
        final File my_dir = getCurrentDir();
        if ( my_dir != null ) {
            _writetopdf_filechooser.setCurrentDirectory( my_dir );
        }
        final int result = _writetopdf_filechooser.showSaveDialog( _contentpane );
        File file = _writetopdf_filechooser.getSelectedFile();
        setCurrentDir( _writetopdf_filechooser.getCurrentDirectory() );
        if ( ( file != null ) && ( result == JFileChooser.APPROVE_OPTION ) ) {
            if ( !file.toString().toLowerCase().endsWith( ".pdf" ) ) {
                file = new File( file.toString() + ".pdf" );
            }
            if ( file.exists() ) {
                final int i = JOptionPane.showConfirmDialog( this,
                                                             file + " already exists. Overwrite?",
                                                             "WARNING",
                                                             JOptionPane.OK_CANCEL_OPTION,
                                                             JOptionPane.WARNING_MESSAGE );
                if ( i != JOptionPane.OK_OPTION ) {
                    return;
                }
            }
            printPhylogenyToPdf( file.toString() );
        }
    }

    public static MainFrameApplication createInstance( final Phylogeny[] phys, final Configuration config ) {
        return new MainFrameApplication( phys, config );
    }

    public static MainFrame createInstance( final Phylogeny[] phys,
                                            final Configuration config,
                                            final String title,
                                            final File current_dir ) {
        return new MainFrameApplication( phys, config, title, current_dir );
    }

    static MainFrame createInstance( final Phylogeny[] phys, final Configuration config, final String title ) {
        return new MainFrameApplication( phys, config, title );
    }

    static MainFrame createInstance( final Phylogeny[] phys, final String config_file_name, final String title ) {
        return new MainFrameApplication( phys, config_file_name, title );
    }

    static void setTextForGraphicsSizeChooserMenuItem( final JMenuItem mi, final Options o ) {
        mi.setText( "Enter Default Size for Graphics Export... (current: " + o.getPrintSizeX() + ", "
                + o.getPrintSizeY() + ")" );
    }

    static void setTextForPdfLineWidthChooserMenuItem( final JMenuItem mi, final Options o ) {
        mi.setText( "Enter Default Line Width for PDF Export... (current: " + o.getPrintLineWidth() + ")" );
    }

    static void warnIfNotPhyloXmlValidation( final Configuration c ) {
        if ( !c.isValidatePhyloXmlAgainstSchema() ) {
            JOptionPane
                    .showMessageDialog( null,
                                        ForesterUtil
                                                .wordWrap( "phyloXML XSD-based validation is turned off [enable with line 'validate_against_phyloxml_xsd_schem: true' in configuration file]",
                                                           80 ),
                                        "Warning",
                                        JOptionPane.WARNING_MESSAGE );
        }
    }
} // MainFrameApplication.

class DefaultFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith( ".nh" ) || file_name.endsWith( ".newick" ) || file_name.endsWith( ".phy" )
                || file_name.endsWith( ".nwk" ) || file_name.endsWith( ".phb" ) || file_name.endsWith( ".ph" )
                || file_name.endsWith( ".tr" ) || file_name.endsWith( ".dnd" ) || file_name.endsWith( ".tree" )
                || file_name.endsWith( ".nhx" ) || file_name.endsWith( ".xml" ) || file_name.endsWith( ".phyloxml" )
                || file_name.endsWith( "phylo.xml" ) || file_name.endsWith( ".pxml" ) || file_name.endsWith( ".nexus" )
                || file_name.endsWith( ".nx" ) || file_name.endsWith( ".nex" ) || file_name.endsWith( ".tre" )
                || file_name.endsWith( ".zip" ) || file_name.endsWith( ".tol" ) || file_name.endsWith( ".tolxml" )
                || file_name.endsWith( ".con" ) || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "All supported files (*.xml, *.phyloxml, *phylo.xml, *.nhx, *.nh, *.newick, *.nex, *.nexus, *.phy, *.tre, *.tree, *.tol, ...)";
    }
}

class GraphicsFileFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith( ".jpg" ) || file_name.endsWith( ".jpeg" ) || file_name.endsWith( ".png" )
                || file_name.endsWith( ".gif" ) || file_name.endsWith( ".bmp" ) || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "Image files (*.jpg, *.jpeg, *.png, *.gif, *.bmp)";
    }
}

class MsaFileFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith( ".msa" ) || file_name.endsWith( ".aln" ) || file_name.endsWith( ".fasta" )
                || file_name.endsWith( ".fas" ) || file_name.endsWith( ".fa" ) || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "Multiple sequence alignment files (*.msa, *.aln, *.fasta, *.fa, *.fas)";
    }
}

class NexusFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith( ".nex" ) || file_name.endsWith( ".nexus" ) || file_name.endsWith( ".nx" )
                || file_name.endsWith( ".tre" ) || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "Nexus files (*.nex, *.nexus, *.nx, *.tre)";
    }
} // NexusFilter

class NHFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith( ".nh" ) || file_name.endsWith( ".newick" ) || file_name.endsWith( ".phy" )
                || file_name.endsWith( ".tr" ) || file_name.endsWith( ".tree" ) || file_name.endsWith( ".dnd" )
                || file_name.endsWith( ".ph" ) || file_name.endsWith( ".phb" ) || file_name.endsWith( ".nwk" )
                || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "New Hampshire - Newick files (*.nh, *.newick, *.phy, *.tree, *.dnd, *.tr, *.ph, *.phb, *.nwk)";
    }
} // NHFilter

class NHXFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith( ".nhx" ) || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "NHX files (*.nhx) [deprecated]";
    }
}

class PdfFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        return f.getName().trim().toLowerCase().endsWith( ".pdf" ) || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "PDF files (*.pdf)";
    }
} // PdfFilter

class SequencesFileFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith( ".fasta" ) || file_name.endsWith( ".fa" ) || file_name.endsWith( ".fas" )
                || file_name.endsWith( ".seqs" ) || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "Sequences files (*.fasta, *.fa, *.fas, *.seqs )";
    }
}

class TolFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        final String file_name = f.getName().trim().toLowerCase();
        return ( file_name.endsWith( ".tol" ) || file_name.endsWith( ".tolxml" ) || file_name.endsWith( ".zip" ) || f
                .isDirectory() ) && ( !file_name.endsWith( ".xml.zip" ) );
    }

    @Override
    public String getDescription() {
        return "Tree of Life files (*.tol, *.tolxml)";
    }
} // TolFilter

class XMLFilter extends FileFilter {

    @Override
    public boolean accept( final File f ) {
        final String file_name = f.getName().trim().toLowerCase();
        return file_name.endsWith( ".xml" ) || file_name.endsWith( ".phyloxml" ) || file_name.endsWith( "phylo.xml" )
                || file_name.endsWith( ".pxml" ) || file_name.endsWith( ".zip" ) || f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "phyloXML files (*.xml, *.phyloxml, *phylo.xml, *.pxml, *.zip)";
    }
} // XMLFilter
