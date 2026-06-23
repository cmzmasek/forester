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
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.JTabbedPane;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.tools.SequenceAndTaxonomyDataObtainer;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.tol.TolParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.ws.seqdb.NcbiTaxonomyLineageService;

public final class MainFrameApplication extends MainFrame {

    private final static int FRAME_X_SIZE = 900;
    private final static int FRAME_Y_SIZE = 900;
    // Filters for the file-open dialog (classes defined in this file)
    private static final long serialVersionUID = -799735726778865234L;
    private static final boolean PREPROCESS_TREES = false;
    private final JFileChooser _open_filechooser;
    private final JFileChooser _open_filechooser_for_species_tree;
    // Application-only print menu items
    private JMenuItem _collapse_below_threshold;
    private JMenuItem _collapse_below_branch_length;
    // Others:
    double _min_not_collapse = AptxConstants.MIN_NOT_COLLAPSE_DEFAULT;
    double _min_not_collapse_bl = 0.001;

    private MainFrameApplication(final Phylogeny[] phys, final Configuration config) {
        _configuration = config;
        if (_configuration == null) {
            throw new IllegalArgumentException("configuration is null");
        }
        setVisible(false);
        setOptions(Options.createInstance(_configuration));
        _mainpanel = new MainPanel(_configuration, this);
        installTabContextMenu();
        _open_filechooser = null;
        _open_filechooser_for_species_tree = null;
        _save_filechooser = null;
        _writetopdf_filechooser = null;
        _writetographics_filechooser = null;
        _jmenubar = new JMenuBar();
        buildFileMenu();
        buildTypeMenu();
        _contentpane = getContentPane();
        _contentpane.setLayout(new BorderLayout());
        _contentpane.add(_mainpanel, BorderLayout.CENTER);
        // App is this big
        setSize(MainFrameApplication.FRAME_X_SIZE, MainFrameApplication.FRAME_Y_SIZE);
        // The window listener
        setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
        addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(final WindowEvent e) {
                exit();
            }
        });
        //   setVisible( true );
        if ((phys != null) && (phys.length > 0)) {
            AptxUtil.addPhylogeniesToTabs(phys, "", null, _configuration, _mainpanel);
            validate();
            getMainPanel().getControlPanel().showWholeAll();
            getMainPanel().getControlPanel().showWhole();
            // offer label extraction for a command-line / initially-loaded header tree too, but only after
            // the frame is realized (deferred), so the dialog never appears mid-construction
            final Phylogeny[] loaded = phys;
            javax.swing.SwingUtilities.invokeLater(() -> offerLabelExtraction(loaded));
        }
        // align the tree canvas with the resolved (light/dark) theme at startup
        updateTreeCanvasColors(getConfiguration().getUi());
        //activateSaveAllIfNeeded();
        // ...and its children
        _contentpane.repaint();
    }

    private MainFrameApplication(final Phylogeny[] phys, final Configuration config, final String title) {
        this(phys, config, title, null);
    }

    private MainFrameApplication(final Phylogeny[] phys,
                                 final Configuration config,
                                 final String title,
                                 final File current_dir) {
        super();
        _configuration = config;
        if (_configuration == null) {
            throw new IllegalArgumentException("configuration is null");
        }
        installLookAndFeel(_configuration.getUi());
        // the export/save choosers were created in the super-constructor, before the
        // look-and-feel above was installed; refresh them, so they are not left with the
        // platform default (e.g. the native macOS file dialog instead of FlatLaf).
        refreshFileChoosersLookAndFeel();
        if ((current_dir != null) && current_dir.canRead() && current_dir.isDirectory()) {
            setCurrentDir(current_dir);
        }
        // hide until everything is ready
        setVisible(false);
        setOptions(Options.createInstance(_configuration));
        // set title
        setTitle(AptxConstants.PRG_NAME + " " + AptxConstants.VERSION + " (" + AptxConstants.PRG_DATE + ")");
        _mainpanel = new MainPanel(_configuration, this);
        installTabContextMenu();
        // The file dialogs
        _open_filechooser = new JFileChooser();
        _open_filechooser.setMultiSelectionEnabled(true);
        _open_filechooser.addChoosableFileFilter(MainFrame.xmlfilter);
        _open_filechooser.addChoosableFileFilter(MainFrame.nhxfilter);
        _open_filechooser.addChoosableFileFilter(MainFrame.nhfilter);
        _open_filechooser.addChoosableFileFilter(MainFrame.nexusfilter);
        _open_filechooser.addChoosableFileFilter(MainFrame.tolfilter);
        _open_filechooser.addChoosableFileFilter(_open_filechooser.getAcceptAllFileFilter());
        _open_filechooser.setFileFilter(MainFrame.defaultfilter);
        _open_filechooser_for_species_tree = new JFileChooser();
        _open_filechooser_for_species_tree.setMultiSelectionEnabled(false);
        _open_filechooser_for_species_tree.addChoosableFileFilter(MainFrame.xmlfilter);
        _open_filechooser_for_species_tree.addChoosableFileFilter(MainFrame.tolfilter);
        _open_filechooser_for_species_tree.setFileFilter(MainFrame.xmlfilter);
        try {
            final String home_dir = System.getProperty("user.home");
            _open_filechooser.setCurrentDirectory(new File(home_dir));
            _open_filechooser_for_species_tree.setCurrentDirectory(new File(home_dir));
        } catch (final Exception e) {
            e.printStackTrace();
            // Do nothing. Not important.
        }
        // build the menu bar
        _jmenubar = new JMenuBar();
        if (_configuration.isApplyCustomGuiColors()) {
            _jmenubar.setBackground(getConfiguration().getGuiMenuBackgroundColor());
        }
        buildFileMenu();
        buildAnalysisMenu();
        buildToolsMenu();
        buildViewMenu();
        buildOptionItems();
        buildTypeMenu();
        buildSettingsMenu();
        buildHelpMenu();
        setJMenuBar(_jmenubar);
        _jmenubar.add(_help_jmenu);
        _contentpane = getContentPane();
        _contentpane.setLayout(new BorderLayout());
        _contentpane.add(_mainpanel, BorderLayout.CENTER);
        // App is this big
        setSize(MainFrameApplication.FRAME_X_SIZE, MainFrameApplication.FRAME_Y_SIZE);
        //        addWindowFocusListener( new WindowAdapter() {
        //
        //            @Override
        //            public void windowGainedFocus( WindowEvent e ) {
        //                requestFocusInWindow();
        //            }
        //        } );
        // The window listener
        setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
        addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(final WindowEvent e) {
                if (isUnsavedDataPresent()) {
                    final int r = JOptionPane.showConfirmDialog(null,
                            "Exit despite potentially unsaved changes?",
                            "Exit?",
                            JOptionPane.YES_NO_OPTION);
                    if (r != JOptionPane.YES_OPTION) {
                        return;
                    }
                } else {
                    final int r = JOptionPane
                            .showConfirmDialog(null, "Exit Archaeopteryx?", "Exit?", JOptionPane.YES_NO_OPTION);
                    if (r != JOptionPane.YES_OPTION) {
                        return;
                    }
                }
                exit();
            }
        });
        // The component listener
        addComponentListener(new ComponentAdapter() {

            @Override
            public void componentResized(final ComponentEvent e) {
                if (_mainpanel.getCurrentTreePanel() != null) {
                    _mainpanel.getCurrentTreePanel()
                            .calcParametersForPainting(_mainpanel.getCurrentTreePanel().getWidth(),
                                    _mainpanel.getCurrentTreePanel().getHeight());
                }
            }
        });
        requestFocusInWindow();
        // addKeyListener( this );
        setVisible(true);
        if ((phys != null) && (phys.length > 0)) {
            AptxUtil.addPhylogeniesToTabs(phys, title, null, _configuration, _mainpanel);
            validate();
            getMainPanel().getControlPanel().showWholeAll();
            getMainPanel().getControlPanel().showWhole();
        }
        // align the tree canvas with the resolved (light/dark) theme at startup
        updateTreeCanvasColors(getConfiguration().getUi());
        activateSaveAllIfNeeded();
        // ...and its children
        _contentpane.repaint();
        System.gc();
        // warm the persistent taxonomy cache off the EDT now, so the first colorize/fetch is snappy
        NcbiTaxonomyLineageService.getShared().primeAsync();
    }

    @Override
    public void actionPerformed(final ActionEvent e) {
        try {
            super.actionPerformed(e);
            final Object o = e.getSource();
            // Handle app-specific actions here:
            if (o == _open_item) {
                readPhylogeniesFromFile();
            }
            if (o == _new_item) {
                newTree();
            } else if (o == _close_item) {
                closeCurrentPane();
            } else if (o == _load_species_tree_item) {
                readSpeciesTreeFromFile();
            } else if (o == _obtain_seq_and_tax_information_jmi) {
                if (isSubtreeDisplayed()) {
                    return;
                }
                obtainSequenceAndTaxonomicInformation();
            } else if (o == _extract_label_data_jmi) {
                if (isSubtreeDisplayed()) {
                    return;
                }
                extractLabelData();
            } else if (o == _internal_number_are_confidence_for_nh_parsing_cbmi) {
                updateOptions(getOptions());
            } else if (o == _replace_underscores_cbmi) {
                if ((_extract_taxonomy_no_rbmi != null) && !_extract_taxonomy_no_rbmi.isSelected()) {
                    _extract_taxonomy_no_rbmi.setSelected(true);
                }
                updateOptions(getOptions());
            } else if (o == _allow_errors_in_distance_to_parent_cbmi) {
                updateOptions(getOptions());
            } else if (o == _collapse_below_threshold) {
                if (isSubtreeDisplayed()) {
                    return;
                }
                collapseBelowThreshold();
            } else if (o == _collapse_below_branch_length) {
                if (isSubtreeDisplayed()) {
                    return;
                }
                collapseBelowBranchLengthThreshold();
            } else if ((o == _extract_taxonomy_pfam_strict_rbmi) || (o == _extract_taxonomy_pfam_relaxed_rbmi)
                    || (o == _extract_taxonomy_agressive_rbmi)) {
                if (_replace_underscores_cbmi != null) {
                    _replace_underscores_cbmi.setSelected(false);
                }
                updateOptions(getOptions());
            } else if (o == _extract_taxonomy_no_rbmi) {
                updateOptions(getOptions());
            }
            _contentpane.repaint();
        } catch (final Exception ex) {
            AptxUtil.unexpectedException(ex);
        } catch (final Error err) {
            AptxUtil.unexpectedError(err);
        }
    }

    public void end() {
        _mainpanel.terminate();
        _contentpane.removeAll();
        setVisible(false);
        dispose();
    }

    @Override
    public MainPanel getMainPanel() {
        return _mainpanel;
    }

    private void closeCurrentPane() {
        if (getMainPanel().getCurrentTreePanel() != null) {
            if (getMainPanel().getCurrentTreePanel().isEdited()) {
                final int r = JOptionPane.showConfirmDialog(this,
                        "Close tab despite potentially unsaved changes?",
                        "Close Tab?",
                        JOptionPane.YES_NO_OPTION);
                if (r != JOptionPane.YES_OPTION) {
                    return;
                }
            }
            getMainPanel().closeCurrentPane();
            activateSaveAllIfNeeded();
        }
    }

    /** Browser-style: right-clicking a tree tab opens a context menu with a "Close Tab" item. */
    private void installTabContextMenu() {
        final JTabbedPane tabs = getMainPanel().getTabbedPane();
        if (tabs == null) {
            return;
        }
        tabs.addMouseListener(new MouseAdapter() {

            @Override
            public void mousePressed(final MouseEvent e) {
                maybeShowPopup(e);
            }

            @Override
            public void mouseReleased(final MouseEvent e) {
                maybeShowPopup(e);
            }

            private void maybeShowPopup(final MouseEvent e) {
                if (!e.isPopupTrigger()) {
                    return;
                }
                final int index = tabs.indexAtLocation(e.getX(), e.getY());
                if (index < 0) {
                    return; // not on a tab label
                }
                createTabPopupMenu(index).show(tabs, e.getX(), e.getY());
            }
        });
    }

    /** The context menu shown when right-clicking the tree tab at {@code index}. */
    JPopupMenu createTabPopupMenu(final int index) {
        final JPopupMenu popup = new JPopupMenu();
        final JMenuItem close = new JMenuItem("Close Tab");
        close.addActionListener(ae -> closeTabAt(index));
        popup.add(close);
        return popup;
    }

    /**
     * Selects the tab at {@code index} (making it the current pane) and closes it, with the same
     * unsaved-changes confirmation as File &gt; Close Tab.
     */
    void closeTabAt(final int index) {
        final JTabbedPane tabs = getMainPanel().getTabbedPane();
        if ((tabs == null) || (index < 0) || (index >= tabs.getTabCount())) {
            return;
        }
        tabs.setSelectedIndex(index);
        closeCurrentPane();
    }

    private void collapseBelowThreshold(final Phylogeny phy) {
        final double threshold = getMinNotCollapseConfidenceValue();
        final List<PhylogenyNode> candidates = AptxUtil.branchesToCollapseByConfidence(phy, threshold);
        if (candidates.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "No branches have a confidence value below " + threshold + " — nothing to collapse.",
                    "No branches collapsed",
                    JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        final int confirm = JOptionPane.showConfirmDialog(this,
                "Permanently collapse " + candidates.size() + " branch(es) with confidence below " + threshold
                        + " into polytomies?\nThis cannot be undone.",
                "Collapse " + candidates.size() + " branch(es)?",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.WARNING_MESSAGE);
        if (confirm != JOptionPane.OK_OPTION) {
            return;
        }
        for (final PhylogenyNode n : candidates) {
            PhylogenyMethods.removeNode(n, phy);
        }
        refreshAfterBranchCollapse(phy);
        JOptionPane.showMessageDialog(this,
                "Collapsed " + candidates.size() + " branch(es) with confidence below " + threshold + ".",
                "Collapsed " + candidates.size() + " branch(es)",
                JOptionPane.INFORMATION_MESSAGE);
    }

    /** Recomputes derived tree state and repaints after branches were removed by a collapse tool. */
    private void refreshAfterBranchCollapse(final Phylogeny phy) {
        phy.externalNodesHaveChanged();
        phy.clearHashIdToNodeMap();
        phy.recalculateNumberOfExternalDescendants(true);
        getCurrentTreePanel().resetNodeIdToDistToLeafMap();
        getCurrentTreePanel().updateSetOfCollapsedExternalNodes();
        getCurrentTreePanel().calculateLongestExtNodeInfo();
        getCurrentTreePanel().setNodeInPreorderToNull();
        getCurrentTreePanel().recalculateMaxDistanceToRoot();
        getCurrentTreePanel().resetPreferredSize();
        getCurrentTreePanel().setEdited(true);
        getCurrentTreePanel().repaint();
        repaint();
    }

    private void collapseBelowBranchLengthThreshold() {
        if (getCurrentTreePanel() != null) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ((phy != null) && !phy.isEmpty()) {
                final String s = (String) JOptionPane
                        .showInputDialog(this,
                                "Please enter the minimum branch length value\n",
                                "Minimal Branch Length Value",
                                JOptionPane.QUESTION_MESSAGE,
                                null,
                                null,
                                getMinNotCollapseBlValue());
                if (!ForesterUtil.isEmpty(s)) {
                    boolean success = true;
                    double m = 0.0;
                    final String m_str = s.trim();
                    if (!ForesterUtil.isEmpty(m_str)) {
                        try {
                            m = Double.parseDouble(m_str);
                        } catch (final Exception ex) {
                            success = false;
                        }
                    } else {
                        success = false;
                    }
                    if (success && (m >= 0.0)) {
                        setMinNotCollapseBlValue(m);
                        collapseBl(phy);
                    }
                }
            }
        }
    }

    private void collapseBelowThreshold() {
        if (getCurrentTreePanel() != null) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ((phy != null) && !phy.isEmpty()) {
                final String s = (String) JOptionPane.showInputDialog(this,
                        "Please enter the minimum confidence value\n",
                        "Minimal Confidence Value",
                        JOptionPane.QUESTION_MESSAGE,
                        null,
                        null,
                        getMinNotCollapseConfidenceValue());
                if (!ForesterUtil.isEmpty(s)) {
                    boolean success = true;
                    double m = 0.0;
                    final String m_str = s.trim();
                    if (!ForesterUtil.isEmpty(m_str)) {
                        try {
                            m = Double.parseDouble(m_str);
                        } catch (final Exception ex) {
                            success = false;
                        }
                    } else {
                        success = false;
                    }
                    if (success && (m >= 0.0)) {
                        setMinNotCollapseConfidenceValue(m);
                        collapseBelowThreshold(phy);
                    }
                }
            }
        }
    }

    private void collapseBl(final Phylogeny phy) {
        final double threshold = getMinNotCollapseBlValue();
        final List<PhylogenyNode> candidates = AptxUtil.branchesToCollapseByBranchLength(phy, threshold);
        if (candidates.isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "No branches have a branch length below " + threshold + " — nothing to collapse.",
                    "No branches collapsed",
                    JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        final int confirm = JOptionPane.showConfirmDialog(this,
                "Permanently collapse " + candidates.size() + " branch(es) shorter than " + threshold
                        + " into polytomies?\nThis cannot be undone.",
                "Collapse " + candidates.size() + " branch(es)?",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.WARNING_MESSAGE);
        if (confirm != JOptionPane.OK_OPTION) {
            return;
        }
        for (final PhylogenyNode n : candidates) {
            PhylogenyMethods.removeNode(n, phy);
        }
        refreshAfterBranchCollapse(phy);
        JOptionPane.showMessageDialog(this,
                "Collapsed " + candidates.size() + " branch(es) shorter than " + threshold + ".",
                "Collapsed " + candidates.size() + " branch(es)",
                JOptionPane.INFORMATION_MESSAGE);
    }

    private PhyloXmlParser createPhyloXmlParser() {
        PhyloXmlParser xml_parser = null;
        if (getConfiguration().isValidatePhyloXmlAgainstSchema()) {
            try {
                xml_parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
            } catch (final Exception e) {
                JOptionPane.showMessageDialog(this,
                        e.getLocalizedMessage(),
                        "failed to create validating XML parser",
                        JOptionPane.WARNING_MESSAGE);
            }
        }
        if (xml_parser == null) {
            xml_parser = PhyloXmlParser.createPhyloXmlParser();
        }
        return xml_parser;
    }

    private double getMinNotCollapseBlValue() {
        return _min_not_collapse_bl;
    }

    private double getMinNotCollapseConfidenceValue() {
        return _min_not_collapse;
    }

    private boolean isUnsavedDataPresent() {
        final List<TreePanel> tps = getMainPanel().getTreePanels();
        for (final TreePanel tp : tps) {
            if (tp.isEdited()) {
                return true;
            }
        }
        return false;
    }

    private void newTree() {
        final Phylogeny[] phys = new Phylogeny[1];
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode node = new PhylogenyNode();
        phy.setRoot(node);
        phy.setRooted(true);
        phys[0] = phy;
        AptxUtil.addPhylogeniesToTabs(phys, "", "", getConfiguration(), getMainPanel());
        _mainpanel.getControlPanel().showWhole();
        _mainpanel.getCurrentTreePanel().setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
        _mainpanel.getOptions().setPhylogenyGraphicsType(PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
        getMainPanel().getMainFrame().setSelectedTypeInTypeMenu(PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
        activateSaveAllIfNeeded();
        System.gc();
    }


    private void obtainSequenceAndTaxonomicInformation() {
        if (getCurrentTreePanel() != null) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ((phy != null) && !phy.isEmpty()) {
                final SequenceAndTaxonomyDataObtainer t = new SequenceAndTaxonomyDataObtainer(this,
                        _mainpanel.getCurrentTreePanel(),
                        phy.copy());
                new Thread(t).start();
            }
        }
    }

    private void preProcessTreesUponReading(final Phylogeny[] phys) {
        for (final Phylogeny phy : phys) {
            if ((phy != null) && !phy.isEmpty()) {
                for (final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if (n.isExternal()) {
                        if (n.getNodeData().isHasSequence()) {
                            final Sequence s = n.getNodeData().getSequence();
                            if (ForesterUtil.isEmpty(s.getGeneName()) || s.getGeneName().startsWith("LOC")) {
                                if ((s.getAccession() != null)
                                        && !ForesterUtil.isEmpty(s.getAccession().getValue())) {
                                    s.setGeneName(s.getAccession().getValue());
                                } else if (!ForesterUtil.isEmpty(n.getName())) {
                                    s.setGeneName(n.getName());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void readPhylogeniesFromFile() {
        boolean exception = false;
        Phylogeny[] phys = null;
        // Set an initial directory if none set yet
        final File my_dir = getCurrentDir();
        // Open file-open dialog and set current directory
        if (my_dir != null) {
            _open_filechooser.setCurrentDirectory(my_dir);
        }
        final int result = _open_filechooser.showOpenDialog(_contentpane);
        // All done: get the file
        final File[] files = _open_filechooser.getSelectedFiles();
        setCurrentDir(_open_filechooser.getCurrentDirectory());
        boolean nhx_or_nexus = false;
        if ((files != null) && (files.length > 0) && (result == JFileChooser.APPROVE_OPTION)) {
            for (final File file : files) {
                if ((file != null) && !file.isDirectory()) {
                    if (_mainpanel.getCurrentTreePanel() != null) {
                        _mainpanel.getCurrentTreePanel().setWaitCursor();
                    } else {
                        _mainpanel.setWaitCursor();
                    }
                    if ((_open_filechooser.getFileFilter() == MainFrame.nhfilter)
                            || (_open_filechooser.getFileFilter() == MainFrame.nhxfilter)) {
                        try {
                            final NHXParser nhx = new NHXParser();
                            setSpecialOptionsForNhxParser(nhx);
                            phys = PhylogenyMethods.readPhylogenies(nhx, file);
                            nhx_or_nexus = true;
                        } catch (final Exception e) {
                            exception = true;
                            exceptionOccuredDuringOpenFile(e);
                        }
                    } else if (_open_filechooser.getFileFilter() == MainFrame.xmlfilter) {
                        try {
                            final PhyloXmlParser xml_parser = createPhyloXmlParser();
                            phys = PhylogenyMethods.readPhylogenies(xml_parser, file);
                        } catch (final Exception e) {
                            exception = true;
                            exceptionOccuredDuringOpenFile(e);
                        }
                    } else if (_open_filechooser.getFileFilter() == MainFrame.tolfilter) {
                        try {
                            phys = PhylogenyMethods.readPhylogenies(new TolParser(), file);
                        } catch (final Exception e) {
                            exception = true;
                            exceptionOccuredDuringOpenFile(e);
                        }
                    } else if (_open_filechooser.getFileFilter() == MainFrame.nexusfilter) {
                        try {
                            final NexusPhylogeniesParser nex = new NexusPhylogeniesParser();
                            setSpecialOptionsForNexParser(nex);
                            phys = PhylogenyMethods.readPhylogenies(nex, file);
                            nhx_or_nexus = true;
                        } catch (final Exception e) {
                            exception = true;
                            exceptionOccuredDuringOpenFile(e);
                        }
                    }
                    // "*.*":
                    else {
                        try {
                            final PhylogenyParser parser = ParserUtils
                                    .createParserDependingOnFileType(file,
                                            getConfiguration()
                                                    .isValidatePhyloXmlAgainstSchema());
                            if (parser instanceof NexusPhylogeniesParser nex) {
                                setSpecialOptionsForNexParser(nex);
                                nhx_or_nexus = true;
                            } else if (parser instanceof NHXParser nhx) {
                                setSpecialOptionsForNhxParser(nhx);
                                nhx_or_nexus = true;
                            }
                            phys = PhylogenyMethods.readPhylogenies(parser, file);
                        } catch (final Exception e) {
                            exception = true;
                            exceptionOccuredDuringOpenFile(e);
                        }
                    }
                    if (_mainpanel.getCurrentTreePanel() != null) {
                        _mainpanel.getCurrentTreePanel().setArrowCursor();
                    } else {
                        _mainpanel.setArrowCursor();
                    }
                    if (!exception && (phys.length > 0)) {
                        boolean one_desc = false;
                        if (nhx_or_nexus) {
                            for (final Phylogeny phy : phys) {
                                if (getOptions().isInternalNumberAreConfidenceForNhParsing()) {
                                    PhylogenyMethods.transferInternalNodeNamesToConfidence(phy, "");
                                }
                                if (PhylogenyMethods.getMinimumDescendentsPerInternalNodes(phy) == 1) {
                                    one_desc = true;
                                    break;
                                }
                            }
                        }
                        if (PREPROCESS_TREES) {
                            preProcessTreesUponReading(phys);
                        }
                        AptxUtil.addPhylogeniesToTabs(phys,
                                file.getName(),
                                file.getAbsolutePath(),
                                getConfiguration(),
                                getMainPanel());
                        _mainpanel.getControlPanel().showWhole();
                        if (nhx_or_nexus && one_desc) {
                            JOptionPane
                                    .showMessageDialog(this,
                                            "One or more trees contain (a) node(s) with one descendant, "
                                                    + ForesterUtil.LINE_SEPARATOR
                                                    + "possibly indicating illegal parentheses within node names.",
                                            "Warning: Possible Error in New Hampshire Formatted Data",
                                            JOptionPane.WARNING_MESSAGE);
                        }
                        offerLabelExtraction(phys);
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
        _open_filechooser_for_species_tree.setSelectedFile(new File(""));
        if (my_dir != null) {
            _open_filechooser_for_species_tree.setCurrentDirectory(my_dir);
        }
        final int result = _open_filechooser_for_species_tree.showOpenDialog(_contentpane);
        final File file = _open_filechooser_for_species_tree.getSelectedFile();
        if ((file != null) && (result == JFileChooser.APPROVE_OPTION)) {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            try {
                t = factory.create(file, ParserUtils.createParserDependingOnFileType(file, true))[0];
            } catch (final Exception e) {
                exception = true;
                exceptionOccuredDuringOpenFile(e);
            }
            if (!exception && (t != null) && !t.isRooted()) {
                exception = true;
                t = null;
                JOptionPane.showMessageDialog(this,
                        "Species tree is not rooted",
                        "Species tree not loaded",
                        JOptionPane.ERROR_MESSAGE);
            }
            if (!exception && (t != null)) {
                final Set<Taxonomy> tax_set = new HashSet<>();
                for (final PhylogenyNodeIterator it = t.iteratorExternalForward(); it.hasNext(); ) {
                    final PhylogenyNode node = it.next();
                    if (!node.getNodeData().isHasTaxonomy()) {
                        exception = true;
                        t = null;
                        JOptionPane
                                .showMessageDialog(this,
                                        "Species tree contains external node(s) without taxonomy information",
                                        "Species tree not loaded",
                                        JOptionPane.ERROR_MESSAGE);
                        break;
                    } else {
                        if (tax_set.contains(node.getNodeData().getTaxonomy())) {
                            exception = true;
                            t = null;
                            JOptionPane
                                    .showMessageDialog(this,
                                            "Taxonomy [" + node.getNodeData().getTaxonomy().asSimpleText()
                                                    + "] is not unique in species tree",
                                            "Species tree not loaded",
                                            JOptionPane.ERROR_MESSAGE);
                            break;
                        } else {
                            tax_set.add(node.getNodeData().getTaxonomy());
                        }
                    }
                }
            }
            if (!exception && (t != null)) {
                setSpeciesTree(t);
                JOptionPane.showMessageDialog(this,
                        "Species tree successfully loaded",
                        "Species tree loaded",
                        JOptionPane.INFORMATION_MESSAGE);
            }
            _contentpane.repaint();
            System.gc();
        }
    }


    private void setMinNotCollapseBlValue(final double min_not_collapse_bl) {
        _min_not_collapse_bl = min_not_collapse_bl;
    }

    private void setMinNotCollapseConfidenceValue(final double min_not_collapse) {
        _min_not_collapse = min_not_collapse;
    }

    private void setSpecialOptionsForNexParser(final NexusPhylogeniesParser nex) {
        nex.setReplaceUnderscores(getOptions().isReplaceUnderscoresInNhParsing());
        nex.setTaxonomyExtraction(getOptions().getTaxonomyExtraction());
        nex.setParseBeastStyleExtendedTags(getOptions().isParseBeastStyleExtendedNexusTags());
    }

    private void setSpecialOptionsForNhxParser(final NHXParser nhx) {
        nhx.setReplaceUnderscores(getOptions().isReplaceUnderscoresInNhParsing());
        nhx.setTaxonomyExtraction(getOptions().getTaxonomyExtraction());
        nhx.setAllowErrorsInDistanceToParent(getOptions().isAllowErrorsInDistanceToParent());
        nhx.setParseBeastStyleExtendedTags(getOptions().isParseBeastStyleExtendedNexusTags());
    }

    void buildAnalysisMenu() {
        _analysis_menu = MainFrame.createMenu("Analysis", getConfiguration());
        _analysis_menu.setToolTipText("Reconcile gene and species trees, and infer lineages");
        _analysis_menu.add(_gsdi_item = new JMenuItem("GSDI (Generalized Speciation Duplication Inference)"));
        _analysis_menu.add(_gsdir_item = new JMenuItem("GSDIR (GSDI with re-rooting)"));
        _analysis_menu.add(_load_species_tree_item = new JMenuItem("Load Species Tree..."));
        customizeJMenuItem(_gsdi_item);
        customizeJMenuItem(_gsdir_item);
        customizeJMenuItem(_load_species_tree_item);
        _analysis_menu.addSeparator();
        _analysis_menu.add(_lineage_inference = new JMenuItem(INFER_ANCESTOR_TAXONOMIES));
        customizeJMenuItem(_lineage_inference);
        _lineage_inference.setToolTipText("Inference of ancestor taxonomies/lineages");
        _jmenubar.add(_analysis_menu);
    }

    @Override
    void refreshFileChoosersLookAndFeel() {
        super.refreshFileChoosersLookAndFeel();
        // these choosers are cached and standalone, so a runtime theme switch does not
        // reach them through the window tree; refresh them explicitly as well
        for (final JFileChooser fc : new JFileChooser[] { _open_filechooser, _open_filechooser_for_species_tree }) {
            if (fc != null) {
                SwingUtilities.updateComponentTreeUI(fc);
            }
        }
    }

    @Override
    void buildFileMenu() {
        _file_jmenu = MainFrame.createMenu("File", getConfiguration());
        _file_jmenu.setToolTipText("Read, save, and export trees; close tabs or exit");
        _file_jmenu.add(_open_item = new JMenuItem("Read Tree from File..."));
        if (getConfiguration().isEditable()) {
            _file_jmenu.addSeparator();
            _file_jmenu.add(_new_item = new JMenuItem("New"));
            _new_item.setToolTipText("to create a new tree with one node, as source for manual tree construction");
        }
        _file_jmenu.addSeparator();
        _file_jmenu.add(_save_item = new JMenuItem("Save Tree As..."));
        _file_jmenu.add(_save_all_item = new JMenuItem("Save All Trees As..."));
        _save_all_item.setToolTipText("Write all phylogenies to one file.");
        _save_all_item.setEnabled(false);
        _file_jmenu.addSeparator();
        _file_jmenu.add(_write_to_pdf_item = new JMenuItem("Export to PDF file ..."));
        _file_jmenu.add(_write_to_svg_item = new JMenuItem("Export to SVG file..."));
        _write_to_svg_item.setToolTipText("Scalable vector graphics for publication (edit in Illustrator/Inkscape)");
        _file_jmenu.add(_write_to_eps_item = new JMenuItem("Export to EPS file..."));
        _write_to_eps_item.setToolTipText("Encapsulated PostScript vector graphics for publication");
        if (AptxUtil.canWriteFormat("tif") || AptxUtil.canWriteFormat("tiff")
                || AptxUtil.canWriteFormat("TIF")) {
            _file_jmenu.add(_write_to_tif_item = new JMenuItem("Export to TIFF file..."));
        }
        _file_jmenu.add(_write_to_png_item = new JMenuItem("Export to PNG file..."));
        _file_jmenu.add(_write_to_jpg_item = new JMenuItem("Export to JPG file..."));
        _file_jmenu.addSeparator();
        _file_jmenu.add(_close_item = new JMenuItem("Close Tab"));
        _close_item.setToolTipText("To close the current pane.");
        _close_item.setEnabled(true);
        _file_jmenu.addSeparator();
        _file_jmenu.add(_exit_item = new JMenuItem("Exit"));
        customizeJMenuItem(_open_item);
          customizeJMenuItem(_save_item);
        if (getConfiguration().isEditable()) {
            customizeJMenuItem(_new_item);
        }
        customizeJMenuItem(_close_item);
        customizeJMenuItem(_save_all_item);
        customizeJMenuItem(_write_to_pdf_item);
        customizeJMenuItem(_write_to_svg_item);
        customizeJMenuItem(_write_to_eps_item);
        customizeJMenuItem(_write_to_png_item);
        customizeJMenuItem(_write_to_jpg_item);
        customizeJMenuItem(_write_to_tif_item);
        customizeJMenuItem(_exit_item);
        _jmenubar.add(_file_jmenu);
    }

    /**
     * Builds the checkbox/radio "menu items" that exist only as the backing model for the Settings
     * dialog (which drives them via doClick()). They are no longer shown in any menu -- the old,
     * never-displayed "Options" menu has been retired -- so they are created detached here. The
     * customize* calls below set each item's initial state and, crucially, wire its action listener.
     */
    void buildOptionItems() {
        _ext_node_dependent_cladogram_rbmi = new JRadioButtonMenuItem(MainFrame.NONUNIFORM_CLADOGRAMS_LABEL);
        _non_lined_up_cladograms_rbmi = new JRadioButtonMenuItem(NON_LINED_UP_CLADOGRAMS_LABEL);
        ButtonGroup _radio_group_1 = new ButtonGroup();
        _radio_group_1.add(_ext_node_dependent_cladogram_rbmi);
        _radio_group_1.add(_non_lined_up_cladograms_rbmi);
        _show_overview_cbmi = new JCheckBoxMenuItem(SHOW_OVERVIEW_LABEL);
        _show_scale_cbmi = new JCheckBoxMenuItem(DISPLAY_SCALE_LABEL);
        _show_default_node_shapes_internal_cbmi = new JCheckBoxMenuItem(DISPLAY_NODE_BOXES_LABEL_INT);
        _internal_labels_above_branch_rbmi = new JRadioButtonMenuItem(MainFrame.INTERNAL_LABELS_ABOVE_BRANCH_LABEL);
        _internal_labels_above_branch_rbmi.setToolTipText(MainFrame.INTERNAL_LABELS_ABOVE_BRANCH_TIP);
        _internal_labels_right_of_node_rbmi = new JRadioButtonMenuItem(MainFrame.INTERNAL_LABELS_RIGHT_OF_NODE_LABEL);
        _internal_labels_right_of_node_rbmi.setToolTipText(MainFrame.INTERNAL_LABELS_RIGHT_OF_NODE_TIP);
        final ButtonGroup _radio_group_internal_labels = new ButtonGroup();
        _radio_group_internal_labels.add(_internal_labels_above_branch_rbmi);
        _radio_group_internal_labels.add(_internal_labels_right_of_node_rbmi);
        _show_default_node_shapes_external_cbmi = new JCheckBoxMenuItem(DISPLAY_NODE_BOXES_LABEL_EXT);
        _show_default_node_shapes_for_marked_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_NODE_BOXES_LABEL_MARKED);
        _collapsed_with_average_height_cbmi = new JCheckBoxMenuItem("Proportional Height of Collapsed Subtrees");
        _show_abbreviated_labels_for_collapsed_nodes_cbmi = new JCheckBoxMenuItem("Add Abbreviated Labels to Collapsed Subtrees");
        _line_up_renderable_data_cbmi = new JCheckBoxMenuItem(MainFrame.LINE_UP_RENDERABLE_DATA);
        if (getConfiguration().doDisplayOption(Configuration.show_domain_architectures)) {
            _right_line_up_domains_cbmi = new JCheckBoxMenuItem(MainFrame.RIGHT_LINE_UP_DOMAINS);
            _show_domain_labels = new JCheckBoxMenuItem(MainFrame.SHOW_DOMAIN_LABELS_LABEL);
        }
        _show_confidence_stddev_cbmi = new JCheckBoxMenuItem(SHOW_CONF_STDDEV_LABEL);
        _show_mad_confidence_cbmi = new JCheckBoxMenuItem(MainFrame.SHOW_MAD_CONF_LABEL);
        _color_labels_same_as_parent_branch = new JCheckBoxMenuItem(COLOR_LABELS_LABEL);
        _color_labels_same_as_parent_branch.setToolTipText(MainFrame.COLOR_LABELS_TIP);
        _abbreviate_scientific_names = new JCheckBoxMenuItem(ABBREV_SN_LABEL);
        _use_italic_scientific_names_cbmi = new JCheckBoxMenuItem(MainFrame.ITALIC_SN_LABEL);
        _use_italic_scientific_names_cbmi.setToolTipText(MainFrame.ITALIC_SN_TIP);
        _outline_fonts_in_vector_export_cbmi = new JCheckBoxMenuItem(MainFrame.OUTLINE_FONTS_VECTOR_LABEL);
        _outline_fonts_in_vector_export_cbmi.setToolTipText(MainFrame.OUTLINE_FONTS_VECTOR_TIP);
        _label_direction_cbmi = new JCheckBoxMenuItem(LABEL_DIRECTION_LABEL);
        _label_direction_cbmi.setToolTipText(LABEL_DIRECTION_TIP);
        _color_all_found_nodes_when_coloring_subtree_cbmi = new JCheckBoxMenuItem("Colorize All Found Nodes When Colorizing Subtree(s)");
        _antialias_print_cbmi = new JCheckBoxMenuItem("Antialias (export)");
        _print_black_and_white_cbmi = new JCheckBoxMenuItem("Export in Black and White");
        _graphics_export_visible_only_cbmi = new JCheckBoxMenuItem("Limit to Visible ('Screenshot') for PNG and JPG export");
        _internal_number_are_confidence_for_nh_parsing_cbmi = new JCheckBoxMenuItem("Internal Node Names are Confidence Values");
        _replace_underscores_cbmi = new JCheckBoxMenuItem("Replace Underscores with Spaces");
        _parse_beast_style_extended_nexus_tags_cbmi = new JCheckBoxMenuItem("Parse BEAST-style extended Newick/Nexus tags");
        _parse_beast_style_extended_nexus_tags_cbmi
                .setToolTipText("to parse elements in the form of \"[&!color=#800080]\" in Newick/Nexus formatted trees");
        _allow_errors_in_distance_to_parent_cbmi = new JCheckBoxMenuItem("Ignore Distance Values Format Errors");
        _extract_taxonomy_no_rbmi = new JRadioButtonMenuItem("No Taxonomy Extraction");
        _extract_taxonomy_pfam_strict_rbmi = new JRadioButtonMenuItem("Extract Taxonomy Codes/Ids from Pfam-style Node Names");
        _extract_taxonomy_pfam_relaxed_rbmi = new JRadioButtonMenuItem("Extract Taxonomy Codes/Ids from Pfam-style like Node Names");
        _extract_taxonomy_agressive_rbmi = new JRadioButtonMenuItem("Extract Taxonomy Codes/Ids/Scientific Names from Node Names");
        _extract_taxonomy_pfam_strict_rbmi
                .setToolTipText("To extract taxonomy codes/ids from node names in the form of e.g. \"BCL2_MOUSE/123-304\" or \"BCL2_10090/123-304\"");
        _extract_taxonomy_pfam_relaxed_rbmi
                .setToolTipText("To extract taxonomy codes/ids from node names in the form of e.g. \"bax_MOUSE\" or \"bax_10090\"");
        _extract_taxonomy_agressive_rbmi
                .setToolTipText("To extract taxonomy codes/ids or scientific names from node names in the form of e.g. \"MOUSE\" or \"10090\" or \"xyz_Nematostella_vectensis\"");
        ButtonGroup _radio_group_2 = new ButtonGroup();
        _radio_group_2.add(_extract_taxonomy_no_rbmi);
        _radio_group_2.add(_extract_taxonomy_pfam_strict_rbmi);
        _radio_group_2.add(_extract_taxonomy_pfam_relaxed_rbmi);
        _radio_group_2.add(_extract_taxonomy_agressive_rbmi);
        _use_brackets_for_conf_in_nh_export_cbmi = new JCheckBoxMenuItem(USE_BRACKETS_FOR_CONF_IN_NH_LABEL);
        _use_brackets_for_conf_in_nh_export_cbmi
                .setToolTipText("e.g. \"0.1[90]\" for a branch with support 90 and a length of 0.1");
        _use_internal_names_for_conf_in_nh_export_cbmi = new JCheckBoxMenuItem(USE_INTERNAL_NAMES_FOR_CONF_IN_NH_LABEL);
        customizeCheckBoxMenuItem(_show_default_node_shapes_external_cbmi,
                getOptions().isShowDefaultNodeShapesExternal());
        customizeCheckBoxMenuItem(_show_default_node_shapes_internal_cbmi,
                getOptions().isShowDefaultNodeShapesInternal());
        customizeRadioButtonMenuItem(_internal_labels_above_branch_rbmi,
                getOptions().isInternalLabelsAboveBranch());
        customizeRadioButtonMenuItem(_internal_labels_right_of_node_rbmi,
                !getOptions().isInternalLabelsAboveBranch());
        customizeCheckBoxMenuItem(_show_default_node_shapes_for_marked_cbmi,
                getOptions().isShowDefaultNodeShapesForMarkedNodes());
        customizeCheckBoxMenuItem(_color_labels_same_as_parent_branch,
                getOptions().isColorLabelsSameAsParentBranch());
        customizeCheckBoxMenuItem(_show_domain_labels, getOptions().isShowDomainLabels());
        customizeCheckBoxMenuItem(_abbreviate_scientific_names, getOptions().isAbbreviateScientificTaxonNames());
        customizeCheckBoxMenuItem(_use_italic_scientific_names_cbmi, getOptions().isUseItalicScientificNames());
        customizeCheckBoxMenuItem(_outline_fonts_in_vector_export_cbmi,
                getOptions().isOutlineFontsInVectorExport());
        customizeCheckBoxMenuItem(_show_scale_cbmi, getOptions().isShowScale());
        customizeCheckBoxMenuItem(_collapsed_with_average_height_cbmi, getOptions().isCollapsedWithAverageHeigh());
        customizeCheckBoxMenuItem(_show_abbreviated_labels_for_collapsed_nodes_cbmi,
                getOptions().isShowAbbreviatedLabelsForCollapsedNodes());
        customizeRadioButtonMenuItem(_non_lined_up_cladograms_rbmi,
                getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP);
        customizeRadioButtonMenuItem(_ext_node_dependent_cladogram_rbmi,
                getOptions().getCladogramType() == CLADOGRAM_TYPE.LINED_UP);
        customizeCheckBoxMenuItem(_show_overview_cbmi, getOptions().isShowOverview());
        customizeCheckBoxMenuItem(_label_direction_cbmi,
                getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL);
        customizeCheckBoxMenuItem(_antialias_print_cbmi, getOptions().isAntialiasPrint());
        customizeCheckBoxMenuItem(_print_black_and_white_cbmi, getOptions().isPrintBlackAndWhite());
        customizeCheckBoxMenuItem(_internal_number_are_confidence_for_nh_parsing_cbmi,
                getOptions().isInternalNumberAreConfidenceForNhParsing());
        customizeRadioButtonMenuItem(_extract_taxonomy_no_rbmi,
                getOptions().getTaxonomyExtraction() == TAXONOMY_EXTRACTION.NO);
        customizeRadioButtonMenuItem(_extract_taxonomy_pfam_strict_rbmi,
                getOptions().getTaxonomyExtraction() == TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT);
        customizeRadioButtonMenuItem(_extract_taxonomy_pfam_relaxed_rbmi,
                getOptions().getTaxonomyExtraction() == TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED);
        customizeRadioButtonMenuItem(_extract_taxonomy_agressive_rbmi,
                getOptions().getTaxonomyExtraction() == TAXONOMY_EXTRACTION.AGGRESSIVE);
        customizeCheckBoxMenuItem(_replace_underscores_cbmi, getOptions().isReplaceUnderscoresInNhParsing());
        customizeCheckBoxMenuItem(_allow_errors_in_distance_to_parent_cbmi,
                getOptions().isAllowErrorsInDistanceToParent());
        customizeCheckBoxMenuItem(_color_all_found_nodes_when_coloring_subtree_cbmi,
                getOptions().isColorAllFoundNodesWhenColoringSubtree());
        customizeCheckBoxMenuItem(_parse_beast_style_extended_nexus_tags_cbmi,
                getOptions().isParseBeastStyleExtendedNexusTags());
        customizeCheckBoxMenuItem(_graphics_export_visible_only_cbmi, getOptions().isGraphicsExportVisibleOnly());
        customizeCheckBoxMenuItem(_show_confidence_stddev_cbmi, getOptions().isShowConfidenceStddev());
        customizeCheckBoxMenuItem(_show_mad_confidence_cbmi, getOptions().isShowMadConfidence());
        customizeCheckBoxMenuItem(_use_brackets_for_conf_in_nh_export_cbmi,
                getOptions()
                        .getNhConversionSupportValueStyle() == NH_CONVERSION_SUPPORT_VALUE_STYLE.IN_SQUARE_BRACKETS);
        customizeCheckBoxMenuItem(_use_internal_names_for_conf_in_nh_export_cbmi,
                getOptions()
                        .getNhConversionSupportValueStyle() == NH_CONVERSION_SUPPORT_VALUE_STYLE.AS_INTERNAL_NODE_NAMES);
        customizeCheckBoxMenuItem(_line_up_renderable_data_cbmi, getOptions().isLineUpRendarableNodeData());
        customizeCheckBoxMenuItem(_right_line_up_domains_cbmi, getOptions().isRightLineUpDomains());
    }

    void buildSettingsMenu() {
        _settings_jmenu = createMenu("Settings", getConfiguration());
        _settings_jmenu.setToolTipText("Display, node, export, file and theme settings");
        // The "Settings" entry acts as a one-click launcher: selecting it opens the Settings dialog
        // instead of dropping down an (empty) menu.
        _settings_jmenu.addMenuListener(new javax.swing.event.MenuListener() {

            @Override
            public void menuSelected(final javax.swing.event.MenuEvent e) {
                javax.swing.MenuSelectionManager.defaultManager().clearSelectedPath();
                SwingUtilities.invokeLater(() -> new SettingsDialog(MainFrameApplication.this).setVisible(true));
            }

            @Override
            public void menuDeselected(final javax.swing.event.MenuEvent e) {
            }

            @Override
            public void menuCanceled(final javax.swing.event.MenuEvent e) {
            }
        });
        _jmenubar.add(_settings_jmenu);
    }

    void buildToolsMenu() {
        _tools_menu = createMenu("Tools", getConfiguration());
        _tools_menu.setToolTipText("Root, prune, colorize, collapse, and fetch data for the current tree");
        // Rooting
        _tools_menu.add(_mad_root_item = new JMenuItem("MAD-Root"));
        _mad_root_item.setToolTipText("Root by Minimal Ancestor Deviation (Tria et al. 2017); requires branch lengths");
        customizeJMenuItem(_mad_root_item);
        _tools_menu.add(_midpoint_root_item = new JMenuItem("Midpoint-Root"));
        customizeJMenuItem(_midpoint_root_item);
        _tools_menu.addSeparator();
        // Pruning
        _tools_menu.add(_delete_selected_nodes_item = new JMenuItem("Delete Selected Nodes"));
        _delete_selected_nodes_item.setToolTipText("To delete all selected external nodes");
        customizeJMenuItem(_delete_selected_nodes_item);
        _tools_menu.add(_delete_not_selected_nodes_item = new JMenuItem("Retain Selected Nodes"));
        _delete_not_selected_nodes_item.setToolTipText("To delete all not selected external nodes");
        customizeJMenuItem(_delete_not_selected_nodes_item);
        _tools_menu.addSeparator();
        // Colorizing
        _tools_menu.add(_color_rank_jmi = new JMenuItem("Colorize Subtrees via Taxonomic Rank"));
        customizeJMenuItem(_color_rank_jmi);
        _color_rank_jmi.setToolTipText("for example, at \"Class\" level, colorize mammal specific subtree red");
        _tools_menu.add(_clade_bands_jmi = new JMenuItem("Annotate Clades by Rank…"));
        customizeJMenuItem(_clade_bands_jmi);
        _clade_bands_jmi.setToolTipText("mark clades at a chosen rank with shaded boxes or right-edge bars + labels");
        _tools_menu.addSeparator();
        // Clearing styles & colors
        _tools_menu.add(_remove_visual_styles_item = new JMenuItem("Delete All Visual Styles From Nodes"));
        _remove_visual_styles_item
                .setToolTipText("To remove all node visual styles (fonts, colors) from the current phylogeny");
        customizeJMenuItem(_remove_visual_styles_item);
        _tools_menu.add(_remove_branch_color_item = new JMenuItem("Delete All Colors From Branches"));
        _remove_branch_color_item.setToolTipText("To remove all branch color values from the current phylogeny");
        customizeJMenuItem(_remove_branch_color_item);
        _tools_menu.addSeparator();
        // Collapsing
        final JMenu collapse_menu = createMenu("Collapse Branches", getConfiguration());
        collapse_menu.setFont(MainFrame.menu_font); // match the sibling items (createMenu sets the font only in custom-colors mode)
        collapse_menu.setToolTipText("Permanently collapse weakly-supported or very short branches into polytomies");
        collapse_menu.add(_collapse_below_threshold = new JMenuItem("Collapse Weakly-Supported Branches…"));
        customizeJMenuItem(_collapse_below_threshold);
        _collapse_below_threshold.setToolTipText(
                "Permanently collapse internal branches whose confidence is below a threshold into polytomies (multifurcations). Cannot be undone.");
        collapse_menu.add(_collapse_below_branch_length = new JMenuItem("Collapse Short Branches…"));
        customizeJMenuItem(_collapse_below_branch_length);
        _collapse_below_branch_length.setToolTipText(
                "Permanently collapse internal branches shorter than a threshold into polytomies (multifurcations). Cannot be undone.");
        _tools_menu.add(collapse_menu);
        _tools_menu.addSeparator();
        // Data retrieval
        _tools_menu
                .add(_obtain_seq_and_tax_information_jmi = new JMenuItem(OBTAIN_SEQUENCE_AND_TAXONOMIC_INFORMATION));
        customizeJMenuItem(_obtain_seq_and_tax_information_jmi);
        _obtain_seq_and_tax_information_jmi
                .setToolTipText("To add additional sequence information and detailed taxonomic information (from UniProt/EMBL-GenBank and UniProt Taxonomy)");
        _tools_menu.add(_extract_label_data_jmi = new JMenuItem("Extract Data from Labels…"));
        customizeJMenuItem(_extract_label_data_jmi);
        _extract_label_data_jmi.setToolTipText(
                "Parse UniProt or GenBank/RefSeq FASTA-header node names into accession, description, gene and taxonomy fields (offline, no network); only empty fields are filled");
        _jmenubar.add(_tools_menu);
    }

    @Override
    void close() {
        if (isUnsavedDataPresent()) {
            final int r = JOptionPane.showConfirmDialog(this,
                    "Exit despite potentially unsaved changes?",
                    "Exit?",
                    JOptionPane.YES_NO_OPTION);
            if (r != JOptionPane.YES_OPTION) {
                return;
            }
        }
        exit();
    }

    void exit() {
        removeAllTextFrames();
        _mainpanel.terminate();
        _contentpane.removeAll();
        setVisible(false);
        dispose();
        // System.exit( 0 ); //TODO reconfirm that this is OK, then remove.
    }



    public static MainFrameApplication createInstance(final Phylogeny[] phys, final Configuration config) {
        return new MainFrameApplication(phys, config);
    }

    public static MainFrame createInstance(final Phylogeny[] phys,
                                           final Configuration config,
                                           final String title,
                                           final File current_dir) {
        return new MainFrameApplication(phys, config, title, current_dir);
    }

    static MainFrame createInstance(final Phylogeny[] phys, final Configuration config, final String title) {
        return new MainFrameApplication(phys, config, title);
    }
} // MainFrameApplication.
