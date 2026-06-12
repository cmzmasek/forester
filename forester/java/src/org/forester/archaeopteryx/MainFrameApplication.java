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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.JTabbedPane;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.forester.analysis.TaxonomyDataManager;
import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.tools.SequenceDataRetriver;
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
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

public final class MainFrameApplication extends MainFrame {

    private final static int FRAME_X_SIZE = 900;
    private final static int FRAME_Y_SIZE = 900;
    // Filters for the file-open dialog (classes defined in this file)
    private static final long serialVersionUID = -799735726778865234L;
    private static final boolean PREPROCESS_TREES = false;
    private final JFileChooser _values_filechooser;
    private final JFileChooser _open_filechooser;
    private final JFileChooser _open_filechooser_for_species_tree;
    private final JFileChooser _annotations_filechooser;
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
        _values_filechooser = null;
        _annotations_filechooser = null;
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
        // Expression
        _values_filechooser = new JFileChooser();
        _values_filechooser.setMultiSelectionEnabled(false);
        // Annotations
        _annotations_filechooser = new JFileChooser();
        _annotations_filechooser.setMultiSelectionEnabled(false);
        try {
            final String home_dir = System.getProperty("user.home");
            _open_filechooser.setCurrentDirectory(new File(home_dir));
            _open_filechooser_for_species_tree.setCurrentDirectory(new File(home_dir));
            _values_filechooser.setCurrentDirectory(new File(home_dir));
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
        buildFontSizeMenu();
        buildOptionsMenu();
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
    }

    private MainFrameApplication(final Phylogeny[] phys, final String config_file, final String title) {
        // Reads the config file (false, false => not url, not applet):
        this(phys, new Configuration(config_file, false, false, true), title);
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
            } else if (o == _replace_names_item) {
                replaceNodeNames();
            } else if (o == _close_item) {
                closeCurrentPane();
            } else if (o == _load_species_tree_item) {
                readSpeciesTreeFromFile();
            } else if (o == _obtain_detailed_taxonomic_information_jmi) {
                if (isSubtreeDisplayed()) {
                    return;
                }
                obtainDetailedTaxonomicInformation();
            } else if (o == _obtain_detailed_taxonomic_information_deleting_jmi) {
                if (isSubtreeDisplayed()) {
                    return;
                }
                obtainDetailedTaxonomicInformationDelete();
            } else if (o == _obtain_seq_information_jmi) {
                obtainSequenceInformation();
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
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        final List<PhylogenyNode> to_be_removed = new ArrayList<>();
        double min_support = Double.MAX_VALUE;
        boolean conf_present = false;
        while (it.hasNext()) {
            final PhylogenyNode n = it.next();
            if (!n.isExternal() && !n.isRoot()) {
                final List<Confidence> c = n.getBranchData().getConfidences();
                if ((c != null) && (c.size() > 0)) {
                    conf_present = true;
                    double max = 0;
                    for (final Confidence confidence : c) {
                        if (confidence.getValue() > max) {
                            max = confidence.getValue();
                        }
                    }
                    if (max < getMinNotCollapseConfidenceValue()) {
                        to_be_removed.add(n);
                    }
                    if (max < min_support) {
                        min_support = max;
                    }
                }
            }
        }
        if (conf_present) {
            for (final PhylogenyNode node : to_be_removed) {
                PhylogenyMethods.removeNode(node, phy);
            }
            if (to_be_removed.size() > 0) {
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
            if (to_be_removed.size() > 0) {
                JOptionPane.showMessageDialog(this,
                        "Collapsed " + to_be_removed.size()
                                + " branches with\nconfidence values below "
                                + getMinNotCollapseConfidenceValue(),
                        "Collapsed " + to_be_removed.size() + " branches",
                        JOptionPane.INFORMATION_MESSAGE);
            } else {
                JOptionPane.showMessageDialog(this,
                        "No branch collapsed,\nminimum confidence value per branch is "
                                + min_support,
                        "No branch collapsed",
                        JOptionPane.INFORMATION_MESSAGE);
            }
        } else {
            JOptionPane.showMessageDialog(this,
                    "No branch collapsed because no confidence values present",
                    "No confidence values present",
                    JOptionPane.INFORMATION_MESSAGE);
        }
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
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        final List<PhylogenyNode> to_be_removed = new ArrayList<>();
        double min_bl = Double.MAX_VALUE;
        boolean bl_present = false;
        while (it.hasNext()) {
            final PhylogenyNode n = it.next();
            if (!n.isExternal() && !n.isRoot()) {
                final double bl = n.getDistanceToParent();
                if (bl != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT) {
                    bl_present = true;
                    if (bl < getMinNotCollapseBlValue()) {
                        to_be_removed.add(n);
                    }
                    if (bl < min_bl) {
                        min_bl = bl;
                    }
                }
            }
        }
        if (bl_present) {
            for (final PhylogenyNode node : to_be_removed) {
                PhylogenyMethods.removeNode(node, phy);
            }
            if (to_be_removed.size() > 0) {
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
            if (to_be_removed.size() > 0) {
                JOptionPane.showMessageDialog(this,
                        "Collapsed " + to_be_removed.size()
                                + " branches with\nbranch length values below "
                                + getMinNotCollapseBlValue(),
                        "Collapsed " + to_be_removed.size() + " branches",
                        JOptionPane.INFORMATION_MESSAGE);
            } else {
                JOptionPane.showMessageDialog(this,
                        "No branch collapsed,\nminimum branch length is " + min_bl,
                        "No branch collapsed",
                        JOptionPane.INFORMATION_MESSAGE);
            }
        } else {
            JOptionPane.showMessageDialog(this,
                    "No branch collapsed because no branch length values present",
                    "No branch length values present",
                    JOptionPane.INFORMATION_MESSAGE);
        }
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


    private void replaceNodeNames() {

        if ((getCurrentTreePanel() == null) || (getCurrentTreePanel().getPhylogeny() == null)) {
            JOptionPane.showMessageDialog(this,
                    "Need to load phylogenetic tree first",
                    "Can Not Replace Node Names",
                    JOptionPane.WARNING_MESSAGE);
            return;
        }
        if (isSubtreeDisplayed()) {
            JOptionPane.showMessageDialog(this,
                    "Cannot replace node names in sub-tree",
                    "Can Not Replace Node Names In Sub-Tree",
                    JOptionPane.WARNING_MESSAGE);
            return;
        }
        final File my_dir = getCurrentDir();
        if (my_dir != null) {
            _values_filechooser.setCurrentDirectory(my_dir);
        }
        final int result = _values_filechooser.showOpenDialog(_contentpane);
        final File file = _values_filechooser.getSelectedFile();
        if ((file != null) && (file.length() > 0) && (result == JFileChooser.APPROVE_OPTION)) {
            BasicTable<String> t;
            try {
                t = BasicTableParser.parse(file, '\t');
                if (t.getNumberOfColumns() < 2) {
                    t = BasicTableParser.parse(file, ',');
                }
            } catch (final IOException e) {
                JOptionPane.showMessageDialog(this,
                        e.getMessage(),
                        "Could Not Read Mapping File",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }
            if (t.getNumberOfColumns() != 2) {
                JOptionPane.showMessageDialog(this,
                        "Table contains " + t.getNumberOfColumns() + " column(s)",
                        "Problem with Mapping File",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if (t.getNumberOfRows() != phy.getNumberOfExternalNodes()) {
                JOptionPane.showMessageDialog(this,
                        "Mapping file contains " + t.getNumberOfRows()
                                + " rows, tree has "
                                + phy.getNumberOfExternalNodes() + " external nodes",
                        "Warning",
                        JOptionPane.WARNING_MESSAGE);
            }


            int replaced_names = 0;
            for (int row = 0; row < t.getNumberOfRows(); ++row) {
                final String key = t.getValueAsString(0, row);
                PhylogenyNode found_node = null;
                String new_name = null;
                int count = 0;
                for (final PhylogenyNodeIterator iter = phy.iteratorExternalForward(); iter.hasNext(); ) {
                    final PhylogenyNode node = iter.next();
                    final String node_name = node.getName();
                    if (!ForesterUtil.isEmpty(node_name)) {
                        final Pattern pattern = Pattern.compile("\\b" + key + "\\b");
                        final Matcher m = pattern.matcher(node_name);
                        if (m.find()) {
                            found_node = node;
                            new_name = t.getValueAsString(1, row);
                            ++count;
                        }
                    }
                }
                if (count == 1 && found_node != null) {
                    found_node.setName(new_name);
                    ++replaced_names;
                }
            }
            if (replaced_names > 0) {
                JOptionPane.showMessageDialog(this,
                        "Successfully replaced " + replaced_names + " external nodes (out of " + phy.getNumberOfExternalNodes() + ")",
                        "Success",
                        JOptionPane.INFORMATION_MESSAGE);
            } else {
                JOptionPane.showMessageDialog(this,
                        "Failed to replace any node names",
                        "Failed to replace any node names",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
    }


    private void obtainDetailedTaxonomicInformation() {
        if (getCurrentTreePanel() != null) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ((phy != null) && !phy.isEmpty()) {
                final TaxonomyDataManager t = new TaxonomyDataManager(this,
                        _mainpanel.getCurrentTreePanel(),
                        phy.copy(),
                        false,
                        true);
                new Thread(t).start();
            }
        }
    }

    private void obtainDetailedTaxonomicInformationDelete() {
        if (getCurrentTreePanel() != null) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ((phy != null) && !phy.isEmpty()) {
                final TaxonomyDataManager t = new TaxonomyDataManager(this,
                        _mainpanel.getCurrentTreePanel(),
                        phy.copy(),
                        true,
                        true);
                new Thread(t).start();
            }
        }
    }

    private void obtainSequenceInformation() {
        if (getCurrentTreePanel() != null) {
            final Phylogeny phy = getCurrentTreePanel().getPhylogeny();
            if ((phy != null) && !phy.isEmpty()) {
                final SequenceDataRetriver u = new SequenceDataRetriver(this,
                        _mainpanel.getCurrentTreePanel(),
                        phy.copy());
                new Thread(u).start();
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
                        warnIfNotPhyloXmlValidation(getConfiguration());
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
                            } else if (parser instanceof PhyloXmlParser) {
                                warnIfNotPhyloXmlValidation(getConfiguration());
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
        for (final JFileChooser fc : new JFileChooser[] { _open_filechooser, _open_filechooser_for_species_tree,
                _values_filechooser, _annotations_filechooser }) {
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
        if (AptxUtil.canWriteFormat("tif") || AptxUtil.canWriteFormat("tiff")
                || AptxUtil.canWriteFormat("TIF")) {
            _file_jmenu.add(_write_to_tif_item = new JMenuItem("Export to TIFF file..."));
        }
        _file_jmenu.add(_write_to_png_item = new JMenuItem("Export to PNG file..."));
        _file_jmenu.add(_write_to_jpg_item = new JMenuItem("Export to JPG file..."));
        if (AptxUtil.canWriteFormat("gif")) {
            _file_jmenu.add(_write_to_gif_item = new JMenuItem("Export to GIF file..."));
        }
        if (AptxUtil.canWriteFormat("bmp")) {
            _file_jmenu.add(_write_to_bmp_item = new JMenuItem("Export to BMP file..."));
        }
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
        customizeJMenuItem(_write_to_png_item);
        customizeJMenuItem(_write_to_jpg_item);
        customizeJMenuItem(_write_to_gif_item);
        customizeJMenuItem(_write_to_tif_item);
        customizeJMenuItem(_write_to_bmp_item);
        customizeJMenuItem(_exit_item);
        _jmenubar.add(_file_jmenu);
    }

    void buildOptionsMenu() {
        _options_jmenu = MainFrame.createMenu(OPTIONS_HEADER, getConfiguration());
        _options_jmenu.addChangeListener(new ChangeListener() {

            @Override
            public void stateChanged(final ChangeEvent e) {
                MainFrame.setOvPlacementColorChooseMenuItem(_overview_placment_mi, getOptions());
                MainFrame.setTextColorChooseMenuItem(_switch_colors_mi, getCurrentTreePanel());
                MainFrame.setTextMinSupportMenuItem(_choose_minimal_confidence_mi,
                        getOptions(),
                        getCurrentTreePanel());
                MainFrame.setTextForFontChooserMenuItem(_choose_font_mi,
                        MainFrame.createCurrentFontDesc(getMainPanel()
                                .getTreeFontSet()));

                MainFrame.setTextForPdfLineWidthChooserMenuItem(_choose_pdf_width_mi, getOptions());
                MainFrame.setCycleNodeFillMenuItem(_cycle_node_fill_mi, getOptions());
                MainFrame.setCycleNodeShapeMenuItem(_cycle_node_shape_mi, getOptions());
                MainFrame.setCycleDataReturnMenuItem(_cycle_data_return, getOptions());
                MainFrame.setTextNodeSizeMenuItem(_choose_node_size_mi, getOptions());
                MainFrame.setDefaultBranchWidthMenuItem(_choose_default_branch_width_mi, getOptions());
                try {
                    getMainPanel().getControlPanel().setVisibilityOfDomainStrucureCB();
                    getMainPanel().getControlPanel().setVisibilityOfX();
                } catch (final Exception ignore) {
                    // do nothing, not important.
                }
            }
        });
        _options_jmenu.add(customizeMenuItemAsLabel(new JMenuItem(DISPLAY_SUBHEADER), getConfiguration()));
        _options_jmenu
                .add(_ext_node_dependent_cladogram_rbmi = new JRadioButtonMenuItem(MainFrame.NONUNIFORM_CLADOGRAMS_LABEL));
        _options_jmenu.add(_non_lined_up_cladograms_rbmi = new JRadioButtonMenuItem(NON_LINED_UP_CLADOGRAMS_LABEL));
        ButtonGroup _radio_group_1 = new ButtonGroup();
        _radio_group_1.add(_ext_node_dependent_cladogram_rbmi);
        _radio_group_1.add(_non_lined_up_cladograms_rbmi);
        _options_jmenu.add(_show_overview_cbmi = new JCheckBoxMenuItem(SHOW_OVERVIEW_LABEL));
        _options_jmenu.add(_show_scale_cbmi = new JCheckBoxMenuItem(DISPLAY_SCALE_LABEL));
        _options_jmenu
                .add(_show_default_node_shapes_internal_cbmi = new JCheckBoxMenuItem(DISPLAY_NODE_BOXES_LABEL_INT));
        _options_jmenu
                .add(_show_default_node_shapes_external_cbmi = new JCheckBoxMenuItem(DISPLAY_NODE_BOXES_LABEL_EXT));
        _options_jmenu
                .add(_show_default_node_shapes_for_marked_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_NODE_BOXES_LABEL_MARKED));
        _options_jmenu
                .add(_collapsed_with_average_height_cbmi = new JCheckBoxMenuItem("Proportional Height of Collapsed Subtrees"));
        _options_jmenu
                .add(_show_abbreviated_labels_for_collapsed_nodes_cbmi = new JCheckBoxMenuItem("Add Abbreviated Labels to Collapsed Subtrees"));
        _options_jmenu
                .add(_line_up_renderable_data_cbmi = new JCheckBoxMenuItem(MainFrame.LINE_UP_RENDERABLE_DATA));
        if (getConfiguration().doDisplayOption(Configuration.show_domain_architectures)) {
            _options_jmenu
                    .add(_right_line_up_domains_cbmi = new JCheckBoxMenuItem(MainFrame.RIGHT_LINE_UP_DOMAINS));
            _options_jmenu.add(_show_domain_labels = new JCheckBoxMenuItem(MainFrame.SHOW_DOMAIN_LABELS_LABEL));
        }
        _options_jmenu.add(_show_confidence_stddev_cbmi = new JCheckBoxMenuItem(SHOW_CONF_STDDEV_LABEL));
        _options_jmenu.add(_color_labels_same_as_parent_branch = new JCheckBoxMenuItem(COLOR_LABELS_LABEL));
        _color_labels_same_as_parent_branch.setToolTipText(MainFrame.COLOR_LABELS_TIP);
        _options_jmenu.add(_abbreviate_scientific_names = new JCheckBoxMenuItem(ABBREV_SN_LABEL));
        _options_jmenu.add(_label_direction_cbmi = new JCheckBoxMenuItem(LABEL_DIRECTION_LABEL));
        _label_direction_cbmi.setToolTipText(LABEL_DIRECTION_TIP);
        _options_jmenu.add(_screen_antialias_cbmi = new JCheckBoxMenuItem(SCREEN_ANTIALIAS_LABEL));
        _options_jmenu.add(_background_gradient_cbmi = new JCheckBoxMenuItem(BG_GRAD_LABEL));
        _options_jmenu.add(_cycle_node_shape_mi = new JMenuItem(MainFrame.CYCLE_NODE_SHAPE_LABEL));
        _options_jmenu.add(_cycle_node_fill_mi = new JMenuItem(MainFrame.CYCLE_NODE_FILL_LABEL));
        _options_jmenu.add(_choose_default_branch_width_mi = new JMenuItem(MainFrame.CHOOSE_BRANCH_WIDTH_LABEL));
        _options_jmenu.add(_choose_node_size_mi = new JMenuItem(MainFrame.CHOOSE_NODE_SIZE_LABEL));
        _options_jmenu.add(_choose_minimal_confidence_mi = new JMenuItem(""));
        _options_jmenu.add(_overview_placment_mi = new JMenuItem(""));
        _options_jmenu.add(_switch_colors_mi = new JMenuItem(""));
        _options_jmenu.add(_choose_font_mi = new JMenuItem(""));
        _options_jmenu.addSeparator();
        _options_jmenu.add(_cycle_data_return = new JMenuItem("Cycle Data Return"));
        _options_jmenu.addSeparator();
        _options_jmenu.add(customizeMenuItemAsLabel(new JMenuItem(SEARCH_SUBHEADER), getConfiguration()));
        // The five search-behaviour options (Match Case / Words / Regex / Inverse / Search Properties)
        // now live in the left control panel, right next to the search fields, where users can see them.
        _options_jmenu
                .add(_color_all_found_nodes_when_coloring_subtree_cbmi = new JCheckBoxMenuItem("Colorize All Found Nodes When Colorizing Subtree(s)"));
        _options_jmenu.addSeparator();
        _options_jmenu
                .add(customizeMenuItemAsLabel(new JMenuItem("Graphics Export & Printing:"), getConfiguration()));
        _options_jmenu.add(_antialias_print_cbmi = new JCheckBoxMenuItem("Antialias"));
        _options_jmenu.add(_print_black_and_white_cbmi = new JCheckBoxMenuItem("Export in Black and White"));
        _options_jmenu
                .add(_graphics_export_visible_only_cbmi = new JCheckBoxMenuItem("Limit to Visible ('Screenshot') for PNG, JPG, and GIF export"));
        _options_jmenu.add(_choose_pdf_width_mi = new JMenuItem(""));
        _options_jmenu.addSeparator();
        _options_jmenu.add(customizeMenuItemAsLabel(new JMenuItem("Newick/NHX/Nexus Read:"), getConfiguration()));
        _options_jmenu
                .add(_internal_number_are_confidence_for_nh_parsing_cbmi = new JCheckBoxMenuItem("Internal Node Names are Confidence Values"));
        _options_jmenu.add(_replace_underscores_cbmi = new JCheckBoxMenuItem("Replace Underscores with Spaces"));
        _options_jmenu
                .add(_parse_beast_style_extended_nexus_tags_cbmi = new JCheckBoxMenuItem("Parse BEAST-style extended Newick/Nexus tags"));
        _parse_beast_style_extended_nexus_tags_cbmi
                .setToolTipText("to parse elements in the form of \"[&!color=#800080]\" in Newick/Nexus formatted trees");
        _options_jmenu
                .add(_allow_errors_in_distance_to_parent_cbmi = new JCheckBoxMenuItem("Ignore Distance Values Format Errors"));
        _options_jmenu.add(_extract_taxonomy_no_rbmi = new JRadioButtonMenuItem("No Taxonomy Extraction"));
        _options_jmenu
                .add(_extract_taxonomy_pfam_strict_rbmi = new JRadioButtonMenuItem("Extract Taxonomy Codes/Ids from Pfam-style Node Names"));
        _options_jmenu
                .add(_extract_taxonomy_pfam_relaxed_rbmi = new JRadioButtonMenuItem("Extract Taxonomy Codes/Ids from Pfam-style like Node Names"));
        _options_jmenu
                .add(_extract_taxonomy_agressive_rbmi = new JRadioButtonMenuItem("Extract Taxonomy Codes/Ids/Scientific Names from Node Names"));
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
        _options_jmenu.add(customizeMenuItemAsLabel(new JMenuItem("Newick/Nexus Save:"), getConfiguration()));
        _options_jmenu
                .add(_use_brackets_for_conf_in_nh_export_cbmi = new JCheckBoxMenuItem(USE_BRACKETS_FOR_CONF_IN_NH_LABEL));
        _use_brackets_for_conf_in_nh_export_cbmi
                .setToolTipText("e.g. \"0.1[90]\" for a branch with support 90 and a length of 0.1");
        _options_jmenu
                .add(_use_internal_names_for_conf_in_nh_export_cbmi = new JCheckBoxMenuItem(USE_INTERNAL_NAMES_FOR_CONF_IN_NH_LABEL));
        customizeJMenuItem(_choose_font_mi);
        customizeJMenuItem(_choose_minimal_confidence_mi);
        customizeJMenuItem(_switch_colors_mi);
        customizeJMenuItem(_choose_pdf_width_mi);
        customizeJMenuItem(_overview_placment_mi);
        customizeCheckBoxMenuItem(_show_default_node_shapes_external_cbmi,
                getOptions().isShowDefaultNodeShapesExternal());
        customizeCheckBoxMenuItem(_show_default_node_shapes_internal_cbmi,
                getOptions().isShowDefaultNodeShapesInternal());
        customizeCheckBoxMenuItem(_show_default_node_shapes_for_marked_cbmi,
                getOptions().isShowDefaultNodeShapesForMarkedNodes());
        customizeJMenuItem(_cycle_node_shape_mi);
        customizeJMenuItem(_choose_default_branch_width_mi);
        customizeJMenuItem(_cycle_node_fill_mi);
        customizeJMenuItem(_choose_node_size_mi);
        customizeJMenuItem(_cycle_data_return);
        customizeCheckBoxMenuItem(_color_labels_same_as_parent_branch,
                getOptions().isColorLabelsSameAsParentBranch());
        customizeCheckBoxMenuItem(_screen_antialias_cbmi, getOptions().isAntialiasScreen());
        customizeCheckBoxMenuItem(_background_gradient_cbmi, getOptions().isBackgroundColorGradient());
        customizeCheckBoxMenuItem(_show_domain_labels, getOptions().isShowDomainLabels());
        customizeCheckBoxMenuItem(_abbreviate_scientific_names, getOptions().isAbbreviateScientificTaxonNames());
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
                getOptions().isReplaceUnderscoresInNhParsing());
        customizeCheckBoxMenuItem(_color_all_found_nodes_when_coloring_subtree_cbmi,
                getOptions().isColorAllFoundNodesWhenColoringSubtree());
        customizeCheckBoxMenuItem(_parse_beast_style_extended_nexus_tags_cbmi,
                getOptions().isParseBeastStyleExtendedNexusTags());
        customizeCheckBoxMenuItem(_graphics_export_visible_only_cbmi, getOptions().isGraphicsExportVisibleOnly());
        customizeCheckBoxMenuItem(_show_confidence_stddev_cbmi, getOptions().isShowConfidenceStddev());
        customizeCheckBoxMenuItem(_use_brackets_for_conf_in_nh_export_cbmi,
                getOptions()
                        .getNhConversionSupportValueStyle() == NH_CONVERSION_SUPPORT_VALUE_STYLE.IN_SQUARE_BRACKETS);
        customizeCheckBoxMenuItem(_use_internal_names_for_conf_in_nh_export_cbmi,
                getOptions()
                        .getNhConversionSupportValueStyle() == NH_CONVERSION_SUPPORT_VALUE_STYLE.AS_INTERNAL_NODE_NAMES);
        customizeCheckBoxMenuItem(_line_up_renderable_data_cbmi, getOptions().isLineUpRendarableNodeData());
        customizeCheckBoxMenuItem(_right_line_up_domains_cbmi, getOptions().isRightLineUpDomains());
        // _options_jmenu is not added to the menu bar; its items are folded into the Settings dialog.
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
        _tools_menu.setToolTipText("Colorize branches and subtrees, and replace node names");
        if (getConfiguration().isEditable()) {
            _tools_menu.add(_replace_names_item = new JMenuItem("Replace Node Names"));
            _replace_names_item.setToolTipText("to replace external node names using a tab separated mapping file");
            _tools_menu.addSeparator();
            customizeJMenuItem(_replace_names_item);
        }
        _tools_menu.add(_confcolor_item = new JMenuItem("Colorize Branches Depending on Confidence"));
        customizeJMenuItem(_confcolor_item);
        _tools_menu.add(_color_rank_jmi = new JMenuItem("Colorize Subtrees via Taxonomic Rank"));
        customizeJMenuItem(_color_rank_jmi);
        _color_rank_jmi.setToolTipText("for example, at \"Class\" level, colorize mammal specific subtree red");
        _tools_menu.add(_taxcolor_item = new JMenuItem("Taxonomy Colorize Branches"));
        customizeJMenuItem(_taxcolor_item);
        _tools_menu.addSeparator();
        _tools_menu.add(_remove_visual_styles_item = new JMenuItem("Delete All Visual Styles From Nodes"));
        _remove_visual_styles_item
                .setToolTipText("To remove all node visual styles (fonts, colors) from the current phylogeny");
        customizeJMenuItem(_remove_visual_styles_item);
        _tools_menu.add(_remove_branch_color_item = new JMenuItem("Delete All Colors From Branches"));
        _remove_branch_color_item.setToolTipText("To remove all branch color values from the current phylogeny");
        customizeJMenuItem(_remove_branch_color_item);
        _tools_menu.addSeparator();
        _tools_menu.add(_mad_root_item = new JMenuItem("MAD-Root"));
        _mad_root_item.setToolTipText("Root by Minimal Ancestor Deviation (Tria et al. 2017); requires branch lengths");
        customizeJMenuItem(_mad_root_item);
        _tools_menu.add(_midpoint_root_item = new JMenuItem("Midpoint-Root"));
        customizeJMenuItem(_midpoint_root_item);
        _tools_menu.addSeparator();
        _tools_menu.add(_delete_selected_nodes_item = new JMenuItem("Delete Selected Nodes"));
        _delete_selected_nodes_item.setToolTipText("To delete all selected external nodes");
        customizeJMenuItem(_delete_selected_nodes_item);
        _tools_menu.add(_delete_not_selected_nodes_item = new JMenuItem("Retain Selected Nodes"));
        _delete_not_selected_nodes_item.setToolTipText("To delete all not selected external nodes");
        customizeJMenuItem(_delete_not_selected_nodes_item);
        _tools_menu.addSeparator();
        _tools_menu.add(_collapse_species_specific_subtrees = new JMenuItem("Collapse Single Taxonomy-Subtrees"));
        customizeJMenuItem(_collapse_species_specific_subtrees);
        _collapse_species_specific_subtrees
                .setToolTipText("To (reversibly) collapse subtrees associated with only one taxonomy (such as species specific subtrees)");
        _tools_menu
                .add(_collapse_below_threshold = new JMenuItem("Collapse Branches with Confidence Below Threshold into Multifurcations"));
        customizeJMenuItem(_collapse_below_threshold);
        _collapse_below_threshold
                .setToolTipText("To (permanently) collapse branches with confidence values below a threshold into multifurcations (in the case of multiple confidences per branch: without at least one confidence value above a threshold)");
        //
        _tools_menu
                .add(_collapse_below_branch_length = new JMenuItem("Collapse Branches with Branch Lengths Below Threshold into Multifurcations"));
        customizeJMenuItem(_collapse_below_branch_length);
        _collapse_below_branch_length
                .setToolTipText("To (permanently) collapse branches with branches with branch lengths below a threshold into multifurcations");
        //
        _tools_menu.addSeparator();
        _tools_menu.add(_obtain_seq_information_jmi = new JMenuItem("Obtain Sequence Information"));
        customizeJMenuItem(_obtain_seq_information_jmi);
        _obtain_seq_information_jmi.setToolTipText("To add additional sequence information");
        _tools_menu
                .add(_obtain_detailed_taxonomic_information_jmi = new JMenuItem(OBTAIN_DETAILED_TAXONOMIC_INFORMATION));
        customizeJMenuItem(_obtain_detailed_taxonomic_information_jmi);
        _obtain_detailed_taxonomic_information_jmi
                .setToolTipText("To add additional taxonomic information (from UniProt Taxonomy)");
        _tools_menu
                .add(_obtain_detailed_taxonomic_information_deleting_jmi = new JMenuItem("Obtain Detailed Taxonomic Information (deletes nodes!)"));
        customizeJMenuItem(_obtain_detailed_taxonomic_information_deleting_jmi);
        _obtain_detailed_taxonomic_information_deleting_jmi
                .setToolTipText("To add additional taxonomic information, deletes nodes for which taxonomy cannot found (from UniProt Taxonomy)");
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

    static MainFrame createInstance(final Phylogeny[] phys, final String config_file_name, final String title) {
        return new MainFrameApplication(phys, config_file_name, title);
    }

    static void warnIfNotPhyloXmlValidation(final Configuration c) {
        if (!c.isValidatePhyloXmlAgainstSchema()) {
            JOptionPane
                    .showMessageDialog(null,
                            ForesterUtil
                                    .wordWrap("phyloXML XSD-based validation is turned off [enable with line 'validate_against_phyloxml_xsd_schem: true' in configuration file]",
                                            80),
                            "Warning",
                            JOptionPane.WARNING_MESSAGE);
        }
    }
} // MainFrameApplication.
