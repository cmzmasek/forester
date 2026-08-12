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
import java.awt.Component;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
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

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButton;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.tools.NodeDataExporter;
import org.forester.archaeopteryx.tools.RepresentativeTipSelector;
import org.forester.archaeopteryx.tools.SequenceAndTaxonomyDataObtainer;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
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
    private JMenuItem _select_representative_tips_jmi;
    // Others:
    double _min_not_collapse = AptxConstants.MIN_NOT_COLLAPSE_DEFAULT;
    double _min_not_collapse_bl = 0.001;
    // "Select Representative Tips" remembered inputs (this session)
    private double _repsel_cutoff = 0.05;
    private int _repsel_target = 100;
    private boolean _repsel_by_cutoff = true;
    private int _repsel_pick_index = 0;

    private MainFrameApplication(final Phylogeny[] phys, final Configuration config) {
        _configuration = config;
        if (_configuration == null) {
            throw new IllegalArgumentException("configuration is null");
        }
        setVisible(false);
        setOptions(optionsWithSavedPreferences());
        registerMacOsQuitHandler();
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
            _startup_dir = current_dir; // a transient default for dialogs; NOT persisted (only user choices are)
        }
        // restore each dialog's last-used directory from the previous session (a saved dir that no longer exists is
        // ignored later, lazily, by getCurrentDir -- which then falls back to the startup dir / home).
        new DirectoryPreferences().applyTo(_current_dirs);
        // hide until everything is ready
        setVisible(false);
        setOptions(optionsWithSavedPreferences());
        registerMacOsQuitHandler();
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
        buildEditMenu();
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
            } else if (o == _select_representative_tips_jmi) {
                selectRepresentativeTips();
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
        final JMenuItem edit = new JMenuItem(EDIT_TREE_INFO_LABEL);
        edit.addActionListener(ae -> {
            final JTabbedPane tabs = getMainPanel().getTabbedPane();
            if ((tabs == null) || (index < 0) || (index >= tabs.getTabCount())) {
                return; // the clicked tab was closed/moved since the popup opened -- don't edit the wrong tree
            }
            tabs.setSelectedIndex(index); // act on the clicked tab, then edit its (now current) tree
            showTreeInfoDialog();
        });
        popup.add(edit);
        popup.addSeparator();
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
                "Collapse " + candidates.size() + " branch(es) with confidence below " + threshold
                        + " into polytomies?",
                "Collapse " + candidates.size() + " branch(es)?",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.QUESTION_MESSAGE);
        if (confirm != JOptionPane.OK_OPTION) {
            return;
        }
        getCurrentTreePanel().pushUndoCheckpoint("Collapse Branches");
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
                "Collapse " + candidates.size() + " branch(es) shorter than " + threshold + " into polytomies?",
                "Collapse " + candidates.size() + " branch(es)?",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.QUESTION_MESSAGE);
        if (confirm != JOptionPane.OK_OPTION) {
            return;
        }
        getCurrentTreePanel().pushUndoCheckpoint("Collapse Branches");
        for (final PhylogenyNode n : candidates) {
            PhylogenyMethods.removeNode(n, phy);
        }
        refreshAfterBranchCollapse(phy);
        JOptionPane.showMessageDialog(this,
                "Collapsed " + candidates.size() + " branch(es) shorter than " + threshold + ".",
                "Collapsed " + candidates.size() + " branch(es)",
                JOptionPane.INFORMATION_MESSAGE);
    }

    private void selectRepresentativeTips() {
        if (getCurrentTreePanel() == null) {
            return;
        }
        final TreePanel tp = getCurrentTreePanel();
        final Phylogeny phy = tp.getPhylogeny();
        if ((phy == null) || phy.isEmpty() || (phy.getNumberOfExternalNodes() < 3)) {
            JOptionPane.showMessageDialog(this,
                    "The tree needs at least three external nodes to select representatives.",
                    "Too Few Tips",
                    JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        // the external tips the user has selected (search a/b + manual clicks; branch-click a clade to protect
        // its tips) can be protected from removal -- captured now, before the run overwrites the found highlight
        final Set<Long> protected_ids = new HashSet<>();
        for (final PhylogenyNode t : NodeDataExporter.selectedExternalTips(phy,
                tp.getFoundNodesAsListOfPhylogenyNodes())) {
            protected_ids.add(t.getId());
        }
        // --- input dialog (distance-cutoff mode is unavailable without branch lengths) ---
        final boolean has_bl = RepresentativeTipSelector.hasUsableBranchLengths(phy);
        final RepSelInputs in = buildRepSelInputs(has_bl, _repsel_by_cutoff, _repsel_cutoff, _repsel_target,
                _repsel_pick_index, protected_ids.size());
        final JOptionPane pane = new JOptionPane(in._panel, JOptionPane.PLAIN_MESSAGE, JOptionPane.OK_CANCEL_OPTION);
        final JDialog dialog = pane.createDialog(this, "Select Representative Tips");
        dialog.setResizable(true);
        // The plain showConfirmDialog is fixed-size and centered on the main window, so on a shorter/scaled
        // display (or a larger UI font) its bottom -- the OK/Cancel row -- can be clipped or pushed under the
        // Dock. Pack once to measure the real chrome (title bar + button row), and if the dialog exceeds the
        // usable screen, scroll-cap just the form so the whole dialog fits; then clamp it fully on-screen.
        final Rectangle usable = GraphicsEnvironment.getLocalGraphicsEnvironment().getMaximumWindowBounds();
        if (dialog.getHeight() > usable.height) {
            final int chrome = dialog.getHeight() - in._panel.getHeight(); // measured title + buttons + insets
            final Component fitted = fitToScreen(in._panel, usable.height - chrome);
            if (fitted != in._panel) {
                pane.setMessage(fitted);
                dialog.pack();
            }
        }
        dialog.setLocation(clampFullyOnScreen(dialog.getLocation(), dialog.getSize(), usable));
        dialog.setVisible(true);
        dialog.dispose();
        final Object pane_value = pane.getValue();
        if (!(pane_value instanceof Integer) || (((Integer) pane_value).intValue() != JOptionPane.OK_OPTION)) {
            return;
        }
        // --- parse & validate ---
        final boolean by_cutoff = in._cutoff_rb.isSelected();
        final String s = in._value_tf.getText().trim();
        double cutoff = 0.0;
        int target = 0;
        if (by_cutoff) {
            try {
                cutoff = Double.parseDouble(s);
            } catch (final NumberFormatException ex) {
                JOptionPane.showMessageDialog(this, "Please enter a numeric distance cutoff (for example 0.05).",
                        "Invalid Cutoff", JOptionPane.WARNING_MESSAGE);
                return;
            }
            if (Double.isNaN(cutoff) || Double.isInfinite(cutoff) || (cutoff <= 0.0)) {
                JOptionPane.showMessageDialog(this, "The distance cutoff must be a positive number.",
                        "Invalid Cutoff", JOptionPane.WARNING_MESSAGE);
                return;
            }
        }
        else {
            try {
                target = Integer.parseInt(s);
            } catch (final NumberFormatException ex) {
                JOptionPane.showMessageDialog(this,
                        "Please enter a whole number of representatives (for example 100).",
                        "Invalid Target", JOptionPane.WARNING_MESSAGE);
                return;
            }
            if (target < 1) {
                JOptionPane.showMessageDialog(this, "The target number of representatives must be at least 1.",
                        "Invalid Target", JOptionPane.WARNING_MESSAGE);
                return;
            }
        }
        // persist the inputs for next time; keep the remembered mode when cutoff was forced off by a
        // branch-length-less tree, so it does not overwrite the user's real preference
        if (has_bl) {
            _repsel_by_cutoff = by_cutoff;
        }
        _repsel_pick_index = in._pick_box.getSelectedIndex();
        if (by_cutoff) {
            _repsel_cutoff = cutoff;
        }
        else {
            _repsel_target = target;
        }
        final RepresentativeTipSelector.RepresentativePick pick = (_repsel_pick_index == 1)
                ? RepresentativeTipSelector.RepresentativePick.LONGEST_BRANCH
                : RepresentativeTipSelector.RepresentativePick.MEDOID;
        final Set<Long> protect = in._protect_cb.isSelected() ? protected_ids : new HashSet<>();
        // --- compute on the displayed tree & highlight ---
        final RepresentativeTipSelector.SelectionResult result = by_cutoff
                ? RepresentativeTipSelector.selectByCutoff(phy, cutoff, pick, protect)
                : RepresentativeTipSelector.selectByTargetCount(phy, target, pick, protect);
        tp.setFoundNodes0(new HashSet<>(result.representativeIds()));
        _mainpanel.getControlPanel().displayedPhylogenyMightHaveChanged(true);
        // --- report & offer extraction into a new tab ---
        final int kept = result.getKeptCount();
        final int choice = JOptionPane.showConfirmDialog(this,
                result.summary() + "\n\nCreate a new tab containing only these " + kept
                        + (kept == 1 ? " tip?" : " tips?"),
                "Representative Tips Selected", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
        if (choice != JOptionPane.YES_OPTION) {
            return;
        }
        // Build the extracted tree by index-mapping the kept tips onto a copy (copy() preserves external-node
        // order), so we never depend on node ids matching across the copy.
        final Set<Long> keep_ids = result.representativeIds();
        final List<PhylogenyNode> orig_ext = phy.getExternalNodes();
        final Set<Integer> keep_indices = new HashSet<>();
        for (int i = 0; i < orig_ext.size(); ++i) {
            if (keep_ids.contains(orig_ext.get(i).getId())) {
                keep_indices.add(i);
            }
        }
        final Phylogeny copy = phy.copy();
        final List<PhylogenyNode> copy_ext = copy.getExternalNodes();
        final Set<Long> to_delete = new HashSet<>();
        for (int i = 0; i < copy_ext.size(); ++i) {
            if (!keep_indices.contains(i)) {
                to_delete.add(copy_ext.get(i).getId());
            }
        }
        PhylogenyMethods.deleteExternalNodesNegativeSelection(to_delete, copy);
        // Base the extracted tree's name on the parent's DISPLAYED name (its tab title, e.g. "mammals"),
        // which is what the user sees; phy.getName() is often empty for loaded trees (the tab title then
        // comes from the file name).
        final int parent_tab = _mainpanel.getTabbedPane().getSelectedIndex();
        final String parent_title = (parent_tab >= 0) ? _mainpanel.getTabbedPane().getTitleAt(parent_tab) : null;
        final String parent_display = resolveParentTreeName(phy.getName(), parent_title);
        final int n_reps = copy.getNumberOfExternalNodes();
        copy.setName(representativeTreeName(parent_display, n_reps));
        // Record how this tree was produced in its description (provenance), preserving any inherited one.
        final String provenance = representativeTreeDescription(by_cutoff, cutoff, target, pick, n_reps,
                stripShortExtension(parent_display), phy.getNumberOfExternalNodes());
        final String existing_desc = copy.getDescription();
        copy.setDescription(ForesterUtil.isEmpty(existing_desc) ? provenance : existing_desc + " " + provenance);
        addDerivedPhylogenyInNewTab(copy);
        showWhole();
        getCurrentTreePanel().setEdited(true);
        JOptionPane.showMessageDialog(this,
                "Created a new tab with " + copy.getNumberOfExternalNodes()
                        + (copy.getNumberOfExternalNodes() == 1 ? " tip." : " tips."),
                "Representatives Extracted", JOptionPane.INFORMATION_MESSAGE);
    }

    /**
     * Opens a tree derived from the currently displayed one (e.g. the representative-tips extraction) in a new
     * tab, having the new tab inherit the source tab's phylogram/cladogram display type. Without this the new
     * tab falls back to the global config default (see {@link ControlPanel#phylogenyAdded(Configuration)}),
     * which shows the derived tree of a phylogram parent as a cladogram -- confusing, since the derived tree
     * keeps the parent's branch lengths. A phylogram type is only honored when the derived tree actually carries
     * branch lengths (otherwise a phylogram would be degenerate, so the cladogram default is kept). Must be
     * called while the source tab is still the selected one. Package-visible for testing.
     */
    void addDerivedPhylogenyInNewTab(final Phylogeny derived) {
        final Options.PHYLOGENY_DISPLAY_TYPE source_type = _mainpanel.getControlPanel().getTreeDisplayType();
        _mainpanel.addPhylogenyInNewTab(derived, getConfiguration(), derived.getName(), null);
        final boolean source_is_phylogram = (source_type == Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM)
                || (source_type == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM);
        if (!source_is_phylogram || ((_mainpanel.getCurrentTreePanel() != null)
                && _mainpanel.getCurrentTreePanel().isPhyHasBranchLengths())) {
            _mainpanel.getControlPanel().setTreeDisplayType(source_type);
        }
    }

    /** The input controls for the "Select Representative Tips" dialog, exposed so they can be unit-tested. */
    static final class RepSelInputs {

        final JPanel _panel;
        final JRadioButton _cutoff_rb;
        final JRadioButton _target_rb;
        final JTextField _value_tf;
        final JComboBox<String> _pick_box;
        final JCheckBox _protect_cb;

        RepSelInputs(final JPanel panel, final JRadioButton cutoff_rb, final JRadioButton target_rb,
                     final JTextField value_tf, final JComboBox<String> pick_box, final JCheckBox protect_cb) {
            _panel = panel;
            _cutoff_rb = cutoff_rb;
            _target_rb = target_rb;
            _value_tf = value_tf;
            _pick_box = pick_box;
            _protect_cb = protect_cb;
        }
    }

    /**
     * Builds the "Select Representative Tips" input panel. When {@code has_branch_lengths} is false the
     * distance-cutoff option is disabled and the target-count option is preselected (a note explains why),
     * because patristic distance — and therefore a distance cutoff — is meaningless without branch lengths.
     */
    static RepSelInputs buildRepSelInputs(final boolean has_branch_lengths, final boolean prefer_cutoff,
                                          final double cutoff_default, final int target_default,
                                          final int pick_index, final int selection_count) {
        final boolean start_cutoff = has_branch_lengths && prefer_cutoff;
        final JRadioButton cutoff_rb = new JRadioButton("By distance cutoff", start_cutoff);
        final JRadioButton target_rb = new JRadioButton("By target number of representatives", !start_cutoff);
        final ButtonGroup bg = new ButtonGroup();
        bg.add(cutoff_rb);
        bg.add(target_rb);
        if (!has_branch_lengths) {
            cutoff_rb.setEnabled(false);
        }
        final JTextField value_tf = new JTextField(8);
        value_tf.setText(start_cutoff ? String.valueOf(cutoff_default) : String.valueOf(target_default));
        final JLabel hint = new JLabel();
        final Runnable sync_hint = () -> hint.setText(cutoff_rb.isSelected()
                ? "<html><i>Maximum distance between any two tips kept in the same group.</i></html>"
                : "<html><i>Approximate number of representative tips to keep.</i></html>");
        cutoff_rb.addItemListener(e -> {
            if (cutoff_rb.isSelected()) {
                value_tf.setText(String.valueOf(cutoff_default));
                sync_hint.run();
            }
        });
        target_rb.addItemListener(e -> {
            if (target_rb.isSelected()) {
                value_tf.setText(String.valueOf(target_default));
                sync_hint.run();
            }
        });
        sync_hint.run();
        final JComboBox<String> pick_box = new JComboBox<>(
                new String[] { "Most central (medoid)", "Most divergent (longest branch)" });
        pick_box.setSelectedIndex(pick_index);
        final JCheckBox protect_cb;
        if (selection_count > 0) {
            protect_cb = new JCheckBox("Keep the " + selection_count + " selected/found "
                    + (selection_count == 1 ? "tip" : "tips") + " (never drop them)", true);
        }
        else {
            protect_cb = new JCheckBox("Keep selected/found tips (none are selected)");
            protect_cb.setEnabled(false);
        }
        // A vertical BoxLayout (not GridLayout) so each row keeps its own natural height: a GridLayout(0,1)
        // stretches every row -- labels and blank spacers included -- to the tallest child (the combo box),
        // inflating the dialog by ~100px and, on a smaller/scaled display, pushing the OK/Cancel row off the
        // bottom. Small struts give the group spacing the old blank-label rows provided, far more compactly.
        final JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        value_tf.setMaximumSize(new Dimension(Integer.MAX_VALUE, value_tf.getPreferredSize().height));
        pick_box.setMaximumSize(new Dimension(Integer.MAX_VALUE, pick_box.getPreferredSize().height));
        addLeft(panel, new JLabel("Select one representative per group of similar tips:"));
        addLeft(panel, cutoff_rb);
        addLeft(panel, target_rb);
        if (!has_branch_lengths) {
            addLeft(panel, new JLabel(
                    "<html><i>This tree has no branch lengths, so only \"target number\" is available.</i></html>"));
        }
        panel.add(Box.createVerticalStrut(10));
        addLeft(panel, new JLabel("Value:"));
        addLeft(panel, value_tf);
        addLeft(panel, hint);
        panel.add(Box.createVerticalStrut(10));
        addLeft(panel, new JLabel("Representative:"));
        addLeft(panel, pick_box);
        panel.add(Box.createVerticalStrut(10));
        addLeft(panel, protect_cb);
        return new RepSelInputs(panel, cutoff_rb, target_rb, value_tf, pick_box, protect_cb);
    }

    /** Adds a component to a vertical-BoxLayout panel, left-aligned so the column does not zig-zag. */
    private static void addLeft(final JPanel panel, final JComponent c) {
        c.setAlignmentX(Component.LEFT_ALIGNMENT);
        panel.add(c);
    }

    /**
     * Caps a dialog's message component to a maximum height. When {@code form}'s preferred height exceeds
     * {@code max_form_height} it is wrapped in a borderless vertical scroll pane sized to that height so the
     * surrounding dialog's OK/Cancel row cannot be pushed off-screen; otherwise the form is returned
     * unchanged. The caller supplies the available height (usable screen minus the dialog's actual title +
     * button-row chrome, measured after a first pack), so this stays pure and headless-safe for unit tests.
     */
    static Component fitToScreen(final JComponent form, final int max_form_height) {
        final Dimension pref = form.getPreferredSize();
        if (pref.height <= max_form_height) {
            return form;
        }
        final JScrollPane scroller = new JScrollPane(form, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scroller.setBorder(BorderFactory.createEmptyBorder());
        scroller.getViewport().setPreferredSize(new Dimension(pref.width + 24, Math.max(120, max_form_height)));
        return scroller;
    }

    /**
     * Clamps a window location so a window of {@code size} sits fully inside {@code usable} (the screen area
     * excluding the menu bar and Dock). Prevents a dialog centered on a large parent from spilling its bottom
     * (the OK/Cancel row) under the Dock or off the screen edge. If the window is larger than the usable area
     * in a dimension it is pinned to that edge (top/left) so at least the start is reachable. Pure/testable.
     */
    static Point clampFullyOnScreen(final Point location, final Dimension size, final Rectangle usable) {
        final int max_x = Math.max(usable.x, (usable.x + usable.width) - size.width);
        final int max_y = Math.max(usable.y, (usable.y + usable.height) - size.height);
        final int x = Math.min(Math.max(location.x, usable.x), max_x);
        final int y = Math.min(Math.max(location.y, usable.y), max_y);
        return new Point(x, y);
    }

    /**
     * The parent tree's display name for {@link #representativeTreeName}: its tab title (what the user sees,
     * typically the file name for a loaded tree), preferred over the often-empty internal {@code phy.getName()}.
     * The {@code [N]} placeholder that an unnamed tree's tab gets is treated as no name. Returns {@code null}
     * when neither is usable (the caller then falls back to {@code tree}).
     */
    static String resolveParentTreeName(final String phylogeny_name, final String tab_title) {
        if (!ForesterUtil.isEmpty(tab_title) && !tab_title.matches("\\[\\d+\\]")) {
            return tab_title;
        }
        if (!ForesterUtil.isEmpty(phylogeny_name)) {
            return phylogeny_name;
        }
        return null;
    }

    /**
     * Strips a trailing file-type suffix -- a dot plus 1 to 5 non-dot characters, e.g. {@code .xml},
     * {@code .nwk}, {@code .nexus} -- from a tree name, so the extracted tree reads {@code mammals_233reps}
     * rather than {@code mammals.xml_233reps}. Only the last such suffix is removed; a name without one (or
     * with a longer suffix) is returned unchanged. {@code null} in, {@code null} out.
     */
    static String stripShortExtension(final String name) {
        if (name == null) {
            return null;
        }
        return name.replaceFirst("\\.[^.]{1,5}$", "");
    }

    /**
     * Name for the tree extracted by "Select Representative Tips": the parent tree's name (with any file
     * extension stripped) followed by the representative count, e.g. {@code mammals_233reps} (singular
     * {@code _1rep}). Underscore-joined so it stays clean as an export filename; falls back to {@code tree}
     * when the parent has no name.
     */
    static String representativeTreeName(final String parent_name, final int n_representatives) {
        final String stripped = stripShortExtension(parent_name);
        final String base = ForesterUtil.isEmpty(stripped) ? "tree" : stripped;
        return base + "_" + n_representatives + (n_representatives == 1 ? "rep" : "reps");
    }

    /**
     * Provenance sentence stored in the extracted tree's description, e.g. <i>Used the distance-cutoff
     * (maximum distance 0.05, medoid representative) algorithm to select 233 representative tips from tree
     * named "mammals" with 1000 tips.</i> Records the selection mode + its value, the representative-pick
     * rule, how many tips were kept, and the parent tree's name and size.
     */
    static String representativeTreeDescription(final boolean by_cutoff, final double cutoff, final int target,
                                                final RepresentativeTipSelector.RepresentativePick pick,
                                                final int n_representatives, final String parent_name,
                                                final int parent_tip_count) {
        final String pick_desc = (pick == RepresentativeTipSelector.RepresentativePick.LONGEST_BRANCH)
                ? "longest-branch" : "medoid";
        final String algorithm = by_cutoff
                ? "distance-cutoff (maximum distance " + cutoff + ", " + pick_desc + " representative)"
                : "target-count (target " + target + ", " + pick_desc + " representative)";
        final String base = ForesterUtil.isEmpty(parent_name) ? "tree" : parent_name;
        return "Used the " + algorithm + " algorithm to select " + n_representatives + " representative "
                + (n_representatives == 1 ? "tip" : "tips") + " from tree named \"" + base + "\" with "
                + parent_tip_count + (parent_tip_count == 1 ? " tip." : " tips.");
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
        final File my_dir = getCurrentDir(DirectoryPreferences.Category.OPEN);
        // Open file-open dialog and set current directory
        if (my_dir != null) {
            _open_filechooser.setCurrentDirectory(my_dir);
        }
        final int result = _open_filechooser.showOpenDialog(_contentpane);
        // All done: get the file
        final File[] files = _open_filechooser.getSelectedFiles();
        setCurrentDir(DirectoryPreferences.Category.OPEN, _open_filechooser.getCurrentDirectory());
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
                        if (nhx_or_nexus) {
                            offerInternalNamesAsConfidence(phys);
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
        final File my_dir = getCurrentDir(DirectoryPreferences.Category.OPEN);
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
        _file_jmenu.add(_copy_image_to_clipboard_item = new JMenuItem("Copy Image to Clipboard"));
        _copy_image_to_clipboard_item
                .setToolTipText("Copy the current tree as an image, ready to paste into a document or slide");
        _file_jmenu.addSeparator();
        // Most users never read the manual, so the menu tooltips advertise how to narrow the export.
        final String scope_hint = "<br><i>Use sub-tree display, search, and/or node selection to restrict "
                + "the export.</i>";
        _file_jmenu.add(_export_seqs_fasta_item = new JMenuItem("Export Sequences (FASTA)..."));
        _export_seqs_fasta_item.setToolTipText("<html>Write the tip molecular sequences to a FASTA file."
                + scope_hint + "</html>");
        _file_jmenu.add(_export_node_data_item = new JMenuItem("Export Node Data (TSV)..."));
        _export_node_data_item.setToolTipText("<html>Write the tip data (names, taxonomy, sequence, branch "
                + "length, properties) to a tab-separated table." + scope_hint + "</html>");
        _file_jmenu.add(_import_annotations_item = new JMenuItem("Import Annotations (CSV/TSV)..."));
        _import_annotations_item.setToolTipText("<html>Read a CSV or TSV table and join its columns onto the tips: "
                + "match a chosen key column against the tip name, sequence accession, or taxonomy, with a preview of "
                + "the match before committing.<br><i>Recognized columns fill taxonomy/sequence fields; any other "
                + "column becomes a node property you can color by or show as an annotation column.</i></html>");
        _file_jmenu.add(_import_annotations_url_item = new JMenuItem("Import Annotations from URL..."));
        _import_annotations_url_item.setToolTipText("<html>Fetch a CSV/TSV from a URL (e.g. a Google Sheet published "
                + "to the web as CSV) and run the same import dialog.</html>");
        _file_jmenu.add(_reimport_annotations_item = new JMenuItem("Re-import Annotations"));
        _reimport_annotations_item.setToolTipText("<html>Re-fetch this tree's last annotation source (file or URL) and "
                + "re-apply the same key column, match attribute, and column mapping with one click.<br><i>Edit your "
                + "sheet/file, then pull the changes in without walking the dialog again.</i></html>");
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
        customizeJMenuItem(_copy_image_to_clipboard_item);
        customizeJMenuItem(_write_to_tif_item);
        customizeJMenuItem(_export_seqs_fasta_item);
        customizeJMenuItem(_export_node_data_item);
        customizeJMenuItem(_import_annotations_item);
        customizeJMenuItem(_import_annotations_url_item);
        customizeJMenuItem(_reimport_annotations_item);
        customizeJMenuItem(_exit_item);
        // Keyboard shortcuts for the common File actions. The platform menu-shortcut key is Cmd on macOS, Ctrl
        // elsewhere. Copy Image uses Shift too, deliberately NOT plain Cmd-C, so it can't hijack text copy in the
        // search field. Zoom already has +/- keys (TreePanel), so it needs none here.
        final int shortcut = Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx();
        _open_item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, shortcut));
        _save_item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, shortcut));
        _close_item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_W, shortcut));
        _copy_image_to_clipboard_item
                .setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C, shortcut | InputEvent.SHIFT_DOWN_MASK));
        if (_new_item != null) {
            _new_item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, shortcut));
        }
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
        _show_tree_name_cbmi = new JCheckBoxMenuItem(SHOW_TREE_NAME_LABEL);
        _show_tree_name_cbmi.setToolTipText("show the tree's name on the canvas (lower-left, or lower-right when the scale is shown)");
        _show_scale_grid_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_SCALE_GRID_LABEL);
        _show_scale_grid_cbmi.setToolTipText(MainFrame.DISPLAY_SCALE_GRID_TIP);
        _show_scale_axis_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_SCALE_AXIS_LABEL);
        _show_scale_axis_cbmi.setToolTipText(MainFrame.DISPLAY_SCALE_AXIS_TIP);
        _show_hpd_bars_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_HPD_BARS_LABEL);
        _show_hpd_bars_cbmi.setToolTipText(MainFrame.DISPLAY_HPD_BARS_TIP);
        _show_zebra_stripes_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_ZEBRA_STRIPES_LABEL);
        _show_zebra_stripes_cbmi.setToolTipText(MainFrame.DISPLAY_ZEBRA_STRIPES_TIP);
        _show_internal_taxonomy_key_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_INTERNAL_TAXONOMY_KEY_LABEL);
        _show_internal_taxonomy_key_cbmi.setToolTipText(MainFrame.DISPLAY_INTERNAL_TAXONOMY_KEY_TIP);
        _tip_labels_below_columns_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_TIP_LABELS_BELOW_COLUMNS_LABEL);
        _tip_labels_below_columns_cbmi.setToolTipText(MainFrame.DISPLAY_TIP_LABELS_BELOW_COLUMNS_TIP);
        _reverse_tip_order_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_REVERSE_TIP_ORDER_LABEL);
        _reverse_tip_order_cbmi.setToolTipText(MainFrame.DISPLAY_REVERSE_TIP_ORDER_TIP);
        _bold_found_labels_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_BOLD_FOUND_LABELS_LABEL);
        _bold_found_labels_cbmi.setToolTipText(MainFrame.DISPLAY_BOLD_FOUND_LABELS_TIP);
        _dim_non_matches_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_DIM_NON_MATCHES_LABEL);
        _dim_non_matches_cbmi.setToolTipText(MainFrame.DISPLAY_DIM_NON_MATCHES_TIP);
        _pulse_found_nodes_cbmi = new JCheckBoxMenuItem(MainFrame.DISPLAY_PULSE_FOUND_NODES_LABEL);
        _pulse_found_nodes_cbmi.setToolTipText(MainFrame.DISPLAY_PULSE_FOUND_NODES_TIP);
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
        _transparent_export_background_cbmi = new JCheckBoxMenuItem(MainFrame.TRANSPARENT_BG_LABEL);
        _transparent_export_background_cbmi.setToolTipText(MainFrame.TRANSPARENT_BG_TIP);
        _graphics_export_white_background_cbmi = new JCheckBoxMenuItem(MainFrame.WHITE_BG_LABEL);
        _graphics_export_white_background_cbmi.setToolTipText(MainFrame.WHITE_BG_TIP);
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
        // The "Taxonomy Extraction from Node Names" GUI controls were retired (the TAXONOMY_EXTRACTION
        // enum + ParserUtils stay for the CLI/library and config); the GUI now reads with no extraction.
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
        customizeCheckBoxMenuItem(_transparent_export_background_cbmi,
                getOptions().isTransparentExportBackground());
        customizeCheckBoxMenuItem(_graphics_export_white_background_cbmi,
                getOptions().isGraphicsExportWhiteBackground());
        customizeCheckBoxMenuItem(_show_scale_cbmi, getOptions().isShowScale());
        customizeCheckBoxMenuItem(_show_tree_name_cbmi, getOptions().isShowTreeName());
        customizeCheckBoxMenuItem(_show_scale_grid_cbmi, getOptions().isShowScaleGrid());
        customizeCheckBoxMenuItem(_show_scale_axis_cbmi, getOptions().isShowScaleAxis());
        customizeCheckBoxMenuItem(_show_hpd_bars_cbmi, getOptions().isShowHpdBars());
        customizeCheckBoxMenuItem(_show_zebra_stripes_cbmi, getOptions().isShowZebraStripes());
        customizeCheckBoxMenuItem(_show_internal_taxonomy_key_cbmi, getOptions().isShowInternalTaxonomyKey());
        customizeCheckBoxMenuItem(_tip_labels_below_columns_cbmi, getOptions().isTipLabelsBelowColumns());
        customizeCheckBoxMenuItem(_reverse_tip_order_cbmi, getOptions().isReverseTipOrder());
        customizeCheckBoxMenuItem(_bold_found_labels_cbmi, getOptions().isBoldFoundLabels());
        customizeCheckBoxMenuItem(_dim_non_matches_cbmi, getOptions().isDimNonMatches());
        customizeCheckBoxMenuItem(_pulse_found_nodes_cbmi, getOptions().isPulseFoundNodes());
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
        _tools_menu.add(_delete_selected_nodes_item = new JMenuItem("Delete Selected Tips"));
        _delete_selected_nodes_item.setToolTipText("To delete all selected external nodes");
        customizeJMenuItem(_delete_selected_nodes_item);
        _tools_menu.add(_delete_not_selected_nodes_item = new JMenuItem("Retain Selected Tips"));
        _delete_not_selected_nodes_item.setToolTipText("To delete all not selected external nodes");
        customizeJMenuItem(_delete_not_selected_nodes_item);
        _tools_menu.addSeparator();
        _tools_menu.add(_select_representative_tips_jmi = new JMenuItem("Select Representative Tips…"));
        customizeJMenuItem(_select_representative_tips_jmi);
        _select_representative_tips_jmi.setToolTipText(
                "Reduce redundancy: group tips that are close in evolutionary distance and keep one representative per group, either within a distance cutoff or to reach a target number. Highlights the representatives and can extract them into a new tree; the current tree is left unchanged.");
        _tools_menu.addSeparator();
        // Colorizing
        _tools_menu.add(_color_rank_jmi = new JMenuItem("Colorize Subtrees via Taxonomic Rank"));
        customizeJMenuItem(_color_rank_jmi);
        _color_rank_jmi.setToolTipText("for example, at \"Class\" level, colorize mammal specific subtree red");
        _tools_menu.add(_clade_bands_jmi = new JMenuItem("Annotate Clades by Rank…"));
        customizeJMenuItem(_clade_bands_jmi);
        _clade_bands_jmi.setToolTipText("mark clades at a chosen rank with shaded boxes or right-edge bars + labels");
        _tools_menu.add(_annotation_columns_jmi = new JMenuItem("Annotation Columns…"));
        customizeJMenuItem(_annotation_columns_jmi);
        _annotation_columns_jmi.setToolTipText(
                "show node annotation fields as tip-aligned columns (color strip, heat map, bar, or text) to the right of the tree");
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
                "Collapse internal branches whose confidence is below a threshold into polytomies (multifurcations). Can be undone (Edit ▸ Undo).");
        collapse_menu.add(_collapse_below_branch_length = new JMenuItem("Collapse Short Branches…"));
        customizeJMenuItem(_collapse_below_branch_length);
        _collapse_below_branch_length.setToolTipText(
                "Collapse internal branches shorter than a threshold into polytomies (multifurcations). Can be undone (Edit ▸ Undo).");
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

    /** Fresh Options seeded from the configuration, then overlaid with the user's persisted display toggles
     *  (so the view the user last chose is restored). Paired with the saveFrom(...) in {@link #exit()}. */
    private Options optionsWithSavedPreferences() {
        final Options options = Options.createInstance(_configuration);
        new GuiPreferences().applyTo(options);
        return options;
    }

    /** {@link #exit()} saves the display toggles on the normal window-close / File&gt;Exit paths. On macOS, Cmd-Q
     *  and the app-menu Quit can bypass those, so route the app-quit request through {@link #close()} (same
     *  unsaved-data confirmation; exit() then saves the settings). macOS-only and best-effort. Unlike a JVM
     *  shutdown hook this fires ONLY on a real user quit -- never during (headless or standalone) test runs, and
     *  it's a per-JVM singleton, so there is no concurrent-write race and no risk of overwriting the developer's
     *  real ~/.archaeopteryx from a test. */
    // Set true ONLY by Archaeopteryx.main -- i.e. when this JVM IS the standalone Archaeopteryx application, not
    // an embedder (msa_compactor, forester_applications) that opens a tree view via Archaeopteryx.createApplication
    // and keeps its own JVM working. Gates the System.exit(0) in exit(); also means a test (which never sets it)
    // can never force-kill the suite JVM even if it reached exit().
    private static boolean _launched_as_standalone_application = false;

    /** Marks this JVM as running the standalone Archaeopteryx application; called only from
     *  {@link Archaeopteryx#main(String[])}. Lets {@link #exit()} force-terminate the JVM on quit. */
    static void setLaunchedAsStandaloneApplication() {
        _launched_as_standalone_application = true;
    }

    private void registerMacOsQuitHandler() {
        try {
            if (Desktop.isDesktopSupported()
                    && Desktop.getDesktop().isSupported(Desktop.Action.APP_QUIT_HANDLER)) {
                Desktop.getDesktop().setQuitHandler((event, response) -> {
                    response.cancelQuit(); // we drive the quit ourselves via close() (which may prompt / dispose)
                    close();
                });
            }
        } catch (final Throwable t) {
            // best-effort: exit()'s save still covers the normal window-close / File>Exit paths
        }
    }

    void exit() {
        new GuiPreferences().saveFrom(getOptions()); // persist the display toggles for the next session
        new DirectoryPreferences().saveFrom(_current_dirs); // and each dialog's last-used directory
        removeAllTextFrames();
        _mainpanel.terminate();
        _contentpane.removeAll();
        setVisible(false);
        dispose();
        if (_launched_as_standalone_application) {
            // Standalone app only: force JVM termination. Disposing the frame is not enough -- a heavyweight
            // rollover Popup (TreePanel/PopupFactory) leaves a cached, still-displayable native window, and the
            // Fetch/Infer/taxonomy worker threads are NOT stopped by terminate(); both keep AWT's toolkit thread
            // (and thus the JVM) alive, so relying on AWT auto-shutdown hangs. NOT done when embedded (a tree view
            // opened via Archaeopteryx.createApplication, e.g. msa_compactor): there a host owns the JVM and
            // closing a window must not kill it or its sibling windows. Runs after the prefs save above.
            System.exit(0);
        }
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
