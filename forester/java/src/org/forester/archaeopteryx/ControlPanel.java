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
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.InputEvent;
import javax.swing.AbstractAction;
import javax.swing.KeyStroke;
import javax.swing.undo.UndoManager;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.UIManager;
import javax.swing.ButtonGroup;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JSlider;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollBar;
import javax.swing.JTextField;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.util.TypomaticJButton;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

final class ControlPanel extends JPanel implements ActionListener {

    enum NodeClickAction {
        ADD_NEW_NODE,
        COLLAPSE,
        COLOR_SUBTREE,
        COPY_SUBTREE,
        CUT_SUBTREE,
        DELETE_NODE_OR_SUBTREE,
        EDIT_NODE_DATA,
        OPEN_SEQ_WEB,
        OPEN_TAX_WEB,
        PASTE_SUBTREE,
        REROOT,
        SELECT_NODES,
        SHOW_DATA,
        SUBTREE,
        SWAP,
        NODE_STYLE,
        UNCOLLAPSE_ALL,
        ORDER_SUBTREE;
    }

    final static Font jcb_bold_font = new Font(Configuration
            .getDefaultFontFamilyName(), Font.BOLD, Configuration.getGuiFontSize());
    final static Font jcb_font = new Font(Configuration
            .getDefaultFontFamilyName(), Font.PLAIN, Configuration.getGuiFontSize());
    final static Font js_font = new Font(Configuration
            .getDefaultFontFamilyName(), Font.PLAIN, Configuration.getGuiFontSize());
    // Two sub-tree navigation buttons: "R" jumps all the way back to the complete tree,
    // "R1" moves up by a single level to the immediate super-tree.
    private static final String RETURN_TO_WHOLE_TREE_TEXT  = "R";
    private static final String RETURN_UP_ONE_LEVEL_TEXT   = "R1";
    // Radial (circular/unrooted) re-labels of the X-/X+ zoom buttons (rotate) and the "W" fit button (label flip).
    // Unicode circular arrows; swap here for a plain-text fallback if a platform's button font lacks the glyphs.
    static final String ROTATE_CCW_LABEL                   = "↺"; // anticlockwise open circle arrow
    static final String ROTATE_CW_LABEL                    = "↻"; // clockwise open circle arrow
    // Plain-text fallback used when the actual button font cannot display the circular-arrow glyph (never on
    // macOS/Windows/mainstream Linux fonts, but keeps a minimal-font platform from showing a missing-glyph box).
    static final String ROTATE_CCW_FALLBACK                = "CCW";
    static final String ROTATE_CW_FALLBACK                 = "CW";
    static final String LABEL_DIRECTION_BUTTON_LABEL       = "L";
    // The zoom cross (Y+, X-, F, E, X+, Y-) holds commonly-used functions, so give those buttons a
    // taller click target than the other small control-panel buttons.
    private static final int    ZOOM_BUTTON_HEIGHT        = 24;
    // The full-width Y+/Y- buttons can be a little shorter than the cross-row buttons.
    private static final int    ZOOM_Y_BUTTON_HEIGHT      = 18;
    // Display checkboxes are packed tightly together (no extra gap between consecutive ones).
    private static final int    CHECKBOX_GAP              = 0;
    private static final String SEARCH_TIP_TEXT = "Enter text to search for. Use ',' for logical OR and '+' for logical AND (not used in this manner for regular expression searches).";
    private static final long serialVersionUID = -8463483932821545633L;
    private NodeClickAction _action_when_node_clicked;
    private int _add_new_node_item;
    private JComboBox<String> _click_to_combobox;
    private JLabel _click_to_label;
    private List<String> _click_to_names;
    private int _collapse_cb_item;
    private int _uncollapse_all_cb_item;
    private int _order_subtree_cb_item;
    private JComboBox<String> _color_by_property_cb;
    private JComboBox<String> _size_by_property_cb;
    // "Ancestral pie:" -- shown adaptively (label + combo hidden together) only when the tree carries BEAST
    // discrete/geographic traits (a `beast:<trait>_set_prob` property); see populateAncestralPieBox().
    private JComboBox<String> _ancestral_pie_property_cb;
    private JLabel            _ancestral_pie_label;
    // "Branch lengths:" (Time / Divergence) -- shown adaptively only for an Auspice/Nextstrain tree that carries BOTH
    // a date and a nextstrain:div property; a reversible display mode. See populateBranchLengthsControl().
    private JComboBox<String> _branch_lengths_cb;
    private JLabel            _branch_lengths_label;
    private static final String COLOR_BY_PROPERTY_NONE = "None";
    private boolean _color_branches;
    private JCheckBox _use_visual_styles_cb;
    private int _color_subtree_cb_item;
    // The settings from the conf file
    private final Configuration _configuration;
    private int _copy_subtree_item;
    private int _cut_subtree_item;
    private JButton _decr_domain_structure_evalue_thr;
    private int _delete_node_or_subtree_item;
    private JRadioButton _display_as_unaligned_phylogram_rb;
    private JRadioButton _display_as_aligned_phylogram_rb;
    private JRadioButton _display_as_cladogram_rb;
    private ButtonGroup _display_as_buttongroup;
    private JRadioButton _light_mode_rb;
    private JRadioButton _dark_mode_rb;
    // Tree checkboxes
    private JCheckBox _display_internal_data;
    private JCheckBox _display_external_data;
    private JLabel _domain_display_label;
    private JTextField _domain_structure_evalue_thr_tf;
    private List<Options.PHYLOGENY_DISPLAY_TYPE> _tree_display_types;
    private JCheckBox _dynamically_hide_data;
    private int _edit_node_data_item;
    private JButton _incr_domain_structure_evalue_thr;
    private final MainPanel _mainpanel;
    private JCheckBox _node_desc_popup_cb;
    private int _open_seq_web_item;
    private int _open_tax_web_item;
    private int _node_style_item;
    private JButton _order;
    private int _paste_subtree_item;
    private int _reroot_cb_item;
    private JButton _return_to_whole_tree;
    private JButton _return_to_super_tree;
    // Search
    private JLabel _search_found_label_0;
    private JLabel _search_found_label_1;
    private JButton _search_reset_button_0;
    private JButton _search_reset_button_1;
    private JTextField _search_tf_0;
    private JTextField _search_tf_1;
    private JButton _search_prev_button; // next/previous step-through of the search/selection hits
    private JButton _search_next_button;
    private JLabel _search_nav_label;    // "k / N"
    private JPanel _search_nav_panel;    // whole row; hidden while there are no hits
    private int _select_nodes_item;
    // Indices for the click-to options in the combo box
    private int _show_data_item;
    private JCheckBox _show_domain_architectures;
    private JCheckBox _write_branch_length_values;
    private JCheckBox _show_events;
    private JCheckBox _show_gene_names;
    private JCheckBox _show_node_names;
    private JCheckBox _shorten_labels_cb;
    private JCheckBox _show_properties_cb;
    private JCheckBox _show_seq_names;
    private JCheckBox _show_seq_symbols;
    private JCheckBox _show_sequence_acc;
    private JCheckBox _show_taxo_code;
    private JCheckBox _show_taxo_rank;
    private JCheckBox _show_taxo_common_names;
    private JCheckBox _show_taxo_scientific_names;
    private JButton _show_whole;
    private JButton _expand_y;
    private JButton _fit_width;
    private Map<String, Color> _species_colors;
    private int _subtree_cb_item;
    private int _swap_cb_item;
    private JButton _uncollapse_all;
    private JCheckBox _width_branches;
    private JCheckBox _write_confidence;
    private JButton _zoom_in_domain_structure;
    private JButton _zoom_in_x;
    private JButton _zoom_in_y;
    private JLabel _zoom_label;
    private JButton _zoom_out_domain_structure;
    private JButton _zoom_out_x;
    private JButton _zoom_out_y;
    // tip-label font size (replaces the old "Font Size" menu); the slider is the user-chosen size
    private JLabel  _font_size_label;
    private JSlider _font_size_slider;
    private boolean _font_slider_is_being_set; // guard so programmatic slider updates don't re-apply
    // "Display Data" checkboxes are shown only when the current tree actually carries the
    // corresponding data (so the panel is not a wall of toggles that silently do nothing). This
    // maps each option constant to its wrapper panel (so a whole row can be collapsed/shown) and
    // caches the last presence scan, keyed by tree identity so pure navigation reuses it.
    private final Map<DisplayOption, JPanel> _checkbox_panels = new EnumMap<DisplayOption, JPanel>(DisplayOption.class);
    private Set<DisplayOption>         _data_presence;
    private Phylogeny                  _data_presence_for;
    // The data-dependent checkboxes (Node Name is intentionally NOT here: it is always shown). Their
    // visibility tracks scanForDataPresence over the whole loaded tree.
    private final static DisplayOption[] DATA_GATED_CHECKBOXES   = {
            DisplayOption.SHOW_TAX_CODE, DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES,
            DisplayOption.SHOW_TAXONOMY_COMMON_NAMES, DisplayOption.SHOW_TAX_RANK, DisplayOption.SHOW_SEQ_NAMES,
            DisplayOption.SHOW_GENE_NAMES, DisplayOption.SHOW_SEQ_SYMBOLS, DisplayOption.SHOW_SEQUENCE_ACC,
            DisplayOption.WRITE_CONFIDENCE_VALUES, DisplayOption.WRITE_BRANCH_LENGTH_VALUES,
            DisplayOption.WRITE_EVENTS, DisplayOption.SHOW_DOMAIN_ARCHITECTURES,
            DisplayOption.SHOW_PROPERTIES, DisplayOption.USE_STYLE,
            DisplayOption.WIDTH_BRANCHES, DisplayOption.SHORTEN_LABELS };

    ControlPanel(final MainPanel ap, final Configuration configuration) {
        init();
        _mainpanel = ap;
        _configuration = configuration;
        if (_configuration.isApplyCustomGuiColors()) {
            setBackground(getConfiguration().getGuiBackgroundColor());
            setBorder(BorderFactory.createRaisedBevelBorder());
        } else {
            // modern look-and-feels: a little breathing room instead of the legacy bevel
            setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));
        }
        setLayout(new GridBagLayout());
        setupControls();
    }

    // ---- Vertical layout of the control-panel rows ------------------------------------------
    // Rows are stacked in a single column. Each row added via add(...) gets a top gap; the gap
    // for the NEXT row can be overridden with nextRowGap(...) to separate visual groups. Unlike
    // the old GridLayout, a hidden row (e.g. the domain-display controls when the tree has no
    // domains) takes up no vertical space.
    private static final int ROW_GAP     = 2;   // default gap between rows within a group
    private static final int SECTION_GAP = 10;  // gap above a new visual group
    private static final int TIGHT_GAP   = 4;   // a deliberately small (but non-zero) gap
    private int              _next_row_top_gap = ROW_GAP + 2; // top margin above the very first row

    /** One-shot: sets the vertical gap above the next row added to the control panel. */
    private void nextRowGap(final int px) {
        _next_row_top_gap = px;
    }

    // Give every plain add(component) call a single full-width column cell with the pending top
    // gap, so the existing setup code keeps working while we control inter-group spacing.
    @Override
    protected void addImpl(final Component comp, final Object constraints, final int index) {
        if ((constraints == null) && (getLayout() instanceof GridBagLayout)) {
            final GridBagConstraints gbc = new GridBagConstraints();
            gbc.gridx = 0;
            gbc.gridy = GridBagConstraints.RELATIVE;
            gbc.weightx = 1.0;
            gbc.fill = GridBagConstraints.HORIZONTAL;
            gbc.anchor = GridBagConstraints.PAGE_START;
            gbc.insets = new Insets(_next_row_top_gap, 0, 0, 0);
            _next_row_top_gap = ROW_GAP;
            super.addImpl(comp, gbc, index);
        } else {
            super.addImpl(comp, constraints, index);
        }
    }

    /** A vertical filler added last so the rows stay top-aligned if the panel is taller than its content. */
    private void addControlPanelGlue() {
        final GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = GridBagConstraints.RELATIVE;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.VERTICAL;
        add(new JLabel(""), gbc);
    }

    /**
     * Handle an action.
     */
    @Override
    public void actionPerformed(final ActionEvent e) {
        try {
            final TreePanel tp = getMainPanel().getCurrentTreePanel();
            if (tp == null) {
                return;
            }
            if (e.getSource() == _color_by_property_cb) {
                final Object sel = _color_by_property_cb.getSelectedItem();
                tp.setColorByPropertyRef(
                        ((sel == null) || COLOR_BY_PROPERTY_NONE.equals(sel)) ? null : sel.toString());
                tp.repaint();
            } else if (e.getSource() == _size_by_property_cb) {
                final Object sel = _size_by_property_cb.getSelectedItem();
                tp.setSizeByPropertyRef(
                        ((sel == null) || COLOR_BY_PROPERTY_NONE.equals(sel)) ? null : sel.toString());
                tp.repaint();
            } else if (e.getSource() == _ancestral_pie_property_cb) {
                final Object sel = _ancestral_pie_property_cb.getSelectedItem();
                tp.setAncestralPieTrait(
                        ((sel == null) || COLOR_BY_PROPERTY_NONE.equals(sel)) ? null : sel.toString());
                tp.repaint();
            } else if (e.getSource() == _branch_lengths_cb) {
                final Object sel = _branch_lengths_cb.getSelectedItem();
                // setNextstrainBranchMode is self-contained: it rewrites the branch lengths, swaps the distance unit +
                // time axis, and re-fits to the viewport (a reversible display mode -- no setEdited / no undo).
                tp.setNextstrainBranchMode(
                        TreePanel.NEXTSTRAIN_BRANCH_MODE.DIVERGENCE.label().equals(sel)
                                ? TreePanel.NEXTSTRAIN_BRANCH_MODE.DIVERGENCE
                                : TreePanel.NEXTSTRAIN_BRANCH_MODE.TIME);
            } else if (e.getSource() == _click_to_combobox) {
                setClickToAction(_click_to_combobox.getSelectedIndex());
                getCurrentTreePanel().repaint();
            } else if (e.getSource() == _show_domain_architectures) {
                reRunSearches();
                // When the user switches domains ON, re-fit the (now wider) tree horizontally so the
                // domains actually become visible -- otherwise it looks like nothing happened. The
                // horizontal-only fit keeps the user's vertical zoom; turning domains OFF just repaints.
                if (_show_domain_architectures.isSelected()) {
                    final MainFrame mf = getMainFrame();
                    if (mf != null) {
                        // radial layout: turning domains on auto-enables "Radial Labels" so they actually show
                        mf.enableRadialLabelsIfDomainsInRadialLayout();
                    }
                    fitWidth();
                } else {
                    displayedPhylogenyMightHaveChanged(true);
                }
            } else if ((tp != null) && (tp.getPhylogeny() != null)) {
                // A user P/A/C click both sets the current tab's display type AND records it as the persisted global
                // DEFAULT (so it survives a restart and opens future branch-length trees the same way). This runs
                // only for a real click -- setTreeDisplayType's setSelected(true) does not fire actionPerformed --
                // so the internal callers (load auto-detect / tabChanged / reset) can't clobber the saved preference.
                if (e.getSource() == getDisplayAsUnalignedPhylogramRb()) {
                    setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM);
                    getOptions().setPhylogenyDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM);
                    showWhole();
                }
                if (e.getSource() == getDisplayAsAlignedPhylogramRb()) {
                    setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM);
                    getOptions().setPhylogenyDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM);
                    showWhole();
                }
                if (e.getSource() == getDisplayAsCladogramRb()) {
                    setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM);
                    getOptions().setPhylogenyDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM);
                    showWhole();
                }
                // Zoom buttons
                // Zoom buttons are labeled by SCREEN direction (X = horizontal, Y = vertical); the screen-oriented
                // helpers flip depth<->breadth automatically in a vertical orientation.
                else if (e.getSource() == _zoom_in_x) {
                    // In a radial layout X+ becomes "rotate clockwise" (X/Y zoom both drive the one radial diameter,
                    // so the second axis is free); otherwise it is the horizontal zoom-in.
                    if (isRadialLayout()) {
                        tp.rotateRadial(true);
                    } else {
                        zoomInScreenX(AptxConstants.BUTTON_ZOOM_IN_FACTOR,
                                AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR);
                        displayedPhylogenyMightHaveChanged(false);
                    }
                } else if (e.getSource() == _zoom_in_y) {
                    zoomInScreenY(AptxConstants.BUTTON_ZOOM_IN_FACTOR, AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR);
                    displayedPhylogenyMightHaveChanged(false);
                } else if (e.getSource() == _zoom_out_x) {
                    // In a radial layout X- becomes "rotate counter-clockwise" (see _zoom_in_x above).
                    if (isRadialLayout()) {
                        tp.rotateRadial(false);
                    } else {
                        zoomOutScreenX(AptxConstants.BUTTON_ZOOM_OUT_FACTOR,
                                AptxConstants.BUTTON_ZOOM_OUT_X_CORRECTION_FACTOR);
                        displayedPhylogenyMightHaveChanged(false);
                    }
                } else if (e.getSource() == _zoom_out_y) {
                    zoomOutScreenY(AptxConstants.BUTTON_ZOOM_OUT_FACTOR,
                            AptxConstants.BUTTON_ZOOM_OUT_X_CORRECTION_FACTOR);
                    displayedPhylogenyMightHaveChanged(false);
                } else if (e.getSource() == _show_whole) {
                    displayedPhylogenyMightHaveChanged(true);
                    showWhole();
                } else if (e.getSource() == _expand_y) {
                    expandYToFitLabels();
                } else if (e.getSource() == _fit_width) {
                    // "W" is fit-width in a horizontal rectangular tree, fit-height ("H") in a vertical one, and the
                    // node-label-direction flip ("L") in a radial layout (where fit-width duplicates "F").
                    if (isRadialLayout()) {
                        toggleNodeLabelDirection();
                    } else if (isVerticalOrientation()) {
                        fitHeight();
                    } else {
                        fitWidth();
                    }
                } else if (e.getSource() == _return_to_whole_tree) {
                    returnedToWholeTreePressed();
                } else if (e.getSource() == _return_to_super_tree) {
                    returnedToSuperTreePressed();
                } else if (e.getSource() == _order) {
                    orderPressed(tp);
                } else if (e.getSource() == _uncollapse_all) {
                    uncollapseAll(tp);
                    displayedPhylogenyMightHaveChanged(false);
                } else if (e.getSource() == _zoom_in_domain_structure) {
                    _mainpanel.getCurrentTreePanel().zoomInDomainStructure();
                    displayedPhylogenyMightHaveChanged(true);
                } else if (e.getSource() == _zoom_out_domain_structure) {
                    _mainpanel.getCurrentTreePanel().zoomOutDomainStructure();
                    displayedPhylogenyMightHaveChanged(true);
                } else if (e.getSource() == _decr_domain_structure_evalue_thr) {
                    _mainpanel.getCurrentTreePanel().decreaseDomainStructureEvalueThresholdExp();
                    reRunSearches();
                    displayedPhylogenyMightHaveChanged(true);
                } else if (e.getSource() == _incr_domain_structure_evalue_thr) {
                    _mainpanel.getCurrentTreePanel().increaseDomainStructureEvalueThresholdExp();
                    reRunSearches();
                    displayedPhylogenyMightHaveChanged(true);
                } else if (e.getSource() == _search_tf_0) {
                    search0();
                } else if (e.getSource() == _search_tf_1) {
                    search1();
                } else if ((_dynamically_hide_data != null) && (e.getSource() == _dynamically_hide_data)
                        && !_dynamically_hide_data.isSelected()) {
                    setDynamicHidingIsOn(false);
                    displayedPhylogenyMightHaveChanged(true);
                } else {
                    displayedPhylogenyMightHaveChanged(true);
                }
            }
            tp.requestFocus();
            tp.requestFocusInWindow();
            tp.requestFocus();
        } catch (final Exception ex) {
            AptxUtil.unexpectedException(ex);
        } catch (final Error err) {
            AptxUtil.unexpectedError(err);
        }
    }

    void orderPressed(final TreePanel tp) {
        DESCENDANT_SORT_PRIORITY pri = DESCENDANT_SORT_PRIORITY.NODE_NAME;
        if (isShowTaxonomyScientificNames() || isShowTaxonomyCode()) {
            pri = DESCENDANT_SORT_PRIORITY.TAXONOMY;
        } else if (isShowSeqNames() || isShowSeqSymbols() || isShowGeneNames()) {
            pri = DESCENDANT_SORT_PRIORITY.SEQUENCE;
        }
        PhylogenyMethods.orderAppearanceX(tp.getPhylogeny().getRoot(), true, pri);
        tp.setNodeInPreorderToNull();
        tp.getPhylogeny().externalNodesHaveChanged();
        tp.getPhylogeny().clearHashIdToNodeMap();
        tp.getPhylogeny().recalculateNumberOfExternalDescendants(true);
        tp.resetNodeIdToDistToLeafMap();
        tp.setEdited(true);
        displayedPhylogenyMightHaveChanged(true);
    }

    void returnedToSuperTreePressed() {
        // "R1": climb exactly one branch toward the root (not one navigation-stack frame, which
        // can span many levels when the user descended by clicking a leaf).
        getCurrentTreePanel().superTreeOneLevel();
    }

    void returnedToWholeTreePressed() {
        final TreePanel tp = getCurrentTreePanel();
        boolean changed = false;
        while (tp.isCurrentTreeIsSubtree()) {
            tp.superTree();
            changed = true;
        }
        if (changed) {
            showWhole();
        }
    }

    public JRadioButton getDisplayAsCladogramRb() {
        return _display_as_cladogram_rb;
    }

    public JRadioButton getDisplayAsAlignedPhylogramRb() {
        return _display_as_aligned_phylogram_rb;
    }

    public JRadioButton getDisplayAsUnalignedPhylogramRb() {
        return _display_as_unaligned_phylogram_rb;
    }

    public JCheckBox getDynamicallyHideData() {
        return _dynamically_hide_data;
    }

    public JCheckBox getNodeDescPopupCb() {
        return _node_desc_popup_cb;
    }

    public JCheckBox getShowEventsCb() {
        return _show_events;
    }

    public JCheckBox getUseVisualStylesCb() {
        return _use_visual_styles_cb;
    }

    public JCheckBox getWriteConfidenceCb() {
        return _write_confidence;
    }

    public boolean isShowProperties() {
        return ((_show_properties_cb != null) && _show_properties_cb.isSelected());
    }

    /** Adds {@code option} to the click-to dropdown (at combo position {@code cb_index}) and returns {@code cb_index}. */
    private int addClickToOption(final int cb_index, final ClickToOption option) {
        _click_to_combobox.addItem(option.title());
        _click_to_names.add(option.title());
        if (_configuration.isApplyCustomGuiColors()) {
            _click_to_combobox.setBackground(getConfiguration().getGuiButtonBackgroundColor());
            _click_to_combobox.setForeground(getConfiguration().getGuiButtonTextColor());
        }
        return cb_index;
    }

    private List<Options.PHYLOGENY_DISPLAY_TYPE> getTreeDisplayTypes() {
        return _tree_display_types;
    }

    final private MainFrame getMainFrame() {
        return getMainPanel().getMainFrame();
    }

    private void init() {
        _tree_display_types = new ArrayList<Options.PHYLOGENY_DISPLAY_TYPE>();
        setSpeciesColors(new HashMap<String, Color>());
    }

    private Options.PHYLOGENY_DISPLAY_TYPE getTreeDisplayType(final int index) {
        return getTreeDisplayTypes().get(index);
    }

    private void search0(final MainPanel main_panel, final Phylogeny tree, final String query_str) {
        getSearchFoundCountsLabel0().setVisible(true);
        getSearchResetButton0().setEnabled(true);
        getSearchResetButton0().setVisible(true);
        final Set<Long> nodes = runSearch(tree, (SearchField) _search_field_0.getSelectedItem(),
                (SearchMode) _search_mode_0.getSelectedItem(), query_str, _search_range_tf_0.getText());
        if ((nodes != null) && !nodes.isEmpty()) {
            // Hand the finished set to the panel in one call: setFoundNodes0 is the chokepoint that resets the
            // step-through position (only when the hit set actually changed) and refreshes the "k / N" navigator.
            main_panel.getCurrentTreePanel().setFoundNodes0(nodes);
            setSearchFoundCountsOnLabel0(nodes.size());
        } else {
            setSearchFoundCountsOnLabel0(0);
            searchReset0();
        }
    }

    private void search1(final MainPanel main_panel, final Phylogeny tree, final String query_str) {
        getSearchFoundCountsLabel1().setVisible(true);
        getSearchResetButton1().setEnabled(true);
        getSearchResetButton1().setVisible(true);
        final Set<Long> nodes = runSearch(tree, (SearchField) _search_field_1.getSelectedItem(),
                (SearchMode) _search_mode_1.getSelectedItem(), query_str, _search_range_tf_1.getText());
        if ((nodes != null) && !nodes.isEmpty()) {
            main_panel.getCurrentTreePanel().setFoundNodes1(nodes); // see search0: single chokepoint for the reset+nav
            setSearchFoundCountsOnLabel1(nodes.size());
        } else {
            setSearchFoundCountsOnLabel1(0);
            searchReset1();
        }
    }

    /**
     * Runs one search box against {@code tree} using its chosen {@link SearchField} and {@link SearchMode} plus
     * the two shared modifiers (Match Case, Inverse). For a NUMERIC field, {@code value} is the operand (the low
     * bound for {@link SearchMode#RANGE}) and {@code range_high} the range upper bound. For a TEXT field the
     * legacy multi-term combinators apply to {@code value}: ',' = OR and '+' = AND. Inverse is applied ONCE at the
     * end as the complement over data-bearing nodes, but only for an actual query (a separator-only text query or
     * an incomplete range produces no terms and resets, rather than selecting the whole tree).
     */
    private Set<Long> runSearch(final Phylogeny tree, final SearchField field, final SearchMode mode, String value,
                                final String range_high) {
        // guard: field and mode must be the same kind; a transient mismatch (should be reconciled first) must
        // reset the search, never build a mismatched SearchSpec (which throws -> an unexpected-exception dialog).
        if ((field == null) || (mode == null) || (field.isNumeric() != mode.isNumeric())) {
            return null;
        }
        final boolean inverse = (_search_inverse_cb != null) && _search_inverse_cb.isSelected();
        value = (value == null) ? "" : value.replaceAll("\\s+", " ").trim();
        final Set<Long> nodes = new HashSet<>();
        boolean any_term = false; // whether an actual query was run (so Inverse selects the complement, not all)
        if (field.isNumeric()) {
            final String hi = (range_high == null) ? "" : range_high.trim();
            final Double lo_num = SearchMatcher.parseFiniteDouble(value);
            final Double hi_num = SearchMatcher.parseFiniteDouble(hi);
            // an actual numeric query needs a PARSEABLE operand (both bounds for a range); an unparseable operand
            // is a no-op that must reset -- with Inverse on it must NOT select the whole tree (the Phase-2 rule).
            any_term = (mode == SearchMode.RANGE) ? ((lo_num != null) && (hi_num != null)) : (lo_num != null);
            if (any_term) {
                nodes.addAll(SearchMatcher.search(
                        new SearchSpec(field, mode, value, hi.isEmpty() ? null : hi, false, false), tree));
            }
        } else {
            final boolean case_sensitive = (_search_case_sensitive_cb != null)
                    && _search_case_sensitive_cb.isSelected();
            // ',' OR and '+' AND only apply to plain text matching, not to regular expressions (which can contain
            // those characters), consistent with the legacy behaviour.
            final boolean splittable = (mode != SearchMode.REGEX);
            final String[] or_terms = (splittable && (value.indexOf(',') >= 0)) ? value.split(",+")
                    : new String[] { value };
            for (String or_term : or_terms) {
                or_term = or_term.trim();
                if (ForesterUtil.isEmpty(or_term)) {
                    continue;
                }
                any_term = true;
                if (splittable && (or_term.indexOf('+') > 0)) {
                    nodes.addAll(searchLogicalAnd(tree, field, mode, or_term.split("\\++"), case_sensitive));
                } else {
                    nodes.addAll(SearchMatcher
                            .search(new SearchSpec(field, mode, or_term, null, case_sensitive, false), tree));
                }
            }
        }
        return (inverse && any_term) ? complement(tree, field, nodes) : nodes;
    }

    /** The nodes matching EVERY '+'-separated term (positive) -- the intersection of the per-term matches. */
    private static Set<Long> searchLogicalAnd(final Phylogeny tree, final SearchField field, final SearchMode mode,
                                              final String[] and_terms, final boolean case_sensitive) {
        Set<Long> acc = null;
        for (String term : and_terms) {
            term = term.trim();
            if (ForesterUtil.isEmpty(term)) {
                continue;
            }
            final Set<Long> matches = new HashSet<>(
                    SearchMatcher.search(new SearchSpec(field, mode, term, null, case_sensitive, false), tree));
            if (acc == null) {
                acc = matches;
            } else {
                acc.retainAll(matches);
            }
        }
        return (acc == null) ? new HashSet<>() : acc;
    }

    /** The complement of {@code matched} over the nodes that CARRY {@code field} (the "Inverse" modifier): the
     *  inverse of a field-scoped search is over the nodes the field applies to, not every data-bearing node -- so
     *  e.g. "inverse of scientific-name = X" no longer selects tips that have no scientific name at all. Applies to
     *  text and numeric fields alike (see {@link SearchField#carries}). */
    private static Set<Long> complement(final Phylogeny tree, final SearchField field, final Set<Long> matched) {
        final Set<Long> out = new HashSet<>();
        for (final PhylogenyNode n : PhylogenyMethods.obtainAllNodesAsList(tree)) {
            if (!matched.contains(n.getId()) && field.carries(n)) {
                out.add(n.getId());
            }
        }
        return out;
    }

    /** The distinct existing values the value-autocomplete popup offers for one search box: the selected field's
     *  values across the current tree -- but only for a SPECIFIC text field in a non-regex mode. Empty (so no popup)
     *  for "Any text field" (a mixed value set, not a meaningful pick list), a numeric field (continuous), regex mode
     *  (a pattern, not a value), or no tree. Recomputed on each popup open (the popup does not cache), so it can never
     *  go stale against a tree/data edit. */
    private List<String> autocompleteValues(final boolean box_a) {
        final JComboBox<SearchField> fc = box_a ? _search_field_0 : _search_field_1;
        final JComboBox<SearchMode> mc = box_a ? _search_mode_0 : _search_mode_1;
        if ((fc == null) || (mc == null)) {
            return new ArrayList<String>();
        }
        final SearchField f = (SearchField) fc.getSelectedItem();
        final SearchMode m = (SearchMode) mc.getSelectedItem();
        if ((f == null) || (m == null) || f.isNumeric() || (f.kind() == SearchField.Kind.ANY_TEXT)
                || (m == SearchMode.REGEX)) {
            return new ArrayList<String>();
        }
        final Phylogeny phy = (getMainPanel() == null) ? null : getMainPanel().getCurrentPhylogeny();
        return SearchField.distinctValues(phy, f);
    }

    /** Whether a key does NOT modify the query text -- arrows / caret navigation, and the autocomplete popup's
     *  Down/Up/Escape/Tab. The value fields re-run the search on keyReleased; skipping these keeps arrowing through the
     *  suggestion popup (or moving the caret) from re-running the whole search on every keypress. ENTER is deliberately
     *  NOT here -- it still runs the search (a plain Enter, or after accepting a suggestion). */
    static boolean isNonEditingKey(final int key_code) {
        switch (key_code) {
            case KeyEvent.VK_UP:
            case KeyEvent.VK_DOWN:
            case KeyEvent.VK_LEFT:
            case KeyEvent.VK_RIGHT:
            case KeyEvent.VK_HOME:
            case KeyEvent.VK_END:
            case KeyEvent.VK_PAGE_UP:
            case KeyEvent.VK_PAGE_DOWN:
            case KeyEvent.VK_ESCAPE:
            case KeyEvent.VK_TAB:
            case KeyEvent.VK_SHIFT:
            case KeyEvent.VK_CONTROL:
            case KeyEvent.VK_ALT:
            case KeyEvent.VK_META:
                return true;
            default:
                return false;
        }
    }

    /** Closes a box's value-autocomplete popup (and drops its cached session values), e.g. when its field or mode
     *  changed so the value set / eligibility differs; the next focus/keystroke recomputes from the live tree. */
    private void dismissAutocomplete(final boolean box_a) {
        final SearchValueAutocomplete ac = box_a ? _search_autocomplete_0 : _search_autocomplete_1;
        if (ac != null) {
            ac.endSession();
        }
    }

    /** Test hook: the value-autocomplete candidate list for a box (see {@link #autocompleteValues(boolean)}). */
    List<String> autocompleteValuesForTest(final boolean box_a) {
        return autocompleteValues(box_a);
    }

    private void setTreeDisplayType(final int index, final Options.PHYLOGENY_DISPLAY_TYPE t) {
        getTreeDisplayTypes().set(index, t);
    }

    private void setupClickToOptions() {
        // Every click-to action is always offered; the tree-editing subset (Edit Node Data, Delete, Add,
        // Cut/Copy/Paste) is gated per-item by Options.isEditable(), so on a non-editable tree those rows
        // are simply skipped and the surrounding always-shown actions close up around them.
        final ClickToOption default_option = _configuration.getDefaultDisplayClicktoOption();
        int selected_index = 0;
        int cb_index = 0;
        final boolean editable = getOptions().isEditable();
        _show_data_item = addClickToOption(cb_index++, ClickToOption.DISPLAY_NODE_DATA);
        selected_index = clickToSelectedIndex(selected_index, _show_data_item, default_option, ClickToOption.DISPLAY_NODE_DATA);
        _subtree_cb_item = addClickToOption(cb_index++, ClickToOption.SUBTREE);
        selected_index = clickToSelectedIndex(selected_index, _subtree_cb_item, default_option, ClickToOption.SUBTREE);
        _select_nodes_item = addClickToOption(cb_index++, ClickToOption.SELECT_NODES);
        selected_index = clickToSelectedIndex(selected_index, _select_nodes_item, default_option, ClickToOption.SELECT_NODES);
        _collapse_cb_item = addClickToOption(cb_index++, ClickToOption.COLLAPSE_UNCOLLAPSE);
        selected_index = clickToSelectedIndex(selected_index, _collapse_cb_item, default_option, ClickToOption.COLLAPSE_UNCOLLAPSE);
        _uncollapse_all_cb_item = addClickToOption(cb_index++, ClickToOption.UNCOLLAPSE_ALL);
        selected_index = clickToSelectedIndex(selected_index, _uncollapse_all_cb_item, default_option, ClickToOption.UNCOLLAPSE_ALL);
        _swap_cb_item = addClickToOption(cb_index++, ClickToOption.SWAP);
        selected_index = clickToSelectedIndex(selected_index, _swap_cb_item, default_option, ClickToOption.SWAP);
        _order_subtree_cb_item = addClickToOption(cb_index++, ClickToOption.ORDER_SUBTREE);
        selected_index = clickToSelectedIndex(selected_index, _order_subtree_cb_item, default_option, ClickToOption.ORDER_SUBTREE);
        _open_tax_web_item = addClickToOption(cb_index++, ClickToOption.OPEN_TAX_WEB);
        selected_index = clickToSelectedIndex(selected_index, _open_tax_web_item, default_option, ClickToOption.OPEN_TAX_WEB);
        _open_seq_web_item = addClickToOption(cb_index++, ClickToOption.OPEN_SEQ_WEB);
        selected_index = clickToSelectedIndex(selected_index, _open_seq_web_item, default_option, ClickToOption.OPEN_SEQ_WEB);
        _node_style_item = addClickToOption(cb_index++, ClickToOption.NODE_STYLE);
        selected_index = clickToSelectedIndex(selected_index, _node_style_item, default_option, ClickToOption.NODE_STYLE);
        _color_subtree_cb_item = addClickToOption(cb_index++, ClickToOption.COLOR_SUBTREE);
        selected_index = clickToSelectedIndex(selected_index, _color_subtree_cb_item, default_option, ClickToOption.COLOR_SUBTREE);
        if (editable) {
            _edit_node_data_item = addClickToOption(cb_index++, ClickToOption.EDIT_NODE_DATA);
            selected_index = clickToSelectedIndex(selected_index, _edit_node_data_item, default_option, ClickToOption.EDIT_NODE_DATA);
        }
        _reroot_cb_item = addClickToOption(cb_index++, ClickToOption.REROOT);
        selected_index = clickToSelectedIndex(selected_index, _reroot_cb_item, default_option, ClickToOption.REROOT);
        if (editable) {
            _delete_node_or_subtree_item = addClickToOption(cb_index++, ClickToOption.DELETE_SUBTREE_OR_NODE);
            selected_index = clickToSelectedIndex(selected_index, _delete_node_or_subtree_item, default_option, ClickToOption.DELETE_SUBTREE_OR_NODE);
            _add_new_node_item = addClickToOption(cb_index++, ClickToOption.ADD_NEW_NODE);
            selected_index = clickToSelectedIndex(selected_index, _add_new_node_item, default_option, ClickToOption.ADD_NEW_NODE);
            _cut_subtree_item = addClickToOption(cb_index++, ClickToOption.CUT_SUBTREE);
            selected_index = clickToSelectedIndex(selected_index, _cut_subtree_item, default_option, ClickToOption.CUT_SUBTREE);
            _copy_subtree_item = addClickToOption(cb_index++, ClickToOption.COPY_SUBTREE);
            selected_index = clickToSelectedIndex(selected_index, _copy_subtree_item, default_option, ClickToOption.COPY_SUBTREE);
            _paste_subtree_item = addClickToOption(cb_index++, ClickToOption.PASTE_SUBTREE);
            selected_index = clickToSelectedIndex(selected_index, _paste_subtree_item, default_option, ClickToOption.PASTE_SUBTREE);
        }
        // Set default selection and its action
        _click_to_combobox.setSelectedIndex(selected_index);
        setClickToAction(selected_index);
    }

    /** Returns {@code candidate_index} when {@code option} is the default click-to action, else keeps {@code current}. */
    private static int clickToSelectedIndex(final int current, final int candidate_index, final ClickToOption default_option,
                                            final ClickToOption option) {
        return (default_option == option) ? candidate_index : current;
    }

    private void setupDisplayCheckboxes() {
        addDisplayCheckbox(DisplayOption.DYNAMICALLY_HIDE_DATA);
        addDisplayCheckbox(DisplayOption.NODE_DATA_POPUP);
        addDisplayCheckbox(DisplayOption.DISPLAY_INTERNAL_DATA);
        addDisplayCheckbox(DisplayOption.DISPLAY_EXTERNAL_DATA);
        addDisplayCheckbox(DisplayOption.USE_STYLE);
        addDisplayCheckbox(DisplayOption.WIDTH_BRANCHES);
        final JLabel label = new JLabel("Display Data:");
        label.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        add(label);
        addDisplayCheckbox(DisplayOption.SHOW_NODE_NAMES);
        addDisplayCheckbox(DisplayOption.SHORTEN_LABELS);
        addDisplayCheckbox(DisplayOption.SHOW_TAX_CODE);
        addDisplayCheckbox(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES);
        addDisplayCheckbox(DisplayOption.SHOW_TAXONOMY_COMMON_NAMES);
        addDisplayCheckbox(DisplayOption.SHOW_TAX_RANK);
        addDisplayCheckbox(DisplayOption.SHOW_SEQ_NAMES);
        addDisplayCheckbox(DisplayOption.SHOW_GENE_NAMES);
        addDisplayCheckbox(DisplayOption.SHOW_SEQ_SYMBOLS);
        addDisplayCheckbox(DisplayOption.SHOW_SEQUENCE_ACC);
        addDisplayCheckbox(DisplayOption.WRITE_CONFIDENCE_VALUES);
        addDisplayCheckbox(DisplayOption.WRITE_BRANCH_LENGTH_VALUES);
        addDisplayCheckbox(DisplayOption.SHOW_DOMAIN_ARCHITECTURES);
        addDisplayCheckbox(DisplayOption.WRITE_EVENTS);
        addDisplayCheckbox(DisplayOption.SHOW_PROPERTIES);
        // Data-dependent checkboxes start hidden; the first scan of the loaded tree reveals only the
        // ones whose data is actually present (Node Name is intentionally always shown).
        for (final DisplayOption which : DATA_GATED_CHECKBOXES) {
            final JPanel p = _checkbox_panels.get(which);
            if (p != null) {
                p.setVisible(false);
            }
        }
    }

    // Build a "Display Data" checkbox and set its initial checked state. Whether it is actually
    // visible is then governed by updateDataCheckboxVisibility (data presence in the current tree).
    private void addDisplayCheckbox(final DisplayOption which) {
        addCheckbox(which, which.title());
        setCheckbox(which, which.isCheckedByDefault());
    }

    private void setVisibilityOfDomainStrucureControls() {
        if (_zoom_in_domain_structure != null) {
            final MainFrame mf = getMainFrame();
            if (mf != null) {
                if (isShowDomainArchitectures()) {
                    _domain_display_label.setVisible(true);
                    _zoom_in_domain_structure.setVisible(true);
                    _zoom_out_domain_structure.setVisible(true);
                    _decr_domain_structure_evalue_thr.setVisible(true);
                    _incr_domain_structure_evalue_thr.setVisible(true);
                    _domain_structure_evalue_thr_tf.setVisible(true);
                } else {
                    _domain_display_label.setVisible(false);
                    _zoom_in_domain_structure.setVisible(false);
                    _zoom_out_domain_structure.setVisible(false);
                    _decr_domain_structure_evalue_thr.setVisible(false);
                    _incr_domain_structure_evalue_thr.setVisible(false);
                    _domain_structure_evalue_thr_tf.setVisible(false);
                }
            }
        }
    }

    void activateButtonsToReturnToSuperTree() {
        _return_to_whole_tree.setForeground(getConfiguration().getGuiCheckboxAndButtonActiveColor());
        _return_to_whole_tree.setEnabled(true);
        _return_to_super_tree.setForeground(getConfiguration().getGuiCheckboxAndButtonActiveColor());
        _return_to_super_tree.setEnabled(true);
    }

    void activateButtonToUncollapseAll() {
        _uncollapse_all.setForeground(getConfiguration().getGuiCheckboxAndButtonActiveColor());
        _uncollapse_all.setEnabled(true);
    }

    /**
     * Add zoom and quick edit buttons. (Last modified 8/9/04)
     */
    void addButtons() {
        final JPanel x_panel = new JPanel(new GridLayout(1, 1, 0, 0));
        final JPanel y_panel = new JPanel(new GridLayout(1, 5, 0, 0));
        final JPanel z_panel = new JPanel(new GridLayout(1, 1, 0, 0));
        final JPanel o_panel = new JPanel(new GridLayout(1, 4, 0, 0));
        if (getConfiguration().isApplyCustomGuiColors()) {
            x_panel.setBackground(getBackground());
            y_panel.setBackground(getBackground());
            z_panel.setBackground(getBackground());
            o_panel.setBackground(getBackground());
        }
        nextRowGap(SECTION_GAP);
        add(_zoom_label = new JLabel("Zoom:"));
        customizeLabel(_zoom_label, getConfiguration());
        add(x_panel);
        add(y_panel);
        add(z_panel);
        if (getConfiguration().isUseNativeUI()) {
            _zoom_in_x = new TypomaticJButton("+");
            _zoom_out_x = new TypomaticJButton("-");
        } else {
            _zoom_in_x = new TypomaticJButton("X+");
            _zoom_out_x = new TypomaticJButton("X-");
        }
        _zoom_in_y = new TypomaticJButton("Y+");
        _zoom_out_y = new TypomaticJButton("Y-");
        _show_whole = new JButton("F");
        _show_whole.setToolTipText("fit and center tree display [Alt+C, Home, or Esc]");
        _expand_y = new JButton("E");
        _expand_y.setToolTipText("expand the tree in vertical direction so labels do not overlap at the current font size [Alt+E]");
        _fit_width = new JButton("W");
        _fit_width.setToolTipText("fit the tree to the window width, keeping the current vertical zoom [Alt+W]");
        _zoom_in_x.setToolTipText("zoom in horizontally [Alt+Right or Shift+Alt+mousewheel]");
        _zoom_in_y.setToolTipText("zoom in vertically [Alt+Up or Shift+mousewheel]");
        _zoom_out_x.setToolTipText("zoom out horizontally [Alt+Left or Shift+Alt+mousewheel]");
        _zoom_out_y.setToolTipText("zoom out vertically [Alt+Down or Shift+mousewheel]");
        if (getConfiguration().isUseNativeUI() && ForesterUtil.isMac()) {
            _zoom_out_x.setPreferredSize(new Dimension(55, ZOOM_BUTTON_HEIGHT));
            _zoom_in_x.setPreferredSize(new Dimension(55, ZOOM_BUTTON_HEIGHT));
        } else {
            _zoom_out_x.setPreferredSize(new Dimension(10, ZOOM_BUTTON_HEIGHT));
            _zoom_in_x.setPreferredSize(new Dimension(10, ZOOM_BUTTON_HEIGHT));
        }
        _zoom_out_y.setPreferredSize(new Dimension(10, ZOOM_Y_BUTTON_HEIGHT));
        _zoom_in_y.setPreferredSize(new Dimension(10, ZOOM_Y_BUTTON_HEIGHT));
        _show_whole.setPreferredSize(new Dimension(10, ZOOM_BUTTON_HEIGHT));
        _expand_y.setPreferredSize(new Dimension(10, ZOOM_BUTTON_HEIGHT));
        _fit_width.setPreferredSize(new Dimension(10, ZOOM_BUTTON_HEIGHT));
        // The middle zoom row now holds five buttons (X- F E W X+); trim the default button
        // padding so the two-character X-/X+ labels still fit instead of being clipped to "...".
        final Insets tight = new Insets(2, 1, 2, 1);
        _zoom_out_x.setMargin(tight);
        _show_whole.setMargin(tight);
        _expand_y.setMargin(tight);
        _fit_width.setMargin(tight);
        _zoom_in_x.setMargin(tight);
        _return_to_whole_tree = new JButton(RETURN_TO_WHOLE_TREE_TEXT);
        _return_to_whole_tree.setToolTipText("return all the way to the complete tree (if in a sub-tree) [Alt+Shift+R]");
        _return_to_whole_tree.setEnabled(false);
        _return_to_super_tree = new JButton(RETURN_UP_ONE_LEVEL_TEXT);
        _return_to_super_tree.setToolTipText("move up by one level towards the complete tree (if in a sub-tree) [Alt+R]");
        _return_to_super_tree.setEnabled(false);
        _order = new JButton("O");
        _order.setToolTipText("order all [Alt+O]");
        _uncollapse_all = new JButton("U");
        _uncollapse_all.setToolTipText("uncollapse all [Alt+U]");
        // Four buttons share the bottom row (O R R1 U); trim the default padding so the
        // two-character "R1" label is not clipped to "..." under FlatLaf.
        _order.setMargin(tight);
        _return_to_whole_tree.setMargin(tight);
        _return_to_super_tree.setMargin(tight);
        _uncollapse_all.setMargin(tight);
        addJButton(_zoom_in_y, x_panel);
        addJButton(_zoom_out_x, y_panel);
        addJButton(_show_whole, y_panel);
        addJButton(_expand_y, y_panel);
        addJButton(_fit_width, y_panel);
        addJButton(_zoom_in_x, y_panel);
        addJButton(_zoom_out_y, z_panel);
        // tip-label font-size slider (replaces the retired "Font Size" menu)
        nextRowGap(SECTION_GAP);
        add(_font_size_label = new JLabel(fontSizeLabelText(AptxConstants.DEFAULT_TREE_FONT_SIZE)));
        customizeLabel(_font_size_label, getConfiguration());
        _font_size_slider = new JSlider(TreeFontSet.MIN_FONT_SIZE, TreeFontSet.MAX_FONT_SIZE,
                AptxConstants.DEFAULT_TREE_FONT_SIZE);
        _font_size_slider.setToolTipText("tip-label font size [PageUp/PageDown, Shift +/-, or Ctrl+mousewheel]");
        // a small preferred width so the slider does not widen the control panel; the full-width GridBag row
        // stretches it to the panel width (set by the comboboxes), like the zoom buttons do
        _font_size_slider.setPreferredSize(new Dimension(10, _font_size_slider.getPreferredSize().height));
        if (!getConfiguration().isUseNativeUI()) {
            _font_size_slider.setBackground(getBackground());
        }
        _font_size_slider.addChangeListener(e -> fontSizeSliderChanged());
        add(_font_size_slider);
        nextRowGap(SECTION_GAP);
        add(o_panel);
        addJButton(_order, o_panel);
        addJButton(_return_to_whole_tree, o_panel);
        addJButton(_return_to_super_tree, o_panel);
        addJButton(_uncollapse_all, o_panel);
        nextRowGap(SECTION_GAP);
        setUpControlsForDomainStrucures();
        setVisibilityOfDomainStrucureControls();
    }

    /** Applies the slider value as the user font size (on release / discrete change); keeps the label live. */
    private void fontSizeSliderChanged() {
        if (_font_slider_is_being_set) {
            return;
        }
        final int size = _font_size_slider.getValue();
        _font_size_label.setText(fontSizeLabelText(size)); // live feedback while dragging
        if (_font_size_slider.getValueIsAdjusting()) {
            return; // defer the (re)layout until the drag settles
        }
        getMainPanel().getTreeFontSet().setUserFontSize(size);
        displayedPhylogenyMightHaveChanged(true);
    }

    /**
     * Syncs the font-size slider + label to the font size currently <i>displayed</i> ({@code getLargeFont}),
     * so it tracks the overlap auto-fit too -- e.g. after F/E/W or a horizontal zoom shrink/grow the labels.
     */
    void updateFontSizeSlider() {
        if (_font_size_slider == null) {
            return;
        }
        final TreeFontSet tfs = getMainPanel().getTreeFontSet();
        if (tfs == null) {
            return;
        }
        final int size = Math.max(TreeFontSet.MIN_FONT_SIZE,
                Math.min(TreeFontSet.MAX_FONT_SIZE, tfs.getLargeFont().getSize()));
        if (_font_size_slider.getValue() != size) {
            _font_slider_is_being_set = true; // don't let setValue() re-trigger an apply
            _font_size_slider.setValue(size);
            _font_slider_is_being_set = false;
        }
        _font_size_label.setText(fontSizeLabelText(size));
    }

    private static String fontSizeLabelText(final int size) {
        return "Font size: " + size;
    }

    /** For tests: the current slider value (-1 if the slider is not built). */
    int getFontSizeSliderValue() {
        return (_font_size_slider == null) ? -1 : _font_size_slider.getValue();
    }

    /** For tests: the property refs offered in the "Color by" dropdown (excluding the "None" entry). */
    java.util.List<String> colorByPropertyRefs() {
        final java.util.List<String> refs = new java.util.ArrayList<>();
        if (_color_by_property_cb != null) {
            for (int i = 0; i < _color_by_property_cb.getItemCount(); ++i) {
                final String item = _color_by_property_cb.getItemAt(i);
                if (!COLOR_BY_PROPERTY_NONE.equals(item)) {
                    refs.add(item);
                }
            }
        }
        return refs;
    }

    /** For tests: the "Size by" dropdown's selected item as a string ("None" when off / not built). */
    String getSizeByPropertySelection() {
        if (_size_by_property_cb == null) {
            return COLOR_BY_PROPERTY_NONE;
        }
        final Object sel = _size_by_property_cb.getSelectedItem();
        return (sel == null) ? COLOR_BY_PROPERTY_NONE : sel.toString();
    }

    /** For tests: move the slider (fires the change listener, so it applies the size like a user drag-release). */
    void setFontSizeSliderValue(final int value) {
        if (_font_size_slider != null) {
            _font_size_slider.setValue(value);
        }
    }

    void addCheckbox(final DisplayOption which, final String title) {
        nextRowGap(CHECKBOX_GAP); // pack the display checkboxes tightly together
        final JPanel ch_panel = new JPanel(new BorderLayout(0, 0));
        switch (which) {
            case DISPLAY_INTERNAL_DATA:
                _display_internal_data = new JCheckBox(title);
                _display_internal_data.setToolTipText("To allow or disallow display of internal labels");
                addJCheckBox(_display_internal_data, ch_panel);
                add(ch_panel);
                break;
            case DISPLAY_EXTERNAL_DATA:
                _display_external_data = new JCheckBox(title);
                _display_external_data.setToolTipText("To allow or disallow display of external labels");
                addJCheckBox(_display_external_data, ch_panel);
                add(ch_panel);
                break;
            case SHOW_NODE_NAMES:
                _show_node_names = new JCheckBox(title);
                addJCheckBox(_show_node_names, ch_panel);
                add(ch_panel);
                break;
            case SHORTEN_LABELS:
                _shorten_labels_cb = new JCheckBox(title);
                _shorten_labels_cb.setToolTipText(
                        "Display over-long external labels (e.g. full UniProt/NCBI descriptions) shortened with an ellipsis; the underlying names are not changed");
                addJCheckBox(_shorten_labels_cb, ch_panel);
                add(ch_panel);
                break;
            case SHOW_TAXONOMY_SCIENTIFIC_NAMES:
                _show_taxo_scientific_names = new JCheckBox(title);
                addJCheckBox(_show_taxo_scientific_names, ch_panel);
                add(ch_panel);
                break;
            case SHOW_TAXONOMY_COMMON_NAMES:
                _show_taxo_common_names = new JCheckBox(title);
                addJCheckBox(_show_taxo_common_names, ch_panel);
                add(ch_panel);
                break;
            case SHOW_TAX_CODE:
                _show_taxo_code = new JCheckBox(title);
                addJCheckBox(_show_taxo_code, ch_panel);
                add(ch_panel);
                break;
            case SHOW_TAX_RANK:
                _show_taxo_rank = new JCheckBox(title);
                addJCheckBox(_show_taxo_rank, ch_panel);
                add(ch_panel);
                break;
            case WRITE_CONFIDENCE_VALUES:
                _write_confidence = new JCheckBox(title);
                addJCheckBox(getWriteConfidenceCb(), ch_panel);
                add(ch_panel);
                break;
            case WRITE_EVENTS:
                _show_events = new JCheckBox(title);
                addJCheckBox(getShowEventsCb(), ch_panel);
                add(ch_panel);
                break;
            case USE_STYLE:
                _use_visual_styles_cb = new JCheckBox(title);
                getUseVisualStylesCb()
                        .setToolTipText("To use visual styles (node colors, fonts) and branch colors, if present");
                addJCheckBox(getUseVisualStylesCb(), ch_panel);
                add(ch_panel);
                break;
            case WIDTH_BRANCHES:
                _width_branches = new JCheckBox(title);
                _width_branches.setToolTipText("To use branch width values, if present");
                addJCheckBox(_width_branches, ch_panel);
                add(ch_panel);
                break;
            case WRITE_BRANCH_LENGTH_VALUES:
                _write_branch_length_values = new JCheckBox(title);
                addJCheckBox(_write_branch_length_values, ch_panel);
                add(ch_panel);
                break;
            case SHOW_DOMAIN_ARCHITECTURES:
                _show_domain_architectures = new JCheckBox(title);
                addJCheckBox(_show_domain_architectures, ch_panel);
                add(ch_panel);
                break;
            case SHOW_SEQ_NAMES:
                _show_seq_names = new JCheckBox(title);
                addJCheckBox(_show_seq_names, ch_panel);
                add(ch_panel);
                break;
            case SHOW_GENE_NAMES:
                _show_gene_names = new JCheckBox(title);
                addJCheckBox(_show_gene_names, ch_panel);
                add(ch_panel);
                break;
            case SHOW_SEQ_SYMBOLS:
                _show_seq_symbols = new JCheckBox(title);
                addJCheckBox(_show_seq_symbols, ch_panel);
                add(ch_panel);
                break;
            case SHOW_SEQUENCE_ACC:
                _show_sequence_acc = new JCheckBox(title);
                addJCheckBox(_show_sequence_acc, ch_panel);
                add(ch_panel);
                break;
            case DYNAMICALLY_HIDE_DATA:
                _dynamically_hide_data = new JCheckBox(title);
                getDynamicallyHideData().setToolTipText("To hide labels depending on expected visibility");
                addJCheckBox(getDynamicallyHideData(), ch_panel);
                add(ch_panel);
                break;
            case NODE_DATA_POPUP:
                _node_desc_popup_cb = new JCheckBox(title);
                getNodeDescPopupCb().setToolTipText("To enable mouse rollover display of basic node data");
                addJCheckBox(getNodeDescPopupCb(), ch_panel);
                add(ch_panel);
                break;
            case SHOW_PROPERTIES:
                _show_properties_cb = new JCheckBox(title);
                addJCheckBox(_show_properties_cb, ch_panel);
                add(ch_panel);
                break;
            default:
                throw new RuntimeException("unknown checkbox: " + which);
        }
        _checkbox_panels.put(which, ch_panel);
    }// addCheckbox

    void addJButton(final JButton jb, final JPanel p) {
        jb.setFocusPainted(false);
        jb.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            jb.setBorder(BorderFactory.createLineBorder(getConfiguration().getGuiButtonBorderColor()));
            jb.setBackground(getConfiguration().getGuiButtonBackgroundColor());
            jb.setForeground(getConfiguration().getGuiButtonTextColor());
        }
        p.add(jb);
        jb.addActionListener(this);
    }

    void addJCheckBox(final JCheckBox jcb, final JPanel p) {
        jcb.setFocusPainted(false);
        jcb.setFont(ControlPanel.jcb_font);
        jcb.setMargin(new Insets(0, 0, 0, 0)); // trim vertical padding so the checkboxes pack tightly
        if (_configuration.isApplyCustomGuiColors()) {
            jcb.setBackground(getConfiguration().getGuiBackgroundColor());
            jcb.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        p.add(jcb, "Center");
        jcb.addActionListener(this);
    }

    private final void setupJRadioButton(final JRadioButton rb) {
        rb.setFocusPainted(false);
        rb.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            rb.setBackground(getConfiguration().getGuiBackgroundColor());
            rb.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        rb.addActionListener(this);
    }

    void addJTextField(final JTextField tf, final JPanel p) {
        if (_configuration.isApplyCustomGuiColors()) {
            tf.setForeground(getConfiguration().getGuiBackgroundColor());
            tf.setFont(ControlPanel.jcb_font);
        }
        p.add(tf);
        tf.addActionListener(this);
    }

    void deactivateButtonsToReturnToSuperTree() {
        _return_to_whole_tree.setForeground(getConfiguration().getGuiButtonTextColor());
        _return_to_whole_tree.setEnabled(false);
        _return_to_super_tree.setForeground(getConfiguration().getGuiButtonTextColor());
        _return_to_super_tree.setEnabled(false);
    }

    void deactivateButtonToUncollapseAll() {
        _uncollapse_all.setForeground(getConfiguration().getGuiButtonTextColor());
        _uncollapse_all.setEnabled(false);
    }

    /**
     * Shows only the "Display Data" checkboxes whose data is actually present somewhere in the whole
     * loaded tree, and collapses the rest. The presence scan is cached by tree identity, so pure
     * navigation (same tree) reuses it; pass {@code force_rescan == true} after the tree's data may
     * have changed (load, tab switch, node/annotation edit, add/remove node). Node Name is never
     * gated here -- it is always shown.
     */
    void updateDataCheckboxVisibility(final boolean force_rescan) {
        if ((_mainpanel == null) || (_mainpanel.getCurrentPhylogeny() == null)
                || _mainpanel.getCurrentPhylogeny().isEmpty()) {
            return;
        }
        final Phylogeny phy = _mainpanel.getCurrentPhylogeny();
        final Set<DisplayOption> prev_presence = _data_presence;
        if (force_rescan || (phy != _data_presence_for) || (_data_presence == null)) {
            _data_presence = AptxUtil.scanForDataPresence(phy);
            _data_presence_for = phy;
        }
        boolean changed = false;
        for (final DisplayOption which : DATA_GATED_CHECKBOXES) {
            final JPanel p = _checkbox_panels.get(which);
            if (p != null) {
                final boolean vis = _data_presence.contains(which);
                if (p.isVisible() != vis) {
                    p.setVisible(vis);
                    changed = true;
                }
            }
        }
        if (changed) {
            revalidate();
            repaint();
        }
        // Tailor the search field selectors to this tree, but force the (O(n)) rebuild ONLY when the tree's set of
        // present data types actually CHANGED -- i.e. an in-place data edit added/removed a field (Extract Data from
        // Labels, a node-data edit; Import forces its own rebuild). A structural/display event (collapse, prune, swap,
        // reorder, a display-checkbox toggle) keeps the same data-presence set, so it takes the identity-guarded
        // rebuildSearchFields(false) path and pays nothing; a NEW tree still rebuilds via that method's phy-identity
        // guard. This keeps the frequent structural/display events off the extra full-tree availableFields scan.
        final boolean presence_changed = !java.util.Objects.equals(prev_presence, _data_presence);
        rebuildSearchFields(presence_changed);
    }

    // For tests: whether the "Display Data" checkbox row for the given option constant is showing.
    boolean isDisplayDataCheckboxVisible(final DisplayOption which) {
        final JPanel p = _checkbox_panels.get(which);
        return (p != null) && p.isVisible();
    }

    /** Re-seeds the always-visible control-panel controls that hold their OWN state (theme radios; the two search
     *  modifier checkboxes; the per-box field/mode selectors) from the current Configuration/Options, for Reset to
     *  Defaults. Without this they stay stale after a reset -- and worse, the search checkboxes write their state
     *  back to Options on the next click, which would silently clobber the reset. Uses setSelected (fires no
     *  ActionListener) for the checkboxes; the combos are re-seeded under the _search_controls_adjusting guard so
     *  their ActionListeners don't launch a search mid-reset. */
    void resyncFromOptions() {
        if ( ( _light_mode_rb != null ) && ( _dark_mode_rb != null ) ) {
            final boolean dark = getConfiguration().getUi() == Configuration.UI.FLAT_DARK;
            _light_mode_rb.setSelected( !dark );
            _dark_mode_rb.setSelected( dark );
        }
        final Options o = getOptions();
        if ( o != null ) {
            if ( _search_case_sensitive_cb != null ) {
                _search_case_sensitive_cb.setSelected( o.isSearchCaseSensitive() );
            }
            if ( _search_inverse_cb != null ) {
                _search_inverse_cb.setSelected( o.isInverseSearchResult() );
            }
        }
        _search_controls_adjusting = true;
        try {
            // clear the remember-last state so the reset really returns to the defaults (Any text / contains / =)
            _last_field_label_0 = null;
            _last_field_label_1 = null;
            _last_string_mode_0 = SearchMode.CONTAINS;
            _last_string_mode_1 = SearchMode.CONTAINS;
            _last_numeric_mode_0 = SearchMode.EQ;
            _last_numeric_mode_1 = SearchMode.EQ;
            if ( ( _search_field_0 != null ) && ( _search_field_0.getItemCount() > 0 ) ) {
                _search_field_0.setSelectedIndex( 0 ); // "Any text field"
            }
            if ( ( _search_field_1 != null ) && ( _search_field_1.getItemCount() > 0 ) ) {
                _search_field_1.setSelectedIndex( 0 );
            }
            // Any text is a string field, so make sure the mode combos hold the string set, then default them to
            // "contains" and hide the range fields.
            reconcileModeCombo( true );
            reconcileModeCombo( false );
            if ( ( _search_mode_0 != null ) && ( _search_mode_0.getItemCount() > 0 ) ) {
                _search_mode_0.setSelectedIndex( 0 );
            }
            if ( ( _search_mode_1 != null ) && ( _search_mode_1.getItemCount() > 0 ) ) {
                _search_mode_1.setSelectedIndex( 0 );
            }
            updateRangeFieldVisibility( true );
            updateRangeFieldVisibility( false );
            if ( _search_range_tf_0 != null ) {
                _search_range_tf_0.setText( "" );
            }
            if ( _search_range_tf_1 != null ) {
                _search_range_tf_1.setText( "" );
            }
            if ( _search_combine_combo != null ) {
                _search_combine_combo.setSelectedIndex( 0 ); // back to independent
            }
            _search_was_combined = false;
        }
        finally {
            _search_controls_adjusting = false;
        }
        updateCombineControlVisibility();
    }

    /** Resets the "Color by" property dropdown to None (for Reset to Defaults). The per-tab coloring state is
     *  cleared separately on each TreePanel; this just brings the shared dropdown back in line for the shown tab. */
    void setColorByPropertySelectionToNone() {
        if ( _color_by_property_cb != null ) {
            _color_by_property_cb.setSelectedItem( COLOR_BY_PROPERTY_NONE );
        }
    }

    void displayedPhylogenyMightHaveChanged(final boolean recalc_longest_ext_node_info) {
        if ((_mainpanel != null)
                && ((_mainpanel.getCurrentPhylogeny() != null) && !_mainpanel.getCurrentPhylogeny().isEmpty())) {
            if (recalc_longest_ext_node_info) {
                _mainpanel.getCurrentTreePanel().initNodeData();
                _mainpanel.getCurrentTreePanel().calculateLongestExtNodeInfo();
            }
            updateDataCheckboxVisibility(recalc_longest_ext_node_info);
            if (getOptions().isShowOverview()) {
                _mainpanel.getCurrentTreePanel().updateOvSizes();
            }
            _mainpanel.getCurrentTreePanel().recalculateMaxDistanceToRoot();
            _mainpanel.getCurrentTreePanel().rebuildPropertyDisplays();
            _mainpanel.getCurrentTreePanel().rebuildAnnotationColumns();
            setVisibilityOfDomainStrucureControls();
            updateDomainStructureEvaluethresholdDisplay();
            getMainPanel().getControlPanel();
            _mainpanel.getCurrentTreePanel().updateButtonToUncollapseAll();
            _mainpanel.getCurrentTreePanel().calculateScaleDistance();
            _mainpanel.getCurrentTreePanel().calcMaxDepth();
            _mainpanel.adjustJScrollPane();
            _mainpanel.getCurrentTreePanel().repaint();
            // _mainpanel.getCurrentTreePanel().setUpUrtFactor();
            updateFontSizeSlider(); // keep the slider in sync with keyboard/wheel nudges and the font dialog
            // Safety net for the "k / N" navigator: a tree edit (prune/delete-subtree) can drop found nodes without
            // touching the found-set objects, so the setFoundNodes0/1 chokepoint would miss the shrunk count.
            updateSearchHitNavigation();
        }
    }

    void endClickToOptions() {
        _click_to_combobox.addActionListener(this);
    }

    /**
     * Indicates what action should be execute when a node is clicked
     *
     * @return the click-on action
     */
    NodeClickAction getActionWhenNodeClicked() {
        return _action_when_node_clicked;
    }

    Configuration getConfiguration() {
        return _configuration;
    }

    TreePanel getCurrentTreePanel() {
        return getMainPanel().getCurrentTreePanel();
    }

    MainPanel getMainPanel() {
        return _mainpanel;
    }

    Options getOptions() {
        return getMainPanel().getOptions();
    }

    JLabel getSearchFoundCountsLabel0() {
        return _search_found_label_0;
    }

    JLabel getSearchFoundCountsLabel1() {
        return _search_found_label_1;
    }

    JButton getSearchResetButton0() {
        return _search_reset_button_0;
    }

    JButton getSearchResetButton1() {
        return _search_reset_button_1;
    }

    JTextField getSearchTextField0() {
        return _search_tf_0;
    }

    JTextField getSearchTextField1() {
        return _search_tf_1;
    }

    List<String> getSingleClickToNames() {
        return _click_to_names;
    }

    Map<String, Color> getSpeciesColors() {
        return _species_colors;
    }

    boolean isDrawPhylogram() {
        final Options.PHYLOGENY_DISPLAY_TYPE t = getTreeDisplayType(getMainPanel().getCurrentTabIndex());
        return ((t == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM)
                || (t == Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM));
    }

    boolean isDynamicallyHideData() {
        return ((getDynamicallyHideData() != null) && getDynamicallyHideData().isSelected());
    }

    boolean isEvents() {
        return ((getShowEventsCb() != null) && getShowEventsCb().isSelected());
    }

    boolean isNodeDescPopup() {
        return ((getNodeDescPopupCb() != null) && getNodeDescPopupCb().isSelected());
    }

    boolean isShowConfidenceValues() {
        return ((getWriteConfidenceCb() != null) && getWriteConfidenceCb().isSelected());
    }

    boolean isShowDomainArchitectures() {
        return ((_show_domain_architectures != null) && _show_domain_architectures.isSelected());
    }

    /** Test hook: force the "Show Domain Architectures" checkbox on/off. */
    void setShowDomainArchitecturesForTest(final boolean selected) {
        if (_show_domain_architectures != null) {
            _show_domain_architectures.setSelected(selected);
        }
    }

    /**
     * Gives a text field its own Undo/Redo (Cmd/Ctrl-Z, Cmd/Ctrl-Shift-Z) bound at WHEN_FOCUSED scope, so
     * those keystrokes edit the field's text instead of falling through to the tree-level Undo menu
     * accelerator -- which would otherwise revert the tree while the user is typing a search query.
     */
    private static void installTextUndo(final JTextField tf) {
        final UndoManager um = new UndoManager();
        tf.getDocument().addUndoableEditListener(um);
        final int shortcut = Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx();
        tf.getInputMap().put(KeyStroke.getKeyStroke(KeyEvent.VK_Z, shortcut), "text-undo");
        tf.getInputMap().put(KeyStroke.getKeyStroke(KeyEvent.VK_Z, shortcut | InputEvent.SHIFT_DOWN_MASK),
                "text-redo");
        tf.getActionMap().put("text-undo", new AbstractAction() {

            @Override
            public void actionPerformed(final ActionEvent e) {
                if (um.canUndo()) {
                    um.undo();
                }
            }
        });
        tf.getActionMap().put("text-redo", new AbstractAction() {

            @Override
            public void actionPerformed(final ActionEvent e) {
                if (um.canRedo()) {
                    um.redo();
                }
            }
        });
    }

    /**
     * If the "Domain Architectures" checkbox is available, switch it on and fit the tree -- now wider
     * because of the domain rows -- to the screen. Called on load for a tree that contains domain
     * architectures so the domains are shown immediately (most users never find the checkbox).
     */
    // On load, open the domain architectures at ~this fraction of the viewport width (the widest architecture), so
    // they read at a useful size out of the box instead of a fixed ~90px. The user's domain zoom (+/-) scales from here.
    private static final double INITIAL_DOMAIN_WIDTH_SCREEN_FRACTION = 0.25;

    void showDomainArchitecturesFitted() {
        if ((_show_domain_architectures != null) && _show_domain_architectures.isVisible()) {
            _show_domain_architectures.setSelected(true); // no-op (and fires nothing) if already on
            // Set THIS tree's initial domain width regardless of whether the checkbox was already selected by a
            // previously-open domain tree: _domain_structure_width is per-TreePanel and a freshly-loaded tree has
            // no user domain-zoom to clobber. (This method runs only from the per-load property scan, never on a
            // manual checkbox toggle, so it can't reset the width out from under a user who is zooming domains.)
            getCurrentTreePanel().fitDomainWidthToScreen(INITIAL_DOMAIN_WIDTH_SCREEN_FRACTION,
                    getMainPanel().getSizeOfViewport().width);
            final MainFrame mf = getMainFrame();
            if (mf != null) {
                // a tree that opens in a radial layout (persisted graphics type): auto-enable "Radial Labels" so its
                // domains show right away (they are suppressed under horizontal labels)
                mf.enableRadialLabelsIfDomainsInRadialLayout();
            }
            fitWidth(); // re-fit the (now wider) tree horizontally so the domains show; keep vertical zoom
        }
    }

    boolean isShowGeneNames() {
        return ((_show_gene_names != null) && _show_gene_names.isSelected());
    }

    boolean isShowInternalData() {
        return ((_display_internal_data == null) || _display_internal_data.isSelected());
    }

    boolean isShowExternalData() {
        return ((_display_external_data == null) || _display_external_data.isSelected());
    }

    boolean isShowNodeNames() {
        return ((_show_node_names != null) && _show_node_names.isSelected());
    }

    boolean isShortenLabels() {
        return ((_shorten_labels_cb != null) && _shorten_labels_cb.isSelected());
    }

    boolean isShowSeqNames() {
        return ((_show_seq_names != null) && _show_seq_names.isSelected());
    }

    boolean isShowSeqSymbols() {
        return ((_show_seq_symbols != null) && _show_seq_symbols.isSelected());
    }

    boolean isShowSequenceAcc() {
        return ((_show_sequence_acc != null) && _show_sequence_acc.isSelected());
    }

    boolean isShowTaxonomyCode() {
        return ((_show_taxo_code != null) && _show_taxo_code.isSelected());
    }

    boolean isShowTaxonomyRank() {
        return ((_show_taxo_rank != null) && _show_taxo_rank.isSelected());
    }

    boolean isShowTaxonomyCommonNames() {
        return ((_show_taxo_common_names != null) && _show_taxo_common_names.isSelected());
    }

    boolean isShowTaxonomyScientificNames() {
        return ((_show_taxo_scientific_names != null) && _show_taxo_scientific_names.isSelected());
    }

    boolean isUseVisualStyles() {
        return (((getUseVisualStylesCb() != null) && getUseVisualStylesCb().isSelected())
                || ((getUseVisualStylesCb() == null) && _color_branches));
    }

    boolean isWidthBranches() {
        return ((_width_branches != null) && _width_branches.isSelected());
    }

    boolean isWriteBranchLengthValues() {
        return ((_write_branch_length_values != null) && _write_branch_length_values.isSelected());
    }

    void phylogenyAdded(final Configuration configuration) {
        if (configuration.isDrawAsPhylogram()) {
            getTreeDisplayTypes().add(Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM);
        } else {
            getTreeDisplayTypes().add(Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM);
        }
    }

    void phylogenyRemoved(final int index) {
        getTreeDisplayTypes().remove(index);
    }

    void search0() {
        runSearchBox(true);
    }

    void search1() {
        runSearchBox(false);
    }

    /**
     * Runs a search box, honouring the "Combine A & B" mode. When it is AND/OR and BOTH boxes carry a query, the two
     * boxes are combined into a SINGLE found set (A &cap; B / A &cup; B) in found-set 0; otherwise the two boxes search
     * INDEPENDENTLY (the default A/B dual-highlight). Leaving combine mode (a box cleared, or the mode set back to
     * independent) re-runs BOTH boxes independently so their separate highlights are restored.
     */
    private void runSearchBox(final boolean box_a) {
        final MainPanel main_panel = getMainPanel();
        final Phylogeny tree = main_panel.getCurrentPhylogeny();
        if ((tree == null) || tree.isEmpty()) {
            return;
        }
        if (combineActive()) {
            runCombinedSearch(tree);
            _search_was_combined = true;
        } else if (_search_was_combined) {
            _search_was_combined = false;
            independentSearch(true);  // restore BOTH boxes' own found sets (the combined result had cleared box B)
            independentSearch(false);
        } else {
            independentSearch(box_a);
        }
        updateCombineControlVisibility();
        updateSearchFieldValidity(box_a);
        repaintTreeAfterSearch(); // light: only the highlight changed (no full layout recalc)
    }

    /** The independent (single-box) search: compute the box's own found set from its field/mode/value -- the default,
     *  dual-highlight behaviour. */
    private void independentSearch(final boolean box_a) {
        final MainPanel main_panel = getMainPanel();
        final Phylogeny tree = main_panel.getCurrentPhylogeny();
        if ((tree == null) || tree.isEmpty()) {
            return;
        }
        String query = (box_a ? getSearchTextField0() : getSearchTextField1()).getText();
        if (query != null) {
            query = query.trim();
        }
        if (!ForesterUtil.isEmpty(query)) {
            if (box_a) {
                search0(main_panel, tree, query);
            } else {
                search1(main_panel, tree, query);
            }
        } else {
            final JLabel lbl = box_a ? getSearchFoundCountsLabel0() : getSearchFoundCountsLabel1();
            final JButton btn = box_a ? getSearchResetButton0() : getSearchResetButton1();
            lbl.setVisible(false);
            btn.setEnabled(false);
            btn.setVisible(false);
            if (box_a) {
                searchReset0();
            } else {
                searchReset1();
            }
        }
    }

    /** The combine selection: 0 = independent (default) / 1 = AND (A &cap; B) / 2 = OR (A &cup; B). */
    private int combineModeIndex() {
        return (_search_combine_combo == null) ? 0 : _search_combine_combo.getSelectedIndex();
    }

    /** Whether BOTH search boxes currently hold a non-blank query. */
    private boolean bothQueriesActive() {
        return (_search_tf_0 != null) && (_search_tf_1 != null)
                && !ForesterUtil.isEmptyTrimmed(_search_tf_0.getText())
                && !ForesterUtil.isEmptyTrimmed(_search_tf_1.getText());
    }

    /** Whether the two boxes should be COMBINED into one found set: a non-independent mode AND both boxes have a query. */
    private boolean combineActive() {
        return (combineModeIndex() > 0) && bothQueriesActive();
    }

    /** A description of the active combine mode for the menu-bar counter ("A AND B" / "A OR B"), or null when the two
     *  boxes are independent -- so the counter labels a combined result correctly instead of misreporting it as "A". */
    String searchCombineDescription() {
        return combineActive() ? ((combineModeIndex() == 1) ? "A AND B" : "A OR B") : null;
    }

    /** Re-runs BOTH search boxes (after a display or shared-option change). In combine mode a single search0() already
     *  recomputes the combined A&cap;B/A&cup;B for both boxes, so search1() is skipped -- otherwise the full combined
     *  search would run twice. */
    private void reRunSearches() {
        search0();
        if (!combineActive()) {
            search1();
        }
    }

    /** Shows the "Combine A & B" control only when both boxes have a query (adaptive -- zero panel space for the
     *  common single-search case, like the step-through navigator). */
    private void updateCombineControlVisibility() {
        if (_search_combine_panel == null) {
            return;
        }
        final boolean show = bothQueriesActive();
        if (_search_combine_panel.isVisible() != show) {
            _search_combine_panel.setVisible(show);
            revalidate();
            repaint();
        }
    }

    /** The raw found set of one box (its own field/mode/value + the ',' OR / '+' AND / Inverse handling), WITHOUT
     *  installing it on the panel -- used to combine the two boxes. */
    private Set<Long> rawFoundSet(final boolean box_a, final Phylogeny tree) {
        final String q = (box_a ? _search_tf_0 : _search_tf_1).getText();
        if (ForesterUtil.isEmptyTrimmed(q)) {
            return new HashSet<>();
        }
        final SearchField f = (SearchField) (box_a ? _search_field_0 : _search_field_1).getSelectedItem();
        final SearchMode m = (SearchMode) (box_a ? _search_mode_0 : _search_mode_1).getSelectedItem();
        final String range = (box_a ? _search_range_tf_0 : _search_range_tf_1).getText();
        final Set<Long> r = runSearch(tree, f, m, q.trim(), range);
        return (r == null) ? new HashSet<>() : r;
    }

    /** Runs the combined A/B search: A &cap; B (AND) or A &cup; B (OR), placing the single result in found-set 0
     *  (found-set 1 cleared) so the highlight, count, step-through, the menu-bar counter and export all reflect it. */
    private void runCombinedSearch(final Phylogeny tree) {
        final TreePanel tp = getMainPanel().getCurrentTreePanel();
        if (tp == null) {
            return;
        }
        final Set<Long> combined = new HashSet<>(rawFoundSet(true, tree));
        if (combineModeIndex() == 1) {
            combined.retainAll(rawFoundSet(false, tree)); // AND
        } else {
            combined.addAll(rawFoundSet(false, tree));     // OR
        }
        tp.setFoundNodes1(null); // the combined result lives in found-set 0 only
        if (!combined.isEmpty()) {
            tp.setFoundNodes0(combined);
            setSearchFoundCountsOnLabel0(combined.size());
        } else {
            searchReset0();
            setSearchFoundCountsOnLabel0(0);
        }
        getSearchFoundCountsLabel0().setVisible(true);
        getSearchResetButton0().setEnabled(true);
        getSearchResetButton0().setVisible(true);
        // box B's own count / reset are not meaningful while the two boxes are combined into one result
        getSearchFoundCountsLabel1().setVisible(false);
        getSearchResetButton1().setEnabled(false);
        getSearchResetButton1().setVisible(false);
    }

    /** The "Combine A & B" control (independent / AND / OR), hidden until both boxes carry a query. */
    void setupSearchCombine() {
        _search_combine_panel = new JPanel(new BorderLayout(4, 0));
        _search_combine_panel.setBackground(getBackground());
        final JLabel lbl = new JLabel("Combine:");
        lbl.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            lbl.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _search_combine_combo = new JComboBox<>(new String[] { "independent", "A AND B", "A OR B" });
        _search_combine_combo.setSelectedIndex(0);
        _search_combine_combo.setToolTipText("How to combine the two search boxes: independent A/B highlights, or one "
                + "result set -- A AND B (in both) / A OR B (in either)");
        _search_combine_combo.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            _search_combine_combo.setBackground(getConfiguration().getGuiButtonBackgroundColor());
            _search_combine_combo.setForeground(getConfiguration().getGuiButtonTextColor());
        }
        _search_combine_combo.setPreferredSize(new Dimension(10, _search_combine_combo.getPreferredSize().height));
        _search_combine_combo.addActionListener(e -> {
            if (!_search_controls_adjusting) {
                search0(); // recompute with the new combine mode (routes through runSearchBox)
            }
        });
        _search_combine_panel.add(lbl, BorderLayout.WEST);
        _search_combine_panel.add(_search_combine_combo, BorderLayout.CENTER);
        _search_combine_panel.setVisible(false); // adaptive: shown only when both boxes have a query
        add(_search_combine_panel);
    }

    /** Test hook: set the combine mode (0 independent / 1 AND / 2 OR); fires the real recompute. */
    void setSearchCombineForTest(final int index) {
        if (_search_combine_combo != null) {
            _search_combine_combo.setSelectedIndex(index);
        }
    }

    /** Test hook: whether the adaptive "Combine A & B" control is currently showing. */
    boolean isSearchCombineControlVisibleForTest() {
        return (_search_combine_panel != null) && _search_combine_panel.isVisible();
    }

    /** Test hook: the combine control's laid-out height (0 when hidden/not yet laid out). */
    int searchCombinePanelHeightForTest() {
        return (_search_combine_panel == null) ? 0 : _search_combine_panel.getHeight();
    }

    void searchReset0() {
        if (getMainPanel().getCurrentTreePanel() != null) {
            getMainPanel().getCurrentTreePanel().setFoundNodes0(null);
        }
    }

    void searchReset1() {
        if (getMainPanel().getCurrentTreePanel() != null) {
            getMainPanel().getCurrentTreePanel().setFoundNodes1(null);
        }
    }

    void setActionWhenNodeClicked(final NodeClickAction action) {
        _action_when_node_clicked = action;
    }

    void setCheckbox(final DisplayOption which, final boolean state) {
        switch (which) {
            case DISPLAY_AS_PHYLOGRAM:
                if (getDisplayAsUnalignedPhylogramRb() != null) {
                    getDisplayAsUnalignedPhylogramRb().setSelected(state);
                    getDisplayAsAlignedPhylogramRb().setSelected(!state);
                    getDisplayAsCladogramRb().setSelected(!state);
                }
                break;
            case DISPLAY_INTERNAL_DATA:
                if (_display_internal_data != null) {
                    _display_internal_data.setSelected(state);
                }
                break;
            case DISPLAY_EXTERNAL_DATA:
                if (_display_external_data != null) {
                    _display_external_data.setSelected(state);
                }
                break;
            case SHOW_NODE_NAMES:
                if (_show_node_names != null) {
                    _show_node_names.setSelected(state);
                }
                break;
            case SHORTEN_LABELS:
                if (_shorten_labels_cb != null) {
                    _shorten_labels_cb.setSelected(state);
                }
                break;
            case SHOW_TAXONOMY_SCIENTIFIC_NAMES:
                if (_show_taxo_scientific_names != null) {
                    _show_taxo_scientific_names.setSelected(state);
                }
                break;
            case SHOW_TAXONOMY_COMMON_NAMES:
                if (_show_taxo_common_names != null) {
                    _show_taxo_common_names.setSelected(state);
                }
                break;
            case SHOW_TAX_CODE:
                if (_show_taxo_code != null) {
                    _show_taxo_code.setSelected(state);
                }
                break;
            case SHOW_TAX_RANK:
                if (_show_taxo_rank != null) {
                    _show_taxo_rank.setSelected(state);
                }
                break;
            case WRITE_CONFIDENCE_VALUES:
                if (getWriteConfidenceCb() != null) {
                    getWriteConfidenceCb().setSelected(state);
                }
                break;
            case WRITE_EVENTS:
                if (getShowEventsCb() != null) {
                    getShowEventsCb().setSelected(state);
                }
                break;
            case USE_STYLE:
                if (getUseVisualStylesCb() != null) {
                    getUseVisualStylesCb().setSelected(state);
                }
                break;
            case WIDTH_BRANCHES:
                if (_width_branches != null) {
                    _width_branches.setSelected(state);
                }
                break;
            case SHOW_DOMAIN_ARCHITECTURES:
                if (_show_domain_architectures != null) {
                    _show_domain_architectures.setSelected(state);
                }
                break;
            case WRITE_BRANCH_LENGTH_VALUES:
                if (_write_branch_length_values != null) {
                    _write_branch_length_values.setSelected(state);
                }
                break;
            case SHOW_SEQ_NAMES:
                if (_show_seq_names != null) {
                    _show_seq_names.setSelected(state);
                }
                break;
            case SHOW_GENE_NAMES:
                if (_show_gene_names != null) {
                    _show_gene_names.setSelected(state);
                }
                break;
            case SHOW_SEQ_SYMBOLS:
                if (_show_seq_symbols != null) {
                    _show_seq_symbols.setSelected(state);
                }
                break;
            case SHOW_PROPERTIES:
                if (_show_properties_cb != null) {
                    _show_properties_cb.setSelected(state);
                }
                break;
            case SHOW_SEQUENCE_ACC:
                if (_show_sequence_acc != null) {
                    _show_sequence_acc.setSelected(state);
                }
                break;
            case DYNAMICALLY_HIDE_DATA:
                if (getDynamicallyHideData() != null) {
                    getDynamicallyHideData().setSelected(state);
                }
                break;
            case NODE_DATA_POPUP:
                if (getNodeDescPopupCb() != null) {
                    getNodeDescPopupCb().setSelected(state);
                }
                break;
            default:
                throw new AssertionError("unknown checkbox: " + which);
        }
    }

    /**
     * Set this checkbox state. Not all checkboxes have been instantiated
     * depending on the config.
     */
    void setCheckbox(final JCheckBox cb, final boolean state) {
        if (cb != null) {
            cb.setSelected(state);
        }
    }

    void setClickToAction(final int action) {
        // Set click-to action
        if (action == _show_data_item) {
            setActionWhenNodeClicked(NodeClickAction.SHOW_DATA);
        } else if (action == _collapse_cb_item) {
            setActionWhenNodeClicked(NodeClickAction.COLLAPSE);
        } else if (action == _reroot_cb_item) {
            setActionWhenNodeClicked(NodeClickAction.REROOT);
        } else if (action == _subtree_cb_item) {
            setActionWhenNodeClicked(NodeClickAction.SUBTREE);
        } else if (action == _swap_cb_item) {
            setActionWhenNodeClicked(NodeClickAction.SWAP);
        } else if (action == _color_subtree_cb_item) {
            setActionWhenNodeClicked(NodeClickAction.COLOR_SUBTREE);
        } else if (action == _open_seq_web_item) {
            setActionWhenNodeClicked(NodeClickAction.OPEN_SEQ_WEB);
        } else if (action == _open_tax_web_item) {
            setActionWhenNodeClicked(NodeClickAction.OPEN_TAX_WEB);
        } else if (action == _cut_subtree_item) {
            setActionWhenNodeClicked(NodeClickAction.CUT_SUBTREE);
        } else if (action == _copy_subtree_item) {
            setActionWhenNodeClicked(NodeClickAction.COPY_SUBTREE);
        } else if (action == _delete_node_or_subtree_item) {
            setActionWhenNodeClicked(NodeClickAction.DELETE_NODE_OR_SUBTREE);
        } else if (action == _paste_subtree_item) {
            setActionWhenNodeClicked(NodeClickAction.PASTE_SUBTREE);
        } else if (action == _add_new_node_item) {
            setActionWhenNodeClicked(NodeClickAction.ADD_NEW_NODE);
        } else if (action == _edit_node_data_item) {
            setActionWhenNodeClicked(NodeClickAction.EDIT_NODE_DATA);
        } else if (action == _select_nodes_item) {
            setActionWhenNodeClicked(NodeClickAction.SELECT_NODES);
        } else if (action == _node_style_item) {
            setActionWhenNodeClicked(NodeClickAction.NODE_STYLE);
        } else if (action == _uncollapse_all_cb_item) {
            setActionWhenNodeClicked(NodeClickAction.UNCOLLAPSE_ALL);
        } else if (action == _order_subtree_cb_item) {
            setActionWhenNodeClicked(NodeClickAction.ORDER_SUBTREE);
        } else {
            throw new RuntimeException("unknown action: " + action);
        }
        // make sure drop down is displaying the correct action
        // in case this was called from outside the class
        _click_to_combobox.setSelectedIndex(action);
    }

    void setColorBranches(final boolean color_branches) {
        _color_branches = color_branches;
    }

    void setTreeDisplayType(final Options.PHYLOGENY_DISPLAY_TYPE t) {
        switch (t) {
            case UNALIGNED_PHYLOGRAM:
                getDisplayAsUnalignedPhylogramRb().setSelected(true);
                break;
            case ALIGNED_PHYLOGRAM:
                getDisplayAsAlignedPhylogramRb().setSelected(true);
                break;
            case CLADOGRAM:
                getDisplayAsCladogramRb().setSelected(true);
                break;
        }
        setTreeDisplayType(getMainPanel().getCurrentTabIndex(), t);
    }

    void setDrawPhylogramEnabled(final boolean b) {
        if (getDisplayAsAlignedPhylogramRb() != null && getDisplayAsUnalignedPhylogramRb() != null
                && getDisplayAsCladogramRb() != null) {
            getDisplayAsAlignedPhylogramRb().setEnabled(b);
            getDisplayAsUnalignedPhylogramRb().setEnabled(b);
            getDisplayAsCladogramRb().setEnabled(b);
        }
    }

    void setDynamicHidingIsOn(final boolean is_on) {
        if (is_on) {
            getDynamicallyHideData().setForeground(getConfiguration().getGuiCheckboxAndButtonActiveColor());
        } else {
            if (_configuration.isApplyCustomGuiColors()) {
                getDynamicallyHideData().setForeground(getConfiguration().getGuiButtonTextColor());
            } else {
                // reset to the look-and-feel default so the label stays visible in dark themes
                getDynamicallyHideData().setForeground(UIManager.getColor("CheckBox.foreground"));
            }
        }
    }

    void setSearchFoundCountsOnLabel0(final int counts) {
        getSearchFoundCountsLabel0().setText("Found: " + counts);
    }

    void setSearchFoundCountsOnLabel1(final int counts) {
        getSearchFoundCountsLabel1().setText("Found: " + counts);
    }

    void setShowEvents(final boolean show_events) {
        if (getShowEventsCb() == null) {
            _show_events = new JCheckBox("");
        }
        getShowEventsCb().setSelected(show_events);
    }

    void setSpeciesColors(final Map<String, Color> species_colors) {
        _species_colors = species_colors;
    }

    /**
     * A "Color by:" dropdown that colors the leaves on the fly by the value of a chosen
     * phyloXML property (e.g. host). It is (re)populated per displayed tree; see
     * {@link #populateColorByPropertyBox()}.
     */
    void setupColorByProperty() {
        final JLabel label = new JLabel("Color by:");
        label.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _color_by_property_cb = new JComboBox<String>();
        _color_by_property_cb.setFont(ControlPanel.js_font);
        _color_by_property_cb.setToolTipText("color leaves by the value of a phyloXML property");
        _color_by_property_cb.addItem(COLOR_BY_PROPERTY_NONE);
        // show a friendly property name (no namespace prefix, underscores as spaces,
        // capitalized) in the dropdown -- the underlying ref is unchanged
        _color_by_property_cb.setRenderer(new DefaultListCellRenderer() {
            @Override
            public Component getListCellRendererComponent(final JList<?> list, final Object value, final int index,
                    final boolean is_selected, final boolean has_focus) {
                super.getListCellRendererComponent(list, value, index, is_selected, has_focus);
                if (value instanceof String) {
                    setText(PropertyColorScheme.displayName((String) value));
                }
                return this;
            }
        });
        _color_by_property_cb.addActionListener(this);
        add(label);
        add(_color_by_property_cb);
    }

    /** The "Size by:" dropdown: scale each tip symbol by the value of a chosen NUMERIC phyloXML property (the size
     *  counterpart of "Color by:"). Only numeric refs appear -- see {@link #populateSizeByPropertyBox()}. */
    void setupSizeByProperty() {
        final JLabel label = new JLabel("Size by:");
        label.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _size_by_property_cb = new JComboBox<String>();
        _size_by_property_cb.setFont(ControlPanel.js_font);
        _size_by_property_cb.setToolTipText("scale the tip symbols by the value of a numeric phyloXML property");
        _size_by_property_cb.addItem(COLOR_BY_PROPERTY_NONE);
        _size_by_property_cb.setRenderer(new DefaultListCellRenderer() {
            @Override
            public Component getListCellRendererComponent(final JList<?> list, final Object value, final int index,
                    final boolean is_selected, final boolean has_focus) {
                super.getListCellRendererComponent(list, value, index, is_selected, has_focus);
                if (value instanceof String) {
                    setText(PropertyColorScheme.displayName((String) value));
                }
                return this;
            }
        });
        _size_by_property_cb.addActionListener(this);
        add(label);
        add(_size_by_property_cb);
    }

    /** Repopulate the "Size by:" dropdown from the currently displayed tree's NUMERIC properties. */
    void populateSizeByPropertyBox() {
        if (_size_by_property_cb == null) {
            return;
        }
        final TreePanel tp = getMainPanel().getCurrentTreePanel();
        _size_by_property_cb.removeActionListener(this);
        _size_by_property_cb.removeAllItems();
        _size_by_property_cb.addItem(COLOR_BY_PROPERTY_NONE);
        if ((tp != null) && (tp.getPhylogeny() != null)) {
            for (final String ref : PropertyColorScheme.numericRefs(tp.getPhylogeny())) {
                _size_by_property_cb.addItem(ref);
            }
            _size_by_property_cb.setSelectedItem(
                    (tp.getPropertySizeScale() != null) ? tp.getPropertySizeScale().getRef() : COLOR_BY_PROPERTY_NONE);
        }
        _size_by_property_cb.addActionListener(this);
    }

    /** Resets the "Size by" dropdown to None (for Reset to Defaults); the per-tab scale is cleared on the TreePanel. */
    void setSizeByPropertySelectionToNone() {
        if (_size_by_property_cb != null) {
            _size_by_property_cb.setSelectedItem(COLOR_BY_PROPERTY_NONE);
        }
    }

    /** The "Ancestral pie:" dropdown: draw an ancestral-state pie chart at each node for a chosen BEAST discrete/
     *  geographic trait. Only appears for trees that carry such traits -- see {@link #populateAncestralPieBox()}. */
    void setupAncestralPieProperty() {
        _ancestral_pie_label = new JLabel("Ancestral pie:");
        _ancestral_pie_label.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            _ancestral_pie_label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _ancestral_pie_property_cb = new JComboBox<String>();
        _ancestral_pie_property_cb.setFont(ControlPanel.js_font);
        _ancestral_pie_property_cb
                .setToolTipText("show an ancestral-state pie chart at each node for a discrete/geographic trait");
        _ancestral_pie_property_cb.addItem(COLOR_BY_PROPERTY_NONE);
        _ancestral_pie_property_cb.addActionListener(this);
        add(_ancestral_pie_label);
        add(_ancestral_pie_property_cb);
        // start hidden: no tree is loaded yet, and most trees carry no discrete-trait data (populate reveals it)
        _ancestral_pie_label.setVisible(false);
        _ancestral_pie_property_cb.setVisible(false);
    }

    /** Repopulate the "Ancestral pie:" dropdown from the displayed tree's discrete-trait properties, and show/hide
     *  the whole control (label + combo, so the row collapses) depending on whether the tree has any such trait. */
    void populateAncestralPieBox() {
        if (_ancestral_pie_property_cb == null) {
            return;
        }
        final TreePanel tp = getMainPanel().getCurrentTreePanel();
        _ancestral_pie_property_cb.removeActionListener(this);
        _ancestral_pie_property_cb.removeAllItems();
        _ancestral_pie_property_cb.addItem(COLOR_BY_PROPERTY_NONE);
        boolean has_traits = false;
        if ((tp != null) && (tp.getPhylogeny() != null)) {
            for (final String trait : TreePanelUtil.ancestralStateTraits(tp.getPhylogeny())) {
                _ancestral_pie_property_cb.addItem(trait);
                has_traits = true;
            }
            _ancestral_pie_property_cb.setSelectedItem(
                    (tp.getAncestralPieTrait() != null) ? tp.getAncestralPieTrait() : COLOR_BY_PROPERTY_NONE);
        }
        _ancestral_pie_property_cb.addActionListener(this);
        final boolean changed = _ancestral_pie_label.isVisible() != has_traits;
        _ancestral_pie_label.setVisible(has_traits);
        _ancestral_pie_property_cb.setVisible(has_traits);
        if (changed) {
            revalidate();
            repaint();
        }
    }

    /** Resets the "Ancestral pie" dropdown to None (for Reset to Defaults); the per-tab trait is cleared on the TreePanel. */
    void setAncestralPieSelectionToNone() {
        if (_ancestral_pie_property_cb != null) {
            _ancestral_pie_property_cb.setSelectedItem(COLOR_BY_PROPERTY_NONE);
        }
    }

    /** Test hook: the trait refs currently offered by the "Ancestral pie:" dropdown (excluding the "None" entry). */
    java.util.List<String> ancestralPieTraitRefs() {
        final java.util.List<String> refs = new java.util.ArrayList<>();
        if (_ancestral_pie_property_cb != null) {
            for (int i = 0; i < _ancestral_pie_property_cb.getItemCount(); ++i) {
                final String item = _ancestral_pie_property_cb.getItemAt(i);
                if (!COLOR_BY_PROPERTY_NONE.equals(item)) {
                    refs.add(item);
                }
            }
        }
        return refs;
    }

    /** Test hook: whether the "Ancestral pie:" control (dropdown) is currently visible. */
    boolean isAncestralPieControlVisible() {
        return (_ancestral_pie_property_cb != null) && _ancestral_pie_property_cb.isVisible();
    }

    /** For tests: the "Ancestral pie" dropdown's selected item as a string ("None" when off / not built). */
    String getAncestralPieSelection() {
        if (_ancestral_pie_property_cb == null) {
            return COLOR_BY_PROPERTY_NONE;
        }
        final Object sel = _ancestral_pie_property_cb.getSelectedItem();
        return (sel == null) ? COLOR_BY_PROPERTY_NONE : sel.toString();
    }

    /** The "Branch lengths:" dropdown (Time / Divergence): lay an Auspice/Nextstrain tree out by time (num_date) or by
     *  divergence (nextstrain:div) -- a reversible display mode. Only appears for a tree that carries BOTH signals --
     *  see {@link #populateBranchLengthsControl()}. */
    void setupBranchLengthsControl() {
        _branch_lengths_label = new JLabel("Branch lengths:");
        _branch_lengths_label.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            _branch_lengths_label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _branch_lengths_cb = new JComboBox<String>();
        _branch_lengths_cb.setFont(ControlPanel.js_font);
        _branch_lengths_cb.setToolTipText(
                "lay the tree out by time (num_date) or by divergence (nextstrain:div) -- a reversible display mode");
        _branch_lengths_cb.addItem(TreePanel.NEXTSTRAIN_BRANCH_MODE.TIME.label());
        _branch_lengths_cb.addItem(TreePanel.NEXTSTRAIN_BRANCH_MODE.DIVERGENCE.label());
        _branch_lengths_cb.addActionListener(this);
        add(_branch_lengths_label);
        add(_branch_lengths_cb);
        // start hidden: revealed by populate only for an Auspice/Nextstrain tree carrying both a date and a div
        _branch_lengths_label.setVisible(false);
        _branch_lengths_cb.setVisible(false);
    }

    /** Reseed the "Branch lengths:" dropdown from the current tree and show/hide the whole control (label + combo, so
     *  the row collapses) depending on whether the tree carries both a time and a divergence signal. */
    void populateBranchLengthsControl() {
        if (_branch_lengths_cb == null) {
            return;
        }
        final TreePanel tp = getMainPanel().getCurrentTreePanel();
        final boolean applicable = (tp != null) && (tp.getPhylogeny() != null)
                && tp.isNextstrainTimeDivergenceApplicable();
        _branch_lengths_cb.removeActionListener(this);
        if (applicable) {
            _branch_lengths_cb.setSelectedItem(tp.getNextstrainBranchMode().label());
        }
        _branch_lengths_cb.addActionListener(this);
        final boolean changed = _branch_lengths_label.isVisible() != applicable;
        _branch_lengths_label.setVisible(applicable);
        _branch_lengths_cb.setVisible(applicable);
        if (changed) {
            revalidate();
            repaint();
        }
    }

    /** Reseed the "Branch lengths" dropdown to the Time default WITHOUT firing the listener (for Reset to Defaults; the
     *  per-tab TreePanel model is reset separately via {@code resetNextstrainBranchModeToDefault}). */
    void setBranchLengthsSelectionToTime() {
        if (_branch_lengths_cb != null) {
            _branch_lengths_cb.removeActionListener(this);
            _branch_lengths_cb.setSelectedItem(TreePanel.NEXTSTRAIN_BRANCH_MODE.TIME.label());
            _branch_lengths_cb.addActionListener(this);
        }
    }

    /** Test hook: whether the "Branch lengths:" control is currently visible. */
    boolean isBranchLengthsControlVisible() {
        return (_branch_lengths_cb != null) && _branch_lengths_cb.isVisible();
    }

    /** Test hook: the "Branch lengths" dropdown's selected item as a string. */
    String getBranchLengthsSelection() {
        if (_branch_lengths_cb == null) {
            return TreePanel.NEXTSTRAIN_BRANCH_MODE.TIME.label();
        }
        final Object sel = _branch_lengths_cb.getSelectedItem();
        return (sel == null) ? TreePanel.NEXTSTRAIN_BRANCH_MODE.TIME.label() : sel.toString();
    }

    /** Test hook: select a "Branch lengths" mode the way a user would (fires the real actionPerformed dispatch). */
    void userSelectBranchLengthsForTest(final TreePanel.NEXTSTRAIN_BRANCH_MODE mode) {
        if (_branch_lengths_cb != null) {
            _branch_lengths_cb.setSelectedItem(mode.label());
        }
    }

    /** Repopulate the "Color by:" dropdown from the currently displayed tree's properties. */
    void populateColorByPropertyBox() {
        if (_color_by_property_cb == null) {
            return;
        }
        final TreePanel tp = getMainPanel().getCurrentTreePanel();
        _color_by_property_cb.removeActionListener(this);
        _color_by_property_cb.removeAllItems();
        _color_by_property_cb.addItem(COLOR_BY_PROPERTY_NONE);
        if ((tp != null) && (tp.getPhylogeny() != null)) {
            for (final String ref : PropertyColorScheme.colorableRefs(tp.getPhylogeny())) {
                _color_by_property_cb.addItem(ref);
            }
            // reflect the tree panel's current state
            _color_by_property_cb.setSelectedItem(
                    (tp.getPropertyColorScheme() != null) ? tp.getPropertyColorScheme().getRef()
                            : COLOR_BY_PROPERTY_NONE);
        }
        _color_by_property_cb.addActionListener(this);
    }

    /** Demo gallery: colour the current tree by a property, reflected in the "Color by" dropdown. Populates the box,
     *  selects the ref in it, and drives the render. For a demo tree the ref is always a colourable option (so the
     *  dropdown reflects it); the explicit setter is a fallback ensuring the render even if the combo did not fire. */
    void demoSelectColorByProperty(final String ref) {
        populateColorByPropertyBox();
        if (_color_by_property_cb != null) {
            _color_by_property_cb.setSelectedItem(ref);
        }
        final TreePanel tp = getMainPanel().getCurrentTreePanel();
        if (tp != null) {
            tp.setColorByPropertyRef(ref); // ensure the render even if the combo selection did not change/fire
            tp.repaint();
        }
    }

    /** Demo gallery: turn on ancestral-state pies for a discrete trait, reflected in the "Ancestral pie" dropdown. */
    void demoSelectAncestralPie(final String trait) {
        populateAncestralPieBox();
        if (_ancestral_pie_property_cb != null) {
            _ancestral_pie_property_cb.setSelectedItem(trait);
        }
        final TreePanel tp = getMainPanel().getCurrentTreePanel();
        if (tp != null) {
            tp.setAncestralPieTrait(trait);
            tp.repaint();
        }
    }

    void setupControls() {
        setupThemeButtons();
        nextRowGap(SECTION_GAP); // more space between the Light/Dark row and the P/A/C row
        setupTreeDisplayTypeOptions();
        nextRowGap(SECTION_GAP); // more space between the P/A/C row and "Color by"
        setupColorByProperty();
        setupSizeByProperty();
        setupAncestralPieProperty();
        setupBranchLengthsControl();
        setupDisplayCheckboxes();
        // Click-to options
        nextRowGap(SECTION_GAP);
        startClickToOptions();
        setupClickToOptions();
        endClickToOptions();
        // Zoom and quick edit buttons
        addButtons();
        nextRowGap(SECTION_GAP);
        setupSearchOptions();
        setupSearchTools0();
        nextRowGap(TIGHT_GAP); // less space between Search (A) and Search (B)
        setupSearchTools1();
        nextRowGap(TIGHT_GAP);
        setupSearchCombine(); // adaptive "Combine A & B" control (hidden until both boxes have a query)
        setupSearchNavigation();
        addControlPanelGlue();
    }

    // ---- Search options ---------------------------------------------------------------------
    // The redesigned search tool: each search box (A/B) has its own FIELD selector (what to look at) and
    // MODE selector (how to compare) -- see setupSearchTools0/1. Only two modifiers stay as checkboxes here,
    // because they apply to BOTH boxes and are orthogonal to any field/mode: "Match Case" and "Inverse".
    // (The old "Words"/"Regex" folded into the mode selector; "Properties" into per-property field entries;
    // "Visible" was retired.) The two checkboxes drive the Options state that the searches read.
    private JCheckBox _search_case_sensitive_cb;
    private JCheckBox _search_inverse_cb;
    private JComboBox<SearchField> _search_field_0;
    private JComboBox<SearchField> _search_field_1;
    private JComboBox<SearchMode> _search_mode_0;
    private JComboBox<SearchMode> _search_mode_1;
    private JTextField _search_range_tf_0;    // the RANGE upper bound; shown only in range mode
    private JTextField _search_range_tf_1;
    private JPanel _search_range_panel_0;     // wraps the range field so it collapses to no space when hidden
    private JPanel _search_range_panel_1;
    // value-autocomplete popups: offer the selected field's existing values as you type in the value box
    private SearchValueAutocomplete _search_autocomplete_0;
    private SearchValueAutocomplete _search_autocomplete_1;
    // "Combine A & B" control: 0 = independent (dual-highlight, default) / 1 = AND (A ∩ B) / 2 = OR (A ∪ B). Adaptive:
    // the panel is shown only when BOTH boxes have a query, so it costs no panel space for the common single search.
    private JComboBox<String> _search_combine_combo;
    private JPanel            _search_combine_panel;
    private boolean           _search_was_combined; // the last search produced a combined (not dual) result
    private Phylogeny _search_fields_tree;    // identity guard: the tree the field lists were last built for
    private List<String> _search_fields_sig;  // last-built field signatures (label + numeric-ness), to skip a no-op
    // set while re-seeding the search combos programmatically (e.g. Reset to Defaults / per-tree rebuild), so
    // their listeners don't fire a spurious search mid-adjust.
    private boolean _search_controls_adjusting;
    // remember-last (in-session): each box keeps its chosen field + its last string / numeric mode, so switching
    // field KIND (which repopulates the mode combo) or navigating to another tree doesn't reset the user's choice.
    private String _last_field_label_0;
    private String _last_field_label_1;
    private SearchMode _last_string_mode_0 = SearchMode.CONTAINS;
    private SearchMode _last_string_mode_1 = SearchMode.CONTAINS;
    private SearchMode _last_numeric_mode_0 = SearchMode.EQ;
    private SearchMode _last_numeric_mode_1 = SearchMode.EQ;

    void setupSearchOptions() {
        final JLabel header = new JLabel("Search Options:");
        header.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            header.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        add(header);
        _search_case_sensitive_cb = new JCheckBox(MainFrame.SEARCH_CASE_SENSITIVE_LABEL);
        _search_inverse_cb = new JCheckBox(MainFrame.INVERSE_SEARCH_RESULT_LABEL);
        _search_case_sensitive_cb.setToolTipText("search is case sensitive");
        _search_inverse_cb.setToolTipText("select the nodes that do NOT match");
        final Options o = getOptions();
        if (o != null) {
            _search_case_sensitive_cb.setSelected(o.isSearchCaseSensitive());
            _search_inverse_cb.setSelected(o.isInverseSearchResult());
        }
        final ActionListener l = new ActionListener() {
            @Override
            public void actionPerformed(final ActionEvent e) {
                searchOptionChanged();
            }
        };
        for (final JCheckBox cb : new JCheckBox[] { _search_case_sensitive_cb, _search_inverse_cb }) {
            cb.setFocusPainted(false);
            cb.setFont(ControlPanel.jcb_font);
            cb.setMargin(new Insets(0, 0, 0, 0)); // trim vertical padding so the checkboxes pack tightly
            if (_configuration.isApplyCustomGuiColors()) {
                cb.setBackground(getConfiguration().getGuiBackgroundColor());
                cb.setForeground(getConfiguration().getGuiCheckboxTextColor());
            }
            cb.addActionListener(l);
        }
        nextRowGap(CHECKBOX_GAP);
        add(searchOptionsRow(_search_case_sensitive_cb, _search_inverse_cb));
    }

    /** A field selector (what a query is matched against) for one search box. Seeded with the static text fields;
     *  {@link #rebuildSearchFields(boolean)} then tailors it to the loaded tree (only present fields + one entry per
     *  custom property + numeric built-ins). "Any text field" (index 0) is the default = search-everything. */
    private JComboBox<SearchField> makeSearchFieldCombo(final boolean box_a) {
        final JComboBox<SearchField> combo = new JComboBox<>(SearchField.stringMenuFields());
        combo.setSelectedIndex(0); // "Any text field"
        combo.setToolTipText("which node data to search");
        styleSearchCombo(combo);
        combo.addActionListener(e -> {
            if (!_search_controls_adjusting) {
                onSearchFieldChanged(box_a);
            }
        });
        return combo;
    }

    /** A match-mode selector (how a query is compared) for one search box. Holds the string modes or the numeric
     *  operators depending on the selected field's type (see {@link #reconcileModeCombo(boolean)}); the default is
     *  the first entry ("contains" for text, "=" for numbers). */
    private JComboBox<SearchMode> makeSearchModeCombo(final boolean box_a) {
        final JComboBox<SearchMode> combo = new JComboBox<>(SearchMode.stringModes());
        combo.setSelectedIndex(0); // "contains"
        combo.setToolTipText("how to match");
        styleSearchCombo(combo);
        combo.addActionListener(e -> {
            if (!_search_controls_adjusting) {
                onSearchModeChanged(box_a);
            }
        });
        return combo;
    }

    /** Reacts to a field-selector change: switch the mode set to match the field's type, show/hide the range
     *  field, then re-run that box's search. */
    private void onSearchFieldChanged(final boolean box_a) {
        dismissAutocomplete(box_a); // the value set / eligibility changed with the field; recompute on next open
        final SearchField f = (SearchField) (box_a ? _search_field_0 : _search_field_1).getSelectedItem();
        if (f != null) { // remember this box's field so a later per-tree rebuild re-selects it when available
            if (box_a) {
                _last_field_label_0 = f.label();
            }
            else {
                _last_field_label_1 = f.label();
            }
        }
        reconcileModeCombo(box_a);
        updateRangeFieldVisibility(box_a);
        if (box_a) {
            search0();
        }
        else {
            search1();
        }
    }

    /** Reacts to a mode-selector change: remember it (per box, per kind), show/hide the range field, re-run. */
    private void onSearchModeChanged(final boolean box_a) {
        dismissAutocomplete(box_a); // regex mode suppresses the popup; other modes keep it -- recompute on next open
        final SearchField f = (SearchField) (box_a ? _search_field_0 : _search_field_1).getSelectedItem();
        final SearchMode m = (SearchMode) (box_a ? _search_mode_0 : _search_mode_1).getSelectedItem();
        if ((f != null) && (m != null)) {
            rememberMode(box_a, f.isNumeric(), m);
        }
        updateRangeFieldVisibility(box_a);
        if (box_a) {
            search0();
        }
        else {
            search1();
        }
    }

    private SearchMode rememberedMode(final boolean box_a, final boolean numeric) {
        if (box_a) {
            return numeric ? _last_numeric_mode_0 : _last_string_mode_0;
        }
        return numeric ? _last_numeric_mode_1 : _last_string_mode_1;
    }

    private void rememberMode(final boolean box_a, final boolean numeric, final SearchMode m) {
        if (box_a) {
            if (numeric) {
                _last_numeric_mode_0 = m;
            }
            else {
                _last_string_mode_0 = m;
            }
        }
        else if (numeric) {
            _last_numeric_mode_1 = m;
        }
        else {
            _last_string_mode_1 = m;
        }
    }

    /** Makes the box's mode combo hold the numeric operators when its field is numeric, else the string modes --
     *  repopulating (and resetting to the first mode) only when the kind actually flips. Runs under the adjusting
     *  guard so it fires no search. */
    private void reconcileModeCombo(final boolean box_a) {
        final JComboBox<SearchField> fc = box_a ? _search_field_0 : _search_field_1;
        final JComboBox<SearchMode> mc = box_a ? _search_mode_0 : _search_mode_1;
        if ((fc == null) || (mc == null)) {
            return;
        }
        final SearchField f = (SearchField) fc.getSelectedItem();
        if (f == null) {
            return;
        }
        final SearchMode current = (SearchMode) mc.getSelectedItem();
        final boolean is_numeric_now = (current != null) && current.isNumeric();
        if ((mc.getItemCount() > 0) && (f.isNumeric() == is_numeric_now)) {
            return; // already the right set
        }
        final boolean was_adjusting = _search_controls_adjusting;
        _search_controls_adjusting = true;
        try {
            mc.removeAllItems();
            for (final SearchMode m : (f.isNumeric() ? SearchMode.numericModes() : SearchMode.stringModes())) {
                mc.addItem(m);
            }
            // restore this box's last-used mode of the new kind (remember-last), not always the first entry
            mc.setSelectedItem(rememberedMode(box_a, f.isNumeric()));
        }
        finally {
            _search_controls_adjusting = was_adjusting;
        }
    }

    /** Shows the box's range upper-bound field only when its mode is {@link SearchMode#RANGE}. */
    private void updateRangeFieldVisibility(final boolean box_a) {
        final JComboBox<SearchMode> mc = box_a ? _search_mode_0 : _search_mode_1;
        final JPanel panel = box_a ? _search_range_panel_0 : _search_range_panel_1;
        if ((mc == null) || (panel == null)) {
            return;
        }
        final boolean show = (mc.getSelectedItem() == SearchMode.RANGE);
        if (panel.isVisible() != show) {
            panel.setVisible(show);
            revalidate();
            repaint();
        }
    }

    /**
     * Tailors the two field selectors to the currently displayed tree (see {@link SearchField#availableFields}).
     * Identity-guarded: a {@code false} force is a no-op when the tree hasn't changed (so it costs nothing on the
     * per-search repaint path); pass {@code true} after a data edit on the SAME tree (e.g. importing annotations).
     * Preserves the user's selection by label, and only repopulates the combos when the field list actually
     * differs -- runs under the adjusting guard so it launches no search.
     */
    void rebuildSearchFields(final boolean force) {
        if ((_search_field_0 == null) || (_search_field_1 == null) || (_mainpanel == null)) {
            return;
        }
        final Phylogeny phy = _mainpanel.getCurrentPhylogeny();
        if (!force && (phy == _search_fields_tree)) {
            return;
        }
        final List<SearchField> fields = SearchField.availableFields(phy);
        // the signature includes each field's KIND (numeric vs string), not just its label, so a custom property
        // that flips numeric<->string (e.g. re-import changes its values) repopulates rather than keeping a stale
        // field of the wrong kind.
        final List<String> sig = new ArrayList<String>();
        for (final SearchField f : fields) {
            sig.add(f.label() + "\t" + f.isNumeric());
        }
        _search_fields_tree = phy;
        if (sig.equals(_search_fields_sig)) {
            return; // same fields (label AND kind) -> keep the current selection untouched
        }
        _search_fields_sig = sig;
        _search_controls_adjusting = true;
        try {
            repopulateFieldCombo(_search_field_0, fields, true);
            repopulateFieldCombo(_search_field_1, fields, false);
        }
        finally {
            _search_controls_adjusting = false;
        }
    }

    private void repopulateFieldCombo(final JComboBox<SearchField> combo, final List<SearchField> fields,
                                      final boolean box_a) {
        // prefer this box's remembered field (remember-last), falling back to the current selection; so switching
        // to a tree that lacks it and back restores it rather than sticking at "Any text".
        final SearchField prev = (SearchField) combo.getSelectedItem();
        String want_label = box_a ? _last_field_label_0 : _last_field_label_1;
        if ((want_label == null) && (prev != null)) {
            want_label = prev.label();
        }
        combo.removeAllItems();
        int select = 0;
        for (int i = 0; i < fields.size(); i++) {
            combo.addItem(fields.get(i));
            if (fields.get(i).label().equals(want_label)) {
                select = i;
            }
        }
        if (combo.getItemCount() > 0) {
            combo.setSelectedIndex(select);
        }
        reconcileModeCombo(box_a);          // the (re)selected field may be a different kind than before
        updateRangeFieldVisibility(box_a);
    }

    /** The range upper-bound text field for a box; typing in it re-runs that box's (numeric range) search. */
    private JTextField makeRangeField(final boolean box_a) {
        final JTextField tf = new JTextField(3);
        tf.setFont(ControlPanel.jcb_font);
        tf.setToolTipText("range upper bound");
        installTextUndo(tf);
        if (getConfiguration().isApplyCustomGuiColors()) {
            tf.setForeground(getConfiguration().getGuiMenuBackgroundColor());
            tf.setBackground(getConfiguration().getGuiCheckboxTextColor());
            tf.setBorder(null);
        }
        tf.addKeyListener(new KeyAdapter() {

            @Override
            public void keyReleased(final KeyEvent e) {
                if (isNonEditingKey(e.getKeyCode())) {
                    return; // caret-navigation keys don't change the range bound -> don't re-run the search
                }
                if (box_a) {
                    search0();
                }
                else {
                    search1();
                }
            }
        });
        return tf;
    }

    /** Wraps a range field with a "to" label; the whole row is hidden (and takes no space) unless in range mode. */
    private JPanel makeRangePanel(final JTextField range_tf) {
        final JPanel p = new JPanel(new BorderLayout(2, 0));
        p.setBackground(getBackground());
        final JLabel to = new JLabel("to ");
        to.setFont(ControlPanel.jcb_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            to.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        p.add(to, BorderLayout.WEST);
        p.add(range_tf, BorderLayout.CENTER);
        p.setVisible(false);
        return p;
    }

    /** Shared styling for the search field/mode combos: the small control-panel font, a label renderer (so the
     *  combo shows each field's/mode's friendly label), and a tiny preferred width so the full-width GridBag row
     *  stretches them without the widest item forcing the whole panel wider. */
    private void styleSearchCombo(final JComboBox<?> combo) {
        combo.setFont(ControlPanel.jcb_font);
        combo.setRenderer(SEARCH_LABEL_RENDERER);
        if (_configuration.isApplyCustomGuiColors()) {
            // match the other control-panel combos (Click-on / sequence-relation / color-by), not the text fields
            combo.setBackground(getConfiguration().getGuiButtonBackgroundColor());
            combo.setForeground(getConfiguration().getGuiButtonTextColor());
        }
        combo.setPreferredSize(new Dimension(10, combo.getPreferredSize().height));
    }

    /** Renders a {@link SearchField}/{@link SearchMode} combo entry by its friendly {@code label()}. */
    private static final DefaultListCellRenderer SEARCH_LABEL_RENDERER = new DefaultListCellRenderer() {

        @Override
        public Component getListCellRendererComponent(final JList<?> list, final Object value, final int index,
                                                      final boolean is_selected, final boolean has_focus) {
            final Object text = (value instanceof SearchField) ? ((SearchField) value).label()
                    : (value instanceof SearchMode) ? ((SearchMode) value).label() : value;
            return super.getListCellRendererComponent(list, text, index, is_selected, has_focus);
        }
    };

    // ---- test hooks for the search tool -------------------------------------------------------
    int searchFieldCountForTest() {
        return ( _search_field_0 != null ) ? _search_field_0.getItemCount() : 0;
    }

    int searchModeCountForTest() {
        return ( _search_mode_0 != null ) ? _search_mode_0.getItemCount() : 0;
    }

    SearchField getSearchFieldForTest(final boolean box_a) {
        return (SearchField) ( box_a ? _search_field_0 : _search_field_1 ).getSelectedItem();
    }

    SearchMode getSearchModeForTest(final boolean box_a) {
        return (SearchMode) ( box_a ? _search_mode_0 : _search_mode_1 ).getSelectedItem();
    }

    /** Selects the field whose backing NDF equals {@code ndf} (or the "Any text" field when {@code ndf} is null),
     *  without firing a search (guarded), so a test can set field + text + then call search0()/search1() itself. */
    void setSearchFieldForTest(final boolean box_a, final PhylogenyMethods.NDF ndf) {
        final JComboBox<SearchField> combo = box_a ? _search_field_0 : _search_field_1;
        _search_controls_adjusting = true;
        try {
            for ( int i = 0; i < combo.getItemCount(); i++ ) {
                final SearchField f = combo.getItemAt( i );
                final boolean match = ( ndf == null ) ? ( f.kind() == SearchField.Kind.ANY_TEXT ) : ( f.ndf() == ndf );
                if ( match ) {
                    combo.setSelectedIndex( i );
                    return;
                }
            }
        }
        finally {
            _search_controls_adjusting = false;
        }
    }

    void setSearchModeForTest(final boolean box_a, final SearchMode mode) {
        _search_controls_adjusting = true;
        try {
            ( box_a ? _search_mode_0 : _search_mode_1 ).setSelectedItem( mode );
            updateRangeFieldVisibility( box_a );
        }
        finally {
            _search_controls_adjusting = false;
        }
    }

    List<String> searchFieldLabelsForTest(final boolean box_a) {
        final JComboBox<SearchField> c = box_a ? _search_field_0 : _search_field_1;
        final List<String> out = new ArrayList<String>();
        for ( int i = 0; i < c.getItemCount(); i++ ) {
            out.add( c.getItemAt( i ).label() );
        }
        return out;
    }

    /** Selects a field by its label and reconciles the mode combo to its type (numeric vs string), without firing
     *  a search -- the test then sets the mode + value and calls search0()/search1() itself. */
    void setSearchFieldByLabelForTest(final boolean box_a, final String label) {
        final JComboBox<SearchField> c = box_a ? _search_field_0 : _search_field_1;
        _search_controls_adjusting = true;
        try {
            boolean found = false;
            for ( int i = 0; i < c.getItemCount(); i++ ) {
                if ( c.getItemAt( i ).label().equals( label ) ) {
                    c.setSelectedIndex( i );
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                throw new IllegalArgumentException( "no such search field: " + label );
            }
            reconcileModeCombo( box_a );
            updateRangeFieldVisibility( box_a );
        }
        finally {
            _search_controls_adjusting = false;
        }
    }

    boolean rangeFieldVisibleForTest(final boolean box_a) {
        final JPanel p = box_a ? _search_range_panel_0 : _search_range_panel_1;
        return ( p != null ) && p.isVisible();
    }

    void setRangeHighForTest(final boolean box_a, final String text) {
        ( box_a ? _search_range_tf_0 : _search_range_tf_1 ).setText( text );
    }

    /** Simulates a USER field selection (fires the real listener, so remember-last records it), for tests. */
    void userSelectFieldForTest(final boolean box_a, final String label) {
        final JComboBox<SearchField> c = box_a ? _search_field_0 : _search_field_1;
        for ( int i = 0; i < c.getItemCount(); i++ ) {
            if ( c.getItemAt( i ).label().equals( label ) ) {
                c.setSelectedIndex( i );
                return;
            }
        }
        throw new IllegalArgumentException( "no such search field: " + label );
    }

    /** Simulates a USER mode selection (fires the real listener, so remember-last records it), for tests. */
    void userSelectModeForTest(final boolean box_a, final SearchMode m) {
        ( box_a ? _search_mode_0 : _search_mode_1 ).setSelectedItem( m );
    }

    void setSearchCaseSensitiveForTest(final boolean b) {
        _search_case_sensitive_cb.setSelected( b );
    }

    void setSearchInverseForTest(final boolean b) {
        _search_inverse_cb.setSelected( b );
    }

    /** A cheap repaint of just the tree canvas after a search: the found-set highlight is the only thing that
     *  changed, so a search must NOT run the full layout recalc of displayedPhylogenyMightHaveChanged(true) (which,
     *  per keystroke on a large tree, re-measures every label -- the dominant search-lag source). */
    private void repaintTreeAfterSearch() {
        if ((_mainpanel != null) && (_mainpanel.getCurrentTreePanel() != null)) {
            _mainpanel.getCurrentTreePanel().repaint();
        }
    }

    /** Draws a red "error" outline on a box's query field (and its range field) when the current query can never
     *  match: an uncompilable regular expression, or a non-numeric operand on a numeric field. Because the search
     *  runs on every keystroke (no Enter button), this passive cue -- not a dialog -- is the right signal for a
     *  query that is only transiently invalid while being typed. Uses FlatLaf's {@code JComponent.outline}. */
    private void updateSearchFieldValidity(final boolean box_a) {
        final JComboBox<SearchField> fc = box_a ? _search_field_0 : _search_field_1;
        final JComboBox<SearchMode> mc = box_a ? _search_mode_0 : _search_mode_1;
        if ((fc == null) || (mc == null)) {
            return;
        }
        final SearchField field = (SearchField) fc.getSelectedItem();
        final SearchMode mode = (SearchMode) mc.getSelectedItem();
        final JTextField value_tf = box_a ? _search_tf_0 : _search_tf_1;
        final JTextField range_tf = box_a ? _search_range_tf_0 : _search_range_tf_1;
        setOutlineError(value_tf, isSearchQueryInvalid(field, mode, (value_tf == null) ? null : value_tf.getText()));
        setOutlineError(range_tf, (mode == SearchMode.RANGE)
                && isSearchQueryInvalid(field, mode, (range_tf == null) ? null : range_tf.getText()));
    }

    /** Whether a NON-EMPTY query can never match: an invalid regex (REGEX mode), or a non-number on a numeric
     *  field. An empty query is not "invalid" -- it just clears the search. */
    private static boolean isSearchQueryInvalid(final SearchField field, final SearchMode mode, final String text) {
        if ((field == null) || (mode == null) || ForesterUtil.isEmpty(text) || ForesterUtil.isEmpty(text.trim())) {
            return false;
        }
        if (field.isNumeric()) {
            return SearchMatcher.parseFiniteDouble(text) == null;
        }
        return (mode == SearchMode.REGEX) && !SearchMatcher.isCompilableRegex(text);
    }

    private static void setOutlineError(final JTextField tf, final boolean error) {
        if (tf == null) {
            return;
        }
        final Object want = error ? "error" : null;
        if (!java.util.Objects.equals(tf.getClientProperty("JComponent.outline"), want)) {
            tf.putClientProperty("JComponent.outline", want);
            tf.repaint();
        }
    }

    private JPanel searchOptionsRow(final JCheckBox a, final JCheckBox b) {
        final JPanel p = new JPanel(new GridLayout(1, 2, 0, 0));
        if (_configuration.isApplyCustomGuiColors()) {
            p.setBackground(getBackground());
        }
        p.add(a);
        p.add(b);
        return p;
    }

    /**
     * Pushes the two shared search modifiers (Match Case, Inverse) into {@link Options}, then re-runs BOTH
     * searches and repaints so the change is reflected immediately (the old menu got the repaint for free when
     * it closed). The field/mode selectors are per-box and re-run only their own box from their own listeners.
     */
    private void searchOptionChanged() {
        final Options o = getOptions();
        o.setSearchCaseSensitive(_search_case_sensitive_cb.isSelected());
        o.setInverseSearchResult(_search_inverse_cb.isSelected());
        reRunSearches();
    }

    /**
     * A compact Light/Dark theme switch shown at the very top of the control panel,
     * just above the P/A/C display-type selector. The two radio buttons drive
     * {@link MainFrame#setDarkMode(boolean)} directly.
     */
    void setupThemeButtons() {
        _light_mode_rb = new JRadioButton("Light");
        _dark_mode_rb = new JRadioButton("Dark");
        final boolean dark = getConfiguration().getUi() == Configuration.UI.FLAT_DARK;
        _light_mode_rb.setSelected(!dark);
        _dark_mode_rb.setSelected(dark);
        final ButtonGroup group = new ButtonGroup();
        group.add(_light_mode_rb);
        group.add(_dark_mode_rb);
        _light_mode_rb.setToolTipText("light theme");
        _dark_mode_rb.setToolTipText("dark theme");
        for (final JRadioButton rb : new JRadioButton[] { _light_mode_rb, _dark_mode_rb }) {
            rb.setFocusPainted(false);
            rb.setFont(ControlPanel.jcb_font);
            if (_configuration.isApplyCustomGuiColors()) {
                rb.setBackground(getConfiguration().getGuiBackgroundColor());
                rb.setForeground(getConfiguration().getGuiCheckboxTextColor());
            }
        }
        _light_mode_rb.addActionListener(e -> getMainPanel().getMainFrame().setDarkMode(false));
        _dark_mode_rb.addActionListener(e -> getMainPanel().getMainFrame().setDarkMode(true));
        final JPanel p = new JPanel(new GridLayout(1, 2, 0, 0));
        p.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            p.setBackground(getConfiguration().getGuiBackgroundColor());
            p.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        p.add(_light_mode_rb);
        p.add(_dark_mode_rb);
        add(p);
    }

    void setupTreeDisplayTypeOptions() {
        _display_as_unaligned_phylogram_rb = new JRadioButton("P");
        _display_as_aligned_phylogram_rb = new JRadioButton("A");
        _display_as_cladogram_rb = new JRadioButton("C");
        _display_as_buttongroup = new ButtonGroup();
        _display_as_buttongroup.add(_display_as_unaligned_phylogram_rb);
        _display_as_buttongroup.add(_display_as_aligned_phylogram_rb);
        _display_as_buttongroup.add(_display_as_cladogram_rb);
        getDisplayAsUnalignedPhylogramRb().setToolTipText("(unaligned) phylogram");
        getDisplayAsAlignedPhylogramRb().setToolTipText("aligned phylogram");
        getDisplayAsCladogramRb().setToolTipText("cladogram");
        setupJRadioButton(getDisplayAsUnalignedPhylogramRb());
        setupJRadioButton(getDisplayAsAlignedPhylogramRb());
        setupJRadioButton(getDisplayAsCladogramRb());
        final JPanel p = new JPanel(new GridLayout(1, 3, 0, 0));
        p.setFont(ControlPanel.jcb_font);
        if (_configuration.isApplyCustomGuiColors()) {
            p.setBackground(getConfiguration().getGuiBackgroundColor());
            p.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        p.add(_display_as_unaligned_phylogram_rb);
        p.add(_display_as_aligned_phylogram_rb);
        p.add(_display_as_cladogram_rb);
        add(p);
    }

    void setUpControlsForDomainStrucures() {
        _domain_display_label = new JLabel("Domain Architectures:");
        add(customizeLabel(_domain_display_label, getConfiguration()));
        add(_domain_display_label);
        _zoom_in_domain_structure = new TypomaticJButton("d+");
        _zoom_out_domain_structure = new TypomaticJButton("d-");
        _decr_domain_structure_evalue_thr = new JButton("-");
        _incr_domain_structure_evalue_thr = new JButton("+");
        _zoom_in_domain_structure.setPreferredSize(new Dimension(10, 10));
        _zoom_out_domain_structure.setPreferredSize(new Dimension(10, 10));
        _decr_domain_structure_evalue_thr.setPreferredSize(new Dimension(10, 10));
        _incr_domain_structure_evalue_thr.setPreferredSize(new Dimension(10, 10));
        _incr_domain_structure_evalue_thr.setToolTipText("Increase the E-value threshold by a factor of 10");
        _decr_domain_structure_evalue_thr.setToolTipText("Decrease the E-value threshold by a factor of 10");
        _domain_structure_evalue_thr_tf = new JTextField(3);
        _domain_structure_evalue_thr_tf.setFont(ControlPanel.jcb_font);
        _domain_structure_evalue_thr_tf.setEditable(false);
        if (getConfiguration().isApplyCustomGuiColors()) {
            _domain_structure_evalue_thr_tf.setForeground(getConfiguration().getGuiMenuBackgroundColor());
            _domain_structure_evalue_thr_tf.setBackground(getConfiguration().getGuiCheckboxTextColor());
            _domain_structure_evalue_thr_tf.setBorder(null);
        }
        final JPanel d1_panel = new JPanel(new GridLayout(1, 2, 0, 0));
        final JPanel d2_panel = new JPanel(new GridLayout(1, 3, 0, 0));
        if (_configuration.isApplyCustomGuiColors()) {
            d1_panel.setBackground(getBackground());
            d2_panel.setBackground(getBackground());
        }
        add(d1_panel);
        add(d2_panel);
        addJButton(_zoom_out_domain_structure, d1_panel);
        addJButton(_zoom_in_domain_structure, d1_panel);
        addJButton(_decr_domain_structure_evalue_thr, d2_panel);
        addJTextField(_domain_structure_evalue_thr_tf, d2_panel);
        addJButton(_incr_domain_structure_evalue_thr, d2_panel);
    }


    void setupSearchTools0() {
        final JLabel search_label = new JLabel("Search (A):");
        search_label.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            search_label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        add(search_label);
        search_label.setToolTipText(SEARCH_TIP_TEXT);
        _search_field_0 = makeSearchFieldCombo(true);
        _search_mode_0 = makeSearchModeCombo(true);
        add(_search_field_0);
        add(_search_mode_0);
        _search_found_label_0 = new JLabel();
        getSearchFoundCountsLabel0().setVisible(false);
        _search_found_label_0.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            _search_found_label_0.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _search_tf_0 = new JTextField(3);
        _search_tf_0.setFont(ControlPanel.jcb_font);
        _search_tf_0.setToolTipText(SEARCH_TIP_TEXT);
        installTextUndo(_search_tf_0);
        _search_tf_0.setEditable(true);
        if (getConfiguration().isApplyCustomGuiColors()) {
            _search_tf_0.setForeground(getConfiguration().getGuiMenuBackgroundColor());
            _search_tf_0.setBackground(getConfiguration().getGuiCheckboxTextColor());
            _search_tf_0.setBorder(null);
        }
        _search_reset_button_0 = new JButton();
        getSearchResetButton0().setText("Reset");
        getSearchResetButton0().setEnabled(false);
        getSearchResetButton0().setVisible(false);
        _search_range_tf_0 = makeRangeField(true);
        _search_range_panel_0 = makeRangePanel(_search_range_tf_0);
        final JPanel s_panel_1 = new JPanel(new BorderLayout());
        final JPanel s_panel_2 = new JPanel(new GridLayout(1, 2, 0, 0));
        s_panel_1.setBackground(getBackground());
        add(s_panel_1);
        add(_search_range_panel_0); // hidden unless the mode is "range"
        s_panel_2.setBackground(getBackground());
        add(s_panel_2);
        final KeyAdapter key_adapter = new KeyAdapter() {

            @Override
            public void keyReleased(final KeyEvent key_event) {
                if (isNonEditingKey(key_event.getKeyCode())) {
                    return; // navigation / autocomplete-popup keys don't change the query -> don't re-run the search
                }
                search0();
            }
        };
        final ActionListener action_listener = new ActionListener() {

            @Override
            public void actionPerformed(final ActionEvent e) {
                getSearchTextField0().setText("");
                _search_range_tf_0.setText("");
                // route through search0() so the empty box clears its found set + labels AND the combine state is
                // reconciled (restore box B's independent highlight, drop _search_was_combined, hide the Combine row)
                search0();
            }
        };
        _search_reset_button_0.addActionListener(action_listener);
        _search_tf_0.addKeyListener(key_adapter);
        _search_autocomplete_0 = new SearchValueAutocomplete(_search_tf_0, () -> autocompleteValues(true),
                () -> (_search_case_sensitive_cb != null) && _search_case_sensitive_cb.isSelected(), this::search0);
        addJTextField(_search_tf_0, s_panel_1);
        s_panel_2.add(_search_found_label_0);
        addJButton(_search_reset_button_0, s_panel_2);
    }

    void setupSearchTools1() {
        final JLabel search_label = new JLabel("Search (B):");
        search_label.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            search_label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        add(search_label);
        search_label.setToolTipText(SEARCH_TIP_TEXT);
        _search_field_1 = makeSearchFieldCombo(false);
        _search_mode_1 = makeSearchModeCombo(false);
        add(_search_field_1);
        add(_search_mode_1);
        _search_found_label_1 = new JLabel();
        getSearchFoundCountsLabel1().setVisible(false);
        _search_found_label_1.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            _search_found_label_1.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _search_tf_1 = new JTextField(3);
        _search_tf_1.setFont(ControlPanel.jcb_font);
        _search_tf_1.setToolTipText(SEARCH_TIP_TEXT);
        installTextUndo(_search_tf_1);
        _search_tf_1.setEditable(true);
        if (getConfiguration().isApplyCustomGuiColors()) {
            _search_tf_1.setForeground(getConfiguration().getGuiMenuBackgroundColor());
            _search_tf_1.setBackground(getConfiguration().getGuiCheckboxTextColor());
            _search_tf_1.setBorder(null);
        }
        _search_reset_button_1 = new JButton();
        getSearchResetButton1().setText("Reset");
        getSearchResetButton1().setEnabled(false);
        getSearchResetButton1().setVisible(false);
        _search_range_tf_1 = makeRangeField(false);
        _search_range_panel_1 = makeRangePanel(_search_range_tf_1);
        final JPanel s_panel_1 = new JPanel(new BorderLayout());
        final JPanel s_panel_2 = new JPanel(new GridLayout(1, 2, 0, 0));
        s_panel_1.setBackground(getBackground());
        add(s_panel_1);
        add(_search_range_panel_1); // hidden unless the mode is "range"
        s_panel_2.setBackground(getBackground());
        add(s_panel_2);
        final KeyAdapter key_adapter = new KeyAdapter() {

            @Override
            public void keyReleased(final KeyEvent key_event) {
                if (isNonEditingKey(key_event.getKeyCode())) {
                    return; // navigation / autocomplete-popup keys don't change the query -> don't re-run the search
                }
                search1();
            }
        };
        final ActionListener action_listener = new ActionListener() {

            @Override
            public void actionPerformed(final ActionEvent e) {
                getSearchTextField1().setText("");
                _search_range_tf_1.setText("");
                // route through search1() so the empty box clears its found set + labels AND the combine state is
                // reconciled (restore box A's independent highlight, drop _search_was_combined, hide the Combine row)
                search1();
            }
        };
        _search_reset_button_1.addActionListener(action_listener);
        _search_tf_1.addKeyListener(key_adapter);
        _search_autocomplete_1 = new SearchValueAutocomplete(_search_tf_1, () -> autocompleteValues(false),
                () -> (_search_case_sensitive_cb != null) && _search_case_sensitive_cb.isSelected(), this::search1);
        addJTextField(_search_tf_1, s_panel_1);
        s_panel_2.add(_search_found_label_1);
        addJButton(_search_reset_button_1, s_panel_2);
    }

    /** The "◀  k / N  ▶" step-through row: jumps the viewport from hit to hit. Hidden until a search finds
     *  something ({@link #updateSearchHitNavigation()} shows it). Mirrors the View-menu Find Next/Previous. */
    void setupSearchNavigation() {
        _search_prev_button = makeSearchNavButton("◀", "Center the previous search hit in the view (⌘⇧G)",
                -1);
        _search_next_button = makeSearchNavButton("▶", "Center the next search hit in the view (⌘G)", 1);
        _search_nav_label = new JLabel("", javax.swing.SwingConstants.CENTER);
        _search_nav_label.setFont(ControlPanel.jcb_bold_font);
        _search_nav_label.setToolTipText("Position among the search hits");
        if (getConfiguration().isApplyCustomGuiColors()) {
            _search_nav_label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _search_nav_panel = new JPanel(new BorderLayout());
        _search_nav_panel.setBackground(getBackground());
        _search_nav_panel.add(_search_prev_button, BorderLayout.WEST);
        _search_nav_panel.add(_search_nav_label, BorderLayout.CENTER);
        _search_nav_panel.add(_search_next_button, BorderLayout.EAST);
        _search_nav_panel.setVisible(false);
        add(_search_nav_panel);
    }

    private JButton makeSearchNavButton(final String glyph, final String tip, final int dir) {
        final JButton b = new JButton(glyph);
        b.setFocusPainted(false);
        b.setFont(ControlPanel.jcb_bold_font);
        b.setMargin(new Insets(2, 10, 2, 10)); // roomy enough to match the sibling zoom buttons' height
        b.setToolTipText(tip);
        if (getConfiguration().isApplyCustomGuiColors()) {
            b.setBorder(BorderFactory.createLineBorder(getConfiguration().getGuiButtonBorderColor()));
            b.setBackground(getConfiguration().getGuiButtonBackgroundColor());
            b.setForeground(getConfiguration().getGuiButtonTextColor());
        }
        b.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(final ActionEvent e) {
                stepSearchHit(dir);
            }
        });
        return b;
    }

    /** Test hook: is the "k / N" step-through row currently shown? */
    boolean isSearchNavVisibleForTest() {
        return (_search_nav_panel != null) && _search_nav_panel.isVisible();
    }

    /** Test hook: the current "k / N" navigator text (empty string if the row was never built). */
    String getSearchNavLabelForTest() {
        return (_search_nav_label == null) ? "" : _search_nav_label.getText();
    }

    /** Test hook: the navigator row's laid-out height (0 when hidden/collapsed by the GridBag). Call after forcing
     *  a layout pass -- the revalidate() in updateSearchHitNavigation is otherwise applied asynchronously. */
    int getSearchNavPanelHeightForTest() {
        return (_search_nav_panel == null) ? 0 : _search_nav_panel.getHeight();
    }

    /** Advances to the next (dir=+1) or previous (dir=-1) search hit and re-labels the navigator. */
    void stepSearchHit(final int dir) {
        final TreePanel tp = getCurrentTreePanel();
        if (tp != null) {
            tp.stepToFoundNode(dir);
        }
    }

    /** Refreshes the "k / N" navigator (label text + row visibility) from the current tree's hit set. Cheap when
     *  there is no active search (no hits -> row hidden). Revalidates only when the row actually appears/disappears
     *  -- a bare setVisible() on this GridBag row would otherwise not be laid out until some later event (mirrors
     *  updateDataCheckboxVisibility). */
    void updateSearchHitNavigation() {
        if (_search_nav_panel == null) {
            return;
        }
        final TreePanel tp = getCurrentTreePanel();
        final int count = (tp == null) ? 0 : tp.getSearchHitCount();
        final boolean show = count > 0;
        if (show) {
            final int idx = tp.getSearchHitIndex();
            // clamp: an out-of-band found-set shrink can leave the index past the end until the next step
            final String pos = ((idx < 0) || (idx >= count)) ? "–" : String.valueOf(idx + 1);
            _search_nav_label.setText(pos + " / " + count);
        }
        if (_search_nav_panel.isVisible() != show) {
            _search_nav_panel.setVisible(show);
            revalidate();
            repaint();
        }
        // keep the prominent menu-bar "Found / Selected: N" counter in sync (this method is the single choke point
        // every found-set change funnels through: search, reset, manual selection, prune, undo, tab change)
        final MainFrame mf = getMainFrame();
        if (mf != null) {
            mf.updateFoundSelectedCounter();
        }
    }

    void setVisibilityOfX() {
        final MainFrame mf = getMainFrame();
        if (mf != null) {
            if ((getCurrentTreePanel() != null) && (getCurrentTreePanel().getPhylogeny() != null)) {
                if (AptxUtil.isHasAtLeastOneBranchWithSupportSD(getCurrentTreePanel().getPhylogeny())) {
                    if (mf._show_confidence_stddev_cbmi != null) {
                        mf._show_confidence_stddev_cbmi.setVisible(true);
                    }
                } else {
                    if (mf._show_confidence_stddev_cbmi != null) {
                        mf._show_confidence_stddev_cbmi.setVisible(false);
                    }
                }
                if (AptxUtil.isHasAtLeastOneNodeWithScientificName(getCurrentTreePanel().getPhylogeny())) {
                    if (mf._abbreviate_scientific_names != null) {
                        mf._abbreviate_scientific_names.setVisible(true);
                    }
                } else {
                    if (mf._abbreviate_scientific_names != null) {
                        mf._abbreviate_scientific_names.setVisible(false);
                    }
                }
            }
            if (isDrawPhylogram() || ((getCurrentTreePanel() != null) && ((getCurrentTreePanel()
                    .getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                    || (getCurrentTreePanel().getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)))) {
                if (mf._non_lined_up_cladograms_rbmi != null) {
                    mf._non_lined_up_cladograms_rbmi.setVisible(false);
                }
                if (mf._ext_node_dependent_cladogram_rbmi != null) {
                    mf._ext_node_dependent_cladogram_rbmi.setVisible(false);
                }
            } else {
                if (mf._non_lined_up_cladograms_rbmi != null) {
                    mf._non_lined_up_cladograms_rbmi.setVisible(true);
                }
                if (mf._ext_node_dependent_cladogram_rbmi != null) {
                    mf._ext_node_dependent_cladogram_rbmi.setVisible(true);
                }
            }
            if (isDrawPhylogram()) {
                if (mf._show_scale_cbmi != null) {
                    mf._show_scale_cbmi.setVisible(true);
                }
            } else {
                if (mf._show_scale_cbmi != null) {
                    mf._show_scale_cbmi.setVisible(false);
                }
            }
            if (getCurrentTreePanel() != null) {
                if ((getCurrentTreePanel().getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                        || (getCurrentTreePanel().getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
                    if (mf._label_direction_cbmi != null) {
                        mf._label_direction_cbmi.setVisible(true);
                    }
                } else {
                    if (mf._label_direction_cbmi != null) {
                        mf._label_direction_cbmi.setVisible(false);
                    }
                }
            }
        }
    }

    /**
     * Fit entire tree into window.
     */
    void showWhole() {
        if ((_mainpanel.getCurrentScrollPane() == null)
                || _mainpanel.getCurrentTreePanel().getPhylogeny().isEmpty()) {
            return;
        }
        getCurrentTreePanel().updateSetOfCollapsedExternalNodes();
        displayedPhylogenyMightHaveChanged(true);
        _mainpanel.getCurrentTreePanel().updateOvSettings();
        _mainpanel.getCurrentTreePanel().validate();
        _mainpanel.validate();
        _mainpanel.getCurrentTreePanel().calcParametersForPainting(_mainpanel.getSizeOfViewport().width,
                _mainpanel.getSizeOfViewport().height);
        // radial "fit to window": re-fit the square radial canvas to the viewport (zoom back to whole)
        final PHYLOGENY_GRAPHICS_TYPE gt = _mainpanel.getCurrentTreePanel().getPhylogenyGraphicsType();
        if ((gt == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) || (gt == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
            _mainpanel.getCurrentTreePanel().fitRadialTo(_mainpanel.getSizeOfViewport().width,
                    _mainpanel.getSizeOfViewport().height);
        }
        _mainpanel.getCurrentTreePanel().resetPreferredSize();
        _mainpanel.adjustJScrollPane();
        _mainpanel.getCurrentTreePanel().repaint();
        _mainpanel.getCurrentTreePanel().validate();
        _mainpanel.validate();
        _mainpanel.getCurrentTreePanel().calcParametersForPainting(_mainpanel.getSizeOfViewport().width,
                _mainpanel.getSizeOfViewport().height);
        if ((gt == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) || (gt == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
            _mainpanel.getCurrentTreePanel().fitRadialTo(_mainpanel.getSizeOfViewport().width,
                    _mainpanel.getSizeOfViewport().height);
        }
        _mainpanel.getCurrentTreePanel().resetPreferredSize();
        _mainpanel.adjustJScrollPane();
        _mainpanel.getCurrentTreePanel().repaint();
        _mainpanel.getCurrentTreePanel().updateOvSizes();
    }

    /**
     * Fits the tree to the viewport in the horizontal direction only (so e.g. newly-shown domain
     * architectures become visible), leaving the vertical zoom untouched -- unlike {@link #showWhole()},
     * which re-fits both directions. Works by passing the tree's current vertical extent as the height
     * to {@code calcParametersForPainting()}, which derives the y-distance from it, so the y-distance
     * is left unchanged while the x-direction is re-fitted to the viewport width.
     */
    void fitWidth() {
        if ((_mainpanel.getCurrentScrollPane() == null)
                || _mainpanel.getCurrentTreePanel().getPhylogeny().isEmpty()) {
            return;
        }
        final TreePanel tp = getCurrentTreePanel();
        tp.updateSetOfCollapsedExternalNodes();
        displayedPhylogenyMightHaveChanged(true); // recalc longest-ext-node info (e.g. now with domains)
        tp.updateOvSettings();
        tp.validate();
        _mainpanel.validate();
        tp.resetPreferredSize();
        final int keep_height = tp.getPreferredSize().height; // the current vertical extent
        tp.calcParametersForPainting(_mainpanel.getSizeOfViewport().width, keep_height);
        tp.resetPreferredSize();
        _mainpanel.adjustJScrollPane();
        tp.repaint();
        tp.validate();
        _mainpanel.validate();
        // second pass: a vertical scrollbar may have appeared/disappeared, changing the usable width
        tp.calcParametersForPainting(_mainpanel.getSizeOfViewport().width, keep_height);
        tp.resetPreferredSize();
        _mainpanel.adjustJScrollPane();
        tp.repaint();
        tp.updateOvSizes();
    }

    /** The vertical-orientation analog of {@link #fitWidth()}: fit the depth (branch-length) axis -- which is drawn
     *  vertically in a root-top/bottom tree -- to the window HEIGHT, keeping the current horizontal (tip-spread)
     *  zoom. (calcParametersForPainting swaps the depth/breadth budgets in a vertical orientation, so this is the
     *  transpose of fitWidth.) */
    void fitHeight() {
        if ((_mainpanel.getCurrentScrollPane() == null)
                || _mainpanel.getCurrentTreePanel().getPhylogeny().isEmpty()) {
            return;
        }
        final TreePanel tp = getCurrentTreePanel();
        tp.updateSetOfCollapsedExternalNodes();
        displayedPhylogenyMightHaveChanged(true);
        tp.updateOvSettings();
        tp.validate();
        _mainpanel.validate();
        tp.resetPreferredSize();
        // the current breadth (tip-spread) extent WITHOUT the tilted-label diagonal allowance that resetPreferredSize
        // adds to the preferred width -- feeding the padded width back would grow the breadth zoom on every press
        final int keep_width = tp.logicalBreadthExtent();
        tp.calcParametersForPainting(keep_width, _mainpanel.getSizeOfViewport().height);
        tp.resetPreferredSize();
        _mainpanel.adjustJScrollPane();
        tp.repaint();
        tp.validate();
        _mainpanel.validate();
        // second pass: a horizontal scrollbar may have appeared/disappeared, changing the usable height
        tp.calcParametersForPainting(keep_width, _mainpanel.getSizeOfViewport().height);
        tp.resetPreferredSize();
        _mainpanel.adjustJScrollPane();
        tp.repaint();
        tp.updateOvSizes();
    }

    private boolean isVerticalOrientation() {
        return (getCurrentTreePanel() != null) && getCurrentTreePanel().isVerticalOrientation();
    }

    private boolean isRadialLayout() {
        return (getCurrentTreePanel() != null) && getCurrentTreePanel().isRadialLayout();
    }

    // A zoom re-centers the scroll bar of the SCREEN axis it changes. The depth (x) axis is drawn horizontally
    // normally but VERTICALLY in a vertical orientation, and the breadth (y) axis is the other way round -- so the
    // zoom methods must pick the scroll bar by orientation, else zooming re-centers the wrong axis (or not at all).
    private JScrollBar depthScrollBar() {
        return isVerticalOrientation() ? getMainPanel().getCurrentScrollPane().getVerticalScrollBar()
                : getMainPanel().getCurrentScrollPane().getHorizontalScrollBar();
    }

    private JScrollBar breadthScrollBar() {
        return isVerticalOrientation() ? getMainPanel().getCurrentScrollPane().getHorizontalScrollBar()
                : getMainPanel().getCurrentScrollPane().getVerticalScrollBar();
    }

    // Screen-oriented zoom: the buttons/keyboard/wheel are labeled by SCREEN direction (X = horizontal, Y =
    // vertical). In a vertical orientation the tree's depth/breadth axes are swapped on screen, so "screen X"
    // (horizontal) drives the tip-spread (the zoomY methods) and "screen Y" (vertical) drives the depth (zoomX).
    // The correction factor is only consumed by the depth (zoomX) path.
    void zoomInScreenX(final float factor, final float x_correction_factor) {
        if (isVerticalOrientation()) {
            zoomInY(factor);
        } else {
            zoomInX(factor, x_correction_factor);
        }
    }

    void zoomOutScreenX(final float factor, final float x_correction_factor) {
        if (isVerticalOrientation()) {
            zoomOutY(factor);
        } else {
            zoomOutX(factor, x_correction_factor);
        }
    }

    void zoomInScreenY(final float factor, final float x_correction_factor) {
        if (isVerticalOrientation()) {
            zoomInX(factor, x_correction_factor);
        } else {
            zoomInY(factor);
        }
    }

    void zoomOutScreenY(final float factor, final float x_correction_factor) {
        if (isVerticalOrientation()) {
            zoomOutX(factor, x_correction_factor);
        } else {
            zoomOutY(factor);
        }
    }

    /** The current label of the fit-width/height button ("W" horizontal, "H" vertical). For tests. */
    String getFitButtonText() {
        return (_fit_width == null) ? null : _fit_width.getText();
    }

    // Test hooks for the layout-aware zoom-cluster relabeling (see updateZoomButtonsForLayout): expose the buttons
    // so a test can read their text/enabled state and drive a real doClick gesture through actionPerformed.
    JButton getZoomInXButtonForTest()  { return _zoom_in_x; }
    JButton getZoomOutXButtonForTest() { return _zoom_out_x; }
    JButton getZoomInYButtonForTest()  { return _zoom_in_y; }
    JButton getZoomOutYButtonForTest() { return _zoom_out_y; }
    JButton getExpandButtonForTest()   { return _expand_y; }
    JButton getFitWidthButtonForTest() { return _fit_width; }

    /** Keeps the zoom-cluster buttons in sync with the current LAYOUT. Three cases:
     *  <ul>
     *  <li><b>Radial</b> (circular/unrooted): the two zoom axes collapse to one (both drive the single radial
     *  diameter since 0.11.9), so Y+/Y- become a plain +/- zoom and the now-redundant X-/X+ become rotate
     *  counter-clockwise / clockwise (the same {@code setStartingAngle} path as the S/A keys). "E" (vertical
     *  label spacing) does nothing in a fan so it is greyed out, and "W" (fit-width) -- redundant with "F" (fit
     *  &amp; center) radially -- becomes a node-label-direction flip (horizontal &lt;-&gt; radial, the "Radial
     *  Labels" setting).</li>
     *  <li><b>Vertical</b> (root-top/bottom rectangular): "W" becomes "H" (fit-height) and the E tooltip
     *  describes the now-horizontal tip spacing; the X/Y labels stay screen-relative (their actions flip at
     *  click time via the screen-oriented zoom helpers).</li>
     *  <li><b>Horizontal</b> (default rectangular): the plain X+/X-/Y+/Y-/E/W labels.</li>
     *  </ul>
     *  Call whenever the layout may have changed (panel build, tab change, Type change, Settings apply, Reset). */
    void updateZoomButtonsForLayout() {
        final boolean radial = isRadialLayout();
        final boolean vertical = isVerticalOrientation();
        final boolean native_ui = getConfiguration().isUseNativeUI();
        if (_zoom_in_y != null) {
            _zoom_in_y.setText(radial ? "+" : "Y+");
            _zoom_in_y.setToolTipText(radial ? "zoom in [mouse wheel, +, or Alt+Up]"
                    : "zoom in vertically [Alt+Up or Shift+mousewheel]");
        }
        if (_zoom_out_y != null) {
            _zoom_out_y.setText(radial ? "-" : "Y-");
            _zoom_out_y.setToolTipText(radial ? "zoom out [mouse wheel, -, or Alt+Down]"
                    : "zoom out vertically [Alt+Down or Shift+mousewheel]");
        }
        if (_zoom_in_x != null) {
            _zoom_in_x.setText(radial ? rotateLabel(_zoom_in_x, ROTATE_CW_LABEL, ROTATE_CW_FALLBACK)
                    : (native_ui ? "+" : "X+"));
            _zoom_in_x.setToolTipText(radial ? "rotate clockwise [S or Shift+mousewheel]"
                    : "zoom in horizontally [Alt+Right or Shift+Alt+mousewheel]");
        }
        if (_zoom_out_x != null) {
            _zoom_out_x.setText(radial ? rotateLabel(_zoom_out_x, ROTATE_CCW_LABEL, ROTATE_CCW_FALLBACK)
                    : (native_ui ? "-" : "X-"));
            _zoom_out_x.setToolTipText(radial ? "rotate counter-clockwise [A or Shift+mousewheel]"
                    : "zoom out horizontally [Alt+Left or Shift+Alt+mousewheel]");
        }
        if (_expand_y != null) {
            _expand_y.setEnabled(!radial);
            _expand_y.setToolTipText(radial
                    ? "expand to fit labels -- not available in radial (circular/unrooted) display"
                    : (vertical
                            ? "expand the tree horizontally so labels do not overlap at the current font size [Alt+E]"
                            : "expand the tree in vertical direction so labels do not overlap at the current font size [Alt+E]"));
        }
        if (_fit_width != null) {
            if (radial) {
                _fit_width.setText(LABEL_DIRECTION_BUTTON_LABEL);
                _fit_width.setToolTipText(labelDirectionButtonTooltip());
            }
            else {
                _fit_width.setText(vertical ? "H" : "W");
                _fit_width.setToolTipText(vertical
                        ? "fit the tree to the window height, keeping the current horizontal zoom [Alt+W]"
                        : "fit the tree to the window width, keeping the current vertical zoom [Alt+W]");
            }
        }
    }

    /** The circular-arrow glyph if the button's font can render it, else a plain-text fallback (so a minimal-font
     *  platform never shows a missing-glyph box for the radial rotate buttons). */
    private static String rotateLabel(final JButton button, final String glyph, final String fallback) {
        return ((button.getFont() != null) && button.getFont().canDisplay(glyph.charAt(0))) ? glyph : fallback;
    }

    /** The state-reflecting tooltip for the radial "L" (label-direction flip) button. */
    private String labelDirectionButtonTooltip() {
        final boolean radial_labels = getOptions()
                .getNodeLabelDirection() == Options.NODE_LABEL_DIRECTION.RADIAL;
        return "node labels: " + (radial_labels ? "radial" : "horizontal") + " -- click to flip [Alt+W]";
    }

    /** Flips the node-label direction (horizontal &lt;-&gt; radial) from the radial "L" button, reusing the same
     *  Settings-checkbox path so the option, its menu/dialog state, persistence and Reset stay in sync, then
     *  refreshes the button's state tooltip. Package-visible so the Alt+W keyboard shortcut (TreePanel) drives the
     *  same flip in a radial layout, matching the button. */
    void toggleNodeLabelDirection() {
        getMainPanel().getMainFrame().toggleRadialLabelDirection();
        updateZoomButtonsForLayout();
    }

    void showWholeAll() {
        for (final TreePanel tree_panel : _mainpanel.getTreePanels()) {
            if (tree_panel != null) {
                tree_panel.validate();
                tree_panel.calcParametersForPainting(_mainpanel.getSizeOfViewport().width,
                        _mainpanel.getSizeOfViewport().height);
                tree_panel.resetPreferredSize();
                tree_panel.repaint();
            }
        }
    }

    // Create header for click-to combo box.
    void startClickToOptions() {
        final JLabel spacer = new JLabel("");
        spacer.setFont(ControlPanel.jcb_font);
        add(spacer);
        _click_to_label = new JLabel("Click on Node to:");
        add(customizeLabel(_click_to_label, getConfiguration()));
        _click_to_combobox = new JComboBox<String>();
        _click_to_combobox.setFocusable(false);
        _click_to_combobox.setMaximumRowCount(14);
        _click_to_combobox.setFont(ControlPanel.js_font);
        if (_configuration.isApplyCustomGuiColors()) {
            _click_to_combobox.setBackground(getConfiguration().getGuiBackgroundColor());
        }
        // don't add listener until all items are set (or each one will trigger
        // an event)
        // click_to_list.addActionListener(this);
        add(_click_to_combobox);
        _click_to_names = new ArrayList<String>();
    }

    void tabChanged() {
        if (getMainPanel().getTabbedPane().getTabCount() > 0) {
            // P/A/C (phylogram/cladogram) apply in every layout that can honor branch lengths -- CIRCULAR renders a
            // real phylogram since 0.11.7 (isCircularPhylogram), UNROOTED always has -- so enable on branch-length
            // presence alone and PRESERVE the tab's chosen display type (the old "&& != CIRCULAR" branch force-disabled
            // the radios AND silently forced CLADOGRAM whenever a branch-length tree was viewed circular).
            if (getCurrentTreePanel().isPhyHasBranchLengths()) {
                setDrawPhylogramEnabled(true);
                setTreeDisplayType(getTreeDisplayType(getMainPanel().getCurrentTabIndex()));
            } else {
                setDrawPhylogramEnabled(false);
                setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM);
            }
            getMainPanel().getMainFrame()
                    .setSelectedTypeInTypeMenu(getMainPanel().getCurrentTreePanel().getPhylogenyGraphicsType());
            // re-seed a left-open modeless Settings dialog's per-tab controls (tree style, palette, Time Axis) so it
            // reflects the now-current tab -- AFTER the Type-menu radios above are synced (the tree-style reseeder
            // reads them)
            getMainPanel().getMainFrame().refreshOpenSettingsDialog();
            updateZoomButtonsForLayout(); // relabel the zoom cluster for the current (possibly persisted) layout
            getMainPanel().getCurrentTreePanel().updateSubSuperTreeButton();
            getMainPanel().getCurrentTreePanel().updateButtonToUncollapseAll();
            getMainPanel().getControlPanel().updateDomainStructureEvaluethresholdDisplay();
            updateDataCheckboxVisibility(true);
            populateColorByPropertyBox();
            populateSizeByPropertyBox();
            populateAncestralPieBox();
            populateBranchLengthsControl();
            // run the searches AFTER the field selectors are rebuilt for this tab's tree, so the highlight reflects
            // the field the combo now shows (not the previous tab's field).
            reRunSearches();
            if (getMainPanel().getMainFrame() != null) {
                getMainPanel().getMainFrame().updateEditMenu(); // undo history is per-tab
            }
            getMainPanel().getCurrentTreePanel().repaint();
        }
        else if (getMainPanel().getMainFrame() != null) {
            getMainPanel().getMainFrame().updateEditMenu(); // last tab closed -> disable Undo/Redo
        }
    }

    /**
     * Uncollapse all nodes.
     */
    final void uncollapseAll(final TreePanel tp) {
        final Phylogeny t = tp.getPhylogeny();
        if ((t != null) && !t.isEmpty()) {
            for (final PhylogenyNodeIterator iter = t.iteratorPreorder(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
                node.setCollapse(false);
            }
            tp.resetNodeIdToDistToLeafMap();
            tp.updateSetOfCollapsedExternalNodes();
            t.recalculateNumberOfExternalDescendants(false);
            tp.setNodeInPreorderToNull();
            t.clearHashIdToNodeMap();
            showWhole();
        }
    }

    final void updateDomainStructureEvaluethresholdDisplay() {
        if (_domain_structure_evalue_thr_tf != null) {
            _domain_structure_evalue_thr_tf.setText(
                    "10" + toSuperscript(getMainPanel().getCurrentTreePanel().getDomainStructureEvalueThresholdExp()));
        }
    }

    /** Renders an integer exponent using Unicode superscript glyphs, e.g. -3 becomes "⁻³". */
    private static String toSuperscript(final int exponent) {
        final char[] superscript_digits = { '⁰', '¹', '²', '³', '⁴', '⁵', '⁶',
                '⁷', '⁸', '⁹' };
        final StringBuilder sb = new StringBuilder();
        if (exponent < 0) {
            sb.append('⁻'); // superscript minus
        }
        for (final char digit : Integer.toString(Math.abs(exponent)).toCharArray()) {
            sb.append(superscript_digits[digit - '0']);
        }
        return sb.toString();
    }

    /** Radial zoom: scale the single radial-canvas diameter (decoupled from the rectangular x/y-distance), keeping the
     *  relative scroll position on both axes so the zoom stays roughly centred. Returns true when it handled a radial
     *  layout (so the caller returns), false otherwise. Every radial zoom gesture (X/Y buttons + wheel) routes here. */
    private boolean zoomRadial(final TreePanel tp, final float factor) {
        final PHYLOGENY_GRAPHICS_TYPE t = tp.getPhylogenyGraphicsType();
        if ((t != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) && (t != PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
            return false;
        }
        final JScrollBar hsb = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar();
        final JScrollBar vsb = getMainPanel().getCurrentScrollPane().getVerticalScrollBar();
        final double hc = (hsb.getMaximum() - hsb.getMinimum()) / (hsb.getValue() + (hsb.getVisibleAmount() / 2.0));
        final double vc = (vsb.getMaximum() - vsb.getMinimum()) / (vsb.getValue() + (vsb.getVisibleAmount() / 2.0));
        tp.multiplyRadialDiameter(factor);
        // recompute the layout params against the NEW diameter -- the domain-width cap (effectiveDomainStructureWidth),
        // the domain/label reservation and the label truncation are all radialDiameter-based and are set only here, so
        // without this an incremental +/- zoom leaves them stale (domains sized for the old diameter -> clip/misalign)
        tp.calcParametersForPainting(getMainPanel().getSizeOfViewport().width, getMainPanel().getSizeOfViewport().height);
        getMainPanel().adjustJScrollPane();
        tp.resetPreferredSize();
        getMainPanel().getCurrentScrollPane().getViewport().validate();
        if (hc > 0) {
            hsb.setValue(ForesterUtil.roundToInt(((hsb.getMaximum() - hsb.getMinimum()) / hc)
                    - (hsb.getVisibleAmount() / 2.0)));
        }
        if (vc > 0) {
            vsb.setValue(ForesterUtil.roundToInt(((vsb.getMaximum() - vsb.getMinimum()) / vc)
                    - (vsb.getVisibleAmount() / 2.0)));
        }
        tp.resetPreferredSize();
        tp.updateOvSizes();
        return true;
    }

    final void zoomInX(final float factor, final float x_correction_factor) {
        final JScrollBar sb = depthScrollBar();
        final TreePanel treepanel = getMainPanel().getCurrentTreePanel();
        if (zoomRadial(treepanel, factor)) {
            return;
        }
        treepanel.multiplyUrtFactor(1f);
        if ((treepanel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                || (treepanel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                || isDrawPhylogram(getMainPanel().getCurrentTabIndex())
                || (getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP)) {
            final double x = (sb.getMaximum() - sb.getMinimum())
                    / (sb.getValue() + (sb.getVisibleAmount() / 2.0));
            treepanel.setXdistance((treepanel.getXdistance() * factor));
            treepanel.setXcorrectionFactor((treepanel.getXcorrectionFactor() * x_correction_factor));
            getMainPanel().adjustJScrollPane();
            treepanel.resetPreferredSize();
            getMainPanel().getCurrentScrollPane().getViewport().validate();
            sb.setValue(ForesterUtil
                    .roundToInt(((sb.getMaximum() - sb.getMinimum()) / x) - (sb.getVisibleAmount() / 2.0)));
        } else {
            final int x = sb.getMaximum() - sb.getMinimum() - sb.getVisibleAmount() - sb.getValue();
            treepanel.setXdistance((treepanel.getXdistance() * factor));
            treepanel.setXcorrectionFactor((treepanel.getXcorrectionFactor() * x_correction_factor));
            getMainPanel().adjustJScrollPane();
            treepanel.resetPreferredSize();
            getMainPanel().getCurrentScrollPane().getViewport().validate();
            sb.setValue(sb.getMaximum() - sb.getMinimum() - x - sb.getVisibleAmount());
        }
        treepanel.resetPreferredSize();
        treepanel.updateOvSizes();
    }

    /**
     * Expands the tree vertically just enough that leaf labels no longer overlap at the current
     * font size (so dynamic hiding is no longer needed and the "Dyna Hide" indicator clears on
     * the next repaint). Only ever expands -- if the tree is already spaced out enough, it is left
     * alone. Has no effect for the circular/unrooted layouts, where there are no vertical rows.
     */
    final void expandYToFitLabels() {
        final TreePanel treepanel = getMainPanel().getCurrentTreePanel();
        if ((treepanel == null) || (treepanel.getPhylogeny() == null) || treepanel.getPhylogeny().isEmpty()) {
            return;
        }
        final PHYLOGENY_GRAPHICS_TYPE t = treepanel.getPhylogenyGraphicsType();
        if ((t == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) || (t == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
            return;
        }
        final float target = TreePanelUtil.yDistanceToAvoidLabelOverlap(treepanel.getLargeFontHeight());
        if (treepanel.getYdistance() >= target) {
            return;
        }
        // expand grows the tip-spread (breadth) axis, which is drawn HORIZONTALLY in a vertical orientation, so
        // re-center the breadth scroll bar (vertical normally, horizontal when vertical) -- like the zoom methods
        final JScrollBar sb = breadthScrollBar();
        final double center = (sb.getMaximum() - sb.getMinimum()) / (sb.getValue() + (sb.getVisibleAmount() / 2.0));
        treepanel.setYdistance(target);
        getMainPanel().adjustJScrollPane();
        treepanel.resetPreferredSize();
        getMainPanel().getCurrentScrollPane().getViewport().validate();
        if (center > 0) {
            sb.setValue(ForesterUtil
                    .roundToInt(((sb.getMaximum() - sb.getMinimum()) / center) - (sb.getVisibleAmount() / 2.0)));
        }
        treepanel.resetPreferredSize();
        treepanel.updateOvSizes();
        treepanel.repaint();
    }

    final void zoomInY(final float factor) {
        final JScrollBar sb = breadthScrollBar();
        final TreePanel treepanel = getMainPanel().getCurrentTreePanel();
        if (zoomRadial(treepanel, factor)) {
            return;
        }
        treepanel.multiplyUrtFactor(1.1f);
        final double x = (sb.getMaximum() - sb.getMinimum()) / (sb.getValue() + (sb.getVisibleAmount() / 2.0));
        treepanel.setYdistance((treepanel.getYdistance() * factor));
        getMainPanel().adjustJScrollPane();
        treepanel.resetPreferredSize();
        getMainPanel().getCurrentScrollPane().getViewport().validate();
        sb.setValue(ForesterUtil
                .roundToInt(((sb.getMaximum() - sb.getMinimum()) / x) - (sb.getVisibleAmount() / 2.0)));
        treepanel.resetPreferredSize();
        treepanel.updateOvSizes();
    }

    final void zoomOutX(final float factor, final float x_correction_factor) {
        final TreePanel treepanel = getMainPanel().getCurrentTreePanel();
        if (zoomRadial(treepanel, factor)) {
            return;
        }
        treepanel.multiplyUrtFactor(1f);
        if ((treepanel.getXdistance() * factor) > 0.0) {
            final JScrollBar sb = depthScrollBar();
            if ((treepanel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                    || (treepanel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                    || isDrawPhylogram(getMainPanel().getCurrentTabIndex())
                    || (getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP)) {
                getMainPanel().adjustJScrollPane();
                treepanel.resetPreferredSize();
                getMainPanel().getCurrentScrollPane().getViewport().validate();
                final double x = (sb.getMaximum() - sb.getMinimum())
                        / (sb.getValue() + (sb.getVisibleAmount() / 2.0));
                treepanel.setXdistance((treepanel.getXdistance() * factor));
                treepanel.setXcorrectionFactor((treepanel.getXcorrectionFactor() * x_correction_factor));
                getMainPanel().adjustJScrollPane();
                treepanel.resetPreferredSize();
                getMainPanel().getCurrentScrollPane().getViewport().validate();
                sb.setValue(ForesterUtil.roundToInt(((sb.getMaximum() - sb.getMinimum()) / x)
                        - (sb.getVisibleAmount() / 2.0)));
            } else {
                final int x = sb.getMaximum() - sb.getMinimum() - sb.getVisibleAmount() - sb.getValue();
                treepanel.setXdistance(treepanel.getXdistance() * factor);
                treepanel.setXcorrectionFactor(treepanel.getXcorrectionFactor() * x_correction_factor);
                if (x > 0) {
                    getMainPanel().adjustJScrollPane();
                    treepanel.resetPreferredSize();
                    getMainPanel().getCurrentScrollPane().getViewport().validate();
                    sb.setValue(sb.getMaximum() - sb.getMinimum() - x - sb.getVisibleAmount());
                }
            }
            treepanel.resetPreferredSize();
            treepanel.updateOvSizes();
        }
    }

    private final boolean isDrawPhylogram(int currentTabIndex) {
        Options.PHYLOGENY_DISPLAY_TYPE t = getTreeDisplayType(currentTabIndex);
        return ((t == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM)
                | (t == Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM));
    }

    final void zoomOutY(final float factor) {
        final TreePanel treepanel = getMainPanel().getCurrentTreePanel();
        if (zoomRadial(treepanel, factor)) {
            return;
        }
        treepanel.multiplyUrtFactor(0.9f);
        if ((treepanel.getYdistance() * factor) > 0.0) {
            final JScrollBar sb = breadthScrollBar();
            final double x = (sb.getMaximum() - sb.getMinimum())
                    / (sb.getValue() + (sb.getVisibleAmount() / 2.0));
            treepanel.setYdistance((treepanel.getYdistance() * factor));
            getMainPanel().adjustJScrollPane();
            treepanel.resetPreferredSize();
            getMainPanel().getCurrentScrollPane().getViewport().validate();
            sb.setValue(ForesterUtil
                    .roundToInt(((sb.getMaximum() - sb.getMinimum()) / x) - (sb.getVisibleAmount() / 2.0)));
            treepanel.resetPreferredSize();
            treepanel.updateOvSizes();
        }
    }

    final static JLabel customizeLabel(final JLabel label, final Configuration configuration) {
        label.setFont(ControlPanel.jcb_bold_font);
        if (configuration.isApplyCustomGuiColors()) {
            label.setForeground(configuration.getGuiCheckboxTextColor());
            label.setBackground(configuration.getGuiBackgroundColor());
        }
        return label;
    }

    final public JCheckBox getUseBranchWidthsCb() {
        return _width_branches;
    }

    public Options.PHYLOGENY_DISPLAY_TYPE getTreeDisplayType() {
        if (_display_as_unaligned_phylogram_rb.isSelected()) {
            return Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM;
        } else if (_display_as_aligned_phylogram_rb.isSelected()) {
            return Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM;
        }
        return Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM;
    }
}
