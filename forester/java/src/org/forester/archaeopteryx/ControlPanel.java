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
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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
import javax.swing.ListCellRenderer;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.util.TypomaticJButton;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.SequenceRelation;
import org.forester.phylogeny.data.SequenceRelation.SEQUENCE_RELATION_TYPE;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

final class ControlPanel extends JPanel implements ActionListener {

    enum NodeClickAction {
        ADD_NEW_NODE,
        BLAST,
        COLLAPSE,
        COLOR_SUBTREE,
        COPY_SUBTREE,
        CUT_SUBTREE,
        DELETE_NODE_OR_SUBTREE,
        EDIT_NODE_DATA,
        OPEN_PDB_WEB,
        OPEN_SEQ_WEB,
        OPEN_TAX_WEB,
        PASTE_SUBTREE,
        REROOT,
        SELECT_NODES,
        SHOW_DATA,
        SORT_DESCENDENTS,
        SUBTREE,
        SWAP,
        COLOR_NODE_FONT,
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
    private Map<Integer, String> _all_click_to_names;
    private int _blast_item;
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
    private int _open_pdb_item;
    private int _open_seq_web_item;
    private int _open_tax_web_item;
    private int _color_node_font_item;
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
    private Sequence _selected_query_seq;
    private JCheckBox _seq_relation_confidence_switch;
    private JComboBox<SEQUENCE_RELATION_TYPE> _sequence_relation_type_box;
    private JCheckBox _show_binary_character_counts;
    private JCheckBox _show_binary_characters;
    // Indices for the click-to options in the combo box
    private int _show_data_item;
    private JCheckBox _show_domain_architectures;
    private JCheckBox _show_mol_seqs;
    private JCheckBox _write_branch_length_values;
    private JCheckBox _show_events;
    private JCheckBox _show_gene_names;
    private JCheckBox _show_node_names;
    private JCheckBox _shorten_labels_cb;
    private JCheckBox _show_properties_cb;
    private JCheckBox _show_seq_names;
    private JCheckBox _show_seq_symbols;
    private JCheckBox _show_sequence_acc;
    private JComboBox<String> _show_sequence_relations;
    private JCheckBox _show_taxo_code;
    private JCheckBox _show_taxo_rank;
    private JCheckBox _show_taxo_common_names;
    private JCheckBox _show_taxo_scientific_names;
    private JCheckBox _show_vector_data_cb;
    private JButton _show_whole;
    private JButton _expand_y;
    private JButton _fit_width;
    private int _sort_descendents_item;
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
    private final Map<Integer, JPanel> _checkbox_panels = new HashMap<Integer, JPanel>();
    private Set<Integer>               _data_presence;
    private Phylogeny                  _data_presence_for;
    // The data-dependent checkboxes (Node Name is intentionally NOT here: it is always shown). Their
    // visibility tracks scanForDataPresence over the whole loaded tree.
    private final static int[]         DATA_GATED_CHECKBOXES     = {
            Configuration.show_tax_code, Configuration.show_taxonomy_scientific_names,
            Configuration.show_taxonomy_common_names, Configuration.show_tax_rank, Configuration.show_seq_names,
            Configuration.show_gene_names, Configuration.show_seq_symbols, Configuration.show_sequence_acc,
            Configuration.write_confidence_values, Configuration.write_branch_length_values,
            Configuration.write_events, Configuration.show_binary_characters,
            Configuration.show_binary_character_counts, Configuration.show_domain_architectures,
            Configuration.show_vector_data, Configuration.show_properties, Configuration.use_style,
            Configuration.width_branches, Configuration.shorten_labels };

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
            } else if (e.getSource() == _click_to_combobox) {
                setClickToAction(_click_to_combobox.getSelectedIndex());
                getCurrentTreePanel().repaint();
            } else if (e.getSource() == _show_binary_characters) {
                if ((_show_binary_character_counts != null) && _show_binary_characters.isSelected()) {
                    _show_binary_character_counts.setSelected(false);
                }
                displayedPhylogenyMightHaveChanged(true);
            } else if (e.getSource() == _show_binary_character_counts) {
                if ((_show_binary_characters != null) && _show_binary_character_counts.isSelected()) {
                    _show_binary_characters.setSelected(false);
                }
                displayedPhylogenyMightHaveChanged(true);
            } else if (e.getSource() == _show_domain_architectures) {
                search0();
                search1();
                // When the user switches domains ON, re-fit the (now wider) tree horizontally so the
                // domains actually become visible -- otherwise it looks like nothing happened. The
                // horizontal-only fit keeps the user's vertical zoom; turning domains OFF just repaints.
                if (_show_domain_architectures.isSelected()) {
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
                    search0();
                    search1();
                    displayedPhylogenyMightHaveChanged(true);
                } else if (e.getSource() == _incr_domain_structure_evalue_thr) {
                    _mainpanel.getCurrentTreePanel().increaseDomainStructureEvalueThresholdExp();
                    search0();
                    search1();
                    displayedPhylogenyMightHaveChanged(true);
                } else if (e.getSource() == _search_tf_0) {
                    search0();
                    displayedPhylogenyMightHaveChanged(true);
                } else if (e.getSource() == _search_tf_1) {
                    search1();
                    displayedPhylogenyMightHaveChanged(true);
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

    public Sequence getSelectedQuerySequence() {
        return _selected_query_seq;
    }

    public JComboBox<String> getSequenceRelationBox() {
        if (_show_sequence_relations == null) {
            _show_sequence_relations = new JComboBox<String>();
            _show_sequence_relations.setFocusable(false);
            _show_sequence_relations.setMaximumRowCount(20);
            _show_sequence_relations.setFont(ControlPanel.js_font);
            if (_configuration.isApplyCustomGuiColors()) {
                _show_sequence_relations.setBackground(getConfiguration().getGuiButtonBackgroundColor());
                _show_sequence_relations.setForeground(getConfiguration().getGuiButtonTextColor());
            }
            _show_sequence_relations.addItem("-----");
            _show_sequence_relations.setToolTipText("To display orthology information for selected query");
        }
        return _show_sequence_relations;
    }

    /* GUILHEM_BEG */
    public JComboBox<SEQUENCE_RELATION_TYPE> getSequenceRelationTypeBox() {
        if (_sequence_relation_type_box == null) {
            _sequence_relation_type_box = new JComboBox<SEQUENCE_RELATION_TYPE>();
            for (final SequenceRelation.SEQUENCE_RELATION_TYPE type : SequenceRelation.SEQUENCE_RELATION_TYPE
                    .values()) {
                _sequence_relation_type_box.addItem(type);
            }
            _sequence_relation_type_box.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(final ActionEvent e) {
                    if (_mainpanel.getCurrentPhylogeny() != null) {
                        setSequenceRelationQueries(getMainPanel().getCurrentPhylogeny().getSequenceRelationQueries());
                    }
                }
            });
        }
        return _sequence_relation_type_box;
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

    public boolean isShowMolSequences() {
        return ((_show_mol_seqs != null) && _show_mol_seqs.isSelected());
    }

    public boolean isShowProperties() {
        return ((_show_properties_cb != null) && _show_properties_cb.isSelected());
    }


    public boolean isShowVectorData() {
        return ((_show_vector_data_cb != null) && _show_vector_data_cb.isSelected());
    }

    public void setSequenceRelationQueries(final Collection<Sequence> sequenceRelationQueries) {
        final JComboBox<String> box = getSequenceRelationBox();
        while (box.getItemCount() > 1) {
            box.removeItemAt(1);
        }
        final HashMap<String, Sequence> sequencesByName = new HashMap<String, Sequence>();
        final SequenceRelation.SEQUENCE_RELATION_TYPE relationType = (SequenceRelation.SEQUENCE_RELATION_TYPE) _sequence_relation_type_box
                .getSelectedItem();
        if (relationType == null) {
            return;
        }
        final ArrayList<String> sequenceNamesToAdd = new ArrayList<String>();
        for (final Sequence seq : sequenceRelationQueries) {
            if (seq.isHasSequenceRelations()) {
                boolean fFoundForCurrentType = false;
                for (final SequenceRelation sq : seq.getSequenceRelations()) {
                    if (sq.getType().equals(relationType)) {
                        fFoundForCurrentType = true;
                        break;
                    }
                }
                if (fFoundForCurrentType) {
                    sequenceNamesToAdd.add(seq.getName());
                    sequencesByName.put(seq.getName(), seq);
                }
            }
        }
        // sort sequences by name before adding them to the combo
        final String[] sequenceNameArray = sequenceNamesToAdd.toArray(new String[sequenceNamesToAdd.size()]);
        Arrays.sort(sequenceNameArray, String.CASE_INSENSITIVE_ORDER);
        for (final String seqName : sequenceNameArray) {
            box.addItem(seqName);
        }
        for (final ItemListener oldItemListener : box.getItemListeners()) {
            box.removeItemListener(oldItemListener);
        }
        box.addItemListener(new ItemListener() {

            @Override
            public void itemStateChanged(final ItemEvent e) {
                _selected_query_seq = sequencesByName.get(e.getItem());
                _mainpanel.getCurrentTreePanel().repaint();
            }
        });
    }

    private void addClickToOption(final int which, final String title) {
        _click_to_combobox.addItem(title);
        _click_to_names.add(title);
        _all_click_to_names.put(Integer.valueOf(which), title);
        if (_configuration.isApplyCustomGuiColors()) {
            _click_to_combobox.setBackground(getConfiguration().getGuiButtonBackgroundColor());
            _click_to_combobox.setForeground(getConfiguration().getGuiButtonTextColor());
        }
    }

    /* GUILHEM_BEG */
    private void addSequenceRelationBlock() {
        final JLabel spacer = new JLabel("");
        spacer.setSize(1, 1);
        add(spacer);
        final JLabel mainLabel = new JLabel("Sequence relations to display");
        final JLabel typeLabel = customizeLabel(new JLabel("(type) "), getConfiguration());
        typeLabel.setFont(ControlPanel.js_font.deriveFont(7));
        getSequenceRelationTypeBox().setFocusable(false);
        _sequence_relation_type_box.setFont(ControlPanel.js_font);
        if (_configuration.isApplyCustomGuiColors()) {
            _sequence_relation_type_box.setBackground(getConfiguration().getGuiButtonBackgroundColor());
            _sequence_relation_type_box.setForeground(getConfiguration().getGuiButtonTextColor());
        }
        _sequence_relation_type_box.setRenderer(new ListCellRenderer<Object>() {

            @Override
            public Component getListCellRendererComponent(final JList<?> list,
                                                          final Object value,
                                                          final int index,
                                                          final boolean isSelected,
                                                          final boolean cellHasFocus) {
                final Component component = new DefaultListCellRenderer()
                        .getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                if ((value != null) && (value instanceof SequenceRelation.SEQUENCE_RELATION_TYPE)) {
                    ((DefaultListCellRenderer) component).setText(SequenceRelation
                            .getPrintableNameByType((SequenceRelation.SEQUENCE_RELATION_TYPE) value));
                }
                return component;
            }
        });
        final GridBagLayout gbl = new GridBagLayout();
        _sequence_relation_type_box.setMinimumSize(new Dimension(115, 17));
        _sequence_relation_type_box.setPreferredSize(new Dimension(115, 20));
        final JPanel horizGrid = new JPanel(gbl);
        horizGrid.setBackground(getBackground());
        horizGrid.add(typeLabel);
        horizGrid.add(_sequence_relation_type_box);
        add(customizeLabel(mainLabel, getConfiguration()));
        add(horizGrid);
        add(getSequenceRelationBox());
        if (_configuration.doDisplayOption(Configuration.show_relation_confidence)) {
            addCheckbox(Configuration.show_relation_confidence,
                    _configuration.getDisplayTitle(Configuration.show_relation_confidence));
            setCheckbox(Configuration.show_relation_confidence,
                    _configuration.doCheckOption(Configuration.show_relation_confidence));
        }
    }// addSequenceRelationBlock

    /* GUILHEM_END */
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
        if ((field == null) || (mode == null)) {
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

    /** The complement of {@code matched} over the nodes that CARRY {@code field} (the "Inverse" modifier): for a
     *  numeric field the nodes that actually have a value for it (branch length / support / a numeric property live
     *  outside node data), for a text field any data-bearing node (the legacy behaviour). */
    private static Set<Long> complement(final Phylogeny tree, final SearchField field, final Set<Long> matched) {
        final Set<Long> out = new HashSet<>();
        for (final PhylogenyNode n : PhylogenyMethods.obtainAllNodesAsList(tree)) {
            if (matched.contains(n.getId())) {
                continue;
            }
            final boolean carries = field.isNumeric() ? (field.numericValues(n).length > 0) : n.isHasNodeData();
            if (carries) {
                out.add(n.getId());
            }
        }
        return out;
    }

    private void setTreeDisplayType(final int index, final Options.PHYLOGENY_DISPLAY_TYPE t) {
        getTreeDisplayTypes().set(index, t);
    }

    private void setupClickToOptions() {
        final int default_option = _configuration.getDefaultDisplayClicktoOption();
        int selected_index = 0;
        int cb_index = 0;
        if (_configuration.doDisplayClickToOption(Configuration.display_node_data)) {
            _show_data_item = cb_index;
            addClickToOption(Configuration.display_node_data,
                    _configuration.getClickToTitle(Configuration.display_node_data));
            if (default_option == Configuration.display_node_data) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.collapse_uncollapse)) {
            _collapse_cb_item = cb_index;
            addClickToOption(Configuration.collapse_uncollapse,
                    _configuration.getClickToTitle(Configuration.collapse_uncollapse));
            if (default_option == Configuration.collapse_uncollapse) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.uncollapse_all)) {
            _uncollapse_all_cb_item = cb_index;
            addClickToOption(Configuration.uncollapse_all,
                    _configuration.getClickToTitle(Configuration.uncollapse_all));
            if (default_option == Configuration.uncollapse_all) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.reroot)) {
            _reroot_cb_item = cb_index;
            addClickToOption(Configuration.reroot, _configuration.getClickToTitle(Configuration.reroot));
            if (default_option == Configuration.reroot) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.subtree)) {
            _subtree_cb_item = cb_index;
            addClickToOption(Configuration.subtree, _configuration.getClickToTitle(Configuration.subtree));
            if (default_option == Configuration.subtree) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.swap)) {
            _swap_cb_item = cb_index;
            addClickToOption(Configuration.swap, _configuration.getClickToTitle(Configuration.swap));
            if (default_option == Configuration.swap) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.order_subtree)) {
            _order_subtree_cb_item = cb_index;
            addClickToOption(Configuration.order_subtree,
                    _configuration.getClickToTitle(Configuration.order_subtree));
            if (default_option == Configuration.order_subtree) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.sort_descendents)) {
            _sort_descendents_item = cb_index;
            addClickToOption(Configuration.sort_descendents,
                    _configuration.getClickToTitle(Configuration.sort_descendents));
            if (default_option == Configuration.sort_descendents) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.color_node_font)) {
            _color_node_font_item = cb_index;
            addClickToOption(Configuration.color_node_font,
                    _configuration.getClickToTitle(Configuration.color_node_font));
            if (default_option == Configuration.color_node_font) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.color_subtree)) {
            _color_subtree_cb_item = cb_index;
            addClickToOption(Configuration.color_subtree,
                    _configuration.getClickToTitle(Configuration.color_subtree));
            if (default_option == Configuration.color_subtree) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.open_seq_web)) {
            _open_seq_web_item = cb_index;
            addClickToOption(Configuration.open_seq_web,
                    _configuration.getClickToTitle(Configuration.open_seq_web));
            if (default_option == Configuration.open_seq_web) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.open_pdb_web)) {
            _open_pdb_item = cb_index;
            addClickToOption(Configuration.open_pdb_web,
                    _configuration.getClickToTitle(Configuration.open_pdb_web));
            if (default_option == Configuration.open_pdb_web) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.open_tax_web)) {
            _open_tax_web_item = cb_index;
            addClickToOption(Configuration.open_tax_web,
                    _configuration.getClickToTitle(Configuration.open_tax_web));
            if (default_option == Configuration.open_tax_web) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.blast)) {
            _blast_item = cb_index;
            addClickToOption(Configuration.blast, _configuration.getClickToTitle(Configuration.blast));
            if (default_option == Configuration.blast) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (_configuration.doDisplayClickToOption(Configuration.select_nodes)) {
            _select_nodes_item = cb_index;
            addClickToOption(Configuration.select_nodes,
                    _configuration.getClickToTitle(Configuration.select_nodes));
            if (default_option == Configuration.select_nodes) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if (getOptions().isEditable()) {
            if (_configuration.doDisplayClickToOption(Configuration.cut_subtree)) {
                _cut_subtree_item = cb_index;
                addClickToOption(Configuration.cut_subtree,
                        _configuration.getClickToTitle(Configuration.cut_subtree));
                if (default_option == Configuration.cut_subtree) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if (_configuration.doDisplayClickToOption(Configuration.copy_subtree)) {
                _copy_subtree_item = cb_index;
                addClickToOption(Configuration.copy_subtree,
                        _configuration.getClickToTitle(Configuration.copy_subtree));
                if (default_option == Configuration.copy_subtree) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if (_configuration.doDisplayClickToOption(Configuration.paste_subtree)) {
                _paste_subtree_item = cb_index;
                addClickToOption(Configuration.paste_subtree,
                        _configuration.getClickToTitle(Configuration.paste_subtree));
                if (default_option == Configuration.paste_subtree) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if (_configuration.doDisplayClickToOption(Configuration.delete_subtree_or_node)) {
                _delete_node_or_subtree_item = cb_index;
                addClickToOption(Configuration.delete_subtree_or_node,
                        _configuration.getClickToTitle(Configuration.delete_subtree_or_node));
                if (default_option == Configuration.delete_subtree_or_node) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if (_configuration.doDisplayClickToOption(Configuration.add_new_node)) {
                _add_new_node_item = cb_index;
                addClickToOption(Configuration.add_new_node,
                        _configuration.getClickToTitle(Configuration.add_new_node));
                if (default_option == Configuration.add_new_node) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if (_configuration.doDisplayClickToOption(Configuration.edit_node_data)) {
                _edit_node_data_item = cb_index;
                addClickToOption(Configuration.edit_node_data,
                        _configuration.getClickToTitle(Configuration.edit_node_data));
                if (default_option == Configuration.edit_node_data) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
        }
        // Set default selection and its action
        _click_to_combobox.setSelectedIndex(selected_index);
        setClickToAction(selected_index);
    }

    private void setupDisplayCheckboxes() {
        addDisplayCheckbox(Configuration.dynamically_hide_data);
        addDisplayCheckbox(Configuration.node_data_popup);
        addDisplayCheckbox(Configuration.display_internal_data);
        addDisplayCheckbox(Configuration.display_external_data);
        addDisplayCheckbox(Configuration.use_style);
        addDisplayCheckbox(Configuration.width_branches);
        final JLabel label = new JLabel("Display Data:");
        label.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        add(label);
        addDisplayCheckbox(Configuration.show_node_names);
        addDisplayCheckbox(Configuration.shorten_labels);
        addDisplayCheckbox(Configuration.show_tax_code);
        addDisplayCheckbox(Configuration.show_taxonomy_scientific_names);
        addDisplayCheckbox(Configuration.show_taxonomy_common_names);
        addDisplayCheckbox(Configuration.show_tax_rank);
        addDisplayCheckbox(Configuration.show_seq_names);
        addDisplayCheckbox(Configuration.show_gene_names);
        addDisplayCheckbox(Configuration.show_seq_symbols);
        addDisplayCheckbox(Configuration.show_sequence_acc);
        addDisplayCheckbox(Configuration.write_confidence_values);
        addDisplayCheckbox(Configuration.write_branch_length_values);
        addDisplayCheckbox(Configuration.show_binary_characters);
        addDisplayCheckbox(Configuration.show_binary_character_counts);
        addDisplayCheckbox(Configuration.show_domain_architectures);
        // NB: "Multiple Seq Alignment" (show_mol_seqs) is deliberately NOT built — there is currently no
        // MSA viewer, so the checkbox would only mislead. Re-add here if/when one exists.
        addDisplayCheckbox(Configuration.write_events);
        addDisplayCheckbox(Configuration.show_vector_data);
        addDisplayCheckbox(Configuration.show_properties);
        // Data-dependent checkboxes start hidden; the first scan of the loaded tree reveals only the
        // ones whose data is actually present (Node Name is intentionally always shown).
        for (final int which : DATA_GATED_CHECKBOXES) {
            final JPanel p = _checkbox_panels.get(which);
            if (p != null) {
                p.setVisible(false);
            }
        }
    }

    // Build a "Display Data" checkbox and set its initial checked state. Whether it is actually
    // visible is then governed by updateDataCheckboxVisibility (data presence in the current tree).
    private void addDisplayCheckbox(final int which) {
        addCheckbox(which, _configuration.getDisplayTitle(which));
        setCheckbox(which, _configuration.doCheckOption(which));
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
                    if (mf._right_line_up_domains_cbmi != null) {
                        mf._right_line_up_domains_cbmi.setVisible(true);
                    }
                    if (mf._show_domain_labels != null) {
                        mf._show_domain_labels.setVisible(true);
                    }
                } else {
                    _domain_display_label.setVisible(false);
                    _zoom_in_domain_structure.setVisible(false);
                    _zoom_out_domain_structure.setVisible(false);
                    _decr_domain_structure_evalue_thr.setVisible(false);
                    _incr_domain_structure_evalue_thr.setVisible(false);
                    _domain_structure_evalue_thr_tf.setVisible(false);
                    if (mf._right_line_up_domains_cbmi != null) {
                        mf._right_line_up_domains_cbmi.setVisible(false);
                    }
                    if (mf._show_domain_labels != null) {
                        mf._show_domain_labels.setVisible(false);
                    }
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
        if (getConfiguration().doDisplayOption(Configuration.show_domain_architectures)) {
            nextRowGap(SECTION_GAP);
            setUpControlsForDomainStrucures();
        }
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

    void addCheckbox(final int which, final String title) {
        nextRowGap(CHECKBOX_GAP); // pack the display checkboxes tightly together
        final JPanel ch_panel = new JPanel(new BorderLayout(0, 0));
        switch (which) {
            case Configuration.display_internal_data:
                _display_internal_data = new JCheckBox(title);
                _display_internal_data.setToolTipText("To allow or disallow display of internal labels");
                addJCheckBox(_display_internal_data, ch_panel);
                add(ch_panel);
                break;
            case Configuration.display_external_data:
                _display_external_data = new JCheckBox(title);
                _display_external_data.setToolTipText("To allow or disallow display of external labels");
                addJCheckBox(_display_external_data, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_node_names:
                _show_node_names = new JCheckBox(title);
                addJCheckBox(_show_node_names, ch_panel);
                add(ch_panel);
                break;
            case Configuration.shorten_labels:
                _shorten_labels_cb = new JCheckBox(title);
                _shorten_labels_cb.setToolTipText(
                        "Display over-long external labels (e.g. full UniProt/NCBI descriptions) shortened with an ellipsis; the underlying names are not changed");
                addJCheckBox(_shorten_labels_cb, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_taxonomy_scientific_names:
                _show_taxo_scientific_names = new JCheckBox(title);
                addJCheckBox(_show_taxo_scientific_names, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_taxonomy_common_names:
                _show_taxo_common_names = new JCheckBox(title);
                addJCheckBox(_show_taxo_common_names, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_tax_code:
                _show_taxo_code = new JCheckBox(title);
                addJCheckBox(_show_taxo_code, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_tax_rank:
                _show_taxo_rank = new JCheckBox(title);
                addJCheckBox(_show_taxo_rank, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_binary_characters:
                _show_binary_characters = new JCheckBox(title);
                addJCheckBox(_show_binary_characters, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_binary_character_counts:
                _show_binary_character_counts = new JCheckBox(title);
                addJCheckBox(_show_binary_character_counts, ch_panel);
                add(ch_panel);
                break;
            case Configuration.write_confidence_values:
                _write_confidence = new JCheckBox(title);
                addJCheckBox(getWriteConfidenceCb(), ch_panel);
                add(ch_panel);
                break;
            case Configuration.write_events:
                _show_events = new JCheckBox(title);
                addJCheckBox(getShowEventsCb(), ch_panel);
                add(ch_panel);
                break;
            case Configuration.use_style:
                _use_visual_styles_cb = new JCheckBox(title);
                getUseVisualStylesCb()
                        .setToolTipText("To use visual styles (node colors, fonts) and branch colors, if present");
                addJCheckBox(getUseVisualStylesCb(), ch_panel);
                add(ch_panel);
                break;
            case Configuration.width_branches:
                _width_branches = new JCheckBox(title);
                _width_branches.setToolTipText("To use branch width values, if present");
                addJCheckBox(_width_branches, ch_panel);
                add(ch_panel);
                break;
            case Configuration.write_branch_length_values:
                _write_branch_length_values = new JCheckBox(title);
                addJCheckBox(_write_branch_length_values, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_domain_architectures:
                _show_domain_architectures = new JCheckBox(title);
                addJCheckBox(_show_domain_architectures, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_mol_seqs:
                _show_mol_seqs = new JCheckBox(title);
                addJCheckBox(_show_mol_seqs, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_seq_names:
                _show_seq_names = new JCheckBox(title);
                addJCheckBox(_show_seq_names, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_gene_names:
                _show_gene_names = new JCheckBox(title);
                addJCheckBox(_show_gene_names, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_seq_symbols:
                _show_seq_symbols = new JCheckBox(title);
                addJCheckBox(_show_seq_symbols, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_sequence_acc:
                _show_sequence_acc = new JCheckBox(title);
                addJCheckBox(_show_sequence_acc, ch_panel);
                add(ch_panel);
                break;
            case Configuration.dynamically_hide_data:
                _dynamically_hide_data = new JCheckBox(title);
                getDynamicallyHideData().setToolTipText("To hide labels depending on expected visibility");
                addJCheckBox(getDynamicallyHideData(), ch_panel);
                add(ch_panel);
                break;
            case Configuration.node_data_popup:
                _node_desc_popup_cb = new JCheckBox(title);
                getNodeDescPopupCb().setToolTipText("To enable mouse rollover display of basic node data");
                addJCheckBox(getNodeDescPopupCb(), ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_relation_confidence:
                _seq_relation_confidence_switch = new JCheckBox(title);
                addJCheckBox(_seq_relation_confidence_switch, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_vector_data:
                _show_vector_data_cb = new JCheckBox(title);
                addJCheckBox(_show_vector_data_cb, ch_panel);
                add(ch_panel);
                break;
            case Configuration.show_properties:
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
        if (force_rescan || (phy != _data_presence_for) || (_data_presence == null)) {
            _data_presence = AptxUtil.scanForDataPresence(phy);
            _data_presence_for = phy;
        }
        boolean changed = false;
        for (final int which : DATA_GATED_CHECKBOXES) {
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
        // tailor the search field selectors to this tree (identity-guarded: a no-op when it hasn't changed, so
        // the per-search repaint path pays nothing).
        rebuildSearchFields(false);
    }

    // For tests: whether the "Display Data" checkbox row for the given option constant is showing.
    boolean isDisplayDataCheckboxVisible(final int which) {
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
        }
        finally {
            _search_controls_adjusting = false;
        }
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

    Map<Integer, String> getAllClickToItems() {
        return _all_click_to_names;
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

    boolean isShowBinaryCharacterCounts() {
        return ((_show_binary_character_counts != null) && _show_binary_character_counts.isSelected());
    }

    boolean isShowBinaryCharacters() {
        return ((_show_binary_characters != null) && _show_binary_characters.isSelected());
    }

    boolean isShowConfidenceValues() {
        return ((getWriteConfidenceCb() != null) && getWriteConfidenceCb().isSelected());
    }

    boolean isShowDomainArchitectures() {
        return ((_show_domain_architectures != null) && _show_domain_architectures.isSelected());
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
    void showDomainArchitecturesFitted() {
        if ((_show_domain_architectures != null) && _show_domain_architectures.isVisible()
                && !_show_domain_architectures.isSelected()) {
            _show_domain_architectures.setSelected(true);
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

    boolean isShowSequenceRelationConfidence() {
        return ((_seq_relation_confidence_switch != null) && (_seq_relation_confidence_switch.isSelected()));
    }

    boolean isShowSequenceRelations() {
        return ((_show_sequence_relations != null) && (_show_sequence_relations.getSelectedIndex() > 0));
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
        final MainPanel main_panel = getMainPanel();
        final Phylogeny tree = main_panel.getCurrentPhylogeny();
        if ((tree == null) || tree.isEmpty()) {
            return;
        }
        String query = getSearchTextField0().getText();
        if (query != null) {
            query = query.trim();
        }
        if (!ForesterUtil.isEmpty(query)) {
            search0(main_panel, tree, query);
        } else {
            getSearchFoundCountsLabel0().setVisible(false);
            getSearchResetButton0().setEnabled(false);
            getSearchResetButton0().setVisible(false);
            searchReset0();
        }
    }

    void search1() {
        final MainPanel main_panel = getMainPanel();
        final Phylogeny tree = main_panel.getCurrentPhylogeny();
        if ((tree == null) || tree.isEmpty()) {
            return;
        }
        String query = getSearchTextField1().getText();
        if (query != null) {
            query = query.trim();
        }
        if (!ForesterUtil.isEmpty(query)) {
            search1(main_panel, tree, query);
        } else {
            getSearchFoundCountsLabel1().setVisible(false);
            getSearchResetButton1().setEnabled(false);
            getSearchResetButton1().setVisible(false);
            searchReset1();
        }
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

    void setCheckbox(final int which, final boolean state) {
        switch (which) {
            case Configuration.display_as_phylogram:
                if (getDisplayAsUnalignedPhylogramRb() != null) {
                    getDisplayAsUnalignedPhylogramRb().setSelected(state);
                    getDisplayAsAlignedPhylogramRb().setSelected(!state);
                    getDisplayAsCladogramRb().setSelected(!state);
                }
                break;
            case Configuration.display_internal_data:
                if (_display_internal_data != null) {
                    _display_internal_data.setSelected(state);
                }
                break;
            case Configuration.display_external_data:
                if (_display_external_data != null) {
                    _display_external_data.setSelected(state);
                }
                break;
            case Configuration.show_node_names:
                if (_show_node_names != null) {
                    _show_node_names.setSelected(state);
                }
                break;
            case Configuration.shorten_labels:
                if (_shorten_labels_cb != null) {
                    _shorten_labels_cb.setSelected(state);
                }
                break;
            case Configuration.show_taxonomy_scientific_names:
                if (_show_taxo_scientific_names != null) {
                    _show_taxo_scientific_names.setSelected(state);
                }
                break;
            case Configuration.show_taxonomy_common_names:
                if (_show_taxo_common_names != null) {
                    _show_taxo_common_names.setSelected(state);
                }
                break;
            case Configuration.show_tax_code:
                if (_show_taxo_code != null) {
                    _show_taxo_code.setSelected(state);
                }
                break;
            case Configuration.show_tax_rank:
                if (_show_taxo_rank != null) {
                    _show_taxo_rank.setSelected(state);
                }
                break;
            case Configuration.show_binary_characters:
                if (_show_binary_characters != null) {
                    _show_binary_characters.setSelected(state);
                }
                break;
            case Configuration.show_binary_character_counts:
                if (_show_binary_character_counts != null) {
                    _show_binary_character_counts.setSelected(state);
                }
                break;
            case Configuration.write_confidence_values:
                if (getWriteConfidenceCb() != null) {
                    getWriteConfidenceCb().setSelected(state);
                }
                break;
            case Configuration.write_events:
                if (getShowEventsCb() != null) {
                    getShowEventsCb().setSelected(state);
                }
                break;
            case Configuration.use_style:
                if (getUseVisualStylesCb() != null) {
                    getUseVisualStylesCb().setSelected(state);
                }
                break;
            case Configuration.width_branches:
                if (_width_branches != null) {
                    _width_branches.setSelected(state);
                }
                break;
            case Configuration.show_domain_architectures:
                if (_show_domain_architectures != null) {
                    _show_domain_architectures.setSelected(state);
                }
                break;
            case Configuration.write_branch_length_values:
                if (_write_branch_length_values != null) {
                    _write_branch_length_values.setSelected(state);
                }
                break;
            case Configuration.show_mol_seqs:
                if (_show_mol_seqs != null) {
                    _show_mol_seqs.setSelected(state);
                }
                break;
            case Configuration.show_seq_names:
                if (_show_seq_names != null) {
                    _show_seq_names.setSelected(state);
                }
                break;
            case Configuration.show_gene_names:
                if (_show_gene_names != null) {
                    _show_gene_names.setSelected(state);
                }
                break;
            case Configuration.show_seq_symbols:
                if (_show_seq_symbols != null) {
                    _show_seq_symbols.setSelected(state);
                }
                break;
            case Configuration.show_vector_data:
                if (_show_vector_data_cb != null) {
                    _show_vector_data_cb.setSelected(state);
                }
                break;
            case Configuration.show_properties:
                if (_show_properties_cb != null) {
                    _show_properties_cb.setSelected(state);
                }
                break;
            case Configuration.show_sequence_acc:
                if (_show_sequence_acc != null) {
                    _show_sequence_acc.setSelected(state);
                }
                break;
            case Configuration.dynamically_hide_data:
                if (getDynamicallyHideData() != null) {
                    getDynamicallyHideData().setSelected(state);
                }
                break;
            case Configuration.node_data_popup:
                if (getNodeDescPopupCb() != null) {
                    getNodeDescPopupCb().setSelected(state);
                }
                break;
            /* GUILHEM_BEG */
            case Configuration.show_relation_confidence:
                if (_seq_relation_confidence_switch != null) {
                    _seq_relation_confidence_switch.setSelected(state);
                }
                break;
            /* GUILHEM_END */
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
        } else if (action == _sort_descendents_item) {
            setActionWhenNodeClicked(NodeClickAction.SORT_DESCENDENTS);
        } else if (action == _blast_item) {
            setActionWhenNodeClicked(NodeClickAction.BLAST);
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
        } else if (action == _open_pdb_item) {
            setActionWhenNodeClicked(NodeClickAction.OPEN_PDB_WEB);
        } else if (action == _color_node_font_item) {
            setActionWhenNodeClicked(NodeClickAction.COLOR_NODE_FONT);
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

    void setupControls() {
        setupThemeButtons();
        nextRowGap(SECTION_GAP); // more space between the Light/Dark row and the P/A/C row
        setupTreeDisplayTypeOptions();
        nextRowGap(SECTION_GAP); // more space between the P/A/C row and "Color by"
        setupColorByProperty();
        setupSizeByProperty();
        setupAncestralPieProperty();
        setupDisplayCheckboxes();
        /* GUILHEM_BEG */
        // The sequence relation query selection combo-box
        if (_configuration.displaySequenceRelations()) {
            addSequenceRelationBlock();
        }
        /* GUILHEM_END */
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
        displayedPhylogenyMightHaveChanged(true);
    }

    /** Reacts to a mode-selector change: remember it (per box, per kind), show/hide the range field, re-run. */
    private void onSearchModeChanged(final boolean box_a) {
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
        displayedPhylogenyMightHaveChanged(true);
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
                if (box_a) {
                    search0();
                }
                else {
                    search1();
                }
                displayedPhylogenyMightHaveChanged(true);
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
        search0();
        search1();
        displayedPhylogenyMightHaveChanged(true);
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
                search0();
                displayedPhylogenyMightHaveChanged(true);
            }
        };
        final ActionListener action_listener = new ActionListener() {

            @Override
            public void actionPerformed(final ActionEvent e) {
                searchReset0();
                setSearchFoundCountsOnLabel0(0);
                getSearchFoundCountsLabel0().setVisible(false);
                getSearchTextField0().setText("");
                _search_range_tf_0.setText("");
                getSearchResetButton0().setEnabled(false);
                getSearchResetButton0().setVisible(false);
                displayedPhylogenyMightHaveChanged(true);
            }
        };
        _search_reset_button_0.addActionListener(action_listener);
        _search_tf_0.addKeyListener(key_adapter);
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
                search1();
                displayedPhylogenyMightHaveChanged(true);
            }
        };
        final ActionListener action_listener = new ActionListener() {

            @Override
            public void actionPerformed(final ActionEvent e) {
                searchReset1();
                setSearchFoundCountsOnLabel1(0);
                getSearchFoundCountsLabel1().setVisible(false);
                getSearchTextField1().setText("");
                _search_range_tf_1.setText("");
                getSearchResetButton1().setEnabled(false);
                getSearchResetButton1().setVisible(false);
                displayedPhylogenyMightHaveChanged(true);
            }
        };
        _search_reset_button_1.addActionListener(action_listener);
        _search_tf_1.addKeyListener(key_adapter);
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
    }

    void setVisibilityOfDomainStrucureCB() {
        try {
            if ((getCurrentTreePanel() != null) && ((getCurrentTreePanel()
                    .getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                    || (getCurrentTreePanel().getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED))) {
                if (getMainPanel().getMainFrame()._right_line_up_domains_cbmi != null) {
                    getMainPanel().getMainFrame()._right_line_up_domains_cbmi.setVisible(false);
                }
                if (getMainPanel().getMainFrame()._show_domain_labels != null) {
                    getMainPanel().getMainFrame()._show_domain_labels.setVisible(false);
                }
            } else if (isShowDomainArchitectures()) {
                if (getMainPanel().getMainFrame()._right_line_up_domains_cbmi != null) {
                    getMainPanel().getMainFrame()._right_line_up_domains_cbmi.setVisible(true);
                }
                if (getMainPanel().getMainFrame()._show_domain_labels != null) {
                    getMainPanel().getMainFrame()._show_domain_labels.setVisible(true);
                }
            } else {
                if (getMainPanel().getMainFrame()._right_line_up_domains_cbmi != null) {
                    getMainPanel().getMainFrame()._right_line_up_domains_cbmi.setVisible(false);
                }
                if (getMainPanel().getMainFrame()._show_domain_labels != null) {
                    getMainPanel().getMainFrame()._show_domain_labels.setVisible(false);
                }
            }
        } catch (final Exception ignore) {
            //not important...
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
        // Correlates option names to titles
        _all_click_to_names = new HashMap<Integer, String>();
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
            updateZoomButtonsForLayout(); // relabel the zoom cluster for the current (possibly persisted) layout
            getMainPanel().getCurrentTreePanel().updateSubSuperTreeButton();
            getMainPanel().getCurrentTreePanel().updateButtonToUncollapseAll();
            getMainPanel().getControlPanel().search0();
            getMainPanel().getControlPanel().search1();
            getMainPanel().getControlPanel().updateDomainStructureEvaluethresholdDisplay();
            updateDataCheckboxVisibility(true);
            populateColorByPropertyBox();
            populateSizeByPropertyBox();
            populateAncestralPieBox();
            if (getMainPanel().getMainFrame() != null) {
                getMainPanel().getMainFrame().updateEditMenu(); // undo history is per-tab
            }
            getSequenceRelationTypeBox().removeAllItems();
            for (final SequenceRelation.SEQUENCE_RELATION_TYPE type : getMainPanel().getCurrentPhylogeny()
                    .getRelevantSequenceRelationTypes()) {
                _sequence_relation_type_box.addItem(type);
            }
            getMainPanel().getCurrentTreePanel().repaint();
            //setSequenceRelationQueries( getMainPanel().getCurrentPhylogeny().getSequenceRelationQueries() );
            // according to GUILHEM the line above can be removed.
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
