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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.EnumSet;
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
        GET_EXT_DESC_DATA,
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
        CHANGE_NODE_FONT,
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
    private static final String COLOR_BY_PROPERTY_NONE = "None";
    private boolean _color_branches;
    private JCheckBox _use_visual_styles_cb;
    private int _color_subtree_cb_item;
    private int _change_node_font_item;
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
    private int _get_ext_desc_data;
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
                if (e.getSource() == getDisplayAsUnalignedPhylogramRb()) {
                    setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM);
                    showWhole();
                }
                if (e.getSource() == getDisplayAsAlignedPhylogramRb()) {
                    setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM);
                    showWhole();
                }
                if (e.getSource() == getDisplayAsCladogramRb()) {
                    setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM);
                    showWhole();
                }
                // Zoom buttons
                else if (e.getSource() == _zoom_in_x) {
                    zoomInX(AptxConstants.BUTTON_ZOOM_IN_FACTOR, AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR);
                    displayedPhylogenyMightHaveChanged(false);
                } else if (e.getSource() == _zoom_in_y) {
                    zoomInY(AptxConstants.BUTTON_ZOOM_IN_FACTOR);
                    displayedPhylogenyMightHaveChanged(false);
                } else if (e.getSource() == _zoom_out_x) {
                    zoomOutX(AptxConstants.BUTTON_ZOOM_OUT_FACTOR, AptxConstants.BUTTON_ZOOM_OUT_X_CORRECTION_FACTOR);
                    displayedPhylogenyMightHaveChanged(false);
                } else if (e.getSource() == _zoom_out_y) {
                    zoomOutY(AptxConstants.BUTTON_ZOOM_OUT_FACTOR);
                    displayedPhylogenyMightHaveChanged(false);
                } else if (e.getSource() == _show_whole) {
                    displayedPhylogenyMightHaveChanged(true);
                    showWhole();
                } else if (e.getSource() == _expand_y) {
                    expandYToFitLabels();
                } else if (e.getSource() == _fit_width) {
                    fitWidth();
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
                    // A Display Data checkbox (or similar) changed. If the "Visible" search option is on,
                    // the set of searched fields changed too, so re-run both searches to reflect it.
                    if ((_search_visible_cb != null) && _search_visible_cb.isSelected()) {
                        search0();
                        search1();
                    }
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
        if (getCurrentTreePanel().isCurrentTreeIsSubtree()) {
            getCurrentTreePanel().superTree();
            showWhole();
        }
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

    private void search0(final MainPanel main_panel, final Phylogeny tree, String query_str) {
        getSearchFoundCountsLabel0().setVisible(true);
        getSearchResetButton0().setEnabled(true);
        getSearchResetButton0().setVisible(true);
        String[] queries = null;
        Set<Long> nodes = null;
        // "Visible" search option: restrict the search to the currently displayed node-data fields.
        final Set<PhylogenyMethods.NDF> restrict = visibleSearchFields();
        final boolean search_props = (restrict != null) ? restrict.contains(PhylogenyMethods.NDF.Properties)
                : getOptions().isSearchProperties();
        query_str = query_str.replaceAll("\\s+", " ");
        if ((query_str.indexOf(',') >= 0) && !getOptions().isSearchWithRegex()) {
            queries = query_str.split(",+");
        } else {
            queries = new String[1];
            queries[0] = query_str.trim();
        }
        if ((queries != null) && (queries.length > 0)) {
            nodes = new HashSet<Long>();
            for (String query : queries) {
                if (ForesterUtil.isEmpty(query)) {
                    continue;
                }
                query = query.trim();
                final TreePanel tp = getMainPanel().getCurrentTreePanel();
                if ((query.indexOf('+') > 0) && !getOptions().isSearchWithRegex()) {
                    nodes.addAll(PhylogenyMethods.searchDataLogicalAnd(query.split("\\++"),
                            tree,
                            search_props,
                            getOptions().isSearchCaseSensitive(),
                            !getOptions().isMatchWholeTermsOnly(),
                            isShowDomainArchitectures(),
                            tp != null
                                    ? Math.pow(10,
                                    tp.getDomainStructureEvalueThresholdExp())
                                    : 0,
                            restrict));
                } else {
                    nodes.addAll(PhylogenyMethods
                            .searchData(query,
                                    tree,
                                    search_props,
                                    getOptions().isSearchCaseSensitive(),
                                    !getOptions().isMatchWholeTermsOnly(),
                                    getOptions().isSearchWithRegex(),
                                    isShowDomainArchitectures(),
                                    tp != null ? Math.pow(10, tp.getDomainStructureEvalueThresholdExp()) : 0,
                                    restrict));
                }
            }
            if (getOptions().isInverseSearchResult()) {
                final List<PhylogenyNode> all = PhylogenyMethods.obtainAllNodesAsList(tree);
                final Set<Long> temp_nodes = nodes;
                nodes = new HashSet<Long>();
                for (final PhylogenyNode n : all) {
                    if ((!temp_nodes.contains(n.getId())) && n.isHasNodeData()) {
                        nodes.add(n.getId());
                    }
                }
            }
        }
        if (restrict != null) { // "Visible": keep only nodes whose data is currently displayed
            nodes = retainDisplayedNodes(nodes, tree);
        }
        if ((nodes != null) && (nodes.size() > 0)) {
            main_panel.getCurrentTreePanel().setFoundNodes0(new HashSet<Long>());
            for (final Long node : nodes) {
                main_panel.getCurrentTreePanel().getFoundNodes0().add(node);
            }
            setSearchFoundCountsOnLabel0(nodes.size());
        } else {
            setSearchFoundCountsOnLabel0(0);
            searchReset0();
        }
    }

    private void search1(final MainPanel main_panel, final Phylogeny tree, String query_str) {
        getSearchFoundCountsLabel1().setVisible(true);
        getSearchResetButton1().setEnabled(true);
        getSearchResetButton1().setVisible(true);
        String[] queries = null;
        Set<Long> nodes = null;
        // "Visible" search option: restrict the search to the currently displayed node-data fields.
        final Set<PhylogenyMethods.NDF> restrict = visibleSearchFields();
        final boolean search_props = (restrict != null) ? restrict.contains(PhylogenyMethods.NDF.Properties)
                : getOptions().isSearchProperties();
        query_str = query_str.replaceAll("\\s+", " ");
        if ((query_str.indexOf(',') >= 0) && !getOptions().isSearchWithRegex()) {
            queries = query_str.split(",+");
        } else {
            queries = new String[1];
            queries[0] = query_str.trim();
        }
        if ((queries != null) && (queries.length > 0)) {
            nodes = new HashSet<Long>();
            for (String query : queries) {
                if (ForesterUtil.isEmpty(query)) {
                    continue;
                }
                query = query.trim();
                final TreePanel tp = getMainPanel().getCurrentTreePanel();
                if ((query.indexOf('+') > 0) && !getOptions().isSearchWithRegex()) {
                    nodes.addAll(PhylogenyMethods.searchDataLogicalAnd(query.split("\\++"),
                            tree,
                            search_props,
                            getOptions().isSearchCaseSensitive(),
                            !getOptions().isMatchWholeTermsOnly(),
                            isShowDomainArchitectures(),
                            tp != null
                                    ? Math.pow(10,
                                    tp.getDomainStructureEvalueThresholdExp())
                                    : 0,
                            restrict));
                } else {
                    nodes.addAll(PhylogenyMethods
                            .searchData(query,
                                    tree,
                                    search_props,
                                    getOptions().isSearchCaseSensitive(),
                                    !getOptions().isMatchWholeTermsOnly(),
                                    getOptions().isSearchWithRegex(),
                                    isShowDomainArchitectures(),
                                    tp != null ? Math.pow(10, tp.getDomainStructureEvalueThresholdExp()) : 0,
                                    restrict));
                }
            }
            if (getOptions().isInverseSearchResult()) {
                final List<PhylogenyNode> all = PhylogenyMethods.obtainAllNodesAsList(tree);
                final Set<Long> temp_nodes = nodes;
                nodes = new HashSet<Long>();
                for (final PhylogenyNode n : all) {
                    if ((!temp_nodes.contains(n.getId())) && n.isHasNodeData()) {
                        nodes.add(n.getId());
                    }
                }
            }
        }
        if (restrict != null) { // "Visible": keep only nodes whose data is currently displayed
            nodes = retainDisplayedNodes(nodes, tree);
        }
        if ((nodes != null) && (nodes.size() > 0)) {
            main_panel.getCurrentTreePanel().setFoundNodes1(new HashSet<Long>());
            for (final Long node : nodes) {
                main_panel.getCurrentTreePanel().getFoundNodes1().add(node);
            }
            setSearchFoundCountsOnLabel1(nodes.size());
        } else {
            setSearchFoundCountsOnLabel1(0);
            searchReset1();
        }
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
        if (_configuration.doDisplayClickToOption(Configuration.change_node_font)) {
            _change_node_font_item = cb_index;
            addClickToOption(Configuration.change_node_font,
                    _configuration.getClickToTitle(Configuration.change_node_font));
            if (default_option == Configuration.change_node_font) {
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
        if (_configuration.doDisplayClickToOption(Configuration.get_ext_desc_data)) {
            _get_ext_desc_data = cb_index;
            if (!ForesterUtil.isEmpty(getConfiguration().getLabelForGetExtDescendentsData())) {
                addClickToOption(Configuration.get_ext_desc_data,
                        getConfiguration().getLabelForGetExtDescendentsData());
            } else {
                addClickToOption(Configuration.get_ext_desc_data,
                        getConfiguration().getClickToTitle(Configuration.get_ext_desc_data));
            }
            if (default_option == Configuration.get_ext_desc_data) {
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
        if (_configuration.doDisplayOption(Configuration.dynamically_hide_data)) {
            addCheckbox(Configuration.dynamically_hide_data,
                    _configuration.getDisplayTitle(Configuration.dynamically_hide_data));
            setCheckbox(Configuration.dynamically_hide_data,
                    _configuration.doCheckOption(Configuration.dynamically_hide_data));
        }
        if (_configuration.doDisplayOption(Configuration.node_data_popup)) {
            addCheckbox(Configuration.node_data_popup,
                    _configuration.getDisplayTitle(Configuration.node_data_popup));
            setCheckbox(Configuration.node_data_popup, _configuration.doCheckOption(Configuration.node_data_popup));
        }
        if (_configuration.doDisplayOption(Configuration.display_internal_data)) {
            addCheckbox(Configuration.display_internal_data,
                    _configuration.getDisplayTitle(Configuration.display_internal_data));
            setCheckbox(Configuration.display_internal_data,
                    _configuration.doCheckOption(Configuration.display_internal_data));
        }
        if (_configuration.doDisplayOption(Configuration.display_external_data)) {
            addCheckbox(Configuration.display_external_data,
                    _configuration.getDisplayTitle(Configuration.display_external_data));
            setCheckbox(Configuration.display_external_data,
                    _configuration.doCheckOption(Configuration.display_external_data));
        }
        if (_configuration.doDisplayOption(Configuration.use_style)) {
            addCheckbox(Configuration.use_style, _configuration.getDisplayTitle(Configuration.use_style));
            setCheckbox(Configuration.use_style, _configuration.doCheckOption(Configuration.use_style));
        }
        if (_configuration.doDisplayOption(Configuration.width_branches)) {
            addCheckbox(Configuration.width_branches, _configuration.getDisplayTitle(Configuration.width_branches));
            setCheckbox(Configuration.width_branches, _configuration.doCheckOption(Configuration.width_branches));
        }
        final JLabel label = new JLabel("Display Data:");
        label.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            label.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        add(label);
        if (_configuration.doDisplayOption(Configuration.show_node_names)) {
            addCheckbox(Configuration.show_node_names,
                    _configuration.getDisplayTitle(Configuration.show_node_names));
            setCheckbox(Configuration.show_node_names, _configuration.doCheckOption(Configuration.show_node_names));
        }
        if (_configuration.doDisplayOption(Configuration.show_tax_code)) {
            addCheckbox(Configuration.show_tax_code, _configuration.getDisplayTitle(Configuration.show_tax_code));
            setCheckbox(Configuration.show_tax_code, _configuration.doCheckOption(Configuration.show_tax_code));
        }
        if (_configuration.doDisplayOption(Configuration.show_taxonomy_scientific_names)) {
            addCheckbox(Configuration.show_taxonomy_scientific_names,
                    _configuration.getDisplayTitle(Configuration.show_taxonomy_scientific_names));
            setCheckbox(Configuration.show_taxonomy_scientific_names,
                    _configuration.doCheckOption(Configuration.show_taxonomy_scientific_names));
        }
        if (_configuration.doDisplayOption(Configuration.show_taxonomy_common_names)) {
            addCheckbox(Configuration.show_taxonomy_common_names,
                    _configuration.getDisplayTitle(Configuration.show_taxonomy_common_names));
            setCheckbox(Configuration.show_taxonomy_common_names,
                    _configuration.doCheckOption(Configuration.show_taxonomy_common_names));
        }
        if (_configuration.doDisplayOption(Configuration.show_tax_rank)) {
            addCheckbox(Configuration.show_tax_rank, _configuration.getDisplayTitle(Configuration.show_tax_rank));
            setCheckbox(Configuration.show_tax_rank, _configuration.doCheckOption(Configuration.show_tax_rank));
        }
        if (_configuration.doDisplayOption(Configuration.show_seq_names)) {
            addCheckbox(Configuration.show_seq_names, _configuration.getDisplayTitle(Configuration.show_seq_names));
            setCheckbox(Configuration.show_seq_names, _configuration.doCheckOption(Configuration.show_seq_names));
        }
        if (_configuration.doDisplayOption(Configuration.show_gene_names)) {
            addCheckbox(Configuration.show_gene_names,
                    _configuration.getDisplayTitle(Configuration.show_gene_names));
            setCheckbox(Configuration.show_gene_names, _configuration.doCheckOption(Configuration.show_gene_names));
        }
        if (_configuration.doDisplayOption(Configuration.show_seq_symbols)) {
            addCheckbox(Configuration.show_seq_symbols,
                    _configuration.getDisplayTitle(Configuration.show_seq_symbols));
            setCheckbox(Configuration.show_seq_symbols,
                    _configuration.doCheckOption(Configuration.show_seq_symbols));
        }
        if (_configuration.doDisplayOption(Configuration.show_sequence_acc)) {
            addCheckbox(Configuration.show_sequence_acc,
                    _configuration.getDisplayTitle(Configuration.show_sequence_acc));
            setCheckbox(Configuration.show_sequence_acc,
                    _configuration.doCheckOption(Configuration.show_sequence_acc));
        }
        if (_configuration.doDisplayOption(Configuration.write_confidence_values)) {
            addCheckbox(Configuration.write_confidence_values,
                    _configuration.getDisplayTitle(Configuration.write_confidence_values));
            setCheckbox(Configuration.write_confidence_values,
                    _configuration.doCheckOption(Configuration.write_confidence_values));
        }
        if (_configuration.doDisplayOption(Configuration.write_branch_length_values)) {
            addCheckbox(Configuration.write_branch_length_values,
                    _configuration.getDisplayTitle(Configuration.write_branch_length_values));
            setCheckbox(Configuration.write_branch_length_values,
                    _configuration.doCheckOption(Configuration.write_branch_length_values));
        }
        if (_configuration.doDisplayOption(Configuration.show_binary_characters)) {
            addCheckbox(Configuration.show_binary_characters,
                    _configuration.getDisplayTitle(Configuration.show_binary_characters));
            setCheckbox(Configuration.show_binary_characters,
                    _configuration.doCheckOption(Configuration.show_binary_characters));
        }
        if (_configuration.doDisplayOption(Configuration.show_binary_character_counts)) {
            addCheckbox(Configuration.show_binary_character_counts,
                    _configuration.getDisplayTitle(Configuration.show_binary_character_counts));
            setCheckbox(Configuration.show_binary_character_counts,
                    _configuration.doCheckOption(Configuration.show_binary_character_counts));
        }
        if (_configuration.doDisplayOption(Configuration.show_domain_architectures)) {
            addCheckbox(Configuration.show_domain_architectures,
                    _configuration.getDisplayTitle(Configuration.show_domain_architectures));
            setCheckbox(Configuration.show_domain_architectures,
                    _configuration.doCheckOption(Configuration.show_domain_architectures));
        }
        if (_configuration.doDisplayOption(Configuration.show_mol_seqs)) {
            addCheckbox(Configuration.show_mol_seqs, _configuration.getDisplayTitle(Configuration.show_mol_seqs));
            setCheckbox(Configuration.show_mol_seqs, _configuration.doCheckOption(Configuration.show_mol_seqs));
        }
        if (_configuration.doDisplayOption(Configuration.write_events)) {
            addCheckbox(Configuration.write_events, _configuration.getDisplayTitle(Configuration.write_events));
            setCheckbox(Configuration.write_events, _configuration.doCheckOption(Configuration.write_events));
        }
        if (_configuration.doDisplayOption(Configuration.show_vector_data)) {
            addCheckbox(Configuration.show_vector_data,
                    _configuration.getDisplayTitle(Configuration.show_vector_data));
            setCheckbox(Configuration.show_vector_data,
                    _configuration.doCheckOption(Configuration.show_vector_data));
        }
        if (_configuration.doDisplayOption(Configuration.show_properties)) {
            addCheckbox(Configuration.show_properties,
                    _configuration.getDisplayTitle(Configuration.show_properties));
            setCheckbox(Configuration.show_properties, _configuration.doCheckOption(Configuration.show_properties));
        }
        if (_configuration.doDisplayOption(Configuration.show_taxonomy_images)) {
            addCheckbox(Configuration.show_taxonomy_images,
                    _configuration.getDisplayTitle(Configuration.show_taxonomy_images));
            setCheckbox(Configuration.show_taxonomy_images,
                    _configuration.doCheckOption(Configuration.show_taxonomy_images));
        }
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
        _return_to_whole_tree.setToolTipText("return all the way to the complete tree (if in a sub-tree) [Alt+R]");
        _return_to_whole_tree.setEnabled(false);
        _return_to_super_tree = new JButton(RETURN_UP_ONE_LEVEL_TEXT);
        _return_to_super_tree.setToolTipText("move up by one level towards the complete tree (if in a sub-tree)");
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

    void displayedPhylogenyMightHaveChanged(final boolean recalc_longest_ext_node_info) {
        if ((_mainpanel != null)
                && ((_mainpanel.getCurrentPhylogeny() != null) && !_mainpanel.getCurrentPhylogeny().isEmpty())) {
            if (recalc_longest_ext_node_info) {
                _mainpanel.getCurrentTreePanel().initNodeData();
                _mainpanel.getCurrentTreePanel().calculateLongestExtNodeInfo();
            }
            if (getOptions().isShowOverview()) {
                _mainpanel.getCurrentTreePanel().updateOvSizes();
            }
            _mainpanel.getCurrentTreePanel().recalculateMaxDistanceToRoot();
            _mainpanel.getCurrentTreePanel().rebuildPropertyColorScheme();
            setVisibilityOfDomainStrucureControls();
            updateDomainStructureEvaluethresholdDisplay();
            getMainPanel().getControlPanel();
            _mainpanel.getCurrentTreePanel().updateButtonToUncollapseAll();
            _mainpanel.getCurrentTreePanel().calculateScaleDistance();
            _mainpanel.getCurrentTreePanel().calcMaxDepth();
            _mainpanel.adjustJScrollPane();
            _mainpanel.getCurrentTreePanel().repaint();
            // _mainpanel.getCurrentTreePanel().setUpUrtFactor();
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

    boolean isAntialiasScreenText() {
        return true;
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
        } else if (action == _get_ext_desc_data) {
            setActionWhenNodeClicked(NodeClickAction.GET_EXT_DESC_DATA);
        } else if (action == _open_pdb_item) {
            setActionWhenNodeClicked(NodeClickAction.OPEN_PDB_WEB);
        } else if (action == _color_node_font_item) {
            setActionWhenNodeClicked(NodeClickAction.COLOR_NODE_FONT);
        } else if (action == _change_node_font_item) {
            setActionWhenNodeClicked(NodeClickAction.CHANGE_NODE_FONT);
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
        addControlPanelGlue();
    }

    // ---- Search options ---------------------------------------------------------------------
    // These five used to live in the Options menu, where users rarely noticed them and then
    // wondered why search behaved unexpectedly; they now sit right above the search boxes. They
    // drive the same Options state that search0()/search1() read.
    private JCheckBox _search_case_sensitive_cb;
    private JCheckBox _search_whole_words_only_cb;
    private JCheckBox _search_regex_cb;
    private JCheckBox _search_inverse_cb;
    private JCheckBox _search_properties_cb;
    private JCheckBox _search_visible_cb;

    void setupSearchOptions() {
        final JLabel header = new JLabel("Search Options:");
        header.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            header.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        add(header);
        _search_case_sensitive_cb = new JCheckBox(MainFrame.SEARCH_CASE_SENSITIVE_LABEL);
        _search_whole_words_only_cb = new JCheckBox(MainFrame.SEARCH_TERMS_ONLY_LABEL);
        _search_regex_cb = new JCheckBox(MainFrame.SEARCH_REGEX_LABEL);
        _search_inverse_cb = new JCheckBox(MainFrame.INVERSE_SEARCH_RESULT_LABEL);
        _search_properties_cb = new JCheckBox("Properties");
        _search_visible_cb = new JCheckBox("Visible");
        _search_case_sensitive_cb.setToolTipText("search is case sensitive");
        _search_whole_words_only_cb.setToolTipText("match complete words/terms only");
        _search_regex_cb.setToolTipText(MainFrame.SEARCH_WITH_REGEX_TIP);
        _search_inverse_cb.setToolTipText("select the nodes that do NOT match");
        _search_properties_cb.setToolTipText("also search in (phyloXML) properties");
        _search_visible_cb.setToolTipText("search only the data that is currently displayed (per the Display Data checkboxes)");
        final Options o = getOptions();
        if (o != null) {
            _search_case_sensitive_cb.setSelected(o.isSearchCaseSensitive());
            _search_whole_words_only_cb.setSelected(o.isMatchWholeTermsOnly());
            _search_regex_cb.setSelected(o.isSearchWithRegex());
            _search_inverse_cb.setSelected(o.isInverseSearchResult());
            _search_properties_cb.setSelected(o.isSearchProperties());
        }
        final ActionListener l = new ActionListener() {
            @Override
            public void actionPerformed(final ActionEvent e) {
                searchOptionChanged((JCheckBox) e.getSource());
            }
        };
        for (final JCheckBox cb : new JCheckBox[] { _search_case_sensitive_cb, _search_whole_words_only_cb,
                _search_regex_cb, _search_inverse_cb, _search_properties_cb, _search_visible_cb }) {
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
        add(searchOptionsRow(_search_case_sensitive_cb, _search_whole_words_only_cb));
        nextRowGap(CHECKBOX_GAP);
        add(searchOptionsRow(_search_regex_cb, _search_inverse_cb));
        nextRowGap(CHECKBOX_GAP);
        add(searchOptionsRow(_search_properties_cb, _search_visible_cb));
    }

    /**
     * When the "Visible" search option is on, the set of node-data fields that are currently displayed (per the
     * Display Data checkboxes), so the search is limited to them; {@code null} when "Visible" is off (search all
     * fields). Fields with no Display Data checkbox (taxonomy synonyms/lineage/identifier, annotations, etc.)
     * are never searched in "Visible" mode.
     */
    private Set<PhylogenyMethods.NDF> visibleSearchFields() {
        if ((_search_visible_cb == null) || !_search_visible_cb.isSelected()) {
            return null;
        }
        final EnumSet<PhylogenyMethods.NDF> fields = EnumSet.noneOf(PhylogenyMethods.NDF.class);
        if (isShowNodeNames()) {
            fields.add(PhylogenyMethods.NDF.NodeName);
        }
        if (isShowTaxonomyCode()) {
            fields.add(PhylogenyMethods.NDF.TaxonomyCode);
        }
        if (isShowTaxonomyScientificNames()) {
            fields.add(PhylogenyMethods.NDF.TaxonomyScientificName);
        }
        if (isShowTaxonomyCommonNames()) {
            fields.add(PhylogenyMethods.NDF.TaxonomyCommonName);
        }
        if (isShowSeqNames()) {
            fields.add(PhylogenyMethods.NDF.SequenceName);
        }
        if (isShowGeneNames()) {
            fields.add(PhylogenyMethods.NDF.GeneName);
        }
        if (isShowSeqSymbols()) {
            fields.add(PhylogenyMethods.NDF.SequenceSymbol);
        }
        if (isShowSequenceAcc()) {
            fields.add(PhylogenyMethods.NDF.SequenceAccession);
        }
        if (isShowDomainArchitectures()) {
            fields.add(PhylogenyMethods.NDF.Domain);
        }
        if (isShowProperties()) {
            fields.add(PhylogenyMethods.NDF.Properties);
        }
        return fields;
    }

    /**
     * For the "Visible" search option: keep only those found nodes whose data is currently displayed --
     * external (leaf) nodes only when "Show External Data" is on, internal nodes only when "Show Internal
     * Data" is on -- so the search never matches data that is not on screen.
     */
    private Set<Long> retainDisplayedNodes(final Set<Long> ids, final Phylogeny tree) {
        final boolean show_external = isShowExternalData();
        final boolean show_internal = isShowInternalData();
        if ((ids == null) || ids.isEmpty() || (show_external && show_internal)) {
            return ids;
        }
        final Set<Long> kept = new HashSet<Long>();
        for (final PhylogenyNode n : PhylogenyMethods.obtainAllNodesAsList(tree)) {
            if (ids.contains(n.getId()) && (n.isExternal() ? show_external : show_internal)) {
                kept.add(n.getId());
            }
        }
        return kept;
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
     * Pushes the search-option checkboxes into {@link Options} (which search0()/search1() read), then
     * re-runs BOTH searches and repaints so the changed options are reflected immediately. "Regex" is
     * mutually exclusive with "Match Case"/"Words", as in the old menu.
     */
    private void searchOptionChanged(final JCheckBox source) {
        if ((source == _search_regex_cb) && _search_regex_cb.isSelected()) {
            _search_case_sensitive_cb.setSelected(false);
            _search_whole_words_only_cb.setSelected(false);
        } else if (((source == _search_case_sensitive_cb) || (source == _search_whole_words_only_cb))
                && source.isSelected()) {
            _search_regex_cb.setSelected(false);
        }
        final Options o = getOptions();
        o.setSearchCaseSensitive(_search_case_sensitive_cb.isSelected());
        o.setMatchWholeTermsOnly(_search_whole_words_only_cb.isSelected());
        o.setSearchWithRegex(_search_regex_cb.isSelected());
        o.setInverseSearchResult(_search_inverse_cb.isSelected());
        o.setSearchProperties(_search_properties_cb.isSelected());
        search0();
        search1();
        // re-run alone only updates the found-node sets and count labels; repaint so the tree's
        // highlighting reflects the new options right away (the old menu got this for free when the
        // menu closed and repainted the canvas).
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
        _search_found_label_0 = new JLabel();
        getSearchFoundCountsLabel0().setVisible(false);
        _search_found_label_0.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            _search_found_label_0.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _search_tf_0 = new JTextField(3);
        _search_tf_0.setFont(ControlPanel.jcb_font);
        _search_tf_0.setToolTipText(SEARCH_TIP_TEXT);
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
        final JPanel s_panel_1 = new JPanel(new BorderLayout());
        final JPanel s_panel_2 = new JPanel(new GridLayout(1, 2, 0, 0));
        s_panel_1.setBackground(getBackground());
        add(s_panel_1);
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
        _search_found_label_1 = new JLabel();
        getSearchFoundCountsLabel1().setVisible(false);
        _search_found_label_1.setFont(ControlPanel.jcb_bold_font);
        if (getConfiguration().isApplyCustomGuiColors()) {
            _search_found_label_1.setForeground(getConfiguration().getGuiCheckboxTextColor());
        }
        _search_tf_1 = new JTextField(3);
        _search_tf_1.setFont(ControlPanel.jcb_font);
        _search_tf_1.setToolTipText(SEARCH_TIP_TEXT);
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
        final JPanel s_panel_1 = new JPanel(new BorderLayout());
        final JPanel s_panel_2 = new JPanel(new GridLayout(1, 2, 0, 0));
        s_panel_1.setBackground(getBackground());
        add(s_panel_1);
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
        _mainpanel.getCurrentTreePanel().resetPreferredSize();
        _mainpanel.adjustJScrollPane();
        _mainpanel.getCurrentTreePanel().repaint();
        _mainpanel.getCurrentTreePanel().validate();
        _mainpanel.validate();
        _mainpanel.getCurrentTreePanel().calcParametersForPainting(_mainpanel.getSizeOfViewport().width,
                _mainpanel.getSizeOfViewport().height);
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
            if (getCurrentTreePanel().isPhyHasBranchLengths()
                    && (getCurrentTreePanel().getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)) {
                setDrawPhylogramEnabled(true);
                setTreeDisplayType(getTreeDisplayType(getMainPanel().getCurrentTabIndex()));
            } else {
                setDrawPhylogramEnabled(false);
                setTreeDisplayType(Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM);
            }
            getMainPanel().getMainFrame()
                    .setSelectedTypeInTypeMenu(getMainPanel().getCurrentTreePanel().getPhylogenyGraphicsType());
            getMainPanel().getCurrentTreePanel().updateSubSuperTreeButton();
            getMainPanel().getCurrentTreePanel().updateButtonToUncollapseAll();
            getMainPanel().getControlPanel().search0();
            getMainPanel().getControlPanel().search1();
            getMainPanel().getControlPanel().updateDomainStructureEvaluethresholdDisplay();
            populateColorByPropertyBox();
            getSequenceRelationTypeBox().removeAllItems();
            for (final SequenceRelation.SEQUENCE_RELATION_TYPE type : getMainPanel().getCurrentPhylogeny()
                    .getRelevantSequenceRelationTypes()) {
                _sequence_relation_type_box.addItem(type);
            }
            getMainPanel().getCurrentTreePanel().repaint();
            //setSequenceRelationQueries( getMainPanel().getCurrentPhylogeny().getSequenceRelationQueries() );
            // according to GUILHEM the line above can be removed.
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

    final void zoomInX(final float factor, final float x_correction_factor) {
        final JScrollBar sb = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar();
        final TreePanel treepanel = getMainPanel().getCurrentTreePanel();
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
        final JScrollBar sb = getMainPanel().getCurrentScrollPane().getVerticalScrollBar();
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
        final JScrollBar sb = getMainPanel().getCurrentScrollPane().getVerticalScrollBar();
        final TreePanel treepanel = getMainPanel().getCurrentTreePanel();
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
        treepanel.multiplyUrtFactor(1f);
        if ((treepanel.getXdistance() * factor) > 0.0) {
            final JScrollBar sb = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar();
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
        treepanel.multiplyUrtFactor(0.9f);
        if ((treepanel.getYdistance() * factor) > 0.0) {
            final JScrollBar sb = getMainPanel().getCurrentScrollPane().getVerticalScrollBar();
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
