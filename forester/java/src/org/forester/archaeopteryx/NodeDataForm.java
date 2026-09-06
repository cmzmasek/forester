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
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.DefaultCellEditor;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import javax.swing.text.DefaultCaret;
import javax.swing.text.JTextComponent;

import org.forester.archaeopteryx.NodeDataDraft.ConfidenceDraft;
import org.forester.archaeopteryx.NodeDataDraft.Problem;
import org.forester.archaeopteryx.NodeDataDraft.PropertyDraft;
import org.forester.archaeopteryx.NodeDataDraft.SequenceDraft;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.NodeData;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.ForesterUtil;
import org.forester.util.TaxonomyUtil;

/**
 * The node-data window's content: one scrolling page of collapsible sections (Basic, Taxonomy, Sequences, Events,
 * Date, Distribution, Reference, Properties) built from a {@link NodeDataDraft}. In {@link Mode#EDIT} every value
 * is a live widget; the form validates continuously (a problem outlines its widget and is reported through
 * {@link #problems()}) and tracks dirtiness against the draft it opened with ({@link #isDirty()}). NOTHING reaches
 * the tree until {@link #write()}, which is atomic (no write while any field is invalid), one undo step, and appends
 * a provenance sentence. In {@link Mode#VIEW} the same page shows only what is present, read-only, plus the
 * read-only extras (annotations, lineage, binary characters, ...) that the editor does not edit.
 * <p>
 * The window chrome (buttons, status line, close confirmation) lives in {@link NodeFrame}; this panel is a plain
 * {@code JPanel} so its logic runs headless in tests.
 */
final class NodeDataForm extends JPanel {

    enum Mode {
        VIEW,
        EDIT
    }

    private static final long serialVersionUID = 1L;
    private static final int  SECTION_INDENT   = 18;
    private static final int  MOL_SEQ_ROWS     = 3;
    /** The preferred width a read-only value reports (it stretches to the row anyway). */
    private static final int  VIEW_VALUE_WIDTH = 120;

    private final PhylogenyNode                 _node;
    private final TreePanel                     _tree_panel;
    private final Mode                          _mode;
    private final boolean                       _internal;
    private NodeDataDraft                       _baseline;
    private final Map<String, JComponent>       _fields         = new HashMap<>();
    private final Map<String, Section>          _sections       = new LinkedHashMap<>();
    private final List<ConfidenceRow>           _confidence_rows = new ArrayList<>();
    private final List<SequenceCard>            _sequence_cards = new ArrayList<>();
    private final List<Runnable>                _change_listeners = new ArrayList<>();
    private final List<JComponent>              _outlined       = new ArrayList<>();
    private JPanel                              _confidence_list;
    private JPanel                              _sequence_list;
    private PropertyTableModel                  _property_model;
    private JTable                              _property_table;
    private JButton                             _remove_property_button;
    private JScrollPane                         _scroll;
    private final int                           _label_width;
    private boolean                             _building       = true;
    /** The draft/problems computed for the CURRENT widget state (null = stale); every keystroke recomputes them
     *  once in {@link #fireChanged}, and isDirty()/problems()/write() reuse them instead of re-reading the form. */
    private NodeDataDraft                       _current;
    private List<Problem>                       _current_problems;
    private final DocumentListener              _doc_listener   = onChange( this::fireChanged );

    /** A document listener that runs {@code r} on every kind of change. */
    private static DocumentListener onChange( final Runnable r ) {
        return new DocumentListener() {

            @Override
            public void insertUpdate( final DocumentEvent e ) {
                r.run();
            }

            @Override
            public void removeUpdate( final DocumentEvent e ) {
                r.run();
            }

            @Override
            public void changedUpdate( final DocumentEvent e ) {
                r.run();
            }
        };
    }

    NodeDataForm( final PhylogenyNode node, final TreePanel tree_panel, final Mode mode ) {
        super( new BorderLayout() );
        _node = node;
        _tree_panel = tree_panel;
        _mode = mode;
        _internal = !node.isExternal();
        _baseline = NodeDataDraft.from( node );
        _label_width = getFontMetrics( getFont() ).stringWidth( "Scientific name" ) + 12;
        add( buildHeader(), BorderLayout.NORTH );
        final JPanel page = new ScrollablePage();
        page.setLayout( new BoxLayout( page, BoxLayout.Y_AXIS ) );
        page.setBorder( BorderFactory.createEmptyBorder( 4, 14, 10, 14 ) );
        buildSections( page, _baseline );
        _scroll = new JScrollPane( page );
        _scroll.setBorder( BorderFactory.createEmptyBorder() );
        _scroll.setHorizontalScrollBarPolicy( JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
        _scroll.getVerticalScrollBar().setUnitIncrement( 16 );
        add( _scroll, BorderLayout.CENTER );
        _building = false;
        applyProblems( problems() );
        SwingUtilities.invokeLater( this::scrollToTop );
    }

    /** Puts the page at its top (a focused or caret-bearing field further down must not win on open). */
    private void scrollToTop() {
        _scroll.getViewport().setViewPosition( new java.awt.Point( 0, 0 ) );
    }

    // ------------------------------------------------------------------ public-ish API (used by NodeFrame + tests)
    boolean isEditable() {
        return _mode == Mode.EDIT;
    }

    /** The draft this form last read from / wrote to the node (what dirtiness is measured against). */
    NodeDataDraft baseline() {
        return _baseline;
    }

    /** Notified whenever a value changes (so dirtiness / validity may have changed). */
    void addChangeListener( final Runnable r ) {
        _change_listeners.add( r );
    }

    /** The current widget values as a draft (cached until the next change). In VIEW mode nothing can change, so
     *  this is simply the draft the node was read into (the read-only widgets are display, not a source). */
    NodeDataDraft collect() {
        if ( !isEditable() ) {
            return _baseline;
        }
        commitTableEdits();
        if ( _current == null ) {
            _current = readWidgets();
        }
        return _current;
    }

    private NodeDataDraft readWidgets() {
        final NodeDataDraft d = new NodeDataDraft();
        d.name = text( NodeDataDraft.NAME );
        d.branchLength = text( NodeDataDraft.BRANCH_LENGTH );
        d.branchWidth = text( NodeDataDraft.BRANCH_WIDTH );
        for( final ConfidenceRow r : _confidence_rows ) {
            d.confidences.add( r.toDraft() );
        }
        d.taxId = text( NodeDataDraft.TAX_ID );
        d.taxProvider = text( NodeDataDraft.TAX_PROVIDER );
        d.taxCode = text( NodeDataDraft.TAX_CODE );
        d.taxSciName = text( NodeDataDraft.TAX_SCI_NAME );
        d.taxAuthority = text( NodeDataDraft.TAX_AUTHORITY );
        d.taxCommonName = text( NodeDataDraft.TAX_COMMON_NAME );
        d.taxSynonyms = text( NodeDataDraft.TAX_SYNONYMS );
        d.taxRank = text( NodeDataDraft.TAX_RANK );
        d.taxUris = text( NodeDataDraft.TAX_URIS );
        for( final SequenceCard c : _sequence_cards ) {
            d.sequences.add( c.toDraft() );
        }
        d.duplications = text( NodeDataDraft.EV_DUPLICATIONS );
        d.speciations = text( NodeDataDraft.EV_SPECIATIONS );
        d.geneLosses = text( NodeDataDraft.EV_GENE_LOSSES );
        d.dateDesc = text( NodeDataDraft.DATE_DESC );
        d.dateValue = text( NodeDataDraft.DATE_VALUE );
        d.dateMin = text( NodeDataDraft.DATE_MIN );
        d.dateMax = text( NodeDataDraft.DATE_MAX );
        d.dateUnit = text( NodeDataDraft.DATE_UNIT );
        d.distDesc = text( NodeDataDraft.DIST_DESC );
        d.distDatum = text( NodeDataDraft.DIST_DATUM );
        d.distLat = text( NodeDataDraft.DIST_LAT );
        d.distLong = text( NodeDataDraft.DIST_LONG );
        d.distAlt = text( NodeDataDraft.DIST_ALT );
        d.distAltUnit = text( NodeDataDraft.DIST_ALT_UNIT );
        d.refDesc = text( NodeDataDraft.REF_DESC );
        d.refDoi = text( NodeDataDraft.REF_DOI );
        if ( _property_model != null ) {
            for( final PropertyDraft p : _property_model.rows() ) {
                d.properties.add( p.copy() );
            }
        }
        return d;
    }

    boolean isDirty() {
        return isEditable() && !collect().equals( _baseline );
    }

    /** Every current validation problem (empty = writable); cached with the draft. Always empty in VIEW mode. */
    List<Problem> problems() {
        if ( !isEditable() ) {
            return java.util.Collections.emptyList();
        }
        collect();
        if ( _current_problems == null ) {
            _current_problems = _current.validate( _internal );
        }
        return _current_problems;
    }

    /**
     * Writes the current values to the node, IF they validate: one undo checkpoint, only the changed fields (see
     * {@link NodeDataDraft#writeTo}), a provenance sentence on the tree, then the tree panel is refreshed. A draft
     * with no changes is a no-op that still returns true; an invalid one returns false and touches nothing.
     */
    boolean write() {
        if ( !isEditable() ) {
            return false;
        }
        final NodeDataDraft draft = collect();
        final List<Problem> ps = problems();
        applyProblems( ps );
        if ( !ps.isEmpty() ) {
            return false;
        }
        final Set<String> sections = draft.changedSections( _baseline, _internal );
        if ( sections.isEmpty() ) {
            return true;
        }
        if ( _tree_panel != null ) {
            _tree_panel.pushUndoCheckpoint( "Edit Node Data" );
        }
        draft.writeTo( _node, _baseline );
        // a new card's sequence now exists on the node -- bind the card to it so a later edit mutates that object
        for( int i = 0; i < _sequence_cards.size(); ++i ) {
            _sequence_cards.get( i ).bindOrigin( draft.sequences.get( i ).origin );
        }
        _baseline = draft.copy();
        if ( _tree_panel != null ) {
            final Phylogeny phy = _tree_panel.getPhylogeny();
            if ( phy != null ) {
                final String prov = NodeDataDraft.provenance( nodeLabel( _node ), sections );
                final String existing = phy.getDescription();
                phy.setDescription( ForesterUtil.isEmpty( existing ) ? prov : existing + " " + prov );
            }
            refreshTreePanel( sections );
        }
        fireChanged();
        return true;
    }

    /** What the window title / header calls this node. */
    static String nodeLabel( final PhylogenyNode n ) {
        if ( !ForesterUtil.isEmpty( n.getName() ) ) {
            return n.getName();
        }
        if ( n.getNodeData().isHasTaxonomy() ) {
            final Taxonomy t = n.getNodeData().getTaxonomy();
            if ( !ForesterUtil.isEmpty( t.getScientificName() ) ) {
                return t.getScientificName();
            }
            if ( !ForesterUtil.isEmpty( t.getTaxonomyCode() ) ) {
                return t.getTaxonomyCode();
            }
        }
        if ( n.getNodeData().isHasSequence() && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getName() ) ) {
            return n.getNodeData().getSequence().getName();
        }
        return "node " + n.getId();
    }

    /** The muted line under the header title: what kind of node this is and where it sits. Pure. */
    static String headerSubtitle( final PhylogenyNode n ) {
        final StringBuilder sb = new StringBuilder();
        if ( n.isRoot() ) {
            sb.append( "Root" );
        }
        else {
            sb.append( n.isExternal() ? "External node" : "Internal node" );
        }
        if ( !n.isExternal() ) {
            final int children = n.getNumberOfDescendants();
            sb.append( " · " ).append( children ).append( children == 1 ? " child" : " children" );
            // ONE pass over the tips: how many, how many distinct taxonomies, how many tips carry none
            final Set<Taxonomy> distinct = new HashSet<>();
            int tips = 0;
            int without = 0;
            for( final PhylogenyNode tip : n.getAllExternalDescendants() ) {
                ++tips;
                if ( tip.getNodeData().isHasTaxonomy() && !tip.getNodeData().getTaxonomy().isEmpty() ) {
                    distinct.add( tip.getNodeData().getTaxonomy() );
                }
                else {
                    ++without;
                }
            }
            sb.append( " · " ).append( tips ).append( tips == 1 ? " tip" : " tips" );
            if ( !distinct.isEmpty() ) {
                sb.append( " · " ).append( distinct.size() )
                        .append( distinct.size() == 1 ? " taxonomy" : " distinct taxonomies" );
                if ( without > 0 ) {
                    sb.append( " (" ).append( without ).append( without == 1 ? " tip without)" : " tips without)" );
                }
            }
        }
        if ( !n.isRoot() ) {
            sb.append( " · depth " ).append( n.calculateDepth() );
            final double d = n.calculateDistanceToRoot();
            if ( d > 0 ) {
                sb.append( " · " ).append( NodeDataDraft.formatNumber( d ) ).append( " from root" );
            }
        }
        return sb.toString();
    }

    // ------------------------------------------------------------------ test hooks
    JComponent fieldForTest( final String key ) {
        return fieldFor( key );
    }

    void setTextForTest( final String key, final String value ) {
        setValue( fieldFor( key ), value );
    }

    boolean isSectionExpandedForTest( final String title ) {
        final Section s = _sections.get( title );
        return ( s != null ) && s.isExpanded();
    }

    boolean hasSectionForTest( final String title ) {
        return _sections.containsKey( title );
    }

    void toggleSectionForTest( final String title ) {
        _sections.get( title ).toggle();
    }

    int sequenceCardCountForTest() {
        return _sequence_cards.size();
    }

    void addSequenceForTest() {
        addSequenceCard( new SequenceDraft(), true );
    }

    void removeSequenceForTest( final int i ) {
        removeSequenceCard( _sequence_cards.get( i ) );
    }

    int confidenceRowCountForTest() {
        return _confidence_rows.size();
    }

    void addConfidenceForTest() {
        addConfidenceRow( new ConfidenceDraft(), true );
    }

    void removeConfidenceForTest( final int i ) {
        removeConfidenceRow( _confidence_rows.get( i ) );
    }

    JTable propertyTableForTest() {
        return _property_table;
    }

    void addPropertyForTest( final PropertyDraft p ) {
        _property_model.add( p );
    }

    void removePropertyForTest( final int row ) {
        _property_model.remove( row );
    }

    boolean isOutlinedForTest( final String key ) {
        final JComponent c = outlineTarget( fieldFor( key ) );
        return ( c != null ) && "error".equals( c.getClientProperty( "JComponent.outline" ) );
    }

    JScrollPane scrollPaneForTest() {
        return _scroll;
    }

    // ------------------------------------------------------------------ building
    private JComponent buildHeader() {
        final JPanel header = new JPanel( new BorderLayout( 0, 2 ) );
        header.setBorder( BorderFactory.createEmptyBorder( 12, 14, 8, 14 ) );
        final JLabel title = new JLabel( nodeLabel( _node ) );
        title.putClientProperty( "FlatLaf.styleClass", "h2" );
        title.setFont( title.getFont().deriveFont( Font.BOLD, title.getFont().getSize2D() * 1.35f ) );
        final JLabel subtitle = new JLabel( headerSubtitle( _node ) );
        subtitle.setForeground( mutedColor() );
        header.add( title, BorderLayout.NORTH );
        header.add( subtitle, BorderLayout.CENTER );
        return header;
    }

    private void buildSections( final JPanel page, final NodeDataDraft d ) {
        final NodeData nd = _node.getNodeData();
        // -- Basic --
        {
            final Grid g = new Grid();
            addText( g, "Name", NodeDataDraft.NAME, d.name, null );
            addText( g, "Branch length", NodeDataDraft.BRANCH_LENGTH, d.branchLength,
                     _node.isRoot() ? "none (root)" : null, "Width", NodeDataDraft.BRANCH_WIDTH, d.branchWidth, "1" );
            if ( isEditable() || !d.confidences.isEmpty() ) {
                _confidence_list = new JPanel();
                _confidence_list.setOpaque( false );
                _confidence_list.setLayout( new BoxLayout( _confidence_list, BoxLayout.Y_AXIS ) );
                for( final ConfidenceDraft c : d.confidences ) {
                    addConfidenceRow( c, false );
                }
                final JPanel holder = new JPanel( new BorderLayout() );
                holder.setOpaque( false );
                holder.add( _confidence_list, BorderLayout.CENTER );
                if ( isEditable() ) {
                    holder.add( addButton( "+ Add confidence", () -> addConfidenceRow( new ConfidenceDraft(), true ) ),
                                BorderLayout.SOUTH );
                }
                g.row( "Confidence", holder, true );
            }
            if ( !isEditable() && _node.getBranchData().isHasBranchColor() ) {
                final Color c = _node.getBranchData().getBranchColor().getValue();
                final JPanel sw = new JPanel();
                sw.setBackground( c );
                sw.setPreferredSize( new Dimension( 14, 14 ) );
                sw.setBorder( BorderFactory.createLineBorder( borderColor() ) );
                final JPanel row = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 0 ) );
                row.setOpaque( false );
                row.add( sw );
                row.add( new JLabel( c.getRed() + ", " + c.getGreen() + ", " + c.getBlue() ) );
                g.row( "Branch color", row, false );
            }
            final boolean has = !d.name.isEmpty() || !d.branchLength.isEmpty() || !d.branchWidth.isEmpty()
                    || !d.confidences.isEmpty() || _node.getBranchData().isHasBranchColor();
            if ( isEditable() || has ) {
                addSection( page, NodeDataDraft.SEC_BASIC, null, g, true );
            }
        }
        // -- Taxonomy --
        {
            final Grid g = new Grid();
            addText( g, "Scientific name", NodeDataDraft.TAX_SCI_NAME, d.taxSciName, "e.g. Homo sapiens" );
            addText( g, "Code", NodeDataDraft.TAX_CODE, d.taxCode, "e.g. HUMAN", "Rank", NodeDataDraft.TAX_RANK,
                     d.taxRank, rankChoices(), true );
            addText( g, "Identifier", NodeDataDraft.TAX_ID, d.taxId, "e.g. 9606", "Provider",
                     NodeDataDraft.TAX_PROVIDER, d.taxProvider, "e.g. ncbi" );
            addText( g, "Common name", NodeDataDraft.TAX_COMMON_NAME, d.taxCommonName, null );
            addText( g, "Authority", NodeDataDraft.TAX_AUTHORITY, d.taxAuthority, "e.g. Linnaeus, 1758" );
            addArea( g, "Synonyms", NodeDataDraft.TAX_SYNONYMS, d.taxSynonyms, "one per line", 2, false );
            addArea( g, "URIs", NodeDataDraft.TAX_URIS, d.taxUris, "one per line, e.g. https://...", 2, false );
            if ( !isEditable() && nd.isHasTaxonomy() && !ForesterUtil.isEmpty( nd.getTaxonomy().getLineage() ) ) {
                viewRow( g, "Lineage", String.join( " > ", nonEmpty( nd.getTaxonomy().getLineage() ) ) );
            }
            if ( isEditable() || d.hasTaxonomy() ) {
                addSection( page, NodeDataDraft.SEC_TAXONOMY, null, g, d.hasTaxonomy() );
            }
        }
        // -- Sequences --
        {
            _sequence_list = new JPanel();
            _sequence_list.setOpaque( false );
            _sequence_list.setLayout( new BoxLayout( _sequence_list, BoxLayout.Y_AXIS ) );
            for( final SequenceDraft s : d.sequences ) {
                addSequenceCard( s, false );
            }
            final JPanel holder = new JPanel( new BorderLayout() );
            holder.setOpaque( false );
            holder.add( _sequence_list, BorderLayout.CENTER );
            if ( isEditable() ) {
                holder.add( addButton( "+ Add sequence", () -> addSequenceCard( new SequenceDraft(), true ) ),
                            BorderLayout.SOUTH );
            }
            if ( isEditable() || !d.sequences.isEmpty() ) {
                addSection( page, NodeDataDraft.SEC_SEQUENCES, sequenceDetail(), holder, !d.sequences.isEmpty() );
            }
        }
        // -- Events (internal nodes only) --
        if ( _internal ) {
            final Grid g = new Grid();
            addText( g, "Duplications", NodeDataDraft.EV_DUPLICATIONS, d.duplications, "0" );
            addText( g, "Speciations", NodeDataDraft.EV_SPECIATIONS, d.speciations, "0" );
            addText( g, "Gene losses", NodeDataDraft.EV_GENE_LOSSES, d.geneLosses, "0" );
            boolean has = d.hasEvents();
            if ( !isEditable() && nd.isHasEvent() ) {
                final Event e = nd.getEvent();
                if ( !e.isUnassigned() ) {
                    viewRow( g, "Type", e.getEventType().toString() );
                    has = true;
                }
                if ( e.getConfidence() != null ) {
                    viewRow( g, "Confidence", e.getConfidence().asText().toString() );
                    has = true;
                }
            }
            if ( isEditable() || has ) {
                addSection( page, NodeDataDraft.SEC_EVENTS, null, g, d.hasEvents() );
            }
        }
        // -- Date --
        {
            final Grid g = new Grid();
            addText( g, "Value", NodeDataDraft.DATE_VALUE, d.dateValue, null, "Unit", NodeDataDraft.DATE_UNIT,
                     d.dateUnit, "e.g. mya" );
            addText( g, "Min", NodeDataDraft.DATE_MIN, d.dateMin, null, "Max", NodeDataDraft.DATE_MAX, d.dateMax,
                     null );
            addText( g, "Description", NodeDataDraft.DATE_DESC, d.dateDesc, null );
            if ( isEditable() || d.hasDate() ) {
                addSection( page, NodeDataDraft.SEC_DATE, null, g, d.hasDate() );
            }
        }
        // -- Distribution --
        {
            final Grid g = new Grid();
            addText( g, "Description", NodeDataDraft.DIST_DESC, d.distDesc, "e.g. Pacific Northwest" );
            addText( g, "Latitude", NodeDataDraft.DIST_LAT, d.distLat, "-90 to 90", "Longitude",
                     NodeDataDraft.DIST_LONG, d.distLong, "-180 to 180" );
            addText( g, "Altitude", NodeDataDraft.DIST_ALT, d.distAlt, null, "Unit", NodeDataDraft.DIST_ALT_UNIT,
                     d.distAltUnit, "e.g. m" );
            addText( g, "Geodetic datum", NodeDataDraft.DIST_DATUM, d.distDatum, "e.g. WGS84" );
            if ( isEditable() || d.hasDistribution() ) {
                addSection( page, NodeDataDraft.SEC_DISTRIBUTION, null, g, d.hasDistribution() );
            }
        }
        // -- Reference --
        {
            final Grid g = new Grid();
            addText( g, "DOI", NodeDataDraft.REF_DOI, d.refDoi, "e.g. 10.1093/bioinformatics/btq243" );
            addText( g, "Description", NodeDataDraft.REF_DESC, d.refDesc, null );
            if ( isEditable() || d.hasReference() ) {
                addSection( page, NodeDataDraft.SEC_REFERENCE, null, g, d.hasReference() );
            }
        }
        // -- Binary characters (view only) --
        if ( !isEditable() && nd.isHasBinaryCharacters() ) {
            final BinaryCharacters bc = nd.getBinaryCharacters();
            final Grid g = new Grid();
            viewRow( g, "Present", bc.getPresentCount() + "  " + bc.getPresentCharactersAsStringBuffer() );
            viewRow( g, "Gained", bc.getGainedCount() + "  " + bc.getGainedCharactersAsStringBuffer() );
            viewRow( g, "Lost", bc.getLostCount() + "  " + bc.getLostCharactersAsStringBuffer() );
            addSection( page, "Binary characters", null, g, true );
        }
        // -- Properties --
        {
            _property_model = new PropertyTableModel( d.properties, isEditable() );
            _property_table = new JTable( _property_model );
            _property_table.setSelectionMode( ListSelectionModel.SINGLE_SELECTION );
            _property_table.setRowHeight( _property_table.getFontMetrics( _property_table.getFont() ).getHeight()
                    + 8 );
            _property_table.setAutoResizeMode( JTable.AUTO_RESIZE_ALL_COLUMNS );
            _property_table.putClientProperty( "terminateEditOnFocusLost", Boolean.TRUE );
            _property_table.setSurrendersFocusOnKeystroke( true );
            _property_table.setShowGrid( false );
            _property_table.setIntercellSpacing( new Dimension( 6, 0 ) );
            final int[] widths = { 150, 140, 90, 100, 80 };
            for( int i = 0; i < widths.length; ++i ) {
                final TableColumn col = _property_table.getColumnModel().getColumn( i );
                col.setPreferredWidth( widths[ i ] );
            }
            if ( isEditable() ) {
                final JComboBox<String> datatype = new JComboBox<>( NodeDataDraft.PROPERTY_DATATYPES
                        .toArray( new String[ 0 ] ) );
                datatype.setEditable( true );
                _property_table.getColumnModel().getColumn( 3 ).setCellEditor( new DefaultCellEditor( datatype ) );
                final JComboBox<AppliesTo> applies = new JComboBox<>( AppliesTo.values() );
                _property_table.getColumnModel().getColumn( 4 ).setCellEditor( new DefaultCellEditor( applies ) );
            }
            _property_table.setDefaultRenderer( Object.class, new ProblemCellRenderer() );
            _property_model.addTableModelListener( e -> {
                final Section sec = _sections.get( NodeDataDraft.SEC_PROPERTIES );
                if ( sec != null ) {
                    sec.setDetail( propertyDetail() );
                }
                // a cell commit by the table itself (Enter, focus lost) fires this while its editor is still
                // installed -- collect() must not try to stop that editor again (see commitTableEdits)
                final boolean was = _committing_table_edit;
                _committing_table_edit = true;
                try {
                    fireChanged();
                }
                finally {
                    _committing_table_edit = was;
                }
            } );
            final JPanel holder = new JPanel( new BorderLayout( 0, 4 ) );
            holder.setOpaque( false );
            final JPanel table_panel = new JPanel( new BorderLayout() );
            table_panel.add( _property_table.getTableHeader(), BorderLayout.NORTH );
            table_panel.add( _property_table, BorderLayout.CENTER );
            table_panel.setBorder( BorderFactory.createLineBorder( borderColor() ) );
            holder.add( table_panel, BorderLayout.CENTER );
            if ( isEditable() ) {
                final JPanel buttons = new JPanel( new FlowLayout( FlowLayout.LEFT, 4, 0 ) );
                buttons.setOpaque( false );
                buttons.add( addButton( "+ Add property", () -> {
                    _property_model.add( new PropertyDraft() );
                    final int row = _property_model.getRowCount() - 1;
                    _property_table.setRowSelectionInterval( row, row );
                    _property_table.editCellAt( row, 0 );
                    final Component ed = _property_table.getEditorComponent();
                    if ( ed != null ) {
                        ed.requestFocusInWindow();
                    }
                } ) );
                _remove_property_button = addButton( "− Remove property", () -> {
                    final int row = _property_table.getSelectedRow();
                    if ( row >= 0 ) {
                        commitTableEdits();
                        _property_model.remove( row );
                    }
                } );
                _remove_property_button.setEnabled( false );
                _property_table.getSelectionModel().addListSelectionListener( e -> _remove_property_button
                        .setEnabled( _property_table.getSelectedRow() >= 0 ) );
                buttons.add( _remove_property_button );
                holder.add( buttons, BorderLayout.SOUTH );
            }
            if ( isEditable() || !d.properties.isEmpty() ) {
                addSection( page, NodeDataDraft.SEC_PROPERTIES, propertyDetail(), holder, !d.properties.isEmpty() );
            }
        }
        page.add( Box.createVerticalGlue() );
    }

    private void addSection( final JPanel page, final String title, final String detail, final JComponent body,
                             final boolean expanded ) {
        final Section s = new Section( title, detail, body, expanded );
        _sections.put( title, s );
        page.add( s );
    }

    private String sequenceDetail() {
        final int n = _sequence_cards.size();
        return ( n == 0 ) ? ( isEditable() ? "none" : "" ) : String.valueOf( n );
    }

    private String propertyDetail() {
        final int n = ( _property_model == null ) ? 0 : _property_model.getRowCount();
        return ( n == 0 ) ? ( isEditable() ? "none" : "" ) : String.valueOf( n );
    }

    private static List<String> rankChoices() {
        final List<String> out = new ArrayList<>();
        out.add( "" );
        out.addAll( TaxonomyUtil.TAXONOMY_RANKS_LIST );
        return out;
    }

    // ---- row factories (EDIT: a live widget; VIEW: a read-only value, and the row is skipped when empty) ----
    private void addText( final Grid g, final String label, final String key, final String value,
                          final String placeholder ) {
        if ( !isEditable() && value.isEmpty() ) {
            return;
        }
        g.row( label, valueComponent( key, value, placeholder ), false );
    }

    /** Two text fields on one row (EDIT); in VIEW mode each becomes its own row, and only when non-empty. */
    private void addText( final Grid g, final String label, final String key, final String value,
                          final String placeholder, final String label2, final String key2, final String value2,
                          final String placeholder2 ) {
        if ( isEditable() ) {
            g.row( label, valueComponent( key, value, placeholder ), label2,
                   valueComponent( key2, value2, placeholder2 ) );
            return;
        }
        addText( g, label, key, value, placeholder );
        addText( g, label2, key2, value2, placeholder2 );
    }

    /** Two fields on one row; the second is a combo when {@code choices2} is given. In VIEW mode each is its own
     *  row (only when non-empty) so the page reads as a clean list. */
    private void addText( final Grid g, final String label, final String key, final String value,
                          final String placeholder, final String label2, final String key2, final String value2,
                          final List<String> choices2, final boolean editable_combo ) {
        if ( isEditable() ) {
            final JComponent second = ( choices2 != null ) ? comboField( key2, value2, choices2, editable_combo )
                    : valueComponent( key2, value2, null );
            g.row( label, valueComponent( key, value, placeholder ), label2, second );
            return;
        }
        if ( !value.isEmpty() ) {
            g.row( label, valueComponent( key, value, placeholder ), false );
        }
        if ( ( choices2 != null ) && !value2.isEmpty() ) {
            g.row( label2, valueComponent( key2, value2, null ), false );
        }
    }

    private void addArea( final Grid g, final String label, final String key, final String value,
                          final String placeholder, final int rows, final boolean monospace ) {
        if ( !isEditable() && value.isEmpty() ) {
            return;
        }
        g.row( label, areaComponent( key, value, placeholder, rows, monospace ), true );
    }

    private void viewRow( final Grid g, final String label, final String value ) {
        if ( ForesterUtil.isEmpty( value ) ) {
            return;
        }
        g.row( label, viewValue( value, false ), false );
    }

    private JComponent valueComponent( final String key, final String value, final String placeholder ) {
        if ( !isEditable() ) {
            final JComponent c = viewValue( value, false );
            register( key, c );
            return c;
        }
        final JTextField tf = editField( value, placeholder );
        tf.getDocument().addDocumentListener( _doc_listener );
        register( key, tf );
        return tf;
    }

    private JComponent comboField( final String key, final String value, final List<String> choices,
                                   final boolean editable ) {
        final JComboBox<String> cb = new JComboBox<>( choices.toArray( new String[ 0 ] ) );
        cb.setEditable( editable );
        cb.setSelectedItem( value );
        cb.addActionListener( e -> fireChanged() );
        if ( editable ) {
            final Component ed = cb.getEditor().getEditorComponent();
            if ( ed instanceof JTextComponent ) {
                ( (JTextComponent) ed ).getDocument().addDocumentListener( _doc_listener );
            }
        }
        register( key, cb );
        return cb;
    }

    private JComponent areaComponent( final String key, final String value, final String placeholder,
                                      final int rows, final boolean monospace ) {
        if ( !isEditable() ) {
            final JComponent c = viewValue( value, true );
            if ( monospace ) {
                c.setFont( monoFont( c.getFont() ) );
            }
            register( key, c );
            return c;
        }
        final JTextArea ta = new JTextArea( value, rows, 20 );
        ta.setLineWrap( monospace );
        ta.setWrapStyleWord( false );
        if ( monospace ) {
            ta.setFont( monoFont( ta.getFont() ) );
        }
        if ( placeholder != null ) {
            ta.setToolTipText( placeholder );
        }
        ta.getDocument().addDocumentListener( _doc_listener );
        register( key, ta );
        return areaScrollPane( ta, monospace );
    }

    /** An editable single-line field: placeholder hint, and a click puts the caret where you click (no select-all). */
    private static JTextField editField( final String value, final String placeholder ) {
        final JTextField tf = new JTextField( value, 10 ); // a column count: the TEXT must not size the grid
        if ( placeholder != null ) {
            tf.putClientProperty( "JTextField.placeholderText", placeholder );
        }
        tf.putClientProperty( "JTextField.selectAllOnFocusPolicy", "never" );
        return tf;
    }

    /** A read-only value that still looks like text (selectable, copyable), not like a disabled field. */
    private static JComponent viewValue( final String value, final boolean multiline ) {
        final JTextComponent tc;
        if ( multiline ) {
            final JTextArea ta = new JTextArea( value ) {

                private static final long serialVersionUID = 1L;

                @Override
                public Dimension getPreferredSize() {
                    // the HEIGHT follows the wrapped text at the current width; the WIDTH must never be the text's
                    // (a wrapped area reports its last laid-out width, which would ratchet the grid wider than the
                    // window -- and GridBagLayout answers an over-wide grid by collapsing every row to its minimum)
                    final Dimension d = super.getPreferredSize();
                    return new Dimension( VIEW_VALUE_WIDTH, d.height );
                }
            };
            ta.setLineWrap( true );
            ta.setWrapStyleWord( true );
            tc = ta;
        }
        else {
            tc = new JTextField( value ) {

                private static final long serialVersionUID = 1L;

                @Override
                public Dimension getPreferredSize() {
                    final Dimension d = super.getPreferredSize();
                    return new Dimension( Math.min( d.width, VIEW_VALUE_WIDTH ), d.height );
                }
            };
        }
        tc.setFocusable( false ); // read-only: no caret, no focus ring; the value is display, not input
        tc.setEditable( false );
        tc.setOpaque( false );
        tc.setBorder( BorderFactory.createEmptyBorder( 2, 0, 2, 0 ) );
        tc.setForeground( UIManager.getColor( "Label.foreground" ) );
        if ( tc.getCaret() instanceof DefaultCaret ) {
            // a read-only value must never scroll the page to itself (a long mol seq would open the window scrolled)
            ( (DefaultCaret) tc.getCaret() ).setUpdatePolicy( DefaultCaret.NEVER_UPDATE );
        }
        tc.setCaretPosition( 0 );
        return tc;
    }

    private JButton addButton( final String text, final Runnable action ) {
        final JButton b = new JButton( text );
        b.putClientProperty( "JButton.buttonType", "toolBarButton" );
        b.setHorizontalAlignment( SwingConstants.LEFT );
        b.setFocusable( false );
        b.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
        b.addActionListener( e -> action.run() );
        return b;
    }

    private JButton removeButton( final String tooltip, final Runnable action ) {
        final JButton b = new JButton( "×" );
        b.putClientProperty( "JButton.buttonType", "toolBarButton" );
        b.setToolTipText( tooltip );
        b.setFocusable( false );
        b.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
        b.addActionListener( e -> action.run() );
        return b;
    }

    // ---- confidences ----
    private void addConfidenceRow( final ConfidenceDraft c, final boolean by_user ) {
        final ConfidenceRow row = new ConfidenceRow( c );
        _confidence_rows.add( row );
        _confidence_list.add( row );
        if ( by_user ) {
            revalidatePage();
            row.focus();
            fireChanged();
        }
    }

    private void removeConfidenceRow( final ConfidenceRow row ) {
        _confidence_rows.remove( row );
        _confidence_list.remove( row );
        revalidatePage();
        fireChanged();
    }

    /** value [type] ± sd [×] on one line (EDIT) / "95 (bootstrap)" (VIEW). */
    private final class ConfidenceRow extends JPanel {

        private static final long serialVersionUID = 1L;
        private final JComponent  _value;
        private final JComponent  _type;
        private final JComponent  _sd;

        ConfidenceRow( final ConfidenceDraft c ) {
            super( new GridBagLayout() );
            setOpaque( false );
            setAlignmentX( LEFT_ALIGNMENT );
            final GridBagConstraints gc = new GridBagConstraints();
            gc.insets = new Insets( 2, 0, 2, 6 );
            gc.anchor = GridBagConstraints.WEST;
            if ( !isEditable() ) {
                final StringBuilder sb = new StringBuilder( c.value );
                if ( !c.type.isEmpty() ) {
                    sb.append( " (" ).append( c.type ).append( ")" );
                }
                if ( !c.sd.isEmpty() ) {
                    sb.append( " ± " ).append( c.sd );
                }
                _value = viewValue( sb.toString(), false );
                _type = null;
                _sd = null;
                gc.weightx = 1;
                gc.fill = GridBagConstraints.HORIZONTAL;
                add( _value, gc );
                return;
            }
            final JTextField value = editField( c.value, "value" );
            value.setColumns( 7 );
            value.getDocument().addDocumentListener( _doc_listener );
            _value = value;
            final JComboBox<String> type = new JComboBox<>( NodeDataDraft.CONFIDENCE_TYPES.toArray( new String[ 0 ] ) );
            type.setEditable( true );
            type.setSelectedItem( c.type );
            type.addActionListener( e -> fireChanged() );
            final Component ed = type.getEditor().getEditorComponent();
            if ( ed instanceof JTextComponent ) {
                ( (JTextComponent) ed ).getDocument().addDocumentListener( _doc_listener );
                ( (JTextComponent) ed ).putClientProperty( "JTextField.placeholderText", "type" );
            }
            _type = type;
            final JTextField sd = editField( c.sd, "sd" );
            sd.setColumns( 5 );
            sd.getDocument().addDocumentListener( _doc_listener );
            _sd = sd;
            gc.gridx = 0;
            add( value, gc );
            gc.gridx = 1;
            gc.weightx = 1;
            gc.fill = GridBagConstraints.HORIZONTAL;
            add( type, gc );
            gc.gridx = 2;
            gc.weightx = 0;
            gc.fill = GridBagConstraints.NONE;
            add( new JLabel( "±" ), gc );
            gc.gridx = 3;
            add( sd, gc );
            gc.gridx = 4;
            gc.insets = new Insets( 2, 0, 2, 0 );
            add( removeButton( "Remove this confidence", () -> removeConfidenceRow( this ) ), gc );
        }

        /** EDIT mode only (a VIEW row shows one formatted text and is never read back). */
        ConfidenceDraft toDraft() {
            return new ConfidenceDraft( valueOf( _value ), valueOf( _type ), valueOf( _sd ) );
        }

        JComponent field( final String name ) {
            switch ( name ) {
                case NodeDataDraft.CONF_VALUE:
                    return _value;
                case NodeDataDraft.CONF_TYPE:
                    return _type;
                case NodeDataDraft.CONF_SD:
                    return _sd;
                default:
                    return null;
            }
        }

        void focus() {
            _value.requestFocusInWindow();
        }
    }

    // ---- sequences ----
    private void addSequenceCard( final SequenceDraft s, final boolean by_user ) {
        final SequenceCard card = new SequenceCard( s );
        _sequence_cards.add( card );
        _sequence_list.add( card );
        _sequence_list.add( Box.createVerticalStrut( 8 ) );
        renumberSequenceCards();
        if ( by_user ) {
            revalidatePage();
            card.focus();
            fireChanged();
        }
    }

    private void removeSequenceCard( final SequenceCard card ) {
        _sequence_cards.remove( card );
        final Component[] comps = _sequence_list.getComponents();
        for( int i = 0; i < comps.length; ++i ) {
            if ( comps[ i ] == card ) {
                _sequence_list.remove( card );
                if ( ( i < comps.length - 1 ) ) {
                    _sequence_list.remove( comps[ i + 1 ] ); // its strut
                }
                break;
            }
        }
        renumberSequenceCards();
        revalidatePage();
        fireChanged();
    }

    private void renumberSequenceCards() {
        for( int i = 0; i < _sequence_cards.size(); ++i ) {
            _sequence_cards.get( i ).setNumber( i + 1, _sequence_cards.size() );
        }
        final Section s = _sections.get( NodeDataDraft.SEC_SEQUENCES );
        if ( s != null ) {
            s.setDetail( sequenceDetail() );
        }
    }

    /** One {@code <sequence>}: a bordered card with its own field grid and a remove button. */
    private final class SequenceCard extends JPanel {

        private static final long            serialVersionUID = 1L;
        private final Map<String, JComponent> _f              = new HashMap<>();
        private final JLabel                 _title;
        private final JLabel                 _name_hint;
        private JCheckBox                    _aligned;
        private JLabel                       _length_hint;
        private Sequence                     _origin;

        SequenceCard( final SequenceDraft s ) {
            super( new BorderLayout( 0, 6 ) );
            _origin = s.origin;
            setAlignmentX( LEFT_ALIGNMENT );
            setBorder( BorderFactory.createCompoundBorder( BorderFactory.createLineBorder( borderColor(), 1, true ),
                                                           BorderFactory.createEmptyBorder( 8, 10, 8, 10 ) ) );
            setOpaque( false );
            final JPanel head = new JPanel( new BorderLayout( 8, 0 ) );
            head.setOpaque( false );
            _title = new JLabel( "Sequence" );
            _title.setFont( _title.getFont().deriveFont( Font.BOLD ) );
            _name_hint = new JLabel( s.name );
            _name_hint.setForeground( mutedColor() );
            final JPanel left = new JPanel( new FlowLayout( FlowLayout.LEFT, 8, 0 ) );
            left.setOpaque( false );
            left.add( _title );
            left.add( _name_hint );
            head.add( left, BorderLayout.CENTER );
            if ( isEditable() ) {
                head.add( removeButton( "Remove this sequence", () -> removeSequenceCard( this ) ),
                          BorderLayout.EAST );
            }
            add( head, BorderLayout.NORTH );
            final Grid g = new Grid();
            addCardText( g, "Name", NodeDataDraft.SEQ_NAME, s.name, null, "Symbol", NodeDataDraft.SEQ_SYMBOL,
                         s.symbol, "e.g. BRCA1" );
            if ( isEditable() ) {
                final JComboBox<String> type = new JComboBox<>( new String[] { "", PhyloXmlUtil.SEQ_TYPE_PROTEIN,
                        PhyloXmlUtil.SEQ_TYPE_DNA, PhyloXmlUtil.SEQ_TYPE_RNA } );
                type.setSelectedItem( s.type );
                type.addActionListener( e -> fireChanged() );
                _f.put( NodeDataDraft.SEQ_TYPE, type );
                _aligned = new JCheckBox( "Aligned", s.aligned );
                _aligned.setOpaque( false );
                _aligned.setToolTipText( "The molecular sequence is part of an alignment (is_aligned)" );
                _aligned.addItemListener( e -> fireChanged() );
                final JPanel type_row = new JPanel( new FlowLayout( FlowLayout.LEFT, 8, 0 ) );
                type_row.setOpaque( false );
                type_row.add( type );
                type_row.add( _aligned );
                g.row( "Gene name", cardValue( NodeDataDraft.SEQ_GENE, s.geneName, null ), "Type", type_row );
            }
            else {
                addCardText( g, "Gene name", NodeDataDraft.SEQ_GENE, s.geneName, null );
                addCardText( g, "Type", NodeDataDraft.SEQ_TYPE, s.type + ( s.aligned ? " (aligned)" : "" ), null );
            }
            addCardText( g, "Accession", NodeDataDraft.SEQ_ACC, s.accession, "e.g. P38398", "Source",
                         NodeDataDraft.SEQ_SOURCE, s.source, "e.g. UniProt" );
            addCardText( g, "Location", NodeDataDraft.SEQ_LOCATION, s.location, "e.g. chr17:43044295-43125483" );
            // molecular sequence, with a live residue count beside the label
            if ( isEditable() || !s.molSeq.isEmpty() ) {
                final JComponent area = cardArea( NodeDataDraft.SEQ_MOL_SEQ, s.molSeq, MOL_SEQ_ROWS, true );
                _length_hint = new JLabel();
                _length_hint.setForeground( mutedColor() );
                _length_hint.setFont( _length_hint.getFont().deriveFont( _length_hint.getFont().getSize2D() * 0.9f ) );
                final JPanel label = new JPanel();
                label.setOpaque( false );
                label.setLayout( new BoxLayout( label, BoxLayout.Y_AXIS ) );
                label.add( new JLabel( "Mol seq" ) );
                label.add( _length_hint );
                g.row( label, area, true );
                updateLengthHint();
            }
            addCardArea( g, "URIs", NodeDataDraft.SEQ_URIS, s.uris, 2 );
            if ( !isEditable() && ( _origin != null ) ) {
                addViewExtras( g, _origin );
            }
            add( g, BorderLayout.CENTER );
        }

        private void addViewExtras( final Grid g, final Sequence seq ) {
            if ( seq.getPrimaryAccession() != null ) {
                viewRow( g, "Primary accession", seq.getPrimaryAccession().asText().toString() );
            }
            final SortedSet<Annotation> anns = seq.getAnnotations();
            if ( ( anns != null ) && !anns.isEmpty() ) {
                final List<String> lines = new ArrayList<>();
                for( final Annotation a : anns ) {
                    final StringBuilder sb = new StringBuilder( a.asText() );
                    final List<String> attrs = new ArrayList<>();
                    if ( !ForesterUtil.isEmpty( a.getSource() ) ) {
                        attrs.add( "source " + a.getSource() );
                    }
                    if ( !ForesterUtil.isEmpty( a.getType() ) ) {
                        attrs.add( "type " + a.getType() );
                    }
                    if ( !ForesterUtil.isEmpty( a.getEvidence() ) ) {
                        attrs.add( "evidence " + a.getEvidence() );
                    }
                    if ( a.getConfidence() != null ) {
                        attrs.add( "confidence " + a.getConfidence().asText() );
                    }
                    if ( !attrs.isEmpty() ) {
                        sb.append( " [" ).append( String.join( ", ", attrs ) ).append( "]" );
                    }
                    lines.add( sb.toString() );
                }
                g.row( "Annotations", viewValue( String.join( "\n", lines ), true ), true );
            }
            final SortedSet<Accession> xrefs = seq.getCrossReferences();
            if ( ( xrefs != null ) && !xrefs.isEmpty() ) {
                final List<String> lines = new ArrayList<>();
                for( final Accession x : xrefs ) {
                    lines.add( x.asText().toString() );
                }
                g.row( "Cross references", viewValue( String.join( "\n", lines ), true ), true );
            }
            final DomainArchitecture da = seq.getDomainArchitecture();
            if ( ( da != null ) && ( da.getNumberOfDomains() > 0 ) ) {
                final List<String> names = new ArrayList<>();
                for( int i = 0; i < da.getNumberOfDomains(); ++i ) {
                    names.add( da.getDomain( i ).getName() );
                }
                viewRow( g, "Domains", da.getNumberOfDomains() + ": " + String.join( ", ", names ) );
            }
        }

        private void addCardText( final Grid g, final String label, final String key, final String value,
                                  final String placeholder ) {
            if ( !isEditable() && value.isEmpty() ) {
                return;
            }
            g.row( label, cardValue( key, value, placeholder ), false );
        }

        private void addCardText( final Grid g, final String label, final String key, final String value,
                                  final String placeholder, final String label2, final String key2,
                                  final String value2, final String placeholder2 ) {
            if ( isEditable() ) {
                g.row( label, cardValue( key, value, placeholder ), label2, cardValue( key2, value2, placeholder2 ) );
                return;
            }
            addCardText( g, label, key, value, placeholder );
            addCardText( g, label2, key2, value2, placeholder2 );
        }

        private void addCardArea( final Grid g, final String label, final String key, final String value,
                                  final int rows ) {
            if ( !isEditable() && value.isEmpty() ) {
                return;
            }
            g.row( label, cardArea( key, value, rows, false ), true );
        }

        private JComponent cardValue( final String key, final String value, final String placeholder ) {
            if ( !isEditable() ) {
                final JComponent c = viewValue( value, false );
                _f.put( key, c );
                return c;
            }
            final JTextField tf = editField( value, placeholder );
            tf.getDocument().addDocumentListener( _doc_listener );
            if ( NodeDataDraft.SEQ_NAME.equals( key ) ) {
                tf.getDocument().addDocumentListener( onChange( () -> _name_hint.setText( tf.getText() ) ) );
            }
            _f.put( key, tf );
            return tf;
        }

        private JComponent cardArea( final String key, final String value, final int rows, final boolean mono ) {
            if ( !isEditable() ) {
                final JComponent c = viewValue( value, true );
                _f.put( key, c );
                if ( !mono ) {
                    return c;
                }
                // a molecular sequence can run to thousands of residues: show it in a capped, scrolling box
                c.setFont( monoFont( c.getFont() ) );
                ( (JTextArea) c ).setRows( rows * 2 );
                return areaScrollPane( (JTextArea) c, true );
            }
            final JTextArea ta = new JTextArea( value, rows, 20 );
            ta.setLineWrap( mono );
            if ( mono ) {
                ta.setFont( monoFont( ta.getFont() ) );
            }
            ta.getDocument().addDocumentListener( _doc_listener );
            if ( NodeDataDraft.SEQ_MOL_SEQ.equals( key ) ) {
                ta.getDocument().addDocumentListener( onChange( this::updateLengthHint ) );
            }
            _f.put( key, ta );
            return areaScrollPane( ta, mono );
        }

        private void updateLengthHint() {
            if ( _length_hint == null ) {
                return;
            }
            final int n = NodeDataDraft.countResidues( valueOf( _f.get( NodeDataDraft.SEQ_MOL_SEQ ) ) );
            _length_hint.setText( ( n == 0 ) ? "" : n + " residues" );
        }

        void setNumber( final int n, final int of ) {
            _title.setText( ( of > 1 ) || isEditable() ? "Sequence " + n : "Sequence" );
        }

        void bindOrigin( final Sequence s ) {
            _origin = s;
        }

        /** EDIT mode only (a VIEW card is never read back). */
        SequenceDraft toDraft() {
            final SequenceDraft d = SequenceDraft.withOrigin( _origin ); // identity only; the values are the widgets'
            d.name = valueOf( _f.get( NodeDataDraft.SEQ_NAME ) );
            d.symbol = valueOf( _f.get( NodeDataDraft.SEQ_SYMBOL ) );
            d.geneName = valueOf( _f.get( NodeDataDraft.SEQ_GENE ) );
            d.type = valueOf( _f.get( NodeDataDraft.SEQ_TYPE ) );
            d.aligned = _aligned.isSelected();
            d.accession = valueOf( _f.get( NodeDataDraft.SEQ_ACC ) );
            d.source = valueOf( _f.get( NodeDataDraft.SEQ_SOURCE ) );
            d.location = valueOf( _f.get( NodeDataDraft.SEQ_LOCATION ) );
            d.molSeq = valueOf( _f.get( NodeDataDraft.SEQ_MOL_SEQ ) );
            d.uris = valueOf( _f.get( NodeDataDraft.SEQ_URIS ) );
            return d;
        }

        JComponent field( final String name ) {
            if ( NodeDataDraft.SEQ_ALIGNED.equals( name ) ) {
                return _aligned;
            }
            return _f.get( name );
        }

        void focus() {
            final JComponent c = _f.get( NodeDataDraft.SEQ_NAME );
            if ( c != null ) {
                c.requestFocusInWindow();
            }
        }
    }

    // ---- properties table ----
    private static final class PropertyTableModel extends AbstractTableModel {

        private static final long   serialVersionUID = 1L;
        private static final String[] COLUMNS        = { "Reference", "Value", "Unit", "Datatype", "Applies to" };
        private final List<PropertyDraft> _rows      = new ArrayList<>();
        private final boolean       _editable;

        PropertyTableModel( final List<PropertyDraft> rows, final boolean editable ) {
            for( final PropertyDraft p : rows ) {
                _rows.add( p.copy() );
            }
            _editable = editable;
        }

        List<PropertyDraft> rows() {
            return _rows;
        }

        void add( final PropertyDraft p ) {
            _rows.add( p );
            fireTableRowsInserted( _rows.size() - 1, _rows.size() - 1 );
        }

        void remove( final int row ) {
            _rows.remove( row );
            fireTableRowsDeleted( row, row );
        }

        @Override
        public int getRowCount() {
            return _rows.size();
        }

        @Override
        public int getColumnCount() {
            return COLUMNS.length;
        }

        @Override
        public String getColumnName( final int col ) {
            return COLUMNS[ col ];
        }

        @Override
        public boolean isCellEditable( final int row, final int col ) {
            return _editable;
        }

        @Override
        public Object getValueAt( final int row, final int col ) {
            final PropertyDraft p = _rows.get( row );
            switch ( col ) {
                case 0:
                    return p.ref;
                case 1:
                    return p.value;
                case 2:
                    return p.unit;
                case 3:
                    return p.datatype;
                default:
                    return p.appliesTo;
            }
        }

        @Override
        public void setValueAt( final Object value, final int row, final int col ) {
            final PropertyDraft p = _rows.get( row );
            final String s = ( value == null ) ? "" : String.valueOf( value ).trim();
            switch ( col ) {
                case 0:
                    p.ref = s;
                    break;
                case 1:
                    p.value = s;
                    break;
                case 2:
                    p.unit = s;
                    break;
                case 3:
                    p.datatype = s;
                    break;
                default:
                    p.appliesTo = ( value instanceof AppliesTo ) ? (AppliesTo) value : AppliesTo.NODE;
            }
            fireTableCellUpdated( row, col );
        }
    }

    /** Which property cells currently have a validation problem ("row:col" -> message); painted by the renderer. */
    private final Map<String, String> _bad_property_cells = new HashMap<>();

    private static int propertyColumn( final String field ) {
        switch ( field ) {
            case NodeDataDraft.PROP_REF:
                return 0;
            case NodeDataDraft.PROP_VALUE:
                return 1;
            case NodeDataDraft.PROP_UNIT:
                return 2;
            case NodeDataDraft.PROP_TYPE:
                return 3;
            default:
                return 4;
        }
    }

    /** Paints a cell with a problem in the error colour and shows the problem as its tooltip. */
    private final class ProblemCellRenderer extends javax.swing.table.DefaultTableCellRenderer {

        private static final long serialVersionUID = 1L;

        @Override
        public Component getTableCellRendererComponent( final JTable table, final Object value,
                                                        final boolean selected, final boolean focus, final int row,
                                                        final int col ) {
            final Component c = super.getTableCellRendererComponent( table, value, selected, focus, row, col );
            final String problem = _bad_property_cells.get( row + ":" + col );
            if ( problem != null ) {
                c.setForeground( errorColor() );
                if ( c instanceof JComponent ) {
                    ( (JComponent) c ).setToolTipText( problem );
                }
            }
            else {
                c.setForeground( selected ? table.getSelectionForeground() : table.getForeground() );
                if ( c instanceof JComponent ) {
                    ( (JComponent) c ).setToolTipText( null );
                }
            }
            return c;
        }
    }

    /** For tests: whether the property cell at {@code row}/{@code field} is currently marked as a problem. */
    boolean isPropertyCellMarkedForTest( final int row, final String field ) {
        return _bad_property_cells.containsKey( row + ":" + propertyColumn( field ) );
    }

    /** Set while a property cell edit is being committed (by us or by the table itself), see commitTableEdits. */
    private boolean _committing_table_edit = false;

    /**
     * Pushes a property cell that is still being edited into the model, so collect() sees what the user typed.
     * NOT re-entrant: JTable.editingStopped calls setValueAt BEFORE it removes the editor, so the model event it
     * fires (-> fireChanged -> collect -> here) arrives while isEditing() is still true -- stopping the editor
     * again from inside that would recurse without end (a StackOverflowError on the first property edit).
     */
    private void commitTableEdits() {
        if ( _committing_table_edit || ( _property_table == null ) || !_property_table.isEditing() ) {
            return;
        }
        _committing_table_edit = true;
        try {
            _property_table.getCellEditor().stopCellEditing();
        }
        finally {
            _committing_table_edit = false;
        }
    }

    // ---- collapsible section ----
    private final class Section extends JPanel {

        private static final long serialVersionUID = 1L;
        private final JComponent  _body;
        private final JButton     _header;
        private final JLabel      _detail;
        private boolean           _expanded;

        Section( final String title, final String detail, final JComponent body, final boolean expanded ) {
            super( new BorderLayout() );
            setOpaque( false );
            setAlignmentX( LEFT_ALIGNMENT );
            _body = body;
            _expanded = expanded;
            _header = new JButton( title );
            _header.putClientProperty( "JButton.buttonType", "toolBarButton" );
            _header.setHorizontalAlignment( SwingConstants.LEFT );
            _header.setFont( _header.getFont().deriveFont( Font.BOLD ) );
            _header.setFocusable( false );
            _header.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
            _header.addActionListener( e -> toggle() );
            _detail = new JLabel( ( detail == null ) ? "" : detail );
            _detail.setForeground( mutedColor() );
            final JPanel head = new JPanel( new BorderLayout( 6, 0 ) );
            head.setOpaque( false );
            head.add( _header, BorderLayout.WEST );
            head.add( _detail, BorderLayout.CENTER );
            head.setBorder( BorderFactory.createCompoundBorder( BorderFactory
                    .createMatteBorder( 1, 0, 0, 0, borderColor() ), BorderFactory.createEmptyBorder( 4, 0, 4, 0 ) ) );
            head.addMouseListener( new MouseAdapter() {

                @Override
                public void mouseClicked( final MouseEvent e ) {
                    toggle();
                }
            } );
            head.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
            add( head, BorderLayout.NORTH );
            final JPanel body_wrap = new JPanel( new BorderLayout() );
            body_wrap.setOpaque( false );
            body_wrap.setBorder( BorderFactory.createEmptyBorder( 2, SECTION_INDENT, 10, 0 ) );
            body_wrap.add( body, BorderLayout.CENTER );
            add( body_wrap, BorderLayout.CENTER );
            body_wrap.setVisible( expanded );
            updateIcon();
        }

        boolean isExpanded() {
            return _expanded;
        }

        void toggle() {
            _expanded = !_expanded;
            ( (JComponent) getComponent( 1 ) ).setVisible( _expanded );
            updateIcon();
            revalidatePage();
        }

        void setDetail( final String detail ) {
            _detail.setText( ( detail == null ) ? "" : detail );
        }

        private void updateIcon() {
            final Icon icon = UIManager.getIcon( _expanded ? "Tree.expandedIcon" : "Tree.collapsedIcon" );
            _header.setIcon( icon );
            if ( icon == null ) {
                _header.setText( ( _expanded ? "▾ " : "▸ " ) + _header.getText().replaceFirst( "^[▾▸] ", "" ) );
            }
        }

        @Override
        public Dimension getMaximumSize() {
            final Dimension d = super.getPreferredSize();
            return new Dimension( Integer.MAX_VALUE, d.height );
        }
    }

    // ---- the label/field grid used inside every section ----
    private final class Grid extends JPanel {

        private static final long serialVersionUID = 1L;
        private int               _row             = 0;

        Grid() {
            super( new GridBagLayout() );
            setOpaque( false );
        }

        void row( final String label, final JComponent field, final boolean top_aligned ) {
            row( makeLabel( label ), field, top_aligned );
        }

        void row( final JComponent label, final JComponent field, final boolean top_aligned ) {
            final GridBagConstraints c = base( top_aligned );
            c.gridx = 0;
            c.weightx = 0;
            add( label, c );
            c.gridx = 1;
            c.gridwidth = 3;
            c.weightx = 1;
            c.fill = GridBagConstraints.HORIZONTAL;
            c.insets = new Insets( 3, 0, 3, 0 );
            add( field, c );
            _row++;
        }

        void row( final String label, final JComponent field, final String label2, final JComponent field2 ) {
            final GridBagConstraints c = base( false );
            c.gridx = 0;
            c.weightx = 0;
            add( makeLabel( label ), c );
            c.gridx = 1;
            c.weightx = 0.6;
            c.fill = GridBagConstraints.HORIZONTAL;
            c.insets = new Insets( 3, 0, 3, 0 );
            add( field, c );
            c.gridx = 2;
            c.weightx = 0;
            c.fill = GridBagConstraints.NONE;
            c.insets = new Insets( 3, 14, 3, 8 );
            final JLabel l2 = new JLabel( label2 );
            add( l2, c );
            c.gridx = 3;
            c.weightx = 0.4;
            c.fill = GridBagConstraints.HORIZONTAL;
            c.insets = new Insets( 3, 0, 3, 0 );
            add( field2, c );
            _row++;
        }

        private GridBagConstraints base( final boolean top_aligned ) {
            final GridBagConstraints c = new GridBagConstraints();
            c.gridy = _row;
            c.insets = new Insets( 3, 0, 3, 8 );
            c.anchor = top_aligned ? GridBagConstraints.NORTHWEST : GridBagConstraints.WEST;
            return c;
        }

        private JLabel makeLabel( final String text ) {
            final JLabel l = new JLabel( text );
            final Dimension d = l.getPreferredSize();
            l.setPreferredSize( new Dimension( Math.max( _label_width, d.width ), d.height ) );
            return l;
        }
    }

    /** The scrolling page: tracks the viewport WIDTH (so nothing scrolls sideways) but not the height. */
    private static final class ScrollablePage extends JPanel implements Scrollable {

        private static final long serialVersionUID = 1L;

        @Override
        public Dimension getPreferredScrollableViewportSize() {
            return getPreferredSize();
        }

        @Override
        public int getScrollableUnitIncrement( final Rectangle r, final int o, final int d ) {
            return 16;
        }

        @Override
        public int getScrollableBlockIncrement( final Rectangle r, final int o, final int d ) {
            return Math.max( 16, r.height - 16 );
        }

        @Override
        public boolean getScrollableTracksViewportWidth() {
            return true;
        }

        @Override
        public boolean getScrollableTracksViewportHeight() {
            return false;
        }
    }

    // ------------------------------------------------------------------ plumbing
    private void register( final String key, final JComponent c ) {
        _fields.put( key, c );
    }

    private JComponent fieldFor( final String key ) {
        if ( key.startsWith( "sequence." ) ) {
            final String[] parts = key.split( "\\.", 3 );
            final int i = Integer.parseInt( parts[ 1 ] );
            return ( i < _sequence_cards.size() ) ? _sequence_cards.get( i ).field( parts[ 2 ] ) : null;
        }
        if ( key.startsWith( "confidence." ) ) {
            final String[] parts = key.split( "\\.", 3 );
            final int i = Integer.parseInt( parts[ 1 ] );
            return ( i < _confidence_rows.size() ) ? _confidence_rows.get( i ).field( parts[ 2 ] ) : null;
        }
        if ( key.startsWith( "property." ) ) {
            return _property_table; // problems on properties are marked per CELL (see applyProblems)
        }
        return _fields.get( key );
    }

    private String text( final String key ) {
        return valueOf( _fields.get( key ) );
    }

    private static String valueOf( final JComponent c ) {
        if ( c == null ) {
            return "";
        }
        if ( c instanceof JTextComponent ) {
            return NodeDataDraft.nn( ( (JTextComponent) c ).getText() );
        }
        if ( c instanceof JComboBox ) {
            final JComboBox<?> cb = (JComboBox<?>) c;
            if ( cb.isEditable() ) {
                return NodeDataDraft.nn( String.valueOf( cb.getEditor().getItem() ) ).trim();
            }
            final Object o = cb.getSelectedItem();
            return ( o == null ) ? "" : String.valueOf( o );
        }
        if ( c instanceof JCheckBox ) {
            return String.valueOf( ( (JCheckBox) c ).isSelected() );
        }
        return "";
    }

    private static void setValue( final JComponent c, final String value ) {
        if ( c instanceof JTextComponent ) {
            ( (JTextComponent) c ).setText( value );
        }
        else if ( c instanceof JComboBox ) {
            ( (JComboBox<?>) c ).setSelectedItem( value );
        }
        else if ( c instanceof AbstractButton ) {
            ( (AbstractButton) c ).setSelected( Boolean.parseBoolean( value ) );
        }
    }

    private void fireChanged() {
        _current = null;
        _current_problems = null;
        if ( _building ) {
            return;
        }
        applyProblems( problems() ); // reads the widgets ONCE; listeners then see the cached draft/problems
        for( final Runnable r : _change_listeners ) {
            r.run();
        }
    }

    /** Outlines each problem's widget (FlatLaf {@code JComponent.outline}) and clears the outlines that healed. */
    private void applyProblems( final List<Problem> ps ) {
        for( final JComponent c : _outlined ) {
            c.putClientProperty( "JComponent.outline", null );
            c.setToolTipText( null );
        }
        _outlined.clear();
        _bad_property_cells.clear();
        for( final Problem p : ps ) {
            if ( p.key.startsWith( "property." ) ) {
                final String[] parts = p.key.split( "\\.", 3 );
                _bad_property_cells.put( parts[ 1 ] + ":" + propertyColumn( parts[ 2 ] ), p.message );
                continue;
            }
            final JComponent target = outlineTarget( fieldFor( p.key ) );
            if ( target == null ) {
                continue;
            }
            if ( target.getClientProperty( "JComponent.outline" ) == null ) {
                target.putClientProperty( "JComponent.outline", "error" );
                target.setToolTipText( p.message );
                _outlined.add( target );
            }
        }
        repaint();
    }

    /** The component that carries the outline for a field: a text area's scroll pane, otherwise the field. */
    private static JComponent outlineTarget( final JComponent field ) {
        if ( field == null ) {
            return null;
        }
        if ( field instanceof JTextArea ) {
            final JScrollPane sp = (JScrollPane) SwingUtilities.getAncestorOfClass( JScrollPane.class, field );
            return ( sp != null ) ? sp : field;
        }
        return field;
    }

    private void revalidatePage() {
        revalidate();
        repaint();
    }

    /** Refreshes the tree panel after a write -- only the parts the CHANGED sections can affect (each of these
     *  is a full-tree scan, so renaming one node must not rebuild every combo box). */
    private void refreshTreePanel( final Set<String> sections ) {
        final TreePanel tp = _tree_panel;
        final ControlPanel cp = tp.getControlPanel();
        final boolean basic = sections.contains( NodeDataDraft.SEC_BASIC );
        final boolean properties = sections.contains( NodeDataDraft.SEC_PROPERTIES );
        if ( basic ) {
            tp.recalculateMaxDistanceToRoot(); // a branch length may have changed
            if ( ( tp.getPhylogeny() != null ) && ( cp != null ) ) {
                AptxUtil.lookAtRealBranchLengthsForAptxControlSettings( tp.getPhylogeny(), cp );
            }
        }
        if ( sections.contains( NodeDataDraft.SEC_DATE ) ) {
            tp.invalidateTimeAxisDerivation();
        }
        if ( cp != null ) {
            // which Display Data checkboxes / search fields exist follows what data the tree now carries (this
            // rebuilds the search fields itself when the set of present data changed)
            cp.updateDataCheckboxVisibility( true );
            if ( properties ) {
                cp.populateColorByPropertyBox(); // a new/removed property ref -> the "Color by" choices
                cp.populateSizeByPropertyBox(); // ... and "Size by"
                cp.populateAncestralPieBox(); // ... and the discrete-trait pies
                cp.rebuildSearchFields( true ); // a property ref is a searchable field even if presence is unchanged
            }
        }
        tp.setEdited( true );
        tp.repaint();
    }

    /** The small scrolling box around a multi-line field: its height is fixed by the area's row count (so a long
     *  value scrolls inside it), it hands the mouse wheel to the page when it cannot scroll itself, and it keeps
     *  that height even if the surrounding grid falls back to minimum sizes. */
    private static JScrollPane areaScrollPane( final JTextArea ta, final boolean wrap ) {
        final JScrollPane sp = new JScrollPane( ta );
        sp.setVerticalScrollBarPolicy( JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED );
        sp.setHorizontalScrollBarPolicy( wrap ? JScrollPane.HORIZONTAL_SCROLLBAR_NEVER
                : JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED );
        final Dimension pref = sp.getPreferredSize();
        sp.setMinimumSize( new Dimension( 60, pref.height ) );
        forwardWheelWhenIdle( sp );
        return sp;
    }

    /** A nested text-area scroll pane hands the wheel to the page when it cannot scroll further itself. */
    private static void forwardWheelWhenIdle( final JScrollPane inner ) {
        inner.addMouseWheelListener( e -> {
            final javax.swing.JScrollBar bar = inner.getVerticalScrollBar();
            final boolean can_scroll = bar.isVisible()
                    && ( ( e.getWheelRotation() < 0 ) ? ( bar.getValue() > bar.getMinimum() )
                            : ( bar.getValue() + bar.getVisibleAmount() < bar.getMaximum() ) );
            if ( !can_scroll ) {
                final Component parent = inner.getParent();
                if ( parent != null ) {
                    final MouseWheelEvent copy = new MouseWheelEvent( parent, e.getID(), e.getWhen(),
                            e.getModifiersEx(), e.getX(), e.getY(), e.getClickCount(), e.isPopupTrigger(),
                            e.getScrollType(), e.getScrollAmount(), e.getWheelRotation() );
                    parent.dispatchEvent( copy );
                    e.consume();
                }
            }
        } );
    }

    private static Color mutedColor() {
        final Color c = UIManager.getColor( "Label.disabledForeground" );
        return ( c != null ) ? c : Color.GRAY;
    }

    private static Color borderColor() {
        final Color c = UIManager.getColor( "Component.borderColor" );
        return ( c != null ) ? c : Color.LIGHT_GRAY;
    }

    /** Error TEXT colour (FlatLaf's action red reads on both themes; the border red is too dim on dark). */
    private static Color errorColor() {
        Color c = UIManager.getColor( "Actions.Red" );
        if ( c == null ) {
            c = UIManager.getColor( "Component.error.focusedBorderColor" );
        }
        return ( c != null ) ? c : new Color( 0xD0342C );
    }

    private static Font monoFont( final Font base ) {
        return new Font( Font.MONOSPACED, Font.PLAIN, base.getSize() );
    }

    private static List<String> nonEmpty( final List<String> in ) {
        final List<String> out = new ArrayList<>();
        for( final String s : in ) {
            if ( !ForesterUtil.isEmpty( s ) ) {
                out.add( s );
            }
        }
        return out;
    }
}
