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

import static org.forester.archaeopteryx.FormWidgets.editField;
import static org.forester.archaeopteryx.FormWidgets.newPage;
import static org.forester.archaeopteryx.FormWidgets.onChange;
import static org.forester.archaeopteryx.FormWidgets.pageScroller;
import static org.forester.archaeopteryx.FormWidgets.viewValue;

import java.awt.BorderLayout;
import java.awt.Point;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.event.DocumentListener;
import javax.swing.text.JTextComponent;

import org.forester.archaeopteryx.FormWidgets.Grid;
import org.forester.archaeopteryx.FormWidgets.Header;
import org.forester.archaeopteryx.FormWidgets.Section;
import org.forester.archaeopteryx.NodeDataDraft.Problem;
import org.forester.archaeopteryx.TreeFacts.Fact;
import org.forester.archaeopteryx.TreeFacts.Group;
import org.forester.phylogeny.Phylogeny;

/**
 * The Tree Properties page: the tree-level metadata a user may edit (name, description, identifier, type,
 * branch-length unit -- see {@link TreePropertiesDraft}) on top, then everything the tree can tell about itself
 * ({@link TreeFacts}: file, structure, branch lengths, support values, annotation coverage, time axis) as
 * read-only sections. Writing goes through {@link #write()}: one undo checkpoint, the changed fields only, the
 * tab label kept in step with the name. The fact sections are recomputed by {@link #refresh()} whenever the tree
 * changes (the panel marks the window stale), and {@link #rebind()} re-reads everything when an undo or redo
 * swapped the tree -- keeping unwritten edits, now measured against the restored tree.
 * <p>
 * Like {@link NodeDataForm} this is a plain {@code JPanel}; the window chrome is {@link EditorFrame}.
 */
final class TreePropertiesForm extends JPanel implements EditorFrame.Form {

    static final String SEC_NAME     = "Name & Description";
    static final String SEC_IDENTITY = "Identifier & Type";
    static final String UNTITLED     = "Untitled tree";

    private static final long             serialVersionUID  = 1L;
    private final TreePanel               _tree_panel;
    private Phylogeny                     _phylogeny;
    private TreePropertiesDraft           _baseline;
    private final Header                  _header;
    private final JPanel                  _page;
    private final JScrollPane             _scroll;
    private final Map<String, JComponent> _fields           = new HashMap<>();
    private final Map<String, Section>    _sections         = new LinkedHashMap<>();
    private final List<Section>           _fact_sections    = new ArrayList<>();
    private final List<Runnable>          _change_listeners = new ArrayList<>();
    private final List<JComponent>        _outlined         = new ArrayList<>();
    private final int                     _label_width;
    private boolean                       _loading          = false;
    private TreePropertiesDraft           _current;
    private List<Problem>                 _current_problems;
    private List<Group>                   _groups           = Collections.emptyList();
    private final DocumentListener        _doc_listener     = onChange( this::fireChanged );

    /** @param tree_panel the panel whose tree is edited (its file, edited flag and time axis feed the facts) */
    TreePropertiesForm( final TreePanel tree_panel ) {
        this( tree_panel.getPhylogeny(), tree_panel );
    }

    /** For tests: a form over {@code phy} with an optional panel (null = no file, no undo, no time axis). */
    TreePropertiesForm( final Phylogeny phy, final TreePanel tree_panel ) {
        super( new BorderLayout() );
        _tree_panel = tree_panel;
        _phylogeny = phy;
        _baseline = TreePropertiesDraft.from( phy );
        _label_width = getFontMetrics( getFont() ).stringWidth( "Branches with support" ) + 12;
        _header = new Header( titleText(), subtitleText() );
        add( _header, BorderLayout.NORTH );
        _page = newPage();
        buildEditableSections( _baseline );
        rebuildFactSections();
        _scroll = pageScroller( _page );
        add( _scroll, BorderLayout.CENTER );
        applyProblems( problems() );
        SwingUtilities.invokeLater( () -> _scroll.getViewport().setViewPosition( new Point( 0, 0 ) ) );
    }

    // ------------------------------------------------------------------ EditorFrame.Form
    @Override
    public JComponent component() {
        return this;
    }

    @Override
    public boolean isEditable() {
        return true;
    }

    @Override
    public void addChangeListener( final Runnable r ) {
        _change_listeners.add( r );
    }

    @Override
    public boolean isDirty() {
        return !collect().normalized().equals( _baseline.normalized() );
    }

    @Override
    public List<Problem> problems() {
        collect();
        if ( _current_problems == null ) {
            _current_problems = _current.validate( _baseline );
        }
        return _current_problems;
    }

    /**
     * Writes the edited fields to the tree if they validate: one undo checkpoint ("Edit Tree Properties"), every
     * field re-applied from the normalized draft, the tab label re-derived when the name changed, the panel
     * marked edited and repainted, and the facts refreshed. No changes: a no-op that returns true. Invalid:
     * returns false and touches nothing.
     */
    @Override
    public boolean write() {
        final TreePropertiesDraft draft = collect();
        final List<Problem> ps = problems();
        applyProblems( ps );
        if ( !ps.isEmpty() ) {
            return false;
        }
        final Set<String> changed = draft.changedFields( _baseline );
        if ( changed.isEmpty() ) {
            return true;
        }
        if ( _tree_panel != null ) {
            _tree_panel.pushUndoCheckpoint( "Edit Tree Properties" ); // Phylogeny.copy() carries all six fields
        }
        draft.writeTo( _phylogeny );
        _baseline = draft.normalized();
        if ( _tree_panel != null ) {
            if ( changed.contains( "name" ) && ( _tree_panel.getMainPanel() != null ) ) {
                _tree_panel.getMainPanel().syncTabTitle( _tree_panel ); // the tab label IS the tree name
            }
            _tree_panel.setEdited( true );
            _tree_panel.repaint();
        }
        fireChanged();
        refresh();
        return true;
    }

    // ------------------------------------------------------------------ API for the frame + tests
    TreePropertiesDraft baseline() {
        return _baseline;
    }

    /** The current widget values as a draft (cached until the next change). */
    TreePropertiesDraft collect() {
        if ( _current == null ) {
            final TreePropertiesDraft d = new TreePropertiesDraft();
            d.name = text( TreePropertiesDraft.NAME );
            d.description = text( TreePropertiesDraft.DESCRIPTION );
            d.idValue = text( TreePropertiesDraft.ID_VALUE );
            d.idProvider = text( TreePropertiesDraft.ID_PROVIDER );
            d.type = text( TreePropertiesDraft.TYPE );
            d.distanceUnit = text( TreePropertiesDraft.DISTANCE_UNIT );
            _current = d;
        }
        return _current;
    }

    /** What the window title calls this tree (its written name, or a placeholder). */
    String titleText() {
        return _baseline.normalized().name.isEmpty() ? UNTITLED : _baseline.normalized().name;
    }

    /** The muted header line: file name, tip count, format, and whether there are unsaved changes. */
    String subtitleText() {
        final StringBuilder sb = new StringBuilder();
        final File f = ( _tree_panel != null ) ? _tree_panel.getTreeFile() : null;
        sb.append( ( f != null ) ? f.getName() : "not saved to a file yet" );
        if ( ( _phylogeny != null ) && !_phylogeny.isEmpty() ) {
            final int tips = _phylogeny.getNumberOfExternalNodes();
            sb.append( " · " ).append( TreeFacts.count( tips ) ).append( tips == 1 ? " tip" : " tips" );
        }
        if ( ( _tree_panel != null ) && _tree_panel.isEdited() ) {
            sb.append( " · unsaved changes" );
        }
        return sb.toString();
    }

    /** Recomputes the read-only fact sections (and the header) for the current tree; edits are untouched. */
    void refresh() {
        _header.setTitle( titleText() );
        _header.setSubtitle( subtitleText() );
        rebuildFactSections();
    }

    /**
     * Re-reads the tree from the panel (an undo/redo may have swapped it): the baseline is the restored tree's
     * metadata; clean widgets are reloaded from it, while unwritten edits are kept and now measured against the
     * restored tree. Then everything is refreshed.
     */
    void rebind() {
        final boolean was_dirty = isDirty();
        if ( ( _tree_panel != null ) && !_tree_panel.isCurrentTreeIsSubtree() ) {
            _phylogeny = _tree_panel.getPhylogeny(); // (a displayed sub-tree is transient: stay on the whole tree)
        }
        _baseline = TreePropertiesDraft.from( _phylogeny );
        if ( !was_dirty ) {
            load( _baseline );
        }
        fireChanged();
        refresh();
    }

    /** The fact groups behind the read-only sections (as last computed). */
    List<Group> groups() {
        return _groups;
    }

    // ---- test hooks ----
    JComponent fieldForTest( final String key ) {
        return _fields.get( key );
    }

    void setTextForTest( final String key, final String value ) {
        final JComponent c = _fields.get( key );
        if ( c instanceof JTextComponent ) {
            ( (JTextComponent) c ).setText( value );
        }
    }

    boolean hasSectionForTest( final String title ) {
        return _sections.containsKey( title );
    }

    boolean isSectionExpandedForTest( final String title ) {
        final Section s = _sections.get( title );
        return ( s != null ) && s.isExpanded();
    }

    void toggleSectionForTest( final String title ) {
        _sections.get( title ).toggle();
    }

    boolean isOutlinedForTest( final String key ) {
        final JComponent c = _fields.get( key );
        return ( c != null ) && "error".equals( c.getClientProperty( "JComponent.outline" ) );
    }

    JScrollPane scrollPaneForTest() {
        return _scroll;
    }

    String headerTitleForTest() {
        return _header.getTitle();
    }

    String headerSubtitleForTest() {
        return _header.getSubtitle();
    }

    // ------------------------------------------------------------------ building
    private void buildEditableSections( final TreePropertiesDraft d ) {
        {
            final Grid g = new Grid( _label_width );
            final JTextField name = editField( d.name, "a short title for the tree" );
            register( TreePropertiesDraft.NAME, name );
            g.row( "Name", name, false );
            final JTextArea desc = new JTextArea( d.description, 4, 20 );
            desc.setLineWrap( true );
            desc.setWrapStyleWord( true );
            desc.setToolTipText( "Free text; tools append a provenance sentence here when they change the tree" );
            desc.getDocument().addDocumentListener( _doc_listener );
            register( TreePropertiesDraft.DESCRIPTION, desc );
            g.row( "Description", FormWidgets.areaScrollPane( desc, true ), true );
            addSection( SEC_NAME, null, g, true );
        }
        {
            final Grid g = new Grid( _label_width );
            final JTextField id = editField( d.idValue, "e.g. 12345" );
            register( TreePropertiesDraft.ID_VALUE, id );
            final JTextField provider = editField( d.idProvider, "e.g. treebase" );
            register( TreePropertiesDraft.ID_PROVIDER, provider );
            g.row( "Identifier", id, "Provider", provider );
            final JTextField type = editField( d.type, "e.g. gene tree, species tree" );
            register( TreePropertiesDraft.TYPE, type );
            final JTextField unit = editField( d.distanceUnit, "e.g. substitutions/site, Ma" );
            register( TreePropertiesDraft.DISTANCE_UNIT, unit );
            g.row( "Type", type, "Branch-length unit", unit );
            final boolean any = !d.idValue.isEmpty() || !d.idProvider.isEmpty() || !d.type.isEmpty()
                    || !d.distanceUnit.isEmpty();
            addSection( SEC_IDENTITY, null, g, any );
        }
    }

    private void addSection( final String title, final String detail, final JComponent body, final boolean expanded ) {
        final Section s = new Section( title, detail, body, expanded );
        _sections.put( title, s );
        _page.add( s );
    }

    /** Drops and rebuilds the read-only sections from fresh {@link TreeFacts}, keeping each title's expanded state. */
    private void rebuildFactSections() {
        final Map<String, Boolean> expanded = new HashMap<>();
        for( final Section s : _fact_sections ) {
            expanded.put( s.getTitle(), s.isExpanded() );
            _sections.remove( s.getTitle() );
            _page.remove( s );
        }
        _fact_sections.clear();
        final File file = ( _tree_panel != null ) ? _tree_panel.getTreeFile() : null;
        final boolean edited = ( _tree_panel != null ) && _tree_panel.isEdited();
        _groups = TreeFacts.compute( _phylogeny, file, edited, timeAxis() );
        for( final Group group : _groups ) {
            final Grid g = new Grid( _label_width );
            for( final Fact f : group.facts ) {
                final JTextComponent v = viewValue( f.value, f.value.length() > 60 );
                g.row( f.key, v, f.value.length() > 60 );
            }
            if ( group.histogram != null ) {
                g.span( new HistogramPanel( group.histogram ) );
            }
            final Boolean was = expanded.get( group.title );
            final Section s = new Section( group.title, group.detail, g, ( was != null ) ? was : true );
            _sections.put( group.title, s );
            _fact_sections.add( s );
            _page.add( s );
        }
        _page.revalidate();
        _page.repaint();
    }

    /** What the panel knows about the displayed tree's time axis; null without a panel. */
    private TreeFacts.TimeAxis timeAxis() {
        if ( ( _tree_panel == null ) || ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return null;
        }
        final boolean dated = AptxUtil.detectTimeTree( _phylogeny ) == AptxUtil.TIME_TREE_KIND.DATED;
        return new TreeFacts.TimeAxis( _tree_panel.effectiveTimeAxisType(), dated,
                                       dated ? AptxUtil.timeTreeUnit( _phylogeny ) : null,
                                       _tree_panel.timeAxisRootAgeMa(), _tree_panel.timeAxisPresentDate() );
    }

    // ------------------------------------------------------------------ plumbing
    private void register( final String key, final JComponent c ) {
        _fields.put( key, c );
        if ( c instanceof JTextComponent ) {
            ( (JTextComponent) c ).getDocument().addDocumentListener( _doc_listener );
        }
    }

    private String text( final String key ) {
        final JComponent c = _fields.get( key );
        return ( c instanceof JTextComponent ) ? ( (JTextComponent) c ).getText() : "";
    }

    /** Puts {@code d}'s values into the widgets without firing per-keystroke change events. */
    private void load( final TreePropertiesDraft d ) {
        _loading = true;
        try {
            setText( TreePropertiesDraft.NAME, d.name );
            setText( TreePropertiesDraft.DESCRIPTION, d.description );
            setText( TreePropertiesDraft.ID_VALUE, d.idValue );
            setText( TreePropertiesDraft.ID_PROVIDER, d.idProvider );
            setText( TreePropertiesDraft.TYPE, d.type );
            setText( TreePropertiesDraft.DISTANCE_UNIT, d.distanceUnit );
        }
        finally {
            _loading = false;
        }
        _current = null;
        _current_problems = null;
    }

    private void setText( final String key, final String value ) {
        final JComponent c = _fields.get( key );
        if ( ( c instanceof JTextComponent ) && !value.equals( ( (JTextComponent) c ).getText() ) ) {
            ( (JTextComponent) c ).setText( value );
            if ( c instanceof JTextComponent ) {
                ( (JTextComponent) c ).setCaretPosition( 0 );
            }
        }
    }

    private void fireChanged() {
        if ( _loading ) {
            return;
        }
        _current = null;
        _current_problems = null;
        applyProblems( problems() );
        for( final Runnable r : _change_listeners ) {
            r.run();
        }
    }

    private void applyProblems( final List<Problem> problems ) {
        for( final JComponent c : _outlined ) {
            c.putClientProperty( "JComponent.outline", null );
        }
        _outlined.clear();
        for( final Problem p : problems ) {
            final JComponent c = _fields.get( p.key );
            if ( c != null ) {
                c.putClientProperty( "JComponent.outline", "error" );
                _outlined.add( c );
            }
        }
    }
}
