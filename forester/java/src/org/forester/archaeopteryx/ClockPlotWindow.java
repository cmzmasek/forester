// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.UIManager;
import javax.swing.WindowConstants;

import org.forester.phylogeny.PhylogenyNode;

/**
 * The CLOCK PLOT window: modeless, linked to the tree of the current tab (Christian, 2026-10-02: "desktop should copy
 * what JS did for the clock plot"; Archaeopteryx.js 3.23.0 {@code showClockPlot}). The drawing is
 * {@link ClockPlotPanel}; under it the two checkboxes, Deselect all, what the line says (Rate, the root date by the
 * line over the one in the tree, R², Tips) and a note.
 * <p>
 * The plot follows the tree on view. {@link #treePainted} is called at the end of every screen paint of the current
 * tab, so it costs next to nothing when nothing changed: the plot is worked out again only for another tab, another
 * tree on display (a subtree view, an undo), a different number of nodes (one deleted), or a change the Time | Div
 * offer is asked again for (a date, property or length written in the node editor); otherwise it only asks whether
 * any point's colour changed. A tab whose tree has no plot closes the window.
 */
final class ClockPlotWindow extends JDialog {

    private static final long  serialVersionUID = 1L;
    static final String        TITLE            = "Clock plot";
    static final String        NOTE_LINE        = "The line is fitted to the tips only. Tips share ancestry, so R² describes the fit and is no test. Click a point to select its node; drag to select the tips in a box.";
    static final String        NOTE_RATES       = "Divergence here is each branch’s time × its clock rate: the points show the model’s rates, not a measurement of their own. ";
    /** Tips only to begin with (Christian, 2026-10-02, in Archaeopteryx.js): the line is theirs, and a deep root crowds
     *  them into a corner. Kept while the program runs; Reset to Defaults puts them back. */
    static final boolean       DEFAULT_LINE     = true;
    static final boolean       DEFAULT_INTERNAL = false;
    private static boolean     _line_option     = DEFAULT_LINE;
    private static boolean     _internal_option = DEFAULT_INTERNAL;

    private final MainFrame    _frame;
    private final ClockPlotPanel _plot;
    private final JCheckBox    _line_cb;
    private final JCheckBox    _internal_cb;
    private final JButton      _deselect;
    private final JPanel       _stats;
    private final JTextArea    _note;
    private TreePanel          _tp;
    private boolean            _closed;
    // what the plot was laid down for
    private Object             _seen_tree;
    private int                _seen_count = -1;
    private long               _seen_epoch = -1;
    private int[]              _seen_colors;
    private List<Object>       _seen_colour_key;
    private int                _colour_passes;
    private final List<String[]> _stat_rows = new ArrayList<String[]>();

    ClockPlotWindow( final MainFrame frame ) {
        super( frame, TITLE, false );
        _frame = frame;
        setDefaultCloseOperation( WindowConstants.DISPOSE_ON_CLOSE );
        _plot = new ClockPlotPanel( new ClockPlotPanel.Host() {

            @Override
            public Color colorOf( final PhylogenyNode node, final boolean tip ) {
                return ( _tp == null ) ? null : _tp.clockPlotColorOf( node, tip );
            }

            @Override
            public boolean isHit( final PhylogenyNode node ) {
                return ( _tp != null ) && _tp.isClockPlotHit( node );
            }

            @Override
            public boolean isSelected( final PhylogenyNode node ) {
                return ( _tp != null ) && _tp.isClockPlotSelected( node );
            }

            @Override
            public void select( final List<PhylogenyNode> nodes, final boolean add ) {
                if ( _tp != null ) {
                    _tp.selectFromClockPlot( nodes, add );
                }
            }

            @Override
            public void light( final List<PhylogenyNode> nodes ) {
                if ( _tp != null ) {
                    _tp.setClockPlotLit( nodes );
                }
            }
        } );
        _plot.setShowLine( _line_option );
        _plot.setShowInternal( _internal_option );
        final JPanel body = new JPanel( new BorderLayout( 0, 6 ) );
        body.setBorder( BorderFactory.createEmptyBorder( 10, 12, 10, 12 ) );
        body.add( _plot, BorderLayout.NORTH );
        final JPanel below = new JPanel();
        below.setLayout( new BoxLayout( below, BoxLayout.Y_AXIS ) );
        final JPanel options = new JPanel();
        options.setLayout( new BoxLayout( options, BoxLayout.X_AXIS ) );
        _line_cb = new JCheckBox( "Regression line", _line_option );
        _line_cb.setToolTipText( "to show/hide the least-squares line through the tips" );
        _line_cb.addActionListener( e -> {
            _line_option = _line_cb.isSelected();
            _plot.setShowLine( _line_option );
        } );
        _internal_cb = new JCheckBox( "Internal nodes", _internal_option );
        _internal_cb.setToolTipText( "to show/hide the internal nodes: their dates were inferred, and they take no part in the line" );
        _internal_cb.addActionListener( e -> {
            _internal_option = _internal_cb.isSelected();
            _plot.setShowInternal( _internal_option );
            _seen_colors = null;
        } );
        _deselect = new JButton( "Deselect all" );
        _deselect.setToolTipText( "deselect every selected node of the tree" );
        _deselect.addActionListener( e -> {
            if ( _tp != null ) {
                _tp.deselectAllFromClockPlot();
            }
        } );
        options.add( _line_cb );
        options.add( Box.createHorizontalStrut( 14 ) );
        options.add( _internal_cb );
        options.add( Box.createHorizontalGlue() );
        options.add( _deselect );
        options.setAlignmentX( LEFT_ALIGNMENT );
        below.add( options );
        below.add( Box.createVerticalStrut( 6 ) );
        _stats = new JPanel( new GridBagLayout() );
        _stats.setAlignmentX( LEFT_ALIGNMENT );
        below.add( _stats );
        below.add( Box.createVerticalStrut( 8 ) );
        _note = new JTextArea();
        _note.setEditable( false );
        _note.setFocusable( false );
        _note.setLineWrap( true );
        _note.setWrapStyleWord( true );
        _note.setOpaque( false );
        _note.setBorder( null );
        _note.setAlignmentX( LEFT_ALIGNMENT );
        below.add( _note );
        body.add( below, BorderLayout.CENTER );
        setContentPane( body );
    }

    /** Closes NOW, whoever closes it (the title bar's close, the button, a tab with no plot): nothing waits for a
     *  {@code windowClosed} event, which arrives later (Archaeopteryx.js 3.23.0 learned this on the live demo: the
     *  button "closed" a window already closed). */
    @Override
    public void dispose() {
        if ( !_closed ) {
            _closed = true;
            closed();
        }
        super.dispose();
    }

    /** Opens the window for {@code tp}'s tree; false (and nothing opened) where it has no plot. */
    boolean showFor( final TreePanel tp ) {
        if ( !rebind( tp ) ) {
            return false;
        }
        pack();
        place();
        setVisible( true );
        return true;
    }

    /** In the tree area's top right corner, inside the screen. */
    private void place() {
        final Rectangle screen = getGraphicsConfiguration().getBounds();
        int right = screen.x + screen.width;
        int top = screen.y;
        if ( ( _tp != null ) && _tp.isShowing() && ( _tp.getParent() != null ) ) {
            final java.awt.Component area = ( _tp.getParent().getParent() != null ) ? _tp.getParent().getParent()
                    : _tp.getParent();
            final Point p = area.getLocationOnScreen();
            right = Math.min( right, p.x + area.getWidth() );
            top = Math.max( top, p.y );
        }
        final Dimension d = getSize();
        setLocation( Math.max( screen.x + 6, right - d.width - 14 ),
                     Math.max( screen.y + 6, Math.min( top + 14, ( screen.y + screen.height ) - d.height - 6 ) ) );
    }

    /** The tab the window plots, worked out afresh; false where its tree has no plot. */
    private boolean rebind( final TreePanel tp ) {
        if ( ( _tp != null ) && ( _tp != tp ) ) {
            _tp.setClockPlotLit( null );
        }
        _tp = tp;
        _seen_tree = null;
        return refresh( true );
    }

    /**
     * After a screen paint of the CURRENT tab's tree: the plot follows the view, the colours and the selection. Closes
     * the window where the tree has lost its plot.
     */
    void treePainted( final TreePanel tp ) {
        if ( _closed ) {
            return;
        }
        if ( tp != _tp ) {
            if ( !rebind( tp ) ) {
                dispose();
            }
            return;
        }
        if ( !refresh( false ) ) {
            dispose();
        }
    }

    /** False where the tree no longer has a plot. */
    private boolean refresh( final boolean force ) {
        if ( _tp == null ) {
            return false;
        }
        final Object tree = _tp.getPhylogeny();
        final int count = ( _tp.getPhylogeny() == null ) ? 0 : _tp.getPhylogeny().getNodeCount();
        final long epoch = _tp.clockPlotEpoch();
        boolean laid = false;
        if ( force || ( tree != _seen_tree ) || ( count != _seen_count ) || ( epoch != _seen_epoch ) ) {
            final ClockPlot.Data data = _tp.clockPlotData();
            if ( data == null ) {
                return false;
            }
            _seen_tree = tree;
            _seen_count = count;
            _seen_epoch = epoch;
            _seen_colors = null;
            _tp.setClockPlotLit( null );
            _plot.setData( data, _tp.isClockPlotCalendar() );
            updateStats( data );
            laid = true;
        }
        // the colours, worked out again only where what they are read from changed
        final List<Object> key = _tp.clockPlotColourKey();
        if ( laid || !key.equals( _seen_colour_key ) ) {
            _seen_colour_key = key;
            recolour();
        }
        relight();
        _deselect.setEnabled( _tp.hasClockPlotSelection() );
        return true;
    }

    /** One pass over the points' colours; the drawing repaints only where one changed. */
    private void recolour() {
        ++_colour_passes;
        final List<ClockPlotPanel.Mark> marks = _plot.marksForTest();
        final int[] colors = new int[ marks.size() ];
        for( int i = 0; i < colors.length; ++i ) {
            final PhylogenyNode n = marks.get( i )._p._node;
            final Color c = _tp.clockPlotColorOf( n, marks.get( i )._p._tip );
            colors[ i ] = ( ( c == null ) ? 0 : ( c.getRGB() | 1 ) ) ^ ( _tp.isClockPlotHit( n ) ? 0x5A5A5A5A : 0 );
        }
        if ( ( _seen_colors == null ) || !java.util.Arrays.equals( colors, _seen_colors ) ) {
            _seen_colors = colors;
            _plot.repaint();
        }
    }

    /** The tree lets go of its light when its structure changes (a collapse); the pointer is still on its dot, so the
     *  dot's nodes are lit again, as the tree now shows them. Repaints the tree only where that changes anything. */
    private void relight() {
        _tp.setClockPlotLit( _plot.hoveredNodes() );
    }

    /** Whether the window plots {@code tp}'s tree. */
    boolean isFor( final TreePanel tp ) {
        return _tp == tp;
    }

    private void updateStats( final ClockPlot.Data data ) {
        _stat_rows.clear();
        final boolean calendar = _plot.isCalendar();
        final ClockPlot.Fit fit = data._fit;
        final String per = " subs/site per "
                + ( calendar ? "year" : ( ( data._unit != null ) && !data._unit.trim().isEmpty() ? data._unit.trim() : "unit of time" ) );
        if ( fit != null ) {
            _stat_rows.add( new String[] { "Rate", ClockPlot.number( fit._rate, 3 ) + per } );
            _stat_rows.add( new String[] { "Root date, by the line",
                    ( fit._root_date == null ) ? "none: divergence does not rise with time"
                            : ClockPlot.date( fit._root_date.doubleValue(), calendar ) } );
            _stat_rows.add( new String[] { "Root date, in the tree", ClockPlot.date( data._root_date, calendar ) } );
            _stat_rows.add( new String[] { "R²",
                    ( fit._r2 == null ) ? "none: every tip has one divergence" : ClockPlot.number( fit._r2.doubleValue(), 3 ) } );
        }
        else {
            _stat_rows.add( new String[] { "Line", "none: it takes three tips, not all on one date" } );
        }
        _stat_rows.add( new String[] { "Tips",
                String.format( java.util.Locale.US, "%,d", data.tipCount() )
                        + ( _tp.isCurrentTreeIsSubtree() ? " (the clade on view)" : "" ) } );
        _stats.removeAll();
        final GridBagConstraints c = new GridBagConstraints();
        c.gridy = 0;
        c.anchor = GridBagConstraints.WEST;
        c.insets = new Insets( 1, 0, 1, 12 );
        final Color muted = UIManager.getColor( "Label.disabledForeground" );
        for( final String[] row : _stat_rows ) {
            c.gridx = 0;
            c.weightx = 0;
            final JLabel k = new JLabel( row[ 0 ] );
            if ( muted != null ) {
                k.setForeground( muted );
            }
            _stats.add( k, c );
            c.gridx = 1;
            c.weightx = 1;
            final JLabel v = new JLabel( row[ 1 ] );
            v.setFont( v.getFont().deriveFont( Font.BOLD ) );
            _stats.add( v, c );
            ++c.gridy;
        }
        _note.setText( ( data._from_rates ? NOTE_RATES : "" ) + NOTE_LINE );
        _note.setColumns( 0 );
        _note.setSize( new Dimension( ClockPlotPanel.W, Short.MAX_VALUE ) );
        _stats.revalidate();
        if ( isVisible() ) {
            pack();
        }
    }

    /** A node pointed at in the tree rings its point. */
    void markNode( final PhylogenyNode node ) {
        _plot.markNode( node );
    }

    private void closed() {
        if ( _tp != null ) {
            _tp.setClockPlotLit( null );
        }
        _frame.clockPlotWindowClosed( this );
    }

    /** Reset to Defaults: the two checkboxes as they open. */
    static void resetOptionsToDefaults() {
        _line_option = DEFAULT_LINE;
        _internal_option = DEFAULT_INTERNAL;
    }

    void applyOptions() {
        _line_cb.setSelected( _line_option );
        _internal_cb.setSelected( _internal_option );
        _plot.setShowLine( _line_option );
        if ( _plot.isShowInternal() != _internal_option ) {
            _plot.setShowInternal( _internal_option );
            _seen_colors = null;
        }
    }

    // ---- test seams ----------------------------------------------------------------------------------------------

    ClockPlotPanel plotForTest() {
        return _plot;
    }

    int colourPassesForTest() {
        return _colour_passes;
    }

    TreePanel treePanelForTest() {
        return _tp;
    }

    List<String[]> statsForTest() {
        return _stat_rows;
    }

    String noteForTest() {
        return _note.getText();
    }

    JCheckBox lineCheckBoxForTest() {
        return _line_cb;
    }

    JCheckBox internalCheckBoxForTest() {
        return _internal_cb;
    }

    JButton deselectButtonForTest() {
        return _deselect;
    }

    static void setOptionsForTest( final boolean line, final boolean internal ) {
        _line_option = line;
        _internal_option = internal;
    }

    static boolean lineOptionForTest() {
        return _line_option;
    }

    static boolean internalOptionForTest() {
        return _internal_option;
    }

}
