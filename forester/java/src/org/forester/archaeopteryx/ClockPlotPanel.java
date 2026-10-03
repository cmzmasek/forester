// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.RoundRectangle2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

/**
 * The drawing of the clock plot, and its link to the tree. What is plotted, which trees have a plot and what the
 * line is fitted to are {@link ClockPlot}'s; this is the drawing and the pointer, ported from the "Clock plot" section
 * of Archaeopteryx.js 3.23.0 ({@code layClockPlot}, {@code colorClockPlot}, {@code clockMarksAt},
 * {@code bindClockPlotPointer}):
 * <ul>
 * <li>a point wears its node's colour in the tree: the selection and the search hits, "Color by", the file's styles
 * ({@link Host#colorOf});</li>
 * <li>pointing at one lights its node in the tree, and pointing at a node in the tree rings its point;</li>
 * <li>a click selects or deselects the node, a drag selects the tips in the box -- the tree's own selection;</li>
 * <li>tips drawn on one spot are one dot, and the dot answers for all of them: named together, lit together, selected
 * together.</li>
 * </ul>
 * The readout is painted on the plot, last, never a popup (as the tree's hover card is).
 */
final class ClockPlotPanel extends JComponent {

    private static final long serialVersionUID = 1L;

    /** What the panel needs from the tree it plots. */
    interface Host {

        /** The colour the node wears in the tree, or null where it wears none and the point takes the panel's ink. */
        Color colorOf( PhylogenyNode node, boolean tip );

        /** Whether the node is a search hit or selected: its point is drawn larger, outlined, over its neighbours. */
        boolean isHit( PhylogenyNode node );

        boolean isSelected( PhylogenyNode node );

        /** Selects ({@code add}) or deselects every one of {@code nodes}. */
        void select( List<PhylogenyNode> nodes, boolean add );

        /** Lights these nodes in the tree; empty puts the light out. */
        void light( List<PhylogenyNode> nodes );
    }

    static final int            W            = 440;
    static final int            H            = 300;
    static final int            PAD_TOP      = 10;
    static final int            PAD_RIGHT    = 14;
    static final int            PAD_BOTTOM   = 38;
    static final int            PAD_LEFT     = 58;
    /** How near the pointer has to be to a point. */
    static final int            HIT_PX       = 7;
    /** Points nearer each other than this are one dot on the screen. */
    static final double         SAME_PX      = 1;
    /** How many of a dot's tips the readout names. */
    static final int            NAMES_MAX    = 5;
    /** A press that moves further than this is a box, not a click. */
    static final int            DRAG_PX      = 4;
    private static final String MONTHS[]     = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct",
            "Nov", "Dec" };
    private static final double TIP_R        = 3;
    private static final double INTERNAL_R   = 2;
    private static final double HIT_R        = 4;
    private static final double RING_R       = 6.5;

    /** One point as laid down. */
    static final class Mark {

        final ClockPlot.Point _p;
        final double          _px;
        final double          _py;

        Mark( final ClockPlot.Point p, final double px, final double py ) {
            _p = p;
            _px = px;
            _py = py;
        }
    }

    private final Host       _host;
    private ClockPlot.Data   _data;
    private boolean          _calendar;
    private boolean          _show_line     = true;
    private boolean          _show_internal = false;
    private List<Mark>       _marks         = Collections.emptyList();
    // the scales: x(date) = _x0 + (date - _d0) * _xk ; y(div) = _y0 - (div - _v0) * _yk
    private double           _dom_x0, _dom_x1, _dom_y0, _dom_y1;
    private List<Mark>       _hover         = null;
    /** The mark of the node pointed at IN THE TREE, ringed; null for none. */
    private Mark             _ringed_from_tree = null;
    private int              _hover_x, _hover_y;
    private int              _press_x       = -1, _press_y = -1;
    private boolean          _boxed         = false;
    private int              _box_x, _box_y;
    private int              _lays          = 0;

    ClockPlotPanel( final Host host ) {
        _host = host;
        setPreferredSize( new Dimension( W, H ) );
        setOpaque( false );
        final MouseAdapter m = new MouseAdapter() {

            @Override
            public void mouseMoved( final MouseEvent e ) {
                pointAt( e.getX(), e.getY() );
            }

            @Override
            public void mousePressed( final MouseEvent e ) {
                if ( SwingUtilities.isLeftMouseButton( e ) ) {
                    press( e.getX(), e.getY() );
                }
            }

            @Override
            public void mouseDragged( final MouseEvent e ) {
                drag( e.getX(), e.getY() );
            }

            @Override
            public void mouseReleased( final MouseEvent e ) {
                if ( SwingUtilities.isLeftMouseButton( e ) ) {
                    release( e.getX(), e.getY() );
                }
            }

            @Override
            public void mouseExited( final MouseEvent e ) {
                if ( _press_x < 0 ) {
                    leave();
                }
            }
        };
        addMouseListener( m );
        addMouseMotionListener( m );
    }

    // ---- what is shown ---------------------------------------------------------------------------------------------

    /** Lays the plot down for {@code data}. */
    void setData( final ClockPlot.Data data, final boolean calendar ) {
        _data = data;
        _calendar = calendar && ( data != null ) && data._forward;
        lay();
    }

    ClockPlot.Data getData() {
        return _data;
    }

    boolean isCalendar() {
        return _calendar;
    }

    void setShowLine( final boolean show ) {
        _show_line = show;
        repaint();
    }

    boolean isShowLine() {
        return _show_line;
    }

    void setShowInternal( final boolean show ) {
        _show_internal = show;
        lay();
    }

    boolean isShowInternal() {
        return _show_internal;
    }

    private void lay() {
        ++_lays;
        _hover = null; // every point may have moved from under the pointer: known at its next move
        _ringed_from_tree = null;
        if ( _data == null ) {
            _marks = Collections.emptyList();
            repaint();
            return;
        }
        double x0 = Double.POSITIVE_INFINITY, x1 = Double.NEGATIVE_INFINITY;
        double y0 = Double.POSITIVE_INFINITY, y1 = Double.NEGATIVE_INFINITY;
        final List<ClockPlot.Point> shown = new ArrayList<ClockPlot.Point>();
        for( final ClockPlot.Point p : _data._points ) {
            if ( p._tip || _show_internal ) {
                shown.add( p );
                x0 = Math.min( x0, p._date );
                x1 = Math.max( x1, p._date );
                y0 = Math.min( y0, p._div );
                y1 = Math.max( y1, p._div );
            }
        }
        double pad_x = ( x1 - x0 ) * 0.04;
        if ( pad_x == 0 || !Double.isFinite( pad_x ) ) {
            pad_x = 0.5;
        }
        double pad_y = ( y1 - y0 ) * 0.06;
        if ( pad_y == 0 || !Double.isFinite( pad_y ) ) {
            pad_y = Math.abs( y1 ) * 0.06;
        }
        if ( pad_y == 0 || !Double.isFinite( pad_y ) ) {
            pad_y = 1e-6;
        }
        _dom_x0 = x0 - pad_x;
        _dom_x1 = x1 + pad_x;
        _dom_y0 = y0 - pad_y;
        _dom_y1 = y1 + pad_y;
        final List<Mark> marks = new ArrayList<Mark>( shown.size() );
        for( final ClockPlot.Point p : shown ) {
            marks.add( new Mark( p, x( p._date ), y( p._div ) ) );
        }
        _marks = marks;
        repaint();
    }

    /** Time runs left to right, whichever way its numbers do. */
    double x( final double date ) {
        final double left = PAD_LEFT;
        final double right = W - PAD_RIGHT;
        final double t = ( date - _dom_x0 ) / ( _dom_x1 - _dom_x0 );
        return ( ( _data == null ) || _data._forward ) ? ( left + ( t * ( right - left ) ) )
                : ( right - ( t * ( right - left ) ) );
    }

    double y( final double div ) {
        final double top = PAD_TOP;
        final double bottom = H - PAD_BOTTOM;
        return bottom - ( ( ( div - _dom_y0 ) / ( _dom_y1 - _dom_y0 ) ) * ( bottom - top ) );
    }

    // ---- painting ------------------------------------------------------------------------------------------------

    private static Color ui( final String key, final Color fallback ) {
        final Color c = UIManager.getColor( key );
        return ( c != null ) ? c : fallback;
    }

    private static boolean dark() {
        final Color bg = ui( "Panel.background", Color.WHITE );
        return ( ( bg.getRed() * 299 ) + ( bg.getGreen() * 587 ) + ( bg.getBlue() * 114 ) ) < 128000;
    }

    private static Color mix( final Color a, final Color b, final double t ) {
        return new Color( (int) Math.round( a.getRed() + ( ( b.getRed() - a.getRed() ) * t ) ),
                          (int) Math.round( a.getGreen() + ( ( b.getGreen() - a.getGreen() ) * t ) ),
                          (int) Math.round( a.getBlue() + ( ( b.getBlue() - a.getBlue() ) * t ) ) );
    }

    private static Color alpha( final Color c, final double a ) {
        return new Color( c.getRed(), c.getGreen(), c.getBlue(), (int) Math.round( a * 255 ) );
    }

    Color ink() {
        return ui( "Label.foreground", Color.BLACK );
    }

    Color ground() {
        final Color bg = ui( "Panel.background", Color.WHITE );
        return dark() ? mix( bg, Color.WHITE, 0.05 ) : mix( bg, Color.WHITE, 0.6 );
    }

    private Color muted() {
        return mix( ink(), ui( "Panel.background", Color.WHITE ), 0.4 );
    }

    private Color line() {
        return mix( ink(), ui( "Panel.background", Color.WHITE ), 0.85 );
    }

    private Color lineStrong() {
        return mix( ink(), ui( "Panel.background", Color.WHITE ), 0.6 );
    }

    /** The panel's font, or the look and feel's label font: a bare component has none until it is parented. */
    private Font baseFont() {
        final Font f = getFont();
        if ( f != null ) {
            return f;
        }
        final Font label = UIManager.getFont( "Label.font" );
        return ( label != null ) ? label : new Font( Font.SANS_SERIF, Font.PLAIN, 12 );
    }

    static Color accent() {
        return TreePanel.uiAccentColor();
    }

    @Override
    protected void paintComponent( final Graphics g0 ) {
        final Graphics2D g = (Graphics2D) g0.create();
        try {
            g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
            g.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            g.setRenderingHint( RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE );
            paintPlot( g );
        }
        finally {
            g.dispose();
        }
    }

    void paintPlot( final Graphics2D g ) {
        final int left = PAD_LEFT;
        final int right = W - PAD_RIGHT;
        final int top = PAD_TOP;
        final int bottom = H - PAD_BOTTOM;
        // an opaque ground under the points
        g.setColor( ground() );
        g.fillRect( left, top, right - left, bottom - top );
        if ( _data == null ) {
            return;
        }
        final Font tick_font = baseFont().deriveFont( Font.PLAIN, 10f );
        final Font title_font = baseFont().deriveFont( Font.PLAIN, 10.5f );
        g.setFont( tick_font );
        final FontMetrics fm = g.getFontMetrics();
        final BasicStroke one = new BasicStroke( 1f );
        g.setStroke( one );
        // --- the grid, the ticks and what each axis measures
        for( final Tick t : ticksX() ) {
            final double px = Math.round( x( t._value ) ) + 0.5;
            if ( ( px < left ) || ( px > right ) ) {
                continue;
            }
            g.setColor( line() );
            g.draw( new Line2D.Double( px, top, px, bottom ) );
            g.setColor( muted() );
            g.drawString( t._label, (float) ( px - ( fm.stringWidth( t._label ) / 2.0 ) ), bottom + 13 );
        }
        for( final double v : ClockPlot.ticks( _dom_y0, _dom_y1, 5 ) ) {
            final double py = Math.round( y( v ) ) + 0.5;
            g.setColor( line() );
            g.draw( new Line2D.Double( left, py, right, py ) );
            g.setColor( muted() );
            final String s = ClockPlot.plain( Double.parseDouble( String.format( java.util.Locale.ROOT, "%.12g", v ) ) );
            g.drawString( s, left - 6 - fm.stringWidth( s ), (float) ( py + 3.5 ) );
        }
        g.setColor( lineStrong() );
        g.draw( new Line2D.Double( left, bottom + 0.5, right, bottom + 0.5 ) );
        g.draw( new Line2D.Double( left + 0.5, top, left + 0.5, bottom ) );
        g.setFont( title_font );
        final FontMetrics tfm = g.getFontMetrics();
        g.setColor( ink() );
        final String x_title = xTitle();
        g.drawString( x_title, (float) ( ( ( left + right ) / 2.0 ) - ( tfm.stringWidth( x_title ) / 2.0 ) ), H - 6 );
        final String y_title = "Divergence (subs/site)";
        final AffineTransform saved = g.getTransform();
        g.translate( 12 + tfm.getAscent() - 9, ( top + bottom ) / 2.0 );
        g.rotate( -Math.PI / 2 );
        g.drawString( y_title, (float) ( -tfm.stringWidth( y_title ) / 2.0 ), 0 );
        g.setTransform( saved );
        // --- inside the frame
        final Shape clip = g.getClip();
        g.clipRect( left, top, right - left, bottom - top );
        // the ancestors under the line, the tips over it, the hits over everything
        final List<Mark> hits = new ArrayList<Mark>();
        for( final Mark m : _marks ) {
            if ( !m._p._tip ) {
                if ( _host.isHit( m._p._node ) ) {
                    hits.add( m );
                }
                else {
                    paintMark( g, m, false );
                }
            }
        }
        final ClockPlot.Fit fit = _data._fit;
        if ( ( fit != null ) && _show_line ) {
            g.setColor( accent() );
            g.setStroke( new BasicStroke( 1.6f ) );
            g.draw( new Line2D.Double( x( _dom_x0 ), y( fit.at( _dom_x0 ) ), x( _dom_x1 ), y( fit.at( _dom_x1 ) ) ) );
            g.setStroke( one );
        }
        for( final Mark m : _marks ) {
            if ( m._p._tip ) {
                if ( _host.isHit( m._p._node ) ) {
                    hits.add( m );
                }
                else {
                    paintMark( g, m, false );
                }
            }
        }
        for( final Mark m : hits ) {
            paintMark( g, m, true );
        }
        g.setClip( clip );
        // --- the ring, the box, the readout
        final Mark ring = ( _hover != null ) ? _hover.get( 0 ) : _ringed_from_tree;
        if ( ring != null ) {
            g.setColor( accent() );
            g.setStroke( new BasicStroke( 1.8f ) );
            g.draw( new Ellipse2D.Double( ring._px - RING_R, ring._py - RING_R, 2 * RING_R, 2 * RING_R ) );
            g.setStroke( one );
        }
        if ( _boxed ) {
            final Rectangle2D box = box( _box_x, _box_y );
            g.setColor( alpha( accent(), 0.16 ) );
            g.fill( box );
            g.setColor( accent() );
            g.draw( box );
        }
        if ( ( _hover != null ) && !_boxed ) {
            paintReadout( g, readout( _hover ) );
        }
    }

    private void paintMark( final Graphics2D g, final Mark m, final boolean hit ) {
        final Color own = _host.colorOf( m._p._node, m._p._tip );
        final double r = hit ? HIT_R : ( m._p._tip ? TIP_R : INTERNAL_R );
        final Ellipse2D dot = new Ellipse2D.Double( m._px - r, m._py - r, 2 * r, 2 * r );
        if ( own != null ) {
            g.setColor( hit ? own : alpha( own, 0.8 ) );
        }
        else if ( m._p._tip ) {
            g.setColor( alpha( ink(), hit ? 1 : 0.8 ) );
        }
        else {
            g.setColor( alpha( muted(), hit ? 1 : 0.55 ) );
        }
        g.fill( dot );
        if ( hit ) {
            g.setColor( ( own != null ) ? own.darker() : ink() );
            g.setStroke( new BasicStroke( 1.2f ) );
            g.draw( dot );
            g.setStroke( new BasicStroke( 1f ) );
        }
    }

    String xTitle() {
        final String unit = ( !ForesterUtil.isEmpty( _data._unit ) && !_calendar ) ? " (" + _data._unit + ")" : "";
        return ( _data._forward ? "Date" : "Age" ) + unit;
    }

    static final class Tick {

        final double _value;
        final String _label;

        Tick( final double value, final String label ) {
            _value = value;
            _label = label;
        }
    }

    /** The date axis's ticks: whole years over three years or more, months over less, else round numbers. */
    List<Tick> ticksX() {
        final double from = Math.min( _dom_x0, _dom_x1 );
        final double to = Math.max( _dom_x0, _dom_x1 );
        final List<Tick> ticks = new ArrayList<Tick>();
        if ( _calendar ) {
            if ( ( to - from ) >= 3 ) {
                for( final Integer y : ClockPlot.calendarTickYears( from, to ) ) {
                    ticks.add( new Tick( y.doubleValue(), String.valueOf( y ) ) );
                }
                return ticks;
            }
            final List<ClockPlot.MonthTick> months = ClockPlot.calendarTickMonths( from, to );
            if ( months.size() >= 2 ) {
                for( final ClockPlot.MonthTick m : months ) {
                    ticks.add( new Tick( m._value, MONTHS[ m._month - 1 ] + " " + m._year ) );
                }
                return ticks;
            }
        }
        for( final double v : ClockPlot.ticks( from, to, 6 ) ) {
            ticks.add( new Tick( v, ClockPlot.plain( Double.parseDouble( String.format( java.util.Locale.ROOT,
                                                                                        "%.12g",
                                                                                        v ) ) ) ) );
        }
        return ticks;
    }

    // ---- the readout ---------------------------------------------------------------------------------------------

    /** What the readout says of the points under the pointer: of one, what it is; of several on one spot, how many,
     *  the names of the first few, and what they share. */
    static List<String> readout( final List<Mark> marks ) {
        final List<String> lines = new ArrayList<String>();
        final ClockPlot.Point p = marks.get( 0 )._p;
        if ( marks.size() == 1 ) {
            lines.add( !ForesterUtil.isEmpty( p._node.getName() ) ? "Name: " + p._node.getName()
                    : ( p._tip ? "Tip" : "Internal node" ) );
            final String date = dateText( p._node );
            if ( date != null ) {
                lines.add( "Date: " + date );
            }
            lines.add( "Divergence: " + ClockPlot.number( p._div, 4 ) );
            if ( p._residual != null ) {
                lines.add( "Off the line: " + ClockPlot.signed( p._residual.doubleValue() ) );
            }
            if ( !p._tip ) {
                lines.add( "Tips below: " + p._node.getNumberOfExternalNodes() );
            }
            return lines;
        }
        lines.add( marks.size() + ( p._tip ? " tips here" : " internal nodes here" ) );
        final List<Mark> named = new ArrayList<Mark>();
        for( final Mark m : marks ) {
            if ( !ForesterUtil.isEmpty( m._p._node.getName() ) ) {
                named.add( m );
            }
        }
        for( int i = 0; i < Math.min( NAMES_MAX, named.size() ); ++i ) {
            lines.add( "Name: " + named.get( i )._p._node.getName() );
        }
        if ( named.size() > NAMES_MAX ) {
            lines.add( "Name: … and " + ( named.size() - NAMES_MAX ) + " more" );
        }
        final String unit = ( p._node.getNodeData().isHasDate()
                && !ForesterUtil.isEmpty( p._node.getNodeData().getDate().getUnit() ) )
                        ? " " + p._node.getNodeData().getDate().getUnit() : "";
        final double[] dates = new double[ marks.size() ];
        final double[] divs = new double[ marks.size() ];
        boolean all_residuals = true;
        final double[] res = new double[ marks.size() ];
        for( int i = 0; i < marks.size(); ++i ) {
            dates[ i ] = marks.get( i )._p._date;
            divs[ i ] = marks.get( i )._p._div;
            if ( marks.get( i )._p._residual == null ) {
                all_residuals = false;
            }
            else {
                res[ i ] = marks.get( i )._p._residual.doubleValue();
            }
        }
        lines.add( "Date: " + span( dates, 0 ) + unit );
        lines.add( "Divergence: " + span( divs, 1 ) );
        if ( all_residuals ) {
            lines.add( "Off the line: " + span( res, 2 ) );
        }
        return lines;
    }

    /** One value where the points agree, "least – greatest" where they do not (points a fraction of a pixel apart
     *  are one dot and may differ a hair). {@code how}: 0 a date as stated, 1 a divergence, 2 a residual. */
    private static String span( final double[] values, final int how ) {
        double lo = Double.POSITIVE_INFINITY, hi = Double.NEGATIVE_INFINITY;
        for( final double v : values ) {
            lo = Math.min( lo, v );
            hi = Math.max( hi, v );
        }
        return ( lo == hi ) ? format( lo, how ) : ( format( lo, how ) + " – " + format( hi, how ) );
    }

    private static String format( final double v, final int how ) {
        if ( how == 0 ) {
            return ClockPlot.plain( v );
        }
        return ( how == 1 ) ? ClockPlot.number( v, 4 ) : ClockPlot.signed( v );
    }

    /** A node's date as the node window states it: value, a genuine interval, its unit, its description. */
    private static String dateText( final PhylogenyNode n ) {
        if ( !n.getNodeData().isHasDate() ) {
            return null;
        }
        final org.forester.phylogeny.data.Date d = n.getNodeData().getDate();
        String s = ( d.getValue() != null ) ? d.getValue().toPlainString() : "";
        if ( ( d.getMin() != null ) && ( d.getMax() != null ) && AptxUtil.hasDateIntervalWidth( d ) ) {
            s += " [" + d.getMin().toPlainString() + " - " + d.getMax().toPlainString() + "]";
        }
        if ( !s.isEmpty() && !ForesterUtil.isEmpty( d.getUnit() ) ) {
            s += " " + d.getUnit();
        }
        s = s.trim();
        if ( !ForesterUtil.isEmpty( d.getDesc() ) ) {
            s = s.isEmpty() ? d.getDesc() : s + " (" + d.getDesc() + ")";
        }
        return s.isEmpty() ? null : s;
    }

    private void paintReadout( final Graphics2D g, final List<String> lines ) {
        g.setFont( baseFont().deriveFont( Font.PLAIN, 11f ) );
        final FontMetrics fm = g.getFontMetrics();
        int w = 0;
        for( final String s : lines ) {
            w = Math.max( w, fm.stringWidth( s ) );
        }
        final int pad = 6;
        final int lh = fm.getHeight();
        final int bw = w + ( 2 * pad );
        final int bh = ( lines.size() * lh ) + ( 2 * pad ) - fm.getLeading();
        int bx = _hover_x + 14;
        int by = _hover_y + 14;
        if ( ( bx + bw ) > ( W - 2 ) ) {
            bx = _hover_x - 14 - bw;
        }
        if ( ( by + bh ) > ( H - 2 ) ) {
            by = _hover_y - 14 - bh;
        }
        bx = Math.max( 2, Math.min( bx, W - 2 - bw ) );
        by = Math.max( 2, Math.min( by, H - 2 - bh ) );
        final RoundRectangle2D card = new RoundRectangle2D.Double( bx, by, bw, bh, 8, 8 );
        g.setColor( alpha( ui( "Panel.background", Color.WHITE ), 0.95 ) );
        g.fill( card );
        g.setColor( lineStrong() );
        g.draw( card );
        g.setColor( ink() );
        int ty = by + pad + fm.getAscent();
        for( final String s : lines ) {
            g.drawString( s, bx + pad, ty );
            ty += lh;
        }
    }

    // ---- the pointer ---------------------------------------------------------------------------------------------

    /**
     * The points under the pointer: the nearest one within reach, and every other of its kind drawn on the same spot.
     * Tips with one date and one divergence -- identical sequences sampled on one day -- are one dot on the screen, and
     * the dot stands for all of them (Christian, 2026-10-02, in Archaeopteryx.js). The SAME SPOT, not everything within
     * reach: in a dense cloud that is hundreds of tips, and a click meant for one dot would take them all.
     *
     * @return the marks in the tree's own order (preorder); empty where the pointer is on none
     */
    List<Mark> marksAt( final double px, final double py ) {
        Mark best = null;
        double best_d = HIT_PX * HIT_PX;
        for( final Mark m : _marks ) {
            final double dx = m._px - px;
            final double dy = m._py - py;
            final double d = ( dx * dx ) + ( dy * dy );
            // a tip wins over an ancestor at the same distance: it is the datum
            if ( ( d < best_d ) || ( ( d == best_d ) && ( best != null ) && m._p._tip && !best._p._tip ) ) {
                best = m;
                best_d = d;
            }
        }
        if ( best == null ) {
            return Collections.emptyList();
        }
        final List<Mark> same = new ArrayList<Mark>();
        for( final Mark m : _marks ) {
            if ( ( m._p._tip == best._p._tip ) && ( Math.abs( m._px - best._px ) < SAME_PX )
                    && ( Math.abs( m._py - best._py ) < SAME_PX ) ) {
                same.add( m );
            }
        }
        return same;
    }

    private static List<PhylogenyNode> nodesOf( final List<Mark> marks ) {
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>( marks.size() );
        for( final Mark m : marks ) {
            nodes.add( m._p._node );
        }
        return nodes;
    }

    void pointAt( final int px, final int py ) {
        final List<Mark> found = marksAt( px, py );
        final List<Mark> marks = found.isEmpty() ? null : found;
        _hover_x = px;
        _hover_y = py;
        if ( !same( marks, _hover ) ) {
            _hover = marks;
            _host.light( ( marks == null ) ? Collections.<PhylogenyNode> emptyList() : nodesOf( marks ) );
            setCursor( ( marks != null ) ? Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) : Cursor.getDefaultCursor() );
        }
        repaint();
    }

    private static boolean same( final List<Mark> a, final List<Mark> b ) {
        if ( ( a == null ) || ( b == null ) ) {
            return a == b;
        }
        return a.equals( b );
    }

    void press( final int px, final int py ) {
        _press_x = px;
        _press_y = py;
        _boxed = false;
    }

    void drag( final int px, final int py ) {
        if ( _press_x < 0 ) {
            return;
        }
        if ( !_boxed && ( Math.hypot( px - _press_x, py - _press_y ) > DRAG_PX ) ) {
            _boxed = true;
            _hover = null;
            _host.light( Collections.<PhylogenyNode> emptyList() );
        }
        _box_x = px;
        _box_y = py;
        repaint();
    }

    private Rectangle2D box( final int px, final int py ) {
        return new Rectangle2D.Double( Math.min( _press_x, px ),
                                       Math.min( _press_y, py ),
                                       Math.abs( px - _press_x ),
                                       Math.abs( py - _press_y ) );
    }

    void release( final int px, final int py ) {
        if ( _press_x < 0 ) {
            return;
        }
        final boolean boxed = _boxed;
        final Rectangle2D box = box( px, py );
        _press_x = -1;
        _press_y = -1;
        _boxed = false;
        repaint();
        if ( boxed ) {
            // the tips in the box join the selection; an ancestor is selected by a click, where it is one node and
            // meant
            final List<PhylogenyNode> add = new ArrayList<PhylogenyNode>();
            for( final Mark m : _marks ) {
                if ( m._p._tip && ( m._px >= box.getMinX() ) && ( m._px <= box.getMaxX() ) && ( m._py >= box.getMinY() )
                        && ( m._py <= box.getMaxY() ) && !_host.isSelected( m._p._node ) ) {
                    add.add( m._p._node );
                }
            }
            if ( !add.isEmpty() ) {
                _host.select( add, true );
            }
            return;
        }
        // a click selects every node of the dot; where all of them are selected already, it deselects them
        final List<Mark> marks = marksAt( px, py );
        if ( !marks.isEmpty() ) {
            boolean all = true;
            for( final Mark m : marks ) {
                if ( !_host.isSelected( m._p._node ) ) {
                    all = false;
                    break;
                }
            }
            _host.select( nodesOf( marks ), !all );
        }
    }

    void leave() {
        if ( _hover != null ) {
            _hover = null;
            _host.light( Collections.<PhylogenyNode> emptyList() );
        }
        setCursor( Cursor.getDefaultCursor() );
        repaint();
    }

    /** A node pointed at IN THE TREE rings its point, where it has one; null takes the ring away. */
    void markNode( final PhylogenyNode node ) {
        Mark mark = null;
        if ( node != null ) {
            for( final Mark m : _marks ) {
                if ( m._p._node == node ) {
                    mark = m;
                    break;
                }
            }
        }
        if ( mark != _ringed_from_tree ) {
            _ringed_from_tree = mark;
            repaint();
        }
    }

    /** The points under the pointer, for the tree's paint: lit again after a render (a click redraws the tree). */
    List<PhylogenyNode> hoveredNodes() {
        return ( _hover == null ) ? Collections.<PhylogenyNode> emptyList() : nodesOf( _hover );
    }

    // ---- test seams ----------------------------------------------------------------------------------------------

    List<Mark> marksForTest() {
        return _marks;
    }

    Mark ringedForTest() {
        return ( _hover != null ) ? _hover.get( 0 ) : _ringed_from_tree;
    }

    int laysForTest() {
        return _lays;
    }
}
