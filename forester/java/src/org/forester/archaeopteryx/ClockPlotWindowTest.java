// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The Clock plot WINDOW and its link to the tree, in a real frame on the demo pair ({@code clock-plot.nex}, and its
 * twin {@code clock-plot-one-date.nex} that has Time | Div but no plot): the button, what the window says, pointing
 * (a dot of four identical samples names and lights all four), the selection both ways (click, box, Deselect all),
 * a node pointed at in the tree ringing its point, the checkboxes, a subtree view, a collapsed clade, a tab without a
 * plot, and closing at once. Headful: headless, every check here is a green no-op.
 */
public final class ClockPlotWindowTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ClockPlotWindowTest: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean line = ClockPlotWindow.lineOptionForTest();
        final boolean internal = ClockPlotWindow.internalOptionForTest();
        try {
            ClockPlotWindow.setOptionsForTest( ClockPlotWindow.DEFAULT_LINE, ClockPlotWindow.DEFAULT_INTERNAL );
            return panelOk() && windowOk();
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return fail( "exception: " + t );
        }
        finally {
            ClockPlotWindow.setOptionsForTest( line, internal );
        }
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [ClockPlotWindowTest] " + m );
        return false;
    }

    private static void fail( final boolean[] ok, final String m ) {
        ok[ 0 ] = false;
        fail( m );
    }

    private static Phylogeny read( final String demo ) throws Exception {
        final Phylogeny[] phys = FigureRenderer.readTrees( new File( System.getProperty( "user.dir" ), "forester/demo/" + demo ) );
        return phys[ 0 ];
    }

    // ---- the panel alone, with a host that records ---------------------------------------------------------------

    private static final class RecordingHost implements ClockPlotPanel.Host {

        final Set<PhylogenyNode>  _selected = new java.util.HashSet<PhylogenyNode>();
        List<PhylogenyNode>       _lit      = new ArrayList<PhylogenyNode>();
        int                       _selects  = 0;

        @Override
        public java.awt.Color colorOf( final PhylogenyNode node, final boolean tip ) {
            return null;
        }

        @Override
        public boolean isHit( final PhylogenyNode node ) {
            return _selected.contains( node );
        }

        @Override
        public boolean isSelected( final PhylogenyNode node ) {
            return _selected.contains( node );
        }

        @Override
        public void select( final List<PhylogenyNode> nodes, final boolean add ) {
            ++_selects;
            if ( add ) {
                _selected.addAll( nodes );
            }
            else {
                _selected.removeAll( nodes );
            }
        }

        @Override
        public void light( final List<PhylogenyNode> nodes ) {
            _lit = new ArrayList<PhylogenyNode>( nodes );
        }
    }

    private static List<ClockPlotPanel.Mark> named( final ClockPlotPanel panel, final String prefix ) {
        final List<ClockPlotPanel.Mark> out = new ArrayList<ClockPlotPanel.Mark>();
        for( final ClockPlotPanel.Mark m : panel.marksForTest() ) {
            if ( ( m._p._node.getName() != null ) && m._p._node.getName().startsWith( prefix ) ) {
                out.add( m );
            }
        }
        return out;
    }

    private static boolean panelOk() throws Exception {
        final Phylogeny phy = read( "clock-plot.nex" );
        final ClockPlot.Data data = ClockPlot.data( phy, BranchLengthLayout.TimeLengths.onScreen( phy ), null );
        final RecordingHost host = new RecordingHost();
        final ClockPlotPanel panel = new ClockPlotPanel( host );
        panel.setSize( ClockPlotPanel.W, ClockPlotPanel.H );
        panel.setData( data, true );
        if ( panel.marksForTest().size() != 20 ) {
            return fail( "tips only to begin with: 20 marks, got " + panel.marksForTest().size() );
        }
        // the four identical samples are ONE dot that stands for all of them, in the tree's order
        final List<ClockPlotPanel.Mark> four = new ArrayList<ClockPlotPanel.Mark>();
        for( final String n : new String[] { "Americas/7/2023", "Americas/8/2023", "Americas/9/2023", "Americas/10/2023" } ) {
            four.addAll( named( panel, n ) );
        }
        final ClockPlotPanel.Mark dot = four.get( 0 );
        final List<ClockPlotPanel.Mark> at = panel.marksAt( dot._px + 2, dot._py - 2 );
        if ( !at.equals( four ) ) {
            return fail( "the dot of four identical samples answers for all four, in tree order: " + at.size() );
        }
        final List<String> lines = ClockPlotPanel.readout( at );
        if ( !lines.get( 0 ).equals( "4 tips here" ) || !lines.get( 1 ).equals( "Name: Americas/7/2023" )
                || !lines.contains( "Date: 2023.25 year" ) || !lines.contains( "Divergence: 0.0128" ) ) {
            return fail( "the readout of the dot: how many, the names, what they share: " + lines );
        }
        // one spot, not everything in reach: a neighbour 3 px off is its own dot
        final ClockPlotPanel.Mark asia6 = named( panel, "Asia/6/2022" ).get( 0 );
        final List<ClockPlotPanel.Mark> one = panel.marksAt( asia6._px, asia6._py );
        if ( ( one.size() != 1 ) || ( one.get( 0 ) != asia6 ) ) {
            return fail( "a lone point is one mark" );
        }
        final List<String> single = ClockPlotPanel.readout( one );
        if ( !single.equals( java.util.Arrays.asList( "Name: Asia/6/2022", "Date: 2022.7 year", "Divergence: 0.0154",
                                                      "Off the line: +0.00365" ) ) ) {
            return fail( "the readout of one tip: " + single );
        }
        if ( !panel.marksAt( ClockPlotPanel.PAD_LEFT + 1, ClockPlotPanel.PAD_TOP + 1 ).isEmpty() ) {
            return fail( "nothing under an empty corner" );
        }
        // pointing lights the dot's nodes; leaving puts the light out
        panel.pointAt( (int) Math.round( dot._px ), (int) Math.round( dot._py ) );
        if ( host._lit.size() != 4 ) {
            return fail( "pointing at the dot lights its four nodes, got " + host._lit.size() );
        }
        panel.leave();
        if ( !host._lit.isEmpty() ) {
            return fail( "leaving puts the light out" );
        }
        // a click selects all four; a second click, all four selected, deselects them; with one of them selected a
        // click selects the rest
        click( panel, dot );
        if ( host._selected.size() != 4 ) {
            return fail( "a click on the dot selects all four, got " + host._selected.size() );
        }
        click( panel, dot );
        if ( !host._selected.isEmpty() ) {
            return fail( "a click on the dot, all four selected, deselects them" );
        }
        host._selected.add( four.get( 2 )._p._node );
        click( panel, dot );
        if ( host._selected.size() != 4 ) {
            return fail( "one of four selected: a click selects them all" );
        }
        host._selected.clear();
        // a press that moves a few pixels is still a click
        panel.press( (int) Math.round( dot._px ), (int) Math.round( dot._py ) );
        panel.drag( (int) Math.round( dot._px ) + 3, (int) Math.round( dot._py ) );
        panel.release( (int) Math.round( dot._px ) + 3, (int) Math.round( dot._py ) );
        if ( host._selected.size() != 4 ) {
            return fail( "a press moved 3 px is a click" );
        }
        host._selected.clear();
        // a box over the right third of the plot selects the tips in it, and only tips, with ancestors shown too
        panel.setShowInternal( true );
        if ( panel.marksForTest().size() != 29 ) {
            return fail( "Internal nodes: 29 marks, got " + panel.marksForTest().size() );
        }
        final int x0 = (int) Math.round( panel.x( 2022.0 ) );
        boxDrag( panel, x0, ClockPlotPanel.PAD_TOP, ClockPlotPanel.W - ClockPlotPanel.PAD_RIGHT, ClockPlotPanel.H - ClockPlotPanel.PAD_BOTTOM );
        int want = 0;
        boolean ancestor_inside = false;
        for( final ClockPlotPanel.Mark m : panel.marksForTest() ) {
            if ( m._px >= x0 ) {
                if ( m._p._tip ) {
                    ++want;
                }
                else {
                    ancestor_inside = true;
                }
            }
        }
        if ( !ancestor_inside ) {
            return fail( "precondition: an ancestor lies inside the box" );
        }
        if ( ( want < 5 ) || ( host._selected.size() != want ) ) {
            return fail( "the box selects its " + want + " tips and no ancestor, got " + host._selected.size() );
        }
        for( final PhylogenyNode n : host._selected ) {
            if ( !n.isExternal() ) {
                return fail( "a box selected an ancestor" );
            }
        }
        panel.setShowInternal( false );
        // a node pointed at in the tree rings its point; one with no point (an ancestor, tips only) rings none
        panel.markNode( asia6._p._node );
        if ( ( panel.ringedForTest() == null ) || ( panel.ringedForTest()._p._node != asia6._p._node ) ) {
            return fail( "a node pointed at in the tree rings its point" );
        }
        panel.markNode( phy.getRoot() );
        if ( panel.ringedForTest() != null ) {
            return fail( "an ancestor with no point shown rings nothing" );
        }
        // the line: the ink under it changes with Regression line
        final int with = accentInk( panel );
        panel.setShowLine( false );
        final int without = accentInk( panel );
        panel.setShowLine( true );
        if ( !( with > ( without + 100 ) ) ) {
            return fail( "Regression line draws the line: accent ink " + with + " with, " + without + " without" );
        }
        // the date axis: years over six years, months over less, round numbers where the dates are not calendar
        if ( !panel.ticksX().get( 0 )._label.equals( "2018" ) ) {
            return fail( "calendar years over six years: " + panel.ticksX().get( 0 )._label );
        }
        panel.setData( data, false );
        if ( !panel.xTitle().equals( "Date (year)" ) ) {
            return fail( "a plain date axis names its unit: " + panel.xTitle() );
        }
        // an ancestor drawn on a tip's spot is not of the tip's dot, and the tip wins the tie: it is the datum
        if ( !ancestorOnATipsSpotOk( data ) ) {
            return false;
        }
        // ages run the other way: time still runs left to right
        return agesOk();
    }

    private static boolean ancestorOnATipsSpotOk( final ClockPlot.Data data ) {
        ClockPlot.Point tip = null;
        ClockPlot.Point ancestor = null;
        for( final ClockPlot.Point p : data._points ) {
            if ( p._tip && "Asia/6/2022".equals( p._node.getName() ) ) {
                tip = p;
            }
            if ( !p._tip && ( ancestor == null ) ) {
                ancestor = p;
            }
        }
        final List<ClockPlot.Point> points = new ArrayList<ClockPlot.Point>();
        // the ancestor FIRST, so it is the nearest found first and only the tie rule can hand the dot to the tip
        points.add( new ClockPlot.Point( ancestor._node, tip._date, tip._div, false ) );
        points.addAll( data._points );
        final ClockPlotPanel panel = new ClockPlotPanel( new RecordingHost() );
        panel.setShowInternal( true );
        panel.setData( new ClockPlot.Data( data._forward, data._unit, data._from_rates, points, data._root, data._root_date,
                                           data._root_div, data._fit ),
                       true );
        ClockPlotPanel.Mark m = null;
        for( final ClockPlotPanel.Mark k : panel.marksForTest() ) {
            if ( k._p == tip ) {
                m = k;
            }
        }
        final List<ClockPlotPanel.Mark> at = panel.marksAt( m._px, m._py );
        if ( ( at.size() != 1 ) || ( at.get( 0 )._p != tip ) ) {
            return fail( "an ancestor on a tip's spot: the dot is the tip's alone, got " + at.size() );
        }
        // the same SPOT, along each axis: tips a few pixels off to the side, or above, are dots of their own
        final double px_per_year = Math.abs( panel.x( tip._date + 1 ) - panel.x( tip._date ) );
        final double px_per_div = Math.abs( panel.y( tip._div + 1e-3 ) - panel.y( tip._div ) ) / 1e-3;
        final List<ClockPlot.Point> near = new ArrayList<ClockPlot.Point>( data._points );
        near.add( new ClockPlot.Point( tip._node, tip._date + ( 3 / px_per_year ), tip._div, true ) );
        near.add( new ClockPlot.Point( tip._node, tip._date, tip._div + ( 3 / px_per_div ), true ) );
        final ClockPlotPanel p2 = new ClockPlotPanel( new RecordingHost() );
        p2.setData( new ClockPlot.Data( data._forward, data._unit, data._from_rates, near, data._root, data._root_date,
                                        data._root_div, data._fit ),
                    true );
        for( final ClockPlotPanel.Mark k : p2.marksForTest() ) {
            if ( k._p == tip ) {
                final List<ClockPlotPanel.Mark> one = p2.marksAt( k._px, k._py );
                if ( ( one.size() != 1 ) || ( one.get( 0 )._p != tip ) ) {
                    return fail( "tips 3 px to the side or above are not of the dot, got " + one.size() );
                }
            }
        }
        return true;
    }

    private static boolean agesOk() throws Exception {
        final Phylogeny t = read( "beast-tip-dates.nex" );
        final BranchLengthLayout.TimeLengths time = BranchLengthLayout.TimeLengths.onScreen( t );
        final ClockPlot.Data d = ClockPlot.data( t, time, null );
        if ( ( d == null ) || !d._forward || !d._from_rates ) {
            return fail( "beast-tip-dates.nex: calendar dates, divergence from the rates" );
        }
        final ClockPlotPanel panel = new ClockPlotPanel( new RecordingHost() );
        panel.setData( d, true );
        if ( !( panel.x( 2010 ) < panel.x( 2012 ) ) ) {
            return fail( "calendar dates: earlier to the left" );
        }
        final ClockPlot.Data ages = new ClockPlot.Data( false, null, true, d._points, d._root, d._root_date, d._root_div, d._fit );
        panel.setData( ages, false );
        if ( !( panel.x( 2012 ) < panel.x( 2010 ) ) || !panel.xTitle().equals( "Age" ) ) {
            return fail( "ages: the older (larger) to the left, the axis called Age" );
        }
        return true;
    }

    private static void click( final ClockPlotPanel panel, final ClockPlotPanel.Mark m ) {
        panel.press( (int) Math.round( m._px ), (int) Math.round( m._py ) );
        panel.release( (int) Math.round( m._px ), (int) Math.round( m._py ) );
    }

    private static void boxDrag( final ClockPlotPanel panel, final int xa, final int ya, final int xb, final int yb ) {
        panel.press( xa, ya );
        panel.drag( ( xa + xb ) / 2, ( ya + yb ) / 2 );
        panel.drag( xb, yb );
        panel.release( xb, yb );
    }

    /** Pixels in the accent's hue (its core, every channel within 40). */
    private static int accentInk( final ClockPlotPanel panel ) {
        final BufferedImage img = new BufferedImage( ClockPlotPanel.W, ClockPlotPanel.H, BufferedImage.TYPE_INT_RGB );
        final java.awt.Graphics2D g = img.createGraphics();
        g.setColor( java.awt.Color.WHITE );
        g.fillRect( 0, 0, ClockPlotPanel.W, ClockPlotPanel.H );
        g.setRenderingHint( java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON );
        panel.paintPlot( g );
        g.dispose();
        final java.awt.Color a = ClockPlotPanel.accent();
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                if ( ( Math.abs( ( ( rgb >> 16 ) & 0xff ) - a.getRed() ) < 40 ) && ( Math.abs( ( ( rgb >> 8 ) & 0xff ) - a.getGreen() ) < 40 )
                        && ( Math.abs( ( rgb & 0xff ) - a.getBlue() ) < 40 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    // ---- in a real frame -----------------------------------------------------------------------------------------

    private static void paint( final TreePanel tp ) {
        tp.paintImmediately( 0, 0, Math.max( 1, tp.getWidth() ), Math.max( 1, tp.getHeight() ) );
    }

    private static boolean windowOk() throws Exception {
        final Phylogeny one_date = read( "clock-plot-one-date.nex" );
        final Phylogeny clock = read( "clock-plot.nex" );
        final MainFrame[] mf = new MainFrame[ 1 ];
        // several trees: the LAST tab is selected
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { one_date, clock }, new Configuration(), "clock" ) );
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> inFrame( mf[ 0 ], ok ) );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> {
                if ( mf[ 0 ].clockPlotWindow() != null ) {
                    mf[ 0 ].clockPlotWindow().dispose();
                }
                mf[ 0 ].dispose();
            } );
        }
        return ok[ 0 ];
    }

    private static void inFrame( final MainFrame frame, final boolean[] ok ) {
        final ControlPanel cp = frame.getMainPanel().getControlPanel();
        final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
        if ( !"clock_plot".equals( tp.getPhylogeny().getName() ) && ( tp.getPhylogeny().getNumberOfExternalNodes() != 20 ) ) {
            fail( ok, "precondition: the clock-plot demo is the current tab" );
            return;
        }
        cp.populateBranchLengthsControl();
        if ( !cp.isBranchLengthsControlVisible() || !cp.clockPlotButtonForTest().isVisible() ) {
            fail( ok, "clock-plot.nex: Time | Div and the Clock plot button are shown" );
        }
        if ( !ControlPanel.CLOCK_PLOT_TIP.equals( cp.clockPlotButtonForTest().getToolTipText() ) ) {
            fail( ok, "the button says what it does" );
        }
        cp.clockPlotButtonForTest().doClick();
        final ClockPlotWindow w = frame.clockPlotWindow();
        if ( ( w == null ) || !w.isVisible() || !cp.clockPlotButtonForTest().isSelected() ) {
            fail( ok, "the button opens the window and stays pressed while it is open" );
            return;
        }
        // what the line says
        final String stats = rows( w );
        if ( !stats.equals( "Rate=0.0021 subs/site per year|Root date, by the line=2017.1|Root date, in the tree=2017|R²=0.942|Tips=20|" ) ) {
            fail( ok, "what the window says: " + stats );
        }
        if ( !w.noteForTest().equals( ClockPlotWindow.NOTE_LINE ) ) {
            fail( ok, "a recorded divergence: no clock-rate note" );
        }
        final ClockPlotPanel plot = w.plotForTest();
        // pointing at the dot of four lights the four in the tree
        final List<ClockPlotPanel.Mark> four = named( plot, "Americas/" );
        ClockPlotPanel.Mark dot = null;
        for( final ClockPlotPanel.Mark m : four ) {
            if ( m._p._date == 2023.25 ) {
                dot = m;
                break;
            }
        }
        plot.pointAt( (int) Math.round( dot._px ), (int) Math.round( dot._py ) );
        if ( tp.clockPlotLitForTest().size() != 4 ) {
            fail( ok, "pointing at the dot lights its four tips in the tree, got " + tp.clockPlotLitForTest().size() );
        }
        // a click selects them in the TREE's selection, whatever Click-to says
        click( plot, dot );
        if ( ( tp.getFoundNodes0() == null ) || ( tp.getFoundNodes0().size() != 4 ) ) {
            fail( ok, "a click on the dot selects its four tips in the tree (found set 0)" );
        }
        paint( tp );
        if ( !w.deselectButtonForTest().isEnabled() ) {
            fail( ok, "Deselect all is enabled while something is selected" );
        }
        w.deselectButtonForTest().doClick();
        paint( tp );
        if ( ( tp.getFoundNodes0() != null ) && !tp.getFoundNodes0().isEmpty() ) {
            fail( ok, "Deselect all clears the tree's selection" );
        }
        if ( w.deselectButtonForTest().isEnabled() ) {
            fail( ok, "Deselect all is disabled with nothing selected" );
        }
        plot.leave();
        if ( !tp.clockPlotLitForTest().isEmpty() ) {
            fail( ok, "leaving the plot puts the light out" );
        }
        // a selection made in the tree asks the plot to repaint; a paint that changes nothing does not. Read off the
        // RepaintManager: the plot paints its colours live, so a picture of it cannot tell whether it was asked
        final javax.swing.RepaintManager rm = javax.swing.RepaintManager.currentManager( plot );
        paint( tp );
        rm.paintDirtyRegions();
        final java.util.Set<Long> found = new java.util.HashSet<Long>();
        found.add( tp.getPhylogeny().getNode( "Europe/2/2018" ).getId() );
        tp.setFoundNodes0( found );
        paint( tp );
        if ( rm.getDirtyRegion( plot ).isEmpty() ) {
            fail( ok, "a selection made in the tree repaints the plot" );
        }
        rm.paintDirtyRegions();
        paint( tp );
        if ( !rm.getDirtyRegion( plot ).isEmpty() ) {
            fail( ok, "a paint that changes nothing repaints nothing: " + rm.getDirtyRegion( plot ) );
        }
        tp.setFoundNodes0( null );
        paint( tp );
        // the colour pass: skipped by a paint that changes nothing, done again after an edit
        final int passes = w.colourPassesForTest();
        paint( tp );
        if ( w.colourPassesForTest() != passes ) {
            fail( ok, "a paint that changes nothing works out no colour" );
        }
        // with Color by OFF, so the edit cannot reach the key through a rebuilt colour scheme: only the edit itself
        tp.setColorByPropertyRef( null );
        paint( tp );
        final int off = w.colourPassesForTest();
        final java.util.List<Object> before = tp.clockPlotColourKey();
        tp.setEdited( true );
        if ( !before.subList( 0, before.size() - 1 ).equals( tp.clockPlotColourKey().subList( 0, before.size() - 1 ) ) ) {
            fail( ok, "precondition: the edit changed nothing in the colour key but the edit count" );
        }
        paint( tp );
        if ( ( off != ( passes + 1 ) ) || ( w.colourPassesForTest() != ( off + 1 ) ) ) {
            fail( ok, "Color by off recolours once; then an edit (a node style may have changed) recolours again: "
                    + passes + " " + off + " " + w.colourPassesForTest() );
        }
        // a node pointed at in the TREE rings its point
        final PhylogenyNode asia6 = tp.getPhylogeny().getNode( "Asia/6/2022" );
        tp.setHoverForTest( asia6, false );
        if ( ( plot.ringedForTest() == null ) || ( plot.ringedForTest()._p._node != asia6 ) ) {
            fail( ok, "pointing at a tip in the tree rings its point" );
        }
        tp.setHoverForTest( null, false );
        if ( plot.ringedForTest() != null ) {
            fail( ok, "pointing away takes the ring off" );
        }
        tp.setHoverForTest( asia6, false );
        tp.setNodeInPreorderToNull(); // what a collapse, a paste, a cut do to the node under the pointer
        if ( plot.ringedForTest() != null ) {
            fail( ok, "a change of the tree's structure takes the ring off with the node pointed at" );
        }
        // a render that changes nothing works nothing out again
        final int lays = plot.laysForTest();
        paint( tp );
        paint( tp );
        if ( plot.laysForTest() != lays ) {
            fail( ok, "a paint that changes nothing lays nothing down again" );
        }
        // ...and one after a change the Time | Div offer is asked again for does
        tp.invalidateBranchLengthToggle();
        paint( tp );
        if ( plot.laysForTest() != ( lays + 1 ) ) {
            fail( ok, "a change in place (the node editor) works the plot out again" );
        }
        // a tip deleted: another node count, the plot is worked out again
        final PhylogenyNode gone = tp.getPhylogeny().getNode( "Europe/1/2018" );
        tp.getPhylogeny().deleteSubtree( gone, true );
        tp.getPhylogeny().externalNodesHaveChanged();
        tp.setNodeInPreorderToNull();
        paint( tp );
        if ( plot.getData().tipCount() != 19 ) {
            fail( ok, "a tip deleted: the plot has 19 tips, got " + plot.getData().tipCount() );
        }
        // the switch to Div changes no number
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        paint( tp );
        final String stats19 = rows( w );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        paint( tp );
        if ( !rows( w ).equals( stats19 ) || !stats19.contains( "Tips=19|" ) ) {
            fail( ok, "in Div the plot says what it says in Time: " + stats19 + " / " + rows( w ) );
        }
        // Internal nodes, kept while the program runs
        w.internalCheckBoxForTest().doClick();
        if ( ( plot.marksForTest().size() != 28 ) || !ClockPlotWindow.internalOptionForTest() ) {
            fail( ok, "Internal nodes adds the 9 ancestors to the 19 tips, and is kept: " + plot.marksForTest().size() );
        }
        w.internalCheckBoxForTest().doClick();
        // a collapsed clade: pointing at a tip hidden in it lights the clade
        final PhylogenyNode americas = tp.getPhylogeny().getNode( "NODE_Americas_2022.9" );
        final ClockPlotPanel.Mark now = plot.marksAt( dot._px, dot._py ).get( 0 );
        plot.pointAt( (int) Math.round( now._px ), (int) Math.round( now._py ) );
        if ( tp.clockPlotLitForTest().size() != 4 ) {
            fail( ok, "precondition: the dot's four tips are lit before the collapse" );
        }
        // the pointer stays on the dot while the clade is collapsed: lit again, as the tree now shows them
        tp.collapse( americas );
        paint( tp );
        if ( ( tp.clockPlotLitForTest().size() != 1 ) || ( tp.clockPlotLitForTest().get( 0 ) != americas ) ) {
            fail( ok, "tips hidden in a collapsed clade light the clade, once: " + tp.clockPlotLitForTest() );
        }
        plot.leave();
        tp.collapse( americas );
        paint( tp );
        // a subtree view: the clade's own plot, its own line; the button stays (the rule is the tree's)
        tp.subTree( tp.getPhylogeny().getNode( "NODE_Asia" ) );
        paint( tp );
        cp.populateBranchLengthsControl();
        if ( w.plotForTest().getData().tipCount() != 6 ) {
            fail( ok, "in the view of the Asia clade, its six tips: " + w.plotForTest().getData().tipCount() );
        }
        if ( !rows( w ).contains( "Tips=6 (the clade on view)" ) || !cp.clockPlotButtonForTest().isVisible() ) {
            fail( ok, "the clade on view is said, and the button stays: " + rows( w ) );
        }
        tp.superTree();
        paint( tp );
        if ( w.plotForTest().getData().tipCount() != 19 ) {
            fail( ok, "back to the whole tree, its 19 tips" );
        }
        // a tab whose tree has no plot closes the window, and the button follows
        frame.getMainPanel().getTabbedPane().setSelectedIndex( 0 );
        final TreePanel other = frame.getMainPanel().getCurrentTreePanel();
        paint( other );
        if ( ( frame.clockPlotWindow() != null ) || w.isDisplayable() ) {
            fail( ok, "a tab with no plot closes the window" );
        }
        cp.populateBranchLengthsControl();
        if ( !cp.isBranchLengthsControlVisible() || cp.clockPlotButtonForTest().isVisible() ) {
            fail( ok, "clock-plot-one-date.nex: Time | Div shown, the Clock plot button not" );
        }
        if ( cp.clockPlotButtonForTest().isSelected() ) {
            fail( ok, "the button is up once the window is gone" );
        }
        // back, open again, and close: closed AT ONCE, not when a later event says so
        frame.getMainPanel().getTabbedPane().setSelectedIndex( 1 );
        cp.populateBranchLengthsControl();
        cp.clockPlotButtonForTest().doClick();
        final ClockPlotWindow w2 = frame.clockPlotWindow();
        if ( w2 == null ) {
            fail( ok, "opened again" );
            return;
        }
        cp.clockPlotButtonForTest().doClick();
        if ( ( frame.clockPlotWindow() != null ) || cp.clockPlotButtonForTest().isSelected() ) {
            fail( ok, "the button closes the window at once" );
        }
        cp.clockPlotButtonForTest().doClick();
        final ClockPlotWindow w3 = frame.clockPlotWindow();
        w3.dispose(); // the title bar's close
        if ( ( frame.clockPlotWindow() != null ) || cp.clockPlotButtonForTest().isSelected() ) {
            fail( ok, "the window's own close is at once too, and the button follows" );
        }
        // the clock-rate tree gets the note
        final Phylogeny beast;
        try {
            beast = read( "beast-tip-dates.nex" );
        }
        catch ( final Exception e ) {
            fail( ok, "beast-tip-dates.nex unreadable" );
            return;
        }
        frame.getMainPanel().addPhylogenyInNewTab( beast, frame.getConfiguration(), "beast", null );
        final TreePanel btp = frame.getMainPanel().getCurrentTreePanel();
        cp.populateBranchLengthsControl();
        cp.clockPlotButtonForTest().doClick();
        final ClockPlotWindow w4 = frame.clockPlotWindow();
        if ( ( w4 == null ) || !w4.noteForTest().equals( ClockPlotWindow.NOTE_RATES + ClockPlotWindow.NOTE_LINE ) ) {
            fail( ok, "a clock-rate tree: the window says its divergence is time x rate" );
        }
        else if ( !rows( w4 ).startsWith( "Rate=0.004 subs/site per year|Root date, by the line=2008.3|Root date, in the tree=2008.3|R²=1|Tips=10" ) ) {
            fail( ok, "beast-tip-dates.nex: " + rows( w4 ) );
        }
        if ( btp != frame.getMainPanel().getCurrentTreePanel() ) {
            fail( ok, "precondition: the BEAST tab is current" );
        }
        // closing the tab closes its plot, at once (no paint of another tab need follow)
        frame.getMainPanel().closeCurrentPane();
        if ( ( frame.clockPlotWindow() != null ) || ( ( w4 != null ) && w4.isDisplayable() ) ) {
            fail( ok, "closing the tab closes its clock plot" );
        }
    }

    private static String rows( final ClockPlotWindow w ) {
        final StringBuilder sb = new StringBuilder();
        for( final String[] r : w.statsForTest() ) {
            sb.append( r[ 0 ] ).append( '=' ).append( r[ 1 ] ).append( '|' );
        }
        return sb.toString();
    }

    private ClockPlotWindowTest() {
    }
}
