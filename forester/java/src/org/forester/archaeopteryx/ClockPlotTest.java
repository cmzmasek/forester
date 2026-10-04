// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import java.io.File;
import java.math.BigDecimal;
import java.util.List;
import java.util.Map;

import org.forester.archaeopteryx.tools.TipDateExtractor;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * {@link ClockPlot}, the pure half of the clock plot: which trees have one, its points and its line, the ticks and the
 * numbers as the window prints them. Ported from Archaeopteryx.js's {@code test/clock_test.js}, and like it, every
 * expected number was computed BY HAND from the five tips below, not read off the code:
 *
 * <pre>
 *   tip   date   divergence
 *   A     2000   0.001
 *   B     2001   0.003
 *   C     2002   0.002
 *   D     2003   0.006
 *   E     2004   0.005
 *
 *   means 2002 and 0.0034; Sxx = 10, Sxy = 0.011, Syy = 17.2e-6
 *   slope 0.0011, intercept 0.0034 - 0.0011 x 2002 = -2.1988
 *   R2 = 0.011^2 / (10 x 17.2e-6) = 121 / 172
 *   the line reaches the root's divergence (0) at 2.1988 / 0.0011 = 1998.9090...
 *   off the line: A -0.0002, B +0.0007, C -0.0014, D +0.0015, E -0.0006
 * </pre>
 *
 * and for the clade of C, D and E (its root Y: 2000.5, 0.001):
 *
 * <pre>
 *   means 2003 and 0.013 / 3; Sxx = 2, Sxy = 0.003, Syy = 26e-6 / 3
 *   slope 0.0015, R2 = 0.003^2 / (2 x 26e-6 / 3) = 27 / 52
 *   the line reaches Y's divergence (0.001) at 2000 + 7 / 9
 * </pre>
 *
 * The tree: ROOT(1999, 0) -> X(1999.5, 0.0005) -> A, B; ROOT -> Y(2000.5, 0.001) -> C, D, E. Its branch lengths are
 * the gaps between its dates, so it is in time, as an Auspice build arrives.
 * <p>
 * The SECOND kind, a divergence tree whose tips alone are dated (Archaeopteryx.js e22361c, 2d7fabe; the fixtures and
 * numbers of its {@code test/clock_test.js}), from the Nexus
 * {@code ((A[&num_date=2000]:0.01,B[&num_date=2001]:0.03)X:0.005,(C[&num_date=2002]:0.02,D[&num_date=2003]:0.06,E[&num_date=2004]:0.04)Y:0.01)}:
 *
 * <pre>
 *   tip   date   path from the root
 *   A     2000   0.005 + 0.01 = 0.015
 *   B     2001   0.005 + 0.03 = 0.035
 *   C     2002   0.01 + 0.02 = 0.03
 *   D     2003   0.01 + 0.06 = 0.07
 *   E     2004   0.01 + 0.04 = 0.05
 *
 *   means 2002 and 0.04; Sxx = 10, Sxy = 0.105, Syy = 0.00175
 *   slope 0.0105, intercept 0.04 - 0.0105 x 2002 = -20.981
 *   R2 = 0.105^2 / (10 x 0.00175) = 0.63
 *   the line reaches the root (0) at 20.981 / 0.0105 = 1998.190476...
 * </pre>
 *
 * and for the clade Y = (C, D, E), itself 0.01 from the root: means 2003 and 0.05; Sxx = 2, Sxy = 0.02, Syy = 0.0008;
 * slope 0.01, intercept -19.98, R2 0.25; the line reaches Y's divergence (0.01) at 1999. With the dates in the tips'
 * NAMES instead (A_2000, B|2001-07-02, ...: a year alone is mid-year, and 2 July 2001 is the middle of a non-leap
 * year exactly) every date shifts by 0.5: the slope and R2 stand, the intercept is -20.98625 and the root by the
 * line 1998.690476...
 */
public final class ClockPlotTest {

    private static final String DIV  = BranchLengthLayout.DIV_PROPERTY_REF;
    private static final String RATE = BranchLengthLayout.RATE_PROPERTY_REF;

    public static void main( final String[] args ) {
        System.out.println( test() ? "ClockPlotTest: OK." : "ClockPlotTest: FAILED." );
    }

    public static boolean test() {
        boolean ok = true;
        ok &= testRule();
        ok &= testWholeTree();
        ok &= testAncestorsTakeNoPart();
        ok &= testClade();
        ok &= testAges();
        ok &= testNoSignal();
        ok &= testRegression();
        ok &= testDivergenceFromRoot();
        ok &= testMonthTicks();
        ok &= testYearTicks();
        ok &= testD3Ticks();
        ok &= testNumbers();
        ok &= testDemoPair();
        ok &= testRootToTipRule();
        ok &= testRootToTipData();
        ok &= testRootToTipAges();
        ok &= testRootToTipNames();
        ok &= testRootToTipFromLabels();
        ok &= testSwineH1CrossRun();
        ok &= testRootToTipDemoPair();
        return ok;
    }

    // ---- fixtures ------------------------------------------------------------------------------------------------

    private static final double[][] TIPS      = { { 2000, 0.001 }, { 2001, 0.003 }, { 2002, 0.002 }, { 2003, 0.006 },
            { 2004, 0.005 } };
    private static final String[]   TIP_NAMES = { "A", "B", "C", "D", "E" };

    private static void date( final PhylogenyNode n, final double v ) {
        final Date d = new Date();
        d.setValue( new BigDecimal( String.valueOf( v ) ) );
        d.setUnit( "year" );
        n.getNodeData().setDate( d );
    }

    private static void prop( final PhylogenyNode n, final String ref, final double v ) {
        if ( n.getNodeData().getProperties() == null ) {
            n.getNodeData().setProperties( new PropertiesList() );
        }
        n.getNodeData().getProperties().addProperty( new Property( ref, String.valueOf( v ), "", "xsd:decimal",
                                                                   AppliesTo.NODE ) );
    }

    private static PhylogenyNode node( final String name, final Double date, final double div, final PhylogenyNode parent ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        if ( date != null ) {
            date( n, date );
        }
        prop( n, DIV, div );
        if ( parent != null ) {
            parent.addAsChild( n );
            if ( ( date != null ) && parent.getNodeData().isHasDate() ) {
                n.setDistanceToParent( date - parent.getNodeData().getDate().getValue().doubleValue() );
            }
            else {
                n.setDistanceToParent( 1.0 );
            }
        }
        return n;
    }

    /** The five-tip build. {@code tips}: [date, div] per tip, a null row takes the tip out; {@code x} and {@code y}:
     *  the two ancestors' [date, div], a null date leaves one undated. */
    private static Phylogeny build( final double[][] tips, final Double[] x, final Double[] y ) {
        final PhylogenyNode root = node( "ROOT", 1999.0, 0, null );
        final Double[] xx = ( x != null ) ? x : new Double[] { 1999.5, 0.0005 };
        final Double[] yy = ( y != null ) ? y : new Double[] { 2000.5, 0.001 };
        final PhylogenyNode nx = node( "X", xx[ 0 ], xx[ 1 ], root );
        final PhylogenyNode ny = node( "Y", yy[ 0 ], yy[ 1 ], root );
        final double[][] t = ( tips != null ) ? tips : TIPS;
        for( int i = 0; i < 5; ++i ) {
            if ( t[ i ] != null ) {
                node( TIP_NAMES[ i ], t[ i ][ 0 ], t[ i ][ 1 ], ( i < 2 ) ? nx : ny );
            }
        }
        final Phylogeny p = new Phylogeny();
        p.setRoot( root );
        p.setRooted( true );
        p.externalNodesHaveChanged();
        return p;
    }

    private static Phylogeny build() {
        return build( null, null, null );
    }

    private static BranchLengthLayout.TimeLengths time( final Phylogeny p ) {
        return BranchLengthLayout.TimeLengths.onScreen( p );
    }

    private static ClockPlot.Data data( final Phylogeny p ) {
        return ClockPlot.data( p, time( p ), null );
    }

    private static boolean close( final double a, final double b ) {
        return close( a, b, 1e-9 );
    }

    private static boolean close( final double a, final double b, final double rel ) {
        return Math.abs( a - b ) <= ( Math.max( Math.abs( a ), Math.abs( b ) ) * rel );
    }

    private static PhylogenyNode named( final Phylogeny p, final String name ) {
        return p.getNode( name );
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [ClockPlotTest] " + m );
        return false;
    }

    // ---- the rule ------------------------------------------------------------------------------------------------

    private static boolean testRule() {
        final Phylogeny whole = build();
        if ( !BranchLengthLayout.isApplicable( whole ) || !ClockPlot.isOffered( whole, time( whole ) ) ) {
            return fail( "the build of five tips on five dates has a plot" );
        }
        // every tip on one date: Time | Div is offered, the plot is not
        final double[][] one = { { 2003, 0.001 }, { 2003, 0.003 }, { 2003, 0.002 }, { 2003, 0.006 }, { 2003, 0.005 } };
        final Phylogeny same = build( one, null, null );
        if ( !BranchLengthLayout.isApplicable( same ) ) {
            return fail( "control: the one-date build is still offered Time | Div" );
        }
        if ( ClockPlot.isOffered( same, time( same ) ) || ( data( same ) != null ) ) {
            return fail( "every tip on one date: no plot" );
        }
        // ...and one tip off that date is enough
        final double[][] off = { { 2002, 0.001 }, { 2003, 0.003 }, { 2003, 0.002 }, { 2003, 0.006 }, { 2003, 0.005 } };
        final Phylogeny one_off = build( off, null, null );
        if ( !ClockPlot.isOffered( one_off, time( one_off ) ) ) {
            return fail( "one tip on another date: a plot" );
        }
        // two tips, on two dates: Time | Div is offered, the plot is not
        final Phylogeny two = build( new double[][] { TIPS[ 0 ], null, TIPS[ 2 ], null, null }, null, null );
        if ( ( two.getNumberOfExternalNodes() != 2 ) || !BranchLengthLayout.isApplicable( two ) ) {
            return fail( "control: two tips, and still offered Time | Div" );
        }
        if ( ClockPlot.isOffered( two, time( two ) ) ) {
            return fail( "two tips: no plot" );
        }
        final Phylogeny three = build( new double[][] { TIPS[ 0 ], TIPS[ 1 ], TIPS[ 2 ], null, null }, null, null );
        if ( ( three.getNumberOfExternalNodes() != 3 ) || !ClockPlot.isOffered( three, time( three ) ) ) {
            return fail( "three tips: a plot" );
        }
        // no Time | Div, no plot, whatever the tips state: one ancestor undated
        final Phylogeny gap = build( null, new Double[] { null, 0.0005 }, null );
        if ( BranchLengthLayout.isApplicable( gap ) ) {
            return fail( "control: one ancestor undated, and no Time | Div" );
        }
        if ( ClockPlot.isOffered( gap, time( gap ) ) || ( data( gap ) != null ) ) {
            return fail( "no Time | Div: no plot" );
        }
        if ( ClockPlot.isOffered( null, null ) || ClockPlot.isOffered( new Phylogeny(), new BranchLengthLayout.TimeLengths() ) ) {
            return fail( "no tree, an empty tree: no plot" );
        }
        return true;
    }

    // ---- the plot ------------------------------------------------------------------------------------------------

    private static boolean testWholeTree() {
        final Phylogeny phy = build();
        final ClockPlot.Data d = data( phy );
        if ( ( d == null ) || !d._forward || !"year".equals( d._unit ) || d._from_rates ) {
            return fail( "what the plot is: forward, unit year, recorded divergence" );
        }
        if ( d._points.size() != 8 ) {
            return fail( "a point per node: five tips and three ancestors, got " + d._points.size() );
        }
        final StringBuilder got = new StringBuilder();
        final Map<String, Double> want = Map.of( "A", -0.0002, "B", 0.0007, "C", -0.0014, "D", 0.0015, "E", -0.0006 );
        ClockPlot.Point x = null;
        for( final ClockPlot.Point p : d._points ) {
            if ( p._tip ) {
                got.append( p._node.getName() ).append( ' ' ).append( ClockPlot.plain( p._date ) ).append( ' ' )
                        .append( p._div ).append( ", " );
                if ( Math.abs( p._residual - want.get( p._node.getName() ) ) > 1e-12 ) {
                    return fail( "off the line, " + p._node.getName() + ": " + p._residual );
                }
            }
            else if ( "X".equals( p._node.getName() ) ) {
                x = p;
            }
        }
        if ( !got.toString().equals( "A 2000 0.001, B 2001 0.003, C 2002 0.002, D 2003 0.006, E 2004 0.005, " ) ) {
            return fail( "the tips, each at its date and divergence: " + got );
        }
        if ( ( x == null ) || ( x._date != 1999.5 ) || ( x._div != 0.0005 ) || ( x._residual != null ) ) {
            return fail( "an ancestor is a point too, at its own date and divergence, with no residual" );
        }
        final ClockPlot.Fit f = d._fit;
        if ( ( f == null ) || ( f._n != 5 ) ) {
            return fail( "the line counts the five tips and nothing else" );
        }
        if ( !close( f._slope, 0.0011 ) || !close( f._intercept, -2.1988 ) || !close( f._r2, 121.0 / 172 ) ) {
            return fail( "slope 0.0011, intercept -2.1988, R2 121/172: " + f._slope + " " + f._intercept + " " + f._r2 );
        }
        if ( !close( f._rate, 0.0011 ) || !close( f._root_date, 2.1988 / 0.0011 ) ) {
            return fail( "rate 0.0011; the root by the line at 1998.909: " + f._rate + " " + f._root_date );
        }
        if ( ( d._root != phy.getRoot() ) || ( d._root_date != 1999 ) || ( d._root_div != 0 ) ) {
            return fail( "the root as the tree states it" );
        }
        if ( d.tipCount() != 5 ) {
            return fail( "five tips" );
        }
        return true;
    }

    /** The ancestors' dates were inferred with a clock: moved anywhere, they leave the line where it was. */
    private static boolean testAncestorsTakeNoPart() {
        final ClockPlot.Fit a = data( build() )._fit;
        final ClockPlot.Fit b = data( build( null, new Double[] { 1999.9, 0.0009 }, new Double[] { 1999.2, 0.0001 } ) )._fit;
        if ( ( a._slope != b._slope ) || ( a._intercept != b._intercept ) || !a._r2.equals( b._r2 ) || ( a._n != b._n ) ) {
            return fail( "the fit moved with the ancestors" );
        }
        final double[][] moved = { TIPS[ 0 ], TIPS[ 1 ], TIPS[ 2 ], TIPS[ 3 ], { 2004, 0.009 } };
        if ( data( build( moved, null, null ) )._fit._slope == a._slope ) {
            return fail( "control: a tip moved, and the line did not" );
        }
        return true;
    }

    /** A subtree view, made the way the tab makes one (a copy of the clade's top over the shared nodes). */
    private static boolean testClade() {
        final Phylogeny phy = build();
        final PhylogenyNode y = named( phy, "Y" );
        final Phylogeny view = TreePanelUtil.subTree( y, phy );
        final ClockPlot.Data d = ClockPlot.data( phy, time( phy ), view );
        if ( ( d._points.size() != 4 ) || ( d._root != view.getRoot() ) || ( d._root_date != 2000.5 ) || ( d._root_div != 0.001 ) ) {
            return fail( "the clade of Y: four points under its own root" );
        }
        final ClockPlot.Fit f = d._fit;
        if ( ( f == null ) || ( f._n != 3 ) || !close( f._slope, 0.0015 ) || !close( f._intercept, ( 0.013 / 3 ) - ( 0.0015 * 2003 ) )
                || !close( f._r2, 27.0 / 52 ) ) {
            return fail( "the clade's own line: slope 0.0015, R2 27/52" );
        }
        // where the line reaches the CLADE's root divergence, not zero
        if ( !close( f._root_date, 2000 + ( 7.0 / 9 ) ) ) {
            return fail( "the clade's root by the line at 2000.777: " + f._root_date );
        }
        // the rule was the tree's: a clade of two tips draws its points and no line
        final Phylogeny phy2 = build();
        final ClockPlot.Data x = ClockPlot.data( phy2, time( phy2 ), TreePanelUtil.subTree( named( phy2, "X" ), phy2 ) );
        if ( ( x == null ) || ( x._points.size() != 3 ) || ( x._fit != null ) ) {
            return fail( "a clade of two tips: three points and no line" );
        }
        return true;
    }

    /**
     * A BEAST clock-model tree stating AGES (0 at the youngest tip, 4 at the root) and one rate, 0.01, on every branch.
     * A tip's divergence is 0.01 x (4 - its age): A 0.04, B 0.03, C 0.035, D 0.02 -- on the line div = 0.04 - 0.01 x
     * age exactly. ((A:2, B:1):2, (C:2.5, D:1):1), node ages 2 and 3.
     */
    private static Phylogeny agesTree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.setName( "R" );
        dateNoUnit( root, 4 );
        final PhylogenyNode p = rated( "P", 2, 2, root );
        final PhylogenyNode q = rated( "Q", 3, 1, root );
        rated( "A", 0, 2, p );
        rated( "B", 1, 1, p );
        rated( "C", 0.5, 2.5, q );
        rated( "D", 2, 1, q );
        final Phylogeny t = new Phylogeny();
        t.setRoot( root );
        t.setRooted( true );
        t.externalNodesHaveChanged();
        return t;
    }

    private static boolean testAges() {
        final Phylogeny t = agesTree();
        if ( !BranchLengthLayout.isApplicable( t ) ) {
            return fail( "control: the ages tree is offered Time | Div" );
        }
        final ClockPlot.Data d = data( t );
        if ( ( d == null ) || d._forward || !d._from_rates || ( d._unit != null ) ) {
            return fail( "ages, their divergence from the rates, no unit" );
        }
        final StringBuilder tips = new StringBuilder();
        for( final ClockPlot.Point pt : d._points ) {
            if ( pt._tip ) {
                tips.append( pt._node.getName() ).append( ' ' ).append( ClockPlot.plain( pt._date ) ).append( ' ' )
                        .append( ClockPlot.number( pt._div, 12 ) ).append( ", " );
            }
        }
        if ( !tips.toString().equals( "A 0 0.04, B 1 0.03, C 0.5 0.035, D 2 0.02, " ) ) {
            return fail( "each tip at its age and at 0.01 x the time above it: " + tips );
        }
        final ClockPlot.Fit f = d._fit;
        // divergence FALLS as the age rises: the slope is negative and the rate is not
        if ( !close( f._slope, -0.01 ) || !close( f._rate, 0.01 ) || !close( f._r2, 1 ) ) {
            return fail( "slope -0.01, rate +0.01, R2 1: " + f._slope + " " + f._rate + " " + f._r2 );
        }
        if ( !close( f._root_date, 4 ) || ( d._root_date != 4 ) ) {
            return fail( "the root by the line at age 4, where the tree has it" );
        }
        return true;
    }

    private static void dateNoUnit( final PhylogenyNode n, final double v ) {
        final Date d = new Date();
        d.setValue( new BigDecimal( String.valueOf( v ) ) );
        n.getNodeData().setDate( d );
    }

    private static PhylogenyNode rated( final String name, final double age, final double length, final PhylogenyNode parent ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        dateNoUnit( n, age );
        n.setDistanceToParent( length );
        prop( n, RATE, 0.01 );
        parent.addAsChild( n );
        return n;
    }

    /** Divergence falling with time: a negative rate, reported as it is, and no date such a line "started" at. */
    private static boolean testNoSignal() {
        final double[][] falling = { { 2000, 0.006 }, { 2001, 0.005 }, { 2002, 0.004 }, { 2003, 0.004 }, { 2004, 0.003 } };
        final ClockPlot.Data d = data( build( falling, new Double[] { 1999.5, 0.0005 }, new Double[] { 2000.5, 0.001 } ) );
        if ( ( d == null ) || ( d._fit == null ) || !( d._fit._rate < 0 ) || ( d._fit._root_date != null ) ) {
            return fail( "a falling line: a negative rate and no root date" );
        }
        return true;
    }

    private static boolean testRegression() {
        if ( ClockPlot.regression( new double[] { 1, 2 }, new double[] { 1, 2 } ) != null ) {
            return fail( "two points make no line" );
        }
        if ( ClockPlot.regression( new double[] { 5, 5, 5 }, new double[] { 1, 2, 3 } ) != null ) {
            return fail( "one date makes no line" );
        }
        if ( ClockPlot.regression( new double[] { 1, 2, 3 }, new double[] { 1, 2 } ) != null ) {
            return fail( "as many dates as divergences" );
        }
        final ClockPlot.Fit flat = ClockPlot.regression( new double[] { 1, 2, 3 }, new double[] { 7, 7, 7 } );
        if ( ( flat == null ) || ( flat._slope != 0 ) || ( flat._intercept != 7 ) || ( flat._r2 != null ) ) {
            return fail( "one divergence: a flat line and no R2" );
        }
        // 1, 3, 2 against 1, 2, 3: slope 1/2, intercept 1, R2 1/4
        final ClockPlot.Fit r = ClockPlot.regression( new double[] { 1, 2, 3 }, new double[] { 1, 3, 2 } );
        if ( !close( r._slope, 0.5 ) || !close( r._intercept, 1 ) || !close( r._r2, 0.25 ) || ( r._n != 3 ) || ( r._mean_x != 2 )
                || ( r._mean_y != 2 ) ) {
            return fail( "slope 1/2, intercept 1, R2 1/4" );
        }
        // dates a day apart in 2020: about the means, dx -1.5e-3..1.5e-3, Sxx 5e-6, Sxy 4.9e-9 -> slope 9.8e-4
        final ClockPlot.Fit near = ClockPlot.regression( new double[] { 2020.001, 2020.002, 2020.003, 2020.004 },
                                                         new double[] { 0.0000010, 0.0000021, 0.0000029, 0.0000040 } );
        if ( !close( near._slope, 9.8e-4, 1e-7 ) ) {
            return fail( "slope 9.8e-4 from dates a day apart: " + near._slope );
        }
        return true;
    }

    /**
     * A clock-rate tree's divergence from the root is the SUM of each branch's as the divergence layout draws it,
     * negative spans at 0. (P at age 2 and its child B dated AFTER it, at age 2.5 -- a negative span of -0.5 -- with
     * rate 0.01: B's divergence is P's, 0.02, not 0.015.)
     */
    private static boolean testDivergenceFromRoot() {
        final PhylogenyNode root = new PhylogenyNode();
        dateNoUnit( root, 4 );
        final PhylogenyNode p = rated( "P", 2, 2, root );
        rated( "A", 0, 2, p );
        rated( "B", 2.5, -0.5, p );
        rated( "C", 1, 3, root );
        final Phylogeny t = new Phylogeny();
        t.setRoot( root );
        t.setRooted( true );
        t.externalNodesHaveChanged();
        final Map<Long, Double> div = BranchLengthLayout.divergenceFromRoot( t, time( t ) );
        if ( div == null ) {
            return fail( "control: the tree states both layouts" );
        }
        if ( div.get( root.getId() ) != 0 ) {
            return fail( "the root at 0" );
        }
        if ( !close( div.get( t.getNode( "A" ).getId() ), 0.04 ) || !close( div.get( t.getNode( "B" ).getId() ), 0.02 )
                || !close( div.get( t.getNode( "C" ).getId() ), 0.03 ) ) {
            return fail( "A 0.04, B 0.02 (its negative span adds nothing), C 0.03: " + div );
        }
        // a recording tree: the recorded value, a negative one included
        final Phylogeny rec = build( new double[][] { TIPS[ 0 ], TIPS[ 1 ], TIPS[ 2 ], TIPS[ 3 ], { 2004, -0.002 } }, null, null );
        final Map<Long, Double> rd = BranchLengthLayout.divergenceFromRoot( rec, time( rec ) );
        if ( ( rd == null ) || ( rd.get( rec.getNode( "E" ).getId() ) != -0.002 ) ) {
            return fail( "a recorded divergence is plotted as recorded, a negative one too" );
        }
        if ( BranchLengthLayout.divergenceFromRoot( build( null, new Double[] { null, 0.0005 }, null ),
                                                    new BranchLengthLayout.TimeLengths() ) != null ) {
            return fail( "a tree the layouts cannot state: none" );
        }
        return true;
    }

    // ---- ticks ---------------------------------------------------------------------------------------------------

    private static String label( final List<ClockPlot.MonthTick> ticks ) {
        final StringBuilder sb = new StringBuilder();
        for( final ClockPlot.MonthTick t : ticks ) {
            sb.append( sb.length() > 0 ? " " : "" ).append( t._year ).append( '-' ).append( t._month );
        }
        return sb.toString();
    }

    private static boolean testMonthTicks() {
        // 29 months: every sixth month is five ticks, every third would be ten
        final List<ClockPlot.MonthTick> lng = ClockPlot.calendarTickMonths( 2019.95, 2022.4 );
        if ( !label( lng ).equals( "2020-1 2020-7 2021-1 2021-7 2022-1" ) ) {
            return fail( "every sixth month over 29 months: " + label( lng ) );
        }
        // 1 July is day 183 of a leap year and day 182 of any other
        if ( ( lng.get( 0 )._value != 2020 ) || !close( lng.get( 1 )._value, 2020 + ( 182.0 / 366 ) )
                || !close( lng.get( 3 )._value, 2021 + ( 181.0 / 365 ) ) ) {
            return fail( "a tick is where its month begins" );
        }
        final List<ClockPlot.MonthTick> sht = ClockPlot.calendarTickMonths( 2020, 2020.2 );
        if ( !label( sht ).equals( "2020-1 2020-2 2020-3" ) || !close( sht.get( 1 )._value, 2020 + ( 31.0 / 366 ) )
                || !close( sht.get( 2 )._value, 2020 + ( 60.0 / 366 ) ) ) {
            return fail( "every month over ten weeks: " + label( sht ) );
        }
        if ( !label( ClockPlot.calendarTickMonths( 2020, 2020.7 ) ).equals( "2020-1 2020-3 2020-5 2020-7 2020-9" ) ) {
            return fail( "every second month over eight months" );
        }
        if ( !label( ClockPlot.calendarTickMonths( 2020.4, 2021.6 ) ).equals( "2020-7 2020-10 2021-1 2021-4 2021-7" ) ) {
            return fail( "every third month over 14 months" );
        }
        if ( !ClockPlot.calendarTickMonths( 2020, 2020 ).isEmpty() || !ClockPlot.calendarTickMonths( 2021, 2020 ).isEmpty()
                || !ClockPlot.calendarTickMonths( 1900, 2000 ).isEmpty()
                || !ClockPlot.calendarTickMonths( Double.NaN, 2000 ).isEmpty() ) {
            return fail( "no span, a backwards one, a century, no number: no ticks" );
        }
        return true;
    }

    private static boolean testYearTicks() {
        // 2017.75 - 2024.35: a span of 6.6, /7 -> 0.94 -> a step of 1
        if ( !ClockPlot.calendarTickYears( 2017.75, 2024.35 ).toString().equals( "[2018, 2019, 2020, 2021, 2022, 2023, 2024]" ) ) {
            return fail( "a tick a year over six years: " + ClockPlot.calendarTickYears( 2017.75, 2024.35 ) );
        }
        // 1940 - 2020: 80 / 7 = 11.4 -> a step of 20
        if ( !ClockPlot.calendarTickYears( 1940, 2020 ).toString().equals( "[1940, 1960, 1980, 2000, 2020]" ) ) {
            return fail( "every twenty years over eighty: " + ClockPlot.calendarTickYears( 1940, 2020 ) );
        }
        if ( ( ClockPlot.niceAxisStep( 0.94 ) != 1 ) || ( ClockPlot.niceAxisStep( 11.4 ) != 20 ) || ( ClockPlot.niceAxisStep( 3 ) != 5 )
                || ( ClockPlot.niceAxisStep( 0 ) != 1 ) ) {
            return fail( "nice steps 1, 20, 5, and 1 for nothing" );
        }
        return true;
    }

    /** d3.ticks, worked out by hand from its tickSpec (step = span / count; power, error, factor 1/2/5/10). */
    private static boolean testD3Ticks() {
        // span 1 / 5 = 0.2: power -1, error 2 >= sqrt 2 -> factor 2, inc 10 / 2 = 5 -> 0, 1/5, ... 5/5
        if ( !java.util.Arrays.toString( ClockPlot.ticks( 0, 1, 5 ) ).equals( "[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]" ) ) {
            return fail( "d3.ticks(0, 1, 5): " + java.util.Arrays.toString( ClockPlot.ticks( 0, 1, 5 ) ) );
        }
        // 0.00049 - 0.0033 over 5: step 5.62e-4, power -4, error 5.62 >= sqrt 10 -> factor 5, inc 1e4 / 5 = 2000;
        // i1 round(0.98) = 1, i2 round(6.6) = 7 -> 7 / 2000 > 0.0033 -> 6
        if ( !java.util.Arrays.toString( ClockPlot.ticks( 0.00049, 0.0033, 5 ) )
                .equals( "[5.0E-4, 0.001, 0.0015, 0.002, 0.0025, 0.003]" ) ) {
            return fail( "d3.ticks(0.00049, 0.0033, 5): " + java.util.Arrays.toString( ClockPlot.ticks( 0.00049, 0.0033, 5 ) ) );
        }
        // 1999.5 - 2004.5 over 6: step 0.833, power -1, error 8.33 >= sqrt 50 -> factor 10, inc 1 -> 2000 .. 2004
        if ( !java.util.Arrays.toString( ClockPlot.ticks( 1999.5, 2004.5, 6 ) ).equals( "[2000.0, 2001.0, 2002.0, 2003.0, 2004.0]" ) ) {
            return fail( "d3.ticks(1999.5, 2004.5, 6): " + java.util.Arrays.toString( ClockPlot.ticks( 1999.5, 2004.5, 6 ) ) );
        }
        if ( ( ClockPlot.ticks( 3, 3, 5 ).length != 1 ) || ( ClockPlot.ticks( 0, 1, 0 ).length != 0 ) ) {
            return fail( "one value for no span, none for no count" );
        }
        return true;
    }

    /** The numbers as Archaeopteryx.js prints them ({@code clockNumber}: toPrecision, toExponential). */
    private static boolean testNumbers() {
        final String[][] cases = { { ClockPlot.number( 0.0011451878920276468, 3 ), "0.00115" },
                { ClockPlot.number( 0.004000001025922829, 3 ), "0.004" }, { ClockPlot.number( 0.9999999999963703, 3 ), "1" },
                { ClockPlot.number( 0.9417217767897014, 3 ), "0.942" }, { ClockPlot.number( 5e-5, 3 ), "5.00e-5" },
                { ClockPlot.number( -1.2345e-5, 3 ), "-1.23e-5" }, { ClockPlot.number( 1.2345e8, 3 ), "1.23e+8" },
                { ClockPlot.number( 0, 3 ), "0" }, { ClockPlot.number( 0.0154, 4 ), "0.0154" },
                { ClockPlot.number( Double.NaN, 3 ), "" }, { ClockPlot.date( 2019.3796527821794, true ), "2019.38" },
                { ClockPlot.date( 2019.9, true ), "2019.9" }, { ClockPlot.date( 2017.0, true ), "2017" },
                { ClockPlot.date( 2.1988 / 0.0011, false ), "1998.9" }, { ClockPlot.signed( 0.0036516 ), "+0.00365" },
                { ClockPlot.signed( -0.0014 ), "-0.0014" }, { ClockPlot.plain( 2023.25 ), "2023.25" } };
        for( final String[] c : cases ) {
            if ( !c[ 0 ].equals( c[ 1 ] ) ) {
                return fail( "printed " + c[ 0 ] + ", want " + c[ 1 ] );
            }
        }
        return true;
    }

    // ---- the demo pair -------------------------------------------------------------------------------------------

    /**
     * forester/demo/clock-plot.nex and its refused twin, read as File > Open reads them. The numbers its README row
     * quotes are pinned here, and were also run on Archaeopteryx.js.
     */
    private static boolean testDemoPair() {
        final Phylogeny phy = read( "clock-plot.nex" );
        final Phylogeny one = read( "clock-plot-one-date.nex" );
        if ( ( phy == null ) || ( one == null ) ) {
            return fail( "the clock-plot demo pair could not be read" );
        }
        if ( BranchLengthLayout.arrivesShowingDivergence( phy ) || !ClockPlot.isOffered( phy, time( phy ) ) ) {
            return fail( "clock-plot.nex arrives in time and has a plot" );
        }
        final ClockPlot.Data d = data( phy );
        if ( ( d.tipCount() != 20 ) || ( d._fit._n != 20 ) ) {
            return fail( "20 tips in the line" );
        }
        if ( !ClockPlot.number( d._fit._rate, 3 ).equals( "0.0021" ) || !ClockPlot.number( d._fit._r2, 3 ).equals( "0.942" )
                || !ClockPlot.date( d._fit._root_date, true ).equals( "2017.1" ) || !ClockPlot.date( d._root_date, true ).equals( "2017" ) ) {
            return fail( "rate 0.0021, R2 0.942, root by the line 2017.1 over 2017 in the tree: " + d._fit._rate + " " + d._fit._r2 + " "
                    + d._fit._root_date );
        }
        // the outlier is the tip furthest off the line, 0.00365 above it
        ClockPlot.Point worst = null;
        int identical = 0;
        for( final ClockPlot.Point p : d._points ) {
            if ( p._tip && ( ( worst == null ) || ( Math.abs( p._residual ) > Math.abs( worst._residual ) ) ) ) {
                worst = p;
            }
            if ( p._tip && ( p._date == 2023.25 ) && ( p._div == 0.0128 ) ) {
                ++identical;
            }
        }
        if ( !"Asia/6/2022".equals( worst._node.getName() ) || !ClockPlot.signed( worst._residual ).equals( "+0.00365" ) ) {
            return fail( "the outlier: Asia/6/2022, +0.00365 off the line" );
        }
        if ( identical != 4 ) {
            return fail( "four identical samples on one spot, got " + identical );
        }
        // the twin: Time | Div offered, no plot
        if ( !BranchLengthLayout.isApplicable( one ) || ClockPlot.isOffered( one, time( one ) ) ) {
            return fail( "clock-plot-one-date.nex: offered Time | Div, refused the plot" );
        }
        return true;
    }

    private static Phylogeny read( final String demo ) {
        return readFile( new File( System.getProperty( "user.dir" ), "forester/demo/" + demo ) );
    }

    private static Phylogeny readFile( final File f ) {
        try {
            final Phylogeny[] phys = FigureRenderer.readTrees( f );
            return ( phys.length == 1 ) ? phys[ 0 ] : null;
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    /** A tree from its text, read the way File > Open reads a file of that {@code suffix} (".nex", ".nwk"). */
    private static Phylogeny parse( final String text, final String suffix ) {
        try {
            final File f = File.createTempFile( "clock_plot_test", suffix );
            f.deleteOnExit();
            java.nio.file.Files.writeString( f.toPath(), text );
            return readFile( f );
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    // ---- the second kind: a divergence tree whose tips alone are dated ------------------------------------------

    private static final String ROOT_TO_TIP_NEXUS = "#NEXUS\nBegin trees;\ntree T = [&R] ((A[&num_date=2000]:0.01,B[&num_date=2001]:0.03)X:0.005,"
            + "(C[&num_date=2002]:0.02,D[&num_date=2003]:0.06,E[&num_date=2004]:0.04)Y:0.01);\nEnd;\n";

    private static Phylogeny rootToTip() {
        return parse( ROOT_TO_TIP_NEXUS, ".nex" );
    }

    private static ClockPlot.Kind kind( final Phylogeny p ) {
        return ClockPlot.kind( p, time( p ) );
    }

    /** Which trees are kind 2: dated tips, undated ancestors, a length on every branch, no clock rates -- and the
     *  line's own rule on top. */
    private static boolean testRootToTipRule() {
        final Phylogeny t = rootToTip();
        if ( t == null ) {
            return fail( "the root-to-tip fixture could not be read" );
        }
        if ( BranchLengthLayout.isApplicable( t ) ) {
            return fail( "control: the ancestors are undated, so no Time | Div" );
        }
        if ( !ClockPlot.isRootToTipTree( t ) || ( kind( t ) != ClockPlot.Kind.ROOT_TO_TIP ) || !ClockPlot.isOffered( t, time( t ) ) ) {
            return fail( "dated tips over divergence lengths: a root-to-tip plot, got " + kind( t ) );
        }
        // one ancestor dated: neither kind (not a time tree, and not tips-only)
        final Phylogeny dated = rootToTip();
        date( named( dated, "X" ), 1999 );
        if ( ( kind( dated ) != null ) || ( data( dated ) != null ) ) {
            return fail( "one ancestor dated: no plot" );
        }
        // one tip undated: refused, not dropped from the plot
        final Phylogeny undated = rootToTip();
        named( undated, "C" ).getNodeData().setDate( null );
        if ( ( kind( undated ) != null ) || ( data( undated ) != null ) ) {
            return fail( "one tip undated: no plot" );
        }
        // one branch without a length
        final Phylogeny gap = parse( ROOT_TO_TIP_NEXUS.replace( "C[&num_date=2002]:0.02", "C[&num_date=2002]" ), ".nex" );
        if ( named( gap, "C" ).getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
            return fail( "fixture: C has no length, got " + named( gap, "C" ).getDistanceToParent() );
        }
        if ( kind( gap ) != null ) {
            return fail( "a branch without a length: no plot" );
        }
        // the line's own rule: three tips, not all on one date
        final Phylogeny two = parse( "#NEXUS\nBegin trees;\ntree T = [&R] (A[&num_date=2000]:0.01,B[&num_date=2001]:0.03);\nEnd;\n", ".nex" );
        if ( !ClockPlot.isRootToTipTree( two ) || ( kind( two ) != null ) ) {
            return fail( "two tips: the kind, and no plot" );
        }
        final Phylogeny same = parse( ROOT_TO_TIP_NEXUS.replaceAll( "num_date=200\\d", "num_date=2000" ), ".nex" );
        if ( !ClockPlot.isRootToTipTree( same ) || ( kind( same ) != null ) ) {
            return fail( "every tip on one date: the kind, and no plot" );
        }
        // a clock-model tree is kind 1 and never this one: its lengths are time
        final Phylogeny ages = agesTree();
        if ( ( kind( ages ) != ClockPlot.Kind.DIVERGENCE ) || ClockPlot.isRootToTipTree( ages ) ) {
            return fail( "a clock-model tree keeps its Time | Div plot" );
        }
        // ...and one with a rate on every branch but no heights -- dated tips, undated ancestors -- is neither kind
        final Phylogeny rated = parse( "#NEXUS\nBegin trees;\ntree T = [&R] ((A[&rate=0.01,num_date=2000]:2,B[&rate=0.01,num_date=2001]:1)[&rate=0.01]:2,"
                + "(C[&rate=0.01,num_date=2002]:2.5,D[&rate=0.01,num_date=2003]:1)[&rate=0.01]:1);\nEnd;\n", ".nex" );
        if ( ( BranchLengthLayout.divergenceSource( rated ) != BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE )
                || BranchLengthLayout.isApplicable( rated ) ) {
            return fail( "fixture: a rate on every branch, and no Time | Div" );
        }
        if ( ClockPlot.isRootToTipTree( rated ) || ( kind( rated ) != null ) ) {
            return fail( "a clock rate on every branch: the lengths are time, so no root-to-tip plot" );
        }
        // a plain tree with no date anywhere, in the file or in a name: no plot
        final Phylogeny plain = parse( "((A:1,B:2):1,(C:1,D:3):2);", ".nwk" );
        if ( ( plain == null ) || ClockPlot.isRootToTipTree( plain ) || ( kind( plain ) != null ) ) {
            return fail( "a plain Newick tree without dates: no plot" );
        }
        return true;
    }

    private static boolean testRootToTipData() {
        final Phylogeny t = rootToTip();
        final ClockPlot.Data d = data( t );
        if ( ( d == null ) || ( d._kind != ClockPlot.Kind.ROOT_TO_TIP ) || !d._forward || !"year".equals( d._unit ) || d._from_rates
                || ( d._div_unit != null ) || !d._calendar || ( d._from_names != 0 ) || ( d._ambiguous_names != 0 ) ) {
            return fail( "what the plot is: root-to-tip, calendar years by the unit, no divergence unit stated, no names read" );
        }
        if ( d._points.size() != 5 ) {
            return fail( "the tips only: no ancestor is dated, got " + d._points.size() );
        }
        final String[][] want = { { "A", "2000", "0.015" }, { "B", "2001", "0.035" }, { "C", "2002", "0.03" }, { "D", "2003", "0.07" },
                { "E", "2004", "0.05" } };
        for( final ClockPlot.Point p : d._points ) {
            if ( !p._tip ) {
                return fail( "an ancestor has a point" );
            }
            boolean found = false;
            for( final String[] w : want ) {
                if ( w[ 0 ].equals( p._node.getName() ) ) {
                    found = true;
                    if ( ( p._date != Double.parseDouble( w[ 1 ] ) ) || !close( p._div, Double.parseDouble( w[ 2 ] ) ) ) {
                        return fail( "tip " + w[ 0 ] + " at " + w[ 1 ] + ", " + w[ 2 ] + " from the root: " + p._date + " " + p._div );
                    }
                }
            }
            if ( !found ) {
                return fail( "an unexpected point: " + p._node.getName() );
            }
        }
        if ( ( d._root != t.getRoot() ) || d.rootDated() || ( d._root_div != 0 ) ) {
            return fail( "the root: no date of its own, divergence 0" );
        }
        final ClockPlot.Fit f = d._fit;
        if ( ( f == null ) || ( f._n != 5 ) || !close( f._slope, 0.0105 ) || !close( f._intercept, -20.981 ) || !close( f._r2, 0.63 ) ) {
            return fail( "slope 0.0105, intercept -20.981, R2 0.63: " + f._slope + " " + f._intercept + " " + f._r2 );
        }
        if ( !close( f._rate, 0.0105 ) || !close( f._root_date, 20.981 / 0.0105 ) ) {
            return fail( "rate 0.0105; the root by the line at 1998.190: " + f._rate + " " + f._root_date );
        }
        // E: the line at 2004 is -20.981 + 21.042 = 0.061, and E is 0.05
        for( final ClockPlot.Point p : d._points ) {
            if ( "E".equals( p._node.getName() ) && ( Math.abs( p._residual - ( -0.011 ) ) > 1e-12 ) ) {
                return fail( "E lies 0.011 under the line: " + p._residual );
            }
        }
        // the view of Y: its own line, and its divergence is its path from the TREE's root
        final Phylogeny view = TreePanelUtil.subTree( named( t, "Y" ), t );
        final ClockPlot.Data c = ClockPlot.data( t, time( t ), view );
        if ( ( c == null ) || ( c._points.size() != 3 ) || ( c._root != view.getRoot() ) || !close( c._root_div, 0.01 ) || c.rootDated() ) {
            return fail( "the view of Y: three tips, Y 0.01 from the root and undated" );
        }
        final ClockPlot.Fit cf = c._fit;
        if ( ( cf == null ) || ( cf._n != 3 ) || !close( cf._slope, 0.01 ) || !close( cf._intercept, -19.98 ) || !close( cf._r2, 0.25 ) ) {
            return fail( "the clade's line: slope 0.01, intercept -19.98, R2 0.25" );
        }
        if ( !close( cf._root_date, 1999 ) ) {
            return fail( "the clade's root by the line at 1999: " + cf._root_date );
        }
        // the divergence's unit is the tree's own where it states one; the first kind's is subs/site
        t.setDistanceUnit( "substitutions/site" );
        if ( !"substitutions/site".equals( data( t )._div_unit ) ) {
            return fail( "a stated branch length unit is the divergence's unit" );
        }
        t.setDistanceUnit( "  " );
        if ( data( t )._div_unit != null ) {
            return fail( "a blank unit is none" );
        }
        final ClockPlot.Data k1 = data( build() );
        if ( !"subs/site".equals( k1._div_unit ) || !k1._calendar || ( k1._from_names != 0 ) || ( k1._ambiguous_names != 0 ) ) {
            return fail( "kind 1: subs/site, calendar years by the unit, no names read" );
        }
        return true;
    }

    /** Tips dated by AGE (a geologic unit, no dated pair to read the direction off): time still runs toward the
     *  youngest tip, so the slope is negative and the rate is reported positive. A..E at 4, 3, 2, 1, 0 million
     *  years: 2000..2004 reflected, so the same sums with Sxy's sign flipped; intercept 0.04 + 0.0105 x 2 = 0.061. */
    private static boolean testRootToTipAges() {
        final Phylogeny t = rootToTip();
        final String[] names = { "A", "B", "C", "D", "E" };
        for( int i = 0; i < names.length; ++i ) {
            final Date d = new Date();
            d.setValue( new BigDecimal( String.valueOf( 4 - i ) ) );
            d.setUnit( "mya" );
            named( t, names[ i ] ).getNodeData().setDate( d );
        }
        final ClockPlot.Data d = data( t );
        if ( ( d == null ) || ( d._kind != ClockPlot.Kind.ROOT_TO_TIP ) || d._forward || !"mya".equals( d._unit ) || d._calendar ) {
            return fail( "ages: the dates run the other way, and are no calendar" );
        }
        final ClockPlot.Fit f = d._fit;
        if ( !close( f._slope, -0.0105 ) || !close( f._rate, 0.0105 ) || !close( f._r2, 0.63 ) || !close( f._intercept, 0.061 ) ) {
            return fail( "slope -0.0105, rate +0.0105, R2 0.63, intercept 0.061: " + f._slope + " " + f._rate + " " + f._r2 + " " + f._intercept );
        }
        if ( !close( f._root_date, 0.061 / 0.0105 ) ) {
            return fail( "the root by the line at 5.81 Ma: " + f._root_date );
        }
        return true;
    }

    /** Christian, 2026-10-03: a tip dated only by its NAME counts (the class comment's numbers). */
    private static boolean testRootToTipNames() {
        final Phylogeny t = parse( "((A_2000:0.01,B|2001-07-02:0.03)X:0.005,(C_2002:0.02,D_2003:0.06,E_2004:0.04)Y:0.01);", ".nwk" );
        if ( ( t == null ) || ( t.getNumberOfExternalNodes() != 5 ) ) {
            return fail( "the names fixture could not be read" );
        }
        if ( kind( t ) != ClockPlot.Kind.ROOT_TO_TIP ) {
            return fail( "dates in the tip names: a root-to-tip plot, got " + kind( t ) );
        }
        final ClockPlot.Data d = data( t );
        if ( ( d._from_names != 5 ) || ( d._ambiguous_names != 0 ) || !d._calendar || !d._forward || ( d._unit != null ) ) {
            return fail( "five tips dated by their names, none ambiguous, calendar years, no stated unit: " + d._from_names + " "
                    + d._ambiguous_names + " " + d._calendar + " " + d._unit );
        }
        final java.util.List<String> dates = new java.util.ArrayList<String>();
        for( final ClockPlot.Point p : d._points ) {
            dates.add( p._node.getName().charAt( 0 ) + " " + ClockPlot.plain( p._date ) );
        }
        java.util.Collections.sort( dates );
        if ( !dates.toString().equals( "[A 2000.5, B 2001.5, C 2002.5, D 2003.5, E 2004.5]" ) ) {
            return fail( "a year alone is mid-year; 2 July 2001 is mid-year exactly: " + dates );
        }
        final ClockPlot.Fit f = d._fit;
        if ( ( f == null ) || !close( f._slope, 0.0105 ) || !close( f._r2, 0.63 ) || !close( f._intercept, -20.98625 )
                || !close( f._root_date, 20.98625 / 0.0105 ) ) {
            return fail( "slope 0.0105, R2 0.63, intercept -20.98625, the root by the line at 1998.690: " + f._intercept + " " + f._root_date );
        }
        // a date the file states wins over the one in the name
        final Phylogeny stated = parse( "#NEXUS\nBegin trees;\ntree T = [&R] ((A_2000[&num_date=1999]:0.01,B_2001:0.03)X:0.005,"
                + "(C_2002:0.02,D_2003:0.06,E_2004:0.04)Y:0.01);\nEnd;\n", ".nex" );
        final ClockPlot.Data s = data( stated );
        ClockPlot.Point a = null;
        for( final ClockPlot.Point p : s._points ) {
            if ( "A_2000".equals( p._node.getName() ) ) {
                a = p;
            }
        }
        if ( ( s._from_names != 4 ) || ( a == null ) || ( a._date != 1999 ) ) {
            return fail( "A: the stated 1999, not the name's 2000.5; the other four by their names: " + s._from_names + " " + ( a == null ? "-" : a._date ) );
        }
        // an ancestor's name dates nothing
        final Phylogeny xn = parse( "((A_2000:0.01,B_2001:0.03)X_1999:0.005,(C_2002:0.02,D_2003:0.06,E_2004:0.04)Y:0.01);", ".nwk" );
        if ( ( kind( xn ) != ClockPlot.Kind.ROOT_TO_TIP ) || ( data( xn )._points.size() != 5 ) ) {
            return fail( "an ancestor named with a year is not dated by it" );
        }
        // one tip with no date in the file and none in its name: refused, not dropped
        final Phylogeny gap = parse( "((A_2000:0.01,B_2001:0.03)X:0.005,(C:0.02,D_2003:0.06,E_2004:0.04)Y:0.01);", ".nwk" );
        if ( ClockPlot.isRootToTipTree( gap ) || ( kind( gap ) != null ) ) {
            return fail( "one tip undated by file and name: no plot" );
        }
        // a day/month-ambiguous numeric date is read day first and counted: 03/04/2002 is 3 April, day 93 of a
        // non-leap year
        final Phylogeny amb = parse( "((A_2000:0.01,B_2001:0.03)X:0.005,(C_03/04/2002:0.02,D_2003:0.06,E_2004:0.04)Y:0.01);", ".nwk" );
        final ClockPlot.Data ad = data( amb );
        ClockPlot.Point c = null;
        for( final ClockPlot.Point p : ad._points ) {
            if ( "C_03/04/2002".equals( p._node.getName() ) ) {
                c = p;
            }
        }
        if ( ( ad._from_names != 5 ) || ( ad._ambiguous_names != 1 ) || ( c == null ) || !close( c._date, 2002 + ( 92.5 / 365 ) ) ) {
            return fail( "one ambiguous name, read day first: 3 April 2002: " + ad._from_names + " " + ad._ambiguous_names
                    + " " + ( c == null ? "-" : c._date ) );
        }
        // the one-date check reads the same dates: every name the same year, no plot; the kind stands
        final Phylogeny one = parse( "((A_2000:0.01,B_2000:0.03)X:0.005,(C_2000:0.02,D_2000:0.06,E_2000:0.04)Y:0.01);", ".nwk" );
        if ( !ClockPlot.isRootToTipTree( one ) || ( kind( one ) != null ) ) {
            return fail( "every name on one year: the kind, and no plot" );
        }
        // ...and stated dates that differ make the line though every name says one year: the file's date wins there too
        final Phylogeny differ = parse( "#NEXUS\nBegin trees;\ntree T = [&R] ((A_2000[&num_date=2000]:0.01,B_2000[&num_date=2001]:0.03)X:0.005,"
                + "(C_2000[&num_date=2002]:0.02,D_2000[&num_date=2003]:0.06,E_2000[&num_date=2004]:0.04)Y:0.01);\nEnd;\n", ".nex" );
        if ( ( kind( differ ) != ClockPlot.Kind.ROOT_TO_TIP ) || ( data( differ )._from_names != 0 ) ) {
            return fail( "every name one year, the stated dates five: a plot, by the file" );
        }
        return true;
    }

    /** A plain tree whose tips get their dates from Extract Dates from Labels is the same kind before (by the names)
     *  and after (by the file, which then wins), with the same plot either way. */
    private static boolean testRootToTipFromLabels() {
        final Phylogeny phy = read( "date-in-labels.xml" );
        if ( phy == null ) {
            return fail( "date-in-labels.xml could not be read" );
        }
        if ( ( kind( phy ) != ClockPlot.Kind.ROOT_TO_TIP ) || ( data( phy )._from_names != 10 ) ) {
            return fail( "date-in-labels.xml: ten tips dated by their names, a root-to-tip plot" );
        }
        final ClockPlot.Fit by_names = data( phy )._fit;
        for( final TipDateExtractor.TipDate row : TipDateExtractor.preview( phy, TipDateExtractor.DayMonthOrder.DAY_FIRST ) ) {
            if ( !TipDateExtractor.applyToNode( row.node(), row.match() ) ) {
                return fail( "fixture: " + row.label() + " has a date to extract" );
            }
        }
        final ClockPlot.Data after = data( phy );
        if ( ( after == null ) || ( after._kind != ClockPlot.Kind.ROOT_TO_TIP ) || ( after._from_names != 0 ) || !after._calendar ) {
            return fail( "dates written onto the tips: still a root-to-tip plot, now by the file" );
        }
        // the extractor rounds a decimal year to four places: the line it fits is the same to that
        if ( !close( after._fit._slope, by_names._slope, 1e-3 ) || !close( after._fit._root_date, by_names._root_date, 1e-5 ) ) {
            return fail( "the same line by the file as by the names: " + after._fit._slope + " / " + by_names._slope );
        }
        return true;
    }

    /**
     * JOINT CROSS-RUN with Archaeopteryx.js (2d7fabe): its {@code docs/data/swH1-HA1.nexus} (fifty swine H1 HA1
     * strains named ".../YYYY|clade", lengths stated, no unit, no date annotation anywhere) is a root-to-tip tree
     * dated by its names. The numbers are THEIRS, from their run of their code, held to 1e-12; ours are sent back.
     */
    private static boolean testSwineH1CrossRun() {
        final Phylogeny t = readFile( new File( System.getProperty( "user.dir" ), "forester/test_data/clock/swH1-HA1.nexus" ) );
        if ( t == null ) {
            return fail( "forester/test_data/clock/swH1-HA1.nexus could not be read" );
        }
        final ClockPlot.Data d = data( t );
        if ( ( d == null ) || ( d._kind != ClockPlot.Kind.ROOT_TO_TIP ) || ( d._fit == null ) ) {
            return fail( "swH1-HA1.nexus: a root-to-tip plot with a line" );
        }
        final ClockPlot.Fit f = d._fit;
        if ( ( f._n != 50 ) || ( d._from_names != 50 ) || ( d._ambiguous_names != 0 ) || !d._calendar || ( d._unit != null )
                || ( d._div_unit != null ) || d.rootDated() || ( d._root_div != 0 ) ) {
            return fail( "swH1-HA1.nexus: 50 tips by their names, calendar, no units, an undated root at 0: " + f._n + " "
                    + d._from_names + " " + d._ambiguous_names + " " + d._calendar + " " + d._unit + " " + d._div_unit );
        }
        if ( !close( f._slope, 0.00434714238019987, 1e-12 ) || !close( f._rate, 0.00434714238019987, 1e-12 )
                || !close( f._intercept, -8.493052180996681, 1e-12 ) || !close( f._r2, 0.5915897595418261, 1e-12 )
                || !close( f._root_date, 1953.7092273950764, 1e-12 ) ) {
            return fail( "swH1-HA1.nexus against Archaeopteryx.js: slope " + f._slope + " intercept " + f._intercept + " R2 " + f._r2
                    + " root by the line " + f._root_date );
        }
        return true;
    }

    /** forester/demo/root-to-tip.xml and its refused twin, read as File > Open reads them: the numbers the README row
     *  quotes, computed by hand in DemoTreeGenerator. */
    private static boolean testRootToTipDemoPair() {
        final Phylogeny phy = read( "root-to-tip.xml" );
        final Phylogeny twin = read( "root-to-tip-tip-undated.xml" );
        if ( ( phy == null ) || ( twin == null ) ) {
            return fail( "the root-to-tip demo pair could not be read" );
        }
        final ClockPlot.Data d = data( phy );
        if ( ( d == null ) || ( d._kind != ClockPlot.Kind.ROOT_TO_TIP ) || ( d.tipCount() != 16 ) || ( d._fit._n != 16 ) || d.rootDated()
                || !d._calendar || ( d._from_names != 0 ) || !"substitutions/site".equals( d._div_unit ) ) {
            return fail( "root-to-tip.xml: 16 tips in the line, no root date, calendar years, substitutions/site" );
        }
        if ( !close( d._fit._slope, 0.0039386126, 1e-7 ) || !ClockPlot.number( d._fit._rate, 3 ).equals( "0.00394" )
                || !ClockPlot.number( d._fit._r2, 3 ).equals( "0.975" ) || !ClockPlot.date( d._fit._root_date, true ).equals( "2017.94" ) ) {
            return fail( "rate 0.00394, R2 0.975, root by the line 2017.94: " + d._fit._rate + " " + d._fit._r2 + " " + d._fit._root_date );
        }
        // the tip furthest off the line
        ClockPlot.Point worst = null;
        for( final ClockPlot.Point p : d._points ) {
            if ( ( worst == null ) || ( Math.abs( p._residual ) > Math.abs( worst._residual ) ) ) {
                worst = p;
            }
        }
        if ( !"R-1".equals( worst._node.getName() ) || !ClockPlot.signed( worst._residual ).equals( "+0.00159" ) ) {
            return fail( "the tip furthest off the line: R-1, +0.00159: " + worst._node.getName() + " " + worst._residual );
        }
        // the clade of S, 0.0128 from the root: its own line
        final ClockPlot.Data c = ClockPlot.data( phy, time( phy ), TreePanelUtil.subTree( named( phy, "S" ), phy ) );
        if ( ( c == null ) || ( c._points.size() != 4 ) || !close( c._root_div, 0.0128 ) || c.rootDated()
                || !close( c._fit._slope, 0.0042135922, 1e-7 ) || !ClockPlot.number( c._fit._r2, 3 ).equals( "0.933" )
                || !ClockPlot.date( c._fit._root_date, true ).equals( "2021.37" ) ) {
            return fail( "the clade of S: four tips, 0.0128 from the root, R2 0.933, its top by the line at 2021.37" );
        }
        // the twin: S-3 undated, the plot refused
        if ( named( twin, "S-3" ).getNodeData().isHasDate() || ( ClockPlot.rootToTipDate( named( twin, "S-3" ) ) != null ) ) {
            return fail( "fixture: S-3 has no date in the file and none in its name" );
        }
        if ( ClockPlot.isRootToTipTree( twin ) || ( kind( twin ) != null ) || ( data( twin ) != null ) ) {
            return fail( "root-to-tip-tip-undated.xml: refused the plot" );
        }
        return true;
    }

    private ClockPlotTest() {
    }
}
