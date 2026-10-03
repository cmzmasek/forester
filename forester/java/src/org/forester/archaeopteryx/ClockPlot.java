// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * The CLOCK PLOT, its pure half: every node as a point, its DATE against its DIVERGENCE from the root, and a straight
 * line through the tips. The slope is the rate the tips accumulated divergence at, and where the line comes down to
 * the root's divergence is the date it puts the root on: the root-to-tip regression. Ported from Archaeopteryx.js
 * 3.23.0 ({@code forester.clockPlotKind}, {@code clockPlotData}, {@code clockRegression},
 * {@code calendarTickMonths}); Christian, 2026-10-02: "desktop should copy what JS did for the clock plot". The panel
 * is {@link ClockPlotPanel}.
 * <p>
 * <b>Which trees have one.</b> A tree offered Time | Div ({@link BranchLengthLayout#isApplicable}): each node states
 * a date and a divergence from the root, recorded on it or summed from each branch's time x its clock rate
 * ({@link BranchLengthLayout#divergenceFromRoot}). Which of the two layouts is on screen does not matter. And its LINE
 * can be drawn: three tips or more, not all on one date -- a tree sampled at one moment has no slope to estimate.
 * Asked of the WHOLE tree, never of the view: in the view of a clade too small or too uniform for a line the points
 * are drawn and the line is not. (A second kind, a divergence tree whose tips alone are dated, is described by
 * Archaeopteryx.js and built by neither side.)
 * <p>
 * <b>The fit</b> is ordinary least squares over the TIPS alone. A tip's date is an observation; an ancestor's was
 * inferred, usually with a clock, and counting it would have the estimate confirm itself. The line is not forced
 * through the root, so the date it reaches the root's divergence at can be held against the date the tree states.
 */
final class ClockPlot {

    /** The fewest tips a line is drawn through. */
    static final int MIN_TIPS = 3;

    private ClockPlot() {
    }

    /** One node of the plot. */
    static final class Point {

        final PhylogenyNode _node;
        final double        _date;
        final double        _div;
        final boolean       _tip;
        /** A tip's divergence less the line's at its date; null for an ancestor, or where there is no line. */
        Double              _residual;

        Point( final PhylogenyNode node, final double date, final double div, final boolean tip ) {
            _node = node;
            _date = date;
            _div = div;
            _tip = tip;
        }
    }

    /** Ordinary least squares of y on x. */
    static final class Fit {

        final int    _n;
        final double _slope;
        final double _intercept;
        /** Null where every y is the same: no variance to explain. */
        final Double _r2;
        final double _mean_x;
        final double _mean_y;
        /** The slope in the direction time runs: divergence per unit of time, negative where the tips' divergence
         *  falls with time. Set by {@link ClockPlot#data}. */
        double       _rate;
        /** The date the line reaches the view root's divergence at; null unless the rate is positive. */
        Double       _root_date;

        Fit( final int n, final double slope, final double intercept, final Double r2, final double mean_x,
             final double mean_y ) {
            _n = n;
            _slope = slope;
            _intercept = intercept;
            _r2 = r2;
            _mean_x = mean_x;
            _mean_y = mean_y;
        }

        double at( final double x ) {
            return _intercept + ( _slope * x );
        }
    }

    /** The clock plot of a tree, or of the clade on view. */
    static final class Data {

        /** Whether the dates increase toward the tips (calendar dates); false for ages. */
        final boolean       _forward;
        /** The dates' unit as the tree states it, or null. */
        final String        _unit;
        /** Whether the divergence is time x clock rate (a BEAST clock-model tree), not recorded. */
        final boolean       _from_rates;
        /** Every node of the view, in preorder. */
        final List<Point>   _points;
        /** The top of the view. */
        final PhylogenyNode _root;
        final double        _root_date;
        final double        _root_div;
        /** The line through the view's tips, or null where it has none (fewer than three, or one date). */
        final Fit           _fit;

        Data( final boolean forward, final String unit, final boolean from_rates, final List<Point> points,
              final PhylogenyNode root, final double root_date, final double root_div, final Fit fit ) {
            _forward = forward;
            _unit = unit;
            _from_rates = from_rates;
            _points = Collections.unmodifiableList( points );
            _root = root;
            _root_date = root_date;
            _root_div = root_div;
            _fit = fit;
        }

        int tipCount() {
            int n = 0;
            for( final Point p : _points ) {
                if ( p._tip ) {
                    ++n;
                }
            }
            return n;
        }
    }

    /**
     * Whether the tree has a clock plot (see the class comment).
     *
     * @param whole the WHOLE tree a tab holds
     * @param time the lengths it has in time: on screen, or kept
     */
    static boolean isOffered( final Phylogeny whole, final BranchLengthLayout.TimeLengths time ) {
        return ( whole != null ) && !whole.isEmpty() && ( time != null ) && BranchLengthLayout.isApplicable( whole, time )
                && tipsMakeALine( whole );
    }

    /** The half of {@link #isOffered} that is the plot's own: three tips or more, not all on one date. (The other
     *  half, Time | Div, a tab caches.) */
    static boolean tipsMakeALine( final Phylogeny whole ) {
        int tips = 0;
        Double first = null;
        boolean differ = false;
        // preorder, not the external-node iterator: that one reads a cache a subtree view of the tree leaves stale
        for( final PhylogenyNodeIterator it = whole.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() ) {
                continue;
            }
            final Double d = BranchLengthLayout.dateValue( n );
            ++tips;
            if ( first == null ) {
                first = d;
            }
            else if ( ( d != null ) && ( d.doubleValue() != first.doubleValue() ) ) {
                differ = true;
            }
        }
        return ( tips >= MIN_TIPS ) && differ;
    }

    /**
     * Ordinary least squares of y on x, computed about the means: a date is near 2020 and its spread a few months,
     * and sums of raw squares would spend every digit on the 2020.
     *
     * @return null with fewer than three points, as many x as y not given, or every x the same
     */
    static Fit regression( final double[] xs, final double[] ys ) {
        final int n = xs.length;
        if ( ( n < MIN_TIPS ) || ( ys.length != n ) ) {
            return null;
        }
        double mx = 0;
        double my = 0;
        for( int i = 0; i < n; ++i ) {
            mx += xs[ i ];
            my += ys[ i ];
        }
        mx /= n;
        my /= n;
        double sxx = 0;
        double sxy = 0;
        double syy = 0;
        for( int i = 0; i < n; ++i ) {
            final double dx = xs[ i ] - mx;
            final double dy = ys[ i ] - my;
            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
        }
        if ( !( sxx > 0 ) || !Double.isFinite( sxx ) || !Double.isFinite( sxy ) || !Double.isFinite( syy ) ) {
            return null;
        }
        final double slope = sxy / sxx;
        return new Fit( n, slope, my - ( slope * mx ), ( syy > 0 ) ? Double.valueOf( ( sxy * sxy ) / ( sxx * syy ) )
                : null, mx, my );
    }

    /**
     * The clock plot of the tree, or of the clade on view.
     *
     * @param whole the WHOLE tree a tab holds: the rule is asked of it, and the divergence from ITS root is plotted
     * @param time the lengths it has in time: on screen, or kept
     * @param view the tree on display (the whole tree, or a subtree view of it); the whole tree where null
     * @return null where the tree has no clock plot
     */
    static Data data( final Phylogeny whole, final BranchLengthLayout.TimeLengths time, final Phylogeny view ) {
        if ( !isOffered( whole, time ) ) {
            return null;
        }
        final Map<Long, Double> div_of = BranchLengthLayout.divergenceFromRoot( whole, time );
        final Phylogeny shown = ( ( view == null ) || view.isEmpty() ) ? whole : view;
        final List<Point> points = new ArrayList<Point>();
        final List<Double> xs = new ArrayList<Double>();
        final List<Double> ys = new ArrayList<Double>();
        for( final PhylogenyNodeIterator it = shown.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final Double date = BranchLengthLayout.dateValue( n );
            final Double div = div_of.get( n.getId() );
            if ( ( date == null ) || ( div == null ) ) {
                continue; // a guard, not a rule: every node of a tree offered Time | Div states both
            }
            final Point p = new Point( n, date.doubleValue(), div.doubleValue(), n.isExternal() );
            points.add( p );
            if ( p._tip ) {
                xs.add( date );
                ys.add( div );
            }
        }
        final boolean forward = BranchLengthLayout.datesIncreaseTowardTips( whole );
        final PhylogenyNode top = shown.getRoot();
        final Double top_date = BranchLengthLayout.dateValue( top );
        final Double top_div = div_of.get( top.getId() );
        final Fit fit = regression( unbox( xs ), unbox( ys ) );
        if ( fit != null ) {
            fit._rate = forward ? fit._slope : -fit._slope;
            fit._root_date = ( ( fit._rate > 0 ) && ( top_div != null ) )
                    ? Double.valueOf( ( top_div.doubleValue() - fit._intercept ) / fit._slope ) : null;
            for( final Point p : points ) {
                if ( p._tip ) {
                    p._residual = Double.valueOf( p._div - fit.at( p._date ) );
                }
            }
        }
        return new Data( forward,
                         AptxUtil.timeTreeUnit( whole ),
                         BranchLengthLayout.divergenceSource( whole ) == BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE,
                         points,
                         top,
                         ( top_date == null ) ? Double.NaN : top_date.doubleValue(),
                         ( top_div == null ) ? Double.NaN : top_div.doubleValue(),
                         fit );
    }

    private static double[] unbox( final List<Double> values ) {
        final double[] a = new double[ values.size() ];
        for( int i = 0; i < a.length; ++i ) {
            a[ i ] = values.get( i ).doubleValue();
        }
        return a;
    }

    // ---- ticks ---------------------------------------------------------------------------------------------------

    /** A month tick: the decimal year at which that month begins. */
    static final class MonthTick {

        final double _value;
        final int    _year;
        /** 1-12 */
        final int    _month;

        MonthTick( final double value, final int year, final int month ) {
            _value = value;
            _year = year;
            _month = month;
        }
    }

    private static int yearLength( final int y ) {
        return ( ( ( ( y % 4 ) == 0 ) && ( ( y % 100 ) != 0 ) ) || ( ( y % 400 ) == 0 ) ) ? 366 : 365;
    }

    private static int monthLength( final int y, final int m ) {
        if ( m == 2 ) {
            return ( yearLength( y ) == 366 ) ? 29 : 28;
        }
        return ( ( m == 4 ) || ( m == 6 ) || ( m == 9 ) || ( m == 11 ) ) ? 30 : 31;
    }

    private static int dayOfYear( final int y, final int m, final int d ) {
        int doy = d;
        for( int i = 1; i < m; ++i ) {
            doy += monthLength( y, i );
        }
        return doy;
    }

    /**
     * Month ticks over [from, to] in decimal years, for a span too short for whole years (an outbreak sampled over
     * some months): the first of every 1st, 2nd, 3rd, 6th or 12th month, the smallest step leaving six ticks or
     * fewer; none where even a tick a year is too many. As Archaeopteryx.js's {@code forester.calendarTickMonths}.
     */
    static List<MonthTick> calendarTickMonths( final double from, final double to ) {
        final double span = to - from;
        if ( !( span > 0 ) || !Double.isFinite( span ) || !Double.isFinite( from ) ) {
            return Collections.emptyList();
        }
        final int[] steps = { 1, 2, 3, 6, 12 };
        for( final int step : steps ) {
            final List<MonthTick> ticks = new ArrayList<MonthTick>();
            final long year = (long) Math.floor( from );
            // month index = year * 12 + (month - 1); ticks where it divides
            for( long k = year * 12; ticks.size() <= 6; ++k ) {
                final int y = (int) Math.floorDiv( k, 12L );
                final int m = (int) ( k - ( y * 12L ) ) + 1;
                final double v = y + ( ( dayOfYear( y, m, 1 ) - 1 ) / (double) yearLength( y ) );
                if ( v > ( to + 1e-9 ) ) {
                    break;
                }
                if ( ( v >= ( from - 1e-9 ) ) && ( Math.floorMod( k, (long) step ) == 0 ) ) {
                    ticks.add( new MonthTick( v, y, m ) );
                }
            }
            if ( ticks.size() <= 6 ) {
                return ticks;
            }
        }
        return Collections.emptyList();
    }

    /** The smallest of 1, 2, 5, 10 x a power of ten at or above {@code target} (Archaeopteryx.js's
     *  {@code forester.niceAxisStep}). */
    static double niceAxisStep( final double target ) {
        if ( !( target > 0 ) || !Double.isFinite( target ) ) {
            return 1;
        }
        final double mag = Math.pow( 10, Math.floor( Math.log10( target ) ) );
        for( final double c : new double[] { 1, 2, 5, 10 } ) {
            final double s = c * mag;
            if ( s >= ( target - 1e-12 ) ) {
                return s;
            }
        }
        return 10 * mag;
    }

    /** Whole-year ticks over [from, to], about seven of them (Archaeopteryx.js's {@code forester.calendarTickYears}). */
    static List<Integer> calendarTickYears( final double from, final double to ) {
        final double span = to - from;
        if ( !( span > 0 ) || !Double.isFinite( span ) || !Double.isFinite( from ) ) {
            return Collections.emptyList();
        }
        final long step = Math.max( 1, Math.round( niceAxisStep( span / 7 ) ) );
        final List<Integer> years = new ArrayList<Integer>();
        for( long y = (long) Math.ceil( from / step ) * step; y <= ( to + 1e-9 ); y += step ) {
            years.add( Integer.valueOf( (int) y ) );
        }
        return years;
    }

    /**
     * About {@code count} round values over [start, stop], as d3's {@code d3.ticks} gives them (d3-array 3:
     * {@code tickSpec}), which is what Archaeopteryx.js marks the plot's axes with.
     */
    static double[] ticks( final double start, final double stop, final int count ) {
        if ( !( count > 0 ) || !Double.isFinite( start ) || !Double.isFinite( stop ) ) {
            return new double[ 0 ];
        }
        if ( start == stop ) {
            return new double[] { start };
        }
        final boolean reverse = stop < start;
        final double lo = reverse ? stop : start;
        final double hi = reverse ? start : stop;
        final double step = ( hi - lo ) / count;
        final double power = Math.floor( Math.log10( step ) );
        final double error = step / Math.pow( 10, power );
        final double factor = ( error >= Math.sqrt( 50 ) ) ? 10 : ( error >= Math.sqrt( 10 ) ) ? 5
                : ( error >= Math.sqrt( 2 ) ) ? 2 : 1;
        double i1, i2, inc;
        if ( power < 0 ) {
            inc = Math.pow( 10, -power ) / factor;
            i1 = Math.round( lo * inc );
            i2 = Math.round( hi * inc );
            if ( ( i1 / inc ) < lo ) {
                ++i1;
            }
            if ( ( i2 / inc ) > hi ) {
                --i2;
            }
            inc = -inc;
        }
        else {
            inc = Math.pow( 10, power ) * factor;
            i1 = Math.round( lo / inc );
            i2 = Math.round( hi / inc );
            if ( ( i1 * inc ) < lo ) {
                ++i1;
            }
            if ( ( i2 * inc ) > hi ) {
                --i2;
            }
        }
        if ( ( i2 < i1 ) && ( 0.5 <= count ) && ( count < 2 ) ) {
            return ticks( start, stop, count * 2 );
        }
        if ( !( i2 >= i1 ) ) {
            return new double[ 0 ];
        }
        final int n = (int) ( i2 - i1 + 1 );
        final double[] t = new double[ n ];
        for( int i = 0; i < n; ++i ) {
            final double k = reverse ? ( i2 - i ) : ( i1 + i );
            t[ i ] = ( inc < 0 ) ? ( k / -inc ) : ( k * inc );
        }
        return t;
    }

    // ---- numbers as the panel prints them ------------------------------------------------------------------------

    /** A number at {@code digits} significant digits (3 by default); very small or very large in exponent form, which
     *  is how a rate is usually written. As Archaeopteryx.js's {@code clockNumber}: JS's {@code toExponential} and
     *  {@code String(Number(toPrecision))}. */
    static String number( final double v, final int digits ) {
        if ( !Double.isFinite( v ) ) {
            return "";
        }
        final int sig = ( digits > 0 ) ? digits : 3;
        final double a = Math.abs( v );
        if ( ( a != 0 ) && ( ( a < 1e-4 ) || ( a >= 1e7 ) ) ) {
            final String s = String.format( java.util.Locale.ROOT, "%." + ( sig - 1 ) + "e", v );
            // JS writes 1.10e-3 as "1.10e-3", Java as "1.10e-03"
            final int e = s.indexOf( 'e' );
            final String mantissa = s.substring( 0, e );
            final int exp = Integer.parseInt( s.substring( e + 1 ) );
            return mantissa + "e" + ( ( exp < 0 ) ? "-" : "+" ) + Math.abs( exp );
        }
        return plain( new java.math.BigDecimal( v ).round( new java.math.MathContext( sig ) ).doubleValue() );
    }

    /** A double the way JavaScript's {@code String(number)} writes one in the plain range: no trailing ".0". */
    static String plain( final double v ) {
        if ( v == Math.rint( v ) && ( Math.abs( v ) < 1e15 ) ) {
            return Long.toString( (long) v );
        }
        return new java.math.BigDecimal( Double.toString( v ) ).stripTrailingZeros().toPlainString();
    }

    /** A date of the plot: a calendar year to two decimals (about four days), anything else at five digits. */
    static String date( final double v, final boolean calendar ) {
        if ( !Double.isFinite( v ) ) {
            return "";
        }
        return calendar ? plain( Math.round( v * 100.0 ) / 100.0 ) : number( v, 5 );
    }

    /** A residual with its sign. */
    static String signed( final double v ) {
        return ( ( v > 0 ) ? "+" : "" ) + number( v, 3 );
    }
}
