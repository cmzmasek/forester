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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.forester.archaeopteryx.tools.TipDateExtractor;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/**
 * The CLOCK PLOT, its pure half: every node as a point, its DATE against its DIVERGENCE from the root, and a straight
 * line through the tips. The slope is the rate the tips accumulated divergence at, and where the line comes down to
 * the root's divergence is the date it puts the root on: the root-to-tip regression. Ported from Archaeopteryx.js
 * 3.23.0 ({@code forester.clockPlotKind}, {@code clockPlotData}, {@code clockRegression},
 * {@code calendarTickMonths}); Christian, 2026-10-02: "desktop should copy what JS did for the clock plot". The panel
 * is {@link ClockPlotPanel}.
 * <p>
 * <b>Which trees have one</b> ({@link Kind}, THE RULE, joint with Archaeopteryx.js; the second kind built there in
 * e22361c on Christian's word, 2026-10-03, and ported here):
 * <ol>
 * <li>DIVERGENCE: a tree offered Time | Div ({@link BranchLengthLayout#isApplicable}): each node states a date and a
 * divergence from the root, recorded on it or summed from each branch's time x its clock rate
 * ({@link BranchLengthLayout#divergenceFromRoot}). Every node is a point; which of the two layouts is on screen does
 * not matter.</li>
 * <li>ROOT_TO_TIP: a tree whose branch lengths are its divergence and whose TIPS alone are dated, with no date on
 * any ancestor (it is not a time tree): a tip's divergence is the sum of the stated lengths from the root. Tips
 * only. This is the tree a clock plot is wanted for most -- is there a signal worth dating? -- and the one whose
 * answer depends on its root: re-rooting moves every point ({@link #isRootToTipTree}). A tip is dated by what the
 * file states for it, or else by the date its NAME carries (Christian, 2026-10-03: names count; JS 2d7fabe;
 * {@link #rootToTipDate}).</li>
 * </ol>
 * The first kind is asked first; the two cannot both hold. In both the plot is offered when its LINE can be drawn:
 * three tips or more, not all on one date -- a tree sampled at one moment has no slope to estimate. Asked of the
 * WHOLE tree, never of the view: in the view of a clade too small or too uniform for a line the points are drawn
 * and the line is not.
 * <p>
 * <b>The fit</b> is ordinary least squares over the TIPS alone. A tip's date is an observation; an ancestor's was
 * inferred, usually with a clock, and counting it would have the estimate confirm itself. The line is not forced
 * through the root, so the date it reaches the root's divergence at can be held against the date the tree states --
 * where the tree states one; a root-to-tip tree does not, and there the line's date is the only estimate.
 */
final class ClockPlot {

    /** The fewest tips a line is drawn through. */
    static final int MIN_TIPS = 3;

    private ClockPlot() {
    }

    /** Which clock plot a tree has (the class comment). */
    enum Kind {
        /** A tree with Time | Div: every node dated and with a divergence from the root. */
        DIVERGENCE,
        /** A divergence tree whose tips alone are dated: a tip's divergence is its path from the root. */
        ROOT_TO_TIP
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

        final Kind          _kind;
        /** Whether the dates increase toward the tips (calendar dates); false for ages. */
        final boolean       _forward;
        /** The dates' unit as the tree states it, or null. */
        final String        _unit;
        /** The divergence's unit: "subs/site" on a tree with Time | Div, a root-to-tip tree's stated branch length
         *  unit, or null where it states none. */
        final String        _div_unit;
        /** Whether the divergence is time x clock rate (a BEAST clock-model tree), not recorded. */
        final boolean       _from_rates;
        /** Whether the dates are calendar years: their unit says so, or they were read off the tips' names, which
         *  carry no other kind of date. The model's own answer, as Archaeopteryx.js's: never the axis on screen. */
        final boolean       _calendar;
        /** How many tips are dated by their names (a root-to-tip tree; 0 otherwise). */
        final int           _from_names;
        /** How many of those names carry a day/month-ambiguous numeric date, read day first. */
        final int           _ambiguous_names;
        /** Every node of the view, in preorder; on a root-to-tip tree the tips only, no ancestor being dated. */
        final List<Point>   _points;
        /** The top of the view. */
        final PhylogenyNode _root;
        /** NaN where the top of the view states no date: on a root-to-tip tree always. */
        final double        _root_date;
        final double        _root_div;
        /** The line through the view's tips, or null where it has none (fewer than three, or one date). */
        final Fit           _fit;

        Data( final Kind kind, final boolean forward, final String unit, final String div_unit, final boolean from_rates,
              final boolean calendar, final int from_names, final int ambiguous_names, final List<Point> points,
              final PhylogenyNode root, final double root_date, final double root_div, final Fit fit ) {
            _kind = kind;
            _forward = forward;
            _unit = unit;
            _div_unit = div_unit;
            _from_rates = from_rates;
            _calendar = calendar;
            _from_names = from_names;
            _ambiguous_names = ambiguous_names;
            _points = Collections.unmodifiableList( points );
            _root = root;
            _root_date = root_date;
            _root_div = root_div;
            _fit = fit;
        }

        /** Whether the top of the view states a date of its own, to hold the line's against. */
        boolean rootDated() {
            return !Double.isNaN( _root_date );
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

    /** A tip's date on a root-to-tip tree, and where it came from. */
    static final class TipDate {

        final double  _value;
        /** Read off the tip's name, no date being stated for it. */
        final boolean _from_name;
        /** A numeric date whose day and month could be either way round, read day first. */
        final boolean _ambiguous;

        TipDate( final double value, final boolean from_name, final boolean ambiguous ) {
            _value = value;
            _from_name = from_name;
            _ambiguous = ambiguous;
        }
    }

    /**
     * A tip's date on a root-to-tip tree: the one the file states, else the one its NAME carries, read by the
     * label-date parser both programs anchor BEAST heights with ({@link TipDateExtractor#parse}, as Archaeopteryx.js's
     * {@code forester.parseTipLabelDate}): the rightmost date in the label, a year alone taken as mid-year
     * (2020 -> 2020.5), a day/month-ambiguous numeric date read DAY FIRST and counted. Null for no date at all. A
     * stated date wins over the name; an ancestor is never dated by its name.
     */
    static TipDate rootToTipDate( final PhylogenyNode tip ) {
        final Double stated = BranchLengthLayout.dateValue( tip );
        if ( stated != null ) {
            return new TipDate( stated.doubleValue(), false, false );
        }
        final TipDateExtractor.DateMatch m = ForesterUtil.isEmpty( tip.getName() ) ? null
                : TipDateExtractor.parse( tip.getName(), TipDateExtractor.DayMonthOrder.DAY_FIRST );
        return ( m == null ) ? null : new TipDate( m.decimalYear(), true, m.ambiguous() );
    }

    /**
     * Kind 2 of the rule, before the tips are counted (Archaeopteryx.js's {@code isRootToTipTree}): every tip is
     * dated, by the file or by its name ({@link #rootToTipDate}); no ancestor states a date (the root included; a
     * name dates no ancestor); every branch states a length -- a length of 0, or a negative one, is a length the
     * file stated; one that is not there is not -- and no clock rate on every branch: a clock-model tree's lengths
     * are time, so it keeps its Time | Div plot or has none. Every tip must be dated and every branch must state a
     * length: a tip with either missing has no point, and a plot that drops tips silently answers a different
     * question than it was asked, so the tree is refused the plot, not plotted with tips missing.
     */
    static boolean isRootToTipTree( final Phylogeny whole ) {
        if ( ( whole == null ) || whole.isEmpty() ) {
            return false;
        }
        int branches = 0;
        for( final PhylogenyNodeIterator it = whole.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isExternal() ? ( rootToTipDate( n ) == null ) : ( BranchLengthLayout.dateValue( n ) != null ) ) {
                return false;
            }
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            ++branches;
            if ( n.getDistanceToParent() == PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                return false;
            }
        }
        return ( branches > 0 )
                && ( BranchLengthLayout.divergenceSource( whole ) != BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE );
    }

    /** Each node's distance from the TREE's root along the stated branch lengths, by node id: a root-to-tip tree's
     *  divergence (Archaeopteryx.js's {@code rootToTipDistances}). A branch without a length adds nothing -- a guard:
     *  {@link #isRootToTipTree} refuses such a tree. */
    static Map<Long, Double> rootToTipDistances( final Phylogeny whole ) {
        final Map<Long, Double> of = new HashMap<Long, Double>();
        for( final PhylogenyNodeIterator it = whole.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                of.put( n.getId(), Double.valueOf( 0.0 ) );
            }
            else {
                final double len = n.getDistanceToParent();
                of.put( n.getId(), Double.valueOf( of.get( n.getParent().getId() ).doubleValue()
                        + ( ( ( len != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) && Double.isFinite( len ) ) ? len : 0 ) ) );
            }
        }
        return of;
    }

    /**
     * Which clock plot the tree has, or null for none (the class comment): DIVERGENCE where it is offered Time | Div,
     * else ROOT_TO_TIP where its tips alone are dated over stated lengths -- and in either only where its tips make
     * a line. Asked of the WHOLE tree a tab holds, never of the view.
     *
     * @param time_div whether the tree is offered Time | Div (a tab caches that answer)
     */
    static Kind kind( final Phylogeny whole, final boolean time_div ) {
        if ( ( whole == null ) || whole.isEmpty() ) {
            return null;
        }
        final Kind kind = time_div ? Kind.DIVERGENCE : ( isRootToTipTree( whole ) ? Kind.ROOT_TO_TIP : null );
        return ( ( kind != null ) && tipsMakeALine( whole, kind == Kind.ROOT_TO_TIP ) ) ? kind : null;
    }

    /**
     * {@link #kind(Phylogeny, boolean)}, asking the tree itself whether it is offered Time | Div.
     *
     * @param time the lengths it has in time: on screen, or kept
     */
    static Kind kind( final Phylogeny whole, final BranchLengthLayout.TimeLengths time ) {
        if ( ( whole == null ) || whole.isEmpty() || ( time == null ) ) {
            return null;
        }
        return kind( whole, BranchLengthLayout.isApplicable( whole, time ) );
    }

    /**
     * Whether the tree has a clock plot of either kind.
     *
     * @param whole the WHOLE tree a tab holds
     * @param time the lengths it has in time: on screen, or kept
     */
    static boolean isOffered( final Phylogeny whole, final BranchLengthLayout.TimeLengths time ) {
        return kind( whole, time ) != null;
    }

    /** The half of {@link #kind} that is the line's own, the same for both kinds: three tips or more, not all on one
     *  date -- the dates the kind reads: on a root-to-tip tree ({@code names}) a tip's name dates it too. */
    static boolean tipsMakeALine( final Phylogeny whole, final boolean names ) {
        int tips = 0;
        Double first = null;
        boolean differ = false;
        // preorder, not the external-node iterator: that one reads a cache a subtree view of the tree leaves stale
        for( final PhylogenyNodeIterator it = whole.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() ) {
                continue;
            }
            final Double d;
            if ( names ) {
                final TipDate t = rootToTipDate( n );
                d = ( t == null ) ? null : Double.valueOf( t._value );
            }
            else {
                d = BranchLengthLayout.dateValue( n );
            }
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
        final Kind kind = kind( whole, time );
        if ( kind == null ) {
            return null;
        }
        final boolean root_to_tip = kind == Kind.ROOT_TO_TIP;
        // a root-to-tip tree's divergence is the path from the TREE's root, in the view of a clade too: the clade's
        // top then has one as well, and the line's root date is where it reaches that
        final Map<Long, Double> div_of = root_to_tip ? rootToTipDistances( whole )
                : BranchLengthLayout.divergenceFromRoot( whole, time );
        final Phylogeny shown = ( ( view == null ) || view.isEmpty() ) ? whole : view;
        final List<Point> points = new ArrayList<Point>();
        final List<Double> xs = new ArrayList<Double>();
        final List<Double> ys = new ArrayList<Double>();
        int from_names = 0;
        int ambiguous_names = 0;
        for( final PhylogenyNodeIterator it = shown.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            Double date = BranchLengthLayout.dateValue( n );
            if ( root_to_tip && n.isExternal() && ( date == null ) ) {
                final TipDate named = rootToTipDate( n ); // every tip has one: the kind says so
                date = ( named == null ) ? null : Double.valueOf( named._value );
                if ( ( named != null ) && named._from_name ) {
                    ++from_names;
                    if ( named._ambiguous ) {
                        ++ambiguous_names;
                    }
                }
            }
            final Double div = div_of.get( n.getId() );
            if ( ( date == null ) || ( div == null ) ) {
                continue; // an undated ancestor (root-to-tip); else a guard, not a rule: a tree offered Time | Div states both
            }
            final Point p = new Point( n, date.doubleValue(), div.doubleValue(), n.isExternal() );
            points.add( p );
            if ( p._tip ) {
                xs.add( date );
                ys.add( div );
            }
        }
        // a root-to-tip tree has no dated pair to read the direction off: its dates are calendar dates unless their
        // unit says ages (the first stated unit here; a majority of them in Archaeopteryx.js)
        final Options.TIME_AXIS_TYPE axis = AptxUtil.deriveTimeAxisType( whole );
        final boolean forward = root_to_tip ? ( axis != Options.TIME_AXIS_TYPE.GEOLOGIC )
                : BranchLengthLayout.datesIncreaseTowardTips( whole );
        // calendar years: the dates' unit says so, or they were read off the tips' names, which carry no other kind
        // of date
        final boolean calendar = forward && ( ( axis == Options.TIME_AXIS_TYPE.CALENDAR )
                || ( root_to_tip && ( axis == Options.TIME_AXIS_TYPE.NONE ) && ( from_names > 0 ) ) );
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
        final String length_unit = ForesterUtil.isEmpty( whole.getDistanceUnit() ) ? null : whole.getDistanceUnit().trim();
        return new Data( kind,
                         forward,
                         AptxUtil.timeTreeUnit( whole ),
                         root_to_tip ? ( ( ( length_unit != null ) && !length_unit.isEmpty() ) ? length_unit : null )
                                 : BranchLengthLayout.distanceUnit( BranchLengthLayout.MODE.DIVERGENCE, null ),
                         !root_to_tip && ( BranchLengthLayout.divergenceSource( whole ) == BranchLengthLayout.DIVERGENCE_SOURCE.CLOCK_RATE ),
                         calendar,
                         from_names,
                         ambiguous_names,
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
