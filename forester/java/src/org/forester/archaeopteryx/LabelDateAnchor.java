// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;

/**
 * Where a tree's tip labels put its time anchor -- the one calendar date from which every other date on the tree
 * follows. Shared by the two ways a tree can state time without saying what its numbers mean:
 * <ul>
 * <li>{@link HeightDateConverter}: BEAST HEIGHTS, unit-less ages before the youngest tip. Each tip offers
 * {@code label + height}, and the anchor is the date of height 0.</li>
 * <li>{@link BranchLengthDateConverter}: BRANCH LENGTHS with no dates at all. Each tip offers
 * {@code label - distance from the root}, and the anchor is the date of the root.</li>
 * </ul>
 * One rule serves both because the question is the same: do the tip labels' sampling dates, shifted by what the tree
 * says about each tip's position in time, all point at one date? If they do, the tree's numbers are years.
 * <p>
 * A tip offers not a date but the RANGE its label states ({@code 2012} is all of 2012, an ISO day is that day), widened
 * by {@link #TOLERANCE_YEARS} at each end, which covers the day-either-way differences between programs' decimal-year
 * conventions and stays far below the gaps a wrong unit opens up. The anchor is then the date allowed by the most tips
 * -- required to be at least {@link #MIN_AGREEING_NUM}/{@link #MIN_AGREEING_DEN} of them, with no separate stretch of
 * dates allowed by as many -- placed by the median of what the most precisely dated agreeing tips say, and kept inside
 * the dates those tips actually state (the tolerance decides agreement; it does not place the anchor). Two of the agreeing tips must have been sampled at different times: tips from one year
 * agree with any unit and so prove nothing.
 */
final class LabelDateAnchor {

    /** How far a tip may miss the anchor, in years (about four days). */
    static final double TOLERANCE_YEARS  = 0.01;
    /** At least this share (19/20) of the compared tips must agree on the anchor. */
    static final int    MIN_AGREEING_NUM = 19;
    static final int    MIN_AGREEING_DEN = 20;

    /** The anchor date (decimal year), how many of the compared tips agree on it, and how many were compared. */
    record Anchor(BigDecimal value, int agreeing, int compared) {
    }

    /**
     * One compared tip: the calendar range its label states, the OFFSET the tree adds to that range to point at the
     * anchor (a height, or minus a distance from the root), and the label's own decimal year -- the point estimate
     * used when the precise tips place the anchor.
     */
    record TipOffer(double rangeStart, double rangeEnd, double offset, double labelDate) {
    }

    /**
     * The anchor the offers agree on, or null when they do not: fewer than a strict majority of the tree's {@code tips}
     * made an offer, too few of them agree, a second stretch of dates is equally well supported, or the agreeing tips
     * were all sampled at the same time. Pure.
     */
    static Anchor infer( final List<TipOffer> offers, final int tips ) {
        final int compared = offers.size();
        if ( ( compared * 2 ) <= tips ) {
            return null;
        }
        final double[] best = mostAllowedDates( offers );
        if ( best == null ) {
            return null; // two separate stretches of dates are equally well supported
        }
        final double mid = ( best[ 0 ] + best[ 1 ] ) / 2;
        final List<TipOffer> agreeing_tips = new ArrayList<>();
        double latest_start = -Double.MAX_VALUE;
        double earliest_end = Double.MAX_VALUE;
        for( final TipOffer o : offers ) {
            if ( allows( o, mid ) ) {
                agreeing_tips.add( o );
                latest_start = Math.max( latest_start, o.rangeStart() );
                earliest_end = Math.min( earliest_end, o.rangeEnd() );
            }
        }
        final int agreeing = agreeing_tips.size();
        if ( ( agreeing * MIN_AGREEING_DEN ) < ( compared * MIN_AGREEING_NUM ) ) {
            return null;
        }
        if ( latest_start <= earliest_end ) {
            return null; // every agreeing label range overlaps every other: no two samples are known to differ in time
        }
        // The tolerance decides AGREEMENT; it must not place the anchor. Clamp into the dates the agreeing tips
        // actually state -- their un-widened intersection -- and where they agree only thanks to the tolerance (an
        // empty intersection) let the precise median stand. Measured against ground truth: three Nextstrain time
        // trees whose .nexus siblings carry the real dates were out by 4.2, 4.0 and 2.9 days when the anchor was
        // clamped to the widened stretch, and by 0.5, 0.0 and 0.1 days once it was not. The four BEAST trees are
        // unaffected (their medians already lie inside).
        double lo = -Double.MAX_VALUE;
        double hi = Double.MAX_VALUE;
        for( final TipOffer o : agreeing_tips ) {
            lo = Math.max( lo, o.rangeStart() + o.offset() );
            hi = Math.min( hi, o.rangeEnd() + o.offset() );
        }
        final double median = precisestMedian( agreeing_tips );
        final double anchor = ( lo <= hi ) ? Math.max( lo, Math.min( hi, median ) ) : median;
        return new Anchor( rounded( BigDecimal.valueOf( anchor ) ), agreeing, compared );
    }

    /** Where the most precisely dated tips put the anchor: the median of label date plus offset over the tips whose
     *  label range is at most twice as wide as the narrowest. The middle of the stretch all agreeing tips allow would
     *  be pulled about by the coarse labels -- influenza.tree mixes {@code 1993.11} with {@code 1997} (meaning
     *  1997.00), and that middle, 2005.2525, showed the tip labelled 1994.1 as 1994.1025. */
    private static double precisestMedian( final List<TipOffer> tips ) {
        double narrowest = Double.MAX_VALUE;
        for( final TipOffer o : tips ) {
            narrowest = Math.min( narrowest, o.rangeEnd() - o.rangeStart() );
        }
        final List<Double> values = new ArrayList<>();
        for( final TipOffer o : tips ) {
            if ( ( o.rangeEnd() - o.rangeStart() ) <= ( 2 * narrowest ) ) {
                values.add( o.labelDate() + o.offset() );
            }
        }
        values.sort( null );
        final int n = values.size();
        return ( ( n % 2 ) == 1 ) ? values.get( n / 2 )
                : ( ( values.get( ( n / 2 ) - 1 ) + values.get( n / 2 ) ) / 2 );
    }

    /** Whether a tip allows the anchor to sit at calendar date {@code x}. */
    private static boolean allows( final TipOffer o, final double x ) {
        return ( x >= ( ( o.rangeStart() + o.offset() ) - TOLERANCE_YEARS ) )
                && ( x <= ( o.rangeEnd() + o.offset() + TOLERANCE_YEARS ) );
    }

    /** The stretch of calendar dates allowed by the most tips, as {start, end}; null when a separate stretch is allowed
     *  by as many. A sweep over the tips' allowed intervals, closed at both ends. */
    private static double[] mostAllowedDates( final List<TipOffer> offers ) {
        final List<double[]> events = new ArrayList<>(); // {x, +1 start / -1 end}
        for( final TipOffer o : offers ) {
            events.add( new double[] { ( o.rangeStart() + o.offset() ) - TOLERANCE_YEARS, 1 } );
            events.add( new double[] { o.rangeEnd() + o.offset() + TOLERANCE_YEARS, -1 } );
        }
        // starts before ends at the same x: intervals that touch overlap
        events.sort( ( a, b ) -> ( a[ 0 ] != b[ 0 ] ) ? Double.compare( a[ 0 ], b[ 0 ] )
                : Double.compare( b[ 1 ], a[ 1 ] ) );
        int depth = 0;
        int best = 0;
        double start = 0;
        double end = 0;
        boolean open = false;
        boolean tied = false;
        for( final double[] e : events ) {
            if ( e[ 1 ] > 0 ) {
                depth++;
                if ( depth > best ) {
                    best = depth;
                    start = e[ 0 ];
                    open = true;
                    tied = false;
                }
                else if ( ( depth == best ) && !open ) {
                    tied = true;
                }
            }
            else {
                if ( open ) {
                    end = e[ 0 ];
                    open = false;
                }
                depth--;
            }
        }
        return ( ( best == 0 ) || tied ) ? null : new double[] { start, end };
    }

    /** Five decimals (about five minutes of a year), so a converted date reads as the label wrote it rather than
     *  carrying the float noise of the arithmetic that produced it. */
    static BigDecimal rounded( final BigDecimal x ) {
        return new BigDecimal( x.setScale( 5, RoundingMode.HALF_UP ).stripTrailingZeros().toPlainString() );
    }

    private LabelDateAnchor() {
    }
}
