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

import org.forester.archaeopteryx.tools.TipDateExtractor;
import org.forester.archaeopteryx.tools.TipDateExtractor.DateMatch;
import org.forester.archaeopteryx.tools.TipDateExtractor.DayMonthOrder;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/**
 * Turns the node HEIGHTS of a tip-dated time tree into calendar DATES when the tip labels say which calendar date
 * height 0 is -- so a BEAST / TreeAnnotator MCC tree opens on the Calendar axis, like a Nextstrain tree does.
 * <p>
 * A BEAST tree states every node's age as a {@code height}: time before the youngest tip, with no unit. That alone
 * cannot give a time axis ({@link AptxUtil#deriveTimeAxisType} will not guess a unit from the size of the numbers).
 * But a tip-dated analysis almost always names its tips with their sampling dates ({@code A_duck_Guangdong_12_2000},
 * {@code NewYork_705_1994.1}, {@code EBOV|KR817226|2014-06-10}), and those give the missing evidence: when the
 * heights are years, every tip's label date plus its height is the same calendar date -- the date of height 0.
 * Heights in months or days, or a label that is a strain number and not a date, break that agreement.
 * <p>
 * The rule, all of it required:
 * <ol>
 * <li>no node's date carries a unit (a unit says what the numbers are; this never overrides one);</li>
 * <li>the tree is a time tree ({@link AptxUtil#isTimeTree}: its internal nodes are dated);</li>
 * <li>a strict majority of the tips carry both a height and a date in the label ({@link TipDateExtractor});</li>
 * <li>each such tip allows height 0 to lie anywhere in its label's calendar range plus its height, widened by
 * {@link #TOLERANCE_YEARS} on each side; the calendar date allowed by the most tips must be allowed by at least
 * {@link #MIN_AGREEING_NUM}/{@link #MIN_AGREEING_DEN} of them, and no separate stretch of dates may be allowed by as
 * many;</li>
 * <li>the agreeing tips must have been sampled at different times: two of their label ranges must not overlap. Tips
 * all from one year agree with heights in any unit, and so prove nothing.</li>
 * </ol>
 * Height 0 is then the median of label date plus height over the most precisely dated agreeing tips, kept inside the
 * dates all of them allow, and each node's date is that minus its height (the older HPD bound becomes the earlier
 * date), in years, rounded to 5 decimals. The heights are not kept: the sentence appended to the tree description
 * names the date of height 0, from which each one follows. (Keeping them as a {@code beast:height} property was
 * measured and rejected: the PearTree Ebola tree, which has no other Color-by field, would open coloured by it.)
 * <p>
 * Measured 2026-09-17 on every tree file in the forester, Archaeopteryx.js and Downloads corpora: the four real BEAST
 * trees (influenza.tree, HA_discrete_MCC, HA_continuous_MCC, the PearTree Ebola example) agree on every tip, and no
 * other tree has unit-less dates.
 */
final class HeightDateConverter {

    /** How far a tip's label date plus its height may miss, in years: covers the one-day differences between the
     *  decimal-year conventions of different programs with room to spare, and is far below the gaps heights in months
     *  or days open up. */
    static final double TOLERANCE_YEARS  = 0.01;
    /** At least this share (19/20) of the tips carrying a height and a label date must agree on height 0. */
    static final int    MIN_AGREEING_NUM = 19;
    static final int    MIN_AGREEING_DEN = 20;
    static final String YEAR_UNIT        = "year";

    /** The calendar date (decimal year) of height 0, how many of the compared tips agree on it, and how many tips were
     *  compared (those with a height and a date in the label). */
    record Anchor(BigDecimal present, int agreeing, int compared) {
    }

    /** Converts every tree of a load that qualifies, appending the provenance sentence to each converted tree's
     *  description. Returns the number of trees converted. Called by every load path, next to
     *  {@link AptxUtil#applyInternalLabelPolicy}. */
    static int convertHeightsToDates( final Phylogeny[] phys ) {
        if ( phys == null ) {
            return 0;
        }
        int converted = 0;
        for( final Phylogeny phy : phys ) {
            final Anchor anchor = inferAnchor( phy );
            if ( anchor == null ) {
                continue;
            }
            convert( phy, anchor.present() );
            final String prov = provenanceSentence( anchor, phy.getName(), phy.getNumberOfExternalNodes() );
            final String existing = phy.getDescription();
            phy.setDescription( ForesterUtil.isEmpty( existing ) ? prov : ( existing + " " + prov ) );
            converted++;
        }
        return converted;
    }

    /** The calendar date of height 0 by the rule in the class comment, or null when the tree does not qualify. Pure. */
    static Anchor inferAnchor( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return null;
        }
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final Date d = it.next().getNodeData().getDate();
            if ( ( d != null ) && !ForesterUtil.isEmpty( d.getUnit() ) ) {
                return null;
            }
        }
        if ( !AptxUtil.isTimeTree( phy ) ) {
            return null;
        }
        int tips = 0;
        final List<double[]> ranges = new ArrayList<>(); // per compared tip: label start, end, height, label date
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode tip = it.next();
            tips++;
            // not isHasDate(): it reads a value of 0 -- the youngest tip's height -- as no date
            final Date d = tip.getNodeData().getDate();
            if ( ( d == null ) || ( d.getValue() == null ) ) {
                continue;
            }
            final DateMatch m = TipDateExtractor.parse( tip.getName(), DayMonthOrder.DAY_FIRST );
            if ( m != null ) {
                ranges.add( new double[] { m.rangeStart(), m.rangeEnd(), d.getValue().doubleValue(),
                        m.decimalYear() } );
            }
        }
        final int compared = ranges.size();
        if ( ( compared * 2 ) <= tips ) {
            return null;
        }
        final double[] best = mostAllowedDates( ranges );
        if ( best == null ) {
            return null; // two separate stretches of dates are equally well supported
        }
        final double mid = ( best[ 0 ] + best[ 1 ] ) / 2;
        final List<double[]> agreeing_tips = new ArrayList<>();
        double latest_start = -Double.MAX_VALUE;
        double earliest_end = Double.MAX_VALUE;
        for( final double[] r : ranges ) {
            if ( allows( r, mid ) ) {
                agreeing_tips.add( r );
                latest_start = Math.max( latest_start, r[ 0 ] );
                earliest_end = Math.min( earliest_end, r[ 1 ] );
            }
        }
        final int agreeing = agreeing_tips.size();
        if ( ( agreeing * MIN_AGREEING_DEN ) < ( compared * MIN_AGREEING_NUM ) ) {
            return null;
        }
        if ( latest_start <= earliest_end ) {
            return null; // every agreeing label range overlaps every other: no two samples are known to differ in time
        }
        final double present = Math.max( best[ 0 ], Math.min( best[ 1 ], precisestMedian( agreeing_tips ) ) );
        return new Anchor( rounded( BigDecimal.valueOf( present ) ), agreeing, compared );
    }

    /** Where the most precisely dated tips put height 0: the median of label date plus height over the tips whose
     *  label range is at most twice as wide as the narrowest. The middle of the stretch all agreeing tips allow would
     *  be pulled about by the coarse labels -- influenza.tree mixes {@code 1993.11} with {@code 1997} (meaning
     *  1997.00), and that middle, 2005.2525, showed the tip labelled 1994.1 as 1994.1025. */
    private static double precisestMedian( final List<double[]> tips ) {
        double narrowest = Double.MAX_VALUE;
        for( final double[] r : tips ) {
            narrowest = Math.min( narrowest, r[ 1 ] - r[ 0 ] );
        }
        final List<Double> values = new ArrayList<>();
        for( final double[] r : tips ) {
            if ( ( r[ 1 ] - r[ 0 ] ) <= ( 2 * narrowest ) ) {
                values.add( r[ 3 ] + r[ 2 ] );
            }
        }
        values.sort( null );
        final int n = values.size();
        return ( ( n % 2 ) == 1 ) ? values.get( n / 2 )
                : ( ( values.get( ( n / 2 ) - 1 ) + values.get( n / 2 ) ) / 2 );
    }

    /** Rewrites every node date as {@code present - height}, the interval bounds swapped (the older bound is the
     *  earlier date), in years. A date without a point value keeps an empty unit, as the parser writes it. */
    static void convert( final Phylogeny phy, final BigDecimal present ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final Date d = n.getNodeData().getDate();
            if ( d == null ) {
                continue;
            }
            final BigDecimal value = ( d.getValue() == null ) ? null : rounded( present.subtract( d.getValue() ) );
            final BigDecimal min = ( d.getMax() == null ) ? null : rounded( present.subtract( d.getMax() ) );
            final BigDecimal max = ( d.getMin() == null ) ? null : rounded( present.subtract( d.getMin() ) );
            n.getNodeData().setDate( new Date( d.getDesc(), value, min, max, ( value != null ) ? YEAR_UNIT : "" ) );
        }
    }

    /** e.g. <i>Converted the node heights of tree named "TREE1" with 190 tips to calendar dates: the sampling dates in
     *  190 of 190 tip labels put height 0 at 2005.5, so each date is 2005.5 minus the height.</i> */
    static String provenanceSentence( final Anchor anchor, final String tree_name, final int num_ext_nodes ) {
        final String present = anchor.present().toPlainString();
        return "Converted the node heights of " + TreePanelUtil.provenanceTreePhrase( tree_name, num_ext_nodes )
                + " to calendar dates: the sampling dates in " + anchor.agreeing() + " of " + anchor.compared()
                + " tip labels put height 0 at " + present + ", so each date is " + present + " minus the height.";
    }

    /** Whether a compared tip ({label start, label end, height, label date}) allows height 0 at calendar date
     *  {@code x}. */
    private static boolean allows( final double[] r, final double x ) {
        return ( x >= ( ( r[ 0 ] + r[ 2 ] ) - TOLERANCE_YEARS ) ) && ( x <= ( r[ 1 ] + r[ 2 ] + TOLERANCE_YEARS ) );
    }

    /** The stretch of calendar dates allowed by the most tips, as {start, end}; null when a separate stretch is allowed
     *  by as many. A sweep over the tips' allowed intervals, closed at both ends. */
    private static double[] mostAllowedDates( final List<double[]> ranges ) {
        final List<double[]> events = new ArrayList<>(); // {x, +1 start / -1 end}
        for( final double[] r : ranges ) {
            events.add( new double[] { ( r[ 0 ] + r[ 2 ] ) - TOLERANCE_YEARS, 1 } );
            events.add( new double[] { r[ 1 ] + r[ 2 ] + TOLERANCE_YEARS, -1 } );
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

    private static BigDecimal rounded( final BigDecimal x ) {
        return new BigDecimal( x.setScale( 5, RoundingMode.HALF_UP ).stripTrailingZeros().toPlainString() );
    }

    private HeightDateConverter() {
    }
}
