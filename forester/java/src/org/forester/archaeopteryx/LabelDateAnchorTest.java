// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.forester.archaeopteryx.LabelDateAnchor.Anchor;
import org.forester.archaeopteryx.LabelDateAnchor.TipOffer;

/**
 * {@link LabelDateAnchor} at the offers level -- the shared rule behind {@link HeightDateConverter} and
 * {@link BranchLengthDateConverter}, tested without building a tree, so the boundaries are stated as numbers: the
 * strict-majority gate, the 19-in-20 agreement floor, a rival stretch, tips that were all sampled at the same time,
 * and where the anchor lands (the precise tips' median, kept inside what the agreeing tips actually state).
 */
public final class LabelDateAnchorTest {

    public static void main( final String[] args ) {
        System.out.println( "LabelDateAnchor: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        try {
            return gates() && agreementFloor() && rivalStretch() && anchorPlacement();
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    /** No offers, too few offers for the tree's tip count, and tips that all state the same time. */
    private static boolean gates() {
        if ( LabelDateAnchor.infer( new ArrayList<>(), 4 ) != null ) {
            return fail( "no offers, no anchor" );
        }
        // two offers from a four-tip tree is not a strict majority, even though they agree perfectly
        final List<TipOffer> two = Arrays.asList( year( 2000, 5 ), year( 2005, 0 ) );
        if ( LabelDateAnchor.infer( two, 4 ) != null ) {
            return fail( "half the tips is not a strict majority" );
        }
        if ( LabelDateAnchor.infer( two, 3 ) == null ) {
            return fail( "two of three tips is" );
        }
        // every label states the same year: they would agree with the offsets in any unit, so they prove nothing
        final List<TipOffer> one_year = Arrays.asList( year( 2020, 0 ), year( 2020, 0 ), year( 2020, 0 ) );
        if ( LabelDateAnchor.infer( one_year, 3 ) != null ) {
            return fail( "tips sampled at one time state nothing about the unit" );
        }
        return true;
    }

    /** 19 of 20 agreeing converts, 18 does not -- and the count is of the OFFERS, not of the tree's tips. */
    private static boolean agreementFloor() {
        for( final int wrong : new int[] { 1, 2 } ) {
            final List<TipOffer> offers = new ArrayList<>();
            for( int i = 0; i < 20; ++i ) {
                final int year = 1996 + ( i % 10 );
                offers.add( year( year, ( 2005 - year ) + ( ( i < wrong ) ? 3 : 0 ) ) );
            }
            final Anchor a = LabelDateAnchor.infer( offers, 20 );
            if ( ( wrong == 1 ) && ( ( a == null ) || ( a.agreeing() != 19 ) || ( a.compared() != 20 ) ) ) {
                return fail( "19 of 20 agreeing is enough, got " + a );
            }
            if ( ( wrong == 2 ) && ( a != null ) ) {
                return fail( "18 of 20 is not" );
            }
        }
        return true;
    }

    /** Two separate stretches of dates, each allowed by as many tips: which one is the anchor is undecided. */
    private static boolean rivalStretch() {
        final List<TipOffer> offers = new ArrayList<>();
        for( int i = 0; i < 18; ++i ) {
            final int y = 1996 + ( i % 10 );
            offers.add( year( y, 2005 - y ) ); // all allow every date in [2005, 2006]
        }
        offers.add( new TipOffer( 2005.0384, 2005.0411, 0, 2005.0397 ) ); // 15 Jan 2005
        offers.add( new TipOffer( 2005.8849, 2005.8877, 0, 2005.8863 ) ); // 20 Nov 2005
        if ( LabelDateAnchor.infer( offers, 20 ) != null ) {
            return fail( "19 tips allow each of two separate days and 18 the time between: undecided" );
        }
        // CONTROL: with only one of the two rivals the same offers DO anchor, so it is the pair that refuses them
        final List<TipOffer> one_rival = offers.subList( 0, 19 );
        if ( LabelDateAnchor.infer( one_rival, 19 ) == null ) {
            return fail( "control: one rival day anchors, or the refusal above is not the rival stretch" );
        }
        return true;
    }

    /** The precise tips place the anchor; the tolerance decides agreement but must not place it. The offers are the
     *  ones the tree-level tests build: labels like {@code 1993.11} (precise, +/- half its last digit) and
     *  {@code 1997} (a whole year), each shifted by its node's height. */
    private static boolean anchorPlacement() {
        // two precise labels nearly agree (2005.250 and 2005.252, median 2005.251); the year-labelled tip states
        // nothing before 2005.253, and all three overlap from there
        final List<TipOffer> overlapping = Arrays.asList( decimals( 1993.11, 12.14 ), decimals( 1993.13, 12.122 ),
                                                          year( 1997, 8.253 ) );
        final Anchor clamped = LabelDateAnchor.infer( overlapping, 3 );
        if ( ( clamped == null ) || ( clamped.value().compareTo( new BigDecimal( "2005.253" ) ) != 0 ) ) {
            return fail( "the anchor is kept inside the dates the agreeing tips state (2005.253), got " + clamped );
        }
        // now no two of them overlap: they agree only within the tolerance, so the precise median stands rather than
        // being pushed to an edge that no tip states
        final List<TipOffer> tolerance_only = Arrays.asList( decimals( 1993.11, 12.14 ), decimals( 2005.25, 0.008 ),
                                                             year( 1997, 8.265 ), year( 1999, 6.0 ) );
        final Anchor median = LabelDateAnchor.infer( tolerance_only, 4 );
        if ( ( median == null ) || ( median.value().compareTo( new BigDecimal( "2005.254" ) ) != 0 ) ) {
            return fail( "with no date in common the precise median stands (2005.254), got " + median );
        }
        return true;
    }

    /** A label written to two decimals, which states its last digit either way (1993.11 is 1993.105 to 1993.115). */
    private static TipOffer decimals( final double label, final double offset ) {
        return new TipOffer( label - 0.005, label + 0.005, offset, label );
    }

    /** A whole-year label at {@code year}, offering {@code year + offset} as the anchor. */
    private static TipOffer year( final int year, final double offset ) {
        return new TipOffer( year, year + 1, offset, year + 0.5 );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [LabelDateAnchorTest] " + msg );
        return false;
    }

    private LabelDateAnchorTest() {
    }
}
