// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import java.math.BigDecimal;
import java.time.LocalDate;
import java.time.Year;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * {@link HeightDateConverter}: BEAST heights become calendar dates only when the tip labels' sampling dates agree with
 * them. Covers the conversion itself (values, swapped interval bounds, rounding away float noise, descriptions kept,
 * the unit, the provenance sentence), where height 0 is put (the precisest tips, clamped), and every refusal: a unit,
 * no dated internal nodes, heights in months or substitutions, contemporaneous tips, too few labelled tips, too little
 * agreement, two equally supported anchors, the tolerance edge.
 */
public final class HeightDateConverterTest {

    /** HA-like: bare-year labels, integer heights, height 0 on c (with a date="" description), float noise on a. */
    private static final String BARE_YEARS = "((a_1996[&height=9.000000000000002,"
            + "height_95%_HPD={9.0,9.000000000000004}]:1,"
            + "b_2000[&height=5]:5)[&height=10,height_95%_HPD={9.5,10.8}]:1,"
            + "(c_2005[&height=0.0,date=\"2005-06-01\"]:4,d_2003[&height=2]:2)[&height=4]:7)[&height=11];";

    public static void main( final String[] args ) {
        System.out.println( "HeightDateConverter: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        try {
            return convertsBareYears() && anchorsOnPrecisestTips() && clampsAnchorIntoAgreement()
                    && dayDatesInBeastConvention() && refusesAUnit() && refusesUndatedInternalNodes()
                    && refusesOtherHeightUnits() && refusesContemporaneousTips() && needsMostTipsLabelled()
                    && needsNineteenOfTwenty() && refusesTwoEquallySupportedAnchors() && toleranceEdge()
                    && heightZeroCounts() && idempotentAndNullSafe() && intervalWithoutValue();
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    private static boolean convertsBareYears() throws Exception {
        final Phylogeny p = nh( BARE_YEARS );
        p.setName( "TREE1" );
        if ( HeightDateConverter.convertHeightsToDates( new Phylogeny[] { p } ) != 1 ) {
            return fail( "a bare-year BEAST tree whose labels agree with its heights must convert" );
        }
        final String sentence = "Converted the node heights of tree named \"TREE1\" with 4 tips to calendar dates: "
                + "the sampling dates in 4 of 4 tip labels put height 0 at 2005.5, so each date is 2005.5 minus the "
                + "height.";
        if ( !sentence.equals( p.getDescription() ) ) {
            return fail( "provenance sentence, got: " + p.getDescription() );
        }
        if ( !dateIs( p.getRoot(), "1994.5", null, null ) ) {
            return false;
        }
        final PhylogenyNode ab = p.getNode( "a_1996" ).getParent();
        if ( !dateIs( ab, "1995.5", "1994.7", "1996" ) ) {
            return fail( "the OLDER HPD bound (height 10.8) must become the EARLIER date" );
        }
        // 2005.5 - 9.000000000000002 rounds to 1996.5; the {9.0, 9.000000000000004} noise interval becomes {d,d}
        if ( !dateIs( p.getNode( "a_1996" ), "1996.5", "1996.5", "1996.5" ) ) {
            return fail( "float noise must round away" );
        }
        final Date c = p.getNode( "c_2005" ).getNodeData().getDate();
        if ( !dateIs( p.getNode( "c_2005" ), "2005.5", null, null ) || !"2005-06-01".equals( c.getDesc() ) ) {
            return fail( "a height of 0 converts, and the date description is kept" );
        }
        if ( AptxUtil.deriveTimeAxisType( p ) != Options.TIME_AXIS_TYPE.CALENDAR ) {
            return fail( "a converted tree must derive the Calendar axis" );
        }
        if ( !AptxUtil.isTimeTree( p ) ) {
            return fail( "a converted tree stays a time tree" );
        }
        if ( AptxUtil.isHasFossilRanges( p ) || AptxUtil.isHasSampledTipWithDateUncertainty( p ) ) {
            return fail( "a rounded-away noise interval is neither a fossil range nor a sampling uncertainty" );
        }
        return true;
    }

    /** influenza.tree's mix: 1993.11 and 2005.25 (to 0.01), 1994.1 (to 0.1), 1997 meaning 1997.00. The middle of the
     *  stretch all four allow is 2005.2525; the two precisest labels say 2005.25. */
    private static boolean anchorsOnPrecisestTips() throws Exception {
        final Phylogeny p = ladder( new String[] { "NY_1_1994.1", "NY_2_1993.11", "NY_3_1997", "NY_4_2005.25" },
                                    new double[] { 11.15, 12.14, 8.25, 0 } );
        final HeightDateConverter.Anchor a = HeightDateConverter.inferAnchor( p );
        if ( ( a == null ) || ( a.present().compareTo( new BigDecimal( "2005.25" ) ) != 0 ) || ( a.agreeing() != 4 ) ) {
            return fail( "height 0 must be where the precisest labels put it (2005.25), got " + a );
        }
        HeightDateConverter.convertHeightsToDates( new Phylogeny[] { p } );
        if ( !dateIs( p.getNode( "NY_1_1994.1" ), "1994.1", null, null ) ) {
            return fail( "a tip labelled 1994.1 must read 1994.1" );
        }
        return true;
    }

    /** The precisest tips say 2005.254 (median of 2005.25 and 2005.258), but a year-labelled tip allows nothing before
     *  2005.255: height 0 is kept where every agreeing tip allows it. */
    private static boolean clampsAnchorIntoAgreement() throws Exception {
        final Phylogeny p = ladder( new String[] { "p_1993.11", "q_2005.25", "r_1997", "s_1999" },
                                    new double[] { 12.14, 0.008, 8.265, 6.0 } );
        final HeightDateConverter.Anchor a = HeightDateConverter.inferAnchor( p );
        if ( ( a == null ) || ( a.present().compareTo( new BigDecimal( "2005.255" ) ) != 0 ) ) {
            return fail( "height 0 must be clamped into the stretch every agreeing tip allows (2005.255), got " + a );
        }
        return true;
    }

    /** Day-precision labels with heights computed as BEAST does ((day of year - 1) / days in year): every tip's date
     *  lands inside its own labelled day. */
    private static boolean dayDatesInBeastConvention() throws Exception {
        final String[] days = { "2014-03-17", "2014-06-10", "2014-11-02", "2015-01-29", "2015-05-12", "2015-10-24" };
        final String[] names = new String[ days.length ];
        final double[] heights = new double[ days.length ];
        final double youngest = beastDecimal( days[ days.length - 1 ] );
        for( int i = 0; i < days.length; ++i ) {
            names[ i ] = "EBOV|KR8172" + i + "|" + days[ i ];
            heights[ i ] = youngest - beastDecimal( days[ i ] );
        }
        final Phylogeny p = ladder( names, heights );
        if ( HeightDateConverter.convertHeightsToDates( new Phylogeny[] { p } ) != 1 ) {
            return fail( "day-dated labels with BEAST-convention heights must convert" );
        }
        for( int i = 0; i < days.length; ++i ) {
            final double v = p.getNode( names[ i ] ).getNodeData().getDate().getValue().doubleValue();
            final LocalDate d = LocalDate.parse( days[ i ] );
            final double len = Year.of( d.getYear() ).length();
            if ( ( v < ( d.getYear() + ( ( d.getDayOfYear() - 1 ) / len ) - 1e-5 ) )
                    || ( v > ( d.getYear() + ( d.getDayOfYear() / len ) + 1e-5 ) ) ) {
                return fail( names[ i ] + " converted to " + v + ", outside its labelled day" );
            }
        }
        return true;
    }

    /** A unit ANYWHERE stops it -- on the root, on an inner node, or on a single tip. One node is enough: a file that
     *  says what its numbers are is never overridden, and a tree cannot be half in years. */
    private static boolean refusesAUnit() throws Exception {
        for( final String unit : new String[] { "year", "mya", "days" } ) {
            for( final String where : new String[] { "root", "inner", "tip" } ) {
                final Phylogeny p = nh( BARE_YEARS );
                p.setDescription( "kept" );
                nodeAt( p, where ).getNodeData().getDate().setUnit( unit );
                if ( ( HeightDateConverter.convertHeightsToDates( new Phylogeny[] { p } ) != 0 )
                        || !"kept".equals( p.getDescription() )
                        || !"11".equals( p.getRoot().getNodeData().getDate().getValue().toPlainString() ) ) {
                    return fail( "a date unit (" + unit + " on the " + where + ") says what the numbers are: never converted" );
                }
            }
        }
        return true;
    }

    private static PhylogenyNode nodeAt( final Phylogeny p, final String where ) {
        switch ( where ) {
            case "root":
                return p.getRoot();
            case "inner":
                return p.getNode( "a_1996" ).getParent();
            default:
                return p.getNode( "d_2003" );
        }
    }

    private static boolean refusesUndatedInternalNodes() throws Exception {
        final Phylogeny p = nh( "((a_1996[&height=9]:1,b_2000[&height=5]:5):1,"
                + "(c_2005[&height=0]:4,d_2003[&height=2]:2):7);" );
        if ( HeightDateConverter.inferAnchor( p ) != null ) {
            return fail( "dated tips alone are not a time tree" );
        }
        return true;
    }

    /** The same tree with its heights in months, and in substitutions per site: the labels no longer agree. */
    private static boolean refusesOtherHeightUnits() throws Exception {
        final String[] names = { "a_1996", "b_2000", "c_2005", "d_2003" };
        final double[] years = { 9, 5, 0, 2 };
        for( final double factor : new double[] { 12, 0.001 } ) {
            final double[] h = new double[ years.length ];
            for( int i = 0; i < h.length; ++i ) {
                h[ i ] = years[ i ] * factor;
            }
            if ( HeightDateConverter.inferAnchor( ladder( names, h ) ) != null ) {
                return fail( "heights " + factor + " x years must not be read as years" );
            }
        }
        if ( HeightDateConverter.inferAnchor( ladder( names, years ) ) == null ) {
            return fail( "control: the same tree in years converts" );
        }
        return true;
    }

    /** Four tips all from 2020 at height 0 agree with any unit: no evidence. */
    private static boolean refusesContemporaneousTips() throws Exception {
        final Phylogeny p = ladder( new String[] { "a_2020", "b_2020", "c_2020", "d_2020" },
                                    new double[] { 0, 0, 0, 0 } );
        if ( HeightDateConverter.inferAnchor( p ) != null ) {
            return fail( "tips from one year prove nothing about the unit" );
        }
        return true;
    }

    private static boolean needsMostTipsLabelled() throws Exception {
        if ( HeightDateConverter.inferAnchor( ladder( new String[] { "a_1996", "b_2000", "c_x", "d_y" },
                                                      new double[] { 9, 5, 0, 2 } ) ) != null ) {
            return fail( "half the tips labelled is not a strict majority" );
        }
        if ( HeightDateConverter.inferAnchor( ladder( new String[] { "a_1996", "b_2000", "c_2005", "d_y" },
                                                      new double[] { 9, 5, 0, 2 } ) ) == null ) {
            return fail( "three of four tips labelled is a strict majority" );
        }
        return true;
    }

    /** 20 tips: 19 agreeing converts, 18 does not. */
    private static boolean needsNineteenOfTwenty() throws Exception {
        for( final int wrong : new int[] { 1, 2 } ) {
            final String[] names = new String[ 20 ];
            final double[] h = new double[ 20 ];
            for( int i = 0; i < 20; ++i ) {
                final int year = 1996 + ( i % 10 );
                names[ i ] = "t" + i + "_" + year;
                h[ i ] = ( 2005 - year ) + ( ( i < wrong ) ? 3 : 0 ); // a wrong tip's label is 3 years off
            }
            final HeightDateConverter.Anchor a = HeightDateConverter.inferAnchor( ladder( names, h ) );
            if ( ( wrong == 1 ) && ( ( a == null ) || ( a.agreeing() != 19 ) || ( a.compared() != 20 ) ) ) {
                return fail( "19 of 20 tips agreeing must convert, got " + a );
            }
            if ( ( wrong == 2 ) && ( a != null ) ) {
                return fail( "18 of 20 tips agreeing must not convert" );
            }
        }
        return true;
    }

    /** 18 year-labelled tips allow all of 2005; two day-labelled tips at height 0 disagree with each other (15 Jan vs
     *  20 Nov). 19 tips allow each of the two days, 18 the time between: which day is height 0 is not decided. */
    private static boolean refusesTwoEquallySupportedAnchors() throws Exception {
        final String[] names = new String[ 20 ];
        final double[] h = new double[ 20 ];
        for( int i = 0; i < 18; ++i ) {
            final int year = 1996 + ( i % 10 );
            names[ i ] = "t" + i + "_" + year;
            h[ i ] = 2005 - year;
        }
        names[ 18 ] = "x_2005-01-15";
        names[ 19 ] = "y_2005-11-20";
        if ( HeightDateConverter.inferAnchor( ladder( names, h ) ) != null ) {
            return fail( "two separate, equally supported dates for height 0 must not convert" );
        }
        return true;
    }

    /** d_2003's allowed dates start 0.019 years after the others end, less than twice the tolerance: agrees. At
     *  0.021 it does not, and 3 of 4 is too few. */
    private static boolean toleranceEdge() throws Exception {
        final String[] names = { "a_1996", "b_2000", "c_2005", "d_2003" };
        if ( HeightDateConverter.inferAnchor( ladder( names, new double[] { 9, 5, 0, 3.019 } ) ) == null ) {
            return fail( "a miss of 0.019 years is within the tolerance on both sides" );
        }
        if ( HeightDateConverter.inferAnchor( ladder( names, new double[] { 9, 5, 0, 3.021 } ) ) != null ) {
            return fail( "a miss of 0.021 years is beyond it" );
        }
        return true;
    }

    /** Two youngest tips at height exactly 0 (which NodeData.isHasDate reads as no date) must count. */
    private static boolean heightZeroCounts() throws Exception {
        final HeightDateConverter.Anchor a = HeightDateConverter.inferAnchor( ladder( new String[] { "a_1996", "b_2005",
                "c_2005", "d_2000" }, new double[] { 9, 0, 0, 5 } ) );
        if ( ( a == null ) || ( a.compared() != 4 ) ) {
            return fail( "tips at height 0 are compared, got " + a );
        }
        return true;
    }

    private static boolean idempotentAndNullSafe() throws Exception {
        final Phylogeny p = nh( BARE_YEARS );
        p.setDescription( "An H5N1 tree." );
        final Phylogeny[] phys = { null, new Phylogeny(), p };
        if ( ( HeightDateConverter.convertHeightsToDates( phys ) != 1 )
                || ( HeightDateConverter.convertHeightsToDates( phys ) != 0 )
                || ( HeightDateConverter.convertHeightsToDates( null ) != 0 ) ) {
            return fail( "converts once, then the unit \"year\" stops it; null and empty trees are skipped" );
        }
        if ( !p.getDescription().startsWith( "An H5N1 tree. Converted the node heights of a tree with 4 tips" )
                || ( p.getDescription().indexOf( "Converted" ) != p.getDescription().lastIndexOf( "Converted" ) ) ) {
            return fail( "the sentence is appended once, after the existing description: " + p.getDescription() );
        }
        return true;
    }

    /** A node with only an HPD interval (no point height) converts its bounds and keeps an empty unit, as the parser
     *  writes such a date; an undated internal node stays undated. */
    private static boolean intervalWithoutValue() throws Exception {
        final Phylogeny p = nh( "(((a_1996[&height=9]:1,b_2000[&height=5]:5)[&height_95%_HPD={9.5,10.8}]:1,"
                + "(c_2005[&height=0]:4,d_2003[&height=2]:2)[&height=4]:7)[&height=11]:1,"
                + "(e_2004[&height=1]:3,f_2001[&height=4]:0))[&height=12];" );
        if ( HeightDateConverter.convertHeightsToDates( new Phylogeny[] { p } ) != 1 ) {
            return fail( "control: this tree converts" );
        }
        final Date d = p.getNode( "a_1996" ).getParent().getNodeData().getDate();
        if ( ( d.getValue() != null ) || !"".equals( d.getUnit() )
                || ( d.getMin().compareTo( new BigDecimal( "1994.7" ) ) != 0 )
                || ( d.getMax().compareTo( new BigDecimal( "1996" ) ) != 0 ) ) {
            return fail( "an interval-only date: bounds converted, unit empty, got " + d.getMin() + " " + d.getMax()
                    + " '" + d.getUnit() + "'" );
        }
        if ( p.getNode( "e_2004" ).getParent().getNodeData().getDate() != null ) {
            return fail( "an undated node stays undated" );
        }
        return oneBoundOnly();
    }

    /** One bound without the other -- which the node-data editor can save, since its minimum and maximum are separate
     *  fields. The bounds swap sides: a lower AGE bound is the LATER calendar date, so a min-only age becomes a
     *  max-only date. */
    private static boolean oneBoundOnly() throws Exception {
        final Phylogeny p = nh( BARE_YEARS );
        final PhylogenyNode inner = p.getNode( "a_1996" ).getParent();
        // keep each node's height: stripping the values would leave too few dated internal nodes and the tree would
        // stop being a time tree -- the control below would then fail for the wrong reason
        inner.getNodeData().setDate( new Date( "", new BigDecimal( "10" ), new BigDecimal( "9.5" ), null, "" ) );
        final PhylogenyNode other = p.getNode( "c_2005" ).getParent();
        other.getNodeData().setDate( new Date( "", new BigDecimal( "4" ), null, new BigDecimal( "4.5" ), "" ) );
        if ( HeightDateConverter.convertHeightsToDates( new Phylogeny[] { p } ) != 1 ) {
            return fail( "control: this tree converts" );
        }
        final Date a = inner.getNodeData().getDate();
        if ( ( a.getMin() != null ) || ( a.getMax() == null ) || ( a.getMax().compareTo( new BigDecimal( "1996" ) ) != 0 ) ) {
            return fail( "a min-only age becomes a max-only date (2005.5 - 9.5 = 1996), got " + a.getMin() + ".."
                    + a.getMax() );
        }
        final Date b = other.getNodeData().getDate();
        if ( ( b.getMax() != null ) || ( b.getMin() == null ) || ( b.getMin().compareTo( new BigDecimal( "2001" ) ) != 0 ) ) {
            return fail( "a max-only age becomes a min-only date (2005.5 - 4.5 = 2001), got " + b.getMin() + ".."
                    + b.getMax() );
        }
        return true;
    }

    // ---- helpers ----------------------------------------------------------------------------------------------

    private static double beastDecimal( final String iso ) {
        final LocalDate d = LocalDate.parse( iso );
        return d.getYear() + ( ( d.getDayOfYear() - 1.0 ) / Year.of( d.getYear() ).length() );
    }

    /** A ladder tree: tip 0 and tip 1 join first, then each further tip; every internal node dated one year older than
     *  its oldest child. */
    private static Phylogeny ladder( final String[] names, final double[] heights ) throws Exception {
        String s = names[ 0 ] + "[&height=" + plain( heights[ 0 ] ) + "]";
        double oldest = heights[ 0 ];
        for( int i = 1; i < names.length; ++i ) {
            oldest = Math.max( oldest, heights[ i ] ) + 1;
            s = "(" + s + "," + names[ i ] + "[&height=" + plain( heights[ i ] ) + "])[&height=" + plain( oldest )
                    + "]";
        }
        return nh( s + ";" );
    }

    private static String plain( final double d ) {
        return BigDecimal.valueOf( d ).toPlainString();
    }

    private static Phylogeny nh( final String s ) throws Exception {
        return ParserBasedPhylogenyFactory.getInstance().create( s, new NHXParser() )[ 0 ];
    }

    /** The node's date equals value/min/max (null = absent) as numbers AND as written (no float noise), unit "year". */
    private static boolean dateIs( final PhylogenyNode n, final String value, final String min, final String max ) {
        final Date d = n.getNodeData().getDate();
        final String got = ( d == null ) ? "no date"
                : ( d.getValue() + " [" + d.getMin() + ", " + d.getMax() + "] '" + d.getUnit() + "'" );
        if ( ( d == null ) || !same( d.getValue(), value ) || !same( d.getMin(), min ) || !same( d.getMax(), max ) ) {
            return fail( "node " + n.getName() + ": expected " + value + " [" + min + ", " + max + "], got " + got );
        }
        if ( ( d.getValue() != null ) && !"year".equals( d.getUnit() ) ) {
            return fail( "node " + n.getName() + ": a converted value is in years, got " + got );
        }
        return true;
    }

    private static boolean same( final BigDecimal got, final String expected ) {
        return ( expected == null ) ? ( got == null ) : ( ( got != null ) && got.toPlainString().equals( expected ) );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [HeightDateConverterTest] " + msg );
        return false;
    }

    private HeightDateConverterTest() {
    }
}
