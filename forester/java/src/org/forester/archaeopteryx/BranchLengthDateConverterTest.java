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
 * {@link BranchLengthDateConverter}: a Newick time tree with NO dates is dated from its tip labels, when its branch
 * lengths are years. Covers the dating itself (root date, every node at root + its distance, the unit, the Calendar
 * axis, the provenance sentence) and the refusals that keep a DIVERGENCE tree -- the same file shape, lengths in
 * substitutions -- from being read as time: precise labels that cannot agree, year-only labels one year apart (which
 * touch at a single date), a tree that already carries dates, one without branch lengths, too few dated labels, and
 * tips all sampled at the same time.
 */
public final class BranchLengthDateConverterTest {

    public static void main( final String[] args ) {
        System.out.println( "BranchLengthDateConverter: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        try {
            return datesATimeTree() && refusesADivergenceTree() && refusesTwoAdjacentYears() && refusesADatedTree()
                    && refusesWithoutBranchLengths() && needsMostTipsLabelled() && unmeasuredTipDoesNotVote()
                    && refusesContemporaneousTips()
                    && idempotentAndNullSafe();
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    /** Six day-labelled tips whose distances from the root are the years between the root and each sample. */
    private static boolean datesATimeTree() throws Exception {
        final String[] days = { "2014-03-17", "2014-06-10", "2014-11-02", "2015-01-29", "2015-05-12", "2015-10-24" };
        final double root = decimalYear( "2013-08-01" );
        final Phylogeny phy = timeTree( days, root );
        phy.setName( "timetree" );
        if ( BranchLengthDateConverter.dateTreesFromBranchLengths( new Phylogeny[] { phy } ) != 1 ) {
            return fail( "a Newick time tree whose labels agree with its branch lengths must be dated" );
        }
        for( final String day : days ) {
            final PhylogenyNode tip = phy.getNode( "EBOV|" + day );
            final Date d = tip.getNodeData().getDate();
            final LocalDate ld = LocalDate.parse( day );
            final double len = Year.of( ld.getYear() ).length();
            if ( ( d == null ) || !"year".equals( d.getUnit() )
                    || ( d.getValue().doubleValue() < ( ld.getYear() + ( ( ld.getDayOfYear() - 1 ) / len ) - 1e-5 ) )
                    || ( d.getValue().doubleValue() > ( ld.getYear() + ( ld.getDayOfYear() / len ) + 1e-5 ) ) ) {
                return fail( day + " must be dated inside its own labelled day, got " + d );
            }
        }
        final double dated_root = phy.getRoot().getNodeData().getDate().getValue().doubleValue();
        if ( Math.abs( dated_root - root ) > ( 1.0 / 365 ) ) {
            return fail( "the root must be dated where the tips put it (" + root + "), got " + dated_root );
        }
        if ( AptxUtil.deriveTimeAxisType( phy ) != Options.TIME_AXIS_TYPE.CALENDAR ) {
            return fail( "a dated tree must derive the Calendar axis" );
        }
        if ( !AptxUtil.isTimeTree( phy ) ) {
            return fail( "its internal nodes are dated, so it is a time tree (and must refuse re-rooting)" );
        }
        if ( !phy.getDescription().startsWith( "Dated tree named \"timetree\" with 6 tips from the sampling dates in 6 "
                + "of 6 tip labels: the branch lengths are years, and the root is at " ) ) {
            return fail( "provenance sentence, got: " + phy.getDescription() );
        }
        return true;
    }

    /** THE case this rule exists to refuse: the same tree with its lengths in substitutions per site. Precise labels
     *  spread over years cannot agree on a root date when every tip sits ~0.001 from it. */
    private static boolean refusesADivergenceTree() throws Exception {
        final String[] days = { "2014-03-17", "2014-06-10", "2014-11-02", "2015-01-29", "2015-05-12", "2015-10-24" };
        final Phylogeny div = divergenceTree( days, 0.001 );
        if ( BranchLengthDateConverter.inferRootDate( div ) != null ) {
            return fail( "substitutions per site are not years" );
        }
        return true;
    }

    /** The subtler half: a divergence tree whose labels carry only a YEAR. Two adjacent years' ranges touch at exactly
     *  one date, so the tips can agree -- and the different-sampling-times rule refuses it, because touching ranges
     *  overlap. This is the shape that made a naive root-to-tip check accept 2535 of 2537 tips of a real SARS-CoV-2
     *  divergence tree. */
    private static boolean refusesTwoAdjacentYears() throws Exception {
        final String[] years = { "iso_2020", "iso_2020", "iso_2021", "iso_2021", "iso_2020", "iso_2021" };
        if ( BranchLengthDateConverter.inferRootDate( divergenceTree( years, 0.001 ) ) != null ) {
            return fail( "year-only labels one year apart touch at a single date and prove nothing" );
        }
        // three years apart cannot overlap at all, so it fails on agreement instead
        final String[] three = { "iso_2019", "iso_2019", "iso_2020", "iso_2021", "iso_2021", "iso_2020" };
        if ( BranchLengthDateConverter.inferRootDate( divergenceTree( three, 0.001 ) ) != null ) {
            return fail( "year-only labels three years apart cannot agree on a root date" );
        }
        // CONTROL: the same year-only labels on a tree whose lengths ARE years must date, or the refusals above would
        // pass on a fixture that can never convert whatever its unit
        final String[] names = new String[ three.length ];
        final double[] lengths = new double[ three.length ];
        for( int i = 0; i < three.length; ++i ) {
            names[ i ] = "seq" + i + "_" + three[ i ].substring( 4 );
            lengths[ i ] = ( Integer.parseInt( three[ i ].substring( 4 ) ) + 0.5 ) - 2018.0; // root at 2018.0
        }
        // (to within the 1e-4 trunk this ladder puts between its branching points)
        final LabelDateAnchor.Anchor dated = BranchLengthDateConverter.inferRootDate( ladder( names, lengths ) );
        if ( ( dated == null ) || ( Math.abs( dated.value().doubleValue() - 2018.0 ) > 1e-3 ) ) {
            return fail( "control: the same labels with the lengths in YEARS must date the root 2018, got " + dated );
        }
        return true;
    }

    /** A tree that carries dates is not this case -- whatever they are, and even on one node. */
    private static boolean refusesADatedTree() throws Exception {
        final String[] days = { "2014-03-17", "2014-06-10", "2014-11-02", "2015-01-29", "2015-05-12", "2015-10-24" };
        final Phylogeny phy = timeTree( days, decimalYear( "2013-08-01" ) );
        phy.getFirstExternalNode().getNodeData().setDate( new Date( "", new BigDecimal( "3" ), null, null, "" ) );
        if ( BranchLengthDateConverter.inferRootDate( phy ) != null ) {
            return fail( "one node with a date is enough: the file has its own idea of time" );
        }
        return true;
    }

    private static boolean refusesWithoutBranchLengths() throws Exception {
        final Phylogeny phy = nh( "((a_2014-03-17,b_2014-06-10),(c_2015-01-29,d_2015-10-24));" );
        if ( BranchLengthDateConverter.inferRootDate( phy ) != null ) {
            return fail( "a cladogram says nothing about time" );
        }
        return true;
    }

    private static boolean needsMostTipsLabelled() throws Exception {
        final String[] half = { "2014-03-17", "2014-06-10", "2014-11-02" };
        final Phylogeny phy = timeTree( half, decimalYear( "2013-08-01" ) );
        // rename two of the three tips so only one carries a date
        phy.getNode( "EBOV|2014-03-17" ).setName( "no_date_here" );
        phy.getNode( "EBOV|2014-06-10" ).setName( "nor_here" );
        if ( BranchLengthDateConverter.inferRootDate( phy ) != null ) {
            return fail( "one dated label of three is not a strict majority" );
        }
        return true;
    }

    /** A tip whose path to the root carries no branch length sits at distance 0 and knows nothing about time, so it
     *  must not vote. Here its label is 24 years off: counted, it would drag the agreement below the floor and the
     *  tree -- which the other six tips date perfectly well -- would be refused. */
    private static boolean unmeasuredTipDoesNotVote() throws Exception {
        final String[] days = { "2014-03-17", "2014-06-10", "2014-11-02", "2015-01-29", "2015-05-12", "2015-10-24" };
        final Phylogeny phy = timeTree( days, decimalYear( "2013-08-01" ) );
        final PhylogenyNode stray = new PhylogenyNode(); // no branch length: distance from the root is 0
        stray.setName( "EBOV|1990-01-01" );
        phy.getRoot().addAsChild( stray );
        phy.externalNodesHaveChanged();
        final LabelDateAnchor.Anchor a = BranchLengthDateConverter.inferRootDate( phy );
        if ( ( a == null ) || ( a.compared() != 6 ) ) {
            return fail( "a tip with no branch length must be left out of the vote, got " + a );
        }
        return true;
    }

    private static boolean refusesContemporaneousTips() throws Exception {
        final String[] same = { "2014-03-17", "2014-03-17", "2014-03-17", "2014-03-17" };
        if ( BranchLengthDateConverter.inferRootDate( timeTree( same, decimalYear( "2013-08-01" ) ) ) != null ) {
            return fail( "tips sampled at one moment fit branch lengths in any unit" );
        }
        return true;
    }

    private static boolean idempotentAndNullSafe() throws Exception {
        final String[] days = { "2014-03-17", "2014-06-10", "2014-11-02", "2015-01-29", "2015-05-12", "2015-10-24" };
        final Phylogeny phy = timeTree( days, decimalYear( "2013-08-01" ) );
        phy.setDescription( "An Ebola tree." );
        final Phylogeny[] phys = { null, new Phylogeny(), phy };
        if ( ( BranchLengthDateConverter.dateTreesFromBranchLengths( phys ) != 1 )
                || ( BranchLengthDateConverter.dateTreesFromBranchLengths( phys ) != 0 )
                || ( BranchLengthDateConverter.dateTreesFromBranchLengths( null ) != 0 ) ) {
            return fail( "dates once, then the dates it wrote stop it; null and empty trees are skipped" );
        }
        if ( !phy.getDescription().startsWith( "An Ebola tree. Dated a tree with 6 tips" )
                || ( phy.getDescription().indexOf( "Dated" ) != phy.getDescription().lastIndexOf( "Dated" ) ) ) {
            return fail( "the sentence is appended once: " + phy.getDescription() );
        }
        return true;
    }

    // ---- helpers ----------------------------------------------------------------------------------------------

    /** A date as a decimal year, counting a day from its start (as the programs that write these trees do). */
    private static double decimalYear( final String iso ) {
        final LocalDate d = LocalDate.parse( iso );
        return d.getYear() + ( ( d.getDayOfYear() - 1.0 ) / Year.of( d.getYear() ).length() );
    }

    /** A ladder whose branch lengths are the YEARS between the root and each labelled tip. */
    private static Phylogeny timeTree( final String[] days, final double root ) throws Exception {
        final String[] names = new String[ days.length ];
        final double[] lengths = new double[ days.length ];
        for( int i = 0; i < days.length; ++i ) {
            names[ i ] = "EBOV|" + days[ i ];
            lengths[ i ] = decimalYear( days[ i ] ) - root;
        }
        return ladder( names, lengths );
    }

    /** The same shape with the lengths in substitutions per site: the tips sit {@code scale} x their time from the
     *  root, which is what a divergence tree looks like. Labels are used verbatim. */
    private static Phylogeny divergenceTree( final String[] labels, final double scale ) throws Exception {
        final String[] names = new String[ labels.length ];
        final double[] lengths = new double[ labels.length ];
        for( int i = 0; i < labels.length; ++i ) {
            names[ i ] = labels[ i ].startsWith( "iso_" ) ? ( "seq" + i + "_" + labels[ i ].substring( 4 ) )
                    : ( "EBOV|" + labels[ i ] );
            lengths[ i ] = ( 1 + ( i * 0.37 ) ) * scale;
        }
        return ladder( names, lengths );
    }

    /** A ladder: every tip hangs off the trunk, at its own distance from the root; no dates anywhere. */
    private static Phylogeny ladder( final String[] names, final double[] distances ) throws Exception {
        final StringBuilder sb = new StringBuilder();
        double trunk = 0;
        for( int i = 0; i < names.length; ++i ) {
            final double step = ( i == 0 ) ? 0 : 1e-4; // a little trunk between the branching points
            trunk += step;
            final String tip = names[ i ] + ":" + plain( Math.max( 1e-6, distances[ i ] - trunk ) );
            sb.insert( 0, ( i == 0 ) ? tip : ( "(" ) );
            if ( i > 0 ) {
                sb.append( "," ).append( tip ).append( "):" ).append( plain( step ) );
            }
        }
        return nh( "(" + sb + ");" );
    }

    private static String plain( final double d ) {
        return BigDecimal.valueOf( d ).toPlainString();
    }

    private static Phylogeny nh( final String s ) throws Exception {
        return ParserBasedPhylogenyFactory.getInstance().create( s, new NHXParser() )[ 0 ];
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [BranchLengthDateConverterTest] " + msg );
        return false;
    }

    private BranchLengthDateConverterTest() {
    }
}
