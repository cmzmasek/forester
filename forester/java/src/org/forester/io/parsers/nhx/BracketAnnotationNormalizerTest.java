// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.io.parsers.nhx;

import java.util.List;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Headless tests for {@link BracketAnnotationNormalizer}: the three per-tree decisions (whose annotations, whether an
 * Auspice num_date stands, whether a TreeTime date= is promoted), their joint tolerances and their run order. The
 * fixtures mirror the real pairs that make the rules necessary: TreeTime writes the SAME date= on timetree.nexus and
 * divergence_tree.nexus, Auspice the SAME num_date on its time and divergence exports.
 */
public final class BracketAnnotationNormalizerTest {

    // TreeTime's shape: the blob AFTER the length. Years: every parent-to-child date difference IS the length.
    private static final String TT_YEARS      = "((A:1.00[&mutations=\"A1G,C5T\",date=2003.00],B:2.50[&date=2004.50])"
            + "NODE_0000001:1.00[&mutations=\"G7A\",date=2002.00],C:2.25[&date=2003.25])NODE_0000000:0.0[&date=2001.00];";
    // the SAME dates over substitution lengths
    private static final String TT_DIVERGENCE = "((A:0.0012[&mutations=\"A1G,C5T\",date=2003.00],B:0.0031[&date=2004.50])"
            + "NODE_0000001:0.0011[&mutations=\"G7A\",date=2002.00],C:0.0027[&date=2003.25])NODE_0000000:0.0[&date=2001.00];";
    // Auspice's shape: the blob BEFORE the length
    private static final String AU_YEARS      = "((A[&num_date=2003.0,num_date_CI={2002.9,2003.1},country=Ghana]:1.0,"
            + "B[&num_date=2004.5]:2.5)NODE_0000001[&num_date=2002.0,num_date_CI={2001.5,2002.5}]:1.0,"
            + "C[&num_date=2003.25]:2.25)NODE_0000000[&num_date=2001.0];";
    private static final String AU_DIVERGENCE = "((A[&num_date=2003.0,num_date_CI={2002.9,2003.1},country=Ghana]:0.0012,"
            + "B[&num_date=2004.5]:0.0031)NODE_0000001[&num_date=2002.0,num_date_CI={2001.5,2002.5}]:0.0011,"
            + "C[&num_date=2003.25,date=2003-04-01]:0.0027)NODE_0000000[&num_date=2001.0];";

    public static boolean test() {
        return testProducer() && testDatePromotion() && testPromotionThresholds() && testTolerances()
                && testNumDateSettles() && testNumDateBurdenOfProof() && testRunOrder() && testOptionOff()
                && testJointNumDateCases() && testDenseTreesUninformativePairs();
    }

    /** "Mutations present, age absent" -> treetime:, for EVERY beast: property of the tree. */
    private static boolean testProducer() {
        try {
            final Phylogeny tt = parse( "((A:1[&mutations=\"A1G\",mcc=3],B:1[&region=\"east\"]):1,C:1);" );
            if ( ( prop( tt.getNode( "A" ), "treetime:mutations" ) == null ) || ( prop( tt.getNode( "A" ), "treetime:mcc" ) == null )
                    || ( prop( tt.getNode( "B" ), "treetime:region" ) == null ) ) {
                return fail( "a mutations-only tree is TreeTime's: every beast: property becomes treetime:" );
            }
            if ( countPrefix( tt, "beast:" ) != 0 ) {
                return fail( "no beast: property may be left on a TreeTime tree" );
            }
            final Property mcc = prop( tt.getNode( "A" ), "treetime:mcc" );
            if ( !"3".equals( mcc.getValue() ) || !"xsd:string".equals( mcc.getDataType() ) ) {
                return fail( "a renamed property keeps its value and datatype" );
            }
            // a tree carrying BOTH mutations and an age stays beast:
            final Phylogeny both = parse( "((A:1[&mutations=\"A1G\"],B:1):1[&height=3.1],C:1);" );
            if ( ( prop( both.getNode( "A" ), "beast:mutations" ) == null ) || ( countPrefix( both, "treetime:" ) != 0 ) ) {
                return fail( "mutations beside a node age: not TreeTime's own, stays beast:" );
            }
            // an HPD alone states an age too
            final Phylogeny hpd = parse( "((A:1[&mutations=\"A1G\"],B:1):1[&height_95%_HPD={2.9,3.4}],C:1);" );
            if ( countPrefix( hpd, "treetime:" ) != 0 ) {
                return fail( "an age interval alone states an age: stays beast:" );
            }
            // ...and so does an age of exactly 0, all a BEAST tip ever states (NodeData.isHasDate() reads it as no date)
            final Phylogeny zero = parse( "((A:1[&mutations=\"A1G\",height=0.0],B:1[&height=0.0]):1,C:1[&height=0.0]);" );
            if ( countPrefix( zero, "treetime:" ) != 0 ) {
                return fail( "height=0.0 states an age: stays beast:" );
            }
            // bare mugration output can be attributed to nobody
            final Phylogeny mug = parse( "((A:1[&region=\"east\"],B:1[&region=\"west\"]):1[&region=\"east\"],C:1);" );
            if ( ( countPrefix( mug, "treetime:" ) != 0 ) || ( prop( mug.getNode( "A" ), "beast:region" ) == null ) ) {
                return fail( "a bare trait with no mutations anywhere keeps the generic namespace" );
            }
            // nextstrain: is never renamed
            final Phylogeny ns = parse( "((A:1[&mutations=\"A1G\",div=0.1],B:1):1,C:1);" );
            if ( ( prop( ns.getNode( "A" ), "nextstrain:div" ) == null ) || ( prop( ns.getNode( "A" ), "treetime:mutations" ) == null ) ) {
                return fail( "only beast: refs are renamed; nextstrain:div stays" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "producer parse threw: " + e );
        }
    }

    /** The headline pair: identical date= values, promoted on the years tree and NOT on the substitutions tree. */
    private static boolean testDatePromotion() {
        try {
            final Phylogeny years = parse( TT_YEARS );
            for( final PhylogenyNodeIterator it = years.iteratorPreorder(); it.hasNext(); ) {
                final Date d = it.next().getNodeData().getDate();
                if ( ( d == null ) || ( d.getValue() == null ) || !"year".equals( d.getUnit() ) ) {
                    return fail( "years tree: every date= must be promoted to a value with unit year, got " + d );
                }
                if ( d.getValue().doubleValue() != Double.parseDouble( d.getDesc() ) ) {
                    return fail( "years tree: the value is the desc's number and the desc STAYS, got " + d );
                }
            }
            if ( !"2003.00".equals( years.getNode( "A" ).getNodeData().getDate().getDesc() ) ) {
                return fail( "the desc is kept as written" );
            }
            final Phylogeny div = parse( TT_DIVERGENCE );
            for( final PhylogenyNodeIterator it = div.iteratorPreorder(); it.hasNext(); ) {
                final Date d = it.next().getNodeData().getDate();
                if ( ( d == null ) || ( d.getValue() != null ) || !isEmpty( d.getUnit() ) || isEmpty( d.getDesc() ) ) {
                    return fail( "divergence tree: the SAME dates must stay descs only, got " + d );
                }
            }
            // a calendar STRING is never a candidate, on any tree
            final Phylogeny str = parse( TT_YEARS.replace( "date=2003.00", "date=2003-01-01" ) );
            final Date ds = str.getNode( "A" ).getNodeData().getDate();
            if ( ( ds.getValue() != null ) || !"2003-01-01".equals( ds.getDesc() ) ) {
                return fail( "a non-numeric date= stays a desc even on a promoted tree, got " + ds );
            }
            if ( str.getNode( "B" ).getNodeData().getDate().getValue() == null ) {
                return fail( "...while the numeric ones around it are still promoted (3 of 3 remaining pairs agree)" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "date-promotion parse threw: " + e );
        }
    }

    /** Evidence FOR is needed: two comparable pairs at least, and a STRICT majority of them agreeing. */
    private static boolean testPromotionThresholds() {
        try {
            // ONE comparable pair (only A and its parent are dated)
            final Phylogeny one = parse( "((A:1.0[&date=2003.0],B:2.5):1.0[&date=2002.0],C:2.25);" );
            if ( one.getNode( "A" ).getNodeData().getDate().getValue() != null ) {
                return fail( "one comparable pair is no evidence: not promoted" );
            }
            // a dated node whose DIRECT parent is undated is no pair: here 0 pairs although 3 nodes are dated
            final Phylogeny gaps = parse( "((A:1.0[&date=2003.0],B:2.5[&date=2004.5]):1.0,C:2.25[&date=2003.25]);" );
            if ( gaps.getNode( "A" ).getNodeData().getDate().getValue() != null ) {
                return fail( "pairs are node + DIRECT parent; with none, nothing is promoted" );
            }
            // ...and a GRANDPARENT is no parent: A and B would each match the dated root two branches up (their own
            // length equals that date difference), but only C has a dated direct parent -> one pair -> not promoted
            final Phylogeny skip = parse( "((A:2.0[&date=2003.0],B:3.5[&date=2004.5]):1.0,C:2.25[&date=2003.25]):0.0[&date=2001.0];" );
            if ( skip.getNode( "A" ).getNodeData().getDate().getValue() != null ) {
                return fail( "a pair never skips a generation" );
            }
            // 4 pairs, exactly 2 agree: no STRICT majority
            final Phylogeny half = parse( TT_YEARS.replace( "A:1.00", "A:0.10" ).replace( "C:2.25", "C:0.25" ) );
            if ( half.getNode( "B" ).getNodeData().getDate().getValue() != null ) {
                return fail( "2 of 4 pairs agreeing is not a strict majority: not promoted" );
            }
            // 4 pairs, 3 agree: promoted, ALL of them (the disagreeing node too)
            final Phylogeny most = parse( TT_YEARS.replace( "A:1.00", "A:0.10" ) );
            if ( ( most.getNode( "B" ).getNodeData().getDate().getValue() == null )
                    || ( most.getNode( "A" ).getNodeData().getDate().getValue() == null ) ) {
                return fail( "3 of 4 pairs agreeing is a strict majority: every numeric date is promoted" );
            }
            // a pair needs a branch length
            final Phylogeny no_len = parse( "((A[&date=2003.0],B[&date=2004.5])[&date=2002.0],C[&date=2003.25])[&date=2001.0];" );
            if ( no_len.getNode( "A" ).getNodeData().getDate().getValue() != null ) {
                return fail( "without branch lengths there is nothing to confirm the dates: not promoted" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "threshold parse threw: " + e );
        }
    }

    /** JOINT tolerances: |dDate - length| <= 0.02 + 0.01 x |length|. */
    private static boolean testTolerances() {
        if ( ( BracketAnnotationNormalizer.NUMERIC_DATE_ABS_TOL != 0.02 )
                || ( BracketAnnotationNormalizer.NUMERIC_DATE_REL_TOL != 0.01 ) ) {
            return fail( "the tolerances are JOINT constants (0.02 abs, 0.01 rel): neither side retunes alone" );
        }
        try {
            // every pair is off by 0.029 on a length of 1.0: inside 0.02 + 0.01 x 1.0 = 0.03
            final Phylogeny inside = parse( "((A:1.0[&date=2003.029],B:1.0[&date=2003.029]):1.0[&date=2002.0],C:1.0[&date=2001.971])"
                    + ":0.0[&date=2001.0];" );
            if ( inside.getNode( "A" ).getNodeData().getDate().getValue() == null ) {
                return fail( "an error of 0.029 on a length of 1.0 is inside the tolerance" );
            }
            // off by 0.04: outside (A, B and C disagree; only the inner node agrees)
            final Phylogeny outside = parse( "((A:1.0[&date=2003.04],B:1.0[&date=2003.04]):1.0[&date=2002.0],C:1.0[&date=2001.96])"
                    + ":0.0[&date=2001.0];" );
            if ( outside.getNode( "A" ).getNodeData().getDate().getValue() != null ) {
                return fail( "an error of 0.04 on a length of 1.0 is outside the tolerance" );
            }
            // the relative part: off by 0.9 on a length of 100.5 is inside 0.02 + 0.01 x 100.5 = 1.025
            // (the tip lengths differ: tips that ALL share one length >= 10 read as bootstrap values and leave the
            // branches -- NHXParser.isBranchLengthsLikeBootstrapValues)
            final Phylogeny rel = parse( "((A:100.5[&date=2201.9],B:80.5[&date=2181.7]):100.5[&date=2100.5],"
                    + "C:120.5[&date=2119.6]):0.0[&date=2000.0];" );
            if ( rel.getNode( "A" ).getNodeData().getDate().getValue() == null ) {
                return fail( "the tolerance grows with the branch length (0.01 relative)" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "tolerance parse threw: " + e );
        }
    }

    /** Auspice writes num_date on its divergence export too (nextstrain_lassa_gpc_tree.nexus: 169 of 2295 pairs
     *  agree). There the date FALLS -- no calendar axis over substitutions, re-rooting not refused -- and nothing
     *  is lost: the year stays as nextstrain:num_date, the interval as nextstrain:num_date_CI. */
    private static boolean testNumDateSettles() {
        try {
            final Phylogeny years = parse( AU_YEARS );
            final Date da = years.getNode( "A" ).getNodeData().getDate();
            if ( ( da == null ) || ( da.getValue() == null ) || ( da.getValue().doubleValue() != 2003.0 ) || !"year".equals( da.getUnit() ) ) {
                return fail( "time-scaled Auspice tree: num_date stands, got " + da );
            }
            // ...whole: a TIP keeps its interval (a sample dated only to the month or year has a real one), as the
            // INTERNAL node does -- what AuspiceJsonParser makes of the same build
            if ( ( da.getMin() == null ) || ( da.getMin().doubleValue() != 2002.9 ) || ( da.getMax() == null )
                    || ( da.getMax().doubleValue() != 2003.1 ) ) {
                return fail( "a standing TIP keeps its num_date_CI interval, got " + da.getMin() + " .. " + da.getMax() );
            }
            final Date di = years.getNode( "NODE_0000001" ).getNodeData().getDate();
            if ( ( di == null ) || ( di.getMin() == null ) || ( di.getMin().doubleValue() != 2001.5 ) || ( di.getMax() == null )
                    || ( di.getMax().doubleValue() != 2002.5 ) ) {
                return fail( "a standing INTERNAL node keeps its num_date_CI interval, got " + di );
            }
            if ( ( prop( years.getNode( "A" ), "nextstrain:num_date_CI" ) != null )
                    || ( prop( years.getNode( "NODE_0000001" ), "nextstrain:num_date_CI" ) != null ) ) {
                return fail( "where the dates stand an interval lives IN the date, never as a property as well" );
            }
            final Phylogeny div = parse( AU_DIVERGENCE );
            final PhylogenyNode a = div.getNode( "A" );
            if ( a.getNodeData().isHasDate() ) {
                return fail( "divergence Auspice tree: num_date must fall, got " + a.getNodeData().getDate() );
            }
            final Property nd = prop( a, "nextstrain:num_date" );
            final Property ci = prop( a, "nextstrain:num_date_CI" );
            if ( ( nd == null ) || !"2003.0".equals( nd.getValue() ) || !"xsd:decimal".equals( nd.getDataType() ) ) {
                return fail( "a fallen num_date stays as the numeric nextstrain:num_date, got " + nd );
            }
            if ( ( ci == null ) || !"{2002.9,2003.1}".equals( ci.getValue() ) || !"xsd:string".equals( ci.getDataType() ) ) {
                return fail( "a fallen num_date's interval stays as nextstrain:num_date_CI text, got "
                        + ( ci == null ? "null" : ci.getValue() ) );
            }
            if ( prop( div.getNode( "NODE_0000001" ), "nextstrain:num_date_CI" ) == null ) {
                return fail( "a fallen INTERNAL node's interval is re-filed too" );
            }
            if ( prop( div.getNode( "B" ), "nextstrain:num_date_CI" ) != null ) {
                return fail( "a node that had no interval gets no interval property" );
            }
            if ( prop( a, "beast:country" ) == null ) {
                return fail( "the node's other properties are untouched" );
            }
            // a date= desc beside a fallen num_date survives as the desc it always was
            final Date dc = div.getNode( "C" ).getNodeData().getDate();
            if ( ( dc == null ) || ( dc.getValue() != null ) || !isEmpty( dc.getUnit() ) || !"2003-04-01".equals( dc.getDesc() ) ) {
                return fail( "a desc beside a fallen num_date stays, got " + dc );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "num_date parse threw: " + e );
        }
    }

    /** The burden of proof is REVERSED from date=: a num_date is a date by name, so it stands unless there is
     *  evidence AGAINST -- a tree too small to judge keeps its dates; exactly half agreeing is no majority. */
    private static boolean testNumDateBurdenOfProof() {
        try {
            // ONE comparable pair, and it disagrees: too small to judge -> stands
            final Phylogeny one = parse( "((A[&num_date=2003.0]:0.001,B:0.002)[&num_date=2002.0]:0.001,C:0.002);" );
            if ( !one.getNode( "A" ).getNodeData().isHasDate() || ( one.getNode( "A" ).getNodeData().getDate().getValue() == null ) ) {
                return fail( "one disagreeing pair is no evidence against: the num_date stands" );
            }
            // a single dated tip: stands
            final Phylogeny single = parse( "(A[&num_date=2003.0]:0.001,B:0.002);" );
            if ( !single.getNode( "A" ).getNodeData().isHasDate() ) {
                return fail( "no pair at all: the num_date stands" );
            }
            // 4 pairs, exactly 2 agree: NOT a strict majority -> falls
            final Phylogeny half = parse( AU_YEARS.replace( "country=Ghana]:1.0", "country=Ghana]:0.1" ).replace( "]:2.25", "]:0.25" ) );
            if ( half.getNode( "B" ).getNodeData().isHasDate() ) {
                return fail( "2 of 4 pairs agreeing is no strict majority: the num_dates fall" );
            }
            // 4 pairs, 3 agree -> stand, ALL of them
            final Phylogeny most = parse( AU_YEARS.replace( "country=Ghana]:1.0", "country=Ghana]:0.1" ) );
            if ( !most.getNode( "A" ).getNodeData().isHasDate() || !most.getNode( "B" ).getNodeData().isHasDate() ) {
                return fail( "3 of 4 pairs agreeing: the num_dates stand, the disagreeing node's too" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "burden-of-proof parse threw: " + e );
        }
    }

    /** The producer test runs FIRST: a num_date says "not TreeTime's own" while it is still provisional -- even
     *  one that is about to fall. And a BEAST height is never a num_date: it is not settled, whatever the lengths. */
    private static boolean testRunOrder() {
        try {
            final Phylogeny p = parse( AU_DIVERGENCE.replace( "country=Ghana", "mutations=\"A1G\"" ) );
            if ( p.getNode( "A" ).getNodeData().isHasDate() ) {
                return fail( "fixture: the num_dates of this tree must have fallen" );
            }
            if ( ( prop( p.getNode( "A" ), "beast:mutations" ) == null ) || ( countPrefix( p, "treetime:" ) != 0 ) ) {
                return fail( "a tree with num_date is Auspice's, not TreeTime's own -- even where the num_date falls" );
            }
            final Phylogeny beast = parse( "((A:0.001[&height=0.0,height_95%_HPD={0.0,0.2}],B:0.002[&height=0.0]):0.001[&height=1.5,height_95%_HPD={1.0,2.0}],"
                    + "C:0.002[&height=0.0])[&height=3.0];" );
            // (4 comparable pairs, none agreeing: exactly what makes a num_date fall)
            final Date d = beast.getNode( "A" ).getParent().getNodeData().getDate();
            if ( ( d == null ) || ( d.getValue() == null ) || ( d.getValue().doubleValue() != 1.5 ) || ( d.getMin() == null )
                    || !isEmpty( d.getUnit() ) ) {
                return fail( "a BEAST height is an age, not a num_date: never settled, never given a unit, got " + d );
            }
            final Date tip = beast.getNode( "A" ).getNodeData().getDate();
            if ( ( tip == null ) || ( tip.getValue() == null ) || ( tip.getValue().doubleValue() != 0.0 ) ) {
                return fail( "a BEAST tip's height of 0.0 stays its date value, got " + tip );
            }
            if ( ( tip.getMax() == null ) || ( tip.getMax().doubleValue() != 0.2 ) ) {
                return fail( "a BEAST tip's height HPD is a sampled tip date and is kept" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "run-order parse threw: " + e );
        }
    }

    /** The JOINT acceptance cases, as Archaeopteryx.js runs them (its testNumDateOnlyOnTimeScaledTrees): one tree,
     *  four branch lengths; year differences are A 2, B 3.5, ab 1, C 6 -- four comparable pairs. */
    private static boolean testJointNumDateCases() {
        final String[][] lengths = { { "2", "3.5", "1", "6" }, { "0.002", "0.0035", "0.001", "0.006" },
                { "2", "3.5", "0.001", "0.006" }, { "2", "3.5", "1", "0.006" } };
        final int[] expected_dated = { 5, 0, 0, 5 };
        try {
            for( int i = 0; i < lengths.length; ++i ) {
                final String[] l = lengths[ i ];
                final Phylogeny phy = parse( "((A:" + l[ 0 ] + "[&num_date=2003,num_date_CI={2002.5,2003.5}],B:" + l[ 1 ]
                        + "[&num_date=2004.5])ab:" + l[ 2 ] + "[&num_date=2001,num_date_CI={2000.5,2001.5}],C:" + l[ 3 ]
                        + "[&num_date=2006])root[&num_date=2000];" );
                int dated = 0;
                int num_date = 0;
                int ci = 0;
                for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if ( n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getValue() != null ) ) {
                        ++dated;
                    }
                    num_date += ( prop( n, "nextstrain:num_date" ) != null ) ? 1 : 0;
                    ci += ( prop( n, "nextstrain:num_date_CI" ) != null ) ? 1 : 0;
                }
                if ( ( dated != expected_dated[ i ] ) || ( num_date != 5 ) || ( ci != ( ( expected_dated[ i ] == 0 ) ? 2 : 0 ) ) ) {
                    return fail( "joint case " + i + ": expected " + expected_dated[ i ] + " dated, 5 nextstrain:num_date, got "
                            + dated + " dated, " + num_date + " num_date, " + ci + " num_date_CI" );
                }
                final PhylogenyNode a = phy.getNode( "A" );
                if ( expected_dated[ i ] == 5 ) {
                    final Date d = a.getNodeData().getDate();
                    if ( ( d.getValue().doubleValue() != 2003 ) || !"year".equals( d.getUnit() ) || ( d.getMin() == null )
                            || ( d.getMin().doubleValue() != 2002.5 ) || ( d.getMax().doubleValue() != 2003.5 ) ) {
                        return fail( "joint case " + i + ": tip A is {2003, year} with min 2002.5 / max 2003.5" );
                    }
                    final Date dab = phy.getNode( "ab" ).getNodeData().getDate();
                    if ( ( dab.getMin() == null ) || ( dab.getMin().doubleValue() != 2000.5 ) || ( dab.getMax().doubleValue() != 2001.5 ) ) {
                        return fail( "joint case " + i + ": internal ab has min 2000.5 / max 2001.5" );
                    }
                }
                else if ( a.getNodeData().isHasDate() ) {
                    return fail( "joint case " + i + ": A has NO date element" );
                }
            }
            final Phylogeny two = parse( "(A:0.002[&num_date=2003],B:0.004[&num_date=2004])root[&num_date=2000];" );
            if ( two.getNode( "A" ).getNodeData().isHasDate() || two.getNode( "B" ).getNodeData().isHasDate() ) {
                return fail( "joint: two pairs, neither agrees -> 0 date values" );
            }
            final Phylogeny one = parse( "(A:0.002[&num_date=2003],B:0.004)root[&num_date=2000];" );
            if ( !one.getNode( "A" ).getNodeData().isHasDate() || !one.getRoot().getNodeData().isHasDate() ) {
                return fail( "joint: ONE pair is not evidence -> both dated" );
            }
            final Date lone_tip = parse( "(A:1[&num_date=2001.5,num_date_CI={2000.1,2002.9}],B:1);" ).getNode( "A" ).getNodeData().getDate();
            if ( ( lone_tip == null ) || ( lone_tip.getValue().doubleValue() != 2001.5 ) || ( lone_tip.getMin() == null )
                    || ( lone_tip.getMin().doubleValue() != 2000.1 ) || ( lone_tip.getMax().doubleValue() != 2002.9 ) ) {
                return fail( "joint: a lone dated tip keeps min 2000.1 / max 2002.9" );
            }
            final Date lone_internal = parse( "((X:1,Y:1)A:1[&num_date=2001.5,num_date_CI={2000.1,2002.9}],B:1);" ).getNode( "A" )
                    .getNodeData().getDate();
            if ( ( lone_internal == null ) || ( lone_internal.getMin() == null ) || ( lone_internal.getMin().doubleValue() != 2000.1 )
                    || ( lone_internal.getMax().doubleValue() != 2002.9 ) ) {
                return fail( "joint: a lone dated INTERNAL node keeps min 2000.1 / max 2002.9" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "joint-case parse threw: " + e );
        }
    }

    /** A DENSELY sampled tree: most parent-to-child date differences are below the 0.02 absolute tolerance. On
     *  substitution lengths such a pair "agrees" by accident, so it is not counted at all -- only a pair whose date
     *  difference or length exceeds the tolerance can tell years from substitutions. Here 4 of 6 pairs are
     *  uninformative (and would agree), 2 are informative. */
    private static boolean testDenseTreesUninformativePairs() {
        // {A, B, N1, C, D, N2} lengths: divergence (substitutions) and time (the date differences themselves)
        final String div = "((A[&D=2020.02]:0.00001,B[&D=2020.51]:0.0001)N1[&D=2020.01]:0.00002,"
                + "(C[&D=2020.52]:0.00011,D[&D=2020.02]:0.00003)N2[&D=2020.01]:0.00004)R[&D=2020.00];";
        final String time = "((A[&D=2020.02]:0.01,B[&D=2020.51]:0.5)N1[&D=2020.01]:0.011,"
                + "(C[&D=2020.52]:0.51,D[&D=2020.02]:0.0105)N2[&D=2020.01]:0.0101)R[&D=2020.00];";
        try {
            // Auspice num_date: on the divergence tree the dates FALL (2 informative pairs, 0 agree) -- counted
            // naively 4 of 6 would agree, a majority, and the tree would be dated over substitutions
            final Phylogeny ns_div = parse( div.replace( "[&D=", "[&num_date=" ) );
            if ( ns_div.getNode( "A" ).getNodeData().isHasDate() ) {
                return fail( "dense divergence tree: uninformative pairs must not keep the num_dates standing" );
            }
            final Phylogeny ns_time = parse( time.replace( "[&D=", "[&num_date=" ) );
            if ( !ns_time.getNode( "A" ).getNodeData().isHasDate() ) {
                return fail( "dense TIME tree: its informative pairs agree, the num_dates stand" );
            }
            // TreeTime date=: not promoted on the divergence tree, promoted on the time tree
            final Phylogeny tt_div = parse( div.replace( "[&D=", "[&mutations=\"A1G\",date=" ) );
            if ( tt_div.getNode( "B" ).getNodeData().getDate().getValue() != null ) {
                return fail( "dense divergence tree: uninformative pairs must not promote date=" );
            }
            final Phylogeny tt_time = parse( time.replace( "[&D=", "[&mutations=\"A1G\",date=" ) );
            if ( tt_time.getNode( "B" ).getNodeData().getDate().getValue() == null ) {
                return fail( "dense TIME tree: date= is promoted on its informative pairs" );
            }
            // a tree with NO informative pair at all: date= needs evidence FOR (not promoted), num_date needs evidence
            // AGAINST (stands) -- the joint burden of proof, unchanged
            final String flat = "((A[&D=2020.02]:0.00001,B[&D=2020.03]:0.00002)N1[&D=2020.01]:0.00003,C[&D=2020.015]:0.00004)R[&D=2020.00];";
            if ( parse( flat.replace( "[&D=", "[&mutations=\"A1G\",date=" ) ).getNode( "A" ).getNodeData().getDate().getValue() != null ) {
                return fail( "no informative pair: date= stays a description" );
            }
            if ( !parse( flat.replace( "[&D=", "[&num_date=" ) ).getNode( "A" ).getNodeData().isHasDate() ) {
                return fail( "no informative pair: a num_date stands" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "dense-tree parse threw: " + e );
        }
    }

    private static boolean testOptionOff() {
        try {
            final NHXParser p = new NHXParser();
            p.setParseBeastStyleExtendedTags( false );
            p.setSource( TT_YEARS );
            final Phylogeny phy = p.parse()[ 0 ];
            if ( phy.getNode( "A" ).getNodeData().isHasDate() || ( countPrefix( phy, "treetime:" ) != 0 ) ) {
                return fail( "with the bracket-annotation option OFF nothing is normalized" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "option-off parse threw: " + e );
        }
    }

    private static Phylogeny parse( final String nh ) throws Exception {
        final NHXParser p = new NHXParser();
        p.setParseBeastStyleExtendedTags( true );
        p.setSource( nh );
        return p.parse()[ 0 ];
    }

    private static int countPrefix( final Phylogeny phy, final String prefix ) {
        int c = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().getProperties() != null ) {
                for( final Property p : n.getNodeData().getProperties().getProperties() ) {
                    if ( p.getRef().startsWith( prefix ) ) {
                        ++c;
                    }
                }
            }
        }
        return c;
    }

    private static Property prop( final PhylogenyNode node, final String ref ) {
        if ( node.getNodeData().getProperties() == null ) {
            return null;
        }
        final List<Property> ps = node.getNodeData().getProperties().getProperties( ref );
        return ps.isEmpty() ? null : ps.get( 0 );
    }

    private static boolean isEmpty( final String s ) {
        return ( s == null ) || s.isEmpty();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "BracketAnnotationNormalizer test failed: " + msg );
        return false;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }
}
