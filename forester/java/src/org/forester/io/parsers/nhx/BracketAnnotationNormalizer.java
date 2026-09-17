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

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.json.AuspiceJsonParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/**
 * What {@link BeastAnnotationParser} cannot decide one node at a time, decided once the whole tree has been read.
 * Each question is put to the TREE, never guessed from a file name or sniffed from one annotation -- two producers
 * write the SAME annotation on two different trees:
 * <ul>
 * <li><b>Whose annotations are these?</b> A TreeTime tree carries {@code mutations=} and no node age at all, where
 *     every BEAST / MrBayes run states one. "Mutations present, age absent" re-namespaces the tree's
 *     {@code beast:} properties to {@code treetime:}. A bare mugration trait ({@code [&region="east"]}) can be
 *     attributed to nobody and keeps the generic namespace.</li>
 * <li><b>Is an Auspice {@code num_date} a date here?</b> Auspice's "download Nexus" writes it on its time tree
 *     (branch lengths in years) AND on its divergence tree (substitutions). A date value is what makes a tree a
 *     time tree -- calendar axis, re-rooting refused -- so it stands only where the tree does not contradict it.</li>
 * <li><b>Is a bare numeric {@code date=} a date value?</b> TreeTime writes the same {@code date=2003.84} on
 *     timetree.nexus and on divergence_tree.nexus. It becomes a value only where the tree confirms it.</li>
 * </ul>
 * The test is the same for both dates: do the parent-to-child date differences reproduce the branch lengths? Only
 * the burden of proof differs. A {@code date=} may not be a value at all, so it needs evidence FOR (two comparable
 * pairs or more, a strict majority agreeing); a {@code num_date} is a date by its very name, so it stands unless
 * there is evidence AGAINST (two comparable pairs or more, and NO strict majority agreeing) -- a tree too small to
 * say anything keeps its dates. Nothing is lost either way: the desc always stays, and a fallen {@code num_date}
 * stays as the numeric {@code nextstrain:num_date} with its interval as {@code nextstrain:num_date_CI}. Where the
 * num_dates stand they stand whole -- a TIP keeps its interval too (the sampling-date uncertainty of a sample dated
 * only to the month or the year), exactly as {@link AuspiceJsonParser} keeps it in the same build's JSON, so the two
 * downloads of one dataset open alike.
 * <p>
 * JOINT with Archaeopteryx.js ({@code renameTreeTimeProperties}, {@code settleNumDates},
 * {@code promoteTimeScaledDates} in forester.js), including the run order and both tolerances: neither side
 * retunes alone.
 */
public final class BracketAnnotationNormalizer {

    /** JOINT: {@code date=} is written to two decimals, so a difference of two carries about 0.01 of error. */
    static final double         NUMERIC_DATE_ABS_TOL = 0.02;
    /** JOINT. */
    static final double         NUMERIC_DATE_REL_TOL = 0.01;
    static final String         TREETIME_PREFIX      = "treetime:";
    private static final String BEAST_PREFIX         = "beast:";
    private static final String MUTATIONS_REF        = BEAST_PREFIX + "mutations";
    static final String         NUM_DATE_CI_REF      = AuspiceJsonParser.PREFIX + "num_date_CI";

    private BracketAnnotationNormalizer() {
        // pure utility
    }

    /** All three, in the one order that works: the producer test must see a provisional {@code num_date} (it says
     *  "not TreeTime's own") before that is settled, and both must run before a {@code date=} is promoted (until
     *  then a date VALUE can only have come from a height or a num_date). */
    public static void normalize( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return;
        }
        renameTreeTimeProperties( phy );
        settleNumDates( phy );
        promoteTimeScaledDates( phy );
    }

    static void renameTreeTimeProperties( final Phylogeny phy ) {
        boolean mutations = false;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            // (not isHasDate(): that reads a value of exactly 0 -- every BEAST tip's height -- as "no date")
            final Date d = n.getNodeData().getDate();
            if ( d != null ) {
                if ( ( d.getValue() != null ) || ( d.getMin() != null ) || ( d.getMax() != null ) ) {
                    return; // a node states an age (or a num_date): not TreeTime's own Nexus
                }
            }
            if ( !mutations && ( n.getNodeData().getProperties() != null )
                    && !n.getNodeData().getProperties().getProperties( MUTATIONS_REF ).isEmpty() ) {
                mutations = true;
            }
        }
        if ( !mutations ) {
            return;
        }
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final PropertiesList old = n.getNodeData().getProperties();
            if ( ( old == null ) || ( old.size() == 0 ) ) {
                continue;
            }
            final PropertiesList renamed = new PropertiesList(); // a Property's ref is final: rebuild the list
            for( final Property p : old.getProperties() ) {
                renamed.addProperty( p.getRef().startsWith( BEAST_PREFIX )
                        ? new Property( TREETIME_PREFIX + p.getRef().substring( BEAST_PREFIX.length() ),
                                        p.getValue(),
                                        p.getUnit(),
                                        p.getDataType(),
                                        p.getAppliesTo(),
                                        p.getIdRef() )
                        : p );
            }
            n.getNodeData().setProperties( renamed );
        }
    }

    /** A {@code num_date} date (the only source of unit "year" before {@link #promoteTimeScaledDates}) falls where
     *  the tree is measurably NOT time-scaled. */
    static void settleNumDates( final Phylogeny phy ) {
        final List<PhylogenyNode> marked = new ArrayList<PhylogenyNode>();
        final int[] pairs_agree = { 0, 0 };
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final Double year = numDateYear( n );
            if ( year != null ) {
                marked.add( n );
                count( n, year, n.isRoot() ? null : numDateYear( n.getParent() ), pairs_agree );
            }
        }
        if ( ( pairs_agree[ 0 ] < 2 ) || ( ( pairs_agree[ 1 ] * 2 ) > pairs_agree[ 0 ] ) ) {
            // Too small to judge, or time-scaled: the dates stand, intervals and all -- a TIP's too. A sample dated only
            // to the month or the year carries a genuine interval, and AuspiceJsonParser keeps the same build's JSON
            // the same way (Christian, 2026-09-17; both readers used to drop tip intervals, to keep the geologic
            // Fossil Range Bars off viral trees -- the display now tells the two apart instead).
            return;
        }
        for( final PhylogenyNode n : marked ) {
            final Date d = n.getNodeData().getDate();
            if ( ( d.getMin() != null ) && ( d.getMax() != null ) ) {
                PropertiesList pl = n.getNodeData().getProperties();
                if ( pl == null ) {
                    pl = new PropertiesList();
                    n.getNodeData().setProperties( pl );
                }
                pl.addProperty( new Property( NUM_DATE_CI_REF,
                                              "{" + d.getMin() + "," + d.getMax() + "}",
                                              "",
                                              "xsd:string",
                                              AppliesTo.NODE ) );
            }
            if ( ForesterUtil.isEmpty( d.getDesc() ) ) {
                n.getNodeData().setDate( null );
            }
            else {
                n.getNodeData().setDate( new Date( d.getDesc() ) );
            }
        }
    }

    /** A bare numeric {@code date=} desc becomes the date VALUE (unit "year") where the tree confirms it. The desc
     *  stays. */
    static void promoteTimeScaledDates( final Phylogeny phy ) {
        final List<PhylogenyNode> dated = new ArrayList<PhylogenyNode>();
        final int[] pairs_agree = { 0, 0 };
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final Double year = numericDateDesc( n );
            if ( year != null ) {
                dated.add( n );
                count( n, year, n.isRoot() ? null : numericDateDesc( n.getParent() ), pairs_agree );
            }
        }
        if ( ( pairs_agree[ 0 ] < 2 ) || ( ( pairs_agree[ 1 ] * 2 ) <= pairs_agree[ 0 ] ) ) {
            return;
        }
        for( final PhylogenyNode n : dated ) {
            final Date d = n.getNodeData().getDate();
            d.setValue( new BigDecimal( d.getDesc().trim() ) );
            d.setUnit( BeastAnnotationParser.YEAR_UNIT );
        }
    }

    /** One comparable pair = a dated node whose DIRECT parent is dated too and which has a branch length, and which
     *  can tell years from substitutions at all: its date difference or its length exceeds the absolute tolerance. A
     *  pair with both below it "agrees" whatever the branch lengths measure -- on a densely sampled DIVERGENCE tree
     *  most pairs are like that (rebuilt from a real Auspice H5N1 export, 39% of all pairs agreed on substitution
     *  lengths), which came close to dating a tree measured in substitutions. Informative pairs separate the two
     *  cleanly: that H5N1 divergence tree 1 of 5586, its time tree 5586 of 5586. A pair agrees when the date difference
     *  reproduces the length within the joint tolerances. JOINT (Christian, 2026-09-17). */
    private static void count( final PhylogenyNode n,
                               final Double year,
                               final Double parent_year,
                               final int[] pairs_agree ) {
        final double length = n.getDistanceToParent();
        if ( ( parent_year == null ) || ( length == PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) ) {
            return;
        }
        if ( Math.max( Math.abs( year.doubleValue() - parent_year.doubleValue() ), Math.abs( length ) ) <= NUMERIC_DATE_ABS_TOL ) {
            return; // uninformative: agrees on any scale
        }
        ++pairs_agree[ 0 ];
        final double tol = NUMERIC_DATE_ABS_TOL + ( NUMERIC_DATE_REL_TOL * Math.abs( length ) );
        if ( Math.abs( ( year.doubleValue() - parent_year.doubleValue() ) - length ) <= tol ) {
            ++pairs_agree[ 1 ];
        }
    }

    private static Double numDateYear( final PhylogenyNode n ) {
        final Date d = n.getNodeData().getDate();
        if ( d == null ) {
            return null;
        }
        return ( ( d.getValue() != null ) && BeastAnnotationParser.YEAR_UNIT.equals( d.getUnit() ) )
                ? Double.valueOf( d.getValue().doubleValue() ) : null;
    }

    private static Double numericDateDesc( final PhylogenyNode n ) {
        final Date d = n.getNodeData().getDate();
        if ( ( d == null ) || ( d.getValue() != null ) || ForesterUtil.isEmpty( d.getDesc() ) ) {
            return null;
        }
        return BeastAnnotationParser.parseNumber( d.getDesc().trim() );
    }
}
