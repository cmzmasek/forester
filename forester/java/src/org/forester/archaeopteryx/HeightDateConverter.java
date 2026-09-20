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
 * <li>the tips agree on where height 0 is, by the shared rule in {@link LabelDateAnchor}: each offers its label's
 * calendar range plus its height, at least 19 in 20 must allow one date, no rival stretch may be as well supported,
 * and the agreeing tips must have been sampled at different times.</li>
 * </ol>
 * Height 0 is then the median of label date plus height over the most precisely dated agreeing tips, kept inside the
 * dates all of them allow, and each node's date is that minus its height (the older HPD bound becomes the earlier
 * date), in years, rounded to 5 decimals. The heights are not kept: the sentence appended to the tree description
 * names the date of height 0, from which each one follows. (Keeping them as a {@code beast:height} property was
 * measured and rejected: the third-party Ebola example tree, which has no other Color-by field, would open
 * coloured by it.)
 * <p>
 * Measured 2026-09-17 on every tree file in the forester, Archaeopteryx.js and Downloads corpora: the four real BEAST
 * trees (influenza.tree, HA_discrete_MCC, HA_continuous_MCC, and the third-party Ebola example named in the
 * working notes) agree on every tip, and no
 * other tree has unit-less dates.
 */
final class HeightDateConverter {

    static final String YEAR_UNIT = "year";

    /** Converts every tree of a load that qualifies, appending the provenance sentence to each converted tree's
     *  description. Returns the number of trees converted. Called by every load path, next to
     *  {@link AptxUtil#applyInternalLabelPolicy}. */
    static int convertHeightsToDates( final Phylogeny[] phys ) {
        if ( phys == null ) {
            return 0;
        }
        int converted = 0;
        for( final Phylogeny phy : phys ) {
            final LabelDateAnchor.Anchor anchor = inferAnchor( phy );
            if ( anchor == null ) {
                continue;
            }
            convert( phy, anchor.value() );
            final String prov = provenanceSentence( anchor, phy.getName(), phy.getNumberOfExternalNodes() );
            final String existing = phy.getDescription();
            phy.setDescription( ForesterUtil.isEmpty( existing ) ? prov : ( existing + " " + prov ) );
            converted++;
        }
        return converted;
    }

    /** The calendar date of height 0 by the rule in the class comment, or null when the tree does not qualify. Pure. */
    static LabelDateAnchor.Anchor inferAnchor( final Phylogeny phy ) {
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
        final List<LabelDateAnchor.TipOffer> offers = new ArrayList<>();
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
                // a height is time BEFORE the anchor, so the tip's label plus its height points at the anchor
                offers.add( new LabelDateAnchor.TipOffer( m.rangeStart(), m.rangeEnd(), d.getValue().doubleValue(),
                                                          m.decimalYear() ) );
            }
        }
        return LabelDateAnchor.infer( offers, tips );
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
            final BigDecimal value = ( d.getValue() == null ) ? null
                    : LabelDateAnchor.rounded( present.subtract( d.getValue() ) );
            final BigDecimal min = ( d.getMax() == null ) ? null
                    : LabelDateAnchor.rounded( present.subtract( d.getMax() ) );
            final BigDecimal max = ( d.getMin() == null ) ? null
                    : LabelDateAnchor.rounded( present.subtract( d.getMin() ) );
            n.getNodeData().setDate( new Date( d.getDesc(), value, min, max, ( value != null ) ? YEAR_UNIT : "" ) );
        }
    }

    /** e.g. <i>Converted the node heights of tree named "TREE1" with 190 tips to calendar dates: the sampling dates in
     *  190 of 190 tip labels put height 0 at 2005.5, so each date is 2005.5 minus the height.</i> */
    static String provenanceSentence( final LabelDateAnchor.Anchor anchor, final String tree_name,
                                      final int num_ext_nodes ) {
        final String present = anchor.value().toPlainString();
        return "Converted the node heights of " + TreePanelUtil.provenanceTreePhrase( tree_name, num_ext_nodes )
                + " to calendar dates: the sampling dates in " + anchor.agreeing() + " of " + anchor.compared()
                + " tip labels put height 0 at " + present + ", so each date is " + present + " minus the height.";
    }

    private HeightDateConverter() {
    }
}
