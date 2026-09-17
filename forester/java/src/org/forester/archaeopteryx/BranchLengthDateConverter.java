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
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/**
 * Dates a time tree that carries NO dates at all, when its BRANCH LENGTHS are years and its tip labels say so -- the
 * dateless sibling of {@link HeightDateConverter}.
 * <p>
 * TreeTime and Nextstrain both export a time tree as a plain Newick: the branch lengths are years, but nothing in the
 * file says that, and the same programs write divergence trees, in substitutions per site, in exactly the same shape.
 * The Nexus exports of those runs carry dates and open on the Calendar axis; their {@code .nwk} siblings used to open
 * as a bare distance tree. The tip labels settle it, as they do for BEAST heights: if the branch lengths are years,
 * then every tip's label date MINUS its distance from the root is the same calendar date -- the date of the root.
 * Where the labels agree ({@link LabelDateAnchor}), each node is dated {@code root date + its distance from the root},
 * in years, and the tree opens on the Calendar axis.
 * <p>
 * A divergence tree cannot pass. Its tips sit ~0.001 substitutions from the root while their labels span years, so
 * precise labels put every tip's offer somewhere different and nothing agrees. Coarse labels are caught by the rest of
 * the rule rather than by this one: year-only labels one year apart touch at a single date, which the
 * different-sampling-times test refuses, and labels three or more years apart cannot overlap at all. Measured over
 * every tree file in the forester, Archaeopteryx.js and Downloads corpora (2026-09-17), no divergence tree converts.
 * <p>
 * The dates this writes are DERIVED from the branch lengths, not measured, which is why the tree says so: the sentence
 * appended to the description names the root date, and every node's date follows from it and the branch lengths that
 * were already on screen. Because internal nodes end up dated, the tree is then a time tree -- so re-rooting is
 * refused, which is the point: a new root would contradict the dates the branch lengths imply.
 */
final class BranchLengthDateConverter {

    /** Dates every tree of a load that qualifies, appending the provenance sentence to each one's description. Returns
     *  the number of trees dated. Called by every load path, after {@link HeightDateConverter}. */
    static int dateTreesFromBranchLengths( final Phylogeny[] phys ) {
        if ( phys == null ) {
            return 0;
        }
        int dated = 0;
        for( final Phylogeny phy : phys ) {
            final LabelDateAnchor.Anchor anchor = inferRootDate( phy );
            if ( anchor == null ) {
                continue;
            }
            dateNodes( phy, anchor.value() );
            final String prov = provenanceSentence( anchor, phy.getName(), phy.getNumberOfExternalNodes() );
            final String existing = phy.getDescription();
            phy.setDescription( ForesterUtil.isEmpty( existing ) ? prov : ( existing + " " + prov ) );
            dated++;
        }
        return dated;
    }

    /** The calendar date of the ROOT, or null when the tree does not qualify: it already carries dates (any date at
     *  all -- the file has its own idea of time, and {@link HeightDateConverter} has had its turn), it has no branch
     *  lengths to read, or the tip labels do not agree. Pure. */
    static LabelDateAnchor.Anchor inferRootDate( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return null;
        }
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().getDate() != null ) {
                return null; // a dated tree is not this case, whatever its unit
            }
        }
        int tips = 0;
        final List<LabelDateAnchor.TipOffer> offers = new ArrayList<>();
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode tip = it.next();
            tips++;
            final double distance = distanceToRoot( tip );
            if ( !( distance > 0 ) ) {
                continue; // no branch lengths on this tip's path: it says nothing about time
            }
            final DateMatch m = TipDateExtractor.parse( tip.getName(), DayMonthOrder.DAY_FIRST );
            if ( m != null ) {
                // a distance is time AFTER the root, so the tip's label minus its distance points at the root
                offers.add( new LabelDateAnchor.TipOffer( m.rangeStart(), m.rangeEnd(), -distance, m.decimalYear() ) );
            }
        }
        return LabelDateAnchor.infer( offers, tips );
    }

    /** Dates every node {@code root_date + its distance from the root}, in years. */
    static void dateNodes( final Phylogeny phy, final BigDecimal root_date ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final BigDecimal value = LabelDateAnchor
                    .rounded( root_date.add( BigDecimal.valueOf( distanceToRoot( n ) ) ) );
            n.getNodeData().setDate( new Date( "", value, null, null, HeightDateConverter.YEAR_UNIT ) );
        }
    }

    /** e.g. <i>Dated tree named "timetree" with 99 tips from the sampling dates in 99 of 99 tip labels: the branch
     *  lengths are years, and the root is at 2012.976, so each node's date is that plus its distance from the
     *  root.</i> */
    static String provenanceSentence( final LabelDateAnchor.Anchor anchor, final String tree_name,
                                      final int num_ext_nodes ) {
        final String root = anchor.value().toPlainString();
        return "Dated " + TreePanelUtil.provenanceTreePhrase( tree_name, num_ext_nodes ) + " from the sampling dates in "
                + anchor.agreeing() + " of " + anchor.compared() + " tip labels: the branch lengths are years, and the "
                + "root is at " + root + ", so each node's date is that plus its distance from the root.";
    }

    /** Distance from the root, counting only branches whose length is set; 0 at the root, and 0 (i.e. "says nothing")
     *  for a tip whose path to the root carries no length at all. */
    private static double distanceToRoot( final PhylogenyNode node ) {
        double d = 0;
        PhylogenyNode n = node;
        while ( !n.isRoot() && ( n.getParent() != null ) ) {
            final double b = n.getDistanceToParent();
            if ( b != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                d += b;
            }
            n = n.getParent();
        }
        return d;
    }

    private BranchLengthDateConverter() {
    }
}
