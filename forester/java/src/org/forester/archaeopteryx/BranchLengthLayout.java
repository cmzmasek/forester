// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Laying a dated tree out by TIME or by DIVERGENCE, for the trees that carry both.
 * <p>
 * A time-calibrated tree usually knows two different things about every branch: how much TIME it spans, and how much
 * genetic change happened along it. Only one can be the x-axis at a time, so this switches between them. It is a pure
 * DISPLAY mode: both quantities stay retained on the tree, nothing is edited, and switching back is exact.
 * <p>
 * Two shapes are recognised, and the difference between them is not cosmetic:
 * <ul>
 * <li>{@link DIVERGENCE_SOURCE#STORED} -- the file RECORDS divergence per node, as a cumulative value (Auspice /
 * Nextstrain {@code nextstrain:div}). Branch divergence is the successive difference.</li>
 * <li>{@link DIVERGENCE_SOURCE#CLOCK_RATE} -- the file records a per-branch clock RATE (BEAST {@code beast:rate}) and
 * divergence is DERIVED as rate x the branch's time span. This is an inference from the model that produced the tree,
 * not a measurement in the file, and the UI says so rather than presenting the two as the same kind of number.
 * Only when EVERY branch states a rate (finite, not negative): a branch without one has no divergence to draw, and
 * drawing it at 0 would state "no change along this branch", which the file never said. Joint with
 * Archaeopteryx.js ({@code everyBranchClockRate}).</li>
 * </ul>
 * TIME is computed the same way for both: the difference between a node's date and its parent's. That is what makes
 * one implementation serve both formats -- the dates are already there, whether they came from an Auspice
 * {@code num_date} or a BEAST {@code height}.
 * <p>
 * <b>The joint rule</b> (with Archaeopteryx.js; Christian, 2026-09-29): <i>Time | Div is offered only when both
 * layouts can state every branch: each non-root node and its parent state what the time layout reads, and each
 * non-root node states what the divergence layout reads. A value that is stated is stated, whether zero or negative;
 * only an absent one is missing. Otherwise the switch is not offered and the tree stays in the layout it arrived
 * in.</i> A branch drawn at 0 for want of a number says "nothing happened here", which no file said.
 * <p>
 * <b>A span may be negative.</b> Real files date a child before its parent (a summary tree's medians: 35 of the
 * 1372 branches of one influenza tree). TIME keeps the sign, in the direction the tree's dates run, so the lengths
 * from the root to a node add up to that node's own date; DIVERGENCE states 0, a negative amount of change meaning
 * nothing. This is about the VALUE -- what is shown as a branch length and written to a file. What is DRAWN for a
 * negative length is the painter's business, and the painter draws every negative length of every tree at 0
 * ({@code TreePanel.calculateBranchLengthToParent}).
 */
final class BranchLengthLayout {

    /** Which quantity the branch lengths currently express. */
    enum MODE {
        TIME( "Time" ),
        DIVERGENCE( "Divergence" );

        private final String _label;

        MODE( final String label ) {
            _label = label;
        }

        String label() {
            return _label;
        }

        @Override
        public String toString() {
            return _label;
        }
    }

    /** Where a tree's divergence numbers come from -- and whether they are recorded or inferred. */
    enum DIVERGENCE_SOURCE {
        STORED, CLOCK_RATE, NONE
    }

    /** Auspice / Nextstrain: cumulative divergence from the root, per node. */
    final static String DIV_PROPERTY_REF  = "nextstrain:div";
    /** BEAST: the clock rate on the branch leading to this node (substitutions/site per unit time). */
    final static String RATE_PROPERTY_REF = "beast:rate";

    private BranchLengthLayout() {
    }

    private static Double numericProperty( final PhylogenyNode node, final String ref ) {
        if ( ( node.getNodeData() == null ) || ( node.getNodeData().getProperties() == null ) ) {
            return null;
        }
        for( final Property p : node.getNodeData().getProperties().getProperties() ) {
            if ( ref.equals( p.getRef() ) ) {
                try {
                    final double d = Double.parseDouble( p.getValue().trim() );
                    return Double.isFinite( d ) ? Double.valueOf( d ) : null;
                }
                catch ( final Exception e ) {
                    return null;
                }
            }
        }
        return null;
    }

    /** The node's date VALUE, or null when it states none. Asks whether the value is THERE, never whether it is
     *  zero: 0 is the height of every contemporaneous BEAST tip. */
    private static Double dateValue( final PhylogenyNode node ) {
        if ( ( node.getNodeData() == null ) || ( node.getNodeData().getDate() == null )
                || ( node.getNodeData().getDate().getValue() == null ) ) {
            return null;
        }
        return Double.valueOf( node.getNodeData().getDate().getValue().doubleValue() );
    }

    /**
     * Which way the tree's dates run: true when they INCREASE toward the tips (calendar dates), false when they
     * decrease (ages, heights -- largest at the root). Measured on the tree, not read off a unit: the majority of the
     * parent-child pairs whose two dates differ decides (as Archaeopteryx.js measures it). No such pair, or as many
     * one way as the other: ages.
     */
    static boolean datesIncreaseTowardTips( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return false;
        }
        int up = 0;
        int down = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            final Double d = dateValue( n );
            final Double dp = dateValue( n.getParent() );
            if ( ( d == null ) || ( dp == null ) ) {
                continue;
            }
            if ( d.doubleValue() > dp.doubleValue() ) {
                ++up;
            }
            else if ( d.doubleValue() < dp.doubleValue() ) {
                ++down;
            }
        }
        return up > down;
    }

    /**
     * The time span of the branch leading to {@code node}, from the dates, or null when either end is undated.
     * SIGNED, in the direction the tree's dates run: negative when the node is dated BEFORE its parent. Summed from
     * the root, signed spans add up to each node's own date; their absolute values do not.
     */
    private static Double timeSpan( final PhylogenyNode node, final boolean increase_toward_tips ) {
        if ( node.isRoot() || ( node.getParent() == null ) ) {
            return null;
        }
        final Double d = dateValue( node );
        final Double dp = dateValue( node.getParent() );
        if ( ( d == null ) || ( dp == null ) ) {
            return null;
        }
        return Double.valueOf( increase_toward_tips ? ( d.doubleValue() - dp.doubleValue() )
                : ( dp.doubleValue() - d.doubleValue() ) );
    }

    /**
     * Whether the TIME layout can state every branch: every non-root node states a date, and so does its parent --
     * that is, every node of the tree. (It was a strict majority until the joint rule: a tree offered the switch on
     * a majority has branches the layout cannot state, and drew them at 0.)
     */
    static boolean isTimeDerivable( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return false;
        }
        int branches = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            ++branches;
            if ( ( dateValue( n ) == null ) || ( dateValue( n.getParent() ) == null ) ) {
                return false;
            }
        }
        return branches > 0;
    }

    /**
     * Where divergence comes from for this tree, when the divergence layout can state EVERY branch with it.
     * <ul>
     * <li>RECORDED ({@code nextstrain:div}): every node states one, the root included -- a branch is the difference
     * of two. Stated on some nodes only, the tree has no divergence layout: a recorded source wins over a derivable
     * one, so there is no falling back on rates.</li>
     * <li>CLOCK RATE: no node states a recorded divergence, and every non-root node states a rate that is finite and
     * not negative (the root has no branch, so its rate is never asked for). Its other factor, the time span, is
     * what {@link #isTimeDerivable} asks for.</li>
     * </ul>
     */
    static DIVERGENCE_SOURCE divergenceSource( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return DIVERGENCE_SOURCE.NONE;
        }
        boolean any_recorded = false;
        boolean every_node_recorded = true;
        boolean every_branch_rated = true;
        int branches = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( numericProperty( n, DIV_PROPERTY_REF ) != null ) {
                any_recorded = true;
            }
            else {
                every_node_recorded = false;
            }
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            ++branches;
            if ( every_branch_rated ) {
                final Double rate = numericProperty( n, RATE_PROPERTY_REF );
                every_branch_rated = ( rate != null ) && ( rate.doubleValue() >= 0 );
            }
        }
        if ( branches < 1 ) {
            return DIVERGENCE_SOURCE.NONE;
        }
        if ( any_recorded ) {
            return every_node_recorded ? DIVERGENCE_SOURCE.STORED : DIVERGENCE_SOURCE.NONE;
        }
        return every_branch_rated ? DIVERGENCE_SOURCE.CLOCK_RATE : DIVERGENCE_SOURCE.NONE;
    }

    /** Whether BOTH layouts can state every branch, and so the tree is offered the toggle. */
    static boolean isApplicable( final Phylogeny phy ) {
        return isTimeDerivable( phy ) && ( divergenceSource( phy ) != DIVERGENCE_SOURCE.NONE );
    }

    /**
     * Branch lengths = the time each branch spans, from the node dates, SIGNED (see {@link #timeSpan}). Never
     * refused: it is the way back from divergence. On a tree offered the switch every branch has a span; a branch
     * that has LOST a date since (the node editor, while divergence was on screen) has nothing to be laid out by and
     * goes to 0, because the alternative is the divergence length it still holds, in the wrong unit.
     */
    static void applyTime( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return;
        }
        final boolean up = datesIncreaseTowardTips( phy );
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            final Double t = timeSpan( n, up );
            n.setDistanceToParent( ( t == null ) ? 0.0 : t.doubleValue() );
        }
    }

    /** {@link #applyTime(Phylogeny)} for the branches BELOW {@code clade_root} only: its own branch, and every
     *  branch outside the clade, is left as it is. For a clade a subtree view rewrote, inside a tree that has no
     *  second layout of its own and so must not be rewritten as a whole. {@code increase_toward_tips} is the way
     *  the dates of the tree AROUND the clade run ({@link #datesIncreaseTowardTips}). */
    static void applyTimeBelow( final PhylogenyNode clade_root, final boolean increase_toward_tips ) {
        if ( clade_root == null ) {
            return;
        }
        for( final PhylogenyNode n : clade_root.getDescendants() ) {
            final Double t = timeSpan( n, increase_toward_tips );
            n.setDistanceToParent( ( t == null ) ? 0.0 : t.doubleValue() );
            applyTimeBelow( n, increase_toward_tips );
        }
    }

    /**
     * Branch lengths = genetic change along each branch: the successive difference of a recorded cumulative
     * divergence, or rate x time span when a clock rate is what the file states. A tree that is NOT offered the
     * switch is left exactly as it was -- no divergence source, or a branch a layout cannot state.
     * <p>
     * Nothing is absent here (that is what being offered the switch means), so there is no null to check for. What is
     * left to decide is a NEGATIVE amount: a recorded divergence that falls along a branch, or a span that runs
     * backwards. Change cannot be negative, so the length is 0 -- and max() also turns the negative zero of a rate
     * written "-0.0" into a plain one, so no length is ever written "-0.0".
     */
    static void applyDivergence( final Phylogeny phy ) {
        if ( !isApplicable( phy ) ) {
            return;
        }
        final DIVERGENCE_SOURCE source = divergenceSource( phy );
        final boolean up = datesIncreaseTowardTips( phy );
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            if ( source == DIVERGENCE_SOURCE.STORED ) {
                n.setDistanceToParent( Math.max( 0.0, numericProperty( n, DIV_PROPERTY_REF ).doubleValue()
                        - numericProperty( n.getParent(), DIV_PROPERTY_REF ).doubleValue() ) );
            }
            else {
                n.setDistanceToParent( Math.max( 0.0, numericProperty( n, RATE_PROPERTY_REF ).doubleValue()
                        * timeSpan( n, up ).doubleValue() ) );
            }
        }
    }

    /**
     * What to call the mode in the UI. A DERIVED divergence says so: on a BEAST tree the number is rate x time, an
     * inference from the clock model rather than something the file recorded, and a reader deserves to know which
     * they are looking at before they quote it.
     */
    static String label( final MODE mode, final DIVERGENCE_SOURCE source ) {
        if ( mode == MODE.DIVERGENCE ) {
            return ( source == DIVERGENCE_SOURCE.CLOCK_RATE ) ? "Divergence (from clock rate)" : "Divergence";
        }
        return MODE.TIME.label();
    }

    /** The distance unit to stamp on the tree for a mode, so the scale bar and the tree-properties window agree. */
    static String distanceUnit( final MODE mode, final String date_unit ) {
        if ( mode == MODE.DIVERGENCE ) {
            return "subs/site";
        }
        return ( ( date_unit == null ) || date_unit.trim().isEmpty() ) ? "time" : date_unit.trim();
    }
}
