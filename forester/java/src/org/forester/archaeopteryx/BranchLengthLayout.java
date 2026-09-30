// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

package org.forester.archaeopteryx;

import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Laying a dated tree out by TIME or by DIVERGENCE, for the trees that carry both.
 * <p>
 * A time-calibrated tree usually knows two different things about every branch: how much TIME it spans, and how much
 * genetic change happened along it. Only one can be the x-axis at a time, so this switches between them. It is a pure
 * DISPLAY mode: both quantities stay retained, nothing is edited, and switching back is exact.
 * <p>
 * <b>TIME is what the tree shows when it is not showing divergence</b>: the branch lengths it arrived with (the
 * lengths its file states; for an Auspice build, which states none, the ones its reader made of its dates). They are
 * KEPT when the tree leaves time ({@link TimeLengths}) and given back when it returns -- the picture the tree opened
 * with, to the digit. It used to be recomputed from the node dates, which is a different tree wherever a file's
 * lengths and its dates disagree: one summary tree states a branch of 1.407 between two nodes dated 1.098 apart
 * (its lengths and its median heights are two summaries of one posterior). Joint with Archaeopteryx.js (Christian,
 * 2026-09-29).
 * <p>
 * DIVERGENCE comes in two shapes, and the difference between them is not cosmetic:
 * <ul>
 * <li>{@link DIVERGENCE_SOURCE#STORED} -- the file RECORDS divergence per node, as a cumulative value (Auspice /
 * Nextstrain {@code nextstrain:div}). Branch divergence is the successive difference.</li>
 * <li>{@link DIVERGENCE_SOURCE#CLOCK_RATE} -- the file records a per-branch clock RATE (BEAST {@code beast:rate}) and
 * divergence is DERIVED as rate x the branch's length in time -- for a branch that took in a deleted node, each
 * piece at its own node's rate ({@link TimeLengths#divergence}, joint with Archaeopteryx.js). This is an inference
 * from the model that produced the tree, not a measurement in the file, and the UI says so rather than presenting
 * the two as the same kind of number.</li>
 * </ul>
 * <b>A tree may ARRIVE showing divergence</b>: one saved while divergence was on screen carries its dates, its
 * rates or recorded divergence, and branch lengths that are its DIVERGENCE. Taken for time, those lengths would be
 * what Time gives back, and what a rate is multiplied by. So it is asked once, when the tree arrives
 * ({@link #arrivesShowingDivergence}); such a tree starts in divergence with no length kept, and its time is laid
 * out from its dates.
 * <p>
 * <b>The joint rule</b> (with Archaeopteryx.js; Christian, 2026-09-29): <i>Time | Div is offered only when both
 * layouts can state every branch: each non-root node and its parent state what the time layout reads, and each
 * non-root node states what the divergence layout reads. A value that is stated is stated, whether zero or negative;
 * only an absent one is missing. Otherwise the switch is not offered and the tree stays in the layout it arrived
 * in.</i> A branch drawn at 0 for want of a number says "nothing happened here", which no file said. Here: a date on
 * every node and a length on every branch; a recorded divergence on every node, or a rate on every branch. And
 * each picture must have some DEPTH, measured at the tips as Archaeopteryx.js measures it: a tip whose divergence
 * from the root is above 0, and a tip whose date is not the root's. A tree whose divergence is 0 along every branch
 * is not offered a picture of nothing.
 * <p>
 * <b>One tree, one layout.</b> The switch acts on the WHOLE tree a tab holds and is judged on the whole tree, also
 * while a subtree of it is on view (as in Archaeopteryx.js): a subtree view shares its nodes with the tree it is a
 * view of, so a switch that acted on the view alone left the whole tree in two units at once.
 * <p>
 * <b>A length may be negative.</b> Real files state a node before its parent (35 of the 1372 branches of one
 * influenza tree). TIME gives the length back as stated; DIVERGENCE states 0, a negative amount of change meaning
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

    /**
     * The lengths a tree's branches have IN TIME: what the time layout reads, and what Time gives back.
     * <p>
     * Kept by node id, which is what a node has in common with its copies: an undo snapshot is a copy of the tree, and
     * the root of a subtree view is a copy of the clade's top node -- each with the id of its original. Each length is
     * kept WITH THE ID OF THE PARENT it was measured to, because the tree can change while divergence is on screen:
     * Delete Node, Cut Subtree and the prunes remove a node and add its length to its child's. A branch whose parent
     * is no longer the one its length was kept to spans, in time, its own kept length and the kept lengths of the
     * nodes removed between them; a branch whose way up to its parent is not all kept has no kept length, and is
     * completed from its dates. Immutable once handed out, so an undo snapshot can hold one.
     */
    static final class TimeLengths {

        private static final class Kept {

            final double _length;
            final long   _parent;
            /** The clock rate the node stated for its branch when this was kept, or null: what a piece of a
             *  merged branch goes on being multiplied by after its node is gone. */
            final Double _rate;

            Kept( final double length, final long parent, final Double rate ) {
                _length = length;
                _parent = parent;
                _rate = rate;
            }
        }

        private final Map<Long, Kept> _of = new HashMap<Long, Kept>();

        /** The lengths ON SCREEN now, of trees that are showing time. Where two trees hold a node each with the same
         *  id, the tree given LAST is the one that counts. A length that is not stated is not kept. */
        static TimeLengths onScreen( final Phylogeny... trees ) {
            final TimeLengths t = new TimeLengths();
            for( final Phylogeny phy : trees ) {
                t.keep( phy );
            }
            return t;
        }

        /** These lengths, with {@code other}'s for every branch of {@code phy} (not its root's): an undo snapshot's lengths for the
         *  tree it restores, over what the tab has for the rest of its trees. */
        TimeLengths over( final TimeLengths other, final Phylogeny phy ) {
            final TimeLengths t = new TimeLengths();
            t._of.putAll( _of );
            if ( ( other != null ) && ( phy != null ) && !phy.isEmpty() ) {
                for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if ( n.isRoot() || ( n.getParent() == null ) ) {
                        continue; // no branch IN phy: the root of a view stands for a node of the tree beneath it
                    }
                    final Long id = Long.valueOf( n.getId() );
                    if ( other._of.containsKey( id ) ) {
                        t._of.put( id, other._of.get( id ) );
                    }
                }
            }
            return t;
        }

        private void keep( final Phylogeny phy ) {
            if ( ( phy == null ) || phy.isEmpty() ) {
                return;
            }
            for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( n.isRoot() || ( n.getParent() == null ) ) {
                    continue;
                }
                final double d = n.getDistanceToParent();
                if ( ( d != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) && !Double.isNaN( d ) && !Double.isInfinite( d ) ) {
                    _of.put( Long.valueOf( n.getId() ), new Kept( d, n.getParent().getId(), rate( n ) ) );
                }
            }
        }

        /** The PIECES of the branch leading to {@code node}, the node's own first: its own kept entry when its
         *  parent is the one it was kept to, else the kept entries along the way up to its parent (the nodes removed
         *  since). Null when the way up is not all kept. */
        private java.util.List<Kept> pieces( final PhylogenyNode node ) {
            if ( node.getParent() == null ) {
                return null;
            }
            final long parent = node.getParent().getId();
            final java.util.List<Kept> pieces = new java.util.ArrayList<Kept>( 2 );
            long id = node.getId();
            for( int step = 0; step <= _of.size(); ++step ) {
                final Kept k = _of.get( Long.valueOf( id ) );
                if ( k == null ) {
                    return null;
                }
                pieces.add( k );
                if ( k._parent == parent ) {
                    return pieces;
                }
                id = k._parent;
            }
            return null;
        }

        /** The length in time of the branch leading to {@code node}, or null when none is kept for it: the kept
         *  lengths of its pieces, summed, with their signs. */
        Double of( final PhylogenyNode node ) {
            final java.util.List<Kept> pieces = pieces( node );
            if ( pieces == null ) {
                return null;
            }
            double sum = 0;
            for( final Kept k : pieces ) {
                sum += k._length;
            }
            return Double.valueOf( sum );
        }

        /**
         * Rate x time along the branch leading to {@code node}, PIECE BY PIECE: its own piece at the rate it states
         * now ({@code rate}), a piece a removed node left at the rate that node stated (at {@code rate} when it
         * stated none) -- the change along a branch does not become another change because a node on it was
         * deleted (as Archaeopteryx.js keeps it, per segment). Null when no length is kept for the branch.
         */
        Double divergence( final PhylogenyNode node, final double rate ) {
            final java.util.List<Kept> pieces = pieces( node );
            if ( pieces == null ) {
                return null;
            }
            double sum = rate * pieces.get( 0 )._length;
            for( int i = 1; i < pieces.size(); ++i ) {
                final Kept k = pieces.get( i );
                sum += ( ( k._rate != null ) ? k._rate.doubleValue() : rate ) * k._length;
            }
            return Double.valueOf( sum );
        }

        /**
         * These lengths brought up to what {@code trees} (showing TIME) show now: every stated length on screen is
         * kept, except where a branch on screen is still the SUM of the pieces it was made of by a delete -- then
         * the pieces stay, with the rates of the nodes that are gone. What a tab in time remembers before each change
         * to its tree, so that a node deleted while time is on screen keeps its rate too.
         */
        TimeLengths refreshedBy( final Phylogeny... trees ) {
            final TimeLengths t = new TimeLengths();
            t._of.putAll( _of );
            for( final Phylogeny phy : trees ) {
                if ( ( phy == null ) || phy.isEmpty() ) {
                    continue;
                }
                for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if ( n.isRoot() || ( n.getParent() == null ) ) {
                        continue;
                    }
                    final double d = n.getDistanceToParent();
                    if ( ( d == PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) || Double.isNaN( d ) || Double.isInfinite( d ) ) {
                        t._of.remove( Long.valueOf( n.getId() ) ); // a branch that states none on screen has none
                        continue;
                    }
                    final java.util.List<Kept> pieces = pieces( n );
                    if ( ( pieces != null ) && ( pieces.size() > 1 ) ) {
                        double sum = 0;
                        for( final Kept k : pieces ) {
                            sum += k._length;
                        }
                        if ( Math.abs( sum - d ) <= ( 1e-9 * Math.max( 1.0, Math.abs( d ) ) ) ) {
                            continue; // still the sum of its pieces: they stay
                        }
                    }
                    t._of.put( Long.valueOf( n.getId() ), new Kept( d, n.getParent().getId(), rate( n ) ) );
                }
            }
            return t;
        }

        /**
         * These lengths, and for every branch of {@code trees} that has none among them the GAP BETWEEN ITS DATES,
         * signed the way the dates run ({@code up}: they increase toward the tips). For a tree that is showing
         * divergence: a branch no length was kept for -- every branch of a tree that ARRIVED showing divergence, a
         * node added since -- has its dates to say what it spans. A branch with neither stays without a length.
         */
        TimeLengths completedFromDates( final boolean up, final Phylogeny... trees ) {
            final TimeLengths t = new TimeLengths();
            t._of.putAll( _of );
            for( final Phylogeny phy : trees ) {
                if ( ( phy == null ) || phy.isEmpty() ) {
                    continue;
                }
                for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if ( ( n.getParent() != null ) && ( of( n ) == null ) ) {
                        final Double gap = dateGap( n, up );
                        if ( gap != null ) {
                            t._of.put( Long.valueOf( n.getId() ), new Kept( gap.doubleValue(), n.getParent().getId(), rate( n ) ) );
                        }
                    }
                }
            }
            return t;
        }
    }

    /** Auspice / Nextstrain: cumulative divergence from the root, per node. */
    final static String          DIV_PROPERTY_REF  = "nextstrain:div";
    /** BEAST: the clock rate on the branch leading to this node (substitutions/site per unit time). */
    final static String          RATE_PROPERTY_REF = "beast:rate";
    /** JOINT with Archaeopteryx.js: a rate AND a recorded divergence are written as a plain decimal number, with or
     *  without an exponent (Christian, 2026-09-29, both). {@code Double.parseDouble} alone also reads "0.005d", "1f"
     *  and "0x1p-8", which no other reader of the file takes for that number. */
    private final static Pattern PLAIN_DECIMAL     = Pattern.compile( "[+-]?(\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+)?" );

    private BranchLengthLayout() {
    }

    private static String property( final PhylogenyNode node, final String ref ) {
        if ( ( node.getNodeData() == null ) || ( node.getNodeData().getProperties() == null ) ) {
            return null;
        }
        for( final Property p : node.getNodeData().getProperties().getProperties() ) {
            if ( ref.equals( p.getRef() ) ) {
                return ( p.getValue() == null ) ? null : p.getValue().trim();
            }
        }
        return null;
    }

    /** The number {@code text} states when it is a plain decimal and finite, else null. */
    private static Double plainDecimal( final String text ) {
        return ( ( text == null ) || !PLAIN_DECIMAL.matcher( text ).matches() ) ? null : finite( text );
    }

    private static Double finite( final String number ) {
        if ( number == null ) {
            return null;
        }
        try {
            final double d = Double.parseDouble( number );
            return Double.isFinite( d ) ? Double.valueOf( d ) : null;
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    /** The node's recorded cumulative divergence, or null when it records none that can be read as one: a plain
     *  decimal, finite. A NEGATIVE one is a value (Christian, 2026-09-29, joint): it is what the file records, and a
     *  branch along which it falls is drawn at 0 ({@link #divergence}). */
    private static Double recordedDivergence( final PhylogenyNode node ) {
        return plainDecimal( property( node, DIV_PROPERTY_REF ) );
    }

    /** The clock rate the node states for its branch, or null when it states none that can be read as one: a plain
     *  decimal, finite. Whether it may be NEGATIVE is for the caller to say. */
    private static Double rate( final PhylogenyNode node ) {
        return plainDecimal( property( node, RATE_PROPERTY_REF ) );
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
     * The gap between the dates of a node and its parent, or null when either is undated. SIGNED, in the direction
     * the tree's dates run: negative when the node is dated BEFORE its parent. What a branch spans in time where
     * no length was kept for it ({@link TimeLengths#completedFromDates}).
     */
    private static Double dateGap( final PhylogenyNode node, final boolean increase_toward_tips ) {
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
     * Whether the TIME layout can state every branch: every node of the tree states a date, the root included, and
     * every branch has a length in time.
     *
     * @param time the lengths the tree has in time: on screen, or kept
     */
    static boolean isTimeDerivable( final Phylogeny phy, final TimeLengths time ) {
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
            if ( ( dateValue( n ) == null ) || ( dateValue( n.getParent() ) == null ) || ( time.of( n ) == null ) ) {
                return false;
            }
        }
        return branches > 0;
    }

    /** {@link #isTimeDerivable(Phylogeny, TimeLengths)} of a tree that is showing time. */
    static boolean isTimeDerivable( final Phylogeny phy ) {
        return isTimeDerivable( phy, TimeLengths.onScreen( phy ) );
    }

    /**
     * Where divergence comes from for this tree, when the divergence layout can state EVERY branch with it.
     * <ul>
     * <li>RECORDED ({@code nextstrain:div}): every node states one, the root included -- a branch is the difference
     * of two. Stated on some nodes only, the tree has no divergence layout: a recorded source wins over a derivable
     * one, so there is no falling back on rates.</li>
     * <li>CLOCK RATE: no node states a recorded divergence, and every non-root node states a rate that is a plain
     * decimal, finite and not negative (the root has no branch, so its rate is never asked for). Its other factor,
     * the length in time, is what {@link #isTimeDerivable} asks for.</li>
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
            if ( recordedDivergence( n ) != null ) {
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
                final Double rate = rate( n );
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

    /**
     * The divergence along the branch leading to {@code n}, of a tree both layouts can state every branch of --
     * nothing is absent here, so there is no null to check for. What is left to decide is a NEGATIVE amount: a
     * recorded divergence that falls along a branch, or a rate times a length that is negative. Change cannot be
     * negative, so it is 0 -- and max() also turns the negative zero of a rate written "-0.0" into a plain one, so
     * no length is ever written "-0.0".
     */
    private static double divergence( final PhylogenyNode n, final DIVERGENCE_SOURCE source, final TimeLengths time ) {
        if ( source == DIVERGENCE_SOURCE.STORED ) {
            return Math.max( 0.0, recordedDivergence( n ).doubleValue() - recordedDivergence( n.getParent() ).doubleValue() );
        }
        return Math.max( 0.0, time.divergence( n, rate( n ).doubleValue() ).doubleValue() );
    }

    /**
     * Whether BOTH layouts can state every branch, and both pictures have some depth, and so the tree is offered
     * the toggle. JOINT, the depth as Archaeopteryx.js measures it, at the TIPS: the divergence from the root to
     * some tip is above 0 (divergence is never negative, so that is any branch above 0), and some tip is dated
     * differently from the root (NOT "some branch spans time": branches of +1 and -1 leave a tip on the root's
     * date).
     *
     * @param time the lengths the tree has in time: on screen, or kept
     */
    static boolean isApplicable( final Phylogeny phy, final TimeLengths time ) {
        final DIVERGENCE_SOURCE source = everyBranchSource( phy, time );
        if ( source == DIVERGENCE_SOURCE.NONE ) {
            return false;
        }
        final double root_date = dateValue( phy.getRoot() ).doubleValue();
        boolean divergence_has_depth = false;
        boolean time_has_depth = false;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            if ( !divergence_has_depth && ( divergence( n, source, time ) > 0 ) ) {
                divergence_has_depth = true;
            }
            if ( n.isExternal() && ( dateValue( n ).doubleValue() != root_date ) ) {
                time_has_depth = true;
            }
            if ( divergence_has_depth && time_has_depth ) {
                return true;
            }
        }
        return false;
    }

    /** The divergence source when BOTH layouts can state every branch of {@code phy} (no depth asked), else NONE. */
    private static DIVERGENCE_SOURCE everyBranchSource( final Phylogeny phy, final TimeLengths time ) {
        return isTimeDerivable( phy, time ) ? divergenceSource( phy ) : DIVERGENCE_SOURCE.NONE;
    }

    /** {@link #isApplicable(Phylogeny, TimeLengths)} of a tree that is showing time. */
    static boolean isApplicable( final Phylogeny phy ) {
        return isApplicable( phy, TimeLengths.onScreen( phy ) );
    }

    /**
     * Branch lengths = the lengths the tree has in time. Never refused: it is the way back from divergence. A
     * branch that has NO length in time (none kept, and no dates to complete it from) goes to 0, because the
     * alternative is the divergence length it still holds, in the wrong unit.
     */
    static void applyTime( final Phylogeny phy, final TimeLengths time ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return;
        }
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            final Double t = time.of( n );
            n.setDistanceToParent( ( t == null ) ? 0.0 : t.doubleValue() );
        }
    }

    /**
     * Whether the tree, as it ARRIVES, is showing divergence -- its branch lengths are its divergence, not its time.
     * Asked only of a tree that states everything both layouts read; any other tree shows what its file states and
     * is not asked.
     * <ul>
     * <li>A tree that RECORDS its divergence is showing time when its lengths are the gaps between its dates, signed
     * the way its dates run: 19 branches in 20, each to a millionth of the larger of the two (JOINT, the test of
     * Archaeopteryx.js). Otherwise it is showing divergence.</li>
     * <li>A CLOCK-RATE tree's lengths are not the gaps between its dates even in time (they are two summaries of one
     * posterior, and differ by a quarter on some branches), so that test cannot be asked of it. It is showing
     * divergence when its file SAYS so: the unit of its branch lengths is the one the divergence layout stamps.</li>
     * </ul>
     */
    static boolean arrivesShowingDivergence( final Phylogeny phy ) {
        if ( !isTimeDerivable( phy ) ) {
            return false;
        }
        final DIVERGENCE_SOURCE source = divergenceSource( phy );
        if ( source == DIVERGENCE_SOURCE.NONE ) {
            return false;
        }
        if ( source == DIVERGENCE_SOURCE.CLOCK_RATE ) {
            return ( phy.getDistanceUnit() != null )
                    && distanceUnit( MODE.DIVERGENCE, null ).equals( phy.getDistanceUnit().trim() );
        }
        final boolean up = datesIncreaseTowardTips( phy );
        int branches = 0;
        int are_the_gap = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            ++branches;
            final double length = n.getDistanceToParent();
            final double gap = dateGap( n, up ).doubleValue();
            if ( Math.abs( gap - length ) <= ( Math.max( Math.max( Math.abs( length ), Math.abs( gap ) ), 1e-9 ) * 1e-6 ) ) {
                ++are_the_gap;
            }
        }
        return ( are_the_gap * 20 ) < ( branches * 19 );
    }

    /**
     * Branch lengths = genetic change along each branch: the successive difference of a recorded cumulative
     * divergence, or rate x the branch's length in time when a clock rate is what the file states. A tree that is
     * NOT offered the switch is left exactly as it was.
     *
     * @param time the lengths the tree has in time. NOT read off the tree: two trees of a tab share their nodes,
     *            and the second to be laid out would take the divergence the first has just written for time
     */
    static void applyDivergence( final Phylogeny phy, final TimeLengths time ) {
        if ( isApplicable( phy, time ) ) {
            layOutDivergence( phy, time, divergenceSource( phy ) );
        }
    }

    /**
     * {@link #applyDivergence}, for the trees of a tab that was judged BEFORE its mode was set: the whole tree, the
     * clade of a subtree view, or a copy of one an undo put on display. It is asked only whether both layouts can
     * state its every branch, not whether its pictures have depth -- a clade of identical sequences has none, and is
     * a part of a tree in divergence all the same, not a tree in time stamped with the unit of divergence. A tree
     * that cannot state every branch is left as it was.
     */
    static void applyDivergenceToPart( final Phylogeny phy, final TimeLengths time ) {
        final DIVERGENCE_SOURCE source = everyBranchSource( phy, time );
        if ( source != DIVERGENCE_SOURCE.NONE ) {
            layOutDivergence( phy, time, source );
        }
    }

    private static void layOutDivergence( final Phylogeny phy, final TimeLengths time, final DIVERGENCE_SOURCE source ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isRoot() && ( n.getParent() != null ) ) {
                n.setDistanceToParent( divergence( n, source, time ) );
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
