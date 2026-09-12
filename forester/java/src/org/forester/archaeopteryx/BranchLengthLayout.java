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
 * not a measurement in the file, and the UI says so rather than presenting the two as the same kind of number.</li>
 * </ul>
 * TIME is computed the same way for both: the difference between a node's date and its parent's. That is what makes
 * one implementation serve both formats -- the dates are already there, whether they came from an Auspice
 * {@code num_date} or a BEAST {@code height}.
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

    private static Double dateValue( final PhylogenyNode node ) {
        if ( ( node.getNodeData() == null ) || !node.getNodeData().isHasDate()
                || ( node.getNodeData().getDate().getValue() == null ) ) {
            return null;
        }
        return Double.valueOf( node.getNodeData().getDate().getValue().doubleValue() );
    }

    /** The time span of the branch leading to {@code node}, from the dates, or null when either end is undated. */
    private static Double timeLength( final PhylogenyNode node ) {
        if ( node.isRoot() || ( node.getParent() == null ) ) {
            return null;
        }
        final Double d = dateValue( node );
        final Double dp = dateValue( node.getParent() );
        return ( ( d == null ) || ( dp == null ) ) ? null : Double.valueOf( Math.abs( dp.doubleValue() - d.doubleValue() ) );
    }

    /**
     * Whether a TIME layout can be computed: a strict majority of the non-root nodes must have a date AND a dated
     * parent. A majority rather than all, because a real tree can carry a stray undated node; but not "any", or two
     * stray dates on an undated tree would offer a toggle that lays almost every branch out at zero.
     */
    static boolean isTimeDerivable( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return false;
        }
        int derivable = 0;
        int non_root = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            ++non_root;
            if ( timeLength( n ) != null ) {
                ++derivable;
            }
        }
        return ( non_root > 0 ) && ( ( derivable * 2 ) > non_root );
    }

    /** Where divergence would come from for this tree. A RECORDED source wins over a derivable one. */
    static DIVERGENCE_SOURCE divergenceSource( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return DIVERGENCE_SOURCE.NONE;
        }
        boolean any_rate = false;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( numericProperty( n, DIV_PROPERTY_REF ) != null ) {
                return DIVERGENCE_SOURCE.STORED;
            }
            if ( !any_rate && !n.isRoot() && ( numericProperty( n, RATE_PROPERTY_REF ) != null ) ) {
                any_rate = true;
            }
        }
        return any_rate ? DIVERGENCE_SOURCE.CLOCK_RATE : DIVERGENCE_SOURCE.NONE;
    }

    /** Whether the tree can be laid out BOTH ways, and so should be offered the toggle. */
    static boolean isApplicable( final Phylogeny phy ) {
        return isTimeDerivable( phy ) && ( divergenceSource( phy ) != DIVERGENCE_SOURCE.NONE );
    }

    /** Branch lengths = the time each branch spans, from the node dates. An underivable branch gets 0 rather than a
     *  stale cross-scale length left over from the other mode. */
    static void applyTime( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return;
        }
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            final Double t = timeLength( n );
            n.setDistanceToParent( ( t == null ) ? 0.0 : t.doubleValue() );
        }
    }

    /** Branch lengths = genetic change along each branch: the successive difference of a recorded cumulative
     *  divergence, or rate x time when only a clock rate is recorded. Same zero-rather-than-stale rule. */
    static void applyDivergence( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return;
        }
        final DIVERGENCE_SOURCE source = divergenceSource( phy );
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || ( n.getParent() == null ) ) {
                continue;
            }
            double len = 0.0;
            if ( source == DIVERGENCE_SOURCE.STORED ) {
                final Double d = numericProperty( n, DIV_PROPERTY_REF );
                final Double dp = numericProperty( n.getParent(), DIV_PROPERTY_REF );
                if ( ( d != null ) && ( dp != null ) ) {
                    len = Math.max( 0.0, d.doubleValue() - dp.doubleValue() ); // clamp a spurious negative to 0
                }
            }
            else if ( source == DIVERGENCE_SOURCE.CLOCK_RATE ) {
                final Double rate = numericProperty( n, RATE_PROPERTY_REF );
                final Double t = timeLength( n );
                if ( ( rate != null ) && ( t != null ) ) {
                    len = Math.max( 0.0, rate.doubleValue() * t.doubleValue() );
                }
            }
            n.setDistanceToParent( len );
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
