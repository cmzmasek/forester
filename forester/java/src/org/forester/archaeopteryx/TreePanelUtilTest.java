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

package org.forester.archaeopteryx;

import java.awt.Color;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Unit tests for {@link TreePanelUtil}. Lives in the {@code org.forester.archaeopteryx} package
 * because the methods under test are package-private. Run standalone via {@link #main(String[])},
 * or as part of the suite via {@link #test()} (called from {@code org.forester.test.Test}).
 */
public final class TreePanelUtilTest {

    public static void main( final String[] args ) {
        System.out.println( "TreePanelUtil: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        return testYDistanceToAvoidLabelOverlap() && testSupportSymbolMath() && testDetectConfidenceScaleMax()
                && testCapEntries() && testTaxonomyLabel();
    }

    /**
     * taxonomyLabel picks scientific name, else common name, else taxonomy code -- so a rank
     * colorization legend gets a row even for taxa identified only by a common name or code (not just
     * a scientific name).
     */
    private static boolean testTaxonomyLabel() {
        final Taxonomy sci_and_common = new Taxonomy();
        sci_and_common.setScientificName( "Diptera" );
        sci_and_common.setCommonName( "flies" );
        if ( !"Diptera".equals( TreePanelUtil.taxonomyLabel( sci_and_common ) ) ) {
            return fail( "scientific name must win over common name" );
        }
        final Taxonomy common_only = new Taxonomy();
        common_only.setCommonName( "flies" );
        if ( !"flies".equals( TreePanelUtil.taxonomyLabel( common_only ) ) ) {
            return fail( "common name must be used when no scientific name (the legend-empty bug fix)" );
        }
        if ( !"".equals( TreePanelUtil.taxonomyLabel( new Taxonomy() ) )
                || !"".equals( TreePanelUtil.taxonomyLabel( null ) ) ) {
            return fail( "an empty or null taxonomy must yield an empty label" );
        }
        return true;
    }

    /** capEntries bounds a legend to its first N entries (iteration order), preserving keys+colors. */
    private static boolean testCapEntries() {
        if ( !TreePanelUtil.capEntries( null, 5 ).isEmpty() ) {
            return fail( "null input must yield an empty cap" );
        }
        final Map<String, Color> in = new LinkedHashMap<String, Color>();
        in.put( "a", Color.RED );
        in.put( "b", Color.GREEN );
        in.put( "c", Color.BLUE );
        // cap >= size keeps everything, in order
        final Map<String, Color> all = TreePanelUtil.capEntries( in, 5 );
        if ( !new ArrayList<String>( all.keySet() ).equals( new ArrayList<String>( in.keySet() ) ) ) {
            return fail( "cap >= size must keep all entries in iteration order; got " + all.keySet() );
        }
        if ( all.get( "a" ) != Color.RED ) {
            return fail( "cap must preserve the colors" );
        }
        // cap < size keeps exactly the first N entries
        final Map<String, Color> capped = TreePanelUtil.capEntries( in, 2 );
        if ( ( capped.size() != 2 ) || !capped.containsKey( "a" ) || !capped.containsKey( "b" )
                || capped.containsKey( "c" ) ) {
            return fail( "cap < size must keep exactly the first N entries; got " + capped.keySet() );
        }
        // cap of 0 -> empty (the legend's "+N more" footer then accounts for all of them)
        if ( !TreePanelUtil.capEntries( in, 0 ).isEmpty() ) {
            return fail( "cap of 0 must yield an empty map" );
        }
        return true;
    }

    /**
     * The node-symbol support math (see TreePanel.paintNodeSupportSymbol): scale detection picks the
     * absolute ceiling (1 or 100) so symbols mean the same across trees; the support fraction clamps
     * to 0..1; size interpolates min..max by that fraction; the threshold test compares the fraction
     * to the cutoff.
     */
    private static boolean testSupportSymbolMath() {
        // scale ceiling: anything above 1 implies the 0..100 family, otherwise 0..1
        if ( TreePanelUtil.confidenceScaleMaxFor( 0.0 ) != 1.0 ) {
            return fail( "empty/zero support implies the 0..1 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 0.95 ) != 1.0 ) {
            return fail( "max 0.95 (posterior probability) implies the 0..1 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 1.0 ) != 1.0 ) {
            return fail( "max exactly 1.0 stays on the 0..1 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 1.01 ) != 100.0 ) {
            return fail( "any value above 1 implies the 0..100 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 70.0 ) != 100.0 ) {
            return fail( "max 70 (bootstrap) implies the 0..100 scale" );
        }
        // support fraction: clamped to 0..1, scale-relative
        if ( TreePanelUtil.supportFraction( 50.0, 100.0 ) != 0.5 ) {
            return fail( "50 on a 0..100 scale is fraction 0.5" );
        }
        if ( TreePanelUtil.supportFraction( 0.8, 1.0 ) != 0.8 ) {
            return fail( "0.8 on a 0..1 scale is fraction 0.8" );
        }
        if ( TreePanelUtil.supportFraction( 150.0, 100.0 ) != 1.0 ) {
            return fail( "fraction must clamp to 1.0 above the scale ceiling" );
        }
        if ( TreePanelUtil.supportFraction( -1.0, 100.0 ) != 0.0 ) {
            return fail( "negative support clamps to fraction 0.0" );
        }
        if ( TreePanelUtil.supportFraction( 50.0, 0.0 ) != 0.0 ) {
            return fail( "a non-positive scale must not divide; fraction 0.0" );
        }
        // size interpolation: min at 0, max at full support, monotonic
        final float min = 2.0f;
        final float max = 8.0f;
        if ( TreePanelUtil.supportSymbolSize( 0.0, 100.0, min, max ) != min ) {
            return fail( "zero support gives the minimum symbol size" );
        }
        if ( TreePanelUtil.supportSymbolSize( 100.0, 100.0, min, max ) != max ) {
            return fail( "full support gives the maximum symbol size" );
        }
        if ( TreePanelUtil.supportSymbolSize( 50.0, 100.0, min, max ) != 5.0f ) {
            return fail( "half support gives the midpoint symbol size" );
        }
        if ( TreePanelUtil.supportSymbolSize( 30.0, 100.0, min, max ) >= TreePanelUtil
                .supportSymbolSize( 60.0, 100.0, min, max ) ) {
            return fail( "symbol size must grow with support" );
        }
        // threshold test (cutoff is a fraction of the scale)
        if ( !TreePanelUtil.isSupportAtOrAboveThreshold( 95.0, 100.0, 0.95 ) ) {
            return fail( "95/100 must meet the 0.95 cutoff" );
        }
        if ( TreePanelUtil.isSupportAtOrAboveThreshold( 94.0, 100.0, 0.95 ) ) {
            return fail( "94/100 must fall below the 0.95 cutoff" );
        }
        if ( !TreePanelUtil.isSupportAtOrAboveThreshold( 0.96, 1.0, 0.95 ) ) {
            return fail( "0.96 on a 0..1 scale must meet the 0.95 cutoff (scale-independent)" );
        }
        return true;
    }

    /** detectConfidenceScaleMax scans only internal-node confidences and infers the absolute scale. */
    private static boolean testDetectConfidenceScaleMax() {
        // a bootstrap tree (support 90 at the internal node) -> 0..100 scale
        if ( TreePanelUtil.detectConfidenceScaleMax( treeWithInternalConfidence( 90.0, "bootstrap" ) ) != 100.0 ) {
            return fail( "a tree with bootstrap support 90 must be detected as the 0..100 scale" );
        }
        // a Bayesian tree (posterior probability 0.99) -> 0..1 scale
        if ( TreePanelUtil.detectConfidenceScaleMax( treeWithInternalConfidence( 0.99, "posterior" ) ) != 1.0 ) {
            return fail( "a tree with posterior probability 0.99 must be detected as the 0..1 scale" );
        }
        // no confidences anywhere -> defaults to the 0..1 scale (harmless; nothing is drawn)
        final Phylogeny bare = new Phylogeny();
        final PhylogenyNode r = new PhylogenyNode();
        r.addAsChild( new PhylogenyNode() );
        r.addAsChild( new PhylogenyNode() );
        bare.setRoot( r );
        bare.externalNodesHaveChanged();
        if ( TreePanelUtil.detectConfidenceScaleMax( bare ) != 1.0 ) {
            return fail( "a tree without confidences defaults to the 0..1 scale" );
        }
        return true;
    }

    /** (internal:conf, leaf, leaf) under a root -- the confidence sits on an internal branch. */
    private static Phylogeny treeWithInternalConfidence( final double value, final String type ) {
        final Phylogeny p = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode internal = new PhylogenyNode();
        internal.getBranchData().addConfidence( new Confidence( value, type ) );
        internal.addAsChild( new PhylogenyNode() );
        internal.addAsChild( new PhylogenyNode() );
        root.addAsChild( internal );
        root.addAsChild( new PhylogenyNode() );
        p.setRoot( root );
        p.externalNodesHaveChanged();
        return p;
    }

    /**
     * The y-distance returned for a given label height must (a) space adjacent leaf rows
     * (2 * y-distance apart) at least one label-height apart so labels do not overlap, and
     * (b) drive the dynamic-hiding factor -- the same formula TreePanel.calcDynamicHidingFactor
     * uses -- down to <= 1 so the "Dyna Hide" indicator clears.
     */
    private static boolean testYDistanceToAvoidLabelOverlap() {
        final int[] heights = { 2, 8, 10, 11, 12, 14, 16, 20, 27, 40 };
        float previous = -1.0f;
        for( final int h : heights ) {
            final float y_dist = TreePanelUtil.yDistanceToAvoidLabelOverlap( h );
            if ( y_dist <= 0.0f ) {
                return fail( "y-distance must be positive for height " + h + " (got " + y_dist + ")" );
            }
            // (a) leaf rows are 2 * y-distance apart; that must be >= the label height
            if ( ( 2.0f * y_dist ) < h ) {
                return fail( "labels would overlap at height " + h + ": spacing " + ( 2.0f * y_dist ) + " < " + h );
            }
            // (b) same formula as TreePanel.calcDynamicHidingFactor: round( h / (1.5 * y-distance) )
            final int hiding_factor = (int) ( 0.5 + ( h / ( 1.5 * y_dist ) ) );
            if ( hiding_factor > 1 ) {
                return fail( "dynamic-hiding factor should be <= 1 at height " + h + " (got " + hiding_factor + ")" );
            }
            // monotonic in the label height (taller labels never need less spacing)
            if ( y_dist < previous ) {
                return fail( "y-distance should grow with label height; broke at height " + h );
            }
            previous = y_dist;
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [TreePanelUtilTest] " + message );
        return false;
    }

    private TreePanelUtilTest() {
        // not instantiable
    }
}
