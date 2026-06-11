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
        return testYDistanceToAvoidLabelOverlap();
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
