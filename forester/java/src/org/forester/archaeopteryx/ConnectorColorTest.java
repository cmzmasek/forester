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

/**
 * Unit tests for {@link TreePanel#connectorColor()}: the lined-up-data connector (the dashed guide
 * line drawn in ALIGNED_PHYLOGRAM mode between a branch tip and its right-aligned label) must use a
 * fixed, neutral guide color. The bug it guards against: {@code drawConnection} used to draw with
 * the {@link java.awt.Graphics2D}'s current color, which is the node's label color. For a
 * search-found node that color is the bright red/green highlight; on screen the sub-pixel dashed
 * stroke antialiased it away, but in a vector PDF/PNG export it rendered as a fully saturated
 * red/green "branch". The connector color must therefore be neutral and must NOT equal any of the
 * search-found highlight colors. In {@code org.forester.archaeopteryx} so it can reach the
 * package-private method; run via the suite ({@link #test()}).
 */
public final class ConnectorColorTest {

    public static void main( final String[] args ) {
        System.out.println( "ConnectorColor: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        final Color c = TreePanel.connectorColor();
        if ( c == null ) {
            return fail( "connector color is null" );
        }
        // neutral: a true gray (no hue) -- never a colored highlight
        if ( ( c.getRed() != c.getGreen() ) || ( c.getGreen() != c.getBlue() ) ) {
            return fail( "connector color is not a neutral gray: " + c );
        }
        // faint, but visible: a light guide, not black and not white
        if ( ( c.getRed() <= 0 ) || ( c.getRed() >= 255 ) ) {
            return fail( "connector color should be a light-but-visible gray: " + c );
        }
        // crucially, must not be any of the search-found highlight colors (red / green / yellow / blue)
        final TreeColorSet tcs = TreeColorSet.createInstance( new Configuration( null, false, false, false ) );
        if ( c.equals( tcs.getFoundColor0() ) || c.equals( tcs.getFoundColor1() )
                || c.equals( tcs.getFoundColor0and1() ) ) {
            return fail( "connector color must not be a search-found highlight color: " + c );
        }
        // and it must be stable / independent of any caller state: same value every call
        if ( !c.equals( TreePanel.connectorColor() ) ) {
            return fail( "connector color must be stable across calls" );
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [ConnectorColorTest] " + message );
        return false;
    }

    private ConnectorColorTest() {
    }
}
