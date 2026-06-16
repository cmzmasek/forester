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
 * Unit tests for {@link TreeColorSet}: Archaeopteryx now has exactly two tree color schemes -- Dark
 * (index 0) and Light (index 1) -- matching the light/dark UI theme. There is no scheme cycling or
 * chooser. Headless; runs in the suite.
 */
public final class TreeColorSchemeTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreeColorScheme: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        // exactly two schemes, named Dark and Light
        if ( ( TreeColorSet.SCHEME_NAMES.length != 2 ) || !"Dark".equals( TreeColorSet.SCHEME_NAMES[ 0 ] )
                || !"Light".equals( TreeColorSet.SCHEME_NAMES[ 1 ] ) ) {
            return false;
        }
        final TreeColorSet tcs = TreeColorSet.createInstance();
        // every scheme row must hold exactly one color per COLOR_FIELDS entry, otherwise setColorSchema's
        // positional [0..N-1] mapping would silently shift colors or throw ArrayIndexOutOfBounds
        if ( tcs._color_schemes.length != TreeColorSet.SCHEME_NAMES.length ) {
            return false;
        }
        for ( final Color[] row : tcs._color_schemes ) {
            if ( row.length != TreeColorSet.COLOR_FIELDS.length ) {
                return false;
            }
        }
        // the default (scheme 0) is Dark: the modern dark-gray background, light branches
        if ( tcs.getCurrentColorScheme() != 0 ) {
            return false;
        }
        tcs.setColorSchema( 0 );
        if ( !new Color( 43, 43, 43 ).equals( tcs.getBackgroundColor() )
                || !"Dark".equals( tcs.getCurrentColorSchemeName() ) ) {
            return false;
        }
        if ( !Color.WHITE.equals( tcs.getBranchColor() ) ) {
            return false; // Dark scheme draws branches white
        }
        // scheme 1 is Light: white background, black branches/text
        tcs.setColorSchema( 1 );
        if ( !Color.WHITE.equals( tcs.getBackgroundColor() ) || !Color.BLACK.equals( tcs.getBranchColor() )
                || !Color.BLACK.equals( tcs.getSequenceColor() ) || !"Light".equals( tcs.getCurrentColorSchemeName() ) ) {
            return false;
        }
        return true;
    }

    private TreeColorSchemeTest() {
    }
}
