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
        // NO RAW PRIMARY may sit in a scheme row. The six dead slots that held them -- the three MATCHING_NODES_*,
        // NODE_BOX, BINARY_DOMAIN_COMBINATIONS and ANNOTATION -- were removed once their features were; this keeps
        // (255,0,0), (0,255,0), (0,0,255), (255,255,0) and chartreuse from creeping back into a row. Greys and the
        // black/white ends are fine: it is SATURATION that is garish, not extremity.
        // What makes a colour garish is not extremity but sitting on the OUTER SKIN of the RGB cube: a channel
        // blown out to 255 while another is near zero. That catches (255,0,0), (0,255,0), (0,0,255), (255,255,0),
        // chartreuse (173,255,47) and royal blue (65,105,255) -- every value the removed slots held -- while
        // leaving the deliberate Okabe-Ito tones alone, since vermillion (213,94,0) and amber (230,159,0) top out
        // BELOW 255. Greys and pure white/black are min == max, so they never trip it.
        for( final Color[] row : tcs._color_schemes ) {
            for( final Color c : row ) {
                final int max = Math.max( c.getRed(), Math.max( c.getGreen(), c.getBlue() ) );
                final int min = Math.min( c.getRed(), Math.min( c.getGreen(), c.getBlue() ) );
                if ( ( max == 255 ) && ( min <= 110 ) ) {
                    System.out.println( "  [TreeColorSchemeTest] raw primary in a scheme row: " + c );
                    return false;
                }
            }
        }
        // ...and every slot a row carries must actually be READ by setColorSchema -- one name per colour, no
        // dead slots left for a future getter to resurrect
        if ( tcs._color_schemes[ 0 ].length != TreeColorSet.COLOR_FIELDS.length ) {
            System.out.println( "  [TreeColorSchemeTest] a scheme row and COLOR_FIELDS must agree: "
                    + tcs._color_schemes[ 0 ].length + " vs " + TreeColorSet.COLOR_FIELDS.length );
            return false;
        }
        // ---- found-node highlights: UNIFIED across themes (A = red, "both" = teal, same hue per role in Light and
        //      Dark -- no longer flipped), with Search B = the "Found/Selected Colors" choice ----
        tcs.setFoundColor( Options.FOUND_COLOR.ELECTRIC_VIOLET );
        tcs.setColorSchema( 1 ); // Light: A red, B electric violet, both teal
        if ( !new Color( 0xE0, 0x00, 0x00 ).equals( tcs.getFoundColor0() )
                || !new Color( 0x9D, 0x00, 0xFF ).equals( tcs.getFoundColor1() )
                || !new Color( 0x00, 0x8B, 0x8B ).equals( tcs.getFoundColor0and1() ) ) {
            return false;
        }
        tcs.setColorSchema( 0 ); // Dark: SAME hues, brightened for the dark ground
        if ( !new Color( 0xFF, 0x40, 0x40 ).equals( tcs.getFoundColor0() )
                || !new Color( 0xB4, 0x4D, 0xFF ).equals( tcs.getFoundColor1() )
                || !new Color( 0x00, 0xC8, 0xC8 ).equals( tcs.getFoundColor0and1() ) ) {
            return false;
        }
        // A stays red-dominant and "both" stays teal (green+blue over red) in BOTH themes -- the anti-flip guarantee
        for ( final int s : new int[] { 0, 1 } ) {
            tcs.setColorSchema( s );
            final Color a = tcs.getFoundColor0();
            final Color both = tcs.getFoundColor0and1();
            if ( ( a.getRed() <= a.getGreen() ) || ( a.getRed() <= a.getBlue() ) ) {
                return false; // A must be red
            }
            if ( ( both.getRed() >= both.getGreen() ) || ( both.getRed() >= both.getBlue() ) ) {
                return false; // "both" must be teal
            }
        }
        // only Search B changes with the choice; A and "both" do not
        tcs.setColorSchema( 1 );
        tcs.setFoundColor( Options.FOUND_COLOR.NEON_MAGENTA );
        if ( !new Color( 0xFF, 0x00, 0xAA ).equals( tcs.getFoundColor1() )
                || !new Color( 0xE0, 0x00, 0x00 ).equals( tcs.getFoundColor0() )
                || !new Color( 0x00, 0x8B, 0x8B ).equals( tcs.getFoundColor0and1() ) ) {
            return false;
        }
        tcs.setFoundColor( Options.FOUND_COLOR.EMERALD_GREEN );
        if ( !new Color( 0x00, 0xA0, 0x00 ).equals( tcs.getFoundColor1() ) ) {
            return false;
        }
        // ---- reconciliation event colors: Okabe-Ito (colorblind-safe), same in both themes; NOT the old pure
        //      red/green/yellow primaries (duplication-vs-speciation was the classic red-green confusion pair) ----
        for ( final int s : new int[] { 0, 1 } ) {
            tcs.setColorSchema( s );
            if ( !new Color( 0xD5, 0x5E, 0x00 ).equals( tcs.getDuplicationBoxColor() )      // vermillion
                    || !new Color( 0x00, 0x9E, 0x73 ).equals( tcs.getSpecBoxColor() )        // bluish-green
                    || !new Color( 0xE6, 0x9F, 0x00 ).equals( tcs.getDuplicationOrSpeciationColor() ) ) { // amber
                return false;
            }
            // duplication and speciation must be distinguishable under red-green color blindness: the old pair
            // (pure red 255,0,0 vs pure green 0,255,0) shared a blue channel of 0 and differed only in red-vs-green;
            // vermillion vs bluish-green differ strongly in the BLUE channel too, which colorblind vision retains
            if ( Math.abs( tcs.getDuplicationBoxColor().getBlue() - tcs.getSpecBoxColor().getBlue() ) < 60 ) {
                return false;
            }
        }
        return true;
    }

    private TreeColorSchemeTest() {
    }
}
