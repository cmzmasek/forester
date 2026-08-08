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

import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Properties;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;

/**
 * Guards the removal of the "Curved" and "Convex" tree styles: the two enum constants are gone (while the surviving
 * styles remain), and an old settings file that stored a removed style is tolerated on load (the unknown value is
 * ignored, while a still-valid stored style IS applied). Headless.
 */
public final class TreeStyleRemovalTest {

    // the styles that must survive the removal (asserted individually, so adding a NEW style later does not
    // spuriously fail this removal-specific test)
    private static final String[] SURVIVING = { "RECTANGULAR", "EURO_STYLE", "ROUNDED", "TRIANGULAR", "UNROOTED",
                                                "CIRCULAR" };

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreeStyleRemoval: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return enumRemoved() && persistence();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected: " + e );
        }
    }

    /** CURVED and CONVEX are gone; the surviving styles are still present. */
    private static boolean enumRemoved() {
        for( final String gone : new String[] { "CURVED", "CONVEX" } ) {
            try {
                PHYLOGENY_GRAPHICS_TYPE.valueOf( gone );
                return fail( gone + " should have been removed from PHYLOGENY_GRAPHICS_TYPE" );
            }
            catch ( final IllegalArgumentException expected ) {
                // good -- the constant is gone
            }
        }
        for( final String live : SURVIVING ) {
            try {
                PHYLOGENY_GRAPHICS_TYPE.valueOf( live );
            }
            catch ( final IllegalArgumentException e ) {
                return fail( "surviving style " + live + " is unexpectedly missing" );
            }
        }
        return true;
    }

    /** An old settings file that stored the now-removed style is ignored on load (the option keeps its value); a
     *  positive control (a still-valid style IS applied) proves the pref is actually wired, so the ignore is meaningful. */
    private static boolean persistence() throws Exception {
        return appliesStored( "CURVED", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR )
                && appliesStored( "ROUNDED", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, PHYLOGENY_GRAPHICS_TYPE.ROUNDED );
    }

    private static boolean appliesStored( final String stored, final PHYLOGENY_GRAPHICS_TYPE initial,
                                          final PHYLOGENY_GRAPHICS_TYPE expected ) throws Exception {
        final Path file = Files.createTempFile( "aptx-prefs", ".properties" );
        final Properties p = new Properties();
        p.setProperty( "phylogeny_graphics_type", stored );
        try ( final OutputStream o = Files.newOutputStream( file ) ) {
            p.store( o, "test" );
        }
        final Options opts = Options.createInstance( new Configuration() );
        opts.setPhylogenyGraphicsType( initial );
        new GuiPreferences( file ).applyTo( opts ); // must not throw
        if ( opts.getPhylogenyGraphicsType() != expected ) {
            return fail( "stored graphics type '" + stored + "' -> expected " + expected + ", got "
                    + opts.getPhylogenyGraphicsType() );
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [TreeStyleRemovalTest] " + msg );
        return false;
    }

    private TreeStyleRemovalTest() {
    }
}
