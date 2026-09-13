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

package org.forester.archaeopteryx.phylogeny.data;

/**
 * Where a protein domain is drawn along its architecture's track: residue r covers {@code [start + (r-1) f, start + r f]},
 * so a domain {@code from..to} spans {@code [start + (from-1) f, start + to f]} and a domain {@code 1..L} covers the
 * backbone exactly. JOINT rule with Archaeopteryx.js (Christian, 2026-09-12); the desktop used to place a domain one
 * residue to the right.
 */
public final class RenderableDomainArchitectureTest {

    private static final float EPS = 1e-3f;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RenderableDomainArchitecture: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final float start = 100f;
        final float f = 0.129808f; // the apaf.xml scale at a 1200 px viewport: 300 / 2080 x 0.9
        // a domain 1..L covers the backbone [start, start + L f] exactly -- no residue of offset at either end
        final float[] whole = RenderableDomainArchitecture.domainExtent( start, 1, 1249, f );
        if ( ( whole[ 0 ] != start ) || !near( whole[ 1 ], start + ( 1249 * f ) ) ) {
            return fail( "a domain 1..L must cover the backbone exactly, got [" + whole[ 0 ] + ", " + whole[ 1 ] + "]" );
        }
        // the spec's worked numbers for 22_MOUSE: CARD 6-90 at x offset 0.65, width 11.03
        final float[] card = RenderableDomainArchitecture.domainExtent( start, 6, 90, f );
        if ( !near( card[ 0 ] - start, 5 * f ) || !near( card[ 1 ] - card[ 0 ], 85 * f ) ) {
            return fail( "CARD 6-90 must start 5 residues in and span 85, got offset " + ( card[ 0 ] - start ) + " width "
                    + ( card[ 1 ] - card[ 0 ] ) );
        }
        // a single residue r is exactly one residue wide, starting (r-1) residues in
        final float[] one = RenderableDomainArchitecture.domainExtent( start, 1, 1, f );
        if ( ( one[ 0 ] != start ) || !near( one[ 1 ] - one[ 0 ], f ) ) {
            return fail( "residue 1 must occupy [start, start + f]" );
        }
        // adjacent domains a..b and b+1..c meet without a gap or an overlap
        final float[] left = RenderableDomainArchitecture.domainExtent( start, 605, 643, f );
        final float[] right = RenderableDomainArchitecture.domainExtent( start, 644, 685, f );
        if ( !near( left[ 1 ], right[ 0 ] ) ) {
            return fail( "adjacent domains must meet, got " + left[ 1 ] + " and " + right[ 0 ] );
        }
        // a domain ending at the last residue ends at the backbone's end
        final float[] last = RenderableDomainArchitecture.domainExtent( start, 1168, 1249, f );
        if ( !near( last[ 1 ], whole[ 1 ] ) ) {
            return fail( "a domain ending at residue L must end where the backbone ends" );
        }
        return true;
    }

    private static boolean near( final float a, final float b ) {
        return Math.abs( a - b ) < EPS;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [RenderableDomainArchitectureTest] " + message );
        return false;
    }

    private RenderableDomainArchitectureTest() {
        // not instantiable
    }
}
