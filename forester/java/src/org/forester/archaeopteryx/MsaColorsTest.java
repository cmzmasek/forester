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

/** Headless tests for the alignment residue coloring ({@link MsaColors}). */
public final class MsaColorsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MsaColors: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            // type detection: protein vs nucleotide (over the ACTUAL residues)
            if ( MsaColors.isNucleotide( "MKLVWEHD" ) ) {
                return fail( "a protein sequence must not read as nucleotide" );
            }
            if ( !MsaColors.isNucleotide( "ACGTACGTTT" ) ) {
                return fail( "a DNA sequence must read as nucleotide" );
            }
            if ( !MsaColors.isNucleotide( "ACGUACGUUU" ) ) {
                return fail( "an RNA sequence must read as nucleotide" );
            }

            // a gap has no color, in either scheme
            if ( ( MsaColors.colorFor( '-', false ) != null ) || ( MsaColors.colorFor( '-', true ) != null )
                    || ( MsaColors.colorFor( '.', false ) != null ) ) {
                return fail( "a gap must have no color" );
            }

            // amino acids: same physico-chemical group shares a color; different groups differ
            final Color asp = MsaColors.colorFor( 'D', false );
            final Color glu = MsaColors.colorFor( 'E', false ); // both negative
            final Color lys = MsaColors.colorFor( 'K', false ); // positive
            if ( ( asp == null ) || !asp.equals( glu ) ) {
                return fail( "D and E (both negative) must share a color" );
            }
            if ( asp.equals( lys ) ) {
                return fail( "a negative and a positive residue must not share a color" );
            }
            if ( MsaColors.colorFor( 'i', false ) == null || !MsaColors.colorFor( 'i', false )
                    .equals( MsaColors.colorFor( 'I', false ) ) ) {
                return fail( "coloring must be case-insensitive" );
            }

            // nucleotides: A/C/G/T distinct; U colors as T
            final Color a = MsaColors.colorFor( 'A', true );
            final Color c = MsaColors.colorFor( 'C', true );
            final Color g = MsaColors.colorFor( 'G', true );
            final Color t = MsaColors.colorFor( 'T', true );
            if ( ( a == null ) || a.equals( c ) || a.equals( g ) || a.equals( t ) || c.equals( g ) ) {
                return fail( "A/C/G/T must be visually distinct" );
            }
            if ( !MsaColors.colorFor( 'U', true ).equals( t ) ) {
                return fail( "U must color the same as T" );
            }

            // letter contrast: light bg -> dark ink, dark bg -> light ink
            if ( !MsaColors.letterColor( new Color( 240, 240, 240 ) ).equals( Color.BLACK )
                    || !MsaColors.letterColor( new Color( 20, 20, 90 ) ).equals( Color.WHITE ) ) {
                return fail( "letter color must contrast the cell background" );
            }
            return true;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [MsaColorsTest] " + message );
        return false;
    }

    private MsaColorsTest() {
        // not instantiable
    }
}
