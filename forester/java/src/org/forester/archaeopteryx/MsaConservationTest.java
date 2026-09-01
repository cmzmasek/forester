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

import java.util.Arrays;
import java.util.List;

import org.forester.archaeopteryx.MsaConservation.Measure;

/**
 * The conservation measures, against values worked out by hand.
 * <p>
 * A conservation bar is read off a published figure, so the numbers behind it have to be exactly what the
 * documentation says they are -- an off-by-a-gap denominator would be invisible on screen and wrong in a caption.
 * Every expected value below is derived in the comment beside it.
 */
public final class MsaConservationTest {

    private static final double EPS = 1e-6;
    /** ln(20) and ln(4): the amino-acid and nucleotide alphabet sizes the information content normalises by. */
    private static final double LN20 = Math.log( 20 );
    private static final double LN4  = Math.log( 4 );

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MsaConservation: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return identity() && information() && gapsAndEdges() && conventions();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [MsaConservationTest] " + msg );
        return false;
    }

    private static boolean near( final String what, final double got, final double expected ) {
        if ( Math.abs( got - expected ) > EPS ) {
            return fail( what + ": expected " + expected + ", got " + got );
        }
        return true;
    }

    private static MsaConservation protein( final String... rows ) {
        int len = 0;
        for( final String r : rows ) {
            len = Math.max( len, r.length() );
        }
        return MsaConservation.compute( Arrays.asList( rows ), len, false );
    }

    private static boolean identity() {
        final MsaConservation c = protein( "MKTA", "MKTA", "MRTA", "MKSG" );
        if ( c.rows() != 4 ) {
            return fail( "rows: " + c.rows() );
        }
        if ( c.length() != 4 ) {
            return fail( "length: " + c.length() );
        }
        if ( !near( "col 0 identity (M,M,M,M = 4/4)", c.identityAt( 0 ), 1.0 ) ) {
            return false;
        }
        if ( !near( "col 1 identity (K,K,R,K = 3/4)", c.identityAt( 1 ), 0.75 ) ) {
            return false;
        }
        if ( ( c.consensusAt( 0 ) != 'M' ) || ( c.consensusAt( 1 ) != 'K' ) || ( c.consensusAt( 2 ) != 'T' )
                || ( c.consensusAt( 3 ) != 'A' ) ) {
            return fail( "consensus: " + c.consensusAt( 0 ) + c.consensusAt( 1 ) + c.consensusAt( 2 )
                    + c.consensusAt( 3 ) );
        }
        return true;
    }

    private static boolean information() {
        // K,K,R,K -> p = 3/4, 1/4; H = -(0.75 ln0.75 + 0.25 ln0.25) = 0.5623351 nats
        final double h = -( ( 0.75 * Math.log( 0.75 ) ) + ( 0.25 * Math.log( 0.25 ) ) );
        final MsaConservation c = protein( "MKTA", "MKTA", "MRTA", "MKSG" );
        if ( !near( "col 1 information", c.informationAt( 1 ), 1.0 - ( h / LN20 ) ) ) {
            return false;
        }
        if ( !near( "a fully conserved column is 1.0", c.informationAt( 0 ), 1.0 ) ) {
            return false;
        }
        // THE reason both measures are offered: two columns with the SAME majority fraction, different diversity.
        // Eight rows, majority 4/8 in both -> identity ties at 0.5, information separates them.
        final MsaConservation two_way = protein( "A", "A", "A", "A", "C", "C", "C", "C" );
        final MsaConservation four_way = protein( "A", "A", "A", "A", "C", "C", "D", "E" );
        if ( !near( "two-way identity", two_way.identityAt( 0 ), 0.5 )
                || !near( "four-way identity", four_way.identityAt( 0 ), 0.5 ) ) {
            return false;
        }
        final double h2 = Math.log( 2 ); // p = 1/2, 1/2
        // p = 1/2, 1/4, 1/8, 1/8
        final double h4 = -( ( 0.5 * Math.log( 0.5 ) ) + ( 0.25 * Math.log( 0.25 ) )
                + ( 2 * ( 0.125 * Math.log( 0.125 ) ) ) );
        if ( !near( "two-way information", two_way.informationAt( 0 ), 1.0 - ( h2 / LN20 ) )
                || !near( "four-way information", four_way.informationAt( 0 ), 1.0 - ( h4 / LN20 ) ) ) {
            return false;
        }
        if ( two_way.informationAt( 0 ) <= four_way.informationAt( 0 ) ) {
            return fail( "information must rank a two-way split above a four-way one at equal identity" );
        }
        // the alphabet size matters: the same column normalises against log(4), not log(20), for nucleotides
        final MsaConservation nuc = MsaConservation
                .compute( Arrays.asList( "A", "A", "C", "C" ), 1, true );
        if ( !near( "nucleotide information", nuc.informationAt( 0 ), 1.0 - ( h2 / LN4 ) ) ) {
            return false;
        }
        final MsaConservation as_protein = protein( "A", "A", "C", "C" );
        if ( nuc.informationAt( 0 ) >= as_protein.informationAt( 0 ) ) {
            return fail( "a 2-of-4-letter split is less surprising than a 2-of-20 one, so it must score LOWER" );
        }
        return true;
    }

    private static boolean gapsAndEdges() {
        // gaps count in the DENOMINATOR: one residue among four rows scores 0.25 under both measures, even though
        // the residues present are perfectly conserved
        final MsaConservation c = protein( "AA", "A-", "A-", "A-" );
        if ( !near( "gappy column identity", c.identityAt( 1 ), 0.25 )
                || !near( "gappy column information", c.informationAt( 1 ), 0.25 ) ) {
            return false;
        }
        if ( c.consensusAt( 1 ) != 'A' ) {
            return fail( "a gappy column still names its consensus residue, got " + (int) c.consensusAt( 1 ) );
        }
        // an all-gap column scores zero and names nothing
        final MsaConservation all_gap = protein( "A-", "A-", "A." );
        if ( ( all_gap.identityAt( 1 ) != 0 ) || ( all_gap.informationAt( 1 ) != 0 )
                || ( all_gap.consensusAt( 1 ) != 0 ) ) {
            return fail( "an all-gap column must score 0 with no consensus" );
        }
        // a row SHORTER than the alignment counts as gapped past its end -- never as absent from the denominator
        final MsaConservation ragged = protein( "AAAA", "AA" );
        if ( !near( "short row counts as a gap", ragged.identityAt( 2 ), 0.5 ) ) {
            return false;
        }
        // out of range and empty input are answers, not exceptions
        if ( ( c.scoreAt( -1, Measure.IDENTITY ) != 0 ) || ( c.scoreAt( 99, Measure.INFORMATION ) != 0 )
                || ( c.consensusAt( 99 ) != 0 ) ) {
            return fail( "out-of-range columns must score 0" );
        }
        final MsaConservation empty = MsaConservation.compute( null, 5, false );
        if ( ( empty.rows() != 0 ) || ( empty.length() != 5 ) || ( empty.identityAt( 0 ) != 0 ) ) {
            return fail( "no rows -> no scores, no exception" );
        }
        final List<String> none = Arrays.asList();
        if ( MsaConservation.compute( none, 0, false ).length() != 0 ) {
            return fail( "zero columns -> zero length" );
        }
        return true;
    }

    private static boolean conventions() {
        // case-insensitive: the alignment is drawn and coloured that way, so it must be scored that way
        final MsaConservation mixed = protein( "a", "A", "A" );
        if ( !near( "case-insensitive identity", mixed.identityAt( 0 ), 1.0 ) || ( mixed.consensusAt( 0 ) != 'A' ) ) {
            return fail( "lower and upper case are the same residue" );
        }
        // every gap character the display recognises counts as a gap here too (ONE definition, MsaColors.isGap)
        final MsaConservation gaps = protein( "A", "-", ".", "~" );
        if ( !near( "all gap characters count as gaps", gaps.identityAt( 0 ), 0.25 ) ) {
            return false;
        }
        // ties go to the alphabetically first residue, so the same alignment always yields the same figure
        final MsaConservation tie = protein( "C", "A" );
        if ( tie.consensusAt( 0 ) != 'A' ) {
            return fail( "a tie must resolve to the alphabetically first residue, got " + tie.consensusAt( 0 ) );
        }
        // more distinct symbols than the alphabet (ambiguity codes) would drive the raw value below zero: clamped,
        // and low -- never silently treated as conserved
        final String[] many = new String[ 25 ];
        for( int i = 0; i < many.length; i++ ) {
            many[ i ] = String.valueOf( (char) ( 'A' + i ) );
        }
        final MsaConservation over = protein( many );
        if ( ( over.informationAt( 0 ) < 0 ) || ( over.informationAt( 0 ) > 0.01 ) ) {
            return fail( "25 distinct residues must clamp to ~0, got " + over.informationAt( 0 ) );
        }
        // scoreAt dispatches to the measure it is given
        final MsaConservation c = protein( "MKTA", "MKTA", "MRTA", "MKSG" );
        if ( ( c.scoreAt( 1, Measure.IDENTITY ) != c.identityAt( 1 ) )
                || ( c.scoreAt( 1, Measure.INFORMATION ) != c.informationAt( 1 ) ) ) {
            return fail( "scoreAt must dispatch on the measure" );
        }
        return true;
    }

    private MsaConservationTest() {
    }
}
