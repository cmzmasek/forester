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
 * The alignment residue readout: the Kyte-Doolittle table, the ungapped position, and the tooltip text.
 * <p>
 * The 20 hydropathy values are SHIPPED SCIENTIFIC CONSTANTS, so they are pinned here one by one rather than spot
 * checked -- a transposed pair would be invisible on screen and wrong in a figure caption.
 */
public final class ResidueInfoTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ResidueInfo: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return hydropathyTable() && names() && classes() && gapAgreement() && ungapped() && tooltip();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ResidueInfoTest] " + msg );
        return false;
    }

    /** Kyte J, Doolittle RF (1982) J Mol Biol 157(1):105-132 -- all 20, cross-checked against two tabulations. */
    private static boolean hydropathyTable() {
        final Object[][] expected = { { 'I', 4.5 }, { 'V', 4.2 }, { 'L', 3.8 }, { 'F', 2.8 }, { 'C', 2.5 },
                                      { 'M', 1.9 }, { 'A', 1.8 }, { 'G', -0.4 }, { 'T', -0.7 }, { 'S', -0.8 },
                                      { 'W', -0.9 }, { 'Y', -1.3 }, { 'P', -1.6 }, { 'H', -3.2 }, { 'E', -3.5 },
                                      { 'Q', -3.5 }, { 'D', -3.5 }, { 'N', -3.5 }, { 'K', -3.9 }, { 'R', -4.5 } };
        for( final Object[] e : expected ) {
            final char c = (Character) e[ 0 ];
            final double want = (Double) e[ 1 ];
            if ( Math.abs( ResidueInfo.hydropathy( c ) - want ) > 1e-9 ) {
                return fail( "hydropathy(" + c + ") should be " + want + ", got " + ResidueInfo.hydropathy( c ) );
            }
            if ( Math.abs( ResidueInfo.hydropathy( Character.toLowerCase( c ) ) - want ) > 1e-9 ) {
                return fail( "hydropathy must be case-insensitive, failed for " + c );
            }
        }
        // the extremes of the scale, so a wholesale sign flip or rescale cannot pass
        if ( ResidueInfo.hydropathy( 'I' ) <= ResidueInfo.hydropathy( 'R' ) ) {
            return fail( "Ile must be the most hydrophobic and Arg the most hydrophilic" );
        }
        // a letter the scale does not define must yield NaN -- NOT 0, which would read as 'neutral'
        for( final char c : new char[] { 'X', 'B', 'Z', 'O', '-', '.', '*', '1' } ) {
            if ( !Double.isNaN( ResidueInfo.hydropathy( c ) ) ) {
                return fail( "the scale has no value for '" + c + "'; it must be NaN, not "
                        + ResidueInfo.hydropathy( c ) );
            }
        }
        return true;
    }

    private static boolean names() {
        if ( !"Tryptophan".equals( ResidueInfo.fullName( 'W', false ) )
                || !"Tryptophan".equals( ResidueInfo.fullName( 'w', false ) )
                || !"Aspartic acid".equals( ResidueInfo.fullName( 'D', false ) )
                || !"Selenocysteine".equals( ResidueInfo.fullName( 'U', false ) ) ) {
            return fail( "amino-acid names are wrong" );
        }
        // 'U' is selenocysteine in a protein but URACIL in a nucleotide alignment -- the flag must decide
        if ( !"Uracil".equals( ResidueInfo.fullName( 'U', true ) )
                || !"Adenine".equals( ResidueInfo.fullName( 'A', true ) ) ) {
            return fail( "nucleotide names are wrong" );
        }
        if ( ResidueInfo.fullName( '-', false ) != null ) {
            return fail( "a gap has no name" );
        }
        return true;
    }

    /** The names must come from the SAME grouping the alignment is coloured by, or the tooltip contradicts the
     *  colours: every residue sharing a colour must share a class name. */
    private static boolean classes() {
        final String[][] same_colour = { { "ILVAM" }, { "FWY" }, { "KRH" }, { "DE" }, { "STNQ" }, { "PG" } };
        for( final String[] grp : same_colour ) {
            final String first = MsaColors.className( grp[ 0 ].charAt( 0 ), false );
            for( final char c : grp[ 0 ].toCharArray() ) {
                if ( !first.equals( MsaColors.className( c, false ) ) ) {
                    return fail( "'" + c + "' is coloured with '" + grp[ 0 ].charAt( 0 )
                            + "' but named differently: " + MsaColors.className( c, false ) + " vs " + first );
                }
            }
        }
        if ( MsaColors.className( '-', false ) != null ) {
            return fail( "a gap has no class" );
        }
        if ( !"purine".equals( MsaColors.className( 'G', true ) )
                || !"pyrimidine".equals( MsaColors.className( 'T', true ) ) ) {
            return fail( "nucleotide classes are wrong" );
        }
        // the SAME letter means different things in the two alphabets
        if ( MsaColors.className( 'G', true ).equals( MsaColors.className( 'G', false ) ) ) {
            return fail( "'G' must not be described identically as a base and as an amino acid" );
        }
        // For NUCLEOTIDES the class is deliberately NOT the colour grouping: the paint gives A/C/G/T four distinct
        // colours while purine/pyrimidine is coarser. Pinned so nobody "fixes" the javadoc back to claiming it is.
        if ( MsaColors.colorFor( 'A', true ).equals( MsaColors.colorFor( 'G', true ) ) ) {
            return fail( "precondition: A and G are painted in DIFFERENT colours" );
        }
        if ( !MsaColors.className( 'A', true ).equals( MsaColors.className( 'G', true ) ) ) {
            return fail( "A and G are both purines -- the base NAME is what distinguishes them, not the class" );
        }
        if ( ResidueInfo.fullName( 'A', true ).equals( ResidueInfo.fullName( 'G', true ) ) ) {
            return fail( "...so the readout must still name the base itself" );
        }
        return true;
    }

    /** The gap set is ONE definition shared by the paint and the readout. When they diverged, a '~' (a terminal
     *  gap in PIR/NBRF and MView output) was painted as a filled grey cell yet reported as "gap" -- and every later
     *  residue's ungapped number was computed as if those columns were not there. */
    private static boolean gapAgreement() {
        for( final char g : new char[] { '-', '.', ' ', '~' } ) {
            if ( !ResidueInfo.isGap( g ) || !MsaColors.isGap( g ) ) {
                return fail( "'" + g + "' must count as a gap" );
            }
            if ( ( MsaColors.colorFor( g, false ) != null ) || ( MsaColors.colorFor( g, true ) != null ) ) {
                return fail( "'" + g + "' reads as a gap but the PAINT still fills its cell" );
            }
            if ( ( MsaColors.className( g, false ) != null ) || ( ResidueInfo.fullName( g, false ) != null ) ) {
                return fail( "'" + g + "' must have no class and no name" );
            }
        }
        // ...and a real residue is neither a gap nor unfilled
        for( final char r : new char[] { 'W', 'X', 'B' } ) {
            if ( ResidueInfo.isGap( r ) || ( MsaColors.colorFor( r, false ) == null ) ) {
                return fail( "'" + r + "' is a residue: it must be filled and must not read as a gap" );
            }
        }
        return true;
    }

    private static boolean ungapped() {
        //            column: 0123456
        final String aligned = "M--KL-A";
        if ( ResidueInfo.ungappedPosition( aligned, 0 ) != 1 ) {
            return fail( "the first residue is number 1" );
        }
        if ( ResidueInfo.ungappedPosition( aligned, 3 ) != 2 ) {
            return fail( "gaps must not be counted: K is residue 2, got "
                    + ResidueInfo.ungappedPosition( aligned, 3 ) );
        }
        if ( ResidueInfo.ungappedPosition( aligned, 4 ) != 3 ) {
            return fail( "L is residue 3, got " + ResidueInfo.ungappedPosition( aligned, 4 ) );
        }
        if ( ResidueInfo.ungappedPosition( aligned, 6 ) != 4 ) {
            return fail( "A is residue 4, got " + ResidueInfo.ungappedPosition( aligned, 6 ) );
        }
        // a gap column has no residue number at all
        for( final int gap_col : new int[] { 1, 2, 5 } ) {
            if ( ResidueInfo.ungappedPosition( aligned, gap_col ) != 0 ) {
                return fail( "a gap column has no residue number, got "
                        + ResidueInfo.ungappedPosition( aligned, gap_col ) + " at " + gap_col );
            }
        }
        if ( ( ResidueInfo.ungappedPosition( aligned, -1 ) != 0 )
                || ( ResidueInfo.ungappedPosition( aligned, 99 ) != 0 )
                || ( ResidueInfo.ungappedPosition( null, 0 ) != 0 ) ) {
            return fail( "out-of-range and null must be handled, not thrown" );
        }
        // every gap convention the display honours: in "A.B~C D" the residues are A, B, C, D -- so D is number 4
        if ( ResidueInfo.ungappedPosition( "A.B~C D", 6 ) != 4 ) {
            return fail( "'.', '~' and ' ' are gaps too, got " + ResidueInfo.ungappedPosition( "A.B~C D", 6 ) );
        }
        if ( ResidueInfo.ungappedPosition( "A.B~C D", 3 ) != 0 ) {
            return fail( "'~' is a gap column, got " + ResidueInfo.ungappedPosition( "A.B~C D", 3 ) );
        }
        return true;
    }

    private static boolean tooltip() {
        final String t = ResidueInfo.describeCell( "seq_1", "M--KL-A", 4, false );
        if ( ( t == null ) || !t.startsWith( "<html>" ) || !t.endsWith( "</html>" ) ) {
            return fail( "the tooltip should be HTML, got " + t );
        }
        if ( !t.contains( "seq_1" ) || !t.contains( "Alignment column 5" )
                || !t.contains( "Residue 3 of this sequence" ) || !t.contains( "Leucine" )
                || !t.contains( "aliphatic" ) || !t.contains( "3.8" ) ) {
            return fail( "the tooltip is missing something: " + t );
        }
        // columns are 1-based on screen but 0-based in the code -- an off-by-one here misreports every position
        if ( t.contains( "column 4" ) ) {
            return fail( "the alignment column must be shown 1-based: " + t );
        }
        final String gap = ResidueInfo.describeCell( "seq_1", "M--KL-A", 1, false );
        if ( ( gap == null ) || !gap.contains( "gap" ) || gap.contains( "Residue" ) || gap.contains( "Hydropathy" ) ) {
            return fail( "a gap cell says gap and nothing else: " + gap );
        }
        // a nucleotide alignment must NOT claim a hydropathy (the scale is for amino acids)
        final String nuc = ResidueInfo.describeCell( "seq_1", "ACGT", 3, true );
        if ( ( nuc == null ) || !nuc.contains( "Thymine" ) || !nuc.contains( "pyrimidine" )
                || nuc.contains( "Hydropathy" ) ) {
            return fail( "a base has no hydropathy: " + nuc );
        }
        // an ambiguity code has a name and a class but no hydropathy
        final String x = ResidueInfo.describeCell( "s", "X", 0, false );
        if ( ( x == null ) || !x.contains( "Any amino acid" ) || x.contains( "Hydropathy" ) ) {
            return fail( "X has no hydropathy value: " + x );
        }
        // the tip name is user data rendered as HTML
        final String esc = ResidueInfo.describeCell( "a<b>&c", "A", 0, false );
        if ( ( esc == null ) || esc.contains( "<b>&c" ) || !esc.contains( "a&lt;b&gt;&amp;c" ) ) {
            return fail( "the tip name must be HTML-escaped: " + esc );
        }
        if ( ResidueInfo.describeCell( "s", "ACGT", 9, false ) != null ) {
            return fail( "a column past the end of the sequence has no tooltip" );
        }
        return true;
    }

    private ResidueInfoTest() {
    }
}
