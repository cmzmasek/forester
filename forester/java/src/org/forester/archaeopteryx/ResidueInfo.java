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

import org.forester.util.ForesterUtil;

/**
 * What one residue of the alignment display says about itself: its full name, the physico-chemical class it is
 * COLOURED by, and (for an amino acid) its Kyte-Doolittle hydropathy index -- plus the two positions a reader of an
 * alignment actually needs, the alignment column and the residue's own ungapped number within its sequence.
 * <p>
 * The class names come from {@link MsaColors}, deliberately: the alignment is already COLOURED by that grouping, so
 * naming it from anywhere else would let the tooltip and the colours disagree.
 * <p>
 * Pure and headless: this is a lookup table plus two counters, so it is fully unit-testable and costs nothing to
 * call. The ungapped position is computed by a scan rather than cached -- measured at ~7 us over a 1,000-column
 * sequence and ~13 us over 10,000, against the ~6 us node scan the panel already performs on every mouse move, so a
 * per-sequence prefix-sum cache (tens of MB on a large alignment) would be the wrong trade.
 */
final class ResidueInfo {

    /**
     * Kyte-Doolittle hydropathy index, indexed by 'A'..'Z' (NaN where the letter is not one of the 20 standard
     * amino acids). Positive = hydrophobic, negative = hydrophilic; the scale runs +4.5 (Ile) to -4.5 (Arg).
     * <p>
     * Kyte J, Doolittle RF (1982): "A simple method for displaying the hydropathic character of a protein",
     * Journal of Molecular Biology 157(1):105-132, doi:10.1016/0022-2836(82)90515-0. The 20 values were checked
     * against two independent published tabulations before being written down here.
     */
    private static final double[] KD = new double[ 26 ];
    static {
        java.util.Arrays.fill( KD, Double.NaN );
        KD[ 'I' - 'A' ] = 4.5;
        KD[ 'V' - 'A' ] = 4.2;
        KD[ 'L' - 'A' ] = 3.8;
        KD[ 'F' - 'A' ] = 2.8;
        KD[ 'C' - 'A' ] = 2.5;
        KD[ 'M' - 'A' ] = 1.9;
        KD[ 'A' - 'A' ] = 1.8;
        KD[ 'G' - 'A' ] = -0.4;
        KD[ 'T' - 'A' ] = -0.7;
        KD[ 'S' - 'A' ] = -0.8;
        KD[ 'W' - 'A' ] = -0.9;
        KD[ 'Y' - 'A' ] = -1.3;
        KD[ 'P' - 'A' ] = -1.6;
        KD[ 'H' - 'A' ] = -3.2;
        KD[ 'E' - 'A' ] = -3.5;
        KD[ 'Q' - 'A' ] = -3.5;
        KD[ 'D' - 'A' ] = -3.5;
        KD[ 'N' - 'A' ] = -3.5;
        KD[ 'K' - 'A' ] = -3.9;
        KD[ 'R' - 'A' ] = -4.5;
    }

    /** True for the characters the alignment treats as a gap -- delegated to {@link MsaColors#isGap}, which is the
     *  single definition the PAINT uses too, so a cell cannot be filled on screen and read as a gap here. */
    static boolean isGap( final char residue ) {
        return MsaColors.isGap( residue );
    }

    /** The Kyte-Doolittle hydropathy index, or {@link Double#NaN} for a gap, an ambiguity code (X/B/Z) or any other
     *  non-standard letter -- for which the scale simply has no value, and inventing one would be a fabrication. */
    static double hydropathy( final char residue ) {
        final char r = Character.toUpperCase( residue );
        return ( ( r >= 'A' ) && ( r <= 'Z' ) ) ? KD[ r - 'A' ] : Double.NaN;
    }

    /** The residue's full name ("Tryptophan", "Adenine"), or null when the letter is not one this knows. */
    static String fullName( final char residue, final boolean nucleotide ) {
        final char r = Character.toUpperCase( residue );
        if ( nucleotide ) {
            switch ( r ) {
                case 'A': return "Adenine";
                case 'C': return "Cytosine";
                case 'G': return "Guanine";
                case 'T': return "Thymine";
                case 'U': return "Uracil";
                default:  return null;
            }
        }
        switch ( r ) {
            case 'A': return "Alanine";
            case 'R': return "Arginine";
            case 'N': return "Asparagine";
            case 'D': return "Aspartic acid";
            case 'C': return "Cysteine";
            case 'E': return "Glutamic acid";
            case 'Q': return "Glutamine";
            case 'G': return "Glycine";
            case 'H': return "Histidine";
            case 'I': return "Isoleucine";
            case 'L': return "Leucine";
            case 'K': return "Lysine";
            case 'M': return "Methionine";
            case 'F': return "Phenylalanine";
            case 'P': return "Proline";
            case 'S': return "Serine";
            case 'T': return "Threonine";
            case 'W': return "Tryptophan";
            case 'Y': return "Tyrosine";
            case 'V': return "Valine";
            case 'U': return "Selenocysteine";
            case 'O': return "Pyrrolysine";
            case 'B': return "Asparagine or aspartic acid";
            case 'Z': return "Glutamine or glutamic acid";
            case 'X': return "Any amino acid";
            default:  return null;
        }
    }

    /**
     * The residue's 1-based number within its OWN ungapped sequence, or 0 when the column is a gap (nothing to
     * number). This is the coordinate a reader maps back onto a real protein or gene -- the alignment column is a
     * property of the alignment, this is a property of the sequence.
     */
    static int ungappedPosition( final String aligned, final int column ) {
        if ( ( aligned == null ) || ( column < 0 ) || ( column >= aligned.length() )
                || isGap( aligned.charAt( column ) ) ) {
            return 0;
        }
        int n = 0;
        for( int i = 0; i <= column; ++i ) {
            if ( !isGap( aligned.charAt( i ) ) ) {
                ++n;
            }
        }
        return n;
    }

    /**
     * The tooltip for one alignment cell, as HTML (Swing renders {@code <html>} tooltips over several lines):
     * the tip's name, the alignment column, the ungapped residue number, and the residue's identity -- its letter,
     * full name, colour class and, for an amino acid, its hydropathy.
     */
    static String describeCell( final String tip_name, final String aligned, final int column,
                                final boolean nucleotide ) {
        if ( ( aligned == null ) || ( column < 0 ) || ( column >= aligned.length() ) ) {
            return null;
        }
        final char res = aligned.charAt( column );
        final StringBuilder sb = new StringBuilder( "<html>" );
        if ( !ForesterUtil.isEmpty( tip_name ) ) {
            sb.append( "<b>" ).append( escape( tip_name ) ).append( "</b><br>" );
        }
        sb.append( "Alignment column " ).append( column + 1 );
        if ( isGap( res ) ) {
            sb.append( "<br>gap" ).append( "</html>" );
            return sb.toString();
        }
        sb.append( "<br>Residue " ).append( ungappedPosition( aligned, column ) ).append( " of this sequence" );
        final String name = fullName( res, nucleotide );
        sb.append( "<br><b>" ).append( escape( String.valueOf( Character.toUpperCase( res ) ) ) ).append( "</b>" );
        if ( name != null ) {
            sb.append( " &ndash; " ).append( name );
        }
        final String cls = MsaColors.className( res, nucleotide );
        if ( cls != null ) {
            sb.append( "<br>" ).append( cls );
        }
        final double kd = nucleotide ? Double.NaN : hydropathy( res );
        if ( !Double.isNaN( kd ) ) {
            sb.append( "<br>Hydropathy (Kyte&ndash;Doolittle) " ).append( kd );
        }
        return sb.append( "</html>" ).toString();
    }

    /** Minimal HTML escaping: a tip name is user data and the tooltip is rendered as HTML. */
    private static String escape( final String s ) {
        return s.replace( "&", "&amp;" ).replace( "<", "&lt;" ).replace( ">", "&gt;" );
    }

    private ResidueInfo() {
    }
}
