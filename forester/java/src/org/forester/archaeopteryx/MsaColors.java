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

import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterUtil;

/**
 * Per-residue colors for the sequence-alignment display. Amino acids use a fixed physico-chemical (Zappo-style)
 * scheme -- aliphatic / aromatic / positive / negative / hydrophilic / conformationally-special / cysteine; nucleotides
 * use the standard A/C/G/T-U convention. A gap (or blank) has no color (drawn as background). The scheme is a fixed,
 * residue-identity mapping (not conservation-dependent), so it is O(1) per cell -- fast enough for wide alignments.
 */
final class MsaColors {

    static final char GAP = MolecularSequence.GAP; // '-'

    // ----- amino-acid (Zappo-style physico-chemical) -----
    private static final Color AA_ALIPHATIC = new Color( 240, 170, 170 ); // I L V A M
    private static final Color AA_AROMATIC  = new Color( 240, 190, 90 );  // F W Y
    private static final Color AA_POSITIVE  = new Color( 120, 130, 240 ); // K R H
    private static final Color AA_NEGATIVE  = new Color( 230, 100, 100 ); // D E
    private static final Color AA_HYDROPHIL = new Color( 120, 200, 120 ); // S T N Q
    private static final Color AA_SPECIAL   = new Color( 220, 130, 220 ); // P G  (conformationally special)
    private static final Color AA_CYSTEINE  = new Color( 235, 220, 110 ); // C

    // ----- nucleotide -----
    private static final Color NUC_A = new Color( 120, 200, 120 ); // green
    private static final Color NUC_C = new Color( 120, 130, 240 ); // blue
    private static final Color NUC_G = new Color( 230, 185, 80 );  // amber
    private static final Color NUC_T = new Color( 230, 110, 110 ); // red (also U)

    /** A muted grey for a recognized-but-uncategorized residue (X, B, Z, *, ...), so it is visible but not loud. */
    private static final Color UNKNOWN = new Color( 205, 205, 205 );

    /**
     * Whether {@code sample} residues read as nucleotide (DNA/RNA) rather than amino acid. Uses
     * {@link ForesterUtil#guessMolecularSequenceType(String)} on the ACTUAL residues -- NOT {@code Msa.getType()},
     * which the FASTA parser hard-codes to AA. An undecidable/empty sample defaults to amino acid (false).
     */
    static boolean isNucleotide( final String sample ) {
        final MolecularSequence.TYPE t = ForesterUtil.guessMolecularSequenceType( sample == null ? "" : sample );
        return ( t == MolecularSequence.TYPE.DNA ) || ( t == MolecularSequence.TYPE.RNA );
    }

    /** The fill color for {@code residue}, or {@code null} for a gap / blank (which is drawn as background). */
    static Color colorFor( final char residue, final boolean nucleotide ) {
        final char r = Character.toUpperCase( residue );
        if ( ( r == GAP ) || ( r == '.' ) || ( r == ' ' ) ) {
            return null; // a gap has no fill
        }
        return nucleotide ? nucleotideColor( r ) : aminoAcidColor( r );
    }

    private static Color aminoAcidColor( final char r ) {
        switch ( r ) {
            case 'I': case 'L': case 'V': case 'A': case 'M':
                return AA_ALIPHATIC;
            case 'F': case 'W': case 'Y':
                return AA_AROMATIC;
            case 'K': case 'R': case 'H':
                return AA_POSITIVE;
            case 'D': case 'E':
                return AA_NEGATIVE;
            case 'S': case 'T': case 'N': case 'Q':
                return AA_HYDROPHIL;
            case 'P': case 'G':
                return AA_SPECIAL;
            case 'C':
                return AA_CYSTEINE;
            default:
                return UNKNOWN; // X, B, Z, *, U (selenocysteine), etc.
        }
    }

    private static Color nucleotideColor( final char r ) {
        switch ( r ) {
            case 'A':
                return NUC_A;
            case 'C':
                return NUC_C;
            case 'G':
                return NUC_G;
            case 'T': case 'U':
                return NUC_T;
            default:
                return UNKNOWN; // N, R, Y, ambiguity codes
        }
    }

    /** Black or white, whichever reads better on {@code bg} -- for drawing the residue letter over its cell. */
    static Color letterColor( final Color bg ) {
        if ( bg == null ) {
            return Color.DARK_GRAY;
        }
        final double luminance = ( 0.299 * bg.getRed() ) + ( 0.587 * bg.getGreen() ) + ( 0.114 * bg.getBlue() );
        return ( luminance < 140 ) ? Color.WHITE : Color.BLACK;
    }

    private MsaColors() {
        // not instantiable
    }
}
