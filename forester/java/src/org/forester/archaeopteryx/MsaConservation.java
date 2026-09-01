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

import java.util.List;

/**
 * Per-column conservation and consensus for an alignment, for the track drawn under the alignment display.
 * <p>
 * Two measures, both on <b>[0,1]</b> so either can drive the same bar:
 * <ul>
 * <li><b>{@link Measure#IDENTITY Consensus identity}</b> -- the fraction of ROWS carrying the column's most common
 * residue. Gaps are counted in the denominator, so a column that is half gaps cannot score above 0.5. This is the
 * same definition as {@code MsaMethods.calculateIdentityRatio}, which the msa_compactor tools use.</li>
 * <li><b>{@link Measure#INFORMATION Information content}</b> -- the sequence-logo measure (Schneider &amp; Stephens
 * 1990), {@code (log2(K) - H) / log2(K)} where {@code H} is the Shannon entropy of the residues actually observed
 * and {@code K} is the alphabet size (4 for nucleotide, 20 for amino acid), multiplied by the non-gap fraction so
 * gaps count against it exactly as they do for identity. It separates cases identity cannot: a column split evenly
 * between two residues scores higher than one split four ways, though both have the same majority fraction.</li>
 * </ul>
 * <p>
 * <b>Deliberately NOT Jalview's conservation</b> (Livingstone &amp; Barton 1993): that scores conserved
 * physico-chemical property sets and is defined for amino acids only, while this display serves nucleotide
 * alignments equally.
 * <p>
 * Conventions, all of which matter for reading a figure:
 * <ul>
 * <li>Residues are compared case-insensitively, matching how the alignment is coloured and drawn.</li>
 * <li>The gap set is {@link MsaColors#isGap} -- one definition shared with the paint, the residue readout and this.</li>
 * <li>A row that is SHORTER than the alignment counts as gapped in the columns past its end (a ragged alignment is
 * malformed; treating the missing tail as gaps is the conservative reading, never the flattering one).</li>
 * <li>The consensus is the most common NON-gap residue, so a gappy column still names its residue -- the bar
 * already carries the gappiness. Ties go to the alphabetically first residue, so a figure is reproducible.</li>
 * <li>Ambiguity codes (N, X, B, Z...) are counted as ordinary distinct symbols. A column full of them therefore
 * scores LOW rather than being silently treated as conserved; the information content is clamped to [0,1] since
 * more than {@code K} distinct symbols can push the raw value below zero.</li>
 * </ul>
 */
final class MsaConservation {

    /** Which measure the bar height shows. {@code toString} is what the Settings dropdown displays. */
    enum Measure {
        IDENTITY( "Consensus identity" ),
        INFORMATION( "Information content" );

        private final String _label;

        Measure( final String label ) {
            _label = label;
        }

        @Override
        public String toString() {
            return _label;
        }
    }

    /** Alphabet sizes used to normalise the information content. */
    private static final int K_NUCLEOTIDE = 4;
    private static final int K_AMINO_ACID = 20;

    private final char[]   _consensus;   // 0 where the column holds no residue at all
    private final double[] _identity;
    private final double[] _information;
    private final int      _rows;

    private MsaConservation( final char[] consensus, final double[] identity, final double[] information,
                             final int rows ) {
        _consensus = consensus;
        _identity = identity;
        _information = information;
        _rows = rows;
    }

    /**
     * Scores {@code rows} (the aligned sequences, one string each) over {@code length} columns.
     *
     * @param nucleotide whether to normalise the information content against a 4-letter or a 20-letter alphabet
     */
    static MsaConservation compute( final List<String> rows, final int length, final boolean nucleotide ) {
        final int n_cols = Math.max( 0, length );
        final int n_rows = ( rows == null ) ? 0 : rows.size();
        final char[] consensus = new char[ n_cols ];
        final double[] identity = new double[ n_cols ];
        final double[] information = new double[ n_cols ];
        if ( ( n_rows == 0 ) || ( n_cols == 0 ) ) {
            return new MsaConservation( consensus, identity, information, n_rows );
        }
        final int k = nucleotide ? K_NUCLEOTIDE : K_AMINO_ACID;
        final double log_k = Math.log( k );
        final int[] counts = new int[ 128 ]; // ASCII is enough: residues are letters
        for( int col = 0; col < n_cols; col++ ) {
            java.util.Arrays.fill( counts, 0 );
            int non_gap = 0;
            int wide = 0; // residues outside ASCII -- counted as present, but each as its own symbol
            for( final String row : rows ) {
                if ( ( row == null ) || ( col >= row.length() ) ) {
                    continue; // past this row's end: a gap (see the class comment)
                }
                final char c = Character.toUpperCase( row.charAt( col ) );
                if ( MsaColors.isGap( c ) ) {
                    continue;
                }
                non_gap++;
                if ( c < counts.length ) {
                    counts[ c ]++;
                }
                else {
                    wide++;
                }
            }
            if ( non_gap == 0 ) {
                continue; // consensus 0, both scores 0
            }
            int best_count = 0;
            char best = 0;
            for( char c = 0; c < counts.length; c++ ) {
                if ( counts[ c ] > best_count ) { // strictly greater -> the alphabetically first of a tie wins
                    best_count = counts[ c ];
                    best = c;
                }
            }
            consensus[ col ] = best;
            identity[ col ] = (double) best_count / n_rows;
            // Shannon entropy over the residues ACTUALLY OBSERVED (p sums to 1 over the non-gap rows), then
            // normalised by log2(K) and scaled by the non-gap fraction so gaps cost the same as they do above.
            double h = 0;
            for( char c = 0; c < counts.length; c++ ) {
                if ( counts[ c ] > 0 ) {
                    final double p = (double) counts[ c ] / non_gap;
                    h -= p * Math.log( p );
                }
            }
            if ( wide > 0 ) { // treat each non-ASCII residue as its own symbol, the least flattering reading
                final double p = 1.0 / non_gap;
                h -= wide * ( p * Math.log( p ) );
            }
            final double normalized = 1.0 - ( h / log_k );
            information[ col ] = clamp01( normalized ) * ( (double) non_gap / n_rows );
        }
        return new MsaConservation( consensus, identity, information, n_rows );
    }

    private static double clamp01( final double d ) {
        if ( Double.isNaN( d ) || ( d < 0 ) ) {
            return 0;
        }
        return ( d > 1 ) ? 1 : d;
    }

    /**
     * How the track names itself in the figure: the measure, and how many sequences it was scored over. The count
     * is the load-bearing half -- the profile covers the tips CURRENTLY DISPLAYED, so without it a reader cannot
     * tell whether a bar describes the whole alignment or the six tips left after a clade was collapsed.
     */
    static String label( final Measure measure, final int rows ) {
        return measure + " (n = " + rows + ")";
    }

    /** Number of columns scored. */
    int length() {
        return _identity.length;
    }

    /** Number of rows the scores were computed over -- the tips that were on screen. */
    int rows() {
        return _rows;
    }

    /** The column's most common non-gap residue (upper case), or 0 when the column holds no residue. */
    char consensusAt( final int col ) {
        return ( ( col < 0 ) || ( col >= _consensus.length ) ) ? 0 : _consensus[ col ];
    }

    /** The column's score under {@code measure}, in [0,1]. Out-of-range columns score 0. */
    double scoreAt( final int col, final Measure measure ) {
        if ( ( col < 0 ) || ( col >= _identity.length ) ) {
            return 0;
        }
        return ( measure == Measure.INFORMATION ) ? _information[ col ] : _identity[ col ];
    }

    double identityAt( final int col ) {
        return scoreAt( col, Measure.IDENTITY );
    }

    double informationAt( final int col ) {
        return scoreAt( col, Measure.INFORMATION );
    }
}
