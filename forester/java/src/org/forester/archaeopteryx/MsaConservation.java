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
 * Two measures, both on <b>[0,1]</b> so either can drive the same bar -- and a third choice,
 * {@link Measure#LOGO}, which is the second of them drawn as letters rather than as a bar:
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
 * <li>A residue outside ASCII (malformed input) is counted BY CHARACTER, exactly like any other residue -- two
 * occurrences of the same one are two of the same symbol. Until 2026-09-25 each occurrence was counted as its own
 * distinct symbol, described here as "the least flattering reading". It was simply wrong, and measurably so: a
 * three-row column reading A plus two identical non-ASCII residues scored 2.736966 bits where it holds 3.403632,
 * and the logo's letter stack summed to a third of the height the bar drew for the same column.</li>
 * </ul>
 */
final class MsaConservation {

    /** What the track shows. {@code toString} is what the Settings dropdown displays. */
    enum Measure {
        IDENTITY( "Consensus identity" ),
        INFORMATION( "Information content" ),
        /** The information content drawn as a SEQUENCE LOGO -- a stack of letters per column, each letter's height
         *  its share of the column's information. Scores the same number as {@link #INFORMATION}; only the drawing
         *  differs, which is why it lives in this enum rather than in a second on/off switch the user would have to
         *  reconcile with this one. */
        LOGO( "Sequence logo" );

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

    private final char[]     _consensus;   // 0 where the column holds no residue at all
    private final double[]   _identity;
    private final double[]   _information;
    private final char[][]   _stack_residues;  // per column, residues ordered as the logo stacks them
    private final double[][] _stack_fractions; // per column, each residue's share of the BAND (sums to _information)
    private final int        _rows;

    private MsaConservation( final char[] consensus, final double[] identity, final double[] information,
                             final char[][] stack_residues, final double[][] stack_fractions, final int rows ) {
        _consensus = consensus;
        _identity = identity;
        _information = information;
        _stack_residues = stack_residues;
        _stack_fractions = stack_fractions;
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
        final char[][] stack_residues = new char[ n_cols ][];
        final double[][] stack_fractions = new double[ n_cols ][];
        if ( ( n_rows == 0 ) || ( n_cols == 0 ) ) {
            return new MsaConservation( consensus, identity, information, stack_residues, stack_fractions, n_rows );
        }
        final int k = nucleotide ? K_NUCLEOTIDE : K_AMINO_ACID;
        final double log_k = Math.log( k );
        final int[] counts = new int[ 128 ]; // the fast path: a residue is a letter, and letters are ASCII
        for( int col = 0; col < n_cols; col++ ) {
            java.util.Arrays.fill( counts, 0 );
            int non_gap = 0;
            // Residues outside ASCII (malformed input). Counted BY CHARACTER like any other residue, in a map
            // allocated only when one actually appears, so a normal alignment pays nothing. A TreeMap because its
            // ascending key order is the alphabetical tie-break the stack needs, for free.
            java.util.TreeMap<Character, Integer> wide = null;
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
                    if ( wide == null ) {
                        wide = new java.util.TreeMap<Character, Integer>();
                    }
                    final Integer prev = wide.get( Character.valueOf( c ) );
                    wide.put( Character.valueOf( c ),
                              Integer.valueOf( ( prev == null ) ? 1 : ( prev.intValue() + 1 ) ) );
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
            if ( wide != null ) { // ...and the same rule across the two, ASCII sorting first
                for( final java.util.Map.Entry<Character, Integer> e : wide.entrySet() ) {
                    if ( e.getValue().intValue() > best_count ) {
                        best_count = e.getValue().intValue();
                        best = e.getKey().charValue();
                    }
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
            if ( wide != null ) {
                for( final Integer n : wide.values() ) {
                    final double p = n.intValue() / (double) non_gap;
                    h -= p * Math.log( p );
                }
            }
            final double normalized = 1.0 - ( h / log_k );
            information[ col ] = clamp01( normalized ) * ( (double) non_gap / n_rows );
            buildStack( counts, wide, non_gap, information[ col ], col, stack_residues, stack_fractions );
        }
        return new MsaConservation( consensus, identity, information, stack_residues, stack_fractions, n_rows );
    }

    /**
     * The logo stack for one column: every residue present, ordered by count DESCENDING with alphabetical order
     * within a tie, each carrying its share of the column's band height ({@code p * information}).
     * <p>
     * The order is part of the joint contract with archaeopteryx.js -- four residues at equal count must come out
     * in the same order in both viewers, or the same alignment reads differently in each. It falls out of a STABLE
     * sort over {@code counts}, which is already in ascending character order.
     * <p>
     * Built from the ASCII counts, so a residue outside ASCII (malformed input; see the class comment) contributes
     * to the column's height through the entropy above but gets no letter of its own. The stack is then SHORTER
     * than the column scored -- which is the honest drawing, since there is no residue to name.
     */
    private static void buildStack( final int[] counts, final java.util.TreeMap<Character, Integer> wide,
                                    final int non_gap, final double column_height, final int col,
                                    final char[][] stack_residues, final double[][] stack_fractions ) {
        int distinct = ( wide == null ) ? 0 : wide.size();
        for( int c = 0; c < counts.length; c++ ) {
            if ( counts[ c ] > 0 ) {
                distinct++;
            }
        }
        if ( distinct == 0 ) {
            return;
        }
        final char[] residues = new char[ distinct ];
        final int[] found = new int[ distinct ];
        int at = 0;
        for( char c = 0; c < counts.length; c++ ) {
            if ( counts[ c ] > 0 ) {
                residues[ at ] = c;
                found[ at ] = counts[ c ];
                at++;
            }
        }
        if ( wide != null ) { // after the ASCII residues, themselves in ascending character order
            for( final java.util.Map.Entry<Character, Integer> e : wide.entrySet() ) {
                residues[ at ] = e.getKey().charValue();
                found[ at ] = e.getValue().intValue();
                at++;
            }
        }
        // insertion sort by count descending -- STABLE, so the ascending-character order above survives a tie
        for( int i = 1; i < distinct; i++ ) {
            final char r = residues[ i ];
            final int n = found[ i ];
            int j = i - 1;
            while( ( j >= 0 ) && ( found[ j ] < n ) ) {
                residues[ j + 1 ] = residues[ j ];
                found[ j + 1 ] = found[ j ];
                j--;
            }
            residues[ j + 1 ] = r;
            found[ j + 1 ] = n;
        }
        final double[] fractions = new double[ distinct ];
        for( int i = 0; i < distinct; i++ ) {
            fractions[ i ] = ( (double) found[ i ] / non_gap ) * column_height;
        }
        stack_residues[ col ] = residues;
        stack_fractions[ col ] = fractions;
    }

    /** The residues of {@code col}'s logo stack, most frequent FIRST (the drawing stacks from the bottom up, so
     *  this one ends on top). Empty when the column holds no residue at all. */
    char[] stackResiduesAt( final int col ) {
        final char[] r = ( ( col < 0 ) || ( col >= _stack_residues.length ) ) ? null : _stack_residues[ col ];
        return ( r == null ) ? NO_RESIDUES : r;
    }

    /** Each stacked residue's share of the band, in the same order as {@link #stackResiduesAt}. The sum is the
     *  column's information content, so a gappy or variable column stacks shorter. */
    double[] stackFractionsAt( final int col ) {
        final double[] f = ( ( col < 0 ) || ( col >= _stack_fractions.length ) ) ? null : _stack_fractions[ col ];
        return ( f == null ) ? NO_FRACTIONS : f;
    }

    private static final char[]   NO_RESIDUES  = new char[ 0 ];
    private static final double[] NO_FRACTIONS = new double[ 0 ];

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
        // LOGO is INFORMATION drawn as letters, so it must score the identical number: were it to fall through to
        // identity here, the letter stack and the value the rest of the program reports for the same column would
        // disagree, and the stack heights would be built from one measure while the readout named another.
        return ( ( measure == Measure.INFORMATION ) || ( measure == Measure.LOGO ) ) ? _information[ col ]
                : _identity[ col ];
    }

    double identityAt( final int col ) {
        return scoreAt( col, Measure.IDENTITY );
    }

    double informationAt( final int col ) {
        return scoreAt( col, Measure.INFORMATION );
    }
}
