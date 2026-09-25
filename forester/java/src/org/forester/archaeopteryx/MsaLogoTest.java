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
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;

/**
 * The sequence logo drawn under the alignment: a stack of letters per column, each letter's height its share of the
 * column's information content.
 * <p>
 * <b>The numbers are a JOINT CONTRACT with archaeopteryx.js</b> and are checked here against the values that session
 * published from its own implementation, to six decimal places, in ITS units (bits) rather than ours. Both sides
 * must agree on the information formula, the deliberate ABSENCE of a small-sample correction, what counts as a gap,
 * how gaps scale a column, and -- the one that is invisible until two viewers are put side by side -- the order
 * letters stack in when their counts TIE.
 * <p>
 * The drawing is checked against PIXELS, not against the code's own idea of where it drew. Two traps live there and
 * both are silent: scaling a letter by the font's LINE height instead of its own INK makes every letter slightly
 * wrong, and putting the baseline at the letter's bottom instead of one scaled DESCENT above it makes every G, Q and
 * J hang through the ruler underneath. A test that asked the layout where the letters went would agree with both.
 */
public final class MsaLogoTest {

    private final static int    W        = 1400;
    private final static int    H        = 900;
    private final static double TOL      = 1e-6;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MsaLogo: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        try {
            jointNumbers( ok );
            gapsAndOrder( ok );
            glyphInkBox( ok );
            if ( !GraphicsEnvironment.isHeadless() ) {
                drawnStack( ok );
            }
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
        return ok[ 0 ];
    }

    // ---- the joint numbers, in archaeopteryx.js's units --------------------------------------------------------

    /**
     * Every case the archaeopteryx.js session published from its own implementation. Ours stores each letter's share
     * of the BAND (a fraction of the maximum information), theirs stores it in BITS, so the comparison multiplies by
     * {@code log2(K)} -- which is itself part of what is being checked, since a wrong alphabet size would show up
     * here as every protein number being off by the same factor.
     */
    private static void jointNumbers( final boolean[] ok ) {
        final double max_nuc = Math.log( 4 ) / Math.log( 2 );   // 2
        final double max_aa = Math.log( 20 ) / Math.log( 2 );   // 4.321928
        // rows = ACGT-A / ACGT-A / ACGA-C / TCGA-G, nucleotide
        final MsaConservation main = MsaConservation
                .compute( Arrays.asList( "ACGT-A", "ACGT-A", "ACGA-C", "TCGA-G" ), 6, true );
        bits( ok, main, 0, max_nuc, "col1", new char[] { 'A', 'T' }, new double[] { 0.891541, 0.297180 }, 1.188722 );
        bits( ok, main, 1, max_nuc, "col2", new char[] { 'C' }, new double[] { 2.000000 }, 2.000000 );
        bits( ok, main, 2, max_nuc, "col3", new char[] { 'G' }, new double[] { 2.000000 }, 2.000000 );
        bits( ok, main, 3, max_nuc, "col4", new char[] { 'A', 'T' }, new double[] { 0.500000, 0.500000 }, 1.000000 );
        bits( ok, main, 4, max_nuc, "col5 (all gaps)", new char[] {}, new double[] {}, 0.000000 );
        bits( ok, main, 5, max_nuc, "col6", new char[] { 'A', 'C', 'G' },
              new double[] { 0.250000, 0.125000, 0.125000 }, 0.500000 );
        // partial gaps: full conservation at HALF height -- the occupancy rule, and the case that separates
        // "frequencies over non-gap rows" from "frequencies over all rows"
        bits( ok, MsaConservation.compute( Arrays.asList( "A", "A", "-", "-" ), 1, true ), 0, max_nuc,
              "partial gaps", new char[] { 'A' }, new double[] { 1.000000 }, 1.000000 );
        // a single sequence: H = 0, so full information. NOT special-cased on either side, and worth pinning
        // because a small-sample correction -- which we deliberately do not apply -- would blank it.
        final MsaConservation one = MsaConservation.compute( Arrays.asList( "AC" ), 2, true );
        bits( ok, one, 0, max_nuc, "single seq col1", new char[] { 'A' }, new double[] { 2.000000 }, 2.000000 );
        bits( ok, one, 1, max_nuc, "single seq col2", new char[] { 'C' }, new double[] { 2.000000 }, 2.000000 );
        // protein: the maximum is log2(20), so a conserved column is 4.321928 bits and never 2
        bits( ok, MsaConservation.compute( Arrays.asList( "L", "L", "L" ), 1, false ), 0, max_aa, "protein LLL",
              new char[] { 'L' }, new double[] { 4.321928 }, 4.321928 );
        // A residue outside ASCII is malformed input, but it is still ONE symbol however often it occurs. Counting
        // each occurrence separately -- which this did until 2026-09-25 -- under-reports the column AND leaves the
        // logo's stack summing to less than the height the bar draws beside it. Numbers from the JS session, run
        // against their implementation. The escape keeps this source file ASCII.
        final char omega = '\u03A9';
        final MsaConservation wide = MsaConservation
                .compute( Arrays.asList( "AA" + omega, "AA" + omega, "AAA" ), 3, false );
        bits( ok, wide, 2, max_aa, "non-ASCII residue", new char[] { omega, 'A' },
              new double[] { 2.269088, 1.134544 }, 3.403632 );
        double wide_sum = 0;
        for( final double f : wide.stackFractionsAt( 2 ) ) {
            wide_sum += f;
        }
        if ( Math.abs( wide_sum - wide.informationAt( 2 ) ) > TOL ) {
            fail( ok, "the stack must sum to the column height even for a residue outside ASCII: stack "
                    + wide_sum + " against column " + wide.informationAt( 2 ) );
        }
        // THE TIE-BREAK, made visible: four residues at equal count come out alphabetical. Every height is zero, so
        // only the ORDER is under test -- and the order is what makes two viewers stack the same column alike.
        bits( ok, MsaConservation.compute( Arrays.asList( "G", "C", "A", "T" ), 1, true ), 0, max_nuc, "4-way tie",
              new char[] { 'A', 'C', 'G', 'T' }, new double[] { 0, 0, 0, 0 }, 0.000000 );
    }

    /** Compares one column against the published residues and their heights IN BITS. */
    private static void bits( final boolean[] ok, final MsaConservation c, final int col, final double max_bits,
                              final String what, final char[] want_res, final double[] want_bits,
                              final double want_height ) {
        final char[] got_res = c.stackResiduesAt( col );
        final double[] got_f = c.stackFractionsAt( col );
        if ( !Arrays.equals( got_res, want_res ) ) {
            fail( ok, what + ": stack must be " + Arrays.toString( want_res ) + ", got "
                    + Arrays.toString( got_res ) );
            return;
        }
        if ( Math.abs( ( c.informationAt( col ) * max_bits ) - want_height ) > 1e-5 ) {
            fail( ok, what + ": column height must be " + want_height + " bits, got "
                    + ( c.informationAt( col ) * max_bits ) );
        }
        for( int i = 0; i < want_bits.length; i++ ) {
            if ( Math.abs( ( got_f[ i ] * max_bits ) - want_bits[ i ] ) > 1e-5 ) {
                fail( ok, what + ": " + want_res[ i ] + " must be " + want_bits[ i ] + " bits, got "
                        + ( got_f[ i ] * max_bits ) );
            }
        }
    }

    // ---- what is a gap, and what the stack sums to --------------------------------------------------------------

    private static void gapsAndOrder( final boolean[] ok ) {
        // ? and X are UNKNOWN RESIDUES, not gaps -- they take part in the stack and cost the column information.
        // Reading them as gaps would make a column of them look perfectly conserved.
        final MsaConservation q = MsaConservation.compute( Arrays.asList( "A", "A", "?", "X" ), 1, false );
        final char[] res = q.stackResiduesAt( 0 );
        if ( ( res.length != 3 ) || ( res[ 0 ] != 'A' ) ) {
            fail( ok, "? and X must stack as residues alongside A, got " + Arrays.toString( res ) );
        }
        if ( !( q.informationAt( 0 ) < 1.0 ) ) {
            fail( ok, "a column of A A ? X is not fully conserved, yet it scored " + q.informationAt( 0 ) );
        }
        // every gap character, on the other hand, is absent from the stack entirely
        final MsaConservation g = MsaConservation.compute( Arrays.asList( "A", "-", ".", "~", " " ), 1, false );
        if ( !Arrays.equals( g.stackResiduesAt( 0 ), new char[] { 'A' } ) ) {
            fail( ok, "- . ~ and space are gaps and must not stack, got "
                    + Arrays.toString( g.stackResiduesAt( 0 ) ) );
        }
        // Two residues that are the SAME SYMBOL must stack as one letter. Case first: the alignment is drawn and
        // coloured case-insensitively, so it must be stacked that way too.
        final MsaConservation mixed = MsaConservation.compute( Arrays.asList( "a", "A", "a", "A" ), 1, false );
        if ( !Arrays.equals( mixed.stackResiduesAt( 0 ), new char[] { 'A' } ) ) {
            fail( ok, "lower and upper case are one residue and must stack as one letter, got "
                    + Arrays.toString( mixed.stackResiduesAt( 0 ) ) );
        }
        // THE STACK IS THE COLUMN: its parts must add up to the height the bar draws beside them, on every
        // alignment, not just a tidy one. The earlier version of this ran on a single clean four-row fixture --
        // which is exactly why a column holding a residue outside ASCII could sum to a THIRD of its own height and
        // still pass. The cases below are the ones that broke it or could: a row shorter than the alignment, an
        // all-gap column, a repeated non-ASCII residue, mixed case, and a single row.
        final char omega = '\u03A9';
        final List<List<String>> corpus = new ArrayList<List<String>>();
        corpus.add( Arrays.asList( "AC-A", "AG-A", "AT-C", "CT-G" ) );
        corpus.add( Arrays.asList( "ACGT", "AC", "ACG", "A" ) );                    // ragged: short rows are gaps
        corpus.add( Arrays.asList( "A-" + omega + "A", "a-" + omega + "C", "A-A-" ) ); // non-ASCII, case, gaps
        corpus.add( Arrays.asList( "----" ) );                                       // nothing but gaps
        corpus.add( Arrays.asList( "ACGT" ) );                                       // one row
        for( int f = 0; f < corpus.size(); f++ ) {
            final MsaConservation m = MsaConservation.compute( corpus.get( f ), 4, false );
            for( int col = 0; col < 4; col++ ) {
                double sum = 0;
                for( final double frac : m.stackFractionsAt( col ) ) {
                    sum += frac;
                }
                if ( Math.abs( sum - m.informationAt( col ) ) > TOL ) {
                    fail( ok, "alignment " + f + " column " + col + ": the stack sums to " + sum
                            + " but the column scores " + m.informationAt( col )
                            + " -- the letters would be drawn a different height from the bar for the same data" );
                }
                // ...and no residue may be stacked twice: a symbol counted under two identities inflates the
                // letter count while the heights still add up, which the sum alone cannot see.
                final char[] stacked = m.stackResiduesAt( col );
                for( int i = 1; i < stacked.length; i++ ) {
                    if ( stacked[ i ] == stacked[ i - 1 ] ) {
                        fail( ok, "alignment " + f + " column " + col + ": " + stacked[ i ]
                                + " is stacked twice -- one symbol, one letter" );
                    }
                }
                // LOGO is INFORMATION drawn as letters, so it must report the identical score
                if ( Math.abs( m.scoreAt( col, MsaConservation.Measure.LOGO )
                        - m.scoreAt( col, MsaConservation.Measure.INFORMATION ) ) > TOL ) {
                    fail( ok, "alignment " + f + " column " + col
                            + ": LOGO and INFORMATION must score the same number" );
                }
            }
        }
    }

    // ---- the ink box, which is not the line box -----------------------------------------------------------------

    private static void glyphInkBox( final boolean[] ok ) {
        if ( GraphicsEnvironment.isHeadless() ) {
            return; // no font rasteriser to measure against
        }
        MsaLogoGlyphs.clearCacheForTest();
        final Font f = new Font( Font.MONOSPACED, Font.BOLD, 24 );
        final float[] cap_a = MsaLogoGlyphs.metrics( f, 'A' );
        final float[] cap_q = MsaLogoGlyphs.metrics( f, 'Q' );
        if ( ( cap_a[ 1 ] <= 0 ) || ( cap_q[ 1 ] <= 0 ) ) {
            fail( ok, "A and Q must report ink above the baseline, got " + cap_a[ 1 ] + " / " + cap_q[ 1 ] );
            return;
        }
        // THE distinction. A capital A sits ON the baseline and has no ink below it; Q's tail does. The font's own
        // descent is the same number for both, which is exactly why it cannot be used here.
        if ( cap_a[ 2 ] > 1.0f ) {
            fail( ok, "a capital A has no descender, yet its ink descent is " + cap_a[ 2 ]
                    + " -- that is the LINE box, not the ink" );
        }
        if ( !( cap_q[ 2 ] > 1.0f ) ) {
            fail( ok, "Q's tail falls below the baseline, yet its ink descent is " + cap_q[ 2 ] );
        }
        // a space paints nothing, so it must report no ink at all rather than a line-box height
        if ( MsaLogoGlyphs.inkHeight( f, ' ' ) != 0 ) {
            fail( ok, "a space paints nothing, yet it reports " + MsaLogoGlyphs.inkHeight( f, ' ' ) + " of ink" );
        }
        // the cache must hand back the same measurement, not re-derive a different one
        if ( !Arrays.equals( cap_q, MsaLogoGlyphs.metrics( f, 'Q' ) ) ) {
            fail( ok, "a second measurement of the same glyph disagreed with the first" );
        }
    }

    // ---- the drawn stack, read off the pixels --------------------------------------------------------------------

    /**
     * Renders a designed alignment and measures the INK in two columns that are each 100% conserved -- one on a
     * letter with a descender ({@code Q}) and one without ({@code A}).
     * <p>
     * Both are the same height by construction, so whatever is true of one must be true of the other. That makes the
     * descender trap catchable without knowing a single font metric: if the baseline is placed at the bottom of the
     * ink instead of one scaled descent above it, Q's tail drops below A's and below the band.
     */
    private static void drawnStack( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        final Phylogeny tree = logoTree();
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { tree }, new Configuration(), "msa-logo" ) );
        final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
        final BufferedImage[] img = new BufferedImage[ 1 ];
        final int[][] band = new int[ 1 ][];
        SwingUtilities.invokeAndWait( () -> {
            tp.getOptions().setShowOverview( false );
            tp.setOvOn( false );
            tp.setShowMsa( true );
            tp.getOptions().setShowMsaConservation( true );
            tp.getOptions().setMsaColumnWidth( 26 );
            tp.getOptions().setMsaConservationMeasure( MsaConservation.Measure.LOGO );
            tp.setSize( W, H );
            tp.calcParametersForPainting( W, H );
            final BufferedImage i = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
            final Graphics2D g = i.createGraphics();
            g.setColor( tp.getTreeColorSet().getBackgroundColor() );
            g.fillRect( 0, 0, W, H );
            // EXPORT geometry: on screen the band floats at the VIEWPORT bottom, and an offscreen panel's visible
            // rect is not its size, so a screen render would place the band somewhere this test cannot predict.
            tp.paintPhylogeny( g, false, true, W, H, 0, 0 );
            g.dispose();
            img[ 0 ] = i;
            band[ 0 ] = tp.msaLogoBandForTest();
        } );
        if ( band[ 0 ] == null ) {
            fail( ok, "no logo band was drawn at all" );
            dispose( mf );
            return;
        }
        final int origin_x = band[ 0 ][ 0 ];
        final int baseline = band[ 0 ][ 1 ];
        final int band_h = band[ 0 ][ 2 ];
        final int cw = band[ 0 ][ 3 ];
        final int bg = tp.getTreeColorSet().getBackgroundColor().getRGB();
        // column 0 is all A, column 1 is all Q -- both fully conserved, so both fill the band
        final int[] a_ink = inkExtent( img[ 0 ], bg, origin_x, cw, baseline, band_h );
        final int[] q_ink = inkExtent( img[ 0 ], bg + 0, origin_x + cw, cw, baseline, band_h );
        if ( ( a_ink == null ) || ( q_ink == null ) ) {
            fail( ok, "a fully conserved column drew no ink in the band (A: " + Arrays.toString( a_ink ) + ", Q: "
                    + Arrays.toString( q_ink ) + ")" );
            dispose( mf );
            return;
        }
        // a column at full information fills the band. Measured to just ABOVE the baseline (see inkExtent), so
        // one pixel short of the band is the expected answer, not a discrepancy.
        if ( Math.abs( ( a_ink[ 1 ] - a_ink[ 0 ] ) - band_h ) > 5 ) {
            fail( ok, "a fully conserved column must fill the " + band_h + " px band, its ink is "
                    + ( a_ink[ 1 ] - a_ink[ 0 ] ) + " px" );
        }
        // THE descender check, done on the TOP of the ink. Both columns are 100% conserved, so both stacks start at
        // the top of the band whatever letter they are made of. Place a descender by its baseline instead of by its
        // ink and the whole letter shifts DOWN by its scaled descent -- about 9 px here -- which shows up as Q
        // starting lower than A, and its tail leaving the band at the other end.
        if ( Math.abs( q_ink[ 0 ] - a_ink[ 0 ] ) > 2 ) {
            fail( ok, "Q's ink starts at y=" + q_ink[ 0 ] + " but A's at y=" + a_ink[ 0 ]
                    + " -- both columns are fully conserved, so a letter with a descender is being placed by its "
                    + "baseline instead of by its ink, and its tail hangs below the band" );
        }
        if ( a_ink[ 0 ] < ( ( baseline - band_h ) - 2 ) ) {
            fail( ok, "the stack reaches y=" + a_ink[ 0 ] + ", above the band's top at y=" + ( baseline - band_h ) );
        }
        // ...and a less conserved column is SHORTER IN PROPORTION to its information. Column 3 is an even A/C
        // split. Its expected height is deliberately NOT written here as a number: in a 20-letter alphabet an even
        // split still carries 1 - 1/log2(20) = 0.77 of the maximum, not the half that the nucleotide intuition
        // suggests, and hard-coding the wrong intuition is how the first version of this check failed a correct
        // drawing. The claim is the one that matters anyway -- that the PIXELS follow the model.
        final int[] half = inkExtent( img[ 0 ], bg, origin_x + ( 3 * cw ), cw, baseline, band_h );
        if ( half == null ) {
            fail( ok, "the evenly split column drew nothing" );
        }
        else {
            final double want_ratio = tp.msaLogoFractionsForTest( 3 )[ 0 ] + tp.msaLogoFractionsForTest( 3 )[ 1 ];
            final double got_ratio = ( half[ 1 ] - half[ 0 ] ) / (double) ( a_ink[ 1 ] - a_ink[ 0 ] );
            if ( !( want_ratio < 0.95 ) ) {
                fail( ok, "precondition: the split column must be less conserved than the full one, it scored "
                        + want_ratio );
            }
            if ( Math.abs( got_ratio - want_ratio ) > 0.06 ) {
                fail( ok, "stack heights must be proportional to the information: column 3 scores "
                        + String.format( "%.3f", Double.valueOf( want_ratio ) ) + " of a full column but was drawn "
                        + String.format( "%.3f", Double.valueOf( got_ratio ) ) + " of one ("
                        + ( half[ 1 ] - half[ 0 ] ) + " px against " + ( a_ink[ 1 ] - a_ink[ 0 ] ) + ")" );
            }
        }
        // an all-gap column draws nothing at all -- "no data" must not look like "no conservation"
        if ( inkExtent( img[ 0 ], bg, origin_x + ( 4 * cw ), cw, baseline, band_h ) != null ) {
            fail( ok, "an all-gap column must draw no letter" );
        }
        // THE ORDER, read off the pixels. Column 6 is 3 parts A to 1 part C, and the two residues take different
        // colours from the very palette the alignment cells use, so "which letter is on top" is answerable without
        // recognising a glyph: the majority residue's ink must sit ABOVE the minority's. Drawn the other way up, a
        // reader takes the bottom letter for the consensus.
        final double a_y = meanInkY( img[ 0 ], MsaColors.colorFor( 'A', false ), origin_x + ( 6 * cw ), cw, baseline,
                                     band_h );
        final double c_y = meanInkY( img[ 0 ], MsaColors.colorFor( 'C', false ), origin_x + ( 6 * cw ), cw, baseline,
                                     band_h );
        if ( Double.isNaN( a_y ) || Double.isNaN( c_y ) ) {
            fail( ok, "the 3:1 column must draw both residues (A found: " + !Double.isNaN( a_y ) + ", C found: "
                    + !Double.isNaN( c_y ) + ")" );
        }
        else if ( !( a_y < c_y ) ) {
            fail( ok, "the stack is upside down: A is 3 of 4 rows and C is 1, so A must be drawn ABOVE C, but A's "
                    + "ink centres at y=" + Math.round( a_y ) + " and C's at y=" + Math.round( c_y ) );
        }
        dispose( mf );
    }

    /**
     * {@code { topY, bottomY }} of the ink in one band column, or null when the column is blank.
     * <p>
     * Stops one pixel ABOVE the baseline on purpose. The band draws a baseline RULE across its whole width, and the
     * column ruler sits immediately under it -- include either and every column reports ink, blank ones included,
     * which is exactly how the first version of this test "passed" an all-gap column and measured every stack at
     * the same 74 px.
     */
    private static int[] inkExtent( final BufferedImage img, final int bg, final int x0, final int w,
                                    final int baseline, final int band_h ) {
        // ...and starts exactly AT the band top, not above it: the row naming the measure is drawn immediately
        // over the band, and a few pixels of slack upward pulled that label into every column's measurement.
        final int y_from = Math.max( 0, baseline - band_h );
        final int y_to = Math.min( img.getHeight() - 1, baseline - 1 );
        int top = Integer.MAX_VALUE;
        int bottom = Integer.MIN_VALUE;
        for( int y = y_from; y <= y_to; y++ ) {
            for( int x = Math.max( 0, x0 + 2 ); x < Math.min( img.getWidth(), ( x0 + w ) - 2 ); x++ ) {
                if ( differs( img.getRGB( x, y ), bg ) ) {
                    if ( y < top ) {
                        top = y;
                    }
                    if ( y > bottom ) {
                        bottom = y;
                    }
                }
            }
        }
        return ( top == Integer.MAX_VALUE ) ? null : new int[] { top, bottom };
    }

    /** The mean y of the pixels painted in {@code want} inside one band column, or NaN when there are none. */
    private static double meanInkY( final BufferedImage img, final Color want, final int x0, final int w,
                                    final int baseline, final int band_h ) {
        long sum = 0;
        long n = 0;
        for( int y = Math.max( 0, baseline - band_h ); y <= ( baseline - 1 ); y++ ) {
            for( int x = Math.max( 0, x0 + 2 ); x < Math.min( img.getWidth(), ( x0 + w ) - 2 ); x++ ) {
                final Color got = new Color( img.getRGB( x, y ) );
                // an exact match only: an antialiased edge pixel is a blend and belongs to neither letter
                if ( ( got.getRed() == want.getRed() ) && ( got.getGreen() == want.getGreen() )
                        && ( got.getBlue() == want.getBlue() ) ) {
                    sum += y;
                    n++;
                }
            }
        }
        return ( n == 0 ) ? Double.NaN : ( sum / (double) n );
    }

    private static boolean differs( final int rgb, final int bg ) {
        final Color a = new Color( rgb );
        final Color b = new Color( bg );
        return ( Math.abs( a.getRed() - b.getRed() ) + Math.abs( a.getGreen() - b.getGreen() )
                + Math.abs( a.getBlue() - b.getBlue() ) ) > 60;
    }

    /**
     * 24 tips over a designed 6-column alignment: col0 all A, col1 all Q (a descender), col2 all G, col3 an even
     * A/C split, col4 all gap, col5 all W. Repeating three row patterns keeps the frequencies exact while giving
     * enough tips that the rows are thin.
     */
    private static Phylogeny logoTree() throws org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException {
        // col0 all A, col1 all Q (a descender), col2 all G, col3 an even A/C split, col4 all gap, col5 all W,
        // col6 a 3:1 A/C split -- the last one exists so the stack has two letters of DIFFERENT heights, which is
        // what makes the ORDER visible. With only the even split, drawing the stack upside down changed nothing.
        final String[] patterns = { "AQGA-WA", "AQGA-WA", "AQGC-WA", "AQGC-WC" };
        final String[] names = { "alpha", "bravo", "cobra", "delta", "eagle", "fjord" };
        final PhylogenyNode root = new PhylogenyNode();
        final List<PhylogenyNode> tips = new ArrayList<PhylogenyNode>();
        for( int i = 0; i < 24; i++ ) {
            final PhylogenyNode t = new PhylogenyNode();
            t.setDistanceToParent( 0.1 );
            t.setName( names[ i % names.length ] + ( i / names.length ) );
            final Sequence s = new Sequence();
            s.setName( t.getName() );
            s.setMolecularSequence( patterns[ i % patterns.length ] );
            s.setMolecularSequenceAligned( true );
            s.setType( "protein" );
            t.getNodeData().addSequence( s );
            root.addAsChild( t );
            tips.add( t );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        phy.setName( "msa-logo" );
        return phy;
    }

    private static void dispose( final MainFrame[] mf ) throws Exception {
        SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [MsaLogoTest] " + message );
        ok[ 0 ] = false;
    }

    private MsaLogoTest() {
        // not instantiable
    }
}
