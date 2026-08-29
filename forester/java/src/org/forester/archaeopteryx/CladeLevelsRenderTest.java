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

import java.awt.GraphicsEnvironment;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.TreePanel.CLADE_LABEL_ANGLE;
import org.forester.archaeopteryx.TreePanel.CLADE_VIS;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * "Annotate Clades by Rank" at up to three ranks at once, drawn as NESTED bar/bracket columns. Runs on the bat
 * demo, whose internal nodes carry family / subfamily / genus offline, and asserts the geometry contract the
 * nesting rests on:
 * <ul>
 * <li>the levels come out finest-rank-first, whatever order they were asked for;</li>
 * <li>each level gets its own column, marching strictly OUTWARD, with enough room between columns for the mark
 * and its labels -- overlapping columns would draw one rank's bars through another's names;</li>
 * <li>adding outer levels does not move the innermost column, so the finest rank sits in the same place
 * regardless of how many broader ranks are drawn beside it;</li>
 * <li>a wider label angle really does claim more room (the reserve is angle-aware, which is what lets three
 * levels fit at all);</li>
 * <li>BOXES stay single-level -- nesting alpha-composited washes would double-darken the inner clade;</li>
 * <li>and the marks actually paint, in a rectangular AND in the circular layout.</li>
 * </ul>
 * Guarded to a no-op on a headless box.
 */
public final class CladeLevelsRenderTest {

    private static final String DEMO = "forester/demo/bat-phylogeny.xml";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CladeLevelsRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final File file = new File( System.getProperty( "user.dir" ), DEMO );
            if ( !file.exists() ) {
                System.out.println( "  [CladeLevelsRenderTest] missing demo tree " + file.getAbsolutePath() );
                return false;
            }
            final Phylogeny[] phys = ParserBasedPhylogenyFactory.getInstance()
                    .create( file, PhyloXmlParser.createPhyloXmlParser() );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( phys, new Configuration(), "clade-levels" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    exercise( mf[ 0 ], ok );
                }
                finally {
                    ( ( JFrame ) mf[ 0 ] ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void exercise( final MainFrame frame, final boolean[] ok ) {
        frame.setSize( 1200, 900 );
        ( ( JFrame ) frame ).validate();
        final TreePanel tp = frame.getCurrentTreePanel();
        frame.getOptions().setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );

        // ---- one level: the baseline the nested versions must not disturb -------------------------------
        final int one = tp.setCladeLevels( specs( "genus" ), CLADE_VIS.BARS, true );
        if ( one <= 0 ) {
            fail( ok, "the bat demo should place genus-level clades offline, got " + one + " marks" );
            return;
        }
        final float[] one_start = tp.cladeLevelStartsForTest();
        if ( one_start.length != 1 ) {
            fail( ok, "one rank must give one column, got " + one_start.length );
            return;
        }

        // ---- three levels, asked for in the WRONG order: they must come back finest-first ---------------
        final int three = tp.setCladeLevels(
                specs( "family", CLADE_LABEL_ANGLE.VERTICAL, "genus", CLADE_LABEL_ANGLE.VERTICAL, "subfamily",
                       CLADE_LABEL_ANGLE.VERTICAL ),
                CLADE_VIS.BARS, true );
        if ( !tp.cladeLevelRanks().equals( Arrays.asList( "genus", "subfamily", "family" ) ) ) {
            fail( ok, "levels must draw finest-first regardless of pick order, got " + tp.cladeLevelRanks() );
        }
        if ( three <= one ) {
            fail( ok, "three levels must draw more marks than one, got " + three + " vs " + one );
        }
        final float[] starts = tp.cladeLevelStartsForTest();
        if ( starts.length != 3 ) {
            fail( ok, "three ranks must give three columns, got " + starts.length );
            return;
        }
        // strictly outward, and far enough apart that a column's mark + labels clear the next column
        final int label_h = tp.getLargeFontHeight(); // a VERTICAL level's label occupies exactly one line height
        for ( int i = 1; i < starts.length; ++i ) {
            final float gap = starts[ i ] - starts[ i - 1 ];
            if ( gap <= 0 ) {
                fail( ok, "clade column " + i + " must sit OUTSIDE column " + ( i - 1 ) + " (starts "
                        + Arrays.toString( starts ) + ")" );
            }
            else if ( gap < label_h ) {
                fail( ok, "clade columns " + ( i - 1 ) + " and " + i + " are only " + gap
                        + " px apart -- not room for a mark plus its labels (line height " + label_h + ")" );
            }
        }
        // the innermost column is anchored to the tips, so adding outer levels must not shift it
        if ( Math.abs( starts[ 0 ] - one_start[ 0 ] ) > 0.5f ) {
            fail( ok, "adding outer levels moved the innermost column: " + one_start[ 0 ] + " -> " + starts[ 0 ] );
        }

        // ---- a wider label angle must claim more room, or three levels could never be made to fit -------
        tp.setCladeLevels( specs( "genus", CLADE_LABEL_ANGLE.HORIZONTAL, "family", CLADE_LABEL_ANGLE.VERTICAL ),
                           CLADE_VIS.BARS, true );
        final float wide_gap = tp.cladeLevelStartsForTest()[ 1 ] - tp.cladeLevelStartsForTest()[ 0 ];
        tp.setCladeLevels( specs( "genus", CLADE_LABEL_ANGLE.VERTICAL, "family", CLADE_LABEL_ANGLE.VERTICAL ),
                           CLADE_VIS.BARS, true );
        final float narrow_gap = tp.cladeLevelStartsForTest()[ 1 ] - tp.cladeLevelStartsForTest()[ 0 ];
        if ( wide_gap <= narrow_gap ) {
            fail( ok, "a HORIZONTAL level must reserve more depth than a VERTICAL one, got " + wide_gap + " vs "
                    + narrow_gap );
        }

        // ---- BOXES stay single-level (nested translucent washes would double-darken the inner clade) ----
        tp.setCladeLevels( specs( "genus", "family", "subfamily" ), CLADE_VIS.BOXES, true );
        if ( tp.cladeLevelStartsForTest().length != 1 ) {
            fail( ok, "BOXES must draw a single level, got " + tp.cladeLevelStartsForTest().length );
        }

        // ---- a rank that places NOTHING must not claim a column, nor be reported as drawn ---------------
        // ("mystery" is not a rank any tree resolves, so its level comes back empty)
        tp.setCladeLevels( specs( "genus", "mystery" ), CLADE_VIS.BARS, true );
        if ( tp.cladeLevelRanks().contains( "mystery" ) ) {
            fail( ok, "a rank that placed no clade must not be reported as drawn, got " + tp.cladeLevelRanks() );
        }
        final float[] one_real = tp.cladeLevelStartsForTest();
        tp.setCladeLevels( specs( "genus" ), CLADE_VIS.BARS, true );
        if ( one_real.length != tp.cladeLevelStartsForTest().length ) {
            fail( ok, "an empty level must not reserve a column: " + one_real.length + " vs "
                    + tp.cladeLevelStartsForTest().length );
        }

        // ---- the vertical orientations: an upright label must clear the bar it names ---------------------
        verticalLabelsClearTheirBands( frame, ok );
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );

        // ---- and it all actually paints, rectangular AND circular --------------------------------------
        for ( final CLADE_VIS mode : new CLADE_VIS[] { CLADE_VIS.BARS, CLADE_VIS.BRACKETS } ) {
            for ( final Options.PHYLOGENY_GRAPHICS_TYPE gt : new Options.PHYLOGENY_GRAPHICS_TYPE[] {
                    Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR } ) {
                tp.setPhylogenyGraphicsType( gt );
                tp.setCladeLevels( specs( "genus", "subfamily", "family" ), mode, true );
                try {
                    if ( AptxUtil.renderPhylogenyToImage( 900, 900, tp, frame.getOptions(), false, 1,
                                                          false ) == null ) {
                        fail( ok, "rendering three " + mode + " levels in " + gt + " produced no image" );
                    }
                }
                catch ( final Exception e ) {
                    fail( ok, "rendering three " + mode + " levels in " + gt + " threw " + e );
                }
            }
        }
    }

    /**
     * In a root-top / root-bottom layout a clade label is drawn UPRIGHT, re-anchored to the un-rotated frame at the
     * point just past its mark. It used to be CENTRED on that anchor, which left the lower half of every name
     * printed inside the coloured bar it labels ("Ctenophora" sitting halfway into the blue band). The whole glyph
     * box must fall on the far side of the anchor instead: below it with the root at the top, above it with the
     * root at the bottom -- the same rule an upright tip label follows.
     */
    private static void verticalLabelsClearTheirBands( final MainFrame frame, final boolean[] ok ) {
        final int ascent = 11, descent = 3; // representative metrics; the rule must hold for any
        for ( final double anchor : new double[] { 0, 37.5, 400 } ) {
            // root at the BOTTOM: the label sits ABOVE the anchor, so its lowest pixel may not pass it
            final double up = TreePanel.uprightCladeLabelBaseline( anchor, ascent, descent, true );
            if ( ( up + descent ) > anchor ) {
                fail( ok, "root-bottom: the label's bottom (" + ( up + descent ) + ") must not cross the mark at "
                        + anchor );
            }
            // root at the TOP: the label sits BELOW the anchor, so its highest pixel may not pass it
            final double down = TreePanel.uprightCladeLabelBaseline( anchor, ascent, descent, false );
            if ( ( down - ascent ) < anchor ) {
                fail( ok, "root-top: the label's top (" + ( down - ascent ) + ") must not cross the mark at "
                        + anchor );
            }
            // and the two must fall on OPPOSITE sides -- one rule flipped by the orientation, not one offset
            if ( up >= down ) {
                fail( ok, "the root-bottom baseline must sit above the root-top baseline at anchor " + anchor );
            }
        }
        // it also has to survive a real render in both vertical orientations
        final TreePanel tp = frame.getCurrentTreePanel();
        for ( final Options.TREE_ORIENTATION o : new Options.TREE_ORIENTATION[] {
                Options.TREE_ORIENTATION.ROOT_BOTTOM, Options.TREE_ORIENTATION.ROOT_TOP } ) {
            frame.getOptions().setTreeOrientation( o );
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp.setCladeLevels( specs( "phylum", "class" ), CLADE_VIS.BARS, false );
            try {
                if ( AptxUtil.renderPhylogenyToImage( 1100, 900, tp, frame.getOptions(), false, 1, false ) == null ) {
                    fail( ok, "nested clade bars produced no image in " + o );
                }
            }
            catch ( final Exception e ) {
                fail( ok, "nested clade bars in " + o + " threw " + e );
            }
        }
        frame.getOptions().setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT ); // leave it as we found it
    }

    private static List<CladeLevel.Spec> specs( final String... ranks ) {
        final List<CladeLevel.Spec> out = new ArrayList<>();
        for ( final String r : ranks ) {
            out.add( new CladeLevel.Spec( r, CLADE_LABEL_ANGLE.VERTICAL ) );
        }
        return out;
    }

    private static List<CladeLevel.Spec> specs( final Object... rank_then_angle ) {
        final List<CladeLevel.Spec> out = new ArrayList<>();
        for ( int i = 0; i < rank_then_angle.length; i += 2 ) {
            out.add( new CladeLevel.Spec( ( String ) rank_then_angle[ i ],
                                          ( CLADE_LABEL_ANGLE ) rank_then_angle[ i + 1 ] ) );
        }
        return out;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [CladeLevelsRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private CladeLevelsRenderTest() {
    }
}
