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
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
 * <li>shaded boxes go with up to two bar/bracket columns (Christian, 2026-11-04, "go, box + two"): the box level is
 * not a column, its names (beside the boxes) take their own slot before the first column, a column at the boxes'
 * own rank names each box, the pair wears one colour, the legend is keyed per level, the figure keeps it all and an
 * older figure still opens as it did;</li>
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
                    oneTipBoxesDraw( mf[ 0 ], ok );
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
        frame.getCurrentTreePanel().setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
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
        if ( ( tp.cladeLevelStartsForTest().length != 0 ) || !tp.cladeLevelRanks().equals( Arrays.asList( "genus" ) )
                || ( tp.drawnCladeColumnCount() != 0 ) || ( tp.drawnCladeBoxCount() <= 0 ) ) {
            fail( ok, "BOXES must draw a single level, the finest, and no column; got columns "
                    + tp.cladeLevelStartsForTest().length + ", ranks " + tp.cladeLevelRanks() );
        }

        // ---- shaded boxes WITH bar/bracket columns --------------------------------------------------------
        boxesWithColumns( frame, ok, one_start[ 0 ] );
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        frame.getCurrentTreePanel().setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );

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
            frame.getCurrentTreePanel().setTreeOrientation( o );
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
        frame.getCurrentTreePanel().setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT ); // leave it as we found it
    }

    /**
     * Shaded boxes plus columns. {@code bars_start}: where a lone bar column at genus starts -- the boxes, named in
     * the legend, take no column room, so a family column past genus boxes starts exactly there.
     */
    private static void boxesWithColumns( final MainFrame frame, final boolean[] ok, final float bars_start ) {
        final TreePanel tp = frame.getCurrentTreePanel();
        final int label_h = tp.getLargeFontHeight();
        tp.setCladeBoxLabels( false );
        // genus boxes inside family bars: two kinds of mark, two ranks
        final List<CladeLevel.Spec> nested = new ArrayList<>();
        nested.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BOXES ) );
        nested.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BARS ) );
        final int n = tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
        if ( !"genus".equals( tp.cladeBoxRank() ) || !tp.cladeColumnRanks().equals( Arrays.asList( "family" ) )
                || !tp.cladeLevelRanks().equals( Arrays.asList( "genus", "family" ) ) ) {
            fail( ok, "genus boxes + family bars: box rank " + tp.cladeBoxRank() + ", columns " + tp.cladeColumnRanks()
                    + ", ranks " + tp.cladeLevelRanks() );
            return;
        }
        if ( ( tp.drawnCladeBoxCount() <= 0 ) || ( tp.drawnCladeColumnCount() <= 0 )
                || ( n != ( tp.drawnCladeBoxCount() + tp.drawnCladeColumnCount() ) ) ) {
            fail( ok, "the drawn count is the boxes plus the column marks: " + n + " vs " + tp.drawnCladeBoxCount()
                    + " + " + tp.drawnCladeColumnCount() );
        }
        if ( tp.getCladeBandsMode() != CLADE_VIS.BARS ) {
            fail( ok, "the mode a figure records is the columns' style, got " + tp.getCladeBandsMode() );
        }
        // every box draws (the one-tip case, which the bat demo has no genus for, is oneTipBoxesDraw's)
        if ( tp.drawnCladeBoxCount() != tp.cladeLevelBandsForTest().get( 0 ).size() ) {
            fail( ok, "every box draws: drawn " + tp.drawnCladeBoxCount() + " of "
                    + tp.cladeLevelBandsForTest().get( 0 ).size() );
        }
        final int unnamed_reserve = tp.cladeBandRightReserveForTest();
        float[] starts = tp.cladeLevelStartsForTest();
        if ( starts.length != 1 ) {
            fail( ok, "boxes + one column must give one column, got " + starts.length );
            return;
        }
        if ( Math.abs( starts[ 0 ] - bars_start ) > 0.5f ) {
            fail( ok, "boxes named in the legend take no column room: the family column must start where a lone column"
                    + " starts (" + bars_start + "), got " + starts[ 0 ] );
        }
        // the legend keys both: the boxes' genera and the bars' families, in two sections
        final Map<String, Color> rows = tp.rankLegendColorsForTest();
        if ( !tp.hasRankLegend() || !rows.containsKey( "Myotis" ) || !rows.containsKey( "Molossidae" ) ) {
            fail( ok, "boxes in the legend + bars: a legend with a genus and a family row, got " + rows.keySet() );
        }
        // the boxes' names beside them take their own slot BEFORE the first column
        tp.setCladeBoxLabels( true );
        tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
        if ( !tp.cladeBoxNamesShownForTest() ) {
            fail( ok, "named beside each box, with a column at another rank: the names are drawn" );
        }
        final float[] named = tp.cladeLevelStartsForTest();
        if ( ( named.length != 1 ) || ( named[ 0 ] < ( starts[ 0 ] + label_h ) ) ) {
            fail( ok, "the column must move past the boxes' names: " + starts[ 0 ] + " -> "
                    + ( ( named.length == 1 ) ? named[ 0 ] : -1 ) + " (a line is " + label_h + ")" );
        }
        if ( !tp.hasRankLegend() || tp.rankLegendColorsForTest().containsKey( "Myotis" )
                || !tp.rankLegendColorsForTest().containsKey( "Molossidae" ) ) {
            fail( ok, "named boxes + bars: the legend keys the bars only, got " + tp.rankLegendColorsForTest().keySet() );
        }
        if ( tp.cladeBandRightReserveForTest() < ( unnamed_reserve + label_h ) ) {
            fail( ok, "the names beside the boxes take room past the labels before the column: reserve "
                    + unnamed_reserve + " -> " + tp.cladeBandRightReserveForTest() );
        }
        // named boxes + brackets at another rank: nothing needs a key, so no legend at all
        nested.set( 1, new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BRACKETS ) );
        tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
        if ( tp.hasRankLegend() || !tp.cladeBoxNamesShownForTest() ) {
            fail( ok, "named boxes + brackets: the names and the brackets are the keys, no legend; got legend "
                    + tp.hasRankLegend() + ", names " + tp.cladeBoxNamesShownForTest() );
        }
        nested.set( 1, new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BARS ) );
        // a column at the boxes' OWN rank names each box: no names beside, one rank once, the boxes' colours the key
        // even under brackets (which carry none)
        final List<CladeLevel.Spec> same = new ArrayList<>();
        same.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BOXES ) );
        same.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.HORIZONTAL, CLADE_VIS.BRACKETS ) );
        tp.setCladeLevels( same, CLADE_VIS.BOXES, true ); // the names choice is still "beside"
        if ( tp.cladeBoxNamesShownForTest() || tp.cladeBoxNamesDrawnForTest().size() > 0 ) {
            fail( ok, "a bracket at the boxes' rank names each box: no names beside the boxes" );
        }
        if ( !tp.cladeLevelRanks().equals( Arrays.asList( "genus" ) ) || ( tp.cladeLevelStartsForTest().length != 1 )
                || ( tp.getCladeBandsMode() != CLADE_VIS.BRACKETS ) ) {
            fail( ok, "genus boxes + genus brackets: the rank once, one column, mode BRACKETS; got "
                    + tp.cladeLevelRanks() + " / " + tp.cladeLevelStartsForTest().length + " / " + tp.getCladeBandsMode() );
        }
        if ( Math.abs( tp.cladeLevelStartsForTest()[ 0 ] - bars_start ) > 0.5f ) {
            fail( ok, "with the names not drawn the column starts where a lone column starts, got "
                    + tp.cladeLevelStartsForTest()[ 0 ] );
        }
        if ( !tp.hasRankLegend() || !tp.rankLegendColorsForTest().containsKey( "Myotis" ) ) {
            fail( ok, "boxes + brackets at one rank: the boxes' colours are the key (a legend), got "
                    + tp.rankLegendColorsForTest().keySet() );
        }
        // ...and with bars at that rank the bar wears its box's colour, taxon by taxon (one pair of marks per clade)
        same.set( 1, new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BARS ) );
        tp.setCladeLevels( same, CLADE_VIS.BOXES, true );
        final List<List<CladeBand>> levels = tp.cladeLevelBandsForTest();
        if ( levels.size() != 2 ) {
            fail( ok, "genus boxes + genus bars: two levels drawn, got " + levels.size() );
        }
        else {
            final Map<String, Color> box_colors = new HashMap<>();
            for ( final CladeBand b : levels.get( 0 ) ) {
                box_colors.put( b.getTaxon(), b.getColor() );
            }
            for ( final CladeBand b : levels.get( 1 ) ) {
                if ( !b.getColor().equals( box_colors.get( b.getTaxon() ) ) ) {
                    fail( ok, "the bar of " + b.getTaxon() + " must wear its box's colour: " + b.getColor() + " vs "
                            + box_colors.get( b.getTaxon() ) );
                    break;
                }
            }
            if ( levels.get( 0 ).size() != levels.get( 1 ).size() ) {
                fail( ok, "one rank, one set of clades for both marks: " + levels.get( 0 ).size() + " vs "
                        + levels.get( 1 ).size() );
            }
        }
        // a box rank the tree cannot place leaves the next box rank to draw (as several ranks in the old BOXES mode
        // did): subspecies is finer than genus and the bat demo names none -- genus draws, and with the colours it
        // has alone (the unplaced box rank is no level of the figure, so it must not hue-band the drawn one)
        tp.setCladeBoxLabels( false );
        tp.setCladeLevels( specs( "genus" ), CLADE_VIS.BOXES, true );
        final Map<String, Color> alone = new HashMap<>();
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            alone.put( b.getTaxon(), b.getColor() );
        }
        tp.setCladeLevels( specs( "subspecies", "genus" ), CLADE_VIS.BOXES, true );
        if ( !tp.cladeLevelRanks().equals( Arrays.asList( "genus" ) ) || ( tp.drawnCladeBoxCount() <= 0 )
                || ( tp.cladeLevelSpecCountForTest() != 2 ) ) {
            fail( ok, "an unplaceable finer box rank must leave genus boxes drawn (2 specs kept): ranks "
                    + tp.cladeLevelRanks() + ", boxes " + tp.drawnCladeBoxCount() + ", specs "
                    + tp.cladeLevelSpecCountForTest() );
        }
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            if ( !b.getColor().equals( alone.get( b.getTaxon() ) ) ) {
                fail( ok, "an unplaceable box rank must not change the drawn boxes' colours: " + b.getTaxon() + " "
                        + b.getColor() + " vs " + alone.get( b.getTaxon() ) );
                break;
            }
        }
        // the box at the BROADER rank: the first-level hook reads the boxes, the dialog reopens on their rank
        final List<CladeLevel.Spec> broad = new ArrayList<>();
        broad.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BOXES ) );
        broad.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BARS ) );
        tp.setCladeLevels( broad, CLADE_VIS.BOXES, true );
        boolean has_family = false;
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            has_family |= "Molossidae".equals( b.getTaxon() );
        }
        if ( !has_family || !"family".equals( MainFrame.cladeRankToPreselect( tp ) ) ) {
            fail( ok, "family boxes + genus bars: the first-level hook reads the boxes and the dialog preselects family;"
                    + " got boxes " + has_family + ", preselect " + MainFrame.cladeRankToPreselect( tp ) );
        }
        tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
        if ( !"genus".equals( MainFrame.cladeRankToPreselect( tp ) ) ) {
            fail( ok, "genus boxes + family bars: the dialog preselects the boxes' rank, got "
                    + MainFrame.cladeRankToPreselect( tp ) );
        }
        final List<CladeLevel.Spec> family_bars = new ArrayList<>();
        family_bars.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BARS ) );
        tp.setCladeLevels( family_bars, CLADE_VIS.BARS, true );
        if ( !"family".equals( MainFrame.cladeRankToPreselect( tp ) ) ) {
            fail( ok, "bars alone: the dialog preselects the finest column, got " + MainFrame.cladeRankToPreselect( tp ) );
        }
        // a colour picked on a legend row is stored under the rank the legend KEYS: named genus boxes + family bars
        // key family alone (no section blocks), so a pick on "Molossidae" must go to family, not to the finest drawn
        // rank -- where the bars would never find it
        tp.setCladeBoxLabels( true );
        tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
        if ( !"family".equals( tp.rankLegendRankForTest( "Molossidae" ) ) ) {
            fail( ok, "named genus boxes + family bars: a legend pick on Molossidae is stored under family, got "
                    + tp.rankLegendRankForTest( "Molossidae" ) );
        }
        tp.setCladeBoxLabels( false );
        tp.setCladeLevels( nested, CLADE_VIS.BOXES, true ); // both keyed: the section blocks decide, as before
        if ( !"family".equals( tp.rankLegendRankForTest( "Molossidae" ) )
                || !"genus".equals( tp.rankLegendRankForTest( "Myotis" ) ) ) {
            fail( ok, "genus boxes in the legend + family bars: each row under its own rank, got "
                    + tp.rankLegendRankForTest( "Molossidae" ) + " / " + tp.rankLegendRankForTest( "Myotis" ) );
        }

        // the figure keeps the per-level styles; an older figure (no styles) opens as it always did
        tp.setCladeBoxLabels( false );
        tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
        final String saved = FigureSpec.capture( tp ).toPropertyValue();
        final FigureSpec captured = FigureSpec.parse( saved ); // the stored text escapes ~ and |; read it back
        if ( !"genus~VERTICAL~BOXES|family~VERTICAL~BARS".equals( captured.get( "cladelevels" ) )
                || !"BARS".equals( captured.get( "clademode" ) ) ) {
            fail( ok, "the figure must name each level's style and the columns' mode: " + captured.get( "cladelevels" )
                    + " / " + captured.get( "clademode" ) );
        }
        tp.clearCladeBands();
        FigureSpec.parse( saved ).applyTo( tp );
        if ( !"genus".equals( tp.cladeBoxRank() ) || !tp.cladeColumnRanks().equals( Arrays.asList( "family" ) )
                || ( tp.getCladeBandsMode() != CLADE_VIS.BARS ) ) {
            fail( ok, "the saved figure must bring the boxes and the bars back: box " + tp.cladeBoxRank() + ", columns "
                    + tp.cladeColumnRanks() );
        }
        FigureSpec.parse( "v1;clademode=BOXES;cladeskip=true;cladelevels=genus~VERTICAL|family~VERTICAL" ).applyTo( tp );
        if ( !"genus".equals( tp.cladeBoxRank() ) || !tp.cladeColumnRanks().isEmpty() ) {
            fail( ok, "an older BOXES figure of two ranks must open as one box level, the finest: box "
                    + tp.cladeBoxRank() + ", columns " + tp.cladeColumnRanks() );
        }
        FigureSpec.parse( "v1;clademode=BRACKETS;cladeskip=true;cladelevels=genus~VERTICAL|family~VERTICAL" ).applyTo( tp );
        if ( ( tp.cladeBoxRank() != null ) || !tp.cladeColumnRanks().equals( Arrays.asList( "genus", "family" ) ) ) {
            fail( ok, "an older BRACKETS figure must open as two bracket columns: box " + tp.cladeBoxRank()
                    + ", columns " + tp.cladeColumnRanks() );
        }
        // the wash really is painted with a column beside it: at a point just past a two-child genus root, on its
        // midline (where no branch runs), the boxes + bars render differs from the bars-only render, and the latter
        // is the plain background there
        tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
        paint( tp );
        final BufferedImage with_boxes = render( tp );
        final List<int[]> probes = new ArrayList<>();
        for ( final CladeBand b : tp.cladeLevelBandsForTest().get( 0 ) ) {
            if ( b.getRoot().getNumberOfDescendants() == 2 ) {
                probes.add( new int[] { Math.round( b.getRoot().getXcoord() ) + 6, Math.round( b.getRoot().getYcoord() ) } );
            }
        }
        final List<CladeLevel.Spec> bars_only = new ArrayList<>();
        bars_only.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BARS ) );
        tp.setCladeLevels( bars_only, CLADE_VIS.BARS, true );
        paint( tp );
        final BufferedImage without = render( tp );
        final Color bg = tp.getTreeColorSet().getBackgroundColor();
        int washed = 0;
        for ( final int[] p : probes ) {
            if ( ( p[ 0 ] < 0 ) || ( p[ 1 ] < 0 ) || ( p[ 0 ] >= with_boxes.getWidth() ) || ( p[ 1 ] >= with_boxes.getHeight() ) ) {
                continue;
            }
            final Color a = new Color( with_boxes.getRGB( p[ 0 ], p[ 1 ] ) );
            final Color b = new Color( without.getRGB( p[ 0 ], p[ 1 ] ) );
            if ( near( b, bg, 8 ) && !near( a, b, 8 ) ) {
                ++washed;
            }
        }
        if ( probes.isEmpty() || ( washed == 0 ) ) {
            fail( ok, "with a bar column beside them the boxes' wash must still be painted: " + washed + " of "
                    + probes.size() + " probes washed" );
        }

        // circular: boxes named in the legend keep their sectors within the labels and reserve NOTHING past them --
        // alone (the rule since they exist) and beside a column, whose ring then starts where it starts alone; named
        // beside the boxes, the names are reserved before the ring
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        tp.setCladeBoxLabels( false );
        tp.setCladeLevels( bars_only, CLADE_VIS.BARS, true );
        final int ring_alone = tp.circularCladeBandReserveForTest();
        tp.setCladeLevels( specs( "genus" ), CLADE_VIS.BOXES, true );
        final int boxes_alone = tp.circularCladeBandReserveForTest();
        tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
        final int both = tp.circularCladeBandReserveForTest();
        tp.setCladeBoxLabels( true );
        tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
        final int named_both = tp.circularCladeBandReserveForTest();
        tp.setCladeBoxLabels( false );
        if ( ( boxes_alone != 0 ) || ( both != ring_alone ) || ( ring_alone <= 0 ) || ( named_both < ( both + label_h ) ) ) {
            fail( ok, "circular reserve: boxes in the legend 0, with a column the column's own, named boxes more; got "
                    + boxes_alone + " / " + both + " (column alone " + ring_alone + ") / " + named_both );
        }
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );

        // and the combination paints in all five layouts, boxes named in the legend and beside them
        for ( final boolean beside : new boolean[] { false, true } ) {
            tp.setCladeBoxLabels( beside );
            for ( final Object[] layout : new Object[][] {
                    { Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, Options.TREE_ORIENTATION.ROOT_LEFT },
                    { Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, Options.TREE_ORIENTATION.ROOT_TOP },
                    { Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, Options.TREE_ORIENTATION.ROOT_BOTTOM },
                    { Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, Options.TREE_ORIENTATION.ROOT_LEFT },
                    { Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED, Options.TREE_ORIENTATION.ROOT_LEFT } } ) {
                frame.getCurrentTreePanel().setTreeOrientation( ( Options.TREE_ORIENTATION ) layout[ 1 ] );
                tp.setPhylogenyGraphicsType( ( Options.PHYLOGENY_GRAPHICS_TYPE ) layout[ 0 ] );
                tp.setCladeLevels( nested, CLADE_VIS.BOXES, true );
                try {
                    if ( AptxUtil.renderPhylogenyToImage( 1000, 900, tp, frame.getOptions(), false, 1,
                                                          false ) == null ) {
                        fail( ok, "boxes + bars in " + layout[ 0 ] + "/" + layout[ 1 ] + " produced no image" );
                    }
                }
                catch ( final Exception e ) {
                    fail( ok, "boxes + bars in " + layout[ 0 ] + "/" + layout[ 1 ] + " threw " + e );
                }
            }
        }
        tp.setCladeBoxLabels( false );
    }

    /**
     * The carnivore demo, opened the way File &gt; Open opens it: its saved figure draws genus boxes inside family
     * brackets, and the five one-tip genera keep their wash -- a shaded box is never skipped as a single-member clade,
     * whatever "Skip single-member clades" says -- while every family (two tips or more) gets its bracket. The counts
     * are the ones the demo's description and README row quote.
     */
    private static void oneTipBoxesDraw( final MainFrame frame, final boolean[] ok ) {
        final File file = new File( System.getProperty( "user.dir" ), "forester/demo/clade-boxes-with-brackets.xml" );
        if ( !file.exists() ) {
            fail( ok, "missing demo tree " + file.getAbsolutePath() );
            return;
        }
        final Phylogeny phy;
        try {
            phy = ParserBasedPhylogenyFactory.getInstance().create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
        }
        catch ( final Exception e ) {
            fail( ok, "the carnivore demo could not be read: " + e );
            return;
        }
        AptxUtil.addPhylogeniesToTabs( new Phylogeny[] { phy }, "", "", frame.getConfiguration(), frame.getMainPanel() );
        final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
        if ( !"genus".equals( tp.cladeBoxRank() ) || !tp.cladeColumnRanks().equals( Arrays.asList( "family" ) )
                || ( tp.getCladeBandsMode() != CLADE_VIS.BRACKETS ) || !tp.isCladeBandsSkipSingletons() ) {
            fail( ok, "the carnivore demo's figure must open as genus boxes inside family brackets, singletons skipped;"
                    + " got box " + tp.cladeBoxRank() + ", columns " + tp.cladeColumnRanks() + ", mode "
                    + tp.getCladeBandsMode() + ", skip " + tp.isCladeBandsSkipSingletons() );
            return;
        }
        final List<CladeBand> genera = tp.cladeLevelBandsForTest().get( 0 );
        int one_tip = 0;
        for ( final CladeBand b : genera ) {
            one_tip += b.getRoot().isExternal() ? 1 : 0;
        }
        if ( ( genera.size() != 12 ) || ( one_tip != 5 ) ) {
            fail( ok, "precondition: 12 genera, 5 of them one tip; got " + genera.size() + " / " + one_tip );
            return;
        }
        if ( ( tp.drawnCladeBoxCount() != 12 ) || ( tp.drawnCladeColumnCount() != 5 ) ) {
            fail( ok, "every genus box draws, the one-tip ones included, and every family bracket: drawn "
                    + tp.drawnCladeBoxCount() + " boxes, " + tp.drawnCladeColumnCount() + " brackets" );
        }
        if ( tp.cladeBoxNamesShownForTest() || !tp.hasRankLegend() || !tp.rankLegendColorsForTest().containsKey( "Lutra" )
                || tp.rankLegendColorsForTest().containsKey( "Mustelidae" ) ) {
            fail( ok, "named in the legend: the genera (the one-tip Lutra too) are legend rows, the brackets' families "
                    + "are not; got " + tp.rankLegendColorsForTest().keySet() );
        }
        try {
            if ( AptxUtil.renderPhylogenyToImage( 1000, 900, tp, frame.getOptions(), false, 1, false ) == null ) {
                fail( ok, "the carnivore demo produced no image" );
            }
        }
        catch ( final Exception e ) {
            fail( ok, "rendering the carnivore demo threw " + e );
        }
        // the legend's "no value" row counts the KEYED rank's own unplaced tips: take the genus off two tips (they keep
        // their family) -- genus boxes alone report 2, genus boxes + family brackets still 2 (the family's placed
        // tips must not hide them), and two keyed ranks report none, as sectioned legends always did
        int stripped = 0;
        for ( final org.forester.phylogeny.PhylogenyNode tip : phy.getExternalNodes() ) {
            if ( stripped < 2 ) {
                stripped += tip.getNodeData().getProperties().getProperties()
                        .removeIf( p -> "ncbi:genus".equals( p.getRef() ) ) ? 1 : 0;
            }
        }
        if ( stripped != 2 ) {
            fail( ok, "precondition: two tips stripped of their genus, got " + stripped );
            return;
        }
        final List<CladeLevel.Spec> boxes_only = new ArrayList<>();
        boxes_only.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BOXES ) );
        tp.setCladeLevels( boxes_only, CLADE_VIS.BOXES, true );
        final int alone_missing = tp.rankLegendMissingForTest();
        final List<CladeLevel.Spec> with_brackets = new ArrayList<>( boxes_only );
        with_brackets.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BRACKETS ) );
        tp.setCladeLevels( with_brackets, CLADE_VIS.BOXES, true );
        final int with_missing = tp.rankLegendMissingForTest();
        final List<CladeLevel.Spec> with_bars = new ArrayList<>( boxes_only );
        with_bars.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BARS ) );
        tp.setCladeLevels( with_bars, CLADE_VIS.BOXES, true );
        final int sectioned_missing = tp.rankLegendMissingForTest();
        if ( ( alone_missing != 2 ) || ( with_missing != 2 ) || ( sectioned_missing != 0 ) ) {
            fail( ok, "the no-value row: genus boxes alone 2, with family brackets 2, with family bars (two keyed) 0; got "
                    + alone_missing + " / " + with_missing + " / " + sectioned_missing );
        }
    }

    private static void paint( final TreePanel tp ) {
        tp.paintImmediately( 0, 0, tp.getWidth(), tp.getHeight() );
        tp.paintImmediately( 0, 0, tp.getWidth(), tp.getHeight() );
    }

    private static BufferedImage render( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( tp.getWidth(), tp.getHeight(), BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            tp.printAll( g );
        }
        finally {
            g.dispose();
        }
        return img;
    }

    private static boolean near( final Color a, final Color b, final int tol ) {
        return ( Math.abs( a.getRed() - b.getRed() ) <= tol ) && ( Math.abs( a.getGreen() - b.getGreen() ) <= tol )
                && ( Math.abs( a.getBlue() - b.getBlue() ) <= tol );
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
