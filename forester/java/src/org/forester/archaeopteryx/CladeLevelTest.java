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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.forester.archaeopteryx.TreePanel.CLADE_LABEL_ANGLE;
import org.forester.archaeopteryx.TreePanel.CLADE_VIS;

/**
 * "Annotate Clades by Rank" can draw up to three ranks at once as nested bar/bracket columns, one of them instead a
 * level of shaded boxes. Which rank ends up where is DERIVED, never chosen: {@link CladeLevel#order} sorts
 * finest-rank-first so the finest is drawn nearest the tips and the broadest outermost, which is the only order in
 * which the nesting reads correctly -- a broad clade's mark has to span the several finer marks inside it. If this
 * ordering breaks, the figure silently lies about the hierarchy, so it is pinned down here. The styles: each spec
 * carries its own, a spec without one takes the mode, every box rank is kept (the finest that places anything is
 * the one drawn, decided where the bands are known), the columns are capped so the drawn levels never exceed three,
 * and a rank may appear once as boxes and once as a column (the bracket that names the shaded clade) but never as
 * two columns. Pure + headless.
 */
public final class CladeLevelTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CladeLevel: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return rankDepth() && ordering() && styles() && clearAffordance();
    }

    /**
     * The level-1 rank chooser's "remove the clade annotations" entry. Levels 2 and 3 always had a "(none)";
     * level 1 did not, so once clades were annotated there was no way to take them off -- {@code clearCladeBands()}
     * had no caller outside the tests. The entry is offered ONLY when there is something to remove, and it must not
     * disturb the pre-selection of the last-used rank, whose index it shifts by one.
     */
    private static boolean clearAffordance() {
        final String[] ranks = { "genus (12) (98%)", "family (5) (100%)", "order (2) (100%)" };
        if ( RankSourcePanel.withFirst( ranks, null ) != ranks ) {
            return fail( "with nothing drawn, the chooser should offer the ranks unchanged" );
        }
        final String[] with = RankSourcePanel.withFirst( ranks, MainFrame.CLADE_LEVEL_CLEAR );
        if ( ( with.length != ranks.length + 1 ) || !MainFrame.CLADE_LEVEL_CLEAR.equals( with[ 0 ] ) ) {
            return fail( "with bands drawn, the \"remove\" entry should come FIRST, got "
                    + java.util.Arrays.toString( with ) );
        }
        for ( int i = 0; i < ranks.length; i++ ) {
            if ( !ranks[ i ].equals( with[ i + 1 ] ) ) {
                return fail( "the ranks must keep their order after the \"remove\" entry, got "
                        + java.util.Arrays.toString( with ) );
            }
        }
        // the pre-selection looks the rank up IN THE MODEL SHOWN, so it has to survive the extra leading entry
        if ( AptxUtil.indexOfRank( ranks, "family" ) != 1 ) {
            return fail( "precondition: 'family' should be at index 1 without the remove entry" );
        }
        if ( AptxUtil.indexOfRank( with, "family" ) != 2 ) {
            return fail( "'family' should shift to index 2 once the remove entry is prepended, got "
                    + AptxUtil.indexOfRank( with, "family" ) );
        }
        // ...and the remove entry itself must never read as a rank (it starts with '(' so it has no bare prefix)
        if ( AptxUtil.indexOfRank( with, "none" ) >= 0 ) {
            return fail( "the \"remove\" entry must not be matchable as a rank" );
        }
        // A DESTRUCTIVE entry must never be what the dialog opens on. With the remove entry at index 0, a chooser
        // that simply fell back to "index 0" would greet the user with "remove the clade annotations" pre-selected
        // -- press OK expecting to confirm a rank and the annotations are gone.
        if ( MainFrame.cladeRankPreselectIndex( with, null, null ) != 1 ) {
            return fail( "with nothing remembered, the chooser must open on the first real RANK, not the remove "
                    + "entry, got index " + MainFrame.cladeRankPreselectIndex( with, null, null ) );
        }
        if ( MainFrame.cladeRankPreselectIndex( with, "family", null ) != 2 ) {
            return fail( "the remembered rank should win, got index "
                    + MainFrame.cladeRankPreselectIndex( with, "family", null ) );
        }
        if ( MainFrame.cladeRankPreselectIndex( with, null, "order" ) != 3 ) {
            return fail( "with nothing remembered, the rank currently DRAWN should be pre-selected, got index "
                    + MainFrame.cladeRankPreselectIndex( with, null, "order" ) );
        }
        if ( MainFrame.cladeRankPreselectIndex( ranks, null, null ) != 0 ) {
            return fail( "without the remove entry, index 0 is a perfectly good default" );
        }
        // degenerate: the remove entry with no ranks behind it must not index past the end
        if ( MainFrame.cladeRankPreselectIndex( new String[] { MainFrame.CLADE_LEVEL_CLEAR }, null, null ) != 0 ) {
            return fail( "a one-entry model must not index out of bounds" );
        }
        return true;
    }

    private static boolean rankDepth() {
        // finer rank => strictly greater depth, right down the chain a user actually picks
        final String[] broad_to_fine = { "phylum", "class", "order", "family", "subfamily", "genus", "species" };
        for ( int i = 1; i < broad_to_fine.length; i++ ) {
            final int prev = CladeLevel.rankDepth( broad_to_fine[ i - 1 ] );
            final int cur = CladeLevel.rankDepth( broad_to_fine[ i ] );
            if ( cur <= prev ) {
                return fail( broad_to_fine[ i ] + " must be FINER than " + broad_to_fine[ i - 1 ] + ", got depths "
                        + cur + " vs " + prev );
            }
        }
        if ( CladeLevel.rankDepth( "GENUS" ) != CladeLevel.rankDepth( "genus" ) ) {
            return fail( "rank depth must be case-insensitive" );
        }
        // "unknown"/"other" sit at the END of the underlying rank list but carry no real depth; treating them as
        // the finest ranks would wedge a rank of unknown breadth up against the tips
        for ( final String no_depth : new String[] { "unknown", "other", "not-a-rank", "", null } ) {
            if ( CladeLevel.rankDepth( no_depth ) != -1 ) {
                return fail( "\"" + no_depth + "\" has no meaningful depth and must report -1, got "
                        + CladeLevel.rankDepth( no_depth ) );
            }
        }
        return true;
    }

    private static boolean ordering() {
        // the headline behavior: whatever order the user picks, the FINEST rank is drawn first (nearest the tips)
        for ( final List<String> picked : Arrays.asList( Arrays.asList( "family", "genus", "order" ),
                                                         Arrays.asList( "order", "family", "genus" ),
                                                         Arrays.asList( "genus", "order", "family" ) ) ) {
            final List<String> got = ranks( CladeLevel.order( specs( picked ), CLADE_VIS.BARS ) );
            if ( !got.equals( Arrays.asList( "genus", "family", "order" ) ) ) {
                return fail( "picked " + picked + " must draw as [genus, family, order], got " + got );
            }
        }
        // blanks are dropped, and a repeated rank cannot claim two columns (it would draw the same marks twice)
        final List<CladeLevel.Spec> messy = new ArrayList<>();
        messy.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL ) );
        messy.add( null );
        messy.add( new CladeLevel.Spec( "", CLADE_LABEL_ANGLE.VERTICAL ) );
        messy.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL ) );
        messy.add( new CladeLevel.Spec( "FAMILY", CLADE_LABEL_ANGLE.VERTICAL ) ); // same rank, other case
        final List<String> cleaned = ranks( CladeLevel.order( messy, CLADE_VIS.BARS ) );
        if ( !cleaned.equals( Arrays.asList( "genus", "family" ) ) ) {
            return fail( "blanks/nulls/duplicates must be dropped, got " + cleaned );
        }
        // never more columns than the layout can carry
        final List<String> many = Arrays.asList( "species", "genus", "subfamily", "family", "order", "class" );
        final List<CladeLevel.Spec> capped = CladeLevel.order( specs( many ), CLADE_VIS.BARS );
        if ( capped.size() != CladeLevel.MAX_LEVELS ) {
            return fail( "at most " + CladeLevel.MAX_LEVELS + " levels may survive, got " + capped.size() );
        }
        if ( !ranks( capped ).equals( Arrays.asList( "species", "genus", "subfamily" ) ) ) {
            return fail( "the cap must keep the FINEST ranks (they are the ones nearest the tips), got "
                    + ranks( capped ) );
        }
        // a rank of unknown depth goes OUTERMOST (the safe end) rather than being taken for the finest
        final List<String> with_unknown = ranks( CladeLevel.order( specs( Arrays.asList( "mystery", "genus" ) ), CLADE_VIS.BARS ) );
        if ( !with_unknown.equals( Arrays.asList( "genus", "mystery" ) ) ) {
            return fail( "an unrecognized rank must sort outermost, got " + with_unknown );
        }
        // equal-depth (here: both unknown) specs keep the order the user gave them -- a stable sort
        final List<String> stable = ranks( CladeLevel.order( specs( Arrays.asList( "zeta", "alpha" ) ), CLADE_VIS.BARS ) );
        if ( !stable.equals( Arrays.asList( "zeta", "alpha" ) ) ) {
            return fail( "equal-depth specs must keep the user's order, got " + stable );
        }
        // the per-level label angle survives ordering (it is the whole point of the Spec)
        final List<CladeLevel.Spec> angled = new ArrayList<>();
        angled.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.HORIZONTAL ) );
        angled.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL ) );
        final List<CladeLevel.Spec> ordered = CladeLevel.order( angled, CLADE_VIS.BARS );
        if ( ( ordered.get( 0 ).getLabelAngle() != CLADE_LABEL_ANGLE.VERTICAL )
                || ( ordered.get( 1 ).getLabelAngle() != CLADE_LABEL_ANGLE.HORIZONTAL ) ) {
            return fail( "each level's label angle must travel with its rank through the re-ordering" );
        }
        if ( !CladeLevel.order( null, CLADE_VIS.BARS ).isEmpty() ) {
            return fail( "a null spec list must yield no levels" );
        }
        return true;
    }

    /**
     * The per-level styles (Christian, 2026-11-04, "go, box + two"): a spec without a style takes the mode; a box
     * and a column may share a rank (one pair of marks for the same clades) while two columns may not; a second
     * box level is dropped, the finest kept (what a tab of several box levels drew before styles were per level);
     * the cap counts the box; the style travels with its rank through the re-ordering.
     */
    private static boolean styles() {
        // the mode fills in what a spec does not say -- every caller from before, and every figure before 0.11.200
        final List<CladeLevel.Spec> plain = CladeLevel.order( specs( Arrays.asList( "genus", "family" ) ),
                                                              CLADE_VIS.BRACKETS );
        for ( final CladeLevel.Spec s : plain ) {
            if ( s.getStyle() != CLADE_VIS.BRACKETS ) {
                return fail( "a spec without a style must take the mode, got " + plain );
            }
        }
        // a box and a column at one rank: both kept, the box first (before a column at its own rank, whichever was picked first)
        final List<CladeLevel.Spec> pair = new ArrayList<>();
        pair.add( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BOXES ) );
        pair.add( new CladeLevel.Spec( "Genus", CLADE_LABEL_ANGLE.HORIZONTAL, CLADE_VIS.BRACKETS ) );
        pair.add( new CladeLevel.Spec( "family", CLADE_LABEL_ANGLE.VERTICAL, CLADE_VIS.BRACKETS ) );
        final List<CladeLevel.Spec> both = CladeLevel.order( pair, CLADE_VIS.BOXES );
        if ( ( both.size() != 3 ) || !both.get( 0 ).isBox() || !"genus".equals( both.get( 0 ).getRank() )
                || both.get( 1 ).isBox() || !"Genus".equals( both.get( 1 ).getRank() )
                || ( both.get( 1 ).getStyle() != CLADE_VIS.BRACKETS )
                || ( both.get( 1 ).getLabelAngle() != CLADE_LABEL_ANGLE.HORIZONTAL )
                || !"family".equals( both.get( 2 ).getRank() ) ) {
            return fail( "a box and a column at one rank must both survive, the box first, got " + both );
        }
        // ...but two columns at one rank are one column; every BOX rank is kept (which one draws is decided where the
        // bands are known: the finest that places anything, so a box rank the tree cannot place leaves the next)
        final List<CladeLevel.Spec> twice = new ArrayList<>();
        twice.add( new CladeLevel.Spec( "family", null, CLADE_VIS.BOXES ) );
        twice.add( new CladeLevel.Spec( "genus", null, CLADE_VIS.BARS ) );
        twice.add( new CladeLevel.Spec( "genus", null, CLADE_VIS.BARS ) );
        twice.add( new CladeLevel.Spec( "genus", null, CLADE_VIS.BOXES ) );
        final List<CladeLevel.Spec> once = CladeLevel.order( twice, CLADE_VIS.BARS );
        if ( ( once.size() != 3 ) || !once.get( 0 ).isBox() || !"genus".equals( once.get( 0 ).getRank() )
                || once.get( 1 ).isBox() || !"genus".equals( once.get( 1 ).getRank() )
                || !once.get( 2 ).isBox() || !"family".equals( once.get( 2 ).getRank() ) ) {
            return fail( "two columns at one rank are one; every box rank is kept, finest first; got " + once );
        }
        // the old BOXES-mode tab of three ranks: three box specs, finest first (one of them draws)
        final List<CladeLevel.Spec> old = CladeLevel.order( specs( Arrays.asList( "order", "genus", "family" ) ),
                                                            CLADE_VIS.BOXES );
        if ( ( old.size() != 3 ) || !ranks( old ).equals( Arrays.asList( "genus", "family", "order" ) )
                || !old.get( 0 ).isBox() || !old.get( 1 ).isBox() || !old.get( 2 ).isBox() ) {
            return fail( "several ranks in BOXES mode keep every box rank, finest first, got " + old );
        }
        // the cap counts the box and never drops it: a box and three columns are the box and the two finest columns
        final List<CladeLevel.Spec> four = new ArrayList<>();
        four.add( new CladeLevel.Spec( "species", null, CLADE_VIS.BOXES ) );
        four.add( new CladeLevel.Spec( "order", null, CLADE_VIS.BARS ) );
        four.add( new CladeLevel.Spec( "genus", null, CLADE_VIS.BARS ) );
        four.add( new CladeLevel.Spec( "family", null, CLADE_VIS.BARS ) );
        final List<CladeLevel.Spec> capped = CladeLevel.order( four, CLADE_VIS.BARS );
        if ( ( capped.size() != CladeLevel.MAX_LEVELS ) || !ranks( capped ).equals( Arrays.asList( "species", "genus",
                                                                                                   "family" ) )
                || !capped.get( 0 ).isBox() || capped.get( 1 ).isBox() || capped.get( 2 ).isBox() ) {
            return fail( "the cap counts the box level and keeps the finest columns, got " + capped );
        }
        // ...even when the box is the BROADEST rank: the broadest COLUMN goes, the box stays
        final List<CladeLevel.Spec> broad_box = new ArrayList<>();
        broad_box.add( new CladeLevel.Spec( "species", null, CLADE_VIS.BARS ) );
        broad_box.add( new CladeLevel.Spec( "genus", null, CLADE_VIS.BARS ) );
        broad_box.add( new CladeLevel.Spec( "family", null, CLADE_VIS.BARS ) );
        broad_box.add( new CladeLevel.Spec( "order", null, CLADE_VIS.BOXES ) );
        final List<CladeLevel.Spec> kept_box = CladeLevel.order( broad_box, CLADE_VIS.BARS );
        if ( ( kept_box.size() != CladeLevel.MAX_LEVELS )
                || !ranks( kept_box ).equals( Arrays.asList( "species", "genus", "order" ) )
                || !kept_box.get( 2 ).isBox() ) {
            return fail( "the cap must drop the broadest COLUMN, never the box, got " + kept_box );
        }
        // three columns and no box: the cap is three columns
        if ( CladeLevel.order( specs( Arrays.asList( "species", "genus", "family" ) ), CLADE_VIS.BARS ).size() != 3 ) {
            return fail( "three columns without a box all survive" );
        }
        // a level cannot be built from an unresolved spec: the paints switch on the style
        try {
            new CladeLevel( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL ), new ArrayList<CladeBand>() );
            return fail( "a CladeLevel from a spec without a style must be refused" );
        }
        catch ( final IllegalArgumentException expected ) {
            // the guard
        }
        // the report's word for the columns: by the specs' styles, both when a figure mixed them
        final List<CladeLevel.Spec> mixed = new ArrayList<>();
        mixed.add( new CladeLevel.Spec( "genus", null, CLADE_VIS.BOXES ) );
        mixed.add( new CladeLevel.Spec( "family", null, CLADE_VIS.BARS ) );
        if ( !"bar(s)".equals( MainFrame.cladeColumnKind( mixed ) ) ) {
            return fail( "boxes + bars report bar(s), got " + MainFrame.cladeColumnKind( mixed ) );
        }
        mixed.set( 1, new CladeLevel.Spec( "family", null, CLADE_VIS.BRACKETS ) );
        if ( !"bracket(s)".equals( MainFrame.cladeColumnKind( mixed ) ) ) {
            return fail( "boxes + brackets report bracket(s), got " + MainFrame.cladeColumnKind( mixed ) );
        }
        mixed.add( new CladeLevel.Spec( "order", null, CLADE_VIS.BARS ) );
        if ( !"bar(s) and bracket(s)".equals( MainFrame.cladeColumnKind( mixed ) ) ) {
            return fail( "mixed columns report both, got " + MainFrame.cladeColumnKind( mixed ) );
        }
        // a null mode is not a crash: boxes, the dialog's default
        final List<CladeLevel.Spec> none = CladeLevel.order( specs( Arrays.asList( "genus" ) ), null );
        if ( ( none.size() != 1 ) || !none.get( 0 ).isBox() ) {
            return fail( "a null mode must resolve to boxes, got " + none );
        }
        // a resolved spec keeps its own style over the mode
        final CladeLevel.Spec styled = new CladeLevel.Spec( "genus", null, CLADE_VIS.BARS );
        if ( styled.resolved( CLADE_VIS.BOXES ) != styled ) {
            return fail( "a spec with a style must come back unchanged from resolved()" );
        }
        return true;
    }

    private static List<CladeLevel.Spec> specs( final List<String> ranks ) {
        final List<CladeLevel.Spec> out = new ArrayList<>();
        for ( final String r : ranks ) {
            out.add( new CladeLevel.Spec( r, CLADE_LABEL_ANGLE.VERTICAL ) );
        }
        return out;
    }

    private static List<String> ranks( final List<CladeLevel.Spec> specs ) {
        final List<String> out = new ArrayList<>();
        for ( final CladeLevel.Spec s : specs ) {
            out.add( s.getRank() );
        }
        return out;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [CladeLevelTest] " + msg );
        return false;
    }

    private CladeLevelTest() {
    }
}
