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
import java.util.HashSet;
import java.util.Set;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Coverage for the menu-bar "Found / Selected: N" counter ({@link FoundSelectedCounter}) and its tree-side inputs.
 * The pure text formatters and the component's show/hide + two-pass-sweep lifecycle are tested without a display; the
 * end-to-end wiring (a found-set change refreshes the menu-bar counter, and the tree-validated breakdown/total) needs
 * a realized frame, so it is skipped when headless.
 */
public final class FoundSelectedCounterTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "FoundSelectedCounter: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            staticsOk();
            componentOk();
        }
        catch ( final AssertionError e ) {
            System.out.println( "  [FoundSelectedCounterTest] " + e.getMessage() );
            return false;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return wiringOk();
    }

    // ---- pure text (headless) ----------------------------------------------------------------------------------

    private static void staticsOk() {
        ck( FoundSelectedCounter.labelText( 42 ).equals( "Found / Selected: 42" ), "label text" );
        final String tip = FoundSelectedCounter.tooltipText( 10, 7, 5, 2, true, null );
        ck( tip.contains( "10 highlighted nodes" ), "tooltip shows the distinct total: " + tip );
        ck( tip.contains( "A: 7" ) && tip.contains( "B: 5" ) && tip.contains( "2 in both" ),
            "tooltip shows the A/B/both breakdown: " + tip );
        final String solo = FoundSelectedCounter.tooltipText( 3, 3, 0, 0, true, null );
        ck( solo.contains( "A: 3" ) && !solo.contains( "B:" ), "a solo-A tooltip omits B: " + solo );
        ck( FoundSelectedCounter.tooltipText( 1, 1, 0, 0, true, null ).contains( "1 highlighted node (" ),
            "singular wording" );
        // found-set 0 is a manual SELECTION (Search A has no query) -> labelled "Selected", not "A"
        final String sel = FoundSelectedCounter.tooltipText( 4, 3, 1, 0, false, null );
        ck( sel.contains( "Selected: 3" ) && !sel.contains( "A:" ), "a manual selection is labelled Selected: " + sel );
        ck( sel.contains( "B: 1" ), "Search B is still labelled B alongside a selection: " + sel );
        // when the two boxes are COMBINED the tooltip shows the combine description, NOT the per-box breakdown
        final String comb = FoundSelectedCounter.tooltipText( 6, 6, 0, 0, true, "A AND B" );
        ck( comb.contains( "6 highlighted nodes" ) && comb.contains( "A AND B" ),
            "a combined tooltip names the combine mode: " + comb );
        ck( !comb.contains( "A: 6" ), "a combined tooltip omits the per-box breakdown: " + comb );
    }

    // ---- component lifecycle (headless -- no painting needed) --------------------------------------------------

    private static void componentOk() {
        final FoundSelectedCounter c = new FoundSelectedCounter();
        ck( !c.isShowingForTest(), "counter starts hidden" );
        c.setCounts( 5, 3, 4, 2, true, null );
        ck( c.isShowingForTest(), "counter shows when total > 0" );
        ck( c.getTotalForTest() == 5, "total is set" );
        ck( c.isSweepingForTest(), "a change starts the sweep" );
        ck( ( c.getTooltipForTest() != null ) && c.getTooltipForTest().contains( "5" ), "tooltip carries the total" );
        c.runSweepToEndForTest();
        ck( !c.isSweepingForTest(), "the sweep settles static after its passes" );
        c.setCounts( 5, 3, 4, 2, true, null );
        ck( !c.isSweepingForTest(), "an unchanged total does NOT restart the sweep" );
        c.setCounts( 7, 7, 0, 0, true, null );
        ck( c.isSweepingForTest() && ( c.getTotalForTest() == 7 ), "a changed total re-sweeps + updates" );
        c.setCounts( 0, 0, 0, 0, true, null );
        ck( !c.isShowingForTest(), "counter hides at 0" );
        ck( !c.isSweepingForTest(), "the sweep is stopped when hidden" );

        // the sweep is TWO passes (with a rest between them): running it out takes clearly more frames than one pass
        final FoundSelectedCounter c2 = new FoundSelectedCounter();
        c2.setFont( new java.awt.Font( java.awt.Font.SANS_SERIF, java.awt.Font.PLAIN, 13 ) );
        c2.setSize( 160, 20 );
        c2.setCounts( 9, 9, 0, 0, true, null );
        final int ticks = c2.runSweepToEndForTest();
        // at STEP_PX=8 over a 160px width one pass is ~27 frames; two passes + the rest gap is ~60+
        ck( ticks > 45, "the sweep must run TWO passes, not one (took only " + ticks + " frames)" );
    }

    // ---- end-to-end wiring (headful) ---------------------------------------------------------------------------

    private static boolean wiringOk() {
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        final PhylogenyNode[] leaves = new PhylogenyNode[ 3 ];
        try {
            final Phylogeny phy = buildTree( leaves );
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "counter" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    final FoundSelectedCounter counter = mf[ 0 ]._found_selected_counter;
                    ck( ok, counter != null, "the counter is installed in the menu bar" );
                    ck( ok, !counter.isShowingForTest(), "counter hidden when nothing is highlighted" );

                    final PhylogenyNode a = leaves[ 0 ], b = leaves[ 1 ], c = leaves[ 2 ];
                    // A = {a, b}, B = {b, c} -> a in A only, b in both, c in B only
                    tp.setFoundNodes0( new HashSet<>( java.util.Arrays.asList( a.getId(), b.getId() ) ) );
                    tp.setFoundNodes1( new HashSet<>( java.util.Arrays.asList( b.getId(), c.getId() ) ) );

                    final int[] br = tp.foundSelectedBreakdown();
                    ck( ok, ( br[ 0 ] == 2 ) && ( br[ 1 ] == 2 ) && ( br[ 2 ] == 1 ),
                        "breakdown should be A=2, B=2, both=1 (got " + br[ 0 ] + "/" + br[ 1 ] + "/" + br[ 2 ] + ")" );
                    ck( ok, tp.getSearchHitCount() == 3, "the distinct union total should be 3" );
                    ck( ok, counter.isShowingForTest(), "the counter shows after a found-set change" );
                    ck( ok, counter.getTotalForTest() == 3, "the counter shows the union total (3), got "
                            + counter.getTotalForTest() );
                    // box A carries no query here, so found-set-0's nodes are a manual SELECTION -> "Selected", not "A"
                    ck( ok, ( counter.getTooltipForTest() != null ) && counter.getTooltipForTest().contains( "Selected: 2" )
                            && !counter.getTooltipForTest().contains( "A:" ),
                        "with no Search-A query the breakdown labels found-set 0 as Selected: "
                                + counter.getTooltipForTest() );

                    // clearing the found sets hides the counter again
                    tp.setFoundNodes0( null );
                    tp.setFoundNodes1( null );
                    ck( ok, !counter.isShowingForTest(), "the counter hides when the found sets are cleared" );

                    ( (javax.swing.JFrame) mf[ 0 ] ).dispose();
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
                }
            } );
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    private static Phylogeny buildTree( final PhylogenyNode[] out ) {
        final PhylogenyNode a = named( "a" ), b = named( "b" ), c = named( "c" );
        final PhylogenyNode mid = new PhylogenyNode();
        mid.addAsChild( a );
        mid.addAsChild( b );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( mid );
        root.addAsChild( c );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        out[ 0 ] = a;
        out[ 1 ] = b;
        out[ 2 ] = c;
        return phy;
    }

    private static PhylogenyNode named( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static void ck( final boolean cond, final String msg ) {
        if ( !cond ) {
            throw new AssertionError( msg );
        }
    }

    private static void ck( final boolean[] ok, final boolean cond, final String msg ) {
        if ( !cond ) {
            System.out.println( "  [FoundSelectedCounterTest] " + msg );
            ok[ 0 ] = false;
        }
    }

    private FoundSelectedCounterTest() {
    }
}
