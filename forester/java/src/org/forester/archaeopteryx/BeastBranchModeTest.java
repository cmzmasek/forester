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

import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.TIME_AXIS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * The Time | Div switch on a BEAST clock-model tree, driven the way a user drives it, on the demo PAIR:
 * {@code beast-annotations.nex} states a clock rate on every branch and is offered the switch;
 * {@code beast-rate-missing.nex} is the same tree with ONE rate taken out and is not (joint with Archaeopteryx.js:
 * a branch without a rate has no divergence to draw, and 0 would say "no change here").
 * <p>
 * Both files date their five tips at height 0, which is what this also guards: a date of 0 is a date, so the tree
 * has a time layout at all, no tip branch is laid out at 0, and the way back -- by the switch AND by Reset to
 * Defaults -- restores the lengths the file stated. The tree states HEIGHTS (ages, largest at the root), the shape
 * the Auspice parser's child - parent routine lays out at 0 on every branch.
 * <p>
 * Headful; a green no-op when headless.
 */
public final class BeastBranchModeTest {

    private static final String RATED       = "beast-annotations.nex";
    private static final String ONE_UNRATED = "beast-rate-missing.nex";
    /** {@link #RATED} with isolate_B's height taken out: the time layout cannot state its branch. */
    private static final String ONE_UNDATED = "beast-date-missing.nex";
    /** {@link #RATED} with the (D,E) node dated 0.05 BEFORE its parent: a span that runs backwards. */
    private static final String BACKWARDS   = "beast-negative-span.nex";
    /** A recorded divergence (Nextstrain), and its twin with one tip's div taken out. */
    private static final String RECORDED    = "nextstrain-nexus.nex";
    private static final String ONE_WITHOUT_DIV = "nextstrain-div-missing.nex";
    /** A BEAST tree on CALENDAR time (its heights convert on opening), for what the time axis does. */
    private static final String CALENDAR    = "beast-tip-dates.nex";
    private static final String RATE        = BranchLengthLayout.RATE_PROPERTY_REF;
    /** The sum of the eight branch lengths the file states: 1.2 + 1.2 + 0.9 + 0.8 + 0.5 + 0.5 + 0.3 + 1.3. */
    private static final double TIME_SUM    = 6.7;
    /** The sum of rate x time over the eight branches, from the numbers in the file. */
    private static final double DIV_SUM     = 0.02005;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "BEAST branch mode: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            return fixturesOk() && inFrame( RATED, BeastBranchModeTest::switchAndBack )
                    && inFrame( RATED, BeastBranchModeTest::resetFromDivergence )
                    && inFrame( RATED, BeastBranchModeTest::rateLostWhileInTime )
                    && inFrame( RATED, BeastBranchModeTest::rateLostWhileInDivergence )
                    && inFrame( RATED, BeastBranchModeTest::resetAfterTreeReplaced )
                    && inFrame( RATED, BeastBranchModeTest::datesLostThenTime )
                    && inFrame( RATED, BeastBranchModeTest::datesLostThenReset )
                    && inFrame( RATED, BeastBranchModeTest::undoAcrossTheSwitch )
                    && inFrame( RATED, ( f, tp, cp, ok ) -> undoIntoDivergence( tp, cp, ok, false ) )
                    && inFrame( RATED, ( f, tp, cp, ok ) -> undoIntoDivergence( tp, cp, ok, true ) )
                    && inFrame( RATED, BeastBranchModeTest::redoFilesTheLayoutItLeaves )
                    && inFrame( CALENDAR, BeastBranchModeTest::undoBringsTheAxisBack )
                    && inFrame( RATED, BeastBranchModeTest::divergenceInsideASubtree )
                    && inFrame( RATED, BeastBranchModeTest::timeInsideASubtree )
                    && inFrame( ONE_UNRATED, BeastBranchModeTest::subtreeOfATreeWithoutTheSwitch )
                    && inFrame( RATED, BeastBranchModeTest::navigationAloneRewritesNothing )
                    && inFrame( RATED, BeastBranchModeTest::pressedAndPressedBackInASubtree )
                    && inFrame( RATED, BeastBranchModeTest::descendedInDivergence )
                    && inFrame( RATED, BeastBranchModeTest::rateLostInsideASubtree )
                    && inFrame( CALENDAR, BeastBranchModeTest::calendarCladeBackInTime )
                    && inFrame( ONE_UNRATED, BeastBranchModeTest::editorDecides )
                    && inFrame( RATED, BeastBranchModeTest::editorTakesTheDates )
                    && inFrame( ONE_UNRATED, BeastBranchModeTest::notOffered )
                    && inFrame( ONE_UNDATED, ( f, tp, cp, ok ) -> refused( ONE_UNDATED, tp, cp, ok ) )
                    && inFrame( ONE_WITHOUT_DIV, ( f, tp, cp, ok ) -> refused( ONE_WITHOUT_DIV, tp, cp, ok ) )
                    && inFrame( RECORDED, BeastBranchModeTest::recordedTwinIsOffered )
                    && inFrame( BACKWARDS, BeastBranchModeTest::aSpanThatRunsBackwards );
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return fail( "exception: " + t );
        }
    }

    private interface Scenario {

        void run( MainFrame frame, TreePanel tp, ControlPanel cp, boolean[] ok );
    }

    private static boolean inFrame( final String demo, final Scenario scenario ) throws Exception {
        final Phylogeny phy = read( demo );
        if ( phy == null ) {
            return fail( demo + " demo missing / unreadable" );
        }
        final Configuration conf = new Configuration();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "beast" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            try {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
                if ( ( tp == null ) || ( cp == null ) ) {
                    fail( ok, "no tree panel / control panel" );
                    return;
                }
                scenario.run( mf[ 0 ], tp, cp, ok );
            }
            catch ( final Throwable t ) {
                t.printStackTrace();
                fail( ok, "exception: " + t );
            }
            finally {
                mf[ 0 ].dispose();
            }
        } );
        return ok[ 0 ];
    }

    /** The pair is what it is said to be: the same tree, one rate apart, tips at 0 -- read from the files. */
    private static boolean fixturesOk() {
        final Phylogeny rated = read( RATED );
        final Phylogeny unrated = read( ONE_UNRATED );
        if ( ( rated == null ) || ( unrated == null ) ) {
            return fail( "the demo pair is missing / unreadable" );
        }
        if ( ( branches( rated ) != 8 ) || ( ratedBranches( rated ) != 8 ) ) {
            return fail( RATED + " must state a rate on all 8 branches, got " + ratedBranches( rated ) + " of "
                    + branches( rated ) );
        }
        if ( ( branches( unrated ) != 8 ) || ( ratedBranches( unrated ) != 7 ) ) {
            return fail( ONE_UNRATED + " must state a rate on 7 of 8 branches, got " + ratedBranches( unrated )
                    + " of " + branches( unrated ) );
        }
        int zero_tips = 0;
        for( final PhylogenyNode tip : rated.getExternalNodes() ) {
            if ( ( tip.getNodeData().getDate() != null ) && ( tip.getNodeData().getDate().getValue() != null )
                    && ( tip.getNodeData().getDate().getValue().signum() == 0 ) ) {
                ++zero_tips;
            }
        }
        if ( zero_tips != 5 ) {
            return fail( RATED + " must date its 5 tips at exactly 0, got " + zero_tips );
        }
        if ( !near( sum( rated ), TIME_SUM ) || !near( sum( unrated ), TIME_SUM ) ) {
            return fail( "both files must state branch lengths summing to " + TIME_SUM + ", got " + sum( rated )
                    + " and " + sum( unrated ) );
        }
        return true;
    }

    private static void switchAndBack( final MainFrame frame,
                                       final TreePanel tp,
                                       final ControlPanel cp,
                                       final boolean[] ok ) {
        if ( !tp.isBranchLengthToggleApplicable() ) {
            fail( ok, RATED + " states dates and a rate on every branch: the switch must apply" );
        }
        if ( !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "the Time | Div control must be visible for " + RATED );
        }
        final String tip = cp.branchLengthDivTooltipForTest();
        if ( ( tip == null ) || !tip.contains( "DERIVED" ) ) {
            fail( ok, "the Div tooltip must say the divergence is DERIVED from the clock rate, got: " + tip );
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !"Time".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "the default mode must be TIME" );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        final Phylogeny phy = tp.getPhylogeny();
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) {
            fail( ok, "selecting Div must switch the mode" );
        }
        // rate x time, from the numbers the file states: 1.2 x 0.0031, 0.8 x 0.0026, (0.8 - 0.5) x 0.0034
        if ( !near( length( phy, "isolate_A" ), 0.00372 ) || !near( length( phy, "isolate_C" ), 0.00208 )
                || !near( parentLength( phy, "isolate_D" ), 0.00102 ) ) {
            fail( ok, "divergence must be rate x time; got A=" + length( phy, "isolate_A" ) + " C="
                    + length( phy, "isolate_C" ) + " (D,E)=" + parentLength( phy, "isolate_D" ) );
        }
        if ( !near( sum( phy ), DIV_SUM ) ) {
            fail( ok, "the divergence lengths must sum to " + DIV_SUM + ", got " + sum( phy ) );
        }
        if ( zeroBranches( phy ) != 0 ) {
            fail( ok, "no branch of " + RATED + " has zero divergence, got " + zeroBranches( phy )
                    + " (a tip dated 0 read as undated is laid out at 0)" );
        }
        if ( !"subs/site".equals( phy.getDistanceUnit() ) ) {
            fail( ok, "the divergence unit must be subs/site, got " + phy.getDistanceUnit() );
        }
        if ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) {
            fail( ok, "the divergence view has no time axis, got " + tp.effectiveTimeAxisType() );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "selecting Time must switch back" );
        }
        if ( !near( length( phy, "isolate_A" ), 1.2 ) || !near( length( phy, "isolate_C" ), 0.8 )
                || !near( parentLength( phy, "isolate_D" ), 0.3 ) || !near( sum( phy ), TIME_SUM ) ) {
            fail( ok, "switching back must restore the time lengths; got A=" + length( phy, "isolate_A" ) + " C="
                    + length( phy, "isolate_C" ) + " (D,E)=" + parentLength( phy, "isolate_D" ) + " sum="
                    + sum( phy ) );
        }
        if ( !"time".equals( phy.getDistanceUnit() ) ) {
            fail( ok, "heights state no unit: the time unit must read \"time\", got " + phy.getDistanceUnit() );
        }
        if ( !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "the control must still be there after a round trip" );
        }
    }

    private static void resetFromDivergence( final MainFrame frame,
                                             final TreePanel tp,
                                             final ControlPanel cp,
                                             final boolean[] ok ) {
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        final Phylogeny phy = tp.getPhylogeny();
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !near( sum( phy ), DIV_SUM ) ) {
            fail( ok, "reset: the tree must BE in the divergence view first" );
            return;
        }
        tp.resetBranchLengthModeToDefault();
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "Reset must return the mode to TIME" );
        }
        if ( !near( sum( phy ), TIME_SUM ) || !near( length( phy, "isolate_A" ), 1.2 ) ) {
            fail( ok, "Reset must restore the time lengths of a tree dated in heights; got sum=" + sum( phy )
                    + " A=" + length( phy, "isolate_A" ) );
        }
        if ( !"time".equals( phy.getDistanceUnit() ) ) {
            fail( ok, "Reset must not call unit-less heights years; got " + phy.getDistanceUnit() );
        }
    }

    private static void rateLostWhileInTime( final MainFrame frame,
                                             final TreePanel tp,
                                             final ControlPanel cp,
                                             final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        if ( !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "rate lost: the control must be there to press" );
            return;
        }
        final String unit_before = phy.getDistanceUnit();
        if ( !takeRateAway( phy, "isolate_C" ) ) {
            fail( ok, "rate lost: the rate was not taken away" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "a tree that has lost a rate must refuse the divergence view" );
        }
        if ( !near( sum( phy ), TIME_SUM ) || !near( length( phy, "isolate_C" ), 0.8 ) ) {
            fail( ok, "a refused switch must leave the lengths alone; got sum=" + sum( phy ) + " C="
                    + length( phy, "isolate_C" ) );
        }
        if ( !same( unit_before, phy.getDistanceUnit() ) ) {
            fail( ok, "a refused switch must leave the unit alone; was " + unit_before + ", is "
                    + phy.getDistanceUnit() );
        }
        if ( cp.isBranchLengthsControlVisible() ) {
            fail( ok, "a refused switch must take the control away" );
        }
        if ( !"Time".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "a refused switch must not leave Div pressed" );
        }
    }

    private static void rateLostWhileInDivergence( final MainFrame frame,
                                                   final TreePanel tp,
                                                   final ControlPanel cp,
                                                   final boolean[] ok ) {
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        final Phylogeny phy = tp.getPhylogeny();
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !near( sum( phy ), DIV_SUM ) ) {
            fail( ok, "stranded: the tree must BE in the divergence view first" );
            return;
        }
        if ( !takeRateAway( phy, "isolate_C" ) ) {
            fail( ok, "stranded: the rate was not taken away" );
            return;
        }
        // the panel learns of it (as after an edit, a tab change): the row is the way back, so it stays
        tp.invalidateBranchLengthToggle();
        cp.populateBranchLengthsControl();
        if ( tp.isBranchLengthToggleApplicable() ) {
            fail( ok, "stranded: the panel must know the tree has lost its second layout" );
        }
        if ( !cp.isBranchLengthsControlVisible() || !"Divergence".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "while divergence is on screen the control must stay, with Div pressed" );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "the way back to time needs only the dates: the tree must not be stranded in divergence" );
        }
        if ( !near( sum( phy ), TIME_SUM ) ) {
            fail( ok, "the way back must restore the time lengths, got sum=" + sum( phy ) );
        }
        if ( cp.isBranchLengthsControlVisible() ) {
            fail( ok, "back in time on a tree with no divergence layout: the control must go" );
        }
    }

    /** Divergence on screen, a rate gone, and the tree object REPLACED (what a redo does with its snapshot), so the
     *  panel knows the tree has no second layout and the switch is not there to press: Reset to Defaults is the way
     *  out, and it needs only the dates. */
    private static void resetAfterTreeReplaced( final MainFrame frame,
                                                final TreePanel tp,
                                                final ControlPanel cp,
                                                final boolean[] ok ) {
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || !takeRateAway( tp.getPhylogeny(), "isolate_C" ) ) {
            fail( ok, "replaced: the tree must be in the divergence view, and then lose a rate" );
            return;
        }
        final Phylogeny snapshot = tp.getPhylogeny().copy();
        tp.setTree( snapshot );
        if ( tp.isBranchLengthToggleApplicable() || !near( sum( snapshot ), DIV_SUM ) ) {
            fail( ok, "replaced: the new tree must be refused the switch while showing divergence; applicable="
                    + tp.isBranchLengthToggleApplicable() + " sum=" + sum( snapshot ) );
            return;
        }
        cp.populateBranchLengthsControl();
        if ( !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "replaced: divergence is on screen, so the control must be" );
        }
        tp.resetBranchLengthModeToDefault();
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "replaced: Reset must return the mode to TIME" );
        }
        if ( !near( sum( snapshot ), TIME_SUM ) ) {
            fail( ok, "replaced: Reset must restore the time lengths from the dates alone, got sum=" + sum( snapshot ) );
        }
    }

    /** Divergence on screen and the tree's DATES gone (its three inner nodes'): no branch has a time span left.
     *  The way back is still not refused -- time is laid out from what the dates say, which is nothing, and no
     *  branch keeps the divergence length it had. */
    private static boolean loseTheDates( final TreePanel tp, final ControlPanel cp, final boolean[] ok ) {
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        final Phylogeny phy = tp.getPhylogeny();
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || ( zeroBranches( phy ) != 0 ) ) {
            fail( ok, "dates lost: the tree must BE in the divergence view first, no branch at 0" );
            return false;
        }
        phy.getNode( "isolate_A" ).getParent().getNodeData().setDate( null );
        phy.getNode( "isolate_D" ).getParent().getNodeData().setDate( null );
        phy.getNode( "isolate_D" ).getParent().getParent().getNodeData().setDate( null );
        if ( BranchLengthLayout.isTimeDerivable( phy ) ) {
            fail( ok, "dates lost: the tree must have no time layout left" );
            return false;
        }
        return true;
    }

    private static void timeWithoutDates( final String what,
                                          final TreePanel tp,
                                          final ControlPanel cp,
                                          final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, what + " must return the mode to TIME even with the dates gone" );
        }
        if ( zeroBranches( phy ) != 8 ) {
            fail( ok, what + ": a branch without dates is laid out at 0, never at the divergence it had; "
                    + zeroBranches( phy ) + " of 8 at 0, sum=" + sum( phy ) );
        }
        if ( !"time".equals( phy.getDistanceUnit() ) ) {
            fail( ok, what + " must take the divergence unit away, got " + phy.getDistanceUnit() );
        }
    }

    private static void datesLostThenTime( final MainFrame frame,
                                           final TreePanel tp,
                                           final ControlPanel cp,
                                           final boolean[] ok ) {
        if ( !loseTheDates( tp, cp, ok ) ) {
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        timeWithoutDates( "the switch", tp, cp, ok );
        if ( cp.isBranchLengthsControlVisible() ) {
            fail( ok, "back in time on a tree with no second layout: the control must go" );
        }
    }

    private static void datesLostThenReset( final MainFrame frame,
                                            final TreePanel tp,
                                            final ControlPanel cp,
                                            final boolean[] ok ) {
        if ( !loseTheDates( tp, cp, ok ) ) {
            return;
        }
        tp.resetBranchLengthModeToDefault();
        timeWithoutDates( "Reset", tp, cp, ok );
    }

    /** Undo and redo bring back a TREE, and its branch lengths are in the layout they were captured in: the mode
     *  follows the tree. An edit made in time, the switch pressed, the edit undone -- time lengths are back, and
     *  the panel must say Time, not Div. */
    private static void undoAcrossTheSwitch( final MainFrame frame,
                                             final TreePanel tp,
                                             final ControlPanel cp,
                                             final boolean[] ok ) {
        tp.pushUndoCheckpoint( "Rename" );
        tp.getPhylogeny().getNode( "isolate_A" ).setName( "isolate_A_renamed" );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !near( sum( tp.getPhylogeny() ), DIV_SUM ) ) {
            fail( ok, "undo: the tree must BE in the divergence view first" );
            return;
        }
        if ( !tp.undo() ) {
            fail( ok, "undo: there must be something to undo" );
            return;
        }
        Phylogeny phy = tp.getPhylogeny();
        if ( !near( sum( phy ), TIME_SUM ) || !hasNode( phy, "isolate_A" ) ) {
            fail( ok, "undo: fixture -- the restored tree must be the one from before the edit, in time lengths" );
            return;
        }
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "undo brought time lengths back: the mode must be TIME" );
        }
        if ( !tp.isBranchLengthTimeCalibrated() ) {
            fail( ok, "undo brought time lengths back: the node-age bars must be allowed again" );
        }
        if ( !cp.isBranchLengthsControlVisible() || !"Time".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "undo brought time lengths back: the control must show Time pressed" );
        }
        if ( !tp.redo() ) {
            fail( ok, "redo: there must be something to redo" );
            return;
        }
        phy = tp.getPhylogeny();
        if ( !near( sum( phy ), DIV_SUM ) || !hasNode( phy, "isolate_A_renamed" ) ) {
            fail( ok, "redo: fixture -- the restored tree must be the edited one, in divergence lengths" );
            return;
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || tp.isBranchLengthTimeCalibrated() ) {
            fail( ok, "redo brought divergence lengths back: the mode must be DIVERGENCE" );
        }
        if ( !"Divergence".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "redo brought divergence lengths back: the control must show Div pressed" );
        }
        // and the button still works from there
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( sum( phy ), TIME_SUM ) ) {
            fail( ok, "after a redo into divergence, Time must still lay the tree out by time" );
        }
    }

    // ---- the subtree view. It SHARES its nodes with the tree it came from, so the switch pressed down there
    //      rewrites that clade up here, and nothing else: back on the whole tree, ONE layout must hold.

    private static PhylogenyNode cladeOf( final Phylogeny phy, final String tip ) {
        return phy.getNode( tip ).getParent();
    }

    /** Div pressed in the view of (C,(D,E)): back on the whole tree, ALL eight branches are in divergence -- not
     *  four of them beside four in time, under a unit that is blank. */
    private static void divergenceInsideASubtree( final MainFrame frame,
                                                  final TreePanel tp,
                                                  final ControlPanel cp,
                                                  final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        tp.subTree( cladeOf( whole, "isolate_C" ) );
        if ( !tp.isCurrentTreeIsSubtree() || ( tp.getPhylogeny().getNodeCount() != 5 ) ) {
            fail( ok, "subtree: the view must be of (C,(D,E)), got " + tp.getPhylogeny().getNodeCount() + " nodes" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || !near( length( whole, "isolate_C" ), 0.00208 ) || !near( length( whole, "isolate_A" ), 1.2 ) ) {
            fail( ok, "subtree: fixture -- Div in the view must rewrite the clade and only the clade; C="
                    + length( whole, "isolate_C" ) + " A=" + length( whole, "isolate_A" ) );
            return;
        }
        tp.superTree();
        if ( ( tp.getPhylogeny() != whole ) || tp.isCurrentTreeIsSubtree() ) {
            fail( ok, "subtree: the whole tree must be back on display" );
            return;
        }
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) {
            fail( ok, "back on a tree that can show divergence, the mode stays DIVERGENCE" );
        }
        if ( !near( sum( whole ), DIV_SUM ) || !near( length( whole, "isolate_A" ), 0.00372 ) ) {
            fail( ok, "back on the whole tree EVERY branch must be in divergence; A=" + length( whole, "isolate_A" )
                    + " C=" + length( whole, "isolate_C" ) + " sum=" + sum( whole ) );
        }
        if ( !"subs/site".equals( whole.getDistanceUnit() ) ) {
            fail( ok, "the whole tree must say its lengths are subs/site, got \"" + whole.getDistanceUnit() + "\"" );
        }
        if ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) {
            fail( ok, "divergence has no time axis, got " + tp.effectiveTimeAxisType() );
        }
        if ( !cp.isBranchLengthsControlVisible() || !"Divergence".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "the control must show Div pressed" );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( !near( sum( whole ), TIME_SUM ) ) {
            fail( ok, "and Time must bring every branch back, got sum=" + sum( whole ) );
        }
    }

    /** The other way round: the whole tree in divergence, Time pressed in the view of the clade. */
    private static void timeInsideASubtree( final MainFrame frame,
                                            final TreePanel tp,
                                            final ControlPanel cp,
                                            final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        tp.subTree( cladeOf( whole, "isolate_C" ) );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( length( whole, "isolate_C" ), 0.8 )
                || !near( length( whole, "isolate_A" ), 0.00372 ) ) {
            fail( ok, "subtree: fixture -- Time in the view must rewrite the clade and only the clade; C="
                    + length( whole, "isolate_C" ) + " A=" + length( whole, "isolate_A" ) );
            return;
        }
        tp.superTree();
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "back on the whole tree the mode stays TIME" );
        }
        if ( !near( sum( whole ), TIME_SUM ) || !near( length( whole, "isolate_A" ), 1.2 ) ) {
            fail( ok, "back on the whole tree EVERY branch must be in time; A=" + length( whole, "isolate_A" )
                    + " sum=" + sum( whole ) );
        }
        if ( !"time".equals( whole.getDistanceUnit() ) ) {
            fail( ok, "the whole tree must say its lengths are time, got \"" + whole.getDistanceUnit() + "\"" );
        }
        if ( !"Time".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "the control must show Time pressed" );
        }
    }

    /** A clade that can be laid out both ways inside a tree that cannot: (A,B) of the tree whose isolate_C has no
     *  rate. Div pressed in the view; back on the whole tree there is no divergence layout to show, so the clade
     *  goes back to time -- and NOTHING outside it is rewritten (isolate_D is given a length that is not its date
     *  gap, as a file may state one: it must survive). */
    private static void subtreeOfATreeWithoutTheSwitch( final MainFrame frame,
                                                        final TreePanel tp,
                                                        final ControlPanel cp,
                                                        final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        whole.getNode( "isolate_D" ).setDistanceToParent( 0.47 );
        final String unit_before = whole.getDistanceUnit();
        tp.subTree( cladeOf( whole, "isolate_A" ) );
        cp.populateBranchLengthsControl(); // what a change of tab does
        if ( !tp.isBranchLengthToggleApplicable() || !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "clade: (A,B) states dates and a rate on both branches, its view must be offered the switch" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || !near( length( whole, "isolate_A" ), 0.00372 ) ) {
            fail( ok, "clade: fixture -- Div must take in the view of (A,B)" );
            return;
        }
        tp.superTree();
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "back on a tree with no divergence layout the mode must be TIME" );
        }
        if ( !tp.isBranchLengthTimeCalibrated() ) {
            fail( ok, "back in time: the node-age bars must be allowed again" );
        }
        if ( !near( length( whole, "isolate_A" ), 1.2 ) || !near( length( whole, "isolate_B" ), 1.2 ) ) {
            fail( ok, "the clade must be back in time; A=" + length( whole, "isolate_A" ) + " B="
                    + length( whole, "isolate_B" ) );
        }
        if ( !near( length( whole, "isolate_D" ), 0.47 ) || !near( length( whole, "isolate_C" ), 0.8 )
                || !near( parentLength( whole, "isolate_A" ), 0.9 ) ) {
            fail( ok, "no branch OUTSIDE the clade may be rewritten; D=" + length( whole, "isolate_D" ) + " C="
                    + length( whole, "isolate_C" ) + " (A,B)=" + parentLength( whole, "isolate_A" ) );
        }
        if ( !same( unit_before, whole.getDistanceUnit() ) ) {
            fail( ok, "the whole tree's unit was never changed and must not be now; was " + unit_before + ", is "
                    + whole.getDistanceUnit() );
        }
        if ( cp.isBranchLengthsControlVisible() || !"Time".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "back on a tree with no divergence layout, showing time: the control must go" );
        }
    }

    /** Divergence on screen all the way, no button pressed in the view -- but a rate taken away down there. Back on
     *  the whole tree there is no divergence layout left to show: time, on every branch. */
    private static void rateLostInsideASubtree( final MainFrame frame,
                                                final TreePanel tp,
                                                final ControlPanel cp,
                                                final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        tp.subTree( cladeOf( whole, "isolate_C" ) );
        if ( !tp.isCurrentTreeIsSubtree() || !takeRateAway( whole, "isolate_D" ) ) {
            fail( ok, "rate lost below: the view must descend and the rate must go" );
            return;
        }
        tp.superTree();
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "back on a tree that has lost its divergence layout the mode must be TIME" );
        }
        if ( !near( sum( whole ), TIME_SUM ) ) {
            fail( ok, "...and every branch in time, got sum=" + sum( whole ) );
        }
        if ( cp.isBranchLengthsControlVisible() ) {
            fail( ok, "...and the control gone" );
        }
    }

    /** The clade-only way back, on CALENDAR time, where the axis can be seen to return: a cherry that can be laid out
     *  both ways inside a tree that cannot (a tip outside it has lost its rate). */
    private static void calendarCladeBackInTime( final MainFrame frame,
                                                 final TreePanel tp,
                                                 final ControlPanel cp,
                                                 final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        PhylogenyNode cherry = null;
        for( final PhylogenyNodeIterator it = whole.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isRoot() && !n.isExternal() && ( n.getNumberOfDescendants() == 2 ) && n.getChildNode( 0 ).isExternal()
                    && n.getChildNode( 1 ).isExternal() ) {
                cherry = n;
                break;
            }
        }
        if ( cherry == null ) {
            fail( ok, "calendar: " + CALENDAR + " must have a cherry" );
            return;
        }
        PhylogenyNode outside = null;
        for( final PhylogenyNode tip : whole.getExternalNodes() ) {
            if ( tip.getParent() != cherry ) {
                outside = tip;
            }
        }
        final PhylogenyNode in = cherry.getChildNode( 0 );
        final double in_before = in.getDistanceToParent();
        final double out_before = outside.getDistanceToParent();
        final double depth_before = tp.getMaxDistanceToRootForTest();
        if ( !takeRateAway( whole, outside.getName() ) || ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) ) {
            fail( ok, "calendar: the tree must lose its second layout and be on the calendar axis" );
            return;
        }
        tp.subTree( cherry );
        cp.populateBranchLengthsControl(); // what a change of tab does
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) || !( in.getDistanceToParent() < ( in_before / 10 ) ) ) {
            fail( ok, "calendar: fixture -- Div must take in the view of the cherry; length " + in.getDistanceToParent()
                    + " was " + in_before );
            return;
        }
        tp.superTree();
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "calendar: back on a tree with no divergence layout the mode must be TIME" );
        }
        if ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) {
            fail( ok, "back in time the calendar axis must be back, got " + tp.effectiveTimeAxisType() );
        }
        if ( Math.abs( in.getDistanceToParent() - in_before ) > 1e-6 ) {
            fail( ok, "the cherry must be back in years; " + in.getDistanceToParent() + " was " + in_before );
        }
        if ( outside.getDistanceToParent() != out_before ) {
            fail( ok, "no branch outside the cherry may be rewritten; " + outside.getDistanceToParent() + " was "
                    + out_before );
        }
        if ( Math.abs( tp.getMaxDistanceToRootForTest() - depth_before ) > 1e-6 ) {
            fail( ok, "the depth the tree is drawn to must be the whole tree's again; "
                    + tp.getMaxDistanceToRootForTest() + " was " + depth_before );
        }
    }

    /** A length the file stated stays the length the file stated: going in and out of a subtree must not replace
     *  it by a date gap, inside the clade visited (A) or outside it (C). */
    private static void navigationAloneRewritesNothing( final MainFrame frame,
                                                        final TreePanel tp,
                                                        final ControlPanel cp,
                                                        final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        whole.getNode( "isolate_A" ).setDistanceToParent( 1.17 ); // not its date gap (1.2)
        whole.getNode( "isolate_C" ).setDistanceToParent( 0.77 ); // not its date gap (0.8)
        tp.subTree( cladeOf( whole, "isolate_A" ) );
        if ( !tp.isCurrentTreeIsSubtree() ) {
            fail( ok, "navigation: the view must descend" );
            return;
        }
        tp.superTree();
        if ( ( tp.getPhylogeny() != whole ) || !near( length( whole, "isolate_A" ), 1.17 )
                || !near( length( whole, "isolate_C" ), 0.77 ) ) {
            fail( ok, "navigation alone must rewrite no branch length; A=" + length( whole, "isolate_A" ) + " C="
                    + length( whole, "isolate_C" ) );
        }
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "navigation alone must not change the mode" );
        }
    }

    /** Div pressed in the view of (C,(D,E)) and pressed back: the clade is in time again, as the switch lays time
     *  out. The branches OUTSIDE it were in time all along and are not rewritten -- isolate_A keeps the length the
     *  file stated, which is not its date gap. */
    private static void pressedAndPressedBackInASubtree( final MainFrame frame,
                                                         final TreePanel tp,
                                                         final ControlPanel cp,
                                                         final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        whole.getNode( "isolate_A" ).setDistanceToParent( 1.17 ); // not its date gap (1.2)
        tp.subTree( cladeOf( whole, "isolate_C" ) );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || !near( length( whole, "isolate_C" ), 0.00208 ) ) {
            fail( ok, "pressed back: fixture -- Div must take in the view first" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        tp.superTree();
        if ( ( tp.getPhylogeny() != whole ) || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) ) {
            fail( ok, "pressed back: the whole tree must be back, in TIME" );
        }
        if ( !near( length( whole, "isolate_C" ), 0.8 ) || !near( length( whole, "isolate_D" ), 0.5 ) ) {
            fail( ok, "the clade must be in time; C=" + length( whole, "isolate_C" ) + " D="
                    + length( whole, "isolate_D" ) );
        }
        if ( !near( length( whole, "isolate_A" ), 1.17 ) ) {
            fail( ok, "a branch outside the clade was in time all along and must keep its length; A="
                    + length( whole, "isolate_A" ) );
        }
        if ( !cp.isBranchLengthsControlVisible() || !"Time".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "the tree still has both layouts: the control stays, Time pressed" );
        }
    }

    /** Two levels down, and the layout each level was LEFT in is what counts on the way back. The whole tree is left
     *  in divergence; in the view of (C,(D,E)) Time is pressed; one level further down a rate is taken away, so no
     *  tree on the way up has a second layout any more. The middle view was left in time: nothing to do there. The
     *  whole tree was left in DIVERGENCE: its branches outside the clade still are, and must go back to time. */
    private static void descendedInDivergence( final MainFrame frame,
                                               final TreePanel tp,
                                               final ControlPanel cp,
                                               final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        tp.subTree( cladeOf( whole, "isolate_C" ) );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        tp.subTree( cladeOf( whole, "isolate_D" ) );
        if ( ( tp.getPhylogeny().getNodeCount() != 3 ) || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME )
                || !near( length( whole, "isolate_A" ), 0.00372 ) || !near( length( whole, "isolate_D" ), 0.5 ) ) {
            fail( ok, "two levels: fixture -- the view of (D,E), in time, under a whole tree left in divergence" );
            return;
        }
        if ( !takeRateAway( whole, "isolate_D" ) ) {
            fail( ok, "two levels: the rate was not taken away" );
            return;
        }
        tp.superTree();
        if ( ( tp.getPhylogeny().getNodeCount() != 5 ) || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME )
                || !near( length( whole, "isolate_C" ), 0.8 ) || !near( length( whole, "isolate_D" ), 0.5 )
                || !near( parentLength( whole, "isolate_D" ), 0.3 ) ) {
            fail( ok, "the middle view was left in time and must still be; C=" + length( whole, "isolate_C" )
                    + " D=" + length( whole, "isolate_D" ) + " (D,E)=" + parentLength( whole, "isolate_D" ) );
        }
        if ( !near( length( whole, "isolate_A" ), 0.00372 ) ) {
            fail( ok, "two levels: fixture -- the whole tree's other clade must still be in divergence here; A="
                    + length( whole, "isolate_A" ) );
            return;
        }
        tp.superTree();
        if ( ( tp.getPhylogeny() != whole ) || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) ) {
            fail( ok, "two levels: the whole tree must be back, in TIME" );
        }
        if ( !near( sum( whole ), TIME_SUM ) || !near( length( whole, "isolate_A" ), 1.2 ) ) {
            fail( ok, "the whole tree was left in divergence: its branches outside the clade must go back to time; A="
                    + length( whole, "isolate_A" ) + " sum=" + sum( whole ) );
        }
        if ( !"time".equals( whole.getDistanceUnit() ) ) {
            fail( ok, "the whole tree must no longer say subs/site, got \"" + whole.getDistanceUnit() + "\"" );
        }
        if ( cp.isBranchLengthsControlVisible() ) {
            fail( ok, "in time on a tree that has lost its second layout: the control must go" );
        }
    }

    /** The other direction: an edit made while DIVERGENCE is on screen is captured in divergence, so undoing it
     *  from the time view brings divergence back -- lengths and mode. */
    private static void undoIntoDivergence( final TreePanel tp,
                                            final ControlPanel cp,
                                            final boolean[] ok,
                                            final boolean pre_captured ) {
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( pre_captured ) {
            // the other way in: a state captured first and filed once the operation is known to have changed
            // something (an import does this)
            tp.pushUndoSnapshot( tp.getPhylogeny().copy(), false, "Rename" );
        }
        else {
            tp.pushUndoCheckpoint( "Rename" );
        }
        tp.getPhylogeny().getNode( "isolate_B" ).setName( "isolate_B_renamed" );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( sum( tp.getPhylogeny() ), TIME_SUM ) ) {
            fail( ok, "undo into divergence: the tree must BE in the time view first" );
            return;
        }
        if ( !tp.undo() ) {
            fail( ok, "undo into divergence: there must be something to undo" );
            return;
        }
        final Phylogeny phy = tp.getPhylogeny();
        if ( !near( sum( phy ), DIV_SUM ) || !hasNode( phy, "isolate_B" ) ) {
            fail( ok, "undo into divergence: fixture -- the restored tree must be in divergence lengths" );
            return;
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || tp.isBranchLengthTimeCalibrated() ) {
            fail( ok, "undo brought divergence lengths back: the mode must be DIVERGENCE" );
        }
        if ( !"Divergence".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "undo brought divergence lengths back: the control must show Div pressed" );
        }
    }

    /** Redo puts the tree it REPLACES on the undo stack, in the layout that tree is in at that moment -- here
     *  divergence, switched to between the undo and the redo. */
    private static void redoFilesTheLayoutItLeaves( final MainFrame frame,
                                                    final TreePanel tp,
                                                    final ControlPanel cp,
                                                    final boolean[] ok ) {
        tp.pushUndoCheckpoint( "Rename" );
        tp.getPhylogeny().getNode( "isolate_E" ).setName( "isolate_E_renamed" );
        if ( !tp.undo() ) {
            fail( ok, "redo files: there must be something to undo" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE ); // the un-edited tree, in divergence
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !tp.canRedo() ) {
            fail( ok, "redo files: the switch must take, and must leave the redo in place" );
            return;
        }
        if ( !tp.redo() ) {
            fail( ok, "redo files: there must be something to redo" );
            return;
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( sum( tp.getPhylogeny() ), TIME_SUM )
                || !hasNode( tp.getPhylogeny(), "isolate_E_renamed" ) ) {
            fail( ok, "redo brought back the edited tree, captured in time: mode " + tp.getBranchLengthMode()
                    + ", sum " + sum( tp.getPhylogeny() ) );
            return;
        }
        if ( !tp.undo() ) {
            fail( ok, "redo files: the redo must be undoable" );
            return;
        }
        if ( !near( sum( tp.getPhylogeny() ), DIV_SUM ) || !hasNode( tp.getPhylogeny(), "isolate_E" ) ) {
            fail( ok, "redo files: fixture -- the tree redo replaced was in divergence lengths" );
            return;
        }
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) {
            fail( ok, "the tree redo replaced was showing divergence: undoing the redo must say so" );
        }
    }

    /** The time axis follows too. On calendar time the time view has a CALENDAR axis and divergence has none. */
    private static void undoBringsTheAxisBack( final MainFrame frame,
                                               final TreePanel tp,
                                               final ControlPanel cp,
                                               final boolean[] ok ) {
        if ( !tp.isBranchLengthToggleApplicable() || ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) ) {
            fail( ok, CALENDAR + " must be offered the switch and open on the calendar axis, got "
                    + tp.effectiveTimeAxisType() );
            return;
        }
        tp.pushUndoCheckpoint( "Edit" );
        tp.getPhylogeny().setDescription( "edited" );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) {
            fail( ok, "axis: divergence must have no time axis first" );
            return;
        }
        if ( !tp.undo() ) {
            fail( ok, "axis: there must be something to undo" );
            return;
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME )
                || ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) ) {
            fail( ok, "undo brought time lengths back: the calendar axis must come back with them, got "
                    + tp.effectiveTimeAxisType() );
        }
        if ( !tp.redo() ) {
            fail( ok, "axis: there must be something to redo" );
            return;
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) ) {
            fail( ok, "redo brought divergence lengths back: the time axis must go with them, got "
                    + tp.effectiveTimeAxisType() );
        }
    }

    /** A DATE edit decides too, in both directions and at once. The date value of ONE inner node cleared in the
     *  editor, and the time layout cannot state three of the eight branches (its own, and its two tips'): the
     *  control goes with that Write -- five of eight was a majority once, and was enough. Written back, it returns. */
    private static void editorTakesTheDates( final MainFrame frame,
                                             final TreePanel tp,
                                             final ControlPanel cp,
                                             final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        if ( !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "editor dates: the control must be there first" );
            return;
        }
        final PhylogenyNode ab = phy.getNode( "isolate_A" ).getParent();
        final NodeDataForm form = new NodeDataForm( ab, tp, NodeDataForm.Mode.EDIT );
        form.setTextForTest( NodeDataDraft.DATE_VALUE, "" );
        if ( !form.write() || ( ab.getNodeData().getDate() != null && ab.getNodeData().getDate().getValue() != null ) ) {
            fail( ok, "editor dates: the date value was not cleared; problems: " + form.problems() );
            return;
        }
        if ( BranchLengthLayout.isTimeDerivable( phy ) ) {
            fail( ok, "one node undated: the time layout cannot state every branch" );
        }
        if ( tp.isBranchLengthToggleApplicable() || cp.isBranchLengthsControlVisible() ) {
            fail( ok, "one node undated: the control must go with the Write" );
        }
        form.setTextForTest( NodeDataDraft.DATE_VALUE, "1.2" );
        if ( !form.write() || !BranchLengthLayout.isTimeDerivable( phy ) ) {
            fail( ok, "editor dates: the date value was not written back; problems: " + form.problems() );
            return;
        }
        if ( !tp.isBranchLengthToggleApplicable() || !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "every node dated again: the control must return with the Write" );
        }
    }

    /** The node editor decides, in both directions and at once: the missing rate written onto isolate_C, and the
     *  control is there; taken off again, and it is gone. Driven through the editor's own Write. */
    private static void editorDecides( final MainFrame frame,
                                       final TreePanel tp,
                                       final ControlPanel cp,
                                       final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        if ( tp.isBranchLengthToggleApplicable() || cp.isBranchLengthsControlVisible() ) {
            fail( ok, "editor: " + ONE_UNRATED + " must start without the control" );
            return;
        }
        final NodeDataForm form = new NodeDataForm( phy.getNode( "isolate_C" ), tp, NodeDataForm.Mode.EDIT );
        form.addPropertyForTest( new NodeDataDraft.PropertyDraft( RATE, "0.0026", "", "xsd:decimal",
                                                                  org.forester.phylogeny.data.Property.AppliesTo.NODE ) );
        if ( !form.write() || ( ratedBranches( phy ) != 8 ) ) {
            fail( ok, "editor: the rate was not written; problems: " + form.problems() );
            return;
        }
        if ( !tp.isBranchLengthToggleApplicable() ) {
            fail( ok, "a rate on every branch now: the panel must know at once" );
        }
        if ( !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "a rate on every branch now: the control must appear with the Write" );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !near( sum( phy ), DIV_SUM ) ) {
            fail( ok, "the completed tree must lay out like its twin, got sum=" + sum( phy ) );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        final javax.swing.JTable table = form.propertyTableForTest();
        int row = -1;
        for( int i = 0; i < table.getRowCount(); ++i ) {
            if ( RATE.equals( table.getValueAt( i, 0 ) ) ) {
                row = i;
            }
        }
        if ( row < 0 ) {
            fail( ok, "editor: the form must list the rate it wrote" );
            return;
        }
        form.removePropertyForTest( row );
        if ( !form.write() || ( ratedBranches( phy ) != 7 ) ) {
            fail( ok, "editor: the rate was not taken off again" );
            return;
        }
        if ( tp.isBranchLengthToggleApplicable() || cp.isBranchLengthsControlVisible() ) {
            fail( ok, "a branch without a rate again: the control must go with the Write" );
        }
    }

    private static void notOffered( final MainFrame frame,
                                    final TreePanel tp,
                                    final ControlPanel cp,
                                    final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        if ( !BranchLengthLayout.isTimeDerivable( phy ) ) {
            fail( ok, ONE_UNRATED + " must have a time layout: what it lacks is a rate" );
        }
        if ( tp.isBranchLengthToggleApplicable() ) {
            fail( ok, ONE_UNRATED + " has a branch without a rate: the switch must not apply" );
        }
        if ( cp.isBranchLengthsControlVisible() ) {
            fail( ok, "the Time | Div control must not be shown for " + ONE_UNRATED );
        }
        final String unit_before = phy.getDistanceUnit();
        tp.setBranchLengthMode( BranchLengthLayout.MODE.DIVERGENCE );
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, "asked directly, the panel must still refuse the divergence view" );
        }
        if ( !near( sum( phy ), TIME_SUM ) || ( zeroBranches( phy ) != 0 ) ) {
            fail( ok, "a refused tree keeps the lengths its file states; got sum=" + sum( phy ) + ", "
                    + zeroBranches( phy ) + " at 0" );
        }
        if ( !same( unit_before, phy.getDistanceUnit() ) ) {
            fail( ok, "a refused tree keeps its unit; was " + unit_before + ", is " + phy.getDistanceUnit() );
        }
    }

    /** A tree one of the two layouts cannot state every branch of: no control, the panel refuses when asked
     *  directly, and the tree stays in the layout it arrived in -- every length, and the unit. */
    private static void refused( final String demo, final TreePanel tp, final ControlPanel cp, final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        final java.util.List<Double> arrived = lengths( phy );
        final String unit_before = phy.getDistanceUnit();
        if ( tp.isBranchLengthToggleApplicable() ) {
            fail( ok, demo + ": a layout cannot state every branch, the switch must not apply" );
        }
        if ( cp.isBranchLengthsControlVisible() ) {
            fail( ok, demo + ": the Time | Div control must not be shown" );
        }
        tp.setBranchLengthMode( BranchLengthLayout.MODE.DIVERGENCE );
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
            fail( ok, demo + ": asked directly, the panel must still refuse the divergence view" );
        }
        if ( !arrived.equals( lengths( phy ) ) ) {
            fail( ok, demo + ": the tree must stay in the layout it arrived in; " + arrived + " -> " + lengths( phy ) );
        }
        if ( !same( unit_before, phy.getDistanceUnit() ) ) {
            fail( ok, demo + ": a refused tree keeps its unit; was " + unit_before + ", is " + phy.getDistanceUnit() );
        }
    }

    /** ...and the twin that records a divergence on EVERY node is offered the switch, or the refusal proves nothing. */
    private static void recordedTwinIsOffered( final MainFrame frame,
                                               final TreePanel tp,
                                               final ControlPanel cp,
                                               final boolean[] ok ) {
        if ( !tp.isBranchLengthToggleApplicable() || !cp.isBranchLengthsControlVisible() ) {
            fail( ok, RECORDED + " records a divergence and a date on every node: it must be offered the switch" );
        }
        final String tip = cp.branchLengthDivTooltipForTest();
        if ( ( tip == null ) || !tip.contains( "recorded in the file" ) ) {
            fail( ok, RECORDED + ": the Div tooltip must say the divergence is recorded, got: " + tip );
        }
    }

    /** The (D,E) node is dated 0.05 BEFORE its parent. Time keeps the sign, so the lengths from the root add up to
     *  each node's own date; Div states 0; and back in time it is -0.05 again. About the VALUES: the window draws a
     *  negative length at 0, as it always has. */
    private static void aSpanThatRunsBackwards( final MainFrame frame,
                                                final TreePanel tp,
                                                final ControlPanel cp,
                                                final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        if ( !cp.isBranchLengthsControlVisible() || !near( parentLength( phy, "isolate_D" ), -0.05 ) ) {
            fail( ok, BACKWARDS + " must be offered the switch and open with (D,E) at -0.05, got "
                    + parentLength( phy, "isolate_D" ) );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) {
            fail( ok, "backwards: selecting Div must switch the mode" );
            return;
        }
        if ( Double.doubleToRawLongBits( parentLength( phy, "isolate_D" ) ) != 0L ) {
            fail( ok, "divergence draws a span that runs backwards at 0, got " + parentLength( phy, "isolate_D" ) );
        }
        // every other branch at rate x span: D is 0.85 x 0.0035, C is 0.8 x 0.0026
        if ( !near( length( phy, "isolate_D" ), 0.85 * 0.0035 ) || !near( length( phy, "isolate_C" ), 0.00208 )
                || ( zeroBranches( phy ) != 1 ) ) {
            fail( ok, "...and only that one; D=" + length( phy, "isolate_D" ) + " C=" + length( phy, "isolate_C" ) + ", "
                    + zeroBranches( phy ) + " branches at 0" );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( !near( parentLength( phy, "isolate_D" ), -0.05 ) || !near( length( phy, "isolate_D" ), 0.85 ) ) {
            fail( ok, "back in time the span is -0.05 again; got (D,E)=" + parentLength( phy, "isolate_D" ) + " D="
                    + length( phy, "isolate_D" ) );
        }
        // every tip is dated 0 under a root dated 2.1: the signed lengths from the root add up to 2.1
        // (NOT calculateDistanceToRoot(): it leaves negative lengths out of the sum)
        for( final PhylogenyNode tip : phy.getExternalNodes() ) {
            double signed = 0;
            for( PhylogenyNode n = tip; !n.isRoot(); n = n.getParent() ) {
                signed += n.getDistanceToParent();
            }
            if ( !near( signed, 2.1 ) ) {
                fail( ok, "the lengths from the root to " + tip.getName() + " must add up to 2.1, got " + signed );
            }
        }
    }

    private static java.util.List<Double> lengths( final Phylogeny phy ) {
        final java.util.List<Double> l = new java.util.ArrayList<>();
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( !node.isRoot() ) {
                l.add( node.getDistanceToParent() );
            }
        }
        return l;
    }

    /** Takes the clock rate off the named node IN PLACE, as the node editor would; true when the tree then has no
     *  divergence source (the drive took). */
    private static boolean takeRateAway( final Phylogeny phy, final String name ) {
        final PhylogenyNode n = phy.getNode( name );
        if ( ( n.getNodeData().getProperties() == null )
                || !n.getNodeData().getProperties().getProperties().removeIf( p -> RATE.equals( p.getRef() ) ) ) {
            return false;
        }
        return BranchLengthLayout.divergenceSource( phy ) == BranchLengthLayout.DIVERGENCE_SOURCE.NONE;
    }

    private static boolean hasNode( final Phylogeny phy, final String name ) {
        return phy.getNodes( name ).size() == 1;
    }

    private static int branches( final Phylogeny phy ) {
        return phy.getNodeCount() - 1;
    }

    private static int ratedBranches( final Phylogeny phy ) {
        int n = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( !node.isRoot() && ( node.getNodeData().getProperties() != null )
                    && !node.getNodeData().getProperties().getProperties( RATE ).isEmpty() ) {
                ++n;
            }
        }
        return n;
    }

    private static int zeroBranches( final Phylogeny phy ) {
        int n = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( !node.isRoot() && ( node.getDistanceToParent() == 0.0 ) ) {
                ++n;
            }
        }
        return n;
    }

    private static double sum( final Phylogeny phy ) {
        double s = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( !node.isRoot() ) {
                s += node.getDistanceToParent();
            }
        }
        return s;
    }

    private static double length( final Phylogeny phy, final String name ) {
        return phy.getNode( name ).getDistanceToParent();
    }

    private static double parentLength( final Phylogeny phy, final String name ) {
        return phy.getNode( name ).getParent().getDistanceToParent();
    }

    private static boolean near( final double a, final double b ) {
        return Math.abs( a - b ) < 1e-12;
    }

    private static boolean same( final String a, final String b ) {
        return ( a == null ) ? ( b == null ) : a.equals( b );
    }

    /** A demo read exactly as File > Open reads it. */
    private static Phylogeny read( final String demo ) {
        final File f = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
        if ( !f.exists() ) {
            return null;
        }
        try {
            final Phylogeny[] phys = FigureRenderer.readTrees( f );
            return ( phys.length == 1 ) ? phys[ 0 ] : null;
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [BeastBranchModeTest] " + m );
        return false;
    }

    private static void fail( final boolean[] ok, final String m ) {
        ok[ 0 ] = false;
        fail( m );
    }

    private BeastBranchModeTest() {
    }
}
