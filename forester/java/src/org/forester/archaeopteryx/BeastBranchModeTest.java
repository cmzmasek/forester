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
    /** {@link #RATED} with isolate_A stating 1.407, its dates 1.2 apart: a file's length that is not its date gap. */
    private static final String NOT_HEIGHTS = "beast-lengths-not-heights.nex";
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
                    && inFrame( RATED, BeastBranchModeTest::theFilesLengthNotItsDateGap )
                    && inFrame( RATED, BeastBranchModeTest::editorTakesALength )
                    && inFrame( RATED, BeastBranchModeTest::rateLostWhileInTime )
                    && inFrame( RATED, BeastBranchModeTest::rateLostWhileInDivergence )
                    && inFrame( RATED, BeastBranchModeTest::resetAfterTreeReplaced )
                    && inFrame( RATED, BeastBranchModeTest::datesLostThenTime )
                    && inFrame( RATED, BeastBranchModeTest::datesLostThenReset )
                    && inFrame( RATED, BeastBranchModeTest::undoAcrossTheSwitch )
                    && inFrame( RATED, ( f, tp, cp, ok ) -> undoIntoDivergence( tp, cp, ok, false ) )
                    && inFrame( RATED, ( f, tp, cp, ok ) -> undoIntoDivergence( tp, cp, ok, true ) )
                    && inFrame( RATED, BeastBranchModeTest::redoFilesTheLayoutItLeaves )
                    && inFrame( RATED, BeastBranchModeTest::undoIntoDivergenceTakesItsOwnTime )
                    && inFrame( RATED, BeastBranchModeTest::redoIntoDivergenceTakesItsOwnTime )
                    && inFrame( RATED, BeastBranchModeTest::redoFilesTheTimeOfTheTreeItLeaves )
                    && inFrame( RATED, BeastBranchModeTest::deleteNodeInDivergence )
                    && inFrame( RATED, BeastBranchModeTest::deleteSubtreeInDivergence )
                    && inFrame( RATED, BeastBranchModeTest::deleteNodeInTime )
                    && inFrame( NOT_HEIGHTS, BeastBranchModeTest::aTabMadeFromATimeTab )
                    && inFrame( NOT_HEIGHTS, BeastBranchModeTest::aTabMadeFromADivergenceTab )
                    && aBackgroundTabArrivesBeforeItsAxisIsRead()
                    && inFrame( CALENDAR, BeastBranchModeTest::undoBringsTheAxisBack )
                    && inFrame( RATED, BeastBranchModeTest::divergenceInsideASubtree )
                    && inFrame( RATED, BeastBranchModeTest::timeInsideASubtree )
                    && inFrame( ONE_UNRATED, BeastBranchModeTest::subtreeOfATreeWithoutTheSwitch )
                    && inFrame( RATED, BeastBranchModeTest::navigationAloneRewritesNothing )
                    && inFrame( RATED, BeastBranchModeTest::twoLevelsDown )
                    && inFrame( RATED, BeastBranchModeTest::rateLostInsideASubtree )
                    && inFrame( BACKWARDS, BeastBranchModeTest::aViewRunsTheWayItsTreeDoes )
                    && inFrame( CALENDAR, BeastBranchModeTest::aCherryOfACalendarTree )
                    && inFrame( RATED, BeastBranchModeTest::undoInsideASubtree )
                    && inFrame( RATED, BeastBranchModeTest::undoIntoDivergenceInsideASubtree )
                    && inFrame( RATED, BeastBranchModeTest::aCopyOfACladeWithoutChange )
                    && inFrame( ONE_UNRATED, BeastBranchModeTest::editorDecides )
                    && inFrame( RATED, BeastBranchModeTest::editorTakesTheDates )
                    && inFrame( ONE_UNRATED, BeastBranchModeTest::notOffered )
                    && inFrame( ONE_UNDATED, ( f, tp, cp, ok ) -> refused( ONE_UNDATED, tp, cp, ok ) )
                    && inFrame( ONE_WITHOUT_DIV, ( f, tp, cp, ok ) -> refused( ONE_WITHOUT_DIV, tp, cp, ok ) )
                    && inFrame( "beast-length-missing.nex", ( f, tp, cp, ok ) -> refused( "beast-length-missing.nex", tp, cp, ok ) )
                    && inFrame( "beast-rates-zero.nex", ( f, tp, cp, ok ) -> refused( "beast-rates-zero.nex", tp, cp, ok ) )
                    && inFrame( "beast-rate-spelling.nex", ( f, tp, cp, ok ) -> refused( "beast-rate-spelling.nex", tp, cp, ok ) )
                    && arrivesShowingDivergence( RATED, TIME_SUM, DIV_SUM )
                    && arrivesShowingDivergence( RECORDED, 22.45, 0.06735 )
                    && aViewOfATreeThatKeptNothing()
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
        if ( ( tip == null ) || !tip.contains( "DERIVED" ) || !tip.contains( "its length in time" ) ) {
            fail( ok, "the Div tooltip must say the divergence is DERIVED from the rate and the branch's length in "
                    + "time, got: " + tip );
        }
        // what the button SAYS it does is what theFilesLengthNotItsDateGap measures it doing
        final String time_tip = cp.branchLengthTimeTooltipForTest();
        if ( ( time_tip == null ) || !time_tip.contains( "stated length in time" )
                || time_tip.contains( "difference between" ) ) {
            fail( ok, "the Time tooltip must say Time is the lengths the tree states, and must not say it is the "
                    + "difference between the dates, got: " + time_tip );
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !"Time".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "the default mode must be TIME" );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        final Phylogeny phy = tp.getPhylogeny();
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) {
            fail( ok, "selecting Div must switch the mode" );
        }
        if ( !"Divergence".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "the control must show Div pressed" );
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
        // the deepest tip in divergence is D: 1.3 x 0.0029 + 0.3 x 0.0034 + 0.5 x 0.0035 (A is at 0.00642, E at
        // 0.00644, C at 0.00585) -- in time all five are 2.1 deep
        if ( Math.abs( tp.getMaxDistanceToRootForTest() - 0.00654 ) > 1e-9 ) {
            fail( ok, "the depth the tree is drawn to must be the divergence tree's, 0.00654; got "
                    + tp.getMaxDistanceToRootForTest() );
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

    /** Where the length a file states is not the gap between its two dates (one summary tree states 1.407 between
     *  nodes dated 1.098 apart): divergence is the rate x the STATED length, and Time gives the stated length back
     *  -- the tree that was opened, not another one. By the switch and by Reset. */
    private static void theFilesLengthNotItsDateGap( final MainFrame frame,
                                                     final TreePanel tp,
                                                     final ControlPanel cp,
                                                     final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        phy.getNode( "isolate_A" ).setDistanceToParent( 1.407 ); // its dates are 1.2 apart
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( !near( length( phy, "isolate_A" ), 1.407 * 0.0031 ) || !near( length( phy, "isolate_B" ), 0.00336 ) ) {
            fail( ok, "divergence is the rate x the length the file states; got A=" + length( phy, "isolate_A" ) + " B="
                    + length( phy, "isolate_B" ) );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( length( phy, "isolate_A" ) != 1.407 ) || ( length( phy, "isolate_B" ) != 1.2 ) ) {
            fail( ok, "Time gives back the length the file states, to the digit; got A=" + length( phy, "isolate_A" )
                    + " B=" + length( phy, "isolate_B" ) );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        tp.resetBranchLengthModeToDefault();
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || ( length( phy, "isolate_A" ) != 1.407 ) ) {
            fail( ok, "...and so does Reset; got A=" + length( phy, "isolate_A" ) );
        }
        // a length changed WHILE time is on screen is the length the tree has in time from then on
        phy.getNode( "isolate_A" ).setDistanceToParent( 1.3 );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( length( phy, "isolate_A" ) != 1.3 ) {
            fail( ok, "what is kept is what was on screen when the tree LEFT time, each time; got A="
                    + length( phy, "isolate_A" ) );
        }
    }

    /** A branch LENGTH decides too: taken out in the editor, the time layout cannot state that branch and the
     *  control goes with the Write; written back, it returns. */
    private static void editorTakesALength( final MainFrame frame,
                                            final TreePanel tp,
                                            final ControlPanel cp,
                                            final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        if ( !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "editor length: the control must be there first" );
            return;
        }
        final PhylogenyNode c = phy.getNode( "isolate_C" );
        final NodeDataForm form = new NodeDataForm( c, tp, NodeDataForm.Mode.EDIT );
        form.setTextForTest( NodeDataDraft.BRANCH_LENGTH, "" );
        if ( !form.write() || ( c.getDistanceToParent() != org.forester.phylogeny.data.PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) ) {
            fail( ok, "editor length: the length was not taken out (it is " + c.getDistanceToParent() + "); problems: "
                    + form.problems() );
            return;
        }
        if ( tp.isBranchLengthToggleApplicable() || cp.isBranchLengthsControlVisible() ) {
            fail( ok, "one branch without a length: the control must go with the Write" );
        }
        form.setTextForTest( NodeDataDraft.BRANCH_LENGTH, "0.8" );
        if ( !form.write() || ( c.getDistanceToParent() != 0.8 ) ) {
            fail( ok, "editor length: the length was not written back; problems: " + form.problems() );
            return;
        }
        if ( !tp.isBranchLengthToggleApplicable() || !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "every branch with a length again: the control must return with the Write" );
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
        if ( Math.abs( tp.getMaxDistanceToRootForTest() - 2.1 ) > 1e-9 ) {
            fail( ok, "after Reset the depth the tree is drawn to must be the time tree's, 2.1; got "
                    + tp.getMaxDistanceToRootForTest() );
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

    /** Divergence on screen and the tree's DATES gone (its three inner nodes'): it has no second layout any more.
     *  The way back is still not refused, and it gives back the lengths the tree had in time -- they were kept when
     *  it left time, and are not read off the dates. */
    private static boolean loseTheDates( final TreePanel tp, final ControlPanel cp, final boolean[] ok ) {
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        final Phylogeny phy = tp.getPhylogeny();
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !near( sum( phy ), DIV_SUM ) ) {
            fail( ok, "dates lost: the tree must BE in the divergence view first" );
            return false;
        }
        phy.getNode( "isolate_A" ).getParent().getNodeData().setDate( null );
        phy.getNode( "isolate_D" ).getParent().getNodeData().setDate( null );
        phy.getNode( "isolate_D" ).getParent().getParent().getNodeData().setDate( null );
        tp.invalidateBranchLengthToggle();
        if ( tp.isBranchLengthToggleApplicable() ) {
            fail( ok, "dates lost: the tree must have no second layout left" );
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
        if ( !near( sum( phy ), TIME_SUM ) || !near( length( phy, "isolate_A" ), 1.2 ) || ( zeroBranches( phy ) != 0 ) ) {
            fail( ok, what + " must give back the lengths the tree had in time; sum=" + sum( phy ) + ", "
                    + zeroBranches( phy ) + " of 8 at 0" );
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
        tp.getPhylogeny().getNode( "isolate_B" ).setDistanceToParent( 1.17 ); // as a file may state it: not its date gap, 1.2
        // ...so B's divergence is 1.17 x 0.0028, and the eight add up to that much less
        final double DIV_SUM_B = DIV_SUM + ( 0.1 * 0.0028 ); // the edited tree: B is 1.3 long, not 1.2
        tp.pushUndoCheckpoint( "Rename" );
        tp.getPhylogeny().getNode( "isolate_A" ).setName( "isolate_A_renamed" );
        // the edit also gives B a length of 1.3: what the tab KEEPS when it leaves time is 1.3, what the snapshot
        // holds is 1.17 -- a snapshot laid out again from what was kept since would come back with 1.3
        tp.getPhylogeny().getNode( "isolate_B" ).setDistanceToParent( 1.3 );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !near( sum( tp.getPhylogeny() ), DIV_SUM_B ) ) {
            fail( ok, "undo: the tree must BE in the divergence view first" );
            return;
        }
        if ( !tp.undo() ) {
            fail( ok, "undo: there must be something to undo" );
            return;
        }
        Phylogeny phy = tp.getPhylogeny();
        if ( !hasNode( phy, "isolate_A" ) || !near( length( phy, "isolate_A" ), 1.2 ) ) {
            fail( ok, "undo: fixture -- the restored tree must be the one from before the edit, in time lengths" );
            return;
        }
        if ( !near( length( phy, "isolate_B" ), 1.17 ) || !near( sum( phy ), TIME_SUM - 0.03 ) ) {
            fail( ok, "the restored tree is left AS IT WAS CAPTURED: B stated 1.17, not its date gap; got B="
                    + length( phy, "isolate_B" ) + " sum=" + sum( phy ) );
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
        if ( !near( sum( phy ), DIV_SUM_B ) || !hasNode( phy, "isolate_A_renamed" ) ) {
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
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || ( length( phy, "isolate_B" ) != 1.3 )
                || !near( sum( phy ), TIME_SUM + 0.1 ) ) {
            fail( ok, "after a redo into divergence, Time must give back the lengths the EDITED tree had in time; B="
                    + length( phy, "isolate_B" ) + " sum=" + sum( phy ) );
        }
    }

    // ---- the subtree view. ONE TREE, ONE LAYOUT: the switch acts on the WHOLE tree a tab holds and is judged on
    //      the whole tree, also while a subtree of it is on view (as in Archaeopteryx.js). A view shares its nodes
    //      with the tree it is a view of; a switch that acted on the view alone left the whole tree in two units.

    private static PhylogenyNode cladeOf( final Phylogeny phy, final String tip ) {
        return phy.getNode( tip ).getParent();
    }

    /** Div pressed in the view of (C,(D,E)): ALL eight branches of the whole tree are in divergence AT ONCE -- the
     *  four outside the view too -- and both trees say subs/site. Going back up rewrites nothing more. */
    private static void divergenceInsideASubtree( final MainFrame frame,
                                                  final TreePanel tp,
                                                  final ControlPanel cp,
                                                  final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        tp.subTree( cladeOf( whole, "isolate_C" ) );
        final Phylogeny view = tp.getPhylogeny();
        if ( !tp.isCurrentTreeIsSubtree() || ( view.getNodeCount() != 5 ) ) {
            fail( ok, "subtree: the view must be of (C,(D,E)), got " + view.getNodeCount() + " nodes" );
            return;
        }
        if ( !tp.isBranchLengthToggleApplicable() || !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "the whole tree has both layouts: the view of its clade is offered the switch" );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) {
            fail( ok, "subtree: Div must take in the view" );
            return;
        }
        if ( !near( length( whole, "isolate_C" ), 0.00208 ) || !near( length( whole, "isolate_A" ), 0.00372 )
                || !near( sum( whole ), DIV_SUM ) ) {
            fail( ok, "pressed in the view, the WHOLE tree is in divergence, the branches outside the view too; C="
                    + length( whole, "isolate_C" ) + " A=" + length( whole, "isolate_A" ) + " sum=" + sum( whole ) );
        }
        if ( !"subs/site".equals( whole.getDistanceUnit() ) || !"subs/site".equals( view.getDistanceUnit() ) ) {
            fail( ok, "the whole tree and the view must both say subs/site; got \"" + whole.getDistanceUnit()
                    + "\" and \"" + view.getDistanceUnit() + "\"" );
        }
        final java.util.List<Double> pressed = lengths( whole );
        tp.superTree();
        if ( ( tp.getPhylogeny() != whole ) || tp.isCurrentTreeIsSubtree() ) {
            fail( ok, "subtree: the whole tree must be back on display" );
            return;
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !pressed.equals( lengths( whole ) ) ) {
            fail( ok, "going back up changes nothing: the mode stays, and every length" );
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
        final Phylogeny view = tp.getPhylogeny();
        if ( !"Divergence".equals( cp.getBranchLengthsSelection() ) || !"subs/site".equals( view.getDistanceUnit() ) ) {
            fail( ok, "the view of a tree in divergence is in divergence, and says so" );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( length( whole, "isolate_C" ), 0.8 )
                || !near( length( whole, "isolate_A" ), 1.2 ) || !near( sum( whole ), TIME_SUM ) ) {
            fail( ok, "pressed in the view, the WHOLE tree is in time; C=" + length( whole, "isolate_C" ) + " A="
                    + length( whole, "isolate_A" ) + " sum=" + sum( whole ) );
        }
        if ( !"time".equals( whole.getDistanceUnit() ) || !"time".equals( view.getDistanceUnit() ) ) {
            fail( ok, "the whole tree and the view must both say time; got \"" + whole.getDistanceUnit() + "\" and \""
                    + view.getDistanceUnit() + "\"" );
        }
        tp.superTree();
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( sum( whole ), TIME_SUM )
                || !"Time".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "back on the whole tree: time, as it was left" );
        }
    }

    /** The switch is judged on the WHOLE tree. (A,B) states dates and a rate on both of its branches, but it is a
     *  clade of the tree whose isolate_C has no rate: its view is not offered the switch, asked in any way, and not
     *  a length of the tree is rewritten. */
    private static void subtreeOfATreeWithoutTheSwitch( final MainFrame frame,
                                                        final TreePanel tp,
                                                        final ControlPanel cp,
                                                        final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        final java.util.List<Double> arrived = lengths( whole );
        tp.subTree( cladeOf( whole, "isolate_A" ) );
        final Phylogeny view = tp.getPhylogeny();
        if ( !BranchLengthLayout.isApplicable( view ) ) {
            fail( ok, "clade: fixture -- taken by itself the view of (A,B) has both layouts" );
            return;
        }
        cp.populateBranchLengthsControl(); // what a change of tab does
        if ( tp.isBranchLengthToggleApplicable() || cp.isBranchLengthsControlVisible() ) {
            fail( ok, "a view of a tree that is not offered the switch is not offered it either" );
        }
        tp.setBranchLengthMode( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !arrived.equals( lengths( whole ) ) ) {
            fail( ok, "asked directly in the view, the panel must refuse, and rewrite nothing; A="
                    + length( whole, "isolate_A" ) );
        }
        tp.superTree();
        if ( !arrived.equals( lengths( whole ) ) || cp.isBranchLengthsControlVisible() ) {
            fail( ok, "back on the whole tree nothing has changed" );
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
        final String unit_before = whole.getDistanceUnit();
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
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !same( unit_before, whole.getDistanceUnit() ) ) {
            fail( ok, "navigation alone must change neither the mode nor the unit" );
        }
    }

    /** Two levels down, the switch pressed at the bottom: every tree on the way up -- the view of (D,E), the view of
     *  (C,(D,E)) it was descended from, and the whole tree -- is in the one layout, and says so. */
    private static void twoLevelsDown( final MainFrame frame,
                                       final TreePanel tp,
                                       final ControlPanel cp,
                                       final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        tp.subTree( cladeOf( whole, "isolate_C" ) );
        final Phylogeny middle = tp.getPhylogeny();
        tp.subTree( cladeOf( whole, "isolate_D" ) );
        final Phylogeny bottom = tp.getPhylogeny();
        if ( ( middle.getNodeCount() != 5 ) || ( bottom.getNodeCount() != 3 ) ) {
            fail( ok, "two levels: fixture -- views of 5 and of 3 nodes, got " + middle.getNodeCount() + " and "
                    + bottom.getNodeCount() );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( !near( sum( whole ), DIV_SUM ) || !near( length( whole, "isolate_A" ), 0.00372 ) ) {
            fail( ok, "pressed two levels down, the WHOLE tree is in divergence; sum=" + sum( whole ) );
        }
        for( final Phylogeny tree : new Phylogeny[] { whole, middle, bottom } ) {
            if ( !"subs/site".equals( tree.getDistanceUnit() ) ) {
                fail( ok, "every tree of the tab must say subs/site; the one of " + tree.getNodeCount() + " nodes says \""
                        + tree.getDistanceUnit() + "\"" );
            }
        }
        // the root of each view states the length of the node IT stands for: (C,(D,E)) is 1.3 x 0.0029, (D,E) is
        // 0.3 x 0.0034
        if ( !near( middle.getRoot().getDistanceToParent(), 0.00377 ) || !near( bottom.getRoot().getDistanceToParent(), 0.00102 ) ) {
            fail( ok, "the roots of the two views must state 0.00377 and 0.00102; got "
                    + middle.getRoot().getDistanceToParent() + " and " + bottom.getRoot().getDistanceToParent() );
        }
        tp.superTree();
        tp.superTree();
        if ( ( tp.getPhylogeny() != whole ) || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || !near( sum( whole ), DIV_SUM ) ) {
            fail( ok, "two levels up again: the whole tree, in divergence, as the press left it" );
        }
    }

    /** Divergence on screen, and a rate taken away in a view, no button pressed. Going back up rewrites nothing:
     *  divergence is still what every branch holds, so the mode still says so and the control stays -- it is the
     *  way back. Pressed, Time lays the whole tree out by time, and the control goes. */
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
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !near( sum( whole ), DIV_SUM ) ) {
            fail( ok, "going back up rewrites nothing: divergence, every branch; mode " + tp.getBranchLengthMode()
                    + ", sum=" + sum( whole ) );
        }
        tp.invalidateBranchLengthToggle();
        cp.populateBranchLengthsControl();
        if ( tp.isBranchLengthToggleApplicable() || !cp.isBranchLengthsControlVisible() ) {
            fail( ok, "the tree has lost its second layout, but divergence is on screen: the control stays" );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( sum( whole ), TIME_SUM ) ) {
            fail( ok, "the way back: time, every branch; sum=" + sum( whole ) );
        }
        if ( cp.isBranchLengthsControlVisible() ) {
            fail( ok, "...and the control gone" );
        }
    }

    /** A view is laid out the way ITS TREE's dates run. The view of (C,(D,E)) of the tree whose (D,E) node is dated
     *  before its parent: taken by itself, of its four spans one runs backwards -- but the sign of each is the one
     *  the whole tree gives it, whichever tree is on view when the switch is pressed. */
    private static void aViewRunsTheWayItsTreeDoes( final MainFrame frame,
                                                    final TreePanel tp,
                                                    final ControlPanel cp,
                                                    final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        // the node of the WHOLE tree, taken before the view descends: in the view, D's parent is the view's root, a
        // copy of this node
        final PhylogenyNode de = cladeOf( whole, "isolate_D" );
        tp.subTree( de ); // the view of (D,E): its own root stands for the node dated backwards
        final PhylogenyNode view_root = tp.getPhylogeny().getRoot();
        if ( ( view_root == de ) || ( whole.getNode( "isolate_D" ).getParent() != view_root ) ) {
            fail( ok, "backwards in a view: fixture -- the view's root must be a copy, and D's parent while on view" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || ( Double.doubleToRawLongBits( de.getDistanceToParent() ) != 0L )
                || !near( length( whole, "isolate_D" ), 0.85 * 0.0035 ) ) {
            fail( ok, "pressed in the view of (D,E): the span that runs backwards is 0 in divergence, D is 0.85 x"
                    + " 0.0035; got (D,E)=" + de.getDistanceToParent() + " D=" + length( whole, "isolate_D" ) );
        }
        if ( Double.doubleToRawLongBits( view_root.getDistanceToParent() ) != 0L ) {
            fail( ok, "the view's root states the length of the node it stands for, got "
                    + view_root.getDistanceToParent() );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( !near( de.getDistanceToParent(), -0.05 ) || !near( length( whole, "isolate_D" ), 0.85 )
                || !near( length( whole, "isolate_A" ), 1.2 ) ) {
            fail( ok, "and back in time it is -0.05, D 0.85, A 1.2; got (D,E)=" + de.getDistanceToParent() + " D="
                    + length( whole, "isolate_D" ) + " A=" + length( whole, "isolate_A" ) );
        }
        if ( !near( view_root.getDistanceToParent(), -0.05 ) ) {
            fail( ok, "...and the view's root with it, got " + view_root.getDistanceToParent() );
        }
    }

    /** A cherry of a tree on CALENDAR time, one of its two tips stated 0.25 BEFORE their parent. Pressed in the view
     *  of the cherry: that tip is 0 in divergence and the other is not, and back in time both are what they were. */
    private static void aCherryOfACalendarTree( final MainFrame frame,
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
            fail( ok, "cherry: " + CALENDAR + " must have a cherry" );
            return;
        }
        final PhylogenyNode before = cherry.getChildNode( 0 );
        final PhylogenyNode after = cherry.getChildNode( 1 );
        final double after_length = after.getDistanceToParent();
        before.setDistanceToParent( -0.25 );
        if ( !( after_length > 0 ) || !tp.isBranchLengthToggleApplicable() ) {
            fail( ok, "cherry: fixture -- the other tip must have a length above 0, and the tree be offered the switch" );
            return;
        }
        tp.subTree( cherry );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( Double.doubleToRawLongBits( before.getDistanceToParent() ) != 0L ) || !( after.getDistanceToParent() > 0 )
                || !( after.getDistanceToParent() < ( after_length / 10 ) ) ) {
            fail( ok, "in divergence the tip stated before its parent is at 0 and the other is rate x its length; got "
                    + before.getDistanceToParent() + " and " + after.getDistanceToParent() );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( before.getDistanceToParent() != -0.25 ) || ( after.getDistanceToParent() != after_length ) ) {
            fail( ok, "back in time both tips have the lengths they had; got " + before.getDistanceToParent() + " and "
                    + after.getDistanceToParent() + ", were -0.25 and " + after_length );
        }
    }

    /** An undo INSIDE a view puts a copy on display, and brings back the layout it was captured in. The trees the
     *  view was descended from follow it: back on the whole tree, one layout, the one the panel says. */
    private static void undoInsideASubtree( final MainFrame frame,
                                            final TreePanel tp,
                                            final ControlPanel cp,
                                            final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        tp.subTree( cladeOf( whole, "isolate_C" ) );
        tp.pushUndoCheckpoint( "Rename" ); // captured in time
        tp.getPhylogeny().getNode( "isolate_E" ).setName( "isolate_E_renamed" );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !near( sum( whole ), DIV_SUM ) ) {
            fail( ok, "undo in a view: the whole tree must BE in divergence first" );
            return;
        }
        if ( !tp.undo() ) {
            fail( ok, "undo in a view: there must be something to undo" );
            return;
        }
        final Phylogeny restored = tp.getPhylogeny();
        if ( ( restored.getNodeCount() != 5 ) || !near( restored.getNode( "isolate_C" ).getDistanceToParent(), 0.8 ) ) {
            fail( ok, "undo in a view: fixture -- the restored view must be the one captured, in time lengths" );
            return;
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !"Time".equals( cp.getBranchLengthsSelection() ) ) {
            fail( ok, "undo brought time lengths back: the mode must be TIME" );
        }
        if ( !near( sum( whole ), TIME_SUM ) || !near( length( whole, "isolate_A" ), 1.2 )
                || !"time".equals( whole.getDistanceUnit() ) ) {
            fail( ok, "...and the tree the view was descended from must follow: time, every branch; sum="
                    + sum( whole ) + ", unit \"" + whole.getDistanceUnit() + "\"" );
        }
        // the switch pressed now reaches the copy on display AND the whole tree under it
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( !near( restored.getNode( "isolate_C" ).getDistanceToParent(), 0.00208 ) || !near( sum( whole ), DIV_SUM ) ) {
            fail( ok, "pressed after the undo, the copy on display and the whole tree are both in divergence; C="
                    + restored.getNode( "isolate_C" ).getDistanceToParent() + " sum=" + sum( whole ) );
        }
    }

    /** ...and the other way round, inside a view: an undo that brings DIVERGENCE back while the tab shows time. The
     *  trees beneath the view are showing time, and what they show is what is kept -- isolate_A, outside the view,
     *  was given a length of 1.17 while time was on screen, AFTER the tab had last left time. */
    /**
     * The copy an undo puts on display inside a view, of a clade with NO change in it (D and E at a rate of 0: the
     * whole tree has depth, the clade by itself none). Div lays the whole tree out, and the copy with it: a part of
     * a tree in divergence is in divergence, depth or none -- not left in years under the unit of divergence.
     */
    private static void aCopyOfACladeWithoutChange( final MainFrame frame,
                                                    final TreePanel tp,
                                                    final ControlPanel cp,
                                                    final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        for( final String tip : new String[] { "isolate_D", "isolate_E" } ) {
            final PhylogenyNode n = whole.getNode( tip );
            n.getNodeData().getProperties().getProperties().removeIf( p -> RATE.equals( p.getRef() ) );
            n.getNodeData().getProperties().addProperty( new org.forester.phylogeny.data.Property( RATE, "0", "",
                    "xsd:decimal", org.forester.phylogeny.data.Property.AppliesTo.NODE ) );
        }
        tp.invalidateBranchLengthToggle();
        tp.subTree( cladeOf( whole, "isolate_D" ) );
        tp.pushUndoCheckpoint( "Rename" );
        tp.getPhylogeny().getNode( "isolate_D" ).setName( "isolate_D_renamed" );
        if ( !tp.undo() ) {
            fail( ok, "copy without change: there must be something to undo" );
            return;
        }
        final Phylogeny copy = tp.getPhylogeny();
        // the rename reached the whole tree (a view shares its nodes); the copy is from before it
        if ( !hasNode( whole, "isolate_D_renamed" ) || !hasNode( copy, "isolate_D" ) || !near( length( copy, "isolate_D" ), 0.5 )
                || BranchLengthLayout.isApplicable( copy ) || !tp.isBranchLengthToggleApplicable() ) {
            fail( ok, "copy without change: fixture -- a COPY on display, in time, with no depth of its own, of a"
                    + " tree that is offered the switch" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !near( sum( whole ), DIV_SUM - ( 0.5 * 0.0035 ) - ( 0.5 * 0.0033 ) ) ) {
            fail( ok, "copy without change: fixture -- the whole tree must be in divergence; sum=" + sum( whole ) );
            return;
        }
        if ( ( length( copy, "isolate_D" ) != 0.0 ) || ( length( copy, "isolate_E" ) != 0.0 )
                || !"subs/site".equals( copy.getDistanceUnit() ) ) {
            fail( ok, "a copy of a clade without change is laid out in divergence with its tree: D and E 0 subs/site;"
                    + " got D=" + length( copy, "isolate_D" ) + " E=" + length( copy, "isolate_E" ) + " \""
                    + copy.getDistanceUnit() + "\"" );
        }
    }

    private static void undoIntoDivergenceInsideASubtree( final MainFrame frame,
                                                          final TreePanel tp,
                                                          final ControlPanel cp,
                                                          final boolean[] ok ) {
        final Phylogeny whole = tp.getPhylogeny();
        tp.subTree( cladeOf( whole, "isolate_C" ) );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        tp.pushUndoCheckpoint( "Rename" ); // captured in divergence
        tp.getPhylogeny().getNode( "isolate_E" ).setName( "isolate_E_renamed" );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        whole.getNode( "isolate_A" ).setDistanceToParent( 1.17 );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( sum( whole ), TIME_SUM - 0.03 ) ) {
            fail( ok, "undo to divergence in a view: the whole tree must BE in time first, A at 1.17" );
            return;
        }
        if ( !tp.undo() ) {
            fail( ok, "undo to divergence in a view: there must be something to undo" );
            return;
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || !near( tp.getPhylogeny().getNode( "isolate_C" ).getDistanceToParent(), 0.00208 ) ) {
            fail( ok, "undo to divergence in a view: fixture -- the restored view must be in divergence" );
            return;
        }
        if ( !near( length( whole, "isolate_A" ), 1.17 * 0.0031 ) || !near( length( whole, "isolate_B" ), 0.00336 ) ) {
            fail( ok, "the tree beneath the view follows into divergence, from the lengths it was SHOWING; A="
                    + length( whole, "isolate_A" ) + " B=" + length( whole, "isolate_B" ) );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( length( whole, "isolate_A" ) != 1.17 ) || ( length( whole, "isolate_B" ) != 1.2 ) ) {
            fail( ok, "...and Time gives those back: A 1.17, B 1.2; got " + length( whole, "isolate_A" ) + " and "
                    + length( whole, "isolate_B" ) );
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

    /**
     * An undo INTO a snapshot taken in divergence takes the lengths in time OF THAT MOMENT, not the ones the tab
     * kept since: B is lengthened in time to 1.3 after the snapshot (the edit is taken back by the undo), and Div is
     * pressed over it, which keeps 1.3. Twice: once with the undo crossing from time into divergence, once from
     * divergence into divergence, where the mode does not change at all.
     */
    private static void undoIntoDivergenceTakesItsOwnTime( final MainFrame frame,
                                                           final TreePanel tp,
                                                           final ControlPanel cp,
                                                           final boolean[] ok ) {
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        tp.pushUndoCheckpoint( "Rename" ); // S1: in divergence, B 1.2 in time
        tp.getPhylogeny().getNode( "isolate_A" ).setName( "isolate_A_renamed" );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        tp.pushUndoCheckpoint( "Length" ); // S2: in time, B 1.2
        tp.getPhylogeny().getNode( "isolate_B" ).setDistanceToParent( 1.3 );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE ); // keeps B at 1.3
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( length( tp.getPhylogeny(), "isolate_B" ) != 1.3 ) || !tp.undo() || !tp.undo() ) {
            fail( ok, "own time: fixture -- B 1.3 in time, and two undos to take" );
            return;
        }
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || !hasNode( tp.getPhylogeny(), "isolate_A" ) ) {
            fail( ok, "own time: fixture -- the second undo must restore S1, in divergence" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || ( length( tp.getPhylogeny(), "isolate_B" ) != 1.2 ) ) {
            fail( ok, "an undo into divergence takes the lengths in time the snapshot was taken with: B 1.2, not the"
                    + " 1.3 of an edit the undo took back; got " + length( tp.getPhylogeny(), "isolate_B" ) );
            return;
        }
        // divergence into divergence: the mode does not change, the lengths in time still do
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        tp.pushUndoCheckpoint( "Rename" ); // S3: in divergence, B 1.2 in time
        tp.getPhylogeny().getNode( "isolate_C" ).setName( "isolate_C_renamed" );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        tp.getPhylogeny().getNode( "isolate_B" ).setDistanceToParent( 1.3 );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE ); // keeps B at 1.3
        if ( !tp.undo() || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || !hasNode( tp.getPhylogeny(), "isolate_C" ) ) {
            fail( ok, "own time: fixture -- the undo must restore S3, in divergence" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( length( tp.getPhylogeny(), "isolate_B" ) != 1.2 ) {
            fail( ok, "an undo from divergence into divergence also takes the snapshot's lengths in time: B 1.2; got "
                    + length( tp.getPhylogeny(), "isolate_B" ) );
        }
    }

    /**
     * A REDO into divergence takes the lengths in time the tab had when the undo filed it: B was lengthened to 1.3
     * in time as part of the edit, then Div pressed; the undo goes back to before the edit (B 1.2 in time), and
     * the redo must bring back 1.3 -- a redo that filed no lengths would take B's time from its dates, 1.2.
     */
    private static void redoIntoDivergenceTakesItsOwnTime( final MainFrame frame,
                                                           final TreePanel tp,
                                                           final ControlPanel cp,
                                                           final boolean[] ok ) {
        tp.pushUndoCheckpoint( "Edit" ); // in time, B 1.2
        tp.getPhylogeny().getNode( "isolate_A" ).setName( "isolate_A_renamed" );
        tp.getPhylogeny().getNode( "isolate_B" ).setDistanceToParent( 1.3 );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE ); // keeps B at 1.3
        if ( !tp.undo() || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME )
                || ( length( tp.getPhylogeny(), "isolate_B" ) != 1.2 ) ) {
            fail( ok, "redo own time: fixture -- the undo must go back to time, B 1.2" );
            return;
        }
        if ( !tp.redo() || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || !hasNode( tp.getPhylogeny(), "isolate_A_renamed" ) ) {
            fail( ok, "redo own time: fixture -- the redo must bring back the edited tree, in divergence" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( length( tp.getPhylogeny(), "isolate_B" ) != 1.3 ) {
            fail( ok, "a redo into divergence takes the lengths in time it was filed with: B 1.3; got "
                    + length( tp.getPhylogeny(), "isolate_B" ) );
        }
    }

    /**
     * A redo FILES the tree it replaces for undo, with that tree's lengths in time: B is 1.25 in time when Div is
     * pressed over the undone state, the redo leaves it, and undoing the redo must give 1.25 back -- filed with no
     * lengths, its time would come from its dates, 1.2.
     */
    private static void redoFilesTheTimeOfTheTreeItLeaves( final MainFrame frame,
                                                           final TreePanel tp,
                                                           final ControlPanel cp,
                                                           final boolean[] ok ) {
        tp.pushUndoCheckpoint( "Rename" );
        tp.getPhylogeny().getNode( "isolate_A" ).setName( "isolate_A_renamed" );
        if ( !tp.undo() ) {
            fail( ok, "redo files time: fixture -- there must be something to undo" );
            return;
        }
        tp.getPhylogeny().getNode( "isolate_B" ).setDistanceToParent( 1.25 );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE ); // keeps B at 1.25
        if ( !tp.canRedo() || !tp.redo() || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) ) {
            fail( ok, "redo files time: fixture -- the redo must bring back the renamed tree, in time" );
            return;
        }
        if ( !tp.undo() || ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || !hasNode( tp.getPhylogeny(), "isolate_A" ) ) {
            fail( ok, "redo files time: fixture -- undoing the redo must restore the tree it left, in divergence" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( length( tp.getPhylogeny(), "isolate_B" ) != 1.25 ) {
            fail( ok, "the tree a redo left comes back with its lengths in time: B 1.25; got "
                    + length( tp.getPhylogeny(), "isolate_B" ) );
        }
    }

    /**
     * Delete Node while divergence is on screen adds the removed node's length to its child's, in divergence; Time
     * must do the same with the lengths KEPT, or the child spans only its own time. D is given 0.55 (its date gap
     * 0.5), so what it spans after the delete, 0.55 + 0.3 = 0.85, is neither its own kept length nor the gap to its
     * new parent (0.8).
     */
    private static void deleteNodeInDivergence( final MainFrame frame,
                                                final TreePanel tp,
                                                final ControlPanel cp,
                                                final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        phy.getNode( "isolate_D" ).setDistanceToParent( 0.55 );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) {
            fail( ok, "delete node: fixture -- the tree must be in divergence" );
            return;
        }
        tp.deleteNodeOrSubtreeConfirmed( phy.getNode( "isolate_D" ).getParent(), true );
        if ( phy.getNode( "isolate_D" ).getParent() != phy.getNode( "isolate_C" ).getParent() ) {
            fail( ok, "delete node: fixture -- D and E must now hang from C's parent" );
            return;
        }
        final double on_screen = length( phy, "isolate_D" );
        final double pieces = ( 0.0035 * 0.55 ) + ( 0.0034 * 0.3 );
        if ( !near( on_screen, pieces ) ) {
            fail( ok, "delete node: fixture -- in divergence D is its piece and the removed one's, " + pieces + "; got "
                    + on_screen );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( length( phy, "isolate_D" ), 0.85 )
                || !near( length( phy, "isolate_E" ), 0.8 ) ) {
            fail( ok, "a node deleted in divergence: its children span its time and their own, D 0.85 and E 0.8; got D="
                    + length( phy, "isolate_D" ) + " E=" + length( phy, "isolate_E" ) );
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        if ( !near( length( phy, "isolate_D" ), pieces ) ) {
            fail( ok, "Div, Time, Div after a delete gives back the same divergence, piece by piece, each at its own"
                    + " rate: " + pieces + "; got " + length( phy, "isolate_D" ) );
        }
    }

    /**
     * Delete Node while TIME is on screen: the merged branch is 0.8 in time, and its divergence is still its two
     * pieces at their own rates, 0.0035 x 0.5 + 0.0034 x 0.3 -- not D's rate over the whole, 0.0035 x 0.8. Both
     * ways round, twice.
     */
    private static void deleteNodeInTime( final MainFrame frame,
                                          final TreePanel tp,
                                          final ControlPanel cp,
                                          final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        tp.deleteNodeOrSubtreeConfirmed( phy.getNode( "isolate_D" ).getParent(), true );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near( length( phy, "isolate_D" ), 0.8 ) ) {
            fail( ok, "delete in time: fixture -- in time, D now 0.8 under C's parent" );
            return;
        }
        final double pieces = ( 0.0035 * 0.5 ) + ( 0.0034 * 0.3 );
        for( int round = 1; round <= 2; ++round ) {
            cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
            if ( !near( length( phy, "isolate_D" ), pieces ) ) {
                fail( ok, "a node deleted in time: Div (round " + round + ") gives D its pieces at their own rates, "
                        + pieces + ", not " + ( 0.0035 * 0.8 ) + "; got " + length( phy, "isolate_D" ) );
                return;
            }
            cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
            if ( !near( length( phy, "isolate_D" ), 0.8 ) ) {
                fail( ok, "a node deleted in time: Time (round " + round + ") gives D 0.8; got " + length( phy, "isolate_D" ) );
                return;
            }
        }
    }

    /**
     * A tab made from a tab in TIME, its copy pruned (B goes, and A's parent with it): in Div, A is its own piece at
     * its rate and the removed parent's at the parent's, 0.0031 x 1.407 + 0.0030 x 0.9.
     */
    private static void aTabMadeFromATimeTab( final MainFrame frame,
                                              final TreePanel tp,
                                              final ControlPanel cp,
                                              final boolean[] ok ) {
        final Phylogeny copy = tp.getPhylogeny().copy();
        copy.deleteSubtree( copy.getNode( "isolate_B" ), true );
        copy.externalNodesHaveChanged();
        copy.clearHashIdToNodeMap();
        if ( ( copy.getNode( "isolate_A" ).getParent() != copy.getRoot() ) || !near( length( copy, "isolate_A" ), 2.307 ) ) {
            fail( ok, "time tab: fixture -- A hangs from the root of the pruned copy at 1.407 + 0.9" );
            return;
        }
        ( ( MainFrameApplication ) frame ).addDerivedPhylogenyInNewTab( copy );
        final TreePanel derived = frame.getMainPanel().getCurrentTreePanel();
        if ( ( derived == tp ) || ( derived.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) ) {
            fail( ok, "a tab made from a tab in time is in time" );
            return;
        }
        frame.getMainPanel().getControlPanel().userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        final double pieces = ( 0.0031 * 1.407 ) + ( 0.0030 * 0.9 );
        if ( !near( length( copy, "isolate_A" ), pieces ) ) {
            fail( ok, "the new tab's Div: A is its pieces at their own rates, " + pieces + "; got " + length( copy, "isolate_A" ) );
        }
    }

    /**
     * Delete Subtree while divergence is on screen: C goes, its parent is left with one child and goes too, and the
     * (D,E) node now hangs from the root, spanning its own kept 0.35 (its date gap 0.3) and its old parent's 1.3.
     */
    private static void deleteSubtreeInDivergence( final MainFrame frame,
                                                   final TreePanel tp,
                                                   final ControlPanel cp,
                                                   final boolean[] ok ) {
        final Phylogeny phy = tp.getPhylogeny();
        final PhylogenyNode de = phy.getNode( "isolate_D" ).getParent();
        de.setDistanceToParent( 0.35 );
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        tp.deleteNodeOrSubtreeConfirmed( phy.getNode( "isolate_C" ), false );
        if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || hasNode( phy, "isolate_C" )
                || ( phy.getNode( "isolate_D" ).getParent().getParent() != phy.getRoot() ) ) {
            fail( ok, "delete subtree: fixture -- in divergence, C gone, (D,E) now a child of the root" );
            return;
        }
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( !near( phy.getNode( "isolate_D" ).getParent().getDistanceToParent(), 1.65 )
                || !near( length( phy, "isolate_D" ), 0.5 ) ) {
            fail( ok, "a subtree deleted in divergence: the node left joins its parent's time to its own, 0.35 + 1.3;"
                    + " got " + phy.getNode( "isolate_D" ).getParent().getDistanceToParent() );
        }
    }

    /**
     * A tab made from a tab in divergence (Select Representative Tips) starts in divergence, with the source tab's
     * lengths in time: A states 1.407 (its dates 1.2 apart). The copy is pruned (B goes, and A's parent with it), so
     * A spans 1.407 + 0.9 -- a tab that took its time from its dates would give it 2.1.
     */
    private static void aTabMadeFromADivergenceTab( final MainFrame frame,
                                                    final TreePanel tp,
                                                    final ControlPanel cp,
                                                    final boolean[] ok ) {
        cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
        final Phylogeny copy = tp.getPhylogeny().copy();
        copy.deleteSubtree( copy.getNode( "isolate_B" ), true );
        copy.externalNodesHaveChanged();
        copy.clearHashIdToNodeMap();
        if ( ( copy.getNode( "isolate_A" ).getParent() != copy.getRoot() ) || !( frame instanceof MainFrameApplication ) ) {
            fail( ok, "derived tab: fixture -- A must hang from the root of the pruned copy" );
            return;
        }
        ( ( MainFrameApplication ) frame ).addDerivedPhylogenyInNewTab( copy );
        final TreePanel derived = frame.getMainPanel().getCurrentTreePanel();
        if ( ( derived == tp ) || ( derived.getPhylogeny() != copy ) ) {
            fail( ok, "derived tab: fixture -- the copy must be in a new, current tab" );
            return;
        }
        if ( ( derived.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE )
                || ( derived.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) ) {
            fail( ok, "a tab made from a tab in divergence is in divergence, with no time axis" );
        }
        frame.getMainPanel().getControlPanel().userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
        if ( ( derived.getBranchLengthMode() != BranchLengthLayout.MODE.TIME )
                || !near( length( copy, "isolate_A" ), 2.307 ) ) {
            fail( ok, "the new tab's Time gives back the source tab's lengths in time, A 1.407 + 0.9; got "
                    + length( copy, "isolate_A" ) );
        }
    }

    /**
     * A tree saved from the Div view, opened in a tab that is NOT the current one (two trees: the last is selected):
     * its time axis read FIRST is none -- the arrival is settled before the axis is answered, not whenever some
     * other accessor happens to ask. On calendar time, where the axis it would otherwise derive is CALENDAR.
     */
    private static boolean aBackgroundTabArrivesBeforeItsAxisIsRead() throws Exception {
        final File saved = File.createTempFile( "saved_in_div", ".xml" );
        saved.deleteOnExit();
        final boolean[] written = { false };
        if ( !inFrame( RECORDED, ( f, tp, cp, ok ) -> {
            cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
            try {
                new org.forester.io.writers.PhylogenyWriter().toPhyloXML( saved, tp.getPhylogeny(), 0 );
                written[ 0 ] = tp.getBranchLengthMode() == BranchLengthLayout.MODE.DIVERGENCE;
            }
            catch ( final Exception e ) {
                fail( ok, "background tab: could not be saved: " + e );
            }
        } ) || !written[ 0 ] ) {
            return fail( "background tab: " + RECORDED + " must be saved from the Div view" );
        }
        final Phylogeny reopened = FigureRenderer.readTrees( saved )[ 0 ];
        if ( AptxUtil.deriveTimeAxisType( reopened ) != TIME_AXIS_TYPE.CALENDAR ) {
            return fail( "background tab: fixture -- the tree by itself derives a CALENDAR axis" );
        }
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { reopened, read( RATED ) }, new Configuration(), "saved" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            try {
                final TreePanel background = mf[ 0 ].getMainPanel().getTreePanels().get( 0 );
                if ( ( background.getPhylogeny() != reopened ) || ( mf[ 0 ].getMainPanel().getCurrentTreePanel() == background ) ) {
                    fail( ok, "background tab: fixture -- the saved tree must be in the tab that is not current" );
                    return;
                }
                if ( background.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) {
                    fail( ok, "a tree that arrives in divergence has no time axis, also when the axis is the first"
                            + " thing asked; got " + background.effectiveTimeAxisType() );
                }
                if ( background.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) {
                    fail( ok, "background tab: the tree arrives in divergence" );
                }
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

    /**
     * A tree SAVED while divergence is on screen, and opened again. It carries its dates, its rates or recorded
     * divergence, and lengths that are its DIVERGENCE: the tab must say so -- Div pressed, no time axis -- and Time
     * must lay it out by time, from its dates. Taken for a time tree it said Time over a picture of divergence, gave
     * the divergence back as its time, and multiplied its rates by it.
     */
    private static boolean arrivesShowingDivergence( final String demo, final double time_sum, final double div_sum )
            throws Exception {
        final File saved = File.createTempFile( "saved_in_div", ".xml" );
        saved.deleteOnExit();
        final boolean[] written = { false };
        if ( !inFrame( demo, ( f, tp, cp, ok ) -> {
            if ( !near9( sum( tp.getPhylogeny() ), time_sum ) ) {
                fail( ok, demo + ": fixture -- it must open in time, its lengths adding up to " + time_sum + "; got "
                        + sum( tp.getPhylogeny() ) );
                return;
            }
            cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
            if ( !near9( sum( tp.getPhylogeny() ), div_sum ) ) {
                fail( ok, demo + ": fixture -- in divergence its lengths must add up to " + div_sum + "; got "
                        + sum( tp.getPhylogeny() ) );
                return;
            }
            try {
                new org.forester.io.writers.PhylogenyWriter().toPhyloXML( saved, tp.getPhylogeny(), 0 );
                written[ 0 ] = true;
            }
            catch ( final Exception e ) {
                fail( ok, demo + ": could not be saved: " + e );
            }
        } ) || !written[ 0 ] ) {
            return false;
        }
        final Phylogeny reopened = FigureRenderer.readTrees( saved )[ 0 ];
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { reopened }, new Configuration(), "saved" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            try {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
                final Phylogeny phy = tp.getPhylogeny();
                if ( !near9( sum( phy ), div_sum ) ) {
                    fail( ok, demo + " reopened: fixture -- the file must hold the divergence lengths, got " + sum( phy ) );
                    return;
                }
                if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) {
                    fail( ok, demo + " saved from the Div view arrives SHOWING DIVERGENCE: the mode must say so" );
                }
                if ( !cp.isBranchLengthsControlVisible() || !"Divergence".equals( cp.getBranchLengthsSelection() ) ) {
                    fail( ok, demo + " reopened: the control must be there, Div pressed" );
                }
                if ( tp.isBranchLengthTimeCalibrated() || ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) ) {
                    fail( ok, demo + " reopened: divergence has no time axis and no node-age bars" );
                }
                if ( !tp.isBranchLengthToggleApplicable() ) {
                    fail( ok, demo + " reopened: it has both layouts, its time from its dates" );
                }
                cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
                if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) || !near9( sum( phy ), time_sum ) ) {
                    fail( ok, demo + " reopened: Time must lay it out by time, " + time_sum + "; got " + sum( phy ) );
                }
                cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
                if ( !near9( sum( phy ), div_sum ) ) {
                    fail( ok, demo + " reopened: and Div must bring its divergence back, " + div_sum + "; got " + sum( phy ) );
                }
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

    /**
     * A tree that arrived showing divergence has kept no length: its time is read off its dates, signed the way its
     * dates run -- the WHOLE tree's, also while a view of it is up. The view here is a cherry on calendar time, one
     * tip dated BEFORE their parent: by itself a tie, which reads as ages and would give each of the two the other
     * sign.
     */
    private static boolean aViewOfATreeThatKeptNothing() throws Exception {
        final File saved = File.createTempFile( "saved_in_div", ".xml" );
        saved.deleteOnExit();
        final boolean[] written = { false };
        if ( !inFrame( RECORDED, ( f, tp, cp, ok ) -> {
            cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.DIVERGENCE );
            try {
                new org.forester.io.writers.PhylogenyWriter().toPhyloXML( saved, tp.getPhylogeny(), 0 );
                written[ 0 ] = tp.getBranchLengthMode() == BranchLengthLayout.MODE.DIVERGENCE;
            }
            catch ( final Exception e ) {
                fail( ok, "kept nothing: could not be saved: " + e );
            }
        } ) || !written[ 0 ] ) {
            return fail( "kept nothing: " + RECORDED + " must be saved from the Div view" );
        }
        final Phylogeny reopened = FigureRenderer.readTrees( saved )[ 0 ];
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { reopened }, new Configuration(), "saved" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            try {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
                final Phylogeny whole = tp.getPhylogeny();
                final PhylogenyNode before = whole.getNode( "A/Abidjan/3/2015" );
                final PhylogenyNode after = whole.getNode( "A/Dakar/8/2016" );
                final PhylogenyNode cherry = before.getParent();
                if ( ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.DIVERGENCE ) || ( after.getParent() != cherry )
                        || ( cherry.getNumberOfDescendants() != 2 ) ) {
                    fail( ok, "kept nothing: fixture -- the tree must arrive showing divergence, the two tips a cherry" );
                    return;
                }
                final double parent_date = cherry.getNodeData().getDate().getValue().doubleValue();
                final double after_gap = after.getNodeData().getDate().getValue().doubleValue() - parent_date;
                before.getNodeData().getDate().setValue( new java.math.BigDecimal( String.valueOf( parent_date - 0.25 ) ) );
                tp.subTree( cherry );
                if ( !BranchLengthLayout.datesIncreaseTowardTips( whole )
                        || BranchLengthLayout.datesIncreaseTowardTips( tp.getPhylogeny() ) || !( after_gap > 0 ) ) {
                    fail( ok, "kept nothing: fixture -- the whole tree on calendar time, the view by itself a tie" );
                    return;
                }
                cp.userSelectBranchLengthsForTest( BranchLengthLayout.MODE.TIME );
                if ( tp.getBranchLengthMode() != BranchLengthLayout.MODE.TIME ) {
                    fail( ok, "kept nothing: Time must take in the view" );
                }
                if ( ( Math.abs( before.getDistanceToParent() - ( -0.25 ) ) > 1e-9 )
                        || ( Math.abs( after.getDistanceToParent() - after_gap ) > 1e-9 ) ) {
                    fail( ok, "the dates are read the way the WHOLE tree's run: the tip dated before its parent -0.25, the"
                            + " other +" + after_gap + "; got " + before.getDistanceToParent() + " and "
                            + after.getDistanceToParent() );
                }
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

    private static boolean near9( final double a, final double b ) {
        return Math.abs( a - b ) < 1e-9;
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
        if ( ( parentLength( phy, "isolate_D" ) != -0.05 ) || ( length( phy, "isolate_D" ) != 0.85 ) ) {
            fail( ok, "back in time the length is -0.05 again, as the file states it; got (D,E)="
                    + parentLength( phy, "isolate_D" ) + " D=" + length( phy, "isolate_D" ) );
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
