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
import java.awt.Point;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Verifies the "next / previous search hit" step-through (#5.3) on the search-emphasis demo:
 * <ul>
 *   <li>the ordered hit list and "k / N" count match the found set;</li>
 *   <li>stepping forward/backward advances and wraps in both directions; the first forward step from the
 *       unpositioned state lands on the first hit (last, stepping backward);</li>
 *   <li>a hit hidden under a collapsed clade targets that clade's (drawn) triangle rather than the invisible node;</li>
 *   <li>re-assigning an <em>equal</em> found set preserves the position (the regression that made ⌘G re-run by the
 *       search box's keyReleased snap back to hit 0), while a <em>changed</em> set restarts it;</li>
 *   <li>{@code centerOnNode} actually scrolls the viewport to the target;</li>
 *   <li>manual node selection reveals/updates the navigator, and both search boxes drive it;</li>
 *   <li>the ControlPanel navigator row appears/disappears (a laid-out visibility flip, not just a flag) and the
 *       View-menu Find Next / Find Previous (⌘G / ⌘⇧G) dispatch steps.</li>
 * </ul>
 * Headful; a green no-op when headless.
 */
public final class SearchStepThroughTest {

    private static final int W = 640, H = 560;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SearchStepThrough: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return stepOk();
    }

    private static boolean stepOk() {
        final File file = new File( System.getProperty( "user.dir" ), "forester/demo/search-emphasis.xml" );
        if ( !file.exists() ) {
            return fail( "demo tree missing: " + file.getAbsolutePath() );
        }
        final MainFrame[] mf = new MainFrame[ 1 ];
        final boolean[] ok = { true };
        try {
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "step" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                final ControlPanel cp = tp.getControlPanel();
                cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                frame.showWhole();
                tp.setSize( W, H );
                tp.calcParametersForPainting( W, H );

                // the *_kinase external nodes, in the same preorder the step-through walks -- the expected hit list
                final List<PhylogenyNode> expected = new ArrayList<>();
                final Set<Long> found = new HashSet<>();
                for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if ( n.isExternal() && n.getName().contains( "kinase" ) ) {
                        expected.add( n );
                        found.add( n.getId() );
                    }
                }
                if ( expected.size() < 3 ) {
                    fail( ok, "demo must contain at least 3 *_kinase tips, got " + expected.size() );
                    return;
                }
                tp.setFoundNodes0( found );

                // (1) count matches the found set, and nothing is positioned until the first step
                if ( tp.getSearchHitCount() != expected.size() ) {
                    fail( ok, "hit count " + tp.getSearchHitCount() + " != found " + expected.size() );
                }
                if ( tp.getSearchHitIndex() != -1 ) {
                    fail( ok, "should start unpositioned, index=" + tp.getSearchHitIndex() );
                }

                // (2) forward walk: index 0..N-1, target = the matching hit in preorder
                for( int i = 0; i < expected.size(); i++ ) {
                    tp.stepToFoundNode( 1 );
                    if ( tp.getSearchHitIndex() != i ) {
                        fail( ok, "forward step " + i + " -> index " + tp.getSearchHitIndex() );
                    }
                    if ( tp.getLastStepTargetForTest() != expected.get( i ) ) {
                        fail( ok, "forward step " + i + " targeted the wrong node" );
                    }
                }
                // (3) one more forward step wraps back to the first hit
                tp.stepToFoundNode( 1 );
                if ( tp.getSearchHitIndex() != 0 ) {
                    fail( ok, "forward wrap should return to index 0, got " + tp.getSearchHitIndex() );
                }

                // (4) B1 regression: re-assigning an EQUAL found set (what the search box's keyReleased does when the
                //     user presses ⌘G without editing the query) must NOT move the position; a CHANGED set restarts.
                resetPosition( tp, found );
                tp.stepToFoundNode( 1 );
                tp.stepToFoundNode( 1 ); // index 1
                tp.setFoundNodes0( new HashSet<>( found ) ); // equal contents, different object
                if ( tp.getSearchHitIndex() != 1 ) {
                    fail( ok, "re-setting an equal found set must preserve the position, got "
                            + tp.getSearchHitIndex() );
                }
                final Set<Long> subset = new HashSet<>();
                subset.add( found.iterator().next() );
                tp.setFoundNodes0( subset ); // different contents -> restart
                if ( tp.getSearchHitIndex() != -1 ) {
                    fail( ok, "a changed found set must restart the step-through, got " + tp.getSearchHitIndex() );
                }

                // (5) backward from the unpositioned state lands on the LAST hit, then wraps below 0
                resetPosition( tp, found );
                tp.stepToFoundNode( -1 );
                if ( tp.getSearchHitIndex() != ( expected.size() - 1 ) ) {
                    fail( ok, "first backward step should land on the last hit, got " + tp.getSearchHitIndex() );
                }
                tp.stepToFoundNode( 1 ); // -> wraps to 0
                tp.stepToFoundNode( -1 ); // 0 - 1 -> wraps to last
                if ( tp.getSearchHitIndex() != ( expected.size() - 1 ) ) {
                    fail( ok, "backward wrap below 0 should reach the last hit, got " + tp.getSearchHitIndex() );
                }

                // (6) a hit hidden under a collapsed clade steps to that clade's triangle, not the invisible node
                final PhylogenyNode first_hit = expected.get( 0 );
                final PhylogenyNode collapsed = first_hit.getParent();
                collapsed.setCollapse( true );
                resetPosition( tp, found );
                tp.stepToFoundNode( 1 ); // -> index 0 == first_hit, which is now hidden
                if ( tp.getSearchHitIndex() != 0 ) {
                    fail( ok, "collapsed-clade step should still be index 0, got " + tp.getSearchHitIndex() );
                }
                if ( tp.getLastStepTargetForTest() != collapsed ) {
                    fail( ok, "a hidden hit must target its collapsed ancestor (the drawn triangle)" );
                }
                if ( tp.getLastStepTargetForTest() == first_hit ) {
                    fail( ok, "a hidden hit must NOT target the invisible node itself" );
                }
                collapsed.setCollapse( false );

                // (7) centerOnNode really scrolls: with the tree taller than the viewport, stepping to different hits
                //     moves the viewport (not all positions equal) and at least one lands scrolled away from the top.
                if ( !scrollsToHits( frame, tp, found, expected.size(), ok ) ) {
                    // scrollsToHits already recorded the specific failure
                }

                // (8) an empty found set clears the navigator
                tp.setFoundNodes0( null );
                if ( tp.getSearchHitCount() != 0 ) {
                    fail( ok, "clearing the found set should leave 0 hits, got " + tp.getSearchHitCount() );
                }
                tp.stepToFoundNode( 1 );
                if ( tp.getSearchHitIndex() != -1 ) {
                    fail( ok, "stepping with no hits must leave the position unset, got " + tp.getSearchHitIndex() );
                }

                // (9) the ControlPanel glue: a real search (box A) shows the "k / N" row (a real laid-out flip, not
                //     just the visible flag), stepping via the panel updates it, clearing the box hides it again
                cp.getSearchTextField0().setText( "" );
                cp.search0();
                if ( cp.isSearchNavVisibleForTest() ) {
                    fail( ok, "the step-through row should be hidden before any search" );
                }
                cp.getSearchTextField0().setText( "kinase" );
                cp.search0();
                if ( !cp.isSearchNavVisibleForTest() ) {
                    fail( ok, "a search with hits should reveal the step-through row" );
                }
                // the revalidate() in updateSearchHitNavigation lays out asynchronously; force the pass the live EDT
                // would run, then confirm the now-visible GridBag row actually gets a non-zero height (not 0x0)
                forceLayout( cp );
                if ( cp.getSearchNavPanelHeightForTest() <= 0 ) {
                    fail( ok, "the revealed navigator row must be laid out (height>0), got "
                            + cp.getSearchNavPanelHeightForTest() );
                }
                if ( !( "– / " + expected.size() ).equals( cp.getSearchNavLabelForTest() ) ) {
                    fail( ok, "unstepped navigator should read '– / " + expected.size() + "', got '"
                            + cp.getSearchNavLabelForTest() + "'" );
                }
                cp.stepSearchHit( 1 );
                if ( !( "1 / " + expected.size() ).equals( cp.getSearchNavLabelForTest() ) ) {
                    fail( ok, "after one step the navigator should read '1 / " + expected.size() + "', got '"
                            + cp.getSearchNavLabelForTest() + "'" );
                }

                // (10) the View-menu Find Next / Find Previous (⌘G / ⌘⇧G) dispatch a step through MainFrame
                frame.getFindNextHitItemForTest().doClick();
                if ( tp.getSearchHitIndex() != 1 ) {
                    fail( ok, "Find Next menu item should advance to index 1, got " + tp.getSearchHitIndex() );
                }
                frame.getFindPreviousHitItemForTest().doClick();
                if ( tp.getSearchHitIndex() != 0 ) {
                    fail( ok, "Find Previous menu item should return to index 0, got " + tp.getSearchHitIndex() );
                }

                // clearing the search box hides the row
                cp.getSearchTextField0().setText( "" );
                cp.search0();
                if ( cp.isSearchNavVisibleForTest() ) {
                    fail( ok, "clearing the search box should hide the step-through row" );
                }

                // (11) search box B also drives the navigator
                cp.getSearchTextField1().setText( "kinase" );
                cp.search1();
                if ( !cp.isSearchNavVisibleForTest() || ( tp.getSearchHitCount() != expected.size() ) ) {
                    fail( ok, "search box B should reveal the navigator with " + expected.size() + " hits, got "
                            + cp.getSearchNavLabelForTest() );
                }
                cp.getSearchTextField1().setText( "" );
                cp.search1();

                // (12) manual node selection reveals and updates the navigator (README advertises selection hits)
                if ( cp.isSearchNavVisibleForTest() ) {
                    fail( ok, "navigator should be hidden with no search and no selection" );
                }
                final PhylogenyNode a = expected.get( 0 ), b = expected.get( 1 );
                tp.selectNode( a );
                tp.selectNode( b );
                if ( !cp.isSearchNavVisibleForTest() ) {
                    fail( ok, "selecting nodes should reveal the step-through navigator" );
                }
                if ( tp.getSearchHitCount() != 2 ) {
                    fail( ok, "two selected nodes should be 2 hits, got " + tp.getSearchHitCount() );
                }
                tp.selectNode( a ); // toggle a back off
                tp.selectNode( b ); // toggle b back off -> empty -> navigator hidden
                if ( cp.isSearchNavVisibleForTest() ) {
                    fail( ok, "deselecting all nodes should hide the navigator" );
                }
            } );
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            fail( ok, "unexpected: " + t );
        }
        finally {
            if ( mf[ 0 ] != null ) {
                try {
                    SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
                }
                catch ( final Exception ignore ) {
                    // ignore teardown failure
                }
            }
        }
        return ok[ 0 ];
    }

    /** Puts the step-through back to the unpositioned (-1) state without relying on a same-set assignment (which,
     *  post-fix, no longer resets): null then the set is a genuine change. */
    private static void resetPosition( final TreePanel tp, final Set<Long> found ) {
        tp.setFoundNodes0( null );
        tp.setFoundNodes0( new HashSet<>( found ) );
    }

    /** Runs the layout pass the live EDT would run after a revalidate(), so a just-shown/hidden GridBag row gets
     *  its real bounds within a single synchronous invokeAndWait block. */
    private static void forceLayout( final ControlPanel cp ) {
        cp.setSize( cp.getPreferredSize() );
        cp.doLayout();
    }

    /** Makes the tree taller than its scroll viewport, steps to every hit, and asserts centerOnNode moved the
     *  viewport: the recorded scroll positions are not all identical, and at least one hit scrolled off the top. */
    private static boolean scrollsToHits( final MainFrame frame, final TreePanel tp, final Set<Long> found,
                                          final int n, final boolean[] ok ) {
        final JScrollPane sp = frame.getMainPanel().getCurrentScrollPane();
        if ( sp == null ) {
            fail( ok, "no scroll pane to test centering" );
            return false;
        }
        final int vp_w = 300, vp_h = 220, tree_h = 2400;
        sp.setSize( vp_w, vp_h );
        sp.getViewport().setSize( vp_w, vp_h );
        tp.setSize( vp_w, tree_h ); // taller than the viewport -> a real vertical scroll range
        tp.calcParametersForPainting( vp_w, tree_h );
        sp.validate();
        resetPosition( tp, found );
        final Set<Integer> ys = new HashSet<>();
        int max_y = 0;
        for( int i = 0; i < n; i++ ) {
            tp.stepToFoundNode( 1 );
            final Point p = sp.getViewport().getViewPosition();
            ys.add( p.y );
            max_y = Math.max( max_y, p.y );
        }
        // restore the normal painting size for the remaining sub-tests
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        if ( ys.size() < 2 ) {
            fail( ok, "centerOnNode should move the viewport to different y for different hits, got " + ys );
            return false;
        }
        if ( max_y <= 0 ) {
            fail( ok, "centering on a lower hit should scroll the viewport down (max y=" + max_y + ")" );
            return false;
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [SearchStepThroughTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [SearchStepThroughTest] " + msg );
        ok[ 0 ] = false;
    }

    private SearchStepThroughTest() {
    }
}
