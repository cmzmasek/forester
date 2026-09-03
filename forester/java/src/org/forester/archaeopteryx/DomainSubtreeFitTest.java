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

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Guards the fix for "domains clipped after entering a subtree and pressing F (fit)". A phylogram draws its root
 * shifted right by the root's own branch length, and getMaxDistanceToRoot()/the aligned domain column both fold that
 * branch in -- so a SUBTREE (whose root inherits a nonzero branch from its former parent) shifted the aligned domain
 * column right by rootBranch*xcorr and clipped the domains off the near edge. The fix draws a subtree root as a fixed
 * short stub (not to scale) and corrects the domain-column double-count -- and, since it rides {@code _orientation_R},
 * it holds in all three rectangular orientations (the domain column X is a logical value). Headful; a green no-op when
 * headless. Dogfoods forester/demo/domain-architectures.xml (tips renamed long so the aligned domain column reaches
 * the right edge).
 */
public final class DomainSubtreeFitTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DomainSubtreeFit: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        // break OFF (the reported case) AND break ON (breakCappedHeight also folds in the root branch, so the clip
        // must not resurface); plus a full-tree genuine root edge (the anchor double-count, which the stub can't cover)
        return subtreeDomainsFitOk( false ) && subtreeDomainsFitOk( true ) && fullTreeRootEdgeOk()
                && domainZoomGrowsTheViewOk();
    }

    /** Load the domain demo and give every tip a long label so the aligned domain column reaches the right edge. */
    private static Phylogeny loadLongLabelledDomainTree() throws Exception {
        final File file = new File( System.getProperty( "user.dir" ), "forester/demo/domain-architectures.xml" );
        if ( !file.exists() ) {
            throw new IllegalStateException( "demo tree missing: " + file.getAbsolutePath() );
        }
        final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
        final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
        for ( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.setName( "very_long_isolate_label_that_fills_the_available_width_XXXXXXXXXXXX_" + n.getName() );
        }
        return phy;
    }

    /** Entering a subtree whose root carries an inherited branch length: the root draws as a fixed stub, the depth
     *  cache excludes that branch, and the aligned domain column fits the preferred width. Run with Break Long
     *  Branches OFF (the reported case) and ON (breakCappedHeight also folds in the root branch). */
    private static boolean subtreeDomainsFitOk( final boolean break_long_branches ) {
        try {
            final Phylogeny phy = loadLongLabelledDomainTree();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "domsub" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainFrame frame = mf[ 0 ];
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    // pin the display state (a standalone run must not inherit the developer's persisted ~/.archaeopteryx
                    // graphics type / orientation -- the documented headful-test gotcha)
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    o.setBreakLongBranches( break_long_branches );
                    if ( !tp.getControlPanel().isShowDomainArchitectures() ) {
                        fail( ok, "the demo tree should auto-enable Domain Architectures" );
                    }
                    PhylogenyNode clade = null;
                    for ( final PhylogenyNode c : tp.getPhylogeny().getRoot().getDescendants() ) {
                        if ( !c.isExternal() ) {
                            clade = c;
                            break;
                        }
                    }
                    if ( clade == null ) {
                        fail( ok, "no internal clade to enter as a subtree" );
                        return;
                    }
                    tp.subTree( clade ); // sets up the subtree (its root inherits a branch length) + fits
                    final double root_branch = tp.getPhylogeny().getRoot().getDistanceToParent();
                    if ( !( root_branch > 0 ) ) {
                        fail( ok, "the subtree root must carry an inherited branch length to reproduce, got " + root_branch );
                        return;
                    }
                    final int w = 820, h = 600;
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    tp.resetPreferredSize();
                    AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ); // paint -> assign rootX + coords

                    final String tag = "[break=" + break_long_branches + "] ";
                    // 1. the subtree root draws as a FIXED short stub (MOVE + xdist), NOT scaled by its inherited branch
                    if ( Math.abs( tp.rootXcoordForTest() - tp.stubRootXForTest() ) > 2f ) {
                        fail( ok, tag + "subtree root must draw as a fixed stub at " + tp.stubRootXForTest()
                                + ", got rootX=" + tp.rootXcoordForTest() );
                    }
                    // 2. the depth cache getMaxDistanceToRoot() EXCLUDES the inherited (stubbed) root branch (break OFF)
                    if ( !break_long_branches ) {
                        final double dmax = PhylogenyMethods.calculateMaxDistanceToRoot( tp.getPhylogeny() );
                        if ( Math.abs( tp.getMaxDistanceToRootForTest() - dmax ) > 1e-4 ) {
                            fail( ok, tag + "subtree getMaxDistanceToRoot must exclude the root branch (=" + dmax
                                    + "), got " + tp.getMaxDistanceToRootForTest() );
                        }
                    }
                    // 3. the aligned DOMAIN column fits within the preferred width (the reported clip -- it must NOT
                    //    resurface with Break Long Branches on, where breakCappedHeight also folds in the root branch)
                    final int pref_w = tp.getPreferredSize().width;
                    if ( tp.alignedDomainColumnRightEdgeForTest() > ( pref_w + 1 ) ) {
                        fail( ok, tag + "domains must fit: aligned column right edge "
                                + tp.alignedDomainColumnRightEdgeForTest() + " exceeds preferred width " + pref_w );
                    }
                    // 4. the depth height feeding the extent excludes the stubbed root branch (no under-fill), break OFF
                    if ( !break_long_branches ) {
                        final double h_excl = tp.getPhylogeny().calculateHeight( !o.isCollapsedWithAverageHeigh() )
                                - root_branch;
                        if ( Math.abs( tp.displayedTreeHeightForTest() - h_excl ) > 1e-4 ) {
                            fail( ok, tag + "subtree depth height must exclude the root branch (" + h_excl + "), got "
                                    + tp.displayedTreeHeightForTest() );
                        }
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "threw: " + t );
                }
            } );
            SwingUtilities.invokeAndWait( () -> mf[ 0 ].dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            return fail( "threw: " + t );
        }
    }

    /** A FULL tree whose ROOT carries a genuine branch length (a root edge, drawn to scale) plus domains: the aligned
     *  column must still fit -- guards the domain-column double-count fix, which the subtree stub does NOT cover. */
    /**
     * The "+" domain zoom must GROW the scrollable view, not just the drawing.
     * <p>
     * Widening the domain track widens the label/track reservation, and that reservation IS the panel's width --
     * so if the preferred size does not follow it, the domains simply extend past the right edge with nothing to
     * scroll to. Reported by a user: "the size of the view does not get updated and i can not scroll along the
     * domains". Drives the REAL "+" button, because that is the path the user takes.
     */
    private static boolean domainZoomGrowsTheViewOk() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            final Phylogeny phy = loadLongLabelledDomainTree();
            if ( ( phy == null ) || phy.isEmpty() ) {
                return fail( "could not load domain-architectures.xml" );
            }
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "domzoom" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainPanel mp = mf[ 0 ].getMainPanel();
                final ControlPanel cp = mp.getControlPanel();
                final TreePanel tp = mp.getCurrentTreePanel();
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.setSize( 900, 460 );
                mf[ 0 ].showWhole();
                AptxUtil.renderPhylogenyToImage( 900, 460, tp, tp.getOptions(), false, 1, false );
                final double start_w = tp.domainStructureWidthForTest();
                final int start_pref = tp.getPreferredSize().width;
                for( int i = 1; i <= 4; ++i ) {
                    cp.actionPerformed( new java.awt.event.ActionEvent( cp.zoomInDomainButtonForTest(), 0, "+" ) );
                    AptxUtil.renderPhylogenyToImage( 900, 460, tp, tp.getOptions(), false, 1, false );
                    final int pref_w = tp.getPreferredSize().width;
                    if ( tp.alignedDomainColumnRightEdgeForTest() > ( pref_w + 1 ) ) {
                        fail( ok, "after " + i + " domain zoom-in click(s) the domains reach "
                                + Math.round( tp.alignedDomainColumnRightEdgeForTest() - pref_w )
                                + "px past the scrollable width (" + pref_w + ") -- the view did not grow with them" );
                        return;
                    }
                }
                if ( tp.domainStructureWidthForTest() <= start_w ) {
                    fail( ok, "precondition: the domain zoom did not widen the track at all" );
                }
                if ( tp.getPreferredSize().width <= start_pref ) {
                    fail( ok, "the scrollable width must GROW with the domain track (" + start_pref + " -> "
                            + tp.getPreferredSize().width + ")" );
                }
            } );
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    private static boolean fullTreeRootEdgeOk() {
        try {
            final Phylogeny phy = loadLongLabelledDomainTree();
            phy.getRoot().setDistanceToParent( 0.5 ); // a genuine (large) root edge on the full tree
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "domedge" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainFrame frame = mf[ 0 ];
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    tp.recalculateMaxDistanceToRoot();
                    final int w = 820, h = 600;
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    tp.resetPreferredSize();
                    AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    // the root edge IS drawn to scale here (not a subtree), so rootX is shifted right; the aligned
                    // domain column must still fit (the anchor subtracts the folded root branch once)
                    final int pref_w = tp.getPreferredSize().width;
                    if ( tp.alignedDomainColumnRightEdgeForTest() > ( pref_w + 1 ) ) {
                        fail( ok, "full-tree root-edge: domains must fit: right edge "
                                + tp.alignedDomainColumnRightEdgeForTest() + " exceeds preferred width " + pref_w );
                    }
                    // the numeric scale axis reaches EXACTLY the deepest tip (rootX + span*corr), not a root-branch
                    // past it -- so it agrees with the scale grid (both subtract the folded root edge)
                    double deepest_tip_x = 0;
                    for ( final PhylogenyNodeIterator it = tp.getPhylogeny().iteratorExternalForward(); it.hasNext(); ) {
                        deepest_tip_x = Math.max( deepest_tip_x, it.next().getXcoord() );
                    }
                    final double axis_reach = tp.rootXcoordForTest()
                            + ( tp.numericScaleAxisMaxDist() * tp.getXcorrectionFactor() );
                    if ( Math.abs( axis_reach - deepest_tip_x ) > 2.0 ) {
                        fail( ok, "full-tree root-edge: scale axis reach " + axis_reach + " must equal the deepest tip "
                                + deepest_tip_x + " (no root-branch overshoot)" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "threw: " + t );
                }
            } );
            SwingUtilities.invokeAndWait( () -> mf[ 0 ].dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            return fail( "threw: " + t );
        }
    }

    private static boolean fail( final boolean[] ok, final String msg ) {
        ok[ 0 ] = false;
        return fail( msg );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [DomainSubtreeFitTest] " + msg );
        return false;
    }

    private DomainSubtreeFitTest() {
    }
}
