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
import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.TIME_AXIS_TYPE;
import org.forester.io.parsers.json.AuspiceJsonParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The Auspice/Nextstrain time&harr;divergence branch-length toggle (Increment 2). Loads the synthetic nCoV Auspice v2
 * demo, then drives the ControlPanel "Branch lengths:" dropdown the way a user would and asserts: the control appears
 * (the tree carries both a date and a nextstrain:div); the default TIME view lays out by year with a CALENDAR axis;
 * switching to DIVERGENCE re-lays out by the (tiny) nextstrain:div with the calendar axis OFF and the unit set to
 * substitutions/site; switching back is exactly reversible; and Reset returns to the TIME view. Headful; a green no-op
 * when headless. Dogfoods {@code nextstrain-ncov.json}.
 */
public final class NextstrainBranchModeTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "Nextstrain branch mode: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = parseAuspice( "nextstrain-ncov.json" );
            if ( phy == null ) {
                return fail( "nextstrain-ncov.json demo missing / unreadable" );
            }
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "ns" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final ControlPanel cp = frame.getMainPanel().getControlPanel();
                    if ( ( tp == null ) || ( cp == null ) ) {
                        fail( ok, "no tree panel / control panel" );
                        return;
                    }
                    // the tree carries both signals -> the toggle applies and its control is shown (via tabChanged)
                    if ( !tp.isNextstrainTimeDivergenceApplicable() ) {
                        fail( ok, "the demo must carry both a date and a nextstrain:div (toggle applicable)" );
                    }
                    if ( !cp.isBranchLengthsControlVisible() ) {
                        fail( ok, "the 'Branch lengths:' control must be visible for an Auspice tree" );
                    }
                    // default = TIME: laid out by year (root-to-tip span > 1) with a CALENDAR axis
                    if ( tp.getNextstrainBranchMode() != TreePanel.NEXTSTRAIN_BRANCH_MODE.TIME ) {
                        fail( ok, "default branch mode must be TIME" );
                    }
                    if ( !"Time".equals( cp.getBranchLengthsSelection() ) ) {
                        fail( ok, "the control must show 'Time' by default, got " + cp.getBranchLengthsSelection() );
                    }
                    // TIME lengths are num_date deltas (year-scale, ~0.8 yr for the nCoV demo); DIVERGENCE lengths are
                    // the tiny nextstrain:div deltas (subs/site, ~0.0024) -- a ~hundreds-fold separation.
                    final double time_span = maxRootDist( tp.getPhylogeny() );
                    if ( time_span < 0.1 ) {
                        fail( ok, "TIME view root-to-tip span must be year-scale (>= 0.1), got " + time_span );
                    }
                    if ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) {
                        fail( ok, "TIME view must derive a CALENDAR axis, got " + tp.effectiveTimeAxisType() );
                    }
                    // capture the TIME layout for a render diff below
                    final Options o = frame.getOptions();
                    final int w = 800, h = 500;
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage img_time = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );

                    // switch to DIVERGENCE via the real dropdown dispatch
                    cp.userSelectBranchLengthsForTest( TreePanel.NEXTSTRAIN_BRANCH_MODE.DIVERGENCE );
                    if ( tp.getNextstrainBranchMode() != TreePanel.NEXTSTRAIN_BRANCH_MODE.DIVERGENCE ) {
                        fail( ok, "selecting Divergence must switch the mode" );
                    }
                    final double div_span = maxRootDist( tp.getPhylogeny() );
                    if ( div_span >= ( time_span / 10.0 ) ) {
                        fail( ok, "DIVERGENCE view must re-lay-out by the (tiny) div, got span " + div_span
                                + " vs time " + time_span );
                    }
                    if ( !"subs/site".equals( tp.getPhylogeny().getDistanceUnit() ) ) {
                        fail( ok, "DIVERGENCE view distance unit must be subs/site, got "
                                + tp.getPhylogeny().getDistanceUnit() );
                    }
                    if ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.NONE ) {
                        fail( ok, "DIVERGENCE view must turn the calendar axis OFF (NONE), got "
                                + tp.effectiveTimeAxisType() );
                    }
                    // the re-layout must be visible: the DIVERGENCE render differs substantially from the TIME render
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage img_div = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    if ( diffPixels( img_time, img_div ) < 500 ) {
                        fail( ok, "the DIVERGENCE layout must differ visibly from the TIME layout" );
                    }

                    // switch back to TIME -> exactly reversible (span + axis + unit restored)
                    cp.userSelectBranchLengthsForTest( TreePanel.NEXTSTRAIN_BRANCH_MODE.TIME );
                    if ( Math.abs( maxRootDist( tp.getPhylogeny() ) - time_span ) > 1e-6 ) {
                        fail( ok, "switching back to TIME must restore the year span exactly" );
                    }
                    if ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) {
                        fail( ok, "switching back to TIME must restore the CALENDAR axis" );
                    }
                    if ( !"year".equals( tp.getPhylogeny().getDistanceUnit() ) ) {
                        fail( ok, "switching back to TIME must restore the year unit" );
                    }

                    // Reset to Defaults path: from DIVERGENCE, reset returns to TIME
                    cp.userSelectBranchLengthsForTest( TreePanel.NEXTSTRAIN_BRANCH_MODE.DIVERGENCE );
                    tp.resetNextstrainBranchModeToDefault();
                    if ( ( tp.getNextstrainBranchMode() != TreePanel.NEXTSTRAIN_BRANCH_MODE.TIME )
                            || ( Math.abs( maxRootDist( tp.getPhylogeny() ) - time_span ) > 1e-6 )
                            || ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) ) {
                        fail( ok, "Reset must return the branch mode to the TIME view" );
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    fail( ok, "exception: " + t );
                }
                finally {
                    frame.dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return fail( "exception: " + t );
        }
    }

    /** Count of pixels that differ between two equal-size renders (the layout change is visible when this is large). */
    private static int diffPixels( final BufferedImage a, final BufferedImage b ) {
        final int w = Math.min( a.getWidth(), b.getWidth() );
        final int h = Math.min( a.getHeight(), b.getHeight() );
        int n = 0;
        for ( int y = 0; y < h; ++y ) {
            for ( int x = 0; x < w; ++x ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Max root-to-tip distance over the external nodes (sum of positive branch lengths), independent of TreePanel. */
    private static double maxRootDist( final Phylogeny phy ) {
        double max = 0;
        for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
            max = Math.max( max, ext.calculateDistanceToRoot() );
        }
        return max;
    }

    private static Phylogeny parseAuspice( final String demo ) {
        final File f = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
        if ( !f.exists() ) {
            return null;
        }
        try {
            final AuspiceJsonParser parser = new AuspiceJsonParser();
            parser.setSource( f );
            final Phylogeny[] phys = parser.parse();
            return ( ( phys == null ) || ( phys.length == 0 ) ) ? null : phys[ 0 ];
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [NextstrainBranchModeTest] " + m );
        return false;
    }

    private static void fail( final boolean[] ok, final String m ) {
        ok[ 0 ] = false;
        fail( m );
    }

    private NextstrainBranchModeTest() {
    }
}
