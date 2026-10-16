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

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * A tree opened straight into the circular layout -- the saved default display type is CIRCULAR, as a user who
 * last closed the program in that layout has it -- is painted BEFORE calcParametersForPainting or the control
 * panel ever ran. Until 0.11.182 that first paint divided by the never-computed radial depth, 0, which put every
 * internal node at infinity; the branch-obstacle grid then registered each fork arc chord by chord, two billion
 * chords at infinity per arc -- the OutOfMemoryError of 2026-10-16 (0.11.180) and the hang that replaced it when
 * the grid went linear (0.11.181). The depth is now computed on first use, and this test holds the exact state
 * of that first paint: the real start-up path with the saved default, nothing computed yet, a cladogram.
 * <p>
 * Also pins that the depth the first paint uses is the one every later layout computation derives (both
 * calcParametersForPainting and the control panel's calcMaxDepth: the deepest displayed tip's depth plus one,
 * a collapsed clade counted as a tip), where until 0.11.182 the two derivations differed and a circular
 * cladogram's internal nodes moved between a resize and a control-panel change.
 */
public final class CircularFirstPaintTest {

    public static void main( final String[] args ) {
        System.out.println( test() ? "CircularFirstPaintTest: OK." : "CircularFirstPaintTest: FAILED." );
    }

    public static boolean test() {
        if ( java.awt.GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final String dir = System.getProperty( GuiPreferences.DIR_PROPERTY );
        if ( ( dir == null ) || dir.isEmpty() ) {
            // never write into the user's own ~/.archaeopteryx: the suite isolates the directory (Test.main), a
            // standalone run must pass -Darchaeopteryx.cache.dir=<tmp>
            return fail( "the settings directory must be isolated: set -D" + GuiPreferences.DIR_PROPERTY );
        }
        final File settings = new File( dir, GuiPreferences.SETTINGS_FILE );
        byte[] previous = null;
        try {
            if ( settings.exists() ) {
                previous = Files.readAllBytes( settings.toPath() );
            }
            Files.createDirectories( settings.toPath().getParent() );
            Files.write( settings.toPath(), "phylogeny_graphics_type=CIRCULAR\n".getBytes( StandardCharsets.ISO_8859_1 ) );
            return firstPaint();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected: " + e );
        }
        finally {
            try {
                if ( previous == null ) {
                    settings.delete();
                }
                else {
                    Files.write( settings.toPath(), previous );
                }
            }
            catch ( final Exception e ) {
                e.printStackTrace();
            }
        }
    }

    /** The exact state of the first paint of a tab opened circular, then the paint itself. */
    private static boolean firstPaint() throws Exception {
        final Phylogeny phy = fixture();
        final MainFrame[] mf = new MainFrame[ 1 ];
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> {
                mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(), "first paint" );
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                // the state of the first paint: the saved default opened the tab circular, the tab lays out as a
                // cladogram (the branch the division is in), and nothing has computed the depth. The construction
                // computes it before any paint can run -- since 0.11.183 on every path, the window being built on the
                // event thread whoever asks (raceAsStarted); until then Archaeopteryx.main built it on its own thread
                // and the first paint came first (measured 2026-10-16). The never-computed state is put back
                // deliberately here: the lazy depth is the second line, for any paint that meets it
                if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                    ok[ 0 ] = fail( "fixture: the saved default must open the tab circular, got " + tp.getPhylogenyGraphicsType() );
                }
                if ( tp.isPhyHasBranchLengths() || tp.getControlPanel().isDrawPhylogram() ) {
                    ok[ 0 ] = fail( "fixture: the tab must lay out as a cladogram" );
                }
                if ( ( tp.getWidth() <= 0 ) || ( tp.getHeight() <= 0 ) ) {
                    ok[ 0 ] = fail( "fixture: the panel must have a size to paint, got " + tp.getSize() );
                }
                tp.forgetCircMaxDepthForTest();
                if ( tp.storedCircMaxDepthForTest() != 0 ) {
                    ok[ 0 ] = fail( "fixture: the depth must read as never computed, got " + tp.storedCircMaxDepthForTest() );
                }
                if ( !ok[ 0 ] ) {
                    return;
                }
                final PhylogenyNode inner = phy.getNode( "abcd" );
                final double f = tp.circularRadiusFractionForTest( inner );
                if ( !( ( f >= 0 ) && ( f <= 1 ) ) ) { // a NaN fails too
                    ok[ 0 ] = fail( "an internal node's ring fraction before any layout computation must lie in [0, 1], got " + f );
                    return; // before 0.11.182 the paint below ran for minutes
                }
                final int expected = PhylogenyMethods.calculateMaxDepthConsiderCollapsed( phy ) + 1;
                if ( tp.storedCircMaxDepthForTest() != expected ) {
                    ok[ 0 ] = fail( "asking for the fraction computes the depth, the deepest tip's plus one: expected " + expected
                            + ", got " + tp.storedCircMaxDepthForTest() );
                }
                final double expected_f = inner.calculateDepth() / (double) expected;
                if ( Math.abs( f - expected_f ) > 1e-12 ) {
                    ok[ 0 ] = fail( "the fraction is depth / that: expected " + expected_f + ", got " + f );
                }
                // the paint itself, in that state: every node lands on the canvas
                tp.forgetCircMaxDepthForTest();
                tp.paintImmediately( 0, 0, tp.getWidth(), tp.getHeight() );
                if ( !onCanvas( tp, phy ) ) {
                    ok[ 0 ] = false;
                }
                // one depth, however it is derived: the layout computation and the control panel agree with the
                // first paint -- on a tree with collapsed clades, where the old derivations differed
                tp.calcParametersForPainting( tp.getWidth(), tp.getHeight() );
                if ( tp.storedCircMaxDepthForTest() != expected ) {
                    ok[ 0 ] = fail( "the layout computation must derive the depth the first paint used: " + expected + ", got "
                            + tp.storedCircMaxDepthForTest() );
                }
                // two collapsed clades on the deepest path: with one, the collapse-aware depth plus one equals the
                // plain depth, and a derivation using the plain depth would pass (a mutant did, 2026-10-16)
                phy.getNode( "abc" ).setCollapse( true );
                phy.getNode( "abcde" ).setCollapse( true );
                final int collapsed = PhylogenyMethods.calculateMaxDepthConsiderCollapsed( phy ) + 1;
                final int plain = PhylogenyMethods.calculateMaxDepth( phy );
                if ( ( collapsed == plain ) || ( collapsed == ( plain + 1 ) ) ) {
                    ok[ 0 ] = fail( "precondition: the expected depth " + collapsed + " must differ from the plain depth " + plain
                            + " and from it plus one" );
                }
                tp.calcMaxDepth();
                final int by_panel = tp.storedCircMaxDepthForTest();
                tp.calcParametersForPainting( tp.getWidth(), tp.getHeight() );
                if ( ( by_panel != collapsed ) || ( tp.storedCircMaxDepthForTest() != collapsed ) ) {
                    ok[ 0 ] = fail( "with collapsed clades the control panel's depth (" + by_panel + ") and the layout computation's ("
                            + tp.storedCircMaxDepthForTest() + ") must both be the deepest displayed tip's plus one, " + collapsed );
                }
            } );
        }
        finally {
            dispose( mf[ 0 ] );
        }
        return ok[ 0 ] & raceAsStarted();
    }

    /** The real path: Archaeopteryx.main asks for the window from ITS thread. Until 0.11.183 the window was built
     *  there, became visible while its construction was still running, and the event thread's first paint raced the
     *  layout computation -- and won, on 2026-10-16. Asked from off the event thread here, the window must have been
     *  CONSTRUCTED on it (the construction then completes before any paint), the event thread must be free again
     *  within seconds, and every node must be on the canvas. firstPaint holds the lazy depth, the second line. */
    private static boolean raceAsStarted() throws Exception {
        if ( SwingUtilities.isEventDispatchThread() ) {
            return fail( "fixture: the window must be asked for off the event thread" );
        }
        final Phylogeny phy = fixture();
        final MainFrame mf = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(), "as started" );
        boolean ok = true;
        try {
            if ( !( (MainFrameApplication) mf ).builtOnEventThreadForTest() ) {
                ok = fail( "a window asked for off the event thread must be constructed on it" );
            }
            final Thread waiter = new Thread( () -> {
                try {
                    SwingUtilities.invokeAndWait( () -> {} );
                }
                catch ( final Exception e ) {
                    throw new RuntimeException( e );
                }
            }, "first-paint-watch" );
            waiter.setDaemon( true );
            waiter.start();
            waiter.join( 30_000 );
            if ( waiter.isAlive() ) {
                return fail( "the first paint of a window asked for the way Archaeopteryx.main asks hangs the event thread" );
            }
            final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
            final boolean[] on = { true };
            SwingUtilities.invokeAndWait( () -> {
                if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                    on[ 0 ] = fail( "fixture: the saved default must open the tab circular, got " + tp.getPhylogenyGraphicsType() );
                }
                tp.paintImmediately( 0, 0, tp.getWidth(), tp.getHeight() );
                on[ 0 ] &= onCanvas( tp, phy );
                if ( tp.storedCircMaxDepthForTest() <= 0 ) {
                    on[ 0 ] = fail( "after the first paints the depth must be computed, got " + tp.storedCircMaxDepthForTest() );
                }
            } );
            ok &= on[ 0 ];
        }
        finally {
            dispose( mf );
        }
        return ok;
    }

    /** Deeper on one side, no branch lengths: a cladogram, whose internal nodes sit by DEPTH. */
    private static Phylogeny fixture() throws Exception {
        return ParserBasedPhylogenyFactory.getInstance()
                .create( "(((((a,b)ab,c)abc,d)abcd,e)abcde,(f,g)fg)root;", new NHXParser() )[ 0 ];
    }

    /** Every node within reach of the ring after a paint (an infinite or NaN coordinate fails). */
    private static boolean onCanvas( final TreePanel tp, final Phylogeny phy ) {
        boolean ok = true;
        final double reach = 4.0 * tp.radialDiameter();
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !( ( Math.abs( n.getXcoord() ) <= reach ) && ( Math.abs( n.getYcoord() ) <= reach ) ) ) {
                ok = fail( "after the first paint node " + n.getName() + " must be on the canvas, is at " + n.getXcoord() + ", "
                        + n.getYcoord() );
            }
        }
        return ok;
    }

    private static void dispose( final MainFrame mf ) throws Exception {
        if ( mf != null ) {
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf ).dispose() );
        }
    }

    private static boolean fail( final String what ) {
        System.out.println( "  [CircularFirstPaintTest] " + what );
        return false;
    }

    private CircularFirstPaintTest() {
    }
}
