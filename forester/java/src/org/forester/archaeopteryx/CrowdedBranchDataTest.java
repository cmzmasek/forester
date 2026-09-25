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

import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;

/**
 * Tests "auto-hide crowded data" for BRANCH-anchored marks -- support values, branch-length values and support
 * symbols -- and the diagnostic paint-time readout that shares this corner of the code.
 * <p>
 * The rule is an occupancy test: a mark is drawn unless something already drawn is in its way. What that has to
 * get right is both directions -- a nest of near-zero branches must lose the numbers that would land on top of
 * each other, and a lone short branch in a sparse tree must KEEP its number, because "zero is a value, not an
 * absence" and nothing is in the way.
 */
public final class CrowdedBranchDataTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CrowdedBranchData: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            crowded( ok );
            sparse( ok );
            fpsCounter( ok );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
        return ok[ 0 ];
    }

    // ---- a nest of near-zero branches: the numbers that would collide are dropped ------------------------------
    private static void crowded( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { nest() }, new Configuration(), "crowded" ) );
        final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
        final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
        SwingUtilities.invokeAndWait( () -> {
            cp.setCheckbox( DisplayOption.WRITE_CONFIDENCE_VALUES, true );
            cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
            tp.getOptions().setSupportVisualization( Options.SUPPORT_VISUALIZATION.SIZE_SCALED );
            tp.getOptions().setShowOverview( false );
            tp.setOvOn( false );
            paint( tp, 900, 600 );
        } );
        if ( ( tp.numbersDrawnForTest() + tp.numbersSuppressedForTest() ) < 10 ) {
            fail( ok, "precondition: the nest must offer plenty of numbers, got "
                    + ( tp.numbersDrawnForTest() + tp.numbersSuppressedForTest() ) );
            dispose( mf );
            return;
        }
        if ( tp.numbersSuppressedForTest() < 1 ) {
            fail( ok, "a nest of near-zero branches must lose the numbers that would collide -- none were" );
        }
        if ( tp.numbersDrawnForTest() < 1 ) {
            fail( ok, "...but not ALL of them: the first claim on a free space must still be drawn" );
        }
        if ( tp.symbolsSuppressedForTest() < 1 ) {
            fail( ok, "the support symbols pile up in the same nest and must be thinned too" );
        }
        final int tp_numbers = tp.numbersDrawnForTest();
        final int tp_symbols = tp.symbolsDrawnForTest();

        // Repainting an unchanged view must give the SAME answer. This is what a per-paint reset buys: without
        // it the second paint meets a map still full of the first one's marks, drops everything, and the tree
        // flickers as it redraws.
        final int[] again = new int[ 2 ];
        SwingUtilities.invokeAndWait( () -> {
            paint( tp, 900, 600 );
            again[ 0 ] = tp.numbersDrawnForTest();
            again[ 1 ] = tp.symbolsDrawnForTest();
        } );
        if ( ( again[ 0 ] != tp_numbers ) || ( again[ 1 ] != tp_symbols ) ) {
            fail( ok, "repainting the same view must draw the same marks: " + tp_numbers + "/" + tp_symbols
                    + " then " + again[ 0 ] + "/" + again[ 1 ] );
        }

        // switching auto-hide OFF draws everything again -- the answer to "my figure is missing values"
        final int[] all = new int[ 2 ];
        SwingUtilities.invokeAndWait( () -> {
            cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, false );
            paint( tp, 900, 600 );
            all[ 0 ] = tp.numbersSuppressedForTest();
            all[ 1 ] = tp.symbolsSuppressedForTest();
        } );
        if ( ( all[ 0 ] != 0 ) || ( all[ 1 ] != 0 ) ) {
            fail( ok, "with auto-hide off nothing may be suppressed, got " + all[ 0 ] + " numbers / " + all[ 1 ]
                    + " symbols" );
        }
        dispose( mf );
    }

    // ---- a sparse tree: a zero-length branch KEEPS its number ---------------------------------------------------
    private static void sparse( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { sparseWithAZero() }, new Configuration(), "sparse" ) );
        final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
        final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
        SwingUtilities.invokeAndWait( () -> {
            cp.setCheckbox( DisplayOption.WRITE_CONFIDENCE_VALUES, true );
            cp.setCheckbox( DisplayOption.WRITE_BRANCH_LENGTH_VALUES, true );
            cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
            tp.getOptions().setShowOverview( false );
            tp.setOvOn( false );
            paint( tp, 900, 600 );
        } );
        // THE case a branch-length threshold got wrong: the branch is drawn 0 px, but nothing is near it, so the
        // number reads perfectly and hiding it would claim there is no value where there is one.
        if ( tp.numbersSuppressedForTest() != 0 ) {
            fail( ok, "in a sparse tree nothing collides, so nothing may be hidden -- "
                    + tp.numbersSuppressedForTest() + " number(s) were, including on a zero-length branch" );
        }
        if ( tp.numbersDrawnForTest() < 2 ) {
            fail( ok, "precondition: the sparse tree must actually draw some numbers, got "
                    + tp.numbersDrawnForTest() );
        }
        dispose( mf );
    }

    // ---- the diagnostic readout --------------------------------------------------------------------------------
    private static void fpsCounter( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { sparseWithAZero() }, new Configuration(), "fps" ) );
        final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
        // The SHIPPED default, asked of a fresh Options -- not of this frame's, which legitimately carries
        // whatever the user last chose (the setting is persisted). Reading the live value here would be testing
        // the developer's own settings file.
        if ( Options.createDefaultInstance().isShowFps() ) {
            fail( ok, "the FPS counter must be OFF in a fresh install's defaults" );
        }
        final BufferedImage[] img = new BufferedImage[ 3 ];
        SwingUtilities.invokeAndWait( () -> {
            tp.getOptions().setShowOverview( false );
            tp.setOvOn( false );
            img[ 0 ] = paint( tp, 900, 600 ); // off
            tp.getOptions().setShowFps( true );
            paint( tp, 900, 600 ); // one frame so there is a timing to report
            img[ 1 ] = paint( tp, 900, 600 );
        } );
        if ( tp.fpsReadout() == null ) {
            fail( ok, "after painting, there must be a paint time to report" );
        }
        if ( tp.averagePaintMillis() <= 0 ) {
            fail( ok, "the measured paint time must be positive, got " + tp.averagePaintMillis() );
        }
        if ( !differ( img[ 0 ], img[ 1 ] ) ) {
            fail( ok, "switching the FPS counter on must change what is drawn" );
        }
        // EVERY layout, not just the rectangular family. The counter first shipped inside the block that draws
        // the "Time tree" badge, which is gated on the graphics type -- so circular and unrooted silently had no
        // readout at all. It reports on the PAINT, and every layout paints.
        for( final Options.PHYLOGENY_GRAPHICS_TYPE type : Options.PHYLOGENY_GRAPHICS_TYPE.values() ) {
            final boolean[] shows = { false };
            SwingUtilities.invokeAndWait( () -> {
                tp.setPhylogenyGraphicsType( type );
                tp.getOptions().setShowFps( false );
                final BufferedImage without = paint( tp, 900, 600 );
                tp.getOptions().setShowFps( true );
                paint( tp, 900, 600 );
                shows[ 0 ] = differ( without, paint( tp, 900, 600 ) );
            } );
            if ( !shows[ 0 ] ) {
                fail( ok, "the FPS counter must be drawn in the " + type + " layout too" );
            }
        }
        SwingUtilities.invokeAndWait( () -> {
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp.getOptions().setShowFps( true );
        } );

        // ...and it must never reach an export: a figure cannot ship with a diagnostic on it
        final int[] export_diff = new int[ 1 ];
        SwingUtilities.invokeAndWait( () -> {
            final BufferedImage on = exportImage( tp, 900, 600 );
            tp.getOptions().setShowFps( false );
            final BufferedImage off = exportImage( tp, 900, 600 );
            export_diff[ 0 ] = differ( on, off ) ? 1 : 0;
        } );
        if ( export_diff[ 0 ] != 0 ) {
            fail( ok, "the FPS counter must not appear in an exported figure" );
        }
        dispose( mf );
    }

    // ---- helpers ------------------------------------------------------------------------------------------------
    private static BufferedImage paint( final TreePanel tp, final int w, final int h ) {
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.paintPhylogeny( g, false, false, w, h, 0, 0 );
        g.dispose();
        return img;
    }

    /** The same paint down the EXPORT path (to_graphics_file), which is what an exported figure gets. */
    private static BufferedImage exportImage( final TreePanel tp, final int w, final int h ) {
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.paintPhylogeny( g, false, true, w, h, 0, 0 );
        g.dispose();
        return img;
    }

    private static boolean differ( final BufferedImage a, final BufferedImage b ) {
        if ( ( a == null ) || ( b == null ) ) {
            return true;
        }
        for( int x = 0; x < a.getWidth(); ++x ) {
            for( int y = 0; y < a.getHeight(); ++y ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    /** A caterpillar of 40 internal nodes on near-zero branches, each carrying support -- the "nest". */
    private static Phylogeny nest() {
        final PhylogenyNode root = new PhylogenyNode();
        PhylogenyNode cursor = root;
        for( int i = 0; i < 40; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "t" + i );
            tip.setDistanceToParent( 0.4 );
            final PhylogenyNode inner = new PhylogenyNode();
            inner.setDistanceToParent( 0.0001 ); // the nest: essentially no length
            inner.getBranchData().addConfidence( new Confidence( 50 + ( i % 50 ), "bootstrap" ) );
            cursor.addAsChild( tip );
            cursor.addAsChild( inner );
            cursor = inner;
        }
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "last_a" );
        a.setDistanceToParent( 0.4 );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "last_b" );
        b.setDistanceToParent( 0.4 );
        cursor.addAsChild( a );
        cursor.addAsChild( b );
        return wrap( root );
    }

    /** Four tips, well spread, one of them on a zero-length branch whose parent carries support. */
    private static Phylogeny sparseWithAZero() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode left = new PhylogenyNode();
        left.setDistanceToParent( 0.3 );
        left.getBranchData().addConfidence( new Confidence( 0, "bootstrap" ) );
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "A" );
        a.setDistanceToParent( 0.0 ); // a real, measured zero
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "B" );
        b.setDistanceToParent( 0.2 );
        left.addAsChild( a );
        left.addAsChild( b );
        final PhylogenyNode right = new PhylogenyNode();
        right.setDistanceToParent( 0.6 );
        right.getBranchData().addConfidence( new Confidence( 95, "bootstrap" ) );
        final PhylogenyNode c = new PhylogenyNode();
        c.setName( "C" );
        c.setDistanceToParent( 0.4 );
        final PhylogenyNode d = new PhylogenyNode();
        d.setName( "D" );
        d.setDistanceToParent( 0.5 );
        right.addAsChild( c );
        right.addAsChild( d );
        root.addAsChild( left );
        root.addAsChild( right );
        return wrap( root );
    }

    private static Phylogeny wrap( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static void dispose( final MainFrame[] mf ) throws Exception {
        SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [CrowdedBranchDataTest] " + message );
        ok[ 0 ] = false;
    }

    private CrowdedBranchDataTest() {
        // not instantiable
    }
}
