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

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Two time-tree UX refinements. (b) The data-driven global overlay toggles (Node Age HPD bars, Fossil Range bars) are
 * reconciled to whether ANY open tree carries the interval data — so they enable for a dated/fossil tree, a still-open
 * dated tree keeps them on, and a session with none clears them (bidirectional, non-regressing). (a) A dated tree
 * shown as a CLADOGRAM draws a faint screen-only hint (the time axis needs a phylogram) rather than silently omitting
 * the axis. Headful; a green no-op when headless.
 */
public final class TimeOverlayUxTest {

    private static final int W = 760, H = 520;

    public static void main( final String[] args ) {
        System.out.println( "Time overlay UX: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        if ( java.awt.GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny hpd = parse( "node-hpd-bars.xml" );      // dated, internal HPD intervals (unit mya)
            final Phylogeny fossil = parse( "fossil-range-bars.xml" ); // fossil TIP ranges (external intervals)
            final Phylogeny plain = parse( "color-by-property.xml" );  // no dates
            if ( ( hpd == null ) || ( fossil == null ) || ( plain == null ) ) {
                return fail( "demo trees missing" );
            }
            final boolean[] ok = { true };

            // (b) reconcile: a dated-interval tree -> HPD on, fossil off
            withFrame( new Phylogeny[] { copy( hpd ) }, ok, f -> {
                if ( !f.getMainPanel().getOptions().isShowHpdBars()
                        || f.getMainPanel().getOptions().isShowFossilRangeBars() ) {
                    fail( ok, "a node-HPD tree must enable HPD bars only" );
                }
            } );
            // a fossil-range tree -> fossil on, HPD off
            withFrame( new Phylogeny[] { copy( fossil ) }, ok, f -> {
                if ( !f.getMainPanel().getOptions().isShowFossilRangeBars()
                        || f.getMainPanel().getOptions().isShowHpdBars() ) {
                    fail( ok, "a fossil-range tree must enable Fossil bars only" );
                }
            } );
            // dated + plain -> HPD still on (a plain tab does NOT turn it off = non-regression)
            withFrame( new Phylogeny[] { copy( hpd ), copy( plain ) }, ok, f -> {
                if ( !f.getMainPanel().getOptions().isShowHpdBars() ) {
                    fail( ok, "a plain second tab must not clear HPD while a dated tab is open" );
                }
            } );
            // bidirectional CLEAR: a stuck-on toggle is cleared when no open tree has the data (the fix vs. ON-only)
            withFrame( new Phylogeny[] { copy( plain ) }, ok, f -> {
                f.getMainPanel().getOptions().setShowHpdBars( true );
                f.getMainPanel().getOptions().setShowFossilRangeBars( true );
                AptxUtil.reconcileDataDrivenOverlayToggles( f.getMainPanel() );
                if ( f.getMainPanel().getOptions().isShowHpdBars()
                        || f.getMainPanel().getOptions().isShowFossilRangeBars() ) {
                    fail( ok, "reconcile must CLEAR the toggles when no open tree has interval data" );
                }
            } );

            // (a) the cladogram hint: a dated tree shown as a CLADOGRAM with the auto (Geologic) axis draws the hint;
            // the same cladogram with the axis explicitly OFF does not -- the ONLY difference is the hint text
            withFrame( new Phylogeny[] { copy( hpd ) }, ok, f -> {
                final TreePanel tp = f.getMainPanel().getCurrentTreePanel();
                tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM );
                tp.setTimeAxisType( null ); // auto -> GEOLOGIC (dated) -> hint fires on a cladogram
                final BufferedImage with_hint = renderScreen( tp );
                tp.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE ); // axis off -> no hint
                final BufferedImage no_hint = renderScreen( tp );
                if ( diffPixels( with_hint, no_hint ) < 30 ) {
                    fail( ok, "a dated CLADOGRAM must draw the 'display as a phylogram' hint (diff "
                            + diffPixels( with_hint, no_hint ) + ")" );
                }
            } );

            // (a, circular) the hint ALSO shows in the CIRCULAR layout: a dated circular cladogram has no ring axis
            withFrame( new Phylogeny[] { copy( hpd ) }, ok, f -> {
                final TreePanel tp = f.getMainPanel().getCurrentTreePanel();
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM );
                tp.setTimeAxisType( null ); // auto -> GEOLOGIC -> hint fires on a circular cladogram
                final BufferedImage with_hint = renderScreen( tp );
                tp.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE ); // axis off -> no hint
                final BufferedImage no_hint = renderScreen( tp );
                if ( diffPixels( with_hint, no_hint ) < 30 ) {
                    fail( ok, "a dated CIRCULAR cladogram must draw the hint (diff "
                            + diffPixels( with_hint, no_hint ) + ")" );
                }
            } );

            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return fail( "exception: " + t );
        }
    }

    private interface FrameCheck {
        void run( MainFrame f ) throws Exception;
    }

    private static void withFrame( final Phylogeny[] phys, final boolean[] ok, final FrameCheck check )
            throws Exception {
        final Configuration conf = new Configuration();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication.createInstance( phys, conf, "tou" ) );
        SwingUtilities.invokeAndWait( () -> {
            try {
                check.run( mf[ 0 ] );
            }
            catch ( final Throwable t ) {
                t.printStackTrace();
                fail( ok, "exception: " + t );
            }
            finally {
                mf[ 0 ].dispose();
            }
        } );
    }

    private static BufferedImage renderScreen( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, W, H );
        tp.paintPhylogeny( g, false, false, 0, 0, 0, 0 ); // screen path -- the hint is screen-only
        g.dispose();
        return img;
    }

    private static int diffPixels( final BufferedImage a, final BufferedImage b ) {
        int n = 0;
        for ( int y = 0; y < H; ++y ) {
            for ( int x = 0; x < W; ++x ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static Phylogeny copy( final Phylogeny phy ) {
        return phy.copy();
    }

    private static Phylogeny parse( final String demo ) {
        final File f = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
        if ( !f.exists() ) {
            return null;
        }
        try {
            return ParserBasedPhylogenyFactory.getInstance().create( f, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [TimeOverlayUxTest] " + m );
        return false;
    }

    private static void fail( final boolean[] ok, final String m ) {
        ok[ 0 ] = false;
        fail( m );
    }

    private TimeOverlayUxTest() {
    }
}
