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

import org.forester.archaeopteryx.Options.PHYLOGENY_DISPLAY_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Picking the Triangular tree style nudges the tab to a CLADOGRAM (triangular draws straight chords to a clade's
 * extreme tips, so it only reads cleanly with aligned tips / a near-clock tree). The user stays free to switch back to
 * a phylogram (P), and switching to a NON-triangular style does not force a cladogram. Headful; green no-op headless.
 */
public final class TriangularNudgeTest {

    public static void main( final String[] args ) {
        System.out.println( "Triangular nudge: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phylo = parseXml( "scale-axis.xml" ); // has branch lengths -> opens as a phylogram
            final Phylogeny nobl = parseNh( "((a,b),(c,d));" );    // no branch lengths -> always a cladogram
            if ( ( phylo == null ) || ( nobl == null ) ) {
                return fail( "test trees missing" );
            }
            final boolean[] ok = { true };

            withFrame( phylo.copy(), ok, f -> {
                final ControlPanel cp = f.getMainPanel().getControlPanel();
                final TreePanel tp = f.getMainPanel().getCurrentTreePanel();
                // pin a NON-triangular start style so the test is independent of the developer's persisted
                // ~/.archaeopteryx graphics type (a standalone run inherits it; the full suite isolates the cache dir)
                tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                cp.setTreeDisplayType( PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM ); // start as a phylogram
                if ( !cp.isDrawPhylogram() ) {
                    fail( ok, "the branch-length tree must start as a phylogram" );
                }
                // picking Triangular nudges to a cladogram
                f.typeChanged( f._triangular_type_cbmi );
                if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR ) {
                    fail( ok, "the style must become TRIANGULAR" );
                }
                if ( cp.isDrawPhylogram() ) {
                    fail( ok, "Triangular must nudge the tab to a CLADOGRAM" );
                }
                // the user remains free to switch back to a phylogram
                cp.setTreeDisplayType( PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                if ( !cp.isDrawPhylogram() ) {
                    fail( ok, "the user must be able to switch back to P after the nudge" );
                }
                // re-selecting the ALREADY-current Triangular style must NOT re-clobber that deliberate P choice
                // (the nudge is gated on a real transition to Triangular, not merely on the current style)
                f.typeChanged( f._triangular_type_cbmi );
                if ( !cp.isDrawPhylogram() ) {
                    fail( ok, "re-selecting the current Triangular style must not force a cladogram again" );
                }
                // switching to a NON-triangular style must NOT force a cladogram (stays a phylogram)
                f.typeChanged( f._rectangular_type_cbmi );
                if ( !cp.isDrawPhylogram() ) {
                    fail( ok, "switching to Rectangular must NOT force a cladogram" );
                }
            } );

            // a no-branch-length tree is already a cladogram: Triangular is a harmless no-op (no crash)
            withFrame( nobl.copy(), ok, f -> {
                final ControlPanel cp = f.getMainPanel().getControlPanel();
                f.typeChanged( f._triangular_type_cbmi );
                if ( cp.isDrawPhylogram() ) {
                    fail( ok, "a no-branch-length tree must remain a cladogram under Triangular" );
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

    private static void withFrame( final Phylogeny phy, final boolean[] ok, final FrameCheck check ) throws Exception {
        final Configuration conf = new Configuration();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "tri" ) );
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

    private static Phylogeny parseXml( final String demo ) {
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

    private static Phylogeny parseNh( final String nh ) {
        try {
            return ParserBasedPhylogenyFactory.getInstance().create( nh, new NHXParser() )[ 0 ];
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [TriangularNudgeTest] " + m );
        return false;
    }

    private static void fail( final boolean[] ok, final String m ) {
        ok[ 0 ] = false;
        fail( m );
    }

    private TriangularNudgeTest() {
    }
}
