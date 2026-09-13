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
import java.awt.event.ActionEvent;
import java.io.File;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The domain-track width buttons (d+ / d-) work in the circular and unrooted layouts. The radial layouts keep a width
 * of their own, starting at the smaller of the rectangular width and a fifth of the radius (so the first sight is
 * unchanged), and a press steps whichever width the current layout draws with -- the rectangular width is untouched by
 * a radial press and vice versa. Until 2026-09-13 the rectangular width was capped at that fifth on every redraw, so
 * from a fresh load (255 px against a cap of 71) three d+ presses changed nothing and d- did nothing for nine presses.
 * Christian: "follow the JS fix" (Archaeopteryx.js 1439b5a).
 */
public final class DomainRadialWidthTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DomainRadialWidth: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // needs a display
        }
        try {
            return buttonsStepTheRadialWidth( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR )
                    && buttonsStepTheRadialWidth( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
        }
        catch ( final Throwable e ) {
            e.printStackTrace( System.out );
            return false;
        }
    }

    private static boolean buttonsStepTheRadialWidth( final Options.PHYLOGENY_GRAPHICS_TYPE layout ) throws Exception {
        final File f = new File( System.getProperty( "user.dir" ), "forester/demo/domain-architectures.xml" );
        final Phylogeny phy = ParserUtils.readPhylogenies( f )[ 0 ];
        final MainFrame[] mf = new MainFrame[ 1 ];
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(), "radialwidth" );
            mf[ 0 ].setSize( 1200, 800 );
            mf[ 0 ].validate();
        } );
        try {
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final ControlPanel cp = tp.getControlPanel();
                cp.showDomainArchitecturesFitted(); // the load: the rectangular width = a quarter of the viewport
                final double rect_w = tp.domainStructureWidthForTest();
                mf[ 0 ].getOptions().setPhylogenyGraphicsType( layout );
                tp.setPhylogenyGraphicsType( layout );
                mf[ 0 ].enableRadialLabelsIfDomainsInRadialLayout();
                cp.displayedPhylogenyMightHaveChanged( true );
                // the first sight is unchanged: the radial width starts at min(rectangular width, a fifth of the radius)
                final double cap = 0.2 * tp.radialDiameter() / 2.0;
                final double start = tp.effectiveDomainStructureWidthForTest();
                check( ok, layout + ": precondition, the rectangular width (" + rect_w + ") exceeds the radial cap (" + cap
                        + ") so the old cap would have bitten", rect_w > cap );
                check( ok, layout + ": the radial width starts at the cap, got " + start + " vs " + cap,
                       Math.abs( start - cap ) < 0.01 );
                final float f0 = factor( tp );
                // d+ grows the DRAWN scale by 1.2, three presses in a row
                press( cp, cp.zoomInDomainButtonForTest() );
                check( ok, layout + ": one d+ must grow the drawn scale by 1.2, got " + factor( tp ) + " from " + f0,
                       near( factor( tp ), f0 * 1.2f ) );
                press( cp, cp.zoomInDomainButtonForTest() );
                press( cp, cp.zoomInDomainButtonForTest() );
                check( ok, layout + ": three d+ presses must compound to 1.728, got " + factor( tp ) + " from " + f0,
                       near( factor( tp ), f0 * 1.728f ) );
                check( ok, layout + ": the radial width follows: " + tp.effectiveDomainStructureWidthForTest(),
                       Math.abs( tp.effectiveDomainStructureWidthForTest() - ( start * 1.728 ) ) < 0.05 );
                // d- shrinks it again
                press( cp, cp.zoomOutDomainButtonForTest() );
                check( ok, layout + ": d- must shrink the drawn scale by 0.8, got " + factor( tp ),
                       near( factor( tp ), f0 * 1.728f * 0.8f ) );
                // the rectangular width is untouched by radial presses...
                check( ok, layout + ": radial presses must leave the rectangular width alone, got "
                        + tp.domainStructureWidthForTest() + " vs " + rect_w, tp.domainStructureWidthForTest() == rect_w );
                // ...and a rectangular press leaves the radial width alone
                final double radial_w = tp.radialDomainStructureWidthForTest();
                mf[ 0 ].getOptions().setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                cp.displayedPhylogenyMightHaveChanged( true );
                press( cp, cp.zoomInDomainButtonForTest() );
                check( ok, layout + ": a rectangular d+ steps the rectangular width, got " + tp.domainStructureWidthForTest(),
                       Math.abs( tp.domainStructureWidthForTest() - ( rect_w * 1.2 ) ) < 0.05 );
                check( ok, layout + ": a rectangular d+ must leave the radial width alone, got "
                        + tp.radialDomainStructureWidthForTest() + " vs " + radial_w,
                       tp.radialDomainStructureWidthForTest() == radial_w );
                // back in the radial layout the radial width is still the stepped one, not re-derived from the cap
                mf[ 0 ].getOptions().setPhylogenyGraphicsType( layout );
                tp.setPhylogenyGraphicsType( layout );
                cp.displayedPhylogenyMightHaveChanged( true );
                check( ok, layout + ": returning to the layout keeps its stepped width, got "
                        + tp.effectiveDomainStructureWidthForTest() + " vs " + radial_w,
                       tp.effectiveDomainStructureWidthForTest() == radial_w );
                // the lower limit: shrinking stops once the width is at or under 20
                for ( int i = 0; i < 12; ++i ) {
                    press( cp, cp.zoomOutDomainButtonForTest() );
                }
                final double floor = tp.effectiveDomainStructureWidthForTest();
                press( cp, cp.zoomOutDomainButtonForTest() );
                check( ok, layout + ": at the floor (" + floor + ") a further d- is refused, got "
                        + tp.effectiveDomainStructureWidthForTest(),
                       ( floor <= 20.0 ) && ( tp.effectiveDomainStructureWidthForTest() == floor ) );
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    private static void press( final ControlPanel cp, final Object button ) {
        cp.actionPerformed( new ActionEvent( button, 0, "press" ) );
    }

    /** The drawn scale (px per residue) of the first architecture: what the width actually does to the tracks. */
    private static float factor( final TreePanel tp ) {
        for ( final PhylogenyNode n : tp.getPhylogeny().getExternalNodes() ) {
            if ( n.getNodeData().isHasSequence()
                    && ( n.getNodeData().getSequence().getDomainArchitecture() instanceof RenderableDomainArchitecture ) ) {
                return ( (RenderableDomainArchitecture) n.getNodeData().getSequence().getDomainArchitecture() )
                        .getRenderingFactorWidth();
            }
        }
        throw new IllegalStateException( "no renderable architecture" );
    }

    private static boolean near( final float a, final float b ) {
        return Math.abs( a - b ) < ( 0.002f + ( 0.002f * Math.abs( b ) ) );
    }

    private static void check( final boolean[] ok, final String what, final boolean cond ) {
        if ( !cond ) {
            System.out.println( "  [DomainRadialWidthTest] FAILED: " + what );
            ok[ 0 ] = false;
        }
    }

    private DomainRadialWidthTest() {
        // not instantiable
    }
}
