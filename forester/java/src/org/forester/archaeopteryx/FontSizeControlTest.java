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

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;

/**
 * Integration test for the control-panel font-size slider that replaced the "Font Size" menu, and for the
 * single-source-of-truth {@link TreeFontSet} user-size API. Builds a real {@link MainFrame} (needs FlatLaf
 * + a toolkit) and exercises both directions: the {@code TreeFontSet} size methods, and that the slider
 * reflects the size (sync) and applies it (apply). A no-op when headless, like the other GUI tests here.
 */
public final class FontSizeControlTest {

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a toolkit + FlatLaf
        }
        try {
            final Phylogeny phy = Phylogeny.createInstanceFromNhxString( "((Apis,Bombus)x,(Felis,Canis)y)root" );
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "font size test" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreeFontSet tfs = mf[ 0 ].getMainPanel().getTreeFontSet();
                final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();

                // the default is the tuned base size; small font is large minus 2
                if ( tfs.getUserFontSize() != AptxConstants.DEFAULT_TREE_FONT_SIZE ) {
                    ok[ 0 ] = false;
                    System.out.println( "  default user font size wrong: " + tfs.getUserFontSize() );
                }
                // setUserFontSize re-derives large/small coherently (the single source of truth)
                tfs.setUserFontSize( 18 );
                if ( ( tfs.getUserFontSize() != 18 ) || ( tfs.getLargeFont().getSize() != 18 )
                        || ( tfs.getSmallFont().getSize() != 16 ) ) {
                    ok[ 0 ] = false;
                    System.out.println( "  setUserFontSize(18) -> large/small wrong: " + tfs.getLargeFont().getSize()
                            + "/" + tfs.getSmallFont().getSize() );
                }
                // increase/decrease step by 1
                tfs.increaseUserFontSize();
                tfs.decreaseUserFontSize();
                if ( tfs.getUserFontSize() != 18 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  increase/decrease did not round-trip: " + tfs.getUserFontSize() );
                }
                // clamping to the slider bounds
                tfs.setUserFontSize( 10000 );
                if ( tfs.getUserFontSize() != TreeFontSet.MAX_FONT_SIZE ) {
                    ok[ 0 ] = false;
                    System.out.println( "  upper clamp failed: " + tfs.getUserFontSize() );
                }
                tfs.setUserFontSize( -9 );
                if ( tfs.getUserFontSize() != TreeFontSet.MIN_FONT_SIZE ) {
                    ok[ 0 ] = false;
                    System.out.println( "  lower clamp failed: " + tfs.getUserFontSize() );
                }
                // sync: the slider tracks the user size
                tfs.setUserFontSize( 14 );
                cp.updateFontSizeSlider();
                if ( cp.getFontSizeSliderValue() != 14 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  slider did not sync to the user size: " + cp.getFontSizeSliderValue() );
                }
                // apply: moving the slider sets the user size (the change-listener path, incl. the re-entry guard)
                cp.setFontSizeSliderValue( 9 );
                if ( tfs.getUserFontSize() != 9 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  slider did not apply the size: " + tfs.getUserFontSize() );
                }
                // ---- the NODE-SIZE slider, directly below: same shape, on Options.defaultNodeShapeSize ----
                mf[ 0 ].getOptions().setDefaultNodeShapeSize( (short) 7 );
                cp.updateNodeSizeSlider();
                if ( cp.getNodeSizeSliderValue() != 7 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  node slider did not sync from Options: " + cp.getNodeSizeSliderValue() );
                }
                cp.setNodeSizeSliderValue( 12 );
                if ( mf[ 0 ].getOptions().getDefaultNodeShapeSize() != 12 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  node slider did not apply the size: "
                            + mf[ 0 ].getOptions().getDefaultNodeShapeSize() );
                }
                // A value set from Settings (range 0..100) beyond the slider's everyday range must be SHOWN, not
                // silently clamped -- a clamped slider would misreport the size and change it on the next touch.
                mf[ 0 ].getOptions().setDefaultNodeShapeSize( (short) 40 );
                cp.updateNodeSizeSlider();
                if ( cp.getNodeSizeSliderValue() != 40 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  node slider must stretch to show an out-of-range value, got "
                            + cp.getNodeSizeSliderValue() );
                }
                if ( mf[ 0 ].getOptions().getDefaultNodeShapeSize() != 40 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  ...without changing it: " + mf[ 0 ].getOptions().getDefaultNodeShapeSize() );
                }
                ( (javax.swing.JFrame) mf[ 0 ] ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private FontSizeControlTest() {
    }
}
