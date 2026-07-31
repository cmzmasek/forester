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

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Integration test for the "White Image Background" export option. With a DARK tree color scheme active, an
 * opaque raster render (the path shared by file export and Copy Image to Clipboard) must, when the option is on,
 * come out DOCUMENT-READY: a white background AND dark, legible ink (not just white-behind-light-ink, which would
 * be an invisible white-on-white figure). When off it must be WYSIWYG (the dark theme background). Also checks
 * the Settings checkbox drives the option through updateOptions. No-op when headless (needs FlatLaf).
 */
public final class ExportBackgroundTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ExportBackground: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, conf, "bg" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                final Options opts = frame.getOptions();
                final int w = 400, h = 300;
                frame.showWhole();
                tp.setSize( w, h );

                // force the DARK scheme (branches/text are light) so the two cases are clearly distinguishable
                tp.getTreeColorSet().setColorSchema( TreeColorSet.DARK_COLOR_SCHEME );
                final int theme_bg = tp.getTreeColorSet().getBackgroundColor().getRGB() & 0xFFFFFF;
                if ( theme_bg == 0xFFFFFF ) {
                    fail( ok, "test setup: dark scheme background should not be white" );
                }

                // option ON -> white background AND dark, VISIBLE ink (document-ready), despite the dark theme
                opts.setGraphicsExportWhiteBackground( true );
                tp.getTreeColorSet().setColorSchema( TreeColorSet.DARK_COLOR_SCHEME );
                tp.calcParametersForPainting( w, h );
                final BufferedImage on = AptxUtil.renderPhylogenyToImage( w, h, tp, opts, false, 1, false );
                if ( corner( on ) != 0xFFFFFF ) {
                    fail( ok, "white-background ON must give a white background corner, got " + hex( corner( on ) ) );
                }
                if ( darkPixelCount( on ) < 15 ) {
                    fail( ok, "white-background ON must keep the ink DARK/legible on white (light-on-white would be "
                            + "invisible); dark pixels=" + darkPixelCount( on ) );
                }
                // the transient light theme must be RESTORED after the render (no lasting side effect on the panel)
                if ( tp.getTreeColorSet().getCurrentColorScheme() != TreeColorSet.DARK_COLOR_SCHEME ) {
                    fail( ok, "the export must restore the on-screen color scheme afterwards" );
                }

                // option OFF -> the (dark) theme background, NOT white (WYSIWYG)
                opts.setGraphicsExportWhiteBackground( false );
                tp.getTreeColorSet().setColorSchema( TreeColorSet.DARK_COLOR_SCHEME );
                tp.calcParametersForPainting( w, h );
                final int corner_off = corner( AptxUtil.renderPhylogenyToImage( w, h, tp, opts, false, 1, false ) );
                if ( corner_off == 0xFFFFFF ) {
                    fail( ok, "white-background OFF in a dark theme must NOT give a white corner (WYSIWYG expected)" );
                }
                if ( corner_off != theme_bg ) {
                    fail( ok, "white-background OFF must use the theme background " + hex( theme_bg ) + ", got "
                            + hex( corner_off ) );
                }

                // the Settings checkbox must drive the option through actionPerformed -> updateOptions
                final javax.swing.JCheckBoxMenuItem cb = frame._graphics_export_white_background_cbmi;
                if ( cb == null ) {
                    fail( ok, "the White Image Background checkbox was not created" );
                }
                else {
                    final boolean cb_before = cb.isSelected();
                    cb.doClick();
                    if ( cb.isSelected() == cb_before ) {
                        fail( ok, "doClick should toggle the checkbox" );
                    }
                    if ( frame.getOptions().isGraphicsExportWhiteBackground() != cb.isSelected() ) {
                        fail( ok, "updateOptions must sync the option to the checkbox state" );
                    }
                }

                ( (JFrame) frame ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** A background pixel near the top-left corner (the whole canvas is bg-filled before the tree is drawn). */
    private static int corner( final BufferedImage img ) {
        return ( img == null ) ? -1 : ( img.getRGB( 2, 2 ) & 0xFFFFFF );
    }

    /** Count of clearly-dark pixels (ink) -- proves the tree content is visible, not white-on-white. */
    private static int darkPixelCount( final BufferedImage img ) {
        if ( img == null ) {
            return 0;
        }
        int n = 0;
        for( int x = 0; x < img.getWidth(); ++x ) {
            for( int y = 0; y < img.getHeight(); ++y ) {
                final int rgb = img.getRGB( x, y );
                if ( ( ( ( rgb >> 16 ) & 0xFF ) + ( ( rgb >> 8 ) & 0xFF ) + ( rgb & 0xFF ) ) < 150 ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static String hex( final int rgb ) {
        return String.format( "#%06X", rgb );
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ExportBackgroundTest] " + msg );
        ok[ 0 ] = false;
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 5; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private ExportBackgroundTest() {
    }
}
