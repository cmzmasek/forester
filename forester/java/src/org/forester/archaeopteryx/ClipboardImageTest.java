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
import java.awt.Image;
import java.awt.Transparency;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.image.BufferedImage;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Integration test for the "Copy Image to Clipboard" plumbing on a real {@link MainFrame}/{@link TreePanel}:
 * {@link AptxUtil#renderPhylogenyToImage} must return an opaque, actually-painted raster, and
 * {@link AptxUtil#copyPhylogenyImageToClipboard} must place an {@code imageFlavor} transferable on the clipboard.
 * A PRIVATE {@link Clipboard} is used so the test never clobbers the user's real system clipboard. Guarded to a
 * no-op when headless (needs FlatLaf via {@code createInstance}).
 */
public final class ClipboardImageTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ClipboardImage: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, conf, "clip" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try { // dispose the frame even if an assertion path throws (symmetry with the NO_TREE block below)
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                final Options opts = frame.getOptions();
                opts.setGraphicsExportVisibleOnly( false ); // deterministic: render the whole 800x600 canvas
                opts.setRasterExportScale( 1 );
                frame.showWhole();
                final int w = 800, h = 600;
                tp.setSize( w, h );
                tp.calcParametersForPainting( w, h );

                // renderPhylogenyToImage (opaque, scale 1, whole figure): a right-sized, actually-painted image
                final BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, opts, false, 1, false );
                if ( img == null ) {
                    fail( ok, "renderPhylogenyToImage returned null for a real tree" );
                }
                else {
                    if ( ( img.getWidth() != w ) || ( img.getHeight() != h ) ) {
                        fail( ok, "image size " + img.getWidth() + "x" + img.getHeight() + " != requested " + w + "x"
                                + h );
                    }
                    if ( img.getTransparency() != Transparency.OPAQUE ) {
                        fail( ok, "a clipboard/opaque image must have no alpha channel" );
                    }
                    if ( distinctColorCount( img ) < 2 ) {
                        fail( ok, "the rendered image is a single flat color -- the tree was not painted" );
                    }
                }

                // copyPhylogenyImageToClipboard: writes an image transferable onto a PRIVATE clipboard
                final Clipboard clip = new Clipboard( "aptx-test" );
                final boolean copied = AptxUtil.copyPhylogenyImageToClipboard( tp, opts, clip, null );
                if ( !copied ) {
                    fail( ok, "copyPhylogenyImageToClipboard returned false for a real tree" );
                }
                else {
                    try {
                        final Transferable t = clip.getContents( null );
                        if ( ( t == null ) || !t.isDataFlavorSupported( DataFlavor.imageFlavor ) ) {
                            fail( ok, "the clipboard does not hold an image after the copy" );
                        }
                        else if ( !( t.getTransferData( DataFlavor.imageFlavor ) instanceof Image ) ) {
                            fail( ok, "the clipboard image flavor did not yield an Image" );
                        }
                    }
                    catch ( final Exception e ) {
                        fail( ok, "reading the clipboard image back failed: " + e.getMessage() );
                    }
                }

                // the clipboard render now honors the user's raster-export scale (publication quality), so a high
                // export scale must produce a correspondingly larger bitmap -- NOT a screen-scale one. (cw x ch here
                // is small, so cw*4 x ch*4 stays well under renderPhylogenyToImage's 100 MP cap; scale is not clamped.)
                opts.setRasterExportScale( 4 );
                final int cw = tp.getWidth();
                final int ch = tp.getHeight();
                final Clipboard scale_clip = new Clipboard( "aptx-scale-test" );
                if ( !AptxUtil.copyPhylogenyImageToClipboard( tp, opts, scale_clip, null ) ) {
                    fail( ok, "copyPhylogenyImageToClipboard returned false with export scale 4" );
                }
                else {
                    try {
                        final BufferedImage cbi = (BufferedImage) scale_clip.getContents( null )
                                .getTransferData( DataFlavor.imageFlavor );
                        if ( ( cbi.getWidth() != cw * 4 ) || ( cbi.getHeight() != ch * 4 ) ) {
                            fail( ok, "clipboard image must honor the 4x export scale (" + ( cw * 4 ) + "x" + ( ch * 4 )
                                    + "), got " + cbi.getWidth() + "x" + cbi.getHeight() );
                        }
                    }
                    catch ( final Exception e ) {
                        fail( ok, "reading the scale-test clipboard image failed: " + e.getMessage() );
                    }
                }

                // MainFrame.copyImageToClipboard(Clipboard) core: COPIED for a real tree, image lands on the clipboard
                final Clipboard mf_clip = new Clipboard( "aptx-mf-test" );
                final MainFrame.ClipboardCopyResult res = frame.copyImageToClipboard( mf_clip );
                if ( res != MainFrame.ClipboardCopyResult.COPIED ) {
                    fail( ok, "copyImageToClipboard(core) should return COPIED for a real tree, got " + res );
                }
                else if ( ( mf_clip.getContents( null ) == null )
                        || !mf_clip.getContents( null ).isDataFlavorSupported( DataFlavor.imageFlavor ) ) {
                    fail( ok, "copyImageToClipboard(core) did not put an image on the clipboard" );
                }

                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            // NO_TREE: a frame with an empty/absent current tree must report NO_TREE and write nothing
            SwingUtilities.invokeAndWait( () -> {
                MainFrame empty = null;
                try {
                    empty = MainFrameApplication.createInstance( new Phylogeny[] { new Phylogeny() }, conf, "empty" );
                    final Clipboard c = new Clipboard( "aptx-empty-test" );
                    if ( empty.copyImageToClipboard( c ) != MainFrame.ClipboardCopyResult.NO_TREE ) {
                        fail( ok, "an empty/absent current tree must yield NO_TREE" );
                    }
                    if ( c.getContents( null ) != null ) {
                        fail( ok, "NO_TREE must not write anything to the clipboard" );
                    }
                }
                catch ( final Throwable t ) {
                    // if this toolkit refuses an empty-tree frame, skip rather than fail the suite
                }
                finally {
                    if ( empty instanceof JFrame ) {
                        ( (JFrame) empty ).dispose();
                    }
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ClipboardImageTest] " + msg );
        ok[ 0 ] = false;
    }

    /** Distinct RGB values over a coarse sample of the image (enough to tell "blank" from "something drawn"). */
    private static int distinctColorCount( final BufferedImage img ) {
        final Set<Integer> colors = new HashSet<>();
        for( int x = 0; x < img.getWidth(); x += 4 ) {
            for( int y = 0; y < img.getHeight(); y += 4 ) {
                colors.add( img.getRGB( x, y ) & 0x00FFFFFF );
                if ( colors.size() >= 2 ) {
                    return colors.size();
                }
            }
        }
        return colors.size();
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 6; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "tip" + i );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private ClipboardImageTest() {
    }
}
