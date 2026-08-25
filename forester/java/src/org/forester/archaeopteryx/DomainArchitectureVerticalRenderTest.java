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

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the domain-architecture demo (forester/demo/domain-architectures.xml) as a phylogram and asserts the
 * colored domain boxes appear in a vertical (root-top) orientation as well as the horizontal one -- i.e. each tip's
 * renderable domain architecture rides the canvas rotation into a vertical track. Headful; a green no-op when
 * headless. Dogfoods the demo.
 */
public final class DomainArchitectureVerticalRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DomainArchitectureVerticalRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return domainsRenderOk();
    }

    private static boolean domainsRenderOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/domain-architectures.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "dom" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true ); // predictable white bg so the colored boxes read
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    if ( !tp.getControlPanel().isShowDomainArchitectures() ) {
                        fail( ok, "the demo tree should auto-enable Domain Architectures" );
                    }
                    final int w = 760, h = 620;
                    // horizontal baseline: the domain boxes are drawn as colored (non-gray) pixels
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final int horiz = countColorful( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( horiz < 800 ) {
                        fail( ok, "the domain architectures should paint many colored pixels horizontally, got " + horiz );
                    }
                    // vertical: the boxes ride R into per-tip vertical tracks -- a comparable amount of colored ink
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage vimg = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int vert = countColorful( vimg );
                    if ( vert < ( horiz / 2 ) ) {
                        fail( ok, "the domain architectures should draw in a vertical orientation (vertical=" + vert
                                + " horizontal=" + horiz + ")" );
                    }
                    // ALIGNMENT: on this UNALIGNED phylogram the tips sit at RAGGED depths, but the domain tracks must
                    // still line up in a COMMON column (like the horizontal layout) -- so every track's TOP is at the
                    // same depth and the per-column top-of-domain spread is small. A per-tip (ragged) layout spreads
                    // the tops across the tree height. Measured over the left 60% (clear of the bottom-right legend).
                    final int spread = domainTrackTopSpread( vimg, w, h );
                    if ( spread > 50 ) {
                        fail( ok, "vertical domain tracks must line up in a common column, not hang ragged off each tip"
                                + " (top-of-domain spread=" + spread + "px)" );
                    }
                    // and root-at-bottom too
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_BOTTOM );
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final int vert_b = countColorful( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( vert_b < ( horiz / 2 ) ) {
                        fail( ok, "the domain architectures should draw in a root-at-bottom orientation (got " + vert_b
                                + " horizontal=" + horiz + ")" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Count "colored" (clearly non-gray, not near-white/black) pixels -- the pastel domain boxes. Branches (black),
     *  labels (dark gray text) and the white background are excluded, so on this tree the count is the domain ink. */
    private static int countColorful( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int max = Math.max( r, Math.max( g, b ) ), min = Math.min( r, Math.min( g, b ) );
                if ( ( ( max - min ) > 40 ) && ( max > 90 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Spread (max - min) of the per-column TOP of domain ink in a root-top render, over the left 60% of the width
     *  (clear of the bottom-right legend). Small when the tracks line up in a common column; large when each hangs
     *  ragged off its own ragged-depth tip. Only the domain boxes are "colorful"; branches and text are not. */
    private static int domainTrackTopSpread( final BufferedImage img, final int w, final int h ) {
        final int x_end = ( w * 60 ) / 100;
        int min_top = Integer.MAX_VALUE, max_top = Integer.MIN_VALUE;
        for( int x = 0; x < x_end; ++x ) {
            for( int y = 0; y < h; ++y ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int max = Math.max( r, Math.max( g, b ) ), min = Math.min( r, Math.min( g, b ) );
                if ( ( ( max - min ) > 40 ) && ( max > 90 ) ) {
                    min_top = Math.min( min_top, y ); // first (top) colorful pixel in this column
                    max_top = Math.max( max_top, y );
                    break;
                }
            }
        }
        return ( max_top < min_top ) ? 0 : ( max_top - min_top );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [DomainArchitectureVerticalRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [DomainArchitectureVerticalRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private DomainArchitectureVerticalRenderTest() {
    }
}
