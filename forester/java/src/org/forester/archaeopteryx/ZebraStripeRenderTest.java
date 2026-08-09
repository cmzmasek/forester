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
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the wide demo tree (forester/demo/zebra-stripes.xml) with "Zebra Stripes" ON and asserts faint alternating
 * row bands appear over a white background (and none when off), and that the bands alternate leaf row by leaf row.
 * Headful; a green no-op when headless. Dogfoods the demo.
 */
public final class ZebraStripeRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ZebraStripeRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return stripesRenderOk();
    }

    private static boolean stripesRenderOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/zebra-stripes.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "zebra" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true ); // predictable white background so the faint tint reads
                    // blank the tip labels so the ONLY faint-gray pixels come from the stripes, not antialiased text
                    // edges -- makes the pixel thresholds robust across fonts / rasterizers / HiDPI
                    for( final PhylogenyNode tip : phy.getExternalNodes() ) {
                        tip.setName( "" );
                    }
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final int w = 900, h = 620;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    o.setShowZebraStripes( false );
                    final int off = countFaintGray( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    o.setShowZebraStripes( true );
                    final BufferedImage on_img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int on = countFaintGray( on_img );
                    if ( on <= ( off + 3000 ) ) {
                        fail( ok, "Zebra Stripes should add many faint-gray pixels (on=" + on + " off=" + off + ")" );
                    }
                    if ( off >= 1500 ) {
                        fail( ok, "no faint bands should appear when Zebra Stripes is off, got " + off );
                    }
                    // alternation: external rows in y-order must alternate striped / not-striped
                    final List<PhylogenyNode> tips = new ArrayList<>( phy.getExternalNodes() );
                    tips.sort( Comparator.comparingDouble( PhylogenyNode::getYcoord ) );
                    Boolean prev = null;
                    int striped = 0;
                    for( final PhylogenyNode tip : tips ) {
                        final boolean s = rowIsStriped( on_img, Math.round( tip.getYcoord() ) );
                        if ( s ) {
                            ++striped;
                        }
                        if ( ( prev != null ) && ( prev.booleanValue() == s ) ) {
                            fail( ok, "adjacent leaf rows must alternate striped/not-striped (row y=" + Math
                                    .round( tip.getYcoord() ) + ")" );
                            break;
                        }
                        prev = Boolean.valueOf( s );
                    }
                    if ( ( striped < ( ( tips.size() / 2 ) - 1 ) ) || ( striped > ( ( tips.size() / 2 ) + 1 ) ) ) {
                        fail( ok, "about half the leaf rows should be striped, got " + striped + " of " + tips.size() );
                    }
                    // VERTICAL PARITY: the row bands ride R into vertical column bands over alternate tips in a
                    // root-top/bottom orientation (the cross-tree extent becomes the depth/height axis). Render in
                    // ROOT_TOP and confirm the faint bands still draw.
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    final int vertical_on = countFaintGray( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1,
                            false ) );
                    if ( vertical_on <= ( off + 3000 ) ) {
                        fail( ok, "Zebra Stripes should draw in a vertical orientation (on=" + vertical_on + " off="
                                + off + ")" );
                    }
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    // a transparent-PNG export must NOT bake in the full-width bands (they'd defeat the clean
                    // cut-out) -- turning zebra on must add no semi-transparent stripe pixels vs. zebra off
                    o.setTransparentExportBackground( true );
                    o.setShowZebraStripes( false );
                    final int t_off = countSemiTransparent( AptxUtil.renderPhylogenyToImage( w, h, tp, o, true, 1,
                            false ) );
                    o.setShowZebraStripes( true );
                    final int t_on = countSemiTransparent( AptxUtil.renderPhylogenyToImage( w, h, tp, o, true, 1,
                            false ) );
                    o.setTransparentExportBackground( false );
                    if ( t_on > ( t_off + 200 ) ) {
                        fail( ok, "zebra stripes must be suppressed on a transparent-PNG export (on=" + t_on + " off="
                                + t_off + ")" );
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

    /** Count semi-transparent pixels (0 < alpha < 100) in an ARGB export -- the alpha-16 stripe bands over the
     *  otherwise-transparent canvas fall in this range, while the fully-opaque branches and fully-transparent
     *  background do not; antialiased branch edges are the same with zebra on or off, so they cancel in the
     *  on-vs-off comparison. */
    private static int countSemiTransparent( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int a = ( img.getRGB( x, y ) >>> 24 ) & 0xFF;
                if ( ( a > 0 ) && ( a < 100 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** A striped row's full width is tinted faint gray (the white gaps become ~239), so its row band carries far more
     *  faint-gray pixels than an un-striped row (which only has a little antialiased text). */
    private static boolean rowIsStriped( final BufferedImage img, final int y ) {
        int n = 0;
        for( int yy = Math.max( 0, y - 1 ); yy < Math.min( img.getHeight(), y + 2 ); ++yy ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( isFaintGray( img.getRGB( x, yy ) ) ) {
                    ++n;
                }
            }
        }
        return n > ( img.getWidth() / 4 );
    }

    private static int countFaintGray( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( isFaintGray( img.getRGB( x, y ) ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** A faint neutral gray -- the translucent stripe over the white background (~239), NOT pure white, black
     *  branches, or the darker/colored antialiased text. */
    private static boolean isFaintGray( final int rgb ) {
        final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
        return ( r >= 232 ) && ( r <= 250 ) && ( Math.abs( r - g ) <= 5 ) && ( Math.abs( g - b ) <= 5 );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ZebraStripeRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ZebraStripeRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private ZebraStripeRenderTest() {
    }
}
