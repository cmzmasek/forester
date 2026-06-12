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
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Integration test for the long-standing "search highlight bleeds into PDF branches" bug. iText's
 * {@code PdfGraphics2D} draws a bold label (a search-found node's label is bold) by *stroking* the
 * glyphs, setting the PDF stroke color directly while bypassing its own stroke-color cache; a
 * following stroke whose color matches that now-stale cache is then skipped by iText and silently
 * inherits the label's red/green color. So in a regular phylogram the *branches* of the next clade
 * came out red/green, and in an aligned phylogram the lined-up *connector* did -- in the PDF only;
 * raster exports (PNG/TIFF) and the screen were always fine. {@code resyncPdfStrokeColor} re-syncs
 * the cache after each found label, so nothing past the found node's own box+label is colored.
 *
 * <p>This drives the real iText PDF path ({@link PdfExporter}) -- for both a regular and an aligned
 * phylogram -- and rasterizes the result with {@code sips} (macOS), then asserts that color appears
 * only on the found nodes' own rows. It needs FlatLaf + a display + {@code sips}, so it is a
 * standalone test, not part of the headless suite, and it is a no-op (returns true) wherever any of
 * those is unavailable.
 */
public final class FoundBranchExportTest {

    private static final int FOUND_COUNT = 2; // one node in each highlight set (green + red)

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "FoundBranchExport: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI + PDF integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration( null, false, false, true );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, conf, "found" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    ok[ 0 ] = run( mf[ 0 ] );
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
                }
                ( (JFrame) mf[ 0 ] ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean run( final MainFrame mf ) throws Exception {
        final MainPanel mp = mf.getMainPanel();
        final TreePanel tp = mp.getCurrentTreePanel();
        tp.setSize( 900, 700 );
        tp.validate();
        mp.getOptions().setLineUpRendarableNodeData( true ); // draws aligned connectors when aligned

        // mark two leaves as search-found (one in each highlight set: green and red)
        final Set<Long> found0 = new HashSet<>();
        final Set<Long> found1 = new HashSet<>();
        int leaf = 0;
        for ( final PhylogenyNodeIterator it = tp.getPhylogeny().iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isExternal() ) {
                leaf++;
                if ( leaf == 2 ) {
                    found0.add( n.getId() );
                }
                if ( leaf == 5 ) {
                    found1.add( n.getId() );
                }
            }
        }
        tp.setFoundNodes0( found0 );
        tp.setFoundNodes1( found1 );

        // a regular phylogram (next clade's branches used to leak) and an aligned one (connector leaked)
        return checkExport( tp, Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM, "regular phylogram" )
                && checkExport( tp, Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM, "aligned phylogram" );
    }

    private static boolean checkExport( final TreePanel tp, final Options.PHYLOGENY_DISPLAY_TYPE type,
                                        final String label ) throws Exception {
        tp.getControlPanel().setTreeDisplayType( type );
        tp.printAll( new BufferedImage( 900, 700, BufferedImage.TYPE_INT_RGB ).createGraphics() ); // lay out coords

        final File pdf = File.createTempFile( "aptx_found_export", ".pdf" );
        final File png = File.createTempFile( "aptx_found_export", ".png" );
        try {
            PdfExporter.writePhylogenyToPdf( pdf.getAbsolutePath(), tp, tp.getWidth(), tp.getHeight() );
            if ( !rasterize( pdf, png ) ) {
                System.out.println( "  [FoundBranchExportTest] skipped: could not rasterize PDF (sips unavailable?)" );
                return true;
            }
            final BufferedImage img = ImageIO.read( png );
            if ( img == null ) {
                System.out.println( "  [FoundBranchExportTest] skipped: rasterized image unreadable" );
                return true;
            }
            final int bands = countColoredRowBands( img );
            if ( bands > FOUND_COUNT ) {
                System.out.println( "  [FoundBranchExportTest] the " + label + " export leaks the search-found color "
                        + "past the found nodes: " + bands + " colored row-bands but only " + FOUND_COUNT + " found "
                        + "nodes (branches or connectors are being colored)" );
                return false;
            }
            return true;
        }
        finally {
            pdf.delete();
            png.delete();
        }
    }

    private static boolean rasterize( final File pdf, final File png ) {
        try {
            final Process p = new ProcessBuilder( "sips", "-s", "format", "png", "--resampleHeight", "700",
                                                  pdf.getAbsolutePath(), "--out", png.getAbsolutePath() )
                    .redirectErrorStream( true )
                    .redirectOutput( ProcessBuilder.Redirect.DISCARD )
                    .start();
            if ( !p.waitFor( 30, TimeUnit.SECONDS ) ) {
                p.destroyForcibly();
                return false;
            }
            return ( p.exitValue() == 0 ) && png.isFile() && ( png.length() > 0 );
        }
        catch ( final Exception e ) {
            return false; // sips not present / not this platform -> skip
        }
    }

    // Number of distinct horizontal bands (groups of consecutive rows) that contain fully-saturated
    // found colors (red 255,0,0 or green 0,255,0). After the fix this equals the number of found
    // nodes (only their own box+label rows are colored); a leaked connector or branch adds bands.
    private static int countColoredRowBands( final BufferedImage img ) {
        final int w = img.getWidth();
        final int h = img.getHeight();
        int bands = 0;
        int gap = 99;
        for ( int y = 0; y < h; y++ ) {
            int colored = 0;
            for ( int x = 0; x < w; x++ ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xff;
                final int g = ( rgb >> 8 ) & 0xff;
                final int b = rgb & 0xff;
                if ( ( ( g > 150 ) && ( r < 90 ) && ( b < 90 ) ) || ( ( r > 150 ) && ( g < 90 ) && ( b < 90 ) ) ) {
                    colored++;
                }
            }
            if ( colored >= 3 ) {
                if ( gap >= 2 ) {
                    bands++; // start of a new band (separated from the previous by >= 2 blank rows)
                }
                gap = 0;
            }
            else {
                gap++;
            }
        }
        return bands;
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        for ( int i = 0; i < 3; i++ ) {
            final PhylogenyNode inter = new PhylogenyNode();
            inter.setDistanceToParent( 0.5 );
            for ( int j = 0; j < 3; j++ ) {
                final PhylogenyNode leaf = new PhylogenyNode();
                leaf.setName( "leaf_" + i + "_" + j );
                leaf.setDistanceToParent( 0.3 + ( 0.1 * j ) );
                inter.addAsChild( leaf );
            }
            root.addAsChild( inter );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private FoundBranchExportTest() {
    }
}
