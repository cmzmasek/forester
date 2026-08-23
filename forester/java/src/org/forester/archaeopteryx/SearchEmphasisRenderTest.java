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

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Verifies the two WYSIWYG search-emphasis toggles on the demo tree (forester/demo/search-emphasis.xml): with a tip
 * marked as found, "Bold Found Labels" adds ink to the (red) found label, and "Dim Non-Matches" fades a NON-found
 * (black) label toward the white background. Headful; a green no-op when headless. Dogfoods the demo.
 */
public final class SearchEmphasisRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SearchEmphasisRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return emphasisRendersOk();
    }

    private static boolean emphasisRendersOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/search-emphasis.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "emph" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true ); // LIGHT theme -> found label is red, others black on white
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final PhylogenyNode found = tipNamed( phy, "AKT1_kinase" );
                    final PhylogenyNode other = tipNamed( phy, "actin" );
                    if ( ( found == null ) || ( other == null ) ) {
                        fail( ok, "demo must contain AKT1_kinase and actin" );
                        return;
                    }
                    final Set<Long> found_set = new HashSet<>();
                    found_set.add( found.getId() );
                    tp.setFoundNodes0( found_set );
                    final int w = 620, h = 560;
                    // BASELINE: both emphases off
                    final BufferedImage base = render( frame, tp, o, false, false, w, h );
                    final int red_base = redInk( base, found );
                    final int dark_base = darkInk( base, other );
                    if ( red_base < 40 ) {
                        fail( ok, "the found label should already be red (found color); got " + red_base + " red px" );
                    }
                    if ( dark_base < 40 ) {
                        fail( ok, "the non-found label should have dark ink at baseline; got " + dark_base + " px" );
                    }
                    // BOLD ON: the found (red) label gains ink
                    final BufferedImage bold = render( frame, tp, o, true, false, w, h );
                    final int red_bold = redInk( bold, found );
                    if ( red_bold <= ( red_base + 20 ) ) {
                        fail( ok, "Bold Found Labels should thicken the found label (red px bold=" + red_bold + " base="
                                + red_base + ")" );
                    }
                    // DIM ON: the non-found (black) label fades toward the white background (loses its near-black ink)
                    final BufferedImage dim = render( frame, tp, o, false, true, w, h );
                    final int dark_dim = darkInk( dim, other );
                    if ( dark_dim >= ( dark_base / 2 ) ) {
                        fail( ok, "Dim Non-Matches should fade the non-found label (dark px dim=" + dark_dim + " base="
                                + dark_base + ")" );
                    }
                    // the found label must NOT be dimmed: still fully red with dim on
                    if ( redInk( dim, found ) < ( red_base - 20 ) ) {
                        fail( ok, "a found label must NOT be dimmed (red px with dim on="
                                + redInk( dim, found ) + " base=" + red_base + ")" );
                    }
                    // the NUMBERS dim too: with branch-length + confidence values shown, a non-hit tip's branch-length
                    // number fades when Dim is on (it routes through the same dimNonMatch helper as the labels).
                    tp.getControlPanel().setCheckbox( DisplayOption.WRITE_BRANCH_LENGTH_VALUES, true );
                    tp.getControlPanel().setCheckbox( DisplayOption.WRITE_CONFIDENCE_VALUES, true );
                    final PhylogenyNode ins = tipNamed( phy, "insulin" ); // a non-hit tip
                    final int bl_off = branchLenInk( render( frame, tp, o, false, false, w, h ), ins );
                    final int bl_on = branchLenInk( render( frame, tp, o, false, true, w, h ), ins );
                    if ( bl_off < 12 ) {
                        fail( ok, "the non-hit branch-length number should have ink to begin with; got " + bl_off );
                    }
                    if ( bl_on >= ( bl_off / 2 ) ) {
                        fail( ok, "a non-hit branch-length number should fade when Dim is on (dark px off=" + bl_off
                                + " on=" + bl_on + ")" );
                    }
                    // and a non-hit internal node's CONFIDENCE value fades too (paintConfidenceValues path). "insulin"'s
                    // parent clade carries bootstrap 76; measure its number just below the branch (away from the line).
                    final PhylogenyNode inode = ins.getParent();
                    final int cf_off = confInk( render( frame, tp, o, false, false, w, h ), inode );
                    final int cf_on = confInk( render( frame, tp, o, false, true, w, h ), inode );
                    if ( cf_off < 8 ) {
                        fail( ok, "the non-hit confidence value should have ink to begin with; got " + cf_off );
                    }
                    if ( cf_on >= ( cf_off / 2 ) ) {
                        fail( ok, "a non-hit confidence value should fade when Dim is on (dark px off=" + cf_off
                                + " on=" + cf_on + ")" );
                    }
                    // Dim must NOT engage when NO hit is VISIBLE (the only hit hidden under a collapse) -- else the
                    // whole tree washes out with nothing emphasised. "insulin" (a different clade) stays visible.
                    final PhylogenyNode other2 = tipNamed( phy, "insulin" );
                    if ( other2 == null ) {
                        fail( ok, "demo must contain insulin" );
                        return;
                    }
                    final int insulin_dim = darkInk( render( frame, tp, o, false, true, w, h ), other2 ); // hit visible
                    tp.collapse( found.getParent() ); // hide the only found tip (AKT1_kinase) under a collapsed clade
                    final int insulin_free = darkInk( render( frame, tp, o, false, true, w, h ), other2 ); // no vis. hit
                    if ( insulin_free <= ( insulin_dim + 30 ) ) {
                        fail( ok, "dim must lift when the only hit is hidden (insulin dark px hidden-hit=" + insulin_free
                                + " visible-hit=" + insulin_dim + ")" );
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

    private static BufferedImage render( final MainFrame frame, final TreePanel tp, final Options o,
                                         final boolean bold, final boolean dim, final int w, final int h ) {
        o.setBoldFoundLabels( bold );
        o.setDimNonMatches( dim );
        frame.showWhole();
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        return AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ); // scale 1 -> coords are pixel coords
    }

    /** Red ink (r dominant) in the label band just right of the node -- the found label is drawn in the found color
     *  (red in the light theme). */
    private static int redInk( final BufferedImage img, final PhylogenyNode tip ) {
        return countBand( img, tip, true );
    }

    /** Near-black ink in the label band just right of the node -- a normal (non-found, non-dimmed) label. */
    private static int darkInk( final BufferedImage img, final PhylogenyNode tip ) {
        return countBand( img, tip, false );
    }

    /** Near-black ink in a non-hit tip's branch-length number box (drawn above the branch at parent_x+3,
     *  node_y-descent) -- interior of the tree, away from tip labels. */
    private static int branchLenInk( final BufferedImage img, final PhylogenyNode tip ) {
        final int px = Math.round( tip.getParent().getXcoord() );
        final int x0 = Math.max( 0, px + 3 ), x1 = Math.min( img.getWidth() - 1, px + 55 );
        final int y0 = Math.max( 0, Math.round( tip.getYcoord() ) - 15 );
        final int y1 = Math.min( img.getHeight() - 1, Math.round( tip.getYcoord() ) - 1 );
        int n = 0;
        for( int y = y0; y <= y1; ++y ) {
            for( int x = x0; x <= x1; ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( r <= 90 ) && ( g <= 90 ) && ( b <= 90 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Near-black ink in an internal node's confidence-value box: centred between parent_x and node_x, just BELOW
     *  the branch line (drawn at node_y + ascent) so the horizontal branch itself is excluded. */
    private static int confInk( final BufferedImage img, final PhylogenyNode inode ) {
        final int px = Math.round( inode.getParent().getXcoord() );
        final int nx = Math.round( inode.getXcoord() );
        final int x0 = Math.max( 0, px + 4 ), x1 = Math.min( img.getWidth() - 1, nx - 4 );
        final int y0 = Math.max( 0, Math.round( inode.getYcoord() ) + 1 );
        final int y1 = Math.min( img.getHeight() - 1, Math.round( inode.getYcoord() ) + 15 );
        int n = 0;
        for( int y = y0; y <= y1; ++y ) {
            for( int x = x0; x <= x1; ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( r <= 90 ) && ( g <= 90 ) && ( b <= 90 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static int countBand( final BufferedImage img, final PhylogenyNode tip, final boolean red ) {
        final int x0 = Math.max( 0, Math.round( tip.getXcoord() ) + 8 ); // skip the node dot
        final int x1 = Math.min( img.getWidth() - 1, Math.round( tip.getXcoord() ) + 190 );
        final int y0 = Math.max( 0, Math.round( tip.getYcoord() ) - 8 );
        final int y1 = Math.min( img.getHeight() - 1, Math.round( tip.getYcoord() ) + 8 );
        int n = 0;
        for( int y = y0; y <= y1; ++y ) {
            for( int x = x0; x <= x1; ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( red ) {
                    if ( ( r >= 150 ) && ( g <= 110 ) && ( b <= 110 ) ) {
                        ++n;
                    }
                }
                else if ( ( r <= 90 ) && ( g <= 90 ) && ( b <= 90 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static PhylogenyNode tipNamed( final Phylogeny phy, final String name ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        return null;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [SearchEmphasisRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [SearchEmphasisRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private SearchEmphasisRenderTest() {
    }
}
