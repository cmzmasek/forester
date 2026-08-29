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

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;

/**
 * Integration test for the SEQUENCE-ALIGNMENT track: on a real {@link MainFrame}/{@link TreePanel} it renders a tree
 * whose tips carry aligned molecular sequences and checks the alignment draws colored residue cells right of the labels
 * in the rectangular root-left layout, draws NONE in the vertical/circular layouts (approved parity exception), reserves
 * horizontal width so the pane scrolls, and embeds the cells in a VECTOR export. No-op on a headless box.
 */
public final class MsaTrackRenderTest {

    private static final String[] ALIGN = { "MKQLEDPFGH-WYVAST", "MKQIEDPFGY-WYVAST", "LRQMEDANGH-WFVCST",
            "MK-LEDPFGH-WYVAGT" };

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MsaTrackRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = alignedTree();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "msa" ) );
            final boolean[] ok = { true };
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            SwingUtilities.invokeAndWait( () -> {
                tp.getOptions().setGraphicsExportWhiteBackground( true );
                tp.getOptions().setShowOverview( false );
                tp.setOvOn( false );
                tp.getOptions().setMsaColumnWidth( 10 );
            } );

            // (1) alignment colors appear right of the labels in root-left, and NOT in vertical / circular
            checkLayout( mf[ 0 ], tp, ok, "root-left", Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR,
                    Options.TREE_ORIENTATION.ROOT_LEFT, 760, 460, true );
            checkLayout( mf[ 0 ], tp, ok, "root-top", Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR,
                    Options.TREE_ORIENTATION.ROOT_TOP, 600, 600, false );
            checkLayout( mf[ 0 ], tp, ok, "circular", Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                    Options.TREE_ORIENTATION.ROOT_LEFT, 600, 600, false );

            // (2) showing the alignment reserves horizontal width (so the JScrollPane can scroll to it)
            SwingUtilities.invokeAndWait( () -> {
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.setSize( 760, 460 );
                tp.getOptions().setShowMsa( false );
                tp.resetPreferredSize();
                final int off = tp.getPreferredSize().width;
                tp.getOptions().setShowMsa( true );
                tp.resetPreferredSize();
                final int on = tp.getPreferredSize().width;
                if ( on <= off ) {
                    fail( ok, "the alignment must widen the reserved/preferred width (" + off + " -> " + on + ")" );
                }
                // the column ruler reserves a bottom band in root-left, none when off
                if ( tp.msaRulerReserveForTest() <= 0 ) {
                    fail( ok, "the MSA column ruler must reserve a bottom band when shown in root-left" );
                }
                tp.getOptions().setShowMsa( false );
                if ( tp.msaRulerReserveForTest() != 0 ) {
                    fail( ok, "the MSA column ruler must reserve nothing when the alignment is hidden" );
                }
                tp.getOptions().setShowMsa( true );
            } );
            // nice tick step: ~5-8 ticks across the visible span (1/2/5 x 10^k)
            if ( ( TreePanel.niceColumnStepForTest( 30 ) != 5 ) || ( TreePanel.niceColumnStepForTest( 300 ) != 50 )
                    || ( TreePanel.niceColumnStepForTest( 6 ) != 1 ) ) {
                fail( ok, "niceColumnStep should give ~5-8 ticks (30->5, 300->50, 6->1), got "
                        + TreePanel.niceColumnStepForTest( 30 ) + "/" + TreePanel.niceColumnStepForTest( 300 ) + "/"
                        + TreePanel.niceColumnStepForTest( 6 ) );
            }

            // (2b) an alignment wider than its window is SCROLLABLE, and panning the column offset changes what renders
            SwingUtilities.invokeAndWait( () -> {
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.getOptions().setShowMsa( true );
                tp.getOptions().setMsaColumnWidth( 40 ); // wide cells -> the 34-col alignment exceeds its window
                tp.setPreferredSize( new java.awt.Dimension( 560, 460 ) );
                tp.setSize( 560, 460 );
                mf[ 0 ].showWhole();
                tp.calcParametersForPainting( 560, 460 );
                if ( !tp.isMsaScrollableForTest() ) {
                    fail( ok, "a wide-cell alignment must be scrollable (window narrower than the alignment)" );
                }
                if ( tp.getMsaVisibleColumnsForTest() >= ALIGN[ 0 ].length() ) {
                    fail( ok, "the window must show fewer columns than the alignment total" );
                }
                tp.setMsaColumnOffset( 0 );
                final BufferedImage a = AptxUtil.renderPhylogenyToImage( 560, 460, tp, tp.getOptions(), false, 1, false );
                tp.setMsaColumnOffset( 18 );
                final BufferedImage b = AptxUtil.renderPhylogenyToImage( 560, 460, tp, tp.getOptions(), false, 1, false );
                if ( !imagesDiffer( a, b ) ) {
                    fail( ok, "panning the alignment offset must change which columns render" );
                }
                tp.getOptions().setMsaColumnWidth( 10 );
                tp.setMsaColumnOffset( 0 );
            } );

            // (3) a VECTOR export embeds the residue cells (<rect> elements)
            SwingUtilities.invokeAndWait( () -> {
                try {
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    tp.getOptions().setShowMsa( true );
                    tp.setPreferredSize( new java.awt.Dimension( 900, 460 ) );
                    tp.setSize( 900, 460 );
                    mf[ 0 ].showWhole();
                    tp.calcParametersForPainting( 900, 460 );
                    final File dir = java.nio.file.Files.createTempDirectory( "aptx-msa" ).toFile();
                    dir.deleteOnExit();
                    final File svg = new File( dir, "out.svg" );
                    VectorGraphicsExporter.writePhylogenyToVectorGraphicsFile( svg.getAbsolutePath(), tp, 900, 460,
                            VectorGraphicsExporter.Format.SVG, true, true );
                    final String content = new String( java.nio.file.Files.readAllBytes( svg.toPath() ), "UTF-8" );
                    final int rects = content.split( "<rect", -1 ).length - 1;
                    if ( rects < 20 ) {
                        fail( ok, "the SVG export must embed the alignment cells (>=20 <rect>), got " + rects );
                    }
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                    fail( ok, "SVG export threw: " + e.getMessage() );
                }
            } );

            // (4) a gap draws a faint line (not blank), so the alignment's extent is visible: the all-gap column in
            // ALIGN (index 10) adds faint-gray ink when the alignment is shown
            final int[] gap_off = { 0 }, gap_on = { 0 };
            SwingUtilities.invokeAndWait( () -> {
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.getOptions().setMsaColumnWidth( 12 );
                tp.setPreferredSize( new java.awt.Dimension( 760, 460 ) );
                tp.setSize( 760, 460 );
                mf[ 0 ].showWhole();
                tp.calcParametersForPainting( 760, 460 );
                tp.getOptions().setShowMsa( false );
                gap_off[ 0 ] = countFaintGray(
                        AptxUtil.renderPhylogenyToImage( 760, 460, tp, tp.getOptions(), false, 1, false ) );
                tp.getOptions().setShowMsa( true );
                gap_on[ 0 ] = countFaintGray(
                        AptxUtil.renderPhylogenyToImage( 760, 460, tp, tp.getOptions(), false, 1, false ) );
            } );
            if ( gap_on[ 0 ] <= ( gap_off[ 0 ] + 30 ) ) {
                fail( ok, "gaps must draw a faint line (" + gap_off[ 0 ] + " -> " + gap_on[ 0 ] + ")" );
            }

            // (5) the real start/end boundary lines show ONLY at the true edges, so a true edge is distinguishable
            // from a scroll cutoff (fits -> both; offset 0 -> start only; a middle offset -> neither; the end -> end only)
            SwingUtilities.invokeAndWait( () -> {
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.getOptions().setShowMsa( true );
                final int total = ALIGN[ 0 ].length();
                // (a) the whole alignment FITS -> both boundaries bracket it
                tp.getOptions().setMsaColumnWidth( 6 );
                tp.setPreferredSize( new java.awt.Dimension( 900, 460 ) );
                tp.setSize( 900, 460 );
                mf[ 0 ].showWhole();
                tp.calcParametersForPainting( 900, 460 );
                tp.setMsaColumnOffset( 0 );
                if ( !tp.isMsaScrollableForTest()
                        && ( !tp.msaStartBoundaryVisibleForTest() || !tp.msaEndBoundaryVisibleForTest() ) ) {
                    fail( ok, "a fully-visible alignment must show BOTH start and end boundary lines" );
                }
                // (b) narrow the window so it scrolls, then check the boundary logic at the edges vs the middle
                tp.getOptions().setMsaColumnWidth( 40 );
                tp.setSize( 700, 460 );
                mf[ 0 ].showWhole();
                tp.calcParametersForPainting( 700, 460 );
                if ( tp.isMsaScrollableForTest() ) {
                    final int vis = tp.getMsaVisibleColumnsForTest();
                    final int maxoff = total - vis;
                    tp.setMsaColumnOffset( 0 );
                    if ( !tp.msaStartBoundaryVisibleForTest() || tp.msaEndBoundaryVisibleForTest() ) {
                        fail( ok, "at the start only the START boundary must show" );
                    }
                    tp.setMsaColumnOffset( maxoff );
                    if ( tp.msaStartBoundaryVisibleForTest() || !tp.msaEndBoundaryVisibleForTest() ) {
                        fail( ok, "at the end only the END boundary must show" );
                    }
                    if ( maxoff >= 2 ) { // a genuine middle position leaves BOTH edges out of the window
                        tp.setMsaColumnOffset( 1 );
                        if ( tp.msaStartBoundaryVisibleForTest() || tp.msaEndBoundaryVisibleForTest() ) {
                            fail( ok, "a middle-scrolled window must show NEITHER boundary (both edges are cutoffs)" );
                        }
                    }
                }
                tp.getOptions().setMsaColumnWidth( 10 );
                tp.setMsaColumnOffset( 0 );
            } );

            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void checkLayout( final MainFrame mf, final TreePanel tp, final boolean[] ok, final String name,
                                     final Options.PHYLOGENY_GRAPHICS_TYPE type, final Options.TREE_ORIENTATION orient,
                                     final int w, final int h, final boolean expect_msa ) throws Exception {
        final int[] off = { 0 }, on = { 0 };
        SwingUtilities.invokeAndWait( () -> {
            tp.setPhylogenyGraphicsType( type );
            tp.setTreeOrientation( orient );
            tp.setPreferredSize( new java.awt.Dimension( w, h ) );
            tp.setSize( w, h );
            mf.showWhole();
            tp.calcParametersForPainting( w, h );
            tp.getOptions().setShowMsa( false );
            off[ 0 ] = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, tp.getOptions(), false, 1, false ) );
            tp.getOptions().setShowMsa( true );
            on[ 0 ] = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, tp.getOptions(), false, 1, false ) );
        } );
        if ( expect_msa ) {
            if ( on[ 0 ] <= ( off[ 0 ] + 300 ) ) {
                fail( ok, "the alignment must add colored cells in " + name + " (" + off[ 0 ] + " -> " + on[ 0 ] + ")" );
            }
        }
        else if ( on[ 0 ] > ( off[ 0 ] + 300 ) ) {
            fail( ok, "the alignment must NOT draw in " + name + " (" + off[ 0 ] + " -> " + on[ 0 ] + ")" );
        }
    }

    /** Whether two same-size renders differ in a non-trivial number of pixels (the alignment window shows other columns). */
    private static boolean imagesDiffer( final BufferedImage a, final BufferedImage b ) {
        if ( ( a.getWidth() != b.getWidth() ) || ( a.getHeight() != b.getHeight() ) ) {
            return true;
        }
        int diff = 0;
        for( int x = 0; x < a.getWidth(); ++x ) {
            for( int y = 0; y < a.getHeight(); ++y ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    if ( ++diff > 200 ) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /** Light-gray, near-neutral pixels -- the faint gap lines land here; the tree's antialiasing is identical on vs
     *  off so it cancels in the delta, isolating the gap lines. */
    private static int countFaintGray( final BufferedImage img ) {
        int n = 0;
        for( int x = 0; x < img.getWidth(); ++x ) {
            for( int y = 0; y < img.getHeight(); ++y ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, gc = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( ( Math.max( r, Math.max( gc, b ) ) - Math.min( r, Math.min( gc, b ) ) ) <= 12 ) && ( r >= 180 )
                        && ( r <= 245 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Colored (non-grayscale) pixels -- the tree is grayscale, so this isolates the colored alignment cells. */
    private static int countSaturated( final BufferedImage img ) {
        int n = 0;
        for( int x = 0; x < img.getWidth(); ++x ) {
            for( int y = 0; y < img.getHeight(); ++y ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, gc = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( Math.max( r, Math.max( gc, b ) ) - Math.min( r, Math.min( gc, b ) ) ) > 40 ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static Phylogeny alignedTree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < ALIGN.length; ++i ) {
            root.addAsChild( alignedTip( "seq_" + ( i + 1 ), ALIGN[ i ] ) );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode alignedTip( final String name, final String seq ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final Sequence s = new Sequence();
        s.setMolecularSequence( seq );
        s.setMolecularSequenceAligned( true );
        n.getNodeData().addSequence( s );
        return n;
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [MsaTrackRenderTest] " + message );
        ok[ 0 ] = false;
    }

    private MsaTrackRenderTest() {
        // not instantiable
    }
}
