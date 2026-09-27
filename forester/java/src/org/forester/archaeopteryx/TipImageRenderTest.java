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

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.imageio.ImageIO;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;

/**
 * Integration test for TIP IMAGES: on a real {@link MainFrame}/{@link TreePanel} it puts a distinctly-coloured local
 * image at each ragged tip and checks the pictures render (add coloured ink over the grayscale tree) in EVERY layout
 * -- rectangular root-left, root-top (vertical), circular, unrooted -- that the label is shifted to make room, that a
 * broken reference doesn't crash, and that a VECTOR export embeds the images. Guarded to a no-op on a headless box
 * (needs FlatLaf via {@code createInstance}); run standalone or as part of the non-headless suite.
 */
public final class TipImageRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TipImageRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final File dir = java.nio.file.Files.createTempDirectory( "aptx-tipimgrender" ).toFile();
            dir.deleteOnExit();
            // five distinctly-coloured, ragged-aspect local images
            writeSolid( new File( dir, "r.png" ), 90, 50, Color.RED );
            writeSolid( new File( dir, "g.png" ), 60, 60, Color.GREEN );
            writeSolid( new File( dir, "b.png" ), 120, 45, Color.BLUE );
            writeSolid( new File( dir, "c.png" ), 55, 80, Color.CYAN );
            final Phylogeny phy = raggedTree();

            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "tipimg" ) );
            final boolean[] ok = { true };
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            SwingUtilities.invokeAndWait( () -> {
                tp.setTreeFile( new File( dir, "tree.xml" ) ); // base dir for the relative image paths
                tp.getOptions().setGraphicsExportWhiteBackground( true );
                tp.getOptions().setShowOverview( false );
                tp.setOvOn( false );
                tp.getOptions().setTipImageSize( 44 );
            } );
            tp.preloadTipImagesForTest( 8000 ); // wait (off-EDT) for the local images to load, so renders are deterministic

            // (1) the label advance is > 0 for an imaged tip when tip images are on, and 0 when off
            SwingUtilities.invokeAndWait( () -> {
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                final PhylogenyNode imaged = tp.getPhylogeny().getExternalNodes().get( 0 );
                tp.getOptions().setShowTipImages( false );
                if ( tp.tipImageAdvanceForTest( imaged ) != 0 ) {
                    fail( ok, "with tip images OFF the label advance must be 0" );
                }
                tp.getOptions().setShowTipImages( true );
                if ( tp.tipImageAdvanceForTest( imaged ) <= 0 ) {
                    fail( ok, "with tip images ON an imaged tip must reserve a positive label advance" );
                }
            } );

            // (1b) the reserved label reach GROWS when tip images are on, so annotation columns / clade bands anchor
            // past the images instead of overprinting the shifted labels
            SwingUtilities.invokeAndWait( () -> {
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.setSize( 760, 480 );
                tp.getOptions().setShowTipImages( false );
                tp.calculateLongestExtNodeInfo();
                final int off = tp.getLongestExtNodeInfo();
                tp.getOptions().setShowTipImages( true );
                tp.calculateLongestExtNodeInfo();
                final int on = tp.getLongestExtNodeInfo();
                if ( on <= off ) {
                    fail( ok, "tip images must widen the reserved label reach (" + off + " -> " + on + ")" );
                }
            } );

            // (2) images render (add coloured ink) in every layout -- off vs on
            checkLayout( mf[ 0 ], tp, ok, "root-left", Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR,
                    Options.TREE_ORIENTATION.ROOT_LEFT, 760, 480 );
            checkLayout( mf[ 0 ], tp, ok, "root-top", Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR,
                    Options.TREE_ORIENTATION.ROOT_TOP, 600, 600 );
            checkLayout( mf[ 0 ], tp, ok, "circular", Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                    Options.TREE_ORIENTATION.ROOT_LEFT, 600, 600 );
            checkLayout( mf[ 0 ], tp, ok, "unrooted", Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED,
                    Options.TREE_ORIENTATION.ROOT_LEFT, 600, 600 );

            // (2a) an image is one mark with its label: on a fan too dense for its names, a tip whose label the rule
            // hides shows NO image either, and every image drawn sits inside the box the paint recorded for it --
            // an orphan image under a neighbour's name would read as that neighbour's (a review find, 2026-09-27:
            // the image used to be painted before the label's claim). Saturated pixels are the images' alone.
            final int[] star_hidden = { 0 }, star_outside = { 0 }, star_saturated = { 0 };
            SwingUtilities.invokeAndWait( () -> {
                final PhylogenyNode star_root = new PhylogenyNode();
                final String[] imgs = { "r.png", "g.png", "b.png", "c.png" };
                for( int i = 0; i < 48; ++i ) {
                    star_root.addAsChild( imgTip( "lineage_" + i, imgs[ i % 4 ], ( ( i % 2 ) == 0 ) ? 0.2 : 0.6 ) );
                }
                final Phylogeny star = new Phylogeny();
                star.setRoot( star_root );
                star.setRooted( true );
                star.recalculateNumberOfExternalDescendants( false );
                tp.setTree( star );
                tp.setTreeFile( new File( dir, "tree.xml" ) );
                tp.recalculateMaxDistanceToRoot();
                // everything this block sets is put back in its finally: (2b) and (3) below run under the settings
                // they find, so leaving these behind would make them depend on this block having run
                final Options.PHYLOGENY_GRAPHICS_TYPE type0 = tp.getPhylogenyGraphicsType();
                final Options.NODE_LABEL_DIRECTION dir0 = tp.getOptions().getNodeLabelDirection();
                final boolean hide0 = tp.shows( DisplayOption.DYNAMICALLY_HIDE_DATA );
                final boolean imgs0 = tp.getOptions().isShowTipImages();
                final boolean names0 = tp.shows( DisplayOption.SHOW_NODE_NAMES );
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                tp.getOptions().setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.RADIAL );
                tp.getControlPanel().setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                tp.getOptions().setShowTipImages( true );
                tp.setShows( DisplayOption.SHOW_NODE_NAMES, true );
                final int[] prior = tp.layoutForExportSize( 600, 600 );
                try {
                    tp.setRecordLabelBoxesForTest( true );
                    final BufferedImage img = AptxUtil.renderPhylogenyToImage( 600, 600, tp, tp.getOptions(), false, 1,
                            false );
                    star_hidden[ 0 ] = tp.labelsHiddenForTest();
                    final java.util.List<double[]> boxes = tp.labelBoxesForTest();
                    // the recorded boxes -- names and images alike -- never overlap one another: an image that was
                    // drawn but not RESERVED would let a later name be placed across it, and the seam alone would
                    // not show it (the image's box is recorded whether or not it was reserved)
                    for( int i = 0; i < boxes.size(); ++i ) {
                        for( int j = i + 1; j < boxes.size(); ++j ) {
                            final double[] a = boxes.get( i ), b = boxes.get( j );
                            if ( OrientedOccupancy.overlap( a[ 0 ], a[ 1 ], a[ 2 ], a[ 3 ], Math.cos( a[ 4 ] ),
                                    Math.sin( a[ 4 ] ), b[ 0 ], b[ 1 ], b[ 2 ], b[ 3 ], Math.cos( b[ 4 ] ),
                                    Math.sin( b[ 4 ] ) ) ) {
                                ++star_outside[ 0 ]; // counted with the stray pixels: either way something is drawn across something
                            }
                        }
                    }
                    for( int x = 0; x < img.getWidth(); ++x ) {
                        for( int y = 0; y < img.getHeight(); ++y ) {
                            final int rgb = img.getRGB( x, y );
                            final int r = ( rgb >> 16 ) & 0xFF, gc = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                            if ( ( Math.max( r, Math.max( gc, b ) ) - Math.min( r, Math.min( gc, b ) ) ) <= 40 ) {
                                continue;
                            }
                            ++star_saturated[ 0 ];
                            boolean in = false;
                            for( final double[] bx : boxes ) {
                                final double dx = ( x + 0.5 ) - bx[ 0 ], dy = ( y + 0.5 ) - bx[ 1 ];
                                final double u = ( dx * Math.cos( bx[ 4 ] ) ) + ( dy * Math.sin( bx[ 4 ] ) );
                                final double v = ( -dx * Math.sin( bx[ 4 ] ) ) + ( dy * Math.cos( bx[ 4 ] ) );
                                if ( ( Math.abs( u ) <= ( bx[ 2 ] + 1.5 ) ) && ( Math.abs( v ) <= ( bx[ 3 ] + 1.5 ) ) ) {
                                    in = true;
                                    break;
                                }
                            }
                            if ( !in ) {
                                ++star_outside[ 0 ];
                            }
                        }
                    }
                }
                finally {
                    tp.setRecordLabelBoxesForTest( false );
                    tp.restoreLayoutAfterExport( prior );
                    tp.setTree( phy );
                    tp.setTreeFile( new File( dir, "tree.xml" ) );
                    tp.recalculateMaxDistanceToRoot();
                    tp.setPhylogenyGraphicsType( type0 );
                    tp.getOptions().setNodeLabelDirection( dir0 );
                    tp.getControlPanel().setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, hide0 );
                    tp.getOptions().setShowTipImages( imgs0 );
                    tp.setShows( DisplayOption.SHOW_NODE_NAMES, names0 );
                }
            } );
            if ( star_hidden[ 0 ] == 0 ) {
                fail( ok, "precondition: the 48-tip imaged fan must make the label rule hide some names (0 hidden)" );
            }
            else if ( star_saturated[ 0 ] < 200 ) {
                fail( ok, "precondition: some images must be drawn on the fan (" + star_saturated[ 0 ] + " saturated px)" );
            }
            else if ( star_outside[ 0 ] > 0 ) {
                fail( ok, "an image is drawn outside every recorded label/image box -- a hidden tip's image, or an "
                        + "image not reserved with its label (" + star_outside[ 0 ] + " of " + star_saturated[ 0 ]
                        + " saturated px)" );
            }

            // (2b) "Show External Data" off suppresses the tip images (parity with the radial layout)
            final int[] shown = { 0 }, hidden = { 0 };
            SwingUtilities.invokeAndWait( () -> {
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.setPreferredSize( new java.awt.Dimension( 760, 480 ) );
                tp.setSize( 760, 480 );
                mf[ 0 ].showWhole();
                tp.calcParametersForPainting( 760, 480 );
                tp.getOptions().setShowTipImages( true );
                tp.getControlPanel().setShowExternalDataForTest( true );
                shown[ 0 ] = countSaturated(
                        AptxUtil.renderPhylogenyToImage( 760, 480, tp, tp.getOptions(), false, 1, false ) );
                tp.getControlPanel().setShowExternalDataForTest( false );
                hidden[ 0 ] = countSaturated(
                        AptxUtil.renderPhylogenyToImage( 760, 480, tp, tp.getOptions(), false, 1, false ) );
                tp.getControlPanel().setShowExternalDataForTest( true );
            } );
            if ( hidden[ 0 ] >= ( shown[ 0 ] - 200 ) ) {
                fail( ok, "with 'Show External Data' off the tip images must not draw (" + shown[ 0 ] + " -> "
                        + hidden[ 0 ] + ")" );
            }

            // (3) a VECTOR export embeds the images (SVG <image> elements), i.e. drawImage works through the export path
            SwingUtilities.invokeAndWait( () -> {
                try {
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    tp.getOptions().setShowTipImages( true );
                    tp.setPreferredSize( new java.awt.Dimension( 760, 480 ) );
                    tp.setSize( 760, 480 );
                    mf[ 0 ].showWhole();
                    tp.calcParametersForPainting( 760, 480 );
                    final File svg = new File( dir, "out.svg" );
                    VectorGraphicsExporter.writePhylogenyToVectorGraphicsFile( svg.getAbsolutePath(), tp, 760, 480,
                            VectorGraphicsExporter.Format.SVG, true, true );
                    final String content = new String( java.nio.file.Files.readAllBytes( svg.toPath() ), "UTF-8" );
                    final int images = content.split( "<image", -1 ).length - 1;
                    if ( images < 4 ) {
                        fail( ok, "the SVG export must embed the tip images (>=4 <image> elements), got " + images );
                    }
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                    fail( ok, "SVG export threw: " + e.getMessage() );
                }
            } );

            // (4) a broken image reference must not crash the render (a faint placeholder is drawn instead)
            SwingUtilities.invokeAndWait( () -> {
                tp.getPhylogeny().getExternalNodes().get( 1 ).getNodeData().getProperties()
                        .addProperty( new Property( "data:image", "no_such_file.png", "", "xsd:string", AppliesTo.NODE ) );
            } );
            try {
                AptxUtil.renderPhylogenyToImage( 600, 400, tp, tp.getOptions(), false, 1, false );
            }
            catch ( final Throwable t ) {
                t.printStackTrace();
                fail( ok, "a broken image reference must not crash the render: " + t );
            }

            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Renders {@code layout} with tip images OFF then ON and asserts ON adds coloured (saturated) ink -- the tree is
     *  grayscale, so only the coloured images (and any placeholder) contribute. */
    private static void checkLayout( final MainFrame mf, final TreePanel tp, final boolean[] ok, final String name,
                                     final Options.PHYLOGENY_GRAPHICS_TYPE type, final Options.TREE_ORIENTATION orient,
                                     final int w, final int h ) throws Exception {
        final int[] off = { 0 }, on = { 0 };
        SwingUtilities.invokeAndWait( () -> {
            tp.setPhylogenyGraphicsType( type );
            tp.setTreeOrientation( orient );
            tp.setPreferredSize( new java.awt.Dimension( w, h ) );
            tp.setSize( w, h );
            mf.showWhole();
            // lay the tree out FOR the export canvas, as a real fixed-size export does: the radial ring is a square
            // of its own diameter, and one fitted to the window is larger than a 600 px canvas -- with every circular
            // phylogram's labels (and images) now on that ring, an unfitted ring puts them all off the canvas (the
            // ragged tree used to survive this because its images sat at the tips, some of which were inside)
            final int[] prior = tp.layoutForExportSize( w, h );
            try {
                tp.getOptions().setShowTipImages( false );
                off[ 0 ] = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, tp.getOptions(), false, 1, false ) );
                tp.getOptions().setShowTipImages( true );
                on[ 0 ] = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, tp.getOptions(), false, 1, false ) );
            }
            finally {
                tp.restoreLayoutAfterExport( prior );
            }
        } );
        if ( on[ 0 ] <= ( off[ 0 ] + 200 ) ) {
            fail( ok, "tip images must add coloured ink in the " + name + " layout (" + off[ 0 ] + " -> " + on[ 0 ] + ")" );
        }
    }

    /** Coloured (non-grayscale) pixels -- the tree is grayscale, so this isolates the coloured tip images. */
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

    private static Phylogeny raggedTree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( clade( 0.1, imgTip( "T_rex", "r.png", 0.9 ), imgTip( "Raptor", "g.png", 0.4 ) ) );
        root.addAsChild( clade( 0.1, imgTip( "Archie", "b.png", 0.2 ), imgTip( "Trike", "c.png", 0.7 ) ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode imgTip( final String name, final String img, final double bl ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( bl );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( "data:image", img, "", "xsd:string", AppliesTo.NODE ) );
        n.getNodeData().setProperties( pl );
        return n;
    }

    private static PhylogenyNode clade( final double bl, final PhylogenyNode... kids ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setDistanceToParent( bl );
        for( final PhylogenyNode k : kids ) {
            n.addAsChild( k );
        }
        return n;
    }

    private static void writeSolid( final File f, final int w, final int h, final Color c ) throws Exception {
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( c );
        g.fillRect( 0, 0, w, h );
        g.dispose();
        ImageIO.write( img, "png", f );
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [TipImageRenderTest] " + message );
        ok[ 0 ] = false;
    }

    private TipImageRenderTest() {
        // not instantiable
    }
}
