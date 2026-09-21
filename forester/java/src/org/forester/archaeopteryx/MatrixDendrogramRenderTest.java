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
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.nio.file.Files;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.tools.NodeDataImporter;
import org.forester.phylogeny.Phylogeny;

/**
 * The clustergram's own rendering, on the fixture that exercises it: the COLUMN dendrogram drawn above the matrix,
 * and the seamless tiling of the matrix cells themselves. Both live here because this is the panel that shows them
 * -- 18 columns over 50 tips, cells small enough for a hairline between two of them to be visible.
 * <p>
 * The structure itself -- merges, heights and leaf order -- is pinned against R in {@link MatrixColumnOrderTest}.
 * What is tested here is what reaches the screen:
 * <ul>
 * <li>it is drawn in the two Clustered modes and in NO other, and the band it needs is 0 px in the others, so a
 * mode without a clustering behind it pays no space for one;</li>
 * <li>it disappears the moment the drawn order stops being the clustering's -- dragging a column is the ordinary
 * way to get there, and a dendrogram over a hand-reordered matrix would draw connectors that cross;</li>
 * <li>it is really INK, measured in the band it claims, and the band is really RESERVED -- the first tip moves down
 * by it rather than the lines landing on the headers;</li>
 * <li>it is part of the figure, so it is in an EXPORT too (unlike the drag marker, which is screen-only);</li>
 * <li>not in CIRCULAR -- a dendrogram over concentric rings has no sensible geometry (Christian, 2026-09-21:
 * "circular does not need it") -- and not in UNROOTED, which draws no columns at all.</li>
 * </ul>
 * Headful; a green no-op when headless.
 */
public final class MatrixDendrogramRenderTest {

    private static final int W = 1300;
    private static final int H = 900;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MatrixDendrogramRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        heightMapping( ok );
        if ( GraphicsEnvironment.isHeadless() ) {
            return ok[ 0 ];
        }
        return dendrogramOk() && ok[ 0 ];
    }

    /**
     * The height-to-y mapping, on its own: LINEAR in the height. Even spacing by merge rank would look tidier and
     * would claim a structure the data does not have, so it is pinned here rather than left to the eye -- the render
     * below can only see that SOMETHING was drawn in the band.
     */
    private static void heightMapping( final boolean[] ok ) {
        final double base = 100;
        final int band = 40;
        if ( TreePanel.dendrogramMergeY( base, band, 0, 10 ) != 100 ) {
            fail( ok, "a merge at height 0 sits on the baseline" );
        }
        if ( TreePanel.dendrogramMergeY( base, band, 10, 10 ) != 60 ) {
            fail( ok, "the tallest merge sits at the top of the band" );
        }
        if ( TreePanel.dendrogramMergeY( base, band, 2.5, 10 ) != 90 ) {
            fail( ok, "a quarter-height merge must be a QUARTER of the way up (linear), got "
                    + TreePanel.dendrogramMergeY( base, band, 2.5, 10 ) );
        }
        if ( TreePanel.dendrogramMergeY( base, band, Double.POSITIVE_INFINITY, 10 ) != 60 ) {
            fail( ok, "a merge with no finite height goes to the top, above every measurable one" );
        }
        if ( TreePanel.dendrogramMergeY( base, band, 5, 0 ) != 60 ) {
            fail( ok, "with no positive maximum there is nothing to scale against: the top" );
        }
    }

    private static boolean fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [MatrixDendrogramRenderTest] " + msg );
        ok[ 0 ] = false;
        return false;
    }

    /** The demo tree with its table joined on, as File &gt; Demo Trees leaves it: 18 matrix columns. */
    private static Phylogeny demo() throws Exception {
        final File dir = new File( System.getProperty( "user.dir" ), "forester/demo" );
        final Phylogeny phy = FigureRenderer.readTrees( new File( dir, "sparse-accessory-genome.xml" ) )[ 0 ];
        final NodeDataImporter.Table table = NodeDataImporter
                .parseTable( Files.readString( new File( dir, "sparse-accessory-genome.tsv" ).toPath() ) );
        NodeDataImporter.apply( phy, table, table.defaultKeyColumn(), NodeDataImporter.MatchBy.TIP_NAME );
        return phy;
    }

    private static boolean dendrogramOk() {
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        try {
            final Phylogeny phy = demo();
            SwingUtilities.invokeAndWait( () -> {
                mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(),
                                                               "column dendrogram" );
                ( ( JFrame ) mf[ 0 ] ).setLocation( -32000, -32000 );
                ( ( JFrame ) mf[ 0 ] ).setSize( W, H );
                ( ( JFrame ) mf[ 0 ] ).setVisible( true );
            } );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    mf[ 0 ].applyClustergramPreset();
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    // FIRST, on the pristine clustergram: the checks below deliberately mutate the panel (one of
                    // them rewrites every cell value), and the seam detector reads the real data's cell colours
                    seamlessWhenScaled( ok, tp );
                    modes( ok, tp );
                    ink( ok, tp );
                    leavesUnderTheirColumns( ok, tp );
                    reserved( ok, tp );
                    staleCache( ok, tp );
                    vertical( ok, tp );
                    vanishesOnDrag( ok, tp );
                    notRadial( ok, tp );
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                    t.printStackTrace();
                }
                finally {
                    ( ( JFrame ) mf[ 0 ] ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    /** Only the two Clustered modes have a clustering behind them -- and only they may pay for a band. */
    private static void modes( final boolean[] ok, final TreePanel tp ) {
        for( final MatrixColumnOrder.Mode m : MatrixColumnOrder.Mode.values() ) {
            tp.setMatrixColumnOrder( m );
            paint( tp );
            final boolean clustered = ( m == MatrixColumnOrder.Mode.CLUSTERED )
                    || ( m == MatrixColumnOrder.Mode.CLUSTERED_PRESENCE );
            final MatrixColumnOrder.Dendrogram d = tp.matrixColumnDendrogram();
            if ( ( d != null ) != clustered ) {
                fail( ok, m.label() + ": dendrogram present = " + ( d != null ) + ", expected " + clustered );
            }
            if ( ( tp.matrixDendrogramBandHeight() > 0 ) != clustered ) {
                fail( ok, m.label() + ": reserves " + tp.matrixDendrogramBandHeight()
                        + " px for a dendrogram it does not draw" );
            }
            if ( clustered && ( d.stages() != 17 ) ) {
                fail( ok, m.label() + ": 18 columns must give 17 merges, got " + d.stages() );
            }
        }
        // the two clusterings are different trees, or the mode switch would be drawing the same picture twice
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        paint( tp );
        final double[] euclid = tp.matrixColumnDendrogram().height().clone();
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED_PRESENCE );
        paint( tp );
        if ( java.util.Arrays.equals( euclid, tp.matrixColumnDendrogram().height() ) ) {
            fail( ok, "the two Clustered modes must draw different dendrograms on this fixture" );
        }
    }

    /** There is really ink in the band, and really none there when the mode has no clustering behind it. */
    private static void ink( final boolean[] ok, final TreePanel tp ) {
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        final Rectangle band = bandOnScreen( tp );
        if ( band == null ) {
            fail( ok, "the clustered clustergram must report a dendrogram band" );
            return;
        }
        final int before = tp.matrixDendrogramPaintsForTest();
        final int drawn = inkIn( render( tp ), band, screenBg( tp ) );
        if ( tp.matrixDendrogramPaintsForTest() <= before ) {
            fail( ok, "the dendrogram painter was never reached on a clustered clustergram" );
        }
        if ( drawn < 50 ) {
            fail( ok, "only " + drawn + " marked pixels in the dendrogram band " + band + " -- nothing was drawn" );
        }
        // an EXPORT is the same figure: unlike the drag marker, this belongs in the file
        final int exported = inkIn( exportRender( tp ), band, Color.WHITE.getRGB() );
        if ( exported < 50 ) {
            fail( ok, "the dendrogram must be in an export too, got " + exported + " marked pixels" );
        }
        // ...drawn at the FIGURE's line weight, not a hardcoded one. Only a weight measurement sees this: a fixed
        // stroke still puts ink in the band, it just comes out at twice the branch weight in a default PDF.
        final float saved = tp.getOptions().getPdfLineWidth();
        try {
            tp.getOptions().setPdfLineWidth( 0.5f );
            final int thin = inkIn( exportRender( tp ), band, Color.WHITE.getRGB() );
            tp.getOptions().setPdfLineWidth( 4f );
            final int thick = inkIn( exportRender( tp ), band, Color.WHITE.getRGB() );
            if ( thick <= ( thin * 1.5 ) ) {
                fail( ok, "the dendrogram must follow the export's line width: " + thin + " px at 0.5 vs " + thick
                        + " at 4.0 -- it is drawn at a fixed weight" );
            }
        }
        finally {
            tp.getOptions().setPdfLineWidth( saved );
        }
        // The neighbouring case: a mode with no clustering. Its strip has to be measured in ITS OWN layout -- the
        // band is a layout reserve, so dropping it moves the whole tree up and a rectangle from the other mode would
        // land on the headers (measured: 558 header pixels, which is what this check used to "catch").
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.TABLE );
        paint( tp );
        final int quiet = tp.matrixDendrogramPaintsForTest();
        render( tp );
        if ( tp.matrixDendrogramPaintsForTest() != quiet ) {
            fail( ok, "the dendrogram painter must not run at all in a mode with no clustering behind it" );
        }
        final Rectangle same_strip = tp.matrixHeaderStripForTest( band.height );
        if ( same_strip == null ) {
            fail( ok, "the strip above the headers must be measurable in Same as Table too" );
        }
        else {
            final int none = inkIn( render( tp ), same_strip, screenBg( tp ) );
            if ( none > 0 ) {
                fail( ok, "Same as Table draws no dendrogram, yet " + none + " pixels are marked above its headers in "
                        + same_strip );
            }
        }
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        paint( tp );
    }

    /**
     * Every leaf stem stands under its own column. The x positions are checked against
     * {@link TreePanel#annotationColumnGrabPointForTest} -- the geometry a header CLICK uses, which the dendrogram
     * knows nothing about -- so a leaf mapping that drifted (counting the colour strips as leaves, say, or reading
     * the clustering's leaf order the wrong way round) cannot pass by moving the band and the stems together.
     */
    private static void leavesUnderTheirColumns( final boolean[] ok, final TreePanel tp ) {
        // The demo is all MATRIX columns, so a leaf mapping that counted the OTHER column types would pass unnoticed
        // (measured). Turn the first field into a colour strip: 17 leaves among 18 drawn columns, and the strip must
        // get none of them.
        final java.util.List<AnnotationColumns.ColumnSpec> original =
                new java.util.ArrayList<AnnotationColumns.ColumnSpec>( tp.getAnnotationColumnSpecs() );
        final java.util.List<AnnotationColumns.ColumnSpec> mixed =
                new java.util.ArrayList<AnnotationColumns.ColumnSpec>( original );
        mixed.set( 0, new AnnotationColumns.ColumnSpec( mixed.get( 0 )._ref,
                                                        AnnotationColumns.Type.COLOR_STRIP ) );
        tp.setAnnotationColumns( mixed );
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        paint( tp );
        final Rectangle band = tp.matrixDendrogramBandForTest();
        final int n = MatrixColumnOrder.matrixRefs( tp.getAnnotationColumnSpecs() ).size();
        if ( ( band == null ) || ( n != 17 ) ) {
            fail( ok, "fixture: 17 matrix columns beside one colour strip expected, got " + n );
            return;
        }
        final MatrixColumnOrder.Dendrogram d = tp.matrixColumnDendrogram();
        if ( ( d == null ) || ( d.stages() < 1 ) ) {
            fail( ok, "fixture: a dendrogram is needed here" );
            return;
        }
        final BufferedImage img = render( tp );
        final int bg = screenBg( tp );
        int missing = 0;
        for( int c = 0; c < n; ++c ) {
            final java.awt.Point header = tp.annotationColumnGrabPointForTest( matrixColumnToDrawn( tp, c ) );
            boolean stem = false;
            for( int dx = -2; ( dx <= 2 ) && !stem; ++dx ) {
                for( int dy = 1; ( dy <= 4 ) && !stem; ++dy ) { // just above the baseline: the leaf stems
                    final int x = header.x + dx;
                    final int y = ( band.y + band.height ) - dy;
                    if ( ( x >= 0 ) && ( x < img.getWidth() ) && ( y >= 0 ) && ( y < img.getHeight() )
                            && ( img.getRGB( x, y ) != bg ) ) {
                        stem = true;
                    }
                }
            }
            if ( !stem ) {
                ++missing;
            }
        }
        if ( missing > 0 ) {
            fail( ok, missing + " of " + n + " columns have no dendrogram leaf standing over them" );
        }
        // ...and the leaves are the RIGHT way round. A stem over every column says nothing about which leaf is
        // where -- a mirrored mapping still puts one over each (measured). The FIRST merge always joins two
        // singletons, so where its horizontal bar sits is a statement about which two columns the clustering says
        // are closest, checked against the header geometry the dendrogram knows nothing about.
        final int k1 = positionOf( d, -d.left()[ 0 ] - 1 );
        final int k2 = positionOf( d, -d.right()[ 0 ] - 1 );
        if ( ( d.left()[ 0 ] >= 0 ) || ( d.right()[ 0 ] >= 0 ) || ( k1 < 0 ) || ( k2 < 0 ) ) {
            fail( ok, "the first merge must join two single columns, got " + d.left()[ 0 ] + " and " + d.right()[ 0 ] );
            return;
        }
        double max_h = 0;
        for( final double h : d.height() ) {
            if ( Double.isFinite( h ) && ( h > max_h ) ) {
                max_h = h;
            }
        }
        final int x1 = tp.annotationColumnGrabPointForTest( matrixColumnToDrawn( tp, k1 ) ).x;
        final int x2 = tp.annotationColumnGrabPointForTest( matrixColumnToDrawn( tp, k2 ) ).x;
        final int bar_y = ( int ) Math.round( TreePanel.dendrogramMergeY( band.y + band.height, band.height,
                                                                          d.height()[ 0 ], max_h ) );
        final int mid = ( x1 + x2 ) / 2;
        boolean bar = false;
        for( int dy = -2; ( dy <= 2 ) && !bar; ++dy ) {
            bar |= ( img.getRGB( mid, Math.max( 0, Math.min( img.getHeight() - 1, bar_y + dy ) ) ) != bg );
        }
        if ( !bar ) {
            fail( ok, "the first merge's bar must span the two columns the clustering joined first (x " + x1 + ".."
                    + x2 + " at y " + bar_y + ") -- a mirrored leaf mapping puts it elsewhere" );
        }
        tp.setAnnotationColumns( original ); // hand the fixture back: the checks after this one assume all-MATRIX
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        paint( tp );
    }

    /** Where leaf {@code leaf} (an index into the clustering) sits among the drawn matrix columns, or -1. */
    private static int positionOf( final MatrixColumnOrder.Dendrogram d, final int leaf ) {
        for( int k = 0; k < d.order().length; ++k ) {
            if ( d.order()[ k ] == leaf ) {
                return k;
            }
        }
        return -1;
    }

    /** The drawn-column index of the {@code k}-th MATRIX column (the other types are not leaves). */
    private static int matrixColumnToDrawn( final TreePanel tp, final int k ) {
        int seen = 0;
        final java.util.List<AnnotationColumns.ColumnSpec> specs = tp.getAnnotationColumnSpecs();
        for( int i = 0; i < specs.size(); ++i ) {
            if ( specs.get( i )._type == AnnotationColumns.Type.MATRIX ) {
                if ( seen == k ) {
                    return i;
                }
                ++seen;
            }
        }
        return -1;
    }

    /** The band is RESERVED, not overlaid: turning the dendrogram on pushes the first tip down by about its height. */
    private static void reserved( final boolean[] ok, final TreePanel tp ) {
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.TABLE );
        paint( tp );
        final float without = firstTipY( tp );
        final float spread_without = tipSpread( tp );
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        paint( tp );
        final float with = firstTipY( tp );
        final float spread_with = tipSpread( tp );
        final int band = tp.matrixDendrogramBandHeight();
        // Reserving the band is TWO things, and each has to be asserted: the tree is SHIFTED down so it starts below
        // the band (the root's Ycoord), and it is COMPRESSED so it still ends inside the canvas (the ydist budget).
        // Checking only the shift passes while the tree overflows the bottom by the height of the band -- measured,
        // by removing the compression and watching this test stay green.
        if ( with <= without ) {
            fail( ok, "reserving the band must push the tips down: first tip at " + without + " without it, " + with
                    + " with it" );
        }
        else if ( ( with - without ) < ( band / 2.0 ) ) {
            fail( ok, "the tips moved only " + ( with - without ) + " for a band of " + band
                    + " -- the band is being overlaid rather than reserved" );
        }
        if ( ( spread_without - spread_with ) < ( band / 2.0 ) ) {
            fail( ok, "the tip spread must SHRINK by about the band (" + band + "): " + spread_without + " -> "
                    + spread_with + " -- the tree is shifted down without being compressed, so it now runs past the"
                    + " bottom of the canvas" );
        }
    }

    /** The canvas background the screen render really uses -- the tree colour set's, not a corner pixel. */
    private static int screenBg( final TreePanel tp ) {
        return tp.getTreeColorSet().getBackgroundColor().getRGB();
    }

    private static float tipSpread( final TreePanel tp ) {
        float min = Float.MAX_VALUE;
        float max = -Float.MAX_VALUE;
        for( final org.forester.phylogeny.PhylogenyNode t : tp.getPhylogeny().getExternalNodes() ) {
            min = Math.min( min, t.getYcoord() );
            max = Math.max( max, t.getYcoord() );
        }
        return max - min;
    }

    /** Move a column by hand and the drawn order is no longer the clustering's: the dendrogram has to go. */
    private static void vanishesOnDrag( final boolean[] ok, final TreePanel tp ) {
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        paint( tp );
        final Rectangle band = bandOnScreen( tp );
        if ( ( band == null ) || ( tp.matrixColumnDendrogram() == null ) ) {
            fail( ok, "fixture: a dendrogram must be showing before the drag that has to remove it" );
            return;
        }
        if ( !tp.moveAnnotationColumn( 0, 5 ) ) {
            fail( ok, "fixture: the column move did not happen" );
            return;
        }
        paint( tp );
        if ( tp.matrixColumnDendrogram() != null ) {
            fail( ok, "a hand-moved column leaves a matrix the clustering no longer describes: no dendrogram" );
        }
        if ( tp.matrixDendrogramBandHeight() != 0 ) {
            fail( ok, "...and no band either, got " + tp.matrixDendrogramBandHeight() );
        }
        final Rectangle after = tp.matrixHeaderStripForTest( band.height ); // in the layout the drag left behind
        if ( ( after != null ) && ( inkIn( render( tp ), after, screenBg( tp ) ) > 0 ) ) {
            fail( ok, "...and no ink above the headers where the dendrogram used to be, in " + after );
        }
    }

    /**
     * The clustering is cached on (refs, mode), which a change to the DATA leaves untouched -- a subtree, an undo, a
     * tip deletion, a re-import. The self-check that would catch a stale one runs only on a cache MISS, so without
     * an explicit invalidation the old dendrogram is handed back and drawn over the new numbers. Every one of those
     * paths goes through rebuildAnnotationColumns(), which is where the cache is dropped, and this is that.
     */
    private static void staleCache( final boolean[] ok, final TreePanel tp ) {
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        paint( tp );
        final MatrixColumnOrder.Dendrogram before = tp.matrixColumnDendrogram();
        if ( before == null ) {
            fail( ok, "fixture: a dendrogram is needed before the data is changed under it" );
            return;
        }
        // change what the columns are clustered FROM, leaving the refs and the mode exactly as they were
        final java.util.List<String> refs = MatrixColumnOrder.matrixRefs( tp.getAnnotationColumnSpecs() );
        int touched = 0;
        for( final org.forester.phylogeny.PhylogenyNode t : tp.getPhylogeny().getExternalNodes() ) {
            for( final String ref : refs ) {
                for( final org.forester.phylogeny.data.Property pr : t.getNodeData().getProperties()
                        .getProperties( ref ) ) {
                    pr.setValue( String.valueOf( ( touched++ % 5 ) ) );
                }
            }
        }
        if ( touched < 100 ) {
            fail( ok, "fixture: too few cells changed (" + touched + ") to change the clustering" );
            return;
        }
        tp.rebuildAnnotationColumns();
        paint( tp );
        final MatrixColumnOrder.Dendrogram after = tp.matrixColumnDendrogram();
        if ( ( after != null ) && after.equals( before ) ) {
            fail( ok, "the dendrogram of the OLD data is still being drawn over the new: " + after );
        }
    }

    /**
     * A scaled EXPORT tiles its cells with no seam between them.
     * <p>
     * On screen the cells are integer rectangles at scale 1, so two neighbours share a pixel boundary exactly. An
     * export scales user space by a fractional factor: the shared edge then lands in the middle of a device pixel
     * and BOTH neighbours antialias against the background instead of against each other -- a pale hairline grid
     * across the whole matrix in a PDF, while the screen looked right.
     * <p>
     * A seam pixel is one whose neighbours two rows above and below are the SAME cell colour while it is not, so it
     * can only be a gap and never a real boundary between two different cells. The comparison is RELATIVE, against
     * the same picture drawn the screen way at the same scale, for two measured reasons: the detector has a noise
     * floor of ~200 px where antialiased text and branches cross the grid, and the un-bled fill produces tens of
     * thousands (measured here: 7330 against 257). At WHOLE-number scales both collapse to the noise floor, which is
     * why the scale below is deliberately not one -- and the screen render is the fixture check: it must seam, or
     * the detector is blind and the export assertion would pass on any picture at all.
     */
    private static void seamlessWhenScaled( final boolean[] ok, final TreePanel tp ) {
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        paint( tp );
        final double scale = 1.37; // deliberately not a whole number: that is what puts an edge mid-pixel
        final int screen = seamPixels( scaledRender( tp, scale, false ) );
        final int export = seamPixels( scaledRender( tp, scale, true ) );
        if ( screen < 500 ) {
            fail( ok, "fixture: the un-bled fill must seam heavily at a fractional scale (" + screen
                    + " found), or the export check below proves nothing" );
        }
        else if ( export > ( screen / 5 ) ) {
            fail( ok, export + " seam pixels between cells in a scaled export against " + screen
                    + " drawn the screen way -- the cells abut instead of overlapping, so the background shows"
                    + " through every shared edge" );
        }
    }

    /** The panel painted through the export path (or the screen one) into a raster at a fractional scale. */
    private static BufferedImage scaledRender( final TreePanel tp, final double scale, final boolean exporting ) {
        final int w = ( int ) Math.ceil( tp.getWidth() * scale );
        final int h = ( int ) Math.ceil( tp.getHeight() * scale );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        g.setRenderingHint( java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON );
        g.scale( scale, scale );
        tp.paintPhylogeny( g, false, exporting, tp.getWidth(), tp.getHeight(), 0, 0 );
        g.dispose();
        return img;
    }

    /** Pixels sitting between two rows of the SAME cell colour without being it: only a gap can do that. */
    private static int seamPixels( final BufferedImage img ) {
        int n = 0;
        for( int y = 2; y < ( img.getHeight() - 2 ); ++y ) {
            for( int x = 0; x < img.getWidth(); x += 2 ) {
                final int above = img.getRGB( x, y - 2 );
                if ( ( above == img.getRGB( x, y + 2 ) ) && isCellColor( above ) && ( img.getRGB( x, y ) != above ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** A saturated cell colour -- not the white background, not black text or branches. */
    private static boolean isCellColor( final int rgb ) {
        final int r = ( rgb >> 16 ) & 0xFF;
        final int g = ( rgb >> 8 ) & 0xFF;
        final int b = rgb & 0xFF;
        return ( Math.max( r, Math.max( g, b ) ) > 60 ) && ( ( r + g + b ) < 700 )
                && ( ( Math.abs( r - g ) + Math.abs( g - b ) + Math.abs( r - b ) ) > 40 );
    }

    /**
     * The vertical clustergram (root-top / root-bottom) draws it too -- the painter is wired into that branch, and
     * every geometry input there is different: the band reserve maps onto the other screen axis, the lines ride the
     * rotation R, and the headers are anchored upright rather than rotated. Without this the two orientations the
     * code deliberately paints were the only ones nothing checked.
     */
    private static void vertical( final boolean[] ok, final TreePanel tp ) {
        for( final Options.TREE_ORIENTATION o : new Options.TREE_ORIENTATION[] { Options.TREE_ORIENTATION.ROOT_TOP,
                Options.TREE_ORIENTATION.ROOT_BOTTOM } ) {
            tp.setTreeOrientation( o );
            tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
            paint( tp );
            final Rectangle logical = tp.matrixDendrogramBandForTest();
            if ( ( logical == null ) || ( tp.matrixColumnDendrogram() == null ) ) {
                fail( ok, o + ": the clustergram must still have a dendrogram" );
                continue;
            }
            final int before = tp.matrixDendrogramPaintsForTest();
            final BufferedImage img = render( tp );
            if ( tp.matrixDendrogramPaintsForTest() <= before ) {
                fail( ok, o + ": the dendrogram painter was never reached" );
                continue;
            }
            // the hook reports LOGICAL coordinates; in a vertical orientation the picture is a quarter turn away,
            // so the band has to be mapped through the same transform the paint rides
            final Rectangle device = toDevice( tp, logical );
            final int drawn = inkIn( img, device, screenBg( tp ) );
            if ( drawn < 50 ) {
                fail( ok, o + ": only " + drawn + " marked pixels in the dendrogram band " + device );
            }
        }
        tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
        paint( tp );
    }

    /** A logical rectangle's device bounding box (a quarter turn keeps it axis-aligned). */
    private static Rectangle toDevice( final TreePanel tp, final Rectangle r ) {
        final java.awt.geom.Point2D.Double a = tp.screenPoint( r.x, r.y );
        final java.awt.geom.Point2D.Double b = tp.screenPoint( r.x + r.width, r.y + r.height );
        final int x = ( int ) Math.floor( Math.min( a.x, b.x ) );
        final int y = ( int ) Math.floor( Math.min( a.y, b.y ) );
        return new Rectangle( x, y, ( int ) Math.ceil( Math.abs( b.x - a.x ) ),
                              ( int ) Math.ceil( Math.abs( b.y - a.y ) ) );
    }

    /** Circular is an approved exception; unrooted draws no columns at all. */
    private static void notRadial( final boolean[] ok, final TreePanel tp ) {
        tp.setMatrixColumnOrder( MatrixColumnOrder.Mode.CLUSTERED );
        for( final PHYLOGENY_GRAPHICS_TYPE t : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
            tp.setPhylogenyGraphicsType( t );
            paint( tp );
            if ( ( tp.matrixColumnDendrogram() != null ) || ( tp.matrixDendrogramBandHeight() != 0 ) ) {
                fail( ok, t + " must draw no column dendrogram and reserve nothing for one, got band "
                        + tp.matrixDendrogramBandHeight() );
            }
        }
        tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        paint( tp );
    }

    // ---- measuring -----------------------------------------------------------------------------------------------

    private static float firstTipY( final TreePanel tp ) {
        float min = Float.MAX_VALUE;
        for( final org.forester.phylogeny.PhylogenyNode t : tp.getPhylogeny().getExternalNodes() ) {
            min = Math.min( min, t.getYcoord() );
        }
        return min;
    }

    /** The band in DEVICE coordinates (the hook reports logical; root-left is the same space). */
    private static Rectangle bandOnScreen( final TreePanel tp ) {
        paint( tp );
        return tp.matrixDendrogramBandForTest();
    }

    private static void paint( final TreePanel tp ) {
        tp.getMainPanel().getControlPanel().showWhole();
        tp.validate();
        render( tp );
    }

    private static BufferedImage render( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( Math.max( tp.getWidth(), 100 ), Math.max( tp.getHeight(), 100 ),
                                                     BufferedImage.TYPE_INT_RGB );
        tp.printAll( img.getGraphics() );
        return img;
    }

    /** The same panel through the EXPORT path, so the figure -- not the screen -- is what is measured. */
    private static BufferedImage exportRender( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( Math.max( tp.getWidth(), 100 ), Math.max( tp.getHeight(), 100 ),
                                                     BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, img.getWidth(), img.getHeight() );
        tp.paintPhylogeny( g, false, true, img.getWidth(), img.getHeight(), 0, 0 );
        g.dispose();
        return img;
    }

    /**
     * Pixels in {@code r} that differ from {@code bg}. The background is passed IN rather than read from a corner:
     * the overview is on by default and placed upper-left, so pixel (0,0) is not reliably canvas background, and a
     * day when it is not, every empty pixel counts as ink and the positive checks below pass unconditionally.
     */
    private static int inkIn( final BufferedImage img, final Rectangle r, final int bg ) {
        int n = 0;
        for( int x = Math.max( 0, r.x ); x < Math.min( img.getWidth(), r.x + r.width ); ++x ) {
            for( int y = Math.max( 0, r.y ); y < Math.min( img.getHeight(), r.y + r.height ); ++y ) {
                if ( img.getRGB( x, y ) != bg ) {
                    ++n;
                }
            }
        }
        return n;
    }
}
