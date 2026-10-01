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
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * Heat-map cell borders (Christian, 2026-09-30, as Archaeopteryx.js a72365b): a HEATMAP or MATRIX cell whose smaller
 * side is at least 5 panel pixels, on a heat map of at most 60,000 cells, is outlined 0.75 px wide in its own colour
 * darkened as d3's darker(0.8) does it. Smaller cells, bigger heat maps and categorical colour strips are drawn as
 * before.
 * <p>
 * The colour is pinned to values computed by d3 ITSELF (the d3.v7.min.js Archaeopteryx.js ships, run in Node), not
 * to this implementation's formula restated. The renders measure BORDER INK: pixels on the line from a cell's colour
 * to its border colour (antialiasing toward the background or toward a differently coloured neighbour runs in another
 * direction), so equal neighbours that used to merge into a bar show the edge between them. Headful; the render parts
 * are a green no-op when headless.
 */
public final class HeatmapCellBorderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "HeatmapCellBorder: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( !borderColorOk() || !thresholdOk() ) {
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return rendersOk() && sizeCapOk();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [HeatmapCellBorderTest] " + msg );
        return false;
    }

    // ---- the colour: d3's own numbers ---------------------------------------------------------------------------
    /** {fill, d3.color(fill).darker(0.8).formatHex()}, computed with Archaeopteryx.js's d3.v7.min.js in Node. */
    private static final String[][] D3_DARKER_08 = { { "#440154", "#33013f" }, { "#fde725", "#beae1c" },
            { "#21918c", "#196d69" }, { "#3b528b", "#2c3e68" }, { "#5ec962", "#47974a" }, { "#c8643c", "#964b2d" },
            { "#ffffff", "#c0c0c0" }, { "#000000", "#000000" }, { "#010203", "#010202" } };

    private static boolean borderColorOk() {
        for( final String[] p : D3_DARKER_08 ) {
            final Color got = TreePanel.heatmapBorderColor( Color.decode( p[ 0 ] ) );
            final String hex = String.format( "#%02x%02x%02x", got.getRed(), got.getGreen(), got.getBlue() );
            if ( !p[ 1 ].equals( hex ) ) {
                return fail( "border of " + p[ 0 ] + " must be d3's " + p[ 1 ] + ", got " + hex );
            }
        }
        if ( TreePanel.heatmapBorderColor( new Color( 253, 231, 37, 128 ) ).getAlpha() != 128 ) {
            return fail( "the border keeps the fill's alpha" );
        }
        return true;
    }

    private static boolean thresholdOk() {
        if ( !TreePanel.heatmapCellsBordered( 5, 60000 ) ) {
            return fail( "a 5 px cell on a 60,000-cell heat map is bordered (both limits are inclusive)" );
        }
        if ( TreePanel.heatmapCellsBordered( 4.99, 10 ) ) {
            return fail( "a cell under 5 px is not bordered" );
        }
        if ( TreePanel.heatmapCellsBordered( 100, 60001 ) ) {
            return fail( "a heat map of more than 60,000 cells is not bordered, however big its cells" );
        }
        return true;
    }

    // ---- renders --------------------------------------------------------------------------------------------------
    private static final String[] LAYOUTS = { "ROOT_LEFT", "ROOT_TOP", "ROOT_BOTTOM", "CIRCULAR" };

    /** Border ink an UNBORDERED render may still show: antialiasing that happens to lie on the line (measured: at
     *  most 39 px, circular; a bordered render showed at least 794). */
    private static final int UNBORDERED_MAX = 150;

    /** Viridis' top end: the colour of the fixture's value 9, on a 0..9 scale. */
    private static final Color HIGH = Color.decode( "#fde725" );

    private static boolean rendersOk() {
        final boolean[] ok = { true };
        // BIG cells: 12 tips; s1 is 9 on EVERY tip -- a column of equal neighbours, the case borders are for
        runIn( tree( 12, 4 ), AnnotationColumns.Type.MATRIX, 4, ( frame, tp ) -> {
            for( final String layout : LAYOUTS ) {
                final BufferedImage img = render( frame, tp, layout, 900, 900 );
                final int ink = borderInk( img, HIGH );
                if ( !tp.heatmapBorderedForTest() ) {
                    ok[ 0 ] = fail( layout + ": cells this big must be bordered" );
                }
                if ( ink < 300 ) {
                    ok[ 0 ] = fail( layout + ": bordered cells must show border ink between equal neighbours, got "
                            + ink + " px" );
                }
                // the count the cap is weighed against is tips x heat-map columns in EVERY layout (the cap's own
                // render pair below runs rectangular only: 61 columns of 1,000 tips is too big an image elsewhere)
                if ( tp.heatmapCellCountForTest() != ( 12 * 4 ) ) {
                    ok[ 0 ] = fail( layout + ": the cell count must be 12 tips x 4 columns = 48, got "
                            + tp.heatmapCellCountForTest() );
                }
                // THIN: a 0.75 px border centred on an edge touches at most 2 pixels across it
                if ( !"CIRCULAR".equals( layout ) ) {
                    final int[] runs = inkRunsBetweenEqualCells( img, HIGH );
                    if ( runs[ 0 ] < 5 ) {
                        ok[ 0 ] = fail( layout + ": fixture: the scan must cross at least 5 edges between equal"
                                + " neighbours, crossed " + runs[ 0 ] );
                    }
                    else if ( runs[ 1 ] > 2 ) {
                        ok[ 0 ] = fail( layout + ": a border is " + runs[ 1 ] + " px thick across an edge; 0.75 px"
                                + " touches at most 2" );
                    }
                }
            }
        } );
        // SMALL cells: 400 tips in the same window -- rows and arcs under 5 px, drawn as before
        runIn( tree( 400, 4 ), AnnotationColumns.Type.MATRIX, 4, ( frame, tp ) -> {
            for( final String layout : LAYOUTS ) {
                final BufferedImage img = render( frame, tp, layout, 900, 900 );
                final int ink = borderInk( img, HIGH );
                if ( tp.heatmapBorderedForTest() ) {
                    ok[ 0 ] = fail( layout + ": cells under 5 px must not be bordered" );
                }
                if ( ink > UNBORDERED_MAX ) {
                    ok[ 0 ] = fail( layout + ": unbordered cells must show no border ink, got " + ink + " px" );
                }
            }
        } );
        // a categorical COLOR_STRIP is not a heat map: big cells, still no border, in any layout
        runIn( tree( 12, 4 ), AnnotationColumns.Type.COLOR_STRIP, 4, ( frame, tp ) -> {
            for( final String layout : LAYOUTS ) {
                final BufferedImage img = render( frame, tp, layout, 900, 900 );
                final Color strip = commonestBrightColor( img );
                if ( strip == null ) {
                    ok[ 0 ] = fail( layout + ": fixture: the colour strip drew no cells" );
                    continue;
                }
                // the commonest strip colour (s1: one value on every tip, equal neighbours) AND Viridis' top: a strip
                // of numbers is drawn on the gradient, so its 9s are the same yellow a heat map's are
                for( final Color fill : new Color[] { strip, HIGH } ) {
                    final int ink = borderInk( img, fill );
                    if ( ink > UNBORDERED_MAX ) {
                        ok[ 0 ] = fail( layout + ": a colour strip must not be bordered, got " + ink
                                + " px of border ink around " + fill );
                    }
                }
            }
        } );
        return ok[ 0 ];
    }

    /** The 60,000-cell cap, as a pair that differs only in one column: 1,000 tips x 60 is bordered, x 61 is not. Rows
     *  are 7 px in both, so the cell size cannot be the reason. */
    private static boolean sizeCapOk() {
        final boolean[] ok = { true };
        for( final int cols : new int[] { 60, 61 } ) {
            runIn( tree( 1000, cols ), AnnotationColumns.Type.MATRIX, cols, ( frame, tp ) -> {
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                final int w = 1400, h = 7300;
                frame.showWhole();
                tp.setSize( w, h );
                tp.calcParametersForPainting( w, h );
                final BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, frame.getOptions(), false, 1,
                                                                           false );
                final double row = 2.0 * tp.getYdistance();
                if ( row < 5 ) {
                    ok[ 0 ] = fail( "fixture: rows must be at least 5 px so only the count decides, got " + row );
                }
                final boolean want = ( cols == 60 );
                if ( tp.heatmapBorderedForTest() != want ) {
                    ok[ 0 ] = fail( "1,000 tips x " + cols + " heat-map columns must " + ( want ? "" : "NOT " )
                            + "be bordered" );
                }
                final int ink = borderInk( img, HIGH );
                if ( want ? ( ink < 1000 ) : ( ink > UNBORDERED_MAX ) ) {
                    ok[ 0 ] = fail( "1,000 x " + cols + ": border ink " + ink + " px does not match bordered=" + want );
                }
            } );
        }
        return ok[ 0 ];
    }

    // ---- helpers --------------------------------------------------------------------------------------------------
    interface PanelCheck {

        void check( MainFrame frame, TreePanel tp ) throws Exception;
    }

    private static void runIn( final Phylogeny phy, final AnnotationColumns.Type type, final int cols,
                               final PanelCheck check ) {
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "borders" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    tp.setColorByPropertyRef( null );
                    mf[ 0 ].getOptions().setGraphicsExportWhiteBackground( true );
                    final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<>();
                    for( int j = 1; j <= cols; ++j ) {
                        specs.add( new AnnotationColumns.ColumnSpec( "data:s" + j, type ) );
                    }
                    tp.setAnnotationColumns( specs );
                    check.check( mf[ 0 ], tp );
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    fail( "unexpected: " + t );
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            fail( "unexpected: " + e );
        }
    }

    private static BufferedImage render( final MainFrame frame, final TreePanel tp, final String layout, final int w,
                                         final int h ) {
        if ( "CIRCULAR".equals( layout ) ) {
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        }
        else {
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp.setTreeOrientation( Options.TREE_ORIENTATION.valueOf( layout ) );
        }
        frame.showWhole();
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        return AptxUtil.renderPhylogenyToImage( w, h, tp, frame.getOptions(), false, 1, false );
    }

    /** A star-ish tree of {@code tips} tips (pairs under one root), each with s1..s{cols}: s1 = 9 on every tip (equal
     *  neighbours), the other columns 9 or 0 in a pattern, so the scale runs 0..9 and 9 is Viridis' top. */
    private static Phylogeny tree( final int tips, final int cols ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        PhylogenyNode pair = null;
        for( int i = 0; i < tips; ++i ) {
            if ( ( i % 2 ) == 0 ) {
                pair = new PhylogenyNode();
                pair.setDistanceToParent( 1 );
                root.addAsChild( pair );
            }
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "t" + i );
            n.setDistanceToParent( 1 );
            final PropertiesList pl = new PropertiesList();
            for( int j = 1; j <= cols; ++j ) {
                final String v = ( ( j == 1 ) || ( ( ( i + j ) % 3 ) == 0 ) ) ? "9" : "0";
                pl.addProperty( new Property( "data:s" + j, v, "", "xsd:decimal", AppliesTo.NODE ) );
            }
            n.getNodeData().setProperties( pl );
            pair.addAsChild( n );
        }
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /**
     * Pixels on the line from {@code fill} to its border colour: p = fill + t (border - fill) for one t in
     * [0.2, 1.05], every channel within 3 of it. Antialiasing toward white, toward another cell colour or toward a
     * label runs in another direction and is not counted.
     */
    static int borderInk( final BufferedImage img, final Color fill ) {
        final Color b = TreePanel.heatmapBorderColor( fill );
        final int[] c = { fill.getRed(), fill.getGreen(), fill.getBlue() };
        final int[] d = { b.getRed() - c[ 0 ], b.getGreen() - c[ 1 ], b.getBlue() - c[ 2 ] };
        int k = 0;
        for( int i = 1; i < 3; ++i ) {
            if ( Math.abs( d[ i ] ) > Math.abs( d[ k ] ) ) {
                k = i;
            }
        }
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int[] p = { ( rgb >> 16 ) & 0xFF, ( rgb >> 8 ) & 0xFF, rgb & 0xFF };
                final double t = ( p[ k ] - c[ k ] ) / (double) d[ k ];
                if ( ( t < 0.2 ) || ( t > 1.05 ) ) {
                    continue;
                }
                boolean on = true;
                for( int i = 0; i < 3; ++i ) {
                    if ( Math.abs( p[ i ] - ( c[ i ] + ( t * d[ i ] ) ) ) > 3 ) {
                        on = false;
                        break;
                    }
                }
                if ( on ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /**
     * Scans across the column (or, in a vertical orientation, the row) of the image that holds the most pixels of
     * exactly {@code fill}, and measures every run of border ink between two stretches of {@code fill}: an edge
     * between EQUAL neighbours. Returns {runs crossed, thickest run}.
     */
    private static int[] inkRunsBetweenEqualCells( final BufferedImage img, final Color fill ) {
        final int f = fill.getRGB() & 0xFFFFFF;
        int best_line = -1, best_n = -1;
        boolean best_is_column = true;
        for( final boolean column : new boolean[] { true, false } ) {
            final int lines = column ? img.getWidth() : img.getHeight();
            final int len = column ? img.getHeight() : img.getWidth();
            for( int l = 0; l < lines; ++l ) {
                int n = 0;
                for( int i = 0; i < len; ++i ) {
                    if ( ( pixel( img, column, l, i ) & 0xFFFFFF ) == f ) {
                        ++n;
                    }
                }
                if ( n > best_n ) {
                    best_n = n;
                    best_line = l;
                    best_is_column = column;
                }
            }
        }
        final int len = best_is_column ? img.getHeight() : img.getWidth();
        int runs = 0, thickest = 0, run = 0;
        boolean after_fill = false;
        for( int i = 0; i < len; ++i ) {
            final int rgb = pixel( img, best_is_column, best_line, i ) & 0xFFFFFF;
            if ( rgb == f ) {
                if ( after_fill && ( run > 0 ) ) {
                    ++runs;
                    thickest = Math.max( thickest, run );
                }
                after_fill = true;
                run = 0;
            }
            else if ( after_fill && onBorderLine( rgb, fill ) ) {
                ++run;
            }
            else {
                after_fill = false;
                run = 0;
            }
        }
        return new int[] { runs, thickest };
    }

    private static int pixel( final BufferedImage img, final boolean column, final int line, final int i ) {
        return column ? img.getRGB( line, i ) : img.getRGB( i, line );
    }

    /** Whether {@code rgb} lies on the line from {@code fill} to its border colour (the test borderInk counts). */
    private static boolean onBorderLine( final int rgb, final Color fill ) {
        final BufferedImage one = new BufferedImage( 1, 1, BufferedImage.TYPE_INT_RGB );
        one.setRGB( 0, 0, rgb );
        return borderInk( one, fill ) == 1;
    }

    /** The commonest colour in the image whose brightest channel is over 150 and which is not grey -- a strip cell. */
    private static Color commonestBrightColor( final BufferedImage img ) {
        final Map<Integer, Integer> counts = new HashMap<>();
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y ) & 0xFFFFFF;
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int max = Math.max( r, Math.max( g, b ) ), min = Math.min( r, Math.min( g, b ) );
                if ( ( max > 150 ) && ( ( max - min ) > 60 ) ) {
                    counts.merge( rgb, 1, Integer::sum );
                }
            }
        }
        int best = -1, bestN = 0;
        for( final Map.Entry<Integer, Integer> e : counts.entrySet() ) {
            if ( e.getValue() > bestN ) {
                best = e.getKey();
                bestN = e.getValue();
            }
        }
        return ( best < 0 ) ? null : new Color( best );
    }

    private HeatmapCellBorderTest() {
    }
}
