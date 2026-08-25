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
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;

/**
 * Integration test for the "Annotation Columns" tool: on a real {@link MainFrame}/{@link TreePanel} it sets
 * annotation columns of every render type, renders offscreen (exercising the paint + width-reservation path),
 * and checks the menu item is present. Guarded to a no-op on a headless box (needs FlatLaf via
 * {@code createInstance}); run standalone or as part of the non-headless suite.
 */
public final class AnnotationColumnsToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "AnnotationColumnsTool: " + ( ok ? "OK." : "FAILED." ) );
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
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { annotatedTree() }, conf,
                                                                        "anncols" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainPanel mp = mf[ 0 ].getMainPanel();
                final TreePanel tp = mp.getCurrentTreePanel();
                // pin the rectangular layout up front so the header-reserve / header-click assertions below are
                // independent of the developer's persisted ~/.archaeopteryx graphics type (a standalone run inherits
                // it; the isolated-cache suite does not) -- the documented standalone-vs-suite gotcha
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );

                // menu item present + tooltip
                final JMenuItem item = toolsItem( mf[ 0 ].getJMenuBar(), "Annotation Columns" );
                if ( item == null ) {
                    System.out.println( "  [AnnotationColumnsToolTest] Tools menu item not found" );
                    ok[ 0 ] = false;
                }
                else if ( ( item.getToolTipText() == null )
                        || !item.getToolTipText().toLowerCase().contains( "column" ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] menu item missing tooltip" );
                    ok[ 0 ] = false;
                }

                // baseline render (bare tree, no columns) for a pixel-count comparison
                mf[ 0 ].showWhole();
                final int w = 700, h = 360;
                final int without = renderOffscreen( tp, w, h );

                // set one column of every render type
                final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<AnnotationColumns.ColumnSpec>();
                specs.add( new AnnotationColumns.ColumnSpec( "data:region", AnnotationColumns.Type.COLOR_STRIP ) );
                specs.add( new AnnotationColumns.ColumnSpec( "data:value", AnnotationColumns.Type.HEATMAP ) );
                specs.add( new AnnotationColumns.ColumnSpec( "data:value", AnnotationColumns.Type.BAR ) );
                specs.add( new AnnotationColumns.ColumnSpec( "data:region", AnnotationColumns.Type.TEXT ) );
                tp.setAnnotationColumns( specs );
                if ( !tp.hasAnnotationColumns() ) {
                    System.out.println( "  [AnnotationColumnsToolTest] columns were not set" );
                    ok[ 0 ] = false;
                }
                final float first_tip_y_without = tp.getPhylogeny().getFirstExternalNode().getYcoord();
                mf[ 0 ].showWhole();
                // render with the columns -- exercises paintAnnotationColumns + the width reservation; the
                // columns must add painted pixels beyond what the bare tree draws (a whole-right-half scan
                // would false-pass on the tree's own aligned tips/labels)
                final int with = renderOffscreen( tp, w, h );
                if ( with <= without ) {
                    System.out.println( "  [AnnotationColumnsToolTest] annotation columns drew no additional pixels" );
                    ok[ 0 ] = false;
                }

                // every column reserves a positive width (exercises the cached width path)
                for( int i = 0; i < specs.size(); ++i ) {
                    if ( tp.annotationColumnWidthForTest( i ) <= 0 ) {
                        System.out.println( "  [AnnotationColumnsToolTest] column " + i + " has non-positive width" );
                        ok[ 0 ] = false;
                    }
                }
                // the rotated headers reserve top space, and that space pushes the first tip down so the
                // headers don't overlap the top cells (item: header/cell overlap on tight trees)
                if ( tp.annotationHeaderTopReserveForTest() <= 0 ) {
                    System.out.println( "  [AnnotationColumnsToolTest] no header top-reserve was computed" );
                    ok[ 0 ] = false;
                }
                if ( tp.getPhylogeny().getFirstExternalNode().getYcoord() <= first_tip_y_without ) {
                    System.out.println( "  [AnnotationColumnsToolTest] header reserve did not push the first tip down" );
                    ok[ 0 ] = false;
                }
                // the preferred (scrollable) height must cover the reserve-shifted tree, or the bottom cells
                // clip when zoomed/scrolled -- resetPreferredSize must include the same reserve the paint shifts by
                tp.resetPreferredSize();
                float last_tip_y = 0f;
                for( final PhylogenyNode t : tp.getPhylogeny().getExternalNodes() ) {
                    last_tip_y = Math.max( last_tip_y, t.getYcoord() );
                }
                final int last_cell_bottom = Math.round( last_tip_y + tp.getYdistance() );
                if ( tp.getPreferredSize().height < last_cell_bottom ) {
                    System.out.println( "  [AnnotationColumnsToolTest] preferred height " + tp.getPreferredSize().height
                            + " clips the bottom annotation cells (need >= " + last_cell_bottom + ")" );
                    ok[ 0 ] = false;
                }

                // clicking a color-column header toggles its legend -- exercises the real hit-test path
                final java.awt.Point hp = tp.annotationHeaderMidpointForTest( 0 ); // COLOR_STRIP
                if ( hp == null ) {
                    System.out.println( "  [AnnotationColumnsToolTest] no header hit-point for column 0" );
                    ok[ 0 ] = false;
                }
                else {
                    final MouseEvent click = new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(),
                                                             0, hp.x, hp.y, 1, false );
                    if ( !tp.handleAnnotationHeaderClick( click ) || !tp.hasFocusedAnnotationColumn() ) {
                        System.out.println( "  [AnnotationColumnsToolTest] header click did not open the legend" );
                        ok[ 0 ] = false;
                    }
                    tp.handleAnnotationHeaderClick( click ); // a second click toggles it off
                    if ( tp.hasFocusedAnnotationColumn() ) {
                        System.out.println( "  [AnnotationColumnsToolTest] second header click did not close the legend" );
                        ok[ 0 ] = false;
                    }
                }

                // focusing a color-coded column shows that column's legend (bounds recorded on render)
                tp.setFocusedAnnotationColumn( 0 ); // COLOR_STRIP
                renderOffscreen( tp, w, h );
                if ( !tp.hasFocusedAnnotationColumn() || ( tp.getPropertyLegendBounds() == null ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] focused strip legend not drawn" );
                    ok[ 0 ] = false;
                }
                tp.setFocusedAnnotationColumn( 1 ); // HEATMAP
                renderOffscreen( tp, w, h );
                if ( !tp.hasFocusedAnnotationColumn() || ( tp.getPropertyLegendBounds() == null ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] focused heat-map legend not drawn" );
                    ok[ 0 ] = false;
                }
                tp.setFocusedAnnotationColumn( 2 ); // BAR -> a length/range legend (not a color key, but a legend)
                renderOffscreen( tp, w, h );
                if ( !tp.hasFocusedAnnotationColumn() || ( tp.getPropertyLegendBounds() == null ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] focused bar legend not drawn" );
                    ok[ 0 ] = false;
                }
                tp.setFocusedAnnotationColumn( 3 ); // TEXT -> no legend
                if ( tp.hasFocusedAnnotationColumn() ) {
                    System.out.println( "  [AnnotationColumnsToolTest] a text column should not focus a legend" );
                    ok[ 0 ] = false;
                }
                tp.setFocusedAnnotationColumn( 0 );
                tp.setFocusedAnnotationColumn( 0 ); // toggling the same column off
                if ( tp.hasFocusedAnnotationColumn() ) {
                    System.out.println( "  [AnnotationColumnsToolTest] a second header click should toggle the legend off" );
                    ok[ 0 ] = false;
                }

                // VERTICAL PARITY: the tip-aligned columns become horizontal BANDS below the tips. Re-fit in a
                // vertical orientation and verify (a) the columns paint (more ink than the bare vertical tree),
                // (b) the header hit-test works, and (c) the focused-column legend still draws.
                tp.setFocusedAnnotationColumn( -1 );
                tp.getOptions().setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                tp.getControlPanel().updateZoomButtonsForLayout();
                // count SATURATED (colored) pixels: the tree/labels/bars/text are grayscale, so only the strip +
                // heat-map CELLS produce colored pixels -- a check that is not confounded by the layout reservation
                // (which compresses the tree even when the columns do not paint)
                tp.clearAnnotationColumns();
                tp.setSize( w, h );
                tp.calcParametersForPainting( w, h ); // fit to the render size so the bands are within the image
                tp.resetPreferredSize();
                final int v_colored_without = coloredPixels( tp, w, h );
                tp.setAnnotationColumns( specs );
                tp.setSize( w, h );
                tp.calcParametersForPainting( w, h );
                tp.resetPreferredSize();
                final int v_colored_with = coloredPixels( tp, w, h );
                if ( v_colored_with <= ( v_colored_without + 50 ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] vertical annotation cells drew no colored "
                            + "pixels (" + v_colored_without + " -> " + v_colored_with + ")" );
                    ok[ 0 ] = false;
                }
                final java.awt.Point vhp = tp.annotationHeaderMidpointForTest( 0 ); // COLOR_STRIP header, in device px
                if ( vhp == null ) {
                    System.out.println( "  [AnnotationColumnsToolTest] no vertical header hit-point for column 0" );
                    ok[ 0 ] = false;
                }
                else {
                    final MouseEvent vclick = new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(),
                                                              0, vhp.x, vhp.y, 1, false );
                    if ( !tp.handleAnnotationHeaderClick( vclick ) || !tp.hasFocusedAnnotationColumn() ) {
                        System.out.println( "  [AnnotationColumnsToolTest] vertical header click did not open the legend" );
                        ok[ 0 ] = false;
                    }
                    tp.handleAnnotationHeaderClick( vclick ); // second click toggles it off
                    if ( tp.hasFocusedAnnotationColumn() ) {
                        System.out.println( "  [AnnotationColumnsToolTest] vertical second header click did not close it" );
                        ok[ 0 ] = false;
                    }
                }
                tp.setFocusedAnnotationColumn( 0 ); // COLOR_STRIP -> a color key
                renderOffscreen( tp, w, h );
                if ( !tp.hasFocusedAnnotationColumn() || ( tp.getPropertyLegendBounds() == null ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] vertical focused legend not drawn" );
                    ok[ 0 ] = false;
                }
                // WYSIWYG: exports go through the same paintPhylogeny funnel, so the vertical bands must appear in an
                // EXPORT render too (not just the screen path)
                tp.setFocusedAnnotationColumn( -1 );
                tp.setSize( w, h );
                tp.calcParametersForPainting( w, h );
                tp.resetPreferredSize();
                final BufferedImage exp = AptxUtil.renderPhylogenyToImage( w, h, tp, tp.getOptions(), false, 1, false );
                if ( countSaturated( exp ) <= 50 ) {
                    System.out.println( "  [AnnotationColumnsToolTest] vertical annotation cells absent from the export" );
                    ok[ 0 ] = false;
                }
                tp.getOptions().setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT ); // restore for the rest
                mf[ 0 ].showWhole();

                // STACKED_BAR: three numeric fields MERGE into ONE segmented bar (each a distinctly-coloured series).
                // The tree/labels are grayscale, so the colored-pixel count isolates the vivid stack segments.
                tp.setFocusedAnnotationColumn( -1 );
                tp.clearAnnotationColumns();
                mf[ 0 ].showWhole();
                final int stack_bare = coloredPixels( tp, w, h );
                final List<AnnotationColumns.ColumnSpec> stacked = new ArrayList<AnnotationColumns.ColumnSpec>();
                stacked.add( new AnnotationColumns.ColumnSpec( "data:seg_a", AnnotationColumns.Type.STACKED_BAR ) );
                stacked.add( new AnnotationColumns.ColumnSpec( "data:seg_b", AnnotationColumns.Type.STACKED_BAR ) );
                stacked.add( new AnnotationColumns.ColumnSpec( "data:seg_c", AnnotationColumns.Type.STACKED_BAR ) );
                tp.setAnnotationColumns( stacked );
                if ( tp.annotationColumnCountForTest() != 1 ) {
                    System.out.println( "  [AnnotationColumnsToolTest] three STACKED_BAR fields should MERGE into 1 "
                            + "column, got " + tp.annotationColumnCountForTest() );
                    ok[ 0 ] = false;
                }
                mf[ 0 ].showWhole();
                final int stack_with = coloredPixels( tp, w, h );
                if ( stack_with <= ( stack_bare + 50 ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] stacked-bar segments drew no colored pixels ("
                            + stack_bare + " -> " + stack_with + ")" );
                    ok[ 0 ] = false;
                }
                // the merged column shows a series-colour legend when focused
                tp.setFocusedAnnotationColumn( 0 );
                renderOffscreen( tp, w, h );
                if ( !tp.hasFocusedAnnotationColumn() || ( tp.getPropertyLegendBounds() == null ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] focused stacked-bar legend not drawn" );
                    ok[ 0 ] = false;
                }
                tp.setFocusedAnnotationColumn( -1 );
                // a NORMALIZED stacked column also renders segments (the fraction math is covered in AnnotationColumnsTest)
                final List<AnnotationColumns.ColumnSpec> norm = new ArrayList<AnnotationColumns.ColumnSpec>();
                norm.add( new AnnotationColumns.ColumnSpec( "data:seg_a", AnnotationColumns.Type.STACKED_BAR, true ) );
                norm.add( new AnnotationColumns.ColumnSpec( "data:seg_b", AnnotationColumns.Type.STACKED_BAR, true ) );
                norm.add( new AnnotationColumns.ColumnSpec( "data:seg_c", AnnotationColumns.Type.STACKED_BAR, true ) );
                tp.setAnnotationColumns( norm );
                mf[ 0 ].showWhole();
                if ( coloredPixels( tp, w, h ) <= ( stack_bare + 50 ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] normalized stacked-bar drew no colored pixels" );
                    ok[ 0 ] = false;
                }
                // CIRCULAR parity: the segments become stacked radial arcs (vivid, distinct from the grayscale disc)
                tp.getOptions().setGraphicsExportWhiteBackground( true );
                tp.getOptions().setShowOverview( false );
                tp.setOvOn( false );
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                final int cw = 700, chh = 700;
                tp.setPreferredSize( new java.awt.Dimension( cw, chh ) );
                tp.setSize( cw, chh );
                tp.clearAnnotationColumns();
                tp.calcParametersForPainting( cw, chh );
                final int circ_bare = countSaturated(
                        AptxUtil.renderPhylogenyToImage( cw, chh, tp, tp.getOptions(), false, 1, false ) );
                tp.setAnnotationColumns( stacked );
                tp.calcParametersForPainting( cw, chh );
                final int circ_with = countSaturated(
                        AptxUtil.renderPhylogenyToImage( cw, chh, tp, tp.getOptions(), false, 1, false ) );
                if ( circ_with <= ( circ_bare + 50 ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] circular stacked-bar arcs drew no colored pixels ("
                            + circ_bare + " -> " + circ_with + ")" );
                    ok[ 0 ] = false;
                }
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.clearAnnotationColumns();
                mf[ 0 ].showWhole();

                // PIE: the same numeric fields as WEDGES of one pie glyph per tip (merge + rectangular + circular)
                final List<AnnotationColumns.ColumnSpec> pies = new ArrayList<AnnotationColumns.ColumnSpec>();
                pies.add( new AnnotationColumns.ColumnSpec( "data:seg_a", AnnotationColumns.Type.PIE ) );
                pies.add( new AnnotationColumns.ColumnSpec( "data:seg_b", AnnotationColumns.Type.PIE ) );
                pies.add( new AnnotationColumns.ColumnSpec( "data:seg_c", AnnotationColumns.Type.PIE ) );
                tp.setAnnotationColumns( pies );
                if ( tp.annotationColumnCountForTest() != 1 ) {
                    System.out.println( "  [AnnotationColumnsToolTest] three PIE fields should MERGE into 1 column, got "
                            + tp.annotationColumnCountForTest() );
                    ok[ 0 ] = false;
                }
                mf[ 0 ].showWhole();
                if ( coloredPixels( tp, w, h ) <= ( stack_bare + 50 ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] pie wedges drew no colored pixels" );
                    ok[ 0 ] = false;
                }
                tp.setFocusedAnnotationColumn( 0 );
                renderOffscreen( tp, w, h );
                if ( !tp.hasFocusedAnnotationColumn() || ( tp.getPropertyLegendBounds() == null ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] focused pie legend not drawn" );
                    ok[ 0 ] = false;
                }
                tp.setFocusedAnnotationColumn( -1 );
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                tp.setPreferredSize( new java.awt.Dimension( cw, chh ) );
                tp.setSize( cw, chh );
                tp.clearAnnotationColumns();
                tp.calcParametersForPainting( cw, chh );
                final int pie_circ_bare = countSaturated(
                        AptxUtil.renderPhylogenyToImage( cw, chh, tp, tp.getOptions(), false, 1, false ) );
                tp.setAnnotationColumns( pies );
                tp.calcParametersForPainting( cw, chh );
                final int pie_circ_with = countSaturated(
                        AptxUtil.renderPhylogenyToImage( cw, chh, tp, tp.getOptions(), false, 1, false ) );
                if ( pie_circ_with <= ( pie_circ_bare + 50 ) ) {
                    System.out.println( "  [AnnotationColumnsToolTest] circular pies drew no colored pixels ("
                            + pie_circ_bare + " -> " + pie_circ_with + ")" );
                    ok[ 0 ] = false;
                }
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.clearAnnotationColumns();
                mf[ 0 ].showWhole();

                // a BAR column whose values are not a numeric gradient (categorical here; or a subtree that
                // collapsed to a single value) must NOT focus a legend -- its min/max range would be degenerate
                final List<AnnotationColumns.ColumnSpec> deg = new ArrayList<AnnotationColumns.ColumnSpec>();
                deg.add( new AnnotationColumns.ColumnSpec( "data:region", AnnotationColumns.Type.BAR ) );
                tp.setAnnotationColumns( deg );
                tp.setFocusedAnnotationColumn( 0 );
                if ( tp.hasFocusedAnnotationColumn() ) {
                    System.out.println( "  [AnnotationColumnsToolTest] a non-gradient bar column should not focus a legend" );
                    ok[ 0 ] = false;
                }

                // clearing removes them
                tp.clearAnnotationColumns();
                if ( tp.hasAnnotationColumns() ) {
                    System.out.println( "  [AnnotationColumnsToolTest] clear did not remove the columns" );
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

    /** Renders the panel offscreen and returns the number of non-white pixels drawn. */
    private static int renderOffscreen( final TreePanel tp, final int w, final int h ) {
        tp.setSize( w, h );
        tp.validate();
        tp.doLayout();
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.printAll( g );
        g.dispose();
        int n = 0;
        for( int x = 0; x < w; ++x ) {
            for( int y = 0; y < h; ++y ) {
                if ( ( img.getRGB( x, y ) & 0x00FFFFFF ) != 0x00FFFFFF ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Renders the panel offscreen (screen path) and counts SATURATED (colored, non-grayscale) pixels -- the
     *  tree/labels/bars/text are grayscale, so this isolates the strip/heat-map annotation CELLS. */
    private static int coloredPixels( final TreePanel tp, final int w, final int h ) {
        tp.setSize( w, h );
        tp.validate();
        tp.doLayout();
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.printAll( g );
        g.dispose();
        return countSaturated( img );
    }

    /** Colored (non-grayscale) pixels in an image -- the annotation strip/heat-map cells against a grayscale tree. */
    private static int countSaturated( final BufferedImage img ) {
        int n = 0;
        for( int x = 0; x < img.getWidth(); ++x ) {
            for( int y = 0; y < img.getHeight(); ++y ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, gc = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int mx = Math.max( r, Math.max( gc, b ) ), mn = Math.min( r, Math.min( gc, b ) );
                if ( ( mx - mn ) > 40 ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static JMenuItem toolsItem( final JMenuBar bar, final String contains ) {
        for( int i = 0; i < bar.getMenuCount(); ++i ) {
            final JMenu m = bar.getMenu( i );
            if ( ( m != null ) && "Tools".equals( m.getText() ) ) {
                for( int j = 0; j < m.getItemCount(); ++j ) {
                    final JMenuItem it = m.getItem( j );
                    if ( ( it != null ) && ( it.getText() != null ) && it.getText().contains( contains ) ) {
                        return it;
                    }
                }
            }
        }
        return null;
    }

    private static Phylogeny annotatedTree() {
        final PhylogenyNode root = new PhylogenyNode();
        final String[] regions = { "Asia", "Asia", "Europe", "Europe", "Africa", "Africa" };
        PhylogenyNode cur = root;
        for( int i = 0; i < 6; ++i ) {
            final PhylogenyNode in = new PhylogenyNode();
            cur.addAsChild( in );
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:region", regions[ i ], "", "xsd:string", AppliesTo.NODE ) );
            pl.addProperty( new Property( "data:value", String.valueOf( ( i + 1 ) * 5 ), "", "xsd:string",
                                          AppliesTo.NODE ) );
            // three compositional numeric fields (varying totals) for the STACKED_BAR checks
            pl.addProperty( new Property( "data:seg_a", String.valueOf( i + 1 ), "", "xsd:string", AppliesTo.NODE ) );
            pl.addProperty( new Property( "data:seg_b", String.valueOf( 6 - i ), "", "xsd:string", AppliesTo.NODE ) );
            pl.addProperty( new Property( "data:seg_c", String.valueOf( ( i % 3 ) + 1 ), "", "xsd:string",
                                          AppliesTo.NODE ) );
            leaf.getNodeData().setProperties( pl );
            in.addAsChild( leaf );
            cur = in;
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }
}
