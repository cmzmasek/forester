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
import java.io.File;
import java.nio.file.Files;
import java.util.Arrays;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The command-line figure renderer behind {@code aptx_render}.
 * <p>
 * The load-bearing property is DETERMINISM: a figure pipeline that produced a different result on a different
 * machine would be worthless, so the same command must yield byte-identical output, and it must not inherit
 * whatever the operator happens to have saved in the GUI. The format routing is the other thing worth pinning --
 * PDF has its own exporter, and handing it to the raster path silently writes nothing at all.
 * <p>
 * The pure parts run anywhere; the rendering parts need a display (see {@link FigureRenderer}) and are skipped
 * when headless.
 */
public final class AptxRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "AptxRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return pureOk() && readTreesOk() && ( GraphicsEnvironment.isHeadless() || ( renderOk() && legendColumnOk() ) );
    }

    /** Format routing and size parsing -- no display needed. */
    private static boolean pureOk() {
        if ( FigureRenderer.formatFor( new File( "a.pdf" ) ) != GraphicsExportType.PDF ) {
            return fail( ".pdf must route to the PDF exporter" );
        }
        if ( FigureRenderer.formatFor( new File( "a.SVG" ) ) != GraphicsExportType.SVG ) {
            return fail( "the extension test must be case-insensitive" );
        }
        if ( FigureRenderer.formatFor( new File( "a.eps" ) ) != GraphicsExportType.EPS ) {
            return fail( ".eps" );
        }
        if ( FigureRenderer.formatFor( new File( "a.png" ) ) != GraphicsExportType.PNG ) {
            return fail( ".png" );
        }
        if ( FigureRenderer.formatFor( new File( "a.jpeg" ) ) != GraphicsExportType.JPG ) {
            return fail( ".jpeg must be accepted as well as .jpg" );
        }
        if ( FigureRenderer.formatFor( new File( "a.tiff" ) ) != GraphicsExportType.TIFF ) {
            return fail( ".tiff" );
        }
        // anything we cannot write must be REFUSED, not guessed at
        if ( FigureRenderer.formatFor( new File( "a.txt" ) ) != null ) {
            return fail( "an unknown extension must not be guessed at" );
        }
        if ( FigureRenderer.formatFor( new File( "noextension" ) ) != null ) {
            return fail( "a name with no extension must not be guessed at" );
        }
        // sizes
        final ExportSizeSpec mm = FigureRenderer.parseSize( "170x120mm", 300 );
        if ( ( mm == null ) || ( mm.unit() != ExportSizeSpec.Unit.MILLIMETERS ) ) {
            return fail( "170x120mm must parse as millimetres" );
        }
        if ( FigureRenderer.parseSize( "8x6in", 300 ).unit() != ExportSizeSpec.Unit.INCHES ) {
            return fail( "8x6in must parse as inches" );
        }
        if ( FigureRenderer.parseSize( "1200x900px", 72 ).unit() != ExportSizeSpec.Unit.PIXELS ) {
            return fail( "1200x900px must parse as pixels" );
        }
        if ( FigureRenderer.parseSize( " 170X120MM ", 300 ) == null ) {
            return fail( "a size should tolerate case and surrounding space" );
        }
        for ( final String bad : new String[] { null, "", "170x120", "170mm", "axbmm", "0x120mm", "-5x10mm",
                                                "170x120cm" } ) {
            if ( FigureRenderer.parseSize( bad, 300 ) != null ) {
                return fail( "\"" + bad + "\" must NOT parse as a size" );
            }
        }
        // The DEFAULT size must be physical (mm/inch), not a pixel count. The tree lays out in point-space, so a
        // pixel default at 300 dpi is only a few inches across and the font swamps it -- figures came out cramped
        // and radial labels truncated to a few characters. This is a design decision worth pinning.
        final String default_size = new FigureRenderer.Spec().size;
        if ( ( default_size == null ) || !( default_size.endsWith( "mm" ) || default_size.endsWith( "in" ) ) ) {
            return fail( "the default figure size must be a PHYSICAL size, got \"" + default_size + "\"" );
        }
        return settingsIsolationOk();
    }

    /**
     * The render must start from documented DEFAULTS, not from whatever the operator last saved in the GUI --
     * otherwise the same command yields a different figure on a different machine, which defeats the point of a
     * figure pipeline. isolateSettings() is what guarantees it, by pointing the settings directory somewhere
     * throwaway before any Archaeopteryx class reads it.
     * <p>
     * Checked directly, because the suite ALREADY isolates the settings directory -- so simply rendering twice
     * cannot detect the isolation being removed.
     */
    private static boolean settingsIsolationOk() {
        final String saved = System.getProperty( "archaeopteryx.cache.dir" );
        try {
            System.clearProperty( "archaeopteryx.cache.dir" );
            FigureRenderer.isolateSettings();
            final String set = System.getProperty( "archaeopteryx.cache.dir" );
            if ( ForesterUtilLocal.isEmpty( set ) ) {
                return fail( "isolateSettings must point the settings directory somewhere throwaway" );
            }
            final File real = new File( System.getProperty( "user.home" ), ".archaeopteryx" );
            if ( new File( set ).getAbsolutePath().equals( real.getAbsolutePath() ) ) {
                return fail( "isolateSettings must NOT use the user's own settings directory" );
            }
            // and it must not clobber a directory the caller chose deliberately
            System.setProperty( "archaeopteryx.cache.dir", "/tmp/chosen-by-caller" );
            FigureRenderer.isolateSettings();
            if ( !"/tmp/chosen-by-caller".equals( System.getProperty( "archaeopteryx.cache.dir" ) ) ) {
                return fail( "isolateSettings must leave an explicitly chosen settings directory alone" );
            }
            return true;
        }
        finally {
            if ( saved != null ) {
                System.setProperty( "archaeopteryx.cache.dir", saved );
            }
            else {
                System.clearProperty( "archaeopteryx.cache.dir" );
            }
        }
    }

    /** Local emptiness helper, so this test needs no extra import. */
    private static final class ForesterUtilLocal {

        static boolean isEmpty( final String s ) {
            return ( s == null ) || ( s.trim().length() < 1 );
        }
    }

    private static boolean renderOk() {
        File dir = null;
        try {
            dir = Files.createTempDirectory( "aptxrender" ).toFile();
            final File nh = new File( dir, "t.nh" );
            Files.writeString( nh.toPath(), "((AAAA:0.1,BBBB:0.2)N1:0.1,(CCCC:0.3,DDDD:0.1)N2:0.1)root;\n" );
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( nh, new NHXParser() )[ 0 ];

            // ---- a figure is actually produced, in each family of format --------------------------------
            final FigureRenderer.Spec spec = new FigureRenderer.Spec();
            spec.size = "120x90mm";
            spec.dpi = 150;
            final File png = new File( dir, "a.png" );
            FigureRenderer.render( phy, png, spec );
            if ( !png.exists() || ( png.length() < 1000 ) ) {
                return fail( "no PNG produced (" + png.length() + " bytes)" );
            }
            // PDF is the one that used to write NOTHING: it has its own exporter, and the raster path
            // accepts the request and silently produces no file.
            final File pdf = new File( dir, "a.pdf" );
            FigureRenderer.render( phy, pdf, spec );
            if ( !pdf.exists() || ( pdf.length() < 1000 ) ) {
                return fail( "no PDF produced -- is PDF still routed to the raster path? (" + pdf.length()
                        + " bytes)" );
            }
            final File svg = new File( dir, "a.svg" );
            FigureRenderer.render( phy, svg, spec );
            if ( !svg.exists() || !Files.readString( svg.toPath() ).contains( "<" ) ) {
                return fail( "no SVG produced" );
            }
            // and the un-sized path must work too (it takes a different branch)
            final File png2 = new File( dir, "b.png" );
            final FigureRenderer.Spec plain = new FigureRenderer.Spec();
            FigureRenderer.render( phy, png2, plain );
            if ( !png2.exists() || ( png2.length() < 1000 ) ) {
                return fail( "no PNG produced without an explicit size" );
            }

            // ---- DETERMINISM: the whole point of a figure pipeline --------------------------------------
            final File d1 = new File( dir, "d1.png" );
            final File d2 = new File( dir, "d2.png" );
            FigureRenderer.render( phy, d1, spec );
            FigureRenderer.render( phy, d2, spec );
            if ( !Arrays.equals( Files.readAllBytes( d1.toPath() ), Files.readAllBytes( d2.toPath() ) ) ) {
                return fail( "the same command must produce a byte-identical figure" );
            }

            // ---- -dpi must actually do something, with or without an explicit size ---------------------
            // The unsized raster path used to size from the GUI's raster-export SCALE and ignore DPI entirely,
            // so -dpi=600 and -dpi=72 produced identical files while the help promised otherwise.
            final FigureRenderer.Spec lo = new FigureRenderer.Spec();
            lo.dpi = 72;
            final FigureRenderer.Spec hi = new FigureRenderer.Spec();
            hi.dpi = 300;
            hi.size = lo.size = "100x75mm";
            final File p_lo = new File( dir, "lo.png" );
            final File p_hi = new File( dir, "hi.png" );
            FigureRenderer.render( phy, p_lo, lo );
            FigureRenderer.render( phy, p_hi, hi );
            if ( p_hi.length() <= p_lo.length() ) {
                return fail( "-dpi must change the raster output (300 dpi gave " + p_hi.length()
                        + " bytes, 72 dpi gave " + p_lo.length() + ")" );
            }

            // ---- a radial layout must FILL the canvas, not sit at the hidden frame's size ---------------
            // Switching to circular invalidates the radial diameter; without an explicit fit it is lazily set
            // from the hidden frame's viewport, so the ring came out small and off-centre in the output.
            final FigureRenderer.Spec circ = new FigureRenderer.Spec();
            circ.style = FigureRenderer.Style.CIRCULAR;
            final File c = new File( dir, "circ.png" );
            FigureRenderer.render( phy, c, circ );
            final java.awt.image.BufferedImage ci = javax.imageio.ImageIO.read( c );
            if ( inkExtent( ci ) < 0.6 ) {
                return fail( "a circular figure must fill its canvas; the drawing spans only "
                        + Math.round( 100.0 * inkExtent( ci ) ) + "% of it" );
            }

            // ---- a -color ref the tree cannot honour must be REFUSED, not silently ignored --------------
            final FigureRenderer.Spec bad_color = new FigureRenderer.Spec();
            bad_color.color_by_ref = "data:no_such_property";
            if ( !refused( phy, new File( dir, "e.png" ), bad_color ) ) {
                return fail( "an unknown -color ref must be refused, not silently produce a plain figure" );
            }

            // ---- rendering must not leave the process-wide 'no dialogs' switch on -----------------------
            if ( MainFrame.isRenderingOnly() ) {
                return fail( "render() must restore the rendering-only flag; leaving it on kills the load-time "
                        + "offers for every tree opened afterwards in this JVM" );
            }

            // ---- and it must REFUSE what it cannot do, rather than write nothing quietly ----------------
            if ( !refused( phy, new File( dir, "a.txt" ), new FigureRenderer.Spec() ) ) {
                return fail( "an unknown output extension must be refused" );
            }
            final FigureRenderer.Spec bad_size = new FigureRenderer.Spec();
            bad_size.size = "enormous";
            if ( !refused( phy, new File( dir, "c.png" ), bad_size ) ) {
                return fail( "an unreadable size must be refused" );
            }
            if ( !refused( new Phylogeny(), new File( dir, "d.png" ), new FigureRenderer.Spec() ) ) {
                return fail( "an empty tree must be refused" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "unexpected exception: " + e );
        }
        finally {
            if ( dir != null ) {
                final File[] files = dir.listFiles(); // null if already gone or unreadable -- must not mask the
                if ( files != null ) {                // real result (or a real exception) with an NPE from here
                    for ( final File f : files ) {
                        f.delete();
                    }
                }
                dir.delete();
            }
        }
    }

    /**
     * How much of the canvas the drawing SPANS (its ink bounding box over the canvas), not how many pixels it
     * covers -- a tree is thin line art, so even a perfectly fitted one inks only a few percent of the page.
     */
    private static double inkExtent( final java.awt.image.BufferedImage img ) {
        final int bg = img.getRGB( 1, 1 );
        int minx = Integer.MAX_VALUE, miny = Integer.MAX_VALUE, maxx = -1, maxy = -1;
        for ( int y = 0; y < img.getHeight(); y++ ) {
            for ( int x = 0; x < img.getWidth(); x++ ) {
                if ( img.getRGB( x, y ) != bg ) {
                    minx = Math.min( minx, x );
                    miny = Math.min( miny, y );
                    maxx = Math.max( maxx, x );
                    maxy = Math.max( maxy, y );
                }
            }
        }
        if ( maxx < 0 ) {
            return 0;
        }
        return Math.min( ( ( maxx - minx ) + 1 ) / ( double ) img.getWidth(),
                         ( ( maxy - miny ) + 1 ) / ( double ) img.getHeight() );
    }

    private static boolean refused( final Phylogeny phy, final File out, final FigureRenderer.Spec spec ) {
        try {
            FigureRenderer.render( phy, out, spec );
            return false;
        }
        catch ( final Exception e ) {
            return true;
        }
    }

    /** In a rendered figure nobody can drag the legend off the tree, and its default corner -- top right -- is where
     *  a root-left tree puts its top tips: aptx_render used to draw the legend ON them (and on the "Time tree" badge).
     *  There the legend gets a column of its own, reserved the way annotation columns are. The window is unchanged. */
    private static boolean legendColumnOk() {
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        try {
            final org.forester.io.parsers.nhx.NHXParser p = new org.forester.io.parsers.nhx.NHXParser();
            // long category names: a WIDE legend, which without its column really does lie over the tips; a numeric trait
            // for Size-by; a discrete-trait distribution for the ancestral pies
            // pie states WIDER than the grp colour legend below, so that counting the pie legend where it does not stand
            // (bottom LEFT, beside a colour legend) would change the column
            final String pie = ",loc.set={Upper Guinean Forest Zone,Sahel Transition Zone},loc.set.prob={0.6,0.4}";
            p.setSource( "((isolate_A[&rate=Democratic Republic of the Congo,rate_n=1,loc=Ghana,grp=x]:1.2,isolate_B[&rate=Central African Republic,rate_n=4,loc=Togo,grp=y]:1.2)"
                    + "[&loc=Ghana" + pie + "]:0.9,"
                    + "(isolate_C[&rate=Democratic Republic of the Congo,rate_n=2,loc=Ghana,grp=x]:0.8,(isolate_D[&rate=Central African Republic,rate_n=8,loc=Togo,grp=y]:0.5,"
                    + "isolate_E[&rate=Sao Tome and Principe,rate_n=3,loc=Ghana,grp=x]:0.5)[&loc=Togo" + pie + "]:0.3)[&loc=Ghana" + pie + "]:1.3)[&loc=Ghana" + pie + "];" );
            final Phylogeny phy = p.parse()[ 0 ];
            javax.swing.SwingUtilities.invokeAndWait( () -> {
                mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(), "legend column" );
                mf[ 0 ].setLocation( -32000, -32000 );
                mf[ 0 ].setVisible( true );
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final int w = 1200;
                final int h = 800;
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                // explicitly: this frame reads the developer's own saved preferences (aptx_render isolates itself
                // from them), and a saved root-top orientation would silently turn every check below into a no-op
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                mf[ 0 ].getMainPanel().getControlPanel().demoSelectColorByProperty( "beast:rate" );
                tp.setSize( w, h );
                // the WINDOW: no column, the tree uses the width
                tp.layoutForExportSize( w, h );
                final double free_tip_x = maxTipX( tp, w, h );
                final java.awt.Dimension free_box = tp.sharedLegendSizeForTest();
                if ( ( free_box == null ) || ( ( free_tip_x + tp.getLongestExtNodeInfo() ) <= ( w - 10 - free_box.width ) ) ) {
                    ok[ 0 ] = fail( "fixture: without a column this legend must overlap the tips, or the checks below prove nothing" );
                }
                if ( tp.legendColumnReserve() != 0 ) {
                    ok[ 0 ] = fail( "outside figure rendering no legend column is reserved: the window is unchanged" );
                }
                // FIGURE rendering
                tp.setReserveLegendColumn( true );
                final java.awt.Dimension box = tp.sharedLegendSizeForTest();
                tp.layoutForExportSize( w, h );
                final int reserve = tp.legendColumnReserve();
                final double figure_tip_x = maxTipX( tp, w, h );
                if ( ( box == null ) || ( reserve < ( box.width + 10 ) ) ) {
                    ok[ 0 ] = fail( "the column must hold the whole legend box and its 10 px inset, got reserve " + reserve
                            + " for a box of " + box );
                }
                // the legend's left edge (default corner: w - 10 - box width) must lie right of every tip AND its label
                final int legend_left = w - 10 - ( ( box == null ) ? 0 : box.width );
                if ( ( figure_tip_x + tp.getLongestExtNodeInfo() ) > legend_left ) {
                    ok[ 0 ] = fail( "tips + labels reach x=" + ( figure_tip_x + tp.getLongestExtNodeInfo() )
                            + ", under a legend starting at x=" + legend_left );
                }
                // (not pixel for pixel: the layout re-fits the first depth step inside the narrower budget as well)
                if ( ( free_tip_x - figure_tip_x ) < ( reserve / 2.0 ) ) {
                    ok[ 0 ] = fail( "the column must actually narrow the tree: tips moved only " + ( free_tip_x - figure_tip_x )
                            + " for a reserve of " + reserve );
                }
                // no column where the default corner does not collide, or is not where the legend is
                tp.setLegendOffsetForTest( new java.awt.Point( 40, 40 ) );
                if ( tp.legendColumnReserve() != 0 ) {
                    ok[ 0 ] = fail( "a legend that was moved needs no column" );
                }
                tp.setLegendOffsetForTest( null );
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                if ( tp.legendColumnReserve() != 0 ) {
                    ok[ 0 ] = fail( "a vertical orientation (not reachable from aptx_render) gets no column" );
                }
                tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                if ( tp.legendColumnReserve() != 0 ) {
                    ok[ 0 ] = fail( "a radial tree is a disc and leaves the corners free: no column" );
                }
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.layoutForExportSize( 200, 800 ); // the legend would take over 40% of this width: squeezing the tree is worse
                if ( tp.legendColumnReserve() != 0 ) {
                    ok[ 0 ] = fail( "no column when it would take more than 40% of the width" );
                }
                tp.layoutForExportSize( w, h );
                if ( tp.legendColumnReserve() == 0 ) {
                    ok[ 0 ] = fail( "...and the column is back at a width that has room for it" );
                }
                mf[ 0 ].getMainPanel().getControlPanel().demoSelectColorByProperty( null );
                if ( tp.legendColumnReserve() != 0 ) {
                    ok[ 0 ] = fail( "no legend, no column" );
                }
                // the Size-by legend stands at the right edge too (top right alone): it gets the column on its own
                tp.setSizeByPropertyRef( "beast:rate_n" );
                tp.layoutForExportSize( w, h );
                final int size_reserve = tp.legendColumnReserve();
                final double size_tip_x = maxTipX( tp, w, h );
                if ( size_reserve == 0 ) {
                    ok[ 0 ] = fail( "a Size-by legend alone must get the legend column" );
                }
                else if ( ( size_tip_x + tp.getLongestExtNodeInfo() ) > ( w - size_reserve + LEGEND_GAP_FOR_TEST ) ) {
                    ok[ 0 ] = fail( "with a Size-by legend alone, tips + labels reach x=" + ( size_tip_x + tp.getLongestExtNodeInfo() )
                            + " inside its column of " + size_reserve );
                }
                // ...and the ancestral-pie legend, when it holds the top right itself
                tp.setSizeByPropertyRef( null );
                mf[ 0 ].getMainPanel().getControlPanel().demoSelectAncestralPie( "loc" );
                tp.layoutForExportSize( w, h );
                if ( !tp.isShowAncestralPies() || ( tp.legendColumnReserve() == 0 ) ) {
                    ok[ 0 ] = fail( "an ancestral-pie legend alone at the top right must get the legend column (pies shown: "
                            + tp.isShowAncestralPies() + ")" );
                }
                // beside a colour legend the pie legend drops to the bottom LEFT: only the colour legend is at the right
                // (a NARROW one here, grp = x / y, so a wrongly counted pie legend would widen the column)
                mf[ 0 ].getMainPanel().getControlPanel().demoSelectColorByProperty( "beast:grp" );
                tp.layoutForExportSize( w, h );
                final java.awt.Dimension colour_only = tp.sharedLegendSizeForTest();
                if ( ( colour_only == null ) || ( tp.legendColumnReserve() != ( colour_only.width + 22 ) ) ) {
                    ok[ 0 ] = fail( "beside a colour legend the pie legend moves left, so the column is the colour legend's: "
                            + tp.legendColumnReserve() + " for a box of " + colour_only );
                }
            } );
        }
        catch ( final Exception e ) {
            return fail( "legend column check threw: " + e );
        }
        finally {
            if ( mf[ 0 ] != null ) {
                javax.swing.SwingUtilities.invokeLater( mf[ 0 ]::dispose );
            }
        }
        return ok[ 0 ];
    }

    /** The x of the right-most tip, from a paint into a scratch image at the panel's size (paint assigns the coords). */
    private static final int LEGEND_GAP_FOR_TEST = 12; // the column's gap between the tree and the legend inset

    private static double maxTipX( final TreePanel tp, final int w, final int h ) {
        final java.awt.image.BufferedImage img = new java.awt.image.BufferedImage( w, h, java.awt.image.BufferedImage.TYPE_INT_RGB );
        final java.awt.Graphics2D g = img.createGraphics();
        try {
            tp.paintPhylogeny( g, false, true, w, h, 0, 0 );
        }
        finally {
            g.dispose();
        }
        double max = 0;
        for( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = tp.getPhylogeny().iteratorExternalForward(); it.hasNext(); ) {
            max = Math.max( max, it.next().getXcoord() );
        }
        return max;
    }

    /** aptx_render reads a file the way the WINDOW does. It used to build a bare parser, so a BEAST tree rendered
     *  with no posterior, no node ages and no traits to colour by (-color=beast:rate was refused as unknown), and
     *  the measles Nextstrain export as an undated tree. */
    private static boolean readTreesOk() {
        try {
            final File nex = File.createTempFile( "aptx_render_read", ".nex" );
            nex.deleteOnExit();
            java.nio.file.Files.write( nex.toPath(),
                                       ( "#NEXUS\nbegin trees;\n\ttree t1 = [&R] ((A[&rate=0.5]:1.0,B[&rate=0.7]:1.0)"
                                               + "[&posterior=0.98,height=1.0,height_95%_HPD={0.8,1.2}]:1.0,C[&rate=0.2]:2.0)[&height=2.0];\nend;\n" )
                                                       .getBytes( java.nio.charset.StandardCharsets.UTF_8 ) );
            final Phylogeny[] phys = FigureRenderer.readTrees( nex );
            if ( phys.length != 1 ) {
                return fail( "readTrees must return the one tree of the file" );
            }
            final org.forester.phylogeny.PhylogenyNode ab = phys[ 0 ].getNode( "A" ).getParent();
            if ( !ab.getBranchData().isHasConfidences() || ( ab.getBranchData().getConfidence( 0 ).getValue() != 0.98 ) ) {
                return fail( "readTrees must read a BEAST posterior, as the window does" );
            }
            if ( ( ab.getNodeData().getDate() == null ) || ( ab.getNodeData().getDate().getMax() == null ) ) {
                return fail( "readTrees must read BEAST node ages and their HPD, as the window does" );
            }
            if ( ( phys[ 0 ].getNode( "A" ).getNodeData().getProperties() == null )
                    || phys[ 0 ].getNode( "A" ).getNodeData().getProperties().getProperties( "beast:rate" ).isEmpty() ) {
                return fail( "readTrees must read BEAST traits: they are what -color=beast:rate colours by" );
            }
            // ...and the internal-label promotion it always applied is still applied
            final File nwk = File.createTempFile( "aptx_render_read", ".nwk" );
            nwk.deleteOnExit();
            java.nio.file.Files.write( nwk.toPath(), "((A:1,B:1)95:1,(C:1,D:1)70:1);".getBytes( java.nio.charset.StandardCharsets.UTF_8 ) );
            final Phylogeny bs = FigureRenderer.readTrees( nwk )[ 0 ];
            if ( !bs.getNode( "A" ).getParent().getBranchData().isHasConfidences() ) {
                return fail( "readTrees must still promote bootstrap labels" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "readTrees threw: " + e );
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [AptxRenderTest] " + msg );
        return false;
    }
}
