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
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.function.Consumer;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.archaeopteryx.VectorGraphicsExporter.Format;
import org.forester.phylogeny.Phylogeny;

/**
 * Tests for {@link VectorGraphicsExporter}.
 *
 * <p>{@link #test()} exercises the pure, GUI-free {@link VectorGraphicsExporter#render} for SVG and
 * EPS (a simple line + label painter) and is wired into the headless suite -- it needs only the
 * VectorGraphics2D library, no display. The full tree-export path
 * ({@link VectorGraphicsExporter#writePhylogenyToVectorGraphicsFile} via
 * {@link AptxUtil#writePhylogenyToGraphicsFile}) drives a real {@link TreePanel} and so needs
 * FlatLaf + a display; that runs only from {@link #main} when not headless (the established pattern
 * for the GUI tests in this package).
 */
public final class VectorGraphicsExporterTest {

    private static final String LABEL = "LEAFXYZ";

    public static void main( final String[] args ) {
        final boolean unit = test();
        System.out.println( "VectorGraphicsExporter (render SVG/EPS): " + ( unit ? "OK." : "FAILED." ) );
        final boolean e2e = headfulTreeExport();
        System.out.println( "VectorGraphicsExporter (tree export)   : " + ( e2e ? "OK/skipped." : "FAILED." ) );
        System.exit( ( unit && e2e ) ? 0 : 1 );
    }

    public static boolean test() {
        try {
            final Consumer<Graphics2D> painter = g -> {
                g.setColor( Color.BLACK );
                g.drawLine( 10, 10, 290, 10 );
                g.setFont( new Font( "SansSerif", Font.PLAIN, 12 ) );
                g.drawString( LABEL, 20, 40 );
            };
            // ---- SVG, outlining ON (default): text is vectorized to glyph outlines (no viewer font
            // substitution), so the literal label string must NOT appear ----
            final String svg = new String( VectorGraphicsExporter.render( 300, 120, Format.SVG, true, painter ),
                                           StandardCharsets.UTF_8 );
            if ( !svg.contains( "<svg" ) || !svg.contains( "</svg>" ) ) {
                return fail( "svg: not a well-formed <svg> document: " + head( svg ) );
            }
            if ( !svg.contains( "viewBox=\"0 0 300 120\"" ) ) {
                return fail( "svg: viewBox should equal the pixel size 300x120: " + head( svg ) );
            }
            if ( svg.contains( LABEL ) ) {
                return fail( "svg(outline): label must be vectorized, not emitted as substitutable text" );
            }
            // ---- SVG, outlining OFF (the opt-out): the label is kept as real, selectable text ----
            final String svg_text = new String( VectorGraphicsExporter.render( 300, 120, Format.SVG, false, painter ),
                                                StandardCharsets.UTF_8 );
            if ( !svg_text.contains( LABEL ) ) {
                return fail( "svg(no-outline): label should be emitted as real text when outlining is off" );
            }
            // Outlining must turn the label into vector path content -> strictly MORE <path> than the
            // text version; a backend that silently dropped the label could not satisfy this.
            if ( countPaths( svg ) <= countPaths( svg_text ) ) {
                return fail( "svg: outlining must add vector paths vs the text version (" + countPaths( svg )
                        + " vs " + countPaths( svg_text ) + ")" );
            }

            // ---- EPS, outlining ON: valid encapsulated PostScript, bbox ~ pixel size; the literal
            // label string must NOT appear (vectorized) ----
            final String eps = new String( VectorGraphicsExporter.render( 300, 120, Format.EPS, true, painter ),
                                           StandardCharsets.ISO_8859_1 );
            if ( !eps.startsWith( "%!PS-Adobe" ) ) {
                return fail( "eps: bad header: " + head( eps ) );
            }
            if ( !eps.contains( "%%BoundingBox" ) ) {
                return fail( "eps: missing %%BoundingBox" );
            }
            if ( eps.contains( LABEL ) ) {
                return fail( "eps: label must be vectorized, not emitted as substitutable text" );
            }
            final int bbox_w = epsBoundingBoxWidth( eps );
            if ( Math.abs( bbox_w - 300 ) > 3 ) {
                return fail( "eps: bounding-box width " + bbox_w + " pt should be ~300 (1 px -> 1 pt)" );
            }

            // ---- EPS content must FIT the page. Regression guard for the "too zoomed in" bug: the
            // pixel-valued commands used to be ~2.83x too big for the mm-sized page, so only a top-left
            // corner of the figure appeared. The fix makes the NET CTM scale (VG2D's mm->pt page scale
            // composed with our px->mm pre-scale) exactly 1.0 -- a pixel coordinate then lands at the same
            // number of points, and since the page is sized to the pixel dimensions in points (the bbox
            // check above), EVERY drawn coordinate stays within the page, regardless of glyph shape. Assert
            // that invariant (the robust guarantee), plus a concrete "the corner content lands on the page"
            // check. A line-only painter keeps exactly the two expected scale operators (no per-glyph
            // scales), so the net-scale product is unambiguous. Without the pre-scale the net scale is ~2.83
            // and the figure overflows. ----
            final int fw = 400;
            final int fh = 200;
            final String eps_fit = new String( VectorGraphicsExporter.render( fw, fh, Format.EPS, true, g -> {
                g.setColor( Color.BLACK );
                g.drawLine( 0, 0, fw - 1, fh - 1 ); // reaches the far corner
            } ), StandardCharsets.ISO_8859_1 );
            final double net_scale = epsCumulativeScale( eps_fit );
            if ( Math.abs( net_scale - 1.0 ) > 0.01 ) {
                return fail( "eps: net CTM scale should be ~1.0 (px->pt identity) so the whole figure fits; got "
                        + net_scale + " (the drawing would be that many times too big for the page)" );
            }
            final double[] max_pt = epsMaxContentPoint( eps_fit );
            if ( ( max_pt[ 0 ] > ( fw + 2 ) ) || ( max_pt[ 1 ] > ( fh + 2 ) ) ) {
                return fail( "eps: drawing overflows the page (content " + Math.round( max_pt[ 0 ] ) + "x"
                        + Math.round( max_pt[ 1 ] ) + " pt vs page " + fw + "x" + fh
                        + "): only part of the figure would be in the file" );
            }
            // and it must actually reach near the far corner, so a no-op/collapsed painter can't pass vacuously
            if ( ( max_pt[ 0 ] < ( fw * 0.9 ) ) || ( max_pt[ 1 ] < ( fh * 0.9 ) ) ) {
                return fail( "eps: content should span nearly the whole page; got " + Math.round( max_pt[ 0 ] ) + "x"
                        + Math.round( max_pt[ 1 ] ) + " pt" );
            }

            // ---- format ids line up with the AptxUtil export-type suffixes ----
            if ( !Format.SVG.id().equals( GraphicsExportType.SVG.toString() )
                    || !Format.EPS.id().equals( GraphicsExportType.EPS.toString() ) ) {
                return fail( "format id / GraphicsExportType suffix mismatch" );
            }

            // ---- Regression: domain-tree SVG/EPS export (Archaeopteryx's only GradientPaint is the domain box).
            // VectorGraphics2D rasterizes a non-Color (gradient) fill into a BufferedImage the size of the shape's
            // ROUNDED bounds -- and it stays ARMED after a gradient (cleared only by a dispose the single-context
            // paint pass never issues), so EVERY later fill is rasterized too. Two failures, both reproduced on real
            // domain trees: (a) a following SUB-PIXEL fill (a second domain's Color shadow, an empty glyph) rounds to
            // 0 -> new BufferedImage(0, ..) -> IllegalArgumentException aborts the whole export; (b) normal content
            // after a gradient (labels, boxes) embeds as a gradient-TINTED <image> instead of vector.
            // GuardedVectorGraphics2D flattens the gradient to a solid Color up front, so the filter is never armed.
            // This painter reproduces the stale-arming chain; exercises BOTH the outlining and non-outlining
            // Graphics2D, for BOTH formats. If render() threw (crash not fixed), the outer catch fails the test.
            final Consumer<Graphics2D> gradient_painter = g -> {
                g.setPaint( new java.awt.GradientPaint( 0, 0, Color.RED, 0, 10, Color.BLUE ) );
                g.fill( new java.awt.geom.Rectangle2D.Float( 20, 20, 40f, 12f ) ); // a domain body (gradient)
                g.setColor( new Color( 8, 18, 21, 90 ) );                          // a following shadow (Color)
                g.fill( new java.awt.geom.Rectangle2D.Float( 5, 5, 0.3f, 8f ) );   // sub-pixel AFTER gradient -> crashed
                g.setColor( Color.BLACK );
                g.drawString( "X", 60, 40 );                                       // a label after gradient -> tinted
            };
            for( final Format f : new Format[] { Format.SVG, Format.EPS } ) {
                for( final boolean outline : new boolean[] { true, false } ) {
                    final byte[] doc = VectorGraphicsExporter.render( 120, 60, f, outline, gradient_painter );
                    if ( doc.length < 80 ) { // a real document, not a trivial stub
                        return fail( f + "(outline=" + outline + "): domain-style gradient export produced no document" );
                    }
                }
            }
            // the gradient content must be CLEAN VECTOR, not a rasterized (and tinted) <image>: flattening the
            // gradient to a solid Color yields a real <rect>, and nothing after it is rasterized.
            final String grad_svg = new String( VectorGraphicsExporter.render( 120, 60, Format.SVG, true,
                                                                               gradient_painter ),
                                                StandardCharsets.UTF_8 );
            if ( grad_svg.contains( "<image" ) ) {
                return fail( "svg: a gradient fill must export as flat vector, not a rasterized <image> (which also "
                        + "tints the content drawn after it): " + head( grad_svg ) );
            }
            if ( !grad_svg.contains( "<rect" ) ) {
                return fail( "svg: the flattened gradient box should render as a <rect>: " + head( grad_svg ) );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( e.toString() );
        }
    }

    /** Full tree -> SVG export through the public dispatch; no-op (returns true) when headless. */
    private static boolean headfulTreeExport() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            final Phylogeny phy = Phylogeny
                    .createInstanceFromNhxString( "(((A,B)x,(C,D)y)p,((E,F)z,(G,H,I)w)q)root" );
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "vg e2e" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final JFrame f = (JFrame) mf[ 0 ];
                    f.setSize( 1100, 800 );
                    f.validate();
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    mp.getControlPanel().showWhole();
                    if ( ( tp.getWidth() < 200 ) || ( tp.getHeight() < 200 ) ) {
                        f.dispose();
                        return; // no usable viewport in this environment; nothing to assert
                    }
                    final File out = new File( System.getProperty( "java.io.tmpdir" ), "aptx_vg_realtree.svg" );
                    // through AptxUtil.writePhylogenyToGraphicsFile -> exercises the enum + dispatch + exporter
                    final String written = AptxUtil.writePhylogenyToGraphicsFile( out.getAbsolutePath(),
                                                                                  tp.getWidth(),
                                                                                  tp.getHeight(),
                                                                                  tp,
                                                                                  mp.getControlPanel(),
                                                                                  GraphicsExportType.SVG,
                                                                                  mp.getOptions() );
                    final String svg = new String( Files.readAllBytes( out.toPath() ), StandardCharsets.UTF_8 );
                    // labels are vectorized to outlines, so they are not literal text; instead require a
                    // path per drawn glyph/branch -- 9 leaf labels + internal labels + branches => many.
                    final int paths = svg.split( "<path", -1 ).length - 1;
                    if ( paths < 9 ) {
                        ok[ 0 ] = false;
                        System.out.println( "  too few vector paths in svg (labels/branches missing?): " + paths );
                    }
                    if ( svg.contains( ">A<" ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  leaf labels should be outlined, not emitted as substitutable text" );
                    }
                    if ( !svg.contains( "<svg" ) || !svg.contains( "</svg>" ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  malformed svg" );
                    }
                    System.out.println( "  wrote " + written + " ; " + out.length() + " bytes -> " + out );
                    f.dispose();
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                    ok[ 0 ] = TestFail.here();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    private static int countPaths( final String svg ) {
        return svg.split( "<path", -1 ).length - 1;
    }

    /**
     * The net CTM scale VectorGraphics2D applies to drawing coordinates = the product of the |x| factors of
     * every {@code scale} operator. VG2D's EPS processor prepends the mm-&gt;pt page scale and (post-fix) our
     * px-&gt;mm pre-scale; their product is ~1.0 by design, so a pixel coordinate lands at the same number of
     * points. Callers should use a line-only render (no per-glyph {@code scale} ops) so the product is exactly
     * the page+pre scales. The x factor suffices -- the two scales are uniform (|x| == |y|).
     */
    private static double epsCumulativeScale( final String eps ) {
        double scale = 1;
        final java.util.regex.Matcher sc = java.util.regex.Pattern.compile( "([-0-9.]+)\\s+([-0-9.]+)\\s+scale" )
                .matcher( eps );
        while ( sc.find() ) {
            scale *= Math.abs( Double.parseDouble( sc.group( 1 ) ) );
        }
        return scale;
    }

    /**
     * The furthest drawing-command coordinate, in PostScript points: {@code max|M/L coord| * net scale}
     * (see {@link #epsCumulativeScale}), per axis. Multiplying the raw command coordinate by the net scale
     * gives where it actually lands on the page. Intended for a line-only render, whose extent is straight
     * moveto/lineto (M/L); it is not a full path parser (it ignores curveto control points), so the
     * page-fit guarantee for arbitrary glyphs rests on the net-scale-1.0 invariant, not on this scan.
     */
    private static double[] epsMaxContentPoint( final String eps ) {
        final double scale = epsCumulativeScale( eps );
        double max_x = 0;
        double max_y = 0;
        // M and L are this EPS's aliases for moveto/lineto (see the prologue "/M /moveto load def").
        final java.util.regex.Matcher cmd = java.util.regex.Pattern.compile( "([-0-9.]+)\\s+([-0-9.]+)\\s+[ML]\\b" )
                .matcher( eps );
        while ( cmd.find() ) {
            max_x = Math.max( max_x, Math.abs( Double.parseDouble( cmd.group( 1 ) ) ) );
            max_y = Math.max( max_y, Math.abs( Double.parseDouble( cmd.group( 2 ) ) ) );
        }
        return new double[] { max_x * scale, max_y * scale };
    }

    private static int epsBoundingBoxWidth( final String eps ) {
        // tolerate any line ending (\n, \r\n, \r) and surrounding whitespace; require the four numbers.
        for( final String raw : eps.split( "\\r\\n|\\r|\\n" ) ) {
            final String line = raw.trim();
            if ( line.startsWith( "%%BoundingBox:" ) ) {
                final String[] p = line.split( "\\s+" ); // %%BoundingBox: x0 y0 x1 y1
                if ( p.length < 5 ) {
                    return -1;
                }
                return Math.round( Float.parseFloat( p[ 3 ] ) );
            }
        }
        return -1;
    }

    private static String head( final String s ) {
        return s.substring( 0, Math.min( 120, s.length() ) ).replace( "\n", " " );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "\nVectorGraphicsExporterTest failed: " + msg );
        return false;
    }

    private VectorGraphicsExporterTest() {
    }
}
