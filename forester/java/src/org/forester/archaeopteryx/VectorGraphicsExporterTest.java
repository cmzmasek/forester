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
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

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
            // ---- SVG: viewBox is the pixel size, content is real vector text + lines ----
            final String svg = new String( VectorGraphicsExporter.render( 300, 120, Format.SVG, g -> {
                g.setColor( Color.BLACK );
                g.drawLine( 10, 10, 290, 10 );
                g.setFont( new Font( "SansSerif", Font.PLAIN, 12 ) );
                g.drawString( LABEL, 20, 40 );
            } ), StandardCharsets.UTF_8 );
            if ( !svg.contains( "<svg" ) || !svg.contains( "</svg>" ) ) {
                return fail( "svg: not a well-formed <svg> document: " + head( svg ) );
            }
            if ( !svg.contains( "viewBox=\"0 0 300 120\"" ) ) {
                return fail( "svg: viewBox should equal the pixel size 300x120: " + head( svg ) );
            }
            if ( !svg.contains( LABEL ) ) {
                return fail( "svg: label must be emitted as real text (not vectorized)" );
            }
            if ( !svg.contains( "<line" ) && !svg.contains( "<path" ) ) {
                return fail( "svg: the drawn line is missing" );
            }

            // ---- EPS: valid encapsulated PostScript, bounding box ~ pixel size (1 px -> 1 pt) ----
            final String eps = new String( VectorGraphicsExporter.render( 300, 120, Format.EPS, g -> {
                g.setColor( Color.BLACK );
                g.drawLine( 10, 10, 290, 10 );
                g.setFont( new Font( "SansSerif", Font.PLAIN, 12 ) );
                g.drawString( LABEL, 20, 40 );
            } ), StandardCharsets.ISO_8859_1 );
            if ( !eps.startsWith( "%!PS-Adobe" ) ) {
                return fail( "eps: bad header: " + head( eps ) );
            }
            if ( !eps.contains( "%%BoundingBox" ) ) {
                return fail( "eps: missing %%BoundingBox" );
            }
            if ( !eps.contains( LABEL ) ) {
                return fail( "eps: label text missing" );
            }
            final int bbox_w = epsBoundingBoxWidth( eps );
            if ( Math.abs( bbox_w - 300 ) > 3 ) {
                return fail( "eps: bounding-box width " + bbox_w + " pt should be ~300 (1 px -> 1 pt)" );
            }

            // ---- format ids line up with the AptxUtil export-type suffixes ----
            if ( !Format.SVG.id().equals( GraphicsExportType.SVG.toString() )
                    || !Format.EPS.id().equals( GraphicsExportType.EPS.toString() ) ) {
                return fail( "format id / GraphicsExportType suffix mismatch" );
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
            final Configuration conf = new Configuration( null, false, false, true );
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
                    for( final String leaf : new String[] { "A", "B", "C", "D", "E", "F", "G", "H", "I" } ) {
                        if ( !svg.contains( ">" + leaf + "<" ) ) {
                            ok[ 0 ] = false;
                            System.out.println( "  missing leaf label in svg: " + leaf );
                        }
                    }
                    if ( !svg.contains( "<path" ) && !svg.contains( "<line" ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  no branches drawn in svg" );
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
                    ok[ 0 ] = false;
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
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
