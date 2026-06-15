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

import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.phylogeny.Phylogeny;

/**
 * Tests for {@link FontResources} -- the bundled Source Sans 3 default figure font.
 *
 * <p>{@link #test()} is headless-suite-safe: it registers the bundled faces and checks they are visible
 * to the default-font selection and resolve (no fallback), plus that the bundled face ships on the
 * classpath (staged by the build's {@code copy_resources}). The full export proof
 * ({@link #headfulExportUsesFont()}, in {@link #main}) drives a real {@link TreePanel} -> SVG and so
 * needs FlatLaf + a display; it is a no-op when headless (the established pattern for GUI tests here).
 */
public final class FontResourcesTest {

    public static void main( final String[] args ) {
        final boolean unit = test();
        System.out.println( "FontResources (registration): " + ( unit ? "OK." : "FAILED." ) );
        final boolean e2e = headfulExportUsesFont();
        System.out.println( "FontResources (export uses font): " + ( e2e ? "OK/skipped." : "FAILED." ) );
        System.exit( ( unit && e2e ) ? 0 : 1 );
    }

    public static boolean test() {
        FontResources.registerBundledFonts(); // idempotent; refreshes AptxUtil's available-font list
        final String[] available = AptxUtil.getAvailableFontFamiliesSorted();
        // all three bundled families must ship, register, and resolve (plain/bold/italic -> real faces)
        for( final FontResources.Preferred p : FontResources.PREFERRED ) {
            if ( FontResources.class.getResourceAsStream( FontResources.FONT_DIR + p.faces[ 0 ] ) == null ) {
                return fail( "bundled font not on the classpath: " + p.family + " -- did 'copy_resources' run?" );
            }
            if ( Arrays.binarySearch( available, p.family ) < 0 ) {
                return fail( p.family + " was not registered / not visible to AptxUtil" );
            }
            if ( !p.family.equals( new Font( p.family, Font.PLAIN, 12 ).getFamily() )
                    || !p.family.equals( new Font( p.family, Font.BOLD, 12 ).getFamily() )
                    || !p.family.equals( new Font( p.family, Font.ITALIC, 12 ).getFamily() ) ) {
                return fail( p.family + ": plain/bold/italic must resolve to the family (not a fallback)" );
            }
        }
        // the default must be the first preferred family
        if ( !FontResources.DEFAULT_FIGURE_FONT.equals( FontResources.PREFERRED[ 0 ].family )
                || !FontResources.DEFAULT_FIGURE_FONT.equals( AptxConstants.DEFAULT_FONT_CHOICES[ 0 ] ) ) {
            return fail( "DEFAULT_FONT_CHOICES[0] / PREFERRED[0] should be the bundled default font" );
        }
        return true;
    }

    /** Full tree -> SVG export must carry the bundled family as the figure font. No-op (true) when headless. */
    private static boolean headfulExportUsesFont() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            FontResources.registerBundledFonts(); // before Configuration loads, so the default flips to it
            final Configuration conf = new Configuration();
            final Phylogeny phy = Phylogeny.createInstanceFromNhxString( "((Apis,Bombus)x,(Felis,Canis)y)root" );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "font e2e" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final JFrame fr = (JFrame) mf[ 0 ];
                    fr.setSize( 1100, 800 );
                    fr.validate();
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    mp.getControlPanel().showWhole();
                    if ( ( tp.getWidth() < 200 ) || ( tp.getHeight() < 200 ) ) {
                        fr.dispose();
                        return; // no usable viewport in this environment
                    }
                    final String dir = System.getProperty( "java.io.tmpdir" );
                    // SVG keeps text as <text>, so we can assert the figure font by name
                    final File svgf = new File( dir, "aptx_font_tree.svg" );
                    AptxUtil.writePhylogenyToGraphicsFile( svgf.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                           mp.getControlPanel(), GraphicsExportType.SVG, mp.getOptions() );
                    final String svg = new String( Files.readAllBytes( svgf.toPath() ), StandardCharsets.UTF_8 );
                    if ( !svg.contains( FontResources.DEFAULT_FIGURE_FONT ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  SVG does not reference the bundled figure font '"
                                + FontResources.DEFAULT_FIGURE_FONT + "'" );
                    }
                    else {
                        System.out.println( "  SVG figure font = '" + FontResources.DEFAULT_FIGURE_FONT + "' -> " + svgf );
                    }
                    // PDF (outlines) and PNG (raster) share the identical paintPhylogeny path; just confirm they
                    // render non-trivially with the font registered (text is not inspectable by name there)
                    final File pdf = new File( dir, "aptx_font_tree.pdf" ); // PDF has its own exporter path
                    PdfExporter.writePhylogenyToPdf( pdf.getAbsolutePath(), tp, tp.getWidth(), tp.getHeight() );
                    final File png = new File( dir, "aptx_font_tree.png" );
                    AptxUtil.writePhylogenyToGraphicsFile( png.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                           mp.getControlPanel(), GraphicsExportType.PNG, mp.getOptions() );
                    if ( ( pdf.length() < 1000L ) || ( png.length() < 1000L ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  PDF/PNG export looks empty: pdf=" + pdf.length() + " png=" + png.length() );
                    }
                    else {
                        System.out.println( "  PDF -> " + pdf + " (" + pdf.length() + " b) ; PNG -> " + png + " ("
                                + png.length() + " b)" );
                    }
                    fr.dispose();
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

    private static boolean fail( final String msg ) {
        System.out.println( "  FontResources test failed: " + msg );
        return false;
    }

    private FontResourcesTest() {
    }
}
