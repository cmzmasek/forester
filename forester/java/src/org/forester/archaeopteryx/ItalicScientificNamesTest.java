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
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * End-to-end render test for italic scientific names ({@link Options#isUseItalicScientificNames()},
 * applied by {@link TreePanel#taxonomyLabel}). Drives a real {@link TreePanel} through the PNG export
 * path and checks two things by comparing the ON vs OFF renders pixel-for-pixel:
 *
 * <ol>
 * <li>with <i>scientific names</i> shown, toggling italics CHANGES the image (italic glyphs differ from
 * roman) -- so the option actually takes effect; and</li>
 * <li>with only the <i>taxonomy code</i> shown (a non-scientific field), toggling italics leaves the
 * image IDENTICAL -- so italics are scoped to the scientific name and never leak onto code/common/rank.</li>
 * </ol>
 *
 * <p>Needs FlatLaf + a display, so {@link #test()} is a no-op (returns true) when headless -- the same
 * pattern the other GUI tests in this package use; it runs for real from {@link #main} or any
 * non-headless invocation.
 */
public final class ItalicScientificNamesTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ItalicScientificNames: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = taxonomyTree();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "italic test" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final JFrame f = (JFrame) mf[ 0 ];
                try {
                    f.setSize( 760, 460 );
                    f.validate();
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    final ControlPanel cp = mp.getControlPanel();
                    cp.setCheckbox( Configuration.display_external_data, true );
                    cp.showWhole();
                    if ( ( tp.getWidth() < 200 ) || ( tp.getHeight() < 200 ) ) {
                        return; // no usable viewport in this environment; nothing to assert
                    }
                    // (1) scientific names shown: italics ON vs OFF must differ
                    cp.setCheckbox( Configuration.show_taxonomy_scientific_names, true );
                    cp.setCheckbox( Configuration.show_tax_code, false );
                    if ( !cp.isShowTaxonomyScientificNames() ) {
                        return; // can't enable scientific-name display here; nothing to assert
                    }
                    if ( pixelDiff( render( tp, mp, true ), render( tp, mp, false ) ) == 0 ) {
                        ok[ 0 ] = false;
                        System.out.println( "  italics had no effect on a scientific name" );
                    }
                    // (2) only the taxonomy code shown: italics ON vs OFF must be identical
                    cp.setCheckbox( Configuration.show_taxonomy_scientific_names, false );
                    cp.setCheckbox( Configuration.show_tax_code, true );
                    if ( cp.isShowTaxonomyCode() ) {
                        final long leak = pixelDiff( render( tp, mp, true ), render( tp, mp, false ) );
                        if ( leak != 0 ) {
                            ok[ 0 ] = false;
                            System.out.println( "  italics leaked onto a non-scientific field (code); diff=" + leak );
                        }
                    }
                    // (3) vector (SVG/EPS) export: VectorGraphics2D embeds no fonts, so every label is
                    // rendered as glyph OUTLINES to stop the viewer substituting the bundled face (EPS ->
                    // Times serif, SVG -> generic sans). So NO label text -- scientific name OR taxonomy
                    // code -- may appear as literal characters, <path> outlines must be present, and
                    // toggling italics must still change the output (italic vs roman species geometry).
                    cp.setCheckbox( Configuration.show_taxonomy_scientific_names, true );
                    cp.setCheckbox( Configuration.show_tax_code, true );
                    if ( cp.isShowTaxonomyScientificNames() ) {
                        final String svg_on = exportSvg( tp, mp, true );
                        final String svg_off = exportSvg( tp, mp, false );
                        if ( svg_on.contains( "Homo sapiens" ) || svg_on.contains( "HUMAN" )
                                || !svg_on.contains( "<path" ) ) {
                            ok[ 0 ] = false;
                            System.out.println( "  SVG export did not outline all label text (font substitution risk)" );
                        }
                        if ( svg_on.equals( svg_off ) ) {
                            ok[ 0 ] = false;
                            System.out.println( "  italics had no effect on the vector export" );
                        }
                    }
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                    ok[ 0 ] = false;
                }
                finally {
                    f.dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    private static BufferedImage render( final TreePanel tp, final MainPanel mp, final boolean italic )
            throws Exception {
        mp.getOptions().setUseItalicScientificNames( italic );
        final File out = File.createTempFile( "aptx_italic_" + italic, ".png" );
        try {
            AptxUtil.writePhylogenyToGraphicsFile( out.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                   mp.getControlPanel(), GraphicsExportType.PNG, mp.getOptions() );
            return ImageIO.read( out );
        }
        finally {
            out.delete();
        }
    }

    private static String exportSvg( final TreePanel tp, final MainPanel mp, final boolean italic )
            throws Exception {
        mp.getOptions().setUseItalicScientificNames( italic );
        final File out = File.createTempFile( "aptx_italic_" + italic, ".svg" );
        try {
            AptxUtil.writePhylogenyToGraphicsFile( out.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                   mp.getControlPanel(), GraphicsExportType.SVG, mp.getOptions() );
            return new String( Files.readAllBytes( out.toPath() ), StandardCharsets.UTF_8 );
        }
        finally {
            out.delete();
        }
    }

    /** Number of pixels that differ between two same-size images (or a large sentinel if sizes differ). */
    private static long pixelDiff( final BufferedImage a, final BufferedImage b ) {
        if ( ( a.getWidth() != b.getWidth() ) || ( a.getHeight() != b.getHeight() ) ) {
            return Long.MAX_VALUE;
        }
        long n = 0;
        for( int y = 0; y < a.getHeight(); y++ ) {
            for( int x = 0; x < a.getWidth(); x++ ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Two tips, each with a scientific (binomial) name and a taxonomy code. */
    private static Phylogeny taxonomyTree() throws Exception {
        final Phylogeny p = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( tip( "Homo sapiens", "HUMAN" ) );
        root.addAsChild( tip( "Mus musculus", "MOUSE" ) );
        p.setRoot( root );
        p.externalNodesHaveChanged();
        return p;
    }

    private static PhylogenyNode tip( final String scientific_name, final String code ) throws Exception {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        t.setTaxonomyCode( code );
        n.getNodeData().addTaxonomy( t );
        return n;
    }

    private ItalicScientificNamesTest() {
    }
}
