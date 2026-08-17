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

import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.forester.archaeopteryx.TanglegramExporter.Format;
import org.forester.archaeopteryx.TanglegramLinker.LinkField;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Tests {@link TanglegramExporter}: each format writes a non-empty file carrying that format's signature (PDF header,
 * SVG/EPS text, PNG magic). Rendering into a supplied Graphics2D needs no display, so this runs headless (it does need
 * the vendored VectorGraphics2D / OpenPDF libraries on the classpath, like VectorGraphicsExporterTest).
 */
public final class TanglegramExporterTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TanglegramExporter: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return exportSizeOk() && exportFilesOk();
    }

    private static boolean exportSizeOk() {
        final TanglegramPanel panel = new TanglegramPanel( abcd(), abcd(), LinkField.NODE_NAME );
        panel.setSize( 100, 100 ); // below the floor -> floored to 600 x 400
        final Dimension floored = panel.exportSize();
        if ( ( floored.width != 600 ) || ( floored.height != 400 ) ) {
            return fail( "a small panel should export at the 600x400 floor, got " + floored.width + "x"
                    + floored.height );
        }
        panel.setSize( 1200, 800 ); // above the floor -> its own size
        final Dimension actual = panel.exportSize();
        if ( ( actual.width != 1200 ) || ( actual.height != 800 ) ) {
            return fail( "a large panel should export at its own size, got " + actual.width + "x" + actual.height );
        }
        return true;
    }

    private static boolean exportFilesOk() {
        final TanglegramPanel panel = new TanglegramPanel( abcd(), abcd(), LinkField.NODE_NAME );
        final int w = 700;
        final int h = 400;
        try {
            for( final Format format : Format.values() ) {
                final File file = File.createTempFile( "tanglegram-export-", "." + format.extension() );
                file.deleteOnExit();
                try {
                    TanglegramExporter.write( file, format, panel, w, h );
                    if ( !file.exists() || ( file.length() <= 0 ) ) {
                        return fail( format + " wrote a missing or empty file" );
                    }
                    if ( !renderedOk( format, file ) ) {
                        return fail( format + " file lacks the format signature or has no rendered content" );
                    }
                }
                finally {
                    file.delete();
                }
            }
            return true;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return fail( "export threw: " + e );
        }
    }

    /** True when the file has the right format signature AND actually contains rendered content (not an empty
     *  container): the vector text formats must carry drawing primitives (SVG {@code <line>}s, EPS {@code stroke}s);
     *  the binary formats must clear a size floor an empty figure could not reach. */
    private static boolean renderedOk( final Format format, final File file ) throws IOException {
        final byte[] bytes = Files.readAllBytes( file.toPath() );
        final String text = new String( bytes, StandardCharsets.ISO_8859_1 );
        switch ( format ) {
            case PDF:
                return text.startsWith( "%PDF" ) && ( bytes.length > 1200 );
            case SVG:
                return text.contains( "<svg" ) && text.contains( "<line" );
            case EPS:
                return text.contains( "%!PS" ) && text.contains( "stroke" );
            case PNG:
                return ( bytes.length > 8 ) && ( ( bytes[ 0 ] & 0xFF ) == 0x89 ) && ( bytes[ 1 ] == 'P' )
                        && ( bytes[ 2 ] == 'N' ) && ( bytes[ 3 ] == 'G' ) && ( bytes.length > 4000 );
            default:
                return false;
        }
    }

    private static boolean fail( final String message ) {
        System.out.println( "TanglegramExporter test failed: " + message );
        return false;
    }

    private static Phylogeny abcd() {
        return tree( clade( clade( leaf( "A" ), leaf( "B" ) ), clade( leaf( "C" ), leaf( "D" ) ) ) );
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode clade( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode child : children ) {
            n.addAsChild( child );
        }
        return n;
    }

    private static Phylogeny tree( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }
}
