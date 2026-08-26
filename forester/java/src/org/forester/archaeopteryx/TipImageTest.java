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
import java.awt.image.BufferedImage;
import java.io.File;

import javax.imageio.ImageIO;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;

/**
 * Unit tests for the tip-images foundation: {@link TipImages} (reference resolution + size math) and
 * {@link TipImageCache} (local-file decode + the async cache). Headless (ImageIO works without a display); run via
 * the suite or {@link #main(String[])}.
 */
public final class TipImageTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TipImage: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return testRefResolution() && testPredicates() && testScaledSize() && testLoadLocal() && testAsyncCache();
    }

    // ---- resolve a tip's image reference from a property or a taxonomy <uri> ----
    private static boolean testRefResolution() {
        // a recognized image PROPERTY (by ref namespace/local name)
        if ( !"figs/a.png".equals( TipImages.imageRefFor( withProperty( "image:src", "figs/a.png" ) ) ) ) {
            return fail( "an image:src property should be the image ref" );
        }
        if ( !"http://x/y.jpg".equals( TipImages.imageRefFor( withProperty( "data:photo", "http://x/y.jpg" ) ) ) ) {
            return fail( "a data:photo property should be the image ref" );
        }
        // a non-image-named property whose VALUE has an image extension still counts
        if ( !"s.gif".equals( TipImages.imageRefFor( withProperty( "data:whatever", "s.gif" ) ) ) ) {
            return fail( "a value ending in an image extension should be recognized" );
        }
        // a plain non-image property is NOT a ref
        if ( TipImages.imageRefFor( withProperty( "data:host", "human" ) ) != null ) {
            return fail( "a non-image property must not be treated as an image ref" );
        }
        // the phyloXML taxonomy <uri type="image_url">
        if ( !"http://x/t.png".equals( TipImages.imageRefFor( withTaxonomyUri( "http://x/t.png", "image_url" ) ) ) ) {
            return fail( "a taxonomy image <uri> should be the image ref" );
        }
        // a non-image uri type without an image extension is ignored
        if ( TipImages.imageRefFor( withTaxonomyUri( "http://x/page", "webpage" ) ) != null ) {
            return fail( "a non-image <uri> must not be treated as an image ref" );
        }
        // trimming (a real property ref always carries a namespace colon, e.g. from Import Annotations: data:image)
        if ( !"a.png".equals( TipImages.imageRefFor( withProperty( "data:image", "  a.png  " ) ) ) ) {
            return fail( "the image ref should be trimmed" );
        }
        return true;
    }

    private static boolean testPredicates() {
        if ( !TipImages.isImageRef( "image:src", "x" ) || !TipImages.isImageRef( "img:foo", "x" )
                || !TipImages.isImageRef( "data:image", "x" ) || !TipImages.isImageRef( "silhouette", "x" ) ) {
            return fail( "recognized image refs should match on namespace or local name" );
        }
        if ( TipImages.isImageRef( "data:host", "human" ) ) {
            return fail( "a non-image ref with a non-image value should not match" );
        }
        if ( !TipImages.hasImageExtension( "a/b/c.PNG" ) || !TipImages.hasImageExtension( "http://x/y.jpg?v=2" )
                || TipImages.hasImageExtension( "a.svg" ) || TipImages.hasImageExtension( "noext" ) ) {
            return fail( "hasImageExtension: raster exts (case-insensitive, query-stripped) yes, svg/none no" );
        }
        if ( !TipImages.hasImageExtension( "http://x/y.png#frag" )
                || !TipImages.hasImageExtension( "http://x/y.jpg?v=2#frag" ) ) {
            return fail( "hasImageExtension: a trailing #fragment must be stripped too" );
        }
        if ( !TipImages.isUrl( "https://a" ) || !TipImages.isUrl( "HTTP://a" ) || TipImages.isUrl( "figs/a.png" )
                || TipImages.isUrl( "/abs/a.png" ) ) {
            return fail( "isUrl should be true only for http(s)" );
        }
        // hasTipImages over a tree
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( withProperty( "data:image", "a.png" ) );
        root.addAsChild( new PhylogenyNode() ); // no image
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        if ( !TipImages.hasTipImages( phy ) ) {
            return fail( "hasTipImages should be true when a tip has an image" );
        }
        final Phylogeny bare = new Phylogeny();
        final PhylogenyNode br = new PhylogenyNode();
        br.addAsChild( new PhylogenyNode() );
        bare.setRoot( br );
        bare.externalNodesHaveChanged();
        if ( TipImages.hasTipImages( bare ) ) {
            return fail( "hasTipImages should be false for a tree with no image tips" );
        }
        return true;
    }

    private static boolean testScaledSize() {
        // scale to a fixed HEIGHT, preserving aspect (a wide silhouette keeps its width)
        final int[] wide = TipImages.scaledSize( 300, 100, 40, 0 ); // 3:1 -> 120 x 40
        if ( ( wide[ 0 ] != 120 ) || ( wide[ 1 ] != 40 ) ) {
            return fail( "3:1 image scaled to height 40 should be 120x40, got " + wide[ 0 ] + "x" + wide[ 1 ] );
        }
        final int[] tall = TipImages.scaledSize( 100, 200, 40, 0 ); // 1:2 -> 20 x 40
        if ( ( tall[ 0 ] != 20 ) || ( tall[ 1 ] != 40 ) ) {
            return fail( "1:2 image scaled to height 40 should be 20x40, got " + tall[ 0 ] + "x" + tall[ 1 ] );
        }
        // clamp a very wide image to max_w (height scales down to keep aspect)
        final int[] clamped = TipImages.scaledSize( 1000, 100, 40, 200 ); // would be 400x40 -> clamp to 200x20
        if ( ( clamped[ 0 ] != 200 ) || ( clamped[ 1 ] != 20 ) ) {
            return fail( "a 10:1 image should clamp to 200x20, got " + clamped[ 0 ] + "x" + clamped[ 1 ] );
        }
        final int[] degen = TipImages.scaledSize( 0, 100, 40, 0 );
        if ( ( degen[ 0 ] != 0 ) || ( degen[ 1 ] != 0 ) ) {
            return fail( "a degenerate image should scale to 0x0" );
        }
        return true;
    }

    // ---- TipImageCache.loadImage: decode a real local PNG; missing / non-image -> null ----
    private static boolean testLoadLocal() {
        try {
            final File dir = Files.tempDir();
            final File png = new File( dir, "square.png" );
            writePng( png, 24, 12, Color.RED );
            // absolute path
            final BufferedImage abs = TipImageCache.loadImage( png.getAbsolutePath(), null );
            if ( ( abs == null ) || ( abs.getWidth() != 24 ) || ( abs.getHeight() != 12 ) ) {
                return fail( "loading an absolute PNG path should give a 24x12 image" );
            }
            // relative path resolved against base_dir
            final BufferedImage rel = TipImageCache.loadImage( "square.png", dir );
            if ( ( rel == null ) || ( rel.getWidth() != 24 ) ) {
                return fail( "a relative path should resolve against base_dir" );
            }
            // missing file
            if ( TipImageCache.loadImage( "does_not_exist.png", dir ) != null ) {
                return fail( "a missing file should load as null" );
            }
            // a file that is not a decodable image
            final File notimg = new File( dir, "notimg.png" );
            java.nio.file.Files.write( notimg.toPath(), "not an image".getBytes( "UTF-8" ) );
            if ( TipImageCache.loadImage( "notimg.png", dir ) != null ) {
                return fail( "a non-image file should load as null" );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "loadLocal threw: " + e.getMessage() );
        }
    }

    // ---- the async cache: get() is null first, then serves the image after the background load, and negative-caches ----
    private static boolean testAsyncCache() {
        try {
            final File dir = Files.tempDir();
            final File png = new File( dir, "c.png" );
            writePng( png, 10, 10, Color.BLUE );
            final TipImageCache cache = new TipImageCache( 8 );
            // first get schedules the load and returns null
            if ( cache.get( "c.png", dir ) != null ) {
                return fail( "the first get() should return null while the image loads" );
            }
            final BufferedImage img = waitFor( cache, "c.png", dir );
            if ( ( img == null ) || ( img.getWidth() != 10 ) ) {
                return fail( "the cache should serve the loaded image after the background load" );
            }
            // a broken reference is negative-cached (isFailed true after the load attempt)
            cache.get( "missing.png", dir );
            long end = System.currentTimeMillis() + 3000;
            while ( !cache.isFailed( "missing.png", dir ) && ( System.currentTimeMillis() < end ) ) {
                Thread.sleep( 20 );
            }
            if ( !cache.isFailed( "missing.png", dir ) ) {
                return fail( "a missing reference should be negative-cached (isFailed)" );
            }
            if ( cache.get( "missing.png", dir ) != null ) {
                return fail( "a negative-cached reference should keep returning null" );
            }
            cache.clear();
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "asyncCache threw: " + e.getMessage() );
        }
    }

    private static BufferedImage waitFor( final TipImageCache cache, final String ref, final File dir )
            throws InterruptedException {
        final long end = System.currentTimeMillis() + 3000;
        BufferedImage img = cache.get( ref, dir );
        while ( ( img == null ) && !cache.isFailed( ref, dir ) && ( System.currentTimeMillis() < end ) ) {
            Thread.sleep( 20 );
            img = cache.get( ref, dir );
        }
        return img;
    }

    // ---------------------------------------------------------------------------------------

    private static PhylogenyNode withProperty( final String ref, final String value ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "t" );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( ref, value, "", "xsd:string", AppliesTo.NODE ) );
        n.getNodeData().setProperties( pl );
        return n;
    }

    private static PhylogenyNode withTaxonomyUri( final String uri, final String type ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "t" );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( "Homo sapiens" );
        t.addUri( new Uri( uri, "", type ) );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static void writePng( final File f, final int w, final int h, final Color c ) throws Exception {
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final java.awt.Graphics2D g = img.createGraphics();
        g.setColor( c );
        g.fillRect( 0, 0, w, h );
        g.dispose();
        ImageIO.write( img, "png", f );
    }

    /** A per-run temp directory (deleted on JVM exit). */
    private static final class Files {

        static File tempDir() throws Exception {
            final File d = java.nio.file.Files.createTempDirectory( "aptx-tipimg" ).toFile();
            d.deleteOnExit();
            return d;
        }
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [TipImageTest] " + message );
        return false;
    }

    private TipImageTest() {
        // not instantiable
    }
}
