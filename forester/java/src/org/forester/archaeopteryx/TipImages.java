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

import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;
import org.forester.util.ForesterUtil;

/**
 * Pure helpers for the tip-images feature: resolve a tip's image REFERENCE (a local file path OR an http(s) URL)
 * from its node data, and small size math. The reference can live in a node PROPERTY (the common "import a table of
 * tip -> image path" case, fed by the Import Annotations tool) or in the phyloXML {@code <taxonomy><uri
 * type="image_url">} (the standard slot that round-trips in the tree file). The actual loading / caching / scaling
 * lives in {@link TipImageCache}; the rendering lives in {@code TreePanel}.
 * <p>
 * v1 handles RASTER images ImageIO can decode (PNG / JPG / GIF / BMP) -- PNG-with-transparency suits silhouettes;
 * SVG is not supported (it needs a rasterizer dependency) and is a deliberate follow-up.
 */
final class TipImages {

    /** Property refs (matched on the namespace OR the local name, case-insensitive) recognized as a tip-image
     *  reference -- so an "image", "img:src", "data:photo", ... column becomes a tip image. */
    private static final Set<String> IMAGE_KEYS = new HashSet<String>( Arrays.asList( "image", "img", "picture",
            "photo", "thumbnail", "image_url", "tip_image", "silhouette" ) );

    /** Raster image file extensions ImageIO can decode (SVG is intentionally absent -- v1 has no rasterizer). */
    private static final Set<String> IMAGE_EXTS = new HashSet<String>( Arrays.asList( "png", "jpg", "jpeg", "gif",
            "bmp" ) );

    private TipImages() {
        // not instantiable
    }

    /** The image reference (local path or URL) for {@code node}, or null: the first recognized image PROPERTY, else
     *  an image {@code <uri>} on the node's taxonomy. */
    static String imageRefFor( final PhylogenyNode node ) {
        if ( ( node == null ) || ( node.getNodeData() == null ) ) {
            return null;
        }
        if ( node.getNodeData().getProperties() != null ) {
            for( final Property p : node.getNodeData().getProperties().getProperties() ) {
                if ( ( p != null ) && !ForesterUtil.isEmpty( p.getValue() ) && isImageRef( p.getRef(), p.getValue() ) ) {
                    return p.getValue().trim();
                }
            }
        }
        if ( node.getNodeData().isHasTaxonomy() ) {
            final Taxonomy t = node.getNodeData().getTaxonomy();
            if ( t.getUris() != null ) {
                for( final Uri u : t.getUris() ) {
                    if ( ( u != null ) && ( u.getValue() != null ) && isImageUri( u ) ) {
                        return u.getValue().toString().trim();
                    }
                }
            }
        }
        return null;
    }

    /** Whether a property (ref + value) denotes a tip image: the ref's namespace or local name is a recognized image
     *  key, OR the value itself ends with a known raster-image extension. */
    static boolean isImageRef( final String ref, final String value ) {
        if ( !ForesterUtil.isEmpty( ref ) ) {
            final String lower = ref.toLowerCase( Locale.ROOT );
            final int colon = lower.lastIndexOf( ':' );
            final String local = ( colon >= 0 ) ? lower.substring( colon + 1 ) : lower;
            final String ns = ( colon >= 0 ) ? lower.substring( 0, colon ) : "";
            if ( IMAGE_KEYS.contains( local ) || ( !ns.isEmpty() && IMAGE_KEYS.contains( ns ) ) ) {
                return true;
            }
        }
        return hasImageExtension( value );
    }

    private static boolean isImageUri( final Uri u ) {
        final String type = ( u.getType() == null ) ? "" : u.getType().toLowerCase( Locale.ROOT );
        return type.contains( "image" ) || hasImageExtension( u.getValue().toString() );
    }

    /** Whether {@code s} ends with a raster-image extension (a trailing {@code ?query} or {@code #fragment} on a URL is
     *  ignored). */
    static boolean hasImageExtension( final String s ) {
        if ( ForesterUtil.isEmpty( s ) ) {
            return false;
        }
        String path = s.trim().toLowerCase( Locale.ROOT );
        final int q = path.indexOf( '?' );
        if ( q >= 0 ) {
            path = path.substring( 0, q );
        }
        final int hash = path.indexOf( '#' );
        if ( hash >= 0 ) {
            path = path.substring( 0, hash );
        }
        final int dot = path.lastIndexOf( '.' );
        return ( dot >= 0 ) && IMAGE_EXTS.contains( path.substring( dot + 1 ) );
    }

    /** Whether {@code ref} is an http(s) URL (else a local file path). */
    static boolean isUrl( final String ref ) {
        if ( ForesterUtil.isEmpty( ref ) ) {
            return false;
        }
        final String l = ref.trim().toLowerCase( Locale.ROOT );
        return l.startsWith( "http://" ) || l.startsWith( "https://" );
    }

    /** True iff any external tip carries an image reference -- drives the on-load auto-enable. */
    static boolean hasTipImages( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return false;
        }
        for( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( imageRefFor( ext ) != null ) {
                return true;
            }
        }
        return false;
    }

    /**
     * The drawn {@code [width, height]} for an image scaled so its HEIGHT is {@code target_h}, preserving aspect --
     * so a wide silhouette keeps its full width. A very wide image is clamped to {@code max_w} (its height then scaled
     * down to keep aspect); {@code max_w <= 0} means no clamp. Degenerate inputs give {@code [0, 0]}.
     */
    static int[] scaledSize( final int img_w, final int img_h, final int target_h, final int max_w ) {
        if ( ( img_w <= 0 ) || ( img_h <= 0 ) || ( target_h <= 0 ) ) {
            return new int[] { 0, 0 };
        }
        int h = target_h;
        int w = Math.max( 1, Math.round( target_h * ( (float) img_w / (float) img_h ) ) );
        if ( ( max_w > 0 ) && ( w > max_w ) ) {
            w = max_w;
            h = Math.max( 1, Math.round( max_w * ( (float) img_h / (float) img_w ) ) );
        }
        return new int[] { w, h };
    }
}
