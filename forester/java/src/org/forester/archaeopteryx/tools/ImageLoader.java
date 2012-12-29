// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx.tools;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

import javax.imageio.ImageIO;

import org.forester.archaeopteryx.AptxUtil;
import org.forester.archaeopteryx.Constants;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class ImageLoader implements Runnable {

    private final TreePanel            _tp;
    private static final BufferedImage PLACEHOLDER = new BufferedImage( 1, 1, BufferedImage.TYPE_INT_RGB );
    private final static boolean       DEBUG       = false;

    public ImageLoader( final TreePanel tp ) {
        _tp = tp;
    }

    private void load() {
        Hashtable<String, BufferedImage> image_map = null;
        if ( _tp.getImageMap() != null ) {
            image_map = _tp.getImageMap();
        }
        else {
            image_map = new Hashtable<String, BufferedImage>();
            _tp.setImageMap( image_map );
        }
        // ImageIO.setUseCache( false );
        for( final PhylogenyNodeIterator it = _tp.getPhylogeny().iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( node.getNodeData().isHasTaxonomy() && ( node.getNodeData().getTaxonomy().getUris() != null )
                    && !node.getNodeData().getTaxonomy().getUris().isEmpty() ) {
                final List<Uri> us = new ArrayList<Uri>();
                for( final Taxonomy t : node.getNodeData().getTaxonomies() ) {
                    for( final Uri uri : t.getUris() ) {
                        us.add( uri );
                    }
                }
                for( final Uri uri : us ) {
                    if ( uri != null ) {
                        final String type = uri.getType().toLowerCase();
                        final String uri_str = uri.getValue().toString().toLowerCase();
                        if ( ( !image_map.containsKey( uri_str ) )
                                && ( type.equals( "image" ) || type.equals( "img" ) || type.equals( "photo" )
                                        || type.equals( "picture" ) || uri_str.endsWith( ".jpg" )
                                        || uri_str.endsWith( ".jpeg" ) || uri_str.endsWith( ".png" )
                                        || uri_str.endsWith( ".gif" ) || uri_str.endsWith( ".bmp" ) ) ) {
                            image_map.put( uri_str, PLACEHOLDER );
                            BufferedImage bi = null;
                            if ( DEBUG ) {
                                System.out.println( "accessing: " + uri );
                            }
                            try {
                                bi = ImageIO.read( uri.getValue().toURL() );
                            }
                            catch ( final MalformedURLException e ) {
                                AptxUtil.printWarningMessage( Constants.PRG_NAME, "\"" + uri.getValue()
                                        + "\": Malformed URL Exception: " + e.getLocalizedMessage() );
                            }
                            catch ( final IOException e ) {
                                AptxUtil.printWarningMessage( Constants.PRG_NAME, "\"" + uri.getValue()
                                        + "\": IO Exception: " + e.getLocalizedMessage() );
                            }
                            if ( bi != null ) {
                                image_map.put( uri_str, bi );
                                _tp.repaint();
                            }
                            else {
                                image_map.remove( uri_str );
                            }
                        }
                    }
                }
            }
        }
    }

    @Override
    public void run() {
        load();
    }
}
