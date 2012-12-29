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

package org.forester.io.parsers.phyloxml.data;

import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Point;
import org.forester.phylogeny.data.Polygon;

public class PolygonParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new PolygonParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private PolygonParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        final List<Point> points = new ArrayList<Point>();
        for( int j = 0; j < element.getNumberOfChildElements(); ++j ) {
            final XmlElement e = element.getChildElement( j );
            if ( e.getQualifiedName().equals( PhyloXmlMapping.POINT ) ) {
                points.add( ( Point ) PointParser.getInstance().parse( e ) );
            }
        }
        return new Polygon( points );
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
