// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Point;
import org.forester.phylogeny.data.Polygon;

public class DistributionParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new DistributionParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private DistributionParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        String desc = "";
        List<Point> points = null;
        List<Polygon> polygons = null;
        for( int i = 0; i < element.getNumberOfChildElements(); ++i ) {
            final XmlElement child_element = element.getChildElement( i );
            if ( child_element.getQualifiedName().equals( PhyloXmlMapping.DISTRIBUTION_DESC ) ) {
                desc = child_element.getValueAsString();
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.POINT ) ) {
                if ( points == null ) {
                    points = new ArrayList<Point>();
                }
                points.add( ( Point ) PointParser.getInstance().parse( child_element ) );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.POLYGON ) ) {
                if ( polygons == null ) {
                    polygons = new ArrayList<Polygon>();
                }
                polygons.add( ( Polygon ) PolygonParser.getInstance().parse( child_element ) );
            }
        }
        return new Distribution( desc, points, polygons );
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
