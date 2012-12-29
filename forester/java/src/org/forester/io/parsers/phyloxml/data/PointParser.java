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

import java.math.BigDecimal;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Point;
import org.forester.util.ForesterUtil;

public class PointParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new PointParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private PointParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        String alt_unit = "";
        String geo_datum = "";
        if ( element.isHasAttribute( PhyloXmlMapping.POINT_ALTITUDE_UNIT_ATTR ) ) {
            alt_unit = element.getAttribute( PhyloXmlMapping.POINT_ALTITUDE_UNIT_ATTR );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.POINT_GEODETIC_DATUM ) ) {
            geo_datum = element.getAttribute( PhyloXmlMapping.POINT_GEODETIC_DATUM );
        }
        String lat_str = null;
        String lon_str = null;
        String alt_str = null;
        for( int j = 0; j < element.getNumberOfChildElements(); ++j ) {
            final XmlElement e = element.getChildElement( j );
            if ( e.getQualifiedName().equals( PhyloXmlMapping.POINT_LATITUDE ) ) {
                lat_str = e.getValueAsString();
            }
            else if ( e.getQualifiedName().equals( PhyloXmlMapping.POINT_LONGITUDE ) ) {
                lon_str = e.getValueAsString();
            }
            else if ( e.getQualifiedName().equals( PhyloXmlMapping.POINT_ALTITUDE ) ) {
                alt_str = e.getValueAsString();
            }
        }
        BigDecimal lat = null;
        BigDecimal lon = null;
        BigDecimal alt = null;
        if ( !ForesterUtil.isEmpty( lat_str ) ) {
            lat = new BigDecimal( lat_str );
        }
        if ( !ForesterUtil.isEmpty( lon_str ) ) {
            lon = new BigDecimal( lon_str );
        }
        if ( !ForesterUtil.isEmpty( alt_str ) ) {
            alt = new BigDecimal( alt_str );
        }
        return new Point( geo_datum, lat, lon, alt, alt_unit );
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
