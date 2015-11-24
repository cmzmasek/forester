// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public class Distribution implements PhylogenyData {

    private final String        _desc;
    private final List<Point>   _points;
    private final List<Polygon> _polygons;

    public Distribution( final String desc ) {
        _desc = desc;
        _points = null;
        _polygons = null;
    }

    public Distribution( final String desc, final List<Point> points ) {
        _desc = null;
        _points = points;
        _polygons = null;
    }

    public Distribution( final String desc, final List<Point> points, final List<Polygon> polygons ) {
        _desc = desc;
        _points = points;
        _polygons = polygons;
    }

    public boolean isEmpty() {
        if ( ForesterUtil.isEmpty( _desc ) && ( ( getPoints() != null ) && ( getPoints().size() == 1 ) )
                && ForesterUtil.isEmpty( _polygons ) ) {
            if ( Point.isSeemsEmpty( getPoints().get( 0 ) ) ) {
                return true;
            }
        }
        return ForesterUtil.isEmpty( _desc ) && ForesterUtil.isEmpty( _points ) && ForesterUtil.isEmpty( _polygons );
    }

    @Override
    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        if ( isEmpty() ) {
            return sb;
        }
        sb.append( "Distribution: " );
        if ( !ForesterUtil.isEmpty( getDesc() ) ) {
            sb.append( ForesterUtil.LINE_SEPARATOR );
            sb.append( " Description: " );
            sb.append( getDesc() );
        }
        int i = 0;
        if ( getPoints() != null ) {
            for( final Point point : getPoints() ) {
                if ( ( point != null ) && !Point.isSeemsEmpty( point ) ) {
                    sb.append( ForesterUtil.LINE_SEPARATOR );
                    sb.append( " Point " + i + ": " );
                    sb.append( point.asSimpleText() );
                    i++;
                }
            }
        }
        i = 0;
        if ( getPolygons() != null ) {
            for( final Polygon polygon : getPolygons() ) {
                if ( polygon != null ) {
                    sb.append( ForesterUtil.LINE_SEPARATOR );
                    sb.append( " Polygon " + i + ":" );
                    sb.append( ForesterUtil.LINE_SEPARATOR );
                    sb.append( polygon.asSimpleText() );
                    i++;
                }
            }
        }
        return sb;
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        List<Point> new_points = null;
        List<Polygon> new_polygons = null;
        if ( getPoints() != null ) {
            new_points = new ArrayList<Point>();
            for( final Point point : getPoints() ) {
                new_points.add( ( Point ) point.copy() );
            }
        }
        if ( getPolygons() != null ) {
            new_polygons = new ArrayList<Polygon>();
            for( final Polygon polygon : getPolygons() ) {
                new_polygons.add( ( Polygon ) polygon.copy() );
            }
        }
        return new Distribution( getDesc(), new_points, new_polygons );
    }

    public String getDesc() {
        return _desc;
    }

    public List<Point> getPoints() {
        return _points;
    }

    public List<Polygon> getPolygons() {
        return _polygons;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( isEmpty() ) {
            return;
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.DISTRIBUTION );
        if ( !ForesterUtil.isEmpty( getDesc() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.DISTRIBUTION_DESC, getDesc(), indentation );
        }
        final String ind = indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE;
        if ( getPoints() != null ) {
            for( final Point point : getPoints() ) {
                if ( ( point != null ) && !Point.isSeemsEmpty( point ) ) {
                    point.toPhyloXML( writer, level, ind );
                }
            }
        }
        if ( getPolygons() != null ) {
            for( final Polygon polygon : getPolygons() ) {
                if ( polygon != null ) {
                    polygon.toPhyloXML( writer, level, ind );
                }
            }
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.DISTRIBUTION );
    }

    @Override
    public String toString() {
        return asSimpleText().toString();
    }
}
