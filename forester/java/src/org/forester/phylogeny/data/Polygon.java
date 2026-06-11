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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public class Polygon implements PhylogenyData {

    private final List<Point> _points;

    public Polygon( final List<Point> points ) {
        _points = points;
    }

    @Override
    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        boolean first = true;
        for( final Point point : getPoints() ) {
            if ( first ) {
                first = false;
            }
            else {
                sb.append( ForesterUtil.LINE_SEPARATOR );
            }
            sb.append( point.asSimpleText() );
        }
        return sb;
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        final List<Point> new_points = new ArrayList<Point>();
        for( final Point point : getPoints() ) {
            new_points.add( ( Point ) point.copy() );
        }
        return new Polygon( new_points );
    }

    public List<Point> getPoints() {
        return _points;
    }

    public boolean isEmpty() {
        return ForesterUtil.isEmpty( _points );
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
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.POLYGON );
        for( final Point point : getPoints() ) {
            point.toPhyloXML( writer, level, indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
            writer.write( indentation );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.POLYGON );
    }
}
