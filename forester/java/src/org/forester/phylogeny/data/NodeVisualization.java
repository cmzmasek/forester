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
// WWW: www.phylosoft.org

package org.forester.phylogeny.data;

import java.awt.Color;
import java.io.IOException;
import java.io.Writer;

public class NodeVisualization implements PhylogenyData {

    public enum NodeFill {
        NONE, GRADIENT, SOLID
    }

    public enum NodeShape {
        CIRCLE, RECTANGLE
    }
    private NodeShape _shape;
    private NodeFill  _fill_type;
    private Color     _border_color;
    private Color     _fill_color;
    private double    _size;

    public NodeVisualization() {
        _shape = NodeShape.CIRCLE;
        _fill_type = NodeFill.SOLID;
        _border_color = null;
        _fill_color = null;
        _size = 0;
    }

    public NodeVisualization( final NodeShape shape,
                              final NodeFill fill_type,
                              final Color border_color,
                              final Color fill_color,
                              final double size ) {
        _shape = shape;
        _fill_type = fill_type;
        _border_color = border_color;
        _fill_color = fill_color;
        _size = size;
    }

    @Override
    public StringBuffer asSimpleText() {
        return asText();
    }

    @Override
    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        return sb;
    }

    @Override
    public PhylogenyData copy() {
        return new NodeVisualization( getShape(),
                                      getFillType(),
                                      getBorderColor() != null ? new Color( getBorderColor().getRed(), getBorderColor()
                                              .getGreen(), getBorderColor().getBlue() ) : null,
                                      getFillColor() != null ? new Color( getFillColor().getRed(), getFillColor()
                                              .getGreen(), getFillColor().getBlue() ) : null,
                                      getSize() );
    }

    public Color getBorderColor() {
        return _border_color;
    }

    public Color getFillColor() {
        return _fill_color;
    }

    public NodeFill getFillType() {
        return _fill_type;
    }

    public NodeShape getShape() {
        return _shape;
    }

    public double getSize() {
        return _size;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    public void setBorderColor( final Color border_color ) {
        _border_color = border_color;
    }

    public void setFillColor( final Color fill_color ) {
        _fill_color = fill_color;
    }

    public void setFillType( final NodeFill fill_type ) {
        _fill_type = fill_type;
    }

    public void setShape( final NodeShape shape ) {
        _shape = shape;
    }

    public void setSize( final double size ) {
        _size = size;
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
