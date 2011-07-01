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
import java.util.ArrayList;
import java.util.List;

import org.forester.phylogeny.data.Property.AppliesTo;

public class NodeVisualization implements PhylogenyData {

    private NodeShape _shape;
    private NodeFill  _fill_type;
    private Color     _border_color;
    private Color     _fill_color;
    private double    _size;
    private double    _transparancy;

    public NodeVisualization() {
        init();
    }

    public NodeVisualization( final NodeShape shape,
                              final NodeFill fill_type,
                              final Color border_color,
                              final Color fill_color,
                              final double size,
                              final double transparancy ) {
        setShape( shape );
        setFillType( fill_type );
        setBorderColor( border_color );
        setFillColor( fill_color );
        setSize( size );
        setTransparancy( transparancy );
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
                                      getSize(),
                                      getTransparancy() );
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

    public double getTransparancy() {
        return _transparancy;
    }

    private void init() {
        setShape( NodeShape.CIRCLE );
        setFillType( NodeFill.SOLID );
        setBorderColor( null );
        setFillColor( null );
        setSize( 0 );
        setTransparancy( 1 );
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

    public void setTransparancy( final double transparancy ) {
        _transparancy = transparancy;
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

    public enum NodeFill {
        NONE, GRADIENT, SOLID
    }

    public enum NodeShape {
        CIRCLE, RECTANGLE
    }

    private List<Property> toProperties() {
        final List<Property> properties = new ArrayList<Property>();
        properties.add( new Property( SIZE_REF, String.valueOf( getSize() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getShape() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getFillType() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getTransparancy() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getFillColor() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        properties.add( new Property( SIZE_REF, String.valueOf( getBorderColor() ), "", SIZE_TYPE, AppliesTo.NODE ) );
        return properties;
    }
    public static final String SIZE_REF  = "aptx_visualiation:node_sise";
    public static final String SIZE_TYPE = "xsd:decimal";
}
