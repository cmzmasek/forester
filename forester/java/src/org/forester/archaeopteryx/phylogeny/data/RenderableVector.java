// $Id:
// $
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
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx.phylogeny.data;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

import org.forester.archaeopteryx.Configuration;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public final class RenderableVector implements RenderablePhylogenyData {

    final static public int         DEFAULT_HEIGHT          = 12;
    final static public int         DEFAULT_WIDTH           = 120;
    private double                  _rendering_factor_width = 1.0;
    private List<Double>            _values;
    private final Rectangle2D       _rectangle              = new Rectangle2D.Float();
    private Configuration           _configuration;
    private double                  _height;
    private double                  _min;
    private double                  _max;
    private double                  _mean;
    private static RenderableVector _instance               = null;

    public static RenderableVector createInstance( final List<Double> values,
                                                   final DescriptiveStatistics stats,
                                                   final Configuration configuration ) {
        if ( _instance == null ) {
            _instance = new RenderableVector();
        }
        _instance.setRenderingHeight( DEFAULT_HEIGHT );
        _instance._values = values;
        _instance._configuration = configuration;
        _instance._min = stats.getMin();
        _instance._max = stats.getMax();
        _instance._mean = stats.arithmeticMean();
        return _instance;
    }

    @Override
    public Object clone() {
        throw new NoSuchMethodError();
    }

    private RenderableVector() {
        _values = null;
        _configuration = null;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( _values.toString() );
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    public double getRenderingFactorWidth() {
        return _rendering_factor_width;
    }

    public int getTotalLength() {
        return ( int ) ( _values.size() * getRenderingHeight() );
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new NoSuchMethodError();
    }

    @Override
    public void render( final double x1,
                        final double y1,
                        final Graphics2D g,
                        final TreePanel tree_panel,
                        final boolean to_pdf ) {
        final double y = y1;
        final double start = x1 + 20.0;
        final double width = ( double ) DEFAULT_WIDTH / _values.size();
        for( int i = 0; i < _values.size(); ++i ) {
            g.setColor( calculateColor( _values.get( i ) ) );
            _rectangle.setFrame( start + ( i * width ), y - 0.5, width, getRenderingHeight() );
            g.fill( _rectangle );
        }
    }

    private Color calculateColor( final double v ) {
        return ForesterUtil.calcColor( v, _min, _max, _mean, Color.MAGENTA, Color.GREEN, Color.WHITE );
    }

    public void setRenderingFactorWidth( final double rendering_factor_width ) {
        _rendering_factor_width = rendering_factor_width;
    }

    @Override
    public StringBuffer toNHX() {
        throw new NoSuchMethodError();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        throw new NoSuchMethodError();
    }

    @Override
    public PhylogenyData copy() {
        throw new NoSuchMethodError();
    }

    @Override
    public Dimension getOriginalSize() {
        return new Dimension( getTotalLength(), ( int ) getRenderingHeight() );
    }

    @Override
    public Object getParameter() {
        return null;
    }

    @Override
    public Dimension getRenderingSize() {
        return getOriginalSize();
    }

    @Override
    public void setParameter( final double parameter ) {
        throw new NoSuchMethodError();
    }

    @Override
    public void setRenderingHeight( final double height ) {
        _height = height;
    }

    private double getRenderingHeight() {
        return _height;
    }
}
