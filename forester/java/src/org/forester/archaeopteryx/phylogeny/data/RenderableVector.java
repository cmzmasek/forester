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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx.phylogeny.data;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

import org.forester.archaeopteryx.AptxUtil;
import org.forester.archaeopteryx.Configuration;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public final class RenderableVector implements RenderablePhylogenyData {

    final static int                VECTOR_DEFAULT_HEIGHT   = 12;
    public final static int         VECTOR_DEFAULT_WIDTH    = 120;
    private double                  _rendering_factor_width = 1.0;
    private List<Double>            _values;
    private final Rectangle2D       _rectangle              = new Rectangle2D.Float();
    private double                  _height                 = VECTOR_DEFAULT_HEIGHT;
    private double                  _min;
    private double                  _max;
    private double                  _mean;
    private Color                   _min_color              = Color.BLUE;
    private Color                   _max_color              = Color.YELLOW;
    private Color                   _mean_color             = Color.WHITE;
    private int                     _width                  = VECTOR_DEFAULT_WIDTH;
    private static RenderableVector _instance               = null;

    private RenderableVector() {
        _values = null;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( _values.toString() );
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public Object clone() {
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

    public double getRenderingFactorWidth() {
        return _rendering_factor_width;
    }

    @Override
    public Dimension getRenderingSize() {
        return getOriginalSize();
    }

    public int getTotalLength() {
        return ( int ) ( _values.size() * getRenderingHeight() );
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new NoSuchMethodError();
    }

    @Override
    public void render( final float x1,
                        final float y1,
                        final Graphics2D g,
                        final TreePanel tree_panel,
                        final boolean to_pdf ) {
        final double y = y1;
        final double start = x1 + 20.0;
        final double width = ( double ) _width / _values.size();
        for( int i = 0; i < _values.size(); ++i ) {
            g.setColor( calculateColor( _values.get( i ) ) );
            _rectangle.setFrame( start + ( i * width ), y - 0.5, width + 1, getRenderingHeight() );
            g.fill( _rectangle );
        }
    }

    @Override
    public void setParameter( final double parameter ) {
        throw new NoSuchMethodError();
    }

    public void setRenderingFactorWidth( final double rendering_factor_width ) {
        _rendering_factor_width = rendering_factor_width;
    }

    @Override
    public void setRenderingHeight( final float height ) {
        _height = height;
    }

    @Override
    public StringBuffer toNHX() {
        throw new NoSuchMethodError();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        throw new NoSuchMethodError();
    }

    private Color calculateColor( final double v ) {
        return ForesterUtil.calcColor( v, _min, _max, _mean, _min_color, _max_color, _mean_color );
    }

    private double getRenderingHeight() {
        return _height;
    }

    public static RenderableVector createInstance( final List<Double> values,
                                                   final DescriptiveStatistics stats,
                                                   final Configuration configuration ) {
        if ( _instance == null ) {
            _instance = new RenderableVector();
        }
        _instance._values = values;
        if ( configuration != null ) {
            _instance._min_color = configuration.getVectorDataMinColor();
            _instance._max_color = configuration.getVectorDataMaxColor();
            _instance._mean_color = configuration.getVectorDataMeanColor();
            _instance._width = configuration.getVectorDataWidth();
            _instance._height = configuration.getVectorDataHeight();
        }
        if ( stats.getN() > 0 ) {
            _instance._min = stats.getMin();
            _instance._max = stats.getMax();
            _instance._mean = stats.arithmeticMean();
        }
        else {
            _instance._min = 0;
            _instance._max = 0;
            _instance._mean = 0;
            AptxUtil.printWarningMessage( "Archaeopteryx", "creating renderable vector with empty statistics" );
        }
        return _instance;
    }
}
