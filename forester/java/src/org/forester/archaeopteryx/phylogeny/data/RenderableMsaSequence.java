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
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public final class RenderableMsaSequence implements RenderablePhylogenyData {

    final static int                VECTOR_DEFAULT_HEIGHT   = 12;
    public final static int         VECTOR_DEFAULT_WIDTH    = 120;
    private double                  _rendering_factor_width = 1.0;
    private String            _seq;
    private final Rectangle2D       _rectangle              = new Rectangle2D.Float();
    private double                  _height                 = VECTOR_DEFAULT_HEIGHT;
    private double                  _min;
    private double                  _max;
    private double                  _mean;
    private Color                   _min_color              = Color.BLUE;
    private Color                   _max_color              = Color.YELLOW;
    private Color                   _mean_color             = Color.WHITE;
    private int                     _width                  = VECTOR_DEFAULT_WIDTH;
    private static RenderableMsaSequence _instance               = null;

    private RenderableMsaSequence() {
        _seq = null;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( _seq );
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
        return _seq.length();
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
        g.drawString( _seq, x1, y1
                                       );
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

    public static RenderableMsaSequence createInstance( final String seq,
                                                   final Configuration configuration ) {
        if ( _instance == null ) {
            _instance = new RenderableMsaSequence();
        }
        _instance._seq = seq;
        if ( configuration != null ) {
          
        }
       
        return _instance;
    }
}
