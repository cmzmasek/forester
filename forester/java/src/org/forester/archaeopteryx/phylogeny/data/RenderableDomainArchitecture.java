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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.io.Writer;
import java.math.BigDecimal;
import java.util.Map;
import java.util.SortedMap;

import org.forester.archaeopteryx.AptxUtil;
import org.forester.archaeopteryx.Constants;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.util.ForesterUtil;

public final class RenderableDomainArchitecture extends DomainArchitecture implements RenderablePhylogenyData {

    final static private int          BRIGHTEN_COLOR_BY             = 200;
    final static private int          E_VALUE_THRESHOLD_EXP_DEFAULT = 0;
    final static private BasicStroke  STROKE_1                      = new BasicStroke( 1f );
    private static Map<String, Color> _domain_colors;
    private final DomainArchitecture  _domain_structure;
    private int                       _e_value_threshold_exp        = E_VALUE_THRESHOLD_EXP_DEFAULT;
    private final Rectangle2D         _rectangle                    = new Rectangle2D.Float();
    private float                    _rendering_factor_width       = 1;
    private float                    _rendering_height             = 0;

    public RenderableDomainArchitecture( final DomainArchitecture domain_structure ) {
        _domain_structure = domain_structure;
    }

    public static void setColorMap( final Map<String, Color> domain_colors ) {
        _domain_colors = domain_colors;
    }

    @Override
    public StringBuffer asSimpleText() {
        return _domain_structure.asSimpleText();
    }

    @Override
    public StringBuffer asText() {
        return _domain_structure.asText();
    }

    @Override
    public PhylogenyData copy() {
        return _domain_structure.copy();
    }

    private final void drawDomain( final double x,
                                   final double y,
                                   final double width,
                                   final double heigth,
                                   final String name,
                                   final Graphics2D g,
                                   final boolean to_pdf ) {
        final double h2 = heigth / 2.0;
        final Color color_one = getColorOne( name );
        final Color color_two = getColorTwo( color_one );
        double step = 1;
        if ( to_pdf ) {
            step = 0.05;
        }
        for( double i = 0; i < heigth; i += step ) {
            g.setColor( org.forester.util.ForesterUtil
                    .calcColor( i >= h2 ? heigth - i : i, 0, h2, color_one, color_two ) );
            _rectangle.setFrame( x, i + y, width, step );
            g.fill( _rectangle );
        }
    }

    private final Color getColorOne( final String name ) {
        Color c = _domain_colors.get( name );
        if ( c == null ) {
            c = AptxUtil.calculateColorFromString( name, false );
            if ( c == null ) {
                throw new IllegalStateException();
            }
            _domain_colors.put( name, c );
        }
        return c;
    }

    private Color getColorTwo( final Color color_one ) {
        final int red = color_one.getRed() + RenderableDomainArchitecture.BRIGHTEN_COLOR_BY;
        final int green = color_one.getGreen() + RenderableDomainArchitecture.BRIGHTEN_COLOR_BY;
        final int blue = color_one.getBlue() + RenderableDomainArchitecture.BRIGHTEN_COLOR_BY;
        return new Color( red > 255 ? 255 : red, green > 255 ? 255 : green, blue > 255 ? 255 : blue );
    }

    @Override
    public ProteinDomain getDomain( final int i ) {
        return _domain_structure.getDomain( i );
    }

    @Override
    public SortedMap<BigDecimal, ProteinDomain> getDomains() {
        return _domain_structure.getDomains();
    }

    @Override
    public int getNumberOfDomains() {
        return _domain_structure.getNumberOfDomains();
    }

    @Override
    public Dimension getOriginalSize() {
        return new Dimension( _domain_structure.getTotalLength(), ForesterUtil.roundToInt( _rendering_height ) );
    }

    @Override
    public Object getParameter() {
        return new Integer( _e_value_threshold_exp );
    }

    public float getRenderingFactorWidth() {
        return _rendering_factor_width;
    }

    @Override
    public Dimension getRenderingSize() {
        return new Dimension( ForesterUtil.roundToInt( _domain_structure.getTotalLength() * getRenderingFactorWidth() ),
                              ForesterUtil.roundToInt( _rendering_height ) );
    }

    @Override
    public int getTotalLength() {
        return _domain_structure.getTotalLength();
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        return _domain_structure.isEqual( data );
    }

    @Override
    public void render( final float x1,
                        final float y1,
                        final Graphics2D g,
                        final TreePanel tree_panel,
                        final boolean to_pdf ) {
        final float f = getRenderingFactorWidth();
        final float y = y1 + ( _rendering_height / 2 );
        final float start = x1 + 20;
        final Stroke s = g.getStroke();
        g.setStroke( STROKE_1 );
        if ( !to_pdf ) {
            g.setColor( tree_panel.getTreeColorSet().getDomainBaseColor() );
        }
        else {
            g.setColor( Constants.DOMAIN_BASE_COLOR_FOR_PDF );
        }
        _rectangle.setFrame( start, y - 0.5, _domain_structure.getTotalLength() * f, 1 );
        g.fill( _rectangle );
        for( int i = 0; i < _domain_structure.getDomains().size(); ++i ) {
            final ProteinDomain d = _domain_structure.getDomain( i );
            if ( d.getConfidence() <= Math.pow( 10, _e_value_threshold_exp ) ) {
                final float xa = start + ( d.getFrom() * f );
                final float xb = xa + ( d.getLength() * f );
                if ( tree_panel.getMainPanel().getOptions().isShowDomainLabels()
                        && ( tree_panel.getMainPanel().getTreeFontSet().getFontMetricsSmall().getHeight() > 4 ) ) {
                    g.setFont( tree_panel.getMainPanel().getTreeFontSet().getSmallFont() );
                    if ( !to_pdf ) {
                        g.setColor( tree_panel.getTreeColorSet().getDomainLabelColor() );
                    }
                    else {
                        g.setColor( Constants.DOMAIN_LABEL_COLOR_FOR_PDF );
                    }
                    g.drawString( d.getName(),
                                  xa,
                                  y1
                                  + tree_panel.getMainPanel().getTreeFontSet().getFontMetricsSmall().getAscent()
                                  + _rendering_height  );
                    
                    
                }
                drawDomain( xa, y1, xb - xa, _rendering_height, d.getName(), g, to_pdf );
            }
        }
        g.setStroke( s );
    }

    @Override
    public void setParameter( final double e_value_threshold_exp ) {
        _e_value_threshold_exp = ( int ) e_value_threshold_exp;
    }

    public void setRenderingFactorWidth( final float rendering_factor_width ) {
        _rendering_factor_width = rendering_factor_width;
    }

    @Override
    public void setRenderingHeight( final float rendering_height ) {
        _rendering_height = rendering_height;
    }

    @Override
    public StringBuffer toNHX() {
        return _domain_structure.toNHX();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        _domain_structure.toPhyloXML( writer, level, indentation );
    }
}
