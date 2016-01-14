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

import org.forester.archaeopteryx.Configuration;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.sequence.MolecularSequence;
import org.forester.sequence.MolecularSequence.TYPE;

public final class RenderableMsaSequence implements RenderablePhylogenyData {

    final static int                     DEFAULT_HEIGHT          = 12;
    final public static int              DEFAULT_WIDTH           = 400;
    private double                       _rendering_factor_width = 1.0;
    private char                         _seq[];
    private final Rectangle2D            _rectangle              = new Rectangle2D.Float();
    private double                       _height                 = DEFAULT_HEIGHT;
    private final float                  _width                  = DEFAULT_WIDTH;
    private MolecularSequence.TYPE       _type;
    private static RenderableMsaSequence _instance               = null;

    private RenderableMsaSequence() {
        _seq = null;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( _seq.toString() );
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
        return _seq.length;
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
        final float y = y1;
        final float start = x1 + 20;
        final float width = _width / _seq.length;
        for( int i = 0; i < _seq.length; ++i ) {
            final char c = _seq[ i ];
            if ( width < 4 ) {
                if ( c != '-' ) {
                    g.setColor( calculateColor( c ) );
                    _rectangle.setFrame( start + ( i * width ), y - 0.5, width + 1, getRenderingHeight() );
                    g.fill( _rectangle );
                }
            }
            else {
                g.setColor( calculateColor( c ) );
                g.drawString( String.valueOf( c ), start + ( i * width ), y - 0.5f );
            }
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

    private Color calculateColor( final char c ) {
        if ( _type == TYPE.AA ) {
            return calculateAAColor( c );
        }
        return calculateNucleotideColor( c );
    }

    private Color calculateNucleotideColor( final char c ) {
        if ( c == 'A' ) {
            return Color.YELLOW;
        }
        if ( ( c == 'T' ) || ( c == 'U' ) ) {
            return Color.ORANGE;
        }
        if ( c == 'G' ) {
            return Color.BLUE;
        }
        if ( c == 'C' ) {
            return Color.CYAN;
        }
        else if ( c == '-' ) {
            return Color.GRAY;
        }
        else {
            return Color.GRAY;
        }
    }

    private Color calculateAAColor( final char c ) {
        if ( ( c == 'G' ) || ( c == 'A' ) || ( c == 'S' ) || ( c == 'T' ) ) {
            return Color.YELLOW;
        }
        else if ( ( c == 'N' ) || ( c == 'Q' ) || ( c == 'H' ) ) {
            return Color.PINK;
        }
        else if ( ( c == 'D' ) || ( c == 'E' ) ) {
            return Color.RED;
        }
        else if ( ( c == 'K' ) || ( c == 'R' ) ) {
            return Color.BLUE;
        }
        else if ( c == '-' ) {
            return Color.GRAY;
        }
        else if ( c == 'X' ) {
            return Color.GRAY;
        }
        else {
            return Color.GREEN;
        }
    }

    private double getRenderingHeight() {
        return _height;
    }

    public static RenderableMsaSequence createInstance( final String seq,
                                                        final String type,
                                                        final Configuration configuration ) {
        if ( _instance == null ) {
            _instance = new RenderableMsaSequence();
        }
        if ( type.equals( "protein" ) ) {
            _instance._type = TYPE.AA;
        }
        else if ( type.equals( "dna" ) ) {
            _instance._type = TYPE.DNA;
        }
        else {
            _instance._type = TYPE.RNA;
        }
        _instance._seq = seq.toCharArray();
        if ( configuration != null ) {
        }
        return _instance;
    }
}
