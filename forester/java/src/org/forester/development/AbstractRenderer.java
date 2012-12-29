// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

package org.forester.development;

import java.awt.Color;
import java.awt.Graphics;

import javax.swing.JComponent;

public abstract class AbstractRenderer extends JComponent {

    /**
     * 
     */
    private static final long serialVersionUID   = 7236434322552764776L;
    static final Color        DEFAULT_COLOR      = new Color( 0, 0, 0 );
    static final Color        MARKED_COLOR       = new Color( 255, 255, 0 );
    static final Color        USER_FLAGGED_COLOR = new Color( 255, 0, 255 );
    static final Color        SELECTED_COLOR     = new Color( 255, 0, 0 );
    int                       _x;
    int                       _y;
    int                       _well_size;
    byte                      _status;

    public AbstractRenderer() {
    }

    abstract MsaRenderer getParentPlateRenderer();

    byte getStatus() {
        return _status;
    }

    int getWellSize() {
        return _well_size;
    }

    @Override
    public int getX() {
        return _x;
    }

    @Override
    public int getY() {
        return _y;
    }

    abstract boolean isSelected();

    @Override
    public abstract void paint( Graphics g );

    abstract void setIsSelected( boolean flag );

    void setStatus( final byte status ) {
        _status = status;
    }

    void setWellSize( final int well_size ) {
        _well_size = well_size;
    }

    void setX( final int x ) {
        _x = x;
    }

    void setY( final int y ) {
        _y = y;
    }
}
