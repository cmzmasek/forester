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
import java.awt.Graphics2D;
import java.awt.RenderingHints;

public class ResidueRenderer extends AbstractRenderer {

    static final Color        EMPTY_COLOR          = new Color( 250, 0, 250 );
    static final Color        POSITIVE_COLOR       = new Color( 250, 0, 250 );
    static final Color        NEGATIVE_COLOR       = new Color( 250, 0, 250 );
    static final Color        NULL_COLOR           = new Color( 50, 50, 50 );
    static final int          DISTANCE_OVAL_BORDER = 1;
    static final int          SIZE_LIMIT           = 7;
    /**
     *
     */
    private static final long serialVersionUID     = -2331160296913478874L;
    private final char        _value;
    private Color             _wellColor;
    private boolean           _isMarked;
    private boolean           _isSelected;
    private final MsaRenderer _parentPlateRenderer;

    public ResidueRenderer( final char value, final MsaRenderer parentPlateRenderer ) {
        _value = value;
        _parentPlateRenderer = parentPlateRenderer;
        setIsSelected( false );
        setIsMarked( false );
        setStatus( ( byte ) 0 );
    }

    private double calcFactor( final double min, final double max ) {
        return ( max - min ) / 255D;
    }

    private Color calcWellColor( double value,
                                 final double min,
                                 final double max,
                                 final Color minColor,
                                 final Color maxColor ) {
        if ( value < min ) {
            value = min;
        }
        if ( value > max ) {
            value = max;
        }
        final double x = ( 255D * ( value - min ) ) / ( max - min );
        final int red = ( int ) ( minColor.getRed() + ( x * calcFactor( minColor.getRed(), maxColor.getRed() ) ) );
        final int green = ( int ) ( minColor.getGreen() + ( x * calcFactor( minColor.getGreen(), maxColor.getGreen() ) ) );
        final int blue = ( int ) ( minColor.getBlue() + ( x * calcFactor( minColor.getBlue(), maxColor.getBlue() ) ) );
        return new Color( red, green, blue );
    }

    private Color calcWellColor( double value,
                                 final double min,
                                 final double max,
                                 final double mean,
                                 final Color minColor,
                                 final Color maxColor,
                                 final Color meanColor ) {
        //  if ( isEmpty() ) {
        //     return ResidueRenderer.NULL_COLOR;
        // }
        if ( meanColor == null ) {
            return calcWellColor( value, min, max, minColor, maxColor );
        }
        if ( value < min ) {
            value = min;
        }
        if ( value > max ) {
            value = max;
        }
        if ( value < mean ) {
            final double x = ( 255D * ( value - min ) ) / ( mean - min );
            final int red = ( int ) ( minColor.getRed() + ( x * calcFactor( minColor.getRed(), meanColor.getRed() ) ) );
            final int green = ( int ) ( minColor.getGreen() + ( x * calcFactor( minColor.getGreen(),
                                                                                meanColor.getGreen() ) ) );
            final int blue = ( int ) ( minColor.getBlue() + ( x * calcFactor( minColor.getBlue(), meanColor.getBlue() ) ) );
            return new Color( red, green, blue );
        }
        if ( value > mean ) {
            final double x = ( 255D * ( value - mean ) ) / ( max - mean );
            final int red = ( int ) ( meanColor.getRed() + ( x * calcFactor( meanColor.getRed(), maxColor.getRed() ) ) );
            final int green = ( int ) ( meanColor.getGreen() + ( x * calcFactor( meanColor.getGreen(),
                                                                                 maxColor.getGreen() ) ) );
            final int blue = ( int ) ( meanColor.getBlue() + ( x * calcFactor( meanColor.getBlue(), maxColor.getBlue() ) ) );
            return new Color( red, green, blue );
        }
        else {
            return meanColor;
        }
    }

    public double getDataValue() {
        return _value;
    }

    @Override
    public MsaRenderer getParentPlateRenderer() {
        return _parentPlateRenderer;
    }

    public Color getWellColor() {
        return _wellColor;
    }

    public boolean isMarked() {
        return _isMarked;
    }

    @Override
    public boolean isSelected() {
        return _isSelected;
    }

    @Override
    public void paint( final Graphics g ) {
        final int width = getWellSize() - 1;
        final int width_ = width - 1;
        final int width__ = ( width_ - 1 ) + 1;
        final int width__s = width__ - 2;
        final int x_ = getX() + 1;
        final int y_ = getY() + 1;
        //     final PlateDisplayPanel hmp = getParentPlateRenderer()
        //             .getPlateDisplayPanel();
        //     boolean draw_circle = hmp.isDrawCircle()
        //            || ( !hmp.isDrawCircle() && !hmp.isDrawSquare() && ( width > 7 ) );
        //  final boolean show_user_flags = _parentPlateRenderer
        //          .getPlateDisplayPanel().showUserFlagsCheckBox.isSelected();
        //  final boolean show_outlier_flags = _parentPlateRenderer
        //          .getPlateDisplayPanel().showOutliersCheckBox.isSelected();
        // final boolean show_hit_picks = _parentPlateRenderer
        // .getPlateDisplayPanel().showHitPicksCheckBox.isSelected();
        // final boolean show_confirmed_hits = _parentPlateRenderer
        // .getPlateDisplayPanel().showConfirmedHitsCheckBox.isSelected();
        final Graphics2D g2 = ( Graphics2D ) g;
        g2.setRenderingHint( RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED );
        if ( isMarked() ) {
            g2.setColor( AbstractRenderer.MARKED_COLOR );
        }
        // else if ( !draw_circle && isSelected() ) {
        //     g2.setColor( AbstractRenderer.SELECTED_COLOR );
        // }
        else {
            g2.setColor( AbstractRenderer.DEFAULT_COLOR );
        }
        g2.drawRect( getX(), getY(), width, width );
        //        if ( draw_circle ) {
        //            if ( isSelected() && isMarked() ) {
        //                g2.setColor( AbstractRenderer.MARKED_COLOR );
        //            }
        //            else if ( isSelected() ) {
        //                g2.setColor( AbstractRenderer.SELECTED_COLOR );
        //            }
        //            else {
        //                g2.setColor( AbstractRenderer.DEFAULT_COLOR );
        //            }
        //            g2.fillRect( x_, y_, width_, width_ );
        //        }
        g2.setColor( getWellColor() );
        //        if ( draw_circle ) {
        //            g2.setRenderingHint( RenderingHints.KEY_ANTIALIASING,
        //                    RenderingHints.VALUE_ANTIALIAS_ON );
        //            if ( isSelected() && ( width > 6 ) ) {
        //                g2.fillOval( getX() + 2, getY() + 2, width__s, width__s );
        //            }
        //            else if ( width < 5 ) {
        //                g2.fillOval( ( getX() + 1 ) - 1, ( getY() + 1 ) - 1,
        //                        width__ + 2, width__ + 2 );
        //            }
        //            else {
        //                g2.fillOval( getX() + 1, getY() + 1, width__, width__ );
        //            }
        //        }
        if ( isMarked() || isSelected() ) {
            g2.fillRect( getX() + 1, getY() + 1, width_, width_ );
        }
        else {
            g2.fillRect( ( getX() + 1 ) - 1, ( getY() + 1 ) - 1, width_ + 2, width_ + 2 );
        }
    }

    public void resetWellColor( final double min, final double max, final Color minColor, final Color maxColor ) {
        setWellColor( calcWellColor( getDataValue(), min, max, minColor, maxColor ) );
    }

    public void resetWellColor( final double min,
                                final double max,
                                final double mean,
                                final Color minColor,
                                final Color maxColor,
                                final Color meanColor ) {
        setWellColor( calcWellColor( getDataValue(), min, max, mean, minColor, maxColor, meanColor ) );
    }

    public void setIsMarked( final boolean isMarked ) {
        _isMarked = isMarked;
    }

    @Override
    public void setIsSelected( final boolean isSelected ) {
        _isSelected = isSelected;
    }

    private void setWellColor( final Color wellColor ) {
        _wellColor = wellColor;
    }
}
