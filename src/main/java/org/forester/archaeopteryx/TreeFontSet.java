// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// Copyright (C) 2003-2007 Ethalinda K.S. Cannon
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

package org.forester.archaeopteryx;

import java.awt.Font;
import java.awt.FontMetrics;

/*
 * Maintains the fonts for drawing a tree.
 */
public final class TreeFontSet {

    static final int            BOLD_AND_ITALIC           = Font.BOLD + Font.ITALIC;
    final static float          FONT_SIZE_CHANGE_STEP     = 1.0f;
    final static float          SMALL_FONTS_BASE          = 8;
    private final static String DEFAULT_FONT              = "Verdana";
    private Font                _base_font;
    private boolean             _decreased_size_by_system = false;
    private FontMetrics         _fm_large;
    // Handy holders for font metrics
    private FontMetrics         _fm_small;
    private Font                _large_font;
    private Font                _large_font_memory;
    private Font                _large_font_system;
    private final int           _max;
    private final int           _min;
    // the owner (needed to get font metrics)
    private final MainPanel     _owner;
    // The fonts
    private Font                _small_font;
    private Font                _small_font_memory;
    private Font                _small_font_system;
    private int                 _small_max_ascent         = 0;
    // hold font measurements
    private int                 _small_max_descent        = 0;

    TreeFontSet( final MainPanel owner ) {
        _owner = owner;
        _min = _owner.getConfiguration().getMinBaseFontSize();
        _max = _owner.getConfiguration().getMinBaseFontSize();
        setBaseFont( new Font( DEFAULT_FONT, Font.PLAIN, 10 ) );
    }

    public FontMetrics getFontMetricsLarge() {
        return _fm_large;
    }

    public FontMetrics getFontMetricsSmall() {
        return _fm_small;
    }

    public Font getSmallFont() {
        return _small_font;
    }

    public int getSmallMaxAscent() {
        return _small_max_ascent;
    }

    public int getSmallMaxDescent() {
        return _small_max_descent;
    }

    private Font getLargeFontSystem() {
        return _large_font_system;
    }

    private void intializeFonts() {
        final int small_size = getBaseFont().getSize() - 2;
        int italic = Font.ITALIC;
        if ( getBaseFont().getStyle() == Font.BOLD ) {
            italic = italic + Font.BOLD;
        }
        _small_font = new Font( getBaseFont().getFontName(), getBaseFont().getStyle(), small_size );
        _large_font = new Font( getBaseFont().getFontName(), getBaseFont().getStyle(), getBaseFont().getSize() );
        _small_font_system = new Font( getBaseFont().getFontName(), getBaseFont().getStyle(), small_size );
        _large_font_system = new Font( getBaseFont().getFontName(), getBaseFont().getStyle(), getBaseFont().getSize() );
        _small_font_memory = _small_font;
        _large_font_memory = _large_font;
        setupFontMetrics();
    }

    private void setDecreasedSizeBySystem( final boolean decreased_size_by_system ) {
        _decreased_size_by_system = decreased_size_by_system;
    }

    private void setupFontMetrics() {
        _fm_small = _owner.getFontMetrics( _small_font );
        _fm_large = _owner.getFontMetrics( _large_font );
        _small_max_descent = _fm_small.getMaxDescent();
        _small_max_ascent = _fm_small.getMaxAscent() + 1;
    }

    void decreaseFontSize( final int min, final boolean decreased_size_by_system ) {
        if ( decreased_size_by_system && !isDecreasedSizeBySystem() ) {
            _small_font_memory = _small_font;
            _large_font_memory = _large_font;
        }
        setDecreasedSizeBySystem( decreased_size_by_system );
        if ( _large_font.getSize() >= min ) {
            _small_font = _small_font.deriveFont( _small_font.getSize() - FONT_SIZE_CHANGE_STEP );
            _large_font = _large_font.deriveFont( _large_font.getSize() - FONT_SIZE_CHANGE_STEP );
            setupFontMetrics();
        }
    }

    Font getBaseFont() {
        return _base_font;
    }

    Font getLargeFont() {
        return _large_font;
    }

    Font getLargeFontMemory() {
        return _large_font_memory;
    }

    Font getSmallFontSystem() {
        return _small_font_system;
    }

    void increaseFontSize() {
        _small_font = _small_font.deriveFont( _small_font.getSize() + FONT_SIZE_CHANGE_STEP );
        _large_font = _large_font.deriveFont( _large_font.getSize() + FONT_SIZE_CHANGE_STEP );
        setupFontMetrics();
    }

    boolean isDecreasedSizeBySystem() {
        return _decreased_size_by_system;
    }

    void largeFonts() {
        setDecreasedSizeBySystem( false );
        _small_font = _small_font.deriveFont( 12f );
        _large_font = _large_font.deriveFont( 14f );
        setupFontMetrics();
    }

    void mediumFonts() {
        setDecreasedSizeBySystem( false );
        _small_font = _small_font.deriveFont( 8f );
        _large_font = _large_font.deriveFont( 10f );
        setupFontMetrics();
    }

    void reset() {
        _large_font_system = _large_font;
    }

    void setBaseFont( final Font base_font ) {
        _base_font = base_font;
        intializeFonts();
    }

    void smallFonts() {
        setDecreasedSizeBySystem( false );
        _small_font = _small_font.deriveFont( SMALL_FONTS_BASE - 2 );
        _large_font = _large_font.deriveFont( SMALL_FONTS_BASE );
        setupFontMetrics();
    }

    void superTinyFonts() {
        setDecreasedSizeBySystem( false );
        _small_font = _small_font.deriveFont( 2f );
        _large_font = _large_font.deriveFont( 4f );
        setupFontMetrics();
    }

    void tinyFonts() {
        setDecreasedSizeBySystem( false );
        _small_font = _small_font.deriveFont( 4f );
        _large_font = _large_font.deriveFont( 6f );
        setupFontMetrics();
    }
}
