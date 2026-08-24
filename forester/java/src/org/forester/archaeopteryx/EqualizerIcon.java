// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.archaeopteryx;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import javax.swing.Icon;
import javax.swing.UIManager;

/**
 * A small animated "equalizer" activity indicator for the menu bar: a few short bars that gently bounce in a
 * phased sine wave, in the FlatLaf theme accent color. Used on the "processes running" menu to draw the eye to
 * a background task without a (too-large, dated) spinning circle. The paint is a pure function of {@link
 * #advance() the phase}; a Swing {@code Timer} advances the phase and repaints the host component (see
 * {@code MainFrame}). Gentle by design -- a low frame rate + small amplitude.
 */
final class EqualizerIcon implements Icon {

    // used only if the current look-and-feel exposes no accent color (e.g. a non-FlatLaf/headless context)
    private static final Color  FALLBACK_ACCENT = new Color( 0x26, 0x75, 0xbf );
    // radians the whole pattern advances per animation tick (small -> a gentle bounce at a low frame rate)
    static final double         PHASE_STEP      = 0.42;
    // phase offset between adjacent bars, so they bounce out of step (the classic equalizer look)
    private static final double BAR_OFFSET      = 0.9;

    private final int _width;
    private final int _height;
    private final int _bars;
    private double    _phase = 0.0;

    EqualizerIcon( final int width, final int height, final int bars ) {
        _width = width;
        _height = height;
        _bars = Math.max( 1, bars );
    }

    /** Advance the animation by one step (called on each timer tick). */
    void advance() {
        _phase += PHASE_STEP;
    }

    /** Test hook: set the phase directly so two frames can be compared. */
    void setPhaseForTest( final double phase ) {
        _phase = phase;
    }

    @Override
    public int getIconWidth() {
        return _width;
    }

    @Override
    public int getIconHeight() {
        return _height;
    }

    @Override
    public void paintIcon( final Component c, final Graphics g, final int x, final int y ) {
        final Graphics2D g2 = (Graphics2D) g.create();
        try {
            g2.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
            Color accent = UIManager.getColor( "Component.accentColor" );
            if ( accent == null ) {
                accent = FALLBACK_ACCENT;
            }
            g2.setColor( accent );
            final int gap = 1;
            final int bar_w = Math.max( 1, ( _width - ( ( _bars - 1 ) * gap ) ) / _bars );
            final int min_h = Math.max( 2, _height / 4 );
            for ( int i = 0; i < _bars; i++ ) {
                final double s = 0.5 + ( 0.5 * Math.sin( _phase + ( i * BAR_OFFSET ) ) );
                final int bar_h = min_h + ( int ) Math.round( s * ( _height - min_h ) );
                final int bx = x + ( i * ( bar_w + gap ) );
                final int by = y + ( _height - bar_h );
                g2.fillRoundRect( bx, by, bar_w, bar_h, 2, 2 );
            }
        }
        finally {
            g2.dispose();
        }
    }
}
