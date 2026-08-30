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

import java.awt.BasicStroke;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;

import javax.swing.Icon;

/**
 * The glyph on the control panel's single light/dark theme toggle: either a SUN (a small disc with eight rays)
 * or a crescent MOON. The button offers the theme it will switch TO -- a moon while the light theme is active,
 * a sun while the dark one is -- so exactly one icon is showing at any time (this replaced the former pair of
 * "Light"/"Dark" radio buttons).
 * <p>
 * Drawn with Java2D rather than a bitmap so it stays crisp at any {@code flatlaf.uiScale} / GUI font size, and
 * painted in the host component's foreground color, so it follows the FlatLaf theme and greys out by itself on
 * a disabled button.
 */
final class ThemeToggleIcon implements Icon {

    /** Number of rays around the sun's disc. */
    private static final int    SUN_RAYS        = 8;
    /** Sun disc radius, as a fraction of the icon size. */
    private static final double SUN_DISC_R      = 0.20;
    /** Inner / outer end of a sun ray, as a fraction of the icon size. */
    private static final double SUN_RAY_INNER_R = 0.29;
    private static final double SUN_RAY_OUTER_R = 0.43;
    /** Moon: the full disc, and the offset disc that is subtracted from it to carve the crescent. */
    private static final double MOON_R          = 0.44;
    private static final double MOON_BITE_R     = 0.40;
    private static final double MOON_BITE_DIST  = 0.33;
    /** Direction (up and to the right) of the subtracted disc, so the crescent's horns point up-left/down-left. */
    private static final double MOON_BITE_ANGLE = Math.toRadians( -35 );

    private final boolean _sun;
    private final int     _size;

    /**
     * @param sun true for the sun glyph (shown while the dark theme is active), false for the crescent moon
     * @param size the icon's width and height in pixels
     */
    ThemeToggleIcon( final boolean sun, final int size ) {
        _sun = sun;
        _size = Math.max( 4, size );
    }

    /** True if this icon paints the sun (i.e. the button offers the LIGHT theme). */
    boolean isSun() {
        return _sun;
    }

    /** True if this icon paints the crescent moon (i.e. the button offers the DARK theme). */
    boolean isMoon() {
        return !_sun;
    }

    @Override
    public int getIconWidth() {
        return _size;
    }

    @Override
    public int getIconHeight() {
        return _size;
    }

    @Override
    public void paintIcon( final Component c, final Graphics g, final int x, final int y ) {
        final Graphics2D g2 = (Graphics2D) g.create();
        try {
            g2.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
            g2.setColor( ControlButtonIcon.inkFor( c ) );
            final double cx = x + ( _size / 2.0 );
            final double cy = y + ( _size / 2.0 );
            if ( _sun ) {
                paintSun( g2, cx, cy );
            }
            else {
                paintMoon( g2, cx, cy );
            }
        }
        finally {
            g2.dispose();
        }
    }

    private void paintSun( final Graphics2D g2, final double cx, final double cy ) {
        final double disc_r = SUN_DISC_R * _size;
        g2.fill( new Ellipse2D.Double( cx - disc_r, cy - disc_r, 2 * disc_r, 2 * disc_r ) );
        g2.setStroke( new BasicStroke( Math.max( 1.1f, ( float ) ( _size * 0.09 ) ), BasicStroke.CAP_ROUND,
                                       BasicStroke.JOIN_ROUND ) );
        final double inner = SUN_RAY_INNER_R * _size;
        final double outer = SUN_RAY_OUTER_R * _size;
        for ( int i = 0; i < SUN_RAYS; i++ ) {
            final double a = ( i * 2 * Math.PI ) / SUN_RAYS;
            final double dx = Math.cos( a );
            final double dy = Math.sin( a );
            g2.draw( new java.awt.geom.Line2D.Double( cx + ( dx * inner ), cy + ( dy * inner ),
                                                      cx + ( dx * outer ), cy + ( dy * outer ) ) );
        }
    }

    private void paintMoon( final Graphics2D g2, final double cx, final double cy ) {
        final double r = MOON_R * _size;
        final double bite_r = MOON_BITE_R * _size;
        final double bx = cx + ( Math.cos( MOON_BITE_ANGLE ) * MOON_BITE_DIST * _size );
        final double by = cy + ( Math.sin( MOON_BITE_ANGLE ) * MOON_BITE_DIST * _size );
        final Area crescent = new Area( new Ellipse2D.Double( cx - r, cy - r, 2 * r, 2 * r ) );
        crescent.subtract( new Area( new Ellipse2D.Double( bx - bite_r, by - bite_r, 2 * bite_r, 2 * bite_r ) ) );
        g2.fill( crescent );
    }
}
