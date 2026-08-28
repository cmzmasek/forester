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
import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;

import javax.swing.Icon;

/**
 * The glyphs on the control panel's layout row: one small tree silhouette per PRIMARY DISPLAY TYPE. There are
 * exactly five, matching the project's "all five display types at parity" rule -- the three rectangular
 * orientations (root at left / top / bottom) plus the two radial layouts (circular, unrooted). The row is a
 * single exclusive group, so one click reaches any layout from any other.
 * <p>
 * The three rectangular glyphs are the SAME mini tree rotated, which is exactly what the orientations are, so
 * they read as a family. The circular glyph is a two-level circular dendrogram (centre hub, three arcs, six
 * radial tips) -- deliberately arc-bearing so it cannot be confused with the sun on the theme toggle. The
 * unrooted glyph is a free-form star with irregular angles and no hub, the way an unrooted tree actually draws.
 * <p>
 * Vector-drawn (crisp at any {@code flatlaf.uiScale}) and painted in the host component's foreground, so it
 * follows the FlatLaf theme and greys out by itself on a disabled button.
 */
final class LayoutIcon implements Icon {

    /** The five primary display types, in the order they appear on the control panel. */
    enum Kind {
        ROOT_LEFT, ROOT_TOP, ROOT_BOTTOM, CIRCULAR, UNROOTED
    }

    // ---- the rectangular mini tree, in a normalized [0,1] x [0,1] box, drawn root-at-LEFT --------------
    // A root stub, a spine that forks into an upper and a lower branch, and the upper branch forking again:
    // three tips, which is the fewest that still reads as a tree rather than a bracket.
    private static final double[][] RECT_SEGMENTS = {
            { 0.06, 0.50, 0.24, 0.50 }, // root stub
            { 0.24, 0.24, 0.24, 0.76 }, // spine
            { 0.24, 0.24, 0.52, 0.24 }, // upper branch
            { 0.24, 0.76, 0.94, 0.76 }, // lower branch -> tip 3
            { 0.52, 0.10, 0.52, 0.38 }, // sub-spine
            { 0.52, 0.10, 0.94, 0.10 }, // tip 1
            { 0.52, 0.38, 0.94, 0.38 }, // tip 2
    };

    // ---- the circular mini tree ----------------------------------------------------------------------
    // The unmistakable "circular" cue at 16 px is a RING, so the tip circle is drawn as one big arc and the
    // branches are spokes hanging inside it. The arc is deliberately left OPEN (a wedge missing on the right,
    // where a real circular tree's root sits), which is also what keeps it from reading as the theme toggle's
    // sun -- the sun's rays are detached from its disc and run the whole way round.
    private static final double   CIRC_HUB_R      = 0.075; // filled centre dot (the root)
    private static final double   CIRC_SPOKE_IN_R = 0.13;  // where a spoke leaves the hub
    private static final double   CIRC_RIM_R      = 0.43;  // the tip circle
    private static final double   CIRC_RIM_START  = 32;    // arc start, in screen degrees
    private static final double   CIRC_RIM_EXTENT = 296;   // arc sweep -- the remainder is the root wedge
    private static final double[] CIRC_SPOKES     = { 75, 150, 225, 300 };

    // ---- the unrooted mini tree ----------------------------------------------------------------------
    // Irregular angles and lengths on purpose: an unrooted tree has no centre and no symmetry to it.
    private static final double[] UNROOTED_ANGLES = { -105, 20, 145 };
    private static final double[] UNROOTED_LENGTH = { 0.26, 0.21, 0.24 };
    private static final double   UNROOTED_FORK   = 0.20; // length of the two tips off each edge
    private static final double   UNROOTED_SPREAD = 38;   // half-angle between the two tips, in degrees

    private final Kind _kind;
    private final int  _size;

    LayoutIcon( final Kind kind, final int size ) {
        _kind = kind;
        _size = Math.max( 6, size );
    }

    Kind getKind() {
        return _kind;
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
            g2.setColor( ( c != null ) ? c.getForeground() : Color.BLACK );
            g2.setStroke( new BasicStroke( strokeWidth(), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND ) );
            switch ( _kind ) {
                case ROOT_LEFT:
                    paintRectangular( g2, x, y, 0 );
                    break;
                case ROOT_TOP:
                    paintRectangular( g2, x, y, 90 ); // right -> down, so the root ends up at the top
                    break;
                case ROOT_BOTTOM:
                    paintRectangular( g2, x, y, -90 ); // right -> up, so the root ends up at the bottom
                    break;
                case CIRCULAR:
                    paintCircular( g2, x, y );
                    break;
                case UNROOTED:
                    paintUnrooted( g2, x, y );
                    break;
            }
        }
        finally {
            g2.dispose();
        }
    }

    private float strokeWidth() {
        return Math.max( 1.0f, ( float ) ( _size * 0.075 ) );
    }

    /** The root-left mini tree, optionally rotated about the icon's centre to give the other two orientations. */
    private void paintRectangular( final Graphics2D g2, final int x, final int y, final double degrees ) {
        if ( degrees != 0 ) {
            g2.rotate( Math.toRadians( degrees ), x + ( _size / 2.0 ), y + ( _size / 2.0 ) );
        }
        for( final double[] s : RECT_SEGMENTS ) {
            g2.draw( new Line2D.Double( x + ( s[ 0 ] * _size ), y + ( s[ 1 ] * _size ), x + ( s[ 2 ] * _size ),
                                        y + ( s[ 3 ] * _size ) ) );
        }
    }

    private void paintCircular( final Graphics2D g2, final int x, final int y ) {
        final double cx = x + ( _size / 2.0 );
        final double cy = y + ( _size / 2.0 );
        final double hub = CIRC_HUB_R * _size;
        g2.fill( new Ellipse2D.Double( cx - hub, cy - hub, 2 * hub, 2 * hub ) );
        final double r = CIRC_RIM_R * _size;
        // Arc2D measures degrees counter-clockwise from 3 o'clock, while everything else here uses the screen
        // convention (y down, so clockwise) -- hence the negated start and extent.
        g2.draw( new Arc2D.Double( cx - r, cy - r, 2 * r, 2 * r, -CIRC_RIM_START, -CIRC_RIM_EXTENT, Arc2D.OPEN ) );
        for( final double a : CIRC_SPOKES ) {
            radial( g2, cx, cy, a, CIRC_SPOKE_IN_R, CIRC_RIM_R );
        }
    }

    private void paintUnrooted( final Graphics2D g2, final int x, final int y ) {
        final double cx = x + ( _size / 2.0 );
        final double cy = y + ( _size / 2.0 );
        for( int i = 0; i < UNROOTED_ANGLES.length; i++ ) {
            final double a = Math.toRadians( UNROOTED_ANGLES[ i ] );
            final double ex = cx + ( Math.cos( a ) * UNROOTED_LENGTH[ i ] * _size );
            final double ey = cy + ( Math.sin( a ) * UNROOTED_LENGTH[ i ] * _size );
            g2.draw( new Line2D.Double( cx, cy, ex, ey ) );
            for( final int sign : new int[] { -1, 1 } ) {
                final double fa = a + ( sign * Math.toRadians( UNROOTED_SPREAD ) );
                g2.draw( new Line2D.Double( ex, ey, ex + ( Math.cos( fa ) * UNROOTED_FORK * _size ),
                                            ey + ( Math.sin( fa ) * UNROOTED_FORK * _size ) ) );
            }
        }
    }

    /** A line along the given angle (degrees, screen convention) between two radii given as size fractions. */
    private void radial( final Graphics2D g2, final double cx, final double cy, final double degrees,
                         final double from_r, final double to_r ) {
        final double a = Math.toRadians( degrees );
        g2.draw( new Line2D.Double( cx + ( Math.cos( a ) * from_r * _size ), cy + ( Math.sin( a ) * from_r * _size ),
                                    cx + ( Math.cos( a ) * to_r * _size ), cy + ( Math.sin( a ) * to_r * _size ) ) );
    }
}
