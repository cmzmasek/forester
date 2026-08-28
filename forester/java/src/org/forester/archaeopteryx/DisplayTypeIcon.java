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
import java.awt.geom.Line2D;

import javax.swing.Icon;

/**
 * The glyphs on the control panel's phylogram/cladogram row -- the former "P" / "A" / "C" letter buttons.
 * <p>
 * Each glyph is a mini root-left tree PLUS the tip labels it produces, drawn as short ticks past the tips. The
 * labels are the point: what actually changes on screen between the three settings is not just the branches but
 * where the reader's eye finds the names, and a glyph that shows only branches leaves the most useful half of
 * the difference out. So the three are told apart by two independent cues:
 * <table>
 * <tr><th></th><th>branches</th><th>labels</th></tr>
 * <tr><td>{@code PHYLOGRAM}</td><td>ragged (lengths to scale)</td><td>ragged -- each label sits at its own tip</td></tr>
 * <tr><td>{@code ALIGNED_PHYLOGRAM}</td><td>ragged (lengths to scale)</td><td>lined up in one column</td></tr>
 * <tr><td>{@code CLADOGRAM}</td><td>flush (lengths ignored)</td><td>lined up in one column</td></tr>
 * </table>
 * Adjacent pairs therefore differ in exactly one cue, which is what makes the row readable at a glance: P vs A
 * changes the labels, A vs C changes the branches.
 * <p>
 * The glyph is WIDER than it is tall, because a tree plus a label column is a wide thing and the three buttons
 * are the widest on the control panel; a square icon would have squeezed it for no reason.
 * <p>
 * Vector-drawn (crisp at any {@code flatlaf.uiScale}) and painted in the host component's foreground, so it
 * follows the FlatLaf theme and greys out by itself on a disabled button.
 */
final class DisplayTypeIcon implements Icon {

    enum Kind {
        PHYLOGRAM, ALIGNED_PHYLOGRAM, CLADOGRAM
    }

    /** Width as a multiple of height -- room for the tree AND the label column past it. */
    private static final double ASPECT     = 1.6;

    // Normalized geometry: x as a fraction of the icon WIDTH, y as a fraction of its HEIGHT. Root at left.
    private static final double ROOT_X     = 0.02;
    private static final double SPINE_X    = 0.13;
    // The fork carrying the upper two tips. Every RAGGED end must sit to the RIGHT of the fork its branch leaves
    // (SUB_X for the upper two, SPINE_X for the lowest) -- a "shorter" tip drawn back past its own fork is not a
    // short branch, it is a backwards one. That is why the SHORTEST tip here is the bottom one: it hangs off the
    // spine, so it has the room to be dramatically short.
    private static final double SUB_X      = 0.30;
    private static final double[] TIP_Y    = { 0.13, 0.45, 0.85 };
    private static final double UPPER_Y    = 0.29;  // where the upper two tips' branch leaves the spine
    private static final double ROOT_Y     = 0.57;  // midway down the spine
    /**
     * Where the three branches END when lengths are drawn to scale. PHYLOGRAM and ALIGNED_PHYLOGRAM draw the
     * very same tree from these -- they differ ONLY in where the labels go -- so the raggedness is pushed well
     * past what a real tree would show: at 17 px a subtle length difference is no difference at all.
     */
    private static final double[] RAGGED   = { 0.62, 0.42, 0.26 };
    /** ... and where they all end when lengths are ignored (a cladogram): flush with the longest of them. */
    private static final double FLUSH      = 0.62;
    private static final double LABEL_GAP  = 0.11;  // clear air between a tip and its label
    private static final double LABEL_LEN  = 0.22;
    /** The common label column: the longest branch plus the gap, so no label can collide with a branch. */
    private static final double LABEL_COL  = FLUSH + LABEL_GAP;

    private final Kind _kind;
    private final int  _height;
    private final int  _width;

    /** @param height the icon's height in pixels; its width follows from {@link #ASPECT}. */
    DisplayTypeIcon( final Kind kind, final int height ) {
        _kind = kind;
        _height = Math.max( 6, height );
        _width = ( int ) Math.round( _height * ASPECT );
    }

    Kind getKind() {
        return _kind;
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
            g2.setColor( ( c != null ) ? c.getForeground() : Color.BLACK );
            g2.setStroke( new BasicStroke( Math.max( 1.0f, ( float ) ( _height * 0.085 ) ), BasicStroke.CAP_BUTT,
                                           BasicStroke.JOIN_MITER ) );
            final boolean flush = ( _kind == Kind.CLADOGRAM );
            final double[] ends = { flush ? FLUSH : RAGGED[ 0 ], flush ? FLUSH : RAGGED[ 1 ],
                                    flush ? FLUSH : RAGGED[ 2 ] };
            // the tree
            line( g2, x, y, ROOT_X, ROOT_Y, SPINE_X, ROOT_Y );      // root stub
            line( g2, x, y, SPINE_X, UPPER_Y, SPINE_X, TIP_Y[ 2 ] ); // spine
            line( g2, x, y, SPINE_X, UPPER_Y, SUB_X, UPPER_Y );      // branch to the upper clade
            line( g2, x, y, SUB_X, TIP_Y[ 0 ], SUB_X, TIP_Y[ 1 ] );  // that clade's own fork
            line( g2, x, y, SUB_X, TIP_Y[ 0 ], ends[ 0 ], TIP_Y[ 0 ] );
            line( g2, x, y, SUB_X, TIP_Y[ 1 ], ends[ 1 ], TIP_Y[ 1 ] );
            line( g2, x, y, SPINE_X, TIP_Y[ 2 ], ends[ 2 ], TIP_Y[ 2 ] );
            // the tip labels: at each tip in a plain phylogram, in one column otherwise
            final boolean aligned = ( _kind != Kind.PHYLOGRAM );
            for ( int i = 0; i < TIP_Y.length; i++ ) {
                final double label_x = aligned ? LABEL_COL : ( ends[ i ] + LABEL_GAP );
                line( g2, x, y, label_x, TIP_Y[ i ], label_x + LABEL_LEN, TIP_Y[ i ] );
            }
        }
        finally {
            g2.dispose();
        }
    }

    /** Test hook: the pixel rows the three tip branches (and their labels) sit on, at a given icon height. */
    static int[] tipRowsForTest( final int height ) {
        final int h = Math.max( 6, height );
        final int[] rows = new int[ TIP_Y.length ];
        for ( int i = 0; i < TIP_Y.length; i++ ) {
            rows[ i ] = ( int ) Math.round( TIP_Y[ i ] * h );
        }
        return rows;
    }

    private void line( final Graphics2D g2, final int x, final int y, final double x1, final double y1,
                       final double x2, final double y2 ) {
        g2.draw( new Line2D.Double( x + ( x1 * _width ), y + ( y1 * _height ), x + ( x2 * _width ),
                                    y + ( y2 * _height ) ) );
    }
}
