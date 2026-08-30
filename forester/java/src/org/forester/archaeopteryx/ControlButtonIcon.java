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
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;

import javax.swing.Icon;

/**
 * Glyphs for the control panel's action buttons that used to be bare letters (or, for the two radial rotate
 * buttons, a Unicode circular arrow with a "CW"/"CCW" text fallback for fonts that lacked the character).
 * <p>
 * Drawing them removes the font dependency entirely -- there is nothing left to fall back from -- and lets them
 * tint with the FlatLaf theme and scale with the GUI font like every other icon on the panel.
 * <p>
 * Two deliberate choices, both about not colliding with glyphs the panel already shows:
 * <ul>
 * <li>the rotate arcs carry NO centre dot: a hub would make them near-twins of the circular <em>layout</em>
 * button one row above;</li>
 * <li>"return to the whole tree" is an arrow to a BAR and "up one level" the same arrow without it, so the two
 * read as one step versus all the way. Drawing the whole tree instead would have reproduced the layout row's
 * root-left glyph exactly, and a house would not have paired with anything.</li>
 * </ul>
 * Both arrows point LEFT, toward the root: in the classic root-left tree that is literally the direction of the
 * parent clade, so the glyph agrees with the tree on screen instead of borrowing the file-manager "up a level"
 * metaphor. The pair is then the familiar "back one" / "back to the start" idiom.
 */
final class ControlButtonIcon implements Icon {

    enum Kind {
        /** Rotate the radial layout clockwise (the old "↻" / "CW"). */
        ROTATE_CW,
        /** Rotate the radial layout counter-clockwise (the old "↺" / "CCW"). */
        ROTATE_CCW,
        /** Back to the complete tree from a subtree view (the old "R"). */
        WHOLE_TREE,
        /** Up to the immediate super-tree (the old "R1"). */
        UP_ONE_LEVEL,
        /** Uncollapse every collapsed clade (the old "U"). */
        UNCOLLAPSE_ALL,
        /** Fit the tree to the window WIDTH, keeping the vertical zoom (the old "W"). */
        FIT_WIDTH,
        /** Fit the tree to the window HEIGHT, keeping the horizontal zoom (the old "H"). */
        FIT_HEIGHT,
        /** Fit and centre the whole tree (the old "F"). */
        FIT_ALL,
        /** Spread the label rows apart vertically so they stop overlapping (the old "E", root-left layout). */
        EXPAND_VERTICAL,
        /** The same, horizontally (the root-top/bottom orientations). */
        EXPAND_HORIZONTAL,
        /** Switch the radial layouts' node labels to RIDE THE SPOKE (shown while they are horizontal). */
        LABELS_RADIAL,
        /** Switch them back to horizontal (shown while they ride the spoke). */
        LABELS_HORIZONTAL
    }

    private final Kind _kind;
    private final int  _size;

    ControlButtonIcon( final Kind kind, final int size ) {
        _kind = kind;
        _size = Math.max( 8, size );
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
            g2.setColor( inkFor( c ) );
            switch ( _kind ) {
                case ROTATE_CW:
                    paintRotate( g2, x, y, true );
                    break;
                case ROTATE_CCW:
                    paintRotate( g2, x, y, false );
                    break;
                case WHOLE_TREE:
                    paintBackArrow( g2, x, y, true );
                    break;
                case UP_ONE_LEVEL:
                    paintBackArrow( g2, x, y, false );
                    break;
                case FIT_WIDTH:
                    paintFit( g2, x, y, true, false );
                    break;
                case FIT_HEIGHT:
                    paintFit( g2, x, y, false, false );
                    break;
                case FIT_ALL:
                    paintFit( g2, x, y, false, true );
                    break;
                case EXPAND_VERTICAL:
                    paintExpand( g2, x, y, true );
                    break;
                case EXPAND_HORIZONTAL:
                    paintExpand( g2, x, y, false );
                    break;
                case LABELS_RADIAL:
                    paintLabelDirection( g2, x, y, true );
                    break;
                case LABELS_HORIZONTAL:
                    paintLabelDirection( g2, x, y, false );
                    break;
                case UNCOLLAPSE_ALL:
                    paintUncollapse( g2, x, y );
                    break;
                default:
                    // every Kind gets its own case, so a new one fails loudly here instead of silently
                    // painting whichever glyph happened to sit in the default arm
                    throw new IllegalStateException( "unhandled control-button glyph: " + _kind );
            }
        }
        finally {
            g2.dispose();
        }
    }

    /**
     * The colour a glyph paints in on {@code c}: its foreground -- FADED when the component is DISABLED.
     * Swing never fades it for us: FlatLaf's {@code getDisabledIcon()} synthesizes a greyed image only for an
     * {@code ImageIcon}, and for a custom {@link Icon} it returns null and paints the normal icon full-strength.
     * So without this, a disabled glyph button (R/R1 at startup, say) looks exactly like an enabled one -- the
     * affordance the old greyed text labels gave for free. Alpha rather than a fixed grey, so it mutes correctly
     * against both the light and the dark theme.
     */
    static Color inkFor( final Component c ) {
        final Color fg = ( c != null ) ? c.getForeground() : Color.BLACK;
        if ( ( c == null ) || c.isEnabled() ) {
            return fg;
        }
        // A SOLID blend toward the button's background, not the foreground at reduced alpha: translucent ink
        // COMPOUNDS wherever two strokes of the same glyph overlap (a shaft under its arrowhead, a bar crossing
        // a line), leaving dark spots in a supposedly faded icon. A solid colour paints the same however often
        // it is overdrawn -- and blending toward the actual background mutes correctly on both themes.
        final Color bg = ( c.getBackground() != null ) ? c.getBackground() : Color.GRAY;
        return new Color( blend( fg.getRed(), bg.getRed() ), blend( fg.getGreen(), bg.getGreen() ),
                          blend( fg.getBlue(), bg.getBlue() ) );
    }

    /** How much of the FOREGROUND survives in a disabled glyph's ink. */
    private static final double DISABLED_KEEP = 0.38;

    private static int blend( final int fg, final int bg ) {
        return ( int ) Math.round( ( fg * DISABLED_KEEP ) + ( bg * ( 1.0 - DISABLED_KEEP ) ) );
    }

    /**
     * The "fit" family (the old F / W / H letters): a rounded window FRAME with arrows pushing outward against
     * its edges. The frame's shape carries the axis -- LANDSCAPE with two horizontal arrows for fit-width,
     * PORTRAIT with two vertical arrows for fit-height, SQUARE with a two-headed diagonal (the established
     * fit-to-screen idiom, and calmer at 15 px than four arrows) for fit-everything.
     */
    private void paintFit( final Graphics2D g2, final int x, final int y, final boolean horizontal,
                           final boolean both ) {
        g2.setStroke( new BasicStroke( Math.max( 1.0f, _size * 0.065f ), BasicStroke.CAP_BUTT,
                                       BasicStroke.JOIN_ROUND ) );
        final double cx = x + ( _size / 2.0 );
        final double cy = y + ( _size / 2.0 );
        final double fw = _size * ( both ? 0.92 : ( horizontal ? 0.96 : 0.66 ) );
        final double fh = _size * ( both ? 0.92 : ( horizontal ? 0.66 : 0.96 ) );
        final double corner = _size * 0.18;
        g2.draw( new java.awt.geom.RoundRectangle2D.Double( cx - ( fw / 2 ), cy - ( fh / 2 ), fw, fh, corner,
                                                            corner ) );
        final double inset = _size * 0.10;
        final double gap = _size * 0.05;
        final double head = _size * 0.20;
        if ( both ) {
            // one two-headed diagonal, corner to corner
            fitArrow( g2, cx + gap, cy + gap, ( cx + ( fw / 2 ) ) - ( inset * 1.2 ),
                      ( cy + ( fh / 2 ) ) - ( inset * 1.2 ), head );
            fitArrow( g2, cx - gap, cy - gap, ( cx - ( fw / 2 ) ) + ( inset * 1.2 ),
                      ( cy - ( fh / 2 ) ) + ( inset * 1.2 ), head );
        }
        else if ( horizontal ) {
            fitArrow( g2, cx - gap, cy, ( cx - ( fw / 2 ) ) + inset, cy, head );
            fitArrow( g2, cx + gap, cy, ( cx + ( fw / 2 ) ) - inset, cy, head );
        }
        else {
            fitArrow( g2, cx, cy - gap, cx, ( cy - ( fh / 2 ) ) + inset, head );
            fitArrow( g2, cx, cy + gap, cx, ( cy + ( fh / 2 ) ) - inset, head );
        }
    }

    /**
     * "Expand to fit labels" (the old "E"): the LABEL ROWS as three short parallel lines, with an arrow beyond
     * each outer row pushing the stack apart -- the increase-line-spacing idiom, which is exactly what the
     * button does. Deliberately FRAMELESS: the fit family owns the window frame, and this is not a fit -- it
     * grows the layout so nothing has to be hidden. Vertical form for root-left (rows are horizontal, spread
     * vertically), rotated for the root-top/bottom orientations.
     */
    private void paintExpand( final Graphics2D g2, final int x, final int y, final boolean vertical ) {
        g2.setStroke( new BasicStroke( Math.max( 1.0f, _size * 0.075f ), BasicStroke.CAP_BUTT,
                                       BasicStroke.JOIN_ROUND ) );
        final double cx = x + ( _size / 2.0 );
        final double cy = y + ( _size / 2.0 );
        final double line_half = _size * 0.28;  // half-length of a label row's line
        final double row_off = _size * 0.16;    // the outer rows' offset from the centre row
        final double head = _size * 0.19;
        final double tip_at = _size * 0.485;    // the arrows' tips, just inside the icon edge
        final double tail_at = row_off + ( _size * 0.06 ); // arrows set off just past the outer rows
        if ( vertical ) {
            for ( final double off : new double[] { -row_off, 0, row_off } ) {
                g2.draw( new Line2D.Double( cx - line_half, cy + off, cx + line_half, cy + off ) );
            }
            fitArrow( g2, cx, cy - tail_at, cx, cy - tip_at, head );
            fitArrow( g2, cx, cy + tail_at, cx, cy + tip_at, head );
        }
        else {
            for ( final double off : new double[] { -row_off, 0, row_off } ) {
                g2.draw( new Line2D.Double( cx + off, cy - line_half, cx + off, cy + line_half ) );
            }
            fitArrow( g2, cx - tail_at, cy, cx - tip_at, cy, head );
            fitArrow( g2, cx + tail_at, cy, cx + tip_at, cy, head );
        }
    }

    /**
     * The radial label-direction flip (the old "L"): the tree's centre as a dot, one spoke running out to the
     * upper right, and a thick short BAR -- the label -- at the spoke's end: riding ALONG the spoke for the
     * radial-labels state, lying flat for the horizontal one. Like the theme toggle, the button shows the state
     * it will switch TO, so exactly one of the two forms is ever visible.
     */
    private void paintLabelDirection( final Graphics2D g2, final int x, final int y, final boolean radial ) {
        final double cx = x + ( _size * 0.22 );
        final double cy = y + ( _size * 0.78 );
        final double dot = _size * 0.10;
        g2.fill( new java.awt.geom.Ellipse2D.Double( cx - dot, cy - dot, 2 * dot, 2 * dot ) );
        g2.setStroke( new BasicStroke( Math.max( 1.0f, _size * 0.075f ), BasicStroke.CAP_BUTT,
                                       BasicStroke.JOIN_ROUND ) );
        final double a = Math.toRadians( -45 ); // the spoke, up and to the right
        final double spoke_r = _size * 0.52;
        final double sx = cx + ( Math.cos( a ) * spoke_r );
        final double sy = cy + ( Math.sin( a ) * spoke_r );
        g2.draw( new Line2D.Double( cx + ( Math.cos( a ) * dot * 1.6 ), cy + ( Math.sin( a ) * dot * 1.6 ), sx,
                                    sy ) );
        // the label: a thick bar just past the spoke's end -- along it, or lying flat
        g2.setStroke( new BasicStroke( Math.max( 2.0f, _size * 0.17f ), BasicStroke.CAP_BUTT,
                                       BasicStroke.JOIN_ROUND ) );
        final double bar = _size * 0.34;
        if ( radial ) {
            g2.draw( new Line2D.Double( sx + ( Math.cos( a ) * _size * 0.06 ),
                                        sy + ( Math.sin( a ) * _size * 0.06 ),
                                        sx + ( Math.cos( a ) * ( ( _size * 0.06 ) + bar ) ),
                                        sy + ( Math.sin( a ) * ( ( _size * 0.06 ) + bar ) ) ) );
        }
        else {
            g2.draw( new Line2D.Double( sx + ( _size * 0.04 ), sy, sx + ( _size * 0.04 ) + bar, sy ) );
        }
    }

    /** A short arrow from (x1,y1) whose head TIP lands exactly on (x2,y2). */
    private static void fitArrow( final Graphics2D g2, final double x1, final double y1, final double x2,
                                  final double y2, final double head ) {
        final double a = Math.atan2( y2 - y1, x2 - x1 );
        g2.draw( new Line2D.Double( x1, y1, x2 - ( Math.cos( a ) * head * 0.85 ),
                                    y2 - ( Math.sin( a ) * head * 0.85 ) ) );
        final Path2D h = new Path2D.Double();
        h.moveTo( x2, y2 );
        h.lineTo( x2 - ( Math.cos( a - 0.62 ) * head ), y2 - ( Math.sin( a - 0.62 ) * head ) );
        h.lineTo( x2 - ( Math.cos( a + 0.62 ) * head ), y2 - ( Math.sin( a + 0.62 ) * head ) );
        h.closePath();
        g2.fill( h );
    }

    /**
     * The ORBIT form (the user's chosen reference): a small open circle as the pivot, with a thin arrow curving
     * around it and the arrowhead capping the arc's TERMINUS, tip pointing INTO the gap along the direction of
     * travel. Semantically exact for this button -- the tree stays put and the view rotates around it -- and
     * the pivot ring plus arrowhead keep it from reading as the circular LAYOUT button or the theme sun.
     * The head placement is the part earlier drafts got wrong twice: the line must visibly END in the head like
     * a real arrow. A head at the arc's start (tip pointing along where the line continues) runs the stroke
     * straight through the head, so it reads as a lump partway along the ring.
     * Clockwise: the gap at the lower left, the head at the bottom pointing left (which IS the clockwise
     * direction at 6 o'clock). Counter-clockwise is the mirror image.
     */
    private void paintRotate( final Graphics2D g2, final int x, final int y, final boolean cw ) {
        g2.setStroke( new BasicStroke( Math.max( 1.0f, _size * 0.075f ), BasicStroke.CAP_BUTT,
                                       BasicStroke.JOIN_ROUND ) );
        final double cx = x + ( _size / 2.0 );
        final double cy = y + ( _size / 2.0 );
        // the pivot: a small OPEN circle at the centre -- the thing that stays put while the view turns
        final double pr = _size * 0.12;
        g2.draw( new java.awt.geom.Ellipse2D.Double( cx - pr, cy - pr, 2 * pr, 2 * pr ) );
        // The orbit arc. Screen angles (y down, so increasing = clockwise on screen); Arc2D measures the other
        // way round, hence the negations. The 60-degree gap sits at the lower LEFT for clockwise (lower right
        // for counter-clockwise); the arc sets off from the gap's far edge, runs the long way round, and
        // TERMINATES at the gap's near edge -- where the head caps it, tip into the gap.
        final double r = _size * 0.37;
        final double start = cw ? 170 : 10;               // the gap's far edge, where the arc sets off
        final double sweep = cw ? 300 : -300;             // the long way round
        final double terminus = start + sweep;            // 110 (CW) / 70 (CCW): just off 6 o'clock
        g2.draw( new Arc2D.Double( cx - r, cy - r, 2 * r, 2 * r, -start, -sweep, Arc2D.OPEN ) );
        final double a = Math.toRadians( terminus );
        final double px = cx + ( r * Math.cos( a ) );
        final double py = cy + ( r * Math.sin( a ) );
        final double tangent = a + ( cw ? ( Math.PI / 2 ) : ( -Math.PI / 2 ) ); // the way it is travelling
        // the head's BASE sits on the terminus, square to the travel; only the tip reaches past it, into the gap
        final double len = _size * 0.26;
        final double half = _size * 0.13;
        final double nx = -Math.sin( tangent );
        final double ny = Math.cos( tangent );
        final Path2D head = new Path2D.Double();
        head.moveTo( px + ( Math.cos( tangent ) * len ), py + ( Math.sin( tangent ) * len ) ); // the tip
        head.lineTo( px + ( nx * half ), py + ( ny * half ) );
        head.lineTo( px - ( nx * half ), py - ( ny * half ) );
        head.closePath();
        g2.fill( head );
    }

    /**
     * An arrow pointing LEFT, toward the root; {@code to_bar} adds the bar it stops against, i.e. "all the way
     * back" rather than "one step back".
     */
    private void paintBackArrow( final Graphics2D g2, final int x, final int y, final boolean to_bar ) {
        g2.setStroke( new BasicStroke( Math.max( 1.1f, _size * 0.11f ), BasicStroke.CAP_BUTT,
                                       BasicStroke.JOIN_MITER ) );
        final double cy = y + ( _size / 2.0 );
        final double tip = x + ( _size * ( to_bar ? 0.30 : 0.14 ) );
        final double tail = x + ( _size * 0.88 );
        g2.draw( new Line2D.Double( tip + ( _size * 0.16 ), cy, tail, cy ) );
        final Path2D head = new Path2D.Double();
        head.moveTo( tip, cy );
        head.lineTo( tip + ( _size * 0.26 ), cy - ( _size * 0.22 ) );
        head.lineTo( tip + ( _size * 0.26 ), cy + ( _size * 0.22 ) );
        head.closePath();
        g2.fill( head );
        if ( to_bar ) {
            g2.draw( new Line2D.Double( x + ( _size * 0.14 ), y + ( _size * 0.18 ), x + ( _size * 0.14 ),
                                        y + ( _size * 0.82 ) ) );
        }
    }

    /** The collapsed-clade triangle Archaeopteryx actually draws, with the tips it opens back out into. */
    private void paintUncollapse( final Graphics2D g2, final int x, final int y ) {
        g2.setStroke( new BasicStroke( Math.max( 1.0f, _size * 0.085f ), BasicStroke.CAP_ROUND,
                                       BasicStroke.JOIN_ROUND ) );
        final Path2D tri = new Path2D.Double(); // apex at the parent, base toward the tips (as on screen)
        tri.moveTo( x + ( _size * 0.06 ), y + ( _size * 0.50 ) );
        tri.lineTo( x + ( _size * 0.44 ), y + ( _size * 0.16 ) );
        tri.lineTo( x + ( _size * 0.44 ), y + ( _size * 0.84 ) );
        tri.closePath();
        g2.fill( tri );
        for ( final double ty : new double[] { 0.20, 0.50, 0.80 } ) { // the tips it fans back out into
            g2.draw( new Line2D.Double( x + ( _size * 0.56 ), y + ( _size * ty ), x + ( _size * 0.94 ),
                                        y + ( _size * ty ) ) );
        }
    }
}
