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

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.font.TextAttribute;
import java.awt.geom.RoundRectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import javax.swing.UIManager;

import org.forester.archaeopteryx.NodeHoverText.Row;

/**
 * The node hover card, drawn straight onto the tree canvas as its last overlay: a rounded, slightly translucent
 * card in the UI theme's colours with a hairline border and a soft shadow, small uppercase muted section headings,
 * and key/value rows with a muted key column -- the Archaeopteryx.js tooltip, disc for disc. It is NOT a popup
 * window: it is state the canvas paints, so it can be nowhere but on the canvas and dies with every repaint that
 * does not draw it -- the structural cure for the old popup that lingered on the desktop. Layout is measured once
 * per node ({@link #NodeHoverCard}); {@link #paint} is then a handful of fills.
 */
final class NodeHoverCard {

    /** The JS card's proportions at a 13 px base font; everything scales with the font. */
    static final int   MAX_WIDTH_PX      = 320;
    static final float ALPHA             = 0.97f;
    static final int   FADE_IN_MS        = 150;
    /** Offset of the card's top-left corner from the pointer (the JS uses +14/+14). */
    static final int   POINTER_OFFSET    = 14;
    private static final int   CORNER    = 8;
    private static final int   PAD_X     = 10;
    private static final int   PAD_Y     = 7;
    private static final int   KEY_GAP   = 10;
    private static final float KEY_FRACTION = 0.42f;
    private static final int   ROW_GAP   = 2;
    private static final int   HEAD_GAP  = 6;

    /** One laid-out text line (a heading, a key, or one wrapped fragment of a value). */
    private static final class Line {

        final String text;
        final int    x;
        final int    baseline;
        final Font   font;
        final Color  color;

        Line( final String text, final int x, final int baseline, final Font font, final Color color ) {
            this.text = text;
            this.x = x;
            this.baseline = baseline;
            this.font = font;
            this.color = color;
        }
    }

    private final List<Line> _lines = new ArrayList<>();
    private final int        _width;
    private final int        _height;
    private final Color      _background;
    private final Color      _border;
    private final float      _scale;

    /**
     * Lays the rows out for {@code base} (the UI label font) using {@code fm_of} to measure -- a component's
     * {@code getFontMetrics}, so this works on and off screen. {@code dark} picks the palette.
     */
    NodeHoverCard( final List<Row> rows, final Font base, final java.util.function.Function<Font, FontMetrics> fm_of,
                   final boolean dark ) {
        _scale = base.getSize2D() / 13f;
        final Font value_font = base.deriveFont( Font.PLAIN, Math.max( 9f, base.getSize2D() * 0.9f ) );
        final Map<TextAttribute, Object> attrs = new HashMap<>();
        attrs.put( TextAttribute.TRACKING, 0.07f );
        final Font head_font = base.deriveFont( Font.BOLD, Math.max( 8f, base.getSize2D() * 0.72f ) )
                .deriveFont( attrs );
        final Color ink = dark ? new Color( 0xE8, 0xED, 0xF2 ) : new Color( 0x1E, 0x2A, 0x35 );
        final Color muted = dark ? new Color( 0xA8, 0xB4, 0xC0 ) : new Color( 0x6B, 0x7A, 0x89 );
        final Color faint = dark ? new Color( 0x8A, 0x98, 0xA6 ) : new Color( 0x93, 0xA3, 0xB2 );
        _background = dark ? new Color( 0x2B, 0x2F, 0x35 ) : Color.WHITE;
        _border = dark ? new Color( 0x4A, 0x52, 0x5B ) : new Color( 0xCA, 0xD6, 0xE1 );
        final FontMetrics vfm = fm_of.apply( value_font );
        final FontMetrics hfm = fm_of.apply( head_font );
        final int max_w = Math.round( MAX_WIDTH_PX * _scale );
        final int pad_x = Math.round( PAD_X * _scale );
        final int pad_y = Math.round( PAD_Y * _scale );
        final int key_gap = Math.round( KEY_GAP * _scale );
        // key column: the widest key, capped at a fraction of the card
        int key_w = 0;
        int longest_value = 0;
        for( final Row r : rows ) {
            if ( !r.isHeading() ) {
                key_w = Math.max( key_w, vfm.stringWidth( r.key ) );
                longest_value = Math.max( longest_value, vfm.stringWidth( r.value ) );
            }
        }
        final int inner_max = max_w - ( 2 * pad_x );
        key_w = Math.min( key_w, Math.round( inner_max * KEY_FRACTION ) );
        final int value_x = pad_x + key_w + key_gap;
        final int value_max_w = inner_max - key_w - key_gap;
        // width: as wide as the content needs, up to the cap (the JS "width:max-content; max-width")
        int inner_w = 0;
        for( final Row r : rows ) {
            if ( r.isHeading() ) {
                inner_w = Math.max( inner_w, hfm.stringWidth( r.key.toUpperCase( Locale.ROOT ) ) );
            }
            else {
                inner_w = Math.max( inner_w, key_w + key_gap + Math.min( vfm.stringWidth( r.value ), value_max_w ) );
            }
        }
        inner_w = Math.min( inner_w, inner_max );
        int y = pad_y;
        boolean first = true;
        for( final Row r : rows ) {
            if ( r.isHeading() ) {
                y += first ? 0 : Math.round( HEAD_GAP * _scale );
                y += hfm.getAscent();
                _lines.add( new Line( r.key.toUpperCase( Locale.ROOT ), pad_x, y, head_font, faint ) );
                y += hfm.getDescent() + Math.round( ROW_GAP * _scale );
            }
            else {
                final List<String> frags = wrap( r.value, vfm, value_max_w );
                y += vfm.getAscent();
                _lines.add( new Line( clip( r.key, vfm, key_w ), pad_x, y, value_font, muted ) );
                for( int i = 0; i < frags.size(); ++i ) {
                    if ( i > 0 ) {
                        y += vfm.getHeight();
                    }
                    _lines.add( new Line( frags.get( i ), value_x, y, value_font, ink ) );
                }
                y += vfm.getDescent() + Math.round( ROW_GAP * _scale );
            }
            first = false;
        }
        _width = inner_w + ( 2 * pad_x );
        _height = y - Math.round( ROW_GAP * _scale ) + pad_y;
    }

    int getWidth() {
        return _width;
    }

    int getHeight() {
        return _height;
    }

    /** Where to put the card for a pointer at {@code pointer} so it stays inside {@code visible}: below-right of
     *  the pointer by default, flipped above / to the left when that would leave the visible area. Pure. */
    Rectangle placement( final int pointer_x, final int pointer_y, final Rectangle visible ) {
        final int off = Math.round( POINTER_OFFSET * _scale );
        int x = pointer_x + off;
        int y = pointer_y + off;
        if ( x + _width > visible.x + visible.width ) {
            x = pointer_x - off - _width;
        }
        if ( y + _height > visible.y + visible.height ) {
            y = pointer_y - off - _height;
        }
        x = Math.max( visible.x, Math.min( x, visible.x + visible.width - _width ) );
        y = Math.max( visible.y, Math.min( y, visible.y + visible.height - _height ) );
        return new Rectangle( x, y, _width, _height );
    }

    /** Paints the card with its top-left corner at ({@code x},{@code y}), at {@code alpha} (the fade-in). */
    void paint( final Graphics2D g0, final int x, final int y, final float alpha ) {
        final Graphics2D g = (Graphics2D) g0.create();
        try {
            g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
            g.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            g.setRenderingHint( RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON );
            final float a = Math.max( 0f, Math.min( 1f, alpha ) );
            final int corner = Math.round( CORNER * _scale );
            // shadow: three soft rings, offset downwards, faintest outermost (the JS box-shadow)
            final int[] spread = { 6, 3, 1 };
            final float[] shade = { 0.05f, 0.07f, 0.10f };
            for( int i = 0; i < spread.length; ++i ) {
                final int s = Math.round( spread[ i ] * _scale );
                g.setComposite( AlphaComposite.getInstance( AlphaComposite.SRC_OVER, shade[ i ] * a ) );
                g.setColor( Color.BLACK );
                g.fill( new RoundRectangle2D.Float( x - s, y - s + Math.round( 3 * _scale ), _width + ( 2 * s ),
                                                    _height + ( 2 * s ), corner + s, corner + s ) );
            }
            g.setComposite( AlphaComposite.getInstance( AlphaComposite.SRC_OVER, ALPHA * a ) );
            g.setColor( _background );
            final RoundRectangle2D card = new RoundRectangle2D.Float( x, y, _width, _height, corner, corner );
            g.fill( card );
            g.setComposite( AlphaComposite.getInstance( AlphaComposite.SRC_OVER, a ) );
            g.setColor( _border );
            g.draw( card );
            for( final Line l : _lines ) {
                g.setFont( l.font );
                g.setColor( l.color );
                g.drawString( l.text, x + l.x, y + l.baseline );
            }
        }
        finally {
            g.dispose();
        }
    }

    /** Greedy word wrap into lines no wider than {@code max_w}; a single over-long word is clipped with "…". */
    static List<String> wrap( final String text, final FontMetrics fm, final int max_w ) {
        final List<String> out = new ArrayList<>();
        if ( ( text == null ) || text.isEmpty() ) {
            out.add( "" );
            return out;
        }
        if ( fm.stringWidth( text ) <= max_w ) {
            out.add( text );
            return out;
        }
        final StringBuilder line = new StringBuilder();
        for( final String word : text.split( " " ) ) {
            final String candidate = ( line.length() == 0 ) ? word : line + " " + word;
            if ( fm.stringWidth( candidate ) <= max_w ) {
                line.setLength( 0 );
                line.append( candidate );
            }
            else {
                if ( line.length() > 0 ) {
                    out.add( line.toString() );
                }
                line.setLength( 0 );
                line.append( clip( word, fm, max_w ) );
            }
            if ( out.size() >= 6 ) { // a card is a glance, not a document (the node window has the rest)
                out.set( out.size() - 1, clip( out.get( out.size() - 1 ) + " …", fm, max_w ) );
                return out;
            }
        }
        if ( line.length() > 0 ) {
            out.add( line.toString() );
        }
        return out;
    }

    /** {@code s} shortened with "…" to fit {@code max_w}. */
    static String clip( final String s, final FontMetrics fm, final int max_w ) {
        if ( fm.stringWidth( s ) <= max_w ) {
            return s;
        }
        final String ell = "…";
        int end = s.length();
        while ( ( end > 0 ) && ( fm.stringWidth( s.substring( 0, end ) + ell ) > max_w ) ) {
            --end;
        }
        return s.substring( 0, end ) + ell;
    }

    /** Whether the UI theme is dark (drives the card palette). */
    static boolean isDarkTheme() {
        final Color bg = UIManager.getColor( "Panel.background" );
        if ( bg == null ) {
            return false;
        }
        return ( ( 0.299 * bg.getRed() ) + ( 0.587 * bg.getGreen() ) + ( 0.114 * bg.getBlue() ) ) < 128;
    }
}
