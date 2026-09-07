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
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.MouseEvent;

import javax.swing.JComponent;

import org.forester.archaeopteryx.TreeFacts.Histogram;

/**
 * A small painted bar histogram (the Tree Properties window's branch-length and support-value distributions):
 * accent-coloured bars on a baseline, the range's minimum and maximum in muted small text underneath, and a
 * per-bar tooltip with the bin's range and count. Replaces the old ASCII-art histogram.
 */
final class HistogramPanel extends JComponent {

    private static final long serialVersionUID = 1L;
    static final int          BAR_AREA_HEIGHT  = 48;
    static final int          GAP              = 2;
    /** Bars wider than this look like blocks, not a distribution. */
    static final int          MAX_WIDTH        = 360;
    private final Histogram   _h;

    HistogramPanel( final Histogram h ) {
        _h = h;
        setOpaque( false );
        final Font base = javax.swing.UIManager.getFont( "Label.font" );
        setFont( ( base != null ) ? base : new Font( Font.SANS_SERIF, Font.PLAIN, 13 ) ); // a bare JComponent has none
        setToolTipText( "" ); // registers with the tooltip manager; the text comes from getToolTipText(MouseEvent)
        final int label_h = getFontMetrics( smallFont() ).getHeight();
        setPreferredSize( new Dimension( 240, BAR_AREA_HEIGHT + 4 + label_h ) );
        setMinimumSize( new Dimension( 120, BAR_AREA_HEIGHT + 4 + label_h ) );
        setMaximumSize( new Dimension( MAX_WIDTH, BAR_AREA_HEIGHT + 4 + label_h ) );
    }

    Histogram histogramForTest() {
        return _h;
    }

    private Font smallFont() {
        final Font f = getFont();
        return f.deriveFont( Math.max( 9f, f.getSize2D() - 2f ) );
    }

    /** The bin under x, or -1. */
    int binAt( final int x ) {
        final int n = _h.counts.length;
        final double bar_w = ( Math.min( getWidth(), MAX_WIDTH ) - ( GAP * ( n - 1 ) ) ) / (double) n;
        if ( ( bar_w <= 0 ) || ( x < 0 ) || ( x >= Math.min( getWidth(), MAX_WIDTH ) ) ) {
            return -1;
        }
        final int i = (int) ( x / ( bar_w + GAP ) );
        return ( i < n ) ? i : -1;
    }

    @Override
    public String getToolTipText( final MouseEvent e ) {
        final int i = binAt( e.getX() );
        if ( i < 0 ) {
            return null;
        }
        final String lo = TreeFacts.number( _h.edge( i ) );
        final String hi = TreeFacts.number( _h.edge( i + 1 ) );
        final int c = _h.counts[ i ];
        return lo + " – " + hi + ": " + TreeFacts.count( c ) + ( c == 1 ? " branch" : " branches" );
    }

    @Override
    protected void paintComponent( final Graphics g ) {
        final Graphics2D g2 = (Graphics2D) g.create();
        try {
            g2.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
            g2.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            final int n = _h.counts.length;
            final int w = Math.min( getWidth(), MAX_WIDTH );
            final double bar_w = ( w - ( GAP * ( n - 1 ) ) ) / (double) n;
            final int max = Math.max( 1, _h.maxCount() );
            final Color bar = FormWidgets.accentColor();
            final Color muted = FormWidgets.mutedColor();
            for( int i = 0; i < n; ++i ) {
                final int x = (int) Math.round( i * ( bar_w + GAP ) );
                final int bw = Math.max( 1, (int) Math.round( bar_w ) );
                final int bh = (int) Math.round( ( BAR_AREA_HEIGHT - 1 ) * ( _h.counts[ i ] / (double) max ) );
                g2.setColor( _h.counts[ i ] == 0 ? muted : bar );
                if ( _h.counts[ i ] == 0 ) {
                    g2.fillRect( x, BAR_AREA_HEIGHT - 1, bw, 1 ); // an empty bin: just its baseline tick
                }
                else {
                    g2.fillRoundRect( x, BAR_AREA_HEIGHT - bh, bw, bh, 3, 3 );
                }
            }
            g2.setColor( FormWidgets.borderColor() );
            g2.fillRect( 0, BAR_AREA_HEIGHT, w, 1 );
            g2.setFont( smallFont() );
            g2.setColor( muted );
            final FontMetrics fm = g2.getFontMetrics();
            final int ty = BAR_AREA_HEIGHT + 3 + fm.getAscent();
            final String lo = TreeFacts.number( _h.min );
            final String hi = TreeFacts.number( _h.max );
            g2.drawString( lo, 0, ty );
            g2.drawString( hi, w - fm.stringWidth( hi ), ty );
        }
        finally {
            g2.dispose();
        }
    }
}
