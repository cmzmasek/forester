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
import java.awt.LinearGradientPaint;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.Point2D;

import javax.swing.JComponent;
import javax.swing.Timer;
import javax.swing.UIManager;

/**
 * A prominent, menu-bar counter of the currently highlighted nodes -- the combined "Found / Selected" pool. In
 * Archaeopteryx a search hit (box A / box B) and a manual selection are the SAME thing (selection reuses found-set 0),
 * so this shows ONE number: the count of distinct highlighted nodes. It hides itself when that count is zero, and on
 * each change it plays a brief two-pass "run-through" shimmer sweep across the text (clipped to the glyphs) to draw the
 * eye to the update, then settles static. Display-only chrome; the count is supplied by {@code MainFrame} from the
 * current tree panel's found sets.
 */
final class FoundSelectedCounter extends JComponent {

    private static final int    PAD_X   = 10;                  // horizontal padding around the text
    private static final int    BAND    = 26;                  // shimmer band width (px)
    private static final int    TICK_MS = 25;                  // animation frame interval
    private static final int    STEP_PX = 8;                   // shimmer advance per frame (lower = slower sweep)
    private static final int    PASSES  = 2;                   // run-throughs per change (the "2-pass sweep")
    private static final int    GAP     = 10;                  // ticks the text rests between passes (~250ms)
    private static final String PREFIX  = "Found / Selected: ";

    private int          _total;
    private int          _a;
    private int          _b;
    private int          _both;
    private boolean      _showing;
    private final Timer  _timer;
    private double       _sweep_x;
    private int          _pass;
    private int          _gap;      // remaining rest ticks between two passes
    private boolean      _sweeping;

    FoundSelectedCounter() {
        setOpaque( false );
        setVisible( false );
        _timer = new Timer( TICK_MS, e -> onTick() );
        _timer.setCoalesce( true );
    }

    /**
     * Sets the highlighted-node counts and (when they changed, or the counter was hidden) triggers the sweep.
     * {@code total} = distinct highlighted nodes (union of the two found sets); {@code a}/{@code b}/{@code both} are
     * the per-set breakdown for the tooltip. {@code a_is_search} tells whether found-set 0 currently holds Search-A
     * hits ({@code true}) or a manual SELECTION ({@code false}, i.e. Search A has no query) -- Archaeopteryx reuses
     * found-set 0 for both, so this decides whether that count is labelled "A" or "Selected". When {@code combine_label}
     * is non-null (the two boxes are combined via "Combine A & B"), the tooltip shows that description (e.g. "A AND B")
     * instead of the per-box breakdown, since the highlight is one combined set. A {@code total} of 0 hides the counter.
     */
    void setCounts( final int total, final int a, final int b, final int both, final boolean a_is_search,
                    final String combine_label ) {
        if ( total <= 0 ) {
            _total = 0;
            _a = 0;
            _b = 0;
            _both = 0;
            if ( _showing ) {
                _showing = false;
                stopSweep();
                setVisible( false );
            }
            return;
        }
        final boolean changed = ( total != _total ) || !_showing;
        _total = total;
        _a = a;
        _b = b;
        _both = both;
        setToolTipText( tooltipText( total, a, b, both, a_is_search, combine_label ) );
        _showing = true;
        if ( !isVisible() ) {
            setVisible( true );
        }
        revalidate(); // the text width may have changed -> let the menu bar re-lay-out
        if ( changed ) {
            startSweep();
        }
        repaint();
    }

    // ---- text (pure, testable) ---------------------------------------------------------------------------------

    /** The label text, e.g. {@code "Found / Selected: 42"}. */
    static String labelText( final int total ) {
        return PREFIX + total;
    }

    /** The tooltip: the distinct total plus the per-set breakdown, stressing that hits and selections are one pool.
     *  Found-set 0 is labelled "A" when Search A holds a query, else "Selected" (a manual selection). Each part is
     *  shown only when non-zero (a solo search / selection shows just its own count). */
    static String tooltipText( final int total, final int a, final int b, final int both, final boolean a_is_search,
                               final String combine_label ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( total ).append( total == 1 ? " highlighted node" : " highlighted nodes" );
        sb.append( " (search hits + selection)" );
        if ( combine_label != null ) {
            return sb.append( " — " ).append( combine_label ).toString(); // combined result: no per-box breakdown
        }
        if ( ( a > 0 ) || ( b > 0 ) ) {
            sb.append( ": " );
            boolean first = true;
            if ( a > 0 ) {
                sb.append( a_is_search ? "A: " : "Selected: " ).append( a );
                first = false;
            }
            if ( b > 0 ) {
                if ( !first ) {
                    sb.append( ", " );
                }
                sb.append( "B: " ).append( b );
            }
            if ( both > 0 ) {
                sb.append( " (" ).append( both ).append( " in both)" );
            }
        }
        return sb.toString();
    }

    private String fullText() {
        return labelText( _total );
    }

    private Font paintFont() {
        final Font f = getFont();
        return ( f == null ) ? new Font( Font.SANS_SERIF, Font.BOLD, 12 ) : f.deriveFont( Font.BOLD );
    }

    @Override
    public Dimension getPreferredSize() {
        final Font f = paintFont();
        final FontMetrics fm = getFontMetrics( f );
        return new Dimension( fm.stringWidth( fullText() ) + ( 2 * PAD_X ), fm.getHeight() + 4 );
    }

    // ---- animation ---------------------------------------------------------------------------------------------

    private void startSweep() {
        _sweeping = true;
        _pass = 0;
        _gap = 0;
        _sweep_x = -BAND;
        _timer.start();
    }

    private void stopSweep() {
        _sweeping = false;
        _timer.stop();
    }

    private void onTick() {
        if ( _gap > 0 ) {
            _gap--;      // resting between the two passes -- the band sits off-screen, so the text shows static
            repaint();
            return;
        }
        _sweep_x += STEP_PX;
        if ( _sweep_x > ( getWidth() + BAND ) ) {
            _pass++;
            if ( _pass >= PASSES ) {
                stopSweep();
                repaint();
                return;
            }
            _sweep_x = -BAND; // next pass starts off the left edge
            _gap = GAP;       // pause so the two passes read as distinct, not one continuous sweep
        }
        repaint();
    }

    @Override
    public void removeNotify() {
        stopSweep(); // don't leave a timer running after the frame is torn down
        super.removeNotify();
    }

    // ---- painting ----------------------------------------------------------------------------------------------

    @Override
    protected void paintComponent( final Graphics g ) {
        final Graphics2D g2 = (Graphics2D) g.create();
        try {
            g2.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
            g2.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            final Font f = paintFont();
            g2.setFont( f );
            final FontMetrics fm = g2.getFontMetrics();
            final String text = fullText();
            final int x = PAD_X;
            final int baseline = ( ( getHeight() - fm.getHeight() ) / 2 ) + fm.getAscent();
            g2.setColor( baseColor() );
            g2.drawString( text, x, baseline );
            if ( _sweeping ) {
                paintSweep( g2, f, text, x, baseline );
            }
        }
        finally {
            g2.dispose();
        }
    }

    /** The run-through: a bright band, clipped to the text glyphs, swept horizontally. */
    private void paintSweep( final Graphics2D g2, final Font f, final String text, final int x, final int baseline ) {
        final FontRenderContext frc = g2.getFontRenderContext();
        final GlyphVector gv = f.createGlyphVector( frc, text );
        final Shape outline = gv.getOutline( x, baseline );
        final Shape old_clip = g2.getClip();
        g2.clip( outline );
        final float x0 = (float) ( _sweep_x - ( BAND / 2.0 ) );
        final float x1 = (float) ( _sweep_x + ( BAND / 2.0 ) );
        if ( x1 > x0 ) {
            final Color sh = shimmerColor();
            final Color trans = new Color( sh.getRed(), sh.getGreen(), sh.getBlue(), 0 );
            g2.setPaint( new LinearGradientPaint( new Point2D.Float( x0, 0 ), new Point2D.Float( x1, 0 ),
                                                  new float[] { 0f, 0.5f, 1f }, new Color[] { trans, sh, trans } ) );
            g2.fillRect( 0, 0, getWidth(), getHeight() );
        }
        g2.setClip( old_clip );
    }

    private static Color baseColor() {
        Color c = UIManager.getColor( "Component.accentColor" );
        if ( c == null ) {
            c = UIManager.getColor( "Label.foreground" );
        }
        return ( c == null ) ? Color.GRAY : c;
    }

    /** A brighter version of the base colour -- the moving highlight reads as a shine sweeping through the text. */
    private static Color shimmerColor() {
        return blend( baseColor(), Color.WHITE, 0.65f );
    }

    private static Color blend( final Color a, final Color b, final float t ) {
        final float s = 1f - t;
        return new Color( Math.round( ( a.getRed() * s ) + ( b.getRed() * t ) ),
                          Math.round( ( a.getGreen() * s ) + ( b.getGreen() * t ) ),
                          Math.round( ( a.getBlue() * s ) + ( b.getBlue() * t ) ) );
    }

    // ---- test hooks --------------------------------------------------------------------------------------------

    boolean isShowingForTest() {
        return _showing && isVisible();
    }

    boolean isSweepingForTest() {
        return _sweeping;
    }

    int getTotalForTest() {
        return _total;
    }

    String getTooltipForTest() {
        return getToolTipText();
    }

    /** Drives one animation frame (the Timer is EDT-based; a test advances it deterministically). */
    void tickForTest() {
        onTick();
    }

    /** Runs the sweep to completion (both passes) so a test can assert it settles static; returns the number of
     *  animation frames it took (a single-pass regression would be far fewer). */
    int runSweepToEndForTest() {
        int ticks = 0;
        while ( _sweeping && ( ticks < 100000 ) ) {
            onTick();
            ticks++;
        }
        return ticks;
    }
}
