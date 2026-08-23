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

package org.forester.archaeopteryx.phylogeny.data;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Rectangle2D;
import java.awt.geom.RoundRectangle2D;
import java.io.IOException;
import java.io.Writer;
import java.math.BigDecimal;
import java.util.Map;
import java.util.SortedMap;

import org.forester.archaeopteryx.AptxConstants;
import org.forester.archaeopteryx.AptxUtil;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.util.ForesterUtil;

public final class RenderableDomainArchitecture extends DomainArchitecture implements RenderablePhylogenyData {

    private final static int          E_VALUE_THRESHOLD_EXP_DEFAULT = 0;
    private final static BasicStroke  STROKE_1                      = new BasicStroke( 1f );
    // The flat modern look (fixed style): a rounded gradient body, a thin hue border, and a soft drop-shadow. The
    // shadow (and the optional glow) are STEPPED -- a few translucent offset/outset copies rather than a Gaussian blur
    // -- so they render IDENTICALLY on screen and in vector export (PDF / SVG / EPS have no blur). GradientPaint is
    // vector-safe (verified: OpenPDF renders a native shading, VectorGraphics2D tessellates it correctly).
    // Small radius on purpose: domain boxes are SHORT (height capped ~11px), so a larger radius rounds the short sides
    // into a pill/oval. 2px reads as a rounded RECTANGLE at every height (r is still clamped to min(w,h)/2 below).
    private final static float        CORNER_RADIUS                 = 2f;
    private final static float        GRAD_LIGHTEN_TOP              = 0.12f; // top edge toward white
    private final static float        GRAD_DARKEN_BOTTOM            = 0.10f; // bottom edge toward black
    private final static float        BORDER_DARKEN                 = 0.24f;
    // { dx, dy, alpha } per drop-shadow layer -- darkest/closest first, softest/farthest last (a stepped soft shadow).
    private final static float[][]    SHADOW_LAYERS                 = { { 0.4f, 0.7f, 40 }, { 0.9f, 1.5f, 28 },
            { 1.6f, 2.5f, 17 } };
    // { outset, alpha } per glow layer (a hue halo; off by default).
    private final static float[][]    GLOW_LAYERS                   = { { 3.2f, 20 }, { 1.6f, 34 } };
    private static Map<String, Color> _domain_colors;
    private final DomainArchitecture  _domain_structure;
    private int                       _e_value_threshold_exp        = E_VALUE_THRESHOLD_EXP_DEFAULT;
    private final Rectangle2D         _rectangle                    = new Rectangle2D.Float();
    private final RoundRectangle2D    _round                        = new RoundRectangle2D.Float();
    private float                     _rendering_factor_width        = 1;
    private float                     _rendering_height              = 0;

    public RenderableDomainArchitecture( final DomainArchitecture domain_structure ) {
        _domain_structure = domain_structure;
    }

    public static void setColorMap( final Map<String, Color> domain_colors ) {
        _domain_colors = domain_colors;
    }

    @Override
    public StringBuffer asSimpleText() {
        return _domain_structure.asSimpleText();
    }

    @Override
    public StringBuffer asText() {
        return _domain_structure.asText();
    }

    @Override
    public PhylogenyData copy() {
        return _domain_structure.copy();
    }

    /**
     * Draw one domain box in the flat modern style: a soft stepped drop-shadow, an optional hue glow, a rounded body
     * with a mild vertical gradient, and a thin hue border. Pure Graphics2D primitives (translucent rounded rects +
     * a GradientPaint), so it is WYSIWYG on screen and in every vector export.
     */
    private void drawDomainFlat( final double x, final double y, final double width, final double height,
                                 final Color base, final boolean glow, final Graphics2D g ) {
        if ( ( width <= 0 ) || !Double.isFinite( width ) ) { // non-finite guards a degenerate factor (all-zero-length architectures)
            return;
        }
        final double r = Math.min( CORNER_RADIUS, Math.min( width, height ) / 2.0 );
        final double arc = 2 * r;
        // stepped soft drop-shadow (behind the box)
        for( final float[] s : SHADOW_LAYERS ) {
            g.setColor( new Color( 8, 18, 21, (int) s[ 2 ] ) );
            _round.setRoundRect( x + s[ 0 ], y + s[ 1 ], width, height, arc, arc );
            g.fill( _round );
        }
        // optional hue glow (behind the box)
        if ( glow ) {
            for( final float[] gl : GLOW_LAYERS ) {
                final float o = gl[ 0 ];
                g.setColor( new Color( base.getRed(), base.getGreen(), base.getBlue(), (int) gl[ 1 ] ) );
                _round.setRoundRect( x - o, y - o, width + ( 2 * o ), height + ( 2 * o ), arc + ( 2 * o ),
                                     arc + ( 2 * o ) );
                g.fill( _round );
            }
        }
        // rounded body with a mild vertical gradient
        g.setPaint( new GradientPaint( 0f, (float) y, lighten( base, GRAD_LIGHTEN_TOP ), 0f,
                                       (float) ( y + height ), darken( base, GRAD_DARKEN_BOTTOM ) ) );
        _round.setRoundRect( x, y, width, height, arc, arc );
        g.fill( _round );
        // thin hue border
        g.setColor( darken( base, BORDER_DARKEN ) );
        g.setStroke( STROKE_1 );
        g.draw( _round );
    }

    /**
     * The colour the renderer uses for a domain NAME: the palette entry (see
     * {@code AptxUtil.assignDomainPalette}) if present, else a deterministic fallback colour, cached. Public + static
     * so a legend (painted by the {@code TreePanel}) can match the drawn boxes exactly.
     */
    public static Color colorFor( final String name ) {
        if ( name == null ) {
            return Color.GRAY; // a nameless domain (guards calculateColorFromString(null) -> NPE)
        }
        if ( _domain_colors == null ) {
            _domain_colors = new java.util.HashMap<String, Color>();
        }
        Color c = _domain_colors.get( name );
        if ( c == null ) {
            c = AptxUtil.calculateColorFromString( name, false );
            if ( c == null ) {
                c = Color.GRAY;
            }
            _domain_colors.put( name, c );
        }
        return c;
    }

    private static Color lighten( final Color c, final float t ) {
        return new Color( c.getRed() + Math.round( ( 255 - c.getRed() ) * t ),
                          c.getGreen() + Math.round( ( 255 - c.getGreen() ) * t ),
                          c.getBlue() + Math.round( ( 255 - c.getBlue() ) * t ) );
    }

    private static Color darken( final Color c, final float t ) {
        return new Color( Math.round( c.getRed() * ( 1 - t ) ), Math.round( c.getGreen() * ( 1 - t ) ),
                          Math.round( c.getBlue() * ( 1 - t ) ) );
    }

    /** Near-black or white label ink, whichever reads on {@code c} (by relative luminance). */
    private static Color contrastInk( final Color c ) {
        final double lum = ( ( 0.2126 * c.getRed() ) + ( 0.7152 * c.getGreen() ) + ( 0.0722 * c.getBlue() ) ) / 255.0;
        return lum > 0.55 ? new Color( 20, 26, 29 ) : Color.WHITE;
    }

    @Override
    public ProteinDomain getDomain( final int i ) {
        return _domain_structure.getDomain( i );
    }

    @Override
    public SortedMap<BigDecimal, ProteinDomain> getDomains() {
        return _domain_structure.getDomains();
    }

    @Override
    public int getNumberOfDomains() {
        return _domain_structure.getNumberOfDomains();
    }

    @Override
    public Dimension getOriginalSize() {
        return new Dimension( _domain_structure.getTotalLength(), ForesterUtil.roundToInt( _rendering_height ) );
    }

    @Override
    public Object getParameter() {
        return Integer.valueOf( _e_value_threshold_exp );
    }

    public float getRenderingFactorWidth() {
        return _rendering_factor_width;
    }

    @Override
    public Dimension getRenderingSize() {
        return new Dimension( ForesterUtil.roundToInt( _domain_structure.getTotalLength() * getRenderingFactorWidth() ),
                              ForesterUtil.roundToInt( _rendering_height ) );
    }

    @Override
    public int getTotalLength() {
        return _domain_structure.getTotalLength();
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        return _domain_structure.isEqual( data );
    }

    @Override
    public void render( final float x1,
                        final float y1,
                        final Graphics2D g,
                        final TreePanel tree_panel,
                        final boolean to_pdf ) {
        render( x1, y1, g, tree_panel, to_pdf, true );
    }

    /**
     * As {@link #render(float, float, Graphics2D, TreePanel, boolean)}, but {@code draw_labels} can suppress the
     * per-domain name labels. Used by the vertical (root-top/bottom) orientation, where the architecture boxes ride
     * the canvas rotation into a thin vertical track per tip and upright names would collide with neighbouring tips.
     */
    public void render( final float x1,
                        final float y1,
                        final Graphics2D g,
                        final TreePanel tree_panel,
                        final boolean to_pdf,
                        final boolean draw_labels ) {
        final float f = getRenderingFactorWidth();
        final float y = y1 + ( _rendering_height / 2 );
        final float start = x1 + 20;
        final Stroke s = g.getStroke();
        g.setStroke( STROKE_1 );
        g.setColor( to_pdf ? AptxConstants.DOMAIN_BASE_COLOR_FOR_PDF
                           : tree_panel.getTreeColorSet().getDomainBaseColor() );
        _rectangle.setFrame( start, y - 0.5, _domain_structure.getTotalLength() * f, 1 );
        g.fill( _rectangle );
        // "Labels on domains" draws the name centred inside each box; "Legend" and "No labels" draw none here (the
        // legend is a separate, draggable, E-value-aware slot painted by the TreePanel).
        final boolean glow = tree_panel.getMainPanel().getOptions().isShowDomainGlow();
        final boolean on_domain_labels = draw_labels
                && tree_panel.getMainPanel().getOptions().isDomainLabelsOnDomains()
                && ( tree_panel.getMainPanel().getTreeFontSet().getFontMetricsSmall().getHeight() > 4 );
        for( int i = 0; i < _domain_structure.getDomains().size(); ++i ) {
            final ProteinDomain d = _domain_structure.getDomain( i );
            if ( d.getConfidence() <= Math.pow( 10, _e_value_threshold_exp ) ) {
                final float xa = start + ( d.getFrom() * f );
                final float xb = xa + ( d.getLength() * f );
                final Color base = colorFor( d.getName() );
                drawDomainFlat( xa, y1, xb - xa, _rendering_height, base, glow, g );
                if ( on_domain_labels && ( d.getName() != null ) ) { // a nameless domain still draws its box, just no label
                    final FontMetrics fm = tree_panel.getMainPanel().getTreeFontSet().getFontMetricsSmall();
                    final int tw = fm.stringWidth( d.getName() );
                    if ( tw <= ( ( xb - xa ) - 4 ) ) { // fits centred inside the box, else drop it
                        g.setFont( tree_panel.getMainPanel().getTreeFontSet().getSmallFont() );
                        g.setColor( contrastInk( base ) );
                        final float tx = xa + ( ( ( xb - xa ) - tw ) / 2f );
                        final float ty = y1 + ( ( _rendering_height + fm.getAscent() - fm.getDescent() ) / 2f );
                        g.drawString( d.getName(), tx, ty );
                    }
                }
            }
        }
        g.setStroke( s );
    }

    @Override
    public void setParameter( final double e_value_threshold_exp ) {
        _e_value_threshold_exp = ( int ) e_value_threshold_exp;
    }

    public void setRenderingFactorWidth( final float rendering_factor_width ) {
        _rendering_factor_width = rendering_factor_width;
    }

    @Override
    public void setRenderingHeight( final float rendering_height ) {
        _rendering_height = rendering_height;
    }

    @Override
    public StringBuffer toNHX() {
        return _domain_structure.toNHX();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        _domain_structure.toPhyloXML( writer, level, indentation );
    }
}
