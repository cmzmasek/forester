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
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JPanel;
import javax.swing.JViewport;
import javax.swing.Scrollable;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;

/**
 * A dedicated, deliberately minimal renderer for a tanglegram: two rectangular cladograms drawn facing each other
 * (left tree root-left, right tree mirrored root-right) with the tips of one linked to the tips of the other by
 * connector lines. It owns independent deep copies of both trees (so it never mutates the source tabs), computes a
 * simple lined-up cladogram layout (tips flush at the facing edge, internal node = midpoint of its first/last child),
 * and draws matched-tip connectors from {@link TanglegramLinker}. Unmatched tips are greyed and get no connector.
 *
 * <p>This is NOT a {@link TreePanel}: it intentionally lacks domains/MSA/annotations and every editing tool -- a
 * tanglegram is a read-only comparison view. Interactive rotation / auto-untangle are deferred to a later increment.
 */
final class TanglegramPanel extends JPanel implements Scrollable {

    private static final long   serialVersionUID    = 1L;
    private static final int    MARGIN              = 14;
    private static final int    TOP_PAD             = 12;
    private static final int    LABEL_GAP           = 6;
    private static final int    MIN_TREE_WIDTH      = 40;
    private static final int    MIN_CONNECTOR_GAP   = 40;
    private static final int    BASE_WIDTH          = 900;
    private static final int    ROW_PIXELS          = 20;
    private static final double ZOOM_STEP           = 1.25;
    private static final double MIN_ZOOM            = 0.2;
    private static final double MAX_ZOOM            = 20.0;
    private final Phylogeny                _left;
    private final Phylogeny                _right;
    private final List<PhylogenyNode>      _left_tips;
    private final List<PhylogenyNode>      _right_tips;
    private final TanglegramLinker.Result  _result;
    private final Set<PhylogenyNode>       _unmatched;
    private final int                      _crossings;
    private double                         _zoom       = 1.0;
    private int                            _laid_out_w = -1;
    private int                            _laid_out_h = -1;

    TanglegramPanel( final Phylogeny left_source, final Phylogeny right_source, final TanglegramLinker.LinkField field ) {
        _left = left_source.copy();
        _right = right_source.copy();
        _left_tips = TanglegramLinker.externalTipsInDisplayOrder( _left );
        _right_tips = TanglegramLinker.externalTipsInDisplayOrder( _right );
        _result = TanglegramLinker.link( _left, _right, field );
        _unmatched = Collections.newSetFromMap( new IdentityHashMap<PhylogenyNode, Boolean>() );
        _unmatched.addAll( _result.getUnmatchedA() );
        _unmatched.addAll( _result.getUnmatchedB() );
        _crossings = computeCrossings();
        setOpaque( true );
    }

    TanglegramLinker.Result getResult() {
        return _result;
    }

    int getCrossingCount() {
        return _crossings;
    }

    int getLeftTipCount() {
        return _left_tips.size();
    }

    int getRightTipCount() {
        return _right_tips.size();
    }

    int getUnmatchedCount() {
        return _unmatched.size();
    }

    TanglegramLinker.LinkField getField() {
        return _result.getField();
    }

    /** Test hook: the label that would be drawn for the i-th left tip (exercises the identity-field fallback). */
    String leftTipLabelForTest( final int i ) {
        return labelFor( _left_tips.get( i ) );
    }

    void zoomIn() {
        setZoom( _zoom * ZOOM_STEP );
    }

    void zoomOut() {
        setZoom( _zoom / ZOOM_STEP );
    }

    /** Scales the tanglegram so the whole (taller) axis fits the current viewport height; falls back to natural
     *  spacing (zoom 1) when not yet in a viewport. Width already tracks the viewport, so this fits the figure. */
    void fit() {
        if ( getParent() instanceof JViewport ) {
            final int viewport_h = getParent().getHeight();
            if ( viewport_h > 0 ) {
                setZoom( (double) viewport_h / naturalHeight() );
                return;
            }
        }
        setZoom( 1.0 );
    }

    private int naturalHeight() {
        final int max_tips = Math.max( _left_tips.size(), _right_tips.size() );
        return Math.max( ( 2 * TOP_PAD ) + ( ROW_PIXELS * max_tips ), 300 );
    }

    private void setZoom( final double zoom ) {
        final double clamped = Math.max( MIN_ZOOM, Math.min( MAX_ZOOM, zoom ) );
        if ( clamped != _zoom ) {
            _zoom = clamped;
            revalidate();
            repaint();
        }
    }

    private int computeCrossings() {
        final Map<PhylogenyNode, Integer> left_index = indexOf( _left_tips );
        final Map<PhylogenyNode, Integer> right_index = indexOf( _right_tips );
        final List<int[]> pairs = new ArrayList<>();
        for( final TanglegramLinker.Link link : _result.getLinks() ) {
            final Integer a = left_index.get( link.getA() );
            final Integer b = right_index.get( link.getB() );
            if ( ( a != null ) && ( b != null ) ) {
                pairs.add( new int[] { a, b } );
            }
        }
        return TanglegramLinker.countCrossings( pairs );
    }

    private static Map<PhylogenyNode, Integer> indexOf( final List<PhylogenyNode> tips ) {
        final Map<PhylogenyNode, Integer> index = new IdentityHashMap<>();
        for( int i = 0; i < tips.size(); i++ ) {
            index.put( tips.get( i ), i );
        }
        return index;
    }

    @Override
    protected void paintComponent( final Graphics g ) {
        super.paintComponent( g );
        final Graphics2D g2 = (Graphics2D) g.create();
        try {
            g2.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
            g2.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            g2.setFont( getFont() );
            final int w = getWidth();
            final int h = getHeight();
            if ( ( w <= 0 ) || ( h <= 0 ) || _left_tips.isEmpty() || _right_tips.isEmpty() ) {
                return;
            }
            final FontMetrics fm = g2.getFontMetrics();
            final double avail = w - ( 2.0 * MARGIN );
            // reserve space for tip labels, but never let them starve the two trees or the central connector gap
            final double label_cap = Math.max( 0.0, ( avail - ( 2.0 * MIN_TREE_WIDTH ) - MIN_CONNECTOR_GAP ) / 2.0 );
            final double w_left = Math.min( maxLabelWidth( _left_tips, fm ), label_cap );
            final double w_right = Math.min( maxLabelWidth( _right_tips, fm ), label_cap );
            final double tree_w = Math.max( MIN_TREE_WIDTH, ( avail - w_left - w_right - MIN_CONNECTOR_GAP ) / 2.0 );
            final double avail_h = h - ( 2.0 * TOP_PAD );
            // the layout (topology, tip spacing, depth) is a pure function of the panel size, so recompute it only
            // when the size changes -- a plain scroll (SIMPLE_SCROLL_MODE repaints the viewport) reuses the coords
            if ( ( w != _laid_out_w ) || ( h != _laid_out_h ) ) {
                layoutTree( _left, _left_tips, tree_w, avail_h );
                layoutTree( _right, _right_tips, tree_w, avail_h );
                _laid_out_w = w;
                _laid_out_h = h;
            }
            final double left_root_x = MARGIN;
            final double left_tips_x = MARGIN + tree_w;
            final double left_label_end = left_tips_x + w_left;
            final double right_root_x = w - MARGIN;
            final double right_tips_x = w - MARGIN - tree_w;
            final double right_label_start = right_tips_x - w_right;
            final Color branch_color = getForeground();
            final Color connector_color = withAlpha( branch_color, 120 );
            final Color unmatched_color = blend( branch_color, getBackground(), 0.55 );
            drawTree( g2, _left, left_root_x, false, branch_color );
            drawTree( g2, _right, right_root_x, true, branch_color );
            drawTipLabels( g2, _left_tips, left_tips_x + LABEL_GAP, false, w_left, fm, branch_color, unmatched_color );
            drawTipLabels( g2, _right_tips, right_tips_x - LABEL_GAP, true, w_right, fm, branch_color, unmatched_color );
            drawConnectors( g2, left_label_end, right_label_start, connector_color );
        }
        finally {
            g2.dispose();
        }
    }

    // ---- layout ------------------------------------------------------------------------------------------------

    /** Assigns each node an (x,y): x is a lined-up cladogram depth in [0, tree_w] (tips flush at tree_w); y is the
     *  screen y filling the available height, tips evenly spaced and internal nodes at their child block's midpoint. */
    private static void layoutTree( final Phylogeny phy, final List<PhylogenyNode> tips, final double tree_w,
                                    final double avail_h ) {
        final int n = tips.size();
        if ( n == 0 ) {
            return;
        }
        final double half_pitch = avail_h / ( 2.0 * n );
        for( int i = 0; i < n; i++ ) {
            tips.get( i ).setYcoord( (float) ( TOP_PAD + ( ( ( 2 * i ) + 1 ) * half_pitch ) ) );
        }
        final Map<PhylogenyNode, Integer> heights = new IdentityHashMap<>();
        final int tree_height = computeHeight( phy.getRoot(), heights );
        final double x_step = ( tree_height > 0 ) ? ( tree_w / tree_height ) : tree_w;
        for( final PhylogenyNode node : heights.keySet() ) {
            node.setXcoord( (float) ( ( tree_height - heights.get( node ) ) * x_step ) );
        }
        assignInternalY( phy.getRoot() );
    }

    private static int computeHeight( final PhylogenyNode node, final Map<PhylogenyNode, Integer> heights ) {
        int height = 0;
        if ( !node.isExternal() ) {
            for( final PhylogenyNode child : node.getDescendants() ) {
                height = Math.max( height, computeHeight( child, heights ) + 1 );
            }
        }
        heights.put( node, height );
        return height;
    }

    private static double assignInternalY( final PhylogenyNode node ) {
        if ( node.isExternal() ) {
            return node.getYcoord();
        }
        final List<PhylogenyNode> children = node.getDescendants();
        double first_y = 0;
        double last_y = 0;
        for( int i = 0; i < children.size(); i++ ) {
            final double child_y = assignInternalY( children.get( i ) );
            if ( i == 0 ) {
                first_y = child_y;
            }
            if ( i == ( children.size() - 1 ) ) {
                last_y = child_y;
            }
        }
        final double y = ( first_y + last_y ) / 2.0;
        node.setYcoord( (float) y );
        return y;
    }

    // ---- drawing -----------------------------------------------------------------------------------------------

    private static void drawTree( final Graphics2D g2, final Phylogeny phy, final double root_screen_x,
                                  final boolean mirror, final Color color ) {
        g2.setColor( color );
        g2.setStroke( new BasicStroke( 1f ) );
        for( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            final double sx = screenX( node, root_screen_x, mirror );
            final double sy = node.getYcoord();
            if ( node.getParent() != null ) {
                final double px = screenX( node.getParent(), root_screen_x, mirror );
                g2.drawLine( round( px ), round( sy ), round( sx ), round( sy ) );
            }
            if ( !node.isExternal() ) {
                final List<PhylogenyNode> children = node.getDescendants();
                final double y_first = children.get( 0 ).getYcoord();
                final double y_last = children.get( children.size() - 1 ).getYcoord();
                g2.drawLine( round( sx ), round( y_first ), round( sx ), round( y_last ) );
            }
        }
    }

    private static double screenX( final PhylogenyNode node, final double root_screen_x, final boolean mirror ) {
        return mirror ? ( root_screen_x - node.getXcoord() ) : ( root_screen_x + node.getXcoord() );
    }

    private void drawTipLabels( final Graphics2D g2, final List<PhylogenyNode> tips, final double edge_x,
                                final boolean mirror, final double max_width, final FontMetrics fm,
                                final Color matched_color, final Color unmatched_color ) {
        final int ascent_offset = ( fm.getAscent() - fm.getDescent() ) / 2;
        for( final PhylogenyNode tip : tips ) {
            final String label = labelFor( tip );
            if ( label.isEmpty() ) {
                continue;
            }
            final String shown = truncateToWidth( label, fm, (int) max_width );
            g2.setColor( _unmatched.contains( tip ) ? unmatched_color : matched_color );
            final int y = round( tip.getYcoord() ) + ascent_offset;
            final int x = mirror ? ( round( edge_x ) - fm.stringWidth( shown ) ) : round( edge_x );
            g2.drawString( shown, x, y );
        }
    }

    private void drawConnectors( final Graphics2D g2, final double left_x, final double right_x, final Color color ) {
        g2.setColor( color );
        final Stroke old = g2.getStroke();
        g2.setStroke( new BasicStroke( 1.2f ) );
        for( final TanglegramLinker.Link link : _result.getLinks() ) {
            g2.drawLine( round( left_x ), round( link.getA().getYcoord() ), round( right_x ),
                         round( link.getB().getYcoord() ) );
        }
        g2.setStroke( old );
    }

    private static double maxLabelWidth( final List<PhylogenyNode> tips, final FontMetrics fm ) {
        int max = 0;
        for( final PhylogenyNode tip : tips ) {
            max = Math.max( max, fm.stringWidth( labelFor( tip ) ) );
        }
        return max + ( 2.0 * LABEL_GAP );
    }

    /** The displayed tip label: the first non-empty of node name, scientific name, taxonomy code, common name, NCBI
     *  id, sequence name -- so a tip linked on ANY of those fields still shows a readable label (not a blank row). */
    private static String labelFor( final PhylogenyNode tip ) {
        if ( notEmpty( tip.getName() ) ) {
            return tip.getName().trim();
        }
        if ( tip.getNodeData().isHasTaxonomy() ) {
            if ( notEmpty( tip.getNodeData().getTaxonomy().getScientificName() ) ) {
                return tip.getNodeData().getTaxonomy().getScientificName().trim();
            }
            if ( notEmpty( tip.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
                return tip.getNodeData().getTaxonomy().getTaxonomyCode().trim();
            }
            if ( notEmpty( tip.getNodeData().getTaxonomy().getCommonName() ) ) {
                return tip.getNodeData().getTaxonomy().getCommonName().trim();
            }
        }
        if ( notEmpty( PhylogenyMethods.getTaxonomyIdentifier( tip ) ) ) {
            return PhylogenyMethods.getTaxonomyIdentifier( tip ).trim();
        }
        if ( tip.getNodeData().isHasSequence() && notEmpty( tip.getNodeData().getSequence().getName() ) ) {
            return tip.getNodeData().getSequence().getName().trim();
        }
        return "";
    }

    private static boolean notEmpty( final String s ) {
        return ( s != null ) && !s.trim().isEmpty();
    }

    private static String truncateToWidth( final String s, final FontMetrics fm, final int max_width ) {
        if ( max_width <= 0 ) {
            // no column room at all (an extremely narrow window) -- hide the label rather than overrun the trees
            return "";
        }
        if ( fm.stringWidth( s ) <= max_width ) {
            return s;
        }
        final String ellipsis = "…";
        final int ell_w = fm.stringWidth( ellipsis );
        int end = s.length();
        while ( ( end > 0 ) && ( ( fm.stringWidth( s.substring( 0, end ) ) + ell_w ) > max_width ) ) {
            end--;
        }
        return ( end <= 0 ) ? ellipsis : ( s.substring( 0, end ) + ellipsis );
    }

    private static int round( final double v ) {
        return (int) Math.round( v );
    }

    private static Color withAlpha( final Color c, final int alpha ) {
        return new Color( c.getRed(), c.getGreen(), c.getBlue(), alpha );
    }

    private static Color blend( final Color a, final Color b, final double fraction ) {
        final double f = Math.max( 0.0, Math.min( 1.0, fraction ) );
        final int r = (int) Math.round( ( a.getRed() * ( 1 - f ) ) + ( b.getRed() * f ) );
        final int g = (int) Math.round( ( a.getGreen() * ( 1 - f ) ) + ( b.getGreen() * f ) );
        final int bl = (int) Math.round( ( a.getBlue() * ( 1 - f ) ) + ( b.getBlue() * f ) );
        return new Color( r, g, bl );
    }

    // ---- Scrollable: fill the viewport width, and the height too until the tree needs to scroll ------------------

    // Zoom scales the VERTICAL axis only (tip-row spacing) -- that is the crowded axis of a tanglegram; the width
    // always tracks the viewport (see getScrollableTracksViewportWidth) so both trees stay spread across the window.
    @Override
    public Dimension getPreferredSize() {
        return new Dimension( BASE_WIDTH, (int) Math.round( naturalHeight() * _zoom ) );
    }

    @Override
    public Dimension getPreferredScrollableViewportSize() {
        return getPreferredSize();
    }

    @Override
    public int getScrollableUnitIncrement( final Rectangle visible, final int orientation, final int direction ) {
        return 16;
    }

    @Override
    public int getScrollableBlockIncrement( final Rectangle visible, final int orientation, final int direction ) {
        return ( orientation == javax.swing.SwingConstants.VERTICAL ) ? visible.height : visible.width;
    }

    @Override
    public boolean getScrollableTracksViewportWidth() {
        return true;
    }

    @Override
    public boolean getScrollableTracksViewportHeight() {
        return ( getParent() instanceof JViewport ) && ( getPreferredSize().height <= getParent().getHeight() );
    }
}
