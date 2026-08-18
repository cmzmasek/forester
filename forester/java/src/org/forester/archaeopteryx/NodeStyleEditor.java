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
import java.util.List;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.util.ForesterUtil;

/**
 * Pure, Swing-free applier for per-node visual style ({@link NodeVisualData}: font style/size/colour + node
 * shape/fill/size/colour). A {@link Spec} carries ONLY the attributes the user chose to change -- each null field
 * means "leave this attribute unchanged" -- so a bulk edit over many nodes writes just the ticked attributes and
 * never clobbers the rest. Reached from the "Node Style" click-to action (one node) and the Tools "Node Style for
 * Selected Nodes" action (the found/selected set). Rendering + phyloXML round-trip already exist; this is the editor
 * side. Headless-testable.
 */
final class NodeStyleEditor {

    /** The attributes to write. A null field = "leave that attribute unchanged". Font STYLE and SIZE need a font
     *  FAMILY for the renderer's {@link NodeVisualData#getFont()} to build a Font (it returns null without a name),
     *  so {@link #apply} supplies a default family when the node has none. */
    static final class Spec {

        private final Integer   _font_style; // java.awt.Font PLAIN / BOLD / ITALIC / BOLD+ITALIC, or null
        private final Integer   _font_size;  // points, or null
        private final Color     _font_color;
        private final NodeShape _shape;
        private final NodeFill  _fill;
        private final Float     _node_size;
        private final Color     _node_color;

        Spec( final Integer font_style,
              final Integer font_size,
              final Color font_color,
              final NodeShape shape,
              final NodeFill fill,
              final Float node_size,
              final Color node_color ) {
            _font_style = font_style;
            _font_size = font_size;
            _font_color = font_color;
            _shape = shape;
            _fill = fill;
            _node_size = node_size;
            _node_color = node_color;
        }

        boolean isEmpty() {
            return ( _font_style == null ) && ( _font_size == null ) && ( _font_color == null ) && ( _shape == null )
                    && ( _fill == null ) && ( _node_size == null ) && ( _node_color == null );
        }

        /** Whether the spec sets the per-node FONT (style or size) -- which additionally requires a font family. */
        boolean touchesFont() {
            return ( _font_style != null ) || ( _font_size != null );
        }

        Integer fontStyle() {
            return _font_style;
        }

        Integer fontSize() {
            return _font_size;
        }

        Color fontColor() {
            return _font_color;
        }

        NodeShape shape() {
            return _shape;
        }

        NodeFill fill() {
            return _fill;
        }

        Float nodeSize() {
            return _node_size;
        }

        Color nodeColor() {
            return _node_color;
        }
    }

    /**
     * Applies {@code spec} to every node in {@code targets}, creating a {@link NodeVisualData} where absent and
     * writing ONLY the spec's non-null attributes. Setting a per-node font style/size needs a family AND a size for
     * the renderer's {@link NodeVisualData#getFont()} to build a valid Font ({@code new Font(name, style, size)}
     * returns nothing without a name and renders at size 0 without a size), so where the node has none the tree's
     * current {@code default_font_family} / {@code default_font_size} fill in (invisible to the user, per the "core
     * set" choice that hides the font family). Returns the number of nodes actually modified.
     */
    static int apply( final List<PhylogenyNode> targets, final Spec spec, final String default_font_family,
                      final int default_font_size ) {
        if ( ( targets == null ) || ( spec == null ) || spec.isEmpty() ) {
            return 0;
        }
        int n = 0;
        for ( final PhylogenyNode node : targets ) {
            if ( node == null ) {
                continue;
            }
            NodeVisualData vis = node.getNodeData().getNodeVisualData();
            if ( vis == null ) {
                vis = new NodeVisualData();
                node.getNodeData().setNodeVisualData( vis );
            }
            if ( spec.touchesFont() && ForesterUtil.isEmpty( vis.getFontName() ) ) {
                vis.setFontName( default_font_family );
            }
            if ( spec.fontStyle() != null ) {
                vis.setFontStyle( spec.fontStyle().intValue() );
            }
            if ( spec.fontSize() != null ) {
                vis.setFontSize( spec.fontSize().intValue() );
            }
            if ( spec.touchesFont() && ( vis.getFontSize() <= 0 ) ) {
                // a per-node font needs a real size, else getFont() renders at size 0; keep an existing size, else
                // fall back to the tree's current size
                vis.setFontSize( default_font_size );
            }
            if ( spec.fontColor() != null ) {
                vis.setFontColor( spec.fontColor() );
            }
            if ( spec.shape() != null ) {
                vis.setShape( spec.shape() );
            }
            if ( spec.fill() != null ) {
                vis.setFillType( spec.fill() );
            }
            if ( spec.nodeSize() != null ) {
                vis.setSize( spec.nodeSize().floatValue() );
            }
            if ( spec.nodeColor() != null ) {
                vis.setNodeColor( spec.nodeColor() );
            }
            n++;
        }
        return n;
    }

    /** A human-readable provenance sentence naming the attributes changed, for the tree description. */
    static String provenance( final Spec spec, final int count ) {
        final StringBuilder attrs = new StringBuilder();
        appendAttr( attrs, spec.fontStyle() != null, "font style" );
        appendAttr( attrs, spec.fontSize() != null, "font size" );
        appendAttr( attrs, spec.fontColor() != null, "font color" );
        appendAttr( attrs, spec.shape() != null, "node shape" );
        appendAttr( attrs, spec.fill() != null, "node fill" );
        appendAttr( attrs, spec.nodeSize() != null, "node size" );
        appendAttr( attrs, spec.nodeColor() != null, "node color" );
        return "Set the visual style (" + attrs + ") of " + count + ( count == 1 ? " node." : " nodes." );
    }

    private static void appendAttr( final StringBuilder sb, final boolean present, final String label ) {
        if ( present ) {
            if ( sb.length() > 0 ) {
                sb.append( ", " );
            }
            sb.append( label );
        }
    }

    private NodeStyleEditor() {
    }
}
