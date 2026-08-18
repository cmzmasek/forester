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
import java.awt.Font;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.FontType;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;

/**
 * Headless tests for {@link NodeStyleEditor}: the per-node visual-style applier writes ONLY the spec's ticked
 * attributes (so a bulk edit never clobbers the others), creates a {@link NodeVisualData} where absent, fills in a
 * font family + size when a per-node font is set (so it actually renders), and produces a sensible provenance
 * sentence.
 */
public final class NodeStyleEditorTest {

    public static boolean test() {
        try {
            return perAttributeOk() && bulkOk() && leavesOthersOk() && fontFallbackOk() && keepsExistingSizeOk()
                    && emptySpecOk() && provenanceOk();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** A single-node edit writes exactly the ticked attributes; the rest stay at their defaults. */
    private static boolean perAttributeOk() {
        final PhylogenyNode node = new PhylogenyNode();
        final NodeStyleEditor.Spec spec = new NodeStyleEditor.Spec( null, null, Color.RED, NodeShape.CIRCLE, null,
                                                                    null, null );
        final int n = NodeStyleEditor.apply( Arrays.asList( node ), spec, "Source Sans 3", 12 );
        ck( n == 1, "one node modified" );
        final NodeVisualData vis = node.getNodeData().getNodeVisualData();
        ck( vis != null, "visual data created" );
        ck( ( vis.getFontColor() != null ) && vis.getFontColor().equals( Color.RED ), "font color written" );
        ck( vis.getShape() == NodeShape.CIRCLE, "shape written" );
        // un-ticked attributes stay default
        ck( vis.getFillType() == NodeFill.DEFAULT, "fill left default" );
        ck( vis.getNodeColor() == null, "node color left unset" );
        ck( vis.getFontStyle() == FontType.PLAIN, "font style left default" );
        return true;
    }

    /** The same spec applies to every node in a bulk target list. */
    private static boolean bulkOk() {
        final List<PhylogenyNode> nodes = new ArrayList<>();
        for ( int i = 0; i < 3; ++i ) {
            nodes.add( new PhylogenyNode() );
        }
        final NodeStyleEditor.Spec spec = new NodeStyleEditor.Spec( null, null, null, null, NodeFill.SOLID, null,
                                                                    Color.BLUE );
        final int n = NodeStyleEditor.apply( nodes, spec, "Source Sans 3", 12 );
        ck( n == 3, "three nodes modified" );
        for ( final PhylogenyNode node : nodes ) {
            final NodeVisualData vis = node.getNodeData().getNodeVisualData();
            ck( vis.getFillType() == NodeFill.SOLID, "fill written on each node" );
            ck( ( vis.getNodeColor() != null ) && vis.getNodeColor().equals( Color.BLUE ), "node color on each node" );
        }
        return true;
    }

    /** A per-attribute edit does NOT overwrite a node's OTHER, already-set attributes (the anti-clobber guarantee). */
    private static boolean leavesOthersOk() {
        final PhylogenyNode node = new PhylogenyNode();
        final NodeVisualData pre = new NodeVisualData();
        pre.setShape( NodeShape.RECTANGLE );
        pre.setSize( 9f );
        node.getNodeData().setNodeVisualData( pre );
        // apply only a font color
        NodeStyleEditor.apply( Arrays.asList( node ),
                               new NodeStyleEditor.Spec( null, null, Color.GREEN, null, null, null, null ),
                               "Source Sans 3", 12 );
        final NodeVisualData vis = node.getNodeData().getNodeVisualData();
        ck( ( vis.getFontColor() != null ) && vis.getFontColor().equals( Color.GREEN ), "font color added" );
        ck( vis.getShape() == NodeShape.RECTANGLE, "pre-existing shape NOT clobbered" );
        ck( vis.getSize() == 9f, "pre-existing node size NOT clobbered" );
        return true;
    }

    /** Setting a font STYLE with no family/size fills both from the tree defaults, so getFont() renders a real font. */
    private static boolean fontFallbackOk() {
        final PhylogenyNode node = new PhylogenyNode();
        NodeStyleEditor.apply( Arrays.asList( node ),
                               new NodeStyleEditor.Spec( Integer.valueOf( Font.BOLD ), null, null, null, null, null,
                                                         null ),
                               "Source Sans 3", 14 );
        final Font f = node.getNodeData().getNodeVisualData().getFont();
        ck( f != null, "getFont() builds a real font (family filled in)" );
        ck( f.getStyle() == Font.BOLD, "font style is bold" );
        ck( f.getSize() == 14, "font size filled from the tree default (14)" );
        return true;
    }

    /** When only the style changes, an EXISTING per-node font size is kept, not overwritten with the tree default. */
    private static boolean keepsExistingSizeOk() {
        final PhylogenyNode node = new PhylogenyNode();
        final NodeVisualData pre = new NodeVisualData();
        pre.setFontName( "Source Sans 3" );
        pre.setFontSize( 20 );
        node.getNodeData().setNodeVisualData( pre );
        NodeStyleEditor.apply( Arrays.asList( node ),
                               new NodeStyleEditor.Spec( Integer.valueOf( Font.ITALIC ), null, null, null, null, null,
                                                         null ),
                               "Source Sans 3", 12 );
        final Font f = node.getNodeData().getNodeVisualData().getFont();
        ck( f.getSize() == 20, "existing font size (20) preserved, not reset to the default 12" );
        ck( f.getStyle() == Font.ITALIC, "font style updated to italic" );
        return true;
    }

    /** An empty spec is a no-op and creates no NodeVisualData. */
    private static boolean emptySpecOk() {
        final PhylogenyNode node = new PhylogenyNode();
        final int n = NodeStyleEditor.apply( Arrays.asList( node ),
                                             new NodeStyleEditor.Spec( null, null, null, null, null, null, null ),
                                             "Source Sans 3", 12 );
        ck( n == 0, "empty spec modifies nothing" );
        ck( node.getNodeData().getNodeVisualData() == null, "empty spec creates no visual data" );
        return true;
    }

    private static boolean provenanceOk() {
        final NodeStyleEditor.Spec spec = new NodeStyleEditor.Spec( null, null, Color.RED, NodeShape.CIRCLE, null,
                                                                    null, null );
        final String p = NodeStyleEditor.provenance( spec, 3 );
        ck( p.contains( "font color" ) && p.contains( "node shape" ), "provenance names the changed attributes: " + p );
        ck( p.contains( "3 nodes" ), "provenance has the plural count: " + p );
        ck( !p.contains( "node size" ) && !p.contains( "node fill" ), "provenance omits un-changed attributes: " + p );
        final String one = NodeStyleEditor.provenance(
                new NodeStyleEditor.Spec( null, null, Color.RED, null, null, null, null ), 1 );
        ck( one.contains( "1 node." ), "singular count: " + one );
        return true;
    }

    private static void ck( final boolean cond, final String msg ) {
        if ( !cond ) {
            throw new RuntimeException( "NodeStyleEditorTest: " + msg );
        }
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private NodeStyleEditorTest() {
    }
}
