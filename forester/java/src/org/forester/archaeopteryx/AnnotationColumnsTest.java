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
import java.util.ArrayList;
import java.util.List;

import org.forester.archaeopteryx.AnnotationColumns.Type;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;

/**
 * Unit tests for {@link AnnotationColumns} (the pure per-cell model behind the tip-aligned annotation
 * columns). Same package as the class under test (package-private). Headless; run via the suite or
 * {@link #main(String[])}.
 */
public final class AnnotationColumnsTest {

    public static void main( final String[] args ) {
        System.out.println( "AnnotationColumns: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        return testTypeSuggestion() && testCells() && testHeaders() && testLabels() && testMatrix() && testSymbol();
    }

    // ---- SYMBOL column: a shape glyph colored like a strip, with presence-driven fill ----
    private static boolean testSymbol() {
        // isFalsy: the explicit false/absent tokens (case-insensitive, trimmed) render as a hollow OUTLINE
        if ( !AnnotationColumns.isFalsy( "no" ) || !AnnotationColumns.isFalsy( "No" )
                || !AnnotationColumns.isFalsy( " 0 " ) || !AnnotationColumns.isFalsy( "FALSE" )
                || !AnnotationColumns.isFalsy( "absent" ) || !AnnotationColumns.isFalsy( "-" )
                || !AnnotationColumns.isFalsy( "n" ) ) {
            return fail( "isFalsy should accept the false/absent tokens (case-insensitive, trimmed)" );
        }
        if ( AnnotationColumns.isFalsy( "yes" ) || AnnotationColumns.isFalsy( "human" )
                || AnnotationColumns.isFalsy( "1" ) || AnnotationColumns.isFalsy( "positive" )
                || AnnotationColumns.isFalsy( "" ) || AnnotationColumns.isFalsy( null ) ) {
            return fail( "isFalsy should reject present/other values (blank/null are handled as NONE, not falsy)" );
        }
        // a SYMBOL is offered for a categorical field (with COLOR_STRIP + TEXT) but NOT for a numeric one
        final Phylogeny phy = tree();
        if ( !AnnotationColumns.allowedTypes( phy, "repseq:host" ).contains( Type.SYMBOL ) ) {
            return fail( "a categorical field should allow SYMBOL" );
        }
        if ( AnnotationColumns.allowedTypes( phy, "data:score" ).contains( Type.SYMBOL ) ) {
            return fail( "a numeric field should NOT allow SYMBOL" );
        }
        // symbolFill: present -> FILLED, explicit-false -> OUTLINE, missing -> NONE
        final Phylogeny bphy = binaryTree();
        final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<AnnotationColumns.ColumnSpec>();
        specs.add( new AnnotationColumns.ColumnSpec( "data:present", Type.SYMBOL ) );      // 0
        specs.add( new AnnotationColumns.ColumnSpec( "data:present", Type.COLOR_STRIP ) ); // 1 (non-SYMBOL)
        final AnnotationColumns ac = new AnnotationColumns( bphy, specs );
        if ( ac.symbolFill( tip( bphy, "yes_tip" ), 0 ) != AnnotationColumns.Fill.FILLED ) {
            return fail( "a present value should be a FILLED symbol" );
        }
        if ( ac.symbolFill( tip( bphy, "no_tip" ), 0 ) != AnnotationColumns.Fill.OUTLINE ) {
            return fail( "an explicitly-false value should be a hollow OUTLINE symbol" );
        }
        if ( ac.symbolFill( tip( bphy, "blank_tip" ), 0 ) != AnnotationColumns.Fill.NONE ) {
            return fail( "a missing value should draw NONE" );
        }
        // a SYMBOL column carries a (categorical) cell color, so it colors like a strip and gets a legend
        if ( ac.cellColor( tip( bphy, "yes_tip" ), 0 ) == null ) {
            return fail( "a SYMBOL column should carry a categorical cell color" );
        }
        // symbolFill is always NONE for a non-SYMBOL column
        if ( ac.symbolFill( tip( bphy, "yes_tip" ), 1 ) != AnnotationColumns.Fill.NONE ) {
            return fail( "symbolFill on a non-SYMBOL column must be NONE" );
        }
        return true;
    }

    private static Phylogeny binaryTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( binaryLeaf( "yes_tip", "yes" ) );
        root.addAsChild( binaryLeaf( "no_tip", "no" ) );
        root.addAsChild( binaryLeaf( "blank_tip", null ) ); // no property -> a missing value
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode binaryLeaf( final String name, final String present ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        if ( present != null ) {
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:present", present, "", "xsd:string", AppliesTo.NODE ) );
            n.getNodeData().setProperties( pl );
        }
        return n;
    }

    // ---- heat-map MATRIX: all matrix columns share ONE color scale ----
    private static boolean testMatrix() {
        // two numeric columns with DIFFERENT per-column ranges: s1 in [0,10], s2 in [0,100].
        final Phylogeny phy = matrixTree();
        final List<AnnotationColumns.ColumnSpec> ms = new ArrayList<AnnotationColumns.ColumnSpec>();
        ms.add( new AnnotationColumns.ColumnSpec( "data:s1", Type.MATRIX ) ); // 0
        ms.add( new AnnotationColumns.ColumnSpec( "data:s2", Type.MATRIX ) ); // 1
        final AnnotationColumns matrix = new AnnotationColumns( phy, ms );
        // SHARED scale: the value 10 in s1 and the value 10 in s2 map to the SAME color (both 10 on the global 0..100)
        final Color m_s1_10 = matrix.cellColor( tip( phy, "p" ), 0 );
        final Color m_s2_10 = matrix.cellColor( tip( phy, "p" ), 1 );
        if ( ( m_s1_10 == null ) || !m_s1_10.equals( m_s2_10 ) ) {
            return fail( "MATRIX: equal values in different columns must share a color (shared scale): " + m_s1_10
                    + " vs " + m_s2_10 );
        }
        // PER-COLUMN HEATMAP on the SAME data: 10 is the MAX of s1 (top of its scale) but low in s2 -> DIFFERENT colors
        final List<AnnotationColumns.ColumnSpec> hs = new ArrayList<AnnotationColumns.ColumnSpec>();
        hs.add( new AnnotationColumns.ColumnSpec( "data:s1", Type.HEATMAP ) );
        hs.add( new AnnotationColumns.ColumnSpec( "data:s2", Type.HEATMAP ) );
        final AnnotationColumns heat = new AnnotationColumns( phy, hs );
        final Color h_s1_10 = heat.cellColor( tip( phy, "p" ), 0 );
        final Color h_s2_10 = heat.cellColor( tip( phy, "p" ), 1 );
        if ( ( h_s1_10 == null ) || h_s1_10.equals( h_s2_10 ) ) {
            return fail( "per-column HEATMAP: the value 10 should color differently in s1 (its max) vs s2 (low): "
                    + h_s1_10 + " vs " + h_s2_10 );
        }
        // a MATRIX column whose values are ALL identical would, ALONE, auto-detect as CATEGORICAL (one distinct
        // value) -- it must STILL color on the shared gradient, not a categorical palette color, or that one column
        // silently breaks the shared scale. Here s2 is constant (50) while s1 spans 0..100 (shared range [0,100]);
        // the constant column's cell for 50 must match the SAME value's cell in the spanning column.
        final Phylogeny phy2 = new Phylogeny();
        final PhylogenyNode root2 = new PhylogenyNode();
        root2.addAsChild( matrixLeaf( "a", "0", "50" ) );
        root2.addAsChild( matrixLeaf( "b", "100", "50" ) );
        root2.addAsChild( matrixLeaf( "c", "50", "50" ) );
        phy2.setRoot( root2 );
        phy2.externalNodesHaveChanged();
        final List<AnnotationColumns.ColumnSpec> ds = new ArrayList<AnnotationColumns.ColumnSpec>();
        ds.add( new AnnotationColumns.ColumnSpec( "data:s1", Type.MATRIX ) ); // spans 0..100
        ds.add( new AnnotationColumns.ColumnSpec( "data:s2", Type.MATRIX ) ); // constant 50 -> categorical alone
        final AnnotationColumns degen = new AnnotationColumns( phy2, ds );
        final Color d_spanning_50 = degen.cellColor( tip( phy2, "c" ), 0 ); // 50 on the shared [0,100] gradient
        final Color d_constant_50 = degen.cellColor( tip( phy2, "c" ), 1 ); // 50 in the all-equal column
        if ( ( d_constant_50 == null ) || !d_constant_50.equals( d_spanning_50 ) ) {
            return fail( "MATRIX: an all-equal column must still use the shared gradient (value 50 -> same color as "
                    + "the spanning column), not a categorical palette color: " + d_constant_50 + " vs "
                    + d_spanning_50 );
        }
        return true;
    }

    private static Phylogeny matrixTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( matrixLeaf( "p", "10", "10" ) );
        root.addAsChild( matrixLeaf( "q", "10", "100" ) );
        root.addAsChild( matrixLeaf( "r", "0", "0" ) );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode matrixLeaf( final String name, final String s1, final String s2 ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( "data:s1", s1, "", "xsd:string", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:s2", s2, "", "xsd:string", AppliesTo.NODE ) );
        n.getNodeData().setProperties( pl );
        return n;
    }

    // ---- friendly render-type labels ----
    private static boolean testLabels() {
        if ( !"Color strip".equals( AnnotationColumns.label( Type.COLOR_STRIP ) )
                || !"Symbol".equals( AnnotationColumns.label( Type.SYMBOL ) )
                || !"Heat map".equals( AnnotationColumns.label( Type.HEATMAP ) )
                || !"Bar".equals( AnnotationColumns.label( Type.BAR ) )
                || !"Text".equals( AnnotationColumns.label( Type.TEXT ) ) ) {
            return fail( "render-type labels are wrong" );
        }
        return true;
    }

    // ---- default/allowed render types follow the data (categorical vs numeric) ----
    private static boolean testTypeSuggestion() {
        final Phylogeny phy = tree();
        if ( AnnotationColumns.defaultType( phy, "repseq:host" ) != Type.COLOR_STRIP ) {
            return fail( "a categorical field should default to COLOR_STRIP" );
        }
        if ( AnnotationColumns.defaultType( phy, "data:score" ) != Type.HEATMAP ) {
            return fail( "a numeric field should default to HEATMAP" );
        }
        final List<Type> host_types = AnnotationColumns.allowedTypes( phy, "repseq:host" );
        if ( ( host_types.size() != 3 ) || !host_types.contains( Type.COLOR_STRIP )
                || !host_types.contains( Type.SYMBOL ) || !host_types.contains( Type.TEXT )
                || host_types.contains( Type.HEATMAP ) ) {
            return fail( "a categorical field should allow COLOR_STRIP + SYMBOL + TEXT" );
        }
        final List<Type> score_types = AnnotationColumns.allowedTypes( phy, "data:score" );
        if ( ( score_types.size() != 4 ) || !score_types.contains( Type.HEATMAP )
                || !score_types.contains( Type.MATRIX ) || !score_types.contains( Type.BAR )
                || !score_types.contains( Type.TEXT ) || score_types.contains( Type.COLOR_STRIP )
                || score_types.contains( Type.SYMBOL ) ) {
            return fail( "a numeric field should allow HEATMAP + MATRIX + BAR + TEXT" );
        }
        return true;
    }

    // ---- per-cell color / bar fraction / text ----
    private static boolean testCells() {
        final Phylogeny phy = tree();
        final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<AnnotationColumns.ColumnSpec>();
        specs.add( new AnnotationColumns.ColumnSpec( "repseq:host", Type.COLOR_STRIP ) ); // 0
        specs.add( new AnnotationColumns.ColumnSpec( "data:score", Type.HEATMAP ) );       // 1
        specs.add( new AnnotationColumns.ColumnSpec( "data:score", Type.BAR ) );           // 2
        specs.add( new AnnotationColumns.ColumnSpec( "repseq:host", Type.TEXT ) );          // 3
        final AnnotationColumns ac = new AnnotationColumns( phy, specs );
        if ( ac.size() != 4 ) {
            return fail( "expected 4 columns" );
        }
        // color-strip: same value -> same color; different value -> different; all non-null
        final Color a = ac.cellColor( tip( phy, "A" ), 0 ); // cat
        final Color c = ac.cellColor( tip( phy, "C" ), 0 ); // cat
        final Color b = ac.cellColor( tip( phy, "B" ), 0 ); // dog
        if ( ( a == null ) || ( b == null ) || !a.equals( c ) || a.equals( b ) ) {
            return fail( "color-strip: cat==cat, cat!=dog, none null" );
        }
        // heat-map: numeric -> non-null, min and max differ
        final Color h_min = ac.cellColor( tip( phy, "A" ), 1 ); // score 10
        final Color h_max = ac.cellColor( tip( phy, "D" ), 1 ); // score 40
        if ( ( h_min == null ) || ( h_max == null ) || h_min.equals( h_max ) ) {
            return fail( "heat-map: min/max colors should differ and be non-null" );
        }
        // bar / text columns carry no cell color
        if ( ( ac.cellColor( tip( phy, "A" ), 2 ) != null ) || ( ac.cellColor( tip( phy, "A" ), 3 ) != null ) ) {
            return fail( "bar/text columns should have no cell color" );
        }
        // bar fraction spans 0..1 across the numeric range; NaN for non-bar columns
        if ( Math.abs( ac.barFraction( tip( phy, "A" ), 2 ) - 0.0 ) > 1e-9 ) {
            return fail( "min value should be bar fraction 0" );
        }
        if ( Math.abs( ac.barFraction( tip( phy, "D" ), 2 ) - 1.0 ) > 1e-9 ) {
            return fail( "max value should be bar fraction 1" );
        }
        if ( Math.abs( ac.barFraction( tip( phy, "B" ), 2 ) - ( 1.0 / 3.0 ) ) > 1e-9 ) {
            return fail( "score 20 should be bar fraction 1/3" );
        }
        if ( !Double.isNaN( ac.barFraction( tip( phy, "A" ), 1 ) )
                || !Double.isNaN( ac.barFraction( tip( phy, "A" ), 0 ) ) ) {
            return fail( "non-bar columns should have NaN bar fraction" );
        }
        // text column: the raw value
        if ( !"cat".equals( ac.cellText( tip( phy, "A" ), 3 ) ) || !"dog".equals( ac.cellText( tip( phy, "B" ), 3 ) ) ) {
            return fail( "text column should show the raw value" );
        }
        return true;
    }

    // ---- headers are the prettified field names ----
    private static boolean testHeaders() {
        final Phylogeny phy = tree();
        final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<AnnotationColumns.ColumnSpec>();
        specs.add( new AnnotationColumns.ColumnSpec( "repseq:host", Type.COLOR_STRIP ) );
        specs.add( new AnnotationColumns.ColumnSpec( "data:score", Type.HEATMAP ) );
        final AnnotationColumns ac = new AnnotationColumns( phy, specs );
        if ( !"Host".equals( ac.getColumn( 0 ).getHeader() ) ) {
            return fail( "header should be 'Host', got " + ac.getColumn( 0 ).getHeader() );
        }
        if ( !"Score".equals( ac.getColumn( 1 ).getHeader() ) ) {
            return fail( "header should be 'Score', got " + ac.getColumn( 1 ).getHeader() );
        }
        return true;
    }

    // ---------------------------------------------------------------------------------------

    private static Phylogeny tree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( leaf( "A", "cat", "10" ) );
        root.addAsChild( leaf( "B", "dog", "20" ) );
        root.addAsChild( leaf( "C", "cat", "30" ) );
        root.addAsChild( leaf( "D", "fish", "40" ) );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode leaf( final String name, final String host, final String score ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( "repseq:host", host, "", "xsd:string", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:score", score, "", "xsd:string", AppliesTo.NODE ) );
        n.getNodeData().setProperties( pl );
        return n;
    }

    private static PhylogenyNode tip( final Phylogeny phy, final String name ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        return null;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [AnnotationColumnsTest] " + message );
        return false;
    }

    private AnnotationColumnsTest() {
        // not instantiable
    }
}
