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
        return testTypeSuggestion() && testCells() && testHeaders() && testLabels() && testMatrix() && testSymbol()
                && testStackedBar() && testPie();
    }

    // ---- PIE column: several numeric fields MERGE into one pie glyph, wedge angle = the tip's own proportion ----
    private static boolean testPie() {
        final Phylogeny phy = stackTree();
        if ( !"Pie".equals( AnnotationColumns.label( Type.PIE ) ) ) {
            return fail( "PIE label should be 'Pie'" );
        }
        if ( !AnnotationColumns.isMergedType( Type.PIE ) || !AnnotationColumns.isMergedType( Type.STACKED_BAR )
                || AnnotationColumns.isMergedType( Type.BAR ) || AnnotationColumns.isMergedType( Type.COLOR_STRIP ) ) {
            return fail( "isMergedType should be true only for STACKED_BAR and PIE" );
        }
        if ( !AnnotationColumns.allowedTypes( phy, "data:a" ).contains( Type.PIE ) ) {
            return fail( "a numeric field should allow PIE" );
        }
        // three PIE fields MERGE into ONE pie column
        final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<AnnotationColumns.ColumnSpec>();
        specs.add( new AnnotationColumns.ColumnSpec( "data:a", Type.PIE ) );
        specs.add( new AnnotationColumns.ColumnSpec( "data:b", Type.PIE ) );
        specs.add( new AnnotationColumns.ColumnSpec( "data:c", Type.PIE ) );
        final AnnotationColumns pie = new AnnotationColumns( phy, specs );
        if ( pie.size() != 1 ) {
            return fail( "three PIE fields must MERGE into one column, got " + pie.size() );
        }
        if ( ( pie.getColumn( 0 ).getType() != Type.PIE ) || !"Pie".equals( pie.getColumn( 0 ).getHeader() ) ) {
            return fail( "the merged column should be a PIE with the generic 'Pie' header" );
        }
        if ( ( pie.stackColors( 0 ).size() != 3 ) || ( pie.stackHeaders( 0 ).size() != 3 )
                || !"A".equals( pie.stackHeaders( 0 ).get( 0 ) ) || !"C".equals( pie.stackHeaders( 0 ).get( 2 ) ) ) {
            return fail( "a 3-series pie should expose 3 colours + headers in series order, got " + pie.stackHeaders( 0 ) );
        }
        // wedge fractions = the tip's OWN proportions (a pie is inherently normalized), summing to 1
        final double[] pf = pie.stackFractions( tip( phy, "P" ), 0 );
        if ( !approx( pf, new double[] { 1.0 / 4, 2.0 / 4, 1.0 / 4 } ) || ( Math.abs( sum( pf ) - 1.0 ) > 1e-9 ) ) {
            return fail( "pie P wedges should be the tip's own proportions (sum 1): " + java.util.Arrays.toString( pf ) );
        }
        if ( sum( pie.stackFractions( tip( phy, "R" ), 0 ) ) != 0.0 ) {
            return fail( "a tip missing all series should have no pie (all-zero fractions)" );
        }
        if ( pie.stackFractions( tip( phy, "N" ), 0 )[ 0 ] != 0.0 ) {
            return fail( "a negative series value should count as 0 in a pie" );
        }
        // STACKED_BAR and PIE fields form TWO SEPARATE merged columns (in first-seen-per-type order)
        final List<AnnotationColumns.ColumnSpec> mixed = new ArrayList<AnnotationColumns.ColumnSpec>();
        mixed.add( new AnnotationColumns.ColumnSpec( "data:a", Type.STACKED_BAR ) );
        mixed.add( new AnnotationColumns.ColumnSpec( "data:b", Type.STACKED_BAR ) );
        mixed.add( new AnnotationColumns.ColumnSpec( "data:a", Type.PIE ) );
        mixed.add( new AnnotationColumns.ColumnSpec( "data:b", Type.PIE ) );
        final AnnotationColumns both = new AnnotationColumns( phy, mixed );
        if ( ( both.size() != 2 ) || ( both.getColumn( 0 ).getType() != Type.STACKED_BAR )
                || ( both.getColumn( 1 ).getType() != Type.PIE ) ) {
            return fail( "STACKED_BAR + PIE fields should form 2 merged columns [STACKED_BAR, PIE], got size "
                    + both.size() );
        }
        return true;
    }

    // ---- STACKED_BAR column: several numeric fields MERGE into one segmented bar, absolute or normalized ----
    private static boolean testStackedBar() {
        final Phylogeny phy = stackTree();
        if ( !"Stacked bar".equals( AnnotationColumns.label( Type.STACKED_BAR ) ) ) {
            return fail( "STACKED_BAR label should be 'Stacked bar'" );
        }
        // offered only for a numeric field
        if ( !AnnotationColumns.allowedTypes( phy, "data:a" ).contains( Type.STACKED_BAR ) ) {
            return fail( "a numeric field should allow STACKED_BAR" );
        }
        // three STACKED_BAR fields MERGE into ONE column, at the first stacked field's position
        final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<AnnotationColumns.ColumnSpec>();
        specs.add( new AnnotationColumns.ColumnSpec( "data:a", Type.STACKED_BAR ) );
        specs.add( new AnnotationColumns.ColumnSpec( "data:b", Type.STACKED_BAR ) );
        specs.add( new AnnotationColumns.ColumnSpec( "data:c", Type.STACKED_BAR ) );
        final AnnotationColumns abs = new AnnotationColumns( phy, specs );
        if ( abs.size() != 1 ) {
            return fail( "three STACKED_BAR fields must MERGE into one column, got " + abs.size() );
        }
        if ( ( abs.getColumn( 0 ).getType() != Type.STACKED_BAR )
                || !"Stacked bar".equals( abs.getColumn( 0 ).getHeader() ) ) {
            return fail( "the merged column should be a STACKED_BAR with the generic 'Stacked bar' header" );
        }
        // series colours + headers: 3, distinct, parallel and in series order
        if ( ( abs.stackColors( 0 ).size() != 3 ) || ( abs.stackHeaders( 0 ).size() != 3 ) ) {
            return fail( "a 3-series stacked bar should expose 3 colours and 3 headers" );
        }
        if ( abs.stackColors( 0 ).get( 0 ).equals( abs.stackColors( 0 ).get( 1 ) )
                || abs.stackColors( 0 ).get( 1 ).equals( abs.stackColors( 0 ).get( 2 ) )
                || abs.stackColors( 0 ).get( 0 ).equals( abs.stackColors( 0 ).get( 2 ) ) ) {
            return fail( "the stacked-bar series should have distinct colours" );
        }
        if ( !"A".equals( abs.stackHeaders( 0 ).get( 0 ) ) || !"B".equals( abs.stackHeaders( 0 ).get( 1 ) )
                || !"C".equals( abs.stackHeaders( 0 ).get( 2 ) ) ) {
            return fail( "series headers should be the prettified field names in series order, got "
                    + abs.stackHeaders( 0 ) );
        }
        // ABSOLUTE fractions: value / the LARGEST per-tip total (max total = 8 at Q). P totals 4 -> half-width bar.
        if ( !approx( abs.stackFractions( tip( phy, "P" ), 0 ), new double[] { 1.0 / 8, 2.0 / 8, 1.0 / 8 } ) ) {
            return fail( "absolute P fractions wrong: "
                    + java.util.Arrays.toString( abs.stackFractions( tip( phy, "P" ), 0 ) ) );
        }
        final double[] qf = abs.stackFractions( tip( phy, "Q" ), 0 );
        if ( !approx( qf, new double[] { 4.0 / 8, 4.0 / 8, 0.0 } ) || ( Math.abs( sum( qf ) - 1.0 ) > 1e-9 ) ) {
            return fail( "the largest-total tip's absolute bar should fill the width (sum 1): "
                    + java.util.Arrays.toString( qf ) );
        }
        // a tip missing ALL series -> all-zero (draws nothing), distinct from a full bar
        if ( sum( abs.stackFractions( tip( phy, "R" ), 0 ) ) != 0.0 ) {
            return fail( "a tip missing all series should have all-zero stacked fractions" );
        }
        // a NEGATIVE series value counts as 0 (a stacked bar is non-negative)
        if ( abs.stackFractions( tip( phy, "N" ), 0 )[ 0 ] != 0.0 ) {
            return fail( "a negative series value should count as 0 in a stacked bar" );
        }
        // NORMALIZED: value / the tip's OWN total -> every tip WITH data fills the width (sum 1)
        final List<AnnotationColumns.ColumnSpec> nspecs = new ArrayList<AnnotationColumns.ColumnSpec>();
        nspecs.add( new AnnotationColumns.ColumnSpec( "data:a", Type.STACKED_BAR, true ) );
        nspecs.add( new AnnotationColumns.ColumnSpec( "data:b", Type.STACKED_BAR, true ) );
        nspecs.add( new AnnotationColumns.ColumnSpec( "data:c", Type.STACKED_BAR, true ) );
        final AnnotationColumns norm = new AnnotationColumns( phy, nspecs );
        final double[] npf = norm.stackFractions( tip( phy, "P" ), 0 );
        if ( !approx( npf, new double[] { 1.0 / 4, 2.0 / 4, 1.0 / 4 } ) || ( Math.abs( sum( npf ) - 1.0 ) > 1e-9 ) ) {
            return fail( "normalized P fractions should be the tip's own proportions (sum 1): "
                    + java.util.Arrays.toString( npf ) );
        }
        if ( Math.abs( sum( norm.stackFractions( tip( phy, "Q" ), 0 ) ) - 1.0 ) > 1e-9 ) {
            return fail( "normalized Q bar should also fill the width" );
        }
        if ( sum( norm.stackFractions( tip( phy, "R" ), 0 ) ) != 0.0 ) {
            return fail( "normalized: a tip with no data should still draw nothing" );
        }
        // a NON-stacked column exposes no stacked fractions/colours/headers
        final List<AnnotationColumns.ColumnSpec> bar_only = new ArrayList<AnnotationColumns.ColumnSpec>();
        bar_only.add( new AnnotationColumns.ColumnSpec( "data:a", Type.BAR ) );
        final AnnotationColumns bar = new AnnotationColumns( phy, bar_only );
        if ( ( bar.stackFractions( tip( phy, "P" ), 0 ).length != 0 ) || !bar.stackColors( 0 ).isEmpty()
                || !bar.stackHeaders( 0 ).isEmpty() ) {
            return fail( "a non-stacked column should expose no stacked fractions/colours/headers" );
        }
        return true;
    }

    private static Phylogeny stackTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( stackLeaf( "P", "1", "2", "1" ) );  // total 4
        root.addAsChild( stackLeaf( "Q", "4", "4", "0" ) );  // total 8 (the max)
        root.addAsChild( stackLeaf( "R", null, null, null ) ); // missing all
        root.addAsChild( stackLeaf( "N", "-1", "2", "1" ) ); // negative a -> counts as 0 (total 3)
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode stackLeaf( final String name, final String a, final String b, final String c ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        if ( ( a != null ) || ( b != null ) || ( c != null ) ) {
            final PropertiesList pl = new PropertiesList();
            if ( a != null ) {
                pl.addProperty( new Property( "data:a", a, "", "xsd:string", AppliesTo.NODE ) );
            }
            if ( b != null ) {
                pl.addProperty( new Property( "data:b", b, "", "xsd:string", AppliesTo.NODE ) );
            }
            if ( c != null ) {
                pl.addProperty( new Property( "data:c", c, "", "xsd:string", AppliesTo.NODE ) );
            }
            n.getNodeData().setProperties( pl );
        }
        return n;
    }

    private static boolean approx( final double[] a, final double[] b ) {
        if ( ( a == null ) || ( b == null ) || ( a.length != b.length ) ) {
            return false;
        }
        for( int i = 0; i < a.length; ++i ) {
            if ( Math.abs( a[ i ] - b[ i ] ) > 1e-9 ) {
                return false;
            }
        }
        return true;
    }

    private static double sum( final double[] a ) {
        double s = 0;
        for( final double v : a ) {
            s += v;
        }
        return s;
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
        // shape picker: labels; default CIRCLE (2-arg + null coerced); the chosen shape rides into the column
        if ( !"Circle".equals( AnnotationColumns.shapeLabel( AnnotationColumns.SymbolShape.CIRCLE ) )
                || !"Square".equals( AnnotationColumns.shapeLabel( AnnotationColumns.SymbolShape.SQUARE ) )
                || !"Diamond".equals( AnnotationColumns.shapeLabel( AnnotationColumns.SymbolShape.DIAMOND ) )
                || !"Triangle".equals( AnnotationColumns.shapeLabel( AnnotationColumns.SymbolShape.TRIANGLE ) ) ) {
            return fail( "symbol-shape labels are wrong" );
        }
        if ( ( new AnnotationColumns.ColumnSpec( "data:present", Type.SYMBOL )._shape
                != AnnotationColumns.SymbolShape.CIRCLE )
                || ( new AnnotationColumns.ColumnSpec( "data:present", Type.SYMBOL, null )._shape
                        != AnnotationColumns.SymbolShape.CIRCLE ) ) {
            return fail( "a SYMBOL column should default to a CIRCLE glyph (2-arg / null shape)" );
        }
        final List<AnnotationColumns.ColumnSpec> shaped = new ArrayList<AnnotationColumns.ColumnSpec>();
        shaped.add( new AnnotationColumns.ColumnSpec( "data:present", Type.SYMBOL,
                AnnotationColumns.SymbolShape.TRIANGLE ) );
        if ( new AnnotationColumns( bphy, shaped ).symbolShape( 0 ) != AnnotationColumns.SymbolShape.TRIANGLE ) {
            return fail( "the chosen symbol shape should ride into the resolved column" );
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
                || !"Stacked bar".equals( AnnotationColumns.label( Type.STACKED_BAR ) )
                || !"Pie".equals( AnnotationColumns.label( Type.PIE ) )
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
                || host_types.contains( Type.HEATMAP ) || host_types.contains( Type.STACKED_BAR )
                || host_types.contains( Type.PIE ) ) {
            return fail( "a categorical field should allow COLOR_STRIP + SYMBOL + TEXT (no STACKED_BAR/PIE)" );
        }
        final List<Type> score_types = AnnotationColumns.allowedTypes( phy, "data:score" );
        if ( ( score_types.size() != 6 ) || !score_types.contains( Type.HEATMAP )
                || !score_types.contains( Type.MATRIX ) || !score_types.contains( Type.BAR )
                || !score_types.contains( Type.STACKED_BAR ) || !score_types.contains( Type.PIE )
                || !score_types.contains( Type.TEXT ) || score_types.contains( Type.COLOR_STRIP )
                || score_types.contains( Type.SYMBOL ) ) {
            return fail( "a numeric field should allow HEATMAP + MATRIX + BAR + STACKED_BAR + PIE + TEXT" );
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
