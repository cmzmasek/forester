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
        return testTypeSuggestion() && testCells() && testHeaders() && testLabels();
    }

    // ---- friendly render-type labels ----
    private static boolean testLabels() {
        if ( !"Color strip".equals( AnnotationColumns.label( Type.COLOR_STRIP ) )
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
        if ( ( host_types.size() != 2 ) || !host_types.contains( Type.COLOR_STRIP )
                || !host_types.contains( Type.TEXT ) || host_types.contains( Type.HEATMAP ) ) {
            return fail( "a categorical field should allow only COLOR_STRIP + TEXT" );
        }
        final List<Type> score_types = AnnotationColumns.allowedTypes( phy, "data:score" );
        if ( ( score_types.size() != 3 ) || !score_types.contains( Type.HEATMAP )
                || !score_types.contains( Type.BAR ) || !score_types.contains( Type.TEXT )
                || score_types.contains( Type.COLOR_STRIP ) ) {
            return fail( "a numeric field should allow HEATMAP + BAR + TEXT" );
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
