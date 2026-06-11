// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// Contact: phylosoft @ gmail . com

package org.forester.archaeopteryx;

import java.awt.Color;
import java.util.Map;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;

/**
 * Unit tests for {@link PropertyColorScheme}. Lives in the {@code org.forester.archaeopteryx}
 * package because the class under test (and its methods) are package-private. Run standalone
 * via {@link #main(String[])}, or as part of the suite via {@link #test()} (called from
 * {@code org.forester.test.Test}).
 */
public final class PropertyColorSchemeTest {

    public static void main( final String[] args ) {
        System.out.println( "PropertyColorScheme: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        return testDisplayName() && testCategoricalGrouping() && testCountryGrouping() && testYearGradient()
                && testAbsentAndEmpty() && testCollapseExcludesHiddenLeaves() && testCollapseRescalesGradient();
    }

    // ---- displayName: namespace strip, '_' -> space, capitalize, acronyms preserved ----
    private static boolean testDisplayName() {
        if ( !eq( "Protein Names", PropertyColorScheme.displayName( "repseq:protein_names" ), "displayName protein_names" ) ) {
            return false;
        }
        if ( !eq( "Host", PropertyColorScheme.displayName( "repseq:host" ), "displayName host" ) ) {
            return false;
        }
        if ( !eq( "RNA Type", PropertyColorScheme.displayName( "RNA_type" ), "displayName acronym" ) ) {
            return false;
        }
        if ( !eq( "Year", PropertyColorScheme.displayName( "year" ), "displayName no-prefix" ) ) {
            return false;
        }
        if ( !eq( "No Namespace Here", PropertyColorScheme.displayName( "no_namespace_here" ), "displayName no-namespace" ) ) {
            return false;
        }
        if ( !eq( "None", PropertyColorScheme.displayName( "None" ), "displayName already-clean" ) ) {
            return false;
        }
        // null / empty are returned as-is (must not throw)
        if ( PropertyColorScheme.displayName( null ) != null ) {
            return fail( "displayName(null) should be null" );
        }
        if ( !eq( "", PropertyColorScheme.displayName( "" ), "displayName empty" ) ) {
            return false;
        }
        return true;
    }

    // ---- categorical grouping: case / whitespace / underscore variants share a color ----
    private static boolean testCategoricalGrouping() {
        final String ref = "repseq:host";
        // "Human" x3, "human" x1, "homo_sapiens" x1, "Homo sapiens" x2, "man" x1
        final Phylogeny phy = treeWith( ref, "Human", "Human", "Human", "human", "homo_sapiens", "Homo sapiens",
                                        "Homo sapiens", "man" );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        if ( s.isGradient() ) {
            return fail( "host should not be a gradient" );
        }
        // three groups: {Human, human}, {homo_sapiens, Homo sapiens}, {man}
        if ( s.getValueColors().size() != 3 ) {
            return fail( "host groups expected 3, got " + s.getValueColors().size() );
        }
        // trivial variants -> same color
        if ( !sameColor( s, phy, ref, "Human", "human" ) ) {
            return fail( "Human/human should share a color" );
        }
        if ( !sameColor( s, phy, ref, "homo_sapiens", "Homo sapiens" ) ) {
            return fail( "homo_sapiens/Homo sapiens should share a color" );
        }
        // semantically equal but lexically different -> intentionally NOT merged
        if ( sameColor( s, phy, ref, "Human", "man" ) ) {
            return fail( "Human/man must NOT be merged (out of scope)" );
        }
        if ( sameColor( s, phy, ref, "Human", "homo_sapiens" ) ) {
            return fail( "Human/homo_sapiens must NOT be merged" );
        }
        // legend shows the most frequent spelling per group
        final Map<String, Color> legend = s.getValueColors();
        if ( !legend.containsKey( "Human" ) || legend.containsKey( "human" ) ) {
            return fail( "legend should show 'Human' (most frequent), not 'human'" );
        }
        if ( !legend.containsKey( "Homo sapiens" ) || legend.containsKey( "homo_sapiens" ) ) {
            return fail( "legend should show 'Homo sapiens' (most frequent), not 'homo_sapiens'" );
        }
        if ( !legend.containsKey( "man" ) ) {
            return fail( "legend should contain 'man'" );
        }
        return true;
    }

    // ---- country: group by the part before the first ':' (USA:CA == USA:IL) ----
    private static boolean testCountryGrouping() {
        final String ref = "repseq:country";
        final Phylogeny phy = treeWith( ref, "USA:CA", "USA:IL", "usa:ny", "Canada:ON", "Canada", "Brazil:SP" );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        if ( s.getValueColors().size() != 3 ) {
            return fail( "country groups expected 3 (USA, Canada, Brazil), got " + s.getValueColors().size() );
        }
        if ( !sameColor( s, phy, ref, "USA:CA", "USA:IL" ) || !sameColor( s, phy, ref, "USA:CA", "usa:ny" ) ) {
            return fail( "USA:* subdivisions should share a color" );
        }
        if ( !sameColor( s, phy, ref, "Canada:ON", "Canada" ) ) {
            return fail( "Canada:ON / Canada should share a color" );
        }
        if ( sameColor( s, phy, ref, "USA:CA", "Brazil:SP" ) ) {
            return fail( "USA and Brazil must differ" );
        }
        // legend labels carry no ':' subdivision
        for( final String label : s.getValueColors().keySet() ) {
            if ( label.indexOf( ':' ) >= 0 ) {
                return fail( "country legend label should not contain ':' -> " + label );
            }
        }
        return true;
    }

    // ---- year: continuous gradient over the numeric range; non-numeric/missing -> null ----
    private static boolean testYearGradient() {
        final String ref = "repseq:year";
        final Phylogeny phy = treeWith( ref, "1927", "2000", "2025", null, "n/a" );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        if ( !s.isGradient() ) {
            return fail( "year should be a gradient" );
        }
        if ( !eq( "1927", s.getGradientMinLabel(), "year min label" ) ) {
            return false;
        }
        if ( !eq( "2025", s.getGradientMaxLabel(), "year max label" ) ) {
            return false;
        }
        // min-year leaf gets the low (t=0) color, max-year leaf the high (t=1) color
        final Color c_min = colorForValue( s, phy, ref, "1927" );
        final Color c_max = colorForValue( s, phy, ref, "2025" );
        if ( ( c_min == null ) || !c_min.equals( s.gradientColorAt( 0.0 ) ) ) {
            return fail( "1927 should map to gradientColorAt(0)" );
        }
        if ( ( c_max == null ) || !c_max.equals( s.gradientColorAt( 1.0 ) ) ) {
            return fail( "2025 should map to gradientColorAt(1)" );
        }
        if ( c_min.equals( c_max ) ) {
            return fail( "min and max year colors should differ" );
        }
        // non-numeric and missing -> no color
        if ( colorForValue( s, phy, ref, "n/a" ) != null ) {
            return fail( "non-numeric year should have no color" );
        }
        if ( colorForMissing( s, phy, ref ) != null ) {
            return fail( "leaf without a year should have no color" );
        }
        // gradient legend is not a categorical value list
        if ( !s.getValueColors().isEmpty() ) {
            return fail( "gradient scheme should have no categorical legend entries" );
        }
        return true;
    }

    // ---- absent property / empty tree -> empty scheme ----
    private static boolean testAbsentAndEmpty() {
        final Phylogeny phy = treeWith( "repseq:host", "cat", "dog" );
        final PropertyColorScheme absent = new PropertyColorScheme( phy, "repseq:not_present" );
        if ( !absent.isEmpty() ) {
            return fail( "scheme for an absent property should be empty" );
        }
        final PropertyColorScheme empty_year = new PropertyColorScheme( phy, "repseq:year" );
        if ( !empty_year.isEmpty() ) {
            return fail( "year gradient with no numeric values should be empty" );
        }
        if ( new PropertyColorScheme( null, "repseq:host" ).isEmpty() != true ) {
            return fail( "scheme over a null phylogeny should be empty" );
        }
        return true;
    }

    // ---- collapsing a clade drops its (now hidden) leaves from the categorical legend ----
    private static boolean testCollapseExcludesHiddenLeaves() {
        final String ref = "repseq:host";
        //        root
        //       /    \
        //   cladeA   cladeB
        //   /   \     /   \
        // cat  dog  fish  bird
        final PhylogenyNode clade_b = internal( leaf( "b1", ref, "fish" ), leaf( "b2", ref, "bird" ) );
        final Phylogeny phy = treeOf( internal( leaf( "a1", ref, "cat" ), leaf( "a2", ref, "dog" ) ), clade_b );
        // nothing collapsed: all four values present
        if ( new PropertyColorScheme( phy, ref ).getValueColors().size() != 4 ) {
            return fail( "expected 4 host groups while nothing is collapsed" );
        }
        // collapse cladeB: its leaves (fish, bird) are hidden and must drop out
        clade_b.setCollapse( true );
        final PropertyColorScheme collapsed = new PropertyColorScheme( phy, ref );
        if ( collapsed.getValueColors().size() != 2 ) {
            return fail( "expected 2 host groups after collapsing a clade, got " + collapsed.getValueColors().size() );
        }
        final Map<String, Color> legend = collapsed.getValueColors();
        if ( !legend.containsKey( "cat" ) || !legend.containsKey( "dog" ) ) {
            return fail( "visible leaves (cat, dog) should remain in the legend" );
        }
        if ( legend.containsKey( "fish" ) || legend.containsKey( "bird" ) ) {
            return fail( "collapsed-away leaves (fish, bird) should be gone from the legend" );
        }
        // a hidden leaf gets no color
        if ( colorForValue( collapsed, phy, ref, "fish" ) != null ) {
            return fail( "a collapsed-away leaf should have no color" );
        }
        // uncollapsing restores all four
        clade_b.setCollapse( false );
        if ( new PropertyColorScheme( phy, ref ).getValueColors().size() != 4 ) {
            return fail( "uncollapsing should restore all 4 host groups" );
        }
        return true;
    }

    // ---- collapsing a clade rescales the year gradient to the still-visible range ----
    private static boolean testCollapseRescalesGradient() {
        final String ref = "repseq:year";
        final PhylogenyNode older = internal( leaf( "b1", ref, "1950" ), leaf( "b2", ref, "1960" ) );
        final Phylogeny phy = treeOf( internal( leaf( "a1", ref, "2000" ), leaf( "a2", ref, "2010" ) ), older );
        final PropertyColorScheme full = new PropertyColorScheme( phy, ref );
        if ( !eq( "1950", full.getGradientMinLabel(), "year min uncollapsed" )
                || !eq( "2010", full.getGradientMaxLabel(), "year max uncollapsed" ) ) {
            return false;
        }
        // collapse the older clade: the gradient rescales to the visible 2000..2010
        older.setCollapse( true );
        final PropertyColorScheme collapsed = new PropertyColorScheme( phy, ref );
        if ( !eq( "2000", collapsed.getGradientMinLabel(), "year min after collapse" )
                || !eq( "2010", collapsed.getGradientMaxLabel(), "year max after collapse" ) ) {
            return false;
        }
        return true;
    }

    // ---------------------------------------------------------------------------------------

    /** A leaf node named {@code name}; a {@code null} value means "no property". */
    private static PhylogenyNode leaf( final String name, final String ref, final String value ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        if ( value != null ) {
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( ref, value, "", "xsd:string", AppliesTo.NODE ) );
            n.getNodeData().setProperties( pl );
        }
        return n;
    }

    /** An internal node with the given children. */
    private static PhylogenyNode internal( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode c : children ) {
            n.addAsChild( c );
        }
        return n;
    }

    /** A tree rooted at a new node with the given (internal) clades as children. */
    private static Phylogeny treeOf( final PhylogenyNode... clades ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final PhylogenyNode c : clades ) {
            root.addAsChild( c );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /** A flat tree with one external node per value; a {@code null} value means "no property". */
    private static Phylogeny treeWith( final String ref, final String... values ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        int i = 0;
        for( final String v : values ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "leaf" + ( i++ ) );
            if ( v != null ) {
                final PropertiesList pl = new PropertiesList();
                pl.addProperty( new Property( ref, v, "", "xsd:string", AppliesTo.NODE ) );
                n.getNodeData().setProperties( pl );
            }
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean sameColor( final PropertyColorScheme s, final Phylogeny phy, final String ref,
                                      final String raw_a, final String raw_b ) {
        final Color a = colorForValue( s, phy, ref, raw_a );
        final Color b = colorForValue( s, phy, ref, raw_b );
        return ( a != null ) && a.equals( b );
    }

    private static Color colorForValue( final PropertyColorScheme s, final Phylogeny phy, final String ref,
                                        final String raw ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( raw.equals( valueOf( n, ref ) ) ) {
                return s.colorFor( n );
            }
        }
        return null;
    }

    private static Color colorForMissing( final PropertyColorScheme s, final Phylogeny phy, final String ref ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( valueOf( n, ref ) == null ) {
                return s.colorFor( n );
            }
        }
        return null;
    }

    private static String valueOf( final PhylogenyNode n, final String ref ) {
        if ( ( n.getNodeData() == null ) || ( n.getNodeData().getProperties() == null ) ) {
            return null;
        }
        final java.util.List<Property> ps = n.getNodeData().getProperties().getProperties( ref );
        return ps.isEmpty() ? null : ps.get( 0 ).getValue();
    }

    private static boolean eq( final String expected, final String actual, final String what ) {
        if ( ( expected == null ) ? ( actual != null ) : !expected.equals( actual ) ) {
            return fail( what + ": expected [" + expected + "] got [" + actual + "]" );
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [PropertyColorSchemeTest] " + message );
        return false;
    }

    private PropertyColorSchemeTest() {
        // not instantiable
    }
}
