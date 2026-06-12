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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
        return testDisplayName() && testCategoricalGrouping() && testHumanSynonym() && testCountryGrouping()
                && testHostQualifierGrouping() && testYearGradient() && testAbsentAndEmpty()
                && testCollapseExcludesHiddenLeaves() && testCollapseRescalesGradient()
                && testFrequencyColorsAndLegend();
    }

    // ---- colors assigned by frequency (distinct for the most common values); legend = top-N most
    //      frequent, re-sorted alphabetically ----
    private static boolean testFrequencyColorsAndLegend() {
        final String ref = "repseq:host";
        // 26 distinct values "x00".."x25"; value x_i occurs (i+1) times, so frequency INCREASES with
        // i (x25 most frequent, x00 least) -- the reverse of alphabetical order.
        final List<String> vals = new ArrayList<String>();
        for( int i = 0; i <= 25; ++i ) {
            final String name = String.format( "x%02d", i );
            for( int c = 0; c <= i; ++c ) {
                vals.add( name );
            }
        }
        final Phylogeny phy = treeWith( ref, vals.toArray( new String[ 0 ] ) );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        if ( s.numberOfValues() != 26 ) {
            return fail( "expected 26 distinct values, got " + s.numberOfValues() );
        }
        // the legend shows the 20 MOST FREQUENT values (x06..x25), re-sorted alphabetically
        final Map<String, Color> legend = s.legendValues( 20 );
        final List<String> keys = new ArrayList<String>( legend.keySet() );
        if ( keys.size() != 20 ) {
            return fail( "legend should hold 20 entries, got " + keys.size() );
        }
        for( int k = 0; k < 20; ++k ) {
            final String expected = String.format( "x%02d", k + 6 ); // x06..x25, in alphabetical order
            if ( !expected.equals( keys.get( k ) ) ) {
                return fail( "legend entry " + k + " expected " + expected + " got " + keys.get( k ) );
            }
        }
        // the 24 most frequent values (x02..x25) must all have distinct colors (no palette cycling)
        final Set<Color> colors = new HashSet<Color>();
        for( int k = 2; k <= 25; ++k ) {
            colors.add( colorForValue( s, phy, ref, String.format( "x%02d", k ) ) );
        }
        if ( colors.size() != 24 ) {
            return fail( "the 24 most frequent values should have 24 distinct colors, got " + colors.size() );
        }
        // per-value leaf counts (for the legend): value x_i occurs (i+1) times
        final Map<String, Integer> ct = s.getValueCounts();
        if ( ct.size() != 26 ) {
            return fail( "expected 26 per-value counts, got " + ct.size() );
        }
        for( int k = 0; k <= 25; ++k ) {
            final String name = String.format( "x%02d", k );
            if ( ( ct.get( name ) == null ) || ( ct.get( name ).intValue() != ( k + 1 ) ) ) {
                return fail( "count for " + name + " expected " + ( k + 1 ) + " got " + ct.get( name ) );
            }
        }
        return true;
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
        // "Mouse" x3, "mouse" x1, "mus_musculus" x1, "Mus musculus" x2, "rat" x1
        final Phylogeny phy = treeWith( ref, "Mouse", "Mouse", "Mouse", "mouse", "mus_musculus", "Mus musculus",
                                        "Mus musculus", "rat" );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        if ( s.isGradient() ) {
            return fail( "host should not be a gradient" );
        }
        // three groups: {Mouse, mouse}, {mus_musculus, Mus musculus}, {rat}
        if ( s.getValueColors().size() != 3 ) {
            return fail( "host groups expected 3, got " + s.getValueColors().size() );
        }
        // trivial variants -> same color
        if ( !sameColor( s, phy, ref, "Mouse", "mouse" ) ) {
            return fail( "Mouse/mouse should share a color" );
        }
        if ( !sameColor( s, phy, ref, "mus_musculus", "Mus musculus" ) ) {
            return fail( "mus_musculus/Mus musculus should share a color" );
        }
        // semantically equal but lexically different -> intentionally NOT merged
        if ( sameColor( s, phy, ref, "Mouse", "rat" ) ) {
            return fail( "Mouse/rat must NOT be merged (out of scope)" );
        }
        if ( sameColor( s, phy, ref, "Mouse", "mus_musculus" ) ) {
            return fail( "Mouse/mus_musculus must NOT be merged" );
        }
        // legend shows the most frequent spelling per group
        final Map<String, Color> legend = s.getValueColors();
        if ( !legend.containsKey( "Mouse" ) || legend.containsKey( "mouse" ) ) {
            return fail( "legend should show 'Mouse' (most frequent), not 'mouse'" );
        }
        if ( !legend.containsKey( "Mus musculus" ) || legend.containsKey( "mus_musculus" ) ) {
            return fail( "legend should show 'Mus musculus' (most frequent), not 'mus_musculus'" );
        }
        if ( !legend.containsKey( "rat" ) ) {
            return fail( "legend should contain 'rat'" );
        }
        return true;
    }

    // ---- the one deliberate synonym fold: human/Human/humans -> Homo sapiens ----
    private static boolean testHumanSynonym() {
        final String ref = "repseq:host";
        final Phylogeny phy = treeWith( ref, "Human", "human", "Homo sapiens", "homo_sapiens", "man" );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        // two groups: the merged { Human, human, Homo sapiens, homo_sapiens } and { man }
        if ( s.getValueColors().size() != 2 ) {
            return fail( "human synonym: expected 2 groups, got " + s.getValueColors().size() );
        }
        if ( !sameColor( s, phy, ref, "Human", "Homo sapiens" ) || !sameColor( s, phy, ref, "human", "homo_sapiens" ) ) {
            return fail( "human/Human should fold into Homo sapiens" );
        }
        if ( sameColor( s, phy, ref, "Human", "man" ) ) {
            return fail( "'man' must NOT fold into Homo sapiens (only 'human' is special-cased)" );
        }
        // the legend shows the canonical "Homo sapiens", not "Human"/"human"
        final Map<String, Color> legend = s.getValueColors();
        if ( !legend.containsKey( "Homo sapiens" ) || legend.containsKey( "Human" ) || legend.containsKey( "human" ) ) {
            return fail( "legend should show 'Homo sapiens', not 'Human'/'human'" );
        }
        // the merged group counts all four folded/variant leaves
        final Integer count = s.getValueCounts().get( "Homo sapiens" );
        if ( ( count == null ) || ( count.intValue() != 4 ) ) {
            return fail( "Homo sapiens count should be 4, got " + count );
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

    // ---- host: drop the qualifier after the first ';' (Homo sapiens; male == Homo sapiens) ----
    private static boolean testHostQualifierGrouping() {
        final String ref = "repseq:host";
        final Phylogeny phy = treeWith( ref, "Homo sapiens; male 35", "Homo sapiens; female old", "Homo sapiens",
                                        "homo_sapiens; juvenile", "Mus musculus; female", "Gallus gallus" );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        // three groups: {the four Homo sapiens}, {Mus musculus}, {Gallus gallus}
        if ( s.getValueColors().size() != 3 ) {
            return fail( "host groups expected 3 (Homo sapiens, Mus musculus, Gallus gallus), got "
                    + s.getValueColors().size() );
        }
        // same base host with different qualifiers -- and underscore/case variants -- share a color
        if ( !sameColor( s, phy, ref, "Homo sapiens; male 35", "Homo sapiens; female old" )
                || !sameColor( s, phy, ref, "Homo sapiens; male 35", "Homo sapiens" )
                || !sameColor( s, phy, ref, "Homo sapiens; male 35", "homo_sapiens; juvenile" ) ) {
            return fail( "Homo sapiens with different qualifiers should share a color" );
        }
        // different base hosts stay distinct
        if ( sameColor( s, phy, ref, "Homo sapiens; male 35", "Mus musculus; female" ) ) {
            return fail( "Homo sapiens and Mus musculus must differ" );
        }
        // legend labels carry no ';' qualifier and show the most frequent spelling
        for( final String label : s.getValueColors().keySet() ) {
            if ( label.indexOf( ';' ) >= 0 ) {
                return fail( "host legend label should not contain ';' -> " + label );
            }
        }
        if ( !s.getValueColors().containsKey( "Homo sapiens" ) ) {
            return fail( "host legend should show 'Homo sapiens' (most frequent spelling)" );
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
