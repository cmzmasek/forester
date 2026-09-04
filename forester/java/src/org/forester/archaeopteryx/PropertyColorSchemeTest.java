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
import java.util.HashMap;
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
                && testHostQualifierGrouping() && testYearGradient() && testNumericGradientGeneralization()
                && testColorableRefs() && testAbsentAndEmpty() && testCollapseExcludesHiddenLeaves()
                && testCollapseRescalesGradient() && testFrequencyColorsAndLegend() && testColorOverrides()
                && testPalettes() && testOrderLegendEntriesEdges() && testMissingCount()
                && testColorIdentityMemory();
    }

    // ---- orderLegendEntries edge cases (formerly covered via legendValues/capEntries) ----
    private static boolean testOrderLegendEntriesEdges() {
        final Map<String, Color> colors = new HashMap<String, Color>();
        colors.put( "b", Color.RED );
        colors.put( "a", Color.GREEN );
        colors.put( "c", Color.BLUE );
        // null counts -> every key ranks as 0 -> the cap keeps the alphabetically-first max, displayed A-Z
        final List<String> nc = new ArrayList<String>( TreePanel.orderLegendEntries( colors, null, 2, false ).keySet() );
        if ( ( nc.size() != 2 ) || !"a".equals( nc.get( 0 ) ) || !"b".equals( nc.get( 1 ) ) ) {
            return fail( "null counts, max 2: keep a,b A-Z; got " + nc );
        }
        // max >= size keeps everything
        if ( TreePanel.orderLegendEntries( colors, null, 99, false ).size() != 3 ) {
            return fail( "max >= size must keep all entries" );
        }
        // max == 0 -> empty
        if ( !TreePanel.orderLegendEntries( colors, null, 0, false ).isEmpty() ) {
            return fail( "max 0 must yield an empty legend" );
        }
        // partial counts: only "c" has a count, so by_count keeps it and ranks it first
        final Map<String, Integer> partial = new HashMap<String, Integer>();
        partial.put( "c", 5 );
        final List<String> byc = new ArrayList<String>(
                TreePanel.orderLegendEntries( colors, partial, 2, true ).keySet() );
        if ( ( byc.size() != 2 ) || !"c".equals( byc.get( 0 ) ) ) {
            return fail( "partial counts, by_count: c (count 5) must lead; got " + byc );
        }
        return true;
    }

    // ---- colorableRefs: numeric all-unique columns stay colorable (gradient); categorical all-unique don't ----
    private static boolean testColorableRefs() {
        // a numeric column where every tip has a distinct value is still colorable (renders as a gradient)
        if ( !PropertyColorScheme.colorableRefs( treeWith( "data:year", "2015", "2016", "2020", "2024" ) )
                .contains( "data:year" ) ) {
            return fail( "an all-distinct numeric column should be colorable (gradient)" );
        }
        // a categorical column where every tip is unique is NOT colorable (one color per tip is useless)
        if ( PropertyColorScheme.colorableRefs( treeWith( "data:id", "a", "b", "c", "d" ) ).contains( "data:id" ) ) {
            return fail( "an all-distinct categorical column should not be colorable" );
        }
        // a repeated categorical column is colorable
        if ( !PropertyColorScheme.colorableRefs( treeWith( "repseq:host", "cat", "cat", "dog", "dog" ) )
                .contains( "repseq:host" ) ) {
            return fail( "a repeated categorical column should be colorable" );
        }
        // a constant column (one distinct value) is not colorable
        if ( !PropertyColorScheme.colorableRefs( treeWith( "data:const", "X", "X", "X", "X" ) ).isEmpty() ) {
            return fail( "a constant column should not be colorable" );
        }
        return true;
    }

    // ---- gradient generalizes beyond "year": any predominantly-numeric column is colored by a gradient ----
    private static boolean testNumericGradientGeneralization() {
        // a non-"year" numeric column -> gradient over its range
        final String age = "data:age";
        final Phylogeny p1 = treeWith( age, "1", "2", "3", "4", "5" );
        final PropertyColorScheme s1 = new PropertyColorScheme( p1, age );
        if ( !s1.isGradient() ) {
            return fail( "a numeric column (age) should be a gradient" );
        }
        if ( !eq( "1", s1.getGradientMinLabel(), "age min" ) || !eq( "5", s1.getGradientMaxLabel(), "age max" ) ) {
            return false;
        }
        // predominantly numeric with a non-numeric sentinel (3 of 4) -> still a gradient; sentinel uncolored
        final String pct = "data:pct_identity";
        final Phylogeny p2 = treeWith( pct, "88.5", "91.0", "97.2", "unknown" );
        final PropertyColorScheme s2 = new PropertyColorScheme( p2, pct );
        if ( !s2.isGradient() ) {
            return fail( "a 3-of-4 numeric column should still be a gradient" );
        }
        if ( colorForValue( s2, p2, pct, "unknown" ) != null ) {
            return fail( "the non-numeric sentinel should have no gradient color" );
        }
        // mostly textual with a stray number (1 of 4) -> categorical, one color per distinct value
        final String mixed = "data:mixed";
        final Phylogeny p3 = treeWith( mixed, "cat", "dog", "fish", "7" );
        final PropertyColorScheme s3 = new PropertyColorScheme( p3, mixed );
        if ( s3.isGradient() ) {
            return fail( "a mostly-textual column should stay categorical" );
        }
        if ( s3.getValueColors().size() != 4 ) {
            return fail( "mixed column should have 4 categorical groups, got " + s3.getValueColors().size() );
        }
        // exactly half numeric (2 of 4) is not a strict majority -> categorical
        final String half = "data:half";
        if ( new PropertyColorScheme( treeWith( half, "1", "2", "cat", "dog" ), half ).isGradient() ) {
            return fail( "a half-numeric column should not be a gradient (needs a strict majority)" );
        }
        // a single distinct number (no range) -> categorical, not a degenerate gradient
        final String flat = "data:flat";
        if ( new PropertyColorScheme( treeWith( flat, "5", "5", "5", "cat" ), flat ).isGradient() ) {
            return fail( "a column with only one distinct number should not be a gradient" );
        }
        return true;
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
        // the legend keeps the 20 MOST FREQUENT values (x06..x25); A-Z display order (by_count=false)
        final Map<String, Color> legend = TreePanel.orderLegendEntries( s.getValueColors(), s.getValueCounts(), 20,
                                                                        false );
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
        // by_count=true display order: same 20 kept, but most-frequent first (x25, x24, ..., x06)
        final List<String> by_count = new ArrayList<String>(
                TreePanel.orderLegendEntries( s.getValueColors(), s.getValueCounts(), 20, true ).keySet() );
        for( int k = 0; k < 20; ++k ) {
            final String expected = String.format( "x%02d", 25 - k ); // x25 down to x06
            if ( !expected.equals( by_count.get( k ) ) ) {
                return fail( "by-count legend entry " + k + " expected " + expected + " got " + by_count.get( k ) );
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

    // ---- selectable categorical palettes ----
    private static boolean testPalettes() {
        if ( !PropertyColorScheme.paletteNames().contains( "Default" )
                || !PropertyColorScheme.paletteNames().contains( "Colorblind-friendly" ) ) {
            return fail( "expected Default and Colorblind-friendly palettes" );
        }
        final String ref = "repseq:host";
        // "common" is the most frequent value, so it gets the first color of the chosen palette
        final Phylogeny phy = treeWith( ref, "common", "common", "common", "rare" );
        if ( !new Color( 0x4E79A7 ).equals( colorForValue( new PropertyColorScheme( phy, ref, null, "Default" ), phy,
                                                           ref, "common" ) ) ) {
            return fail( "Default palette: 'common' should get the first default (Tableau 10) color" );
        }
        if ( !new Color( 0xE69F00 ).equals( colorForValue( new PropertyColorScheme( phy, ref, null,
                                                                                    "Colorblind-friendly" ),
                                                           phy, ref, "common" ) ) ) {
            return fail( "Colorblind palette: 'common' should get the first colorblind color" );
        }
        // an unknown palette name falls back to the default
        if ( !new Color( 0x4E79A7 ).equals( colorForValue( new PropertyColorScheme( phy, ref, null, "Nonexistent" ),
                                                           phy, ref, "common" ) ) ) {
            return fail( "unknown palette name should fall back to Default" );
        }
        return true;
    }

    // ---- user-assigned per-value color overrides (keyed by group key) ----
    private static boolean testColorOverrides() {
        final String ref = "repseq:host";
        final Phylogeny phy = treeWith( ref, "cat", "cat", "dog", "fish" );
        final Map<String, Color> overrides = new HashMap<String, Color>();
        overrides.put( "cat", new Color( 0x123456 ) ); // keyed by the group key (lower-cased "cat")
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref, overrides );
        if ( !new Color( 0x123456 ).equals( colorForValue( s, phy, ref, "cat" ) ) ) {
            return fail( "cat should use the override color" );
        }
        // dog and fish keep distinct automatic palette colors, different from the override
        final Color dog = colorForValue( s, phy, ref, "dog" );
        final Color fish = colorForValue( s, phy, ref, "fish" );
        if ( new Color( 0x123456 ).equals( dog ) || new Color( 0x123456 ).equals( fish ) || dog.equals( fish ) ) {
            return fail( "dog/fish should keep distinct automatic colors" );
        }
        // getValueKeys maps a representative label to its (stable) group key
        if ( !"cat".equals( s.getValueKeys().get( "cat" ) ) ) {
            return fail( "getValueKeys should map 'cat' -> 'cat'" );
        }
        // without overrides, cat is automatic again
        if ( new Color( 0x123456 ).equals( colorForValue( new PropertyColorScheme( phy, ref, null ), phy, ref, "cat" ) ) ) {
            return fail( "without overrides, cat should use an automatic color" );
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

    // ---- the "no value" legend row's count: visible tips that draw NO mark under the scheme ----
    private static boolean testMissingCount() {
        // categorical: two of five tips carry no host at all -> missing 2
        final Phylogeny phy = treeWith( "repseq:host", "cat", "dog", "cat", null, null );
        final PropertyColorScheme cat = new PropertyColorScheme( phy, "repseq:host" );
        if ( ( cat.visibleTipCount() != 5 ) || ( cat.missingCount() != 2 ) ) {
            return fail( "categorical: expected 5 visible / 2 missing, got " + cat.visibleTipCount() + " / "
                    + cat.missingCount() );
        }
        // full coverage -> 0 (the row must not appear)
        final PropertyColorScheme full = new PropertyColorScheme( treeWith( "repseq:host", "cat", "dog" ),
                                                                  "repseq:host" );
        if ( full.missingCount() != 0 ) {
            return fail( "full coverage must count 0 missing, got " + full.missingCount() );
        }
        // gradient: an UNPARSEABLE value ("n/a") draws no mark, so it counts as missing too -- the row must
        // agree with what the tree shows, not with mere value presence
        final PropertyColorScheme grad = new PropertyColorScheme(
                treeWith( "repseq:year", "2001", "2005", "2010", "n/a", null ), "repseq:year" );
        if ( !grad.isGradient() || ( grad.missingCount() != 2 ) ) {
            return fail( "gradient: expected 2 missing (one n/a + one absent), got " + grad.missingCount() );
        }
        // collapsed-away tips are not counted either way: the row describes the DISPLAYED tree
        final PhylogenyNode hidden = internal( leaf( "b1", "repseq:host", null ),
                                               leaf( "b2", "repseq:host", null ) );
        final Phylogeny phy2 = treeOf( internal( leaf( "a1", "repseq:host", "cat" ),
                                                 leaf( "a2", "repseq:host", "dog" ) ),
                                       hidden );
        hidden.setCollapse( true );
        final PropertyColorScheme collapsed = new PropertyColorScheme( phy2, "repseq:host" );
        if ( ( collapsed.visibleTipCount() != 2 ) || ( collapsed.missingCount() != 0 ) ) {
            return fail( "collapsed-away value-less tips must not count as missing, got "
                    + collapsed.missingCount() + " of " + collapsed.visibleTipCount() );
        }
        return true;
    }

    // ---- value-color IDENTITY memory (JS parity): colors survive view changes, new values extend ----
    private static boolean testColorIdentityMemory() {
        final String ref = "repseq:host";
        final Map<String, Color> memory = new HashMap<String, Color>();
        final int[] next = new int[ 1 ];
        // view 1 (the "launch" view): frequency-ordered assignment, remembered
        final PropertyColorScheme v1 = new PropertyColorScheme(
                treeWith( ref, "cat", "cat", "cat", "dog", "dog", "bird" ), ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        final Color c_cat = v1.getValueColors().get( "cat" );
        final Color c_dog = v1.getValueColors().get( "dog" );
        final Color c_bird = v1.getValueColors().get( "bird" );
        if ( ( c_cat == null ) || c_cat.equals( c_dog ) || c_dog.equals( c_bird ) || c_cat.equals( c_bird ) ) {
            return fail( "launch view should assign three distinct colors" );
        }
        // view 2 = a "subtree" where the frequencies FLIP (bird now beats dog, cat gone): without memory the
        // frequency-sorted palette re-spreads; with it, every surviving value keeps its color
        final Phylogeny sub = treeWith( ref, "dog", "bird", "bird" );
        final PropertyColorScheme legacy = new PropertyColorScheme( sub, ref ); // no memory = old behavior
        if ( !c_cat.equals( legacy.getValueColors().get( "bird" ) ) ) {
            return fail( "precondition lost its teeth: the legacy re-spread should hand bird cat's old color" );
        }
        final PropertyColorScheme v2 = new PropertyColorScheme( sub, ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        if ( !c_dog.equals( v2.getValueColors().get( "dog" ) )
                || !c_bird.equals( v2.getValueColors().get( "bird" ) ) ) {
            return fail( "a value must keep its color across a view change (dog/bird re-spread)" );
        }
        // view 3: a value met for the FIRST time takes the next free slot -- never a color already handed out
        final PropertyColorScheme v3 = new PropertyColorScheme( treeWith( ref, "fish", "dog" ), ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        final Color c_fish = v3.getValueColors().get( "fish" );
        if ( c_fish.equals( c_cat ) || c_fish.equals( c_dog ) || c_fish.equals( c_bird ) ) {
            return fail( "a new value must take a FREE palette slot, not collide with a remembered one" );
        }
        final PropertyColorScheme v3b = new PropertyColorScheme( treeWith( ref, "fish", "dog" ), ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        if ( !c_fish.equals( v3b.getValueColors().get( "fish" ) ) ) {
            return fail( "a newly-met value must be REMEMBERED too" );
        }
        // an override wins over the memory but is never stored in it
        final Map<String, Color> ov = new HashMap<String, Color>();
        ov.put( "dog", Color.MAGENTA );
        final PropertyColorScheme with_ov = new PropertyColorScheme( sub, ref, ov,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        if ( !Color.MAGENTA.equals( with_ov.getValueColors().get( "dog" ) ) ) {
            return fail( "an override must win over the identity memory" );
        }
        final PropertyColorScheme after_ov = new PropertyColorScheme( sub, ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        if ( !c_dog.equals( after_ov.getValueColors().get( "dog" ) ) ) {
            return fail( "clearing an override must return the REMEMBERED automatic color, not the override" );
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
