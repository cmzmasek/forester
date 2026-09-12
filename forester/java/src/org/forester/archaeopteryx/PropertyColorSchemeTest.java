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
import org.forester.phylogeny.data.Sequence;

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
                && testColorIdentityMemory() && testSynonymDictionary() && testElementSlots()
                && testCandidateOrdering() && testModeBands() && testExcludedRefs()
                && testUniquenessTiersAndOpening() && testNumericGrammarAndFolding() && testMultiValuedAndAppliesTo()
                && testViewSummarizesNotReclassifies() && testEditKeepsValuedChoice() && testFoldToEmptyIsNoValue()
                && testContractFixtures();
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

    // ---- colorableRefs: an all-distinct NUMERIC column is offered; an all-distinct CATEGORICAL one is refused ----
    private static boolean testColorableRefs() {
        // a measurement is naturally one value per sample, so a numeric field is never refused for uniqueness
        if ( !PropertyColorScheme.colorableRefs( treeWith( "data:year", "2015", "2016", "2020", "2024" ) )
                .contains( "data:year" ) ) {
            return fail( "an all-distinct numeric column must be offered" );
        }
        // JS-authoritative (Christian, 2026-09-12): a CATEGORICAL column with as many distinct values as the tree
        // has tips is an identifier, and is REFUSED
        if ( PropertyColorScheme.colorableRefs( treeWith( "data:strain", "a", "b", "c", "d" ) )
                .contains( "data:strain" ) ) {
            return fail( "an all-distinct categorical column must be refused" );
        }
        // ...counted against ALL tips, not the covered ones: unique over the 3 of 4 tips that carry it is offered.
        // Deliberate, and not to be "improved" without asking (see the classifier's comment).
        if ( !PropertyColorScheme.colorableRefs( treeWith( "data:strain", "a", "b", "c", null ) )
                .contains( "data:strain" ) ) {
            return fail( "distinct < TOTAL tips must be offered even when every covered tip is unique" );
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
            final String name = String.format( "X%02d", i );
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
            final String expected = String.format( "X%02d", k + 6 ); // x06..x25, in alphabetical order
            if ( !expected.equals( keys.get( k ) ) ) {
                return fail( "legend entry " + k + " expected " + expected + " got " + keys.get( k ) );
            }
        }
        // by_count=true display order: same 20 kept, but most-frequent first (x25, x24, ..., x06)
        final List<String> by_count = new ArrayList<String>(
                TreePanel.orderLegendEntries( s.getValueColors(), s.getValueCounts(), 20, true ).keySet() );
        for( int k = 0; k < 20; ++k ) {
            final String expected = String.format( "X%02d", 25 - k ); // x25 down to x06
            if ( !expected.equals( by_count.get( k ) ) ) {
                return fail( "by-count legend entry " + k + " expected " + expected + " got " + by_count.get( k ) );
            }
        }
        // the 24 most frequent values (x02..x25) must all have distinct colors (no palette cycling)
        final Set<Color> colors = new HashSet<Color>();
        for( int k = 2; k <= 25; ++k ) {
            colors.add( colorForValue( s, phy, ref, String.format( "X%02d", k ) ) );
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
            final String name = String.format( "X%02d", k );
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
        // dictionary-NEUTRAL values (fox/wolf are not synonyms of anything), so this tests the pure
        // spelling fold in isolation: "Fox" x3, "fox" x1, "red_fox" x1, "Red fox" x2, "wolf" x1
        final Phylogeny phy = treeWith( ref, "Fox", "Fox", "Fox", "fox", "red_fox", "Red fox",
                                        "Red fox", "wolf" );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        if ( s.isGradient() ) {
            return fail( "host should not be a gradient" );
        }
        // three groups: {Fox, fox}, {red_fox, Red fox}, {wolf}
        if ( s.getValueColors().size() != 3 ) {
            return fail( "host groups expected 3, got " + s.getValueColors().size() );
        }
        // trivial variants -> same color
        if ( !sameColor( s, phy, ref, "Fox", "fox" ) ) {
            return fail( "Fox/fox should share a color" );
        }
        if ( !sameColor( s, phy, ref, "red_fox", "Red fox" ) ) {
            return fail( "red_fox/Red fox should share a color" );
        }
        // different values -> NOT merged (and never by substring: "Fox" is not "Red fox")
        if ( sameColor( s, phy, ref, "Fox", "wolf" ) ) {
            return fail( "Fox/wolf must NOT be merged" );
        }
        if ( sameColor( s, phy, ref, "Fox", "red_fox" ) ) {
            return fail( "Fox/red_fox must NOT be merged (no substring matching)" );
        }
        // legend shows the most frequent spelling per group, first character uppercased
        final Map<String, Color> legend = s.getValueColors();
        if ( !legend.containsKey( "Fox" ) || legend.containsKey( "fox" ) ) {
            return fail( "legend should show 'Fox' (most frequent), not 'fox'" );
        }
        if ( !legend.containsKey( "Red fox" ) || legend.containsKey( "red_fox" ) ) {
            return fail( "legend should show 'Red fox' (most frequent), not 'red_fox'" );
        }
        if ( !legend.containsKey( "Wolf" ) || legend.containsKey( "wolf" ) ) {
            return fail( "a lowercase-only group's legend label gets its first char uppercased: 'Wolf'" );
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
        if ( !"cat".equals( s.getValueKeys().get( "Cat" ) ) ) {
            return fail( "getValueKeys should map 'cat' -> 'cat'" );
        }
        // without overrides, cat is automatic again
        if ( new Color( 0x123456 ).equals( colorForValue( new PropertyColorScheme( phy, ref, null ), phy, ref, "cat" ) ) ) {
            return fail( "without overrides, cat should use an automatic color" );
        }
        return true;
    }

    // ---- the synonym dictionary folds human/Homo sapiens/h. sapiens into the canonical "Human" ----
    private static boolean testHumanSynonym() {
        final String ref = "repseq:host";
        final Phylogeny phy = treeWith( ref, "Human", "human", "Homo sapiens", "homo_sapiens", "man" );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        // two groups: the merged { Human, human, Homo sapiens, homo_sapiens } and { man }
        if ( s.getValueColors().size() != 2 ) {
            return fail( "human synonym: expected 2 groups, got " + s.getValueColors().size() );
        }
        if ( !sameColor( s, phy, ref, "Human", "Homo sapiens" ) || !sameColor( s, phy, ref, "human", "homo_sapiens" ) ) {
            return fail( "human/Homo sapiens variants should fold into one group" );
        }
        if ( sameColor( s, phy, ref, "Human", "man" ) ) {
            return fail( "'man' must NOT fold into Human (whole-value dictionary, no fuzziness)" );
        }
        // the legend shows the CANONICAL COMMON NAME "Human" (the JS-parity convention), never the variants
        final Map<String, Color> legend = s.getValueColors();
        if ( !legend.containsKey( "Human" ) || legend.containsKey( "Homo sapiens" ) || legend.containsKey( "human" ) ) {
            return fail( "legend should show the canonical 'Human', got " + legend.keySet() );
        }
        // the merged group counts all four folded/variant leaves
        final Integer count = s.getValueCounts().get( "Human" );
        if ( ( count == null ) || ( count.intValue() != 4 ) ) {
            return fail( "Human count should be 4, got " + count );
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
        // three groups, shown under their canonical common names: {Human x4}, {Mouse}, {Chicken}
        if ( s.getValueColors().size() != 3 ) {
            return fail( "host groups expected 3 (Human, Mouse, Chicken), got "
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
        if ( !s.getValueColors().containsKey( "Human" ) || !s.getValueColors().containsKey( "Mouse" )
                || !s.getValueColors().containsKey( "Chicken" ) ) {
            return fail( "host legend should show the canonical names, got " + s.getValueColors().keySet() );
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
        if ( !legend.containsKey( "Cat" ) || !legend.containsKey( "Dog" ) ) {
            return fail( "visible leaves (Cat, Dog) should remain in the legend" );
        }
        if ( legend.containsKey( "Fish" ) || legend.containsKey( "Bird" ) ) {
            return fail( "collapsed-away leaves (Fish, Bird) should be gone from the legend" );
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
        final Color c_cat = v1.getValueColors().get( "Cat" );
        final Color c_dog = v1.getValueColors().get( "Dog" );
        final Color c_bird = v1.getValueColors().get( "Bird" );
        if ( ( c_cat == null ) || c_cat.equals( c_dog ) || c_dog.equals( c_bird ) || c_cat.equals( c_bird ) ) {
            return fail( "launch view should assign three distinct colors" );
        }
        // view 2 = a "subtree" where the frequencies FLIP (bird now beats dog, cat gone): without memory the
        // frequency-sorted palette re-spreads; with it, every surviving value keeps its color
        final Phylogeny sub = treeWith( ref, "dog", "bird", "bird" );
        final PropertyColorScheme legacy = new PropertyColorScheme( sub, ref ); // no memory = old behavior
        if ( !c_cat.equals( legacy.getValueColors().get( "Bird" ) ) ) {
            return fail( "precondition lost its teeth: the legacy re-spread should hand bird cat's old color" );
        }
        final PropertyColorScheme v2 = new PropertyColorScheme( sub, ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        if ( !c_dog.equals( v2.getValueColors().get( "Dog" ) )
                || !c_bird.equals( v2.getValueColors().get( "Bird" ) ) ) {
            return fail( "a value must keep its color across a view change (dog/bird re-spread)" );
        }
        // view 3: a value met for the FIRST time takes the next free slot -- never a color already handed out
        final PropertyColorScheme v3 = new PropertyColorScheme( treeWith( ref, "fish", "dog" ), ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        final Color c_fish = v3.getValueColors().get( "Fish" );
        if ( c_fish.equals( c_cat ) || c_fish.equals( c_dog ) || c_fish.equals( c_bird ) ) {
            return fail( "a new value must take a FREE palette slot, not collide with a remembered one" );
        }
        final PropertyColorScheme v3b = new PropertyColorScheme( treeWith( ref, "fish", "dog" ), ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        if ( !c_fish.equals( v3b.getValueColors().get( "Fish" ) ) ) {
            return fail( "a newly-met value must be REMEMBERED too" );
        }
        // an override wins over the memory but is never stored in it
        final Map<String, Color> ov = new HashMap<String, Color>();
        ov.put( "dog", Color.MAGENTA );
        final PropertyColorScheme with_ov = new PropertyColorScheme( sub, ref, ov,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        if ( !Color.MAGENTA.equals( with_ov.getValueColors().get( "Dog" ) ) ) {
            return fail( "an override must win over the identity memory" );
        }
        final PropertyColorScheme after_ov = new PropertyColorScheme( sub, ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, memory, next );
        if ( !c_dog.equals( after_ov.getValueColors().get( "Dog" ) ) ) {
            return fail( "clearing an override must return the REMEMBERED automatic color, not the override" );
        }
        return true;
    }

    // ---- the synonym dictionary + parenthesis repair: the CROSS-IMPLEMENTATION CONTRACT with
    //      Archaeopteryx.js (VIS_SYNONYMS in forester.js) -- cases lifted from its spec ----
    private static boolean testSynonymDictionary() {
        final String ref = "repseq:host";
        // dictionary folds, across naming styles: common name, scientific, abbreviated scientific
        final Phylogeny phy = treeWith( ref, "swine", "Sus scrofa", "porcine", "equine", "Bos taurus (cattle)",
                                        "broiler chicken", "gallus_gallus", "ferret badger", "42-day-old pig" );
        final PropertyColorScheme s = new PropertyColorScheme( phy, ref );
        if ( !sameColor( s, phy, ref, "swine", "Sus scrofa" ) || !sameColor( s, phy, ref, "swine", "porcine" ) ) {
            return fail( "swine/Sus scrofa/porcine should all fold to Pig" );
        }
        final Map<String, Color> legend = s.getValueColors();
        if ( !legend.containsKey( "Pig" ) || !legend.containsKey( "Horse" ) || !legend.containsKey( "Chicken" ) ) {
            return fail( "expected canonical Pig/Horse/Chicken rows, got " + legend.keySet() );
        }
        // a miss with a trailing parenthetical is retried once with it stripped: "Bos taurus (cattle)" -> Cow
        if ( !legend.containsKey( "Cow" ) ) {
            return fail( "'Bos taurus (cattle)' should fold to Cow via the trailing-parenthetical retry" );
        }
        // WHOLE-VALUE only, never substring: these keep their own groups
        if ( sameColor( s, phy, ref, "ferret badger", "42-day-old pig" )
                || legend.containsKey( "Ferret" ) && !legend.containsKey( "Ferret badger" ) ) {
            return fail( "'ferret badger' must NOT fold to Ferret (no substring matching)" );
        }
        if ( !legend.containsKey( "Ferret badger" ) || !legend.containsKey( "42-day-old pig" ) ) {
            return fail( "non-dictionary values keep their own (uppercased) rows, got " + legend.keySet() );
        }
        // the qualifier cut can land inside a parenthetical; repair trims back to before the unmatched '('
        final Phylogeny phy2 = treeWith( ref, "Saimiri boliviensis (squirrel monkey; voucher: SBB04)",
                                         "Saimiri boliviensis" );
        final PropertyColorScheme s2 = new PropertyColorScheme( phy2, ref );
        if ( s2.getValueColors().size() != 1 ) {
            return fail( "parenthesis repair should merge the vouchered spelling: got "
                    + s2.getValueColors().keySet() );
        }
        if ( !s2.getValueColors().containsKey( "Saimiri boliviensis" ) ) {
            return fail( "repaired label should read 'Saimiri boliviensis', got " + s2.getValueColors().keySet() );
        }
        // a dictionary MISS keeps its trailing parenthetical in the DISPLAY (only the lookup strips it)
        final PropertyColorScheme s3 = new PropertyColorScheme( treeWith( ref, "red fox (wild)" ), ref );
        if ( !s3.getValueColors().containsKey( "Red fox (wild)" ) ) {
            return fail( "a dictionary miss keeps its parenthetical in the display, got "
                    + s3.getValueColors().keySet() );
        }
        // a value that folds to NOTHING ("_") forms no group -- and counts as missing (it draws no mark)
        final PropertyColorScheme s4 = new PropertyColorScheme( treeWith( ref, "_", "cat" ), ref );
        if ( ( s4.getValueColors().size() != 1 ) || ( s4.missingCount() != 1 ) ) {
            return fail( "a fold-to-empty value must form no group and count as missing: groups "
                    + s4.getValueColors().keySet() + ", missing " + s4.missingCount() );
        }
        return true;
    }

    // ---- element slots (JS parity): taxonomy/sequence fields color-able under the JS's reserved ids,
    //      values VERBATIM (no dictionary/fold), per-tip-unique slots refused like properties ----
    private static boolean testElementSlots() {
        // four tips: codes repeat (2x MOUSE, 2x CHICK); common names would be dictionary bait; sequence
        // names are per-tip-unique (identifier-like -> must NOT be offered)
        final Phylogeny phy = treeOf( internal( taxLeaf( "a", "MOUSE", "Mus musculus", "swine", "seqA" ),
                                                taxLeaf( "b", "MOUSE", "Mus musculus", "swine", "seqB" ) ),
                                      internal( taxLeaf( "c", "CHICK", "Gallus gallus", "chicken", "seqC" ),
                                                taxLeaf( "d", "CHICK", "Gallus gallus", "chicken", "seqD" ) ) );
        final PhylogenyNode any = phy.getFirstExternalNode();
        if ( !"MOUSE".equals( PropertyColorScheme.valueFor( any, PropertyColorScheme.TAX_CODE_REF ) )
                || !"Mus musculus".equals( PropertyColorScheme.valueFor( any, PropertyColorScheme.TAX_SCI_NAME_REF ) )
                || !"swine".equals( PropertyColorScheme.valueFor( any, PropertyColorScheme.TAX_COMMON_NAME_REF ) )
                || !"seqA".equals( PropertyColorScheme.valueFor( any, PropertyColorScheme.SEQ_NAME_REF ) ) ) {
            return fail( "element-slot valueFor should read the taxonomy/sequence fields" );
        }
        if ( PropertyColorScheme.valueFor( any, PropertyColorScheme.SEQ_SYMBOL_REF ) != null ) {
            return fail( "an absent sequence symbol should read null" );
        }
        final List<String> refs = PropertyColorScheme.colorableRefs( phy );
        if ( !refs.contains( PropertyColorScheme.TAX_CODE_REF )
                || !refs.contains( PropertyColorScheme.TAX_SCI_NAME_REF )
                || !refs.contains( PropertyColorScheme.TAX_COMMON_NAME_REF ) ) {
            return fail( "repeating taxonomy slots should be offered, got " + refs );
        }
        // a per-tip-unique sequence-name slot is an identifier: REFUSED, exactly as a property would be
        if ( refs.contains( PropertyColorScheme.SEQ_NAME_REF ) ) {
            return fail( "a per-tip-unique sequence-name slot must be refused, got " + refs );
        }
        // element values are used VERBATIM: the common name "swine" keeps its own row -- it must NOT fold
        // to the dictionary's "Pig" when the field being colored IS the taxonomy common name
        final PropertyColorScheme s = new PropertyColorScheme( phy, PropertyColorScheme.TAX_COMMON_NAME_REF );
        if ( !s.getValueColors().containsKey( "swine" ) || s.getValueColors().containsKey( "Pig" ) ) {
            return fail( "element-slot values must stay verbatim (no dictionary), got "
                    + s.getValueColors().keySet() );
        }
        if ( s.missingCount() != 0 ) {
            return fail( "all four tips carry a common name; missing should be 0, got " + s.missingCount() );
        }
        // ...and the JS display labels
        if ( !"Taxonomy Code".equals( PropertyColorScheme.displayName( PropertyColorScheme.TAX_CODE_REF ) )
                || !"Gene Name".equals( PropertyColorScheme.displayName( PropertyColorScheme.SEQ_GENE_NAME_REF ) ) ) {
            return fail( "element slots should carry the JS Color-menu display labels" );
        }
        return true;
    }

    /** A leaf with taxonomy (code/sci/common) and a named sequence. */
    private static PhylogenyNode taxLeaf( final String name, final String code, final String sci,
                                          final String common, final String seq_name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final org.forester.phylogeny.data.Taxonomy t = new org.forester.phylogeny.data.Taxonomy();
        try {
            t.setTaxonomyCode( code );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        t.setScientificName( sci );
        t.setCommonName( common );
        n.getNodeData().setTaxonomy( t );
        final org.forester.phylogeny.data.Sequence q = new org.forester.phylogeny.data.Sequence();
        q.setName( seq_name );
        n.getNodeData().setSequence( q );
        return n;
    }

    // ---- candidate ordering (JS parity): tiers, entropy score, sparse ranked LAST not refused ----
    private static boolean testCandidateOrdering() {
        // 12 tips; four fields of known character:
        //   data:balanced  categorical, full coverage, 3 even groups   -> tier 0, HIGH score
        //   data:skewed    categorical, full coverage, 11-vs-1 split   -> tier 0, LOW score
        //   data:years     numeric, full coverage, 12 distinct         -> tier 1 (gradient default)
        //   data:sparse    categorical, 3 of 12 tips                   -> tier 4 (sparse, ranked last)
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 12; ++i ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "t" + i );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:balanced", "g" + ( i % 3 ), "", "xsd:string", AppliesTo.NODE ) );
            pl.addProperty( new Property( "data:skewed", ( i == 0 ) ? "rare" : "common", "", "xsd:string",
                                          AppliesTo.NODE ) );
            pl.addProperty( new Property( "data:years", Integer.toString( 2000 + i ), "", "xsd:decimal",
                                          AppliesTo.NODE ) );
            if ( i < 3 ) {
                pl.addProperty( new Property( "data:sparse", "s" + i % 2, "", "xsd:string", AppliesTo.NODE ) );
            }
            n.getNodeData().setProperties( pl );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        final List<String> refs = PropertyColorScheme.colorableRefs( phy );
        final int i_bal = refs.indexOf( "data:balanced" );
        final int i_skew = refs.indexOf( "data:skewed" );
        final int i_years = refs.indexOf( "data:years" );
        final int i_sparse = refs.indexOf( "data:sparse" );
        if ( ( i_bal < 0 ) || ( i_skew < 0 ) || ( i_years < 0 ) || ( i_sparse < 0 ) ) {
            return fail( "all four fields should be OFFERED (sparse ranked, not refused): " + refs );
        }
        if ( !( ( i_bal < i_skew ) && ( i_skew < i_years ) && ( i_years < i_sparse ) ) ) {
            return fail( "expected balanced < skewed < years < sparse, got " + refs );
        }
        // the tree opens with the first candidate that is not wide...
        if ( !"data:balanced".equals( PropertyColorScheme.autoColorCandidate( phy ) ) ) {
            return fail( "auto-color should pick the top candidate, got "
                    + PropertyColorScheme.autoColorCandidate( phy ) );
        }
        // ...so a sparse field DOES open a tree when nothing above it exists. JS-authoritative 2026-09-12: the old
        // "only tiers 0 and 1 open a tree" is retired.
        final Phylogeny sparse_only = new Phylogeny();
        final PhylogenyNode r2 = new PhylogenyNode();
        for( int i = 0; i < 12; ++i ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "u" + i );
            if ( i < 3 ) {
                final PropertiesList pl = new PropertiesList();
                pl.addProperty( new Property( "data:sparse", "s" + i % 2, "", "xsd:string", AppliesTo.NODE ) );
                n.getNodeData().setProperties( pl );
            }
            r2.addAsChild( n );
        }
        sparse_only.setRoot( r2 );
        sparse_only.externalNodesHaveChanged();
        if ( !"data:sparse".equals( PropertyColorScheme.autoColorCandidate( sparse_only ) ) ) {
            return fail( "a sparse field that is the only candidate must open the tree; got "
                    + PropertyColorScheme.autoColorCandidate( sparse_only ) );
        }
        return true;
    }

    // ---- the numeric three-band system + the forced-mode constructor ----
    private static boolean testModeBands() {
        final String ref = "repseq:n";
        final Phylogeny five = numericTree( ref, 5, 3 );   // 5 distinct numbers over 15 tips
        final Phylogeny fifteen = numericTree( ref, 15, 1 );
        final Phylogeny thirty = numericTree( ref, 30, 1 );
        if ( PropertyColorScheme.colorModeBand( five.getExternalNodes(), ref )
                != PropertyColorScheme.ModeBand.SMALL ) {
            return fail( "5 distinct numbers should band SMALL" );
        }
        if ( PropertyColorScheme.colorModeBand( fifteen.getExternalNodes(), ref )
                != PropertyColorScheme.ModeBand.MEDIUM ) {
            return fail( "15 distinct numbers should band MEDIUM" );
        }
        if ( PropertyColorScheme.colorModeBand( thirty.getExternalNodes(), ref )
                != PropertyColorScheme.ModeBand.LARGE ) {
            return fail( "30 distinct numbers should band LARGE" );
        }
        if ( PropertyColorScheme.colorModeBand( treeWith( ref, "cat", "dog" ).getExternalNodes(), ref )
                != PropertyColorScheme.ModeBand.NOT_NUMERIC ) {
            return fail( "a text field should band NOT_NUMERIC" );
        }
        if ( !PropertyColorScheme.ModeBand.SMALL.isSwitchable()
                || !PropertyColorScheme.ModeBand.MEDIUM.isSwitchable()
                || PropertyColorScheme.ModeBand.LARGE.isSwitchable()
                || PropertyColorScheme.ModeBand.SMALL.defaultsToGradient()
                || !PropertyColorScheme.ModeBand.MEDIUM.defaultsToGradient() ) {
            return fail( "band switchability/defaults wrong" );
        }
        // the forced-mode constructor: the same numeric data as COLORS or as a GRADIENT
        final PropertyColorScheme colors = new PropertyColorScheme( five, ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, null, null, Boolean.FALSE );
        if ( colors.isGradient() || ( colors.getValueColors().size() != 5 ) ) {
            return fail( "forced FALSE should color 5 distinct numbers individually" );
        }
        final PropertyColorScheme grad = new PropertyColorScheme( five, ref, null,
                PropertyColorScheme.DEFAULT_PALETTE_NAME, null, null, Boolean.TRUE );
        if ( !grad.isGradient() ) {
            return fail( "forced TRUE should be a gradient" );
        }
        // the UN-forced default is unchanged (plain numeric detection -- annotation columns rely on it)
        if ( !new PropertyColorScheme( five, ref ).isGradient() ) {
            return fail( "the legacy (null-forced) default must stay the plain numeric detection" );
        }
        return true;
    }

    /** A tree with {@code distinct * per} tips carrying numeric values 0..distinct-1, each {@code per} times. */
    private static Phylogeny numericTree( final String ref, final int distinct, final int per ) {
        final List<String> vals = new ArrayList<String>();
        for( int v = 0; v < distinct; ++v ) {
            for( int k = 0; k < per; ++k ) {
                vals.add( Integer.toString( v ) );
            }
        }
        return treeWith( ref, vals.toArray( new String[ 0 ] ) );
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

    /**
     * Refs that must never be offered for visualization, and the displayed name the rules match against.
     * CROSS-IMPLEMENTATION CONTRACT with Archaeopteryx.js, agreed 2026-09-11 -- neither side adds an exclusion
     * or changes the prettifier alone.
     */
    private static boolean testExcludedRefs() {
        try {
            // --- JS-authoritative (forester.js visExcludedRef): refs that describe the RECORD, not the organism ---
            // Measured on a 13,246-tip BV-BRC tree: all three of these were offered in Color-by before this rule.
            final String[][] excluded = { { "vipr:BVBRC_Accession", "accession" }, { "vipr:GB_Accession", "accession" },
                    { "vipr:NCBI_Taxon_Id", "taxon id" }, { "ncbi_taxid", "taxon id" }, { "taxonomy_id", "taxon id" },
                    { "style:branch_color", "style" }, { "x:Abbr_Authors", "author" }, { "x:Region Set", "set" },
                    { "x:Data-Use_Terms", "data use" }, { "dataUseTerms", "data use" },
                    { "x:strain_identifier", "identifier" }, { "x:accessions", "accession" },
                    // no separator and no lowercase->uppercase boundary: the displayed name is ONE word, and the
                    // accession / identifier rules are SUFFIX rules, so they still catch it
                    { "x:GBAccession", "accession" }, { "x:GBIdentifier", "identifier" },
                    { "x:XYZaccession", "accession" }, { "x:XYZidentifier", "identifier" },
                    // the WORD "id" (JS e66ccff): BV-BRC's genome_id / patric_id pass every statistical test
                    { "BVBRC:genome_id", "id" }, { "BVBRC:patric_id", "id" }, { "x:GenomeID", "id" },
                    { "x:genomeId", "id" }, { "x:Feature_ID", "id" }, { "x:genome_ids", "id" }, { "x:id", "id" },
                    // "restricted until" in every spelling (JS f4fa9e3): a data-use embargo date. The space is
                    // optional so the one-word "restricteduntil", which prettifying cannot split, is caught too
                    { "x:Restricted Until", "restricted until" }, { "x:restricted_until", "restricted until" },
                    { "x:restricted-until", "restricted until" }, { "x:restrictedUntil", "restricted until" },
                    { "x:RESTRICTED_UNTIL", "restricted until" }, { "x:restricteduntil", "restricted until" },
                    { "BVBRC:restricted_until", "restricted until" } };
            for( final String[] e : excluded ) {
                if ( !e[ 1 ].equals( PropertyColorScheme.excludedRefReason( e[ 0 ] ) ) ) {
                    System.out.println( "[PropertyColorSchemeTest] " + e[ 0 ] + " must be excluded as '" + e[ 1 ]
                            + "'; got " + PropertyColorScheme.excludedRefReason( e[ 0 ] ) );
                    return false;
                }
            }
            // ...and the controls that keep the word matching honest: these must SURVIVE
            for( final String keep : new String[] { "x:Dataset", "x:Subset", "x:Authority", "x:Within_Group",
                    "vipr:Host", "repseq:protein_names", "data:year",
                    // the accession rule matches a SUFFIX, so a trailing qualifier keeps the field: these read as
                    // "which number/column", not "the record's accession"
                    "x:Accession Number", "x:accession_number", "x:AccessionX",
                    // "id" is a WORD rule: two letters end a great many ordinary words
                    "x:Plasmid", "x:Hybrid", "x:Lipid", "x:Nucleic Acid", "x:Idea",
                    // "restricted until" is word-anchored too
                    "x:Restricted", "x:Unrestricted", "x:Restriction",
                    // strain and genome name are NOT excluded by name (Christian, 2026-09-12): they group sequences
                    // wherever an isolate contributes several
                    "x:strain", "BVBRC:genome_name" } ) {
                if ( PropertyColorScheme.isExcludedRef( keep ) ) {
                    System.out.println( "[PropertyColorSchemeTest] " + keep + " must NOT be excluded; got "
                            + PropertyColorScheme.excludedRefReason( keep ) );
                    return false;
                }
            }
            // ...and the rules must actually be APPLIED to candidacy, not merely exist. A sabotage that removed the
            // filter from colorableRefs PASSED the assertions above, which is what caught this gap: they pinned the
            // predicate and nothing pinned its use.
            {
                final Phylogeny p = Phylogeny.createInstanceFromNhxString( "((a,b),(c,d));" );
                int i = 0;
                for( final PhylogenyNode n : p.getExternalNodes() ) {
                    ++i;
                    final PropertiesList pl = new PropertiesList();
                    pl.addProperty( new Property( "vipr:NCBI_Taxon_Id", String.valueOf( 100 + ( i % 2 ) ), "",
                                                  "xsd:string", AppliesTo.NODE ) );
                    pl.addProperty( new Property( "vipr:Host", ( ( i % 2 ) == 0 ) ? "human" : "avian", "",
                                                  "xsd:string", AppliesTo.NODE ) );
                    n.getNodeData().setProperties( pl );
                }
                final List<String> refs = PropertyColorScheme.colorableRefs( p );
                if ( refs.contains( "vipr:NCBI_Taxon_Id" ) ) {
                    System.out.println( "[PropertyColorSchemeTest] an excluded ref must not be a colorable candidate" );
                    return false;
                }
                if ( !refs.contains( "vipr:Host" ) ) {
                    System.out.println( "[PropertyColorSchemeTest] an ordinary ref must still be offered; got " + refs );
                    return false;
                }
            }
            // --- the displayed name, which the exclusion rules above match against ---
            final String[][] pretty = { { "geographic_group", "Geographic Group" },
                    { "dataUseTerms", "Data Use Terms" }, { "FluSeason", "Flu Season" },
                    { "NCBI_Taxon_Id", "NCBI Taxon Id" }, { "GlobalH1Clade", "Global H1 Clade" },
                    { "host", "Host" }, { "lowerUPPER", "Lower UPPER" }, { "PANGO", "PANGO" },
                    // H5N1 must NOT become "H5 N1": a digit followed by an uppercase is not a word break, or every
                    // flu subtype in the data is mangled
                    { "H5N1", "H5N1" }, { "HA", "HA" },
                    // A HYPHEN IS NOT A WORD SEPARATOR -- it is preserved, and a hyphenated word counts as ONE
                    // word for capitalisation. Diffed byte-for-byte against the Archaeopteryx.js contract
                    // fixture (141/141), so neither side may reintroduce a hyphen-to-space normalisation alone.
                    { "In-Group", "In-Group" }, { "Out-Groups", "Out-Groups" }, { "IN-GROUP", "IN-GROUP" },
                    { "strain-identifier", "Strain-identifier" }, { "Data-Use Terms", "Data-Use Terms" },
                    { "Within-Group", "Within-Group" },
                    // forester.js prettifyVisLabel VERBATIM: a word is capitalised when its FIRST character is a
                    // lowercase letter, whatever follows (the desktop used to leave any word carrying a capital
                    // alone, so "h5N1" stayed "h5N1"); runs of spaces are kept; the namespace ends at the FIRST colon
                    { "h5N1", "H5N1" }, { "iPhone", "I Phone" }, { "two  spaces", "Two  Spaces" },
                    { "a:b:c", "B:c" } };
            for( final String[] p : pretty ) {
                if ( !p[ 1 ].equals( PropertyColorScheme.displayName( p[ 0 ] ) ) ) {
                    System.out.println( "[PropertyColorSchemeTest] displayName(" + p[ 0 ] + ") must be '" + p[ 1 ]
                            + "'; got '" + PropertyColorScheme.displayName( p[ 0 ] ) + "'" );
                    return false;
                }
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    /**
     * The three-way uniqueness rule, the tiers and what OPENS a tree -- JS-authoritative (Christian, 2026-09-12;
     * forester.js 79dd9de). A categorical field with as many distinct values as the tree has TIPS is refused; above
     * 20 distinct values it is WIDE (tier 2), and also NEAR-UNIQUE (tier 5) when distinct/covered > 3/5; a numeric
     * field is never refused for being unique and is always tier 1. The tree opens with the FIRST CANDIDATE THAT IS
     * NOT WIDE, so tiers 0, 1, 3 and 4 can open one and 2 and 5 never do. This replaces the desktop's 2026-09-11 rule,
     * which ranked an all-distinct categorical instead of refusing it.
     */
    private static boolean testUniquenessTiersAndOpening() {
        try {
            // fully annotated, per-tip-unique categorical: an identifier, REFUSED
            if ( cand( propTree( 30, 30, true ), "x:field" ) != null ) {
                return fail( "a fully annotated all-distinct categorical must be refused" );
            }
            // the FluSeason shape: 2 of 6 tips, 2 values. distinct (2) < TOTAL (6), so offered -- sparse, tier 4 --
            // and, as the only candidate, it OPENS the tree
            final Phylogeny sample_of_two = propTree( 6, 2, true );
            final PropertyColorScheme.VisCandidate fs = cand( sample_of_two, "x:field" );
            if ( ( fs == null ) || !fs._sparse || ( fs.tier() != 4 ) ) {
                return fail( "two of six tips with two values must be offered as sparse (tier 4)" );
            }
            if ( !"x:field".equals( PropertyColorScheme.autoColorCandidate( sample_of_two ) ) ) {
                return fail( "a sparse field that is the only candidate must open the tree" );
            }
            // coverage at the bar: 20 of 30 is exactly 2/3 and NOT sparse; 19 of 30 is
            if ( cand( fieldTree( 30, refs( "x:F" ), i -> ( i < 20 ) ? ( ( ( i % 2 ) == 1 ) ? "A" : "B" ) : null ),
                       "x:F" )._sparse
                    || !cand( fieldTree( 30, refs( "x:F" ),
                                         i -> ( i < 19 ) ? ( ( ( i % 2 ) == 1 ) ? "A" : "B" ) : null ), "x:F" )._sparse ) {
                return fail( "sparse means covered/total < 2/3: 20 of 30 is not sparse, 19 of 30 is" );
            }
            // WIDE at exactly 3/5: 21 values over 35 tips -> wide, NOT near-unique, tier 2, and never opens
            final Phylogeny wide = fieldTree( 35, refs( "x:F" ), i -> "v" + ( i % 21 ) );
            final PropertyColorScheme.VisCandidate w = cand( wide, "x:F" );
            if ( ( w == null ) || !w._wide || w._near_unique || ( w.tier() != 2 ) ) {
                return fail( "21 values over 35 tips (exactly 3/5) must be wide, tier 2" );
            }
            if ( PropertyColorScheme.autoColorCandidate( wide ) != null ) {
                return fail( "a tree whose only candidate is wide must open uncoloured" );
            }
            // NEAR-UNIQUE just past 3/5: 21 values over 34 tips -> tier 5
            final PropertyColorScheme.VisCandidate nu = cand( fieldTree( 34, refs( "x:F" ), i -> "v" + ( i % 21 ) ),
                                                              "x:F" );
            if ( ( nu == null ) || !nu._wide || !nu._near_unique || ( nu.tier() != 5 ) ) {
                return fail( "21 values over 34 tips must be near-unique, tier 5" );
            }
            // exactly 20 distinct values is still a clean category
            final PropertyColorScheme.VisCandidate twenty = cand( fieldTree( 22, refs( "x:F" ),
                                                                             i -> "v" + ( ( i < 20 ) ? i : 0 ) ),
                                                                  "x:F" );
            if ( ( twenty == null ) || twenty._wide || ( twenty.tier() != 0 ) ) {
                return fail( "20 distinct values must be a clean category, tier 0" );
            }
            // OPENS PAST WIDE: a wide field above an In-Group does not stop the In-Group from opening the tree
            final Phylogeny past = fieldTree( 36, refs( "x:Wide", "x:In-Group" ), i -> "v" + ( i % 21 ),
                                              i -> ( ( i % 2 ) == 1 ) ? "in" : "out" );
            if ( !PropertyColorScheme.colorableRefs( past ).equals( java.util.Arrays.asList( "x:Wide", "x:In-Group" ) ) ) {
                return fail( "the wide field (tier 2) must precede the In-Group (tier 3); got "
                        + PropertyColorScheme.colorableRefs( past ) );
            }
            if ( !"x:In-Group".equals( PropertyColorScheme.autoColorCandidate( past ) ) ) {
                return fail( "a wide field must not block the In-Group from opening the tree; got "
                        + PropertyColorScheme.autoColorCandidate( past ) );
            }
            // a LOPSIDED category (29 cat / 1 dog) opens past a near-unique field (25 distinct over 30 tips). The
            // companion is lopsided on purpose so it scores far below: only the tier can put it first.
            final Phylogeny near = fieldTree( 30, refs( "x:wide", "x:host" ),
                                              i -> ( i < 20 ) ? ( "u" + i ) : ( "r" + ( i % 5 ) ),
                                              i -> ( i == 0 ) ? "dog" : "cat" );
            if ( !"x:host".equals( PropertyColorScheme.colorableRefs( near ).get( 0 ) )
                    || !"x:host".equals( PropertyColorScheme.autoColorCandidate( near ) ) ) {
                return fail( "a lopsided category must lead and open past a near-unique field; got "
                        + PropertyColorScheme.colorableRefs( near ) );
            }
            // NUMERIC: eight tips, eight distinct numbers -- offered, tier 1, and it opens the tree (the shape of
            // size-by-property.xml)
            final Phylogeny nums = fieldTree( 8, refs( "x:read_count" ), i -> String.valueOf( 100 + ( i * 7 ) ) );
            final PropertyColorScheme.VisCandidate rc = cand( nums, "x:read_count" );
            if ( ( rc == null ) || !rc._numeric || ( rc.tier() != 1 )
                    || !"x:read_count".equals( PropertyColorScheme.autoColorCandidate( nums ) ) ) {
                return fail( "an all-distinct numeric field must be offered at tier 1 and open the tree" );
            }
            // ...never refused however unique: 30 distinct numbers over 30 tips is a range with no switch
            final PropertyColorScheme.VisCandidate r30 = cand( fieldTree( 30, refs( "x:F" ),
                                                                          i -> String.valueOf( ( i * 7 ) + 3 ) ),
                                                               "x:F" );
            if ( ( r30 == null ) || !"range".equals( r30._color_mode ) || r30._switchable ) {
                return fail( "30 distinct numbers must be offered as a range, not switchable" );
            }
            // TIER beats SCORE: a numeric field (8 even values) ranks below a lopsided clean category (25/75)
            final Phylogeny tiers = fieldTree( 40, refs( "x:Numeric", "x:Clean" ), i -> String.valueOf( i % 8 ),
                                               i -> ( ( i % 4 ) == 0 ) ? "P" : "Q" );
            if ( !PropertyColorScheme.colorableRefs( tiers ).equals( java.util.Arrays.asList( "x:Clean", "x:Numeric" ) ) ) {
                return fail( "a clean category (tier 0) must precede a numeric field (tier 1); got "
                        + PropertyColorScheme.colorableRefs( tiers ) );
            }
            // SHAPE: seven distinct values or fewer
            if ( !cand( fieldTree( 14, refs( "x:F" ), i -> "v" + ( i % 7 ) ), "x:F" )._shape
                    || cand( fieldTree( 16, refs( "x:F" ), i -> "v" + ( i % 8 ) ), "x:F" )._shape ) {
                return fail( "shape means <= 7 distinct values" );
            }
            // LABEL COLLISION: two refs whose prettified labels coincide are BOTH shown as their full ref
            final Phylogeny coll = fieldTree( 10, refs( "a:host", "b:host" ), i -> ( ( i % 2 ) == 0 ) ? "cat" : "dog",
                                              i -> ( ( i % 3 ) == 0 ) ? "x" : "y" );
            if ( !"a:host".equals( cand( coll, "a:host" )._label ) || !"b:host".equals( cand( coll, "b:host" )._label ) ) {
                return fail( "colliding labels must both read as the full ref" );
            }
            if ( !"Host".equals( cand( fieldTree( 10, refs( "a:host" ), i -> ( ( i % 2 ) == 0 ) ? "cat" : "dog" ),
                                       "a:host" )._label ) ) {
                return fail( "an uncontested label must stay prettified" );
            }
            // DEPRIORITIZED: in/out-group is offered but ranks at tier 3. Paired with a field of the SAME distribution
            // so the scores tie and only the tier decides -- and the companion is "Zone" on purpose, since "in group"
            // sorts BEFORE "zone", so this fails if the tier is removed.
            for( final String ig : new String[] { "x:In-Group", "x:InGroup", "x:In Group", "x:in_group",
                                                  "x:ingroup", "x:IN-GROUP", "x:Out-Group", "x:OutGroup",
                                                  "x:out_group", "x:outgroup", "x:In-Groups", "x:Outgroups" } ) {
                final Phylogeny t = fieldTree( 8, refs( ig, "x:Zone" ), i -> ( i < 4 ) ? "Alpha" : "Beta",
                                               i -> ( i < 4 ) ? "Alpha" : "Beta" );
                final List<String> refs = PropertyColorScheme.colorableRefs( t );
                if ( !refs.contains( ig ) ) {
                    return fail( ig + " must still be OFFERED, only demoted" );
                }
                if ( !"x:Zone".equals( refs.get( 0 ) ) || !"x:Zone".equals( PropertyColorScheme.autoColorCandidate( t ) ) ) {
                    return fail( ig + " must not lead or open the tree; got " + refs );
                }
            }
            // ...and the word anchoring is not decoration: "Within Group" CONTAINS "in group" and must keep leading
            if ( PropertyColorScheme.isDeprioritizedRef( "x:Within-Group" )
                    || PropertyColorScheme.isDeprioritizedRef( "x:Within Group" )
                    || PropertyColorScheme.isDeprioritizedRef( "x:Ingredient" )
                    || PropertyColorScheme.isDeprioritizedRef( "x:Grouping" ) ) {
                return fail( "the in/out-group rule must be anchored as WORDS of the displayed label" );
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    /**
     * "Numeric" is a pinned GRAMMAR, not a host language's parser, and spellings of one number are ONE value --
     * JS-authoritative (decision 3, 2026-09-12). Java's parseDouble takes "Infinity" and "NaN" and JavaScript's
     * Number() takes "0x1A" and "0b101", so a field of either would be a gradient in one program and a category in
     * the other.
     */
    private static boolean testNumericGrammarAndFolding() {
        for( final String num : new String[] { "+5", ".5", "5.", "1e3", "-2E-1", "2009", "007", "-0", "1.5e+10" } ) {
            if ( !PropertyColorScheme.isVisNumber( num ) ) {
                return fail( "'" + num + "' must be a number under the grammar" );
            }
        }
        for( final String word : new String[] { "0x1A", "0b101", "Infinity", "NaN", "1,5", "1_000", "", "e3", ".",
                                                "+", "5f", "5d", "1e", "--5", " 5" } ) {
            if ( PropertyColorScheme.isVisNumber( word ) ) {
                return fail( "'" + word + "' must NOT be a number under the grammar" );
            }
        }
        // the numeric_grammar fixture tree, built directly
        final Phylogeny t = fieldTree( 20, refs( "x:Fold", "x:Hex", "x:Inf", "x:Forms" ),
                                       i -> new String[] { "1", "1.0", "2", "+2", "3" }[ i % 5 ],
                                       i -> ( ( i % 4 ) != 0 ) ? String.valueOf( i % 4 ) : "0x1A",
                                       i -> ( ( i % 4 ) != 0 ) ? String.valueOf( i % 4 ) : "Infinity",
                                       i -> new String[] { "+5", ".5", "5.", "1e3", "7" }[ i % 5 ] );
        final PropertyColorScheme.VisCandidate fold = cand( t, "x:Fold" );
        if ( ( fold == null ) || !fold._numeric || !fold._values.equals( java.util.Arrays.asList( "1", "2", "3" ) ) ) {
            return fail( "x:Fold must be numeric with 3 folded values; got " + ( ( fold == null ) ? null : fold._values ) );
        }
        if ( !Integer.valueOf( 8 ).equals( fold._counts.get( "1" ) ) || !Integer.valueOf( 8 ).equals( fold._counts.get( "2" ) )
                || !Integer.valueOf( 4 ).equals( fold._counts.get( "3" ) ) ) {
            return fail( "folded counts must add up (1: 8, 2: 8, 3: 4); got " + fold._counts );
        }
        final PropertyColorScheme.VisCandidate hex = cand( t, "x:Hex" );
        final PropertyColorScheme.VisCandidate inf = cand( t, "x:Inf" );
        if ( ( hex == null ) || hex._numeric || ( hex._values.size() != 4 ) || ( inf == null ) || inf._numeric
                || ( inf._values.size() != 4 ) ) {
            return fail( "a field containing 0x1A or Infinity must be categorical with 4 values" );
        }
        // every decimal spelling is a number; "+5" and "5." fold to the SHORTEST spelling, ties by code point
        final PropertyColorScheme.VisCandidate forms = cand( t, "x:Forms" );
        if ( ( forms == null ) || !forms._numeric
                || !forms._values.equals( java.util.Arrays.asList( ".5", "+5", "7", "1e3" ) ) ) {
            return fail( "x:Forms must fold +5/5. and sort numerically; got " + ( ( forms == null ) ? null : forms._values ) );
        }
        // a node carrying "1.0" reads "1": the same colour and the same legend row
        final PropertyColorScheme s = new PropertyColorScheme( t, "x:Fold", null, PropertyColorScheme.DEFAULT_PALETTE_NAME,
                                                               null, null, Boolean.FALSE, fold );
        if ( !s.getValueColors().keySet().equals( new java.util.HashSet<String>( java.util.Arrays.asList( "1", "2", "3" ) ) )
                || !Integer.valueOf( 8 ).equals( s.getValueCounts().get( "1" ) ) ) {
            return fail( "the scheme must show one row per number; got " + s.getValueCounts() );
        }
        if ( !sameColor( s, t, "x:Fold", "1", "1.0" ) || !sameColor( s, t, "x:Fold", "2", "+2" ) ) {
            return fail( "spellings of one number must share a colour" );
        }
        // ...and as a gradient every spelling is placed
        final PropertyColorScheme g = new PropertyColorScheme( t, "x:Fold", null, PropertyColorScheme.DEFAULT_PALETTE_NAME,
                                                               null, null, Boolean.TRUE, fold );
        if ( !g.isGradient() || !"1".equals( g.getGradientMinLabel() ) || !"3".equals( g.getGradientMaxLabel() )
                || ( colorForValue( g, t, "x:Fold", "+2" ) == null ) ) {
            return fail( "the folded field as a gradient must span 1..3 and colour every spelling" );
        }
        // folding can leave ONE value, which is then refused like any other
        if ( cand( fieldTree( 10, refs( "x:F" ), i -> ( ( i % 2 ) == 0 ) ? "1" : "1.0" ), "x:F" ) != null ) {
            return fail( "a field whose spellings fold to one number must be refused" );
        }
        // folding applies to a numeric ELEMENT SLOT too
        final Phylogeny slot = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final String[] genes = { "1", "1.0", "2", "2.0", "3", "3" };
        for( final String gene : genes ) {
            final PhylogenyNode n = new PhylogenyNode();
            final Sequence q = new Sequence();
            q.setGeneName( gene );
            n.getNodeData().setSequence( q );
            root.addAsChild( n );
        }
        slot.setRoot( root );
        slot.externalNodesHaveChanged();
        final PropertyColorScheme.VisCandidate gn = cand( slot, PropertyColorScheme.SEQ_GENE_NAME_REF );
        if ( ( gn == null ) || !gn._numeric || !gn._values.equals( java.util.Arrays.asList( "1", "2", "3" ) ) ) {
            return fail( "a numeric element slot must fold its spellings; got " + ( ( gn == null ) ? null : gn._values ) );
        }
        if ( !"1".equals( PropertyColorScheme.visualizationNodeValue( slot.getExternalNodes().get( 1 ), gn ) ) ) {
            return fail( "a slot value \"1.0\" must read as \"1\"" );
        }
        return true;
    }

    /**
     * A ref carried TWICE by one tip refuses the field (this was a DESKTOP DEFECT: it offered the field and silently
     * coloured that tip by its first value), and only applies_to "node" and "clade" are node data.
     */
    private static boolean testMultiValuedAndAppliesTo() {
        final Phylogeny multi = fieldTree( 30, refs( "x:F", "x:G" ),
                                           i -> ( i == 0 ) ? new String[] { "A", "B" } : ( ( ( i % 2 ) == 1 ) ? "A" : "B" ),
                                           i -> ( ( i % 2 ) == 1 ) ? "A" : "B" );
        if ( cand( multi, "x:F" ) != null ) {
            return fail( "a ref carried twice by one tip must refuse the whole field" );
        }
        if ( cand( multi, "x:G" ) == null ) {
            return fail( "the same distribution carried once must be offered" );
        }
        // the SAME value twice is still carried twice
        if ( cand( fieldTree( 30, refs( "x:F" ),
                              i -> ( i == 0 ) ? new String[] { "A", "A" } : ( ( ( i % 2 ) == 1 ) ? "A" : "B" ) ),
                   "x:F" ) != null ) {
            return fail( "a ref carried twice with one value must still refuse the field" );
        }
        // applies_to: node and clade are the node's own data; parent_branch is the branch above it
        final Phylogeny at = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 30; ++i ) {
            final PhylogenyNode n = new PhylogenyNode();
            final String v = ( ( i % 2 ) == 1 ) ? "A" : "B";
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "x:OnNode", v, "", "xsd:string", AppliesTo.NODE ) );
            pl.addProperty( new Property( "x:OnClade", v, "", "xsd:string", AppliesTo.CLADE ) );
            pl.addProperty( new Property( "x:OnBranch", v, "", "xsd:string", AppliesTo.PARENT_BRANCH ) );
            // the same ref as a node property AND a parent_branch property: one value, not two
            pl.addProperty( new Property( "x:Both", v, "", "xsd:string", AppliesTo.NODE ) );
            pl.addProperty( new Property( "x:Both", "Z", "", "xsd:string", AppliesTo.PARENT_BRANCH ) );
            n.getNodeData().setProperties( pl );
            root.addAsChild( n );
        }
        at.setRoot( root );
        at.externalNodesHaveChanged();
        if ( ( cand( at, "x:OnNode" ) == null ) || ( cand( at, "x:OnClade" ) == null ) || ( cand( at, "x:OnBranch" ) != null ) ) {
            return fail( "applies_to node and clade must be offered, parent_branch refused; got "
                    + PropertyColorScheme.colorableRefs( at ) );
        }
        final PropertyColorScheme.VisCandidate both = cand( at, "x:Both" );
        if ( ( both == null ) || ( both._values.size() != 2 )
                || "Z".equals( PropertyColorScheme.visualizationNodeValue( at.getExternalNodes().get( 0 ), both ) ) ) {
            return fail( "a parent_branch property must be neither a second value nor the node's value" );
        }
        return true;
    }

    /**
     * A VIEW never re-decides candidacy -- it only re-summarizes (JS-authoritative, decision 1, 2026-09-12). Measured
     * on the JS side: re-running the classifier per subtree and applying its refusals to the chosen field dropped the
     * colouring in 61% of the corpus's clades, because a clade is by nature a set of tips sharing a value and one
     * value is refused. The desktop already rebuilt its scheme over the visible tips; this pins that the scheme, the
     * band and the offered fields come from the TREE's candidate, not from classifying the view.
     */
    private static boolean testViewSummarizesNotReclassifies() {
        // two clades of six: Genus is one value per clade; Year is 12 distinct numbers on the tree and 6 in a clade;
        // Hex is categorical on the tree ("0x1A" lives in clade B only) but all digits in clade A; Host only in B
        final Phylogeny t = twoCladeTree( 6, refs( "x:Genus", "x:Year", "x:Hex", "x:Host" ),
                                          i -> ( i < 6 ) ? "Mastadenovirus" : "Aviadenovirus",
                                          i -> String.valueOf( 1990 + i ),
                                          i -> ( i < 6 ) ? String.valueOf( 1 + ( i % 3 ) ) : "0x1A",
                                          i -> ( i < 6 ) ? null : ( ( ( i % 2 ) == 0 ) ? "duck" : "chicken" ) );
        final List<PropertyColorScheme.VisCandidate> tree_cands = PropertyColorScheme.visualizationCandidates( t );
        final PropertyColorScheme.VisCandidate genus = PropertyColorScheme.findCandidate( tree_cands, "x:Genus" );
        final PropertyColorScheme.VisCandidate year = PropertyColorScheme.findCandidate( tree_cands, "x:Year" );
        final PropertyColorScheme.VisCandidate hex = PropertyColorScheme.findCandidate( tree_cands, "x:Hex" );
        final PropertyColorScheme.VisCandidate host = PropertyColorScheme.findCandidate( tree_cands, "x:Host" );
        if ( ( genus == null ) || ( year == null ) || ( hex == null ) || ( host == null ) || hex._numeric ) {
            return fail( "Genus, Year, Hex (categorical) and Host must all be tree candidates; got " + refsOf( tree_cands ) );
        }
        // the VIEW: collapse clade B, so only clade A is on screen (a collapse is a view, like a subtree)
        final PhylogenyNode clade_b = t.getRoot().getChildNode( 1 );
        clade_b.setCollapse( true );
        final List<PhylogenyNode> view = PropertyColorScheme.visibleExternalNodes( t );
        if ( view.size() != 6 ) {
            return fail( "the view should hold clade A's six tips, got " + view.size() );
        }
        // the one-value clade is SUMMARIZED, not refused: one legend row, and it is coloured rather than greyed
        final PropertyColorScheme.VisSummary gs = PropertyColorScheme.visualizationSummary( genus, view );
        if ( ( gs._distinct != 1 ) || ( gs._coverage != 6 ) || !"Mastadenovirus".equals( gs._values.get( 0 ) ) ) {
            return fail( "a one-genus clade must summarize to its one value" );
        }
        final PropertyColorScheme gsch = new PropertyColorScheme( t, "x:Genus", null,
                                                                  PropertyColorScheme.DEFAULT_PALETTE_NAME, null, null,
                                                                  Boolean.FALSE, genus );
        if ( gsch.isEmpty() || ( gsch.getValueColors().size() != 1 ) || ( gsch.missingCount() != 0 ) ) {
            return fail( "the chosen Genus must colour its one-value clade" );
        }
        // a NUMERIC candidate's band is the VIEW's: 12 years on the tree (range, switchable), 6 in the clade (colours)
        if ( ( PropertyColorScheme.colorModeBand( year, PropertyColorScheme.allExternalNodes( t ) )
                != PropertyColorScheme.ModeBand.MEDIUM )
                || ( PropertyColorScheme.colorModeBand( year, view ) != PropertyColorScheme.ModeBand.SMALL ) ) {
            return fail( "a numeric candidate's band must follow the distinct count in the view" );
        }
        // a CATEGORICAL candidate stays a category in every view, even where its visible values are all digits...
        if ( PropertyColorScheme.colorModeBand( hex, view ) != PropertyColorScheme.ModeBand.NOT_NUMERIC ) {
            return fail( "a categorical candidate must keep its mode in a view that happens to look numeric" );
        }
        // ...where classifying the view itself would have called it numeric (the control)
        if ( PropertyColorScheme.colorModeBand( view, "x:Hex" ) == PropertyColorScheme.ModeBand.NOT_NUMERIC ) {
            return fail( "control: per-view detection should read clade A's Hex values as numbers" );
        }
        // a chosen field ABSENT from the view stays chosen: an empty scheme whose every visible tip is "no value"
        final PropertyColorScheme hs = new PropertyColorScheme( t, "x:Host", null,
                                                                PropertyColorScheme.DEFAULT_PALETTE_NAME, null, null,
                                                                Boolean.FALSE, host );
        if ( !hs.isEmpty() || ( hs.missingCount() != 6 ) ) {
            return fail( "a chosen field absent from the view must be an empty scheme with 6 'no value' tips" );
        }
        // and the view changed nothing about what the TREE offers
        if ( !refsOf( PropertyColorScheme.visualizationCandidates( t ) ).equals( refsOf( tree_cands ) ) ) {
            return fail( "a collapse must not change the tree's candidates" );
        }
        clade_b.setCollapse( false );
        return true;
    }

    /**
     * An EDIT re-derives the candidates but KEEPS the fields the user chose for as long as they still carry a value on
     * a remaining tip. THE REFUSAL RULES DECIDE WHAT IS OFFERED, NEVER WHAT IS ALREADY CHOSEN (JS-authoritative,
     * decision 1). On the desktop an edit also includes a node-data write, which Archaeopteryx.js cannot do.
     */
    private static boolean testEditKeepsValuedChoice() {
        final Phylogeny t = twoCladeTree( 6, refs( "x:Host", "x:Zone" ),
                                          i -> ( i < 6 ) ? ( ( ( i % 2 ) == 0 ) ? "cat" : "dog" ) : "duck",
                                          i -> ( ( i % 2 ) == 0 ) ? "N" : "S" );
        final PropertyColorScheme.VisCandidate host = cand( t, "x:Host" );
        final PropertyColorScheme.VisCandidate zone = cand( t, "x:Zone" );
        if ( ( host == null ) || ( zone == null ) ) {
            return fail( "Host and Zone must be offered before the edit" );
        }
        // delete clade A: one host is left, so the rules now refuse Host...
        t.deleteSubtree( t.getRoot().getChildNode( 0 ), true );
        if ( cand( t, "x:Host" ) != null ) {
            return fail( "control: with one host left the field must not be OFFERED" );
        }
        // ...but the CHOSEN Host is kept, appended after the offered fields and flagged
        final List<PropertyColorScheme.VisCandidate> kept = PropertyColorScheme
                .visualizationCandidatesKeeping( t, java.util.Arrays.asList( host ) );
        final PropertyColorScheme.VisCandidate k = PropertyColorScheme.findCandidate( kept, "x:Host" );
        if ( ( k == null ) || !k._kept || ( kept.indexOf( k ) != ( kept.size() - 1 ) ) ) {
            return fail( "a chosen field with a value left must be KEPT, last; got " + refsOf( kept ) );
        }
        if ( !k._values.equals( java.util.Arrays.asList( "Duck" ) ) || ( k._coverage != 6 ) || ( k._total != 6 ) ) {
            return fail( "a kept field must be re-summarized over what is left; got " + k._values + " " + k._coverage
                    + "/" + k._total );
        }
        // a field that was NOT chosen is not kept
        if ( PropertyColorScheme.findCandidate( PropertyColorScheme.visualizationCandidatesKeeping( t, null ), "x:Host" ) != null ) {
            return fail( "an unchosen refused field must not be kept" );
        }
        // a chosen field that is still OFFERED comes back re-derived, not as a kept copy
        final PropertyColorScheme.VisCandidate z = PropertyColorScheme.findCandidate( PropertyColorScheme
                .visualizationCandidatesKeeping( t, java.util.Arrays.asList( zone, host ) ), "x:Zone" );
        if ( ( z == null ) || z._kept || ( z == zone ) ) {
            return fail( "a still-offered chosen field must be the re-derived candidate" );
        }
        // a NODE-DATA edit (desktop only) that removes the chosen field from every tip drops it: only a field with no
        // value left anywhere falls back to the default
        for( final PhylogenyNode n : PropertyColorScheme.allExternalNodes( t ) ) {
            final PropertiesList pl = new PropertiesList();
            for( final Property p : n.getNodeData().getProperties().getProperties() ) {
                if ( !"x:Host".equals( p.getRef() ) ) {
                    pl.addProperty( p );
                }
            }
            n.getNodeData().setProperties( pl );
        }
        if ( PropertyColorScheme.findCandidate( PropertyColorScheme.visualizationCandidatesKeeping( t, java.util.Arrays.asList( k ) ),
                                                "x:Host" ) != null ) {
            return fail( "a chosen field with no value left anywhere must be dropped" );
        }
        return true;
    }

    /**
     * THE ACCEPTANCE TEST: both Archaeopteryx.js contract fixtures, run through the desktop's own phyloXML parser and
     * classifier -- the NAMES (vis-contract.tsv: verdict, label, deprioritization) and the DATA (vis-trees.tsv +
     * vis-trees/*.xml: every column). Copied from the JS repo's test/fixtures (at 635741b) into
     * forester/test_data/vis_contract/. When they regenerate, recopy; never edit a row to make this pass. The row
     * counts are pinned so a truncated copy cannot pass.
     */
    private static boolean testContractFixtures() {
        final java.io.File dir = new java.io.File( System.getProperty( "user.dir" ), "forester/test_data/vis_contract" );
        try {
            int name_rows = 0;
            for( final String line : java.nio.file.Files.readAllLines( new java.io.File( dir, "vis-contract.tsv" ).toPath() ) ) {
                if ( line.isEmpty() || line.startsWith( "#" ) ) {
                    continue;
                }
                final String[] f = line.split( "\t", -1 );
                ++name_rows;
                final String ref = f[ 0 ];
                final String verdict = f[ 1 ];
                final String label = ( f.length > 2 ) ? f[ 2 ] : "";
                final boolean excluded = PropertyColorScheme.isExcludedRef( ref );
                if ( excluded != "EXCLUDED".equals( verdict ) ) {
                    return fail( "name fixture: " + ref + " js=" + verdict + " desktop excluded=" + excluded );
                }
                if ( !label.equals( PropertyColorScheme.displayName( ref ) ) ) {
                    return fail( "name fixture label: " + ref + " js='" + label + "' desktop='"
                            + PropertyColorScheme.displayName( ref ) + "'" );
                }
                if ( !excluded && ( PropertyColorScheme.isDeprioritizedRef( ref ) != "DEPRIORITIZED".equals( verdict ) ) ) {
                    return fail( "name fixture deprioritization: " + ref + " js=" + verdict );
                }
            }
            if ( name_rows != 182 ) {
                return fail( "name fixture: expected 182 rows, read " + name_rows );
            }
            final Map<String, List<String[]>> by_tree = new java.util.LinkedHashMap<String, List<String[]>>();
            for( final String line : java.nio.file.Files.readAllLines( new java.io.File( dir, "vis-trees.tsv" ).toPath() ) ) {
                if ( line.isEmpty() || line.startsWith( "#" ) ) {
                    continue;
                }
                final String[] f = line.split( "\t", -1 );
                if ( !by_tree.containsKey( f[ 0 ] ) ) {
                    by_tree.put( f[ 0 ], new ArrayList<String[]>() );
                }
                by_tree.get( f[ 0 ] ).add( f );
            }
            int data_rows = 0;
            for( final Map.Entry<String, List<String[]>> e : by_tree.entrySet() ) {
                final java.io.File xml = new java.io.File( new java.io.File( dir, "vis-trees" ), e.getKey() + ".xml" );
                // STRICT, XSD-validating parse for every tree: a fixture that is not valid phyloXML fails here (JS c757959
                // fixed the one that was -- a ref containing a space)
                final Phylogeny phy = org.forester.phylogeny.factories.ParserBasedPhylogenyFactory.getInstance()
                        .create( xml, org.forester.io.parsers.util.ParserUtils.createParserDependingOnFileType( xml, true ) )[ 0 ];
                final List<PropertyColorScheme.VisCandidate> cands = PropertyColorScheme.visualizationCandidates( phy );
                final PropertyColorScheme.VisCandidate first = PropertyColorScheme.openingVisualization( cands );
                final String opens = ( first == null ) ? "-" : first._ref;
                for( final String[] want : e.getValue() ) {
                    ++data_rows;
                    final PropertyColorScheme.VisCandidate c = PropertyColorScheme.findCandidate( cands, want[ 1 ] );
                    final String got = ( c == null )
                            ? String.join( "\t", e.getKey(), want[ 1 ], "refused", "-", "-", "-", "-", "-", "-", "-", "-",
                                           "-", opens, "-" )
                            : String.join( "\t", e.getKey(), want[ 1 ], "offered", String.valueOf( c.tier() ),
                                           String.valueOf( cands.indexOf( c ) ), bit( c._numeric ), bit( c._wide ),
                                           bit( c._near_unique ), bit( c._sparse ), bit( c._deprioritized ),
                                           c._color_mode, bit( c._shape ), opens, String.valueOf( c._values.size() ) );
                    if ( !got.equals( String.join( "\t", want ) ) ) {
                        return fail( "data fixture row differs:\n    js      " + String.join( "\t", want )
                                + "\n    desktop " + got );
                    }
                }
                for( final PropertyColorScheme.VisCandidate c : cands ) {
                    boolean declared = false;
                    for( final String[] want : e.getValue() ) {
                        declared |= want[ 1 ].equals( c._ref );
                    }
                    if ( !declared ) {
                        return fail( "data fixture: " + e.getKey() + " offers undeclared " + c._ref );
                    }
                }
            }
            if ( data_rows != 39 ) {
                return fail( "data fixture: expected 39 rows, read " + data_rows );
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    /**
     * A property value that folds to NOTHING -- "_", "___", a host that is only its ";" qualifier -- is NO VALUE,
     * everywhere, exactly like an empty string (Christian, 2026-09-12; forester.js 635741b): the tip is not covered, no
     * group or legend row is made, it does not make the ref "carried twice", and a field of nothing else is no
     * candidate.
     */
    private static boolean testFoldToEmptyIsNoValue() {
        final Phylogeny t = fieldTree( 30, refs( "x:Blanks", "x:Host", "x:Void", "x:Beside" ),
                                       i -> ( i < 12 ) ? ( ( ( i % 2 ) == 0 ) ? "_" : "___" )
                                               : ( ( ( i % 2 ) == 0 ) ? "A" : "B" ),
                                       i -> ( i < 12 ) ? "; cell culture" : ( ( ( i % 2 ) == 0 ) ? "Human; male" : "Pig" ),
                                       i -> "_",
                                       i -> ( i == 0 ) ? new String[] { "_", "Human" }
                                               : ( ( ( i % 2 ) == 0 ) ? "Human" : "Pig" ) );
        final PropertyColorScheme.VisCandidate blanks = cand( t, "x:Blanks" );
        if ( ( blanks == null ) || ( blanks._coverage != 18 ) || !blanks._sparse || ( blanks._values.size() != 2 ) ) {
            return fail( "12 of 30 tips of \"_\"/\"___\" must leave 18 covered, sparse, 2 values; got "
                    + ( ( blanks == null ) ? null : ( blanks._coverage + " " + blanks._values ) ) );
        }
        final PropertyColorScheme.VisCandidate host = cand( t, "x:Host" );
        if ( ( host == null ) || ( host._coverage != 18 ) || !host._sparse
                || !host._values.equals( java.util.Arrays.asList( "Human", "Pig" ) ) ) {
            return fail( "a host that is only its ';' qualifier must be no value; got "
                    + ( ( host == null ) ? null : ( host._coverage + " " + host._values ) ) );
        }
        if ( cand( t, "x:Void" ) != null ) {
            return fail( "a field whose every value folds to nothing must not be a candidate" );
        }
        final PropertyColorScheme.VisCandidate beside = cand( t, "x:Beside" );
        if ( ( beside == null ) || ( beside._coverage != 30 ) ) {
            return fail( "\"_\" beside \"Human\" on one tip must not make the ref carried twice" );
        }
        PhylogenyNode t0 = null;
        for( final PhylogenyNode n : t.getExternalNodes() ) {
            if ( "t0".equals( n.getName() ) ) {
                t0 = n;
            }
        }
        if ( !"Human".equals( PropertyColorScheme.visualizationNodeValue( t0, beside ) ) ) {
            return fail( "the node value must skip \"_\" and read the tip's real value; got "
                    + PropertyColorScheme.visualizationNodeValue( t0, beside ) );
        }
        if ( PropertyColorScheme.visualizationNodeValue( t0, blanks ) != null ) {
            return fail( "a tip carrying only \"_\" must read as no value" );
        }
        final PropertyColorScheme.VisSummary sum = PropertyColorScheme
                .visualizationSummary( blanks, PropertyColorScheme.allExternalNodes( t ) );
        if ( ( sum._coverage != 18 ) || sum._counts.containsKey( "" ) ) {
            return fail( "the summary must not cover a fold-to-empty tip or give it a row; got " + sum._counts );
        }
        final PropertyColorScheme s = new PropertyColorScheme( t, "x:Blanks", null, PropertyColorScheme.DEFAULT_PALETTE_NAME,
                                                               null, null, Boolean.FALSE, blanks );
        if ( ( s.getValueColors().size() != 2 ) || ( s.missingCount() != 12 ) ) {
            return fail( "the legend must have 2 rows and 12 'no value' tips; got " + s.getValueCounts() + ", missing "
                    + s.missingCount() );
        }
        return true;
    }

    private static String bit( final boolean b ) {
        return b ? "1" : "0";
    }

    /** A tip's value for a field: null = absent, a String[] = the ref carried once per element. */
    private interface V {

        Object at( int i );
    }

    private static String[] refs( final String... r ) {
        return r;
    }

    private static PropertyColorScheme.VisCandidate cand( final Phylogeny p, final String ref ) {
        return PropertyColorScheme.findCandidate( PropertyColorScheme.visualizationCandidates( p ), ref );
    }

    private static List<String> refsOf( final List<PropertyColorScheme.VisCandidate> cs ) {
        final List<String> r = new ArrayList<String>();
        for( final PropertyColorScheme.VisCandidate c : cs ) {
            r.add( c._ref );
        }
        return r;
    }

    private static PhylogenyNode fieldTip( final int i, final String[] refs, final V[] fns ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "t" + i );
        final PropertiesList pl = new PropertiesList();
        for( int k = 0; k < refs.length; ++k ) {
            final Object v = fns[ k ].at( i );
            if ( v instanceof String[] ) {
                for( final String one : (String[]) v ) {
                    pl.addProperty( new Property( refs[ k ], one, "", "xsd:string", AppliesTo.NODE ) );
                }
            }
            else if ( v != null ) {
                pl.addProperty( new Property( refs[ k ], (String) v, "", "xsd:string", AppliesTo.NODE ) );
            }
        }
        n.getNodeData().setProperties( pl );
        return n;
    }

    /** A star tree of {@code n} tips; tip i carries each ref's value from its function. */
    private static Phylogeny fieldTree( final int n, final String[] refs, final V... fns ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < n; ++i ) {
            root.addAsChild( fieldTip( i, refs, fns ) );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /** Two clades of {@code per} tips under the root; tips are numbered 0.. across both clades. */
    private static Phylogeny twoCladeTree( final int per, final String[] refs, final V... fns ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int c = 0; c < 2; ++c ) {
            final PhylogenyNode clade = new PhylogenyNode();
            for( int i = 0; i < per; ++i ) {
                clade.addAsChild( fieldTip( ( c * per ) + i, refs, fns ) );
            }
            root.addAsChild( clade );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /** A star tree of {@code tips} leaves, the first {@code annotated} of them carrying {@code x:field}. */
    private static Phylogeny propTree( final int tips, final int annotated, final boolean unique ) throws Exception {
        final StringBuilder sb = new StringBuilder( "(" );
        for( int i = 0; i < tips; ++i ) {
            if ( i > 0 ) {
                sb.append( ',' );
            }
            sb.append( 't' ).append( i );
        }
        sb.append( ");" );
        final Phylogeny p = Phylogeny.createInstanceFromNhxString( sb.toString() );
        int i = 0;
        for( final PhylogenyNode n : p.getExternalNodes() ) {
            if ( i < annotated ) {
                final PropertiesList pl = new PropertiesList();
                pl.addProperty( new Property( "x:field", unique ? ( "v" + i ) : ( "v" + ( i % 2 ) ), "", "xsd:string",
                                              AppliesTo.NODE ) );
                n.getNodeData().setProperties( pl );
            }
            ++i;
        }
        return p;
    }

}
