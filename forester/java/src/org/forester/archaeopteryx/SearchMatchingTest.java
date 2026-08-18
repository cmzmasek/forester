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

import java.util.List;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.NDF;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headless coverage for the redesigned search backend ({@link SearchField} / {@link SearchMode} /
 * {@link SearchSpec} / {@link SearchMatcher}). The legacy search engine had essentially no matching tests
 * (one field-restriction case in {@code org.forester.test.Test}); this pins the new semantics: every string
 * mode, case sensitivity, the "Any text" field, per-property (key-scoped) matching, all numeric operators
 * incl. range boundaries, the inverse flag, and a whole-tree search.
 */
public final class SearchMatchingTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SearchMatching: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            stringModes();
            caseSensitivity();
            anyTextField();
            perPropertyScoping();
            numericOperators();
            rangeAndValuelessNumeric();
            numericDatatypeClassifier();
            inverse();
            specValidation();
            parseDouble();
            wholeWordBoundaries();
            regexValidity();
            wholeTreeSearch();
            availableFieldsScan();
            structuralFields();
            return true;
        }
        catch ( final AssertionError e ) {
            System.out.println( "  [SearchMatchingTest] " + e.getMessage() );
            return false;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    // ---- string modes ---------------------------------------------------------------------------------------

    private static void stringModes() {
        final PhylogenyNode n = new PhylogenyNode();
        n.getNodeData().setTaxonomy( sci( "Homo sapiens" ) );
        final SearchField f = SearchField.ofNdf( NDF.TaxonomyScientificName );
        ck( pos( f, SearchMode.CONTAINS, "sapi", n ), "CONTAINS should match a substring" );
        ck( !pos( f, SearchMode.CONTAINS, "xyz", n ), "CONTAINS should not match an absent substring" );
        ck( pos( f, SearchMode.STARTS_WITH, "Homo", n ), "STARTS_WITH should match the prefix" );
        ck( !pos( f, SearchMode.STARTS_WITH, "sapiens", n ), "STARTS_WITH should not match a non-prefix" );
        ck( pos( f, SearchMode.ENDS_WITH, "sapiens", n ), "ENDS_WITH should match the suffix" );
        ck( !pos( f, SearchMode.ENDS_WITH, "Homo", n ), "ENDS_WITH should not match a non-suffix" );
        ck( pos( f, SearchMode.WHOLE_WORD, "Homo", n ), "WHOLE_WORD should match a whole term" );
        ck( !pos( f, SearchMode.WHOLE_WORD, "Hom", n ), "WHOLE_WORD should not match a partial term" );
        ck( pos( f, SearchMode.REGEX, "Ho.*ens", n ), "REGEX should match a valid pattern" );
        // an invalid regex must fail closed, never throw (mutant: let PatternSyntaxException escape)
        ck( !pos( f, SearchMode.REGEX, "[", n ), "an invalid regex should return false, not throw" );
    }

    private static void caseSensitivity() {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "Kinase" );
        final SearchField f = SearchField.ofNdf( NDF.NodeName );
        ck( SearchMatcher.matchesPositive( new SearchSpec( f, SearchMode.CONTAINS, "kinase", null, false, false ), n ),
            "case-insensitive CONTAINS should match regardless of case" );
        ck( !SearchMatcher.matchesPositive( new SearchSpec( f, SearchMode.CONTAINS, "kinase", null, true, false ), n ),
            "case-sensitive CONTAINS should respect case" );
        // regex honours case sensitivity too
        ck( SearchMatcher.matchesPositive( new SearchSpec( f, SearchMode.REGEX, "kin.se", null, false, false ), n ),
            "case-insensitive REGEX should match" );
        ck( !SearchMatcher.matchesPositive( new SearchSpec( f, SearchMode.REGEX, "kin.se", null, true, false ), n ),
            "case-sensitive REGEX should respect case" );
    }

    // ---- Any Text -------------------------------------------------------------------------------------

    private static void anyTextField() {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "leaf1" );
        n.getNodeData().setTaxonomy( sci( "Felis catus" ) );
        addProp( n, "data:host", "human", "xsd:string" );             // a user-visible property
        addProp( n, "aptx:import_profile", "v1;secret", "xsd:string" ); // an internal (aptx:) property
        final SearchField any = SearchField.anyText();
        ck( pos( any, SearchMode.CONTAINS, "catus", n ), "Any-text should search the scientific name" );
        ck( pos( any, SearchMode.CONTAINS, "leaf1", n ), "Any-text should search the node name" );
        ck( pos( any, SearchMode.CONTAINS, "human", n ), "Any-text should search user-visible property values" );
        ck( !pos( any, SearchMode.CONTAINS, "secret", n ),
            "Any-text must NOT search internal (aptx:) property values" );
    }

    // ---- per-property (key-scoped) matching -----------------------------------------------------------------

    private static void perPropertyScoping() {
        final PhylogenyNode n = new PhylogenyNode();
        addProp( n, "data:host", "human", "xsd:string" );
        addProp( n, "data:diet", "omnivore", "xsd:string" );
        ck( pos( SearchField.property( "data:host", false ), SearchMode.CONTAINS, "human", n ),
            "a property field should match its own value" );
        // the legacy search lumped all property VALUES together and ignored the key; the redesign is key-scoped
        ck( !pos( SearchField.property( "data:diet", false ), SearchMode.CONTAINS, "human", n ),
            "a property field must not match a value living in a different property (key-scoped)" );
    }

    // ---- numeric operators ----------------------------------------------------------------------------------

    private static void numericOperators() {
        final PhylogenyNode n = new PhylogenyNode();
        n.setDistanceToParent( 0.5 );
        PhylogenyMethods.setConfidence( n, 95 );
        addProp( n, "aptx:reads", "1000", "xsd:decimal" );
        final SearchField bl = SearchField.branchLength();
        ck( pos( bl, SearchMode.EQ, "0.5", n ), "EQ should match the exact branch length" );
        ck( !pos( bl, SearchMode.EQ, "0.6", n ), "EQ should not match a different value" );
        ck( pos( bl, SearchMode.NE, "0.6", n ), "NE should match a differing value" );
        ck( !pos( bl, SearchMode.NE, "0.5", n ), "NE should not match an equal value" );
        ck( pos( bl, SearchMode.LT, "1", n ), "LT should match a smaller-than value" );
        ck( !pos( bl, SearchMode.LT, "0.5", n ), "LT is strict" );
        ck( pos( bl, SearchMode.LE, "0.5", n ), "LE is inclusive" );
        ck( pos( bl, SearchMode.GT, "0.4", n ), "GT should match a greater-than value" );
        ck( !pos( bl, SearchMode.GT, "0.5", n ), "GT is strict" );
        ck( pos( bl, SearchMode.GE, "0.5", n ), "GE is inclusive" );
        // confidence + numeric property share the numeric path
        ck( pos( SearchField.confidence(), SearchMode.GE, "90", n ), "confidence GE should match" );
        ck( !pos( SearchField.confidence(), SearchMode.LT, "90", n ), "confidence LT should not match here" );
        ck( pos( SearchField.property( "aptx:reads", true ), SearchMode.GT, "500", n ),
            "a numeric property GT should match" );
        ck( !pos( SearchField.property( "aptx:reads", true ), SearchMode.LT, "500", n ),
            "a numeric property LT should not match here" );
    }

    private static void rangeAndValuelessNumeric() {
        final PhylogenyNode n = new PhylogenyNode();
        n.setDistanceToParent( 0.5 );
        final SearchField bl = SearchField.branchLength();
        ck( range( bl, "0.4", "0.6", n ), "RANGE should contain an interior value" );
        ck( !range( bl, "0.6", "0.9", n ), "RANGE should exclude an out-of-range value" );
        ck( range( bl, "0.5", "0.6", n ), "RANGE is inclusive at the low bound" );
        ck( range( bl, "0.4", "0.5", n ), "RANGE is inclusive at the high bound" );
        // swapped bounds are normalised (min/max), so order does not matter
        ck( range( bl, "0.6", "0.4", n ), "RANGE should normalise swapped bounds" );
        // a value-less node (default branch length) matches no numeric predicate, including NE
        final PhylogenyNode empty = new PhylogenyNode();
        ck( !pos( bl, SearchMode.EQ, "0", empty ), "a value-less numeric field should not match EQ" );
        ck( !pos( bl, SearchMode.NE, "0", empty ), "a value-less numeric field should not match NE" );
        ck( !range( bl, "0", "10", empty ), "a value-less numeric field should not match RANGE" );
    }

    private static void numericDatatypeClassifier() {
        ck( SearchField.datatypeIsNumeric( "xsd:decimal" ), "xsd:decimal is numeric" );
        ck( SearchField.datatypeIsNumeric( "xsd:integer" ), "xsd:integer is numeric" );
        ck( SearchField.datatypeIsNumeric( "xsd:double" ), "xsd:double is numeric" );
        ck( SearchField.datatypeIsNumeric( "xsd:long" ), "xsd:long is numeric" );
        ck( !SearchField.datatypeIsNumeric( "xsd:string" ), "xsd:string is not numeric" );
        ck( !SearchField.datatypeIsNumeric( null ), "a null datatype is not numeric" );
        ck( !SearchField.datatypeIsNumeric( "" ), "an empty datatype is not numeric" );
        // a custom datatype whose local name merely CONTAINS a numeric-type substring must not be misclassified
        ck( !SearchField.datatypeIsNumeric( "x:footprint" ), "x:footprint (contains 'int') is not numeric" );
        ck( !SearchField.datatypeIsNumeric( "data:print" ), "data:print (contains 'int') is not numeric" );
    }

    // ---- inverse --------------------------------------------------------------------------------------------

    private static void inverse() {
        final PhylogenyNode n = new PhylogenyNode();
        n.getNodeData().setTaxonomy( sci( "Homo sapiens" ) );
        final SearchField f = SearchField.ofNdf( NDF.TaxonomyScientificName );
        ck( !SearchMatcher.matches( new SearchSpec( f, SearchMode.CONTAINS, "sapi", null, false, true ), n ),
            "inverse should flip a positive match to false" );
        ck( SearchMatcher.matches( new SearchSpec( f, SearchMode.CONTAINS, "xyz", null, false, true ), n ),
            "inverse should flip a negative match to true" );
    }

    // ---- SearchSpec validation ------------------------------------------------------------------------------

    private static void specValidation() {
        // a numeric mode on a string field (and vice versa) is a programming error
        ck( throwsIae( () -> new SearchSpec( SearchField.ofNdf( NDF.NodeName ), SearchMode.GT, "5" ) ),
            "a numeric mode on a string field should be rejected" );
        ck( throwsIae( () -> new SearchSpec( SearchField.branchLength(), SearchMode.CONTAINS, "x" ) ),
            "a string mode on a numeric field should be rejected" );
    }

    private static void parseDouble() {
        ck( ( SearchMatcher.parseFiniteDouble( "1.5" ) != null )
                && ( SearchMatcher.parseFiniteDouble( "1.5" ).doubleValue() == 1.5 ), "\"1.5\" should parse to 1.5" );
        ck( SearchMatcher.parseFiniteDouble( "abc" ) == null, "non-numeric should parse to null" );
        ck( SearchMatcher.parseFiniteDouble( "" ) == null, "empty should parse to null" );
        ck( SearchMatcher.parseFiniteDouble( "NaN" ) == null, "NaN should be rejected" );
        ck( SearchMatcher.parseFiniteDouble( "Infinity" ) == null, "Infinity should be rejected" );
        // a comma is accepted as the decimal separator when unambiguous (European locales)
        ck( eqD( SearchMatcher.parseFiniteDouble( "0,5" ), 0.5 ), "\"0,5\" should parse to 0.5" );
        ck( eqD( SearchMatcher.parseFiniteDouble( "1,5" ), 1.5 ), "\"1,5\" should parse to 1.5" );
        ck( eqD( SearchMatcher.parseFiniteDouble( "-2,25" ), -2.25 ), "\"-2,25\" should parse to -2.25" );
        ck( SearchMatcher.parseFiniteDouble( "1,5,5" ) == null, "two commas should not parse" );
        ck( SearchMatcher.parseFiniteDouble( "1.5,5" ) == null, "a mixed '.' and ',' should not parse" );
        // the US thousands-grouping pattern (a comma followed by EXACTLY three digits) is ambiguous with a European
        // decimal, so it is REJECTED (null, as before the comma support) rather than misparsed as a fraction
        ck( SearchMatcher.parseFiniteDouble( "12,500" ) == null, "\"12,500\" (US thousands) must be null, not 12.5" );
        ck( SearchMatcher.parseFiniteDouble( "1,000" ) == null, "\"1,000\" (ambiguous thousands) must be null" );
        ck( eqD( SearchMatcher.parseFiniteDouble( "1,25" ), 1.25 ), "\"1,25\" (2 decimals) should still parse to 1.25" );
        ck( eqD( SearchMatcher.parseFiniteDouble( "1,2345" ), 1.2345 ),
            "\"1,2345\" (4 decimals, not a thousands group) should parse to 1.2345" );
    }

    private static boolean eqD( final Double d, final double v ) {
        return ( d != null ) && ( d.doubleValue() == v );
    }

    // ---- whole-word / regex ---------------------------------------------------------------------------------

    private static void wholeWordBoundaries() {
        final PhylogenyNode n = new PhylogenyNode();
        final SearchField f = SearchField.ofNdf( NDF.NodeName );
        n.setName( "Homo_sapiens" );
        ck( pos( f, SearchMode.WHOLE_WORD, "Homo", n ), "whole-word should match across an underscore (Homo_sapiens)" );
        ck( pos( f, SearchMode.WHOLE_WORD, "sapiens", n ), "whole-word should match the trailing token" );
        n.setName( "Rot1-kinase" );
        ck( pos( f, SearchMode.WHOLE_WORD, "kinase", n ), "whole-word should match across a hyphen (Rot1-kinase)" );
        n.setName( "kinase," );
        ck( pos( f, SearchMode.WHOLE_WORD, "kinase", n ), "whole-word should match a token before a comma" );
        n.setName( "kinases" );
        ck( !pos( f, SearchMode.WHOLE_WORD, "kinase", n ), "whole-word should NOT match a longer word (kinases)" );
        n.setName( "prokinase" );
        ck( !pos( f, SearchMode.WHOLE_WORD, "kinase", n ), "whole-word should NOT match a suffix (prokinase)" );
        // accented / non-ASCII letters are part of the word (Unicode-aware boundary), NOT a token boundary
        n.setName( "Bäcker" );
        ck( !pos( f, SearchMode.WHOLE_WORD, "B", n ), "whole-word 'B' should NOT match before a diacritic (Bäcker)" );
        n.setName( "naïve" );
        ck( !pos( f, SearchMode.WHOLE_WORD, "na", n ), "whole-word 'na' should NOT split at a diacritic (naïve)" );
        ck( pos( f, SearchMode.WHOLE_WORD, "naïve", n ), "whole-word should match the full accented word (naïve)" );
    }

    private static void regexValidity() {
        ck( SearchMatcher.isCompilableRegex( "Ho.*ens" ), "a valid regex is compilable" );
        ck( !SearchMatcher.isCompilableRegex( "[" ), "an invalid regex is not compilable" );
        ck( !SearchMatcher.isCompilableRegex( "" ), "an empty query is not 'compilable'" );
    }

    // ---- whole-tree search ----------------------------------------------------------------------------------

    private static void wholeTreeSearch() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode leaf1 = new PhylogenyNode();
        leaf1.getNodeData().setTaxonomy( sci( "Homo sapiens" ) );
        final PhylogenyNode leaf2 = new PhylogenyNode();
        leaf2.getNodeData().setTaxonomy( sci( "Felis catus" ) );
        root.addAsChild( leaf1 );
        root.addAsChild( leaf2 );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();

        final SearchField any = SearchField.anyText();
        final List<Long> hits = SearchMatcher.search( new SearchSpec( any, SearchMode.CONTAINS, "catus" ), phy );
        ck( ( hits.size() == 1 ) && hits.contains( leaf2.getId() ),
            "whole-tree search should return exactly the matching leaf" );
        // inverse over the tree: root (no data) and leaf1 do NOT contain "catus", so they are the complement
        final List<Long> inv = SearchMatcher
                .search( new SearchSpec( any, SearchMode.CONTAINS, "catus", null, false, true ), phy );
        ck( ( inv.size() == 2 ) && inv.contains( root.getId() ) && inv.contains( leaf1.getId() )
                && !inv.contains( leaf2.getId() ), "inverse whole-tree search should return the complement" );
    }

    // ---- available (per-tree) field discovery ---------------------------------------------------------------

    private static void availableFieldsScan() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode leaf = new PhylogenyNode();
        leaf.setName( "leaf1" );
        leaf.getNodeData().setTaxonomy( sci( "Homo sapiens" ) );
        leaf.setDistanceToParent( 0.5 );
        PhylogenyMethods.setConfidence( leaf, 95 );
        addProp( leaf, "data:host", "human", "xsd:string" );
        addProp( leaf, "data:reads", "1000", "xsd:decimal" );
        addProp( leaf, "x:num", "5", "" );   // no datatype, value parses -> numeric
        addProp( leaf, "x:txt", "abc", "" ); // no datatype, value does not parse -> string
        addProp( leaf, "aptx:import_profile", "v1;x", "xsd:string" ); // internal -> must NOT be offered as a field
        root.addAsChild( leaf );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();

        final List<SearchField> fields = SearchField.availableFields( phy );
        ck( byLabel( fields, "Any Text" ) != null, "availableFields should always include Any text" );
        ck( byLabel( fields, "Node Name" ) != null, "availableFields should always include Node Name" );
        ck( byLabel( fields, "Taxonomy Scientific" ) != null, "a present scientific name should be offered" );
        ck( byLabel( fields, "Branch Length" ) != null, "a present branch length should be offered" );
        ck( byLabel( fields, "Support / Confidence" ) != null, "a present confidence should be offered" );
        ck( byLabel( fields, "Gene Name" ) == null, "an absent gene name should not be offered" );
        ck( byLabel( fields, "Domain" ) == null, "an absent domain should not be offered" );
        final SearchField host = byLabel( fields, "data:host" );
        final SearchField reads = byLabel( fields, "data:reads" );
        final SearchField x_num = byLabel( fields, "x:num" );
        final SearchField x_txt = byLabel( fields, "x:txt" );
        ck( byLabel( fields, "aptx:import_profile" ) == null,
            "an internal (aptx:) property must NOT be offered as a searchable field" );
        ck( ( host != null ) && !host.isNumeric(), "data:host (xsd:string) should be a string field" );
        ck( ( reads != null ) && reads.isNumeric(), "data:reads (xsd:decimal) should be a numeric field" );
        ck( ( x_num != null ) && x_num.isNumeric(), "x:num (no datatype, numeric values) should be numeric" );
        ck( ( x_txt != null ) && !x_txt.isNumeric(), "x:txt (no datatype, non-numeric value) should be string" );
        ck( SearchField.availableFields( new Phylogeny() ).size() == 2,
            "an empty tree should offer only Any text + Node Name" );
    }

    // ---- topological / structural fields --------------------------------------------------------------------

    private static void structuralFields() {
        // root -> ( mid -> (a, b), c )
        final PhylogenyNode a = named( "a" );
        final PhylogenyNode b = named( "b" );
        final PhylogenyNode c = named( "c" );
        final PhylogenyNode mid = new PhylogenyNode();
        mid.addAsChild( a );
        mid.addAsChild( b );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( mid );
        root.addAsChild( c );
        mid.setDistanceToParent( 1.0 ); // branch lengths, so distance-from-root is meaningful/offered
        a.setDistanceToParent( 2.0 );
        b.setDistanceToParent( 3.0 );
        c.setDistanceToParent( 0.5 );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        phy.recalculateNumberOfExternalDescendants( false ); // populate the subtree tip counts

        final SearchField size = SearchField.cladeSize();
        final SearchField depth = SearchField.depth();
        ck( size.isNumeric() && depth.isNumeric(), "clade size and depth are numeric fields" );

        // clade size (subtree tips): root=3, mid=2, a/b/c=1
        ck( pos( size, SearchMode.EQ, "3", root ), "root clade size should be 3" );
        ck( pos( size, SearchMode.EQ, "2", mid ), "mid clade size should be 2" );
        ck( pos( size, SearchMode.GT, "1", mid ) && pos( size, SearchMode.GT, "1", root ),
            "internal nodes have clade size > 1" );
        ck( !pos( size, SearchMode.GT, "1", a ), "a leaf's clade size (1) is not > 1" );

        // number of DIRECT children: root=2, mid=2, a/b/c=0
        final SearchField kids = SearchField.numChildren();
        ck( pos( kids, SearchMode.EQ, "2", root ), "root has 2 children" );
        ck( pos( kids, SearchMode.EQ, "0", a ), "a leaf has 0 children" );
        ck( !pos( kids, SearchMode.GT, "0", c ), "a leaf's child count (0) is not > 0" );

        // depth from root (edges): root=0, mid=1, c=1, a/b=2
        ck( pos( depth, SearchMode.EQ, "0", root ), "root depth should be 0" );
        ck( pos( depth, SearchMode.EQ, "1", mid ), "mid depth should be 1" );
        ck( pos( depth, SearchMode.EQ, "2", a ), "a depth should be 2" );
        ck( pos( depth, SearchMode.LT, "2", root ), "root depth is < 2" );

        // distance from root (branch lengths): root=0, mid=1, c=0.5, a=3, b=4
        final SearchField dist = SearchField.distanceToRoot();
        ck( pos( dist, SearchMode.EQ, "0", root ), "root distance should be 0" );
        ck( pos( dist, SearchMode.EQ, "3", a ), "a distance should be 3 (1+2)" );
        ck( pos( dist, SearchMode.GT, "3.5", b ), "b distance (4) should be > 3.5" );
        ck( !pos( dist, SearchMode.GT, "3.5", a ), "a distance (3) is not > 3.5" );

        // they are always offered (distance only when the tree has branch lengths)
        final List<SearchField> fields = SearchField.availableFields( phy );
        ck( byLabel( fields, "Structure: Clade Size (tips)" ) != null, "clade-size field should be offered" );
        ck( byLabel( fields, "Structure: Number of Children" ) != null, "num-children field should be offered" );
        ck( byLabel( fields, "Structure: Depth from Root (edges)" ) != null, "depth field should be offered" );
        ck( byLabel( fields, "Structure: Distance from Root" ) != null,
            "distance-from-root should be offered when branch lengths are present" );

        // ... but distance-from-root is NOT offered for a tree without branch lengths
        final PhylogenyNode r2 = new PhylogenyNode();
        r2.addAsChild( named( "x" ) );
        final Phylogeny no_bl = new Phylogeny();
        no_bl.setRoot( r2 );
        no_bl.externalNodesHaveChanged();
        ck( byLabel( SearchField.availableFields( no_bl ), "Structure: Distance from Root" ) == null,
            "distance-from-root should NOT be offered without branch lengths" );
    }

    private static PhylogenyNode named( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static SearchField byLabel( final List<SearchField> fields, final String label ) {
        for ( final SearchField f : fields ) {
            if ( f.label().equals( label ) ) {
                return f;
            }
        }
        return null;
    }

    // ---- helpers --------------------------------------------------------------------------------------------

    private static boolean pos( final SearchField f, final SearchMode m, final String v, final PhylogenyNode n ) {
        return SearchMatcher.matchesPositive( new SearchSpec( f, m, v ), n );
    }

    private static boolean range( final SearchField f, final String lo, final String hi, final PhylogenyNode n ) {
        return SearchMatcher.matchesPositive( new SearchSpec( f, SearchMode.RANGE, lo, hi, false, false ), n );
    }

    private static Taxonomy sci( final String s ) {
        final Taxonomy t = new Taxonomy();
        t.setScientificName( s );
        return t;
    }

    private static void addProp( final PhylogenyNode n, final String ref, final String value, final String datatype ) {
        PropertiesList pl = n.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            n.getNodeData().setProperties( pl );
        }
        pl.addProperty( new Property( ref, value, "", datatype, AppliesTo.NODE ) );
    }

    private static boolean throwsIae( final Runnable r ) {
        try {
            r.run();
            return false;
        }
        catch ( final IllegalArgumentException e ) {
            return true;
        }
    }

    private static void ck( final boolean cond, final String msg ) {
        if ( !cond ) {
            throw new AssertionError( msg );
        }
    }

    private SearchMatchingTest() {
    }
}
