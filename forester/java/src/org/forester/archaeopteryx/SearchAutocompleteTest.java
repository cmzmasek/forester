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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

import javax.swing.JTextField;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods.NDF;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;

/**
 * Headless coverage for the search value-autocomplete popup ({@link SearchValueAutocomplete}) and the categorical
 * {@link SearchField#nodeType()} field it reintroduces. Exercises the pure filter / display-model logic (substring,
 * prefix-first ordering, case honouring, cap + hint row), the distinct-value extraction the popup feeds on (empty for
 * "Any text" / numeric fields), and an accept (fills the field + runs the search) -- all without realising a window,
 * so it runs for real even in the headless suite.
 */
public final class SearchAutocompleteTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SearchAutocomplete: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            filterOrdering();
            capAndHint();
            distinctValuesAndNodeType();
            acceptFillsAndSearches();
            return true;
        }
        catch ( final AssertionError e ) {
            System.out.println( "  [SearchAutocompleteTest] " + e.getMessage() );
            return false;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    // ---- pure filter / display-model logic ------------------------------------------------------------------

    private static void filterOrdering() {
        // sorted like SearchField.distinctValues (a TreeSet): [Homo sapiens, Pan, Sapporo virus, sapiens]
        final List<String> all = sorted( "Homo sapiens", "Pan", "Sapporo virus", "sapiens" );
        // case-INSENSITIVE "sap": prefix (starts-with) matches first, then other substring matches, alphabetical
        // within each group -> [Sapporo virus, sapiens] ++ [Homo sapiens]
        ck( SearchValueAutocomplete.filterValues( all, "sap", false )
                .equals( Arrays.asList( "Sapporo virus", "sapiens", "Homo sapiens" ) ),
            "case-insensitive filter should rank prefix matches first" );
        // case-SENSITIVE "sap": "Sapporo" (capital S) no longer matches at all; "sapiens" is the only prefix
        ck( SearchValueAutocomplete.filterValues( all, "sap", true )
                .equals( Arrays.asList( "sapiens", "Homo sapiens" ) ),
            "case-sensitive filter should honour case" );
        // empty query -> the whole browse list, returned as-is (no copy)
        ck( SearchValueAutocomplete.filterValues( all, "", false ) == all, "empty query returns the browse list" );
        ck( SearchValueAutocomplete.filterValues( all, "  ", false ) == all, "a blank query returns the browse list" );
        ck( SearchValueAutocomplete.filterValues( all, "zzz", false ).isEmpty(), "no match -> empty" );
    }

    private static void capAndHint() {
        final List<String> five = sorted( "a1", "a2", "a3", "a4", "a5" );
        // more matches than the cap -> the first N + a non-selectable hint row
        final List<String> capped = SearchValueAutocomplete.displayModel( five, "", false, 3 );
        ck( capped.size() == 4, "display model should be capped at 3 + a hint row" );
        ck( capped.subList( 0, 3 ).equals( Arrays.asList( "a1", "a2", "a3" ) ), "the cap keeps the first N matches" );
        ck( SearchValueAutocomplete.HINT_TEXT.equals( capped.get( 3 ) ), "a truncated list ends with the hint row" );
        // within the cap -> no hint row
        final List<String> full = SearchValueAutocomplete.displayModel( five, "", false, 10 );
        ck( full.equals( five ), "an un-truncated list is shown verbatim" );
        ck( !full.contains( SearchValueAutocomplete.HINT_TEXT ), "an un-truncated list has no hint row" );
    }

    // ---- distinct values + the node-type categorical field --------------------------------------------------

    private static void distinctValuesAndNodeType() {
        // root -> ( mid -> (a, b), c ); a/b/c named leaves, mid + root unnamed
        final PhylogenyNode a = named( "a" );
        final PhylogenyNode b = named( "b" );
        final PhylogenyNode c = named( "c" );
        final PhylogenyNode mid = new PhylogenyNode();
        mid.addAsChild( a );
        mid.addAsChild( b );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( mid );
        root.addAsChild( c );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();

        // node-name distinct values: only the three named leaves (mid/root are unnamed -> skipped), sorted
        ck( SearchField.distinctValues( phy, SearchField.ofNdf( NDF.NodeName ) )
                .equals( Arrays.asList( "a", "b", "c" ) ), "node-name distinct values should be the named leaves" );
        // node-type distinct values: the categorical set present in the tree, sorted
        ck( SearchField.distinctValues( phy, SearchField.nodeType() )
                .equals( Arrays.asList( "internal", "leaf", "root" ) ),
            "node-type distinct values should be internal/leaf/root" );
        // "Any text" and numeric fields are not meaningful pick lists -> empty
        ck( SearchField.distinctValues( phy, SearchField.anyText() ).isEmpty(),
            "\"Any text\" should yield no autocomplete candidates" );
        ck( SearchField.distinctValues( phy, SearchField.cladeSize() ).isEmpty(),
            "a numeric field should yield no autocomplete candidates" );
        // a molecular sequence is (near-)unique per tip and can be huge -> no autocomplete pick-list even when present
        final Sequence seq = new Sequence();
        seq.setMolecularSequence( "ACGTACGTACGTACGTACGT" );
        a.getNodeData().setSequence( seq );
        ck( SearchField.distinctValues( phy, SearchField.ofNdf( NDF.MolecularSequence ) ).isEmpty(),
            "molecular sequence should have no autocomplete pick-list (even with a sequence present)" );

        // node type is offered in the per-tree field list, is a text field, and matches
        ck( byLabel( SearchField.availableFields( phy ), "Structure: Node Type" ) != null,
            "node type should be offered by availableFields" );
        final SearchField nt = SearchField.nodeType();
        ck( !nt.isNumeric(), "node type is a text (categorical) field" );
        ck( SearchField.nodeTypeLabel( a ).equals( "leaf" ), "a tip is a leaf" );
        ck( SearchField.nodeTypeLabel( mid ).equals( "internal" ), "a non-root internal node is internal" );
        ck( SearchField.nodeTypeLabel( root ).equals( "root" ), "the root is root" );
        ck( SearchMatcher.matchesPositive( new SearchSpec( nt, SearchMode.WHOLE_WORD, "leaf" ), a ),
            "node type 'leaf' should match a tip" );
        ck( !SearchMatcher.matchesPositive( new SearchSpec( nt, SearchMode.WHOLE_WORD, "leaf" ), root ),
            "node type 'leaf' should not match the root" );
        ck( SearchMatcher.matchesPositive( new SearchSpec( nt, SearchMode.WHOLE_WORD, "internal" ), mid ),
            "node type 'internal' should match an internal node" );
    }

    // ---- accept (fills the field + runs the search); N/A supplier -> no popup -------------------------------

    private static void acceptFillsAndSearches() {
        final JTextField tf = new JTextField();
        final boolean[] ran = { false };
        final SearchValueAutocomplete ac = new SearchValueAutocomplete( tf,
                () -> sorted( "Africa", "Asia", "Europe" ), () -> false, () -> ran[ 0 ] = true );
        tf.setText( "a" ); // the field is unfocused, so this must NOT trigger the popup (guarded document listener)
        final List<String> shown = ac.displayModelForTest();
        ck( shown.contains( "Africa" ) && shown.contains( "Asia" ) && !shown.contains( "Europe" ),
            "typing 'a' should filter to the a-matches: " + shown );
        ac.acceptValueForTest( "Europe" );
        ck( tf.getText().equals( "Europe" ), "accepting a suggestion fills the field with the exact value" );
        ck( ran[ 0 ], "accepting a suggestion runs the search" );
        // an empty value supplier (N/A field: Any-text / numeric / regex) offers nothing -> no popup
        final SearchValueAutocomplete none = new SearchValueAutocomplete( new JTextField(),
                () -> Collections.<String>emptyList(), () -> false, () -> {} );
        ck( none.displayModelForTest().isEmpty(), "an empty value supplier yields no suggestions" );
    }

    // ---- helpers --------------------------------------------------------------------------------------------

    private static List<String> sorted( final String... values ) {
        return new ArrayList<>( new TreeSet<>( Arrays.asList( values ) ) );
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

    private static void ck( final boolean cond, final String msg ) {
        if ( !cond ) {
            throw new AssertionError( msg );
        }
    }

    private SearchAutocompleteTest() {
    }
}
