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

import java.awt.GraphicsEnvironment;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods.NDF;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;

/**
 * The search value box's suggestion list, a PORT of the Archaeopteryx.js list (b7e1bc8): only the term being typed
 * is matched and completed (',' is OR, '+' is AND), values are matched the way the box's mode will match them and
 * always ignoring case, at most ten rows in the values' sorted order with the matched part in bold accent, then
 * "N more — keep typing"; nothing when nothing matches or the one match is what is typed; the arrows open the list and
 * wrap, Enter picks, Escape and Tab close. Christian, 2026-09-13: search is a strong point of both programs, and the
 * two must behave exactly alike.
 */
public final class SearchAutocompleteTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SearchAutocomplete: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            termOfAQuery();
            matchingFollowsTheMode();
            tenRowsThenMore();
            nothingToOffer();
            matchedPartIsMarked();
            pickReplacesOnlyTheTerm();
            valuesSortedAndCapped();
            distinctValuesAndNodeType();
            acceptFillsAndSearches();
            if ( !GraphicsEnvironment.isHeadless() ) {
                keyboardOnARealBox();
            }
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

    // ---- the pure rules ---------------------------------------------------------------------------------------

    private static void termOfAQuery() {
        ck( SearchValueAutocomplete.termStart( "abc" ) == 0, "a plain query is one term" );
        ck( SearchValueAutocomplete.termStart( "cat,do" ) == 4, "the term starts after the last ','" );
        ck( SearchValueAutocomplete.termStart( "a+b" ) == 2, "the term starts after the last '+'" );
        ck( SearchValueAutocomplete.termStart( "a,b+c" ) == 4, "the LAST separator of either kind wins" );
        ck( SearchValueAutocomplete.termStart( "x," ) == 2, "a trailing separator leaves an empty term" );
        ck( SearchValueAutocomplete.currentTerm( "cat, do " ).equals( "do" ), "the term is trimmed" );
        ck( SearchValueAutocomplete.currentTerm( "" ).isEmpty() && SearchValueAutocomplete.currentTerm( null ).isEmpty(),
            "an empty box has an empty term" );
    }

    private static void matchingFollowsTheMode() {
        final List<String> all = Arrays.asList( "Homo sapiens", "Pan", "Sapporo virus", "sapiens" );
        // contains (and whole word / exact): a substring, case ignored, in the values' own order -- no re-ranking
        ck( SearchValueAutocomplete.matches( all, "sap", SearchMode.CONTAINS )
                .equals( Arrays.asList( "Homo sapiens", "Sapporo virus", "sapiens" ) ),
            "contains should admit every substring match, in the list's order" );
        ck( SearchValueAutocomplete.matches( all, "SAP", SearchMode.WHOLE_WORD )
                .equals( Arrays.asList( "Homo sapiens", "Sapporo virus", "sapiens" ) ),
            "matching always ignores case" );
        ck( SearchValueAutocomplete.matches( all, "sap", SearchMode.STARTS_WITH )
                .equals( Arrays.asList( "Sapporo virus", "sapiens" ) ), "starts-with should admit prefixes only" );
        ck( SearchValueAutocomplete.matches( all, "ens", SearchMode.ENDS_WITH )
                .equals( Arrays.asList( "Homo sapiens", "sapiens" ) ), "ends-with should admit suffixes only" );
        ck( SearchValueAutocomplete.matches( all, "", SearchMode.STARTS_WITH ).equals( all ),
            "an empty term admits everything" );
        ck( SearchValueAutocomplete.matches( all, "zzz", SearchMode.CONTAINS ).isEmpty(), "no match -> empty" );
    }

    private static void tenRowsThenMore() {
        final List<String> many = new ArrayList<>();
        for ( int i = 1; i <= 25; ++i ) {
            many.add( String.format( "v%02d", i ) );
        }
        final SearchValueAutocomplete.Model m = SearchValueAutocomplete.model( many, "v", SearchMode.CONTAINS );
        ck( m._rows.size() == SearchValueAutocomplete.MAX_ROWS, "at most ten rows are shown" );
        ck( m._rows.equals( many.subList( 0, 10 ) ), "the ten are the first matches in order" );
        ck( m._more == 15, "the rest are counted: 25 - 10 = 15, got " + m._more );
        ck( SearchValueAutocomplete.moreText( 15 ).equals( "15 more — keep typing" ),
            "the last row says how many more, with an em dash" );
        final SearchValueAutocomplete.Model few = SearchValueAutocomplete.model( many, "v1", SearchMode.CONTAINS );
        ck( ( few._rows.size() == 10 ) && ( few._more == 0 ), "v1 matches exactly v10..v19: ten rows and none more" );
        final SearchValueAutocomplete.Model exact = SearchValueAutocomplete.model( many, "v2", SearchMode.CONTAINS );
        ck( ( exact._rows.size() == 6 ) && ( exact._more == 0 ), "within ten rows there is no more row" );
        // the whole browse list when the term is empty
        final SearchValueAutocomplete.Model browse = SearchValueAutocomplete.model( many, "", SearchMode.CONTAINS );
        ck( ( browse._rows.size() == 10 ) && ( browse._more == 15 ), "an empty term browses the first ten" );
    }

    private static void nothingToOffer() {
        final List<String> all = Arrays.asList( "Africa", "Asia", "Europe" );
        ck( SearchValueAutocomplete.model( all, "zzz", SearchMode.CONTAINS ).isEmpty(), "no match -> nothing" );
        ck( SearchValueAutocomplete.model( all, "europe", SearchMode.CONTAINS ).isEmpty(),
            "the one match being what is typed (case ignored) -> nothing" );
        ck( SearchValueAutocomplete.model( all, "Asia, EUROPE", SearchMode.CONTAINS ).isEmpty(),
            "...and that is judged on the term being typed, not the whole box" );
        ck( SearchValueAutocomplete.model( all, "a", SearchMode.CONTAINS )._rows.equals( Arrays.asList( "Africa", "Asia" ) ),
            "two matches are offered even if one equals the term" );
        ck( SearchValueAutocomplete.model( Collections.<String>emptyList(), "", SearchMode.CONTAINS ).isEmpty(),
            "no values -> nothing" );
    }

    private static void matchedPartIsMarked() {
        ck( Arrays.equals( SearchValueAutocomplete.matchSpan( "Homo sapiens", "SAP" ), new int[] { 5, 3 } ),
            "the matched part is the term's first occurrence, case ignored" );
        ck( Arrays.equals( SearchValueAutocomplete.matchSpan( "sapiens sapiens", "sap" ), new int[] { 0, 3 } ),
            "the FIRST occurrence, even when a later one exists" );
        ck( SearchValueAutocomplete.matchSpan( "Pan", "sap" ) == null, "no occurrence -> no mark" );
        ck( SearchValueAutocomplete.matchSpan( "Pan", "" ) == null, "an empty term marks nothing" );
        final String html = SearchValueAutocomplete.rowHtml( "Homo sapiens", "sap", "#2675bf" );
        ck( html.equals( "<html>Homo <b><font color=\"#2675bf\">sap</font></b>iens" ),
            "the row marks the matched part in bold accent, got " + html );
        ck( SearchValueAutocomplete.rowHtml( "a<b>&c", "b", "#000000" )
                .equals( "<html>a&lt;<b><font color=\"#000000\">b</font></b>&gt;&amp;c" ), "values are HTML-escaped" );
        ck( SearchValueAutocomplete.rowHtml( "Pan", "sap", "#000000" ).equals( "<html>Pan" ),
            "an unmarked row is the plain escaped value" );
    }

    private static void pickReplacesOnlyTheTerm() {
        ck( SearchValueAutocomplete.pickInto( "cat,do", "dog" ).equals( "cat,dog" ), "a pick keeps the OR before it" );
        ck( SearchValueAutocomplete.pickInto( "cat+ do", "dog" ).equals( "cat+dog" ),
            "a pick keeps the AND before it and drops the term's leading space" );
        ck( SearchValueAutocomplete.pickInto( "do", "dog" ).equals( "dog" ), "a single term is replaced whole" );
        ck( SearchValueAutocomplete.pickInto( "", "dog" ).equals( "dog" ), "an empty box takes the value" );
        ck( SearchValueAutocomplete.pickInto( "a,b,", "c" ).equals( "a,b,c" ), "a trailing separator is kept" );
    }

    private static void valuesSortedAndCapped() {
        final List<String> l = Arrays.asList( "cherry", "Banana", "apple", "Apple" );
        Collections.sort( l, SearchField.SUGGESTION_ORDER );
        ck( l.equals( Arrays.asList( "apple", "Apple", "Banana", "cherry" ) ),
            "values sort case-insensitively first, lowercase before uppercase on a tie (like localeCompare): " + l );
        // ...and the VALUE LIST uses that order (code-unit order would put "Banana" first): a sabotage that sorted the
        // list plainly survived while only the comparator was pinned
        final PhylogenyNode r = new PhylogenyNode();
        for ( final String n : new String[] { "cherry", "Banana", "apple" } ) {
            r.addAsChild( named( n ) );
        }
        final Phylogeny small = new Phylogeny();
        small.setRoot( r );
        small.setRooted( true );
        small.externalNodesHaveChanged();
        ck( SearchField.distinctValues( small, SearchField.ofNdf( NDF.NodeName ) )
                .equals( Arrays.asList( "apple", "Banana", "cherry" ) ), "the value list is sorted case-insensitively" );
        // the cap: a tree with more distinct names than the cap offers the first CAP in sorted order
        final PhylogenyNode root = new PhylogenyNode();
        for ( int i = 0; i < SearchField.AUTOCOMPLETE_VALUE_CAP + 5; ++i ) {
            root.addAsChild( named( String.format( "n%05d", i ) ) );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        final List<String> vals = SearchField.distinctValues( phy, SearchField.ofNdf( NDF.NodeName ) );
        ck( vals.size() == SearchField.AUTOCOMPLETE_VALUE_CAP, "the value list is capped at " + SearchField.AUTOCOMPLETE_VALUE_CAP
                + ", got " + vals.size() );
        ck( vals.get( 0 ).equals( "n00000" ) && vals.get( vals.size() - 1 ).equals( "n01999" ),
            "the cap keeps the first values in sorted order" );
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

    // ---- accept (replaces the term + runs the search); N/A supplier -> no popup -------------------------------

    private static void acceptFillsAndSearches() {
        final JTextField tf = new JTextField();
        final boolean[] ran = { false };
        final SearchValueAutocomplete ac = new SearchValueAutocomplete( tf,
                () -> Arrays.asList( "Africa", "Asia", "Europe" ), () -> SearchMode.CONTAINS, () -> ran[ 0 ] = true );
        tf.setText( "Europe,a" ); // the field is unfocused, so this must NOT trigger the popup (guarded document listener)
        final SearchValueAutocomplete.Model shown = ac.modelForTest();
        ck( shown._rows.equals( Arrays.asList( "Africa", "Asia" ) ),
            "typing 'a' after an OR should offer the a-matches of the TERM: " + shown._rows );
        ac.acceptValueForTest( "Asia" );
        ck( tf.getText().equals( "Europe,Asia" ), "accepting a suggestion replaces only the term, got " + tf.getText() );
        ck( ran[ 0 ], "accepting a suggestion runs the search" );
        // an empty value supplier (N/A field: Any-text / numeric / regex) offers nothing -> no popup
        final SearchValueAutocomplete none = new SearchValueAutocomplete( new JTextField(),
                () -> Collections.<String>emptyList(), () -> SearchMode.CONTAINS, () -> {} );
        ck( none.modelForTest().isEmpty(), "an empty value supplier yields no suggestions" );
    }

    // ---- the keys, on a real (shown) box ----------------------------------------------------------------------

    private static void keyboardOnARealBox() throws Exception {
        final JFrame[] frame = new JFrame[ 1 ];
        final Throwable[] failure = { null };
        SwingUtilities.invokeAndWait( () -> {
            try {
                final JTextField tf = new JTextField( 20 );
                final SearchValueAutocomplete ac = new SearchValueAutocomplete( tf,
                        () -> Arrays.asList( "Africa", "Antarctica", "Asia", "Europe" ), () -> SearchMode.CONTAINS,
                        () -> {} );
                frame[ 0 ] = new JFrame( "suggest" );
                final JPanel p = new JPanel();
                p.add( tf );
                frame[ 0 ].getContentPane().add( p );
                frame[ 0 ].pack();
                frame[ 0 ].setLocation( 40, 40 );
                frame[ 0 ].setVisible( true );
                tf.setText( "Europe,a" );
                // Down opens the closed list and highlights the first row
                key( tf, KeyEvent.VK_DOWN );
                ck( ac.isShowingForTest(), "Down should open the list" );
                ck( ac.shownRowsForTest().equals( Arrays.asList( "Africa", "Antarctica", "Asia" ) ),
                    "the list should hold the term's matches, got " + ac.shownRowsForTest() );
                ck( ac.selectedIndexForTest() == 0, "Down on a fresh list highlights the first row" );
                key( tf, KeyEvent.VK_DOWN );
                key( tf, KeyEvent.VK_DOWN );
                ck( ac.selectedIndexForTest() == 2, "Down moves down" );
                key( tf, KeyEvent.VK_DOWN );
                ck( ac.selectedIndexForTest() == 0, "Down past the last row wraps to the first" );
                key( tf, KeyEvent.VK_UP );
                ck( ac.selectedIndexForTest() == 2, "Up from the first row wraps to the last" );
                // Escape closes the list only; the text is untouched
                key( tf, KeyEvent.VK_ESCAPE );
                ck( !ac.isShowingForTest() && tf.getText().equals( "Europe,a" ), "Escape closes the list and keeps the text" );
                // Up on a closed list opens it and highlights the last row
                key( tf, KeyEvent.VK_UP );
                ck( ac.isShowingForTest() && ( ac.selectedIndexForTest() == 2 ), "Up on a fresh list highlights the last row" );
                // Tab is a focus-traversal key: Swing moves the focus before any key listener sees it, so what Tab does
                // to the list is what LEAVING THE BOX does -- the focus listener calls endSession(), which is exercised
                // here directly (a synthetic FOCUS_LOST is swallowed by the focus manager when the box is not the real
                // focus owner, and real focus is not something a test can rely on): it closes, without picking
                ac.endSession();
                ck( !ac.isShowingForTest() && tf.getText().equals( "Europe,a" ),
                    "leaving the box (Tab) closes the list without picking" );
                // Enter on a highlighted row replaces only the term
                key( tf, KeyEvent.VK_DOWN );
                key( tf, KeyEvent.VK_DOWN );
                ck( ac.selectedIndexForTest() == 1, "two Downs highlight the second row" );
                key( tf, KeyEvent.VK_ENTER );
                ck( tf.getText().equals( "Europe,Antarctica" ), "Enter picks the highlighted row into the term, got "
                        + tf.getText() );
                ck( !ac.isShowingForTest(), "a pick closes the list" );
                // more than ten matches: ten rows and the "N more" row, which the arrows never land on
                final JTextField tf2 = new JTextField( 20 );
                final List<String> many = new ArrayList<>();
                for ( int i = 1; i <= 12; ++i ) {
                    many.add( String.format( "v%02d", i ) );
                }
                final SearchValueAutocomplete ac2 = new SearchValueAutocomplete( tf2, () -> many, () -> SearchMode.CONTAINS,
                                                                                 () -> {} );
                p.add( tf2 );
                frame[ 0 ].pack();
                tf2.setText( "v" );
                key( tf2, KeyEvent.VK_DOWN );
                ck( ac2.shownRowsForTest().size() == 11, "ten rows plus the more row, got " + ac2.shownRowsForTest() );
                ck( ac2.shownRowsForTest().get( 10 ).equals( "2 more — keep typing" ),
                    "the last row counts the rest, got " + ac2.shownRowsForTest().get( 10 ) );
                key( tf2, KeyEvent.VK_UP ); // from the first row
                ck( ac2.selectedIndexForTest() == 9, "Up wraps to the last VALUE row, never the more row, got "
                        + ac2.selectedIndexForTest() );
            }
            catch ( final Throwable t ) {
                failure[ 0 ] = t;
            }
            finally {
                if ( frame[ 0 ] != null ) {
                    frame[ 0 ].dispose();
                }
            }
        } );
        if ( failure[ 0 ] instanceof AssertionError ) {
            throw (AssertionError) failure[ 0 ];
        }
        if ( failure[ 0 ] != null ) {
            throw new RuntimeException( failure[ 0 ] );
        }
    }

    private static void key( final JTextField tf, final int code ) {
        tf.dispatchEvent( new KeyEvent( tf, KeyEvent.KEY_PRESSED, System.currentTimeMillis(), 0, code,
                                        KeyEvent.CHAR_UNDEFINED ) );
    }

    // ---- helpers --------------------------------------------------------------------------------------------

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
