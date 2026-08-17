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
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headful coverage for the redesigned search UI. Phase 2: the two search boxes each carry a field + match-mode
 * selector, defaulting to "Any text field" + "contains"; any-text vs field-scoped, the string modes, Match Case,
 * Inverse, A/B independence. Phase 3: the field selector is tailored to the loaded tree (only present fields + one
 * entry per custom property + numeric built-ins), the mode selector switches STRING&harr;NUMERIC by field type,
 * the numeric operators and the on-demand range field work, and per-property matching is key-scoped. A green no-op
 * when headless.
 */
public final class SearchToolTest {

    private static final String ANY  = "Any text field";
    private static final String NODE = "Node name";
    private static final String SCI  = "Taxonomy: scientific name";
    private static final String BL   = "Branch length";
    private static final String HOST = "aptx:host";
    private static final String READS = "data:reads";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SearchTool: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        final PhylogenyNode[] leaves = new PhylogenyNode[ 4 ];
        try {
            final Phylogeny phy = buildTree( leaves );
            final Configuration conf = new Configuration();
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "search" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    final ControlPanel cp = tp.getControlPanel();
                    cp.rebuildSearchFields( true ); // ensure the field lists reflect the loaded tree
                    final PhylogenyNode a = leaves[ 0 ]; // alpha / Homo sapiens / bl 0.5 / host human / reads 1000
                    final PhylogenyNode b = leaves[ 1 ]; // beta / Felis catus / bl 1.5 / host cat / reads 200
                    final PhylogenyNode c = leaves[ 2 ]; // gamma_kinase / Mus musculus / bl 2.5 / host mouse / reads 50
                    final PhylogenyNode mid = leaves[ 3 ]; // internal: bl 0.2, but NO node data

                    // defaults + dynamic field list
                    ck( ok, cp.getSearchFieldForTest( true ).kind() == SearchField.Kind.ANY_TEXT,
                        "the default field should be Any text" );
                    ck( ok, cp.getSearchModeForTest( true ) == SearchMode.CONTAINS, "the default mode should be contains" );
                    ck( ok, cp.searchModeCountForTest() == 5, "a string field should offer 5 modes" );
                    final List<String> fields = cp.searchFieldLabelsForTest( true );
                    for ( final String want : new String[] { ANY, NODE, SCI, BL, HOST, READS } ) {
                        ck( ok, fields.contains( want ), "field selector should list \"" + want + "\", has " + fields );
                    }
                    ck( ok, !fields.contains( "Gene name" ), "field selector should NOT list absent fields (Gene name)" );

                    // --- string: any-text vs field-scoped ---
                    ck( ok, same( runA( cp, tp, ANY, SearchMode.CONTAINS, false, false, "catus", "" ), b ),
                        "Any-text 'catus' should find Felis catus" );
                    ck( ok, same( runA( cp, tp, ANY, SearchMode.CONTAINS, false, false, "alpha", "" ), a ),
                        "Any-text 'alpha' should find alpha by node name" );
                    ck( ok, runA( cp, tp, NODE, SearchMode.CONTAINS, false, false, "Homo", "" ).isEmpty(),
                        "Node-name 'Homo' should find nothing" );
                    ck( ok, same( runA( cp, tp, SCI, SearchMode.CONTAINS, false, false, "Homo", "" ), a ),
                        "Scientific-name 'Homo' should find alpha" );
                    ck( ok, same( runA( cp, tp, NODE, SearchMode.STARTS_WITH, false, false, "gam", "" ), c ),
                        "starts-with 'gam' should find gamma_kinase" );
                    ck( ok, same( runA( cp, tp, NODE, SearchMode.ENDS_WITH, false, false, "kinase", "" ), c ),
                        "ends-with 'kinase' should find gamma_kinase" );
                    ck( ok, same( runA( cp, tp, SCI, SearchMode.CONTAINS, false, false, "homo", "" ), a ),
                        "case-insensitive 'homo' should find Homo sapiens" );
                    ck( ok, runA( cp, tp, SCI, SearchMode.CONTAINS, true, false, "homo", "" ).isEmpty(),
                        "case-sensitive 'homo' should find nothing" );
                    final Set<Long> inv = runA( cp, tp, NODE, SearchMode.CONTAINS, false, true, "alpha", "" );
                    ck( ok, !inv.contains( a.getId() ) && inv.contains( b.getId() ) && inv.contains( c.getId() )
                            && ( inv.size() == 2 ), "inverse of node-name 'alpha' should be {beta, gamma_kinase}" );
                    ck( ok, runA( cp, tp, ANY, SearchMode.CONTAINS, false, true, ",", "" ).isEmpty(),
                        "a separator-only query with Inverse should reset, not select the whole tree" );

                    // --- per-property (key-scoped) ---
                    ck( ok, same( runA( cp, tp, HOST, SearchMode.CONTAINS, false, false, "cat", "" ), b ),
                        "aptx:host 'cat' should find beta" );
                    ck( ok, runA( cp, tp, HOST, SearchMode.CONTAINS, false, false, "1000", "" ).isEmpty(),
                        "aptx:host must not match a value living in data:reads (key-scoped)" );

                    // --- numeric: branch length ---
                    ck( ok, same( runA( cp, tp, BL, SearchMode.GT, false, false, "1", "" ), b, c ),
                        "branch length > 1 should be {beta, gamma_kinase}" );
                    ck( ok, same( runA( cp, tp, BL, SearchMode.LT, false, false, "1", "" ), a, mid ),
                        "branch length < 1 should be {alpha, mid}" );
                    ck( ok, same( runA( cp, tp, BL, SearchMode.EQ, false, false, "0.5", "" ), a ),
                        "branch length = 0.5 should be {alpha}" );
                    ck( ok, same( runA( cp, tp, BL, SearchMode.RANGE, false, false, "1", "2" ), b ),
                        "branch length in [1,2] should be {beta}" );
                    // inverse of a numeric query complements over nodes that CARRY the field (branch length), not
                    // just data-bearing nodes: the internal 'mid' has a branch length (0.2) but no node data.
                    final Set<Long> bl_inv = runA( cp, tp, BL, SearchMode.GT, false, true, "0.5", "" );
                    ck( ok, bl_inv.contains( a.getId() ) && bl_inv.contains( mid.getId() ) && ( bl_inv.size() == 2 ),
                        "inverse of branch length > 0.5 should be {alpha, mid} (mid has a branch length, no data)" );
                    // a non-parseable numeric operand with Inverse on must reset, not select the whole tree
                    ck( ok, runA( cp, tp, BL, SearchMode.GT, false, true, "kinase", "" ).isEmpty(),
                        "a non-numeric operand on a numeric field (Inverse on) should reset, not select all" );

                    // --- numeric: a decimal property ---
                    ck( ok, same( runA( cp, tp, READS, SearchMode.GT, false, false, "500", "" ), a ),
                        "data:reads > 500 should be {alpha}" );
                    ck( ok, same( runA( cp, tp, READS, SearchMode.LT, false, false, "100", "" ), c ),
                        "data:reads < 100 should be {gamma_kinase}" );

                    // --- mode-set switching + range field visibility ---
                    cp.setSearchFieldByLabelForTest( true, BL );
                    ck( ok, cp.searchModeCountForTest() == 7, "a numeric field should offer 7 operators, got "
                            + cp.searchModeCountForTest() );
                    ck( ok, cp.getSearchModeForTest( true ).isNumeric(), "a numeric field should make the mode combo numeric" );
                    cp.setSearchModeForTest( true, SearchMode.RANGE );
                    ck( ok, cp.rangeFieldVisibleForTest( true ), "range mode should reveal the range field" );
                    cp.setSearchModeForTest( true, SearchMode.GT );
                    ck( ok, !cp.rangeFieldVisibleForTest( true ), "a non-range mode should hide the range field" );
                    cp.setSearchFieldByLabelForTest( true, ANY );
                    ck( ok, cp.searchModeCountForTest() == 5, "switching back to a string field should restore 5 modes" );
                    ck( ok, !cp.getSearchModeForTest( true ).isNumeric(), "a string field should make the mode combo string" );

                    // --- a custom property that flips numeric -> string must repopulate (label unchanged) ---
                    ck( ok, cp.searchFieldLabelsForTest( true ).contains( "data:n" ), "data:n should be a field" );
                    cp.setSearchFieldByLabelForTest( true, "data:n" );
                    ck( ok, cp.getSearchModeForTest( true ).isNumeric(),
                        "data:n starts numeric (datatype-less, all values parse)" );
                    setPropValue( a, "data:n", "x" ); // its values are no longer numeric
                    setPropValue( b, "data:n", "y" );
                    setPropValue( c, "data:n", "z" );
                    cp.rebuildSearchFields( true );
                    cp.setSearchFieldByLabelForTest( true, "data:n" );
                    ck( ok, !cp.getSearchModeForTest( true ).isNumeric(),
                        "data:n should become a STRING field once its values stop being numeric" );

                    // --- A/B independence ---
                    cp.setSearchInverseForTest( false );
                    cp.setSearchCaseSensitiveForTest( false );
                    final Set<Long> found_a = runA( cp, tp, ANY, SearchMode.CONTAINS, false, false, "catus", "" );
                    final Set<Long> found_b = runB( cp, tp, ANY, SearchMode.CONTAINS, "musculus" );
                    ck( ok, same( found_a, b ) && same( found_b, c ),
                        "A ('catus') and B ('musculus') should hold independent found sets" );
                    ck( ok, same( nonNull( tp.getFoundNodes0() ), b ), "running B must not disturb A's found set" );

                    // --- remember-last (in-session): a chosen mode survives a field-KIND switch ---
                    cp.userSelectFieldForTest( true, NODE );
                    cp.userSelectModeForTest( true, SearchMode.STARTS_WITH );
                    cp.userSelectFieldForTest( true, BL );   // numeric detour repopulates the mode combo
                    cp.userSelectFieldForTest( true, NODE ); // back to a string field
                    ck( ok, cp.getSearchModeForTest( true ) == SearchMode.STARTS_WITH,
                        "the string mode should be remembered across a numeric-field detour" );
                    cp.userSelectFieldForTest( true, BL );
                    cp.userSelectModeForTest( true, SearchMode.GT );
                    cp.userSelectFieldForTest( true, NODE );
                    cp.userSelectFieldForTest( true, BL );
                    ck( ok, cp.getSearchModeForTest( true ) == SearchMode.GT,
                        "the numeric mode should be remembered across a string-field detour" );
                    // Reset clears the remembered modes back to the defaults
                    cp.userSelectFieldForTest( true, NODE );
                    cp.userSelectModeForTest( true, SearchMode.ENDS_WITH );
                    cp.resyncFromOptions();
                    cp.userSelectFieldForTest( true, BL );
                    cp.userSelectFieldForTest( true, NODE );
                    ck( ok, cp.getSearchModeForTest( true ) == SearchMode.CONTAINS,
                        "Reset should clear the remembered mode (contains, not the pre-reset ends-with)" );

                    ( (javax.swing.JFrame) mf[ 0 ] ).dispose();
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
                }
            } );
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
        return ok[ 0 ];
    }

    // ---- helpers --------------------------------------------------------------------------------------------

    private static Set<Long> runA( final ControlPanel cp, final TreePanel tp, final String field_label,
                                   final SearchMode mode, final boolean case_sensitive, final boolean inverse,
                                   final String value, final String range_high ) {
        cp.setSearchFieldByLabelForTest( true, field_label );
        cp.setSearchModeForTest( true, mode );
        cp.setSearchCaseSensitiveForTest( case_sensitive );
        cp.setSearchInverseForTest( inverse );
        cp.getSearchTextField0().setText( value );
        cp.setRangeHighForTest( true, range_high );
        cp.search0();
        return nonNull( tp.getFoundNodes0() );
    }

    private static Set<Long> runB( final ControlPanel cp, final TreePanel tp, final String field_label,
                                   final SearchMode mode, final String value ) {
        cp.setSearchFieldByLabelForTest( false, field_label );
        cp.setSearchModeForTest( false, mode );
        cp.getSearchTextField1().setText( value );
        cp.search1();
        return nonNull( tp.getFoundNodes1() );
    }

    private static Set<Long> nonNull( final Set<Long> s ) {
        return ( s == null ) ? Collections.emptySet() : s;
    }

    private static boolean same( final Set<Long> got, final PhylogenyNode... expected ) {
        final Set<Long> want = new HashSet<>();
        for ( final PhylogenyNode n : expected ) {
            want.add( n.getId() );
        }
        return got.equals( want );
    }

    private static Phylogeny buildTree( final PhylogenyNode[] out ) {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode mid = new PhylogenyNode();
        mid.setDistanceToParent( 0.2 ); // a branch length, but no node data
        final PhylogenyNode a = leaf( "alpha", "Homo sapiens", 0.5, "human", "1000", "1" );
        final PhylogenyNode b = leaf( "beta", "Felis catus", 1.5, "cat", "200", "2" );
        final PhylogenyNode c = leaf( "gamma_kinase", "Mus musculus", 2.5, "mouse", "50", "3" );
        mid.addAsChild( a );
        mid.addAsChild( b );
        root.addAsChild( mid );
        root.addAsChild( c );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.setName( "searchtool" );
        phy.externalNodesHaveChanged();
        out[ 0 ] = a;
        out[ 1 ] = b;
        out[ 2 ] = c;
        out[ 3 ] = mid;
        return phy;
    }

    private static PhylogenyNode leaf( final String name, final String sci, final double bl, final String host,
                                       final String reads, final String n_val ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( bl );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        n.getNodeData().setTaxonomy( t );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( HOST, host, "", "xsd:string", AppliesTo.NODE ) );
        pl.addProperty( new Property( READS, reads, "", "xsd:decimal", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:n", n_val, "", "", AppliesTo.NODE ) ); // datatype-less numeric
        n.getNodeData().setProperties( pl );
        return n;
    }

    private static void setPropValue( final PhylogenyNode n, final String ref, final String val ) {
        for ( final Property p : n.getNodeData().getProperties().getProperties() ) {
            if ( ref.equals( p.getRef() ) ) {
                p.setValue( val );
            }
        }
    }

    private static void ck( final boolean[] ok, final boolean cond, final String msg ) {
        if ( !cond ) {
            System.out.println( "  [SearchToolTest] " + msg );
            ok[ 0 ] = false;
        }
    }

    private SearchToolTest() {
    }
}
