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
import java.util.Set;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods.NDF;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headful coverage for the redesigned search UI (Phase 2): the two search boxes each carry a field selector and
 * a match-mode selector, defaulting to "Any text field" + "contains" (today's behaviour). Drives real searches
 * through {@link ControlPanel#search0()}/{@link ControlPanel#search1()} and asserts the found sets: any-text vs
 * field-scoped search, the string modes, Match Case, Inverse, and A/B independence. A green no-op when headless.
 */
public final class SearchToolTest {

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
        final PhylogenyNode[] leaves = new PhylogenyNode[ 3 ];
        try {
            final Phylogeny phy = buildTree( leaves );
            final Configuration conf = new Configuration();
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "search" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    final ControlPanel cp = tp.getControlPanel();
                    final PhylogenyNode a = leaves[ 0 ]; // "alpha" / Homo sapiens
                    final PhylogenyNode b = leaves[ 1 ]; // "beta" / Felis catus
                    final PhylogenyNode c = leaves[ 2 ]; // "gamma_kinase" / Mus musculus

                    // the combos exist with the expected contents + defaults
                    ck( ok, cp.searchFieldCountForTest() == 18, "field combo should list 18 fields, got "
                            + cp.searchFieldCountForTest() );
                    ck( ok, cp.searchModeCountForTest() == 5, "mode combo should list 5 string modes, got "
                            + cp.searchModeCountForTest() );
                    ck( ok, cp.getSearchFieldForTest( true ).kind() == SearchField.Kind.ANY_TEXT,
                        "the default field should be Any text" );
                    ck( ok, cp.getSearchModeForTest( true ) == SearchMode.CONTAINS, "the default mode should be contains" );

                    // Any-text searches every text field (node name AND scientific name)
                    ck( ok, same( runA( cp, tp, null, SearchMode.CONTAINS, false, false, "catus" ), b ),
                        "Any-text 'catus' should find the Felis catus leaf" );
                    ck( ok, same( runA( cp, tp, null, SearchMode.CONTAINS, false, false, "alpha" ), a ),
                        "Any-text 'alpha' should find the alpha leaf by node name" );

                    // field-scoped: a scientific-name value is NOT found via the Node-name field, and vice versa
                    ck( ok, runA( cp, tp, NDF.NodeName, SearchMode.CONTAINS, false, false, "Homo" ).isEmpty(),
                        "Node-name 'Homo' should find nothing (Homo is a scientific name)" );
                    ck( ok, same( runA( cp, tp, NDF.TaxonomyScientificName, SearchMode.CONTAINS, false, false, "Homo" ),
                                  a ),
                        "Scientific-name 'Homo' should find the alpha leaf" );

                    // string modes on the node name
                    ck( ok, same( runA( cp, tp, NDF.NodeName, SearchMode.STARTS_WITH, false, false, "gam" ), c ),
                        "starts-with 'gam' should find gamma_kinase" );
                    ck( ok, same( runA( cp, tp, NDF.NodeName, SearchMode.ENDS_WITH, false, false, "kinase" ), c ),
                        "ends-with 'kinase' should find gamma_kinase" );
                    ck( ok, runA( cp, tp, NDF.NodeName, SearchMode.STARTS_WITH, false, false, "kinase" ).isEmpty(),
                        "starts-with 'kinase' should find nothing" );

                    // Match Case
                    ck( ok, same( runA( cp, tp, NDF.TaxonomyScientificName, SearchMode.CONTAINS, false, false, "homo" ),
                                  a ),
                        "case-insensitive 'homo' should find Homo sapiens" );
                    ck( ok, runA( cp, tp, NDF.TaxonomyScientificName, SearchMode.CONTAINS, true, false, "homo" )
                            .isEmpty(), "case-sensitive 'homo' should find nothing" );

                    // Inverse: complement over data-bearing nodes (the data-less root is excluded)
                    final Set<Long> inv = runA( cp, tp, NDF.NodeName, SearchMode.CONTAINS, false, true, "alpha" );
                    ck( ok, !inv.contains( a.getId() ) && inv.contains( b.getId() ) && inv.contains( c.getId() )
                            && ( inv.size() == 2 ), "inverse of node-name 'alpha' should be {beta, gamma_kinase}" );

                    // a separator-only query with Inverse must RESET the search, not select the whole tree
                    ck( ok, runA( cp, tp, null, SearchMode.CONTAINS, false, true, "," ).isEmpty(),
                        "a separator-only query with Inverse should reset, not select the whole tree" );

                    // A and B are independent (dual-highlight, not combined)
                    cp.setSearchInverseForTest( false );
                    cp.setSearchCaseSensitiveForTest( false );
                    final Set<Long> found_a = runA( cp, tp, null, SearchMode.CONTAINS, false, false, "catus" );
                    final Set<Long> found_b = runB( cp, tp, null, SearchMode.CONTAINS, "musculus" );
                    ck( ok, same( found_a, b ) && same( found_b, c ),
                        "A ('catus') and B ('musculus') should hold independent found sets" );
                    ck( ok, same( nonNull( tp.getFoundNodes0() ), b ),
                        "running B must not disturb A's found set" );

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

    private static Set<Long> runA( final ControlPanel cp, final TreePanel tp, final NDF field, final SearchMode mode,
                                   final boolean case_sensitive, final boolean inverse, final String text ) {
        cp.setSearchFieldForTest( true, field );
        cp.setSearchModeForTest( true, mode );
        cp.setSearchCaseSensitiveForTest( case_sensitive );
        cp.setSearchInverseForTest( inverse );
        cp.getSearchTextField0().setText( text );
        cp.search0();
        return nonNull( tp.getFoundNodes0() );
    }

    private static Set<Long> runB( final ControlPanel cp, final TreePanel tp, final NDF field, final SearchMode mode,
                                   final String text ) {
        cp.setSearchFieldForTest( false, field );
        cp.setSearchModeForTest( false, mode );
        cp.getSearchTextField1().setText( text );
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
        final PhylogenyNode a = leaf( "alpha", "Homo sapiens" );
        final PhylogenyNode b = leaf( "beta", "Felis catus" );
        final PhylogenyNode c = leaf( "gamma_kinase", "Mus musculus" );
        root.addAsChild( a );
        root.addAsChild( b );
        root.addAsChild( c );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.setName( "searchtool" );
        phy.externalNodesHaveChanged();
        out[ 0 ] = a;
        out[ 1 ] = b;
        out[ 2 ] = c;
        return phy;
    }

    private static PhylogenyNode leaf( final String name, final String sci ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        n.getNodeData().setTaxonomy( t );
        return n;
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
