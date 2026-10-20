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

import java.util.Arrays;
import java.util.List;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * Tests for {@link RankSourcePanel}: the chooser entries state what the chosen source places, a source change
 * rebuilds them keeping the selected rank, and the first / "(none)" entries survive. Headless (the Swing
 * components are built, never shown).
 */
public final class RankSourcePanelTest {

    public static void main( final String[] args ) {
        System.out.println( "RankSourcePanel: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        return testChoices() && testPanel() && testEntryAndBareRank();
    }

    /** order Carnivora over family Felidae over (cat, lion); order Rodentia over (rat, mouse): four species tips --
     *  four, so that order's two taxa repeat enough to be offered (the joint rule offers a rank when at most three
     *  taxa fall on five placed tips; two orders on three tips would not). */
    private static Phylogeny fixture() {
        final PhylogenyNode carnivora = annotated( "order", "Carnivora" );
        final PhylogenyNode felidae = annotated( "family", "Felidae" );
        felidae.addAsChild( tip( "cat", "Felis catus" ) );
        felidae.addAsChild( tip( "lion", "Panthera leo" ) );
        carnivora.addAsChild( felidae );
        final PhylogenyNode rodentia = annotated( "order", "Rodentia" );
        rodentia.addAsChild( tip( "rat", "Rattus norvegicus" ) );
        rodentia.addAsChild( tip( "mouse", "Mus musculus" ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( carnivora );
        root.addAsChild( rodentia );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean testChoices() {
        final Phylogeny phy = fixture();
        // the cache knows the rat's family (and nothing else), so the family share depends on the source
        final TaxonomicLineageService svc = new RatService();
        final List<String> tree = Arrays.asList( RankSourcePanel.choices( phy, svc, TaxonSource.TREE ) );
        // the joint offer rule (0.11.186): a rank is offered when it places a tip and its taxa repeat (at most three
        // taxa per five placed tips); species names four taxa on four tips and genus places nothing: neither is
        // offered, and the tree source never offers a rank nothing carries
        if ( !tree.equals( Arrays.asList( "order (2 annotated; 100% placed)", "family (1 annotated; 50% placed)" ) ) ) {
            return fail( "TREE entries: the offered ranks in canonical order, each with its annotation count and what the tree places; got " + tree );
        }
        for( final String e : tree ) {
            if ( e.contains( "look up" ) ) {
                return fail( "the tree source never looks anything up: " + e );
            }
        }
        final List<String> db = Arrays.asList( RankSourcePanel.choices( phy, svc, TaxonSource.DATABASE ) );
        // placed NOW: cat and lion by their ancestors (the cache cannot place them), the rat by the cache, the mouse
        // by neither at family (Rodentia names no family) and by its ancestor at order; every uncached tip is to be
        // looked up; an online source also offers the major ranks nothing carries (genus) for what it would look up
        if ( !db.contains( "family (1 annotated; 75% placed now; 3 to look up)" )
                || !db.contains( "order (2 annotated; 100% placed now; 3 to look up)" )
                || !db.contains( "genus (not annotated; 0% placed now; 3 to look up)" ) ) {
            return fail( "DATABASE entries count what is placed NOW and every tip it would look up; got " + db );
        }
        final List<String> gaps = Arrays.asList( RankSourcePanel.choices( phy, svc, TaxonSource.TREE_THEN_DATABASE ) );
        if ( !gaps.contains( "family (1 annotated; 75% placed now; 1 to look up)" )
                || !gaps.contains( "order (2 annotated; 100% placed now; 0 to look up)" ) ) {
            return fail( "TREE_THEN_DATABASE looks up only the tips the tree leaves unplaced; got " + gaps );
        }
        if ( RankSourcePanel.choices( new Phylogeny(), svc, TaxonSource.TREE ).length != 0
                || RankSourcePanel.choices( null, svc, TaxonSource.TREE ).length != 0 ) {
            return fail( "an empty or null tree has no entries" );
        }
        return true;
    }

    private static boolean testPanel() {
        final Phylogeny phy = fixture();
        final RankSourcePanel p = new RankSourcePanel( phy, new EmptyService(), TaxonSource.TREE, "FIRST", 2 );
        if ( !"FIRST".equals( p.rankBox().getItemAt( 0 ) ) || ( p.selectedSource() != TaxonSource.TREE )
                || !p.radioForTest( TaxonSource.TREE ).isSelected() ) {
            return fail( "the first entry leads the main chooser and the initial source is selected" );
        }
        if ( ( p.extraRankBoxes().size() != 2 ) || !RankSourcePanel.NONE.equals( p.extraRankBoxes().get( 0 ).getSelectedItem() )
                || ( p.selectedExtraRank( 0 ) != null ) ) {
            return fail( "extra choosers open on (none)" );
        }
        p.selectRank( "family" );
        p.extraRankBoxes().get( 1 ).setSelectedItem( "order (2 annotated; 100% placed)" );
        if ( !"family".equals( p.selectedRank() ) || !"order".equals( p.selectedExtraRank( 1 ) ) ) {
            return fail( "selecting by bare rank" );
        }
        final int[] rebuilt = { 0 };
        p.addSourceListener( () -> rebuilt[ 0 ]++ );
        p.setSelectedSource( TaxonSource.DATABASE );
        if ( ( rebuilt[ 0 ] != 1 ) || ( p.selectedSource() != TaxonSource.DATABASE ) ) {
            return fail( "a source change rebuilds once and is reported" );
        }
        if ( !"family".equals( p.selectedRank() ) || !"order".equals( p.selectedExtraRank( 1 ) )
                || ( p.selectedExtraRank( 0 ) != null ) || !"FIRST".equals( p.rankBox().getItemAt( 0 ) ) ) {
            return fail( "a rebuild keeps every chooser's selection and the first entry; got " + p.rankBox().getSelectedItem()
                    + " / " + p.extraRankBoxes().get( 1 ).getSelectedItem() );
        }
        if ( !String.valueOf( p.rankBox().getSelectedItem() ).contains( "to look up" ) ) {
            return fail( "after switching to an online source the entries say what would be looked up; got "
                    + p.rankBox().getSelectedItem() );
        }
        p.setSelectedIndex( 0 );
        if ( !p.isFirstEntrySelected() || ( p.selectedRank() != null ) ) {
            return fail( "the first entry is not a rank" );
        }
        p.setSelectedSource( TaxonSource.TREE );
        if ( !p.isFirstEntrySelected() ) {
            return fail( "a rebuild keeps the first entry selected" );
        }
        p.selectRank( "no such rank" );
        if ( !p.isFirstEntrySelected() ) {
            return fail( "selecting a rank the chooser does not offer changes nothing" );
        }
        // the write-lineages checkbox: on by default, live only for an online source, and never "on" for the tree
        final RankSourcePanel w = new RankSourcePanel( phy, new EmptyService(), TaxonSource.TREE, null, 0 );
        if ( !w.writeLineagesBox().isSelected() || w.writeLineagesBox().isEnabled() || w.isWriteLineages() ) {
            return fail( "for the tree source the checkbox is checked but disabled, and no lineages are written" );
        }
        w.setSelectedSource( TaxonSource.DATABASE );
        if ( !w.writeLineagesBox().isEnabled() || !w.isWriteLineages() ) {
            return fail( "an online source enables the checkbox, and it is on by default" );
        }
        w.writeLineagesBox().setSelected( false );
        if ( w.isWriteLineages() ) {
            return fail( "unchecked: no lineages are written" );
        }
        w.setSelectedSource( TaxonSource.TREE_THEN_DATABASE );
        if ( !w.writeLineagesBox().isEnabled() || w.isWriteLineages() ) {
            return fail( "the other online source enables it too; the user's unchecking is kept" );
        }
        w.setSelectedSource( TaxonSource.TREE );
        if ( w.writeLineagesBox().isEnabled() ) {
            return fail( "back to the tree source: disabled again" );
        }
        // no first entry: the ranks start at index 0 and nothing is a first entry
        final RankSourcePanel q = new RankSourcePanel( phy, new EmptyService(), TaxonSource.DATABASE, null, 0 );
        if ( q.isFirstEntrySelected() || !q.rankBox().getItemAt( 0 ).startsWith( "domain" ) || !q.extraRankBoxes().isEmpty() ) {
            return fail( "without a first entry the chooser starts with the ranks; got " + q.rankBox().getItemAt( 0 ) );
        }
        return true;
    }

    private static boolean testEntryAndBareRank() {
        if ( !"family (19 annotated; 29% placed)".equals( RankSourcePanel.entry( "family", 19, 55, 191, 136, false ) )
                || !"family (19 annotated; 29% placed now; 136 to look up)".equals( RankSourcePanel.entry( "family", 19, 55, 191, 136, true ) )
                || !"order (not annotated; 0% placed now; 191 to look up)".equals( RankSourcePanel.entry( "order", 0, 0, 191, 191, true ) ) ) {
            return fail( "entry format" );
        }
        if ( !"family".equals( RankSourcePanel.bareRank( "family (19 annotated; 29% placed)" ) )
                || !"family".equals( RankSourcePanel.bareRank( "family" ) ) || ( RankSourcePanel.bareRank( RankSourcePanel.NONE ) != null )
                || ( RankSourcePanel.bareRank( "" ) != null ) || ( RankSourcePanel.bareRank( null ) != null ) ) {
            return fail( "bareRank" );
        }
        if ( ( TaxonSource.defaultFor( true ) != TaxonSource.TREE ) || ( TaxonSource.defaultFor( false ) != TaxonSource.DATABASE )
                || TaxonSource.TREE.isOnline() || !TaxonSource.DATABASE.isOnline() || !TaxonSource.TREE_THEN_DATABASE.isOnline() ) {
            return fail( "a tree with ranks opens offline, a bare tree online; only the database sources are online" );
        }
        final String[] with = RankSourcePanel.withFirst( new String[] { "a", "b" }, "x" );
        if ( ( with.length != 3 ) || !"x".equals( with[ 0 ] ) || !"b".equals( with[ 2 ] ) ) {
            return fail( "withFirst" );
        }
        return true;
    }

    private static PhylogenyNode annotated( final String rank, final String sci ) {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        try {
            t.setRank( rank );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode tip( final String name, final String sci ) {
        final PhylogenyNode n = annotated( "species", sci );
        n.setName( name );
        return n;
    }

    /** A service that has nothing cached and must not be fetched from. */
    private static class EmptyService implements TaxonomicLineageService {

        @Override
        public TaxonLineage lineageOf( final String taxon ) {
            return null;
        }

        @Override
        public TaxonLineage fetch( final String taxon ) {
            throw new UnsupportedOperationException( "the chooser must never fetch" );
        }
    }

    /** Knows the rat's family only. */
    private static final class RatService extends EmptyService {

        @Override
        public TaxonLineage lineageOf( final String taxon ) {
            return "Rattus norvegicus".equalsIgnoreCase( taxon ) ? new TaxonLineage( null, null, null, null,
                    java.util.Collections.singletonList( new TaxonLineage.Ancestor( "Muridae", "family", null ) ) ) : null;
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  RankSourcePanel test failed: " + msg );
        return false;
    }

    private RankSourcePanelTest() {
    }
}
