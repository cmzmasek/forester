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

package org.forester.archaeopteryx.tools;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonLineage.Ancestor;

/**
 * Headless unit test for {@link TaxonomySpeciesTreeBuilder}: the induced NCBI-taxonomy species tree merges tip
 * lineages correctly (shared higher taxa become internal nodes), keys by tax-id (homonyms don't merge), and
 * handles the single-species / no-lineage edge cases.
 */
public final class TaxonomySpeciesTreeBuilderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TaxonomySpeciesTreeBuilder: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        boolean ok = true;

        // three species: human + chimp (both Hominidae) and mouse (Rodentia) -- the induced tree must group the two
        // primates together under Hominidae and put mouse on its own, all under Mammalia.
        final Map<PhylogenyNode, TaxonLineage> tips = new LinkedHashMap<PhylogenyNode, TaxonLineage>();
        tips.put( new PhylogenyNode(), lineage( "9606", "species", "Homo sapiens",
                anc( "Mammalia", "class", "40674" ), anc( "Primates", "order", "9443" ),
                anc( "Hominidae", "family", "9604" ), anc( "Homo", "genus", "9605" ) ) );
        tips.put( new PhylogenyNode(), lineage( "9598", "species", "Pan troglodytes",
                anc( "Mammalia", "class", "40674" ), anc( "Primates", "order", "9443" ),
                anc( "Hominidae", "family", "9604" ), anc( "Pan", "genus", "9596" ) ) );
        tips.put( new PhylogenyNode(), lineage( "10090", "species", "Mus musculus",
                anc( "Mammalia", "class", "40674" ), anc( "Rodentia", "order", "9989" ),
                anc( "Muridae", "family", "10066" ), anc( "Mus", "genus", "10088" ) ) );

        final TaxonomySpeciesTreeBuilder.Result r = TaxonomySpeciesTreeBuilder.build( tips );
        final Phylogeny st = r.getSpeciesTree();

        ok &= eq( "species count", r.getSpeciesCount(), 3 );
        ok &= eq( "tips without lineage", r.getTipsWithoutLineage(), 0 );
        // Hominidae groups exactly the two primates (a shared ancestor node, merged from both lineages)
        ok &= sameSet( "Hominidae descendants", leavesUnder( st, "Hominidae" ), "Homo sapiens", "Pan troglodytes" );
        // Mammalia is the common ancestor of all three
        ok &= sameSet( "Mammalia descendants", leavesUnder( st, "Mammalia" ),
                       "Homo sapiens", "Pan troglodytes", "Mus musculus" );
        // Mus musculus is NOT under Primates
        ok &= !leavesUnder( st, "Primates" ).contains( "Mus musculus" )
                && leavesUnder( st, "Primates" ).size() == 2;
        if ( leavesUnder( st, "Primates" ).contains( "Mus musculus" ) ) {
            System.out.println( "  [TaxonomySpeciesTreeBuilderTest] mouse must not be under Primates" );
        }
        // the shared ancestors are single nodes (Mammalia/Hominidae merged, not duplicated per tip)
        ok &= eq( "one Mammalia node", countByTaxon( st, "Mammalia" ), 1 );
        ok &= eq( "one Hominidae node", countByTaxon( st, "Hominidae" ), 1 );
        // the tree is a single rooted tree
        ok &= ( st.getRoot() != null ) && st.isRooted();

        // tax-id keying: two lineages whose family shares a NAME but has DIFFERENT tax-ids must NOT merge into one
        // node (a homonym), so the two species land under separate family nodes -> 2 "Homonymidae" nodes.
        final Map<PhylogenyNode, TaxonLineage> homo = new LinkedHashMap<PhylogenyNode, TaxonLineage>();
        homo.put( new PhylogenyNode(),
                  lineage( "1", "species", "Alpha one", anc( "Homonymidae", "family", "111" ) ) );
        homo.put( new PhylogenyNode(),
                  lineage( "2", "species", "Beta two", anc( "Homonymidae", "family", "222" ) ) );
        final Phylogeny hn = TaxonomySpeciesTreeBuilder.build( homo ).getSpeciesTree();
        ok &= eq( "homonym family NOT merged (distinct tax-ids)", countByTaxon( hn, "Homonymidae" ), 2 );
        // ...but the SAME family tax-id DOES merge (one node grouping both species)
        final Map<PhylogenyNode, TaxonLineage> same = new LinkedHashMap<PhylogenyNode, TaxonLineage>();
        same.put( new PhylogenyNode(),
                  lineage( "1", "species", "Alpha one", anc( "Realidae", "family", "111" ) ) );
        same.put( new PhylogenyNode(),
                  lineage( "2", "species", "Beta two", anc( "Realidae", "family", "111" ) ) );
        ok &= eq( "same family tax-id merged", countByTaxon(
                TaxonomySpeciesTreeBuilder.build( same ).getSpeciesTree(), "Realidae" ), 1 );

        // MIXED PROVENANCE: the same taxon can arrive with a tax-id (an online-resolved lineage) in one tip and
        // NAME-ONLY (a tip's stored in-tree lineage) in another -- it must still MERGE into one node, not split.
        // human's ancestors carry ids; chimp's shared ancestors are name-only (null tax-id).
        final Map<PhylogenyNode, TaxonLineage> mix = new LinkedHashMap<PhylogenyNode, TaxonLineage>();
        mix.put( new PhylogenyNode(), lineage( "9606", "species", "Homo sapiens",
                anc( "Mammalia", "class", "40674" ), anc( "Primates", "order", "9443" ),
                anc( "Hominidae", "family", "9604" ) ) );
        mix.put( new PhylogenyNode(), lineage( null, "species", "Pan troglodytes",
                anc( "Mammalia", "class", null ), anc( "Primates", "order", null ),
                anc( "Hominidae", "family", null ) ) );
        final Phylogeny mx = TaxonomySpeciesTreeBuilder.build( mix ).getSpeciesTree();
        ok &= eq( "mixed-provenance Mammalia merged", countByTaxon( mx, "Mammalia" ), 1 );
        ok &= eq( "mixed-provenance Hominidae merged", countByTaxon( mx, "Hominidae" ), 1 );
        ok &= sameSet( "mixed Hominidae groups both", leavesUnder( mx, "Hominidae" ),
                       "Homo sapiens", "Pan troglodytes" );
        // duplicate-leaf guard: the SAME species with an id in one tip and name-only in another -> ONE leaf
        // (else GSDIR rejects a non-unique species-tree taxonomy and the reconciliation hard-fails)
        final Map<PhylogenyNode, TaxonLineage> dup = new LinkedHashMap<PhylogenyNode, TaxonLineage>();
        dup.put( new PhylogenyNode(), lineage( "9606", "species", "Homo sapiens", anc( "Homo", "genus", "9605" ) ) );
        dup.put( new PhylogenyNode(), lineage( null, "species", "Homo sapiens", anc( "Homo", "genus", null ) ) );
        ok &= eq( "same species (id + name-only) -> one leaf",
                  TaxonomySpeciesTreeBuilder.build( dup ).getSpeciesTree().getNumberOfExternalNodes(), 1 );

        // single species -> a one-leaf tree; empty lineage -> counted as without-lineage, not placed
        final Map<PhylogenyNode, TaxonLineage> mixed = new LinkedHashMap<PhylogenyNode, TaxonLineage>();
        mixed.put( new PhylogenyNode(), lineage( "9606", "species", "Homo sapiens", anc( "Homo", "genus", "9605" ) ) );
        mixed.put( new PhylogenyNode(), TaxonLineage.EMPTY );
        final TaxonomySpeciesTreeBuilder.Result rm = TaxonomySpeciesTreeBuilder.build( mixed );
        ok &= eq( "single-species leaf count", rm.getSpeciesCount(), 1 );
        ok &= eq( "empty lineage counted", rm.getTipsWithoutLineage(), 1 );

        return ok;
    }

    private static TaxonLineage lineage( final String tax_id, final String rank, final String sci,
                                         final Ancestor... ancestors ) {
        final List<Ancestor> a = new ArrayList<Ancestor>();
        for( final Ancestor an : ancestors ) {
            a.add( an );
        }
        return new TaxonLineage( tax_id, rank, sci, null, a );
    }

    private static Ancestor anc( final String name, final String rank, final String tax_id ) {
        return new Ancestor( name, rank, tax_id );
    }

    /** Scientific names of the external descendants of the (single) node whose taxonomy scientific name is {@code taxon}. */
    private static List<String> leavesUnder( final Phylogeny phy, final String taxon ) {
        final List<String> names = new ArrayList<String>();
        final PhylogenyNode n = firstByTaxon( phy, taxon );
        if ( n != null ) {
            for( final PhylogenyNode ext : n.getAllExternalDescendants() ) {
                names.add( scientificName( ext ) );
            }
        }
        return names;
    }

    private static PhylogenyNode firstByTaxon( final Phylogeny phy, final String taxon ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( taxon.equals( scientificName( n ) ) ) {
                return n;
            }
        }
        return null;
    }

    private static int countByTaxon( final Phylogeny phy, final String taxon ) {
        int c = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            if ( taxon.equals( scientificName( it.next() ) ) ) {
                ++c;
            }
        }
        return c;
    }

    private static String scientificName( final PhylogenyNode n ) {
        if ( n.getNodeData().isHasTaxonomy() && ( n.getNodeData().getTaxonomy().getScientificName() != null ) ) {
            return n.getNodeData().getTaxonomy().getScientificName();
        }
        return null;
    }

    private static boolean eq( final String name, final int a, final int b ) {
        if ( a != b ) {
            System.out.println( "  [TaxonomySpeciesTreeBuilderTest] " + name + ": got " + a + ", expected " + b );
            return false;
        }
        return true;
    }

    private static boolean sameSet( final String name, final List<String> got, final String... expected ) {
        boolean ok = ( got.size() == expected.length );
        for( final String e : expected ) {
            if ( !got.contains( e ) ) {
                ok = false;
            }
        }
        if ( !ok ) {
            System.out.println( "  [TaxonomySpeciesTreeBuilderTest] " + name + ": got " + got );
        }
        return ok;
    }

    private TaxonomySpeciesTreeBuilderTest() {
    }
}
