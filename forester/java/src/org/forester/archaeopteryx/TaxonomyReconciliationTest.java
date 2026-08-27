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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.forester.archaeopteryx.tools.TaxonomySpeciesTreeBuilder;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sdi.GSDIR;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonLineage.Ancestor;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * Headless end-to-end test of the "reconcile using NCBI taxonomy" flow (the pieces {@code MainFrame.
 * executeGSDIRwithTaxonomySpeciesTree} wires together): resolve each gene-tip's lineage via the service
 * ({@link TreePanelUtil#tipLineages}), build the induced taxonomy species tree
 * ({@link TaxonomySpeciesTreeBuilder}), and run {@link GSDIR}. A paralog duplication that predates the
 * human/chimp split must be recovered as a duplication. Uses a FAKE lineage service, so it needs no network.
 */
public final class TaxonomyReconciliationTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TaxonomyReconciliation: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            // gene family duplicated BEFORE the human/chimp split: two paralog clades (A and B), each with a human
            // and a chimp gene, plus a mouse ortholog. The node joining the two clades spans the SAME species on
            // both sides -> a duplication.
            final Phylogeny gene_tree = Phylogeny
                    .createInstanceFromNhxString( "(((humanA,chimpA),(humanB,chimpB)),mouse);" );
            final Map<String, String> tip_species = new HashMap<String, String>();
            tip_species.put( "humanA", "Homo sapiens" );
            tip_species.put( "humanB", "Homo sapiens" );
            tip_species.put( "chimpA", "Pan troglodytes" );
            tip_species.put( "chimpB", "Pan troglodytes" );
            tip_species.put( "mouse", "Mus musculus" );
            for( final PhylogenyNodeIterator it = gene_tree.iteratorExternalForward(); it.hasNext(); ) {
                final PhylogenyNode tip = it.next();
                final Taxonomy t = new Taxonomy();
                t.setScientificName( tip_species.get( tip.getName() ) );
                tip.getNodeData().setTaxonomy( t );
            }

            // resolve lineages via the (fake) service exactly as the GUI does, then build the species tree
            final TaxonomicLineageService service = new FakeService();
            final Map<PhylogenyNode, TaxonLineage> tip_lineages = TreePanelUtil.tipLineages( gene_tree, service );
            final TaxonomySpeciesTreeBuilder.Result r = TaxonomySpeciesTreeBuilder.build( tip_lineages );
            final Phylogeny species_tree = r.getSpeciesTree();
            if ( r.getSpeciesCount() != 3 ) {
                return fail( "species tree should have 3 species (human/chimp/mouse), got " + r.getSpeciesCount() );
            }

            // run the SAME reconciler the GUI uses
            final GSDIR gsdir = new GSDIR( gene_tree.copy(), species_tree.copy(), true, true, true );
            if ( gsdir.getMinDuplicationsSum() < 1 ) {
                return fail( "the pre-split paralog duplication must be recovered; duplications = "
                        + gsdir.getMinDuplicationsSum() );
            }
            if ( gsdir.getStrippedExternalGeneTreeNodes().size() != 0 ) {
                return fail( "all gene tips should map to a species; stripped = "
                        + gsdir.getStrippedExternalGeneTreeNodes().size() );
            }
            // the reconciled gene tree must carry at least one duplication event
            final Phylogeny result = gsdir.getMinDuplicationsSumGeneTree();
            int dups = 0;
            for( final PhylogenyNodeIterator it = result.iteratorPostorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( n.getNodeData().isHasEvent() && n.getNodeData().getEvent().isDuplication() ) {
                    ++dups;
                }
            }
            if ( dups < 1 ) {
                return fail( "the reconciled gene tree must carry a duplication event, got " + dups );
            }

            // an all-orthologs gene tree (one gene per species, matching the species topology) has NO duplication
            final Phylogeny orthologs = Phylogeny.createInstanceFromNhxString( "((human,chimp),mouse);" );
            final Map<String, String> osp = new HashMap<String, String>();
            osp.put( "human", "Homo sapiens" );
            osp.put( "chimp", "Pan troglodytes" );
            osp.put( "mouse", "Mus musculus" );
            for( final PhylogenyNodeIterator it = orthologs.iteratorExternalForward(); it.hasNext(); ) {
                final PhylogenyNode tip = it.next();
                final Taxonomy t = new Taxonomy();
                t.setScientificName( osp.get( tip.getName() ) );
                tip.getNodeData().setTaxonomy( t );
            }
            final Phylogeny osp_tree = TaxonomySpeciesTreeBuilder
                    .build( TreePanelUtil.tipLineages( orthologs, service ) ).getSpeciesTree();
            final GSDIR og = new GSDIR( orthologs.copy(), osp_tree.copy(), true, true, true );
            if ( og.getMinDuplicationsSum() != 0 ) {
                return fail( "an all-orthologs tree must have 0 duplications, got " + og.getMinDuplicationsSum() );
            }
            return true;
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    /** Fake lineage service: returns a fixed NCBI-like lineage for the three demo species (no network). */
    private static final class FakeService implements TaxonomicLineageService {

        @Override
        public TaxonLineage lineageOf( final String taxon ) {
            if ( "Homo sapiens".equals( taxon ) ) {
                return lin( "9606", "Homo sapiens", anc( "Mammalia", "class", "40674" ),
                            anc( "Primates", "order", "9443" ), anc( "Hominidae", "family", "9604" ),
                            anc( "Homo", "genus", "9605" ) );
            }
            if ( "Pan troglodytes".equals( taxon ) ) {
                return lin( "9598", "Pan troglodytes", anc( "Mammalia", "class", "40674" ),
                            anc( "Primates", "order", "9443" ), anc( "Hominidae", "family", "9604" ),
                            anc( "Pan", "genus", "9596" ) );
            }
            if ( "Mus musculus".equals( taxon ) ) {
                return lin( "10090", "Mus musculus", anc( "Mammalia", "class", "40674" ),
                            anc( "Rodentia", "order", "9989" ), anc( "Muridae", "family", "10066" ),
                            anc( "Mus", "genus", "10088" ) );
            }
            return null;
        }

        @Override
        public TaxonLineage fetch( final String taxon ) throws IOException {
            return lineageOf( taxon );
        }
    }

    private static TaxonLineage lin( final String tax_id, final String sci, final Ancestor... ancestors ) {
        final List<Ancestor> a = new ArrayList<Ancestor>();
        for( final Ancestor an : ancestors ) {
            a.add( an );
        }
        return new TaxonLineage( tax_id, "species", sci, null, a );
    }

    private static Ancestor anc( final String name, final String rank, final String tax_id ) {
        return new Ancestor( name, rank, tax_id );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [TaxonomyReconciliationTest] " + msg );
        return false;
    }

    private TaxonomyReconciliationTest() {
    }
}
