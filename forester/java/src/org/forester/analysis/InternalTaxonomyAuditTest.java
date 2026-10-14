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

package org.forester.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.forester.analysis.InternalTaxonomyAudit.Finding;
import org.forester.analysis.InternalTaxonomyAudit.Verdict;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonLineage.Ancestor;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * Headless tests for {@link InternalTaxonomyAudit}: every verdict is reached from a fixture built like the real
 * case (a bacterial backbone whose "Bacteria" nodes were resolved to the stick-insect genus), the proposals are the
 * right ones, {@code apply} writes exactly them, and the lookup list and summary are right. Tip lineages are given
 * directly; the database is an in-memory cache.
 */
public final class InternalTaxonomyAuditTest {

    public static boolean test() {
        return testVerdicts() && testApply() && testLookupsSummaryAndEdges() && testLevelOf();
    }

    /** The fixture: see each node's comment. {@code map} receives the tip lineages. */
    private static final class Fixture {

        final Map<PhylogenyNode, TaxonLineage> map = new IdentityHashMap<PhylogenyNode, TaxonLineage>();
        final FakeService                       svc = new FakeService();
        PhylogenyNode                           stick_insect, too_deep, wrong_rank, wrong_name, fine, unknown, unchecked, no_proposal;
        Phylogeny                               phy;

        Fixture() {
            // tips: two E. coli and a Bacillus subtilis (all Bacteria [2]), two cats (Carnivora > Felidae) and a dog
            final PhylogenyNode ecoli1 = tip( "ecoli1", map, bacterium( "Escherichia", "561", "Enterobacterales", "91347" ) );
            final PhylogenyNode ecoli2 = tip( "ecoli2", map, bacterium( "Escherichia", "561", "Enterobacterales", "91347" ) );
            final PhylogenyNode bsub = tip( "bsub", map, bacterium( "Bacillus", "1386", "Bacillales", "1385" ) );
            final PhylogenyNode cat1 = tip( "cat1", map, mammal( "Felidae", "9681" ) );
            final PhylogenyNode cat2 = tip( "cat2", map, mammal( "Felidae", "9681" ) );
            final PhylogenyNode dog = tip( "dog", map, mammal( "Canidae", "9608" ) );
            final PhylogenyNode orphan = tip( "orphan", map, null ); // no lineage at all
            // the stick insect: "Bacteria Latreille et al. 1825" [629395], genus -- over bacteria. WRONG_TAXON,
            // proposed: the deepest taxon all three share = Bacteria [2], superkingdom
            stick_insect = internal( "Bacteria Latreille et al. 1825", "genus", "629395" );
            final PhylogenyNode enterobacterales = internal( "Enterobacterales", "order", "91347" ); // fine (ids)
            enterobacterales.addAsChild( ecoli1 );
            enterobacterales.addAsChild( ecoli2 );
            stick_insect.addAsChild( enterobacterales );
            stick_insect.addAsChild( bsub );
            fine = enterobacterales;
            // Felidae over cats AND the dog: NOT_SHARED_BY_ALL_TIPS, proposed Carnivora [33554]
            too_deep = internal( "Felidae", "family", "9681" );
            // genus Carnivora over the cats: the taxon is right, the rank is not: WRONG_RANK -> order (from the tips)
            wrong_rank = internal( "Carnivora", "genus", null );
            wrong_rank.addAsChild( cat1 );
            wrong_rank.addAsChild( cat2 );
            too_deep.addAsChild( wrong_rank );
            too_deep.addAsChild( dog );
            // id 9608 named wrongly, over the dog alone: the database (cached) names it Canidae: WRONG_NAME
            wrong_name = internal( "Canidae-typo", "family", "9608" );
            svc.know( "9608", new TaxonLineage( "9608", "family", "Canidae", null, Collections.singletonList( new Ancestor( "Carnivora", "order", "33554" ) ) ) );
            // an unknown name over the orphan only: nothing cached but asked -> UNKNOWN_TO_DATABASE (EMPTY cached)
            unknown = internal( "Nonexistus", null, null );
            svc.empty( "Nonexistus" );
            unknown.addAsChild( orphan );
            // a name nobody looked up over the orphan: NOT_CHECKED
            unchecked = internal( "Maybeus", null, null );
            unchecked.addAsChild( tip( "orphan2", map, null ) );
            // a wrong taxon whose tips share nothing: a finding without a proposal
            no_proposal = internal( "Felidae", "family", "9681" );
            no_proposal.addAsChild( tip( "x", map, bacterium( "Escherichia", "561", "Enterobacterales", "91347" ) ) );
            no_proposal.addAsChild( tip( "y", map, new TaxonLineage( "7", "species", "Alienus", null, new ArrayList<Ancestor>() ) ) );
            final PhylogenyNode root = new PhylogenyNode();
            final PhylogenyNode wrong_name_clade = new PhylogenyNode();
            wrong_name_clade.addAsChild( wrong_name );
            wrong_name.addAsChild( tip( "dog2", map, mammal( "Canidae", "9608" ) ) );
            for( final PhylogenyNode n : new PhylogenyNode[] { stick_insect, too_deep, wrong_name_clade, unknown, unchecked, no_proposal } ) {
                root.addAsChild( n );
            }
            phy = new Phylogeny();
            phy.setRoot( root );
            phy.externalNodesHaveChanged();
            phy.setName( "audit" );
        }
    }

    private static boolean testVerdicts() {
        final Fixture f = new Fixture();
        final List<Finding> findings = InternalTaxonomyAudit.audit( f.phy, f.map, f.svc );
        final Map<PhylogenyNode, Finding> by_node = new IdentityHashMap<PhylogenyNode, Finding>();
        for( final Finding x : findings ) {
            by_node.put( x.getNode(), x );
        }
        if ( findings.size() != 8 ) {
            return fail( "eight annotated internal nodes, got " + findings.size() );
        }
        final Finding stick = by_node.get( f.stick_insect );
        if ( ( stick.getVerdict() != Verdict.WRONG_TAXON ) || !"superkingdom Bacteria [2]".equals( stick.getProposed() )
                || !stick.getDetail().contains( "none of the 3 resolved tips" ) ) {
            return fail( "the stick insect is a wrong taxon, replaced by the Bacteria all tips share; got " + stick.getVerdict() + " / " + stick.getProposed() + " / " + stick.getDetail() );
        }
        if ( !"genus Bacteria Latreille et al. 1825 [629395]".equals( stick.getCurrent() ) ) {
            return fail( "the current annotation reads rank name [id]; got " + stick.getCurrent() );
        }
        final Finding deep = by_node.get( f.too_deep );
        if ( ( deep.getVerdict() != Verdict.NOT_SHARED_BY_ALL_TIPS ) || !"order Carnivora [33554]".equals( deep.getProposed() )
                || !deep.getDetail().contains( "only 2 of the 3 resolved tips" ) ) {
            return fail( "Felidae over cats and a dog is too deep, Carnivora proposed; got " + deep.getVerdict() + " / " + deep.getProposed() );
        }
        final Finding rank = by_node.get( f.wrong_rank );
        if ( ( rank.getVerdict() != Verdict.WRONG_RANK ) || !"order Carnivora [33554]".equals( rank.getProposed() )
                || !rank.getDetail().contains( "tips' lineages give the rank order" ) ) {
            return fail( "genus Carnivora has the wrong rank; the tips say order (and give the id); got " + rank.getVerdict() + " / " + rank.getProposed() + " / " + rank.getDetail() );
        }
        final Finding name = by_node.get( f.wrong_name );
        if ( ( name.getVerdict() != Verdict.WRONG_NAME ) || !"family Canidae [9608]".equals( name.getProposed() )
                || !name.getDetail().contains( "NCBI names id 9608 \"Canidae\"" ) ) {
            return fail( "id 9608 is Canidae; got " + name.getVerdict() + " / " + name.getProposed() + " / " + name.getDetail() );
        }
        final Finding fine = by_node.get( f.fine );
        if ( ( fine.getVerdict() != Verdict.OK ) || fine.hasProposal() || !fine.getDetail().contains( "all 2 resolved tips" ) ) {
            return fail( "Enterobacterales over two E. coli is fine; got " + fine.getVerdict() + " / " + fine.getDetail() );
        }
        if ( by_node.get( f.unknown ).getVerdict() != Verdict.UNKNOWN_TO_DATABASE ) {
            return fail( "a cached-empty name over an unresolved tip is unknown to the database; got " + by_node.get( f.unknown ).getVerdict() );
        }
        if ( by_node.get( f.unchecked ).getVerdict() != Verdict.NOT_CHECKED ) {
            return fail( "nothing cached, no tip lineage: not checked; got " + by_node.get( f.unchecked ).getVerdict() );
        }
        final Finding none = by_node.get( f.no_proposal );
        if ( ( none.getVerdict() != Verdict.WRONG_TAXON ) || none.hasProposal() || !"".equals( none.getProposed() )
                || !none.getDetail().contains( "share no taxon" ) ) {
            return fail( "a wrong taxon over tips that share nothing has no proposal; got " + none.getVerdict() + " / " + none.getProposed() );
        }
        // a correct annotation whose rank NCBI calls "domain" and the tree "superkingdom" is not a rank mismatch
        if ( !InternalTaxonomyAudit.sameRank( "domain", "superkingdom" ) || !InternalTaxonomyAudit.sameRank( "Order", "order" )
                || InternalTaxonomyAudit.sameRank( "order", "family" ) ) {
            return fail( "sameRank aliases domain and superkingdom, nothing else" );
        }
        // preorder, and a node without name and id is skipped
        if ( findings.get( 0 ).getNode() != f.stick_insect ) {
            return fail( "findings come in preorder" );
        }
        return true;
    }

    private static boolean testApply() {
        final Fixture f = new Fixture();
        final List<Finding> findings = InternalTaxonomyAudit.audit( f.phy, f.map, f.svc );
        // keep a field the in-place corrections must not lose
        f.wrong_rank.getNodeData().getTaxonomy().setCommonName( "carnivores" );
        final int changed = InternalTaxonomyAudit.apply( findings, EnumSet.allOf( Verdict.class ) );
        if ( changed != 4 ) {
            return fail( "four findings carry a proposal; applied " + changed );
        }
        final Taxonomy stick = f.stick_insect.getNodeData().getTaxonomy();
        if ( !"Bacteria".equals( stick.getScientificName() ) || !"superkingdom".equals( stick.getRank() )
                || ( stick.getIdentifier() == null ) || !"2".equals( stick.getIdentifier().getValue() ) ) {
            return fail( "the wrong taxon is REPLACED by the shared one; got " + stick );
        }
        final Taxonomy deep = f.too_deep.getNodeData().getTaxonomy();
        if ( !"Carnivora".equals( deep.getScientificName() ) || !"order".equals( deep.getRank() ) ) {
            return fail( "the too-deep taxon is replaced too; got " + deep );
        }
        final Taxonomy rank = f.wrong_rank.getNodeData().getTaxonomy();
        if ( !"order".equals( rank.getRank() ) || !"Carnivora".equals( rank.getScientificName() )
                || !"carnivores".equals( rank.getCommonName() ) || ( rank.getIdentifier() == null )
                || !"33554".equals( rank.getIdentifier().getValue() ) ) {
            return fail( "a wrong rank is corrected IN the taxonomy, other fields kept, the id filled in; got " + rank + " / " + rank.getCommonName() );
        }
        final Taxonomy name = f.wrong_name.getNodeData().getTaxonomy();
        if ( !"Canidae".equals( name.getScientificName() ) || !"family".equals( name.getRank() ) ) {
            return fail( "a wrong name is corrected in place; got " + name );
        }
        if ( !"Felidae".equals( f.no_proposal.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "a finding without a proposal changes nothing" );
        }
        // a second audit of the corrected tree finds nothing left to correct
        final List<Finding> again = InternalTaxonomyAudit.audit( f.phy, f.map, f.svc );
        for( final Finding x : again ) {
            if ( x.hasProposal() && ( x.getNode() != f.no_proposal ) ) {
                return fail( "after applying, nothing but the un-proposable finding remains: " + x.getNode() + " " + x.getVerdict() + " " + x.getProposed() );
            }
        }
        // apply restricted to one verdict touches only that verdict
        final Fixture g = new Fixture();
        final int only_ranks = InternalTaxonomyAudit.apply( InternalTaxonomyAudit.audit( g.phy, g.map, g.svc ),
                                                            EnumSet.of( Verdict.WRONG_RANK ) );
        if ( ( only_ranks != 1 ) || !"Bacteria Latreille et al. 1825".equals( g.stick_insect.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "apply honours the chosen verdicts; applied " + only_ranks );
        }
        if ( ( InternalTaxonomyAudit.apply( null, null ) != 0 ) || ( InternalTaxonomyAudit.apply( new ArrayList<Finding>(), null ) != 0 ) ) {
            return fail( "nothing to apply" );
        }
        return true;
    }

    private static boolean testLookupsSummaryAndEdges() {
        final Fixture f = new Fixture();
        final java.util.Set<String> keys = InternalTaxonomyAudit.internalTaxaToLookUp( f.phy, f.svc );
        // by id where there is one, else by name; the cached ones (9608, Nonexistus) excluded
        if ( !keys.equals( new java.util.TreeSet<String>( java.util.Arrays.asList( "629395", "91347", "9681", "Carnivora", "Maybeus" ) ) ) ) {
            return fail( "the internal taxa to look up: " + keys );
        }
        if ( !InternalTaxonomyAudit.internalTaxaToLookUp( new Phylogeny(), f.svc ).isEmpty() ) {
            return fail( "an empty tree has nothing to look up" );
        }
        final List<Finding> findings = InternalTaxonomyAudit.audit( f.phy, f.map, f.svc );
        final String summary = InternalTaxonomyAudit.summary( findings );
        for( final String piece : new String[] { "Checked 8 internal annotations", "1 consistent", "2 wrong taxon", "1 not shared by all tips",
                "1 wrong rank", "1 name differs from NCBI's", "1 unknown to the database", "1 not checked", "4 corrections proposed." } ) {
            if ( !summary.contains( piece ) ) {
                return fail( "summary lacks \"" + piece + "\": " + summary );
            }
        }
        final Map<Verdict, Integer> counts = InternalTaxonomyAudit.countByVerdict( findings );
        if ( ( counts.get( Verdict.WRONG_TAXON ) != 2 ) || ( counts.get( Verdict.OK ) != 1 ) || ( counts.size() != Verdict.values().length ) ) {
            return fail( "countByVerdict: " + counts );
        }
        final String prov = InternalTaxonomyAudit.provenance( f.phy, 4 );
        if ( !prov.startsWith( "Used the internal-taxonomy check (tip lineages and NCBI Taxonomy) to correct 4 internal node taxonomies" )
                || !prov.endsWith( "in tree named \"audit\" with 11 tips." ) ) {
            return fail( "provenance: " + prov );
        }
        if ( !InternalTaxonomyAudit.audit( new Phylogeny(), f.map, f.svc ).isEmpty()
                || !InternalTaxonomyAudit.audit( null, f.map, f.svc ).isEmpty() ) {
            return fail( "an empty or null tree has no findings" );
        }
        // no service at all: the tips still judge; the database cannot
        final List<Finding> offline = InternalTaxonomyAudit.audit( f.phy, f.map, null );
        if ( offline.size() != 8 ) {
            return fail( "offline: every annotated node is still examined" );
        }
        for( final Finding x : offline ) {
            if ( x.getNode() == f.wrong_name && ( x.getVerdict() != Verdict.OK ) ) {
                return fail( "offline, the misnamed id over its own dog reads consistent (the tips carry id 9608); got " + x.getVerdict() );
            }
            if ( x.getNode() == f.unknown && ( x.getVerdict() != Verdict.NOT_CHECKED ) ) {
                return fail( "offline, an unknown name is merely not checked" );
            }
        }
        if ( !"order Carnivora [33554]".equals( InternalTaxonomyAudit.describe( "Carnivora", "order", "33554" ) )
                || !"Carnivora".equals( InternalTaxonomyAudit.describe( "Carnivora", null, "" ) )
                || !"(unnamed) [1]".equals( InternalTaxonomyAudit.describe( null, null, "1" ) ) ) {
            return fail( "describe" );
        }
        return true;
    }

    /** By id when both carry one -- a same name with another id is NOT this taxon -- else by name. */
    private static boolean testLevelOf() {
        final TaxonLineage tl = mammal( "Felidae", "9681" );
        if ( InternalTaxonomyAudit.levelOf( tl, "Carnivora", "33554" ) == null
                || InternalTaxonomyAudit.levelOf( tl, "Carnivora", null ) == null
                || InternalTaxonomyAudit.levelOf( tl, "carnivora", null ) == null
                || InternalTaxonomyAudit.levelOf( tl, "Felidae", "9681" ) == null ) {
            return fail( "levelOf must find an ancestor or the taxon itself by id or by name" );
        }
        if ( InternalTaxonomyAudit.levelOf( tl, "Carnivora", "99999" ) != null ) {
            return fail( "a same-named level with ANOTHER id is not this taxon" );
        }
        if ( InternalTaxonomyAudit.levelOf( tl, "Rodentia", null ) != null || InternalTaxonomyAudit.levelOf( tl, "Rodentia", "9989" ) != null ) {
            return fail( "a taxon not in the lineage" );
        }
        return true;
    }

    // ---- fixture helpers

    private static TaxonLineage bacterium( final String genus, final String genus_id, final String order, final String order_id ) {
        final List<Ancestor> anc = new ArrayList<Ancestor>();
        anc.add( new Ancestor( "cellular organisms", "no rank", "131567" ) );
        anc.add( new Ancestor( "Bacteria", "superkingdom", "2" ) );
        anc.add( new Ancestor( order, "order", order_id ) );
        return new TaxonLineage( genus_id, "genus", genus, null, anc );
    }

    private static TaxonLineage mammal( final String family, final String family_id ) {
        final List<Ancestor> anc = new ArrayList<Ancestor>();
        anc.add( new Ancestor( "Eukaryota", "superkingdom", "2759" ) );
        anc.add( new Ancestor( "Mammalia", "class", "40674" ) );
        anc.add( new Ancestor( "Carnivora", "order", "33554" ) );
        return new TaxonLineage( family_id, "family", family, null, anc );
    }

    private static PhylogenyNode tip( final String name, final Map<PhylogenyNode, TaxonLineage> map, final TaxonLineage tl ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        if ( tl != null ) {
            map.put( n, tl );
        }
        return n;
    }

    private static PhylogenyNode internal( final String sci, final String rank, final String id ) {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        if ( rank != null ) {
            try {
                t.setRank( rank );
            }
            catch ( final Exception e ) {
                throw new RuntimeException( e );
            }
        }
        if ( id != null ) {
            t.setIdentifier( new Identifier( id, "ncbi" ) );
        }
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    /** Cache-only. */
    private static final class FakeService implements TaxonomicLineageService {

        private final Map<String, TaxonLineage> _cache = new HashMap<String, TaxonLineage>();

        void know( final String key, final TaxonLineage tl ) {
            _cache.put( key.toLowerCase( Locale.ROOT ), tl );
        }

        void empty( final String key ) {
            _cache.put( key.toLowerCase( Locale.ROOT ), TaxonLineage.EMPTY );
        }

        @Override
        public TaxonLineage lineageOf( final String taxon ) {
            return ( taxon == null ) ? null : _cache.get( taxon.trim().toLowerCase( Locale.ROOT ) );
        }

        @Override
        public TaxonLineage fetch( final String taxon ) {
            throw new UnsupportedOperationException( "the audit never fetches" );
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "InternalTaxonomyAudit test failed: " + msg );
        return false;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }
}
