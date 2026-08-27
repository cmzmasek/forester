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
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;
import org.forester.util.TaxonomyUtil;
import org.forester.ws.seqdb.TaxonLineage;

/**
 * Builds an induced NCBI-taxonomy "species tree" from a set of already-resolved tip lineages -- for APPROXIMATE
 * gene-tree/species-tree reconciliation (GSDIR) when the user has no curated species tree. Each tip's lineage
 * (root &rarr; ... &rarr; its own taxon) is merged into a trie: shared higher taxa become internal nodes and the
 * tips' own taxa become the leaves. Nodes are keyed by NCBI tax-id where present (else by lower-cased scientific
 * name), so two distinct taxa that happen to share a name at some rank are NOT merged. Because the input lineages
 * can come from two sources with different id-completeness -- online-resolved lineages carry tax-ids, a tip's stored
 * in-tree lineage is name-only -- the builder first learns a name&rarr;tax-id mapping so a taxon that is name-only in
 * one tip still merges with the SAME taxon that carries a tax-id in another (rather than splitting into two nodes).
 * Relationships the taxonomy leaves unresolved appear as POLYTOMIES, which GSDIR tolerates and which keep
 * duplication calls conservative.
 *
 * <p><b>This is an approximation:</b> the NCBI taxonomy is a classification, not a phylogeny (polytomies, and the
 * occasional non-monophyletic rank), so the reconciliation it feeds is a quick, zero-setup estimate rather than a
 * substitute for a curated species tree -- the caller surfaces that caveat.
 *
 * <p>Pure and GUI-free: it takes an already-resolved {@code {tip -> TaxonLineage}} map (built cache-only by
 * {@code TreePanelUtil.tipLineages}), so it is unit-testable with no network. Each node's {@code <taxonomy>} is
 * built via {@link TaxonomyUtil#buildNcbiTaxonomy} (Linnaean rank kept, non-Linnaean dropped, NCBI id when known),
 * so a gene-tree tip links to its species-tree leaf on the same taxonomic basis GSDI already uses.
 */
public final class TaxonomySpeciesTreeBuilder {

    /** The built species tree plus how many gene tips had no usable lineage (so were left out of it). */
    public static final class Result {

        private final Phylogeny _species_tree;
        private final int       _species_count;
        private final int       _tips_without_lineage;

        Result( final Phylogeny species_tree, final int species_count, final int tips_without_lineage ) {
            _species_tree = species_tree;
            _species_count = species_count;
            _tips_without_lineage = tips_without_lineage;
        }

        /** The induced taxonomy species tree (may be empty if no tip had a usable lineage). */
        public Phylogeny getSpeciesTree() {
            return _species_tree;
        }

        /** Number of distinct species (external nodes) in the built tree. */
        public int getSpeciesCount() {
            return _species_count;
        }

        /** Gene tips whose lineage was empty/unusable -- they are not represented in the species tree. */
        public int getTipsWithoutLineage() {
            return _tips_without_lineage;
        }
    }

    private TaxonomySpeciesTreeBuilder() {
        // not instantiable
    }

    /** Builds the induced taxonomy species tree from the per-tip lineages. Tips whose lineage is empty/unusable are
     *  skipped (counted in {@link Result#getTipsWithoutLineage()}). */
    public static Result build( final Map<PhylogenyNode, TaxonLineage> tip_lineages ) {
        // Pass 1: collect the usable lineage paths and learn a name -> tax-id mapping. The lineages can arrive from
        // TWO sources with different id-completeness -- online-resolved lineages carry NCBI tax-ids, a tip's stored
        // in-tree lineage is name-only -- so the SAME taxon (e.g. Mammalia) may appear with an id in one tip and
        // name-only in another. Learning the id here lets the name-only occurrence key on the same tax-id, so a shared
        // taxon MERGES into one node instead of splitting (which would corrupt the reconciliation / duplicate a leaf).
        final List<List<Step>> paths = new ArrayList<List<Step>>();
        final Map<String, String> name_to_tax_id = new HashMap<String, String>();
        int without_lineage = 0;
        if ( tip_lineages != null ) {
            for( final TaxonLineage lineage : tip_lineages.values() ) {
                final List<Step> path = pathOf( lineage );
                if ( path.isEmpty() ) {
                    ++without_lineage;
                    continue;
                }
                for( final Step step : path ) {
                    if ( !ForesterUtil.isEmpty( step._tax_id ) && !name_to_tax_id.containsKey( step._name_key ) ) {
                        name_to_tax_id.put( step._name_key, step._tax_id );
                    }
                }
                paths.add( path );
            }
        }
        // Pass 2: merge each path into a trie. Key each taxon by its tax-id -- its own, or the one another occurrence
        // of the same name supplied in pass 1 -- else by the (lower-cased) name; so two DISTINCT taxa sharing a name
        // but with different tax-ids (a homonym) still stay separate. Insertion-ordered for a deterministic tree.
        final Map<String, PhylogenyNode> by_key = new LinkedHashMap<String, PhylogenyNode>();
        for( final List<Step> path : paths ) {
            PhylogenyNode parent = null;
            final Set<String> seen_in_path = new HashSet<String>();
            for( final Step step : path ) {
                final String tax_id = ForesterUtil.isEmpty( step._tax_id ) ? name_to_tax_id.get( step._name_key )
                                                                           : step._tax_id;
                final String key = ForesterUtil.isEmpty( tax_id ) ? ( "n:" + step._name_key ) : tax_id;
                // a taxon can appear only once per lineage; ignore a (malformed) repeat -- no re-attach, no cycle
                if ( !seen_in_path.add( key ) ) {
                    continue;
                }
                PhylogenyNode node = by_key.get( key );
                if ( node == null ) {
                    node = new PhylogenyNode();
                    node.getNodeData().setTaxonomy( TaxonomyUtil.buildNcbiTaxonomy( step._name, step._rank, tax_id ) );
                    node.setName( step._name ); // also as the node name, so it shows even with taxonomy display off
                    by_key.put( key, node );
                }
                // attach to its parent exactly once (the NCBI taxonomy is a tree, so a taxon has one parent)
                if ( ( parent != null ) && ( node != parent ) && ( node.getParent() == null ) ) {
                    parent.addAsChild( node );
                }
                parent = node;
            }
        }
        final Phylogeny phy = new Phylogeny();
        // the roots are the created nodes that never got a parent; join several disjoint tops (e.g. Viruses vs
        // cellular organisms) under one synthetic root so the species tree is a single rooted tree for GSDIR
        final List<PhylogenyNode> roots = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNode n : by_key.values() ) {
            if ( n.getParent() == null ) {
                roots.add( n );
            }
        }
        if ( roots.size() == 1 ) {
            phy.setRoot( roots.get( 0 ) );
        }
        else if ( roots.size() > 1 ) {
            final PhylogenyNode synthetic = new PhylogenyNode();
            synthetic.setName( "root" );
            for( final PhylogenyNode r : roots ) {
                synthetic.addAsChild( r );
            }
            phy.setRoot( synthetic );
        }
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        final int species = phy.isEmpty() ? 0 : phy.getNumberOfExternalNodes();
        return new Result( phy, species, without_lineage );
    }

    /** One taxon on a lineage path: name, rank, tax-id, plus the lower-cased name used to unify name-only occurrences
     *  with id-bearing ones of the same taxon (see {@link #build}). */
    private static final class Step {

        private final String _name;
        private final String _rank;
        private final String _tax_id;
        private final String _name_key;

        private Step( final String name, final String rank, final String tax_id ) {
            _name = name;
            _rank = rank;
            _tax_id = tax_id;
            _name_key = name.toLowerCase( Locale.ROOT );
        }
    }

    /** The lineage as an ordered root&rarr;taxon path of named steps; empty if the lineage carries no usable taxon. */
    private static List<Step> pathOf( final TaxonLineage lineage ) {
        final List<Step> path = new ArrayList<Step>();
        if ( ( lineage == null ) || lineage.isEmpty() ) {
            return path;
        }
        for( final TaxonLineage.Ancestor a : lineage.getAncestors() ) {
            if ( !ForesterUtil.isEmpty( a.getName() ) ) {
                path.add( new Step( a.getName(), a.getRank(), a.getTaxId() ) );
            }
        }
        if ( !ForesterUtil.isEmpty( lineage.getScientificName() ) ) {
            path.add( new Step( lineage.getScientificName(), lineage.getRank(), lineage.getTaxId() ) );
        }
        return path;
    }
}
