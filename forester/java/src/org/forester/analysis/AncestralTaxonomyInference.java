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
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonLineage.Ancestor;

/**
 * Assigns each internal node the deepest taxon shared by ALL its descendant tips, computed from the tips'
 * resolved {@link TaxonLineage}s (from the unified taxonomy service / cache -- see Spine A/B). The inferred
 * internal taxon carries a scientific name, a rank, and (when known) an NCBI tax-id, so it is a first-class
 * taxon that "Colorize / Annotate Clades by Rank" can use directly.
 *
 * <p>The computation is a pure, headless-testable function of the tree shape + a tip&rarr;lineage map. It is
 * <b>n-ary</b> (a multifurcation intersects all its children, not just the first two) and <b>conservative</b>:
 * a node with even one unresolved descendant tip infers nothing (an empty shared path intersects to empty),
 * which is why the caller resolves the tips online first. Existing (hand-curated) internal taxonomies are
 * preserved unless {@code overwrite} is set.
 *
 * <p>Two deliberate properties: (1) it does NOT clear a child's label when it equals its parent's (the old code
 * did, but that mutated the child's data) -- nested clades that share their deepest taxon simply show it on both,
 * which is honest and left to the renderer to de-duplicate; (2) with {@code overwrite} on, an existing internal
 * taxon whose shared path is empty (an un-inferrable node) is left in place rather than blanked.
 */
public final class AncestralTaxonomyInference {

    private AncestralTaxonomyInference() {
        // utility class
    }

    /**
     * Writes an inferred taxonomy onto each internal node of {@code phy}, in place.
     *
     * @param phy          the (rooted) tree to enrich
     * @param tip_lineages each external node's resolved lineage; an absent/{@link TaxonLineage#isEmpty() empty}
     *                     entry marks a tip that could not be resolved (its ancestors then infer nothing)
     * @param overwrite    when false, an internal node that already carries a taxonomy is left untouched; when
     *                     true, every internal node with a shared path is (re)written
     * @return counts of what happened (never null)
     */
    public static InferenceResult inferInternalTaxonomies( final Phylogeny phy,
                                                           final Map<PhylogenyNode, TaxonLineage> tip_lineages,
                                                           final boolean overwrite ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return new InferenceResult( 0, 0, 0, 0 );
        }
        // Each node's "shared path" = the levels (root->deepest) common to ALL of its descendant tips.
        // Computed once postorder (children are visited before their parent).
        final Map<PhylogenyNode, List<Ancestor>> shared = new HashMap<PhylogenyNode, List<Ancestor>>();
        int resolved_tips = 0;
        int unresolved_tips = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( node.isExternal() ) {
                final TaxonLineage tl = ( tip_lineages == null ) ? null : tip_lineages.get( node );
                shared.put( node, tipPath( tl ) );
                // "resolved" == has at least one ANCESTOR to intersect; an own-only taxon (no ancestors) can seed no
                // shared ancestor for its parent, so it is reported as unresolved (resolve it online to fix that)
                if ( ( tl == null ) || tl.getAncestors().isEmpty() ) {
                    unresolved_tips++;
                }
                else {
                    resolved_tips++;
                }
            }
            else {
                List<Ancestor> common = null;
                for( final PhylogenyNode child : node.getDescendants() ) {
                    final List<Ancestor> child_path = shared.get( child );
                    common = ( common == null ) ? child_path : commonPrefix( common, child_path );
                    if ( common.isEmpty() ) {
                        break; // an unresolved / disagreeing descendant -> this node (and its ancestors) infer nothing
                    }
                }
                shared.put( node, ( common == null ) ? Collections.<Ancestor> emptyList() : common );
            }
        }
        int assigned = 0;
        int skipped_existing = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( node.isExternal() ) {
                continue;
            }
            final List<Ancestor> path = shared.get( node );
            if ( ( path == null ) || path.isEmpty() ) {
                continue;
            }
            if ( node.getNodeData().isHasTaxonomy() && !overwrite ) {
                skipped_existing++; // preserve hand-curated internal taxa unless overwriting
                continue;
            }
            node.getNodeData().setTaxonomy( toTaxonomy( path ) );
            assigned++;
        }
        return new InferenceResult( assigned, skipped_existing, resolved_tips, unresolved_tips );
    }

    /** A tip's lineage as an ordered root&rarr;deepest path: its ancestors (root&rarr;parent) then its own taxon. */
    private static List<Ancestor> tipPath( final TaxonLineage tl ) {
        final List<Ancestor> path = new ArrayList<Ancestor>();
        if ( ( tl == null ) || tl.isEmpty() ) {
            return path;
        }
        path.addAll( tl.getAncestors() ); // root -> parent
        if ( !ForesterUtil.isEmpty( tl.getScientificName() ) ) {
            path.add( new Ancestor( tl.getScientificName(), tl.getRank(), tl.getTaxId() ) ); // the taxon itself
        }
        return path;
    }

    /** The common leading run of two root&rarr;deepest paths (levels matched by {@link #sameLevel}). */
    private static List<Ancestor> commonPrefix( final List<Ancestor> a, final List<Ancestor> b ) {
        final List<Ancestor> out = new ArrayList<Ancestor>();
        final int n = Math.min( a.size(), b.size() );
        for( int i = 0; i < n; i++ ) {
            if ( sameLevel( a.get( i ), b.get( i ) ) ) {
                out.add( a.get( i ) );
            }
            else {
                break;
            }
        }
        return out;
    }

    /** Two lineage levels are the same taxon: by NCBI tax-id when both carry one, else by name (case-insensitive). */
    private static boolean sameLevel( final Ancestor x, final Ancestor y ) {
        final String xi = x.getTaxId();
        final String yi = y.getTaxId();
        if ( !ForesterUtil.isEmpty( xi ) && !ForesterUtil.isEmpty( yi ) ) {
            return xi.equals( yi );
        }
        final String xn = x.getName();
        return !ForesterUtil.isEmpty( xn ) && xn.equalsIgnoreCase( y.getName() );
    }

    /** Build the internal-node taxonomy from the DEEPEST shared level (+ the full shared path as its lineage). */
    private static Taxonomy toTaxonomy( final List<Ancestor> path ) {
        final Ancestor deepest = path.get( path.size() - 1 );
        final Taxonomy tax = new Taxonomy();
        tax.setScientificName( deepest.getName() );
        if ( !ForesterUtil.isEmpty( deepest.getRank() ) ) {
            try {
                tax.setRank( deepest.getRank().toLowerCase( Locale.ROOT ) );
            }
            catch ( final PhyloXmlDataFormatException e ) {
                // a non-Linnaean rank ("no rank" / "clade" / unknown) -- keep the name, leave the rank unset
            }
        }
        if ( !ForesterUtil.isEmpty( deepest.getTaxId() ) ) {
            tax.setIdentifier( new Identifier( deepest.getTaxId(), "ncbi" ) );
        }
        final List<String> names = new ArrayList<String>();
        for( final Ancestor a : path ) {
            names.add( a.getName() );
        }
        tax.setLineage( names );
        return tax;
    }

    /** Immutable summary of an inference pass. */
    public static final class InferenceResult {

        private final int _assigned;
        private final int _skipped_existing;
        private final int _resolved_tips;
        private final int _unresolved_tips;

        InferenceResult( final int assigned,
                         final int skipped_existing,
                         final int resolved_tips,
                         final int unresolved_tips ) {
            _assigned = assigned;
            _skipped_existing = skipped_existing;
            _resolved_tips = resolved_tips;
            _unresolved_tips = unresolved_tips;
        }

        /** Internal nodes that received an inferred taxonomy. */
        public int getAssigned() {
            return _assigned;
        }

        /** Internal nodes left untouched because they already carried a taxonomy and {@code overwrite} was false. */
        public int getSkippedExisting() {
            return _skipped_existing;
        }

        /** Tips that carried an ancestor lineage the inference could intersect. */
        public int getResolvedTips() {
            return _resolved_tips;
        }

        /** Tips with no ancestor lineage (empty or own-taxon-only) -- each blanks the inference of every node
         *  above it. Distinct from {@code TreePanelUtil.tipsWithoutLineage}, which is the set to fetch online. */
        public int getUnresolvedTips() {
            return _unresolved_tips;
        }
    }
}
