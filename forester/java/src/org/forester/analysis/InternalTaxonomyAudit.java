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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeSet;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.util.TaxonomyUtil;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonLineage.Ancestor;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * Checks every INTERNAL node's {@code <taxonomy>} against two references and proposes a correction for each
 * finding: (1) the lineages of the tips below it -- an annotation must name a taxon those lineages contain, else it
 * is a wrong taxon for that clade (the real case: a tree whose "Bacteria" nodes had been resolved to the stick-insect
 * genus <i>Bacteria</i> Latreille, 17 of them), and the tips also tell the taxon's RANK; (2) NCBI Taxonomy by the
 * node's own tax-id, when cached -- which tells the taxon's name and rank authoritatively. The proposal for a wrong
 * taxon is the deepest taxon ALL tips below share ({@link AncestralTaxonomyInference#commonTaxonOf}); for a wrong
 * rank or name, the node's taxon with the reference's rank or name. Pure: reads the tree and the service's cache,
 * writes nothing until {@link #apply}, which the caller wraps in an undo checkpoint and a provenance sentence.
 */
public final class InternalTaxonomyAudit {

    /** What the check found for one internal annotation. */
    public enum Verdict {
        /** The taxon is in every tip's lineage (or the database confirms it) with the same rank and name. */
        OK( "consistent", false ),
        /** The taxon is not in the lineage of ANY tip below the node: a wrong taxon for this clade. */
        WRONG_TAXON( "wrong taxon: no tip below it belongs to it", true ),
        /** Some tips below belong to the taxon, others do not: the annotation sits too deep in the tree. */
        NOT_SHARED_BY_ALL_TIPS( "not shared by all tips below it", true ),
        /** The taxon is right, its rank is not (the tips' lineages or NCBI give another). */
        WRONG_RANK( "wrong rank", true ),
        /** The taxon is right by id, its scientific name is not NCBI's. */
        WRONG_NAME( "name differs from NCBI's", true ),
        /** The database does not know the taxon (an unknown or ambiguous name). */
        UNKNOWN_TO_DATABASE( "unknown to the database", false ),
        /** Neither the tips below nor the database could say anything (nothing cached, no tip lineages). */
        NOT_CHECKED( "not checked: no tip lineage and nothing cached", false );

        private final String  _label;
        private final boolean _correctable;

        Verdict( final String label, final boolean correctable ) {
            _label = label;
            _correctable = correctable;
        }

        public String getLabel() {
            return _label;
        }

        /** True when the finding carries a proposed correction. */
        public boolean isCorrectable() {
            return _correctable;
        }
    }

    /** One internal annotation, its verdict, the evidence, and the proposed taxonomy (null when none). */
    public static final class Finding {

        private final PhylogenyNode _node;
        private final String        _current;
        private final Verdict       _verdict;
        private final String        _detail;
        private final String        _proposed_name;
        private final String        _proposed_rank;
        private final String        _proposed_id;

        Finding( final PhylogenyNode node, final String current, final Verdict verdict, final String detail,
                 final String proposed_name, final String proposed_rank, final String proposed_id ) {
            _node = node;
            _current = current;
            _verdict = verdict;
            _detail = detail;
            _proposed_name = proposed_name;
            _proposed_rank = proposed_rank;
            _proposed_id = proposed_id;
        }

        public PhylogenyNode getNode() {
            return _node;
        }

        /** The annotation as it is: "rank name [id]". */
        public String getCurrent() {
            return _current;
        }

        public Verdict getVerdict() {
            return _verdict;
        }

        public String getDetail() {
            return _detail;
        }

        /** The proposed annotation, "rank name [id]", or "" when there is none. */
        public String getProposed() {
            return ( _proposed_name == null ) ? "" : describe( _proposed_name, _proposed_rank, _proposed_id );
        }

        public boolean hasProposal() {
            return _proposed_name != null;
        }

        String getProposedName() {
            return _proposed_name;
        }

        String getProposedRank() {
            return _proposed_rank;
        }

        String getProposedId() {
            return _proposed_id;
        }
    }

    /** "order Carnivora [33554]"; the rank or id left out when absent. */
    public static String describe( final String name, final String rank, final String id ) {
        final StringBuilder sb = new StringBuilder();
        if ( !ForesterUtil.isEmpty( rank ) ) {
            sb.append( rank ).append( ' ' );
        }
        sb.append( ForesterUtil.isEmpty( name ) ? "(unnamed)" : name );
        if ( !ForesterUtil.isEmpty( id ) ) {
            sb.append( " [" ).append( id ).append( ']' );
        }
        return sb.toString();
    }

    /**
     * Audits every internal node that carries a taxonomy. {@code tip_lineages} is each tip's resolved lineage (from
     * the cache or the tree, as {@code TreePanelUtil.tipLineages} builds it; an absent or empty entry = unresolved);
     * {@code service} is consulted cache-only for the node's own taxon (may be null). Findings come in preorder.
     */
    public static List<Finding> audit( final Phylogeny phy, final Map<PhylogenyNode, TaxonLineage> tip_lineages,
                                       final TaxonomicLineageService service ) {
        final List<Finding> out = new ArrayList<Finding>();
        if ( ( phy == null ) || phy.isEmpty() ) {
            return out;
        }
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( node.isExternal() || !node.getNodeData().isHasTaxonomy() ) {
                continue;
            }
            final Taxonomy tax = node.getNodeData().getTaxonomy();
            final String name = tax.getScientificName();
            final String id = ncbiId( tax );
            if ( ForesterUtil.isEmpty( name ) && ForesterUtil.isEmpty( id ) ) {
                continue; // nothing to check (a code- or common-name-only annotation is not resolvable here)
            }
            out.add( check( node, tax, name, id, tip_lineages, service ) );
        }
        return out;
    }

    private static Finding check( final PhylogenyNode node, final Taxonomy tax, final String name, final String id,
                                  final Map<PhylogenyNode, TaxonLineage> tip_lineages,
                                  final TaxonomicLineageService service ) {
        final String current = describe( name, tax.getRank(), id );
        // --- the database's word on the node's own taxon (cache-only): by id first, else by name
        TaxonLineage db = null;
        boolean db_asked = false;
        if ( service != null ) {
            if ( !ForesterUtil.isEmpty( id ) ) {
                db = service.lineageOf( id );
                db_asked = true;
            }
            if ( ( db == null ) && !ForesterUtil.isEmpty( name ) ) {
                db = service.lineageOf( name );
                db_asked = true;
            }
        }
        final boolean db_known = ( db != null ) && !db.isEmpty();
        // --- the tips' word: which of the resolved tips below carry the taxon in their lineage, and at what rank
        int resolved = 0;
        int carrying = 0;
        String rank_from_tips = null;
        String id_from_tips = null;
        for( final PhylogenyNode tip : node.getAllExternalDescendants() ) {
            final TaxonLineage tl = ( tip_lineages == null ) ? null : tip_lineages.get( tip );
            if ( ( tl == null ) || tl.isEmpty() ) {
                continue;
            }
            resolved++;
            final Ancestor level = levelOf( tl, name, id );
            if ( level != null ) {
                carrying++;
                if ( ( rank_from_tips == null ) && isLinnaean( level.getRank() ) ) {
                    rank_from_tips = level.getRank().toLowerCase( Locale.ROOT );
                }
                if ( ( id_from_tips == null ) && !ForesterUtil.isEmpty( level.getTaxId() ) ) {
                    id_from_tips = level.getTaxId();
                }
            }
        }
        if ( ( resolved > 0 ) && ( carrying == 0 ) ) {
            // no tip below belongs to this taxon: the clade is annotated with the wrong taxon. Propose the deepest
            // taxon ALL resolved tips share (null when they share none or some tip is unresolved).
            final Ancestor common = AncestralTaxonomyInference.commonTaxonOf( node, tip_lineages );
            final String detail = "none of the " + resolved + " resolved tip" + ( resolved == 1 ? "" : "s" ) + " below it belongs to it"
                    + ( db_known ? ( "; NCBI knows " + current + " as " + describeDb( db ) ) : "" )
                    + ( ( common == null ) ? "; the tips share no taxon that could replace it" : "" );
            return new Finding( node, current, Verdict.WRONG_TAXON, detail, ( common == null ) ? null : common.getName(),
                                ( common == null ) ? null : linnaeanOrNull( common.getRank() ),
                                ( common == null ) ? null : common.getTaxId() );
        }
        if ( ( resolved > 0 ) && ( carrying < resolved ) ) {
            final Ancestor common = AncestralTaxonomyInference.commonTaxonOf( node, tip_lineages );
            final String detail = "only " + carrying + " of the " + resolved + " resolved tips below it belong to it"
                    + ( ( common == null ) ? "; the tips share no taxon that could replace it" : "" );
            return new Finding( node, current, Verdict.NOT_SHARED_BY_ALL_TIPS, detail, ( common == null ) ? null : common.getName(),
                                ( common == null ) ? null : linnaeanOrNull( common.getRank() ),
                                ( common == null ) ? null : common.getTaxId() );
        }
        // the taxon is right for the clade (every resolved tip carries it), or no tip could say. Now rank and name.
        final String reference_rank = db_known ? db.ownRankOrNull() : rank_from_tips;
        final String reference_name = ( db_known && !ForesterUtil.isEmpty( id ) ) ? db.getScientificName() : null;
        final String reference_id = ForesterUtil.isEmpty( id ) ? ( db_known ? db.getTaxId() : id_from_tips ) : id;
        final String rank = ForesterUtil.isEmpty( tax.getRank() ) ? null : tax.getRank().toLowerCase( Locale.ROOT );
        if ( ( reference_rank != null ) && ( rank != null ) && !sameRank( rank, reference_rank ) ) {
            return new Finding( node, current, Verdict.WRONG_RANK,
                                "the " + ( db_known ? "database" : "tips' lineages" ) + " give the rank " + reference_rank,
                                name, reference_rank, reference_id );
        }
        if ( ( reference_name != null ) && !ForesterUtil.isEmpty( name ) && !reference_name.equalsIgnoreCase( name ) ) {
            return new Finding( node, current, Verdict.WRONG_NAME, "NCBI names id " + id + " \"" + reference_name + "\"",
                                reference_name, ( rank == null ) ? reference_rank : rank, id );
        }
        if ( ( resolved == 0 ) && !db_known ) {
            return new Finding( node, current, db_asked && ( db != null ) ? Verdict.UNKNOWN_TO_DATABASE : Verdict.NOT_CHECKED,
                                db_asked && ( db != null ) ? "the database has no taxon for it" : null, null, null, null );
        }
        return new Finding( node, current, Verdict.OK, ( resolved > 0 ) ? ( "in the lineage of all " + resolved + " resolved tip"
                + ( resolved == 1 ? "" : "s" ) + " below it" ) : "confirmed by the database", null, null, null );
    }

    private static String describeDb( final TaxonLineage db ) {
        return describe( db.getScientificName(), db.ownRankOrNull(), db.getTaxId() );
    }

    /** The level of {@code tl} (an ancestor or the taxon itself) that is this taxon: by id when both carry one,
     *  else by name; null when the lineage does not contain it. */
    static Ancestor levelOf( final TaxonLineage tl, final String name, final String id ) {
        final List<Ancestor> levels = new ArrayList<Ancestor>( tl.getAncestors() );
        if ( !ForesterUtil.isEmpty( tl.getScientificName() ) ) {
            levels.add( new Ancestor( tl.getScientificName(), tl.getRank(), tl.getTaxId() ) );
        }
        for( final Ancestor a : levels ) {
            if ( !ForesterUtil.isEmpty( id ) && !ForesterUtil.isEmpty( a.getTaxId() ) ) {
                if ( id.equals( a.getTaxId() ) ) {
                    return a;
                }
                continue; // both have ids and they differ: not this taxon, whatever the names say
            }
            if ( !ForesterUtil.isEmpty( name ) && name.equalsIgnoreCase( a.getName() ) ) {
                return a;
            }
        }
        return null;
    }

    private static boolean isLinnaean( final String rank ) {
        return !ForesterUtil.isEmpty( rank ) && !"no rank".equalsIgnoreCase( rank ) && !"clade".equalsIgnoreCase( rank );
    }

    private static String linnaeanOrNull( final String rank ) {
        return isLinnaean( rank ) ? rank.toLowerCase( Locale.ROOT ) : null;
    }

    /** NCBI's "superkingdom" and phyloXML's "domain" are the same level. */
    static boolean sameRank( final String a, final String b ) {
        final String x = a.toLowerCase( Locale.ROOT );
        final String y = b.toLowerCase( Locale.ROOT );
        if ( x.equals( y ) ) {
            return true;
        }
        return ( x.equals( "domain" ) && y.equals( "superkingdom" ) ) || ( x.equals( "superkingdom" ) && y.equals( "domain" ) );
    }

    private static String ncbiId( final Taxonomy tax ) {
        if ( ( tax.getIdentifier() != null ) && !ForesterUtil.isEmpty( tax.getIdentifier().getValue() )
                && "ncbi".equalsIgnoreCase( tax.getIdentifier().getProvider() ) ) {
            return tax.getIdentifier().getValue();
        }
        return null;
    }

    /**
     * The lookup keys (an NCBI id, else a scientific name) of the internal taxa the service has NOT cached -- what a
     * caller fetches off the EDT before {@link #audit} so the database can have its word. Sorted, distinct.
     */
    public static TreeSet<String> internalTaxaToLookUp( final Phylogeny phy, final TaxonomicLineageService service ) {
        final TreeSet<String> keys = new TreeSet<String>();
        if ( ( phy == null ) || phy.isEmpty() ) {
            return keys;
        }
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( node.isExternal() || !node.getNodeData().isHasTaxonomy() ) {
                continue;
            }
            final Taxonomy tax = node.getNodeData().getTaxonomy();
            final String key = !ForesterUtil.isEmpty( ncbiId( tax ) ) ? ncbiId( tax ) : tax.getScientificName();
            if ( ForesterUtil.isEmpty( key ) ) {
                continue;
            }
            if ( ( service != null ) && ( service.lineageOf( key ) != null ) ) {
                continue;
            }
            keys.add( key );
        }
        return keys;
    }

    /**
     * Writes the proposed correction of every finding whose verdict is in {@code which} onto its node: a wrong
     * rank or name is corrected IN the node's taxonomy (its other fields kept); a wrong or too-deep taxon is REPLACED
     * by the proposed one (name + rank + id, via {@link TaxonomyUtil#buildNcbiTaxonomy}). Findings without a proposal
     * are skipped. Returns the number of nodes changed. The caller checkpoints undo first.
     */
    public static int apply( final List<Finding> findings, final java.util.Set<Verdict> which ) {
        int changed = 0;
        if ( findings == null ) {
            return 0;
        }
        for( final Finding f : findings ) {
            if ( !f.hasProposal() || ( which != null && !which.contains( f.getVerdict() ) ) ) {
                continue;
            }
            final Taxonomy tax = f.getNode().getNodeData().getTaxonomy();
            switch ( f.getVerdict() ) {
                case WRONG_RANK:
                    if ( setRank( tax, f.getProposedRank() ) ) {
                        if ( ForesterUtil.isEmpty( ncbiId( tax ) ) && !ForesterUtil.isEmpty( f.getProposedId() ) ) {
                            tax.setIdentifier( new Identifier( f.getProposedId(), "ncbi" ) );
                        }
                        changed++;
                    }
                    break;
                case WRONG_NAME:
                    tax.setScientificName( f.getProposedName() );
                    changed++;
                    break;
                case WRONG_TAXON:
                case NOT_SHARED_BY_ALL_TIPS:
                    f.getNode().getNodeData().setTaxonomy( TaxonomyUtil.buildNcbiTaxonomy( f.getProposedName(),
                                                                                             f.getProposedRank(),
                                                                                             f.getProposedId() ) );
                    changed++;
                    break;
                default:
                    break;
            }
        }
        return changed;
    }

    private static boolean setRank( final Taxonomy tax, final String rank ) {
        if ( ForesterUtil.isEmpty( rank ) ) {
            return false;
        }
        try {
            tax.setRank( rank.toLowerCase( Locale.ROOT ) );
            return true;
        }
        catch ( final PhyloXmlDataFormatException e ) {
            return false; // a rank outside the phyloXML vocabulary: leave it
        }
    }

    /** The findings per verdict, in enum order, for a summary. */
    public static Map<Verdict, Integer> countByVerdict( final List<Finding> findings ) {
        final Map<Verdict, Integer> m = new LinkedHashMap<Verdict, Integer>();
        for( final Verdict v : Verdict.values() ) {
            m.put( v, 0 );
        }
        for( final Finding f : findings ) {
            m.merge( f.getVerdict(), 1, Integer::sum );
        }
        return m;
    }

    /** "Checked 190 internal annotations: 170 consistent, 10 wrong taxon, ..." -- the verdicts with a count,
     *  then how many corrections are proposed. */
    public static String summary( final List<Finding> findings ) {
        final StringBuilder sb = new StringBuilder( "Checked " ).append( findings.size() ).append( " internal annotation" )
                .append( findings.size() == 1 ? "" : "s" );
        int proposals = 0;
        final List<String> parts = new ArrayList<String>();
        for( final Map.Entry<Verdict, Integer> e : countByVerdict( findings ).entrySet() ) {
            if ( e.getValue() > 0 ) {
                parts.add( e.getValue() + " " + e.getKey().getLabel() );
            }
        }
        for( final Finding f : findings ) {
            if ( f.hasProposal() ) {
                proposals++;
            }
        }
        if ( !parts.isEmpty() ) {
            sb.append( ": " ).append( String.join( "; ", parts ) );
        }
        sb.append( ". " ).append( proposals ).append( " correction" ).append( proposals == 1 ? "" : "s" ).append( " proposed." );
        return sb.toString();
    }

    /** The provenance sentence a caller appends after {@link #apply}. */
    public static String provenance( final Phylogeny phy, final int changed ) {
        final String name = ForesterUtil.isEmpty( phy.getName() ) ? "" : phy.getName();
        return "Used the internal-taxonomy check (tip lineages and NCBI Taxonomy) to correct " + changed + " internal node taxonom"
                + ( changed == 1 ? "y" : "ies" ) + " (wrong taxa replaced by the deepest taxon their tips share; ranks and names"
                + " corrected) in tree named \"" + name + "\" with " + phy.getNumberOfExternalNodes() + " tips.";
    }

    private InternalTaxonomyAudit() {
    }
}
