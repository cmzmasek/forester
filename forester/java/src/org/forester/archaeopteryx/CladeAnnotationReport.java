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

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.forester.archaeopteryx.TreePanelUtil.RankTaxon;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * The audit of a rank annotation: for every rank asked for, one row per TIP -- its taxon at that rank, the
 * {@link TaxonPlacement.Source source} that placed it (or why nothing did) and a note (the ancestor consulted, the
 * query sent, the candidates of an ambiguous name, the tree annotation the database overrode) -- plus one row per
 * INTERNAL annotation at that rank that the grouping had to IGNORE because tips under it are placed elsewhere, and
 * a summary line per rank. This is what makes an incomplete or inconsistent annotation visible instead of silently
 * merged. Pure (cache-only through the service), so it can be built before anything is drawn or written, and
 * exported as TSV.
 */
final class CladeAnnotationReport implements TableReport {

    /** One line of the report. {@code internal} marks an ignored internal annotation rather than a tip. */
    static final class Row {

        private final String  _rank;
        private final String  _node;
        private final String  _taxon;
        private final String  _source;
        private final String  _note;
        private final boolean _internal;

        Row( final String rank, final String node, final String taxon, final String source, final String note,
             final boolean internal ) {
            _rank = rank;
            _node = node;
            _taxon = ( taxon == null ) ? "" : taxon;
            _source = source;
            _note = ( note == null ) ? "" : note;
            _internal = internal;
        }

        String getRank() {
            return _rank;
        }

        String getNode() {
            return _node;
        }

        String getTaxon() {
            return _taxon;
        }

        String getSource() {
            return _source;
        }

        String getNote() {
            return _note;
        }

        boolean isInternal() {
            return _internal;
        }
    }

    static final String   IGNORED_SOURCE = "internal annotation ignored";
    static final String[] COLUMNS        = { "Rank", "Node", "Taxon", "Source", "Note" };

    private final List<Row>    _rows;
    private final List<String> _summary;
    private final int          _conflicts;
    private final int          _ignored;
    private final int          _unplaced;

    private CladeAnnotationReport( final List<Row> rows, final List<String> summary, final int conflicts,
                                   final int ignored, final int unplaced ) {
        _rows = Collections.unmodifiableList( rows );
        _summary = Collections.unmodifiableList( summary );
        _conflicts = conflicts;
        _ignored = ignored;
        _unplaced = unplaced;
    }

    /**
     * Builds the report for {@code ranks} (in that order) from {@code source}, reading only the service's cache.
     * Tips come in the tree's order; the ignored internal annotations of a rank follow its tips.
     */
    static CladeAnnotationReport of( final Phylogeny phy, final List<String> ranks,
                                     final TaxonomicLineageService service, final TaxonSource source ) {
        final List<Row> rows = new ArrayList<Row>();
        final List<String> summary = new ArrayList<String>();
        int conflicts = 0;
        int ignored = 0;
        int unplaced = 0;
        if ( ( phy == null ) || phy.isEmpty() || ( ranks == null ) ) {
            return new CladeAnnotationReport( rows, summary, 0, 0, 0 );
        }
        for( final String rank : ranks ) {
            if ( ForesterUtil.isEmpty( rank ) ) {
                continue;
            }
            final Map<PhylogenyNode, TaxonPlacement> placements = TreePanelUtil.placeNodes( phy, rank, service, source );
            final Map<TaxonPlacement.Source, Integer> by_source = new LinkedHashMap<TaxonPlacement.Source, Integer>();
            int tips = 0;
            int placed = 0;
            int rank_conflicts = 0;
            for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
                final PhylogenyNode tip = it.next();
                tips++;
                final TaxonPlacement p = placements.get( tip );
                if ( p == null ) {
                    continue; // cannot happen: every tip gets a placement
                }
                by_source.merge( p.getSource(), 1, Integer::sum );
                if ( p.isPlaced() ) {
                    placed++;
                }
                String note = p.getDetail();
                if ( p.isConflict() ) {
                    rank_conflicts++;
                    note = "the tree's ancestor says \"" + p.getTreeSays().getName() + "\""
                            + ( ForesterUtil.isEmpty( note ) ? "" : ( "; " + note ) );
                }
                rows.add( new Row( rank, nodeLabel( tip ), p.isPlaced() ? p.getTaxon().getName() : null,
                                   p.getSource().getLabel(), note, false ) );
            }
            // internal annotations at this rank that tips under them contradict
            final List<String> ignored_names = new ArrayList<String>();
            for( final Entry<PhylogenyNode, TaxonPlacement> e : placements.entrySet() ) {
                final PhylogenyNode node = e.getKey();
                final TaxonPlacement p = e.getValue();
                if ( node.isExternal() || ( ( p.getSource() != TaxonPlacement.Source.OWN_ANNOTATION )
                        && ( p.getSource() != TaxonPlacement.Source.OWN_METADATA ) ) ) {
                    continue; // an internal node names a taxon by its own taxonomy or its own field
                }
                final Map<String, Integer> elsewhere = new LinkedHashMap<String, Integer>();
                for( final PhylogenyNode tip : node.getAllExternalDescendants() ) {
                    final TaxonPlacement tp = placements.get( tip );
                    if ( ( tp != null ) && tp.isPlaced() && !p.getTaxon().sameTaxon( tp.getTaxon() ) ) {
                        elsewhere.merge( tp.getTaxon().getName(), 1, Integer::sum );
                    }
                }
                if ( !elsewhere.isEmpty() ) {
                    final StringBuilder sb = new StringBuilder( "tips under it are placed as " );
                    boolean first = true;
                    for( final Entry<String, Integer> x : elsewhere.entrySet() ) {
                        sb.append( first ? "" : ", " ).append( x.getKey() ).append( " (" ).append( x.getValue() ).append( ')' );
                        first = false;
                    }
                    rows.add( new Row( rank, "[internal] " + nodeLabel( node ), p.getTaxon().getName(), IGNORED_SOURCE,
                                       sb.toString(), true ) );
                    ignored_names.add( p.getTaxon().getName() );
                }
            }
            final int rank_unplaced = tips - placed;
            conflicts += rank_conflicts;
            ignored += ignored_names.size();
            unplaced += rank_unplaced;
            summary.add( summaryLine( rank, tips, placed, by_source, rank_conflicts, ignored_names ) );
        }
        return new CladeAnnotationReport( rows, summary, conflicts, ignored, unplaced );
    }

    /** "order: 18 of 18 tips placed (ancestor's annotation 18); 0 unplaced; the database overrode the tree on 0;
     *  1 internal annotation ignored (Rodentia)". */
    static String summaryLine( final String rank, final int tips, final int placed,
                               final Map<TaxonPlacement.Source, Integer> by_source, final int conflicts,
                               final List<String> ignored_names ) {
        final StringBuilder sb = new StringBuilder( rank ).append( ": " ).append( placed ).append( " of " ).append( tips )
                .append( ( tips == 1 ) ? " tip placed" : " tips placed" );
        final StringBuilder sources = new StringBuilder();
        for( final Entry<TaxonPlacement.Source, Integer> e : by_source.entrySet() ) {
            if ( e.getKey().isPlaced() ) {
                sources.append( ( sources.length() == 0 ) ? "" : ", " ).append( e.getKey().getLabel() ).append( ' ' )
                        .append( e.getValue() );
            }
        }
        if ( sources.length() > 0 ) {
            sb.append( " (" ).append( sources ).append( ')' );
        }
        sb.append( "; " ).append( tips - placed ).append( " unplaced" );
        final StringBuilder why = new StringBuilder();
        for( final Entry<TaxonPlacement.Source, Integer> e : by_source.entrySet() ) {
            if ( !e.getKey().isPlaced() ) {
                why.append( ( why.length() == 0 ) ? "" : ", " ).append( e.getKey().getLabel().replace( "unplaced: ", "" ) )
                        .append( ' ' ).append( e.getValue() );
            }
        }
        if ( why.length() > 0 ) {
            sb.append( " (" ).append( why ).append( ')' );
        }
        sb.append( "; the database overrode the tree on " ).append( conflicts );
        sb.append( "; " ).append( ignored_names.size() ).append( ( ignored_names.size() == 1 ) ? " internal annotation ignored"
                : " internal annotations ignored" );
        if ( !ignored_names.isEmpty() ) {
            final Map<String, Integer> grouped = new LinkedHashMap<String, Integer>();
            for( final String n : ignored_names ) {
                grouped.merge( n, 1, Integer::sum );
            }
            final StringBuilder names = new StringBuilder();
            for( final Entry<String, Integer> e : grouped.entrySet() ) {
                names.append( ( names.length() == 0 ) ? "" : ", " ).append( e.getKey() )
                        .append( ( e.getValue() > 1 ) ? ( " \u00d7" + e.getValue() ) : "" );
            }
            sb.append( " (" ).append( names ).append( ')' );
        }
        return sb.toString();
    }

    /** A tip's or node's label for the report: its name, else its taxonomy label, else "(unnamed)". */
    static String nodeLabel( final PhylogenyNode node ) {
        if ( !ForesterUtil.isEmpty( node.getName() ) ) {
            return node.getName();
        }
        if ( node.getNodeData().isHasTaxonomy() ) {
            final String t = TreePanelUtil.taxonomyLabel( node.getNodeData().getTaxonomy() );
            if ( !ForesterUtil.isEmpty( t ) ) {
                return t;
            }
        }
        return "(unnamed)";
    }

    List<Row> getRows() {
        return _rows;
    }

    @Override
    public String[] getColumns() {
        return COLUMNS.clone();
    }

    @Override
    public List<String[]> getCells() {
        final List<String[]> out = new ArrayList<String[]>();
        for( final Row r : _rows ) {
            out.add( new String[] { r.getRank(), r.getNode(), r.getTaxon(), r.getSource(), r.getNote() } );
        }
        return out;
    }

    /** One line per rank. */
    List<String> getSummaryLines() {
        return _summary;
    }

    @Override
    public String getSummary() {
        return String.join( "\n", _summary );
    }

    /** Tips the database placed somewhere other than the tree's ancestor annotation, over all ranks. */
    int getConflictCount() {
        return _conflicts;
    }

    /** Internal annotations contradicted by the tips under them, over all ranks. */
    int getIgnoredInternalCount() {
        return _ignored;
    }

    /** Tip rows without a taxon, over all ranks. */
    int getUnplacedCount() {
        return _unplaced;
    }

    @Override
    public String toTsv() {
        return TableReport.tsvOf( COLUMNS, getCells() );
    }

    static String tsvCell( final String s ) {
        return TableReport.tsvCell( s );
    }

    /** The taxon of the tip-level rows at {@code rank}, keyed by node label -- a test convenience. */
    Map<String, String> taxaByNodeForTest( final String rank ) {
        final Map<String, String> out = new LinkedHashMap<String, String>();
        for( final Row r : _rows ) {
            if ( !r.isInternal() && r.getRank().equals( rank ) ) {
                out.put( r.getNode(), r.getTaxon() );
            }
        }
        return out;
    }
}
