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
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.ForesterUtil;

/**
 * Annotation import: read a tab-separated, tip-keyed table (the round-trip of {@link NodeDataExporter})
 * and write its columns onto the matching external nodes. Pure (no GUI).
 *
 * <p>The first row is the header. The <b>key column</b> is {@code node_id} when present, otherwise
 * {@code name} -- inverting the export, which only emits {@code node_id} when the tip names are not a usable
 * key. Rows are matched to tips by that key (node ids are stable only within the loaded tree, so a
 * {@code node_id} table is meant to be re-imported into the same session, not across a save/reload).
 * Unmatched rows and tips absent from the table are reported, never an error.
 *
 * <p>The table carries values only: the export does not record an identifier provider, accession source, or
 * property datatype, so those are not restored (a numeric {@code taxonomy_id} defaults to the NCBI provider;
 * a custom column becomes an {@code xsd:string} property).
 *
 * <p>Recognized columns ({@code taxonomy_scientific_name}, {@code sequence_accession}, ...) are written to
 * the corresponding model fields, creating a {@link Taxonomy}/{@link Sequence} on the tip when absent. Any
 * <b>other</b> column becomes a node {@link Property} keyed on the column header (a {@code data:} namespace
 * is prepended when the header carries no {@code ':'}); these custom columns are what column-coloring
 * consumes. A non-empty cell overwrites the existing value; an empty cell is left untouched, so a sparse
 * table only adds/updates the fields it fills. {@code branch_length} is recognized but intentionally not
 * applied -- an annotation import never changes branch lengths or tree geometry.
 */
public final class NodeDataImporter {

    private static final String TAB         = "\t";
    /** Namespace prepended to a custom column header that is not already a valid (contains {@code ':'}) property ref. */
    private static final String DATA_PREFIX = "data:";
    /** Column headers (lower-cased) that map to model fields rather than to a custom property. */
    private static final Set<String> RESERVED = Set.of( "node_id", "name", "taxonomy_scientific_name",
                                                        "taxonomy_common_name", "taxonomy_code", "taxonomy_id",
                                                        "taxonomy_rank", "sequence_name", "gene_name",
                                                        "sequence_symbol", "sequence_accession", "sequence_type",
                                                        "branch_length" );

    /**
     * Outcome of an {@link #apply} run: how many tips were annotated, which table rows matched no tip, how many
     * tree tips were absent from the table, the custom (property) columns applied (the color-by-column payoff),
     * and any per-cell warnings (values rejected by the model, e.g. an illegal taxonomy rank).
     */
    public static final class ImportResult {

        private final int          tips_annotated;
        private final int          rows_matched;
        private final List<String> unmatched_row_keys;
        private final int          tips_not_in_table;
        private final List<String> property_columns;
        private final List<String> warnings;

        private ImportResult( final int tips_annotated,
                              final int rows_matched,
                              final List<String> unmatched_row_keys,
                              final int tips_not_in_table,
                              final List<String> property_columns,
                              final List<String> warnings ) {
            this.tips_annotated = tips_annotated;
            this.rows_matched = rows_matched;
            // defensive, immutable copies so a caller cannot mutate the result's internals
            this.unmatched_row_keys = List.copyOf( unmatched_row_keys );
            this.tips_not_in_table = tips_not_in_table;
            this.property_columns = List.copyOf( property_columns );
            this.warnings = List.copyOf( warnings );
        }

        public int getTipsAnnotated() {
            return tips_annotated;
        }

        public int getRowsMatched() {
            return rows_matched;
        }

        public List<String> getUnmatchedRowKeys() {
            return unmatched_row_keys;
        }

        public int getTipsNotInTable() {
            return tips_not_in_table;
        }

        /** The custom (non-reserved) columns that were written as node properties -- usable for column coloring. */
        public List<String> getPropertyColumns() {
            return property_columns;
        }

        public List<String> getWarnings() {
            return warnings;
        }

        /** A human-readable, multi-line summary for the completion dialog. */
        public String summary() {
            final StringBuilder sb = new StringBuilder();
            sb.append( "Annotated " ).append( tips_annotated ).append( tips_annotated == 1 ? " tip" : " tips" );
            sb.append( " from " ).append( rows_matched ).append( rows_matched == 1 ? " matched row." : " matched rows." );
            if ( !property_columns.isEmpty() ) {
                sb.append( "\n\nColumns available for coloring: " ).append( String.join( ", ", property_columns ) ).append( "." );
            }
            if ( tips_not_in_table > 0 ) {
                sb.append( "\n\n" ).append( tips_not_in_table )
                  .append( tips_not_in_table == 1 ? " tip in the tree was not in the table."
                                                  : " tips in the tree were not in the table." );
            }
            if ( !unmatched_row_keys.isEmpty() ) {
                sb.append( "\n\n" ).append( unmatched_row_keys.size() )
                  .append( unmatched_row_keys.size() == 1 ? " table row matched no tip"
                                                          : " table rows matched no tip" );
                sb.append( " (" ).append( examples( unmatched_row_keys ) ).append( ")." );
            }
            if ( !warnings.isEmpty() ) {
                sb.append( "\n\n" ).append( warnings.size() )
                  .append( warnings.size() == 1 ? " value was skipped:\n" : " values were skipped:\n" );
                sb.append( examplesLines( warnings ) );
            }
            return sb.toString();
        }

        private static String examples( final List<String> items ) {
            final int show = Math.min( 5, items.size() );
            final String head = String.join( ", ", items.subList( 0, show ) );
            return ( items.size() > show ) ? head + ", and " + ( items.size() - show ) + " more" : head;
        }

        private static String examplesLines( final List<String> items ) {
            final int show = Math.min( 5, items.size() );
            final StringBuilder sb = new StringBuilder();
            for( int i = 0; i < show; i++ ) {
                sb.append( "  • " ).append( items.get( i ) ).append( '\n' );
            }
            if ( items.size() > show ) {
                sb.append( "  • ...and " ).append( items.size() - show ).append( " more" );
            }
            return sb.toString();
        }
    }

    /**
     * Apply the tab-separated table {@code tsv} to the external nodes of {@code phy}.
     *
     * @throws IllegalArgumentException if the table has no header or no usable key column ({@code name}/{@code node_id})
     */
    public static ImportResult apply( final Phylogeny phy, final String tsv ) {
        if ( phy == null ) {
            throw new IllegalArgumentException( "no tree to annotate" );
        }
        String text = ( tsv == null ) ? "" : tsv;
        if ( text.startsWith( "﻿" ) ) {
            text = text.substring( 1 ); // strip a UTF-8 byte-order mark (Excel-exported TSV often carries one)
        }
        final String[] lines = text.split( "\\R", -1 );
        if ( ( lines.length == 0 ) || ForesterUtil.isEmptyTrimmed( lines[ 0 ] ) ) {
            throw new IllegalArgumentException( "the table is empty (no header row)" );
        }
        final String[] headers = lines[ 0 ].split( TAB, -1 );
        final String[] lc_headers = new String[ headers.length ];
        for( int i = 0; i < headers.length; i++ ) {
            headers[ i ] = headers[ i ].trim();
            lc_headers[ i ] = headers[ i ].toLowerCase();
        }
        final int key_index = keyColumnIndex( lc_headers );
        if ( key_index < 0 ) {
            throw new IllegalArgumentException( "the table has no \"name\" or \"node_id\" key column" );
        }
        final boolean key_is_node_id = "node_id".equals( lc_headers[ key_index ] );
        final Map<String, List<PhylogenyNode>> by_key = indexTips( phy, key_is_node_id );

        final Set<Long> annotated_ids = new LinkedHashSet<>();
        final Set<String> matched_keys = new LinkedHashSet<>();
        final List<String> unmatched_row_keys = new ArrayList<>();
        final LinkedHashSet<String> property_columns = new LinkedHashSet<>();
        final List<String> warnings = new ArrayList<>();
        int rows_matched = 0;

        for( int r = 1; r < lines.length; r++ ) {
            if ( ForesterUtil.isEmpty( lines[ r ] ) ) {
                continue; // skip blank lines (e.g. a trailing newline)
            }
            final String[] cells = lines[ r ].split( TAB, -1 );
            final String key = ( key_index < cells.length ) ? cells[ key_index ].trim() : "";
            final List<PhylogenyNode> targets = key.isEmpty() ? null : by_key.get( key );
            if ( ( targets == null ) || targets.isEmpty() ) {
                unmatched_row_keys.add( key.isEmpty() ? "(blank)" : key );
                continue;
            }
            rows_matched++;
            matched_keys.add( key );
            for( int c = 0; c < headers.length; c++ ) {
                if ( ( c == key_index ) || ( c >= cells.length ) ) {
                    continue;
                }
                final String value = cells[ c ].trim();
                if ( value.isEmpty() ) {
                    continue; // never clobber an existing value with a blank cell
                }
                final boolean reserved = RESERVED.contains( lc_headers[ c ] );
                if ( reserved && "branch_length".equals( lc_headers[ c ] ) ) {
                    continue; // recognized but intentionally not applied (annotation import leaves geometry alone)
                }
                final String prop_ref = reserved ? null : propertyRef( headers[ c ] );
                for( final PhylogenyNode n : targets ) {
                    try {
                        if ( reserved ) {
                            applyReservedField( n, lc_headers[ c ], value );
                        }
                        else {
                            setProperty( n, prop_ref, value );
                        }
                        annotated_ids.add( n.getId() );
                    }
                    catch ( final PhyloXmlDataFormatException e ) {
                        warnings.add( "row \"" + key + "\", column \"" + headers[ c ] + "\": " + e.getMessage() );
                    }
                }
                if ( !reserved ) {
                    property_columns.add( prop_ref );
                }
            }
        }

        int tips_not_in_table = 0;
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            final String key = tipKey( n, key_is_node_id );
            if ( key.isEmpty() || !matched_keys.contains( key ) ) {
                tips_not_in_table++;
            }
        }
        return new ImportResult( annotated_ids.size(), rows_matched, unmatched_row_keys, tips_not_in_table,
                                 new ArrayList<>( property_columns ), warnings );
    }

    private static int keyColumnIndex( final String[] lc_headers ) {
        int name_index = -1;
        for( int i = 0; i < lc_headers.length; i++ ) {
            if ( "node_id".equals( lc_headers[ i ] ) ) {
                return i; // node_id wins when present (it is the export's unique key for non-unique names)
            }
            if ( ( name_index < 0 ) && "name".equals( lc_headers[ i ] ) ) {
                name_index = i;
            }
        }
        return name_index;
    }

    /** Map each tip's key (its name, or its node id as text) to the tip(s) that carry it; blank keys are skipped. */
    private static Map<String, List<PhylogenyNode>> indexTips( final Phylogeny phy, final boolean key_is_node_id ) {
        final Map<String, List<PhylogenyNode>> by_key = new TreeMap<>();
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            final String key = tipKey( n, key_is_node_id );
            if ( !key.isEmpty() ) {
                by_key.computeIfAbsent( key, k -> new ArrayList<>() ).add( n );
            }
        }
        return by_key;
    }

    /**
     * A tip's match key. The row key cells are trimmed (see {@link #apply}), so the tip name is trimmed too --
     * otherwise a tip named {@code " A "} would never match a {@code "A"} row. Node ids carry no whitespace.
     */
    private static String tipKey( final PhylogenyNode n, final boolean key_is_node_id ) {
        if ( key_is_node_id ) {
            return String.valueOf( n.getId() );
        }
        return ForesterUtil.isEmpty( n.getName() ) ? "" : n.getName().trim();
    }

    private static void applyReservedField( final PhylogenyNode n, final String col, final String value )
            throws PhyloXmlDataFormatException {
        switch ( col ) {
            case "name":
                n.setName( value );
                break;
            case "taxonomy_scientific_name":
                writeTaxonomy( n, t -> t.setScientificName( value ) );
                break;
            case "taxonomy_common_name":
                writeTaxonomy( n, t -> t.setCommonName( value ) );
                break;
            case "taxonomy_code":
                writeTaxonomy( n, t -> t.setTaxonomyCode( value ) );
                break;
            case "taxonomy_id":
                writeTaxonomy( n, t -> t.setIdentifier( isAllDigits( value ) ? new Identifier( value, "ncbi" )
                                                                             : new Identifier( value ) ) );
                break;
            case "taxonomy_rank":
                writeTaxonomy( n, t -> t.setRank( value.toLowerCase() ) );
                break;
            case "sequence_name":
                writeSequence( n, s -> s.setName( value ) );
                break;
            case "gene_name":
                writeSequence( n, s -> s.setGeneName( value ) );
                break;
            case "sequence_symbol":
                writeSequence( n, s -> s.setSymbol( value ) );
                break;
            case "sequence_accession":
                writeSequence( n, s -> s.setAccession( new Accession( value ) ) );
                break;
            case "sequence_type":
                writeSequence( n, s -> s.setType( value ) );
                break;
            default:
                break; // node_id / branch_length are never applied here
        }
    }

    @FunctionalInterface
    private interface TaxonomyWrite {
        void apply( Taxonomy t ) throws PhyloXmlDataFormatException;
    }

    @FunctionalInterface
    private interface SequenceWrite {
        void apply( Sequence s ) throws PhyloXmlDataFormatException;
    }

    /**
     * Apply a taxonomy field, attaching a freshly-created {@link Taxonomy} only when the write succeeds. The
     * validating setters ({@code setRank}/{@code setTaxonomyCode}) throw before mutating, so a rejected cell
     * never leaves a phantom empty taxonomy behind (which would flip {@code isHasTaxonomy()} on with no data).
     */
    private static void writeTaxonomy( final PhylogenyNode n, final TaxonomyWrite write )
            throws PhyloXmlDataFormatException {
        final boolean had = n.getNodeData().isHasTaxonomy();
        final Taxonomy t = had ? n.getNodeData().getTaxonomy() : new Taxonomy();
        write.apply( t ); // throws on an invalid value before we attach a newly-created taxonomy
        if ( !had ) {
            n.getNodeData().addTaxonomy( t );
        }
    }

    /** As {@link #writeTaxonomy} but for the tip's {@link Sequence} ({@code setSymbol}/{@code setType} validate). */
    private static void writeSequence( final PhylogenyNode n, final SequenceWrite write )
            throws PhyloXmlDataFormatException {
        final boolean had = n.getNodeData().isHasSequence();
        final Sequence s = had ? n.getNodeData().getSequence() : new Sequence();
        write.apply( s );
        if ( !had ) {
            n.getNodeData().addSequence( s );
        }
    }

    /** ASCII-digits check (matching the prior {@code \d+} test) used to default a numeric taxonomy id to NCBI. */
    private static boolean isAllDigits( final String s ) {
        if ( s.isEmpty() ) {
            return false;
        }
        for( int i = 0; i < s.length(); i++ ) {
            final char c = s.charAt( i );
            if ( ( c < '0' ) || ( c > '9' ) ) {
                return false;
            }
        }
        return true;
    }

    /** Set (overwrite) the custom property {@code ref} on the tip, replacing any existing property with that ref. */
    private static void setProperty( final PhylogenyNode n, final String ref, final String value ) {
        PropertiesList props = n.getNodeData().getProperties();
        if ( props == null ) {
            props = new PropertiesList();
            n.getNodeData().setProperties( props );
        }
        final Iterator<Property> it = props.getProperties().iterator();
        while ( it.hasNext() ) {
            if ( ref.equals( it.next().getRef() ) ) {
                it.remove();
            }
        }
        props.addProperty( new Property( ref, value, "", "xsd:string", AppliesTo.NODE ) );
    }

    /** A legal property ref for a custom column: the header verbatim if it already namespaces with {@code ':'}, else prefixed. */
    private static String propertyRef( final String header ) {
        return ( header.indexOf( ':' ) >= 1 ) ? header : DATA_PREFIX + header;
    }

    private NodeDataImporter() {
    }
}
