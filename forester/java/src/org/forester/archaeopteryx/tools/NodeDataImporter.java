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
 * Annotation import: read a tip-keyed <b>CSV or TSV</b> table and write its columns onto the matching external
 * nodes (the "Import Annotations" tool). Pure (no GUI).
 *
 * <p>The first row is the header. Rows are joined to tips by a <b>keyed join</b>: a chosen key column matched
 * against a chosen tip attribute (see {@link MatchBy} -- tip name, sequence accession, taxonomy id, or taxonomy
 * scientific name). {@link #parseTable(String)} auto-detects the delimiter (tab, else comma) and unquotes CSV
 * fields; {@link #dryRun} reports the match counts before committing; {@link #apply(Phylogeny, Table, int, MatchBy)}
 * writes the columns. Unmatched rows and tips absent from the table are reported, never an error.
 *
 * <p>Recognized columns ({@code taxonomy_scientific_name}, {@code sequence_accession}, ...) are written to the
 * corresponding model fields, creating a {@link Taxonomy}/{@link Sequence} on the tip when absent. Any <b>other</b>
 * column becomes a node {@link Property} keyed on the column header (a {@code data:} namespace is prepended when
 * the header carries no {@code ':'}); these custom columns are what column-coloring consumes. A non-empty cell
 * overwrites the existing value; an empty cell is left untouched, so a sparse table only adds/updates the fields
 * it fills. {@code branch_length} is recognized but intentionally not applied -- an annotation import never
 * changes branch lengths or tree geometry.
 *
 * <p>The legacy {@link #apply(Phylogeny, String)} TSV entry point is kept: it parses the table, picks the
 * {@code node_id}/{@code name} key column, and delegates -- the round-trip of {@link NodeDataExporter}.
 */
public final class NodeDataImporter {

    /** Namespace prepended to a custom column header that is not already a valid (contains {@code ':'}) property ref. */
    private static final String DATA_PREFIX = "data:";
    /** Column headers (lower-cased) that map to model fields rather than to a custom property. */
    private static final Set<String> RESERVED = Set.of( "node_id", "name", "taxonomy_scientific_name",
                                                        "taxonomy_common_name", "taxonomy_code", "taxonomy_id",
                                                        "taxonomy_rank", "sequence_name", "gene_name",
                                                        "sequence_symbol", "sequence_accession", "sequence_type",
                                                        "branch_length" );

    /** The tip attribute a table's key column is matched against for the keyed join. */
    public enum MatchBy {
        /** The tip's node id as text -- internal, stable only within the loaded session (the {@link NodeDataExporter} round-trip). */
        NODE_ID( "Node id" ),
        /** The tip's name ({@link PhylogenyNode#getName()}) -- the default, and the common case (names = labels/accessions). */
        TIP_NAME( "Tip name" ),
        /** The tip's sequence accession value. */
        SEQUENCE_ACCESSION( "Sequence accession" ),
        /** The tip's taxonomy identifier value (e.g. an NCBI tax id). */
        TAXONOMY_ID( "Taxonomy id" ),
        /** The tip's taxonomy scientific name. */
        TAXONOMY_SCIENTIFIC_NAME( "Taxonomy scientific name" );

        private final String _label;

        MatchBy( final String label ) {
            _label = label;
        }

        @Override
        public String toString() {
            return _label;
        }

        /** This tip's key for this match attribute, trimmed; {@code ""} when the tip lacks it (never matched). */
        String nodeKey( final PhylogenyNode n ) {
            switch ( this ) {
                case NODE_ID:
                    return String.valueOf( n.getId() );
                case TIP_NAME:
                    return trimmed( n.getName() );
                case SEQUENCE_ACCESSION:
                    if ( n.getNodeData().isHasSequence() && ( n.getNodeData().getSequence().getAccession() != null ) ) {
                        return trimmed( n.getNodeData().getSequence().getAccession().getValue() );
                    }
                    return "";
                case TAXONOMY_ID:
                    if ( n.getNodeData().isHasTaxonomy() && ( n.getNodeData().getTaxonomy().getIdentifier() != null ) ) {
                        return trimmed( n.getNodeData().getTaxonomy().getIdentifier().getValue() );
                    }
                    return "";
                case TAXONOMY_SCIENTIFIC_NAME:
                    if ( n.getNodeData().isHasTaxonomy() ) {
                        return trimmed( n.getNodeData().getTaxonomy().getScientificName() );
                    }
                    return "";
                default:
                    return "";
            }
        }
    }

    /** The user-selectable match attributes (NODE_ID is internal, used only by the legacy TSV entry point). */
    public static MatchBy[] userMatchOptions() {
        return new MatchBy[] { MatchBy.TIP_NAME, MatchBy.SEQUENCE_ACCESSION, MatchBy.TAXONOMY_ID,
                MatchBy.TAXONOMY_SCIENTIFIC_NAME };
    }

    /** A parsed delimited table: trimmed headers + the raw data rows (each cell trimmed at use, as the TSV path did). */
    public static final class Table {

        private final String[]       _headers;    // trimmed
        private final String[]       _lc_headers; // lower-cased (for reserved-column / key lookup)
        private final List<String[]> _rows;       // data rows; a row may be shorter than the header (a short/ragged line)
        private final char           _delimiter;

        private Table( final String[] headers, final String[] lc_headers, final List<String[]> rows, final char delimiter ) {
            _headers = headers;
            _lc_headers = lc_headers;
            _rows = rows;
            _delimiter = delimiter;
        }

        public String[] getHeaders() {
            return _headers.clone();
        }

        public int getColumnCount() {
            return _headers.length;
        }

        public int getRowCount() {
            return _rows.size();
        }

        /** Cell {@code (row, col)}, trimmed, or {@code ""} when the row is too short (a ragged line). */
        public String getCell( final int row, final int col ) {
            final String[] cells = _rows.get( row );
            return ( col < cells.length ) ? cells[ col ].trim() : "";
        }

        public char getDelimiter() {
            return _delimiter;
        }

        /**
         * The best default key column for the dialog: a {@code name} column, else a {@code node_id} column, else 0.
         * The dialog matches by a tip ATTRIBUTE (name/accession/taxonomy) and cannot match by node_id, so {@code name}
         * is preferred over {@code node_id} here (the legacy TSV round-trip uses {@link #keyColumnIndex} directly).
         */
        public int defaultKeyColumn() {
            int node_id = -1;
            for( int i = 0; i < _lc_headers.length; i++ ) {
                if ( "name".equals( _lc_headers[ i ] ) ) {
                    return i;
                }
                if ( ( node_id < 0 ) && "node_id".equals( _lc_headers[ i ] ) ) {
                    node_id = i;
                }
            }
            return ( node_id >= 0 ) ? node_id : 0;
        }
    }

    /**
     * Outcome of a {@link #dryRun}: the row-join counts and the columns that would be imported -- for a preview
     * shown to the user before committing.
     */
    public static final class MatchReport {

        private final int          _total_rows;
        private final int          _rows_matched;   // rows whose key hit >= 1 tip
        private final int          _rows_ambiguous; // rows whose key hit > 1 tip (a subset of matched)
        private final int          _rows_unmatched; // rows whose key hit 0 tips (or a blank key)
        private final int          _total_tips;
        private final int          _tips_without_row;
        private final List<String> _property_columns; // become data:* properties (color-able)
        private final List<String> _reserved_columns; // fill model fields (taxonomy_*, sequence_*, name)

        private MatchReport( final int total_rows, final int rows_matched, final int rows_ambiguous,
                             final int rows_unmatched, final int total_tips, final int tips_without_row,
                             final List<String> property_columns, final List<String> reserved_columns ) {
            _total_rows = total_rows;
            _rows_matched = rows_matched;
            _rows_ambiguous = rows_ambiguous;
            _rows_unmatched = rows_unmatched;
            _total_tips = total_tips;
            _tips_without_row = tips_without_row;
            _property_columns = List.copyOf( property_columns );
            _reserved_columns = List.copyOf( reserved_columns );
        }

        public int getTotalRows() { return _total_rows; }
        public int getRowsMatched() { return _rows_matched; }
        public int getRowsAmbiguous() { return _rows_ambiguous; }
        public int getRowsUnmatched() { return _rows_unmatched; }
        public int getTotalTips() { return _total_tips; }
        public int getTipsWithoutRow() { return _tips_without_row; }
        public List<String> getPropertyColumns() { return _property_columns; }
        public List<String> getReservedColumns() { return _reserved_columns; }

        /** A one-line preview for the dialog, e.g. "42/45 rows match · 3 unmatched · 2 of 44 tips have no row · imports: host, reads". */
        public String summaryLine() {
            final StringBuilder sb = new StringBuilder();
            sb.append( _rows_matched ).append( '/' ).append( _total_rows ).append( _total_rows == 1 ? " row" : " rows" )
              .append( " match" );
            if ( _rows_ambiguous > 0 ) {
                sb.append( " (" ).append( _rows_ambiguous ).append( " ambiguous)" );
            }
            if ( _rows_unmatched > 0 ) {
                sb.append( " · " ).append( _rows_unmatched ).append( " unmatched" );
            }
            if ( _tips_without_row > 0 ) {
                sb.append( " · " ).append( _tips_without_row ).append( " of " ).append( _total_tips )
                  .append( " tips have no row" );
            }
            if ( !_property_columns.isEmpty() ) {
                sb.append( " · imports: " ).append( String.join( ", ", _property_columns ) );
            }
            return sb.toString();
        }
    }

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

    /** Parse CSV/TSV text into a {@link Table} with the delimiter auto-detected. */
    public static Table parseTable( final String raw ) {
        return parseTable( raw, null );
    }

    /**
     * Parse CSV/TSV text into a {@link Table}: strip a UTF-8 BOM, then tokenize into records. The delimiter is
     * {@code forced_delimiter} when given, else auto-detected from the first line (TAB if it carries a tab, else
     * comma, else TAB). A comma table is parsed RFC-4180-style: double-quoted fields may contain the delimiter,
     * doubled {@code ""} quotes, AND embedded newlines (a value may span lines). A tab table is split literally
     * (quotes are ordinary characters, as TSV values rarely quote).
     *
     * @throws IllegalArgumentException if the text has no non-blank header row
     */
    public static Table parseTable( final String raw, final Character forced_delimiter ) {
        String text = ( raw == null ) ? "" : raw;
        if ( text.startsWith( "﻿" ) ) {
            text = text.substring( 1 ); // strip a UTF-8 byte-order mark (Excel-exported files often carry one)
        }
        final char delim = ( forced_delimiter != null ) ? forced_delimiter.charValue() : detectDelimiter( text );
        final List<String[]> records = ( delim == ',' ) ? tokenizeCsv( text ) : tokenizeTsv( text );
        // the header is the FIRST non-blank record (so a leading blank line does not break the parse); the data rows
        // follow it, blank records skipped
        int header_idx = 0;
        while ( ( header_idx < records.size() ) && isBlankRecord( records.get( header_idx ) ) ) {
            header_idx++;
        }
        if ( header_idx >= records.size() ) {
            throw new IllegalArgumentException( "the table is empty (no header row)" );
        }
        final String[] headers = records.get( header_idx );
        final String[] lc_headers = new String[ headers.length ];
        for( int i = 0; i < headers.length; i++ ) {
            headers[ i ] = headers[ i ].trim();
            lc_headers[ i ] = headers[ i ].toLowerCase();
        }
        final List<String[]> rows = new ArrayList<>();
        for( int r = header_idx + 1; r < records.size(); r++ ) {
            if ( !isBlankRecord( records.get( r ) ) ) {
                rows.add( records.get( r ) ); // skip fully-blank records (e.g. a trailing newline)
            }
        }
        return new Table( headers, lc_headers, rows, delim );
    }

    /** Auto-detect the delimiter from the first NON-BLANK line: TAB if it carries a tab, else comma, else TAB. */
    private static char detectDelimiter( final String text ) {
        for( final String line : text.split( "\\R", -1 ) ) {
            if ( !ForesterUtil.isEmpty( line ) ) {
                if ( line.indexOf( '\t' ) >= 0 ) {
                    return '\t';
                }
                return ( line.indexOf( ',' ) >= 0 ) ? ',' : '\t';
            }
        }
        return '\t';
    }

    private static boolean isBlankRecord( final String[] record ) {
        for( final String c : record ) {
            if ( !ForesterUtil.isEmpty( c ) ) {
                return false;
            }
        }
        return true;
    }

    /** Tab-delimited: one record per line (any Unicode break), each split literally on tabs (no quote handling). */
    private static List<String[]> tokenizeTsv( final String text ) {
        final List<String[]> records = new ArrayList<>();
        for( final String line : text.split( "\\R", -1 ) ) {
            records.add( line.split( "\t", -1 ) );
        }
        return records;
    }

    /**
     * RFC-4180 comma tokenizer over the WHOLE text: a double-quoted field may contain commas, doubled {@code ""}
     * quotes, and newlines (so a value can span lines); an unquoted CR/LF/CRLF ends the record.
     */
    private static List<String[]> tokenizeCsv( final String text ) {
        final List<String[]> records = new ArrayList<>();
        List<String> record = new ArrayList<>();
        final StringBuilder cur = new StringBuilder();
        boolean in_quotes = false;
        int i = 0;
        final int n = text.length();
        while ( i < n ) {
            final char c = text.charAt( i );
            if ( in_quotes ) {
                if ( c == '"' ) {
                    if ( ( ( i + 1 ) < n ) && ( text.charAt( i + 1 ) == '"' ) ) {
                        cur.append( '"' ); // a doubled "" is a literal quote
                        i += 2;
                    }
                    else {
                        in_quotes = false;
                        i++;
                    }
                }
                else {
                    cur.append( c ); // commas and newlines are literal inside quotes
                    i++;
                }
            }
            else if ( c == '"' ) {
                in_quotes = true;
                i++;
            }
            else if ( c == ',' ) {
                record.add( cur.toString() );
                cur.setLength( 0 );
                i++;
            }
            else if ( ( c == '\r' ) || ( c == '\n' ) ) {
                record.add( cur.toString() );
                cur.setLength( 0 );
                records.add( record.toArray( new String[ 0 ] ) );
                record = new ArrayList<>();
                i++;
                if ( ( c == '\r' ) && ( i < n ) && ( text.charAt( i ) == '\n' ) ) {
                    i++; // consume the LF of a CRLF pair
                }
            }
            else {
                cur.append( c );
                i++;
            }
        }
        if ( ( cur.length() > 0 ) || !record.isEmpty() ) {
            record.add( cur.toString() ); // flush the last field/record when there is no trailing newline
            records.add( record.toArray( new String[ 0 ] ) );
        }
        return records;
    }

    /** A read-only preview of the join importing every non-key column (see {@link #dryRun(Phylogeny, Table, int, MatchBy, ColumnPlan)}). */
    public static MatchReport dryRun( final Phylogeny phy, final Table table, final int key_col, final MatchBy match_by ) {
        return dryRun( phy, table, key_col, match_by, ColumnPlan.importAll( table ) );
    }

    /**
     * A read-only preview of the keyed join: how many rows match / are ambiguous / are unmatched, how many tips
     * have no row, and the columns that would be imported (honoring the {@code plan}'s include/rename choices).
     * Does not mutate the tree.
     */
    public static MatchReport dryRun( final Phylogeny phy, final Table table, final int key_col, final MatchBy match_by,
                                      final ColumnPlan plan ) {
        if ( phy == null ) {
            throw new IllegalArgumentException( "no tree to annotate" );
        }
        if ( ( key_col < 0 ) || ( key_col >= table.getColumnCount() ) ) {
            throw new IllegalArgumentException( "the key column is out of range" );
        }
        final Map<String, List<PhylogenyNode>> by_key = indexTips( phy, match_by );
        int matched = 0, ambiguous = 0, unmatched = 0;
        final Set<String> matched_keys = new LinkedHashSet<>();
        for( int r = 0; r < table.getRowCount(); r++ ) {
            final String key = table.getCell( r, key_col );
            final List<PhylogenyNode> targets = key.isEmpty() ? null : by_key.get( key );
            if ( ( targets == null ) || targets.isEmpty() ) {
                unmatched++;
            }
            else {
                matched++;
                if ( targets.size() > 1 ) {
                    ambiguous++;
                }
                matched_keys.add( key );
            }
        }
        int total_tips = 0, tips_without_row = 0;
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            total_tips++;
            final String key = match_by.nodeKey( n );
            if ( key.isEmpty() || !matched_keys.contains( key ) ) {
                tips_without_row++;
            }
        }
        final List<String> property_cols = new ArrayList<>();
        final List<String> reserved_cols = new ArrayList<>();
        for( int c = 0; c < table.getColumnCount(); c++ ) {
            if ( ( c == key_col ) || !plan.isIncluded( c ) ) {
                continue;
            }
            final String header = plan.header( c ).trim();
            if ( header.isEmpty() || isNeverApplied( header.toLowerCase() ) ) {
                continue; // a blanked rename has no ref; node_id / branch_length are never written as annotations
            }
            if ( RESERVED.contains( header.toLowerCase() ) ) {
                reserved_cols.add( header ); // recognized model field
            }
            else {
                property_cols.add( header );
            }
        }
        return new MatchReport( table.getRowCount(), matched, ambiguous, unmatched, total_tips, tips_without_row,
                                property_cols, reserved_cols );
    }

    /**
     * Which columns to import and under what (possibly renamed) header. A rename changes how the column is
     * classified: renaming a column to a reserved name (e.g. {@code taxonomy_id}) fills that model field; any other
     * name becomes a {@code data:} property. The key column is never imported regardless of its flag.
     */
    public static final class ColumnPlan {

        private final boolean[] _included;
        private final String[]  _header; // effective (possibly renamed) header per column

        private ColumnPlan( final boolean[] included, final String[] header ) {
            _included = included;
            _header = header;
        }

        /** The default plan for a table: import every column under its own header. */
        public static ColumnPlan importAll( final Table table ) {
            final int n = table.getColumnCount();
            final boolean[] inc = new boolean[ n ];
            java.util.Arrays.fill( inc, true );
            return new ColumnPlan( inc, table.getHeaders() );
        }

        public void setIncluded( final int col, final boolean included ) {
            _included[ col ] = included;
        }

        public void setHeader( final int col, final String header ) {
            _header[ col ] = ( header == null ) ? "" : header;
        }

        public boolean isIncluded( final int col ) {
            return _included[ col ];
        }

        public String header( final int col ) {
            return _header[ col ];
        }
    }

    /**
     * A remembered import, so it can be repeated with one click after the source (file or URL) changes -- the
     * annotation "profile". Keyed by column HEADER NAME (not index) + the delimiter used, so a re-import survives
     * the source gaining/reordering rows or columns: a new column is default-included, a missing key column falls
     * back to the table's default.
     */
    public static final class ImportProfile {

        private final String              _source;           // the file path or URL to re-fetch
        private final boolean             _is_url;
        private final Character           _delimiter;        // the delimiter that was used (reproduce the same parse)
        private final MatchBy             _match_by;
        private final String              _key_header;       // the key column's header name
        private final Set<String>         _excluded_headers; // original headers the user unchecked
        private final Map<String, String> _renames;          // original header -> effective (renamed) header

        private ImportProfile( final String source, final boolean is_url, final Character delimiter,
                final MatchBy match_by, final String key_header, final Set<String> excluded_headers,
                final Map<String, String> renames ) {
            _source = source;
            _is_url = is_url;
            _delimiter = delimiter;
            _match_by = match_by;
            _key_header = key_header;
            _excluded_headers = excluded_headers;
            _renames = renames;
        }

        /** Capture the profile from a completed import choice. */
        public static ImportProfile from( final Table table, final int key_col, final MatchBy match_by,
                final ColumnPlan plan, final String source, final boolean is_url ) {
            final String[] headers = table.getHeaders();
            final Set<String> excluded = new LinkedHashSet<>();
            final Map<String, String> renames = new java.util.LinkedHashMap<>();
            for( int c = 0; c < headers.length; c++ ) {
                if ( c == key_col ) {
                    continue;
                }
                if ( !plan.isIncluded( c ) ) {
                    excluded.add( headers[ c ] );
                }
                else if ( !headers[ c ].equals( plan.header( c ) ) ) {
                    renames.put( headers[ c ], plan.header( c ) );
                }
            }
            return new ImportProfile( source, is_url, Character.valueOf( table.getDelimiter() ), match_by,
                                      headers[ key_col ], excluded, renames );
        }

        public String getSource() {
            return _source;
        }

        public boolean isUrl() {
            return _is_url;
        }

        public Character getDelimiter() {
            return _delimiter;
        }

        public MatchBy getMatchBy() {
            return _match_by;
        }

        /** The key column index in {@code table} by the remembered header name, or the table's default if absent. */
        public int keyColumn( final Table table ) {
            final String[] headers = table.getHeaders();
            for( int c = 0; c < headers.length; c++ ) {
                if ( _key_header.equals( headers[ c ] ) ) {
                    return c;
                }
            }
            return table.defaultKeyColumn();
        }

        /** Rebuild the ColumnPlan against a freshly-parsed {@code table}: exclude the remembered headers, apply the
         *  remembered renames, and default-include any NEW column. */
        public ColumnPlan columnPlan( final Table table ) {
            final ColumnPlan plan = ColumnPlan.importAll( table );
            final int key_col = keyColumn( table );
            final String[] headers = table.getHeaders();
            for( int c = 0; c < headers.length; c++ ) {
                plan.setIncluded( c, ( c != key_col ) && !_excluded_headers.contains( headers[ c ] ) );
                plan.setHeader( c, _renames.getOrDefault( headers[ c ], headers[ c ] ) );
            }
            return plan;
        }
    }

    /** Apply the parsed {@code table} importing every non-key column (see {@link #apply(Phylogeny, Table, int, MatchBy, ColumnPlan)}). */
    public static ImportResult apply( final Phylogeny phy, final Table table, final int key_col, final MatchBy match_by ) {
        return apply( phy, table, key_col, match_by, ColumnPlan.importAll( table ) );
    }

    /**
     * Apply the parsed {@code table} to {@code phy}, matching its {@code key_col} against each tip's
     * {@code match_by} attribute and importing the columns the {@code plan} includes (under their possibly-renamed
     * headers). Non-empty cells overwrite; blank cells never clobber; {@code branch_length} is ignored. Returns the
     * accounting; a value the model rejects becomes a warning, never a thrown error.
     */
    public static ImportResult apply( final Phylogeny phy, final Table table, final int key_col, final MatchBy match_by,
                                      final ColumnPlan plan ) {
        if ( phy == null ) {
            throw new IllegalArgumentException( "no tree to annotate" );
        }
        if ( ( key_col < 0 ) || ( key_col >= table.getColumnCount() ) ) {
            throw new IllegalArgumentException( "the key column is out of range" );
        }
        final Map<String, List<PhylogenyNode>> by_key = indexTips( phy, match_by );

        final Set<Long> annotated_ids = new LinkedHashSet<>();
        final Set<String> matched_keys = new LinkedHashSet<>();
        final List<String> unmatched_row_keys = new ArrayList<>();
        final LinkedHashSet<String> property_columns = new LinkedHashSet<>();
        final List<String> warnings = new ArrayList<>();
        int rows_matched = 0;

        for( int r = 0; r < table.getRowCount(); r++ ) {
            final String key = table.getCell( r, key_col );
            final List<PhylogenyNode> targets = key.isEmpty() ? null : by_key.get( key );
            if ( ( targets == null ) || targets.isEmpty() ) {
                unmatched_row_keys.add( key.isEmpty() ? "(blank)" : key );
                continue;
            }
            rows_matched++;
            matched_keys.add( key );
            for( int c = 0; c < table.getColumnCount(); c++ ) {
                if ( ( c == key_col ) || !plan.isIncluded( c ) ) {
                    continue; // the key column and de-selected columns are never imported
                }
                final String value = table.getCell( r, c );
                if ( value.isEmpty() ) {
                    continue; // never clobber an existing value with a blank cell
                }
                final String header = plan.header( c ).trim();
                final String lc = header.toLowerCase();
                if ( header.isEmpty() || isNeverApplied( lc ) ) {
                    // a blanked rename has no ref (avoid a malformed "data:" property); node_id / branch_length are
                    // identity/geometry columns never written as annotations (and would otherwise falsely count the
                    // tip as annotated -> a spurious undo checkpoint + provenance for a no-op import)
                    continue;
                }
                final boolean reserved = RESERVED.contains( lc );
                final String prop_ref = reserved ? null : propertyRef( header );
                for( final PhylogenyNode n : targets ) {
                    try {
                        if ( reserved ) {
                            applyReservedField( n, lc, value );
                        }
                        else {
                            setProperty( n, prop_ref, value );
                        }
                        annotated_ids.add( n.getId() );
                    }
                    catch ( final PhyloXmlDataFormatException e ) {
                        warnings.add( "row \"" + key + "\", column \"" + header + "\": " + e.getMessage() );
                    }
                }
                if ( !reserved ) {
                    property_columns.add( prop_ref );
                }
            }
        }

        int tips_not_in_table = 0;
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            final String key = match_by.nodeKey( n );
            if ( key.isEmpty() || !matched_keys.contains( key ) ) {
                tips_not_in_table++;
            }
        }
        return new ImportResult( annotated_ids.size(), rows_matched, unmatched_row_keys, tips_not_in_table,
                                 new ArrayList<>( property_columns ), warnings );
    }

    /**
     * Legacy TSV entry point (the {@link NodeDataExporter} round-trip): parse the table, pick the
     * {@code node_id}/{@code name} key column, and apply.
     *
     * @throws IllegalArgumentException if the table has no header or no usable key column ({@code name}/{@code node_id})
     */
    public static ImportResult apply( final Phylogeny phy, final String tsv ) {
        if ( phy == null ) {
            throw new IllegalArgumentException( "no tree to annotate" );
        }
        final Table table = parseTable( tsv );
        final String[] lc_headers = new String[ table.getColumnCount() ];
        final String[] headers = table.getHeaders();
        for( int i = 0; i < headers.length; i++ ) {
            lc_headers[ i ] = headers[ i ].toLowerCase();
        }
        final int key_index = keyColumnIndex( lc_headers );
        if ( key_index < 0 ) {
            throw new IllegalArgumentException( "the table has no \"name\" or \"node_id\" key column" );
        }
        final MatchBy match_by = "node_id".equals( lc_headers[ key_index ] ) ? MatchBy.NODE_ID : MatchBy.TIP_NAME;
        return apply( phy, table, key_index, match_by );
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

    /** Map each tip's {@code match_by} key to the tip(s) that carry it; blank keys are skipped (a list, so a shared key annotates all). */
    private static Map<String, List<PhylogenyNode>> indexTips( final Phylogeny phy, final MatchBy match_by ) {
        final Map<String, List<PhylogenyNode>> by_key = new TreeMap<>();
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            final String key = match_by.nodeKey( n );
            if ( !key.isEmpty() ) {
                by_key.computeIfAbsent( key, k -> new ArrayList<>() ).add( n );
            }
        }
        return by_key;
    }

    private static String trimmed( final String s ) {
        return ForesterUtil.isEmpty( s ) ? "" : s.trim();
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

    /** Reserved columns that are recognized but NEVER written by an import: {@code node_id} (identity key, no model
     *  field) and {@code branch_length} (geometry). Skipping them also keeps a no-op column from falsely counting a
     *  tip as annotated. */
    private static boolean isNeverApplied( final String lc_header ) {
        return "node_id".equals( lc_header ) || "branch_length".equals( lc_header );
    }

    private NodeDataImporter() {
    }
}
