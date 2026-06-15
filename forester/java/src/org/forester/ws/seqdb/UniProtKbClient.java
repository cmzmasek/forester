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

package org.forester.ws.seqdb;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.forester.util.ForesterUtil;

/**
 * A clean UniProtKB client over the UniProt REST API, returning a {@link SequenceEntry} for an
 * accession. It requests a single <b>TSV</b> row with an explicit, fixed field list (no JSON library
 * needed) and parses it positionally. Network I/O goes through {@link WsHttp} (timeout-bounded +
 * throttled). Results -- including a definitive "not found" -- are cached for the life of the process.
 */
public final class UniProtKbClient {

    // Requested fields, in a fixed order; the TSV response has one column per field in this order.
    private static final String[]            FIELDS = { "accession", "id", "protein_name", "gene_primary",
            "gene_names", "sequence", "organism_name", "organism_id", "lineage", "go_id", "xref_refseq", "xref_pdb" };
    private static final String              BASE   = "https://rest.uniprot.org/uniprotkb/search?format=tsv&size=1&fields="
            + String.join( ",", FIELDS ) + "&query=accession:";
    private final Map<String, SequenceEntry> _cache = Collections
            .synchronizedMap( new HashMap<String, SequenceEntry>() );

    /**
     * The UniProtKB entry for {@code accession}, or {@link SequenceEntry#EMPTY} when not found. Caches
     * the result (a not-found accession is a definitive negative). Does network I/O -- call off the EDT.
     *
     * @throws IOException on a connection/transport failure.
     */
    public SequenceEntry fetchEntry( final String accession ) throws IOException {
        if ( ForesterUtil.isEmpty( accession ) ) {
            return SequenceEntry.EMPTY;
        }
        final String k = accession.trim().toUpperCase( Locale.ROOT );
        synchronized ( _cache ) {
            if ( _cache.containsKey( k ) ) {
                return _cache.get( k );
            }
        }
        final String tsv = WsHttp.httpGet( BASE + WsHttp.encode( accession.trim() ) );
        final SequenceEntry parsed = parseTsv( tsv );
        if ( !parsed.isEmpty() ) {
            _cache.put( k, parsed );
            return parsed;
        }
        // Cache a negative ONLY for a well-formed "no match" (the header row, no data row). An empty,
        // truncated, or non-TSV body (e.g. a transient 200 error page) is left uncached so a later
        // attempt can retry rather than stranding a resolvable accession for the whole process.
        if ( isHeaderOnly( tsv ) ) {
            _cache.put( k, SequenceEntry.EMPTY );
        }
        return SequenceEntry.EMPTY;
    }

    /** A UniProt "no match" response is exactly the header row (which contains tabs) with no data row. */
    private static boolean isHeaderOnly( final String tsv ) {
        if ( ForesterUtil.isEmpty( tsv ) ) {
            return false;
        }
        String first = null;
        int nonblank = 0;
        for( final String line : tsv.split( "\n" ) ) {
            if ( !line.trim().isEmpty() ) {
                if ( first == null ) {
                    first = line;
                }
                ++nonblank;
            }
        }
        return ( nonblank == 1 ) && ( first != null ) && first.contains( "\t" );
    }

    /**
     * Parses a UniProt REST TSV response (a header row + zero/one data row) into a {@link SequenceEntry}.
     * Columns are read positionally in the order of {@link #FIELDS}. UniProtKB entries are proteins.
     * Pure -- no I/O. Returns {@link SequenceEntry#EMPTY} when there is no data row.
     */
    static SequenceEntry parseTsv( final String tsv ) {
        if ( ForesterUtil.isEmpty( tsv ) ) {
            return SequenceEntry.EMPTY;
        }
        final String[] lines = tsv.split( "\n" );
        if ( lines.length < 2 ) {
            return SequenceEntry.EMPTY; // header only (or nothing) -> no entry
        }
        // limit -1 keeps trailing empty columns (e.g. an empty last xref field)
        final String[] c = lines[ 1 ].split( "\t", -1 );
        final SequenceEntry d = new SequenceEntry( col( c, 0 ), // accession
                                                   col( c, 1 ), // id (entry name)
                                                   col( c, 2 ), // protein_name -> sequence name
                                                   geneName( col( c, 3 ), col( c, 4 ) ), // gene_primary | gene_names
                                                   col( c, 5 ), // sequence
                                                   SequenceEntry.MoleculeType.PROTEIN, // UniProtKB entries are proteins
                                                   organismName( col( c, 6 ) ), // organism_name (strip "(Common name)")
                                                   col( c, 7 ), // organism_id
                                                   lineageNames( col( c, 8 ) ), // lineage (strip each "(rank)")
                                                   splitList( col( c, 9 ), ";" ), // go_id
                                                   crossRefs( col( c, 10 ), col( c, 11 ) ), // refseq, pdb
                                                   null ); // no source-feature properties from UniProt
        return d.isEmpty() ? SequenceEntry.EMPTY : d;
    }

    private static String col( final String[] cols, final int i ) {
        return ( i < cols.length ) ? emptyToNull( cols[ i ] ) : null;
    }

    private static String emptyToNull( final String s ) {
        if ( s == null ) {
            return null;
        }
        final String t = s.trim();
        return t.isEmpty() ? null : t;
    }

    /**
     * UniProt's lineage cell is "Name (rank), Name (rank), ..." (e.g. "Eukaryota (domain), Metazoa
     * (kingdom)"); split on ',' and drop each trailing "(rank)" so the node lineage holds clean names.
     */
    private static List<String> lineageNames( final String cell ) {
        final List<String> out = new ArrayList<String>();
        if ( !ForesterUtil.isEmpty( cell ) ) {
            for( final String part : cell.split( "," ) ) {
                final String name = organismName( part ); // strips " (...)" and trims; null if empty
                if ( name != null ) {
                    out.add( name );
                }
            }
        }
        return out;
    }

    /** The scientific name: UniProt's organism field is "Scientific name (Common name)"; drop the parenthetical. */
    private static String organismName( final String organism ) {
        if ( organism == null ) {
            return null;
        }
        final int paren = organism.indexOf( '(' );
        final String sci = ( paren > 0 ) ? organism.substring( 0, paren ).trim() : organism.trim();
        return sci.isEmpty() ? null : sci;
    }

    /** Primary gene name when present, else the first token of the (space-separated) gene-names list. */
    private static String geneName( final String gene_primary, final String gene_names ) {
        if ( !ForesterUtil.isEmpty( gene_primary ) ) {
            return gene_primary;
        }
        if ( !ForesterUtil.isEmpty( gene_names ) ) {
            final String first = gene_names.trim().split( "\\s+" )[ 0 ];
            return first.isEmpty() ? null : first;
        }
        return null;
    }

    /** Splits a delimited UniProt cell into a trimmed, non-empty list (e.g. lineage on ",", GO on ";"). */
    private static List<String> splitList( final String cell, final String delim ) {
        final List<String> out = new ArrayList<String>();
        if ( !ForesterUtil.isEmpty( cell ) ) {
            for( final String part : cell.split( java.util.regex.Pattern.quote( delim ) ) ) {
                final String t = part.trim();
                if ( !t.isEmpty() ) {
                    out.add( t );
                }
            }
        }
        return out;
    }

    /** RefSeq + PDB cross-references as "db:id" strings (each cell is ";"-separated, often with a trailing ";"). */
    private static List<String> crossRefs( final String refseq, final String pdb ) {
        final List<String> out = new ArrayList<String>();
        for( final String id : splitList( refseq, ";" ) ) {
            out.add( "refseq:" + id );
        }
        for( final String id : splitList( pdb, ";" ) ) {
            out.add( "pdb:" + id );
        }
        return out;
    }
}
