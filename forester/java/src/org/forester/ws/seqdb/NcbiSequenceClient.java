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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.forester.util.ForesterUtil;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A clean NCBI sequence client (GenBank/RefSeq/EMBL/GI accessions) over E-utilities efetch in GBSeq
 * <b>XML</b>, returning a {@link SequenceEntry}. This is the non-UniProt path of "Fetch Sequence &amp;
 * Taxonomic Data"; UniProt accessions go through {@link UniProtKbClient}. Network I/O goes through
 * {@link WsHttp} (timeout-bounded + throttled). Results -- including "not found" -- are cached.
 */
public final class NcbiSequenceClient {

    private static final String              EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?tool=Archaeopteryx&rettype=gb&retmode=xml&db=";
    // GenBank "source" feature qualifiers -> phyloXML property refs, so existing property-color schemes
    // keep working. NCBI renamed "country" to "geo_loc_name"; accept both. Refs mirror the prior tool.
    private static final Map<String, String> SOURCE_PROP_REFS = sourcePropRefs();
    private final Map<String, SequenceEntry> _cache = Collections
            .synchronizedMap( new HashMap<String, SequenceEntry>() );

    private static Map<String, String> sourcePropRefs() {
        final Map<String, String> m = new HashMap<String, String>();
        m.put( "strain", "taxonomy:strain" );
        m.put( "host", "source:host" );
        m.put( "country", "source:country" );
        m.put( "geo_loc_name", "source:country" );
        m.put( "isolate", "source:isolate" );
        m.put( "segment", "genome:viral_segment" );
        m.put( "collection_date", "source:collection_date" );
        m.put( "isolation_source", "source:isolation_source" );
        return m;
    }

    /**
     * The GenBank entry for {@code accession} ({@code protein} selects the {@code protein} vs
     * {@code nuccore} efetch database), or {@link SequenceEntry#EMPTY} when not found. Caches the
     * result. Does network I/O -- call off the EDT.
     *
     * @throws IOException on a connection/transport failure.
     */
    public SequenceEntry fetchEntry( final String accession, final boolean protein ) throws IOException {
        if ( ForesterUtil.isEmpty( accession ) ) {
            return SequenceEntry.EMPTY;
        }
        final String db = protein ? "protein" : "nuccore";
        final String k = db + ":" + accession.trim().toUpperCase( Locale.ROOT );
        synchronized ( _cache ) {
            if ( _cache.containsKey( k ) ) {
                return _cache.get( k );
            }
        }
        final String xml = WsHttp.httpGet( EFETCH + db + "&id=" + WsHttp.encode( accession.trim() ) );
        final SequenceEntry parsed = parseGbSeqXml( xml );
        if ( !parsed.isEmpty() ) {
            _cache.put( k, parsed );
            return parsed;
        }
        // Cache a negative ONLY when the response actually parsed (valid XML, just no usable GBSeq = a
        // real "not in this db"); a parse failure (truncated/garbled/HTML body) is likely transient, so
        // leave it uncached rather than strand a resolvable accession for the whole process.
        if ( WsHttp.parseXml( xml ) != null ) {
            _cache.put( k, SequenceEntry.EMPTY );
        }
        return SequenceEntry.EMPTY;
    }

    /**
     * Parses an NCBI efetch GBSeq XML response into a {@link SequenceEntry}: accession/locus,
     * definition (sequence name), molecule type, organism + lineage, the molecular sequence, and -- from
     * the feature-table qualifiers -- the gene name, organism tax-id ({@code db_xref taxon:N}) and other
     * {@code db_xref} cross-references. Pure -- no I/O. {@link SequenceEntry#EMPTY} when nothing usable.
     */
    static SequenceEntry parseGbSeqXml( final String xml ) {
        final Document doc = WsHttp.parseXml( xml );
        if ( doc == null ) {
            return SequenceEntry.EMPTY;
        }
        final Element gbseq = WsHttp.firstChildElement( doc.getDocumentElement(), "GBSeq" );
        if ( gbseq == null ) {
            return SequenceEntry.EMPTY;
        }
        final String locus = trimOrNull( WsHttp.text( WsHttp.firstChildElement( gbseq, "GBSeq_locus" ) ) );
        String accession = trimOrNull( WsHttp.text( WsHttp.firstChildElement( gbseq, "GBSeq_primary-accession" ) ) );
        if ( accession == null ) {
            accession = locus;
        }
        final String definition = trimOrNull( WsHttp.text( WsHttp.firstChildElement( gbseq, "GBSeq_definition" ) ) );
        final String organism = trimOrNull( WsHttp.text( WsHttp.firstChildElement( gbseq, "GBSeq_organism" ) ) );
        final String moltype = WsHttp.text( WsHttp.firstChildElement( gbseq, "GBSeq_moltype" ) );
        final String seq = upperOrNull( WsHttp.text( WsHttp.firstChildElement( gbseq, "GBSeq_sequence" ) ) );
        final SequenceEntry.MoleculeType type = moleculeType( moltype );

        // lineage = the semicolon-separated GBSeq_taxonomy, then the organism itself appended
        final List<String> lineage = new ArrayList<String>();
        final String taxonomy = WsHttp.text( WsHttp.firstChildElement( gbseq, "GBSeq_taxonomy" ) );
        if ( !ForesterUtil.isEmpty( taxonomy ) ) {
            for( final String part : taxonomy.split( ";" ) ) {
                final String t = part.trim();
                if ( !t.isEmpty() ) {
                    lineage.add( t );
                }
            }
        }
        if ( organism != null ) {
            lineage.add( organism );
        }

        // gene name, organism tax-id, cross-references and source-feature metadata from the qualifiers
        String gene = null;
        String organism_id = null;
        final List<String> xrefs = new ArrayList<String>();
        final Map<String, String> props = new LinkedHashMap<String, String>(); // source qualifiers -> property refs
        final Element feature_table = WsHttp.firstChildElement( gbseq, "GBSeq_feature-table" );
        if ( feature_table != null ) {
            for( final Element feature : WsHttp.childElements( feature_table, "GBFeature" ) ) {
                final Element quals = WsHttp.firstChildElement( feature, "GBFeature_quals" );
                if ( quals == null ) {
                    continue;
                }
                for( final Element q : WsHttp.childElements( quals, "GBQualifier" ) ) {
                    final String name = WsHttp.text( WsHttp.firstChildElement( q, "GBQualifier_name" ) );
                    final String value = trimOrNull( WsHttp.text( WsHttp.firstChildElement( q, "GBQualifier_value" ) ) );
                    if ( ( name == null ) || ( value == null ) ) {
                        continue;
                    }
                    if ( ( gene == null ) && "gene".equals( name ) ) {
                        gene = value;
                    }
                    else if ( "db_xref".equals( name ) ) {
                        final int colon = value.indexOf( ':' );
                        if ( colon > 0 ) {
                            final String xdb = value.substring( 0, colon ).trim().toLowerCase( Locale.ROOT );
                            final String xid = value.substring( colon + 1 ).trim();
                            if ( "taxon".equals( xdb ) ) {
                                if ( organism_id == null ) {
                                    organism_id = xid;
                                }
                            }
                            else if ( !xid.isEmpty() ) {
                                xrefs.add( xdb + ":" + xid );
                            }
                        }
                    }
                    else if ( SOURCE_PROP_REFS.containsKey( name ) ) {
                        // strain/host/country/isolate/segment/collection_date/isolation_source (source feature)
                        props.putIfAbsent( SOURCE_PROP_REFS.get( name ), value );
                    }
                }
            }
        }
        final SequenceEntry e = new SequenceEntry( accession, locus, definition, gene, seq, type, organism,
                                                   organism_id, lineage, Collections.<String> emptyList(), xrefs,
                                                   props );
        return e.isEmpty() ? SequenceEntry.EMPTY : e;
    }

    private static SequenceEntry.MoleculeType moleculeType( final String moltype ) {
        if ( ForesterUtil.isEmpty( moltype ) ) {
            return null;
        }
        final String m = moltype.toUpperCase( Locale.ROOT );
        if ( m.contains( "AA" ) ) {
            return SequenceEntry.MoleculeType.PROTEIN;
        }
        if ( m.contains( "RNA" ) ) {
            return SequenceEntry.MoleculeType.RNA;
        }
        return SequenceEntry.MoleculeType.DNA;
    }

    private static String trimOrNull( final String s ) {
        if ( s == null ) {
            return null;
        }
        final String t = s.trim();
        return t.isEmpty() ? null : t;
    }

    private static String upperOrNull( final String s ) {
        final String t = trimOrNull( s );
        return ( t == null ) ? null : t.toUpperCase( Locale.ROOT );
    }
}
