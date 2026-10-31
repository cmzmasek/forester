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

import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;

/**
 * Where a value in the node window ("Display Node Data") links to, and where the node menu's "Open Taxonomy DB"
 * goes. Pure: each method returns a URL, or null when the data does not name a page we know how to address. Only
 * ever a LINK, never a fetch: Archaeopteryx asks no server for any of this.
 * <p>
 * The rules (Christian, 2026-10-30: links go to "the provider's database"):
 * <ul>
 * <li>A taxonomy IDENTIFIER opens the database its provider names: NCBI Taxonomy for ncbi (and for no provider,
 * which {@link Identifier#normalizedProvider} reads as ncbi), UniProt Taxonomy for uniprot. An identifier that is
 * itself a web address opens that address. Any other provider links nothing.</li>
 * <li>A SCIENTIFIC NAME opens its identifier's page when there is one, else NCBI Taxonomy's lookup by name.</li>
 * <li>A taxonomy CODE is a UniProt mnemonic ("HUMAN") and opens UniProt's taxonomy search for it.</li>
 * <li>A COMMON NAME links nothing ("bat", "E. coli": no page is named), so the node menu never opens one either.</li>
 * <li>A sequence ACCESSION opens the database its SOURCE names, and only when the value is shaped like that
 * database's accessions. With no source -- or the placeholders forester writes for none, "unknown" (phyloXML)
 * and "?" (NHX) -- only a value that IS a UniProt, RefSeq or GenBank accession links. A source we do not know
 * links nothing: a BV-BRC "fig|1000.9.peg.1" is not guessed at.</li>
 * <li>At NCBI the accession's FORMAT picks Protein or Nucleotide, not the sequence's declared type. Measured
 * 2026-10-30 on 725 trees (Downloads, the family reference trees, the demos): 2,287 sequences declared "protein"
 * carry a genome accession (KJ473821, NC_009021). The declared type decides only a format we do not recognise.</li>
 * <li>A DOI opens doi.org. A URI, or a metadata value, opens itself when it is an http(s) address.</li>
 * </ul>
 * Every pattern was opened once by hand (2026-10-30) and answered with the page it names.
 */
final class NodeDataLinks {

    static final String NCBI_TAXONOMY_BY_ID     = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=";
    static final String NCBI_TAXONOMY_BY_NAME   = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name=";
    static final String UNIPROT_TAXONOMY        = "https://www.uniprot.org/taxonomy/";
    static final String UNIPROT_TAXONOMY_BY_CODE = "https://www.uniprot.org/taxonomy?query=mnemonic:";
    static final String NCBI_GENE               = "https://www.ncbi.nlm.nih.gov/gene/";
    static final String ENA                     = "https://www.ebi.ac.uk/ena/browser/view/";
    static final String ENSEMBL                 = "https://www.ensembl.org/id/";
    static final String INTERPRO_ENTRY          = "https://www.ebi.ac.uk/interpro/entry/InterPro/";
    static final String PATHOPLEXUS             = "https://pathoplexus.org/seq/";
    static final String DOI                     = "https://doi.org/";

    /** The leading run of characters an accession is made of; what follows must be a separator (":NS5"). */
    private static final Pattern TOKEN              = Pattern.compile( "[A-Za-z0-9_.]+" );
    private static final Pattern SEPARATOR          = Pattern.compile( "[\\s:|/;,]" );
    /** RefSeq: one of NCBI's prefixes and an underscore (a Pathoplexus "PP_..." is not one); the second letter P
     *  is a protein (AP_, NP_, WP_, XP_, YP_). */
    private static final Pattern REFSEQ             = Pattern
            .compile( "(AC|AP|NC|NG|NM|NP|NR|NT|NW|NZ|WP|XM|XP|XR|YP)_[A-Z0-9]+(\\.\\d+)?" );
    /** GenBank protein: three letters and five or seven digits. */
    private static final Pattern GENBANK_PROTEIN    = Pattern.compile( "[A-Z]{3}\\d{5}(\\d{2})?(\\.\\d+)?" );
    /** GenBank nucleotide: 1+5, 2+6, 2+8, and the WGS/TSA/MGA shapes (4+8..10, 5+7, 6+9..11). */
    private static final Pattern GENBANK_NUCLEOTIDE = Pattern
            .compile( "([A-Z]\\d{5}|[A-Z]{2}\\d{6}|[A-Z]{2}\\d{8}|[A-Z]{4}\\d{8,10}|[A-Z]{5}\\d{7}|[A-Z]{6}\\d{9,11})"
                    + "(\\.\\d+)?" );
    private static final Pattern UNIPROT            = Pattern
            .compile( SequenceAccessionTools.UNIPROT_KB_BASE_PATTERN_STR + "(\\.\\d+)?" );
    /** A UniProt entry name: a mnemonic or an accession, an underscore, a species mnemonic (BRCA1_HUMAN). */
    private static final Pattern UNIPROT_ENTRY_NAME = Pattern.compile( "[A-Z0-9]{1,10}_[A-Z0-9]{3,5}" );
    /** An accession of a shape we cannot place: letters, then digits. Only the declared type can place it. */
    private static final Pattern ACCESSION_SHAPE    = Pattern.compile( "[A-Z]{1,6}_?\\d{5,}(\\.\\d+)?" );
    private static final Pattern DIGITS             = Pattern.compile( "\\d+" );
    private static final Pattern INTERPRO           = Pattern.compile( "IPR\\d{6}" );
    private static final Pattern ENSEMBL_ID         = Pattern.compile( "(ENS[A-Z]*\\d+)(\\.\\d+)?" );
    private static final Pattern PATHOPLEXUS_ID     = Pattern.compile( "PP_[A-Z0-9]+(\\.\\d+)?" );
    private static final Pattern DOI_PATTERN        = Pattern.compile( "10\\.\\d{4,9}/\\S+" );
    private static final Pattern WEB                = Pattern.compile( "(?i)https?://\\S+" );
    /** What a URI may hold as written, besides a "%XX" escape (RFC 3986's unreserved and reserved characters, but
     *  '[' and ']', which only an IPv6 host may hold). */
    private static final String  URI_CHARS          = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
            + "-._~:/?#@!$&'()*+,;=";
    private static final char[]  HEX                = "0123456789ABCDEF".toCharArray();

    private NodeDataLinks() {
    }

    // ------------------------------------------------------------------ taxonomy
    /** The page a taxonomy identifier names, or null (see the class comment). */
    static String taxonomyIdUrl( final String value, final String provider ) {
        final String v = trim( value );
        if ( v.isEmpty() ) {
            return null;
        }
        final String web = webUrl( v );
        if ( web != null ) {
            return web;
        }
        if ( !DIGITS.matcher( v ).matches() ) {
            return null;
        }
        final String p = trim( provider ).toLowerCase( Locale.ROOT );
        if ( p.equals( "uniprot" ) || p.equals( "uniprot.taxonomy" ) ) {
            return UNIPROT_TAXONOMY + v;
        }
        return "ncbi".equals( Identifier.normalizedProvider( p ) ) ? NCBI_TAXONOMY_BY_ID + v : null;
    }

    /** A scientific name's page: its identifier's ({@code id_url}) when there is one, else NCBI's lookup by name. */
    static String scientificNameUrl( final String name, final String id_url ) {
        if ( trim( name ).isEmpty() ) {
            return null;
        }
        return ( id_url != null ) ? id_url : NCBI_TAXONOMY_BY_NAME + encodeQuery( trim( name ) );
    }

    /** A UniProt taxonomy mnemonic's search page, or null when the code is not shaped like one. */
    static String taxonomyCodeUrl( final String code ) {
        final String c = trim( code );
        // phyloXML's taxonomy code IS the UniProt mnemonic: one pattern for both
        return PhyloXmlUtil.TAXOMONY_CODE_PATTERN.matcher( c ).matches() ? UNIPROT_TAXONOMY_BY_CODE + c : null;
    }

    /**
     * The node menu's "Open Taxonomy DB": the first of the identifier's page, the scientific name's and the code's
     * -- always a page the node window links -- or null when none names a page.
     */
    static String taxonomyUrl( final Taxonomy t ) {
        if ( t == null ) {
            return null;
        }
        final Identifier id = t.getIdentifier();
        final String id_url = ( id == null ) ? null : taxonomyIdUrl( id.getValue(), id.getProvider() );
        if ( id_url != null ) {
            return id_url;
        }
        final String name_url = scientificNameUrl( t.getScientificName(), null );
        return ( name_url != null ) ? name_url : taxonomyCodeUrl( t.getTaxonomyCode() );
    }

    // ------------------------------------------------------------------ sequences
    /**
     * The page an accession names, from its {@code source}; {@code seq_type} (the sequence's declared type, may be
     * null) places only an NCBI accession whose format we do not recognise. Null when nothing names a page.
     */
    static String accessionUrl( final String value, final String source, final String seq_type ) {
        final String v = trim( value );
        if ( v.isEmpty() ) {
            return null;
        }
        final String s = trim( source ).toLowerCase( Locale.ROOT );
        if ( s.isEmpty() || s.equals( "unknown" ) || s.equals( "?" ) ) { // what forester itself writes for none
            return accessionUrlWithoutSource( v );
        }
        final String token = token( v );
        if ( token == null ) {
            return null;
        }
        switch ( s ) {
            case "ncbi":
            case "genbank":
            case "refseq":
            case "insdc":
            case "ddbj":
                return ncbiUrl( token, seq_type, token.equals( v ) );
            case "gi":
                return DIGITS.matcher( token ).matches() ? ForesterUtil.NCBI_GI + token : null;
            case "geneid":
            case "entrez":
            case "entrezgene":
            case "ncbigene":
            case "ncbi_gene":
                return DIGITS.matcher( token ).matches() ? NCBI_GENE + token : null;
            case "embl":
            case "ena":
                return ACCESSION_SHAPE.matcher( token ).matches() ? ENA + token : null;
            case "ensembl": {
                final Matcher m = ENSEMBL_ID.matcher( token );
                return m.matches() ? ENSEMBL + m.group( 1 ) : null;
            }
            case "interpro":
                return INTERPRO.matcher( token ).matches() ? INTERPRO_ENTRY + token + "/" : null;
            case "pfam": {
                final String acc = NodeHoverText.pfamAccession( token );
                return ( acc != null ) ? TreePanel.pfamUrl( acc ) : null;
            }
            case "pathoplexus":
                return PATHOPLEXUS_ID.matcher( token ).matches() ? PATHOPLEXUS + token : null;
            default:
                if ( s.startsWith( "uniprot" ) || s.equals( "swissprot" ) || s.equals( "swiss-prot" )
                        || s.equals( "trembl" ) || s.equals( "sp" ) || s.equals( "tr" ) ) {
                    // an accession or an entry name; not the "sp" of "sp|P38398|BRCA1_HUMAN"
                    return ( UNIPROT.matcher( token ).matches() || UNIPROT_ENTRY_NAME.matcher( token ).matches() )
                            ? ForesterUtil.UNIPROT_KB + withoutUniProtVersion( token ) : null;
                }
                return null;
        }
    }

    /** No source: only a value that IS a UniProt, RefSeq or GenBank accession links (UniProt first, as the
     *  node menu's "Open Sequence DB" reads one). */
    private static String accessionUrlWithoutSource( final String v ) {
        if ( UNIPROT.matcher( v ).matches() ) {
            return ForesterUtil.UNIPROT_KB + withoutUniProtVersion( v );
        }
        if ( REFSEQ.matcher( v ).matches() || GENBANK_PROTEIN.matcher( v ).matches()
                || GENBANK_NUCLEOTIDE.matcher( v ).matches() ) {
            return ncbiUrl( v, null, true );
        }
        return null;
    }

    /** NCBI Protein or Nucleotide by the accession's format; the declared type only for a shape we cannot place,
     *  and then only when the token is the whole value. */
    private static String ncbiUrl( final String token, final String seq_type, final boolean whole_value ) {
        if ( REFSEQ.matcher( token ).matches() ) {
            return ( ( token.charAt( 1 ) == 'P' ) ? ForesterUtil.NCBI_PROTEIN : ForesterUtil.NCBI_NUCCORE ) + token;
        }
        if ( GENBANK_PROTEIN.matcher( token ).matches() ) {
            return ForesterUtil.NCBI_PROTEIN + token;
        }
        final boolean nucleotide = GENBANK_NUCLEOTIDE.matcher( token ).matches();
        final boolean uniprot = UNIPROT.matcher( token ).matches();
        final String type = trim( seq_type ).toLowerCase( Locale.ROOT );
        if ( nucleotide && uniprot ) { // [OPQ] and five digits: either; the declared type decides, protein if none
            return ( type.equals( "dna" ) || type.equals( "rna" ) ? ForesterUtil.NCBI_NUCCORE
                    : ForesterUtil.NCBI_PROTEIN ) + token;
        }
        if ( nucleotide ) {
            return ForesterUtil.NCBI_NUCCORE + token;
        }
        if ( uniprot ) { // a Swiss-Prot accession under an NCBI source: NCBI Protein carries Swiss-Prot
            return ForesterUtil.NCBI_PROTEIN + token;
        }
        if ( whole_value && ACCESSION_SHAPE.matcher( token ).matches() ) {
            if ( type.equals( "protein" ) ) {
                return ForesterUtil.NCBI_PROTEIN + token;
            }
            if ( type.equals( "dna" ) || type.equals( "rna" ) ) {
                return ForesterUtil.NCBI_NUCCORE + token;
            }
        }
        return null;
    }

    /** The accession at the start of {@code v}: the leading run of accession characters, provided the value
     *  ends there or goes on after a separator ("QEA08541.1:NS5" gives QEA08541.1). Null otherwise. */
    private static String token( final String v ) {
        final Matcher m = TOKEN.matcher( v );
        if ( !m.lookingAt() ) {
            return null;
        }
        if ( ( m.end() < v.length() ) && !SEPARATOR.matcher( v.substring( m.end(), m.end() + 1 ) ).matches() ) {
            return null;
        }
        String t = v.substring( 0, m.end() );
        while ( t.endsWith( "." ) ) {
            t = t.substring( 0, t.length() - 1 );
        }
        return t.isEmpty() ? null : t;
    }

    /** UniProt addresses an entry without a sequence version: "Q5UQR0.1" is Q5UQR0. */
    private static String withoutUniProtVersion( final String token ) {
        final Matcher m = UNIPROT.matcher( token );
        if ( m.matches() && ( m.group( m.groupCount() ) != null ) ) {
            return token.substring( 0, token.length() - m.group( m.groupCount() ).length() );
        }
        return token;
    }

    // ------------------------------------------------------------------ the rest
    /** doi.org for a DOI, written bare or as "doi:..." or as a doi.org address; null for anything else. */
    static String doiUrl( final String doi ) {
        String d = trim( doi );
        final String lower = d.toLowerCase( Locale.ROOT );
        for( final String prefix : new String[] { "https://doi.org/", "http://doi.org/", "https://dx.doi.org/",
                "http://dx.doi.org/", "doi:" } ) {
            if ( lower.startsWith( prefix ) ) {
                d = d.substring( prefix.length() ).trim();
                break;
            }
        }
        if ( !DOI_PATTERN.matcher( d ).matches() ) {
            return null;
        }
        // the DOI is the PATH: a '?' or '#' in it is part of the name, not a query or a fragment
        final String path = escape( d, "?#" );
        return ( path == null ) ? null : parsedWebAddress( DOI + path );
    }

    /**
     * A text that is one http(s) address and nothing else, as an address the browser can be handed; null
     * otherwise. Characters a URI may not hold ('|', '{', non-ASCII) are escaped; the address's own "%XX" escapes
     * are kept as written (re-quoting would turn "%20" into "%2520").
     */
    static String webUrl( final String text ) {
        final String t = trim( text );
        return WEB.matcher( t ).matches() ? parsedWebAddress( escape( t, "" ) ) : null;
    }

    /** {@code escaped} as a URI with a host, or null (no host; not well-formed Unicode; not a URI). */
    private static String parsedWebAddress( final String escaped ) {
        if ( escaped == null ) {
            return null;
        }
        try {
            final URI uri = new URI( escaped );
            return ( uri.getHost() != null ) ? uri.toASCIIString() : null;
        }
        catch ( final URISyntaxException e ) {
            return null;
        }
    }

    /**
     * {@code s} with every character a URI may not hold, and those in {@code also}, percent-encoded as UTF-8; an
     * existing "%XX" escape is kept, a '#' after the first is escaped. Null when {@code s} is not well-formed
     * Unicode (an unpaired surrogate, which {@code URI} cannot encode at all).
     */
    private static String escape( final String s, final String also ) {
        final StringBuilder sb = new StringBuilder( s.length() + 16 );
        boolean seen_hash = false;
        for( int i = 0; i < s.length(); ) {
            final int cp = s.codePointAt( i );
            final int len = Character.charCount( cp );
            if ( ( cp >= Character.MIN_SURROGATE ) && ( cp <= Character.MAX_SURROGATE ) ) {
                return null;
            }
            if ( ( cp == '%' ) && ( ( i + 2 ) < s.length() ) && isHex( s.charAt( i + 1 ) ) && isHex( s.charAt( i + 2 ) ) ) {
                sb.append( s, i, i + 3 );
                i += 3;
                continue;
            }
            final boolean keep = ( cp < 0x80 ) && ( URI_CHARS.indexOf( cp ) >= 0 ) && ( also.indexOf( cp ) < 0 )
                    && !( ( cp == '#' ) && seen_hash );
            if ( keep ) {
                sb.append( (char) cp );
                seen_hash |= ( cp == '#' );
            }
            else {
                for( final byte b : new String( Character.toChars( cp ) ).getBytes( StandardCharsets.UTF_8 ) ) {
                    sb.append( '%' ).append( HEX[ ( b >> 4 ) & 0xF ] ).append( HEX[ b & 0xF ] );
                }
            }
            i += len;
        }
        return sb.toString();
    }

    private static boolean isHex( final char c ) {
        return ( ( c >= '0' ) && ( c <= '9' ) ) || ( ( c >= 'A' ) && ( c <= 'F' ) ) || ( ( c >= 'a' ) && ( c <= 'f' ) );
    }

    private static String encodeQuery( final String s ) {
        try {
            return URLEncoder.encode( s, ForesterConstants.UTF_8 );
        }
        catch ( final UnsupportedEncodingException e ) {
            throw new IllegalStateException( e ); // UTF-8 is always there
        }
    }

    private static String trim( final String s ) {
        return ( s == null ) ? "" : s.trim();
    }
}
