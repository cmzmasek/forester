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

/**
 * Unit tests for the pure XML parsers of {@link NcbiTaxonomyLineageService} -- driven entirely by
 * captured fixtures, so they exercise the lineage/rank extraction with no network access. Lives in
 * the same package to reach the package-private parse methods.
 */
public final class NcbiTaxonomyLineageServiceTest {

    // A realistic esearch response (taxonomy db) for a scientific-name query.
    private static final String ESEARCH = """
            <?xml version="1.0" encoding="UTF-8"?>
            <eSearchResult><Count>1</Count><RetMax>1</RetMax><RetStart>0</RetStart>
            <IdList><Id>9682</Id></IdList></eSearchResult>""";

    // A realistic efetch taxonomy response for the genus Felis. Note the DOCTYPE referencing a remote
    // DTD: parsing must NOT reach out for it (success here proves external-DTD loading is disabled),
    // and the trailing "no rank" entry must be skipped.
    private static final String EFETCH_FELIS = """
            <?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE TaxaSet PUBLIC "-//NLM//DTD Taxon, 14th January 2002//EN" "https://www.ncbi.nlm.nih.gov/entrez/query/DTD/taxon.dtd">
            <TaxaSet>
              <Taxon>
                <TaxId>9682</TaxId>
                <ScientificName>Felis</ScientificName>
                <OtherNames><GenbankCommonName>cats</GenbankCommonName></OtherNames>
                <Rank>genus</Rank>
                <LineageEx>
                  <Taxon><TaxId>131567</TaxId><ScientificName>cellular organisms</ScientificName><Rank>no rank</Rank></Taxon>
                  <Taxon><TaxId>2759</TaxId><ScientificName>Eukaryota</ScientificName><Rank>superkingdom</Rank></Taxon>
                  <Taxon><TaxId>33208</TaxId><ScientificName>Metazoa</ScientificName><Rank>kingdom</Rank></Taxon>
                  <Taxon><TaxId>7711</TaxId><ScientificName>Chordata</ScientificName><Rank>phylum</Rank></Taxon>
                  <Taxon><TaxId>40674</TaxId><ScientificName>Mammalia</ScientificName><Rank>class</Rank></Taxon>
                  <Taxon><TaxId>33554</TaxId><ScientificName>Carnivora</ScientificName><Rank>order</Rank></Taxon>
                  <Taxon><TaxId>9681</TaxId><ScientificName>Felidae</ScientificName><Rank>family</Rank></Taxon>
                  <Taxon><TaxId>338152</TaxId><ScientificName>Felinae</ScientificName><Rank>subfamily</Rank></Taxon>
                </LineageEx>
              </Taxon>
            </TaxaSet>""";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NcbiTaxonomyLineageService: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        // esearch id extraction
        if ( !"9682".equals( NcbiTaxonomyLineageService.parseEsearchFirstId( ESEARCH ) ) ) {
            return fail( "esearch first id not extracted" );
        }
        if ( NcbiTaxonomyLineageService.parseEsearchFirstId(
                "<eSearchResult><IdList></IdList></eSearchResult>" ) != null ) {
            return fail( "empty IdList must yield no id" );
        }
        if ( NcbiTaxonomyLineageService.parseEsearchFirstId( "" ) != null ) {
            return fail( "empty esearch must yield no id" );
        }

        // efetch extraction -> ONE TaxonLineage serving every consumer (colorizer + Fetch tool)
        final TaxonLineage rl = NcbiTaxonomyLineageService.parseEfetchFull( EFETCH_FELIS );
        if ( ( rl == null ) || rl.isEmpty() ) {
            return fail( "Felis lineage did not parse (was external DTD fetched?)" );
        }
        // rank -> name lookup (the core of the colorize-by-rank fix: Felis -> order Carnivora)
        if ( !"Carnivora".equals( rl.at( "order" ) ) ) {
            return fail( "expected order Carnivora, got " + rl.at( "order" ) );
        }
        if ( !"Carnivora".equals( rl.at( "ORDER" ) ) ) {
            return fail( "rank lookup must be case-insensitive" );
        }
        if ( !"Felis".equals( rl.at( "genus" ) ) ) { // the queried taxon's own rank/name
            return fail( "expected genus Felis, got " + rl.at( "genus" ) );
        }
        if ( !"Felidae".equals( rl.at( "family" ) ) || !"Mammalia".equals( rl.at( "class" ) )
                || !"Metazoa".equals( rl.at( "kingdom" ) ) || !"Eukaryota".equals( rl.at( "superkingdom" ) ) ) {
            return fail( "an expected ancestor rank was missing/wrong" );
        }
        if ( rl.at( "species" ) != null ) {
            return fail( "Felis has no species rank in this lineage" );
        }
        // the "no rank" cellular-organisms entry must be skipped, not stored under some rank
        if ( rl.namesByRank().containsValue( "cellular organisms" ) ) {
            return fail( "'no rank' entry should have been skipped" );
        }
        // rank -> tax-id lookup (Spine A's new capability): ancestors' ids, the taxon's OWN id, and the alias
        if ( !"33554".equals( rl.taxIdAt( "order" ) ) || !"40674".equals( rl.taxIdAt( "class" ) )
                || !"9682".equals( rl.taxIdAt( "genus" ) ) ) { // 9682 is Felis's own id
            return fail( "per-rank tax-ids not plumbed: order=" + rl.taxIdAt( "order" ) + " genus="
                    + rl.taxIdAt( "genus" ) );
        }
        // domain<->superkingdom alias: NCBI labels the domain level "superkingdom", so a "domain" query must
        // resolve via the superkingdom entry (Eukaryota / 2759) for BOTH name and tax-id, and vice versa
        if ( !"Eukaryota".equals( rl.at( "domain" ) ) || !"2759".equals( rl.taxIdAt( "domain" ) ) ) {
            return fail( "at/taxIdAt(\"domain\") must alias to the superkingdom entry; got " + rl.at( "domain" ) );
        }
        final TaxonLineage domain_only = new TaxonLineage( null, null, null, null,
                java.util.Collections.singletonList( new TaxonLineage.Ancestor( "Bacteria", "domain", "2" ) ) );
        if ( !"Bacteria".equals( domain_only.at( "superkingdom" ) ) || !"2".equals( domain_only.taxIdAt( "superkingdom" ) ) ) {
            return fail( "at(\"superkingdom\") must alias to a domain entry when present" );
        }

        // Fetch-tool facet on the SAME object: own scientific name/rank/tax-id/common + the COMPLETE lineage of
        // scientific names (incl. the no-rank ancestor), ending in the taxon itself
        if ( !"Felis".equals( rl.getScientificName() ) || !"genus".equals( rl.ownRankOrNull() )
                || !"9682".equals( rl.getTaxId() ) || !"cats".equals( rl.getCommonName() ) ) {
            return fail( "own-taxon detail fields wrong: " + rl );
        }
        if ( !rl.lineageNames().contains( "Carnivora" )
                || !"Felis".equals( rl.lineageNames().get( rl.lineageNames().size() - 1 ) )
                || !rl.lineageNames().contains( "cellular organisms" ) ) {
            return fail( "lineageNames must be the full name path ending in Felis; got " + rl.lineageNames() );
        }

        // edge cases: empty / no-Taxon / malformed efetch all yield an empty TaxonLineage, never throw
        if ( !NcbiTaxonomyLineageService.parseEfetchFull( "" ).isEmpty()
                || !NcbiTaxonomyLineageService.parseEfetchFull( "<TaxaSet></TaxaSet>" ).isEmpty()
                || !NcbiTaxonomyLineageService.parseEfetchFull( "<broken" ).isEmpty() ) {
            return fail( "empty/no-Taxon/malformed efetch must yield an empty TaxonLineage, not throw" );
        }

        // request throttling (keeps the clients under NCBI's ~3 req/sec limit -> avoids HTTP 429)
        if ( WsHttp.throttleDelayMs( 1000L, 1100L, 350L ) != 250L ) {
            return fail( "throttle must wait the remaining interval (250ms)" );
        }
        if ( WsHttp.throttleDelayMs( 1000L, 1350L, 350L ) != 0L ) {
            return fail( "throttle must not wait once the interval has elapsed (boundary)" );
        }
        if ( WsHttp.throttleDelayMs( 1000L, 9999L, 350L ) != 0L ) {
            return fail( "throttle must never return a negative wait" );
        }
        if ( WsHttp.throttleDelayMs( 0L, 0L, 350L ) != 350L ) {
            return fail( "the first request waits the full interval from a zero baseline" );
        }

        // tax-id detection: a bare numeric query is efetched directly (no esearch); anything else is a name
        if ( !NcbiTaxonomyLineageService.isTaxId( "9606" ) || !NcbiTaxonomyLineageService.isTaxId( " 9682 " ) ) {
            return fail( "a run of digits (optionally padded) must be recognized as a tax-id" );
        }
        if ( NcbiTaxonomyLineageService.isTaxId( "Felis" ) || NcbiTaxonomyLineageService.isTaxId( "P12345" )
                || NcbiTaxonomyLineageService.isTaxId( "9606A" ) || NcbiTaxonomyLineageService.isTaxId( "" ) ) {
            return fail( "names / accessions / mixed strings must not be treated as a tax-id" );
        }
        return testEsearchIds() && testEmbeddedTaxId() && testStagedSearch();
    }

    /** Every id of an IdList, in order; none for an empty list, a malformed response or an error-only response. */
    private static boolean testEsearchIds() {
        if ( !java.util.Arrays.asList( "55087", "1386" ).equals( NcbiTaxonomyLineageService.parseEsearchIds( esearchXml( "55087", "1386" ) ) ) ) {
            return fail( "parseEsearchIds must return every id in NCBI's order" );
        }
        if ( !NcbiTaxonomyLineageService.parseEsearchIds( esearchXml() ).isEmpty()
                || !NcbiTaxonomyLineageService.parseEsearchIds( "<broken" ).isEmpty()
                || !NcbiTaxonomyLineageService.parseEsearchIds( "" ).isEmpty() ) {
            return fail( "an empty, malformed or blank esearch must yield no ids" );
        }
        // an <Id> outside the IdList (e.g. inside an ErrorList) is not a hit
        if ( !NcbiTaxonomyLineageService.parseEsearchIds(
                "<eSearchResult><Count>0</Count><IdList/><ErrorList><Id>1</Id></ErrorList></eSearchResult>" ).isEmpty() ) {
            return fail( "only ids under IdList are hits" );
        }
        return true;
    }

    /** The UniProt "OX=" tax-id carried in a label: at the start, after whitespace or after a '|'; never part of
     *  another token, never a partial number. */
    private static boolean testEmbeddedTaxId() {
        if ( !"47803".equals( NcbiTaxonomyLineageService.embeddedTaxId(
                "tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Radical SAM protein OS=Fusarium phyllophilum OX=47803 GN=FPHYL_5758 PE=4 SV=1" ) ) ) {
            return fail( "the OX= tax-id of a UniProt header must be read" );
        }
        if ( !"9606".equals( NcbiTaxonomyLineageService.embeddedTaxId( "OX=9606" ) )
                || !"9606".equals( NcbiTaxonomyLineageService.embeddedTaxId( "P12345|OX=9606" ) ) ) {
            return fail( "OX= at the start or after a '|' must be read" );
        }
        if ( ( NcbiTaxonomyLineageService.embeddedTaxId( "BOX=9606" ) != null )
                || ( NcbiTaxonomyLineageService.embeddedTaxId( "OX=96a" ) != null )
                || ( NcbiTaxonomyLineageService.embeddedTaxId( "Felis catus" ) != null )
                || ( NcbiTaxonomyLineageService.embeddedTaxId( "" ) != null )
                || ( NcbiTaxonomyLineageService.embeddedTaxId( null ) != null ) ) {
            return fail( "a token merely ending in OX=, a non-numeric value, a plain name or nothing carry no tax-id" );
        }
        return true;
    }

    /**
     * The staged, quoted name search and the "never guess" rule, driven by an in-memory transport: the first stage
     * with a hit decides; one hit resolves, several leave the name unplaced with its candidates on record; a tax-id
     * in the label skips the search entirely; a name's own brackets travel inside the quotes.
     */
    private static boolean testStagedSearch() {
        final FakeHttp http = new FakeHttp();
        http.esearch( "\"Felis\"[Scientific Name]", "9682" );
        http.esearch( "\"Bacillus\"[Scientific Name]", "55087", "1386" );
        http.esearch( "\"Bacillus\"[All Names]", "55087", "1386", "99999" ); // must never be asked
        http.esearch( "\"human\"[All Names]", "9606" );
        http.esearch( "Bacillus subtilis (strain 168)", "1423" ); // only the unrestricted stage knows this form
        http.esearch( "\"[Haemophilus] ducreyi\"[Scientific Name]", "730" );
        http.efetch( "9682", EFETCH_FELIS );
        http.efetch( "55087,1386", efetchOf( taxonXml( "55087", "Bacillus", "genus", "Arthropoda" ),
                                            taxonXml( "1386", "Bacillus", "genus", "Bacillota" ) ) );
        http.efetch( "9606", efetchOf( taxonXml( "9606", "Homo sapiens", "species", "Chordata" ) ) );
        http.efetch( "1423", efetchOf( taxonXml( "1423", "Bacillus subtilis", "species", "Bacillota" ) ) );
        http.efetch( "47803", efetchOf( taxonXml( "47803", "Fusarium phyllophilum", "species", "Ascomycota" ) ) );
        http.efetch( "730", efetchOf( taxonXml( "730", "[Haemophilus] ducreyi", "species", "Pseudomonadota" ) ) );
        final NcbiTaxonomyLineageService svc;
        try {
            svc = new NcbiTaxonomyLineageService( http,
                    new TaxonomyDiskCache( java.nio.file.Files.createTempDirectory( "aptx-taxonomy-test" ) ) );
        }
        catch ( final IOException e ) {
            return fail( "could not create a temporary cache directory: " + e );
        }
        try {
            // (1) the scientific-name stage resolves Felis; the later stages are never asked
            if ( !"Carnivora".equals( svc.fetch( "Felis" ).at( "order" ) ) ) {
                return fail( "Felis must resolve at the scientific-name stage" );
            }
            if ( http.asked( "\"Felis\"[All Names]" ) || http.asked( "Felis" ) ) {
                return fail( "a hit at the first stage must end the search" );
            }
            // (2) a homonym is AMBIGUOUS: left unplaced (an EMPTY, cached), its candidates fetched in ONE efetch and
            //     exposed, no later stage asked, no second query ever
            if ( !svc.fetch( "Bacillus" ).isEmpty() ) {
                return fail( "a name several taxa carry must NOT resolve to the first of them" );
            }
            if ( http.asked( "\"Bacillus\"[All Names]" ) || http.asked( "Bacillus" ) ) {
                return fail( "an ambiguous stage must not fall through to a looser one" );
            }
            final java.util.List<TaxonLineage> cands = svc.ambiguousCandidates( " bacillus " );
            if ( ( cands.size() != 2 ) || !"Arthropoda".equals( cands.get( 0 ).at( "phylum" ) )
                    || !"Bacillota".equals( cands.get( 1 ).at( "phylum" ) ) || !"55087".equals( cands.get( 0 ).getTaxId() ) ) {
                return fail( "the candidates of an ambiguous name must be on record, each with its lineage; got " + cands );
            }
            if ( ( svc.lineageOf( "Bacillus" ) == null ) || !svc.lineageOf( "Bacillus" ).isEmpty() ) {
                return fail( "an ambiguous name must be cached as unplaced so it is not re-prompted" );
            }
            final int calls = http.calls();
            svc.fetch( "Bacillus" );
            if ( http.calls() != calls ) {
                return fail( "a cached ambiguity must not be queried again" );
            }
            if ( !svc.ambiguousCandidates( "Felis" ).isEmpty() ) {
                return fail( "a resolved name has no candidates" );
            }
            // (3) a common name is not a scientific name: it falls to the All Names stage
            if ( !"Homo sapiens".equals( svc.fetch( "human" ).getScientificName() ) ) {
                return fail( "a common name must resolve at the All Names stage" );
            }
            if ( !http.asked( "\"human\"[Scientific Name]" ) ) {
                return fail( "the scientific-name stage must be tried first" );
            }
            // (4) a form NCBI lists differently falls to the unrestricted stage, after both quoted stages
            if ( !"Bacillus subtilis".equals( svc.fetch( "Bacillus subtilis (strain 168)" ).getScientificName() ) ) {
                return fail( "a strain-suffixed name must resolve at the unrestricted stage" );
            }
            if ( !http.asked( "\"Bacillus subtilis (strain 168)\"[Scientific Name]" )
                    || !http.asked( "\"Bacillus subtilis (strain 168)\"[All Names]" ) ) {
                return fail( "both quoted stages must be tried before the unrestricted one" );
            }
            // (5) nothing at any stage: a cached negative, and NOT ambiguous
            if ( !svc.fetch( "Nonexistus fictus" ).isEmpty() || !svc.ambiguousCandidates( "Nonexistus fictus" ).isEmpty() ) {
                return fail( "an unknown name is unplaced and has no candidates" );
            }
            if ( !http.asked( "Nonexistus fictus" ) ) {
                return fail( "all three stages must be tried before giving up" );
            }
            // (6) a tax-id carried in the label is used as is: no name search at all, result cached under the label
            final int searches = http.searches();
            final String header = "tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Radical SAM protein OS=Fusarium phyllophilum OX=47803 GN=X PE=4 SV=1";
            if ( !"Fusarium phyllophilum".equals( svc.fetch( header ).getScientificName() ) ) {
                return fail( "a label carrying OX= must resolve by that tax-id" );
            }
            if ( http.searches() != searches ) {
                return fail( "an embedded tax-id must skip the name search" );
            }
            if ( ( svc.lineageOf( header ) == null ) || svc.lineageOf( header ).isEmpty() ) {
                return fail( "the lineage must be cached under the label it was asked for" );
            }
            // (7) a name's own brackets travel INSIDE the quotes, not as search syntax
            if ( !"730".equals( svc.fetch( "[Haemophilus] ducreyi" ).getTaxId() ) ) {
                return fail( "a bracketed name must be sent quoted and resolve" );
            }
        }
        catch ( final IOException e ) {
            return fail( "the in-memory transport must not throw: " + e );
        }
        return true;
    }

    private static String esearchXml( final String... ids ) {
        final StringBuilder sb = new StringBuilder( "<?xml version=\"1.0\"?><eSearchResult><Count>" ).append( ids.length )
                .append( "</Count><IdList>" );
        for( final String id : ids ) {
            sb.append( "<Id>" ).append( id ).append( "</Id>" );
        }
        return sb.append( "</IdList></eSearchResult>" ).toString();
    }

    private static String taxonXml( final String id, final String sci, final String rank, final String phylum ) {
        return "<Taxon><TaxId>" + id + "</TaxId><ScientificName>" + sci + "</ScientificName><Rank>" + rank
                + "</Rank><LineageEx><Taxon><TaxId>1</TaxId><ScientificName>" + phylum
                + "</ScientificName><Rank>phylum</Rank></Taxon></LineageEx></Taxon>";
    }

    private static String efetchOf( final String... taxa ) {
        return "<?xml version=\"1.0\"?><TaxaSet>" + String.join( "", taxa ) + "</TaxaSet>";
    }

    /** An in-memory transport: esearch answers keyed by the DECODED term, efetch answers by the id list; records
     *  what was asked. An unknown term finds nothing; an unknown id fetches nothing. */
    private static final class FakeHttp implements NcbiTaxonomyLineageService.Http {

        private final java.util.Map<String, String> _esearch = new java.util.HashMap<String, String>();
        private final java.util.Map<String, String> _efetch  = new java.util.HashMap<String, String>();
        private final java.util.List<String>        _terms   = new java.util.ArrayList<String>();
        private int                                 _calls;

        void esearch( final String term, final String... ids ) {
            _esearch.put( term, esearchXml( ids ) );
        }

        void efetch( final String ids, final String xml ) {
            _efetch.put( ids, xml );
        }

        boolean asked( final String term ) {
            return _terms.contains( term );
        }

        int calls() {
            return _calls;
        }

        int searches() {
            return _terms.size();
        }

        @Override
        public String get( final String url ) throws IOException {
            ++_calls;
            try {
                final int t = url.indexOf( "term=" );
                if ( t >= 0 ) {
                    final String term = java.net.URLDecoder.decode( url.substring( t + 5 ), "UTF-8" );
                    _terms.add( term );
                    return _esearch.containsKey( term ) ? _esearch.get( term ) : esearchXml();
                }
                final int i = url.indexOf( "id=" );
                final String ids = java.net.URLDecoder.decode( url.substring( i + 3 ), "UTF-8" );
                return _efetch.containsKey( ids ) ? _efetch.get( ids ) : efetchOf();
            }
            catch ( final java.io.UnsupportedEncodingException e ) {
                throw new IOException( e );
            }
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  NcbiTaxonomyLineageService test failed: " + msg );
        return false;
    }

    private NcbiTaxonomyLineageServiceTest() {
    }
}
