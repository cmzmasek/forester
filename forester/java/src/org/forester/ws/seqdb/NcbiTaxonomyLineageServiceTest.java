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

        // efetch ranked-lineage extraction (the core of the colorize-by-rank fix: Felis -> order Carnivora)
        final RankedLineage rl = NcbiTaxonomyLineageService.parseEfetchTaxonomyXml( EFETCH_FELIS );
        if ( ( rl == null ) || rl.isEmpty() ) {
            return fail( "Felis lineage did not parse (was external DTD fetched?)" );
        }
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
        if ( rl.asMap().containsValue( "cellular organisms" ) ) {
            return fail( "'no rank' entry should have been skipped" );
        }
        // domain<->superkingdom alias: NCBI labels the domain level "superkingdom", so at("domain") must
        // resolve via the superkingdom entry (Eukaryota), and vice versa
        if ( !"Eukaryota".equals( rl.at( "domain" ) ) ) {
            return fail( "at(\"domain\") must alias to the superkingdom entry; got " + rl.at( "domain" ) );
        }
        final java.util.Map<String, String> dm = new java.util.LinkedHashMap<String, String>();
        dm.put( "domain", "Bacteria" );
        if ( !"Bacteria".equals( new RankedLineage( dm ).at( "superkingdom" ) ) ) {
            return fail( "at(\"superkingdom\") must alias to a domain entry when present" );
        }

        // edge cases
        if ( !NcbiTaxonomyLineageService.parseEfetchTaxonomyXml( "" ).isEmpty() ) {
            return fail( "empty efetch must yield empty lineage" );
        }
        if ( !NcbiTaxonomyLineageService.parseEfetchTaxonomyXml( "<TaxaSet></TaxaSet>" ).isEmpty() ) {
            return fail( "efetch with no Taxon must yield empty lineage" );
        }
        if ( !NcbiTaxonomyLineageService.parseEfetchTaxonomyXml( "<broken" ).isEmpty() ) {
            return fail( "malformed efetch must yield empty lineage, not throw" );
        }

        // full taxonomy detail (the "Fetch Sequence & Taxonomic Data" path): scientific name, rank,
        // tax-id, common name, and the COMPLETE lineage of scientific names (incl. the no-rank ancestor)
        final ResolvedTaxonomy rt = NcbiTaxonomyLineageService.parseEfetchTaxonomyDetail( EFETCH_FELIS );
        if ( ( rt == null ) || rt.isEmpty() ) {
            return fail( "Felis detail did not parse" );
        }
        if ( !"Felis".equals( rt.getScientificName() ) || !"genus".equals( rt.getRank() )
                || !"9682".equals( rt.getTaxId() ) || !"cats".equals( rt.getCommonName() ) ) {
            return fail( "detail fields wrong: " + rt );
        }
        if ( !rt.getLineage().contains( "Carnivora" )
                || !"Felis".equals( rt.getLineage().get( rt.getLineage().size() - 1 ) )
                || !rt.getLineage().contains( "cellular organisms" ) ) {
            return fail( "detail lineage must be the full name path ending in Felis; got " + rt.getLineage() );
        }
        if ( !NcbiTaxonomyLineageService.parseEfetchTaxonomyDetail( "" ).isEmpty()
                || !NcbiTaxonomyLineageService.parseEfetchTaxonomyDetail( "<broken" ).isEmpty() ) {
            return fail( "empty/malformed efetch detail must yield an empty ResolvedTaxonomy, not throw" );
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
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  NcbiTaxonomyLineageService test failed: " + msg );
        return false;
    }

    private NcbiTaxonomyLineageServiceTest() {
    }
}
