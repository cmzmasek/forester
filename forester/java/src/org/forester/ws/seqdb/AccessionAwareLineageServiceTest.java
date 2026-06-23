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
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

import org.forester.phylogeny.data.Accession;

/**
 * Unit tests for {@link AccessionAwareLineageService} -- the bridge that lets the rank colorizer /
 * "Annotate Clades by Rank" place tips identified by a UniProt/SwissProt/RefSeq sequence accession.
 * Driven by in-memory fakes (no network): an {@link OrganismSource} that maps an accession to an
 * organism tax-id, and a {@link TaxonomicLineageService} delegate that maps a tax-id/name to a ranked
 * lineage. Lives in the same package to use the package-private seams.
 */
public final class AccessionAwareLineageServiceTest {

    public static void main( final String[] args ) {
        System.out.println( "AccessionAwareLineageService: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        final FakeDelegate delegate = new FakeDelegate();
        delegate.know( "9606", "order", "Primates" );        // human, keyed by NCBI tax-id
        delegate.know( "10090", "order", "Rodentia" );       // mouse, keyed by NCBI tax-id
        delegate.know( "Felis", "order", "Carnivora" );      // a plain scientific name
        delegate.know( "Mus musculus", "order", "Rodentia" );// organism scientific name (id-less entry)
        final FakeOrganisms seqs = new FakeOrganisms();
        seqs.know( "P12345", "9606", "Homo sapiens" );        // UniProt accession -> human
        seqs.know( "RL7_HUMAN", "9606", "Homo sapiens" );     // SwissProt entry name -> human
        seqs.know( "Q12340", null, "Mus musculus" );          // id-less organism: only a scientific name
        final AccessionAwareLineageService svc = new AccessionAwareLineageService( delegate, seqs );

        try {
            // (1) a UniProt accession bridges acc -> tax-id 9606 -> ranked lineage; lineageOf hits afterward
            if ( svc.lineageOf( "P12345" ) != null ) {
                return fail( "lineageOf must be null before the accession is fetched" );
            }
            if ( !"Primates".equals( svc.fetch( "P12345" ).at( "order" ) ) ) {
                return fail( "UniProt accession must resolve to its organism's order (Primates)" );
            }
            if ( ( svc.lineageOf( "P12345" ) == null ) || !"Primates".equals( svc.lineageOf( "P12345" ).at( "order" ) ) ) {
                return fail( "after fetch, the cache-only lineageOf must return the bridged lineage under the accession key" );
            }

            // (2) a SwissProt entry name (mnemonic) resolves the same way
            if ( !"Primates".equals( svc.fetch( "RL7_HUMAN" ).at( "order" ) ) ) {
                return fail( "a SwissProt entry name must resolve via the sequence DB too" );
            }

            // (3) an id-less entry falls back to resolving the organism scientific name
            if ( !"Rodentia".equals( svc.fetch( "Q12340" ).at( "order" ) ) ) {
                return fail( "an entry without a tax-id must fall back to its organism name (Mus musculus -> Rodentia)" );
            }

            // (4) a plain scientific name is NOT treated as an accession -- it passes straight to the delegate
            final int seq_calls_before = seqs.calls();
            if ( !"Carnivora".equals( svc.fetch( "Felis" ).at( "order" ) ) ) {
                return fail( "a plain scientific name must resolve through the taxonomy delegate" );
            }
            if ( seqs.calls() != seq_calls_before ) {
                return fail( "a plain name must not hit the organism source" );
            }

            // (5) a not-found accession caches a negative: lineageOf returns EMPTY (non-null) and a re-fetch
            //     does not hit the organism source again
            final RankedLineage miss = svc.fetch( "P99988" );
            if ( ( miss == null ) || !miss.isEmpty() ) {
                return fail( "an unresolvable accession must return an empty (not null) lineage" );
            }
            if ( ( svc.lineageOf( "P99988" ) == null ) || !svc.lineageOf( "P99988" ).isEmpty() ) {
                return fail( "a not-found accession must cache a negative so it is not re-fetched/re-prompted" );
            }
            final int calls_after_miss = seqs.calls();
            svc.fetch( "P99988" );
            if ( seqs.calls() != calls_after_miss ) {
                return fail( "a cached-negative accession must not be fetched from the organism source again" );
            }

            // (6) empty / blank queries are safe
            if ( ( svc.fetch( "" ) != RankedLineage.EMPTY ) || ( svc.lineageOf( null ) != null ) ) {
                return fail( "empty/blank queries must be handled without a fetch" );
            }
        }
        catch ( final IOException e ) {
            return fail( "fakes must not throw: " + e );
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  AccessionAwareLineageService test failed: " + msg );
        return false;
    }

    /** In-memory taxonomy delegate: {@code lineageOf} is cache-only; {@code fetch} copies from the "DB". */
    private static final class FakeDelegate implements TaxonomicLineageService {

        private final Map<String, RankedLineage> _db    = new HashMap<String, RankedLineage>();
        private final Map<String, RankedLineage> _cache = new HashMap<String, RankedLineage>();

        void know( final String key, final String rank, final String name ) {
            final Map<String, String> m = new LinkedHashMap<String, String>();
            m.put( rank, name );
            _db.put( key.toLowerCase( Locale.ROOT ), new RankedLineage( m ) );
        }

        @Override
        public RankedLineage lineageOf( final String taxon ) {
            return ( taxon == null ) ? null : _cache.get( taxon.trim().toLowerCase( Locale.ROOT ) );
        }

        @Override
        public RankedLineage fetch( final String taxon ) {
            final String k = taxon.trim().toLowerCase( Locale.ROOT );
            final RankedLineage rl = _db.containsKey( k ) ? _db.get( k ) : RankedLineage.EMPTY;
            _cache.put( k, rl );
            return rl;
        }
    }

    /** In-memory organism source: maps an accession value to an organism (tax-id + scientific name); counts calls. */
    private static final class FakeOrganisms implements OrganismSource {

        private final Map<String, Organism> _db = new HashMap<String, Organism>();
        private int                         _calls;

        void know( final String accession, final String tax_id, final String scientific_name ) {
            _db.put( accession.toUpperCase( Locale.ROOT ), new Organism( tax_id, scientific_name ) );
        }

        @Override
        public Organism organismOf( final Accession acc ) {
            ++_calls;
            final Organism o = _db.get( acc.getValue().toUpperCase( Locale.ROOT ) );
            return ( o == null ) ? Organism.EMPTY : o;
        }

        int calls() {
            return _calls;
        }
    }

    private AccessionAwareLineageServiceTest() {
    }
}
