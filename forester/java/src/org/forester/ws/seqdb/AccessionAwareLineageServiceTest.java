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
            final TaxonLineage miss = svc.fetch( "P99988" );
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
            if ( ( svc.fetch( "" ) != TaxonLineage.EMPTY ) || ( svc.lineageOf( null ) != null ) ) {
                return fail( "empty/blank queries must be handled without a fetch" );
            }

            // (7) a tax-id carried IN the label (the UniProt header's "OX=") goes straight to the taxonomy delegate --
            //     no organism lookup, although the header also carries an accession -- and is remembered under the label
            delegate.know( "47803", "order", "Hypocreales" );
            final String header = "tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Radical SAM protein OS=Fusarium phyllophilum OX=47803 GN=F PE=4 SV=1";
            final int before_ox = seqs.calls();
            if ( !"Hypocreales".equals( svc.fetch( header ).at( "order" ) ) ) {
                return fail( "a label carrying OX= must resolve by that tax-id through the delegate" );
            }
            if ( seqs.calls() != before_ox ) {
                return fail( "an embedded tax-id must not cost a sequence-database call" );
            }
            if ( ( svc.lineageOf( header ) == null ) || !"Hypocreales".equals( svc.lineageOf( header ).at( "order" ) ) ) {
                return fail( "the lineage must be cached under the label for the cache-only lineageOf" );
            }

            // (8) ambiguity is the delegate's to report; an accession is never ambiguous
            delegate.ambiguous( "Bacillus", new TaxonLineage( "1386", "genus", "Bacillus", null, null ),
                                new TaxonLineage( "55087", "genus", "Bacillus", null, null ) );
            if ( svc.ambiguousCandidates( "Bacillus" ).size() != 2 ) {
                return fail( "ambiguousCandidates must come from the taxonomy delegate" );
            }
            if ( !svc.ambiguousCandidates( "P12345" ).isEmpty() || !svc.ambiguousCandidates( null ).isEmpty() ) {
                return fail( "an accession or a blank query has no candidates" );
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

        private final Map<String, TaxonLineage>                _db        = new HashMap<String, TaxonLineage>();
        private final Map<String, TaxonLineage>                _cache     = new HashMap<String, TaxonLineage>();
        private final Map<String, java.util.List<TaxonLineage>> _ambiguous = new HashMap<String, java.util.List<TaxonLineage>>();

        void know( final String key, final String rank, final String name ) {
            // the delegate resolves this taxon as its own rank/name (no ancestors); TaxonLineage.at(rank) -> name
            _db.put( key.toLowerCase( Locale.ROOT ), new TaxonLineage( null, rank, name, null, null ) );
        }

        void ambiguous( final String name, final TaxonLineage... candidates ) {
            _ambiguous.put( name.toLowerCase( Locale.ROOT ), java.util.Arrays.asList( candidates ) );
        }

        @Override
        public java.util.List<TaxonLineage> ambiguousCandidates( final String taxon ) {
            final java.util.List<TaxonLineage> c = ( taxon == null ) ? null
                    : _ambiguous.get( taxon.trim().toLowerCase( Locale.ROOT ) );
            return ( c == null ) ? java.util.Collections.<TaxonLineage> emptyList() : c;
        }

        @Override
        public TaxonLineage lineageOf( final String taxon ) {
            return ( taxon == null ) ? null : _cache.get( taxon.trim().toLowerCase( Locale.ROOT ) );
        }

        @Override
        public TaxonLineage fetch( final String taxon ) {
            final String k = taxon.trim().toLowerCase( Locale.ROOT );
            final TaxonLineage rl = _db.containsKey( k ) ? _db.get( k ) : TaxonLineage.EMPTY;
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
