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
import java.util.Collections;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import org.forester.phylogeny.data.Accession;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;

/**
 * A {@link TaxonomicLineageService} that also resolves <b>sequence accessions</b> -- UniProt/SwissProt
 * (e.g. {@code P12345}, {@code RL7_HUMAN}), RefSeq, GenBank, GI -- which a bare taxonomy database cannot
 * place because they name a <i>protein</i>, not an organism. Trees whose tips are UniProt and/or mixed
 * NCBI/UniProt identifiers are very common, so the rank colorizer and "Annotate Clades by Rank" need to
 * place them too.
 *
 * <p>It decorates a taxonomy {@code delegate} (NCBI) with a {@link SequenceFetcher}: when a query looks
 * like a sequence accession, it fetches that entry, reads the organism's NCBI tax-id, and resolves the
 * ranked lineage for that tax-id through the delegate (an exact efetch -- the dominant case is many
 * proteins of a few organisms, so the per-organism lineage is fetched once and reused). The result is
 * cached under the <i>accession's own query string</i> so the cache-only {@link #lineageOf} the EDT
 * uses hits exactly the key the colorizer looks up. Plain scientific names (and bare tax-ids) pass
 * straight through to the delegate.
 *
 * <p>Network I/O happens in {@link #fetch} only -- call it off the EDT; {@link #lineageOf} is cache-only.
 * The seam design (delegate + fetcher) makes it unit-testable with in-memory fakes.
 */
public final class AccessionAwareLineageService implements TaxonomicLineageService {

    private final TaxonomicLineageService            _delegate;
    private final SequenceFetcher                    _sequences;
    // accession query string -> ranked lineage (including EMPTY negatives), so a later cache-only
    // lineageOf(accession) hits and the accession is never re-fetched.
    private final Map<String, RankedLineage>         _alias = Collections
            .synchronizedMap( new HashMap<String, RankedLineage>() );

    public AccessionAwareLineageService( final TaxonomicLineageService delegate, final SequenceFetcher sequences ) {
        _delegate = delegate;
        _sequences = sequences;
    }

    @Override
    public RankedLineage lineageOf( final String taxon ) {
        if ( ForesterUtil.isEmpty( taxon ) ) {
            return null;
        }
        final RankedLineage rl = _delegate.lineageOf( taxon );
        if ( rl != null ) {
            return rl; // a plain name already resolved through the delegate
        }
        return _alias.get( key( taxon ) ); // an accession we bridged earlier (null if not yet fetched)
    }

    @Override
    public RankedLineage fetch( final String taxon ) throws IOException {
        if ( ForesterUtil.isEmpty( taxon ) ) {
            return RankedLineage.EMPTY;
        }
        final RankedLineage cached = lineageOf( taxon );
        if ( cached != null ) {
            return cached;
        }
        // UniProt-priority parsing matches Aptx's own accession handling: a UniProt accession such as
        // P12345 also fits the single-letter+5-digit GenBank shape, and protein trees are UniProt-first.
        final Accession acc = SequenceAccessionTools.parseAccessorFromString_UniProtPriority( taxon );
        if ( acc != null ) {
            return fetchViaAccession( taxon, acc );
        }
        return _delegate.fetch( taxon ); // a plain scientific name or bare tax-id
    }

    /**
     * Resolves {@code acc} (the accession parsed from {@code query}) to its organism's ranked lineage via
     * the sequence DB + the taxonomy delegate, caching the result under {@code query}. A genuine miss is
     * cached as a negative so the accession is not re-fetched (and the user is not re-prompted); a
     * transport failure propagates so the caller can stop the run cleanly.
     */
    private RankedLineage fetchViaAccession( final String query, final Accession acc ) throws IOException {
        final SequenceEntry entry = _sequences.fetch( acc );
        RankedLineage rl = RankedLineage.EMPTY;
        if ( ( entry != null ) && !entry.isEmpty() ) {
            if ( !ForesterUtil.isEmpty( entry.getOrganismId() ) ) {
                rl = _delegate.fetch( entry.getOrganismId() ); // exact efetch by NCBI tax-id
            }
            if ( ( ( rl == null ) || rl.isEmpty() ) && !ForesterUtil.isEmpty( entry.getOrganismName() ) ) {
                rl = _delegate.fetch( entry.getOrganismName() ); // fall back to the organism scientific name
            }
        }
        if ( rl == null ) {
            rl = RankedLineage.EMPTY;
        }
        _alias.put( key( query ), rl );
        return rl;
    }

    private static String key( final String s ) {
        return s.trim().toLowerCase( Locale.ROOT );
    }
}
