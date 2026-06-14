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
 * Resolves a taxon's full {@link RankedLineage} from a taxonomy database, with a clean split
 * between a <b>blocking</b> {@link #fetch(String)} (does network I/O; call off the EDT) and a
 * <b>non-blocking</b> {@link #lineageOf(String)} (cache-only; safe on the EDT). This lets a GUI
 * resolve everything it can from the cache instantly, then fetch only the misses in the background
 * and re-read the cache -- never blocking the event thread on the network.
 *
 * <p>Implementations cache results (including negatives) so repeated questions about the same taxon
 * cost nothing. This interface is the seam that lets the colorizer be unit-tested with an in-memory
 * fake instead of the network.
 */
public interface TaxonomicLineageService {

    /**
     * The cached lineage for {@code taxon}, or {@code null} if it has not been {@link #fetch fetched}
     * yet (or was fetched and not found). Never does network I/O -- safe to call on the EDT.
     */
    RankedLineage lineageOf( String taxon );

    /**
     * Resolves {@code taxon} against the database (network I/O) and caches the result, returning it.
     * Returns {@link RankedLineage#EMPTY} (cached) when the taxon cannot be resolved, so it is not
     * re-queried. Call off the EDT.
     *
     * @throws IOException on a connection/transport failure (a genuinely-not-found taxon is not an
     *             error -- it returns {@link RankedLineage#EMPTY}).
     */
    RankedLineage fetch( String taxon ) throws IOException;
}
