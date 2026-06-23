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

import org.forester.phylogeny.data.Accession;

/**
 * Resolves which {@link Organism} a sequence accession belongs to -- the narrow seam the rank-colorizer
 * bridge ({@link AccessionAwareLineageService}) needs. It is intentionally taxonomy-only: callers want
 * the organism, not the protein, so an implementation should request and retain only the organism (never
 * the full record / molecular sequence). The production implementation is {@link WebOrganismSource}; the
 * seam lets the bridge be unit-tested with an in-memory fake.
 */
public interface OrganismSource {

    /**
     * The organism for {@code acc}, or {@link Organism#EMPTY} when nothing resolves. Does network I/O --
     * call off the EDT.
     *
     * @throws IOException on a connection/transport failure.
     */
    Organism organismOf( Accession acc ) throws IOException;
}
