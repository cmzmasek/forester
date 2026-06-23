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
import org.forester.phylogeny.data.Accession.Source;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;

/**
 * The production {@link OrganismSource}: routes a UniProt accession to {@link UniProtKbClient#organism}
 * (a minimal, taxonomy-only UniProt request that caches nothing heavy), and any other accession
 * (RefSeq/EMBL/GI/NCBI) to {@link NcbiSequenceClient}, keeping only the organism from the result.
 *
 * <p>The UniProt path -- the dominant case for protein trees, and the one the rank-colorizer bridge most
 * needs -- never pulls or caches the full protein record. The NCBI efetch path reuses the shared
 * sequence client (so it does cache the GBSeq record); that less-common path could be made lightweight
 * later, but the organism is all this class exposes.
 */
public final class WebOrganismSource implements OrganismSource {

    private final UniProtKbClient    _uniprot = new UniProtKbClient();
    private final NcbiSequenceClient _ncbi    = new NcbiSequenceClient();

    @Override
    public Organism organismOf( final Accession acc ) throws IOException {
        if ( ( acc == null ) || ForesterUtil.isEmpty( acc.getValue() ) ) {
            return Organism.EMPTY;
        }
        final String value = acc.getValue();
        if ( Source.UNIPROT.toString().equalsIgnoreCase( acc.getSource() ) ) {
            return _uniprot.organism( value );
        }
        // RefSeq / EMBL / NCBI / GI -> NCBI efetch; the protein-vs-nucleotide shape guess can be wrong, so
        // if the guessed database has no entry, try the other before giving up.
        final boolean protein = SequenceAccessionTools.isProteinDbQuery( value );
        SequenceEntry e = _ncbi.fetchEntry( value, protein );
        if ( e.isEmpty() ) {
            e = _ncbi.fetchEntry( value, !protein );
        }
        return e.isEmpty() ? Organism.EMPTY : new Organism( e.getOrganismId(), e.getOrganismName() );
    }
}
