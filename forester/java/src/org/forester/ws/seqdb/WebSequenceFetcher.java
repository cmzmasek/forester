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
 * The production {@link SequenceFetcher}: routes a UniProt accession to {@link UniProtKbClient} (UniProt
 * REST TSV) and any other accession (RefSeq/EMBL/GI/NCBI) to {@link NcbiSequenceClient} (NCBI efetch
 * GBSeq XML), choosing the protein vs. nucleotide database by accession shape. Both clients are
 * timeout-bounded, throttled and cached via {@link WsHttp}.
 */
public final class WebSequenceFetcher implements SequenceFetcher {

    private final UniProtKbClient    _uniprot = new UniProtKbClient();
    private final NcbiSequenceClient _ncbi    = new NcbiSequenceClient();

    @Override
    public SequenceEntry fetch( final Accession acc ) throws IOException {
        if ( ( acc == null ) || ForesterUtil.isEmpty( acc.getValue() ) ) {
            return SequenceEntry.EMPTY;
        }
        final String value = acc.getValue();
        if ( Source.UNIPROT.toString().equalsIgnoreCase( acc.getSource() ) ) {
            return _uniprot.fetchEntry( value );
        }
        // RefSeq / EMBL / NCBI / GI -> NCBI efetch; pick protein vs nucleotide by accession shape, but the
        // shape guess can be wrong (a protein GI is bare digits; AP_/YP_ RefSeq are proteins) -- so if the
        // guessed database has no entry, try the other one before giving up.
        final boolean protein = SequenceAccessionTools.isProteinDbQuery( value );
        final SequenceEntry e = _ncbi.fetchEntry( value, protein );
        return e.isEmpty() ? _ncbi.fetchEntry( value, !protein ) : e;
    }
}
