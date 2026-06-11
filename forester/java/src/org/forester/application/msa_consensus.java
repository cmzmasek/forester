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

package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.msa.DeleteableMsa;
import org.forester.msa.MsaMethods;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;

public final class msa_consensus {

    public static void main( final String args[] ) {
        final boolean special = false;
        try {
            final File infile = new File( args[ 0 ] );
            File outfile = new File( args[ 1 ] );
            final String name = args[ 2 ];
            DeleteableMsa msa = null;
            final FileInputStream is = new FileInputStream( infile );
            if ( FastaParser.isLikelyFasta( infile ) ) {
                msa = DeleteableMsa.createInstance( FastaParser.parseMsa( is ) );
            }
            else {
                msa = DeleteableMsa.createInstance( GeneralMsaParser.parseMsa( is ) );
            }
            final String cons = MsaMethods.calculateMajorityConsensusSequence( msa );
            final List<MolecularSequence> seqs = new ArrayList<>();
            seqs.add( BasicSequence.createGeneralSequence( name, cons ) );
            if ( special ) {
                final boolean all_identical = MsaMethods.isAllSequencesIdentical( msa );
                if ( !all_identical ) {
                    final String out = outfile.toString();
                    outfile = new File( out.substring( 0, out.length() - 6 ) + "_TO_BLAST.fasta" );
                }
            }
            SequenceWriter.writeSeqs( seqs, outfile, SEQ_FORMAT.FASTA, 80 );
        }
        catch ( final FileNotFoundException e ) {
            e.printStackTrace();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }
}
