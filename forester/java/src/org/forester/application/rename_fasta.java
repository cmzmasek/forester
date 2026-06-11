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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.FastaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class rename_fasta {

    

    public static void main( final String args[] ) {
        try {
            final File infile = new File( args[ 0 ] );
            final File outfile = new File( args[ 1 ] );
            List<MolecularSequence> seqs;
            seqs = FastaParser.parse( new FileInputStream( infile ) );
            for( MolecularSequence seq : seqs ) {
                BasicSequence bseq = ( BasicSequence ) seq;
                final int i = bseq.getIdentifier().lastIndexOf( '_' );
                if ( i > 0 ) {
                    bseq.setIdentifier( bseq.getIdentifier().substring( i + 1 ) );
                }
            }
            SequenceWriter.writeSeqs( seqs, outfile, SEQ_FORMAT.FASTA, 60 );
        }
        catch ( FileNotFoundException e ) {
            e.printStackTrace();
        }
        catch ( IOException e ) {
            e.printStackTrace();
        }
    }

  
}
