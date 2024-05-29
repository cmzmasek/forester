// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester
//
//
// "java -Xmx1024m -cp path\to\forester.jar org.forester.application.fasta_split
//
//

package org.forester.application;

import org.forester.io.parsers.FastaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.ViralUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class rename_fasta_acc {

    public final static Pattern PATTERN_GB = Pattern
            .compile("\\|([A-Z][A-Z0-9.]{5,})");

    public static void main( final String args[] ) {
        try {
            final File infile = new File( args[ 0 ] );
            final File outfile = new File( args[ 1 ] );
            List<MolecularSequence> seqs;
            seqs = FastaParser.parse( new FileInputStream( infile ) );
            for( MolecularSequence seq : seqs ) {
                final BasicSequence bseq = ( BasicSequence ) seq;
                final Matcher mg = PATTERN_GB.matcher(bseq.getIdentifier());
                if (mg.find()) {
                    bseq.setIdentifier( mg.group(1).trim() );
                }
                else {
                    System.out.println( bseq.getIdentifier() + " did not match");
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
