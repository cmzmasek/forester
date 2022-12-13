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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public final class blast_pars {

    public static void main( final String args[] ) {
        try {
            final File infile = new File( args[ 0 ] );
            String desc = "";
            String acc = "";
            boolean saw_description = false;
            boolean saw_query = false;
            try (BufferedReader br = new BufferedReader( new FileReader( infile ) )) {
                String line;
                while ( ( line = br.readLine() ) != null ) {
                    line = line.trim();
                    if ( line.length() > 0 ) {
                        if ( line.startsWith( "Description" ) ) {
                            saw_description = true;
                            saw_query = false;
                            desc = "";
                            acc = "";
                        }
                        else if ( saw_description && ( desc.length() == 0 ) ) {
                            final String[] s = line.split( "\\s{2,}" );
                            desc = s[ 0 ];
                            acc = s[ s.length - 1 ];
                            System.out.println( ">" + desc + "|" + acc );
                        }
                        else if ( saw_description && ( desc.length() > 0 ) && line.startsWith( "Query" ) ) {
                            saw_query = true;
                            final String[] s = line.split( "\\s+" );
                            final String seq = s[ 2 ];
                            System.out.println( seq );
                        }
                        else if ( saw_query && line.startsWith( ">" ) ) {
                            saw_description = false;
                            saw_query = false;
                        }
                    }
                }
            }
            //            final String name = args[ 2 ];
            //            DeleteableMsa msa = null;
            //            final FileInputStream is = new FileInputStream( infile );
            //            if ( FastaParser.isLikelyFasta( infile ) ) {
            //                msa = DeleteableMsa.createInstance( FastaParser.parseMsa( is ) );
            //            }
            //            else {
            //                msa = DeleteableMsa.createInstance( GeneralMsaParser.parseMsa( is ) );
            //            }
            //            final String cons = MsaMethods.calculateMajorityConsensusSequence( msa );
            //            final List<MolecularSequence> seqs = new ArrayList<>();
            //            seqs.add( BasicSequence.createGeneralSequence( name, cons ) );
            //            if ( special ) {
            //                final boolean all_identical = MsaMethods.isAllSequencesIdentical( msa );
            //                if ( !all_identical ) {
            //                    final String out = outfile.toString();
            //                    outfile = new File( out.substring( 0, out.length() - 6 ) + "_TO_BLAST.fasta" );
            //                }
            //            }
            //            SequenceWriter.writeSeqs( seqs, outfile, SEQ_FORMAT.FASTA, 80 );
        }
        catch ( final FileNotFoundException e ) {
            e.printStackTrace();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }
}
