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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public final class blast_pars {

    private final static String VERSION = "1.0.0";
    public static void main( final String args[] ) {
        try {
            File infile = null;
            String annot = "";
            if ( args.length == 1 ) {
                infile = new File( args[ 0 ] );
            }
            else if ( args.length == 2 ) {
                annot = args[ 0 ];
                infile = new File( args[ 1 ] );
            }
            else {
                System.err.println( "Usage: blast_pars [annotation tag] <infile>" );
                System.exit( -1 );
            }
            String desc = "";
            String acc = "";
            boolean saw_description = false;
            boolean saw_Sbjct = false;
            try (BufferedReader br = new BufferedReader( new FileReader( infile ) )) {
                String line;
                while ( ( line = br.readLine() ) != null ) {
                    line = line.trim();
                    if ( line.length() > 0 ) {
                        if ( line.startsWith( "Description" ) ) {
                            saw_description = true;
                            saw_Sbjct = false;
                            desc = "";
                            acc = "";
                        }
                        else if ( saw_description && ( desc.length() == 0 ) ) {
                            final String[] s = line.split( "\\s{2,}" );
                            desc = s[ 0 ];
                            acc = s[ s.length - 1 ];
                            if ( annot.length() > 0 ) {
                                System.out.println( ">" + desc + "|" + acc + "|" + annot );
                            }
                            else {
                                System.out.println( ">" + desc + "|" + acc );
                            }
                        }
                        else if ( saw_description && ( desc.length() > 0 ) && line.startsWith( "Sbjct" ) ) {
                            saw_Sbjct = true;
                            final String[] s = line.split( "\\s+" );
                            final String seq = s[ 2 ];
                            System.out.println( seq );
                        }
                        else if ( saw_Sbjct && line.startsWith( ">" ) ) {
                            saw_description = false;
                            saw_Sbjct = false;
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
