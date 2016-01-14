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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class table2fasta {

    final static private String PRG_NAME    = "table2fasta";
    final static private String PRG_VERSION = "1.00";
    final static private String PRG_DATE    = "150327";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( table2fasta.PRG_NAME, table2fasta.PRG_VERSION, table2fasta.PRG_DATE );
        System.out.println();
        if ( ( args.length != 3 ) ) {
            table2fasta.argumentsError();
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        final int position = Integer.parseInt( cla.getName( 0 ) );
        final File intable = cla.getFile( 1 );
        final File outfile = cla.getFile( 2 );
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, outfile + " already exists" );
        }
        if ( !intable.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, intable + " does not exist" );
        }
        BasicTable<String> t = null;
        try {
            t = BasicTableParser.parse( intable, '\t' );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        final List<MolecularSequence> seqs = new ArrayList<MolecularSequence>();
        System.out.println( "Number of rows: " + t.getNumberOfRows() );
        for( int r = 0; r < t.getNumberOfRows(); ++r ) {
            String seq = null;
            final StringBuilder id = new StringBuilder();
            for( int c = 0; c < t.getNumberOfColumns(); ++c ) {
                if ( c == position ) {
                    seq = t.getValue( c, r );
                }
                else if ( ( c == 0 ) || ( c == 1 ) ) {
                    id.append( t.getValue( c, r ) );
                    id.append( " " );
                }
            }
            if ( id.length() < 2 ) {
                ForesterUtil.fatalError( PRG_NAME, "row " + r + " id is empty" );
            }
            String id_str = id.toString().trim();
            if ( id_str.startsWith( ">" ) ) {
                id_str = id_str.substring( 1 );
            }
            if ( ForesterUtil.isEmpty( seq ) ) {
                seq = t.getValue( position - 1, r );
                if ( ForesterUtil.isEmpty( seq ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "row " + r + " seq is empty" );
                }
            }
            MolecularSequence s = null;
            try {
                s = BasicSequence.createAaSequence( id_str, seq );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            }
            if ( s != null ) {
                seqs.add( s );
            }
        }
        try {
            SequenceWriter.writeSeqs( seqs, outfile, SEQ_FORMAT.FASTA, 60 );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    private static void argumentsError() {
        System.out.println( PRG_NAME + " <position> <infile> <outfile>" );
        System.out.println();
        System.exit( -1 );
    }
}
