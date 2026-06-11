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

package org.forester.applications;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;

public class get_distances {

    public static void main( final String args[] ) {
        if ( args.length != 3 ) {
            System.out.println( "\nget_distances: Wrong number of arguments.\n" );
            System.out.println( "Usage: \"get_distances <phylogeny file> <file with node names> <outfile>\"\n" );
            System.exit( -1 );
        }
        final File phylogeny_infile = new File( args[ 0 ] );
        final File names_infile = new File( args[ 1 ] );
        final File outfile = new File( args[ 2 ] );
        Phylogeny p = null;
        try {
            final PhylogenyParser pp = org.forester.io.parsers.util.ParserUtils
                    .createParserDependingOnFileType( phylogeny_infile, true );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            p = factory.create( phylogeny_infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + phylogeny_infile + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        String line = "";
        try {
            final BufferedReader in = new BufferedReader( new FileReader( names_infile ) );
            final BufferedWriter out = new BufferedWriter( new FileWriter( outfile ) );
            while ( ( line = in.readLine() ) != null ) {
                if ( line.length() < 3 ) {
                    continue;
                }
                final StringTokenizer st = new StringTokenizer( line );
                if ( st.countTokens() < 2 ) {
                    continue;
                }
                final double d = PhylogenyMethods.calculateDistance( p.getNode( st.nextToken() ),
                                                                     p.getNode( st.nextToken() ) );
                out.write( line + " " + d );
                out.newLine();
            }
            out.flush();
            out.close();
            in.close();
        }
        catch ( final IOException e ) {
            System.out.println( "\nError during processing of \"" + names_infile + "\" [" + e.getMessage()
                                + "] at line \"" + line + "\"\n" );
            System.exit( -1 );
        }
        System.out.println( "\nDone.\n" );
    }
}
