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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

public class pfamacc2pfamid {

    final static private String PRG_NAME = "pfamacc2pfamid";

    public static void main( final String args[] ) {
        if ( args.length != 2 ) {
            printHelp();
            System.exit( -1 );
        }
        BufferedReader br = null;
        try {
            br = new BufferedReader( new FileReader( args[ 0 ] ) );
        }
        catch ( final FileNotFoundException e ) {
            printHelp();
            e.printStackTrace();
        }
        String line;
        final Map<String, String> acc_id = new HashMap<String, String>();
        String id = null;
        try {
            while ( ( line = br.readLine() ) != null ) {
                if ( line.startsWith( "#=GF ID" ) ) {
                    if ( id != null ) {
                        System.err.println( "illegal format" );
                        System.exit( -1 );
                    }
                    id = line.substring( 7 ).trim();
                }
                else if ( line.startsWith( "#=GF AC" ) ) {
                    if ( id == null ) {
                        System.err.println( "illegal format" );
                        System.exit( -1 );
                    }
                    String acc = line.substring( 7 ).trim();
                    if ( acc.indexOf( '.' ) > 0 ) {
                        acc = acc.substring( 0, acc.indexOf( '.' ) );
                    }
                    acc_id.put( acc, id );
                    id = null;
                }
                else if ( line.startsWith( "//" ) ) {
                    if ( id != null ) {
                        System.err.println( "illegal format" );
                        System.exit( -1 );
                    }
                }
            }
        }
        catch ( final Exception e ) {
            printHelp();
            e.printStackTrace();
        }
        try {
            br = new BufferedReader( new FileReader( args[ 1 ] ) );
        }
        catch ( final FileNotFoundException e ) {
            printHelp();
            e.printStackTrace();
        }
        int not_found = 0;
        try {
            while ( ( line = br.readLine() ) != null ) {
                line = line.trim();
                if ( ( line.length() > 0 ) && !line.startsWith( "#" ) ) {
                    String[] pfam_accs = null;
                    if ( line.contains( "," ) ) {
                        pfam_accs = line.split( "," );
                    }
                    else {
                        pfam_accs = new String[ 1 ];
                        pfam_accs[ 0 ] = line;
                    }
                    for( final String pfam_acc : pfam_accs ) {
                        if ( acc_id.containsKey( pfam_acc ) ) {
                            System.out.println( acc_id.get( pfam_acc ) );
                        }
                        else {
                            not_found++;
                        }
                    }
                }
            }
        }
        catch ( final Exception e ) {
            printHelp();
            e.printStackTrace();
        }
        System.err.println( "# not found: " + not_found );
    }

    private static void printHelp() {
        System.out.println();
        System.out.println( PRG_NAME + " <Pfam full> <file with pfam accessors, newline and/or comma separated>" );
        System.out.println();
    }
}
