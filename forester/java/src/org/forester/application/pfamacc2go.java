// $Id:
//
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

package org.forester.application;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import org.forester.go.PfamToGoMapping;
import org.forester.go.PfamToGoParser;

public class pfamacc2go {

    final static private String PRG_NAME = "pfamacc2go";

    public static void main( final String args[] ) {
        if ( args.length != 2 ) {
            printHelp();
            System.exit( -1 );
        }
        final PfamToGoParser p = new PfamToGoParser( new File( args[ 0 ] ) );
        p.setUseAccessors( true );
        List<PfamToGoMapping> pfam2go = null;
        try {
            pfam2go = p.parse();
        }
        catch ( final IOException e ) {
            printHelp();
            e.printStackTrace();
        }
        BufferedReader br = null;
        try {
            br = new BufferedReader( new FileReader( args[ 1 ] ) );
        }
        catch ( final FileNotFoundException e ) {
            printHelp();
            e.printStackTrace();
        }
        String line;
        int total_pfam_ids = 0;
        int mapped_pfam_ids = 0;
        try {
            while ( ( line = br.readLine() ) != null ) {
                line = line.trim();
                if ( ( line.length() > 0 ) && !line.startsWith( "#" ) ) {
                    String[] pfam_ids = null;
                    if ( line.contains( "," ) ) {
                        pfam_ids = line.split( "," );
                    }
                    else {
                        pfam_ids = new String[ 1 ];
                        pfam_ids[ 0 ] = line;
                    }
                    for( final String pfam_id : pfam_ids ) {
                        total_pfam_ids++;
                        boolean mapped = false;
                        for( final PfamToGoMapping pfam_to_go_mapping : pfam2go ) {
                            if ( pfam_to_go_mapping.getKey().getId().equals( pfam_id ) ) {
                                mapped = true;
                                System.out.println( pfam_to_go_mapping.getValue().toString() );
                            }
                        }
                        if ( mapped ) {
                            mapped_pfam_ids++;
                        }
                    }
                }
            }
        }
        catch ( final Exception e ) {
            printHelp();
            e.printStackTrace();
        }
        System.out.println( "# total pfam ids : " + total_pfam_ids );
        System.out.println( "# pfam ids mapped: " + mapped_pfam_ids );
    }

    private static void printHelp() {
        System.out.println();
        System.out.println( PRG_NAME
                + " <pfam2go mapping file> <file with pfam accessors, newline and/or comma separated>" );
        System.out.println();
    }
}
