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
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

public final class ortholog_annotate_tab {

    private final static String VERSION = "1.0.0";
    public static void main( final String args[] ) {
        try {
            File infile = null;
            final String annot = "";
            infile = new File( args[ 0 ] );
            final SortedMap<String, SortedMap<String, SortedSet<String>>> ortho_to_protein_ids_map = new TreeMap<>();
            final SortedSet<String> all_genomes = new TreeSet<>();
            try (BufferedReader br = new BufferedReader( new FileReader( infile ) )) {
                String line;
                while ( ( line = br.readLine() ) != null ) {
                    line = line.trim();
                    if ( line.length() > 0 ) {
                        final String[] x = line.split( "\t" );
                        if ( x.length != 3 ) {
                            System.err.println( "Unexpected line: " + line );
                            System.exit( -1 );
                        }
                        final String protein_id = x[ 0 ];
                        final String da = x[ 1 ];
                        final String da_name = x[ 2 ];
                        final String key = da + "@" + da_name;
                        final String[] y = protein_id.split( "\\." );
                        final String genome = y[ 0 ];
                        all_genomes.add( genome );
                        if ( ortho_to_protein_ids_map.containsKey( key ) ) {
                            final SortedMap<String, SortedSet<String>> protein_ids = ortho_to_protein_ids_map
                                    .get( key );
                            if ( protein_ids.containsKey( genome ) ) {
                                protein_ids.get( genome ).add( protein_id );
                            }
                            else {
                                final SortedSet<String> s = new TreeSet<>();
                                s.add( protein_id );
                                protein_ids.put( genome, s );
                            }
                        }
                        else {
                            final SortedMap<String, SortedSet<String>> protein_ids = new TreeMap<>();
                            final SortedSet<String> s = new TreeSet<>();
                            s.add( protein_id );
                            protein_ids.put( genome, s );
                            ortho_to_protein_ids_map.put( key, protein_ids );
                        }
                    }
                }
            }
            System.out.print( '\t' );
            for( final String genome : all_genomes ) {
                System.out.print( '\t' );
                System.out.print( genome );
            }
            System.out.println();
            for( final Entry<String, SortedMap<String, SortedSet<String>>> entry : ortho_to_protein_ids_map
                    .entrySet() ) {
                final String ortho = entry.getKey();
                final SortedMap<String, SortedSet<String>> protein_ids = entry.getValue();
                final String[] o = ortho.split( "@" );
                System.out.print( o[ 1 ] );
                System.out.print( '\t' );
                System.out.print( o[ 0 ] );
                System.out.print( '\t' );
                for( final String genome : all_genomes ) {
                    if ( protein_ids.containsKey( genome ) ) {
                        final SortedSet<String> actual_protein_ids = protein_ids.get( genome );
                        boolean first = true;
                        for( final String protein_id : actual_protein_ids ) {
                            if ( !first ) {
                                System.out.print( ", " );
                            }
                            first = false;
                            System.out.print( protein_id );
                        }
                    }
                    else {
                        System.out.print( "-" );
                    }
                    System.out.print( '\t' );
                }
                System.out.println();
            }
        }
        catch ( final FileNotFoundException e ) {
            e.printStackTrace();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }
}
