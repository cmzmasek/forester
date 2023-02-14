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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.stream.Collectors;

import org.forester.evoinference.distance.PairwiseDistanceCalculator;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.DeleteableMsa;

public final class msa2pwd_table {

    public static void main( final String args[] ) {
        try {
            final File indir = new File( args[ 0 ] );
            final String target = args[ 1 ];
            final File[] list_of_files = indir.listFiles();
            Arrays.sort( list_of_files );
            for( final File msa_file : list_of_files ) {
                if ( msa_file.isFile() && msa_file.toString().endsWith( ".fasta" )
                        && ( msa_file.toString().toLowerCase().indexOf( "mafft" ) > 0 ) ) {
                    doit( target, msa_file );
                }
            }
        }
        catch ( final FileNotFoundException e ) {
            e.printStackTrace();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }

    private static void doit( final String target, final File msa_file ) throws FileNotFoundException, IOException {
        System.out.print( msa_file.getName().substring( 0, msa_file.getName().toLowerCase().indexOf( "mafft" ) - 1 ) );
        //final File outfile = new File( args[ 1 ] );
        DeleteableMsa msa = null;
        final FileInputStream is = new FileInputStream( msa_file );
        if ( FastaParser.isLikelyFasta( msa_file ) ) {
            msa = DeleteableMsa.createInstance( FastaParser.parseMsa( is ) );
        }
        else {
            msa = DeleteableMsa.createInstance( GeneralMsaParser.parseMsa( is ) );
        }
        final BasicSymmetricalDistanceMatrix pwd = PairwiseDistanceCalculator.calcFractionalSimilarities( msa );
        int target_index = -1;
        boolean first = true;
        for( int i = 0; i < pwd.getSize(); i++ ) {
            if ( pwd.getIdentifier( i ).indexOf( target ) > -1 ) {
                target_index = i;
                final Map<String, Double> id_to_score = new HashMap<>();
                for( int col = 0; col < pwd.getSize(); ++col ) {
                    if ( col != target_index ) {
                        final double val = pwd.getValue( col, target_index );
                        id_to_score.put( pwd.getIdentifier( col ), val );
                    }
                }
                final Map<String, Double> sorted_id_to_score = id_to_score.entrySet().stream()
                        .sorted( Map.Entry.comparingByValue( Comparator.reverseOrder() ) )
                        .collect( Collectors.toMap( Map.Entry::getKey,
                                                    Map.Entry::getValue,
                                                    ( e1, e2 ) -> e1,
                                                    LinkedHashMap::new ) );
                if ( first ) {
                    first = false;
                }
                System.out.print( "\t" );
                System.out.print( pwd.getIdentifier( target_index ) );
                boolean first_inner = true;
                for( final Map.Entry<String, Double> entry : sorted_id_to_score.entrySet() ) {
                    if ( first_inner ) {
                        first_inner = false;
                        System.out.print( "\t" );
                    }
                    else {
                        System.out.print( "\t" );
                        System.out.print( "\t" );
                    }
                    System.out.println( entry.getKey() + "\t" + entry.getValue() );
                }
            }
        }
        if ( target_index > -1 ) {
        }
        else {
            //dosomething
        }
        // BufferedWriter w = ForesterUtil.createBufferedWriter( outfile );
        // pwd.write( w);
        //w.close();
    }
}
