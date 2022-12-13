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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.DeleteableMsa;
import org.forester.util.ForesterUtil;

public final class epi_t {

    private final static Pattern GAP_C_TERM_ONLY = Pattern.compile( "[^\\-]+\\-+" );
    private final static Pattern GAP_N_TERM_ONLY = Pattern.compile( "\\-+[^\\-]+" );
    private final static boolean HEATMAP         = false;
    private final static boolean TEST_1          = false;
    public static void main( final String args[] ) {
        // System.out.println( x("a", "abc", 0) );
        // System.out.println( x("b", "abc", 0) );
        // System.out.println( x("c", "abc", 0) );
        // System.out.println( x("d", "abc", 0) );
        // System.out.println( x("ab", "abc", 0) );
        // System.out.println( x("bc", "abc", 0) );
        // System.out.println( x("abc", "abc", 0) );
        // System.out.println( x("d", "abc", 0) );
        // System.out.println( x("ef_hi", "abcdefghijklmnopqrstuvw", 0) );
        // System.out.println( x("ef_hi", "abcdefghijklmnopqrstuvw", 1) );
        //
        // System.out.println( x("ef__i", "abcdefghijklmnopqrstuvw", 0) );
        // System.out.println( x("ef__i", "abcdefghijklmnopqrstuvw", 1) );
        // System.out.println( x("ef__i", "abcdefghijklmnopqrstuvw", 2) );
        // System.out.println( x("_f__i", "abcdefghijklmnopqrstuvw", 3) );
        try {
            final File msa_file = new File( args[ 0 ] );
            final File peptide_seqs_file = new File( args[ 1 ] );
            final File outfile = new File( args[ 2 ] );
            DeleteableMsa msa = null;
            final FileInputStream is = new FileInputStream( msa_file );
            if ( FastaParser.isLikelyFasta( msa_file ) ) {
                msa = DeleteableMsa.createInstance( FastaParser.parseMsa( is ) );
            }
            else {
                msa = DeleteableMsa.createInstance( GeneralMsaParser.parseMsa( is ) );
            }
            //
            final List<String> peptide_seqs = new ArrayList<>();
            final List<String> peptide_ids = new ArrayList<>();
            final List<String> taxonomy = new ArrayList<>();
            final List<String> cd = new ArrayList<>();
            final List<String> megapool = new ArrayList<>();
            final List<String> start = new ArrayList<>();
            final List<String> num_times_tested = new ArrayList<>();
            final List<String> num_times_positive = new ArrayList<>();
            final List<String> total_magnitude = new ArrayList<>();
            final List<String> genome_accs = new ArrayList<>();
            final BufferedReader reader = new BufferedReader( new FileReader( peptide_seqs_file ) );
            String line = reader.readLine();
            while ( line != null ) {
                line = line.trim();
                if ( line.length() > 0 ) {
                    final String[] s = line.split( "\t" );
                    taxonomy.add( s[ 0 ] );
                    cd.add( s[ 1 ] );
                    peptide_seqs.add( s[ 2 ] );
                    //megapool.add( s[ 1 ] );
                    //start.add( s[ 2 ] );
                    // peptide_seqs.add( s[ 3 ] );
                    // peptide_ids.add( s[ 4 ] );
                    // num_times_tested.add( s[ 5 ] );
                    //  num_times_positive.add( s[ 6 ] );
                    // total_magnitude.add( s[ 7 ] );
                    //genome_accs.add( s[ 8 ] );
                }
                line = reader.readLine();
            }
            reader.close();
            System.out.print( "TAXO" );
            System.out.print( "\t" );
            System.out.print( "CD" );
            System.out.print( "\t" );
            System.out.print( "SEQ" );
            System.out.print( "\t" );
            //
            //            System.out.print( "#CD4+" );
            //            System.out.print( "\t" );
            //            System.out.print( "Megapool" );
            //            System.out.print( "\t" );
            //            System.out.print( "Start" );
            //            System.out.print( "\t" );
            //            System.out.print( "Sequence" );
            //            System.out.print( "\t" );
            //            System.out.print( "Peptide_ID" );
            //            System.out.print( "\t" );
            //            System.out.print( "Num_times_tested" );
            //            System.out.print( "\t" );
            //            System.out.print( "Num_times_positiv" );
            //            System.out.print( "\t" );
            //            System.out.print( "Total_magnitude" );
            //            System.out.print( "\t" );
            //            System.out.print( "Genome_accession" );
            //            System.out.print( "\t" );
            //            System.out.print( "Match first" );
            //            System.out.print( "\t" );
            //            System.out.print( "Match last" );
            //            System.out.print( "\t" );
            for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                final String current_seq_name = msa.getIdentifier( row );
                System.out.print( current_seq_name );
                System.out.print( "\t" );
                if ( !HEATMAP ) {
                    System.out.print( "" );
                    System.out.print( "\t" );
                }
            }
            System.out.println();
            for( int p = 0; p < peptide_seqs.size(); ++p ) {
                final String peptide_seq = peptide_seqs.get( p );
                boolean found = false;
                int first = -1;
                int last = -1;
                for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                    final String current_seq_str = msa.getSequenceAsString( row ).toString();
                    final int i = current_seq_str.indexOf( peptide_seq );
                    if ( i > -1 ) {
                        first = i;
                        last = ( first + peptide_seq.length() ) - 1;
                        found = true;
                        break;
                    }
                }
                if ( !found ) {
                    final int max_t = ( peptide_seq.length() / 2 ) - 1;
                    T: for( int t = 0; t < max_t; ++t ) {
                        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                            final String current_seq_str = msa.getSequenceAsString( row ).toString();
                            final List<Object> match_result = match( peptide_seq, current_seq_str, t );
                            if ( match_result != null ) {
                                found = true;
                                if ( TEST_1 ) {
                                    System.out.println( t + ")  " + peptide_seq + " -> " + match_result.get( 0 ) );
                                }
                                first = ( int ) match_result.get( 1 );
                                last = ( int ) match_result.get( 2 ) - 1;
                                break T;
                            }
                        }
                    }
                }
                if ( !found ) {
                    if ( TEST_1 ) {
                        System.out.println( "STILL NOT FOUND: " + peptide_seq );
                    }
                }
                if ( found ) {
                    System.out.print( taxonomy.get( p ) );
                    System.out.print( "\t" );
                    System.out.print( cd.get( p ) );
                    System.out.print( "\t" );
                    System.out.print( peptide_seq );
                    System.out.print( "\t" );
                    //                    System.out.print( megapool.get( p ) );
                    //                    System.out.print( "\t" );
                    //                    System.out.print( start.get( p ) );
                    //                    System.out.print( "\t" );
                    //                    System.out.print( peptide_seq );
                    //                    System.out.print( "\t" );
                    //                    System.out.print( peptide_id );
                    //                    System.out.print( "\t" );
                    //                    System.out.print( num_times_tested.get( p ) );
                    //                    System.out.print( "\t" );
                    //                    System.out.print( num_times_positive.get( p ) );
                    //                    System.out.print( "\t" );
                    //                    System.out.print( total_magnitude.get( p ) );
                    //                    System.out.print( "\t" );
                    //                    System.out.print( genome_accs.get( p ) );
                    //                    System.out.print( "\t" );
                    //                    System.out.print( first );
                    //                    System.out.print( "\t" );
                    //                    System.out.print( last );
                    //System.out.print( "\t" );
                    for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                        final String current_seq_str = msa.getSequenceAsString( row ).toString();
                        //final String current_seq_name = msa.getIdentifier( row );
                        //System.out.print( current_seq_name );
                        //System.out.print( "\t" );
                        String positional_homolog = current_seq_str.substring( first, last + 1 );
                        if ( positional_homolog.indexOf( '-' ) > -1 ) {
                            final Matcher ma_n = GAP_N_TERM_ONLY.matcher( positional_homolog );
                            final Matcher ma_c = GAP_C_TERM_ONLY.matcher( positional_homolog );
                            // System.out.println( "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" );
                            // System.out.println( ">>>>>> " + positional_homolog );
                            if ( ma_n.matches() ) {
                                // System.out.println( "N-TERM GAP" );
                                String new_positional_homolog = positional_homolog.replace( "-", "" );
                                int e = 1;
                                final int orig_length = positional_homolog.length();
                                while ( ( ( new_positional_homolog.indexOf( '-' ) > -1 )
                                        || ( new_positional_homolog.length() < orig_length ) )
                                        && ( ( first - e ) > 0 ) ) {
                                    new_positional_homolog = current_seq_str.substring( first - e++, last + 1 );
                                    new_positional_homolog = new_positional_homolog.replace( "-", "" );
                                }
                                // System.out.println( "--> " + new_positional_homolog );
                                positional_homolog = new_positional_homolog;
                            }
                            else if ( ma_c.matches() ) {
                                //  System.out.println( "C-TERM GAP" );
                                String new_positional_homolog = positional_homolog.replace( "-", "" );
                                int e = 1;
                                final int orig_length = positional_homolog.length();
                                while ( ( ( new_positional_homolog.indexOf( '-' ) > -1 )
                                        || ( new_positional_homolog.length() < orig_length ) )
                                        && ( ( last + 1 + e ) < ( current_seq_str.length() - 0 ) ) ) {
                                    new_positional_homolog = current_seq_str.substring( first, last + 1 + e++ );
                                    new_positional_homolog = new_positional_homolog.replace( "-", "" );
                                }
                                // System.out.println( "--> " + new_positional_homolog );
                                positional_homolog = new_positional_homolog;
                            }
                            else {
                                //System.out.println( "GAPS EVERYWHERE!!!" );
                                boolean done_n = false;
                                boolean done_c = false;
                                String new_positional_homolog = positional_homolog.replace( "-", "" );
                                int e_n = 0;
                                int e_c = 0;
                                final int orig_length = positional_homolog.length();
                                while ( !done_n || !done_c ) {
                                    if ( ( new_positional_homolog.indexOf( '-' ) < 0 )
                                            && ( new_positional_homolog.length() >= orig_length ) ) {
                                        break;
                                    }
                                    if ( ( first - e_n ) > 0 ) {
                                        e_n++;
                                        new_positional_homolog = current_seq_str.substring( first - e_n,
                                                                                            last + 1 + e_c );
                                        new_positional_homolog = new_positional_homolog.replace( "-", "" );
                                    }
                                    else {
                                        done_n = true;
                                    }
                                    if ( ( new_positional_homolog.indexOf( '-' ) < 0 )
                                            && ( new_positional_homolog.length() >= orig_length ) ) {
                                        break;
                                    }
                                    if ( ( last + 1 + e_c ) < ( current_seq_str.length() - 0 ) ) {
                                        e_c++;
                                        new_positional_homolog = current_seq_str.substring( first - e_n,
                                                                                            last + 1 + e_c );
                                        new_positional_homolog = new_positional_homolog.replace( "-", "" );
                                    }
                                    else {
                                        done_c = true;
                                    }
                                }
                                // System.out.println( "--> " + new_positional_homolog );
                                positional_homolog = new_positional_homolog;
                            }
                            // System.out.println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" );
                        }
                        if ( !HEATMAP ) {
                            System.out.print( positional_homolog );
                            System.out.print( "\t" );
                        }
                        if ( peptide_seq.length() != positional_homolog.length() ) {
                            // System.out.println( "      !! WARNING " + peptide_seq + " != " + positional_homolog );
                            // System.out.println( "                 " + taxonomy.get( p ) );
                        }
                        final double diss = calcDissimilarity( peptide_seq, positional_homolog );
                        System.out.print( ForesterUtil.round( diss, 4 ) );
                        System.out.print( "\t" );
                    }
                    System.out.println();
                }
            }
            //
        }
        catch ( final FileNotFoundException e ) {
            e.printStackTrace();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }

    private final static List<Object> match( final String query, final String target, final int tolerance ) {
        final int target_index_max = ( target.length() - query.length() ) + 1;
        for( int target_index = 0; target_index < target_index_max; target_index++ ) {
            int missed = 0;
            for( int query_index = 0; query_index < query.length(); ++query_index ) {
                final char target_char = target.charAt( target_index + query_index );
                final char query_char = query.charAt( query_index );
                if ( target_char != query_char ) {
                    missed++;
                }
                if ( missed > tolerance ) {
                    break;
                }
            }
            if ( missed <= tolerance ) {
                return Arrays.asList( target.substring( target_index, target_index + query.length() ),
                                      target_index,
                                      target_index + query.length() );
            }
        }
        return null;
    }

    final static double calcDissimilarity( final String s1, final String s2 ) {
        final int l1 = s1.length();
        final int l2 = s2.length();
        int shorter = 0;
        int longer = 0;
        if ( l1 > l2 ) {
            longer = l1;
            shorter = l2;
        }
        else {
            longer = l2;
            shorter = l1;
        }
        int d = 0;
        for( int i = 0; i < shorter; ++i ) {
            if ( s1.charAt( i ) != s2.charAt( i ) ) {
                ++d;
            }
        }
        return ( ( double ) d ) / longer;
    }
}
