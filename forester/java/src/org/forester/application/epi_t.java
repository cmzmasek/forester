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
import org.forester.msa.MsaMethods;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.ForesterUtil;

public final class epi_t {

    private final static String  VERSION         = "1.0.0";
    private final static Pattern GAP_C_TERM_ONLY = Pattern.compile( "[^\\-]+\\-+" );
    private final static Pattern GAP_N_TERM_ONLY = Pattern.compile( "\\-+[^\\-]+" );
    private final static Pattern GAP_ONLY        = Pattern.compile( "\\-+" );
    private final static boolean TEST_1          = false;
    private final static boolean TEST_2          = false;
    public static void main( final String args[] ) {
        if ( args.length != 5 ) {
            System.err.println( "Usage: epi_t <'h' or 's'> <'k' or 'f'> <'i' or 'n'> <ms> <peptides file>" );
            System.exit( -1 );
        }
        try {
            final String o_hs = args[ 0 ];
            boolean heatmap = false;
            if ( o_hs.equals( "h" ) ) {
                heatmap = true;
            }
            else if ( o_hs.equals( "s" ) ) {
                heatmap = false;
            }
            else {
                System.err.println( "use 'h' for heatmap or 's' for sequences" );
                System.exit( -1 );
            }
            final String o_kcgm = args[ 1 ];
            boolean keep_complete_gap_matches = false;
            if ( o_kcgm.equals( "k" ) ) {
                keep_complete_gap_matches = true;
            }
            else if ( o_kcgm.equals( "f" ) ) {
                keep_complete_gap_matches = false;
            }
            else {
                System.err.println( "use 'k' the keep completely gapped matches or 'f' to fill in with non-gap chars" );
                System.exit( -1 );
            }
            final String o_mpi = args[ 2 ];
            boolean mat_pep_inference = false;
            if ( o_mpi.equals( "i" ) ) {
                mat_pep_inference = true;
            }
            else if ( o_mpi.equals( "n" ) ) {
                mat_pep_inference = false;
            }
            else {
                System.err.println( "use 'i' to infer mature peptides or 'n' to not infer" );
                System.exit( -1 );
            }
            System.out.println( "Version:" + VERSION );
            final File msa_file = new File( args[ 3 ] );
            final File peptide_seqs_file = new File( args[ 4 ] );
            DeleteableMsa msa = null;
            final FileInputStream is = new FileInputStream( msa_file );
            if ( FastaParser.isLikelyFasta( msa_file ) ) {
                msa = DeleteableMsa.createInstance( FastaParser.parseMsa( is ) );
            }
            else {
                msa = DeleteableMsa.createInstance( GeneralMsaParser.parseMsa( is ) );
            }
            final List<String> peptide_seqs = new ArrayList<>();
            final List<String> taxonomy = new ArrayList<>();
            final List<String> pop = new ArrayList<>();
            final List<String> source = new ArrayList<>();
            final List<String> epitope_names = new ArrayList<>();
            final List<String> epitope_accs = new ArrayList<>();
            final BufferedReader reader = new BufferedReader( new FileReader( peptide_seqs_file ) );
            String line = reader.readLine();
            while ( line != null ) {
                if ( line.length() > 0 ) {
                    if ( line.startsWith( "\t" ) ) {
                        line = "NA" + line;
                    }
                    if ( line.endsWith( "\t" ) ) {
                        line = line + " ";
                    }
                    final String[] s = line.split( "\t" );
                    taxonomy.add( s[ 0 ] );
                    source.add( s[ 1 ] );
                    pop.add( s[ 2 ] );
                    peptide_seqs.add( s[ 3 ] );
                    epitope_names.add( s[ 6 ] );
                    epitope_accs.add( s[ 7 ] );
                }
                line = reader.readLine();
            }
            reader.close();
            System.out.print( "TAXO" );
            System.out.print( "\t" );
            System.out.print( "SOURCE" );
            System.out.print( "\t" );
            System.out.print( "POP" );
            System.out.print( "\t" );
            System.out.print( "EPITOPE" );
            System.out.print( "\t" );
            if ( mat_pep_inference ) {
                System.out.print( "NAME INFERRED" );
                System.out.print( "\t" );
            }
            System.out.print( "NAME" );
            System.out.print( "\t" );
            System.out.print( "ACC" );
            System.out.print( "\t" );
            System.out.print( "FIRST (MSA)" );
            System.out.print( "\t" );
            System.out.print( "LAST (MSA)" );
            System.out.print( "\t" );
            System.out.print( "MEDIAN CONS" );
            System.out.print( "\t" );
            System.out.print( "IQR" );
            System.out.print( "\t" );
            System.out.print( "MIN CONS" );
            System.out.print( "\t" );
            System.out.print( "MAX CONS" );
            System.out.print( "\t" );
            System.out.print( "SHANNON ENT" );
            System.out.print( "\t" );
            for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                final String current_seq_name = msa.getIdentifier( row );
                System.out.print( current_seq_name );
                System.out.print( "\t" );
                if ( !heatmap ) {
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
                String inferred_name = "na";
                for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                    final String current_seq_str = msa.getSequenceAsString( row ).toString();
                    final int i = current_seq_str.indexOf( peptide_seq );
                    if ( i > -1 ) {
                        first = i;
                        last = ( first + peptide_seq.length() ) - 1;
                        found = true;
                        if ( mat_pep_inference ) {
                            inferred_name = inferNSPname( peptide_seq, i );
                        }
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
                                if ( mat_pep_inference ) {
                                    inferred_name = inferNSPname( peptide_seq, first );
                                }
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
                    System.out.print( source.get( p ) );
                    System.out.print( "\t" );
                    System.out.print( pop.get( p ) );
                    System.out.print( "\t" );
                    System.out.print( peptide_seq );
                    System.out.print( "\t" );
                    if ( mat_pep_inference ) {
                        System.out.print( inferred_name );
                        System.out.print( "\t" );
                    }
                    System.out.print( epitope_names.get( p ) );
                    System.out.print( "\t" );
                    System.out.print( epitope_accs.get( p ) );
                    System.out.print( "\t" );
                    final BasicDescriptiveStatistics stats = new BasicDescriptiveStatistics();
                    final StringBuilder data_sb = new StringBuilder();
                    for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                        final String current_seq_str = msa.getSequenceAsString( row ).toString();
                        String positional_homolog = current_seq_str.substring( first, last + 1 );
                        final String orig_positional_homolog = positional_homolog;
                        if ( positional_homolog.indexOf( '-' ) > -1 ) {
                            final Matcher ma_n = GAP_N_TERM_ONLY.matcher( positional_homolog );
                            final Matcher ma_c = GAP_C_TERM_ONLY.matcher( positional_homolog );
                            final Matcher ma_go = GAP_ONLY.matcher( positional_homolog );
                            final int orig_length = positional_homolog.length();
                            if ( TEST_2 ) {
                                System.out
                                        .println( "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" );
                                System.out.println( ">>>>>> " + positional_homolog );
                            }
                            if ( ma_n.matches() ) {
                                if ( TEST_2 ) {
                                    System.out.println( "N-TERM GAP" );
                                }
                                String new_positional_homolog = positional_homolog.replace( "-", "" );
                                int e = 1;
                                while ( ( ( new_positional_homolog.indexOf( '-' ) > -1 )
                                        || ( new_positional_homolog.length() < orig_length ) )
                                        && ( ( first - e ) > 0 ) ) {
                                    new_positional_homolog = current_seq_str.substring( first - e++, last + 1 );
                                    new_positional_homolog = new_positional_homolog.replace( "-", "" );
                                }
                                if ( TEST_2 ) {
                                    System.out.println( "--> " + new_positional_homolog );
                                }
                                positional_homolog = new_positional_homolog;
                            }
                            else if ( ma_c.matches() ) {
                                if ( TEST_2 ) {
                                    System.out.println( "C-TERM GAP" );
                                }
                                String new_positional_homolog = positional_homolog.replace( "-", "" );
                                int e = 1;
                                while ( ( ( new_positional_homolog.indexOf( '-' ) > -1 )
                                        || ( new_positional_homolog.length() < orig_length ) )
                                        && ( ( last + 1 + e ) < ( current_seq_str.length() - 0 ) ) ) {
                                    new_positional_homolog = current_seq_str.substring( first, last + 1 + e++ );
                                    new_positional_homolog = new_positional_homolog.replace( "-", "" );
                                }
                                if ( TEST_2 ) {
                                    System.out.println( "--> " + new_positional_homolog );
                                }
                                positional_homolog = new_positional_homolog;
                            }
                            if ( ( !( keep_complete_gap_matches && ma_go.matches() )
                                    && ( !ma_c.matches() && !ma_c.matches() ) )
                                    || ( positional_homolog.length() < orig_length ) ) {
                                boolean done_n = false;
                                boolean done_c = false;
                                String new_positional_homolog = positional_homolog.replace( "-", "" );
                                int e_n = 0;
                                int e_c = 0;
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
                                if ( TEST_2 ) {
                                    System.out.println( "--> " + new_positional_homolog );
                                }
                                positional_homolog = new_positional_homolog;
                            }
                            if ( TEST_2 ) {
                                System.out.println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" );
                            }
                        }
                        if ( !heatmap ) {
                            data_sb.append( positional_homolog );
                            data_sb.append( "\t" );
                        }
                        final double sim = calcSimilarity( peptide_seq, orig_positional_homolog );
                        stats.addValue( sim );
                        data_sb.append( ForesterUtil.round( sim, 4 ) );
                        data_sb.append( "\t" );
                    }
                    System.out.print( first );
                    System.out.print( "\t" );
                    System.out.print( last );
                    System.out.print( "\t" );
                    System.out.print( ForesterUtil.round( stats.median(), 4 ) );
                    System.out.print( "\t" );
                    System.out.print( ForesterUtil.round( stats.interquartileRange(), 4 ) );
                    System.out.print( "\t" );
                    System.out.print( ForesterUtil.round( stats.getMin(), 4 ) );
                    System.out.print( "\t" );
                    System.out.print( ForesterUtil.round( stats.getMax(), 4 ) );
                    System.out.print( "\t" );
                    System.out.print( ForesterUtil
                            .round( MsaMethods.calcAvgNormalizedShannonsEntropy( 21, msa, first, last ), 4 ) );
                    System.out.print( "\t" );
                    System.out.print( data_sb.toString() );
                    System.out.println();
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

    private static String inferNSPname( final String peptide_seq, final int i ) {
        String inferred_name = "na";
        int ii = i + 1;
        final int offset = peptide_seq.length() / 2;
        if ( ii > offset ) {
            ii -= offset;
        }
        if ( ( ii >= 1 ) && ( ii <= 269 ) ) {
            inferred_name = "NSP1 (Host translation inhibitor)";
        }
        else if ( ( ii >= 270 ) && ( ii <= 1022 ) ) {
            inferred_name = "NSP2";
        }
        else if ( ( ii >= 1023 ) && ( ii <= 3651 ) ) {
            inferred_name = "NSP3 (Papain-like protease)";
        }
        else if ( ( ii >= 3652 ) && ( ii <= 4164 ) ) {
            inferred_name = "NSP4";
        }
        else if ( ( ii >= 4165 ) && ( ii <= 4492 ) ) {
            inferred_name = "NSP5 (3C-like proteinase)";
        }
        else if ( ( ii >= 4493 ) && ( ii <= 4795 ) ) {
            inferred_name = "NSP6";
        }
        else if ( ( ii >= 4796 ) && ( ii <= 4884 ) ) {
            inferred_name = "NSP7";
        }
        else if ( ( ii >= 4885 ) && ( ii <= 5085 ) ) {
            inferred_name = "NSP8";
        }
        else if ( ( ii >= 5086 ) && ( ii <= 5198 ) ) {
            inferred_name = "NSP9 (RNA-capping enzyme subunit)";
        }
        else if ( ( ii >= 5199 ) && ( ii <= 5338 ) ) {
            inferred_name = "NSP10";
        }
        else if ( ( ii >= 5339 ) && ( ii <= 5370 ) ) {
            inferred_name = "NSP11";
        }
        else if ( ( ii >= 5371 ) && ( ii <= 6275 ) ) {
            inferred_name = "NSP12 (RNA-directed RNA polymerase)";
        }
        else if ( ( ii >= 6276 ) && ( ii <= 6882 ) ) {
            inferred_name = "NSP13 (Helicase)";
        }
        else if ( ( ii >= 6883 ) && ( ii <= 7416 ) ) {
            inferred_name = "NSP14 (Guanine-N7 methyltransferase)";
        }
        else if ( ( ii >= 7417 ) && ( ii <= 7810 ) ) {
            inferred_name = "NSP15 (Uridylate-specific endoribonuclease)";
        }
        else if ( ii >= 7811 ) {
            inferred_name = "NSP16 (2'-O-methyltransferase)";
        }
        return inferred_name;
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

    final static double calcSimilarity( final String s1, final String s2 ) {
        final int l = s1.length();
        if ( l != s2.length() ) {
            throw new IllegalArgumentException( "Unequal length: " + s1 + ", " + s2 );
        }
        int s = 0;
        for( int i = 0; i < l; ++i ) {
            if ( s1.charAt( i ) == s2.charAt( i ) ) {
                ++s;
            }
        }
        return ( ( double ) s ) / l;
    }
}
