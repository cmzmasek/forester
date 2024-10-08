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
    // add to tolerance: Spike: 12, PP: 3, N: 4, ORF7a: 2, ORF8: 3

    private final static String  VERSION         = "1.0.1";
    private final static Pattern GAP_C_TERM_ONLY = Pattern.compile( "[^\\-]+\\-+" );
    private final static Pattern GAP_N_TERM_ONLY = Pattern.compile( "\\-+[^\\-]+" );
    private final static Pattern GAP_ONLY        = Pattern.compile( "\\-+" );
    private final static String  SARBECO_TAG     = "|Sarbeco";
    private final static boolean TEST_1          = false;
    private final static boolean TEST_2          = false;
    public static void main( final String args[] ) {
        if ( args.length != 7 ) {
            System.err.println( "Usage: epi_t <h|s> <k|f> <i|n> <s|n> <add to tolerance> <msa> <peptides file>" );
            System.err.println( "Usage: epi_t h f n s 4 B_Membrane_BLAST_results_mafft.fasta Membrane.txt > m_01.tsv" );
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
            final String o_s = args[ 3 ];
            boolean separate_sarbeco_vs_non_sarbeco_stats = false;
            if ( o_s.equals( "s" ) ) {
                separate_sarbeco_vs_non_sarbeco_stats = true;
            }
            else if ( o_s.equals( "n" ) ) {
                separate_sarbeco_vs_non_sarbeco_stats = false;
            }
            else {
                System.err.println( "use 's' for separate Sarbeco vs non-Sarbeco statistics, 'n' otherwise" );
                System.exit( -1 );
            }
            final String add_to_tolerance_str = args[ 4 ];
            final int add_to_tolerance = Integer.parseInt( add_to_tolerance_str );
            System.out.println( "Version: " + VERSION );
            System.out.println( "Add to tolerance: " + add_to_tolerance );
            final File msa_file = new File( args[ 5 ] );
            final File peptide_seqs_file = new File( args[ 6 ] );
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
                    if ( line.endsWith( "\t" ) ) {
                        line = line + " ";
                    }
                    final String[] s = line.split( "\t" );
                    if ( s.length < 8 ) {
                        System.err.println( "error: unexpected format: " + line );
                        System.exit( -1 );
                    }
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
            System.out.print( "MEDIAN CNSV" );
            System.out.print( "\t" );
            System.out.print( "IQR CNSV" );
            System.out.print( "\t" );
            System.out.print( "MIN CNSV" );
            System.out.print( "\t" );
            System.out.print( "MAX CNSV" );
            System.out.print( "\t" );
            if ( separate_sarbeco_vs_non_sarbeco_stats ) {
                System.out.print( "MEDIAN CNSV SARB" );
                System.out.print( "\t" );
                System.out.print( "IQR CNSV" );
                System.out.print( "\t" );
                System.out.print( "MIN CNSV" );
                System.out.print( "\t" );
                System.out.print( "MAX CNSV" );
                System.out.print( "\t" );
                System.out.print( "MEDIAN CNSV NON-SARB" );
                System.out.print( "\t" );
                System.out.print( "IQR CNSV" );
                System.out.print( "\t" );
                System.out.print( "MIN CNSV" );
                System.out.print( "\t" );
                System.out.print( "MAX CNSV" );
                System.out.print( "\t" );
            }
            System.out.print( "SHANNON ENT" );
            System.out.print( "\t" );
            for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                System.out.print( msa.getIdentifier( row ) );
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
                    final int max_t = ( peptide_seq.length() / 2 ) + add_to_tolerance;
                    T: for( int t = 0; t < max_t; ++t ) {
                        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                            final String current_seq_str = msa.getSequenceAsString( row ).toString();
                            final List<Object> match_result = match( peptide_seq, current_seq_str, t );
                            if ( match_result != null ) {
                                found = true;
                                if ( TEST_1 ) {
                                    System.err.println( t + ")  " + peptide_seq + " -> " + match_result.get( 0 ) );
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
                    System.err.println( "WARNING: NOT FOUND: " + peptide_seq );
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
                    final BasicDescriptiveStatistics stats_non_sarbeco = new BasicDescriptiveStatistics();
                    final BasicDescriptiveStatistics stats_sarbeco = new BasicDescriptiveStatistics();
                    final StringBuilder data_sb = new StringBuilder();
                    for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
                        final String current_seq_str = msa.getSequenceAsString( row ).toString();
                        final String current_seq_name = msa.getIdentifier( row );
                        String positional_homolog = current_seq_str.substring( first, last + 1 );
                        final String orig_positional_homolog = positional_homolog;
                        if ( positional_homolog.indexOf( '-' ) > -1 ) {
                            final Matcher ma_n = GAP_N_TERM_ONLY.matcher( positional_homolog );
                            final Matcher ma_c = GAP_C_TERM_ONLY.matcher( positional_homolog );
                            final Matcher ma_go = GAP_ONLY.matcher( positional_homolog );
                            final int orig_length = positional_homolog.length();
                            if ( TEST_2 ) {
                                System.err
                                        .println( "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" );
                                System.err.println( ">>>>>> " + positional_homolog );
                            }
                            if ( ma_n.matches() ) {
                                if ( TEST_2 ) {
                                    System.err.println( "N-TERM GAP" );
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
                                    System.err.println( "--> " + new_positional_homolog );
                                }
                                positional_homolog = new_positional_homolog;
                            }
                            else if ( ma_c.matches() ) {
                                if ( TEST_2 ) {
                                    System.err.println( "C-TERM GAP" );
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
                                    System.err.println( "--> " + new_positional_homolog );
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
                                    System.err.println( "--> " + new_positional_homolog );
                                }
                                positional_homolog = new_positional_homolog;
                            }
                            if ( TEST_2 ) {
                                System.err.println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" );
                            }
                        }
                        if ( !heatmap ) {
                            data_sb.append( positional_homolog );
                            data_sb.append( "\t" );
                        }
                        final double sim = calcSimilarity( peptide_seq, orig_positional_homolog );
                        stats.addValue( sim );
                        if ( separate_sarbeco_vs_non_sarbeco_stats ) {
                            if ( current_seq_name.indexOf( SARBECO_TAG ) > -1 ) {
                                stats_sarbeco.addValue( sim );
                            }
                            else {
                                stats_non_sarbeco.addValue( sim );
                            }
                        }
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
                    if ( separate_sarbeco_vs_non_sarbeco_stats ) {
                        if ( ( stats_sarbeco.getN() > 0 ) && ( stats_non_sarbeco.getN() > 0 ) ) {
                            System.out.print( ForesterUtil.round( stats_sarbeco.median(), 4 ) );
                            System.out.print( "\t" );
                            System.out.print( ForesterUtil.round( stats_sarbeco.interquartileRange(), 4 ) );
                            System.out.print( "\t" );
                            System.out.print( ForesterUtil.round( stats_sarbeco.getMin(), 4 ) );
                            System.out.print( "\t" );
                            System.out.print( ForesterUtil.round( stats_sarbeco.getMax(), 4 ) );
                            System.out.print( "\t" );
                            System.out.print( ForesterUtil.round( stats_non_sarbeco.median(), 4 ) );
                            System.out.print( "\t" );
                            System.out.print( ForesterUtil.round( stats_non_sarbeco.interquartileRange(), 4 ) );
                            System.out.print( "\t" );
                            System.out.print( ForesterUtil.round( stats_non_sarbeco.getMin(), 4 ) );
                            System.out.print( "\t" );
                            System.out.print( ForesterUtil.round( stats_non_sarbeco.getMax(), 4 ) );
                            System.out.print( "\t" );
                        }
                        else {
                            System.out.print( "" );
                            System.out.print( "\t" );
                            System.out.print( "" );
                            System.out.print( "\t" );
                            System.out.print( "" );
                            System.out.print( "\t" );
                            System.out.print( "" );
                            System.out.print( "\t" );
                            System.out.print( "" );
                            System.out.print( "\t" );
                            System.out.print( "" );
                            System.out.print( "\t" );
                            System.out.print( "" );
                            System.out.print( "\t" );
                            System.out.print( "" );
                            System.out.print( "\t" );
                        }
                    }
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

    /*
     *
     * Positions for Polyprotein NSPs for MSA "Beta_PP_02_mafft" 2022/12/17
     *
     */
    private final static String inferNSPname( final String peptide_seq, final int i ) {
        String inferred_name = "na";
        int ii = i + 1;
        final int offset = peptide_seq.length() / 2;
        if ( ii > offset ) {
            ii -= offset;
        }
        if ( ( ii >= 1 ) && ( ii <= 269 ) ) {
            inferred_name = "NSP1 (Host translation inhibitor)";
        }
        else if ( ( ii >= 270 ) && ( ii <= 985 ) ) {
            inferred_name = "NSP2";
        }
        else if ( ( ii >= 986 ) && ( ii <= 3625 ) ) {
            inferred_name = "NSP3 (Papain-like protease)";
        }
        else if ( ( ii >= 3626 ) && ( ii <= 4138 ) ) {
            inferred_name = "NSP4";
        }
        else if ( ( ii >= 4139 ) && ( ii <= 4466 ) ) {
            inferred_name = "NSP5 (3C-like proteinase)";
        }
        else if ( ( ii >= 4467 ) && ( ii <= 4769 ) ) {
            inferred_name = "NSP6";
        }
        else if ( ( ii >= 4770 ) && ( ii <= 4858 ) ) {
            inferred_name = "NSP7";
        }
        else if ( ( ii >= 4859 ) && ( ii <= 5059 ) ) {
            inferred_name = "NSP8";
        }
        else if ( ( ii >= 5060 ) && ( ii <= 5172 ) ) {
            inferred_name = "NSP9 (RNA-capping enzyme subunit)";
        }
        else if ( ( ii >= 5173 ) && ( ii <= 5312 ) ) {
            inferred_name = "NSP10";
        }
        else if ( ( ii >= 5313 ) && ( ii <= 5344 ) ) {
            inferred_name = "NSP11";
        }
        else if ( ( ii >= 5345 ) && ( ii <= 6252 ) ) {
            inferred_name = "NSP12 (RNA-directed RNA polymerase)";
        }
        else if ( ( ii >= 6253 ) && ( ii <= 6858 ) ) {
            inferred_name = "NSP13 (Helicase)";
        }
        else if ( ( ii >= 6859 ) && ( ii <= 7390 ) ) {
            inferred_name = "NSP14 (Guanine-N7 methyltransferase)";
        }
        else if ( ( ii >= 7391 ) && ( ii <= 7781 ) ) {
            inferred_name = "NSP15 (Uridylate-specific endoribonuclease)";
        }
        else if ( ii >= 7782 ) {
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
