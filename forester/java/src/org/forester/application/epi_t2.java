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

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.msa.DeleteableMsa;
import org.forester.msa.MsaMethods;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.ForesterUtil;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class epi_t2 {

    private final static String VERSION = "1.0.1";
    private final static Pattern GAP_C_TERM_ONLY = Pattern.compile("[^\\-]+\\-+");
    private final static Pattern GAP_N_TERM_ONLY = Pattern.compile("\\-+[^\\-]+");
    private final static Pattern GAP_ONLY = Pattern.compile("\\-+");
    private final static boolean TEST_1 = false;
    private final static boolean TEST_2 = false;

    public static void main(final String[] args) {
        if (args.length != 5) {
            System.err.println("Usage  : epi_t2 <h|s> <k|f> <add to tolerance> <msa> <peptides file> (h: heatmap, s: sequences, k: keep gaps, f: fill-in gaps)");
            System.err.println("Example: epi_t2 s f 2 Mammarenavirus_L_protein_all_mafft.fasta lassa_all_proteins_L.tsv > Mammarenavirus_L_protein_all_SEQUENCES.tsv");
            System.exit(-1);
        }
        try {
            final String o_hs = args[0];
            boolean heatmap = false;
            if (o_hs.equals("h")) {
                heatmap = true;
            } else if (o_hs.equals("s")) {
                heatmap = false;
            } else {
                System.err.println("use 'h' for heatmap or 's' for sequences");
                System.exit(-1);
            }
            final String o_kcgm = args[1];
            boolean keep_complete_gap_matches = false;
            if (o_kcgm.equals("k")) {
                keep_complete_gap_matches = true;
            } else if (o_kcgm.equals("f")) {
                keep_complete_gap_matches = false;
            } else {
                System.err.println("use 'k' the keep completely gapped matches or 'f' to fill in with non-gap chars");
                System.exit(-1);
            }


            final String add_to_tolerance_str = args[2];
            final int add_to_tolerance = Integer.parseInt(add_to_tolerance_str);
            System.out.println("Version: " + VERSION);
            System.out.println("Add to tolerance: " + add_to_tolerance);
            final File msa_file = new File(args[3]);
            final File peptide_seqs_file = new File(args[4]);
            DeleteableMsa msa = null;
            final FileInputStream is = new FileInputStream(msa_file);
            if (FastaParser.isLikelyFasta(msa_file)) {
                msa = DeleteableMsa.createInstance(FastaParser.parseMsa(is));
            } else {
                msa = DeleteableMsa.createInstance(GeneralMsaParser.parseMsa(is));
            }
            final List<String> peptide_seqs = new ArrayList<>();
            final List<Integer> peptide_start = new ArrayList<>();
            final List<Integer> peptide_stop = new ArrayList<>();
            final List<Integer> peptide_length = new ArrayList<>();


            final BufferedReader reader = new BufferedReader(new FileReader(peptide_seqs_file));
            String line = reader.readLine();
            while (line != null) {
                if (line.length() > 0) {
                    if (line.endsWith("\t")) {
                        line = line + " ";
                    }
                    final String[] s = line.split("\t");
                    if (s.length != 4) {
                        System.err.println("error: unexpected format: " + line);
                        System.exit(-1);
                    }
                    peptide_start.add(Integer.valueOf(s[0]));
                    peptide_stop.add(Integer.valueOf(s[1]));
                    peptide_length.add(Integer.valueOf(s[2]));
                    peptide_seqs.add(s[3]);
                }
                line = reader.readLine();
            }
            reader.close();
            System.out.print("START");
            System.out.print("\t");
            System.out.print("STOP");
            System.out.print("\t");
            System.out.print("LENGTH");
            System.out.print("\t");
            System.out.print("EPITOPE");
            System.out.print("\t");
            System.out.print("FIRST (MSA)");
            System.out.print("\t");
            System.out.print("LAST (MSA)");
            System.out.print("\t");
            System.out.print("MEDIAN CNSV");
            System.out.print("\t");
            System.out.print("IQR CNSV");
            System.out.print("\t");
            System.out.print("MIN CNSV");
            System.out.print("\t");
            System.out.print("MAX CNSV");
            System.out.print("\t");
            System.out.print("SHANNON ENT");
            System.out.print("\t");
            for (int row = 0; row < msa.getNumberOfSequences(); ++row) {
                System.out.print(msa.getIdentifier(row));
                System.out.print("\t");
                if (!heatmap) {
                    System.out.print("");
                    System.out.print("\t");
                }
            }
            System.out.println();
            for (int p = 0; p < peptide_seqs.size(); ++p) {
                final String peptide_seq = peptide_seqs.get(p);
                boolean found = false;
                int first = -1;
                int last = -1;
                for (int row = 0; row < msa.getNumberOfSequences(); ++row) {
                    final String current_seq_str = msa.getSequenceAsString(row).toString();
                    final int i = current_seq_str.indexOf(peptide_seq);
                    if (i > -1) {
                        first = i;
                        last = (first + peptide_seq.length()) - 1;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    final int max_t = (peptide_seq.length() / 2) + add_to_tolerance;
                    T:
                    for (int t = 0; t < max_t; ++t) {
                        for (int row = 0; row < msa.getNumberOfSequences(); ++row) {
                            final String current_seq_str = msa.getSequenceAsString(row).toString();
                            final List<Object> match_result = match(peptide_seq, current_seq_str, t);
                            if (match_result != null) {
                                found = true;
                                if (TEST_1) {
                                    System.err.println(t + ")  " + peptide_seq + " -> " + match_result.get(0));
                                }
                                first = (int) match_result.get(1);
                                last = (int) match_result.get(2) - 1;

                                if ( first < peptide_start.get(p) ) {
                                    System.err.println("WARNING: PROBLEM WITH: " + peptide_seq + " (\"add to tolerance\" likely too high)");
                                }


                                break T;
                            }
                        }
                    }
                }
                if (!found) {
                    System.err.println("WARNING: NOT FOUND: " + peptide_seq);
                }
                if (found) {
                    System.out.print(peptide_start.get(p));
                    System.out.print("\t");
                    System.out.print(peptide_stop.get(p));
                    System.out.print("\t");
                    System.out.print(peptide_length.get(p));
                    System.out.print("\t");
                    System.out.print(peptide_seq);
                    System.out.print("\t");

                    final BasicDescriptiveStatistics stats = new BasicDescriptiveStatistics();
                    final StringBuilder data_sb = new StringBuilder();
                    for (int row = 0; row < msa.getNumberOfSequences(); ++row) {
                        final String current_seq_str = msa.getSequenceAsString(row).toString();
                        String positional_homolog = current_seq_str.substring(first, last + 1);
                        final String orig_positional_homolog = positional_homolog;
                        if (positional_homolog.indexOf('-') > -1) {
                            final Matcher ma_n = GAP_N_TERM_ONLY.matcher(positional_homolog);
                            final Matcher ma_c = GAP_C_TERM_ONLY.matcher(positional_homolog);
                            final Matcher ma_go = GAP_ONLY.matcher(positional_homolog);
                            final int orig_length = positional_homolog.length();
                            if (TEST_2) {
                                System.err
                                        .println("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                                System.err.println(">>>>>> " + positional_homolog);
                            }
                            if (ma_n.matches()) {
                                if (TEST_2) {
                                    System.err.println("N-TERM GAP");
                                }
                                String new_positional_homolog = positional_homolog.replace("-", "");
                                int e = 1;
                                while (((new_positional_homolog.indexOf('-') > -1)
                                        || (new_positional_homolog.length() < orig_length))
                                        && ((first - e) > 0)) {
                                    new_positional_homolog = current_seq_str.substring(first - e++, last + 1);
                                    new_positional_homolog = new_positional_homolog.replace("-", "");
                                }
                                if (TEST_2) {
                                    System.err.println("--> " + new_positional_homolog);
                                }
                                positional_homolog = new_positional_homolog;
                            } else if (ma_c.matches()) {
                                if (TEST_2) {
                                    System.err.println("C-TERM GAP");
                                }
                                String new_positional_homolog = positional_homolog.replace("-", "");
                                int e = 1;
                                while (((new_positional_homolog.indexOf('-') > -1)
                                        || (new_positional_homolog.length() < orig_length))
                                        && ((last + 1 + e) < (current_seq_str.length()))) {
                                    new_positional_homolog = current_seq_str.substring(first, last + 1 + e++);
                                    new_positional_homolog = new_positional_homolog.replace("-", "");
                                }
                                if (TEST_2) {
                                    System.err.println("--> " + new_positional_homolog);
                                }
                                positional_homolog = new_positional_homolog;
                            }
                            if ((!(keep_complete_gap_matches && ma_go.matches())
                                    && (!ma_c.matches() && !ma_c.matches()))
                                    || (positional_homolog.length() < orig_length)) {
                                boolean done_n = false;
                                boolean done_c = false;
                                String new_positional_homolog = positional_homolog.replace("-", "");
                                int e_n = 0;
                                int e_c = 0;
                                while (!done_n || !done_c) {
                                    if ((new_positional_homolog.indexOf('-') < 0)
                                            && (new_positional_homolog.length() >= orig_length)) {
                                        break;
                                    }
                                    if ((first - e_n) > 0) {
                                        e_n++;
                                        new_positional_homolog = current_seq_str.substring(first - e_n,
                                                last + 1 + e_c);
                                        new_positional_homolog = new_positional_homolog.replace("-", "");
                                    } else {
                                        done_n = true;
                                    }
                                    if ((new_positional_homolog.indexOf('-') < 0)
                                            && (new_positional_homolog.length() >= orig_length)) {
                                        break;
                                    }
                                    if ((last + 1 + e_c) < (current_seq_str.length())) {
                                        e_c++;
                                        new_positional_homolog = current_seq_str.substring(first - e_n,
                                                last + 1 + e_c);
                                        new_positional_homolog = new_positional_homolog.replace("-", "");
                                    } else {
                                        done_c = true;
                                    }
                                }
                                if (TEST_2) {
                                    System.err.println("--> " + new_positional_homolog);
                                }
                                positional_homolog = new_positional_homolog;
                            }
                            if (TEST_2) {
                                System.err.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                            }
                        }
                        if (!heatmap) {
                            data_sb.append(positional_homolog);
                            data_sb.append("\t");
                        }
                        final double sim = calcSimilarity(peptide_seq, orig_positional_homolog);
                        stats.addValue(sim);

                        data_sb.append(ForesterUtil.round(sim, 4));
                        data_sb.append("\t");
                    }
                    System.out.print(first);
                    System.out.print("\t");
                    System.out.print(last);
                    System.out.print("\t");
                    System.out.print(ForesterUtil.round(stats.median(), 4));
                    System.out.print("\t");
                    System.out.print(ForesterUtil.round(stats.interquartileRange(), 4));
                    System.out.print("\t");
                    System.out.print(ForesterUtil.round(stats.getMin(), 4));
                    System.out.print("\t");
                    System.out.print(ForesterUtil.round(stats.getMax(), 4));
                    System.out.print("\t");

                    System.out.print(ForesterUtil
                            .round(MsaMethods.calcAvgNormalizedShannonsEntropy(21, msa, first, last), 4));
                    System.out.print("\t");
                    System.out.print(data_sb);
                    System.out.println();
                }
            }
        } catch (final FileNotFoundException e) {
            e.printStackTrace();
        } catch (final IOException e) {
            e.printStackTrace();
        }
    }

    private static List<Object> match(final String query, final String target, final int tolerance) {
        final int target_index_max = (target.length() - query.length()) + 1;
        for (int target_index = 0; target_index < target_index_max; target_index++) {
            int missed = 0;
            for (int query_index = 0; query_index < query.length(); ++query_index) {
                final char target_char = target.charAt(target_index + query_index);
                final char query_char = query.charAt(query_index);
                if (target_char != query_char) {
                    missed++;
                }
                if (missed > tolerance) {
                    break;
                }
            }
            if (missed <= tolerance) {
                return Arrays.asList(target.substring(target_index, target_index + query.length()),
                        target_index,
                        target_index + query.length());
            }
        }
        return null;
    }

    static double calcSimilarity(final String s1, final String s2) {
        final int l = s1.length();
        if (l != s2.length()) {
            throw new IllegalArgumentException("Unequal length: " + s1 + ", " + s2);
        }
        int s = 0;
        for (int i = 0; i < l; ++i) {
            if (s1.charAt(i) == s2.charAt(i)) {
                ++s;
            }
        }
        return ((double) s) / l;
    }
}
