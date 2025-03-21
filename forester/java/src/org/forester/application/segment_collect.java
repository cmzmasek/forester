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
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.ForesterUtil;
import org.forester.util.ViralUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;

public final class segment_collect {


    private final static String PRG_NAME = "segment_collect";
    private static final String PRG_DATE = "2024-08-16";
    private static final String PRG_VERSION = "0.0.1";

    private final static int NUMBER_OF_SEGMENTS = 8;

    private final static int LOW_Q_THRESHOLD = 15;
    private final static String SEGMENT_OUTFILE_BASE = "sc_segment_";
    private final static String ALL_OUTFILE = "sc_all.fasta";

    private final static String SEP = "/";

    private static final String UNKNOWN = "unknown";

    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        if (args.length != 1) {
            System.out.println("\nWrong number of arguments, expected: <infile>\n");
            System.exit(-1);
        }

        final File infile = new File(args[0]);
        final File outfile_all = new File(ALL_OUTFILE);

        if (!infile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile + "] does not exist");
        }
        if (outfile_all.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + outfile_all + "] already exists");
        }
        for (int i = 0; i < NUMBER_OF_SEGMENTS; ++i) {
            final File seq_outfile = new File(SEGMENT_OUTFILE_BASE + (i + 1) + ".fasta");
            if (seq_outfile.exists()) {
                ForesterUtil.fatalError(PRG_NAME, "[" + seq_outfile + "] already exists");
            }
        }

        List<MolecularSequence> in_seqs = null;
        try {
            in_seqs = FastaParser.parse(new FileInputStream(infile));
        } catch (IOException e) {
            ForesterUtil.fatalError(PRG_NAME, "\nCould not read \"" + infile + "\" [" + e.getMessage() + "]");
        }

        final BasicDescriptiveStatistics infile_length_stats = new BasicDescriptiveStatistics();
        final BasicDescriptiveStatistics outfile_all_length_stats = new BasicDescriptiveStatistics();

        final SortedMap<String, Integer> strain_to_counts = new TreeMap();

        for (final MolecularSequence seq : in_seqs) {
            final BasicSequence bseq = (BasicSequence) seq;
            final String name = bseq.getIdentifier();

            infile_length_stats.addValue(bseq.getLength());

            final String new_name = makeNewName(name);
            if (new_name == null) continue;

            if (strain_to_counts.containsKey(new_name)) {
                strain_to_counts.put(new_name, strain_to_counts.get(new_name) + 1);
            } else {
                strain_to_counts.put(new_name, 1);
            }

        }

        final SortedMap<String, MolecularSequence[]> strain_to_segments = new TreeMap();

        int c = 0;
        int input_seqs_without_segment = 0;
        int input_seqs_failed_to_match = 0;
        int input_seqs_too_short = 0;
        int input_seqs_low_qual = 0;
        int not_eight = 0;
        for (final MolecularSequence seq : in_seqs) {

            final BasicSequence bseq = (BasicSequence) seq;
            final String name = bseq.getIdentifier();

            if (getIndCount(bseq) > LOW_Q_THRESHOLD) {
                ++input_seqs_low_qual;

                continue;
            }

            String new_name = makeNewName(name);
            if (new_name == null) {

                ++input_seqs_failed_to_match;
                continue;
            }

            if (strain_to_counts.get(new_name) == NUMBER_OF_SEGMENTS) {
                String name_lc = name.toLowerCase();
                final MolecularSequence[] current_seg_list;
                if (!strain_to_segments.containsKey(new_name)) {
                    MolecularSequence[] l = new MolecularSequence[8];
                    strain_to_segments.put(new_name, l);
                    current_seg_list = l;
                } else {
                    current_seg_list = strain_to_segments.get(new_name);
                }

                if (name_lc.indexOf("segment_1") > 0 || name_lc.indexOf("segment 1") > 0 ||
                        name_lc.indexOf("001_a_pb2_") > 0 || name_lc.indexOf(" pb2 ") > 0 ||
                        name_lc.indexOf("(pb2)") > 0 || name_lc.indexOf("_pb2|") > 0) {

                    if (seq.getLength() > 1500) {
                        current_seg_list[0] = seq;
                    } else {
                        ++input_seqs_too_short;
                    }
                } else if (name_lc.indexOf("segment_2") > 0 || name_lc.indexOf("segment 2") > 0 ||
                        name_lc.indexOf("001_a_pb1_") > 0 || name_lc.indexOf(" pb1 ") > 0 ||
                        name_lc.indexOf(" pb1,") > 0 || name_lc.indexOf("(pb1)") > 0 || name_lc.indexOf("_pb1|") > 0) {
                    if (seq.getLength() > 1500) {
                        current_seg_list[1] = seq;
                    } else {
                        ++input_seqs_too_short;
                    }
                } else if (name_lc.indexOf("segment_3") > 0 || name_lc.indexOf("segment 3") > 0 ||
                        name_lc.indexOf("001_a_pa_") > 0 || name_lc.indexOf(" pa ") > 0 ||
                        name_lc.indexOf("(pa)") > 0 || name_lc.indexOf("_pa|") > 0) {
                    if (seq.getLength() > 1500) {
                        current_seg_list[2] = seq;
                    } else {
                        ++input_seqs_too_short;
                    }
                } else if (name_lc.indexOf("segment_4") > 0 || name_lc.indexOf("segment 4") > 0 ||
                        name_lc.indexOf("001_a_ha_") > 0 || name_lc.indexOf(" ha ") > 0 ||
                        name_lc.indexOf("(ha)") > 0 || name_lc.indexOf("_ha|") > 0) {
                    if (seq.getLength() > 1200) {
                        current_seg_list[3] = seq;
                    } else {
                        ++input_seqs_too_short;
                    }
                } else if (name_lc.indexOf("segment_5") > 0 || name_lc.indexOf("segment 5") > 0 ||
                        name_lc.indexOf("001_a_np_") > 0 || name_lc.indexOf(" np ") > 0 ||
                        name_lc.indexOf("(np)") > 0 || name_lc.indexOf("_np|") > 0) {
                    if (seq.getLength() > 1000) {
                        current_seg_list[4] = seq;
                    } else {
                        ++input_seqs_too_short;
                    }
                } else if (name_lc.indexOf("segment_6") > 0 || name_lc.indexOf("segment 6") > 0 ||
                        name_lc.indexOf("001_a_na_") > 0 || name_lc.indexOf(" na ") > 0 ||
                        name_lc.indexOf("(na)") > 0 || name_lc.indexOf("_na|") > 0) {
                    if (seq.getLength() > 1000) {
                        current_seg_list[5] = seq;
                    } else {
                        ++input_seqs_too_short;
                    }
                } else if (name_lc.indexOf("segment_7") > 0 || name_lc.indexOf("segment 7") > 0 ||
                        name_lc.indexOf("001_a_mp_") > 0 || name_lc.indexOf(" mp ") > 0 || name_lc.indexOf("(mp)") > 0
                        || name_lc.indexOf("(m1)") > 0 || name_lc.indexOf("(m2)") > 0 ||
                        name_lc.indexOf("m1, m2") > 0 || name_lc.indexOf("_mp|") > 0) {
                    if (seq.getLength() > 900) {
                        current_seg_list[6] = seq;
                    } else {
                        ++input_seqs_too_short;
                    }
                } else if (name_lc.indexOf("segment_8") > 0 || name_lc.indexOf("segment 8") > 0
                        || name_lc.indexOf("001_a_ns_") > 0 || name_lc.indexOf(" ns ") > 0
                        || name_lc.indexOf("(ns)") > 0 || name_lc.indexOf("(nep)") > 0
                        || name_lc.indexOf("(ns1)") > 0 || name_lc.indexOf("ns1, ns2") > 0 || name_lc.indexOf("_ns|") > 0) {
                    if (seq.getLength() > 800) {
                        current_seg_list[7] = seq;
                    } else {
                        ++input_seqs_too_short;
                    }
                } else {
                    ++input_seqs_without_segment;
                    System.out.println("Warning: could not obtain segment from \"" + name + "\"");
                }
            } else {
                ++not_eight;
            }

        }

        final List<MolecularSequence> outseqs_all = new ArrayList();

        final List<MolecularSequence>[] outseqs_individual_segments = new List[NUMBER_OF_SEGMENTS];
        for (int i = 0; i < NUMBER_OF_SEGMENTS; ++i) {
            outseqs_individual_segments[i] = new ArrayList<>();
        }

        int input_strain_with_less_than_eight_segments = 0;

        F:
        for (final Map.Entry<String, MolecularSequence[]> entry : strain_to_segments.entrySet()) {
            final String strain = entry.getKey();
            final MolecularSequence[] s = entry.getValue();

            for (int i = 0; i < NUMBER_OF_SEGMENTS; ++i) {
                if (s[i] == null) {
                    ++input_strain_with_less_than_eight_segments;
                    continue F;
                }
            }

            final StringBuilder sb = new StringBuilder(0);
            for (int i = 0; i < NUMBER_OF_SEGMENTS; ++i) {
                // outseqs_individual_segments[i].add(new BasicSequence(strain + "_segment_" + (i + 1), s[i].getMolecularSequenceAsString(), MolecularSequence.TYPE.DNA));
                outseqs_individual_segments[i].add(new BasicSequence("(" + strain + ")", s[i].getMolecularSequenceAsString(), MolecularSequence.TYPE.DNA));

                sb.append(s[i].getMolecularSequenceAsString());
            }
            final MolecularSequence newseq = new BasicSequence(strain, sb.toString(), MolecularSequence.TYPE.DNA);
            outfile_all_length_stats.addValue(newseq.getLength());
            outseqs_all.add(newseq);
        }


        try {
            SequenceWriter.writeSeqs(outseqs_all, outfile_all, SEQ_FORMAT.FASTA, 80);
        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, "failed to write to [" + outfile_all + "]: " + e.getMessage());
        }
        System.out.println("[" + PRG_NAME + "] wrote " + outseqs_all.size() + " sequences to " + outfile_all);

        for (int i = 0; i < NUMBER_OF_SEGMENTS; ++i) {
            final File seq_outfile = new File(SEGMENT_OUTFILE_BASE + (i + 1) + ".fasta");
            try {
                SequenceWriter.writeSeqs(outseqs_individual_segments[i], seq_outfile, SEQ_FORMAT.FASTA, 80);
            } catch (final IOException e) {
                ForesterUtil.fatalError(PRG_NAME, "failed to write to [" + seq_outfile + "]: " + e.getMessage());
            }
            final BasicDescriptiveStatistics stats = new BasicDescriptiveStatistics();
            for (int j = 0; j < outseqs_individual_segments[i].size(); ++j) {
                stats.addValue(outseqs_individual_segments[i].get(j).getLength());
            }

            System.out.println("[" + PRG_NAME + "] wrote " + outseqs_individual_segments[i].size() + " sequences to " + seq_outfile);
            System.out.println();
            System.out.println("Length statistics for segment " + (i + 1));
            System.out.println("N     : " + stats.getN());
            System.out.println("Min   : " + stats.getMin());
            System.out.println("Max   : " + stats.getMax());
            System.out.println("IQR   : " + stats.interquartileRange());
            System.out.println("Median: " + stats.median());
            System.out.println("Mean  : " + stats.arithmeticMean());
            System.out.println("SD    : " + stats.sampleStandardDeviation());
            System.out.println();
        }


        System.out.println();
        System.out.println("Input length statistics:");
        System.out.println("N     : " + infile_length_stats.getN());
        System.out.println("Min   : " + infile_length_stats.getMin());
        System.out.println("Max   : " + infile_length_stats.getMax());
        System.out.println("IQR   : " + infile_length_stats.interquartileRange());
        System.out.println("Median: " + infile_length_stats.median());
        System.out.println("Mean  : " + infile_length_stats.arithmeticMean());
        System.out.println("SD    : " + infile_length_stats.sampleStandardDeviation());
        System.out.println();
        System.out.println("Output all length statistics:");
        System.out.println("N     : " + outfile_all_length_stats.getN());
        System.out.println("Min   : " + outfile_all_length_stats.getMin());
        System.out.println("Max   : " + outfile_all_length_stats.getMax());
        System.out.println("IQR   : " + outfile_all_length_stats.interquartileRange());
        System.out.println("Median: " + outfile_all_length_stats.median());
        System.out.println("Mean  : " + outfile_all_length_stats.arithmeticMean());
        System.out.println("SD    : " + outfile_all_length_stats.sampleStandardDeviation());
        System.out.println();
        System.out.println("Number of sequences in input  : " + in_seqs.size());
        System.out.println("Number of distinct isolates   : " + strain_to_counts.size());
        System.out.println("   low quality                : " + input_seqs_low_qual);
        System.out.println("   named failed to match      : " + input_seqs_failed_to_match);
        System.out.println("   without segment information: " + input_seqs_without_segment);
        System.out.println("   not eight segments         : " + not_eight);
        System.out.println("   too short                  : " + input_seqs_too_short);
        System.out.println("Distinct isolates after qc    : " + strain_to_segments.size());
        System.out.println("With missing segment after qc : " + input_strain_with_less_than_eight_segments);


        System.out.println();
        System.out.println("Number of sequences in output : " + outseqs_all.size());

        System.out.println("[" + PRG_NAME + "] OK");
        System.out.println();

    }

    private static int getIndCount(final BasicSequence bseq) {
        int ind_count = 0;
        for (int i = 0; i < bseq.getLength(); ++i) {
            char r = bseq.getResidueAt(i);
            if (r != 'a' && r != 'c' && r != 'g' && r != 't' && r != 'A' && r != 'C' && r != 'G' && r != 'T') {
                ++ind_count;
            }
        }
        return ind_count;
    }

    private static String makeNewName(final String name) {
        final Matcher m0 = ViralUtils.PATTERN_0.matcher(name);
        final Matcher m1 = ViralUtils.PATTERN_1.matcher(name);
        final Matcher m2 = ViralUtils.PATTERN_2.matcher(name);
        final Matcher m1b = ViralUtils.PATTERN_1b.matcher(name);
        final Matcher m3 = ViralUtils.PATTERN_3.matcher(name);
        final Matcher mg = ViralUtils.PATTERN_G.matcher(name);
        String type = UNKNOWN;
        String host = UNKNOWN;
        String location = UNKNOWN;
        String strain_number = UNKNOWN;
        String year = UNKNOWN;
        String subtype = UNKNOWN;

        if (m0.find()) {
            // 1. type
            // 2. host
            // 3. country/state
            // 4. strain number
            // 5. year
            // 6. subtype
            type = m0.group(1).trim();
            host = m0.group(2).trim();
            location = m0.group(3).trim();
            strain_number = m0.group(4).trim();
            year = m0.group(5).trim();
            subtype = m0.group(6).trim().toUpperCase();
        } else if (m1.find()) {
            // 1. type
            // 2. host
            // 3. country/state
            // 4. number
            // 5. year
            // 6. subtype
            type = m1.group(1).trim();
            host = m1.group(2).trim();
            location = m1.group(3).trim();
            strain_number = m1.group(4).trim();
            year = m1.group(5).trim();
            subtype = m1.group(6).trim().toUpperCase();
        } else if (m1b.find()) {
            // 1. type
            // 2. host
            // 3. country/state
            // 4. number
            // 5. year
            // 6. subtype
            type = m1b.group(1).trim();
            host = m1b.group(2).trim();
            location = m1b.group(3).trim();
            strain_number = m1b.group(4).trim();
            year = m1b.group(5).trim();
            subtype = m1b.group(6).trim().toUpperCase();

        } else if (m2.find()) {
            // 1. type
            // 2. country/state
            // 3. number
            // 4. year
            // 5. subtype
            host = "human";
            type = m2.group(1).trim();
            location = m2.group(2).trim();
            strain_number = m2.group(3).trim();
            year = m2.group(4).trim();
            subtype = m2.group(5).trim().toUpperCase();
        } else if (m3.find()) {
            // 1. type
            // 2. host
            // 3. country/state
            // 4. number
            // 5. year
            type = "A";
            host = m3.group(2).trim();
            location = m3.group(3).trim();
            strain_number = m3.group(4).trim();
            year = m3.group(5).trim();
            subtype = "H5N1";
        } else if (mg.find()) {
            type = "A";
            host = "human";
            location = mg.group(2).trim();
            strain_number = mg.group(3).trim();
            year = mg.group(4).trim();
            subtype = "H5N1";

        } else {
            System.out.println("Warning: name \"" + name + "\" could not be matched");
            return null;
        }

        host = host.replace('_', ' ');
        location = location.replace('_', ' ');

        host = host.replaceAll("\\s+", " ");
        location = location.replaceAll("\\s+", " ");

        host = ViralUtils.cleanHost(host);

        host = ViralUtils.cleanHostString(host);
        location = ViralUtils.cleanLocationString(location);

        String new_name = type + SEP + host + SEP + location + SEP + strain_number + SEP + year + "(" + subtype + ")";

        new_name = new_name.replaceAll("\\s+", "_");
        return new_name;
    }
}
