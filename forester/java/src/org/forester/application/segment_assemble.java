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
import java.util.regex.Pattern;


// Example:
// segment_assemble "(.+?)\|.*?\|.*?\|.*?\|.*?\|([^|]{3,}?)\|" segment_S_600.fasta segment_M_3000.fasta segment_L_5000.fasta
public final class segment_assemble {


    private final static String PRG_NAME = "segment_assemble";
    private static final String PRG_DATE = "2026-01-07";
    private static final String PRG_VERSION = "0.0.1";

    private final static String SEGMENT_OUTFILE_BASE = "sa_segment_";
    private final static String ALL_OUTFILE = "sa_all.fasta";

    private static final String UNKNOWN = "unknown";

    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        if (args.length != 4 && args.length != 9) {
            System.out.println("\nWrong number of arguments, expected: <segment 1,S> <segment 2,M> <segment 3,L> OR <segment 1 - 8>\n");
            System.exit(-1);
        }
        final String regex = args[0];
        final Pattern pattern = Pattern.compile(regex);

        ForesterUtil.programMessage(PRG_NAME, "regex: " + pattern);
        final int number_of_segments = args.length - 1;
        ForesterUtil.programMessage(PRG_NAME, "number of segments: " + number_of_segments);
        final List<File> infiles = new ArrayList<File>();
        for (int i = 0; i < number_of_segments; ++i) {
            infiles.add(new File(args[i + 1]));
            if (!infiles.get(i).exists()) {
                final String e = ForesterUtil.isReadableFile(infiles.get(i));
                if (!ForesterUtil.isEmpty(e)) {
                    ForesterUtil.fatalError(PRG_NAME, e);
                }
            }
            ForesterUtil.programMessage(PRG_NAME, "infile " + i + ": " + infiles.get(i));
        }

        final File outfile_all = new File(ALL_OUTFILE);


        if (outfile_all.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + outfile_all + "] already exists");
        }
        for (int i = 0; i < number_of_segments; ++i) {
            final File seq_outfile = new File(SEGMENT_OUTFILE_BASE + (i + 1) + ".fasta");
            if (seq_outfile.exists()) {
                ForesterUtil.fatalError(PRG_NAME, "[" + seq_outfile + "] already exists");
            }
        }

        List<List<MolecularSequence>> in_segments = new ArrayList<>();
        for (int i = 0; i < number_of_segments; ++i) {
            List<MolecularSequence> in_seqs = null;
            try {
                in_seqs = FastaParser.parse(new FileInputStream(infiles.get(i)));
            } catch (IOException e) {
                ForesterUtil.fatalError(PRG_NAME, "\nCould not read \"" + infiles.get(i) + "\" [" + e.getMessage() + "]");
            }
            in_segments.add(in_seqs);
        }

        printInputSegmentsLengthStats(number_of_segments, in_segments);

        final SortedMap<String, Integer> strain_to_counts = countSegmentsPerStrain(pattern, number_of_segments, in_segments);

        System.out.println(strain_to_counts);

        final SortedMap<String, List<MolecularSequence>> perstrain = organizeSeqsPerStrain(pattern, number_of_segments, in_segments, strain_to_counts);

       // System.out.println(perstrain);

        System.out.println("[" + PRG_NAME + "] OK");
        System.out.println();

    }

    private static SortedMap<String, Integer> countSegmentsPerStrain(final Pattern pattern,
                                                                     final int number_of_segments,
                                                                     final List<List<MolecularSequence>> in_segments) {
        final SortedMap<String, Integer> strain_to_counts = new TreeMap();
        for (int i = 0; i < number_of_segments; ++i) {
            final List<MolecularSequence> seqs = in_segments.get(i);
            for (final MolecularSequence seq : seqs) {
                final String strain = extractStrain(seq.getIdentifier(), pattern);
                if (!ForesterUtil.isEmpty(strain)) {
                    if (strain_to_counts.containsKey(strain)) {
                        strain_to_counts.put(strain, strain_to_counts.get(strain) + 1);
                    } else {
                        strain_to_counts.put(strain, 1);
                    }
                }
            }
        }
        return strain_to_counts;
    }

    private static SortedMap<String, List<MolecularSequence>> organizeSeqsPerStrain(final Pattern pattern,
                                                                                    final int number_of_segments,
                                                                                    final List<List<MolecularSequence>> in_segments,
                                                                                    final SortedMap<String, Integer> counts) {
        final SortedMap<String, List<MolecularSequence>> strain_to_seqs = new TreeMap();
        for (int seg = 0; seg < number_of_segments; ++seg) {
            final List<MolecularSequence> seqs = in_segments.get(seg);
            for (final MolecularSequence seq : seqs) {
                final String strain = extractStrain(seq.getIdentifier(), pattern);
                if (!ForesterUtil.isEmpty(strain)) {
                    if (counts.get(strain) == number_of_segments) {
                        if (!strain_to_seqs.containsKey(strain)) {
                            strain_to_seqs.put(strain, new ArrayList());
                        }
                        List<MolecularSequence> x = strain_to_seqs.get(strain);
                        x.add(seq);
                    }
                }
            }
        }
        return strain_to_seqs;
    }

    private static void printInputSegmentsLengthStats(final int number_of_segments,
                                                      final List<List<MolecularSequence>> in_segments) {
        for (int i = 0; i < number_of_segments; ++i) {
            final BasicDescriptiveStatistics stats = new BasicDescriptiveStatistics();
            final List<MolecularSequence> seqs = in_segments.get(i);
            for (final MolecularSequence seq : seqs) {
                stats.addValue(seq.getLength());
            }
            System.out.println("Length statistics for input segment " + (i + 1));
            System.out.println("N     : " + stats.getN());
            System.out.println("Min   : " + stats.getMin());
            System.out.println("Max   : " + stats.getMax());
            System.out.println("IQR   : " + stats.interquartileRange());
            System.out.println("Median: " + stats.median());
            System.out.println("Mean  : " + stats.arithmeticMean());
            System.out.println("SD    : " + stats.sampleStandardDeviation());
        }
    }


    private static String extractStrain(final String name, final Pattern p) {
        final Matcher m = p.matcher(name);
        if (m.find()) {
            return m.group(2).trim();
        } else {
            return null;
        }
    }

    private static String extractSeqId(final String name, final Pattern p) {
        final Matcher m = p.matcher(name);
        if (m.find()) {
            return m.group(1).trim();
        } else {
            return null;
        }
    }

}
