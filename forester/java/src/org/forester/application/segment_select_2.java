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
import org.forester.util.ForesterUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class segment_select_2 {

    private static final String PRG_DATE = "2025-10-20";
    private static final String PRG_VERSION = "2.0.0";

    private static final boolean VERBOSE = true;

    private final static String PRG_NAME = "segment_select_2";

    // public final static Pattern P0 = Pattern
    //    .compile("(\\(A/.+?)[_\\s]segment[_\\s]");

    //  public final static Pattern P0 = Pattern
    //          .compile("(\\(A/.+?\\)\\))");

    public final static Pattern P0 = Pattern
            .compile("(A/.+?\\))");


    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);

        if (args.length != 3 && args.length != 4) {
            System.out.println("\nWrong number of arguments, expected: [regex] <fasta file to be trimmed> <\"guide\" fasta file> <outfile>\n");
            System.exit(-1);
        }
        String regex;
        final File infile_to_be_trimmed;
        final File infile_guide;
        final File outfile;

        if (args.length == 4) {
            regex = args[0];
            infile_to_be_trimmed = new File(args[1]);
            infile_guide = new File(args[2]);
            outfile = new File(args[3]);
        } else {
            regex = null;
            infile_to_be_trimmed = new File(args[0]);
            infile_guide = new File(args[1]);
            outfile = new File(args[2]);
        }
        if (!infile_to_be_trimmed.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile_to_be_trimmed + "] does not exist");
        }

        if (!infile_guide.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile_guide + "] does not exist");
        }
        if (outfile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + outfile + "] already exists");
        }
        final Pattern compiled_regex;
        if (regex != null) {
            compiled_regex = Pattern.compile(regex);
        } else {
            compiled_regex = P0;
        }

        Map<String, MolecularSequence> outseqs = new HashMap<>();

        List<MolecularSequence> seqs_to_be_trimmed;

        try {
            seqs_to_be_trimmed = FastaParser.parse(new FileInputStream(infile_to_be_trimmed));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        List<MolecularSequence> seqs_guide;

        try {
            seqs_guide = FastaParser.parse(new FileInputStream(infile_guide));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        final SortedSet matches  = new TreeSet<String>();

        for (MolecularSequence s_trimmed : seqs_to_be_trimmed) {
            final Matcher m0 = compiled_regex.matcher(s_trimmed.getIdentifier());

            String virus_id_tr = "";
            if (m0.find()) {
                virus_id_tr = m0.group(1).trim().toLowerCase().replaceAll("\\s+", "_");
            } else {
                System.out.println("Error: Sequence \"" + s_trimmed.getIdentifier() + "\" could not be matched");
                //System.exit(-1);
                continue;
            }
            if (matches.contains(virus_id_tr)) {
                System.out.println("Error: Match \"" +virus_id_tr + "\" is not specific");
                System.out.println("     " + s_trimmed.getIdentifier());
                System.exit(-1);
            }
            else {
                matches.add(virus_id_tr);
            }



            int counter = 0;
            boolean found = false;
            for (MolecularSequence sguide : seqs_guide) {
                final Matcher m1 = compiled_regex.matcher(sguide.getIdentifier());
                ++counter;
                String virus_id_g = "";
                if (m1.find()) {
                    virus_id_g = m1.group(1).trim().toLowerCase().replaceAll("\\s+", "_");
                } else {
                    System.out.println("Error: Guide sequence " + counter + " \"" + sguide.getIdentifier() + "\" could not be matched");
                    System.exit(-1);
                }
                // if (virus_id_g.split("/").length > 4 ) {
                if (virus_id_g.equals(virus_id_tr)) {
                    if ( found) {
                        System.out.println("Error: " + virus_id_g + " is not specific");
                        System.out.println("Guide: " + sguide.getIdentifier());
                        System.out.println("     : " + s_trimmed.getIdentifier());
                        System.exit(-1);
                    }
                    found = true;
                    System.out.println(" == " + virus_id_g);
                    outseqs.put(s_trimmed.getIdentifier(), s_trimmed);
                }
                // }
            }
        }

        List<MolecularSequence> outseqlist = new ArrayList<MolecularSequence>(outseqs.values());

        try {
            SequenceWriter.writeSeqs(outseqlist, outfile, SEQ_FORMAT.FASTA, 60);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        System.out.println();

        System.out.println("[" + PRG_NAME + "] orig    : " + seqs_to_be_trimmed.size() + " sequences");
        System.out.println("[" + PRG_NAME + "] guide   : " + seqs_guide.size() + " sequences");
        System.out.println("[" + PRG_NAME + "] retained: " + outseqlist.size() + " sequences");
        System.out.println("[" + PRG_NAME + "] wrote: [" + outfile + "]");
        System.out.println("[" + PRG_NAME + "] OK");
        System.out.println();
    }
}
