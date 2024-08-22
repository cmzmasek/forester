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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class segment_select_2 {

    private static final String PRG_DATE = "2024-05-22";
    private static final String PRG_VERSION = "1.0.1";

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

        if (args.length != 3) {
            System.out.println("\nWrong number of arguments, expected: <fasta file to be trimmed> <fasta file> <outfile>\n");
            System.exit(-1);
        }

        final File infile_to_be_trimmed = new File(args[0]);
        final File infile_guide = new File(args[1]);
        final File outfile = new File(args[2]);

        if (!infile_to_be_trimmed.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile_to_be_trimmed + "] does not exist");
        }

        if (!infile_guide.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile_guide + "] does not exist");
        }
        if (outfile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + outfile + "] already exists");
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

        for (MolecularSequence s_trimmed : seqs_to_be_trimmed) {
            final Matcher m0 = P0.matcher(s_trimmed.getIdentifier());

            String virus_id_tr = "";
            if (m0.find()) {
                virus_id_tr = m0.group(1).trim().toLowerCase().replaceAll( "\\s+", "_");
            } else {
                System.out.println("Error: MSA name \"" + s_trimmed.getIdentifier() + "\" could not be matched");
                //System.exit(-1);
                continue;
            }

            I: for (MolecularSequence sguide : seqs_guide) {
                final Matcher m1 = P0.matcher(sguide.getIdentifier());

                String virus_id_g = "";
                if (m1.find()) {
                    virus_id_g = m1.group(1).trim().toLowerCase().replaceAll( "\\s+", "_");
                } else {
                    System.out.println("Error: MSA name \"" + sguide.getIdentifier() + "\" could not be matched");
                    System.exit(-1);
                }

                if (virus_id_g.equals(virus_id_tr)) {
                    outseqs.put(s_trimmed.getIdentifier(), s_trimmed);
                    break I;
                }
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
