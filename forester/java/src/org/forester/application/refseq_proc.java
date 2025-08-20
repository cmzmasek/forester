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
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class refseq_proc {

    private final static String PRG_NAME = "refseq_proc";
    private static final String PRG_DATE = "2025-08-20";
    private static final String PRG_VERSION = "1.0.0";
    public final static Pattern PATTERN_GENE = Pattern
            .compile("\\[gene=(.+?)\\]");
    public final static Pattern PATTERN_PROTEIN = Pattern
            .compile("\\[protein=(.+?)\\]");
    public final static Pattern PATTERN_PROTEIN_ID = Pattern
            .compile("\\[protein_id=(.+?)\\]");
    public final static Pattern PATTERN_LOCATION = Pattern
            .compile("\\[location=.+?(\\d+?)\\]");
    public final static Pattern PATTERN_GBKEY = Pattern
            .compile("\\[gbkey=(.+?)\\]");


    public static void main(final String args[]) {

        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        if (args.length != 2) {
            System.out.println("\nWrong number of arguments, expected: <infile> <outfile>\n");
            System.exit(-1);
        }

        try {
            final File infile = new File(args[0]);
            final File outfile = new File(args[1]);
            final List<MolecularSequence> seqs = FastaParser.parse(new FileInputStream(infile));
            for (final MolecularSequence seq : seqs) {
                final BasicSequence bseq = (BasicSequence) seq;
                final String id = bseq.getIdentifier();

                String gene = "";
                final Matcher mg = PATTERN_GENE.matcher(id);
                if (mg.find()) {
                    gene = mg.group(1).trim();
                    System.out.println(gene);
                } else {
                    System.out.println(bseq.getIdentifier() + " did not match");
                }

                String protein = "";
                final Matcher mp = PATTERN_PROTEIN.matcher(id);
                if (mp.find()) {
                    protein = mp.group(1).trim();
                    System.out.println(protein);
                } else {
                    System.out.println(bseq.getIdentifier() + " did not match");
                }

                String protein_id = "";
                final Matcher mpi = PATTERN_PROTEIN_ID.matcher(id);
                if (mpi.find()) {
                    protein_id = mpi.group(1).trim();
                    System.out.println(protein_id);
                } else {
                    System.out.println(bseq.getIdentifier() + " did not match");
                }

                String loc = "";
                final Matcher ml = PATTERN_LOCATION.matcher(id);
                if (ml.find()) {
                    loc = ml.group(1).trim();
                    System.out.println(loc);
                } else {
                    System.out.println(bseq.getIdentifier() + " did not match");
                }

                String gbkey = "";
                final Matcher mk = PATTERN_GBKEY.matcher(id);
                if (mk.find()) {
                    gbkey = mk.group(1).trim();
                    System.out.println(gbkey);
                } else {
                    System.out.println(bseq.getIdentifier() + " did not match");
                }


            }
            SequenceWriter.writeSeqs(seqs, outfile, SEQ_FORMAT.FASTA, 60);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
