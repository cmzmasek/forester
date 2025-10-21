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
import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class rem_dup_regex {

    private static final String PRG_DATE = "2025-10-20";
    private static final String PRG_VERSION = "1.0.0";
    private static final boolean VERBOSE = true;
    private final static String PRG_NAME = "rem_dup_regex";

    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);

        if (args.length != 3) {
            System.out.println("\nWrong number of arguments, expected: [regex] <infile> <outfile>\n");
            System.exit(-1);
        }

        final String regex = args[0];
        final File infile = new File(args[1]);
        final File outfile = new File(args[2]);

        if (!infile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile + "] does not exist");
        }
        if (outfile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + outfile + "] already exists");
        }
        final Pattern compiled_regex = Pattern.compile(regex);
        List<MolecularSequence> seqs;

        try {
            seqs = FastaParser.parse(new FileInputStream(infile));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        final List<String> names = new ArrayList<>();
        for (final MolecularSequence s : seqs) {
            final Matcher m_i = compiled_regex.matcher(s.getIdentifier());
            if (m_i.find()) {
                names.add(m_i.group(1).trim().toLowerCase().replaceAll("\\s+", "_"));
            }
        }

        int no_match = 0;
        int multiples = 0;
        final List<MolecularSequence> outseqlist = new ArrayList<>();
        for (final MolecularSequence s : seqs) {
            final String name_o = s.getIdentifier();
            final Matcher m_o = compiled_regex.matcher(name_o);
            String id_o = "";
            if (m_o.find()) {
                id_o = m_o.group(1).trim().toLowerCase().replaceAll("\\s+", "_");
            } else {
                System.out.println("\"" + s.getIdentifier() + "\" could not be matched");
                ++no_match;
                continue;
            }
            int count = 0;
            for (final String name : names) {
                if (name.equals(id_o)) {
                    ++count;
                }
            }
            if (count > 1) {
                ++multiples;
                System.out.println("\"" + name_o + "\":  " + count);
            } else if (count == 1) {
                outseqlist.add(s);
            } else {
                ForesterUtil.fatalError(PRG_NAME, "This should never have happened.");
            }
        }

        try {
            SequenceWriter.writeSeqs(outseqlist, outfile, SEQ_FORMAT.FASTA, 60);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        System.out.println();
        System.out.println("[" + PRG_NAME + "] version   : " + PRG_VERSION);
        System.out.println("[" + PRG_NAME + "] infile    : " + infile);
        System.out.println("[" + PRG_NAME + "] regex     : " + regex);
        System.out.println("[" + PRG_NAME + "] input     : " + seqs.size() + " sequences");
        System.out.println("[" + PRG_NAME + "] no match  : " + no_match);
        System.out.println("[" + PRG_NAME + "] unspecific: " + multiples);
        System.out.println("[" + PRG_NAME + "] retained  : " + outseqlist.size() + " sequences");
        System.out.println("[" + PRG_NAME + "] wrote     : [" + outfile + "]");
        System.out.println("[" + PRG_NAME + "] OK");
        System.out.println();
    }
}
