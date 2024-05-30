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
import org.forester.util.ViralUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;

public final class rename_fasta_vipr_x4 {


    private final static String PRG_NAME = "rename_fasta_vipr_x4";
    private static final String PRG_DATE = "2024-05-30";
    private static final String PRG_VERSION = "1.0.0";

    private static final String UNKNOWN = "unknown";

    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        if (args.length != 2) {
            System.out.println("\nWrong number of arguments, expected: <infile> <outfile>\n");
            System.exit(-1);
        }

        final File infile = new File(args[0]);
        final File outfile = new File(args[1]);

        if (!infile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile + "] does not exist");
        }
        if (outfile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + outfile + "] already exists");
        }

        List<MolecularSequence> seqs = null;
        try {
            seqs = FastaParser.parse(new FileInputStream(infile));
        } catch (IOException e) {
            System.out.println("\nCould not read \"" + infile + "\" [" + e.getMessage() + "]\n");
            System.exit(-1);
        }
        for (final MolecularSequence seq : seqs) {
            final BasicSequence bseq = (BasicSequence) seq;
            final String name = bseq.getIdentifier();

            final Matcher mg = ViralUtils.PATTERN_GB.matcher(name);
            final Matcher m0 = ViralUtils.PATTERN_0.matcher(name);
            final Matcher m1 = ViralUtils.PATTERN_1.matcher(name);
            final Matcher m2 = ViralUtils.PATTERN_2.matcher(name);
            final Matcher m1b = ViralUtils.PATTERN_1b.matcher(name);
            String type = UNKNOWN;
            String host = UNKNOWN;
            String location = UNKNOWN;
            String strain_number = UNKNOWN;
            String year = UNKNOWN;
            String subtype = UNKNOWN;
            String genbank_acc = "";
            if (mg.find()) {
                genbank_acc = mg.group(1).trim();
            }

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
                // 7. acc
                type = m1b.group(1).trim();
                host = m1b.group(2).trim();
                location = m1b.group(3).trim();
                strain_number = m1b.group(4).trim();
                year = m1b.group(5).trim();
                subtype = m1b.group(6).trim().toUpperCase();
                genbank_acc = m1b.group(7).trim();

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
            } else {
                System.out.println("ERROR: name \"" + name + "\" could not be matched");
                System.exit(-1);
            }

            host = host.replace('_', ' ');
            location = location.replace('_', ' ');

            host = ViralUtils.cleanHost(host);

            host = ViralUtils.cleanHostOrLocationString(host);
            location = ViralUtils.cleanHostOrLocationString(location);

            year = ViralUtils.checkYear(year);

            String new_name = type + "/" + host + "/" + location + "/" + strain_number + "/" + year + "|" + subtype;

            if (genbank_acc.length() > 5) {
                new_name += "|" + genbank_acc;
            }
            new_name = new_name.replace(' ', '_');
            bseq.setIdentifier(new_name);

        }
        try {
            SequenceWriter.writeSeqs(seqs, outfile, SEQ_FORMAT.FASTA, 60);
        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, "failed to write to [" + outfile + "]: " + e.getMessage());
        }

        System.out.println();
        System.out.println("[" + PRG_NAME + "] wrote: [" + outfile + "]");
        System.out.println("[" + PRG_NAME + "] OK");
        System.out.println();

    }


}
