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
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterUtil;
import org.forester.util.ViralUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;

public final class segment_select {

    private static final String PRG_DATE = "2024-05-17";
    private static final String PRG_VERSION = "1.0.0";

    private static final boolean VERBOSE = false;

    private final static String PRG_NAME = "segment_select";


    public static void main(final String args[]) {


        if (args.length != 3) {
            System.out.println("\nWrong number of arguments, expected: <fasta file> <intree> <outfile>\n");
            System.exit(-1);
        }

        final File infile = new File(args[0]);
        final File intree = new File(args[1]);
        final File outfile = new File(args[2]);

        if (!infile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile + "] does not exist");
        }

        if (!intree.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + intree + "] does not exist");
        }
        if (outfile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + outfile + "] already exists");
        }

        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(intree, true);
            p = factory.create(intree, pp)[0];
        } catch (final Exception e) {
            System.out.println("\nCould not read \"" + infile + "\" [" + e.getMessage() + "]\n");
            System.exit(-1);
        }
        final List<PhylogenyNode> ext_nodes = p.getExternalNodes();

        final ArrayList<String> strain_numbers = new ArrayList<>();
        final ArrayList<String> hosts = new ArrayList<>();
        final ArrayList<String> countries = new ArrayList<>();
        final ArrayList<String> years = new ArrayList<>();

        List<MolecularSequence> seqs;
        Map<String, MolecularSequence> outseqs = new HashMap<>();
        try {
            seqs = FastaParser.parse(new FileInputStream(infile));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        for (MolecularSequence s : seqs) {
            BasicSequence seq = (BasicSequence) s;
            final String seq_name = seq.getIdentifier();

            final Matcher m0 = ViralUtils.PATTERN_0.matcher(seq_name);
            final Matcher m1 = ViralUtils.PATTERN_1.matcher(seq_name);
            final Matcher m2 = ViralUtils.PATTERN_2.matcher(seq_name);
            final Matcher m1b = ViralUtils.PATTERN_1b.matcher(seq_name);

            String host = "";
            String location = "";
            String strain_number = "";
            String year = "";

            if (m0.find()) {
                host = m0.group(2).trim();
                location = m0.group(3).trim();
                strain_number = m0.group(4).trim();
                year = m0.group(5).trim();

            } else if (m1.find()) {
                host = m1.group(2).trim();
                location = m1.group(3).trim();
                strain_number = m1.group(4).trim();
                year = m1.group(5).trim();

            } else if (m1b.find()) {
                host = m1b.group(2).trim();
                location = m1b.group(3).trim();
                strain_number = m1b.group(4).trim();
                year = m1b.group(5).trim();

            } else if (m2.find()) {
                host = "human";
                location = m2.group(2).trim();
                strain_number = m2.group(3).trim();
                year = m2.group(4).trim();

            } else {
                System.out.println("Warning: MSA name \"" + seq_name + "\" could not be matched");
                continue;
            }

            host = ViralUtils.cleanHost(host);

            host = ViralUtils.cleanHostOrLocationString(host);
            location = ViralUtils.cleanHostOrLocationString(location);

            final String country = ViralUtils.determineCountry(location);

            if (VERBOSE) {
                System.out.println();
                System.out.println();
                System.out.println("Name    : " + seq_name);
                System.out.println("Host    : " + host);
                System.out.println("Location: " + location);
                System.out.println("Country : " + country);
                System.out.println("Number  : " + strain_number);
                System.out.println("Year    : " + year);
            }

            year = ViralUtils.checkYear(year);

            for (PhylogenyNode node : ext_nodes) {
                final String node_name = node.getName();
                final Matcher m1b2 = ViralUtils.PATTERN_1b.matcher(node_name);

                if (m1b2.find()) {
                    final String node_host = m1b2.group(2).trim();
                    final String node_location = m1b2.group(3).trim();
                    final String node_strain_number = m1b2.group(4).trim();
                    final String node_year = m1b2.group(5).trim();

                    if (node_host.equalsIgnoreCase(host)
                            && node_year.equalsIgnoreCase(year)
                            && node_strain_number.equalsIgnoreCase(strain_number)
                            && node_location.equalsIgnoreCase(location)
                    ) {
                        outseqs.put(seq_name, seq);
                    }
                } else {
                    System.out.println("ERROR: node name \"" + node_name + "\" could not be matched");
                    System.exit(-1);
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
        System.out.println("[" + PRG_NAME + "] retained " + outseqlist.size() + " sequences");
        System.out.println("[" + PRG_NAME + "] wrote: [" + outfile + "]");
        System.out.println("[" + PRG_NAME + "] OK");
        System.out.println();
    }
}
