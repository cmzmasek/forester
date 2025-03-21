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

package org.forester.application;

import org.forester.io.parsers.FastaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;
import org.forester.util.ViralUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class name_cons {

    final static private String PRG_NAME = "name_coms";
    final static private String PRG_VERSION = "1.00";
    final static private String PRG_DATE = "250319";

    public final static Pattern P1 = Pattern.compile("^.+?\\s+(.+?)\\s+\\[");

    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(name_cons.PRG_NAME, name_cons.PRG_VERSION, name_cons.PRG_DATE);
        System.out.println();
        if ((args.length != 2)) {
            name_cons.argumentsError();
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments(args);
        } catch (final Exception e) {
            ForesterUtil.fatalError(PRG_NAME, e.getMessage());
        }
        final File indir = cla.getFile(0);
        final File outdir = cla.getFile(1);
        if (!indir.isDirectory()) {
            ForesterUtil.fatalError(PRG_NAME, indir + " is not a directory");
        }
        if (!outdir.isDirectory()) {
            ForesterUtil.fatalError(PRG_NAME, outdir + " is not a directory");
        }
        final File[] list_of_files = indir.listFiles();
        final List<File> infiles = new ArrayList<File>();
        for (final File file : list_of_files) {
            if (file.isFile()
                    && file.canRead()
                    && (file.toString().toLowerCase().endsWith(".fasta") || file.toString().toLowerCase()
                    .endsWith(".fas"))) {
                infiles.add(file);
            }
        }
        Collections.sort(infiles);
        int c = 0;
        for (final File infile : infiles) {
            //System.out.println(++c + "/" + infiles.size() + ": " + infile);
            execute(outdir, infile);
        }
    }

    private static void execute(final File outdir, final File infile) {
        final File outfile = new File(outdir.getAbsolutePath().toString() + "/" + infile.getName());
        if (outfile.exists()) {
            System.out.println(outfile + " already exists");
        } else {
            try {
                final List<MolecularSequence> seqs = FastaParser.parse(new FileInputStream(infile));
                final Map<String, Integer> names = new HashMap<String, Integer>();

                for (final MolecularSequence seq : seqs) {


                    final String id = seq.getIdentifier();
                    //System.out.println(id);
                    final Matcher m = P1.matcher(id);
                    if (m.find()) {
                        final String name = m.group(1).trim().toLowerCase();
                        //System.out.println("-> " + name );
                        if (!name.startsWith("hypothetical")) {
                            if (names.containsKey(name)) {
                                names.put(name, names.get(name) + 1);
                            } else {
                                names.put(name, 1);
                            }
                        }

                    }

                }
                //System.out.println(names);
                List<Map.Entry<String, Integer>> list = new ArrayList<>(names.entrySet());
                list.sort(Map.Entry.comparingByValue());
                //System.out.println(list);
                System.out.print(infile.getName().replace(".fasta", ""));
                System.out.print('\t');
                if (list.size() > 0 ) {
                    String top_name = list.get(list.size() - 1).getKey();
                    top_name = top_name.replace("rna", "RNA");
                    top_name = top_name.replace("dna", "DNA");
                    final String s1 = top_name.substring(0, 1).toUpperCase();
                    final String top_name_cap = s1 +top_name.substring(1);
                    System.out.print(top_name_cap);
                }
                else {
                    System.out.print("");
                }
                System.out.print('\t');
                if (list.size() > 0 ) {
                    System.out.print(list.get(list.size() - 1).getValue() + "/" + seqs.size());
                }
                else {
                    System.out.print("0/" + seqs.size());
                }
                System.out.println();

            } catch (final IOException e) {
                ForesterUtil.fatalError(PRG_NAME, e.getMessage());
            }
        }
    }


    private static void argumentsError() {
        System.out.println(PRG_NAME + " <indir> <outdir>");
        System.out.println();
        System.exit(-1);
    }
}
