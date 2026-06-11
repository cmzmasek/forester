// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.application;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class vipr_xc {

    private static final String PRG_DATE = "2025-01-27";
    private static final String PRG_VERSION = "1.1.0";
    private final static String PRG_NAME = "vipr_xc";
    private static final String XSD_STRING = "xsd:string";
    private static final String CLADE = "vipr:H3_clade";
    private static final String SHORTCLADE = "vipr:H3_shortclade";
    private static final String LEGACYCLADE = "vipr:H3_legacyclade";


    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        if (args.length != 3) {
            System.out.println("\nWrong number of arguments, expected: <intree> <clade map> <outtree>\n");
            System.exit(-1);
        }

        final File intree = new File(args[0]);
        final File clade_infile = new File(args[1]);
        final File outfile = new File(args[2]);
        if (!intree.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + intree + "] does not exist");
        }
        if (!clade_infile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + clade_infile + "] does not exist");
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
            System.out.println("\nCould not read \"" + intree + "\" [" + e.getMessage() + "]\n");
            System.exit(-1);
        }
        BasicTable<String> isolate_to_clade = null;
        try {
            isolate_to_clade = BasicTableParser.parse(clade_infile, '\t', false, false);
        } catch (final Exception e) {
            ForesterUtil.fatalError(PRG_NAME, "failed to read [" + clade_infile + "] [" + e.getMessage() + "]");
        }


        final SortedMap<String, String> isolate_to_clade_map = new TreeMap<String, String>();
        final SortedMap<String, String> isolate_to_shortclade_map = new TreeMap<String, String>();
        final SortedMap<String, String> isolate_to_legacyclade_map = new TreeMap<String, String>();
        for (int row = 0; row < isolate_to_clade.getNumberOfRows(); ++row) {
            final String key = isolate_to_clade.getValue(1, row);
            final String clade = isolate_to_clade.getValue(2, row);
            final String shortclade = isolate_to_clade.getValue(4, row);
            final String legacyclade = isolate_to_clade.getValue(5, row);
            isolate_to_clade_map.put(key, clade);
            isolate_to_shortclade_map.put(key, shortclade);
            isolate_to_legacyclade_map.put(key, legacyclade);
        }

        final List<PhylogenyNode> ext_nodes = p.getExternalNodes();
        int ext_nodes_count = 0;
        int mapped = 0;
        int not_mapped = 0;
        for (final PhylogenyNode ext_node : ext_nodes) {
            ++ext_nodes_count;
            String name = ext_node.getName();
            if (isolate_to_clade_map.containsKey(name)) {
                ++mapped;
                final PropertiesList custom_data;
                if (ext_node.isHasNodeData() && ext_node.getNodeData().isHasProperties()) {
                    custom_data = ext_node.getNodeData().getProperties();
                } else {
                    custom_data = new PropertiesList();
                    ext_node.getNodeData().setProperties(custom_data);
                }
                custom_data.addProperty(new Property(CLADE, isolate_to_clade_map.get(name), "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(SHORTCLADE, isolate_to_shortclade_map.get(name), "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(LEGACYCLADE, isolate_to_legacyclade_map.get(name), "", XSD_STRING, AppliesTo.NODE));

            } else {
                ++not_mapped;
                ForesterUtil.printWarningMessage(PRG_NAME, "No mapping for: " + name);
            }
        }


        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML(p, 0, outfile);
        } catch (final IOException e) {
            System.out.println("\nFailure to write output [" + e.getMessage() + "]\n");
            System.exit(-1);
        }
        System.out.println();
        System.out.println("External nodes                    : " + ext_nodes_count);
        System.out.println("External nodes with clade added   : " + mapped);
        System.out.println("External nodes without clade added: " + not_mapped);
        System.out.println();
        System.out.println();
        System.out.println("[" + PRG_NAME + "] wrote: [" + outfile + "]");
        System.out.println("[" + PRG_NAME + "] OK");
        System.out.println();

    }
}