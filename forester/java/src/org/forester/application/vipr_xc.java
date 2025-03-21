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

    private static final String PRG_DATE = "2024-06-11";
    private static final String PRG_VERSION = "1.0.2";
    private final static String PRG_NAME = "vipr_xc";
    private static final String XSD_STRING = "xsd:string";
    private static final String H5_CLADE = "vipr:H5_clade";


    private final static Pattern PATTERN_GB = Pattern.compile("\\|([A-Z][A-Z0-9.]{4,10}?)(_|$)");

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


        /*final SortedMap<String, String> isolate_to_clade_map = new TreeMap<String, String>();
        for (int row = 0; row < isolate_to_clade.getNumberOfRows(); ++row) {
            final String keys = isolate_to_clade.getValue(0, row);
            final String value = isolate_to_clade.getValue(1, row);
            final Matcher mg = PATTERN_GB.matcher(keys);
            while(mg.find()) {
                final String key_acc = mg.group(1);
                if ((key_acc != null) && (value != null)) {
                    if (isolate_to_clade_map.containsKey(key_acc)) {
                        throw new IllegalArgumentException("attempt to use non-unique table value as key [" + key_acc + "]");
                    }
                    isolate_to_clade_map.put(key_acc, value);
                }
            }
        }*/

        final List<PhylogenyNode> ext_nodes = p.getExternalNodes();
        int ext_nodes_count = 0;
        int mapped = 0;
        int not_mapped = 0;
        for (final PhylogenyNode ext_node : ext_nodes) {
            ++ext_nodes_count;
            String name = ext_node.getName();
            final String x[] = name.split("\\|");
            //final String name_acc = x[x.length - 1];
            name = x[0];
            if (!ForesterUtil.isEmpty(name)) {
                name = name.replaceAll(" ", "_").replace(",", "").replaceAll("'", "");
                boolean could_map = false;
                F:
                for (int row = 0; row < isolate_to_clade.getNumberOfRows(); ++row) {
                    final String names = isolate_to_clade.getValue(0, row);
                    final String clade = isolate_to_clade.getValue(1, row);
                    if (names.indexOf(name) > -1 && clade.indexOf("cannot") < 0) {
                        PropertiesList custom_data = ext_node.getNodeData().getProperties();
                        custom_data.addProperty(new Property(H5_CLADE, clade, "", XSD_STRING, AppliesTo.NODE));
                        could_map = true;
                        break F;
                    }
                }
                if (could_map) {
                    ++mapped;
                } else {
                    System.out.println( "Could not map: " + name);
                    ++not_mapped;
                }


                /*if (isolate_to_clade_map.containsKey(name_acc)) {
                    final String clade = isolate_to_clade_map.get(name_acc);
                    if (clade.length() > 0 && clade.indexOf("cannot") < 0) {
                        PropertiesList custom_data = ext_node.getNodeData().getProperties();
                        custom_data.addProperty(new Property(H5_CLADE, clade, "", XSD_STRING, AppliesTo.NODE));
                        ++mapped;
                    } else {
                        ++not_mapped;
                    }
                } else {
                    System.out.println("Error: not found in map: " + name);
                    System.exit(-1);
                }*/
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