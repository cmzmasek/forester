
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

public class tree_label_interpreter {

    private static final String PRG_DATE = "2024-04-15";
    private static final String PRG_VERSION = "1.0.0";
    private static final String PRG_NAME = "tree_label_interpreter";
    private static final String XSD_STRING = "xsd:string";

    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        System.out.println();
        if ((args.length != 3)) {
            System.out.println(PRG_NAME + ": Wrong number of arguments.\n");
            System.out.println("Usage  : " + PRG_NAME + " <legend> <in-tree> <out-tree>\n");
            System.out.println("Example: " + PRG_NAME + " labels.txt p1.xml p2.xml\n");
            System.exit(-1);
        }
        final File legend = new File(args[0]);
        final File intree = new File(args[1]);
        final File outtree = new File(args[2]);

        final String e0 = ForesterUtil.isWritableFile(outtree);
        if (!ForesterUtil.isEmpty(e0)) {
            ForesterUtil.fatalError(PRG_NAME, e0);
        }
        final String e1 = ForesterUtil.isReadableFile(intree);
        if (!ForesterUtil.isEmpty(e1)) {
            ForesterUtil.fatalError(PRG_NAME, e1);
        }
        final String e2 = ForesterUtil.isReadableFile(legend);
        if (!ForesterUtil.isEmpty(e2)) {
            ForesterUtil.fatalError(PRG_NAME, e2);
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(intree, true);
            p = factory.create(intree, pp)[0];
        } catch (final Exception e) {
            ForesterUtil.fatalError(PRG_NAME, "Could not read \"" + intree + "\" [" + e.getMessage() + "]");
        }
        ForesterUtil
                .programMessage(PRG_NAME,
                        "Successfully read in tree with " + p.getNumberOfExternalNodes() + " external nodes");
        BasicTable<String> legend_table = null;
        try {
            legend_table = BasicTableParser.parse(legend, '\t');
        } catch (final Exception e) {
            ForesterUtil.fatalError(PRG_NAME, "failed to read \"" + legend + "\" [" + e.getMessage() + "]");
        }
        if (legend_table.getNumberOfColumns() != 2) {
            ForesterUtil.fatalError(PRG_NAME,
                    "table has " + legend_table.getNumberOfColumns() + " columns, expected two");
        }
        if (legend_table.getNumberOfRows() < 1) {
            ForesterUtil.fatalError(PRG_NAME, "table has no rows (i.e. is empty)");
        }
        ForesterUtil.programMessage(PRG_NAME,
                "Successfully read in table with " + legend_table.getNumberOfColumns() + " columns and "
                        + legend_table.getNumberOfRows() + " rows");


        final SortedMap<String, String> legend_map = legend_table.getColumnsAsMap(0, 1);

        final String sep;
        if (legend_map.containsKey("separator")) {
            sep = legend_map.get("separator");
        } else {
            sep = "/";
        }

        final SortedMap<String, Pattern> property_to_regex_map = new TreeMap<String, Pattern>();
        for (int i = 0; i < legend_table.getNumberOfRows(); ++i) {
            if (legend_table.getValueAsString(0, i).startsWith("/")) {
                String pattern_str = legend_table.getValueAsString(0, i);

                pattern_str = pattern_str.substring(1,pattern_str.length() - 1);

                property_to_regex_map.put(legend_table.getValueAsString(1, i), Pattern.compile(pattern_str));

                System.out.println( );

                ForesterUtil.programMessage(PRG_NAME, "Pattern " + pattern_str + ": " + legend_table.getValueAsString(1, i));

            }
        }

        ForesterUtil.programMessage(PRG_NAME, "Separator: " + sep);

        int properties_added = 0;
        for (final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();

            if (!ForesterUtil.isEmpty(node.getName())) {
                final String name = node.getName();

                String[] elements = name.split(sep);
                for (int i = 0; i < elements.length; ++i) {
                    final String i_str = String.valueOf(i);
                    if (legend_map.containsKey(i_str)) {
                        final String element = elements[i].trim().replaceAll("\\s+", " ");
                        final String ref = legend_map.get(i_str);
                        if (ref.indexOf(':') < 0) {
                            ForesterUtil
                                    .fatalError(PRG_NAME,
                                            "Property reference " + "\"" + ref + "\" lacks required :");
                        }
                        PropertiesList custom_data = node.getNodeData().getProperties();
                        if (custom_data == null) {
                            custom_data = new PropertiesList();
                            node.getNodeData().setProperties(custom_data);
                        }
                        custom_data.addProperty(new Property(ref, element, "", XSD_STRING, AppliesTo.NODE));
                        ++properties_added;

                    }
                }

                for (Map.Entry<String, Pattern> entry : property_to_regex_map.entrySet()) {
                    final Pattern pa = entry.getValue();
                    final String ref = entry.getKey();
                    final Matcher matcher = pa.matcher(name);
                    if (matcher.find() && matcher.groupCount() == 1) {
                        final String g = matcher.group(1);
                        System.out.println(g);
                        PropertiesList custom_data = node.getNodeData().getProperties();
                        if (custom_data == null) {
                            custom_data = new PropertiesList();
                            node.getNodeData().setProperties(custom_data);
                        }
                        custom_data.addProperty(new Property(ref, g, "", XSD_STRING, AppliesTo.NODE));
                        ++properties_added;
                    }
                }

            }
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML(p, 0, outtree);
        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, "failed to write to [" + outtree + "]: " + e.getLocalizedMessage());
        }
        System.out.println();
        System.out.println("Sum of properties added: " + properties_added);
        System.out.println("Wrote outtree to       : " + outtree);
        System.out.println();
    }
}
