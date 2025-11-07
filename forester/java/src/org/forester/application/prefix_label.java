
package org.forester.application;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class prefix_label {

    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";

    final static private String PRG_NAME = "prefix_label";
    final static private String PRG_DESC = "add taxonomy based prefixes";
    final static private String PRG_VERSION = "1.0.0";
    final static private String PRG_DATE = "2025.11.07";
    final static private String E_MAIL = "";
    final static private String WWW = "https://github.com/cmzmasek";


    final static private boolean DEBUG = false;

    public static void main(final String[] args) {

        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments(args);
        } catch (IOException e) {
            ForesterUtil.fatalError(PRG_NAME, e.getMessage());
        }

        if (cla.isOptionSet(HELP_OPTION_1) || cla.isOptionSet(HELP_OPTION_2)) {
            printHelp();
            System.exit(0);
        }

        if (args.length != 4) {
            printHelp();
            System.exit(-1);
        }

        final String pattern_str = args[0];
        final File intree = new File(args[1]);
        final File labelsfile = new File(args[2]);
        final File outfile = new File(args[3]);
        Pattern genome_id_pattern = null;
        try {
            genome_id_pattern = Pattern.compile(pattern_str);
        } catch (final Exception e) {
            ForesterUtil.fatalError(PRG_NAME, "Pattern \"" + pattern_str + "\" cannot be compiled [" + e.getMessage() + "]");
        }
        System.out.println("[" + PRG_NAME + "] Pattern: " + genome_id_pattern);

        if (!intree.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + intree + "] does not exist");
        }
        if (!labelsfile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + labelsfile + "] does not exist");
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
            ForesterUtil.fatalError(PRG_NAME, "Could not read \"" + intree + "\" [\" + e.getMessage() + \"]");
        }

        System.out.println("[" + PRG_NAME + "] read in: " + intree);

        String str = null;
        try {
            str = Files.readString(Path.of(labelsfile.toString()));
        } catch (final Exception e) {
            ForesterUtil.fatalError(PRG_NAME, "Could not read \"" + labelsfile + "\" [\" + e.getMessage() + \"]");
        }

        if (ForesterUtil.isEmpty(str)) {
            ForesterUtil.fatalError(PRG_NAME, "\"" + labelsfile + "\" is empty");
        }

        System.out.println("[" + PRG_NAME + "] read in: " + labelsfile);
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(outfile));
        } catch (final Exception e) {
            ForesterUtil.fatalError(PRG_NAME, "Could not open \"" + outfile + "\" for writing  [\" + e.getMessage() + \"]");
        }

        final String[] labels = str.split("\t");

        System.out.println("[" + PRG_NAME + "] split into " + labels.length + " labels");

        try {
            extracted(genome_id_pattern, p, labels, writer);
            writer.close();

        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, e.getMessage());
        }
        System.out.println("[" + PRG_NAME + "] wrote: [" + outfile + "]");
        System.out.println("[" + PRG_NAME + "] OK");
        System.out.println();
    }

    private static void extracted(final Pattern genome_id_pattern, final Phylogeny p, final String[] labels, final BufferedWriter writer) throws IOException {
        for (final String label : labels) {
            writer.write(label);
            writer.write("\t");
        }
        writer.write("\n");

        for (final String label : labels) {
            if (DEBUG) {
                System.out.println(label);
            }
            final Matcher m = genome_id_pattern.matcher(label);
            String genome_id = null;
            if (m.find()) {
                genome_id = m.group(1);
            }
            if (ForesterUtil.isEmpty(genome_id)) {
                writer.write("??????_" + label);
                writer.write("\t");
                System.out.println("WARNING: no match: " + label);
            } else {
                if (DEBUG) {
                    System.out.println("=> " + genome_id);
                }
                final List<PhylogenyNode> x = p.getNodesPartialMatch(genome_id);
                if (x.size() < 1) {
                    System.out.println("WARNING: not found in tree: " + genome_id);
                    writer.write("??????_" + label);
                    writer.write("\t");
                } else if (x.size() > 1) {
                    ForesterUtil.fatalError(PRG_NAME, genome_id + "is not unique");
                } else {
                    PhylogenyNode n = x.get(0);
                    String prefix = findPrefix(n);
                    if (!ForesterUtil.isEmpty(prefix)) {
                        writer.write(prefix + "_" + label);
                        writer.write("\t");
                    } else {
                        System.out.println("WARNING: not prefix found in tree for: " + label);
                        writer.write("??????_" + label);
                        writer.write("\t");
                    }
                }
            }
        }
        writer.write("\n");
    }

    private static String findPrefix(final PhylogenyNode node) {
        PhylogenyNode n = node;
        String prefix = null;
        if (!n.isRoot()) {
            n = n.getParent();
            while (prefix == null) {
                if (n.isHasNodeData() && n.getNodeData().isHasTaxonomy() &&
                        !ForesterUtil.isEmpty(n.getNodeData().getTaxonomy().getScientificName())) {
                    prefix = n.getNodeData().getTaxonomy().getScientificName();
                } else if (!ForesterUtil.isEmpty(n.getName())) {
                    prefix = n.getName();
                }
                if (n.isRoot()) {
                    break;
                } else {
                    n = n.getParent();
                }
            }
        }
        return prefix;
    }

    private static void printHelp() {
        ForesterUtil.printProgramInformation(PRG_NAME,
                PRG_DESC,
                PRG_VERSION,
                PRG_DATE,
                E_MAIL,
                WWW,
                ForesterUtil.getForesterLibraryInformation());
        System.out.println("Usage:   " + PRG_NAME + " <pattern> <intree> <labels> <outfile>");
        System.out.println();
        System.out.println("Example: " + PRG_NAME + " \"\\|(.+?)_prot\" beta_cov.xml labels.tsv out.tsv");
        System.out.println();
    }
}
