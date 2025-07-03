package org.forester.application;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

public class treesort_reformate {

    private static final String PRG_DATE = "2025-07-03";
    private static final String PRG_VERSION = "0.0.1";
    private static final String PRG_NAME = "treesort_reformate";

    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        System.out.println();
        if ((args.length != 2)) {
            System.out.println(PRG_NAME + ": Wrong number of arguments.\n");
            System.out.println("Usage: " + PRG_NAME + " <in-tree> <out-file> \n");
            System.exit(-1);
        }
        final File intree = new File(args[0]);
        final File outfile = new File(args[1]);
        final String e0 = ForesterUtil.isWritableFile(outfile);
        if (!ForesterUtil.isEmpty(e0)) {
            ForesterUtil.fatalError(PRG_NAME, e0);
        }
        final String e1 = ForesterUtil.isReadableFile(intree);
        if (!ForesterUtil.isEmpty(e1)) {
            ForesterUtil.fatalError(PRG_NAME, e1);
        }
        execute(intree, outfile);
    }

    private static void execute(File intree, File outfile) {
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            p = factory.create(intree, new NHXParser())[0];
        } catch (final Exception e) {
            ForesterUtil.fatalError(PRG_NAME, "Could not read \"" + intree + "\" [" + e.getMessage() + "]");
        }
        ForesterUtil.programMessage(PRG_NAME, "Successfully read in tree with " + p.getNumberOfExternalNodes() + " external nodes");
        Writer w = null;
        int counter = 0;
        int counter_all = 0;
        try {
            w = ForesterUtil.createBufferedWriter(outfile);
            for (final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
                if (node.isExternal()) {
                    ++counter_all;
                    final PropertiesList properties = node.getNodeData().getProperties();
                    if (properties != null && properties.size() > 0) {
                        List<Property> c = properties.getPropertiesWithGivenReferencePrefix(ForesterConstants.NH_COMMENT);
                        if (c != null && c.size() > 0) {
                            w.write(node.getName());
                            w.write(", ");
                            w.write(c.get(0).getValue());
                            w.write('\n');
                            ++counter;
                        }
                    }
                }
            }
            w.close();
        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, e.getLocalizedMessage());
        }
        ForesterUtil.programMessage(PRG_NAME, "Wrote output for " + counter + "/" + counter_all + " external nodes to " + outfile);
    }
}