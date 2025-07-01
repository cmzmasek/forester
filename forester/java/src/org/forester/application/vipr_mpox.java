
package org.forester.application;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class vipr_mpox {

    private static final String PRG_DATE = "2025-06-13";
    private static final String PRG_VERSION = "1.0.0";
    private static final String PRG_NAME = "vipr_mpox";
    private static final String XSD_STRING = "xsd:string";
    private static final String CLADE = "vipr:Clade";
    private static final String COUNTRY = "vipr:Country";
    private static final String YEAR = "vipr:Year";


    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        System.out.println();
        if ((args.length != 3)) {
            System.out.println(PRG_NAME + ": Wrong number of arguments.\n");
            System.out.println("Usage: " + PRG_NAME + " <in-tree> <out-tree> <description> \n");
            System.exit(-1);
        }
        final File intree = new File(args[0]);
        final File outtree = new File(args[1]);
        final String desc = args[2].trim();
        final String e0 = ForesterUtil.isWritableFile(outtree);
        if (!ForesterUtil.isEmpty(e0)) {
            ForesterUtil.fatalError(PRG_NAME, e0);
        }
        final String e1 = ForesterUtil.isReadableFile(intree);
        if (!ForesterUtil.isEmpty(e1)) {
            ForesterUtil.fatalError(PRG_NAME, e1);
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


        p.setName(desc);

        reRoot(p, "KJ642618", "KJ642615");
        p.setRerootable(false);

        assignClade(p, "I", "KJ642618", "PQ221902");
        assignClade(p, "I.a", "KJ642618", "OQ729808");
        assignClade(p, "I.b", "PQ221902", "PP601218");
        assignClade(p, "II", "KJ642615", "KJ642614");
        assignClade(p, "II.a", "MN346692", "KJ642614");
        assignClade(p, "II.b", "KJ642617", "OP535319");
        assignClade(p, "II.b.A", "PP852954", "OR499970");
        assignClade(p, "II.b.A.1", "OR499970", "MT903341");
        assignClade(p, "II.b.A.2", "PP852972", "PP852967");

        for (final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            procNode(node);
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML(p, 0, outtree);
        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, "failed to write to [" + outtree + "]: " + e.getLocalizedMessage());
        }

        System.out.println();
        System.out.println("Wrote outtree to       : " + outtree);
        System.out.println();
    }

    private static void assignClade(final
                                    Phylogeny p, final String clade, final String node1, final String node2) {
        final PhylogenyNode lca = PhylogenyMethods.calculateLCA(p.getNodePartialMatch(node1), p.getNodePartialMatch(node2));
        lca.setName(clade);
    }

    private static void reRoot(final Phylogeny p, final String node1, final String node2) {
        final PhylogenyNode lca = PhylogenyMethods.calculateLCA(p.getNodePartialMatch(node1), p.getNodePartialMatch(node2));
        p.reRoot(lca);
    }

    private static void procNode(final PhylogenyNode node) {
        if (node.isExternal()) {
            String clade = "";
            PhylogenyNode n = node;
            while (!n.isRoot()) {
                n = n.getParent();
                if (n.isInternal() && !ForesterUtil.isEmpty(n.getName())) {
                    clade = n.getName();
                    break;
                }
            }
            String name = node.getName();
            if (name.endsWith("|")) {
                name = name.substring(0, name.length() - 1);
            }
            if (node.getNodeData().isHasProperties()) {
                if (name.length() < 15) {
                    List<Property> c = node.getNodeData().getProperties().getPropertiesWithGivenReferencePrefix(COUNTRY);
                    if (c != null && c.size() == 1) {
                        name = name + "|" + c.get(0).getValue().replaceAll("\\s+", "_");
                    }
                    List<Property> y = node.getNodeData().getProperties().getPropertiesWithGivenReferencePrefix(YEAR);
                    if (y != null && y.size() == 1) {
                        name = name + "|" + y.get(0).getValue();
                    }
                }
            }
            node.setName(name + "/" + clade);

            if (!node.getNodeData().isHasProperties()) {
                node.getNodeData().setProperties(new PropertiesList());
            }
            node.getNodeData().getProperties().addProperty(new Property(CLADE, clade, "", XSD_STRING, Property.AppliesTo.NODE));
        }
    }

}
