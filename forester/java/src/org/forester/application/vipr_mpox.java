
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

    private static final String PRG_DATE = "2024-12-02";
    private static final String PRG_VERSION = "0.0.1";
    private static final String PRG_NAME = "vipr_mpox";
    private static final String XSD_STRING = "xsd:string";
    private static final String CLADE = "vipr:Clade";


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
        final int properties_added = 0;

        p.setName(desc);

        reRoot(p, "KJ642618", "KJ642615");
        p.setRerootable(false);

        assignClade(p, "I", "KJ642618", "PQ305773");
        assignClade(p, "I.a", "KJ642618", "OQ729808");
        assignClade(p, "I.b", "PQ305773", "PP601218");
        assignClade(p, "II", "KJ642615", "KJ642614");
        assignClade(p, "II.a", "MN346692", "KJ642614");
        assignClade(p, "II.b", "KJ642617", "PQ207094");
        assignClade(p, "II.b.A", "PP852954", "PQ207094");
        assignClade(p, "II.b.A.1", "PQ207094", "MT903341");
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
        System.out.println("Sum of properties added: " + properties_added);
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
            node.setName(node.getName() + "|" + clade);

            if (!node.getNodeData().isHasProperties()) {
                node.getNodeData().setProperties(new PropertiesList());
            }
            node.getNodeData().getProperties().addProperty(new Property(CLADE, clade, "", XSD_STRING, Property.AppliesTo.NODE));
        }
    }

}
