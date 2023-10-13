
package org.forester.application;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sequence.BasicSequence;
import org.forester.util.ForesterUtil;

public final class ancestor_seq_x {

    private final static String PRG_NAME = "ancestor_seq_x";
    private static final String PRG_DATE = "2021-11-29";
    private static final String PRG_VERSION = "1.0.2";
    private static final Pattern P1 = Pattern.compile("(\\d+)\\.\\s+(.+?):\\s*(.+)");
    private static final Pattern P2 = Pattern.compile("\\(\\s*(\\d+)\\s*\\.\\s*(\\d+)\\s*\\)");
    private static final String SEQ_NAME = "S";

    public static void main(final String[] args) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        if (args.length != 3) {
            System.out.println(PRG_NAME + ": Wrong number of arguments.\n");
            System.out.println("Usage: " + PRG_NAME + " <in-tree> <most probable sequence file> <out-tree>\n");
            System.out.println("Example: " + PRG_NAME + " tree.xml most_prob_seqs.txt tree_with_anc_seqs.xml\n");
            System.exit(-1);
        }
        final File intree = new File(args[0]);
        final File mega_most_prob_seqs = new File(args[1]);
        final File outtree = new File(args[2]);
        final String error0 = ForesterUtil.isReadableFile(intree);
        if (!ForesterUtil.isEmpty(error0)) {
            ForesterUtil.fatalError(PRG_NAME, error0);
        }
        final String error1 = ForesterUtil.isReadableFile(mega_most_prob_seqs);
        if (!ForesterUtil.isEmpty(error1)) {
            ForesterUtil.fatalError(PRG_NAME, error1);
        }
        final String error2 = ForesterUtil.isWritableFile(outtree);
        if (!ForesterUtil.isEmpty(error2)) {
            ForesterUtil.fatalError(PRG_NAME, error2);
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
        final SortedMap<String, IdSeq> map = new TreeMap<>();
        final SortedMap<String, PhylogenyNode> number_to_node_map = new TreeMap<>();
        readAncestralSeqsFile(mega_most_prob_seqs, p, map, number_to_node_map);
        addSeqsToNodes(map, number_to_node_map);
        int int_nodes_total = 0;
        int int_nodes_with_seq = 0;
        int int_nodes_without_seq = 0;
        int seq_length = -1;
        for (final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if (node.isInternal()) {
                ++int_nodes_total;
                if (node.isHasNodeData() && node.getNodeData().isHasSequence()
                        && !ForesterUtil.isEmpty(node.getNodeData().getSequence().getMolecularSequence())) {
                    ++int_nodes_with_seq;
                    if (seq_length < 0) {
                        seq_length = node.getNodeData().getSequence().getMolecularSequence().length();
                    } else if (seq_length != node.getNodeData().getSequence().getMolecularSequence().length()) {
                        ForesterUtil.fatalError(PRG_NAME, "sequences of unequal length detected");
                    }
                } else {
                    ++int_nodes_without_seq;
                }
            }
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML(p, 0, outtree);
        } catch (final IOException e) {
            System.out.println("\nFailure to write output [" + e.getMessage() + "]\n");
            System.exit(-1);
        }
        System.out.println("Wrote outtree to               : " + outtree);
        System.out.println("Sequence length                : " + seq_length);
        System.out.println("Internal nodes total           : " + int_nodes_total);
        System.out.println("Internal nodes with sequence   : " + int_nodes_with_seq);
        System.out.println("Internal nodes without sequence: " + int_nodes_without_seq + "\n");
    }

    private static void readAncestralSeqsFile(final File mega_most_prob_seqs,
                                              final Phylogeny p,
                                              final SortedMap<String, IdSeq> map,
                                              final SortedMap<String, PhylogenyNode> number_to_node_map) {
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(mega_most_prob_seqs));
            String line;
            boolean lca_is_root = false;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.length() > 0) {
                    final Matcher m1 = P1.matcher(line);
                    if (m1.find()) {
                        final String number = m1.group(1);
                        final String id = m1.group(2);
                        final String seq = m1.group(3);
                        if (map.containsKey(number)) {
                            if (!map.get(number).getId().equals(id)) {
                                ForesterUtil.fatalError(PRG_NAME,
                                        "Error: Ids are not equal: " + map.get(number).getId()
                                                + " != " + id);
                            }
                            map.get(number).getSeq().append(seq);
                        } else {
                            map.put(number, new IdSeq(id, seq));
                        }
                        final Matcher m2 = P2.matcher(id);
                        if (m2.find()) {
                            final String number_1 = m2.group(1);
                            final String number_2 = m2.group(2);
                            final PhylogenyNode node_1 = number_to_node_map.get(number_1);
                            final PhylogenyNode node_2 = number_to_node_map.get(number_2);
                            if (node_1 == null) {
                                System.out.println(number_to_node_map);
                                System.out.println(line);
                                ForesterUtil.fatalError(PRG_NAME, "no node for " + number_1 + " found in map");
                            }
                            if (node_2 == null) {
                                System.out.println(number_to_node_map);
                                System.out.println(line);
                                ForesterUtil.fatalError(PRG_NAME, "no node for " + number_2 + " found in map");
                            }
                            final PhylogenyNode lca = PhylogenyMethods.calculateLCA(node_1, node_2);
                            if (lca == null) {
                                ForesterUtil.fatalError(PRG_NAME, "LCA is null");
                            }
                            if (lca.isRoot()) {
                                if (lca_is_root) {
                                    ForesterUtil
                                            .fatalError(PRG_NAME,
                                                    "More than one LCAs are root, trees probably don't match");
                                }
                                lca_is_root = true;
                            }
                            number_to_node_map.put(number, lca);
                        } else {
                            lca_is_root = false;
                            final PhylogenyNode node = findNodeBySecAcc(p, id);
                            if (node != null) {
                                number_to_node_map.put(number, node);
                            }
                        }
                    }
                }
            }
            reader.close();
        } catch (final IOException e) {
            e.printStackTrace();
        }
    }

    private static void addSeqsToNodes(final SortedMap<String, IdSeq> map,
                                       final SortedMap<String, PhylogenyNode> number_to_node_map) {
        for (final Entry<String, IdSeq> entry : map.entrySet()) {
            final String number = entry.getKey();
            final String id = entry.getValue().getId();
            final String seq = entry.getValue().getSeq().toString();
            //System.out.println( number + ":" + id + ": " + seq );
            final PhylogenyNode node = number_to_node_map.get(number);
            if (node == null) {
                System.out.println("node found node: " + number);
                //System.exit( -1 );
            } else {
                if (node.isHasNodeData() && node.getNodeData().isHasSequence()
                        && !node.getNodeData().getSequence().isEmpty()) {
                    if ((node.getNodeData().getSequence().getAccession() != null)
                            && !ForesterUtil.isEmpty(node.getNodeData().getSequence().getAccession().getValue())) {
                        final String myacc = node.getNodeData().getSequence().getAccession().getValue();
                        if (!myacc.equals(id)) {
                            System.out.println("ERROR");
                            System.exit(-1);
                        }
                        node.getNodeData().getSequence().setMolecularSequence(seq);
                        node.getNodeData().getSequence().setMolecularSequenceAligned(true);
                    }
                } else {
                    node.getNodeData().addSequence(new Sequence(BasicSequence.createAaSequence(SEQ_NAME, seq)));
                    node.getNodeData().getSequence().setMolecularSequenceAligned(true);
                }
            }
        }
    }

    private static PhylogenyNode findNodeBySecAcc(final Phylogeny p, final String acc) {
        for (final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if (node.isHasNodeData() && node.getNodeData().isHasSequence()
                    && !node.getNodeData().getSequence().isEmpty()) {
                if ((node.getNodeData().getSequence().getAccession() != null)
                        && !ForesterUtil.isEmpty(node.getNodeData().getSequence().getAccession().getValue())) {
                    final String myacc = node.getNodeData().getSequence().getAccession().getValue();
                    if (acc.equals(myacc)) {
                        return node;
                    }
                }
            }
        }
        return null;
    }

    final static class IdSeq {

        final String _id;
        final StringBuilder _seq;

        IdSeq(final String id, final String seq) {
            _id = id;
            _seq = new StringBuilder(seq);
        }

        public String getId() {
            return _id;
        }

        public StringBuilder getSeq() {
            return _seq;
        }
    }
}
