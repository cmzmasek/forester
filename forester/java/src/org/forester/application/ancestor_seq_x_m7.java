
package org.forester.application;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
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

public final class ancestor_seq_x_m7 {

    private final static String  PRG_NAME    = "ancestor_seq_x_m7";
    private static final String  PRG_DATE    = "2021-06-01";
    private static final String  PRG_VERSION = "1.0.0";
    private static final Pattern P1          = Pattern.compile( "(\\d+)\\.\\s+?([A-Za-z0-9_]+?):(.+)" );
    private static final Pattern P2          = Pattern.compile( "(\\d+)\\.\\s+?\\(\\s*(\\d+)\\s*\\.\\s*(\\d+)\\s*\\)" );
    private static final String  SEQ_NAME    = "S";
    public static void main( final String[] args ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        if ( args.length != 3 ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree> <most probable sequence file> <out-tree>\n" );
            System.out.println( "Example: " + PRG_NAME + " tree.xml most_prob_seqs.txt tree_with_anc_seqs.xml\n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final File mega_most_prob_seqs = new File( args[ 1 ] );
        final File outtree = new File( args[ 2 ] );
        final String error0 = ForesterUtil.isReadableFile( intree );
        if ( !ForesterUtil.isEmpty( error0 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error0 );
        }
        final String error1 = ForesterUtil.isReadableFile( mega_most_prob_seqs );
        if ( !ForesterUtil.isEmpty( error1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error1 );
        }
        final String error2 = ForesterUtil.isWritableFile( outtree );
        if ( !ForesterUtil.isEmpty( error2 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error2 );
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intree, true );
            p = factory.create( intree, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + intree + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        final SortedMap<String, IdSeq> map = new TreeMap<>();
        final SortedMap<String, PhylogenyNode> number_to_node_map = new TreeMap<>();
        readAncestralSeqsFileF( mega_most_prob_seqs, p, map, number_to_node_map );
        System.exit( 0 );
        addSeqsToNodes( map, number_to_node_map );
        int int_nodes_total = 0;
        int int_nodes_with_seq = 0;
        int int_nodes_without_seq = 0;
        int seq_length = -1;
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isInternal() ) {
                ++int_nodes_total;
                if ( node.isHasNodeData() && node.getNodeData().isHasSequence()
                        && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() ) ) {
                    ++int_nodes_with_seq;
                    if ( seq_length < 0 ) {
                        seq_length = node.getNodeData().getSequence().getMolecularSequence().length();
                    }
                    else if ( seq_length != node.getNodeData().getSequence().getMolecularSequence().length() ) {
                        ForesterUtil.fatalError( PRG_NAME, "sequences of unequal length detected" );
                    }
                }
                else {
                    ++int_nodes_without_seq;
                }
            }
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( p, 0, outtree );
        }
        catch ( final IOException e ) {
            System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        System.out.println( "Wrote outtree to               : " + outtree );
        System.out.println( "Sequence length                : " + seq_length );
        System.out.println( "Internal nodes total           : " + int_nodes_total );
        System.out.println( "Internal nodes with sequence   : " + int_nodes_with_seq );
        System.out.println( "Internal nodes without sequence: " + int_nodes_without_seq + "\n" );
    }

    private static void readAncestralSeqsFileF( final File mega_most_prob_seqs,
                                                final Phylogeny p,
                                                final SortedMap<String, IdSeq> map,
                                                final SortedMap<String, PhylogenyNode> number_to_node_map ) {
        BufferedReader reader;
        final SortedMap<String, String> number_to_seqacc_map = new TreeMap<>();
        final SortedMap<String, NumberPair> number_to_numberpair_map = new TreeMap<>();
        try {
            reader = new BufferedReader( new FileReader( mega_most_prob_seqs ) );
            String line;
            boolean done = false;
            while ( ( line = reader.readLine() ) != null ) {
                if ( !done ) {
                    line.trim();
                    if ( line.length() > 0 ) {
                        final Matcher m1 = P1.matcher( line );
                        final Matcher m2 = P2.matcher( line );
                        if ( m1.find() ) {
                            System.out.println( "P1 matches: " + line );
                            number_to_seqacc_map.put( m1.group( 1 ), m1.group( 2 ) );
                        }
                        if ( m2.find() ) {
                            System.out.println( "P2 matches: " + line );
                            number_to_numberpair_map.put( m2.group( 1 ),
                                                          new NumberPair( m2.group( 2 ), m2.group( 3 ) ) );
                        }
                        else if ( ( number_to_numberpair_map.size() > 0 ) ) {
                            done = true;
                        }
                    }
                }
            }
            reader.close();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        //////
        final Iterator<Entry<String, String>> it = number_to_seqacc_map.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<String, String> e = it.next();
            final String number = e.getKey();
            final String seqacc = e.getValue();
            final PhylogenyNode node = findNodeBySecAcc( p, seqacc );
            if ( node == null ) {
                ForesterUtil.fatalError( PRG_NAME, "node with sequence accession " + seqacc + " not found" );
            }
            number_to_node_map.put( number, node );
        }
        /////
        /////
      
        boolean not_done = true;
        while ( not_done ) {
            boolean all_found = true;
            final Iterator<Entry<String, NumberPair>> it2 = number_to_numberpair_map.entrySet().iterator();
            System.out.println(  );
            while ( it2.hasNext() ) {
                final Map.Entry<String, NumberPair> e = it2.next();
                final String number = e.getKey();
                final NumberPair nn = e.getValue();
                final String nu1 = nn.getNumber1();
                final String nu2 = nn.getNumber2();
                if ( number_to_node_map.containsKey( nu1 ) && number_to_node_map.containsKey( nu2 ) ) {
                    final PhylogenyNode node1 = number_to_node_map.get( nu1 );
                    final PhylogenyNode node2 = number_to_node_map.get( nu2 );
                    final PhylogenyNode lca = PhylogenyMethods.calculateLCA( node1, node2 );
                    number_to_node_map.put( number, lca);
                    System.out.println( number + " => " + lca );
                }
                else {
                    all_found = false;
                }
            }
            if (all_found) {
                not_done = false;
            }
        }
        
        System.out.println( number_to_node_map );
        /////
    }

    private static void readAncestralSeqsFile( final File mega_most_prob_seqs,
                                               final Phylogeny p,
                                               final SortedMap<String, IdSeq> map,
                                               final SortedMap<String, PhylogenyNode> number_to_node_map ) {
        BufferedReader reader;
        try {
            reader = new BufferedReader( new FileReader( mega_most_prob_seqs ) );
            String line;
            while ( ( line = reader.readLine() ) != null ) {
                line.trim();
                if ( line.length() > 0 ) {
                    final Matcher m1 = P1.matcher( line );
                    if ( m1.find() ) {
                        final String number = m1.group( 1 );
                        final String id = m1.group( 2 );
                        final String seq = m1.group( 3 );
                        if ( map.containsKey( number ) ) {
                            if ( !map.get( number ).getId().equals( id ) ) {
                                System.out.println( "Error: Ids are not equal: " + map.get( number ).getId() + " != "
                                        + id );
                                System.exit( -1 );
                            }
                            map.get( number ).getSeq().append( seq );
                        }
                        else {
                            map.put( number, new IdSeq( id, seq ) );
                        }
                        final Matcher m2 = P2.matcher( id );
                        if ( m2.find() ) {
                            final String number_1 = m2.group( 1 );
                            final String number_2 = m2.group( 2 );
                            final PhylogenyNode node_1 = number_to_node_map.get( number_1 );
                            final PhylogenyNode node_2 = number_to_node_map.get( number_2 );
                            if ( ( node_1 != null ) && ( node_2 != null ) ) {
                                final PhylogenyNode lca = PhylogenyMethods.calculateLCA( node_1, node_2 );
                                number_to_node_map.put( number, lca );
                            }
                            else {
                            }
                        }
                        else {
                            final PhylogenyNode node = findNodeBySecAcc( p, id );
                            if ( node != null ) {
                                number_to_node_map.put( number, node );
                            }
                        }
                    }
                }
            }
            reader.close();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }

    private static void addSeqsToNodes( final SortedMap<String, IdSeq> map,
                                        final SortedMap<String, PhylogenyNode> number_to_node_map ) {
        for( final Entry<String, IdSeq> entry : map.entrySet() ) {
            final String number = entry.getKey();
            final String id = entry.getValue().getId();
            final String seq = entry.getValue().getSeq().toString();
            //System.out.println( number + ":" + id + ": " + seq );
            final PhylogenyNode node = number_to_node_map.get( number );
            if ( node == null ) {
                System.out.println( "node found node: " + number );
                //System.exit( -1 );
            }
            else {
                if ( node.isHasNodeData() && node.getNodeData().isHasSequence()
                        && !node.getNodeData().getSequence().isEmpty() ) {
                    if ( ( node.getNodeData().getSequence().getAccession() != null )
                            && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() ) ) {
                        final String myacc = node.getNodeData().getSequence().getAccession().getValue();
                        if ( !myacc.equals( id ) ) {
                            System.out.println( "ERROR" );
                            System.exit( -1 );
                        }
                        node.getNodeData().getSequence().setMolecularSequence( seq );
                        node.getNodeData().getSequence().setMolecularSequenceAligned( true );
                    }
                }
                else {
                    node.getNodeData().addSequence( new Sequence( BasicSequence.createAaSequence( SEQ_NAME, seq ) ) );
                    node.getNodeData().getSequence().setMolecularSequenceAligned( true );
                }
            }
        }
    }

    private static PhylogenyNode findNodeBySecAcc( final Phylogeny p, final String acc ) {
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isHasNodeData() && node.getNodeData().isHasSequence()
                    && !node.getNodeData().getSequence().isEmpty() ) {
                if ( ( node.getNodeData().getSequence().getAccession() != null )
                        && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() ) ) {
                    final String myacc = node.getNodeData().getSequence().getAccession().getValue();
                    if ( acc.equals( myacc ) ) {
                        return node;
                    }
                }
            }
        }
        return null;
    }
    final static class IdSeq {

        final String        _id;
        final StringBuilder _seq;
        IdSeq( final String id, final String seq ) {
            _id = id;
            _seq = new StringBuilder( seq );
        }

        public String getId() {
            return _id;
        }

        public StringBuilder getSeq() {
            return _seq;
        }
    }

    final static class NumberPair {

        final String _n1;
        final String _n2;
        NumberPair( final String n1, final String n2 ) {
            _n1 = n1;
            _n2 = n2;
        }

        public String getNumber1() {
            return _n1;
        }

        public String getNumber2() {
            return _n2;
        }
    }
}
