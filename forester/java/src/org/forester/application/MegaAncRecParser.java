
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

public final class MegaAncRecParser {

    private final static String  PRG_NAME = "MegaAncRecParser";
    private static final Pattern P1       = Pattern.compile( "(\\d+)\\.\\s+?(.+?):(.+)" );
    private static final Pattern P2       = Pattern.compile( "\\(\\s*(\\d+)\\s*\\.\\s*(\\d+)\\s*\\)" );
    public static void main( final String[] args ) {
        // final File intree = new File( args[ 0 ] );
        final File intree = new File( "/home/lambda/Dropbox/WORK/JCVI/ANC_REC/SARS2_4_14_21_29400_09999_pango_3_MAFFT_05_GTR_fastme_fme_mp_pdvxvm.xml" );
        if ( !intree.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + intree + "] does not exist" );
        }
        ///////////////////////
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
        /////////////////////
        final SortedMap<String, IdSeq> map = new TreeMap<>();
        final SortedMap<String, PhylogenyNode> number_to_node_map = new TreeMap<>();
        BufferedReader reader;
        try {
            reader = new BufferedReader( new FileReader( "/home/lambda/Dropbox/WORK/JCVI/ANC_REC/most_prob_seqs.txt" ) );
            String line;
            while ( ( line = reader.readLine() ) != null ) {
                // System.out.println( line );
                line.trim();
                if ( line.length() > 0 ) {
                    final Matcher m = P1.matcher( line );
                    if ( m.find() ) {
                        final String number = m.group( 1 );
                        final String id = m.group( 2 );
                        final String seq = m.group( 3 );
                        //  System.out.print( number );
                        //  System.out.print( " --> " );
                        //  System.out.print( id );
                        //  System.out.print( " --> " );
                        //  System.out.println( seq );
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
                        ////
                        final Matcher m4 = P2.matcher( id );
                        if ( m4.find() ) {
                            final String number_1 = m4.group( 1 );
                            final String number_2 = m4.group( 2 );
                            final PhylogenyNode node_1 = number_to_node_map.get( number_1 );
                            final PhylogenyNode node_2 = number_to_node_map.get( number_2 );
                            final PhylogenyNode lca = PhylogenyMethods.calculateLCA( node_1, node_2 );
                            number_to_node_map.put( number, lca );
                        }
                        else {
                            final PhylogenyNode node = findNodeBySecAcc( p, id );
                            if ( node != null ) {
                                number_to_node_map.put( number, node );
                            }
                        }
                        /////
                    }
                }
            }
            reader.close();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        final SortedMap<String, String> number_to_acc_map = new TreeMap<>();
        for( final Entry<String, IdSeq> entry : map.entrySet() ) {
            final String number = entry.getKey();
            final String id = entry.getValue().getId();
            final String seq = entry.getValue().getSeq().toString();
            System.out.println( number + ":" + id + ": " + seq );
            final PhylogenyNode node = number_to_node_map.get( number );
            if ( node == null ) {
                System.out.println( "------------------------------------------------> " + number );
                //System.exit( -1 );
            }
            else {
                if ( node.isHasNodeData() && node.getNodeData().isHasSequence()
                        && !node.getNodeData().getSequence().isEmpty() ) {
                    if ( ( node.getNodeData().getSequence().getAccession() != null )
                            && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() ) ) {
                        final String myacc = node.getNodeData().getSequence().getAccession().getValue();
                        if ( !myacc.equals( id ) ) {
                            System.out.println( "ERROR " );
                            System.exit( -1 );
                        }
                        node.getNodeData().getSequence().setMolecularSequence( seq );
                    }
                }
                else {
                    node.getNodeData().addSequence( new Sequence( BasicSequence.createAaSequence( id, seq ) ) );
                }
            }
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( p, 0, new File( "/home/lambda/Dropbox/WORK/JCVI/ANC_REC/MegaAncRecPars_out.xml" ) );
        }
        catch ( final IOException e ) {
            System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
            System.exit( -1 );
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
}
