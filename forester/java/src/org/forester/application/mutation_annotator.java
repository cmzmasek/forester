
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.SortedMap;
import java.util.TreeMap;

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
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class mutation_annotator {

    private static final String PRG_NAME   = "mutation_annotator";
    private static final String XSD_STRING = "xsd:string";
    //TODO add option to delete seqs afterwards...
    public static void main( final String args[] ) {
        if ( ( args.length != 4 ) && ( args.length != 6 ) ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME
                    + " <in-tree> <out-tree> <prefix> <property reference> [name of node with reference seq] [property reference for mutations vs reference seq]\n" );
            System.out
                    .println( "Examples : " + PRG_NAME + " tree_with_anc_seqs.xml outtree.xml S aptx:branch_event\n" );
            System.out.println( "         : " + PRG_NAME
                    + " tree_with_anc_seqs.xml outtree.xml S aptx:branch_event whu1 vipr:blah\n" );
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
        final String prefix = args[ 2 ];
        final String property_ref = args[ 3 ];
        String reference_seqe_node_name = null;
        String mut_vs_reference_seq_property_ref = null;
        if ( args.length == 6 ) {
            reference_seqe_node_name = args[ 4 ];
            mut_vs_reference_seq_property_ref = args[ 5 ];
        }
        checkForOutputFileWriteability( outfile );
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( infile, true );
            p = factory.create( infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "Could not read \"" + infile + "\" [" + e.getMessage() + "]" );
        }
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isHasNodeData() && ( node.getNodeData().getSequence() != null )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() ) ) {
            }
            else {
                System.out.println( "no sequence for: " + node );
                if ( node.isInternal()) {
                    System.out.println( "    child 1: " + node.getChildNode1() );
                    System.out.println( "    child 2: " + node.getChildNode2() );
                }
            }
        }
        int total_nodes = 0;
        int updated_nodes = 0;
        int total_mutations = 0;
        int branches_without_seqs = 0;
        int branches_with_seqs = 0;
        final BasicDescriptiveStatistics stats = new BasicDescriptiveStatistics();
        int seq_length = -1;
        final SortedMap<String, Integer> branch_mutations = new TreeMap<>();
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( !node.isRoot() ) {
                ++total_nodes;
                final PhylogenyNode parent = node.getParent();
                if ( node.isHasNodeData() && parent.isHasNodeData() && ( node.getNodeData().getSequence() != null )
                        && ( parent.getNodeData().getSequence() != null )
                        && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() )
                        && !ForesterUtil.isEmpty( parent.getNodeData().getSequence().getMolecularSequence() ) ) {
                    ++branches_with_seqs;
                    final String child_seq = node.getNodeData().getSequence().getMolecularSequence().trim()
                            .toUpperCase();
                    final String parent_seq = parent.getNodeData().getSequence().getMolecularSequence().trim()
                            .toUpperCase();
                    if ( seq_length < 0 ) {
                        seq_length = parent_seq.length();
                    }
                    if ( child_seq.length() != seq_length ) {
                        ForesterUtil.fatalError( PRG_NAME,
                                                 "node with unequal sequence length found: length found is "
                                                         + child_seq.length() + " expected " + seq_length );
                    }
                    if ( parent_seq.length() != seq_length ) {
                        ForesterUtil.fatalError( PRG_NAME,
                                                 "node with unequal sequence length found: length found is "
                                                         + parent_seq.length() + " expected " + seq_length );
                    }
                    int mutations = 0;
                    for( int x = 0; x < seq_length; ++x ) {
                        if ( parent_seq.charAt( x ) != child_seq.charAt( x ) ) {
                            ++mutations;
                            ++total_mutations;
                            final String mutation = prefix + ":" + parent_seq.charAt( x ) + ( x + 1 )
                                    + child_seq.charAt( x );
                            PropertiesList custom_data = node.getNodeData().getProperties();
                            if ( custom_data == null ) {
                                custom_data = new PropertiesList();
                            }
                            custom_data.addProperty( new Property( property_ref,
                                                                   mutation,
                                                                   "",
                                                                   XSD_STRING,
                                                                   AppliesTo.PARENT_BRANCH ) );
                            if ( branch_mutations.containsKey( mutation ) ) {
                                branch_mutations.put( mutation, ( branch_mutations.get( mutation ) + 1 ) );
                            }
                            else {
                                branch_mutations.put( mutation, 1 );
                            }
                        }
                    }
                    if ( mutations > 0 ) {
                        ++updated_nodes;
                    }
                    stats.addValue( mutations );
                }
                else {
                    ++branches_without_seqs;
                }
            }
        }
        System.out.println();
        System.out.println( "INFERRED MUTATIONS ON BRANCHES:" );
        branch_mutations.entrySet().forEach( entry -> {
            System.out.println( entry.getKey() + "\t" + entry.getValue() );
        } );
        System.out.println();
        if ( !ForesterUtil.isEmpty( reference_seqe_node_name )
                && !ForesterUtil.isEmpty( mut_vs_reference_seq_property_ref ) ) {
            externalMutations( prefix, reference_seqe_node_name, mut_vs_reference_seq_property_ref, p );
        }
        p.setRerootable( false );
        p.setRooted( true );
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( p, 0, outfile );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outfile + "]: " + e.getLocalizedMessage() );
        }
        System.out.println( "Wrote outtree to: " + outfile );
        System.out.println( "Total nodes (excluding root): " + total_nodes );
        System.out.println( "Branches with sequences     : " + branches_with_seqs );
        System.out.println( "Branches without sequences  : " + branches_without_seqs );
        System.out.println( "Updated nodes               : " + updated_nodes );
        System.out.println( "Sequence length             : " + seq_length );
        System.out.println( "Sum of mutations            : " + total_mutations );
        System.out.println( "Minimum mutations per branch: " + stats.getMin() );
        System.out.println( "Maximum mutations per branch: " + stats.getMax() );
        System.out.println( "Mean mutations per branch   : " + stats.arithmeticMean() );
        System.out.println( "Median mutations per branch : " + stats.median() );
    }

    private static void externalMutations( final String prefix,
                                           final String reference_seqe_node_name,
                                           final String mut_vs_reference_seq_property_ref,
                                           final Phylogeny p ) {
        final SortedMap<String, Integer> external_mutations = new TreeMap<>();
        PhylogenyNode ref_seq_node = null;
        try {
            ref_seq_node = p.getNode( reference_seqe_node_name );
        }
        catch ( final IllegalArgumentException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
        }
        String ref_seq = null;
        if ( ref_seq_node.isHasNodeData() && ( ref_seq_node.getNodeData().getSequence() != null )
                && !ForesterUtil.isEmpty( ref_seq_node.getNodeData().getSequence().getMolecularSequence() ) ) {
            ref_seq = ref_seq_node.getNodeData().getSequence().getMolecularSequence().trim().toUpperCase();
        }
        else {
            ForesterUtil.fatalError( PRG_NAME, "node sequence associated with node: " + reference_seqe_node_name );
        }
        final int ref_seq_length = ref_seq.length();
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                if ( node.isHasNodeData() && ( node.getNodeData().getSequence() != null )
                        && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() ) ) {
                    final String seq = node.getNodeData().getSequence().getMolecularSequence().trim().toUpperCase();
                    if ( seq.length() != ref_seq_length ) {
                        ForesterUtil.fatalError( PRG_NAME,
                                                 "node with unequal sequence length found: length found is "
                                                         + seq.length() + " expected " + ref_seq_length );
                    }
                    final int mutations = 0;
                    for( int x = 0; x < ref_seq_length; ++x ) {
                        if ( seq.charAt( x ) != ref_seq.charAt( x ) ) {
                            // ++mutations;
                            // ++total_mutations;
                            final String mutation = prefix + ":" + ref_seq.charAt( x ) + ( x + 1 ) + seq.charAt( x );
                            PropertiesList custom_data = node.getNodeData().getProperties();
                            if ( custom_data == null ) {
                                custom_data = new PropertiesList();
                            }
                            custom_data.addProperty( new Property( mut_vs_reference_seq_property_ref,
                                                                   mutation,
                                                                   "",
                                                                   XSD_STRING,
                                                                   AppliesTo.NODE ) );
                            if ( external_mutations.containsKey( mutation ) ) {
                                external_mutations.put( mutation, ( external_mutations.get( mutation ) + 1 ) );
                            }
                            else {
                                external_mutations.put( mutation, 1 );
                            }
                        }
                    }
                }
            }
        }
        System.out.println();
        System.out.println( "EXTERNAL MUTATIONS:" );
        external_mutations.entrySet().forEach( entry -> {
            System.out.println( entry.getKey() + "\t" + entry.getValue() );
        } );
        System.out.println();
    }

    private static void checkForOutputFileWriteability( final File outfile ) {
        final String error = ForesterUtil.isWritableFile( outfile );
        if ( !ForesterUtil.isEmpty( error ) ) {
            ForesterUtil.fatalError( PRG_NAME, error );
        }
    }
}
