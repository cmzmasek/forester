
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

    private static final String PRG_DATE    = "2021-05-27";
    private static final String PRG_VERSION = "1.0.1";
    private static final String PRG_NAME    = "mutation_annotator";
    private static final String XSD_STRING  = "xsd:string";
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        if ( ( args.length != 5 ) && ( args.length != 6 ) && ( args.length != 7 ) ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME
                    + " <in-tree> <out-tree with seqs removed> <out-tree with seqs retained> <prefix> <property reference> [name of node with reference seq] [property reference for mutations vs reference seq]\n" );
            System.out.println( "Examples : " + PRG_NAME
                    + " tree_with_anc_seqs.xml outtree.xml outtree_keep_seqs.xml S aptx:branch_event\n" );
            System.out.println( "           " + PRG_NAME
                    + " tree_with_anc_seqs.xml outtree.xml outtree_keep_seqs.xml S aptx:branch_event 2019_nCoV_WHU01 vipr:Mutation\n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final File outtree_rem_seqs = new File( args[ 1 ] );
        final File outtree_keep_seqs = new File( args[ 2 ] );
        final String prefix = args[ 3 ];
        final String property_ref = args[ 4 ];
        String reference_seqe_node_name = null;
        String mut_vs_reference_seq_property_ref = null;
        if ( args.length == 6 ) {
            mut_vs_reference_seq_property_ref = args[ 5 ];
        }
        if ( args.length == 7 ) {
            reference_seqe_node_name = args[ 5 ];
            mut_vs_reference_seq_property_ref = args[ 6 ];
        }
        checkForOutputFileWriteability( outtree_rem_seqs );
        checkForOutputFileWriteability( outtree_keep_seqs );
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intree, true );
            p = factory.create( intree, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "Could not read \"" + intree + "\" [" + e.getMessage() + "]" );
        }
        System.out.println();
        printNodesWithoutSequences( p );
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
                                node.getNodeData().setProperties( custom_data );
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
        System.out.println();
        System.out.println();
        System.out.println();
        System.out.println( "INFERRED MUTATIONS ON BRANCHES, CONVERGENT/PARALLEL ONLY:" );
        branch_mutations.entrySet().forEach( entry -> {
            if ( entry.getValue() > 1 ) {
                System.out.println( entry.getKey() + "\t" + entry.getValue() );
            }
        } );
        System.out.println();
        System.out.println();
        branch_mutations.entrySet().forEach( entry -> {
            if ( entry.getValue() > 1 ) {
                System.out.print( "'" + entry.getKey() + "', " );
            }
        } );
        
        System.out.println();
        System.out.println();
        System.out.println();
        System.out.println();
        System.out.println( "INFERRED MUTATIONS ON BRANCHES, NON-PARALLEL ONLY:" );
        branch_mutations.entrySet().forEach( entry -> {
            if ( entry.getValue() == 1 ) {
                System.out.println( entry.getKey() + "\t" + entry.getValue() );
            }
        } );
        System.out.println();
        System.out.println();
        branch_mutations.entrySet().forEach( entry -> {
            if ( entry.getValue() == 1 ) {
                System.out.print( "'" + entry.getKey() + "', " );
            }
        } );
        
        System.out.println();
        System.out.println();
        if ( !ForesterUtil.isEmpty( reference_seqe_node_name )
                && !ForesterUtil.isEmpty( mut_vs_reference_seq_property_ref ) ) {
            externalMutations( prefix, reference_seqe_node_name, mut_vs_reference_seq_property_ref, p );
        }
        else if ( ForesterUtil.isEmpty( reference_seqe_node_name )
                && !ForesterUtil.isEmpty( mut_vs_reference_seq_property_ref ) ) {
            externalMutations( prefix, null, mut_vs_reference_seq_property_ref, p );
        }
        p.setRerootable( false );
        p.setRooted( true );
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( p, 0, outtree_keep_seqs );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME,
                                     "failed to write to [" + outtree_keep_seqs + "]: " + e.getLocalizedMessage() );
        }
        removeMolecularSequences( p );
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( p, 0, outtree_rem_seqs );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME,
                                     "failed to write to [" + outtree_rem_seqs + "]: " + e.getLocalizedMessage() );
        }
        System.out.println( "Wrote outtree to            : " + outtree_rem_seqs );
        System.out.println( "Wrote outtree to            : " + outtree_keep_seqs );
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

    private static void removeMolecularSequences( final Phylogeny p ) {
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                if ( node.getNodeData().getSequence() != null ) {
                    node.getNodeData().getSequence().setMolecularSequence( null );
                    node.getNodeData().getSequence().setMolecularSequenceAligned( false );
                    node.getNodeData().getSequence().setGeneName( null );
                    node.getNodeData().getSequence().setName( null );
                }
            }
            else {
                node.getNodeData().setSequence( null );
            }
        }
    }

    private static void printNodesWithoutSequences( final Phylogeny p ) {
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( !node.isHasMolecularSequence() ) {
                if ( node.isInternal() ) {
                    System.out.println( "no sequence for internal node: " + node );
                    if ( node.getChildNode1().isExternal() ) {
                        System.out.println( "    external child 1: " + node.getChildNode1() );
                    }
                    else {
                        System.out.println( "    intrenal child 1: " + node.getChildNode1() );
                    }
                    if ( node.getChildNode2().isExternal() ) {
                        System.out.println( "    external child 2: " + node.getChildNode2() );
                    }
                    else {
                        System.out.println( "    internal child 2: " + node.getChildNode2() );
                    }
                }
                else {
                    System.out.println( "no sequence for external node: " + node );
                }
                System.out.println();
            }
        }
    }

    private static void externalMutations( final String prefix,
                                           final String reference_seqe_node_name,
                                           final String mut_vs_reference_seq_property_ref,
                                           final Phylogeny p ) {
        final SortedMap<String, Integer> external_mutations = new TreeMap<>();
        String ref_seq = null;
        if ( reference_seqe_node_name != null ) {
            PhylogenyNode ref_seq_node = null;
            try {
                ref_seq_node = p.getNode( reference_seqe_node_name );
            }
            catch ( final IllegalArgumentException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
            }
            if ( ref_seq_node.isHasMolecularSequence() ) {
                ref_seq = ref_seq_node.getNodeData().getSequence().getMolecularSequence().trim().toUpperCase();
            }
            else {
                ForesterUtil.fatalError( PRG_NAME, "node sequence associated with node: " + reference_seqe_node_name );
            }
        }
        else {
            ref_seq = getReferenceSequence();
        }
        final int ref_seq_length = ref_seq.length();
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                if ( node.isHasMolecularSequence() ) {
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
                                node.getNodeData().setProperties( custom_data );
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

    private static String getReferenceSequence() {
        return "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFD"
                + "NPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVY"
                + "SSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQT"
                + "LLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRV"
                + "QPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSF"
                + "VIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPC"
                + "NGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFL"
                + "PFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGS"
                + "NVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTI"
                + "SVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGF"
                + "NFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAG"
                + "TITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALN"
                + "TLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRV"
                + "DFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNT"
                + "FVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDL"
                + "QELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT";
    }
}
