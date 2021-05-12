
package org.forester.application;

import java.io.File;
import java.io.IOException;

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
        if ( args.length != 4 ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree> <out-tree> <suffix> <property reference>\n" );
            System.out.println( "Example: " + PRG_NAME + " tree_with_anc_seqs.xml outtree.xml S aptx:branch_event\n" );
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
        final String suffix = args[ 2 ];
        final String property_ref = args[ 3 ];
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
        int total_nodes = 0;
        int updated_nodes = 0;
        int total_mutations = 0;
        int branches_without_seqs = 0;
        int branches_with_seqs = 0;
        final BasicDescriptiveStatistics stats = new BasicDescriptiveStatistics();
        int seq_length = -1;
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
                                                 "node with unequal sequence lenght found: length found is "
                                                         + child_seq.length() + " expected " + seq_length );
                    }
                    if ( parent_seq.length() != seq_length ) {
                        ForesterUtil.fatalError( PRG_NAME,
                                                 "node with unequal sequence lenght found: length found is "
                                                         + parent_seq.length() + " expected " + seq_length );
                    }
                    int mutations = 0;
                    for( int x = 0; x < seq_length; ++x ) {
                        if ( parent_seq.charAt( x ) != child_seq.charAt( x ) ) {
                            ++mutations;
                            ++total_mutations;
                            final String mutation = suffix + ":" + parent_seq.charAt( x ) + ( x + 1 )
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

    private static void checkForOutputFileWriteability( final File outfile ) {
        final String error = ForesterUtil.isWritableFile( outfile );
        if ( !ForesterUtil.isEmpty( error ) ) {
            ForesterUtil.fatalError( PRG_NAME, error );
        }
    }
}
