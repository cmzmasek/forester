// java -Xmx2048m -cp
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.application.decoratorX
// RRMa_ALL_plus_RRMa_ee3_50_hmmalign_05_40_fme_with_seqs_2.phylo.xml
// nature12311-s3_cz_4.txt x

package org.forester.application;

import java.io.File;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;

public class decoratorX {

    private static final int SEQ_COLUMN    = 3;
    private static final int TARGET_COLUMN = 4;

    public static void main( final String args[] ) {
        File intree = null;
        File outtree = null;
        File intable = null;
        try {
            CommandLineArguments cla = null;
            cla = new CommandLineArguments( args );
            intree = cla.getFile( 0 );
            intable = cla.getFile( 1 );
            outtree = cla.getFile( 2 );
            if ( outtree.exists() ) {
                System.out.println( outtree + " already exists" );
                System.exit( -1 );
            }
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = new PhyloXmlParser();
            final Phylogeny phy = factory.create( intree, xml_parser )[ 0 ];
            final BasicTable<String> t = BasicTableParser.parse( intable, '\t' );
            //  System.out.println( t.toString() );
            final PhylogenyNodeIterator it = phy.iteratorExternalForward();
            int i = 0;
            while ( it.hasNext() ) {
                final PhylogenyNode node = it.next();
                processNode( node, t );
                i++;
            }
            final PhylogenyNodeIterator it2 = phy.iteratorExternalForward();
            while ( it2.hasNext() ) {
                final PhylogenyNode node = it2.next();
                processNode2( node, phy );
            }
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( outtree, phy, 0 );
        }
        catch ( final Exception e ) {
            System.out.println( e.getLocalizedMessage() );
            e.printStackTrace();
            System.exit( -1 );
        }
    }

    private static void processNode( final PhylogenyNode node, final BasicTable<String> t ) throws Exception {
        final String node_seq = node.getNodeData().getSequence().getMolecularSequence().toUpperCase();
        boolean found = false;
        for( int col = 0; col < t.getNumberOfRows(); ++col ) {
            final String table_seq = t.getValueAsString( SEQ_COLUMN, col ).toUpperCase();
            if ( table_seq.contains( node_seq ) ) {
                if ( found ) {
                    // throw new Exception( "Sequence from node " + node + " is not unique: " + node_seq );
                }
                found = true;
                final Annotation annotation = new Annotation( "target", t.getValueAsString( TARGET_COLUMN, col ) );
                node.getNodeData().getSequence().addAnnotation( annotation );
                System.out.println( node + "->" + annotation );
            }
        }
        // if ( !found ) {
        //     throw new Exception( "Sequence from node " + node + " not found: " + node_seq );
        // }
    }

    private static void processNode2( final PhylogenyNode node, final Phylogeny t ) {
        if ( ( node.getNodeData().getSequence().getAnnotations() == null )
                || node.getNodeData().getSequence().getAnnotations().isEmpty() ) {
            t.deleteSubtree( node, true );
        }
    }
}
