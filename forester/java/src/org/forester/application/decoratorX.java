// java -Xmx2048m -cp
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.application.decoratorX
// RRMa_ALL_plus_RRMa_ee3_50_hmmalign_05_40_fme_with_seqs_2.phylo.xml
// nature12311-s3_cz_4.txt x1 x2

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

    private static final int SEQ_NAME_COLUMN = 1;
    private static final int SPECIES_COLUMN  = 2;
    private static final int SEQ_COLUMN      = 3;
    private static final int TARGET_COLUMN   = 4;

    public static void main( final String args[] ) {
        File intree = null;
        File outtree1 = null;
        File outtree2 = null;
        File intable = null;
        try {
            CommandLineArguments cla = null;
            cla = new CommandLineArguments( args );
            intree = cla.getFile( 0 );
            intable = cla.getFile( 1 );
            outtree1 = cla.getFile( 2 );
            outtree2 = cla.getFile( 3 );
            if ( outtree1.exists() ) {
                System.out.println( outtree1 + " already exists" );
                System.exit( -1 );
            }
            if ( outtree2.exists() ) {
                System.out.println( outtree2 + " already exists" );
                System.exit( -1 );
            }
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
            final Phylogeny phy = factory.create( intree, xml_parser )[ 0 ];
            final BasicTable<String> t = BasicTableParser.parse( intable, '\t' );
            final PhylogenyNodeIterator it = phy.iteratorExternalForward();
            int i = 0;
            while ( it.hasNext() ) {
                final PhylogenyNode node = it.next();
                processNode( node, t );
                i++;
            }
            final PhylogenyWriter writer1 = new PhylogenyWriter();
            writer1.toPhyloXML( outtree1, phy, 0 );
            final PhylogenyNodeIterator it2 = phy.iteratorExternalForward();
            while ( it2.hasNext() ) {
                final PhylogenyNode node = it2.next();
                processNode2( node, phy );
            }
            final PhylogenyWriter writer2 = new PhylogenyWriter();
            writer2.toPhyloXML( outtree2, phy, 0 );
        }
        catch ( final Exception e ) {
            System.out.println( e.getLocalizedMessage() );
            System.exit( -1 );
        }
    }

    private static void processNode( final PhylogenyNode node, final BasicTable<String> t ) throws Exception {
        final String node_seq = node.getNodeData().getSequence().getMolecularSequence().toUpperCase();
        boolean found = false;
        String found_row = "";
        String found_protein_name = "";
        String found_species = "";
        for( int row = 0; row < t.getNumberOfRows(); ++row ) {
            final String table_seq = t.getValueAsString( SEQ_COLUMN, row ).toUpperCase();
            if ( table_seq.contains( node_seq ) ) {
                if ( found ) {
                    if ( !found_protein_name.equals( t.getValueAsString( SEQ_NAME_COLUMN, row ) )
                            || !found_species.equals( t.getValueAsString( SPECIES_COLUMN, row ) ) ) {
                        throw new Exception( "Sequence from node " + node + " is not unique: " + node_seq + "\n"
                                + "Already found in row " + found_row );
                    }
                }
                else {
                    found = true;
                    found_row = t.getRowAsString( row, ", " );
                    found_protein_name = t.getValueAsString( SEQ_NAME_COLUMN, row );
                    found_species = t.getValueAsString( SPECIES_COLUMN, row );
                }
                final Annotation annotation = new Annotation( "target", t.getValueAsString( TARGET_COLUMN, row ) );
                node.getNodeData().getSequence().addAnnotation( annotation );
                System.out.println( node + "->" + annotation );
            }
        }
    }

    private static void processNode2( final PhylogenyNode node, final Phylogeny t ) {
        if ( ( node.getNodeData().getSequence().getAnnotations() == null )
                || node.getNodeData().getSequence().getAnnotations().isEmpty() ) {
            t.deleteSubtree( node, true );
        }
    }
}
