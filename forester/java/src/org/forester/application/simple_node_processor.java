
package org.forester.application;

import java.io.File;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.CommandLineArguments;

public class simple_node_processor {

    private final static String BASE = "b_";

    public static void main( final String args[] ) {
        File in = null;
        File out = null;
        try {
            CommandLineArguments cla = null;
            cla = new CommandLineArguments( args );
            in = cla.getFile( 0 );
            out = cla.getFile( 1 );
            if ( out.exists() ) {
                System.out.println( out + " already exists" );
                System.exit( -1 );
            }
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = new PhyloXmlParser();
            final Phylogeny[] phylogenies_0 = factory.create( in, xml_parser );
            final Phylogeny phylogeny_0 = phylogenies_0[ 0 ];
            // final PhylogenyNodeIterator it = phylogeny_0.iteratorPostorder();
            final PhylogenyNodeIterator it = phylogeny_0.iteratorExternalForward();
            int i = 0;
            while ( it.hasNext() ) {
                final PhylogenyNode node = it.next();
                processNode( node, i, out.toString() );
                i++;
            }
            //   final PhylogenyWriter writer = new PhylogenyWriter();
            //   writer.toPhyloXML( out, phylogeny_0, 0 );
        }
        catch ( final Exception e ) {
            System.out.println( e.getLocalizedMessage() );
            e.printStackTrace();
            System.exit( -1 );
        }
    }

    //    private static void processNode( final PhylogenyNode node, final int i ) {
    //        node.setDistanceToParent( PhylogenyNode.DISTANCE_DEFAULT );
    //        if ( !node.isExternal() ) {
    //            if ( ( node.getName() == null ) || node.getName().isEmpty() ) {
    //                node.setName( BASE + i );
    //            }
    //        }
    //    }
    private static void processNode( final PhylogenyNode node, final int i, final String label ) {
        //if ( node.isExternal() ) {
        //    final String c = "" + node.getNodeData().getBinaryCharacters().getPresentCount();
        //    final String s = node.getNodeData().getTaxonomy().getScientificName();
        //    System.out.println( s + "\t" + c );
        //}
        //        if ( !node.isExternal() ) {
        //            if ( !node.getNodeData().isHasTaxonomy() ) {
        //                if ( !ForesterUtil.isEmpty( node.getName() ) ) {
        //                    if ( ( node.getName().indexOf( "_" ) < 0 ) && ( node.getName().indexOf( "&" ) < 0 )
        //                            && ( node.getName().indexOf( " " ) < 0 ) ) {
        //                        Taxonomy t = new Taxonomy();
        //                        t.setScientificName( node.getName() );
        //                        node.getNodeData().addTaxonomy( t );
        //                        node.setName( "" );
        //                    }
        //                }
        //            }
        //        }
        if ( node.isExternal() ) {
            if ( node.getNodeData().isHasTaxonomy() ) {
                //  final Taxonomy t = node.getNodeData().getTaxonomy();
                //  t.setIdentifier( null );
                //if ( !ForesterUtil.isEmpty( t.getTaxonomyCode() ) && t.getTaxonomyCode().length() == 5 ) {
                //    if ( node.getName().equalsIgnoreCase( t.getTaxonomyCode() ) ) {
                //        node.setName( "" );
                //    }
                //}
                // node.setName( "" );
                final Taxonomy t = node.getNodeData().getTaxonomy();
                System.out.println( t.getTaxonomyCode() + "\t" + t.getScientificName() + "\t" + t.getCommonName()
                        + "\t" + label );
            }
            else {
                //System.out.println( "node " + node + " has not tax" );
            }
        }
    }
}