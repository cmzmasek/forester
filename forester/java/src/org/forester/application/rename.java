
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class rename {

    public static void main( final String args[] ) {
        final File infile = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( infile, true );
            p = factory.create( infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + infile + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            /*final String node_name = node.getName();
            if ( node.isExternal() && !ForesterUtil.isEmpty( node_name ) ) {
                final int i = node_name.lastIndexOf( '_' );
                if ( i > 0 ) {
                    node.setName( node_name.substring( i + 1 ) );
                }
            }*/
            if ( node.isExternal() ) {
                final PropertiesList custom_data = node.getNodeData().getProperties();
                final List<Property> pl = custom_data.getProperties( "BVBRC:strain" );
                node.setName( pl.get( 0 ).getValue() );
            }
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
           // w.toNewHampshire( p, true, outfile );
            w.toPhyloXML( p, 0, outfile );
        }
        catch ( final IOException e ) {
            System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
    }
}
