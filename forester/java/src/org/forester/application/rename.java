
package org.forester.application;

import java.io.File;
import java.io.IOException;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

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
            final String node_name = node.getName();
            if ( !ForesterUtil.isEmpty( node_name ) ) {
                final int i = node_name.lastIndexOf( '_' );
                if ( i > 0 ) {
                    node.setName( node_name.substring( i + 1 )  );
                }
            }
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toNewHampshire( p, true, outfile );
        }
        catch ( final IOException e ) {
            System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
    }
}
