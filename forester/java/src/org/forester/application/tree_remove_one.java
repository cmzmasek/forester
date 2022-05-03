
package org.forester.application;

import java.io.File;
import java.io.IOException;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterUtil;

public class tree_remove_one {

    private static final String PRG_NAME = "tree_remove_one";
    public static void main( final String args[] ) {
        if ( ( args.length != 1 ) ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree>\n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final String e1 = ForesterUtil.isReadableFile( intree );
        if ( !ForesterUtil.isEmpty( e1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e1 );
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intree, true );
            p = factory.create( intree, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "Could not read \"" + intree + "\" [" + e.getMessage() + "]" );
        }
        final String[] ext_node_names = p.getAllExternalNodeNames();
        for( final String node_name : ext_node_names ) {
            final Phylogeny p_copy = p.copy();
            p_copy.deleteSubtree( p_copy.getNode( node_name ), false );
            final PhylogenyWriter w = PhylogenyWriter.createPhylogenyWriter();
            final File out_file = new File( intree.getAbsolutePath()
                    .substring( 0, intree.getAbsolutePath().indexOf( '.' ) ) + "_" + node_name + ".nwk" );
            try {
                w.toNewHampshire( p_copy, true, out_file );
            }
            catch ( final IOException e ) {
                e.printStackTrace();
            }
        }
    }
}
