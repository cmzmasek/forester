
package org.forester.application;

import java.io.File;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class tree_a {

    private static final String PRG_DATE    = "2021-12-06";
    private static final String PRG_VERSION = "0.0.2";
    private static final String PRG_NAME    = "tree_a";
    public static void main( final String args[] ) {
        if ( ( args.length != 2 ) ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree> <dir>\n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final File dir = new File( args[ 1 ] );
        final String e1 = ForesterUtil.isReadableFile( intree );
        if ( !ForesterUtil.isEmpty( e1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e1 );
        }
        final String query_seq = dir.getName().substring( 6 ).replace( '~', '|' );
        
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intree, true );
            p = factory.create( intree, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "Could not read \"" + intree + "\" [" + e.getMessage() + "]" );
        }
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() && !node.isRoot() ) {
                if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                    final String name = node.getName();
                    if ( name.indexOf( '|' ) > 0 ) {
                        final String[] x = name.split( "\\|" );
                        final String subspecies = x[ x.length - 1 ];
                        System.out.print( query_seq );
                        System.out.print( "\t" );
                        
                        System.out.print( subspecies );
                        System.out.print( "\t" );
                        
                        final PhylogenyNode parent = node.getParent();
                        final List<String> ds = parent.getAllExternalDescendantsNames();
                        System.out.print( ds.size() - 1 );
                        System.out.print( "\t" );
                        for( final String desc : ds ) {
                            if ( desc.indexOf( '|' ) < 0 ) {
                                System.out.print( desc );
                                System.out.print( "\t" );
                            }
                        }
                        System.out.println();
                    }
                }
            }
        }
    }
}
