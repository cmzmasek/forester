
package org.forester.application;

import java.io.File;
import java.io.IOException;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class cladinator_tree_prepare {

    private static final String PRG_DATE    = "2022-05-22";
    private static final String PRG_VERSION = "0.0.1";
    private static final String PRG_NAME    = "cladinator_tree_prepare";
    private static final String REF         = "subspecies:clade";
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        System.out.println();
        if ( ( args.length != 2 ) ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree> <out-tree> \n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final File outtree = new File( args[ 1 ] );
        final String e0 = ForesterUtil.isWritableFile( outtree );
        if ( !ForesterUtil.isEmpty( e0 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e0 );
        }
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
        if ( !p.isRooted() ) {
            ForesterUtil.fatalError( PRG_NAME, "\"" + intree + "\" is not rooted" );
        }
        p.setRerootable( false );
        ForesterUtil
                .programMessage( PRG_NAME,
                                 "Successfully read in tree with " + p.getNumberOfExternalNodes() + " external nodes" );
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                if ( node.isHasNodeData() && ( node.getNodeData().getProperties() != null )
                        && ( node.getNodeData().getProperties().size() > 0 )
                        && ( PhylogenyMethods.getNodePropertyValues( node, REF ).size() == 1 ) ) {
                    node.setName( PhylogenyMethods.getNodePropertyValues( node, REF ).get( 0 ) );
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "No annotation found for node " + node.getName() );
                }
            }
            else {
                node.setName( "" );
            }
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toNewHampshire( p, true, outtree );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outtree + "]: " + e.getLocalizedMessage() );
        }
        System.out.println();
        System.out.println( "Wrote outtree to: " + outtree );
        System.out.println();
    }
}
