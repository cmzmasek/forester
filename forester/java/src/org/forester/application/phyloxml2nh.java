
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
import org.forester.util.ForesterUtil;

public class phyloxml2nh {

    private static final String PRG_DATE    = "2024-05-28";
    private static final String PRG_VERSION = "1.0.1";
    private static final String PRG_NAME    = "phyloxml2nh";
    public static void main( final String args[] ) {
        if ( ( args.length != 3 && args.length != 2 ) ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree> <out-tree> [ref]\n" );
            System.exit( -1 );
        }

        final File intree = new File( args[ 0 ] );
        final File outtree = new File( args[ 1 ] );

        final String ref;
        if  ( args.length == 3 ) {
            ref = args[ 1 ];
        }
        else {
            ref = null;
        }

        final String e1 = ForesterUtil.isReadableFile( intree );
        if ( !ForesterUtil.isEmpty( e1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e1 );
        }
        final String e2 = ForesterUtil.isWritableFile( outtree );
        if ( !ForesterUtil.isEmpty( e2 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e2 );
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

        if ( ref != null ) {
            int counter = 0;
            for (final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
                final PropertiesList custom_data = node.getNodeData().getProperties();
                if (custom_data != null) {
                    final List<Property> props = custom_data.getProperties(ref);
                    if ((props.size() > 0) && (props.get(0) != null)) {
                        node.setName(props.get(0).getValue());
                        ++counter;
                    }
                }
            }
            System.out.println("Replaced " + counter + " names");
        }
        final PhylogenyWriter writer = new PhylogenyWriter();
        try {
            writer.toNewHampshire( p, true, PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE.AS_INTERNAL_NODE_NAMES,outtree );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.toString() );
        }
        System.out.println( "Wrote: " + outtree );
    }
}
