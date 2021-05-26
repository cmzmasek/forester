
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
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

public class tree_annotator {

    private static final String PRG_DATE    = "2021-05-26";
    private static final String PRG_VERSION = "0.0.1";
    private static final String PRG_NAME    = "tree_annotator";
    private static final String XSD_STRING  = "xsd:string";
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        System.out.println();
        if ( ( args.length != 3 ) && ( args.length != 4 ) ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree> <table> <out-tree> [prefix]\n" );
            System.out.println( "Examples : " + PRG_NAME + " tree.nh map.txt outtree.xml\n" );
            System.out.println( "           " + PRG_NAME + " tree.nh map.txt outtree.xml bvbrc\n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final File intable = new File( args[ 1 ] );
        final File outtree = new File( args[ 2 ] );
        String prefix = null;
        if ( args.length == 4 ) {
            prefix = args[ 3 ];
        }
        final String e0 = ForesterUtil.isWritableFile( outtree );
        if ( !ForesterUtil.isEmpty( e0 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e0 );
        }
        final String e1 = ForesterUtil.isReadableFile( intree );
        if ( !ForesterUtil.isEmpty( e1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e1 );
        }
        final String e2 = ForesterUtil.isReadableFile( intable );
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
        ForesterUtil
                .programMessage( PRG_NAME,
                                 "Successfully read in tree with " + p.getNumberOfExternalNodes() + " external nodes" );
        BasicTable<String> table = null;
        try {
            table = BasicTableParser.parse( intable, '\t' );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read \"" + intable + "\" [" + e.getMessage() + "]" );
        }
        if ( table.getNumberOfColumns() < 3 ) {
            ForesterUtil.fatalError( PRG_NAME,
                                     "table has " + table.getNumberOfColumns() + " columns, expected at least 3" );
        }
        if ( table.getNumberOfRows() < 1 ) {
            ForesterUtil.fatalError( PRG_NAME, "table has no rows (i.e. is empty)" );
        }
        ForesterUtil.programMessage( PRG_NAME,
                                     "Successfully read in table with " + table.getNumberOfColumns() + " columns and "
                                             + table.getNumberOfRows() + " rows" );
        int properties_added = 0;
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                    final String name = node.getName();
                    for( int i = 0; i < table.getNumberOfRows(); ++i ) {
                        if ( table.getValue( 0, i ).equals( name ) ) {
                            String ref = table.getValueAsString( 1, i );
                            final String val = table.getValueAsString( 2, i );
                            if ( ForesterUtil.isEmpty( ref ) ) {
                                ForesterUtil.fatalError( PRG_NAME, "property reference empty for node " + node );
                            }
                            if ( ForesterUtil.isEmpty( val ) ) {
                                ForesterUtil.fatalError( PRG_NAME, "property value empty for node " + node );
                            }
                            if ( !ForesterUtil.isEmpty( prefix ) ) {
                                ref = prefix + ":" + ref;
                            }
                            else {
                                if ( ref.indexOf( ':' ) < 0 ) {
                                    ForesterUtil
                                            .fatalError( PRG_NAME,
                                                         "property reference " + "\"" + ref + "\" lacks required :" );
                                }
                            }
                            PropertiesList custom_data = node.getNodeData().getProperties();
                            if ( custom_data == null ) {
                                custom_data = new PropertiesList();
                                node.getNodeData().setProperties( custom_data );
                            }
                            custom_data.addProperty( new Property( ref, val, "", XSD_STRING, AppliesTo.NODE ) );
                            ++properties_added;
                        }
                    }
                }
            }
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( p, 0, outtree );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outtree + "]: " + e.getLocalizedMessage() );
        }
        System.out.println();
        System.out.println( "Sum of properties added: " + properties_added );
        System.out.println( "Wrote outtree to       : " + outtree );
        System.out.println();
    }
}
