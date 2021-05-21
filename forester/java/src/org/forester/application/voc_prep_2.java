
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.SortedMap;

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

public class voc_prep_2 {

    private final static String PRG_NAME      = "voc_prep_2";
    private static final String PREFIX        = "S:";
    private static final String XSD_STRING    = "xsd:string";
    private static final String VIPR_MUTATION = "vipr:Mutation";
    public static void main( final String args[] ) {
        final File intree = new File( args[ 0 ] );
        final File genome_to_mut_infile = new File( args[ 1 ] );
        final File outfile = new File( args[ 2 ] );
        if ( !intree.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + intree + "] does not exist" );
        }
        if ( !genome_to_mut_infile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + genome_to_mut_infile + "] does not exist" );
        }
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intree, true );
            p = factory.create( intree, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + intree + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        BasicTable<String> genome_to_mut = null;
        try {
            genome_to_mut = BasicTableParser.parse( genome_to_mut_infile, '\t', false, false );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME,
                                     "failed to read [" + genome_to_mut_infile + "] [" + e.getMessage() + "]" );
        }
        final SortedMap<String, String> genome_to_mut_map = genome_to_mut.getColumnsAsMap( 0, 1 );
        
        
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isHasNodeData() && node.getNodeData().isHasSequence()
                    && !node.getNodeData().getSequence().isEmpty() ) {
                if ( ( node.getNodeData().getSequence().getAccession() != null )
                        && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() ) ) {
                    final String acc = node.getNodeData().getSequence().getAccession().getValue();
                    if ( genome_to_mut_map.containsKey( acc ) ) {
                        final String mut = genome_to_mut_map.get( acc );
                        final String muts[] = mut.split( "," );
                        for( final String m : muts ) {
                            PropertiesList custom_data = node.getNodeData().getProperties();
                            if ( custom_data == null ) {
                                custom_data = new PropertiesList();
                            }
                            final String mut_str = PREFIX + m;
                            custom_data.addProperty( new Property( VIPR_MUTATION,
                                                                   mut_str,
                                                                   "",
                                                                   XSD_STRING,
                                                                   AppliesTo.NODE ) );
                            node.getNodeData().setProperties( custom_data );
                            System.out.println( acc + "->" + mut_str );
                        }
                    }
                }
            }
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( p, 0, outfile );
        }
        catch ( final IOException e ) {
            System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
    }
}