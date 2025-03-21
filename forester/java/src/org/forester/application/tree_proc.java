
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class tree_proc {

    private static final String PRG_DATE     = "2021-10-08";
    private static final String PRG_VERSION  = "0.0.1";
    private static final String PRG_NAME     = "tree_proc";
    private static final String XSD_STRING   = "xsd:string";
    private static final String VIPR_HOST    = "vipr:Hosts";
    private static final String VIPR_SPECIES = "vipr:Species";
    private static final String VIPR_YEAR    = "vipr:Year";
    private static final String VIPR_EC      = "vipr:EC";
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
        ForesterUtil
                .programMessage( PRG_NAME,
                                 "Successfully read in tree with " + p.getNumberOfExternalNodes() + " external nodes" );
        final int properties_added = 0;
        final Pattern species_p = Pattern.compile( "\\[(\\S+\\s+\\S+)" );
        final Pattern tax_p = Pattern.compile( "\\[(\\S+\\s+[^\\\\|]+)" );
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                node.getNodeData().setSequence(null);

                /*if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                    final String name = node.getName();
                    String species = "";
                    final Matcher m = species_p.matcher( name );
                    if ( m.find() ) {
                        species = m.group( 1 );
                    }
                    else {
                        System.out.println( name );
                    }
                    final Matcher m2 = tax_p.matcher( name );
                    String tax = "";
                    if ( m2.find() ) {
                        tax = m2.group( 1 );
                        System.out.println( tax );
                        final Taxonomy t = new Taxonomy();
                        t.setScientificName( tax );
                        node.getNodeData().setTaxonomy( t );
                    }
                    final Sequence s = new Sequence();
                    s.setName( "DNA polymerase III epsilon subunit" );
                    final Annotation a = new Annotation( "EC", "2.7.7.7" );
                    s.addAnnotation( a );
                    node.getNodeData().addSequence( s );
                    ///////////
                    String host = "Human";
                    final String ec = "2.7.7.7";
                    if ( tax.indexOf( "felis" ) > -1 ) {
                        host = "Human, Cat, Flea";
                    }
                    if ( tax.indexOf( "endosymbiont of Proechinophthirus fluctus" ) > -1 ) {
                        host = "Lice";
                    }
                    if ( tax.indexOf( "endosymbiont of Cimex lectularius" ) > -1 ) {
                        host = "Bed bug";
                    }
                    if ( tax.indexOf( "endosymbiont of Ixodes pacificus" ) > -1 ) {
                        host = "Tick";
                    }
                    node.setName( tax + "|DNA polymerase III epsilon subunit");
                    PropertiesList custom_data = node.getNodeData().getProperties();
                    if ( custom_data == null ) {
                        custom_data = new PropertiesList();
                        node.getNodeData().setProperties( custom_data );
                    }
                    if ( !ForesterUtil.isEmpty( host ) ) {
                        custom_data.addProperty( new Property( VIPR_HOST, host, "", XSD_STRING, AppliesTo.NODE ) );
                    }
                    if ( !ForesterUtil.isEmpty( ec ) ) {
                        custom_data.addProperty( new Property( VIPR_EC, ec, "", XSD_STRING, AppliesTo.NODE ) );
                    }
                    if ( !ForesterUtil.isEmpty( species ) ) {
                        custom_data
                                .addProperty( new Property( VIPR_SPECIES, species, "", XSD_STRING, AppliesTo.NODE ) );
                    }
                    ///////////
                }
                */

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
