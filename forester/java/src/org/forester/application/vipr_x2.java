
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterUtil;

public class vipr_x2 {

    private static int           OPERATION      = 1;
    private static final String  XSD_STRING     = "xsd:string";
    private static final String  VIPR_HOST      = "vipr:Host";
    private static final String  VIPR_COUNTRY   = "vipr:Country";
    private static final String  VIPR_YEAR      = "vipr:Year";
    private final static String  PRG_NAME       = "vipr_x2";
    private static final String  PRG_DATE       = "2021-08-02";
    private static final String  PRG_VERSION    = "1.0.0";
    //  Name: Venezuelan_equine_encephalitis_virus|Tri550|VIPR_ALG4_257129192_1018_2283.1|FJ792614|NA|2007|Mosquito|Trinidad_and_Tobago|E2|E2_envelope_protein
    // 1. species
    // 2. Strain
    // 3. GB protein accession
    // 4. GB genome accession
    // 5. Uniprot acc
    // 6. Date
    // 7. Host
    // 8. Country
    // 9. gene symbol
    // 10. gene prod name
    private final static Pattern VIPR_PATTERN_1 = Pattern
            .compile( "(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)" );
    private final static Pattern YEAR_PATTERN   = Pattern.compile( "(\\d{4}).*" );
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        if ( args.length != 2 ) {
            System.out.println( "\nWrong number of arguments, expected: intree outtree\n" );
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
        if ( !infile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + infile + "] does not exist" );
        }
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
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
        final List<PhylogenyNode> ext_nodes = p.getExternalNodes();
        for( final PhylogenyNode ext_node : ext_nodes ) {
            final String name = ext_node.getName();
            if ( !ForesterUtil.isEmpty( name ) ) {
                if ( OPERATION == 1 ) {
                    final Matcher m = VIPR_PATTERN_1.matcher( name );
                    if ( m.matches() ) {
                        final String species = m.group( 1 );
                        final String strain = m.group( 2 );
                        final String gb_prot_acc = m.group( 3 );
                        final String gb_genome_acc = m.group( 4 );
                        final String uniprot_acc = m.group( 5 );
                        final String date = m.group( 6 );
                        final String host = m.group( 7 );
                        final String country = m.group( 8 );
                        String year = "NA";
                        final Matcher ym = YEAR_PATTERN.matcher( date );
                        if ( ym.matches() ) {
                            year = ym.group( 1 );
                        }
                        System.out.println( name );
                        System.out.println( "species      : " + species );
                        System.out.println( "strain       : " + strain );
                        System.out.println( "host         : " + host );
                        System.out.println( "country      : " + country );
                        System.out.println( "year         : " + year );
                        System.out.println();
                        final PropertiesList custom_data = new PropertiesList();
                        if ( !ForesterUtil.isEmpty( year ) ) {
                            custom_data.addProperty( new Property( VIPR_YEAR, year, "", XSD_STRING, AppliesTo.NODE ) );
                        }
                        custom_data
                                .addProperty( new Property( VIPR_COUNTRY, country, "", XSD_STRING, AppliesTo.NODE ) );
                        custom_data.addProperty( new Property( VIPR_HOST, host, "", XSD_STRING, AppliesTo.NODE ) );
                        ext_node.getNodeData().setProperties( custom_data );
                        Sequence seq = null;
                        if ( ext_node.isHasNodeData() && ext_node.getNodeData().isHasSequence() ) {
                            seq = ext_node.getNodeData().getSequence();
                        }
                        else {
                            seq = new Sequence();
                        }
                        seq.setAccession( new Accession( gb_prot_acc, Accession.Source.NCBI ) );
                        ext_node.getNodeData().addSequence( seq );
                        ext_node.setName( species + "|" + strain + "|" + host + "|" + country + "|" + year );
                    }
                    else {
                        System.out.println( "WARNING name \"" + name + "\" could not be matched" );
                        //System.exit( -1 );
                    }
                }
            }
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( p, 0, outfile );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outfile + "]: " + e.getMessage() );
        }
        System.out.println( "[" + PRG_NAME + "] wrote: [" + outfile + "]" );
        System.out.println( "[" + PRG_NAME + "] OK" );
        System.out.println();
    }
}
