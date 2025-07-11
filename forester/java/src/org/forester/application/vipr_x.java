
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
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
import org.forester.util.ViralUtils;

public class vipr_x {

    private static int           OPERATION         = 2;
    private static final String  XSD_STRING        = "xsd:string";
    private static final String  VIPR_HOST         = "vipr:Host";
    private static final String  VIPR_COUNTRY      = "vipr:Country";
    private static final String  VIPR_YEAR         = "vipr:Year";
    
    private static final String  VIPR_REGION       = "vipr:Region";
    private final static String  PRG_NAME          = "vipr_x";
    private static final String  PRG_DATE          = "2025-05-29";
    private static final String  PRG_VERSION       = "1.0.2";
    //  VP1|VP1_protein|KP322752|US/CA/14_6089|2014_08|Human|USA|NA|Enterovirus_D68
    // 1. gene symbol
    // 2. gene product name
    // 3. GB accession
    // 4. strain name
    // 5. Date
    // 6. Host
    // 7. Country
    // 8. Subtype
    // 9. Virus Type
    private final static Pattern VIPR_PATTERN_1    = Pattern
            .compile( "(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)" );
    // MT451688|SARS_CoV_2/human/AUS/VIC1004/2020|2020_03_31|Human|Australia|NA|Severe_acute_respiratory_syndrome_related_coronavirus
    // MW514307|SARS_CoV_2/human/RUS/Dubrovka/2020|2020_06_04|Human|Russia
    // 1. strain name
    // 2. GB accession
    // 3. Date
    // 4. Host
    // 5. Country
    private final static Pattern VIPR_PATTERN_2    = Pattern.compile( "(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|.+" );
    // Severe acute respiratory syndrome coronavirus 2 |OR469709.1|2023-08-13|Homo sapiens|USA|FL-CDC-LC1050185|
    private final static Pattern NCBI_PATTERN_1    = Pattern
            .compile( "(Severe\\s+acute\\s+respiratory\\s+syndrome\\s+coronavirus\\s+2\\s*)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|" );
    //                  Severe acute respiratory syndrome coronavirus 2                     |ac     |date   |host   |cou    |FL-CDC-LC1050185|
    private final static Pattern YEAR_PATTERN      = Pattern.compile( "(\\d{4}).*" );
    private final static boolean MODIFY_NODE_NAMES = true;
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
                        final String gene_symbol = m.group( 1 );
                        final String gene_product_name = m.group( 2 );
                        final String gb_accession = m.group( 3 );
                        final String strain_name = m.group( 4 );
                        final String date = m.group( 5 );
                        final String host = m.group( 6 );
                        final String country = m.group( 7 );
                        final String subtype = m.group( 8 );
                        final String virus_type = m.group( 9 );
                        String year = "";
                        final Matcher ym = YEAR_PATTERN.matcher( date );
                        if ( ym.matches() ) {
                            year = ym.group( 1 );
                        }
                        System.out.println( name );
                        System.out.println( "gene symbol      : " + gene_symbol );
                        System.out.println( "gene product name: " + gene_product_name );
                        System.out.println( "gb accession     : " + gb_accession );
                        System.out.println( "strain name      : " + strain_name );
                        System.out.println( "date             : " + date );
                        System.out.println( "year             : " + year );
                        System.out.println( "host             : " + host );
                        System.out.println( "country          : " + country );
                        System.out.println( "subtype          : " + subtype );
                        System.out.println( "virus type       : " + virus_type );
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
                        seq.setAccession( new Accession( gb_accession, Accession.Source.NCBI ) );
                        try {
                            seq.setSymbol( gene_symbol );
                        }
                        catch ( final PhyloXmlDataFormatException e ) {
                            // ignore
                        }
                        ext_node.getNodeData().addSequence( seq );
                        ViralUtils.addRegion( country, custom_data, VIPR_REGION );
                    }
                    else {
                        System.out.println( "ERROR: name \"" + name + "\" could not be matched" );
                        System.exit( -1 );
                    }
                }
                else if ( OPERATION == 2 ) {
                    final Matcher m = VIPR_PATTERN_2.matcher( name );
                    final Matcher m2 = NCBI_PATTERN_1.matcher( name );
                    if ( m.matches() ) {
                        final String gb_accession = m.group( 2 );
                        final String strain_name = m.group( 1 );
                        final String date = m.group( 3 );
                        String host = m.group( 4 );
                        final String country = m.group( 5 );
                        String year = "";
                        final Matcher ym = YEAR_PATTERN.matcher( date );
                        if ( ym.matches() ) {
                            year = ym.group( 1 );
                        }
                        if (host.equalsIgnoreCase("Homo sapiens")) {
                            host = "Human";
                        }
                        System.out.println( name );
                        System.out.println( "gb accession     : " + gb_accession );
                        System.out.println( "strain name      : " + strain_name );
                        System.out.println( "date             : " + date );
                        System.out.println( "year             : " + year );
                        System.out.println( "host             : " + host );
                        System.out.println( "country          : " + country );
                        System.out.println();
                        PropertiesList custom_data = ext_node.getNodeData().getProperties();
                        if ( custom_data == null ) {
                            custom_data = new PropertiesList();
                        }
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
                        seq.setAccession( new Accession( gb_accession, Accession.Source.NCBI ) );
                        ext_node.getNodeData().addSequence( seq );
                        ViralUtils.addRegion( country, custom_data, VIPR_REGION );
                    }
                    else if ( m2.matches() ) {
                        final String gb_accession = m2.group( 2 );
                        final String species = m2.group( 1 );
                        final String date = m2.group( 3 );
                        String host = m2.group( 4 );
                        final String country = m2.group( 5 );
                        final String genotype = m2.group( 6 );
                        String year = "";
                        final Matcher ym = YEAR_PATTERN.matcher( date );
                        if ( ym.matches() ) {
                            year = ym.group( 1 );
                        }
                        if (host.equalsIgnoreCase("Homo sapiens")) {
                            host = "Human";
                        }
                        System.out.println( name );
                        System.out.println( "gb accession     : " + gb_accession );
                        System.out.println( "strain name      : " + species );
                        System.out.println( "date             : " + date );
                        System.out.println( "year             : " + year );
                        System.out.println( "host             : " + host );
                        System.out.println( "country          : " + country );
                        System.out.println( "genotype         : " + genotype );
                        System.out.println();
                        PropertiesList custom_data = ext_node.getNodeData().getProperties();
                        if ( custom_data == null ) {
                            custom_data = new PropertiesList();
                        }
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
                        seq.setAccession( new Accession( gb_accession, Accession.Source.NCBI ) );
                        ext_node.getNodeData().addSequence( seq );
                        ViralUtils.addRegion( country, custom_data, VIPR_REGION );
                    }
                    else {
                        ForesterUtil.fatalError( PRG_NAME, "name \"" + name + "\" could not be matched" );
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
