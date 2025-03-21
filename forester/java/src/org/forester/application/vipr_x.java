
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

public class vipr_x {

    private static int           OPERATION         = 2;
    private static final String  XSD_STRING        = "xsd:string";
    private static final String  VIPR_HOST         = "vipr:Host";
    private static final String  VIPR_COUNTRY      = "vipr:Country";
    private static final String  VIPR_YEAR         = "vipr:Year";
    
    private static final String  VIPR_REGION       = "vipr:Region";
    private final static String  PRG_NAME          = "vipr_x";
    private static final String  PRG_DATE          = "2023-09-20";
    private static final String  PRG_VERSION       = "1.0.1";
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
                        addRegion( country, custom_data );
                        addSpecialData1( gb_accession, custom_data, ext_node );
                        addSpecialData12( gb_accession, custom_data, ext_node );
                        addSpecialData2( gb_accession, custom_data, ext_node );
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
                        addRegion( country, custom_data );
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
                        addRegion( country, custom_data );
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

    private static void addRegion( final String country, final PropertiesList custom_data ) {
        String region = "";
        final String c = country.toLowerCase();
        if ( c.equals( "canada" ) || c.equals( "usa" ) || c.equals( "mexico" ) ) {
            region = "North America";
        }
        else if ( c.equals( "peru" ) || c.equals( "ecuador" ) || c.equals( "colombia" ) || c.equals( "chile" )
                || c.equals( "brazil" ) || c.equals( "argentina" ) || c.equals( "guatemala" ) || c.equals( "uruguay" )
                || c.equals( "venezuela" ) ) {
            region = "South America";
        }
        else if ( c.equals( "denmark" ) || c.equals( "finland" ) || c.equals( "france" ) || c.equals( "germany" )
                || c.equals( "netherlands" ) || c.equals( "norway" ) || c.equals( "united_kingdom" )
                || c.equals( "switzerland" ) || c.equals( "austria" ) || c.equals( "estonia" ) || c.equals( "sweden" )
                || c.equals( "belgium" ) ) {
            region = "Western Europe";
        }
        else if ( c.equals( "serbia" ) || c.equals( "greece" ) || c.equals( "malta" ) || c.equals( "italy" )
                || c.equals( "spain" ) || c.equals( "portugal" ) ) {
            region = "Southern Europe";
        }
        else if ( c.equals( "poland" ) ) {
            region = "Central Europe";
        }
        else if ( c.equals( "russia" ) || c.equals( "belarus" ) ) {
            region = "Eastern Europe";
        }
        else if ( c.equals( "japan" ) || c.equals( "taiwan" ) || c.equals( "hong_kong" ) || c.equals( "south_korea" )
                || c.equals( "china" ) ) {
            region = "East Asia";
        }
        else if ( c.equals( "kazakhstan" ) || c.equals( "uzbekistan" ) || c.equals( "armenia" ) ) {
            region = "Central Asia";
        }
        else if ( c.equals( "kuwait" ) || c.equals( "jordan" ) || c.equals( "bahrain" ) || c.equals( "iraq" )
                || c.equals( "saudi arabia" ) || c.equals( "saudi_arabia" ) || c.equals( "turkey" ) || c.equals( "egypt" ) || c.equals( "israel" )
                || c.equals( "west_bank" ) || c.equals( "iran" ) ) {
            region = "West Asia";
        }
        else if ( c.equals( "india" )
                || ( ( c.equals( "pakistan" ) | c.equals( "bangladesh" ) ) || c.equals( "sri_lanka" ) ) ) {
            region = "South Asia";
        }
        else if ( c.equals( "laos" ) || c.equals( "cambodia" ) || c.equals( "thailand" ) || c.equals( "malaysia" )
                || c.equals( "philippines" ) || ( c.equals( "viet_nam" ) | c.equals( "viet nam" ) )
                || c.equals( "myanmar" ) || c.equals( "timor_leste" ) ) {
            region = "Southeast Asia";
        }
        else if ( c.equals( "morocco" ) || c.equals( "gambia" ) || c.equals( "kenya" ) || c.equals( "senegal" )
                || c.equals( "south_africa" ) || c.equals( "tanzania" ) || c.equals( "ghana" ) || c.equals( "benin" )
                || c.equals( "tunisia" ) || c.equals( "nigeria" ) || c.equals( "libya" ) || c.equals( "djibouti" )
                || c.equals( "sierra_leone" ) || c.equals( "guinea" ) || c.equals( "botswana" )
                || c.equals( "ethiopia" ) || c.equals( "malawi" ) || c.equals( "mali" ) ) {
            region = "Africa";
        }
        else if ( c.equals( "australia" ) || c.equals( "new_zealand" ) ) {
            region = "Oceania";
        }
        else if ( c.equals( "dominican_republic" ) || c.equals( "puerto_rico" ) || c.equals( "jamaica" )
                || c.equals( "belize" ) ) {
            region = "Caribbean";
        }
        else if ( c.equals( "na" ) ) {
            region = "";
        }
        else {
            System.out.println( "ERROR: unknown country \"" + c + "\"" );
            System.exit( -1 );
        }
        if ( !ForesterUtil.isEmpty( region ) ) {
            custom_data.addProperty( new Property( VIPR_REGION, region, "", XSD_STRING, AppliesTo.NODE ) );
        }
    }

    private static void addSpecialData1( final String gb_accession,
                                         final PropertiesList custom_data,
                                         final PhylogenyNode node ) {
        if ( gb_accession.equalsIgnoreCase( "KP126912" ) || gb_accession.equalsIgnoreCase( "KP100796" )
                || gb_accession.equalsIgnoreCase( "KP744827" ) || gb_accession.equalsIgnoreCase( "KP126911" )
                || gb_accession.equalsIgnoreCase( "KP100794" ) || gb_accession.equalsIgnoreCase( "KP100792" )
                || gb_accession.equalsIgnoreCase( "KP126910" ) || gb_accession.equalsIgnoreCase( "KP100793" )
                || gb_accession.equalsIgnoreCase( "KP322752" ) || gb_accession.equalsIgnoreCase( "KY358059" )
                || gb_accession.equalsIgnoreCase( "KX685078" ) || gb_accession.equalsIgnoreCase( "KX675263" )
                || gb_accession.equalsIgnoreCase( "KX675261" ) || gb_accession.equalsIgnoreCase( "KX675262" ) ) {
            custom_data.addProperty( new Property( "vipr:AFM_from_Lit", "true", "", XSD_STRING, AppliesTo.NODE ) );
            if ( MODIFY_NODE_NAMES ) {
                node.setName( node.getName() + "*AFM_from_Lit" );
            }
        }
    }

    private static void addSpecialData12( final String gb_accession,
                                          final PropertiesList custom_data,
                                          final PhylogenyNode node ) {
        if ( gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13964_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13965_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13966_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13967_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13968_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82282_1_2347_3273.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82283_1_2352_3278.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82284_1_2347_3273.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13999_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14000_1_1_825.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14002_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14003_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14004_1_1_828.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14005_1_2354_3280.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14006_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14007_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82280_1_2350_3276.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82281_1_2348_3274.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13980_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13981_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13982_1_2354_2537.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13983_1_2354_3280.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13994_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13995_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13996_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13997_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13969_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13970_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13973_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13975_1_2342_3268.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13976_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13977_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13978_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13979_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13957_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13959_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13962_1_2346_3272.1" ) ) {
            custom_data.addProperty( new Property( "vipr:AFM_from_Lit", "true", "", XSD_STRING, AppliesTo.NODE ) );
            if ( MODIFY_NODE_NAMES ) {
                node.setName( node.getName() + "*AFM_from_Lit_2" );
            }
        }
    }

    private static void addSpecialData2( final String gb_accession,
                                         final PropertiesList custom_data,
                                         final PhylogenyNode node ) {
        if ( gb_accession.equalsIgnoreCase( "VIPR_ALG4_694265759_2332_3258.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_694265750_2353_3279.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_694265773_1905_2831.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_694265778_2367_3296.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_217316401_2309_3187.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_217316471_2313_3182.1" )
                || gb_accession.equalsIgnoreCase( "NP_740518.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_348549271_2240_3148.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_333980908_2253_3137.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_348549277_2043_2906.1" ) ) {
            custom_data.addProperty( new Property( "vipr:AFM_EXP", "positive", "", XSD_STRING, AppliesTo.NODE ) );
            if ( MODIFY_NODE_NAMES ) {
                node.setName( node.getName() + "*AFM_EXP_P" );
            }
        }
        if ( gb_accession.equalsIgnoreCase( "VIPR_ALG4_913075276_2342_3265.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_930628106_2365_3291.1" ) ) {
            custom_data.addProperty( new Property( "vipr:AFM_EXP", "negative", "", XSD_STRING, AppliesTo.NODE ) );
            if ( MODIFY_NODE_NAMES ) {
                node.setName( node.getName() + "*AFM_EXP_N" );
            }
        }
    }
}
