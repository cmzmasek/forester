
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

public class node_name_pars {

    private static int           OPERATION          = 1;
    private static final boolean MODIFY_NODE_NAMES  = true;
    private static final Pattern ARENAVIRUS_PATTERN = Pattern.compile( "^([a-zA-Z\\s]+?virus)" );
    private static final String  XSD_STRING         = "xsd:string";
    private static final String  HOST               = "annotation:Host";
    private static final String  COUNTRY            = "annotation:Country";
    private static final String  YEAR               = "annotation:Year";
    private static final String  STRAIN             = "annotation:Strain";
    private static final String  SPECIES            = "annotation:Species";
    private static final String  GENUS              = "annotation:Genus";
    private static final String  REGION             = "annotation:Region";
    private static final String  PRG_NAME           = "node_name_pars";
    private static final String  PRG_DATE           = "2023-04-13";
    private static final String  PRG_VERSION        = "1.0.1";
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
                    final String[] x = name.split( "\\|", -1 );
                    // genome|species|strain|segment|genbankacc|year|country|host
                    if ( x.length != 8 ) {
                        System.out.println( "ERROR: name \"" + name + "\" unexpected format" + x );
                        System.out.println();
                        System.exit( -1 );
                    }
                    final String genome = x[ 0 ];
                    final String genus = x[ 1 ];
                    String strain = x[ 2 ];
                    final String segment = x[ 3 ];
                    final String genbankacc = x[ 4 ];
                    String year = x[ 5 ];
                    String country = x[ 6 ];
                    String host = x[ 7 ];
                    String species = "";
                    if ( strain.equalsIgnoreCase( "not applicable" ) ) {
                        strain = "";
                    }
                    if ( year.equalsIgnoreCase( "not applicable" ) ) {
                        year = "";
                    }
                    if ( country.equalsIgnoreCase( "not applicable" ) ) {
                        country = "";
                    }
                    if ( host.equalsIgnoreCase( "not applicable" ) ) {
                        host = "";
                    }
                    final Matcher m = ARENAVIRUS_PATTERN.matcher( name );
                    if ( m.find() ) {
                        species = m.group( 1 );
                    }
                    System.out.println( name );
                    System.out.println( "genome           : " + genome );
                    System.out.println( "genus            : " + genus );
                    System.out.println( "species          : " + species );
                    System.out.println( "strain           : " + strain );
                    System.out.println( "segment          : " + segment );
                    System.out.println( "genbankacc       : " + genbankacc );
                    System.out.println( "year             : " + year );
                    System.out.println( "country          : " + country );
                    System.out.println( "host             : " + host );
                    System.out.println();
                    final PropertiesList custom_data = new PropertiesList();
                    custom_data.addProperty( new Property( GENUS, genus, "", XSD_STRING, AppliesTo.NODE ) );
                    custom_data.addProperty( new Property( SPECIES, species, "", XSD_STRING, AppliesTo.NODE ) );
                    custom_data.addProperty( new Property( STRAIN, strain, "", XSD_STRING, AppliesTo.NODE ) );
                    custom_data.addProperty( new Property( YEAR, year, "", XSD_STRING, AppliesTo.NODE ) );
                    custom_data.addProperty( new Property( COUNTRY, country, "", XSD_STRING, AppliesTo.NODE ) );
                    custom_data.addProperty( new Property( HOST, host, "", XSD_STRING, AppliesTo.NODE ) );
                    ext_node.getNodeData().setProperties( custom_data );
                    Sequence seq = null;
                    if ( ext_node.isHasNodeData() && ext_node.getNodeData().isHasSequence() ) {
                        seq = ext_node.getNodeData().getSequence();
                    }
                    else {
                        seq = new Sequence();
                    }
                    seq.setAccession( new Accession( genbankacc, Accession.Source.NCBI ) );
                    ext_node.getNodeData().addSequence( seq );
                    addRegion( country, custom_data );
                    if ( MODIFY_NODE_NAMES ) {
                        ext_node.setName( genome );
                    }
                }
            }
            else {
                System.out.println( "WARNING: no name for node" );
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
                || c.equals( "brazil" ) || c.equals( "argentina" ) || c.equals( "uruguay" ) || c.equals( "venezuela" )
                || c.equals( "bolivia" ) || c.equals( "paraguay" ) || c.equals( "french guiana" ) ) {
            region = "South America";
        }
        else if ( c.equals( "costa rica" ) || c.equals( "guatemala" ) ) {
            region = "Central American";
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
        else if ( c.equals( "poland" ) || c.equals( "slovakia" ) || c.equals( "czech republic" )
                || c.equals( "hungary" ) ) {
            region = "Central Europe";
        }
        else if ( c.equals( "bulgaria" ) ) {
            region = "Southeast Europe";
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
        else if ( c.equals( "jordan" ) || c.equals( "bahrain" ) || c.equals( "iraq" ) || c.equals( "saudi_arabia" )
                || c.equals( "turkey" ) || c.equals( "egypt" ) || c.equals( "israel" ) || c.equals( "west_bank" )
                || c.equals( "iran" ) ) {
            region = "West Asia";
        }
        else if ( c.equals( "india" )
                || ( ( c.equals( "pakistan" ) | c.equals( "bangladesh" ) ) || c.equals( "sri_lanka" ) ) ) {
            region = "South Asia";
        }
        else if ( c.equals( "cambodia" ) || c.equals( "thailand" ) || c.equals( "malaysia" )
                || c.equals( "philippines" ) || c.equals( "viet_nam" ) || c.equals( "myanmar" )
                || c.equals( "timor_leste" ) ) {
            region = "Southeast Asia";
        }
        else if ( c.equals( "morocco" ) || c.equals( "gambia" ) || c.equals( "kenya" ) || c.equals( "senegal" )
                || c.equals( "south_africa" ) || c.equals( "tanzania" ) || c.equals( "ghana" ) || c.equals( "benin" )
                || c.equals( "tunisia" ) || c.equals( "nigeria" ) || c.equals( "libya" ) || c.equals( "djibouti" )
                || c.equals( "sierra_leone" ) || c.equals( "guinea" ) || c.equals( "botswana" )
                || c.equals( "ethiopia" ) || c.equals( "malawi" ) || c.equals( "mali" ) || c.equals( "zambia" )
                || c.equals( "cameroon" ) || c.equals( "angola" ) || c.equals( "namibia" ) || c.equals( "south africa" )
                || c.equals( "sierra leone" ) || c.equals( "liberia" ) || c.equals( "togo" )
                || c.equals( "central african republic" ) || c.equals( "mozambique" ) ) {
            region = "Africa";
        }
        else if ( c.equals( "australia" ) || c.equals( "new_zealand" ) ) {
            region = "Oceania";
        }
        else if ( c.equals( "dominican_republic" ) || c.equals( "puerto_rico" ) || c.equals( "jamaica" )
                || c.equals( "belize" )  || c.equals( "trinidad and tobago" ) ) {
            region = "Caribbean";
        }
        else if ( c.equals( "na" ) || ( c.length() < 1 ) ) {
            region = "";
        }
        else {
            System.out.println( "ERROR: unknown country \"" + c + "\"" );
            System.exit( -1 );
        }
        if ( !ForesterUtil.isEmpty( region ) ) {
            custom_data.addProperty( new Property( REGION, region, "", XSD_STRING, AppliesTo.NODE ) );
        }
    }
}
