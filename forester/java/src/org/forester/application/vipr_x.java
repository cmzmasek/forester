
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

    private static final String  XSD_STRING   = "xsd:string";
    private static final String  VIPR_HOST    = "vipr:Host";
    private static final String  VIPR_COUNTRY = "vipr:Country";
    private static final String  VIPR_YEAR    = "vipr:Year";
    private static final String  VIPR_REGION  = "vipr:Region";
    private final static String  PRG_NAME     = "vipr_x";
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
    private final static Pattern VIPR_PATTERN = Pattern
            .compile( "(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)" );
    private final static Pattern YEAR_PATTERN = Pattern.compile( "(\\d{4}).*" );

    public static void main( final String args[] ) {
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
                final Matcher m = VIPR_PATTERN.matcher( name );
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
                    custom_data.addProperty( new Property( VIPR_COUNTRY, country, "", XSD_STRING, AppliesTo.NODE ) );
                    custom_data.addProperty( new Property( VIPR_HOST, host, "", XSD_STRING, AppliesTo.NODE ) );
                    ext_node.getNodeData().setProperties( custom_data );
                    final Sequence seq = new Sequence();
                    seq.setAccession( new Accession( gb_accession, Accession.Source.NCBI ) );
                    try {
                        seq.setSymbol( gene_symbol );
                    }
                    catch ( final PhyloXmlDataFormatException e ) {
                        // ignore
                    }
                    ext_node.getNodeData().addSequence( seq );
                    String region = "";
                    final String c = country.toLowerCase();
                    if ( c.equals( "canada" ) || c.equals( "usa" ) || c.equals( "mexico" ) ) {
                        region = "North America";
                    }
                    else if ( c.equals( "denmark" ) || c.equals( "finland" ) || c.equals( "france" )
                            || c.equals( "germany" ) || c.equals( "italy" ) || c.equals( "netherlands" )
                            || c.equals( "norway" ) || c.equals( "spain" ) || c.equals( "united_kingdom" ) ) {
                        region = "Western Europe";
                    }
                    else if ( c.equals( "japan" ) || c.equals( "taiwan" ) || c.equals( "hong_kong" )
                            || c.equals( "china" ) ) {
                        region = "East Asia";
                    }
                    else if ( c.equals( "india" ) ) {
                        region = "South Asia";
                    }
                    else if ( c.equals( "malaysia" ) || c.equals( "philippines" ) || c.equals( "viet_nam" ) ) {
                        region = "Southeast Asia";
                    }
                    else if ( c.equals( "gambia" ) || c.equals( "kenya" ) || c.equals( "senegal" )
                            || c.equals( "south_africa" ) || c.equals( "tanzania" ) ) {
                        region = "Africa";
                    }
                    else if ( c.equals( "australia" ) || c.equals( "new_zealand" ) ) {
                        region = "Oceania";
                    }
                    else if ( c.equals( "na" ) ) {
                        region = "";
                    }
                    else {
                        System.out.println( "ERROR: unkknown country \"" + c + "\"" );
                        System.exit( -1 );
                    }
                    if ( !ForesterUtil.isEmpty( region ) ) {
                        custom_data.addProperty( new Property( VIPR_REGION, region, "", XSD_STRING, AppliesTo.NODE ) );
                    }
                    if ( gb_accession.equalsIgnoreCase( "KP126912" ) || gb_accession.equalsIgnoreCase( "KP100796" )
                            || gb_accession.equalsIgnoreCase( "KP744827" )
                            || gb_accession.equalsIgnoreCase( "KP126911" )
                            || gb_accession.equalsIgnoreCase( "KP100794" )
                            || gb_accession.equalsIgnoreCase( "KP100792" )
                            || gb_accession.equalsIgnoreCase( "KP126910" )
                            || gb_accession.equalsIgnoreCase( "KP100793" )
                            || gb_accession.equalsIgnoreCase( "KP322752" )
                            || gb_accession.equalsIgnoreCase( "KY358059" )
                            || gb_accession.equalsIgnoreCase( "KX685078" )
                            || gb_accession.equalsIgnoreCase( "KX675263" )
                            || gb_accession.equalsIgnoreCase( "KX675261" )
                            || gb_accession.equalsIgnoreCase( "KX675262" ) ) {
                        custom_data.addProperty( new Property( "vipr:AFM_from_Lit",
                                                               "true",
                                                               "",
                                                               XSD_STRING,
                                                               AppliesTo.NODE ) );
                    }
                }
                else {
                    System.out.println( "ERROR: name \"" + name + "\" could not be matched" );
                    System.exit( -1 );
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
