
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
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class voc_prep {

    private static final Pattern YEAR_MONTH_PATTERN = Pattern.compile( "\\W(\\d{4})_(\\d{1,2})_\\d{1,2}\\W" );
    private static final String  XSD_STRING         = "xsd:string";
    private static final String  VIPR_YEAR_MONTH    = "vipr:Year_Month";
    private final static String  PRG_NAME           = "voc_prep";
    private static final String  PRG_DATE           = "2021-05-24";
    private static final String  PRG_VERSION        = "1.0.0";
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        final File infile = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
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
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final String node_name = node.getName();
            if ( node.isExternal() && !ForesterUtil.isEmpty( node_name ) ) {
                // MT412292|SARS_CoV_2/human/USA/WA_UW_6170/2020|2020_04_01|Human|USA|NA|Severe_acute_respiratory_syndrome_related_coronavirus
                final String[] x = node_name.split( "\\|" );
                if ( x.length > 1 ) {
                    final String new_node_name = x[ 0 ];
                    System.out.println( node_name + " -> " + new_node_name );
                    node.setName( new_node_name );
                }
                else {
                    System.out.println( "WARNING: Odd format: " + node_name );
                }
                final Matcher m = YEAR_MONTH_PATTERN.matcher( node_name );
                if ( m.find() ) {
                    final int year = ( Integer.valueOf( m.group( 1 ) ) );
                    final int month = ( Integer.valueOf( m.group( 2 ) ) );
                    if ( ( month < 1 ) || ( month > 12 ) ) {
                        System.out.println( "Error: Illegal month: " + month );
                        System.exit( -1 );
                    }
                    if ( ( year < 1000 ) || ( year > 2022 ) ) {
                        System.out.println( "Error: Illegal year: " + year );
                        System.exit( -1 );
                    }
                    System.out.println( "year " + year );
                    System.out.println( "month " + month );
                    PropertiesList custom_data = node.getNodeData().getProperties();
                    if ( custom_data == null ) {
                        custom_data = new PropertiesList();
                    }
                    final StringBuilder sb = new StringBuilder();
                    sb.append( year );
                    if ( month < 10 ) {
                        sb.append( "_0" );
                    }
                    else {
                        sb.append( "_" );
                    }
                    sb.append( month );
                    custom_data.addProperty( new Property( VIPR_YEAR_MONTH,
                                                           sb.toString(),
                                                           "",
                                                           XSD_STRING,
                                                           AppliesTo.NODE ) );
                    node.getNodeData().setProperties( custom_data );
                }
                else {
                    System.out.println( "WARNING: No date information in: " + node_name );
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
