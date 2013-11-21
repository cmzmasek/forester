// javac -cp ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/subtree_feature_count.java
//
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.subtree_feature_count

package org.forester.applications;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class subtree_feature_count {

    final static private String MIN_DISTANCE_TO_ROOT_OPTION = "d";
    final static private String E_MAIL                      = "phylosoft@gmail.com";
    final static private String HELP_OPTION_1               = "help";
    final static private String HELP_OPTION_2               = "h";
    final static private String PRG_DATE                    = "131120";
    final static private String PRG_DESC                    = "";
    final static private String PRG_NAME                    = "subtree_feature_count";
    final static private String PRG_VERSION                 = "0.90";
    final static private String WWW                         = "sites.google.com/site/cmzmasek/home/software/forester";

    public static void main( final String args[] ) {
        try {
            final CommandLineArguments cla = new CommandLineArguments( args );
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length < 2 ) ) {
                printHelp();
                System.exit( 0 );
            }
            final List<String> allowed_options = new ArrayList<String>();
            allowed_options.add( MIN_DISTANCE_TO_ROOT_OPTION );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            final double min_distance_to_root = cla.getOptionValueAsDouble( MIN_DISTANCE_TO_ROOT_OPTION );
            if ( min_distance_to_root <= 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "attempt to use min distance to root of: " + min_distance_to_root );
            }
            final File intree_file = cla.getFile( 0 );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny phy = factory.create( intree_file, PhyloXmlParser.createPhyloXmlParserXsdValidating() )[ 0 ];
            execute( phy, min_distance_to_root );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    private final static void execute( final Phylogeny phy, final double min_distance_to_root ) {
        final List<List<PhylogenyNode>> ll = PhylogenyMethods.divideIntoSubTrees( phy, min_distance_to_root );
        for( final List<PhylogenyNode> l : ll ) {
            int xray = 0;
            int nmr = 0;
            int model = 0;
            for( final PhylogenyNode node : l ) {
                if ( node.getNodeData().isHasSequence() ) {
                    final Sequence seq = node.getNodeData().getSequence();
                    final SortedSet<Accession> xrefs = seq.getCrossReferences();
                    if ( !ForesterUtil.isEmpty( xrefs ) ) {
                        for( final Accession xref : xrefs ) {
                            if ( xref.getSource().equalsIgnoreCase( "pdb" ) ) {
                                if ( xref.getComment().equalsIgnoreCase( "x-ray" )
                                        || xref.getComment().equalsIgnoreCase( "xray" ) ) {
                                    ++xray;
                                }
                                if ( xref.getComment().equalsIgnoreCase( "nmr" ) ) {
                                    ++nmr;
                                }
                                if ( xref.getComment().equalsIgnoreCase( "model" ) ) {
                                    ++model;
                                }
                            }
                        }
                    }
                }
            }
            final int n = l.size();
            final double xray_p = ForesterUtil.round( 100.0 * xray / n, 1 );
            final double nmr_p = ForesterUtil.round( 100.0 * nmr / n, 1 );
            final double model_p = ForesterUtil.round( 100.0 * model / n, 1 );
            final StringBuilder sb = new StringBuilder();
            sb.append( String.valueOf( n ) );
            sb.append( "\t" );
            sb.append( String.valueOf( xray ) );
            sb.append( "\t" );
            sb.append( String.valueOf( nmr ) );
            sb.append( "\t" );
            sb.append( String.valueOf( model ) );
            sb.append( "\t" );
            sb.append( String.valueOf( xray_p ) );
            sb.append( "\t" );
            sb.append( String.valueOf( nmr_p ) );
            sb.append( "\t" );
            sb.append( String.valueOf( model_p ) );
            System.out.println( sb );
        }
    }

    private static void printHelp() {
        ForesterUtil.printProgramInformation( PRG_NAME,
                                              PRG_DESC,
                                              PRG_VERSION,
                                              PRG_DATE,
                                              E_MAIL,
                                              WWW,
                                              ForesterUtil.getForesterLibraryInformation() );
        System.out.print( "Usage: " );
        System.out.println( PRG_NAME + " -d=<min distance to root> <intree>" );
        System.out.println();
        System.out.println();
    }
}
