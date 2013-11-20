
package org.forester.application;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class subtree_feature_count {

    final static private String DEPTH_OPTION  = "d";
    final static private String E_MAIL        = "phylosoft@gmail.com";
    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";
    final static private String PRG_DATE      = "131120";
    final static private String PRG_DESC      = "";
    final static private String PRG_NAME      = "subtree_feature_count";
    final static private String PRG_VERSION   = "0.90";
    final static private String WWW           = "sites.google.com/site/cmzmasek/home/software/forester";

    public static void main( final String args[] ) {
        try {
            final CommandLineArguments cla = new CommandLineArguments( args );
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length < 2 ) ) {
                printHelp();
                System.exit( 0 );
            }
            final List<String> allowed_options = new ArrayList<String>();
            allowed_options.add( DEPTH_OPTION );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            final double depth = cla.getOptionValueAsDouble( DEPTH_OPTION );
            final File intree_file = cla.getFile( 0 );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny phy = factory.create( intree_file, PhyloXmlParser.createPhyloXmlParserXsdValidating() )[ 0 ];
            execute( phy, depth );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    private static StringBuilder analyzeSubtree( final PhylogenyNode n, final double depth ) {
        final PhylogenyNode node = moveUp( n, depth );
        final List<PhylogenyNode> ext_descs = node.getAllExternalDescendants();
        for( final PhylogenyNode ext : ext_descs ) {
            if ( ext.getIndicator() != 0 ) {
                throw new RuntimeException( "this should not have happened" );
            }
            ext.setIndicator( ( byte ) 1 );
        }
        int xray = 0;
        int nmr = 0;
        int model = 0;
        boolean is_first = true;
        PhylogenyNode first_node = null;
        PhylogenyNode last_node = null;
        int c = 0;
        for( final PhylogenyNode ext : ext_descs ) {
            if ( is_first ) {
                first_node = ext;
                is_first = false;
            }
            last_node = ext;
            ++c;
            if ( ext.getNodeData().isHasSequence() ) {
                final Sequence seq = ext.getNodeData().getSequence();
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
        final StringBuilder sb = new StringBuilder();
        sb.append( String.valueOf( c ) );
        sb.append( "\t" );
        sb.append( first_node.getName() );
        sb.append( "\t" );
        sb.append( last_node.getName() );
        sb.append( "\t" );
        sb.append( String.valueOf( xray ) );
        sb.append( "\t" );
        sb.append( String.valueOf( nmr ) );
        sb.append( "\t" );
        sb.append( String.valueOf( model ) );
        return sb;
    }

    private static void execute( final Phylogeny phy, final double depth ) {
        setAllIndicatorsToZero( phy );
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getIndicator() != 0 ) {
                continue;
            }
            final StringBuilder s = analyzeSubtree( n, depth );
            System.out.println( s.toString() );
        }
    }

    private static PhylogenyNode moveUp( final PhylogenyNode node, final double depth ) {
        PhylogenyNode n = node;
        PhylogenyNode prev = node;
        while ( depth < n.calculateDistanceToRoot() ) {
            prev = n;
            n = n.getParent();
        }
        return prev;
    }

    private static void printHelp() {
        ForesterUtil.printProgramInformation( PRG_NAME,
                                              PRG_DESC,
                                              PRG_VERSION,
                                              PRG_DATE,
                                              E_MAIL,
                                              WWW,
                                              ForesterUtil.getForesterLibraryInformation() );
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + "" );
        System.out.println();
        System.out.println( " example: " );
        System.out.println();
        System.out.println();
        System.out.println();
    }

    private static void setAllIndicatorsToZero( final Phylogeny phy ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
            it.next().setIndicator( ( byte ) 0 );
        }
    }
}
