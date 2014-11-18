
package org.forester.application;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class dom_dup {

    // HUMAN MOUSE
    // ARATH SOYBN VOLCA CYAME PARTE THAPS EMIHU NAEGR
    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";
    final static private String PRG_NAME      = "dom_dup";
    final static private String PRG_DESC      = "";
    final static private String PRG_VERSION   = "0.90";
    final static private String PRG_DATE      = "2013.03.12";
    final static private String E_MAIL        = "phylosoft@gmail.com";
    final static private String WWW           = "sites.google.com/site/cmzmasek/home/software/forester";

    public static void main( final String args[] ) {
        try {
            final CommandLineArguments cla = new CommandLineArguments( args );
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( cla.getNumberOfNames() != 3 ) ) {
                printHelp();
                System.exit( 0 );
            }
            final String pattern_str = cla.getName( 0 );
            final File intree_file = cla.getFile( 2 );
            final File species_groups_file = cla.getFile( 1 );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny phy = factory.create( intree_file, PhyloXmlParser.createPhyloXmlParserXsdValidating() )[ 0 ];
            ForesterUtil.programMessage( PRG_NAME, "Pattern string: " + pattern_str );
            final Pattern pattern = Pattern.compile( pattern_str );
            ForesterUtil.programMessage( PRG_NAME, "Pattern is: " + pattern );
            final SortedSet<String> set_a = new TreeSet<String>();
            final SortedSet<String> set_b = new TreeSet<String>();
            read( species_groups_file, set_a, set_b );
            print_set( set_a, "Set a:" );
            print_set( set_b, "Set b:" );
            final SortedSet<String> matching_names = obtainMatchingNames( phy, pattern );
            ForesterUtil.programMessage( PRG_NAME, "Found names: " );
            final SortedMap<String, List<String>> pairs = obtainPairs( matching_names );
            int lca_counter = 0;
            int non_lca_counter = 0;
            int missing_counter = 0;
            int total_counter = 0;
            final Iterator<Entry<String, List<String>>> it = pairs.entrySet().iterator();
            while ( it.hasNext() ) {
                final Map.Entry<String, List<String>> x = it.next();
                total_counter++;
                if ( x.getValue().size() == 2 ) {
                    final String a = x.getValue().get( 0 );
                    final String b = x.getValue().get( 1 );
                    System.out.print( a + " - " + b );
                    final PhylogenyNode lca = PhylogenyMethods.calculateLCA( phy.getNode( a ), phy.getNode( b ) );
                    final List<PhylogenyNode> external_descs = lca.getAllExternalDescendants();
                    boolean in_a = false;
                    boolean in_b = false;
                    for( final PhylogenyNode external_desc : external_descs ) {
                        final String tc = external_desc.getNodeData().getTaxonomy().getTaxonomyCode();
                        if ( set_a.contains( tc ) ) {
                            in_a = true;
                        }
                        if ( set_b.contains( tc ) ) {
                            in_b = true;
                        }
                    }
                    if ( in_a && in_b ) {
                        System.out.print( " => LCA " );
                        lca_counter++;
                    }
                    else {
                        non_lca_counter++;
                    }
                    System.out.println();
                }
                else if ( x.getValue().size() == 1 ) {
                    System.out.println( x.getValue().get( 0 ) + " => no partner in current tree!" );
                    missing_counter++;
                }
                else {
                    System.out.println( "error" );
                    System.exit( -1 );
                }
            }
            System.out.println( "Total       : " + total_counter );
            System.out.println( "LCA         : " + lca_counter );
            System.out.println( "Non-LCA     : " + non_lca_counter );
            System.out.println( "With missing: " + missing_counter );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    private static SortedMap<String, List<String>> obtainPairs( final SortedSet<String> matching_names ) {
        final SortedMap<String, List<String>> pairs = new TreeMap<String, List<String>>();
        for( final String m : matching_names ) {
            final String short_m = m.substring( 0, m.indexOf( '~' ) );
            if ( !pairs.containsKey( short_m ) ) {
                final List<String> p = new ArrayList<String>();
                p.add( m );
                pairs.put( short_m, p );
            }
            else {
                pairs.get( short_m ).add( m );
            }
        }
        return pairs;
    }

    private static SortedSet<String> obtainMatchingNames( final Phylogeny phy, final Pattern pattern ) {
        final SortedSet<String> matching_names = new TreeSet<String>();
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final Matcher m = pattern.matcher( n.getName() );
            if ( m.find() ) {
                matching_names.add( n.getName() );
            }
        }
        return matching_names;
    }

    private static void print_set( final Set<String> set_a, final String l ) {
        ForesterUtil.programMessage( PRG_NAME, l );
        for( final String s : set_a ) {
            System.out.print( s + " " );
        }
        System.out.println();
    }

    private static void read( final File species_groups_file, final Set<String> set_a, final Set<String> set_b )
            throws IOException {
        final BufferedReader reader = ForesterUtil.obtainReader( species_groups_file );
        String line;
        boolean first_line = true;
        while ( ( line = reader.readLine() ) != null ) {
            line = line.trim();
            if ( !ForesterUtil.isEmpty( line ) ) {
                final String s[] = line.split( " " );
                for( final String name : s ) {
                    if ( first_line ) {
                        set_a.add( name );
                    }
                    else {
                        set_b.add( name );
                    }
                }
                if ( first_line ) {
                    first_line = false;
                }
            }
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
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + "" );
        System.out.println();
        System.out.println( " example: " );
        System.out.println();
        System.out
        .println( "dom_dup \"HUMAN~[12]-2\" groups.txt RRMa_ALL_plus_RRMa_ee3_50_hmmalign_05_40_fme_gsdi.phylo.xml" );
        System.out.println();
        System.out.println();
    }
}
