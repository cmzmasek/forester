
package org.forester.application;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class rst2phyloxml {

    final static private Pattern p1            = Pattern.compile( "Branch\\s+(\\d+):\\s+" );
    final static private Pattern p2            = Pattern
            .compile( "\\s+(\\d+)\\s+([A-Z-_]+)\\s+(\\S+)\\s->\\s+([A-Z-_]+)" );
    private static final String  XSD_STRING    = "xsd:string";
    private static final String  VIPR_MUTATION = "vipr:Mutation";
    public static void main( final String args[] ) {
        final File infile = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
        final StringBuilder sb1 = new StringBuilder();
        final StringBuilder sb2 = new StringBuilder();
        final Map<String, List<String[]>> changes = new HashMap<>();
        String current_branch = null;
        try (BufferedReader br = new BufferedReader( new FileReader( infile ) )) {
            String line;
            boolean readtree1 = false;
            boolean readtree2 = false;
            while ( ( line = br.readLine() ) != null ) {
                if ( line.length() > 0 ) {
                    final Matcher m1 = p1.matcher( line );
                    final Matcher m2 = p2.matcher( line );
                    if ( line.contains( "Ancestral reconstruction" ) ) {
                        readtree1 = true;
                        readtree2 = false;
                    }
                    else if ( line.contains( "tree with node labels for Rod Page" ) ) {
                        readtree1 = false;
                        readtree2 = true;
                    }
                    else if ( line.startsWith( "Nodes " ) ) {
                        readtree1 = false;
                        readtree2 = false;
                    }
                    else if ( readtree1 ) {
                        sb1.append( line );
                        readtree1 = false;
                    }
                    else if ( readtree2 ) {
                        sb2.append( line );
                        readtree2 = false;
                    }
                    else if ( m1.find() ) {
                        readtree1 = false;
                        readtree2 = false;
                        final String branch = m1.group( 1 );
                        System.out.println( line );
                        System.out.println( "  Branch : " + branch );
                        current_branch = branch;
                    }
                    else if ( m2.find() ) {
                        readtree1 = false;
                        readtree2 = false;
                        final String position = m2.group( 1 );
                        final String from = m2.group( 2 );
                        final String support = m2.group( 3 );
                        final String to = m2.group( 4 );
                        System.out.println( line );
                        System.out.println( "  Position: " + position );
                        System.out.println( "  From    : " + from );
                        System.out.println( "  Support : " + support );
                        System.out.println( "  To      : " + to );
                        System.out.println();
                        if ( !changes.containsKey( current_branch ) ) {
                            changes.put( current_branch, new ArrayList<String[]>() );
                        }
                        changes.get( current_branch ).add( new String[] { position, from, support, to } );
                    }
                }
            }
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        Phylogeny p1 = null;
        Phylogeny p2 = null;
        try {
            p1 = factory.create( sb1.toString(), new NHXParser() )[ 0 ];
        }
        catch ( final IOException e1 ) {
            e1.printStackTrace();
        }
        try {
            p2 = factory.create( sb2.toString(), new NHXParser() )[ 0 ];
        }
        catch ( final IOException e1 ) {
            e1.printStackTrace();
        }
        X: for( final PhylogenyNodeIterator iter1 = p1.iteratorPreorder(); iter1.hasNext(); ) {
            final PhylogenyNode n1 = iter1.next();
            if ( n1.getDistanceToParent() >= 0 ) {
                final Set<String> s1 = new HashSet<>();
                s1.addAll( n1.getAllExternalDescendantsNames() );
                for( final PhylogenyNodeIterator iter2 = p2.iteratorPreorder(); iter2.hasNext(); ) {
                    final PhylogenyNode n2 = iter2.next();
                    final Set<String> s2 = new HashSet<>();
                    final List<String> e2 = n2.getAllExternalDescendantsNames();
                    for( final String name : e2 ) {
                        s2.add( name.substring( name.indexOf( '_' ) + 1, name.length() ) );
                    }
                    if ( s2.containsAll( s1 ) && s1.containsAll( s2 ) ) {
                        n2.setDistanceToParent( n1.getDistanceToParent() );
                        continue X;
                    }
                }
                System.out.println( "\nFailure find\n" + s1 );
                System.exit( -1 );
            }
        }
        for( final PhylogenyNodeIterator iter2 = p2.iteratorPreorder(); iter2.hasNext(); ) {
            final PhylogenyNode n2 = iter2.next();
            if ( changes.containsKey( n2.getName() ) ) {
                final List<String[]> c = changes.get( n2.getName() );
                for( final String[] cn : c ) {
                    PropertiesList custom_data = n2.getNodeData().getProperties();
                    if ( custom_data == null ) {
                        custom_data = new PropertiesList();
                    }
                    final String c_str = cn[ 1 ] + cn[ 0 ] + cn[ 3 ];
                    custom_data.addProperty( new Property( VIPR_MUTATION, c_str, "", XSD_STRING, AppliesTo.NODE ) );
                    n2.getNodeData().setProperties( custom_data );
                }
            }
        }
        if ( p1 != null ) {
            try {
                final PhylogenyWriter w = new PhylogenyWriter();
                w.toPhyloXML( p1, 0, new File( outfile.getAbsoluteFile() + "_1.xml" ) );
            }
            catch ( final IOException e ) {
                System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
                System.exit( -1 );
            }
        }
        if ( p2 != null ) {
            try {
                final PhylogenyWriter w = new PhylogenyWriter();
                w.toPhyloXML( p2, 0, new File( outfile.getAbsoluteFile() + "_2.xml" ) );
            }
            catch ( final IOException e ) {
                System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
                System.exit( -1 );
            }
        }
    }
}
