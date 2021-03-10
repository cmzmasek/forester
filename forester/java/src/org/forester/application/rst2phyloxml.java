
package org.forester.application;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;

public class rst2phyloxml {

    public static void main( final String args[] ) {
        final File infile = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
        final StringBuilder sb1 = new StringBuilder();
        final StringBuilder sb2 = new StringBuilder();
        try (BufferedReader br = new BufferedReader( new FileReader( infile ) )) {
            String line;
            boolean readtree1 = false;
            boolean readtree2 = false;
            while ( ( line = br.readLine() ) != null ) {
                line = line.trim();
                if ( line.length() > 0 ) {
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
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }
        try {
            p2 = factory.create( sb2.toString(), new NHXParser() )[ 0 ];
        }
        catch ( final IOException e1 ) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }
        /*
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final String node_name = node.getName();
            if ( node.isExternal() && !ForesterUtil.isEmpty( node_name ) ) {
                final int i = node_name.lastIndexOf( '_' );
                if ( i > 0 ) {
                    node.setName( node_name.substring( i + 1 ) );
                }
            }
        }
        */
        if ( p1 != null ) {
            try {
                final PhylogenyWriter w = new PhylogenyWriter();
                w.toPhyloXML( p1, 0, new File( outfile.getAbsoluteFile() + "1" ) );
            }
            catch ( final IOException e ) {
                System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
                System.exit( -1 );
            }
        }
        if ( p2 != null ) {
            try {
                final PhylogenyWriter w = new PhylogenyWriter();
                w.toPhyloXML( p2, 0, new File( outfile.getAbsoluteFile() + "2" ) );
            }
            catch ( final IOException e ) {
                System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
                System.exit( -1 );
            }
        }
    }
}
