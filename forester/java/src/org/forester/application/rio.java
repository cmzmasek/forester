// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: www.phylosoft.org/forester

package org.forester.application;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Vector;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.sdi.RIO;
import org.forester.util.ForesterUtil;

public class rio {

    final static private String  PRG_NAME                         = "RIO";
    final static private String  PRG_VERSION                      = "2.03 ALPHA";
    final static private String  PRG_DATE                         = "2010.01.15";
    final static private String  E_MAIL                           = "czmasek@burnham.org";
    final static private String  WWW                              = "www.phylosoft.org/forester/";
    final static private boolean TIME                             = true;
    final static private boolean VERBOSE                          = true;
    // For method getDistances -- calculation of distances.
    final static private boolean MINIMIZE_COST                    = false;
    // For method getDistances -- calculation of distances.
    final static private boolean MINIMIZE_DUPS                    = true;
    // For method getDistances -- calculation of distances.
    final static private boolean MINIMIZE_HEIGHT                  = true;
    final static private int     WARN_NO_ORTHOS_DEFAULT           = 2;
    final static private int
                                 // How many sd away from mean to root.
                                 WARN_MORE_THAN_ONE_ORTHO_DEFAULT = 2;
    // How many sd away from mean to LCA of orthos.
    final static private double  THRESHOLD_ULTRA_PARALOGS_DEFAULT = 50;
    // How many sd away from mean to LCA of orthos.
    final static private double  WARN_ONE_ORTHO_DEFAULT           = 2;

    // Factor between the two distances to their LCA
    // (larger/smaller).
    // Factor between the two distances to their LCA
    // (larger/smaller).
    /**
     * Calculates the mean and standard deviation of all nodes of Phylogeny t
     * which have a bootstrap values zero or more. Returns null in case of
     * failure (e.g t has no bootstrap values, or just one).
     * <p>
     * 
     * @param t
     *            reference to a tree with bootstrap values
     * @return Array of doubles, [0] is the mean, [1] the standard deviation
     */
    private static double[] calculateMeanBoostrapValue( final Phylogeny t ) {
        double b = 0;
        int n = 0;
        long sum = 0;
        double x = 0.0, mean = 0.0;
        final double[] da = new double[ 2 ];
        final Vector<Double> bv = new Vector<Double>();
        PhylogenyNode node = null;
        PreorderTreeIterator i = null;
        i = new PreorderTreeIterator( t );
        // Calculates the mean.
        while ( i.hasNext() ) {
            node = i.next();
            if ( !( ( node.getParent() != null ) && node.getParent().isRoot()
                    && ( PhylogenyMethods.getConfidenceValue( node.getParent().getChildNode1() ) > 0 )
                    && ( PhylogenyMethods.getConfidenceValue( node.getParent().getChildNode2() ) > 0 ) && ( node
                    .getParent().getChildNode2() == node ) ) ) {
                b = PhylogenyMethods.getConfidenceValue( node );
                if ( b > 0 ) {
                    sum += b;
                    bv.addElement( new Double( b ) );
                    n++;
                }
            }
            // i.next();
        }
        if ( n < 2 ) {
            return null;
        }
        mean = ( double ) sum / n;
        // Calculates the standard deviation.
        sum = 0;
        for( int j = 0; j < n; ++j ) {
            b = ( bv.elementAt( j ) ).intValue();
            x = b - mean;
            sum += ( x * x );
        }
        da[ 0 ] = mean;
        da[ 1 ] = java.lang.Math.sqrt( sum / ( n - 1.0 ) );
        return da;
    }

    private final static void errorInCommandLine() {
        System.out.println( "\nrio: Error in command line.\n" );
        printHelp();
        System.exit( -1 );
    }

    public static void main( final String[] args ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE, E_MAIL, WWW );
        File species_tree_file = null;
        File multiple_trees_file = null;
        File outfile = null;
        String seq_name = "";
        String arg = "";
        boolean output_ultraparalogs = false;
        ArrayList<String> orthologs_al_for_dc = null;
        double t_orthologs = 0.0;
        double t_orthologs_dc = 0.0;
        double threshold_ultra_paralogs = 0.0;
        double[] bs_mean_sd = null;
        int sort = 13;
        Phylogeny species_tree = null;
        RIO rio_instance = null;
        PrintWriter out = null;
        long time = 0;
        if ( args.length < 2 ) {
            printHelp();
            System.exit( 0 );
        }
        else if ( ( args.length < 3 ) || ( args.length > 18 ) ) {
            errorInCommandLine();
        }
        for( final String arg2 : args ) {
            if ( arg2.trim().charAt( 0 ) != 'p' ) {
                if ( arg2.trim().length() < 3 ) {
                    errorInCommandLine();
                }
                else {
                    arg = arg2.trim().substring( 2 );
                }
            }
            try {
                switch ( arg2.trim().charAt( 0 ) ) {
                    case 'M':
                        multiple_trees_file = new File( arg );
                        break;
                    case 'N':
                        seq_name = arg;
                        break;
                    case 'S':
                        species_tree_file = new File( arg );
                        break;
                    case 'O':
                        outfile = new File( arg );
                        break;
                    case 'p':
                        output_ultraparalogs = true;
                        break;
                    case 'P':
                        sort = Integer.parseInt( arg );
                        if ( ( sort < 0 ) || ( sort > 17 ) ) {
                            errorInCommandLine();
                        }
                        break;
                    case 'L':
                        t_orthologs = Double.parseDouble( arg );
                        break;
                    case 'v':
                        threshold_ultra_paralogs = Double.parseDouble( arg );
                        break;
                    default:
                        errorInCommandLine();
                }
            }
            catch ( final Exception e ) {
                errorInCommandLine();
            }
        }
        if ( ( seq_name == "" ) || ( species_tree_file == null ) || ( multiple_trees_file == null )
                || ( outfile == null ) ) {
            errorInCommandLine();
        }
        if ( ( sort < 0 ) || ( sort > 2 ) ) {
            errorInCommandLine();
        }
        if ( VERBOSE ) {
            System.out.println( "\nMultiple trees file:                          " + multiple_trees_file );
            System.out.println( "Seq name:                                     " + seq_name );
            System.out.println( "Species tree file:                            " + species_tree_file );
            System.out.println( "Outfile:                                      " + outfile );
            System.out.println( "Sort:                                         " + sort );
            System.out.println( "Threshold orthologs:                          " + t_orthologs );
            System.out.println( "Threshold orthologs for distance calc.:       " + t_orthologs_dc );
            if ( output_ultraparalogs ) {
                System.out.println( "Threshold ultra paralogs:                     " + threshold_ultra_paralogs );
            }
        }
        if ( TIME && VERBOSE ) {
            time = System.currentTimeMillis();
        }
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            species_tree = factory.create( species_tree_file, new PhyloXmlParser() )[ 0 ];
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
        if ( !species_tree.isRooted() ) {
            ForesterUtil.printErrorMessage( PRG_NAME, "Species tree is not rooted" );
            System.exit( -1 );
        }
        if ( !species_tree.isCompletelyBinary() ) {
            ForesterUtil.printErrorMessage( PRG_NAME, "Species tree is not completely binary" );
            System.exit( -1 );
        }
        rio_instance = new RIO();
        final StringBuffer output = new StringBuffer();
        try {
            rio_instance.inferOrthologs( multiple_trees_file, species_tree.copy(), seq_name );
            output.append( rio_instance.inferredOrthologsToString( seq_name, sort, t_orthologs ) );
            if ( output_ultraparalogs ) {
                output.append( "\n\nUltra paralogs:\n" );
                output.append( rio_instance.inferredUltraParalogsToString( seq_name, threshold_ultra_paralogs ) );
            }
            output.append( "\n\nSort priority: " + RIO.getOrder( sort ) );
            output.append( "\nExt nodes    : " + rio_instance.getExtNodesOfAnalyzedGeneTrees() );
            output.append( "\nSamples      : " + rio_instance.getNumberOfSamples() + "\n" );
            out = new PrintWriter( new FileWriter( outfile ), true );
        }
        catch ( final Exception e ) {
            ForesterUtil.printErrorMessage( PRG_NAME, e.getLocalizedMessage() );
            e.printStackTrace();
            System.exit( -1 );
        }
        out.println( output );
        out.close();
        ForesterUtil.programMessage( PRG_NAME, "wrote results to \"" + outfile + "\"" );
        if ( TIME && VERBOSE ) {
            time = System.currentTimeMillis() - time;
            ForesterUtil.programMessage( PRG_NAME, "time: " + time + "ms" );
        }
        ForesterUtil.programMessage( PRG_NAME, "OK." );
        System.exit( 0 );
    }

    private final static void printHelp() {
        System.out.println( "M= (String) Multiple gene tree file (mandatory)" );
        System.out.println( "N= (String) Query sequence name (mandatory)" );
        System.out.println( "S= (String) Species tree file (mandatory)" );
        System.out.println( "O= (String) Output file name -- overwritten without warning! (mandatory)" );
        System.out.println( "P= (int)    Sort priority" );
        System.out.println( "L= (double) Threshold orthologs for output" );
        System.out.println( " Sort priority (\"P=\"):" );
        System.out.println( RIO.getOrderHelp().toString() );
        System.out.println();
        System.out
                .println( " Example: \"rio M=gene_trees.xml N=bcl2_NEMVE S=species_tree.xml D=distances P=13 p O=out\"" );
        System.out.println();
    }
}
