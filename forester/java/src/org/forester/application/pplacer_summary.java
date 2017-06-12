
package org.forester.application;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.BasicDescriptiveStatistics;

public class pplacer_summary {

    public static void main( final String args[] ) {
        final File indir = new File( "." );
        final File[] list_of_files = indir.listFiles();
        final List<File> infiles = new ArrayList<>();
        for( final File file : list_of_files ) {
            if ( file.isFile() && file.canRead() && file.toString().endsWith( ".sing.tre" ) ) {
                infiles.add( file );
            }
        }
        Collections.sort( infiles );
        final BasicDescriptiveStatistics non_unique_placements_stats = new BasicDescriptiveStatistics();
        final BasicDescriptiveStatistics unexpected_top_placements_stats = new BasicDescriptiveStatistics();
        final BasicDescriptiveStatistics unexpected_unique_top_placements_stats = new BasicDescriptiveStatistics();
        for( final File infile : infiles ) {
            Phylogeny phys[] = null;
            try {
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                phys = factory.create( infile,
                                       org.forester.io.parsers.util.ParserUtils
                                               .createParserDependingOnFileType( infile, true ) );
            }
            catch ( final Exception e ) {
                System.out.println( "Could not read " + infile + ":" + e );
                System.exit( -1 );
            }
            final Pattern p = Pattern.compile( "_Q_#\\d+_M=(.+)" );
            int total_trees = 0;
            int total_placements = 0;
            int total_expected_top_placements = 0;
            int total_unexpected_top_placements = 0;
            int total_unexpected_unique_top_placements = 0;
            int unique_placements = 0;
            int non_unique_placements = 0;
            int m1_placements = 0;
            int non_m1_placements = 0;
            for( final Phylogeny phy : phys ) {
                ++total_trees;
                final List<PhylogenyNode> nodes = phy.getNodes( p );
                if ( nodes.isEmpty() ) {
                    System.out.println();
                    System.out.println( "not found!" );
                    System.exit( -1 );
                }
                else {
                    total_placements++;
                    if ( nodes.size() == 1 ) {
                        unique_placements++;
                    }
                    else {
                        non_unique_placements++;
                    }
                    final PhylogenyNode n = nodes.get( 0 );
                    if ( n.isExternal() && !n.isRoot() ) {
                        final Matcher m = p.matcher( n.getName() );
                        if ( m.find() ) {
                            final double M = Double.parseDouble( m.group( 1 ) );
                            if ( M > 0.999999 ) {
                                m1_placements++;
                            }
                            else {
                                non_m1_placements++;
                            }
                        }
                        else {
                            System.out.println();
                            System.out.println( "no match!" );
                            System.exit( -1 );
                        }
                        String sib = null;
                        if ( n.getChildNodeIndex() == 0 ) {
                            sib = n.getParent().getChildNode2().getName();
                        }
                        else if ( n.getChildNodeIndex() == 1 ) {
                            sib = n.getParent().getChildNode1().getName();
                        }
                        else {
                            System.out.println();
                            System.out.println( "more than two children!" );
                            System.exit( -1 );
                        }
                        //  System.out.println( n.getName() + "->" + sib );
                        if ( n.getName().startsWith( sib ) ) {
                            total_expected_top_placements++;
                        }
                        else {
                            total_unexpected_top_placements++;
                            if ( nodes.size() == 1 ) {
                                total_unexpected_unique_top_placements++;
                            }
                        }
                    }
                }
            }
            System.out.println();
            System.out.println( infile.getName() );
            final Pattern pa = Pattern.compile( "(\\d+)-\\d+" );
            final Matcher m = pa.matcher( infile.getName() );
            if ( m.find() ) {
                final int start = Integer.parseInt( m.group( 1 ) );
                //      System.out.println( start + "\t"
                //              + ( ( ( double ) total_unexpected_top_placements ) / total_placements ) );
            }
            System.out.println( "total trees" + "\t" + total_trees );
            System.out.println( "total placements" + "\t" + total_placements );
            System.out.println( "total expected top placements" + "\t" + total_expected_top_placements );
            System.out.println( "total un-expected top placements" + "\t" + total_unexpected_top_placements );
            System.out.println( "total un-expected unique placements" + "\t" + total_unexpected_unique_top_placements );
            System.out.println( "unique placements" + "\t" + unique_placements );
            System.out.println( "non unique placements" + "\t" + non_unique_placements );
            System.out.println( "m1 placements" + "\t" + m1_placements );
            System.out.println( "non m1 placements" + "\t" + non_m1_placements );
            non_unique_placements_stats.addValue( ( ( double ) non_unique_placements ) / total_placements );
            unexpected_top_placements_stats
                    .addValue( ( ( double ) total_unexpected_top_placements ) / total_placements );
            unexpected_unique_top_placements_stats
                    .addValue( ( ( double ) total_unexpected_unique_top_placements ) / total_placements );
        }
        System.out.println( "Non-unique placements: Mean\t" + non_unique_placements_stats.arithmeticMean() );
        System.out.println( "Non-unique placements: Min\t" + non_unique_placements_stats.getMin() );
        System.out.println( "Non-unique placements: Max\t" + non_unique_placements_stats.getMax() );
        System.out.println( "Unexpected top-placements: Mean\t" + unexpected_top_placements_stats.arithmeticMean() );
        System.out.println( "Unexpected top-placements: Min\t" + unexpected_top_placements_stats.getMin() );
        System.out.println( "Unexpected top-placements: Max\t" + unexpected_top_placements_stats.getMax() );
        System.out.println( "Unexpected unique placements: Mean\t"
                + unexpected_unique_top_placements_stats.arithmeticMean() );
        System.out.println( "Unexpected unique placements: Min\t" + unexpected_unique_top_placements_stats.getMin() );
        System.out.println( "Unexpected unique placements: Max\t" + unexpected_unique_top_placements_stats.getMax() );
    }
}
