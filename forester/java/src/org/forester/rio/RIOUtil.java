
package org.forester.rio;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.datastructures.IntMatrix;
import org.forester.io.parsers.IteratingPhylogenyParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.rio.RIO.REROOTING;
import org.forester.sdi.SDIException;
import org.forester.sdi.SDIutil.ALGORITHM;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.EasyWriter;
import org.forester.util.ForesterUtil;

public final class RIOUtil {

    public static final void executeAnalysis( final File gene_trees_file,
                                              final File species_tree_file,
                                              final File orthology_outtable,
                                              final File orthology_groups_outfile,
                                              final File logfile,
                                              final String outgroup,
                                              final REROOTING rerooting,
                                              final int gt_first,
                                              final int gt_last,
                                              final File return_species_tree,
                                              final File return_min_dup_gene_tree,
                                              final File return_median_dup_gene_tree,
                                              final boolean transfer_taxonomy,
                                              final ALGORITHM algorithm,
                                              final boolean use_gene_trees_dir,
                                              final EasyWriter log,
                                              final double ortholog_group_cutoff ) {
        try {
            final RIO rio;
            boolean iterating = false;
            final PhylogenyParser p = ParserUtils.createParserDependingOnFileType( gene_trees_file, true );
            if ( p instanceof PhyloXmlParser ) {
                rio = RIO.executeAnalysis( gene_trees_file,
                                           species_tree_file,
                                           algorithm,
                                           rerooting,
                                           outgroup,
                                           gt_first,
                                           gt_last,
                                           logfile != null,
                                           true,
                                           transfer_taxonomy );
            }
            else {
                iterating = true;
                if ( p instanceof NHXParser ) {
                    final NHXParser nhx = ( NHXParser ) p;
                    nhx.setReplaceUnderscores( false );
                    nhx.setIgnoreQuotes( true );
                    nhx.setTaxonomyExtraction( TAXONOMY_EXTRACTION.AGGRESSIVE );
                }
                else if ( p instanceof NexusPhylogeniesParser ) {
                    final NexusPhylogeniesParser nex = ( NexusPhylogeniesParser ) p;
                    nex.setReplaceUnderscores( false );
                    nex.setIgnoreQuotes( true );
                    nex.setTaxonomyExtraction( TAXONOMY_EXTRACTION.AGGRESSIVE );
                }
                else {
                    throw new RuntimeException( "unknown parser type: " + p );
                }
                final IteratingPhylogenyParser ip = ( IteratingPhylogenyParser ) p;
                ip.setSource( gene_trees_file );
                rio = RIO.executeAnalysis( ip,
                                           species_tree_file,
                                           algorithm,
                                           rerooting,
                                           outgroup,
                                           gt_first,
                                           gt_last,
                                           logfile != null,
                                           !use_gene_trees_dir,
                                           transfer_taxonomy );
            }
            if ( !use_gene_trees_dir ) {
                if ( algorithm == ALGORITHM.GSDIR ) {
                    System.out.println( "Taxonomy linking based on           :\t" + rio.getGSDIRtaxCompBase() );
                }
            }
            ///
            ////
            final IntMatrix m;
            if ( iterating ) {
                m = rio.getOrthologTable();
            }
            else {
                m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            }
            final BasicDescriptiveStatistics stats = rio.getDuplicationsStatistics();
            writeTable( orthology_outtable, stats.getN(), m, !use_gene_trees_dir );
            final int ortholog_groups = writeOrtologGroups( orthology_groups_outfile,
                                                            ortholog_group_cutoff,
                                                            stats.getN(),
                                                            m,
                                                            !use_gene_trees_dir,
                                                            false );
            final int ortholog_groups_005 = writeOrtologGroups( null, 0.05, stats.getN(), m, false, true );
            final int ortholog_groups_025 = writeOrtologGroups( null, 0.25, stats.getN(), m, false, true );
            final int ortholog_groups_05 = writeOrtologGroups( null, 0.5, stats.getN(), m, false, true );
            final int ortholog_groups_075 = writeOrtologGroups( null, 0.75, stats.getN(), m, false, true );
            final int ortholog_groups_095 = writeOrtologGroups( null, 0.95, stats.getN(), m, false, true );
            if ( ( algorithm != ALGORITHM.SDIR ) && ( logfile != null ) ) {
                writeLogFile( logfile,
                              rio,
                              species_tree_file,
                              gene_trees_file,
                              orthology_outtable,
                              org.forester.application.rio.PRG_NAME,
                              org.forester.application.rio.PRG_VERSION,
                              org.forester.application.rio.PRG_DATE,
                              ForesterUtil.getForesterLibraryInformation(),
                              !use_gene_trees_dir );
            }
            if ( return_species_tree != null ) {
                writeTree( rio.getSpeciesTree(),
                           return_species_tree,
                           use_gene_trees_dir ? null : "Wrote (stripped) species tree to    :\t" );
            }
            if ( return_min_dup_gene_tree != null && rio.getMinDuplicationsGeneTree() != null ) {
                final int min = ( int ) rio.getDuplicationsStatistics().getMin();
                writeTree( rio.getMinDuplicationsGeneTree(),
                           new File( return_min_dup_gene_tree.toString() + min + ".xml" ),
                           use_gene_trees_dir ? null : "Wrote one min duplication gene tree :\t" );
            }
            if ( return_median_dup_gene_tree != null && rio.getDuplicationsToTreeMap() != null ) {
                final int med = ( int ) rio.getDuplicationsStatistics().median();
                writeTree( rio.getDuplicationsToTreeMap().get( med ),
                           new File( return_median_dup_gene_tree.toString() + med + ".xml" ),
                           use_gene_trees_dir ? null : "Wrote one med duplication gene tree :\t" );
            }
            final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.##" );
            final int min = ( int ) stats.getMin();
            final int max = ( int ) stats.getMax();
            final int median = ( int ) stats.median();
            int min_count = 0;
            int max_count = 0;
            int median_count = 0;
            for( double d : stats.getData() ) {
                if ( ( ( int ) d ) == min ) {
                    ++min_count;
                }
                if ( ( ( int ) d ) == max ) {
                    ++max_count;
                }
                if ( ( ( int ) d ) == median ) {
                    ++median_count;
                }
            }
            final double min_count_percentage = ( 100.0 * min_count ) / stats.getN();
            final double max_count_percentage = ( 100.0 * max_count ) / stats.getN();
            final double median_count_percentage = ( 100.0 * median_count ) / stats.getN();
            if ( use_gene_trees_dir ) {
                String name = gene_trees_file.getName();
                if ( name.indexOf( "." ) > 0 ) {
                    name = name.substring( 0, name.lastIndexOf( "." ) );
                }
                log.print( name );
                log.print( "\t" );
                log.print( Integer.toString( rio.getExtNodesOfAnalyzedGeneTrees() ) );
                log.print( "\t" );
                log.print( Integer.toString( ortholog_groups ) );
                //
                log.print( "\t" );
                log.print( Integer.toString( ortholog_groups_005 ) );
                log.print( "\t" );
                log.print( Integer.toString( ortholog_groups_025 ) );
                log.print( "\t" );
                log.print( Integer.toString( ortholog_groups_05 ) );
                log.print( "\t" );
                log.print( Integer.toString( ortholog_groups_075 ) );
                log.print( "\t" );
                log.print( Integer.toString( ortholog_groups_095 ) );
                //
                log.print( "\t" );
                if ( stats.getN() > 3 ) {
                    log.print( df.format( median ) );
                }
                else {
                    log.print( "" );
                }
                log.print( "\t" );
                log.print( df.format( stats.arithmeticMean() ) );
                log.print( "\t" );
                if ( stats.getN() > 3 ) {
                    log.print( df.format( stats.sampleStandardDeviation() ) );
                }
                else {
                    log.print( "" );
                }
                log.print( "\t" );
                log.print( Integer.toString( min ) );
                log.print( "\t" );
                log.print( Integer.toString( max ) );
                log.print( "\t" );
                log.print( Integer.toString( rio.getRemovedGeneTreeNodes().size() ) );
                log.print( "\t" );
                log.print( Integer.toString( stats.getN() ) );
                log.println();
            }
            else {
                System.out.println( "Gene tree internal nodes            :\t" + rio.getIntNodesOfAnalyzedGeneTrees() );
                System.out.println( "Gene tree external nodes            :\t" + rio.getExtNodesOfAnalyzedGeneTrees() );
                System.out.println( "Mean number of duplications         :\t" + df.format( stats.arithmeticMean() )
                        + "\t" + df.format( ( 100.0 * stats.arithmeticMean() ) / rio.getIntNodesOfAnalyzedGeneTrees() )
                        + "%\t(sd: " + df.format( stats.sampleStandardDeviation() ) + ")" );
                if ( stats.getN() > 3 ) {
                    System.out.println( "Median number of duplications       :\t" + df.format( median ) + "\t"
                            + df.format( ( 100.0 * median ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%" );
                }
                System.out.println( "Minimum duplications                :\t" + min + "\t"
                        + df.format( ( 100.0 * min ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%" );
                System.out.println( "Maximum duplications                :\t" + ( int ) max + "\t"
                        + df.format( ( 100.0 * max ) / rio.getIntNodesOfAnalyzedGeneTrees() ) + "%" );
                System.out.println( "Gene trees with median duplications :\t" + median_count + "\t"
                        + df.format( median_count_percentage ) + "%" );
                System.out.println( "Gene trees with minimum duplications:\t" + min_count + "\t"
                        + df.format( min_count_percentage ) + "%" );
                System.out.println( "Gene trees with maximum duplications:\t" + max_count + "\t"
                        + df.format( max_count_percentage ) + "%" );
                if ( algorithm == ALGORITHM.GSDIR ) {
                    System.out.println( "Removed ext gene tree nodes         :\t"
                            + rio.getRemovedGeneTreeNodes().size() );
                }
            }
        }
        catch ( final RIOException e ) {
            ForesterUtil.fatalError( e.getLocalizedMessage() );
        }
        catch ( final SDIException e ) {
            ForesterUtil.fatalError( e.getLocalizedMessage() );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( e.getLocalizedMessage() );
        }
        catch ( final OutOfMemoryError e ) {
            ForesterUtil.outOfMemoryError( e );
        }
        catch ( final Exception e ) {
            ForesterUtil.unexpectedFatalError( e );
        }
        catch ( final Error e ) {
            ForesterUtil.unexpectedFatalError( e );
        }
    }

    private static final void writeTable( final File table_outfile,
                                          final int gene_trees_analyzed,
                                          final IntMatrix m,
                                          final boolean verbose )
            throws IOException {
        final EasyWriter w = ForesterUtil.createEasyWriter( table_outfile );
        final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.####" );
        df.setDecimalSeparatorAlwaysShown( false );
        df.setRoundingMode( RoundingMode.HALF_UP );
        for( int i = 0; i < m.size(); ++i ) {
            w.print( "\t" );
            w.print( m.getLabel( i ) );
        }
        w.println();
        for( int x = 0; x < m.size(); ++x ) {
            w.print( m.getLabel( x ) );
            for( int y = 0; y < m.size(); ++y ) {
                w.print( "\t" );
                if ( x == y ) {
                    if ( m.get( x, y ) != gene_trees_analyzed ) {
                        ForesterUtil.unexpectedFatalError( "diagonal value is off" );
                    }
                    w.print( "-" );
                }
                else {
                    w.print( df.format( ( ( double ) m.get( x, y ) ) / gene_trees_analyzed ) );
                }
            }
            w.println();
        }
        w.close();
        if ( verbose ) {
            System.out.println( "Wrote table to                      :\t" + table_outfile.getCanonicalPath() );
        }
    }

    private static final int writeOrtologGroups( final File outfile,
                                                 final double cutoff,
                                                 final int gene_trees_analyzed,
                                                 final IntMatrix m,
                                                 final boolean verbose,
                                                 final boolean calc_conly )
            throws IOException {
        List<SortedSet<String>> groups = new ArrayList<SortedSet<String>>();
        BasicDescriptiveStatistics stats = new BasicDescriptiveStatistics();
        int below_075 = 0;
        int below_05 = 0;
        int below_025 = 0;
        for( int x = 1; x < m.size(); ++x ) {
            final String a = m.getLabel( x );
            for( int y = 0; y < x; ++y ) {
                final String b = m.getLabel( y );
                final double s = ( ( double ) m.get( x, y ) ) / gene_trees_analyzed;
                stats.addValue( s );
                if ( s < 0.75 ) {
                    below_075++;
                    if ( s < 0.5 ) {
                        below_05++;
                        if ( s < 0.25 ) {
                            below_025++;
                        }
                    }
                }
                if ( s >= cutoff ) {
                    boolean found = false;
                    for( final SortedSet<String> group : groups ) {
                        if ( group.contains( a ) ) {
                            group.add( b );
                            found = true;
                        }
                        if ( group.contains( b ) ) {
                            group.add( a );
                            found = true;
                        }
                    }
                    if ( !found ) {
                        final SortedSet<String> new_group = new TreeSet<String>();
                        new_group.add( a );
                        new_group.add( b );
                        groups.add( new_group );
                    }
                }
            }
        }
        //Deal with singlets:
        for( int x = 0; x < m.size(); ++x ) {
            final String a = m.getLabel( x );
            boolean found = false;
            for( final SortedSet<String> group : groups ) {
                if ( group.contains( a ) ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                final SortedSet<String> new_group = new TreeSet<String>();
                new_group.add( a );
                groups.add( new_group );
            }
        }
        if ( calc_conly ) {
            return groups.size();
        }
        final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.####" );
        df.setDecimalSeparatorAlwaysShown( false );
        df.setRoundingMode( RoundingMode.HALF_UP );
        final EasyWriter w = ForesterUtil.createEasyWriter( outfile );
        int counter = 1;
        for( final SortedSet<String> group : groups ) {
            w.print( Integer.toString( counter++ ) );
            for( final String s : group ) {
                w.print( "\t" );
                w.print( s );
            }
            w.println();
        }
        w.println();
        w.println( "# Cutoff\t" + df.format( cutoff ) );
        w.println();
        w.println( "# Orthology support statistics:" );
        if ( stats.getN() > 3 ) {
            w.println( "# Median\t" + df.format( stats.median() ) );
        }
        w.println( "# Mean\t" + df.format( stats.arithmeticMean() ) );
        if ( stats.getN() > 3 ) {
            w.println( "# SD\t" + df.format( stats.sampleStandardDeviation() ) );
        }
        w.println( "# Min\t" + df.format( stats.getMin() ) );
        w.println( "# Max\t" + df.format( stats.getMax() ) );
        w.println( "# Total\t" + df.format( stats.getN() ) );
        w.println( "# Below 0.75\t" + below_075 + "\t" + df.format( ( 100.0 * below_075 / stats.getN() ) ) + "%" );
        w.println( "# Below 0.5\t" + below_05 + "\t" + df.format( ( 100.0 * below_05 / stats.getN() ) ) + "%" );
        w.println( "# Below 0.25\t" + below_025 + "\t" + df.format( ( 100.0 * below_025 / stats.getN() ) ) + "%" );
        w.close();
        if ( verbose ) {
            System.out.println( "Number of ortholog groups           :\t" + groups.size() );
            System.out.println( "Wrote orthologs groups table to     :\t" + outfile.getCanonicalPath() );
        }
        return groups.size();
    }

    private static void writeTree( final Phylogeny p, final File f, final String comment ) throws IOException {
        final PhylogenyWriter writer = new PhylogenyWriter();
        writer.toPhyloXML( f, p, 0 );
        if ( comment != null ) {
            System.out.println( comment + f.getCanonicalPath() );
        }
    }

    private static void writeLogFile( final File logfile,
                                      final RIO rio,
                                      final File species_tree_file,
                                      final File gene_trees_file,
                                      final File outtable,
                                      final String prg_name,
                                      final String prg_v,
                                      final String prg_date,
                                      final String f,
                                      final boolean verbose )
            throws IOException {
        final EasyWriter out = ForesterUtil.createEasyWriter( logfile );
        out.println( "# " + prg_name );
        out.println( "# version : " + prg_v );
        out.println( "# date    : " + prg_date );
        out.println( "# based on: " + f );
        out.println( "# ----------------------------------" );
        out.println( "Gene trees                          :\t" + gene_trees_file.getCanonicalPath() );
        out.println( "Species tree                        :\t" + species_tree_file.getCanonicalPath() );
        out.println( "All vs all orthology table          :\t" + outtable.getCanonicalPath() );
        out.flush();
        out.println( rio.getLog().toString() );
        out.close();
        if ( verbose ) {
            System.out.println( "Wrote log to                        :\t" + logfile.getCanonicalPath() );
        }
    }

    private final static Map<String, String> obtainMapping( final File dir, final String prefix, final String suffix )
            throws IOException {
        if ( !dir.exists() ) {
            throw new IOException( "[" + dir + "] does not exist" );
        }
        if ( !dir.isDirectory() ) {
            throw new IOException( "[" + dir + "] is not a directory" );
        }
        final File mapping_files[] = dir.listFiles( new FilenameFilter() {

            @Override
            public boolean accept( final File dir, final String name ) {
                return ( name.endsWith( suffix ) );
            }
        } );
        String my_suffix = suffix;
        boolean done = false;
        do {
            int matches = 0;
            for( File file : mapping_files ) {
                if ( file.getName().equals( my_suffix ) ) {
                    matches++;
                }
            }
            if ( matches == 1) {
                done = true;
            }
            else {
                my_suffix = my_suffix.substring( 0, my_suffix.length() - 1);
            }
        } while (!done );
        
        
        if ( mapping_files.length == 0 ) {
            throw new IOException( "file with prefix \"" + prefix + "\" and suffix \"" + suffix + "\" not found in ["
                    + dir + "] " );
        }
        if ( mapping_files.length > 1 ) {
            throw new IOException( "file with prefix \"" + prefix + "\" and suffix \"" + suffix + "\" not unique in ["
                    + dir + "] " );
        }
        final BasicTable<String> t = BasicTableParser.parse( mapping_files[ 0 ], '\t' );
        return t.getColumnsAsMap( 0, 1 );
    }
}
