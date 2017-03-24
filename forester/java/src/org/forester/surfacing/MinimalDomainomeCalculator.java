
package org.forester.surfacing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.application.surfacing;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.protein.Domain;
import org.forester.protein.Protein;
import org.forester.species.BasicSpecies;
import org.forester.species.Species;
import org.forester.util.ForesterUtil;

public final class MinimalDomainomeCalculator {

    static final public void calcDomainome( final Phylogeny tre,
                                            final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                            final double ie_cutoff ) {
        if ( protein_lists_per_species == null || tre == null ) {
            throw new IllegalArgumentException( "argument is null" );
        }
        if ( protein_lists_per_species.size() < 2 ) {
            throw new IllegalArgumentException( "not enough genomes" );
        }
        for( final PhylogenyNodeIterator iter = tre.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isInternal() ) {
                System.out.println();
                if ( node.getNodeData().isHasTaxonomy() ) {
                    System.out.println( node.getNodeData().getTaxonomy().getScientificName() + ":" );
                }
                else {
                    System.out.println( node.getName() + ":" );
                }
                final List<PhylogenyNode> external_descs = node.getAllExternalDescendants();
                final List<Set<String>> domains_per_genome_list = new ArrayList<Set<String>>();
                for( final PhylogenyNode external_desc : external_descs ) {
                    final String code = external_desc.getNodeData().getTaxonomy().getTaxonomyCode();
                    System.out.print( code + " " );
                    final List<Protein> proteins_per_species = protein_lists_per_species
                            .get( new BasicSpecies( code ) );
                    if ( proteins_per_species != null ) {
                        final SortedSet<String> domains_per_genome = new TreeSet<String>();
                        for( final Protein protein : proteins_per_species ) {
                            List<Domain> domains = protein.getProteinDomains();
                            for( final Domain domain : domains ) {
                                if ( ( domain.getPerDomainEvalue() <= ie_cutoff ) || ( ie_cutoff <= -1 ) ) {
                                    domains_per_genome.add( domain.getDomainId() );
                                }
                            }
                        }
                        if ( domains_per_genome.size() > 0 ) {
                            domains_per_genome_list.add( domains_per_genome );
                        }
                    }
                }
                System.out.println();
                if ( domains_per_genome_list.size() > 0 ) {
                    Set<String> intersection = calcIntersection( domains_per_genome_list );
                    System.out.println( intersection );
                }
            }
        }
    }

    static final public void calcOme( final boolean use_domain_architectures,
                                      final Phylogeny tre,
                                      final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                      final String separator,
                                      final double ie_cutoff,
                                      final String outfile_base )
            throws IOException {
        final SortedMap<String, SortedSet<String>> species_to_das_map = new TreeMap<String, SortedSet<String>>();
        if ( protein_lists_per_species == null || tre == null ) {
            throw new IllegalArgumentException( "argument is null" );
        }
        if ( protein_lists_per_species.size() < 2 ) {
            throw new IllegalArgumentException( "not enough genomes" );
        }
        final String x;
        if ( use_domain_architectures ) {
            x = "DA";
        }
        else {
            x = "domain";
        }
        final File outfile = new File( outfile_base + "_minimal_" + x + "ome.tsv" );
        final File outfile_table = new File( outfile_base + "_minimal_" + x + "ome_matrix.tsv" );
        SurfacingUtil.checkForOutputFileWriteability( outfile );
        SurfacingUtil.checkForOutputFileWriteability( outfile_table );
        final BufferedWriter out = new BufferedWriter( new FileWriter( outfile ) );
        final BufferedWriter out_table = new BufferedWriter( new FileWriter( outfile_table ) );
        out.write( "SPECIES\tCOMMON NAME\tCODE\tRANK\t#EXT NODES\tEXT NODE CODES\t#DA\tDA" );
        out.write( ForesterUtil.LINE_SEPARATOR );
        for( final PhylogenyNodeIterator iter = tre.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final String species_name = node.getNodeData().isHasTaxonomy()
                    ? node.getNodeData().getTaxonomy().getScientificName() : node.getName();
            final String common = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy().getCommonName()
                    : "";
            final String tcode = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy().getTaxonomyCode()
                    : "";
            final String rank = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy().getRank() : "";
            out.write( species_name );
            if ( !ForesterUtil.isEmpty( common ) ) {
                out.write( "\t" + common );
            }
            else {
                out.write( "\t" );
            }
            if ( !ForesterUtil.isEmpty( tcode ) ) {
                out.write( "\t" + tcode );
            }
            else {
                out.write( "\t" );
            }
            if ( !ForesterUtil.isEmpty( rank ) ) {
                out.write( "\t" + rank );
            }
            else {
                out.write( "\t" );
            }
            final List<PhylogenyNode> external_descs = node.getAllExternalDescendants();
            if ( node.isInternal() ) {
                out.write( "\t" + external_descs.size() + "\t" );
            }
            else {
                out.write( "\t\t" );
            }
            final List<Set<String>> das_per_genome_list = new ArrayList<Set<String>>();
            boolean first = true;
            for( final PhylogenyNode external_desc : external_descs ) {
                final String code = external_desc.getNodeData().getTaxonomy().getTaxonomyCode();
                if ( node.isInternal() ) {
                    if ( first ) {
                        first = false;
                    }
                    else {
                        out.write( ", " );
                    }
                    out.write( code );
                }
                final List<Protein> proteins_per_species = protein_lists_per_species.get( new BasicSpecies( code ) );
                if ( proteins_per_species != null ) {
                    final SortedSet<String> das_per_genome = new TreeSet<String>();
                    for( final Protein protein : proteins_per_species ) {
                        if ( use_domain_architectures ) {
                            final String da = protein.toDomainArchitectureString( separator, ie_cutoff );
                            das_per_genome.add( da );
                        }
                        else {
                            List<Domain> domains = protein.getProteinDomains();
                            for( final Domain domain : domains ) {
                                if ( ( ie_cutoff <= -1 ) || ( domain.getPerDomainEvalue() <= ie_cutoff ) ) {
                                    das_per_genome.add( domain.getDomainId() );
                                }
                            }
                        }
                    }
                    if ( das_per_genome.size() > 0 ) {
                        das_per_genome_list.add( das_per_genome );
                    }
                }
            }
            if ( das_per_genome_list.size() > 0 ) {
                SortedSet<String> intersection = calcIntersection( das_per_genome_list );
                out.write( "\t" + intersection.size() + "\t" );
                first = true;
                for( final String s : intersection ) {
                    if ( first ) {
                        first = false;
                    }
                    else {
                        out.write( ", " );
                    }
                    out.write( s );
                }
                out.write( ForesterUtil.LINE_SEPARATOR );
                species_to_das_map.put( species_name, intersection );
            }
        }
        final SortedSet<String> all_species_names = new TreeSet<String>();
        final SortedSet<String> all_das = new TreeSet<String>();
        for( final Entry<String, SortedSet<String>> e : species_to_das_map.entrySet() ) {
            all_species_names.add( e.getKey() );
            for( final String das : e.getValue() ) {
                all_das.add( das );
            }
        }
        out_table.write( '\t' );
        boolean first = true;
        for( final String species_name : all_species_names ) {
            if ( first ) {
                first = false;
            }
            else {
                out_table.write( '\t' );
            }
            out_table.write( species_name );
        }
        out_table.write( ForesterUtil.LINE_SEPARATOR );
        for( final String das : all_das ) {
            out_table.write( das );
            out_table.write( '\t' );
            first = true;
            for( final String species_name : all_species_names ) {
                if ( first ) {
                    first = false;
                }
                else {
                    out_table.write( '\t' );
                }
                if ( species_to_das_map.get( species_name ).contains( das ) ) {
                    out_table.write( '1' );
                }
                else {
                    out_table.write( '0' );
                }
            }
            out_table.write( ForesterUtil.LINE_SEPARATOR );
        }
        out.flush();
        out.close();
        out_table.flush();
        out_table.close();
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote minimal DAome data to           : " + outfile );
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote minimal DAome data to (as table): " + outfile_table );
    }

    static final public void calcDAome( final Phylogeny tre,
                                        final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                        final String separator,
                                        final double ie_cutoff,
                                        final String outfile_base )
            throws IOException {
        final SortedMap<String, SortedSet<String>> species_to_das_map = new TreeMap<String, SortedSet<String>>();
        if ( protein_lists_per_species == null || tre == null ) {
            throw new IllegalArgumentException( "argument is null" );
        }
        if ( protein_lists_per_species.size() < 2 ) {
            throw new IllegalArgumentException( "not enough genomes" );
        }
        final File outfile = new File( outfile_base + "_minimal_daome.txt" );
        final File outfile_table = new File( outfile_base + "_minimal_daome.tsv" );
        SurfacingUtil.checkForOutputFileWriteability( outfile );
        SurfacingUtil.checkForOutputFileWriteability( outfile_table );
        final BufferedWriter out = new BufferedWriter( new FileWriter( outfile ) );
        final BufferedWriter out_table = new BufferedWriter( new FileWriter( outfile_table ) );
        out.write( "SPECIES\tCOMMON NAME\tCODE\tRANK\t#EXT NODES\tEXT NODE CODES\t#DA\tDA" );
        out.write( ForesterUtil.LINE_SEPARATOR );
        for( final PhylogenyNodeIterator iter = tre.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final String species_name = node.getNodeData().isHasTaxonomy()
                    ? node.getNodeData().getTaxonomy().getScientificName() : node.getName();
            final String common = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy().getCommonName()
                    : "";
            final String tcode = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy().getTaxonomyCode()
                    : "";
            final String rank = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy().getRank() : "";
            out.write( species_name );
            if ( !ForesterUtil.isEmpty( common ) ) {
                out.write( "\t" + common );
            }
            else {
                out.write( "\t" );
            }
            if ( !ForesterUtil.isEmpty( tcode ) ) {
                out.write( "\t" + tcode );
            }
            else {
                out.write( "\t" );
            }
            if ( !ForesterUtil.isEmpty( rank ) ) {
                out.write( "\t" + rank );
            }
            else {
                out.write( "\t" );
            }
            final List<PhylogenyNode> external_descs = node.getAllExternalDescendants();
            if ( node.isInternal() ) {
                out.write( "\t" + external_descs.size() + "\t" );
            }
            else {
                out.write( "\t\t" );
            }
            final List<Set<String>> das_per_genome_list = new ArrayList<Set<String>>();
            boolean first = true;
            for( final PhylogenyNode external_desc : external_descs ) {
                final String code = external_desc.getNodeData().getTaxonomy().getTaxonomyCode();
                if ( node.isInternal() ) {
                    if ( first ) {
                        first = false;
                    }
                    else {
                        out.write( ", " );
                    }
                    out.write( code );
                }
                final List<Protein> proteins_per_species = protein_lists_per_species.get( new BasicSpecies( code ) );
                if ( proteins_per_species != null ) {
                    final SortedSet<String> das_per_genome = new TreeSet<String>();
                    for( final Protein protein : proteins_per_species ) {
                        final String da = protein.toDomainArchitectureString( separator, ie_cutoff );
                        das_per_genome.add( da );
                    }
                    if ( das_per_genome.size() > 0 ) {
                        das_per_genome_list.add( das_per_genome );
                    }
                }
            }
            if ( das_per_genome_list.size() > 0 ) {
                SortedSet<String> intersection = calcIntersection( das_per_genome_list );
                out.write( "\t" + intersection.size() + "\t" );
                first = true;
                for( final String s : intersection ) {
                    if ( first ) {
                        first = false;
                    }
                    else {
                        out.write( ", " );
                    }
                    out.write( s );
                }
                out.write( ForesterUtil.LINE_SEPARATOR );
                species_to_das_map.put( species_name, intersection );
            }
        }
        final SortedSet<String> all_species_names = new TreeSet<String>();
        final SortedSet<String> all_das = new TreeSet<String>();
        for( final Entry<String, SortedSet<String>> e : species_to_das_map.entrySet() ) {
            all_species_names.add( e.getKey() );
            for( final String das : e.getValue() ) {
                all_das.add( das );
            }
        }
        out_table.write( '\t' );
        boolean first = true;
        for( final String species_name : all_species_names ) {
            if ( first ) {
                first = false;
            }
            else {
                out_table.write( '\t' );
            }
            out_table.write( species_name );
        }
        out_table.write( ForesterUtil.LINE_SEPARATOR );
        for( final String das : all_das ) {
            out_table.write( das );
            out_table.write( '\t' );
            first = true;
            for( final String species_name : all_species_names ) {
                if ( first ) {
                    first = false;
                }
                else {
                    out_table.write( '\t' );
                }
                if ( species_to_das_map.get( species_name ).contains( das ) ) {
                    out_table.write( '1' );
                }
                else {
                    out_table.write( '0' );
                }
            }
            out_table.write( ForesterUtil.LINE_SEPARATOR );
        }
        out.flush();
        out.close();
        out_table.flush();
        out_table.close();
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote minimal DAome data to           : " + outfile );
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote minimal DAome data to (as table): " + outfile_table );
    }

    private final static SortedSet<String> calcIntersection( final List<Set<String>> features_per_genome_list ) {
        final Set<String> first = features_per_genome_list.get( 0 );
        final SortedSet<String> my_first = new TreeSet<String>();
        for( final String s : first ) {
            my_first.add( s );
        }
        for( int i = 1; i < features_per_genome_list.size(); ++i ) {
            my_first.retainAll( features_per_genome_list.get( i ) );
        }
        return my_first;
    }

    public static void main( final String[] args ) {
        Set<String> a = new HashSet<String>();
        Set<String> b = new HashSet<String>();
        Set<String> c = new HashSet<String>();
        Set<String> d = new HashSet<String>();
        a.add( "x" );
        a.add( "b" );
        a.add( "c" );
        b.add( "a" );
        b.add( "b" );
        b.add( "c" );
        c.add( "a" );
        c.add( "b" );
        c.add( "c" );
        c.add( "c" );
        c.add( "f" );
        d.add( "a" );
        d.add( "c" );
        d.add( "d" );
        List<Set<String>> domains_per_genome_list = new ArrayList<Set<String>>();
        domains_per_genome_list.add( a );
        domains_per_genome_list.add( b );
        domains_per_genome_list.add( c );
        domains_per_genome_list.add( d );
        Set<String> x = calcIntersection( domains_per_genome_list );
        System.out.println( x );
    }
}
