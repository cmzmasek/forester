
package org.forester.surfacing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.application.surfacing;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.protein.Domain;
import org.forester.protein.Protein;
import org.forester.species.BasicSpecies;
import org.forester.species.Species;
import org.forester.surfacing.SurfacingUtil.DomainComparator;
import org.forester.util.ForesterUtil;

public final class MinimalDomainomeCalculator {

    public final static void calc( final boolean use_domain_architectures,
                                   final Phylogeny tre,
                                   final int target_level,
                                   final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                   final String separator,
                                   final double ie_cutoff,
                                   final String outfile_base,
                                   final boolean write_protein_files )
            throws IOException {
        final SortedMap<String, SortedSet<String>> species_to_features_map = new TreeMap<String, SortedSet<String>>();
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
        out.write( "SPECIES\tCOMMON NAME\tCODE\tRANK\t#EXT NODES\tEXT NODE CODES\t#" + x + "\t" + x + "" );
        out.write( ForesterUtil.LINE_SEPARATOR );
        ///////////
        //////////
        SortedMap<String, List<Protein>> protein_lists_per_quasi_species = null;
        if ( target_level >= 1 ) {
            protein_lists_per_quasi_species = makeProteinListsPerQuasiSpecies( tre,
                                                                               target_level,
                                                                               protein_lists_per_species );
           
        }
        /////////
        ///////////
        for( final PhylogenyNodeIterator iter = tre.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final int node_level = PhylogenyMethods.calculateLevel( node );
            final String species_name = node.getNodeData().isHasTaxonomy()
                    ? node.getNodeData().getTaxonomy().getScientificName() : node.getName();
            final String common = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy().getCommonName()
                    : "";
            final String tcode = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy().getTaxonomyCode()
                    : "";
            final String rank = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy().getRank() : "";
            final List<PhylogenyNode> external_descs = node.getAllExternalDescendants();
            if ( ( target_level < 1 ) || ( node_level >= target_level ) ) {
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
                if ( node.isInternal() ) {
                    out.write( "\t" + external_descs.size() + "\t" );
                }
                else {
                    out.write( "\t\t" );
                }
            }
            final List<Set<String>> features_per_genome_list = new ArrayList<Set<String>>();
            boolean first = true;
            if ( target_level >= 1 ) {
                ////////////
                ////////////
                if ( node_level >= target_level ) {
                    final List<PhylogenyNode> given_level_descs = PhylogenyMethods
                            .getAllDescendantsOfGivenLevel( node, target_level );
                    for( final PhylogenyNode given_level_desc : given_level_descs ) {
                        final String spec_name = given_level_desc.getNodeData().isHasTaxonomy()
                                ? given_level_desc.getNodeData().getTaxonomy().getScientificName()
                                : given_level_desc.getName();
                        if ( node.isInternal() ) {
                            if ( first ) {
                                first = false;
                            }
                            else {
                                out.write( ", " );
                            }
                            out.write( "sp_n=" + spec_name );
                        }
                        final List<Protein> proteins_per_species = protein_lists_per_quasi_species.get( spec_name );
                        if ( proteins_per_species != null ) {
                            final SortedSet<String> features_per_genome = new TreeSet<String>();
                            for( final Protein protein : proteins_per_species ) {
                                if ( use_domain_architectures ) {
                                    final String da = protein.toDomainArchitectureString( separator, ie_cutoff );
                                    features_per_genome.add( da );
                                }
                                else {
                                    List<Domain> domains = protein.getProteinDomains();
                                    for( final Domain domain : domains ) {
                                        if ( ( ie_cutoff <= -1 ) || ( domain.getPerDomainEvalue() <= ie_cutoff ) ) {
                                            features_per_genome.add( domain.getDomainId() );
                                        }
                                    }
                                }
                            }
                            System.out.println( ">>>>>>>>>>>>>> features_per_genome.size()=" + features_per_genome.size() );
                            if ( features_per_genome.size() > 0 ) {
                                features_per_genome_list.add( features_per_genome );
                            }
                            else {
                                System.out.println( "error!" );
                                System.exit( -1 );
                            }
                        }
                        else {
                            System.out.println( "error!" );
                            System.exit( -1 );
                        }
                    }
                }
                ///////////
                ///////////
            }
            else {
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
                    final List<Protein> proteins_per_species = protein_lists_per_species
                            .get( new BasicSpecies( code ) );
                    if ( proteins_per_species != null ) {
                        final SortedSet<String> features_per_genome = new TreeSet<String>();
                        for( final Protein protein : proteins_per_species ) {
                            if ( use_domain_architectures ) {
                                final String da = protein.toDomainArchitectureString( separator, ie_cutoff );
                                features_per_genome.add( da );
                            }
                            else {
                                List<Domain> domains = protein.getProteinDomains();
                                for( final Domain domain : domains ) {
                                    if ( ( ie_cutoff <= -1 ) || ( domain.getPerDomainEvalue() <= ie_cutoff ) ) {
                                        features_per_genome.add( domain.getDomainId() );
                                    }
                                }
                            }
                        }
                        if ( features_per_genome.size() > 0 ) {
                            features_per_genome_list.add( features_per_genome );
                        }
                    }
                } // for( final PhylogenyNode external_desc : external_descs )
            } // else
            if ( features_per_genome_list.size() > 0 ) {
                SortedSet<String> intersection = calcIntersection( features_per_genome_list );
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
                species_to_features_map.put( species_name, intersection );
            }
        }
        final SortedSet<String> all_species_names = new TreeSet<String>();
        final SortedSet<String> all_features = new TreeSet<String>();
        for( final Entry<String, SortedSet<String>> e : species_to_features_map.entrySet() ) {
            all_species_names.add( e.getKey() );
            for( final String f : e.getValue() ) {
                all_features.add( f );
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
        for( final String das : all_features ) {
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
                if ( species_to_features_map.get( species_name ).contains( das ) ) {
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
        if ( write_protein_files ) {
            final String protdirname;
            final String a;
            final String b;
            if ( use_domain_architectures ) {
                a = "_DA";
                b = "domain architectures (DAs)";
                protdirname = "_DAS";
            }
            else {
                a = "_domain";
                b = "domains";
                protdirname = "_DOMAINS";
            }
            final File prot_dir = new File( outfile_base + protdirname );
            final boolean success = prot_dir.mkdir();
            if ( !success ) {
                throw new IOException( "failed to create dir " + prot_dir );
            }
            int total = 0;
            final String dir = outfile_base + protdirname + "/";
            for( final String feat : all_features ) {
                final File extract_outfile = new File( dir + feat + a + surfacing.SEQ_EXTRACT_SUFFIX );
                SurfacingUtil.checkForOutputFileWriteability( extract_outfile );
                final Writer proteins_file_writer = new BufferedWriter( new FileWriter( extract_outfile ) );
                final int counter = extractProteinFeatures( use_domain_architectures,
                                                            protein_lists_per_species,
                                                            feat,
                                                            proteins_file_writer,
                                                            ie_cutoff,
                                                            separator );
                if ( counter < 1 ) {
                    ForesterUtil.printWarningMessage( "surfacing", feat + " not present (in " + b + " extraction)" );
                }
                total += counter;
                proteins_file_writer.close();
            }
            ForesterUtil.programMessage( "surfacing",
                                         "Wrote " + total + " individual " + b + " from a total of "
                                                 + all_features.size() + " into: " + dir );
        }
    }

    private final static SortedMap<String, List<Protein>> makeProteinListsPerQuasiSpecies( final Phylogeny tre,
                                                                                           final int level,
                                                                                           final SortedMap<Species, List<Protein>> protein_lists_per_species ) {
        final SortedMap<String, List<Protein>> protein_lists_per_quasi_species = new TreeMap<String, List<Protein>>();
        System.out.println( "---------------------------------" );
        System.out.println( "level=" + level );
        for( final PhylogenyNodeIterator iter = tre.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final int node_level = PhylogenyMethods.calculateLevel( node );
            if ( node_level == level ) {
                System.out.println( "level=" + level );
                final List<PhylogenyNode> external_descs = node.getAllExternalDescendants();
                final List<Protein> protein_list_per_quasi_species = new ArrayList<Protein>();
                for( final PhylogenyNode external_desc : external_descs ) {
                    final String code = external_desc.getNodeData().getTaxonomy().getTaxonomyCode();
                    final List<Protein> proteins_per_species = protein_lists_per_species
                            .get( new BasicSpecies( code ) );
                    //System.out.println( code );
                    for( Protein protein : proteins_per_species ) {
                        protein_list_per_quasi_species.add( protein );
                    }
                }
                final String species_name = node.getNodeData().isHasTaxonomy()
                        ? node.getNodeData().getTaxonomy().getScientificName() : node.getName();
                System.out.println( "species_name=" + species_name );
                protein_lists_per_quasi_species.put( species_name, protein_list_per_quasi_species );
                System.out.println( ">>>>" + protein_list_per_quasi_species.size() );
            }
        }
      
        return protein_lists_per_quasi_species;
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

    private final static int extractProteinFeatures( final boolean use_domain_architectures,
                                                     final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                                     final String domain_id,
                                                     final Writer out,
                                                     final double ie_cutoff,
                                                     final String domain_separator )
            throws IOException {
        int counter = 0;
        final String separator_for_output = "\t";
        for( final Species species : protein_lists_per_species.keySet() ) {
            final List<Protein> proteins_per_species = protein_lists_per_species.get( species );
            for( final Protein protein : proteins_per_species ) {
                if ( use_domain_architectures ) {
                    if ( domain_id.equals( protein.toDomainArchitectureString( domain_separator, ie_cutoff ) ) ) {
                        int from = Integer.MAX_VALUE;
                        int to = -1;
                        for( final Domain d : protein.getProteinDomains() ) {
                            if ( ( ie_cutoff <= -1 ) || ( d.getPerDomainEvalue() <= ie_cutoff ) ) {
                                if ( d.getFrom() < from ) {
                                    from = d.getFrom();
                                }
                                if ( d.getTo() > to ) {
                                    to = d.getTo();
                                }
                            }
                        }
                        out.write( protein.getSpecies().getSpeciesId() );
                        out.write( separator_for_output );
                        out.write( protein.getProteinId().getId() );
                        out.write( separator_for_output );
                        out.write( domain_id );
                        out.write( separator_for_output );
                        out.write( "/" );
                        out.write( from + "-" + to );
                        out.write( "/" );
                        out.write( SurfacingConstants.NL );
                        ++counter;
                    }
                }
                else {
                    final List<Domain> domains = protein.getProteinDomains( domain_id );
                    if ( domains.size() > 0 ) {
                        out.write( protein.getSpecies().getSpeciesId() );
                        out.write( separator_for_output );
                        out.write( protein.getProteinId().getId() );
                        out.write( separator_for_output );
                        out.write( domain_id );
                        out.write( separator_for_output );
                        for( final Domain domain : domains ) {
                            if ( ( ie_cutoff < 0 ) || ( domain.getPerDomainEvalue() <= ie_cutoff ) ) {
                                out.write( "/" );
                                out.write( domain.getFrom() + "-" + domain.getTo() );
                            }
                        }
                        out.write( "/" );
                        out.write( separator_for_output );
                        final List<Domain> domain_list = new ArrayList<Domain>();
                        for( final Domain domain : protein.getProteinDomains() ) {
                            if ( ( ie_cutoff < 0 ) || ( domain.getPerDomainEvalue() <= ie_cutoff ) ) {
                                domain_list.add( domain );
                            }
                        }
                        final Domain domain_ary[] = new Domain[ domain_list.size() ];
                        for( int i = 0; i < domain_list.size(); ++i ) {
                            domain_ary[ i ] = domain_list.get( i );
                        }
                        Arrays.sort( domain_ary, new DomainComparator( true ) );
                        out.write( "{" );
                        boolean first = true;
                        for( final Domain domain : domain_ary ) {
                            if ( first ) {
                                first = false;
                            }
                            else {
                                out.write( "," );
                            }
                            out.write( domain.getDomainId().toString() );
                            out.write( ":" + domain.getFrom() + "-" + domain.getTo() );
                            out.write( ":" + domain.getPerDomainEvalue() );
                        }
                        out.write( "}" );
                        if ( !( ForesterUtil.isEmpty( protein.getDescription() )
                                || protein.getDescription().equals( SurfacingConstants.NONE ) ) ) {
                            out.write( protein.getDescription() );
                        }
                        out.write( separator_for_output );
                        if ( !( ForesterUtil.isEmpty( protein.getAccession() )
                                || protein.getAccession().equals( SurfacingConstants.NONE ) ) ) {
                            out.write( protein.getAccession() );
                        }
                        out.write( SurfacingConstants.NL );
                        ++counter;
                    }
                }
            }
        }
        out.flush();
        return counter;
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
