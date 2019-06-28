// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2018 Christian M. Zmasek
// Copyright (C) 2018 J. Craig Venter Institute
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

package org.forester.surfacing;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;
import org.forester.ws.seqdb.UniprotData;
import org.forester.ws.seqdb.UniprotRetrieve;

public final class DaioMethods {

    private final static int     VERBOSITY                     = 1;
    private final static int     MAX_IDS_TO_SEARCH_PER_SPECIES = 10;
    private final static String  NA_SYMBOL                     = "-";
    private final static Pattern MULTIPLE_NAME_PATTERN         = Pattern.compile( "(.+)\\s\\[(\\d+)\\]" );
    private final static Pattern ORF_NAME_PATTERN              = Pattern.compile( "orf\\s*(.+)",
                                                                                  Pattern.CASE_INSENSITIVE );
    private final static Pattern ORF_PROTEIN_NAME_PATTERN      = Pattern.compile( "orf\\s*(.+)\\s+protein",
                                                                                  Pattern.CASE_INSENSITIVE );

    public final static String inferTaxonomicInformation( final Phylogeny species_tree, final List<String> all_species )
            throws IOException {
        final List<PhylogenyNode> all_nodes = new ArrayList<>();
        for( final String taxonomy_code : all_species ) {
            final PhylogenyNode node = species_tree.getNodeViaTaxonomyCode( taxonomy_code );
            all_nodes.add( node );
        }
        PhylogenyNode lca = null;
        for( final PhylogenyNode node : all_nodes ) {
            final PhylogenyNode l = PhylogenyMethods.calculateLCA( all_nodes.get( 0 ), node );
            if ( ( lca == null ) || ( l.calculateDepth() < lca.calculateDepth() ) ) {
                lca = l;
            }
        }
        final boolean all;
        if ( lca.getAllExternalDescendants().size() == all_species.size() ) {
            all = true;
        }
        else if ( lca.getAllExternalDescendants().size() > all_species.size() ) {
            all = false;
        }
        else {
            throw new IllegalStateException( "this should never have happened" );
        }
        String taxonomic_string = null;
        if ( lca.isExternal() ) {
            taxonomic_string = "species_specific";
        }
        else {
            if ( lca.isHasNodeData() && lca.getNodeData().isHasTaxonomy() ) {
                if ( !ForesterUtil.isEmpty( lca.getNodeData().getTaxonomy().getSynonyms() )
                        && !ForesterUtil.isEmptyTrimmed( lca.getNodeData().getTaxonomy().getSynonyms().get( 0 ) ) ) {
                    taxonomic_string = lca.getNodeData().getTaxonomy().getSynonyms().get( 0 );
                }
                else if ( !ForesterUtil.isEmptyTrimmed( lca.getNodeData().getTaxonomy().getScientificName() ) ) {
                    taxonomic_string = lca.getNodeData().getTaxonomy().getScientificName();
                }
            }
            if ( ForesterUtil.isEmptyTrimmed( taxonomic_string ) && !ForesterUtil.isEmptyTrimmed( lca.getName() ) ) {
                taxonomic_string = lca.getName();
            }
            if ( ForesterUtil.isEmptyTrimmed( taxonomic_string ) ) {
                throw new IllegalArgumentException( "found LCA in species tree without (1) taxonomy synonym (preferred), (2) taxonomy scientific name, and (3) node name" );
            }
            if ( all ) {
                taxonomic_string = taxonomic_string.toUpperCase();
            }
            else {
                taxonomic_string = taxonomic_string.toLowerCase();
            }
        }
        return taxonomic_string;
    }

    /**
     * Dealing with "ORF" based names in viral genomes.
     *
     */
    public final static String processOrfNames( String name ) {
        final Matcher matcher_orf_p = ORF_PROTEIN_NAME_PATTERN.matcher( name );
        final Matcher matcher_orf = ORF_NAME_PATTERN.matcher( name );
        if ( matcher_orf_p.matches() ) {
            name = "ORF " + matcher_orf_p.group( 1 ) + " protein";
        }
        else if ( matcher_orf.matches() ) {
            name = "ORF " + matcher_orf.group( 1 ) + " protein";
        }
        return name;
    }

    public final static String[] obtainName( final String da,
                                             final SortedMap<String, SortedSet<String>> species_ids_map,
                                             final int max_ids_to_search_per_species,
                                             final int verbosity )
            throws IOException {
        final List<String> ids_per_da = new ArrayList<>();
        final Iterator<Entry<String, SortedSet<String>>> it = species_ids_map.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<String, SortedSet<String>> e = it.next();
            final SortedSet<String> ids = e.getValue();
            if ( ids.size() > 0 ) {
                int counter = 0;
                for( final String id : ids ) {
                    ids_per_da.add( id );
                    if ( ++counter > max_ids_to_search_per_species ) {
                        break;
                    }
                }
            }
        }
        return accessUniprot( ids_per_da, verbosity );
    }

    public final static void outputSpeciesAndIdsPerDAs( final SortedMap<String, SortedMap<String, SortedSet<String>>> da_species_ids_map,
                                                        final Phylogeny species_tree,
                                                        final Writer writer,
                                                        final String separator,
                                                        final String ids_separator,
                                                        final boolean obtain_names_from_db,
                                                        final File input_da_name_file,
                                                        final Writer output_da_name_writer,
                                                        final Writer suffix_da_name_writer )
            throws IOException {
        final SortedMap<String, String> input_da_name_map;
        final Set<String> all_names_lc = new HashSet<>();
        if ( input_da_name_file != null ) {
            final BasicTable<String> input_da_name_table = BasicTableParser.parse( input_da_name_file, '\t' );
            input_da_name_map = input_da_name_table.getColumnsAsMap( 0, 1 );
        }
        else {
            input_da_name_map = null;
        }
        final Iterator<Entry<String, SortedMap<String, SortedSet<String>>>> it = da_species_ids_map.entrySet()
                .iterator();
        if ( obtain_names_from_db ) {
            System.out
                    .println( "[surfacing] > Obtaining names for domain architecures (DA) from Uniprot (slow step)..." );
        }
        int counter = 0;
        final int total = da_species_ids_map.entrySet().size();
        final SortedSet<String> suffix_da_name_set = new TreeSet<>();
        while ( it.hasNext() ) {
            ++counter;
            final Map.Entry<String, SortedMap<String, SortedSet<String>>> e = it.next();
            final String da = e.getKey();
            final SortedMap<String, SortedSet<String>> species_ids_map = e.getValue();
            // Get all species:
            final List<String> all_species = new ArrayList<>();
            final Iterator<Entry<String, SortedSet<String>>> it_for_all_species = species_ids_map.entrySet().iterator();
            String name = null;
            String count = "";
            String out_of = "";
            if ( VERBOSITY > 0 ) {
                System.out.println();
                System.out.println( "[" + counter + "/" + total + "] DA: " + da + ":" );
            }
            if ( ( input_da_name_map != null ) && input_da_name_map.containsKey( da ) ) {
                name = input_da_name_map.get( da );
                if ( VERBOSITY > 0 ) {
                    System.out.println( "  [from file] " + name );
                }
            }
            if ( ( ForesterUtil.isEmptyTrimmed( name ) || name.equals( NA_SYMBOL ) ) && obtain_names_from_db ) {
                final String[] res = obtainName( da, species_ids_map, MAX_IDS_TO_SEARCH_PER_SPECIES, VERBOSITY );
                if ( res != null ) {
                    name = res[ 0 ];
                    count = res[ 1 ];
                    out_of = res[ 2 ];
                }
                else {
                    name = null;
                }
            }
            if ( !ForesterUtil.isEmpty( name ) && !name.equals( NA_SYMBOL ) ) {
                name = processOrfNames( name );
                String name_lc = name.toLowerCase();
                if ( all_names_lc.contains( name_lc ) ) {
                    while ( all_names_lc.contains( name_lc ) ) {
                        final Matcher matcher = MULTIPLE_NAME_PATTERN.matcher( name );
                        if ( matcher.matches() ) {
                            final int d = 1 + Integer.parseInt( matcher.group( 2 ) );
                            name = matcher.group( 1 ) + " [" + d + "]";
                        }
                        else {
                            name = name + " [2]";
                        }
                        name_lc = name.toLowerCase();
                    }
                    if ( VERBOSITY > 0 ) {
                        System.out.println( "   -> " + name );
                    }
                }
                all_names_lc.add( name_lc );
            }
            while ( it_for_all_species.hasNext() ) {
                final Map.Entry<String, SortedSet<String>> e3 = it_for_all_species.next();
                final SortedSet<String> ids = e3.getValue();
                if ( ids.size() > 0 ) {
                    final String taxonomy_code = e3.getKey();
                    all_species.add( taxonomy_code );
                }
            }
            final Iterator<Entry<String, SortedSet<String>>> it2 = species_ids_map.entrySet().iterator();
            output_da_name_writer.write( da );
            output_da_name_writer.write( "\t" );
            if ( !ForesterUtil.isEmpty( name ) ) {
                output_da_name_writer.write( name );
                output_da_name_writer.write( "\t" );
                output_da_name_writer.write( count );
                output_da_name_writer.write( "\t" );
                output_da_name_writer.write( out_of );
            }
            else {
                output_da_name_writer.write( NA_SYMBOL );
                output_da_name_writer.write( "\t" );
                output_da_name_writer.write( "0" );
                output_da_name_writer.write( "\t" );
                output_da_name_writer.write( "0" );
            }
            output_da_name_writer.write( SurfacingConstants.NL );
            output_da_name_writer.flush();
            while ( it2.hasNext() ) {
                final Map.Entry<String, SortedSet<String>> e2 = it2.next();
                final String species = e2.getKey();
                final SortedSet<String> ids = e2.getValue();
                if ( ids.size() > 0 ) {
                    final String taxonomic_suffix = inferTaxonomicInformation( species_tree, all_species );
                    //
                    final StringBuilder suffix_da_name = new StringBuilder();
                    suffix_da_name.append( taxonomic_suffix );
                    suffix_da_name.append( separator );
                    suffix_da_name.append( da );
                    suffix_da_name.append( separator );
                    if ( !ForesterUtil.isEmpty( name ) ) {
                        suffix_da_name.append( name );
                    }
                    else {
                        suffix_da_name.append( NA_SYMBOL );
                    }
                    suffix_da_name_set.add( suffix_da_name.toString() );
                    //
                    writer.write( da );
                    writer.write( separator );
                    if ( !ForesterUtil.isEmpty( name ) ) {
                        writer.write( name );
                    }
                    else {
                        writer.write( NA_SYMBOL );
                    }
                    writer.write( "_" );
                    writer.write( taxonomic_suffix );
                    writer.write( separator );
                    writer.write( species );
                    writer.write( separator );
                    boolean first = true;
                    for( final String id : ids ) {
                        if ( first ) {
                            first = false;
                        }
                        else {
                            writer.write( ids_separator );
                        }
                        writer.write( id );
                    }
                    writer.write( SurfacingConstants.NL );
                }
                writer.flush();
            }
        }
        writer.flush();
        output_da_name_writer.flush();
        final Iterator<String> it3 = suffix_da_name_set.iterator();
        while ( it3.hasNext() ) {
            suffix_da_name_writer.append( it3.next() );
            suffix_da_name_writer.append( SurfacingConstants.NL );
        }
    }

    private static String[] accessUniprot( final List<String> ids, final int verbosity ) throws IOException {
        final UniprotRetrieve ret = new UniprotRetrieve( false );
        final SortedMap<String, UniprotData> m = ret.retrieve( ids );
        final Iterator<Entry<String, UniprotData>> it = m.entrySet().iterator();
        final List<String> names = new ArrayList<>();
        while ( it.hasNext() ) {
            final Map.Entry<String, UniprotData> pair = it.next();
            if ( verbosity > 1 ) {
                System.out.println( "    " + pair.getKey() + " => " + pair.getValue().getProteinNames() );
            }
            String name = pair.getValue().getProteinNames();
            final int index = name.indexOf( " (" );
            if ( index > 0 ) {
                name = name.substring( 0, index );
            }
            final String name_lc = name.toLowerCase();
            if ( !name_lc.startsWith( "uncharacterized" ) ) {
                names.add( name );
            }
        }
        final Optional<Entry<String, Long>> opt = names.stream()
                .collect( Collectors.groupingBy( s -> s, Collectors.counting() ) ).entrySet().stream()
                .max( Comparator.comparing( Entry::getValue ) );
        if ( opt.isPresent() ) {
            if ( verbosity > 0 ) {
                System.out.println( "  [" + opt.get().getValue() + "/" + names.size() + "] " + opt.get().getKey() );
            }
            return new String[] { opt.get().getKey(), String.valueOf( opt.get().getValue() ),
                    String.valueOf( names.size() ) };
        }
        else {
            if ( verbosity > 0 ) {
                System.out.println( "  n/a" );
            }
        }
        return null;
    }
}
