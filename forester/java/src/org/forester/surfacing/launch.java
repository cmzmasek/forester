// $Id:
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2018 Christian M. Zmasek
// Copyright (C) 2018 J. Craig Venter Institute
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.surfacing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.go.GoId;
import org.forester.go.GoNameSpace;
import org.forester.go.GoTerm;
import org.forester.go.GoUtils;
import org.forester.go.OBOparser;
import org.forester.go.PfamToGoMapping;
import org.forester.go.PfamToGoParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser.INDIVIDUAL_SCORE_CUTOFF;
import org.forester.phylogeny.Phylogeny;
import org.forester.protein.BinaryDomainCombination;
import org.forester.protein.Domain;
import org.forester.protein.Protein;
import org.forester.species.BasicSpecies;
import org.forester.species.Species;
import org.forester.surfacing.DomainSimilarity.DomainSimilarityScoring;
import org.forester.surfacing.DomainSimilarity.PRINT_OPTION;
import org.forester.surfacing.DomainSimilarityCalculator.Detailedness;
import org.forester.surfacing.GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public final class launch {

    private static final int                                        MINIMAL_NUMBER_OF_SIMILARITIES_FOR_SPLITTING                       = 1000;
    private static final String                                     DOMAIN_COMBINITONS_OUTPUT_OPTION_FOR_GRAPH_ANALYSIS                = "graph_analysis_out";
    private static final String                                     DOMAIN_COMBINITONS_COUNTS_OUTPUT_OPTION                            = "dcc";
    private static final String                                     HELP_OPTION_1                                                      = "help";
    private static final String                                     HELP_OPTION_2                                                      = "h";
    private static final String                                     OUTPUT_DIR_OPTION                                                  = "out_dir";
    private static final String                                     SCORING_OPTION                                                     = "scoring";
    private static final DomainSimilarityScoring                    SCORING_DEFAULT                                                    = DomainSimilarity.DomainSimilarityScoring.COMBINATIONS;
    private static final String                                     SCORING_DOMAIN_COUNT_BASED                                         = "domains";
    private static final String                                     SCORING_PROTEIN_COUNT_BASED                                        = "proteins";
    private static final String                                     SCORING_COMBINATION_BASED                                          = "combinations";
    private static final String                                     DETAILEDNESS_OPTION                                                = "detail";
    private static final Detailedness                               DETAILEDNESS_DEFAULT                                               = DomainSimilarityCalculator.Detailedness.PUNCTILIOUS;
    private static final String                                     SPECIES_MATRIX_OPTION                                              = "smatrix";
    private static final String                                     DETAILEDNESS_BASIC                                                 = "basic";
    private static final String                                     DETAILEDNESS_LIST_IDS                                              = "list_ids";
    private static final String                                     DETAILEDNESS_PUNCTILIOUS                                           = "punctilious";
    private static final String                                     DOMAIN_SIMILARITY_SORT_OPTION                                      = "sort";
    private static final DomainSimilarity.DomainSimilaritySortField DOMAIN_SORT_FILD_DEFAULT                                           = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
    private static final String                                     DOMAIN_SIMILARITY_SORT_MIN                                         = "min";
    private static final String                                     DOMAIN_SIMILARITY_SORT_MAX                                         = "max";
    private static final String                                     DOMAIN_SIMILARITY_SORT_SD                                          = "sd";
    private static final String                                     DOMAIN_SIMILARITY_SORT_MEAN                                        = "mean";
    private static final String                                     DOMAIN_SIMILARITY_SORT_DIFF                                        = "diff";
    private static final String                                     DOMAIN_SIMILARITY_SORT_COUNTS_DIFF                                 = "count_diff";
    private static final String                                     DOMAIN_SIMILARITY_SORT_ABS_COUNTS_DIFF                             = "abs_count_diff";
    private static final String                                     DOMAIN_SIMILARITY_SORT_SPECIES_COUNT                               = "species";
    private static final String                                     DOMAIN_SIMILARITY_SORT_ALPHA                                       = "alpha";
    private static final String                                     DOMAIN_SIMILARITY_SORT_BY_SPECIES_COUNT_FIRST_OPTION               = "species_first";
    private static final String                                     DOMAIN_COUNT_SORT_OPTION                                           = "dc_sort";
    private static final GenomeWideCombinableDomainsSortOrder       DOMAINS_SORT_ORDER_DEFAULT                                         = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.ALPHABETICAL_KEY_ID;
    private static final String                                     DOMAIN_COUNT_SORT_ALPHA                                            = "alpha";
    private static final String                                     DOMAIN_COUNT_SORT_KEY_DOMAIN_COUNT                                 = "dom";
    private static final String                                     DOMAIN_COUNT_SORT_KEY_DOMAIN_PROTEINS_COUNT                        = "prot";
    private static final String                                     DOMAIN_COUNT_SORT_COMBINATIONS_COUNT                               = "comb";
    private static final String                                     CUTOFF_SCORE_FILE_OPTION                                           = "cos";
    private static final String                                     NOT_IGNORE_DUFS_OPTION                                             = "dufs";
    private static final String                                     MAX_FS_E_VALUE_OPTION                                              = "fs_e";
    private static final String                                     MAX_I_E_VALUE_OPTION                                               = "ie";
    private static final String                                     MIN_REL_ENV_LENGTH_RATIO_OPTION                                    = "mrel";
    private static final String                                     NO_ENGULFING_OVERLAP_OPTION                                        = "no_eo";
    private static final String                                     IGNORE_COMBINATION_WITH_SAME_OPTION                                = "ignore_self_comb";
    private static final String                                     PERFORM_DC_REGAIN_PROTEINS_STATS_OPTION                            = "dc_regain_stats";
    private static final String                                     DA_ANALYSIS_OPTION                                                 = "da_analyis";
    private static final String                                     USE_LAST_IN_FITCH_OPTION                                           = "last";
    private static final String                                     PAIRWISE_DOMAIN_COMPARISONS_OPTION                                 = "pwc";
    private static final String                                     OUTPUT_FILE_OPTION                                                 = "o";
    private static final String                                     PFAM_TO_GO_FILE_USE_OPTION                                         = "p2g";
    private static final String                                     GO_OBO_FILE_USE_OPTION                                             = "obo";
    private static final String                                     GO_NAMESPACE_LIMIT_OPTION                                          = "go_namespace";
    private static final String                                     GO_NAMESPACE_LIMIT_OPTION_MOLECULAR_FUNCTION                       = "molecular_function";
    private static final String                                     GO_NAMESPACE_LIMIT_OPTION_BIOLOGICAL_PROCESS                       = "biological_process";
    private static final String                                     GO_NAMESPACE_LIMIT_OPTION_CELLULAR_COMPONENT                       = "cellular_component";
    private static final String                                     SECONDARY_FEATURES_PARSIMONY_MAP_FILE                              = "secondary";
    private static final String                                     DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_TAB_DELIMITED                = "simple_tab";
    private static final String                                     DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_HTML                         = "simple_html";
    private static final String                                     DOMAIN_SIMILARITY_PRINT_OPTION_DETAILED_HTML                       = "detailed_html";
    private static final String                                     DOMAIN_SIMILARITY_PRINT_OPTION                                     = "ds_output";
    private static final PRINT_OPTION                               DOMAIN_SIMILARITY_PRINT_OPTION_DEFAULT                             = DomainSimilarity.PRINT_OPTION.HTML;
    private static final String                                     IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_OPTION          = "ignore_singlet_domains";
    private static final String                                     IGNORE_VIRAL_IDS                                                   = "ignore_viral_ids";
    private static final boolean                                    IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_DEFAULT         = false;
    private static final String                                     IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION                      = "ignore_species_specific_domains";
    private static final boolean                                    IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION_DEFAULT              = false;
    private static final String                                     MATRIX_MEAN_SCORE_BASED_GENOME_DISTANCE_SUFFIX                     = "_mean_score.pwd";
    private static final String                                     MATRIX_SHARED_DOMAINS_BASED_GENOME_DISTANCE_SUFFIX                 = "_domains.pwd";
    private static final String                                     MATRIX_SHARED_BIN_COMBINATIONS_BASED_GENOME_DISTANCE_SUFFIX        = "_bin_combinations.pwd";
    private static final String                                     NJ_TREE_MEAN_SCORE_BASED_GENOME_DISTANCE_SUFFIX                    = "_mean_score_NJ"
            + ForesterConstants.PHYLO_XML_SUFFIX;
    private static final String                                     NJ_TREE_SHARED_DOMAINS_BASED_GENOME_DISTANCE_SUFFIX                = "_domains_NJ"
            + ForesterConstants.PHYLO_XML_SUFFIX;
    private static final String                                     NJ_TREE_SHARED_BIN_COMBINATIONS_BASED_GENOME_DISTANCE_SUFFIX       = "_bin_combinations_NJ"
            + ForesterConstants.PHYLO_XML_SUFFIX;
    private static final String                                     FILTER_POSITIVE_OPTION                                             = "pos_filter";
    private static final String                                     FILTER_NEGATIVE_OPTION                                             = "neg_filter";
    private static final String                                     FILTER_NEGATIVE_DOMAINS_OPTION                                     = "neg_dom_filter";
    private static final String                                     INPUT_GENOMES_FILE_OPTION                                          = "genomes";
    private static final String                                     INPUT_SPECIES_TREE_OPTION                                          = "species_tree";
    private static final String                                     SEQ_EXTRACT_OPTION                                                 = "prot_extract";
    private static final boolean                                    IGNORE_DUFS_DEFAULT                                                = true;
    private static final boolean                                    IGNORE_COMBINATION_WITH_SAME_DEFAULLT                              = false;
    private static final double                                     MAX_E_VALUE_DEFAULT                                                = -1;
    private static final String                                     RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION                             = "random_seed";
    private static final String                                     CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS                           = "consider_bdc_direction";
    private static final String                                     CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS_AND_ADJACENCY             = "consider_bdc_adj";
    private static final String                                     OUTPUT_LIST_OF_ALL_PROTEINS_OPTIONS                                = "all_prot";
    private static final String                                     OUTPUT_LIST_OF_ALL_PROTEINS_PER_DOMAIN_E_VALUE_OPTION              = "all_prot_e";
    private static final String                                     OUTPUT_DOMAIN_COMBINATIONS_GAINED_MORE_THAN_ONCE_ANALYSIS_SUFFIX   = "_fitch_dc_gains_counts";
    private static final String                                     OUTPUT_DOMAIN_COMBINATIONS_LOST_MORE_THAN_ONCE_ANALYSIS_SUFFIX     = "_fitch_dc_losses_counts";
    private static final String                                     DOMAIN_LENGTHS_ANALYSIS_SUFFIX                                     = "_domain_lengths_analysis";
    private static final String                                     PERFORM_DOMAIN_LENGTH_ANALYSIS_OPTION                              = "dla";
    private static final String                                     D_PROMISCUITY_FILE_SUFFIX                                          = "_domain_promiscuities";
    private static final String                                     LOG_FILE_SUFFIX                                                    = "_log.txt";
    private static final String                                     DATA_FILE_SUFFIX                                                   = "_domain_combination_data.txt";
    private static final String                                     DATA_FILE_DESC                                                     = "#SPECIES\tPRTEIN_ID\tN_TERM_DOMAIN\tC_TERM_DOMAIN\tN_TERM_DOMAIN_PER_DOMAIN_E_VALUE\tC_TERM_DOMAIN_PER_DOMAIN_E_VALUE\tN_TERM_DOMAIN_COUNTS_PER_PROTEIN\tC_TERM_DOMAIN_COUNTS_PER_PROTEIN";
    private static final String                                     WRITE_TO_NEXUS_OPTION                                              = "nexus";
    private static final String                                     PERFORM_DC_FITCH                                                   = "dc_pars";
    private static final INDIVIDUAL_SCORE_CUTOFF                    INDIVIDUAL_SCORE_CUTOFF_DEFAULT                                    = INDIVIDUAL_SCORE_CUTOFF.FULL_SEQUENCE;                                                                                                                                                      //TODO look at me! change?
    private static final boolean                                    CALC_SIMILARITY_SCORES                                             = false;
    private static final String                                     SEPARATOR_FOR_DA                                                   = "--";
    private static final String                                     DOMAIN_SPECIES_IDS_MAP_NAME                                        = "_DOMAIN_SPECIES_IDS_MAP.txt";
    private static final String                                     WRITE_DA_IDS_NAMES_MAPS_OPTION                                     = "write_DA_maps";
    private static final String                                     INPUT_DA_NAME_FILE_OPTION                                          = "input_DA_name_map";
    private static final String                                     OBTAIN_NAMES_FOR_DAS_FROM_DB_OPTION                                = "obtain_DA_names_from_db";
    private static final String                                     VERBOSITY_OPTION                                                   = "verbosity";
    private static final int                                        VERBOSITY_DEFAULT                                                  = 0;
    private static final String                                     OBTAIN_NAMES_FOR_DAS_FROM_DB_MAX_IDS_TO_SEARCH_PER_SPECIES_OPTION  = "max_ids_to_search_per_species";
    private static final int                                        OBTAIN_NAMES_FOR_DAS_FROM_DB_MAX_IDS_TO_SEARCH_PER_SPECIES_DEFAULT = 20;
    private static final String                                     UNIPROT_PRIORITY_FOR_ACCESSOR_PARSING_OPTION                       = "uniprot_priority";
    private static final boolean                                    UNIPROT_PRIORITY_FOR_ACCESSOR_PARSING_DEFAULT                      = false;
    @SuppressWarnings( "unchecked")
    public static void main( final String args[] ) {
        final long start_time = new Date().getTime();
        final StringBuilder html_desc = new StringBuilder();
        ForesterUtil.printProgramInformation( SurfacingConstants.PRG_NAME,
                                              SurfacingConstants.PRG_VERSION,
                                              SurfacingConstants.PRG_DATE,
                                              SurfacingConstants.E_MAIL,
                                              SurfacingConstants.WWW );
        final String nl = ForesterUtil.LINE_SEPARATOR;
        html_desc.append( "<table>" + nl );
        html_desc.append( "<tr><td>Produced by:</td><td>" + SurfacingConstants.PRG_NAME + "</td></tr>" + nl );
        html_desc.append( "<tr><td>Version:</td><td>" + SurfacingConstants.PRG_VERSION + "</td></tr>" + nl );
        html_desc.append( "<tr><td>Release Date:</td><td>" + SurfacingConstants.PRG_DATE + "</td></tr>" + nl );
        html_desc.append( "<tr><td>Contact:</td><td>" + SurfacingConstants.E_MAIL + "</td></tr>" + nl );
        html_desc.append( "<tr><td>WWW:</td><td>" + SurfacingConstants.WWW + "</td></tr>" + nl );
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( launch.HELP_OPTION_1 ) || cla.isOptionSet( launch.HELP_OPTION_2 ) ) {
            launch.printHelp();
            System.exit( 0 );
        }
        if ( ( args.length < 1 ) ) {
            launch.printHelp();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<>();
        allowed_options.add( NOT_IGNORE_DUFS_OPTION );
        allowed_options.add( MAX_FS_E_VALUE_OPTION );
        allowed_options.add( MAX_I_E_VALUE_OPTION );
        allowed_options.add( MIN_REL_ENV_LENGTH_RATIO_OPTION );
        allowed_options.add( DETAILEDNESS_OPTION );
        allowed_options.add( OUTPUT_FILE_OPTION );
        allowed_options.add( DOMAIN_SIMILARITY_SORT_OPTION );
        allowed_options.add( SPECIES_MATRIX_OPTION );
        allowed_options.add( SCORING_OPTION );
        allowed_options.add( SurfacingConstants.MAX_ALLOWED_OVERLAP_OPTION );
        allowed_options.add( NO_ENGULFING_OVERLAP_OPTION );
        allowed_options.add( DOMAIN_COUNT_SORT_OPTION );
        allowed_options.add( CUTOFF_SCORE_FILE_OPTION );
        allowed_options.add( DOMAIN_SIMILARITY_SORT_BY_SPECIES_COUNT_FIRST_OPTION );
        allowed_options.add( OUTPUT_DIR_OPTION );
        allowed_options.add( IGNORE_COMBINATION_WITH_SAME_OPTION );
        allowed_options.add( PFAM_TO_GO_FILE_USE_OPTION );
        allowed_options.add( GO_OBO_FILE_USE_OPTION );
        allowed_options.add( DOMAIN_SIMILARITY_PRINT_OPTION );
        allowed_options.add( GO_NAMESPACE_LIMIT_OPTION );
        allowed_options.add( PAIRWISE_DOMAIN_COMPARISONS_OPTION );
        allowed_options.add( IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_OPTION );
        allowed_options.add( CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS );
        allowed_options.add( INPUT_SPECIES_TREE_OPTION );
        allowed_options.add( FILTER_POSITIVE_OPTION );
        allowed_options.add( FILTER_NEGATIVE_OPTION );
        allowed_options.add( INPUT_GENOMES_FILE_OPTION );
        allowed_options.add( RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION );
        allowed_options.add( FILTER_NEGATIVE_DOMAINS_OPTION );
        allowed_options.add( IGNORE_VIRAL_IDS );
        allowed_options.add( SEQ_EXTRACT_OPTION );
        allowed_options.add( OUTPUT_LIST_OF_ALL_PROTEINS_PER_DOMAIN_E_VALUE_OPTION );
        allowed_options.add( SECONDARY_FEATURES_PARSIMONY_MAP_FILE );
        allowed_options.add( SurfacingConstants.PLUS_MINUS_ANALYSIS_OPTION );
        allowed_options.add( DOMAIN_COMBINITONS_OUTPUT_OPTION_FOR_GRAPH_ANALYSIS );
        allowed_options.add( DOMAIN_COMBINITONS_COUNTS_OUTPUT_OPTION );
        allowed_options.add( OUTPUT_LIST_OF_ALL_PROTEINS_OPTIONS );
        allowed_options.add( CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS_AND_ADJACENCY );
        allowed_options.add( WRITE_TO_NEXUS_OPTION );
        allowed_options.add( PERFORM_DC_REGAIN_PROTEINS_STATS_OPTION );
        allowed_options.add( DA_ANALYSIS_OPTION );
        allowed_options.add( USE_LAST_IN_FITCH_OPTION );
        allowed_options.add( PERFORM_DC_FITCH );
        allowed_options.add( PERFORM_DOMAIN_LENGTH_ANALYSIS_OPTION );
        allowed_options.add( WRITE_DA_IDS_NAMES_MAPS_OPTION );
        allowed_options.add( INPUT_DA_NAME_FILE_OPTION );
        allowed_options.add( OBTAIN_NAMES_FOR_DAS_FROM_DB_OPTION );
        allowed_options.add( VERBOSITY_OPTION );
        allowed_options.add( OBTAIN_NAMES_FOR_DAS_FROM_DB_MAX_IDS_TO_SEARCH_PER_SPECIES_OPTION );
        allowed_options.add( UNIPROT_PRIORITY_FOR_ACCESSOR_PARSING_OPTION );
        boolean ignore_dufs = launch.IGNORE_DUFS_DEFAULT;
        boolean ignore_combination_with_same = launch.IGNORE_COMBINATION_WITH_SAME_DEFAULLT;
        double fs_e_value_max = launch.MAX_E_VALUE_DEFAULT;
        double ie_value_max = launch.MAX_E_VALUE_DEFAULT;
        double rel_env_length_ratio_cutoff = -1;
        int max_allowed_overlap = SurfacingConstants.MAX_ALLOWED_OVERLAP_DEFAULT;
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        boolean use_last_in_fitch_parsimony = false;
        if ( cla.isOptionSet( USE_LAST_IN_FITCH_OPTION ) ) {
            use_last_in_fitch_parsimony = true;
        }
        boolean write_to_nexus = false;
        if ( cla.isOptionSet( WRITE_TO_NEXUS_OPTION ) ) {
            write_to_nexus = true;
        }
        boolean perform_dc_fich = false;
        if ( cla.isOptionSet( PERFORM_DC_FITCH ) ) {
            perform_dc_fich = true;
        }
        boolean perform_dc_regain_proteins_stats = false;
        if ( cla.isOptionSet( PERFORM_DC_REGAIN_PROTEINS_STATS_OPTION ) ) {
            perform_dc_regain_proteins_stats = true;
        }
        boolean da_analysis = false;
        if ( cla.isOptionSet( DA_ANALYSIS_OPTION ) ) {
            da_analysis = true;
        }
        boolean output_binary_domain_combinationsfor_graph_analysis = false;
        if ( cla.isOptionSet( DOMAIN_COMBINITONS_OUTPUT_OPTION_FOR_GRAPH_ANALYSIS ) ) {
            output_binary_domain_combinationsfor_graph_analysis = true;
        }
        boolean output_binary_domain_combinationsfor_counts = false;
        if ( cla.isOptionSet( DOMAIN_COMBINITONS_COUNTS_OUTPUT_OPTION ) ) {
            output_binary_domain_combinationsfor_counts = true;
        }
        if ( cla.isOptionSet( launch.MAX_FS_E_VALUE_OPTION ) ) {
            try {
                fs_e_value_max = cla.getOptionValueAsDouble( launch.MAX_FS_E_VALUE_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "no acceptable value for E-value maximum" );
            }
        }
        if ( cla.isOptionSet( launch.MIN_REL_ENV_LENGTH_RATIO_OPTION ) ) {
            try {
                rel_env_length_ratio_cutoff = cla.getOptionValueAsDouble( launch.MIN_REL_ENV_LENGTH_RATIO_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no acceptable value for min rel env length ratio" );
            }
        }
        if ( cla.isOptionSet( launch.MAX_I_E_VALUE_OPTION ) ) {
            try {
                ie_value_max = cla.getOptionValueAsDouble( launch.MAX_I_E_VALUE_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "no acceptable value for E-value maximum" );
            }
        }
        if ( cla.isOptionSet( SurfacingConstants.MAX_ALLOWED_OVERLAP_OPTION ) ) {
            try {
                max_allowed_overlap = cla.getOptionValueAsInt( SurfacingConstants.MAX_ALLOWED_OVERLAP_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no acceptable value for maximal allowed domain overlap" );
            }
        }
        boolean no_engulfing_overlaps = false;
        if ( cla.isOptionSet( launch.NO_ENGULFING_OVERLAP_OPTION ) ) {
            no_engulfing_overlaps = true;
        }
        boolean ignore_virus_like_ids = false;
        if ( cla.isOptionSet( launch.IGNORE_VIRAL_IDS ) ) {
            ignore_virus_like_ids = true;
        }
        if ( cla.isOptionSet( launch.NOT_IGNORE_DUFS_OPTION ) ) {
            ignore_dufs = false;
        }
        if ( cla.isOptionSet( launch.IGNORE_COMBINATION_WITH_SAME_OPTION ) ) {
            ignore_combination_with_same = true;
        }
        boolean domain_length_analysis = false;
        if ( cla.isOptionSet( launch.PERFORM_DOMAIN_LENGTH_ANALYSIS_OPTION ) ) {
            domain_length_analysis = true;
        }
        boolean ignore_domains_without_combs_in_all_spec = IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_DEFAULT;
        if ( cla.isOptionSet( launch.IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_OPTION ) ) {
            ignore_domains_without_combs_in_all_spec = true;
        }
        boolean ignore_species_specific_domains = IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION_DEFAULT;
        if ( cla.isOptionSet( launch.IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION ) ) {
            ignore_species_specific_domains = true;
        }
        if ( !cla.isOptionValueSet( launch.INPUT_SPECIES_TREE_OPTION ) ) {
            ForesterUtil
                    .fatalError( SurfacingConstants.PRG_NAME,
                                 "no input species tree file given: " + launch.INPUT_SPECIES_TREE_OPTION + "=<file>" );
        }
        File output_file = null;
        if ( cla.isOptionSet( launch.OUTPUT_FILE_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.OUTPUT_FILE_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for domain combinations similarities output file: -"
                                                 + launch.OUTPUT_FILE_OPTION + "=<file>" );
            }
            output_file = new File( cla.getOptionValue( launch.OUTPUT_FILE_OPTION ) );
            SurfacingUtil.checkForOutputFileWriteability( output_file );
        }
        File cutoff_scores_file = null;
        Map<String, Double> individual_score_cutoffs = null;
        if ( cla.isOptionSet( launch.CUTOFF_SCORE_FILE_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.CUTOFF_SCORE_FILE_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for individual domain score cutoffs file: -"
                                                 + launch.CUTOFF_SCORE_FILE_OPTION + "=<file>" );
            }
            cutoff_scores_file = new File( cla.getOptionValue( launch.CUTOFF_SCORE_FILE_OPTION ) );
            final String error = ForesterUtil.isReadableFile( cutoff_scores_file );
            if ( !ForesterUtil.isEmpty( error ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "cannot read individual domain score cutoffs file: " + error );
            }
            try {
                final BasicTable<String> scores_table = BasicTableParser.parse( cutoff_scores_file, ' ' );
                individual_score_cutoffs = scores_table.getColumnsAsMapDouble( 0, 1 );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "cannot read from individual score cutoffs file: " + e );
            }
        }
        BinaryDomainCombination.DomainCombinationType dc_type = BinaryDomainCombination.DomainCombinationType.BASIC;
        if ( cla.isOptionSet( launch.CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS ) ) {
            dc_type = BinaryDomainCombination.DomainCombinationType.DIRECTED;
        }
        if ( cla.isOptionSet( launch.CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS_AND_ADJACENCY ) ) {
            dc_type = BinaryDomainCombination.DomainCombinationType.DIRECTED_ADJACTANT;
        }
        File out_dir = null;
        if ( cla.isOptionSet( launch.OUTPUT_DIR_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.OUTPUT_DIR_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for output directory: -" + launch.OUTPUT_DIR_OPTION + "=<dir>" );
            }
            out_dir = new File( cla.getOptionValue( launch.OUTPUT_DIR_OPTION ) );
            if ( out_dir.exists() && ( out_dir.listFiles().length > 0 ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "\"" + out_dir + "\" aready exists and is not empty" );
            }
            if ( !out_dir.exists() ) {
                final boolean success = out_dir.mkdir();
                if ( !success || !out_dir.exists() ) {
                    ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "failed to create \"" + out_dir + "\"" );
                }
            }
            if ( !out_dir.canWrite() ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "cannot write to \"" + out_dir + "\"" );
            }
        }
        File positive_filter_file = null;
        File negative_filter_file = null;
        File negative_domains_filter_file = null;
        if ( cla.isOptionSet( launch.FILTER_NEGATIVE_OPTION ) && cla.isOptionSet( launch.FILTER_POSITIVE_OPTION ) ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                     "attempt to use both negative and positive protein filter" );
        }
        if ( cla.isOptionSet( launch.FILTER_NEGATIVE_DOMAINS_OPTION )
                && ( cla.isOptionSet( launch.FILTER_NEGATIVE_OPTION )
                        || cla.isOptionSet( launch.FILTER_POSITIVE_OPTION ) ) ) {
            ForesterUtil
                    .fatalError( SurfacingConstants.PRG_NAME,
                                 "attempt to use both negative or positive protein filter together wirh a negative domains filter" );
        }
        if ( cla.isOptionSet( launch.FILTER_NEGATIVE_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.FILTER_NEGATIVE_OPTION ) ) {
                ForesterUtil
                        .fatalError( SurfacingConstants.PRG_NAME,
                                     "no value for negative filter: -" + launch.FILTER_NEGATIVE_OPTION + "=<file>" );
            }
            negative_filter_file = new File( cla.getOptionValue( launch.FILTER_NEGATIVE_OPTION ) );
            final String msg = ForesterUtil.isReadableFile( negative_filter_file );
            if ( !ForesterUtil.isEmpty( msg ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "can not read from \"" + negative_filter_file + "\": " + msg );
            }
        }
        else if ( cla.isOptionSet( launch.FILTER_POSITIVE_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.FILTER_POSITIVE_OPTION ) ) {
                ForesterUtil
                        .fatalError( SurfacingConstants.PRG_NAME,
                                     "no value for positive filter: -" + launch.FILTER_POSITIVE_OPTION + "=<file>" );
            }
            positive_filter_file = new File( cla.getOptionValue( launch.FILTER_POSITIVE_OPTION ) );
            final String msg = ForesterUtil.isReadableFile( positive_filter_file );
            if ( !ForesterUtil.isEmpty( msg ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "can not read from \"" + positive_filter_file + "\": " + msg );
            }
        }
        else if ( cla.isOptionSet( launch.FILTER_NEGATIVE_DOMAINS_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.FILTER_NEGATIVE_DOMAINS_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for negative domains filter: -"
                                                 + launch.FILTER_NEGATIVE_DOMAINS_OPTION + "=<file>" );
            }
            negative_domains_filter_file = new File( cla.getOptionValue( launch.FILTER_NEGATIVE_DOMAINS_OPTION ) );
            final String msg = ForesterUtil.isReadableFile( negative_domains_filter_file );
            if ( !ForesterUtil.isEmpty( msg ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "can not read from \"" + negative_domains_filter_file + "\": " + msg );
            }
        }
        final List<String> plus_minus_analysis_high_copy_base_species = new ArrayList<>();
        final List<String> plus_minus_analysis_high_copy_target_species = new ArrayList<>();
        final List<String> plus_minus_analysis_high_low_copy_species = new ArrayList<>();
        final List<Object> plus_minus_analysis_numbers = new ArrayList<>();
        SurfacingUtil.processPlusMinusAnalysisOption( cla,
                                                      plus_minus_analysis_high_copy_base_species,
                                                      plus_minus_analysis_high_copy_target_species,
                                                      plus_minus_analysis_high_low_copy_species,
                                                      plus_minus_analysis_numbers );
        File input_genomes_file = null;
        if ( cla.isOptionSet( launch.INPUT_GENOMES_FILE_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.INPUT_GENOMES_FILE_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for input genomes file: -" + launch.INPUT_GENOMES_FILE_OPTION
                                                 + "=<file>" );
            }
            input_genomes_file = new File( cla.getOptionValue( launch.INPUT_GENOMES_FILE_OPTION ) );
            final String msg = ForesterUtil.isReadableFile( input_genomes_file );
            if ( !ForesterUtil.isEmpty( msg ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "can not read from \"" + input_genomes_file + "\": " + msg );
            }
        }
        else {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                     "no input genomes file given: " + launch.INPUT_GENOMES_FILE_OPTION + "=<file>" );
        }
        DomainSimilarity.DomainSimilarityScoring scoring = SCORING_DEFAULT;
        if ( cla.isOptionSet( launch.SCORING_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.SCORING_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for scoring method for domain combinations similarity calculation: -"
                                                 + launch.SCORING_OPTION + "=<" + launch.SCORING_DOMAIN_COUNT_BASED
                                                 + "|" + launch.SCORING_PROTEIN_COUNT_BASED + "|"
                                                 + launch.SCORING_COMBINATION_BASED + ">\"" );
            }
            final String scoring_str = cla.getOptionValue( launch.SCORING_OPTION );
            if ( scoring_str.equals( launch.SCORING_DOMAIN_COUNT_BASED ) ) {
                scoring = DomainSimilarity.DomainSimilarityScoring.DOMAINS;
            }
            else if ( scoring_str.equals( launch.SCORING_COMBINATION_BASED ) ) {
                scoring = DomainSimilarity.DomainSimilarityScoring.COMBINATIONS;
            }
            else if ( scoring_str.equals( launch.SCORING_PROTEIN_COUNT_BASED ) ) {
                scoring = DomainSimilarity.DomainSimilarityScoring.PROTEINS;
            }
            else {
                ForesterUtil
                        .fatalError( SurfacingConstants.PRG_NAME,
                                     "unknown value \"" + scoring_str
                                             + "\" for scoring method for domain combinations similarity calculation: \"-"
                                             + launch.SCORING_OPTION + "=<" + launch.SCORING_DOMAIN_COUNT_BASED + "|"
                                             + launch.SCORING_PROTEIN_COUNT_BASED + "|"
                                             + launch.SCORING_COMBINATION_BASED + ">\"" );
            }
        }
        boolean sort_by_species_count_first = false;
        if ( cla.isOptionSet( launch.DOMAIN_SIMILARITY_SORT_BY_SPECIES_COUNT_FIRST_OPTION ) ) {
            sort_by_species_count_first = true;
        }
        boolean species_matrix = false;
        if ( cla.isOptionSet( launch.SPECIES_MATRIX_OPTION ) ) {
            species_matrix = true;
        }
        boolean output_protein_lists_for_all_domains = false;
        double output_list_of_all_proteins_per_domain_e_value_max = -1;
        if ( cla.isOptionSet( launch.OUTPUT_LIST_OF_ALL_PROTEINS_OPTIONS ) ) {
            output_protein_lists_for_all_domains = true;
            if ( cla.isOptionSet( launch.OUTPUT_LIST_OF_ALL_PROTEINS_PER_DOMAIN_E_VALUE_OPTION ) ) {
                try {
                    output_list_of_all_proteins_per_domain_e_value_max = cla
                            .getOptionValueAsDouble( launch.OUTPUT_LIST_OF_ALL_PROTEINS_PER_DOMAIN_E_VALUE_OPTION );
                }
                catch ( final Exception e ) {
                    ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                             "no acceptable value for per domain E-value maximum" );
                }
            }
        }
        Detailedness detailedness = DETAILEDNESS_DEFAULT;
        if ( cla.isOptionSet( launch.DETAILEDNESS_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.DETAILEDNESS_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for -" + launch.DETAILEDNESS_OPTION + "=<"
                                                 + launch.DETAILEDNESS_BASIC + "|" + launch.DETAILEDNESS_LIST_IDS + "|"
                                                 + launch.DETAILEDNESS_PUNCTILIOUS + ">\"" );
            }
            final String detness = cla.getOptionValue( launch.DETAILEDNESS_OPTION ).toLowerCase();
            if ( detness.equals( launch.DETAILEDNESS_BASIC ) ) {
                detailedness = DomainSimilarityCalculator.Detailedness.BASIC;
            }
            else if ( detness.equals( launch.DETAILEDNESS_LIST_IDS ) ) {
                detailedness = DomainSimilarityCalculator.Detailedness.LIST_COMBINING_DOMAIN_FOR_EACH_SPECIES;
            }
            else if ( detness.equals( launch.DETAILEDNESS_PUNCTILIOUS ) ) {
                detailedness = DomainSimilarityCalculator.Detailedness.PUNCTILIOUS;
            }
            else {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "unknown value \"" + detness + "\" for detailedness: \"-"
                                                 + launch.DETAILEDNESS_OPTION + "=<" + launch.DETAILEDNESS_BASIC + "|"
                                                 + launch.DETAILEDNESS_LIST_IDS + "|" + launch.DETAILEDNESS_PUNCTILIOUS
                                                 + ">\"" );
            }
        }
        String automated_pairwise_comparison_suffix = null;
        boolean perform_pwc = false;
        boolean write_pwc_files = false;
        if ( cla.isOptionSet( launch.PAIRWISE_DOMAIN_COMPARISONS_OPTION ) ) {
            perform_pwc = true;
            if ( !cla.isOptionValueSet( launch.PAIRWISE_DOMAIN_COMPARISONS_OPTION ) ) {
                write_pwc_files = false;
            }
            else {
                write_pwc_files = true;
                automated_pairwise_comparison_suffix = "_"
                        + cla.getOptionValue( launch.PAIRWISE_DOMAIN_COMPARISONS_OPTION );
            }
        }
        String query_domain_ids = null;
        if ( cla.isOptionSet( launch.SEQ_EXTRACT_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.SEQ_EXTRACT_OPTION ) ) {
                ForesterUtil
                        .fatalError( SurfacingConstants.PRG_NAME,
                                     "no domain ids given for sequences with given domains to be extracted : -"
                                             + launch.SEQ_EXTRACT_OPTION
                                             + "=<ordered domain sequences, domain ids separated by '~', sequences separated by '#'>" );
            }
            query_domain_ids = cla.getOptionValue( launch.SEQ_EXTRACT_OPTION );
        }
        DomainSimilarity.DomainSimilaritySortField domain_similarity_sort_field = DOMAIN_SORT_FILD_DEFAULT;
        DomainSimilarity.DomainSimilaritySortField domain_similarity_sort_field_for_automated_pwc = DOMAIN_SORT_FILD_DEFAULT;
        if ( cla.isOptionSet( launch.DOMAIN_SIMILARITY_SORT_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.DOMAIN_SIMILARITY_SORT_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for domain combinations similarities sorting: -"
                                                 + launch.DOMAIN_SIMILARITY_SORT_OPTION + "=<"
                                                 + launch.DOMAIN_SIMILARITY_SORT_ALPHA + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_MAX + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_MIN + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_MEAN + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_DIFF + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_ABS_COUNTS_DIFF + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_COUNTS_DIFF + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_SPECIES_COUNT + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_SD + ">\"" );
            }
            final String sort_str = cla.getOptionValue( launch.DOMAIN_SIMILARITY_SORT_OPTION ).toLowerCase();
            if ( sort_str.equals( launch.DOMAIN_SIMILARITY_SORT_ALPHA ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( launch.DOMAIN_SIMILARITY_SORT_MAX ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MAX;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( launch.DOMAIN_SIMILARITY_SORT_MIN ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MIN;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( launch.DOMAIN_SIMILARITY_SORT_MEAN ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MEAN;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.MEAN;
            }
            else if ( sort_str.equals( launch.DOMAIN_SIMILARITY_SORT_SPECIES_COUNT ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.SPECIES_COUNT;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( launch.DOMAIN_SIMILARITY_SORT_SD ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.SD;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( launch.DOMAIN_SIMILARITY_SORT_DIFF ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MAX_DIFFERENCE;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.MAX_DIFFERENCE;
            }
            else if ( sort_str.equals( launch.DOMAIN_SIMILARITY_SORT_ABS_COUNTS_DIFF ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE;
            }
            else if ( sort_str.equals( launch.DOMAIN_SIMILARITY_SORT_COUNTS_DIFF ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE;
            }
            else {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "unknown value \"" + sort_str
                                                 + "\" for domain combinations similarities sorting: \"-"
                                                 + launch.DOMAIN_SIMILARITY_SORT_OPTION + "=<"
                                                 + launch.DOMAIN_SIMILARITY_SORT_ALPHA + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_MAX + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_MIN + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_MEAN + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_DIFF + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_ABS_COUNTS_DIFF + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_COUNTS_DIFF + "|" + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_SPECIES_COUNT + "|"
                                                 + launch.DOMAIN_SIMILARITY_SORT_SD + ">\"" );
            }
        }
        DomainSimilarity.PRINT_OPTION domain_similarity_print_option = DOMAIN_SIMILARITY_PRINT_OPTION_DEFAULT;
        if ( cla.isOptionSet( launch.DOMAIN_SIMILARITY_PRINT_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.DOMAIN_SIMILARITY_PRINT_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for print option: -"
                                                 + launch.DOMAIN_SIMILARITY_PRINT_OPTION_DETAILED_HTML + "|"
                                                 + launch.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_HTML + "|"
                                                 + launch.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_TAB_DELIMITED + ">\"" );
            }
            final String sort = cla.getOptionValue( launch.DOMAIN_SIMILARITY_PRINT_OPTION ).toLowerCase();
            if ( sort.equals( launch.DOMAIN_SIMILARITY_PRINT_OPTION_DETAILED_HTML ) ) {
                domain_similarity_print_option = DomainSimilarity.PRINT_OPTION.HTML;
            }
            else if ( sort.equals( launch.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_HTML ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "simple HTML output not implemented yet :(" );
            }
            else if ( sort.equals( launch.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_TAB_DELIMITED ) ) {
                domain_similarity_print_option = DomainSimilarity.PRINT_OPTION.SIMPLE_TAB_DELIMITED;
            }
            else {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "unknown value \"" + sort + "\" for print option: -"
                                                 + launch.DOMAIN_SIMILARITY_PRINT_OPTION_DETAILED_HTML + "|"
                                                 + launch.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_HTML + "|"
                                                 + launch.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_TAB_DELIMITED + ">\"" );
            }
        }
        GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder dc_sort_order = DOMAINS_SORT_ORDER_DEFAULT;
        if ( cla.isOptionSet( launch.DOMAIN_COUNT_SORT_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.DOMAIN_COUNT_SORT_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for sorting of domain counts: -" + launch.DOMAIN_COUNT_SORT_OPTION
                                                 + "=<" + launch.DOMAIN_COUNT_SORT_ALPHA + "|"
                                                 + launch.DOMAIN_COUNT_SORT_KEY_DOMAIN_COUNT + "|"
                                                 + launch.DOMAIN_COUNT_SORT_KEY_DOMAIN_PROTEINS_COUNT + "|"
                                                 + launch.DOMAIN_COUNT_SORT_COMBINATIONS_COUNT + ">\"" );
            }
            final String sort = cla.getOptionValue( launch.DOMAIN_COUNT_SORT_OPTION ).toLowerCase();
            if ( sort.equals( launch.DOMAIN_COUNT_SORT_ALPHA ) ) {
                dc_sort_order = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.ALPHABETICAL_KEY_ID;
            }
            else if ( sort.equals( launch.DOMAIN_COUNT_SORT_KEY_DOMAIN_COUNT ) ) {
                dc_sort_order = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.KEY_DOMAIN_COUNT;
            }
            else if ( sort.equals( launch.DOMAIN_COUNT_SORT_KEY_DOMAIN_PROTEINS_COUNT ) ) {
                dc_sort_order = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.KEY_DOMAIN_PROTEINS_COUNT;
            }
            else if ( sort.equals( launch.DOMAIN_COUNT_SORT_COMBINATIONS_COUNT ) ) {
                dc_sort_order = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.COMBINATIONS_COUNT;
            }
            else {
                ForesterUtil
                        .fatalError( SurfacingConstants.PRG_NAME,
                                     "unknown value \"" + sort + "\" for sorting of domain counts: \"-"
                                             + launch.DOMAIN_COUNT_SORT_OPTION + "=<" + launch.DOMAIN_COUNT_SORT_ALPHA
                                             + "|" + launch.DOMAIN_COUNT_SORT_KEY_DOMAIN_COUNT + "|"
                                             + launch.DOMAIN_COUNT_SORT_KEY_DOMAIN_PROTEINS_COUNT + "|"
                                             + launch.DOMAIN_COUNT_SORT_COMBINATIONS_COUNT + ">\"" );
            }
        }
        final String[][] input_file_properties = SurfacingUtil.processInputGenomesFile( input_genomes_file );
        final int number_of_genomes = input_file_properties.length;
        if ( number_of_genomes < 2 ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "cannot analyze less than two files" );
        }
        if ( ( number_of_genomes < 3 ) && perform_pwc ) {
            ForesterUtil
                    .fatalError( SurfacingConstants.PRG_NAME,
                                 "cannot use : -" + launch.PAIRWISE_DOMAIN_COMPARISONS_OPTION
                                         + "=<suffix> to turn on pairwise analyses with less than three input files" );
        }
        SurfacingUtil.checkWriteabilityForPairwiseComparisons( domain_similarity_print_option,
                                                               input_file_properties,
                                                               automated_pairwise_comparison_suffix,
                                                               out_dir );
        for( int i = 0; i < number_of_genomes; i++ ) {
            File dcc_outfile = new File( input_file_properties[ i ][ 1 ]
                    + SurfacingConstants.DOMAIN_COMBINITON_COUNTS_OUTPUTFILE_SUFFIX );
            if ( out_dir != null ) {
                dcc_outfile = new File( out_dir + ForesterUtil.FILE_SEPARATOR + dcc_outfile );
            }
            SurfacingUtil.checkForOutputFileWriteability( dcc_outfile );
        }
        File pfam_to_go_file = new File( "pfam2go.txt" );
        if ( cla.isOptionSet( launch.PFAM_TO_GO_FILE_USE_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.PFAM_TO_GO_FILE_USE_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for Pfam to GO mapping file: -" + launch.PFAM_TO_GO_FILE_USE_OPTION
                                                 + "=<file>" );
            }
            pfam_to_go_file = new File( cla.getOptionValue( launch.PFAM_TO_GO_FILE_USE_OPTION ) );
        }
        final String error1 = ForesterUtil.isReadableFile( pfam_to_go_file );
        if ( !ForesterUtil.isEmpty( error1 ) ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "cannot read Pfam to GO mapping file: " + error1 );
        }
        Map<String, List<GoId>> domain_id_to_go_ids_map = null;
        int domain_id_to_go_ids_count = 0;
        try {
            final PfamToGoParser parser = new PfamToGoParser( pfam_to_go_file );
            final List<PfamToGoMapping> pfam_to_go_mappings = parser.parse();
            domain_id_to_go_ids_map = SurfacingUtil.createDomainIdToGoIdMap( pfam_to_go_mappings );
            if ( parser.getMappingCount() < domain_id_to_go_ids_map.size() ) {
                ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME,
                                                   "parser.getMappingCount() < domain_id_to_go_ids_map.size()" );
            }
            domain_id_to_go_ids_count = parser.getMappingCount();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "cannot read from Pfam to GO mapping file: " + e );
        }
        File go_obo_file = new File( "go.obo" );
        if ( cla.isOptionSet( launch.GO_OBO_FILE_USE_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.GO_OBO_FILE_USE_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for GO OBO file: -" + launch.GO_OBO_FILE_USE_OPTION + "=<file>" );
            }
            go_obo_file = new File( cla.getOptionValue( launch.GO_OBO_FILE_USE_OPTION ) );
        }
        final String error2 = ForesterUtil.isReadableFile( go_obo_file );
        if ( !ForesterUtil.isEmpty( error2 ) ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "cannot read GO OBO file: " + error2 );
        }
        List<GoTerm> go_terms = null;
        try {
            final OBOparser parser = new OBOparser( go_obo_file, OBOparser.ReturnType.BASIC_GO_TERM );
            go_terms = parser.parse();
            if ( parser.getGoTermCount() != go_terms.size() ) {
                ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME,
                                                   "parser.getGoTermCount() != go_terms.size()" );
            }
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "cannot read from GO OBO file: " + e );
        }
        Map<GoId, GoTerm> go_id_to_term_map = null;
        if ( ( ( domain_id_to_go_ids_map != null ) && ( domain_id_to_go_ids_map.size() > 0 ) )
                && ( ( go_terms != null ) && ( go_terms.size() > 0 ) ) ) {
            go_id_to_term_map = GoUtils.createGoIdToGoTermMap( go_terms );
        }
        GoNameSpace go_namespace_limit = null;
        if ( cla.isOptionSet( launch.GO_NAMESPACE_LIMIT_OPTION ) ) {
            if ( ( go_id_to_term_map == null ) || go_id_to_term_map.isEmpty() ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "cannot use GO namespace limit (-" + launch.GO_NAMESPACE_LIMIT_OPTION
                                                 + "=<namespace>) without Pfam to GO mapping file ("
                                                 + launch.PFAM_TO_GO_FILE_USE_OPTION + "=<file>) and GO OBO file (-"
                                                 + launch.GO_OBO_FILE_USE_OPTION + "=<file>)" );
            }
            if ( !cla.isOptionValueSet( launch.GO_NAMESPACE_LIMIT_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for GO namespace limit: \"-" + launch.GO_NAMESPACE_LIMIT_OPTION
                                                 + "=<" + launch.GO_NAMESPACE_LIMIT_OPTION_MOLECULAR_FUNCTION + "|"
                                                 + launch.GO_NAMESPACE_LIMIT_OPTION_BIOLOGICAL_PROCESS + "|"
                                                 + launch.GO_NAMESPACE_LIMIT_OPTION_CELLULAR_COMPONENT + ">\"" );
            }
            final String go_namespace_limit_str = cla.getOptionValue( launch.GO_NAMESPACE_LIMIT_OPTION ).toLowerCase();
            if ( go_namespace_limit_str.equals( launch.GO_NAMESPACE_LIMIT_OPTION_MOLECULAR_FUNCTION ) ) {
                go_namespace_limit = GoNameSpace.createMolecularFunction();
            }
            else if ( go_namespace_limit_str.equals( launch.GO_NAMESPACE_LIMIT_OPTION_BIOLOGICAL_PROCESS ) ) {
                go_namespace_limit = GoNameSpace.createBiologicalProcess();
            }
            else if ( go_namespace_limit_str.equals( launch.GO_NAMESPACE_LIMIT_OPTION_CELLULAR_COMPONENT ) ) {
                go_namespace_limit = GoNameSpace.createCellularComponent();
            }
            else {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "unknown value \"" + go_namespace_limit_str + "\" for GO namespace limit: \"-"
                                                 + launch.GO_NAMESPACE_LIMIT_OPTION + "=<"
                                                 + launch.GO_NAMESPACE_LIMIT_OPTION_MOLECULAR_FUNCTION + "|"
                                                 + launch.GO_NAMESPACE_LIMIT_OPTION_BIOLOGICAL_PROCESS + "|"
                                                 + launch.GO_NAMESPACE_LIMIT_OPTION_CELLULAR_COMPONENT + ">\"" );
            }
        }
        if ( ( domain_similarity_sort_field == DomainSimilarity.DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE )
                && ( number_of_genomes > 2 ) ) {
            domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE;
        }
        File[] intree_files = null;
        Phylogeny[] intrees = null;
        if ( cla.isOptionSet( launch.INPUT_SPECIES_TREE_OPTION ) ) {
            if ( number_of_genomes < 3 ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "cannot infer gains and losses on input species trees (-"
                                                 + launch.INPUT_SPECIES_TREE_OPTION + " without pairwise analyses ("
                                                 + launch.PAIRWISE_DOMAIN_COMPARISONS_OPTION
                                                 + "=<suffix for pairwise comparison output files>)" );
            }
            if ( !cla.isOptionValueSet( launch.INPUT_SPECIES_TREE_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for input tree: -" + launch.INPUT_SPECIES_TREE_OPTION
                                                 + "=<tree file in phyloXML format>" );
            }
            final String intrees_str = cla.getOptionValue( launch.INPUT_SPECIES_TREE_OPTION );
            if ( intrees_str.indexOf( "#" ) > 0 ) {
                final String[] intrees_strs = intrees_str.split( "#" );
                intree_files = new File[ intrees_strs.length ];
                int i = 0;
                for( final String s : intrees_strs ) {
                    intree_files[ i++ ] = new File( s.trim() );
                }
            }
            else {
                intree_files = new File[ 1 ];
                intree_files[ 0 ] = new File( intrees_str );
            }
            intrees = SurfacingUtil
                    .obtainAndPreProcessIntrees( intree_files, number_of_genomes, input_file_properties );
        }
        final Phylogeny intree_0_orig = SurfacingUtil.obtainFirstIntree( intree_files[ 0 ] );
        long random_number_seed_for_fitch_parsimony = 0l;
        boolean radomize_fitch_parsimony = false;
        if ( cla.isOptionSet( launch.RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION ) ) {
            if ( !cla.isOptionValueSet( launch.RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for random number seed: -"
                                                 + launch.RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION + "=<seed>" );
            }
            try {
                random_number_seed_for_fitch_parsimony = cla
                        .getOptionValueAsLong( RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.getMessage() );
            }
            radomize_fitch_parsimony = true;
        }
        SortedSet<String> filter = null;
        if ( ( positive_filter_file != null ) || ( negative_filter_file != null )
                || ( negative_domains_filter_file != null ) ) {
            filter = new TreeSet<>();
            if ( positive_filter_file != null ) {
                SurfacingUtil.processFilter( positive_filter_file, filter );
            }
            else if ( negative_filter_file != null ) {
                SurfacingUtil.processFilter( negative_filter_file, filter );
            }
            else if ( negative_domains_filter_file != null ) {
                SurfacingUtil.processFilter( negative_domains_filter_file, filter );
            }
        }
        final boolean obtain_names_for_das_from_db;
        if ( cla.isOptionSet( launch.OBTAIN_NAMES_FOR_DAS_FROM_DB_OPTION ) ) {
            obtain_names_for_das_from_db = true;
        }
        else {
            obtain_names_for_das_from_db = false;
        }
        final boolean write_da_ids_names_maps;
        if ( cla.isOptionSet( WRITE_DA_IDS_NAMES_MAPS_OPTION ) ) {
            write_da_ids_names_maps = true;
        }
        else {
            write_da_ids_names_maps = false;
        }
        if ( cla.isOptionSet( UNIPROT_PRIORITY_FOR_ACCESSOR_PARSING_OPTION ) ) {
            GlobalOptions.setUniprotPriorityForAccessorParsing( true );
        }
        else {
            GlobalOptions.setUniprotPriorityForAccessorParsing( UNIPROT_PRIORITY_FOR_ACCESSOR_PARSING_DEFAULT );
        }
        if ( cla.isOptionSet( OBTAIN_NAMES_FOR_DAS_FROM_DB_MAX_IDS_TO_SEARCH_PER_SPECIES_OPTION ) ) {
            if ( !cla.isOptionValueSet( OBTAIN_NAMES_FOR_DAS_FROM_DB_MAX_IDS_TO_SEARCH_PER_SPECIES_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value given for: -"
                                                 + OBTAIN_NAMES_FOR_DAS_FROM_DB_MAX_IDS_TO_SEARCH_PER_SPECIES_OPTION
                                                 + "=<int>" );
            }
            try {
                GlobalOptions.setObtainNamesForDasFromDbMaxIdsToSearchPerSpecies( cla
                        .getOptionValueAsInt( OBTAIN_NAMES_FOR_DAS_FROM_DB_MAX_IDS_TO_SEARCH_PER_SPECIES_OPTION ) );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.getMessage() );
            }
        }
        else {
            GlobalOptions
                    .setObtainNamesForDasFromDbMaxIdsToSearchPerSpecies( OBTAIN_NAMES_FOR_DAS_FROM_DB_MAX_IDS_TO_SEARCH_PER_SPECIES_DEFAULT );
        }
        if ( cla.isOptionSet( VERBOSITY_OPTION ) ) {
            if ( !cla.isOptionValueSet( VERBOSITY_OPTION ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value given for: -" + VERBOSITY_OPTION + "=<int>" );
            }
            try {
                GlobalOptions.setVerbosity( cla.getOptionValueAsInt( VERBOSITY_OPTION ) );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.getMessage() );
            }
        }
        else {
            GlobalOptions.setVerbosity( VERBOSITY_DEFAULT );
        }
        File input_da_name_file = null;
        if ( write_da_ids_names_maps ) {
            if ( cla.isOptionSet( launch.INPUT_DA_NAME_FILE_OPTION ) ) {
                if ( !cla.isOptionValueSet( launch.INPUT_DA_NAME_FILE_OPTION ) ) {
                    ForesterUtil
                            .fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for GO OBO file: -" + launch.INPUT_DA_NAME_FILE_OPTION + "=<file>" );
                }
                input_da_name_file = new File( cla.getOptionValue( launch.INPUT_DA_NAME_FILE_OPTION ) );
                final String error4 = ForesterUtil.isReadableFile( input_da_name_file );
                if ( !ForesterUtil.isEmpty( error4 ) ) {
                    ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, "cannot read: " + input_da_name_file );
                }
            }
            if ( ( input_da_name_file == null ) && !obtain_names_for_das_from_db ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "need to obtain names for DAs either from file and/or database" );
            }
        }
        Map<String, Set<String>>[] domain_id_to_secondary_features_maps = null;
        File[] secondary_features_map_files = null;
        final File domain_lengths_analysis_outfile = new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file
                + DOMAIN_LENGTHS_ANALYSIS_SUFFIX );
        if ( domain_length_analysis ) {
            SurfacingUtil.checkForOutputFileWriteability( domain_lengths_analysis_outfile );
        }
        if ( cla.isOptionSet( launch.SECONDARY_FEATURES_PARSIMONY_MAP_FILE ) ) {
            if ( !cla.isOptionValueSet( launch.SECONDARY_FEATURES_PARSIMONY_MAP_FILE ) ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                         "no value for secondary features map file: -"
                                                 + launch.SECONDARY_FEATURES_PARSIMONY_MAP_FILE + "=<file>" );
            }
            final String[] secondary_features_map_files_strs = cla
                    .getOptionValue( launch.SECONDARY_FEATURES_PARSIMONY_MAP_FILE ).split( "#" );
            secondary_features_map_files = new File[ secondary_features_map_files_strs.length ];
            domain_id_to_secondary_features_maps = new Map[ secondary_features_map_files_strs.length ];
            int i = 0;
            for( final String secondary_features_map_files_str : secondary_features_map_files_strs ) {
                secondary_features_map_files[ i ] = new File( secondary_features_map_files_str );
                final String error = ForesterUtil.isReadableFile( secondary_features_map_files[ i ] );
                if ( !ForesterUtil.isEmpty( error ) ) {
                    ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                             "cannot read secondary features map file: " + error );
                }
                try {
                    domain_id_to_secondary_features_maps[ i ] = SurfacingUtil
                            .createDomainIdToSecondaryFeaturesMap( secondary_features_map_files[ i ] );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                             "cannot read secondary features map file: " + e.getMessage() );
                }
                catch ( final Exception e ) {
                    ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                             "problem with contents of features map file ["
                                                     + secondary_features_map_files[ i ] + "]: " + e.getMessage() );
                }
                i++;
            }
        }
        if ( out_dir == null ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                     "no output directory indicated (-" + launch.OUTPUT_DIR_OPTION + "=<dir>)" );
        }
        if ( output_file == null ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                     "no name for (main) output file indicated (-" + launch.OUTPUT_FILE_OPTION
                                             + "=<file>)" );
        }
        if ( ( domain_id_to_go_ids_map == null ) || domain_id_to_go_ids_map.isEmpty() ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                     "no (acceptable) Pfam to GO id mapping file provided ('pfam2go file') (-"
                                             + launch.PFAM_TO_GO_FILE_USE_OPTION + "=<file>)" );
        }
        if ( ( go_id_to_term_map == null ) || go_id_to_term_map.isEmpty() ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                     "no (acceptable) go id to term mapping file provided ('GO OBO file') (-"
                                             + launch.GO_OBO_FILE_USE_OPTION + "=<file>)" );
        }
        System.out.println( "Output directory            : " + out_dir );
        System.out.println( "Input genomes from          : " + input_genomes_file );
        html_desc.append( "<tr><td>Input genomes from:</td><td>" + input_genomes_file + "</td></tr>" + nl );
        if ( positive_filter_file != null ) {
            final int filter_size = filter.size();
            System.out.println( "Positive protein filter     : " + positive_filter_file + " [" + filter_size
                    + " domain ids]" );
            html_desc.append( "<tr><td>Positive protein filter:</td><td>" + positive_filter_file + " [" + filter_size
                    + " domain ids]</td></tr>" + nl );
        }
        if ( negative_filter_file != null ) {
            final int filter_size = filter.size();
            System.out.println( "Negative protein filter     : " + negative_filter_file + " [" + filter_size
                    + " domain ids]" );
            html_desc.append( "<tr><td>Negative protein filter:</td><td>" + negative_filter_file + " [" + filter_size
                    + " domain ids]</td></tr>" + nl );
        }
        if ( negative_domains_filter_file != null ) {
            final int filter_size = filter.size();
            System.out.println( "Negative domain filter      : " + negative_domains_filter_file + " [" + filter_size
                    + " domain ids]" );
            html_desc.append( "<tr><td>Negative domain filter:</td><td>" + negative_domains_filter_file + " ["
                    + filter_size + " domain ids]</td></tr>" + nl );
        }
        if ( plus_minus_analysis_high_copy_base_species.size() > 0 ) {
            String plus0 = "";
            for( final String s : plus_minus_analysis_high_copy_base_species ) {
                plus0 += "+" + s + " ";
            }
            String plus1 = "";
            for( final String s : plus_minus_analysis_high_copy_target_species ) {
                plus1 += "*" + s + " ";
            }
            String minus = "";
            for( final String s : plus_minus_analysis_high_low_copy_species ) {
                minus += "-" + s + " ";
            }
            System.out.println( "Plus-minus analysis         : " + plus1 + "&& " + plus0 + "&& " + minus );
            html_desc.append( "<tr><td>Plus-minus analysis:</td><td>" + plus1 + "&& " + plus0 + "&& " + minus
                    + "</td></tr>" + nl );
        }
        if ( cutoff_scores_file != null ) {
            System.out.println( "Cutoff scores file          : " + cutoff_scores_file );
            html_desc.append( "<tr><td>Cutoff scores file:</td><td>" + cutoff_scores_file + "</td></tr>" + nl );
        }
        if ( ie_value_max >= 0.0 ) {
            System.out.println( "iE-value maximum (incl)     : " + ie_value_max );
            html_desc.append( "<tr><td>iE-value maximum (inclusive):</td><td>" + ie_value_max + "</td></tr>" + nl );
        }
        if ( rel_env_length_ratio_cutoff > 0.0 ) {
            System.out.println( "Rel env length ratio min    : " + rel_env_length_ratio_cutoff );
            html_desc.append( "<tr><td>Relative hmm envelope length ratio min (inclusive):</td><td>"
                    + rel_env_length_ratio_cutoff + "</td></tr>" + nl );
        }
        if ( fs_e_value_max >= 0.0 ) {
            System.out.println( "FS E-value maximum (incl)   : " + fs_e_value_max );
            html_desc.append( "<tr><td>FS E-value maximum (inclusive):</td><td>" + fs_e_value_max + "</td></tr>" + nl );
        }
        if ( output_protein_lists_for_all_domains ) {
            System.out.println( "Domain E-value max          : " + output_list_of_all_proteins_per_domain_e_value_max );
            html_desc.append( "<tr><td>Protein lists: E-value maximum per domain (inclusive):</td><td>"
                    + output_list_of_all_proteins_per_domain_e_value_max + "</td></tr>" + nl );
        }
        System.out.println( "Ignore DUFs                 : " + ignore_dufs );
        if ( ignore_virus_like_ids ) {
            System.out.println( "Ignore virus like ids       : " + ignore_virus_like_ids );
            html_desc.append( "<tr><td>Ignore virus, phage, transposition related ids:</td><td>" + ignore_virus_like_ids
                    + "</td></tr>" + nl );
        }
        html_desc.append( "<tr><td>Ignore DUFs:</td><td>" + ignore_dufs + "</td></tr>" + nl );
        if ( max_allowed_overlap != SurfacingConstants.MAX_ALLOWED_OVERLAP_DEFAULT ) {
            System.out.println( "Max allowed domain overlap  : " + max_allowed_overlap );
            html_desc
                    .append( "<tr><td>Max allowed domain overlap:</td><td>" + max_allowed_overlap + "</td></tr>" + nl );
        }
        if ( no_engulfing_overlaps ) {
            System.out.println( "Ignore engulfed domains     : " + no_engulfing_overlaps );
            html_desc.append( "<tr><td>Ignore (lower confidence) engulfed domains:</td><td>" + no_engulfing_overlaps
                    + "</td></tr>" + nl );
        }
        System.out.println( "Ignore singlet domains      : " + ignore_domains_without_combs_in_all_spec );
        html_desc
                .append( "<tr><td>Ignore singlet domains for domain combination similarity analyses (not for parsimony analyses):</td><td>"
                        + ignore_domains_without_combs_in_all_spec + "</td></tr>" + nl );
        System.out.println( "Ignore species specific doms: " + ignore_species_specific_domains );
        html_desc
                .append( "<tr><td>Ignore species specific domains for domain combination similarity analyses (not for parsimony analyses):</td><td>"
                        + ignore_species_specific_domains + "</td></tr>" + nl );
        System.out.println( "Ignore combination with self: " + ignore_combination_with_same );
        html_desc.append( "<tr><td>Ignore combination with self for domain combination similarity analyses:</td><td>"
                + ignore_combination_with_same + "</td></tr>" + nl );
        System.out.println( "Consider directedness       : "
                + ( dc_type != BinaryDomainCombination.DomainCombinationType.BASIC ) );
        html_desc.append( "<tr><td>Consider directedness of binary domain combinations:</td><td>"
                + ( dc_type != BinaryDomainCombination.DomainCombinationType.BASIC ) + "</td></tr>" + nl );
        if ( dc_type != BinaryDomainCombination.DomainCombinationType.BASIC ) {
            System.out.println( "Consider adjacency          : "
                    + ( dc_type == BinaryDomainCombination.DomainCombinationType.DIRECTED_ADJACTANT ) );
            html_desc.append( "<tr><td>Consider djacency of binary domain combinations:</td><td>"
                    + ( dc_type == BinaryDomainCombination.DomainCombinationType.DIRECTED_ADJACTANT ) + "</td></tr>"
                    + nl );
        }
        System.out.println( "Fitch parsimony of DCs      : " + perform_dc_fich );
        html_desc.append( "<tr><td>Fitch parsimony of DCs:</td><td>" + perform_dc_fich + "</td></tr>" + nl );
        if ( perform_dc_fich ) {
            System.out.println( "Use last in Fitch parsimony : " + use_last_in_fitch_parsimony );
            html_desc.append( "<tr><td>Use last in Fitch parsimony:</td><td>" + use_last_in_fitch_parsimony
                    + "</td></tr>" + nl );
        }
        System.out.println( "Write to Nexus files        : " + write_to_nexus );
        html_desc.append( "<tr><td>Write to Nexus files:</td><td>" + write_to_nexus + "</td></tr>" + nl );
        if ( perform_dc_fich ) {
            System.out.println( "DC regain prot stats        : " + perform_dc_regain_proteins_stats );
            html_desc.append( "<tr><td>DC regain prot stats:</td><td>" + perform_dc_regain_proteins_stats + "</td></tr>"
                    + nl );
        }
        System.out.println( "DA analysis                 : " + da_analysis );
        html_desc.append( "<tr><td>DA analysis :</td><td>" + da_analysis + "</td></tr>" + nl );
        System.out.print( "Domain counts sort order    : " );
        html_desc.append( "<tr><td>Domain counts sort order:</td><td>" );
        switch ( dc_sort_order ) {
            case ALPHABETICAL_KEY_ID:
                System.out.println( "alphabetical" );
                html_desc.append( "alphabetical" + "</td></tr>" + nl );
                break;
            case KEY_DOMAIN_COUNT:
                System.out.println( "domain count" );
                html_desc.append( "domain count" + "</td></tr>" + nl );
                break;
            case KEY_DOMAIN_PROTEINS_COUNT:
                System.out.println( "domain proteins count" );
                html_desc.append( "domain proteins count" + "</td></tr>" + nl );
                break;
            case COMBINATIONS_COUNT:
                System.out.println( "domain combinations count" );
                html_desc.append( "domain combinations count" + "</td></tr>" + nl );
                break;
            default:
                ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME, "unknown value for dc sort order" );
        }
        if ( domain_id_to_go_ids_map != null ) {
            System.out.println( "Pfam to GO mappings from    : " + pfam_to_go_file + " [" + domain_id_to_go_ids_count
                    + " mappings]" );
            html_desc.append( "<tr><td>Pfam to GO mappings from:</td><td>" + pfam_to_go_file + " ["
                    + domain_id_to_go_ids_count + " mappings]" + "</td></tr>" + nl );
        }
        if ( go_terms != null ) {
            System.out.println( "GO terms from               : " + go_obo_file + " [" + go_terms.size() + " terms]" );
            html_desc.append( "<tr><td>GO terms from:</td><td>" + go_obo_file + " [" + go_terms.size() + " terms]"
                    + "</td></tr>" + nl );
        }
        if ( go_namespace_limit != null ) {
            System.out.println( "Limit GO terms to           : " + go_namespace_limit.toString() );
            html_desc.append( "<tr><td>Limit GO terms to</td><td>" + go_namespace_limit + "</td></tr>" + nl );
        }
        if ( perform_pwc ) {
            System.out.println( "Suffix for PWC files        : " + automated_pairwise_comparison_suffix );
            html_desc.append( "<tr><td>Suffix for PWC files</td><td>" + automated_pairwise_comparison_suffix
                    + "</td></tr>" + nl );
        }
        if ( out_dir != null ) {
            System.out.println( "Output directory            : " + out_dir );
        }
        if ( query_domain_ids != null ) {
            System.out.println( "Query domains (ordered)     : " + query_domain_ids );
            html_desc.append( "<tr><td></td><td>" + query_domain_ids + "</td></tr>" + nl );
        }
        System.out.println( "Write similarities to       : " + output_file );
        System.out.print( "  Scoring method            : " );
        html_desc.append( "<tr><td>Scoring method:</td><td>" );
        switch ( scoring ) {
            case COMBINATIONS:
                System.out.println( "domain combinations based" );
                html_desc.append( "domain combinations based" + "</td></tr>" + nl );
                break;
            case DOMAINS:
                System.out.println( "domain counts based" );
                html_desc.append( "domain counts based" + "</td></tr>" + nl );
                break;
            case PROTEINS:
                System.out.println( "domain proteins counts based" );
                html_desc.append( "domain proteins counts based" + "</td></tr>" + nl );
                break;
            default:
                ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME,
                                                   "unknown value for sorting for scoring" );
        }
        System.out.print( "  Sort by                   : " );
        html_desc.append( "<tr><td>Sort by:</td><td>" );
        switch ( domain_similarity_sort_field ) {
            case MIN:
                System.out.print( "score minimum" );
                html_desc.append( "score minimum" );
                break;
            case MAX:
                System.out.print( "score maximum" );
                html_desc.append( "score maximum" );
                break;
            case MEAN:
                System.out.print( "score mean" );
                html_desc.append( "score mean" );
                break;
            case SD:
                System.out.print( "score standard deviation" );
                html_desc.append( "score standard deviation" );
                break;
            case SPECIES_COUNT:
                System.out.print( "species number" );
                html_desc.append( "species number" );
                break;
            case DOMAIN_ID:
                System.out.print( "alphabetical domain identifier" );
                html_desc.append( "alphabetical domain identifier" );
                break;
            case MAX_DIFFERENCE:
                System.out.print( "(maximal) difference" );
                html_desc.append( "(maximal) difference" );
                break;
            case ABS_MAX_COUNTS_DIFFERENCE:
                System.out.print( "absolute (maximal) counts difference" );
                html_desc.append( "absolute (maximal) counts difference" );
                break;
            case MAX_COUNTS_DIFFERENCE:
                System.out.print( "(maximal) counts difference" );
                html_desc.append( "(maximal) counts  difference" );
                break;
            default:
                ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME,
                                                   "unknown value for sorting for similarities" );
        }
        if ( sort_by_species_count_first ) {
            System.out.println( " (sort by species count first)" );
            html_desc.append( " (sort by species count first)" );
        }
        else {
            System.out.println();
        }
        html_desc.append( "</td></tr>" + nl );
        System.out.print( "  Detailedness              : " );
        switch ( detailedness ) {
            case BASIC:
                System.out.println( "basic" );
                break;
            case LIST_COMBINING_DOMAIN_FOR_EACH_SPECIES:
                System.out.println( "list combining domains for each species" );
                break;
            case PUNCTILIOUS:
                System.out.println( "punctilious" );
                break;
            default:
                ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME,
                                                   "unknown value for sorting for detailedness" );
        }
        System.out.print( "  Print option              : " );
        switch ( domain_similarity_print_option ) {
            case HTML:
                System.out.println( "HTML" );
                break;
            case SIMPLE_TAB_DELIMITED:
                System.out.println( "simple tab delimited" );
                break;
            default:
                ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME, "unknown value for print option" );
        }
        System.out.print( "  Species matrix            : " + species_matrix );
        System.out.println();
        final File dc_data_file = new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file + DATA_FILE_SUFFIX );
        System.out.println( "Domain comb data output     : " + dc_data_file );
        html_desc.append( "<tr><td>Domain combination data output:</td><td> " + dc_data_file + " </td></tr>" );
        System.out.println();
        if ( perform_pwc ) {
            System.out.println( "Pairwise comparisons: " );
            html_desc.append( "<tr><td>Pairwise comparisons:</td><td></td></tr>" );
            System.out.print( "  Sort by                   : " );
            html_desc.append( "<tr><td>Sort by:</td><td>" );
            switch ( domain_similarity_sort_field_for_automated_pwc ) {
                case MEAN:
                    System.out.print( "score mean" );
                    html_desc.append( "score mean" );
                    break;
                case DOMAIN_ID:
                    System.out.print( "alphabetical domain identifier" );
                    html_desc.append( "alphabetical domain identifier" );
                    break;
                case MAX_DIFFERENCE:
                    System.out.print( "difference" );
                    html_desc.append( "difference" );
                    break;
                case ABS_MAX_COUNTS_DIFFERENCE:
                    System.out.print( "absolute counts difference" );
                    html_desc.append( "absolute counts difference" );
                    break;
                case MAX_COUNTS_DIFFERENCE:
                    System.out.print( "counts difference" );
                    html_desc.append( "counts difference" );
                    break;
                default:
                    ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME,
                                                       "unknown value for sorting for similarities" );
            }
            System.out.println();
            html_desc.append( "</td></tr>" + nl );
            if ( ( intrees != null ) && ( intrees.length > 0 ) ) {
                for( final File intree_file : intree_files ) {
                    html_desc.append( "<tr><td>Intree for gain/loss parsimony analysis:</td><td>" + intree_file
                            + "</td></tr>" + nl );
                    System.out.println( "  Intree for gain/loss pars.: " + intree_file );
                }
            }
            if ( radomize_fitch_parsimony ) {
                html_desc.append( "<tr><td>    Random number seed for Fitch parsimony analysis:</td><td>"
                        + random_number_seed_for_fitch_parsimony + "</td></tr>" + nl );
                System.out.println( "    Random number seed      : " + random_number_seed_for_fitch_parsimony );
            }
            if ( ( domain_id_to_secondary_features_maps != null )
                    && ( domain_id_to_secondary_features_maps.length > 0 ) ) {
                for( int i = 0; i < secondary_features_map_files.length; i++ ) {
                    html_desc.append( "<tr><td>Secondary features map file:</td><td>"
                            + secondary_features_map_files[ i ] + "</td></tr>" + nl );
                    System.out.println( "Secondary features map file : " + secondary_features_map_files[ i ]
                            + " [mappings for " + domain_id_to_secondary_features_maps[ i ].size() + " domain ids]" );
                    if ( GlobalOptions.getVerbosity() > 0 ) {
                        System.out.println();
                        System.out.println( "Domain ids to secondary features map:" );
                        for( final String domain_id : domain_id_to_secondary_features_maps[ i ].keySet() ) {
                            System.out.print( domain_id );
                            System.out.print( " => " );
                            for( final String sec : domain_id_to_secondary_features_maps[ i ].get( domain_id ) ) {
                                System.out.print( sec );
                                System.out.print( " " );
                            }
                            System.out.println();
                        }
                    }
                }
            }
        } // if ( perform_pwc ) {
        if ( obtain_names_for_das_from_db ) {
            System.out.println( "Obtain DA names from online : true" );
            html_desc.append( "<tr><td>Obtain DA names from online database (UniProtKB):</td><td>true</td></tr>" + nl );
            System.out.println( "Max IDs per species         : "
                    + GlobalOptions.getObtainNamesForDasFromDbMaxIdsToSearchPerSpecies() );
            html_desc.append( "<tr><td>Max IDs per species for DB search :</td><td>"
                    + GlobalOptions.getObtainNamesForDasFromDbMaxIdsToSearchPerSpecies() + "</td></tr>" + nl );
            System.out
                    .println( "UniProt priority            : " + GlobalOptions.isUniprotPriorityForAccessorParsing() );
            html_desc.append( "<tr><td>\"UniProt priority:</td><td>"
                    + GlobalOptions.isUniprotPriorityForAccessorParsing() + "</td></tr>" + nl );
        }
        System.out.println( "Verbosity                   : " + GlobalOptions.getVerbosity() );
        html_desc.append( "<tr><td>Verbosity :</td><td>" + GlobalOptions.getVerbosity() + "</td></tr>" + nl );
        html_desc.append( "<tr><td>Command line:</td><td>" + nl + nl + cla.getCommandLineArgsAsString() + nl + nl
                + "</td></tr>" + nl );
        System.out.println( "Command line                : " + cla.getCommandLineArgsAsString() );
        BufferedWriter[] query_domains_writer_ary = null;
        List<String>[] query_domain_ids_array = null;
        if ( query_domain_ids != null ) {
            final String[] query_domain_ids_str_array = query_domain_ids.split( "#" );
            query_domain_ids_array = new ArrayList[ query_domain_ids_str_array.length ];
            query_domains_writer_ary = new BufferedWriter[ query_domain_ids_str_array.length ];
            for( int i = 0; i < query_domain_ids_str_array.length; i++ ) {
                String query_domain_ids_str = query_domain_ids_str_array[ i ];
                final String[] query_domain_ids_str_ary = query_domain_ids_str.split( "~" );
                final List<String> query = new ArrayList<>();
                for( final String element : query_domain_ids_str_ary ) {
                    query.add( element );
                }
                query_domain_ids_array[ i ] = query;
                query_domain_ids_str = query_domain_ids_str.replace( '~', '_' );
                String protein_names_writer_str = query_domain_ids_str + SurfacingConstants.SEQ_EXTRACT_SUFFIX;
                if ( out_dir != null ) {
                    protein_names_writer_str = out_dir + ForesterUtil.FILE_SEPARATOR + protein_names_writer_str;
                }
                try {
                    query_domains_writer_ary[ i ] = new BufferedWriter( new FileWriter( protein_names_writer_str ) );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( SurfacingConstants.PRG_NAME,
                                             "Could not open [" + protein_names_writer_str + "]: "
                                                     + e.getLocalizedMessage() );
                }
            }
        }
        SortedMap<Species, List<Protein>> protein_lists_per_species = null; //This will only be created if needed.
        boolean need_protein_lists_per_species = false;
        //if ( ( plus_minus_analysis_high_copy_base_species.size() > 0 ) || output_protein_lists_for_all_domains
        //    || true ) { //TODO
        need_protein_lists_per_species = true;
        //  }
        if ( need_protein_lists_per_species ) {
            protein_lists_per_species = new TreeMap<>();
        }
        List<GenomeWideCombinableDomains> gwcd_list = new ArrayList<>( number_of_genomes );
        final SortedSet<String> all_domains_encountered = new TreeSet<>();
        final SortedSet<BinaryDomainCombination> all_bin_domain_combinations_encountered = new TreeSet<>();
        List<BinaryDomainCombination> all_bin_domain_combinations_gained_fitch = null;
        List<BinaryDomainCombination> all_bin_domain_combinations_lost_fitch = null;
        if ( ( intrees != null ) && ( intrees.length == 1 ) ) {
            all_bin_domain_combinations_gained_fitch = new ArrayList<>();
            all_bin_domain_combinations_lost_fitch = new ArrayList<>();
        }
        final File per_genome_domain_promiscuity_statistics_file = new File( out_dir + ForesterUtil.FILE_SEPARATOR
                + output_file + D_PROMISCUITY_FILE_SUFFIX );
        BufferedWriter per_genome_domain_promiscuity_statistics_writer = null;
        try {
            per_genome_domain_promiscuity_statistics_writer = new BufferedWriter( new FileWriter( per_genome_domain_promiscuity_statistics_file ) );
            per_genome_domain_promiscuity_statistics_writer.write( "Species:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "Mean:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "SD:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "Median:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "Min:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "Max:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "N:\t" );
            per_genome_domain_promiscuity_statistics_writer
                    .write( "Max Promiscuous Domains:" + ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final IOException e2 ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e2.getMessage() );
        }
        final File log_file = new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file + LOG_FILE_SUFFIX );
        BufferedWriter log_writer = null;
        try {
            log_writer = new BufferedWriter( new FileWriter( log_file ) );
        }
        catch ( final IOException e2 ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e2.getMessage() );
        }
        BufferedWriter dc_data_writer = null;
        try {
            dc_data_writer = new BufferedWriter( new FileWriter( dc_data_file ) );
            dc_data_writer.write( DATA_FILE_DESC );
            dc_data_writer.write( ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final IOException e2 ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e2.getMessage() );
        }
        DescriptiveStatistics protein_coverage_stats = new BasicDescriptiveStatistics();
        DescriptiveStatistics all_genomes_domains_per_potein_stats = new BasicDescriptiveStatistics();
        final SortedMap<Integer, Integer> all_genomes_domains_per_potein_histo = new TreeMap<>();
        final SortedSet<String> domains_which_are_always_single = new TreeSet<>();
        final SortedSet<String> domains_which_are_sometimes_single_sometimes_not = new TreeSet<>();
        final SortedSet<String> domains_which_never_single = new TreeSet<>();
        BufferedWriter domains_per_potein_stats_writer = null;
        try {
            domains_per_potein_stats_writer = new BufferedWriter( new FileWriter( out_dir + ForesterUtil.FILE_SEPARATOR
                    + output_file + "_domains_per_potein_stats.txt" ) );
            domains_per_potein_stats_writer.write( "Genome" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer.write( "Mean" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer.write( "SD" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer.write( "Median" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer.write( "N" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer.write( "Min" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer.write( "Max" );
            domains_per_potein_stats_writer.write( "\n" );
        }
        catch ( final IOException e3 ) {
            e3.printStackTrace();
        }
        Map<String, DescriptiveStatistics> protein_length_stats_by_dc = null;
        Map<String, DescriptiveStatistics> domain_number_stats_by_dc = null;
        final Map<String, DescriptiveStatistics> domain_length_stats_by_domain = new HashMap<>();
        if ( perform_dc_regain_proteins_stats ) {
            protein_length_stats_by_dc = new HashMap<>();
            domain_number_stats_by_dc = new HashMap<>();
        }
        DomainLengthsTable domain_lengths_table = null;
        if ( domain_length_analysis ) {
            domain_lengths_table = new DomainLengthsTable();
        }
        // Main loop:
        final SortedMap<String, Set<String>> distinct_domain_architecutures_per_genome = new TreeMap<>();
        final SortedMap<String, Integer> distinct_domain_architecuture_counts = new TreeMap<>();
        for( int i = 0; i < number_of_genomes; ++i ) {
            System.out.println();
            System.out.println( ( i + 1 ) + "/" + number_of_genomes );
            SurfacingUtil.log( ( i + 1 ) + "/" + number_of_genomes, log_writer );
            System.out.println( "Processing                                     : " + input_file_properties[ i ][ 1 ]
                    + " [" + input_file_properties[ i ][ 0 ] + "]" );
            SurfacingUtil.log(
                               "Genome                                         : " + input_file_properties[ i ][ 1 ]
                                       + " [" + input_file_properties[ i ][ 0 ] + "]",
                               log_writer );
            HmmscanPerDomainTableParser parser = null;
            INDIVIDUAL_SCORE_CUTOFF ind_score_cutoff = INDIVIDUAL_SCORE_CUTOFF.NONE;
            if ( individual_score_cutoffs != null ) {
                ind_score_cutoff = INDIVIDUAL_SCORE_CUTOFF_DEFAULT;
            }
            if ( ( positive_filter_file != null ) || ( negative_filter_file != null )
                    || ( negative_domains_filter_file != null ) ) {
                HmmscanPerDomainTableParser.FilterType filter_type = HmmscanPerDomainTableParser.FilterType.NONE;
                if ( positive_filter_file != null ) {
                    filter_type = HmmscanPerDomainTableParser.FilterType.POSITIVE_PROTEIN;
                }
                else if ( negative_filter_file != null ) {
                    filter_type = HmmscanPerDomainTableParser.FilterType.NEGATIVE_PROTEIN;
                }
                else if ( negative_domains_filter_file != null ) {
                    filter_type = HmmscanPerDomainTableParser.FilterType.NEGATIVE_DOMAIN;
                }
                parser = new HmmscanPerDomainTableParser( new File( input_file_properties[ i ][ 0 ] ),
                                                          input_file_properties[ i ][ 1 ],
                                                          filter,
                                                          filter_type,
                                                          ind_score_cutoff,
                                                          true );
            }
            else {
                parser = new HmmscanPerDomainTableParser( new File( input_file_properties[ i ][ 0 ] ),
                                                          input_file_properties[ i ][ 1 ],
                                                          ind_score_cutoff,
                                                          true );
            }
            if ( fs_e_value_max >= 0.0 ) {
                parser.setFsEValueMaximum( fs_e_value_max );
            }
            if ( ie_value_max >= 0.0 ) {
                parser.setIEValueMaximum( ie_value_max );
            }
            if ( rel_env_length_ratio_cutoff > 0.0 ) {
                parser.setRelEnvLengthRatioCutoff( rel_env_length_ratio_cutoff );
            }
            parser.setIgnoreDufs( ignore_dufs );
            parser.setIgnoreVirusLikeIds( ignore_virus_like_ids );
            parser.setIgnoreEngulfedDomains( no_engulfing_overlaps );
            if ( max_allowed_overlap != SurfacingConstants.MAX_ALLOWED_OVERLAP_DEFAULT ) {
                parser.setMaxAllowedOverlap( max_allowed_overlap );
            }
            parser.setReturnType( HmmscanPerDomainTableParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            if ( individual_score_cutoffs != null ) {
                parser.setIndividualScoreCutoffs( individual_score_cutoffs );
            }
            List<Protein> protein_list = null;
            try {
                protein_list = parser.parse();
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.getMessage() );
            }
            catch ( final Exception e ) {
                ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME, e.getMessage(), e );
            }
            if ( GlobalOptions.getVerbosity() > 0 ) {
                System.out.println( "Domains ignored due to negative domain filter: " );
                ForesterUtil.printCountingMap( parser.getDomainsIgnoredDueToNegativeDomainFilterCountsMap() );
                System.out.println( "Domains ignored due to virus like id: " );
                ForesterUtil.printCountingMap( parser.getDomainsIgnoredDueToVirusLikeIdCountsMap() );
            }
            final double coverage = ( double ) protein_list.size() / parser.getProteinsEncountered();
            protein_coverage_stats.addValue( coverage );
            int distinct_das = -1;
            if ( da_analysis ) {
                final String genome = input_file_properties[ i ][ 0 ];
                distinct_das = SurfacingUtil.storeDomainArchitectures( genome,
                                                                       distinct_domain_architecutures_per_genome,
                                                                       protein_list,
                                                                       distinct_domain_architecuture_counts );
            }
            System.out.println( "Number of proteins encountered                 : " + parser.getProteinsEncountered() );
            SurfacingUtil.log( "Number of proteins encountered                 : " + parser.getProteinsEncountered(),
                               log_writer );
            System.out.println( "Number of proteins stored                      : " + protein_list.size() );
            SurfacingUtil.log( "Number of proteins stored                      : " + protein_list.size(), log_writer );
            System.out.println( "Coverage                                       : "
                    + ForesterUtil.roundToInt( 100.0 * coverage ) + "%" );
            SurfacingUtil.log(
                               "Coverage                                       : "
                                       + ForesterUtil.roundToInt( 100.0 * coverage ) + "%",
                               log_writer );
            System.out.println( "Domains encountered                            : " + parser.getDomainsEncountered() );
            SurfacingUtil.log( "Domains encountered                            : " + parser.getDomainsEncountered(),
                               log_writer );
            System.out.println( "Domains stored                                 : " + parser.getDomainsStored() );
            SurfacingUtil.log( "Domains stored                                 : " + parser.getDomainsStored(),
                               log_writer );
            System.out.println( "Distinct domains stored                        : "
                    + parser.getDomainsStoredSet().size() );
            SurfacingUtil
                    .log( "Distinct domains stored                        : " + parser.getDomainsStoredSet().size(),
                          log_writer );
            System.out.println( "Domains ignored due to individual score cutoffs: "
                    + parser.getDomainsIgnoredDueToIndividualScoreCutoff() );
            SurfacingUtil.log(
                               "Domains ignored due to individual score cutoffs: "
                                       + parser.getDomainsIgnoredDueToIndividualScoreCutoff(),
                               log_writer );
            System.out.println( "Domains ignored due to FS E-value              : "
                    + parser.getDomainsIgnoredDueToFsEval() );
            SurfacingUtil
                    .log( "Domains ignored due to FS E-value              : " + parser.getDomainsIgnoredDueToFsEval(),
                          log_writer );
            System.out.println( "Domains ignored due to iE-value                : "
                    + parser.getDomainsIgnoredDueToIEval() );
            SurfacingUtil
                    .log( "Domains ignored due to iE-value                : " + parser.getDomainsIgnoredDueToIEval(),
                          log_writer );
            System.out.println( "Domains ignored due to rel env length ratio    : "
                    + parser.getDomainsIgnoredDueToRelEnvLengthRatioCutoff() );
            SurfacingUtil.log(
                               "Domains ignored due to rel env length ratio    : "
                                       + parser.getDomainsIgnoredDueToRelEnvLengthRatioCutoff(),
                               log_writer );
            System.out.println( "Domains ignored due to DUF designation         : "
                    + parser.getDomainsIgnoredDueToDuf() );
            SurfacingUtil.log( "Domains ignored due to DUF designation         : " + parser.getDomainsIgnoredDueToDuf(),
                               log_writer );
            if ( ignore_virus_like_ids ) {
                System.out.println( "Domains ignored due virus like ids             : "
                        + parser.getDomainsIgnoredDueToVirusLikeIds() );
                SurfacingUtil.log(
                                   "Domains ignored due virus like ids             : "
                                           + parser.getDomainsIgnoredDueToVirusLikeIds(),
                                   log_writer );
            }
            System.out.println( "Domains ignored due negative domain filter     : "
                    + parser.getDomainsIgnoredDueToNegativeDomainFilter() );
            SurfacingUtil.log(
                               "Domains ignored due negative domain filter     : "
                                       + parser.getDomainsIgnoredDueToNegativeDomainFilter(),
                               log_writer );
            System.out.println( "Domains ignored due to overlap                 : "
                    + parser.getDomainsIgnoredDueToOverlap() );
            SurfacingUtil
                    .log( "Domains ignored due to overlap                 : " + parser.getDomainsIgnoredDueToOverlap(),
                          log_writer );
            if ( negative_filter_file != null ) {
                System.out.println( "Proteins ignored due to negative filter        : "
                        + parser.getProteinsIgnoredDueToFilter() );
                SurfacingUtil.log(
                                   "Proteins ignored due to negative filter        : "
                                           + parser.getProteinsIgnoredDueToFilter(),
                                   log_writer );
            }
            if ( positive_filter_file != null ) {
                System.out.println( "Proteins ignored due to positive filter        : "
                        + parser.getProteinsIgnoredDueToFilter() );
                SurfacingUtil.log(
                                   "Proteins ignored due to positive filter        : "
                                           + parser.getProteinsIgnoredDueToFilter(),
                                   log_writer );
            }
            if ( da_analysis ) {
                System.out.println( "Distinct domain architectures stored           : " + distinct_das );
                SurfacingUtil.log( "Distinct domain architectures stored           : " + distinct_das, log_writer );
            }
            System.out.println( "Time for processing                            : " + parser.getTime() + "ms" );
            SurfacingUtil.log( "", log_writer );
            try {
                int count = 0;
                for( final Protein protein : protein_list ) {
                    dc_data_writer
                            .write( SurfacingUtil.proteinToDomainCombinations( protein, count + "", "\t" ).toString() );
                    ++count;
                    for( final Domain d : protein.getProteinDomains() ) {
                        final String d_str = d.getDomainId().toString();
                        if ( !domain_length_stats_by_domain.containsKey( d_str ) ) {
                            domain_length_stats_by_domain.put( d_str, new BasicDescriptiveStatistics() );
                        }
                        domain_length_stats_by_domain.get( d_str ).addValue( d.getLength() );
                    }
                }
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.toString() );
            }
            SurfacingUtil.domainsPerProteinsStatistics( input_file_properties[ i ][ 1 ],
                                                        protein_list,
                                                        all_genomes_domains_per_potein_stats,
                                                        all_genomes_domains_per_potein_histo,
                                                        domains_which_are_always_single,
                                                        domains_which_are_sometimes_single_sometimes_not,
                                                        domains_which_never_single,
                                                        domains_per_potein_stats_writer );
            if ( domain_length_analysis ) {
                domain_lengths_table.addLengths( protein_list );
            }
            if ( !da_analysis ) {
                gwcd_list.add( BasicGenomeWideCombinableDomains
                        .createInstance( protein_list,
                                         ignore_combination_with_same,
                                         new BasicSpecies( input_file_properties[ i ][ 1 ] ),
                                         domain_id_to_go_ids_map,
                                         dc_type,
                                         protein_length_stats_by_dc,
                                         domain_number_stats_by_dc ) );
                if ( gwcd_list.get( i ).getSize() > 0 ) {
                    if ( output_binary_domain_combinationsfor_counts ) {
                        SurfacingUtil
                                .writeDomainCombinationsCountsFile( input_file_properties,
                                                                    out_dir,
                                                                    per_genome_domain_promiscuity_statistics_writer,
                                                                    gwcd_list.get( i ),
                                                                    i,
                                                                    dc_sort_order );
                    }
                    if ( output_binary_domain_combinationsfor_graph_analysis ) {
                        SurfacingUtil.writeBinaryDomainCombinationsFileForGraphAnalysis( input_file_properties,
                                                                                         out_dir,
                                                                                         gwcd_list.get( i ),
                                                                                         i,
                                                                                         dc_sort_order );
                    }
                    SurfacingUtil.addAllDomainIdsToSet( gwcd_list.get( i ), all_domains_encountered );
                    SurfacingUtil.addAllBinaryDomainCombinationToSet( gwcd_list.get( i ),
                                                                      all_bin_domain_combinations_encountered );
                }
            }
            if ( query_domains_writer_ary != null ) {
                for( int j = 0; j < query_domain_ids_array.length; j++ ) {
                    try {
                        SurfacingUtil.extractProteinNames( protein_list,
                                                           query_domain_ids_array[ j ],
                                                           query_domains_writer_ary[ j ],
                                                           "\t",
                                                           SurfacingConstants.LIMIT_SPEC_FOR_PROT_EX );
                        query_domains_writer_ary[ j ].flush();
                    }
                    catch ( final IOException e ) {
                        e.printStackTrace();
                    }
                }
            }
            if ( need_protein_lists_per_species ) {
                protein_lists_per_species.put( new BasicSpecies( input_file_properties[ i ][ 1 ] ), protein_list );
            }
            try {
                log_writer.flush();
            }
            catch ( final IOException e2 ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e2.getLocalizedMessage() );
            }
            System.gc();
        } // for( int i = 0; i < number_of_genomes; ++i ) {
        ForesterUtil
                .programMessage( SurfacingConstants.PRG_NAME,
                                 "Wrote domain promiscuities to: " + per_genome_domain_promiscuity_statistics_file );
        final int LEVEL = 0;
        try {
            MinimalDomainomeCalculator.calc( false,
                                             intrees[ 0 ],
                                             LEVEL,
                                             protein_lists_per_species,
                                             SEPARATOR_FOR_DA,
                                             -1,
                                             out_dir.toString() + "/" + output_file,
                                             true,
                                             false,
                                             false,
                                             null );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.getLocalizedMessage() );
        }
        try {
            MinimalDomainomeCalculator.calc( true,
                                             intrees[ 0 ],
                                             LEVEL,
                                             protein_lists_per_species,
                                             SEPARATOR_FOR_DA,
                                             -1,
                                             out_dir.toString() + "/" + output_file,
                                             true,
                                             obtain_names_for_das_from_db,
                                             write_da_ids_names_maps,
                                             input_da_name_file );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.getLocalizedMessage() );
        }
        if ( da_analysis ) {
            SurfacingUtil.performDomainArchitectureAnalysis( distinct_domain_architecutures_per_genome,
                                                             distinct_domain_architecuture_counts,
                                                             10,
                                                             new File( out_dir.toString() + "/" + output_file
                                                                     + "_DA_counts.txt" ),
                                                             new File( out_dir.toString() + "/" + output_file
                                                                     + "_unique_DAs.txt" ) );
            distinct_domain_architecutures_per_genome.clear();
            distinct_domain_architecuture_counts.clear();
            System.gc();
        }
        try {
            domains_per_potein_stats_writer.write( "ALL" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer.write( all_genomes_domains_per_potein_stats.arithmeticMean() + "" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer
                    .write( all_genomes_domains_per_potein_stats.sampleStandardDeviation() + "" );
            domains_per_potein_stats_writer.write( "\t" );
            if ( all_genomes_domains_per_potein_stats.getN() <= 300 ) {
                domains_per_potein_stats_writer.write( all_genomes_domains_per_potein_stats.median() + "" );
                domains_per_potein_stats_writer.write( "\t" );
            }
            domains_per_potein_stats_writer.write( all_genomes_domains_per_potein_stats.getN() + "" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer.write( all_genomes_domains_per_potein_stats.getMin() + "" );
            domains_per_potein_stats_writer.write( "\t" );
            domains_per_potein_stats_writer.write( all_genomes_domains_per_potein_stats.getMax() + "" );
            domains_per_potein_stats_writer.write( "\n" );
            domains_per_potein_stats_writer.close();
            all_genomes_domains_per_potein_stats = null;
            SurfacingUtil.printOutPercentageOfMultidomainProteins( all_genomes_domains_per_potein_histo, log_writer );
            ForesterUtil.map2file(
                                   new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file
                                           + "_all_genomes_domains_per_potein_histo.txt" ),
                                   all_genomes_domains_per_potein_histo,
                                   "\t",
                                   "\n" );
            ForesterUtil.collection2file(
                                          new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file
                                                  + "_domains_always_single_.txt" ),
                                          domains_which_are_always_single,
                                          "\n" );
            ForesterUtil.collection2file(
                                          new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file
                                                  + "_domains_single_or_combined.txt" ),
                                          domains_which_are_sometimes_single_sometimes_not,
                                          "\n" );
            ForesterUtil.collection2file(
                                          new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file
                                                  + "_domains_always_combined.txt" ),
                                          domains_which_never_single,
                                          "\n" );
            ForesterUtil.programMessage( SurfacingConstants.PRG_NAME,
                                         "Average of proteins with a least one domain assigned: "
                                                 + ( 100 * protein_coverage_stats.arithmeticMean() ) + "% (+/-"
                                                 + ( 100 * protein_coverage_stats.sampleStandardDeviation() ) + "%)" );
            ForesterUtil.programMessage( SurfacingConstants.PRG_NAME,
                                         "Range of proteins with a least one domain assigned: "
                                                 + ( 100 * protein_coverage_stats.getMin() ) + "%-"
                                                 + ( 100 * protein_coverage_stats.getMax() ) + "%" );
            SurfacingUtil.log(
                               "Average of prot with a least one dom assigned  : "
                                       + ( 100 * protein_coverage_stats.arithmeticMean() ) + "% (+/-"
                                       + ( 100 * protein_coverage_stats.sampleStandardDeviation() ) + "%)",
                               log_writer );
            SurfacingUtil.log(
                               "Range of prot with a least one dom assigned    : "
                                       + ( 100 * protein_coverage_stats.getMin() ) + "%-"
                                       + ( 100 * protein_coverage_stats.getMax() ) + "%",
                               log_writer );
            protein_coverage_stats = null;
        }
        catch ( final IOException e2 ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e2.getLocalizedMessage() );
        }
        if ( query_domains_writer_ary != null ) {
            for( int j = 0; j < query_domain_ids_array.length; j++ ) {
                try {
                    query_domains_writer_ary[ j ].close();
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.toString() );
                }
            }
        }
        try {
            per_genome_domain_promiscuity_statistics_writer.close();
            dc_data_writer.close();
            log_writer.close();
        }
        catch ( final IOException e2 ) {
            ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e2.getLocalizedMessage() );
        }
        if ( domain_length_analysis ) {
            try {
                SurfacingUtil.executeDomainLengthAnalysis( input_file_properties,
                                                           number_of_genomes,
                                                           domain_lengths_table,
                                                           domain_lengths_analysis_outfile );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e1.toString() );
            }
            System.out.println();
            ForesterUtil.programMessage( SurfacingConstants.PRG_NAME,
                                         "Wrote domain length data to: " + domain_lengths_analysis_outfile );
            System.out.println();
        }
        domain_lengths_table = null;
        final long analysis_start_time = new Date().getTime();
        PairwiseDomainSimilarityCalculator pw_calc = null;
        final DomainSimilarityCalculator calc = new BasicDomainSimilarityCalculator( domain_similarity_sort_field,
                                                                                     sort_by_species_count_first,
                                                                                     number_of_genomes == 2,
                                                                                     CALC_SIMILARITY_SCORES,
                                                                                     true );
        switch ( scoring ) {
            case COMBINATIONS:
                pw_calc = new CombinationsBasedPairwiseDomainSimilarityCalculator();
                break;
            case DOMAINS:
                pw_calc = new DomainCountsBasedPairwiseSimilarityCalculator();
                break;
            case PROTEINS:
                pw_calc = new ProteinCountsBasedPairwiseDomainSimilarityCalculator();
                break;
            default:
                ForesterUtil.unexpectedFatalError( SurfacingConstants.PRG_NAME,
                                                   "unknown value for sorting for scoring" );
        }
        DomainSimilarityCalculator.GoAnnotationOutput go_annotation_output = DomainSimilarityCalculator.GoAnnotationOutput.NONE;
        if ( domain_id_to_go_ids_map != null ) {
            go_annotation_output = DomainSimilarityCalculator.GoAnnotationOutput.ALL;
        }
        final SortedSet<DomainSimilarity> similarities = calc
                .calculateSimilarities( pw_calc,
                                        gwcd_list,
                                        ignore_domains_without_combs_in_all_spec,
                                        ignore_species_specific_domains );
        SurfacingUtil.decoratePrintableDomainSimilarities( similarities, detailedness );
        final Map<String, Integer> tax_code_to_id_map = SurfacingUtil.createTaxCodeToIdMap( intrees[ 0 ] );
        try {
            String my_outfile = output_file.toString();
            Map<Character, Writer> split_writers = null;
            Writer writer = null;
            if ( similarities.size() > MINIMAL_NUMBER_OF_SIMILARITIES_FOR_SPLITTING ) {
                if ( my_outfile.endsWith( ".html" ) ) {
                    my_outfile = my_outfile.substring( 0, my_outfile.length() - 5 );
                }
                split_writers = new HashMap<>();
                SurfacingUtil.createSplitWriters( out_dir, my_outfile, split_writers );
            }
            else if ( !my_outfile.endsWith( ".html" ) ) {
                my_outfile += ".html";
                writer = new BufferedWriter( new FileWriter( out_dir + ForesterUtil.FILE_SEPARATOR + my_outfile ) );
            }
            List<Species> species_order = null;
            if ( species_matrix ) {
                species_order = new ArrayList<>();
                for( int i = 0; i < number_of_genomes; i++ ) {
                    species_order.add( new BasicSpecies( input_file_properties[ i ][ 1 ] ) );
                }
            }
            html_desc.append( "<tr><td>Sum of all distinct binary combinations:</td><td>"
                    + all_bin_domain_combinations_encountered.size() + "</td></tr>" + nl );
            html_desc.append( "<tr><td>Sum of all distinct domains:</td><td>" + all_domains_encountered.size()
                    + "</td></tr>" + nl );
            html_desc.append( "<tr><td>Analysis date/time:</td><td>"
                    + new java.text.SimpleDateFormat( "yyyy.MM.dd HH:mm:ss" ).format( new java.util.Date() )
                    + "</td></tr>" + nl );
            html_desc.append( "</table>" + nl );
            final Writer simple_tab_writer = new BufferedWriter( new FileWriter( out_dir + ForesterUtil.FILE_SEPARATOR
                    + my_outfile.replaceFirst( ".html", ".tsv" ) ) );
            final String domain_species_seqid_map_writer_name = my_outfile.replaceFirst( ".html",
                                                                                         DOMAIN_SPECIES_IDS_MAP_NAME );
            final Writer domain_species_seqid_map_writer = new BufferedWriter( new FileWriter( out_dir
                    + ForesterUtil.FILE_SEPARATOR + domain_species_seqid_map_writer_name ) );
            SurfacingUtil.writeDomainSimilaritiesToFile( html_desc,
                                                         new StringBuilder( number_of_genomes + " genomes" ),
                                                         simple_tab_writer,
                                                         writer,
                                                         domain_species_seqid_map_writer,
                                                         split_writers,
                                                         similarities,
                                                         number_of_genomes == 2,
                                                         species_order,
                                                         domain_similarity_print_option,
                                                         scoring,
                                                         true,
                                                         tax_code_to_id_map,
                                                         intree_0_orig,
                                                         positive_filter_file != null ? filter : null );
            simple_tab_writer.close();
            domain_species_seqid_map_writer.flush();
            domain_species_seqid_map_writer.close();
            ForesterUtil.programMessage( SurfacingConstants.PRG_NAME,
                                         "Wrote domain-species-ids map to       : "
                                                 + domain_species_seqid_map_writer_name );
            ForesterUtil.programMessage( SurfacingConstants.PRG_NAME,
                                         "Wrote main output (includes domain similarities) to: \""
                                                 + ( out_dir == null ? my_outfile
                                                         : out_dir + ForesterUtil.FILE_SEPARATOR + my_outfile )
                                                 + "\"" );
        }
        catch ( final IOException e ) {
            ForesterUtil
                    .fatalError( SurfacingConstants.PRG_NAME,
                                 "Failed to write similarites to: \"" + output_file + "\" [" + e.getMessage() + "]" );
        }
        System.out.println();
        final Species[] species = new Species[ number_of_genomes ];
        for( int i = 0; i < number_of_genomes; ++i ) {
            species[ i ] = new BasicSpecies( input_file_properties[ i ][ 1 ] );
        }
        List<Phylogeny> inferred_trees = null;
        if ( ( number_of_genomes > 2 ) && perform_pwc ) {
            final PairwiseGenomeComparator pwgc = new PairwiseGenomeComparator();
            pwgc.performPairwiseComparisons( html_desc,
                                             sort_by_species_count_first,
                                             detailedness,
                                             ignore_domains_without_combs_in_all_spec,
                                             ignore_species_specific_domains,
                                             domain_similarity_sort_field_for_automated_pwc,
                                             domain_similarity_print_option,
                                             scoring,
                                             domain_id_to_go_ids_map,
                                             go_id_to_term_map,
                                             go_namespace_limit,
                                             species,
                                             number_of_genomes,
                                             gwcd_list,
                                             pw_calc,
                                             automated_pairwise_comparison_suffix,
                                             true,
                                             SurfacingConstants.PAIRWISE_DOMAIN_COMPARISONS_PREFIX,
                                             SurfacingConstants.PRG_NAME,
                                             out_dir,
                                             write_pwc_files,
                                             tax_code_to_id_map,
                                             CALC_SIMILARITY_SCORES,
                                             intree_0_orig );
            String matrix_output_file = new String( output_file.toString() );
            if ( matrix_output_file.indexOf( '.' ) > 1 ) {
                matrix_output_file = matrix_output_file.substring( 0, matrix_output_file.indexOf( '.' ) );
            }
            if ( out_dir != null ) {
                matrix_output_file = out_dir + ForesterUtil.FILE_SEPARATOR + matrix_output_file;
                output_file = new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file );
            }
            SurfacingUtil.writeMatrixToFile(
                                             new File( matrix_output_file
                                                     + launch.MATRIX_MEAN_SCORE_BASED_GENOME_DISTANCE_SUFFIX ),
                                             pwgc.getDomainDistanceScoresMeans() );
            SurfacingUtil.writeMatrixToFile( new File( matrix_output_file
                    + launch.MATRIX_SHARED_BIN_COMBINATIONS_BASED_GENOME_DISTANCE_SUFFIX ),
                                             pwgc.getSharedBinaryCombinationsBasedDistances() );
            SurfacingUtil.writeMatrixToFile(
                                             new File( matrix_output_file
                                                     + launch.MATRIX_SHARED_DOMAINS_BASED_GENOME_DISTANCE_SUFFIX ),
                                             pwgc.getSharedDomainsBasedDistances() );
            final Phylogeny nj_gd = SurfacingUtil
                    .createNjTreeBasedOnMatrixToFile( new File( matrix_output_file
                            + launch.NJ_TREE_MEAN_SCORE_BASED_GENOME_DISTANCE_SUFFIX ),
                                                      pwgc.getDomainDistanceScoresMeans().get( 0 ) );
            final Phylogeny nj_bc = SurfacingUtil
                    .createNjTreeBasedOnMatrixToFile( new File( matrix_output_file
                            + launch.NJ_TREE_SHARED_BIN_COMBINATIONS_BASED_GENOME_DISTANCE_SUFFIX ),
                                                      pwgc.getSharedBinaryCombinationsBasedDistances().get( 0 ) );
            final Phylogeny nj_d = SurfacingUtil
                    .createNjTreeBasedOnMatrixToFile( new File( matrix_output_file
                            + launch.NJ_TREE_SHARED_DOMAINS_BASED_GENOME_DISTANCE_SUFFIX ),
                                                      pwgc.getSharedDomainsBasedDistances().get( 0 ) );
            inferred_trees = new ArrayList<>();
            inferred_trees.add( nj_gd );
            inferred_trees.add( nj_bc );
            inferred_trees.add( nj_d );
        } // if ( ( output_file != null ) && ( number_of_genomes > 2 ) && !isEmpty( automated_pairwise_comparison_suffix ) )
        if ( ( out_dir != null ) && ( !perform_pwc ) ) {
            output_file = new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file );
        }
        if ( write_to_nexus ) {
            SurfacingUtil.writePresentToNexus( output_file, positive_filter_file, filter, gwcd_list );
        }
        if ( ( ( intrees != null ) && ( intrees.length > 0 ) ) && ( number_of_genomes > 2 ) ) {
            final StringBuilder parameters_sb = SurfacingUtil.createParametersAsString( ignore_dufs,
                                                                                        ie_value_max,
                                                                                        fs_e_value_max,
                                                                                        max_allowed_overlap,
                                                                                        no_engulfing_overlaps,
                                                                                        cutoff_scores_file,
                                                                                        dc_type );
            String s = "_";
            if ( radomize_fitch_parsimony ) {
                s += random_number_seed_for_fitch_parsimony + "_";
            }
            int i = 0;
            for( final Phylogeny intree : intrees ) {
                final String outfile_name = ForesterUtil.removeSuffix( output_file.toString() ) + s
                        + ForesterUtil.removeSuffix( intree_files[ i ].toString() );
                final DomainParsimonyCalculator domain_parsimony = DomainParsimonyCalculator
                        .createInstance( intree, gwcd_list );
                SurfacingUtil.executeParsimonyAnalysis( random_number_seed_for_fitch_parsimony,
                                                        radomize_fitch_parsimony,
                                                        outfile_name,
                                                        domain_parsimony,
                                                        intree,
                                                        domain_id_to_go_ids_map,
                                                        go_id_to_term_map,
                                                        go_namespace_limit,
                                                        parameters_sb.toString(),
                                                        domain_id_to_secondary_features_maps,
                                                        positive_filter_file == null ? null : filter,
                                                        output_binary_domain_combinationsfor_graph_analysis,
                                                        all_bin_domain_combinations_gained_fitch,
                                                        all_bin_domain_combinations_lost_fitch,
                                                        dc_type,
                                                        protein_length_stats_by_dc,
                                                        domain_number_stats_by_dc,
                                                        domain_length_stats_by_domain,
                                                        tax_code_to_id_map,
                                                        write_to_nexus,
                                                        use_last_in_fitch_parsimony,
                                                        perform_dc_fich );
                // Listing of all domain combinations gained is only done if only one input tree is used.
                if ( ( domain_id_to_secondary_features_maps != null )
                        && ( domain_id_to_secondary_features_maps.length > 0 ) ) {
                    int j = 0;
                    for( final Map<String, Set<String>> domain_id_to_secondary_features_map : domain_id_to_secondary_features_maps ) {
                        final Map<Species, MappingResults> mapping_results_map = new TreeMap<>();
                        final DomainParsimonyCalculator secondary_features_parsimony = DomainParsimonyCalculator
                                .createInstance( intree, gwcd_list, domain_id_to_secondary_features_map );
                        SurfacingUtil.executeParsimonyAnalysisForSecondaryFeatures( outfile_name + "_"
                                + secondary_features_map_files[ j++ ],
                                                                                    secondary_features_parsimony,
                                                                                    intree,
                                                                                    parameters_sb.toString(),
                                                                                    mapping_results_map,
                                                                                    use_last_in_fitch_parsimony );
                        if ( i == 0 ) {
                            System.out.println();
                            System.out.println( "Mapping to secondary features:" );
                            for( final Species spec : mapping_results_map.keySet() ) {
                                final MappingResults mapping_results = mapping_results_map.get( spec );
                                final int total_domains = mapping_results.getSumOfFailures()
                                        + mapping_results.getSumOfSuccesses();
                                System.out.print( spec + ":" );
                                System.out.print( " mapped domains = " + mapping_results.getSumOfSuccesses() );
                                System.out.print( ", not mapped domains = " + mapping_results.getSumOfFailures() );
                                if ( total_domains > 0 ) {
                                    System.out.println( ", mapped ratio = "
                                            + ( ( 100 * mapping_results.getSumOfSuccesses() ) / total_domains ) + "%" );
                                }
                                else {
                                    System.out.println( ", mapped ratio = n/a (total domains = 0 )" );
                                }
                            }
                        }
                    }
                }
                i++;
            } // for( final Phylogeny intree : intrees ) {
        }
        if ( plus_minus_analysis_high_copy_base_species.size() > 0 ) {
            SurfacingUtil.executePlusMinusAnalysis( output_file,
                                                    plus_minus_analysis_high_copy_base_species,
                                                    plus_minus_analysis_high_copy_target_species,
                                                    plus_minus_analysis_high_low_copy_species,
                                                    gwcd_list,
                                                    protein_lists_per_species,
                                                    domain_id_to_go_ids_map,
                                                    go_id_to_term_map,
                                                    plus_minus_analysis_numbers );
        }
        if ( output_protein_lists_for_all_domains ) {
            SurfacingUtil.writeProteinListsForAllSpecies( out_dir,
                                                          protein_lists_per_species,
                                                          gwcd_list,
                                                          output_list_of_all_proteins_per_domain_e_value_max,
                                                          positive_filter_file != null ? filter : null );
        }
        gwcd_list = null;
        if ( all_bin_domain_combinations_gained_fitch != null ) {
            try {
                SurfacingUtil.executeFitchGainsAnalysis( new File( output_file
                        + launch.OUTPUT_DOMAIN_COMBINATIONS_GAINED_MORE_THAN_ONCE_ANALYSIS_SUFFIX ),
                                                         all_bin_domain_combinations_gained_fitch,
                                                         all_domains_encountered.size(),
                                                         all_bin_domain_combinations_encountered,
                                                         true );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.getLocalizedMessage() );
            }
        }
        if ( all_bin_domain_combinations_lost_fitch != null ) {
            try {
                SurfacingUtil.executeFitchGainsAnalysis( new File( output_file
                        + launch.OUTPUT_DOMAIN_COMBINATIONS_LOST_MORE_THAN_ONCE_ANALYSIS_SUFFIX ),
                                                         all_bin_domain_combinations_lost_fitch,
                                                         all_domains_encountered.size(),
                                                         all_bin_domain_combinations_encountered,
                                                         false );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( SurfacingConstants.PRG_NAME, e.getLocalizedMessage() );
            }
        }
        final Runtime rt = java.lang.Runtime.getRuntime();
        final long free_memory = rt.freeMemory() / 1000000;
        final long total_memory = rt.totalMemory() / 1000000;
        ForesterUtil.programMessage( SurfacingConstants.PRG_NAME,
                                     "Time for analysis : " + ( new Date().getTime() - analysis_start_time ) + "ms" );
        ForesterUtil.programMessage( SurfacingConstants.PRG_NAME,
                                     "Total running time: " + ( new Date().getTime() - start_time ) + "ms " );
        ForesterUtil
                .programMessage( SurfacingConstants.PRG_NAME,
                                 "Free memory       : " + free_memory + "MB, total memory: " + total_memory + "MB" );
        ForesterUtil.programMessage( SurfacingConstants.PRG_NAME,
                                     "If this application is useful to you, please cite:" );
        ForesterUtil.programMessage( SurfacingConstants.PRG_NAME, SurfacingConstants.WWW );
        ForesterUtil
                .programMessage( SurfacingConstants.PRG_NAME,
                                 "[next step for phylogenomic analysis pipeline (example, in \"DAS\" dir): % mse.rb .prot . FL_seqs DA_seqs path/to/genome_locations.txt]" );
        ForesterUtil.programMessage( SurfacingConstants.PRG_NAME, "OK" );
        System.out.println();
    }

    private static void printHelp() {
        System.out.println();
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( "% java -Xms256m -Xmx512m -cp forester.jar org.forester.applications."
                + SurfacingConstants.PRG_NAME + " [options] <phylogen(y|ies) infile>" );
        System.out.println();
        System.out.println( " Note: This software might need a significant amount of memory (heap space);" );
        System.out
                .println( "       hence use \"-Xms128m -Xmx512m\" (or more) to prevent a \"java.lang.OutOfMemoryError\"." );
        System.out.println();
        System.out.println( " Options: " );
        System.out.println( launch.DETAILEDNESS_OPTION + ": level of detail for similarities output file (default:"
                + DETAILEDNESS_DEFAULT + ")" );
        System.out.println( launch.IGNORE_COMBINATION_WITH_SAME_OPTION
                + ": to ignore combinations with self (default: not to ignore)" );
        System.out.println( launch.IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_OPTION
                + ": to ignore domains without combinations in any species (for similarity calc purposes, not for parsimony analyses) (default: not to ignore)" );
        System.out.println( launch.IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION
                + ": to ignore domains specific to one species (for similarity calc purposes, not for parsimony analyses) (default: not to ignore)" );
        System.out.println( launch.NOT_IGNORE_DUFS_OPTION
                + ": to _not_ ignore DUFs (domains with unknown function) (default: ignore DUFs)" );
        System.out.println( launch.IGNORE_VIRAL_IDS
                + ": to ignore domains with ids containing 'vir', 'retro', 'transpos', 'phage', or starting with 'rv' or 'gag_'" );
        System.out.println( launch.DOMAIN_SIMILARITY_SORT_OPTION + ": sorting for similarities (default: "
                + DOMAIN_SORT_FILD_DEFAULT + ")" );
        System.out.println( launch.OUTPUT_FILE_OPTION + ": name for (main) output file (mandatory)" );
        System.out.println( launch.MAX_I_E_VALUE_OPTION + ": max (inclusive) iE-value" );
        System.out.println( launch.MAX_FS_E_VALUE_OPTION + ": max (inclusive) FS E-value" );
        System.out
                .println( launch.MIN_REL_ENV_LENGTH_RATIO_OPTION + ": min (inclusive) relative envelope length ratio" );
        System.out.println( SurfacingConstants.MAX_ALLOWED_OVERLAP_OPTION + ": maximal allowed domain overlap" );
        System.out.println( launch.NO_ENGULFING_OVERLAP_OPTION + ": to ignore engulfed lower confidence domains" );
        System.out.println( launch.SPECIES_MATRIX_OPTION + ": species matrix" );
        System.out.println( launch.SCORING_OPTION + ": scoring (default:" + SCORING_DEFAULT + ")" );
        System.out.println( launch.DOMAIN_COUNT_SORT_OPTION + ": sorting for domain counts (default:"
                + DOMAINS_SORT_ORDER_DEFAULT + ")" );
        System.out.println( launch.DOMAIN_SIMILARITY_PRINT_OPTION + ": domain similarity print option (default:"
                + DOMAIN_SIMILARITY_PRINT_OPTION_DEFAULT + ")" );
        System.out.println( launch.CUTOFF_SCORE_FILE_OPTION + ": cutoff score file" );
        System.out.println( launch.DOMAIN_SIMILARITY_SORT_BY_SPECIES_COUNT_FIRST_OPTION
                + ": sort by species count first" );
        System.out.println( launch.OUTPUT_DIR_OPTION + ": output directory" );
        System.out.println( launch.PFAM_TO_GO_FILE_USE_OPTION + ": Pfam to GO mapping file" );
        System.out.println( launch.GO_OBO_FILE_USE_OPTION + ": GO terms file (OBO format)" );
        System.out.println( launch.GO_NAMESPACE_LIMIT_OPTION + ": limit GO term to one GO namespace" );
        System.out.println( launch.PAIRWISE_DOMAIN_COMPARISONS_OPTION
                + "[=<suffix for pairwise comparison output files>]: to perform pairwise comparison based analyses" );
        System.out.println( launch.INPUT_SPECIES_TREE_OPTION
                + ": species tree, to perform (Dollo, Fitch) parismony analyses" );
        System.out.println( launch.INPUT_SPECIES_TREE_OPTION
                + "=<treefiles in phyloXML format, separated by #>: to infer domain/binary domain combination gains/losses on given species trees" );
        System.out.println( launch.FILTER_POSITIVE_OPTION
                + "=<file>: to filter out proteins not containing at least one domain listed in <file>" );
        System.out.println( launch.FILTER_NEGATIVE_OPTION
                + "=<file>: to filter out proteins containing at least one domain listed in <file>" );
        System.out.println( launch.FILTER_NEGATIVE_DOMAINS_OPTION
                + "=<file>: to filter out (ignore) domains listed in <file>" );
        System.out.println( launch.INPUT_GENOMES_FILE_OPTION + "=<file>: to read input files from <file>" );
        System.out.println( launch.RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION
                + "=<seed>: seed for random number generator for Fitch Parsimony analysis (type: long, default: no randomization - given a choice, prefer absence" );
        System.out.println( launch.CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS
                + ": to consider directedness in binary combinations: e.g. A-B != B-A" );
        System.out.println( launch.CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS_AND_ADJACENCY
                + ": to consider directedness and adjacency in binary combinations" );
        System.out.println( launch.SEQ_EXTRACT_OPTION
                + "=<domain ids (Pfam names)>: to extract sequence names of sequences containing matching domains and/or domain-sequences (order N to C) (domain separator: '~', domain sequences speparator: '#', e.g. 'NACHT#BIR~CARD')" );
        System.out.println( launch.SECONDARY_FEATURES_PARSIMONY_MAP_FILE
                + "=<file>: to perfom parsimony analysis on secondary features" );
        System.out.println( SurfacingConstants.PLUS_MINUS_ANALYSIS_OPTION
                + "=<file>: to presence/absence genome analysis" );
        System.out.println( launch.DOMAIN_COMBINITONS_COUNTS_OUTPUT_OPTION
                + ": to output binary domain counts (as individual files)" );
        System.out.println( launch.DOMAIN_COMBINITONS_OUTPUT_OPTION_FOR_GRAPH_ANALYSIS
                + ": to output binary domain combinations for (downstream) graph analysis" );
        System.out.println( launch.OUTPUT_LIST_OF_ALL_PROTEINS_OPTIONS + ": to output all proteins per domain" );
        System.out.println( launch.OUTPUT_LIST_OF_ALL_PROTEINS_PER_DOMAIN_E_VALUE_OPTION
                + ": e value max per domain for output of all proteins per domain" );
        System.out.println( launch.USE_LAST_IN_FITCH_OPTION + ": to use last in Fitch parsimony" );
        System.out.println( launch.WRITE_TO_NEXUS_OPTION + ": to output in Nexus format" );
        System.out.println( PERFORM_DC_FITCH + ": to perform DC Fitch parsimony" );
        System.out.println( PERFORM_DC_REGAIN_PROTEINS_STATS_OPTION + ": to perform DC regain protein statistics" );
        System.out.println( DA_ANALYSIS_OPTION + ": to perform DA analysis" );
        System.out.println( PERFORM_DOMAIN_LENGTH_ANALYSIS_OPTION + ": to perform domain length analysis" );
        System.out.println( WRITE_DA_IDS_NAMES_MAPS_OPTION + ": to write DA-name-seq IDs mapping files" );
        System.out.println( INPUT_DA_NAME_FILE_OPTION + "=<file>: file to obtain DA names from" );
        System.out.println( OBTAIN_NAMES_FOR_DAS_FROM_DB_OPTION
                + ": to obtain DA names from online database (UniProtKB)" );
        System.out.println( OBTAIN_NAMES_FOR_DAS_FROM_DB_MAX_IDS_TO_SEARCH_PER_SPECIES_OPTION
                + "=<n>: max IDs per species for DB search when obtaining DA names from online database (UniProtKB)" );
        System.out.println( UNIPROT_PRIORITY_FOR_ACCESSOR_PARSING_OPTION
                + ": to use UniProt priority for DB accessor parsing (expert option)" );
        System.out.println( VERBOSITY_OPTION + "=<n>: verbosity" );
        System.out.println();
        System.out.println();
        System.out
                .println( "Example 1: surfacing -p2g=pfam2go.txt -obo=go.obo -species_tree=tol_156.xml -no_eo -ie=0.01 -dufs -genomes=genomes_all.txt -pos_filter=tf_1.txt -out_dir=_tf1 -o=tf1" );
        System.out.println();
        System.out
                .println( "Example 2: surfacing -p2g=pfam2go.txt -obo=go.obo -species_tree=tol_156.xml -last -ignore_viral_ids -no_eo -ie=0.1 -dufs -genomes=genomes_all.txt -pos_filter=tf_1.txt -all_prot -all_prot_e=0.1 -out_dir=_tf1_e01_ape01 -o=tf1_e01_ape01" );
        System.out.println();
        System.out
                .println( "Example 3: surfacing -species_tree=master_tree.xml -no_eo -ie=1e-6 -mrel=0.5 -mo=10 -dufs -genomes=genomes.txt -out_dir=a605 -o=a605" );
        System.out.println();
        System.out.println( "Example 4: surfacing -species_tree=Nidovirales_2.xml -no_eo -ie=1e-6 -mrel=0.4 -mo=5 -dufs -genomes=genomes_coronaviridae1.txt -out_dir=_y6 -o=y6 -write_DA_maps -obtain_DA_names_from_db -max_ids_to_search_per_species=20 -verbosity=3" );
        System.out.println();
        System.out
                .println( "[next step for phylogenomic analysis pipeline (example, in \"DAS\" dir): % mse.rb .prot . FL_seqs DA_seqs path/to/genome_locations.txt]" );
        System.out.println();
    }
}
