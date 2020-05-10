// $Id:
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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

import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

final class SurfacingConstants {

    static final String  PRG_NAME                                                                      = "surfacing";
    static final String  PRG_VERSION                                                                   = "2.700";
    static final String  PRG_DATE                                                                      = "200509";
    static final String  E_MAIL                                                                        = "phyloxml@gmail.com";
    static final String  WWW                                                                           = "https://sites.google.com/site/cmzmasek/home/software/forester/surfacing";
    static final String  AMIGO_LINK                                                                    = "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&query=";
    static final String  NL                                                                            = ForesterUtil.LINE_SEPARATOR;
    static final String  NONE                                                                          = "[none]";
    static final boolean PRINT_MORE_DOM_SIMILARITY_INFO                                                = false;
    static final boolean SECONDARY_FEATURES_ARE_SCOP                                                   = true;
    static final String  SECONDARY_FEATURES_SCOP_LINK                                                  = "http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?key=";
    static final String  PLUS_MINUS_DOM_SUFFIX                                                         = "_plus_minus_dom.txt";
    static final String  PLUS_MINUS_DOM_SUFFIX_HTML                                                    = "_plus_minus_dom.html";
    static final String  PLUS_MINUS_DC_SUFFIX_HTML                                                     = "_plus_minus_dc.html";
    static final int     PLUS_MINUS_ANALYSIS_MIN_DIFF_DEFAULT                                          = 0;
    static final double  PLUS_MINUS_ANALYSIS_FACTOR_DEFAULT                                            = 1.0;
    static final String  PLUS_MINUS_ALL_GO_IDS_DOM_SUFFIX                                              = "_plus_minus_go_ids_all.txt";
    static final String  PLUS_MINUS_PASSING_GO_IDS_DOM_SUFFIX                                          = "_plus_minus_go_ids_passing.txt";
    static final String  ALL_PFAMS_ENCOUNTERED_SUFFIX                                                  = "_all_encountered_pfams";
    static final String  ALL_PFAMS_ENCOUNTERED_WITH_GO_ANNOTATION_SUFFIX                               = "_all_encountered_pfams_with_go_annotation";
    static final String  ENCOUNTERED_PFAMS_SUMMARY_SUFFIX                                              = "_encountered_pfams_summary";
    static final String  ALL_PFAMS_GAINED_AS_DOMAINS_SUFFIX                                            = "_all_pfams_gained_as_domains";
    static final String  ALL_PFAMS_LOST_AS_DOMAINS_SUFFIX                                              = "_all_pfams_lost_as_domains";
    static final String  ALL_PFAMS_GAINED_AS_DC_SUFFIX                                                 = "_all_pfams_gained_as_dc";
    static final String  ALL_PFAMS_LOST_AS_DC_SUFFIX                                                   = "_all_pfams_lost_as_dc";
    static final String  BASE_DIRECTORY_PER_NODE_DOMAIN_GAIN_LOSS_FILES                                = "PER_NODE_EVENTS";
    static final String  BASE_DIRECTORY_PER_SUBTREE_DOMAIN_GAIN_LOSS_FILES                             = "PER_SUBTREE_EVENTS";
    static final String  INDEPENDENT_DC_GAINS_FITCH_PARS_COUNTS_OUTPUT_SUFFIX                          = "_indep_dc_gains_fitch_counts.txt";
    static final String  INDEPENDENT_DC_GAINS_FITCH_PARS_DC_OUTPUT_SUFFIX                              = "_indep_dc_gains_fitch_lists.txt";
    static final String  INDEPENDENT_DC_GAINS_FITCH_PARS_DC_FOR_GO_MAPPING_OUTPUT_SUFFIX               = "_indep_dc_gains_fitch_lists_for_go_mapping.txt";
    static final String  INDEPENDENT_DC_GAINS_FITCH_PARS_DC_FOR_GO_MAPPING_OUTPUT_UNIQUE_SUFFIX        = "_indep_dc_gains_fitch_lists_for_go_mapping_unique.txt";
    static final String  BINARY_DOMAIN_COMBINATIONS_PARSIMONY_TREE_OUTPUT_SUFFIX_FITCH_MAPPED          = "_dc_MAPPED_secondary_features_fitch"
            + ForesterConstants.PHYLO_XML_SUFFIX;
    static final String  INDEPENDENT_DC_GAINS_FITCH_PARS_COUNTS_MAPPED_OUTPUT_SUFFIX                   = "_indep_dc_gains_fitch_counts_MAPPED.txt";
    static final String  INDEPENDENT_DC_GAINS_FITCH_PARS_DC_MAPPED_OUTPUT_SUFFIX                       = "_indep_dc_gains_fitch_lists_MAPPED.txt";
    static final String  INDEPENDENT_DC_GAINS_FITCH_PARS_DC_FOR_GO_MAPPING_MAPPED_OUTPUT_SUFFIX        = "_indep_dc_gains_fitch_lists_for_go_mapping_MAPPED.txt";
    static final String  INDEPENDENT_DC_GAINS_FITCH_PARS_DC_FOR_GO_MAPPING_MAPPED_OUTPUT_UNIQUE_SUFFIX = "_indep_dc_gains_fitch_lists_for_go_mapping_unique_MAPPED.txt";
    static final String  DA_SPECIES_IDS_MAP_NAME                                                       = "_DA_SPECIES_IDS_MAP.txt";
    static final String  SUFFIX_DA_NAME_NAME                                                           = "_SUFFIX_DA_NAME_MAP.txt";
    static final String  DA_NAME_MAP_NAME                                                              = "_DA_NAME_MAP.txt";
    static final String  DOMAIN_COMBINITONS_OUTPUTFILE_SUFFIX_FOR_GRAPH_ANALYSIS                       = "_dc.dot";
    static final String  PARSIMONY_OUTPUT_FITCH_PRESENT_BC_OUTPUTFILE_SUFFIX_FOR_GRAPH_ANALYSIS        = "_fitch_present_dc.dot";
    static final String  PARSIMONY_OUTPUT_GL_SUFFIX_DOLLO_DOMAINS                                      = "_dollo_gl_d";
    static final String  PARSIMONY_OUTPUT_GL_SUFFIX_DOLLO_BINARY_COMBINATIONS                          = "_dollo_gl_dc";
    static final String  PARSIMONY_OUTPUT_GL_SUFFIX_FITCH_DOMAINS                                      = "_fitch_gl_d";
    static final String  PARSIMONY_OUTPUT_GL_SUFFIX_FITCH_BINARY_COMBINATIONS                          = "_fitch_gl_dc";
    // gain/loss counts:
    static final String  PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_DOLLO_DOMAINS                               = "_dollo_glc_d";
    static final String  PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_DOLLO_BINARY_COMBINATIONS                   = "_dollo_glc_dc";
    static final String  PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_FITCH_DOMAINS                               = "_fitch_glc_d";
    static final String  PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_FITCH_BINARY_COMBINATIONS                   = "_fitch_glc_dc";
    // tables:
    static final String  PARSIMONY_OUTPUT_FITCH_GAINS_BC                                               = "_fitch_gains_dc";
    static final String  PARSIMONY_OUTPUT_FITCH_GAINS_HTML_BC                                          = "_fitch_gains_dc.html";
    static final String  PARSIMONY_OUTPUT_FITCH_LOSSES_BC                                              = "_fitch_losses_dc";
    static final String  PARSIMONY_OUTPUT_FITCH_LOSSES_HTML_BC                                         = "_fitch_losses_dc.html";
    static final String  PARSIMONY_OUTPUT_FITCH_PRESENT_BC                                             = "_fitch_present_dc";
    static final String  PARSIMONY_OUTPUT_FITCH_PRESENT_HTML_BC                                        = "_fitch_present_dc.html";
    static final String  PARSIMONY_OUTPUT_DOLLO_GAINS_D                                                = "_dollo_gains_d";
    static final String  PARSIMONY_OUTPUT_DOLLO_GAINS_HTML_D                                           = "_dollo_gains_d.html";
    static final String  PARSIMONY_OUTPUT_DOLLO_LOSSES_D                                               = "_dollo_losses_d";
    static final String  PARSIMONY_OUTPUT_DOLLO_LOSSES_HTML_D                                          = "_dollo_losses_d.html";
    static final String  PARSIMONY_OUTPUT_DOLLO_PRESENT_D                                              = "_dollo_present_d";
    static final String  PARSIMONY_OUTPUT_DOLLO_PRESENT_HTML_D                                         = "_dollo_present_d.html";
    static final String  DOMAINS_PRESENT_NEXUS                                                         = "_dom.nex";
    static final String  BDC_PRESENT_NEXUS                                                             = "_dc.nex";
    static final String  DOMAINS_PARSIMONY_TREE_OUTPUT_SUFFIX_DOLLO                                    = "_d_dollo"
            + ForesterConstants.PHYLO_XML_SUFFIX;
    static final String  DOMAINS_PARSIMONY_TREE_OUTPUT_SUFFIX_FITCH                                    = "_d_fitch"
            + ForesterConstants.PHYLO_XML_SUFFIX;
    static final String  BINARY_DOMAIN_COMBINATIONS_PARSIMONY_TREE_OUTPUT_SUFFIX_DOLLO                 = "_dc_dollo"
            + ForesterConstants.PHYLO_XML_SUFFIX;
    static final String  BINARY_DOMAIN_COMBINATIONS_PARSIMONY_TREE_OUTPUT_SUFFIX_FITCH                 = "_dc_fitch"
            + ForesterConstants.PHYLO_XML_SUFFIX;
    static final String  NEXUS_EXTERNAL_DOMAINS                                                        = "_dom.nex";
    static final String  NEXUS_EXTERNAL_DOMAIN_COMBINATIONS                                            = "_dc.nex";
    static final String  NEXUS_SECONDARY_FEATURES                                                      = "_secondary_features.nex";
    static final String  PARSIMONY_OUTPUT_GL_SUFFIX_DOLLO_SECONDARY_FEATURES                           = "_dollo_gl_secondary_features";
    static final String  PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_DOLLO_SECONDARY_FEATURES                    = "_dollo_glc_secondary_features";
    static final String  PARSIMONY_OUTPUT_DOLLO_GAINS_SECONDARY_FEATURES                               = "_dollo_gains_secondary_features";
    static final String  PARSIMONY_OUTPUT_DOLLO_LOSSES_SECONDARY_FEATURES                              = "_dollo_losses_secondary_features";
    static final String  PARSIMONY_OUTPUT_DOLLO_PRESENT_SECONDARY_FEATURES                             = "_dollo_present_secondary_features";
    static final String  SECONDARY_FEATURES_PARSIMONY_TREE_OUTPUT_SUFFIX_DOLLO                         = "_secondary_features_dollo"
            + ForesterConstants.PHYLO_XML_SUFFIX;
    static final String  PARSIMONY_OUTPUT_DOLLO_ALL_GOID_D_ALL_NAMESPACES                              = "_dollo_goid_d";
    static final String  PARSIMONY_OUTPUT_FITCH_ALL_GOID_BC_ALL_NAMESPACES                             = "_fitch_goid_dc";
    static final String  LIMIT_SPEC_FOR_PROT_EX                                                        = null;                                                                                            // e.g. "HUMAN"; set to null for not using this feature (default).
    static final String  SEQ_EXTRACT_SUFFIX                                                            = ".prot";
    static final String  PAIRWISE_DOMAIN_COMPARISONS_PREFIX                                            = "pwc_";
    static final String  DOMAIN_COMBINITON_COUNTS_OUTPUTFILE_SUFFIX                                    = ".dcc";
    static final String  MAX_ALLOWED_OVERLAP_OPTION                                                    = "mo";
    static final int     MAX_ALLOWED_OVERLAP_DEFAULT                                                   = -1;
    static final String  PLUS_MINUS_ANALYSIS_OPTION                                                    = "plus_minus";
}
