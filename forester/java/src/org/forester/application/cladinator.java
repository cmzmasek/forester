// $Id:
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
// Contact: phyloxml @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.forester.clade_analysis.AnalysisMulti;
import org.forester.clade_analysis.Prefix;
import org.forester.clade_analysis.ResultMulti;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.EasyWriter;
import org.forester.util.ForesterUtil;
import org.forester.util.UserException;

public final class cladinator {

    final static private String PRG_NAME = "cladinator";
    final static private String PRG_VERSION = "1.1.1.a1";
    final static private String PRG_DATE = "2025-07-24";
    final static private String PRG_DESC = "clades within clades of annotated labels -- analysis of pplacer-type outputs";
    final static private String E_MAIL = "czmasek AT jcvi DOT org";
    final static private String WWW = "https://github.com/cmzmasek";
    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";
    final static private String SEP_OPTION = "s";
    final static private String QUERY_PATTERN_OPTION = "q";
    final static private String MAPPING_FILE_OPTION = "m";
    final static private String EXTRA_PROCESSING_OPTION1 = "x";
    final static private String EXTRA_PROCESSING1_SEP_OPTION = "xs";
    final static private String EXTRA_PROCESSING1_KEEP_EXTRA_OPTION = "xk";
    final static private String QUIET_OPTION = "Q";
    final static private String SPECIAL_PROCESSING_OPTION = "S";
    final static private String VERBOSE_OPTION = "v";
    final static private String REMOVE_ANNOT_SEP_OPTION = "rs";
    final static private String SIMPLE_OUTPUT_OPTION = "sim";
    final static private String SIMPLE_OUTPUT_CUTOFF_OPTION = "simc";
    final static private double SPECIFICS_CUTOFF_DEFAULT = 0.7;
    final static private String SEP_DEFAULT = ".";
    final static private Pattern QUERY_PATTERN_DEFAULT = AnalysisMulti.DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE;
    final static private String EXTRA_PROCESSING1_SEP_DEFAULT = "|";
    final static private boolean EXTRA_PROCESSING1_KEEP_EXTRA_DEFAULT = false;
    private final static DecimalFormat df = new DecimalFormat("0.0###");
    final static private boolean TWO_COLUMNS_IN_SIMPLE_OUTPUT = true;
    final static private boolean ALLOW_TO_SPLIT_QUERY = true;

    public static void main(final String[] args) {
        try {
            ForesterUtil.printProgramInformation(PRG_NAME,
                    PRG_DESC,
                    PRG_VERSION,
                    PRG_DATE,
                    E_MAIL,
                    WWW,
                    ForesterUtil.getForesterLibraryInformation());
            CommandLineArguments cla = null;
            try {
                cla = new CommandLineArguments(args);
            } catch (final Exception e) {
                ForesterUtil.fatalError(PRG_NAME, e.getMessage());
            }
            if (cla.isOptionSet(HELP_OPTION_1) || cla.isOptionSet(HELP_OPTION_2)) {
                System.out.println();
                print_help();
                System.exit(0);
            }
            if ((cla.getNumberOfNames() != 1) && (cla.getNumberOfNames() != 2)) {
                print_help();
                System.exit(-1);
            }
            final List<String> allowed_options = new ArrayList<>();
            allowed_options.add(SEP_OPTION);
            allowed_options.add(QUERY_PATTERN_OPTION);
            allowed_options.add(MAPPING_FILE_OPTION);
            allowed_options.add(EXTRA_PROCESSING_OPTION1);
            allowed_options.add(EXTRA_PROCESSING1_SEP_OPTION);
            allowed_options.add(EXTRA_PROCESSING1_KEEP_EXTRA_OPTION);
            allowed_options.add(SPECIAL_PROCESSING_OPTION);
            allowed_options.add(VERBOSE_OPTION);
            allowed_options.add(QUIET_OPTION);
            allowed_options.add(REMOVE_ANNOT_SEP_OPTION);
            allowed_options.add(SIMPLE_OUTPUT_OPTION);
            allowed_options.add(SIMPLE_OUTPUT_CUTOFF_OPTION);
            final String dissallowed_options = cla.validateAllowedOptionsAsString(allowed_options);
            if (dissallowed_options.length() > 0) {
                ForesterUtil.fatalError(PRG_NAME, "unknown option(s): " + dissallowed_options);
            }

            String separator = SEP_DEFAULT;
            if (cla.isOptionSet(SEP_OPTION)) {
                if (cla.isOptionValueSet(SEP_OPTION)) {
                    separator = cla.getOptionValue(SEP_OPTION);
                } else {
                    ForesterUtil.fatalError(PRG_NAME, "no value for separator option");
                }
            }
            Pattern compiled_query = null;
            if (cla.isOptionSet(QUERY_PATTERN_OPTION)) {
                if (cla.isOptionValueSet(QUERY_PATTERN_OPTION)) {
                    final String query_str = cla.getOptionValue(QUERY_PATTERN_OPTION);
                    try {
                        compiled_query = Pattern.compile(query_str);
                    } catch (final PatternSyntaxException e) {
                        ForesterUtil.fatalError(PRG_NAME,
                                "error in regular expression: " + query_str + ": " + e.getMessage());
                    }
                } else {
                    ForesterUtil.fatalError(PRG_NAME, "no value for query pattern option");
                }
            }
            File mapping_file = null;
            if (cla.isOptionSet(MAPPING_FILE_OPTION)) {
                if (cla.isOptionValueSet(MAPPING_FILE_OPTION)) {
                    final String mapping_file_str = cla.getOptionValue(MAPPING_FILE_OPTION);
                    final String error = ForesterUtil.isReadableFile(mapping_file_str);
                    if (!ForesterUtil.isEmpty(error)) {
                        ForesterUtil.fatalError(PRG_NAME, error);
                    }
                    mapping_file = new File(mapping_file_str);
                } else {
                    ForesterUtil.fatalError(PRG_NAME, "no value for mapping file");
                }
            }
            final Pattern pattern = (compiled_query != null) ? compiled_query : QUERY_PATTERN_DEFAULT;
            final File intreefile = cla.getFile(0);
            final String error_intreefile = ForesterUtil.isReadableFile(intreefile);
            if (!ForesterUtil.isEmpty(error_intreefile)) {
                ForesterUtil.fatalError(PRG_NAME, error_intreefile);
            }
            final File outtablefile;
            if (cla.getNumberOfNames() > 1) {
                outtablefile = cla.getFile(1);
                final String error_outtablefile = ForesterUtil.isWritableFile(outtablefile);
                if (!ForesterUtil.isEmpty(error_outtablefile)) {
                    ForesterUtil.fatalError(PRG_NAME, error_outtablefile);
                }
            } else {
                outtablefile = null;
            }
            final boolean simple_output;
            double simple_output_cutoff = 0.0;
            simple_output = cla.isOptionSet(SIMPLE_OUTPUT_OPTION);
            if (cla.isOptionSet(SIMPLE_OUTPUT_CUTOFF_OPTION)) {
                if (!simple_output) {
                    ForesterUtil.fatalError(PRG_NAME,
                            "must set simple output option to use cutoff for simple output");
                }
                if (cla.isOptionValueSet(SIMPLE_OUTPUT_CUTOFF_OPTION)) {
                    simple_output_cutoff = cla.getOptionValueAsDouble(SIMPLE_OUTPUT_CUTOFF_OPTION);
                    if ((simple_output_cutoff < 0) || (simple_output_cutoff > 1)) {
                        ForesterUtil.fatalError(PRG_NAME, "cutoff must be between 0.0 and 1.0");
                    }
                } else {
                    ForesterUtil.fatalError(PRG_NAME, "no value for cutoff");
                }
            }
            final BasicTable<String> t;
            final SortedMap<String, String> map;
            if (mapping_file != null) {
                t = BasicTableParser.parse(mapping_file, '\t');
                if (t.getNumberOfColumns() != 2) {
                    ForesterUtil.fatalError(PRG_NAME,
                            "mapping file needs to have 2 tab-separated columns, not "
                                    + t.getNumberOfColumns());
                }
                map = t.getColumnsAsMap(0, 1);
            } else {
                t = null;
                map = null;
            }
            final boolean extra_processing1;
            extra_processing1 = cla.isOptionSet(EXTRA_PROCESSING_OPTION1);
            String extra_processing1_sep = EXTRA_PROCESSING1_SEP_DEFAULT;
            if (cla.isOptionSet(EXTRA_PROCESSING1_SEP_OPTION)) {
                if (!extra_processing1) {
                    ForesterUtil.fatalError(PRG_NAME,
                            "extra processing is not enabled, cannot set -"
                                    + EXTRA_PROCESSING1_SEP_OPTION + " option");
                }
                if (cla.isOptionValueSet(EXTRA_PROCESSING1_SEP_OPTION)) {
                    extra_processing1_sep = cla.getOptionValue(EXTRA_PROCESSING1_SEP_OPTION);
                } else {
                    ForesterUtil.fatalError(PRG_NAME, "no value for extra processing separator");
                }
            }
            if ((extra_processing1_sep != null) && extra_processing1_sep.equals(separator)) {
                ForesterUtil.fatalError(PRG_NAME,
                        "extra processing separator must not be the same the annotation-separator");
            }
            boolean extra_processing1_keep = EXTRA_PROCESSING1_KEEP_EXTRA_DEFAULT;
            if (cla.isOptionSet(EXTRA_PROCESSING1_KEEP_EXTRA_OPTION)) {
                if (!extra_processing1) {
                    ForesterUtil.fatalError(PRG_NAME,
                            "extra processing is not enabled, cannot set -"
                                    + EXTRA_PROCESSING1_KEEP_EXTRA_OPTION + " option");
                }
                extra_processing1_keep = true;
            }
            Pattern special_pattern = null;
            boolean special_processing = false;
            if (cla.isOptionSet(SPECIAL_PROCESSING_OPTION)) {
                if (extra_processing1) {
                    ForesterUtil
                            .fatalError(PRG_NAME,
                                    "extra processing cannot be used together with special processing pattern");
                }
                if (cla.isOptionValueSet(SPECIAL_PROCESSING_OPTION)) {
                    final String str = cla.getOptionValue(SPECIAL_PROCESSING_OPTION);
                    try {
                        special_pattern = Pattern.compile(str);
                    } catch (final PatternSyntaxException e) {
                        ForesterUtil
                                .fatalError(PRG_NAME,
                                        "error in special processing pattern: " + str + ": " + e.getMessage());
                    }
                    special_processing = true;
                } else {
                    ForesterUtil.fatalError(PRG_NAME, "no value for special processing pattern");
                }
            }
            final boolean remove_annotation_sep;
            remove_annotation_sep = cla.isOptionSet(REMOVE_ANNOT_SEP_OPTION);
            final boolean verbose;
            verbose = cla.isOptionSet(VERBOSE_OPTION);
            final boolean quit;
            quit = cla.isOptionSet(QUIET_OPTION);
            System.out.println("Input tree                 : " + intreefile);
            if (mapping_file != null) {
                System.out.println("Mapping file               : " + mapping_file + " (" + t.getNumberOfRows()
                        + " rows)");
            }
            System.out.println("Annotation-separator       : " + separator);
            if (remove_annotation_sep) {
                System.out.println("Remove anno.-sep. in output: " + remove_annotation_sep);
            }
            System.out.println("Query pattern              : " + pattern);
            if (extra_processing1) {
                System.out.println("Extra processing           : " + extra_processing1);
                System.out.println("Extra processing separator : " + extra_processing1_sep);
                System.out.println("Keep extra annotations     : " + extra_processing1_keep);
            }
            if (special_processing) {
                System.out.println("Special processing         : " + special_processing);
                System.out.println("Special processing pattern : " + special_pattern);
            }
            if (outtablefile != null) {
                System.out.println("Output table               : " + outtablefile);
            }
            if (simple_output) {
                System.out.println("Simple output              : true");
                System.out.println("Cutoff for simple output   : " + simple_output_cutoff);
            }
            Phylogeny[] phys = null;
            try {
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(intreefile, true);
                phys = factory.create(intreefile, pp);
            } catch (final IOException e) {
                ForesterUtil.fatalError(PRG_NAME, "Could not read \"" + intreefile + "\" [" + e.getMessage() + "]");
            }
            if (phys.length == 0) {
                ForesterUtil.fatalError(PRG_NAME, "\"" + intreefile + "\" does not contain any trees");
            }
            System.out.println("Number of input trees      : " + phys.length);
            if (phys.length == 1) {
                System.out.println("Ext. nodes in input tree   : " + phys[0].getNumberOfExternalNodes());
            } else {
                System.out.println("Ext. nodes in input tree 1 : " + phys[0].getNumberOfExternalNodes());
            }
            final EasyWriter outtable_writer;
            if (outtablefile != null) {
                outtable_writer = ForesterUtil.createEasyWriter(outtablefile);
                if (!TWO_COLUMNS_IN_SIMPLE_OUTPUT) {
                    outtable_writer.print("#" + PRG_NAME + " " + PRG_VERSION + " " + PRG_DATE);
                    outtable_writer.print(" Input tree: " + intreefile);
                }
            } else {
                outtable_writer = null;
            }
            int counter = 0;
            for (final Phylogeny phy : phys) {
                if (map != null) {
                    AnalysisMulti.performMapping(pattern, map, phy, verbose);
                }
                if (extra_processing1) {
                    AnalysisMulti.performExtraProcessing1(pattern,
                            phy,
                            extra_processing1_sep,
                            extra_processing1_keep,
                            separator,
                            verbose);
                } else if (special_processing) {
                    AnalysisMulti.performSpecialProcessing1(pattern, phy, separator, special_pattern, verbose);
                }

                final boolean e = AnalysisMulti.likelyProblematicQuery(phy, pattern, 2);
                if (e) {
                    System.out.println("Error: Tree #" + counter);
                    System.out.println(" likely non-homolgous query sequence");
                    continue;
                }
                final ResultMulti res = AnalysisMulti.execute(phy, pattern, separator);
                if (!quit) {
                    if (phys.length == 1) {
                        printResult(res, -1, remove_annotation_sep);
                    } else {
                        printResult(res, counter, remove_annotation_sep);
                    }
                }
                if (outtable_writer != null) {
                    if (simple_output) {
                        writeResultToTableSimple(res,
                                outtable_writer,
                                simple_output_cutoff,
                                TWO_COLUMNS_IN_SIMPLE_OUTPUT,
                                ALLOW_TO_SPLIT_QUERY);
                    } else {
                        writeResultToTable(res, outtable_writer, remove_annotation_sep);
                    }
                    outtable_writer.flush();
                }
                ++counter;
            }
            if (outtable_writer != null) {
                outtable_writer.flush();
                outtable_writer.close();
            }
        } catch (final UserException e) {
            ForesterUtil.fatalError(PRG_NAME, e.getMessage());
        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, e.getMessage());
        } catch (final Exception e) {
            e.printStackTrace();
            ForesterUtil.fatalError(PRG_NAME, "Unexpected error!");
        }
    }

    private static void printResult(final ResultMulti res,
                                    final int counter,
                                    final boolean remove_annotation_sep) {
        System.out.println();
        if (counter == -1) {
            System.out.println("Result for " + res.getQueryNamePrefix());
        } else {
            System.out.println("Result for " + res.getQueryNamePrefix() + " [tree " + counter + "]");
        }
        if ((res.getAllMultiHitPrefixes() == null) | (res.getAllMultiHitPrefixes().size() < 1)) {
            System.out.println(" No match to query pattern!");
        } else {
            System.out.println(" Matching Clade(s):");
            for (final Prefix prefix : res.getCollapsedMultiHitPrefixes()) {
                if (remove_annotation_sep) {
                    System.out.println(" " + prefix.toStringRemoveSeparator());
                } else {
                    System.out.println(" " + prefix);
                }
            }

            if (!ForesterUtil.isEmpty(res.getAllMultiHitPrefixesDown())) {
                System.out.println();
                System.out.println(" Matching Down-tree Bracketing Clade(s):");
                for (final Prefix prefix : res.getCollapsedMultiHitPrefixesDown()) {
                    if (remove_annotation_sep) {
                        System.out.println(" " + prefix.toStringRemoveSeparator());
                    } else {
                        System.out.println(" " + prefix);
                    }
                }
            }
            if (!ForesterUtil.isEmpty(res.getAllMultiHitPrefixesUp())) {
                System.out.println();
                System.out.println(" Matching Up-tree Bracketing Clade(s):");
                for (final Prefix prefix : res.getCollapsedMultiHitPrefixesUp()) {
                    if (remove_annotation_sep) {
                        System.out.println(" " + prefix.toStringRemoveSeparator());
                    } else {
                        System.out.println(" " + prefix);
                    }
                }
            }
            System.out.println();
            System.out.println(" Total Number of Matches: " + res.getNumberOfMatches() + "/"
                    + res.getReferenceTreeNumberOfExternalNodes());
        }
        System.out.println();
    }

    private static void writeResultToTableSimple(final ResultMulti res,
                                                 final EasyWriter w,
                                                 final double cutoff,
                                                 final boolean two_columns,
                                                 final boolean split_query)
            throws IOException {
        if ((res.getAllMultiHitPrefixes() == null) | (res.getAllMultiHitPrefixes().size() < 1)) {
            w.print(res.getQueryNamePrefix());
            w.print("\t");
            w.println("No match to query pattern!");
        } else {
            boolean done = false;
            for (final Prefix prefix : res.getCollapsedMultiHitPrefixes()) {
                if ((prefix.getConfidence() >= cutoff) && !prefix.getPrefix().equals("?")) {
                    if (split_query) {
                        final String[] queries = res.getQueryNamePrefix().split("_");
                        for (final String query : queries) {
                            w.print(query);
                            if (!two_columns) {
                                w.print("\t");
                                w.print("Matching Clades");
                            }
                            w.print("\t");
                            w.print(prefix.getPrefix());
                            if (!two_columns) {
                                w.print("\t");
                                w.print(df.format(prefix.getConfidence()));
                            }
                            w.println();
                        }
                    } else {
                        w.print(res.getQueryNamePrefix());
                        if (!two_columns) {
                            w.print("\t");
                            w.print("Matching Clades");
                        }
                        w.print("\t");
                        w.print(prefix.getPrefix());
                        if (!two_columns) {
                            w.print("\t");
                            w.print(df.format(prefix.getConfidence()));
                        }
                        w.println();
                    }
                    done = true;
                    break;
                }
            }
            if (!done) {
                if (!ForesterUtil.isEmpty(res.getAllMultiHitPrefixesDown())) {
                    for (final Prefix prefix : res.getCollapsedMultiHitPrefixesDown()) {
                        if ((prefix.getConfidence() >= cutoff) && !prefix.getPrefix().equals("?")) {
                            if (split_query) {
                                final String[] queries = res.getQueryNamePrefix().split("_");
                                for (final String query : queries) {
                                    w.print(query);
                                    if (!two_columns) {
                                        w.print("\t");
                                        w.print("Matching Down-tree Bracketing Clades");
                                    }
                                    w.print("\t");
                                    w.print(prefix.getPrefix() + "-like");
                                    if (!two_columns) {
                                        w.print("\t");
                                        w.print(df.format(prefix.getConfidence()));
                                    }
                                    w.println();
                                }
                            } else {
                                w.print(res.getQueryNamePrefix());
                                if (!two_columns) {
                                    w.print("\t");
                                    w.print("Matching Down-tree Bracketing Clades");
                                }
                                w.print("\t");
                                w.print(prefix.getPrefix() + "-like");
                                if (!two_columns) {
                                    w.print("\t");
                                    w.print(df.format(prefix.getConfidence()));
                                }
                                w.println();
                            }
                            done = true;
                            break;
                        }
                    }
                }
            }
            if (!done) {
                if (!ForesterUtil.isEmpty(res.getAllMultiHitPrefixesUp())) {
                    for (final Prefix prefix : res.getCollapsedMultiHitPrefixesUp()) {
                        if ((prefix.getConfidence() >= cutoff) && !prefix.getPrefix().equals("?")) {
                            if (split_query) {
                                final String[] queries = res.getQueryNamePrefix().split("_");
                                for (final String query : queries) {
                                    w.print(query);
                                    if (!two_columns) {
                                        w.print("\t");
                                        w.print("Matching Up-tree Bracketing Clades");
                                    }
                                    w.print("\t");
                                    w.print(prefix.getPrefix() + "-like");
                                    if (!two_columns) {
                                        w.print("\t");
                                        w.print(df.format(prefix.getConfidence()));
                                    }
                                    w.println();
                                }
                            } else {
                                w.print(res.getQueryNamePrefix());
                                if (!two_columns) {
                                    w.print("\t");
                                    w.print("Matching Up-tree Bracketing Clades");
                                }
                                w.print("\t");
                                w.print(prefix.getPrefix() + "-like");
                                if (!two_columns) {
                                    w.print("\t");
                                    w.print(df.format(prefix.getConfidence()));
                                }
                                w.println();
                            }
                            done = true;
                            break;
                        }
                    }
                }
            }
            if (!done) {
                if (split_query) {
                    final String[] queries = res.getQueryNamePrefix().split("_");
                    for (final String query : queries) {
                        w.print(query);
                        w.print("\t");
                        w.println("?");
                    }
                } else {
                    w.print(res.getQueryNamePrefix());
                    w.print("\t");
                    w.println("?");
                }
            }
        }
    }

    private static void writeResultToTable(final ResultMulti res,
                                           final EasyWriter w,
                                           final boolean remove_annotation_sep)
            throws IOException {
        if ((res.getAllMultiHitPrefixes() == null) | (res.getAllMultiHitPrefixes().size() < 1)) {
            w.print(res.getQueryNamePrefix());
            w.print("\t");
            w.println("No match to query pattern!");
        } else {
            for (final Prefix prefix : res.getCollapsedMultiHitPrefixes()) {
                w.print(res.getQueryNamePrefix());
                w.print("\t");
                w.print("Matching Clades");
                w.print("\t");
                if (remove_annotation_sep) {
                    w.print(prefix.getPrefixRemoveSeparator());
                } else {
                    w.print(prefix.getPrefix());
                }
                w.print("\t");
                w.print(df.format(prefix.getConfidence()));
                w.print("\t");
                w.print(String.valueOf(res.getNumberOfMatches()));
                w.print("\t");
                w.print(String.valueOf(res.getReferenceTreeNumberOfExternalNodes()));
                w.println();
            }

            if (!ForesterUtil.isEmpty(res.getAllMultiHitPrefixesDown())) {
                for (final Prefix prefix : res.getCollapsedMultiHitPrefixesDown()) {
                    w.print(res.getQueryNamePrefix());
                    w.print("\t");
                    w.print("Matching Down-tree Bracketing Clades");
                    w.print("\t");
                    if (remove_annotation_sep) {
                        w.print(prefix.getPrefixRemoveSeparator());
                    } else {
                        w.print(prefix.getPrefix());
                    }
                    w.print("\t");
                    w.print(df.format(prefix.getConfidence()));
                    w.print("\t");
                    w.print(String.valueOf(res.getNumberOfMatches()));
                    w.print("\t");
                    w.print(String.valueOf(res.getReferenceTreeNumberOfExternalNodes()));
                    w.println();
                }
            }
            if (!ForesterUtil.isEmpty(res.getAllMultiHitPrefixesUp())) {
                for (final Prefix prefix : res.getCollapsedMultiHitPrefixesUp()) {
                    w.print(res.getQueryNamePrefix());
                    w.print("\t");
                    w.print("Matching Up-tree Bracketing Clades");
                    w.print("\t");
                    if (remove_annotation_sep) {
                        w.print(prefix.getPrefixRemoveSeparator());
                    } else {
                        w.print(prefix.getPrefix());
                    }
                    w.print("\t");
                    w.print(df.format(prefix.getConfidence()));
                    w.print("\t");
                    w.print(String.valueOf(res.getNumberOfMatches()));
                    w.print("\t");
                    w.print(String.valueOf(res.getReferenceTreeNumberOfExternalNodes()));
                    w.println();
                }
            }
        }
    }

    private static void print_help() {
        System.out.println("Usage:");
        System.out.println();
        System.out.println(PRG_NAME + " [options] <input tree(s) file> [output table file]");
        System.out.println();
        System.out.println(" options:");
        System.out.println("  -" + SEP_OPTION + "=<separator>     : the annotation-separator to be used (default: \""
                + SEP_DEFAULT + "\")");
        System.out.println("  -" + MAPPING_FILE_OPTION
                + "=<mapping table> : to map node names to appropriate annotations (tab-separated, two columns) (default: no mapping)");
        System.out.println("  -" + EXTRA_PROCESSING_OPTION1
                + "                 : to enable extra processing of annotations (e.g. \"Q16611|A.1.1\" becomes \"A.1.1\")");
        System.out.println("  -" + EXTRA_PROCESSING1_SEP_OPTION
                + "=<separator>    : the separator for extra annotations (default: \"" + EXTRA_PROCESSING1_SEP_DEFAULT
                + "\")");
        System.out.println("  -" + EXTRA_PROCESSING1_KEEP_EXTRA_OPTION
                + "                : to keep extra annotations (e.g. \"Q16611|A.1.1\" becomes \"A.1.1.Q16611\")");
        System.out.println("  -" + SPECIAL_PROCESSING_OPTION
                + "=<pattern>       : special processing with pattern (e.g. \"(\\d+)([a-z]+)_.+\" for changing \"6q_EF42\" to \"6.q\")");
        System.out.println("  -" + REMOVE_ANNOT_SEP_OPTION
                + "                : to remove the annotation-separator in the output (e.g. the \"" + SEP_DEFAULT
                + "\")");
        System.out.println("  -" + VERBOSE_OPTION + "                 : verbose");
        System.out.println("  -" + QUIET_OPTION
                + "                 : quiet (no output to console, for when used in a pipeline)");
        System.out.println("  --" + QUERY_PATTERN_OPTION
                + "=<pattern>      : expert option: the regular expression pattern for the query (default: \""
                + QUERY_PATTERN_DEFAULT + "\" for pplacer output)");
        System.out.println("  -" + SIMPLE_OUTPUT_OPTION + "               : simple output");
        System.out.println("  -" + SIMPLE_OUTPUT_CUTOFF_OPTION
                + "              : probability cutoff for simple output (0.0-1.0, default 0.0)");
        System.out.println();
        System.out.println("Examples:");
        System.out.println();
        System.out.println(" " + PRG_NAME + " pp_out_tree.sing.tre result.tsv");
        System.out.println(" " + PRG_NAME + " -s=. pp_out_tree.sing.tre result.tsv");
        System.out.println(" " + PRG_NAME + " -s=_ -m=map.tsv pp_out_trees.sing.tre result.tsv");
        System.out.println(" " + PRG_NAME + " -x -xs=& -xk pp_out_trees.sing.tre result.tsv");
        System.out.println(" " + PRG_NAME + " -x -xs=\"|\" pp_out_trees.sing.tre result.tsv");
        System.out.println(" " + PRG_NAME + " -x -xk -m=map.tsv pp_out_trees.sing.tre result.tsv");
        System.out.println(" " + PRG_NAME + " -m=map.tsv -S='(\\d+)([a-z?]*)_.+' pp_out_trees.sing.tre result.tsv");
        System.out.println();
    }
}
