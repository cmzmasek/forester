// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2017 Christian M. Zmasek
// Copyright (C) 2017 J. Craig Venter Institute
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
// --------------------
// TODO
// * Multiple "hits" with different "M" values
// * More tests (including multiple children per node), especially on edge cases
// * Utilize relevant support values for warnings

package org.forester.clade_analysis;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.ForesterUtil;
import org.forester.util.UserException;

public final class AnalysisMulti {

    private final static String UNKNOWN = "?";
    public final static double DEFAULT_CUTOFF_FOR_SPECIFICS = 0.5;
    public final static String DEFAULT_SEPARATOR = ".";
    public final static Pattern DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE = Pattern.compile("_#\\d+_M=(.+)");

    public static ResultMulti execute(final Phylogeny p) throws UserException {
        return execute(p, DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE, DEFAULT_SEPARATOR, DEFAULT_CUTOFF_FOR_SPECIFICS);
    }

    public static ResultMulti execute(final Phylogeny p, final String separator) throws UserException {
        return execute(p, DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE, separator, DEFAULT_CUTOFF_FOR_SPECIFICS);
    }

    public static ResultMulti execute(final Phylogeny p, final String separator, final double cutoff_for_specifics)
            throws UserException {
        return execute(p, DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE, separator, cutoff_for_specifics);
    }

    public static ResultMulti execute(final Phylogeny p, final double cutoff_for_specifics) throws UserException {
        return execute(p, DEFAULT_QUERY_PATTERN_FOR_PPLACER_TYPE, DEFAULT_SEPARATOR, cutoff_for_specifics);
    }

    public static String likelyProblematicQuery(final Phylogeny p,
                               final Pattern query, final double factor
    ) {
        final List<PhylogenyNode> qnodes = p.getNodes(query);
        BasicDescriptiveStatistics s = new BasicDescriptiveStatistics();
        for (final PhylogenyNodeIterator it = p.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode i = it.next();
            if (!qnodes.contains(i)) {
                s.addValue(i.calculateDistanceToRoot());
            }
        }
        //final double mean_distance_to_root = s.arithmeticMean();
        final double max_distance_to_root = s.getMax();
        //final double sd = s.sampleStandardDeviation();
        boolean all_outliers = true;
        for (final PhylogenyNode q : qnodes) {
            if (q.calculateDistanceToRoot() < factor * max_distance_to_root) {
                all_outliers = false;
                break;
            }
        }
        if ( all_outliers ) {
            return "Error: Possible non homologous query sequence";
        }
        else {
            return null;
        }
    }


    public static ResultMulti execute(final Phylogeny p,
                                      final Pattern query,
                                      final String separator,
                                      final double cutoff_for_specifics)
            throws UserException {
        if (ForesterUtil.isEmpty(separator)) {
            throw new UserException("separator must not be null or empty");
        }
        cleanUpExternalNames(p, separator);
        final List<PhylogenyNode> qnodes = p.getNodes(query);
        final ResultMulti res = new ResultMulti();
        res.setQueryNamePrefix(obtainQueryPrefix(query, qnodes));
        res.setTotalNumberOfMatches(qnodes.size());
        res.setReferenceTreeNumberOfExternalNodes(p.getNumberOfExternalNodes() - qnodes.size());
        for (int i = 0; i < qnodes.size(); ++i) {
            final PhylogenyNode qnode = qnodes.get(i);
            if (qnode.isRoot()) {
                throw new UserException("query " + query + " is root");
            }
            if (qnode.getParent().isRoot()) {
                throw new UserException("parent of query " + query + " is root");
            }
            PhylogenyNode qnode_p = qnode.getParent();
            PhylogenyNode qnode_pp = qnode.getParent().getParent();
            //This is to deal with internal nodes with 1 descendant.
            while (qnode_p.getNumberOfDescendants() == 1) {
                qnode_p = qnode_p.getParent();
            }
            while (qnode_pp.getNumberOfDescendants() == 1) {
                qnode_pp = qnode_pp.getParent();
            }
            final List<String> qnode_ext_nodes_names = new ArrayList<String>();
            for (final PhylogenyNode qnode_ext_node : qnode_pp.getAllExternalDescendants()) {
                final String name = qnode_ext_node.getName();
                final Matcher m = query.matcher(name);
                if (!m.find()) {
                    qnode_ext_nodes_names.add(name);
                }
            }
            final String greatest_common_prefix = ForesterUtil.greatestCommonPrefix(qnode_ext_nodes_names, separator);
            final Matcher matcher = query.matcher(qnode.getName());
            String conf_str = null;
            if (matcher.find()) {
                conf_str = matcher.group(1);
            } else {
                throw new IllegalStateException("ERROR: query pattern does not match [this should have never happened!]");
            }
            final double conf;
            try {
                conf = Double.parseDouble(conf_str);
            }
            catch (final NumberFormatException ex) {
                throw new UserException("ERROR: Could not parse confidence from \"" + conf_str +"\" from query node name \"" + qnode.getName() +"\"" );
            }
            if (!ForesterUtil.isEmpty(greatest_common_prefix)) {
                res.addGreatestCommonPrefix(greatest_common_prefix, conf);
            } else {
                res.addGreatestCommonPrefix(UNKNOWN, conf);
            }
            final String greatest_common_prefix_up = analyzeSiblings(qnode_p, qnode_pp, separator, query);
            if (!ForesterUtil.isEmpty(greatest_common_prefix_up)) {
                res.addGreatestCommonPrefixUp(greatest_common_prefix_up, conf);
            } else {
                res.addGreatestCommonPrefixUp(UNKNOWN, conf);
            }
            final String greatest_common_prefix_down = analyzeSiblings(qnode, qnode_p, separator, query);
            if (!ForesterUtil.isEmpty(greatest_common_prefix_down)) {
                res.addGreatestCommonPrefixDown(greatest_common_prefix_down, conf);
            } else {
                res.addGreatestCommonPrefixDown(UNKNOWN, conf);
            }
        }
        res.analyze(cutoff_for_specifics);
        return res;
    }

    private final static String obtainQueryPrefix(final Pattern query, final List<PhylogenyNode> qnodes)
            throws UserException {
        String query_name_prefix = null;
        for (final PhylogenyNode n : qnodes) {
            final String name = n.getName();
            final Matcher matcher = query.matcher(name);
            if (matcher.find()) {
                final String prefix = name.substring(0, matcher.start());
                if (ForesterUtil.isEmpty(prefix)) {
                    throw new UserException("query nodes with empty label prefix found: \"" + prefix + "\"");
                }
                if (query_name_prefix == null) {
                    query_name_prefix = prefix;
                } else if (!query_name_prefix.equals(prefix)) {
                    throw new UserException("query nodes with different label prefixes found: \"" + query_name_prefix
                            + "\" and \"" + prefix + "\"");
                }
            }
        }
        return query_name_prefix;
    }

    private final static void cleanUpExternalNames(final Phylogeny p, final String separator) throws UserException {
        final Pattern pattern1 = Pattern.compile("\\Q" + separator + "\\E" + "\\s+");
        final Pattern pattern2 = Pattern.compile("\\s+" + "\\Q" + separator + "\\E");
        final Pattern pattern3 = Pattern.compile("\\Q" + separator + separator + "\\E");
        final PhylogenyNodeIterator it = p.iteratorExternalForward();
        while (it.hasNext()) {
            final PhylogenyNode node = it.next();
            final String name = node.getName().trim();
            if (ForesterUtil.isEmpty(name)) {
                throw new UserException("external node(s) with empty annotation found");
            }
            if (name.endsWith(separator)) {
                throw new UserException("illegally formatted annotation found: annotations cannot end with separator: "
                        + name);
            }
            if (name.startsWith(separator)) {
                throw new UserException("illegally formatted annotation found: annotations cannot start with separator: "
                        + name);
            }
            if (pattern1.matcher(name).find()) {
                throw new UserException("illegally formatted annotation found: separator followed by whitespace: "
                        + name);
            }
            if (pattern2.matcher(name).find()) {
                throw new UserException("illegally formatted annotation found: whitespace followed by separator: "
                        + name);
            }
            if (pattern3.matcher(name).find()) {
                throw new UserException("illegally formatted annotation found: empty annotation level: " + name);
            }
            node.setName(name.replaceAll("\\s+", " "));
        }
    }

    private final static String analyzeSiblings(final PhylogenyNode child,
                                                final PhylogenyNode parent,
                                                final String separator,
                                                final Pattern query) {
        final int child_index = child.getChildNodeIndex();
        final List<String> ext_nodes_names = new ArrayList<String>();
        final List<PhylogenyNode> descs = parent.getDescendants();
        for (int i = 0; i < descs.size(); ++i) {
            if (i != child_index) {
                final PhylogenyNode d = descs.get(i);
                for (final PhylogenyNode n : d.getAllExternalDescendants()) {
                    final String name = n.getName();
                    final Matcher m = query.matcher(name);
                    if (!m.find()) {
                        ext_nodes_names.add(name);
                    }
                }
            }
        }
        final String greatest_common_prefix = ForesterUtil.greatestCommonPrefix(ext_nodes_names, separator);
        return greatest_common_prefix;
    }

    public final static void performMapping(final Pattern pattern,
                                            final SortedMap<String, String> map,
                                            final Phylogeny p,
                                            final boolean verbose)
            throws UserException {
        if (verbose) {
            System.out.println();
            System.out.println("Id to annotation mapping:");
        }
        final PhylogenyNodeIterator it = p.iteratorExternalForward();
        while (it.hasNext()) {
            final PhylogenyNode node = it.next();
            final String name = node.getName().trim();
            if (ForesterUtil.isEmpty(name)) {
                throw new UserException("external node with empty name found");
            }
            final Matcher m = pattern.matcher(name);
            if (!m.find()) {
                if (!map.containsKey(name)) {
                    throw new UserException("no mapping for \"" + name + "\" found");
                }
                node.setName(map.get(name).trim());
                if (verbose) {
                    System.out.println(name + " -> " + node.getName());
                }
            }
        }
        if (verbose) {
            System.out.println();
        }
    }

    public final static void performExtraProcessing1(final Pattern query_pattern,
                                                     final Phylogeny p,
                                                     final String extra_sep,
                                                     final boolean keep,
                                                     final String annotation_sep,
                                                     final boolean verbose)
            throws UserException {
        if (verbose) {
            System.out.println();
            System.out.println("Extra annotation processing:");
        }
        final PhylogenyNodeIterator it = p.iteratorExternalForward();
        while (it.hasNext()) {
            final PhylogenyNode node = it.next();
            final String name = node.getName().trim();
            if (ForesterUtil.isEmpty(name)) {
                throw new UserException("external node with empty name found");
            }
            if (!query_pattern.matcher(name).find()) {
                final StringBuilder sb = new StringBuilder();
                final int last_index = name.lastIndexOf(extra_sep);
                if (last_index >= 0) {
                    final String annotation = name.substring(last_index + 1).trim();
                    if (ForesterUtil.isEmptyTrimmed(annotation)) {
                        throw new UserException("llegally formatted annotation: " + name);
                    }
                    if (keep) {
                        final String extra = name.substring(0, last_index).trim();
                        sb.append(annotation);
                        if (!ForesterUtil.isEmpty(extra)) {
                            sb.append(annotation_sep);
                            sb.append(extra);
                        }
                    } else {
                        sb.append(annotation);
                    }
                    node.setName(sb.toString());
                    if (verbose) {
                        System.out.println(name + " -> " + node.getName());
                    }
                }
            }
        }
        if (verbose) {
            System.out.println();
        }
    }

    public final static void performSpecialProcessing1(final Pattern query_pattern,
                                                       final Phylogeny p,
                                                       final String annotation_sep,
                                                       final Pattern special_pattern,
                                                       final boolean verbose)
            throws UserException {
        if (verbose) {
            System.out.println();
            System.out.println("Special annotation processing:");
        }
        final PhylogenyNodeIterator it = p.iteratorExternalForward();
        while (it.hasNext()) {
            final PhylogenyNode node = it.next();
            final String name = node.getName().trim();
            if (ForesterUtil.isEmpty(name)) {
                throw new UserException("external node with empty name found");
            }
            if (!query_pattern.matcher(name).find()) {
                final Matcher special_m = special_pattern.matcher(name);
                if (special_m.matches()) {
                    final int c = special_m.groupCount();
                    if (c < 1) {
                        throw new UserException("illegal special pattern: " + special_pattern
                                + " (need at least one capturing group)");
                    }
                    final StringBuilder sb = new StringBuilder();
                    for (int i = 1; i <= c; ++i) {
                        final String g = special_m.group(i);
                        if (!ForesterUtil.isEmpty(g)) {
                            if (i > 1) {
                                sb.append(annotation_sep);
                            }
                            sb.append(special_m.group(i));
                        }
                    }
                    node.setName(sb.toString());
                    if (verbose) {
                        System.out.println(name + " -> " + node.getName());
                    }
                } else {
                    throw new UserException("illegally formatted annotation for special processing: " + name
                            + " (expected pattern: " + special_pattern + ")");
                }
            }
        }
        if (verbose) {
            System.out.println();
        }
    }
}
