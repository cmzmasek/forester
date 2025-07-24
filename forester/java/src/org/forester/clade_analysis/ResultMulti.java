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

package org.forester.clade_analysis;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.util.ForesterUtil;
import org.forester.util.UserException;

public final class ResultMulti {

    private final static double MIN_DIFF = 1E-5;
    private final String _separator;
    private final List<Prefix> _greatest_common_prefixes = new ArrayList<Prefix>();
    private final List<Prefix> _greatest_common_prefixes_up = new ArrayList<Prefix>();
    private final List<Prefix> _greatest_common_prefixes_down = new ArrayList<Prefix>();
    private List<Prefix> _all = null;
    private List<Prefix> _collapsed = null;


    private List<Prefix> _all_up = null;
    private List<Prefix> _collapsed_up = null;

    private List<Prefix> _all_down = null;
    private List<Prefix> _collapsed_down = null;
    private int _matches = 0;
    private int _ref_tree_ext_nodes = 0;
    private String _query_name_prefix = "";
    private final static DecimalFormat df = new DecimalFormat("0.0###");

    ResultMulti(final String separator) {
        if (ForesterUtil.isEmpty(separator)) {
            throw new IllegalArgumentException("separator must not be null or empty");
        }
        _separator = separator;
        reset();
    }

    ResultMulti() {
        _separator = AnalysisMulti.DEFAULT_SEPARATOR;
        reset();
    }

    public List<Prefix> getAllMultiHitPrefixesUp() {
        return _all_up;
    }

    public List<Prefix> getCollapsedMultiHitPrefixesUp() {
        return _collapsed_up;
    }


    public List<Prefix> getAllMultiHitPrefixesDown() {
        return _all_down;
    }

    public List<Prefix> getCollapsedMultiHitPrefixesDown() {
        return _collapsed_down;
    }


    public List<Prefix> getAllMultiHitPrefixes() {
        return _all;
    }

    public List<Prefix> getCollapsedMultiHitPrefixes() {
        return _collapsed;
    }


    public String getQueryNamePrefix() {
        return _query_name_prefix;
    }

    public int getNumberOfMatches() {
        return _matches;
    }

    public int getReferenceTreeNumberOfExternalNodes() {
        return _ref_tree_ext_nodes;
    }


    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append("Query: ");
        sb.append(getQueryNamePrefix());
        sb.append(ForesterUtil.LINE_SEPARATOR);
        sb.append("Matching Clade(s):");
        sb.append(ForesterUtil.LINE_SEPARATOR);
        for (final Prefix prefix : _collapsed) {
            sb.append(prefix);
            sb.append(ForesterUtil.LINE_SEPARATOR);
        }

        if (!ForesterUtil.isEmpty(_all_down)) {
            sb.append(ForesterUtil.LINE_SEPARATOR);
            sb.append("Matching Down-tree Bracketing Clade(s):");
            sb.append(ForesterUtil.LINE_SEPARATOR);
            for (final Prefix prefix : _collapsed_down) {
                sb.append(prefix);
                sb.append(ForesterUtil.LINE_SEPARATOR);
            }
        }
        if (!ForesterUtil.isEmpty(_all_up)) {
            sb.append(ForesterUtil.LINE_SEPARATOR);
            sb.append("Matching Up-tree Bracketing Clade(s):");
            sb.append(ForesterUtil.LINE_SEPARATOR);
            for (final Prefix prefix : _collapsed_up) {
                sb.append(prefix);
                sb.append(ForesterUtil.LINE_SEPARATOR);
            }
        }
        sb.append(ForesterUtil.LINE_SEPARATOR);
        sb.append("Total Number of Matches: " + getNumberOfMatches() + "/" + getReferenceTreeNumberOfExternalNodes());
        sb.append(ForesterUtil.LINE_SEPARATOR);
        return sb.toString();
    }

    void addGreatestCommonPrefix(final String prefix, final double confidence) {
        _greatest_common_prefixes.add(new Prefix(prefix, confidence, _separator));
    }

    void addGreatestCommonPrefixUp(final String prefix_up, final double confidence) {
        _greatest_common_prefixes_up.add(new Prefix(prefix_up, confidence, _separator));
    }

    void addGreatestCommonPrefixDown(final String prefix_down, final double confidence) {
        _greatest_common_prefixes_down.add(new Prefix(prefix_down, confidence, _separator));
    }

    void setQueryNamePrefix(final String query_name_prefix) {
        if (!ForesterUtil.isEmpty(_query_name_prefix)) {
            throw new IllegalStateException("illegal attempt to change the query name prefix");
        }
        _query_name_prefix = query_name_prefix;
    }

    void setTotalNumberOfMatches(final int matches) {
        if (_matches > 0) {
            throw new IllegalStateException("illegal attempt to change the number of matches");
        }
        _matches = matches;
    }

    public void setReferenceTreeNumberOfExternalNodes(final int ext_nodes) {
        if (_ref_tree_ext_nodes > 0) {
            throw new IllegalStateException("illegal attempt to change the number of external nodes");
        }
        _ref_tree_ext_nodes = ext_nodes;
    }

    void analyze() throws UserException {
        reset();
        analyzeGreatestCommonPrefixes(_greatest_common_prefixes, _separator);
        analyzeGreatestCommonPrefixesUp(_greatest_common_prefixes_up, _separator);
        analyzeGreatestCommonPrefixesDown(_greatest_common_prefixes_down, _separator);
    }

    private void reset() {
        _all = new ArrayList<Prefix>();
        _collapsed = new ArrayList<Prefix>();
        _all_up = new ArrayList<Prefix>();
        _collapsed_up = new ArrayList<Prefix>();
        _all_down = new ArrayList<Prefix>();
        _collapsed_down = new ArrayList<Prefix>();
    }

    private void analyzeGreatestCommonPrefixes(final List<Prefix> greatest_common_prefixes,
                                               final String separator)
            throws UserException {
        final List<Prefix> l = obtainAllPrefixes(greatest_common_prefixes, separator);
        if (!ForesterUtil.isEmpty(l)) {
            sortPrefixesAccordingToConfidence(l);
            _all = removeLessSpecificPrefixes(l, separator);
            _collapsed = collapse(_all);

        }
    }

    private void analyzeGreatestCommonPrefixesUp(final List<Prefix> greatest_common_prefixes_up,
                                                 final String separator)
            throws UserException {
        final List<Prefix> l = obtainAllPrefixes(greatest_common_prefixes_up, separator);
        if (!ForesterUtil.isEmpty(l)) {
            sortPrefixesAccordingToConfidence(l);
            _all_up = removeLessSpecificPrefixes(l, separator);
            _collapsed_up = collapse(_all_up);

        }
    }

    void analyzeGreatestCommonPrefixesDown(final List<Prefix> greatest_common_prefixes_down,
                                           final String separator)
            throws UserException {
        final List<Prefix> l = obtainAllPrefixes(greatest_common_prefixes_down, separator);
        if (!ForesterUtil.isEmpty(l)) {
            sortPrefixesAccordingToConfidence(l);
            _all_down = removeLessSpecificPrefixes(l, separator);
            _collapsed_down = collapse(_all_down);

        }
    }

    static List<Prefix> obtainSpecifics(final double cutoff,
                                        final List<Prefix> cleaned,
                                        final List<Prefix> collapsed,
                                        final String separator) {
        final List<Prefix> cleaned_spec = new ArrayList<Prefix>();
        final Set<String> collapsed_set = new HashSet<String>();
        for (final Prefix prefix : collapsed) {
            collapsed_set.add(prefix.getPrefix());
        }
        final List<Prefix> spec = new ArrayList<Prefix>();
        for (final Prefix prefix : cleaned) {
            if ((prefix.getConfidence() >= cutoff) && !collapsed_set.contains(prefix.getPrefix())) {
                spec.add(prefix);
            }
        }
        for (final Prefix o : spec) {
            boolean ok = true;
            for (final Prefix i : spec) {
                if ((!o.getPrefix().equals(i.getPrefix()))
                        && (ForesterUtil.isContainsPrefix(i.getPrefix(), o.getPrefix(), separator))) {
                    ok = false;
                    break;
                }

            }
            if (ok) {
                cleaned_spec.add(o);
            }
        }
        return cleaned_spec;
    }

    private static List<Prefix> collapse(final List<Prefix> cleaned) throws UserException {
        final List<Prefix> collapsed = new ArrayList<Prefix>();
        final Set<String> firsts = new HashSet<String>();
        double confidence_sum = 0;
        for (final Prefix prefix : cleaned) {
            final String f = prefix.getPrefixFirstElement();
            if (!firsts.contains(f)) {
                firsts.add(f);
                collapsed.add(prefix);
                confidence_sum += prefix.getConfidence();
            }
        }
        if (!ForesterUtil.isEqual(confidence_sum, 1.0, MIN_DIFF)) {
            throw new UserException("ERROR: confidences add up to " + confidence_sum + " instead of 1.0");
        }
        return collapsed;
    }

    /*
     * This replaces (by way of example)
     * A.1.1 0.9
     * A.1   0.9
     * with
     * A.1.1 0.9
     *
     * I.e. it removes less specific prefixes.
     *
     */
    private static List<Prefix> removeLessSpecificPrefixes(final List<Prefix> l, final String separator) {
        final List<Prefix> cleaned = new ArrayList<Prefix>();
        for (final Prefix o : l) {
            boolean ok = true;
            for (final Prefix i : l) {

                if ((!o.getPrefix().equals(i.getPrefix()))
                        && (ForesterUtil.isContainsPrefix(i.getPrefix(), o.getPrefix(), separator))
                        && ForesterUtil.isEqual(i.getConfidence(),
                        o.getConfidence())) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                cleaned.add(o);
            }
        }
        return cleaned;
    }

    private static void sortPrefixesAccordingToConfidence(final List<Prefix> l) {
        Collections.sort(l, new Comparator<Prefix>() {

            @Override
            public int compare(final Prefix x, final Prefix y) {
                return compare(x.getConfidence(), y.getConfidence());
            }

            private int compare(final double a, final double b) {
                return a > b ? -1 : a > b ? 1 : 0;
            }
        });
    }

    private static List<Prefix> obtainAllPrefixes(final List<Prefix> greatest_common_prefixes,
                                                  final String separator) {
        final SortedMap<String, Double> map = new TreeMap<String, Double>();
        for (final Prefix prefix : greatest_common_prefixes) {
            final List<String> prefixes = ForesterUtil.spliIntoPrefixes(prefix.getPrefix(), separator);
            for (final String p : prefixes) {
                map.put(p, 0.0);
            }
        }
        for (final String key : map.keySet()) {
            for (final Prefix prefix : greatest_common_prefixes) {
                if (ForesterUtil.isContainsPrefix(prefix.getPrefix(), key, separator)) {
                    map.put(key, map.get(key) + prefix.getConfidence());
                }
            }
        }
        final List<Prefix> l = new ArrayList<Prefix>();
        for (final Entry<String, Double> entry : map.entrySet()) {
            l.add(new Prefix(entry.getKey(), entry.getValue(), separator));
        }
        return l;
    }
}
