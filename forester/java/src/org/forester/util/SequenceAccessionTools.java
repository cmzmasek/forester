// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.util;

import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Accession.Source;
import org.forester.phylogeny.data.Sequence;

public final class SequenceAccessionTools {

    // The format for GenBank Accession numbers are:
    // Nucleotide: 1 letter + 5 numerals OR 2 letters + 6 numerals
    // Protein: 3 letters + 5 numerals
    // http://www.ncbi.nlm.nih.gov/Sequin/acc.html
    public final static Pattern GENBANK_NUC_PATTERN_1 = Pattern
            .compile("(?:\\A|.*[^a-zA-Z0-9])([A-Z]\\d{5}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)");
    public final static Pattern GENBANK_NUC_PATTERN_2 = Pattern
            .compile("(?:\\A|.*[^a-zA-Z0-9])([A-Z]{2}\\d{6}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)");
    public final static Pattern GENBANK_PROT_PATTERN = Pattern
            .compile("(?:\\A|.*[^a-zA-Z0-9])([A-Z]{3}\\d{5}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)");

    //public final static Pattern GENBANK_PROT_PATTERN_YP = Pattern
    //        .compile("(?:\\A|.*[^a-zA-Z0-9])(YP_\\d{9}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)");
    public final static Pattern GENBANK_PROT_PATTERN_2 = Pattern
            .compile("(?:\\A|.*[^a-zA-Z0-9])([A-Z]{3}\\d{5}(?:\\.\\d+))(?:[^a-zA-Z0-9]|\\Z)");
    public final static Pattern GI_PATTERN = Pattern
            .compile("(?:\\b|_)(?:GI|gi)[|_=:](\\d+)(?:\\b|_)");
    public final static String UNIPROT_KB_BASE_PATTERN_STR = "((?:[OPQ][0-9][A-Z0-9]{3}[0-9])|(?:[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}))";
    public final static Pattern UNIPROT_KB_PATTERN_0 = Pattern
            .compile("(?:\\b|_)" + UNIPROT_KB_BASE_PATTERN_STR + "(?:\\b|_)");
    public final static Pattern UNIPROT_KB_PATTERN_1 = Pattern
            .compile("(?:\\b|_)(?:sp|tr)[\\.|\\-_=/\\\\]" + UNIPROT_KB_BASE_PATTERN_STR + "(?:\\b|_)");
    public final static Pattern UNIPROT_KB_PATTERN_2 = Pattern.compile("(?:\\b|_)(?:[A-Z0-9]{2,5}|"
            + UNIPROT_KB_BASE_PATTERN_STR + ")_(([A-Z9][A-Z]{2}[A-Z0-9]{2})|RAT|PIG|PEA)(?:\\b|_)");
    public final static Pattern ENSEMBL_PATTERN = Pattern.compile("(?:\\b|_)(ENS[A-Z]*[0-9]+)(?:\\b|_)");
    // RefSeq accession numbers can be distinguished from GenBank accessions
    // by their distinct prefix format of 2 characters followed by an
    // underscore character ('_'). For example, a RefSeq protein accession is NP_015325.
    private final static Pattern REFSEQ_PATTERN = Pattern
            .compile("(?:\\A|.*[^a-zA-Z0-9])([A-Z]{2}_\\d{6,}(\\.\\d)?)(?:[^a-zA-Z0-9]|\\Z)");
    //.compile("(?:\\A|.*[^a-zA-Z0-9])(YP_\\d{9}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)");


    private SequenceAccessionTools() {
        // Hiding the constructor.
    }

    public final static boolean isProteinDbQuery(final String query) {
        final String r1 = parseRefSeqAccessorFromString(query);
        if (!ForesterUtil.isEmpty(r1) && (r1.charAt(1) == 'P')) {
            return true;
        }
        final String r2 = parseUniProtAccessorFromString(query);
        if (!ForesterUtil.isEmpty(r2)) {
            return true;
        }
        return GENBANK_PROT_PATTERN.matcher(query).lookingAt();
    }

    public final static Accession obtainAccessorFromSequence(final Sequence seq) {
        String a = obtainUniProtAccessorFromSequence(seq);
        if (!ForesterUtil.isEmpty(a)) {
            return new Accession(a, Source.UNIPROT);
        }
        a = obtainGenbankAccessorFromSequence(seq);
        if (!ForesterUtil.isEmpty(a)) {
            return new Accession(a, Source.NCBI);
        }
        a = obtainRefSeqAccessorFromSequence(seq);
        if (!ForesterUtil.isEmpty(a)) {
            return new Accession(a, Source.REFSEQ);
        }
        a = obtainGiNumberFromSequence(seq);
        if (!ForesterUtil.isEmpty(a)) {
            return new Accession(a, Source.GI);
        }
        return null;
    }

    public final static String obtainUniProtAccessorFromSequence(final Sequence seq) {
        if (seq == null) {
            return null;
        }
        String a = null;
        if (!ForesterUtil.isEmpty(seq.getSymbol())) {
            a = parseUniProtAccessorFromString(seq.getSymbol());
        }
        if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getName())) {
            a = parseUniProtAccessorFromString(seq.getName());
        }
        if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getGeneName())) {
            a = parseUniProtAccessorFromString(seq.getGeneName());
        }
        if (ForesterUtil.isEmpty(a) && (seq.getAccession() != null)
                && !ForesterUtil.isEmpty(seq.getAccession().getValue())) {
            a = parseUniProtAccessorFromString(seq.getAccession().getValue());
        }
        return a;
    }

    public final static String obtainGenbankAccessorFromSequence(final Sequence seq) {
        if (seq == null) {
            return null;
        }
        String a = null;
        if (!ForesterUtil.isEmpty(seq.getSymbol())) {
            a = parseGenbankAccessorFromString(seq.getSymbol());
        }
        if (!ForesterUtil.isEmpty(seq.getGeneName())) {
            a = parseGenbankAccessorFromString(seq.getGeneName());
        }
        if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getName())) {
            a = parseGenbankAccessorFromString(seq.getName());
        }
        if (ForesterUtil.isEmpty(a) && (seq.getAccession() != null)
                && !ForesterUtil.isEmpty(seq.getAccession().getValue())) {
            a = parseGenbankAccessorFromString(seq.getAccession().getValue());
        }
        return a;
    }

    public final static String obtainRefSeqAccessorFromSequence(final Sequence seq) {
        if (seq == null) {
            return null;
        }
        String a = null;
        if (!ForesterUtil.isEmpty(seq.getSymbol())) {
            a = parseRefSeqAccessorFromString(seq.getSymbol());
        }
        if (!ForesterUtil.isEmpty(seq.getGeneName())) {
            a = parseRefSeqAccessorFromString(seq.getGeneName());
        }
        if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getName())) {
            a = parseRefSeqAccessorFromString(seq.getName());
        }
        if (ForesterUtil.isEmpty(a) && (seq.getAccession() != null)
                && !ForesterUtil.isEmpty(seq.getAccession().getValue())) {
            a = parseRefSeqAccessorFromString(seq.getAccession().getValue());
        }
        return a;
    }

    public final static String obtainGiNumberFromSequence(final Sequence seq) {
        if (seq == null) {
            return null;
        }
        String a = null;
        if (!ForesterUtil.isEmpty(seq.getName())) {
            a = parseGInumberFromString(seq.getName());
        }
        if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getGeneName())) {
            a = parseGInumberFromString(seq.getGeneName());
        }
        if (ForesterUtil.isEmpty(a) && (seq.getAccession() != null)
                && !ForesterUtil.isEmpty(seq.getAccession().getValue())) {
            a = parseGInumberFromString(seq.getAccession().getValue());
        }
        return a;
    }

    public final static Accession obtainAccessorFromDataFields(final PhylogenyNode n) {
        String a = obtainUniProtAccessorFromDataFields(n);
        if (!ForesterUtil.isEmpty(a)) {
            return new Accession(a, Source.UNIPROT);
        }
        a = obtainGenbankAccessorFromDataFields(n);
        if (!ForesterUtil.isEmpty(a)) {
            return new Accession(a, Source.NCBI);
        }
        a = obtainRefSeqAccessorFromDataFields(n);
        if (!ForesterUtil.isEmpty(a)) {
            return new Accession(a, Source.REFSEQ);
        }
        a = obtainGiNumberFromDataFields(n);
        if (!ForesterUtil.isEmpty(a)) {
            return new Accession(a, Source.GI);
        }
        return null;
    }

    public final static Accession obtainFromSeqAccession(final PhylogenyNode n) {
        if (n.getNodeData().isHasSequence() && (n.getNodeData().getSequence().getAccession() != null)
                && !ForesterUtil.isEmpty(n.getNodeData().getSequence().getAccession().getSource())
                && !ForesterUtil.isEmpty(n.getNodeData().getSequence().getAccession().getValue())) {
            final String source = n.getNodeData().getSequence().getAccession().getSource().toLowerCase();
            final String value = n.getNodeData().getSequence().getAccession().getValue();
            if ((source.startsWith("uniprot") || source.equals("swissprot") || source.equals("trembl")
                    || source.equals("sp"))) {
                return new Accession(value, Source.UNIPROT);
            } else if (source.equals("embl") || source.equals("ebi")) {
                return new Accession(value, Source.EMBL);
            } else if (source.equals("ncbi") || source.equals("genbank")) {
                return new Accession(value, Source.NCBI);
            } else if (source.equals("refseq")) {
                return new Accession(value, Source.REFSEQ);
            } else if (source.equals("gi")) {
                return new Accession(value, Source.GI);
            }
        }
        return null;
    }

    public final static String obtainGenbankAccessorFromDataFields(final PhylogenyNode n) {
        String a = null;
        if (n.getNodeData().isHasSequence()) {
            final Sequence seq = n.getNodeData().getSequence();
            if (!ForesterUtil.isEmpty(seq.getSymbol())) {
                a = parseGenbankAccessorFromString(seq.getSymbol());
            }
            if (!ForesterUtil.isEmpty(seq.getGeneName())) {
                a = parseGenbankAccessorFromString(seq.getGeneName());
            }
            if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getName())) {
                a = parseGenbankAccessorFromString(seq.getName());
            }
            if (ForesterUtil.isEmpty(a) && (n.getNodeData().getSequence().getAccession() != null)
                    && !ForesterUtil.isEmpty(seq.getAccession().getValue())) {
                a = parseGenbankAccessorFromString(seq.getAccession().getValue());
            }
        }
        if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(n.getName())) {
            a = parseGenbankAccessorFromString(n.getName());
        }
        return a;
    }

    public final static String obtainGiNumberFromDataFields(final PhylogenyNode n) {
        String a = null;
        if (n.getNodeData().isHasSequence()) {
            final Sequence seq = n.getNodeData().getSequence();
            if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getName())) {
                a = parseGInumberFromString(seq.getName());
            }
            if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getGeneName())) {
                a = parseGInumberFromString(seq.getGeneName());
            }
            if (ForesterUtil.isEmpty(a) && (n.getNodeData().getSequence().getAccession() != null)
                    && !ForesterUtil.isEmpty(seq.getAccession().getValue())) {
                a = parseGInumberFromString(seq.getAccession().getValue());
            }
        }
        if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(n.getName())) {
            a = parseGInumberFromString(n.getName());
        }
        return a;
    }

    public final static String obtainRefSeqAccessorFromDataFields(final PhylogenyNode n) {
        String a = null;
        if (n.getNodeData().isHasSequence()) {
            final Sequence seq = n.getNodeData().getSequence();
            if (!ForesterUtil.isEmpty(seq.getSymbol())) {
                a = parseRefSeqAccessorFromString(seq.getSymbol());
            }
            if (!ForesterUtil.isEmpty(seq.getGeneName())) {
                a = parseRefSeqAccessorFromString(seq.getGeneName());
            }
            if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getName())) {
                a = parseRefSeqAccessorFromString(seq.getName());
            }
            if (ForesterUtil.isEmpty(a) && (n.getNodeData().getSequence().getAccession() != null)
                    && !ForesterUtil.isEmpty(seq.getAccession().getValue())) {
                a = parseRefSeqAccessorFromString(seq.getAccession().getValue());
            }
        }
        if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(n.getName())) {
            a = parseRefSeqAccessorFromString(n.getName());
        }
        return a;
    }

    public final static String obtainUniProtAccessorFromDataFields(final PhylogenyNode n) {
        String a = "";
        if (n.getNodeData().isHasSequence()) {
            final Sequence seq = n.getNodeData().getSequence();
            if (!ForesterUtil.isEmpty(seq.getSymbol())) {
                a = SequenceAccessionTools.parseUniProtAccessorFromString(seq.getSymbol());
            }
            if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getName())) {
                a = SequenceAccessionTools.parseUniProtAccessorFromString(seq.getName());
            }
            if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(seq.getGeneName())) {
                a = SequenceAccessionTools.parseUniProtAccessorFromString(seq.getGeneName());
            }
            if (ForesterUtil.isEmpty(a) && (n.getNodeData().getSequence().getAccession() != null)
                    && !ForesterUtil.isEmpty(seq.getAccession().getValue())) {
                a = SequenceAccessionTools.parseUniProtAccessorFromString(seq.getAccession().getValue());
            }
        }
        if (ForesterUtil.isEmpty(a) && !ForesterUtil.isEmpty(n.getName())) {
            a = SequenceAccessionTools.parseUniProtAccessorFromString(n.getName());
        }
        return a;
    }

    public final static Accession parseAccessorFromString(final String s) {
        if (!ForesterUtil.isEmpty(s)) {
            String v = parseGenbankAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.NCBI);
            }
            v = parseRefSeqAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.REFSEQ);
            }
            v = parseGInumberFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.GI);
            }
            v = parseEnsemlAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.ENSEMBL);
            }
            v = parseUniProtAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.UNIPROT);
            }
        }
        return null;
    }

    public final static Accession parseAccessorFromString_GenbankProteinPriority(final String s) {
        if (!ForesterUtil.isEmpty(s)) {
            String v = parseGenbankProteinAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.NCBI);
            }
            v = parseRefSeqAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.REFSEQ);
            }
            v = parseGInumberFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.GI);
            }
            v = parseEnsemlAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.ENSEMBL);
            }
            v = parseUniProtAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.UNIPROT);
            }
            v = parseGenbankAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.NCBI);
            }
        }
        return null;
    }

    public final static Accession parseAccessorFromString_UniProtPriority(final String s) {
        if (!ForesterUtil.isEmpty(s)) {
            String v = parseUniProtAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.UNIPROT);
            }
            v = parseRefSeqAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.REFSEQ);
            }
            v = parseGInumberFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.GI);
            }
            v = parseEnsemlAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.ENSEMBL);
            }
            v = parseGenbankAccessorFromString(s);
            if (!ForesterUtil.isEmpty(v)) {
                return new Accession(v, Source.NCBI);
            }
        }
        return null;
    }

    public final static String parseGenbankAccessorFromString(final String s) {
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // 2020-05-08:
        // Changed order: before it was GENBANK_NUC_PATTERN_1 first, and GENBANK_PROT_PATTERN last.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Matcher m = GENBANK_PROT_PATTERN.matcher(s);
        if (m.lookingAt()) {
            return m.group(1);
        } else {
            m = GENBANK_NUC_PATTERN_2.matcher(s);
            if (m.lookingAt()) {
                return m.group(1);
            } else {
                m = GENBANK_NUC_PATTERN_1.matcher(s);
                if (m.lookingAt()) {
                    return m.group(1);
                } else {
                    return null;
                }
            }
        }
    }

    public final static String parseGenbankProteinAccessorFromString(final String s) {
        final Matcher m = GENBANK_PROT_PATTERN.matcher(s);
        if (m.lookingAt()) {
            return m.group(1);
        } else {
            return null;
        }
    }

    public final static String parseGInumberFromString(final String s) {
        final Matcher m = GI_PATTERN.matcher(s);
        if (m.find()) {
            return m.group(1);
        }
        return null;
    }

    public final static String parseEnsemlAccessorFromString(final String s) {
        final Matcher m = ENSEMBL_PATTERN.matcher(s);
        if (m.find()) {
            return m.group(1);
        }
        return null;
    }

    public final static String parseRefSeqAccessorFromString(final String s) {
        final Matcher m = REFSEQ_PATTERN.matcher(s);
        if (m.lookingAt()) {
            return m.group(1);
        }
        return null;
    }

    public final static String parseUniProtAccessorFromString(final String s) {
        Matcher m = UNIPROT_KB_PATTERN_1.matcher(s);
        if (m.find()) {
            return m.group(1);
        }
        m = UNIPROT_KB_PATTERN_2.matcher(s);
        if (m.find()) {
            return m.group();
        }
        m = UNIPROT_KB_PATTERN_0.matcher(s);
        if (m.find()) {
            return m.group(1);
        }
        return null;
    }

    public static void main(String[] args) {
        System.out.println(parseGenbankAccessorFromString("LT898418.1 Porcine epidemic diarrhea virus isolate PEDV_AUSTRIA_L01065-M10_15-04_2015 genome assembly, complete genome: monopartite"));
    }

}
