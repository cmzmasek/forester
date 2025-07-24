
package org.forester.clade_analysis;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterUtil;

import java.io.File;
import java.util.regex.Pattern;

public class CladeAnalysisTest {

    private final static String PATH_TO_TEST_DATA = System.getProperty("user.dir") + ForesterUtil.getFileSeparator()
            + "forester/test_data" + ForesterUtil.getFileSeparator();

    public static void main(final String[] args) {
        boolean failed = false;
        if (!testCladeAnalysis1()) {
            System.out.println("Clade analysis 1 failed");
            failed = true;
        }
        if (!testCladeAnalysis2()) {
            System.out.println("Clade analysis 2 failed");
            failed = true;
        }
        if (!testCladeAnalysis3()) {
            System.out.println("Clade analysis 3 failed");
            failed = true;
        }
        if (!testCladeAnalysis4()) {
            System.out.println("Clade analysis 4 failed");
            failed = true;
        }
        if (!testCladeAnalysis5()) {
            System.out.println("Clade analysis 5 failed");
            failed = true;
        }
        if (!testCladeAnalysis6()) {
            System.out.println("Clade analysis 6 failed");
            failed = true;
        }
        if (!failed) {
            System.out.println("OK");
        } else {
            System.out.println("NOT OK");
        }
    }

    public static boolean test() {
        if (!testCladeAnalysis1()) {
            return false;
        }
        if (!testCladeAnalysis2()) {
            return false;
        }
        if (!testCladeAnalysis3()) {
            return false;
        }
        if (!testCladeAnalysis4()) {
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis1() {
        try {
            final File intreefile1 = new File(PATH_TO_TEST_DATA + "clade_analysis_test_1.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(intreefile1, true);
            final Phylogeny p1 = factory.create(intreefile1, pp)[0];

            ResultSingle res = AnalysisSingle.execute(p1, "A.1.1.1", ".");
            if (!res.getGreatestCommonPrefix().equals("A.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A.1.2.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 4) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.1.1.2", ".");
            if (!res.getGreatestCommonPrefix().equals("A.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A.1.2.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 4) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.1.1.3", ".");
            if (!res.getGreatestCommonPrefix().equals("A.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A.1.2.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 4) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.1.1.4", ".");
            if (!res.getGreatestCommonPrefix().equals("A.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A.1.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 3) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.1.2.1", ".");
            if (!res.getGreatestCommonPrefix().equals("A")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 17) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.2.1.1", ".");
            if (!res.getGreatestCommonPrefix().equals("A")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.2.1.2")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 17) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.2.1.2", ".");
            if (!res.getGreatestCommonPrefix().equals("A")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.2.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 17) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.3.1.1", ".");
            if (!res.getGreatestCommonPrefix().equals("A.3")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.3.1.2")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A.3.2.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 2) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.3.1.2", ".");
            if (!res.getGreatestCommonPrefix().equals("A.3")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.3.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A.3.2.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 2) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.3.2.1", ".");
            if (!res.getGreatestCommonPrefix().equals("A.3")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.3.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A.3.3.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 3) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.3.3.1", ".");
            if (!res.getGreatestCommonPrefix().equals("A")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.3")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 10) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.4.1.1", ".");
            if (!res.getGreatestCommonPrefix().equals("A.4.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.4.1.1.a")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A.4.1.2")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 3) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.4.1.1.a", ".");
            if (!res.getGreatestCommonPrefix().equals("A.4.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.4.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A.4.1.2")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 3) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.4.1.2", ".");
            res = AnalysisSingle.execute(p1, "A.4.1.2.a", ".");
            res = AnalysisSingle.execute(p1, "A.5.1.1", ".");
            if (!res.getGreatestCommonPrefix().equals("A")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.5.1.2")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 10) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.5.1.2", ".");
            if (!res.getGreatestCommonPrefix().equals("A")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A.5.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 10) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "A.6.3.12", ".");
            if (!res.getGreatestCommonPrefix().equals("A")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("A")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("A")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 17) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "B.1.1.1", ".");
            if (!res.getGreatestCommonPrefix().equals("")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("B.1.234.3")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 25) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 2) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "B.1.234.3", ".");
            if (!res.getGreatestCommonPrefix().equals("")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("B.1.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 25) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 2) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "C.1.1.1", ".");
            if (!res.getGreatestCommonPrefix().equals("C.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("C.1.1.2")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("C.1.2.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 2) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "C.1.1.2", ".");
            if (!res.getGreatestCommonPrefix().equals("C.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("C.1.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("C.1.2.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 2) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "C.1.2.1", ".");
            if (!res.getGreatestCommonPrefix().equals("C")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("C.1.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("C.2.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 3) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "C.2.1", ".");
            if (!res.getGreatestCommonPrefix().equals("C")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("C.1")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("C.3")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 4) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "C.3", ".");
            if (!res.getGreatestCommonPrefix().equals("")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("C")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("QE.1.1.1.2.1")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 5) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 1) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "QE.1.1.1.2.1", ".");
            if (!res.getGreatestCommonPrefix().equals("")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("C")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 25) {
                return false;
            }
            if (res.getTreeSize() != 25) {
                return false;
            }
            if (res.getWarnings().size() != 2) {
                return false;
            }
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis2() {
        try {
            final File intreefile1 = new File(PATH_TO_TEST_DATA + "clade_analysis_test_2.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(intreefile1, true);
            final Phylogeny p1 = factory.create(intreefile1, pp)[0];
            ResultSingle res = AnalysisSingle.execute(p1, "6_DQ278891", null);
            if (!res.getGreatestCommonPrefix().equals("6_")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("6_DQ278893")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("6_JX183550")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 2) {
                return false;
            }
            if (res.getTreeSize() != 219) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "6xa_EU408330", null);
            if (!res.getGreatestCommonPrefix().equals("6xa_EU40833")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("6xa_EU408331")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("6xa_EU408332")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 2) {
                return false;
            }
            if (res.getTreeSize() != 219) {
                return false;
            }
            if (res.getWarnings().size() != 0) {
                return false;
            }
            res = AnalysisSingle.execute(p1, "7a_EF108306", null);
            if (!res.getGreatestCommonPrefix().equals("")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixDown().equals("2")) {
                return false;
            }
            if (!res.getGreatestCommonPrefixUp().equals("")) {
                return false;
            }
            if (res.getLeastEncompassingCladeSize() != 219) {
                return false;
            }
            if (res.getTreeSize() != 219) {
                return false;
            }
            if (res.getWarnings().size() != 2) {
                return false;
            }
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis3() {
        try {
            final ResultMulti res1 = new ResultMulti();
            res1.addGreatestCommonPrefix("A.1.1", 0.3);
            res1.addGreatestCommonPrefix("A.1.2", 0.3);
            res1.addGreatestCommonPrefix("A.1.3", 0.3);
            res1.addGreatestCommonPrefix("B.1", 0.1);
            res1.analyze();
            System.out.print(res1.toString());
            System.out.println("------------------------- ");
            System.out.println();
            final ResultMulti res2 = new ResultMulti(".");
            res2.addGreatestCommonPrefix("A.1.1.1", 0.1);
            res2.addGreatestCommonPrefix("A.1", 0.7);
            res2.addGreatestCommonPrefix("A.1.2", 0.1);
            res2.addGreatestCommonPrefix("B.1", 0.1);
            res2.analyze();
            System.out.print(res2.toString());
            System.out.println("------------------------- ");
            System.out.println();
            final ResultMulti res3 = new ResultMulti(".");
            res3.addGreatestCommonPrefix("A.1.1.1", 0.1);
            res3.addGreatestCommonPrefix("A.1.1.1.1", 0.6);
            res3.addGreatestCommonPrefix("A.1", 0.1);
            res3.addGreatestCommonPrefix("A.1.2", 0.1);
            res3.addGreatestCommonPrefix("B.1", 0.1);
            res3.analyze();
            System.out.print(res3.toString());
            System.out.println("------------------------- ");
            System.out.println();
            final ResultMulti res33 = new ResultMulti(".");
            res33.addGreatestCommonPrefix("A.1.1.1", 0.1);
            res33.addGreatestCommonPrefix("A.1.1.1.1", 0.3);
            res33.addGreatestCommonPrefix("A.1", 0.1);
            res33.addGreatestCommonPrefix("A.1.2", 0.1);
            res33.addGreatestCommonPrefix("B.1", 0.1);
            res33.addGreatestCommonPrefix("B.1.1.1", 0.3);
            res33.analyze();
            System.out.print(res33.toString());
            System.out.println("------------------------- ");
            System.out.println();
            final ResultMulti res4 = new ResultMulti();
            res4.addGreatestCommonPrefix("A.1.1.1.1", 0.35);
            res4.addGreatestCommonPrefix("A.1.1.1.2", 0.35);
            res4.addGreatestCommonPrefix("A.1", 0.1);
            res4.addGreatestCommonPrefix("A.1.2", 0.1);
            res4.addGreatestCommonPrefix("B.1", 0.1);
            res4.analyze();
            System.out.print(res4.toString());
            System.out.println("------------------------- ");
            System.out.println();
            final ResultMulti res5 = new ResultMulti();
            res5.addGreatestCommonPrefix("A.1.1.1.1", 0.2);
            res5.addGreatestCommonPrefix("C.2.3", 0.2);
            res5.addGreatestCommonPrefix("A.1.5", 0.1);
            res5.addGreatestCommonPrefix("A.3.1.4", 0.2);
            res5.addGreatestCommonPrefix("B.1.1", 0.2);
            res5.addGreatestCommonPrefix("B.1.2", 0.09);
            res5.addGreatestCommonPrefix("D.1.1.1.1", 0.01);
            res5.analyze();
            System.out.print(res5.toString());
            System.out.println("------------------------- ");
            System.out.println();
            final ResultMulti res6 = new ResultMulti();
            res6.addGreatestCommonPrefix("A.1.1.1", 0.05);
            res6.addGreatestCommonPrefix("A.1.1.1.1", 0.65);
            res6.addGreatestCommonPrefix("A.1", 0.1);
            res6.addGreatestCommonPrefix("A.1.2", 0.1);
            res6.addGreatestCommonPrefix("B.1", 0.1);
            res6.analyze();
            System.out.print(res6.toString());
            System.out.println("------------------------- ");
            System.out.println();
            final ResultMulti res7 = new ResultMulti();
            res7.addGreatestCommonPrefix("A.1.1.1", 0.07);
            res7.addGreatestCommonPrefix("A.1.1.1.1", 0.9);
            res7.addGreatestCommonPrefix("A.1", 0.01);
            res7.addGreatestCommonPrefix("A.1.2", 0.01);
            res7.addGreatestCommonPrefix("B.1", 0.01);
            res7.analyze();
            System.out.print(res7.toString());
            System.out.println("------------------------- ");
            System.out.println();
            final ResultMulti res8 = new ResultMulti("_/_");
            res8.addGreatestCommonPrefix("AA_/_abc_/_def", 0.07);
            res8.addGreatestCommonPrefix("AA_/_abc_/_sfc", 0.9);
            res8.addGreatestCommonPrefix("AA_/_abc_/_xcd", 0.01);
            res8.addGreatestCommonPrefix("AA_/_abc_/_memr", 0.01);
            res8.addGreatestCommonPrefix("AA_/_abc_/_fkem_/_odem", 0.01);
            res8.analyze();
            System.out.print(res8.toString());
            System.out.println("------------------------- ");
            System.out.println();
            final ResultMulti res9 = new ResultMulti("_/_");
            res9.addGreatestCommonPrefix("AA_/_abc_/_def", 0.07);
            res9.addGreatestCommonPrefix("AA_/_abc_/_sfc", 0.6);
            res9.addGreatestCommonPrefix("AA_/_abc_/_xcd", 0.01);
            res9.addGreatestCommonPrefix("AA_/_abc_/_memr", 0.01);
            res9.addGreatestCommonPrefix("AA_/_abc_/_fkem_/_odem", 0.01);
            res9.addGreatestCommonPrefix("BB_/_fke_/_dme_/_nx2", 0.3);
            res9.analyze();
            System.out.print(res9.toString());
            System.out.println("------------------------- ");
            System.out.println();
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis4() {
        try {
            final File intreefile1 = new File(PATH_TO_TEST_DATA + "pplacer_2.tre");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(intreefile1, true);
            final Phylogeny p1 = factory.create(intreefile1, pp)[0];
            final ResultMulti res2 = AnalysisMulti.execute(p1);
            res2.analyze();
            System.out.print(res2.toString());
            System.out.println("------------------------- ");
            System.out.println();
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis5() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final String t1s = "(((((A.1.1,Q_#1_M=1),A.1.2),(A.2.1,A.2.2)),((A.3.1,A.3.2),(A.4.1,A.4.2))),(((B.1,B.2),B.3),(C.1,C.2)))";
            final Phylogeny t1 = factory.create(t1s, new NHXParser())[0];
            final ResultMulti res1 = AnalysisMulti.execute(t1, ".");
            res1.analyze();
            System.out.print(res1.toString());
            System.out.println("------------------------- ");
            System.out.println();


        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis6() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final String t1s = "(((((A.1.1,A.1.2),Q_#0_M=0.5),((A.2.1,A.2.2),Q_#1_M=0.5)),((A.3.1,A.3.2),(A.4.1,A.4.2))),(((B.1,B.2),B.3),(C.1,C.2)))";
            final Phylogeny t1 = factory.create(t1s, new NHXParser())[0];
            final ResultMulti res1 = AnalysisMulti.execute(t1, ".");
            res1.analyze();
            System.out.print(res1.toString());
            System.out.println("------------------------- ");
            System.out.println();


        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }


}
