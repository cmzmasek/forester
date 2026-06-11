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

package org.forester.clade_analysis;

import java.io.File;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterUtil;

public class CladeAnalysisDemo {

    private final static String PATH_TO_TEST_DATA = System.getProperty("user.dir") + ForesterUtil.getFileSeparator()
            + "test_data" + ForesterUtil.getFileSeparator();

    public static void main(final String[] args) {
        boolean failed = false;

        if (!testCladeAnalysis1()) {
            System.out.println("Demo 1 failed");
            failed = true;
        }
        if (!testCladeAnalysis2()) {
            System.out.println("Demo 2 failed");
            failed = true;
        }

        if (!testCladeAnalysis3()) {
            System.out.println("Demo 3 failed");
            failed = true;
        }

        if (!testCladeAnalysis4()) {
            System.out.println("Demo 4 failed");
            failed = true;
        }

        if (!testCladeAnalysis5()) {
            System.out.println("Demo 5 failed");
            failed = true;
        }

        if (!testCladeAnalysis6()) {
            System.out.println("Demo 6 failed");
            failed = true;
        }

        if (!testCladeAnalysis7()) {
            System.out.println("Demo 7 failed");
            failed = true;
        }

        if (!testCladeAnalysis8()) {
            System.out.println("Demo 8 failed");
            failed = true;
        }

        if (!testCladeAnalysis9()) {
            System.out.println("Demo 9 failed");
            failed = true;
        }


        if (!failed) {
            System.out.println("OK");
        } else {
            System.out.println("NOT OK");
        }
    }


    private static boolean testCladeAnalysis1() {
        try {
            final File in = new File(PATH_TO_TEST_DATA + "cladinator_demo_1.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(in, true);
            final Phylogeny p1 = factory.create(in, pp)[0];
            ResultMulti res = AnalysisMulti.execute(p1);

            System.out.println("DEMO 1:");
            System.out.println("+++++++");
            System.out.print(res.toString());
            System.out.println("------------------------- ");
            System.out.println();
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis2() {
        try {
            final File in = new File(PATH_TO_TEST_DATA + "cladinator_demo_2.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(in, true);
            final Phylogeny p1 = factory.create(in, pp)[0];
            ResultMulti res = AnalysisMulti.execute(p1);

            System.out.println("DEMO 2:");
            System.out.println("+++++++");
            System.out.print(res.toString());
            System.out.println("------------------------- ");
            System.out.println();
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis3() {
        try {
            final File in = new File(PATH_TO_TEST_DATA + "cladinator_demo_3.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(in, true);
            final Phylogeny p1 = factory.create(in, pp)[0];
            ResultMulti res = AnalysisMulti.execute(p1);

            System.out.println("DEMO 3:");
            System.out.println("+++++++");
            System.out.print(res.toString());
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
            final File in = new File(PATH_TO_TEST_DATA + "cladinator_demo_4.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(in, true);
            final Phylogeny p1 = factory.create(in, pp)[0];
            ResultMulti res = AnalysisMulti.execute(p1);

            System.out.println("DEMO 4:");
            System.out.println("+++++++");
            System.out.print(res.toString());
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
            final File in = new File(PATH_TO_TEST_DATA + "cladinator_demo_5.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(in, true);
            final Phylogeny p1 = factory.create(in, pp)[0];
            ResultMulti res = AnalysisMulti.execute(p1);

            System.out.println("DEMO 5:");
            System.out.println("+++++++");
            System.out.print(res.toString());
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
            final File in = new File(PATH_TO_TEST_DATA + "cladinator_demo_6.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(in, true);
            final Phylogeny p1 = factory.create(in, pp)[0];
            ResultMulti res = AnalysisMulti.execute(p1);

            System.out.println("DEMO 6:");
            System.out.println("+++++++");
            System.out.print(res.toString());
            System.out.println("------------------------- ");
            System.out.println();
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis7() {
        try {
            final File in = new File(PATH_TO_TEST_DATA + "cladinator_demo_7.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(in, true);
            final Phylogeny p1 = factory.create(in, pp)[0];
            ResultMulti res = AnalysisMulti.execute(p1);

            System.out.println("DEMO 7:");
            System.out.println("+++++++");
            System.out.print(res.toString());
            System.out.println("------------------------- ");
            System.out.println();
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis8() {
        try {
            final File in = new File(PATH_TO_TEST_DATA + "cladinator_demo_8.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(in, true);
            final Phylogeny p1 = factory.create(in, pp)[0];
            ResultMulti res = AnalysisMulti.execute(p1);

            System.out.println("DEMO 8:");
            System.out.println("+++++++");
            System.out.print(res.toString());
            System.out.println("------------------------- ");
            System.out.println();
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis9() {
        try {
            final File in = new File(PATH_TO_TEST_DATA + "cladinator_demo_9.xml");
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(in, true);
            final Phylogeny p1 = factory.create(in, pp)[0];
            ResultMulti res = AnalysisMulti.execute(p1);

            System.out.println("DEMO 9:");
            System.out.println("+++++++");
            System.out.print(res.toString());
            System.out.println("------------------------- ");
            System.out.println();
        } catch (final Exception e) {
            e.printStackTrace(System.out);
            return false;
        }
        return true;
    }


}
