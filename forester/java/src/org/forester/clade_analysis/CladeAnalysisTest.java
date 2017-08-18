
package org.forester.clade_analysis;

import java.io.File;
import java.util.List;
import java.util.regex.Pattern;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterUtil;

public class CladeAnalysisTest {

    private final static String PATH_TO_TEST_DATA = System.getProperty( "user.dir" ) + ForesterUtil.getFileSeparator()
            + "test_data" + ForesterUtil.getFileSeparator();

    public static void main( final String[] args ) {
        boolean failed = false;
        if ( !testCladeAnalysis1() ) {
            System.out.println( "Clade analysis 1 failed" );
            failed = true;
        }
        if ( !testCladeAnalysis2() ) {
            System.out.println( "Clade analysis 2 failed" );
            failed = true;
        }
        if ( !testCladeAnalysis3() ) {
            System.out.println( "Clade analysis 3 failed" );
            failed = true;
        }
        if ( !testCladeAnalysis4() ) {
            System.out.println( "Clade analysis 3 failed" );
            failed = true;
        }
        if ( !failed ) {
            System.out.println( "OK" );
        }
    }

    public static boolean test() {
        if ( !testCladeAnalysis1() ) {
            return false;
        }
        if ( !testCladeAnalysis2() ) {
            return false;
        }
        if ( !testCladeAnalysis3() ) {
            return false;
        }
        if ( !testCladeAnalysis4() ) {
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis1() {
        try {
            final File intreefile1 = new File( PATH_TO_TEST_DATA + "clade_analysis_test_1.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intreefile1, true );
            final Phylogeny p1 = factory.create( intreefile1, pp )[ 0 ];
            Result res = Analysis.execute( p1, "A.1.1.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.1.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 4 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.1.1.2", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.1.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 4 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.1.1.3", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.1.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 4 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.1.1.4", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.1.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 3 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.1.2.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 17 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.2.1.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.2.1.2" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 17 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.2.1.2", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.2.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 17 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.3.1.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A.3" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.3.1.2" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.3.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 2 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.3.1.2", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A.3" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.3.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.3.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 2 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.3.2.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A.3" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.3.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.3.3.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 3 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.3.3.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.3" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 10 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.4.1.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A.4.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.4.1.1.a" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.4.1.2" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 3 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.4.1.1.a", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A.4.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.4.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.4.1.2" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 3 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.4.1.2", "." );
            res = Analysis.execute( p1, "A.4.1.2.a", "." );
            res = Analysis.execute( p1, "A.5.1.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.5.1.2" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 10 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.5.1.2", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.5.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 10 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "A.6.3.12", "." );
            if ( !res.getGreatestCommonPrefix().equals( "A" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 17 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "B.1.1.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "B.1.234.3" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 25 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 2 ) {
                return false;
            }
            res = Analysis.execute( p1, "B.1.234.3", "." );
            if ( !res.getGreatestCommonPrefix().equals( "" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "B.1.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 25 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 2 ) {
                return false;
            }
            res = Analysis.execute( p1, "C.1.1.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "C.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "C.1.1.2" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "C.1.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 2 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "C.1.1.2", "." );
            if ( !res.getGreatestCommonPrefix().equals( "C.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "C.1.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "C.1.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 2 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "C.1.2.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "C" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "C.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "C.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 3 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "C.2.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "C" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "C.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "C.3" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 4 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "C.3", "." );
            if ( !res.getGreatestCommonPrefix().equals( "" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "C" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "QE.1.1.1.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 5 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 1 ) {
                return false;
            }
            res = Analysis.execute( p1, "QE.1.1.1.2.1", "." );
            if ( !res.getGreatestCommonPrefix().equals( "" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "C" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 25 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 2 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis2() {
        try {
            final File intreefile1 = new File( PATH_TO_TEST_DATA + "clade_analysis_test_2.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intreefile1, true );
            final Phylogeny p1 = factory.create( intreefile1, pp )[ 0 ];
            Result res = Analysis.execute( p1, "6_DQ278891", null );
            if ( !res.getGreatestCommonPrefix().equals( "6_" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "6_DQ278893" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "6_JX183550" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 2 ) {
                return false;
            }
            if ( res.getTreeSize() != 219 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "6xa_EU408330", null );
            if ( !res.getGreatestCommonPrefix().equals( "6xa_EU40833" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "6xa_EU408331" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "6xa_EU408332" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 2 ) {
                return false;
            }
            if ( res.getTreeSize() != 219 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }
            res = Analysis.execute( p1, "7a_EF108306", null );
            if ( !res.getGreatestCommonPrefix().equals( "" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "2" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 219 ) {
                return false;
            }
            if ( res.getTreeSize() != 219 ) {
                return false;
            }
            if ( res.getWarnings().size() != 2 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testCladeAnalysis3() {
        try {
            final Result2 res1 = new Result2();
            res1.addGreatestCommonPrefix( "A.1.1", 0.3 );
            res1.addGreatestCommonPrefix( "A.1.2", 0.3 );
            res1.addGreatestCommonPrefix( "A.1.3", 0.3 );
            res1.addGreatestCommonPrefix( "B.1", 0.1 );
            res1.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res1.toString());
            System.out.println( "------------------------- ");
            System.out.println();
            
            final Result2 res2 = new Result2( "." );
            res2.addGreatestCommonPrefix( "A.1.1.1", 0.1 );
            res2.addGreatestCommonPrefix( "A.1", 0.7 );
            res2.addGreatestCommonPrefix( "A.1.2", 0.1 );
            res2.addGreatestCommonPrefix( "B.1", 0.1 );
            res2.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res2.toString());
            System.out.println( "------------------------- ");
            System.out.println();
           
            final Result2 res3 = new Result2( "." );
            res3.addGreatestCommonPrefix( "A.1.1.1", 0.1 );
            res3.addGreatestCommonPrefix( "A.1.1.1.1", 0.6 );
            res3.addGreatestCommonPrefix( "A.1", 0.1 );
            res3.addGreatestCommonPrefix( "A.1.2", 0.1 );
            res3.addGreatestCommonPrefix( "B.1", 0.1 );
            res3.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res3.toString());
            System.out.println( "------------------------- ");
            System.out.println();
            
            final Result2 res33 = new Result2( "." );
            res33.addGreatestCommonPrefix( "A.1.1.1", 0.1 );
            res33.addGreatestCommonPrefix( "A.1.1.1.1", 0.3 );
            res33.addGreatestCommonPrefix( "A.1", 0.1 );
            res33.addGreatestCommonPrefix( "A.1.2", 0.1 );
            res33.addGreatestCommonPrefix( "B.1", 0.1 );
            res33.addGreatestCommonPrefix( "B.1.1.1", 0.3 );
            res33.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res33.toString());
            System.out.println( "------------------------- ");
            System.out.println();
            
            final Result2 res4 = new Result2();
            res4.addGreatestCommonPrefix( "A.1.1.1.1", 0.35 );
            res4.addGreatestCommonPrefix( "A.1.1.1.2", 0.35 );
            res4.addGreatestCommonPrefix( "A.1", 0.1 );
            res4.addGreatestCommonPrefix( "A.1.2", 0.1 );
            res4.addGreatestCommonPrefix( "B.1", 0.1 );
            res4.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res4.toString());
            System.out.println( "------------------------- ");
            System.out.println();
            
            final Result2 res5 = new Result2();
            res5.addGreatestCommonPrefix( "A.1.1.1.1", 0.2 );
            res5.addGreatestCommonPrefix( "C.2.3", 0.2 );
            res5.addGreatestCommonPrefix( "A.1.5", 0.1 );
            res5.addGreatestCommonPrefix( "A.3.1.4", 0.2 );
            res5.addGreatestCommonPrefix( "B.1.1", 0.2 );
            res5.addGreatestCommonPrefix( "B.1.2", 0.09 );
            res5.addGreatestCommonPrefix( "D.1.1.1.1", 0.01 );
            res5.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res5.toString());
            System.out.println( "------------------------- ");
            System.out.println();
           
            final Result2 res6 = new Result2();
            res6.addGreatestCommonPrefix( "A.1.1.1", 0.05 );
            res6.addGreatestCommonPrefix( "A.1.1.1.1", 0.65 );
            res6.addGreatestCommonPrefix( "A.1", 0.1 );
            res6.addGreatestCommonPrefix( "A.1.2", 0.1 );
            res6.addGreatestCommonPrefix( "B.1", 0.1 );
            res6.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res6.toString());
            System.out.println( "------------------------- ");
            System.out.println();
           
            final Result2 res7 = new Result2();
            res7.addGreatestCommonPrefix( "A.1.1.1", 0.07 );
            res7.addGreatestCommonPrefix( "A.1.1.1.1", 0.9 );
            res7.addGreatestCommonPrefix( "A.1", 0.01 );
            res7.addGreatestCommonPrefix( "A.1.2", 0.01 );
            res7.addGreatestCommonPrefix( "B.1", 0.01 );
            res7.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res7.toString());
            System.out.println( "------------------------- ");
            System.out.println();
            
            final Result2 res8 = new Result2( "_/_" );
            res8.addGreatestCommonPrefix( "AA_/_abc_/_def", 0.07 );
            res8.addGreatestCommonPrefix( "AA_/_abc_/_sfc", 0.9 );
            res8.addGreatestCommonPrefix( "AA_/_abc_/_xcd", 0.01 );
            res8.addGreatestCommonPrefix( "AA_/_abc_/_memr", 0.01 );
            res8.addGreatestCommonPrefix( "AA_/_abc_/_fkem_/_odem", 0.01 );
            res8.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res8.toString());
            System.out.println( "------------------------- ");
            System.out.println();
           
            final Result2 res9 = new Result2( "_/_" );
            res9.addGreatestCommonPrefix( "AA_/_abc_/_def", 0.07 );
            res9.addGreatestCommonPrefix( "AA_/_abc_/_sfc", 0.6 );
            res9.addGreatestCommonPrefix( "AA_/_abc_/_xcd", 0.01 );
            res9.addGreatestCommonPrefix( "AA_/_abc_/_memr", 0.01 );
            res9.addGreatestCommonPrefix( "AA_/_abc_/_fkem_/_odem", 0.01 );
            res9.addGreatestCommonPrefix( "BB_/_fke_/_dme_/_nx2", 0.3 );
            res9.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res9.toString());
            System.out.println( "------------------------- ");
            System.out.println();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
    
    private static boolean testCladeAnalysis4() {
        try {
            final File intreefile1 = new File( PATH_TO_TEST_DATA + "pplacer_2.tre" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intreefile1, true );
            final Phylogeny p1 = factory.create( intreefile1, pp )[ 0 ];
            Pattern query = Pattern.compile(".+#\\d+_M=(.+)");
            Result2 res = Analysis2.execute( p1, query, "." );
            
            res.analyzeGreatestCommonPrefixes( 0.3 );
            System.out.print( res.toString());
            System.out.println( "------------------------- ");
            System.out.println();
            
           // Result res = Analysis.execute( p1, "A.1.1.1", "." );
           /* if ( !res.getGreatestCommonPrefix().equals( "A.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixDown().equals( "A.1.1" ) ) {
                return false;
            }
            if ( !res.getGreatestCommonPrefixUp().equals( "A.1.2.1" ) ) {
                return false;
            }
            if ( res.getLeastEncompassingCladeSize() != 4 ) {
                return false;
            }
            if ( res.getTreeSize() != 25 ) {
                return false;
            }
            if ( res.getWarnings().size() != 0 ) {
                return false;
            }*/
          
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
}
