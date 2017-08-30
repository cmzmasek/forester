
package org.forester.clade_analysis;

import java.io.File;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterUtil;

public class CladeAnalysisDemo {

    private final static String PATH_TO_TEST_DATA = System.getProperty( "user.dir" ) + ForesterUtil.getFileSeparator()
            + "test_data" + ForesterUtil.getFileSeparator();

    public static void main( final String[] args ) {
        boolean failed = false;
       
        if ( !testCladeAnalysis1() ) {
            System.out.println( "Demo 1 failed" );
            failed = true;
        }
        if ( !testCladeAnalysis2() ) {
            System.out.println( "Demo 2 failed" );
            failed = true;
        }
        
        if ( !testCladeAnalysis3() ) {
            System.out.println( "Demo 3 failed" );
            failed = true;
        }
        
        if ( !testCladeAnalysis4() ) {
            System.out.println( "Demo 4 failed" );
            failed = true;
        }
        
        if ( !testCladeAnalysis5() ) {
            System.out.println( "Demo 5 failed" );
            failed = true;
        }
        
        if ( !testCladeAnalysis6() ) {
            System.out.println( "Demo 6 failed" );
            failed = true;
        }
        
        if ( !testCladeAnalysis7() ) {
            System.out.println( "Demo 7 failed" );
            failed = true;
        }
        
        if ( !testCladeAnalysis8() ) {
            System.out.println( "Demo 8 failed" );
            failed = true;
        }
        
        if ( !testCladeAnalysis9() ) {
            System.out.println( "Demo 9 failed" );
            failed = true;
        }
     
       
        if ( !failed ) {
            System.out.println( "OK" );
        }
        else {
            System.out.println( "NOT OK" );
        }
    }


    
    

    

    private static boolean testCladeAnalysis1() {
        try {
            final File in = new File( PATH_TO_TEST_DATA + "cladinator_demo_1.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny p1 = factory.create( in, pp )[ 0 ];
            ResultMulti res = AnalysisMulti.execute( p1, 0.5 );
            
            System.out.println( "DEMO 1:" );
            System.out.println( "+++++++" );
            System.out.print( res.toString() );
            System.out.println( "------------------------- " );
            System.out.println();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
    
    private static boolean testCladeAnalysis2() {
        try {
            final File in = new File( PATH_TO_TEST_DATA + "cladinator_demo_2.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny p1 = factory.create( in, pp )[ 0 ];
            ResultMulti res = AnalysisMulti.execute( p1, 0.5 );
            
            System.out.println( "DEMO 2:" );
            System.out.println( "+++++++" );
            System.out.print( res.toString() );
            System.out.println( "------------------------- " );
            System.out.println();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
    
    private static boolean testCladeAnalysis3() {
        try {
            final File in = new File( PATH_TO_TEST_DATA + "cladinator_demo_3.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny p1 = factory.create( in, pp )[ 0 ];
            ResultMulti res = AnalysisMulti.execute( p1, 0.5 );
            
            System.out.println( "DEMO 3:" );
            System.out.println( "+++++++" );
            System.out.print( res.toString() );
            System.out.println( "------------------------- " );
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
            final File in = new File( PATH_TO_TEST_DATA + "cladinator_demo_4.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny p1 = factory.create( in, pp )[ 0 ];
            ResultMulti res = AnalysisMulti.execute( p1, 0.5 );
            
            System.out.println( "DEMO 4:" );
            System.out.println( "+++++++" );
            System.out.print( res.toString() );
            System.out.println( "------------------------- " );
            System.out.println();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
    
    private static boolean testCladeAnalysis5() {
        try {
            final File in = new File( PATH_TO_TEST_DATA + "cladinator_demo_5.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny p1 = factory.create( in, pp )[ 0 ];
            ResultMulti res = AnalysisMulti.execute( p1, 0.5 );
            
            System.out.println( "DEMO 5:" );
            System.out.println( "+++++++" );
            System.out.print( res.toString() );
            System.out.println( "------------------------- " );
            System.out.println();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
    
    private static boolean testCladeAnalysis6() {
        try {
            final File in = new File( PATH_TO_TEST_DATA + "cladinator_demo_6.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny p1 = factory.create( in, pp )[ 0 ];
            ResultMulti res = AnalysisMulti.execute( p1, 0.5 );
            
            System.out.println( "DEMO 6:" );
            System.out.println( "+++++++" );
            System.out.print( res.toString() );
            System.out.println( "------------------------- " );
            System.out.println();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
    
    private static boolean testCladeAnalysis7() {
        try {
            final File in = new File( PATH_TO_TEST_DATA + "cladinator_demo_7.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny p1 = factory.create( in, pp )[ 0 ];
            ResultMulti res = AnalysisMulti.execute( p1, 0.5 );
            
            System.out.println( "DEMO 7:" );
            System.out.println( "+++++++" );
            System.out.print( res.toString() );
            System.out.println( "------------------------- " );
            System.out.println();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
    
    private static boolean testCladeAnalysis8() {
        try {
            final File in = new File( PATH_TO_TEST_DATA + "cladinator_demo_8.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny p1 = factory.create( in, pp )[ 0 ];
            ResultMulti res = AnalysisMulti.execute( p1, 0.5 );
            
            System.out.println( "DEMO 8:" );
            System.out.println( "+++++++" );
            System.out.print( res.toString() );
            System.out.println( "------------------------- " );
            System.out.println();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
    
    private static boolean testCladeAnalysis9() {
        try {
            final File in = new File( PATH_TO_TEST_DATA + "cladinator_demo_9.xml" );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny p1 = factory.create( in, pp )[ 0 ];
            ResultMulti res = AnalysisMulti.execute( p1, 0.5 );
            
            System.out.println( "DEMO 9:" );
            System.out.println( "+++++++" );
            System.out.print( res.toString() );
            System.out.println( "------------------------- " );
            System.out.println();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

 
}
