
package org.forester.rio;

import org.forester.datastructures.IntMatrix;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.PhylogenyNodeField;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.rio.RIO.REROOTING;
import org.forester.sdi.SDIutil.ALGORITHM;
import org.forester.sdi.SDIutil.TaxonomyComparisonBase;
import org.forester.util.ForesterUtil;

public final class TestRIO {

    private final static String PATH_TO_TEST_DATA = System.getProperty( "user.dir" ) + ForesterUtil.getFileSeparator()
                                                          + "test_data" + ForesterUtil.getFileSeparator();

    public static boolean test() {
        if ( !testRIO_GSDIR() ) {
            return false;
        }
        return true;
    }

    private static boolean testRIO_GSDIR() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final NHXParser nhx = new NHXParser();
            nhx.setReplaceUnderscores( false );
            nhx.setIgnoreQuotes( true );
            nhx.setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.YES );
            final String gene_trees_1_str = "(((((MOUSE,RAT),HUMAN),CAEEL),YEAST),ARATH);"
                    + "((((MOUSE,RAT),HUMAN),(ARATH,YEAST)),CAEEL);" + "((MOUSE,RAT),(((ARATH,YEAST),CAEEL),HUMAN));"
                    + "(((((MOUSE,HUMAN),RAT),CAEEL),YEAST),ARATH);" + "((((HUMAN,MOUSE),RAT),(ARATH,YEAST)),CAEEL);";
            final Phylogeny[] gene_trees_1 = factory.create( gene_trees_1_str, nhx );
            final String species_trees_1_str = "(((((MOUSE,RAT),HUMAN),CAEEL),YEAST),ARATH);";
            final Phylogeny species_tree_1 = factory.create( species_trees_1_str, new NHXParser() )[ 0 ];
            species_tree_1.setRooted( true );
            PhylogenyMethods.transferNodeNameToField( species_tree_1, PhylogenyNodeField.TAXONOMY_CODE, true );
            //Archaeopteryx.createApplication( species_trees_1 );
            RIO rio = RIO.executeAnalysis( gene_trees_1,
                                           species_tree_1,
                                           ALGORITHM.GSDIR,
                                           REROOTING.BY_ALGORITHM,
                                           "",
                                           true,
                                           false );
            if ( rio.getAnalyzedGeneTrees().length != 5 ) {
                return false;
            }
            if ( rio.getExtNodesOfAnalyzedGeneTrees() != 6 ) {
                return false;
            }
            if ( rio.getGSDIRtaxCompBase() != TaxonomyComparisonBase.CODE ) {
                return false;
            }
            if ( rio.getRemovedGeneTreeNodes().size() != 0 ) {
                return false;
            }
            IntMatrix m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            // System.out.println( m.toString() );
            if ( !m.getRowAsString( 0, ',' ).equals( "ARATH,5,5,5,5,5,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 1, ',' ).equals( "CAEEL,5,5,5,5,5,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 2, ',' ).equals( "HUMAN,5,5,5,5,3,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 3, ',' ).equals( "MOUSE,5,5,5,5,3,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 4, ',' ).equals( "RAT,5,5,3,3,5,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 5, ',' ).equals( "YEAST,5,5,5,5,5,5" ) ) {
                return false;
            }
            //
            final Phylogeny[] gene_trees_2 = factory.create( gene_trees_1_str, nhx );
            final String species_trees_2_str = "((((MOUSE,RAT,HUMAN),CAEEL),YEAST),ARATH);";
            final Phylogeny species_tree_2 = factory.create( species_trees_2_str, new NHXParser() )[ 0 ];
            species_tree_2.setRooted( true );
            PhylogenyMethods.transferNodeNameToField( species_tree_2, PhylogenyNodeField.TAXONOMY_CODE, true );
            rio = RIO.executeAnalysis( gene_trees_2, species_tree_2 );
            m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            // System.out.println( m.toString() );
            if ( !m.getRowAsString( 0, ',' ).equals( "ARATH,5,5,5,5,5,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 1, ',' ).equals( "CAEEL,5,5,5,5,5,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 2, ',' ).equals( "HUMAN,5,5,5,5,5,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 3, ',' ).equals( "MOUSE,5,5,5,5,5,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 4, ',' ).equals( "RAT,5,5,5,5,5,5" ) ) {
                return false;
            }
            if ( !m.getRowAsString( 5, ',' ).equals( "YEAST,5,5,5,5,5,5" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    public static void main( final String[] args ) {
        if ( !testRIO_GSDIR() ) {
            System.out.println( "testRIO GSDIR failed" );
        }
        else {
            System.out.println( "OK" );
        }
    }
}