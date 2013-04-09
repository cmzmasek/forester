
package org.forester.rio;

import java.io.File;

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

    public static void main( final String[] args ) {
        if ( !testRIO_GSDIR() ) {
            System.out.println( "testRIO GSDIR failed" );
        }
        if ( !testRIO_GSDIR_Iterating() ) {
            System.out.println( "testRIO GSDIR iterating failed" );
        }
        else {
            System.out.println( "OK" );
        }
    }

    public static boolean test() {
        if ( !testRIO_GSDIR() ) {
            return false;
        }
        if ( !testRIO_GSDIR_Iterating() ) {
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
            nhx.setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            //
            final String gene_trees_00_str = "(MOUSE,RAT);(MOUSE,RAT);(MOUSE,RAT);(RAT,MOUSE);";
            final Phylogeny[] gene_trees_00 = factory.create( gene_trees_00_str, nhx );
            final String species_trees_00_str = "(MOUSE,RAT);";
            final Phylogeny species_tree_00 = factory.create( species_trees_00_str, new NHXParser() )[ 0 ];
            species_tree_00.setRooted( true );
            PhylogenyMethods.transferNodeNameToField( species_tree_00, PhylogenyNodeField.TAXONOMY_CODE, true );
            RIO rio = RIO.executeAnalysis( gene_trees_00,
                                           species_tree_00,
                                           ALGORITHM.GSDIR,
                                           REROOTING.BY_ALGORITHM,
                                           "",
                                           true,
                                           false,
                                           true );
            if ( rio.getAnalyzedGeneTrees().length != 4 ) {
                return false;
            }
            if ( rio.getExtNodesOfAnalyzedGeneTrees() != 2 ) {
                return false;
            }
            if ( rio.getGSDIRtaxCompBase() != TaxonomyComparisonBase.CODE ) {
                return false;
            }
            if ( rio.getRemovedGeneTreeNodes().size() != 0 ) {
                return false;
            }
            IntMatrix m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            if ( !m.getRowAsString( 0, ',' ).equals( "MOUSE,4,4" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            if ( !m.getRowAsString( 1, ',' ).equals( "RAT,4,4" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            final String gene_trees_000_str = "(MOUSE1[&&NHX:S=MOUSE],MOUSE2[&&NHX:S=MOUSE]);(MOUSE1[&&NHX:S=MOUSE],MOUSE2[&&NHX:S=MOUSE])";
            final Phylogeny[] gene_trees_000 = factory.create( gene_trees_000_str, nhx );
            final String species_trees_000_str = "[&&NHX:S=MOUSE];";
            final Phylogeny species_tree_000 = factory.create( species_trees_000_str, new NHXParser() )[ 0 ];
            species_tree_000.setRooted( true );
            rio = RIO.executeAnalysis( gene_trees_000,
                                       species_tree_000,
                                       ALGORITHM.GSDIR,
                                       REROOTING.BY_ALGORITHM,
                                       "",
                                       true,
                                       false,
                                       true );
            if ( rio.getAnalyzedGeneTrees().length != 2 ) {
                return false;
            }
            if ( rio.getExtNodesOfAnalyzedGeneTrees() != 2 ) {
                return false;
            }
            if ( rio.getGSDIRtaxCompBase() != TaxonomyComparisonBase.SCIENTIFIC_NAME ) {
                return false;
            }
            if ( rio.getRemovedGeneTreeNodes().size() != 0 ) {
                return false;
            }
            m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            if ( !m.getRowAsString( 0, ',' ).equals( "MOUSE1,2,0" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            if ( !m.getRowAsString( 1, ',' ).equals( "MOUSE2,0,2" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            //
            //
            final String gene_trees_0000_str = "(MOUSE1[&&NHX:S=MOUSE],MOUSE2[&&NHX:S=MOUSE]);(MOUSE1[&&NHX:S=MOUSE],MOUSE2[&&NHX:S=MOUSE]);(MOUSE1[&&NHX:S=MOUSE],MOUSE2[&&NHX:S=MOUSE])";
            final Phylogeny[] gene_trees_0000 = factory.create( gene_trees_0000_str, nhx );
            final String species_trees_0000_str = "([&&NHX:S=MOUSE]);";
            final Phylogeny species_tree_0000 = factory.create( species_trees_0000_str, new NHXParser() )[ 0 ];
            species_tree_0000.setRooted( true );
            rio = RIO.executeAnalysis( gene_trees_0000,
                                       species_tree_0000,
                                       ALGORITHM.GSDIR,
                                       REROOTING.BY_ALGORITHM,
                                       "",
                                       true,
                                       false,
                                       true );
            if ( rio.getAnalyzedGeneTrees().length != 3 ) {
                return false;
            }
            if ( rio.getExtNodesOfAnalyzedGeneTrees() != 2 ) {
                return false;
            }
            if ( rio.getGSDIRtaxCompBase() != TaxonomyComparisonBase.SCIENTIFIC_NAME ) {
                return false;
            }
            if ( rio.getRemovedGeneTreeNodes().size() != 0 ) {
                return false;
            }
            m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            if ( !m.getRowAsString( 0, ',' ).equals( "MOUSE1,3,0" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            if ( !m.getRowAsString( 1, ',' ).equals( "MOUSE2,0,3" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            //
            final String gene_trees_x_str = "(MOUSE1[&&NHX:S=MOUSE],MOUSE2[&&NHX:S=MOUSE])";
            final Phylogeny[] gene_trees_x = factory.create( gene_trees_x_str, nhx );
            final String species_trees_x_str = "[&&NHX:S=MOUSE];";
            final Phylogeny species_tree_x = factory.create( species_trees_x_str, new NHXParser() )[ 0 ];
            species_tree_x.setRooted( true );
            rio = RIO.executeAnalysis( gene_trees_x,
                                       species_tree_x,
                                       ALGORITHM.GSDIR,
                                       REROOTING.BY_ALGORITHM,
                                       "",
                                       true,
                                       false,
                                       true );
            if ( rio.getAnalyzedGeneTrees().length != 1 ) {
                return false;
            }
            if ( rio.getExtNodesOfAnalyzedGeneTrees() != 2 ) {
                return false;
            }
            if ( rio.getGSDIRtaxCompBase() != TaxonomyComparisonBase.SCIENTIFIC_NAME ) {
                return false;
            }
            if ( rio.getRemovedGeneTreeNodes().size() != 0 ) {
                return false;
            }
            m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            if ( !m.getRowAsString( 0, ',' ).equals( "MOUSE1,1,0" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            if ( !m.getRowAsString( 1, ',' ).equals( "MOUSE2,0,1" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            final String gene_trees_xx_str = "(MOUSE1[&&NHX:S=MOUSE],RAT1[&&NHX:S=RAT])";
            final Phylogeny[] gene_trees_xx = factory.create( gene_trees_xx_str, nhx );
            final String species_trees_xx_str = "([&&NHX:S=MOUSE],[&&NHX:S=RAT]);";
            final Phylogeny species_tree_xx = factory.create( species_trees_xx_str, new NHXParser() )[ 0 ];
            species_tree_xx.setRooted( true );
            rio = RIO.executeAnalysis( gene_trees_xx,
                                       species_tree_xx,
                                       ALGORITHM.GSDIR,
                                       REROOTING.BY_ALGORITHM,
                                       "",
                                       true,
                                       false,
                                       true );
            if ( rio.getAnalyzedGeneTrees().length != 1 ) {
                return false;
            }
            if ( rio.getExtNodesOfAnalyzedGeneTrees() != 2 ) {
                return false;
            }
            if ( rio.getGSDIRtaxCompBase() != TaxonomyComparisonBase.SCIENTIFIC_NAME ) {
                return false;
            }
            if ( rio.getRemovedGeneTreeNodes().size() != 0 ) {
                return false;
            }
            m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
            if ( !m.getRowAsString( 0, ',' ).equals( "MOUSE1,1,1" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            if ( !m.getRowAsString( 1, ',' ).equals( "RAT1,1,1" ) ) {
                System.out.println( m.toString() );
                return false;
            }
            final String gene_trees_1_str = "(((((MOUSE,RAT),HUMAN),CAEEL),YEAST),ARATH);"
                    + "((((MOUSE,RAT),HUMAN),(ARATH,YEAST)),CAEEL);" + "((MOUSE,RAT),(((ARATH,YEAST),CAEEL),HUMAN));"
                    + "(((((MOUSE,HUMAN),RAT),CAEEL),YEAST),ARATH);" + "((((HUMAN,MOUSE),RAT),(ARATH,YEAST)),CAEEL);";
            final Phylogeny[] gene_trees_1 = factory.create( gene_trees_1_str, nhx );
            final String species_trees_1_str = "(((((MOUSE,RAT),HUMAN),CAEEL),YEAST),ARATH);";
            final Phylogeny species_tree_1 = factory.create( species_trees_1_str, new NHXParser() )[ 0 ];
            species_tree_1.setRooted( true );
            PhylogenyMethods.transferNodeNameToField( species_tree_1, PhylogenyNodeField.TAXONOMY_CODE, true );
            rio = RIO.executeAnalysis( gene_trees_1,
                                       species_tree_1,
                                       ALGORITHM.GSDIR,
                                       REROOTING.BY_ALGORITHM,
                                       "",
                                       true,
                                       false,
                                       true );
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
            m = RIO.calculateOrthologTable( rio.getAnalyzedGeneTrees(), true );
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
            //
            RIO r0 = RIO.executeAnalysis( new File( PATH_TO_TEST_DATA + "rio_mb_taxcode.run1.t" ),
                                          new File( PATH_TO_TEST_DATA + "rio_tol_1.xml" ),
                                          ALGORITHM.GSDIR,
                                          REROOTING.BY_ALGORITHM,
                                          "",
                                          -1,
                                          -1,
                                          true,
                                          false,
                                          true );
            if ( r0.getGSDIRtaxCompBase() != TaxonomyComparisonBase.CODE ) {
                return false;
            }
            if ( r0.getAnalyzedGeneTrees().length != 201 ) {
                return false;
            }
            if ( r0.getExtNodesOfAnalyzedGeneTrees() != 6 ) {
                System.out.println( r0.getExtNodesOfAnalyzedGeneTrees() );
                return false;
            }
            if ( r0.getIntNodesOfAnalyzedGeneTrees() != 5 ) {
                return false;
            }
            if ( r0.getRemovedGeneTreeNodes().size() != 0 ) {
                return false;
            }
            if ( ForesterUtil.roundToInt( r0.getDuplicationsStatistics().median() ) != 1 ) {
                return false;
            }
            m = RIO.calculateOrthologTable( r0.getAnalyzedGeneTrees(), true );
            if ( !m.getRowAsString( 0, ',' ).equals( "A7SHU1_NEMVE,201,201,200,200,200,200" ) ) {
                System.out.println( m.getRowAsString( 0, ',' ) );
                return false;
            }
            if ( !m.getRowAsString( 1, ',' ).equals( "BCDO2_HUMAN,201,201,200,200,200,43" ) ) {
                System.out.println( m.getRowAsString( 1, ',' ) );
                return false;
            }
            if ( !m.getRowAsString( 2, ',' ).equals( "BCDO2_MOUSE,200,200,201,201,201,43" ) ) {
                System.out.println( m.getRowAsString( 2, ',' ) );
                return false;
            }
            if ( !m.getRowAsString( 3, ',' ).equals( "H2ZH97_CIOSA,200,200,201,201,201,201" ) ) {
                System.out.println( m.getRowAsString( 3, ',' ) );
                return false;
            }
            if ( !m.getRowAsString( 4, ',' ).equals( "Q1RLW1_DANRE,200,200,201,201,201,43" ) ) {
                System.out.println( m.getRowAsString( 4, ',' ) );
                return false;
            }
            if ( !m.getRowAsString( 5, ',' ).equals( "Q6DIN7_XENTR,200,43,43,201,43,201" ) ) {
                System.out.println( m.getRowAsString( 5, ',' ) );
                return false;
            }
            r0 = RIO.executeAnalysis( new File( PATH_TO_TEST_DATA + "rio_mb_taxid.run1.t" ),
                                      new File( PATH_TO_TEST_DATA + "rio_tol_1.xml" ),
                                      ALGORITHM.GSDIR,
                                      REROOTING.BY_ALGORITHM,
                                      "",
                                      -1,
                                      -1,
                                      true,
                                      false,
                                      true );
            if ( r0.getGSDIRtaxCompBase() != TaxonomyComparisonBase.ID ) {
                return false;
            }
            if ( r0.getAnalyzedGeneTrees().length != 201 ) {
                return false;
            }
            if ( r0.getExtNodesOfAnalyzedGeneTrees() != 6 ) {
                return false;
            }
            if ( r0.getIntNodesOfAnalyzedGeneTrees() != 5 ) {
                return false;
            }
            if ( r0.getRemovedGeneTreeNodes().size() != 0 ) {
                return false;
            }
            if ( ForesterUtil.roundToInt( r0.getDuplicationsStatistics().median() ) != 1 ) {
                return false;
            }
            m = RIO.calculateOrthologTable( r0.getAnalyzedGeneTrees(), true );
            if ( !m.getRowAsString( 0, ',' ).equals( "A7SHU1_45351,201,200,201,200,200,200" ) ) {
                System.out.println( m.getRowAsString( 0, ',' ) );
                return false;
            }
            //mouse
            if ( !m.getRowAsString( 1, ',' ).equals( "BCDO2_10090,200,201,200,201,201,43" ) ) {
                System.out.println( m.getRowAsString( 1, ',' ) );
                return false;
            }
            //human
            if ( !m.getRowAsString( 2, ',' ).equals( "BCDO2_9606,201,200,201,200,200,43" ) ) {
                System.out.println( m.getRowAsString( 2, ',' ) );
                return false;
            }
            if ( !m.getRowAsString( 3, ',' ).equals( "H2ZH97_51511,200,201,200,201,201,201" ) ) {
                System.out.println( m.getRowAsString( 3, ',' ) );
                return false;
            }
            if ( !m.getRowAsString( 4, ',' ).equals( "Q1RLW1_7955,200,201,200,201,201,43" ) ) {
                System.out.println( m.getRowAsString( 4, ',' ) );
                return false;
            }
            if ( !m.getRowAsString( 5, ',' ).equals( "Q6DIN7_8364,200,43,43,201,43,201" ) ) {
                System.out.println( m.getRowAsString( 5, ',' ) );
                return false;
            }
            r0 = RIO.executeAnalysis( new File( PATH_TO_TEST_DATA + "rio_mb_taxsn.run1.t" ),
                                      new File( PATH_TO_TEST_DATA + "rio_tol_1.xml" ),
                                      ALGORITHM.GSDIR,
                                      REROOTING.BY_ALGORITHM,
                                      "",
                                      -1,
                                      -1,
                                      true,
                                      false,
                                      true );
            if ( r0.getGSDIRtaxCompBase() != TaxonomyComparisonBase.SCIENTIFIC_NAME ) {
                return false;
            }
            if ( r0.getAnalyzedGeneTrees().length != 201 ) {
                return false;
            }
            if ( r0.getExtNodesOfAnalyzedGeneTrees() != 6 ) {
                System.out.println( r0.getExtNodesOfAnalyzedGeneTrees() );
                return false;
            }
            //            if ( r0.getIntNodesOfAnalyzedGeneTrees() != 5 ) {
            //                return false;
            //            }
            //            if ( r0.getRemovedGeneTreeNodes().size() != 0 ) {
            //                return false;
            //            }
            //            if ( ForesterUtil.roundToInt( r0.getDuplicationsStatistics().median() ) != 1 ) {
            //                return false;
            //            }
            //            m = RIO.calculateOrthologTable( r0.getAnalyzedGeneTrees(), true );
            //            if ( !m.getRowAsString( 0, ',' ).equals( "A7SHU1_Nematostella_vectensis,201,201,200,200,200,200" ) ) {
            //                System.out.println( m.getRowAsString( 0, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 1, ',' ).equals( "BCDO2_Homo_sapiens,201,201,200,200,200,43" ) ) {
            //                System.out.println( m.getRowAsString( 1, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 2, ',' ).equals( "BCDO2_Mus_musculus,200,200,201,201,201,43" ) ) {
            //                System.out.println( m.getRowAsString( 2, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 3, ',' ).equals( "H2ZH97_Ciona_savignyi,200,200,201,201,201,201" ) ) {
            //                System.out.println( m.getRowAsString( 3, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 4, ',' ).equals( "Q1RLW1_Danio_rerio,200,200,201,201,201,43" ) ) {
            //                System.out.println( m.getRowAsString( 4, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 5, ',' ).equals( "Q6DIN7_Xenopus_tropicalis,200,43,43,201,43,201" ) ) {
            //                System.out.println( m.getRowAsString( 5, ',' ) );
            //                return false;
            //            }
            //
            //            r0 = RIO.executeAnalysis( new File( PATH_TO_TEST_DATA + "rio_mb_taxsn.run1.t" ),
            //                                      new File( PATH_TO_TEST_DATA + "rio_tol_1.xml" ),
            //                                      ALGORITHM.GSDIR,
            //                                      REROOTING.MIDPOINT,
            //                                      "",
            //                                      -1,
            //                                      -1,
            //                                      true,
            //                                      false,
            //                                      true );
            //            if ( r0.getGSDIRtaxCompBase() != TaxonomyComparisonBase.SCIENTIFIC_NAME ) {
            //                return false;
            //            }
            //            if ( r0.getAnalyzedGeneTrees().length != 201 ) {
            //                return false;
            //            }
            //            if ( r0.getExtNodesOfAnalyzedGeneTrees() != 6 ) {
            //                return false;
            //            }
            //            if ( r0.getIntNodesOfAnalyzedGeneTrees() != 5 ) {
            //                return false;
            //            }
            //            if ( r0.getRemovedGeneTreeNodes().size() != 0 ) {
            //                return false;
            //            }
            //            if ( ForesterUtil.roundToInt( r0.getDuplicationsStatistics().median() ) != 2 ) {
            //                return false;
            //            }
            //            m = RIO.calculateOrthologTable( r0.getAnalyzedGeneTrees(), true );
            //            if ( !m.getRowAsString( 0, ',' ).equals( "A7SHU1_Nematostella_vectensis,201,94,93,160,93,93" ) ) {
            //                System.out.println( m.getRowAsString( 0, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 1, ',' ).equals( "BCDO2_Homo_sapiens,94,201,200,53,200,43" ) ) {
            //                System.out.println( m.getRowAsString( 1, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 2, ',' ).equals( "BCDO2_Mus_musculus,93,200,201,53,201,43" ) ) {
            //                System.out.println( m.getRowAsString( 2, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 3, ',' ).equals( "H2ZH97_Ciona_savignyi,160,53,53,201,53,53" ) ) {
            //                System.out.println( m.getRowAsString( 3, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 4, ',' ).equals( "Q1RLW1_Danio_rerio,93,200,201,53,201,43" ) ) {
            //                System.out.println( m.getRowAsString( 4, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 5, ',' ).equals( "Q6DIN7_Xenopus_tropicalis,93,43,43,53,43,201" ) ) {
            //                System.out.println( m.getRowAsString( 5, ',' ) );
            //                return false;
            //            }
            //
            //            r0 = RIO.executeAnalysis( new File( PATH_TO_TEST_DATA + "rio_mb_taxsn.run1.t" ),
            //                                      new File( PATH_TO_TEST_DATA + "rio_tol_1.xml" ),
            //                                      ALGORITHM.GSDIR,
            //                                      REROOTING.OUTGROUP,
            //                                      "H2ZH97_Ciona_savignyi",
            //                                      -1,
            //                                      -1,
            //                                      true,
            //                                      false,
            //                                      true );
            //            if ( r0.getGSDIRtaxCompBase() != TaxonomyComparisonBase.SCIENTIFIC_NAME ) {
            //                return false;
            //            }
            //            if ( r0.getAnalyzedGeneTrees().length != 201 ) {
            //                return false;
            //            }
            //            if ( r0.getExtNodesOfAnalyzedGeneTrees() != 6 ) {
            //                return false;
            //            }
            //            if ( r0.getIntNodesOfAnalyzedGeneTrees() != 5 ) {
            //                return false;
            //            }
            //            if ( r0.getRemovedGeneTreeNodes().size() != 0 ) {
            //                return false;
            //            }
            //            if ( ForesterUtil.roundToInt( r0.getDuplicationsStatistics().median() ) != 2 ) {
            //                return false;
            //            }
            //            m = RIO.calculateOrthologTable( r0.getAnalyzedGeneTrees(), true );
            //            if ( !m.getRowAsString( 0, ',' ).equals( "A7SHU1_Nematostella_vectensis,201,201,200,0,200,200" ) ) {
            //                System.out.println( m.getRowAsString( 0, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 1, ',' ).equals( "BCDO2_Homo_sapiens,201,201,200,0,200,43" ) ) {
            //                System.out.println( m.getRowAsString( 1, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 2, ',' ).equals( "BCDO2_Mus_musculus,200,200,201,0,201,43" ) ) {
            //                System.out.println( m.getRowAsString( 2, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 3, ',' ).equals( "H2ZH97_Ciona_savignyi,0,0,0,201,0,0" ) ) {
            //                System.out.println( m.getRowAsString( 3, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 4, ',' ).equals( "Q1RLW1_Danio_rerio,200,200,201,0,201,43" ) ) {
            //                System.out.println( m.getRowAsString( 4, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 5, ',' ).equals( "Q6DIN7_Xenopus_tropicalis,200,43,43,0,43,201" ) ) {
            //                System.out.println( m.getRowAsString( 5, ',' ) );
            //                return false;
            //            }
            //
            //
            //            r0 = RIO.executeAnalysis( new File( PATH_TO_TEST_DATA + "rio_mb_taxsn.run1.t" ),
            //                                      new File( PATH_TO_TEST_DATA + "rio_tol_1.xml" ),
            //                                      ALGORITHM.GSDIR,
            //                                      REROOTING.NONE,
            //                                      null,
            //                                      10,
            //                                      19,
            //                                      true,
            //                                      false,
            //                                      true );
            //            if ( r0.getGSDIRtaxCompBase() != TaxonomyComparisonBase.SCIENTIFIC_NAME ) {
            //                return false;
            //            }
            //            if ( r0.getAnalyzedGeneTrees().length != 10 ) {
            //                return false;
            //            }
            //            if ( r0.getExtNodesOfAnalyzedGeneTrees() != 6 ) {
            //                return false;
            //            }
            //            if ( r0.getIntNodesOfAnalyzedGeneTrees() != 5 ) {
            //                return false;
            //            }
            //            if ( r0.getRemovedGeneTreeNodes().size() != 0 ) {
            //                return false;
            //            }
            //            if ( ForesterUtil.roundToInt( r0.getDuplicationsStatistics().median() ) != 4 ) {
            //                return false;
            //            }
            //            m = RIO.calculateOrthologTable( r0.getAnalyzedGeneTrees(), true );
            //            if ( !m.getRowAsString( 0, ',' ).equals( "A7SHU1_Nematostella_vectensis,10,0,0,10,0,0" ) ) {
            //                System.out.println( m.getRowAsString( 0, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 1, ',' ).equals( "BCDO2_Homo_sapiens,0,10,0,0,0,0" ) ) {
            //                System.out.println( m.getRowAsString( 1, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 2, ',' ).equals( "BCDO2_Mus_musculus,0,0,10,0,0,0" ) ) {
            //                System.out.println( m.getRowAsString( 2, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 3, ',' ).equals( "H2ZH97_Ciona_savignyi,10,0,0,10,0,0" ) ) {
            //                System.out.println( m.getRowAsString( 3, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 4, ',' ).equals( "Q1RLW1_Danio_rerio,0,0,0,0,10,0" ) ) {
            //                System.out.println( m.getRowAsString( 4, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 5, ',' ).equals( "Q6DIN7_Xenopus_tropicalis,0,0,0,0,0,10" ) ) {
            //                System.out.println( m.getRowAsString( 5, ',' ) );
            //                return false;
            //            }
            //
            //            r0 = RIO.executeAnalysis( new File( PATH_TO_TEST_DATA + "rio_mb_taxcode_1.run1.t" ),
            //                                      new File( PATH_TO_TEST_DATA + "rio_tol_1.xml" ),
            //                                      ALGORITHM.GSDIR,
            //                                      REROOTING.BY_ALGORITHM,
            //                                      "",
            //                                      -1,
            //                                      -1,
            //                                      true,
            //                                      false,
            //                                      true );
            //            if ( r0.getGSDIRtaxCompBase() != TaxonomyComparisonBase.CODE ) {
            //                return false;
            //            }
            //            if ( r0.getAnalyzedGeneTrees().length != 201 ) {
            //                return false;
            //            }
            //            if ( r0.getExtNodesOfAnalyzedGeneTrees() != 3 ) {
            //                return false;
            //            }
            //            if ( r0.getIntNodesOfAnalyzedGeneTrees() != 2 ) {
            //                return false;
            //            }
            //            if ( r0.getRemovedGeneTreeNodes().size() != 3 ) {
            //                return false;
            //            }
            //            if ( ForesterUtil.roundToInt( r0.getDuplicationsStatistics().median() ) != 0 ) {
            //                return false;
            //            }
            //            m = RIO.calculateOrthologTable( r0.getAnalyzedGeneTrees(), true );
            //            if ( !m.getRowAsString( 0, ',' ).equals( "BCDO2_HUMAN,201,201,201" ) ) {
            //                System.out.println( m.getRowAsString( 0, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 1, ',' ).equals( "Q1RLW1_DANRE,201,201,201" ) ) {
            //                System.out.println( m.getRowAsString( 1, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 2, ',' ).equals( "Q6DIN7_XENTR,201,201,201" ) ) {
            //                System.out.println( m.getRowAsString( 2, ',' ) );
            //                return false;
            //            }
            //
            //
            //            r0 = RIO.executeAnalysis( new File( PATH_TO_TEST_DATA + "rio_mb_taxcode_2.run1.t" ),
            //                                      new File( PATH_TO_TEST_DATA + "rio_tol_1.xml" ),
            //                                      ALGORITHM.GSDIR,
            //                                      REROOTING.BY_ALGORITHM,
            //                                      "",
            //                                      -1,
            //                                      -1,
            //                                      true,
            //                                      false,
            //                                      true );
            //            if ( r0.getGSDIRtaxCompBase() != TaxonomyComparisonBase.CODE ) {
            //                return false;
            //            }
            //            if ( r0.getAnalyzedGeneTrees().length != 201 ) {
            //                return false;
            //            }
            //            if ( r0.getExtNodesOfAnalyzedGeneTrees() != 6 ) {
            //                return false;
            //            }
            //            if ( r0.getIntNodesOfAnalyzedGeneTrees() != 5 ) {
            //                return false;
            //            }
            //            if ( r0.getRemovedGeneTreeNodes().size() != 0 ) {
            //                return false;
            //            }
            //            if ( ForesterUtil.roundToInt( r0.getDuplicationsStatistics().median() ) != 1 ) {
            //                return false;
            //            }
            //            m = RIO.calculateOrthologTable( r0.getAnalyzedGeneTrees(), true );
            //            if ( !m.getRowAsString( 0, ',' ).equals( "A7SHU1_NEMVE&1,201,201,200,200,200,200" ) ) {
            //                System.out.println( m.getRowAsString( 0, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 1, ',' ).equals( "BCDO2_HUMAN+,201,201,200,200,200,43" ) ) {
            //                System.out.println( m.getRowAsString( 1, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 2, ',' ).equals( "BCDO2_MOUSE,200,200,201,201,201,43" ) ) {
            //                System.out.println( m.getRowAsString( 2, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 3, ',' ).equals( "CIOSA,200,200,201,201,201,201" ) ) {
            //                System.out.println( m.getRowAsString( 3, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 4, ',' ).equals( "Q1RLW1_DANRE/12-45,200,200,201,201,201,43" ) ) {
            //                System.out.println( m.getRowAsString( 4, ',' ) );
            //                return false;
            //            }
            //            if ( !m.getRowAsString( 5, ',' ).equals( "Q6DIN7_XENTR-LOUSE,200,43,43,201,43,201" ) ) {
            //                System.out.println( m.getRowAsString( 5, ',' ) );
            //                return false;
            //            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testRIO_GSDIR_Iterating() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final NHXParser nhx = new NHXParser();
            nhx.setReplaceUnderscores( false );
            nhx.setIgnoreQuotes( true );
            nhx.setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            final String gene_trees_1_str = "(((((MOUSE,RAT),HUMAN),CAEEL),YEAST),ARATH);"
                    + "((((MOUSE,RAT),HUMAN),(ARATH,YEAST)),CAEEL);" + "((MOUSE,RAT),(((ARATH,YEAST),CAEEL),HUMAN));"
                    + "(((((MOUSE,HUMAN),RAT),CAEEL),YEAST),ARATH);" + "((((HUMAN,MOUSE),RAT),(ARATH,YEAST)),CAEEL);";
            nhx.setSource( gene_trees_1_str );
            final String species_trees_1_str = "(((((MOUSE,RAT),HUMAN),CAEEL),YEAST),ARATH);";
            final Phylogeny species_tree_1 = factory.create( species_trees_1_str, new NHXParser() )[ 0 ];
            species_tree_1.setRooted( true );
            PhylogenyMethods.transferNodeNameToField( species_tree_1, PhylogenyNodeField.TAXONOMY_CODE, true );
            //Archaeopteryx.createApplication( species_trees_1 );
            RIO rio = RIO.executeAnalysis( nhx,
                                           species_tree_1,
                                           ALGORITHM.GSDIR,
                                           REROOTING.BY_ALGORITHM,
                                           "",
                                           true,
                                           false,
                                           true );
            if ( rio.getExtNodesOfAnalyzedGeneTrees() != 6 ) {
                return false;
            }
            if ( rio.getGSDIRtaxCompBase() != TaxonomyComparisonBase.CODE ) {
                return false;
            }
            if ( rio.getRemovedGeneTreeNodes().size() != 0 ) {
                return false;
            }
            IntMatrix m = rio.getOrthologTable();
            //System.out.println( m.toString() );
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
            final String species_trees_2_str = "((((MOUSE,RAT,HUMAN),CAEEL),YEAST),ARATH);";
            final Phylogeny species_tree_2 = factory.create( species_trees_2_str, new NHXParser() )[ 0 ];
            species_tree_2.setRooted( true );
            PhylogenyMethods.transferNodeNameToField( species_tree_2, PhylogenyNodeField.TAXONOMY_CODE, true );
            rio = RIO.executeAnalysis( nhx,
                                       species_tree_2,
                                       ALGORITHM.GSDIR,
                                       REROOTING.BY_ALGORITHM,
                                       "",
                                       true,
                                       false,
                                       true );
            m = rio.getOrthologTable();
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
}