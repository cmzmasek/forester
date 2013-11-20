// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.test;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.SortedSet;

import org.forester.application.support_transfer;
import org.forester.archaeopteryx.TreePanelUtil;
import org.forester.development.DevelopmentTools;
import org.forester.evoinference.TestPhylogenyReconstruction;
import org.forester.evoinference.matrix.character.CharacterStateMatrix;
import org.forester.evoinference.matrix.character.CharacterStateMatrix.BinaryStates;
import org.forester.go.TestGo;
import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser.INDIVIDUAL_SCORE_CUTOFF;
import org.forester.io.parsers.nexus.NexusBinaryStatesMatrixParser;
import org.forester.io.parsers.nexus.NexusCharactersParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.tol.TolParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.io.writers.SequenceWriter;
import org.forester.msa.BasicMsa;
import org.forester.msa.Mafft;
import org.forester.msa.Msa;
import org.forester.msa.MsaInferrer;
import org.forester.msa.MsaMethods;
import org.forester.pccx.TestPccx;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyBranch;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Accession.Source;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Polygon;
import org.forester.phylogeny.data.PropertiesMap;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.protein.BasicDomain;
import org.forester.protein.BasicProtein;
import org.forester.protein.Domain;
import org.forester.protein.Protein;
import org.forester.protein.ProteinId;
import org.forester.rio.TestRIO;
import org.forester.sdi.SDI;
import org.forester.sdi.SDIR;
import org.forester.sdi.TestGSDI;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.Sequence;
import org.forester.species.BasicSpecies;
import org.forester.species.Species;
import org.forester.surfacing.TestSurfacing;
import org.forester.tools.ConfidenceAssessor;
import org.forester.tools.SupportCount;
import org.forester.tools.TreeSplitMatrix;
import org.forester.util.AsciiHistogram;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.GeneralTable;
import org.forester.util.SequenceAccessionTools;
import org.forester.ws.seqdb.SequenceDatabaseEntry;
import org.forester.ws.seqdb.SequenceDbWsTools;
import org.forester.ws.seqdb.UniProtTaxonomy;
import org.forester.ws.wabi.TxSearch;
import org.forester.ws.wabi.TxSearch.RANKS;
import org.forester.ws.wabi.TxSearch.TAX_NAME_CLASS;
import org.forester.ws.wabi.TxSearch.TAX_RANK;

@SuppressWarnings( "unused")
public final class Test {

    private final static boolean PERFORM_DB_TESTS          = true;
    private final static double  ZERO_DIFF                 = 1.0E-9;
    private final static String  PATH_TO_TEST_DATA         = System.getProperty( "user.dir" )
                                                                   + ForesterUtil.getFileSeparator() + "test_data"
                                                                   + ForesterUtil.getFileSeparator();
    private final static String  PATH_TO_RESOURCES         = System.getProperty( "user.dir" )
                                                                   + ForesterUtil.getFileSeparator() + "resources"
                                                                   + ForesterUtil.getFileSeparator();
    private final static boolean USE_LOCAL_PHYLOXML_SCHEMA = true;
    private static final String  PHYLOXML_REMOTE_XSD       = ForesterConstants.PHYLO_XML_LOCATION + "/"
                                                                   + ForesterConstants.PHYLO_XML_VERSION + "/"
                                                                   + ForesterConstants.PHYLO_XML_XSD;
    private static final String  PHYLOXML_LOCAL_XSD        = PATH_TO_RESOURCES + "phyloxml_schema/"
                                                                   + ForesterConstants.PHYLO_XML_VERSION + "/"
                                                                   + ForesterConstants.PHYLO_XML_XSD;

    public static boolean testOverlapRemoval() {
        try {
            final Domain d0 = new BasicDomain( "d0", ( short ) 2, ( short ) 5, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d1 = new BasicDomain( "d1", ( short ) 7, ( short ) 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d2 = new BasicDomain( "d2", ( short ) 0, ( short ) 20, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d3 = new BasicDomain( "d3", ( short ) 9, ( short ) 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d4 = new BasicDomain( "d4", ( short ) 7, ( short ) 8, ( short ) 1, ( short ) 1, 0.1, 1 );
            final List<Boolean> covered = new ArrayList<Boolean>();
            covered.add( true ); // 0
            covered.add( false ); // 1
            covered.add( true ); // 2
            covered.add( false ); // 3
            covered.add( true ); // 4
            covered.add( true ); // 5
            covered.add( false ); // 6
            covered.add( true ); // 7
            covered.add( true ); // 8
            if ( ForesterUtil.calculateOverlap( d0, covered ) != 3 ) {
                return false;
            }
            if ( ForesterUtil.calculateOverlap( d1, covered ) != 2 ) {
                return false;
            }
            if ( ForesterUtil.calculateOverlap( d2, covered ) != 6 ) {
                return false;
            }
            if ( ForesterUtil.calculateOverlap( d3, covered ) != 0 ) {
                return false;
            }
            if ( ForesterUtil.calculateOverlap( d4, covered ) != 2 ) {
                return false;
            }
            final Domain a = new BasicDomain( "a", ( short ) 2, ( short ) 5, ( short ) 1, ( short ) 1, 0.01, 1 );
            final Domain b = new BasicDomain( "b", ( short ) 2, ( short ) 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Protein ab = new BasicProtein( "ab", "varanus", 0 );
            ab.addProteinDomain( a );
            ab.addProteinDomain( b );
            final Protein ab_s0 = ForesterUtil.removeOverlappingDomains( 3, false, ab );
            if ( ab.getNumberOfProteinDomains() != 2 ) {
                return false;
            }
            if ( ab_s0.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( !ab_s0.getProteinDomain( 0 ).getDomainId().equals( "a" ) ) {
                return false;
            }
            final Protein ab_s1 = ForesterUtil.removeOverlappingDomains( 4, false, ab );
            if ( ab.getNumberOfProteinDomains() != 2 ) {
                return false;
            }
            if ( ab_s1.getNumberOfProteinDomains() != 2 ) {
                return false;
            }
            final Domain c = new BasicDomain( "c", ( short ) 20000, ( short ) 20500, ( short ) 1, ( short ) 1, 10, 1 );
            final Domain d = new BasicDomain( "d",
                                              ( short ) 10000,
                                              ( short ) 10500,
                                              ( short ) 1,
                                              ( short ) 1,
                                              0.0000001,
                                              1 );
            final Domain e = new BasicDomain( "e", ( short ) 5000, ( short ) 5500, ( short ) 1, ( short ) 1, 0.0001, 1 );
            final Protein cde = new BasicProtein( "cde", "varanus", 0 );
            cde.addProteinDomain( c );
            cde.addProteinDomain( d );
            cde.addProteinDomain( e );
            final Protein cde_s0 = ForesterUtil.removeOverlappingDomains( 0, false, cde );
            if ( cde.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( cde_s0.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            final Domain f = new BasicDomain( "f", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 10, 1 );
            final Domain g = new BasicDomain( "g", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.01, 1 );
            final Domain h = new BasicDomain( "h", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.0001, 1 );
            final Domain i = new BasicDomain( "i", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.5, 1 );
            final Domain i2 = new BasicDomain( "i", ( short ) 5, ( short ) 30, ( short ) 1, ( short ) 1, 0.5, 10 );
            final Protein fghi = new BasicProtein( "fghi", "varanus", 0 );
            fghi.addProteinDomain( f );
            fghi.addProteinDomain( g );
            fghi.addProteinDomain( h );
            fghi.addProteinDomain( i );
            fghi.addProteinDomain( i );
            fghi.addProteinDomain( i );
            fghi.addProteinDomain( i2 );
            final Protein fghi_s0 = ForesterUtil.removeOverlappingDomains( 10, false, fghi );
            if ( fghi.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            if ( fghi_s0.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( !fghi_s0.getProteinDomain( 0 ).getDomainId().equals( "h" ) ) {
                return false;
            }
            final Protein fghi_s1 = ForesterUtil.removeOverlappingDomains( 11, false, fghi );
            if ( fghi.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            if ( fghi_s1.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            final Domain j = new BasicDomain( "j", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 10, 1 );
            final Domain k = new BasicDomain( "k", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.01, 1 );
            final Domain l = new BasicDomain( "l", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.0001, 1 );
            final Domain m = new BasicDomain( "m", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 4, 0.5, 1 );
            final Domain m0 = new BasicDomain( "m", ( short ) 10, ( short ) 20, ( short ) 2, ( short ) 4, 0.5, 1 );
            final Domain m1 = new BasicDomain( "m", ( short ) 10, ( short ) 20, ( short ) 3, ( short ) 4, 0.5, 1 );
            final Domain m2 = new BasicDomain( "m", ( short ) 5, ( short ) 30, ( short ) 4, ( short ) 4, 0.5, 10 );
            final Protein jklm = new BasicProtein( "jklm", "varanus", 0 );
            jklm.addProteinDomain( j );
            jklm.addProteinDomain( k );
            jklm.addProteinDomain( l );
            jklm.addProteinDomain( m );
            jklm.addProteinDomain( m0 );
            jklm.addProteinDomain( m1 );
            jklm.addProteinDomain( m2 );
            final Protein jklm_s0 = ForesterUtil.removeOverlappingDomains( 10, false, jklm );
            if ( jklm.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            if ( jklm_s0.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( !jklm_s0.getProteinDomain( 0 ).getDomainId().equals( "l" ) ) {
                return false;
            }
            final Protein jklm_s1 = ForesterUtil.removeOverlappingDomains( 11, false, jklm );
            if ( jklm.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            if ( jklm_s1.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            final Domain only = new BasicDomain( "only", ( short ) 5, ( short ) 30, ( short ) 4, ( short ) 4, 0.5, 10 );
            final Protein od = new BasicProtein( "od", "varanus", 0 );
            od.addProteinDomain( only );
            final Protein od_s0 = ForesterUtil.removeOverlappingDomains( 0, false, od );
            if ( od.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( od_s0.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    public static boolean testEngulfingOverlapRemoval() {
        try {
            final Domain d0 = new BasicDomain( "d0", 0, 8, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d1 = new BasicDomain( "d1", 0, 1, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d2 = new BasicDomain( "d2", 0, 2, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d3 = new BasicDomain( "d3", 7, 8, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d4 = new BasicDomain( "d4", 7, 9, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d5 = new BasicDomain( "d4", 0, 9, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d6 = new BasicDomain( "d4", 4, 5, ( short ) 1, ( short ) 1, 0.1, 1 );
            final List<Boolean> covered = new ArrayList<Boolean>();
            covered.add( true ); // 0
            covered.add( false ); // 1
            covered.add( true ); // 2
            covered.add( false ); // 3
            covered.add( true ); // 4
            covered.add( true ); // 5
            covered.add( false ); // 6
            covered.add( true ); // 7
            covered.add( true ); // 8
            if ( ForesterUtil.isEngulfed( d0, covered ) ) {
                return false;
            }
            if ( ForesterUtil.isEngulfed( d1, covered ) ) {
                return false;
            }
            if ( ForesterUtil.isEngulfed( d2, covered ) ) {
                return false;
            }
            if ( !ForesterUtil.isEngulfed( d3, covered ) ) {
                return false;
            }
            if ( ForesterUtil.isEngulfed( d4, covered ) ) {
                return false;
            }
            if ( ForesterUtil.isEngulfed( d5, covered ) ) {
                return false;
            }
            if ( !ForesterUtil.isEngulfed( d6, covered ) ) {
                return false;
            }
            final Domain a = new BasicDomain( "a", 0, 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain b = new BasicDomain( "b", 8, 20, ( short ) 1, ( short ) 1, 0.2, 1 );
            final Domain c = new BasicDomain( "c", 15, 16, ( short ) 1, ( short ) 1, 0.3, 1 );
            final Protein abc = new BasicProtein( "abc", "nemve", 0 );
            abc.addProteinDomain( a );
            abc.addProteinDomain( b );
            abc.addProteinDomain( c );
            final Protein abc_r1 = ForesterUtil.removeOverlappingDomains( 3, false, abc );
            final Protein abc_r2 = ForesterUtil.removeOverlappingDomains( 3, true, abc );
            if ( abc.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( abc_r1.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( abc_r2.getNumberOfProteinDomains() != 2 ) {
                return false;
            }
            if ( !abc_r2.getProteinDomain( 0 ).getDomainId().equals( "a" ) ) {
                return false;
            }
            if ( !abc_r2.getProteinDomain( 1 ).getDomainId().equals( "b" ) ) {
                return false;
            }
            final Domain d = new BasicDomain( "d", 0, 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain e = new BasicDomain( "e", 8, 20, ( short ) 1, ( short ) 1, 0.3, 1 );
            final Domain f = new BasicDomain( "f", 15, 16, ( short ) 1, ( short ) 1, 0.2, 1 );
            final Protein def = new BasicProtein( "def", "nemve", 0 );
            def.addProteinDomain( d );
            def.addProteinDomain( e );
            def.addProteinDomain( f );
            final Protein def_r1 = ForesterUtil.removeOverlappingDomains( 5, false, def );
            final Protein def_r2 = ForesterUtil.removeOverlappingDomains( 5, true, def );
            if ( def.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( def_r1.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( def_r2.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( !def_r2.getProteinDomain( 0 ).getDomainId().equals( "d" ) ) {
                return false;
            }
            if ( !def_r2.getProteinDomain( 1 ).getDomainId().equals( "f" ) ) {
                return false;
            }
            if ( !def_r2.getProteinDomain( 2 ).getDomainId().equals( "e" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    public static boolean isEqual( final double a, final double b ) {
        return ( ( Math.abs( a - b ) ) < Test.ZERO_DIFF );
    }

    public static void main( final String[] args ) {
        System.out.println( "[Java version: " + ForesterUtil.JAVA_VERSION + " " + ForesterUtil.JAVA_VENDOR + "]" );
        System.out.println( "[OS: " + ForesterUtil.OS_NAME + " " + ForesterUtil.OS_ARCH + " " + ForesterUtil.OS_VERSION
                + "]" );
        Locale.setDefault( Locale.US );
        System.out.println( "[Locale: " + Locale.getDefault() + "]" );
        int failed = 0;
        int succeeded = 0;
        System.out.print( "[Test if directory with files for testing exists/is readable: " );
        if ( Test.testDir( PATH_TO_TEST_DATA ) ) {
            System.out.println( "OK.]" );
        }
        else {
            System.out.println( "could not find/read from directory \"" + PATH_TO_TEST_DATA + "\".]" );
            System.out.println( "Testing aborted." );
            System.exit( -1 );
        }
        System.out.print( "[Test if resources directory exists/is readable: " );
        if ( testDir( PATH_TO_RESOURCES ) ) {
            System.out.println( "OK.]" );
        }
        else {
            System.out.println( "could not find/read from directory \"" + Test.PATH_TO_RESOURCES + "\".]" );
            System.out.println( "Testing aborted." );
            System.exit( -1 );
        }
        final long start_time = new Date().getTime();
        System.out.print( "Basic node methods: " );
        if ( Test.testBasicNodeMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Protein id: " );
        if ( !testProteinId() ) {
            System.out.println( "failed." );
            failed++;
        }
        else {
            succeeded++;
        }
        System.out.println( "OK." );
        System.out.print( "Species: " );
        if ( !testSpecies() ) {
            System.out.println( "failed." );
            failed++;
        }
        else {
            succeeded++;
        }
        System.out.println( "OK." );
        System.out.print( "Basic domain: " );
        if ( !testBasicDomain() ) {
            System.out.println( "failed." );
            failed++;
        }
        else {
            succeeded++;
        }
        System.out.println( "OK." );
        System.out.print( "Basic protein: " );
        if ( !testBasicProtein() ) {
            System.out.println( "failed." );
            failed++;
        }
        else {
            succeeded++;
        }
        System.out.println( "OK." );
        System.out.print( "Sequence writer: " );
        if ( testSequenceWriter() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Sequence id parsing: " );
        if ( testSequenceIdParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "UniProtKB id extraction: " );
        if ( Test.testExtractUniProtKbProteinSeqIdentifier() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Sequence DB tools 1: " );
        if ( testSequenceDbWsTools1() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        if ( PERFORM_DB_TESTS ) {
            System.out.print( "Ebi Entry Retrieval: " );
            if ( Test.testEbiEntryRetrieval() ) {
                System.out.println( "OK." );
                succeeded++;
            }
            else {
                System.out.println( "failed." );
                failed++;
            }
        }
        // System.exit( 0 );
        if ( PERFORM_DB_TESTS ) {
            System.out.print( "Sequence DB tools 2: " );
            if ( testSequenceDbWsTools2() ) {
                System.out.println( "OK." );
                succeeded++;
            }
            else {
                System.out.println( "failed." );
                failed++;
                System.exit( -1 );
            }
        }
        // System.exit( 0 );
        System.out.print( "Hmmscan output parser: " );
        if ( testHmmscanOutputParser() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        //
        System.out.print( "Overlap removal: " );
        if ( !org.forester.test.Test.testOverlapRemoval() ) {
            System.out.println( "failed." );
            failed++;
        }
        else {
            succeeded++;
        }
        System.out.println( "OK." );
        System.out.print( "Engulfing overlap removal: " );
        if ( !Test.testEngulfingOverlapRemoval() ) {
            System.out.println( "failed." );
            failed++;
        }
        else {
            succeeded++;
        }
        System.out.println( "OK." );
        //
        System.out.print( "Taxonomy code extraction: " );
        if ( Test.testExtractTaxonomyCodeFromNodeName() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "SN extraction: " );
        if ( Test.testExtractSNFromNodeName() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Taxonomy extraction (general): " );
        if ( Test.testTaxonomyExtraction() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Uri for Aptx web sequence accession: " );
        if ( Test.testCreateUriForSeqWeb() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic node construction and parsing of NHX (node level): " );
        if ( Test.testNHXNodeParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "NHX parsing iterating: " );
        if ( Test.testNHParsingIter() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "NH parsing: " );
        if ( Test.testNHParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Conversion to NHX (node level): " );
        if ( Test.testNHXconversion() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "NHX parsing: " );
        if ( Test.testNHXParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "NHX parsing with quotes: " );
        if ( Test.testNHXParsingQuotes() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "NHX parsing (MrBayes): " );
        if ( Test.testNHXParsingMB() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Nexus characters parsing: " );
        if ( Test.testNexusCharactersParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Nexus tree parsing iterating: " );
        if ( Test.testNexusTreeParsingIterating() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Nexus tree parsing: " );
        if ( Test.testNexusTreeParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Nexus tree parsing (translating): " );
        if ( Test.testNexusTreeParsingTranslating() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Nexus matrix parsing: " );
        if ( Test.testNexusMatrixParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic phyloXML parsing: " );
        if ( Test.testBasicPhyloXMLparsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic phyloXML parsing (validating against schema): " );
        if ( testBasicPhyloXMLparsingValidating() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Roundtrip phyloXML parsing (validating against schema): " );
        if ( Test.testBasicPhyloXMLparsingRoundtrip() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "phyloXML Distribution Element: " );
        if ( Test.testPhyloXMLparsingOfDistributionElement() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Tol XML parsing: " );
        if ( Test.testBasicTolXMLparsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Copying of node data: " );
        if ( Test.testCopyOfNodeData() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Tree copy: " );
        if ( Test.testTreeCopy() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic tree methods: " );
        if ( Test.testBasicTreeMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Tree methods: " );
        if ( Test.testTreeMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Postorder Iterator: " );
        if ( Test.testPostOrderIterator() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Preorder Iterator: " );
        if ( Test.testPreOrderIterator() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Levelorder Iterator: " );
        if ( Test.testLevelOrderIterator() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Re-id methods: " );
        if ( Test.testReIdMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Methods on last external nodes: " );
        if ( Test.testLastExternalNodeMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Methods on external nodes: " );
        if ( Test.testExternalNodeRelatedMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Deletion of external nodes: " );
        if ( Test.testDeletionOfExternalNodes() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Subtree deletion: " );
        if ( Test.testSubtreeDeletion() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Phylogeny branch: " );
        if ( Test.testPhylogenyBranch() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Rerooting: " );
        if ( Test.testRerooting() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Mipoint rooting: " );
        if ( Test.testMidpointrooting() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Node removal: " );
        if ( Test.testNodeRemoval() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Support count: " );
        if ( Test.testSupportCount() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Support transfer: " );
        if ( Test.testSupportTransfer() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Finding of LCA: " );
        if ( Test.testGetLCA() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Finding of LCA 2: " );
        if ( Test.testGetLCA2() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Calculation of distance between nodes: " );
        if ( Test.testGetDistance() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Descriptive statistics: " );
        if ( Test.testDescriptiveStatistics() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Data objects and methods: " );
        if ( Test.testDataObjects() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Properties map: " );
        if ( Test.testPropertiesMap() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "SDIse: " );
        if ( Test.testSDIse() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "SDIunrooted: " );
        if ( Test.testSDIunrooted() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "GSDI: " );
        if ( TestGSDI.test() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "RIO: " );
        if ( TestRIO.test() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Phylogeny reconstruction:" );
        System.out.println();
        if ( TestPhylogenyReconstruction.test( new File( PATH_TO_TEST_DATA ) ) ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Analysis of domain architectures: " );
        System.out.println();
        if ( TestSurfacing.test( new File( PATH_TO_TEST_DATA ) ) ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "GO: " );
        System.out.println();
        if ( TestGo.test( new File( PATH_TO_TEST_DATA ) ) ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Modeling tools: " );
        if ( TestPccx.test() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Split Matrix strict: " );
        if ( Test.testSplitStrict() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Split Matrix: " );
        if ( Test.testSplit() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Confidence Assessor: " );
        if ( Test.testConfidenceAssessor() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic table: " );
        if ( Test.testBasicTable() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "General table: " );
        if ( Test.testGeneralTable() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Amino acid sequence: " );
        if ( Test.testAminoAcidSequence() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "General MSA parser: " );
        if ( Test.testGeneralMsaParser() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Fasta parser for msa: " );
        if ( Test.testFastaParser() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Creation of balanced phylogeny: " );
        if ( Test.testCreateBalancedPhylogeny() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Genbank accessor parsing: " );
        if ( Test.testGenbankAccessorParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        if ( PERFORM_DB_TESTS ) {
            System.out.print( "Uniprot Entry Retrieval: " );
            if ( Test.testUniprotEntryRetrieval() ) {
                System.out.println( "OK." );
                succeeded++;
            }
            else {
                System.out.println( "failed." );
                failed++;
            }
        }
        if ( PERFORM_DB_TESTS ) {
            System.out.print( "Uniprot Taxonomy Search: " );
            if ( Test.testUniprotTaxonomySearch() ) {
                System.out.println( "OK." );
                succeeded++;
            }
            else {
                System.out.println( "failed." );
                failed++;
            }
        }
        //----
        String path = "";
        final String os = ForesterUtil.OS_NAME.toLowerCase();
        if ( ( os.indexOf( "mac" ) >= 0 ) && ( os.indexOf( "os" ) > 0 ) ) {
            path = "/usr/local/bin/mafft";
        }
        else if ( os.indexOf( "win" ) >= 0 ) {
            path = "C:\\Program Files\\mafft-win\\mafft.bat";
        }
        else {
            path = "/home/czmasek/bin/mafft";
        }
        if ( !MsaInferrer.isInstalled( path ) ) {
            path = "mafft";
        }
        if ( !MsaInferrer.isInstalled( path ) ) {
            path = "/usr/local/bin/mafft";
        }
        if ( MsaInferrer.isInstalled( path ) ) {
            System.out.print( "MAFFT (external program): " );
            if ( Test.testMafft( path ) ) {
                System.out.println( "OK." );
                succeeded++;
            }
            else {
                System.out.println( "failed [will not count towards failed tests]" );
            }
        }
        //----
        System.out.print( "Next nodes with collapsed: " );
        if ( Test.testNextNodeWithCollapsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Simple MSA quality: " );
        if ( Test.testMsaQualityMethod() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.println();
        final Runtime rt = java.lang.Runtime.getRuntime();
        final long free_memory = rt.freeMemory() / 1000000;
        final long total_memory = rt.totalMemory() / 1000000;
        System.out.println( "Running time    : " + ( new Date().getTime() - start_time ) + "ms " + "(free memory: "
                + free_memory + "MB, total memory: " + total_memory + "MB)" );
        System.out.println();
        System.out.println( "Successful tests: " + succeeded );
        System.out.println( "Failed     tests: " + failed );
        System.out.println();
        if ( failed < 1 ) {
            System.out.println( "OK." );
        }
        else {
            System.out.println( "Not OK." );
        }
    }

    private final static Phylogeny createPhylogeny( final String nhx ) throws IOException {
        final Phylogeny p = ParserBasedPhylogenyFactory.getInstance().create( nhx, new NHXParser() )[ 0 ];
        return p;
    }

    private final static Event getEvent( final Phylogeny p, final String n1, final String n2 ) {
        return PhylogenyMethods.calculateLCA( p.getNode( n1 ), p.getNode( n2 ) ).getNodeData().getEvent();
    }

    private static boolean testAminoAcidSequence() {
        try {
            final Sequence aa1 = BasicSequence.createAaSequence( "aa1", "aAklm-?xX*z$#" );
            if ( aa1.getLength() != 13 ) {
                return false;
            }
            if ( aa1.getResidueAt( 0 ) != 'A' ) {
                return false;
            }
            if ( aa1.getResidueAt( 2 ) != 'K' ) {
                return false;
            }
            if ( !new String( aa1.getMolecularSequence() ).equals( "AAKLM-XXX*ZXX" ) ) {
                return false;
            }
            final Sequence aa2 = BasicSequence.createAaSequence( "aa3", "ARNDCQEGHILKMFPSTWYVX*-BZOJU" );
            if ( !new String( aa2.getMolecularSequence() ).equals( "ARNDCQEGHILKMFPSTWYVX*-BZXXU" ) ) {
                return false;
            }
            final Sequence dna1 = BasicSequence.createDnaSequence( "dna1", "ACGTUX*-?RYMKWSN" );
            if ( !new String( dna1.getMolecularSequence() ).equals( "ACGTNN*-NRYMKWSN" ) ) {
                return false;
            }
            final Sequence rna1 = BasicSequence.createRnaSequence( "rna1", "..ACGUTX*-?RYMKWSN" );
            if ( !new String( rna1.getMolecularSequence() ).equals( "--ACGUNN*-NRYMKWSN" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testBasicDomain() {
        try {
            final Domain pd = new BasicDomain( "id", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            if ( !pd.getDomainId().equals( "id" ) ) {
                return false;
            }
            if ( pd.getNumber() != 1 ) {
                return false;
            }
            if ( pd.getTotalCount() != 4 ) {
                return false;
            }
            if ( !pd.equals( new BasicDomain( "id", 22, 111, ( short ) 1, ( short ) 4, 0.2, -12 ) ) ) {
                return false;
            }
            final Domain a1 = new BasicDomain( "a", 1, 10, ( short ) 1, ( short ) 4, 0.1, -12 );
            final BasicDomain a1_copy = new BasicDomain( "a", 1, 10, ( short ) 1, ( short ) 4, 0.1, -12 );
            final BasicDomain a1_equal = new BasicDomain( "a", 524, 743994, ( short ) 1, ( short ) 300, 3.0005, 230 );
            final BasicDomain a2 = new BasicDomain( "a", 1, 10, ( short ) 2, ( short ) 4, 0.1, -12 );
            final BasicDomain a3 = new BasicDomain( "A", 1, 10, ( short ) 1, ( short ) 4, 0.1, -12 );
            if ( !a1.equals( a1 ) ) {
                return false;
            }
            if ( !a1.equals( a1_copy ) ) {
                return false;
            }
            if ( !a1.equals( a1_equal ) ) {
                return false;
            }
            if ( !a1.equals( a2 ) ) {
                return false;
            }
            if ( a1.equals( a3 ) ) {
                return false;
            }
            if ( a1.compareTo( a1 ) != 0 ) {
                return false;
            }
            if ( a1.compareTo( a1_copy ) != 0 ) {
                return false;
            }
            if ( a1.compareTo( a1_equal ) != 0 ) {
                return false;
            }
            if ( a1.compareTo( a2 ) != 0 ) {
                return false;
            }
            if ( a1.compareTo( a3 ) == 0 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicNodeMethods() {
        try {
            if ( PhylogenyNode.getNodeCount() != 0 ) {
                return false;
            }
            final PhylogenyNode n1 = new PhylogenyNode();
            final PhylogenyNode n2 = PhylogenyNode
                    .createInstanceFromNhxString( "", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            final PhylogenyNode n3 = PhylogenyNode
                    .createInstanceFromNhxString( "n3", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            final PhylogenyNode n4 = PhylogenyNode
                    .createInstanceFromNhxString( "n4:0.01", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( n1.isHasAssignedEvent() ) {
                return false;
            }
            if ( PhylogenyNode.getNodeCount() != 4 ) {
                return false;
            }
            if ( n3.getIndicator() != 0 ) {
                return false;
            }
            if ( n3.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( !n3.isExternal() ) {
                return false;
            }
            if ( !n3.isRoot() ) {
                return false;
            }
            if ( !n4.getName().equals( "n4" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicPhyloXMLparsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny[] phylogenies_0 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t1.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_0.length != 4 ) {
                return false;
            }
            final Phylogeny t1 = phylogenies_0[ 0 ];
            final Phylogeny t2 = phylogenies_0[ 1 ];
            final Phylogeny t3 = phylogenies_0[ 2 ];
            final Phylogeny t4 = phylogenies_0[ 3 ];
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( !t1.isRooted() ) {
                return false;
            }
            if ( t1.isRerootable() ) {
                return false;
            }
            if ( !t1.getType().equals( "gene_tree" ) ) {
                return false;
            }
            if ( t2.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "node a" ).getDistanceToParent(), 1.0 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "node b" ).getDistanceToParent(), 2.0 ) ) {
                return false;
            }
            if ( t2.getNode( "node a" ).getNodeData().getTaxonomies().size() != 2 ) {
                return false;
            }
            if ( !t2.getNode( "node a" ).getNodeData().getTaxonomy( 0 ).getCommonName().equals( "some parasite" ) ) {
                return false;
            }
            if ( !t2.getNode( "node a" ).getNodeData().getTaxonomy( 1 ).getCommonName().equals( "the host" ) ) {
                return false;
            }
            if ( t2.getNode( "node a" ).getNodeData().getSequences().size() != 2 ) {
                return false;
            }
            if ( !t2.getNode( "node a" ).getNodeData().getSequence( 0 ).getMolecularSequence()
                    .startsWith( "actgtgggggt" ) ) {
                return false;
            }
            if ( !t2.getNode( "node a" ).getNodeData().getSequence( 1 ).getMolecularSequence()
                    .startsWith( "ctgtgatgcat" ) ) {
                return false;
            }
            if ( t3.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            if ( !t1.getName().equals( "t1" ) ) {
                return false;
            }
            if ( !t2.getName().equals( "t2" ) ) {
                return false;
            }
            if ( !t3.getName().equals( "t3" ) ) {
                return false;
            }
            if ( !t4.getName().equals( "t4" ) ) {
                return false;
            }
            if ( !t3.getIdentifier().getValue().equals( "1-1" ) ) {
                return false;
            }
            if ( !t3.getIdentifier().getProvider().equals( "treebank" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getType().equals( "protein" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getName()
                    .equals( "Apoptosis facilitator Bcl-2-like 14 protein" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getSymbol().equals( "BCL2L14" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getAccession().getValue().equals( "Q9BZR8" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getAccession().getSource().equals( "UniProtKB" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getDesc()
                    .equals( "apoptosis" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getRef()
                    .equals( "GO:0006915" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getSource()
                    .equals( "UniProtKB" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getEvidence()
                    .equals( "experimental" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getType()
                    .equals( "function" ) ) {
                return false;
            }
            if ( ( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getConfidence()
                    .getValue() != 1 ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getConfidence()
                    .getType().equals( "ml" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getDesc()
                    .equals( "apoptosis" ) ) {
                return false;
            }
            if ( ( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getAppliesTo() != AppliesTo.ANNOTATION ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getDataType().equals( "xsd:double" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getRef().equals( "AFFY:expression" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getUnit().equals( "AFFY:x" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getValue().equals( "0.2" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "MED:disease" ).getValue().equals( "lymphoma" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 1 ) ).getRef()
                    .equals( "GO:0005829" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) ).getDesc()
                    .equals( "intracellular organelle" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getUri( 0 ).getType().equals( "source" ) ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getUri( 0 ).getDescription()
                    .equals( "UniProt link" ) ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getLocation().equals( "12p13-p12" ) ) ) {
                return false;
            }
            final SortedSet<Accession> x = t3.getNode( "root node" ).getNodeData().getSequence().getCrossReferences();
            if ( x.size() != 4 ) {
                return false;
            }
            int c = 0;
            for( final Accession acc : x ) {
                if ( c == 0 ) {
                    if ( !acc.getSource().equals( "KEGG" ) ) {
                        return false;
                    }
                    if ( !acc.getValue().equals( "hsa:596" ) ) {
                        return false;
                    }
                }
                c++;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicPhyloXMLparsingRoundtrip() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = PhyloXmlParser.createPhyloXmlParser();
            if ( USE_LOCAL_PHYLOXML_SCHEMA ) {
                xml_parser.setValidateAgainstSchema( PHYLOXML_LOCAL_XSD );
            }
            else {
                xml_parser.setValidateAgainstSchema( PHYLOXML_REMOTE_XSD );
            }
            final Phylogeny[] phylogenies_0 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t1.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_0.length != 4 ) {
                return false;
            }
            final StringBuffer t1_sb = new StringBuffer( phylogenies_0[ 0 ].toPhyloXML( 0 ) );
            final Phylogeny[] phylogenies_t1 = factory.create( t1_sb, xml_parser );
            if ( phylogenies_t1.length != 1 ) {
                return false;
            }
            final Phylogeny t1_rt = phylogenies_t1[ 0 ];
            if ( !t1_rt.getDistanceUnit().equals( "cc" ) ) {
                return false;
            }
            if ( !t1_rt.isRooted() ) {
                return false;
            }
            if ( t1_rt.isRerootable() ) {
                return false;
            }
            if ( !t1_rt.getType().equals( "gene_tree" ) ) {
                return false;
            }
            final StringBuffer t2_sb = new StringBuffer( phylogenies_0[ 1 ].toPhyloXML( 0 ) );
            final Phylogeny[] phylogenies_t2 = factory.create( t2_sb, xml_parser );
            final Phylogeny t2_rt = phylogenies_t2[ 0 ];
            if ( t2_rt.getNode( "node a" ).getNodeData().getTaxonomies().size() != 2 ) {
                return false;
            }
            if ( !t2_rt.getNode( "node a" ).getNodeData().getTaxonomy( 0 ).getCommonName().equals( "some parasite" ) ) {
                return false;
            }
            if ( !t2_rt.getNode( "node a" ).getNodeData().getTaxonomy( 1 ).getCommonName().equals( "the host" ) ) {
                return false;
            }
            if ( t2_rt.getNode( "node a" ).getNodeData().getSequences().size() != 2 ) {
                return false;
            }
            if ( !t2_rt.getNode( "node a" ).getNodeData().getSequence( 0 ).getMolecularSequence()
                    .startsWith( "actgtgggggt" ) ) {
                return false;
            }
            if ( !t2_rt.getNode( "node a" ).getNodeData().getSequence( 1 ).getMolecularSequence()
                    .startsWith( "ctgtgatgcat" ) ) {
                return false;
            }
            final StringBuffer t3_sb_0 = new StringBuffer( phylogenies_0[ 2 ].toPhyloXML( 0 ) );
            final Phylogeny[] phylogenies_1_0 = factory.create( t3_sb_0, xml_parser );
            final StringBuffer t3_sb = new StringBuffer( phylogenies_1_0[ 0 ].toPhyloXML( 0 ) );
            final Phylogeny[] phylogenies_1 = factory.create( t3_sb, xml_parser );
            if ( phylogenies_1.length != 1 ) {
                return false;
            }
            final Phylogeny t3_rt = phylogenies_1[ 0 ];
            if ( !t3_rt.getName().equals( "t3" ) ) {
                return false;
            }
            if ( t3_rt.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            if ( !t3_rt.getIdentifier().getValue().equals( "1-1" ) ) {
                return false;
            }
            if ( !t3_rt.getIdentifier().getProvider().equals( "treebank" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getType().equals( "protein" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getName()
                    .equals( "Apoptosis facilitator Bcl-2-like 14 protein" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getSymbol().equals( "BCL2L14" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getAccession().getValue().equals( "Q9BZR8" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getAccession().getSource()
                    .equals( "UniProtKB" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getDesc()
                    .equals( "apoptosis" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getRef()
                    .equals( "GO:0006915" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getSource()
                    .equals( "UniProtKB" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getEvidence()
                    .equals( "experimental" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getType()
                    .equals( "function" ) ) {
                return false;
            }
            if ( ( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getConfidence()
                    .getValue() != 1 ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getConfidence()
                    .getType().equals( "ml" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getDesc()
                    .equals( "apoptosis" ) ) {
                return false;
            }
            if ( ( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getAppliesTo() != AppliesTo.ANNOTATION ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getDataType().equals( "xsd:double" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getRef().equals( "AFFY:expression" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getUnit().equals( "AFFY:x" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "AFFY:expression" ).getValue().equals( "0.2" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getProperties()
                    .getProperty( "MED:disease" ).getValue().equals( "lymphoma" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 1 ) ).getRef()
                    .equals( "GO:0005829" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) ).getDesc()
                    .equals( "intracellular organelle" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getUri( 0 ).getType().equals( "source" ) ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getUri( 0 ).getDescription()
                    .equals( "UniProt link" ) ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getLocation().equals( "12p13-p12" ) ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getReference().getDoi().equals( "10.1038/387489a0" ) ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getReference().getDescription()
                    .equals( "Aguinaldo, A. M. A.; J. M. Turbeville, L. S. Linford, M. C. Rivera, J. R. Garey, R. A. Raff, & J. A. Lake (1997). \"Evidence for a clade of nematodes, arthropods and other moulting animals\". Nature 387 (6632): 489–493." ) ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getTaxonomyCode().equals( "ECDYS" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getScientificName().equals( "ecdysozoa" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getCommonName().equals( "molting animals" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getIdentifier().getValue().equals( "1" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getIdentifier().getProvider()
                    .equals( "ncbi" ) ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getTotalLength() != 124 ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 )
                    .getName().equals( "B" ) ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 )
                    .getFrom() != 21 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 ).getTo() != 44 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 )
                    .getLength() != 24 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 )
                    .getConfidence() != 2144 ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 ).getId()
                    .equals( "pfam" ) ) {
                return false;
            }
            if ( t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().getGainedCharacters().size() != 3 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().getPresentCharacters().size() != 2 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().getLostCharacters().size() != 1 ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().getType().equals( "domains" ) ) {
                return false;
            }
            final Taxonomy taxbb = t3_rt.getNode( "node bb" ).getNodeData().getTaxonomy();
            if ( !taxbb.getAuthority().equals( "Stephenson, 1935" ) ) {
                return false;
            }
            if ( !taxbb.getCommonName().equals( "starlet sea anemone" ) ) {
                return false;
            }
            if ( !taxbb.getIdentifier().getProvider().equals( "EOL" ) ) {
                return false;
            }
            if ( !taxbb.getIdentifier().getValue().equals( "704294" ) ) {
                return false;
            }
            if ( !taxbb.getTaxonomyCode().equals( "NEMVE" ) ) {
                return false;
            }
            if ( !taxbb.getScientificName().equals( "Nematostella vectensis" ) ) {
                return false;
            }
            if ( taxbb.getSynonyms().size() != 2 ) {
                return false;
            }
            if ( !taxbb.getSynonyms().contains( "Nematostella vectensis Stephenson1935" ) ) {
                return false;
            }
            if ( !taxbb.getSynonyms().contains( "See Anemone" ) ) {
                return false;
            }
            if ( !taxbb.getUri( 0 ).getDescription().equals( "EOL" ) ) {
                return false;
            }
            if ( !taxbb.getUri( 0 ).getType().equals( "linkout" ) ) {
                return false;
            }
            if ( !taxbb.getUri( 0 ).getValue().toString().equals( "http://www.eol.org/pages/704294" ) ) {
                return false;
            }
            if ( ( ( BinaryCharacters ) t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().copy() )
                    .getLostCount() != BinaryCharacters.COUNT_DEFAULT ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getGainedCount() != 1 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getGainedCharacters().size() != 1 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getLostCount() != 3 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getLostCharacters().size() != 3 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getPresentCount() != 2 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getPresentCharacters().size() != 2 ) {
                return false;
            }
            if ( !t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getType().equals( "characters" ) ) {
                return false;
            }
            //
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getDesc().equals( "Silurian" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getValue().toPlainString()
                    .equalsIgnoreCase( "435" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getMin().toPlainString().equalsIgnoreCase( "416" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getMax().toPlainString()
                    .equalsIgnoreCase( "443.7" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getUnit().equals( "mya" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bb" ).getNodeData().getDate().getDesc().equals( "Triassic" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bc" ).getNodeData().getDate().getValue().toPlainString()
                    .equalsIgnoreCase( "433" ) ) {
                return false;
            }
            final SortedSet<Accession> x = t3_rt.getNode( "root node" ).getNodeData().getSequence()
                    .getCrossReferences();
            if ( x.size() != 4 ) {
                return false;
            }
            int c = 0;
            for( final Accession acc : x ) {
                if ( c == 0 ) {
                    if ( !acc.getSource().equals( "KEGG" ) ) {
                        return false;
                    }
                    if ( !acc.getValue().equals( "hsa:596" ) ) {
                        return false;
                    }
                }
                c++;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicPhyloXMLparsingValidating() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            PhyloXmlParser xml_parser = null;
            try {
                xml_parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
            }
            catch ( final Exception e ) {
                // Do nothing -- means were not running from jar.
            }
            if ( xml_parser == null ) {
                xml_parser = PhyloXmlParser.createPhyloXmlParser();
                if ( USE_LOCAL_PHYLOXML_SCHEMA ) {
                    xml_parser.setValidateAgainstSchema( PHYLOXML_LOCAL_XSD );
                }
                else {
                    xml_parser.setValidateAgainstSchema( PHYLOXML_REMOTE_XSD );
                }
            }
            final Phylogeny[] phylogenies_0 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t1.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_0.length != 4 ) {
                return false;
            }
            final Phylogeny t1 = phylogenies_0[ 0 ];
            final Phylogeny t2 = phylogenies_0[ 1 ];
            final Phylogeny t3 = phylogenies_0[ 2 ];
            final Phylogeny t4 = phylogenies_0[ 3 ];
            if ( !t1.getName().equals( "t1" ) ) {
                return false;
            }
            if ( !t2.getName().equals( "t2" ) ) {
                return false;
            }
            if ( !t3.getName().equals( "t3" ) ) {
                return false;
            }
            if ( !t4.getName().equals( "t4" ) ) {
                return false;
            }
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( t2.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            if ( t3.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            final String x2 = Test.PATH_TO_TEST_DATA + "phyloxml_test_t1.xml";
            final Phylogeny[] phylogenies_1 = factory.create( x2, xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( "errors:" );
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_1.length != 4 ) {
                return false;
            }
            final Phylogeny[] phylogenies_2 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t3.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( "errors:" );
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_2.length != 1 ) {
                return false;
            }
            if ( phylogenies_2[ 0 ].getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            final Phylogeny[] phylogenies_3 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t4.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_3.length != 2 ) {
                return false;
            }
            final Phylogeny a = phylogenies_3[ 0 ];
            if ( !a.getName().equals( "tree 4" ) ) {
                return false;
            }
            if ( a.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !a.getNode( "node b1" ).getNodeData().getSequence().getName().equals( "b1 gene" ) ) {
                return false;
            }
            if ( !a.getNode( "node b1" ).getNodeData().getTaxonomy().getCommonName().equals( "b1 species" ) ) {
                return false;
            }
            final Phylogeny[] phylogenies_4 = factory.create( Test.PATH_TO_TEST_DATA + "special_characters.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_4.length != 1 ) {
                return false;
            }
            final Phylogeny s = phylogenies_4[ 0 ];
            if ( s.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            s.getNode( "first" );
            s.getNode( "<>" );
            s.getNode( "\"<a'b&c'd\">\"" );
            s.getNode( "'''\"" );
            s.getNode( "\"\"\"" );
            s.getNode( "dick & doof" );
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicProtein() {
        try {
            final BasicProtein p0 = new BasicProtein( "p0", "owl", 0 );
            final Domain a = new BasicDomain( "a", 1, 10, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain b = new BasicDomain( "b", 11, 20, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain c = new BasicDomain( "c", 9, 23, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain d = new BasicDomain( "d", 15, 30, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain e = new BasicDomain( "e", 60, 70, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain x = new BasicDomain( "x", 100, 110, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain y = new BasicDomain( "y", 100, 110, ( short ) 1, ( short ) 5, 0.1, -12 );
            p0.addProteinDomain( y );
            p0.addProteinDomain( e );
            p0.addProteinDomain( b );
            p0.addProteinDomain( c );
            p0.addProteinDomain( d );
            p0.addProteinDomain( a );
            p0.addProteinDomain( x );
            if ( !p0.toDomainArchitectureString( "~" ).equals( "a~b~c~d~e~x~y" ) ) {
                return false;
            }
            if ( !p0.toDomainArchitectureString( "~", 3, "=" ).equals( "a~b~c~d~e~x~y" ) ) {
                return false;
            }
            //
            final BasicProtein aa0 = new BasicProtein( "aa", "owl", 0 );
            final Domain a1 = new BasicDomain( "a", 1, 10, ( short ) 1, ( short ) 5, 0.1, -12 );
            aa0.addProteinDomain( a1 );
            if ( !aa0.toDomainArchitectureString( "~" ).equals( "a" ) ) {
                return false;
            }
            if ( !aa0.toDomainArchitectureString( "~", 3, "" ).equals( "a" ) ) {
                return false;
            }
            //
            final BasicProtein aa1 = new BasicProtein( "aa", "owl", 0 );
            final Domain a11 = new BasicDomain( "a", 1, 10, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain a12 = new BasicDomain( "a", 2, 20, ( short ) 1, ( short ) 5, 0.1, -12 );
            aa1.addProteinDomain( a11 );
            aa1.addProteinDomain( a12 );
            if ( !aa1.toDomainArchitectureString( "~" ).equals( "a~a" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 3, "" ).equals( "a~a" ) ) {
                return false;
            }
            aa1.addProteinDomain( new BasicDomain( "a", 20, 30, ( short ) 1, ( short ) 5, 0.1, -12 ) );
            if ( !aa1.toDomainArchitectureString( "~" ).equals( "a~a~a" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 3, "" ).equals( "aaa" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 4, "" ).equals( "a~a~a" ) ) {
                return false;
            }
            aa1.addProteinDomain( new BasicDomain( "a", 30, 40, ( short ) 1, ( short ) 5, 0.1, -12 ) );
            if ( !aa1.toDomainArchitectureString( "~" ).equals( "a~a~a~a" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 3, "" ).equals( "aaa" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 4, "" ).equals( "aaa" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 5, "" ).equals( "a~a~a~a" ) ) {
                return false;
            }
            aa1.addProteinDomain( new BasicDomain( "b", 32, 40, ( short ) 1, ( short ) 5, 0.1, -12 ) );
            if ( !aa1.toDomainArchitectureString( "~" ).equals( "a~a~a~a~b" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 3, "" ).equals( "aaa~b" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 4, "" ).equals( "aaa~b" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 5, "" ).equals( "a~a~a~a~b" ) ) {
                return false;
            }
            aa1.addProteinDomain( new BasicDomain( "c", 1, 2, ( short ) 1, ( short ) 5, 0.1, -12 ) );
            if ( !aa1.toDomainArchitectureString( "~" ).equals( "c~a~a~a~a~b" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 3, "" ).equals( "c~aaa~b" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 4, "" ).equals( "c~aaa~b" ) ) {
                return false;
            }
            if ( !aa1.toDomainArchitectureString( "~", 5, "" ).equals( "c~a~a~a~a~b" ) ) {
                return false;
            }
            //
            final BasicProtein p00 = new BasicProtein( "p0", "owl", 0 );
            final Domain a0 = new BasicDomain( "a", 1, 10, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain b0 = new BasicDomain( "b", 11, 20, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain c0 = new BasicDomain( "c", 9, 23, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain d0 = new BasicDomain( "d", 15, 30, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain e0 = new BasicDomain( "e", 60, 70, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain e1 = new BasicDomain( "e", 61, 71, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain e2 = new BasicDomain( "e", 62, 72, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain e3 = new BasicDomain( "e", 63, 73, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain e4 = new BasicDomain( "e", 64, 74, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain e5 = new BasicDomain( "e", 65, 75, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain x0 = new BasicDomain( "x", 100, 110, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain y0 = new BasicDomain( "y", 100, 110, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain y1 = new BasicDomain( "y", 120, 130, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain y2 = new BasicDomain( "y", 140, 150, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain y3 = new BasicDomain( "y", 160, 170, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain z0 = new BasicDomain( "z", 200, 210, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain z1 = new BasicDomain( "z", 300, 310, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain z2 = new BasicDomain( "z", 400, 410, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain zz0 = new BasicDomain( "Z", 500, 510, ( short ) 1, ( short ) 5, 0.1, -12 );
            final Domain zz1 = new BasicDomain( "Z", 600, 610, ( short ) 1, ( short ) 5, 0.1, -12 );
            p00.addProteinDomain( y0 );
            p00.addProteinDomain( e0 );
            p00.addProteinDomain( b0 );
            p00.addProteinDomain( c0 );
            p00.addProteinDomain( d0 );
            p00.addProteinDomain( a0 );
            p00.addProteinDomain( x0 );
            p00.addProteinDomain( y1 );
            p00.addProteinDomain( y2 );
            p00.addProteinDomain( y3 );
            p00.addProteinDomain( e1 );
            p00.addProteinDomain( e2 );
            p00.addProteinDomain( e3 );
            p00.addProteinDomain( e4 );
            p00.addProteinDomain( e5 );
            p00.addProteinDomain( z0 );
            p00.addProteinDomain( z1 );
            p00.addProteinDomain( z2 );
            p00.addProteinDomain( zz0 );
            p00.addProteinDomain( zz1 );
            if ( !p00.toDomainArchitectureString( "~", 3, "" ).equals( "a~b~c~d~eee~x~yyy~zzz~Z~Z" ) ) {
                return false;
            }
            if ( !p00.toDomainArchitectureString( "~", 4, "" ).equals( "a~b~c~d~eee~x~yyy~z~z~z~Z~Z" ) ) {
                return false;
            }
            if ( !p00.toDomainArchitectureString( "~", 5, "" ).equals( "a~b~c~d~eee~x~y~y~y~y~z~z~z~Z~Z" ) ) {
                return false;
            }
            if ( !p00.toDomainArchitectureString( "~", 6, "" ).equals( "a~b~c~d~eee~x~y~y~y~y~z~z~z~Z~Z" ) ) {
                return false;
            }
            if ( !p00.toDomainArchitectureString( "~", 7, "" ).equals( "a~b~c~d~e~e~e~e~e~e~x~y~y~y~y~z~z~z~Z~Z" ) ) {
                return false;
            }
            // A0  A10  B15  A20  B25  A30  B35  B40  C50  A60  C70  D80
            final Domain A0 = new BasicDomain( "A", 0, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain A10 = new BasicDomain( "A", 10, 11, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain B15 = new BasicDomain( "B", 11, 16, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain A20 = new BasicDomain( "A", 20, 100, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain B25 = new BasicDomain( "B", 25, 26, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain A30 = new BasicDomain( "A", 30, 31, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain B35 = new BasicDomain( "B", 31, 40, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain B40 = new BasicDomain( "B", 40, 600, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain C50 = new BasicDomain( "C", 50, 59, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain A60 = new BasicDomain( "A", 60, 395, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain C70 = new BasicDomain( "C", 70, 71, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain D80 = new BasicDomain( "D", 80, 81, ( short ) 1, ( short ) 4, 0.1, -12 );
            final BasicProtein p = new BasicProtein( "p", "owl", 0 );
            p.addProteinDomain( B15 );
            p.addProteinDomain( C50 );
            p.addProteinDomain( A60 );
            p.addProteinDomain( A30 );
            p.addProteinDomain( C70 );
            p.addProteinDomain( B35 );
            p.addProteinDomain( B40 );
            p.addProteinDomain( A0 );
            p.addProteinDomain( A10 );
            p.addProteinDomain( A20 );
            p.addProteinDomain( B25 );
            p.addProteinDomain( D80 );
            List<String> domains_ids = new ArrayList<String>();
            domains_ids.add( "A" );
            domains_ids.add( "B" );
            domains_ids.add( "C" );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids.add( "X" );
            if ( p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<String>();
            domains_ids.add( "A" );
            domains_ids.add( "C" );
            domains_ids.add( "D" );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<String>();
            domains_ids.add( "A" );
            domains_ids.add( "D" );
            domains_ids.add( "C" );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<String>();
            domains_ids.add( "A" );
            domains_ids.add( "A" );
            domains_ids.add( "B" );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<String>();
            domains_ids.add( "A" );
            domains_ids.add( "A" );
            domains_ids.add( "A" );
            domains_ids.add( "B" );
            domains_ids.add( "B" );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<String>();
            domains_ids.add( "A" );
            domains_ids.add( "A" );
            domains_ids.add( "B" );
            domains_ids.add( "A" );
            domains_ids.add( "B" );
            domains_ids.add( "B" );
            domains_ids.add( "A" );
            domains_ids.add( "B" );
            domains_ids.add( "C" );
            domains_ids.add( "A" );
            domains_ids.add( "C" );
            domains_ids.add( "D" );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicTable() {
        try {
            final BasicTable<String> t0 = new BasicTable<String>();
            if ( t0.getNumberOfColumns() != 0 ) {
                return false;
            }
            if ( t0.getNumberOfRows() != 0 ) {
                return false;
            }
            t0.setValue( 3, 2, "23" );
            t0.setValue( 10, 1, "error" );
            t0.setValue( 10, 1, "110" );
            t0.setValue( 9, 1, "19" );
            t0.setValue( 1, 10, "101" );
            t0.setValue( 10, 10, "1010" );
            t0.setValue( 100, 10, "10100" );
            t0.setValue( 0, 0, "00" );
            if ( !t0.getValue( 3, 2 ).equals( "23" ) ) {
                return false;
            }
            if ( !t0.getValue( 10, 1 ).equals( "110" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 1, 10 ).equals( "101" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 10, 10 ).equals( "1010" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 100, 10 ).equals( "10100" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 9, 1 ).equals( "19" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( t0.getNumberOfColumns() != 101 ) {
                return false;
            }
            if ( t0.getNumberOfRows() != 11 ) {
                return false;
            }
            if ( t0.getValueAsString( 49, 4 ) != null ) {
                return false;
            }
            final String l = ForesterUtil.getLineSeparator();
            final StringBuffer source = new StringBuffer();
            source.append( "" + l );
            source.append( "# 1 1 1 1 1 1 1 1" + l );
            source.append( " 00 01 02 03" + l );
            source.append( "   10 11 12 13  " + l );
            source.append( "20 21 22 23 " + l );
            source.append( "    30  31    32 33" + l );
            source.append( "40 41 42 43" + l );
            source.append( "  # 1 1 1 1 1 " + l );
            source.append( "50 51 52 53 54" + l );
            final BasicTable<String> t1 = BasicTableParser.parse( source.toString(), ' ' );
            if ( t1.getNumberOfColumns() != 5 ) {
                return false;
            }
            if ( t1.getNumberOfRows() != 6 ) {
                return false;
            }
            if ( !t1.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( 1, 0 ).equals( "01" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( 3, 0 ).equals( "03" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( 4, 5 ).equals( "54" ) ) {
                return false;
            }
            final StringBuffer source1 = new StringBuffer();
            source1.append( "" + l );
            source1.append( "# 1; 1; 1; 1 ;1 ;1; 1 ;1;" + l );
            source1.append( " 00; 01 ;02;03" + l );
            source1.append( "   10; 11; 12; 13  " + l );
            source1.append( "20; 21; 22; 23 " + l );
            source1.append( "    30;  31;    32; 33" + l );
            source1.append( "40;41;42;43" + l );
            source1.append( "  # 1 1 1 1 1 " + l );
            source1.append( ";;;50  ;  ;52; 53;;54   " + l );
            final BasicTable<String> t2 = BasicTableParser.parse( source1.toString(), ';' );
            if ( t2.getNumberOfColumns() != 5 ) {
                return false;
            }
            if ( t2.getNumberOfRows() != 6 ) {
                return false;
            }
            if ( !t2.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 1, 0 ).equals( "01" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 3, 0 ).equals( "03" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 3, 3 ).equals( "33" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 3, 5 ).equals( "53" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 1, 5 ).equals( "" ) ) {
                return false;
            }
            final StringBuffer source2 = new StringBuffer();
            source2.append( "" + l );
            source2.append( "comment: 1; 1; 1; 1 ;1 ;1; 1 ;1;" + l );
            source2.append( " 00; 01 ;02;03" + l );
            source2.append( "   10; 11; 12; 13  " + l );
            source2.append( "20; 21; 22; 23 " + l );
            source2.append( "                     " + l );
            source2.append( "    30;  31;    32; 33" + l );
            source2.append( "40;41;42;43" + l );
            source2.append( "  comment: 1 1 1 1 1 " + l );
            source2.append( ";;;50  ;   52; 53;;54   " + l );
            final List<BasicTable<String>> tl = BasicTableParser.parse( source2.toString(),
                                                                        ';',
                                                                        false,
                                                                        false,
                                                                        "comment:",
                                                                        false );
            if ( tl.size() != 2 ) {
                return false;
            }
            final BasicTable<String> t3 = tl.get( 0 );
            final BasicTable<String> t4 = tl.get( 1 );
            if ( t3.getNumberOfColumns() != 4 ) {
                return false;
            }
            if ( t3.getNumberOfRows() != 3 ) {
                return false;
            }
            if ( t4.getNumberOfColumns() != 4 ) {
                return false;
            }
            if ( t4.getNumberOfRows() != 3 ) {
                return false;
            }
            if ( !t3.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( !t4.getValueAsString( 0, 0 ).equals( "30" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicTolXMLparsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final TolParser parser = new TolParser();
            final Phylogeny[] phylogenies_0 = factory.create( Test.PATH_TO_TEST_DATA + "tol_2484.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_0.length != 1 ) {
                return false;
            }
            final Phylogeny t1 = phylogenies_0[ 0 ];
            if ( t1.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            if ( !t1.isRooted() ) {
                return false;
            }
            if ( !t1.getRoot().getNodeData().getTaxonomy().getScientificName().equals( "Mesozoa" ) ) {
                return false;
            }
            if ( !t1.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "2484" ) ) {
                return false;
            }
            if ( !t1.getRoot().getChildNode( 0 ).getNodeData().getTaxonomy().getScientificName().equals( "Rhombozoa" ) ) {
                return false;
            }
            if ( t1.getRoot().getChildNode( 0 ).getNumberOfDescendants() != 3 ) {
                return false;
            }
            final Phylogeny[] phylogenies_1 = factory.create( Test.PATH_TO_TEST_DATA + "tol_2.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_1.length != 1 ) {
                return false;
            }
            final Phylogeny t2 = phylogenies_1[ 0 ];
            if ( t2.getNumberOfExternalNodes() != 664 ) {
                return false;
            }
            if ( !t2.isRooted() ) {
                return false;
            }
            if ( !t2.getRoot().getNodeData().getTaxonomy().getScientificName().equals( "Eubacteria" ) ) {
                return false;
            }
            if ( !t2.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "2" ) ) {
                return false;
            }
            if ( t2.getRoot().getNumberOfDescendants() != 24 ) {
                return false;
            }
            if ( t2.getRoot().getNumberOfDescendants() != 24 ) {
                return false;
            }
            if ( !t2.getRoot().getChildNode( 0 ).getNodeData().getTaxonomy().getScientificName().equals( "Aquificae" ) ) {
                return false;
            }
            if ( !t2.getRoot().getChildNode( 0 ).getChildNode( 0 ).getNodeData().getTaxonomy().getScientificName()
                    .equals( "Aquifex" ) ) {
                return false;
            }
            final Phylogeny[] phylogenies_2 = factory.create( Test.PATH_TO_TEST_DATA + "tol_5.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_2.length != 1 ) {
                return false;
            }
            final Phylogeny t3 = phylogenies_2[ 0 ];
            if ( t3.getNumberOfExternalNodes() != 184 ) {
                return false;
            }
            if ( !t3.getRoot().getNodeData().getTaxonomy().getScientificName().equals( "Viruses" ) ) {
                return false;
            }
            if ( !t3.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "5" ) ) {
                return false;
            }
            if ( t3.getRoot().getNumberOfDescendants() != 6 ) {
                return false;
            }
            final Phylogeny[] phylogenies_3 = factory.create( Test.PATH_TO_TEST_DATA + "tol_4567.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_3.length != 1 ) {
                return false;
            }
            final Phylogeny t4 = phylogenies_3[ 0 ];
            if ( t4.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( !t4.getRoot().getNodeData().getTaxonomy().getScientificName().equals( "Marpissa decorata" ) ) {
                return false;
            }
            if ( !t4.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "4567" ) ) {
                return false;
            }
            if ( t4.getRoot().getNumberOfDescendants() != 0 ) {
                return false;
            }
            final Phylogeny[] phylogenies_4 = factory.create( Test.PATH_TO_TEST_DATA + "tol_16299.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_4.length != 1 ) {
                return false;
            }
            final Phylogeny t5 = phylogenies_4[ 0 ];
            if ( t5.getNumberOfExternalNodes() != 13 ) {
                return false;
            }
            if ( !t5.getRoot().getNodeData().getTaxonomy().getScientificName().equals( "Hominidae" ) ) {
                return false;
            }
            if ( !t5.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "16299" ) ) {
                return false;
            }
            if ( t5.getRoot().getNumberOfDescendants() != 2 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicTreeMethods() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t1 = factory.create();
            if ( !t1.isEmpty() ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "((A:1,B:2)AB:1,(C:3,D:5)CD:3)ABCD:0.5", new NHXParser() )[ 0 ];
            if ( t2.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            if ( t2.getHeight() != 8.5 ) {
                return false;
            }
            if ( !t2.isCompletelyBinary() ) {
                return false;
            }
            if ( t2.isEmpty() ) {
                return false;
            }
            final Phylogeny t3 = factory.create( "((A:1,B:2,C:10)ABC:1,(D:3,E:5)DE:3)", new NHXParser() )[ 0 ];
            if ( t3.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            if ( t3.getHeight() != 11 ) {
                return false;
            }
            if ( t3.isCompletelyBinary() ) {
                return false;
            }
            final PhylogenyNode n = t3.getNode( "ABC" );
            final Phylogeny t4 = factory.create( "((A:1,B:2,C:10)ABC:1,(D:3,E:5)DE:3,(F,G,H,I))", new NHXParser() )[ 0 ];
            if ( t4.getNumberOfExternalNodes() != 9 ) {
                return false;
            }
            if ( t4.getHeight() != 11 ) {
                return false;
            }
            if ( t4.isCompletelyBinary() ) {
                return false;
            }
            final StringBuffer sb5 = new StringBuffer( "(((A11:2)A1:2,(A21:1,A22:2,A23)A2:11,A3:2)A:2,B:10,C:3,D:8)" );
            final Phylogeny t5 = factory.create( sb5, new NHXParser() )[ 0 ];
            if ( t5.getNumberOfExternalNodes() != 8 ) {
                return false;
            }
            if ( t5.getHeight() != 15 ) {
                return false;
            }
            final StringBuffer sb6 = new StringBuffer( "(X,Y,Z,(((A111)A11:2)A1:2,(X,Y,Z,A21:1,A22:2,A23)A2:11,A3:2)A:2,B:10,C:3,D:8)" );
            final Phylogeny t6 = factory.create( sb6, new NHXParser() )[ 0 ];
            if ( t6.getHeight() != 15 ) {
                return false;
            }
            final StringBuffer sb7 = new StringBuffer( "(((A11:2)A1:2,(A21:1,A22:2,A23)A2:11,A3:2)A:2,B:10,C:15,D:8)" );
            final Phylogeny t7 = factory.create( sb7, new NHXParser() )[ 0 ];
            if ( t7.getHeight() != 15 ) {
                return false;
            }
            final StringBuffer sb8 = new StringBuffer( "(((A11:11)A1:2,(A21:2,A22:2,A23,A24,AA:)A2:11,A3:2)A:2,B:15,C:15,D:15)" );
            final Phylogeny t8 = factory.create( sb8, new NHXParser() )[ 0 ];
            if ( t8.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( t8.getHeight() != 15 ) {
                return false;
            }
            final char[] a9 = new char[] { 'a' };
            final Phylogeny t9 = factory.create( a9, new NHXParser() )[ 0 ];
            if ( t9.getHeight() != 0 ) {
                return false;
            }
            final char[] a10 = new char[] { 'a', ':', '6' };
            final Phylogeny t10 = factory.create( a10, new NHXParser() )[ 0 ];
            if ( t10.getHeight() != 6 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testConfidenceAssessor() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((((A,B)ab,C)abc,D)abcd,E)abcde", new NHXParser() )[ 0 ];
            final Phylogeny[] ev0 = factory
                    .create( "((((A,B),C),D),E);((((A,B),C),D),E);((((A,B),C),D),E);((((A,B),C),D),E);",
                             new NHXParser() );
            ConfidenceAssessor.evaluate( "bootstrap", ev0, t0, false, 1, 0, 2 );
            if ( !isEqual( t0.getNode( "ab" ).getBranchData().getConfidence( 0 ).getValue(), 3 ) ) {
                return false;
            }
            if ( !isEqual( t0.getNode( "abc" ).getBranchData().getConfidence( 0 ).getValue(), 3 ) ) {
                return false;
            }
            final Phylogeny t1 = factory.create( "((((A,B)ab[&&NHX:B=50],C)abc,D)abcd,E)abcde", new NHXParser() )[ 0 ];
            final Phylogeny[] ev1 = factory
                    .create( "((((A,B),C),D),E);((A,B),((E,D),C));(((A,B),C),(E,D));(A,(((E,D),C),B));(B,(A,((E,D),C)));(C,((E,D),(A,B)));(D,(E,((A,B),C)));",
                             new NHXParser() );
            ConfidenceAssessor.evaluate( "bootstrap", ev1, t1, false, 1 );
            if ( !isEqual( t1.getNode( "ab" ).getBranchData().getConfidence( 1 ).getValue(), 7 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "abc" ).getBranchData().getConfidence( 0 ).getValue(), 7 ) ) {
                return false;
            }
            final Phylogeny t_b = factory.create( "((((A,C)ac,D)acd,E)acde,B)abcde", new NHXParser() )[ 0 ];
            final Phylogeny[] ev_b = factory
                    .create( "((A,C),X);((A,X),C);(A,C);((((A,B),C),D),E);((A,B),((E,D),C));(((A,B),C),(E,D));(A,(((E,D),C),B));(B,(A,((E,D),C)));(C,((E,D),(A,B)));(D,(E,((A,B),C)));((((A,C)ac,D)acd,E)acde,B)abcd",
                             new NHXParser() );
            ConfidenceAssessor.evaluate( "bootstrap", ev_b, t_b, false, 1 );
            if ( !isEqual( t_b.getNode( "ac" ).getBranchData().getConfidence( 0 ).getValue(), 4 ) ) {
                return false;
            }
            if ( !isEqual( t_b.getNode( "acd" ).getBranchData().getConfidence( 0 ).getValue(), 1 ) ) {
                return false;
            }
            //
            final Phylogeny t1x = factory.create( "((((A,B)ab,C)abc,D)abcd,E)abcde", new NHXParser() )[ 0 ];
            final Phylogeny[] ev1x = factory
                    .create( "((((A,B),C),D),E);((A,B),((E,D),C));(((A,B),C),(E,D));(A,(((E,D),C),B));(B,(A,((E,D),C)));(C,((E,D),(A,B)));(D,(E,((A,B),C)));",
                             new NHXParser() );
            ConfidenceAssessor.evaluate( "bootstrap", ev1x, t1x, true, 1 );
            if ( !isEqual( t1x.getNode( "ab" ).getBranchData().getConfidence( 0 ).getValue(), 7 ) ) {
                return false;
            }
            if ( !isEqual( t1x.getNode( "abc" ).getBranchData().getConfidence( 0 ).getValue(), 7 ) ) {
                return false;
            }
            final Phylogeny t_bx = factory.create( "((((A,C)ac,D)acd,E)acde,B)abcde", new NHXParser() )[ 0 ];
            final Phylogeny[] ev_bx = factory
                    .create( "((((A,B),C),D),E);((A,B),((E,D),C));(((A,B),C),(E,D));(A,(((E,D),C),B));(B,(A,((E,D),C)));(C,((E,D),(A,B)));(D,(E,((A,B),C)));((((A,C)ac,D)acd,E)acde,B)abcd",
                             new NHXParser() );
            ConfidenceAssessor.evaluate( "bootstrap", ev_bx, t_bx, true, 1 );
            if ( !isEqual( t_bx.getNode( "ac" ).getBranchData().getConfidence( 0 ).getValue(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t_bx.getNode( "acd" ).getBranchData().getConfidence( 0 ).getValue(), 1 ) ) {
                return false;
            }
            //
            final Phylogeny[] t2 = factory
                    .create( "((((a,b),c),d),e);(((a,b),c),(d,e));(((((a,b),c),d),e),f);((((a,b),c),(d,e)),f);(((a,b),c),d,e);((a,b,c),d,e);",
                             new NHXParser() );
            final Phylogeny[] ev2 = factory
                    .create( "((((a,b),c),d),e);((((a,b),c),d),e);((((a,b),e),d),c);((((a,b),e),d),c);(((a,b),(c,d)),e);((a,b),x);((a,b),(x,y));(a,b);(a,e);(a,b,c);",
                             new NHXParser() );
            for( final Phylogeny target : t2 ) {
                ConfidenceAssessor.evaluate( "bootstrap", ev2, target, false, 1 );
            }
            //
            final Phylogeny t4 = factory.create( "((((((A,B)ab,C)abc,D)abcd,E)abcde,F)abcdef,G)abcdefg",
                                                 new NHXParser() )[ 0 ];
            final Phylogeny[] ev4 = factory.create( "(((A,B),C),(X,Y));((F,G),((A,B,C),(D,E)))", new NHXParser() );
            ConfidenceAssessor.evaluate( "bootstrap", ev4, t4, false, 1 );
            if ( !isEqual( t4.getNode( "ab" ).getBranchData().getConfidence( 0 ).getValue(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t4.getNode( "abc" ).getBranchData().getConfidence( 0 ).getValue(), 2 ) ) {
                return false;
            }
            if ( !isEqual( t4.getNode( "abcde" ).getBranchData().getConfidence( 0 ).getValue(), 1 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testCopyOfNodeData() {
        try {
            final PhylogenyNode n1 = PhylogenyNode
                    .createInstanceFromNhxString( "n5:0.1[&&NHX:S=Ecoli:E=1.1.1.1:D=Y:Co=Y:B=56:T=1:O=22:SO=33:SN=44:W=2:C=10.20.30:XN=S=tag1=value1=unit1]" );
            final PhylogenyNode n2 = n1.copyNodeData();
            if ( !n1.toNewHampshireX().equals( n2.toNewHampshireX() ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testTreeCopy() {
        try {
            final String str_0 = "((((a,b),c),d)[&&NHX:S=lizards],e[&&NHX:S=reptiles])r[&&NHX:S=animals]";
            final Phylogeny t0 = Phylogeny.createInstanceFromNhxString( str_0 );
            final Phylogeny t1 = t0.copy();
            if ( !t1.toNewHampshireX().equals( t0.toNewHampshireX() ) ) {
                return false;
            }
            if ( !t1.toNewHampshireX().equals( str_0 ) ) {
                return false;
            }
            t0.deleteSubtree( t0.getNode( "c" ), true );
            t0.deleteSubtree( t0.getNode( "a" ), true );
            t0.getRoot().getNodeData().getTaxonomy().setScientificName( "metazoa" );
            t0.getNode( "b" ).setName( "Bee" );
            if ( !t0.toNewHampshireX().equals( "((Bee,d)[&&NHX:S=lizards],e[&&NHX:S=reptiles])r[&&NHX:S=metazoa]" ) ) {
                return false;
            }
            if ( !t1.toNewHampshireX().equals( str_0 ) ) {
                return false;
            }
            t0.deleteSubtree( t0.getNode( "e" ), true );
            t0.deleteSubtree( t0.getNode( "Bee" ), true );
            t0.deleteSubtree( t0.getNode( "d" ), true );
            if ( !t1.toNewHampshireX().equals( str_0 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testCreateBalancedPhylogeny() {
        try {
            final Phylogeny p0 = DevelopmentTools.createBalancedPhylogeny( 6, 5 );
            if ( p0.getRoot().getNumberOfDescendants() != 5 ) {
                return false;
            }
            if ( p0.getNumberOfExternalNodes() != 15625 ) {
                return false;
            }
            final Phylogeny p1 = DevelopmentTools.createBalancedPhylogeny( 2, 10 );
            if ( p1.getRoot().getNumberOfDescendants() != 10 ) {
                return false;
            }
            if ( p1.getNumberOfExternalNodes() != 100 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testCreateUriForSeqWeb() {
        try {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "tr|B3RJ64" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.UNIPROT_KB + "B3RJ64" ) ) {
                return false;
            }
            n.setName( "B0LM41_HUMAN" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.UNIPROT_KB + "B0LM41_HUMAN" ) ) {
                return false;
            }
            n.setName( "NP_001025424" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.NCBI_PROTEIN + "NP_001025424" ) ) {
                return false;
            }
            n.setName( "_NM_001030253-" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.NCBI_NUCCORE + "NM_001030253" ) ) {
                return false;
            }
            n.setName( "XM_002122186" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.NCBI_NUCCORE + "XM_002122186" ) ) {
                return false;
            }
            n.setName( "dgh_AAA34956_gdg" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.NCBI_PROTEIN + "AAA34956" ) ) {
                return false;
            }
            n.setName( "AAA34956" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.NCBI_PROTEIN + "AAA34956" ) ) {
                return false;
            }
            n.setName( "GI:394892" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.NCBI_GI + "394892" ) ) {
                System.out.println( TreePanelUtil.createUriForSeqWeb( n, null, null ) );
                return false;
            }
            n.setName( "gi_394892" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.NCBI_GI + "394892" ) ) {
                System.out.println( TreePanelUtil.createUriForSeqWeb( n, null, null ) );
                return false;
            }
            n.setName( "gi6335_gi_394892_56635_Gi_43" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.NCBI_GI + "394892" ) ) {
                System.out.println( TreePanelUtil.createUriForSeqWeb( n, null, null ) );
                return false;
            }
            n.setName( "P12345" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.UNIPROT_KB + "P12345" ) ) {
                System.out.println( TreePanelUtil.createUriForSeqWeb( n, null, null ) );
                return false;
            }
            n.setName( "gi_fdgjmn-3jk5-243 mnefmn fg023-0 P12345 4395jtmnsrg02345m1ggi92450jrg890j4t0j240" );
            if ( !TreePanelUtil.createUriForSeqWeb( n, null, null ).equals( ForesterUtil.UNIPROT_KB + "P12345" ) ) {
                System.out.println( TreePanelUtil.createUriForSeqWeb( n, null, null ) );
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDataObjects() {
        try {
            final Confidence s0 = new Confidence();
            final Confidence s1 = new Confidence();
            if ( !s0.isEqual( s1 ) ) {
                return false;
            }
            final Confidence s2 = new Confidence( 0.23, "bootstrap" );
            final Confidence s3 = new Confidence( 0.23, "bootstrap" );
            if ( s2.isEqual( s1 ) ) {
                return false;
            }
            if ( !s2.isEqual( s3 ) ) {
                return false;
            }
            final Confidence s4 = ( Confidence ) s3.copy();
            if ( !s4.isEqual( s3 ) ) {
                return false;
            }
            s3.asSimpleText();
            s3.asText();
            // Taxonomy
            // ----------
            final Taxonomy t1 = new Taxonomy();
            final Taxonomy t2 = new Taxonomy();
            final Taxonomy t3 = new Taxonomy();
            final Taxonomy t4 = new Taxonomy();
            final Taxonomy t5 = new Taxonomy();
            t1.setIdentifier( new Identifier( "ecoli" ) );
            t1.setTaxonomyCode( "ECOLI" );
            t1.setScientificName( "E. coli" );
            t1.setCommonName( "coli" );
            final Taxonomy t0 = ( Taxonomy ) t1.copy();
            if ( !t1.isEqual( t0 ) ) {
                return false;
            }
            t2.setIdentifier( new Identifier( "ecoli" ) );
            t2.setTaxonomyCode( "OTHER" );
            t2.setScientificName( "what" );
            t2.setCommonName( "something" );
            if ( !t1.isEqual( t2 ) ) {
                return false;
            }
            t2.setIdentifier( new Identifier( "nemve" ) );
            if ( t1.isEqual( t2 ) ) {
                return false;
            }
            t1.setIdentifier( null );
            t3.setTaxonomyCode( "ECOLI" );
            t3.setScientificName( "what" );
            t3.setCommonName( "something" );
            if ( !t1.isEqual( t3 ) ) {
                return false;
            }
            t1.setIdentifier( null );
            t1.setTaxonomyCode( "" );
            t4.setScientificName( "E. ColI" );
            t4.setCommonName( "something" );
            if ( !t1.isEqual( t4 ) ) {
                return false;
            }
            t4.setScientificName( "B. subtilis" );
            t4.setCommonName( "something" );
            if ( t1.isEqual( t4 ) ) {
                return false;
            }
            t1.setIdentifier( null );
            t1.setTaxonomyCode( "" );
            t1.setScientificName( "" );
            t5.setCommonName( "COLI" );
            if ( !t1.isEqual( t5 ) ) {
                return false;
            }
            t5.setCommonName( "vibrio" );
            if ( t1.isEqual( t5 ) ) {
                return false;
            }
            // Identifier
            // ----------
            final Identifier id0 = new Identifier( "123", "pfam" );
            final Identifier id1 = ( Identifier ) id0.copy();
            if ( !id1.isEqual( id1 ) ) {
                return false;
            }
            if ( !id1.isEqual( id0 ) ) {
                return false;
            }
            if ( !id0.isEqual( id1 ) ) {
                return false;
            }
            id1.asSimpleText();
            id1.asText();
            // ProteinDomain
            // ---------------
            final ProteinDomain pd0 = new ProteinDomain( "abc", 100, 200 );
            final ProteinDomain pd1 = ( ProteinDomain ) pd0.copy();
            if ( !pd1.isEqual( pd1 ) ) {
                return false;
            }
            if ( !pd1.isEqual( pd0 ) ) {
                return false;
            }
            pd1.asSimpleText();
            pd1.asText();
            final ProteinDomain pd2 = new ProteinDomain( pd0.getName(), pd0.getFrom(), pd0.getTo(), "id" );
            final ProteinDomain pd3 = ( ProteinDomain ) pd2.copy();
            if ( !pd3.isEqual( pd3 ) ) {
                return false;
            }
            if ( !pd2.isEqual( pd3 ) ) {
                return false;
            }
            if ( !pd0.isEqual( pd3 ) ) {
                return false;
            }
            pd3.asSimpleText();
            pd3.asText();
            // DomainArchitecture
            // ------------------
            final ProteinDomain d0 = new ProteinDomain( "domain0", 10, 20 );
            final ProteinDomain d1 = new ProteinDomain( "domain1", 30, 40 );
            final ProteinDomain d2 = new ProteinDomain( "domain2", 50, 60 );
            final ProteinDomain d3 = new ProteinDomain( "domain3", 70, 80 );
            final ProteinDomain d4 = new ProteinDomain( "domain4", 90, 100 );
            final ArrayList<PhylogenyData> domains0 = new ArrayList<PhylogenyData>();
            domains0.add( d2 );
            domains0.add( d0 );
            domains0.add( d3 );
            domains0.add( d1 );
            final DomainArchitecture ds0 = new DomainArchitecture( domains0, 110 );
            if ( ds0.getNumberOfDomains() != 4 ) {
                return false;
            }
            final DomainArchitecture ds1 = ( DomainArchitecture ) ds0.copy();
            if ( !ds0.isEqual( ds0 ) ) {
                return false;
            }
            if ( !ds0.isEqual( ds1 ) ) {
                return false;
            }
            if ( ds1.getNumberOfDomains() != 4 ) {
                return false;
            }
            final ArrayList<PhylogenyData> domains1 = new ArrayList<PhylogenyData>();
            domains1.add( d1 );
            domains1.add( d2 );
            domains1.add( d4 );
            domains1.add( d0 );
            final DomainArchitecture ds2 = new DomainArchitecture( domains1, 200 );
            if ( ds0.isEqual( ds2 ) ) {
                return false;
            }
            ds1.asSimpleText();
            ds1.asText();
            ds1.toNHX();
            final DomainArchitecture ds3 = new DomainArchitecture( "120>30>40>0.9>b>50>60>0.4>c>10>20>0.1>a" );
            if ( !ds3.toNHX().toString().equals( ":DS=120>10>20>0.1>a>30>40>0.9>b>50>60>0.4>c" ) ) {
                System.out.println( ds3.toNHX() );
                return false;
            }
            if ( ds3.getNumberOfDomains() != 3 ) {
                return false;
            }
            // Event
            // -----
            final Event e1 = new Event( Event.EventType.fusion );
            if ( e1.isDuplication() ) {
                return false;
            }
            if ( !e1.isFusion() ) {
                return false;
            }
            if ( !e1.asText().toString().equals( "fusion" ) ) {
                return false;
            }
            if ( !e1.asSimpleText().toString().equals( "fusion" ) ) {
                return false;
            }
            final Event e11 = new Event( Event.EventType.fusion );
            if ( !e11.isEqual( e1 ) ) {
                return false;
            }
            if ( !e11.toNHX().toString().equals( "" ) ) {
                return false;
            }
            final Event e2 = new Event( Event.EventType.speciation_or_duplication );
            if ( e2.isDuplication() ) {
                return false;
            }
            if ( !e2.isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !e2.asText().toString().equals( "speciation_or_duplication" ) ) {
                return false;
            }
            if ( !e2.asSimpleText().toString().equals( "?" ) ) {
                return false;
            }
            if ( !e2.toNHX().toString().equals( ":D=?" ) ) {
                return false;
            }
            if ( e11.isEqual( e2 ) ) {
                return false;
            }
            final Event e2c = ( Event ) e2.copy();
            if ( !e2c.isEqual( e2 ) ) {
                return false;
            }
            Event e3 = new Event( 1, 2, 3 );
            if ( e3.isDuplication() ) {
                return false;
            }
            if ( e3.isSpeciation() ) {
                return false;
            }
            if ( e3.isGeneLoss() ) {
                return false;
            }
            if ( !e3.asText().toString().equals( "duplications [1] speciations [2] gene-losses [3]" ) ) {
                return false;
            }
            final Event e3c = ( Event ) e3.copy();
            final Event e3cc = ( Event ) e3c.copy();
            if ( !e3c.asSimpleText().toString().equals( "D2S3L" ) ) {
                return false;
            }
            e3 = null;
            if ( !e3c.isEqual( e3cc ) ) {
                return false;
            }
            Event e4 = new Event( 1, 2, 3 );
            if ( !e4.asText().toString().equals( "duplications [1] speciations [2] gene-losses [3]" ) ) {
                return false;
            }
            if ( !e4.asSimpleText().toString().equals( "D2S3L" ) ) {
                return false;
            }
            final Event e4c = ( Event ) e4.copy();
            e4 = null;
            final Event e4cc = ( Event ) e4c.copy();
            if ( !e4cc.asText().toString().equals( "duplications [1] speciations [2] gene-losses [3]" ) ) {
                return false;
            }
            if ( !e4c.isEqual( e4cc ) ) {
                return false;
            }
            final Event e5 = new Event();
            if ( !e5.isUnassigned() ) {
                return false;
            }
            if ( !e5.asText().toString().equals( "unassigned" ) ) {
                return false;
            }
            if ( !e5.asSimpleText().toString().equals( "" ) ) {
                return false;
            }
            final Event e6 = new Event( 1, 0, 0 );
            if ( !e6.asText().toString().equals( "duplication" ) ) {
                return false;
            }
            if ( !e6.asSimpleText().toString().equals( "D" ) ) {
                return false;
            }
            final Event e7 = new Event( 0, 1, 0 );
            if ( !e7.asText().toString().equals( "speciation" ) ) {
                return false;
            }
            if ( !e7.asSimpleText().toString().equals( "S" ) ) {
                return false;
            }
            final Event e8 = new Event( 0, 0, 1 );
            if ( !e8.asText().toString().equals( "gene-loss" ) ) {
                return false;
            }
            if ( !e8.asSimpleText().toString().equals( "L" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDeletionOfExternalNodes() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "A", new NHXParser() )[ 0 ];
            final PhylogenyWriter w = new PhylogenyWriter();
            if ( t0.isEmpty() ) {
                return false;
            }
            if ( t0.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t0.deleteSubtree( t0.getNode( "A" ), false );
            if ( t0.getNumberOfExternalNodes() != 0 ) {
                return false;
            }
            if ( !t0.isEmpty() ) {
                return false;
            }
            final Phylogeny t1 = factory.create( "(A,B)r", new NHXParser() )[ 0 ];
            if ( t1.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "A" ), false );
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( !t1.getNode( "B" ).getName().equals( "B" ) ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "B" ), false );
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "r" ), false );
            if ( !t1.isEmpty() ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "((A,B),C)", new NHXParser() )[ 0 ];
            if ( t2.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "B" ), false );
            if ( t2.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            t2.toNewHampshireX();
            PhylogenyNode n = t2.getNode( "A" );
            if ( !n.getNextExternalNode().getName().equals( "C" ) ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "A" ), false );
            if ( t2.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "C" ), true );
            if ( t2.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            final Phylogeny t3 = factory.create( "((A,B),(C,D))", new NHXParser() )[ 0 ];
            if ( t3.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            t3.deleteSubtree( t3.getNode( "B" ), true );
            if ( t3.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            n = t3.getNode( "A" );
            if ( !n.getNextExternalNode().getName().equals( "C" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNextExternalNode().getName().equals( "D" ) ) {
                return false;
            }
            t3.deleteSubtree( t3.getNode( "A" ), true );
            if ( t3.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            n = t3.getNode( "C" );
            if ( !n.getNextExternalNode().getName().equals( "D" ) ) {
                return false;
            }
            t3.deleteSubtree( t3.getNode( "C" ), true );
            if ( t3.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t3.deleteSubtree( t3.getNode( "D" ), true );
            if ( t3.getNumberOfExternalNodes() != 0 ) {
                return false;
            }
            final Phylogeny t4 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            if ( t4.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            t4.deleteSubtree( t4.getNode( "B2" ), true );
            if ( t4.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            String s = w.toNewHampshire( t4, false, true ).toString();
            if ( !s.equals( "((A,(B11,B12)),(C,D));" ) ) {
                return false;
            }
            t4.deleteSubtree( t4.getNode( "B11" ), true );
            if ( t4.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            t4.deleteSubtree( t4.getNode( "C" ), true );
            if ( t4.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            n = t4.getNode( "A" );
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "B12" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "D" ) ) {
                return false;
            }
            s = w.toNewHampshire( t4, false, true ).toString();
            if ( !s.equals( "((A,B12),D);" ) ) {
                return false;
            }
            final Phylogeny t5 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t5.deleteSubtree( t5.getNode( "A" ), true );
            if ( t5.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t5, false, true ).toString();
            if ( !s.equals( "(((B11,B12),B2),(C,D));" ) ) {
                return false;
            }
            final Phylogeny t6 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t6.deleteSubtree( t6.getNode( "B11" ), true );
            if ( t6.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t6, false, false ).toString();
            if ( !s.equals( "((A,(B12,B2)),(C,D));" ) ) {
                return false;
            }
            final Phylogeny t7 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t7.deleteSubtree( t7.getNode( "B12" ), true );
            if ( t7.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t7, false, true ).toString();
            if ( !s.equals( "((A,(B11,B2)),(C,D));" ) ) {
                return false;
            }
            final Phylogeny t8 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t8.deleteSubtree( t8.getNode( "B2" ), true );
            if ( t8.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t8, false, false ).toString();
            if ( !s.equals( "((A,(B11,B12)),(C,D));" ) ) {
                return false;
            }
            final Phylogeny t9 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t9.deleteSubtree( t9.getNode( "C" ), true );
            if ( t9.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t9, false, true ).toString();
            if ( !s.equals( "((A,((B11,B12),B2)),D);" ) ) {
                return false;
            }
            final Phylogeny t10 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t10.deleteSubtree( t10.getNode( "D" ), true );
            if ( t10.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t10, false, true ).toString();
            if ( !s.equals( "((A,((B11,B12),B2)),C);" ) ) {
                return false;
            }
            final Phylogeny t11 = factory.create( "(A,B,C)", new NHXParser() )[ 0 ];
            t11.deleteSubtree( t11.getNode( "A" ), true );
            if ( t11.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            s = w.toNewHampshire( t11, false, true ).toString();
            if ( !s.equals( "(B,C);" ) ) {
                return false;
            }
            t11.deleteSubtree( t11.getNode( "C" ), true );
            if ( t11.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            s = w.toNewHampshire( t11, false, false ).toString();
            if ( !s.equals( "B;" ) ) {
                return false;
            }
            final Phylogeny t12 = factory.create( "((A1,A2,A3),(B1,B2,B3),(C1,C2,C3))", new NHXParser() )[ 0 ];
            t12.deleteSubtree( t12.getNode( "B2" ), true );
            if ( t12.getNumberOfExternalNodes() != 8 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A1,A2,A3),(B1,B3),(C1,C2,C3));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "B3" ), true );
            if ( t12.getNumberOfExternalNodes() != 7 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A1,A2,A3),B1,(C1,C2,C3));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "C3" ), true );
            if ( t12.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A1,A2,A3),B1,(C1,C2));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "A1" ), true );
            if ( t12.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A2,A3),B1,(C1,C2));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "B1" ), true );
            if ( t12.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A2,A3),(C1,C2));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "A3" ), true );
            if ( t12.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "(A2,(C1,C2));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "A2" ), true );
            if ( t12.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "(C1,C2);" ) ) {
                return false;
            }
            final Phylogeny t13 = factory.create( "(A,B,C,(D:1.0,E:2.0):3.0)", new NHXParser() )[ 0 ];
            t13.deleteSubtree( t13.getNode( "D" ), true );
            if ( t13.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            s = w.toNewHampshire( t13, false, true ).toString();
            if ( !s.equals( "(A,B,C,E:5.0);" ) ) {
                return false;
            }
            final Phylogeny t14 = factory.create( "((A,B,C,(D:0.1,E:0.4):1.0),F)", new NHXParser() )[ 0 ];
            t14.deleteSubtree( t14.getNode( "E" ), true );
            if ( t14.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t14, false, true ).toString();
            if ( !s.equals( "((A,B,C,D:1.1),F);" ) ) {
                return false;
            }
            final Phylogeny t15 = factory.create( "((A1,A2,A3,A4),(B1,B2,B3,B4),(C1,C2,C3,C4))", new NHXParser() )[ 0 ];
            t15.deleteSubtree( t15.getNode( "B2" ), true );
            if ( t15.getNumberOfExternalNodes() != 11 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "B1" ), true );
            if ( t15.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "B3" ), true );
            if ( t15.getNumberOfExternalNodes() != 9 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "B4" ), true );
            if ( t15.getNumberOfExternalNodes() != 8 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "A1" ), true );
            if ( t15.getNumberOfExternalNodes() != 7 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "C4" ), true );
            if ( t15.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDescriptiveStatistics() {
        try {
            final DescriptiveStatistics dss1 = new BasicDescriptiveStatistics();
            dss1.addValue( 82 );
            dss1.addValue( 78 );
            dss1.addValue( 70 );
            dss1.addValue( 58 );
            dss1.addValue( 42 );
            if ( dss1.getN() != 5 ) {
                return false;
            }
            if ( !Test.isEqual( dss1.getMin(), 42 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.getMax(), 82 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.arithmeticMean(), 66 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.sampleStandardDeviation(), 16.24807680927192 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.median(), 70 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.midrange(), 62 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.sampleVariance(), 264 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.pearsonianSkewness(), -0.7385489458759964 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.coefficientOfVariation(), 0.24618298195866547 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.sampleStandardUnit( 66 - 16.24807680927192 ), -1.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.getValue( 1 ), 78 ) ) {
                return false;
            }
            dss1.addValue( 123 );
            if ( !Test.isEqual( dss1.arithmeticMean(), 75.5 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.getMax(), 123 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.standardErrorOfMean(), 11.200446419674531 ) ) {
                return false;
            }
            final DescriptiveStatistics dss2 = new BasicDescriptiveStatistics();
            dss2.addValue( -1.85 );
            dss2.addValue( 57.5 );
            dss2.addValue( 92.78 );
            dss2.addValue( 57.78 );
            if ( !Test.isEqual( dss2.median(), 57.64 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss2.sampleStandardDeviation(), 39.266984753946495 ) ) {
                return false;
            }
            final double[] a = dss2.getDataAsDoubleArray();
            if ( !Test.isEqual( a[ 3 ], 57.78 ) ) {
                return false;
            }
            dss2.addValue( -100 );
            if ( !Test.isEqual( dss2.sampleStandardDeviation(), 75.829111296388 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss2.sampleVariance(), 5750.05412 ) ) {
                return false;
            }
            final double[] ds = new double[ 14 ];
            ds[ 0 ] = 34;
            ds[ 1 ] = 23;
            ds[ 2 ] = 1;
            ds[ 3 ] = 32;
            ds[ 4 ] = 11;
            ds[ 5 ] = 2;
            ds[ 6 ] = 12;
            ds[ 7 ] = 33;
            ds[ 8 ] = 13;
            ds[ 9 ] = 22;
            ds[ 10 ] = 21;
            ds[ 11 ] = 35;
            ds[ 12 ] = 24;
            ds[ 13 ] = 31;
            final int[] bins = BasicDescriptiveStatistics.performBinning( ds, 0, 40, 4 );
            if ( bins.length != 4 ) {
                return false;
            }
            if ( bins[ 0 ] != 2 ) {
                return false;
            }
            if ( bins[ 1 ] != 3 ) {
                return false;
            }
            if ( bins[ 2 ] != 4 ) {
                return false;
            }
            if ( bins[ 3 ] != 5 ) {
                return false;
            }
            final double[] ds1 = new double[ 9 ];
            ds1[ 0 ] = 10.0;
            ds1[ 1 ] = 19.0;
            ds1[ 2 ] = 9.999;
            ds1[ 3 ] = 0.0;
            ds1[ 4 ] = 39.9;
            ds1[ 5 ] = 39.999;
            ds1[ 6 ] = 30.0;
            ds1[ 7 ] = 19.999;
            ds1[ 8 ] = 30.1;
            final int[] bins1 = BasicDescriptiveStatistics.performBinning( ds1, 0, 40, 4 );
            if ( bins1.length != 4 ) {
                return false;
            }
            if ( bins1[ 0 ] != 2 ) {
                return false;
            }
            if ( bins1[ 1 ] != 3 ) {
                return false;
            }
            if ( bins1[ 2 ] != 0 ) {
                return false;
            }
            if ( bins1[ 3 ] != 4 ) {
                return false;
            }
            final int[] bins1_1 = BasicDescriptiveStatistics.performBinning( ds1, 0, 40, 3 );
            if ( bins1_1.length != 3 ) {
                return false;
            }
            if ( bins1_1[ 0 ] != 3 ) {
                return false;
            }
            if ( bins1_1[ 1 ] != 2 ) {
                return false;
            }
            if ( bins1_1[ 2 ] != 4 ) {
                return false;
            }
            final int[] bins1_2 = BasicDescriptiveStatistics.performBinning( ds1, 1, 39, 3 );
            if ( bins1_2.length != 3 ) {
                return false;
            }
            if ( bins1_2[ 0 ] != 2 ) {
                return false;
            }
            if ( bins1_2[ 1 ] != 2 ) {
                return false;
            }
            if ( bins1_2[ 2 ] != 2 ) {
                return false;
            }
            final DescriptiveStatistics dss3 = new BasicDescriptiveStatistics();
            dss3.addValue( 1 );
            dss3.addValue( 1 );
            dss3.addValue( 1 );
            dss3.addValue( 2 );
            dss3.addValue( 3 );
            dss3.addValue( 4 );
            dss3.addValue( 5 );
            dss3.addValue( 5 );
            dss3.addValue( 5 );
            dss3.addValue( 6 );
            dss3.addValue( 7 );
            dss3.addValue( 8 );
            dss3.addValue( 9 );
            dss3.addValue( 10 );
            dss3.addValue( 10 );
            dss3.addValue( 10 );
            final AsciiHistogram histo = new AsciiHistogram( dss3 );
            histo.toStringBuffer( 10, '=', 40, 5 );
            histo.toStringBuffer( 3, 8, 10, '=', 40, 5, null );
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDir( final String file ) {
        try {
            final File f = new File( file );
            if ( !f.exists() ) {
                return false;
            }
            if ( !f.isDirectory() ) {
                return false;
            }
            if ( !f.canRead() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            return false;
        }
        return true;
    }

    private static boolean testGenbankAccessorParsing() {
        //The format for GenBank Accession numbers are:
        //Nucleotide: 1 letter + 5 numerals OR 2 letters + 6 numerals
        //Protein:    3 letters + 5 numerals
        //http://www.ncbi.nlm.nih.gov/Sequin/acc.html
        if ( !SequenceAccessionTools.parseGenbankAccessorFromString( "AY423861" ).equals( "AY423861" ) ) {
            return false;
        }
        if ( !SequenceAccessionTools.parseGenbankAccessorFromString( ".AY423861.2" ).equals( "AY423861.2" ) ) {
            return false;
        }
        if ( !SequenceAccessionTools.parseGenbankAccessorFromString( "345_.AY423861.24_345" ).equals( "AY423861.24" ) ) {
            return false;
        }
        if ( SequenceAccessionTools.parseGenbankAccessorFromString( "AAY423861" ) != null ) {
            return false;
        }
        if ( SequenceAccessionTools.parseGenbankAccessorFromString( "AY4238612" ) != null ) {
            return false;
        }
        if ( SequenceAccessionTools.parseGenbankAccessorFromString( "AAY4238612" ) != null ) {
            return false;
        }
        if ( SequenceAccessionTools.parseGenbankAccessorFromString( "Y423861" ) != null ) {
            return false;
        }
        if ( !SequenceAccessionTools.parseGenbankAccessorFromString( "S12345" ).equals( "S12345" ) ) {
            return false;
        }
        if ( !SequenceAccessionTools.parseGenbankAccessorFromString( "|S12345|" ).equals( "S12345" ) ) {
            return false;
        }
        if ( SequenceAccessionTools.parseGenbankAccessorFromString( "|S123456" ) != null ) {
            return false;
        }
        if ( SequenceAccessionTools.parseGenbankAccessorFromString( "ABC123456" ) != null ) {
            return false;
        }
        if ( !SequenceAccessionTools.parseGenbankAccessorFromString( "ABC12345" ).equals( "ABC12345" ) ) {
            return false;
        }
        if ( !SequenceAccessionTools.parseGenbankAccessorFromString( "&ABC12345&" ).equals( "ABC12345" ) ) {
            return false;
        }
        if ( SequenceAccessionTools.parseGenbankAccessorFromString( "ABCD12345" ) != null ) {
            return false;
        }
        return true;
    }

    private static boolean testExternalNodeRelatedMethods() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t1 = factory.create( "((A,B),(C,D))", new NHXParser() )[ 0 ];
            PhylogenyNode n = t1.getNode( "A" );
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "B" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "C" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "D" ) ) {
                return false;
            }
            n = t1.getNode( "B" );
            while ( !n.isLastExternalNode() ) {
                n = n.getNextExternalNode();
            }
            final Phylogeny t2 = factory.create( "(((A,B),C),D)", new NHXParser() )[ 0 ];
            n = t2.getNode( "A" );
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "B" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "C" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "D" ) ) {
                return false;
            }
            n = t2.getNode( "B" );
            while ( !n.isLastExternalNode() ) {
                n = n.getNextExternalNode();
            }
            final Phylogeny t3 = factory.create( "(((A,B),(C,D)),((E,F),(G,H)))", new NHXParser() )[ 0 ];
            n = t3.getNode( "A" );
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "B" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "C" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "D" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "E" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "F" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "G" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getName().equals( "H" ) ) {
                return false;
            }
            n = t3.getNode( "B" );
            while ( !n.isLastExternalNode() ) {
                n = n.getNextExternalNode();
            }
            final Phylogeny t4 = factory.create( "((A,B),(C,D))", new NHXParser() )[ 0 ];
            for( final PhylogenyNodeIterator iter = t4.iteratorExternalForward(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
            }
            final Phylogeny t5 = factory.create( "(((A,B),(C,D)),((E,F),(G,H)))", new NHXParser() )[ 0 ];
            for( final PhylogenyNodeIterator iter = t5.iteratorExternalForward(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
            }
            final Phylogeny t6 = factory.create( "((((((A))),(((B))),((C)),((((D)))),E)),((F)))", new NHXParser() )[ 0 ];
            final PhylogenyNodeIterator iter = t6.iteratorExternalForward();
            if ( !iter.next().getName().equals( "A" ) ) {
                return false;
            }
            if ( !iter.next().getName().equals( "B" ) ) {
                return false;
            }
            if ( !iter.next().getName().equals( "C" ) ) {
                return false;
            }
            if ( !iter.next().getName().equals( "D" ) ) {
                return false;
            }
            if ( !iter.next().getName().equals( "E" ) ) {
                return false;
            }
            if ( !iter.next().getName().equals( "F" ) ) {
                return false;
            }
            if ( iter.hasNext() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testExtractSNFromNodeName() {
        try {
            if ( !ParserUtils.extractScientificNameFromNodeName( "BCDO2_Mus_musculus" ).equals( "Mus musculus" ) ) {
                return false;
            }
            if ( !ParserUtils.extractScientificNameFromNodeName( "BCDO2_Mus_musculus_musculus" )
                    .equals( "Mus musculus musculus" ) ) {
                return false;
            }
            if ( !ParserUtils.extractScientificNameFromNodeName( "BCDO2_Mus_musculus_musculus-12" )
                    .equals( "Mus musculus musculus" ) ) {
                return false;
            }
            if ( !ParserUtils.extractScientificNameFromNodeName( " -XS12_Mus_musculus-12" ).equals( "Mus musculus" ) ) {
                return false;
            }
            if ( !ParserUtils.extractScientificNameFromNodeName( " -1234_Mus_musculus-12 affrre e" )
                    .equals( "Mus musculus" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testExtractTaxonomyCodeFromNodeName() {
        try {
            if ( ParserUtils.extractTaxonomyCodeFromNodeName( "MOUSE", TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ) != null ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "SOYBN", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( " ARATH ", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "ARATH" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( " ARATH ", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "ARATH" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "RAT", TAXONOMY_EXTRACTION.AGGRESSIVE ).equals( "RAT" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "RAT", TAXONOMY_EXTRACTION.AGGRESSIVE ).equals( "RAT" ) ) {
                return false;
            }
            if ( ParserUtils.extractTaxonomyCodeFromNodeName( "RAT1", TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ) != null ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( " _SOYBN", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "SOYBN", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "qwerty SOYBN", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "qwerty_SOYBN", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "ABCD_SOYBN ", TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "SOYBN", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( ",SOYBN,", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "xxx,SOYBN,xxx", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( ParserUtils.extractTaxonomyCodeFromNodeName( "xxxSOYBNxxx", TAXONOMY_EXTRACTION.AGGRESSIVE ) != null ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "-SOYBN~", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "SOYBN" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "NNN8_ECOLI/1-2:0.01",
                                                               TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT ).equals( "ECOLI" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "blag_9YX45-blag", TAXONOMY_EXTRACTION.AGGRESSIVE )
                    .equals( "9YX45" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_MOUSE function = 23445",
                                                               TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED )
                    .equals( "MOUSE" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_MOUSE+function = 23445",
                                                               TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED )
                    .equals( "MOUSE" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_MOUSE|function = 23445",
                                                               TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED )
                    .equals( "MOUSE" ) ) {
                return false;
            }
            if ( ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_MOUSEfunction = 23445",
                                                              TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ) != null ) {
                return false;
            }
            if ( ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_MOUSEFunction = 23445",
                                                              TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ) != null ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_RAT function = 23445",
                                                               TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ).equals( "RAT" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_RAT function = 23445",
                                                               TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ).equals( "RAT" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_RAT|function = 23445",
                                                               TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ).equals( "RAT" ) ) {
                return false;
            }
            if ( ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_RATfunction = 23445",
                                                              TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ) != null ) {
                return false;
            }
            if ( ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_RATFunction = 23445",
                                                              TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ) != null ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_RAT/1-3", TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED )
                    .equals( "RAT" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_PIG/1-3", TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT )
                    .equals( "PIG" ) ) {
                return false;
            }
            if ( !ParserUtils
                    .extractTaxonomyCodeFromNodeName( "BCL2_MOUSE/1-3", TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED )
                    .equals( "MOUSE" ) ) {
                return false;
            }
            if ( !ParserUtils.extractTaxonomyCodeFromNodeName( "BCL2_MOUSE/1-3", TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT )
                    .equals( "MOUSE" ) ) {
                return false;
            }
            if ( ParserUtils.extractTaxonomyCodeFromNodeName( "_MOUSE ", TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED ) != null ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testExtractUniProtKbProteinSeqIdentifier() {
        try {
            PhylogenyNode n = new PhylogenyNode();
            n.setName( "tr|B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "tr.B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "tr=B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "tr-B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "tr/B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "tr\\B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "tr_B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( " tr|B3RJ64 " );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "-tr|B3RJ64-" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "-tr=B3RJ64-" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "_tr=B3RJ64_" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( " tr_tr|B3RJ64_sp|123 " );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "sp|B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "sp|B3RJ64C" );
            if ( SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ) != null ) {
                return false;
            }
            n.setName( "sp B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n.setName( "sp|B3RJ6X" );
            if ( SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ) != null ) {
                return false;
            }
            n.setName( "sp|B3RJ6" );
            if ( SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ) != null ) {
                return false;
            }
            n.setName( "K1PYK7_CRAGI" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "K1PYK7_CRAGI" ) ) {
                return false;
            }
            n.setName( "K1PYK7_PEA" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "K1PYK7_PEA" ) ) {
                return false;
            }
            n.setName( "K1PYK7_RAT" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "K1PYK7_RAT" ) ) {
                return false;
            }
            n.setName( "K1PYK7_PIG" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "K1PYK7_PIG" ) ) {
                return false;
            }
            n.setName( "~K1PYK7_PIG~" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "K1PYK7_PIG" ) ) {
                return false;
            }
            n.setName( "123456_ECOLI-K1PYK7_CRAGI-sp" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "K1PYK7_CRAGI" ) ) {
                return false;
            }
            n.setName( "K1PYKX_CRAGI" );
            if ( SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ) != null ) {
                return false;
            }
            n.setName( "XXXXX_CRAGI" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "XXXXX_CRAGI" ) ) {
                return false;
            }
            n.setName( "tr|H3IB65|H3IB65_STRPU~2-2" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "H3IB65" ) ) {
                return false;
            }
            n.setName( "jgi|Lacbi2|181470|Lacbi1.estExt_GeneWisePlus_human.C_10729~2-3" );
            if ( SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ) != null ) {
                return false;
            }
            n.setName( "sp|Q86U06|RBM23_HUMAN~2-2" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "Q86U06" ) ) {
                return false;
            }
            n = new PhylogenyNode();
            org.forester.phylogeny.data.Sequence seq = new org.forester.phylogeny.data.Sequence();
            seq.setSymbol( "K1PYK7_CRAGI" );
            n.getNodeData().addSequence( seq );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "K1PYK7_CRAGI" ) ) {
                return false;
            }
            seq.setSymbol( "tr|B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n = new PhylogenyNode();
            seq = new org.forester.phylogeny.data.Sequence();
            seq.setName( "K1PYK7_CRAGI" );
            n.getNodeData().addSequence( seq );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "K1PYK7_CRAGI" ) ) {
                return false;
            }
            seq.setName( "tr|B3RJ64" );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            n = new PhylogenyNode();
            seq = new org.forester.phylogeny.data.Sequence();
            seq.setAccession( new Accession( "K1PYK8_CRAGI", "?" ) );
            n.getNodeData().addSequence( seq );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "K1PYK8_CRAGI" ) ) {
                return false;
            }
            n = new PhylogenyNode();
            seq = new org.forester.phylogeny.data.Sequence();
            seq.setAccession( new Accession( "tr|B3RJ64", "?" ) );
            n.getNodeData().addSequence( seq );
            if ( !SequenceAccessionTools.obtainUniProtAccessorFromDataFields( n ).equals( "B3RJ64" ) ) {
                return false;
            }
            //
            n = new PhylogenyNode();
            n.setName( "ACP19736" );
            if ( !SequenceAccessionTools.obtainGenbankAccessorFromDataFields( n ).equals( "ACP19736" ) ) {
                return false;
            }
            n = new PhylogenyNode();
            n.setName( "|ACP19736|" );
            if ( !SequenceAccessionTools.obtainGenbankAccessorFromDataFields( n ).equals( "ACP19736" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testFastaParser() {
        try {
            if ( !FastaParser.isLikelyFasta( new FileInputStream( PATH_TO_TEST_DATA + "fasta_0.fasta" ) ) ) {
                return false;
            }
            if ( FastaParser.isLikelyFasta( new FileInputStream( PATH_TO_TEST_DATA + "msa_3.txt" ) ) ) {
                return false;
            }
            final Msa msa_0 = FastaParser.parseMsa( new FileInputStream( PATH_TO_TEST_DATA + "fasta_0.fasta" ) );
            if ( !msa_0.getSequenceAsString( 0 ).toString().equalsIgnoreCase( "ACGTGKXFMFDMXEXXXSFMFMF" ) ) {
                return false;
            }
            if ( !msa_0.getIdentifier( 0 ).equals( "one dumb" ) ) {
                return false;
            }
            if ( !msa_0.getSequenceAsString( 1 ).toString().equalsIgnoreCase( "DKXASDFXSFXFKFKSXDFKSLX" ) ) {
                return false;
            }
            if ( !msa_0.getSequenceAsString( 2 ).toString().equalsIgnoreCase( "SXDFKSXLFSFPWEXPRXWXERR" ) ) {
                return false;
            }
            if ( !msa_0.getSequenceAsString( 3 ).toString().equalsIgnoreCase( "AAAAAAAAAAAAAAAAAAAAAAA" ) ) {
                return false;
            }
            if ( !msa_0.getSequenceAsString( 4 ).toString().equalsIgnoreCase( "DDDDDDDDDDDDDDDDDDDDAXF" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testGeneralMsaParser() {
        try {
            final String msa_str_0 = "seq1 abcd\n\nseq2 efgh\n";
            final Msa msa_0 = GeneralMsaParser.parse( new ByteArrayInputStream( msa_str_0.getBytes() ) );
            final String msa_str_1 = "seq1 abc\nseq2 ghi\nseq1 def\nseq2 jkm\n";
            final Msa msa_1 = GeneralMsaParser.parse( new ByteArrayInputStream( msa_str_1.getBytes() ) );
            final String msa_str_2 = "seq1 abc\nseq2 ghi\n\ndef\njkm\n";
            final Msa msa_2 = GeneralMsaParser.parse( new ByteArrayInputStream( msa_str_2.getBytes() ) );
            final String msa_str_3 = "seq1 abc\n def\nseq2 ghi\n jkm\n";
            final Msa msa_3 = GeneralMsaParser.parse( new ByteArrayInputStream( msa_str_3.getBytes() ) );
            if ( !msa_1.getSequenceAsString( 0 ).toString().equalsIgnoreCase( "abcdef" ) ) {
                return false;
            }
            if ( !msa_1.getSequenceAsString( 1 ).toString().equalsIgnoreCase( "ghixkm" ) ) {
                return false;
            }
            if ( !msa_1.getIdentifier( 0 ).toString().equals( "seq1" ) ) {
                return false;
            }
            if ( !msa_1.getIdentifier( 1 ).toString().equals( "seq2" ) ) {
                return false;
            }
            if ( !msa_2.getSequenceAsString( 0 ).toString().equalsIgnoreCase( "abcdef" ) ) {
                return false;
            }
            if ( !msa_2.getSequenceAsString( 1 ).toString().equalsIgnoreCase( "ghixkm" ) ) {
                return false;
            }
            if ( !msa_2.getIdentifier( 0 ).toString().equals( "seq1" ) ) {
                return false;
            }
            if ( !msa_2.getIdentifier( 1 ).toString().equals( "seq2" ) ) {
                return false;
            }
            if ( !msa_3.getSequenceAsString( 0 ).toString().equalsIgnoreCase( "abcdef" ) ) {
                return false;
            }
            if ( !msa_3.getSequenceAsString( 1 ).toString().equalsIgnoreCase( "ghixkm" ) ) {
                return false;
            }
            if ( !msa_3.getIdentifier( 0 ).toString().equals( "seq1" ) ) {
                return false;
            }
            if ( !msa_3.getIdentifier( 1 ).toString().equals( "seq2" ) ) {
                return false;
            }
            final Msa msa_4 = GeneralMsaParser.parse( new FileInputStream( PATH_TO_TEST_DATA + "msa_1.txt" ) );
            if ( !msa_4.getSequenceAsString( 0 ).toString().equalsIgnoreCase( "abcdefeeeeeeeexx" ) ) {
                return false;
            }
            if ( !msa_4.getSequenceAsString( 1 ).toString().equalsIgnoreCase( "efghixffffffffyy" ) ) {
                return false;
            }
            if ( !msa_4.getSequenceAsString( 2 ).toString().equalsIgnoreCase( "klmnxphhhhhhhhzz" ) ) {
                return false;
            }
            final Msa msa_5 = GeneralMsaParser.parse( new FileInputStream( PATH_TO_TEST_DATA + "msa_2.txt" ) );
            if ( !msa_5.getSequenceAsString( 0 ).toString().equalsIgnoreCase( "abcdefxx" ) ) {
                return false;
            }
            if ( !msa_5.getSequenceAsString( 1 ).toString().equalsIgnoreCase( "efghixyy" ) ) {
                return false;
            }
            if ( !msa_5.getSequenceAsString( 2 ).toString().equalsIgnoreCase( "klmnxpzz" ) ) {
                return false;
            }
            final Msa msa_6 = GeneralMsaParser.parse( new FileInputStream( PATH_TO_TEST_DATA + "msa_3.txt" ) );
            if ( !msa_6.getSequenceAsString( 0 ).toString().equalsIgnoreCase( "abcdefeeeeeeeexx" ) ) {
                return false;
            }
            if ( !msa_6.getSequenceAsString( 1 ).toString().equalsIgnoreCase( "efghixffffffffyy" ) ) {
                return false;
            }
            if ( !msa_6.getSequenceAsString( 2 ).toString().equalsIgnoreCase( "klmnxphhhhhhhhzz" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testGeneralTable() {
        try {
            final GeneralTable<Integer, String> t0 = new GeneralTable<Integer, String>();
            t0.setValue( 3, 2, "23" );
            t0.setValue( 10, 1, "error" );
            t0.setValue( 10, 1, "110" );
            t0.setValue( 9, 1, "19" );
            t0.setValue( 1, 10, "101" );
            t0.setValue( 10, 10, "1010" );
            t0.setValue( 100, 10, "10100" );
            t0.setValue( 0, 0, "00" );
            if ( !t0.getValue( 3, 2 ).equals( "23" ) ) {
                return false;
            }
            if ( !t0.getValue( 10, 1 ).equals( "110" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 1, 10 ).equals( "101" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 10, 10 ).equals( "1010" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 100, 10 ).equals( "10100" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 9, 1 ).equals( "19" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 49, 4 ).equals( "" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 22349, 3434344 ).equals( "" ) ) {
                return false;
            }
            final GeneralTable<String, String> t1 = new GeneralTable<String, String>();
            t1.setValue( "3", "2", "23" );
            t1.setValue( "10", "1", "error" );
            t1.setValue( "10", "1", "110" );
            t1.setValue( "9", "1", "19" );
            t1.setValue( "1", "10", "101" );
            t1.setValue( "10", "10", "1010" );
            t1.setValue( "100", "10", "10100" );
            t1.setValue( "0", "0", "00" );
            t1.setValue( "qwerty", "zxcvbnm", "asdef" );
            if ( !t1.getValue( "3", "2" ).equals( "23" ) ) {
                return false;
            }
            if ( !t1.getValue( "10", "1" ).equals( "110" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "1", "10" ).equals( "101" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "10", "10" ).equals( "1010" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "100", "10" ).equals( "10100" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "9", "1" ).equals( "19" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "0", "0" ).equals( "00" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "qwerty", "zxcvbnm" ).equals( "asdef" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "49", "4" ).equals( "" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "22349", "3434344" ).equals( "" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testGetDistance() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "(((A:1,B:2,X:100)ab:3,C:4)abc:5,(D:7,(E:9,F:10)ef:8)def:6)r",
                                                 new NHXParser() )[ 0 ];
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "C" ), p1.getNode( "C" ) ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "def" ), p1.getNode( "def" ) ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "ef" ), p1.getNode( "ef" ) ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "r" ), p1.getNode( "r" ) ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "A" ) ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "B" ) ) != 3 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "B" ), p1.getNode( "A" ) ) != 3 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "C" ) ) != 8 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "C" ), p1.getNode( "A" ) ) != 8 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "D" ) ) != 22 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "E" ) ) != 32 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "E" ), p1.getNode( "A" ) ) != 32 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "F" ) ) != 33 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "F" ), p1.getNode( "A" ) ) != 33 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "ab" ) ) != 1 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "ab" ), p1.getNode( "A" ) ) != 1 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "abc" ) ) != 4 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "abc" ), p1.getNode( "A" ) ) != 4 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "r" ) ) != 9 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "r" ), p1.getNode( "A" ) ) != 9 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "def" ) ) != 15 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "def" ), p1.getNode( "A" ) ) != 15 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "A" ), p1.getNode( "ef" ) ) != 23 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "ef" ), p1.getNode( "A" ) ) != 23 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "ef" ), p1.getNode( "def" ) ) != 8 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "def" ), p1.getNode( "ef" ) ) != 8 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "ef" ), p1.getNode( "r" ) ) != 14 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "ef" ), p1.getNode( "abc" ) ) != 19 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "ef" ), p1.getNode( "ab" ) ) != 22 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "ab" ), p1.getNode( "ef" ) ) != 22 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p1.getNode( "def" ), p1.getNode( "abc" ) ) != 11 ) {
                return false;
            }
            final Phylogeny p2 = factory.create( "((A:4,B:5,C:6)abc:1,(D:7,E:8,F:9)def:2,(G:10,H:11,I:12)ghi:3)r",
                                                 new NHXParser() )[ 0 ];
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "A" ), p2.getNode( "B" ) ) != 9 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "A" ), p2.getNode( "C" ) ) != 10 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "A" ), p2.getNode( "D" ) ) != 14 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "A" ), p2.getNode( "ghi" ) ) != 8 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "A" ), p2.getNode( "I" ) ) != 20 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "G" ), p2.getNode( "ghi" ) ) != 10 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "r" ), p2.getNode( "r" ) ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "r" ), p2.getNode( "G" ) ) != 13 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "G" ), p2.getNode( "r" ) ) != 13 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "G" ), p2.getNode( "H" ) ) != 21 ) {
                return false;
            }
            if ( PhylogenyMethods.calculateDistance( p2.getNode( "G" ), p2.getNode( "I" ) ) != 22 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testGetLCA() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "((((((A,B)ab,C)abc,D)abcd,E)abcde,F)abcdef,(G,H)gh)abcdefgh",
                                                 new NHXParser() )[ 0 ];
            final PhylogenyNode A = PhylogenyMethods.calculateLCA( p1.getNode( "A" ), p1.getNode( "A" ) );
            if ( !A.getName().equals( "A" ) ) {
                return false;
            }
            final PhylogenyNode gh = PhylogenyMethods.calculateLCA( p1.getNode( "gh" ), p1.getNode( "gh" ) );
            if ( !gh.getName().equals( "gh" ) ) {
                return false;
            }
            final PhylogenyNode ab = PhylogenyMethods.calculateLCA( p1.getNode( "A" ), p1.getNode( "B" ) );
            if ( !ab.getName().equals( "ab" ) ) {
                return false;
            }
            final PhylogenyNode ab2 = PhylogenyMethods.calculateLCA( p1.getNode( "B" ), p1.getNode( "A" ) );
            if ( !ab2.getName().equals( "ab" ) ) {
                return false;
            }
            final PhylogenyNode gh2 = PhylogenyMethods.calculateLCA( p1.getNode( "H" ), p1.getNode( "G" ) );
            if ( !gh2.getName().equals( "gh" ) ) {
                return false;
            }
            final PhylogenyNode gh3 = PhylogenyMethods.calculateLCA( p1.getNode( "G" ), p1.getNode( "H" ) );
            if ( !gh3.getName().equals( "gh" ) ) {
                return false;
            }
            final PhylogenyNode abc = PhylogenyMethods.calculateLCA( p1.getNode( "C" ), p1.getNode( "A" ) );
            if ( !abc.getName().equals( "abc" ) ) {
                return false;
            }
            final PhylogenyNode abc2 = PhylogenyMethods.calculateLCA( p1.getNode( "A" ), p1.getNode( "C" ) );
            if ( !abc2.getName().equals( "abc" ) ) {
                return false;
            }
            final PhylogenyNode abcd = PhylogenyMethods.calculateLCA( p1.getNode( "A" ), p1.getNode( "D" ) );
            if ( !abcd.getName().equals( "abcd" ) ) {
                return false;
            }
            final PhylogenyNode abcd2 = PhylogenyMethods.calculateLCA( p1.getNode( "D" ), p1.getNode( "A" ) );
            if ( !abcd2.getName().equals( "abcd" ) ) {
                return false;
            }
            final PhylogenyNode abcdef = PhylogenyMethods.calculateLCA( p1.getNode( "A" ), p1.getNode( "F" ) );
            if ( !abcdef.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcdef2 = PhylogenyMethods.calculateLCA( p1.getNode( "F" ), p1.getNode( "A" ) );
            if ( !abcdef2.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcdef3 = PhylogenyMethods.calculateLCA( p1.getNode( "ab" ), p1.getNode( "F" ) );
            if ( !abcdef3.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcdef4 = PhylogenyMethods.calculateLCA( p1.getNode( "F" ), p1.getNode( "ab" ) );
            if ( !abcdef4.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcde = PhylogenyMethods.calculateLCA( p1.getNode( "A" ), p1.getNode( "E" ) );
            if ( !abcde.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode abcde2 = PhylogenyMethods.calculateLCA( p1.getNode( "E" ), p1.getNode( "A" ) );
            if ( !abcde2.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode r = PhylogenyMethods.calculateLCA( p1.getNode( "abcdefgh" ), p1.getNode( "abcdefgh" ) );
            if ( !r.getName().equals( "abcdefgh" ) ) {
                return false;
            }
            final PhylogenyNode r2 = PhylogenyMethods.calculateLCA( p1.getNode( "A" ), p1.getNode( "H" ) );
            if ( !r2.getName().equals( "abcdefgh" ) ) {
                return false;
            }
            final PhylogenyNode r3 = PhylogenyMethods.calculateLCA( p1.getNode( "H" ), p1.getNode( "A" ) );
            if ( !r3.getName().equals( "abcdefgh" ) ) {
                return false;
            }
            final PhylogenyNode abcde3 = PhylogenyMethods.calculateLCA( p1.getNode( "E" ), p1.getNode( "abcde" ) );
            if ( !abcde3.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode abcde4 = PhylogenyMethods.calculateLCA( p1.getNode( "abcde" ), p1.getNode( "E" ) );
            if ( !abcde4.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode ab3 = PhylogenyMethods.calculateLCA( p1.getNode( "ab" ), p1.getNode( "B" ) );
            if ( !ab3.getName().equals( "ab" ) ) {
                return false;
            }
            final PhylogenyNode ab4 = PhylogenyMethods.calculateLCA( p1.getNode( "B" ), p1.getNode( "ab" ) );
            if ( !ab4.getName().equals( "ab" ) ) {
                return false;
            }
            final Phylogeny p2 = factory.create( "(a,b,(((c,d)cd,e)cde,f)cdef)r", new NHXParser() )[ 0 ];
            final PhylogenyNode cd = PhylogenyMethods.calculateLCA( p2.getNode( "c" ), p2.getNode( "d" ) );
            if ( !cd.getName().equals( "cd" ) ) {
                return false;
            }
            final PhylogenyNode cd2 = PhylogenyMethods.calculateLCA( p2.getNode( "d" ), p2.getNode( "c" ) );
            if ( !cd2.getName().equals( "cd" ) ) {
                return false;
            }
            final PhylogenyNode cde = PhylogenyMethods.calculateLCA( p2.getNode( "c" ), p2.getNode( "e" ) );
            if ( !cde.getName().equals( "cde" ) ) {
                return false;
            }
            final PhylogenyNode cde2 = PhylogenyMethods.calculateLCA( p2.getNode( "e" ), p2.getNode( "c" ) );
            if ( !cde2.getName().equals( "cde" ) ) {
                return false;
            }
            final PhylogenyNode cdef = PhylogenyMethods.calculateLCA( p2.getNode( "c" ), p2.getNode( "f" ) );
            if ( !cdef.getName().equals( "cdef" ) ) {
                return false;
            }
            final PhylogenyNode cdef2 = PhylogenyMethods.calculateLCA( p2.getNode( "d" ), p2.getNode( "f" ) );
            if ( !cdef2.getName().equals( "cdef" ) ) {
                return false;
            }
            final PhylogenyNode cdef3 = PhylogenyMethods.calculateLCA( p2.getNode( "f" ), p2.getNode( "d" ) );
            if ( !cdef3.getName().equals( "cdef" ) ) {
                return false;
            }
            final PhylogenyNode rt = PhylogenyMethods.calculateLCA( p2.getNode( "c" ), p2.getNode( "a" ) );
            if ( !rt.getName().equals( "r" ) ) {
                return false;
            }
            final Phylogeny p3 = factory
                    .create( "((((a,(b,c)bc)abc,(d,e)de)abcde,f)abcdef,(((g,h)gh,(i,j)ij)ghij,k)ghijk,l)",
                             new NHXParser() )[ 0 ];
            final PhylogenyNode bc_3 = PhylogenyMethods.calculateLCA( p3.getNode( "b" ), p3.getNode( "c" ) );
            if ( !bc_3.getName().equals( "bc" ) ) {
                return false;
            }
            final PhylogenyNode ac_3 = PhylogenyMethods.calculateLCA( p3.getNode( "a" ), p3.getNode( "c" ) );
            if ( !ac_3.getName().equals( "abc" ) ) {
                return false;
            }
            final PhylogenyNode ad_3 = PhylogenyMethods.calculateLCA( p3.getNode( "a" ), p3.getNode( "d" ) );
            if ( !ad_3.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode af_3 = PhylogenyMethods.calculateLCA( p3.getNode( "a" ), p3.getNode( "f" ) );
            if ( !af_3.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode ag_3 = PhylogenyMethods.calculateLCA( p3.getNode( "a" ), p3.getNode( "g" ) );
            if ( !ag_3.getName().equals( "" ) ) {
                return false;
            }
            if ( !ag_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode al_3 = PhylogenyMethods.calculateLCA( p3.getNode( "a" ), p3.getNode( "l" ) );
            if ( !al_3.getName().equals( "" ) ) {
                return false;
            }
            if ( !al_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode kl_3 = PhylogenyMethods.calculateLCA( p3.getNode( "k" ), p3.getNode( "l" ) );
            if ( !kl_3.getName().equals( "" ) ) {
                return false;
            }
            if ( !kl_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode fl_3 = PhylogenyMethods.calculateLCA( p3.getNode( "f" ), p3.getNode( "l" ) );
            if ( !fl_3.getName().equals( "" ) ) {
                return false;
            }
            if ( !fl_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode gk_3 = PhylogenyMethods.calculateLCA( p3.getNode( "g" ), p3.getNode( "k" ) );
            if ( !gk_3.getName().equals( "ghijk" ) ) {
                return false;
            }
            final Phylogeny p4 = factory.create( "(a,b,c)r", new NHXParser() )[ 0 ];
            final PhylogenyNode r_4 = PhylogenyMethods.calculateLCA( p4.getNode( "b" ), p4.getNode( "c" ) );
            if ( !r_4.getName().equals( "r" ) ) {
                return false;
            }
            final Phylogeny p5 = factory.create( "((a,b),c,d)root", new NHXParser() )[ 0 ];
            final PhylogenyNode r_5 = PhylogenyMethods.calculateLCA( p5.getNode( "a" ), p5.getNode( "c" ) );
            if ( !r_5.getName().equals( "root" ) ) {
                return false;
            }
            final Phylogeny p6 = factory.create( "((a,b),c,d)rot", new NHXParser() )[ 0 ];
            final PhylogenyNode r_6 = PhylogenyMethods.calculateLCA( p6.getNode( "c" ), p6.getNode( "a" ) );
            if ( !r_6.getName().equals( "rot" ) ) {
                return false;
            }
            final Phylogeny p7 = factory.create( "(((a,b)x,c)x,d,e)rott", new NHXParser() )[ 0 ];
            final PhylogenyNode r_7 = PhylogenyMethods.calculateLCA( p7.getNode( "a" ), p7.getNode( "e" ) );
            if ( !r_7.getName().equals( "rott" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testGetLCA2() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p_a = factory.create( "(a)", new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p_a );
            final PhylogenyNode p_a_1 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p_a.getNode( "a" ),
                                                                                              p_a.getNode( "a" ) );
            if ( !p_a_1.getName().equals( "a" ) ) {
                return false;
            }
            final Phylogeny p_b = factory.create( "((a)b)", new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p_b );
            final PhylogenyNode p_b_1 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p_b.getNode( "b" ),
                                                                                              p_b.getNode( "a" ) );
            if ( !p_b_1.getName().equals( "b" ) ) {
                return false;
            }
            final PhylogenyNode p_b_2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p_b.getNode( "a" ),
                                                                                              p_b.getNode( "b" ) );
            if ( !p_b_2.getName().equals( "b" ) ) {
                return false;
            }
            final Phylogeny p_c = factory.create( "(((a)b)c)", new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p_c );
            final PhylogenyNode p_c_1 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p_c.getNode( "b" ),
                                                                                              p_c.getNode( "a" ) );
            if ( !p_c_1.getName().equals( "b" ) ) {
                return false;
            }
            final PhylogenyNode p_c_2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p_c.getNode( "a" ),
                                                                                              p_c.getNode( "c" ) );
            if ( !p_c_2.getName().equals( "c" ) ) {
                System.out.println( p_c_2.getName() );
                System.exit( -1 );
                return false;
            }
            final PhylogenyNode p_c_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p_c.getNode( "a" ),
                                                                                              p_c.getNode( "b" ) );
            if ( !p_c_3.getName().equals( "b" ) ) {
                return false;
            }
            final PhylogenyNode p_c_4 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p_c.getNode( "c" ),
                                                                                              p_c.getNode( "a" ) );
            if ( !p_c_4.getName().equals( "c" ) ) {
                return false;
            }
            final Phylogeny p1 = factory.create( "((((((A,B)ab,C)abc,D)abcd,E)abcde,F)abcdef,(G,H)gh)abcdefgh",
                                                 new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p1 );
            final PhylogenyNode A = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "A" ),
                                                                                          p1.getNode( "A" ) );
            if ( !A.getName().equals( "A" ) ) {
                return false;
            }
            final PhylogenyNode gh = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "gh" ),
                                                                                           p1.getNode( "gh" ) );
            if ( !gh.getName().equals( "gh" ) ) {
                return false;
            }
            final PhylogenyNode ab = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "A" ),
                                                                                           p1.getNode( "B" ) );
            if ( !ab.getName().equals( "ab" ) ) {
                return false;
            }
            final PhylogenyNode ab2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "B" ),
                                                                                            p1.getNode( "A" ) );
            if ( !ab2.getName().equals( "ab" ) ) {
                return false;
            }
            final PhylogenyNode gh2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "H" ),
                                                                                            p1.getNode( "G" ) );
            if ( !gh2.getName().equals( "gh" ) ) {
                return false;
            }
            final PhylogenyNode gh3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "G" ),
                                                                                            p1.getNode( "H" ) );
            if ( !gh3.getName().equals( "gh" ) ) {
                return false;
            }
            final PhylogenyNode abc = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "C" ),
                                                                                            p1.getNode( "A" ) );
            if ( !abc.getName().equals( "abc" ) ) {
                return false;
            }
            final PhylogenyNode abc2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "A" ),
                                                                                             p1.getNode( "C" ) );
            if ( !abc2.getName().equals( "abc" ) ) {
                return false;
            }
            final PhylogenyNode abcd = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "A" ),
                                                                                             p1.getNode( "D" ) );
            if ( !abcd.getName().equals( "abcd" ) ) {
                return false;
            }
            final PhylogenyNode abcd2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "D" ),
                                                                                              p1.getNode( "A" ) );
            if ( !abcd2.getName().equals( "abcd" ) ) {
                return false;
            }
            final PhylogenyNode abcdef = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "A" ),
                                                                                               p1.getNode( "F" ) );
            if ( !abcdef.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcdef2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "F" ),
                                                                                                p1.getNode( "A" ) );
            if ( !abcdef2.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcdef3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "ab" ),
                                                                                                p1.getNode( "F" ) );
            if ( !abcdef3.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcdef4 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "F" ),
                                                                                                p1.getNode( "ab" ) );
            if ( !abcdef4.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcde = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "A" ),
                                                                                              p1.getNode( "E" ) );
            if ( !abcde.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode abcde2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "E" ),
                                                                                               p1.getNode( "A" ) );
            if ( !abcde2.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode r = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "abcdefgh" ),
                                                                                          p1.getNode( "abcdefgh" ) );
            if ( !r.getName().equals( "abcdefgh" ) ) {
                return false;
            }
            final PhylogenyNode r2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "A" ),
                                                                                           p1.getNode( "H" ) );
            if ( !r2.getName().equals( "abcdefgh" ) ) {
                return false;
            }
            final PhylogenyNode r3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "H" ),
                                                                                           p1.getNode( "A" ) );
            if ( !r3.getName().equals( "abcdefgh" ) ) {
                return false;
            }
            final PhylogenyNode abcde3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "E" ),
                                                                                               p1.getNode( "abcde" ) );
            if ( !abcde3.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode abcde4 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "abcde" ),
                                                                                               p1.getNode( "E" ) );
            if ( !abcde4.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode ab3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "ab" ),
                                                                                            p1.getNode( "B" ) );
            if ( !ab3.getName().equals( "ab" ) ) {
                return false;
            }
            final PhylogenyNode ab4 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p1.getNode( "B" ),
                                                                                            p1.getNode( "ab" ) );
            if ( !ab4.getName().equals( "ab" ) ) {
                return false;
            }
            final Phylogeny p2 = factory.create( "(a,b,(((c,d)cd,e)cde,f)cdef)r", new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p2 );
            final PhylogenyNode cd = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p2.getNode( "c" ),
                                                                                           p2.getNode( "d" ) );
            if ( !cd.getName().equals( "cd" ) ) {
                return false;
            }
            final PhylogenyNode cd2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p2.getNode( "d" ),
                                                                                            p2.getNode( "c" ) );
            if ( !cd2.getName().equals( "cd" ) ) {
                return false;
            }
            final PhylogenyNode cde = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p2.getNode( "c" ),
                                                                                            p2.getNode( "e" ) );
            if ( !cde.getName().equals( "cde" ) ) {
                return false;
            }
            final PhylogenyNode cde2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p2.getNode( "e" ),
                                                                                             p2.getNode( "c" ) );
            if ( !cde2.getName().equals( "cde" ) ) {
                return false;
            }
            final PhylogenyNode cdef = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p2.getNode( "c" ),
                                                                                             p2.getNode( "f" ) );
            if ( !cdef.getName().equals( "cdef" ) ) {
                return false;
            }
            final PhylogenyNode cdef2 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p2.getNode( "d" ),
                                                                                              p2.getNode( "f" ) );
            if ( !cdef2.getName().equals( "cdef" ) ) {
                return false;
            }
            final PhylogenyNode cdef3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p2.getNode( "f" ),
                                                                                              p2.getNode( "d" ) );
            if ( !cdef3.getName().equals( "cdef" ) ) {
                return false;
            }
            final PhylogenyNode rt = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p2.getNode( "c" ),
                                                                                           p2.getNode( "a" ) );
            if ( !rt.getName().equals( "r" ) ) {
                return false;
            }
            final Phylogeny p3 = factory
                    .create( "((((a,(b,c)bc)abc,(d,e)de)abcde,f)abcdef,(((g,h)gh,(i,j)ij)ghij,k)ghijk,l)",
                             new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p3 );
            final PhylogenyNode bc_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p3.getNode( "b" ),
                                                                                             p3.getNode( "c" ) );
            if ( !bc_3.getName().equals( "bc" ) ) {
                return false;
            }
            final PhylogenyNode ac_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p3.getNode( "a" ),
                                                                                             p3.getNode( "c" ) );
            if ( !ac_3.getName().equals( "abc" ) ) {
                return false;
            }
            final PhylogenyNode ad_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p3.getNode( "a" ),
                                                                                             p3.getNode( "d" ) );
            if ( !ad_3.getName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode af_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p3.getNode( "a" ),
                                                                                             p3.getNode( "f" ) );
            if ( !af_3.getName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode ag_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p3.getNode( "a" ),
                                                                                             p3.getNode( "g" ) );
            if ( !ag_3.getName().equals( "" ) ) {
                return false;
            }
            if ( !ag_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode al_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p3.getNode( "a" ),
                                                                                             p3.getNode( "l" ) );
            if ( !al_3.getName().equals( "" ) ) {
                return false;
            }
            if ( !al_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode kl_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p3.getNode( "k" ),
                                                                                             p3.getNode( "l" ) );
            if ( !kl_3.getName().equals( "" ) ) {
                return false;
            }
            if ( !kl_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode fl_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p3.getNode( "f" ),
                                                                                             p3.getNode( "l" ) );
            if ( !fl_3.getName().equals( "" ) ) {
                return false;
            }
            if ( !fl_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode gk_3 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p3.getNode( "g" ),
                                                                                             p3.getNode( "k" ) );
            if ( !gk_3.getName().equals( "ghijk" ) ) {
                return false;
            }
            final Phylogeny p4 = factory.create( "(a,b,c)r", new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p4 );
            final PhylogenyNode r_4 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p4.getNode( "b" ),
                                                                                            p4.getNode( "c" ) );
            if ( !r_4.getName().equals( "r" ) ) {
                return false;
            }
            final Phylogeny p5 = factory.create( "((a,b),c,d)root", new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p5 );
            final PhylogenyNode r_5 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p5.getNode( "a" ),
                                                                                            p5.getNode( "c" ) );
            if ( !r_5.getName().equals( "root" ) ) {
                return false;
            }
            final Phylogeny p6 = factory.create( "((a,b),c,d)rot", new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p6 );
            final PhylogenyNode r_6 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p6.getNode( "c" ),
                                                                                            p6.getNode( "a" ) );
            if ( !r_6.getName().equals( "rot" ) ) {
                return false;
            }
            final Phylogeny p7 = factory.create( "(((a,b)x,c)x,d,e)rott", new NHXParser() )[ 0 ];
            PhylogenyMethods.preOrderReId( p7 );
            final PhylogenyNode r_7 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p7.getNode( "a" ),
                                                                                            p7.getNode( "e" ) );
            if ( !r_7.getName().equals( "rott" ) ) {
                return false;
            }
            final PhylogenyNode r_71 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p7.getNode( "e" ),
                                                                                             p7.getNode( "a" ) );
            if ( !r_71.getName().equals( "rott" ) ) {
                return false;
            }
            final PhylogenyNode r_72 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p7.getNode( "e" ),
                                                                                             p7.getNode( "rott" ) );
            if ( !r_72.getName().equals( "rott" ) ) {
                return false;
            }
            final PhylogenyNode r_73 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p7.getNode( "rott" ),
                                                                                             p7.getNode( "a" ) );
            if ( !r_73.getName().equals( "rott" ) ) {
                return false;
            }
            final PhylogenyNode r_74 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p7.getNode( "rott" ),
                                                                                             p7.getNode( "rott" ) );
            if ( !r_74.getName().equals( "rott" ) ) {
                return false;
            }
            final PhylogenyNode r_75 = PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( p7.getNode( "e" ),
                                                                                             p7.getNode( "e" ) );
            if ( !r_75.getName().equals( "e" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testHmmscanOutputParser() {
        final String test_dir = Test.PATH_TO_TEST_DATA;
        try {
            final HmmscanPerDomainTableParser parser1 = new HmmscanPerDomainTableParser( new File( test_dir
                    + ForesterUtil.getFileSeparator() + "hmmscan30b3_output_1" ), "MONBR", INDIVIDUAL_SCORE_CUTOFF.NONE );
            parser1.parse();
            final HmmscanPerDomainTableParser parser2 = new HmmscanPerDomainTableParser( new File( test_dir
                    + ForesterUtil.getFileSeparator() + "hmmscan30b3_output_2" ), "MONBR", INDIVIDUAL_SCORE_CUTOFF.NONE );
            final List<Protein> proteins = parser2.parse();
            if ( parser2.getProteinsEncountered() != 4 ) {
                return false;
            }
            if ( proteins.size() != 4 ) {
                return false;
            }
            if ( parser2.getDomainsEncountered() != 69 ) {
                return false;
            }
            if ( parser2.getDomainsIgnoredDueToDuf() != 0 ) {
                return false;
            }
            if ( parser2.getDomainsIgnoredDueToEval() != 0 ) {
                return false;
            }
            final Protein p1 = proteins.get( 0 );
            if ( p1.getNumberOfProteinDomains() != 15 ) {
                return false;
            }
            if ( p1.getLength() != 850 ) {
                return false;
            }
            final Protein p2 = proteins.get( 1 );
            if ( p2.getNumberOfProteinDomains() != 51 ) {
                return false;
            }
            if ( p2.getLength() != 1291 ) {
                return false;
            }
            final Protein p3 = proteins.get( 2 );
            if ( p3.getNumberOfProteinDomains() != 2 ) {
                return false;
            }
            final Protein p4 = proteins.get( 3 );
            if ( p4.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( !p4.getProteinDomain( 0 ).getDomainId().toString().equals( "DNA_pol_B_new" ) ) {
                return false;
            }
            if ( p4.getProteinDomain( 0 ).getFrom() != 51 ) {
                return false;
            }
            if ( p4.getProteinDomain( 0 ).getTo() != 395 ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getPerDomainEvalue(), 1.2e-39 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getPerDomainScore(), 135.7 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getPerSequenceEvalue(), 8.3e-40 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getPerSequenceScore(), 136.3 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getNumber(), 1 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getTotalCount(), 1 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testLastExternalNodeMethods() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final char[] a0 = { '(', '(', 'A', ',', 'B', ')', ',', '(', 'C', ',', 'D', ')', ')', };
            final Phylogeny t0 = factory.create( a0, new NHXParser() )[ 0 ];
            final PhylogenyNode n1 = t0.getNode( "A" );
            if ( n1.isLastExternalNode() ) {
                return false;
            }
            final PhylogenyNode n2 = t0.getNode( "B" );
            if ( n2.isLastExternalNode() ) {
                return false;
            }
            final PhylogenyNode n3 = t0.getNode( "C" );
            if ( n3.isLastExternalNode() ) {
                return false;
            }
            final PhylogenyNode n4 = t0.getNode( "D" );
            if ( !n4.isLastExternalNode() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testLevelOrderIterator() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((A,B)ab,(C,D)cd)r", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it0;
            for( it0 = t0.iteratorLevelOrder(); it0.hasNext(); ) {
                it0.next();
            }
            for( it0.reset(); it0.hasNext(); ) {
                it0.next();
            }
            final PhylogenyNodeIterator it = t0.iteratorLevelOrder();
            if ( !it.next().getName().equals( "r" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "ab" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "cd" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "A" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "B" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "C" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "D" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "(((1,2,(a,(X,Y,Z)b)3,4,5,6)A,B,C)abc,(D,E,(f1,(f21)f2,f3)F,G)defg)r",
                                                 new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it2;
            for( it2 = t2.iteratorLevelOrder(); it2.hasNext(); ) {
                it2.next();
            }
            for( it2.reset(); it2.hasNext(); ) {
                it2.next();
            }
            final PhylogenyNodeIterator it3 = t2.iteratorLevelOrder();
            if ( !it3.next().getName().equals( "r" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "abc" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "defg" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "A" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "B" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "C" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "D" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "E" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "F" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "G" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "1" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "2" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "3" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "4" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "5" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "6" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "f1" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "f2" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "f3" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "a" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "b" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "f21" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "X" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "Y" ) ) {
                return false;
            }
            if ( !it3.next().getName().equals( "Z" ) ) {
                return false;
            }
            if ( it3.hasNext() ) {
                return false;
            }
            final Phylogeny t4 = factory.create( "((((D)C)B)A)r", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it4;
            for( it4 = t4.iteratorLevelOrder(); it4.hasNext(); ) {
                it4.next();
            }
            for( it4.reset(); it4.hasNext(); ) {
                it4.next();
            }
            final PhylogenyNodeIterator it5 = t4.iteratorLevelOrder();
            if ( !it5.next().getName().equals( "r" ) ) {
                return false;
            }
            if ( !it5.next().getName().equals( "A" ) ) {
                return false;
            }
            if ( !it5.next().getName().equals( "B" ) ) {
                return false;
            }
            if ( !it5.next().getName().equals( "C" ) ) {
                return false;
            }
            if ( !it5.next().getName().equals( "D" ) ) {
                return false;
            }
            final Phylogeny t5 = factory.create( "A", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it6;
            for( it6 = t5.iteratorLevelOrder(); it6.hasNext(); ) {
                it6.next();
            }
            for( it6.reset(); it6.hasNext(); ) {
                it6.next();
            }
            final PhylogenyNodeIterator it7 = t5.iteratorLevelOrder();
            if ( !it7.next().getName().equals( "A" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testMafft( final String path ) {
        try {
            final List<String> opts = new ArrayList<String>();
            opts.add( "--maxiterate" );
            opts.add( "1000" );
            opts.add( "--localpair" );
            opts.add( "--quiet" );
            Msa msa = null;
            final MsaInferrer mafft = Mafft.createInstance( path );
            msa = mafft.infer( new File( PATH_TO_TEST_DATA + "ncbi_sn.fasta" ), opts );
            if ( ( msa == null ) || ( msa.getLength() < 20 ) || ( msa.getNumberOfSequences() != 19 ) ) {
                return false;
            }
            if ( !msa.getIdentifier( 0 ).toString().equals( "a" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testMidpointrooting() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "(A:1,B:4,C:2,D:2,E:6,F:1,G:1,H:1)", new NHXParser() )[ 0 ];
            PhylogenyMethods.midpointRoot( t0 );
            if ( !isEqual( t0.getNode( "E" ).getDistanceToParent(), 5 ) ) {
                return false;
            }
            if ( !isEqual( t0.getNode( "B" ).getDistanceToParent(), 4 ) ) {
                return false;
            }
            if ( !isEqual( PhylogenyMethods.calculateLCA( t0.getNode( "F" ), t0.getNode( "G" ) ).getDistanceToParent(),
                           1 ) ) {
                return false;
            }
            final Phylogeny t1 = factory.create( "((A:1,B:2)AB:1[&&NHX:B=55],(C:3,D:4)CD:3[&&NHX:B=10])ABCD:0.5",
                                                 new NHXParser() )[ 0 ];
            if ( !t1.isRooted() ) {
                return false;
            }
            PhylogenyMethods.midpointRoot( t1 );
            if ( !isEqual( t1.getNode( "A" ).getDistanceToParent(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "B" ).getDistanceToParent(), 2 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "C" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "D" ).getDistanceToParent(), 4 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "CD" ).getDistanceToParent(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "AB" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
            t1.reRoot( t1.getNode( "A" ) );
            PhylogenyMethods.midpointRoot( t1 );
            if ( !isEqual( t1.getNode( "A" ).getDistanceToParent(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "B" ).getDistanceToParent(), 2 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "C" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "D" ).getDistanceToParent(), 4 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "CD" ).getDistanceToParent(), 1 ) ) {
                System.exit( -1 );
                return false;
            }
            if ( !isEqual( t1.getNode( "AB" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testMsaQualityMethod() {
        try {
            final Sequence s0 = BasicSequence.createAaSequence( "a", "ABAXEFGHIJ" );
            final Sequence s1 = BasicSequence.createAaSequence( "b", "ABBXEFGHIJ" );
            final Sequence s2 = BasicSequence.createAaSequence( "c", "AXCXEFGHIJ" );
            final Sequence s3 = BasicSequence.createAaSequence( "d", "AXDDEFGHIJ" );
            final List<Sequence> l = new ArrayList<Sequence>();
            l.add( s0 );
            l.add( s1 );
            l.add( s2 );
            l.add( s3 );
            final Msa msa = BasicMsa.createInstance( l );
            if ( !isEqual( 1, MsaMethods.calculateIdentityRatio( msa, 0 ) ) ) {
                return false;
            }
            if ( !isEqual( 0.5, MsaMethods.calculateIdentityRatio( msa, 1 ) ) ) {
                return false;
            }
            if ( !isEqual( 0.25, MsaMethods.calculateIdentityRatio( msa, 2 ) ) ) {
                return false;
            }
            if ( !isEqual( 0.75, MsaMethods.calculateIdentityRatio( msa, 3 ) ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNextNodeWithCollapsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            PhylogenyNode n;
            List<PhylogenyNode> ext = new ArrayList<PhylogenyNode>();
            final StringBuffer sb0 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h))fgh)cdefgh)abcdefgh" );
            final Phylogeny t0 = factory.create( sb0, new NHXParser() )[ 0 ];
            t0.getNode( "cd" ).setCollapse( true );
            t0.getNode( "cde" ).setCollapse( true );
            n = t0.getFirstExternalNode();
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( !ext.get( 0 ).getName().equals( "a" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "b" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "cde" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "f" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "g" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "h" ) ) {
                return false;
            }
            ext.clear();
            final StringBuffer sb1 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h))fgh)cdefgh)abcdefgh" );
            final Phylogeny t1 = factory.create( sb1, new NHXParser() )[ 0 ];
            t1.getNode( "ab" ).setCollapse( true );
            t1.getNode( "cd" ).setCollapse( true );
            t1.getNode( "cde" ).setCollapse( true );
            n = t1.getNode( "ab" );
            ext = new ArrayList<PhylogenyNode>();
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( !ext.get( 0 ).getName().equals( "ab" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "cde" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "f" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "g" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "h" ) ) {
                return false;
            }
            //
            //
            ext.clear();
            final StringBuffer sb2 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h)gh)fgh)cdefgh)abcdefgh" );
            final Phylogeny t2 = factory.create( sb2, new NHXParser() )[ 0 ];
            t2.getNode( "ab" ).setCollapse( true );
            t2.getNode( "cd" ).setCollapse( true );
            t2.getNode( "cde" ).setCollapse( true );
            t2.getNode( "c" ).setCollapse( true );
            t2.getNode( "d" ).setCollapse( true );
            t2.getNode( "e" ).setCollapse( true );
            t2.getNode( "gh" ).setCollapse( true );
            n = t2.getNode( "ab" );
            ext = new ArrayList<PhylogenyNode>();
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( !ext.get( 0 ).getName().equals( "ab" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "cde" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "f" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "gh" ) ) {
                return false;
            }
            //
            //
            ext.clear();
            final StringBuffer sb3 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h)gh)fgh)cdefgh)abcdefgh" );
            final Phylogeny t3 = factory.create( sb3, new NHXParser() )[ 0 ];
            t3.getNode( "ab" ).setCollapse( true );
            t3.getNode( "cd" ).setCollapse( true );
            t3.getNode( "cde" ).setCollapse( true );
            t3.getNode( "c" ).setCollapse( true );
            t3.getNode( "d" ).setCollapse( true );
            t3.getNode( "e" ).setCollapse( true );
            t3.getNode( "gh" ).setCollapse( true );
            t3.getNode( "fgh" ).setCollapse( true );
            n = t3.getNode( "ab" );
            ext = new ArrayList<PhylogenyNode>();
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( !ext.get( 0 ).getName().equals( "ab" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "cde" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "fgh" ) ) {
                return false;
            }
            //
            //
            ext.clear();
            final StringBuffer sb4 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h)gh)fgh)cdefgh)abcdefgh" );
            final Phylogeny t4 = factory.create( sb4, new NHXParser() )[ 0 ];
            t4.getNode( "ab" ).setCollapse( true );
            t4.getNode( "cd" ).setCollapse( true );
            t4.getNode( "cde" ).setCollapse( true );
            t4.getNode( "c" ).setCollapse( true );
            t4.getNode( "d" ).setCollapse( true );
            t4.getNode( "e" ).setCollapse( true );
            t4.getNode( "gh" ).setCollapse( true );
            t4.getNode( "fgh" ).setCollapse( true );
            t4.getNode( "abcdefgh" ).setCollapse( true );
            n = t4.getNode( "abcdefgh" );
            if ( n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes() != null ) {
                return false;
            }
            //
            //
            final StringBuffer sb5 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h))fgh)cdefgh)abcdefgh" );
            final Phylogeny t5 = factory.create( sb5, new NHXParser() )[ 0 ];
            ext.clear();
            n = t5.getFirstExternalNode();
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 8 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "a" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "b" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "c" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "d" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "f" ) ) {
                return false;
            }
            if ( !ext.get( 6 ).getName().equals( "g" ) ) {
                return false;
            }
            if ( !ext.get( 7 ).getName().equals( "h" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb6 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h))fgh)cdefgh)abcdefgh" );
            final Phylogeny t6 = factory.create( sb6, new NHXParser() )[ 0 ];
            ext.clear();
            t6.getNode( "ab" ).setCollapse( true );
            n = t6.getNode( "ab" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 7 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "ab" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "c" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "d" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "f" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "g" ) ) {
                return false;
            }
            if ( !ext.get( 6 ).getName().equals( "h" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb7 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h))fgh)cdefgh)abcdefgh" );
            final Phylogeny t7 = factory.create( sb7, new NHXParser() )[ 0 ];
            ext.clear();
            t7.getNode( "cd" ).setCollapse( true );
            n = t7.getNode( "a" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 7 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "a" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "b" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "cd" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "f" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "g" ) ) {
                return false;
            }
            if ( !ext.get( 6 ).getName().equals( "h" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb8 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h))fgh)cdefgh)abcdefgh" );
            final Phylogeny t8 = factory.create( sb8, new NHXParser() )[ 0 ];
            ext.clear();
            t8.getNode( "cd" ).setCollapse( true );
            t8.getNode( "c" ).setCollapse( true );
            t8.getNode( "d" ).setCollapse( true );
            n = t8.getNode( "a" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 7 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "a" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "b" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "cd" ) ) {
                System.out.println( "2 fail" );
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "f" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "g" ) ) {
                return false;
            }
            if ( !ext.get( 6 ).getName().equals( "h" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb9 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h)gh)fgh)cdefgh)abcdefgh" );
            final Phylogeny t9 = factory.create( sb9, new NHXParser() )[ 0 ];
            ext.clear();
            t9.getNode( "gh" ).setCollapse( true );
            n = t9.getNode( "a" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 7 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "a" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "b" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "c" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "d" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "f" ) ) {
                return false;
            }
            if ( !ext.get( 6 ).getName().equals( "gh" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb10 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h)gh)fgh)cdefgh)abcdefgh" );
            final Phylogeny t10 = factory.create( sb10, new NHXParser() )[ 0 ];
            ext.clear();
            t10.getNode( "gh" ).setCollapse( true );
            t10.getNode( "g" ).setCollapse( true );
            t10.getNode( "h" ).setCollapse( true );
            n = t10.getNode( "a" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 7 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "a" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "b" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "c" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "d" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "f" ) ) {
                return false;
            }
            if ( !ext.get( 6 ).getName().equals( "gh" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb11 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h)gh)fgh)cdefgh)abcdefgh" );
            final Phylogeny t11 = factory.create( sb11, new NHXParser() )[ 0 ];
            ext.clear();
            t11.getNode( "gh" ).setCollapse( true );
            t11.getNode( "fgh" ).setCollapse( true );
            n = t11.getNode( "a" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 6 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "a" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "b" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "c" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "d" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "fgh" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb12 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h)gh)fgh)cdefgh)abcdefgh" );
            final Phylogeny t12 = factory.create( sb12, new NHXParser() )[ 0 ];
            ext.clear();
            t12.getNode( "gh" ).setCollapse( true );
            t12.getNode( "fgh" ).setCollapse( true );
            t12.getNode( "g" ).setCollapse( true );
            t12.getNode( "h" ).setCollapse( true );
            t12.getNode( "f" ).setCollapse( true );
            n = t12.getNode( "a" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 6 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "a" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "b" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "c" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "d" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "fgh" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb13 = new StringBuffer( "((a,b)ab,(((c,d)cd,e)cde,(f,(g,h)gh)fgh)cdefgh)abcdefgh" );
            final Phylogeny t13 = factory.create( sb13, new NHXParser() )[ 0 ];
            ext.clear();
            t13.getNode( "ab" ).setCollapse( true );
            t13.getNode( "b" ).setCollapse( true );
            t13.getNode( "fgh" ).setCollapse( true );
            t13.getNode( "gh" ).setCollapse( true );
            n = t13.getNode( "ab" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 5 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "ab" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "c" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "d" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "fgh" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb14 = new StringBuffer( "((a,b,0)ab,(((c,d)cd,e)cde,(f,(g,h,1,2)gh,0)fgh)cdefgh)abcdefgh" );
            final Phylogeny t14 = factory.create( sb14, new NHXParser() )[ 0 ];
            ext.clear();
            t14.getNode( "ab" ).setCollapse( true );
            t14.getNode( "a" ).setCollapse( true );
            t14.getNode( "fgh" ).setCollapse( true );
            t14.getNode( "gh" ).setCollapse( true );
            n = t14.getNode( "ab" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 5 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "ab" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "c" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "d" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "fgh" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb15 = new StringBuffer( "((a,b,0)ab,(((c,d)cd,e)cde,x,(f,(g,h,1,2)gh,0)fgh)cdefgh)abcdefgh" );
            final Phylogeny t15 = factory.create( sb15, new NHXParser() )[ 0 ];
            ext.clear();
            t15.getNode( "ab" ).setCollapse( true );
            t15.getNode( "a" ).setCollapse( true );
            t15.getNode( "fgh" ).setCollapse( true );
            t15.getNode( "gh" ).setCollapse( true );
            n = t15.getNode( "ab" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 6 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "ab" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "c" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "d" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "e" ) ) {
                return false;
            }
            if ( !ext.get( 4 ).getName().equals( "x" ) ) {
                return false;
            }
            if ( !ext.get( 5 ).getName().equals( "fgh" ) ) {
                return false;
            }
            //
            //
            final StringBuffer sb16 = new StringBuffer( "((a,b,0)ab,(((c,d)cd,e)cde,x,(f,(g,h,1,2)gh,0)fgh)cdefgh)abcdefgh" );
            final Phylogeny t16 = factory.create( sb16, new NHXParser() )[ 0 ];
            ext.clear();
            t16.getNode( "ab" ).setCollapse( true );
            t16.getNode( "a" ).setCollapse( true );
            t16.getNode( "fgh" ).setCollapse( true );
            t16.getNode( "gh" ).setCollapse( true );
            t16.getNode( "cd" ).setCollapse( true );
            t16.getNode( "cde" ).setCollapse( true );
            t16.getNode( "d" ).setCollapse( true );
            t16.getNode( "x" ).setCollapse( true );
            n = t16.getNode( "ab" );
            while ( n != null ) {
                ext.add( n );
                n = n.getNextExternalNodeWhileTakingIntoAccountCollapsedNodes();
            }
            if ( ext.size() != 4 ) {
                return false;
            }
            if ( !ext.get( 0 ).getName().equals( "ab" ) ) {
                return false;
            }
            if ( !ext.get( 1 ).getName().equals( "cde" ) ) {
                return false;
            }
            if ( !ext.get( 2 ).getName().equals( "x" ) ) {
                return false;
            }
            if ( !ext.get( 3 ).getName().equals( "fgh" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNexusCharactersParsing() {
        try {
            final NexusCharactersParser parser = new NexusCharactersParser();
            parser.setSource( new File( Test.PATH_TO_TEST_DATA + "nexus_test_7.nex" ) );
            parser.parse();
            String[] labels = parser.getCharStateLabels();
            if ( labels.length != 7 ) {
                return false;
            }
            if ( !labels[ 0 ].equals( "14-3-3" ) ) {
                return false;
            }
            if ( !labels[ 1 ].equals( "2-Hacid_dh" ) ) {
                return false;
            }
            if ( !labels[ 2 ].equals( "2-Hacid_dh_C" ) ) {
                return false;
            }
            if ( !labels[ 3 ].equals( "2-oxoacid_dh" ) ) {
                return false;
            }
            if ( !labels[ 4 ].equals( "2OG-FeII_Oxy" ) ) {
                return false;
            }
            if ( !labels[ 5 ].equals( "3-HAO" ) ) {
                return false;
            }
            if ( !labels[ 6 ].equals( "3_5_exonuc" ) ) {
                return false;
            }
            parser.setSource( new File( Test.PATH_TO_TEST_DATA + "nexus_test_8.nex" ) );
            parser.parse();
            labels = parser.getCharStateLabels();
            if ( labels.length != 7 ) {
                return false;
            }
            if ( !labels[ 0 ].equals( "14-3-3" ) ) {
                return false;
            }
            if ( !labels[ 1 ].equals( "2-Hacid_dh" ) ) {
                return false;
            }
            if ( !labels[ 2 ].equals( "2-Hacid_dh_C" ) ) {
                return false;
            }
            if ( !labels[ 3 ].equals( "2-oxoacid_dh" ) ) {
                return false;
            }
            if ( !labels[ 4 ].equals( "2OG-FeII_Oxy" ) ) {
                return false;
            }
            if ( !labels[ 5 ].equals( "3-HAO" ) ) {
                return false;
            }
            if ( !labels[ 6 ].equals( "3_5_exonuc" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNexusMatrixParsing() {
        try {
            final NexusBinaryStatesMatrixParser parser = new NexusBinaryStatesMatrixParser();
            parser.setSource( new File( Test.PATH_TO_TEST_DATA + "nexus_test_9.nex" ) );
            parser.parse();
            final CharacterStateMatrix<BinaryStates> m = parser.getMatrix();
            if ( m.getNumberOfCharacters() != 9 ) {
                return false;
            }
            if ( m.getNumberOfIdentifiers() != 5 ) {
                return false;
            }
            if ( m.getState( 0, 0 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( m.getState( 0, 1 ) != BinaryStates.ABSENT ) {
                return false;
            }
            if ( m.getState( 1, 0 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( m.getState( 2, 0 ) != BinaryStates.ABSENT ) {
                return false;
            }
            if ( m.getState( 4, 8 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( !m.getIdentifier( 0 ).equals( "MOUSE" ) ) {
                return false;
            }
            if ( !m.getIdentifier( 4 ).equals( "ARATH" ) ) {
                return false;
            }
            //            if ( labels.length != 7 ) {
            //                return false;
            //            }
            //            if ( !labels[ 0 ].equals( "14-3-3" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 1 ].equals( "2-Hacid_dh" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 2 ].equals( "2-Hacid_dh_C" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 3 ].equals( "2-oxoacid_dh" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 4 ].equals( "2OG-FeII_Oxy" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 5 ].equals( "3-HAO" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 6 ].equals( "3_5_exonuc" ) ) {
            //                return false;
            //            }
            //            parser.setSource( new File( Test.PATH_TO_TEST_DATA + "nexus_test_8.nex" ) );
            //            parser.parse();
            //            labels = parser.getCharStateLabels();
            //            if ( labels.length != 7 ) {
            //                return false;
            //            }
            //            if ( !labels[ 0 ].equals( "14-3-3" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 1 ].equals( "2-Hacid_dh" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 2 ].equals( "2-Hacid_dh_C" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 3 ].equals( "2-oxoacid_dh" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 4 ].equals( "2OG-FeII_Oxy" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 5 ].equals( "3-HAO" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 6 ].equals( "3_5_exonuc" ) ) {
            //                return false;
            //            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNexusTreeParsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final NexusPhylogeniesParser parser = new NexusPhylogeniesParser();
            Phylogeny[] phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_1.nex", parser );
            if ( phylogenies.length != 1 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 25 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "" ) ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_2.nex", parser );
            if ( phylogenies.length != 1 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "name" ) ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_3.nex", parser );
            if ( phylogenies.length != 1 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "" ) ) {
                return false;
            }
            if ( phylogenies[ 0 ].isRooted() ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_4.nex", parser );
            if ( phylogenies.length != 18 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "tree 0" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getName().equals( "tree 1" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( phylogenies[ 2 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 3 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 4 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 5 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 6 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 7 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 8 ].getName().equals( "tree 8" ) ) {
                return false;
            }
            if ( phylogenies[ 8 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 8 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 9 ].getName().equals( "tree 9" ) ) {
                return false;
            }
            if ( !phylogenies[ 9 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 9 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 10 ].getName().equals( "tree 10" ) ) {
                return false;
            }
            if ( !phylogenies[ 10 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 10 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 11 ].getName().equals( "tree 11" ) ) {
                return false;
            }
            if ( phylogenies[ 11 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 11 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 12 ].getName().equals( "tree 12" ) ) {
                return false;
            }
            if ( !phylogenies[ 12 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 12 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 13 ].getName().equals( "tree 13" ) ) {
                return false;
            }
            if ( !phylogenies[ 13 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 13 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 14 ].getName().equals( "tree 14" ) ) {
                return false;
            }
            if ( !phylogenies[ 14 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 14 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 15 ].getName().equals( "tree 15" ) ) {
                return false;
            }
            if ( phylogenies[ 15 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 15 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 16 ].getName().equals( "tree 16" ) ) {
                return false;
            }
            if ( !phylogenies[ 16 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 16 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 17 ].getName().equals( "tree 17" ) ) {
                return false;
            }
            if ( phylogenies[ 17 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 17 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNexusTreeParsingIterating() {
        try {
            final NexusPhylogeniesParser p = new NexusPhylogeniesParser();
            p.setSource( Test.PATH_TO_TEST_DATA + "nexus_test_1.nex" );
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 25 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy != null ) {
                return false;
            }
            //
            p.reset();
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 25 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy != null ) {
                return false;
            }
            ////
            p.setSource( Test.PATH_TO_TEST_DATA + "nexus_test_2.nex" );
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phy.getName().equals( "name" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy != null ) {
                return false;
            }
            //
            p.reset();
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phy.getName().equals( "name" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy != null ) {
                return false;
            }
            ////
            p.setSource( Test.PATH_TO_TEST_DATA + "nexus_test_3.nex" );
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( phy.isRooted() ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy != null ) {
                return false;
            }
            //
            p.reset();
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy != null ) {
                return false;
            }
            ////
            p.setSource( Test.PATH_TO_TEST_DATA + "nexus_test_4_1.nex" );
            //            if ( phylogenies.length != 18 ) {
            //                return false;
            //            }
            //0
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 0" ) ) {
                return false;
            }
            //1
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 1" ) ) {
                return false;
            }
            //2
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( phy.isRooted() ) {
                return false;
            }
            //3
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //4
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 5 ) {
                System.out.println( phy.getNumberOfExternalNodes() );
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //5
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( phy.isRooted() ) {
                return false;
            }
            //6
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //7
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.toNewHampshire().equals( "((a,b),c);" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //8
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.toNewHampshire().equals( "((AA,BB),CC);" ) ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 8" ) ) {
                return false;
            }
            //9
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.toNewHampshire().equals( "((a,b),cc);" ) ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 9" ) ) {
                return false;
            }
            //10
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.toNewHampshire().equals( "((a,b),c);" ) ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 10" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //11
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.toNewHampshire().equals( "((1,2),3);" ) ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 11" ) ) {
                return false;
            }
            if ( phy.isRooted() ) {
                return false;
            }
            //12
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.toNewHampshire().equals( "((aa,bb),cc);" ) ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 12" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //13
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.toNewHampshire().equals( "((a,b),c);" ) ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 13" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //14
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                System.out.println( phy.getNumberOfExternalNodes() );
                return false;
            }
            if ( !phy
                    .toNewHampshire()
                    .equals( "(1:0.212481,8:0.297838,(9:0.222729,((6:0.201563,7:0.194547):0.282035,(4:1.146091,(3:1.008881,(10:0.384105,(2:0.235682,5:0.353432):0.32368):0.103875):0.41354):0.254687):0.095341):0.079254):0.0;" ) ) {
                System.out.println( phy.toNewHampshire() );
                return false;
            }
            if ( !phy.getName().equals( "tree 14" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //15
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                System.out.println( phy.getNumberOfExternalNodes() );
                return false;
            }
            if ( !phy
                    .toNewHampshire()
                    .equals( "(1:0.212481,8:0.297838,(9:0.222729,((6:0.201563,7:0.194547):0.282035,(4:1.146091,(3:1.008881,(10:0.384105,(2:0.235682,5:0.353432):0.32368):0.103875):0.41354):0.254687):0.095341):0.079254):0.0;" ) ) {
                System.out.println( phy.toNewHampshire() );
                return false;
            }
            if ( !phy.getName().equals( "tree 15" ) ) {
                return false;
            }
            if ( phy.isRooted() ) {
                return false;
            }
            //16
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                System.out.println( phy.getNumberOfExternalNodes() );
                return false;
            }
            if ( !phy
                    .toNewHampshire()
                    .equals( "(1:0.212481,8:0.297838,(9:0.222729,((6:0.201563,7:0.194547):0.282035,(4:1.146091,(3:1.008881,(10:0.384105,(2:0.235682,5:0.353432):0.32368):0.103875):0.41354):0.254687):0.095341):0.079254):0.0;" ) ) {
                System.out.println( phy.toNewHampshire() );
                return false;
            }
            if ( !phy.getName().equals( "tree 16" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //17
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                System.out.println( phy.getNumberOfExternalNodes() );
                return false;
            }
            if ( !phy
                    .toNewHampshire()
                    .equals( "(1:0.212481,8:0.297838,(9:0.222729,((6:0.201563,7:0.194547):0.282035,(4:1.146091,(3:1.008881,(10:0.384105,(2:0.235682,5:0.353432):0.32368):0.103875):0.41354):0.254687):0.095341):0.079254):0.0;" ) ) {
                System.out.println( phy.toNewHampshire() );
                return false;
            }
            if ( !phy.getName().equals( "tree 17" ) ) {
                return false;
            }
            if ( phy.isRooted() ) {
                return false;
            }
            //
            if ( p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy != null ) {
                return false;
            }
            p.reset();
            //0
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 0" ) ) {
                return false;
            }
            //1
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phy.getName().equals( "tree 1" ) ) {
                return false;
            }
            //2
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( phy.isRooted() ) {
                return false;
            }
            //3
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //4
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 5 ) {
                System.out.println( phy.getNumberOfExternalNodes() );
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( !phy.isRooted() ) {
                return false;
            }
            //5
            if ( !p.hasNext() ) {
                return false;
            }
            phy = p.next();
            if ( phy == null ) {
                return false;
            }
            if ( phy.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phy.getName().equals( "" ) ) {
                return false;
            }
            if ( phy.isRooted() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNexusTreeParsingTranslating() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final NexusPhylogeniesParser parser = new NexusPhylogeniesParser();
            Phylogeny[] phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_5.nex", parser );
            if ( phylogenies.length != 1 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "Tree0" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_6.nex", parser );
            if ( phylogenies.length != 3 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "Tree0" ) ) {
                return false;
            }
            if ( phylogenies[ 0 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getName().equals( "Tree1" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNextExternalNode().getName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            if ( phylogenies[ 2 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getName().equals( "Tree2" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNextExternalNode().getName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_7.nex", parser );
            if ( phylogenies.length != 3 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "Tree0" ) ) {
                return false;
            }
            if ( phylogenies[ 0 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getName().equals( "Tree1" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNextExternalNode().getName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            if ( phylogenies[ 2 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getName().equals( "Tree2" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNextExternalNode().getName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHParsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "(A,B1)", new NHXParser() )[ 0 ];
            if ( !p1.toNewHampshireX().equals( "(A,B1)" ) ) {
                return false;
            }
            final NHXParser nhxp = new NHXParser();
            nhxp.setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.NO );
            nhxp.setReplaceUnderscores( true );
            final Phylogeny uc0 = factory.create( "(A__A_,_B_B)", nhxp )[ 0 ];
            if ( !uc0.getRoot().getChildNode( 0 ).getName().equals( "A A " ) ) {
                return false;
            }
            if ( !uc0.getRoot().getChildNode( 1 ).getName().equals( " B B" ) ) {
                return false;
            }
            final Phylogeny p1b = factory
                    .create( "   \n  \t  \b   \r \f   ; (  \n  \t  \b   \r \f; A ;  \n  \t  \b   \r \f,  \n  \t  \b   \r \f; B ;   \n  \t  \b   \r \f 1  \n  \t  \b   \r \f ;  \n  \t  \b   \r \f );;;;; \n  \t  \b   \r \f;;;  \n  \t  \b   \r \f ",
                             new NHXParser() )[ 0 ];
            if ( !p1b.toNewHampshireX().equals( "(';A;',';B;1;')" ) ) {
                return false;
            }
            if ( !p1b.toNewHampshire().equals( "(';A;',';B;1;');" ) ) {
                return false;
            }
            final Phylogeny p2 = factory.create( new StringBuffer( "(A,B2)" ), new NHXParser() )[ 0 ];
            final Phylogeny p3 = factory.create( new char[] { '(', 'A', ',', 'B', '3', ')' }, new NHXParser() )[ 0 ];
            final Phylogeny p4 = factory.create( "(A,B4);", new NHXParser() )[ 0 ];
            final Phylogeny p5 = factory.create( new StringBuffer( "(A,B5);" ), new NHXParser() )[ 0 ];
            final Phylogeny[] p7 = factory.create( "(A,B7);(C,D7)", new NHXParser() );
            final Phylogeny[] p8 = factory.create( "(A,B8) (C,D8)", new NHXParser() );
            final Phylogeny[] p9 = factory.create( "(A,B9)\n(C,D9)", new NHXParser() );
            final Phylogeny[] p10 = factory.create( "(A,B10);(C,D10);", new NHXParser() );
            final Phylogeny[] p11 = factory.create( "(A,B11);(C,D11) (E,F11)\t(G,H11)", new NHXParser() );
            final Phylogeny[] p12 = factory.create( "(A,B12) (C,D12) (E,F12) (G,H12)", new NHXParser() );
            final Phylogeny[] p13 = factory.create( " ; (;A; , ; B ; 1  3 ; \n)\t ( \n ;"
                                                            + " C ; ,; D;13;);;;;;;(;E;,;F;13 ;) ; "
                                                            + "; ; ( \t\n\r\b; G ;, ;H ;1 3; )  ;  ;   ;",
                                                    new NHXParser() );
            if ( !p13[ 0 ].toNewHampshireX().equals( "(';A;',';B;13;')" ) ) {
                return false;
            }
            if ( !p13[ 1 ].toNewHampshireX().equals( "(';C;',';D;13;')" ) ) {
                return false;
            }
            if ( !p13[ 2 ].toNewHampshireX().equals( "(';E;',';F;13;')" ) ) {
                return false;
            }
            if ( !p13[ 3 ].toNewHampshireX().equals( "(';G;',';H;13;')" ) ) {
                return false;
            }
            final Phylogeny[] p14 = factory.create( "(A,B14)ab", new NHXParser() );
            final Phylogeny[] p15 = factory.create( "(A,B15)ab;", new NHXParser() );
            final String p16_S = "((A,B),C)";
            final Phylogeny[] p16 = factory.create( p16_S, new NHXParser() );
            if ( p16.length != 1 ) {
                return false;
            }
            if ( !p16[ 0 ].toNewHampshireX().equals( p16_S ) ) {
                return false;
            }
            final String p17_S = "(C,(A,B))";
            final Phylogeny[] p17 = factory.create( p17_S, new NHXParser() );
            if ( p17.length != 1 ) {
                return false;
            }
            if ( !p17[ 0 ].toNewHampshireX().equals( p17_S ) ) {
                return false;
            }
            final String p18_S = "((A,B),(C,D))";
            final Phylogeny[] p18 = factory.create( p18_S, new NHXParser() );
            if ( p18.length != 1 ) {
                return false;
            }
            if ( !p18[ 0 ].toNewHampshireX().equals( p18_S ) ) {
                return false;
            }
            final String p19_S = "(((A,B),C),D)";
            final Phylogeny[] p19 = factory.create( p19_S, new NHXParser() );
            if ( p19.length != 1 ) {
                return false;
            }
            if ( !p19[ 0 ].toNewHampshireX().equals( p19_S ) ) {
                return false;
            }
            final String p20_S = "(A,(B,(C,D)))";
            final Phylogeny[] p20 = factory.create( p20_S, new NHXParser() );
            if ( p20.length != 1 ) {
                return false;
            }
            if ( !p20[ 0 ].toNewHampshireX().equals( p20_S ) ) {
                return false;
            }
            final String p21_S = "(A,(B,(C,(D,E))))";
            final Phylogeny[] p21 = factory.create( p21_S, new NHXParser() );
            if ( p21.length != 1 ) {
                return false;
            }
            if ( !p21[ 0 ].toNewHampshireX().equals( p21_S ) ) {
                return false;
            }
            final String p22_S = "((((A,B),C),D),E)";
            final Phylogeny[] p22 = factory.create( p22_S, new NHXParser() );
            if ( p22.length != 1 ) {
                return false;
            }
            if ( !p22[ 0 ].toNewHampshireX().equals( p22_S ) ) {
                return false;
            }
            final String p23_S = "(A,(B,(C,(D,E)de)cde)bcde)abcde";
            final Phylogeny[] p23 = factory.create( p23_S, new NHXParser() );
            if ( p23.length != 1 ) {
                System.out.println( "xl=" + p23.length );
                System.exit( -1 );
                return false;
            }
            if ( !p23[ 0 ].toNewHampshireX().equals( p23_S ) ) {
                return false;
            }
            final String p24_S = "((((A,B)ab,C)abc,D)abcd,E)abcde";
            final Phylogeny[] p24 = factory.create( p24_S, new NHXParser() );
            if ( p24.length != 1 ) {
                return false;
            }
            if ( !p24[ 0 ].toNewHampshireX().equals( p24_S ) ) {
                return false;
            }
            final String p241_S1 = "(A,(B,(C,(D,E)de)cde)bcde)abcde";
            final String p241_S2 = "((((A,B)ab,C)abc,D)abcd,E)abcde";
            final Phylogeny[] p241 = factory.create( p241_S1 + p241_S2, new NHXParser() );
            if ( p241.length != 2 ) {
                return false;
            }
            if ( !p241[ 0 ].toNewHampshireX().equals( p241_S1 ) ) {
                return false;
            }
            if ( !p241[ 1 ].toNewHampshireX().equals( p241_S2 ) ) {
                return false;
            }
            final String p25_S = "((((((((((((((A,B)ab,C)abc,D)abcd,E)"
                    + "abcde,(B,(C,(D,E)de)cde)bcde)abcde,(B,((A,(B,(C,(D,"
                    + "E)de)cde)bcde)abcde,(D,E)de)cde)bcde)abcde,B)ab,C)"
                    + "abc,((((A,B)ab,C)abc,D)abcd,E)abcde)abcd,E)abcde,"
                    + "((((A,((((((((A,B)ab,C)abc,((((A,B)ab,C)abc,D)abcd,"
                    + "E)abcde)abcd,E)abcde,((((A,B)ab,C)abc,D)abcd,E)abcde)"
                    + "ab,C)abc,((((A,B)ab,C)abc,D)abcd,E)abcde)abcd,E)abcde"
                    + ")ab,C)abc,D)abcd,E)abcde)ab,C)abc,((((A,B)ab,C)abc,D)" + "abcd,E)abcde)abcd,E)abcde";
            final Phylogeny[] p25 = factory.create( p25_S, new NHXParser() );
            if ( !p25[ 0 ].toNewHampshireX().equals( p25_S ) ) {
                return false;
            }
            final String p26_S = "(A,B)ab";
            final Phylogeny[] p26 = factory.create( p26_S, new NHXParser() );
            if ( !p26[ 0 ].toNewHampshireX().equals( p26_S ) ) {
                return false;
            }
            final String p27_S = "((((A,B)ab,C)abc,D)abcd,E)abcde";
            final Phylogeny[] p27s = factory.create( p27_S, new NHXParser() );
            if ( p27s.length != 1 ) {
                System.out.println( "xxl=" + p27s.length );
                System.exit( -1 );
                return false;
            }
            if ( !p27s[ 0 ].toNewHampshireX().equals( p27_S ) ) {
                System.out.println( p27s[ 0 ].toNewHampshireX() );
                System.exit( -1 );
                return false;
            }
            final Phylogeny[] p27 = factory.create( new File( Test.PATH_TO_TEST_DATA + "phylogeny27.nhx" ),
                                                    new NHXParser() );
            if ( p27.length != 1 ) {
                System.out.println( "yl=" + p27.length );
                System.exit( -1 );
                return false;
            }
            if ( !p27[ 0 ].toNewHampshireX().equals( p27_S ) ) {
                System.out.println( p27[ 0 ].toNewHampshireX() );
                System.exit( -1 );
                return false;
            }
            final String p28_S1 = "((((A,B)ab,C)abc,D)abcd,E)abcde";
            final String p28_S2 = "(A,(B,(C,(D,E)de)cde)bcde)abcde";
            final String p28_S3 = "(A,B)ab";
            final String p28_S4 = "((((A,B),C),D),;E;)";
            final Phylogeny[] p28 = factory.create( new File( Test.PATH_TO_TEST_DATA + "phylogeny28.nhx" ),
                                                    new NHXParser() );
            if ( !p28[ 0 ].toNewHampshireX().equals( p28_S1 ) ) {
                return false;
            }
            if ( !p28[ 1 ].toNewHampshireX().equals( p28_S2 ) ) {
                return false;
            }
            if ( !p28[ 2 ].toNewHampshireX().equals( p28_S3 ) ) {
                return false;
            }
            if ( !p28[ 3 ].toNewHampshireX().equals( "((((A,B),C),D),';E;')" ) ) {
                return false;
            }
            if ( p28.length != 4 ) {
                return false;
            }
            final String p29_S = "((((A:0.01,B:0.684)ab:0.345,C:0.3451)abc:0.3451,D:1.5)abcd:0.134,E:0.32)abcde:0.1345";
            final Phylogeny[] p29 = factory.create( p29_S, new NHXParser() );
            if ( !p29[ 0 ].toNewHampshireX().equals( p29_S ) ) {
                return false;
            }
            final String p30_S = "((((A:0.01,B:0.02):0.93,C:0.04):0.05,D:1.4):0.06,E):0.72";
            final Phylogeny[] p30 = factory.create( p30_S, new NHXParser() );
            if ( !p30[ 0 ].toNewHampshireX().equals( p30_S ) ) {
                return false;
            }
            final String p32_S = " ;   ; 	\n  \t  \b   \f  \r  ;;;;;; ";
            final Phylogeny[] p32 = factory.create( p32_S, new NHXParser() );
            if ( ( p32.length != 0 ) ) {
                return false;
            }
            final String p33_S = "A";
            final Phylogeny[] p33 = factory.create( p33_S, new NHXParser() );
            if ( !p33[ 0 ].toNewHampshireX().equals( p33_S ) ) {
                return false;
            }
            final String p34_S = "B;";
            final Phylogeny[] p34 = factory.create( p34_S, new NHXParser() );
            if ( !p34[ 0 ].toNewHampshireX().equals( "B" ) ) {
                return false;
            }
            final String p35_S = "B:0.2";
            final Phylogeny[] p35 = factory.create( p35_S, new NHXParser() );
            if ( !p35[ 0 ].toNewHampshireX().equals( p35_S ) ) {
                return false;
            }
            final String p36_S = "(A)";
            final Phylogeny[] p36 = factory.create( p36_S, new NHXParser() );
            if ( !p36[ 0 ].toNewHampshireX().equals( p36_S ) ) {
                return false;
            }
            final String p37_S = "((A))";
            final Phylogeny[] p37 = factory.create( p37_S, new NHXParser() );
            if ( !p37[ 0 ].toNewHampshireX().equals( p37_S ) ) {
                return false;
            }
            final String p38_S = "(((((((A:0.2):0.2):0.3):0.4):0.5):0.6):0.7):0.8";
            final Phylogeny[] p38 = factory.create( p38_S, new NHXParser() );
            if ( !p38[ 0 ].toNewHampshireX().equals( p38_S ) ) {
                return false;
            }
            final String p39_S = "(((B,((((A:0.2):0.2):0.3):0.4):0.5):0.6):0.7):0.8";
            final Phylogeny[] p39 = factory.create( p39_S, new NHXParser() );
            if ( !p39[ 0 ].toNewHampshireX().equals( p39_S ) ) {
                return false;
            }
            final String p40_S = "(A,B,C)";
            final Phylogeny[] p40 = factory.create( p40_S, new NHXParser() );
            if ( !p40[ 0 ].toNewHampshireX().equals( p40_S ) ) {
                return false;
            }
            final String p41_S = "(A,B,C,D,E,F,G,H,I,J,K)";
            final Phylogeny[] p41 = factory.create( p41_S, new NHXParser() );
            if ( !p41[ 0 ].toNewHampshireX().equals( p41_S ) ) {
                return false;
            }
            final String p42_S = "(A,B,(X,Y,Z),D,E,F,G,H,I,J,K)";
            final Phylogeny[] p42 = factory.create( p42_S, new NHXParser() );
            if ( !p42[ 0 ].toNewHampshireX().equals( p42_S ) ) {
                return false;
            }
            final String p43_S = "(A,B,C,(AA,BB,CC,(CCC,DDD,EEE,(FFFF,GGGG)x)y,DD,EE,FF,GG,HH),D,E,(EE,FF),F,G,H,(((((5)4)3)2)1),I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,(XX,(YY)),Y,Z)";
            final Phylogeny[] p43 = factory.create( p43_S, new NHXParser() );
            if ( !p43[ 0 ].toNewHampshireX().equals( p43_S ) ) {
                return false;
            }
            final String p44_S = "(((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)))";
            final Phylogeny[] p44 = factory.create( p44_S, new NHXParser() );
            if ( !p44[ 0 ].toNewHampshireX().equals( p44_S ) ) {
                return false;
            }
            final String p45_S = "((((((((((A))))))))),(((((((((B))))))))),(((((((((C))))))))))";
            final Phylogeny[] p45 = factory.create( p45_S, new NHXParser() );
            if ( !p45[ 0 ].toNewHampshireX().equals( p45_S ) ) {
                return false;
            }
            final String p46_S = "";
            final Phylogeny[] p46 = factory.create( p46_S, new NHXParser() );
            if ( p46.length != 0 ) {
                return false;
            }
            final Phylogeny p47 = factory.create( new StringBuffer( "((A,B)ab:2[0.44],C)" ), new NHXParser() )[ 0 ];
            if ( !isEqual( 0.44, p47.getNode( "ab" ).getBranchData().getConfidence( 0 ).getValue() ) ) {
                return false;
            }
            final Phylogeny p48 = factory.create( new StringBuffer( "((A,B)ab:2[88],C)" ), new NHXParser() )[ 0 ];
            if ( !isEqual( 88, p48.getNode( "ab" ).getBranchData().getConfidence( 0 ).getValue() ) ) {
                return false;
            }
            final Phylogeny p49 = factory
                    .create( new StringBuffer( "((A,B)a[comment:a,b;(a)]b:2[0.44][comment(a,b,b);],C)" ),
                             new NHXParser() )[ 0 ];
            if ( !isEqual( 0.44, p49.getNode( "ab" ).getBranchData().getConfidence( 0 ).getValue() ) ) {
                return false;
            }
            final Phylogeny p50 = factory.create( new StringBuffer( "((\"A\",B)ab:2[88],C)" ), new NHXParser() )[ 0 ];
            if ( p50.getNode( "A" ) == null ) {
                return false;
            }
            if ( !p50.toNewHampshire( false, NH_CONVERSION_SUPPORT_VALUE_STYLE.IN_SQUARE_BRACKETS )
                    .equals( "((A,B)ab:2.0[88],C);" ) ) {
                return false;
            }
            if ( !p50.toNewHampshire( false, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE ).equals( "((A,B)ab:2.0,C);" ) ) {
                return false;
            }
            if ( !p50.toNewHampshire( false, NH_CONVERSION_SUPPORT_VALUE_STYLE.AS_INTERNAL_NODE_NAMES )
                    .equals( "((A,B)88:2.0,C);" ) ) {
                return false;
            }
            final Phylogeny p51 = factory.create( new StringBuffer( "((\"A(A\",B)ab:2[88],C)" ), new NHXParser() )[ 0 ];
            if ( p51.getNode( "A(A" ) == null ) {
                return false;
            }
            final Phylogeny p52 = factory.create( new StringBuffer( "(('A(A',B)ab:2[88],C)" ), new NHXParser() )[ 0 ];
            if ( p52.getNode( "A(A" ) == null ) {
                return false;
            }
            final Phylogeny p53 = factory
                    .create( new StringBuffer( "(('A(A',\"B (x (a' ,b) f(x);\"[com])[ment]ab:2[88],C)" ),
                             new NHXParser() )[ 0 ];
            if ( p53.getNode( "B (x (a' ,b) f(x);" ) == null ) {
                return false;
            }
            // 
            final Phylogeny p54 = factory.create( new StringBuffer( "((A,B):[88],C)" ), new NHXParser() )[ 0 ];
            if ( p54.getNode( "A" ) == null ) {
                return false;
            }
            if ( !p54.toNewHampshire( false, NH_CONVERSION_SUPPORT_VALUE_STYLE.IN_SQUARE_BRACKETS )
                    .equals( "((A,B)[88],C);" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHParsingIter() {
        try {
            final String p0_str = "(A,B);";
            final NHXParser p = new NHXParser();
            p.setSource( p0_str );
            if ( !p.hasNext() ) {
                return false;
            }
            final Phylogeny p0 = p.next();
            if ( !p0.toNewHampshire().equals( p0_str ) ) {
                System.out.println( p0.toNewHampshire() );
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            //
            final String p00_str = "(A,B)root;";
            p.setSource( p00_str );
            final Phylogeny p00 = p.next();
            if ( !p00.toNewHampshire().equals( p00_str ) ) {
                System.out.println( p00.toNewHampshire() );
                return false;
            }
            //
            final String p000_str = "A;";
            p.setSource( p000_str );
            final Phylogeny p000 = p.next();
            if ( !p000.toNewHampshire().equals( p000_str ) ) {
                System.out.println( p000.toNewHampshire() );
                return false;
            }
            //
            final String p0000_str = "A";
            p.setSource( p0000_str );
            final Phylogeny p0000 = p.next();
            if ( !p0000.toNewHampshire().equals( "A;" ) ) {
                System.out.println( p0000.toNewHampshire() );
                return false;
            }
            //
            p.setSource( "(A)" );
            final Phylogeny p00000 = p.next();
            if ( !p00000.toNewHampshire().equals( "(A);" ) ) {
                System.out.println( p00000.toNewHampshire() );
                return false;
            }
            //
            final String p1_str = "(A,B)(C,D)(E,F)(G,H)";
            p.setSource( p1_str );
            if ( !p.hasNext() ) {
                return false;
            }
            final Phylogeny p1_0 = p.next();
            if ( !p1_0.toNewHampshire().equals( "(A,B);" ) ) {
                System.out.println( p1_0.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            final Phylogeny p1_1 = p.next();
            if ( !p1_1.toNewHampshire().equals( "(C,D);" ) ) {
                System.out.println( "(C,D) != " + p1_1.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            final Phylogeny p1_2 = p.next();
            if ( !p1_2.toNewHampshire().equals( "(E,F);" ) ) {
                System.out.println( "(E,F) != " + p1_2.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            final Phylogeny p1_3 = p.next();
            if ( !p1_3.toNewHampshire().equals( "(G,H);" ) ) {
                System.out.println( "(G,H) != " + p1_3.toNewHampshire() );
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            //
            final String p2_str = "((1,2,3),B);(C,D) (E,F)root;(G,H); ;(X)";
            p.setSource( p2_str );
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p2_0 = p.next();
            if ( !p2_0.toNewHampshire().equals( "((1,2,3),B);" ) ) {
                System.out.println( p2_0.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p2_1 = p.next();
            if ( !p2_1.toNewHampshire().equals( "(C,D);" ) ) {
                System.out.println( "(C,D) != " + p2_1.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p2_2 = p.next();
            if ( !p2_2.toNewHampshire().equals( "(E,F)root;" ) ) {
                System.out.println( "(E,F)root != " + p2_2.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p2_3 = p.next();
            if ( !p2_3.toNewHampshire().equals( "(G,H);" ) ) {
                System.out.println( "(G,H) != " + p2_3.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p2_4 = p.next();
            if ( !p2_4.toNewHampshire().equals( "(X);" ) ) {
                System.out.println( "(X) != " + p2_4.toNewHampshire() );
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            ////
            p.reset();
            if ( !p.hasNext() ) {
                return false;
            }
            p2_0 = p.next();
            if ( !p2_0.toNewHampshire().equals( "((1,2,3),B);" ) ) {
                System.out.println( p2_0.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            p2_1 = p.next();
            if ( !p2_1.toNewHampshire().equals( "(C,D);" ) ) {
                System.out.println( "(C,D) != " + p2_1.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            p2_2 = p.next();
            if ( !p2_2.toNewHampshire().equals( "(E,F)root;" ) ) {
                System.out.println( "(E,F)root != " + p2_2.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            p2_3 = p.next();
            if ( !p2_3.toNewHampshire().equals( "(G,H);" ) ) {
                System.out.println( "(G,H) != " + p2_3.toNewHampshire() );
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            p2_4 = p.next();
            if ( !p2_4.toNewHampshire().equals( "(X);" ) ) {
                System.out.println( "(X) != " + p2_4.toNewHampshire() );
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            //
            final String p3_str = "((A,B),C)abc";
            p.setSource( p3_str );
            if ( !p.hasNext() ) {
                return false;
            }
            final Phylogeny p3_0 = p.next();
            if ( !p3_0.toNewHampshire().equals( "((A,B),C)abc;" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            //
            final String p4_str = "((A,B)ab,C)abc";
            p.setSource( p4_str );
            if ( !p.hasNext() ) {
                return false;
            }
            final Phylogeny p4_0 = p.next();
            if ( !p4_0.toNewHampshire().equals( "((A,B)ab,C)abc;" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            //
            final String p5_str = "(((A,B)ab,C)abc,D)abcd";
            p.setSource( p5_str );
            if ( !p.hasNext() ) {
                return false;
            }
            final Phylogeny p5_0 = p.next();
            if ( !p5_0.toNewHampshire().equals( "(((A,B)ab,C)abc,D)abcd;" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            //
            final String p6_str = "(A,(B,(C,(D,E)de)cde)bcde)abcde";
            p.setSource( p6_str );
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p6_0 = p.next();
            if ( !p6_0.toNewHampshire().equals( "(A,(B,(C,(D,E)de)cde)bcde)abcde;" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            p.reset();
            if ( !p.hasNext() ) {
                return false;
            }
            p6_0 = p.next();
            if ( !p6_0.toNewHampshire().equals( "(A,(B,(C,(D,E)de)cde)bcde)abcde;" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            //
            final String p7_str = "((((A,B)ab,C)abc,D)abcd,E)abcde";
            p.setSource( p7_str );
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p7_0 = p.next();
            if ( !p7_0.toNewHampshire().equals( "((((A,B)ab,C)abc,D)abcd,E)abcde;" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            p.reset();
            if ( !p.hasNext() ) {
                return false;
            }
            p7_0 = p.next();
            if ( !p7_0.toNewHampshire().equals( "((((A,B)ab,C)abc,D)abcd,E)abcde;" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            //
            final String p8_str = "((((A,B)ab,C)abc,D)abcd,E)abcde ((((a,b)ab,c)abc,d)abcd,e)abcde";
            p.setSource( p8_str );
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p8_0 = p.next();
            if ( !p8_0.toNewHampshire().equals( "((((A,B)ab,C)abc,D)abcd,E)abcde;" ) ) {
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p8_1 = p.next();
            if ( !p8_1.toNewHampshire().equals( "((((a,b)ab,c)abc,d)abcd,e)abcde;" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            p.reset();
            if ( !p.hasNext() ) {
                return false;
            }
            p8_0 = p.next();
            if ( !p8_0.toNewHampshire().equals( "((((A,B)ab,C)abc,D)abcd,E)abcde;" ) ) {
                return false;
            }
            if ( !p.hasNext() ) {
                return false;
            }
            p8_1 = p.next();
            if ( !p8_1.toNewHampshire().equals( "((((a,b)ab,c)abc,d)abcd,e)abcde;" ) ) {
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            p.reset();
            //
            p.setSource( "" );
            if ( p.hasNext() ) {
                return false;
            }
            //
            p.setSource( new File( Test.PATH_TO_TEST_DATA + "phylogeny27.nhx" ) );
            if ( !p.hasNext() ) {
                return false;
            }
            Phylogeny p_27 = p.next();
            if ( !p_27.toNewHampshireX().equals( "((((A,B)ab,C)abc,D)abcd,E)abcde" ) ) {
                System.out.println( p_27.toNewHampshireX() );
                System.exit( -1 );
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
            p.reset();
            if ( !p.hasNext() ) {
                return false;
            }
            p_27 = p.next();
            if ( !p_27.toNewHampshireX().equals( "((((A,B)ab,C)abc,D)abcd,E)abcde" ) ) {
                System.out.println( p_27.toNewHampshireX() );
                System.exit( -1 );
                return false;
            }
            if ( p.hasNext() ) {
                return false;
            }
            if ( p.next() != null ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHXconversion() {
        try {
            final PhylogenyNode n1 = new PhylogenyNode();
            final PhylogenyNode n2 = PhylogenyNode.createInstanceFromNhxString( "" );
            final PhylogenyNode n3 = PhylogenyNode.createInstanceFromNhxString( "n3" );
            final PhylogenyNode n4 = PhylogenyNode.createInstanceFromNhxString( "n4:0.01" );
            final PhylogenyNode n5 = PhylogenyNode
                    .createInstanceFromNhxString( "n5:0.1[&&NHX:S=Ecoli:E=1.1.1.1:D=Y:Co=Y:B=56:T=1]" );
            final PhylogenyNode n6 = PhylogenyNode
                    .createInstanceFromNhxString( "n6:0.000001[&&NHX:S=Ecoli:E=1.1.1.1:D=N:Co=N:B=100:T=1]" );
            if ( !n1.toNewHampshireX().equals( "" ) ) {
                return false;
            }
            if ( !n2.toNewHampshireX().equals( "" ) ) {
                return false;
            }
            if ( !n3.toNewHampshireX().equals( "n3" ) ) {
                return false;
            }
            if ( !n4.toNewHampshireX().equals( "n4:0.01" ) ) {
                return false;
            }
            if ( !n5.toNewHampshireX().equals( "n5:0.1[&&NHX:T=1:S=Ecoli:D=Y:B=56]" ) ) {
                return false;
            }
            if ( !n6.toNewHampshireX().equals( "n6:1.0E-6[&&NHX:T=1:S=Ecoli:D=N:B=100]" ) ) {
                System.out.println( n6.toNewHampshireX() );
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHXNodeParsing() {
        try {
            final PhylogenyNode n1 = new PhylogenyNode();
            final PhylogenyNode n2 = PhylogenyNode.createInstanceFromNhxString( "" );
            final PhylogenyNode n3 = PhylogenyNode.createInstanceFromNhxString( "n3" );
            final PhylogenyNode n4 = PhylogenyNode.createInstanceFromNhxString( "n4:0.01" );
            final PhylogenyNode n5 = PhylogenyNode
                    .createInstanceFromNhxString( "n5:0.1[&&NHX:S=Ecoli:E=1.1.1.1:D=Y:B=56:T=1:On=22:SOn=33:SNn=44:W=2:C=10.20.30:XN=S=tag1=value1=unit1:XN=S=tag3=value3=unit3]" );
            if ( !n3.getName().equals( "n3" ) ) {
                return false;
            }
            if ( n3.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                return false;
            }
            if ( n3.isDuplication() ) {
                return false;
            }
            if ( n3.isHasAssignedEvent() ) {
                return false;
            }
            if ( PhylogenyMethods.getBranchWidthValue( n3 ) != BranchWidth.BRANCH_WIDTH_DEFAULT_VALUE ) {
                return false;
            }
            if ( !n4.getName().equals( "n4" ) ) {
                return false;
            }
            if ( n4.getDistanceToParent() != 0.01 ) {
                return false;
            }
            if ( !n5.getName().equals( "n5" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( n5 ) != 56 ) {
                return false;
            }
            if ( n5.getDistanceToParent() != 0.1 ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n5 ).equals( "Ecoli" ) ) {
                return false;
            }
            if ( !n5.isDuplication() ) {
                return false;
            }
            if ( !n5.isHasAssignedEvent() ) {
                return false;
            }
            final PhylogenyNode n8 = PhylogenyNode
                    .createInstanceFromNhxString( "ABCD_ECOLI/1-2:0.01",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n8.getName().equals( "ABCD_ECOLI/1-2" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n8 ).equals( "ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n9 = PhylogenyNode
                    .createInstanceFromNhxString( "ABCD_ECOLI/1-12:0.01",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n9.getName().equals( "ABCD_ECOLI/1-12" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n9 ).equals( "ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n10 = PhylogenyNode
                    .createInstanceFromNhxString( "n10.ECOLI", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n10.getName().equals( "n10.ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n20 = PhylogenyNode
                    .createInstanceFromNhxString( "ABCD_ECOLI/1-2", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n20.getName().equals( "ABCD_ECOLI/1-2" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n20 ).equals( "ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n20x = PhylogenyNode
                    .createInstanceFromNhxString( "N20_ECOL1/1-2", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n20x.getName().equals( "N20_ECOL1/1-2" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n20x ).equals( "ECOL1" ) ) {
                return false;
            }
            final PhylogenyNode n20xx = PhylogenyNode
                    .createInstanceFromNhxString( "N20_eCOL1/1-2", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n20xx.getName().equals( "N20_eCOL1/1-2" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n20xx ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n20xxx = PhylogenyNode
                    .createInstanceFromNhxString( "n20_ecoli/1-2", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n20xxx.getName().equals( "n20_ecoli/1-2" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n20xxx ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n20xxxx = PhylogenyNode
                    .createInstanceFromNhxString( "n20_Ecoli/1-2", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n20xxxx.getName().equals( "n20_Ecoli/1-2" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n20xxxx ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n21 = PhylogenyNode
                    .createInstanceFromNhxString( "N21_PIG", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n21.getName().equals( "N21_PIG" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n21 ).equals( "PIG" ) ) {
                return false;
            }
            final PhylogenyNode n21x = PhylogenyNode
                    .createInstanceFromNhxString( "n21_PIG", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n21x.getName().equals( "n21_PIG" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n21x ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n22 = PhylogenyNode
                    .createInstanceFromNhxString( "n22/PIG", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n22.getName().equals( "n22/PIG" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n22 ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n23 = PhylogenyNode
                    .createInstanceFromNhxString( "n23/PIG_1", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n23.getName().equals( "n23/PIG_1" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n23 ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode a = PhylogenyNode
                    .createInstanceFromNhxString( "ABCD_ECOLI/1-2", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !a.getName().equals( "ABCD_ECOLI/1-2" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( a ).equals( "ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode c1 = PhylogenyNode
                    .createInstanceFromNhxString( "n10_BOVIN/1000-2000",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !c1.getName().equals( "n10_BOVIN/1000-2000" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( c1 ).equals( "BOVIN" ) ) {
                return false;
            }
            final PhylogenyNode c2 = PhylogenyNode
                    .createInstanceFromNhxString( "N10_Bovin_1/1000-2000",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !c2.getName().equals( "N10_Bovin_1/1000-2000" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( c2 ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode e3 = PhylogenyNode
                    .createInstanceFromNhxString( "n10_RAT~", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !e3.getName().equals( "n10_RAT~" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( e3 ).equals( "RAT" ) ) {
                return false;
            }
            final PhylogenyNode n11 = PhylogenyNode
                    .createInstanceFromNhxString( "N111111_ECOLI/1-2:0.4",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n11.getName().equals( "N111111_ECOLI/1-2" ) ) {
                return false;
            }
            if ( n11.getDistanceToParent() != 0.4 ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n11 ).equals( "ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n12 = PhylogenyNode
                    .createInstanceFromNhxString( "N111111-ECOLI---/jdj:0.4",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n12.getName().equals( "N111111-ECOLI---/jdj" ) ) {
                return false;
            }
            if ( n12.getDistanceToParent() != 0.4 ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n12 ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode o = PhylogenyNode
                    .createInstanceFromNhxString( "ABCD_MOUSE", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !o.getName().equals( "ABCD_MOUSE" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( o ).equals( "MOUSE" ) ) {
                return false;
            }
            if ( n1.getName().compareTo( "" ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( n1 ) != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                return false;
            }
            if ( n1.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                return false;
            }
            if ( n2.getName().compareTo( "" ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( n2 ) != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                return false;
            }
            if ( n2.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                return false;
            }
            final PhylogenyNode n00 = PhylogenyNode
                    .createInstanceFromNhxString( "n7:0.000001[&&NHX:GN=gene_name:AC=accession123:S=Ecoli:D=N:Co=N:B=100:T=1]" );
            if ( !n00.getNodeData().getSequence().getName().equals( "gene_name" ) ) {
                return false;
            }
            if ( !n00.getNodeData().getSequence().getAccession().getValue().equals( "accession123" ) ) {
                return false;
            }
            final PhylogenyNode nx = PhylogenyNode.createInstanceFromNhxString( "n5:0.1[&&NHX:S=Ecoli:GN=gene_1]" );
            if ( !nx.getNodeData().getSequence().getName().equals( "gene_1" ) ) {
                return false;
            }
            final PhylogenyNode n13 = PhylogenyNode
                    .createInstanceFromNhxString( "blah_12345/1-2", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n13.getName().equals( "blah_12345/1-2" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n13 ).equals( "12345" ) ) {
                return false;
            }
            if ( !n13.getNodeData().getTaxonomy().getIdentifier().getValue().equals( "12345" ) ) {
                return false;
            }
            if ( !n13.getNodeData().getTaxonomy().getIdentifier().getProvider().equals( "uniprot" ) ) {
                return false;
            }
            final PhylogenyNode n14 = PhylogenyNode
                    .createInstanceFromNhxString( "BLA1_9QX45/1-2", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n14.getName().equals( "BLA1_9QX45/1-2" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n14 ).equals( "9QX45" ) ) {
                return false;
            }
            final PhylogenyNode n15 = PhylogenyNode
                    .createInstanceFromNhxString( "something_wicked[123]",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n15.getName().equals( "something_wicked" ) ) {
                return false;
            }
            if ( n15.getBranchData().getNumberOfConfidences() != 1 ) {
                return false;
            }
            if ( !isEqual( n15.getBranchData().getConfidence( 0 ).getValue(), 123 ) ) {
                return false;
            }
            final PhylogenyNode n16 = PhylogenyNode
                    .createInstanceFromNhxString( "something_wicked2[9]",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n16.getName().equals( "something_wicked2" ) ) {
                return false;
            }
            if ( n16.getBranchData().getNumberOfConfidences() != 1 ) {
                return false;
            }
            if ( !isEqual( n16.getBranchData().getConfidence( 0 ).getValue(), 9 ) ) {
                return false;
            }
            final PhylogenyNode n17 = PhylogenyNode
                    .createInstanceFromNhxString( "something_wicked3[a]",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !n17.getName().equals( "something_wicked3" ) ) {
                return false;
            }
            if ( n17.getBranchData().getNumberOfConfidences() != 0 ) {
                return false;
            }
            final PhylogenyNode n18 = PhylogenyNode
                    .createInstanceFromNhxString( ":0.5[91]", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( !isEqual( n18.getDistanceToParent(), 0.5 ) ) {
                return false;
            }
            if ( n18.getBranchData().getNumberOfConfidences() != 1 ) {
                return false;
            }
            if ( !isEqual( n18.getBranchData().getConfidence( 0 ).getValue(), 91 ) ) {
                return false;
            }
            final PhylogenyNode n19 = PhylogenyNode
                    .createInstanceFromNhxString( "blah_1-roejojoej", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n19.getNodeData().getTaxonomy().getIdentifier().getValue().equals( "1" ) ) {
                return false;
            }
            if ( !n19.getNodeData().getTaxonomy().getIdentifier().getProvider().equals( "uniprot" ) ) {
                return false;
            }
            final PhylogenyNode n30 = PhylogenyNode
                    .createInstanceFromNhxString( "blah_1234567-roejojoej",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n30.getNodeData().getTaxonomy().getIdentifier().getValue().equals( "1234567" ) ) {
                return false;
            }
            if ( !n30.getNodeData().getTaxonomy().getIdentifier().getProvider().equals( "uniprot" ) ) {
                return false;
            }
            final PhylogenyNode n31 = PhylogenyNode
                    .createInstanceFromNhxString( "blah_12345678-roejojoej",
                                                  NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n31.getNodeData().isHasTaxonomy() ) {
                return false;
            }
            final PhylogenyNode n32 = PhylogenyNode
                    .createInstanceFromNhxString( "sd_12345678", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n32.getNodeData().isHasTaxonomy() ) {
                return false;
            }
            final PhylogenyNode n40 = PhylogenyNode
                    .createInstanceFromNhxString( "bcl2_12345", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n40.getNodeData().getTaxonomy().getIdentifier().getValue().equals( "12345" ) ) {
                return false;
            }
            final PhylogenyNode n41 = PhylogenyNode
                    .createInstanceFromNhxString( "12345", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n41.getNodeData().isHasTaxonomy() ) {
                return false;
            }
            final PhylogenyNode n42 = PhylogenyNode
                    .createInstanceFromNhxString( "12345", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
            if ( n42.getNodeData().isHasTaxonomy() ) {
                return false;
            }
            final PhylogenyNode n43 = PhylogenyNode.createInstanceFromNhxString( "12345",
                                                                                 NHXParser.TAXONOMY_EXTRACTION.NO );
            if ( n43.getNodeData().isHasTaxonomy() ) {
                return false;
            }
            final PhylogenyNode n44 = PhylogenyNode
                    .createInstanceFromNhxString( "12345~1-2", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n44.getNodeData().isHasTaxonomy() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHXParsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "(A     [&&NHX:S=a_species],B1[&&NHX:S=b_species])", new NHXParser() )[ 0 ];
            if ( !p1.toNewHampshireX().equals( "(A[&&NHX:S=a_species],B1[&&NHX:S=b_species])" ) ) {
                return false;
            }
            final String p2_S = "(((((((A:0.2[&&NHX:S=qwerty]):0.2[&&NHX:S=uiop]):0.3[&&NHX:S=asdf]):0.4[&&NHX:S=zxc]):0.5[&&NHX:S=a]):0.6[&&NHX:S=asd]):0.7[&&NHX:S=za]):0.8[&&NHX:S=zaq]";
            final Phylogeny[] p2 = factory.create( p2_S, new NHXParser() );
            if ( !p2[ 0 ].toNewHampshireX().equals( p2_S ) ) {
                return false;
            }
            final String p2b_S = "(((((((A:0.2[&NHX:S=qw,erty]):0.2[&:S=u(io)p]):0.3[&NHX:S=asdf]):0.4[S=zxc]):0.5[]):0.6[&&NH:S=asd]):0.7[&&HX:S=za]):0.8[&&:S=zaq]";
            final Phylogeny[] p2b = factory.create( p2b_S, new NHXParser() );
            if ( !p2b[ 0 ].toNewHampshireX().equals( "(((((((A:0.2):0.2):0.3):0.4):0.5):0.6):0.7):0.8" ) ) {
                return false;
            }
            final Phylogeny[] p3 = factory
                    .create( "[  comment&&NHX,())))](((((((A:0.2[&&NHX:S=qwerty]):0.2[&&NHX:S=uiop]):0.3[&&NHX:S=asdf]):0.4[&&NHX:S=zxc]):0.5[&&NHX:S=a]):0.6[&&NHX:S=asd]):0.7[&&NHX:S=za]):0.8[&&NHX:S=zaq]",
                             new NHXParser() );
            if ( !p3[ 0 ].toNewHampshireX().equals( p2_S ) ) {
                return false;
            }
            final Phylogeny[] p4 = factory
                    .create( "(((((((A:0.2[&&NHX:S=qwerty]):0.2[&&NHX:S=uiop]):0.3[&&NHX:S=asdf]):0.4[&&NHX:S=zxc]):0.5[&&NHX:S=a]):0.6[&&NHX:S=asd]):0.7[&&NHX:S=za]):0.8[&&NHX:S=zaq][comment(]",
                             new NHXParser() );
            if ( !p4[ 0 ].toNewHampshireX().equals( p2_S ) ) {
                return false;
            }
            final Phylogeny[] p5 = factory
                    .create( "[]  (  [][ ][   ]  ([((( &&NHXcomment only![[[[[[]([]((((A:0.2[&&NHX:S=q[comment )))]werty][,,,,))]):0.2[&&NHX:S=uiop]):0.3[&&NHX:S=a[comment,,))]sdf])[comment(((]:0.4[&&NHX:S=zxc][comment(((][comment(((]):0.5[&&NHX:S=a]):0.6[&&NHX:S=a[comment(((]sd]):0.7[&&NHX:S=za]):0.8[&&NHX:S=zaq][comment(((]",
                             new NHXParser() );
            if ( !p5[ 0 ].toNewHampshireX().equals( p2_S ) ) {
                return false;
            }
            final String p6_S_C = "(A[][][][1][22][333][4444][55555][666666][&&NHX:S=Aspecies],B[))],C,(AA,BB,CC,(CCC,DDD,EEE,[comment](FFFF,GGGG)x)y,D[comment]D,EE,FF,GG,HH),D,E,(EE,FF),F,G,H,(((((5)4)3)2)1),I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,(XX,(YY)),Y,Z)";
            final String p6_S_WO_C = "(A[&&NHX:S=Aspecies],B,C,(AA,BB,CC,(CCC,DDD,EEE,(FFFF,GGGG)x)y,DD,EE,FF,GG,HH),D,E,(EE,FF),F,G,H,(((((5)4)3)2)1),I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,(XX,(YY)),Y,Z)";
            final Phylogeny[] p6 = factory.create( p6_S_C, new NHXParser() );
            if ( !p6[ 0 ].toNewHampshireX().equals( p6_S_WO_C ) ) {
                return false;
            }
            final String p7_S_C = "(((A [&&NHX:S=species_a], B [&&NHX:S=Vstorri] , C   , D),(A,B,C,D[comment])[],[c][]([xxx]A[comment],[comment]B[comment][comment],[comment][comment]C[comment][comment],[comment][comment]D[comment][comment])[comment][comment],[comment]   [comment](A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C[comment][comment][comment][comment][comment]    [comment],D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),[comment][comment]((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)))";
            final String p7_S_WO_C = "(((A[&&NHX:S=species_a],B[&&NHX:S=Vstorri],C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)))";
            final Phylogeny[] p7 = factory.create( p7_S_C, new NHXParser() );
            if ( !p7[ 0 ].toNewHampshireX().equals( p7_S_WO_C ) ) {
                return false;
            }
            final String p8_S_C = "[cmt](((([]([))))))](((((A[&&NHX:S= [a comment] a])))))))[too many comments!:)])),(((((((((B[&&NHX[ a comment in a bad place]:S   =b])))))[] []   )))),(((((((((C[&&NHX:S=c])   ))[,,, ])))))))";
            final String p8_S_WO_C = "((((((((((A[&&NHX:S=a]))))))))),(((((((((B[&&NHX:S=b]))))))))),(((((((((C[&&NHX:S=c]))))))))))";
            final Phylogeny[] p8 = factory.create( p8_S_C, new NHXParser() );
            if ( !p8[ 0 ].toNewHampshireX().equals( p8_S_WO_C ) ) {
                return false;
            }
            final Phylogeny p9 = factory.create( "((A:0.2,B:0.3):0.5[91],C:0.1)root:0.1[100]", new NHXParser() )[ 0 ];
            if ( !p9.toNewHampshireX().equals( "((A:0.2,B:0.3):0.5[&&NHX:B=91],C:0.1)root:0.1[&&NHX:B=100]" ) ) {
                return false;
            }
            final Phylogeny p10 = factory
                    .create( " [79]   ( (A [co mment] :0 .2[comment],B:0.3[com])[com ment]: 0. 5 \t[ 9 1 ][ comment],C: 0.1)[comment]root:0.1[100] [comment]",
                             new NHXParser() )[ 0 ];
            if ( !p10.toNewHampshireX().equals( "((A:0.2,B:0.3):0.5[&&NHX:B=91],C:0.1)root:0.1[&&NHX:B=100]" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHXParsingMB() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "(1[&prob=0.9500000000000000e+00,prob_stddev=0.1100000000000000e+00,"
                    + "prob_range={1.000000000000000e+00,1.000000000000000e+00},prob(percent)=\"100\","
                    + "prob+-sd=\"100+-0\"]:4.129000000000000e-02[&length_mean=4.153987461671767e-02,"
                    + "length_median=4.129000000000000e-02,length_95%HPD={3.217800000000000e-02,"
                    + "5.026800000000000e-02}],2[&prob=0.810000000000000e+00,prob_stddev=0.000000000000000e+00,"
                    + "prob_range={1.000000000000000e+00,1.000000000000000e+00},prob(percent)=\"100\","
                    + "prob+-sd=\"100+-0\"]:6.375699999999999e-02[&length_mean=6.395210411945065e-02,"
                    + "length_median=6.375699999999999e-02,length_95%HPD={5.388600000000000e-02,"
                    + "7.369400000000000e-02}])", new NHXParser() )[ 0 ];
            if ( !isEqual( p1.getNode( "1" ).getDistanceToParent(), 4.129e-02 ) ) {
                return false;
            }
            if ( !isEqual( p1.getNode( "1" ).getBranchData().getConfidence( 0 ).getValue(), 0.9500000000000000e+00 ) ) {
                return false;
            }
            if ( !isEqual( p1.getNode( "1" ).getBranchData().getConfidence( 0 ).getStandardDeviation(),
                           0.1100000000000000e+00 ) ) {
                return false;
            }
            if ( !isEqual( p1.getNode( "2" ).getDistanceToParent(), 6.375699999999999e-02 ) ) {
                return false;
            }
            if ( !isEqual( p1.getNode( "2" ).getBranchData().getConfidence( 0 ).getValue(), 0.810000000000000e+00 ) ) {
                return false;
            }
            final Phylogeny p2 = factory
                    .create( "(1[something_else(?)s,prob=0.9500000000000000e+00{}(((,p)rob_stddev=0.110000000000e+00,"
                                     + "prob_range={1.000000000000000e+00,1.000000000000000e+00},prob(percent)=\"100\","
                                     + "prob+-sd=\"100+-0\"]:4.129000000000000e-02[&length_mean=4.153987461671767e-02,"
                                     + "length_median=4.129000000000000e-02,length_95%HPD={3.217800000000000e-02,"
                                     + "5.026800000000000e-02}],2[&prob=0.810000000000000e+00,prob_stddev=0.000000000000000e+00,"
                                     + "prob_range={1.000000000000000e+00,1.000000000000000e+00},prob(percent)=\"100\","
                                     + "prob+-sd=\"100+-0\"]:6.375699999999999e-02[&length_mean=6.395210411945065e-02,"
                                     + "length_median=6.375699999999999e-02,length_95%HPD={5.388600000000000e-02,"
                                     + "7.369400000000000e-02}])",
                             new NHXParser() )[ 0 ];
            if ( p2.getNode( "1" ) == null ) {
                return false;
            }
            if ( p2.getNode( "2" ) == null ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            System.exit( -1 );
            return false;
        }
        return true;
    }

    private static boolean testNHXParsingQuotes() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final NHXParser p = new NHXParser();
            final Phylogeny[] phylogenies_0 = factory.create( new File( Test.PATH_TO_TEST_DATA + "quotes.nhx" ), p );
            if ( phylogenies_0.length != 5 ) {
                return false;
            }
            final Phylogeny phy = phylogenies_0[ 4 ];
            if ( phy.getNumberOfExternalNodes() != 7 ) {
                return false;
            }
            if ( phy.getNodes( "a name in double quotes from tree ((a,b),c)" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "charles darwin 'origin of species'" ).size() != 1 ) {
                return false;
            }
            if ( !phy.getNodes( "charles darwin 'origin of species'" ).get( 0 ).getNodeData().getTaxonomy()
                    .getScientificName().equals( "hsapiens" ) ) {
                return false;
            }
            if ( phy.getNodes( "shouldbetogether single quotes" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "'single quotes' inside double quotes" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "double quotes inside single quotes" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "noquotes" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "A   (  B    C '" ).size() != 1 ) {
                return false;
            }
            final NHXParser p1p = new NHXParser();
            p1p.setIgnoreQuotes( true );
            final Phylogeny p1 = factory.create( "(\"A\",'B1')", p1p )[ 0 ];
            if ( !p1.toNewHampshire().equals( "(A,B1);" ) ) {
                return false;
            }
            final NHXParser p2p = new NHXParser();
            p1p.setIgnoreQuotes( false );
            final Phylogeny p2 = factory.create( "(\"A\",'B1')", p2p )[ 0 ];
            if ( !p2.toNewHampshire().equals( "(A,B1);" ) ) {
                return false;
            }
            final NHXParser p3p = new NHXParser();
            p3p.setIgnoreQuotes( false );
            final Phylogeny p3 = factory.create( "(\"A)\",'B1')", p3p )[ 0 ];
            if ( !p3.toNewHampshire().equals( "('A)',B1);" ) ) {
                return false;
            }
            final NHXParser p4p = new NHXParser();
            p4p.setIgnoreQuotes( false );
            final Phylogeny p4 = factory.create( "(\"A)\",'B(),; x')", p4p )[ 0 ];
            if ( !p4.toNewHampshire().equals( "('A)','B(),; x');" ) ) {
                return false;
            }
            final Phylogeny p10 = factory
                    .create( " [79]   ( (\"A \n\tB \" [co mment] :0 .2[comment],'B':0.3[com])[com ment]: 0. 5 \t[ 9 1 ][ comment],'C (or D?\\//;,))': 0.1)[comment]'\nroot is here (cool,  was! ) ':0.1[100] [comment]",
                             new NHXParser() )[ 0 ];
            final String p10_clean_str = "(('A B':0.2,B:0.3):0.5[&&NHX:B=91],'C (or D?\\//;,))':0.1)'root is here (cool,  was! )':0.1[&&NHX:B=100]";
            if ( !p10.toNewHampshireX().equals( p10_clean_str ) ) {
                return false;
            }
            final Phylogeny p11 = factory.create( p10.toNewHampshireX(), new NHXParser() )[ 0 ];
            if ( !p11.toNewHampshireX().equals( p10_clean_str ) ) {
                return false;
            }
            //
            final Phylogeny p12 = factory
                    .create( " [79]   ( (\"A \n\tB \" [[][] :0 .2[comment][\t&\t&\n N\tH\tX:S=mo\tnkey !],'\tB\t\b\t\n\f\rB B ':0.0\b3[])\t[com ment]: 0. 5 \t[ 9 1 ][ \ncomment],'C\t (or D?\\//;,))': 0.\b1)[comment]'\nroot \tis here (cool, \b\t\n\f\r was! ) ':0.1[100] [comment]",
                             new NHXParser() )[ 0 ];
            final String p12_clean_str = "(('A B':0.2[&&NHX:S=monkey!],'BB B':0.03):0.5[&&NHX:B=91],'C (or D?\\//;,))':0.1)'root is here (cool,  was! )':0.1[&&NHX:B=100]";
            if ( !p12.toNewHampshireX().equals( p12_clean_str ) ) {
                return false;
            }
            final Phylogeny p13 = factory.create( p12.toNewHampshireX(), new NHXParser() )[ 0 ];
            if ( !p13.toNewHampshireX().equals( p12_clean_str ) ) {
                return false;
            }
            final String p12_clean_str_nh = "(('A B':0.2,'BB B':0.03):0.5,'C (or D?\\//;,))':0.1)'root is here (cool,  was! )':0.1;";
            if ( !p13.toNewHampshire().equals( p12_clean_str_nh ) ) {
                return false;
            }
            final Phylogeny p14 = factory.create( p13.toNewHampshire(), new NHXParser() )[ 0 ];
            if ( !p14.toNewHampshire().equals( p12_clean_str_nh ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNodeRemoval() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((a)b)", new NHXParser() )[ 0 ];
            PhylogenyMethods.removeNode( t0.getNode( "b" ), t0 );
            if ( !t0.toNewHampshire().equals( "(a);" ) ) {
                return false;
            }
            final Phylogeny t1 = factory.create( "((a:2)b:4)", new NHXParser() )[ 0 ];
            PhylogenyMethods.removeNode( t1.getNode( "b" ), t1 );
            if ( !t1.toNewHampshire().equals( "(a:6.0);" ) ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "((a,b),c)", new NHXParser() )[ 0 ];
            PhylogenyMethods.removeNode( t2.getNode( "b" ), t2 );
            if ( !t2.toNewHampshire().equals( "((a),c);" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPhylogenyBranch() {
        try {
            final PhylogenyNode a1 = PhylogenyNode.createInstanceFromNhxString( "a" );
            final PhylogenyNode b1 = PhylogenyNode.createInstanceFromNhxString( "b" );
            final PhylogenyBranch a1b1 = new PhylogenyBranch( a1, b1 );
            final PhylogenyBranch b1a1 = new PhylogenyBranch( b1, a1 );
            if ( !a1b1.equals( a1b1 ) ) {
                return false;
            }
            if ( !a1b1.equals( b1a1 ) ) {
                return false;
            }
            if ( !b1a1.equals( a1b1 ) ) {
                return false;
            }
            final PhylogenyBranch a1_b1 = new PhylogenyBranch( a1, b1, true );
            final PhylogenyBranch b1_a1 = new PhylogenyBranch( b1, a1, true );
            final PhylogenyBranch a1_b1_ = new PhylogenyBranch( a1, b1, false );
            if ( a1_b1.equals( b1_a1 ) ) {
                return false;
            }
            if ( a1_b1.equals( a1_b1_ ) ) {
                return false;
            }
            final PhylogenyBranch b1_a1_ = new PhylogenyBranch( b1, a1, false );
            if ( !a1_b1.equals( b1_a1_ ) ) {
                return false;
            }
            if ( a1_b1_.equals( b1_a1_ ) ) {
                return false;
            }
            if ( !a1_b1_.equals( b1_a1 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPhyloXMLparsingOfDistributionElement() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            PhyloXmlParser xml_parser = null;
            try {
                xml_parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
            }
            catch ( final Exception e ) {
                // Do nothing -- means were not running from jar.
            }
            if ( xml_parser == null ) {
                xml_parser = PhyloXmlParser.createPhyloXmlParser();
                if ( USE_LOCAL_PHYLOXML_SCHEMA ) {
                    xml_parser.setValidateAgainstSchema( PHYLOXML_LOCAL_XSD );
                }
                else {
                    xml_parser.setValidateAgainstSchema( PHYLOXML_REMOTE_XSD );
                }
            }
            final Phylogeny[] phylogenies_0 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_distribution.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_0.length != 1 ) {
                return false;
            }
            final Phylogeny t1 = phylogenies_0[ 0 ];
            PhylogenyNode n = null;
            Distribution d = null;
            n = t1.getNode( "root node" );
            if ( !n.getNodeData().isHasDistribution() ) {
                return false;
            }
            if ( n.getNodeData().getDistributions().size() != 1 ) {
                return false;
            }
            d = n.getNodeData().getDistribution();
            if ( !d.getDesc().equals( "Hirschweg 38" ) ) {
                return false;
            }
            if ( d.getPoints().size() != 1 ) {
                return false;
            }
            if ( d.getPolygons() != null ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getAltitude().toString().equals( "472" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getAltiudeUnit().equals( "m" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getGeodeticDatum().equals( "WGS84" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLatitude().toString().equals( "47.48148427110029" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLongitude().toString().equals( "8.768951296806335" ) ) {
                return false;
            }
            n = t1.getNode( "node a" );
            if ( !n.getNodeData().isHasDistribution() ) {
                return false;
            }
            if ( n.getNodeData().getDistributions().size() != 2 ) {
                return false;
            }
            d = n.getNodeData().getDistribution( 1 );
            if ( !d.getDesc().equals( "San Diego" ) ) {
                return false;
            }
            if ( d.getPoints().size() != 1 ) {
                return false;
            }
            if ( d.getPolygons() != null ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getAltitude().toString().equals( "104" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getAltiudeUnit().equals( "m" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getGeodeticDatum().equals( "WGS84" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLatitude().toString().equals( "32.880933" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLongitude().toString().equals( "-117.217543" ) ) {
                return false;
            }
            n = t1.getNode( "node bb" );
            if ( !n.getNodeData().isHasDistribution() ) {
                return false;
            }
            if ( n.getNodeData().getDistributions().size() != 1 ) {
                return false;
            }
            d = n.getNodeData().getDistribution( 0 );
            if ( d.getPoints().size() != 3 ) {
                return false;
            }
            if ( d.getPolygons().size() != 2 ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLatitude().toString().equals( "1" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLongitude().toString().equals( "2" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 1 ).getLatitude().toString().equals( "3" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 1 ).getLongitude().toString().equals( "4" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 2 ).getLatitude().toString().equals( "5" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 2 ).getLongitude().toString().equals( "6" ) ) {
                return false;
            }
            Polygon p = d.getPolygons().get( 0 );
            if ( p.getPoints().size() != 3 ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getLatitude().toString().equals( "0.1" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getLongitude().toString().equals( "0.2" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getAltitude().toString().equals( "10" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 2 ).getLatitude().toString().equals( "0.5" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 2 ).getLongitude().toString().equals( "0.6" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 2 ).getAltitude().toString().equals( "30" ) ) {
                return false;
            }
            p = d.getPolygons().get( 1 );
            if ( p.getPoints().size() != 3 ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getLatitude().toString().equals( "1.49348902489947473" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getLongitude().toString().equals( "2.567489393947847492" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getAltitude().toString().equals( "10" ) ) {
                return false;
            }
            // Roundtrip:
            final StringBuffer t1_sb = new StringBuffer( t1.toPhyloXML( 0 ) );
            final Phylogeny[] rt = factory.create( t1_sb, xml_parser );
            if ( rt.length != 1 ) {
                return false;
            }
            final Phylogeny t1_rt = rt[ 0 ];
            n = t1_rt.getNode( "root node" );
            if ( !n.getNodeData().isHasDistribution() ) {
                return false;
            }
            if ( n.getNodeData().getDistributions().size() != 1 ) {
                return false;
            }
            d = n.getNodeData().getDistribution();
            if ( !d.getDesc().equals( "Hirschweg 38" ) ) {
                return false;
            }
            if ( d.getPoints().size() != 1 ) {
                return false;
            }
            if ( d.getPolygons() != null ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getAltitude().toString().equals( "472" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getAltiudeUnit().equals( "m" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getGeodeticDatum().equals( "WGS84" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLatitude().toString().equals( "47.48148427110029" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLongitude().toString().equals( "8.768951296806335" ) ) {
                return false;
            }
            n = t1_rt.getNode( "node a" );
            if ( !n.getNodeData().isHasDistribution() ) {
                return false;
            }
            if ( n.getNodeData().getDistributions().size() != 2 ) {
                return false;
            }
            d = n.getNodeData().getDistribution( 1 );
            if ( !d.getDesc().equals( "San Diego" ) ) {
                return false;
            }
            if ( d.getPoints().size() != 1 ) {
                return false;
            }
            if ( d.getPolygons() != null ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getAltitude().toString().equals( "104" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getAltiudeUnit().equals( "m" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getGeodeticDatum().equals( "WGS84" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLatitude().toString().equals( "32.880933" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLongitude().toString().equals( "-117.217543" ) ) {
                return false;
            }
            n = t1_rt.getNode( "node bb" );
            if ( !n.getNodeData().isHasDistribution() ) {
                return false;
            }
            if ( n.getNodeData().getDistributions().size() != 1 ) {
                return false;
            }
            d = n.getNodeData().getDistribution( 0 );
            if ( d.getPoints().size() != 3 ) {
                return false;
            }
            if ( d.getPolygons().size() != 2 ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLatitude().toString().equals( "1" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 0 ).getLongitude().toString().equals( "2" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 1 ).getLatitude().toString().equals( "3" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 1 ).getLongitude().toString().equals( "4" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 2 ).getLatitude().toString().equals( "5" ) ) {
                return false;
            }
            if ( !d.getPoints().get( 2 ).getLongitude().toString().equals( "6" ) ) {
                return false;
            }
            p = d.getPolygons().get( 0 );
            if ( p.getPoints().size() != 3 ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getLatitude().toString().equals( "0.1" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getLongitude().toString().equals( "0.2" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getAltitude().toString().equals( "10" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 2 ).getLatitude().toString().equals( "0.5" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 2 ).getLongitude().toString().equals( "0.6" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 2 ).getAltitude().toString().equals( "30" ) ) {
                return false;
            }
            p = d.getPolygons().get( 1 );
            if ( p.getPoints().size() != 3 ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getLatitude().toString().equals( "1.49348902489947473" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getLongitude().toString().equals( "2.567489393947847492" ) ) {
                return false;
            }
            if ( !p.getPoints().get( 0 ).getAltitude().toString().equals( "10" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPostOrderIterator() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((A,B)ab,(C,D)cd)r", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it0;
            for( it0 = t0.iteratorPostorder(); it0.hasNext(); ) {
                it0.next();
            }
            for( it0.reset(); it0.hasNext(); ) {
                it0.next();
            }
            final Phylogeny t1 = factory.create( "(((A,B)ab,(C,D)cd)abcd,((E,F)ef,(G,H)gh)efgh)r", new NHXParser() )[ 0 ];
            final PhylogenyNodeIterator it = t1.iteratorPostorder();
            if ( !it.next().getName().equals( "A" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "B" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "ab" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "C" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "D" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "cd" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "abcd" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "E" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "F" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "ef" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "G" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "H" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "gh" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "efgh" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "r" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPreOrderIterator() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((A,B)ab,(C,D)cd)r", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it0;
            for( it0 = t0.iteratorPreorder(); it0.hasNext(); ) {
                it0.next();
            }
            for( it0.reset(); it0.hasNext(); ) {
                it0.next();
            }
            PhylogenyNodeIterator it = t0.iteratorPreorder();
            if ( !it.next().getName().equals( "r" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "ab" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "A" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "B" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "cd" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "C" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "D" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
            final Phylogeny t1 = factory.create( "(((A,B)ab,(C,D)cd)abcd,((E,F)ef,(G,H)gh)efgh)r", new NHXParser() )[ 0 ];
            it = t1.iteratorPreorder();
            if ( !it.next().getName().equals( "r" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "abcd" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "ab" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "A" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "B" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "cd" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "C" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "D" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "efgh" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "ef" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "E" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "F" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "gh" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "G" ) ) {
                return false;
            }
            if ( !it.next().getName().equals( "H" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPropertiesMap() {
        try {
            final PropertiesMap pm = new PropertiesMap();
            final Property p0 = new Property( "dimensions:diameter", "1", "metric:mm", "xsd:decimal", AppliesTo.NODE );
            final Property p1 = new Property( "dimensions:length", "2", "metric:mm", "xsd:decimal", AppliesTo.NODE );
            final Property p2 = new Property( "something:else",
                                              "?",
                                              "improbable:research",
                                              "xsd:decimal",
                                              AppliesTo.NODE );
            pm.addProperty( p0 );
            pm.addProperty( p1 );
            pm.addProperty( p2 );
            if ( !pm.getProperty( "dimensions:diameter" ).getValue().equals( "1" ) ) {
                return false;
            }
            if ( !pm.getProperty( "dimensions:length" ).getValue().equals( "2" ) ) {
                return false;
            }
            if ( pm.getProperties().size() != 3 ) {
                return false;
            }
            if ( pm.getPropertiesWithGivenReferencePrefix( "dimensions" ).size() != 2 ) {
                return false;
            }
            if ( pm.getPropertiesWithGivenReferencePrefix( "something" ).size() != 1 ) {
                return false;
            }
            if ( pm.getProperties().size() != 3 ) {
                return false;
            }
            pm.removeProperty( "dimensions:diameter" );
            if ( pm.getProperties().size() != 2 ) {
                return false;
            }
            if ( pm.getPropertiesWithGivenReferencePrefix( "dimensions" ).size() != 1 ) {
                return false;
            }
            if ( pm.getPropertiesWithGivenReferencePrefix( "something" ).size() != 1 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testProteinId() {
        try {
            final ProteinId id1 = new ProteinId( "a" );
            final ProteinId id2 = new ProteinId( "a" );
            final ProteinId id3 = new ProteinId( "A" );
            final ProteinId id4 = new ProteinId( "b" );
            if ( !id1.equals( id1 ) ) {
                return false;
            }
            if ( id1.getId().equals( "x" ) ) {
                return false;
            }
            if ( id1.getId().equals( null ) ) {
                return false;
            }
            if ( !id1.equals( id2 ) ) {
                return false;
            }
            if ( id1.equals( id3 ) ) {
                return false;
            }
            if ( id1.hashCode() != id1.hashCode() ) {
                return false;
            }
            if ( id1.hashCode() != id2.hashCode() ) {
                return false;
            }
            if ( id1.hashCode() == id3.hashCode() ) {
                return false;
            }
            if ( id1.compareTo( id1 ) != 0 ) {
                return false;
            }
            if ( id1.compareTo( id2 ) != 0 ) {
                return false;
            }
            if ( id1.compareTo( id3 ) != 0 ) {
                return false;
            }
            if ( id1.compareTo( id4 ) >= 0 ) {
                return false;
            }
            if ( id4.compareTo( id1 ) <= 0 ) {
                return false;
            }
            if ( !id4.getId().equals( "b" ) ) {
                return false;
            }
            final ProteinId id5 = new ProteinId( " C " );
            if ( !id5.getId().equals( "C" ) ) {
                return false;
            }
            if ( id5.equals( id1 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testReIdMethods() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p = factory.create( "((1,2)A,(((X,Y,Z)a,b)3)B,(4,5,6)C)r", new NHXParser() )[ 0 ];
            final long count = PhylogenyNode.getNodeCount();
            p.levelOrderReID();
            if ( p.getNode( "r" ).getId() != count ) {
                return false;
            }
            if ( p.getNode( "A" ).getId() != ( count + 1 ) ) {
                return false;
            }
            if ( p.getNode( "B" ).getId() != ( count + 1 ) ) {
                return false;
            }
            if ( p.getNode( "C" ).getId() != ( count + 1 ) ) {
                return false;
            }
            if ( p.getNode( "1" ).getId() != ( count + 2 ) ) {
                return false;
            }
            if ( p.getNode( "2" ).getId() != ( count + 2 ) ) {
                return false;
            }
            if ( p.getNode( "3" ).getId() != ( count + 2 ) ) {
                return false;
            }
            if ( p.getNode( "4" ).getId() != ( count + 2 ) ) {
                return false;
            }
            if ( p.getNode( "5" ).getId() != ( count + 2 ) ) {
                return false;
            }
            if ( p.getNode( "6" ).getId() != ( count + 2 ) ) {
                return false;
            }
            if ( p.getNode( "a" ).getId() != ( count + 3 ) ) {
                return false;
            }
            if ( p.getNode( "b" ).getId() != ( count + 3 ) ) {
                return false;
            }
            if ( p.getNode( "X" ).getId() != ( count + 4 ) ) {
                return false;
            }
            if ( p.getNode( "Y" ).getId() != ( count + 4 ) ) {
                return false;
            }
            if ( p.getNode( "Z" ).getId() != ( count + 4 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testRerooting() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t1 = factory.create( "((A:1,B:2)AB:1[&&NHX:B=55],(C:3,D:5)CD:3[&&NHX:B=10])ABCD:0.5",
                                                 new NHXParser() )[ 0 ];
            if ( !t1.isRooted() ) {
                return false;
            }
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "CD" ) );
            t1.reRoot( t1.getNode( "A" ) );
            t1.reRoot( t1.getNode( "B" ) );
            t1.reRoot( t1.getNode( "AB" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "C" ) );
            t1.reRoot( t1.getNode( "CD" ) );
            t1.reRoot( t1.getNode( "A" ) );
            t1.reRoot( t1.getNode( "B" ) );
            t1.reRoot( t1.getNode( "AB" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "C" ) );
            t1.reRoot( t1.getNode( "A" ) );
            t1.reRoot( t1.getNode( "B" ) );
            t1.reRoot( t1.getNode( "AB" ) );
            t1.reRoot( t1.getNode( "C" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "CD" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "A" ) );
            t1.reRoot( t1.getNode( "B" ) );
            t1.reRoot( t1.getNode( "AB" ) );
            t1.reRoot( t1.getNode( "C" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "CD" ) );
            t1.reRoot( t1.getNode( "D" ) );
            if ( !isEqual( t1.getNode( "A" ).getDistanceToParent(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "B" ).getDistanceToParent(), 2 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "C" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "D" ).getDistanceToParent(), 2.5 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "CD" ).getDistanceToParent(), 2.5 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "AB" ).getDistanceToParent(), 4 ) ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "(((A:1,B:2)AB:10[&&NHX:B=55],C)ABC:3[&&NHX:B=33],D:5)ABCD:0.5",
                                                 new NHXParser() )[ 0 ];
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "D" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "ABC" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "AB" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "D" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "AB" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "D" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "D" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "ABC" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            final Phylogeny t3 = factory.create( "(A[&&NHX:B=10],B[&&NHX:B=20],C[&&NHX:B=30],D[&&NHX:B=40])",
                                                 new NHXParser() )[ 0 ];
            t3.reRoot( t3.getNode( "B" ) );
            if ( t3.getNode( "B" ).getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getNumberOfDescendants() != 3 ) {
                return false;
            }
            t3.reRoot( t3.getNode( "B" ) );
            if ( t3.getNode( "B" ).getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getNumberOfDescendants() != 3 ) {
                return false;
            }
            t3.reRoot( t3.getRoot() );
            if ( t3.getNode( "B" ).getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getNumberOfDescendants() != 3 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSDIse() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny species1 = factory.create( "[&&NHX:S=yeast]", new NHXParser() )[ 0 ];
            final Phylogeny gene1 = factory.create( "(A1[&&NHX:S=yeast],A2[&&NHX:S=yeast])", new NHXParser() )[ 0 ];
            gene1.setRooted( true );
            species1.setRooted( true );
            final SDI sdi = new SDI( gene1, species1 );
            if ( !gene1.getRoot().isDuplication() ) {
                return false;
            }
            final Phylogeny species2 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene2 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B])ab,[&&NHX:S=C])abc,[&&NHX:S=D])abcd,([&&NHX:S=E],[&&NHX:S=F])ef)r",
                             new NHXParser() )[ 0 ];
            species2.setRooted( true );
            gene2.setRooted( true );
            final SDI sdi2 = new SDI( gene2, species2 );
            if ( sdi2.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( !gene2.getNode( "ab" ).isSpeciation() ) {
                return false;
            }
            if ( !gene2.getNode( "ab" ).isHasAssignedEvent() ) {
                return false;
            }
            if ( !gene2.getNode( "abc" ).isSpeciation() ) {
                return false;
            }
            if ( !gene2.getNode( "abc" ).isHasAssignedEvent() ) {
                return false;
            }
            if ( !gene2.getNode( "r" ).isSpeciation() ) {
                return false;
            }
            if ( !gene2.getNode( "r" ).isHasAssignedEvent() ) {
                return false;
            }
            final Phylogeny species3 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene3 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=A])aa,[&&NHX:S=C])abc,[&&NHX:S=D])abcd,([&&NHX:S=E],[&&NHX:S=F])ef)r",
                             new NHXParser() )[ 0 ];
            species3.setRooted( true );
            gene3.setRooted( true );
            final SDI sdi3 = new SDI( gene3, species3 );
            if ( sdi3.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( !gene3.getNode( "aa" ).isDuplication() ) {
                return false;
            }
            if ( !gene3.getNode( "aa" ).isHasAssignedEvent() ) {
                return false;
            }
            final Phylogeny species4 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene4 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=C])ac,[&&NHX:S=B])abc,[&&NHX:S=D])abcd,([&&NHX:S=E],[&&NHX:S=F])ef)r",
                             new NHXParser() )[ 0 ];
            species4.setRooted( true );
            gene4.setRooted( true );
            final SDI sdi4 = new SDI( gene4, species4 );
            if ( sdi4.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( !gene4.getNode( "ac" ).isSpeciation() ) {
                return false;
            }
            if ( !gene4.getNode( "abc" ).isDuplication() ) {
                return false;
            }
            if ( gene4.getNode( "abcd" ).isDuplication() ) {
                return false;
            }
            if ( species4.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            if ( gene4.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            final Phylogeny species5 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene5 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=D])ad,[&&NHX:S=C])adc,[&&NHX:S=B])abcd,([&&NHX:S=E],[&&NHX:S=F])ef)r",
                             new NHXParser() )[ 0 ];
            species5.setRooted( true );
            gene5.setRooted( true );
            final SDI sdi5 = new SDI( gene5, species5 );
            if ( sdi5.getDuplicationsSum() != 2 ) {
                return false;
            }
            if ( !gene5.getNode( "ad" ).isSpeciation() ) {
                return false;
            }
            if ( !gene5.getNode( "adc" ).isDuplication() ) {
                return false;
            }
            if ( !gene5.getNode( "abcd" ).isDuplication() ) {
                return false;
            }
            if ( species5.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            if ( gene5.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            // Trees from Louxin Zhang 1997 "On a Mirkin-Muchnik-Smith
            // Conjecture for Comparing Molecular Phylogenies"
            // J. of Comput Bio. Vol. 4, No 2, pp.177-187
            final Phylogeny species6 = factory
                    .create( "(((1:[&&NHX:S=1],5:[&&NHX:S=5])1-5,((4:[&&NHX:S=4],6:[&&NHX:S=6])4-6,2:[&&NHX:S=2])4-6-2)1-5-4-6-2,"
                                     + "((9:[&&NHX:S=9],3:[&&NHX:S=3])9-3,(8:[&&NHX:S=8],7:[&&NHX:S=7])8-7)9-3-8-7)",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene6 = factory
                    .create( "(((1:0.1[&&NHX:S=1],2:0.1[&&NHX:S=2])1-2:0.1,3:0.1[&&NHX:S=3])1-2-3:0.1,"
                                     + "((4:0.1[&&NHX:S=4],(5:0.1[&&NHX:S=5],6:0.1[&&NHX:S=6])5-6:0.1)4-5-6:0.1,"
                                     + "(7:0.1[&&NHX:S=7],(8:0.1[&&NHX:S=8],9:0.1[&&NHX:S=9])8-9:0.1)7-8-9:0.1)4-5-6-7-8-9:0.1)r;",
                             new NHXParser() )[ 0 ];
            species6.setRooted( true );
            gene6.setRooted( true );
            final SDI sdi6 = new SDI( gene6, species6 );
            if ( sdi6.getDuplicationsSum() != 3 ) {
                return false;
            }
            if ( !gene6.getNode( "r" ).isDuplication() ) {
                return false;
            }
            if ( !gene6.getNode( "4-5-6" ).isDuplication() ) {
                return false;
            }
            if ( !gene6.getNode( "7-8-9" ).isDuplication() ) {
                return false;
            }
            if ( !gene6.getNode( "1-2" ).isSpeciation() ) {
                return false;
            }
            if ( !gene6.getNode( "1-2-3" ).isSpeciation() ) {
                return false;
            }
            if ( !gene6.getNode( "5-6" ).isSpeciation() ) {
                return false;
            }
            if ( !gene6.getNode( "8-9" ).isSpeciation() ) {
                return false;
            }
            if ( !gene6.getNode( "4-5-6-7-8-9" ).isSpeciation() ) {
                return false;
            }
            sdi6.computeMappingCostL();
            if ( sdi6.computeMappingCostL() != 17 ) {
                return false;
            }
            if ( species6.getNumberOfExternalNodes() != 9 ) {
                return false;
            }
            if ( gene6.getNumberOfExternalNodes() != 9 ) {
                return false;
            }
            final Phylogeny species7 = Test.createPhylogeny( "(((((((" + "([&&NHX:S=a1],[&&NHX:S=a2]),"
                    + "([&&NHX:S=b1],[&&NHX:S=b2])" + "),[&&NHX:S=x]),(" + "([&&NHX:S=m1],[&&NHX:S=m2]),"
                    + "([&&NHX:S=n1],[&&NHX:S=n2])" + ")),(" + "([&&NHX:S=i1],[&&NHX:S=i2]),"
                    + "([&&NHX:S=j1],[&&NHX:S=j2])" + ")),(" + "([&&NHX:S=e1],[&&NHX:S=e2]),"
                    + "([&&NHX:S=f1],[&&NHX:S=f2])" + ")),[&&NHX:S=y]),[&&NHX:S=z])" );
            species7.setRooted( true );
            final Phylogeny gene7_1 = Test
                    .createPhylogeny( "((((((((a1[&&NHX:S=a1],a2[&&NHX:S=a2]),b1[&&NHX:S=b1]),x[&&NHX:S=x]),m1[&&NHX:S=m1]),i1[&&NHX:S=i1]),e1[&&NHX:S=e1]),y[&&NHX:S=y]),z[&&NHX:S=z])" );
            gene7_1.setRooted( true );
            final SDI sdi7 = new SDI( gene7_1, species7 );
            if ( sdi7.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "a2" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "x" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "m1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "i1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "e1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "y" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "z" ).isSpeciation() ) {
                return false;
            }
            final Phylogeny gene7_2 = Test
                    .createPhylogeny( "(((((((((a1[&&NHX:S=a1],a2[&&NHX:S=a2]),b1[&&NHX:S=b1]),x[&&NHX:S=x]),m1[&&NHX:S=m1]),i1[&&NHX:S=i1]),j2[&&NHX:S=j2]),e1[&&NHX:S=e1]),y[&&NHX:S=y]),z[&&NHX:S=z])" );
            gene7_2.setRooted( true );
            final SDI sdi7_2 = new SDI( gene7_2, species7 );
            if ( sdi7_2.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "a2" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "x" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "m1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "i1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "j2" ).isDuplication() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "e1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "y" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "z" ).isSpeciation() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            return false;
        }
        return true;
    }

    private static boolean testSDIunrooted() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p0 = factory.create( "((((A,B)ab,(C1,C2)cc)abc,D)abcd,(E,F)ef)abcdef", new NHXParser() )[ 0 ];
            final List<PhylogenyBranch> l = SDIR.getBranchesInPreorder( p0 );
            final Iterator<PhylogenyBranch> iter = l.iterator();
            PhylogenyBranch br = iter.next();
            if ( !br.getFirstNode().getName().equals( "abcd" ) && !br.getFirstNode().getName().equals( "ef" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "abcd" ) && !br.getSecondNode().getName().equals( "ef" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "abcd" ) && !br.getFirstNode().getName().equals( "abc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "abcd" ) && !br.getSecondNode().getName().equals( "abc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "abc" ) && !br.getFirstNode().getName().equals( "ab" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "abc" ) && !br.getSecondNode().getName().equals( "ab" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "ab" ) && !br.getFirstNode().getName().equals( "A" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ab" ) && !br.getSecondNode().getName().equals( "A" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "ab" ) && !br.getFirstNode().getName().equals( "B" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ab" ) && !br.getSecondNode().getName().equals( "B" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "ab" ) && !br.getFirstNode().getName().equals( "abc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ab" ) && !br.getSecondNode().getName().equals( "abc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "abc" ) && !br.getFirstNode().getName().equals( "cc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "abc" ) && !br.getSecondNode().getName().equals( "cc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "C1" ) && !br.getFirstNode().getName().equals( "cc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "C1" ) && !br.getSecondNode().getName().equals( "cc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "C2" ) && !br.getFirstNode().getName().equals( "cc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "C2" ) && !br.getSecondNode().getName().equals( "cc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "abc" ) && !br.getFirstNode().getName().equals( "cc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "abc" ) && !br.getSecondNode().getName().equals( "cc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "abc" ) && !br.getFirstNode().getName().equals( "abcd" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "abc" ) && !br.getSecondNode().getName().equals( "abcd" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "abcd" ) && !br.getFirstNode().getName().equals( "D" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "abcd" ) && !br.getSecondNode().getName().equals( "D" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "ef" ) && !br.getFirstNode().getName().equals( "abcd" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ef" ) && !br.getSecondNode().getName().equals( "abcd" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "ef" ) && !br.getFirstNode().getName().equals( "E" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ef" ) && !br.getSecondNode().getName().equals( "E" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getName().equals( "ef" ) && !br.getFirstNode().getName().equals( "F" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ef" ) && !br.getSecondNode().getName().equals( "F" ) ) {
                return false;
            }
            if ( iter.hasNext() ) {
                return false;
            }
            final Phylogeny p1 = factory.create( "(C,(A,B)ab)abc", new NHXParser() )[ 0 ];
            final List<PhylogenyBranch> l1 = SDIR.getBranchesInPreorder( p1 );
            final Iterator<PhylogenyBranch> iter1 = l1.iterator();
            br = iter1.next();
            if ( !br.getFirstNode().getName().equals( "ab" ) && !br.getFirstNode().getName().equals( "C" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ab" ) && !br.getSecondNode().getName().equals( "C" ) ) {
                return false;
            }
            br = iter1.next();
            if ( !br.getFirstNode().getName().equals( "ab" ) && !br.getFirstNode().getName().equals( "A" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ab" ) && !br.getSecondNode().getName().equals( "A" ) ) {
                return false;
            }
            br = iter1.next();
            if ( !br.getFirstNode().getName().equals( "ab" ) && !br.getFirstNode().getName().equals( "B" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ab" ) && !br.getSecondNode().getName().equals( "B" ) ) {
                return false;
            }
            if ( iter1.hasNext() ) {
                return false;
            }
            final Phylogeny p2 = factory.create( "((A,B)ab,C)abc", new NHXParser() )[ 0 ];
            final List<PhylogenyBranch> l2 = SDIR.getBranchesInPreorder( p2 );
            final Iterator<PhylogenyBranch> iter2 = l2.iterator();
            br = iter2.next();
            if ( !br.getFirstNode().getName().equals( "ab" ) && !br.getFirstNode().getName().equals( "C" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ab" ) && !br.getSecondNode().getName().equals( "C" ) ) {
                return false;
            }
            br = iter2.next();
            if ( !br.getFirstNode().getName().equals( "ab" ) && !br.getFirstNode().getName().equals( "A" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ab" ) && !br.getSecondNode().getName().equals( "A" ) ) {
                return false;
            }
            br = iter2.next();
            if ( !br.getFirstNode().getName().equals( "ab" ) && !br.getFirstNode().getName().equals( "B" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getName().equals( "ab" ) && !br.getSecondNode().getName().equals( "B" ) ) {
                return false;
            }
            if ( iter2.hasNext() ) {
                return false;
            }
            final Phylogeny species0 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene1 = factory
                    .create( "(((((A:0.6[&&NHX:S=A],B:0.1[&&NHX:S=B])ab:0.1,C:0.1[&&NHX:S=C])abc:0.3,D:1.0[&&NHX:S=D])abcd:0.2,E:0.1[&&NHX:S=E])abcde:0.2,F:0.2[&&NHX:S=F])",
                             new NHXParser() )[ 0 ];
            species0.setRooted( true );
            gene1.setRooted( true );
            final SDIR sdi_unrooted = new SDIR();
            sdi_unrooted.infer( gene1, species0, false, true, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 0 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.4 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 1.0 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != Integer.MAX_VALUE ) {
                return false;
            }
            final Phylogeny gene2 = factory
                    .create( "(((((A:2.6[&&NHX:S=A],B:0.1[&&NHX:S=B])ab:0.1,C:0.1[&&NHX:S=C])abc:0.3,D:1.0[&&NHX:S=D])abcd:0.2,E:0.1[&&NHX:S=E])abcde:0.2,F:0.2[&&NHX:S=F])",
                             new NHXParser() )[ 0 ];
            gene2.setRooted( true );
            sdi_unrooted.infer( gene2, species0, false, false, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 3 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 2.0 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != Integer.MAX_VALUE ) {
                return false;
            }
            final Phylogeny species6 = factory
                    .create( "(((1:[&&NHX:S=1],5:[&&NHX:S=5])1-5,((4:[&&NHX:S=4],6:[&&NHX:S=6])4-6,2:[&&NHX:S=2])4-6-2)1-5-4-6-2,"
                                     + "((9:[&&NHX:S=9],3:[&&NHX:S=3])9-3,(8:[&&NHX:S=8],7:[&&NHX:S=7])8-7)9-3-8-7)",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene6 = factory
                    .create( "((5:0.1[&&NHX:S=5],6:0.1[&&NHX:S=6])5-6:0.05[&&NHX:S=6],(4:0.1[&&NHX:S=4],"
                                     + "(((1:0.1[&&NHX:S=1],2:0.1[&&NHX:S=2])1-2:0.1[&&NHX:S=2],3:0.25[&&NHX:S=3])1-2-3:0.2[&&NHX:S=2],"
                                     + "(7:0.1[&&NHX:S=7],(8:0.1[&&NHX:S=8],"
                                     + "9:0.1[&&NHX:S=9])8-9:0.1[&&NHX:S=9])7-8-9:0.1[&&NHX:S=8])"
                                     + "4-5-6-7-8-9:0.1[&&NHX:S=5])4-5-6:0.05[&&NHX:S=5])",
                             new NHXParser() )[ 0 ];
            species6.setRooted( true );
            gene6.setRooted( true );
            Phylogeny[] p6 = sdi_unrooted.infer( gene6, species6, false, true, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 0.375 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 3 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != Integer.MAX_VALUE ) {
                return false;
            }
            if ( !p6[ 0 ].getRoot().isDuplication() ) {
                return false;
            }
            if ( !p6[ 0 ].getNode( "4-5-6" ).isDuplication() ) {
                return false;
            }
            if ( !p6[ 0 ].getNode( "7-8-9" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "1-2" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "1-2-3" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "5-6" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "8-9" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "4-5-6-7-8-9" ).isDuplication() ) {
                return false;
            }
            p6 = null;
            final Phylogeny species7 = factory
                    .create( "(((1:[&&NHX:S=1],5:[&&NHX:S=5])1-5,((4:[&&NHX:S=4],6:[&&NHX:S=6])4-6,2:[&&NHX:S=2])4-6-2)1-5-4-6-2,"
                                     + "((9:[&&NHX:S=9],3:[&&NHX:S=3])9-3,(8:[&&NHX:S=8],7:[&&NHX:S=7])8-7)9-3-8-7)",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene7 = factory
                    .create( "((5:0.1[&&NHX:S=5],6:0.1[&&NHX:S=6])5-6:0.05[&&NHX:S=6],(4:0.1[&&NHX:S=4],"
                                     + "(((1:0.1[&&NHX:S=1],2:0.1[&&NHX:S=2])1-2:0.1[&&NHX:S=2],3:0.25[&&NHX:S=3])1-2-3:0.2[&&NHX:S=2],"
                                     + "(7:0.1[&&NHX:S=7],(8:0.1[&&NHX:S=8],"
                                     + "9:0.1[&&NHX:S=9])8-9:0.1[&&NHX:S=9])7-8-9:0.1[&&NHX:S=8])"
                                     + "4-5-6-7-8-9:0.1[&&NHX:S=5])4-5-6:0.05[&&NHX:S=5])",
                             new NHXParser() )[ 0 ];
            species7.setRooted( true );
            gene7.setRooted( true );
            Phylogeny[] p7 = sdi_unrooted.infer( gene7, species7, true, true, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 0.375 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 3 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != 17 ) {
                return false;
            }
            if ( !p7[ 0 ].getRoot().isDuplication() ) {
                return false;
            }
            if ( !p7[ 0 ].getNode( "4-5-6" ).isDuplication() ) {
                return false;
            }
            if ( !p7[ 0 ].getNode( "7-8-9" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "1-2" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "1-2-3" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "5-6" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "8-9" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "4-5-6-7-8-9" ).isDuplication() ) {
                return false;
            }
            p7 = null;
            final Phylogeny species8 = factory
                    .create( "(((1:[&&NHX:S=1],5:[&&NHX:S=5])1-5,((4:[&&NHX:S=4],6:[&&NHX:S=6])4-6,2:[&&NHX:S=2])4-6-2)1-5-4-6-2,"
                                     + "((9:[&&NHX:S=9],3:[&&NHX:S=3])9-3,(8:[&&NHX:S=8],7:[&&NHX:S=7])8-7)9-3-8-7)",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene8 = factory
                    .create( "((5:0.1[&&NHX:S=5],6:0.1[&&NHX:S=6])5-6:0.05[&&NHX:S=6],(4:0.1[&&NHX:S=4],"
                                     + "(((1:0.1[&&NHX:S=1],2:0.1[&&NHX:S=2])1-2:0.1[&&NHX:S=2],3:0.25[&&NHX:S=3])1-2-3:0.2[&&NHX:S=2],"
                                     + "(7:0.1[&&NHX:S=7],(8:0.1[&&NHX:S=8],"
                                     + "9:0.1[&&NHX:S=9])8-9:0.1[&&NHX:S=9])7-8-9:0.1[&&NHX:S=8])"
                                     + "4-5-6-7-8-9:0.1[&&NHX:S=5])4-5-6:0.05[&&NHX:S=5])",
                             new NHXParser() )[ 0 ];
            species8.setRooted( true );
            gene8.setRooted( true );
            Phylogeny[] p8 = sdi_unrooted.infer( gene8, species8, false, false, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 0.375 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 3 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != Integer.MAX_VALUE ) {
                return false;
            }
            if ( !p8[ 0 ].getRoot().isDuplication() ) {
                return false;
            }
            if ( !p8[ 0 ].getNode( "4-5-6" ).isDuplication() ) {
                return false;
            }
            if ( !p8[ 0 ].getNode( "7-8-9" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "1-2" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "1-2-3" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "5-6" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "8-9" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "4-5-6-7-8-9" ).isDuplication() ) {
                return false;
            }
            p8 = null;
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSequenceIdParsing() {
        try {
            Accession id = SequenceAccessionTools.parseAccessorFromString( "gb_ADF31344_segmented_worms_" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "ADF31344" ) || !id.getSource().equals( "ncbi" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            //
            id = SequenceAccessionTools.parseAccessorFromString( "segmented worms|gb_ADF31344" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "ADF31344" ) || !id.getSource().equals( "ncbi" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            //
            id = SequenceAccessionTools.parseAccessorFromString( "segmented worms gb_ADF31344 and more" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "ADF31344" ) || !id.getSource().equals( "ncbi" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            // 
            id = SequenceAccessionTools.parseAccessorFromString( "gb_AAA96518_1" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "AAA96518" ) || !id.getSource().equals( "ncbi" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            // 
            id = SequenceAccessionTools.parseAccessorFromString( "gb_EHB07727_1_rodents_" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "EHB07727" ) || !id.getSource().equals( "ncbi" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            // 
            id = SequenceAccessionTools.parseAccessorFromString( "dbj_BAF37827_1_turtles_" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "BAF37827" ) || !id.getSource().equals( "ncbi" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            // 
            id = SequenceAccessionTools.parseAccessorFromString( "emb_CAA73223_1_primates_" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "CAA73223" ) || !id.getSource().equals( "ncbi" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            // 
            id = SequenceAccessionTools.parseAccessorFromString( "mites|ref_XP_002434188_1" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "XP_002434188" ) || !id.getSource().equals( "refseq" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            // 
            id = SequenceAccessionTools.parseAccessorFromString( "mites_ref_XP_002434188_1_bla_XP_12345" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "XP_002434188" ) || !id.getSource().equals( "refseq" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            // 
            id = SequenceAccessionTools.parseAccessorFromString( "P4A123" );
            if ( ( id == null ) || ForesterUtil.isEmpty( id.getValue() ) || ForesterUtil.isEmpty( id.getSource() )
                    || !id.getValue().equals( "P4A123" ) || !id.getSource().equals( "uniprot" ) ) {
                if ( id != null ) {
                    System.out.println( "value   =" + id.getValue() );
                    System.out.println( "provider=" + id.getSource() );
                }
                return false;
            }
            id = SequenceAccessionTools.parseAccessorFromString( "XP_12345" );
            if ( id != null ) {
                System.out.println( "value   =" + id.getValue() );
                System.out.println( "provider=" + id.getSource() );
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSequenceWriter() {
        try {
            final String n = ForesterUtil.LINE_SEPARATOR;
            if ( !SequenceWriter.toFasta( "name", "awes", 5 ).toString().equals( ">name" + n + "awes" ) ) {
                return false;
            }
            if ( !SequenceWriter.toFasta( "name", "awes", 4 ).toString().equals( ">name" + n + "awes" ) ) {
                return false;
            }
            if ( !SequenceWriter.toFasta( "name", "awes", 3 ).toString().equals( ">name" + n + "awe" + n + "s" ) ) {
                return false;
            }
            if ( !SequenceWriter.toFasta( "name", "awes", 2 ).toString().equals( ">name" + n + "aw" + n + "es" ) ) {
                return false;
            }
            if ( !SequenceWriter.toFasta( "name", "awes", 1 ).toString()
                    .equals( ">name" + n + "a" + n + "w" + n + "e" + n + "s" ) ) {
                return false;
            }
            if ( !SequenceWriter.toFasta( "name", "abcdefghij", 3 ).toString()
                    .equals( ">name" + n + "abc" + n + "def" + n + "ghi" + n + "j" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testSpecies() {
        try {
            final Species s1 = new BasicSpecies( "a" );
            final Species s2 = new BasicSpecies( "a" );
            final Species s3 = new BasicSpecies( "A" );
            final Species s4 = new BasicSpecies( "b" );
            if ( !s1.equals( s1 ) ) {
                return false;
            }
            if ( s1.getSpeciesId().equals( "x" ) ) {
                return false;
            }
            if ( s1.getSpeciesId().equals( null ) ) {
                return false;
            }
            if ( !s1.equals( s2 ) ) {
                return false;
            }
            if ( s1.equals( s3 ) ) {
                return false;
            }
            if ( s1.hashCode() != s1.hashCode() ) {
                return false;
            }
            if ( s1.hashCode() != s2.hashCode() ) {
                return false;
            }
            if ( s1.hashCode() == s3.hashCode() ) {
                return false;
            }
            if ( s1.compareTo( s1 ) != 0 ) {
                return false;
            }
            if ( s1.compareTo( s2 ) != 0 ) {
                return false;
            }
            if ( s1.compareTo( s3 ) != 0 ) {
                return false;
            }
            if ( s1.compareTo( s4 ) >= 0 ) {
                return false;
            }
            if ( s4.compareTo( s1 ) <= 0 ) {
                return false;
            }
            if ( !s4.getSpeciesId().equals( "b" ) ) {
                return false;
            }
            final Species s5 = new BasicSpecies( " C " );
            if ( !s5.getSpeciesId().equals( "C" ) ) {
                return false;
            }
            if ( s5.equals( s1 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSplit() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p0 = factory.create( "(((A,B,C),D),(E,(F,G)))R", new NHXParser() )[ 0 ];
            //Archaeopteryx.createApplication( p0 );
            final Set<PhylogenyNode> ex = new HashSet<PhylogenyNode>();
            ex.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            final TreeSplitMatrix s0 = new TreeSplitMatrix( p0, false, ex );
            // System.out.println( s0.toString() );
            //
            Set<PhylogenyNode> query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            /////////
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "A" ) );
            //            query_nodes.add( new PhylogenyNode( "B" ) );
            //            query_nodes.add( new PhylogenyNode( "C" ) );
            //            query_nodes.add( new PhylogenyNode( "D" ) );
            //            query_nodes.add( new PhylogenyNode( "E" ) );
            //            query_nodes.add( new PhylogenyNode( "F" ) );
            //            query_nodes.add( new PhylogenyNode( "G" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "A" ) );
            //            query_nodes.add( new PhylogenyNode( "B" ) );
            //            query_nodes.add( new PhylogenyNode( "C" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "D" ) );
            //            query_nodes.add( new PhylogenyNode( "E" ) );
            //            query_nodes.add( new PhylogenyNode( "F" ) );
            //            query_nodes.add( new PhylogenyNode( "G" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "A" ) );
            //            query_nodes.add( new PhylogenyNode( "B" ) );
            //            query_nodes.add( new PhylogenyNode( "C" ) );
            //            query_nodes.add( new PhylogenyNode( "D" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "E" ) );
            //            query_nodes.add( new PhylogenyNode( "F" ) );
            //            query_nodes.add( new PhylogenyNode( "G" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "F" ) );
            //            query_nodes.add( new PhylogenyNode( "G" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            ///////////////////////////
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "X" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "Y" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testSplitStrict() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p0 = factory.create( "(((A,B,C),D),(E,(F,G)))R", new NHXParser() )[ 0 ];
            final Set<PhylogenyNode> ex = new HashSet<PhylogenyNode>();
            ex.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            ex.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            final TreeSplitMatrix s0 = new TreeSplitMatrix( p0, true, ex );
            Set<PhylogenyNode> query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "F" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "B" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "E" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "D" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "A" ) );
            query_nodes.add( PhylogenyNode.createInstanceFromNhxString( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testSubtreeDeletion() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t1 = factory.create( "((A,B,C)abc,(D,E,F)def)r", new NHXParser() )[ 0 ];
            t1.deleteSubtree( t1.getNode( "A" ), false );
            if ( t1.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            t1.toNewHampshireX();
            t1.deleteSubtree( t1.getNode( "E" ), false );
            if ( t1.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            t1.toNewHampshireX();
            t1.deleteSubtree( t1.getNode( "F" ), false );
            if ( t1.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            t1.toNewHampshireX();
            t1.deleteSubtree( t1.getNode( "D" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "def" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "B" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "C" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "abc" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "r" ), false );
            if ( t1.getNumberOfExternalNodes() != 0 ) {
                return false;
            }
            if ( !t1.isEmpty() ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "(((1,2,3)A,B,C)abc,(D,E,F)def)r", new NHXParser() )[ 0 ];
            t2.deleteSubtree( t2.getNode( "A" ), false );
            t2.toNewHampshireX();
            if ( t2.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "abc" ), false );
            t2.toNewHampshireX();
            if ( t2.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "def" ), false );
            t2.toNewHampshireX();
            if ( t2.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSupportCount() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0_1 = factory.create( "(((A,B),C),(D,E))", new NHXParser() )[ 0 ];
            final Phylogeny[] phylogenies_1 = factory.create( "(((A,B),C),(D,E)) " + "(((C,B),A),(D,E))"
                                                                      + "(((A,B),C),(D,E)) " + "(((A,B),C),(D,E))"
                                                                      + "(((A,B),C),(D,E))" + "(((C,B),A),(D,E))"
                                                                      + "(((E,B),D),(C,A))" + "(((C,B),A),(D,E))"
                                                                      + "(((A,B),C),(D,E))" + "(((A,B),C),(D,E))",
                                                              new NHXParser() );
            SupportCount.count( t0_1, phylogenies_1, true, false );
            final Phylogeny t0_2 = factory.create( "(((((A,B),C),D),E),(F,G))", new NHXParser() )[ 0 ];
            final Phylogeny[] phylogenies_2 = factory.create( "(((((A,B),C),D),E),(F,G))"
                                                                      + "(((((A,B),C),D),E),((F,G),X))"
                                                                      + "(((((A,Y),B),C),D),((F,G),E))"
                                                                      + "(((((A,B),C),D),E),(F,G))"
                                                                      + "(((((A,B),C),D),E),(F,G))"
                                                                      + "(((((A,B),C),D),E),(F,G))"
                                                                      + "(((((A,B),C),D),E),(F,G),Z)"
                                                                      + "(((((A,B),C),D),E),(F,G))"
                                                                      + "((((((A,B),C),D),E),F),G)"
                                                                      + "(((((X,Y),F,G),E),((A,B),C)),D)",
                                                              new NHXParser() );
            SupportCount.count( t0_2, phylogenies_2, true, false );
            final PhylogenyNodeIterator it = t0_2.iteratorPostorder();
            while ( it.hasNext() ) {
                final PhylogenyNode n = it.next();
                if ( !n.isExternal() && ( PhylogenyMethods.getConfidenceValue( n ) != 10 ) ) {
                    return false;
                }
            }
            final Phylogeny t0_3 = factory.create( "(((A,B)ab,C)abc,((D,E)de,F)def)", new NHXParser() )[ 0 ];
            final Phylogeny[] phylogenies_3 = factory.create( "(((A,B),C),((D,E),F))" + "(((A,C),B),((D,F),E))"
                    + "(((C,A),B),((F,D),E))" + "(((A,B),F),((D,E),C))" + "(((((A,B),C),D),E),F)", new NHXParser() );
            SupportCount.count( t0_3, phylogenies_3, true, false );
            t0_3.reRoot( t0_3.getNode( "def" ).getId() );
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "ab" ) ) != 3 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "abc" ) ) != 4 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "def" ) ) != 4 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "de" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "A" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "B" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "C" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "D" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "E" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "F" ) ) != 5 ) {
                return false;
            }
            final Phylogeny t0_4 = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            final Phylogeny[] phylogenies_4 = factory.create( "((((((A,X),C),B),D),E),F) "
                    + "(((A,B,Z),C,Q),(((D,Y),E),F))", new NHXParser() );
            SupportCount.count( t0_4, phylogenies_4, true, false );
            t0_4.reRoot( t0_4.getNode( "F" ).getId() );
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "1" ) ) != 1 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "2" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "3" ) ) != 1 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "4" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "A" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "B" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "C" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "D" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "E" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "F" ) ) != 2 ) {
                return false;
            }
            Phylogeny a = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            final Phylogeny b1 = factory.create( "(((((B,A)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            double d = SupportCount.compare( b1, a, true, true, true );
            if ( !Test.isEqual( d, 5.0 / 5.0 ) ) {
                return false;
            }
            a = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            final Phylogeny b2 = factory.create( "(((((C,B)1,A)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            d = SupportCount.compare( b2, a, true, true, true );
            if ( !Test.isEqual( d, 4.0 / 5.0 ) ) {
                return false;
            }
            a = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            final Phylogeny b3 = factory.create( "(((((F,C)1,A)2,B)3,D)4,E)", new NHXParser() )[ 0 ];
            d = SupportCount.compare( b3, a, true, true, true );
            if ( !Test.isEqual( d, 2.0 / 5.0 ) ) {
                return false;
            }
            a = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)r", new NHXParser() )[ 0 ];
            final Phylogeny b4 = factory.create( "(((((F,C)1,A)2,B)3,D)4,E)r", new NHXParser() )[ 0 ];
            d = SupportCount.compare( b4, a, true, true, false );
            if ( !Test.isEqual( d, 1.0 / 5.0 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSupportTransfer() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "(((A,B)ab:97,C)abc:57,((D,E)de:10,(F,G)fg:50,(H,I)hi:64)defghi)",
                                                 new NHXParser() )[ 0 ];
            final Phylogeny p2 = factory
                    .create( "(((A:0.1,B:0.3)ab:0.4,C)abc:0.5,((D,E)de,(F,G)fg,(H,I)hi:0.59)defghi)", new NHXParser() )[ 0 ];
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "ab" ) ) >= 0.0 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "abc" ) ) >= 0.0 ) {
                return false;
            }
            support_transfer.moveBranchLengthsToBootstrap( p1 );
            support_transfer.transferSupportValues( p1, p2 );
            if ( p2.getNode( "ab" ).getDistanceToParent() != 0.4 ) {
                return false;
            }
            if ( p2.getNode( "abc" ).getDistanceToParent() != 0.5 ) {
                return false;
            }
            if ( p2.getNode( "hi" ).getDistanceToParent() != 0.59 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "ab" ) ) != 97 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "abc" ) ) != 57 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "de" ) ) != 10 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "fg" ) ) != 50 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "hi" ) ) != 64 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testTaxonomyExtraction() {
        try {
            final PhylogenyNode n0 = PhylogenyNode
                    .createInstanceFromNhxString( "sd_12345678", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n0.getNodeData().isHasTaxonomy() ) {
                return false;
            }
            final PhylogenyNode n1 = PhylogenyNode
                    .createInstanceFromNhxString( "sd_12345x", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n1.getNodeData().isHasTaxonomy() ) {
                System.out.println( n1.toString() );
                return false;
            }
            final PhylogenyNode n2x = PhylogenyNode
                    .createInstanceFromNhxString( "12345", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n2x.getNodeData().isHasTaxonomy() ) {
                return false;
            }
            final PhylogenyNode n3 = PhylogenyNode
                    .createInstanceFromNhxString( "blag_12345", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n3.getNodeData().getTaxonomy().getIdentifier().getValue().equals( "12345" ) ) {
                System.out.println( n3.toString() );
                return false;
            }
            final PhylogenyNode n4 = PhylogenyNode
                    .createInstanceFromNhxString( "blag-12345", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n4.getNodeData().isHasTaxonomy() ) {
                System.out.println( n4.toString() );
                return false;
            }
            final PhylogenyNode n5 = PhylogenyNode
                    .createInstanceFromNhxString( "12345-blag", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n5.getNodeData().isHasTaxonomy() ) {
                System.out.println( n5.toString() );
                return false;
            }
            final PhylogenyNode n6 = PhylogenyNode
                    .createInstanceFromNhxString( "blag-12345-blag", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n6.getNodeData().isHasTaxonomy() ) {
                System.out.println( n6.toString() );
                return false;
            }
            final PhylogenyNode n7 = PhylogenyNode
                    .createInstanceFromNhxString( "blag-12345_blag", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n7.getNodeData().isHasTaxonomy() ) {
                System.out.println( n7.toString() );
                return false;
            }
            final PhylogenyNode n8 = PhylogenyNode
                    .createInstanceFromNhxString( "blag_12345-blag", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n8.getNodeData().getTaxonomy().getIdentifier().getValue().equals( "12345" ) ) {
                System.out.println( n8.toString() );
                return false;
            }
            final PhylogenyNode n9 = PhylogenyNode
                    .createInstanceFromNhxString( "blag_12345/blag", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n9.getNodeData().getTaxonomy().getIdentifier().getValue().equals( "12345" ) ) {
                System.out.println( n9.toString() );
                return false;
            }
            final PhylogenyNode n10x = PhylogenyNode
                    .createInstanceFromNhxString( "blag_12X45-blag", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n10x.getNodeData().isHasTaxonomy() ) {
                System.out.println( n10x.toString() );
                return false;
            }
            final PhylogenyNode n10xx = PhylogenyNode
                    .createInstanceFromNhxString( "blag_1YX45-blag", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( n10xx.getNodeData().isHasTaxonomy() ) {
                System.out.println( n10xx.toString() );
                return false;
            }
            final PhylogenyNode n10 = PhylogenyNode
                    .createInstanceFromNhxString( "blag_9YX45-blag", NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_RELAXED );
            if ( !n10.getNodeData().getTaxonomy().getTaxonomyCode().equals( "9YX45" ) ) {
                System.out.println( n10.toString() );
                return false;
            }
            final PhylogenyNode n11 = PhylogenyNode
                    .createInstanceFromNhxString( "BLAG_Mus_musculus", NHXParser.TAXONOMY_EXTRACTION.AGGRESSIVE );
            if ( !n11.getNodeData().getTaxonomy().getScientificName().equals( "Mus musculus" ) ) {
                System.out.println( n11.toString() );
                return false;
            }
            final PhylogenyNode n12 = PhylogenyNode
                    .createInstanceFromNhxString( "BLAG_Mus_musculus_musculus",
                                                  NHXParser.TAXONOMY_EXTRACTION.AGGRESSIVE );
            if ( !n12.getNodeData().getTaxonomy().getScientificName().equals( "Mus musculus musculus" ) ) {
                System.out.println( n12.toString() );
                return false;
            }
            final PhylogenyNode n13 = PhylogenyNode
                    .createInstanceFromNhxString( "BLAG_Mus_musculus1", NHXParser.TAXONOMY_EXTRACTION.AGGRESSIVE );
            if ( n13.getNodeData().isHasTaxonomy() ) {
                System.out.println( n13.toString() );
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testTreeMethods() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((((A,B)ab,C)abc,D)abcd,E)", new NHXParser() )[ 0 ];
            PhylogenyMethods.collapseSubtreeStructure( t0.getNode( "abcd" ) );
            if ( !t0.toNewHampshireX().equals( "((A,B,C,D)abcd,E)" ) ) {
                System.out.println( t0.toNewHampshireX() );
                return false;
            }
            final Phylogeny t1 = factory.create( "((((A:0.1,B)ab:0.2,C)abc:0.3,D)abcd:0.4,E)", new NHXParser() )[ 0 ];
            PhylogenyMethods.collapseSubtreeStructure( t1.getNode( "abcd" ) );
            if ( !isEqual( t1.getNode( "A" ).getDistanceToParent(), 0.6 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "B" ).getDistanceToParent(), 0.5 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "C" ).getDistanceToParent(), 0.3 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSequenceDbWsTools1() {
        try {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "NP_001025424" );
            Accession acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.REFSEQ.toString() )
                    || !acc.getValue().equals( "NP_001025424" ) ) {
                return false;
            }
            n.setName( "340 0559 -- _NP_001025424_dsfdg15 05" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.REFSEQ.toString() )
                    || !acc.getValue().equals( "NP_001025424" ) ) {
                return false;
            }
            n.setName( "NP_001025424.1" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.REFSEQ.toString() )
                    || !acc.getValue().equals( "NP_001025424" ) ) {
                return false;
            }
            n.setName( "NM_001030253" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.REFSEQ.toString() )
                    || !acc.getValue().equals( "NM_001030253" ) ) {
                return false;
            }
            n.setName( "BCL2_HUMAN" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.UNIPROT.toString() )
                    || !acc.getValue().equals( "BCL2_HUMAN" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "P10415" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.UNIPROT.toString() )
                    || !acc.getValue().equals( "P10415" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( " P10415 " );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.UNIPROT.toString() )
                    || !acc.getValue().equals( "P10415" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "_P10415|" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.UNIPROT.toString() )
                    || !acc.getValue().equals( "P10415" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "AY695820" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.NCBI.toString() )
                    || !acc.getValue().equals( "AY695820" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "_AY695820_" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.NCBI.toString() )
                    || !acc.getValue().equals( "AY695820" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "AAA59452" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.NCBI.toString() )
                    || !acc.getValue().equals( "AAA59452" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "_AAA59452_" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.NCBI.toString() )
                    || !acc.getValue().equals( "AAA59452" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "AAA59452.1" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.NCBI.toString() )
                    || !acc.getValue().equals( "AAA59452.1" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "_AAA59452.1_" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.NCBI.toString() )
                    || !acc.getValue().equals( "AAA59452.1" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "GI:94894583" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.GI.toString() )
                    || !acc.getValue().equals( "94894583" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "gi|71845847|1,4-alpha-glucan branching enzyme [Dechloromonas aromatica RCB]" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.GI.toString() )
                    || !acc.getValue().equals( "71845847" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
            n.setName( "gi|71845847|gb|AAZ45343.1| 1,4-alpha-glucan branching enzyme [Dechloromonas aromatica RCB]" );
            acc = SequenceDbWsTools.obtainSeqAccession( n );
            if ( ( acc == null ) || !acc.getSource().equals( Source.NCBI.toString() )
                    || !acc.getValue().equals( "AAZ45343.1" ) ) {
                System.out.println( acc.toString() );
                return false;
            }
        }
        catch ( final Exception e ) {
            return false;
        }
        return true;
    }

    private static boolean testSequenceDbWsTools2() {
        try {
            final PhylogenyNode n1 = new PhylogenyNode( "NP_001025424" );
            SequenceDbWsTools.obtainSeqInformation( n1 );
            if ( !n1.getNodeData().getSequence().getName().equals( "Bcl2" ) ) {
                return false;
            }
            if ( !n1.getNodeData().getTaxonomy().getScientificName().equals( "Danio rerio" ) ) {
                return false;
            }
            if ( !n1.getNodeData().getSequence().getAccession().getSource().equals( Source.REFSEQ.toString() ) ) {
                return false;
            }
            if ( !n1.getNodeData().getSequence().getAccession().getValue().equals( "NP_001025424" ) ) {
                return false;
            }
            final PhylogenyNode n2 = new PhylogenyNode( "NM_001030253" );
            SequenceDbWsTools.obtainSeqInformation( n2 );
            if ( !n2.getNodeData().getSequence().getName()
                    .equals( "Danio rerio B-cell leukemia/lymphoma 2 (bcl2), mRNA" ) ) {
                return false;
            }
            if ( !n2.getNodeData().getTaxonomy().getScientificName().equals( "Danio rerio" ) ) {
                return false;
            }
            if ( !n2.getNodeData().getSequence().getAccession().getSource().equals( Source.REFSEQ.toString() ) ) {
                return false;
            }
            if ( !n2.getNodeData().getSequence().getAccession().getValue().equals( "NM_001030253" ) ) {
                return false;
            }
            final PhylogenyNode n3 = new PhylogenyNode( "NM_184234.2" );
            SequenceDbWsTools.obtainSeqInformation( n3 );
            if ( !n3.getNodeData().getSequence().getName()
                    .equals( "Homo sapiens RNA binding motif protein 39 (RBM39), transcript variant 1, mRNA" ) ) {
                return false;
            }
            if ( !n3.getNodeData().getTaxonomy().getScientificName().equals( "Homo sapiens" ) ) {
                return false;
            }
            if ( !n3.getNodeData().getSequence().getAccession().getSource().equals( Source.REFSEQ.toString() ) ) {
                return false;
            }
            if ( !n3.getNodeData().getSequence().getAccession().getValue().equals( "NM_184234" ) ) {
                return false;
            }
        }
        catch ( final IOException e ) {
            System.out.println();
            System.out.println( "the following might be due to absence internet connection:" );
            e.printStackTrace( System.out );
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testEbiEntryRetrieval() {
        try {
            final SequenceDatabaseEntry entry = SequenceDbWsTools.obtainEntry( "AAK41263" );
            if ( !entry.getAccession().equals( "AAK41263" ) ) {
                System.out.println( entry.getAccession() );
                return false;
            }
            if ( !entry.getTaxonomyScientificName().equals( "Sulfolobus solfataricus P2" ) ) {
                System.out.println( entry.getTaxonomyScientificName() );
                return false;
            }
            if ( !entry.getSequenceName()
                    .equals( "Sulfolobus solfataricus P2 Glycogen debranching enzyme, hypothetical (treX-like)" ) ) {
                System.out.println( entry.getSequenceName() );
                return false;
            }
            // if ( !entry.getSequenceSymbol().equals( "" ) ) {
            //     System.out.println( entry.getSequenceSymbol() );
            //     return false;
            // }
            if ( !entry.getGeneName().equals( "treX-like" ) ) {
                System.out.println( entry.getGeneName() );
                return false;
            }
            if ( !entry.getTaxonomyIdentifier().equals( "273057" ) ) {
                System.out.println( entry.getTaxonomyIdentifier() );
                return false;
            }
            if ( !entry.getAnnotations().first().getRefValue().equals( "3.2.1.33" ) ) {
                System.out.println( entry.getAnnotations().first().getRefValue() );
                return false;
            }
            if ( !entry.getAnnotations().first().getRefSource().equals( "EC" ) ) {
                System.out.println( entry.getAnnotations().first().getRefSource() );
                return false;
            }
            if ( entry.getCrossReferences().size() != 5 ) {
                return false;
            }
            //
            final SequenceDatabaseEntry entry1 = SequenceDbWsTools.obtainEntry( "ABJ16409" );
            if ( !entry1.getAccession().equals( "ABJ16409" ) ) {
                return false;
            }
            if ( !entry1.getTaxonomyScientificName().equals( "Felis catus" ) ) {
                System.out.println( entry1.getTaxonomyScientificName() );
                return false;
            }
            if ( !entry1.getSequenceName().equals( "Felis catus (domestic cat) partial BCL2" ) ) {
                System.out.println( entry1.getSequenceName() );
                return false;
            }
            if ( !entry1.getTaxonomyIdentifier().equals( "9685" ) ) {
                System.out.println( entry1.getTaxonomyIdentifier() );
                return false;
            }
            if ( !entry1.getGeneName().equals( "BCL2" ) ) {
                System.out.println( entry1.getGeneName() );
                return false;
            }
            if ( entry1.getCrossReferences().size() != 6 ) {
                return false;
            }
            //
            final SequenceDatabaseEntry entry2 = SequenceDbWsTools.obtainEntry( "NM_184234" );
            if ( !entry2.getAccession().equals( "NM_184234" ) ) {
                return false;
            }
            if ( !entry2.getTaxonomyScientificName().equals( "Homo sapiens" ) ) {
                System.out.println( entry2.getTaxonomyScientificName() );
                return false;
            }
            if ( !entry2.getSequenceName()
                    .equals( "Homo sapiens RNA binding motif protein 39 (RBM39), transcript variant 1, mRNA" ) ) {
                System.out.println( entry2.getSequenceName() );
                return false;
            }
            if ( !entry2.getTaxonomyIdentifier().equals( "9606" ) ) {
                System.out.println( entry2.getTaxonomyIdentifier() );
                return false;
            }
            if ( !entry2.getGeneName().equals( "RBM39" ) ) {
                System.out.println( entry2.getGeneName() );
                return false;
            }
            if ( entry2.getCrossReferences().size() != 3 ) {
                return false;
            }
            //
            final SequenceDatabaseEntry entry3 = SequenceDbWsTools.obtainEntry( "HM043801" );
            if ( !entry3.getAccession().equals( "HM043801" ) ) {
                return false;
            }
            if ( !entry3.getTaxonomyScientificName().equals( "Bursaphelenchus xylophilus" ) ) {
                System.out.println( entry3.getTaxonomyScientificName() );
                return false;
            }
            if ( !entry3.getSequenceName().equals( "Bursaphelenchus xylophilus RAF gene, complete cds" ) ) {
                System.out.println( entry3.getSequenceName() );
                return false;
            }
            if ( !entry3.getTaxonomyIdentifier().equals( "6326" ) ) {
                System.out.println( entry3.getTaxonomyIdentifier() );
                return false;
            }
            if ( !entry3.getSequenceSymbol().equals( "RAF" ) ) {
                System.out.println( entry3.getSequenceSymbol() );
                return false;
            }
            if ( !ForesterUtil.isEmpty( entry3.getGeneName() ) ) {
                return false;
            }
            if ( entry3.getCrossReferences().size() != 8 ) {
                return false;
            }
            //
            //
            final SequenceDatabaseEntry entry4 = SequenceDbWsTools.obtainEntry( "AAA36557.1" );
            if ( !entry4.getAccession().equals( "AAA36557" ) ) {
                return false;
            }
            if ( !entry4.getTaxonomyScientificName().equals( "Homo sapiens" ) ) {
                System.out.println( entry4.getTaxonomyScientificName() );
                return false;
            }
            if ( !entry4.getSequenceName().equals( "Homo sapiens (human) ras protein" ) ) {
                System.out.println( entry4.getSequenceName() );
                return false;
            }
            if ( !entry4.getTaxonomyIdentifier().equals( "9606" ) ) {
                System.out.println( entry4.getTaxonomyIdentifier() );
                return false;
            }
            if ( !entry4.getGeneName().equals( "ras" ) ) {
                System.out.println( entry4.getGeneName() );
                return false;
            }
            //   if ( !entry4.getChromosome().equals( "ras" ) ) {
            //     System.out.println( entry4.getChromosome() );
            //     return false;
            // }
            // if ( !entry4.getMap().equals( "ras" ) ) {
            //     System.out.println( entry4.getMap() );
            //     return false;
            // }
            //TODO FIXME gi...
            //
            //TODO fails:
            //            final SequenceDatabaseEntry entry5 = SequenceDbWsTools.obtainEntry( "M30539" );
            //            if ( !entry5.getAccession().equals( "HM043801" ) ) {
            //                return false;
            //            }
            final SequenceDatabaseEntry entry5 = SequenceDbWsTools.obtainEntry( "AAZ45343.1" );
            if ( !entry5.getAccession().equals( "AAZ45343" ) ) {
                return false;
            }
            if ( !entry5.getTaxonomyScientificName().equals( "Dechloromonas aromatica RCB" ) ) {
                System.out.println( entry5.getTaxonomyScientificName() );
                return false;
            }
            if ( !entry5.getSequenceName().equals( "Dechloromonas aromatica RCB 1,4-alpha-glucan branching enzyme" ) ) {
                System.out.println( entry5.getSequenceName() );
                return false;
            }
            if ( !entry5.getTaxonomyIdentifier().equals( "159087" ) ) {
                System.out.println( entry5.getTaxonomyIdentifier() );
                return false;
            }
        }
        catch ( final IOException e ) {
            System.out.println();
            System.out.println( "the following might be due to absence internet connection:" );
            e.printStackTrace( System.out );
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testUniprotEntryRetrieval() {
        try {
            final SequenceDatabaseEntry entry = SequenceDbWsTools.obtainUniProtEntry( "P12345", 200 );
            if ( !entry.getAccession().equals( "P12345" ) ) {
                return false;
            }
            if ( !entry.getTaxonomyScientificName().equals( "Oryctolagus cuniculus" ) ) {
                return false;
            }
            if ( !entry.getSequenceName().equals( "Aspartate aminotransferase, mitochondrial" ) ) {
                return false;
            }
            if ( !entry.getSequenceSymbol().equals( "mAspAT" ) ) {
                return false;
            }
            if ( !entry.getGeneName().equals( "GOT2" ) ) {
                return false;
            }
            if ( !entry.getTaxonomyIdentifier().equals( "9986" ) ) {
                return false;
            }
        }
        catch ( final IOException e ) {
            System.out.println();
            System.out.println( "the following might be due to absence internet connection:" );
            e.printStackTrace( System.out );
            return true;
        }
        catch ( final Exception e ) {
            return false;
        }
        return true;
    }

    private static boolean testUniprotTaxonomySearch() {
        try {
            List<UniProtTaxonomy> results = SequenceDbWsTools.getTaxonomiesFromCommonNameStrict( "starlet sea anemone",
                                                                                                 10 );
            if ( results.size() != 1 ) {
                return false;
            }
            if ( !results.get( 0 ).getCode().equals( "NEMVE" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getCommonName().equalsIgnoreCase( "starlet sea anemone" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getId().equalsIgnoreCase( "45351" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getRank().equalsIgnoreCase( "species" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getScientificName().equals( "Nematostella vectensis" ) ) {
                return false;
            }
            results = null;
            results = SequenceDbWsTools.getTaxonomiesFromScientificNameStrict( "Nematostella vectensis", 10 );
            if ( results.size() != 1 ) {
                return false;
            }
            if ( !results.get( 0 ).getCode().equals( "NEMVE" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getCommonName().equalsIgnoreCase( "starlet sea anemone" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getId().equalsIgnoreCase( "45351" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getRank().equalsIgnoreCase( "species" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getScientificName().equals( "Nematostella vectensis" ) ) {
                return false;
            }
            results = null;
            results = SequenceDbWsTools.getTaxonomiesFromId( "45351", 10 );
            if ( results.size() != 1 ) {
                return false;
            }
            if ( !results.get( 0 ).getCode().equals( "NEMVE" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getCommonName().equalsIgnoreCase( "starlet sea anemone" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getId().equalsIgnoreCase( "45351" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getRank().equalsIgnoreCase( "species" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getScientificName().equals( "Nematostella vectensis" ) ) {
                return false;
            }
            results = null;
            results = SequenceDbWsTools.getTaxonomiesFromTaxonomyCode( "NEMVE", 10 );
            if ( results.size() != 1 ) {
                return false;
            }
            if ( !results.get( 0 ).getCode().equals( "NEMVE" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getCommonName().equalsIgnoreCase( "starlet sea anemone" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getId().equalsIgnoreCase( "45351" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getRank().equalsIgnoreCase( "species" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getScientificName().equals( "Nematostella vectensis" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getLineage().get( 1 ).equals( "Eukaryota" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getLineage().get( 2 ).equals( "Metazoa" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getLineage().get( results.get( 0 ).getLineage().size() - 1 )
                    .equals( "Nematostella vectensis" ) ) {
                System.out.println( results.get( 0 ).getLineage() );
                return false;
            }
            //
            results = null;
            results = SequenceDbWsTools.getTaxonomiesFromScientificNameStrict( "Xenopus tropicalis", 10 );
            if ( results.size() != 1 ) {
                return false;
            }
            if ( !results.get( 0 ).getCode().equals( "XENTR" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getCommonName().equalsIgnoreCase( "Western clawed frog" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getId().equalsIgnoreCase( "8364" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getRank().equalsIgnoreCase( "species" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getScientificName().equals( "Xenopus tropicalis" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getLineage().get( results.get( 0 ).getLineage().size() - 1 )
                    .equals( "Xenopus tropicalis" ) ) {
                System.out.println( results.get( 0 ).getLineage() );
                return false;
            }
            //
            results = null;
            results = SequenceDbWsTools.getTaxonomiesFromId( "8364", 10 );
            if ( results.size() != 1 ) {
                return false;
            }
            if ( !results.get( 0 ).getCode().equals( "XENTR" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getCommonName().equalsIgnoreCase( "Western clawed frog" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getId().equalsIgnoreCase( "8364" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getRank().equalsIgnoreCase( "species" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getScientificName().equals( "Xenopus tropicalis" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getLineage().get( results.get( 0 ).getLineage().size() - 1 )
                    .equals( "Xenopus tropicalis" ) ) {
                System.out.println( results.get( 0 ).getLineage() );
                return false;
            }
            //
            results = null;
            results = SequenceDbWsTools.getTaxonomiesFromTaxonomyCode( "XENTR", 10 );
            if ( results.size() != 1 ) {
                return false;
            }
            if ( !results.get( 0 ).getCode().equals( "XENTR" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getCommonName().equalsIgnoreCase( "Western clawed frog" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getId().equalsIgnoreCase( "8364" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getRank().equalsIgnoreCase( "species" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getScientificName().equals( "Xenopus tropicalis" ) ) {
                return false;
            }
            if ( !results.get( 0 ).getLineage().get( results.get( 0 ).getLineage().size() - 1 )
                    .equals( "Xenopus tropicalis" ) ) {
                System.out.println( results.get( 0 ).getLineage() );
                return false;
            }
        }
        catch ( final IOException e ) {
            System.out.println();
            System.out.println( "the following might be due to absence internet connection:" );
            e.printStackTrace( System.out );
            return true;
        }
        catch ( final Exception e ) {
            return false;
        }
        return true;
    }

    private static boolean testWabiTxSearch() {
        try {
            String result = "";
            result = TxSearch.searchSimple( "nematostella" );
            result = TxSearch.getTxId( "nematostella" );
            if ( !result.equals( "45350" ) ) {
                return false;
            }
            result = TxSearch.getTxName( "45350" );
            if ( !result.equals( "Nematostella" ) ) {
                return false;
            }
            result = TxSearch.getTxId( "nematostella vectensis" );
            if ( !result.equals( "45351" ) ) {
                return false;
            }
            result = TxSearch.getTxName( "45351" );
            if ( !result.equals( "Nematostella vectensis" ) ) {
                return false;
            }
            result = TxSearch.getTxId( "Bacillus subtilis subsp. subtilis str. N170" );
            if ( !result.equals( "536089" ) ) {
                return false;
            }
            result = TxSearch.getTxName( "536089" );
            if ( !result.equals( "Bacillus subtilis subsp. subtilis str. N170" ) ) {
                return false;
            }
            final List<String> queries = new ArrayList<String>();
            queries.add( "Campylobacter coli" );
            queries.add( "Escherichia coli" );
            queries.add( "Arabidopsis" );
            queries.add( "Trichoplax" );
            queries.add( "Samanea saman" );
            queries.add( "Kluyveromyces marxianus" );
            queries.add( "Bacillus subtilis subsp. subtilis str. N170" );
            queries.add( "Bornavirus parrot/PDD/2008" );
            final List<RANKS> ranks = new ArrayList<RANKS>();
            ranks.add( RANKS.SUPERKINGDOM );
            ranks.add( RANKS.KINGDOM );
            ranks.add( RANKS.FAMILY );
            ranks.add( RANKS.GENUS );
            ranks.add( RANKS.TRIBE );
            result = TxSearch.searchLineage( queries, ranks );
            result = TxSearch.searchParam( "Homo sapiens", TAX_NAME_CLASS.ALL, TAX_RANK.SPECIES, 10, true );
            result = TxSearch.searchParam( "Samanea saman", TAX_NAME_CLASS.SCIENTIFIC_NAME, TAX_RANK.ALL, 10, true );
        }
        catch ( final Exception e ) {
            System.out.println();
            System.out.println( "the following might be due to absence internet connection:" );
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
}
