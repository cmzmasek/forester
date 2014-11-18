// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

package org.forester.ws.wabi;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * This is to access the Web API for Biology (WABI) at DDBJ.
 * See: http://xml.nig.ac.jp/
 *
 * Service Description:
 * TXSearch is a retrieval system for a Taxonomy Database which
 * was unified by DDBJ, GenBank and EMBL, which is developed by DDBJ.
 * See: http://xml.nig.ac.jp/wabi/Method?serviceName=TxSearch&mode=methodList
 *
 */
public final class TxSearch {

    private static final String TAXONOMIC_RANK                         = "Taxonomic rank: ";
    private static final String FULL_LINEAGE                           = "Full lineage: ";
    private static final String SEARCH_LINEAGE_QUERY_PARAM_NAME        = "query";
    private static final String SEARCH_LINEAGE_RANKS_PARAM_NAME        = "ranks";
    private static final String SEARCH_LINEAGE_SUPERKINGDOM_PARAM_NAME = "superkingdom";
    private final static String GET_TX_ID_METHOD_NAME                  = "getTxId";
    private final static String GET_TX_NAME_METHOD_NAME                = "getTxName";
    private final static String SEARCH_SIMPLE_METHOD_NAME              = "searchSimple";
    private final static String TX_SEARCH_SERVICE_NAME                 = "TxSearch";
    private final static String TX_NAME_PARAM_NAME                     = "tx_Name";
    private final static String TX_ID_PARAM_NAME                       = "tx_Id";
    private final static String SEARCH_LINEAGE_NAME_METHOD_NAME        = "searchLineage";
    private final static String SEARCH_PARAM_METHOD_NAME               = "searchParam";

    public static String[] getLineage( final String result ) throws IOException {
        String[] lineage = null;
        for( String line : result.split( RestUtil.LINE_SEPARATOR ) ) {
            line = line.trim();
            if ( line.startsWith( FULL_LINEAGE ) ) {
                if ( lineage != null ) {
                    throw new IOException( "search result is not unique" );
                }
                lineage = line.substring( FULL_LINEAGE.length() ).split( ";" );
            }
        }
        return lineage;
    }

    public static String getTaxonomicRank( final String result ) throws IOException {
        String rank = null;
        for( String line : result.split( RestUtil.LINE_SEPARATOR ) ) {
            line = line.trim();
            if ( line.startsWith( TAXONOMIC_RANK ) ) {
                if ( rank != null ) {
                    throw new IOException( "search result is not unique" );
                }
                rank = line.substring( TAXONOMIC_RANK.length() ).trim();
            }
        }
        return rank;
    }

    public static String getTxId( final String tx_name ) throws IOException {
        return RestUtil.getResult( TX_SEARCH_SERVICE_NAME,
                                   GET_TX_ID_METHOD_NAME,
                                   TX_NAME_PARAM_NAME + "=" + RestUtil.encode( tx_name ) ).trim();
    }

    public static String getTxName( final String tx_id ) throws IOException {
        return RestUtil.getResult( TX_SEARCH_SERVICE_NAME,
                                   GET_TX_NAME_METHOD_NAME,
                                   TX_ID_PARAM_NAME + "=" + RestUtil.encode( tx_id ) ).trim();
    }

    public static void main( final String[] args ) throws IOException {
        String result = "";
        try {
            result = searchSimple( "SAMSA" );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
        try {
            result = searchSimple( "nematostella" );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        final String[] lineage = getLineage( result );
        for( final String element : lineage ) {
            System.out.println( element );
        }
        System.out.println( getTaxonomicRank( result ) );
        System.out.println( "---------------" );
        try {
            result = getTxId( "nematostella" );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
        try {
            result = getTxName( "45350" );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
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
        //        ranks.add( RANKS.SUPERKINGDOM );
        //        ranks.add( RANKS.KINGDOM );
        //        ranks.add( RANKS.FAMILY );
        //        ranks.add( RANKS.GENUS );
        ranks.add( RANKS.ALL );
        try {
            result = searchLineage( queries, ranks );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
        try {
            result = searchParam( "Homo sapiens", TAX_NAME_CLASS.ALL, TAX_RANK.SPECIES, 10, true );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
        try {
            result = searchParam( "Samanea saman", TAX_NAME_CLASS.SCIENTIFIC_NAME, TAX_RANK.ALL, 10, true );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
        try {
            result = searchParam( "cow", TAX_NAME_CLASS.COMMON_NAME, TAX_RANK.ALL, 10, true );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
        try {
            result = searchParam( "Helicogloea lagerheimii", TAX_NAME_CLASS.SCIENTIFIC_NAME, TAX_RANK.ALL, 10, true );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
        try {
            result = searchParam( "Cronartium ribicola", TAX_NAME_CLASS.SCIENTIFIC_NAME, TAX_RANK.ALL, 10, true );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
        try {
            result = searchParam( "Peridermium harknessii", TAX_NAME_CLASS.SCIENTIFIC_NAME, TAX_RANK.ALL, 10, true );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
        System.out.println( "---------------" );
        try {
            result = searchParam( "Eukaryota", TAX_NAME_CLASS.SCIENTIFIC_NAME, TAX_RANK.ALL, 10, true );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println( result );
    }

    private static String ranksAsString( final List<RANKS> l ) throws UnsupportedEncodingException {
        final StringBuffer sb = new StringBuffer();
        for( final RANKS r : l ) {
            if ( sb.length() > 0 ) {
                sb.append( RestUtil.LIST_SEPARATOR );
            }
            sb.append( RestUtil.encode( r.toString() ) );
        }
        return sb.toString();
    }

    public static String searchLineage( final List<String> queries, final List<RANKS> ranks ) throws IOException {
        return searchLineage( queries, ranks, "" );
    }

    public static String searchLineage( final List<String> queries, final List<RANKS> ranks, final String superkingdom )
            throws IOException {
        return RestUtil.getResult( TX_SEARCH_SERVICE_NAME,
                                   SEARCH_LINEAGE_NAME_METHOD_NAME,
                                   SEARCH_LINEAGE_QUERY_PARAM_NAME + "=" + RestUtil.listAsString( queries ) + "&"
                                           + SEARCH_LINEAGE_RANKS_PARAM_NAME + "=" + ranksAsString( ranks ) + "&"
                                           + SEARCH_LINEAGE_SUPERKINGDOM_PARAM_NAME + "="
                                           + RestUtil.encode( superkingdom ) ).trim();
    }

    public static String searchParam( final String tx_name,
                                      final TAX_NAME_CLASS tx_name_class,
                                      final TAX_RANK tx_rank,
                                      int tx_rmax,
                                      final boolean as_scientific_name ) throws IOException {
        String as_scientific_name_str = "no";
        if ( as_scientific_name ) {
            as_scientific_name_str = "yes";
        }
        if ( tx_rmax < 1 ) {
            tx_rmax = 1;
        }
        return RestUtil.getResult( TX_SEARCH_SERVICE_NAME,
                                   SEARCH_PARAM_METHOD_NAME,
                                   TX_NAME_PARAM_NAME + "=" + RestUtil.encode( tx_name ) + "&tx_Clas="
                                           + RestUtil.encode( tx_name_class.toString() ) + "&tx_Rank="
                                           + RestUtil.encode( tx_rank.toString() ) + "&tx_Rmax=" + tx_rmax
                                           + "&tx_Dcls=" + as_scientific_name_str ).trim();
    }

    public static String searchSimple( final String tx_name ) throws IOException {
        return RestUtil.getResult( TX_SEARCH_SERVICE_NAME,
                                   SEARCH_SIMPLE_METHOD_NAME,
                                   TX_NAME_PARAM_NAME + "=" + RestUtil.encode( tx_name ) ).trim();
    }

    public enum RANKS {
        ALL( "all" ),
        SUPERKINGDOM( "superkingdom" ),
        KINGDOM( "kingdom" ),
        SUBKINGDOM( "subkingdom" ),
        SUPERPHYLUM( "superphylum" ),
        PHYLUM( "phylum" ),
        SUBPHYLUM( "subphylum" ),
        SUPERCLASS( "superclass" ),
        CLASS( "class" ),
        SUBCLASS( "subclass" ),
        INFRACLASS( "infraclass" ),
        SUPERORDER( "superorder" ),
        ORDER( "order" ),
        SUBORDER( "suborder" ),
        INFRAORDER( "infraorder" ),
        PARVORDER( "parvorder" ),
        SUPERFAMILY( "superfamily" ),
        FAMILY( "family" ),
        SUBFAMILY( "subfamily" ),
        TRIBE( "tribe" ),
        SUBTRIBE( "subtribe" ),
        GENUS( "genus" ),
        SPECIES( "species" );

        private final String _str;

        private RANKS( final String name ) {
            _str = name;
        }

        @Override
        public String toString() {
            return _str;
        }
    }

    public enum TAX_NAME_CLASS {
        ALL( "all" ),
        SCIENTIFIC_NAME( "scientific name" ),
        PREFFERED_COMMON_NAME( "preferred common name" ),
        COMMON_NAME( "common name" ),
        SYNONYM( "synonym" );

        private final String _str;

        private TAX_NAME_CLASS( final String name ) {
            _str = name;
        }

        @Override
        public String toString() {
            return _str;
        }
    }

    public enum TAX_RANK {
        ALL( "All" ),
        NO_RANK( "no rank" ),
        SUPERKINGDOM( "superkingdom" ),
        KINGDOM( "kingdom" ),
        SUBKINGDOM( "subkingdom" ),
        SUPERPHYLUM( "superphylum" ),
        PHYLUM( "phylum" ),
        SUBPHYLUM( "subphylum" ),
        SUPERCLASS( "superclass" ),
        CLASS( "class" ),
        SUBCLASS( "subclass" ),
        INFRACLASS( "infraclass" ),
        SUPERORDER( "superorder" ),
        ORDER( "order" ),
        SUBORDER( "suborder" ),
        INFRAORDER( "infraorder" ),
        PARVORDER( "parvorder" ),
        SUPERFAMILY( "superfamily" ),
        FAMILY( "family" ),
        SUBFAMILY( "subfamily" ),
        TRIBE( "tribe" ),
        SUBTRIBE( "subtribe" ),
        GENUS( "genus" ),
        SUBGENUS( "subgenus" ),
        SPECIES_GROUP( "species group" ),
        SPECIES_SUBGROUP( "species subgroup" ),
        SPECIES( "species" ),
        SUBSPECIES( "subspecies" ),
        VARIETAS( "varietas" ),
        FORMA( "forma" );

        private final String _str;

        private TAX_RANK( final String name ) {
            _str = name;
        }

        @Override
        public String toString() {
            return _str;
        }
    }
}
