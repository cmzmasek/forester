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
// WWW: www.phylosoft.org/forester

package org.forester.ws.uniprot;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;

import org.forester.util.ForesterUtil;

public final class UniProtWsTools {

    public final static String   BASE_URL = "http://www.uniprot.org/";
    private final static String  URL_ENC  = "UTF-8";
    private final static boolean DEBUG    = false;

    synchronized private static String encode( final String str ) throws UnsupportedEncodingException {
        return URLEncoder.encode( str.trim(), URL_ENC );
    }

    synchronized public static List<UniProtTaxonomy> getTaxonomiesFromCommonName( final String cn,
                                                                                  final int max_taxonomies_return )
            throws IOException {
        final List<String> result = getTaxonomyStringFromCommonName( cn, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

    synchronized public static List<UniProtTaxonomy> getTaxonomiesFromCommonNameStrict( final String cn,
                                                                                        final int max_taxonomies_return )
            throws IOException {
        final List<UniProtTaxonomy> taxonomies = getTaxonomiesFromCommonName( cn, max_taxonomies_return );
        if ( ( taxonomies != null ) && ( taxonomies.size() > 0 ) ) {
            final List<UniProtTaxonomy> filtered_taxonomies = new ArrayList<UniProtTaxonomy>();
            for( final UniProtTaxonomy taxonomy : taxonomies ) {
                if ( taxonomy.getCommonName().equalsIgnoreCase( cn ) ) {
                    filtered_taxonomies.add( taxonomy );
                }
            }
            return filtered_taxonomies;
        }
        return null;
    }

    synchronized public static List<UniProtTaxonomy> getTaxonomiesFromId( final String id,
                                                                          final int max_taxonomies_return )
            throws IOException {
        final List<String> result = getTaxonomyStringFromId( id, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

    synchronized public static List<UniProtTaxonomy> getTaxonomiesFromScientificName( final String sn,
                                                                                      final int max_taxonomies_return )
            throws IOException {
        // Hack!  Craniata? .. 
        if ( sn.equals( "Drosophila" ) ) {
            return hack( UniProtTaxonomy.DROSOPHILA_GENUS );
        }
        else if ( sn.equals( "Xenopus" ) ) {
            return hack( UniProtTaxonomy.XENOPUS_GENUS );
        }
        final List<String> result = getTaxonomyStringFromScientificName( sn, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

    /**
     * Does not return "sub-types".
     * For example, for "Mus musculus" only returns "Mus musculus"
     * and not "Mus musculus", "Mus musculus bactrianus", ...
     * 
     */
    synchronized public static List<UniProtTaxonomy> getTaxonomiesFromScientificNameStrict( final String sn,
                                                                                            final int max_taxonomies_return )
            throws IOException {
        final List<UniProtTaxonomy> taxonomies = getTaxonomiesFromScientificName( sn, max_taxonomies_return );
        if ( ( taxonomies != null ) && ( taxonomies.size() > 0 ) ) {
            final List<UniProtTaxonomy> filtered_taxonomies = new ArrayList<UniProtTaxonomy>();
            for( final UniProtTaxonomy taxonomy : taxonomies ) {
                if ( taxonomy.getScientificName().equalsIgnoreCase( sn ) ) {
                    filtered_taxonomies.add( taxonomy );
                }
            }
            return filtered_taxonomies;
        }
        return null;
    }

    synchronized public static List<UniProtTaxonomy> getTaxonomiesFromTaxonomyCode( final String code,
                                                                                    final int max_taxonomies_return )
            throws IOException {
        String my_code = new String( code );
        // Hacks!
        if ( my_code.equals( "FUGRU" ) ) {
            my_code = "TAKRU";
        }
        else if ( my_code.equals( "CAP" ) ) {
            return hack( UniProtTaxonomy.CAPITELLA_TELATA_SPECIES );
        }
        final List<String> result = getTaxonomyStringFromTaxonomyCode( my_code, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

    synchronized private static List<String> getTaxonomyStringFromCommonName( final String cn,
                                                                              final int max_lines_to_return )
            throws IOException {
        return queryUniprot( "taxonomy/?query=common%3a%22" + encode( cn ) + "%22&format=tab", max_lines_to_return );
    }

    synchronized private static List<String> getTaxonomyStringFromId( final String id, final int max_lines_to_return )
            throws IOException {
        return queryUniprot( "taxonomy/?query=id%3a%22" + encode( id ) + "%22&format=tab", max_lines_to_return );
    }

    synchronized private static List<String> getTaxonomyStringFromScientificName( final String sn,
                                                                                  final int max_lines_to_return )
            throws IOException {
        return queryUniprot( "taxonomy/?query=scientific%3a%22" + encode( sn ) + "%22&format=tab", max_lines_to_return );
    }

    synchronized private static List<String> getTaxonomyStringFromTaxonomyCode( final String code,
                                                                                final int max_lines_to_return )
            throws IOException {
        return queryUniprot( "taxonomy/?query=mnemonic%3a%22" + encode( code ) + "%22&format=tab", max_lines_to_return );
    }

    synchronized private static List<UniProtTaxonomy> hack( final UniProtTaxonomy tax ) {
        final List<UniProtTaxonomy> l = new ArrayList<UniProtTaxonomy>();
        l.add( tax );
        return l;
    }

    synchronized private static List<UniProtTaxonomy> parseUniProtTaxonomy( final List<String> result )
            throws IOException {
        final List<UniProtTaxonomy> taxonomies = new ArrayList<UniProtTaxonomy>();
        for( final String line : result ) {
            if ( ForesterUtil.isEmpty( line ) ) {
                // Ignore empty lines.
            }
            else if ( line.startsWith( "Taxon" ) ) {
                //TODO next the check format FIXME
            }
            else {
                if ( line.split( "\t" ).length > 4 ) {
                    taxonomies.add( new UniProtTaxonomy( line ) );
                }
            }
        }
        return taxonomies;
    }

    synchronized public static List<String> queryUniprot( final String query, int max_lines_to_return )
            throws IOException {
        if ( ForesterUtil.isEmpty( query ) ) {
            throw new IllegalArgumentException( "illegal attempt to use empty query " );
        }
        if ( max_lines_to_return < 1 ) {
            max_lines_to_return = 1;
        }
        final URL url = new URL( BASE_URL + query );
        if ( DEBUG ) {
            System.out.println( "url: " + url.toString() );
        }
        final URLConnection urlc = url.openConnection();
        final BufferedReader in = new BufferedReader( new InputStreamReader( urlc.getInputStream() ) );
        String line;
        final List<String> result = new ArrayList<String>();
        while ( ( line = in.readLine() ) != null ) {
            result.add( line );
            if ( result.size() > max_lines_to_return ) {
                break;
            }
        }
        in.close();
        return result;
    }
}
