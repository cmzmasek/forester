
package org.forester.ws.seqdb;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.phylogeny.data.Accession;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;

public final class UniprotRetrieve {

    
    private final boolean        _verbose;
    private final static String  BASE_URL = "http://www.uniprot.org/uploadlists/";
    private final static boolean DEBUG    = false;

    public UniprotRetrieve( final boolean verbose ) {
        _verbose = verbose;
    }

    public UniprotRetrieve() {
        _verbose = false;
    }

    /*
     * See: https://www.uniprot.org/help/api_idmapping
     *
     */
    public final SortedMap<String, UniprotData> retrieve( final Set<String> queries ) throws IOException {
        final List<String> queries_list = new ArrayList<>();
        for( final String q : queries ) {
            queries_list.add( q );
        }
        return retrieve( queries_list );
    }

    public final SortedMap<String, UniprotData> retrieve( final List<String> queries ) throws IOException {
        final List<String> ncbi_queries = new ArrayList<>();
        final List<String> refseq_queries = new ArrayList<>();
        for( final String query : queries ) {
            final Accession acc = SequenceAccessionTools.parseAccessorFromString( query );
            if ( acc != null ) {
                if ( acc.getSource().equals( "ncbi" ) ) {
                    ncbi_queries.add( query );
                }
                else if ( acc.getSource().equals( "refseq" ) ) {
                    refseq_queries.add( query );
                }
                else if ( acc.getSource().equals( "ensembl" ) ) {
                    System.out.println( "NOTE: query \"" + query + "\" is Ensembl -- ignored" );
                }
                else if ( acc.getSource().equals( "uniprot" ) ) {
                    System.out.println( "NOTE: query \"" + query + "\" is UniProt -- ignored" );
                }
                else if ( acc.getSource().equals( SequenceAccessionTools.VIPR_SOURCE ) ) {
                    System.out.println( "NOTE: query \"" + query + "\" is ViPR -- ignored" );
                }
                else {
                    System.out.println( "NOTE: query \"" + query + "\" is of unknown source -- ignored" );
                }
            }
            else {
                if ( _verbose ) {
                    System.out.println( "NOTE: query \"" + query + "\" is of unknown source -- ignored" );
                }
            }
        }
        final SortedMap<String, UniprotData> m = new TreeMap<>();
        runQuery( "P_REFSEQ_AC", "ACC", refseq_queries, m );
        runQuery( "EMBL", "ACC", ncbi_queries, m );
        return m;
    }

    private final void runQuery( final String from,
                                 final String to,
                                 final List<String> queries,
                                 final SortedMap<String, UniprotData> m )
            throws IOException {
        if ( ForesterUtil.isEmpty( queries )) {
            return;
        }
        final StringBuilder intitial_url_str = new StringBuilder( BASE_URL + "?from=" + from + "&to=" + to
                + "&format=tab&query=" );
        boolean first = true;
        for( final String query : queries ) {
            if ( first ) {
                first = false;
            }
            else {
                intitial_url_str.append( "+" );
            }
            intitial_url_str.append( query );
        }
        String url_str = intitial_url_str.toString();
        if ( DEBUG ) {
            System.out.println( "\n\n" + url_str + "\n" );
        }
        URL resource_url;
        URL base;
        URL next;
        HttpURLConnection conn;
        String location;
        while ( true ) {
            resource_url = new URL( url_str );
            conn = ( HttpURLConnection ) resource_url.openConnection();
            conn.setConnectTimeout( 50000 );
            conn.setReadTimeout( 50000 );
            conn.setInstanceFollowRedirects( false ); // make the logic below easier to detect redirections
            conn.setRequestProperty( "User-Agent", "Mozilla/5.0..." );
            switch ( conn.getResponseCode() ) {
                case HttpURLConnection.HTTP_MOVED_PERM:
                case HttpURLConnection.HTTP_MOVED_TEMP:
                    location = conn.getHeaderField( "Location" );
                    location = URLDecoder.decode( location, "UTF-8" );
                    location = location.replaceAll( "\\s+", "%20" );
                    base = new URL( url_str );
                    next = new URL( base, location ); // deal with relative URLs
                    url_str = next.toExternalForm();
                    continue;
            }
            break;
        }
        final BufferedReader in = new BufferedReader( new InputStreamReader( conn.getInputStream() ) );
        String line;
        int counter = 0;
        while ( ( line = in.readLine() ) != null ) {
            if ( DEBUG ) {
                System.out.println( line );
            }
            if ( !line.startsWith( "yourlist" ) ) {
                final UniprotData d = new UniprotData( line, conn.getURL().toString() );
                m.put( d.getId(), d );
                if ( _verbose ) {
                    System.out.println( ( counter++ ) + ": " + d.getId() + " -> " + d.getEntryName() );
                }
            }
        }
        in.close();
    }

    public static void main( final String[] args ) throws Exception {
        final List<String> queries = new ArrayList<>();
        queries.add( "YP_009256203.1" );
        queries.add( "YP_009199248.1" );
        queries.add( "AAX85683.1" );
        queries.add( "AGT20977.1" );
        queries.add( "ABG89282.1" );
        queries.add( "ACN89752.1" );
        queries.add( "YP_006908645.1" );
        queries.add( "ACN89728.1" );
        queries.add( "YP_001718613.1" );
        queries.add( "AIV41800.1" );
        final UniprotRetrieve ret = new UniprotRetrieve( true );
        final SortedMap<String, UniprotData> m = ret.retrieve( queries );
        final Iterator<Entry<String, UniprotData>> it = m.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<String, UniprotData> pair = it.next();
            System.out.println( pair.getKey() + " => " + pair.getValue().getProteinNames() );
        }
    }
}
