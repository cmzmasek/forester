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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.List;

/**
 *
 * This is to access the Web API for Biology (WABI) at DDBJ.
 * See: http://xml.nig.ac.jp/
 *
 */
public final class RestUtil {

    final static String         LIST_SEPARATOR = "%0A";
    final static String         LINE_SEPARATOR = "\n";
    private final static String BASE_URL       = "http://xml.nig.ac.jp/rest/Invoke";
    private final static String SERVICE        = "service";
    private final static String METHOD         = "method";
    private final static String URL_ENC        = "UTF-8";

    static String encode( final String str ) throws UnsupportedEncodingException {
        return URLEncoder.encode( str.trim(), URL_ENC );
    }

    /**
     * Method for access REST
     * @param query
     * service name method name and parameter for executing rest
     * @return execution result
     * @throws IOException
     */
    public static String getResult( final String query ) throws IOException {
        final URL url = new URL( BASE_URL );
        final URLConnection urlc = url.openConnection();
        urlc.setDoOutput( true );
        urlc.setAllowUserInteraction( false );
        final PrintStream ps = new PrintStream( urlc.getOutputStream() );
        //System.out.println( "query: " + query );
        ps.print( query.trim() );
        ps.close();
        final BufferedReader br = new BufferedReader( new InputStreamReader( urlc.getInputStream() ) );
        final StringBuffer sb = new StringBuffer();
        String line = null;
        while ( ( line = br.readLine() ) != null ) {
            sb.append( line + LINE_SEPARATOR );
        }
        br.close();
        return sb.toString().trim();
    }

    public static String getResult( final String service_name, final String method_name, final String parameters )
            throws IOException {
        final String service = SERVICE + '=' + encode( service_name );
        final String method = METHOD + '=' + encode( method_name );
        return getResult( service + '&' + method + '&' + parameters.trim() );
    }

    static String listAsString( final List<String> l ) throws UnsupportedEncodingException {
        final StringBuffer sb = new StringBuffer();
        for( final String s : l ) {
            if ( sb.length() > 0 ) {
                sb.append( LIST_SEPARATOR );
            }
            sb.append( encode( s ) );
        }
        return sb.toString();
    }
}
