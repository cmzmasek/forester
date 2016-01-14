
package org.forester.ws.hmmer;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLEncoder;

public class Test {

    public static void main( final String[] args ) {
        try {
            final URL url = new URL( "http://hmmer.janelia.org/search/hmmscan" );
            final HttpURLConnection connection = ( HttpURLConnection ) url.openConnection();
            connection.setDoOutput( true );
            connection.setDoInput( true );
            connection.setInstanceFollowRedirects( false );
            connection.setRequestMethod( "POST" );
            connection.setRequestProperty( "Content-Type", "application/x-www-form-urlencoded" );
            connection.setRequestProperty( "Accept", "application/json" );
            //Add the database and the sequence. Add more options as you wish!
            final String urlParameters = "hmmdb=" + URLEncoder.encode( "pfam", "UTF-8" ) + "&seq="
                    + ">seq\nEMGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPV"
                    + "NSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEG"
                    + "RVYHYRINTASDGKLYVSSESRFNTLAELVHHHSTVADGLITTLHYPAP";
            connection.setRequestProperty( "Content-Length", "" + Integer.toString( urlParameters.getBytes().length ) );
            //Send request
            final DataOutputStream wr = new DataOutputStream( connection.getOutputStream() );
            wr.writeBytes( urlParameters );
            wr.flush();
            wr.close();
            //Now get the redirect URL
            final URL respUrl = new URL( connection.getHeaderField( "Location" ) );
            final HttpURLConnection connection2 = ( HttpURLConnection ) respUrl.openConnection();
            connection2.setRequestMethod( "GET" );
            connection2.setRequestProperty( "Accept", "text/x-yaml" );
            //Get the response and print it to the screen
            final BufferedReader in = new BufferedReader( new InputStreamReader( connection2.getInputStream() ) );
            String inputLine;
            while ( ( inputLine = in.readLine() ) != null ) {
                System.out.println( inputLine );
            }
            in.close();
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
    }
}
