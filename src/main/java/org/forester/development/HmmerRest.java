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

package org.forester.development;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.URL;
import java.net.URLConnection;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class HmmerRest {

    final static String         LIST_SEPARATOR = "%0A";
    final static String         LINE_SEPARATOR = "\n";
    private final static String BASE_URL       = "http://pfam.sanger.ac.uk/search/sequence";

    public static void main( final String[] args ) {
        final String seq = "MASTENNEKDNFMRDTASRSKKSRRRSLWIAAGAVPTAIALSLSLASPA"
                + "AVAQSSFGSSDIIDSGVLDSITRGLTDYLTPRDEALPAGEVTYPAIEGLP"
                + "AGVRVNSAEYVTSHHVVLSIQSAAMPERPIKVQLLLPRDWYSSPDRDFPE"
                + "IWALDGLRAIEKQSGWTIETNIEQFFADKNAIVVLPVGGESSFYTDWNEP"
                + "NNGKNYQWETFLTEELAPILDKGFRSNGERAITGISMGGTAAVNIATHNP"
                + "EMFNFVGSFSGYLDTTSNGMPAAIGAALADAGGYNVNAMWGPAGSERWLE"
                + "NDPKRNVDQLRGKQVYVSAGSGADDYGQDGSVATGPANAAGVGLELISRM"
                + "TSQTFVDAANGAGVNVIANFRPSGVHAWPYWQFEMTQAWPYMADSLGMSR"
                + "EDRGADCVALGAIADATADGSLGSCLNNEYLVANGVGRAQDFTNGRAYWS"
                + "PNTGAFGLFGRINARYSELGGPDSWLGFPKTRELSTPDGRGRYVHFENGS"
                + "IYWSAATGPWEIPGDMFTAWGTQGYEAGGLGYPVGPAKDFNGGLAQEFQG"
                + "GYVLRTPQNRAYWVRGAISAKYMEPGVATTLGFPTGNERLIPGGAFQEFT"
                + "NGNIYWSASTGAHYILRGGIFDAWGAKGYEQGEYGWPTTDQTSIAAGGET" + "ITFQNGTIRQVNGRIEESR";
        final String query = "seq=" + seq + "" + "&" + "output=xml";
        String result = "";
        try {
            result = getResult( query );
        }
        catch ( final IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        System.out.println( result );
        final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        Document dom = null;
        try {
            //Using factory get an instance of document builder
            final DocumentBuilder db = dbf.newDocumentBuilder();
            //parse using builder to get DOM representation of the XML file
            dom = db.parse( new ByteArrayInputStream( result.getBytes() ) );
        }
        catch ( final ParserConfigurationException pce ) {
            pce.printStackTrace();
        }
        catch ( final SAXException se ) {
            se.printStackTrace();
        }
        catch ( final IOException ioe ) {
            ioe.printStackTrace();
        }
        final Element docEle = dom.getDocumentElement();
        final NodeList nl = docEle.getElementsByTagName( "job" );
        String result_url = "";
        for( int i = 0; i < nl.getLength(); i++ ) {
            //System.out.println( nl.item( i ) );
            result_url = getTextValue( ( Element ) nl.item( i ), "result_url" );
        }
        System.out.println( "result url = " + result_url );
        //gettin the result....
        try {
            //do what you want to do before sleeping
            Thread.sleep( 5000 );//sleep for x ms
            //do what you want to do after sleeptig
        }
        catch ( final InterruptedException ie ) {
            ie.printStackTrace();
        }
        try {
            result = getResult( result_url, "" );
        }
        catch ( final IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        System.out.println( result );
    }

    private static String getTextValue( final Element ele, final String tagName ) {
        String textVal = null;
        final NodeList nl = ele.getElementsByTagName( tagName );
        if ( ( nl != null ) && ( nl.getLength() > 0 ) ) {
            final Element el = ( Element ) nl.item( 0 );
            textVal = el.getFirstChild().getNodeValue();
        }
        return textVal;
    }

    public static String getResult( final String base_url, final String query ) throws IOException {
        System.out.println( query );
        final URL url = new URL( base_url );
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

    public static String getResult( final String query ) throws IOException {
        System.out.println( query );
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
}
