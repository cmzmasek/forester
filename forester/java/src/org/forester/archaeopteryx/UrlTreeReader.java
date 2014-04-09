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

package org.forester.archaeopteryx;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Date;

import javax.swing.JOptionPane;

import org.forester.archaeopteryx.webservices.PhylogeniesWebserviceClient;
import org.forester.archaeopteryx.webservices.WebserviceUtil;
import org.forester.archaeopteryx.webservices.WebservicesManager;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.tol.TolParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterUtil;

public class UrlTreeReader implements Runnable {

    private final MainFrame _main_frame;
    private final int       _webservice_client_index;

    UrlTreeReader( final MainFrame mf, final int webservice_client_index ) {
        _main_frame = mf;
        _webservice_client_index = webservice_client_index;
    }

    @Override
    public void run() {
        readPhylogeniesFromWebservice();
    }

    synchronized void readPhylogeniesFromWebservice() {
        final long start_time = new Date().getTime();
        URL url = null;
        Phylogeny[] trees = null;
        final WebservicesManager webservices_manager = WebservicesManager.getInstance();
        final PhylogeniesWebserviceClient client = webservices_manager
                .getAvailablePhylogeniesWebserviceClient( _webservice_client_index );
        String identifier = JOptionPane.showInputDialog( _main_frame, client.getInstructions() + "\n(Reference: "
                + client.getReference() + ")", client.getDescription(), JOptionPane.QUESTION_MESSAGE );
        if ( ( identifier != null ) && ( identifier.trim().length() > 0 ) ) {
            identifier = identifier.trim();
            if ( client.isQueryInteger() ) {
                identifier = identifier.replaceAll( "^\\D+", "" );
                int id = -1;
                try {
                    id = Integer.parseInt( identifier );
                }
                catch ( final NumberFormatException e ) {
                    id = -1;
                }
                if ( id < 1 ) {
                    JOptionPane.showMessageDialog( _main_frame,
                                                   "Identifier is expected to be a number",
                                                   "Can not open URL",
                                                   JOptionPane.ERROR_MESSAGE );
                    return;
                }
                identifier = id + "";
            }
            boolean exception = false;
            try {
                String url_str = client.getUrl();
                url_str = url_str.replaceFirst( PhylogeniesWebserviceClient.QUERY_PLACEHOLDER, identifier );
                url = new URL( url_str );
                PhylogenyParser parser = null;
                switch ( client.getReturnFormat() ) {
                    case TOL_XML_RESPONSE:
                        parser = new TolParser();
                        break;
                    case NEXUS:
                        parser = new NexusPhylogeniesParser();
                        ( ( NexusPhylogeniesParser ) parser ).setReplaceUnderscores( true );
                        break;
                    case TREEBASE_TREE:
                        parser = new NexusPhylogeniesParser();
                        ( ( NexusPhylogeniesParser ) parser ).setReplaceUnderscores( true );
                        ( ( NexusPhylogeniesParser ) parser ).setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.NO );
                        break;
                    case TREEBASE_STUDY:
                        parser = new NexusPhylogeniesParser();
                        ( ( NexusPhylogeniesParser ) parser ).setReplaceUnderscores( true );
                        ( ( NexusPhylogeniesParser ) parser ).setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.NO );
                        break;
                    case NH:
                        parser = new NHXParser();
                        ( ( NHXParser ) parser ).setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.NO );
                        ( ( NHXParser ) parser ).setReplaceUnderscores( true );
                        ( ( NHXParser ) parser ).setGuessRootedness( true );
                        break;
                    case NH_EXTRACT_TAXONOMY:
                        parser = new NHXParser();
                        ( ( NHXParser ) parser ).setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.AGGRESSIVE );
                        ( ( NHXParser ) parser ).setReplaceUnderscores( false );
                        ( ( NHXParser ) parser ).setGuessRootedness( true );
                        break;
                    case PFAM:
                        parser = new NHXParser();
                        ( ( NHXParser ) parser )
                                .setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.PFAM_STYLE_STRICT );
                        ( ( NHXParser ) parser ).setReplaceUnderscores( false );
                        ( ( NHXParser ) parser ).setGuessRootedness( true );
                        break;
                    case NHX:
                        parser = new NHXParser();
                        ( ( NHXParser ) parser ).setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.NO );
                        ( ( NHXParser ) parser ).setReplaceUnderscores( false );
                        ( ( NHXParser ) parser ).setGuessRootedness( true );
                        break;
                    case PHYLOXML:
                        parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
                        break;
                    default:
                        throw new IllegalArgumentException( "unknown format: " + client.getReturnFormat() );
                }
                if ( _main_frame.getMainPanel().getCurrentTreePanel() != null ) {
                    _main_frame.getMainPanel().getCurrentTreePanel().setWaitCursor();
                }
                else {
                    _main_frame.getMainPanel().setWaitCursor();
                }
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                trees = factory.create( url.openStream(), parser );
            }
            catch ( final MalformedURLException e ) {
                exception = true;
                JOptionPane.showMessageDialog( _main_frame,
                                               "Malformed URL: " + url + "\n" + e.getLocalizedMessage(),
                                               "Malformed URL",
                                               JOptionPane.ERROR_MESSAGE );
            }
            catch ( final IOException e ) {
                exception = true;
                JOptionPane.showMessageDialog( _main_frame,
                                               "Could not read from " + url + "\n" + e.getLocalizedMessage(),
                                               "Failed to read tree from " + client.getName() + " for " + identifier,
                                               JOptionPane.ERROR_MESSAGE );
            }
            catch ( final NumberFormatException e ) {
                exception = true;
                JOptionPane.showMessageDialog( _main_frame,
                                               "Could not read from " + url + "\n" + e.getLocalizedMessage(),
                                               "Failed to read tree from " + client.getName() + " for " + identifier,
                                               JOptionPane.ERROR_MESSAGE );
            }
            catch ( final Exception e ) {
                exception = true;
                e.printStackTrace();
                JOptionPane.showMessageDialog( _main_frame,
                                               e.getLocalizedMessage(),
                                               "Unexpected Exception",
                                               JOptionPane.ERROR_MESSAGE );
            }
            finally {
                if ( _main_frame.getCurrentTreePanel() != null ) {
                    _main_frame.getCurrentTreePanel().setArrowCursor();
                }
                else {
                    _main_frame.getMainPanel().setArrowCursor();
                }
            }
            if ( ( trees != null ) && ( trees.length > 0 ) ) {
                for( final Phylogeny phylogeny : trees ) {
                    if ( !phylogeny.isEmpty() ) {
                        if ( client.getName().equals( WebserviceUtil.TREE_FAM_NAME ) ) {
                            phylogeny.setRerootable( false );
                            phylogeny.setRooted( true );
                        }
                        if ( client.getProcessingInstructions() != null ) {
                            try {
                                WebserviceUtil.processInstructions( client, phylogeny );
                            }
                            catch ( final PhyloXmlDataFormatException e ) {
                                JOptionPane.showMessageDialog( _main_frame,
                                                               "Error:\n" + e.getLocalizedMessage(),
                                                               "Error",
                                                               JOptionPane.ERROR_MESSAGE );
                            }
                        }
                        if ( client.getNodeField() != null ) {
                            try {
                                PhylogenyMethods.transferNodeNameToField( phylogeny, client.getNodeField(), false );
                            }
                            catch ( final PhyloXmlDataFormatException e ) {
                                JOptionPane.showMessageDialog( _main_frame,
                                                               "Error:\n" + e.getLocalizedMessage(),
                                                               "Error",
                                                               JOptionPane.ERROR_MESSAGE );
                            }
                        }
                        phylogeny.setIdentifier( new Identifier( identifier, client.getName() ) );
                        _main_frame.getJMenuBar().remove( _main_frame.getHelpMenu() );
                        _main_frame.getMenuBarOfMainFrame().add( _main_frame.getHelpMenu() );
                        _main_frame.getMainPanel().addPhylogenyInNewTab( phylogeny,
                                                                         _main_frame.getConfiguration(),
                                                                         new File( url.getFile() ).getName(),
                                                                         url.toString() );
                        String my_name_for_file = "";
                        if ( !ForesterUtil.isEmpty( phylogeny.getName() ) ) {
                            my_name_for_file = new String( phylogeny.getName() ).replaceAll( " ", "_" );
                        }
                        else if ( phylogeny.getIdentifier() != null ) {
                            final StringBuffer sb = new StringBuffer();
                            if ( !ForesterUtil.isEmpty( phylogeny.getIdentifier().getProvider() ) ) {
                                sb.append( phylogeny.getIdentifier().getProvider() );
                                sb.append( "_" );
                            }
                            sb.append( phylogeny.getIdentifier().getValue() );
                            my_name_for_file = new String( sb.toString().replaceAll( " ", "_" ) );
                        }
                        _main_frame.getMainPanel().getCurrentTreePanel().setTreeFile( new File( my_name_for_file ) );
                        AptxUtil.lookAtSomeTreePropertiesForAptxControlSettings( phylogeny, _main_frame.getMainPanel()
                                .getControlPanel(), _main_frame.getConfiguration() );
                        _main_frame.getMainPanel().getControlPanel().showWhole();
                    }
                }
            }
            else if ( !exception ) {
                JOptionPane.showMessageDialog( null, ForesterUtil.wordWrap( "Failed to read in tree(s) from [" + url
                        + "]", 80 ), "Error", JOptionPane.ERROR_MESSAGE );
            }
            _main_frame.getContentPane().repaint();
            if ( ( trees != null ) && ( trees.length > 0 ) ) {
                try {
                    JOptionPane.showMessageDialog( null,
                                                   ForesterUtil.wordWrap( "Successfully read in " + trees.length
                                                           + " tree(s) from [" + url + "]", 80 ),
                                                   "Success",
                                                   JOptionPane.INFORMATION_MESSAGE );
                }
                catch ( final Exception e ) {
                    // Not important if this fails, do nothing.
                }
                _main_frame.getContentPane().repaint();
            }
        }
        _main_frame.activateSaveAllIfNeeded();
        System.gc();
    }
}
