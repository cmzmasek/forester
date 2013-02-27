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

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.KeyboardFocusManager;
import java.io.File;
import java.net.URL;

import javax.swing.JApplet;
import javax.swing.UIManager;

import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;

public class ArchaeopteryxA extends JApplet {

    private static final long  serialVersionUID      = 2314899014580484146L;
    private final static Color background_color      = new Color( 0, 0, 0 );
    private final static Color font_color            = new Color( 0, 255, 0 );
    private final static Color ex_background_color   = new Color( 0, 0, 0 );
    private final static Color ex_font_color         = new Color( 255, 0, 0 );
    private final static Font  font                  = new Font( Configuration.getDefaultFontFamilyName(), Font.BOLD, 9 );
    private MainFrameApplet    _mainframe_applet;
    private String             _tree_url_str         = "";
    private String             _species_tree_url_str = "";
    private String             _message_1            = "";
    private String             _message_2            = "";
    public final static String NAME                  = "ArchaeopteryxA";

    @Override
    public void destroy() {
        AptxUtil.printAppletMessage( NAME, "going to be destroyed" );
        if ( getMainFrameApplet() != null ) {
            getMainFrameApplet().close();
        }
    }

    /**
     * This method returns the current external node data which
     * has been selected by the user by clicking the "Return ..."
     * menu item. This method is expected to be called from Javascript or
     * something like it.
     * 
     * @return current external node data as String
     */
    public String getCurrentExternalNodesDataBuffer() {
        return getMainFrameApplet().getCurrentTreePanel().getCurrentExternalNodesDataBufferAsString();
    }

    public int getCurrentExternalNodesDataBufferChangeCounter() {
        return getMainFrameApplet().getCurrentTreePanel().getCurrentExternalNodesDataBufferChangeCounter();
    }

    public int getCurrentExternalNodesDataBufferLength() {
        return getMainFrameApplet().getCurrentTreePanel().getCurrentExternalNodesDataBufferAsString().length();
    }

    public String getTreeUrlStr() {
        return _tree_url_str;
    }

    public String getSpeciesTreeUrlStr() {
        return _species_tree_url_str;
    }

    @Override
    public void init() {
        boolean has_exception = false;
        setName( NAME );
        setTreeUrlStr( getParameter( Constants.APPLET_PARAM_NAME_FOR_URL_OF_TREE_TO_LOAD ) );
        setSpeciesTreeUrlStr( getParameter( Constants.APPLET_PARAM_NAME_FOR_URL_OF_SPECIES_TREE_TO_LOAD ) );
        if ( !ForesterUtil.isEmpty( getTreeUrlStr() ) ) {
            AptxUtil.printAppletMessage( NAME, "URL of tree(s) to load: \"" + getTreeUrlStr() + "\"" );
        }
        else {
            ForesterUtil.printErrorMessage( NAME, "no URL for tree(s) to load!" );
            setBackground( ex_background_color );
            setForeground( ex_font_color );
            has_exception = true;
            setMessage1( "no URL for tree(s) to load" );
            repaint();
        }
        if ( !ForesterUtil.isEmpty( getSpeciesTreeUrlStr() ) ) {
            AptxUtil.printAppletMessage( NAME, "URL of species tree to load: \"" + getSpeciesTreeUrlStr() + "\"" );
        }
        setBackground( background_color );
        setForeground( font_color );
        setFont( font );
        repaint();
        String s = null;
        try {
            s = System.getProperty( "java.version" );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( NAME, "minor error: " + e.getLocalizedMessage() );
        }
        if ( ( s != null ) && ( s.length() > 0 ) ) {
            setMessage2( "[Your Java version: " + s + "]" );
            repaint();
        }
        final String config_filename = getParameter( Constants.APPLET_PARAM_NAME_FOR_CONFIG_FILE_URL );
        AptxUtil.printAppletMessage( NAME, "URL for configuration file is: " + config_filename );
        final Configuration configuration = new Configuration( config_filename, true, true, true );
        try {
            if ( configuration.isUseNativeUI() ) {
                UIManager.setLookAndFeel( UIManager.getSystemLookAndFeelClassName() );
            }
            else {
                UIManager.setLookAndFeel( UIManager.getCrossPlatformLookAndFeelClassName() );
            }
            setVisible( false );
            _mainframe_applet = new MainFrameApplet( this, configuration );
            final URL tree_url = new URL( getTreeUrlStr() );
            final Phylogeny[] phys = AptxUtil.readPhylogeniesFromUrl( tree_url, configuration
                    .isValidatePhyloXmlAgainstSchema(), configuration.isReplaceUnderscoresInNhParsing(), configuration
                    .isInternalNumberAreConfidenceForNhParsing(), configuration.getTaxonomyExtraction() );
            AptxUtil.addPhylogeniesToTabs( phys,
                                           new File( tree_url.getFile() ).getName(),
                                           getTreeUrlStr(),
                                           getMainFrameApplet().getConfiguration(),
                                           getMainFrameApplet().getMainPanel() );
            if ( !ForesterUtil.isEmpty( getSpeciesTreeUrlStr() ) ) {
                final URL species_tree_url = new URL( getSpeciesTreeUrlStr() );
                final Phylogeny[] species_trees = AptxUtil
                        .readPhylogeniesFromUrl( species_tree_url,
                                                 configuration.isValidatePhyloXmlAgainstSchema(),
                                                 configuration.isReplaceUnderscoresInNhParsing(),
                                                 false,
                                                 TAXONOMY_EXTRACTION.NO );
                if ( ( species_trees != null ) && ( species_trees.length > 0 ) ) {
                    AptxUtil.printAppletMessage( NAME, "successfully read species tree" );
                    getMainFrameApplet().setSpeciesTree( species_trees[ 0 ] );
                }
            }
            getMainFrameApplet().getMainPanel().getControlPanel().showWholeAll();
            getMainFrameApplet().getMainPanel().getControlPanel().showWhole();
            setVisible( true );
        }
        catch ( final Exception e ) {
            ForesterUtil.printErrorMessage( NAME, e.toString() );
            setBackground( ex_background_color );
            setForeground( ex_font_color );
            has_exception = true;
            setMessage1( "Exception: " + e );
            e.printStackTrace();
            repaint();
        }
        if ( !has_exception ) {
            setMessage1( NAME + " is now ready!" );
            repaint();
            AptxUtil.printAppletMessage( NAME, "successfully initialized" );
        }
        KeyboardFocusManager.getCurrentKeyboardFocusManager().clearGlobalFocusOwner();
        getMainFrameApplet().requestFocus();
        getMainFrameApplet().requestFocusInWindow();
        getMainFrameApplet().requestFocus();
        /* GUILHEM_BEG */
        final String default_relation = getParameter( Constants.APPLET_PARAM_NAME_FOR_DEFAULT_SEQUENCE_RELATION_TYPE );
        if ( default_relation != null ) {
            getMainFrameApplet().getMainPanel().getControlPanel().getSequenceRelationTypeBox()
                    .setSelectedItem( default_relation );
        }
        final String default_sequence = getParameter( Constants.APPLET_PARAM_NAME_FOR_DEFAULT_QUERY_SEQUENCE );
        if ( default_sequence != null ) {
            getMainFrameApplet().getMainPanel().getControlPanel().getSequenceRelationBox()
                    .setSelectedItem( default_sequence );
            /* GUILHEM_END */
        }
    }

    /**
     * Prints message when initialization is finished. Called automatically.
     * 
     * @param g
     *            Graphics
     */
    @Override
    public void paint( final Graphics g ) {
        g.setColor( background_color );
        g.fillRect( 0, 0, 300, 60 );
        g.setColor( font_color );
        g.drawString( getMessage2(), 10, 20 );
        g.drawString( getMessage1(), 10, 40 );
    }

    @Override
    public void start() {
        getMainFrameApplet().getMainPanel().validate();
        getMainFrameApplet().requestFocus();
        getMainFrameApplet().requestFocusInWindow();
        getMainFrameApplet().requestFocus();
        AptxUtil.printAppletMessage( NAME, "started" );
    }

    private MainFrameApplet getMainFrameApplet() {
        return _mainframe_applet;
    }

    private String getMessage1() {
        return _message_1;
    }

    private String getMessage2() {
        return _message_2;
    }

    private void setMessage1( final String message_1 ) {
        _message_1 = message_1;
    }

    private void setMessage2( final String message_2 ) {
        _message_2 = message_2;
    }

    private void setTreeUrlStr( final String url_string ) {
        _tree_url_str = url_string;
    }

    private void setSpeciesTreeUrlStr( final String url_string ) {
        _species_tree_url_str = url_string;
    }
}
