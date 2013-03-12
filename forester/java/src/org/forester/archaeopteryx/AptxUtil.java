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
import java.awt.Component;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URI;
import java.net.URL;
import java.net.URLEncoder;
import java.text.ParseException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;
import javax.swing.JApplet;
import javax.swing.JOptionPane;
import javax.swing.text.MaskFormatter;

import org.forester.analysis.TaxonomyDataManager;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.io.parsers.tol.TolParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.AsciiHistogram;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceIdParser;
import org.forester.ws.seqdb.UniProtTaxonomy;

public final class AptxUtil {

    private final static Pattern  seq_identifier_pattern_1       = Pattern
                                                                         .compile( "^([A-Za-z]{2,5})[|=:]([0-9A-Za-z_\\.]{5,40})\\s*$" );
    private final static Pattern  seq_identifier_pattern_2       = Pattern
                                                                         .compile( "^([A-Za-z]{2,5})[|=:]([0-9A-Za-z_\\.]{5,40})[|,; ].*$" );
    private final static String[] AVAILABLE_FONT_FAMILIES_SORTED = GraphicsEnvironment.getLocalGraphicsEnvironment()
                                                                         .getAvailableFontFamilyNames();
    static {
        Arrays.sort( AVAILABLE_FONT_FAMILIES_SORTED );
    }

    public final static String createUriForSeqWeb( final PhylogenyNode node,
                                                   final Configuration conf,
                                                   final TreePanel tp ) {
        String uri_str = null;
        final String upkb = ForesterUtil.extractUniProtKbProteinSeqIdentifier( node );
        if ( !ForesterUtil.isEmpty( upkb ) ) {
            try {
                uri_str = ForesterUtil.UNIPROT_KB + URLEncoder.encode( upkb, ForesterConstants.UTF8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                showErrorMessage( tp, e.toString() );
                e.printStackTrace();
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = ForesterUtil.extractGenbankAccessor( node );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    if ( SequenceIdParser.isProtein( v ) ) {
                        uri_str = ForesterUtil.NCBI_PROTEIN + URLEncoder.encode( v, ForesterConstants.UTF8 );
                    }
                    else {
                        uri_str = ForesterUtil.NCBI_NUCCORE + URLEncoder.encode( v, ForesterConstants.UTF8 );
                    }
                }
                catch ( final UnsupportedEncodingException e ) {
                    showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = ForesterUtil.extractRefSeqAccessorAccessor( node );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    if ( SequenceIdParser.isProtein( v ) ) {
                        uri_str = ForesterUtil.NCBI_PROTEIN + URLEncoder.encode( v, ForesterConstants.UTF8 );
                    }
                    else {
                        uri_str = ForesterUtil.NCBI_NUCCORE + URLEncoder.encode( v, ForesterConstants.UTF8 );
                    }
                }
                catch ( final UnsupportedEncodingException e ) {
                    showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        if ( ForesterUtil.isEmpty( uri_str ) ) {
            final String v = ForesterUtil.extractGInumber( node );
            if ( !ForesterUtil.isEmpty( v ) ) {
                try {
                    uri_str = ForesterUtil.NCBI_GI + URLEncoder.encode( v, ForesterConstants.UTF8 );
                }
                catch ( final UnsupportedEncodingException e ) {
                    showErrorMessage( tp, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        return uri_str;
    }

    public static MaskFormatter createMaskFormatter( final String s ) {
        MaskFormatter formatter = null;
        try {
            formatter = new MaskFormatter( s );
        }
        catch ( final ParseException e ) {
            throw new IllegalArgumentException( e );
        }
        return formatter;
    }

    final static public boolean isHasAtLeastNodeWithEvent( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            if ( it.next().getNodeData().isHasEvent() ) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns true if at least one branch has a length larger than zero.
     * 
     * 
     * @param phy
     */
    final static public boolean isHasAtLeastOneBranchLengthLargerThanZero( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            if ( it.next().getDistanceToParent() > 0.0 ) {
                return true;
            }
        }
        return false;
    }

    final static public boolean isHasAtLeastOneBranchWithSupportValues( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            if ( it.next().getBranchData().isHasConfidences() ) {
                return true;
            }
        }
        return false;
    }

    final public static void launchWebBrowser( final URI uri,
                                               final boolean is_applet,
                                               final JApplet applet,
                                               final String frame_name ) throws IOException {
        if ( is_applet ) {
            applet.getAppletContext().showDocument( uri.toURL(), frame_name );
        }
        else {
            // This requires Java 1.6:
            // =======================
            // boolean no_desktop = false;
            // try {
            // if ( Desktop.isDesktopSupported() ) {
            // System.out.println( "desktop supported" );
            // final Desktop dt = Desktop.getDesktop();
            // dt.browse( uri );
            // }
            // else {
            // no_desktop = true;
            // }
            // }
            // catch ( final Exception ex ) {
            // ex.printStackTrace();
            // no_desktop = true;
            // }
            // catch ( final Error er ) {
            // er.printStackTrace();
            // no_desktop = true;
            // }
            // if ( no_desktop ) {
            // System.out.println( "desktop not supported" );
            try {
                openUrlInWebBrowser( uri.toString() );
            }
            catch ( final Exception e ) {
                throw new IOException( e );
            }
            // }
        }
    }

    public static Set<Taxonomy> obtainAllDistinctTaxonomies( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        final Set<Taxonomy> tax_set = new HashSet<Taxonomy>();
        for( final PhylogenyNode n : descs ) {
            if ( n.getNodeData().isHasTaxonomy() && !n.getNodeData().getTaxonomy().isEmpty() ) {
                tax_set.add( n.getNodeData().getTaxonomy() );
            }
        }
        return tax_set;
    }

    /**
     * Returns the set of distinct taxonomies of
     * all external nodes of node.
     * If at least one the external nodes has no taxonomy,
     * null is returned.
     * 
     */
    public static Set<Taxonomy> obtainDistinctTaxonomies( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        final Set<Taxonomy> tax_set = new HashSet<Taxonomy>();
        for( final PhylogenyNode n : descs ) {
            if ( !n.getNodeData().isHasTaxonomy() || n.getNodeData().getTaxonomy().isEmpty() ) {
                return null;
            }
            tax_set.add( n.getNodeData().getTaxonomy() );
        }
        return tax_set;
    }

    public final static Accession obtainSequenceAccessionFromName( final String sequence_name ) {
        final String n = sequence_name.trim();
        final Matcher matcher1 = seq_identifier_pattern_1.matcher( n );
        String group1 = "";
        String group2 = "";
        if ( matcher1.matches() ) {
            group1 = matcher1.group( 1 );
            group2 = matcher1.group( 2 );
        }
        else {
            final Matcher matcher2 = seq_identifier_pattern_2.matcher( n );
            if ( matcher2.matches() ) {
                group1 = matcher2.group( 1 );
                group2 = matcher2.group( 2 );
            }
        }
        if ( ForesterUtil.isEmpty( group1 ) || ForesterUtil.isEmpty( group2 ) ) {
            return null;
        }
        return new Accession( group2, group1 );
    }

    public final static void printWarningMessage( final String name, final String message ) {
        System.out.println( "[" + name + "] > " + message );
    }

    final public static void showErrorMessage( final Component parent, final String error_msg ) {
        printAppletMessage( Constants.PRG_NAME, error_msg );
        JOptionPane.showMessageDialog( parent, error_msg, "[" + Constants.PRG_NAME + " " + Constants.VERSION
                + "] Error", JOptionPane.ERROR_MESSAGE );
    }

    public final static void showExtDescNodeDataUserSelectedHelper( final ControlPanel cp,
                                                                    final PhylogenyNode node,
                                                                    final List<String> data ) {
        final StringBuilder sb = new StringBuilder();
        if ( cp.isShowNodeNames() && !ForesterUtil.isEmpty( node.getName() ) ) {
            showExtDescNodeDataUserSelectedHelperHelper( node.getName(), sb );
        }
        if ( cp.isShowGeneNames() && node.getNodeData().isHasSequence()
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getName() ) ) {
            showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence().getName(), sb );
        }
        if ( cp.isShowGeneSymbols() && node.getNodeData().isHasSequence()
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getSymbol() ) ) {
            showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence().getSymbol(), sb );
        }
        if ( cp.isShowSequenceAcc() && node.getNodeData().isHasSequence()
                && ( node.getNodeData().getSequence().getAccession() != null )
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().toString() ) ) {
            showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence().getAccession().toString(), sb );
        }
        if ( cp.isShowTaxonomyCode() && node.getNodeData().isHasTaxonomy()
                && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
            showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getTaxonomy().getTaxonomyCode(), sb );
        }
        if ( cp.isShowTaxonomyScientificNames() && node.getNodeData().isHasTaxonomy()
                && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) ) {
            showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getTaxonomy().getScientificName(), sb );
        }
        if ( ( cp.isShowGeneNames() || cp.isShowGeneSymbols() || cp.isShowSequenceAcc() )
                && node.getNodeData().isHasSequence()
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() ) ) {
            showExtDescNodeDataUserSelectedHelperHelper( node.getNodeData().getSequence().getMolecularSequence(), sb );
        }
        final String s = sb.toString().trim();
        if ( !ForesterUtil.isEmpty( s ) ) {
            data.add( s );
        }
    }

    public final static void showExtDescNodeDataUserSelectedHelperHelper( final String s, final StringBuilder sb ) {
        if ( sb.length() > 0 ) {
            sb.append( "\t" );
        }
        sb.append( s );
    }

    final public static void showInformationMessage( final Component parent, final String title, final String msg ) {
        JOptionPane.showMessageDialog( parent, msg, title, JOptionPane.INFORMATION_MESSAGE );
    }

    public static void writePhylogenyToGraphicsFile( final File intree,
                                                     final File outfile,
                                                     final int width,
                                                     final int height,
                                                     final GraphicsExportType type,
                                                     final Configuration config ) throws IOException {
        final PhylogenyParser parser = ParserUtils.createParserDependingOnFileType( intree, true );
        Phylogeny[] phys = null;
        phys = PhylogenyMethods.readPhylogenies( parser, intree );
        writePhylogenyToGraphicsFile( phys[ 0 ], outfile, width, height, type, config );
    }

    public static void writePhylogenyToGraphicsFile( final Phylogeny phy,
                                                     final File outfile,
                                                     final int width,
                                                     final int height,
                                                     final GraphicsExportType type,
                                                     final Configuration config ) throws IOException {
        final Phylogeny[] phys = new Phylogeny[ 1 ];
        phys[ 0 ] = phy;
        final MainFrameApplication mf = MainFrameApplication.createInstance( phys, config );
        AptxUtil.writePhylogenyToGraphicsFileNonInteractive( outfile, width, height, mf.getMainPanel()
                .getCurrentTreePanel(), mf.getMainPanel().getControlPanel(), type, mf.getOptions() );
        mf.end();
    }

    public final static void writePhylogenyToGraphicsFileNonInteractive( final File outfile,
                                                                         final int width,
                                                                         final int height,
                                                                         final TreePanel tree_panel,
                                                                         final ControlPanel ac,
                                                                         final GraphicsExportType type,
                                                                         final Options options ) throws IOException {
        tree_panel.calcParametersForPainting( width, height, true );
        tree_panel.resetPreferredSize();
        tree_panel.repaint();
        final RenderingHints rendering_hints = new RenderingHints( RenderingHints.KEY_RENDERING,
                                                                   RenderingHints.VALUE_RENDER_QUALITY );
        rendering_hints.put( RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY );
        if ( options.isAntialiasPrint() ) {
            rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        }
        else {
            rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF );
            rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF );
        }
        final Phylogeny phylogeny = tree_panel.getPhylogeny();
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return;
        }
        if ( outfile.isDirectory() ) {
            throw new IOException( "\"" + outfile + "\" is a directory" );
        }
        final BufferedImage buffered_img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g2d = buffered_img.createGraphics();
        g2d.setRenderingHints( rendering_hints );
        tree_panel.paintPhylogeny( g2d, false, true, width, height, 0, 0 );
        if ( type == GraphicsExportType.TIFF ) {
            writeToTiff( outfile, buffered_img );
        }
        else {
            ImageIO.write( buffered_img, type.toString(), outfile );
        }
        g2d.dispose();
    }

    final static void addPhylogeniesToTabs( final Phylogeny[] phys,
                                            final String default_name,
                                            final String full_path,
                                            final Configuration configuration,
                                            final MainPanel main_panel ) {
        if ( phys.length > Constants.MAX_TREES_TO_LOAD ) {
            JOptionPane.showMessageDialog( main_panel, "Attempt to load " + phys.length
                    + " phylogenies,\ngoing to load only the first " + Constants.MAX_TREES_TO_LOAD, Constants.PRG_NAME
                    + " more than " + Constants.MAX_TREES_TO_LOAD + " phylogenies", JOptionPane.WARNING_MESSAGE );
        }
        int i = 1;
        for( final Phylogeny phy : phys ) {
            if ( !phy.isEmpty() ) {
                if ( i <= Constants.MAX_TREES_TO_LOAD ) {
                    String my_name = "";
                    String my_name_for_file = "";
                    if ( phys.length > 1 ) {
                        if ( !ForesterUtil.isEmpty( default_name ) ) {
                            my_name = new String( default_name );
                        }
                        if ( !ForesterUtil.isEmpty( full_path ) ) {
                            my_name_for_file = new String( full_path );
                        }
                        else if ( !ForesterUtil.isEmpty( default_name ) ) {
                            my_name_for_file = new String( default_name );
                        }
                        String suffix = "";
                        if ( my_name_for_file.indexOf( '.' ) > 0 ) {
                            suffix = my_name_for_file.substring( my_name_for_file.lastIndexOf( '.' ),
                                                                 my_name_for_file.length() );
                            my_name_for_file = my_name_for_file.substring( 0, my_name_for_file.lastIndexOf( '.' ) );
                        }
                        if ( !ForesterUtil.isEmpty( my_name_for_file ) ) {
                            my_name_for_file += "_";
                        }
                        if ( !ForesterUtil.isEmpty( phy.getName() ) ) {
                            my_name_for_file += phy.getName().replaceAll( " ", "_" );
                        }
                        else if ( phy.getIdentifier() != null ) {
                            final StringBuffer sb = new StringBuffer();
                            if ( !ForesterUtil.isEmpty( phy.getIdentifier().getProvider() ) ) {
                                sb.append( phy.getIdentifier().getProvider() );
                                sb.append( "_" );
                            }
                            sb.append( phy.getIdentifier().getValue() );
                            my_name_for_file += sb;
                        }
                        else {
                            my_name_for_file += i;
                        }
                        if ( !ForesterUtil.isEmpty( my_name ) && ForesterUtil.isEmpty( phy.getName() )
                                && ( phy.getIdentifier() == null ) ) {
                            my_name = my_name + " [" + i + "]";
                        }
                        if ( !ForesterUtil.isEmpty( suffix ) ) {
                            my_name_for_file += suffix;
                        }
                    }
                    else {
                        if ( !ForesterUtil.isEmpty( default_name ) ) {
                            my_name = new String( default_name );
                        }
                        my_name_for_file = "";
                        if ( !ForesterUtil.isEmpty( full_path ) ) {
                            my_name_for_file = new String( full_path );
                        }
                        else if ( !ForesterUtil.isEmpty( default_name ) ) {
                            my_name_for_file = new String( default_name );
                        }
                        if ( ForesterUtil.isEmpty( my_name_for_file ) ) {
                            if ( !ForesterUtil.isEmpty( phy.getName() ) ) {
                                my_name_for_file = new String( phy.getName() ).replaceAll( " ", "_" );
                            }
                            else if ( phy.getIdentifier() != null ) {
                                final StringBuffer sb = new StringBuffer();
                                if ( !ForesterUtil.isEmpty( phy.getIdentifier().getProvider() ) ) {
                                    sb.append( phy.getIdentifier().getProvider() );
                                    sb.append( "_" );
                                }
                                sb.append( phy.getIdentifier().getValue() );
                                my_name_for_file = new String( sb.toString().replaceAll( " ", "_" ) );
                            }
                        }
                    }
                    main_panel.addPhylogenyInNewTab( phy, configuration, my_name, full_path );
                    main_panel.getCurrentTreePanel().setTreeFile( new File( my_name_for_file ) );
                    lookAtSomeTreePropertiesForAptxControlSettings( phy, main_panel.getControlPanel(), configuration );
                    ++i;
                }
            }
        }
    }

    final static void addPhylogenyToPanel( final Phylogeny[] phys,
                                           final Configuration configuration,
                                           final MainPanel main_panel ) {
        final Phylogeny phy = phys[ 0 ];
        main_panel.addPhylogenyInPanel( phy, configuration );
        lookAtSomeTreePropertiesForAptxControlSettings( phy, main_panel.getControlPanel(), configuration );
    }

    final static Color calculateColorFromString( final String str ) {
        final String species_uc = str.toUpperCase();
        char first = species_uc.charAt( 0 );
        char second = ' ';
        char third = ' ';
        if ( species_uc.length() > 1 ) {
            second = species_uc.charAt( 1 );
            if ( species_uc.length() > 2 ) {
                if ( species_uc.indexOf( " " ) > 0 ) {
                    third = species_uc.charAt( species_uc.indexOf( " " ) + 1 );
                }
                else {
                    third = species_uc.charAt( 2 );
                }
            }
        }
        first = AptxUtil.normalizeCharForRGB( first );
        second = AptxUtil.normalizeCharForRGB( second );
        third = AptxUtil.normalizeCharForRGB( third );
        if ( ( first > 235 ) && ( second > 235 ) && ( third > 235 ) ) {
            first = 0;
        }
        else if ( ( first < 80 ) && ( second < 80 ) && ( third < 80 ) ) {
            second = 255;
        }
        return new Color( first, second, third );
    }

    // Returns true if the specified format name can be written
    final static boolean canWriteFormat( final String format_name ) {
        final Iterator<ImageWriter> iter = ImageIO.getImageWritersByFormatName( format_name );
        return iter.hasNext();
    }

    final static void collapseSpeciesSpecificSubtrees( final Phylogeny phy ) {
        boolean inferred = false;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.isCollapse() && ( n.getNumberOfDescendants() > 1 ) ) {
                final Set<Taxonomy> taxs = obtainDistinctTaxonomies( n );
                if ( ( taxs != null ) && ( taxs.size() == 1 ) ) {
                    AptxUtil.collapseSubtree( n, true );
                    if ( !n.getNodeData().isHasTaxonomy() ) {
                        n.getNodeData().setTaxonomy( ( Taxonomy ) n.getAllExternalDescendants().get( 0 ).getNodeData()
                                .getTaxonomy().copy() );
                    }
                    inferred = true;
                }
                else {
                    n.setCollapse( false );
                }
            }
        }
        if ( inferred ) {
            phy.setRerootable( false );
        }
    }

    final static void collapseSubtree( final PhylogenyNode node, final boolean collapse ) {
        node.setCollapse( collapse );
        if ( node.isExternal() ) {
            return;
        }
        final PhylogenyNodeIterator it = new PreorderTreeIterator( node );
        while ( it.hasNext() ) {
            it.next().setCollapse( collapse );
        }
    }

    final static void colorPhylogenyAccordingToConfidenceValues( final Phylogeny tree, final TreePanel tree_panel ) {
        double max_conf = 0.0;
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.getBranchData().setBranchColor( null );
            if ( n.getBranchData().isHasConfidences() ) {
                final double conf = PhylogenyMethods.getConfidenceValue( n );
                if ( conf > max_conf ) {
                    max_conf = conf;
                }
            }
        }
        if ( max_conf > 0.0 ) {
            final Color bg = tree_panel.getTreeColorSet().getBackgroundColor();
            final Color br = tree_panel.getTreeColorSet().getBranchColor();
            for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( n.getBranchData().isHasConfidences() ) {
                    final double conf = PhylogenyMethods.getConfidenceValue( n );
                    final BranchColor c = new BranchColor( ForesterUtil.calcColor( conf, 0.0, max_conf, bg, br ) );
                    colorizeSubtree( n, c );
                }
            }
        }
    }

    final static void colorPhylogenyAccordingToExternalTaxonomy( final Phylogeny tree, final TreePanel tree_panel ) {
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            it.next().getBranchData().setBranchColor( null );
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.getBranchData().isHasBranchColor() ) {
                final Taxonomy tax = PhylogenyMethods.getExternalDescendantsTaxonomy( n );
                if ( tax != null ) {
                    n.getBranchData().setBranchColor( new BranchColor( tree_panel.calculateTaxonomyBasedColor( tax ) ) );
                    final List<PhylogenyNode> descs = PhylogenyMethods.getAllDescendants( n );
                    for( final PhylogenyNode desc : descs ) {
                        desc.getBranchData()
                                .setBranchColor( new BranchColor( tree_panel.calculateTaxonomyBasedColor( tax ) ) );
                    }
                }
            }
        }
    }

    final static int colorPhylogenyAccordingToRanks( final Phylogeny tree, final String rank, final TreePanel tree_panel ) {
        final Map<String, Color> true_lineage_to_color_map = new HashMap<String, Color>();
        int colorizations = 0;
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy()
                    && ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() )
                            || !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getCommonName() ) || !ForesterUtil
                            .isEmpty( n.getNodeData().getTaxonomy().getTaxonomyCode() ) ) ) {
                if ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getRank() )
                        && n.getNodeData().getTaxonomy().getRank().equalsIgnoreCase( rank ) ) {
                    final BranchColor c = new BranchColor( tree_panel.calculateTaxonomyBasedColor( n.getNodeData()
                            .getTaxonomy() ) );
                    colorizeSubtree( n, c );
                    ++colorizations;
                    if ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                        true_lineage_to_color_map.put( n.getNodeData().getTaxonomy().getScientificName(), c.getValue() );
                    }
                }
            }
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( ( node.getBranchData().getBranchColor() == null ) && node.getNodeData().isHasTaxonomy()
                    && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getLineage() ) ) {
                boolean success = false;
                if ( !true_lineage_to_color_map.isEmpty() ) {
                    for( final String lin : node.getNodeData().getTaxonomy().getLineage() ) {
                        if ( true_lineage_to_color_map.containsKey( lin ) ) {
                            colorizeSubtree( node, new BranchColor( true_lineage_to_color_map.get( lin ) ) );
                            ++colorizations;
                            success = true;
                            break;
                        }
                    }
                }
                if ( !success ) {
                    final Map<String, String> lineage_to_rank_map = MainPanel.getLineageToRankMap();
                    for( final String lin : node.getNodeData().getTaxonomy().getLineage() ) {
                        final Taxonomy temp_tax = new Taxonomy();
                        temp_tax.setScientificName( lin );
                        if ( lineage_to_rank_map.containsKey( lin )
                                && !ForesterUtil.isEmpty( lineage_to_rank_map.get( lin ) )
                                && lineage_to_rank_map.get( lin ).equalsIgnoreCase( rank ) ) {
                            final BranchColor c = new BranchColor( tree_panel.calculateTaxonomyBasedColor( temp_tax ) );
                            colorizeSubtree( node, c );
                            ++colorizations;
                            true_lineage_to_color_map.put( lin, c.getValue() );
                            break;
                        }
                        else {
                            UniProtTaxonomy up = null;
                            try {
                                up = TaxonomyDataManager.obtainUniProtTaxonomy( temp_tax, null, null );
                            }
                            catch ( final Exception e ) {
                                e.printStackTrace();
                            }
                            if ( ( up != null ) && !ForesterUtil.isEmpty( up.getRank() ) ) {
                                lineage_to_rank_map.put( lin, up.getRank() );
                                if ( up.getRank().equalsIgnoreCase( rank ) ) {
                                    final BranchColor c = new BranchColor( tree_panel.calculateTaxonomyBasedColor( temp_tax ) );
                                    colorizeSubtree( node, c );
                                    ++colorizations;
                                    true_lineage_to_color_map.put( lin, c.getValue() );
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        return colorizations;
    }

    final static String createBasicInformation( final Phylogeny phy ) {
        final StringBuilder desc = new StringBuilder();
        if ( ( phy != null ) && !phy.isEmpty() ) {
            if ( !ForesterUtil.isEmpty( phy.getName() ) ) {
                desc.append( "Name: " );
                desc.append( phy.getName() );
                desc.append( "\n" );
            }
            if ( phy.getIdentifier() != null ) {
                desc.append( "Id: " );
                desc.append( phy.getIdentifier().toString() );
                desc.append( "\n" );
            }
            if ( !ForesterUtil.isEmpty( phy.getDescription() ) ) {
                desc.append( "Description: " );
                desc.append( phy.getDescription() );
                desc.append( "\n" );
            }
            if ( !ForesterUtil.isEmpty( phy.getDistanceUnit() ) ) {
                desc.append( "Distance Unit: " );
                desc.append( phy.getDistanceUnit() );
                desc.append( "\n" );
            }
            if ( !ForesterUtil.isEmpty( phy.getType() ) ) {
                desc.append( "Type: " );
                desc.append( phy.getType() );
                desc.append( "\n" );
            }
            desc.append( "Rooted: " );
            desc.append( phy.isRooted() );
            desc.append( "\n" );
            desc.append( "Rerootable: " );
            desc.append( phy.isRerootable() );
            desc.append( "\n" );
            desc.append( "Nodes: " );
            desc.append( phy.getNodeCount() );
            desc.append( "\n" );
            desc.append( "External nodes: " );
            desc.append( phy.getNumberOfExternalNodes() );
            desc.append( "\n" );
            desc.append( "Internal nodes: " );
            desc.append( phy.getNodeCount() - phy.getNumberOfExternalNodes() );
            desc.append( "\n" );
            desc.append( "Internal nodes with polytomies: " );
            desc.append( PhylogenyMethods.countNumberOfPolytomies( phy ) );
            desc.append( "\n" );
            desc.append( "Branches: " );
            desc.append( phy.getNumberOfBranches() );
            desc.append( "\n" );
            desc.append( "Depth: " );
            desc.append( PhylogenyMethods.calculateMaxDepth( phy ) );
            desc.append( "\n" );
            desc.append( "Maximum distance to root: " );
            desc.append( ForesterUtil.round( PhylogenyMethods.calculateMaxDistanceToRoot( phy ), 6 ) );
            desc.append( "\n" );
            final Set<Taxonomy> taxs = obtainAllDistinctTaxonomies( phy.getRoot() );
            if ( taxs != null ) {
                desc.append( "Distinct external taxonomies: " );
                desc.append( taxs.size() );
            }
            desc.append( "\n" );
            final DescriptiveStatistics bs = PhylogenyMethods.calculatBranchLengthStatistics( phy );
            if ( bs.getN() > 3 ) {
                desc.append( "\n" );
                desc.append( "Branch-length statistics: " );
                desc.append( "\n" );
                desc.append( "    Number of branches with non-negative branch-lengths: " + bs.getN() );
                desc.append( "\n" );
                desc.append( "    Median: " + ForesterUtil.round( bs.median(), 6 ) );
                desc.append( "\n" );
                desc.append( "    Mean: " + ForesterUtil.round( bs.arithmeticMean(), 6 ) );
                desc.append( "\n" );
                desc.append( "    SD: " + ForesterUtil.round( bs.sampleStandardDeviation(), 6 ) );
                desc.append( "\n" );
                desc.append( "    Minimum: " + ForesterUtil.round( bs.getMin(), 6 ) );
                desc.append( "\n" );
                desc.append( "    Maximum: " + ForesterUtil.round( bs.getMax(), 6 ) );
                desc.append( "\n" );
                if ( Math.abs( bs.getMax() - bs.getMin() ) > 0.0001 ) {
                    desc.append( "\n" );
                    final AsciiHistogram histo = new AsciiHistogram( bs );
                    desc.append( histo.toStringBuffer( 12, '#', 40, 7, "    " ) );
                }
            }
            final DescriptiveStatistics ds = PhylogenyMethods.calculatNumberOfDescendantsPerNodeStatistics( phy );
            if ( ds.getN() > 2 ) {
                desc.append( "\n" );
                desc.append( "Descendants per node statistics: " );
                desc.append( "\n" );
                desc.append( "    Median: " + ForesterUtil.round( ds.median(), 2 ) );
                desc.append( "\n" );
                desc.append( "    Mean: " + ForesterUtil.round( ds.arithmeticMean(), 2 ) );
                desc.append( "\n" );
                desc.append( "    SD: " + ForesterUtil.round( ds.sampleStandardDeviation(), 2 ) );
                desc.append( "\n" );
                desc.append( "    Minimum: " + ForesterUtil.roundToInt( ds.getMin() ) );
                desc.append( "\n" );
                desc.append( "    Maximum: " + ForesterUtil.roundToInt( ds.getMax() ) );
                desc.append( "\n" );
            }
            List<DescriptiveStatistics> css = null;
            try {
                css = PhylogenyMethods.calculatConfidenceStatistics( phy );
            }
            catch ( final IllegalArgumentException e ) {
                ForesterUtil.printWarningMessage( Constants.PRG_NAME, e.getMessage() );
            }
            if ( ( css != null ) && ( css.size() > 0 ) ) {
                desc.append( "\n" );
                for( int i = 0; i < css.size(); ++i ) {
                    final DescriptiveStatistics cs = css.get( i );
                    if ( ( cs != null ) && ( cs.getN() > 1 ) ) {
                        if ( css.size() > 1 ) {
                            desc.append( "Support statistics " + ( i + 1 ) + ": " );
                        }
                        else {
                            desc.append( "Support statistics: " );
                        }
                        if ( !ForesterUtil.isEmpty( cs.getDescription() ) ) {
                            desc.append( "\n" );
                            desc.append( "    Type: " + cs.getDescription() );
                        }
                        desc.append( "\n" );
                        desc.append( "    Branches with support: " + cs.getN() );
                        desc.append( "\n" );
                        desc.append( "    Median: " + ForesterUtil.round( cs.median(), 6 ) );
                        desc.append( "\n" );
                        desc.append( "    Mean: " + ForesterUtil.round( cs.arithmeticMean(), 6 ) );
                        desc.append( "\n" );
                        if ( cs.getN() > 2 ) {
                            desc.append( "    SD: " + ForesterUtil.round( cs.sampleStandardDeviation(), 6 ) );
                            desc.append( "\n" );
                        }
                        desc.append( "    Minimum: " + ForesterUtil.roundToInt( cs.getMin() ) );
                        desc.append( "\n" );
                        desc.append( "    Maximum: " + ForesterUtil.roundToInt( cs.getMax() ) );
                        desc.append( "\n" );
                    }
                }
            }
        }
        return desc.toString();
    }

    /**
     * Exits with -1.
     * 
     * 
     * @param message
     *            to message to be printed
     */
    final static void dieWithSystemError( final String message ) {
        System.out.println();
        System.out.println( Constants.PRG_NAME + " encountered the following system error: " + message );
        System.out.println( "Please contact the authors." );
        System.out.println( Constants.PRG_NAME + " needs to close." );
        System.out.println();
        System.exit( -1 );
    }

    final static String[] getAllPossibleRanks() {
        final String[] str_array = new String[ PhyloXmlUtil.TAXONOMY_RANKS_LIST.size() - 2 ];
        int i = 0;
        for( final String e : PhyloXmlUtil.TAXONOMY_RANKS_LIST ) {
            if ( !e.equals( PhyloXmlUtil.UNKNOWN ) && !e.equals( PhyloXmlUtil.OTHER ) ) {
                str_array[ i++ ] = e;
            }
        }
        return str_array;
    }

    final static String[] getAllRanks( final Phylogeny tree ) {
        final SortedSet<String> ranks = new TreeSet<String>();
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() && !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getRank() ) ) {
                ranks.add( n.getNodeData().getTaxonomy().getRank() );
            }
        }
        return ForesterUtil.stringSetToArray( ranks );
    }

    final static String[] getAvailableFontFamiliesSorted() {
        return AVAILABLE_FONT_FAMILIES_SORTED;
    }

    final static boolean isHasAssignedEvent( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasEvent() ) {
            return false;
        }
        if ( ( node.getNodeData().getEvent() ).isUnassigned() ) {
            return false;
        }
        return true;
    }

    final static boolean isMac() {
        try {
            final String s = ForesterUtil.OS_NAME.toLowerCase();
            return s.startsWith( "mac" );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "minor error: " + e );
            return false;
        }
    }

    final static boolean isUsOrCanada() {
        try {
            if ( ( Locale.getDefault().equals( Locale.CANADA ) ) || ( Locale.getDefault().equals( Locale.US ) ) ) {
                return true;
            }
        }
        catch ( final Exception e ) {
            return false;
        }
        return false;
    }

    final static boolean isWindows() {
        try {
            final String s = ForesterUtil.OS_NAME.toLowerCase();
            return s.indexOf( "win" ) > -1;
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "minor error: " + e );
            return false;
        }
    }

    final static void lookAtSomeTreePropertiesForAptxControlSettings( final Phylogeny t,
                                                                      final ControlPanel atv_control,
                                                                      final Configuration configuration ) {
        if ( ( t != null ) && !t.isEmpty() ) {
            if ( !AptxUtil.isHasAtLeastOneBranchLengthLargerThanZero( t ) ) {
                atv_control.setDrawPhylogram( false );
                atv_control.setDrawPhylogramEnabled( false );
            }
            if ( configuration.doGuessCheckOption( Configuration.display_as_phylogram ) ) {
                if ( atv_control.getDisplayAsPhylogramCb() != null ) {
                    if ( AptxUtil.isHasAtLeastOneBranchLengthLargerThanZero( t ) ) {
                        atv_control.setDrawPhylogram( true );
                        atv_control.setDrawPhylogramEnabled( true );
                    }
                    else {
                        atv_control.setDrawPhylogram( false );
                    }
                }
            }
            if ( configuration.doGuessCheckOption( Configuration.write_confidence_values ) ) {
                if ( atv_control.getWriteConfidenceCb() != null ) {
                    if ( AptxUtil.isHasAtLeastOneBranchWithSupportValues( t ) ) {
                        atv_control.setCheckbox( Configuration.write_confidence_values, true );
                    }
                    else {
                        atv_control.setCheckbox( Configuration.write_confidence_values, false );
                    }
                }
            }
            if ( configuration.doGuessCheckOption( Configuration.write_events ) ) {
                if ( atv_control.getShowEventsCb() != null ) {
                    if ( AptxUtil.isHasAtLeastNodeWithEvent( t ) ) {
                        atv_control.setCheckbox( Configuration.write_events, true );
                    }
                    else {
                        atv_control.setCheckbox( Configuration.write_events, false );
                    }
                }
            }
        }
    }

    final static void openWebsite( final String url, final boolean is_applet, final JApplet applet ) throws IOException {
        try {
            AptxUtil.launchWebBrowser( new URI( url ), is_applet, applet, Constants.PRG_NAME );
        }
        catch ( final Exception e ) {
            throw new IOException( e );
        }
    }

    final static void outOfMemoryError( final OutOfMemoryError e ) {
        System.err.println();
        System.err.println( "Java memory allocation might be too small, try \"-Xmx2048m\" java command line option" );
        System.err.println();
        e.printStackTrace();
        System.err.println();
        JOptionPane.showMessageDialog( null,
                                       "Java memory allocation might be too small, try \"-Xmx2048m\" java command line option"
                                               + "\n\nError: " + e.getLocalizedMessage(),
                                       "Out of Memory Error [" + Constants.PRG_NAME + " " + Constants.VERSION + "]",
                                       JOptionPane.ERROR_MESSAGE );
        System.exit( -1 );
    }

    final static void printAppletMessage( final String applet_name, final String message ) {
        System.out.println( "[" + applet_name + "] > " + message );
    }

    final static Phylogeny[] readPhylogeniesFromUrl( final URL url,
                                                     final boolean phyloxml_validate_against_xsd,
                                                     final boolean replace_underscores,
                                                     final boolean internal_numbers_are_confidences,
                                                     final TAXONOMY_EXTRACTION taxonomy_extraction,
                                                     final boolean midpoint_reroot ) throws FileNotFoundException,
            IOException {
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final PhylogenyParser parser;
        boolean nhx_or_nexus = false;
        if ( url.getHost().toLowerCase().indexOf( "tolweb" ) >= 0 ) {
            parser = new TolParser();
        }
        else {
            parser = ParserUtils.createParserDependingOnUrlContents( url, phyloxml_validate_against_xsd );
            if ( parser instanceof NHXParser ) {
                nhx_or_nexus = true;
                final NHXParser nhx = ( NHXParser ) parser;
                nhx.setReplaceUnderscores( replace_underscores );
                nhx.setIgnoreQuotes( false );
                nhx.setTaxonomyExtraction( taxonomy_extraction );
            }
            else if ( parser instanceof NexusPhylogeniesParser ) {
                nhx_or_nexus = true;
                final NexusPhylogeniesParser nex = ( NexusPhylogeniesParser ) parser;
                nex.setReplaceUnderscores( replace_underscores );
                nex.setIgnoreQuotes( false );
            }
        }
        final Phylogeny[] phys = factory.create( url.openStream(), parser );
        if ( nhx_or_nexus && internal_numbers_are_confidences ) {
            for( final Phylogeny phy : phys ) {
                PhylogenyMethods.transferInternalNodeNamesToConfidence( phy );
            }
        }
        if ( midpoint_reroot ) {
            for( final Phylogeny phy : phys ) {
                PhylogenyMethods.midpointRoot( phy );
                PhylogenyMethods.orderAppearance( phy.getRoot(), true, true, DESCENDANT_SORT_PRIORITY.NODE_NAME );
            }
        }
        return phys;
    }

    final static void removeBranchColors( final Phylogeny phy ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            it.next().getBranchData().setBranchColor( null );
        }
    }

    final static void unexpectedError( final Error e ) {
        System.err.println();
        e.printStackTrace( System.err );
        System.err.println();
        final StringBuffer sb = new StringBuffer();
        for( final StackTraceElement s : e.getStackTrace() ) {
            sb.append( s + "\n" );
        }
        JOptionPane
                .showMessageDialog( null,
                                    "An unexpected (possibly severe) error has occured - terminating. \nPlease contact: "
                                            + Constants.AUTHOR_EMAIL + " \nError: " + e.getLocalizedMessage() + "\n"
                                            + sb,
                                    "Unexpected Severe Error [" + Constants.PRG_NAME + " " + Constants.VERSION + "]",
                                    JOptionPane.ERROR_MESSAGE );
        System.exit( -1 );
    }

    final static void unexpectedException( final Exception e ) {
        System.err.println();
        e.printStackTrace( System.err );
        System.err.println();
        final StringBuffer sb = new StringBuffer();
        for( final StackTraceElement s : e.getStackTrace() ) {
            sb.append( s + "\n" );
        }
        JOptionPane.showMessageDialog( null,
                                       "An unexpected exception has occured. \nPlease contact: "
                                               + Constants.AUTHOR_EMAIL + " \nException: " + e.getLocalizedMessage()
                                               + "\n" + sb,
                                       "Unexpected Exception [" + Constants.PRG_NAME + Constants.VERSION + "]",
                                       JOptionPane.ERROR_MESSAGE );
    }

    final static String writePhylogenyToGraphicsByteArrayOutputStream( final ByteArrayOutputStream baos,
                                                                       int width,
                                                                       int height,
                                                                       final TreePanel tree_panel,
                                                                       final ControlPanel ac,
                                                                       final GraphicsExportType type,
                                                                       final Options options ) throws IOException {
        if ( !options.isGraphicsExportUsingActualSize() ) {
            if ( options.isGraphicsExportVisibleOnly() ) {
                throw new IllegalArgumentException( "cannot export visible rectangle only without exporting in actual size" );
            }
            tree_panel.calcParametersForPainting( options.getPrintSizeX(), options.getPrintSizeY(), true );
            tree_panel.resetPreferredSize();
            tree_panel.repaint();
        }
        final RenderingHints rendering_hints = new RenderingHints( RenderingHints.KEY_RENDERING,
                                                                   RenderingHints.VALUE_RENDER_QUALITY );
        rendering_hints.put( RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY );
        if ( options.isAntialiasPrint() ) {
            rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        }
        else {
            rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF );
            rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF );
        }
        final Phylogeny phylogeny = tree_panel.getPhylogeny();
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return "";
        }
        Rectangle visible = null;
        if ( !options.isGraphicsExportUsingActualSize() ) {
            width = options.getPrintSizeX();
            height = options.getPrintSizeY();
        }
        else if ( options.isGraphicsExportVisibleOnly() ) {
            visible = tree_panel.getVisibleRect();
            width = visible.width;
            height = visible.height;
        }
        final BufferedImage buffered_img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB );
        Graphics2D g2d = buffered_img.createGraphics();
        g2d.setRenderingHints( rendering_hints );
        int x = 0;
        int y = 0;
        if ( options.isGraphicsExportVisibleOnly() ) {
            g2d = ( Graphics2D ) g2d.create( -visible.x, -visible.y, visible.width, visible.height );
            g2d.setClip( null );
            x = visible.x;
            y = visible.y;
        }
        tree_panel.paintPhylogeny( g2d, false, true, width, height, x, y );
        ImageIO.write( buffered_img, type.toString(), baos );
        g2d.dispose();
        System.gc();
        if ( !options.isGraphicsExportUsingActualSize() ) {
            tree_panel.getMainPanel().getControlPanel().showWhole();
        }
        String msg = baos.toString();
        if ( ( width > 0 ) && ( height > 0 ) ) {
            msg += " [size: " + width + ", " + height + "]";
        }
        return msg;
    }

    final static String writePhylogenyToGraphicsFile( final String file_name,
                                                      int width,
                                                      int height,
                                                      final TreePanel tree_panel,
                                                      final ControlPanel ac,
                                                      final GraphicsExportType type,
                                                      final Options options ) throws IOException {
        if ( !options.isGraphicsExportUsingActualSize() ) {
            if ( options.isGraphicsExportVisibleOnly() ) {
                throw new IllegalArgumentException( "cannot export visible rectangle only without exporting in actual size" );
            }
            tree_panel.calcParametersForPainting( options.getPrintSizeX(), options.getPrintSizeY(), true );
            tree_panel.resetPreferredSize();
            tree_panel.repaint();
        }
        final RenderingHints rendering_hints = new RenderingHints( RenderingHints.KEY_RENDERING,
                                                                   RenderingHints.VALUE_RENDER_QUALITY );
        rendering_hints.put( RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY );
        if ( options.isAntialiasPrint() ) {
            rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        }
        else {
            rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF );
            rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF );
        }
        final Phylogeny phylogeny = tree_panel.getPhylogeny();
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return "";
        }
        final File file = new File( file_name );
        if ( file.isDirectory() ) {
            throw new IOException( "\"" + file_name + "\" is a directory" );
        }
        Rectangle visible = null;
        if ( !options.isGraphicsExportUsingActualSize() ) {
            width = options.getPrintSizeX();
            height = options.getPrintSizeY();
        }
        else if ( options.isGraphicsExportVisibleOnly() ) {
            visible = tree_panel.getVisibleRect();
            width = visible.width;
            height = visible.height;
        }
        final BufferedImage buffered_img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB );
        Graphics2D g2d = buffered_img.createGraphics();
        g2d.setRenderingHints( rendering_hints );
        int x = 0;
        int y = 0;
        if ( options.isGraphicsExportVisibleOnly() ) {
            g2d = ( Graphics2D ) g2d.create( -visible.x, -visible.y, visible.width, visible.height );
            g2d.setClip( null );
            x = visible.x;
            y = visible.y;
        }
        tree_panel.paintPhylogeny( g2d, false, true, width, height, x, y );
        if ( type == GraphicsExportType.TIFF ) {
            writeToTiff( file, buffered_img );
        }
        else {
            ImageIO.write( buffered_img, type.toString(), file );
        }
        g2d.dispose();
        System.gc();
        if ( !options.isGraphicsExportUsingActualSize() ) {
            tree_panel.getMainPanel().getControlPanel().showWhole();
        }
        String msg = file.toString();
        if ( ( width > 0 ) && ( height > 0 ) ) {
            msg += " [size: " + width + ", " + height + "]";
        }
        return msg;
    }

    final static void writeToTiff( final File file, final BufferedImage image ) throws IOException {
        // See: http://log.robmeek.com/2005/08/write-tiff-in-java.html
        ImageWriter writer = null;
        ImageOutputStream ios = null;
        // Find an appropriate writer:
        final Iterator<ImageWriter> it = ImageIO.getImageWritersByFormatName( "TIF" );
        if ( it.hasNext() ) {
            writer = it.next();
        }
        else {
            throw new IOException( "failed to get TIFF image writer" );
        }
        // Setup writer:
        ios = ImageIO.createImageOutputStream( file );
        writer.setOutput( ios );
        final ImageWriteParam image_write_param = new ImageWriteParam( Locale.getDefault() );
        image_write_param.setCompressionMode( ImageWriteParam.MODE_EXPLICIT );
        // see writeParam.getCompressionTypes() for available compression type
        // strings.
        image_write_param.setCompressionType( "PackBits" );
        final String t[] = image_write_param.getCompressionTypes();
        for( final String string : t ) {
            System.out.println( string );
        }
        // Convert to an IIOImage:
        final IIOImage iio_image = new IIOImage( image, null, null );
        writer.write( null, iio_image, image_write_param );
    }

    private static void colorizeSubtree( final PhylogenyNode node, final BranchColor c ) {
        node.getBranchData().setBranchColor( c );
        final List<PhylogenyNode> descs = PhylogenyMethods.getAllDescendants( node );
        for( final PhylogenyNode desc : descs ) {
            desc.getBranchData().setBranchColor( c );
        }
    }

    final private static char normalizeCharForRGB( char c ) {
        c -= 65;
        c *= 10.2;
        c = c > 255 ? 255 : c;
        c = c < 0 ? 0 : c;
        return c;
    }

    final private static void openUrlInWebBrowser( final String url ) throws IOException, ClassNotFoundException,
            SecurityException, NoSuchMethodException, IllegalArgumentException, IllegalAccessException,
            InvocationTargetException, InterruptedException {
        final String os = System.getProperty( "os.name" );
        final Runtime runtime = Runtime.getRuntime();
        if ( os.toLowerCase().startsWith( "win" ) ) {
            Runtime.getRuntime().exec( "rundll32 url.dll,FileProtocolHandler " + url );
        }
        else if ( isMac() ) {
            final Class<?> file_mgr = Class.forName( "com.apple.eio.FileManager" );
            final Method open_url = file_mgr.getDeclaredMethod( "openURL", new Class[] { String.class } );
            open_url.invoke( null, new Object[] { url } );
        }
        else {
            final String[] browsers = { "firefox", "opera", "konqueror", "mozilla", "netscape", "epiphany" };
            String browser = null;
            for( int i = 0; ( i < browsers.length ) && ( browser == null ); ++i ) {
                if ( runtime.exec( new String[] { "which", browsers[ i ] } ).waitFor() == 0 ) {
                    browser = browsers[ i ];
                }
            }
            if ( browser == null ) {
                throw new IOException( "could not find a web browser to open [" + url + "] in" );
            }
            else {
                runtime.exec( new String[] { browser, url } );
            }
        }
    }

    // See: http://www.xml.nig.ac.jp/tutorial/rest/index.html#2.2
    // static void openDDBJRest() throws IOException {
    // //set URL
    // URL url = new URL( "http://xml.nig.ac.jp/rest/Invoke" );
    // //set parameter
    // String query = "service=GetEntry&method=getDDBJEntry&accession=AB000100";
    // //make connection
    // URLConnection urlc = url.openConnection();
    // //use post mode
    // urlc.setDoOutput( true );
    // urlc.setAllowUserInteraction( false );
    // //send query
    // PrintStream ps = new PrintStream( urlc.getOutputStream() );
    // ps.print( query );
    // ps.close();
    // //get result
    // BufferedReader br = new BufferedReader( new InputStreamReader(
    // urlc.getInputStream() ) );
    // String l = null;
    // while ( ( l = br.readLine() ) != null ) {
    // System.out.println( l );
    // }
    // br.close();
    // }
    public static enum GraphicsExportType {
        GIF( "gif" ), JPG( "jpg" ), PDF( "pdf" ), PNG( "png" ), TIFF( "tif" ), BMP( "bmp" );

        private final String _suffix;

        private GraphicsExportType( final String suffix ) {
            _suffix = suffix;
        }

        @Override
        public String toString() {
            return _suffix;
        }
    }
}
