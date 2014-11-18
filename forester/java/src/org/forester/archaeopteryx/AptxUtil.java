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
import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URI;
import java.net.URL;
import java.net.URLConnection;
import java.text.ParseException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;
import javax.swing.JApplet;
import javax.swing.JOptionPane;
import javax.swing.text.MaskFormatter;

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
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.AsciiHistogram;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public final class AptxUtil {

    private final static String[] AVAILABLE_FONT_FAMILIES_SORTED = GraphicsEnvironment.getLocalGraphicsEnvironment()
                                                                         .getAvailableFontFamilyNames();
    static {
        Arrays.sort( AVAILABLE_FONT_FAMILIES_SORTED );
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

    final static public boolean isHasAtLeastOneBranchWithSupportSD( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( n.getBranchData().isHasConfidences() ) {
                final List<Confidence> c = n.getBranchData().getConfidences();
                for( final Confidence confidence : c ) {
                    if ( confidence.getStandardDeviation() > 0 ) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    final static public boolean isHasAtLeastOneNodeWithScientificName( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy()
                    && !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                return true;
            }
        }
        return false;
    }

    final static public boolean isHasAtLeastOneNodeWithSequenceAnnotation( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasSequence()
                    && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getAnnotations() ) ) {
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

    public final static void printWarningMessage( final String name, final String message ) {
        System.out.println( "[" + name + "] > " + message );
    }

    final public static void showErrorMessage( final Component parent, final String error_msg ) {
        printAppletMessage( Constants.PRG_NAME, error_msg );
        JOptionPane.showMessageDialog( parent, error_msg, "[" + Constants.PRG_NAME + " " + Constants.VERSION
                + "] Error", JOptionPane.ERROR_MESSAGE );
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
        tree_panel.calcParametersForPainting( width, height );
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

    // Returns true if the specified format name can be written
    final static boolean canWriteFormat( final String format_name ) {
        final Iterator<ImageWriter> iter = ImageIO.getImageWritersByFormatName( format_name );
        return iter.hasNext();
    }

    final static String createBasicInformation( final Phylogeny phy, final File treefile ) {
        final StringBuilder desc = new StringBuilder();
        if ( ( phy != null ) && !phy.isEmpty() ) {
            String f = null;
            if ( treefile != null ) {
                try {
                    f = treefile.getCanonicalPath();
                }
                catch ( final IOException e ) {
                    //Not important, ignore.
                }
                if ( !ForesterUtil.isEmpty( f ) ) {
                    desc.append( "Path: " );
                    desc.append( f );
                    desc.append( "\n" );
                }
            }
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
            for( final Taxonomy t : taxs ) {
                System.out.println( t.toString() );
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
                desc.append( "    Mean: " + ForesterUtil.round( bs.arithmeticMean(), 6 ) + " (stdev: "
                        + ForesterUtil.round( bs.sampleStandardDeviation(), 6 ) + ")" );
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
                desc.append( "    Median: " + ForesterUtil.round( ds.median(), 6 ) );
                desc.append( "\n" );
                desc.append( "    Mean: " + ForesterUtil.round( ds.arithmeticMean(), 6 ) + " (stdev: "
                        + ForesterUtil.round( ds.sampleStandardDeviation(), 6 ) + ")" );
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
                        if ( cs.getN() > 2 ) {
                            desc.append( " (stdev: " + ForesterUtil.round( cs.sampleStandardDeviation(), 6 ) + ")" );
                        }
                        desc.append( "\n" );
                        desc.append( "    Minimum: " + ForesterUtil.round( cs.getMin(), 6 ) );
                        desc.append( "\n" );
                        desc.append( "    Maximum: " + ForesterUtil.round( cs.getMax(), 6 ) );
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

    final public static Phylogeny[] readPhylogeniesFromUrl( final URL url,
                                                            final boolean phyloxml_validate_against_xsd,
                                                            final boolean replace_underscores,
                                                            final boolean internal_numbers_are_confidences,
                                                            final TAXONOMY_EXTRACTION taxonomy_extraction,
                                                            final boolean midpoint_reroot )
            throws FileNotFoundException, IOException {
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
        AptxUtil.printAppletMessage( "Archaeopteryx", "parser is " + parser.getName() );
        final URLConnection url_connection = url.openConnection();
        url_connection.setDefaultUseCaches( false );
        final InputStream i = url_connection.getInputStream();
        final Phylogeny[] phys = factory.create( i, parser );
        i.close();
        if ( phys != null ) {
            if ( nhx_or_nexus && internal_numbers_are_confidences ) {
                for( final Phylogeny phy : phys ) {
                    PhylogenyMethods.transferInternalNodeNamesToConfidence( phy, "" );
                }
            }
            if ( midpoint_reroot ) {
                for( final Phylogeny phy : phys ) {
                    PhylogenyMethods.midpointRoot( phy );
                    PhylogenyMethods.orderAppearance( phy.getRoot(), true, true, DESCENDANT_SORT_PRIORITY.NODE_NAME );
                }
            }
        }
        return phys;
    }

    final static void removeBranchColors( final Phylogeny phy ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            it.next().getBranchData().setBranchColor( null );
        }
    }

    final static void removeVisualStyles( final Phylogeny phy ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            it.next().getNodeData().setNodeVisualData( null );
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
            tree_panel.calcParametersForPainting( options.getPrintSizeX(), options.getPrintSizeY() );
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
            tree_panel.calcParametersForPainting( options.getPrintSizeX(), options.getPrintSizeY() );
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

    final private static void openUrlInWebBrowser( final String url ) throws IOException, ClassNotFoundException,
            SecurityException, NoSuchMethodException, IllegalArgumentException, IllegalAccessException,
            InvocationTargetException, InterruptedException {
        final String os = System.getProperty( "os.name" );
        final Runtime runtime = Runtime.getRuntime();
        if ( os.toLowerCase().startsWith( "win" ) ) {
            Runtime.getRuntime().exec( "rundll32 url.dll,FileProtocolHandler " + url );
        }
        else if ( ForesterUtil.isMac() ) {
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

    public static enum GraphicsExportType {
        BMP( "bmp" ), GIF( "gif" ), JPG( "jpg" ), PDF( "pdf" ), PNG( "png" ), TIFF( "tif" );

        private final String _suffix;

        private GraphicsExportType( final String suffix ) {
            _suffix = suffix;
        }

        @Override
        public String toString() {
            return _suffix;
        }
    }

    final public static Color calculateColorFromString( final String str, final boolean is_taxonomy ) {
        final String my_str = str.toUpperCase();
        char first = my_str.charAt( 0 );
        char second = ' ';
        char third = ' ';
        if ( my_str.length() > 1 ) {
            if ( is_taxonomy ) {
                second = my_str.charAt( 1 );
            }
            else {
                second = my_str.charAt( my_str.length() - 1 );
            }
            if ( is_taxonomy ) {
                if ( my_str.length() > 2 ) {
                    if ( my_str.indexOf( " " ) > 0 ) {
                        third = my_str.charAt( my_str.indexOf( " " ) + 1 );
                    }
                    else {
                        third = my_str.charAt( 2 );
                    }
                }
            }
            else if ( my_str.length() > 2 ) {
                third = my_str.charAt( ( my_str.length() - 1 ) / 2 );
            }
        }
        first = normalizeCharForRGB( first );
        second = normalizeCharForRGB( second );
        third = normalizeCharForRGB( third );
        if ( ( first > 235 ) && ( second > 235 ) && ( third > 235 ) ) {
            first = 0;
        }
        else if ( ( first < 60 ) && ( second < 60 ) && ( third < 60 ) ) {
            second = 255;
        }
        return new Color( first, second, third );
    }

    final private static char normalizeCharForRGB( char c ) {
        c -= 65;
        c *= 10.2;
        c = c > 255 ? 255 : c;
        c = c < 0 ? 0 : c;
        return c;
    }
}
