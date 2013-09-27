
package org.forester.applications;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.forester.archaeopteryx.AptxUtil;
import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.archaeopteryx.Configuration;
import org.forester.archaeopteryx.Options;
import org.forester.archaeopteryx.TreeColorSet;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class phylo2coloredgraphics {

    public static void main( final String[] args ) {
        try {
            // Reading-in of a tree from a file.
            final File treefile = new File( "/home/czmasek/tol_117_TEST.xml" );
            final PhylogenyParser parser = ParserUtils.createParserDependingOnFileType( treefile, true );
            final Phylogeny phy = PhylogenyMethods.readPhylogenies( parser, treefile )[ 0 ];
            // Creating a node name -> color map.
            final Map<String, Color> colors = new HashMap<String, Color>();
            colors.put( "Primates", new Color( 255, 255, 0 ) );
            colors.put( "PANTR", new Color( 255, 0, 255 ) );
            colors.put( "HUMAN", new Color( 255, 0, 0 ) );
            colors.put( "RAT", new Color( 155, 0, 0 ) );
            colors.put( "MOUSE", new Color( 55, 155, 0 ) );
            colors.put( "CAVPO", new Color( 155, 155, 0 ) );
            colors.put( "LOTGI", new Color( 155, 155, 255 ) );
            // Setting colors.
            for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( colors.containsKey( n.getName() ) ) {
                    n.getBranchData().setBranchColor( new BranchColor( colors.get( n.getName() ) ) );
                    // To make colored subtrees thicker:
                    n.getBranchData().setBranchWidth( new BranchWidth( 4 ) );
                }
            }
            // Setting up a configuration object.
            final Configuration config = new Configuration();
            config.putDisplayColors( TreeColorSet.BACKGROUND, new Color( 255, 255, 255 ) );
            config.putDisplayColors( TreeColorSet.BRANCH, new Color( 0, 0, 0 ) );
            config.putDisplayColors( TreeColorSet.TAXONOMY, new Color( 0, 0, 0 ) );
            config.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            config.setTaxonomyColorize( false );
            config.setColorizeBranches( true );
            config.setUseBranchesWidths( true );
            config.setDisplayTaxonomyCode( false );
            // Writing to a graphics file.
            AptxUtil.writePhylogenyToGraphicsFile( phy,
                                                   new File( "/home/czmasek/000.png" ),
                                                   1300,
                                                   1300,
                                                   GraphicsExportType.PNG,
                                                   config );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }
}
