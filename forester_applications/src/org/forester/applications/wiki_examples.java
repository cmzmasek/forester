
package org.forester.applications;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.forester.archaeopteryx.Archaeopteryx;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class wiki_examples {

    public static void main( final String[] args ) {
        // Reading-in of (a) tree(s) from a file.
        final File treefile = new File( args[ 0 ] );
        PhylogenyParser parser = null;
        try {
            parser = ParserUtils.createParserDependingOnFileType( treefile, true );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        Phylogeny[] phys = null;
        try {
            phys = PhylogenyMethods.readPhylogenies( parser, treefile );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        final Phylogeny phy = phys[ 0 ];
        // Read node->color map into a map
        final Map<String, Color> colors = new HashMap<String, Color>();
        // read it in from file...
        // Iterate over nodes and set colors from 'colors' map
        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
            // if node-name (?) in 'colors' map
            it.next().getBranchData().setBranchColor( new BranchColor( colors.get( "xx" ) ) );
        }
        // For testing, use Aptx...
        Archaeopteryx.createApplication( phy );
        // Finally, create
    }
}