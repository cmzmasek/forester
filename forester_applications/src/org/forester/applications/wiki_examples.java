// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

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