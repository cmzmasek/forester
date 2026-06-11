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

import java.io.File;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

// javac -cp ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/domainloss_replacement.java
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.domainloss_replacement
public class domainloss_replacement {

    public static void main( final String args[] ) {
        try {
            if ( args.length != 2 ) {
                System.out
                .println( "Usage: domainloss_replacement <phylogeny file> <file with replacement characters>" );
                System.exit( -1 );
            }
            final Phylogeny p = ParserUtils.readPhylogenies( args[ 0 ] )[ 0 ];
            final Set<String> replacement_domains = ForesterUtil.file2set( new File( args[ 1 ] ) );
            for( final PhylogenyNodeIterator it = p.iteratorExternalForward(); it.hasNext(); ) {
                PhylogenyNode n = it.next();
                String name = null;
                if ( n.getNodeData().isHasTaxonomy()
                        && !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                    name = n.getNodeData().getTaxonomy().getScientificName();
                }
                else {
                    name = n.getName();
                }
                final SortedSet<String> lost_chars = new TreeSet<String>();
                while ( !n.isRoot() ) {
                    lost_chars.addAll( n.getNodeData().getBinaryCharacters().getLostCharacters() );
                    n = n.getParent();
                }
                final int losses = lost_chars.size();
                lost_chars.retainAll( replacement_domains );
                final int intersection = lost_chars.size();
                final double percentage = ( 100.0 * intersection ) / losses;
                System.out.println( name + "\t" + intersection + "\t" + losses + "\t"
                        + ForesterUtil.round( percentage, 3 ) );
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
    }
}
