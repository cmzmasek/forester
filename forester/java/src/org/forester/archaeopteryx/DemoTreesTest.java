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

package org.forester.archaeopteryx;

import java.io.File;
import java.util.Iterator;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Anti-rot guard for the {@code forester/demo/} gallery: every committed demo tree (see
 * {@link org.forester.demo.DemoTreeGenerator}) must still parse cleanly AND still carry the data its feature needs to
 * be demonstrable. Headless -- pure phyloXML parsing + model inspection, no GUI. If a parser or the property/taxonomy
 * model changes in a way that breaks a demo, this fails in the always-run suite instead of the demo silently rotting.
 *
 * <p>When a new demo tree is added (per the "one demo tree per new feature" convention), add a shape check here.
 */
public final class DemoTreesTest {

    private static final String DEMO_DIR = System.getProperty( "user.dir" ) + File.separator + "forester"
            + File.separator + "demo" + File.separator;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DemoTrees: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        boolean ok = true;
        // size-by: a numeric property to scale by
        ok &= hasNumericRef( "size-by-property.xml", "data:read_count" );
        // color-by: a categorical property (palette) AND a numeric one (gradient)
        ok &= hasCategoricalRef( "color-by-property.xml", "data:host" );
        ok &= hasNumericRef( "color-by-property.xml", "data:year" );
        // annotation columns: categorical + numeric + text properties
        ok &= hasCategoricalRef( "annotation-columns.xml", "data:host" );
        ok &= hasNumericRef( "annotation-columns.xml", "data:viral_load" );
        ok &= hasCategoricalRef( "annotation-columns.xml", "data:clade" );
        // colorize by rank / clade bands: in-tree 'order' rank annotation so the demo works offline
        ok &= hasRank( "colorize-by-rank.xml", "order" );
        ok &= hasRank( "colorize-by-rank.xml", "species" );
        // scale axis: a phylogram with real branch lengths (so the labeled distance axis has ticks)
        ok &= hasBranchLengths( "scale-axis.xml" );
        // node age bars: internal nodes carry a phyloXML <date> with a min/max interval + branch lengths (dated tree)
        ok &= hasBranchLengths( "node-hpd-bars.xml" );
        ok &= hasInternalDateInterval( "node-hpd-bars.xml" );
        return ok;
    }

    private static boolean hasInternalDateInterval( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getMin() != null )
                    && ( n.getNodeData().getDate().getMax() != null ) ) {
                return true;
            }
        }
        return note( file_name + " must carry an INTERNAL-node <date> with a min/max interval (for HPD age bars)" );
    }

    private static boolean hasBranchLengths( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( !( org.forester.phylogeny.PhylogenyMethods.calculateMaxDistanceToRoot( phy ) > 0.0 ) ) {
            return note( file_name + " must carry branch lengths (a positive max distance to root) for a scale axis" );
        }
        return true;
    }

    /** Parses a demo file, asserting it is present, error-free and non-empty; null on any failure (with a message). */
    private static Phylogeny load( final String file_name ) {
        final File file = new File( DEMO_DIR + file_name );
        if ( !file.exists() ) {
            return fail( file_name + " is missing from the demo gallery (" + file.getAbsolutePath() + ")" );
        }
        try {
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny[] phys = ParserBasedPhylogenyFactory.getInstance().create( file, parser );
            if ( parser.getErrorCount() > 0 ) {
                return fail( file_name + " parsed with errors: " + parser.getErrorMessages() );
            }
            if ( ( phys == null ) || ( phys.length != 1 ) || ( phys[ 0 ] == null ) || phys[ 0 ].isEmpty() ) {
                return fail( file_name + " did not yield exactly one non-empty tree" );
            }
            return phys[ 0 ];
        }
        catch ( final Exception e ) {
            return fail( file_name + " could not be read: " + e.getMessage() );
        }
    }

    private static boolean hasNumericRef( final String file_name, final String ref ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( !PropertyColorScheme.numericRefs( phy ).contains( ref ) ) {
            return note( file_name + " must offer numeric property '" + ref + "' for its feature demo" );
        }
        return true;
    }

    private static boolean hasCategoricalRef( final String file_name, final String ref ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        boolean present = false;
        for( final PhylogenyNode leaf : phy.getExternalNodes() ) {
            if ( PropertyColorScheme.valueFor( leaf, ref ) != null ) {
                present = true;
                break;
            }
        }
        if ( !present ) {
            return note( file_name + " must carry property '" + ref + "' on at least one tip" );
        }
        // 'categorical' means it is NOT treated as a numeric gradient
        if ( PropertyColorScheme.numericRefs( phy ).contains( ref ) ) {
            return note( file_name + " property '" + ref + "' should be categorical, not numeric" );
        }
        return true;
    }

    private static boolean hasRank( final String file_name, final String rank ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() && rank.equals( n.getNodeData().getTaxonomy().getRank() ) ) {
                return true;
            }
        }
        return note( file_name + " must carry an in-tree taxonomy at rank '" + rank + "' (offline colorize)" );
    }

    /** null-returning failure note for the parse path (so load(...) can bail with a message). */
    private static Phylogeny fail( final String msg ) {
        note( msg );
        return null;
    }

    /** false-returning failure note for the shape checks. */
    private static boolean note( final String msg ) {
        System.out.println( "  [DemoTreesTest] " + msg );
        return false;
    }

    private DemoTreesTest() {
    }
}
