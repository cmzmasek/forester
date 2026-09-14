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

import java.awt.Color;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Reference;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * The re-rooting rules shared with Archaeopteryx.js ({@link Rerooting}, {@link AptxUtil#isTimeTree}): which trees
 * refuse a re-root, what counts as internal-node data, the prediction of which annotated nodes change clade (checked
 * against real re-roots), the warning text, and the root-free unrooted view.
 */
public final class RerootingTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "Rerooting: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return timeTree() && refusal() && nodeData() && cladeChangesMatchRealReroots() && warningText()
                    && rootFreeRule() && tipsAround();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    // the JS-pinned cases of forester.isTimeTree: a strict majority, and at least two, of the internal nodes (root
    // included) carry a date VALUE
    private static boolean timeTree() throws Exception {
        // root + i1 + i2 + i3 = four internal nodes, five tips
        final String nh = "((((A:1,B:1)i3:1,C:2)i2:1,D:3)i1:1,E:4)";
        if ( AptxUtil.isTimeTree( nhx( nh ) ) ) {
            return TestFail.here( "no dates" );
        }
        final Phylogeny tips = nhx( nh );
        for( final PhylogenyNode t : tips.getExternalNodes() ) {
            date( t );
        }
        if ( AptxUtil.isTimeTree( tips ) || ( AptxUtil.detectTimeTree( tips ) != AptxUtil.TIME_TREE_KIND.DATED ) ) {
            return TestFail.here( "tip dates alone are not a time tree for re-rooting (detectTimeTree still says DATED)" );
        }
        final Phylogeny one = nhx( nh );
        date( one.getRoot() );
        final Phylogeny two = nhx( nh );
        date( two.getRoot() );
        date( two.getNode( "i1" ) );
        final Phylogeny three = nhx( nh );
        date( three.getRoot() );
        date( three.getNode( "i1" ) );
        date( three.getNode( "i2" ) );
        if ( AptxUtil.isTimeTree( one ) || AptxUtil.isTimeTree( two ) || !AptxUtil.isTimeTree( three ) ) {
            return TestFail.here( "1 of 4 and 2 of 4 are not a strict majority, 3 of 4 is" );
        }
        final Phylogeny all = nhx( nh );
        final Phylogeny no_value = nhx( nh );
        for( final PhylogenyNodeIterator it = all.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isInternal() ) {
                date( n );
            }
        }
        for( final PhylogenyNodeIterator it = no_value.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isInternal() ) {
                n.getNodeData().setDate( new Date( "Late Cretaceous" ) ); // a description, no value
            }
        }
        if ( !AptxUtil.isTimeTree( all ) || AptxUtil.isTimeTree( no_value ) ) {
            return TestFail.here( "all internal dated is a time tree; dates without a value do not count" );
        }
        final Phylogeny two_tips = nhx( "(A:1,B:1)" );
        date( two_tips.getRoot() );
        if ( AptxUtil.isTimeTree( two_tips ) ) {
            return TestFail.here( "1 of 1 dated internal node is fewer than two" );
        }
        final Phylogeny balanced = nhx( "((A:1,B:1)x:1,(C:1,D:1)y:1)" );
        date( balanced.getNode( "x" ) );
        date( balanced.getNode( "y" ) );
        if ( !AptxUtil.isTimeTree( balanced ) || AptxUtil.isTimeTree( null ) || AptxUtil.isTimeTree( new Phylogeny() ) ) {
            return TestFail.here( "2 of 3 is a strict majority; null / empty are not time trees" );
        }
        return true;
    }

    private static boolean refusal() throws Exception {
        final String nh = "((A:1,B:1)x:1,(C:1,D:1)y:1)";
        if ( Rerooting.refusal( nhx( nh ) ) != null || Rerooting.refusal( null ) != null
                || Rerooting.refusal( new Phylogeny() ) != null ) {
            return TestFail.here( "an ordinary / null / empty tree may be re-rooted" );
        }
        final Phylogeny locked = nhx( nh );
        locked.setRerootable( false );
        final Phylogeny timed = nhx( nh );
        date( timed.getNode( "x" ) );
        date( timed.getNode( "y" ) );
        final Phylogeny both = timed.copy();
        both.setRerootable( false );
        final Phylogeny tip_dated = nhx( nh );
        for( final PhylogenyNode t : tip_dated.getExternalNodes() ) {
            date( t );
        }
        if ( ( Rerooting.refusal( locked ) != Rerooting.NOT_REROOTABLE )
                || ( Rerooting.refusal( timed ) != Rerooting.TIME_TREE )
                || ( Rerooting.refusal( both ) != Rerooting.NOT_REROOTABLE ) || ( Rerooting.refusal( tip_dated ) != null ) ) {
            return TestFail.here( "rerootable=false first, then time tree; tip dates alone do not refuse" );
        }
        if ( !Rerooting.NOT_REROOTABLE.equals( "This tree is marked as not re-rootable (rerootable=\"false\")." )
                || !Rerooting.TIME_TREE
                        .equals( "Time trees can't be re-rooted: their branch lengths are times measured from this root." ) ) {
            return TestFail.here( "the refusal texts are a joint contract with Archaeopteryx.js" );
        }
        return true;
    }

    private static boolean nodeData() throws Exception {
        // not data: branch length, support, branch colour, node visual style
        final PhylogenyNode bare = internal();
        bare.setDistanceToParent( 0.5 );
        bare.getBranchData().addConfidence( new Confidence( 90, "bootstrap" ) );
        bare.getBranchData().setBranchColor( new BranchColor( Color.RED ) );
        bare.getNodeData().setNodeVisualData( new NodeVisualData() );
        if ( Rerooting.hasNodeData( bare ) ) {
            return TestFail.here( "branch length, support and visual styling are not node data" );
        }
        final PhylogenyNode empty_tax = internal();
        empty_tax.getNodeData().setTaxonomy( new Taxonomy() );
        if ( Rerooting.hasNodeData( empty_tax ) ) {
            return TestFail.here( "an empty taxonomy is not data" );
        }
        final List<PhylogenyNode> data = new ArrayList<>();
        final PhylogenyNode named = internal();
        named.setName( "Mammalia" );
        data.add( named );
        final PhylogenyNode tax = internal();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( "Primates" );
        tax.getNodeData().setTaxonomy( t );
        data.add( tax );
        final PhylogenyNode seq = internal();
        seq.getNodeData().setSequence( new Sequence() );
        data.add( seq );
        final PhylogenyNode ev = internal();
        ev.getNodeData().setEvent( new Event( 1, 0, 0 ) );
        data.add( ev );
        final PhylogenyNode dist = internal();
        dist.getNodeData().setDistribution( new Distribution( "Africa" ) );
        data.add( dist );
        final PhylogenyNode dated = internal();
        dated.getNodeData().setDate( new Date( "Miocene" ) );
        data.add( dated );
        final PhylogenyNode ref = internal();
        ref.getNodeData().setReference( new Reference( "doi:10.1/x" ) );
        data.add( ref );
        final PhylogenyNode bin = internal();
        bin.getNodeData().setBinaryCharacters( new BinaryCharacters() );
        data.add( bin );
        data.add( withProperty( "data:host", AppliesTo.NODE ) );
        data.add( withProperty( "data:host", AppliesTo.CLADE ) );
        for( final PhylogenyNode n : data ) {
            if ( !Rerooting.hasNodeData( n ) ) {
                return TestFail.here( "must count as node data: " + n.getNodeData() );
            }
        }
        if ( Rerooting.hasNodeData( withProperty( "aptx:reimport_profile", AppliesTo.NODE ) )
                || Rerooting.hasNodeData( withProperty( NodeVisualData.APTX_VISUALIZATION_REF + "font", AppliesTo.NODE ) )
                || Rerooting.hasNodeData( withProperty( "data:rate", AppliesTo.PARENT_BRANCH ) ) ) {
            return TestFail.here( "internal aptx:, style: and parent-branch properties are not node data" );
        }
        // internalNodesWithData: internal nodes only, the root included
        final Phylogeny phy = nhx( "((A:1,B:1)x:1,(C:1,D:1):1)" );
        if ( Rerooting.internalNodesWithData( phy ) != 1 ) {
            return TestFail.here( "named tips do not count; one named internal node does" );
        }
        phy.getRoot().setName( "root" );
        if ( ( Rerooting.internalNodesWithData( phy ) != 2 ) || ( Rerooting.internalNodesWithData( null ) != 0 ) ) {
            return TestFail.here( "the root counts" );
        }
        return true;
    }

    /**
     * The shortcut "clade changed == set of children changed, or node gone" must equal what a real re-root does to
     * each annotated node's tip set -- over random trees, every manual re-root position, midpoint and MAD.
     */
    private static boolean cladeChangesMatchRealReroots() {
        final Random r = new Random( 11 );
        int changed_cases = 0;
        int unchanged_cases = 0;
        for( int t = 0; t < 40; ++t ) {
            final Phylogeny phy = randomTree( r, 4 + r.nextInt( 11 ) );
            final List<Long> targets = new ArrayList<>();
            for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isRoot() ) {
                    targets.add( n.getId() );
                }
            }
            final List<Phylogeny> rerooted = new ArrayList<>();
            for( final long id : targets ) {
                final Phylogeny copy = phy.copy();
                copy.reRoot( copy.getNode( id ) );
                rerooted.add( copy );
            }
            final Phylogeny mid = phy.copy();
            PhylogenyMethods.midpointRoot( mid );
            rerooted.add( mid );
            final Phylogeny mad = phy.copy();
            PhylogenyMethods.madRoot( mad );
            rerooted.add( mad );
            for( final Phylogeny after : rerooted ) {
                final int predicted = Rerooting.dataNodesWhoseCladeChanges( phy, after );
                final int actual = annotatedNodesWhoseTipSetChanged( phy, after );
                if ( predicted != actual ) {
                    return TestFail.here( "tree " + t + ": predicted " + predicted + " changed clades, actually " + actual );
                }
                if ( actual > 0 ) {
                    ++changed_cases;
                }
                else {
                    ++unchanged_cases;
                }
            }
        }
        if ( ( changed_cases < 50 ) || ( unchanged_cases < 20 ) ) {
            return TestFail.here( "the fixture must reach both outcomes: " + changed_cases + " / " + unchanged_cases );
        }
        return true;
    }

    // independent of Rerooting's shortcut: compares each annotated internal node's actual set of tip names
    private static int annotatedNodesWhoseTipSetChanged( final Phylogeny before, final Phylogeny after ) {
        final Map<Long, PhylogenyNode> by_id = new HashMap<>();
        for( final PhylogenyNodeIterator it = after.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            by_id.put( n.getId(), n );
        }
        int count = 0;
        for( final PhylogenyNodeIterator it = before.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isInternal() && Rerooting.hasNodeData( n ) ) {
                final PhylogenyNode m = by_id.get( n.getId() );
                if ( ( m == null ) || !tipNames( n ).equals( tipNames( m ) ) ) {
                    ++count;
                }
            }
        }
        return count;
    }

    private static Set<String> tipNames( final PhylogenyNode n ) {
        final Set<String> names = new HashSet<>();
        for( final PhylogenyNode tip : n.getAllExternalDescendants() ) {
            names.add( tip.getName() );
        }
        return names;
    }

    // a random tree with branch lengths; about half the internal nodes (the root included) are named, i.e. carry data;
    // the root has two or three children
    private static Phylogeny randomTree( final Random r, final int n_tips ) {
        final List<PhylogenyNode> active = new ArrayList<>();
        for( int i = 0; i < n_tips; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "T" + i );
            leaf.setDistanceToParent( 0.05 + r.nextDouble() );
            active.add( leaf );
        }
        final int root_degree = 2 + r.nextInt( 2 );
        int k = 0;
        while ( active.size() > root_degree ) {
            final PhylogenyNode x = active.remove( r.nextInt( active.size() ) );
            final PhylogenyNode y = active.remove( r.nextInt( active.size() ) );
            final PhylogenyNode parent = new PhylogenyNode();
            parent.addAsChild( x );
            parent.addAsChild( y );
            parent.setDistanceToParent( 0.05 + r.nextDouble() );
            if ( r.nextBoolean() ) {
                parent.setName( "N" + k++ );
            }
            active.add( parent );
        }
        final PhylogenyNode root = new PhylogenyNode();
        for( final PhylogenyNode n : active ) {
            root.addAsChild( n );
        }
        if ( r.nextInt( 3 ) == 0 ) {
            root.setName( "ROOT" );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean warningText() {
        if ( !Rerooting.dataWarning( 12, 4 ).equals( "This tree has data on 12 internal nodes. Re-rooting changes the "
                + "clade of 4 of them, so their data may no longer describe them." ) ) {
            return TestFail.here( Rerooting.dataWarning( 12, 4 ) );
        }
        if ( !Rerooting.dataWarning( 5, 1 ).equals( "This tree has data on 5 internal nodes. Re-rooting changes the "
                + "clade of 1 of them, so its data may no longer describe it." ) ) {
            return TestFail.here( Rerooting.dataWarning( 5, 1 ) );
        }
        if ( !Rerooting.dataWarning( 1, 1 ).equals( "This tree has data on 1 internal node. Re-rooting changes its "
                + "clade, so its data may no longer describe it." ) ) {
            return TestFail.here( Rerooting.dataWarning( 1, 1 ) );
        }
        return true;
    }

    private static boolean rootFreeRule() throws Exception {
        final Phylogeny declared = nhx( "((A:1,B:1):1,C:2)" );
        declared.setRooted( false );
        declared.setRootednessDeclared( true );
        final Phylogeny newick = nhx( "((A:1,B:1):1,C:2)" ); // unrooted by default, but declares nothing
        final Phylogeny rooted = nhx( "((A:1,B:1):1,C:2)" );
        rooted.setRooted( true );
        rooted.setRootednessDeclared( true );
        final Options.PHYLOGENY_GRAPHICS_TYPE unrooted = Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED;
        if ( !Rerooting.hidesRootDependentValues( unrooted, declared ) ) {
            return TestFail.here( "declared unrooted AND the unrooted layout hides root-dependent values" );
        }
        if ( Rerooting.hidesRootDependentValues( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, declared )
                || Rerooting.hidesRootDependentValues( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, declared )
                || Rerooting.hidesRootDependentValues( unrooted, newick )
                || Rerooting.hidesRootDependentValues( unrooted, rooted )
                || Rerooting.hidesRootDependentValues( unrooted, null ) ) {
            return TestFail.here( "both conditions are needed" );
        }
        return true;
    }

    private static boolean tipsAround() throws Exception {
        final Phylogeny phy = nhx( "((A,B)x,C,D)" );
        final PhylogenyNode x = phy.getNode( "x" );
        if ( !Rerooting.tipsAround( x ).equals( List.of( 1, 1, 2 ) )
                || !Rerooting.tipsAround( phy.getRoot() ).equals( List.of( 1, 1, 2 ) )
                || !Rerooting.tipsAround( phy.getNode( "A" ) ).isEmpty()
                || !Rerooting.tipsAroundText( x ).equals( "1 · 1 · 2" ) ) {
            return TestFail.here( "sides of x: A, B and the rest (C, D); the stored root: its children's sides" );
        }
        final Phylogeny poly = nhx( "((A,B,C)p,(D,E,F,G)q)" );
        if ( !Rerooting.tipsAroundText( poly.getNode( "p" ) ).equals( "1 · 1 · 1 · 4" )
                || !Rerooting.tipsAroundText( poly.getNode( "q" ) ).equals( "1 · 1 · 1 · 1 · 3" ) ) {
            return TestFail.here( "a polytomy has one side per neighbour, ascending" );
        }
        return true;
    }

    private static Phylogeny nhx( final String nh ) throws Exception {
        return ParserBasedPhylogenyFactory.getInstance().create( nh, new NHXParser() )[ 0 ];
    }

    private static void date( final PhylogenyNode n ) {
        n.getNodeData().setDate( new Date( "", BigDecimal.ONE, null, null, "mya" ) );
    }

    private static PhylogenyNode internal() {
        final PhylogenyNode n = new PhylogenyNode();
        n.addAsChild( new PhylogenyNode() );
        return n;
    }

    private static PhylogenyNode withProperty( final String ref, final AppliesTo applies_to ) {
        final PhylogenyNode n = internal();
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( ref, "v", "", "xsd:string", applies_to ) );
        n.getNodeData().setProperties( pl );
        return n;
    }

    private RerootingTest() {
    }
}
