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

import java.util.List;
import java.util.Set;

import org.forester.archaeopteryx.NodeDataDraft.Problem;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;

/**
 * Headless tests for {@link TreePropertiesDraft}: reading a tree's six metadata fields, normalization (name
 * whitespace collapsed, the rest trimmed), validation (a named tree cannot be blanked; a provider needs a value),
 * change detection, and the write-back (a blank identifier removes it, a blank provider is stored as none).
 */
public final class TreePropertiesDraftTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreePropertiesDraft: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return read() && normalization() && validation() && changes() && write() && equality();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "A" );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "B" );
        root.addAsChild( a );
        root.addAsChild( b );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.setName( "My tree" );
        phy.setDescription( "line one\n  line two" );
        phy.setIdentifier( new Identifier( "TB2:Tr1", "treebase" ) );
        phy.setType( "gene tree" );
        phy.setDistanceUnit( "substitutions/site" );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean read() {
        final TreePropertiesDraft d = TreePropertiesDraft.from( tree() );
        if ( !"My tree".equals( d.name ) || !"line one\n  line two".equals( d.description )
                || !"TB2:Tr1".equals( d.idValue ) || !"treebase".equals( d.idProvider ) || !"gene tree".equals( d.type )
                || !"substitutions/site".equals( d.distanceUnit ) ) {
            return TestFail.here( d.toString() );
        }
        // nulls read as "", including a tree without an identifier and a null tree
        final Phylogeny bare = new Phylogeny();
        bare.setRoot( new PhylogenyNode() );
        final TreePropertiesDraft e = TreePropertiesDraft.from( bare );
        if ( !e.name.isEmpty() || !e.idValue.isEmpty() || !e.idProvider.isEmpty() || !e.type.isEmpty()
                || !e.distanceUnit.isEmpty() ) {
            return TestFail.here( e.toString() );
        }
        if ( !TreePropertiesDraft.from( null ).equals( new TreePropertiesDraft() ) ) {
            return TestFail.here();
        }
        // copy() is a deep-equal, independent object
        final TreePropertiesDraft c = d.copy();
        if ( !c.equals( d ) || ( c == d ) ) {
            return TestFail.here();
        }
        c.name = "other";
        if ( c.equals( d ) ) {
            return TestFail.here( "copy must be independent" );
        }
        return true;
    }

    private static boolean normalization() {
        final TreePropertiesDraft d = new TreePropertiesDraft();
        d.name = "  multi   word\tname  ";
        d.description = "  para one\n\n  indented   two  ";
        d.idValue = " id ";
        d.idProvider = " prov ";
        d.type = "  gene\t tree ";
        d.distanceUnit = " Ma ";
        final TreePropertiesDraft n = d.normalized();
        if ( !"multi word name".equals( n.name ) ) {
            return TestFail.here( "name: " + n.name );
        }
        if ( !"para one\n\n  indented   two".equals( n.description ) ) {
            return TestFail.here( "description must keep internal whitespace: " + n.description );
        }
        if ( !"id".equals( n.idValue ) || !"prov".equals( n.idProvider ) || !"gene tree".equals( n.type )
                || !"Ma".equals( n.distanceUnit ) ) {
            return TestFail.here( n.toString() );
        }
        if ( n.equals( d ) ) {
            return TestFail.here( "raw and normalized differ" );
        }
        if ( !n.normalized().equals( n ) ) {
            return TestFail.here( "normalization is idempotent" );
        }
        return true;
    }

    private static boolean validation() {
        final TreePropertiesDraft base = TreePropertiesDraft.from( tree() );
        if ( !base.validate( base ).isEmpty() ) {
            return TestFail.here( "an unchanged draft validates" );
        }
        // blanking the name of a NAMED tree is a problem, keyed to the name field
        final TreePropertiesDraft blank = base.copy();
        blank.name = "   ";
        final List<Problem> ps = blank.validate( base );
        if ( ( ps.size() != 1 ) || !TreePropertiesDraft.NAME.equals( ps.get( 0 ).key ) ) {
            return TestFail.here( ps.toString() );
        }
        // ... but not of an UNNAMED tree (nothing to lose), nor without a baseline
        final TreePropertiesDraft unnamed = new TreePropertiesDraft();
        if ( !blank.validate( unnamed ).isEmpty() || !blank.validate( null ).isEmpty() ) {
            return TestFail.here();
        }
        // a provider without an identifier value
        final TreePropertiesDraft prov = base.copy();
        prov.idValue = " ";
        final List<Problem> pp = prov.validate( base );
        if ( ( pp.size() != 1 ) || !TreePropertiesDraft.ID_PROVIDER.equals( pp.get( 0 ).key ) ) {
            return TestFail.here( pp.toString() );
        }
        // both at once: two problems, name first
        final TreePropertiesDraft both = prov.copy();
        both.name = "";
        final List<Problem> pb = both.validate( base );
        if ( ( pb.size() != 2 ) || !TreePropertiesDraft.NAME.equals( pb.get( 0 ).key )
                || !TreePropertiesDraft.ID_PROVIDER.equals( pb.get( 1 ).key ) ) {
            return TestFail.here( pb.toString() );
        }
        // a value without a provider is fine; type / unit are free text
        final TreePropertiesDraft free = base.copy();
        free.idProvider = "";
        free.type = "anything at all";
        free.distanceUnit = "furlongs";
        if ( !free.validate( base ).isEmpty() ) {
            return TestFail.here();
        }
        return true;
    }

    private static boolean changes() {
        final TreePropertiesDraft base = TreePropertiesDraft.from( tree() );
        if ( !base.changedFields( base ).isEmpty() ) {
            return TestFail.here();
        }
        // whitespace-only differences are NOT changes
        final TreePropertiesDraft ws = base.copy();
        ws.name = "  My   tree ";
        ws.description = base.description + "   ";
        ws.type = "gene  tree";
        if ( !ws.changedFields( base ).isEmpty() ) {
            return TestFail.here( ws.changedFields( base ).toString() );
        }
        // every field, reported with its user-facing label, in field order
        final TreePropertiesDraft all = base.copy();
        all.name = "x";
        all.description = "y";
        all.idProvider = "z";
        all.type = "t";
        all.distanceUnit = "u";
        final Set<String> ch = all.changedFields( base );
        if ( !ch.toString().equals( "[name, description, identifier, type, branch-length unit]" ) ) {
            return TestFail.here( ch.toString() );
        }
        // the identifier counts as one field whichever half changed
        final TreePropertiesDraft idv = base.copy();
        idv.idValue = "other";
        if ( !idv.changedFields( base ).toString().equals( "[identifier]" ) ) {
            return TestFail.here( idv.changedFields( base ).toString() );
        }
        return true;
    }

    private static boolean write() {
        final Phylogeny phy = tree();
        final TreePropertiesDraft d = TreePropertiesDraft.from( phy );
        d.name = "  Renamed   tree ";
        d.description = "  new text  ";
        d.idValue = " 42 ";
        d.idProvider = "  "; // blank provider -> stored as none
        d.type = " species tree ";
        d.distanceUnit = " Ma ";
        d.writeTo( phy );
        if ( !"Renamed tree".equals( phy.getName() ) || !"new text".equals( phy.getDescription() )
                || !"species tree".equals( phy.getType() ) || !"Ma".equals( phy.getDistanceUnit() ) ) {
            return TestFail.here( phy.getName() + "|" + phy.getDescription() + "|" + phy.getType() + "|"
                    + phy.getDistanceUnit() );
        }
        if ( ( phy.getIdentifier() == null ) || !"42".equals( phy.getIdentifier().getValue() )
                || ( phy.getIdentifier().getProvider() != null ) ) {
            return TestFail.here( String.valueOf( phy.getIdentifier() ) );
        }
        // reading back gives the normalized draft
        if ( !TreePropertiesDraft.from( phy ).equals( d.normalized() ) ) {
            return TestFail.here( TreePropertiesDraft.from( phy ).toString() );
        }
        // a blank identifier value removes the identifier entirely; blank type / unit clear them
        final TreePropertiesDraft clear = TreePropertiesDraft.from( phy );
        clear.idValue = "";
        clear.idProvider = "ignored";
        clear.type = "";
        clear.distanceUnit = "";
        clear.writeTo( phy );
        if ( ( phy.getIdentifier() != null ) || !"".equals( phy.getType() ) || !"".equals( phy.getDistanceUnit() ) ) {
            return TestFail.here( String.valueOf( phy.getIdentifier() ) + "|" + phy.getType() );
        }
        // the tree structure is untouched by a metadata write
        if ( phy.getNumberOfExternalNodes() != 2 ) {
            return TestFail.here();
        }
        return true;
    }

    private static boolean equality() {
        final TreePropertiesDraft a = TreePropertiesDraft.from( tree() );
        final TreePropertiesDraft b = TreePropertiesDraft.from( tree() );
        if ( !a.equals( b ) || ( a.hashCode() != b.hashCode() ) ) {
            return TestFail.here();
        }
        b.distanceUnit = "x";
        if ( a.equals( b ) ) {
            return TestFail.here();
        }
        if ( a.equals( null ) || a.equals( "string" ) ) {
            return TestFail.here();
        }
        if ( !a.toString().contains( "My tree" ) ) {
            return TestFail.here( a.toString() );
        }
        return true;
    }

    private TreePropertiesDraftTest() {
    }
}
