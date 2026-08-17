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

import java.awt.GraphicsEnvironment;
import java.awt.Toolkit;
import java.awt.event.KeyEvent;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JComponent;
import javax.swing.KeyStroke;

import org.forester.archaeopteryx.TanglegramLinker.LinkField;
import org.forester.archaeopteryx.TanglegramPanel.ConnectorColorMode;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Tests the tanglegram GUI glue: the tree-index selection logic ({@link MainFrameApplication#tanglegramTreeIndices},
 * headless) and that a {@link TanglegramFrame} constructs its panel with the correct link result and title (headful,
 * a green no-op when headless).
 */
public final class TanglegramToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TanglegramTool: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( !selectionLogicOk() ) {
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return frameOk() && frameWiringOk() && colorSelectorOk() && crossFieldFrameOk() && associationFrameOk()
                && phylogramToggleOk() && fontControlOk();
    }

    /** The toolbar splits into two rows (was too wide), and its font spinner resizes the panel's tip-label font live. */
    private static boolean fontControlOk() {
        final TanglegramFrame frame = new TanglegramFrame( treeABC(), treeCBA(), LinkField.NODE_NAME, "L", "R" );
        try {
            final TanglegramPanel panel = frame.getTanglegramPanel();
            if ( frame.toolbarRowCountForTest() != 2 ) {
                return fail( "the toolbar should have 2 rows, got " + frame.toolbarRowCountForTest() );
            }
            frame.setFontSizeForTest( 22 );
            if ( panel.getLabelFontSize() != 22 ) {
                return fail( "the font spinner should resize the panel label font to 22, got "
                        + panel.getLabelFontSize() );
            }
            frame.setFontSizeForTest( 9 );
            if ( panel.getLabelFontSize() != 9 ) {
                return fail( "the font spinner should update the panel font again to 9, got "
                        + panel.getLabelFontSize() );
            }
            return true;
        }
        finally {
            frame.dispose();
        }
    }

    /** A frame built with an association table links differently-named host/parasite trees, marks the panel as
     *  association-linked, and notes it in the summary; a value join over the same trees links nothing. */
    private static boolean associationFrameOk() {
        final Phylogeny hosts = tree( clade( leaf( "gopherA" ), leaf( "gopherB" ) ) );
        final Phylogeny lice = tree( clade( leaf( "louseB" ), leaf( "louseA" ) ) ); // reversed -> a crossing
        final Map<String, List<String>> assoc = new LinkedHashMap<>();
        assoc.put( "gopherA", Arrays.asList( "louseA" ) );
        assoc.put( "gopherB", Arrays.asList( "louseB" ) );
        final TanglegramFrame frame = new TanglegramFrame( hosts, lice, LinkField.NODE_NAME, LinkField.NODE_NAME, assoc,
                                                           "hosts", "lice" );
        try {
            final TanglegramPanel panel = frame.getTanglegramPanel();
            if ( !panel.isAssociationLinked() ) {
                return fail( "an association-linked panel should report isAssociationLinked()" );
            }
            if ( ( panel.getResult().getLinks().size() != 2 ) || ( panel.getUnmatchedCount() != 0 ) ) {
                return fail( "the association frame should link both pairs (2 links, 0 unmatched), got "
                        + panel.getResult().getLinks().size() + " links, " + panel.getUnmatchedCount() + " unmatched" );
            }
            if ( !frame.summaryTextForTest().contains( "via association file" ) ) {
                return fail( "the summary should note the association file, was: " + frame.summaryTextForTest() );
            }
        }
        finally {
            frame.dispose();
        }
        // the whole point: a value join over the same (differently-named) trees links nothing
        final TanglegramFrame plain = new TanglegramFrame( hosts, lice, LinkField.NODE_NAME, "hosts", "lice" );
        try {
            if ( !plain.getTanglegramPanel().getResult().getLinks().isEmpty() ) {
                return fail( "a value join must not link the differently-named trees" );
            }
        }
        finally {
            plain.dispose();
        }
        return true;
    }

    /** The aligned-phylogram checkbox is enabled only when a tree has branch lengths, and toggling it drives the panel. */
    private static boolean phylogramToggleOk() {
        final TanglegramFrame frame = new TanglegramFrame( branchLengthTree(), branchLengthTree(), LinkField.NODE_NAME,
                                                           "L", "R" );
        try {
            final TanglegramPanel panel = frame.getTanglegramPanel();
            if ( !frame.isPhylogramCheckboxEnabledForTest() ) {
                return fail( "the phylogram checkbox should be enabled when the trees have branch lengths" );
            }
            if ( panel.isPhylogram() ) {
                return fail( "the phylogram should start off (cladogram by default)" );
            }
            frame.clickPhylogramForTest();
            if ( !panel.isPhylogram() ) {
                return fail( "clicking the phylogram checkbox should turn the phylogram on" );
            }
        }
        finally {
            frame.dispose();
        }
        // length-less trees -> the toggle is greyed out (a cladogram is the only option)
        final TanglegramFrame no_bl = new TanglegramFrame( treeABC(), treeCBA(), LinkField.NODE_NAME, "L", "R" );
        try {
            if ( no_bl.isPhylogramCheckboxEnabledForTest() ) {
                return fail( "the phylogram checkbox should be disabled when there are no branch lengths" );
            }
        }
        finally {
            no_bl.dispose();
        }
        return true;
    }

    private static Phylogeny branchLengthTree() {
        final PhylogenyNode a = leaf( "A" );
        a.setDistanceToParent( 1.0 );
        final PhylogenyNode b = leaf( "B" );
        b.setDistanceToParent( 2.0 );
        return tree( clade( a, b ) );
    }

    /** A frame built with a DIFFERENT link field per tree links two trees that store the same value in different
     *  fields (gene tree's scientific name to species tree's node name), and its summary names both fields. */
    private static boolean crossFieldFrameOk() {
        // gene tree: node name = accession, species in the taxonomy scientific name; species tree: species = node name
        final Phylogeny gene = tree( clade( geneLeaf( "acc1", "Homo sapiens" ), geneLeaf( "acc2", "Mus musculus" ) ) );
        final Phylogeny species = tree( clade( leaf( "Mus musculus" ), leaf( "Homo sapiens" ) ) );
        final TanglegramFrame frame = new TanglegramFrame( gene, species, LinkField.SCIENTIFIC_NAME,
                                                           LinkField.NODE_NAME, "gene", "species" );
        try {
            final TanglegramPanel panel = frame.getTanglegramPanel();
            if ( ( panel.getResult().getLinks().size() != 2 ) || ( panel.getUnmatchedCount() != 0 ) ) {
                return fail( "the two-field frame should link both species (2 links, 0 unmatched), got "
                        + panel.getResult().getLinks().size() + " links, " + panel.getUnmatchedCount() + " unmatched" );
            }
            if ( ( panel.getLeftField() != LinkField.SCIENTIFIC_NAME )
                    || ( panel.getRightField() != LinkField.NODE_NAME ) ) {
                return fail( "the panel should carry both link fields" );
            }
            final String summary = frame.summaryTextForTest();
            if ( !summary.contains( LinkField.SCIENTIFIC_NAME.label() + " ↔ " + LinkField.NODE_NAME.label() ) ) {
                return fail( "summary should name both link fields, was: " + summary );
            }
            return true;
        }
        finally {
            frame.dispose();
        }
    }

    private static PhylogenyNode geneLeaf( final String accession, final String scientific_name ) {
        final PhylogenyNode n = leaf( accession ); // node name is the accession, distinct from the species
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static boolean colorSelectorOk() {
        final TanglegramFrame frame = new TanglegramFrame( taxonTree(), taxonTree(), LinkField.SCIENTIFIC_NAME, "L",
                                                           "R" );
        try {
            if ( frame.colorItemCountForTest() != 3 ) { // Uniform, Crossings, Taxonomy: Scientific Name
                return fail( "the colour selector should have 3 items, got " + frame.colorItemCountForTest() );
            }
            final TanglegramPanel panel = frame.getTanglegramPanel();
            frame.selectColorForTest( 1 );
            if ( panel.colorModeForTest() != ConnectorColorMode.CROSSINGS ) {
                return fail( "selector index 1 should set CROSSINGS mode" );
            }
            frame.selectColorForTest( 2 );
            if ( panel.colorModeForTest() != ConnectorColorMode.FIELD ) {
                return fail( "selector index 2 should set a FIELD mode" );
            }
            frame.selectColorForTest( 0 );
            if ( panel.colorModeForTest() != ConnectorColorMode.UNIFORM ) {
                return fail( "selector index 0 should set UNIFORM mode" );
            }
            return true;
        }
        finally {
            frame.dispose();
        }
    }

    private static Phylogeny taxonTree() {
        return tree( clade( taxonLeaf( "Homo sapiens" ), taxonLeaf( "Mus musculus" ) ) );
    }

    private static PhylogenyNode taxonLeaf( final String scientific_name ) {
        final PhylogenyNode n = leaf( scientific_name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static boolean selectionLogicOk() {
        // the picks are honoured directly (the pickers are always shown so a directional association can be reordered)
        if ( !Arrays.equals( new int[] { 0, 1 }, MainFrameApplication.tanglegramTreeIndices( 0, 1 ) ) ) {
            return fail( "picks (0,1) should resolve to {0,1}" );
        }
        // reversible: the first tree drives the association's left column, so (1,0) must be distinct from (0,1)
        if ( !Arrays.equals( new int[] { 1, 0 }, MainFrameApplication.tanglegramTreeIndices( 1, 0 ) ) ) {
            return fail( "picks (1,0) should resolve to {1,0} (order is meaningful for association linking)" );
        }
        if ( !Arrays.equals( new int[] { 2, 3 }, MainFrameApplication.tanglegramTreeIndices( 2, 3 ) ) ) {
            return fail( "picks (2,3) should resolve to {2,3}" );
        }
        // same tree picked twice -> invalid
        if ( MainFrameApplication.tanglegramTreeIndices( 1, 1 ) != null ) {
            return fail( "the same tree picked twice should be rejected (null)" );
        }
        return true;
    }

    private static boolean frameOk() {
        final TanglegramFrame frame = new TanglegramFrame( treeABC(), treeCBA(), LinkField.NODE_NAME, "leftTree",
                                                           "rightTree" );
        try {
            if ( frame.getTanglegramPanel().getResult().getLinks().size() != 3 ) {
                return fail( "the frame's panel should hold 3 links, got "
                        + frame.getTanglegramPanel().getResult().getLinks().size() );
            }
            if ( frame.getTanglegramPanel().getUnmatchedCount() != 0 ) {
                return fail( "the frame's panel should have 0 unmatched tips" );
            }
            if ( !frame.getTitle().contains( "leftTree" ) || !frame.getTitle().contains( "rightTree" ) ) {
                return fail( "the frame title should name both trees, was: " + frame.getTitle() );
            }
            return true;
        }
        finally {
            frame.dispose();
        }
    }

    private static boolean frameWiringOk() {
        final TanglegramFrame frame = new TanglegramFrame( treeABC(), treeCBA(), LinkField.NODE_NAME, "L", "R" );
        try {
            if ( frame.isUndoButtonEnabledForTest() || frame.isRedoButtonEnabledForTest() ) {
                return fail( "undo/redo buttons should start disabled" );
            }
            final TanglegramPanel panel = frame.getTanglegramPanel();
            panel.rotateLeftRootForTest(); // fires the change listener -> frame.refresh()
            if ( !frame.isUndoButtonEnabledForTest() || frame.isRedoButtonEnabledForTest() ) {
                return fail( "after a flip: undo enabled, redo disabled" );
            }
            if ( !frame.summaryTextForTest().contains( panel.getCrossingCount() + " crossings" ) ) {
                return fail( "summary should show the live crossing count, was: " + frame.summaryTextForTest() );
            }
            frame.clickUndoForTest(); // undo via the toolbar button
            if ( frame.isUndoButtonEnabledForTest() || !frame.isRedoButtonEnabledForTest() ) {
                return fail( "after the undo button: redo enabled, undo disabled" );
            }
            frame.clickRedoForTest();
            if ( !frame.isUndoButtonEnabledForTest() ) {
                return fail( "after the redo button: undo should be enabled again" );
            }
            final int mask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx();
            final Object undo_key = frame.getRootPane().getInputMap( JComponent.WHEN_IN_FOCUSED_WINDOW )
                    .get( KeyStroke.getKeyStroke( KeyEvent.VK_Z, mask ) );
            if ( !"tangle-undo".equals( undo_key ) ) {
                return fail( "Cmd/Ctrl+Z should be bound to undo, was: " + undo_key );
            }
            return true;
        }
        finally {
            frame.dispose();
        }
    }

    private static boolean fail( final String message ) {
        System.out.println( "TanglegramTool test failed: " + message );
        return false;
    }

    private static Phylogeny treeABC() {
        return tree( clade( leaf( "A" ), leaf( "B" ), leaf( "C" ) ) );
    }

    private static Phylogeny treeCBA() {
        return tree( clade( leaf( "C" ), leaf( "B" ), leaf( "A" ) ) );
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode clade( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode child : children ) {
            n.addAsChild( child );
        }
        return n;
    }

    private static Phylogeny tree( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }
}
