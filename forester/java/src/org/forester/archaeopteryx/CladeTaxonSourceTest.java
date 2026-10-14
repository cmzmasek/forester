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
import java.util.Collections;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * Headful: a TreePanel keeps the {@link TaxonSource} the clade bands and the branch colorization were made with,
 * and every later rebuild uses it -- the figure must not change under the user because something was cached.
 * Fixture: an "order Rodentia" annotation over a tip the (faked) database places in Carnivora, so the TREE source
 * draws 2 clades and the DATABASE source 3. A green no-op where there is no display.
 */
public final class CladeTaxonSourceTest {

    public static void main( final String[] args ) {
        System.out.println( "CladeTaxonSource: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        TreePanelUtil.setDefaultLineageServiceForTest( new CatService() );
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { fixture() }, new Configuration(), "cladesource" ) );
            final boolean[] ok = { true };
            try {
                SwingUtilities.invokeAndWait( () -> ok[ 0 ] = run( mf[ 0 ].getMainPanel().getCurrentTreePanel() ) );
            }
            finally {
                SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            }
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        finally {
            TreePanelUtil.setDefaultLineageServiceForTest( null );
        }
    }

    private static boolean run( final TreePanel tp ) {
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        final java.util.List<CladeLevel.Spec> order = Collections
                .singletonList( new CladeLevel.Spec( "order", TreePanel.CLADE_LABEL_ANGLE.VERTICAL ) );
        // TREE: Rodentia and Primates, two clades; DATABASE: the cat leaves Rodentia -> three
        final int from_tree = tp.setCladeLevels( order, TreePanel.CLADE_VIS.BOXES, false, TaxonSource.TREE );
        if ( ( from_tree != 2 ) || ( tp.getCladeTaxonSource() != TaxonSource.TREE ) ) {
            return fail( "from the tree: 2 clades and the source remembered; got " + from_tree + " / " + tp.getCladeTaxonSource() );
        }
        tp.rebuildCladeBands(); // what an edit, an undo or a navigation does
        if ( tp.cladeBandCount() != 2 ) {
            return fail( "a rebuild must use the remembered source; got " + tp.cladeBandCount() + " bands" );
        }
        final int from_db = tp.setCladeLevels( order, TreePanel.CLADE_VIS.BOXES, false, TaxonSource.DATABASE );
        if ( ( from_db != 3 ) || ( tp.getCladeTaxonSource() != TaxonSource.DATABASE ) ) {
            return fail( "from the database: 3 clades; got " + from_db + " / " + tp.getCladeTaxonSource() );
        }
        if ( ( tp.setCladeLevels( order, TreePanel.CLADE_VIS.BOXES, false ) != 3 )
                || ( tp.getCladeTaxonSource() != TaxonSource.DATABASE ) ) {
            return fail( "the source-less call is the DATABASE source (the long-standing behaviour)" );
        }
        if ( ( tp.setCladeLevels( order, TreePanel.CLADE_VIS.BOXES, false, null ) != 3 )
                || ( tp.getCladeTaxonSource() != TaxonSource.DATABASE ) ) {
            return fail( "a null source is the DATABASE source" );
        }
        // the branch colorization keeps its own source the same way
        if ( ( tp.colorByRank( "order", TaxonSource.TREE ) != 2 ) || ( tp.getBranchRankColorizeSource() != TaxonSource.TREE ) ) {
            return fail( "colorize from the tree: 2 clades, source remembered; got " + tp.getBranchRankColorizeSource() );
        }
        if ( ( tp.colorByRank( "order", TaxonSource.DATABASE ) != 3 ) || ( tp.colorByRank( "order" ) != 3 )
                || ( tp.getBranchRankColorizeSource() != TaxonSource.DATABASE ) ) {
            return fail( "colorize from the database: 3 clades; the source-less call is DATABASE" );
        }
        // the lineage write-back: the one tip the (faked) database knows gets its ncbi: field, with provenance
        tp.setEdited( false );
        final String before = tp.getPhylogeny().getDescription();
        if ( tp.writeLineageMetadata() != 1 ) {
            return fail( "one tip has a looked-up lineage to keep" );
        }
        final PhylogenyNode cat = tp.getPhylogeny().getNode( "Felis catus" );
        boolean has_field = false;
        for( final org.forester.phylogeny.data.Property p : cat.getNodeData().getProperties().getProperties() ) {
            has_field |= "ncbi:order".equals( p.getRef() ) && "Carnivora".equals( p.getValue() );
        }
        if ( !has_field ) {
            return fail( "the cat must carry ncbi:order = Carnivora; got " + cat.getNodeData().getProperties() );
        }
        if ( !tp.isEdited() || ( tp.getPhylogeny().getDescription() == null )
                || !tp.getPhylogeny().getDescription().contains( "Used NCBI Taxonomy to write the looked-up lineage of 1 tip" )
                || ( ( before != null ) && !tp.getPhylogeny().getDescription().startsWith( before ) ) ) {
            return fail( "the write marks the tree edited and appends its provenance; got " + tp.getPhylogeny().getDescription() );
        }
        // from now on the TREE source sees the cat as a carnivore: three clades without any lookup
        if ( tp.setCladeLevels( order, TreePanel.CLADE_VIS.BOXES, false, TaxonSource.TREE ) != 3 ) {
            return fail( "the written field is the tip's own identity for the tree source" );
        }
        return true;
    }

    private static Phylogeny fixture() {
        final PhylogenyNode rodents = annotated( "order", "Rodentia" );
        rodents.addAsChild( tip( "Felis catus" ) ); // the database knows this one is a carnivore
        rodents.addAsChild( tip( "u" ) );
        final PhylogenyNode primates = annotated( "order", "Primates" );
        primates.addAsChild( tip( "a" ) );
        primates.addAsChild( tip( "b" ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( rodents );
        root.addAsChild( primates );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode annotated( final String rank, final String sci ) {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        try {
            t.setRank( rank );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode tip( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    /** Cache-only: knows the cat's order and nothing else; never fetches. */
    private static final class CatService implements TaxonomicLineageService {

        @Override
        public TaxonLineage lineageOf( final String taxon ) {
            return "Felis catus".equalsIgnoreCase( taxon ) ? new TaxonLineage( null, null, null, null,
                    Collections.singletonList( new TaxonLineage.Ancestor( "Carnivora", "order", "33554" ) ) ) : null;
        }

        @Override
        public TaxonLineage fetch( final String taxon ) {
            throw new UnsupportedOperationException( "cache-only" );
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  CladeTaxonSource test failed: " + msg );
        return false;
    }

    private CladeTaxonSourceTest() {
    }
}
