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
import java.io.File;

import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Headful: the "Check Internal Taxa" tool on its demo -- the report opens with the four findings and an "Apply 3
 * corrections" button, applying writes the proposals onto the LIVE tree with provenance and the edited mark, and
 * Undo takes them back. A green no-op where there is no display.
 */
public final class InternalTaxaCheckToolTest {

    public static void main( final String[] args ) {
        System.out.println( "InternalTaxaCheckTool: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final File f = new File( System.getProperty( "user.dir" ), "forester/demo/internal-taxa-check.xml" );
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( f, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "internaltaxa" ) );
            final boolean[] ok = { true };
            try {
                SwingUtilities.invokeAndWait( () -> ok[ 0 ] = run( mf[ 0 ] ) );
            }
            finally {
                SwingUtilities.invokeAndWait( () -> mf[ 0 ].dispose() );
            }
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean run( final MainFrame mf ) {
        final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
        final Phylogeny phy = tp.getPhylogeny();
        final PhylogenyNode rodentia = annotatedNode( phy, "Rodentia" );
        if ( rodentia == null ) {
            return fail( "precondition: the demo carries the wrong 'Rodentia' annotation" );
        }
        if ( mf._check_internal_taxa == null ) {
            return fail( "the Analysis menu must carry the check" );
        }
        final ReportWindow w = mf.runInternalTaxaCheck( tp, null, null ); // offline: the fields judge
        try {
            if ( !w.isShowing() || !( w.getReport() instanceof InternalTaxonomyAuditReport ) ) {
                return fail( "the report window must show the audit" );
            }
            final InternalTaxonomyAuditReport report = (InternalTaxonomyAuditReport) w.getReport();
            if ( ( report.getFindings().size() != 4 ) || ( report.proposalCount() != 3 )
                    || ( w.tableForTest().getRowCount() != 4 ) ) {
                return fail( "four findings, three proposals; got " + report.getSummary() );
            }
            if ( ( w.actionButtonForTest() == null ) || !"Apply 3 corrections".equals( w.actionButtonForTest().getText() ) ) {
                return fail( "the action button applies the three corrections; got " + w.actionButtonForTest() );
            }
            final boolean was_edited = tp.isEdited();
            tp.setEdited( false );
            final int changed = tp.applyInternalTaxonomyCorrections( report.getFindings() ); // what the button runs
            if ( changed != 3 ) {
                return fail( "three corrections applied; got " + changed );
            }
            if ( annotatedNode( phy, "Rodentia" ) != null || !"Carnivora".equals( rodentia.getNodeData().getTaxonomy().getScientificName() )
                    || !"order".equals( rodentia.getNodeData().getTaxonomy().getRank() ) ) {
                return fail( "the wrong taxon must be replaced on the LIVE tree; got " + rodentia.getNodeData().getTaxonomy() );
            }
            final PhylogenyNode felidae = annotatedNode( phy, "Felidae" );
            if ( ( felidae == null ) || !"family".equals( felidae.getNodeData().getTaxonomy().getRank() ) ) {
                return fail( "the wrong rank must be corrected in place" );
            }
            if ( !tp.isEdited() || ( phy.getDescription() == null )
                    || !phy.getDescription().contains( "Used the internal-taxonomy check (tip lineages and NCBI Taxonomy) to correct 3 internal node taxonomies" ) ) {
                return fail( "the corrections mark the tree edited and append their provenance; got " + phy.getDescription() );
            }
            if ( !tp.canUndo() || !tp.undo() ) {
                return fail( "the corrections must be undoable" );
            }
            if ( annotatedNode( tp.getPhylogeny(), "Rodentia" ) == null ) {
                return fail( "undo must bring the wrong annotation back" );
            }
            tp.setEdited( was_edited );
        }
        finally {
            w.dispose();
        }
        return true;
    }

    private static PhylogenyNode annotatedNode( final Phylogeny phy, final String sci ) {
        for( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && n.getNodeData().isHasTaxonomy()
                    && sci.equals( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                return n;
            }
        }
        return null;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  InternalTaxaCheckTool test failed: " + msg );
        return false;
    }

    private InternalTaxaCheckToolTest() {
    }
}
