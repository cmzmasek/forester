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
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Collections;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headful test for {@link ReportWindow}: it shows the report's summary and every row, exports the
 * report's TSV byte for byte, and closes. A green no-op where there is no display.
 */
public final class ReportWindowTest {

    public static void main( final String[] args ) {
        System.out.println( "ReportWindow: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { false };
        try {
            SwingUtilities.invokeAndWait( () -> ok[ 0 ] = run() );
        }
        catch ( final Exception e ) {
            return fail( "could not run on the EDT: " + e );
        }
        return ok[ 0 ];
    }

    private static boolean run() {
        final PhylogenyNode carnivora = annotated( "order", "Carnivora" );
        carnivora.addAsChild( tip( "cat" ) );
        carnivora.addAsChild( tip( "dog" ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( carnivora );
        root.addAsChild( tip( "lonely" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        final CladeAnnotationReport report = CladeAnnotationReport.of( phy, Collections.singletonList( "order" ), null,
                                                                       TaxonSource.TREE );
        if ( report.getRows().size() != 3 ) {
            return fail( "precondition: three tip rows, got " + report.getRows().size() );
        }
        final File[] remembered = { null };
        final ReportWindow w = new ReportWindow( null, "report test", report, null,
                                                                               dir -> remembered[ 0 ] = dir );
        try {
            w.setVisible( true );
            if ( !w.isShowing() ) {
                return fail( "the window must show" );
            }
            if ( w.tableForTest().getRowCount() != report.getRows().size() ) {
                return fail( "the table must hold every row; got " + w.tableForTest().getRowCount() );
            }
            if ( !"cat".equals( w.tableForTest().getValueAt( 0, 1 ) ) || !"Carnivora".equals( w.tableForTest().getValueAt( 0, 2 ) )
                    || !"".equals( w.tableForTest().getValueAt( 2, 2 ) ) ) {
                return fail( "the table's cells are the rows' node, taxon, source and note" );
            }
            if ( w.tableForTest().isCellEditable( 0, 0 ) ) {
                return fail( "the table is read-only" );
            }
            if ( !report.getSummary().equals( w.summaryForTest() ) || !w.summaryForTest().contains( "order: 2 of 3 tips placed" ) ) {
                return fail( "the summary area shows the report's summary; got " + w.summaryForTest() );
            }
            final File f = File.createTempFile( "clade-report", ".tsv" );
            f.deleteOnExit();
            w.exportTo( f );
            final String written = new String( Files.readAllBytes( f.toPath() ), StandardCharsets.UTF_8 );
            if ( !written.equals( report.toTsv() ) ) {
                return fail( "the export must be the report's TSV byte for byte" );
            }
            if ( remembered[ 0 ] != null ) {
                return fail( "a programmatic export remembers no directory (only the interactive one does)" );
            }
            if ( w.getReport() != report ) {
                return fail( "the window keeps the report it was given" );
            }
            if ( w.actionButtonForTest() != null ) {
                return fail( "a window without an action has no action button" );
            }
            // an action button runs its action once and then disables itself
            final int[] ran = { 0 };
            final ReportWindow a = new ReportWindow( null, "action test", report, null, null, "Apply", () -> ran[ 0 ]++ );
            try {
                a.setVisible( true );
                a.actionButtonForTest().doClick();
                a.actionButtonForTest().doClick(); // disabled: nothing
                if ( ( ran[ 0 ] != 1 ) || a.actionButtonForTest().isEnabled() ) {
                    return fail( "the action runs once and the button disables itself; ran " + ran[ 0 ] );
                }
            }
            finally {
                a.dispose();
            }
        }
        catch ( final Exception e ) {
            return fail( "unexpected: " + e );
        }
        finally {
            w.dispose();
        }
        return true;
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

    private static boolean fail( final String msg ) {
        System.out.println( "  ReportWindow test failed: " + msg );
        return false;
    }

    private ReportWindowTest() {
    }
}
