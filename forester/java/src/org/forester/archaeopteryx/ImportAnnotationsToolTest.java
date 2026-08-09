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
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.tools.NodeDataImporter;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Integration test for the "Import Annotations (CSV/TSV)" tool: dogfoods the import-annotations demo (a plain-named
 * tip tree + its companion CSV) on a real {@link MainFrame}/{@link TreePanel}, drives the UI-free apply+refresh, and
 * asserts the keyed join wrote the columns onto the tips, surfaced them in "Color by", and appended a provenance
 * sentence -- plus a negative control (matching by taxonomy, which these plain tips lack, annotates nothing).
 * Headful; a green no-op when headless.
 */
public final class ImportAnnotationsToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ImportAnnotationsTool: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( !buildColumnPlanOk() || !provenanceOk() ) {
            return false; // pure, headless-safe checks
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final File dir = new File( System.getProperty( "user.dir" ), "forester/demo" );
            final File tree_file = new File( dir, "import-annotations.xml" );
            final File csv_file = new File( dir, "import-annotations.csv" );
            if ( !tree_file.exists() || !csv_file.exists() ) {
                return fail( "demo files missing under " + dir.getAbsolutePath() );
            }
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                    .create( tree_file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
            final String csv = java.nio.file.Files.readString( csv_file.toPath() );
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "import" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Phylogeny live = tp.getPhylogeny();
                    final NodeDataImporter.Table table = NodeDataImporter.parseTable( csv );
                    final int key_col = table.defaultKeyColumn();
                    // NEGATIVE control: matching by taxonomy (which these plain-named tips lack) annotates nothing
                    final NodeDataImporter.ImportResult none = frame.importAnnotationsAndRefit( live, table, key_col,
                            NodeDataImporter.MatchBy.TAXONOMY_ID, NodeDataImporter.ColumnPlan.importAll( table ),
                            "import-annotations.csv" );
                    if ( none.getTipsAnnotated() != 0 ) {
                        fail( ok, "a taxonomy match on plain-named tips should annotate nothing, got "
                                + none.getTipsAnnotated() );
                    }
                    if ( tip( live, "isolate_01" ).getNodeData().getProperties() != null ) {
                        fail( ok, "the failed (taxonomy) match must not have written any property" );
                    }
                    // POSITIVE: match by tip name, with a ColumnPlan that RENAMES the "reads" column -> data:depth
                    // (headers: name,host,country,reads -> column 3), exercising the plan through the refit
                    final NodeDataImporter.ColumnPlan plan = NodeDataImporter.ColumnPlan.importAll( table );
                    plan.setHeader( 3, "depth" );
                    final NodeDataImporter.ImportResult res = frame.importAnnotationsAndRefit( live, table, key_col,
                            NodeDataImporter.MatchBy.TIP_NAME, plan, "import-annotations.csv" );
                    if ( res.getTipsAnnotated() != 12 ) {
                        fail( ok, "a tip-name join should annotate all 12 tips, got " + res.getTipsAnnotated() );
                    }
                    final PhylogenyNode t1 = tip( live, "isolate_01" );
                    if ( !"mosquito".equals( propertyValue( t1, "data:host" ) )
                            || !"1200".equals( propertyValue( t1, "data:depth" ) )
                            || ( propertyValue( t1, "data:reads" ) != null ) ) {
                        fail( ok, "isolate_01 should carry data:host=mosquito and the RENAMED data:depth=1200 (not data:reads)" );
                    }
                    // the quoted CSV field with an embedded comma survived the parse + join
                    if ( !"Congo, DR".equals( propertyValue( tip( live, "isolate_07" ), "data:country" ) ) ) {
                        fail( ok, "the quoted country with an embedded comma should import intact: "
                                + propertyValue( tip( live, "isolate_07" ), "data:country" ) );
                    }
                    // the imported columns are now offered in the "Color by" dropdown (immediately usable)
                    final List<String> refs = tp.getControlPanel().colorByPropertyRefs();
                    if ( !refs.contains( "data:host" ) || !refs.contains( "data:depth" ) ) {
                        fail( ok, "the Color-by dropdown should now offer the imported columns: " + refs );
                    }
                    // a provenance sentence was appended to the tree description
                    final String desc = live.getDescription();
                    if ( ( desc == null ) || !desc.contains( "Imported annotations from table" ) ) {
                        fail( ok, "an import should append a provenance sentence: " + desc );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Pure: MainFrame.buildColumnPlan excludes the key + unchecked columns and applies renames (the one source of
     *  truth the dialog's preview and commit both use). */
    private static boolean buildColumnPlanOk() {
        final NodeDataImporter.Table t = NodeDataImporter.parseTable( "name,host,reads\nA,cat,5\n" );
        final NodeDataImporter.ColumnPlan plan = MainFrame.buildColumnPlan( t, 0, new boolean[] { true, true, false },
                new String[] { "name", "animal", "reads" } );
        if ( plan.isIncluded( 0 ) ) {
            return fail( "buildColumnPlan must exclude the key column" );
        }
        if ( !plan.isIncluded( 1 ) || !"animal".equals( plan.header( 1 ) ) ) {
            return fail( "buildColumnPlan must include + rename a checked column" );
        }
        if ( plan.isIncluded( 2 ) ) {
            return fail( "buildColumnPlan must exclude an unchecked column" );
        }
        return true;
    }

    /** Pure: MainFrame.importProvenance formats the append sentence (plural/singular, source, match-by, columns). */
    private static boolean provenanceOk() {
        final String p = MainFrame.importProvenance( java.util.List.of( "data:host", "data:reads" ), 3, 5, "t.csv",
                NodeDataImporter.MatchBy.TIP_NAME );
        if ( !p.contains( "\"t.csv\"" ) || !p.contains( "onto 3 of 5 tips" ) || !p.contains( "matched by tip name" )
                || !p.contains( "Columns: data:host, data:reads" ) ) {
            return fail( "importProvenance format wrong: " + p );
        }
        final String p1 = MainFrame.importProvenance( java.util.List.of(), 1, 1, "x",
                NodeDataImporter.MatchBy.SEQUENCE_ACCESSION );
        if ( !p1.contains( "onto 1 of 1 tip " ) || p1.contains( "Columns:" ) ) {
            return fail( "importProvenance singular / no-columns case wrong: " + p1 );
        }
        return true;
    }

    private static PhylogenyNode tip( final Phylogeny phy, final String name ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        throw new IllegalStateException( "no tip named " + name );
    }

    private static String propertyValue( final PhylogenyNode n, final String ref ) {
        if ( n.getNodeData().getProperties() == null ) {
            return null;
        }
        final List<Property> ps = n.getNodeData().getProperties().getPropertiesWithGivenRef( ref );
        return ( ( ps != null ) && !ps.isEmpty() ) ? ps.get( 0 ).getValue() : null;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ImportAnnotationsToolTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ImportAnnotationsToolTest] " + msg );
        ok[ 0 ] = false;
    }

    private ImportAnnotationsToolTest() {
    }
}
